"""Build the in-process VC2010/Qt4 plugins against an external compatible SDK."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[1]
HOST_LIBS = ['GsTLAppli_appli', 'GsTLAppli_geostat', 'GsTLAppli_grid',
             'GsTLAppli_utils', 'GsTLAppli_math', 'GsTLAppli_gui',
             'QtCore4', 'QtGui4', 'QtXml4']


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sdk-root', type=Path, default=os.environ.get('HOSC_SDK_ROOT'))
    parser.add_argument('--build-dir', type=Path, default=ROOT/'build/x64')
    parser.add_argument('--check-only', action='store_true')
    parser.add_argument('--with-tests', action='store_true')
    args = parser.parse_args()
    if args.sdk_root is None:
        parser.error('Supply --sdk-root or HOSC_SDK_ROOT; see docs/SDK.md')
    sdk, build = args.sdk_root.resolve(), args.build_dir.resolve()
    vc = sdk/'winsdk71_extracted/vc_stdamd64/Program Files(64)/Microsoft Visual Studio 10.0/VC'
    vcinc = sdk/'winsdk71_extracted/vc_stdx86/Program Files/Microsoft Visual Studio 10.0/VC/include'
    win = sdk/'winsdk71_extracted/WinSDKBuild_amd64/Program Files/Microsoft SDKs/Windows/v7.1'
    qt = sdk/'qt-4.7.3'
    src = ROOT/'src'
    includes = [src/'native', src/'native/vendor',
                sdk/'SGeMS-5d87afac331e12b5cdf6c1c1b2f0e6e20e9485d7',
                sdk/'gstl', sdk/'boost_1_47_0', qt/'include',
                qt/'include/QtCore', qt/'include/QtGui', qt/'include/QtXml',
                qt/'mkspecs/win32-msvc2010']
    extra = ['GsTLAppli_actions', 'GsTLAppli_filters', 'GsTLAppli_extragui',
             'GsTLwidgets', 'QtDesigner4'] if args.with_tests else []
    required = [vc/'bin/amd64/cl.exe', vc/'bin/amd64/link.exe', vcinc/'vector',
                win/'Include/Windows.h', win/'Lib/x64/kernel32.lib',
                qt/'include/QtCore/QObject', *includes,
                *[sdk/'host_imports'/(s+'.lib') for s in HOST_LIBS+extra]]
    missing = [str(p) for p in required if not p.exists()]
    if missing:
        raise SystemExit('Missing SDK inputs:\n'+'\n'.join(missing))
    with (ROOT/'tests/native_core_sha256.csv').open(encoding='utf-8-sig', newline='') as f:
        core = list(csv.DictReader(f))
    for row in core:
        if sha(src/row['path']) != row['sha256']:
            raise SystemExit('Numerical baseline changed: '+row['path'])
    print('SDK present; numerical baseline verified:', len(core), 'files', flush=True)
    if args.check_only:
        return
    if build == ROOT or any(build == ROOT/p or (ROOT/p) in build.parents
                            for p in ('src', 'bin', 'paper', 'examples', 'ui')):
        raise SystemExit('Build output must not overwrite source or release directories')
    build.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env['PATH'] = ';'.join(map(str, [vc/'bin/amd64', sdk/'winsdk71_extracted/vc_stdamd64/Win/System64']))+';'+env['PATH']
    env['INCLUDE'] = ';'.join(map(str, [vcinc, win/'Include']))
    env['LIB'] = ';'.join(map(str, [vc/'lib/amd64', win/'Lib/x64', sdk/'host_imports']))
    flags = ['/nologo', '/c', '/EHsc', '/GR', '/MD', '/O2', '/W3', '/DNDEBUG',
             '/DWIN32', '/D_WIN32_WINNT=0x0601', '/DQT_DLL', '/DQT_SHARED',
             '/DQT_NO_DEBUG', '/D_CRT_SECURE_NO_WARNINGS', '/D_SCL_SECURE_NO_WARNINGS',
             '/D_ITERATOR_DEBUG_LEVEL=0']
    commands, compiled = [], []

    def run(cmd, label):
        cmd = list(map(str, cmd))
        commands.append(cmd)
        proc = subprocess.run(cmd, env=env, cwd=build, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        (build/(label+'.log')).write_bytes(proc.stdout)
        print(label, proc.returncode, flush=True)
        if proc.returncode:
            raise SystemExit(proc.stdout.decode(errors='replace')[-14000:])

    def compile_source(path, label):
        obj = build/(label+'.obj')
        compiled.append(path)
        run([vc/'bin/amd64/cl.exe', *flags, *['/I'+str(p) for p in includes],
             '/Fo'+str(obj), path], label)
        return obj

    objects = [compile_source(p, 'native_'+p.stem) for p in sorted((src/'native').rglob('*.cpp'))]
    objects.append(compile_source(src/'shared/EnginePlugin.cpp', 'EnginePlugin'))
    for name, sub in [('HOSC2008', 'Template/HOSC2008'), ('HOSCWindow', 'MoveWindow/HOSCWindow')]:
        obj = compile_source(src/sub/'Plugin.cpp', name+'_Plugin')
        run([vc/'bin/amd64/link.exe', '/NOLOGO', '/DLL', '/MACHINE:X64', '/INCREMENTAL:NO',
             '/OUT:'+str(build/(name+'.dll')), *objects, obj,
             *[s+'.lib' for s in HOST_LIBS], 'kernel32.lib', 'user32.lib'], name+'_link')
    if args.with_tests:
        for name, libs in [('run_performance_checks', ['kernel32', 'psapi']),
                           ('qt4_ui_smoke', HOST_LIBS+extra+['kernel32', 'user32'])]:
            obj = compile_source(ROOT/'tests'/(name+'.cpp'), name)
            run([vc/'bin/amd64/link.exe', '/NOLOGO', '/MACHINE:X64', '/INCREMENTAL:NO',
                 '/OUT:'+str(build/(name+'.exe')), obj, *[s+'.lib' for s in libs]], name+'_link')
    report = dict(status='built', sdk_root=str(sdk), numerical_baseline_files=len(core),
                  sources={p.relative_to(ROOT).as_posix(): sha(p) for p in compiled},
                  source_tree={p.relative_to(ROOT).as_posix(): sha(p) for p in src.rglob('*') if p.is_file()},
                  outputs={p.name: sha(p) for p in build.iterdir() if p.suffix in ('.dll', '.exe')},
                  commands=commands)
    (build/'build_manifest.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    print('Build complete:', build)


if __name__ == '__main__':
    main()
