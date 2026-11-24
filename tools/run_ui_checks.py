"""Check actual Qt4 Load/Run/Save in an isolated host, with its embedded Python runtime."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'tests'))
from current_host_checks import read_grid, close


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--host', required=True, type=Path)
    parser.add_argument('--suite', required=True, type=Path, help='Successful current_host_checks output')
    parser.add_argument('--exe', type=Path, default=ROOT / 'build/x64/qt4_ui_smoke.exe')
    parser.add_argument('--out', required=True, type=Path)
    args = parser.parse_args()
    host, suite, out, exe = (p.resolve() for p in (args.host, args.suite, args.out, args.exe))
    if not (host / 'HOSC_ISOLATED_TEST_HOST.txt').is_file():
        parser.error('Stage an isolated test host first')
    if out.exists() or host == out or host in out.parents or not str(out).isascii():
        parser.error('--out must be a new ASCII directory outside the host')
    cases = [('HOSC_TEMPLATE', 'template_exact'), ('HOSC_WINDOW', 'window_grid')]
    for name, case in cases:
        if not (suite / case / 'saved.sgems').is_file():
            parser.error('First run tests/current_host_checks.py in the same isolated host')
    if not exe.is_file():
        parser.error('Build with --with-tests first')
    out.mkdir(parents=True)
    for name, case in cases:
        command = (suite / case / 'run.cmd').read_text().splitlines()[1]
        node = ET.fromstring(command.split('::', 2)[2])
        node.find('props').set('value', 'ui_smoke_' + name)
        node.find('audit_output_directory').set('value', str(out))
        weighting = node.find('enable_weighted_moment')
        if name == 'HOSC_TEMPLATE':
            if weighting is None:
                weighting = ET.SubElement(node, 'enable_weighted_moment')
            weighting.set('value', '1')
        ET.ElementTree(node).write(out / (name + '_run.par'), encoding='unicode')
    env = os.environ.copy()
    env.update(GSTLAPPLIHOME=str(host), PYTHONHOME=str(host), PYTHONPATH=str(host))
    env['PATH'] = str(host) + ';' + env['PATH']
    result = subprocess.run([str(exe), str(host), str(out), str(suite)], cwd=host, env=env,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=90,
                            creationflags=subprocess.CREATE_NO_WINDOW)
    (out / 'host.log').write_bytes(result.stdout)
    print(result.stdout.decode(errors='replace')[-3500:])
    if result.returncode:
        raise SystemExit(result.returncode)
    for name, case in cases:
        expected = read_grid(suite / case / 'saved.sgems')['props'][case]
        actual = read_grid(out / (name + '_ui_result.sgems'))['props']['ui_smoke_' + name]
        assert len(actual) == len(expected)
        for a, b in zip(actual, expected):
            close(a, b, name + ' actual panel result')
    report = dict(status='passed', actual_panel_load_run_numerics='passed',
                  screenshots=[p.name for p in out.glob('*.png')])
    (out / 'report.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    print('PASS both actual SGeMS panel results match independently checked host cases')


if __name__ == '__main__':
    main()
