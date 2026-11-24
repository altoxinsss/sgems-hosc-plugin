"""Copy a small isolated SGeMS runtime; never modify the supplied source host."""
import argparse
import csv
import hashlib
from pathlib import Path
import shutil

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--host', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--dll-dir', type=Path, default=ROOT/'bin/win64')
    args = parser.parse_args()
    source, dest = args.host.resolve(), args.output.resolve()
    if dest.exists() or source == dest or source in dest.parents or dest in source.parents:
        parser.error('--output must be a new directory outside the source host')
    if not (source/'sgems-x64.exe').is_file():
        parser.error('Missing sgems-x64.exe')
    with (ROOT/'bin/win64/host_libraries.csv').open(encoding='utf-8-sig') as f:
        for row in csv.DictReader(f):
            if hashlib.sha256((source/row['file']).read_bytes()).hexdigest() != row['sha256']:
                raise SystemExit('Unsupported host library: '+row['file'])
    overlays = [(args.dll_dir/name, name) for name in ('HOSC2008.dll', 'HOSCWindow.dll')]
    overlays += [(ROOT/'ui'/name, name) for name in ('HOSC_TEMPLATE.ui', 'HOSC_WINDOW.ui')]
    if not all(p.is_file() for p, _ in overlays):
        parser.error('Missing plugin DLL or UI')
    dest.mkdir(parents=True)
    for item in source.iterdir():
        # Do not clone a user's experiment archives or other files beside the host.
        runtime = item.name == 'sgems-x64.exe' or item.suffix.lower() in (
            '.dll', '.pyd', '.py', '.bmp', '.png', '.ico', '.wav')
        if item.is_file() and runtime:
            shutil.copy2(item, dest/item.name)
    for name in ('plugins', 'encodings'):
        if (source/name).is_dir():
            shutil.copytree(source/name, dest/name)
    plugin = dest/'plugins/Geostat'
    plugin.mkdir(parents=True, exist_ok=True)
    for path, name in overlays:
        shutil.copy2(path, plugin/name)
    (dest/'HOSC_ISOLATED_TEST_HOST.txt').write_text('Isolated test copy of '+str(source)+'\n', encoding='utf-8')
    print('Isolated host:', dest)


if __name__ == '__main__':
    main()
