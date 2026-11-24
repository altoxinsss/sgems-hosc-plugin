"""Read-only package, numerical-source and evidence checks; no SGeMS or SDK required."""
import argparse
import csv
import hashlib
from pathlib import Path
import re
import struct
from urllib.parse import unquote
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--paper-source', type=Path, help='Optionally compare to the manuscript Data directory')
    args = parser.parse_args()
    errors = []

    def check(good, message):
        if not good:
            errors.append(message)

    for manifest, base, field in [
        ('bin/win64/SHA256SUMS.csv', ROOT / 'bin/win64', 'file'),
        ('tests/native_core_sha256.csv', ROOT / 'src', 'path'),
    ]:
        with (ROOT / manifest).open(encoding='utf-8-sig', newline='') as stream:
            rows = list(csv.DictReader(stream))
        for row in rows:
            path = (base / row[field]).resolve()
            check(base in path.parents and path.is_file() and sha(path) == row['sha256'].lower(),
                  'Checksum: ' + row[field])
        print(manifest, len(rows), 'entries')
    for name in ('HOSC2008.dll', 'HOSCWindow.dll'):
        data = (ROOT / 'bin/win64' / name).read_bytes()
        offset = struct.unpack_from('<I', data, 0x3c)[0]
        check(data[offset:offset+4] == b'PE\0\0' and struct.unpack_from('<H', data, offset+4)[0] == 0x8664,
              'Not an x64 DLL: ' + name)
    for name in ('HOSC_TEMPLATE.ui', 'HOSC_WINDOW.ui'):
        check(sha(ROOT / 'ui' / name) == sha(ROOT / 'bin/win64' / name), 'UI/release mismatch: ' + name)
        check(ET.parse(ROOT / 'ui' / name).getroot().get('version') == '4.0', 'Not a Qt4 form: ' + name)

    data = ROOT / 'paper/Data'
    entries = []
    for line in (data / 'SHA256SUMS.txt').read_text(encoding='utf-8-sig').splitlines():
        digest, relative = line.split('  ', 1)
        relative = relative.replace('\\', '/')
        path = (data / relative).resolve()
        check(data in path.parents and path.is_file() and sha(path) == digest.lower(), 'Data checksum: ' + relative)
        entries.append(relative)
    actual = {p.relative_to(data).as_posix() for p in data.rglob('*') if p.is_file()}
    check(actual == set(entries) | {'SHA256SUMS.txt'}, 'Data manifest coverage differs')
    check(not any(Path(p).suffix.lower() in ('.json', '.py', '.ps1', '.exe') for p in actual), 'Processing files in paper/Data')
    parameters = list((data / 'Parameters').rglob('*.par'))
    check(len(parameters) == 138, 'Expected 138 paper parameters')
    check(len(list((data / 'Results.prj').glob('*.sgems'))) == 17, 'Expected 17 result grids')
    check(len(list((data / 'CSV').glob('*.csv'))) == 14, 'Expected 14 result/performance CSVs')
    for path in parameters:
        node = ET.parse(path).getroot()
        choice = node.find('Algorithm_Button').get('value')
        if 'bell' in path.stem.lower():
            check(choice == 'Bell', 'Bell selection: ' + str(path))
        if 'avl' in path.stem.lower():
            check(choice == 'AVL Based algorithm', 'AVL selection: ' + str(path))
    if args.paper_source:
        source = args.paper_source.resolve()
        other = {p.relative_to(source).as_posix() for p in source.rglob('*') if p.is_file()}
        check(other == actual, 'Manuscript/repository data file lists differ')
        for relative in sorted(actual & other):
            check(sha(data / relative) == sha(source / relative), 'Manuscript data differs: ' + relative)
    print('Paper package:', len(entries), 'checksummed files;', len(parameters), 'parameters')

    links = 0
    docs = [ROOT / 'README.md', ROOT / 'THIRD_PARTY_NOTICES.md',
            *(ROOT / 'docs').rglob('*.md'), *(ROOT / 'examples').rglob('*.md'),
            ROOT / 'tests/README.md', ROOT / 'paper/README.md', ROOT / 'dev/reference/README.md']
    for doc in docs:
        check(doc.is_file(), 'Missing documentation: ' + str(doc))
        if not doc.is_file():
            continue
        for link in re.findall(r'\]\(([^)]+)\)', doc.read_text(encoding='utf-8')):
            target = unquote(link.split('#')[0]).strip('<>')
            if not target or re.match(r'^[a-zA-Z]+:', target):
                continue
            links += 1
            check((doc.parent / target).exists(), str(doc.relative_to(ROOT)) + ': broken link ' + target)
        for tool in re.findall(r'\btools/[\w.-]+\.(?:py|ps1)', doc.read_text(encoding='utf-8')):
            check((ROOT / tool).is_file(), 'Missing documented tool: ' + tool)
    print('Documentation:', links, 'local links checked')
    if errors:
        raise SystemExit('\n'.join(errors))
    print('PASS package, source baseline, Qt4 forms, paper evidence and local links')


if __name__ == '__main__':
    main()
