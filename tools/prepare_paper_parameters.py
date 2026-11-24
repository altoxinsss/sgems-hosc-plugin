"""Make relocatable run copies of the archived SGeMS parameters, without changing science."""
import argparse
from pathlib import Path
import shutil
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'paper/Data'


def prepare(output):
    output = output.resolve()
    if output.exists():
        raise ValueError('Choose a new output directory; existing files are never overwritten.')
    if not str(output).isascii():
        raise ValueError('The Qt4 host requires an ASCII run path.')
    if output == DATA or DATA in output.parents or output in ROOT.parents:
        raise ValueError('Run copies must not overwrite repository inputs.')
    parameters = sorted((DATA / 'Parameters').rglob('*.par'))
    prepared = []
    for path in parameters:
        tree = ET.parse(path)
        section = path.relative_to(DATA / 'Parameters').parts[0]
        record = tree.getroot().find('audit_output_directory')
        if record is None:
            raise ValueError('Missing Record directory field: ' + str(path))
        record.set('value', (output / 'outputs' / section).as_posix())
        offset = tree.getroot().find('origin_offsets_file')
        if offset is not None and offset.get('value', '').strip():
            offset.set('value', (output / 'protocol/origin_offsets_500.csv').as_posix())
        prepared.append((path, tree, section))
    output.mkdir(parents=True, exist_ok=False)
    (output / 'protocol').mkdir()
    shutil.copy2(DATA / 'Parameters/origin_offsets_500.csv', output / 'protocol')
    for path, tree, section in prepared:
        destination = output / path.relative_to(DATA)
        destination.parent.mkdir(parents=True, exist_ok=True)
        (output / 'outputs' / section).mkdir(parents=True, exist_ok=True)
        tree.write(destination, encoding='utf-8', xml_declaration=True)
    supplementary = Path('Parameters/Supp_weighted_attributes/Supp_samples.dat')
    shutil.copy2(DATA / supplementary, output / supplementary)
    print(f'Prepared {len(prepared)} parameter copies in {output}')
    print('Load inputs before running. Original paper/Data files are unchanged.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-root', required=True, type=Path)
    args = parser.parse_args()
    try:
        prepare(args.output_root)
    except ValueError as error:
        parser.error(str(error))
