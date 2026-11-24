"""Run the six introductory examples in an isolated x64 SGeMS host and verify cumulants."""
import argparse
import csv
import itertools
import json
import math
from pathlib import Path
import sys
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'tests'))
from current_host_checks import read_grid, close, run_process


def partitions(items):
    if not items:
        yield []
        return
    for groups in partitions(items[1:]):
        yield [[items[0]]] + groups
        for i in range(len(groups)):
            yield groups[:i] + [[items[0]] + groups[i]] + groups[i+1:]


def joint_cumulant(events, weights):
    total = sum(weights)
    def moment(columns):
        return sum(w * math.prod(event[j] for j in columns) for event, w in zip(events, weights)) / total
    return sum((-1)**(len(groups)-1) * math.factorial(len(groups)-1) *
               math.prod(moment(group) for group in groups) for groups in partitions(list(range(4))))


def expected_values(source, template, external):
    values = source['props']['value']
    weights = source['props']['weight'] if external else [1.0] * len(values)
    if not template:
        return [joint_cumulant([[v]*4 for v in values], weights)]
    if source['kind'] == 'Point_set':
        points = source['points']
    else:
        dims, origin, spacing = source['dims'], source['origin'], source['spacing']
        points = [(origin[0]+x*spacing[0], origin[1]+y*spacing[1], origin[2]+z*spacing[2])
                  for z, y, x in itertools.product(range(dims[2]), range(dims[1]), range(dims[0]))]
    result = []
    for lag in (0, 1):
        events, event_weights = [], []
        for i, point in enumerate(points):
            matches = [i] if lag == 0 else [j for j, other in enumerate(points)
                        if all(abs(a-b) < 1e-6 for a, b in zip(other, (point[0]+lag, point[1], point[2])))]
            for j in matches:
                events.append([values[i], values[j], values[j], values[j]])
                event_weights.append(weights[i] if i == j else math.sqrt(weights[i]*weights[j]))
        result.append(joint_cumulant(events, event_weights))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--host', required=True, type=Path)
    parser.add_argument('--out', required=True, type=Path)
    args = parser.parse_args()
    host, out = args.host.resolve(), args.out.resolve()
    if not (host / 'HOSC_ISOLATED_TEST_HOST.txt').is_file():
        parser.error('Use tools/stage_test_host.py first; the normal user host is not a test target.')
    if out.exists() or out == host or host in out.parents or not str(out).isascii():
        parser.error('--out must be a new ASCII directory outside the isolated host')
    out.mkdir(parents=True)
    report = []
    for path in sorted((ROOT / 'examples/Parameters').glob('*.par')):
        name = path.stem
        case = out / name
        case.mkdir()
        audit = case / 'records'
        audit.mkdir()
        node = ET.parse(path).getroot()
        # SGeMS command files require the entire XML payload on one line.
        for element in node.iter():
            if element.text is not None and not element.text.strip():
                element.text = None
            if element.tail is not None and not element.tail.strip():
                element.tail = None
        ET.SubElement(node, 'audit_output_directory', value=str(audit))
        algorithm = node.find('algorithm').get('name')
        source_name = node.find('Source_Grid').get('value')
        fixture = ROOT / 'examples/Inputs.prj' / (source_name + '.sgems')
        saved = case / 'result.sgems'
        command = case / 'run.cmd'
        command.write_text('\n'.join([
            f'LoadObjectFromFile {fixture}::sgems',
            f'RunGeostatAlgorithm {algorithm}::/GeostatParamUtils/XML::{ET.tostring(node, encoding="unicode")}',
            f'SaveGeostatGrid {algorithm}_{name}::{saved}::sgems',
        ]) + '\n', encoding='ascii')
        run_process([host / 'sgems-x64.exe', command], host, case / 'host.log')
        actual = read_grid(saved)['props'][name]
        reference = expected_values(read_grid(fixture), algorithm == 'HOSC_TEMPLATE', 'external' in name)
        if len(reference) == 1:
            reference *= len(actual)
        assert len(actual) == len(reference), (name, 'output size')
        for i, (a, b) in enumerate(zip(actual, reference)):
            close(a, b, (name, i), atol=2e-7)
        files = list(audit.rglob('performance.csv'))
        assert len(files) == 1, (name, 'performance record')
        with files[0].open(newline='') as stream:
            records = list(csv.DictReader(stream))
        assert len(records) == 1
        perf = records[0]
        assert perf['status'] == 'success' and perf['measurement_status'] == 'ok', perf
        assert float(perf['process_bits']) == 64 and float(perf['elapsed_seconds']) > 0, perf
        report.append(dict(case=name, outputs=len(actual), status='passed',
                           max_absolute_error=max(abs(a-b) for a, b in zip(actual, reference))))
        print(name + ': native x64 result and performance record passed', flush=True)
    (out / 'report.json').write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
    print(f'PASS {len(report)} examples; results in {out}')


if __name__ == '__main__':
    main()
