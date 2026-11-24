"""Independent membership/formula checks, optionally through the real SGeMS DLL."""
import argparse
import configparser
import csv
import itertools
import json
import math
from pathlib import Path
import subprocess
import sys
import xml.etree.ElementTree as ET


def run(argv, cwd, log):
    with log.open('wb') as stream:
        result = subprocess.run(list(map(str, argv)), cwd=cwd, stdout=stream,
                                stderr=subprocess.STDOUT, timeout=120,
                                creationflags=getattr(subprocess, 'CREATE_NO_WINDOW', 0))
    return result.returncode


def read_values(path):
    return [list(map(float, row.split())) for row in path.read_text().splitlines()[6:] if row.strip()]


def equal(a, b, label, tol=2e-9):
    if math.isnan(a) and math.isnan(b):
        return
    assert math.isfinite(a) and abs(a-b) <= 2e-12 + tol*abs(b), (label, a, b)


def rotate(angles):
    a, b, c = map(math.radians, (90-angles[0], -angles[1], angles[2]))
    ca, sa, cb, sb, cc, sc = math.cos(a), math.sin(a), math.cos(b), math.sin(b), math.cos(c), math.sin(c)
    return [(ca*cb, sa*cb, -sb), (ca*sb*sc-sa*cc, sa*sb*sc+ca*cc, cb*sc),
            (ca*sb*cc+sa*sc, sa*sb*cc-ca*sc, cb*cc)]


def oracle(points, params):
    radii = list(map(float, params['neighborhood_radii'].split(',')))
    axes = rotate(list(map(float, params['neighborhood_angles'].split(','))))
    box = params['neighborhood_shape'] == 'box'
    extent = [sum(abs(radii[k]*axes[k][j]) for k in range(3)) if box else
              math.sqrt(sum((radii[k]*axes[k][j])**2 for k in range(3))) for j in range(3)]
    lo = [min(p[k] for p in points) for k in range(3)]
    hi = [max(p[k] for p in points) for k in range(3)]
    if params['grid_bounds_mode'] == 'manual':
        lo = list(map(float, params['grid_bounds_min'].split(',')))
        hi = list(map(float, params['grid_bounds_max'].split(',')))
    start = list(map(float, params['grid_start'].split(','))) if params.get('grid_start') else lo
    steps = list(map(float, params['grid_spacing'].split(',')))
    coords = []
    for k in range(3):
        if k == 2 and steps[k] == 0:
            coords.append([0.0])
        else:
            n = 1+math.floor((hi[k]-start[k]-2*extent[k])/steps[k]+1e-10)
            assert n > 0
            coords.append([start[k]+extent[k]+i*steps[k] for i in range(n)])
    regular = params.get('center_mode','sample_points') == 'regular_grid'
    centers = [(x, y, z) for z, y, x in itertools.product(coords[2], coords[1], coords[0])] if regular else [p[:3] for p in points]
    members = []
    for center_id, c in enumerate(centers):
        ids = []
        for i, p in enumerate(points):
            if not regular and params.get('includes_center','true') == 'false' and i == center_id:
                continue
            if not math.isfinite(p[3]):
                continue
            v = [sum(axes[k][j]*(p[j]-c[j]) for j in range(3)) for k in range(3)]
            fits = all(abs(v[k]) <= 1e-12 if radii[k] == 0 else
                       abs(v[k])/radii[k] <= 1+1e-12 for k in range(3))
            if fits and (box or sum((v[k]/radii[k])**2 for k in range(3) if radii[k]) <= 1+1e-12):
                ids.append(i)
        members.append(ids)
    return centers, members, tuple(map(len, coords)) if regular else None


def verify(output, points, params, order=5, moment=False):
    rows = read_values(output/'result.txt')
    support = list(csv.DictReader((output/'result.window_support.csv').open()))
    weights = [1.0]*len(points)
    if (output/'result.point_weights.csv').exists():
        weights = [float(r['weight']) for r in csv.DictReader((output/'result.point_weights.csv').open())]
    centers, members, dims = oracle(points, params)
    assert len(rows) == len(centers) == len(support)
    for j, (row, center, ids, s) in enumerate(zip(rows, centers, members, support)):
        for a, b in zip(row[:3], center):
            equal(a, b, 'center')
        assert int(s['sample_count']) == len(ids), (j, s, ids)
        valid = len(ids) >= int(params['min_samples'])
        assert bool(int(s['valid'])) == valid
        wanted = math.nan
        if valid:
            total = sum(weights[i] for i in ids)
            p = [weights[i]/total for i in ids]
            vals = [points[i][3] for i in ids]
            mean = sum(w*v for w, v in zip(p, vals))
            if moment:
                wanted = sum(w*v**order for w, v in zip(p, vals))
            elif order == 1:
                wanted = mean
            else:
                mu = {k: sum(w*(v-mean)**k for w, v in zip(p, vals)) for k in (2, 3, 4, 5)}
                wanted = {2: mu[2], 4: mu[4]-3*mu[2]**2, 5: mu[5]-10*mu[3]*mu[2]}[order]
            ess = total**2/sum(weights[i]**2 for i in ids)
            equal(float(s['effective_sample_size']), ess, 'ESS')
        equal(row[3], wanted, ('statistic', j))
    return rows, dims


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--cli', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--host', type=Path)
    parser.add_argument('--sgems-tests', type=Path)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    if args.host:
        sys.path.insert(0, str(args.sgems_tests))
        from current_host_checks import write_grid, read_grid, f32
    report = []
    for dimensions in (2, 3):
        points = []
        for z, y, x in itertools.product(range(5) if dimensions == 3 else [0], range(7), range(7)):
            if 0 < x < 6 and 0 < y < 6 and (x+2*y+z) % 4 == 0:
                continue
            value = (3*x+2*y+z*z+((x*y)%3))/32
            if (x,y,z) == (3,3,0):
                value = math.nan
            points.append((x-2.25, y+1.125, z+2.5 if dimensions == 3 else 0.0, value))
        # Repeated coordinates are independent input records.
        points.append((*points[5][:3], points[5][3]+0.125))
        fixture = args.out/('points'+str(dimensions)+'.dat')
        fixture.write_text('grid-center fixture\n4\nx\ny\nz\nphi\n'+
                           '\n'.join(' '.join(map(str,p)) for p in points)+'\n')
        if args.host:
            write_grid(fixture.with_suffix('.sgems'), 'points', None, None, None,
                       {'phi':[p[3] for p in points]}, points=[p[:3] for p in points])
        methods = ('uniform','cell_declustering','polygonal_declustering','kriging_declustering')
        for shape, method, evaluator in itertools.product(('box','ellipsoid'), methods, ('AVL Based algorithm','Bell')):
            cases = ['auto','manual','sample_centers']
            if method == 'uniform' and evaluator == 'AVL Based algorithm':
                cases += ['rotated','empty','raw','mean','sample_exclude']
            for case in cases:
                name = f'd{dimensions}_{shape}_{method}_{evaluator.split()[0]}_{case}'
                root = args.out/name
                root.mkdir(exist_ok=True)
                zlo, zhi = (2.5,6.5) if dimensions == 3 else (0,0)
                params = dict(center_mode='regular_grid',neighborhood_shape=shape,
                    neighborhood_radii='1,1,'+('1' if dimensions == 3 else '0'),
                    neighborhood_angles='90,0,0',grid_bounds_mode='auto',
                    grid_spacing='1,1,'+('1' if dimensions == 3 else '0'),
                    min_samples='3',includes_center='true')
                if case in ('manual','empty'):
                    params.update(grid_bounds_mode='manual',
                        grid_bounds_min=f'-4.25,-.875,{zlo}',
                        grid_bounds_max=f'4.75,9.125,{zhi}',
                        grid_start=f'-4.0,-.625,{zlo}')
                if case == 'empty':
                    params.update(grid_bounds_min=f'-10.25,-10.875,{zlo}',
                                  grid_start=f'-10.0,-10.625,{zlo}')
                if case == 'rotated':
                    params['neighborhood_angles'] = '60,15,20' if dimensions == 3 else '60,0,0'
                    params['neighborhood_radii'] = '1,.75,'+('.5' if dimensions == 3 else '0')
                if case.startswith('sample_'):
                    params['center_mode']='sample_points'
                    params['includes_center']='false' if case=='sample_exclude' else 'true'
                decluster = dict(cell_size='2,2'+(',2' if dimensions==3 else ''),
                    shifts_per_axis='2',domain_min='-2.5,1'+(',2.25' if dimensions==3 else ''),
                    domain_max='4,7.5'+(',6.75' if dimensions==3 else ''),
                    kriging_model='exponential',kriging_ranges='3,3'+(',3' if dimensions==3 else ''),
                    kriging_nugget_fraction='.8',kriging_integration_cells='3,3'+(',3' if dimensions==3 else ''))
                order = 1 if case == 'mean' else 5
                target = 'Moment' if case == 'raw' else 'Cumulant'
                config = configparser.ConfigParser()
                config['general'] = dict(source_type='point_set',source_grid=str(fixture),
                    source_property='phi',coordinate_dimensions=str(dimensions),
                    target_property='result',calculate_target=target,calculate_algorithm=evaluator)
                config['sliding_window'] = dict(params,moment_cumulant_order=str(order),
                    center_cumulant_input_for_order_ge_2='true')
                config['weighted_moment'] = dict(weight_mode=method,enabled=str(method!='uniform').lower())
                config['declustering'] = decluster
                config['runtime'] = dict(output_dir=str(root/'output'))
                with (root/'config.ini').open('w') as f:
                    config.write(f)
                assert run([args.cli,root/'config.ini'],root,root/'cli.log') == 0, name
                expected,dims = verify(root/'output',points,params,order,target=='Moment')
                if args.host:
                    audit = root/'audit';audit.mkdir(exist_ok=True)
                    xml = ET.Element('parameters');ET.SubElement(xml,'algorithm',name='HOSC_WINDOW')
                    fields = dict(params,**decluster,Source_Grid='points',Source_Property='phi',props='result',
                        coordinate_dimensions=str(dimensions),moment_order=str(order),Target_Button=target,
                        Algorithm_Button=evaluator,weight_mode=method,output_diagnostics='true',
                        output_point_weights='true',audit_output_directory=str(audit))
                    if target=='Cumulant' and order>=2:
                        fields['output_standardized_cumulant']='true'
                    for key,val in fields.items():
                        ET.SubElement(xml,key,value=val)
                    (root/'parameters.par').write_text(ET.tostring(xml,encoding='unicode'))
                    commands = [f'LoadObjectFromFile {fixture.with_suffix(".sgems")}::s-gems',
                        'RunGeostatAlgorithm HOSC_WINDOW::/GeostatParamUtils/XML::'+ET.tostring(xml,encoding='unicode'),
                        f'SaveGeostatGrid HOSC_WINDOW_result::{root/"saved.sgems"}::s-gems']
                    script = root/'run.cmd';script.write_text('\n'.join(commands)+'\n')
                    assert run([args.host/'sgems.exe',script],args.host,root/'host.log') == 0,name
                    completed=list(audit.rglob('sgems_result.txt'));assert len(completed)==1,(name,(root/'host.log').read_text())
                    actual,_ = verify(completed[0].parent/'output',points,params,order,target=='Moment')
                    saved=read_grid(root/'saved.sgems')
                    if dims:
                        assert saved['kind']=='Cgrid' and saved['dims']==dims,(name,saved)
                        for a,b in zip(saved['origin'],expected[0][:3]):equal(a,b,'host origin',tol=2e-6)
                    else:
                        assert saved['kind']=='Point_set'
                    for i,(a,b) in enumerate(zip(actual,expected)):
                        equal(a[3],b[3],'DLL/CLI')
                        equal(saved['props']['result'][i],f32(b[3]),'stored statistic')
                    assert len(saved['props']['result'])==len(expected)
                    if fields.get('output_standardized_cumulant')=='true':
                        _,members,_=oracle(points,params)
                        weights=[1.0]*len(points)
                        if method!='uniform':
                            weights=[float(r['weight']) for r in csv.DictReader((root/'output/result.point_weights.csv').open())]
                        for i,ids in enumerate(members):
                            wanted=math.nan
                            if math.isfinite(expected[i][3]):
                                total=sum(weights[j] for j in ids)
                                mean=sum(weights[j]*points[j][3] for j in ids)/total
                                variance=sum(weights[j]*(points[j][3]-mean)**2 for j in ids)/total
                                if variance>1e-12:wanted=expected[i][3]/variance**(order/2)
                            equal(saved['props']['result_standardized'][i],f32(wanted),'standardized',tol=2e-6)
                report.append(dict(case=name,centers=len(expected),status='PASS'))
    # Invalid geometry must fail, rather than silently use the legacy centers.
    failures = [
        dict(grid_spacing='0,1,1'),dict(grid_spacing='1,1,nan'),
        dict(grid_bounds_mode='manual',grid_bounds_min='4,4,4',grid_bounds_max='2,2,2'),
        dict(grid_start='-100,0,0'),dict(neighborhood_radii='100,100,100'),
        dict(center_mode='unknown'),dict(neighborhood_shape='unknown'),
        dict(includes_center='false'),
        dict(grid_bounds_mode='manual',grid_bounds_min='',grid_bounds_max='1,1,1')]
    for i,patch in enumerate(failures):
        bad=configparser.ConfigParser()
        bad.read(root/'config.ini')
        bad['sliding_window'].update(center_mode='regular_grid',includes_center='true')
        bad['sliding_window'].update(patch)
        file=args.out/f'invalid_{i}.ini'
        with file.open('w') as f:bad.write(f)
        assert run([args.cli,file],args.out,args.out/f'invalid_{i}.log') != 0,patch
    (args.out/'report.json').write_text(json.dumps(dict(passed=len(report),invalid_checks=len(failures),
        host_checked=bool(args.host),cases=report),indent=2))
    print(f'{len(report)} geometry/weight/evaluator cases and {len(failures)} invalid configurations passed.')


if __name__ == '__main__':
    main()
