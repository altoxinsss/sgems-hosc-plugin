"""Window weighting integration checks against explicit k1--k5 identities.

Spatial weight algorithms already have independent fixtures. Here a clipped
Voronoi strip fixture additionally fixes the expected point masses exactly.
"""
import argparse
import configparser
import csv
import itertools
import json
import math
from pathlib import Path
import subprocess


def records(path):
    with path.open(newline='') as stream:
        return list(csv.DictReader(stream))


def equal(actual, expected, label, atol=2e-10):
    actual, expected = float(actual), float(expected)
    if math.isnan(actual) and math.isnan(expected):
        return
    assert math.isfinite(actual) and math.isfinite(expected), (label, actual, expected)
    assert abs(actual-expected) <= atol + 2e-9*abs(expected), (label, actual, expected)


def statistic(values, weights, order, moment=False):
    total = math.fsum(weights)
    mean = math.fsum(w*v for w,v in zip(weights, values))/total
    if moment:
        return math.fsum(w*v**order for w,v in zip(weights, values))/total
    if order == 1:
        return mean
    mu = {n: math.fsum(w*(v-mean)**n for w,v in zip(weights, values))/total for n in range(2, order+1)}
    if order <= 3:
        return mu[order]
    if order == 4:
        return mu[4] - 3*mu[2]**2
    if order == 5:
        return mu[5] - 10*mu[3]*mu[2]
    raise ValueError('This independent fixture covers orders 1 to 5.')


def verify_output(output, expected_support, point_weights, values, order, moment=False):
    rows = [list(map(float, line.split())) for line in (output/'result.txt').read_text().splitlines()[6:] if line.strip()]
    support = records(output/'result.window_support.csv')
    assert len(rows) == len(support) == len(expected_support)
    nonuniform = False
    for i, (ids, valid) in enumerate(expected_support):
        ws = [point_weights[j] for j in ids]
        vs = [values[j] for j in ids]
        expected = statistic(vs, ws, order, moment) if valid else math.nan
        # The Cartesian CLI omits invalid result rows; fixtures passed here
        # retain all point centers or use complete Cartesian windows.
        equal(rows[i][3], expected, ('value', i))
        equal(support[i]['sample_count'], len(ids), ('count', i))
        equal(support[i]['valid'], int(valid), ('valid', i))
        total, squares = math.fsum(ws), math.fsum(w*w for w in ws)
        equal(support[i]['weight_sum'], total, ('sum', i))
        equal(support[i]['weight_squared_sum'], squares, ('squares', i))
        ess = total*total/squares if squares else 0
        equal(support[i]['effective_sample_size'], ess, ('ess', i))
        if valid and len(ws)>1:
            nonuniform |= ess < len(ws) - 1e-5
    return nonuniform


def run(cli, root, name, source, dimensions, settings, method, order, bell=False,
        moment=False, point=True, extra='', general_extra=''):
    case = root/name
    case.mkdir(parents=True, exist_ok=True)
    config = case/'config.ini'
    config.write_text(f'''[general]
source_type={'point_set' if point else 'cartesian'}
source_grid={source}
source_property=value
coordinate_dimensions={dimensions}
target_property=result
calculate_target={'Moment' if moment else 'Cumulant'}
calculate_algorithm={'Bell' if bell else 'AVL Based algorithm'}
{general_extra}
[sliding_window]
moment_cumulant_order={order}
center_cumulant_input_for_order_ge_2=true
{settings}
[weighted_moment]
weight_mode={method}
normalize_weights=true
[declustering]
{extra}
[runtime]
output_dir={case/'output'}
''')
    process = subprocess.run([str(cli), str(config)], capture_output=True, text=True, timeout=60)
    (case/'run.log').write_text(process.stdout+process.stderr)
    assert process.returncode == 0, (name, process.stdout, process.stderr)
    return case/'output'


def write_source(path, points, values):
    path.write_text('window fixture\n4\nx\ny\nz\nvalue\n'+''.join(
        ' '.join(map(str, (*p,v)))+'\n' for p,v in zip(points,values)))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--cli', required=True, type=Path)
    parser.add_argument('--out', required=True, type=Path)
    a = parser.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    methods = ['cell_declustering','polygonal_declustering','kriging_declustering']
    report = []
    for dims in (2,3):
        points = [(x,.5,z) for z in ([0] if dims==2 else [.25,.75]) for x in (.5,.75,.75,2,3.5,2.5)]
        values = [(.2,.3,.9,.4,.8,math.nan)[i%6]+.05*(i//6) for i in range(len(points))]
        source = a.out/f'points{dims}.gslib'
        write_source(source,points,values)
        suffix = '' if dims==2 else ',1'
        params = f'''cell_size=1,1{suffix}
origin=0,0{',0' if dims==3 else ''}
shifts_per_axis=1
domain_min=0,0{',0' if dims==3 else ''}
domain_max=4,1{suffix}
kriging_model=exponential
kriging_ranges=2,1{suffix}
kriging_nugget_fraction=.8
kriging_integration_cells=12,4{',4' if dims==3 else ''}'''
        for method in methods:
            for excludes in (False,True):
                radii = (1.5,1,0 if dims==2 else 1)
                settings = f"neighborhood_radii={','.join(map(str,radii))}\nmin_samples=2\nincludes_center={'false' if excludes else 'true'}"
                expected = []
                for i,p in enumerate(points):
                    ids = [j for j,q in enumerate(points) if math.isfinite(values[j]) and (not excludes or i!=j)
                           and sum(((p[k]-q[k])/radii[k])**2 for k in range(dims)) <= 1+1e-12]
                    expected.append((ids,len(ids)>=2))
                for order, moment in [(1,False),(2,False),(3,False),(4,False),(5,False),(4,True)]:
                    for bell in (False,True):
                        name = f'p{dims}_{method}_{excludes}_{order}_{moment}_{bell}'
                        output = run(a.cli,a.out,name,source,dims,settings,method,order,bell,moment,extra=params)
                        ws = [float(r['weight']) for r in records(output/'result.point_weights.csv')]
                        assert all(ws[i]==0 for i,v in enumerate(values) if not math.isfinite(v))
                        if method == 'polygonal_declustering' and dims==2:
                            for got,want in zip(ws,[.625/4,.75/8,.75/8,1.375/4,1.25/4,0]):
                                equal(got,want,'exact clipped Voronoi strips')
                        assert verify_output(output,expected,ws,values,order,moment), name
                        meta=(output/'result.meta.txt').read_text()
                        assert ('evaluator=univariate_recurrence' if bell else 'evaluator=AVL') in meta or moment
                        report.append(name)
    # Physical grid coordinates must drive declustering; lattice coordinates
    # remain responsible for window membership.
    for dims in (2,3):
        nz=1 if dims==2 else 2
        lattice=[(x,y,z) for z in range(nz) for y in range(3) for x in range(4)]
        values=[.2+.03*x+.04*y+.02*z+.05*((x*y+z)%3) for x,y,z in lattice]
        source=a.out/f'grid{dims}.gslib';write_source(source,lattice,values)
        params='''cell_size=3,4{size}
origin=10,-2{zmin}
shifts_per_axis=2
domain_min=9.9,-2.2{zmin}
domain_max=17.7,5{zmax}
kriging_model=exponential
kriging_ranges=3,4{size}
kriging_nugget_fraction=.8
kriging_integration_cells=8,6{cells}'''.format(size='' if dims==2 else ',1',zmin='' if dims==2 else ',.9',zmax='' if dims==2 else ',1.9',cells='' if dims==2 else ',4')
        settings=f'window_size=3,2,{nz}\nmove_step=1,1,{0 if dims==2 else 1}'
        expected=[([i for i,(x,y,z) in enumerate(lattice) if ax<=x<ax+3 and ay<=y<ay+2],True)
                  for ay in range(2) for ax in range(2)]
        for method, bell, order in itertools.product(methods,(False,True),(1,4,5)):
            name=f'g{dims}_{method}_{bell}_{order}'
            output=run(a.cli,a.out,name,source,3,settings,method,order,bell,point=False,extra=params,
                       general_extra='source_origin=10,-2,1\nsource_spacing=2,3,.5')
            rows=records(output/'result.point_weights.csv');ws=[float(r['weight']) for r in rows]
            for row,(x,y,z) in zip(rows,lattice):
                for key,value in zip(('x','y','z'),(10+2*x,-2+3*y,1+.5*z)):
                    equal(row[key],value,'physical weight coordinates')
            assert verify_output(output,expected,ws,values,order),name
            report.append(name)
    # Same cumulants under a large shift, with all window weights held fixed.
    pts=[(.5,.5,0),(.75,.5,0),(2,.5,0),(3.5,.5,0)]
    shift_values=[.25,.5,.75,1.5]
    shifted=a.out/'shifted.gslib';write_source(shifted,pts,[v+1024 for v in shift_values])
    for bell in (False,True):
        out=run(a.cli,a.out,'shift_'+str(bell),shifted,2,'neighborhood_radii=10,10,0','cell_declustering',5,bell,
                extra='cell_size=1,1\norigin=0,0\nshifts_per_axis=1')
        ws=[float(r['weight']) for r in records(out/'result.point_weights.csv')]
        verify_output(out,[(list(range(4)),True)]*4,ws,shift_values,5)
        report.append('shift_'+str(bell))
    # Unsupported/ambiguous settings must fail, rather than yield uniform output.
    base=(a.out/report[0]/'config.ini').read_text()
    for i,(old,new) in enumerate([('normalize_weights=true','normalize_weights=false'),
            ('weight_mode=cell_declustering','weight_mode=cell_declustering\nenabled=false'),
            ('weight_mode=cell_declustering','weight_mode=cell_declustering\nuse_weighted_moment=false'),
            ('weight_mode=cell_declustering','weight_mode=cell_declustering\nzero_first_order_moment=true'),
            ('cell_size=1,1','cell_size=0,1')]):
        config=a.out/f'reject_{i}.ini';config.write_text(base.replace(old,new))
        p=subprocess.run([str(a.cli),str(config)],capture_output=True,text=True,timeout=60)
        assert p.returncode!=0,(i,p.stdout,p.stderr)
    (a.out/'report.json').write_text(json.dumps(dict(passed=len(report),cases=report,rejected_invalid_configs=5),indent=2))
    print(f'{len(report)} weighted window configurations and 5 invalid configurations passed.')


if __name__=='__main__':
    main()
