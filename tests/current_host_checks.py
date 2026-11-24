"""Real SGeMS-process integration checks; runs only in the supplied isolated host.

Fixtures are deliberately nonzero-mean, fractional/anisotropic and include
missing values. Saved native SGeMS objects are checked after each run.
"""
import argparse
import configparser
import csv
import io
import itertools
import json
import math
import os
from pathlib import Path
import struct
import subprocess
import xml.etree.ElementTree as ET

NODATA = -9966699.0


def f32(v):
    return struct.unpack('>f', struct.pack('>f', v))[0]


def string(s):
    b = s.encode('ascii') + b'\0'
    return struct.pack('>I', len(b)) + b


def write_grid(path, name, dims, origin, spacing, props, points=None):
    b = struct.pack('>I', 0xB211175D) + string('Point_set' if points is not None else 'Cgrid') + string(name) + struct.pack('>I', 100)
    if points is None:
        b += struct.pack('>3I6f', *dims, *spacing, *origin)
    else:
        b += struct.pack('>I', len(points))
    b += struct.pack('>I', len(props)) + b''.join(string(p) for p in props)
    if points is not None:
        b += b''.join(struct.pack('>3f', *p) for p in points)
    for values in props.values():
        b += b''.join(struct.pack('>f', v if math.isfinite(v) else NODATA) for v in values)
    path.write_bytes(b)


def read_grid(path):
    data = io.BytesIO(path.read_bytes())
    def get(fmt):
        return struct.unpack('>' + fmt, data.read(struct.calcsize('>' + fmt)))
    def st():
        return data.read(get('I')[0]).rstrip(b'\0').decode('ascii')
    assert get('I')[0] == 0xB211175D
    kind, name, version = st(), st(), get('I')[0]
    assert version >= 100
    out = dict(kind=kind, name=name)
    if kind == 'Cgrid':
        out['dims'] = get('3I')
        out['spacing'], out['origin'] = get('3f'), get('3f')
        n = math.prod(out['dims'])
    else:
        assert kind == 'Point_set'
        n = get('I')[0]
    names = [st() for _ in range(get('I')[0])]
    if kind == 'Point_set':
        out['points'] = [get('3f') for _ in range(n)]
    out['props'] = {key: [math.nan if v == NODATA else v for v in get(str(n) + 'f')] for key in names}
    return out


def close(a, b, label, rtol=2e-6, atol=1e-8):
    if math.isnan(a) and math.isnan(b):
        return
    assert math.isfinite(a) and math.isfinite(b) and abs(a - b) <= atol + rtol * abs(b), (label, a, b)


def cumulant(values, order):
    mu = sum(values) / len(values)
    if order == 1:
        return mu
    moments = [1.0, 0.0] + [sum((v - mu) ** r for v in values) / len(values) for r in range(2, order + 1)]
    ks = [0.0] * (order + 1)
    for r in range(2, order + 1):
        ks[r] = moments[r] - sum(math.comb(r - 1, j - 1) * ks[j] * moments[r - j] for j in range(1, r))
    return ks[order]


def run_process(argv, cwd, logfile):
    argv = list(argv)
    if Path(argv[0]).name.lower() == 'sgems.exe' and (Path(cwd)/'sgems-x64.exe').exists():
        argv[0] = Path(cwd)/'sgems-x64.exe'
        # The Qt4 Beta renamed the legacy binary filter to sgems.
        command = Path(argv[1])
        command.write_text(command.read_text(encoding='ascii').replace('::s-gems','::sgems'),encoding='ascii')
    env = os.environ.copy()
    env['GSTLAPPLIHOME'] = str(cwd)
    env['PATH'] = str(cwd)+';'+env['PATH']
    with logfile.open('wb') as log:
        result = subprocess.run([str(v) for v in argv], cwd=cwd, stdout=log, stderr=subprocess.STDOUT,
                                timeout=90, env=env, creationflags=subprocess.CREATE_NO_WINDOW)
    assert result.returncode == 0, (argv, result.returncode, logfile.read_text(errors='replace')[-2500:])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--host', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--case', default='all')
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    dims, origin, spacing = (5, 4, 3), (10.25, -3.5, 100.0), (2.5, 1.25, 0.5)
    values = [f32(1.3 + .07*x + .12*y + .2*z + .03*((x*y+z)%3)) for z in range(3) for y in range(4) for x in range(5)]
    missing = list(values); missing[0] = math.nan
    write_grid(args.out/'grid.sgems', 'grid', dims, origin, spacing, dict(phi=values, missing=missing))
    # Coincident records, fractional coordinates, nonzero Z and an uninformed center.
    points = [(f32(x), f32(y), f32(z)) for z in (2.0, 2.5) for y in (0.0, .6, 1.4) for x in (0.0, .4, 1.3)]
    points += [points[0]]
    pvals = [f32(1.2 + .1*i + .03*(i%3)**2) for i in range(len(points))]
    pvals[2] = math.nan
    write_grid(args.out/'points.sgems', 'points', None, None, None, dict(phi=pvals), points)
    points2 = [(x, y, 0.0) for x, y, z in points[:9]]
    write_grid(args.out/'points2.sgems', 'points2', None, None, None, dict(phi=pvals[:9]), points2)
    base = dict(Source_Grid='grid', Source_Property='phi', Target_Button='Cumulant', Algorithm_Button='AVL Based algorithm',
                moment_order='4', windowl='2', windoww='2', windowh='2', stepx='1', stepy='1', stepz='1',
                output_diagnostics='true', output_standardized_cumulant='true')
    template = dict(Source_Grid='grid', Source_Property='phi', Target_Button='Cumulant', Algorithm_Button='AVL Based algorithm',
                    NOD='1,1,1', NLAG='2,3,2', LAG='1,1,1', LTOL='0,0,0', AZM='0,90,0', ATOL='0,0,0',
                    BANDWH='0,0,0', DIP='0,0,90', DTOL='0,0,0', BANDWD='0,0,0', fst_dir='2', snd_dir='1', trd_dir='3',
                    cross_pos='', template_sampling='exhaustive_integer_offsets', integer_direction_vectors='0,1,0;1,0,0;0,0,1',
                    xmn_='3.25', ymn_='-2.5', zmn_='10', output_diagnostics='true', weight_mode='uniform')
    cases = []
    def add(name, algorithm, params, **extra):
        cases.append((name, algorithm, dict(params, **extra)))
    add('window_grid', 'HOSC_WINDOW', base)
    add('window_bell', 'HOSC_WINDOW', base, Algorithm_Button='Bell')
    add('window_missing', 'HOSC_WINDOW', base, Source_Property='missing')
    add('window_display', 'HOSC_WINDOW', base, invert_window_output='true', enable_percentile_clipping='true')
    add('window_points', 'HOSC_WINDOW', base, Source_Grid='points', neighborhood_radii='1,1,1', neighborhood_angles='90,0,0', min_samples='3')
    add('window_points2', 'HOSC_WINDOW', base, Source_Grid='points2', coordinate_dimensions='2', neighborhood_radii='1,1,1', neighborhood_angles='90,0,0', min_samples='3', includes_center='false')
    add('template_exact', 'HOSC_TEMPLATE', template)
    add('template_bell', 'HOSC_TEMPLATE', template, Algorithm_Button='Bell')
    add('template_moment', 'HOSC_TEMPLATE', template, Target_Button='Moment')
    add('template_legacy', 'HOSC_TEMPLATE', template, template_sampling='legacy_search')
    add('template_legacy_moment', 'HOSC_TEMPLATE', template, template_sampling='legacy_search', Target_Button='Moment')
    add('template_missing', 'HOSC_TEMPLATE', template, Source_Property='missing')
    add('template_slice', 'HOSC_TEMPLATE', template, Output_3D_Slice_Enable='Enable', Output_3D_Slice_Plane='xz', Output_3D_Slice_Position='2')
    for method in ['cell_declustering', 'polygonal_declustering', 'kriging_declustering']:
        common = dict(weight_mode=method, cell_size='4,2,1', declustering_origin='10,-4,99.9', shifts_per_axis='2',
                      domain_min='10,-4,99.9', domain_max='22,1,101.1', kriging_model='exponential', kriging_ranges='5,2,1',
                      kriging_nugget_fraction='0.8', kriging_integration_cells='4,3,2')
        add('grid_'+method, 'HOSC_TEMPLATE', template, **common)
    ptemplate = dict(template, Source_Grid='points', coordinate_dimensions='3', template_sampling='point_set_tolerance',
                     NOD='2,2', NLAG='2,2', LAG='.6,.5', LTOL='.5,.3', AZM='90,0', ATOL='60,60', BANDWH='1,1',
                     DIP='0,90', DTOL='60,60', BANDWD='1,1', fst_dir='2', snd_dir='1', trd_dir='0', cross_pos='')
    add('template_points', 'HOSC_TEMPLATE', ptemplate)
    add('template_points_bell', 'HOSC_TEMPLATE', ptemplate, Algorithm_Button='Bell')
    add('template_points_negative', 'HOSC_TEMPLATE', ptemplate, LAG='-.6,.5')
    add('template_grid_physical', 'HOSC_TEMPLATE', ptemplate, Source_Grid='grid', LAG='2.5,.5', LTOL='0,0', ATOL='0,0', BANDWH='0,0', DTOL='0,0', BANDWD='0,0')
    for method in ['cell_declustering', 'polygonal_declustering', 'kriging_declustering']:
        add('points_'+method, 'HOSC_TEMPLATE', ptemplate, weight_mode=method, cell_size='1,1,.5', shifts_per_axis='2',
            domain_min='-.1,-.1,1.9', domain_max='1.5,1.6,2.6', kriging_model='exponential', kriging_ranges='1,1,.5',
            kriging_nugget_fraction='.8', kriging_integration_cells='4,4,3')
        add('points2_'+method, 'HOSC_TEMPLATE', ptemplate, Source_Grid='points2', coordinate_dimensions='2',
            DIP='0,0', weight_mode=method, cell_size='1,1', shifts_per_axis='2',
            domain_min='-.1,-.1', domain_max='1.5,1.6', kriging_model='exponential', kriging_ranges='1,1',
            kriging_nugget_fraction='.8', kriging_integration_cells='4,4')
    results = []
    for name, algorithm, params in cases:
        if args.case != 'all' and args.case != name:
            continue
        case_dir = args.out/name; case_dir.mkdir(exist_ok=True)
        audit = case_dir/'audit'; audit.mkdir(exist_ok=True)
        params = dict(params, props=name, audit_output_directory=str(audit))
        xml = ET.Element('parameters'); ET.SubElement(xml, 'algorithm', name=algorithm)
        for key, val in params.items(): ET.SubElement(xml, key, value=val)
        grid = params['Source_Grid']
        output = case_dir/'saved.sgems'
        commands = [f'LoadObjectFromFile {args.out/(grid+".sgems")}::s-gems',
                    f'RunGeostatAlgorithm {algorithm}::/GeostatParamUtils/XML::{ET.tostring(xml, encoding="unicode")}',
                    f'SaveGeostatGrid {algorithm}_{name}::{output}::s-gems']
        if name == 'template_slice':
            commands.append(f'SaveGeostatGrid HOSC_slice_{name}::{case_dir/"slice.sgems"}::s-gems')
        if params.get('weight_mode','uniform') != 'uniform':
            commands.append(f'SaveGeostatGrid {algorithm}_{name}_point_weights::{case_dir/"weights.sgems"}::s-gems')
        command = case_dir/'run.cmd'; command.write_text('\n'.join(commands)+'\n', encoding='ascii')
        run_process([args.host/'sgems.exe', command], args.host, case_dir/'host.log')
        completed = [p for p in audit.iterdir() if (p/'sgems_result.txt').exists()]
        assert len(completed) == 1, (name, (case_dir/'host.log').read_text(errors='replace')[-2000:])
        run = completed[0]; saved = read_grid(output)
        cfg = configparser.ConfigParser(); cfg.read(run/'config.ini')
        txt = [list(map(float, line.split())) for line in (run/'output/result.txt').read_text().splitlines()[6:] if line.strip()]
        sprops = saved['props']
        if saved['kind'] == 'Point_set':
            expected_pts = points2 if grid == 'points2' else points
            assert saved['points'] == expected_pts
            expected_values = pvals[:9] if grid == 'points2' else pvals
            for i, p in enumerate(expected_pts):
                support = [v for j,(q,v) in enumerate(zip(expected_pts,expected_values)) if math.isfinite(v)
                           and (params.get('includes_center')!='false' or i != j)
                           and sum((a-b)**2 for a,b in zip(p,q)) <= 1+1e-12]
                expected = cumulant(support,4) if len(support)>=3 else math.nan
                close(sprops[name][i], expected, (name,i),atol=2e-7)
                close(sprops[name+'_support'][i],len(support),(name,'support',i))
                if len(support)>=3:
                    var=cumulant(support,2)
                    close(sprops[name+'_standardized'][i],expected/var**2 if var>1e-12 else math.nan,(name,'std',i))
        else:
            sd=saved['dims']; so=saved['origin']; ss=saved['spacing']
            if algorithm=='HOSC_WINDOW':
                assert sd == (4,3,2)
                for k in range(3): close(so[k],origin[k]+spacing[k]/2,(name,'origin',k))
                src=missing if params['Source_Property']=='missing' else values
                for z,y,x in itertools.product(range(2),range(3),range(4)):
                    node=(z*3+y)*4+x
                    vals=[src[((z+dz)*4+y+dy)*5+x+dx] for dz,dy,dx in itertools.product(range(2),repeat=3)]
                    expected=cumulant(vals,4) if all(map(math.isfinite,vals)) else math.nan
                    if name!='window_display': close(sprops[name][node],expected,(name,node),atol=2e-8)
                    var=cumulant(vals,2)
                    close(sprops[name+'_standardized'][node],expected/var**2 if var>1e-12 else math.nan,(name,'std',node))
            else:
                nlag=list(map(int,params['NLAG'].split(','))); lag=list(map(float,params['LAG'].split(',')))
                display=[int(params[k])-1 for k in ['fst_dir','snd_dir','trd_dir'] if int(params[k])]
                assert sd == tuple([nlag[k] for k in display]+[1]*(3-len(display)))
                for k in range(3):
                    expected=float(params[['xmn_','ymn_','zmn_'][k]])
                    if k<len(display) and lag[display[k]]<0: expected+=(sd[k]-1)*lag[display[k]]
                    close(so[k],expected,(name,'origin',k))
                    close(ss[k],abs(lag[display[k]]) if k<len(display) else 1,(name,'spacing',k))
                if name=='template_slice':
                    sl=read_grid(case_dir/'slice.sgems'); assert sl['dims']==(3,1,2)
                    for z,x in itertools.product(range(2),range(3)):close(sl['props'][name][z*3+x],sprops[name][(z*2+1)*3+x],(name,'slice',x,z))
        # Independently verify the mapping of every core row into native SGeMS nodes.
        for i,row in enumerate(txt):
            if saved['kind']=='Point_set': node=i
            else:
                if algorithm=='HOSC_WINDOW': coords=[round(v) for v in row[:3]]
                else:
                    base_origin=saved['origin'] if cfg['general']['source_type']=='point_set' else [0]*3
                    coords=[round((row[k]-base_origin[k])/saved['spacing'][k]) for k in range(3)]
                node=(coords[2]*saved['dims'][1]+coords[1])*saved['dims'][0]+coords[0]
            close(sprops[name][node],row[3],(name,'import',i))
        if params.get('weight_mode','uniform')!='uniform':
            w=read_grid(case_dir/'weights.sgems')
            with (run/'output/result.point_weights.csv').open() as f: rows=list(csv.DictReader(f))
            assert len(rows)==len(w['props']['declustering_weight'])
            for i,r in enumerate(rows): close(w['props']['declustering_weight'][i],float(r['weight']),(name,'weight',i))
        # At zero lag, original fourth raw moment must survive nonzero mean.
        if name in ('template_moment','template_legacy_moment'):
            close(sprops[name][0],sum(v**4 for v in values)/len(values),(name,'raw fourth moment'),rtol=3e-6)
        results.append(dict(case=name, status='passed', run=str(run), nodes=len(sprops[name])))
        print(name+': passed', flush=True)
    (args.out/'report.json').write_text(json.dumps(results,indent=2)+'\n')
    if args.case=='all':
        for left,right in [('window_grid','window_bell'),('template_exact','template_bell'),('template_points','template_points_bell')]:
            a=read_grid(args.out/left/'saved.sgems')['props'][left]
            b=read_grid(args.out/right/'saved.sgems')['props'][right]
            for i,(x,y) in enumerate(zip(a,b)): close(x,y,(left,right,i),rtol=2e-6,atol=2e-8)
        # An independent raw-product / joint-cumulant calculation at every
        # exact Cartesian lag, including a per-variable partition expansion.
        def partitions(items):
            if not items:
                yield []; return
            for parts in partitions(items[1:]):
                yield [[items[0]]]+parts
                for k in range(len(parts)):
                    yield parts[:k]+[[items[0]]+parts[k]]+parts[k+1:]
        for check,target in [('template_exact','Cumulant'),('template_moment','Moment')]:
            actual=read_grid(args.out/check/'saved.sgems')['props'][check]
            for lz,ly,lx in itertools.product(range(2),range(2),range(3)):
                tuples=[]
                for z,y,x in itertools.product(range(3-lz),range(4-ly),range(5-lx)):
                    ids=[(z*4+y)*5+x,(z*4+y+ly)*5+x,(z*4+y)*5+x+lx,((z+lz)*4+y)*5+x]
                    tuples.append([values[i] for i in ids])
                def moment(cols):return sum(math.prod(t[k] for k in cols) for t in tuples)/len(tuples)
                expected=moment(range(4)) if target=='Moment' else sum((-1)**(len(p)-1)*math.factorial(len(p)-1)*math.prod(moment(block) for block in p) for p in partitions(list(range(4))))
                close(actual[(lz*2+ly)*3+lx],expected,(check,'independent partitions',lx,ly,lz),atol=2e-8)
        print('Bell/AVL and independent partition references: passed',flush=True)


if __name__=='__main__': main()
