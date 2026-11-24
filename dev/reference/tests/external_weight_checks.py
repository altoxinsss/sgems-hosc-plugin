"""End-to-end external observation weights: independent identities and round trips."""
import argparse
import csv
import itertools
import json
import math
from pathlib import Path
import subprocess


def close(a, b, label):
    assert len(a) == len(b), (label, len(a), len(b))
    for x, y in zip(a, b):
        assert (math.isnan(x) and math.isnan(y)) or math.isclose(x, y, rel_tol=2e-8, abs_tol=2e-9), (label, x, y)


def source(path, points, values, weights, label='weight'):
    path.write_text('external-weight fixture\n5\nx\ny\nz\nvalue\n'+label+'\n'+''.join(
        ' '.join(map(str, (*p,v,w)))+'\n' for p,v,w in zip(points,values,weights)))


def partitions(n):
    if n == 0:
        yield []
        return
    for blocks in partitions(n-1):
        yield blocks + [[n-1]]
        for k in range(len(blocks)):
            yield [b+[n-1] if j==k else b[:] for j,b in enumerate(blocks)]


def joint(rows, weights, order, moment=False):
    total = math.fsum(weights)
    def m(ids):
        return math.fsum(w*math.prod(row[j] for j in ids) for row,w in zip(rows,weights))/total
    if moment:
        return m(range(order))
    return math.fsum((-1)**(len(p)-1)*math.factorial(len(p)-1)*math.prod(m(b) for b in p)
                     for p in partitions(order))


def output_values(out):
    lines=(out/'result.txt').read_text().splitlines()
    return [float(line.split()[-1]) for line in lines[2+int(lines[1].split()[-1]):] if line.strip()]


def run(exe, root, name, path, dims, point, module, order=4, mode='external',
        prop='weight', extra='', bell=False, moment=False, fail=False):
    case=root/name
    case.mkdir(parents=True,exist_ok=True)
    if module=='template':
        geometry=f'''[template_shape]
template_sampling={'point_set_tolerance' if point else 'exhaustive_integer_offsets'}
evaluate_raw_complete_tuples=true
NOD={order-1}
NLAG=2
LAG=1
LTOL=0
AZM=0
ATOL=0
BANDWH=0
DIP=0
DTOL=0
BANDWD=0
integer_direction_vectors={'1,0,0' if dims==3 else '1,0'}
minimum_template_count=1
[target_cross_sections]
displayed_dirs=1
fixed_positions=
'''
    else:
        geometry=f'''[sliding_window]
moment_cumulant_order={order}
center_cumulant_input_for_order_ge_2=true
window_size=3,2,{2 if dims==3 else 1}
move_step=1,1,{1 if dims==3 else 0}
require_full_window=true
require_complete_window=false
skip_missing_samples=true
minimum_samples=1
neighborhood_shape=box
neighborhood_radii=100,100,100
neighborhood_angles=90,0,0
center_mode=sample_points
min_samples=1
'''
    if module=='window':
        grid_keys=('window_size','move_step','require_full_window','require_complete_window','skip_missing_samples','minimum_samples')
        point_keys=('neighborhood_shape','neighborhood_radii','neighborhood_angles','center_mode','min_samples')
        geometry='\n'.join(line for line in geometry.splitlines() if line.split('=')[0] not in (grid_keys if point else point_keys))+'\n'
    config=f'''[general]
source_type={'point_set' if point else 'cartesian'}
source_grid={path}
source_property=value
coordinate_dimensions={dims if point else 3}
target_property=result
calculate_target={'Moment' if moment else 'Cumulant'}
calculate_algorithm={'Bell' if bell else 'AVL Based algorithm'}
{geometry}
[weighted_moment]
weight_mode={mode}
weight_property={prop}
normalize_weights=true
zero_first_order_moment=false
{extra}
[runtime]
output_dir={case/'output'}
write_template_tuples=true
'''
    (case/'config.ini').write_text(config)
    p=subprocess.run([str(exe),str(case/'config.ini')],capture_output=True,text=True,timeout=60)
    (case/'run.log').write_text(p.stdout+p.stderr)
    if fail:
        assert p.returncode!=0 and ('weight' in (p.stdout+p.stderr).lower() or 'column' in (p.stdout+p.stderr).lower()), (name,p.stdout,p.stderr)
    else:
        assert p.returncode==0,(name,p.stdout,p.stderr)
    return case/'output'


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--template',type=Path,required=True)
    parser.add_argument('--window',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True)
    a=parser.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    passed=[]
    for dims,point in itertools.product((2,3),(False,True)):
        points=list(itertools.product(range(3),range(2),range(2 if dims==3 else 1)))
        # Deliberately not x-fast input order; points have two coincident records.
        if point: points.append(points[0])
        values=[.5 + ((i*7)%11)/4 for i in range(len(points))]
        weights=[1+(i%5)*.5 for i in range(len(points))]
        tag=f'{dims}d_{"points" if point else "grid"}'
        path=a.out/(tag+'.gslib');source(path,points,values,weights)
        for module in ('template','window'):
            exe=getattr(a,module)
            for order in (range(2,5) if module=='template' else range(1,5)):
                for bell,moment in ((False,False),(True,False),(False,True)):
                    name=f'{tag}_{module}_{order}_{bell}_{moment}'
                    out=run(exe,a.out,name,path,dims,point,module,order,bell=bell,moment=moment)
                    actual=output_values(out)
                    if module=='window':
                        expected=[joint([[v]*order for v in values],weights,order,moment)]*(len(points) if point else 1)
                    else:
                        expected=[]
                        for lag in (0,1):
                            events=[];masses=[]
                            for i,p in enumerate(points):
                                ids=[i] if lag==0 else [j for j,q in enumerate(points) if q==(p[0]+1,p[1],p[2])]
                                for j in ids:
                                    events.append([values[i]]+[values[j]]*(order-1))
                                    masses.append(weights[i] if i==j else math.sqrt(weights[i]*weights[j]))
                            expected.append(joint(events,masses,order,moment))
                    close(actual,expected,name);passed.append(name)
            # All-ones equivalence, positive scaling invariance, missing source filtering.
            baseline=output_values(run(exe,a.out,tag+module+'_base',path,dims,point,module))
            for label,ws,mode in [('scale',[w*17 for w in weights],'external'),('ones',[1]*len(weights),'external'),('uniform',[1]*len(weights),'uniform')]:
                p=a.out/(tag+label+'.gslib');source(p,points,values,ws)
                result=output_values(run(exe,a.out,tag+module+label,p,dims,point,module,mode=mode))
                if label=='ones': ones=result
                else: close(result,baseline if label=='scale' else ones,tag+module+label)
                passed.append(tag+module+label)
            masked=values[:];masked[1]=math.nan
            wm=weights[:];wm[1]=math.nan
            p=a.out/(tag+'_missing.gslib');source(p,points,masked,wm)
            result=output_values(run(exe,a.out,tag+module+'_missing',p,dims,point,module))
            if module=='window':
                good=[i for i,v in enumerate(masked) if math.isfinite(v)]
                expected=joint([[values[i]]*4 for i in good],[weights[i] for i in good],4)
                close(result,[expected]*len(result),'missing')
            passed.append(tag+module+'_missing')
            # Remove one record entirely: Cartesian densification must carry
            # remaining weights with their source rows, including audit IDs.
            keep=[i for i in range(len(points)) if i!=1]
            sparse=a.out/(tag+'_sparse.gslib')
            source(sparse,[points[i] for i in keep],[values[i] for i in keep],[weights[i] for i in keep])
            sparse_result=output_values(run(exe,a.out,tag+module+'_sparse',sparse,dims,point,module))
            if module=='template' or not point:
                close(sparse_result,result,'sparse/missing equivalence')
            else:
                close(sparse_result,[expected]*len(keep),'sparse points')
            passed.append(tag+module+'_sparse')
            for bad in ('nan','inf','0','-1','garbage'):
                ws=weights[:];ws[0]=bad;p=a.out/(tag+'_invalid.gslib');source(p,points,values,ws)
                run(exe,a.out,tag+module+'_reject_'+bad,p,dims,point,module,fail=True)
                passed.append(tag+module+'_reject_'+bad)
            for label in ('','absent'):
                run(exe,a.out,tag+module+'_property_'+label,path,dims,point,module,prop=label,fail=True)
                passed.append(tag+module+'_property_'+label)
            duplicate=a.out/(tag+'_duplicate_column.gslib')
            duplicate.write_text('duplicate weight columns\n6\nx\ny\nz\nvalue\nweight\nweight\n'+''.join(
                ' '.join(map(str,(*p,v,w,w)))+'\n' for p,v,w in zip(points,values,weights)))
            run(exe,a.out,tag+module+'_duplicate_column',duplicate,dims,point,module,fail=True)
            passed.append(tag+module+'_duplicate_column')
            if not point:
                duplicate_node=a.out/(tag+'_duplicate_node.gslib')
                source(duplicate_node,points+[points[0]],values+[values[0]],weights+[weights[0]])
                run(exe,a.out,tag+module+'_duplicate_node',duplicate_node,dims,point,module,fail=True)
                passed.append(tag+module+'_duplicate_node')
            # Reimport existing declustering point masses without recomputation.
            suffix=',-1' if dims==3 else ''
            upper=',2' if dims==3 else ''
            for method in ('cell_declustering','polygonal_declustering','kriging_declustering'):
                extra=f'''[declustering]
dimensions={dims}
cell_size={','.join(['2']*dims)}
shifts_per_axis=2
domain_min=-1,-1{suffix}
domain_max=3,2{upper}
kriging_model=exponential
kriging_ranges={','.join(['4']*dims)}
kriging_nugget_fraction=0.5
kriging_integration_cells={','.join(['3']*dims)}
'''
                out=run(exe,a.out,tag+module+method,path,dims,point,module,mode=method,extra=extra)
                with (out/'result.point_weights.csv').open() as f: audit=list(csv.DictReader(f))
                byid={int(r.get('source_id',r['node_id'])):float(r['weight']) for r in audit}
                external=a.out/(tag+module+method+'.gslib')
                source(external,points,values,[byid[i+1] for i in range(len(points))])
                result=output_values(run(exe,a.out,tag+module+method+'_reimport',external,dims,point,module))
                close(result,output_values(out),tag+module+method);passed.append(tag+module+method+'_reimport')
    (a.out/'summary.json').write_text(json.dumps({'passed':len(passed),'checks':passed},indent=2))
    print(f'{len(passed)} external-weight checks passed')


if __name__=='__main__': main()
