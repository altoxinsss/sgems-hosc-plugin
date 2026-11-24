"""Exercise weighted windows through the actual SGeMS DLL and saved objects."""
import argparse
import configparser
import importlib.util
import itertools
import json
import math
from pathlib import Path
import shutil
import xml.etree.ElementTree as ET
from current_host_checks import run_process, read_grid, close, write_grid, f32
from native_reference_checks import compare, records


def main():
    p=argparse.ArgumentParser()
    for key in ('host','suite','cli','out','reference_tests'):
        p.add_argument('--'+key.replace('_','-'),type=Path,required=True)
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    spec=importlib.util.spec_from_file_location('window_reference',a.reference_tests/'window_declustering_checks.py')
    reference=importlib.util.module_from_spec(spec);spec.loader.exec_module(reference)
    methods=['cell_declustering','polygonal_declustering','kriging_declustering']
    # Reuse already independently validated template settings and source objects.
    fixtures={g:read_grid(a.suite/(g+'.sgems')) for g in ('grid','points','points2')}
    g2values=[f32(.2+.03*x+.05*y+.02*((x*y)%3)) for y in range(4) for x in range(5)]
    write_grid(a.out/'grid2.sgems','grid2',(5,4,1),(10.25,-3.5,100),(2.5,1.25,.5),dict(phi=g2values))
    fixtures['grid2']=read_grid(a.out/'grid2.sgems')
    report=[];parameters=a.out/'test_parameters';parameters.mkdir(exist_ok=True)
    shutil.copy2(a.suite/'grid.sgems',parameters/'grid.sgems')
    shutil.copy2(a.suite/'points.sgems',parameters/'points.sgems')
    for grid,method in itertools.product(('grid','grid2','points','points2'),methods):
        base=('grid_' if grid.startswith('grid') else 'points2_' if grid=='points2' else 'points_')+method
        template=ET.fromstring((a.suite/base/'run.cmd').read_text().splitlines()[1].split('::',2)[2])
        common={key:template.find(key).get('value') for key in ('weight_mode','cell_size','shifts_per_axis','domain_min','domain_max','kriging_model','kriging_ranges','kriging_nugget_fraction','kriging_integration_cells')}
        if grid=='grid2':
            for key in ('cell_size','domain_min','domain_max','kriging_ranges','kriging_integration_cells'):
                common[key]=','.join(common[key].split(',')[:2])
        if grid.startswith('grid'):
            common.update(declustering_origin='10,-4'+(',99.9' if grid=='grid' else ''))
        common.update(Source_Grid=grid,Source_Property='phi',weight_mode=method,
            windowl='2',windoww='2',windowh='2' if grid=='grid' else '1',stepx='1',stepy='1',stepz='1' if grid=='grid' else '0',
            coordinate_dimensions='2' if grid.endswith('2') else '3',neighborhood_radii='1,1,1',
            neighborhood_angles='90,0,0',min_samples='3',includes_center='false' if grid=='points2' else 'true',
            output_diagnostics='true',output_point_weights='true')
        for order,moment,bell in [(1,False,False),(2,False,True),(4,False,False),(4,False,True),
                                  (5,False,False),(5,False,True),(4,True,False)]:
            label=f'{grid}_{method}_{order}_{moment}_{bell}'
            case=a.out/label;case.mkdir(exist_ok=True);audit=case/'audit';audit.mkdir(exist_ok=True)
            params=dict(common,props=label,moment_order=str(order),Target_Button='Moment' if moment else 'Cumulant',
                Algorithm_Button='Bell' if bell else 'AVL Based algorithm',output_standardized_cumulant='true' if order>=2 and not moment else 'false',audit_output_directory=str(audit))
            root=ET.Element('parameters');ET.SubElement(root,'algorithm',name='HOSC_WINDOW')
            for key,val in params.items():ET.SubElement(root,key,value=val)
            fixture=(a.out if grid=='grid2' else a.suite)/(grid+'.sgems')
            command=case/'run.cmd';command.write_text('\n'.join([
                f'LoadObjectFromFile {fixture}::s-gems',
                'RunGeostatAlgorithm HOSC_WINDOW::/GeostatParamUtils/XML::'+ET.tostring(root,encoding='unicode'),
                f'SaveGeostatGrid HOSC_WINDOW_{label}::{case/"saved.sgems"}::s-gems',
                f'SaveGeostatGrid HOSC_WINDOW_{label}_point_weights::{case/"weights.sgems"}::s-gems'])+'\n')
            run_process([a.host/'sgems.exe',command],a.host,case/'host.log')
            done=list(audit.rglob('sgems_result.txt'));assert len(done)==1,(label,(case/'host.log').read_text(errors='replace'))
            run=done[0].parent;out=run/'output';ref=case/'reference';ref.mkdir(exist_ok=True)
            for filename in ('source.gslib','config.ini'):shutil.copy2(run/filename,ref/filename)
            run_process([a.cli/'movewindow_cli.exe',ref/'config.ini'],ref,ref/'run.log')
            maximum=compare(out,ref/'output')
            ws=[float(r['weight']) for r in records(out/'result.point_weights.csv')]
            source=fixtures[grid];values=source['props']['phi'];expected=[]
            if grid.startswith('grid'):
                nx,ny,nz=source['dims'];wz=2 if grid=='grid' else 1
                for z,y,x in itertools.product(range(nz-wz+1),range(ny-1),range(nx-1)):
                    ids=[((z+dz)*ny+y+dy)*nx+x+dx for dz,dy,dx in itertools.product(range(wz),range(2),range(2))]
                    expected.append((ids,True))
            else:
                for i,pt in enumerate(source['points']):
                    ids=[j for j,q in enumerate(source['points']) if math.isfinite(values[j]) and
                        (params['includes_center']=='true' or i!=j) and sum((v-w)**2 for v,w in zip(pt,q))<=1+1e-12]
                    expected.append((ids,len(ids)>=3))
            assert reference.verify_output(out,expected,ws,values,order,moment),label
            # Reused source weights must be identical to the template module.
            if grid!='grid2':
                old=list((a.suite/base/'audit').rglob('result.point_weights.csv'));assert len(old)==1
                for x,y in zip(ws,[float(r['weight']) for r in records(old[0])]):reference.equal(x,y,'template/window weight reuse')
            saved=read_grid(case/'saved.sgems');w_saved=read_grid(case/'weights.sgems')
            for x,y in zip(w_saved['props']['declustering_weight'],ws):close(x,y,('stored weight',label))
            for i,(ids,valid) in enumerate(expected):
                vals=[values[j] for j in ids];weights=[ws[j] for j in ids]
                wanted=reference.statistic(vals,weights,order,moment) if valid else math.nan
                close(saved['props'][label][i],wanted,('saved statistic',label,i),atol=2e-8)
                if params['output_standardized_cumulant']=='true':
                    variance=reference.statistic(vals,weights,2) if valid else math.nan
                    standardized=wanted/variance**(order/2) if variance>1e-12 else math.nan
                    close(saved['props'][label+'_standardized'][i],standardized,('standardized',label,i),atol=2e-6)
            if grid=='points' and order==4 and not moment and not bell:
                # Portable import parameters: source object/phi and all weight
                # settings retained, transient audit directory omitted.
                root.find('props').set('value','window_'+method)
                root.remove(root.find('audit_output_directory'))
                ET.indent(root,space='  ')
                (parameters/(method+'.par')).write_text(ET.tostring(root,encoding='unicode')+'\n')
            report.append(dict(case=label,max_native_reference_difference=maximum,status='passed'))
            print(label+': passed',flush=True)
    # Interleave methods in one host. A failed configuration must not leak
    # source weights or evaluator state into the following valid run.
    repeat=a.out/'repeat';repeat.mkdir(exist_ok=True)
    commands=[f'LoadObjectFromFile {a.suite/"points.sgems"}::s-gems']
    for i in range(9):
        method=methods[i%3];root=ET.parse(parameters/(method+'.par')).getroot()
        root.find('props').set('value','repeat'+str(i));ET.SubElement(root,'audit_output_directory',value=str(repeat))
        if i==3:
            bad=ET.fromstring(ET.tostring(root));bad.find('weight_mode').set('value','cell_declustering');bad.find('cell_size').set('value','-1,1,1')
            commands.append('RunGeostatAlgorithm HOSC_WINDOW::/GeostatParamUtils/XML::'+ET.tostring(bad,encoding='unicode').replace('\n',''))
        commands.append('RunGeostatAlgorithm HOSC_WINDOW::/GeostatParamUtils/XML::'+ET.tostring(root,encoding='unicode').replace('\n',''))
        commands.append(f'SaveGeostatGrid HOSC_WINDOW_repeat{i}::{repeat/(str(i)+".sgems")}::s-gems')
    command=repeat/'run.cmd';command.write_text('\n'.join(commands)+'\n')
    run_process([a.host/'sgems.exe',command],a.host,repeat/'host.log')
    assert len(list(repeat.rglob('sgems_result.txt')))==9
    for i in range(9):
        actual=read_grid(repeat/(str(i)+'.sgems'))['props']['repeat'+str(i)]
        label=f'points_{methods[i%3]}_4_False_False'
        expected=read_grid(a.out/label/'saved.sgems')['props'][label]
        for x,y in zip(actual,expected):close(x,y,'method switching after failure')
    (a.out/'report.json').write_text(json.dumps(dict(cases=report,passed=len(report),repeat_calls=9,recovery_after_error=True),indent=2))
    print(f'{len(report)} host configurations, 9 repeated runs, and error recovery passed.')


if __name__=='__main__':main()
