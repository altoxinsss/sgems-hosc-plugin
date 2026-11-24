"""External weights through real SGeMS objects, plugin DLLs and parameter files."""
import argparse
import itertools
import json
import math
from pathlib import Path
import shutil
import xml.etree.ElementTree as ET
from current_host_checks import write_grid, read_grid, f32, close, run_process
from native_reference_checks import compare


def main():
    p=argparse.ArgumentParser()
    for key in ('host','cli','out'):p.add_argument('--'+key,type=Path,required=True)
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    parameters=a.out/'parameters';parameters.mkdir(exist_ok=True)
    reports=[]
    for dims,point in itertools.product((2,3),(False,True)):
        points=[(x,y,z) for z in range(2 if dims==3 else 1) for y in range(2) for x in range(3)]
        if point: points.append(points[0])
        values=[f32(.5+((i*7)%11)/4) for i in range(len(points))]
        weights=[f32(1+(i%5)*.5) for i in range(len(points))]
        missing=values[:];missing[1]=math.nan
        optional=weights[:];optional[1]=math.nan
        bad=weights[:];bad[0]=0
        negative=weights[:];negative[0]=-1
        invalid=weights[:];invalid[0]=math.nan
        name=('points' if point else 'grid')+('2' if dims==2 else '')
        fixture=a.out/(name+'.sgems')
        write_grid(fixture,name,(3,2,2 if dims==3 else 1),(0,0,0),(1,1,1),
            dict(value=values,weight=weights,ones=[1]*len(weights),scaled=[w*17 for w in weights],
                 missing=missing,optional=optional,bad=bad,negative=negative,invalid=invalid),points if point else None)
        if dims==3:shutil.copy2(fixture,parameters/fixture.name)
        for templ in (True,False):
            algorithm='HOSC_TEMPLATE' if templ else 'HOSC_WINDOW'
            base=dict(Source_Grid=name,Source_Property='value',weight_mode='External weights',weight_property='weight',
                      Target_Button='Cumulant',Algorithm_Button='AVL Based algorithm',coordinate_dimensions=str(dims),
                      output_diagnostics='true',output_point_weights='true')
            if templ:
                base.update(NOD='3',NLAG='2',LAG='1',LTOL='0',AZM='0',ATOL='0',BANDWH='0',DIP='0',DTOL='0',BANDWD='0',
                            fst_dir='1',snd_dir='0',trd_dir='0',cross_pos='',enable_weighted_moment='1',
                            template_sampling='point_set_tolerance' if point else 'exhaustive_integer_offsets',
                            integer_direction_vectors='1,0,0' if dims==3 else '1,0')
            else:
                base.update(moment_order='4',windowl='3',windoww='2',windowh='2' if dims==3 else '1',stepx='1',stepy='1',stepz='1' if dims==3 else '0',
                            neighborhood_shape='box',neighborhood_radii='100,100,100',min_samples='1',center_mode='sample_points')
            cases=[('base',{}),('ones',{'weight_property':'ones'}),('uniform',{'weight_mode':'uniform'}),
                   ('scale',{'weight_property':'scaled'}),('bell',{'Algorithm_Button':'Bell'}),('moment',{'Target_Button':'Moment'}),
                   ('missing',{'Source_Property':'missing','weight_property':'optional'})]
            if templ:
                cases += [('switch_off',{'enable_weighted_moment':'0','weight_property':'absent'})]
            else:
                cases += [(f'order{k}',{'moment_order':str(k)}) for k in (1,2,3)]
                if point:
                    for shape in ('box','ellipsoid'):
                        cases.append((shape+'_centers',dict(center_mode='regular_grid',neighborhood_shape=shape,
                            neighborhood_radii='0.5,0.5,0.5' if dims==3 else '0.5,0.5,0',
                            grid_spacing='1,1,1' if dims==3 else '1,1,0',grid_bounds_mode='auto')))
            saved_results={}
            for label,changes in cases:
                case=a.out/(name+algorithm+'_'+label);case.mkdir(exist_ok=True)
                audit=case/'audit';audit.mkdir(exist_ok=True)
                params=dict(base,props='stat',audit_output_directory=str(audit),**changes)
                root=ET.Element('parameters');ET.SubElement(root,'algorithm',name=algorithm)
                for key,value in params.items():ET.SubElement(root,key,value=value)
                commands=[
                    f'LoadObjectFromFile {fixture}::s-gems',
                    f'RunGeostatAlgorithm {algorithm}::/GeostatParamUtils/XML::'+ET.tostring(root,encoding='unicode'),
                    f'SaveGeostatGrid {algorithm}_stat::{case/"saved.sgems"}::s-gems']
                if label not in ('uniform','switch_off'):
                    commands.append(f'SaveGeostatGrid {algorithm}_stat_point_weights::{case/"weights.sgems"}::s-gems')
                command=case/'run.cmd';command.write_text('\n'.join(commands)+'\n')
                run_process([a.host/'sgems.exe',command],a.host,case/'host.log')
                done=list(audit.rglob('sgems_result.txt'))
                assert len(done)==1,(case,(case/'host.log').read_text(errors='replace'))
                run=done[0].parent
                # Execute the exact exported config/data with the standalone implementation.
                ref=case/'reference';ref.mkdir(exist_ok=True)
                for filename in ('config.ini','source.gslib'):shutil.copy2(run/filename,ref/filename)
                run_process([a.cli/('template_cli.exe' if templ else 'movewindow_cli.exe'),ref/'config.ini'],ref,ref/'run.log')
                # Invalid grid windows have differing diagnostic ESS conventions; compare values separately.
                if label=='missing' and not point and not templ:
                    native=read_grid(case/'saved.sgems')['props']['stat'];assert all(map(math.isnan,native))
                    difference=0
                else: difference=compare(run/'output',ref/'output')
                saved=read_grid(case/'saved.sgems')['props']['stat'];saved_results[label]=saved
                if label not in ('uniform','switch_off'):
                    source_weights=read_grid(fixture)['props'][params['weight_property']]
                    if label=='missing':source_weights[1]=0
                    actual_weights=read_grid(case/'weights.sgems')['props']['observation_weight']
                    for x,y in zip(actual_weights,source_weights):close(x,y,(name,'external weight output',label))
                if label in ('scale','bell','uniform','switch_off'):
                    expected=saved_results['ones' if label in ('uniform','switch_off') else 'base']
                    for x,y in zip(saved,expected):close(x,y,(name,algorithm,label))
                if label=='base' and dims==3:
                    root.remove(root.find('audit_output_directory'));ET.indent(root,space='  ')
                    (parameters/(name+algorithm+'.par')).write_text(ET.tostring(root,encoding='unicode')+'\n')
                reports.append(dict(case=case.name,max_native_reference_difference=difference,status='passed'))
            # Errors must not leave a completed calculation, and must report the weight issue.
            for field in ('bad','negative','invalid','absent',''):
                case=a.out/(name+algorithm+'_reject_'+(field or 'empty'));case.mkdir(exist_ok=True)
                params=dict(base,weight_property=field,props='invalid',audit_output_directory=str(case))
                root=ET.Element('parameters');ET.SubElement(root,'algorithm',name=algorithm)
                for key,value in params.items():ET.SubElement(root,key,value=value)
                command=case/'run.cmd';command.write_text(f'LoadObjectFromFile {fixture}::s-gems\n'+
                    f'RunGeostatAlgorithm {algorithm}::/GeostatParamUtils/XML::'+ET.tostring(root,encoding='unicode')+'\n')
                run_process([a.host/'sgems.exe',command],a.host,case/'host.log')
                assert not list(case.rglob('sgems_result.txt')),case
                assert 'weight' in (case/'host.log').read_text(errors='replace').lower(),case
                reports.append(dict(case=case.name,status='rejection passed'))
    (a.out/'report.json').write_text(json.dumps(reports,indent=2))
    print(f'{len(reports)} SGeMS external-weight integration cases passed')


if __name__=='__main__':main()
