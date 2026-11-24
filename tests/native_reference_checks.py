"""Optional native-plugin versus developer-reference comparisons.

Runs in the isolated host; the x64 programs are test references only.
"""
from pathlib import Path
import argparse
import csv
import json
import math
import shutil
import xml.etree.ElementTree as ET
from current_host_checks import run_process, read_grid, close, write_grid, f32


def records(path):
    return list(csv.DictReader(path.open(newline='')))


def compare(left, right, missing_fixture=False):
    errors = []
    def equal(a, b, context, tol=2e-10):
        a, b = float(a), float(b)
        if math.isnan(a) and math.isnan(b):return
        assert math.isfinite(a) and math.isfinite(b), (context,a,b)
        error=abs(a-b)
        assert error <= 2e-12+tol*abs(b), (context,a,b,error)
        errors.append(error)
    def values(p):return [list(map(float,line.split())) for line in p.read_text().splitlines()[6:] if line.strip()]
    a,b=values(left/'result.txt'),values(right/'result.txt')
    assert len(a)==len(b),(left,len(a),len(b))
    for i,(x,y) in enumerate(zip(a,b)):
        for k in range(4):equal(x[k],y[k],('result',i,k))
    for name,fields in [('result.support.csv',('anchor_count','support_count','value','weight_sum','ESS')),
                        ('result.window_support.csv',('sample_count','effective_sample_size','valid')),
                        ('result.point_weights.csv',('x','y','z','value','weight'))]:
        if not (right/name).exists():continue
        a,b=records(left/name),records(right/name)
        assert len(a)==len(b),(name,len(a),len(b))
        for i,(x,y) in enumerate(zip(a,b)):
            if missing_fixture and name=='result.window_support.csv' and i==0:
                # The exported lattice has one NaN in this eight-cell window.
                # The old dense loader counted its record as present. The DLL
                # reports seven informed cells and zero ESS for invalid output.
                assert (x['sample_count'],x['effective_sample_size'],x['valid'])==('7','0','0')
                assert (y['sample_count'],y['effective_sample_size'],y['valid'])==('8','8','0')
                continue
            for field in fields:
                if field in x and field in y:equal(x[field],y[field],(name,i,field))
    return max(errors,default=0)


def main():
    p=argparse.ArgumentParser()
    for arg in ('host','suite','reference','cli','out'):p.add_argument('--'+arg,type=Path,required=True)
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    report=[]
    for entry in json.loads((a.suite/'report.json').read_text()):
        name=entry['case'];old=list((a.reference/name/'audit').rglob('sgems_result.txt'));assert len(old)==1
        maximum=compare(Path(entry['run'])/'output',old[0].parent/'output',missing_fixture=name=='window_missing')
        report.append(dict(case=name,max_absolute_difference=maximum,status='passed'))
    print('Original suite: native DLL / independent reference matched',flush=True)
    cases=[]
    def add(label,base,**changes):cases.append((label,base,changes))
    for order in (1,5,9):
        add('window_order'+str(order),'window_grid',moment_order=str(order),output_standardized_cumulant='false')
    add('window_raw_fifth','window_grid',Target_Button='Moment',moment_order='5',output_standardized_cumulant='false')
    add('window_rotated3d','window_points',neighborhood_radii='1.7,.7,.55',neighborhood_angles='35,25,15',moment_order='5')
    add('window_rotated2d','window_points2',neighborhood_radii='1.5,.8,0',neighborhood_angles='25,0,0',moment_order='5')
    add('template_fifth','template_exact',NOD='2,1,1')
    add('template_legacy_fifth','template_legacy',NOD='2,1,1')
    add('template_hidden','template_exact',fst_dir='3',snd_dir='1',trd_dir='0',cross_pos='2')
    add('template_points_hidden','template_points',fst_dir='1',snd_dir='0',trd_dir='0',cross_pos='2')
    add('template_points_repeated','template_points',AZM='90,90',DIP='0,0',LAG='.6,.6',LTOL='.5,.5',NOD='2,2')
    add('template_negative_offset','template_exact',integer_direction_vectors='0,-1,0;-1,0,0;0,0,1',NOD='2,1,1')
    add('template_ninth','template_exact',NOD='4,2,2',NLAG='1,1,1')
    add('cell_fifth_bell','points_cell_declustering',Algorithm_Button='Bell')
    add('polygonal_fifth_bell','points_polygonal_declustering',Algorithm_Button='Bell')
    add('kriging_spherical','points_kriging_declustering',kriging_model='spherical')
    add('kriging_gaussian','points_kriging_declustering',kriging_model='gaussian')
    add('point_empty_support','template_points',minimum_template_count='100000')
    add('grid_xy_window','window_grid',Source_Grid='grid2',windowl='3',windoww='2',windowh='1',stepx='2',stepy='1',stepz='0',moment_order='5')
    add('grid_xy_template','template_legacy',Source_Grid='grid2',NOD='2,1,1',NLAG='2,2,2',DIP='0,0,0',AZM='0,90,180',trd_dir='0',cross_pos='1')
    values=[f32(.13+.02*x+.015*y+.017*((x*y+2)%4)) for y in range(7) for x in range(8)]
    write_grid(a.out/'grid2.sgems','grid2',(8,7,1),(0,0,0),(2,3,1),dict(phi=values))
    for label,base,changes in cases:
        target=a.out/label;target.mkdir(exist_ok=True);audit=target/'audit';audit.mkdir(exist_ok=True)
        original=(a.suite/base/'run.cmd').read_text().splitlines()[1].split('::',2)[2]
        root=ET.fromstring(original)
        for key,val in dict(changes,props=label,audit_output_directory=str(audit)).items():
            elem=root.find(key)
            if elem is None:elem=ET.SubElement(root,key)
            elem.set('value',val)
        algo=root.find('algorithm').get('name');grid=root.find('Source_Grid').get('value')
        fixture=(a.out if grid=='grid2' else a.suite)/(grid+'.sgems')
        cmd=target/'run.cmd';cmd.write_text(f'LoadObjectFromFile {fixture}::s-gems\nRunGeostatAlgorithm {algo}::/GeostatParamUtils/XML::{ET.tostring(root,encoding="unicode")}\n')
        run_process([a.host/'sgems.exe',cmd],a.host,target/'host.log')
        completed=list(audit.rglob('sgems_result.txt'));assert len(completed)==1,(label,(target/'host.log').read_text(errors='replace'))
        run=completed[0].parent;ref=target/'reference';ref.mkdir(exist_ok=True)
        for name in ('source.gslib','config.ini'):shutil.copy2(run/name,ref/name)
        exe=a.cli/('template_cli.exe' if algo=='HOSC_TEMPLATE' else 'movewindow_cli.exe')
        run_process([exe,ref/'config.ini'],ref,ref/'reference.log')
        maximum=compare(run/'output',ref/'output')
        report.append(dict(case=label,max_absolute_difference=maximum,status='passed'))
        print(label+': matched',flush=True)
    # Legacy sampler globals and evaluator globals must survive interleaving.
    commands=[f'LoadObjectFromFile {a.suite/"grid.sgems"}::s-gems']
    repeat=a.out/'interleaved';repeat.mkdir(exist_ok=True)
    expected=[]
    for i in range(12):
        base='template_legacy' if i%2==0 else 'window_grid'
        root=ET.fromstring((a.suite/base/'run.cmd').read_text().splitlines()[1].split('::',2)[2])
        root.find('props').set('value','mixed'+str(i));root.find('audit_output_directory').set('value',str(repeat))
        algo=root.find('algorithm').get('name')
        commands.append(f'RunGeostatAlgorithm {algo}::/GeostatParamUtils/XML::{ET.tostring(root,encoding="unicode")}')
        commands.append(f'SaveGeostatGrid {algo}_mixed{i}::{repeat/(str(i)+".sgems")}::s-gems')
        expected.append(read_grid(a.suite/base/'saved.sgems')['props'][base])
    cmd=repeat/'run.cmd';cmd.write_text('\n'.join(commands)+'\n');run_process([a.host/'sgems.exe',cmd],a.host,repeat/'host.log')
    for i,ref in enumerate(expected):
        actual=read_grid(repeat/(str(i)+'.sgems'))['props']['mixed'+str(i)]
        for x,y in zip(actual,ref):close(x,y,('interleaved',i))
    print('12 interleaved template/window calls: passed')
    (a.out/'report.json').write_text(json.dumps(dict(comparisons=report,interleaved_calls=12,
        intentional_diagnostic_difference='window_missing: invalid eight-cell window reports seven informed samples and ESS=0; old dense loader counted the NaN record (sample_count=ESS=8). Result and validity are unchanged.'),indent=2)+'\n')


if __name__=='__main__':main()
