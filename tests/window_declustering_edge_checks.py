"""Weighted-grid missing values, display transforms, and retained order behavior."""
import argparse
import csv
import json
import math
from pathlib import Path
import shutil
import xml.etree.ElementTree as ET
from current_host_checks import run_process, read_grid, close
from native_reference_checks import compare, records


def log_series_cumulant(values, weights, order):
    total=math.fsum(weights);mean=math.fsum(w*v for w,v in zip(weights,values))/total
    a=[0]+[math.fsum(w*(v-mean)**n for w,v in zip(weights,values))/total/math.factorial(n) for n in range(1,order+1)]
    power=[1]+[0]*order;coefficient=0
    for k in range(1,order+1):
        power=[sum(power[j]*a[n-j] for j in range(n+1)) for n in range(order+1)]
        coefficient+=(-1)**(k+1)*power[order]/k
    return coefficient*math.factorial(order)


def main():
    p=argparse.ArgumentParser()
    for k in ('host','suite','weighted','cli','out'):p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True);report=[]
    for method in ('cell_declustering','polygonal_declustering','kriging_declustering'):
        label='grid_'+method+'_4_False_False'
        original=ET.fromstring((a.weighted/label/'run.cmd').read_text().splitlines()[1].split('::',2)[2])
        for case_name,changes in [
            ('missing',dict(Source_Property='missing',output_standardized_cumulant='false')),
            ('display',dict(invert_window_output='true',enable_percentile_clipping='true',percentile_vmin='10',percentile_vmax='90')),
            ('avl9',dict(moment_order='9',output_standardized_cumulant='false')),
            ('bell10',dict(moment_order='10',Algorithm_Button='Bell',output_standardized_cumulant='false')),
            ('avl10_rejected',dict(moment_order='10',output_standardized_cumulant='false'))]:
            label=method+'_'+case_name;folder=a.out/label;folder.mkdir(exist_ok=True);audit=folder/'audit';audit.mkdir(exist_ok=True)
            root=ET.fromstring(ET.tostring(original))
            for k,v in dict(changes,props=label,audit_output_directory=str(audit)).items():
                el=root.find(k)
                if el is None:el=ET.SubElement(root,k)
                el.set('value',v)
            commands=[
                f'LoadObjectFromFile {a.suite/"grid.sgems"}::s-gems',
                'RunGeostatAlgorithm HOSC_WINDOW::/GeostatParamUtils/XML::'+ET.tostring(root,encoding='unicode')]
            if case_name!='avl10_rejected':
                commands.append(f'SaveGeostatGrid HOSC_WINDOW_{label}::{folder/"saved.sgems"}::s-gems')
            cmd=folder/'run.cmd';cmd.write_text('\n'.join(commands)+'\n')
            run_process([a.host/'sgems.exe',cmd],a.host,folder/'host.log')
            completed=list(audit.rglob('sgems_result.txt'))
            if case_name=='avl10_rejected':
                assert not completed
                assert 'encoding' in (folder/'host.log').read_text(errors='replace')
                run=next(audit.iterdir())
            else:
                assert len(completed)==1,(label,(folder/'host.log').read_text(errors='replace'))
                run=completed[0].parent
            ref=folder/'reference';ref.mkdir(exist_ok=True)
            for f in ('source.gslib','config.ini'):shutil.copy2(run/f,ref/f)
            if case_name=='avl10_rejected':
                import subprocess
                result=subprocess.run([str(a.cli/'movewindow_cli.exe'),str(ref/'config.ini')],capture_output=True,text=True,timeout=30)
                assert result.returncode!=0 and 'encoding' in result.stderr
            else:
                run_process([a.cli/'movewindow_cli.exe',ref/'config.ini'],ref,ref/'log.txt')
                output=run/'output';support=records(output/'result.window_support.csv')
                if case_name=='missing':
                    other=records(ref/'output/result.window_support.csv')
                    for x,y in zip(support,other):
                        for key in ('sample_count','weight_sum','weight_squared_sum','effective_sample_size','valid'):
                            close(float(x[key]),float(y[key]),(label,key),atol=2e-10)
                    assert support[0]['sample_count']=='7' and support[0]['valid']=='0'
                    assert float(records(output/'result.point_weights.csv')[0]['weight'])==0
                    # CLI omits invalid Cartesian rows; compare the retained
                    # rows against the same locations in the native output.
                    native=[r for r in (output/'result.txt').read_text().splitlines()[6:] if r.strip()]
                    other=(ref/'output/result.txt').read_text().splitlines()[6:]
                    for x,y in zip(native[1:],other):
                        close(float(x.split()[3]),float(y.split()[3]),label,atol=2e-10)
                else:
                    compare(output,ref/'output')
                if case_name=='display':
                    saved=read_grid(folder/'saved.sgems')['props']
                    original_label='grid_'+method+'_4_False_False'
                    orig=read_grid(a.weighted/original_label/'saved.sgems')['props']
                    for x,y in zip(saved[label+'_standardized'],orig[original_label+'_standardized']):
                        close(x,y,'standardization independent of display')
                if case_name in ('avl9','bell10'):
                    weights=records(output/'result.point_weights.csv')
                    ids=[0,1,5,6,20,21,25,26]
                    expected=log_series_cumulant([float(weights[j]['value']) for j in ids],[float(weights[j]['weight']) for j in ids],int(changes['moment_order']))
                    got=float((output/'result.txt').read_text().splitlines()[6].split()[3])
                    close(got,expected,'log generating-function coefficients',atol=2e-12)
            report.append(label)
            print(label+': passed',flush=True)
    (a.out/'report.json').write_text(json.dumps(dict(passed=len(report),cases=report),indent=2))


if __name__=='__main__':main()
