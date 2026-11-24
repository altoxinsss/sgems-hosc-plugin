"""Generate x64 import libraries directly from the installed host exports.

These libraries are locally derived SDK artifacts, not official SDK binaries.
Compatible headers/compiler and real host tests remain mandatory.
"""
from pathlib import Path
import struct
import hashlib
import json
import subprocess
import os
import argparse

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--sdk-root', type=Path, required=True)
parser.add_argument('--host', type=Path, required=True)
args=parser.parse_args()
SDK=args.sdk_root.resolve()
HOST=args.host.resolve()
LIB=SDK/'winsdk71_extracted/vc_stdamd64/Program Files(64)/Microsoft Visual Studio 10.0/VC/bin/amd64/lib.exe'
if not LIB.is_file() or not (HOST/'sgems-x64.exe').is_file():
    parser.error('Missing VC2010 lib.exe or sgems-x64.exe')
ENV=os.environ.copy()
ENV['PATH']=str(SDK/'winsdk71_extracted/vc_stdamd64/Win/System64')+';'+ENV['PATH']
OUT=SDK/'host_imports'
OUT.mkdir(exist_ok=True)

def exports(path):
    b=path.read_bytes()
    u16=lambda p:struct.unpack_from('<H',b,p)[0]
    u32=lambda p:struct.unpack_from('<I',b,p)[0]
    pe=u32(0x3c)
    assert b[pe:pe+4]==b'PE\0\0' and u16(pe+4)==0x8664
    opt=pe+24
    assert u16(opt)==0x20b
    sections=[]
    for i in range(u16(pe+6)):
        p=opt+u16(pe+20)+40*i
        sections.append((u32(p+12),max(u32(p+8),u32(p+16)),u32(p+20),u32(p+36)))
    def off(rva):
        for r,n,p,f in sections:
            if r<=rva<r+n:return p+rva-r
        raise ValueError(('unmapped RVA',path,rva))
    def string(p):return b[p:b.index(b'\0',p)].decode('ascii')
    erva,esize=u32(opt+112),u32(opt+116)
    ex=off(erva)
    base=u32(ex+16);count=u32(ex+24)
    funcs,names,ords=off(u32(ex+28)),off(u32(ex+32)),off(u32(ex+36))
    result=[]
    for i in range(count):
        name=string(off(u32(names+4*i)));ordinal=u16(ords+2*i)
        address=u32(funcs+4*ordinal)
        data=not any(r<=address<r+n and f&0x20000000 for r,n,p,f in sections)
        if erva<=address<erva+esize:data=False
        result.append((name,base+ordinal,data))
    return result

manifest=[]
for name in ['GsTLAppli_actions','GsTLAppli_filters','GsTLAppli_utils','GsTLAppli_appli','GsTLAppli_grid','GsTLAppli_math','GsTLAppli_geostat','GsTLAppli_extragui','GsTLAppli_gui','GsTLwidgets','QtCore4','QtGui4','QtXml4','QtDesigner4']:
    dll=HOST/(name+'.dll')
    entries=exports(dll)
    definition=OUT/(name+'.def')
    definition.write_text('LIBRARY '+name+'.dll\nEXPORTS\n'+''.join(' '+n+' @'+str(o)+(' DATA' if d else '')+'\n' for n,o,d in entries),encoding='ascii')
    library=OUT/(name+'.lib')
    call=subprocess.run([str(LIB),'/nologo','/MACHINE:X64','/DEF:'+str(definition),'/OUT:'+str(library)],capture_output=True,text=True,encoding='utf-8',errors='replace',env=ENV)
    (OUT/(name+'.log')).write_text(call.stdout+call.stderr,encoding='utf-8')
    assert call.returncode==0,(name,call.stdout,call.stderr)
    manifest.append({'host_dll':str(dll),'sha256':hashlib.sha256(dll.read_bytes()).hexdigest(),'exports':len(entries),'generated_library':str(library)})
(OUT/'manifest.json').write_text(json.dumps(manifest,indent=2))
print('Generated import libraries:',len(manifest))
