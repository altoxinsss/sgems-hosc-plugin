"""Make local forwarding headers from the official Qt 4.7.3 source archive."""
from pathlib import Path
import re
import json
import argparse

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--sdk-root', type=Path, required=True)
args = parser.parse_args()
root = args.sdk_root.resolve()/'qt-4.7.3'
if not (root/'src/corelib/global/qglobal.h').is_file():
    parser.error('Qt 4.7.3 source headers not found under --sdk-root')
modules = {'QtCore':'corelib','QtGui':'gui','QtXml':'xml'}
manifest=[]
for module,sub in modules.items():
    output=root/'include'/module
    output.mkdir(parents=True,exist_ok=True)
    for header in (root/'src'/sub).rglob('*.h'):
        if header.name.endswith('_p.h') or header.name.endswith('_p_p.h'):
            continue
        target=output/header.name
        target.write_text('#include "'+header.as_posix()+'"\n',encoding='utf-8')
        source=header.read_text(encoding='utf-8',errors='replace')
        declarations=re.sub(r'//[^\n]*','',source)
        classes=set(re.findall(r'\bclass\s+(?:Q_\w+_EXPORT\s+)?(Q[A-Za-z0-9_]+)\s*[:{]',declarations))
        classes.update(re.findall(r'\bstruct\s+(?:Q_\w+_EXPORT\s+)?(Q[A-Za-z0-9_]+)\s*[:{]',source))
        classes.update(re.findall(r'#pragma\s+qt_class\((\w+)\)',source))
        # Qt public container classes can put a newline after the class name.
        for name in classes:
            (output/name).write_text('#include "'+header.name+'"\n',encoding='utf-8')
        manifest.append({'header':str(header.relative_to(root)),'module':module,'classes':sorted(classes)})
# Well-known typedef exports in Qt headers (no class declaration).
for name,header in {'QWidgetList':'qwindowdefs.h','QStringList':'qstringlist.h','QObjectList':'qobject.h','QtGlobal':'qglobal.h'}.items():
    module='QtGui' if name=='QWidgetList' else 'QtCore'
    (root/'include'/module/name).write_text('#include "'+header+'"\n',encoding='utf-8')
(root/'forwarding_headers_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
uic=root/'src/tools/uic'
if uic.exists():
    (uic/'qclass_lib_map.h').write_text(''.join('QT_CLASS_LIB('+name+', '+m['module']+', '+Path(m['header']).name+')\n' for m in manifest for name in m['classes']),encoding='utf-8')
designer=root/'include/QtDesigner'
designer.mkdir(exist_ok=True)
for sub in ('uilib','sdk'):
    for header in (root/'tools/designer/src/lib'/sub).glob('*.h'):
        (designer/header.name).write_text('#include "'+header.as_posix()+'"\n',encoding='utf-8')
for cls,h in {'QFormBuilder':'formbuilder.h','QAbstractFormBuilder':'abstractformbuilder.h','QDesignerExportWidget':'qdesignerexportwidget.h'}.items():
    (designer/cls).write_text('#include "'+h+'"\n',encoding='utf-8')
uitools=root/'include/QtUiTools';uitools.mkdir(exist_ok=True)
(uitools/'QUiLoader').write_text('#include "'+(root/'tools/designer/src/uitools/quiloader.h').as_posix()+'"\n',encoding='utf-8')
print('Forwarding headers:',len(manifest))
