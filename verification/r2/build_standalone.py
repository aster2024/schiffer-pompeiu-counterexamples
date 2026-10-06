#!/usr/bin/env python3
"""Embed byte-identical local stage sources and the two frozen decimal inputs."""
import base64
import hashlib
import json
from pathlib import Path
import zlib

root=Path(__file__).resolve().parent
names=('spectral_interval.py','analytic_cap.py','boundary_interval.py','boundary_newton.py',
       'identity_checks.py','coverage_checks.py','scalar_second_iv.py','final_gates.py',
       'negative_controls.py')
bundle={name:(root/name).read_text() for name in names}
bundle['candidate.json']=(root/'center_j40_frozen.json').read_text()
bundle['reference.json']=(root/'reference.json').read_text()
raw=json.dumps(bundle,sort_keys=True,separators=(',',':')).encode()
encoded=base64.b64encode(zlib.compress(raw,9)).decode()
constant='EMBEDDED_SHA256 = '+repr(hashlib.sha256(raw).hexdigest())+'\n'
constant+='EMBEDDED_BUNDLE = (\n'+''.join('    '+repr(encoded[i:i+100])+'\n' for i in range(0,len(encoded),100))+')'
driver=(root/'standalone_driver.py').read_text()
if driver.count('# REPLACE_EMBEDDED_BUNDLE')!=1:
    raise RuntimeError('standalone marker mismatch')
target=root/'verify_planar_convex_v2.py'
target.write_text(driver.replace('# REPLACE_EMBEDDED_BUNDLE',constant))
manifest={name:hashlib.sha256(content.encode()).hexdigest() for name,content in bundle.items()}
manifest['verify_planar_convex_v2.py']=hashlib.sha256(target.read_bytes()).hexdigest()
(root/'BUILD_MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n')
print(target.name,manifest[target.name],len(target.read_bytes()))
