#!/usr/bin/env python3
"""Build the one-file R4 verifier from the transparent certificate package."""
import base64
import hashlib
import io
import json
import re
from pathlib import Path
import tarfile
import zlib

import verify_berenstein_r4_package as package

root=Path(__file__).resolve().parent
manifest={'version':1,'files':{name:hashlib.sha256((root/name).read_bytes()).hexdigest()
                               for name in sorted(package.expected_files())}}
(root/'CERTIFICATE_MANIFEST.json').write_text(json.dumps(manifest,sort_keys=True,indent=2))
manifest_hash=hashlib.sha256((root/'CERTIFICATE_MANIFEST.json').read_bytes()).hexdigest()
driver=root/'verify_berenstein_r4_package.py'
driver_source=driver.read_text()
driver_source,count=re.subn(r"^EXPECTED_MANIFEST_SHA256='[^']+'$",
                           "EXPECTED_MANIFEST_SHA256='"+manifest_hash+"'",
                           driver_source,count=1,flags=re.MULTILINE)
if count!=1:raise ValueError('manifest digest assignment missing')
driver.write_text(driver_source)
names=sorted(set(package.expected_files()+['CERTIFICATE_MANIFEST.json',
                                         'verify_berenstein_r4_package.py']))
buf=io.BytesIO()
with tarfile.open(fileobj=buf,mode='w') as tf:
    for name in names:
        data=(root/name).read_bytes()
        item=tarfile.TarInfo(name)
        item.size=len(data);item.mode=0o644;item.mtime=0
        tf.addfile(item,io.BytesIO(data))
payload=zlib.compress(buf.getvalue(),9)
digest=hashlib.sha256(payload).hexdigest()
encoded=base64.b85encode(payload).decode('ascii')
lines='\n'.join('    '+repr(encoded[i:i+120]) for i in range(0,len(encoded),120))
header='''#!/usr/bin/env python3
"""Self-contained R4 Berenstein Arb verifier.

This single file embeds the transparent certificate package in the same
directory. `--recompute` regenerates every bounded batch before PROVED.
Build source: build_standalone_r4.py; extracted driver source:
verify_berenstein_r4_package.py. Requires python-flint 0.9.0.
"""
import base64
import hashlib
import io
import json
import re
import os
from pathlib import Path
import subprocess
import sys
import tarfile
import tempfile
import zlib

PAYLOAD_SHA256='__DIGEST__'
PAYLOAD_B85=(
__PAYLOAD__
)


def main():
    if sys.flags.optimize:
        raise RuntimeError('Python optimisation disables certificate assertions')
    compressed=base64.b85decode(PAYLOAD_B85.encode('ascii'))
    if hashlib.sha256(compressed).hexdigest()!=PAYLOAD_SHA256:
        raise RuntimeError('embedded certificate payload changed')
    archive=zlib.decompress(compressed)
    with tempfile.TemporaryDirectory(prefix='r4_certificate_') as tmp:
        root=Path(tmp)
        with tarfile.open(fileobj=io.BytesIO(archive),mode='r:') as tf:
            for item in tf.getmembers():
                if not item.isfile() or '/' in item.name or item.name.startswith('.'):
                    raise RuntimeError('invalid embedded member')
                f=tf.extractfile(item)
                if f is None:raise RuntimeError('missing embedded member')
                (root/item.name).write_bytes(f.read())
        inner=root/'verify_berenstein_r4_package.py'
        if not inner.is_file():raise RuntimeError('missing embedded verifier')
        env=dict(os.environ)
        env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',
                   PYTHONDONTWRITEBYTECODE='1')
        result=subprocess.run([sys.executable,'-B',str(inner),*sys.argv[1:]],
                              cwd=root,env=env,capture_output=True,text=True)
        sys.stdout.write(result.stdout)
        sys.stderr.write(result.stderr)
        if result.returncode:raise SystemExit(result.returncode)


if __name__=='__main__':main()
'''
source=header.replace('__DIGEST__',digest).replace('__PAYLOAD__',lines)
(root/'verify_berenstein_r4.py').write_text(source)
print('embedded files',len(names),'compressed bytes',len(payload),
      'standalone bytes',len(source),'payload sha256',digest)
