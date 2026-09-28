"""Aggregate complete R14 finite-column receipts and compute Arb Y, ||A||."""

import hashlib
import json
from pathlib import Path
from flint import arb, ctx
from finite_rank2 import load_inverse
import rank2_cap_core as v


RANGES=[(0,20)]+[(a,a+40) for a in range(20,580,40)]+[(580,600)]
RANGES += [(a,a+1) for a in range(600,610)]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    ctx.prec=128
    ctx.threads=1
    d=json.loads(Path('center_r14_M114_S60.json').read_text())
    rows,g,c,p,V,F=v.make_algebra(d)
    N=len(rows)
    v.need(N==610,'unexpected finite dimension')
    rho=arb(11)/10
    A=load_inverse('inverse_r14_M114_S60.npy')
    inv=v.Inverse(d['M'],d['S'],rho,c[1],A)
    cur=0
    defs=[]
    zs=[]
    for a,b in RANGES:
        v.need(a==cur,'finite batch list gap')
        item=json.loads(Path(f'finite_r14_{a:03d}_{b:03d}.json').read_text())
        v.need(item['kind']=='finite_batch' and item['dimension']==14
               and item['start']==a and item['stop']==b and item['N']==N,
               'finite batch receipt mismatch')
        defs.append(arb(item['inverse_defect_max']))
        zs.append(arb(item['Z_finite_max']))
        cur=b
    v.need(cur==N,'finite columns incomplete')
    defect=v.max_upper(defs)
    Zfinite=v.max_upper(zs)
    Y=inv.apply_many([F])[0].upper()
    Anfin=v.max_upper(sum((abs(A[i,j])*inv.cw[i] for i in range(N)),v.Z)/rho**rows[j][0]
                      for j in range(N))
    out={'status':'OPEN_COMPONENT','dimension':14,'N':N,'bits':ctx.prec,
         'M':d['M'],'S':d['S'],'batch_count':len(RANGES),
         'centre_sha256':sha('center_r14_M114_S60.json'),
         'inverse_sha256':sha('inverse_r14_M114_S60.npy'),
         'inverse_defect':defect.str(30),'Z_finite':Zfinite.str(30),
         'Y':Y.str(30),'A_finite':Anfin.str(30),
         'raw_residual_norm':v.norm(F,rho).str(30)}
    Path('finite_bound_r14_M114_S60.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))


if __name__=='__main__':
    main()
