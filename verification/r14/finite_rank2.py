"""Frozen finite inverse and Arb finite-column receipts for R10/R14.

The inverse is computed numerically, then checked by interval arithmetic.
"""

import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from flint import arb, arb_mat, ctx
import rank2_cap_core as v


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_inverse(path):
    raw=np.load(path)
    return arb_mat([[v.dyadic(float(raw[i,j]).hex())
                     for j in range(raw.shape[1])] for i in range(raw.shape[0])])


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--centre',required=True)
    p.add_argument('--inverse',required=True)
    p.add_argument('--output',required=True)
    p.add_argument('--create-inverse',action='store_true')
    args=p.parse_args()
    ctx.prec=128
    ctx.threads=1
    d=json.loads(Path(args.centre).read_text())
    rows,g,c,pcoef,V,F=v.make_algebra(d)
    N=len(rows)
    rho=arb(11)/10
    cols=v.finite_columns(rows,g,c,pcoef,V,d['S'])
    mat=arb_mat([[col.get(key,v.Z) for col in cols] for key in rows])
    if args.create_inverse:
        mid=np.array([[float(mat[i,j].mid()) for j in range(N)] for i in range(N)])
        frozen=np.linalg.inv(mid)
        v.need(np.all(np.isfinite(frozen)),'nonfinite approximate inverse')
        np.save(args.inverse,frozen)
    A=load_inverse(args.inverse)
    inv=v.Inverse(d['M'],d['S'],rho,c[1],A)
    E=arb_mat(N,N)
    for i in range(N):
        E[i,i]=1
    E=E-A*mat
    defect=v.max_upper(sum((abs(E[i,j])*inv.cw[i] for i in range(N)),v.Z)/inv.cw[j]
                       for j in range(N))
    Zfinite=v.max_upper(val/inv.cw[i] for i,val in enumerate(
        inv.apply_many(cols,subtract=[('finite',i) for i in range(N)])))
    Y=inv.apply_many([F])[0].upper()
    Anfin=v.max_upper(sum((abs(A[i,j])*inv.cw[i] for i in range(N)),v.Z)/rho**rows[j][0]
                      for j in range(N))
    out={'status':'OPEN_COMPONENT','dimension':2*d['m']+2,'N':N,'bits':ctx.prec,
         'M':d['M'],'S':d['S'],'centre_sha256':sha(args.centre),
         'inverse_sha256':sha(args.inverse),'inverse_defect':defect.str(30),
         'Z_finite':Zfinite.str(30),'Y':Y.str(30),'A_finite':Anfin.str(30),
         'raw_residual_norm':v.norm(F,rho).str(30)}
    Path(args.output).write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))


if __name__=='__main__':
    main()
