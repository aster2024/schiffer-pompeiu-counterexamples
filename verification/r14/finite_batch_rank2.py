"""Arb finite Jacobian columns in bounded batches for large R14 matrices."""

import argparse
import json
from pathlib import Path
from flint import arb, arb_mat, ctx
from finite_rank2 import load_inverse
import rank2_cap_core as v


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--centre',required=True)
    p.add_argument('--inverse',required=True)
    p.add_argument('--start',type=int,required=True)
    p.add_argument('--stop',type=int,required=True)
    p.add_argument('--output',required=True)
    args=p.parse_args()
    ctx.prec=128
    ctx.threads=1
    d=json.loads(Path(args.centre).read_text())
    rows,g,c,pcoef,V,F=v.make_algebra(d)
    N=len(rows)
    v.need(0<=args.start<args.stop<=N,'invalid finite batch')
    rho=arb(11)/10
    A=load_inverse(args.inverse)
    inv=v.Inverse(d['M'],d['S'],rho,c[1],A)
    ng=len(g)
    js=sorted(c)
    cols=[]
    for i in range(args.start,args.stop):
        if i<ng:
            key=rows[i]
            e={key:arb(1)}
            cols.append(v.add(e,v.unsine(v.abs2(v.sine(v.K(e)),pcoef))))
        else:
            cols.append(v.shape_column(js[i-ng],c,pcoef,V))
    mat=arb_mat([[col.get(key,v.Z) for col in cols] for key in rows])
    image=A*mat
    defect=[]
    for j,i in enumerate(range(args.start,args.stop)):
        column=sum((abs((arb(1) if k==i else v.Z)-image[k,j])*inv.cw[k]
                    for k in range(N)),v.Z)/inv.cw[i]
        defect.append(column)
    zvals=inv.apply_many(cols,subtract=[('finite',i) for i in range(args.start,args.stop)])
    zrel=[x/inv.cw[i] for i,x in zip(range(args.start,args.stop),zvals)]
    out={'kind':'finite_batch','dimension':2*v.M_ROOT+2,
         'start':args.start,'stop':args.stop,'N':N,
         'inverse_defect_max':v.max_upper(defect).str(30),
         'Z_finite_max':v.max_upper(zrel).str(30)}
    Path(args.output).write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))


if __name__=='__main__':
    main()
