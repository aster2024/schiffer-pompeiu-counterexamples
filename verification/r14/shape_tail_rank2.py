"""Arb near-shape columns for scaled R10/R14 CAPs."""

import argparse
import json
from pathlib import Path
from flint import arb, ctx
from finite_rank2 import load_inverse
import rank2_cap_core as v


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--centre',required=True)
    p.add_argument('--inverse',required=True)
    p.add_argument('--first',type=int,required=True)
    p.add_argument('--last',type=int,required=True)
    p.add_argument('--output',required=True)
    args=p.parse_args()
    ctx.prec=128
    ctx.threads=1
    d=json.loads(Path(args.centre).read_text())
    rows,g,c,pcoef,V,F=v.make_algebra(d)
    M,S=d['M'],d['S']
    rho=arb(11)/10
    v.need(args.first>=M+v.M_ROOT+1 and args.first%v.STEP==1
           and args.last>=args.first and args.last%v.STEP==1,
           'invalid shape-tail range')
    A=load_inverse(args.inverse)
    inv=v.Inverse(M,S,rho,c[1],A)
    js=list(range(args.first,args.last+1,v.STEP))
    fs=[v.add(v.shape_column(j,c,pcoef,V),
              v.scale(v.principal_shape(j,c[1]),-arb(1))) for j in js]
    vals=inv.apply_many(fs)
    rel=[val/inv.omega(j) for val,j in zip(vals,js)]
    imax=v.argmax_upper(rel)
    out={'kind':'near','dimension':2*v.M_ROOT+2,
         'first':args.first,'last':args.last,'count':len(js),
         'max':v.max_upper(rel).str(30),'arg':js[imax],
         'last_value':rel[-1].upper().str(30)}
    Path(args.output).write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))


if __name__=='__main__':
    main()
