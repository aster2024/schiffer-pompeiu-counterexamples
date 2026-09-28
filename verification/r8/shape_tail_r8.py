"""Arb near-shape columns of the R8 CAP, with exact principal subtraction."""

import argparse
import json
from flint import arb, ctx
from tail_inverse_r8 import load_inverse
import verify_r8 as v


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--first',type=int,required=True)
    p.add_argument('--last',type=int,required=True)
    p.add_argument('--output')
    args=p.parse_args()
    ctx.prec=128
    ctx.threads=1
    d=json.load(open('center_r8_M45_S24.json'))
    rows,g,c,pcoef,V,F=v.make_algebra(d)
    M,S=d['M'],d['S']
    rho=arb(11)/10
    v.need(args.first>=M+v.M_ROOT+1 and args.first%v.STEP==1
           and args.last>=args.first and args.last%v.STEP==1,
           'invalid shape-tail range')
    A=load_inverse('inverse_r8_M45_S24.npy')
    inv=v.Inverse(M,S,rho,c[1],A)
    js=list(range(args.first,args.last+1,v.STEP))
    fs=[v.add(v.shape_column(j,c,pcoef,V),
              v.scale(v.principal_shape(j,c[1]),-arb(1))) for j in js]
    vals=inv.apply_many(fs)
    rel=[val/inv.omega(j) for val,j in zip(vals,js)]
    imax=v.argmax_upper(rel)
    out={'kind':'near','first':args.first,'last':args.last,'count':len(js),
         'max':v.max_upper(rel).str(30),'arg':js[imax],
         'last_value':rel[-1].upper().str(30)}
    print(json.dumps(out,indent=2))
    if args.output:
        with open(args.output,'w') as f:
            json.dump(out,f,indent=2)
            f.write('\n')


if __name__=='__main__':
    main()
