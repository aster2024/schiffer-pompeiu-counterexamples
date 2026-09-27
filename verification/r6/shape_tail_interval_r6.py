"""Exact Arb columns for the near part of the R6 conformal shape tail."""
import argparse
import json
from flint import arb, ctx
from finite_interval_r6 import load_inverse
from verify_r6 import (Z, Inverse, abs2, add, hol, hol_psi, make_algebra,
                       max_upper, scale, shift, sine, unsine)


def principal(j,b):
    b3=b**3
    return {(j+1,0):b3*(j+2)/2,
            (j-3,0):-b3*(j-2)/2,
            (j-3,1):-b3*j/(j+1),
            (j-3,2):-b3/(j+1)}


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--first',type=int,required=True)
    ap.add_argument('--last',type=int,required=True)
    ap.add_argument('--output')
    args=ap.parse_args()
    ctx.prec=128;ctx.threads=1
    d=json.load(open('center_r6.json'))
    rows,g,c,p,V,F=make_algebra(d)
    M,S=d['M'],d['S'];rho=arb(11)/10;b=c[1]
    assert args.first>=M+3 and args.first%4==1 and args.last>=args.first and args.last%4==1
    Ahat=load_inverse('inverse_r6_M58_S30.npy',len(rows))
    inv=Inverse(M,S,rho,b,Ahat)
    A=shift(hol(sine(V),p,True),args.first-1)
    B=shift(hol_psi(abs2({(1,0):arb(1)},p),c),args.first-1)
    cols=[];js=[]
    for j in range(args.first,args.last+1,4):
        DF=scale(unsine(add(scale(A,arb(j)),scale(B,-arb(1)/2))),arb(2))
        Delta=add(DF,scale(principal(j,b),-arb(1)))
        cols.append(Delta);js.append(j)
        A=shift(A,4);B=shift(B,4)
    vals=inv.apply_many(cols)
    rel=[v/inv.omega(j) for v,j in zip(vals,js)]
    imax=max(range(len(js)),key=lambda i:rel[i].upper())
    out={'first':args.first,'last':args.last,'count':len(js),'max':max_upper(rel).str(30),
         'arg':js[imax], 'last_value':rel[-1].upper().str(30)}
    print(json.dumps(out,indent=2),flush=True)
    if args.output:
        with open(args.output,'w') as f:json.dump(out,f,indent=2)


if __name__=='__main__':main()
