"""Generic Arb principal-tail inverse bound for R10/R14."""

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
    p.add_argument('--output',required=True)
    p.add_argument('--near-count',type=int,default=100)
    args=p.parse_args()
    ctx.prec=128
    ctx.threads=1
    d=json.loads(Path(args.centre).read_text())
    rows,g,c,pcoef,V,F=v.make_algebra(d)
    M,S=d['M'],d['S']
    rho=arb(11)/10
    A=load_inverse(args.inverse)
    inv=v.Inverse(M,S,rho,c[1],A)
    q=inv.radial_coeff(M)
    vec=[v.Z for _ in rows]
    vec[inv.pos[M,0]]=inv.alpha*(M+1)
    for s in range(1,v.M_ROOT+1):
        vec[inv.pos[M,s]]=inv.alpha*(M+v.M_ROOT+1)*q[s]
    image=A*arb_mat([[x] for x in vec])
    H=sum((abs(image[i,0])*inv.cw[i] for i in range(inv.N)),v.Z)
    n0=M+v.STEP*(args.near_count+1)
    ns=list(range(M+v.STEP,n0,v.STEP))
    vals=inv.apply_many([{(n,0):arb(1)} for n in ns])
    ratios=[val/rho**n for n,val in zip(ns,vals)]
    tau=1/(1-rho**-v.STEP)
    far=tau+v.M_ROOT*tau*rho**-v.STEP/(n0+1)+H/(inv.alpha*(n0+1)*rho**n0)
    bound=v.max_upper([arb(1),far]+ratios)
    imax=v.argmax_upper(ratios)
    out={'status':'OPEN_COMPONENT','dimension':2*v.M_ROOT+2,'M':M,'S':S,
         'count':len(ns),'near_max':v.max_upper(ratios).str(30),
         'near_arg':ns[imax],'far_start':n0,'far_upper':far.upper().str(30),
         'A_tail_upper':bound.str(30),'H':H.str(30)}
    Path(args.output).write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))


if __name__=='__main__':
    main()
