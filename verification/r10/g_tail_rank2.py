"""Arb near g-tail columns and analytic far g-tail bounds for R10/R14."""

import argparse
import json
from pathlib import Path
from flint import arb, ctx
from finite_rank2 import load_inverse
import rank2_cap_core as v


def columns(M,S,degree):
    return [(n,s) for n in range(v.M_ROOT,M+degree+1,v.STEP)
            for s in range(1,S+degree+2)
            if not(n<=M and s<=S)]


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--centre',required=True)
    p.add_argument('--inverse',required=True)
    p.add_argument('--finite-receipt',required=True)
    p.add_argument('--tail-receipt',required=True)
    p.add_argument('--start',type=int,default=0)
    p.add_argument('--stop',type=int,default=0)
    p.add_argument('--short-degree',type=int,default=24)
    p.add_argument('--output',required=True)
    args=p.parse_args()
    ctx.prec=128
    ctx.threads=1
    d=json.loads(Path(args.centre).read_text())
    rows,g,c,pcoef,V,F=v.make_algebra(d)
    M,S=d['M'],d['S']
    degree=args.short_degree
    v.need(degree in pcoef and degree%v.STEP==0,'unsupported short degree')
    rho=arb(11)/10
    A=load_inverse(args.inverse)
    inv=v.Inverse(M,S,rho,c[1],A)
    finite=json.loads(Path(args.finite_receipt).read_text())
    tail=json.loads(Path(args.tail_receipt).read_text())
    An=arb(finite['A_finite']).upper()
    Atail=arb(tail['A_tail_upper']).upper()
    keys=columns(M,S,degree)
    ps={a:x for a,x in pcoef.items() if a<=degree}
    Ps=sum((abs(x)*rho**a for a,x in ps.items()),arb(0))
    Pt=sum((abs(x)*rho**a for a,x in pcoef.items() if a>degree),arb(0))
    delta=2*Ps*Pt+Pt*Pt
    n_far=M+degree+v.STEP
    d_far=n_far+2
    far_n=(Atail*Ps*Ps+An*delta)/(d_far*(d_far+2))
    s_far=S+degree+2
    d_s=v.M_ROOT+2*s_far
    far_s=(Ps*Ps+An*delta)/(d_s*(d_s+2))
    if args.stop==0:
        out={'kind':'far','dimension':2*v.M_ROOT+2,'degree':degree,
             'near_count':len(keys),'n_far':n_far,'s_far':s_far,
             'far_n':far_n.upper().str(30),'far_s':far_s.upper().str(30),
             'Ps':Ps.str(30),'Pt':Pt.str(30),'delta':delta.str(30)}
    else:
        v.need(0<=args.start<args.stop<=len(keys),'invalid g-tail batch')
        selected=keys[args.start:args.stop]
        fs=[v.unsine(v.abs2(v.sine(v.K({key:arb(1)})),ps)) for key in selected]
        vals=inv.apply_many(fs)
        rel=[val/rho**n+An*delta/((n+2*s)*(n+2*s+2))
             for (n,s),val in zip(selected,vals)]
        winner=v.argmax_upper(rel)
        out={'kind':'near','dimension':2*v.M_ROOT+2,
             'start':args.start,'stop':args.stop,
             'max':v.max_upper(rel).str(30),'arg':list(selected[winner]),
             'near_count':len(keys)}
    Path(args.output).write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))


if __name__=='__main__':
    main()
