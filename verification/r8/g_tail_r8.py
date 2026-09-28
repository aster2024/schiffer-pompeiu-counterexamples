"""Arb near columns and analytic far columns for the R8 g tail."""

import argparse
import json
from flint import arb, ctx
from tail_inverse_r8 import load_inverse
import verify_r8 as v


def columns(M,S,degree):
    return [(n,s) for n in range(v.M_ROOT,M+degree+1,v.STEP)
            for s in range(1,S+degree+2)
            if not(n<=M and s<=S)]


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--start',type=int,default=0)
    p.add_argument('--stop',type=int,default=0)
    p.add_argument('--short-degree',type=int,default=24)
    p.add_argument('--output')
    args=p.parse_args()
    ctx.prec=128
    ctx.threads=1
    d=json.load(open('center_r8_M45_S24.json'))
    rows,g,c,pcoef,V,F=v.make_algebra(d)
    M,S=d['M'],d['S']
    degree=args.short_degree
    v.need(degree in pcoef and degree%v.STEP==0,'unsupported short degree')
    rho=arb(11)/10
    A=load_inverse('inverse_r8_M45_S24.npy')
    inv=v.Inverse(M,S,rho,c[1],A)
    keys=columns(M,S,degree)
    ps={a:x for a,x in pcoef.items() if a<=degree}
    Ps=sum((abs(x)*rho**a for a,x in ps.items()),arb(0))
    Pt=sum((abs(x)*rho**a for a,x in pcoef.items() if a>degree),arb(0))
    delta=2*Ps*Pt+Pt*Pt
    An=arb(964)
    Atail=arb('2.303')
    n_far=M+degree+v.STEP
    d_far=n_far+2
    far_n=(Atail*Ps*Ps+An*delta)/(d_far*(d_far+2))
    s_far=S+degree+2
    d_s=v.M_ROOT+2*s_far
    far_s=(Ps*Ps+An*delta)/(d_s*(d_s+2))
    if args.stop==0:
        out={'kind':'far','degree':degree,'near_count':len(keys),
             'far_n':far_n.upper().str(30),'far_s':far_s.upper().str(30),
             'Ps':Ps.str(30),'Pt':Pt.str(30),'delta':delta.str(30),
             'n_far':n_far,'s_far':s_far}
    else:
        v.need(0<=args.start<args.stop<=len(keys),'invalid g-tail batch')
        selected=keys[args.start:args.stop]
        fs=[v.unsine(v.abs2(v.sine(v.K({key:arb(1)})),ps)) for key in selected]
        vals=inv.apply_many(fs)
        rel=[val/rho**n+An*delta/((n+2*s)*(n+2*s+2))
             for (n,s),val in zip(selected,vals)]
        winner=v.argmax_upper(rel)
        out={'kind':'near','start':args.start,'stop':args.stop,
             'max':v.max_upper(rel).str(30),'arg':list(selected[winner]),
             'near_count':len(keys)}
    print(json.dumps(out,indent=2),flush=True)
    if args.output:
        with open(args.output,'w') as f:
            json.dump(out,f,indent=2)
            f.write('\n')


if __name__=='__main__':
    main()
