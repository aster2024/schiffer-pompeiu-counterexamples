"""Rigorous near-g-tail columns and analytic far-g-tail estimates."""
import argparse
import json
from flint import arb, ctx
from finite_interval_r6 import load_inverse
from verify_r6 import Inverse, abs2, K, make_algebra, max_upper, norm, sine, unsine


def columns(M,S,degree):
    return [(n,s) for n in range(2,M+degree+1,4)
            for s in range(1,S+degree+2)
            if not(n<=M and s<=S)]


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--start',type=int,default=0)
    ap.add_argument('--stop',type=int,default=0)
    ap.add_argument('--output')
    ap.add_argument('--short-degree',type=int,default=24)
    args=ap.parse_args()
    ctx.prec=128;ctx.threads=1
    d=json.load(open('center_r6.json'))
    rows,g,c,p,V,F=make_algebra(d)
    M,S=d['M'],d['S'];degree=args.short_degree
    assert degree in p
    rho=arb(11)/10;b=c[1]
    Ahat=load_inverse('inverse_r6_M58_S30.npy',len(rows))
    inv=Inverse(M,S,rho,b,Ahat)
    ks=columns(M,S,degree)
    ps={a:v for a,v in p.items() if a<=degree}
    Ps=sum((abs(v)*rho**a for a,v in ps.items()),arb(0))
    Pt=sum((abs(v)*rho**a for a,v in p.items() if a>degree),arb(0))
    delta=2*Ps*Pt+Pt*Pt
    An=arb('1402.314')
    Atail=arb('3.164016')
    far_n=(Atail*Ps*Ps+An*delta)/((M+degree+4+2)*(M+degree+4+4))
    s_far=S+degree+2
    far_s=(Ps*Ps+An*delta)/((2+2*s_far)*(4+2*s_far))
    if args.stop==0:
        out={'kind':'far','degree':degree,'near_count':len(ks),
             'far_n':far_n.upper().str(30),'far_s':far_s.upper().str(30),
             'Ps':Ps.str(30),'Pt':Pt.str(30),'delta':delta.str(30)}
    else:
        assert 0<=args.start<args.stop<=len(ks)
        keyset=ks[args.start:args.stop]
        fs=[unsine(abs2(sine(K({key:arb(1)})),ps)) for key in keyset]
        vals=inv.apply_many(fs)
        rel=[v/rho**n+An*delta/((n+2*s)*(n+2*s+2)) for (n,s),v in zip(keyset,vals)]
        winner=max(range(len(rel)),key=lambda i:rel[i].upper())
        out={'kind':'near','start':args.start,'stop':args.stop,
             'max':max_upper(rel).str(30),'arg':list(keyset[winner]),
             'near_count':len(ks)}
    print(json.dumps(out,indent=2),flush=True)
    if args.output:
        with open(args.output,'w') as f:json.dump(out,f,indent=2)


if __name__=='__main__':main()
