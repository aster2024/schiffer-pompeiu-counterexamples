"""Arb majorant for the far high-shape Schur perturbation.

This checks the infinite column class estimated in Section 6.2.
"""
from proof_guard import require
import argparse
import hashlib
import json
from flint import arb,ctx
from exact_finite_rank2 import center,KD,h_of,q_of,to_char,sine
from exact_finite_r4 import ZERO,poly,conv


def up(n,s):
    return arb((s+1)*(n+s+1))/((n+2*s+1)*(n+2*s+2))


def down(n,s):
    return ZERO if s==0 else arb(s*(n+s))/((n+2*s)*(n+2*s+1))


def upper(x):
    require((x.is_finite()), 'far_schur_bound_rank2.py:23: proof gate failed')
    y=x.upper()
    require((y.is_finite()), 'far_schur_bound_rank2.py:25: proof gate failed')
    return y


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('center')
    ap.add_argument('--J',type=int,required=True)
    ap.add_argument('--rho-num',type=int,default=11)
    ap.add_argument('--rho-den',type=int,default=10)
    ap.add_argument('--bits',type=int,default=256)
    ap.add_argument('--output',default='')
    args=ap.parse_args()
    ctx.prec=args.bits;ctx.threads=1
    d,g,c,B,keys,ms,js=center(args.center)
    m=d['m'];J=args.J;require((J>js[-1] and (J-1)%(2*m)==0), 'far_schur_bound_rank2.py:40: proof gate failed')
    rho=arb(args.rho_num)/args.rho_den
    b=c[1]
    factor=arb(1)/2 if m==2 else arb(1)
    hb=factor*b**m/B
    alpha=hb*hb
    p={j-1:j*v for j,v in c.items()}
    P=sum((j*rho**(j-1)*abs(v) for j,v in c.items()),ZERO)
    dP=sum((j*rho**(j-1)*abs(v) for j,v in c.items() if j>1),ZERO)
    C=sum((rho**(j-1)*abs(v) for j,v in c.items()),ZERO)
    h=h_of(c,m,B)
    hc=to_char(h,m)
    delta_h1=sum(((n//m)*rho**(n-m)*abs(v-(hb if n==m else ZERO))
                  for n,v in hc.items()),ZERO)
    hh=to_char(conv(h,h),m)
    delta_h2=sum((rho**(n-m)*abs(v-(alpha if n==m else ZERO))
                  for n,v in hh.items()),ZERO)
    dh=factor*m*C**(m-1)/B
    dhdiff=factor*m*abs(C**(m-1)-b**(m-1))/B
    dhball=m*hb/b
    a2diff=2*b*dP+dP*dP
    first=(arb(2)/(b*(J+2*m))*(delta_h1*dh*P*P+
          hb*dhdiff*P*P+hb*dhball*a2diff))
    second=arb(2)/b*(P*delta_h2+alpha*dP)
    boundary_delta=first+second

    V=KD(g,m);Q={};S=d['S']
    for n in ms:
        q={S+1:ZERO,S+2:ZERO}
        for s in range(S,-1,-1):
            q[s]=((up(n,s+1)+down(n,s+1))*q[s+1]
                  -down(n,s+2)*q[s+2]-V.get((n,s+1),ZERO))/up(n,s)
            Q[n,s]=q[s]
        for s in range(S+2):
            val=(up(n,s)+down(n,s))*q.get(s,ZERO)
            if s:val-=up(n,s-1)*q[s-1]
            val-=down(n,s+1)*q.get(s+1,ZERO)
            require(((val-V.get((n,s),ZERO)).contains(0)), 'far_schur_bound_rank2.py:77: proof gate failed')
    W=poly(sine(Q,m),p,True)
    kmin=min(n for n,s in W)
    require((J+kmin>0), 'far_schur_bound_rank2.py:80: proof gate failed')
    Hbound=arb(4)*sum((abs(v)*rho**n*
              (2*(s+max(0,-n))+1)/(J+n)
              for (n,s),v in W.items()),ZERO)
    q0=q_of(g,m);qchar=to_char(q0,m)
    Cq=sum(((n//m)*rho**(n-m)*abs(v) for n,v in qchar.items()),ZERO)
    tau=1/(1-rho**(-2*m))
    shapeK=tau/alpha*boundary_delta
    traceK=tau/alpha*(arb(4)*Cq/(b*rho**m))*sum(
        (abs(v)*rho**n/((J+n)*(J+n+1)) for (n,s),v in W.items()),ZERO)
    total=shapeK+traceK
    print('m',m,'J',J,'rho',rho,'kmin',kmin,
          'P',upper(P),'dP',upper(dP),'delta_h1',upper(delta_h1),
          'delta_h2',upper(delta_h2))
    print('Hbound',upper(Hbound),'Cq',upper(Cq),'tau',upper(tau),
          'shapeK',upper(shapeK),'traceK',upper(traceK),
          'far Schur K upper',upper(total))
    require((upper(total)<arb(1)), 'far_schur_bound_rank2.py:97: proof gate failed')
    if args.output:
        value=upper(total)
        man,exp=value.man_exp()
        hman,hexp=upper(Hbound).man_exp()
        payload=dict(center=args.center,
                     center_sha256=hashlib.sha256(open(args.center,'rb').read()).hexdigest(),
                     J=J,rho_num=args.rho_num,rho_den=args.rho_den,
                     bits=args.bits,m=m,kmin=kmin,
                     Hbound_upper_dyadic=[int(hman),int(hexp)],
                     bound_upper_dyadic=[int(man),int(exp)])
        with open(args.output,'w') as f:json.dump(payload,f,separators=(',',':'))
        print('wrote',args.output)
    print('CERTIFIED FAR SCHUR BOUND ONLY')


if __name__=='__main__':
    main()
