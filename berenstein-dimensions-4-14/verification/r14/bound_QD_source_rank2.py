"""Arb uniform bound for Q |p|^2 K_D on the high-source basis.

Uses positive Zernike multiplication, support of the harmonic row, and the
weighted character multiplier bound. This is one derivative block only.
"""
from proof_guard import require
import argparse
from flint import arb,ctx
from exact_finite_rank2 import center,q_of,to_char
from exact_finite_r4 import ZERO


def upper(x):
    require((x.is_finite()), 'bound_QD_source_rank2.py:13: proof gate failed')
    y=x.upper();require((y.is_finite()), 'bound_QD_source_rank2.py:14: proof gate failed')
    return y


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('center')
    ap.add_argument('--bits',type=int,default=160)
    a=ap.parse_args()
    ctx.prec=a.bits;ctx.threads=1
    d,g,c,B,keys,ms,js=center(a.center)
    m,M,S=d['m'],d['M'],d['S'];rho=arb(102)/100;b=c[1]
    p={j-1:j*v for j,v in c.items()}
    bydeg={}
    for ea,va in p.items():
        for eb,vb in p.items():
            degree=ea+eb
            bydeg[degree]=bydeg.get(degree,ZERO)+abs(va)*abs(vb)*rho**degree
    dmax=max(bydeg)
    qc=to_char(q_of(g,m),m)
    Cq=sum(((n//m)*rho**(n-m)*abs(v) for n,v in qc.items()),ZERO)
    factor=Cq/(b*rho**m)
    nlast=m
    while nlast<M+dmax+2*m:nlast+=2*m
    best=arb(0);arg=None;count=0
    for n in range(m,nlast+1,2*m):
        for s in range(dmax+2):
            if n<=M and s<=S:continue
            if s==0:
                z=arb(1)/(4*(n+1)*(n+2))
                Kterms=((0,z),(1,z))
            else:
                D=n+2*s
                Kterms=((s-1,arb(1)/(4*D*(D+1))),
                        (s,arb(1)/(2*D*(D+2))),
                        (s+1,arb(1)/(4*(D+1)*(D+2))))
            val=arb(0)
            for degree,weight in bydeg.items():
                kco=sum((v for t,v in Kterms if t<=degree),ZERO)
                if kco.is_zero():continue
                Nmin=max(m,n-degree)
                val+=weight*kco/(Nmin+1)
            val=upper(factor*val)
            if val>best:best=val;arg=(n,s)
            count+=1
    # For n>nlast, both kappa(n,s) and 1/(max(m,n-degree)+1)
    # decrease; s>dmax+1 cannot feed the harmonic row at all.
    require((best<arb(114)/1000000), 'bound_QD_source_rank2.py:61: proof gate failed')
    print('QD_SOURCE_BOUND_ONLY','m',m,'checked',count,
          'nlast',nlast,'dmax',dmax,'upper',best,'argmax',arg)


if __name__=='__main__':
    main()
