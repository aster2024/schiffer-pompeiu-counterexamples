"""Arb geometry and sign evaluations at a finite rank-two centre.

The radius-ball bounds are checked in verify_radii_geometry_r14.py, Section 6.6."""
from proof_guard import require
import argparse
from fractions import Fraction
from flint import arb,ctx
from exact_finite_rank2 import center,KD
from exact_finite_r4 import ZERO


def jacobi_values(n,S,x):
    P=[arb(1)]
    if S:
        P.append(1+(n+2)*(x-1)/2)
    for s in range(2,S+1):
        c1=2*s*(s+n)*(2*s+n-2)
        c2=(2*s+n-1)*((2*s+n)*(2*s+n-2)*x-n*n)
        c3=2*(s-1)*(s+n-1)*(2*s+n)
        P.append((c2*P[-1]-c3*P[-2])/c1)
    return P


def value(V,ms,m,S,r):
    x=2*r*r-1
    total=arb(0)
    for n in ms:
        P=jacobi_values(n,S+1,x)
        sign=-1 if ((n//m-1)//2)%2 else 1
        total+=sign*r**n*sum((V.get((n,s),ZERO)*P[s] for s in range(S+2)),ZERO)
    return total


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('center')
    ap.add_argument('--point',action='append',default=[])
    ap.add_argument('--bits',type=int,default=192)
    args=ap.parse_args()
    ctx.prec=args.bits;ctx.threads=1
    d,g,c,B,keys,ms,js=center(args.center)
    m=d['m'];b=c[1]
    L=b-sum((j*abs(v) for j,v in c.items() if j>1),ZERO)
    U=sum((j*(j-1)*abs(v) for j,v in c.items() if j>1),ZERO)
    require((L.lower()>0 and (L-U).lower()>0), 'verify_center_geometry_sign_rank2.py:46: proof gate failed')
    outer=arb(1001)/1000
    Louter=b-sum((j*outer**(j-1)*abs(v) for j,v in c.items() if j>1),ZERO)
    require((Louter.lower()>0), 'verify_center_geometry_sign_rank2.py:49: proof gate failed')
    nonball=arb(0)
    for j in js[1:]:
        if abs(c[j]).lower()>0:
            nonball=abs(c[j]);break
    require((nonball.lower()>0), 'verify_center_geometry_sign_rank2.py:54: proof gate failed')
    V=KD(g,m)
    print('CENTER_GEOMETRY_ONLY m',m,'dimension',2*m+2,
          'b',b,'Re p lower',L.lower(),
          'convexity lower',(1-U/L).lower(),
          'outer Re p lower',Louter.lower(),
          'nonball coefficient lower',nonball.lower())
    for point in args.point:
        q=Fraction(point)
        r=arb(q.numerator)/q.denominator
        require((0<r<1), 'verify_center_geometry_sign_rank2.py:64: proof gate failed')
        v=value(V,ms,m,d['S'],r)
        require((v.is_finite()), 'verify_center_geometry_sign_rank2.py:66: proof gate failed')
        print('V at w radius',point,'theta=pi/(2m):',v,
              'lower',v.lower(),'upper',v.upper())


if __name__=='__main__':
    main()
