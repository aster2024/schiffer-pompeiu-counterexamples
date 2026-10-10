#!/usr/bin/env python3
"""Arb signs of the exact R4 solution at two interior meridian points."""
from flint import arb,ctx
from exact_finite_r4 import get_center,KD,ZERO

ctx.prec=192;ctx.threads=1
_,g,c,keys,ms=get_center('finite_center_r4_M31_S24.json')
V=KD(g)


def jacobi(m,smax,x):
    p=[arb(1),arb(1)+(m+2)*(x-1)/2]
    for n in range(2,smax+1):
        a,b=arb(0),arb(m)
        c1=2*n*(n+a+b)*(2*n+a+b-2)
        c2=(2*n+a+b-1)*((2*n+a+b)*(2*n+a+b-2)*x+a*a-b*b)
        c3=2*(n+a-1)*(n+b-1)*(2*n+a+b)
        p.append((c2*p[-1]-c3*p[-2])/c1)
    return p


def value(r):
    x=2*r*r-1;out=arb(0)
    for m in ms:
        ps=jacobi(m,25,x)
        sign=1 if ((m-1)//2)%2==0 else -1
        out+=sum((V.get((m,s),ZERO)*r**m*ps[s]*sign for s in range(26)),ZERO)
    return out


v1=value(arb(1)/4);v2=value(arb(3)/4)
error=(arb(1)/12)*(arb(1)/12500) # ||K_D||*proof-ball radius
assert (v1-error).lower()>10 and (v2+error).upper()<-8
print('CERTIFIED SIGN CHANGE')
print('V(i/4)',v1,'V(3i/4)',v2,'uniform perturbation',error)
