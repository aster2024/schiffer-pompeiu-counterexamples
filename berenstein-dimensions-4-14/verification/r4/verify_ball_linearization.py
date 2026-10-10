#!/usr/bin/env python3
"""Arb enclosure of one Dirichlet-ball near-resonance in R4.

This verifies Bessel near-resonance and U4 moment inequalities.
It is the consistency check described in Section 8.
"""
from flint import arb, ctx
from fractions import Fraction
from math import comb

ctx.prec=192
ctx.threads=1
U4={0:1,1:-12,2:16} # polynomial in t=cos(theta)^2
def mul(a,b):
    out={}
    for i,x in a.items():
        for j,y in b.items():out[i+j]=out.get(i+j,0)+x*y
    return out
def average(a):
    return sum((Fraction(v*comb(2*k,k),(k+1)*4**k) for k,v in a.items()),Fraction(0))
assert average(mul(U4,U4))==1
assert average(mul(mul(U4,U4),U4))==1
lo=arb(701558666980)/arb(10)**11
hi=arb(701558666984)/arb(10)**11
assert lo.bessel_j(1).upper()<0
assert hi.bessel_j(1).lower()>0
R=(lo+hi)/2+arb(0,2.1e-11)
assert (R-lo).lower()<0 and (R-hi).upper()>0
J1p=(R.bessel_j(0)-R.bessel_j(2))/2
assert J1p.lower()>0
J5=R.bessel_j(5)
J5p=(R.bessel_j(4)-R.bessel_j(6))/2
alpha=R*J5p+2*J5
lam=alpha/J5
assert J5.lower()>0 and alpha.lower()>0 and lam.lower()>0
J5pp=-(R*J5p+(R*R-25)*J5)/(R*R)
beta_bracket=J5+2*R*J5p+R*R*J5pp/2+arb(24)*J5/4
predicted=-5*alpha/beta_bracket
assert beta_bracket.upper()<arb(-27504)/10000
assert predicted.lower()>arb(716)/10000 and predicted.upper()<arb(717)/10000
print('CERTIFIED BALL LINEARIZATION')
print('J1(lo)',lo.bessel_j(1))
print('J1(hi)',hi.bessel_j(1))
print('alpha4',alpha)
print('lambda4',lam)
print('quadratic bracket',beta_bracket,'cubic moment ratio',arb(1)/5)
print('quadratic root prediction',predicted)
