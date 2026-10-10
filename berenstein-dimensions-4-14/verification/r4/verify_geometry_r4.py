#!/usr/bin/env python3
"""Arb geometry and sign margins for every conformal map in the proof ball.

This is conditional on an exact zero being certified in that ball.
"""
from flint import arb,ctx
from exact_finite_r4 import get_center,ZERO

ctx.prec=128;ctx.threads=1
d,g,c,keys,ms=get_center('finite_center_r4_M31_S24.json')
b=c[1];rho=arb(11)/10;R=arb(21)/20;r=arb(1)/12500
tail_derivative=r/b**3
tail_second=(r/b**3)*rho/(arb(1).exp()*rho.log())
deriv_R=sum((j*abs(v)*R**(j-1) for j,v in c.items() if j>=3),ZERO)
deriv_1=sum((j*abs(v) for j,v in c.items() if j>=3),ZERO)
second_1=sum((j*(j-1)*abs(v) for j,v in c.items() if j>=3),ZERO)
lower_univalent=b-deriv_R-tail_derivative
lower_p=b-deriv_1-tail_derivative
upper_pp=second_1+tail_second
lower_convex=1-upper_pp/lower_p
s0=sum((abs(v) for j,v in c.items() if j>=3),ZERO)+r/(5*b**3*rho**2)
s1=sum(((j-1)*abs(v) for j,v in c.items() if j>=3),ZERO)+tail_derivative
denom=b-r/(3*b**3)-s0
star_margin=1-s1/denom
c3_margin=abs(c[3])-r/(5*b**3*rho**2)
q_margin=(g[(1,0)]/4
          -sum((m*abs(g[(m,0)])/(2*(m+1)) for m in ms if m>=3),ZERO)
          -r/2)
for name,x in [('univalence derivative',lower_univalent),('p modulus',lower_p),
               ('convexity criterion',lower_convex),('star-shaped criterion',star_margin),
               ('nonball c3',c3_margin),('q sign',q_margin)]:
    if not x.is_finite() or not (x.lower()>0):
        raise RuntimeError(f'{name} failed: {x}')
    print(name,'lower',x.lower())
print('CERTIFIED GEOMETRY MARGINS CONDITIONAL ON EXACT ZERO')
