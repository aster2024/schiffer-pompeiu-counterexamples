#!/usr/bin/env python3
"""Exact rational radii inequalities from outward-rounded Arb majorants."""
from fractions import Fraction as Q

Y=Q(3038,10**9)   # 3.038e-6
Z=Q(463,500)       # 0.926
C2=Q(268)
C3=Q(23,10000)
C4=Q(6,10**8)
r=Q(1,12500)
p=Y+(Z-1)*r+C2*r*r+C3*r**3+C4*r**4
dp=Z-1+2*C2*r+3*C3*r*r+4*C4*r**3
assert p<0 and dp<0
print('RADII INEQUALITIES PASS')
print('radius',float(r),'polynomial',float(p),'derivative',float(dp))
