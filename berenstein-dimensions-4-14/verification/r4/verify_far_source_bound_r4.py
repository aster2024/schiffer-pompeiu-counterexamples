#!/usr/bin/env python3
"""Arb expressions for the far-source bounds in Section 5.2.

These bound the unenumerated source columns of the approximate inverse."""
from flint import arb,ctx
from exact_finite_r4 import get_center,ZERO

ctx.prec=128;ctx.threads=1
d,g,c,keys,ms=get_center('finite_center_r4_M31_S24.json')
M,S=d['M'],d['S']; assert M==31 and S==24
rho=arb(11)/10;b=c[1]
P=sum((j*rho**(j-1)*abs(v) for j,v in c.items()),ZERO)
Cq=sum((j*rho**(j-1)*abs(g[(j,0)])/(2*(j+1)) for j in ms),ZERO)
L=rho*rho/(rho*rho-1)
edge_col=arb(54778)/1000 # checked in verify_finite_inverse_r4.py
edge=edge_col*(M+2)/(M+4)/rho**(M+1)
mstar=2*M+1 # 63
harmonic=P*P/(2*(mstar+1)*(mstar+2))
nonharmonic=P*P/((mstar+2)*(mstar+4))
if not (harmonic.upper()<nonharmonic.lower()):
    raise RuntimeError('far-source maximum case changed')
d_ratio=nonharmonic
# KD input of degree m>=63 retains angular degree m. Multiplication by
# |p|² shifts angular degree by at most M-1=30. The quotient trace Q then
# costs at most 1/[2 rho (m-M+2)]. Multiplication by q0 costs Cq in the
# unnormalised U-character algebra. The tail shape inverse costs L.
far_angular=d_ratio*(1+(1+edge)*L*Cq/(rho*(mstar-M+2)))
D=1+2*(S+32) # m=1,s=56: D=113
far_radial=P*P/(D*(D+2))
assert far_angular.upper()<arb(9)/10
assert far_radial.upper()<arb(7)/1000
print('CERTIFIED FAR SOURCE BOUNDS ONLY')
print('P upper',P.upper(),'Cq upper',Cq.upper())
print('far angular m>=63 upper',far_angular.upper())
print('far radial s>=56 upper',far_radial.upper())
