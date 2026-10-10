#!/usr/bin/env python3
"""Arb polynomial majorants for the R4 fixed-disc nonlinear remainder."""
from flint import arb,ctx
from exact_finite_r4 import (get_center,KD,ZERO,laur_p,laur_h,laur_a2,
                             conv,to_char)

ctx.prec=128;ctx.threads=1
d,g,c,keys,ms=get_center('finite_center_r4_M31_S24.json')
rho=arb(11)/10;b=c[1]
p=laur_p(c);h=laur_h(c);a2=laur_a2(p)
P=sum((rho**e*abs(v) for e,v in p.items()),ZERO)
P2=sum((rho**abs(e)*abs(v) for e,v in a2.items()),ZERO)
hchar=to_char(h)
H0=sum((rho**(j-1)*abs(v) for j,v in hchar.items()),ZERO)
H1=sum((j*rho**(j-1)*abs(v) for j,v in hchar.items()),ZERO)
h2=to_char(conv(h,h))
H20=sum((rho**(j-1)*abs(v) for j,v in h2.items()),ZERO)
V=KD(g);Vnorm=sum((rho**m*abs(v) for (m,s),v in V.items()),ZERO)
K=arb(1)/12
ap=arb(1)/b**3
ah0=arb(1)/(3*b**3);ah1=arb(1)/b**3
aq0=arb(1)/(4*rho);aq1=arb(1)/(2*rho)
Lh=2*min((H1*ah0).upper(),(H0*ah1).upper())
Qh=ah0*ah1
di2=2*P*ap*K+Vnorm*ap**2
di3=K*ap**2
db2=aq0*aq1+H20*ap**2+Lh*(2*P*ap)+Qh*P2
db3=Lh*ap**2+Qh*(2*P*ap)
db4=Qh*ap**2

normAf=arb(1665434)/1000 # certified <1665.434
edge_col=arb(54778)/1000 # certified <54.778
Cq=sum((j*rho**(j-1)*abs(g[(j,0)])/(2*(j+1)) for j in ms),ZERO)
Lt=rho*rho/(rho*rho-1)
edge=edge_col*(ms[-1]+2)/(ms[-1]+4)/rho**(ms[-1]+1)
normA=normAf+1+(1+edge)*Lt*Cq/(2*rho)
assert normA.upper()<arb(2444)
C2=normA*(di2+db2)
C3=normA*(di3+db3)
C4=normA*db4
assert C2.upper()<arb(268)
assert C3.upper()<arb(23)/10000
assert C4.upper()<arb(6)/10**8
print('P',P,'P2',P2,'H0',H0,'H1',H1,'H20',H20,'Vnorm',Vnorm)
print('normA upper',normA.upper())
print('C2 upper',C2.upper())
print('C3 upper',C3.upper())
print('C4 upper',C4.upper())
