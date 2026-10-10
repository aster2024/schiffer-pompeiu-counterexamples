#!/usr/bin/env python3
"""Exact rational verification of V=K_D g=(1-r²)Q for the frozen centre."""
from fractions import Fraction as F
import json

with open('finite_center_r4_M31_S24.json') as f:data=json.load(f)
M,S=data['M'],data['S'];assert (M,S)==(31,24)
ms=list(range(1,M+1,2))
g={(m,s):F.from_float(float.fromhex(v))
   for (m,s),v in zip(((m,s) for m in ms for s in range(S+1)),data['g_hex'])}
V={}
def acc(k,v):V[k]=V.get(k,F(0))+v
for (m,s),v in g.items():
    if s==0:
        t=v/F(4*(m+1)*(m+2))
        acc((m,0),-t);acc((m,1),t)
    else:
        D=m+2*s
        acc((m,s-1),v/F(4*D*(D+1)))
        acc((m,s),-v/F(2*D*(D+2)))
        acc((m,s+1),v/F(4*(D+1)*(D+2)))

def up(m,s):return F((s+1)*(m+s+1),(m+2*s+1)*(m+2*s+2))
def down(m,s):return F(0) if s==0 else F(s*(m+s),(m+2*s)*(m+2*s+1))

Q={}
for m in ms:
    qs={S+1:F(0),S+2:F(0)}
    for s in range(S,-1,-1):
        qs[s]=((up(m,s+1)+down(m,s+1))*qs[s+1]
               -down(m,s+2)*qs[s+2]-V.get((m,s+1),F(0)))/up(m,s)
        Q[m,s]=qs[s]
    for s in range(S+2):
        got=(up(m,s)+down(m,s))*qs.get(s,F(0))
        if s:got-=up(m,s-1)*qs[s-1]
        got-=down(m,s+1)*qs.get(s+1,F(0))
        assert got==V.get((m,s),F(0)),(m,s)
print('EXACT DIRICHLET FACTOR VERIFIED',len(Q),'rational coefficients')
