#!/usr/bin/env python3
"""Arb constants for the analytic far-shape defect bound j>=301.

The formulas used here are proved in Section 5.2. This certifies
one infinite column class; it does not by itself prove the PDE theorem.
"""
from flint import arb,ctx
from exact_finite_r4 import (get_center,KD,ZERO,laur_p,laur_h,laur_a2,
                             laur_char,conv,to_char,scale,add,poly,sine)

ctx.prec=512;ctx.threads=1
d,g,c,keys,ms=get_center('finite_center_r4_M31_S24.json')
M,S=d['M'],d['S']; assert (M,S)==(31,24)
rho=arb(11)/10;b=c[1];J=301
p=laur_p(c);h=laur_h(c);a2=laur_a2(p)


def boundary_shape(j):
    dp={j-1:arb(j)};dh=laur_char(j)
    da2=add(conv(dp,{-k:v for k,v in p.items()}),
            conv(p,{-k:v for k,v in dp.items()}))
    return to_char(add(scale(conv(conv(h,dh),a2),-2),
                       scale(conv(conv(h,h),da2),-1)))


# In the high-mode character fusion regime the boundary shape column is
# affine in j after aligning each output by its offset from j.
j0,j1=125,127
B0,B1=boundary_shape(j0),boundary_shape(j1)
offsets=sorted(set(k-j0 for k,v in B0.items() if not v.is_zero())
               |set(k-j1 for k,v in B1.items() if not v.is_zero()))
assert min(offsets)==-92 and max(offsets)==60
A={r:(B1.get(j1+r,ZERO)-B0.get(j0+r,ZERO))/2 for r in offsets}
B={r:B0.get(j0+r,ZERO)-j0*A[r] for r in offsets}
assert sum(A.values(),ZERO).is_zero()
check=boundary_shape(129)
assert all((check.get(129+r,ZERO)-(129*A[r]+B[r])).is_zero() for r in offsets)
assert all(v.is_zero() for k,v in check.items() if k-129 not in offsets)
check=boundary_shape(301)
assert all((check.get(301+r,ZERO)-(301*A[r]+B[r])).is_zero() for r in offsets)
assert all(v.is_zero() for k,v in check.items() if k-301 not in offsets)

D={r:B[r]-(r+2)*A[r] for r in offsets}
run=arb(0);limit={}
for r in reversed(offsets):
    run+=A[r]/b**3
    limit[r]=-run
limit[0]-=1
L0=sum((rho**r*abs(v) for r,v in limit.items()),ZERO)
factor=sum((abs(r)*rho**r*abs(limit[r]) for r in offsets),ZERO)
derr=sum((rho**v*sum((abs(D[r]) for r in offsets if r>=v),ZERO)/b**3
          for v in offsets),ZERO)
leak=(rho**(min(offsets)-2)/(1-rho**(-2))
      *sum((abs(x) for x in D.values()),ZERO)/b**3)
boundary_bound=(L0+factor/(J+2)
                +(1+max(offsets)/(J+2))*derr/(J+min(offsets)+2)
                +leak/(J+min(offsets)+2))


# The exact finite V=K_D g vanishes at r=1, so V=(1-r²)Q. The Jacobi
# recurrence for r² gives a descending triangular solve for Q.
V=KD(g);Q={}
def up(m,s):return arb((s+1)*(m+s+1))/((m+2*s+1)*(m+2*s+2))
def down(m,s):
    return ZERO if s==0 else arb(s*(m+s))/((m+2*s)*(m+2*s+1))
for m in ms:
    qs={S+1:ZERO,S+2:ZERO}
    for s in range(S,-1,-1):
        qs[s]=((up(m,s+1)+down(m,s+1))*qs[s+1]
               -down(m,s+2)*qs[s+2]-V.get((m,s+1),ZERO))/up(m,s)
        Q[m,s]=qs[s]
    for s in range(S+2):
        val=(up(m,s)+down(m,s))*qs.get(s,ZERO)
        if s:val-=up(m,s-1)*qs[s-1]
        val-=down(m,s+1)*qs.get(s+1,ZERO)
        assert (val-V.get((m,s),ZERO)).contains(0)
W=poly(sine(Q),p,True)
assert min(m for m,s in W)==-61
interior_bound=(arb(4)/b**3)*sum((abs(v)*rho**m
                    *(2*(s+max(0,-m))+1)/(J+m)
                    for (m,s),v in W.items()),ZERO)

Cq=sum((j*rho**(j-1)*abs(g[(j,0)])/(2*(j+1)) for j in ms),ZERO)
Lt=rho*rho/(rho*rho-1)
coupling_ratio=Cq/(rho*(J-61))*interior_bound
coupling_shape=Lt*coupling_ratio

# For j>=301 direct finite rows vanish by angular support. The only finite
# contribution comes from the first tail shape coefficient c_{M+2} carried
# into the last finite boundary row by the upper bidiagonal ball operator.
edge_col=arb(54778)/1000 # certified by verify_finite_inverse_r4.py
sumD=sum((abs(x) for x in D.values()),ZERO)
finite_boundary=(edge_col*(M+2)*sumD
                 /(b**3*(J+min(offsets)+2)*(J+2)*rho**(J-1)))
finite_coupling=(edge_col*(M+2)*rho**(-(J-93))
                 /(J-90)*coupling_ratio)
total=(boundary_bound+interior_bound+coupling_shape
       +finite_boundary+finite_coupling)
assert total.upper()<arb(463)/500
print('CERTIFIED FAR SHAPE BOUND ONLY')
print('j>=',J,'boundary upper',boundary_bound.upper())
print('interior upper',interior_bound.upper())
print('trace coupling upper',coupling_shape.upper())
print('finite cross upper',(finite_boundary+finite_coupling).upper())
print('total upper',total.upper())
