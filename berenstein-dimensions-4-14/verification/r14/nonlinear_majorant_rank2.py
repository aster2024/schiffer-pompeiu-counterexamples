"""Arb multilinear remainder coefficients for Section 6.5.

Optional radius calculations use the supplied inverse and derivative-defect bounds."""
import argparse
import hashlib
import json
from math import comb
from flint import arb,ctx
from exact_finite_rank2 import center,evaluate,norm,KD,h_of,to_char
from exact_finite_r4 import ZERO


def top(v):
    if not v.is_finite():
        raise ValueError('nonfinite Arb majorant')
    y=v.upper()
    if not y.is_finite():
        raise ValueError('nonfinite Arb upper endpoint')
    return y


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('center')
    ap.add_argument('--rho-num',type=int,default=102)
    ap.add_argument('--rho-den',type=int,default=100)
    ap.add_argument('--bits',type=int,default=192)
    ap.add_argument('--hypothetical-A',type=int,default=1000000)
    ap.add_argument('--hypothetical-Z-num',type=int,default=1)
    ap.add_argument('--hypothetical-Z-den',type=int,default=2)
    ap.add_argument('--bound-output',default='')
    a=ap.parse_args()
    ctx.prec=a.bits;ctx.threads=1
    d,g,c,B,keys,ms,js=center(a.center)
    m=d['m'];b=c[1];rho=arb(a.rho_num)/a.rho_den
    P=sum((j*rho**(j-1)*abs(v) for j,v in c.items()),ZERO)
    C=sum((rho**(j-1)*abs(v) for j,v in c.items()),ZERO)
    V=KD(g,m)
    V0=sum((rho**n*abs(v) for (n,s),v in V.items()),ZERO)
    kappa_h=1/(arb(2)*(m+1)*(m+2))
    kappa_n=1/(arb(m+2)*(m+4))
    kappa=kappa_h if kappa_h.upper()>kappa_n.upper() else kappa_n
    theta0=arb(1)/(1+2*m)
    # The multiplier p=psi' acts in rho-weighted G and U norms.  Relative to
    # the shape weight (j+2m)*rho**(j-1), its norm is bounded by
    # sup_j j/(j+2m) = 1.  The old rho**(-(j-1)) factor only bounded an
    # unweighted Wiener norm and was invalid for the nonlinear remainder.
    theta1=arb(1);arg='weighted supremum'
    factor=arb(1)/2 if m==2 else arb(1)
    h=to_char(h_of(c,m,B),m)
    H0=[sum((rho**(n-m)*abs(v) for n,v in h.items()),ZERO)]
    H1=[sum(((n//m)*rho**(n-m)*abs(v) for n,v in h.items()),ZERO)]
    for k in range(1,m+1):
        base=factor*comb(m,k)/B
        H0.append(base*C**(m-k)*theta0**k)
        p1=(arb(m-k)/m*P*C**(m-k-1)*theta0**k) if k<m else ZERO
        p2=arb(k)/m*C**(m-k)*theta1*theta0**(k-1)
        H1.append(base*(p1+p2))
    HH=[]
    for degree in range(2*m+1):
        s=arb(0)
        for k in range(m+1):
            l=degree-k
            if not 0<=l<=m:continue
            x=top(H1[k]*H0[l]);y=top(H0[k]*H1[l])
            s+=x if x<y else y
        HH.append(s)
    aq0=rho**(-m)/(2*(m+1))
    aq1=rho**(-m)/(2*m)
    coeff=[]
    for degree in range(2,2*m+3):
        boundary=(HH[degree]*P*P if degree<len(HH) else ZERO)
        if degree-1<len(HH):boundary+=HH[degree-1]*2*P*theta1
        if degree-2<len(HH):boundary+=HH[degree-2]*theta1*theta1
        boundary/=b
        if degree==2:boundary+=aq0*aq1/b
        interior=ZERO
        if degree==2:interior=2*P*theta1*kappa+V0*theta1*theta1
        if degree==3:interior=kappa*theta1*theta1
        coeff.append(top(interior+boundary))
    F,G=evaluate(g,c,B,m)
    residual=top(norm(F,G,m,rho))
    deriv_sum=top(sum(((q+2)*coeff[q] for q in range(len(coeff))),ZERO))
    if a.bound_output:
        def pair(x):
            man,exp=top(x).man_exp()
            return [int(man),int(exp)]
        payload=dict(center_sha256=hashlib.sha256(open(a.center,'rb').read()).hexdigest(),
                     m=m,rho_num=a.rho_num,rho_den=a.rho_den,
                     bits=a.bits,derivative_sum_upper_dyadic=pair(deriv_sum),
                     coeff_upper_dyadic=[pair(x) for x in coeff],
                     h_degree_upper_dyadic=[pair(x) for x in H1])
        with open(a.bound_output,'w') as f:json.dump(payload,f,separators=(',',':'))
        print('saved nonlinear bounds',a.bound_output)
    A=arb(a.hypothetical_A)
    Z=arb(a.hypothetical_Z_num)/a.hypothetical_Z_den
    print('CONDITIONAL_MAJORANT_ONLY m',m,'rho',rho,
          'theta0',top(theta0),'theta1',top(theta1),'arg',arg,
          'P',top(P),'C',top(C),'V0',top(V0),'kappa',top(kappa),
          'residual',residual)
    print('sum q*dq upper',deriv_sum)
    print('raw nonlinear coefficients',[(k+2,str(v)) for k,v in enumerate(coeff)])
    print('hypothetical A',A,'Z',Z)
    for k in range(4,22):
        r=arb(1)/arb(10)**k
        mapgap=A*residual+(Z-1)*r+sum((A*coeff[q-2]*r**q
                    for q in range(2,2*m+3)),ZERO)
        lipgap=Z-1+sum((q*A*coeff[q-2]*r**(q-1)
                    for q in range(2,2*m+3)),ZERO)
        if top(mapgap)<0 and top(lipgap)<0:
            print('conditional radius',r,'mapgap upper',top(mapgap),
                  'lipgap upper',top(lipgap))


if __name__=='__main__':
    main()
