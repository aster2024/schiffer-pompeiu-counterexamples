"""Arb radii, sign, convexity and nonball gates for the R14 center.

The A<100000 and Z<1/2 hypotheses are proved by companion block scripts;
this file checks their nonlinear and geometric consequences only.
"""
from proof_guard import require
import hashlib
import json
from flint import arb,ctx
from exact_finite_rank2 import center,evaluate,norm,KD,q_of,h_of,to_char
from exact_finite_r4 import ZERO
from verify_center_geometry_sign_rank2 import value


def top(x):
    require((x.is_finite()), 'verify_radii_geometry_r14.py:15: proof gate failed')
    y=x.upper();require((y.is_finite()), 'verify_radii_geometry_r14.py:16: proof gate failed')
    return y


def low(x):
    require((x.is_finite()), 'verify_radii_geometry_r14.py:21: proof gate failed')
    y=x.lower();require((y.is_finite()), 'verify_radii_geometry_r14.py:22: proof gate failed')
    return y


def dyad(pair):
    a,e=pair
    x=arb(int(a))*arb(2)**int(e)
    require((x.is_finite()), 'verify_radii_geometry_r14.py:29: proof gate failed')
    return x


def char_lower(v,m):
    return v.get(m,ZERO)-sum(((n//m)*abs(a) for n,a in v.items() if n>m),ZERO)


def main():
    ctx.prec=192;ctx.threads=1
    name='center_r14_M186_S48_arbrefined2.json'
    d,g,c,B,keys,ms,js=center(name)
    require(((d['m'],d['M'],d['S'])==(6,186,48)), 'verify_radii_geometry_r14.py:41: proof gate failed')
    old,oldg,oldc,oldB,*_=center('center_r14_M114_S36_arbrefined.json')
    require((old['m']==6 and B==oldB), 'verify_radii_geometry_r14.py:43: proof gate failed')
    m=6;rho=arb(51)/50;r=arb(1)/arb(10)**20
    with open('r14_M186_nonlinear_bounds.json') as f:major=json.load(f)
    require((major['center_sha256']==hashlib.sha256(open(name,'rb').read()).hexdigest()), 'verify_radii_geometry_r14.py:46: proof gate failed')
    require(((major['rho_num'],major['rho_den'])==(102,100)), 'verify_radii_geometry_r14.py:47: proof gate failed')
    dq=[dyad(v) for v in major['coeff_upper_dyadic']]
    hk=[dyad(v) for v in major['h_degree_upper_dyadic']]
    require((len(dq)==13 and len(hk)==7), 'verify_radii_geometry_r14.py:50: proof gate failed')
    if not top(dq[0])<arb(19):
        raise ValueError('quadratic nonlinear majorant fails')
    F,G=evaluate(g,c,B,m)
    residual=top(norm(F,G,m,rho))
    require((residual<arb(2)/arb(10)**30), 'verify_radii_geometry_r14.py:55: proof gate failed')
    # evaluate() divides the boundary output by the argument's c_1. The
    # preconditioner and Z proof fix that denominator at the old center.
    # Their residual and all nonlinear coefficients differ by at most this
    # common conservative factor in the sum output norm.
    ratio=c[1]/oldc[1]
    require((ratio.is_finite() and ratio.lower()>0 and ratio.upper()<arb(2)), 'verify_radii_geometry_r14.py:61: proof gate failed')
    factor=arb(1)
    if ratio.upper()>factor:factor=ratio.upper()
    A=arb(100000);Z=arb(1)/2
    Y=A*factor*residual
    mapgap=Y+(Z-1)*r+sum((A*factor*dq[q-2]*r**q for q in range(2,15)),ZERO)
    lipgap=Z-1+sum((q*A*factor*dq[q-2]*r**(q-1) for q in range(2,15)),ZERO)
    require((top(mapgap)<0 and top(lipgap)<0), 'verify_radii_geometry_r14.py:68: proof gate failed')

    b=c[1]
    L=b-sum((j*abs(v) for j,v in c.items() if j>1),ZERO)
    U=sum((j*(j-1)*abs(v) for j,v in c.items() if j>1),ZERO)
    outer=arb(1001)/1000
    Lout=b-sum((j*outer**(j-1)*abs(v) for j,v in c.items() if j>1),ZERO)
    require((low(Lout-r)>0 and low(L-U-51*r)>0), 'verify_radii_geometry_r14.py:75: proof gate failed')
    require((low(abs(c[13])-r/((13+2*m)*rho**12))>0), 'verify_radii_geometry_r14.py:76: proof gate failed')

    qc=to_char(q_of(g,m),m)
    hc=to_char(h_of(c,m,B),m)
    qlow=char_lower(qc,m)
    hlow=char_lower(hc,m)
    dq_sup=r*rho**(-m)/(2*m)
    dh_sup=sum((hk[k]*r**k for k in range(1,m+1)),ZERO)
    require((low(qlow-dq_sup)>0 and low(hlow-dh_sup)>0), 'verify_radii_geometry_r14.py:84: proof gate failed')

    V=KD(g,m)
    minus=value(V,ms,m,d['S'],arb(1)/10)
    plus=value(V,ms,m,d['S'],arb(37)/100)
    kappa=arb(1)/80
    require((top(minus+kappa*r)<0 and low(plus-kappa*r)>0), 'verify_radii_geometry_r14.py:90: proof gate failed')
    print('RADII_AND_GEOMETRY_PASS_ONLY',
          'radius',r,'Y upper',top(Y),
          'fixed-boundary-scale factor upper',top(factor),
          'map gap upper',top(mapgap),'Lipschitz gap upper',top(lipgap),
          'convex numerator margin lower',low(L-U-51*r),
          'q lower',low(qlow-dq_sup),'h lower',low(hlow-dh_sup),
          'sign values',minus,plus)


if __name__=='__main__':
    main()
