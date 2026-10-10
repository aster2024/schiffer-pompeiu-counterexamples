"""Arb norm bound for the triangular approximate inverse of Section 6.3.

Its finite, Schur and cross-block ceilings are checked by the companion programs."""
from proof_guard import require
from flint import arb,ctx
from exact_finite_rank2 import center,KD,q_of,to_char
from exact_finite_r4 import ZERO


def top(x):
    if not x.is_finite():
        raise ValueError('nonfinite Arb preconditioner bound')
    y=x.upper()
    if not y.is_finite():
        raise ValueError('nonfinite Arb preconditioner upper endpoint')
    return y


def main():
    ctx.prec=192;ctx.threads=1
    d,g,c,B,keys,ms,js=center('center_r14_M114_S36_arbrefined.json')
    require(((d['m'],d['M'],d['S'])==(6,114,36)), 'bound_preconditioner_norm_r14.py:23: proof gate failed')
    rho=arb(102)/100;m=d['m'];b=c[1]
    P=sum((j*rho**(j-1)*abs(v) for j,v in c.items()),ZERO)
    V=KD(g,m)
    V0=sum((rho**n*abs(v) for (n,s),v in V.items()),ZERO)
    H=2*P*V0
    qc=to_char(q_of(g,m),m)
    Cq=sum(((n//m)*rho**(n-m)*abs(v) for n,v in qc.items()),ZERO)
    Q=Cq/(b*rho**m*(m+1))
    require((top(Q)<arb(1)), 'bound_preconditioner_norm_r14.py:32: proof gate failed')
    # Ceilings checked by verify_finite_block_rank2,
    # aggregate_schur_rank2, bound_finite_cross_source_rank2,
    # and aggregate_finite_cross_shape_rank2, respectively.
    Af=arb(82859);Cinv=arb(9);Xs=arb(1)/5;Xc=arb(10)
    bound=Af+(1+Xs)+((1+Xs)*H+Xc+1)*Cinv
    require((top(bound)<arb(100000)), 'bound_preconditioner_norm_r14.py:38: proof gate failed')
    print('PRECONDITIONER_NORM_BOUND_ONLY','P',top(P),'V0',top(V0),
          'H upper',top(H),'Q upper',top(Q),
          'A upper',top(bound))


if __name__=='__main__':
    main()
