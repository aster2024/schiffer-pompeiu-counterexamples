"""Arb Neumann bound for the projected high-source block only."""
from proof_guard import require
import argparse
from flint import arb,ctx
from exact_finite_rank2 import center
from exact_finite_r4 import ZERO


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('center')
    ap.add_argument('--rho-num',type=int,default=102)
    ap.add_argument('--rho-den',type=int,default=100)
    ap.add_argument('--bits',type=int,default=160)
    a=ap.parse_args()
    ctx.prec=a.bits;ctx.threads=1
    d,g,c,B,keys,ms,js=center(a.center)
    rho=arb(a.rho_num)/a.rho_den
    m,M,S=d['m'],d['M'],d['S']
    P=sum((j*rho**(j-1)*abs(v) for j,v in c.items()),ZERO)
    dr=m+2*(S+1)
    ka=1/(arb(dr)*(dr+2))
    n=M+2*m
    kh=1/(arb(2)*(n+1)*(n+2))
    dn=n+2
    kn=1/(arb(dn)*(dn+2))
    require((ka.is_finite() and kh.is_finite() and kn.is_finite()), 'verify_source_tail_rank2.py:26: proof gate failed')
    kmax=ka
    if kh.upper()>kmax.upper():kmax=kh
    if kn.upper()>kmax.upper():kmax=kn
    delta=P*P*kmax
    ceiling=(arb(142)/1000 if (m,M,S)==(6,114,36) else arb(1))
    if not delta.is_finite() or not delta.upper()<ceiling:
        raise ValueError('source-tail perturbation exceeds dimension ceiling')
    inv=1/(1-delta)
    require((inv.is_finite()), 'verify_source_tail_rank2.py:34: proof gate failed')
    print('SOURCE_TAIL_NEUMANN_ONLY m',m,'M',M,'S',S,
          'rho',rho,'P upper',P.upper(),
          'delta upper',delta.upper(),'inverse upper',inv.upper())


if __name__=='__main__':
    main()
