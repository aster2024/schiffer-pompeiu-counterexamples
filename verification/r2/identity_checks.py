#!/usr/bin/env python3
"""Exact polynomial sanity checks, supplementary to the universal proof."""
import json
from pathlib import Path
import sympy as s


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def main():
    x, y = s.symbols('x y', real=True)
    z = x+s.I*y
    lap = lambda f: s.diff(f, x, 2)+s.diff(f, y, 2)
    rr = x*x+y*y
    U = 1+(1-rr)**2*(2+x*x-y*y)
    p = 2+z*z/7
    q = s.expand(p*s.conjugate(p))
    E = lap(U)+q*U
    count = 0
    for degree in (1, 3, 5):
        v = s.expand(z**degree)
        vx, vy = s.re(v), s.im(v)
        vp = degree*z**(degree-1)
        eta_prime = s.expand((2*z/7)*v+p*vp)
        dq = s.expand(s.conjugate(p)*eta_prime+p*s.conjugate(eta_prime))
        T = vx*s.diff(U, x)+vy*s.diff(U, y)
        residual = lap(T)+q*T+dq*U-vx*s.diff(E, x)-vy*s.diff(E, y)-2*s.re(vp)*E
        need(s.expand(residual) == 0, 'Cartesian material derivative identity failed')
        normal_T = x*s.diff(T, x)+y*s.diff(T, y)
        boundary_residual = s.expand(normal_T-(E-q)*(x*vx+y*vy))
        need(s.rem(boundary_residual, y*y+x*x-1, y) == 0,
             'clamped Hessian normal-trace identity failed')
        count += 1
    for n in (0, 1, 2, 5):
        h = s.re(s.expand(z**n))
        hd = h*(1-s.Rational(n, 2)*(rr-1))
        hn = h*(rr-1)/2
        need(s.expand(lap(hd)+2*n*(n+1)*h) == 0, 'Dirichlet lift Laplacian failed')
        need(s.expand(lap(hn)-2*(n+1)*h) == 0, 'Neumann lift Laplacian failed')
        nd = s.expand(x*s.diff(hd, x)+y*s.diff(hd, y))
        nn = s.expand(x*s.diff(hn, x)+y*s.diff(hn, y)-h)
        need(s.rem(nd, y*y+x*x-1, y) == 0 and s.rem(nn, y*y+x*x-1, y) == 0,
             'lift boundary derivative failed')
        count += 1
    result = {'status': 'EXACT_POLYNOMIAL_SANITY_CHECKS_PASS', 'cases': count,
              'scope': 'implementation checks only; universal arguments are in PROOF_PACKAGE.md'}
    Path('identity_checks.json').write_text(json.dumps(result, indent=2)+'\n')
    print(result['status'], count)


if __name__ == '__main__':
    main()
