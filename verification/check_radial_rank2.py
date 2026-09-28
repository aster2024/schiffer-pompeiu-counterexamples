"""Exact su(3) matrix-coordinate and root-coordinate checks for B2/G2.

These checks corroborate, but do not replace, the Weyl-Jacobian proof.
"""

import argparse
import math
import random
import sympy as sp


def need(ok, why):
    if not ok:
        raise RuntimeError(why)


def su3():
    x = sp.symbols("x1:9", real=True)
    a = sp.sqrt(2)
    b = sp.sqrt(6)
    I = sp.I
    X = sp.Matrix([
        [x[2]/a+x[7]/b, (x[0]-I*x[1])/a, (x[3]-I*x[4])/a],
        [(x[0]+I*x[1])/a, -x[2]/a+x[7]/b, (x[5]-I*x[6])/a],
        [(x[3]+I*x[4])/a, (x[5]+I*x[6])/a, -2*x[7]/b],
    ])
    p2 = sp.expand(sp.trace(X*X))
    p3 = sp.expand(sp.trace(X*X*X))
    grad2 = [sp.diff(p2, xi) for xi in x]
    grad3 = [sp.diff(p3, xi) for xi in x]
    checks = {
        "p2_norm": p2 - sum(xi**2 for xi in x),
        "lap_p2": sum(sp.diff(p2, xi, 2) for xi in x) - 16,
        "lap_p3": sum(sp.diff(p3, xi, 2) for xi in x),
        "grad_p2_sq": sum(y*y for y in grad2) - 4*p2,
        "grad_cross": sum(y*z for y,z in zip(grad2,grad3)) - 6*p3,
        "grad_p3_sq": sum(y*y for y in grad3) - sp.Rational(3,2)*p2**2,
    }
    for name, expr in checks.items():
        need(sp.simplify(sp.expand(expr)) == 0, name)
    # On the orthonormal Cartan plane (x3,x8), pi is the Vandermonde.
    t, s = x[2], x[7]
    diagonal = [X[i,i] for i in range(3)]
    pi = sp.expand((diagonal[0]-diagonal[1])
                   *(diagonal[0]-diagonal[2])
                   *(diagonal[1]-diagonal[2]))
    need(sp.simplify(sp.diff(pi,t,2)+sp.diff(pi,s,2)) == 0, 'pi not harmonic')
    for p, target in [(p2, 12), (p3, 0)]:
        p_cartan = sp.simplify(p.subs({xi:0 for xi in x if xi not in (t,s)}))
        drift = 2*(sp.diff(pi,t)*sp.diff(p_cartan,t)
                   +sp.diff(pi,s)*sp.diff(p_cartan,s))
        need(sp.simplify(sp.expand(drift-target*pi)) == 0, 'drift identity failed')
    print("SU3_MATRIX_EXACT_PASS", len(checks), "coordinate identities; pi harmonic; drift on p2,p3")


def roots(m, samples=100):
    # Unit normals to the m reflection walls. Root lengths only change the
    # product by a constant and hence cancel from the radial operator.
    normals = [(-math.sin(k*math.pi/m), math.cos(k*math.pi/m))
               for k in range(m)]
    x, y = sp.symbols("x y", real=True)
    pi = sp.im(sp.expand((x+sp.I*y)**m)).expand(complex=True)
    need(sp.expand(sp.diff(pi,x,2)+sp.diff(pi,y,2)) == 0, 'pi not harmonic')
    invariant = (x*x+y*y)**2 + sp.re(sp.expand((x+sp.I*y)**m)).expand(complex=True)
    f = sp.lambdify((x,y), invariant, "math")
    fx = sp.lambdify((x,y), sp.diff(invariant,x), "math")
    fy = sp.lambdify((x,y), sp.diff(invariant,y), "math")
    lap = sp.lambdify((x,y), sp.diff(invariant,x,2)+sp.diff(invariant,y,2), "math")
    p = sp.lambdify((x,y), pi, "math")
    dpx = sp.lambdify((x,y), sp.diff(pi,x), "math")
    dpy = sp.lambdify((x,y), sp.diff(pi,y), "math")
    lap_pi_f = sp.lambdify((x,y), sp.diff(pi*invariant,x,2)+sp.diff(pi*invariant,y,2), "math")
    rng = random.Random(17+m)
    max_rel = 0.0
    accepted = 0
    while accepted < samples:
        xx, yy = rng.uniform(-1.7,1.7), rng.uniform(-1.7,1.7)
        if min(abs(nx*xx+ny*yy) for nx,ny in normals) < 0.08:
            continue
        drift = sum(2*(nx*fx(xx,yy)+ny*fy(xx,yy))/(nx*xx+ny*yy)
                    for nx,ny in normals)
        lhs = lap(xx,yy)+drift
        rhs = lap_pi_f(xx,yy)/p(xx,yy)
        log_rhs = lap(xx,yy)+2*(dpx(xx,yy)*fx(xx,yy)+dpy(xx,yy)*fy(xx,yy))/p(xx,yy)
        rel = max(abs(lhs-rhs), abs(lhs-log_rhs))/(1+abs(rhs))
        need(math.isfinite(rel), 'nonfinite sample')
        max_rel = max(max_rel, rel)
        accepted += 1
    need(max_rel < 1e-10, 'radial-part identity failed')
    print(f"I2({m})_ROOT_NUMERIC_PASS samples={samples} max_relative={max_rel:.3e}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", type=int, default=100)
    args = parser.parse_args()
    su3()
    roots(4,args.samples)
    roots(6,args.samples)
