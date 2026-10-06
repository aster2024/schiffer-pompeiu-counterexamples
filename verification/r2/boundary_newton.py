#!/usr/bin/env python3
"""Generic finite boundary Newton proposals; not an existence certificate."""
import argparse
import hashlib
import json
import time
from pathlib import Path
from flint import arb,acb,arb_mat,ctx
ctx.threads = 1
from spectral_interval import need,upper,umax,exact_decimal
finite_upper=upper

def centre_digest(state):
    return hashlib.sha256(json.dumps({k:state[k] for k in ("m","modes","p","a","scales")},sort_keys=True).encode()).hexdigest()

def interval_traces_chunk(state, grid, start, stop):
    """Arb enclosure of a finite trapezoidal sum and its analytic Jacobian.

    The continuum quadrature error is not added at this stage.  The result
    therefore has the explicit status DISCRETE_INTERVAL_ONLY.
    """
    m = state["m"]
    need(isinstance(grid, int) and grid >= 32 and grid % 2 == 0,
         "invalid grid")
    need(0 <= start < stop <= grid, "invalid chunk")
    size = state["modes"] + 1
    p = [exact_decimal(x) for x in state["p"]]
    a = [exact_decimal(x) for x in state["a"]]
    scales = [exact_decimal(x) for x in state["scales"]]
    need(len(p) == len(a) == len(scales) == size,
         "coefficient count mismatch")
    need(all(x > 0 for x in scales), "nonpositive basis scale")
    dsum = [arb(0) for _ in range(size)]
    nsum = [arb(0) for _ in range(size)]
    jd = [[arb(0) for _ in range(2 * size)] for _ in range(size)]
    jn = [[arb(0) for _ in range(2 * size)] for _ in range(size)]
    imaginary = acb(arb(0), arb(1))
    for i in range(start, stop):
        phi = arb(2) * arb.pi() * i / grid
        cos_phi, sin_phi = phi.cos(), phi.sin()
        root = acb(cos_phi, sin_phi)
        exps = [acb(arb(1))]
        for _ in range(1, size):
            exps.append(exps[-1] * root)
        A, B = acb(arb(0)), acb(arb(0))
        for k in range(size):
            A += p[k] * exps[k]
            B += (k * m + 1) * p[k] * exps[k]
        radius = (A * A.conjugate()).real.sqrt()
        need(radius > 0, "boundary radius not positive")
        rs = (B * A.conjugate()).real / radius
        theta_s = (B / A).imag
        dr = [(exps[k] * A.conjugate()).real / radius
              for k in range(size)]
        delta = [(exps[k] / A).imag for k in range(size)]
        drs = [(((k * m + 1) * exps[k] * A.conjugate()
                 + B * exps[k].conjugate()).real / radius
                - rs * dr[k] / radius) for k in range(size)]
        dtheta = [((k * m + 1) * exps[k] / A
                   - B * exps[k] / A**2).imag for k in range(size)]
        field_d, field_n = arb(0), arb(0)
        shape_d = [arb(0) for _ in range(size)]
        shape_n = [arb(0) for _ in range(size)]
        basis_d, basis_n = [], []
        for j in range(size):
            order = j * m
            J = radius.bessel_j(order) / scales[j]
            if order == 0:
                Jm = -radius.bessel_j(1) / scales[j]
            else:
                Jm = radius.bessel_j(order - 1) / scales[j]
            Jplus = radius.bessel_j(order + 1) / scales[j]
            Jp = (Jm - Jplus) / 2
            Jpp = ((arb(order)**2 / radius**2 - 1) * J
                   - Jp / radius)
            Q = exps[j] * (A / radius)**order
            core_n = Jp * rs + imaginary * order * J * theta_s
            h = (Q * J).real
            v = (Q * core_n).real
            basis_d.append(h)
            basis_n.append(v)
            field_d += a[j] * h
            field_n += a[j] * v
            for k in range(size):
                shape_d[k] += a[j] * (Q * (
                    Jp * dr[k] + imaginary * order * J * delta[k]
                )).real
                shape_n[k] += a[j] * (Q * (
                    Jpp * dr[k] * rs + Jp * drs[k]
                    + imaginary * order * (Jp * dr[k] * theta_s
                                           + J * dtheta[k])
                    + imaginary * order * delta[k] * core_n
                )).real
        for ell in range(size):
            weight = (arb(1) if ell == 0 else arb(2)) * (ell * phi).cos() / grid
            dsum[ell] += weight * field_d
            nsum[ell] += weight * field_n
            for k in range(size):
                jd[ell][k] += weight * shape_d[k]
                jn[ell][k] += weight * shape_n[k]
                jd[ell][size + k] += weight * basis_d[k]
                jn[ell][size + k] += weight * basis_n[k]
    all_values = dsum + nsum + [x for row in jd + jn for x in row]
    for value in all_values:
        finite_upper(value)
    return {
        "stage": "DISCRETE_INTERVAL_ONLY",
        "algorithm": "fourier_bessel_trace_v1",
        "centre_sha256": centre_digest(state),
        "precision_bits": ctx.prec,
        "grid": grid, "start": start, "stop": stop,
        "residual_partial": [x.str(300) for x in dsum + nsum],
        "jacobian_partial": [[x.str(300) for x in row]
                             for row in jd + jn],
    }



def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--state',required=True)
    ap.add_argument('--modes',type=int,required=True)
    ap.add_argument('--grid',type=int,default=256)
    ap.add_argument('--bits',type=int,default=4096)
    ap.add_argument('--digits',type=int,default=150)
    ap.add_argument('--output',required=True)
    args=ap.parse_args()
    ctx.prec=args.bits
    state=json.loads(Path(args.state).read_text())
    old=len(state['p'])
    need(args.modes+1>=old,'cannot drop modes')
    for j in range(old,args.modes+1):
        scale=abs(arb(state['p'][0]).bessel_j(j*state['m']))
        need(scale>0,'Bessel scale unresolved')
        state['p'].append('0')
        state['a'].append('0')
        state['scales'].append(scale.mid().str(args.digits,radius=False))
    state['modes']=args.modes
    n=args.modes+1
    start=time.monotonic()
    data=interval_traces_chunk(state,args.grid,0,args.grid)
    residual=[arb(s) for s in data['residual_partial']]
    residual[0]-=1
    jac=arb_mat([[arb(s) for s in row] for row in data['jacobian_partial']])
    correction=jac.solve(arb_mat([[r] for r in residual]))
    values=[arb(s)-correction[i,0] for i,s in enumerate(state['p']+state['a'])]
    state['p']=[v.mid().str(args.digits,radius=False) for v in values[:n]]
    state['a']=[v.mid().str(args.digits,radius=False) for v in values[n:]]
    state['status']='FINITE_NUMERICAL_ONLY'
    state['grid']=args.grid
    state['bits']=args.bits
    state['steps']=state.get('steps',0)+1
    state['last_residual_upper']=umax(abs(v) for v in residual).str(25)
    state['last_correction_upper']=umax(abs(correction[i,0]) for i in range(2*n)).str(25)
    state['last_wall_seconds']=time.monotonic()-start
    Path(args.output).write_text(json.dumps(state,indent=2)+'\n')
    print('NUMERICAL_NEWTON',state['modes'],state['last_residual_upper'],state['last_correction_upper'],state['last_wall_seconds'],flush=True)

if __name__=='__main__':
    main()
