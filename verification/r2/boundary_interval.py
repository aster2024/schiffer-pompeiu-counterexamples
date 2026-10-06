#!/usr/bin/env python3
"""All boundary modes of the frozen finite Helmholtz field, rigorously."""
import argparse
import json
from math import comb
from pathlib import Path
from flint import arb,acb,ctx
ctx.threads = 1
from spectral_interval import need,upper,umax,sha,dump,shape_digest
from boundary_newton import exact_decimal,centre_digest
from analytic_cap import field_bounds
finite_upper=upper

def boundary_tail_chunk(state, grid, start, stop, lmax):
    """Interval Fourier coefficients of the explicit finite PDE field."""
    m = state["m"]
    size = state["modes"]+1
    need(grid >= 512 and grid % 2 == 0 and 0 <= start < stop <= grid and
         state["modes"] < lmax < grid // 2,
         "invalid boundary tail grid or chunk")
    p = [exact_decimal(x) for x in state["p"]]
    a = [exact_decimal(x) for x in state["a"]]
    scales = [exact_decimal(x) for x in state["scales"]]
    need(len(p) == len(a) == len(scales) == size and
         all(x > 0 for x in scales), "boundary tail data mismatch")
    dsum = [arb(0) for _ in range(lmax + 1)]
    nsum = [arb(0) for _ in range(lmax + 1)]
    imaginary = acb(arb(0), arb(1))
    for i in range(start, stop):
        phi = 2 * arb.pi() * i / grid
        cos_phi, sin_phi = phi.cos(), phi.sin()
        root = acb(cos_phi, sin_phi)
        powers = [acb(arb(1))]
        for _ in range(1, size):
            powers.append(powers[-1] * root)
        A, B = acb(arb(0)), acb(arb(0))
        for j in range(size):
            A += p[j] * powers[j]
            B += (j * m + 1) * p[j] * powers[j]
        radius = (A * A.conjugate()).real.sqrt()
        need(radius > 0, "nonpositive boundary radius")
        rs = (B * A.conjugate()).real / radius
        theta_s = (B / A).imag
        dv, nv = arb(0), arb(0)
        for j in range(size):
            order = j * m
            J = radius.bessel_j(order) / scales[j]
            Jminus = (-radius.bessel_j(1) if j == 0
                      else radius.bessel_j(order - 1)) / scales[j]
            Jplus = radius.bessel_j(order + 1) / scales[j]
            Jp = (Jminus - Jplus) / 2
            Q = powers[j] * (A / radius)**order
            dv += a[j] * (Q * J).real
            nv += a[j] * (Q * (Jp * rs + imaginary * order * J * theta_s)).real
        c0, c1 = arb(1), cos_phi
        for ell in range(lmax + 1):
            if ell == 0:
                cosine = c0
            elif ell == 1:
                cosine = c1
            else:
                c0, c1 = c1, 2 * cos_phi * c1 - c0
                cosine = c1
            weight = (arb(1) if ell == 0 else arb(2)) * cosine / grid
            dsum[ell] += weight * dv
            nsum[ell] += weight * nv
    for value in dsum + nsum:
        finite_upper(value)
    return {
        "stage": "EXPLICIT_FIELD_BOUNDARY_CHUNK",
        "algorithm": "bessel_boundary_v1",
        "centre_sha256": centre_digest(state),
        "precision_bits": ctx.prec,
        "grid": grid, "start": start, "stop": stop, "lmax": lmax,
        "dirichlet_partial": [x.str(300) for x in dsum],
        "neumann_partial": [x.str(300) for x in nsum],
    }


def merge(state, paths):
    ctx.prec = 4096
    pieces = [json.loads(Path(p).read_text()) for p in paths]
    pieces.sort(key=lambda x: x['start'])
    need(bool(pieces), 'missing boundary leaves')
    grid, lmax = pieces[0]['grid'], pieces[0]['lmax']
    ds, ns = [arb(0) for _ in range(lmax+1)], [arb(0) for _ in range(lmax+1)]
    nxt = 0
    for piece in pieces:
        need(piece['stage'] == 'EXPLICIT_FIELD_BOUNDARY_CHUNK'
             and piece['centre_sha256'] == centre_digest(state)
             and piece['algorithm'] == 'bessel_boundary_v1'
             and piece['precision_bits'] >= 4096, 'boundary leaf provenance mismatch')
        need(piece['grid'] == grid and piece['lmax'] == lmax
             and piece['start'] == nxt < piece['stop'] <= grid,
             'boundary coverage gap/overlap')
        nxt = piece['stop']
        need(len(piece['dirichlet_partial']) == len(piece['neumann_partial']) == lmax+1,
             'boundary dimensions wrong')
        for j, (d, n) in enumerate(zip(piece['dirichlet_partial'], piece['neumann_partial'])):
            dd, nn = arb(d), arb(n)
            upper(dd)
            upper(nn)
            ds[j] += dd
            ns[j] += nn
    need(nxt == grid and grid > 2*lmax, 'incomplete boundary coverage')
    fb = field_bounds(state)
    M = fb['U']+2*fb['dU']
    # Strip |Im(phi)| <= log(4): M bounds BOTH complexified traces.
    qerr = 8*M/arb(4)**(grid-lmax)/(1-arb(1)/4**grid)
    error = arb(0, upper(qerr))
    ds[0] -= 1
    ds = [x+error for x in ds]
    ns = [x+error for x in ns]
    m = state['m']
    # Sum_{t>=0} t^ell / 2^t for ell=0,1,2,3.
    geometric_moments = (2, 2, 6, 26)
    first = lmax+1
    D, N, far = [], [], []
    for k in range(4):
        tail = 2*M/arb(2)**first*sum(
            (arb(comb(k, ell))*(m*first+1)**(k-ell)*m**ell*geometric_moments[ell]
             for ell in range(k+1)), arb(0))
        far.append(tail)
        D.append(sum((abs(v)*2**j*(m*j+1)**k for j, v in enumerate(ds)), arb(0))+tail)
        N.append(sum((abs(v)*2**j*(m*j+1)**k for j, v in enumerate(ns)), arb(0))+tail)
    p = [(m*j+1)*arb(x) for j, x in enumerate(state['p'])]
    P = sum((abs(x)*2**j for j, x in enumerate(p)), arb(0))
    P1 = sum((abs(x)*(m*j)*2**j for j, x in enumerate(p)), arb(0))
    R = fb['R']
    lift = D[0]+N[0]/2
    d_lift = 2*R*(D[2]+N[1])
    lap_lift = 2*(D[2]+N[1])
    dlap_lift = 2*R*(D[3]+N[2])
    E0 = lap_lift+P*P*lift
    E1 = dlap_lift+R*P*P1*lift+P*P*d_lift
    quantities = dict(strip_trace_majorant_upper=M,
                      quadrature_error_per_coefficient=qerr,
                      lift_upper=lift, dlift_upper=d_lift,
                      lap_lift_upper=lap_lift, dlap_lift_upper=dlap_lift,
                      E0_upper=E0, E1_upper=E1)
    result = {'stage': 'CLAMPED_FIELD_RESIDUAL_AND_FIRST_DERIVATIVE',
              'shape_sha256': shape_digest(state),
              'center_sha256': centre_digest(state), 'grid': grid, 'lmax': lmax,
              'modes': state['modes'], 'full_schiffer_proof': False,
              'D_moments_upper': [upper(x).str(55) for x in D],
              'N_moments_upper': [upper(x).str(55) for x in N],
              'each_far_moment_upper': [upper(x).str(55) for x in far]}
    result.update({k: upper(v).str(55) for k, v in quantities.items()})
    return result


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--stage', choices=('leaf', 'merge'), required=True)
    ap.add_argument('--state', required=True)
    ap.add_argument('--grid', type=int, default=1024)
    ap.add_argument('--lmax', type=int, default=400)
    ap.add_argument('--start', type=int, default=0)
    ap.add_argument('--stop', type=int)
    ap.add_argument('--leaf', action='append', default=[])
    ap.add_argument('--output', required=True)
    args = ap.parse_args()
    ctx.prec = 4096
    state = json.loads(Path(args.state).read_text())
    if args.stage == 'leaf':
        result = boundary_tail_chunk(state, args.grid, args.start, args.stop, args.lmax)
        dump(args.output, result)
        print('BOUNDARY_LEAF', args.start, args.stop, flush=True)
    else:
        result = merge(state, args.leaf)
        dump(args.output, result)
        print('CLAMPED_RESIDUAL_COMPONENT', result['E0_upper'], result['E1_upper'], flush=True)


if __name__ == '__main__':
    main()
