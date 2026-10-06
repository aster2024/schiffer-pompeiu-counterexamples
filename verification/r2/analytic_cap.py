#!/usr/bin/env python3
"""Fixed-disc radial-sup / angular-Wiener coupled inverse bounds.

The analytic lemmas and their domains must be proved in PROOF_PACKAGE.md.
This module checks their scalar inequalities, not their derivation.
"""
import argparse
import json
from pathlib import Path
from flint import arb, ctx
ctx.threads = 1
from spectral_interval import need, upper, umax, sha, dump, shape_digest, exact_decimal
from boundary_newton import centre_digest


def field_bounds(state):
    m = state['m']
    c = [exact_decimal(x) for x in state['p']]
    a = [exact_decimal(x) for x in state['a']]
    scales = [exact_decimal(x) for x in state['scales']]
    need(len(c) == len(a) == len(scales) == state['modes']+1 and all(s > 0 for s in scales),
         'invalid finite field')
    rho = c[0]
    R = (arb(2).log()/m).exp()
    R4 = (arb(4).log()/m).exp()
    S = sum((abs(cj)*4**j for j, cj in enumerate(c) if j), arb(0))
    need(rho > 2*S, 'complex angular strip reaches a zero of A')
    displacement = (2*rho*S+S*S)/rho
    phase = ((rho+S)/(rho-S)).log()/2
    P4 = sum((abs(cj)*(m*j+1)*4**j for j, cj in enumerate(c)), arb(0))
    U, dU = arb(0), arb(0)
    for j, (aj, scale) in enumerate(zip(a, scales)):
        n = j*m
        alpha = (arb(n)/(rho+displacement)).acosh() if n > rho+displacement else arb(0)
        exponent = -n*alpha+(rho+displacement)*alpha.sinh()+displacement*alpha.cosh()
        majorant = abs(aj)/scale*4**j*(n*phase+exponent).exp()
        U += majorant
        dU += majorant*alpha.exp()
    need(2*U >= 1, 'Dirichlet strip constant not absorbed')
    return {
        'U': 3*U,
        'dU': 3*R*R4*R4*P4*phase.exp()*dU,
        'strip_displacement': displacement,
        'R': R,
    }


def high_sector_gap(high, m, rho):
    gap = arb((high*m)**2)-rho*rho
    need(gap > 0, 'nonpositive high angular gap d_h')
    return gap


def bounds(state, reference, spectral, boundary, radius, require_radii=True):
    ctx.prec = 256
    m = state['m']
    need(m == reference['m'] == spectral['m'], 'symmetry mismatch')
    need(spectral['stage'] == 'INVARIANT_SECTOR_DIRICHLET_INVERSE', 'wrong spectral scope')
    need(spectral['reference_shape_sha256'] == shape_digest(reference),
         'Dirichlet inverse reference shape mismatch')
    if require_radii:
        need(boundary['stage'] == 'CLAMPED_FIELD_RESIDUAL_AND_FIRST_DERIVATIVE'
             and boundary['shape_sha256'] == shape_digest(state)
             and boundary['center_sha256'] == centre_digest(state),
             'boundary residual shape mismatch')
    c = [exact_decimal(x) for x in state['p']]
    p = [(j*m+1)*x for j, x in enumerate(c)]
    p_ref = [(j*m+1)*exact_decimal(x) for j, x in enumerate(reference['p'])]
    rho = c[0]
    P = sum((abs(x)*2**j for j, x in enumerate(p)), arb(0))
    P1 = sum((abs(x)*(j*m)*2**j for j, x in enumerate(p)), arb(0))
    Pun = sum((abs(x) for x in p), arb(0))
    Pref = sum((abs(x) for x in p_ref), arb(0))
    change = sum((abs((p[j] if j < len(p) else arb(0))-
                      (p_ref[j] if j < len(p_ref) else arb(0)))
                  for j in range(max(len(p), len(p_ref)))), arb(0))
    delta_q = (Pun+Pref)*change
    B_ref = arb(spectral['full_Dirichlet_inverse_L2_norm_upper'])
    upper(B_ref)
    need(B_ref > 0 and B_ref*delta_q < 1, 'Dirichlet inverse transfer failed')
    BL2 = B_ref/(1-B_ref*delta_q)
    high = 1
    while high*m <= rho:
        high += 1
    gap = high_sector_gap(high, m, rho)
    potential = P*P-rho*rho
    beta = potential/gap
    need(beta < 1, 'analytic angular high-sector inverse failed')
    # Sum of all complex Fourier weights for |j|<high.
    low_weight = arb(1+2*sum(2**j for j in range(1, high)))
    low_sup = low_weight*(1+P*P*BL2)/(2*arb(2).sqrt())
    BY = (low_sup+1/gap)/(1-beta)
    Btrace = 1+P*P*BY
    f = field_bounds(state)
    R = f['R']
    E0 = arb(boundary['E0_upper'])
    E1 = arb(boundary['E1_upper'])
    lift = arb(boundary['lift_upper'])
    lap_lift = arb(boundary['lap_lift_upper'])
    dlap_lift = arb(boundary['dlap_lift_upper'])
    for value in (E0, E1, lift, lap_lift, dlap_lift):
        upper(value)
    need(all(x >= 0 for x in (E0, E1, lift, lap_lift, dlap_lift)),
         'negative boundary bound')
    D0 = P*P*f['U']+lap_lift
    D1 = R*P*P1*f['U']+P*P*f['dU']+dlap_lift
    # Boundary d=q-E and its reciprocal in the (1+|n|) Wiener norm.
    denom = rho*rho-potential-E0
    need(denom > 0, 'q-E not separated from zero in Wiener algebra')
    I0 = 1/denom
    I1 = I0+I0*I0*(2*P*P1+2*R*E1)
    Bp = (P+P1)*I1*Btrace
    Bg = Btrace+2*(R*D1+D0)*I1*Btrace
    A = Bg+Bp
    # a=delta_psi/(w*p) has W1 norm <= J1 ||delta_p||_W.
    hol_gap = rho-(P-rho)
    need(hol_gap > 0, 'p not separated from zero on analytic disk')
    J0 = 1/hol_gap
    J1 = J0+P1*J0*J0
    Z = A*2*(R*E1+E0)*J1
    Y = A*E0
    kappa = arb(1)/4
    Uc = f['U']+lift
    C2 = A*(2*P*kappa+Uc)
    C3 = A*kappa
    r = exact_decimal(radius)
    need(r > 0, 'nonpositive existence radius')
    poly = Y+(Z-1)*r+C2*r*r+C3*r*r*r
    contraction = Z+2*C2*r+3*C3*r*r
    # The X norm is ||delta_g||_Y + ||delta_p||_W, with p=psi'.
    center_first = sum((abs(x) for x in p[1:]), arb(0))
    center_second = sum((abs(x)*(j*m) for j, x in enumerate(p) if j), arb(0))
    first_margin = rho-center_first-r
    second_upper = center_second+arb(m+1)*r/2
    curvature = 1-second_upper/first_margin
    convex_margin = rho-center_first-center_second-arb(m+1)*r/2
    nondisc = abs(c[1])-r/(2*(m+1))
    analytic_univalence = rho-(P-rho)-r
    need(first_margin > 0 and curvature > 0 and convex_margin > 0
         and nondisc > 0 and analytic_univalence > 0, 'whole-ball geometry failed')
    if require_radii:
        need(poly < 0 and contraction < 1, 'full-space radii polynomial failed')
    quantities = dict(B_L2=BL2, spectral_transfer_delta_q=delta_q,
                      high_angular_gap_lower=gap.lower(), high_angular_start=arb(high),
                      high_angular_beta=beta, B_Y=BY, B_trace=Btrace,
                      U_field=f['U'], dU_field=f['dU'], D0=D0, D1=D1,
                      inverse_W1_boundary=I1, inverse_W1_holomorphic=J1,
                      B_shape=Bp, B_field=Bg, A=A, Y=Y, Z=Z, C2=C2, C3=C3,
                      radius=r, radii_polynomial=poly, contraction=contraction,
                      four_C2_Y=4*C2*Y, analytic_radius_lower=R.lower(),
                      unit_derivative_lower=first_margin.lower(),
                      curvature_lower=curvature.lower(),
                      convexity_margin_lower=convex_margin.lower(),
                      non_disc_c1_lower=nondisc.lower(),
                      analytic_univalence_margin_lower=analytic_univalence.lower())
    for x in quantities.values():
        upper(x)
    return {'stage': 'COUPLED_RADII_GATES' if require_radii else 'CONSTANTS_DIAGNOSTIC',
            'norm': '||g||_{sum 2^|j| C_r}+||p||_{sum 2^j}, p=psi prime',
            'analytic_lemma_version': 'fixed_disc_wiener_clamped_v1',
            'bounds': {k: v.str(50) for k, v in quantities.items()}}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--state', required=True)
    ap.add_argument('--reference', required=True)
    ap.add_argument('--spectral', required=True)
    ap.add_argument('--boundary')
    ap.add_argument('--radius', default='1e-25')
    ap.add_argument('--diagnostic', action='store_true')
    ap.add_argument('--output', required=True)
    args = ap.parse_args()
    load = lambda p: json.loads(Path(p).read_text())
    if args.boundary:
        boundary = load(args.boundary)
    else:
        need(args.diagnostic, 'missing rigorous boundary bounds')
        boundary = {k: '0' for k in ('E0_upper', 'E1_upper', 'lift_upper',
                                     'lap_lift_upper', 'dlap_lift_upper')}
    result = bounds(load(args.state), load(args.reference), load(args.spectral),
                    boundary, args.radius, not args.diagnostic)
    dump(args.output, result)
    print(result['stage'], json.dumps(result['bounds']), flush=True)


if __name__ == '__main__':
    main()
