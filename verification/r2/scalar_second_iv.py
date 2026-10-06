#!/usr/bin/env python3
"""Second-library scalar enclosure, adapted from critic C2 (AX6 audit).

No Arb or author analytic stage import. Inputs inherited from Arb are only the
spectral inverse upper bound and eight boundary moment upper bounds. This does
NOT independently implement the spectral matrices or trace quadrature.
"""
import argparse
import json
import re
from pathlib import Path
from mpmath import iv, mp, mpf


def need(ok, reason):
    if not ok:
        raise RuntimeError(reason)


def decimal(text):
    need(isinstance(text, str) and re.fullmatch(
        r'[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?', text),
        'iv: exact decimal required')
    v = iv.mpf(text)
    need(mp.isfinite(mpf(v.a)) and mp.isfinite(mpf(v.b)), 'iv: finite input required')
    return v


def verify(cen, ref, spec, bdry, radius):
    iv.prec = 400
    mp.prec = 400
    m = cen['m']
    need(m == 182 == ref['m'], "iv: m == 182 == ref['m']")
    
    def up(text):
        """upper endpoint of an Arb string '[mid +/- rad]' or plain decimal, as an iv point interval."""
        t = text.strip()
        if t.startswith('['):
            mo = re.match('\\[\\s*([-+0-9.eE]+)\\s*\\+/-\\s*([-+0-9.eE]+)\\s*\\]', t)
            need(mo, 'iv: mo')
            need(I(mo.group(2)) >= 0, 'iv: negative radius')
            return H(decimal(mo.group(1)) + decimal(mo.group(2)))
        return H(decimal(t))
    
    def hi(x):
        return iv.mpf(x).b
    
    def lo(x):
        return iv.mpf(x).a
    
    def H(x):
        return iv.mpf(hi(x))
    I = iv.mpf
    c = [decimal(s) for s in cen['p']]
    a = [decimal(s) for s in cen['a']]
    sc = [decimal(s) for s in cen['scales']]
    need(all(x > 0 for x in sc), 'iv: positive scales')
    cref = [decimal(s) for s in ref['p']]
    J = len(c) - 1
    need(J == 40 and len(a) == len(sc) == 41, 'iv: J == 40 and len(a) == len(sc) == 41')
    p = [(j * m + 1) * c[j] for j in range(J + 1)]
    pref = [(j * m + 1) * cref[j] for j in range(len(cref))]
    rho = c[0]
    A_ = lambda x: abs(x) if not hasattr(x, 'a') else iv.mpf([0 if x.a <= 0 <= x.b else min(abs(x.a), abs(x.b)), max(abs(x.a), abs(x.b))])
    P = sum((A_(p[j]) * 2 ** j for j in range(J + 1)))
    P1 = sum((A_(p[j]) * (j * m) * 2 ** j for j in range(J + 1)))
    Pun = sum((A_(x) for x in p))
    Pref = sum((A_(x) for x in pref))
    chg = sum((A_((p[j] if j <= J else I(0)) - (pref[j] if j < len(pref) else I(0))) for j in range(max(J + 1, len(pref)))))
    dq = (Pun + Pref) * chg
    Bref = up(spec['full_Dirichlet_inverse_L2_norm_upper'])
    need(hi(Bref * dq) < 1, 'iv: hi(Bref * dq) < 1')
    BL2 = Bref / (1 - Bref * dq)
    h = 1
    while h * m <= lo(rho):
        h += 1
    need(h * m > hi(rho), 'iv: h * m > hi(rho)')
    dh = I(h * m) ** 2 - rho * rho
    need(lo(dh) > 0, 'iv: lo(dh) > 0')
    Vs = P * P - rho * rho
    beta = Vs / dh
    need(hi(beta) < 1, 'iv: hi(beta) < 1')
    Wl = 1 + 2 * sum((2 ** j for j in range(1, h)))
    low = Wl * (1 + P * P * BL2) / (2 * iv.sqrt(2))
    BY = (low + 1 / dh) / (1 - beta)
    BN = 1 + P * P * BY
    R = iv.exp(iv.log(2) / m)
    R4 = iv.exp(iv.log(4) / m)
    S = sum((A_(c[j]) * 4 ** j for j in range(1, J + 1)))
    need(lo(rho - 2 * S) > 0, 'iv: lo(rho - 2 * S) > 0')
    d = (2 * rho * S + S * S) / rho
    eps = iv.log((rho + S) / (rho - S)) / 2
    P4 = sum((A_(c[j]) * (j * m + 1) * 4 ** j for j in range(J + 1)))
    MU = I(0)
    MdU = I(0)
    for j in range(J + 1):
        n = j * m
        if n > hi(rho + d):
            x = I(n) / (rho + d)
            al = iv.log(x + iv.sqrt(x * x - 1))
        else:
            al = I(0)
        ch = (iv.exp(al) + iv.exp(-al)) / 2
        sh = (iv.exp(al) - iv.exp(-al)) / 2
        b = A_(a[j]) / sc[j] * 4 ** j * iv.exp(n * eps - n * al + (rho + d) * sh + d * ch)
        MU += b
        MdU += b * iv.exp(al)
    Ub = 3 * MU
    Db = 3 * R * R4 * R4 * P4 * iv.exp(eps) * MdU
    M = Ub + 2 * Db
    D = [up(s) for s in bdry['D_moments_upper']]
    N = [up(s) for s in bdry['N_moments_upper']]
    L0 = D[0] + N[0] / 2
    L1 = 2 * R * (D[2] + N[1])
    LD = 2 * (D[2] + N[1])
    LD1 = 2 * R * (D[3] + N[2])
    E0 = LD + P * P * L0
    E1 = LD1 + R * P * P1 * L0 + P * P * L1
    D0 = P * P * Ub + LD
    D1 = R * P * P1 * Ub + P * P * Db + LD1
    dstar = rho * rho - Vs - E0
    need(lo(dstar) > 0, 'iv: lo(dstar) > 0')
    I0 = 1 / dstar
    I1 = I0 + I0 * I0 * (2 * P * P1 + 2 * R * E1)
    Bp = (P + P1) * I1 * BN
    Bg = BN + 2 * (R * D1 + D0) * I1 * BN
    Astar = Bg + Bp
    hg = rho - (P - rho)
    need(lo(hg) > 0, 'iv: lo(hg) > 0')
    J0 = 1 / hg
    J1 = J0 + P1 * J0 * J0
    Z = 2 * Astar * (R * E1 + E0) * J1
    Y0 = Astar * E0
    kap = I(1) / 4
    Uc = Ub + L0
    C2 = Astar * (2 * P * kap + Uc)
    C3 = Astar * kap
    r = I(radius)
    poly = Y0 + (Z - 1) * r + C2 * r * r + C3 * r ** 3
    contr = Z + 2 * C2 * r + 3 * C3 * r * r
    first = sum((A_(x) for x in p[1:]))
    second = sum((A_(p[j]) * (j * m) for j in range(1, J + 1)))
    Lg = rho - first - r
    Ug = second + I(m + 1) * r / 2
    curv = 1 - Ug / Lg
    joint = rho - first - second - I(m + 1) * r / 2
    nond = A_(c[1]) - r / (2 * (m + 1))
    univ = rho - (P - rho) - r

    quantities = {'B_L2': BL2, 'B_Y': BY, 'B_trace': BN, 'U_field': Ub,
                  'dU_field': Db, 'E0': E0, 'E1': E1, 'A': Astar,
                  'Y': Y0, 'Z': Z, 'C2': C2, 'C3': C3,
                  'radii_polynomial': poly, 'contraction': contr,
                  'four_C2_Y': 4*C2*Y0, 'curvature_lower': curv,
                  'convexity_margin_lower': joint, 'non_disc_c1_lower': nond,
                  'analytic_univalence_margin_lower': univ}
    for key, value in quantities.items():
        need(mp.isfinite(mpf(lo(value))) and mp.isfinite(mpf(hi(value))), 'iv nonfinite: '+key)
    need(hi(poly) < 0, 'iv radii polynomial')
    need(hi(contr) <= decimal('1e-5'), 'iv contraction target')
    need(hi(4*C2*Y0) <= decimal('1e-10'), 'iv product target')
    need(lo(curv) > decimal('0.2998879') and lo(joint) > decimal('243.3264')
         and lo(nond) > decimal('0.01645') and lo(univ) > 808 and lo(Lg) > 0,
         'iv whole-ball geometry')
    # Endpoint strings below are display values, not interval inputs to another gate.
    return {'status': 'SECOND_LIBRARY_SCALAR_GATES_PASS', 'library': 'mpmath.iv',
            'precision_bits': 400,
            'scope': 'scalar chain only; spectral inverse and boundary moments inherited from Arb',
            'bounds_display': {k: [mp.nstr(mpf(lo(v)), 55), mp.nstr(mpf(hi(v)), 55)]
                               for k, v in quantities.items()}}


def main():
    ap = argparse.ArgumentParser()
    for name in ('candidate', 'reference', 'spectral', 'boundary', 'output'):
        ap.add_argument('--'+name, required=True)
    ap.add_argument('--radius', required=True)
    args = ap.parse_args()
    load = lambda p: json.loads(Path(p).read_text())
    result = verify(load(args.candidate), load(args.reference), load(args.spectral),
                    load(args.boundary), args.radius)
    Path(args.output).write_text(json.dumps(result, indent=2)+'\n')
    print(result['status'], flush=True)


if __name__ == '__main__':
    main()
