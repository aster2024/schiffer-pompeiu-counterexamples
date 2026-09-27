#!/usr/bin/env python3
"""Exact-rational and Arb checks of two identities used by the R6 certificate.

(I)  K g = (1-r^2)^2 Q   (far_shape_r6.py, Section 5 of FINAL.md)
     I.1 exact (fmpq_poly), for every n = 2 mod 4 with n <= --k-nmax and
         1 <= s <= --k-smax: the radial polynomial of K S_{n,s} from (3.1)
         equals (1-q)^2 P_{s-1}^{(2,n)}(2q-1)/(4s(s+1)), q = r^2; it solves
         Delta(K S_{n,s}) = S_{n,s}; and both boundary traces vanish.
     I.2 exact: the Jacobi connection recurrence of far_shape_r6.q_of_g, run
         from its own byte code with exact rationals instead of Arb balls,
         reproduces P_{s-1}^{(2,n)}/(4s(s+1)) in the P^{(0,n)} basis for all
         centre indices n <= 58, s <= 30, with positive coefficients.
     I.3 Arb (128 bits): for the frozen centre g°, every coefficient ball of
         K g° - (1-|w|^2)^2 Q, computed with the verifier's own disk-polynomial
         algebra, contains 0; the weighted norm of the enclosure is reported.

(II) The principal shape column (5.1), for j = 1 mod 4, j >= 5:
     T_j = b^3 {(j+2)S_{j+1,0}/2 - (j-2)S_{j-3,0}/2 - j S_{j-3,1}/(j+1)
                - S_{j-3,2}/(j+1)}  =  DH(0, b w) e_j.
     II.1 exact: r^{m+4} = r^m[(m+1)/(m+3)P_0 + 2/(m+4)P_1
                 + 2/((m+3)(m+4))P_2]^{(0,m)}(2r^2-1) for every m = j-3.
     II.2 exact: the verifier's own column algebra (zw/hol/sine/unsine,
          exactly as in shape_tail_interval_r6.py) at the ball psi = w, g = 0,
          run in exact rationals, equals principal(j, 1) coefficient by
          coefficient (T_j is homogeneous of degree 3 in b).
     II.3 Arb: the same comparison with b = c°_1 (128 bits); each coefficient
          ball of the difference contains 0.
     II.4 Arb pointwise: the disk-polynomial expansion (5.1) agrees with the
          direct formula b^3[j Re(w^{j-1}) Im(w^2) + Im(w^{j+1})] at sample points.

Every check is fail-closed; the script prints IDENTITIES VERIFIED only if all
of them pass, and writes a JSON receipt.
"""
import argparse
import hashlib
import json
import types
from math import comb, factorial
from pathlib import Path

from flint import arb, ctx, fmpq, fmpq_poly

import verify_r6 as V
import far_shape_r6
import shape_tail_interval_r6


def need(ok, reason):
    if not ok:
        raise RuntimeError('IDENTITY CHECK FAILED: ' + reason)


U = fmpq_poly([0, 1])          # u = q - 1 = r^2 - 1
Q_OF_U = fmpq_poly([1, 1])     # q = u + 1
_JAC = {}
_QPOW = [fmpq_poly([1])]


def qpow(k):
    """q^k = (u+1)^k, cached."""
    while len(_QPOW) <= k:
        _QPOW.append(_QPOW[-1] * Q_OF_U)
    return _QPOW[k]


def pochhammer(a, k):
    out = 1
    for i in range(k):
        out *= a + i
    return out


def jac(s, alpha, n):
    """P_s^{(alpha,n)}(2q-1) as an exact polynomial in u = q-1 (DLMF 18.5.7).

    Coefficient of u^l: (alpha+1)_s (s+alpha+n+1)_l / ((s-l)! (alpha+1)_l l!).
    """
    key = (s, alpha, n)
    if key not in _JAC:
        if s < 0:
            _JAC[key] = fmpq_poly([])
        else:
            num = pochhammer(alpha + 1, s)
            coeffs = [fmpq(num * pochhammer(s + alpha + n + 1, l),
                           factorial(s - l) * pochhammer(alpha + 1, l) * factorial(l))
                      for l in range(s + 1)]
            _JAC[key] = fmpq_poly(coeffs)
    return _JAC[key]


def jac_paper_q(s, n):
    """P_s^{(0,n)}(2q-1) from the explicit q-power formula a^{(n)}_{s,k} of the
    R^4 paper (eq. jacobi), re-expressed in u, used as an independent check."""
    out = fmpq_poly([])
    for k in range(s + 1):
        a = fmpq((-1) ** (s - k) * factorial(n + s + k),
                 factorial(k) * factorial(s - k) * factorial(n + k))
        out += a * qpow(k)
    return out


def check_jacobi_formulas(nmax, smax):
    """Two independent formulas agree, the (2,n) family is orthogonal in
    L^2(q^n (1-q)^2 dq) and both families are normalised by P(1)=binom(s+a,s)."""
    count = 0
    for n in range(2, nmax + 1, 4):
        for s in range(0, smax + 1):
            need(jac(s, 0, n) == jac_paper_q(s, n), f'Jacobi formulas differ at n={n}, s={s}')
            need(jac(s, 0, n)(fmpq(0)) == 1 and jac(s, 2, n)(fmpq(0)) == comb(s + 2, s),
                 f'Jacobi normalisation at n={n}, s={s}')
            count += 1
    # Orthogonality of P^{(2,n)}_s against q^i, i < s, in L^2(q^n(1-q)^2 dq).
    for n in (2, 6, 58):
        for s in range(1, 12):
            for i in range(s):
                integrand = qpow(n + i) * U * U * jac(s, 2, n)
                F = integrand.integral()
                need(F(fmpq(0)) - F(fmpq(-1)) == 0, f'P^(2,n) orthogonality n={n}, s={s}, i={i}')
    return count


def radial_K(n, s):
    d = n + 2 * s
    return (jac(s - 1, 0, n) * fmpq(1, 4 * d * (d + 1))
            - jac(s, 0, n) * fmpq(1, 2 * d * (d + 2))
            + jac(s + 1, 0, n) * fmpq(1, 4 * (d + 1) * (d + 2)))


def check_K_exact(nmax, smax):
    count = 0
    for n in range(2, nmax + 1, 4):
        for s in range(1, smax + 1):
            R = radial_K(n, s)
            target = U * U * jac(s - 1, 2, n) * fmpq(1, 4 * s * (s + 1))
            need(R == target, f'(1-r^2)^2 factorisation of K S_({n},{s})')
            # Delta(r^n R(r^2) sin n theta) = r^n [4 q R'' + 4(n+1) R'] sin n theta.
            R1 = R.derivative()
            lap = 4 * Q_OF_U * R1.derivative() + 4 * (n + 1) * R1
            need(lap == jac(s, 0, n), f'Delta K S_({n},{s}) != S_({n},{s})')
            # Traces at r = 1 (u = 0): value R(1) and radial derivative n R(1) + 2 R'(1).
            need(R(fmpq(0)) == 0 and n * R(fmpq(0)) + 2 * R1(fmpq(0)) == 0,
                 f'boundary traces of K S_({n},{s})')
            # Column sum of |coefficients| in (3.1) is 1/(d(d+2)).
            d = n + 2 * s
            need(fmpq(1, 4 * d * (d + 1)) + fmpq(1, 2 * d * (d + 2)) + fmpq(1, 4 * (d + 1) * (d + 2))
                 == fmpq(1, d * (d + 2)), f'column sum of K at ({n},{s})')
            count += 1
    return count


def exact_acc(out, key, value):
    # NB: in python-flint 0.9.0, fmpq(0).is_zero() returns False, so compare
    # with 0 explicitly. (The verifier itself only uses arb.is_zero().)
    if value != 0:
        out[key] = out.get(key, fmpq(0)) + value


def check_q_of_g_exact(M, S):
    """Run far_shape_r6.q_of_g's byte code with exact rationals."""
    glb = dict(far_shape_r6.q_of_g.__globals__)
    glb.update(arb=fmpq, Z=fmpq(0), acc=exact_acc)
    q_exact = types.FunctionType(far_shape_r6.q_of_g.__code__, glb, 'q_of_g_exact')
    count = 0
    for n in range(2, M + 1, 4):
        for s in range(1, S + 1):
            Q = q_exact({(n, s): fmpq(1)}, M, S)
            need(all(m == n for (m, t) in Q), 'q_of_g changed the angular frequency')
            need(all(v > 0 for v in Q.values()), f'nonpositive connection coefficient at ({n},{s})')
            need(max(t for (m, t) in Q) == s - 1, f'q_of_g degree at ({n},{s})')
            poly = fmpq_poly([])
            for (m, t), v in Q.items():
                poly += v * jac(t, 0, n)
            need(poly == jac(s - 1, 2, n) * fmpq(1, 4 * s * (s + 1)),
                 f'q_of_g connection coefficients at ({n},{s})')
            count += 1
    return count


def one_minus_abs2(f):
    """Multiply a full disk-polynomial expansion by (1-|w|^2)."""
    return V.add(f, V.scale(V.zw(V.zw(f), True), -arb(1)))


def check_Kg_arb(data, rho):
    rows, g, c, p, Vc, F = V.make_algebra(data)
    Q = far_shape_r6.q_of_g(g, data['M'], data['S'])
    KQ = V.unsine(one_minus_abs2(one_minus_abs2(V.sine(Q))))
    Kg = V.K(g)
    diff = V.add(Kg, V.scale(KQ, -arb(1)))
    need(set(diff) <= set(Kg) | set(KQ), 'support bookkeeping')
    for key, v in diff.items():
        need(v.is_finite() and v.contains(0), f'Kg - (1-r^2)^2 Q ball excludes 0 at {key}')
    dnorm = V.norm(diff, rho)
    kg_norm = V.norm(Kg, rho)
    need(dnorm.is_finite() and dnorm.upper() < arb('1e-25'), 'Kg - (1-r^2)^2 Q enclosure too wide')
    return {'coefficients': len(diff), 'diff_norm_upper': dnorm.upper().str(6),
            'Kg_norm': kg_norm.str(20), 'Q_norm': V.norm(Q, rho).str(20)}


def check_radial_51(jmax):
    count = 0
    for j in range(5, jmax + 1, 4):
        m = j - 3
        rhs = (fmpq(m + 1, m + 3) * jac(0, 0, m) + fmpq(2, m + 4) * jac(1, 0, m)
               + fmpq(2, (m + 3) * (m + 4)) * jac(2, 0, m))
        need(Q_OF_U * Q_OF_U == rhs, f'radial expansion r^(m+4) at m={m}')
        count += 1
    return count


def ball_columns(b, one, jmax):
    """Yield (j, DF_j) for the ball x_b = (g=0, psi=b w), using exactly the
    column algebra of shape_tail_interval_r6.py."""
    cb = {1: b}
    pb = {0: b}
    Vb = V.hcoeff(cb)
    A = V.shift(V.hol(V.sine(Vb), pb, True), 4)
    B = V.shift(V.hol_psi(V.abs2({(1, 0): one}, pb), cb), 4)
    half = one / 2
    two = one + one
    for j in range(5, jmax + 1, 4):
        DF = V.scale(V.unsine(V.add(V.scale(A, one * j), V.scale(B, -half))), two)
        yield j, DF
        A = V.shift(A, 4)
        B = V.shift(B, 4)


def check_principal_exact(jmax):
    saved = V.Z
    V.Z = fmpq(0)          # the verifier's accumulators start from V.Z
    try:
        count = 0
        for j, DF in ball_columns(fmpq(1), fmpq(1), jmax):
            T = shape_tail_interval_r6.principal(j, fmpq(1))
            diff = V.add(DF, V.scale(T, fmpq(-1)))
            need(all(v == 0 for v in diff.values()), f'principal column (5.1) at j={j}')
            need(set(k for k, v in DF.items() if v != 0) == set(T), f'support of T_{j}')
            count += 1
    finally:
        V.Z = saved
    return count


def check_principal_arb(b, jmax):
    worst = arb(0)
    count = 0
    for j, DF in ball_columns(b, arb(1), jmax):
        T = shape_tail_interval_r6.principal(j, b)
        diff = V.add(DF, V.scale(T, -arb(1)))
        for key, v in diff.items():
            need(v.is_finite() and v.contains(0), f'Arb principal column j={j} at {key}')
            worst = max(worst, abs(v).upper())
        count += 1
    return count, worst


def eval_S(n, s, r, th):
    """S_{n,s}(r, theta) in Arb, from the exact Jacobi polynomial in u = r^2 - 1."""
    u = r * r - 1
    val = arb(0)
    for cf in reversed(jac(s, 0, n).coeffs()):
        val = val * u + arb(cf.p) / arb(cf.q)
    return r ** n * val * (n * th).sin()


def check_principal_pointwise(b, js, points):
    worst = arb(0)
    for j in js:
        T = shape_tail_interval_r6.principal(j, b)
        for r, th in points:
            series = sum((v * eval_S(n, s, r, th) for (n, s), v in T.items()), arb(0))
            direct = b ** 3 * (j * r ** (j - 1) * ((j - 1) * th).cos() * r * r * (2 * th).sin()
                               + r ** (j + 1) * ((j + 1) * th).sin())
            d = abs(series - direct)
            need(d.is_finite() and (series - direct).contains(0), f'pointwise (5.1) j={j}')
            worst = max(worst, d.upper())
    return worst


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--bits', type=int, default=128)
    ap.add_argument('--k-nmax', type=int, default=202)
    ap.add_argument('--k-smax', type=int, default=80)
    ap.add_argument('--jmax-exact', type=int, default=2001)
    ap.add_argument('--jmax-arb', type=int, default=1001)
    ap.add_argument('--output', default='identities_r6.json')
    args = ap.parse_args()
    ctx.prec = args.bits
    ctx.threads = 1
    need(args.bits >= 96, 'at least 96 bits')
    for name in ('center_r6.json', 'far_shape_r6.py', 'shape_tail_interval_r6.py'):
        need(hashlib.sha256(Path(name).read_bytes()).hexdigest() == V.EXPECTED_SHA256[name],
             'input differs from the frozen certificate: ' + name)
    data = json.loads(Path('center_r6.json').read_text())
    need(data['M'] == 58 and data['S'] == 30, 'unexpected centre')
    rho = arb(11) / 10
    b = V.dyadic(data['c'][0])

    out = {'status': None, 'bits': args.bits}
    out['jacobi_formula_pairs'] = check_jacobi_formulas(args.k_nmax, args.k_smax + 1)
    print('Jacobi formulas agree:', out['jacobi_formula_pairs'], 'pairs', flush=True)
    out['K_exact'] = {'n_max': args.k_nmax, 's_max': args.k_smax,
                      'pairs': check_K_exact(args.k_nmax, args.k_smax)}
    print('I.1 exact K S_{n,s} = (1-r^2)^2 r^n P_{s-1}^{(2,n)} sin/(4s(s+1)):',
          out['K_exact']['pairs'], 'pairs', flush=True)
    out['q_of_g_exact_pairs'] = check_q_of_g_exact(data['M'], data['S'])
    print('I.2 exact q_of_g connection coefficients:', out['q_of_g_exact_pairs'], 'pairs', flush=True)
    out['Kg_arb'] = check_Kg_arb(data, rho)
    print('I.3 Arb ||Kg - (1-r^2)^2 Q||_rho <=', out['Kg_arb']['diff_norm_upper'],
          '(||Kg||_rho =', out['Kg_arb']['Kg_norm'][:24] + ')', flush=True)

    out['radial_51_exact'] = {'j_max': args.jmax_exact, 'count': check_radial_51(args.jmax_exact)}
    print('II.1 exact radial expansion for j = 5..', args.jmax_exact, ':',
          out['radial_51_exact']['count'], 'values', flush=True)
    out['principal_exact'] = {'j_max': args.jmax_exact, 'count': check_principal_exact(args.jmax_exact)}
    print('II.2 exact verifier-algebra ball column == (5.1):', out['principal_exact']['count'],
          'values of j', flush=True)
    cnt, worst = check_principal_arb(b, args.jmax_arb)
    out['principal_arb'] = {'j_max': args.jmax_arb, 'count': cnt, 'max_abs_diff_upper': worst.str(6)}
    print('II.3 Arb ball column == (5.1) for', cnt, 'values of j; max |diff| <=', worst.str(6), flush=True)
    pts = [(arb(1), arb(3) / 10), (arb(7) / 10, arb(11) / 10), (arb(99) / 100, arb(2)),
           (arb(1) / 2, arb(-5) / 7)]
    pw = check_principal_pointwise(b, (5, 9, 61, 101, 397, 1001), pts)
    out['principal_pointwise_arb'] = {'j': [5, 9, 61, 101, 397, 1001], 'points': len(pts),
                                      'max_abs_diff_upper': pw.str(6)}
    print('II.4 Arb pointwise (5.1) vs direct derivative: max |diff| <=', pw.str(6), flush=True)
    out['status'] = 'IDENTITIES_VERIFIED'
    Path(args.output).write_text(json.dumps(out, indent=2) + '\n')
    print('IDENTITIES VERIFIED', flush=True)
    print('receipt:', args.output, flush=True)


if __name__ == '__main__':
    main()
