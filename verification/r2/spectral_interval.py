#!/usr/bin/env python3
"""Arb Dirichlet spectral inverse component. No Schiffer existence claim.

Floating numbers only propose exact dyadic centers and a congruence matrix.
Every accepted numerical inequality is rebuilt with Arb.
"""
import argparse
import hashlib
import json
import math
import re
import time
from pathlib import Path

from flint import arb, arb_mat, ctx
ctx.threads = 1


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def exact_decimal(text):
    """Enclose an exact decimal rational, never a float or an Arb ball string."""
    need(isinstance(text, str) and re.fullmatch(
        r'[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?', text),
        'expected finite exact decimal string')
    value = arb(text)
    upper(value)
    return value


def upper(x):
    need(isinstance(x, arb) and x.is_finite(), 'nonfinite/non-Arb bound')
    y = x.upper()
    need(y.is_finite(), 'nonfinite upper endpoint')
    return y


def umax(xs):
    result = None
    for x in xs:
        x = upper(x)
        if result is None or x > result:
            result = x
    need(result is not None, 'empty maximum')
    return result


def dyad(x):
    """Treat a numerical proposal as an exact binary rational."""
    need(math.isfinite(x), 'nonfinite numerical proposal')
    a, b = x.as_integer_ratio()
    return arb(a)/b


def parse_hex(s):
    need(isinstance(s, str) and s.startswith(('0x', '-0x')), 'not a dyadic hex')
    return dyad(float.fromhex(s))


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def shape_digest(state):
    return hashlib.sha256(json.dumps({'m': state['m'], 'p': state['p']},
                                    sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def dump(path, data):
    Path(path).write_text(json.dumps(data, indent=2)+'\n')


def hull(mid, rad):
    out = arb(mid, upper(rad))
    upper(out)
    return out


def bessel_root(n, proposal):
    ctx.prec = 2048
    x = dyad(proposal)
    for _ in range(7):
        v = x.bessel_j(n)
        deriv = -x.bessel_j(1) if n == 0 else x.bessel_j(n-1)-n*v/x
        need(abs(deriv) > 0, 'Bessel Newton derivative contains zero')
        x = (x-v/deriv).mid()
    eps = arb(1)/2**220
    left, right = x-eps, x+eps
    need(left > 1 and right-left < arb(1)/10, 'invalid Bessel bracket width/location')
    need(left.bessel_j(n)*right.bessel_j(n) < 0, 'Bessel root sign bracket failed')
    norm = abs(hull(x.bessel_j(n+1), eps))
    need(norm > 0, 'zero Bessel normalization')
    return {'n': n, 'left': left.str(95), 'right': right.str(95),
            'mid': x.str(100), 'norm': norm.str(70)}


def gauss_grid(count):
    import numpy as np
    from scipy.special import roots_legendre
    proposals, _ = roots_legendre(count)
    ctx.prec = 256
    eps = arb(1)/2**180
    grid = []
    previous = arb(-1)
    for proposal in proposals:
        x = dyad(float(proposal))
        for _ in range(5):
            p = x.legendre_p(count)
            dp = count*(x*p-x.legendre_p(count-1))/(x*x-1)
            need(abs(dp) > 0, 'Legendre derivative contains zero')
            x = (x-p/dp).mid()
        left, right = x-eps, x+eps
        need(previous < left < right < 1, 'Gauss nodes overlap or leave interval')
        need(left.legendre_p(count)*right.legendre_p(count) < 0,
             'Legendre sign bracket failed')
        previous = right
        node = hull(x, eps)
        p = x.legendre_p(count)
        dp = count*(x*p-x.legendre_p(count-1))/(x*x-1)
        # Max |P_N''| on [-1,1] is its endpoint value.
        second = arb((count-1)*count*(count+1)*(count+2))/8
        dp = hull(dp, eps*second)
        weight = 2/((1-node*node)*dp*dp)
        need(weight > 0, 'Gauss weight not positive')
        r = (1+node)/2
        wr = weight*r/2
        grid.append([r.str(65), wr.str(65)])
    need(len(grid) == count, 'Gauss root coverage failed')
    return grid


def prepare(state_path, window, nodes, output):
    import numpy as np
    from scipy.special import jn_zeros
    state = json.loads(Path(state_path).read_text())
    m = state['m']
    need(m > 0 and len(state['p']) > 1, 'bad candidate')
    ctx.prec = 2048
    rho = exact_decimal(state['p'][0])
    lo, hi = rho-window, rho+window
    need(lo > 1 and window > 0, 'bad spectral window')
    states, coverage = [], []
    maxj = int(float(hi)/m)
    for j in range(maxj+1):
        n = m*j
        proposals = jn_zeros(n, int(float(hi)/np.pi)+6)
        below = [k for k, z in enumerate(proposals) if z < float(lo)]
        above = [k for k, z in enumerate(proposals) if z > float(hi)]
        need(bool(above), 'no upper root guard')
        start = below[-1] if below else 0
        stop = above[0]+1
        roots = [bessel_root(n, float(z)) for z in proposals[start:stop]]
        for root in roots:
            left, right = arb(root['left']), arb(root['right'])
            if lo < left and right < hi:
                states.append(dict(root, j=j))
            else:
                need(right < lo or hi < left, 'root straddles spectral cutoff')
        if start == 0:
            need(n > 0, 'no lower zero guard for order zero')
            end = arb(roots[0]['right'])
            qmax = 1-(arb(n*n)-arb(1)/4)/(end*end)
            need(qmax > 0 and end-n < 3*arb.pi()/(2*qmax.sqrt()),
                 'first Bessel zero coverage failed')
        else:
            need(arb(roots[0]['right']) < lo, 'missing lower Bessel guard')
        for previous, following in zip(roots, roots[1:]):
            left, mid = arb(previous['left']), arb(previous['right'])
            left2, right = arb(following['left']), arb(following['right'])
            endpoint = left if n == 0 else right
            qmax = 1-(arb(n*n)-arb(1)/4)/(endpoint*endpoint)
            sep = arb.pi()/qmax.sqrt()
            need(mid < left2 and right-left < 2*sep,
                 'unexcluded Bessel zero between guards')
        need(hi < arb(roots[-1]['left']), 'missing upper Bessel guard')
        coverage.append({'j': j, 'n': n, 'roots': roots,
                         'starts_at_first_root': start == 0})
    need((maxj+1)*m > hi, 'higher angular order not excluded')
    grid = gauss_grid(nodes)
    result = {'stage': 'CERTIFIED_SPECTRAL_INPUTS', 'm': m,
              'state': state, 'source_sha256': sha(state_path),
              'window': window, 'nodes': nodes, 'states': states,
              'coverage': coverage, 'gauss_grid': grid}
    from coverage_checks import check_root_coverage
    result['sector_counts'] = check_root_coverage(result)
    dump(output, result)
    print('SPECTRAL_INPUTS', len(states), 'states', nodes, 'nodes', flush=True)


def basis_leaf(inputs, start, stop, output):
    data = json.loads(Path(inputs).read_text())
    need(data['stage'] == 'CERTIFIED_SPECTRAL_INPUTS', 'wrong input stage')
    need(0 <= start < stop <= len(data['states']), 'bad basis coverage')
    ctx.prec = 2048
    rs = [arb(pair[0]) for pair in data['gauss_grid']]
    rows, errors = [], []
    for state in data['states'][start:stop]:
        n = state['n']
        left, right = arb(state['left']), arb(state['right'])
        a = hull((left+right)/2, (right-left)/2)
        norm = arb(state['norm'])
        row, errs = [], []
        for r in rs:
            x = a*r
            # Real |J_n'|<=1, so evaluate at an exact midpoint and add radius.
            val = arb(2).sqrt()*hull(x.mid().bessel_j(n), x.rad())/norm
            mid = float(val.mid())
            err = upper(abs(val-dyad(mid)))
            row.append(mid.hex())
            errs.append(err)
        err = umax(errs)
        need(err < arb(1)/10**11, 'basis binary64 enclosure too wide')
        rows.append(row)
        errors.append(err.str(50))
    dump(output, {'stage': 'CERTIFIED_BESSEL_BASIS',
                  'input_sha256': sha(inputs), 'start': start, 'stop': stop,
                  'basis': rows, 'row_error': errors})
    print('BASIS', start, stop, 'max error', umax([arb(x) for x in errors]).str(8), flush=True)


def matrix_inf(mat):
    return umax(sum((abs(mat[i, j]) for j in range(mat.ncols())), arb(0))
                for i in range(mat.nrows()))


def assemble(inputs, leaves, output, coverage_path):
    import numpy as np
    from scipy.linalg import eigh
    ctx.prec = 128
    data = json.loads(Path(inputs).read_text())
    need(data['stage'] == 'CERTIFIED_SPECTRAL_INPUTS', 'wrong input stage')
    from coverage_checks import check_root_coverage, check_scan_binding
    counts = check_root_coverage(data)
    scan = json.loads(Path(coverage_path).read_text())
    check_scan_binding(scan, data)
    ctx.prec = 128
    states = data['states']
    size, nodes = len(states), data['nodes']
    pieces = [json.loads(Path(p).read_text()) for p in leaves]
    pieces.sort(key=lambda x: x['start'])
    nxt = 0
    basis = []
    for piece in pieces:
        need(piece['stage'] == 'CERTIFIED_BESSEL_BASIS'
             and piece['input_sha256'] == sha(inputs), 'basis provenance mismatch')
        need(piece['start'] == nxt < piece['stop'] <= size, 'basis coverage gap/overlap')
        count = piece['stop']-piece['start']
        need(len(piece['basis']) == len(piece['row_error']) == count,
             'basis row count mismatch')
        for row, error in zip(piece['basis'], piece['row_error']):
            need(len(row) == nodes, 'basis node count mismatch')
            e = arb(error)
            need(e >= 0 and e < arb(1)/10**11, 'bad basis error')
            basis.append([hull(parse_hex(x), e) for x in row])
        nxt = piece['stop']
    need(nxt == size, 'incomplete basis coverage')
    rho = arb(data['state']['p'][0])
    m = data['m']
    p = [(m*j+1)*arb(c) for j, c in enumerate(data['state']['p'])]
    grid = [[arb(x) for x in row] for row in data['gauss_grid']]
    q = [[arb(0) for _ in range(nodes)] for _ in p]
    for t, (r, _) in enumerate(grid):
        powers = [arb(1), r**m]
        for _ in range(2, 2*len(p)-1):
            powers.append(powers[-1]*powers[1])
        for h in range(len(p)):
            for k in range(len(p)-h):
                if k or h:
                    q[h][t] += p[k+h]*p[k]*powers[2*k+h]
    q2 = [[arb(0) for _ in range(nodes)] for _ in range(2*len(p)-1)]
    for h in range(len(q2)):
        for k in range(-len(q)+1, len(q)):
            if abs(h-k) < len(q):
                for t in range(nodes):
                    q2[h][t] += q[abs(k)][t]*q[abs(h-k)][t]
    # Analytic Gauss error: ellipse parameter 2 in x=2r-1.
    # |r|<=9/8, |Im r|<=3/8 and |J_n(ar)|<=exp(a |Im r|).
    ellipse_r = arb(9)/8
    ellipse_p = sum((abs(v)*ellipse_r**(m*j) for j, v in enumerate(p)), arb(0))
    ellipse_v = ellipse_p*ellipse_p+rho*rho
    norm_min = None
    for state in states:
        v = arb(state['norm']).lower()
        if norm_min is None or v < norm_min:
            norm_min = v
    need(norm_min > 0, 'spectral normalization failed')
    max_zero = umax(arb(s['right']) for s in states)
    majorant = 8*(1+ellipse_r)*(1+ellipse_v*ellipse_v)*(max_zero*arb(3)/4).exp()/norm_min**2
    quadrature_error = 16*majorant/arb(2)**(2*nodes)
    need(quadrature_error < arb(1)/10**20, 'Gauss analytic error too large')
    qe = hull(arb(0), quadrature_error)
    group = {}
    for i, state in enumerate(states):
        group.setdefault(state['j'], []).append(i)
    blocks = {j: arb_mat([basis[i] for i in ii]) for j, ii in group.items()}

    def project(potential):
        result = arb_mat(size, size)
        for j, ii in group.items():
            for k, kk in group.items():
                if k < j:
                    continue
                v = [arb(0) for _ in range(nodes)]
                if j == 0 and k == 0:
                    v = potential[0]
                elif j == 0:
                    if k < len(potential):
                        v = [arb(2).sqrt()*x for x in potential[k]]
                else:
                    if abs(k-j) < len(potential):
                        v = [a+b for a, b in zip(v, potential[abs(k-j)])]
                    if j+k < len(potential):
                        v = [a+b for a, b in zip(v, potential[j+k])]
                weighted = arb_mat([[basis[col][t]*grid[t][1]*v[t]
                                     for col in kk] for t in range(nodes)])
                block = blocks[j]*weighted
                for a, i in enumerate(ii):
                    for b, col in enumerate(kk):
                        value = block[a, b]+qe
                        upper(value)
                        result[i, col] = value
                        result[col, i] = value
        return result

    v = project(q)
    v2 = project(q2)
    print('MATRICES_ENCLOSED', size, flush=True)
    h = -v
    for i, state in enumerate(states):
        a = hull((arb(state['left'])+arb(state['right']))/2,
                 (arb(state['right'])-arb(state['left']))/2)
        h[i, i] += a*a-rho*rho
    hc = np.array([[float(h[i, j].mid()) for j in range(size)] for i in range(size)])
    evals, eigvecs = eigh(hc)
    d = [dyad(float(x)) for x in evals]
    c = arb_mat([[dyad(float(x)) for x in row] for row in eigvecs])
    ct = c.transpose()
    eye = arb_mat([[int(i == j) for j in range(size)] for i in range(size)])
    c_defect = matrix_inf(ct*c-eye)
    need(c_defect < arb(1)/100, 'congruence matrix not certified invertible')
    diag_abs = [abs(x) for x in d]
    need(all(x > 0 for x in diag_abs), 'zero approximate eigenvalue')
    inv_sqrts = [1/x.sqrt() for x in diag_abs]
    e = ct*h*c
    for i in range(size):
        e[i, i] -= d[i]
    for i in range(size):
        for j in range(size):
            e[i, j] *= inv_sqrts[i]*inv_sqrts[j]
    finite_defect = matrix_inf(e)
    gram = v2-v*v
    gc = gram*c
    trace = arb(0)
    for j in range(size):
        trace += sum((c[i, j]*gc[i, j] for i in range(size)), arb(0))/diag_abs[j]
    perturb = sum((abs(x) for x in p[1:]), arb(0))
    potential_bound = 2*rho*perturb+perturb*perturb
    window = data['window']
    gap = window*(2*rho-window)-potential_bound
    need(gap > 0 and trace >= 0, 'nonpositive complement gap or Gram trace')
    eta = upper(trace)/gap.lower()
    need(finite_defect+eta < 1, 'infinite Dirichlet Schur gate failed')
    min_d = diag_abs[0]
    for x in diag_abs[1:]:
        if x < min_d:
            min_d = x
    s_inverse = (1+c_defect)/(min_d*(1-finite_defect-eta))
    full_inverse = (1+(potential_bound/gap)**2)*s_inverse+1/gap
    result = {
        'stage': 'INVARIANT_SECTOR_DIRICHLET_INVERSE',
        'scope': 'real D_m invariant even L2(unit disc), Dirichlet boundary only',
        'full_schiffer_proof': False, 'input_sha256': sha(inputs),
        'reference_shape_sha256': shape_digest(data['state']),
        'source_sha256': sha(__file__), 'basis_sha256': {str(p): sha(p) for p in leaves},
        'm': m, 'window': window, 'states': size, 'gauss_nodes': nodes,
        'sector_counts': counts,
        'independent_coverage_sha256': sha(coverage_path),
        'independent_sector_counts': scan['counts'],
        'coverage_methods': ['root_brackets_and_sturm_guards', 'quarter_grid_certified_signs'],
        'gauss_error_upper': upper(quadrature_error).str(45),
        'congruence_gram_defect_upper': upper(c_defect).str(45),
        'relative_finite_congruence_defect_upper': upper(finite_defect).str(45),
        'weighted_complement_gram_trace_upper': upper(trace).str(45),
        'potential_sup_upper': upper(potential_bound).str(45),
        'complement_gap_lower': gap.lower().str(45),
        'schur_relative_bound_upper': upper(eta).str(45),
        'minimum_abs_proposed_diagonal': min_d.str(45),
        'full_Dirichlet_inverse_L2_norm_upper': upper(full_inverse).str(45),
    }
    dump(output, result)
    print('DIRICHLET_INVERSE_COMPONENT_PASS', upper(full_inverse).str(12), flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--stage', choices=('prepare', 'basis', 'assemble'), required=True)
    ap.add_argument('--state')
    ap.add_argument('--window', type=int, default=160)
    ap.add_argument('--nodes', type=int, default=1536)
    ap.add_argument('--inputs')
    ap.add_argument('--coverage')
    ap.add_argument('--start', type=int, default=0)
    ap.add_argument('--stop', type=int)
    ap.add_argument('--leaf', action='append', default=[])
    ap.add_argument('--output', required=True)
    args = ap.parse_args()
    if args.stage == 'prepare':
        prepare(args.state, args.window, args.nodes, args.output)
    elif args.stage == 'basis':
        basis_leaf(args.inputs, args.start, args.stop, args.output)
    elif args.stage == 'assemble':
        need(args.coverage is not None, 'missing independent spectral coverage')
        assemble(args.inputs, args.leaf, args.output, args.coverage)


if __name__ == '__main__':
    main()
