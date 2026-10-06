#!/usr/bin/env python3
"""Two complementary spectral-window completeness certificates.

The sign scan uses only exact reference decimals and J_n evaluations. It never
uses scipy, root proposals, the bracket list, or Newton's method.
"""
import argparse
import json
from pathlib import Path
from flint import arb, ctx
from spectral_interval import need, upper, exact_decimal, shape_digest, dump
ctx.threads = 1

EXPECTED = (102, 99, 91, 74, 37, 5, 0)


def check_counts(counts, include_high=False):
    target = EXPECTED if include_high else EXPECTED[:-1]
    need(isinstance(counts, list) and len(counts) == len(target), 'sector count dimensions')
    for j, count in enumerate(target):
        need(type(counts[j]) is int and counts[j] == count, f'sector count j={j}')


def check_root_coverage(data):
    """Recheck guards and require exactly the full in-window root list at assembly."""
    ctx.prec = 2048
    need(data['stage'] == 'CERTIFIED_SPECTRAL_INPUTS', 'coverage input stage')
    m = data['m']
    need(m == data['state']['m'] == 182 and data['window'] == 160, 'coverage window contract')
    rho = exact_decimal(data['state']['p'][0])
    lo, hi = rho-160, rho+160
    need(lo > 1 and 6*m > hi, 'coverage angular cutoff')
    coverage = data['coverage']
    need(len(coverage) == 6, 'missing coverage sector')
    all_states, counts = [], []
    for j, sector in enumerate(coverage):
        n = j*m
        need(sector['j'] == j and sector['n'] == n, 'coverage sector identity')
        roots = sector['roots']
        need(len(roots) >= 2, 'missing root guards')
        count = 0
        for root in roots:
            need(root['n'] == n, 'root angular order')
            left, right, norm = (arb(root[k]) for k in ('left', 'right', 'norm'))
            for v in (left, right, norm):
                upper(v)
            need(1 < left < right and right-left < arb(1)/10, 'coverage bracket geometry')
            need(left.bessel_j(n)*right.bessel_j(n) < 0, 'coverage root endpoint signs')
            # Real-axis Lipschitz bound |J_n prime| <= 1 encloses the normalization.
            mid = (left+right)/2
            true_norm = abs(arb(mid.mid().bessel_j(n+1),
                                upper(mid.rad()+(right-left)/2)))
            need(norm > 0 and norm.contains(true_norm), 'coverage normalization enclosure')
            if lo < left and right < hi:
                all_states.append(dict(root, j=j))
                count += 1
            else:
                need(right < lo or hi < left, 'coverage cutoff straddle')
        first = sector['starts_at_first_root']
        need(type(first) is bool, 'first-root flag type')
        if first:
            need(n > 0, 'order zero needs lower guard')
            end = arb(roots[0]['right'])
            qmax = 1-(arb(n*n)-arb(1)/4)/(end*end)
            need(qmax > 0 and end-n < 3*arb.pi()/(2*qmax.sqrt()), 'first-root coverage')
        else:
            need(arb(roots[0]['right']) < lo, 'coverage lower guard')
        for r1, r2 in zip(roots, roots[1:]):
            left, right1 = arb(r1['left']), arb(r1['right'])
            left2, right = arb(r2['left']), arb(r2['right'])
            endpoint = left if n == 0 else right
            qmax = 1-(arb(n*n)-arb(1)/4)/(endpoint*endpoint)
            need(qmax > 0 and right1 < left2 and right-left < 2*arb.pi()/qmax.sqrt(),
                 'coverage skipped or duplicate root')
        need(hi < arb(roots[-1]['left']), 'coverage upper guard')
        counts.append(count)
    check_counts(counts)
    need(data['states'] == all_states, 'states differ from all covered in-window roots')
    if 'sector_counts' in data:
        need(data['sector_counts'] == counts, 'stored coverage counts differ')
    return counts


def certified_sign(n, rho_text, k):
    for precision in (256, 1024, 4096):
        ctx.prec = precision
        # Reconstruct the exact point at each precision; no rounded grid proposal.
        x = exact_decimal(rho_text)-160+arb(k)/4
        value = x.bessel_j(n)
        upper(value)
        if value > 0:
            return 1
        if value < 0:
            return -1
    raise RuntimeError(f'unresolved independent sign n={n} k={k}')


def sign_scan(reference):
    ctx.prec = 256
    need(reference['m'] == 182, 'sign scan symmetry')
    rho = exact_decimal(reference['p'][0])
    spacing = arb.pi()/(arb(5)/4).sqrt()
    need(rho-160 > 1 and spacing > arb(28)/10 and spacing > arb(1)/4,
         'sign scan zero spacing')
    need(1092 > rho+160, 'sign scan higher-order exclusion')
    sectors, counts = [], []
    for j in range(7):
        signs = [certified_sign(182*j, reference['p'][0], k) for k in range(1281)]
        changes = [k for k in range(1280) if signs[k] != signs[k+1]]
        counts.append(len(changes))
        sectors.append({'j': j, 'n': 182*j, 'signs': signs, 'change_cells': changes})
        print('INDEPENDENT_ZERO_COUNT', 182*j, len(changes), flush=True)
    check_counts(counts, include_high=True)
    return {'stage': 'INDEPENDENT_QUARTER_GRID_COVERAGE',
            'method': 'certified_Jn_signs_and_zero_spacing_gt_2.8',
            'reference_shape_sha256': shape_digest(reference),
            'm': 182, 'window': 160, 'grid_cells': 1280, 'step': '0.25',
            'counts': counts, 'total': sum(counts), 'sectors': sectors}


def check_scan_record(scan, reference):
    need(scan['stage'] == 'INDEPENDENT_QUARTER_GRID_COVERAGE'
         and scan['reference_shape_sha256'] == shape_digest(reference), 'independent coverage identity')
    need(scan['m'] == 182 and scan['window'] == 160 and scan['grid_cells'] == 1280
         and scan['step'] == '0.25', 'independent coverage grid contract')
    check_counts(scan['counts'], include_high=True)
    need(scan['total'] == 408 and len(scan['sectors']) == 7, 'independent coverage total/sectors')
    for j, sector in enumerate(scan['sectors']):
        need(sector['j'] == j and sector['n'] == 182*j, 'independent sector identity')
        signs = sector['signs']
        need(len(signs) == 1281 and all(type(s) is int and s in (-1, 1) for s in signs),
             'independent nonzero signs missing')
        cells = [k for k in range(1280) if signs[k] != signs[k+1]]
        need(cells == sector['change_cells'] and len(cells) == scan['counts'][j],
             'independent sign-change count')


def check_scan_binding(scan, data):
    check_scan_record(scan, data['state'])
    ctx.prec = 2048
    lo = exact_decimal(data['state']['p'][0])-160
    # Match each certified root bracket to a distinct sign-changing grid cell.
    for j in range(6):
        roots = [s for s in data['states'] if s['j'] == j]
        cells = scan['sectors'][j]['change_cells']
        need(len(roots) == len(cells), 'independent root count mismatch')
        for root, k in zip(roots, cells):
            need(lo+arb(k)/4 < arb(root['left']) and
                 arb(root['right']) < lo+arb(k+1)/4, 'independent root-cell mismatch')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--reference', required=True)
    ap.add_argument('--output', required=True)
    args = ap.parse_args()
    ref = json.loads(Path(args.reference).read_text())
    result = sign_scan(ref)
    check_scan_record(result, ref)
    dump(args.output, result)


if __name__ == '__main__':
    main()
