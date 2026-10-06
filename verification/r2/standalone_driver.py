#!/usr/bin/env python3
"""Driver template. The released verifier embeds all sources and exact inputs."""
import argparse
import base64
from concurrent.futures import ThreadPoolExecutor, as_completed
import copy
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import time
import zlib

# REPLACE_EMBEDDED_BUNDLE


def need(ok, reason):
    if not ok:
        raise RuntimeError(reason)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load(path):
    def bad(value):
        raise ValueError('nonfinite JSON constant: '+value)
    return json.loads(Path(path).read_text(), parse_constant=bad)


def run_stage(work, label, argv, timeout=7200):
    env = os.environ.copy()
    env.update(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1',
               NUMEXPR_NUM_THREADS='1', PYTHONDONTWRITEBYTECODE='1')
    start = time.monotonic()
    with (work/(label+'.log')).open('w') as log:
        proc = subprocess.run(['nice', '-n', '19', sys.executable, '-O']+argv,
                              cwd=work, env=env, stdout=log, stderr=subprocess.STDOUT,
                              timeout=timeout)
    need(proc.returncode == 0, 'stage failed: '+label+'; see '+str(work/(label+'.log')))
    return {'stage': label, 'elapsed_seconds': time.monotonic()-start,
            'returncode': proc.returncode}



def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--workdir', required=True)
    ap.add_argument('--output')
    ap.add_argument('--jobs', type=int, default=5)
    ap.add_argument('--stage-timeout', type=int, default=7200)
    args = ap.parse_args()
    need(1 <= args.jobs <= 5, 'jobs must be between one and five')
    need(args.stage_timeout > 0, 'nonpositive stage timeout')
    work = Path(args.workdir).resolve()
    need(not work.exists() or (work.is_dir() and not any(work.iterdir())),
         'replay requires a new or empty directory; existing receipts are forbidden')
    out = Path(args.output).resolve() if args.output else work/'VERIFICATION_RECEIPT.json'
    need(not out.exists(), 'refusing to reuse a pre-existing final receipt')
    work.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    bundle_bytes = zlib.decompress(base64.b64decode(EMBEDDED_BUNDLE))
    need(hashlib.sha256(bundle_bytes).hexdigest() == EMBEDDED_SHA256, 'embedded bundle corrupt')
    bundle = json.loads(bundle_bytes)
    for name, content in bundle.items():
        need(Path(name).name == name and name.endswith(('.py', '.json')), 'unsafe embedded path')
        (work/name).write_text(content)
    sys.path.insert(0, str(work))
    from final_gates import gate, RADIUS
    records = []
    def run(label, argv):
        rec = run_stage(work, label, argv, args.stage_timeout)
        print('REPLAY_STAGE_PASS', label, flush=True)
        return rec
    records.append(run('independent_coverage', ['coverage_checks.py', '--reference',
        'reference.json', '--output', 'independent_coverage.json']))
    records.append(run('spectral_prepare', ['spectral_interval.py', '--stage', 'prepare',
        '--state', 'reference.json', '--window', '160', '--nodes', '1536', '--output', 'spectral_inputs.json']))
    n = len(load(work/'spectral_inputs.json')['states'])
    need(n == 408, 'unexpected spectral window dimension')
    basis_names = []
    units = []
    for a in range(0, n, 24):
        b = min(a+24, n)
        name = f'spectral_basis_{a:04d}_{b:04d}.json'
        basis_names.append(name)
        units.append((name[:-5], ['spectral_interval.py', '--stage', 'basis', '--inputs',
            'spectral_inputs.json', '--start', str(a), '--stop', str(b), '--output', name]))
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = [pool.submit(run, *unit) for unit in units]
        for future in as_completed(futures):
            records.append(future.result())
    command = ['spectral_interval.py', '--stage', 'assemble', '--inputs',
               'spectral_inputs.json', '--coverage', 'independent_coverage.json', '--output', 'spectral_inverse.json']
    for name in basis_names:
        command += ['--leaf', name]
    records.append(run('spectral_assemble', command))
    boundary_names = []
    units = []
    for a in range(0, 1024, 128):
        b = a+128
        name = f'boundary_{a:04d}_{b:04d}.json'
        boundary_names.append(name)
        units.append((name[:-5], ['boundary_interval.py', '--stage', 'leaf', '--state',
            'candidate.json', '--grid', '1024', '--lmax', '400', '--start', str(a),
            '--stop', str(b), '--output', name]))
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = [pool.submit(run, *unit) for unit in units]
        for future in as_completed(futures):
            records.append(future.result())
    command = ['boundary_interval.py', '--stage', 'merge', '--state',
               'candidate.json', '--output', 'clamped_residual.json']
    for name in boundary_names:
        command += ['--leaf', name]
    records.append(run('boundary_merge', command))
    records.append(run('coupled_radii', ['analytic_cap.py', '--state', 'candidate.json',
        '--reference', 'reference.json', '--spectral', 'spectral_inverse.json',
        '--boundary', 'clamped_residual.json', '--radius', RADIUS, '--output', 'coupled_radii.json']))
    records.append(run('identities', ['identity_checks.py']))
    spectral, boundary, cap = (load(work/p) for p in
                               ('spectral_inverse.json', 'clamped_residual.json', 'coupled_radii.json'))
    candidate, reference = load(work/'candidate.json'), load(work/'reference.json')
    scan = load(work/'independent_coverage.json')
    need(spectral['independent_coverage_sha256'] == digest(work/'independent_coverage.json'),
         'independent coverage file changed')
    final = gate(spectral, boundary, cap, candidate, reference, scan)
    records.append(run('second_interval_scalar', ['scalar_second_iv.py', '--candidate', 'candidate.json',
        '--reference', 'reference.json', '--spectral', 'spectral_inverse.json',
        '--boundary', 'clamped_residual.json', '--radius', RADIUS, '--output', 'scalar_second_iv.json']))
    from negative_controls import negative_tests
    safety = negative_tests(spectral, boundary, cap, candidate, reference, scan, work)
    (work/'fail_closed_tests.json').write_text(json.dumps(safety, indent=2)+'\n')
    print(safety['status'], len(safety['rejected']), flush=True)
    import flint, numpy, scipy, sympy, mpmath
    artifacts = {p.name: digest(p) for p in sorted(work.iterdir())
                 if p.is_file() and p.suffix in ('.json', '.py', '.log')}
    receipt = {'status': 'NUMERICAL_GATES_PASS', 'full_proof': False,
               'analytic_proof': 'PROOF_PACKAGE.md; analytic arguments are not machine-checked',
               'trust_base': 'CPython, FLINT/Arb/python-flint; second scalar chain in mpmath.iv',
               'claim': 'bounded real-analytic strictly convex non-disc D182 planar Schiffer domain',
               'proof_scope': 'local numerical hypotheses of the written proof; no new AI or external review',
               'verifier_sha256': digest(__file__), 'embedded_bundle_sha256': EMBEDDED_SHA256,
               'candidate_sha256': digest(work/'candidate.json'),
               'reference_sha256': digest(work/'reference.json'),
               'radius': RADIUS, 'norm': cap['norm'], 'bounds': cap['bounds'],
               'final_rechecked_bounds': final,
               'Dirichlet_inverse_L2_upper': spectral['full_Dirichlet_inverse_L2_norm_upper'],
               'residual_E0_upper': boundary['E0_upper'], 'residual_E1_upper': boundary['E1_upper'],
               'coverage_counts': scan['counts'],
               'second_interval_scalar': load(work/'scalar_second_iv.json'),
               'stages': records, 'stage_timeout_seconds': args.stage_timeout, 'jobs': args.jobs, 'wall_seconds': time.monotonic()-started,
               'versions': {'python': platform.python_version(), 'python_flint': flint.__version__,
                            'numpy': numpy.__version__, 'scipy': scipy.__version__, 'sympy': sympy.__version__, 'mpmath': mpmath.__version__},
               'fail_closed_tests': safety, 'artifacts_sha256': artifacts}
    out.write_text(json.dumps(receipt, indent=2)+'\n')
    readback = load(out)
    need(readback['status'] == 'NUMERICAL_GATES_PASS' and readback['full_proof'] is False
         and readback['verifier_sha256'] == digest(__file__), 'receipt readback failed')
    for name, expected in readback['artifacts_sha256'].items():
        need(digest(work/name) == expected, 'artifact readback changed: '+name)
    print('NUMERICAL_GATES_PASS', flush=True)
    print('receipt:', out, flush=True)


if __name__ == '__main__':
    try:
        main()
    except Exception as error:
        print('NUMERICAL_GATES_FAIL:', type(error).__name__+':', str(error), file=sys.stderr, flush=True)
        raise SystemExit(1)
