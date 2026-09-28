#!/usr/bin/env python3
"""Negative tests for the fail-closed behaviour of the rank-two final gates.

The same file is shipped as test_fail_closed_r8.py, test_fail_closed_r10.py
and test_fail_closed_r14.py; it tests the gate of the directory it lives in.

1. The maximum helpers max_upper and argmax_upper must raise on a NaN or
   infinite ball anywhere in the sequence (Python's max() would silently skip
   a NaN that is not the first element).
2. In a scratch copy of this directory, the final gate must print PROVED on
   the released files, and must refuse a receipt that was changed by a single
   field, both with the frozen-bundle checksum in force and with that
   checksum bypassed (so that the NaN, coverage and threshold checks behind
   it are exercised as well).

The released files are never modified. Prints FAIL-CLOSED TESTS PASSED.
"""
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from flint import arb

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

if (HERE / 'verify_r8.py').exists():
    DIM = 8
    import verify_r8 as LIB
    GATE = ("import sys, verify_r8 as V;{bypass}"
            "sys.argv=['verify_r8.py','--stage','final','--centre',"
            "'center_r8_M45_S24.json','--output','certificate_r8.json'];V.main()")
    BYPASS = "V.frozen_digest=lambda: V.FROZEN_BUNDLE_SHA256;"
    CASES = [
        ('g_near_r8_100_120.json', 'max', 'nan'),
        ('g_near_r8_100_120.json', 'max', '[0.1 +/- inf]'),
        ('g_near_r8_100_120.json', 'start', 101),
        ('far_shape_bound_r8.json', 'Z_shape_far', 'nan'),
        ('finite_bound_r8.json', 'Y', '[1e-5 +/- nan]'),
        ('tail_inverse_bound_r8.json', 'A_tail_upper', '2.31'),
        ('shape_near_r8_139_217.json', 'max', '0.44'),
        ('identities_r8.json', 'Q_connection_max_error', 'nan'),
    ]
elif (HERE / 'verify_r10.py').exists():
    DIM = 10
    import rank2_cap_core as LIB
    GATE = "import verify_r10 as V;{bypass}V.main()"
    BYPASS = "V.digest=lambda: V.FROZEN_BUNDLE_SHA256;"
    CASES = [
        ('g_near_r10_100_120.json', 'max', 'nan'),
        ('g_near_r10_100_120.json', 'stop', 119),
        ('far_shape_bound_r10.json', 'Z_shape_far', 'inf'),
        ('tail_inverse_bound_r10.json', 'A_tail_upper', 'nan'),
        ('shape_near_r10_177_217.json', 'max', '0.99'),
        ('finite_bound_r10.json', 'inverse_defect', '[0.5 +/- nan]'),
        ('identities_r10.json', 'ball_principal_max_error', 'nan'),
    ]
elif (HERE / 'verify_r14.py').exists():
    DIM = 14
    import rank2_cap_core as LIB
    GATE = "import verify_r14 as V;{bypass}V.main()"
    BYPASS = "V.digest=lambda: V.FROZEN_BUNDLE_SHA256;"
    CASES = [
        ('finite_r14_300_340.json', 'Z_finite_max', 'nan'),
        ('finite_r14_300_340.json', 'start', 301),
        ('g_near_r14_0520_0620.json', 'max', '[0.2 +/- nan]'),
        ('far_shape_split_r14_j241.json', 'Z_shape_far', 'nan'),
        ('g_far_bound_r14_M114_S60.json', 'far_n', '-inf'),
        ('shape_near_r14_181_181.json', 'max', '0.62'),
        ('identities_r14_M114_S60.json', 'Q_connection_max_error', 'nan'),
    ]
else:
    raise SystemExit('run this test from an r8/, r10/ or r14/ directory')


def expect_raise(label, func, *args):
    try:
        func(*args)
    except RuntimeError:
        return
    raise SystemExit('NOT FAIL-CLOSED: ' + label)


# 1. Maximum helpers.
BAD = {
    'nan after finite': [arb(1), arb('nan')],
    'nan first': [arb('nan'), arb(1)],
    'infinite value': [arb(1), arb('inf')],
    'infinite radius': [arb(1), arb('[0 +/- inf]')],
}
for label, seq in BAD.items():
    expect_raise('max_upper accepted ' + label, LIB.max_upper, seq)
    expect_raise('argmax_upper accepted ' + label, LIB.argmax_upper, seq)
if not (LIB.max_upper([arb(1) / 3, arb(2) / 3]) >= arb(2) / 3
        and LIB.argmax_upper([arb(1) / 3, arb(2) / 3, arb(2) / 3]) == 1
        and LIB.max_upper([]) == 0):
    raise SystemExit('maximum helpers return wrong values')

# 2. Final gate on released and on tampered receipts.
ENV = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1',
           MKL_NUM_THREADS='1', PYTHONDONTWRITEBYTECODE='1')


def gate(workdir, bypass):
    code = GATE.format(bypass=BYPASS if bypass else '')
    run = subprocess.run([sys.executable, '-c', code], cwd=workdir, env=ENV,
                         capture_output=True, text=True)
    proved = run.returncode == 0 and 'PROVED' in run.stdout.split()
    reason = (run.stderr.strip().splitlines() or [''])[-1]
    return proved, reason


def scratch_copy():
    work = Path(tempfile.mkdtemp(prefix=f'r{DIM}-failclosed-'))
    for f in HERE.iterdir():
        if f.is_file() and f.suffix in ('.py', '.json', '.npy'):
            shutil.copy2(f, work / f.name)
    return work


def tamper(work, name, key, value):
    path = work / name
    data = json.loads(path.read_text())
    if key not in data:
        raise SystemExit(f'test case refers to a missing field {name}:{key}')
    data[key] = value
    path.write_text(json.dumps(data, indent=2) + '\n')


work = scratch_copy()
try:
    for bypass in (False, True):
        proved, reason = gate(work, bypass)
        if not proved:
            raise SystemExit('released files not accepted: ' + reason)
    print(f'R{DIM}: released receipts accepted (with and without checksum bypass)')
    for name, key, value in CASES:
        for bypass in (False, True):
            shutil.copy2(HERE / name, work / name)
            tamper(work, name, key, value)
            proved, reason = gate(work, bypass)
            if proved:
                raise SystemExit(f'NOT FAIL-CLOSED: gate accepted {name}:{key}={value!r}')
            print(f'R{DIM}: {name}:{key}={value!r} rejected'
                  f'{" (checksum bypassed)" if bypass else ""}: {reason}')
        shutil.copy2(HERE / name, work / name)
finally:
    shutil.rmtree(work)
print('FAIL-CLOSED TESTS PASSED')
