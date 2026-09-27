#!/usr/bin/env python3
"""Negative tests: verify_r6.max_upper and receipt parsing must reject NaN/inf.

The pre-release max_upper was max(v.upper() for v in seq), which silently
drops a NaN that is not the first element (e.g. max_upper([1, nan]) == 1).
"""
from flint import arb

import verify_r6 as V

CASES = {
    'nan after finite': [arb(1), arb('nan')],
    'nan first': [arb('nan'), arb(1)],
    'infinite value': [arb(1), arb('inf')],
    'infinite radius': [arb(1), arb('[0 +/- inf]')],
}

for label, seq in CASES.items():
    try:
        V.max_upper(seq)
    except RuntimeError:
        continue
    raise SystemExit('NOT FAIL-CLOSED: max_upper accepted ' + label)

for text in ('nan', 'inf', '[1 +/- inf]'):
    try:
        V._interval(text)
    except RuntimeError:
        continue
    raise SystemExit('NOT FAIL-CLOSED: receipt value accepted ' + text)

assert V.max_upper([arb(1) / 3, arb(2) / 3]) >= arb(2) / 3
assert V.max_upper([]) == 0
print('FAIL-CLOSED TESTS PASSED')
