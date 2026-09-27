#!/usr/bin/env bash
# Recompute every Arb receipt of the R6 certificate from source, compare it
# byte for byte with the frozen files, then rerun the final gate, the identity
# checks and the convexity certificate and compare those outputs as well.
# Each worker is one core under nice 19 and capped at 600 seconds (the
# largest batch takes well under a minute on one core of a Linux server).
#
# Interpreter: $PYTHON if set, else the first python3 on PATH.
# python-flint 0.9.0 is required.
set -euo pipefail
cd "$(dirname "$0")"
if [ -n "${PYTHON:-}" ]; then
    PY="$PYTHON"
else
    PY="$(command -v python3)"
fi
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
LIMIT=600s
run() { nice -n 19 timeout "$LIMIT" "$PY" "$@"; }

"$PY" - <<'PY'
import sys, flint
print('interpreter:', sys.executable, 'Python', sys.version.split()[0], 'python-flint', flint.__version__)
if flint.__version__ != '0.9.0':
    raise SystemExit('python-flint 0.9.0 is required, found ' + flint.__version__)
PY

run test_fail_closed_r6.py

scratch=$(mktemp -d "${TMPDIR:-/tmp}/r6-recheck-XXXXXX")
trap 'rm -rf "$scratch"' EXIT

run finite_interval_r6.py --output "$scratch/residual_bound_r6.json" >/dev/null
for pair in '0 32' '32 96' '96 160' '160 224' '224 288' '288 352' '352 416' '416 465'; do
    read -r a b <<< "$pair"
    run finite_interval_r6.py --start "$a" --stop "$b" --output "$scratch/finite_${a}_${b}.json" >/dev/null
done

run tail_inverse_r6.py --output "$scratch/tail_inverse_bound_r6.json" >/dev/null
run g_tail_interval_r6.py --output "$scratch/g_far_bound_r6.json" >/dev/null
for pair in '0 100' '100 200' '200 300' '300 400' '400 500' '500 600' '600 705'; do
    read -r a b <<< "$pair"
    run g_tail_interval_r6.py --start "$a" --stop "$b" --output "$scratch/g_${a}_${b}.json" >/dev/null
done

for pair in '61 61' '65 125' '129 229' '233 297' '301 397'; do
    read -r a b <<< "$pair"
    run shape_tail_interval_r6.py --first "$a" --last "$b" --output "$scratch/shape_${a}_${b}.json" >/dev/null
done
run far_shape_r6.py --jstar 401 --output "$scratch/far_shape_bound_r6.json" >/dev/null

"$PY" - "$scratch" <<'PY'
import hashlib, json, pathlib, sys
p = pathlib.Path(sys.argv[1])
groups = {
    'finite_batches_r6.json': ['finite_0_32.json','finite_32_96.json','finite_96_160.json',
        'finite_160_224.json','finite_224_288.json','finite_288_352.json',
        'finite_352_416.json','finite_416_465.json'],
    'g_tail_batches_r6.json': ['g_0_100.json','g_100_200.json','g_200_300.json',
        'g_300_400.json','g_400_500.json','g_500_600.json','g_600_705.json'],
    'shape_tail_batches_r6.json': ['shape_61_61.json','shape_65_125.json',
        'shape_129_229.json','shape_233_297.json','shape_301_397.json'],
}
for name, files in groups.items():
    rows = [json.loads((p/f).read_text()) for f in files]
    (p/name).write_text(json.dumps(rows,indent=2))
names = ['residual_bound_r6.json','tail_inverse_bound_r6.json','g_far_bound_r6.json',
         'far_shape_bound_r6.json',*groups]
digest = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
for name in names:
    if digest(p/name) != digest(pathlib.Path(name)):
        raise SystemExit(f'REPRODUCTION MISMATCH: {name}')
print('All seven frozen interval receipts reproduced byte for byte.')
PY

# Final gate on the frozen receipts; its certificate must match the shipped one.
run verify_r6.py --stage final --output "$scratch/certificate_r6.json"
cmp -s "$scratch/certificate_r6.json" certificate_r6.json \
    || { echo 'REPRODUCTION MISMATCH: certificate_r6.json'; exit 1; }
echo 'certificate_r6.json reproduced byte for byte.'

# Exact/Arb identity checks: K g = (1-r^2)^2 Q and the principal column (5.1).
run check_identities_r6.py --output "$scratch/identities_r6.json"
cmp -s "$scratch/identities_r6.json" identities_r6.json \
    || { echo 'REPRODUCTION MISMATCH: identities_r6.json'; exit 1; }
echo 'identities_r6.json reproduced byte for byte.'

# Convexity of the meridian domain for every shape in the certified ball.
run verify_convex_r6.py --main-certificate "$scratch/certificate_r6.json" \
    --output "$scratch/convex_r6.json"
cmp -s "$scratch/convex_r6.json" convex_r6.json \
    || { echo 'REPRODUCTION MISMATCH: convex_r6.json'; exit 1; }
echo 'convex_r6.json reproduced byte for byte.'
echo 'REPRODUCTION COMPLETE'
