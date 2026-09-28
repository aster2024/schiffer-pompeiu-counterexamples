#!/usr/bin/env bash
# Full recomputation of the R^8 (su(3)) interval bundle.
#
# Copies the scripts, the frozen dyadic centre and the frozen binary64 inverse
# into a fresh scratch directory, recomputes every interval receipt there,
# reruns the final gate (verify_r8.py) and the convexity gate
# (verify_convex_r8.py), and compares every released JSON file byte for byte
# with its recomputed counterpart.  The released directory is not modified.
#
# Environment: PYTHON (default python3), CAP_TIMEOUT per job in seconds
# (default 900), TMPDIR.  One thread, nice 19.
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
PY=${PYTHON:-python3}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
work=$(mktemp -d "${TMPDIR:-/tmp}/r8-recheck-XXXXXX")
cp "$here"/*.py "$here"/center_r8_M45_S24.json "$here"/inverse_r8_M45_S24.npy "$work"/
cd "$work"
run_cap() { nice -n 19 timeout "${CAP_TIMEOUT:-900}" "$PY" "$@" >/dev/null; }

run_cap finite_interval_r8.py
run_cap tail_inverse_r8.py
run_cap g_tail_r8.py --output g_far_bound_r8.json
for ((batch_start=0; batch_start<396; batch_start+=20)); do
  batch_stop=$((batch_start+20))
  if ((batch_stop>396)); then batch_stop=396; fi
  printf -v receipt 'g_near_r8_%03d_%03d.json' "$batch_start" "$batch_stop"
  run_cap g_tail_r8.py --start "$batch_start" --stop "$batch_stop" --output "$receipt"
done
run_cap shape_tail_r8.py --first 49 --last 133 --output shape_near_r8_049_133.json
run_cap shape_tail_r8.py --first 139 --last 217 --output shape_near_r8_139_217.json
run_cap far_shape_r8.py --jstar 223 --output far_shape_bound_r8.json
run_cap check_identities_r8.py
nice -n 19 "$PY" verify_r8.py --stage final --centre center_r8_M45_S24.json --output certificate_r8.json
nice -n 19 "$PY" verify_convex_r8.py

for f in "$here"/*.json; do
  cmp "$f" "$work/$(basename "$f")"
done
echo "R8: all $(ls "$here"/*.json | wc -l) released JSON files reproduced byte for byte in $work"
