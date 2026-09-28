#!/usr/bin/env bash
# Full recomputation of the R^10 (so(5)) interval bundle.
#
# Copies the scripts, the frozen dyadic centre and the frozen binary64 inverse
# into a fresh scratch directory, recomputes every interval receipt there,
# reruns the final gate (verify_r10.py) and the convexity gate
# (verify_convex_r10.py), and compares every released JSON file byte for byte
# with its recomputed counterpart.  The released directory is not modified.
#
# Environment: PYTHON (default python3), CAP_TIMEOUT per job in seconds
# (default 900), TMPDIR.  One thread, nice 19.
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
PY=${PYTHON:-python3}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
work=$(mktemp -d "${TMPDIR:-/tmp}/r10-recheck-XXXXXX")
cp "$here"/*.py "$here"/center_r10_M52_S36.json "$here"/inverse_r10_M52_S36.npy "$work"/
cd "$work"
run_cap() { nice -n 19 timeout "${CAP_TIMEOUT:-900}" "$PY" "$@" >/dev/null; }
C=(--centre center_r10_M52_S36.json)
I=(--inverse inverse_r10_M52_S36.npy)
R=(--finite-receipt finite_bound_r10.json --tail-receipt tail_inverse_bound_r10.json)

run_cap finite_rank2.py "${C[@]}" "${I[@]}" --output finite_bound_r10.json
run_cap tail_inverse_rank2.py "${C[@]}" "${I[@]}" --output tail_inverse_bound_r10.json
run_cap g_tail_rank2.py "${C[@]}" "${I[@]}" "${R[@]}" --short-degree 24 --output g_far_bound_r10.json
for ((batch_start=0; batch_start<358; batch_start+=20)); do
  batch_stop=$((batch_start+20))
  if ((batch_stop>358)); then batch_stop=358; fi
  printf -v receipt 'g_near_r10_%03d_%03d.json' "$batch_start" "$batch_stop"
  run_cap g_tail_rank2.py "${C[@]}" "${I[@]}" "${R[@]}" --short-degree 24 \
    --start "$batch_start" --stop "$batch_stop" --output "$receipt"
done
for pair in '57 169' '177 217' '225 265' '273 297'; do
  read -r first last <<< "$pair"
  printf -v receipt 'shape_near_r10_%03d_%03d.json' "$first" "$last"
  run_cap shape_tail_rank2.py "${C[@]}" "${I[@]}" --first "$first" --last "$last" --output "$receipt"
done
run_cap far_shape_rank2.py "${C[@]}" --tail-receipt tail_inverse_bound_r10.json --jstar 305 --output far_shape_bound_r10.json
run_cap check_identities_rank2.py "${C[@]}" --output identities_r10.json
nice -n 19 "$PY" verify_r10.py
nice -n 19 "$PY" verify_convex_r10.py

for f in "$here"/*.json; do
  cmp "$f" "$work/$(basename "$f")"
done
echo "R10: all $(ls "$here"/*.json | wc -l) released JSON files reproduced byte for byte in $work"
