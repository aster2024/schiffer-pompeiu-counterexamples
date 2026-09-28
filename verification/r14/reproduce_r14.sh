#!/usr/bin/env bash
# Full recomputation of the R^14 (g_2) interval bundle.
#
# Copies the scripts, the frozen dyadic centre and the frozen binary64 inverse
# into a fresh scratch directory, recomputes every interval receipt there
# (610 finite columns in 26 batches, 1215 near g-tail columns in 13 batches,
# 8 near shape batches, the split far-shape bound), reruns the final gate
# (verify_r14.py) and the convexity gate (verify_convex_r14.py), and compares
# every released JSON file byte for byte with its recomputed counterpart.
# The released directory is not modified.
#
# Environment: PYTHON (default python3), CAP_TIMEOUT per job in seconds
# (default 900), TMPDIR.  One thread, nice 19.
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
PY=${PYTHON:-python3}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
work=$(mktemp -d "${TMPDIR:-/tmp}/r14-recheck-XXXXXX")
cp "$here"/*.py "$here"/center_r14_M114_S60.json "$here"/inverse_r14_M114_S60.npy "$work"/
cd "$work"
run_cap() { nice -n 19 timeout "${CAP_TIMEOUT:-900}" "$PY" "$@" >/dev/null; }
C=(--centre center_r14_M114_S60.json)
I=(--inverse inverse_r14_M114_S60.npy)
R=(--finite-receipt finite_bound_r14_M114_S60.json --tail-receipt tail_inverse_bound_r14_M114_S60.json)

# Finite Jacobian columns 0..609 in the 26 nonoverlapping batches the verifier expects.
batches=('0 20')
for a in 20 60 100 140 180 220 260 300 340 380 420 460 500 540; do batches+=("$a $((a+40))"); done
batches+=('580 600')
for ((a=600; a<610; a++)); do batches+=("$a $((a+1))"); done
for pair in "${batches[@]}"; do
  read -r a b <<< "$pair"
  printf -v receipt 'finite_r14_%03d_%03d.json' "$a" "$b"
  run_cap finite_batch_rank2.py "${C[@]}" "${I[@]}" --start "$a" --stop "$b" --output "$receipt"
done
run_cap aggregate_finite_r14.py
run_cap tail_inverse_rank2.py "${C[@]}" "${I[@]}" --output tail_inverse_bound_r14_M114_S60.json
run_cap g_tail_rank2.py "${C[@]}" "${I[@]}" "${R[@]}" --short-degree 60 --output g_far_bound_r14_M114_S60.json
gbatches=('0 20')
for a in 20 120 220 320 420 520 620 720 820 920 1020; do gbatches+=("$a $((a+100))"); done
gbatches+=('1120 1215')
for pair in "${gbatches[@]}"; do
  read -r a b <<< "$pair"
  printf -v receipt 'g_near_r14_%04d_%04d.json' "$a" "$b"
  run_cap g_tail_rank2.py "${C[@]}" "${I[@]}" "${R[@]}" --short-degree 60 \
    --start "$a" --stop "$b" --output "$receipt"
done
for pair in '121 121' '133 145' '157 169' '181 181' '193 193' '205 205' '217 217' '229 229'; do
  read -r first last <<< "$pair"
  printf -v receipt 'shape_near_r14_%03d_%03d.json' "$first" "$last"
  run_cap shape_tail_rank2.py "${C[@]}" "${I[@]}" --first "$first" --last "$last" --output "$receipt"
done
run_cap far_shape_split_r14.py --jstar 241 --output far_shape_split_r14_j241.json
run_cap check_identities_rank2.py "${C[@]}" --output identities_r14_M114_S60.json
nice -n 19 "$PY" verify_r14.py
nice -n 19 "$PY" verify_convex_r14.py

for f in "$here"/*.json; do
  cmp "$f" "$work/$(basename "$f")"
done
echo "R14: all $(ls "$here"/*.json | wc -l) released JSON files reproduced byte for byte in $work"
