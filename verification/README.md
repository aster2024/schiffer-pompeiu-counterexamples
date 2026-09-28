# Ancillary files

Scripts, frozen data and certificates for the computer-assisted parts of the
paper. Every inequality is checked in Arb ball arithmetic (python-flint); an
order relation between balls is true only if it holds for all points of the
balls, every check raises an exception when it fails, and NaN or infinite
enclosures are rejected.

| Path | Content |
|---|---|
| `r3_convex/` | `R³` convex main example: single-file verifier, full replay receipt, final certificate, centre, geometry and scale-obstruction receipts; see `r3_convex/README.md`. |
| `r3/` | additional nonconvex `R³` example: single-file verifier `verify_r3.py`, its certificate, the certificate `certificate_r3_replay2.json` of a second, independent from-scratch replay (identical certified values), and the exact check `check_nonconvex_r3.py` of non-convexity at the poles; see `r3/README.md`. |
| `verify_r4.py`, `verify_convex.py`, `data/centre.json`, `certificates/` | `R⁴ = u(2)`: as described in the paper (commands below). |
| `r6/` | `R⁶ = su(2) ⊕ su(2)`: see `r6/RELEASE_README.md`. |
| `r8/` | `R⁸ = su(3)` (`m = 3`): see `r8/README.md`. |
| `r10/` | `R¹⁰ = so(5)` (`m = 4`): see `r10/README.md`. |
| `r14/` | `R¹⁴ = g₂` (`m = 6`): see `r14/README.md`. |
| `check_radial_rank2.py` | Cross-check of the radial part of the Laplacian for `su(3)`, `so(5)`, `g₂` (see below). |
| `SHA256SUMS` | SHA-256 of every other file under this directory (paths relative to it). |

File integrity, from this directory: `sha256sum -c SHA256SUMS`. The
subdirectories `r3/`, `r3_convex/`, `r6/`, `r8/`, `r10/`, `r14/` also carry their own
`SHA256SUMS`, to be checked from inside the subdirectory.

## Requirements

Python 3 with **python-flint 0.9.0** and NumPy for all dimensions; SymPy for
the identity checks and `check_radial_rank2.py`; for both `R³` verifiers in addition SciPy,
`gcc` and the python-flint binary wheel from PyPI (see `r3/README.md`). All
runs recorded here used Python 3.10.19 on 64-bit Linux, one thread per
process, at `nice -n 19`.

## Commands and running times

| Dimension | Command (from the directory named) | Output | Time |
|---|---|---|---|
| 3 convex | `r3_convex/`: `python3 -O verify_r3_convex.py --workdir <empty dir> --output receipt.json --jobs 5` | `PROVED` | about 1.8 h with 5 workers |
| 3 additional nonconvex | `r3/`: `python3 -I -O verify_r3.py --workdir <new-dir> --jobs 4 --max-seconds 17000 --output cert.json`; `python3 check_nonconvex_r3.py` | `PROVED`, `NONCONVEX_AT_POLES_VERIFIED` | about 67 min with 4 workers; < 1 s |
| 4 | `./`: the two commands below | `PROVED`, `PROVED_CONVEX` | about 165 s |
| 6 | `r6/`: `bash reproduce_r6.sh` | `REPRODUCTION COMPLETE` | about 3 min |
| 8 | `r8/`: `python3 verify_r8.py --stage final --centre center_r8_M45_S24.json --output certificate_r8.json`; `python3 verify_convex_r8.py` | `PROVED`, `PROVED_CONVEX` | < 1 s |
| 10 | `r10/`: `python3 verify_r10.py`; `python3 verify_convex_r10.py` | `PROVED`, `PROVED_CONVEX` | < 1 s |
| 14 | `r14/`: `python3 verify_r14.py`; `python3 verify_convex_r14.py` | `PROVED`, `PROVED_CONVEX` | about 3 s |

For `R⁸`, `R¹⁰`, `R¹⁴` the final gates read frozen interval receipts; the
full recomputation of all receipts from the frozen centre and inverse is
`bash reproduce_r8.sh` (about 70 s), `bash reproduce_r10.sh` (about 105 s)
and `bash reproduce_r14.sh` (about 740 s), each on one core; each ends
with a byte-for-byte comparison of every released JSON file. The scripts
`test_fail_closed_r8.py`, `test_fail_closed_r10.py`, `test_fail_closed_r14.py`
check that the gates reject NaN or infinite values, coverage gaps, and
receipts exceeding a threshold (in scratch copies; 5 to 30 s each).

`R⁴`, from this directory:

```bash
mkdir -p logs
nice -n 19 env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  MKL_NUM_THREADS=1 python3 verify_r4.py --stage all --bits 128 --jstar 1001 \
  --output logs/certificate_r4_128.json
python3 verify_convex.py --bits 128 \
  --main-certificate logs/certificate_r4_128.json \
  --output logs/convex_128.json
```

The files written to `logs/` are byte-identical to those in `certificates/`.

## Rank-two computations (`r8/`, `r10/`, `r14/`)

For a compact Lie algebra `k` of rank two with Weyl group `I₂(m)`, the
adjoint-invariant Schiffer problem reduces to a planar Helmholtz problem on the
Cartan plane with Cauchy data `(π, ∇π)`, `π = Im z^m` (`m = 1, 2` are the
cases `u(2)` and `su(2) ⊕ su(2)` of `R⁴` and `R⁶`). Each directory proves,
for `D = ψ(𝔻)`, the existence of an exact zero of
`H(g,c) = g + |ψ'|²(Kg + Im ψ^m / B)` near a frozen centre, and that `D` is
strictly convex, contains `0` and is not a disc:

| | `k` | `m` | `(M, S, N)` | radius `r` | `Y <` | `Z <` | `1 − U/L >` | first shape coefficient |
|---|---|---:|---|---:|---:|---:|---:|---|
| `R⁸` | `su(3)` | 3 | `(45, 24, 200)` | `1/30000` | `1.51945·10⁻⁵` | `0.418572` | `0.831531` | `|c₇| > 0.0774840` |
| `R¹⁰` | `so(5)` | 4 | `(52, 36, 259)` | `10⁻⁷` | `9.59187·10⁻⁹` | `0.552897` | `0.973053` | `|c₉| > 0.0138324` |
| `R¹⁴` | `g₂` | 6 | `(114, 60, 610)` | `1/6000` | `4.59220·10⁻⁵` | `0.604948` | `0.894067` | `|c₁₃| > 0.0463086` |

`check_radial_rank2.py` checks the radial-part formula
`(Δ_k u)|_t = π⁻¹ Δ_t(π u|_t)` used in the reduction: exactly in explicit
matrix coordinates of `su(3)` (six coordinate identities, harmonicity of `π`,
the drift on the two basic invariants), and numerically at 100 regular points
for the root systems `B₂` and `G₂`. It is a consistency check; the formula is
proved in the paper. Expected output:

```text
SU3_MATRIX_EXACT_PASS 6 coordinate identities; pi harmonic; drift on p2,p3
I2(4)_ROOT_NUMERIC_PASS samples=100 max_relative=2.190e-15
I2(6)_ROOT_NUMERIC_PASS samples=100 max_relative=1.597e-14
```

### Changes relative to the reviewed rank-two code

The rank-two scripts were reviewed before this release. The following changes
were made afterwards; all 106 interval receipts were then recomputed from the
frozen centres and inverses with `reproduce_r*.sh` and are byte-identical to
the reviewed ones.

1. **Exact `κ`.** `verify_r10.py` and `verify_r14.py` evaluated the norm bound
   `κ = 1/((m+2)(m+4))` of `K` as a binary64 float. For `m = 4`, the float
   `1/48` lies about `1.2·10⁻¹⁸` below `1/48`, so the enclosures of `C₂` and
   `C₃` in `R¹⁰` were not strictly rigorous (for `m = 6` the float `1/80`
   rounds up). Both verifiers now use the Arb enclosure `arb(1)/((m+2)(m+4))`.
   The enclosures move in the 16th or 17th significant digit (`R¹⁰`:
   `C₂ = 395.40417671385730…`, `C₃ = 0.071132109967883297…`) and all rational
   thresholds hold as before.
2. **NaN-safe maxima.** The Arb check of `Kg° = (1−r²)²Q` (`check_q` in
   `far_shape_r8.py`, `far_shape_rank2.py`) and the maxima in
   `check_identities_r8.py`, `check_identities_rank2.py` used Python's
   `max()` on Arb balls, which silently skips a NaN that is not the first
   element. They now use `max_upper`, which takes upper endpoints and rejects
   NaN and infinite balls. The index of the largest column recorded in the
   receipts (`arg`, `near_arg`) is now computed by `argmax_upper`
   (`verify_r8.py`, `rank2_cap_core.py`), which is NaN-rejecting as well.
3. **No `assert`.** The `assert` statements in `check_identities_r8.py`,
   `check_identities_rank2.py`, `far_shape_r8.py`, `g_tail_r8.py` and
   `check_radial_rank2.py` were replaced by explicit checks that also run
   under `python -O`.
4. Comments in `rank2_cap_core.py` and `verify_r8.py` that referred to a
   single case were made generic.
5. Because stage scripts changed, the frozen-bundle digests in
   `verify_r8.py`, `verify_r10.py`, `verify_r14.py` and all certificates were
   regenerated. New files: `reproduce_r*.sh`, `test_fail_closed_r*.py`,
   `README.md` and `SHA256SUMS` in each directory.
