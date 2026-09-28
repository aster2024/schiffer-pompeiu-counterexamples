# R¹⁰ = so(5): computer-assisted proof package

Computer-assisted part of the proof that there is a bounded **convex** domain
`Ω ⊂ so(5) ≅ sp(2) ≅ R¹⁰`, invariant under the adjoint action of `SO(5)`,
with real-analytic boundary diffeomorphic to `S⁹`, not a ball, carrying a
nonconstant solution of

```text
Δu + u = 0 in Ω,    u = 1,  ∇u = 0 on ∂Ω .
```

Via the radial part of the Laplacian on the Cartan plane (Weyl group `I₂(4)`,
product of positive roots `π = Im z⁴` after normalisation) this reduces to a
planar Helmholtz problem on `D = ψ(𝔻)`, `ψ(w) = Σ_{j ≡ 1 (mod 8)} c_j w^j`,
with Cauchy data `(π, ∇π)` on `∂D`, written as
`H(g,c) = g + |ψ'|²(Kg + Im ψ⁴ / B) = 0` in the disk-polynomial sine sector
`n ≡ 4 (mod 8)`. Here `B > 0` is a fixed field scale (an exact dyadic stored
in the centre file, `B ≈ 2.1204473·10⁶`); the physical field is
`F = B·(Kg + Im ψ⁴/B)∘ψ⁻¹`, so `B` does not change the Cauchy data. The
scripts below prove the existence of an exact zero of `H` near a frozen centre
(contraction argument in a weighted ℓ¹ space with weight `ρ = 11/10`) and
certify that `D` is strictly convex, contains `0`, and is not a disc. The
Lie-theoretic reduction, the lifting to `R¹⁰` and the Pompeiu corollary are
proved in the paper. Rank-two notes common to `r8/`, `r10/`, `r14/`,
including the changes made to the reviewed code, are in `../README.md`.

## Requirements

Python 3 (runs used Python 3.10.19) with **python-flint 0.9.0** and NumPy;
SymPy for `check_identities_rank2.py`. One thread. All commands are run from
this directory.

## Files

| File | Role |
|---|---|
| `center_r10_M52_S36.json` | Frozen centre (`m=4`, `M=52`, `S=36`, field scale `B`): 252 coefficients `g_{n,s}`, `n ∈ {4,12,…,52}`, `1 ≤ s ≤ 36`, and 7 shape coefficients `c_1, c_9, …, c_49`, as binary64 hex strings interpreted as exact dyadic rationals. |
| `inverse_r10_M52_S36.npy` | Frozen binary64 approximate inverse `Â` (259×259) of the finite Jacobian, entries interpreted exactly; only its certified defect is used. |
| `rank2_cap_core.py` | Arb algebra for `m ∈ {3,4,6}` (disk-polynomial products, `K`, principal tail inverse), shared with `r14/`. |
| `finite_rank2.py` → `finite_bound_r10.json` | Finite Jacobian, inverse defect `‖I − ÂM_f‖`, `Z` of the 259 finite columns, `Y`, `‖A‖` on the finite block. |
| `tail_inverse_rank2.py` → `tail_inverse_bound_r10.json` | Norm of the principal-tail inverse (100 columns in Arb, analytic bound beyond). |
| `g_tail_rank2.py` → `g_far_bound_r10.json`, `g_near_r10_*_*.json` (18 files) | Analytic bound for the far `g` tail; the 358 near `g` tail columns in 18 batches. |
| `shape_tail_rank2.py` → `shape_near_r10_*_*.json` (4 files) | Near shape columns `j = 57, 65, …, 297`, with exact subtraction of the principal part. |
| `far_shape_rank2.py` → `far_shape_bound_r10.json` | Analytic far-shape bound for all `j ≥ 305`, using `Kg° = (1−r²)²Q`. |
| `check_identities_rank2.py` → `identities_r10.json` | Exact SymPy check of the three-term formula for `K` (12 cases); 256-bit checks of the ball principal shape column, of the radial-power coefficients `q_s` (positive, sum one), and of the identity `Kg° = (1−r²)²Q`. |
| `verify_r10.py` → `certificate_r10.json` | Final gate and main certificate. |
| `verify_convex_r10.py` → `convex_r10.json` | Independent convexity / non-disc gate, bound to the main certificate and verifier by SHA-256. |
| `test_fail_closed_r10.py` | Negative tests of the final gate on altered receipts (scratch copy). |
| `reproduce_r10.sh` | Full recomputation of every receipt in a scratch directory, with byte-for-byte comparison. |
| `SHA256SUMS` | Checksums of all files in this directory. |

## Commands and expected output

**1. File integrity.** `sha256sum -c SHA256SUMS` — `OK` on every line.

**2. Final gate** (reads the frozen centre, inverse and 27 receipts; < 1 s):

```bash
OMP_NUM_THREADS=1 python3 verify_r10.py
```

Expected output (exactly):

```text
PROVED
certificate: certificate_r10.json
```

and `certificate_r10.json` is rewritten byte-identically. Before printing
`PROVED` the gate checks the SHA-256 digest of the frozen bundle (centre,
inverse, the seven scripts listed above from `rank2_cap_core.py` to
`check_identities_rank2.py`, and the 27 receipts,
`0c8b38d38b96138ae2d59cee3d35b3bf1eecdd9c737a0a6bf39c6e1e70b66884`), that the
batches cover every column class without gaps, that every interval is finite,
all rational thresholds, the radii polynomial and its derivative, and the
geometric inequalities.

**3. Convexity gate** (< 1 s): `python3 verify_convex_r10.py` prints

```text
PROVED_CONVEX [0.97305365889140304527 +/- 4.47e-21]
```

and rewrites `convex_r10.json` byte-identically.

**4. Full recomputation** (about 105 s on one core):
`bash reproduce_r10.sh` recomputes all 27 receipts from the centre and `Â` in
a fresh temporary directory, reruns steps 2–3 there, and compares every
released JSON file byte for byte; the last line is
`R10: all 30 released JSON files reproduced byte for byte in <dir>`.

**5. Fail-closed tests** (a few seconds): `python3 test_fail_closed_r10.py`
checks that `max_upper`/`argmax_upper` reject NaN and infinite balls, and
that, in a scratch copy of this directory, the final gate accepts the
released receipts but rejects receipts altered by a single field (NaN or
infinite values, a coverage gap, a value above a threshold), both with the
frozen-bundle checksum in force and with it bypassed. Last line:
`FAIL-CLOSED TESTS PASSED`.

## Certified values

Arb upper (resp. lower) endpoints at 128 bits, rounded outward:

| Quantity | Bound | Threshold in the gate |
|---|---:|---:|
| `Y` | `< 9.59187·10⁻⁹` | `10⁻⁸` |
| `Z` | `< 0.552897` | `0.56` |
| — finite / near `g` / far `g` | `0.0606041 / 0.462762 / 0.393351` | |
| — near shape / far shape | `0.552897 / 0.212795` | |
| `‖A‖` (finite block) / tail | `< 10194.8 / < 2.01446` | `10195 / 2.015` |
| inverse defect `‖I − ÂM_f‖` | `< 5.72581·10⁻¹²` | `1` |
| `C₂, …, C₆` | `< 395.405, 0.0711322, 3.23526·10⁻⁷, 1.74697·10⁻¹¹, 3.72354·10⁻¹⁶` | `396, 0.072, 4·10⁻⁷, 2·10⁻¹¹, 4·10⁻¹⁶` |
| radius `r` | `10⁻⁷` | |
| radii polynomial / derivative | `< −3.39960·10⁻⁸ / < −0.439920` | `0 / 0` |
| univalence margin on `|w| ≤ 21/20` | `> 37.9650` | `0` |
| `1 − U/L` (convexity of `D`) | `> 0.973053` | `0` (gate), `0.97` (convexity gate) |
| `|c₉|` of the exact solution | `> 0.0138324` | `0` |

The constants `C_q` use `κ = 1/48` as an Arb enclosure (see `../README.md`);
the exact values are `C₂ ∈ [395.4041767138573008… ± 2.6·10⁻²⁸]` and
`C₃ ∈ [0.07113210996788329736… ± 4·10⁻³²]`.

## Checksums

| File | SHA-256 |
|---|---|
| `verify_r10.py` | `85a1d94f956d74b9f04589a08a982e0e75989963948026ad1ec9238b3859ff20` |
| `certificate_r10.json` | `7e084736a9fd834ae58b26031999262d5f12706357291e8dc31a925c14ca2375` |
| frozen bundle digest | `0c8b38d38b96138ae2d59cee3d35b3bf1eecdd9c737a0a6bf39c6e1e70b66884` |

All other checksums are in `SHA256SUMS`.
