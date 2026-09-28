# R¹⁴ = g₂: computer-assisted proof package

Computer-assisted part of the proof that there is a bounded **convex** domain
`Ω ⊂ g₂ ≅ R¹⁴`, invariant under the adjoint action of the compact group
`G₂`, with real-analytic boundary diffeomorphic to `S¹³`, not a ball, carrying
a nonconstant solution of

```text
Δu + u = 0 in Ω,    u = 1,  ∇u = 0 on ∂Ω .
```

Via the radial part of the Laplacian on the Cartan plane (Weyl group `I₂(6)`,
product of positive roots `π = Im z⁶` after normalisation) this reduces to a
planar Helmholtz problem on `D = ψ(𝔻)`, `ψ(w) = Σ_{j ≡ 1 (mod 12)} c_j w^j`,
with Cauchy data `(π, ∇π)` on `∂D`, written as
`H(g,c) = g + |ψ'|²(Kg + Im ψ⁶ / B) = 0` in the disk-polynomial sine sector
`n ≡ 6 (mod 12)`. Here `B > 0` is a fixed field scale (an exact dyadic stored
in the centre file, `B ≈ 4.7962817·10¹¹`); the physical field is
`F = B·(Kg + Im ψ⁶/B)∘ψ⁻¹`, so `B` does not change the Cauchy data. The
scripts below prove the existence of an exact zero of `H` near a frozen centre
(contraction argument in a weighted ℓ¹ space with weight `ρ = 11/10`) and
certify that `D` is strictly convex, contains `0`, and is not a disc. The
Lie-theoretic reduction, the lifting to `R¹⁴` and the Pompeiu corollary are
proved in the paper. Rank-two notes common to `r8/`, `r10/`, `r14/`,
including the changes made to the reviewed code, are in `../README.md`.

## Requirements

Python 3 (runs used Python 3.10.19) with **python-flint 0.9.0** and NumPy;
SymPy for `check_identities_rank2.py`; SciPy only for the optional
provenance script `make_inverse_rank2.py`. One thread. All commands are run
from this directory.

## Files

| File | Role |
|---|---|
| `center_r14_M114_S60.json` | Frozen centre (`m=6`, `M=114`, `S=60`, field scale `B`): 600 coefficients `g_{n,s}`, `n ∈ {6,18,…,114}`, `1 ≤ s ≤ 60`, and 10 shape coefficients `c_1, c_13, …, c_109`, as binary64 hex strings interpreted as exact dyadic rationals. |
| `inverse_r14_M114_S60.npy` | Frozen binary64 approximate inverse `Â` (610×610) of the finite Jacobian, entries interpreted exactly; only its certified defect is used. |
| `rank2_cap_core.py` | Arb algebra for `m ∈ {3,4,6}` (disk-polynomial products, `K`, principal tail inverse), shared with `r10/`. |
| `finite_batch_rank2.py` → `finite_r14_*_*.json` (26 files) | The 610 finite Jacobian columns in 26 batches: inverse defect and `Z` per batch. |
| `aggregate_finite_r14.py` → `finite_bound_r14_M114_S60.json` | Aggregates the 26 batches (coverage checked) and computes `Y` and `‖A‖` on the finite block. `finite_rank2.py` supplies the inverse loader. |
| `tail_inverse_rank2.py` → `tail_inverse_bound_r14_M114_S60.json` | Norm of the principal-tail inverse (100 columns in Arb, analytic bound beyond). |
| `g_tail_rank2.py` → `g_far_bound_r14_M114_S60.json`, `g_near_r14_*_*.json` (13 files) | Analytic bound for the far `g` tail; the 1215 near `g` tail columns in 13 batches. |
| `shape_tail_rank2.py` → `shape_near_r14_*_*.json` (8 files) | Near shape columns `j = 121, 133, …, 229`, with exact subtraction of the principal part. |
| `far_shape_split_r14.py` → `far_shape_split_r14_j241.json` | Analytic far-shape bound for all `j ≥ 241`, using `Kg° = (1−r²)²Q`; each output frequency keeps its own denominator, and the small low-frequency part is charged to the full `‖A‖`. Uses `q_of_g`/`check_q` from `far_shape_rank2.py`. |
| `check_identities_rank2.py` → `identities_r14_M114_S60.json` | Exact SymPy check of the three-term formula for `K` (12 cases); 256-bit checks of the ball principal shape column, of the radial-power coefficients `q_s` (positive, sum one), and of the identity `Kg° = (1−r²)²Q`. |
| `verify_r14.py` → `certificate_r14.json` | Final gate and main certificate. |
| `verify_convex_r14.py` → `convex_r14.json` | Independent convexity / non-disc gate, bound to the main certificate and verifier by SHA-256. |
| `test_fail_closed_r14.py` | Negative tests of the final gate on altered receipts (scratch copy). |
| `reproduce_r14.sh` | Full recomputation of every receipt in a scratch directory, with byte-for-byte comparison. |
| `make_inverse_rank2.py`, `disc_rank2.py`, `conformal_r14_D12_M114_S60_scaled.npz` | Provenance only: the binary64 floating-point computation that produced `Â` (`python3 make_inverse_rank2.py --conformal conformal_r14_D12_M114_S60_scaled.npz --output <file>`; the bytes depend on the BLAS). Not a proof input; `make_inverse_rank2.py` is part of the frozen bundle digest. |
| `SHA256SUMS` | Checksums of all files in this directory. |

## Commands and expected output

**1. File integrity.** `sha256sum -c SHA256SUMS` — `OK` on every line.

**2. Final gate** (reads the frozen centre, inverse and 52 receipts; a few
seconds):

```bash
OMP_NUM_THREADS=1 python3 verify_r14.py
```

Expected output (exactly):

```text
PROVED
certificate: certificate_r14.json
```

and `certificate_r14.json` is rewritten byte-identically. Before printing
`PROVED` the gate checks the SHA-256 digest of the frozen bundle (centre,
inverse, the eleven scripts `rank2_cap_core.py`, `make_inverse_rank2.py`,
`finite_rank2.py`, `finite_batch_rank2.py`, `aggregate_finite_r14.py`,
`tail_inverse_rank2.py`, `g_tail_rank2.py`, `shape_tail_rank2.py`,
`far_shape_rank2.py`, `far_shape_split_r14.py`, `check_identities_rank2.py`,
and the 52 receipts,
`3d6973f5391d99816ab1621bc84b1be144dbbc15f730cd097cc31ee5922156d8`), that the
batches cover every column class without gaps or overlaps, that every
interval is finite, all rational thresholds, the radii polynomial and its
derivative, and the geometric inequalities.

**3. Convexity gate** (< 1 s): `python3 verify_convex_r14.py` prints

```text
PROVED_CONVEX [0.89406714320973174750 +/- 4.79e-21]
```

and rewrites `convex_r14.json` byte-identically.

**4. Full recomputation** (about 740 s on one core):
`bash reproduce_r14.sh` recomputes all 52 receipts from the centre and `Â`
in a fresh temporary directory, reruns steps 2–3 there, and compares every
released JSON file byte for byte; the last line is
`R14: all 55 released JSON files reproduced byte for byte in <dir>`.

**5. Fail-closed tests** (about 30 s): `python3 test_fail_closed_r14.py`
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
| `Y` | `< 4.59220·10⁻⁵` | `4.6·10⁻⁵` |
| `Z` | `< 0.604948` | `0.61` |
| — finite / near `g` / far `g` | `0.120542 / 0.576826 / 0.393946` | |
| — near shape / far shape | `0.600896 / 0.604948` | |
| `‖A‖` (finite block) / tail | `< 11766.9 / < 1.66118` | `11767 / 1.662` |
| inverse defect `‖I − ÂM_f‖` | `< 5.07084·10⁻⁹` | `1` |
| `C₂, …, C₈` | `< 210.533, 0.00606095, 1.91343·10⁻⁹, 2.27176·10⁻¹⁴, 1.59420·10⁻¹⁹, 6.13260·10⁻²⁵, 9.99621·10⁻³¹` | `211, 0.0061, 2·10⁻⁹, 3·10⁻¹⁴, 2·10⁻¹⁹, 7·10⁻²⁵, 2·10⁻³⁰` |
| radius `r` | `1/6000` | |
| radii polynomial / derivative | `< −1.31388·10⁻⁵ / < −0.319666` | `0 / 0` |
| univalence margin on `|w| ≤ 21/20` | `> 87.1340` | `0` |
| `1 − U/L` (convexity of `D`) | `> 0.894067` | `0` (gate), `0.89` (convexity gate) |
| `|c₁₃|` of the exact solution | `> 0.0463086` | `0` |

The constants `C_q` use `κ = 1/80` as an Arb enclosure (see `../README.md`).

## Checksums

| File | SHA-256 |
|---|---|
| `verify_r14.py` | `15152fa2368ceb9ade5d96fa9f695c598e43c803e140af07393e34bb6b3603d2` |
| `certificate_r14.json` | `f565538e828fc581b18a7387e668bc51e769f1dbaf2542c07e04710c0db7f27a` |
| frozen bundle digest | `3d6973f5391d99816ab1621bc84b1be144dbbc15f730cd097cc31ee5922156d8` |

All other checksums are in `SHA256SUMS`.
