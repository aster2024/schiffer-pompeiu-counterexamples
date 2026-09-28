# R⁸ = su(3): computer-assisted proof package

Computer-assisted part of the proof that there is a bounded **convex** domain
`Ω ⊂ su(3) ≅ R⁸`, invariant under the adjoint action of `SU(3)`, with
real-analytic boundary diffeomorphic to `S⁷`, not a ball, carrying a
nonconstant solution of

```text
Δu + u = 0 in Ω,    u = 1,  ∇u = 0 on ∂Ω .
```

Via the radial part of the Laplacian on the Cartan plane (Weyl group `I₂(3)`,
product of positive roots `π = Im z³`) this reduces to a planar Helmholtz
problem on `D = ψ(𝔻)`, `ψ(w) = Σ_{j ≡ 1 (mod 6)} c_j w^j`, with Cauchy data
`(π, ∇π)` on `∂D`, written as `H(g,c) = g + |ψ'|²(Kg + Im ψ³) = 0` in the
disk-polynomial sine sector `n ≡ 3 (mod 6)`. The scripts below prove the
existence of an exact zero of `H` near a frozen centre (contraction argument
in a weighted ℓ¹ space with weight `ρ = 11/10`) and certify that `D` is
strictly convex, contains `0`, and is not a disc. The Lie-theoretic reduction,
the lifting to `R⁸` and the Pompeiu corollary are proved in the paper.
Rank-two notes common to `r8/`, `r10/`, `r14/` are in `../README.md`.

## Requirements

Python 3 (runs used Python 3.10.19) with **python-flint 0.9.0** and NumPy;
SymPy for `check_identities_r8.py`. One thread. All commands are run from this
directory.

## Files

| File | Role |
|---|---|
| `center_r8_M45_S24.json` | Frozen centre (`m=3`, `M=45`, `S=24`): 192 coefficients `g_{n,s}`, `n ∈ {3,9,…,45}`, `1 ≤ s ≤ 24`, and 8 shape coefficients `c_1, c_7, …, c_43`, as binary64 hex strings interpreted as exact dyadic rationals. |
| `inverse_r8_M45_S24.npy` | Frozen binary64 approximate inverse `Â` (200×200) of the finite Jacobian, entries interpreted exactly; only its certified defect is used. |
| `verify_r8.py` | Arb algebra (disk-polynomial products, `K`, principal tail inverse) used by all stage scripts, and the final gate (`--stage final`). |
| `finite_interval_r8.py` → `finite_bound_r8.json` | Finite Jacobian, inverse defect `‖I − ÂM_f‖`, `Z` of the 200 finite columns, `Y`, `‖A‖` on the finite block. |
| `tail_inverse_r8.py` → `tail_inverse_bound_r8.json` | Norm of the principal-tail inverse (100 columns in Arb, analytic bound beyond). |
| `g_tail_r8.py` → `g_far_bound_r8.json`, `g_near_r8_*_*.json` (20 files) | Analytic bound for the far `g` tail; the 396 near `g` tail columns in 20 batches. |
| `shape_tail_r8.py` → `shape_near_r8_049_133.json`, `shape_near_r8_139_217.json` | Near shape columns `j = 49, 55, …, 217`, with exact subtraction of the principal part. |
| `far_shape_r8.py` → `far_shape_bound_r8.json` | Analytic far-shape bound for all `j ≥ 223`, using `Kg° = (1−r²)²Q`. |
| `check_identities_r8.py` → `identities_r8.json` | Exact SymPy check of the three-term formula for `K` (12 cases); 256-bit checks of the ball principal shape column, of the radial-power coefficients `q_s` (positive, sum one), and of the identity `Kg° = (1−r²)²Q`. |
| `certificate_r8.json` | Main certificate, written by the final gate. |
| `verify_convex_r8.py` → `convex_r8.json` | Independent convexity / non-disc gate, bound to the main certificate and verifier by SHA-256. |
| `test_fail_closed_r8.py` | Negative tests of the final gate on altered receipts (scratch copy). |
| `reproduce_r8.sh` | Full recomputation of every receipt in a scratch directory, with byte-for-byte comparison. |
| `SHA256SUMS` | Checksums of all files in this directory. |

## Commands and expected output

**1. File integrity.** `sha256sum -c SHA256SUMS` — `OK` on every line.

**2. Final gate** (reads the frozen centre, inverse and 27 receipts; < 1 s):

```bash
OMP_NUM_THREADS=1 python3 verify_r8.py --stage final \
    --centre center_r8_M45_S24.json --output certificate_r8.json
```

Expected output (exactly):

```text
algebra 200 residual coefficients 1319 residual norm [1.50452988436779e-5 +/- 9.63e-21]
PROVED
certificate: certificate_r8.json
```

and `certificate_r8.json` is rewritten byte-identically. Before printing
`PROVED` the gate checks the SHA-256 digest of the frozen bundle (centre,
inverse, the six stage scripts and the 27 receipts,
`c32ce586dafc332fdf8d05d74a90a98babe97a4275d86529946f482c06d5f21e`), that the
batches cover every column class without gaps, that every interval is finite,
all rational thresholds, the radii polynomial and its derivative, and the
geometric inequalities.

**3. Convexity gate** (< 1 s): `python3 verify_convex_r8.py` prints

```text
PROVED_CONVEX [0.83153153907625301655 +/- 2.45e-21]
```

and rewrites `convex_r8.json` byte-identically.

**4. Full recomputation** (about 70 s on one core):
`bash reproduce_r8.sh` recomputes all 27 receipts from the centre and `Â` in
a fresh temporary directory, reruns steps 2–3 there, and compares every
released JSON file byte for byte; the last line is
`R8: all 30 released JSON files reproduced byte for byte in <dir>`.

**5. Fail-closed tests** (a few seconds): `python3 test_fail_closed_r8.py`
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
| `Y` | `< 1.51945·10⁻⁵` | `1.53·10⁻⁵` |
| `Z` | `< 0.418572` | `0.43` |
| — finite / near `g` / far `g` | `0.0522617 / 0.194306 / 0.179735` | |
| — near shape / far shape | `0.284628 / 0.418572` | |
| `‖A‖` (finite block) / tail | `< 963.778 / < 2.30204` | `964 / 2.303` |
| inverse defect `‖I − ÂM_f‖` | `< 2.28573·10⁻¹⁴` | `1` |
| `C₂, C₃, C₄, C₅` | `< 0.00721911, 4.88354·10⁻¹⁰, 1.54908·10⁻¹⁸, 2.11513·10⁻²⁶` | `1/125, 5·10⁻¹⁰, 2·10⁻¹⁸, 3·10⁻²⁶` |
| radius `r` | `1/30000` | |
| radii polynomial / derivative | `< −3.69999·10⁻⁶ / < −0.569999` | `0 / 0` |
| univalence margin on `|w| ≤ 21/20` | `> 20.0371` | `0` |
| `1 − U/L` (convexity of `D`) | `> 0.831531` | `0` (gate), `3/4` (convexity gate) |
| `|c₇|` of the exact solution | `> 0.0774840` | `0` |

## Checksums

| File | SHA-256 |
|---|---|
| `verify_r8.py` | `c1f3b8ae24b484b6a1ea8b4a016030d9b621bc965ed5e2bfaa3990a9bf5fa586` |
| `certificate_r8.json` | `49d813f5aa6d5b9cd8246b5ef9cd5365adaf11cc4d7e944db1a2f415e6b1a7a1` |
| frozen bundle digest | `c32ce586dafc332fdf8d05d74a90a98babe97a4275d86529946f482c06d5f21e` |

All other checksums are in `SHA256SUMS`.
