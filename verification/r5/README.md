# Convex five-dimensional certificate

This directory contains the computer-assisted part of the proof of
Theorem 1.1 for `n = 5` (Section 22 of the paper): a bounded, convex,
non-ball domain in `R⁵`, invariant under `O(4)` acting on the last four
coordinates and under `X₁ ↦ −X₁`, with a nonconstant solution of
`Δu + u = 0`, `u = 1` and `∇u = 0` on the boundary. The contraction
argument is carried out in coefficient spaces graded by the polynomial
degree, with weight `Υ(D) = 1 + D/10⁶`, `R = 21/20` and shape scale
`σ = 30ρ²`.

| File | Content |
|---|---|
| `verify_r5.py` | The unmodified single-file verifier, SHA-256 `b49dc914721b4f62d5e6ddda2b9ee81a872ed0f77d7a0a42b228dc2a5b040851`. It embeds, as a compressed payload with SHA-256 `83816455a7876dfba33d7fd8f7fbbbfe360009d3c2c934cc601465170362ae5d`, the exact dyadic centre (SHA-256 `c672b3fb2c33ee555903f68cfc2d6b2fb3701877caf6fc27ccc604af1a805b40`), the frozen class parameters, the stage programs and the C kernels. Nothing computed is read from disk. |
| `VERIFICATION_RECEIPT.json` | Receipt of the first full run of this verifier (`PROVED`, 2364 s with six single-threaded workers). Only the machine-specific `replay_workdir` string was replaced by `<workdir>`. |
| `params.json` | The frozen class parameters (radial classes `s ≥ S_ℓ` for `ℓ ≤ 160`, angular gap classes for `162 ≤ ℓ ≤ 254`, uniform class `ℓ ≥ 256`, explicit shape range `163 ≤ j ≤ 599`, far class `j ≥ 601` with `N* = 600`, `N₂ = 1000`, `I = 140`, radius `r = 10⁻¹⁰`). It is byte-identical to the copy embedded in the verifier (SHA-256 `9ebeda4b4621cdecc13c170dda559406f375ee9f486a56e830bab9dd3280511c`) and is distributed for reference only. |
| `SHA256SUMS` | SHA-256 of the other files in this directory. |

## Command

From this directory, with a new or empty scratch directory:

```bash
python3 -O verify_r5.py --workdir <empty dir> --output receipt.json --jobs 6
```

Requirements: 64-bit Linux (the long-double kernels check the x87 format
and rounding mode), Python 3, **python-flint 0.9.0**, NumPy and `gcc`.
The recorded run took 2364 s, about 40 minutes, with six single-threaded
workers (`--jobs` is at most 6), and needs about 0.6 GB of disk space in
the work directory. On success the last line printed is `PROVED` and the
receipt is written; any failure prints `NOT_CERTIFIED` and exits with
status 2. The verifier refuses a nonempty work directory. Check file
integrity with `sha256sum -c SHA256SUMS`.

## What is checked

In 23 stages the verifier compiles its kernels and runs exact self-tests,
then recomputes at the exact centre: the full-support residual; all 4698
finite columns (4617 field, 81 shape) of the linearisation, the 59
exceptional field columns and the 219 explicit shape columns
`163 ≤ j ≤ 599`, each with its complete polynomial support in ball
arithmetic; the approximate inverse and its certified finite defect; the
weighted boundary Toeplitz model of the shape tail and the column
profiles of the approximate inverse; the envelopes of all field tail
classes; the all-index bound for the shape columns `j ≥ 601`; the Taylor
constants; the geometry of the whole existence ball; and the coverage of
every input column by exactly one class.

## Certified values (outward roundings of the receipt)

| Quantity | Bound |
|---|---|
| `Y = ‖A F(x°)‖` | `< 3.19608e-12` |
| `Z = ‖I − A DF(x°)‖` | `< 0.887235` (attained by the uniform angular class `ℓ ≥ 256`) |
| finite columns / exceptional field columns | `< 0.0585324` / `< 0.496612` |
| explicit shape columns / far shape class | `< 0.753050` / `< 0.605269` |
| `‖A‖`, `A_t` | `< 104568.72`, `< 16.60796` |
| `c₂`, `c₃`, `c₄` | `< 1.540189`, `< 6.608341e-6`, `< 1.582007e-11` |
| radius `r` | `10⁻¹⁰` (exact) |
| radii polynomial, contraction bound | `< −8.0788e-12`, `< 0.887267` |
| `Re ψ′` on the closed disc of radius 101/100 | `> 0.981143` |
| `Re(1 + wψ″/ψ′)` on the closed unit disc | `> 0.839460` |
| `Re(wψ′/ψ)` on the unit circle | `> 0.984794` |
| modulus of `c₃` (non-ball witness) | `> 4.45696e-4` |

## Relation with an earlier build

An earlier build of the verifier (not distributed) differed only in the
far-shape stage: for `N > N₂ = 1000` its uniform bound started the
geometric tails one power of `ζ = R⁻²` too late, so that edge terms of
total size below `2e-16` were not covered (the uniform bound is `0.41233`,
the far-class maximum `0.60527` and the binding `Z` `0.88723`). In the
distributed build every tail starts at `ζ^{k₀}`, `k₀ = N₂/2`, with weight 2
for the `H` and `D` sums and 4 for the `Z` sums. It was rebuilt and run
from scratch; all mathematical fields of the receipt are unchanged.

A second from-scratch run of the same unmodified verifier in a separate empty directory also printed `PROVED` (2355 s); its receipt `VERIFICATION_RECEIPT_replay2.json` agrees with `VERIFICATION_RECEIPT.json` in every mathematical field (only the run time, stage order and work directory differ).
