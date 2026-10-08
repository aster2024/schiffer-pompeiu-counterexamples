# Dimension 19

Jizhou Guo — mitsuha2021b@gmail.com — 9 October 2026

Together with the paper of this repository (arXiv:[2609.35419](https://arxiv.org/abs/2609.35419), which covers
dimensions 2 to 18, 20 and 21), this gives convex counterexamples in every dimension from 2 to 21.

**Theorem.** There are a bounded strictly convex domain Ω ⊂ ℝ¹⁹, not a ball, with real-analytic boundary and
invariant under O(18) × ℤ₂, and a nonconstant real function u with

    Δu + u = 0 in Ω,   u = 1 and ∂_ν u = 0 on ∂Ω.

Consequently the Fourier transform of the indicator of Ω vanishes on the unit sphere and Ω fails the Pompeiu
property.

## Proof

The proof is the argument of the part "The exact inverse in general dimension" of the paper, with one bound replaced:
the norm of Q − b in the gate for the reciprocal of Q is estimated by the radial bound proved in
[`Q_RADIAL_PROOF.md`](Q_RADIAL_PROOF.md). All other estimates, thresholds and truncations are those of the paper.
The numerical hypotheses are verified by the single-file program `verify_axisym_exact_v4.py`
(SHA-256 `bcdc648e261d054808a6bf3d21f8e63e3b17070e0d28597a1f88c5923db5dab5`), which embeds the dimension parameters,
the stage sources and the centres.

Certified values for n = 19 (existence ball ‖g − g₀‖_Y + σ‖p − p₀‖_W ≤ r):

| Quantity | Value |
|---|---|
| finite modes of the centre | 1217 |
| frequency scale ρ | 21.42848697211536 |
| R, σ, r | 1.03, 10⁴⁶, 10⁻⁸⁰ |
| Y₀ | < 3.991919 × 10⁻⁸⁵ |
| contraction constant | < 2.750390 × 10⁻¹³ |
| 4 C₂ Y₀ | < 2.195867 × 10⁻¹⁷ |
| radii polynomial | < −9.999600 × 10⁻⁸¹ |
| analytic univalence | > 0.8598170 |
| meridian and transverse curvature numerators | > 0.3266036 |
| non-ball witness (absolute w³ coefficient) | > 1.059923 × 10⁻² |
| L² inverse, spectral states | < 8.389601, 33 |

The curvature numerators, the univalence bound and the non-ball witness hold on the whole existence ball.

## Reproducing

Requirements: Python 3.10 with the packages in `REQUIREMENTS.txt`. In a new directory containing only the verifier:

```bash
nice -n 19 env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  python -E -O -B verify_axisym_exact_v4.py --dimension 19 --workdir result
```

Accept only exit status zero, the final line `PROVED dimensions=19`, and status `PROVED` in both
`result/VERIFICATION_RECEIPT.json` and `result/n19/VERIFICATION_RECEIPT.json`. The run takes about 35 minutes on one
core. With `--dimension N` for N in 5, 9, 11, 12, 13, 15, 16, 17, 18, 20, 21 the same file reproduces every
mathematical field of the receipts of the paper.

## Files

| Path | Content |
|---|---|
| `verify_axisym_exact_v4.py` | The verifier |
| `Q_RADIAL_PROOF.md` | Proof of the radial bound |
| `frozen_centers/n19.json.gz` | The centre (exact decimals; `gzip -dc` returns the JSON input, SHA-256 `3ae81c98befb86fa6eeeaf3862a51c8da8ea8d2364ec78260a33c94785b3f564`) |
| `receipts/` | Receipts of the final isolated run for n = 19 |
| `REQUIREMENTS.txt`, `SHA256SUMS` | Package versions; checksums |
