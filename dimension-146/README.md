# Dimension 146

Jizhou Guo — mitsuha2021b@gmail.com — 9 October 2026

A computer-assisted certificate for one explicit dimension in the two-block setting of the second paper of this
repository ([`infinitely-many-dimensions/`](../infinitely-many-dimensions/)).

**Theorem.** There are a bounded strictly convex domain Ω ⊂ ℝ¹⁴⁶, not a ball, with real-analytic boundary and
invariant under O(54) × O(92), and a nonconstant real function u with

    Δu + u = 0 in Ω,   u = 1 and ∂_ν u = 0 on ∂Ω.

Here p = 73 and R = j_{73,17} = 149.02532977559787… is a zero of the Bessel function J₇₃. The boundary is
ϱ = R(1 + f(x)) with x = |X|²/(|X|² + |Y|²) for (X, Y) ∈ ℝ⁵⁴ × ℝ⁹². The solution lies in the ball of radius 10⁻⁴⁶
around the explicit centre of `center.json`, in the norm of the proof.

## Proof

The analytic argument is [`PROOF.md`](PROOF.md); the same text is embedded in the verifier (its SHA-256,
`332a44a90a063de12f059efa891638cc528e354fc600f4230303c1e1e7478d3e`, is recorded in the final receipt). The numerical
hypotheses are verified by the single-file program `verify_inf3_hardened.py`
(SHA-256 `3abd4552628cd60d084490667887eb9e862d63b457bea42f8ba46022ca94f78c`) in Arb ball arithmetic, with exact
rational identities and, for the large sparse products, floating-point computation with an explicit IEEE binary64
error model. The final phase prints `PROVED_INF3_N146`.

Certified values (from the receipts):

| Quantity | Bound |
|---|---|
| Dirichlet and radial Neumann defects of the centre on the whole boundary | < 1.86 × 10⁻⁷⁵, < 7.99 × 10⁻⁷⁶ |
| residual ‖F(z₀)‖ after exact clamping | < 5.60 × 10⁻⁶⁸ |
| inverse of the graph operator on the deformed domain | < 2.40 × 10⁶ |
| inverse in the original coefficient norm | < 4.30 × 10¹⁰ |
| Newton radius | 10⁻⁴⁶ (exact) |
| radii polynomial | < −9.67 × 10⁻⁴⁷ |
| contraction constant | < 4.42 × 10⁻¹⁰ |
| shape ball on which convexity and the non-ball property are verified | radius 10⁻⁸ |
| Dirichlet spectral gap on O(54) × O(92)-invariant functions | > 9.0 × 10⁻⁹ |

The finite head has angular degrees j ≤ 48 and radial indices s < 128 (dimension 6272); the rectangular computation
covers all inputs with j ≤ 90, s < 214; all further inputs and omitted outputs are covered by explicit bounds for
infinite columns.

## Reproducing

Requirements: Python 3.12 with python-flint, NumPy, SciPy, SymPy and mpmath. Three inputs are needed in an empty
directory: `verify_inf3_hardened.py`, `center.json` and `preconditioner.npy`. The last file (315 MB, SHA-256
`b55f5134707cffa0c5fba5cfa5583284d71e2ca9b56e43bc17114eedc07f411e`) is an approximate inverse of the head matrix; the
verifier treats it as an untrusted proposal and bounds the defect of the product rigorously. It is attached to
[release v5.0](https://github.com/aster2024/schiffer-pompeiu-counterexamples/releases/tag/v5.0) of this repository.

Set `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1` and run the 21 phases in this
order, each as `python verify_inf3_hardened.py PHASE ARGUMENTS`, waiting for each to exit before starting a phase
that depends on it:

```text
center --center center.json --output center_moments.json --bits 4096 --L 180 --max-defect 2e-75 --angular-cut 768
components --center center.json --J 90 --N 214 --shape-J 16 --prefix q16
identities --p 73 --a 54 --output identities.json
geometry --center center.json --radius 1e-8 --bits 4096 --output geometry.json
window --p 73 --R 149.025329775597864 --digits 12 --half-width 4 --bits 768 --output ball_window.json
ball-inverse --p 73 --R 149.025329775597864 --digits 12 --N 128 --J 80 --bits 768 --output ball_inverse.json
build --components q16_components.json --center center.json --head-J 48 --head-N 128 --prefix graph --stage basis
build --components q16_components.json --center center.json --head-J 48 --head-N 128 --prefix graph --stage gradient
build --components q16_components.json --center center.json --head-J 48 --head-N 128 --prefix graph --stage mass
build --components q16_components.json --center center.json --head-J 48 --head-N 128 --prefix graph --stage d4
build --components q16_components.json --center center.json --head-J 48 --head-N 128 --prefix graph --stage md2
build --components q16_components.json --center center.json --head-J 48 --head-N 128 --prefix graph --stage P
build --components q16_components.json --center center.json --head-J 48 --head-N 128 --prefix graph --stage d3
build --components q16_components.json --center center.json --head-J 48 --head-N 128 --prefix graph --stage assemble
products --matrix-prefix graph --inverse preconditioner.npy --output products.json
tails --components q16_components.json --center center.json --head-J 48 --head-N 128 --input-J 90 --input-N 214 --output tails.json
arithmetic --components q16_components.json --output arithmetic.json
inverse --components q16_components.json --center center.json --products products.json --tails tails.json --arithmetic arithmetic.json --output operator_inverse.json
newton --moments center_moments.json --B0 1e7 --B2 1e11 --radius 1e-46 --output newton.json
hilbert --inverse operator_inverse.json --geometry geometry.json --moments center_moments.json --radius 1e-46 --output hilbert_gap.json
final --output THEOREM_RECEIPT.json
```

The argument `--R 149.025329775597864 --digits 12` gives the 12-digit bracket in which the zero of J₇₃ is isolated.
Accept only exit status zero of every phase and status `PROVED_INF3_N146` in `THEOREM_RECEIPT.json`. The complete
run takes about one core-hour; no single phase exceeds 40 minutes.

## Files

| Path | Content |
|---|---|
| `verify_inf3_hardened.py` | The verifier |
| `PROOF.md` | The analytic proof (identical to the text embedded in the verifier) |
| `center.json` | The centre (exact rationals; SHA-256 `9437d79159eb43568d75416a0342d2e56e50cef35a471553479ad0b8d914ecb1`) |
| `receipts/` | Receipts of a complete run of the 21 phases, with the final `THEOREM_RECEIPT.json` |
| `SHA256SUMS` | Checksums |
