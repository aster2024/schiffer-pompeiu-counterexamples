# Convex counterexamples to the Schiffer and Pompeiu conjectures in dimensions 3, 4, 6, 8, 10 and 14

[![DOI](https://zenodo.org/badge/1391007796.svg)](https://doi.org/10.5281/zenodo.22999505)

This repository contains the paper and the computer-assisted proofs for the following result.

**Theorem.** For every n ∈ {3, 4, 6, 8, 10, 14} there are a bounded **convex** domain Ω ⊂ ℝⁿ, **not a ball**, whose boundary is a real-analytic hypersurface diffeomorphic to Sⁿ⁻¹, and a nonconstant function u, real analytic on a neighbourhood of Ω̄, such that

    Δu + u = 0 in Ω,   u = 1 and ∇u = 0 on ∂Ω.

**Consequences.**
- By Green's identity the Fourier transform of the indicator of Ω vanishes on the unit sphere, so Ω fails the Pompeiu property.
- The Schiffer and Pompeiu conjectures therefore fail in these dimensions, even among convex domains.
- The paper also contains a second, independent three-dimensional example. It is strictly star-shaped and non-convex.

**Structure.**
- **n = 3:** the domain is axisymmetric (O(2) × ℤ₂), obtained from a conformal parametrisation of its meridian section.
- **n = 4, 6, 8, 10, 14:** the domains are adjoint-invariant domains in the rank-two compact Lie algebras u(2), so(4), su(3), so(5) and g₂. Harish-Chandra's radial-part formula reduces each problem to a planar Helmholtz problem, and Kostant's convexity theorem reduces convexity to the Cartan section.
- **Existence proofs:** in each case existence follows from a Newton–Kantorovich (radii-polynomial) argument in a weighted coefficient space. Finite blocks are verified in Arb ball arithmetic, and every infinite tail is controlled by explicit analytic bounds.

## Contents

| Path | Description |
|---|---|
| `paper/main.pdf`, `paper/main.tex` | The paper (87 pages). |
| `verification/README.md` | Layout, requirements, commands and running times for all dimensions. |
| `verification/` (top level), `verification/r6/` | n = 4 and n = 6. |
| `verification/r8/`, `r10/`, `r14/` | n = 8, 10, 14 (each with `README.md` and `reproduce_r*.sh`). |
| `verification/r3_convex/` | n = 3, convex example: single-file verifier `verify_r3_convex.py` and two independent full-run receipts. |
| `verification/r3/` | n = 3, additional non-convex example: verifier `verify_r3.py`, two certificates and a non-convexity check. |
| `SHA256SUMS`, `verification/**/SHA256SUMS` | Checksums. |

The repository was previously named `schiffer-pompeiu-r4`. Its first releases contained only the four- and six-dimensional results.

## Reproducing

The requirements are Python ≥ 3.10, `python-flint` 0.9.0, NumPy and, for some dimensions, SciPy and gcc. See `verification/README.md` and the per-directory READMEs. Running times:
- n = 4, 6, 8, 10, 14: minutes each;
- n = 3 convex: about 1.8 hours on 5 cores;
- n = 3 non-convex: about 1.1 hours on 4 cores.

## Status

These are computer-assisted proofs. They were independently re-run and audited (mathematics, code and independent numerics) before release, but they have **not yet been peer reviewed**.

## Author and acknowledgement

Jizhou Guo (Dots Studio, Rednote). This work was carried out with the assistance of AI models (Claude Opus 5.5, GPT-6 Astra and GPT-6 Sol). The author takes full responsibility for the content.

## Licence

The code is under the MIT licence (see `LICENSE`).
