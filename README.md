# Convex counterexamples to the Schiffer and Pompeiu conjectures in dimensions two to eighteen

[![DOI](https://zenodo.org/badge/1391007796.svg)](https://doi.org/10.5281/zenodo.22999505) [![arXiv](https://img.shields.io/badge/arXiv-2609.35419-b31b1b.svg)](https://arxiv.org/abs/2609.35419)

This repository contains the paper (arXiv:[2609.35419](https://arxiv.org/abs/2609.35419), version 2) and the computer-assisted proofs for the following results.

**Theorem.** For every n with 2 ≤ n ≤ 18, and for n = 20 and n = 21, there are a bounded **convex** domain Ω ⊂ ℝⁿ, **not a ball**, whose boundary is a real-analytic hypersurface diffeomorphic to Sⁿ⁻¹, and a nonconstant function u, real analytic on a neighbourhood of Ω̄, such that

    Δu + u = 0 in Ω,   u = 1 and ∇u = 0 on ∂Ω.

For n = 2 the domain is strictly convex; it is the first convex planar counterexample.

**Consequences.**
- By Green's identity the Fourier transform of the indicator of Ω vanishes on the unit sphere, so Ω fails the Pompeiu property.
- The Schiffer and Pompeiu conjectures therefore fail in these dimensions, even among convex domains.

**Berenstein's problem.** The paper also constructs a strictly convex non-disc planar domain with real-analytic boundary carrying a sign-changing solution of Δu + u = 0 with u = 0 and ∂ᵥu = 1 on the boundary, the first convex planar counterexample to Berenstein's conjecture. Planar existence for this problem is due to Colbrook, Sadeghi and Stepaniants.

**Structure.**
- **n = 2:** a D₁₈₂-symmetric domain, obtained through the exact two-sided inverse of the residual-corrected linearisation (a Dirichlet Helmholtz inverse, a holomorphic boundary real-part problem and the material-derivative identity).
- **n = 4, 6, 8, 10, 14:** adjoint-invariant domains in the rank-two compact Lie algebras u(2), so(4), su(3), so(5) and g₂. Harish-Chandra's radial-part formula reduces each problem to a planar Helmholtz problem, and Kostant's convexity theorem reduces convexity to the Cartan section.
- **n = 3, 5, 7:** axisymmetric domains (O(n−1) × ℤ₂), by a quartic equation in weighted coefficient spaces.
- **n = 5 (second proof), 9, 11, 12, 13, 15, 16, 17, 18, 20, 21:** axisymmetric domains, by one dimension-parametrised argument based on the exact inverse; one single-file verifier.
- **Existence proofs:** Newton–Kantorovich (radii-polynomial) arguments. Finite blocks are verified in Arb ball arithmetic, and every infinite tail is controlled by explicit analytic bounds.

## Contents

| Path | Description |
|---|---|
| `paper/main.pdf`, `paper/*.tex` | The paper (196 pages). |
| `verification/README.md` | Layout, requirements, commands and running times. |
| `verification/r2/` | n = 2. |
| `verification/` (top level), `verification/r6/` | n = 4 and n = 6. |
| `verification/r8/`, `r10/`, `r14/` | n = 8, 10, 14. |
| `verification/r3_convex/`, `verification/r3/` | n = 3: the convex example, and an additional non-convex example. |
| `verification/r5/`, `verification/r7/` | n = 5 and n = 7. |
| `verification/gd/` | n = 5, 9, 11, 12, 13, 15, 16, 17, 18, 20, 21: `verify_axisym_exact_v2.py`, frozen centres and receipts. |
| `verification/berenstein/`, `verification/diagnostics/berenstein/` | The convex planar domain for Berenstein's problem: `verify_berenstein_planar_convex_v3.py`, centre and receipt; numerical diagnostics. |
| `SHA256SUMS`, `verification/**/SHA256SUMS` | Checksums. |

The directory `verification/` is identical to the ancillary directory `anc/` of the arXiv submission. The repository was previously named `schiffer-pompeiu-r4`; its first releases contained only the four- and six-dimensional results, and release 3.0 the dimensions 3, 4, 6, 8, 10 and 14.

## Reproducing

The requirements are Python ≥ 3.10, `python-flint` 0.9.0, NumPy and, for some dimensions, SciPy, mpmath and gcc. See `verification/README.md` and the per-directory READMEs for commands and running times.

## Status

These are computer-assisted proofs. The paper has not yet been peer reviewed.

## Author and acknowledgement

Jizhou Guo (Dots Studio, Rednote). The author conceived and directed the project and takes full responsibility for the content. The constructions, proofs, verification code and text were produced with AI models (Claude Opus 5.5, Claude Sonnet 5.5, GPT-6 Astra, GPT-6 Sol and GPT-6.1 Sol) under the author's direction, and were audited with AI assistance.

## Licence

The code is under the MIT licence (see `LICENSE`).
