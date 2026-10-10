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

## Second paper: infinitely many dimensions

[`infinitely-many-dimensions/`](infinitely-many-dimensions/) contains the paper *Convex non-ball Schiffer domains in infinitely many dimensions* (114 pages, 8 October 2026; [PDF](infinitely-many-dimensions/paper.pdf)).

**Theorem.** There are infinitely many dimensions n, both odd and even and in every residue class modulo every integer, for which ℝⁿ contains a bounded strictly convex domain, not a ball, with real-analytic boundary, carrying a nonconstant solution of Δu + u = 0 with u = 1 and ∇u = 0 on the boundary. The domains are invariant under O(a) × O(n − a), and along each family the number of pairwise non-similar such domains in dimension n is at least (1/2 − o(1))·n.

The dimensions are n = 2p for the sufficiently large integers and half-integers p close to a zero in the order of J_p(2·sqrt((p+1)(p+2))); the precise families are in the [README](infinitely-many-dimensions/README.md) of that directory. The proofs in this paper are analytic.

## Further certified dimensions

- [`dimension-19/`](dimension-19/): a convex counterexample in ℝ¹⁹ (verifier, centre, receipts and the proof of the one bound that differs from the paper). With the paper this covers every dimension from 2 to 21.
- [`dimension-146/`](dimension-146/): a computer-assisted certificate for a strictly convex non-ball domain in ℝ¹⁴⁶ invariant under O(54) × O(92), in the two-block setting of the second paper (verifier, centre, analytic proof and receipts).

## Berenstein's problem in dimensions four and fourteen

[`berenstein-dimensions-4-14/`](berenstein-dimensions-4-14/) contains the paper *Counterexamples to Berenstein's conjecture in dimensions four and fourteen* (17 pages, 10 October 2026; [PDF](berenstein-dimensions-4-14/paper.pdf)) and its verification files.

**Theorem.** There are a bounded strictly convex domain in ℝ⁴, different from a ball, with real-analytic boundary diffeomorphic to S³, and a bounded convex domain in ℝ¹⁴, different from a ball and invariant under the adjoint action of G₂, each carrying a sign-changing solution of

    Δu + u = 0 in Ω,   u = 0 and ∂ᵥu = c on ∂Ω,   c ≠ 0.

The proofs are computer-assisted, by the fixed-disc method of the planar case; each dimension has a single-file verifier whose last line is `PROVED`.

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
| `infinitely-many-dimensions/` | The second paper (114 pages): PDF, LaTeX source and checksums. |
| `dimension-19/`, `dimension-146/` | Verifiers, centres, proofs and receipts for dimensions 19 and 146. |
| `berenstein-dimensions-4-14/` | Berenstein's problem in dimensions 4 and 14: paper (17 pages), LaTeX source, verifiers, certificate data and checksums. |

The directory `verification/` is identical to the ancillary directory `anc/` of the arXiv submission. The repository was previously named `schiffer-pompeiu-r4`; its first releases contained only the four- and six-dimensional results, and release 3.0 the dimensions 3, 4, 6, 8, 10 and 14.

## Reproducing

The requirements are Python ≥ 3.10, `python-flint` 0.9.0, NumPy and, for some dimensions, SciPy, mpmath and gcc. See `verification/README.md` and the per-directory READMEs for commands and running times.

## Status

The proofs of the first paper are computer-assisted; the proofs of the second paper are analytic.

## Author and acknowledgement

Jizhou Guo. The author conceived and directed the project and takes full responsibility for the content. The constructions, proofs, verification code and text were produced with AI models (Claude Opus 5.5, Claude Sonnet 5.5, GPT-6 Astra, GPT-6 Sol and GPT-6.1 Sol) under the author's direction, and were audited with AI assistance.

## Licence

The code is under the MIT licence (see `LICENSE`).
