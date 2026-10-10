# Counterexamples to Berenstein's conjecture in dimensions four and fourteen

Jizhou Guo — mitsuha2021b@gmail.com — 10 October 2026

[**Paper (PDF, 17 pages)**](paper.pdf) · [LaTeX source](source/) · [Verification](verification/)

## Results

Berenstein's problem asks whether a bounded smooth domain Ω ⊂ ℝⁿ that carries a solution of

    Δu + λu = 0 in Ω,   u = 0 and ∂_ν u = c on ∂Ω,   λ > 0, c ≠ 0,

must be a ball. For a solution of one sign this follows from Serrin's symmetry theorem; the solutions below
change sign.

**Theorem 1.1.** There are a bounded strictly convex domain Ω ⊂ ℝ⁴, different from a ball, with real-analytic
boundary diffeomorphic to S³, and a real-analytic sign-changing function u, analytic on a neighbourhood of the
closure of Ω, such that Δu + u = 0 in Ω, u = 0 and ∂_ν u = 1 on ∂Ω.

**Theorem 1.2.** There is a bounded convex domain Ω ⊂ ℝ¹⁴, different from a ball, carrying a sign-changing
solution of Δu + u = 0 with u = 0 and non-zero constant normal derivative on its boundary. After an isometric
identification ℝ¹⁴ ≅ 𝔤₂ the domain is invariant under the adjoint action of the compact group G₂, and its Cartan
section is strictly convex.

In dimension fourteen strict convexity is asserted for the Cartan section and convexity for the domain.
Dimensions six, eight and ten are not covered by these certificates.

## Proof

- Dimension four: an O(3)-invariant domain; on the meridian plane domain, which is symmetric under both
  coordinate reflections, the function y·u is odd in y and solves a planar Helmholtz problem with Cauchy data
  (0, y).
- Dimension fourteen: an Ad(G₂)-invariant domain; with π(z) = Im z⁶, a non-zero scalar multiple of the product of the
  six positive roots on the Cartan plane, the function π·u solves a planar Helmholtz problem with Cauchy data (0, π). Convexity is
  transferred from the Cartan section by the convexity theorem of Kostant, in the form given by Lewis.
- Both planar problems are pulled back to the unit disc by a conformal map. An exact solution is obtained from
  the contraction criterion of Section 4 of the paper, applied to x ↦ x − A·F(x) in a ball around an explicit
  finite centre, where A is an approximate inverse of the derivative. Finite blocks are enclosed in Arb interval
  arithmetic; the infinite tails are bounded by analytic estimates proved in the paper, whose finite constants
  are enclosed in the same way.

| Quantity | Dimension four | Dimension fourteen |
|---|---|---|
| radius r of the contraction ball | 1/12500 | 10⁻²⁰ |
| residual bound Y | < 3.038 × 10⁻⁶ | < 1.063 × 10⁻²⁵ |
| Z₀ | < 0.926 | < 0.49 |
| Z₁ | 0 | < 2.4 × 10⁻¹¹ |
| Z₂(r) | < 0.042881 | < 4 × 10⁻¹⁴ |
| norm of the approximate inverse | < 2444 | < 10⁵ |

## Verification

Requirements: Python 3 (recorded runs: 3.10.19), `python-flint` 0.9.0, and NumPy for dimension fourteen. Run one
verifier at a time.

    cd verification/r4
    sha256sum -c SHA256SUMS
    python3 -B verify_berenstein_r4.py          # about 10 s; last line: PROVED
    python3 -B check_inequalities.py            # last line: ALL_PAPER_INEQUALITIES_CHECKED

    cd ../r14
    sha256sum -c SHA256SUMS
    python3 -B verify_r14.py                    # about 35 s; last line: PROVED
    python3 -B check_inequalities.py            # last line: ALL_PAPER_INEQUALITIES_CHECKED

Each single-file verifier embeds its programs and certificate data; the same files are also given separately
in the directory. `verify_berenstein_r4.py --recompute` (about 6 minutes) and `verify_r14.py --recompute-all`
(about 9 minutes) regenerate every bounded-column certificate and compare it byte for byte with the supplied
one. The verifiers refuse to run under `python -O`. The READMEs in the two directories give the details.

SHA-256 of the single-file verifiers: `e3479208ba7234e291cdde4d30e8e8c66ece07ae3bdade67daf9661f622a6337`
(dimension four) and `3ea6375202d31df994cb407603dc7d3c7cbee83587b95c91ab4fe71e8a3b5702` (dimension fourteen).

## Related work

- Planar domains for this problem: Colbrook, Sadeghi and Stepaniants
  ([arXiv:2608.08953](https://arxiv.org/abs/2608.08953)), a bounded simply connected non-disc with analytic boundary
  and a sign-changing eigenfunction; and the strictly convex planar domain in the paper of this repository
  ([arXiv:2609.35419](https://arxiv.org/abs/2609.35419), version 2), whose Question 2 asks whether Berenstein's
  conjecture fails in every dimension n ≥ 3, and within the convex class in those dimensions.
- Dai, Sun, Wei and Zhang ([arXiv:2511.19819](https://arxiv.org/abs/2511.19819)): a bounded uniformly convex
  planar domain with connected C^{2,ε} boundary is a disc if the problem has a non-trivial solution corresponding
  to a large eigenvalue.
- Related sign-changing problems with other nonlinearities or other classes of domains: Ruiz (J. Eur. Math. Soc.,
  2025), Dai and Zhang ([arXiv:2304.04525](https://arxiv.org/abs/2304.04525)), Minlend
  ([arXiv:2307.07784](https://arxiv.org/abs/2307.07784)).

## Files

| Path | Description |
|---|---|
| `paper.pdf`, `source/` | The paper and its LaTeX source (`main.tex`, `references.bib`, `main.bbl`). |
| `verification/r4/` | Dimension four: single-file verifier, its programs and certificate data, README, checksums. |
| `verification/r14/` | Dimension fourteen: the same. |
| `SHA256SUMS` | Checksums of all files in this directory. |

## Citation

```bibtex
@misc{Guo2026Berenstein,
  author = {Guo, Jizhou},
  title = {Counterexamples to {B}erenstein's conjecture in dimensions four and fourteen},
  year = {2026},
  month = oct,
  howpublished = {\url{https://github.com/aster2024/schiffer-pompeiu-counterexamples/tree/main/berenstein-dimensions-4-14}},
  note = {Preprint, 10 October 2026}
}
```
