# Convex non-ball Schiffer domains in infinitely many dimensions

Jizhou Guo — mitsuha2021b@gmail.com — 8 October 2026

[**Paper (PDF, 114 pages)**](paper.pdf) · [LaTeX source](source/) · Archive: [doi:10.5281/zenodo.23238160](https://doi.org/10.5281/zenodo.23238160)

## Results

A *Schiffer domain* is a bounded domain Ω ⊂ ℝⁿ carrying a nonconstant solution of

    Δu + u = 0 in Ω,   u = 1 and ∇u = 0 on ∂Ω.

Let J_ν be the Bessel function of the first kind. The *primary crossing orders* p̂_m are the real solutions, for all
large indices m, of

    J_p( 2·sqrt((p+1)(p+2)) ) = 0.

**Theorem A (bounded detuning).** For every sufficiently large integer or half-integer p with
|p − p̂_m| < 31/(10 p̂_m) for some m, there are an integer a ∈ [p/2, 3p/4] and a bounded strictly convex domain
Ω ⊂ ℝ^(2p), not a ball, with real-analytic boundary and invariant under O(a) × O(2p − a), which is a Schiffer domain.
There are infinitely many such p among the integers and among the half-integers. For every modulus L, every residue
class modulo L contains infinitely many dimensions n = 2p with the same conclusion (for even L a residue class
fixes the parity of n; for odd L each class contains infinitely many dimensions of both parities). In each of these
families, all but O(n^(2/3+η)) integer splits a in a fixed interior range give such domains, pairwise non-similar
for distinct splits.

**Theorem B (logarithmic window).** For every sufficiently large integer or half-integer p with
|p − p̂_m| ≤ log p/(8p) for some m, there are an integer split 0 < a < p and a bounded strictly convex Schiffer domain
Ω ⊂ ℝ^(2p), not a ball, with real-analytic boundary and invariant under O(a) × O(2p − a), with the same count of
pairwise non-similar domains. This family contains
sequences of both parities along which the detuning is unbounded.

**Consequences.** By Green's identity the Fourier transform of the indicator of each of these domains vanishes on
the unit sphere, so the domains fail the Pompeiu property. The Schiffer and Pompeiu conjectures therefore fail
within the class of convex domains in infinitely many dimensions, both odd and even.

## Method

The domains are invariant under O(a) × O(2p − a). The block split a is the parameter, at fixed ambient dimension
and fixed radial spectrum of the ball. The construction starts from a quartic angular seed at a near-resonance of
the ball. An exact spectral-flow argument in the split parameter gives a Dirichlet gap at an integer split;
finite-order centres solve the interior equation exactly with small boundary defects; the linearisation is inverted
in one graded Jacobi coefficient space, through bounded-cardinality angular clusters (Theorem A) or a mixed
conditioning estimate (Theorem B); a Newton iteration closes the construction, and strict convexity, analyticity
and non-similarity are read off from the centre. The proofs are analytic. The thresholds in p are not explicit, and
each family of dimensions has density zero.

Convex counterexamples in the plane, in every dimension from three to eighteen and in dimensions twenty and
twenty-one, with computer-assisted proofs, are in the first paper of this repository
(arXiv:[2609.35419](https://arxiv.org/abs/2609.35419)).

## Files

| Path | Content |
|---|---|
| `paper.pdf` | The paper |
| `source/main.tex`, `source/inf_*.tex` | LaTeX source; the bibliography is in `source/inf_references.tex` |
| `SHA256SUMS` | Checksums of the files above |

To compile: `cd source && tectonic main.tex`, or run `pdflatex main.tex` three times.

The paper is also on ResearchGate: <https://www.researchgate.net/publication/415393919_Convex_non-ball_Schiffer_domains_in_infinitely_many_dimensions>.

## Citation

```bibtex
@misc{Guo2026SchifferInfinitelyMany,
  author = {Guo, Jizhou},
  title = {Convex non-ball {S}chiffer domains in infinitely many dimensions},
  year = {2026},
  month = oct,
  howpublished = {\url{https://github.com/aster2024/schiffer-pompeiu-counterexamples/tree/main/infinitely-many-dimensions}},
  note = {Preprint, 8 October 2026. Archived in Zenodo, version 5.0},
  doi = {10.5281/zenodo.23238160}
}
```
