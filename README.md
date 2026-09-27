# Counterexamples to the Schiffer and Pompeiu conjectures in dimensions four and six

[![DOI](https://zenodo.org/badge/1391007796.svg)](https://doi.org/10.5281/zenodo.22999505)

This repository contains the paper and the computer-assisted verification for the following result.

**Theorem.** For n = 4 and n = 6 there exist a bounded domain Ω ⊂ ℝⁿ and a nonconstant function u such that

    Δu + u = 0 in Ω,   u = 1 and ∇u = 0 on ∂Ω.

The domain Ω has the following properties:
- it is convex and strictly star-shaped, and it is **not a ball**;
- its boundary is a real-analytic hypersurface diffeomorphic to Sⁿ⁻¹.

The function u is real analytic on a neighbourhood of Ω̄.

The symmetries are:
- n = 4: invariant under O(3)×ℤ₂;
- n = 6: invariant under O(3)×O(3), including the exchange of the two factors.

**Consequences.** Schiffer's conjecture fails in ℝ⁴ and ℝ⁶. The Fourier transform of the indicator of Ω vanishes on the unit sphere, so Ω fails the Pompeiu property, and the Pompeiu conjecture fails in these dimensions as well.

**Earlier and other dimensions.**
- Planar counterexamples were found in 2026 by Colbrook–Stepaniants (arXiv:2608.01579) and Cao-Labora–de Dios Pont (arXiv:2608.05114).
- The four- and six-dimensional problems reduce to planar problems through the identities Δ₂(y·u) = y·Δ₄u and Δ₂(y₁y₂·u) = y₁y₂·Δ₆u.
- Dimensions three and five are discussed in the paper only as numerical evidence and are **not** proved here.

## Contents

| Path | Description |
|---|---|
| `paper/main.pdf`, `paper/main.tex` | The paper (42 pages). |
| `verification/verify_r4.py`, `verification/verify_convex.py` | Arb (interval-arithmetic) verifiers for n = 4: the main theorem and convexity. The centre is embedded as exact dyadic rationals. The docstring of `verify_r4.py` refers to "FINAL.md"; the corresponding mathematics is the paper. |
| `verification/data/`, `verification/certificates/` | Centre data and reference certificates for n = 4. |
| `verification/r6/` | Everything for n = 6: verifier, stage scripts, centre, approximate inverse, seven interval receipts, certificates, identity checks and convexity check. See `verification/r6/RELEASE_README.md`. |
| `SHA256SUMS`, `verification/r6/SHA256SUMS` | Checksums. |

## Reproducing the verification

The requirements are Python ≥ 3.10 and `python-flint` 0.9.0 (`pip install python-flint==0.9.0`); n = 6 also needs `numpy`. Runs are single-core.

**n = 4** (about 3 minutes):
```bash
mkdir -p out
OMP_NUM_THREADS=1 python3 verification/verify_r4.py --stage all --bits 128 --jstar 1001 --output out/certificate_r4_128.json      # last line: PROVED
OMP_NUM_THREADS=1 python3 verification/verify_convex.py --bits 128 --main-certificate out/certificate_r4_128.json --output out/convex_128.json   # PROVED_CONVEX
```

**n = 6** (about 3 minutes; recomputes all seven receipts and compares them byte for byte):
```bash
cd verification/r6
PYTHON=python3 bash reproduce_r6.sh      # ends with: REPRODUCTION COMPLETE
sha256sum -c SHA256SUMS
```

All certificate files are byte-for-byte reproducible.

## Status

- These are computer-assisted proofs.
- Both certificates were re-run independently and reproduced byte for byte.
- Each dimension was examined by separate internal mathematical and computational audits.
- The work has **not yet been peer reviewed**. Comments and independent checks are very welcome.

## AI statement

This work was carried out with the assistance of AI models (Claude Opus 5.5, GPT-6 Astra and GPT-6 Sol). The author takes responsibility for the content.

## Author

Jizhou Guo, Dots Studio, Rednote. Contact: mitsuha2021b@gmail.com, sjtu18640985163@sjtu.edu.cn

## License

- Code in `verification/`: MIT License (see `LICENSE`).
- Paper: CC BY 4.0.
