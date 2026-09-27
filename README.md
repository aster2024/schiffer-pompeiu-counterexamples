# A counterexample to the Schiffer and Pompeiu conjectures in dimension four

This repository contains the paper and the computer-assisted verification for the following result.

**Theorem.** There exist:
- a bounded domain Ω ⊂ ℝ⁴ that is strictly star-shaped with respect to the origin, invariant under O(3)×ℤ₂, and **not a ball**, whose boundary is a real-analytic hypersurface diffeomorphic to S³;
- a nonconstant function u, real analytic on a neighbourhood of Ω̄;

such that

    Δu + u = 0 in Ω,   u = 1 and ∇u = 0 on ∂Ω.

**Consequences.**
- Schiffer's conjecture fails in ℝ⁴.
- The Fourier transform of the indicator of Ω vanishes on the unit sphere, so Ω fails the Pompeiu property and the Pompeiu conjecture fails in ℝ⁴.
- The domain Ω is moreover **convex**.

Planar counterexamples were found in 2026 by Colbrook–Stepaniants (arXiv:2608.01579) and Cao-Labora–de Dios Pont (arXiv:2608.05114). This work treats dimension four. Dimensions three and five are discussed only as numerical evidence in the paper; they are **not** proved here.

## Contents

| Path | Description |
|---|---|
| `paper/main.pdf`, `paper/main.tex` | The paper: statement, reduction to a planar problem, radii-polynomial proof and verification details. |
| `verification/verify_r4.py` | Standalone interval-arithmetic (Arb) verifier of the main theorem. The approximate solution is embedded as exact dyadic rationals. Its docstring refers to "FINAL.md"; the corresponding mathematics is the paper. |
| `verification/verify_convex.py` | Interval-arithmetic verification of the convexity statement. |
| `verification/data/centre.json` | Provenance copy of the embedded centre data. |
| `verification/certificates/` | Certificates produced by the reference runs, at 128 and 192 bits. |
| `SHA256SUMS` | Checksums of all verification files. |

## Reproducing the verification

The only requirements are Python ≥ 3.10 and `python-flint` (tested with 0.9.0; install with `pip install python-flint==0.9.0`). Each main run takes about 3 minutes on one core.

```bash
mkdir -p out
OMP_NUM_THREADS=1 python3 verification/verify_r4.py --stage all --bits 128 --jstar 1001 --output out/certificate_r4_128.json
# last line of output: PROVED
OMP_NUM_THREADS=1 python3 verification/verify_convex.py --bits 128 --main-certificate out/certificate_r4_128.json --output out/convex_128.json
# first line of output: PROVED_CONVEX [0.41407413991... +/- ...]
sha256sum out/*.json   # compare with SHA256SUMS
```

Use `--bits 192` for the higher-precision run. The certificate JSON files are byte-for-byte reproducible.

## Status

- This is a computer-assisted proof.
- The verification was re-run independently from this repository layout and reproduced the certificates byte for byte.
- The work has **not yet been peer reviewed**. Comments and independent checks are very welcome.

## AI statement

This work was carried out with the assistance of AI models (Claude Opus 5.5 and GPT-6 Astra). The author takes responsibility for the content.

## Author

Jizhou Guo, Dots Studio, Rednote. Contact: mitsuha2021b@gmail.com, sjtu18640985163@sjtu.edu.cn

## License

- Code in `verification/`: MIT License (see `LICENSE`).
- Paper: CC BY 4.0.
