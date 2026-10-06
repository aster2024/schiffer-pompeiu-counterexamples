# Exact inverse in general dimension

This distribution contains the single-file verifier `verify_axisym_exact_v2.py`,
eleven losslessly compressed exact-decimal centres and the overall final
isolated-run receipts for dimensions 5, 9, 11, 12, 13, 15, 16, 17, 18, 20, 21.
The verifier SHA-256 is
`b33fd6c7e9d72a283c9f752ec55b8e1d4dec48b132872759ac4a32553d1e3940`.
It embeds the dimension parameters, stage sources and centres; it embeds no
computed spectral matrix, boundary enclosure or existence receipt.
`REQUIREMENTS.txt` records the versions of the final runs (Python 3.10.19).

`frozen_centers/nN.json.gz` returns the original JSON input byte for byte
under `gzip -dc`; its keys are `n`, `rho`, `c`, `a`, `scales`, `status`.
No decimal is shortened. The driver uses its own hash-bound embedded copies.

With these dependencies installed, choose an absent or empty output directory:

```bash
nice -n 19 env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  python -E -O verify_axisym_exact_v2.py --dimension 12 --workdir replay12
```

Accept only exit zero, the complete final line `PROVED dimensions=12`, and
status `PROVED` in both `replay12/VERIFICATION_RECEIPT.json` and
`replay12/n12/VERIFICATION_RECEIPT.json`. `NOT_PROVED` and `NOT_ALL_PROVED`
are failure statuses. `--dimension all` executes the eleven dimensions
sequentially; its final line lists precisely those dimensions.
The spectral module pins to the lowest allowed CPU; an optional `taskset`
restriction selects one allowed CPU. All recorded runs use one Arb/BLAS
thread and `nice 19`. Spectral precision and radial quadrature are 768 bits
and 384 nodes, except for dimension 16 (1536 bits and 640 nodes).
The dimension-16 run needs several GB of memory.

The retained `receipts/nN/VERIFICATION_RECEIPT.json` files are byte-for-byte
copies of the receipts of the final isolated runs (dimension 21 has one
completed isolated replay of the verifier). Each receipt binds the dimension,
verifier, embedded payload and per-dimension receipt hash. A fresh replay
regenerates the extracted sources, stage outputs, identity receipt, both
final receipts and hash readbacks. Check each fresh receipt against its
own artifact hashes. Run times and run-specific receipt hashes can change;
compare the mathematical enclosures with the outward-rounded paper tables.

| Dimension | Chain seconds | Complete process seconds |
|---|---:|---:|
| 5 | 528.104 | 535.312 |
| 9 | 637.539 | 645.653 |
| 11 | 370.256 | 378.099 |
| 12 | 286.531 | 294.107 |
| 13 | 581.188 | 588.951 |
| 15 | 209.550 | 217.013 |
| 16 | 1243.897 | 1252.022 |
| 17 | 170.385 | 177.179 |
| 18 | 150.032 | 157.956 |
| 20 | 320.510 | 328.003 |
| 21 | 426.939 | 434.442 |

The analytic proof is Part “The exact inverse in general dimension”. The
verifier recomputes its numerical hypotheses. The name `PROOF_PACKAGE.md`
in the `scope` field of a fresh per-dimension receipt refers to the written
argument now given in that Part; that file is not distributed.
The unchanged symbolic stage
runs 125 connection/lift examples over dimensions 5, 9, 11, 12, 13 and the
symbolic-dimension material identity. There are 52 semantic failure controls.
The second `mpmath.iv` computation recombines nine final scalar expressions
from Arb inputs; it does not independently rebuild special functions,
spectral matrices, boundary samples, shape-norm constants, norm-transfer
constants or fixed-field primitives.

Before creating replay outputs, run `sha256sum -c SHA256SUMS` here.
