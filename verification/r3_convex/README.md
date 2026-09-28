# Convex three-dimensional certificate (AX4)

This directory is the local computer-assisted proof package for the convex, non-ball, axisymmetric example in R^3. It was independently audited (mathematics, code and numerics). It is separate from `../r3/`, which proves an additional nonconvex example.

The exact dyadic centre is `candidate.json` (SHA-256 `547dd683610fbe00753ee1af34ed19f6ff5d1bde76b5b7dacdba2031131bcb3b`). The unmodified single-file verifier is `verify_r3_convex.py` (SHA-256 `372815efbd8e8ed337d0d577ebc9aec3f1212c5ea6ef26acaf206686f9a418d2`). `VERIFICATION_RECEIPT.json` is its completed 17-stage replay receipt; only its machine-specific `replay_workdir` string was replaced for distribution. `certificate_r3_convex.json` records the final numerical existence and geometry gates. `geometry_projected_radius.json` is an earlier conditional geometry calculation on the actual projected radius; the final receipt separately certifies existence. `SCALE3_OBSTRUCTION.json` explains why the scale 3b contraction gate fails at shape column j=183.

From this directory, using an empty scratch directory:

```bash
python3 -O verify_r3_convex.py --workdir <empty dir> --output receipt.json --jobs 5
```

Requires 64-bit Linux, Python 3, python-flint 0.9.0, NumPy, SciPy and gcc. The recorded full replay took 6515.9 seconds, about 1.8 hours, with five single-thread workers. A second, independent from-scratch run of the same unmodified verifier in a separate fresh directory also printed `PROVED`; its receipt `VERIFICATION_RECEIPT_replay2.json` agrees with `VERIFICATION_RECEIPT.json` in every mathematical field (only the run time, the stage completion order, the per-run matrix hash and the work directory differ). The verifier refuses nonempty workdirs and does not accept old PASS receipts. Check file integrity with `sha256sum -c SHA256SUMS`.
