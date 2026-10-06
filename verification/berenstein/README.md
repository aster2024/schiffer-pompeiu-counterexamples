# Convex planar Berenstein certificate

The proof is Part “A convex planar domain for Berenstein's problem” of the
paper. The exact domain is specified by the zero enclosed about the frozen
clamped centre, rather than by the finite conformal polynomial alone.

`verify_berenstein_planar_convex_v3.py` embeds the exact-decimal centre and
ten stage modules, with their SHA-256 hashes. It embeds no computed
spectral matrix, boundary enclosure or existence receipt. Its SHA-256 is
`375ea2e68225d3419309d5359c7d5d2a589c07bb0985e71f2cd0ac9bc3fc95fb`.
The separate `centre_frozen.json` exposes the exact inputs.

Install the versions in `REQUIREMENTS.txt` in an environment visible to
isolated Python. From this directory run, with `replay` absent:

```bash
nice -n 19 env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 \
  python -I -O verify_berenstein_planar_convex_v3.py --workdir replay --jobs 4
```

The final success lines begin:

```text
NUMERICAL_GATES_PASS <absolute-replay-directory>/VERIFICATION_RECEIPT.json
RADIUS [1.0000000000000000000000000000000000000000000000000000000000000000e-36 +/- 3e-105] CONTRACTION
```

The second line continues with the enclosed contraction bound. Accept only
exit zero and a fresh receipt with `status=NUMERICAL_GATES_PASS`, matching
the verifier and centre hashes. The verifier recomputes both complete Robin
window enumerations, the spectral and boundary enclosures, the scalar chain
and whole-ball geometry. The parent rechecks the quoted bounds listed in
`parent_gates` and the actual matrix dimensions. Finite interval identity
checks test translations, the residual-corrected material equations and both
clamping traces. Independent Taylor division binds the Robin coefficient to the
physical centre, and Laurent trace products check the assembled boundary matrix.

The second library recomputes 56 scalar quantities in 400-bit `mpmath.iv`
arithmetic from six Arb spectral primitives and eight boundary moments.
It does not reconstruct the special functions, spectral matrices or boundary
moment enclosures. The 403 semantic controls and 31 identity negative controls
must pass; unresolved comparisons, nonfinite values and child failures reject
the run. Output hashes are read back before the receipt is published.

`VERIFICATION_RECEIPT.json` records an isolated replay: Python
3.10.19, python-flint 0.9.0, NumPy 1.26.4, SciPy 1.15.2 and mpmath 1.3.0,
four single-threaded workers, `nice 19`, 173.04210948944092 seconds.
Only its machine-specific Python executable path is replaced by `python`;
all numerical fields and emitted status names are preserved. Fresh run times
and run-specific hashes may differ. Compare the numerical enclosures with
the paper's outward-rounded bounds.

The verifier certifies the numerical hypotheses of the theorem. The analytic
lemmas are proved in the paper; `full_proof=false` records that they are not
machine-checked. The finite identity stage is a regression check of those
identities. Floating-point diagnostic comparisons with earlier work are in
`../diagnostics/berenstein/` and are not inputs to this certificate.

The field `proof_scope` of the receipt and some strings in the verifier refer
to `PROOF_PACKAGE.md`, the written proof from which the paper was prepared; the
paper part named above replaces that file, which is not distributed. A few
strings in the verifier and receipt also carry labels from the development of
the certificate. The verifier is distributed unchanged because it embeds its
sources byte for byte and its hash is fixed.

Before creating replay outputs, run `sha256sum -c SHA256SUMS` here.
