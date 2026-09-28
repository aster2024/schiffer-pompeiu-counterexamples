# R³: computer-assisted proof package

Computer-assisted proof that there is a bounded, strictly star-shaped,
axisymmetric domain `Ω ⊂ R³`, not convex and not a ball, with real-analytic
boundary diffeomorphic to `S²`, carrying a nonconstant solution of

```text
Δu + u = 0 in Ω,    u = 1,  ∇u = 0 on ∂Ω .
```

Hence `Ω` fails the Pompeiu property. The mathematics is in the R³ section of
the paper; this directory contains the verifier and its certificate.

## Files

| File | Content |
|---|---|
| `verify_r3.py` | Self-contained verifier, 286905 bytes, SHA-256 `c2802ee1202ac2dd899f01c4e5794f1a7e5f78872f53becd1728f23e19180cdb`. It embeds (zlib + base85, SHA-256 `0ba7a364b8bed66c4e50531f871762b8a185968a9bb38b59350c4cc9ad749a4d`) all stage scripts, the C accelerator source and the frozen dyadic centre `refined_M180_step1.json` (SHA-256 `1c2e959714280ca4bbd2e60c5a7bc378486f5c96205d858e0de0c561ef26a344`), 26 files in total, each with its own SHA-256. |
| `certificate_r3.json` | Certificate written by a full replay of `verify_r3.py` (see below). |
| `certificate_r3_replay2.json` | Certificate of a second, independent from-scratch replay in a separate directory, with the job limit raised from four to eight (see below); all certified values agree with `certificate_r3.json`. |
| `check_nonconvex_r3.py` | Exact check, in rational arithmetic, of the non-convexity at the poles (see below). Standard library only; it reads the centre from the payload of `verify_r3.py` and does not modify it. |
| `SHA256SUMS` | Checksums of the files in this directory. |

The verifier is released exactly as it was reviewed; it has not been modified.

## Requirements

* 64-bit Linux, Python 3 (the recorded replay used Python 3.10.19), NumPy,
  SciPy, and **python-flint 0.9.0 installed from the PyPI binary wheel**: the
  verifier compiles its C interval accelerator with `gcc` (`-O2 -fno-fast-math
  -ffp-contract=off`) and links it against the FLINT shared library bundled in
  that wheel (`python_flint.libs/libflint-*.so.*`). It checks the
  python-flint version and that this library is unique.
* `gcc`.
* At most four single-thread worker processes (all at nice 19;
  `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`, `MKL_NUM_THREADS` are forced to 1),
  and about 1 GB of free disk space in the work directory.

## Command and expected output

```bash
sha256sum -c SHA256SUMS
python3 -I -O verify_r3.py --workdir <new-or-empty-dir> --jobs 4 \
    --max-seconds 17000 --output cert.json
```

The work directory must be new or empty (existing receipts are never
accepted); without `--workdir` a fresh temporary directory is created. With
`--jobs 4` the recorded replay took 4008 s (about 67 minutes) of wall-clock
time; a second, independent replay with a job limit of eight (see below) took
3578 s. With fewer workers (`--jobs 1` to `3`) the run takes longer, and
`--max-seconds` (default 28800) must be raised accordingly.

The verifier decompresses its payload, checks the payload hash and every
per-file hash, writes the files into the work directory, compiles the
accelerator, and recomputes all 16 stages (`matrix`, `principal`, `taylor`,
`residual`, `finite_matrix`, `g_near_0/20/40/60`, `jac_g`, `axis`,
`g_envelopes`, `shape_far`, `jac_shape`, `shape_all`, `assemble`); it then
re-checks every file hash and requires the assembled status
`ALL_EXISTENCE_GATES_PASSED`. It prints progress lines (`START`, `DONE`,
`RUNNING`) and finally

```text
PROVED
CERTIFICATE <absolute path of cert.json>
```

with exit code 0. Any missing coverage, NaN or infinite interval, hash
mismatch, failed stage or failed inequality prints `NOT_CERTIFIED <reason>`
and exits with code 2.

## Non-convexity at the poles

```bash
python3 check_nonconvex_r3.py
```

checks the SHA-256 of `verify_r3.py`, of its payload and of the centre file,
reads the 91 dyadic shape coefficients `c°_j`, and verifies in exact rational
arithmetic that `ψ°'(1) = Σ j c°_j ∈ (42.81903, 42.81904)` and
`ψ°'(1) + ψ°''(1) = Σ j² c°_j ∈ (−7.758015, −7.758014)`. For every `c` in the
existence ball `Σ j (21/20)^j |c_j − c°_j| ≤ r/(3b)`, `r = 3·10⁻⁷`, it follows
that `ψ_c'(1) > 0`, `ψ_c'(1) + ψ_c''(1) < −7.758` and
`1 + ψ_c''(1)/ψ_c'(1) < −0.181`, so that both principal curvatures of `∂Ω` at
the two poles are negative (paper, Theorem 1.1(iv)). It prints
`NONCONVEX_AT_POLES_VERIFIED` and runs in a fraction of a second.

## The certificate

`certificate_r3.json` is the certificate written by the replay (status
`PROVED`, `replay_seconds` 4008.36, four workers). The only edit is that the
field `replay_workdir`, which records the scratch directory passed to
`--workdir`, has been replaced by the placeholder `<workdir>`. A new run
writes the same enclosures, with its own `replay_workdir` and
`replay_seconds`.

`certificate_r3_replay2.json` is the certificate of a second replay, run from
scratch in a separate, new work directory (status `PROVED`, `replay_seconds`
3577.54). It used a copy of `verify_r3.py` whose only change is line 3381,
the scheduling limit on the number of worker processes:

```text
<     args=ap.parse_args();need(1<=args.jobs<=4,'jobs must be between one and four')
>     args=ap.parse_args();need(1<=args.jobs<=8,'jobs must be between one and eight')
```

The patched file has 286906 bytes and SHA-256
`829076e8b0755aa0a6e63474395e207f8ed6c54a47793fdab06b12e0d3ad2c8d`, which is the
value of `verifier_sha256` in this certificate; its embedded payload (SHA-256
`0ba7a364…749a4d`) is unchanged. As for `certificate_r3.json`, the field
`replay_workdir` has been replaced by `<workdir>`. Compared field by field with
`certificate_r3.json`, all mathematical fields agree exactly (every enclosure,
the radius, the coverage counts, the status, and the checksums of the centre
and of the payload); the only differences are `replay_seconds`,
`verifier_sha256` and the order in which the sixteen stages completed
(`completed_stages`, the same set).

Certified bounds (outward-rounded here; the certificate contains the Arb
enclosures): approximate-inverse norm `‖A‖ ≤ 6000`; `Y < 2.1003·10⁻⁸`;
`Z < 0.858929` (finite field columns `< 0.720601`, 80 exceptional near field
columns `< 0.327488`, far field classes `< 0.858929`, explicit shape columns
`< 0.330128`, far shape columns `< 0.803677`); `c₂ < 26.60384`,
`c₃ < 1.04114·10⁻³`, `c₄ < 6.16335·10⁻⁹`; radius `r = 3·10⁻⁷`; radii
polynomial `< −4.8199·10⁻⁹`; contraction bound `< 0.95921`; derivative margin
`> 35.4391`; star-shapedness margin `> 0.832086`; non-ball margin
`|c₁₃| > 0.231712`. Coverage: 5824 finite field columns, 80 exceptional field
columns, 91 radial and 1 angular infinite field classes, 260 explicit shape
columns, and all shape columns from index 521 on by the far-shape bound.

## Review notes

An internal review of this package (mathematics, code, numerics) found no
gap. It recorded the following cosmetic points; none of them affects the
certificate.

Written proof (taken into account in the paper):

1. The τ-norm of `q`, `q_x`, `q_y/y`, grouped by degree, is taken equal to the
   raw Chebyshev-weighted sum. This uses that the Legendre coefficients of the
   Chebyshev polynomials `U_n` and of the Gegenbauer polynomials `C_m^(2)` are
   nonnegative. This holds because their generating functions are the second
   and fourth powers of the Legendre generating function
   `(1 − 2xt + t²)^(−1/2)`, and products of Legendre polynomials have
   nonnegative Legendre coefficients (Gaunt); the reviewer also checked it
   exactly for `n ≤ 200`.
2. Proof of `h_{ℓ+1} ≤ R h_ℓ`: write `R h_ℓ − h_{ℓ+1} = Σ_j δ_j R^{n_j}`,
   where only the lowest-frequency coefficient `δ_j` is negative and
   `Σ_j δ_j = 0` (the identity at `R = 1`); hence the sum is `≥ 0` for
   `R ≥ 1`. (Checked numerically for `ℓ ≤ 2600`.)
3. Quadrature: the replay uses 697 radial Gauss–Legendre nodes (exact to
   degree 1393) and 542 angular nodes (exact to degree 1083). The integrands
   have radial degree at most 850 + 540 + 2 = 1392 and angular degree at most
   1080, so these coarse degree counts suffice. (An earlier text referred to
   696 radial nodes and to a finer degree count.)
4. Analytic continuation of `u` across `∂Ω`: cite the analytic boundary
   regularity theorem of Morrey–Nirenberg (analytic boundary, constant
   Dirichlet data) rather than a brief Cauchy–Kovalevskaya argument. The main
   theorem only needs `u ∈ C¹(Ω̄)`, analytic in `Ω`.
5. Typesetting of `{\rm even}` in norm definitions of the working notes.

Code (the verifier is hash-frozen, so these were documented, not changed):

6. `coeff3_native.c`: the return value of `arb_set_str` is not checked in
   `vec_linear` and `projection_prepare`. All strings parse correctly: the
   native weights agree with the Python reference for every `ℓ ≤ 698`.
7. Python `max()` is applied to non-exact Arb balls in
   `certify_envelopes.py` (e.g. `max(arb(0), bb-aa)`), `certify_principal.py`,
   `principal.py`, `inverse_apply.py` and `certify_shape_far.py`. When a ball
   straddles the other argument this can return a value too small by at most
   the ball radius (about `10⁻¹¹⁵` here), and a NaN would be dropped. In this
   run all inputs are finite and the deciding comparisons are far from ties.
8. `assemble_certificate.py`, `ball()`: requires the upper endpoint to be
   `≥ 0` but not the lower one. Harmless, because every gated quantity first
   passes a check rejecting non-finite values.
9. The stage `jac_shape` (`interval_M180_shape.json`, interval finite shape
   Jacobian) is computed and coverage-checked but not used in any inequality;
   the finite shape columns are certified by the full defect computation. It
   is a cross-check only.
10. The Taylor data for the term `H` of the far-shape bound come from the
    native code and were not recomputed independently. Its contribution is
    about `6.9·10⁻⁴/3`; even a hundredfold underestimate would keep the far
    shape bound below `0.862 < 1`. The Taylor data for `W`, `D_x` and `Z` were
    recomputed independently and agree with the certificate.

Numerics and geometry:

11. The domain is not convex: `min Re(1 + wψ''/ψ') = −0.181181…` on the unit
    circle, attained at the two poles (the quantity is negative for
    `θ ∈ [0, 0.0553] ∪ [π − 0.0553, π]`); at the poles both principal
    curvatures of `∂Ω` are about `−0.00423`. The verifier itself checks only
    strict star-shapedness; the negativity at the poles, on the whole
    existence ball, is proved in the paper from the exact values checked by
    `check_nonconvex_r3.py` (added after the review; the verifier is
    unchanged).
12. Independent recomputations by the reviewers: field-tail envelope for the
    class `(46,65)` `0.858929`, angular class `0.832893`, output tail of the
    finite columns `0.718703`, far shape bound `0.80361` (certificate:
    `0.80368`); at the centre, the residual `F/b` is about `10⁻¹⁸`, and
    `max|Δu + u| ≈ 7·10⁻¹⁷` over 722 test points, with `max|u| ≈ 48.2`.
