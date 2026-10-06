# Ancillary files

These files accompany the planar part of "Convex counterexamples to the Schiffer and Pompeiu conjectures in dimensions two to eighteen". The subsection "The plane" of the computer-assisted verification section describes them.

## Verifier and inputs

- `verify_planar_convex_v2.py`: the single-file verifier.  It embeds the nine
  stage sources and the two decimal inputs listed below and reads no
  precomputed numerical output.  Its SHA-256 is printed in the planar replay subsection of the
  paper.  Replay command (from the directory that contains `anc/`):

      nice -n 19 env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
        MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
        PYTHONDONTWRITEBYTECODE=1 PYTHONPATH= \
        python -O anc/r2/verify_planar_convex_v2.py \
        --workdir ./replay --jobs 2 \
        --output ./REPLAY_RECEIPT.json

  The working directory must be new or empty and the output file must not
  exist.  `--jobs` accepts 1 to 5.
- `center_j40_frozen.json`: the forty-mode centre.  The array `p` stores the
  coefficients c_j of psi_0(w) = sum_j c_j w^(182 j + 1); the derivative
  coefficient is (182 j + 1) c_j.  The field coefficient is a[j]/scales[j],
  with exact positive decimal divisors.
- `reference.json`: the eight-mode reference shape used for the spectral
  inverse.
- `spectral_interval.py`, `coverage_checks.py`, `boundary_interval.py`,
  `boundary_newton.py`, `analytic_cap.py`, `final_gates.py`,
  `identity_checks.py`, `scalar_second_iv.py`, `negative_controls.py`: the
  stage sources, byte-identical to the copies embedded in the verifier.
- `standalone_driver.py`, `build_standalone.py`, `BUILD_MANIFEST.json`: the
  driver template, the script that assembles the single-file verifier from the
  sources and inputs, and the hashes of the embedded files.
- `REQUIREMENTS.txt`: software versions of the recorded run.

## Recorded results

- `VERIFICATION_RECEIPT.json`: the receipt of the recorded run (five workers).
  It contains the enclosures of the table of planar enclosed constants, the stage list with
  timings, the software versions, the rejected perturbations and the hashes of
  the replay artifacts.  A successful run certifies the numerical inequalities
  used in the paper; the analytic arguments of the paper are not
  machine-checked.  The receipt field `analytic_proof` and some comments in
  the sources refer to `PROOF_PACKAGE.md`, the written proof from which the
  paper was prepared; the three proof sections of Part 1 (The plane) replace that file, which is
  not distributed.  A few source comments also carry labels from the
  development of the certificate.  The sources are distributed unchanged
  because the verifier embeds them byte for byte and its hash is fixed.
- `ROBUSTNESS_CHECK.json`: the radii inequalities of the planar contraction lemma re-evaluated
  with the inverse bound A_* multiplied by 10^8.
- `exact_constants.py`: exact rational evaluation (Python `fractions`) of the
  algebraic constants of the two inputs that are quoted in the planar spectral, inverse and existence subsections of the
  paper (v_*, gamma, P, P_1, V_*, d_h, beta, d_*, j_*, delta_q, S and the
  geometric margins).  Run `python exact_constants.py` in this directory.

## Non-interval diagnostics

`diagnostics/` contains binary64 computations for the forty-mode centre that
are quoted in the introduction and the rigidity section of the paper for description and for the
comparison with rigidity theorems.  They are not interval enclosures and are
not used in the existence proof.

- `geometry.py`: radius range, curvature range, convexity
  factor, width ratio, relative Hausdorff distance to centred discs, and the
  two numbers of the comparison with Kobayashi's theorem.
- `dswz.py`: the two-point amplitude defect, the function T(s),
  the stationary-phase remainder R_D at two frequencies, and the size of the
  fifth derivative of the support function, for the comparison with
  Dai-Sun-Wei-Zhang.

Run each script from any directory, for example
`python diagnostics/geometry.py`; it reads `center_j40_frozen.json` relative
to its own location. The figure of the paper shows the same centre.

## Checksums

`SHA256SUMS` lists the SHA-256 of every other file in this directory and in
`diagnostics/`.  It was regenerated for the present selection of files; the
hashes of the verifier, the inputs, the stage sources and the receipt are
unchanged from the recorded run and agree with `BUILD_MANIFEST.json` and with
the receipt.  Check with `sha256sum -c SHA256SUMS` in this directory.
