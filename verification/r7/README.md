# Convex seven-dimensional certificate

Computer-assisted part of the paper's section **Dimension seven**:
a bounded strictly convex non-ball domain in R⁷ with real-analytic boundary
S⁶, O(6) × Z₂ symmetry, and a nonconstant solution of Δu + u = 0,
u = 1 and ∇u = 0 on the boundary.

The coefficient spaces use Υ(D) = (1 + D/10⁶)², R = 6/5,
σ = 25ρ² and finite block (M, S) = (210, 84), of dimension 9010.
The existence radius is 10⁻⁹ in the combined field/shape norm.

| File | Content |
|---|---|
| `verify_r7.py` | Unmodified self-contained verifier. SHA-256 `f8fd0ae17f4de35e7031b17ae61bde49e7d6c1e98705456e82cd1cca070989ed`. |
| `VERIFICATION_RECEIPT.json` | Original completed full-run receipt (`PROVED`, 6174.778263 seconds with at most eight single-threaded workers). Only the machine-specific `replay_workdir` string was replaced by `<workdir>`; timing, stage order and every mathematical field are preserved. SHA-256 `2cc620a9d4659c5283b17dea9b99cb7c8f6b84a22775142e2c9f77a4763565a1`. |
| `VERIFICATION_RECEIPT_replay2.json` | Second independent full run of the same unmodified verifier in a separate empty directory (`PROVED`, 10628.949744 seconds with six single-threaded workers). It agrees with the first receipt in every mathematical field; only the run time, stage order and working directory differ. The machine-specific `replay_workdir` string was replaced by `<workdir>`. SHA-256 `63238c99afe693a317173da105d20fe175005548b6822cab5b0945c5f9caff8b`. |
| `SHA256SUMS` | Checksums of the verifier, both receipts and this README. |

The compressed payload has SHA-256
`a9335ae79d93fe9977875e6b744e2d3a7176fe413182d3761d3234a693fdada0`.
It embeds the exact dyadic centre (SHA-256
`31de5cfdaebeb1788e50c9e8fbc57e470ac2252900d12549b308ee654f5a3cec`),
frozen class parameters, all stage programs and C kernels.
No development iteration or precomputed output is read as a proof input.
The payload is the authority for the centre numbers; development-file metadata
is not part of the construction.

## Reproduction

From this directory, using a new or empty scratch directory:

```bash
nice -n 19 env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  python3 -O verify_r7.py --workdir <empty-dir> --output receipt.json --jobs 2
```

Requirements: 64-bit Linux with x87 extended-precision long double,
Python 3, NumPy, **python-flint 0.9.0**, and `gcc`.
The recorded run took about 103 minutes with eight single-threaded workers;
`--jobs` accepts 1 through 8. The command above limits concurrency to two;
it is expected to take several hours, and its runtime has not been measured. The verifier refuses a nonempty work
directory. It prints `PROVED` and writes a receipt only after all gates pass;
a failure prints `NOT_CERTIFIED` and exits with status 2.

Check integrity with `sha256sum -c SHA256SUMS`.

## What is checked

The 28 stages rebuild the kernels and run exact self-tests, then recompute
from the exact centre: the full-support residual; all 9010 finite columns,
133 exceptional field columns and 245 explicit shape tail columns
(odd j = 213 through 701); the dyadic finite inverse and its rigorous
matrix-product defect; the weighted boundary Toeplitz inverse, finite fold,
profiles and Neumann remainders; 106 radial envelopes, 88 angular gap
envelopes and the uniform class ℓ ≥ 388; the all-index far shape class
j ≥ 703; the quartic Taylor coefficients, radii polynomial and contraction;
whole-existence-ball geometry; and exact coverage of every input.

The paper proves the all-index suffix-table extension from an even
Kmax ≥ max(2I + 2, M + 2, 2ic) and the geometric-switch inequalities.
Here (Kmax, I, M, ic) = (1202, 140, 210, 100); the existing integer gates
imply those hypotheses. Taylor angular degrees are not extension hypotheses.
Tables include both parities; entries of the unused parity only enlarge
suffix maxima. The verifier and its payload hashes are preserved.

## Certified values (outward roundings)

| Quantity | Bound |
|---|---|
| Y₀ | < 5.950563e-11 |
| Z₀ | < 0.898934481 |
| finite columns / exceptional field defects | < 0.058622927 / < 0.676824352 |
| explicit shape / far shape defects | < 0.651602561 / < 0.305605856 |
| ‖A‖ / A_t | < 342234.383037 / < 15.223634 |
| c₂ / c₃ / c₄ | < 0.745260675 / < 1.946871609e-6 / < 1.925257056e-12 |
| exact radius | 1/1000000000 |
| radii polynomial / contraction | < −4.13048e-11 / < 0.899444589 |
| Re ψ′ on radius 101/100 | > 0.993557571646 |
| Re(1 + wψ″/ψ′) on the unit disc | > 0.951916263887 |
| Re(wψ′/ψ) on the unit circle | > 0.994796343179 |
| |c₃| | > 0.000189913656777 |
| five transverse principal curvatures, unscaled | > 0.934750157145 |

Z₀ is the bound for the uniform angular class ℓ ≥ 388, s ≥ 1, given by
the envelope at (388, 1); all explicit defects are below 0.68. The geometric bounds are triangle-inequality lower bounds
uniform over the existence ball, rather than measured minima at the centre.
Physical principal curvatures are the unscaled values divided by ρ.
