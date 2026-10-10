# Dimension-fourteen verification

This directory verifies the G2-invariant dimension-fourteen theorem.
It asserts convexity of the ambient domain and strict convexity of its
Cartan section. The self-contained verifier embeds all required programs,
exact centres, dyadic inverse matrix and bounded column certificates.
The same files are supplied transparently here.

Requirements: Python 3, python-flint **0.9.0**, NumPy. The recorded replay
used Python 3.10.19 and NumPy 1.26.4 on Linux x86-64.

Run from this directory, with only one verifier running at a time:

~~~sh
sha256sum -c SHA256SUMS
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  nice -n 19 python3 -B verify_r14.py
~~~

Expected exit code: 0. Expected final line:

~~~text
PROVED
~~~

Both nonlinear coefficient files are regenerated in this ordinary run with
the weighted multiplier theta_1 = 1. Finite, near and far column coverage,
centre transfer and signed boundary factors are checked.

For complete serial regeneration and byte comparison of every bounded column:

~~~sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  nice -n 19 python3 -B verify_r14.py --recompute-all
~~~

Expected exit code: 0. Output includes ALL_RECEIPTS_RECOMPUTED;
its final line is PROVED. Complete regeneration of these same mathematical
programs and numerical data took 518.356 seconds.

Check the sharper numerical inequalities stated in the paper:

~~~sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  nice -n 19 python3 -B check_inequalities.py
~~~

Expected exit code: 0; final line: ALL_PAPER_INEQUALITIES_CHECKED.

Both entry points reject optimization; these controls must exit nonzero:

~~~sh
nice -n 19 python3 -O -B verify_r14.py
nice -n 19 python3 -O -B check_inequalities.py
~~~

It reports: optimized Python disables proof assertions.
A changed embedded nonlinear certificate was also rejected after the outer
payload digest was recomputed, because the inner file manifest did not match.

The executable estimates are described in Sections 6 and 8 of the paper.
The embedded explanatory document specifies the coefficient equations and spaces.
The complete analytic proof is in the paper. No separate file is needed to
run the self-contained verifier. Transparent files support inspection and
the additional numerical comparisons.

Recorded packaged single-file wall time: 34.589 seconds; exit code 0 and final line PROVED.
