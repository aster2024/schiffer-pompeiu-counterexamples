# Dimension-four verification

This directory verifies the dimension-four theorem in the accompanying paper.
The self-contained verifier embeds all required programs, coefficient data,
inverse matrix and column certificates. The same files are supplied
transparently here.

Requirements: Python 3, python-flint **0.9.0**. The recorded replay used
Python 3.10.19 on Linux x86-64. No NumPy computation is required by this proof.

Run from this directory, with only one verifier running at a time:

~~~sh
sha256sum -c SHA256SUMS
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  nice -n 19 python3 -B verify_berenstein_r4.py
~~~

Expected exit code: 0. Expected final line:

~~~text
PROVED
~~~

For complete serial regeneration and byte comparison of every bounded column:

~~~sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  nice -n 19 python3 -B verify_berenstein_r4.py --recompute
~~~

Expected exit code: 0; final line: PROVED.
Recorded wall times: 9.514 seconds for ordinary replay and 374.957 seconds
for complete regeneration.

Check the sharper numerical inequalities stated in the paper:

~~~sh
nice -n 19 python3 -B check_inequalities.py
~~~

Expected exit code: 0; final line: ALL_PAPER_INEQUALITIES_CHECKED.
The single-file run enforces the finite inverse defect below 21e-35.
The additional comparison also checks the sharper constants of the paper.

Both entry points reject optimization; these controls must exit nonzero:

~~~sh
nice -n 19 python3 -O -B verify_berenstein_r4.py
nice -n 19 python3 -O -B check_inequalities.py
~~~

It reports that Python optimization disables certificate assertions.
A changed embedded source-tail certificate was also rejected after the
outer payload digest was recomputed, because the inner manifest did not match.

The analytic reduction, support estimates, injectivity, signed boundary
condition, strict convexity and reconstruction are proved in the paper.
The programs enclose their numerical expressions. The single file can be
copied into an otherwise empty directory and run with the same dependencies.
