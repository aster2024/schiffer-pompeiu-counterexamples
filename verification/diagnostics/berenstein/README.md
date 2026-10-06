# Numerical diagnostics for the planar Berenstein domain

These floating-point scripts describe the frozen finite centre and the
comparisons in “Relation with earlier work”. They are separate from the
interval certificate. They use NumPy 1.26.4, SciPy 1.15.2 and mpmath 1.3.0;
`csb_nonconvex.py` also uses python-flint 0.9.0. Run them from this
directory. `geom.py` loads the frozen centre from
`../../berenstein/centre_frozen.json`.

`python geom.py` evaluates curvature and radial ranges; `python index_count.py`
counts disc eigenvalues; `python defect.py` and `python supportfn.py` evaluate
support-function data. `python csb_convex.py` evaluates the appendix centre
of Colbrook–Sadeghi–Stepaniants. `python csb_nonconvex.py` encloses a
finite triangle witness in Arb and transfers it using their pointwise map
error; its output is `csb_nonconvex.json`. `python table.py` reproduces the stationary-phase
comparison from the included sampled records. `remainder.py` and
`remainder_scan.py` generate new sampled records and may take substantially
longer than the table calculation. None of these sampled maxima is an interval
bound for the exact domain.
