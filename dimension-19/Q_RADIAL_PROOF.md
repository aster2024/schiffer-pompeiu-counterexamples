# The radial bound for the reciprocal of Q

This note proves the bound used by `verify_axisym_exact_v4.py` in the gate "Gegenbauer Q reciprocal" when the
configuration of a dimension selects it (dimension 19). Notation, spaces and the remaining proof are those of the
part "The exact inverse in general dimension" of the paper (arXiv:2609.35419v2); nothing else in that proof changes.

## Setting

Let n ≥ 5 be an integer, λ = (n − 2)/2 and R > 1. Put

    G_ℓ(t) = C_ℓ^λ(t) / C_ℓ^λ(1),    G_ℓ(cos θ) = Σ_m D_{mℓ} cos(mθ),    h_ℓ(R) = Σ_m D_{mℓ} R^m.

For continuous radial coefficients f_ℓ(r), with f_ℓ(0) = 0 for ℓ > 0, the space Y has the norm

    ‖f‖_Y = Σ_{ℓ≥0} h_ℓ(R) · sup_{0≤r≤1} |f_ℓ(r)|,        f(r, t) = Σ_ℓ f_ℓ(r) G_ℓ(t).

Let b > 0 and let ψ(w) = b w + Σ_{k=1}^{J} c_k w^(2k+1) be the odd polynomial with real coefficients that defines the
shape; all decimal inputs are exact rationals. Put Q = Im ψ / y. Let U_j be the Chebyshev polynomials of the second
kind and write U_{2k}(t) = Σ_{ℓ=0}^{2k} A_{kℓ} G_ℓ(t). With K = 40, define for 0 ≤ s ≤ 1

    P_ℓ(s) = Σ_{k=1}^{40} c_k A_{kℓ} s^k        (0 ≤ ℓ ≤ 80),

    T_40(R) = κ(R) Σ_{k=41}^{J} |c_k| (1 + 2 Σ_{j=1}^{k} R^(2j)),

where κ(R) is the bound of the paper for the passage from Chebyshev to Gegenbauer coefficients, and

    q_radial = Σ_{ℓ=0}^{80} h_ℓ(R) sup_{0≤s≤1} |P_ℓ(s)| + T_40(R).

**Lemma.** Let M_ℓ ≥ sup_{[0,1]} |P_ℓ| for 0 ≤ ℓ ≤ 80 and q̂ = Σ_ℓ h_ℓ(R) M_ℓ + T_40(R). If q̂ < b, then Q is invertible
in Y and

    ‖Q^(−1)‖_Y ≤ 1/(b − q̂),
    ‖b Q^(e−1)‖_Y ≤ b^e (1 − q̂/b)^(−|e−1|),    ‖Q^(−e)‖_Y ≤ b^(−e) (1 − q̂/b)^(−e),    e = (n − 2)/2.

## Proof

*Step 1: the head.* For w = r e^(iθ) and y = r sin θ,

    Im w^(2k+1) / y = r^(2k) sin((2k+1)θ)/sin θ = r^(2k) U_{2k}(cos θ).

Both sides are polynomials in (x, y), so the identity holds on the axis and at the origin. Hence

    Q − b = Σ_{k=1}^{J} c_k r^(2k) U_{2k}(t),

and the part with k ≤ 40 equals Σ_{ℓ=0}^{80} P_ℓ(r²) G_ℓ(t). Since r ↦ r² maps [0, 1] onto [0, 1], the Y-norm of
this part is exactly Σ_{ℓ≤80} h_ℓ(R) sup_{[0,1]} |P_ℓ|. The powers s^k are kept inside P_ℓ: the supremum is taken
after summing over k, with signs.

*Step 2: the tail.* For k > 40 use U_{2k} = 1 + 2 Σ_{j=1}^{k} T_{2j} and sup_r r^(2k) = 1; the estimate of the paper
for Chebyshev polynomials in Y gives ‖c_k r^(2k) U_{2k}‖_Y ≤ |c_k| κ(R) (1 + 2 Σ_{j≤k} R^(2j)). Summing over k > 40
gives ‖Q − b‖_Y ≤ q_radial ≤ q̂. No cancellation between different k is used in the tail.

*Step 3: the coefficients A_{kℓ}.* The verifier builds G_ℓ and U_j in ℚ[t] from the recurrences

    (ℓ + 2λ) G_{ℓ+1} = 2(ℓ + λ) t G_ℓ − ℓ G_{ℓ−1},        U_{j+1} = 2t U_j − U_{j−1},

obtains A_{kℓ} by elimination from the top degree (the leading coefficient of G_ℓ is non-zero because λ > 0), and
checks the identity Σ_ℓ A_{kℓ} G_ℓ = U_{2k} in exact rational arithmetic for every k ≤ 40.

*Step 4: the suprema.* Write P(s) = Σ_{k=0}^{d} a_k s^k in the Bernstein basis B_{i,d}(s) = C(d,i) s^i (1 − s)^(d−i):

    P = Σ_i β_i B_{i,d},        β_i = Σ_{k≤i} a_k C(i,k)/C(d,k),

which follows from Σ_{i≥k} (C(i,k)/C(d,k)) B_{i,d}(s) = s^k. The B_{i,d} are non-negative on [0, 1] and sum to 1, so
|P(s)| ≤ max_i |β_i| there. De Casteljau subdivision at the midpoint gives the Bernstein coefficients of P on [0, 1/2]
and on [1/2, 1] exactly (for the left half, Σ_i B_{i,d}(u) 2^(−i) Σ_{j≤i} C(i,j) β_j = P(u/2); the right half follows by
applying this to P(1 − s)). After D subdivisions the 2^D closed intervals cover [0, 1], so the maximum of the absolute
values of all Bernstein coefficients of all leaves bounds sup_{[0,1]} |P|. The verifier uses D = 6 and evaluates the
coefficients in Arb ball arithmetic from the exact rational a_k, taking upper endpoints; this gives M_ℓ. No sampled
value of P enters.

*Step 5: the reciprocal.* By the non-negative linearization G_i G_j = Σ_ℓ L_{ij}^ℓ G_ℓ with L_{ij}^ℓ ≥ 0 and
Σ_ℓ L_{ij}^ℓ h_ℓ ≤ h_i h_j (Section 1 of that part of the paper), together with sup|fg| ≤ sup|f| sup|g|, the space Y
is a Banach algebra with unit of norm 1. Put H = (Q − b)/b, so ‖H‖_Y ≤ q̂/b < 1. The series b^(−1) Σ_{j≥0} (−H)^j
converges absolutely in Y to Q^(−1), and ‖Q^(−1)‖_Y ≤ b^(−1) Σ_j (q̂/b)^j = 1/(b − q̂). The two bounds for powers
follow from the binomial series of (1 + H)^(e−1) and (1 + H)^(−e), whose coefficients are dominated in absolute value
by those of (1 − x)^(−|e−1|) and (1 − x)^(−e). ∎

## Where the bound is used

The bound q̂ replaces the previous bound for ‖Q − b‖_Y only where the proof needs the norm of Q − b in Y: in the
reciprocal 1/(b − q) and in the binomial factors of the conjugation multipliers. Every estimate of the proof that
needs a majorant of the coefficients of Q, or its degree moments (the bounds for the derivatives of Q, the
degree-weighted norms, the quartic expansion with its constants C₂, C₃, C₄, the univalence and curvature
numerators), keeps the previous coefficientwise bound q_old and its moments, exactly as in the paper. For the eleven
dimensions of the paper the configuration does not select the new bound, and every previously recorded mathematical
field of their receipts is reproduced by `verify_axisym_exact_v4.py`.

In dimension 19 the certified values are b = 0.99002…, q_old = 1.77065… (which exceeds b) and q̂ = 0.46646… < b.
