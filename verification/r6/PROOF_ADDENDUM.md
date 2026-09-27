# Proof addendum to the R⁶ certificate (release of 2026-09-28)

This addendum accompanies `FINAL.md` (the proof notes, in Chinese) and the
verifier `verify_r6.py`. It records the changes requested by the audit of the
R⁶ package and gives complete proofs of three statements that the audit asked
to have written out:

1. the fail-closed maximum in the verifier (§1);
2. the principal shape-column identity (5.1) for **every** `j ≡ 1 (mod 4)`,
   `j ≥ 5` (§2), and the factorisation `Kg = (1−r²)²Q` used by the far-shape
   bound (§3), both also checked by the shipped script `check_identities_r6.py`;
3. convexity of the meridian domain `D` for every shape in the certified ball,
   with the certified bound `Re(1 + wψ''/ψ') > 0.71569200` on the closed disc,
   and the implication "`D` convex ⇒ `Ω ⊂ R⁶` convex" (§4), checked by
   `verify_convex_r6.py`.

Notation follows `FINAL.md`: `w = re^{iθ}`, the disc polynomials of the sine
sector are

```text
S_{n,s}(r,θ) = rⁿ P_s^{(0,n)}(2r²−1) sin(nθ),   n ≡ 2 (mod 4), s ≥ 0,
```

with `P_s^{(α,β)}` the Jacobi polynomial in the standard normalisation
(`P_s^{(0,n)}(1) = 1`); `q = r²`, `u = q − 1`; the unknowns are
`x = (g, c)`, `ψ_c(w) = Σ_{j≡1 (4)} c_j w^j` with real `c_j`, and

```text
H(g,c) = g + |ψ_c'|² (Kg + h_c),   h_c = Im(ψ_c²)/2,
‖c‖_C = Σ_{j≡1 (4)} ω_j |c_j|,   ω_j = α (j+2) ρ^{j+1},   α = b³/2,
b = c°₁ = 28.864750643300045…,   ρ = 11/10,   r = 1/5000.
```

`verify_r6.py --stage final` proves (status `PROVED`, `certificate_r6.json`)
that `H` has a unique zero `x* = (g*, c*)` with `‖x* − x°‖_X ≤ r`; in
particular `Σ_j ω_j |c*_j − c°_j| ≤ r`.

---

## 1. Fail-closed maximum (audit item 1)

The released `max_upper` of `verify_r6.py` was `max(v.upper() for v in seq)`.
Python's `max` compares with `>`, and every comparison with a NaN is false, so a
NaN that is not the first element was silently dropped:
`max_upper([1, nan])` returned `1`. (It never happened in the certified run:
the audit re-ran every batch with a strict maximum and all values were finite.)

The release version is identical in spirit to `upper_max` of the R⁴ verifier
`verify_r4.py`: every element must satisfy `is_finite()` and have a finite upper
endpoint, otherwise `need(...)` raises and the verifier stops before printing
`PROVED`. In addition, every Arb value read back from a JSON receipt passes
through `_interval`, which rejects NaN, `inf`, and balls of infinite radius.
All other comparisons in the final gate are strict Arb comparisons (`<`, `>`),
which are false for NaN, so they also fail closed. `test_fail_closed_r6.py`
checks these negative cases (`nan` after a finite value, `nan` first,
`inf`, infinite radius; receipt strings `nan`, `inf`, `[1 +/- inf]`).

The interval receipts do not depend on the change (they are reproduced byte for
byte by `reproduce_r6.sh`); the certificate changes only because it now also
lists the hashes of the five stage scripts and of `verify_r6.py` itself.

---

## 2. The principal shape column (5.1) for general j

**Proposition 2.1.** Let `x_b = (0, ψ_b)`, `ψ_b(w) = b w`, be the disc
(`g = 0`). For every `j ≡ 1 (mod 4)`, `j ≥ 5`,

```text
T_j := DH(x_b) e_j
     = b³ { (j+2)/2 · S_{j+1,0} − (j−2)/2 · S_{j−3,0}
            − j/(j+1) · S_{j−3,1} − 1/(j+1) · S_{j−3,2} }.          (5.1)
```

**Lemma 2.2 (radial expansion).** For every integer `m ≥ 0`,

```text
r^{m+4} sin(mθ) = (m+1)/(m+3) S_{m,0} + 2/(m+4) S_{m,1}
                  + 2/((m+3)(m+4)) S_{m,2}.
```

*Proof.* Divide by `r^m sin(mθ)`; the claim is the polynomial identity
`q² = (m+1)/(m+3) P_0 + 2/(m+4) P_1 + 2/((m+3)(m+4)) P_2` with
`P_s = P_s^{(0,m)}(2q−1)`. By the explicit formula (DLMF 18.5.7; equivalently
the formula for `a^{(n)}_{s,k}` in the R⁴ paper) in the variable `u = q − 1`,

```text
P_0 = 1,   P_1 = 1 + (m+2) u,   P_2 = 1 + 2(m+3) u + (m+3)(m+4) u²/2 .
```

The right-hand side is therefore `c₀ + c₁u + c₂u²` with

```text
c₂ = 2/((m+3)(m+4)) · (m+3)(m+4)/2 = 1,
c₁ = 2(m+2)/(m+4) + 4(m+3)/((m+3)(m+4)) = (2m+8)/(m+4) = 2,
c₀ = [(m+1)(m+4) + 2(m+3) + 2] / ((m+3)(m+4)) = (m²+7m+12)/((m+3)(m+4)) = 1,
```

which is `(1+u)² = q²`. ∎

*Proof of Proposition 2.1.* Differentiating `H` in the direction `e_j`
(`ψ ↦ ψ + ε w^j`) gives, for any `x = (g, c)`,

```text
DH(x) e_j = 2 Re( j w^{j−1} \overline{ψ_c'} ) (Kg + h_c) + |ψ_c'|² Im(ψ_c w^j),
```

because `∂_ε |ψ' + ε j w^{j−1}|² = 2 Re(j w^{j−1} \overline{ψ'})` and
`∂_ε Im((ψ + ε w^j)²)/2 = Im(ψ w^j)`. At `x_b` we have `g = 0`, `ψ_b' = b`,
`h_b = Im(b²w²)/2 = b² r² sin(2θ)/2`, hence

```text
T_j = b³ [ j Re(w^{j−1}) Im(w²) + Im(w^{j+1}) ]
    = b³ [ j r^{j+1} cos((j−1)θ) sin(2θ) + r^{j+1} sin((j+1)θ) ].
```

With `2 cos((j−1)θ) sin(2θ) = sin((j+1)θ) − sin((j−3)θ)`,

```text
T_j = b³ [ (j+2)/2 · r^{j+1} sin((j+1)θ) − j/2 · r^{j+1} sin((j−3)θ) ].
```

The first term is `(j+2)/2 · S_{j+1,0}` (a harmonic polynomial, `P_0 = 1`). For
the second apply Lemma 2.2 with `m = j − 3 ≥ 2`:

```text
−j/2 · r^{j+1} sin((j−3)θ)
  = −j/2 · [ (j−2)/j · S_{j−3,0} + 2/(j+1) · S_{j−3,1} + 2/(j(j+1)) · S_{j−3,2} ]
  = −(j−2)/2 · S_{j−3,0} − j/(j+1) · S_{j−3,1} − 1/(j+1) · S_{j−3,2}. ∎
```

Both frequencies `j+1` and `j−3` are `≡ 2 (mod 4)`, so `T_j` lies in the sine
sector. **Remark (consistency with (5.2)).** The harmonic (`s = 0`) part of
`T_j` puts `α(j+2) = α(n+1)` on row `n = j+1` and `−α(j−2) = −α(n+1)` on row
`n = j−3`. Hence a shape vector `η` produces on the harmonic tail row `n` the
value `α(n+1)(η_{n−1} − η_{n+3})`, which is (5.2); the non-harmonic
contributions `−b³ j/(j+1)`, `−b³/(j+1)` on rows `(j−3, 1)` and `(j−3, 2)` are,
with `n = j−3` and after moving them to the right-hand side (hence the sign),
exactly the corrections `+b³(n+3)/(n+4)·η` and `+b³/(n+4)·η` added in
`Inverse.split` of `verify_r6.py`.

**Machine checks (`check_identities_r6.py`, receipt `identities_r6.json`).**

| check | range | result |
|---|---|---|
| II.1 exact (`fmpq_poly`) radial expansion of Lemma 2.2 | `m = j−3`, all `j ≡ 1 (4)`, `5 ≤ j ≤ 2001` | 500 exact identities |
| II.2 the verifier's own column algebra (`zw`/`hol`/`sine`/`unsine`, as in `shape_tail_interval_r6.py`) at the ball with `b = 1`, run in exact rationals, equals `principal(j, 1)` | `5 ≤ j ≤ 2001` | 500 exact equalities (also equal supports) |
| II.3 same with Arb, `b = c°₁`, 128 bits: every coefficient ball of the difference contains 0 | `5 ≤ j ≤ 1001` | 250 values, max |diff| `≤ 2.21·10⁻²⁸` |
| II.4 Arb pointwise: expansion (5.1) vs the direct formula `b³[j Re(w^{j−1}) Im(w²) + Im(w^{j+1})]` | `j ∈ {5,9,61,101,397,1001}`, 4 points in the closed disc | max |diff| `≤ 8.2·10⁻²⁹` |

(`T_j` is homogeneous of degree 3 in `b`, so II.2 with `b = 1` is exact for
every `b`.)

---

## 3. The factorisation `Kg = (1 − r²)² Q`

**Lemma 3.1.** For every integer `n ≥ 0` and `s ≥ 1`, with `d = n + 2s`, the
operator `K` of (3.1) satisfies

```text
(a) Δ K S_{n,s} = S_{n,s};
(b) K S_{n,s} = 0 and ∇ K S_{n,s} = 0 on ∂𝔻;
(c) K S_{n,s} = (1−r²)² rⁿ P_{s−1}^{(2,n)}(2r²−1) sin(nθ) / (4s(s+1));
(d) the absolute column sum of (3.1) is 1/(d(d+2)); in the R⁶ sector d ≥ 4,
    hence ‖K‖ ≤ 1/24.
```

*Proof.* This is Lemma "K" of the R⁴ paper (`paper/main.tex`, §4, stated there
for odd `n`); its proof of (a)–(c) uses only the explicit coefficients
`a^{(n)}_{s,k}` of `P_s^{(0,n)}(2t−1)` and orthogonality in `L²(tⁿ dt)`, never
the parity of `n`, and therefore applies verbatim to `n ≡ 2 (mod 4)`. For
completeness, (c): by (b) the radial factor `R(q)` of `K S_{n,s} = rⁿ R(r²) sin nθ`
is a polynomial of degree `s+1` with `R(1) = 0` and `nR(1) + 2R'(1) = 0`, so
`R = (1−q)² T` with `deg T = s − 1`. Being a combination of
`P_{s−1}, P_s, P_{s+1}` of the family `(0, n)`, `R` is orthogonal in
`L²([0,1], qⁿ dq)` to all polynomials of degree `< s − 1`; hence `T` is
orthogonal to them in `L²([0,1], qⁿ(1−q)² dq)`, i.e. `T` is a multiple of
`P_{s−1}^{(2,n)}(2q−1)`. The leading coefficient of `P_s^{(α,β)}(2q−1)` in `q`
is `binom(2s+α+β, s)`, so comparing leading coefficients,
`binom(n+2s+2, s+1)/(4(d+1)(d+2)) = binom(n+2s, s−1)/(4s(s+1))`, the factor is
`1/(4s(s+1))`. (d): `1/(4d(d+1)) + 1/(2d(d+2)) + 1/(4(d+1)(d+2)) = 1/(d(d+2))`. ∎

**Connection coefficients.** By DLMF 18.9.5,
`(2k+α+β+1) P_k^{(α,β)} = (k+α+β+1) P_k^{(α+1,β)} − (k+β) P_{k−1}^{(α+1,β)}`.
With `β = n` and `α = 0, 1`:

```text
P_k^{(1,n)} = [ (2k+n+1) P_k^{(0,n)} + (k+n) P_{k−1}^{(1,n)} ] / (k+n+1),
P_k^{(2,n)} = [ (2k+n+2) P_k^{(1,n)} + (k+n) P_{k−1}^{(2,n)} ] / (k+n+2).
```

These are exactly the two recurrences in `far_shape_r6.q_of_g` (`one` and
`two`); by induction all connection coefficients
`P_{s−1}^{(2,n)} = Σ_{t ≤ s−1} C^{(n)}_{s−1,t} P_t^{(0,n)}` are **positive**.
Therefore, for the finite centre `g° = Σ g°_{n,s} S_{n,s}` (`n ≤ 58`, `s ≤ 30`),

```text
K g° = (1−r²)² Q,   Q = Σ_{n,s} g°_{n,s}/(4s(s+1)) Σ_t C^{(n)}_{s−1,t} S_{n,t}.
```

**Rigour of the Arb enclosure of `Q`.** `far_shape_r6.py` evaluates exactly this
rational recurrence in Arb ball arithmetic; Arb's containment guarantee then
gives balls containing the exact rational coefficients of `Q`, which is all the
far-shape bound needs (it uses only upper bounds of `|Q_{m,t}|`, together with
the double-zero lemma `‖(1−r²)² w^N Φ_{m,t}‖ ≤ 4(2t+1)(2t+3)ρ^{m+N}/(m+N+1)²`
of the R⁴ paper, whose proof is independent of `n`).

**Machine checks (`check_identities_r6.py`).**

| check | range | result |
|---|---|---|
| Jacobi formulas: DLMF 18.5.7 (in `u`) equals the R⁴ paper's `a^{(n)}_{s,k}` formula; normalisation; orthogonality of `P^{(2,n)}` in `L²(qⁿ(1−q)²dq)` | `n ≤ 202`, `s ≤ 81` (orthogonality `n ∈ {2,6,58}`, `s ≤ 11`) | 4182 exact pairs |
| I.1 exact: radial part of (3.1) equals (c); `4qR'' + 4(n+1)R' = P_s^{(0,n)}` (i.e. (a)); both traces vanish (b); column sum (d) | all `n ≡ 2 (4)`, `n ≤ 202`, `1 ≤ s ≤ 80` | 4080 exact pairs |
| I.2 exact: the byte code of `far_shape_r6.q_of_g`, run with rationals instead of Arb balls, reproduces `P_{s−1}^{(2,n)}/(4s(s+1))` in the `P^{(0,n)}` basis, with positive coefficients | all centre indices `n ≤ 58`, `s ≤ 30` | 450 exact pairs |
| I.3 Arb, 128 bits: every coefficient ball of `K g° − (1−|w|²)² Q` (computed with the verifier's `zw` algebra) contains 0 | 480 coefficients | weighted norm `≤ 2.02·10⁻³¹` (`‖Kg°‖_ρ ≈ 8352.74`) |

---

## 4. Convexity

Let `B_C = { c : Σ_{j≡1 (4)} ω_j |c_j − c°_j| ≤ r }` (real sequences on the
index set `j ≡ 1 (mod 4)`, infinite tail included). By §1–§6 of `FINAL.md`,
`c* ∈ B_C`.

**Proposition 4.1.** For every `c ∈ B_C` and every `|w| ≤ 1`,

```text
|ψ_c'(w)| ≥ L := b − A₁ − E₁ > 27.67533225,
|ψ_c''(w)| ≤ U := A₂ + E₂ < 7.86831783,
Re(1 + w ψ_c''(w)/ψ_c'(w)) ≥ 1 − U/L > 0.71569200,
```

where `A₁ = Σ_{j≥5} j|c°_j| < 1.18941839`, `A₂ = Σ_{j≥5} j(j−1)|c°_j| < 7.86831777`,

```text
E₁ = r · max( 1/ω₁ , 1/(αρ⁶) ) < 9.39·10⁻⁹,     E₂ = r / (αρ e log ρ) < 5.84·10⁻⁸.
```

*Proof.* Put `δ = c − c°`, so `Σ_j ω_j|δ_j| ≤ r`. For `|w| ≤ 1`,
`|ψ_c'(w)| ≥ c₁ − Σ_{j≥5} j|c_j| ≥ b − A₁ − (|δ₁| + Σ_{j≥5} j|δ_j|)` and
`|ψ_c''(w)| ≤ Σ_{j≥5} j(j−1)|c_j| ≤ A₂ + Σ_{j≥5} j(j−1)|δ_j|`.

*First derivative.* `|δ₁| + Σ_{j≥5} j|δ_j| = Σ_j κ_j ω_j|δ_j| ≤ r sup_j κ_j` with
`κ₁ = 1/ω₁ = 1/(3αρ²)` and, for `j ≥ 5`,
`κ_j = j/ω_j = j/((j+2) α ρ^{j+1}) < 1/(αρ^{j+1}) ≤ 1/(αρ⁶)`.

*Second derivative.* For `j ≥ 5`, `j(j−1)/ω_j = j(j−1)/((j+2)αρ^{j+1}) <
j ρ^{−j}/(αρ)`, and the function `x ↦ xρ^{−x}` on `x > 0` has derivative
`ρ^{−x}(1 − x log ρ)`, hence its maximum `1/(e log ρ)` at `x = 1/log ρ`.
Therefore `Σ_{j≥5} j(j−1)|δ_j| ≤ r/(αρ e log ρ) = E₂`. Both suprema are taken
over **all** `j`, so the infinite tail of `c` is covered without truncation.

Finally `Re(1 + wψ''/ψ') ≥ 1 − |w||ψ''|/|ψ'| ≥ 1 − U/L`. The numerical values
are the Arb enclosures in `convex_r6.json` (128 bits):
`L ∈ [27.675332251016402… ± 2·10⁻³⁷]`, `U ∈ [7.868317826804946… ± 5·10⁻³⁸]`,
`1 − U/L ∈ [0.715692019324683… ± 4·10⁻³⁹]`; the script also certifies the
rational bound `1 − U/L > 71569200/10⁸`. ∎

For comparison, the audit's estimate (sup over `j ≤ 2001` only, hence not a
proof for the tail) was `≥ 0.715692019768`; the rigorous bound differs only in
the eighth significant digit because the tail perturbation is `O(10⁻⁸)`.

**Corollary 4.2 (the meridian domain is strictly convex).** Let `ψ = ψ_{c*}`,
`D = ψ(𝔻)`. Then `ψ` is univalent on `𝔻̄` (indeed on `|w| < 21/20` by the
extension-margin gate of `verify_r6.py`), `ψ' ≠ 0` on `𝔻̄` (`L > 0`), and for
the analytic Jordan curve `z(θ) = ψ(e^{iθ})` one has
`z'(θ) = i e^{iθ} ψ'(e^{iθ}) ≠ 0` and

```text
d/dθ arg z'(θ) = Re(1 + wψ''(w)/ψ'(w)) |_{w=e^{iθ}} > 0.71569200 .
```

The winding number of `ψ'(e^{iθ})` about 0 is zero (no zeros of `ψ'` in `𝔻̄`),
so the tangent direction turns monotonically through exactly `2π`; the
curvature of `∂D` is `Re(1+wψ''/ψ')/|ψ'| > 0.71569200/(b + A₁ + E₁) > 0.0238`.
A closed analytic Jordan curve of everywhere positive curvature whose tangent
turns once bounds a strictly convex domain (classical; equivalently Study's
theorem: `Re(1 + wψ''/ψ') > 0` in `𝔻` iff `ψ` maps `𝔻` univalently onto a
convex domain, Duren, *Univalent Functions*, Ch. 2). Hence `D` is strictly
convex. ∎

**Lemma 4.3 (symmetry of D).** For every real `c` supported on `j ≡ 1 (mod 4)`,
`ψ_c(\bar w) = \overline{ψ_c(w)}` and `ψ_c(iw) = i ψ_c(w)` (as `i^j = i`). Hence
`D = ψ_c(𝔻)` is invariant under `z ↦ \bar z` and `z ↦ iz`, i.e. under the
dihedral group of the square. In the coordinates `z = y₁ + i y₂` this contains
the two reflections `σ₁(y₁, y₂) = (−y₁, y₂)` (`z ↦ −\bar z`) and
`σ₂(y₁, y₂) = (y₁, −y₂)` (`z ↦ \bar z`), and the exchange
`(y₁, y₂) ↦ (y₂, y₁)` (`z ↦ i\bar z`). ∎

**Lemma 4.4 (rectangles).** If `D` is convex and invariant under `σ₁, σ₂`, and
`(m₁, m₂) ∈ D` with `m₁, m₂ ≥ 0`, then the closed rectangle
`[−m₁, m₁] × [−m₂, m₂] ⊂ D`.

*Proof.* The four points `(±m₁, ±m₂)` lie in `D` (images of `(m₁, m₂)` under
`id, σ₁, σ₂, σ₁σ₂`). The rectangle is the convex hull of these four points
(each `(t₁, t₂)` with `|t_i| ≤ m_i` is
`Σ_{ε₁,ε₂=±1} λ_{ε₁}μ_{ε₂}(ε₁m₁, ε₂m₂)` with `λ_± = (1 ± t₁/m₁)/2`,
`μ_± = (1 ± t₂/m₂)/2`, reading `λ_± = 1/2` if `m₁ = 0` and likewise for `μ`),
and a convex set contains the convex hull of any of its points. ∎

**Proposition 4.5 (Ω is convex).** Let
`Ω = { x = (x', x'') ∈ R³ × R³ : (|x'|, |x''|) ∈ D }`. Then `Ω` is a bounded open
convex subset of `R⁶`.

*Proof.* `Ω` is the preimage of the open bounded set `D` under the continuous
map `x ↦ (|x'|, |x''|)`, hence open, and bounded because `|x|² = |x'|² + |x''|²`
is bounded on it. Let `x = (x', x'')`, `z = (z', z'') ∈ Ω` and `t ∈ [0, 1]`.
Put `p = (|x'|, |x''|) ∈ D`, `p̃ = (|z'|, |z''|) ∈ D` and

```text
m = (m₁, m₂) = t p + (1−t) p̃ = ( t|x'| + (1−t)|z'| , t|x''| + (1−t)|z''| ).
```

By convexity of `D` (Corollary 4.2), `m ∈ D`, and `m₁, m₂ ≥ 0`. For
`y = t x + (1−t) z = (t x' + (1−t) z', t x'' + (1−t) z'')` the triangle
inequality in each factor `R³` gives

```text
0 ≤ |t x' + (1−t) z'| ≤ t|x'| + (1−t)|z'| = m₁,
0 ≤ |t x'' + (1−t) z''| ≤ t|x''| + (1−t)|z''| = m₂,
```

so `(|y'|, |y''|) ∈ [0, m₁] × [0, m₂] ⊂ D` by Lemma 4.4 (applicable by Lemma 4.3).
Hence `y ∈ Ω`. ∎

Only the two reflections are used; the exchange symmetry is not needed for
convexity. With Proposition 4.5, the R⁶ domain of `FINAL.md` is a bounded
**convex** domain with real-analytic boundary diffeomorphic to `S⁵`, not a
ball, carrying a nonconstant solution of `Δu + u = 0`, `u = 1`, `∇u = 0` on
`∂Ω`; it fails the Pompeiu property. The convexity statement rests on
`verify_r6.py` (`PROVED`) and `verify_convex_r6.py` (`PROVED_CONVEX`).

---

## 5. Audit items and where they are addressed

| audit item | resolution | files |
|---|---|---|
| `max_upper` could drop NaN | fail-closed as in R⁴; receipts rejected if non-finite; negative tests | `verify_r6.py`, `test_fail_closed_r6.py` |
| `reproduce_r6.sh` called bare `python`, 25 s cap | explicit interpreter (`$PYTHON`, else pinned micromamba `python3`, else `python3`), python-flint 0.9.0 check, 600 s cap; also reproduces certificate, identities and convexity receipts byte for byte | `reproduce_r6.sh` |
| re-freeze hashes and certificate | stage-script hashes added to `EXPECTED_SHA256`; certificate records the verifier's own hash; `PROVED` re-confirmed | `verify_r6.py`, `certificate_r6.json` |
| `Kg = (1−r²)²Q` and (5.1) checked only for a few cases | exact + Arb checks for many `n, s, j`; general proofs §2–§3 | `check_identities_r6.py`, `identities_r6.json` |
| convexity estimate was not part of a certificate | rigorous Arb certificate for every shape in the ball, tail included; implication `D` convex ⇒ `Ω` convex proved | `verify_convex_r6.py`, `convex_r6.json` |
