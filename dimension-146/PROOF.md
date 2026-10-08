# INF3: analytic proof behind the finite certificate

The theorem below is asserted only after the final phase of `verify_inf3_hardened.py` prints `PROVED_INF3_N146`. This text supplies the analytic lemmas; numerical inequalities and artifact bindings are checked by the verifier. There is no appeal to an asymptotic threshold or to an old existence receipt.

## 1. Statement and normalization

Let p=73, a=54, b=92, A=a/2, B=b/2. The exact rational lists in center.json specify a finite degree-96 boundary h and a finite regular Helmholtz field. R is the unique zero of J_73 in the rational bracket checked by Arb: R=j_(73,17)=149.0253297755978772482... . The centre in X is the exactly clamped pair (g0,f0), where f0_j=h_j H_j(1)/R and g0=Delta[u_app(R Phi_f0)-ell-1], with ell the exact double-trace lift specified below. Only the input lists h_j,c_j are rational; R, f0 and the lifted g0 need not be rational. The certified distance is ||g-g0||C+||f-f0||E1<=10^-46. The claim is a bounded, real-analytic, strictly convex nonball O(54) x O(92)-invariant domain in R^146 carrying a nonconstant real u with (Delta+1)u=0, u=1 and partial_nu u=0 on its boundary.

Put r=|y|, x=|X|^2/r^2, t=r^2. The quarter-disc operator is partial_ss+(a-1)/s partial_s+partial_tt+(b-1)/t partial_t in the block radii. At either axis it means the regular Euclidean lift. Invariants are polynomials in |X|^2,|Y|^2. Harmonic decomposition gives exactly one invariant harmonic in each even degree 2j and none in odd degree. Let G_j be the Jacobi polynomial for Beta(A,B), normalized by G_j(1)=1, and S_j=r^(2j)G_j(x). The scalar radial angular eigenvalue is 2j(2j+2p-2).

The classical Jacobi positive-linearization theorem is used in its explicit valid parameter range alpha=B-1>=beta=A-1>=0. Its sufficient Gasper margin is positive for p>=12 and a/(2p) in [1/4,3/8], as checked by the polynomial inequality in the verifier. Thus G_i G_j=sum c_(i,j,k)G_k with c>=0, sum c=1, and |i-j|<=k<=i+j. The endpoint bound |G_j|<=1 follows independently from the Sonin function for the Jacobi ODE: its derivative changes sign at most once from negative to positive, so its maximum is at an endpoint; the larger endpoint is x=1.

## 2. Spaces and exact Poisson columns

Fix rho=2 and Omega(D)=(1+D/p^2)^p. Define Psi_(j,s)=S_j P_s^(0,2j+p-1)(2t-1), and

||v||C=sum rho^j Omega(2j+2s)|v_(j,s)|.

For existence and geometry the coefficients are real; complexification is used only for analytic/Cauchy estimates. This weighted l1 completion is a Banach space. It embeds in normalized invariant L2 with norm at most one: orthogonality gives ||Psi_(j,s)||2^2=p/(2j+p+2s)*||G_j||Beta^2<=1. Let C^[k] add the weight (1+2j+2s)^k. The auxiliary second moment used for the inverse adds (1+(2j+2s)/p)^2. Its space equals C^[2], with norms hatC2<=C2<=p^2 hatC2. Consequently a bound Bhat for the inverse in hatC2 gives p^2 Bhat in C2.

The boundary space E1 has norm sum (1+j)rho^j Omega(2j)|f_j|. Positive linearization, the triangle support, Omega(D1+D2)<=Omega(D1)Omega(D2), and 1+i+j<=(1+i)(1+j) prove the boundary algebra estimate. Only solid polynomial multipliers are used on C; C is not assumed to be a field algebra.

Write beta=2j+p-1, d=beta+2s, b0=beta+1, J_s=P_s^(0,beta)(2t-1). The Dirichlet Poisson columns are

K_0=(J_1-J_0)/(4b0(b0+1)),
K_s=J_(s-1)/(4d(d+1))-J_s/(2d(d+2))+J_(s+1)/(4(d+1)(d+2)), s>=1.

For all s>=0 the coefficient of t^m in J_s^beta is a_(s,m)=(-1)^(s-m) binom(s,m) binom(s+beta+m,s). For 0<=m<=s, the ratios a_(s-1,m)/a_(s,m)=-(s-m)/(s+beta+m) and a_(s+1,m)/a_(s,m)=-(s+beta+m+1)/(s+1-m) give the coefficient of K_s as -a_(s,m)/[4(s+1-m)(s+beta+m)]. Indeed the three terms divided by a_(s,m) obey

-(s-m)/[4d(d+1)(s+beta+m)]-1/[2d(d+2)]-(s+beta+m+1)/[4(d+1)(d+2)(s+1-m)] = -1/[4(s+1-m)(s+beta+m)].

The ratio a_(s,m-1)/a_(s,m)=-m(beta+m)/[(s+1-m)(s+beta+m)] therefore proves 4[tK_s''+(beta+1)K_s']=J_s coefficient by coefficient. The top coefficient uses a_(s+1,s+1)/a_(s,s)=(d+1)(d+2)/[(s+1)(s+beta+1)]; the constant of integration is fixed by K_s(1)=0. The separate s=0 formula is (t-1)/[4(beta+1)]. Their boundary trace is zero. The two useful exact identities are

K_s'=(J_s^(beta+1)-J_(s-1)^(beta+1))/(4(d+1)),
beta K_s+tK_s'=(J_(s+1)^(beta-1)-J_s^(beta-1))/(4(d+1))

for s>=1; for s=0 the right sides are 1/(4b0) and J_1^(beta-1)/(4b0). The harmonic column is not doubled. Exact finite controls in the verifier test this essential normalization separately.

Parameter raising and multiplication/lowering are positive mass-one stencils:

J_s^beta=[(s+beta+1)J_s^(beta+1)+s J_(s-1)^(beta+1)]/(2s+beta+1),
t J_s^(beta+1)=[(s+beta+1)J_s^beta+(s+1)J_(s+1)^beta]/(2s+beta+2).

These contiguous relations follow from the same a_(s,m): the two coefficients in the raising formula, divided by a_(s,m), sum to [(s+beta+1)(s+beta+m+1)-s(s-m)]/[(2s+beta+1)(beta+m+1)]=1; the lowering identity reduces to [(s+beta+1)-(s+1)(s+beta+m+1)/(s+1-m)]/(2s+beta+2)=-m/(s+1-m). Leading coefficients handle m=s+1, and the zero lower-index terms handle s=0. Differentiating K_s gives coefficient 1/[4(beta+m+1)] times a_(s,m), exactly the difference in parameter beta+1. Similarly (beta+m)K_(s,m) equals the difference in parameter beta-1 because

[-(beta+m)/(s+1-m)-(beta+m)/(s+beta+m)]/[4(d+1)] = -(beta+m)/[4(s+1-m)(s+beta+m)].

Thus the derivative and saturated identities hold for every index, with their stated s=0 replacements. They justify every solid multiplier conversion, with total output degree bounded by input degree plus multiplier degree.

For completeness the angular recurrence used by every phase has monic coefficients
H_(j,m)=(-1)^(j-m) binom(j,m) (A+m)_(j-m)/(j+p-1+m)_(j-m),
where (z)_n is the rising factorial. The Jacobi hypergeometric polynomial with parameters (B-1,A-1), argument 2x-1, divided by its leading coefficient, has these coefficients. The ratios
H_(j,m-1)/H_(j,m)=-m(A+m-1)/[(j+1-m)(j+p+m-2)],
H_(j-1,m)/H_(j,m)=-(j-m)(2j+p-2)(2j+p-3)/[j(A+j-1)(j+p+m-2)],
H_(j+1,m)/H_(j,m)=-(j+1)(A+j)(j+p+m-1)/[(j+1-m)(2j+p)(2j+p-1)]
prove H_(j+1)=(x-b_j)H_j-a_j H_(j-1) using the literal a_j,b_j of the verifier, for j>=1; its leading coefficient is one and H_0=1,H_1=x-A/p. The endpoint is H_j(1)=(B)_j/(p+j-1)_j, so H_(j+1)(1)/H_j(1)=(B+j)(p+j-1)/[(p+2j-1)(p+2j)]. Dividing the recurrence by endpoints gives exactly the x-multiplication stencils; up+diagonal+down=1 is an exact rational identity, also at j=0. These proofs cover all indices used, rather than extrapolating finite tests. In angular mode j, partial_r K^D v at r=1 is v_(j,0)/(2(2j+p)). Hence ||partial_r K^D||C->E1<=1/2. The closed subspace G of vanishing harmonic source layers v_(j,0)=0 is precisely the compatible clamped source space. Put K=K^D|G and X=G x E1 with the sum norm. All identities extend from polynomials by the C-to-L2 embedding and Poisson uniqueness.

## 3. Radial chart and exact operator formula

For f in E1 let H=P f=sum f_j S_j, q=1+H, d=q+E q, where E=y dot grad. Set Phi_f(y)=q(y)y. On the chart ball, reciprocal multipliers are convergent series because ||q-1||mult<=||f||E1 and ||d-1||mult<=2||f||E1. On real shapes, ||D(Hy)||infinity<=4||f||E1, giving injectivity whenever this is <1.

Let m=grad q dot grad q and k=grad q dot grad E q. The inverse derivative of Phi is q^-1(I-y tensor grad q/d). Applying the Hessian chain rule, including the derivative of this inverse, gives

Delta_Phi u=q^-2 Delta u -2q^-2 d^-1 grad q dot grad E u +m q^-2 d^-2 E(E+1)u +(2m q^-1 d^-3-Cq/d)Eu,
Cq=q^-2[-2(k-m)/d+m E(E-1)q/d^2].

The first-order drift is included. The independently implemented meridian-coordinate check uses the full 2x2 Jacobian and both block drifts. Clearing denominators gives the exact identity

T=q^2 d^4 (Delta_Phi+R^2)K^D
 =d^4 I-2d^3 grad q dot grad E K^D +m d^2 E(E+1)K^D
  +[2m q d+2d^2(k-m)-m d E(E-1)q]E K^D+R^2 q^2 d^4 K^D.

For an angular product output of degree C=2k from degrees A=2i,B0=2j, put h=i+j-k. The exact gradient-Euler column is

A t^h J_s +[4h(h+C+p-1)-2AB0-4A(p-1)]t^h K_s'
 +2B0 h(h+C+p-1)t^(h-1)K_s.

At the saturated lower output C=B0-A it is grouped as A t^A J_s+2AB0 t^(A-1)(beta K_s+tK_s'). This grouping is necessary for legal parameter lowering. Other outputs have enough radial powers for positive conversion. To derive the formula for every index, write w=B0 K_s+2tK_s'. Poisson's equation gives w'=J_s/2-(B0+2p-2)K_s' and w''=J_s'/2-(B0+2p-2)K_s''. In an output harmonic S_k, the identity 2 grad S_i dot grad S_j=Delta(S_i S_j) gives the profile 2h(h+C+p-1)t^(h-1) w+2A t^h w'. Substituting w,w' yields the displayed three terms. When k=j-i and i<=j, h=A and C=B0-A, so the last two terms combine to 2AB0 t^(A-1)(beta K_s+tK_s'). This proves the saturated grouping and the raising count for all indices. The 7A bound is used for the fixed chart with the explicitly gated 2 qmax<=p-1. Exact squared-radius differentiation controls in the identities phase check the implementation separately.

The scalar pull-back is an all-dimensional chain-rule identity. The inverse differential is (DPhi)^-1=q^-1(I-y tensor grad q/d), so physical coordinate differentiation on a pulled-back scalar is D_i=q^-1(partial_i-q_i E/d). Expanding sum_i D_i D_i, using Delta q=0, E d=E q+E^2 q, and sum_i q_i partial_i(E U)=grad q dot grad(E U), gives exactly the displayed second-order and drift terms. This derivation includes the differentiation of q^-1 and d^-1; clearing by q^2 d^4 gives P0 as stated. The identities phase additionally compares both sides with exact inverse-Jacobian jets in squared block coordinates, including the axes and first-order drift.

## 4. Infinite coefficient inverse by finite head and rigorous tails

The core operator T_* is defined by literal dyadic polynomial coefficients and exact Jacobi/Poisson stencils. The component verifier encloses the difference between this model and the actual order-16 chart. The final matrix builder consumes exactly its validated angular and weight proposals. FINITE_ARITHMETIC.md, embedded with this proof, bounds all matrix assembly and product errors; no eigenvalue or LU accuracy is assumed.

Let P select j<=48,s<128, and Q=I-P. The finite candidate C is an arbitrary dyadic matrix, not assumed to be an inverse. Put A0=diag(C,I). The retained polynomial support gives a total-degree band of 84, angular band 42 and radial-index band 84. Thus every Q-input able to reach P is contained in j<=90,s<214. The verifier checks these inequalities from the actual retained supports. Intermediate truncations cannot remove a returning path: raises decrease radial indices before lowers increase them, and the composite support gates cover the K and gradient intermediates.

Here are the uniform column bounds used for all omitted indices. For input degree D=2j+2s, d0=p-1+D, write omega=Omega(2). After adding an auxiliary moment k=0 or 2, a term raising total degree by L costs at most (1+L/(p+D))^k. Before this moment factor,

||K_s||col<=omega/[d0(d0+2)],
||E K_s||col<=2omega/(d0+1),
||E(E+1)K_s||col<=omega[2+(2p-3)/(d0+1)].

These hold for s=0 as well, using its separate column. The gradient-EK absolute column bound is 7A rho^i Omega(A), A=2i: the three nonsaturated terms cost at most A,4A,2A. For s=0 they cost at most A,2A,A. The saturated formula costs at most 2A. Therefore no clamping restriction is needed for this fixed-polynomial inverse calculation.

Applying these explicit bounds term by term to T_* gives a decreasing majorant in the lower input degree. Angular tail inputs have D>=98; radial tail inputs with j<=48 have D>=256. These are bounded separately. For head input leakage, the perturbation T_*-I-R^2 K^D has upward degree band 84. It cannot leave P when D<98-84=14. For all remaining head inputs its full perturbation column is bounded at D=14. The ball K part can leak only from s=127 and its one upward coefficient is enclosed separately. This covers every output row, not just rows retained in a larger matrix.

The finite product verifier encloses C*T_PP-I, C*T_PQ for angular tail inputs, and C*T_PQ for radial tail inputs. Assembly errors are multiplied by ||C|| and added. Let the resulting head residual be e, finite couplings b_a,b_r, full head leakage h, and infinite tail diagonal bounds t_a,t_r. Then

eta=max(e+h,b_a+t_a,b_r+t_r)

bounds ||I-A0 T_*|| on the entire weighted l1 space. If eta<1, A0 T_* is invertible. Since A0 is identity off a finite square block, surjectivity of A0 T_* forces C to be surjective and hence invertible; thus T_* is invertible and ||T_*^-1||<=max(1,||C||)/(1-eta). This is a finite-head/infinite-tail Schur certificate. It is computed in both C and the auxiliary moment norm.

The actual order-96 chart is then included by a second Neumann gate. All shape coefficient tails are enclosed using the same exact R and exact rational center. For either moment norm, let qmax,dmax bound q,d; dq,dd their order-16 differences; F1 and FE1 bound the shape E1 norm and that of its Euler derivative. The solid gradient inequality

||grad H1 dot grad H2||mult<=4(p+1)||f1||E1 ||f2||E1

bounds differences of m and k. Polynomial telescoping then bounds differences of d^4,d^3,m d^2,P0,q^2 d^4, and the gradient term. `verify_inf3_hardened.py inverse` writes out these telescoping bounds. It also adds the independently enclosed core coefficient-export error and the difference in R^2. If B_* delta<1, the full cleared inverse is <=B_*/(1-B_*delta). Multiplying by ||q^2 d^4|| gives the actual J_D=((Delta_Phi+R^2)K^D) inverse. The auxiliary bound multiplied by p^2 is the original C^[2] bound.

No deformed Hilbert inverse is substituted for this coefficient inverse. The full coefficient inverse suffices for the Newton theorem below. The separately enumerated ball spectrum is an independent finite check and is explicitly not a deformed-domain spectral-gap assertion.

## 5. Approximate center, exact clamping, and residual moments

All center decimals are exact rational inputs. The unique R bracket is checked directly with Arb Bessel values. J_(p-1)(R) is nonzero by the recurrence and ODE uniqueness. The finite regular field is a sum of r^(1-p)J_(p-1+2j)(r)H_j(x)/[R^(1-p)J_(p-1)(R)] with the stored amplitudes. It solves the interior Helmholtz equation exactly and has a regular Euclidean lift.

The collar ODE gives Taylor coefficients and a Cauchy remainder on |y|<=eta. The two-component first-order ODE has row norm <=1+(2p-1)/(R_lower-eta)+1+l(l+2p-2)/(R_lower-eta)^2. Its integral equation bounds both f and f' by C=(|f(0)|+|f'(0)|)exp(eta*row). If H bounds the shape multiplier and z=H/eta<1, the tails are C z^L/(1-z) and C z^(L-1)/(1-z). For moment k and initial angular degree j, polynomial degree is at most j+J m. The tail is bounded by C(1+2j+2J L)^k z^L/[1-z exp(k/L)], with L-1 for the derivative. This follows from (1+2j+2J(L+n))^k<=(1+2j+2JL)^k exp(kn/L). Every denominator is checked.

Finite compositions are performed with Arb polynomials, then converted by the positive x-multiplication Jacobi matrix. The remaining monomials are bounded by rho^m Omega(2m)(1+2m)^k. Thus the resulting bounds hold in coefficient norms on the whole boundary, not only at collocation points.

For physical defects D and N, the reference radial defect is Nref=R d_boundary N. Subtract the exact lift

ell=sum_j [D_j S_j(1-j(t-1))+Nref_j S_j(t-1)/2].

Then U0=uapp(R Phi_f)-ell is exactly clamped, g0=Delta(U0-1) belongs to G, and E=F(g0,f)=-B_f ell. The lift Laplacian is harmonic, with coefficients -2j(4j+2p)D_j+(4j+2p)Nref_j. The explicit lift and differentiated bounds are given in CLAMP_AND_MATERIAL_BOUNDS.md, embedded below. They yield E in C, C^[1], C^[2].

Membership does not depend circularly on an inverse. Each of the 97 regular separated fields has an expansion S_j sum_(n>=0) a_(j,n) R^(2j+2n) q^(2j+2n) t^n, with |a_(j,n)| bounded by a fixed j-dependent constant divided by 4^n n! Gamma(n+2j+p). The chart q is a finite solid polynomial. The positive monomial-to-Jacobi conversions and solid multiplier bound give ||S_j t^n||C[k]<=2^j Omega(2j+2n)(1+2j+2n)^k. The multiplier norm of q to a power grows at most geometrically in n at each fixed moment k; its degree grows linearly. Applying Delta or any fixed number of derivatives costs only a polynomial in that degree. The ratio test with the factorial denominator thus gives absolute convergence of the composed field and its derivatives in every fixed C^[k], at this fixed p=73. No large-p geometric-ratio lemma is used. Its boundary Dirichlet and radial-Neumann defects are entire in x; Jacobi orthogonality and analytic polynomial approximation give faster than any fixed geometric decay of their coefficients. Their exact two-layer lift therefore has finite norm at every fixed angular radius and moment. Hence g0=Delta[U0-1] is in every required C^[k] before using g0=J_D^-1(E-R^2) for its numerical norm bound. Exact clamping and the Poisson normal-trace identity put g0 in G.

## 6. Explicit smoothing and nonlinear constants

For clamped inputs s>=1, multiplication by a solid monomial T=S_c t^h of degree D satisfies

||M_T E K||<=200 rho^c Omega(D)/(p+D).

Here is a proof with an explicit constant. In the parameter-raising stencil the down probability is q_s=s/(2s+beta+1). Couple two adjacent indices with one uniform variable. Their separation remains one until they coalesce. With s denoting the upper index of the adjacent pair, at raising step t the meeting hazard is at least 1/(z+t)-2s/(z+t)^2, z=2s+beta+1. For the two differences in K_s the initial upper indices are s and s+1, so z=d+1 and d+3; for K_s' the parameter is beta+1 and z=d+2. Summation gives failure probability <=e^2 z/(z+M), since 2s sum_(t>=0)(z+t)^-2<=2. Hence the l1 distance of the two rows is <=2e^2 z/(z+M).

If D<=8(B0+p), the direct EK column bound and p+D<=9(B0+p)<=9d give constant at most 36. Otherwise angular triangle support gives a legal raising count M>=D/4 before lowering. Decompose K_s into the two adjacent differences in its three-column formula. The coupling estimate gives ||K_s||<=e^2/[d(d+M)], while 2tK_s' costs <=2e^2/(d+M). Therefore EK costs <=3e^2/(d+M)<=12e^2/(p+D). The only extra grade shift is Omega(2)<2, so 24e^2<200 suffices. The s=0 column is excluded exactly because this lemma is used on G.

The gradient identities above give ||grad(P f) dot grad K^D||<=16||f||E1/p (including s=0), and ||grad(P f) dot grad E K||<=14||f||E1. The former actually follows with constant 6; 16 is a safe enlargement. Applying the smoothing lemma to A(A-1)S_i gives ||M_(E(E-1)P f)E K||<=400||f||E1. Applying it to the solid gradient product gives ||M_(grad(P f) dot grad E(P h))E K||<=1600||f||E1||h||E1: the numerator is bounded by 2AB0(p+B0-1), divided by p+A+B0-2, and A,B0>=2 whenever the term is nonzero.

Together with ||EK||<=4/p and ||E(E+1)K||<=9, these bounds control every grouped term of the scalar pullback. On the complex shape ball ||f||E1<=epsilon=1/10, q^-1 and d^-1 have norms <=(1-epsilon)^-1 and (1-2epsilon)^-1, and m has multiplier norm <=4(p+1)epsilon^2. The verifier's explicit formula M bounds T(f)=Delta_Phi K. Cauchy on the remaining shape-ball distance s gives ||DT||<=M/s and ||D2T||<=8M/s^2. Since F(g,f)=(T(f)+R^2 K)g+R^2 is affine in g, L=2||DT||+||D2T||(||g0||+r) is a valid second-derivative bound. No field-field algebra property or unbounded second derivative of the variable shape is used.

## 7. Exact corrected inverse and Newton gate

At the fixed center let Dp=-partial_nunu U0, d=q+E q, fstar=(E U0)/d. Exact clamping gives Dp=(R^2-E_boundary)/(nu^T A nu). The metric numerator, E1 boundary residual and reciprocal are enclosed. Set T_h=(P h)fstar and Q0(delta g,h)=Delta(K delta g-T_h).

For source z put W=K^D z, h=d_boundary partial_nu W/Dp, delta g=z+Delta T_h. The normal trace of T_h is -Dp h/d, so delta g belongs to G and K delta g=W+T_h. This proves Q0(delta g,h)=z. Conversely Poisson uniqueness and the normal trace recover h and delta g, proving the other inverse identity.

Because d is harmonic,

Delta fstar=d^-1(E+2)g0-2d^-2 grad d dot grad E K g0+2d^-3 |grad d|^2 E K g0.

Positive derivative columns give ||E g0||C<=(p+1)||g0||C^[2]. The explicit Bf bound in the verifier follows. Since fstar has zero Dirichlet trace, fstar=K^D Delta fstar, and the gradient-K estimate gives ||Delta[(P h)fstar]||<=(1+32/p)Bf||h||E1. With H=(1/2)||d/Dp||E1, ||Q0^-1||<=1+H[1+(1+32/p)Bf].

The exact material identity is DF(z0)=J_D Q0+E_material, where E_material h=(P h)E E/d. It follows by differentiating the scalar pullback with reference velocity y(P h)/d. Thus A=Q0^-1 J_D^-1 is a two-sided corrected inverse. The verifier bounds its norm alpha, Y=alpha||E||, Z0<=alpha||d^-1||(p+1)||E||C^[2], and L as above.

For the exact radius r=10^-46 the verifier checks

Y+(Z0-1)r+(alpha L/2)r^2<0 and Z0+alpha Lr<1.

The Banach fixed-point theorem for z->z-AF(z) gives a zero in that X-ball. Since A is invertible, it is a zero of F itself, not merely of a preconditioned equation.

## 8. Geometry and identification of the physical solution

The Sonin bound and great-circle Bernstein inequality give ||Hess_S G_j||<=4j^2. The verifier bounds Cgeo=sup_j 4j^2/[(1+j)2^j Omega(2j)] by a finite maximum plus the decreasing tail 4j/2^j. On the entire E1 ball of radius 10^-8 it checks positive radius, Hessian norm below minimum radius, chart perturbation norm below one, and a nonzero j=2 coefficient. The radial second fundamental form is (1+f)^2 g+2df tensor df-(1+f)Hess f, so it is positive definite in every direction, including both transverse block families and the axes. The radial graph is an embedding; the classical global convexity theorem gives a strictly convex body. Group invariance forces any ball center to be the origin, and the j=2 witness excludes that ball.

Exponential angular coefficients give an analytic boundary, including the block axes. Harmonic extension is analytic through the closed ball. For a general zero (g,f), polynomial source truncations converge in C and hence in invariant L2. The Dirichlet Poisson solution K^D g_n converges to K^D g in H2 of the ball. Angular shape truncations converge through two derivatives on the closed ball: the factor 2^j dominates every fixed derivative power of j. Thus the coefficients of the pulled-back differential operator converge uniformly. For polynomial sources and shapes the coefficient formula equals the pointwise chain-rule operator, as the exact identities show. Passing to the limit identifies F=0 with (Delta_Phi+R^2)U=0 almost everywhere, U=1+K^D g in H2. The Dirichlet trace is one, while the normal trace of K^D g is sum_j g_(j,0)G_j/[2(2j+p)]=0 by the Sobolev trace theorem. Tangential derivatives of its constant boundary trace vanish, so its full gradient trace is zero. The analytic chart is a diffeomorphism of the closed ball: |d-1|<=2||f||E1<1 makes each radial map increasing, and the verified Lipschitz perturbation bound also gives injectivity. Transporting U by Phi^-1 gives an H2 solution with both traces on an analytic boundary; analytic Dirichlet elliptic boundary regularity upgrades it to a classical real analytic field. Scaling by R gives (Delta+1)u=0 with u=1 and partial_nu u=0. Harmonic shape extension is analytic for |y|<sqrt(2), so this includes the block axes. The field cannot be constant because its boundary value is one and the Helmholtz parameter is positive.

## 9. Trust boundary and independent ball check

The complete classical prerequisites are: the invariant-polynomial theorem for O(a)xO(b), a,b>=2, and Fischer harmonic decomposition; Jacobi orthogonality, hypergeometric coefficients and endpoint formula; the Jacobi endpoint maximum theorem (equivalently the stated Sonin argument), with B>=A; Gasper positive linearization in the explicitly checked parameter region; Jacobi contiguous and Bessel recurrence identities; simplicity and spacing greater than pi of positive J_nu zeros for nu>1/2, and positivity before the first zero for orders excluded in the window; uniqueness and H2 regularity of the Dirichlet Poisson problem on a ball and the Sobolev trace theorem; harmonic and analytic elliptic Dirichlet regularity through an analytic boundary; analytic polynomial approximation for entire boundary traces; Cauchy estimates for analytic maps between complex Banach spaces and the Banach fixed-point theorem; great-circle Bernstein inequality and the global convexity theorem for a positively curved embedded radial sphere in dimension at least three; and, for the additional invariant L2 corollary only, Green's identity, compact Dirichlet resolvent, and invariant-sector min-max and domain monotonicity. All parameter hypotheses used here hold at p=73 and split (54,92). The general-index Poisson, Jacobi and pull-back identities have the coefficient/chain-rule proofs above. The remaining estimates are explicit above and in the embedded arithmetic/clamping details. This is not a proof-assistant formalization or a claim of community acceptance.

The independent ball check enumerates every invariant Dirichlet state in its rational window: one angular mode per even degree, orders nu>=upper edge excluded by positivity before nu, and a mesh of width <=1/4 together with zero spacing >=pi proves completeness. It is not used as a surrogate for the deformed coefficient inverse. Neither the asymptotic theorem nor any old producer PASS is a dependency.


## 10. Quantitative invariant-sector L2 gap (additional corollary)

Let D0=Phi_f0(B) and mu=R^2. The verified C inverse has norm B0. Since J_D(lambda)=J_D(mu)+(lambda-mu)K^D and ||K^D||C<=3/[p(p+1)], the verifier checks a Neumann margin below 1/2 throughout |lambda-mu|<=1/4000. Thus the coefficient Dirichlet problem is solvable for every source in C throughout this closed interval.

This excludes invariant L2 Dirichlet eigenfunctions in the same interval without a pointwise norm transfer. If v were such an eigenfunction on D0, choose a source f in C with nonzero physical L2 pairing with v after pushforward. Such a source exists: invariant polynomials are dense in invariant L2, and the chart Jacobian is positive and bounded above and below. The coefficient inverse gives W=K^D J_D(lambda)^-1 f in H2 cap H0, and its pushforward solves (Delta+lambda)W=f. Green's identity against v gives a zero pairing, a contradiction. The invariant Dirichlet Laplacian is self-adjoint with compact resolvent, so absence of these eigenvalues gives a center spectral gap at least 1/4000.

For any shape in the E1 ball of exact radius r=1e-46, its radial function differs from the center by at most r. If m0 is the certified positive lower bound for the center radius, epsilon=r/m0 gives (1-epsilon)D0 subset Df subset (1+epsilon)D0. These domains are all group invariant. Zero extension preserves the invariant H0^1 subspaces, so Dirichlet min-max in that sector bounds every ordered eigenvalue by lambda_j(D0)/(1+epsilon)^2 and lambda_j(D0)/(1-epsilon)^2. The verifier checks that these bounds leave a gap at least 1/5000 around mu. After scaling by R, the physical gap around eigenvalue 1 is at least 1/(5000 R^2), and the invariant L2 Dirichlet inverse norm is at most 5000 R^2.

This assertion is ONLY in the O(54)xO(92)-invariant sector. No inverse on the full L2 space is claimed; translation derivatives of a Schiffer field lie outside this sector. This corollary adds Green's identity, compact Dirichlet resolvent and the invariant min-max principle to the explicitly listed classical facts.


# Finite arithmetic lemma for the frozen-table builder

This is the arithmetic lemma used by the frozen-table verifier. No complete inverse follows from this lemma alone.

The model consists of exact binary64 polynomial coefficients, exact binary64 angular product proposals, and exact binary64 input weights. `verify_inf3_hardened.py components` independently encloses the actual angular products and exact weights with Arb. The matrix builder must load those literal proposals; this error budget does not authorize an arbitrary matrix array merely because its dimensions match.

For a fixed angular input j, multiplication by S_i t^h has nonnegative angular coefficients c_(i,j,k) and nonnegative radial parameter conversions. The conversion raises beta by 2(k-j)+h+i+j-k, then lowers h+i+j-k times; total steps equal 2(i+h). Each step is bidiagonal, its exact entries are nonnegative and each column has mass one. Each computed entry uses at most two products and one sum after one division to form a rational stencil entry. Allowing eight floating operations per step therefore encloses the conversion entrywise by the standard positive-product gamma bound. Its absolute weighted column majorant is rho^i Omega(2i+2h). This remains true for moment (1+D/p)^k after multiplying by (1+(2i+2h)/p)^k.

The K, K', and beta K+tK' columns have at most three entries. Their rational denominators and the integer gradient factors are formed exactly before binary64 division; the verifier checks a bound on every involved integer expression below 2^53. The explicit gradient-EK formula uses at most three such converted columns per angular output. For nonharmonic source columns their absolute sums are bounded by A+4A+2A=7A, A=2i. The saturated output is grouped before conversion and has bound 2A. For harmonic source columns the bounds are A+2A+A=4A (saturated: <=3A/2). Thus 7A is valid for all columns, including s=0. These are bounds on the sum of absolute contributions, so subsequent cancellation cannot invalidate the rounding majorant.

Similarly the absolute computations for K, EK and E(E+1)K have global column bounds

`K <= Omega(2)(1+2/p)^k /[(p-1)(p+1)]`,
`EK <= 2 Omega(2)(1+2/p)^k/p`,
`E(E+1)K <= Omega(2)(1+2/p)^k [2+(2p-3)/p]`.

Here EK is evaluated as B K+2tK', and E(E+1)K as B(B+1)K+(6-4p)tK'+tI. These displayed bounds cover the absolute values of the separate evaluated terms.

Let L be the largest stored polynomial degree and M the largest number of stored polynomial terms. For each output entry, accumulation over polynomial terms has at most M contributions. A nontrivial matrix product has inner dimension at most W=(J+1)N and hence at most W products and W accumulations per entry. Each path through the recipe contains the radial conversions, their short K products, one polynomial-times-derivative product, additions of the five operator terms, and two diagonal scalings. The deliberately enlarged count

`Kops = 8W +32(L+M+shape_order+4)+1024`

exceeds the number of rounding operations along every such path, including sparse duplicate coalescing. Zero additions and index sorting add no numerical error. Every matrix stage checks finiteness and canonical sparse format. The angular export is used at most twice along a path (multiplier and gradient); error e_a therefore contributes at most 4 e_a times the same absolute majorant when e_a<1/10.

The exact unweighted radial conversion columns preserve mass. Grouping every absolute contribution in the five cleared-operator terms therefore gives the majorant

`Tabs = ||d^4|| +2||d^3|| Gabs +||m d^2|| EEKabs +||P0|| EKabs +mu ||q^2 d^4|| Kabs`.

All polynomial norms and Gabs=14 sum_i i|q_i| rho^i Omega(2i)(1+2i/p)^k are computed with Arb from exact dyadic coefficients. This majorant bounds both the exact operator and the absolute computation used for the forward-error estimate; it is not a measured condition number.

If each exported coefficient weight differs relatively from its exact value by at most e_w, diagonal conjugation introduces at most 2e_w/(1-e_w)*Tabs. It enlarges the preceding arithmetic error by (1+e_w)/(1-e_w). Thus the full finite matrix error in the exact weighted norm is at most

`[((gamma_Kops+4e_a)(1+e_w)/(1-e_w))+2e_w/(1-e_w)] Tabs +1e-240`.

The error is uniform over all retained input columns, independent of their number. The underflow allowance dominates the binary64 absolute underflow error: W<=30000, weights in [1,1e40], auxiliary integer weights below 2^53, and the stated path count keep the total weighted underflow many orders below 1e-240. The final finite product verifier separately bounds rounding in C*M_PP-I and C*M_PQ; its error must be added to ||C|| times the assembly error.

The input rectangle must contain every tail column that can reach the head and every intermediate row that can return to it. These are explicit support gates in the uniform-tail verifier. Parameter raises only decrease radial index and are performed before the lowering/multiplication steps, which only increase it. Hence truncating an individual conversion matrix cannot remove a path that later re-enters its retained rows. The composite K/gradient and final multiplier support gates additionally cover their intermediate cutoffs.

The argument assumes IEEE binary64 basic arithmetic; the unit roundoff used is 2^-52, twice the usual round-to-nearest value. Arb evaluates gamma and all final inequalities. Library eigenvalue accuracy, LU stability, and sampled PDE checks are not assumptions of the inverse proof: the candidate C is an arbitrary dyadic matrix, whose full residual and couplings are explicitly checked.


# Explicit clamping and material estimates

These analytic estimates are used by the interval verifier. Numerical success is determined only by its all-links receipt. The general-index smoothing and gradient proofs are in VERIFIER_PROOF.md.

Use boundary moments B_k(v)=sum_j rho^j Omega_p(2j)(1+2j)^k |v_j| and corresponding source moments C^[k]. Let D,N be the physical Dirichlet and radial-Neumann defects. Let q=1+P f, d=q+E q. The reference normal defect is Nref=R d_boundary N. Because (1+D1+D2)<= (1+D1)(1+D2), fixed moment weights are submultiplicative as well as Omega.

For the exact lift ell=sum [D_j S_j(1-j(t-1)) + Nref_j S_j(t-1)/2], put U0=uapp(R Phi_f)-ell. Then U0 is exactly clamped and g0=Delta(U0-1) has zero harmonic source layer.

A conservative collection of lift bounds, with omega=Omega_p(2), is

- ||ell||C[k] <= omega*3^k*(3 B_k(D)+B_k(Nref)/p), from its two Jacobi layers.
- ||Delta ell||C[k] <=2(p+1)B_(k+2)(D)+2p B_(k+1)(Nref), from the exact harmonic Laplacians.
- ||E ell||C[k] <=2 omega*3^k*(B_(k+2)(D)+B_(k+1)(Nref)), by the raw radial monomials in the two lifts.
- ||E(E+1)ell||C[k] <=6 omega*3^k*(B_(k+3)(D)+B_(k+2)(Nref)).
- ||grad(P f) dot grad E ell||C[k] <=8(p+1)omega*3^k*F1_k*(B_(k+3)(D)+B_(k+2)(Nref)), where F1_k=sum(1+j)rho^j Omega(2j)(1+2j)^k |f_j|.

The last estimate follows from the exact solid-gradient product: for output index l and h=i+j-l, its coefficient on S_l t^(h+h0-1) is 2h(h+2l+p-1)+4i h0 when the input is S_j t^h0, h0=0 or 1. It is <=4i(2j+p+1). Euler differentiation adds at most 2j+2, and the lift monomial coefficients add at most (1+2j)|D_j|+|Nref_j|. Every radial monomial-to-Jacobi conversion is positive and sums to one.

All fixed multipliers in the scalar pullback can be bounded in moment k. For example M_m<=4(p+1)F1_k^2, M_(grad q dot grad E q)<=4(p+1)F1_k*F_E1_k. q^-1 and d^-1 can be bounded by geometric series if their moment-k perturbation norms are <1. If that test fails, it must be replaced by a moment-aware sum of the base-norm series, not silently discarded. Together with the lift bounds this gives E=-B_f ell in C and C^[2].

## Avoiding pointwise estimates of the known center

If a full J_D inverse is already certified in C and C^[2], then the exact identity

`g0=J_D^-1(E-R^2)`

bounds G0=||g0||C and G2=||g0||C[2]. This identity is valid for the exactly clamped center and cannot be applied to the unlifted collocation field.

Let v=d^-1 and f0=v E U0. d is harmonic because q and E q are harmonic. Thus

`Delta f0 = d^-1(E+2)g0 -2d^-2 grad d dot grad E K g0 +2d^-3 |grad d|^2 E K g0`.

The Euler column bound is <=(p+1)(1+2j+2s)^2, so ||E g0||C <=(p+1)G2. With delta_d=||d-1||E1 and the grouped gradient estimates of the embedded operator identities, a valid bound is

`Bf = ||d^-1||*((p+1)G2+2G0) +28||d^-1||^2 delta_d G0 +32(p+1)/p ||d^-1||^3 delta_d^2 G0`.

This explicitly exposes the need for the source moment G2. An unweighted inverse norm alone is insufficient.

Since f0 has zero Dirichlet trace, f0=K^D Delta f0. The exact product identity then yields

`||Delta[(P h)f0]|| <=(1+32/p) Bf ||h||E1`.

For the material inverse, certify the boundary E1 norm of d/D_p, where D_p=(R^2-E_boundary)/(nu^T A nu). The numerator of D_p^-1 is a fixed-center metric multiplier; the reciprocal of R^2-E_boundary needs an explicit E1 algebra margin. Put H=(1/2)||d/D_p||E1, using the Poisson normal-trace bound. Then

`||Q0^-1|| <=1+H*(1+(1+32/p)Bf)`.

The material correction has the form (P h)/d times Eul E, so its norm is at most ||d^-1||*(p+1)||E||C[2]. This leads to a concrete Z0 once ||A|| is known.

## Numerical gates

The graph inverse in both C and C^[2], the full-center rather than quartic-chart perturbation, the fixed-multiplier arithmetic, explicit D2F bound on a common X-ball, and the final radii/geometry gates. The verifier must discharge these gates; the formulas alone are not an existence certificate.

## Membership before the inverse identity

The inverse identity cannot be used to infer membership circularly. At a fixed p the finite regular Bessel field has a power series in the squared Euclidean radius with factorial-squared denominators. Composition with the finite polynomial chart and multiplication by the finite solid harmonics converge absolutely in every fixed coefficient norm and every fixed degree moment: positive monomial-to-Jacobi multiplication bounds grow only geometrically, while the radial series has factorial decay. The boundary traces are entire functions of x, hence their Jacobi expansions and the exact lift have finite norms at every fixed angular radius and moment. Consequently the explicitly defined g0 already belongs to C^[k] for each fixed k, before invoking J_D^-1 to obtain a useful numerical bound. This is a membership argument; a pointwise-to-coefficient norm estimate is not substituted for it.
