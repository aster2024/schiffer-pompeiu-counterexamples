"""Recompute final scalar bounds in-process and enforce all published thresholds."""
from flint import arb, ctx
from spectral_interval import need, upper, exact_decimal, shape_digest
from boundary_newton import centre_digest
from analytic_cap import bounds
from coverage_checks import check_counts, check_scan_record
ctx.threads = 1
RADIUS = '1e-32'
UPPER_LIMITS = {
    'B_L2': '6.748', 'high_angular_beta': '0.062', 'B_Y': '1.045e8',
    'B_trace': '7.039e13', 'U_field': '131.610', 'dU_field': '110996',
    'A': '6.419e19', 'C2': '3.479e22', 'C3': '1.605e19',
    'four_C2_Y': '8.611e-18', 'contraction': '6.957e-10', 'Y': '6.189e-41', 'Z': '2.122e-38'}
LOWER_LIMITS = {
    'curvature_lower': '0.2998879', 'convexity_margin_lower': '243.3264',
    'non_disc_c1_lower': '0.01645', 'analytic_univalence_margin_lower': '808',
    'analytic_radius_lower': '1.0038', 'unit_derivative_lower': '811',
    'high_angular_gap_lower': '164759'}
SPECTRAL_UPPER = {
    'gauss_error_upper': '1.35e-360', 'congruence_gram_defect_upper': '1e-12',
    'relative_finite_congruence_defect_upper': '1e-9',
    'schur_relative_bound_upper': '0.240915',
    'full_Dirichlet_inverse_L2_norm_upper': '6.748'}


def check_public_bounds(s, b, boundary):
    for key, value in UPPER_LIMITS.items():
        need(b[key] < arb(value), 'public upper: '+key)
    for key, value in LOWER_LIMITS.items():
        need(b[key] > arb(value), 'public lower: '+key)
    for key, value in SPECTRAL_UPPER.items():
        need(s[key] < arb(value), 'public spectral upper: '+key)
    need(b['radii_polynomial'] < arb('-9.99999993e-33'), 'public radii polynomial')
    need(s['complement_gap_lower'] > arb('230022.34'), 'public complement gap')
    need(arb(boundary['E0_upper']) < arb('9.643e-61'), 'public E0 upper')
    need(arb(boundary['E1_upper']) < arb('5.418e-56'), 'public E1 upper')


def gate(spectral, boundary, cap, candidate, reference, scan):
    ctx.prec = 256
    need(spectral['stage'] == 'INVARIANT_SECTOR_DIRICHLET_INVERSE', 'wrong spectral stage')
    need(boundary['stage'] == 'CLAMPED_FIELD_RESIDUAL_AND_FIRST_DERIVATIVE', 'wrong boundary stage')
    need(cap['stage'] == 'COUPLED_RADII_GATES'
         and cap['analytic_lemma_version'] == 'fixed_disc_wiener_clamped_v1', 'wrong coupled stage')
    need(candidate['m'] == reference['m'] == spectral['m'] == 182, 'symmetry mismatch')
    need(candidate['modes'] == boundary['modes'] == 40, 'candidate mode mismatch')
    need(spectral['states'] == 408 and spectral['window'] == 160
         and spectral['gauss_nodes'] == 1536, 'spectral dimensions')
    check_counts(spectral['sector_counts'])
    check_counts(spectral['independent_sector_counts'], include_high=True)
    check_scan_record(scan, reference)
    need(spectral['independent_sector_counts'] == scan['counts'], 'scan/spectral counts mismatch')
    need(boundary['grid'] == 1024 and boundary['lmax'] == 400, 'boundary grid contract')
    need(spectral['reference_shape_sha256'] == shape_digest(reference), 'spectral center mismatch')
    need(boundary['center_sha256'] == centre_digest(candidate)
         and boundary['shape_sha256'] == shape_digest(candidate), 'boundary center mismatch')
    need(exact_decimal('814.45') < exact_decimal(candidate['p'][0]) < exact_decimal('814.46'),
         'public conformal radius')
    spec_fields = ('gauss_error_upper', 'congruence_gram_defect_upper',
                   'relative_finite_congruence_defect_upper',
                   'weighted_complement_gram_trace_upper', 'potential_sup_upper',
                   'complement_gap_lower', 'schur_relative_bound_upper',
                   'minimum_abs_proposed_diagonal', 'full_Dirichlet_inverse_L2_norm_upper')
    s = {k: arb(spectral[k]) for k in spec_fields}
    for value in s.values():
        upper(value)
        need(value >= 0, 'negative spectral bound')
    need(s['gauss_error_upper'] < arb(1)/10**20, 'quadrature gate failed')
    need(s['congruence_gram_defect_upper'] < 1
         and s['relative_finite_congruence_defect_upper']+s['schur_relative_bound_upper'] < 1
         and s['complement_gap_lower'] > 0 and s['minimum_abs_proposed_diagonal'] > 0,
         'spectral inverse gate failed')
    for key in ('strip_trace_majorant_upper', 'quadrature_error_per_coefficient',
                'lift_upper', 'dlift_upper', 'lap_lift_upper', 'dlap_lift_upper', 'E0_upper', 'E1_upper'):
        value = arb(boundary[key])
        upper(value)
        need(value >= 0, 'negative boundary bound')
    for key in ('D_moments_upper', 'N_moments_upper', 'each_far_moment_upper'):
        need(len(boundary[key]) == 4, 'missing differentiated boundary tail')
        for text in boundary[key]:
            value = arb(text)
            upper(value)
            need(value >= 0, 'negative boundary moment')
    # The stage file is cross-checked only: every final constant comes from HERE.
    fresh = bounds(candidate, reference, spectral, boundary, RADIUS)
    need(cap == fresh, 'coupled stage differs from in-process recomputation')
    b = {k: arb(v) for k, v in fresh['bounds'].items()}
    for key, value in b.items():
        upper(value)
        if key != 'radii_polynomial':
            need(value >= 0, 'negative coupled bound: '+key)
    r = exact_decimal(RADIUS)
    need(b['radius'].contains(r), 'incorrect existence radius')
    polynomial = upper(b['Y'])+(upper(b['Z'])-1)*r+upper(b['C2'])*r*r+upper(b['C3'])*r*r*r
    derivative = upper(b['Z'])+2*upper(b['C2'])*r+3*upper(b['C3'])*r*r
    need(polynomial < 0 and derivative < 1 and b['radii_polynomial'] < 0
         and b['contraction'] < 1 and b['high_angular_beta'] < 1, 'Newton radii gate failed')
    check_public_bounds(s, b, boundary)
    return {'radii_polynomial_upper': upper(polynomial).str(55),
            'contraction_upper': upper(derivative).str(55),
            'four_C2_Y_upper': upper(4*upper(b['C2'])*upper(b['Y'])).str(55),
            'recomputed_coupled_bounds': fresh['bounds']}
