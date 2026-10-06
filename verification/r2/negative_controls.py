"""Fault injections for added gates, with positive controls; safe under python -O."""
import copy
import json
from flint import arb, ctx
from spectral_interval import need, exact_decimal, umax
from analytic_cap import field_bounds, high_sector_gap
from coverage_checks import check_counts, check_root_coverage, check_scan_binding, check_scan_record, certified_sign
from final_gates import gate, check_public_bounds, UPPER_LIMITS, LOWER_LIMITS, SPECTRAL_UPPER, RADIUS
ctx.threads=1

def negative_tests(spectral,boundary,cap,candidate,reference,scan,work):
    ctx.prec=256
    rejected=[]
    def reject(name,fn,reason=None):
        try:
            fn()
        except (RuntimeError,ValueError,KeyError,TypeError) as exc:
            if reason:
                need(reason in str(exc),'wrong rejection for '+name+': '+str(exc))
            rejected.append({'test':name,'reason':str(exc)[:240]})
        else:
            raise RuntimeError('negative test accepted: '+name)
    base=[spectral,boundary,cap,candidate,reference,scan]
    def mutate(name,owner,key,value,delete=False,reason=None):
        args=copy.deepcopy(base);obj=args[owner]
        keys=key if isinstance(key,tuple) else (key,)
        for part in keys[:-1]: obj=obj[part]
        if delete: del obj[keys[-1]]
        else: obj[keys[-1]]=value
        reject(name,lambda:gate(*args),reason)
    gate(*base)
    for t in ('nan','inf','-inf'):
        reject('maximum_'+t,lambda v=t:umax([arb(0),arb(v)]))
    for t in ('nan','inf','-inf','+inf','','[1 +/- 0.01]','1/3',1.0):
        reject('decimal_'+repr(t),lambda v=t:exact_decimal(v),'exact decimal')
        for key in ('p','a','scales'):
            v=copy.deepcopy(candidate);v[key][0]=t
            reject('field_'+key+'_'+repr(t),lambda x=v:field_bounds(x),'exact decimal')
    for t in ('0','-1'):
        v=copy.deepcopy(candidate);v['scales'][0]=t
        reject('scale_'+t,lambda x=v:field_bounds(x),'invalid finite field')
    v=copy.deepcopy(candidate);v['a'].pop()
    reject('field_count',lambda:field_bounds(v),'invalid finite field')
    v=copy.deepcopy(candidate);v['a']=['0']*len(v['a'])
    reject('strip_constant_not_absorbed',lambda:field_bounds(v),'strip constant')
    for rho in (arb(910),arb(911),arb('910 +/- 1')):
        reject('d_h_'+str(rho),lambda x=rho:high_sector_gap(5,182,x),'d_h')
    cases=[
        ('missing_E1',1,'E1_upper',None,True),('nan_E0',1,'E0_upper','nan',False),
        ('bad_shape_hash',1,'shape_sha256','0'*64,False),
        ('missing_schur',0,'schur_relative_bound_upper',None,True),
        ('bad_schur',0,'schur_relative_bound_upper','2',False),
        ('inflated_inverse_not_recoupled',0,'full_Dirichlet_inverse_L2_norm_upper','1e6',False),
        ('inflated_E0_not_recoupled',1,'E0_upper','1',False),
        ('inflated_E1_not_recoupled',1,'E1_upper','1',False),
        ('nan_Y',2,('bounds','Y'),'nan',False),('missing_C3',2,('bounds','C3'),None,True),
        ('forged_Y_zero',2,('bounds','Y'),'0',False),
        ('forged_convexity',2,('bounds','curvature_lower'),'1',False),
        ('wrong_radius',2,('bounds','radius'),'1e-23',False),
        ('states_407',0,'states',407,False),('wrong_window',0,'window',159,False),
        ('wrong_gauss_count',0,'gauss_nodes',1535,False),('grid_8',1,'grid',8,False),
        ('wrong_lmax',1,'lmax',399,False),('wrong_modes',1,'modes',30,False),
        ('missing_derivative_tail',1,'N_moments_upper',[],False)]
    for args in cases: mutate(*args)
    for j in range(6):
        counts=list(spectral['sector_counts']);counts[j]-=1;counts[(j+1)%6]+=1
        mutate('sector_'+str(j)+'_total_preserved',0,'sector_counts',counts,reason='sector count')
    for j in range(7):
        counts=list(spectral['independent_sector_counts']);counts[j]+=1
        mutate('independent_sector_'+str(j),0,'independent_sector_counts',counts,reason='sector count')
    reject('sector_dimensions',lambda:check_counts([408]),'dimensions')
    b={k:arb(v) for k,v in cap['bounds'].items()}
    s={k:arb(v) for k,v in spectral.items() if k in SPECTRAL_UPPER or k=='complement_gap_lower'}
    check_public_bounds(s,b,boundary)
    for key,t in UPPER_LIMITS.items():
        v=dict(b);v[key]=arb(t)
        reject('public_upper_'+key,lambda x=v:check_public_bounds(s,x,boundary),'public upper: '+key)
    for key,t in LOWER_LIMITS.items():
        v=dict(b);v[key]=arb(t)
        reject('public_lower_'+key,lambda x=v:check_public_bounds(s,x,boundary),'public lower: '+key)
    for key,t in SPECTRAL_UPPER.items():
        v=dict(s);v[key]=arb(t)
        reject('public_spectral_'+key,lambda x=v:check_public_bounds(x,b,boundary),'public spectral upper: '+key)
    v=dict(b);v['radii_polynomial']=arb('-9.99999993e-33')
    reject('public_radii_polynomial',lambda:check_public_bounds(s,v,boundary),'public radii polynomial')
    v=dict(s);v['complement_gap_lower']=arb('230022.34')
    reject('public_complement',lambda:check_public_bounds(v,b,boundary),'public complement gap')
    for key,t in [('E0_upper','9.643e-61'),('E1_upper','5.418e-56')]:
        v=dict(boundary);v[key]=t
        reject('public_'+key,lambda x=v:check_public_bounds(s,b,x),'public E')
    data=json.loads((work/'spectral_inputs.json').read_text())
    check_scan_binding(scan,data)
    def bad_roots(name,edit,reason=None):
        value=copy.deepcopy(data);edit(value)
        reject(name,lambda:check_root_coverage(value),reason)
    root_cases=[
        ('input_stage',lambda x:x.update(stage='BAD'),'coverage input stage'),
        ('angular_cutoff',lambda x:x['state']['p'].__setitem__(0,'1000'),'coverage angular cutoff'),
        ('missing_sector',lambda x:x['coverage'].pop(),'missing coverage sector'),
        ('wrong_sector',lambda x:x['coverage'][0].update(j=1),'coverage sector identity'),
        ('missing_guards',lambda x:x['coverage'][0].update(roots=[]),'missing root guards'),
        ('wrong_order',lambda x:x['coverage'][0]['roots'][0].update(n=182),'root angular order'),
        ('wide_bracket',lambda x:x['coverage'][0]['roots'][0].update(right='900'),'bracket geometry'),
        ('same_signs',lambda x:x['coverage'][0]['roots'][0].update(left='600',right='600.01'),'endpoint signs'),
        ('bad_norm',lambda x:x['coverage'][0]['roots'][0].update(norm='0'),'normalization'),
        ('bad_first_flag',lambda x:x['coverage'][0].update(starts_at_first_root='yes'),'flag type'),
        ('order_zero_first',lambda x:x['coverage'][0].update(starts_at_first_root=True),'order zero needs lower guard'),
        ('lost_lower_guard',lambda x:x['coverage'][0]['roots'].pop(0),'lower guard'),
        ('lost_upper_guard',lambda x:x['coverage'][0]['roots'].pop(),'upper guard'),
        ('skipped_root',lambda x:x['coverage'][0]['roots'].pop(50),'skipped or duplicate'),
        ('lost_first_zero',lambda x:x['coverage'][4]['roots'].pop(0),'first-root coverage'),
        ('false_first_guard',lambda x:x['coverage'][4].update(starts_at_first_root=False),'lower guard'),
        ('wrong_window',lambda x:x.update(window=159),'window contract'),
        ('forged_counts',lambda x:x.update(sector_counts=[408]*6),'stored coverage counts'),
        ('state_replaced_same_counts',lambda x:x['states'].__setitem__(151,copy.deepcopy(x['states'][150])),'states differ')]
    for name,edit,reason in root_cases: bad_roots('coverage_'+name,edit,reason)
    for name,key,val in [
        ('bad_stage','stage','WRONG'),('wrong_hash','reference_shape_sha256','0'*64),
        ('wrong_grid','grid_cells',1279),('wrong_step','step','0.5'),
        ('wrong_total','total',407),('missing_sector','sectors',scan['sectors'][:-1])]:
        v=copy.deepcopy(scan);v[key]=val
        reject('scan_'+name,lambda x=v:check_scan_record(x,reference))
    for name,edit,reason in [
        ('zero_sign',lambda x:x['sectors'][0]['signs'].__setitem__(1,0),'nonzero signs'),
        ('missing_sign',lambda x:x['sectors'][0]['signs'].pop(),'nonzero signs'),
        ('wrong_sector_id',lambda x:x['sectors'][0].update(n=182),'sector identity'),
        ('wrong_change_cell',lambda x:x['sectors'][0]['change_cells'].__setitem__(0,0),'sign-change count')]:
        v=copy.deepcopy(scan);edit(v)
        reject('scan_'+name,lambda x=v:check_scan_record(x,reference),reason)
    v=copy.deepcopy(data);v['states'][0].update(left='900',right='900.01')
    reject('root_cell_mismatch',lambda:check_scan_binding(scan,v),'root-cell mismatch')
    import coverage_checks
    original=coverage_checks.exact_decimal
    try:
        coverage_checks.exact_decimal=lambda _:arb('nan')
        reject('unresolved_sign_oracle',lambda:certified_sign(0,'814',0),'nonfinite')
    finally:
        coverage_checks.exact_decimal=original
    original=coverage_checks.exact_decimal
    try:
        coverage_checks.exact_decimal=lambda _:arb(160)
        reject('zero_sign_oracle',lambda:certified_sign(182,'814',0),'unresolved independent sign')
    finally:
        coverage_checks.exact_decimal=original
    from coverage_checks import sign_scan
    for label,key,value,reason in [('symmetry','m',181,'sign scan symmetry'),
                                   ('spacing','p',['100'],'sign scan zero spacing'),
                                   ('high_order','p',['1000'],'higher-order exclusion')]:
        v=copy.deepcopy(reference);v[key]=value
        reject('scan_'+label+'_guard',lambda x=v:sign_scan(x),reason)
    from spectral_interval import shape_digest
    from boundary_newton import centre_digest
    args=copy.deepcopy(base);args[3]['p'][0]='800'
    args[1]['shape_sha256']=shape_digest(args[3]);args[1]['center_sha256']=centre_digest(args[3])
    reject('public_conformal_radius',lambda:gate(*args),'public conformal radius')
    from scalar_second_iv import verify as iv_verify
    iv_verify(candidate,reference,spectral,boundary,RADIUS)
    v=copy.deepcopy(boundary);v['D_moments_upper'][0]='1'
    reject('iv_corrupt_moment',lambda:iv_verify(candidate,reference,spectral,v,RADIUS))
    v=copy.deepcopy(candidate);v['scales'][0]='inf'
    reject('iv_infinite_scale',lambda:iv_verify(v,reference,spectral,boundary,RADIUS),'exact decimal')
    reject('iv_bad_radius',lambda:iv_verify(candidate,reference,spectral,boundary,'1e-23'))
    return {'status':'FAIL_CLOSED_TESTS_PASS','positive_controls':4,'count':len(rejected),'rejected':rejected}
