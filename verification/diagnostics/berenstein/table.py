import numpy as np, ast, json, re, os
from defect import at_normal_angle
from geom import m, rho_c
alpha_u = rho_c**2
th = np.linspace(0, 2*np.pi/m, 4001)
h0,hp0,r0,_ = at_normal_angle(th); h1,hp1,r1,_ = at_normal_angle(th+np.pi)
def Bu(s):
    return max(np.max(np.abs(np.sqrt(r0)*np.exp(t*s*hp0)-np.sqrt(r1)*np.exp(-t*s*hp1))) for t in (-1,1))
def T(s): return np.sqrt(2*np.pi)*np.sqrt(alpha_u+s*s)*Bu(s)
G = {}
for line in open('remainder_160.log'):
    mm = re.match(r'^([0-9.e+]+) (\{.*\}) time', line)
    if mm:
        d = ast.literal_eval(mm.group(2)); G[float(mm.group(1))] = d['G']
src = 'remainder_160.log'
if os.path.exists('remainder_scan.json'):
    for mu,g in json.load(open('remainder_scan.json')): G[mu]=max(G.get(mu,0),g)
    src += ' + remainder_scan.json'
if os.path.exists('remainder_32_highmu.json'):
    for mu,d in json.load(open('remainder_32_highmu.json')).items(): G[float(mu)]=max(G.get(float(mu),0),d['G'])
    src += ' + remainder_32_highmu.json'
mus = np.array(sorted(G)); gs = np.array([G[x] for x in mus])
print('G_u(mu) sources:', src, ' (%d mu values, mu from %g to %g)'%(len(mus),mus[0],mus[-1]))
print('overall max G found = %.1f at mu = %g'%(gs.max(), mus[np.argmax(gs)]))
smax = rho_c/np.sqrt(8)
print('s_max (alpha_s >= 8) = rho_c/sqrt(8) = %.4f'%smax)
print('   s        alpha_s=alpha_u/s^2   B_u(s)     T(s)     S(s)=max_{mu>=3s} G    S/T')
for s in [1e-4,1e-2,0.1,1,3,10,20,30,45,60,75,smax]:
    sel = mus>=3*s
    S = gs[sel].max()
    print('%9.4f  %14.6g  %.6f  %8.3f  %10.1f   %6.1f'%(s, alpha_u/s**2, Bu(s), T(s), S, S/T(s)))
# worst case over a fine grid of s
ss = np.linspace(1e-4, smax, 400)
ratio = [gs[mus>=3*s].max()/T(s) for s in ss]
print('min over s in (0,s_max] of S(s)/T(s) = %.2f at s=%.3f' % (min(ratio), ss[int(np.argmin(ratio))]))
# what the claim itself forces at the claimed scale s=1 (unit conformal radius, lambda=sqrt(alpha_u+1))
print('T(1) = sqrt(2pi) lambda B = %.4f  (claim itself forces C_* >= this: direct reading)'%T(1.0))
