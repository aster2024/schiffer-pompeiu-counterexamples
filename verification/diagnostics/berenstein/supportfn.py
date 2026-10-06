import numpy as np
from defect import at_normal_angle
from geom import m, rho_c
N = 1<<13
th = 2*np.pi*np.arange(N)/N
h,hp,rho,_ = at_normal_angle(th)
H = np.fft.fft(h)/N
k = np.fft.fftfreq(N, 1.0/N)
print('unit conformal radius normalisation: h in [%.8f, %.8f], mean %.8f'%(h.min(),h.max(),h.mean()))
print('harmonics of h (angular frequency, |amp|):')
for kk in range(0, 77*8+1, 77):
    print('   k=%4d  |h_k| = %.4e'%(kk, 2*abs(H[kk]) if kk>0 else abs(H[0])))
D = {}
for j in range(0,7):
    dj = np.fft.ifft(H*(1j*k)**j).real*N
    D[j]=dj
    print('sup|h^(%d)| = %.6e'%(j, np.abs(dj).max()))
M = max(np.abs(D[j]).max() for j in range(0,6))
print('nested C^5 norm M = max_j<=5 sup|h^(j)| =', M)
print('check: h+h\'\' == rho: max dev', np.max(np.abs(D[0]+D[2]-rho)))
print('radius of curvature rho range [%.8f, %.8f]; curvature [%.8f, %.8f]'%(rho.min(),rho.max(),1/rho.max(),1/rho.min()))
# physical (eigenvalue 1) scale
sc = rho_c
print('physical scale (eigenvalue 1): factor %.6f'%sc)
print('  curvature range [%.6e, %.6e]; radius of curvature [%.4f, %.4f]; |h-mean| max %.6f; max|h\'| = %.6f'%(1/rho.max()/sc,1/rho.min()/sc, rho.min()*sc, rho.max()*sc, (h-h.mean()).__abs__().max()*sc, np.abs(hp).max()*sc))
# width function w(theta)=h(theta)+h(theta+pi)
hpi = np.roll(h, N//2)     # theta+pi
w = h+hpi
print('width: min %.8f max %.8f  ratio %.8f  (physical: min %.5f max %.5f)'%(w.min(), w.max(), w.min()/w.max(), w.min()*sc, w.max()*sc))
print('Hausdorff dist to centred discs (max h - min h)/2 = %.4e  (relative, unit) ; physical %.5f'%((h.max()-h.min())/2, (h.max()-h.min())/2*sc))
# odd harmonics fraction: h(theta)-h(theta+pi)
print('max |h(theta)-h(theta+pi)| = %.4e ; max |rho(theta)-rho(theta+pi)| = %.4e'%(np.abs(h-hpi).max(), np.abs(rho-np.roll(rho,N//2)).max()))
# DSWZ class constants (unit scale): R, M, kappa0, kappa1
print('DSWZ class (unit scale): R=%.8f  M=%.6e  kappa0=%.6f kappa1=%.6f'%(1.0000137318, M, 1/rho.max(), 1/rho.min()))
