import numpy as np
phi = [1.50795385361929982e-2,5.76689422728655252e-3,1.01557038850826134e-3,2.71368762415077922e-4,5.94916615455076681e-5,1.40376217202844726e-5,3.19906904482506334e-6,7.39433442683568780e-7,1.69810772559526110e-7,3.90984465868492322e-8,8.99107926650920794e-9,2.06853696471111575e-9,4.75778692314195368e-10,1.09440420794361152e-10,2.51724454756957816e-11,5.78995922871089628e-12,1.33169815228079109e-12,3.06233427501141624e-13,7.09000044668921563e-14,1.71384604459278807e-14]
m=13
N=1<<16
s=2*np.pi*np.arange(N)/N
w=np.exp(1j*s)
psi=w.copy(); d=np.ones_like(w); dd=np.zeros_like(w)
for j,c in enumerate(phi,1):
    e=m*j+1
    psi+=c*w**e; d+=c*e*w**(e-1); dd+=c*e*(e-1)*w**(e-2)
cf=(1+w*dd/d).real
kap=cf/np.abs(d)
print('CSB (Berenstein, D_13) centre: min/max Re(1+w phi\'\'/phi\') =',cf.min(),cf.max())
print('curvature min/max',kap.min(),kap.max(),' radial range',np.abs(psi).min(),np.abs(psi).max())
print('number of sign changes of curvature around the boundary:', int(np.sum(np.diff(np.sign(cf))!=0)))
