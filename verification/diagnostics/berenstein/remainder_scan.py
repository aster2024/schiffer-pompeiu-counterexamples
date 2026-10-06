import numpy as np, json, time
from remainder import ID, lead, pow2_ge
from geom import m
thetas = np.linspace(0, np.pi/m, 48)
mus = np.exp(np.linspace(np.log(200), np.log(30000), 140))
out=[]
for mu in mus:
    N = pow2_ge(max(4096, 3*mu+2000))
    I = ID(mu, thetas, N); L = lead(mu, thetas)
    G = mu**1.5*np.abs(I-L)
    out.append((float(mu), float(G.max())))
json.dump(out, open('remainder_scan.json','w'))
a = np.array(out)
print('max G over scan:', a[:,1].max(), 'at mu=', a[np.argmax(a[:,1]),0])
for lo,hi in [(200,500),(500,1000),(1000,3000),(3000,10000),(10000,30000)]:
    sel=(a[:,0]>=lo)&(a[:,0]<hi)
    print('mu in [%d,%d): max G = %.1f, median G = %.1f, min G=%.1f'%(lo,hi,a[sel,1].max(),np.median(a[sel,1]),a[sel,1].min()))
