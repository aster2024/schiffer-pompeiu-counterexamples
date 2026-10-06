# G_u(mu) = mu^{3/2} max_theta |R_D(mu,0,theta)| for the unit-conformal-radius domain (t=0)
import numpy as np, json, sys, time
from geom import cu, maps, m
from defect import at_normal_angle
def pow2_ge(x):
    n=1
    while n<x: n*=2
    return n
def ID(mu, thetas, N):
    s = 2*np.pi*np.arange(N)/N
    w, psi, dpsi, ddpsi = maps(cu, s)
    jac = np.abs(dpsi)
    out = np.zeros(len(thetas), complex)
    for k,th in enumerate(thetas):
        ph = mu*(psi.real*np.cos(th)+psi.imag*np.sin(th))
        out[k] = np.sum(np.exp(1j*ph)*jac)*(2*np.pi/N)
    return out
def lead(mu, thetas):
    h0,hp0,r0,_ = at_normal_angle(thetas)
    h1,hp1,r1,_ = at_normal_angle(thetas+np.pi)
    return np.sqrt(2*np.pi/mu)*(np.sqrt(r0)*np.exp(1j*mu*h0-1j*np.pi/4) + np.sqrt(r1)*np.exp(-1j*mu*h1+1j*np.pi/4))
if __name__ == '__main__':
    nth = int(sys.argv[1]) if len(sys.argv)>1 else 200
    mus = [float(x) for x in sys.argv[2].split(',')] if len(sys.argv)>2 else [3,10,30,100,238,500,1000,3000,1e4,3e4,1e5]
    thetas = np.linspace(0, np.pi/m, nth)
    res = {}
    for mu in mus:
        t0=time.time()
        N = pow2_ge(max(4096, 3*mu+2000))
        I = ID(mu, thetas, N)
        L = lead(mu, thetas)
        R = I - L
        G = mu**1.5*np.abs(R)
        k = int(np.argmax(G))
        # convergence check with 2N for the argmax theta only
        I2 = ID(mu, thetas[k:k+1], 2*N)
        err = abs(I2[0]-I[k])
        res[mu] = dict(N=N, G=float(G.max()), theta=float(thetas[k]), absR=float(abs(R[k])), absLead=float(abs(L[k])), absI=float(abs(I[k])), Nerr=float(err), Gmin=float(G.min()))
        print(mu, res[mu], 'time %.1f'%(time.time()-t0), flush=True)
    json.dump(res, open('remainder_%d.json'%nth,'w'), indent=1)
