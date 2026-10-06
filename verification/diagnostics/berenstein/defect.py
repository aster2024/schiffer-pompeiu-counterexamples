import numpy as np
from geom import cu, maps, m, rho_c
alpha_u = rho_c**2
def at_normal_angle(theta, coef=cu):
    # solve phi(s)=theta by Newton, phi(s)= s + arg psi'(e^{is})
    s = theta.copy()
    for it in range(30):
        w, psi, dpsi, ddpsi = maps(coef, s)
        phi = s + np.angle(dpsi)
        # unwrap relative to theta
        d = (phi - theta + np.pi) % (2*np.pi) - np.pi
        cf = (1 + w*ddpsi/dpsi).real
        s = s - d/cf
        if np.max(np.abs(d)) < 1e-15: break
    w, psi, dpsi, ddpsi = maps(coef, s)
    cf = (1 + w*ddpsi/dpsi).real
    z = psi*np.exp(-1j*theta)
    h, hp = z.real, z.imag
    rho = np.abs(dpsi)/cf
    return h, hp, rho, s
if __name__ == '__main__':
    th = np.linspace(0, 2*np.pi/m, 4001)
    h0,hp0,r0,_ = at_normal_angle(th)
    h1,hp1,r1,_ = at_normal_angle(th+np.pi)
    print('max |sqrt(rho)-sqrt(rho_pi)| (scale s->0):', np.max(np.abs(np.sqrt(r0)-np.sqrt(r1))))
    print('rho range', r0.min(), r0.max(), ' h range', h0.min(), h0.max())
    print('max|h\'|', np.abs(hp0).max())
    print('h(theta)-h(theta+pi) max', np.max(np.abs(h0-h1)))
    for s in [0,1e-3,0.1,1,10,30,60,84.13]:
        B = 0
        for t in (-1,1):
            B = max(B, np.max(np.abs(np.sqrt(r0)*np.exp(t*s*hp0) - np.sqrt(r1)*np.exp(-t*s*hp1))))
        T = np.sqrt(2*np.pi)*np.sqrt(alpha_u+s**2)*B
        print('s=%8.3f  B_u(s)=%.6e   T(s)=sqrt(2pi) sqrt(alpha+s^2) B = %.4f' % (s,B,T))
    # width
    wd = h0+h1
    print('width min/max', wd.min(), wd.max(), ' min/max', wd.min()/wd.max())
