import numpy as np, pylab as pl
from scipy.special import gamma
from scipy.interpolate import splrep, splint,spalde,splev

r = np.logspace(-6,5,2601)
lr=r[:-100]
rein = 1.

def anal(eta):
    num = np.pi**0.5 * lr**(eta-3.) * gamma(1.5-0.5*eta)
    denom = 2.*gamma(2.-eta/2.)
    rho = num/denom * (2.-eta)/(2.*np.pi*rein**(eta-2))
    return rho

def num(eta):
    kappa = 0.5 * (r/rein)**(eta-2.)
    kappamodel = splrep(r,kappa)
    d_kappa = np.array(spalde(r,kappamodel))[:,1]
    dmodel = splrep(r,d_kappa)
    rho = lr*0.
    for i in range(lr.size):
        R = lr[i]
        rr = np.logspace(-5,0.5*np.log10(r[-1]**2-R**2),1001)
        y = (rr**2+R**2)**0.5
        f = splev(y,dmodel)
        model = splrep(rr,-1*f/(rr**2+R**2)**0.5/np.pi)
        rho[i] = splint(rr[0],rr[-1],model)        
    return rho

ETA = 1.5
rho1 = anal(ETA)
rho2 = num(ETA)

pl.plot(lr[100:],rho1[100:])
pl.plot(lr[100:],rho2[100:])
pl.show()

pl.plot(lr[100:],(rho1[100:]-rho2[100:])/rho2[100:])
pl.show()
