import numpy as np, pylab as pl
from scipy.special import gamma, hyp2f1 as hyp, beta
from scipy.interpolate import splrep, splint,spalde,splev
import lenslib

class MassProfile:
    def align_coords(self,xin,yin,offset=True,revert=False):
        from math import cos,sin,pi
        if offset:
            theta = self.theta-pi/2.
        else:
            theta = self.theta
        ctheta = cos(theta)
        stheta = sin(theta)
        if revert:
            X = xin*ctheta-yin*stheta
            Y = yin*ctheta+xin*stheta
            x = X+self.x
            y = Y+self.y
            return x,y
        X = xin-self.x
        Y = yin-self.y
        x = X*ctheta+Y*stheta
        y = Y*ctheta-X*stheta
        return x,y

class PowerLaw(MassProfile):
    def __init__(self,b=None,eta=1.,pa=None,q=None,x=None,y=None,zl=None,zs=None):
        self.b = b
        self.eta = eta
        self.pa = pa
        self.q = q
        self.x = x
        self.y = y
        self.zl = zl
        self.zs = zs

    def kappa(self,R):
        return 0.5 * (R/self.b)**(self.eta-2.)

    def sigma_crit(self):
        self.sig_crit = lenslib.sig_crit(self.zl,self.zs) # solar masses per Mpc^2
        self.sig_crit /= (1e3)**2. # solar masses per kpc^2
        return self.sig_crit

    def sigma(self,R):
        sig_crit = self.sigma_crit()
        kappa = self.kappa(R)
        return sig_crit * kappa

    def rho_num(self,r):
        sigma = self.sigma(r)
        sigmamodel = splrep(r,sigma)
        d_sigma = np.array(spalde(r,sigmamodel))[:,1]
        dmodel = splrep(r,d_sigma)
        lr = r[:-100]
        rho = lr*0.
        for i in range(lr.size):
            R = lr[i]
            rr = np.logspace(-5,0.5*np.log10(r[-1]**2-R**2),1001)
            y = (rr**2+R**2)**0.5
            f = splev(y,dmodel)
            model = splrep(rr,-1*f/(rr**2+R**2)**0.5/np.pi)
            rho[i] = splint(rr[0],rr[-1],model)        
        return lr, rho

    def rho(self,r):
        sig_crit = self.sigma_crit()
        num = (2.-self.eta) * self.b**(2.-self.eta) * r**(self.eta-3.) * gamma(1.5 - self.eta/2.)
        denom = 4. * gamma(2.-self.eta/2.) * np.pi**0.5
        return sig_crit * num/denom


class gNFW(MassProfile):
    ''' currently must be spherical -- normalise like MWA'''
    def __init__(self,b=None,eta=1.,rs=None,pa=None,q=None,x=None,y=None,zl=None,zs=None,):
        self.b = b
        self.eta = eta
        self.rs = rs
        self.pa = pa
        self.q = q
        self.x = x
        self.y = y
        self.zl = zl
        self.zs = zs
        self.ks = None
        self.rho_s = None

    def sigma_crit(self):
        self.sig_crit = lenslib.sig_crit(self.zl,self.zs) # solar masses per Mpc^2
        self.sig_crit /= (1e3)**2. # solar masses per kpc^2
        return self.sig_crit

    def getKs(self):
        B = self.b/self.rs
        rr = np.logspace(-6,0,61)*B
        rr2p1 = 1.+rr**2
        kappa = 2.*hyp(1.,self.eta/2.,1.5,1./rr2p1)/rr2p1
        c = np.isfinite(kappa)
        kappa = kappa[c]
        rr = rr[c]
        model = splrep(rr,kappa*rr)
        ks = B**2/(2*splint(rr[0],B,model))
        self.ks = ks
        return ks

    def rho(self,r):
        rr = r/self.rs
        rho = rr**self.eta * (1. + rr**2.)**(0.5*(3.-self.eta))
        ks = self.getKs()
        sig_crit = self.sigma_crit()
        self.rho_s = ks * sig_crit / self.rs
        rho = self.rho_s / rho
        return rho

    def kappa(self,R):
        rr = R/self.rs
        rr2p1 = 1.+rr**2
        sig_crit = self.sigma_crit()
        kappa = 2.*hyp(1.,self.eta/2., 1.5, 1./rr2p1)/rr2p1
        ks = self.getKs()
        kappa *= ks
        return kappa

    def sigma(self,R):
        rr = R/self.rs
        kappa = self.kappa(R)
        sig_crit = self.sigma_crit()
        return kappa * sig_crit
    

#PL = PowerLaw(x=0,y=0,b=10.,pa=0,q=1.,zl=0.4,zs=0.7,eta=1.)

PL = PowerLaw(x=0,y=0,b=10.,pa=0,q=1.,zl=0.4,zs=0.7,eta=1.)
DPL = gNFW(x=0,y=0,b=10.,pa=0,q=1.,zl=0.4,zs=0.7,eta=2.,rs=1e6)


r = np.logspace(-5,5,2501)
kappa = PL.kappa(r)
rho = PL.rho(r)

rho2 = DPL.rho(r)
kappa2 = DPL.kappa(r)

# test for eta = 1
# RHO = 0.5 * r**-2 * 10. * PL.sig_crit / np.pi
print rho2/rho
print kappa2/kappa

'''pl.figure()
pl.loglog(r,kappa)
pl.loglog(r,kappa2)
pl.figure()
pl.loglog(r,rho)
pl.loglog(r,rho2)
pl.show()'''

pl.figure()
pl.plot(r,kappa2/kappa)
pl.figure()
pl.plot(rho2/rho)
pl.show()
