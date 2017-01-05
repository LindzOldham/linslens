import numpy as np, pylab as pl
from scipy.special import gamma,gammainc, hyp2f1 as hyp, beta
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
        return (1.-0.5*self.eta) * (R/self.b)**-self.eta * self.q**(-self.eta/2.) # right?

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

    '''def rho(self,r):
        sig_crit = self.sigma_crit()
        num = self.eta * self.b**-self.eta * gamma(0.5 + 0.5*self.eta)
        denom = 4. * np.pi**0.5 * gamma(1.+0.5*self.eta)
        return sig_crit * num/denom * r**(-1.-self.eta)'''

    def rho(self,r):
        sig_crit = self.sigma_crit()
        num = (1.-self.eta/2.) * self.q**(-0.5*self.eta) * self.b**self.eta * gamma((self.eta+1.)/2.)
        denom = np.pi**0.5 * gamma(self.eta/2.)
        rho_0 = num/denom
        return sig_crit * rho_0 * r**(-1.-self.eta)

    def mass(self,r):
        sig_crit = self.sigma_crit()
        num = (1.-self.eta/2.) * self.q**(-0.5*self.eta) * self.b**self.eta * gamma((self.eta+1.)/2.)
        denom = np.pi**0.5 * gamma(self.eta/2.)
        rho_0 = num/denom
        M_0 = rho_0 * 4.*np.pi/(2.-self.eta)
        return sig_crit * M_0 * r**(2.-self.eta)

    '''def mass(self,r):
        sig_crit = self.sigma_crit()
        num = self.eta * self.b**(-1.-self.eta) * gamma(0.5 + 0.5*self.eta)
        denom = 4. * np.pi**0.5 * gamma(1.+0.5*self.eta)
        rho_0 = num/denom
        M_0 = rho_0 * 4.*np.pi/(2.-self.eta)
        return sig_crit * M_0 * r**(2.-self.eta)'''

    def mass_proj(self,R):
        sigma = self.sigma(R)
        lr = R[:-100]
        model = splrep(R,sigma*2.*np.pi*R)
        M = [splint(0,lr[i],model) for i in range(lr.size)]
        return M


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

    def kappa_b_1(self,R):
        rr = R/self.rs
        rr2p1 = 1.+rr**2
        sig_crit = self.sigma_crit()
        kappa = 2.*hyp(1.,self.eta/2., 1.5, 1./rr2p1)/rr2p1
        return kappa

    def sigma(self,R):
        rr = R/self.rs
        kappa = self.kappa(R)
        sig_crit = self.sigma_crit()
        return kappa * sig_crit
    

class Sersic(MassProfile):
    def __init__(self,b=None,n=None,pa=None,q=None,x=None,y=None,re=None,zl=None,zs=None):
        self.b = b
        self.n = n
        self.pa = pa
        self.q = q
        self.x = x
        self.y = y
        self.re = re
        self.zl = zl
        self.zs = zs

    def sigma_crit(self):
        self.sig_crit = lenslib.sig_crit(self.zl,self.zs) # solar masses per Mpc^2
        self.sig_crit /= (1e3)**2. # solar masses per kpc^2
        return self.sig_crit

    def getbFromMass(self,mass):
        sig_crit = self.sigma_crit()
        n,re = self.n,self.re
        b = np.logspace(-3,2,501)*re
        k = 2.*n-1./3+4./(405.*n)+46/(25515.*n**2)
        amp = b*b*k**(2*n)/(2*n*re**2*(gamma(2*n)*gammainc(2*n,k*(b/re)**(1/n))))
        m = 2*np.pi*sig_crit*re**2*amp*gamma(2*n)*n/k**(2*n)
        model = splrep(m,b)
        return splev(mass,model)

    def setbFromMass(self,mass):
        self.b = self.getbFromMass(mass)

    def kappa(self,R):
        b = self.b
        n = self.n
        re = self.re
        k = 2.*n-1./3+4./(405.*n)+46/(25515.*n**2)

        amp = b*b*k**(2*n)/(np.exp(k)*2*n*re**2*(gamma(2*n)*gammainc(2*n,k*(b/re)**(1/n))))
        return amp*np.exp(-k*((R/re)**(1./n)-1.))

    def sigma(self,R):
        sig_crit = self.sigma_crit()
        kappa = self.kappa(R)
        return sig_crit * kappa

    def rho(self,r):
        # deproject sigma
        from jeans.makemodel import deproject
        sigma = self.sigma(r)
        lr,rho = deproject(r,sigma)
        return lr,rho

    def mass(self,r):
        from jeans.makemodel import light2mass
        lr,rho = self.rho(r)
        M = light2mass(lr,rho,1.)
        return M

class gNFW_TC(MassProfile):
    ''' normalise like Tom -- this is actually the same as Matt (except the factor of two) and Keeton '''
    def __init__(self,b=None,eta=1.,rs=None,pa=None,q=None,x=None,y=None,zl=None,zs=None):
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
        G = self.eta
        N = 3.00001 
        B = self.b/self.rs
        R0 = self.rs
        drb = 2*R0/B*(beta(((N-3.)/2.),((3.-G)/2.)) 
            - 
            beta(((N-3.)/2.),(3./2.)) 
            * 
            (1+B**2)**((3-N)/2.) 
            * 
            hyp( 
                (N-3.)/2., 
                G/2., 
                N/2., 
                1./(1+B**2) 
                ) 
            ) 
        self.ks=self.b/drb
        ks=self.ks
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
        #sig_crit = self.sigma_crit()
        kappa = 2.*hyp(1.,self.eta/2., 1.5, 1./rr2p1)/rr2p1
        ks = self.getKs()
        kappa *= ks
        return kappa

    def sigma(self,R):
        kappa = self.kappa(R)
        sig_crit = self.sigma_crit()
        return kappa * sig_crit

    def mass(self,r):
        rho = self.rho(r)
        lr = r[:-100]
        model = splrep(r,4.*np.pi*r**2. * rho)
        M = [splint(0,lr[i],model) for i in range(lr.size)]
        return M

    def mass_proj(self,R):
        lr = R[:-100]
        sigma = self.sigma(lr)
        c = np.isfinite(sigma)==True
        model = splrep(lr[c],sigma[c]*2.*np.pi*lr[c])
        M = [splint(0,lr[i],model) for i in range(lr.size)]
        return M

