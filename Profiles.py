import numpy as np, pylab as pl
from scipy.special import gamma
from scipy.interpolate import splrep, splint

def PowerLaw(r,sig_crit, rein, eta,projected=False):
    print 'eta',eta
    g1 = gamma(0.5*(1.+eta))
    g2 = gamma(0.5*eta)
    rho_0 = (2.-eta)/(2.*np.pi**0.5) * sig_crit * rein**eta * g1 / g2
    DM_rho = rho_0 * r**-(eta+1.)
    DMmodel = splrep(r,DM_rho*4.*np.pi*r**2.)
    DM = np.array([splint(r[0],r[ii],DMmodel) for ii in range(r.size)])
    if projected:
        sigma = (2.-eta)/2. * sig_crit * (r/rein)**(-eta)
        sigmod = splrep(r,sigma*2.*np.pi*r)
        cumsig = sigma*0.
        for i in range(r.size):
            cumsig[i] = splint(0,r[i],sigmod)
        return cumsig, sigma
    return DM, DM_rho # mass and density profiles


def gNFW(r,sig_crit,rein,gamma,rs,projected=False):
    # normalise density profile
    
    z = np.logspace(-5,5,3500)
    sig = r*0.
    for i in range(r.size):
        R = (r[i]**2. + z**2.)**0.5
        integrand = (R/rs)**(-gamma) * (1. + (R/rs)**2.)**(0.5*gamma - 1.5)
        model = splrep(z,integrand)
        sig[i] = 2. * splint(z[0],z[-1],model)
    # integrate sigma out to r_ein
    model = splrep(r,sig * 2. * np.pi * r)
    sig_rein = splint(0,rein,model)
    rho_0 = sig_crit * np.pi * rein**2. / sig_rein
                   
    if projected:
        # normalise sigma
        sig = sig * sig_crit * np.pi * rein**2. / sig_rein
        sigmod = splrep(r,sig*2.*np.pi*r)
        cumsig = sig*0.
        for i in range(r.size):
            cumsig[i] = splint(0,r[i],sigmod)
        return cumsig, sig # surface mass density
    rho = rho_0 * (r/rs)**(-gamma) * (1. + (r/rs)**2.)**(0.5*gamma - 1.5)
    integrand = rho * 4. * np.pi * r**2.
    model = splrep(r,integrand)
    M = [splint(r[0],r[i],model) for i in range(r.size)]
    return M, rho # mass and density profiles

 
def dlogrho_dlogr_gNFW(r,gamma,rs):
    dpdr = (gamma-3.) * (r**2.)/(r**2 + rs**2.) - gamma
    return dpdr

def dlogrho_dlogr_PL(r,eta):
    return np.ones(r.size)*(-eta - 1.)
