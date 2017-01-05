import numpy as np
import pylab as pl
import pyfits as py
from stellarpop import tools, distances
import itertools
from scipy.interpolate import splrep, splev, splint

dist = distances.Distance()
chabrier = np.load('/home/mauger/python/stellarpop/chabrier.dat')
ages = chabrier['age']
ages = np.log10(ages)
ages[0]=5.
wave=chabrier['wave']
solarSpectra = chabrier[6]

### need to calculate new synthetic photometry using Matt's SEDs! Nobody knows the provenance of these ones....
def MassK(K,z,age=4.):
    logT = np.log10(age)
    kfilt = tools.filterfromfile('Kp_NIRC2')
    modobs = np.zeros(ages.size)
    modrest = modobs*0.
    for a in range(ages.size):
        modobs[a] = tools.ABFM(kfilt,[wave,solarSpectra[a]],z)
        modrest[a] = tools.ABFM(kfilt,[wave,solarSpectra[a]],0.)
    interp_Krest,interp_Kobs  = splrep(ages,modrest),splrep(ages,modobs)
    K_modrest,K_modobs = splev(logT,interp_Krest),splev(logT,interp_Kobs)
    K_rest = K + K_modrest - K_modobs
    K_sun = 5.19
    DL = dist.luminosity_distance(z)
    DM = 5.*np.log10(DL*1e6) - 5.
    K_restabs = K_rest - DM
    lum = -0.4*(K_restabs - K_sun)
    return lum, K_rest

def MassI(I,z,age=4.):
    logT = np.log10(age)
    ifilt = tools.filterfromfile('F814W_ACS')  
    modobs = np.zeros(ages.size)
    modrest = modobs*0.
    for a in range(ages.size):
        modobs[a] = tools.ABFM(ifilt,[wave,solarSpectra[a]],z)
        modrest[a] = tools.ABFM(ifilt,[wave,solarSpectra[a]],0.)
    interp_Irest,interp_Iobs  = splrep(ages,modrest),splrep(ages,modobs)
    I_modrest,I_modobs = splev(logT,interp_Irest),splev(logT,interp_Iobs)
    I_rest = I + I_modrest - I_modobs
    I_sun = 4.57
    DL = dist.luminosity_distance(z)
    DM = 5.*np.log10(DL*1e6) - 5.
    I_restabs = I_rest - DM
    lum = -0.4*(I_restabs - I_sun)
    return lum, I_rest

def MassV1(V,z,age=4.):
    logT = np.log10(age)
    vfilt = tools.filterfromfile('F606W_ACS')  
    modobs = np.zeros(ages.size)
    modrest = modobs*0.
    for a in range(ages.size):
        modobs[a] = tools.ABFM(vfilt,[wave,solarSpectra[a]],z)
        modrest[a] = tools.ABFM(vfilt,[wave,solarSpectra[a]],0.)
    interp_Vrest,interp_Vobs  = splrep(ages,modrest),splrep(ages,modobs)
    V_modrest,V_modobs = splev(logT,interp_Vrest),splev(logT,interp_Vobs)
    V_rest = V + V_modrest - V_modobs
    V_sun = 4.74
    DL = dist.luminosity_distance(z)
    DM = 5.*np.log10(DL*1e6) - 5.
    V_restabs = V_rest - DM
    lum = -0.4*(V_restabs - V_sun)
    return lum, V_rest


def MassV2(V,z,age=4.):
    logT = np.log10(age)
    vfilt = tools.filterfromfile('F555W_ACS')  
    modobs = np.zeros(ages.size)
    modrest = modobs*0.
    for a in range(ages.size):
        modobs[a] = tools.ABFM(vfilt,[wave,solarSpectra[a]],z)
        modrest[a] = tools.ABFM(vfilt,[wave,solarSpectra[a]],0.)
    interp_Vrest,interp_Vobs  = splrep(ages,modrest),splrep(ages,modobs)
    V_modrest,V_modobs = splev(logT,interp_Vrest),splev(logT,interp_Vobs)
    V_rest = V + V_modrest - V_modobs
    V_sun = 4.83
    DL = dist.luminosity_distance(z)
    DM = 5.*np.log10(DL*1e6) - 5.
    V_restabs = V_rest - DM
    lum = -0.4*(V_restabs - V_sun)
    return lum, V_rest

def MassV1toB(V,z,age=4.):
    logT = np.log10(age)
    vfilt = tools.filterfromfile('F606W_ACS')  
    restfilt = tools.filterfromfile('F435W_ACS')
    modobs = np.zeros(ages.size)
    modrest = modobs*0.
    for a in range(ages.size):
        modobs[a] = tools.ABFM(vfilt,[wave,solarSpectra[a]],z)
        modrest[a] = tools.ABFM(restfilt,[wave,solarSpectra[a]],0.)
    interp_Vrest,interp_Vobs  = splrep(ages,modrest),splrep(ages,modobs)
    V_modrest,V_modobs = splev(logT,interp_Vrest),splev(logT,interp_Vobs)
    V_rest = V + V_modrest - V_modobs
    V_sun = 4.74
    DL = dist.luminosity_distance(z)
    DM = 5.*np.log10(DL*1e6) - 5.
    V_restabs = V_rest - DM
    lum = -0.4*(V_restabs - V_sun)
    return lum, V_rest

def MassV2toB(V,z,age=4.):
    logT = np.log10(age)
    vfilt = tools.filterfromfile('F555W_ACS')  
    restfilt = tools.filterfromfile('F435W_ACS')
    modobs = np.zeros(ages.size)
    modrest = modobs*0.
    for a in range(ages.size):
        modobs[a] = tools.ABFM(vfilt,[wave,solarSpectra[a]],z)
        modrest[a] = tools.ABFM(restfilt,[wave,solarSpectra[a]],0.)
    interp_Vrest,interp_Vobs  = splrep(ages,modrest),splrep(ages,modobs)
    V_modrest,V_modobs = splev(logT,interp_Vrest),splev(logT,interp_Vobs)
    V_rest = V + V_modrest - V_modobs
    V_sun = 4.83
    DL = dist.luminosity_distance(z)
    DM = 5.*np.log10(DL*1e6) - 5.
    V_restabs = V_rest - DM
    lum = -0.4*(V_restabs - V_sun)
    return lum, V_rest
