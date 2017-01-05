from linslens.Profiles import *
from astLib import astCalc
import glob, numpy as np, pylab as pl
import lenslib
from jeans.makemodel import deproject
import GetStellarMass
import cPickle
from scipy.interpolate import splev,splrep,splint

names = ['J0837','J0901','J0913','J1125','J1144','J1218','J1323','J1347','J1446','J1605','J1606','J1619','J2228']
Mvir, dMvir, r_eff, dr_eff = np.load('/data/ljo31/Lens/LensParams/Mvir.npy')
lz = np.load('/data/ljo31/Lens/LensParams/LensRedshiftsUpdated.npy')[()]
sz = np.load('/data/ljo31/Lens/LensParams/SourceRedshiftsUpdated.npy')[()]

name = 'J0901'

dir = '/data/ljo31b/EELs/galsub/emceeruns/'
file = glob.glob(dir+name+'_parametric_*')
file.sort()
f = file[-1]
while 1:
    if 'DPL' in f:
        file = file[:-1]
        f = file[-1]
    else:
        break
print f
result_PL = np.load(f)

file = glob.glob(dir+name+'_parametric_DPL*')
file.sort()
f = file[-1]
print f
result_DPL = np.load(f)

lp_PL,trace_PL,dic_PL,_ = result_PL
lp_DPL,trace_DPL,dic_DPL,_ = result_DPL
a1_PL,a3_PL = np.unravel_index(lp_PL[:,0].argmax(),lp_PL[:,0].shape)
a1_DPL,a3_DPL = np.unravel_index(lp_DPL[:,0].argmax(),lp_DPL[:,0].shape)

## HACK
if name == 'J1446':
    nn = 0
    while nn<10:
        row = lp_PL[-1,0,:]
        lp_PL = lp_PL[:,:,row>row.min()]
        trace_PL = trace_PL[:,:,row>row.min(),:]
        for key in dic_PL.keys():
            dic_PL[key] = dic_PL[key][:,:,row>row.min()]
        nn+=1
    nn = 0
    while nn<15:
        row = lp_PL[0,0,:]
        lp_PL = lp_PL[:,:,row>row.min()]
        trace_PL = trace_PL[:,:,row>row.min(),:]
        for key in dic_PL.keys():
            dic_PL[key] = dic_PL[key][:,:,row>row.min()]
        nn+=1
    nn = 0
    while nn<10:
        row = lp_PL[250,0,:]
        lp_PL = lp_PL[:,:,row>row.min()]
        trace_PL = trace_PL[:,:,row>row.min(),:]
        for key in dic_PL.keys():
            dic_PL[key] = dic_PL[key][:,:,row>row.min()]
        nn+=1
    

    pl.figure()
    pl.plot(lp_PL[:,0])
    pl.show()
                                

## construct power-law sigma and mass
zl,zs = lz[name][0], sz[name][0]
sig_crit = lenslib.sig_crit(zl,zs) # solar masses per Mpc^2
sig_crit /= (1e3)**2. # solar masses per kpc^2
scale = astCalc.da(zl)*np.pi/180./3600 * 1e3
rein = dic_PL['Lens 1 b'][a1_PL,0,a3_PL]*0.05*scale # rein in kpc
eta = dic_PL['Lens 1 eta'][a1_PL,0,a3_PL]
eta_PL = dic_PL['Lens 1 eta'][a1_PL,0,a3_PL]

## get rho_0 from lensing
g1 = gamma(0.5*(1.+eta))
g2 = gamma(0.5*eta)
rho_0_PL = (2.-eta)/(2.*np.pi**0.5) * sig_crit * rein**eta * g1 / g2

r = np.logspace(-7,4,4500)
#lr = r[:-100]

### gNFW sigma and mass
rs, rein = dic_DPL['Lens 1 rs'][a1_DPL,0,a3_DPL]*0.05*scale,dic_DPL['Lens 1 b'][a1_DPL,0,a3_DPL]*0.05*scale 
gamma = dic_DPL['Lens 1 gamma'][a1_DPL,0,a3_DPL]

# get rho_0 from lensing
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
rho_0_DPL = sig_crit * np.pi * rein**2. / sig_rein
print '###'

### also old result -- Einstein radius
cat = np.load('/data/ljo31/Lens/LensParams/Structure_lensgals_2src.npy')[0]
rein = cat[name]['Lens 1 b']*0.05*scale

# add in LM!
dir = '/data/ljo31/Lens/LensModels/twoband/'
try:
    oresult = np.load(dir+name+'_212')
except:
    if name == 'J1347':
        oresult = np.load(dir+name+'_112')
    else:
        oresult = np.load(dir+name+'_211')

Mstar_PL = dic_PL['stellar mass'][a1_PL,0,a3_PL]*100
Mstar_DPL = dic_DPL['stellar mass'][a1_DPL,0,a3_DPL]*100

## finally, add original inference for total mass
rein, eta = cat[name]['Lens 1 b']*0.05*scale, cat[name]['Lens 1 eta']

## get kinematics
sigma_dm_DPL = np.load('/data/ljo31b/EELs/phys_models/models/interpolators/gNFW_aperture_mass_measure_'+name+'.npy')[()]
sigma_star = np.loadtxt('/data/ljo31b/EELs/phys_models/models/sigma_star_'+name+'.dat')[()]
sigma_dm_PL = np.load('/data/ljo31b/EELs/phys_models/models/interpolators/PL_aperture_mass_measure_'+name+'.npy')[()]

# how to normalise them
# stellar mass normalised to 10^10 solar masses
sigma_PL = (Mstar_PL*sigma_star**2. + rho_0_PL*sigma_dm_PL.eval(np.array([eta_PL]))**2.)**0.5
sigma_DPL = (Mstar_DPL*sigma_star**2. + rho_0_DPL*sigma_dm_DPL.eval(np.column_stack((rs,gamma)))**2.)**0.5

sigma_s_PL, sigma_d_PL = Mstar_PL**0.5 *sigma_star, rho_0_PL**0.5 *sigma_dm_PL.eval(np.array([eta_PL]))
sigma_s_DPL, sigma_d_DPL = Mstar_DPL**0.5 *sigma_star, rho_0_DPL**0.5 *sigma_dm_DPL.eval(np.column_stack((rs,gamma)))

print sigma_s_PL, sigma_d_PL, sigma_s_DPL, sigma_d_DPL

print sigma_PL, sigma_DPL
