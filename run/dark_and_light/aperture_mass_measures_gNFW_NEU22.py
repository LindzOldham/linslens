import numpy as np, pylab as pl, pyfits as py
import glob
import pymc
import myEmcee_blobs as myEmcee
import lenslib
from astLib import astCalc
from SampleOpt import AMAOpt
import cPickle
from scipy.special import beta,hyp2f1 as hyp
import sys

name = sys.argv[1] 
print name


## construct simple model
dir3 = '/data/ljo31b/EELs/esi/kinematics/inference/vdfit/NEW/'
result_K = np.load(dir3+name+'_1.00_lens_esi_indous_vdfit_LENS')
lp_K, trace_K, dic_K, _ = result_K
s1 = np.median(dic_K['lens dispersion'])
dsigmas = s1*0.05

# Einstein radius
lz = np.load('/data/ljo31/Lens/LensParams/LensRedshiftsUpdated.npy')[()]
sz = np.load('/data/ljo31/Lens/LensParams/SourceRedshiftsUpdated.npy')[()]
zl,zs = lz[name][0],sz[name][0]
scale = astCalc.da(zl)*np.pi/180./3600 * 1e3

cat = np.load('/data/ljo31/Lens/LensParams/Structure_lensgals_2src.npy')[0]
rein = cat[name]['Lens 1 b']*0.05*scale # in kpc
dr1  = np.load('/data/ljo31/Lens/LensParams/Structure_lensgals_2src.npy')[1][name]['Lens 1 b']
dr2 = np.load('/data/ljo31/Lens/LensParams/Structure_lensgals_2src.npy')[2][name]['Lens 1 b']
drein = np.mean((dr1,dr2))*0.05*scale

sig_crit = lenslib.sig_crit(zl,zs) # solar masses per Mpc^2
sig_crit /= (1e3)**2. # solar masses per kpc^2

## DATAPOINTS
SIGMA,MEIN = s1, sig_crit * np.pi * rein**2. * 1e-10 #units 10^10 Mdot
DSIGMA = SIGMA*0.05
DMEIN = 2.* drein / rein * MEIN #* 10.
print 'DMEIN', DMEIN

## MODEL
# compute sigma within 1 re
ML = pymc.Uniform('ML',0.001,80,0.05) # this is in units of 1e12 solar masses
GAMMA = pymc.Uniform('gamma',0.3,2.9,1) 
B = pymc.Uniform('B',0.05,10)# this is in kpc
RS = pymc.Uniform('RS',0.05,200)#200,50) # also in kpc

pars = [ML,GAMMA,B,RS]
cov = [1.,0.2,1.,5.]

## get kinematics
sigma_star = np.loadtxt('/data/ljo31b/EELs/phys_models/models/sigma_star_'+name+'.dat')[()]
sigma_dm = np.load('/data/ljo31b/EELs/phys_models/models/interpolators/gNFW_aperture_mass_measure_'+name+'.npy')[()]

## and Einstein masses
Sig_ein_dm = np.load('/data/ljo31b/EELs/phys_models/models/interpolators/gNFW_aperture_mass_measure_Sig_einstein_'+name+'.npy')[()]
Sig_ein_star = np.loadtxt('/data/ljo31b/EELs/phys_models/models/sigma_star_sigma'+name+'.dat')[()]

# POWERLAW quantities need to be normalised by b**eta
# gNFW quantities need to be normalised in a more complicated way...maybe I should grid in r/rs?

@pymc.deterministic
def logP(value=0.,p=pars):
    ml,gamma,b,rs = [pars[i].value for i in range(len(pars))]
    B,G = b/rs,gamma
    R0 = rs
    N=3.0001
    drb = 2*R0/B*(beta(((N-3.)/2.),((3.-G)/2.)) - beta(((N-3.)/2.),(3./2.)) * (1+B**2)**((3-N)/2.) * hyp((N-3.)/2.,G/2.,N/2., 1./(1+B**2))) 
    B = 1./rs
    dr1 = 2*R0/B*(beta(((N-3.)/2.),((3.-G)/2.)) - beta(((N-3.)/2.),(3./2.)) * (1+B**2)**((3-N)/2.) * hyp((N-3.)/2.,G/2.,N/2., 1./(1+B**2))) 
    kappa = b * dr1/drb

    s2 = ml * sigma_star**2. + kappa * sigma_dm.eval(np.column_stack((rs,gamma)))
    lp1 = -0.5*(SIGMA-s2**0.5)**2. /DSIGMA**2.
    
    Mmod = ml * Sig_ein_star + kappa * Sig_ein_dm.eval(np.column_stack((rs,gamma)))
    Mmod /= (1e10) 
    lp2 = -0.5*(MEIN-Mmod)**2. / DMEIN**2.

    #print 'VELDISP', ml * sigma_star**2., kappa * sigma_dm.eval(np.column_stack((rs,gamma)))
    #print 'KAPPA', ml * Sig_ein_star, Sig_ein_dm.eval(np.column_stack((rs,gamma)))

    return lp1+lp2

@pymc.observed
def likelihood(value=0.,lp=logP):
    return lp

SS = AMAOpt(pars,[likelihood],[logP],cov=cov)
SS.sample(4000)
lp,trace,det = SS.result()

print 'results from optimisation:'
for i in range(len(pars)):
    pars[i].value = trace[-1,i]
    print "%18s  %8.3f"%(pars[i].__name__,pars[i].value)

outFile = '/data/ljo31b/EELs/aperture_mass_measure_inference_gNFW_'+name
S = myEmcee.PTEmcee(pars+[likelihood],cov=np.array(cov)/3.,nthreads=8,nwalkers=40,ntemps=10)#28)
S.sample(2000)
f = open(outFile,'wb')
cPickle.dump(S.result(),f,2)
f.close()
result = S.result()
result = np.load(outFile)
lp,trace,dic,_ = result
a1,a3 = np.unravel_index(lp[:,0].argmax(),lp[:,0].shape)
for i in range(len(pars)):
    pars[i].value = np.median(trace[200:,0,:,i])
    print "%18s  %8.5f"%(pars[i].__name__,pars[i].value)

pl.figure()
pl.plot(lp[200:,0])

ml,gamma,b,rs = [pars[i].value for i in range(len(pars))]
B,G = b/rs,gamma
R0 = rs
N=3.0001
drb = 2*R0/B*(beta(((N-3.)/2.),((3.-G)/2.)) - beta(((N-3.)/2.),(3./2.)) * (1+B**2)**((3-N)/2.) * hyp((N-3.)/2.,G/2.,N/2., 1./(1+B**2))) 
B = 1./rs
dr1 = 2*R0/B*(beta(((N-3.)/2.),((3.-G)/2.)) - beta(((N-3.)/2.),(3./2.)) * (1+B**2)**((3-N)/2.) * hyp((N-3.)/2.,G/2.,N/2., 1./(1+B**2))) 
kappa = b * dr1/drb

s2 = ml * sigma_star**2. + kappa * sigma_dm.eval(np.column_stack((rs,gamma)))
lp1 = -0.5*(SIGMA-s2**0.5)**2. /DSIGMA**2.

Mmod = ml * Sig_ein_star + kappa * Sig_ein_dm.eval(np.column_stack((rs,gamma)))
Mmod /= (1e10) 
lp2 = -0.5*(MEIN-Mmod)**2. / DMEIN**2.

print 'lp1,lp2', '%.2f,%.2f'%(lp1,lp2)


print MEIN, Mmod

print SIGMA, s2**0.5
pl.show()

for key in dic.keys():
    pl.figure()
    pl.plot(dic[key][:,0])
    pl.title(key)
    pl.show()
