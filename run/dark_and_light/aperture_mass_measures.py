import numpy as np, pylab as pl, pyfits as py
import glob
import pymc
import myEmcee_blobs as myEmcee
import lenslib
from astLib import astCalc
from SampleOpt import AMAOpt
import cPickle

# contour plot the inference for J0901 from my dark+light models
name='J0901'
dir = '/data/ljo31b/EELs/galsub/emceeruns/'
file = glob.glob(dir+name+'_parametric_DPL*')
file.sort()
f = file[-1]
print f
result_DPL = np.load(f)

## construct simple model
fp = np.load('/data/ljo31b/EELs/esi/kinematics/inference/results_0.30_source_indous_vdfit_jul2016_J2228.npy')
l,m,u = fp
d = np.mean((l,u),axis=0)
dvl,dvs,dsigmal,dsigmas = d.T
vl,vs,sigmal,sigmas = m.T
s1 = sigmal[1]
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

## MODEL
# compute sigma within 1 re
ML = pymc.Uniform('ML',0.01,80,1.)
gamma = pymc.Uniform('gamma',0.3,2.9,1) 
rho0 = pymc.Uniform('rho0',1,10)
rs = pymc.Uniform('rs',10,600,300.)

pars = [ML,gamma,rho0,rs]
cov = [1.,0.2,1.,20.]

sigma_dm = np.load('/data/ljo31b/EELs/phys_models/models/interpolators/gNFW_aperture_mass_measure_J0901.npy')[()]
Sig_ein_dm = np.load('/data/ljo31b/EELs/phys_models/models/interpolators/gNFW_aperture_mass_measure_Sig_einstein_J0901.npy')[()]

sigma_star = np.loadtxt('/data/ljo31b/EELs/phys_models/models/sigma_star_J0901.dat')[()]
sigma_star_ein = np.loadtxt('/data/ljo31b/EELs/phys_models/models/sigma_star_sigma_einstein_J0901.dat')[()]

@pymc.deterministic
def logP(value=0.,p=pars):
    RHO0 = 10**rho0.value
    ml = ML.value*1#e10
    s2 = ml * sigma_star + RHO0 * sigma_dm.eval(np.column_stack((rs.value,gamma.value)))
    lp1 = -0.5*(SIGMA-s2**0.5)**2. /DSIGMA**2.
    
    Mmod = ml * sigma_star_ein + RHO0 * Sig_ein_dm.eval(np.column_stack((rs.value,gamma.value)))
    Mmod /= (1e10) # same units as data
    lp2 = -0.5*(MEIN-Mmod)**2. / DMEIN**2.
    #print lp1,lp2,Mmod, MEIN, SIGMA,s2**0.5,RHO0,ML.value,rs.value,gamma.value
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

outFile = '/data/ljo31b/EELs/aperture_mass_measure_inference_2'
S = myEmcee.PTEmcee(pars+[likelihood],cov=np.array(cov)/3.,nthreads=8,nwalkers=40,ntemps=10)#28)
S.sample(10000)
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

RHO0 = 10**rho0.value
ml = ML.value*1e10
s2 = ml * sigma_star + RHO0 * sigma_dm.eval(np.column_stack((rs.value,gamma.value)))
lp1 = (SIGMA-s2**0.5)**2. /DSIGMA**2.
    
Mmod = ml * sigma_star_ein + RHO0 * Sig_ein_dm.eval(np.column_stack((rs.value,gamma.value)))
Mmod /= (1e10) # same units as data
lp2 = (MEIN-Mmod)**2. / DMEIN**2.

print MEIN, Mmod

print SIGMA, s2**0.5
pl.show()

for key in dic.keys():
    pl.figure()
    pl.plot(dic[key][:,0])
    pl.title(key)
    pl.show()
