import numpy as np, pylab as pl
from imageSim import SBObjects
from distances import D
from jeans.makemodel import *
import time
import lenslib
from multiprocessing import Pool
from scipy.interpolate import splrep, splev
import pymc
from scipy.special import beta,hyp2f1 as hyp
import myEmcee_blobs as myEmcee
from SampleOpt import AMAOpt

# lens SL2SJ142059+563007

name = 'SL2SJ142059+563007'
Reff = 1.62
q=0.54
pa=-12.6
n=4.
Rein = 1.4
dRein = 0.04
ql = 0.67
dql = 0.01
pa = -10.4
dpa = 0.3
z = 0.483
zs = 3.12

scale = 1000*D.Da(z)/206265.
Reff*=scale
Rein*=scale
sig_crit = lenslib.sig_crit(z,zs)/1e6

############################
######### DYNAMICS #########
############################

# STELLAR MASS VELOCITY DISPERSION
gal = PROFILES.Sersic(x=0.,y=0.,n=4.,re=Reff,q=1.,pa=0.,b=50.,zl=z,zs=zs)
gal.setbFromMass(1e12)
r = np.logspace(-5,5,2501)
lr = r[:-100]
mass = gal.mass(r)
sb = gal.sigma(r)
sigma_star = veldisp(r,sb,mass,ap=Reff/2.)
print sigma_star**0.5

# DM VELOCITY DISPERSION
try:
    sigma_dm = np.load('/data/ljo31b/EELs/phys_models/models/sigma_gNFW_SL2Slens.npy')[()]
except:
    start = time.time()
    gammagrid = np.arange(0,2.9,0.01)
    idx = list(range(len(gammagrid)))
    p = Pool(8)
    lr,light = deproject(r,sb)
    arr = [[r,10.*Reff,gammagrid[m],z,zs] for m in range(len(gammagrid))]
    Mdm = np.zeros((lr.size,gammagrid.size))
    sigma_dm = gammagrid*0.
    out = p.map(gridgNFW,arr)
    for i in range(len(arr)):
        Mdm[:,idx[i]] = out[i]
    print 'got mass profiles'
    arr = [[r,sb,Mdm[:,idx[i]],Reff/2.] for i in range(len(arr))]
    out = p.map(gridveldisp,arr)
    for i in range(len(arr)):
        sigma_dm[idx[i]] = out[i]
    print 'got velocity dispersions'
    print sigma_dm**0.5
    np.save('/data/ljo31b/EELs/phys_models/models/sigma_gNFW_SL2Slens.npy',sigma_dm)
    end = time.time()
    print 'time:', '%.3f'%(end-start)

try:
    obj = np.load('/data/ljo31b/EELs/phys_models/models/sigmagrid_gNFW_SL2Slens.npy')
except:
    start = time.time()
    ## make grid
    gammagrid = np.arange(0,2.9,0.01)
    obj = splrep(gammagrid,sigma_dm**0.5,k=1)
    np.save('/data/ljo31b/EELs/phys_models/models/sigmagrid_gNFW_SL2Slens.npy',obj)
    end = time.time()
    print 'time:', '%.3f'%(end-start)

    pl.plot(gammagrid,(sigma_dm**0.5-splev(gammagrid,obj))/sigma_dm**0.5)
    pl.show()
            

############################
######### LENSING ##########
############################

## STELLAR MASS EINSTEIN MASS
sigmod = splrep(r,sb*2.*np.pi*r)
Mein_star = splint(0,Rein,sigmod)
print Mein_star

## DM EINSTEIN MASS
try:
    Mein_dm = np.load('/data/ljo31b/EELs/phys_models/models/Mein_gNFW_SL2Slens.npy')
except:
    start = time.time()
    gammagrid = np.arange(0,2.9,0.01)
    idx = list(range(len(gammagrid)))
    p = Pool(8)
    
    arr = [[r,10.*Reff,G,z,zs,Rein] for G in gammagrid]
    Mein_dm = gammagrid*0.
    out = p.map(gridgNFW_sigma,arr)
    for i in range(len(arr)):
        Mein_dm[idx[i]] = out[i]
    print 'got Einstein masses'
    np.save('/data/ljo31b/EELs/phys_models/models/Mein_gNFW_SL2Slens.npy',Mein_dm)
        
    end = time.time()
    print 'time:', '%.3f'%(end-start)

try:
    obj_Mein = np.load('/data/ljo31b/EELs/phys_models/models/Meingrid_gNFW_SL2Slens.npy')[()]
except:
    start = time.time()
    gammagrid = np.arange(0,2.9,0.01)
    obj_Mein = splrep(gammagrid,Mein_dm,k=1)
    np.save('/data/ljo31b/EELs/phys_models/models/Meingrid_gNFW_SL2Slens.npy',obj_Mein)

    end = time.time()
    print 'time:', '%.3f'%(end-start)

    pl.plot(gammagrid,(Mein_dm-splev(gammagrid,obj_Mein))/Mein_dm)  
    pl.show()
    

### NOW DO INFERENCE
### LOAD DATA
SIGMA, DSIGMA = 228.,19.
MEIN = sig_crit * np.pi * Rein**2. * 1e-10
DMEIN = 2.* dRein / Rein * MEIN
print MEIN, DMEIN


ML = pymc.Uniform('ML',0.001,80,0.05) # this is in units of 1e12 solar masses
#ML = pymc.Uniform('ML',10,12,11.5) # this is in units of 1e12 solar masses
GAMMA = pymc.Uniform('gamma',0.3,2.9,1) 
B = pymc.Uniform('B',0.05,10)# this is in kpc

pars = [ML,GAMMA,B]
cov = [1.,0.2,1.]

rs=10.*Reff

@pymc.deterministic
def logP(value=0.,p=pars):
    ml,gamma,b = [pars[i].value for i in range(len(pars))]
    #ml = 10**ml / (1e12)
    B,G = b/rs,gamma
    R0 = rs
    N=3.0001
    drb = 2*R0/B*(beta(((N-3.)/2.),((3.-G)/2.)) - beta(((N-3.)/2.),(3./2.)) * (1+B**2)**((3-N)/2.) * hyp((N-3.)/2.,G/2.,N/2., 1./(1+B**2))) 
    B = 1./rs
    dr1 = 2*R0/B*(beta(((N-3.)/2.),((3.-G)/2.)) - beta(((N-3.)/2.),(3./2.)) * (1+B**2)**((3-N)/2.) * hyp((N-3.)/2.,G/2.,N/2., 1./(1+B**2))) 
    kappa = b * dr1/drb

    s2 = ml * sigma_star + kappa * splev(gamma,obj)**2.
    lp1 = -0.5*(SIGMA-s2**0.5)**2. /DSIGMA**2.

    Mmod = ml * Mein_star + kappa * splev(gamma,obj_Mein)
    Mmod /= (1e10)
    lp2 = -0.5*(MEIN-Mmod)**2. / DMEIN**2.

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

outFile = '/data/ljo31b/EELs/'+name+'_inference'
S = myEmcee.PTEmcee(pars+[likelihood],cov=np.array(cov)/3.,nthreads=8,nwalkers=40,ntemps=10)
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

#pl.figure()
#pl.plot(lp[200:,0])

ml,gamma,b = [pars[i].value for i in range(len(pars))]
B,G = b/rs,gamma
R0 = rs
N=3.0001
drb = 2*R0/B*(beta(((N-3.)/2.),((3.-G)/2.)) - beta(((N-3.)/2.),(3./2.)) * (1+B**2)**((3-N)/2.) * hyp((N-3.)/2.,G/2.,N/2., 1./(1+B**2))) 
B = 1./rs
dr1 = 2*R0/B*(beta(((N-3.)/2.),((3.-G)/2.)) - beta(((N-3.)/2.),(3./2.)) * (1+B**2)**((3-N)/2.) * hyp((N-3.)/2.,G/2.,N/2., 1./(1+B**2))) 
kappa = b * dr1/drb

s2 = ml * sigma_star + kappa * splev(gamma,obj)**2.
lp1 = -0.5*(SIGMA-s2**0.5)**2. /DSIGMA**2.

Mmod = ml * Mein_star + kappa * splev(gamma,obj_Mein)
Mmod /= (1e10)
lp2 = -0.5*(MEIN-Mmod)**2. / DMEIN**2.

print MEIN, Mmod

print SIGMA, s2**0.5
'''pl.show()

for key in dic.keys():
    pl.figure()
    pl.plot(dic[key][:,0])
    pl.title(key)
    pl.show()'''

trace = trace[500:,0].reshape((1500*40,3))
trace[:,0]= np.log10(trace[:,0]*1e12)
import corner_plot
corner_plot.corner_plot(trace,axis_labels=['ML','GAMMA','B'])
pl.show()
