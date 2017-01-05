from tools import solarmag
import numpy as np, pylab as pl, pyfits as py
import cPickle
from jeans.makemodel import *
from astLib import astCalc
from imageSim import SBObjects
from itertools import product
import ndinterp
from multiprocessing import Pool
from linslens import EELsModels_huge as L, EELsKeckModels as K
import lenslib
from scipy.special import beta, hyp2f1 as hyp

def run(gamma,rs,b):

    mass = PROFILES.gNFW_TC(x=0.,y=0.,eta=gamma,rs=rs,pa=0.,q=1.,b=b,zl=zl,zs=zs) 
    sigma = mass.sigma(lr)
    c = np.where(np.isfinite(sigma)==True)
    sigmod = splrep(lr[c],sigma[c]*2.*np.pi*lr[c])
    M_ein = splint(0,Rein_tot,sigmod)

    return M_ein


phot = py.open('/data/ljo31/Lens/LensParams/Phot_1src_huge_new_new.fits')[1].data
names = phot['name']
lz = np.load('/data/ljo31/Lens/LensParams/LensRedshiftsUpdated.npy')[()]
sz = np.load('/data/ljo31/Lens/LensParams/SourceRedshiftsUpdated.npy')[()]
cat = np.load('/data/ljo31/Lens/LensParams/Structure_lensgals_2src.npy')[0]

name = 'J0901'
zl,zs = lz[name][0],sz[name][0]
scale = astCalc.da(zl)*np.pi/180./3600 * 1e3
Rein_tot = cat[name]['Lens 1 b']*0.05*scale

r = np.logspace(-5,5,2501)
lr = r[:-100]

sigma_dm = np.load('/data/ljo31b/EELs/phys_models/models/interpolators/gNFW_aperture_mass_measure_Sig_einstein_'+name+'.npy')[()]

GAMMAS = np.random.uniform(0.01,2.9,100)
RSS = np.random.uniform(2,550.,100)
BS = np.random.uniform(1.,5,100)

s1,s2 = RSS*0.,BS*0.
for ii in range(len(GAMMAS)):
    gg,rr,bb = GAMMAS[ii],RSS[ii],BS[ii]
    bb=1.
    vd1 = run(gg,rr,bb)
    
    arr = np.column_stack((rr,gg))
   
    B,G = bb/rr,gg
    R0 = rr
    N=3.0001
    drb = 2*R0/B*(beta(((N-3.)/2.),((3.-G)/2.)) - beta(((N-3.)/2.),(3./2.)) * (1+B**2)**((3-N)/2.) * hyp((N-3.)/2.,G/2.,N/2., 1./(1+B**2))) 
    B = 1./rr
    dr1 = 2*R0/B*(beta(((N-3.)/2.),((3.-G)/2.)) - beta(((N-3.)/2.),(3./2.)) * (1+B**2)**((3-N)/2.) * hyp((N-3.)/2.,G/2.,N/2., 1./(1+B**2))) 
    kappa = bb * dr1/drb

    vd2 = kappa * sigma_dm.eval(arr)

    print '%f,%f'%(vd1, vd2)
    s1[ii],s2[ii] = vd1, vd2

pl.figure()
pl.scatter(s1,s2,c='SteelBlue',s=40)
pl.show()
