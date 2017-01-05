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

def run(kresult,result,name,gamma,rs,b):

    mass = PROFILES.gNFW_TC(x=0.,y=0.,eta=gamma,rs=rs,pa=0.,q=1.,b=b,zl=zl,zs=zs) 
    Mdm = mass.mass(r)
    
    arr = [r,sb,Mdm,scale]
    sigma = gridveldisp(arr)
   
    return sigma


phot = py.open('/data/ljo31/Lens/LensParams/Phot_1src_huge_new_new.fits')[1].data
names = phot['name']
lz = np.load('/data/ljo31/Lens/LensParams/LensRedshiftsUpdated.npy')[()]
sz = np.load('/data/ljo31/Lens/LensParams/SourceRedshiftsUpdated.npy')[()]
cat = np.load('/data/ljo31/Lens/LensParams/Structure_lensgals_2src.npy')[0]


dir = '/data/ljo31/Lens/LensModels/twoband/'

name = 'J1144'
gg,rr,bb = 1.,10.,10.

try:
    result = np.load(dir+name+'_212')
    kresult = np.load(dir+name+'_Kp_212')
except:
    if name == 'J1347':
        result = np.load(dir+name+'_112')
        kresult = np.load('/data/ljo31/Lens/J1347/twoband_Kp_112_2')
    elif name == 'J1619':
        result = np.load(dir+name+'_212')
        kresult = np.load(dir+name+'_Kp_212_lensandgalon')
        print 'J1619'                                                                                                                                                                                                       
    else:
        result = np.load(dir+name+'_211')
        kresult = np.load(dir+name+'_Kp_211')
        print 'here'

model = L.EELs(result,name)
kmodel = K.EELs(kresult,result,name)

model.Initialise()
kmodel.Initialise()

zl,zs = lz[name][0],sz[name][0]
scale = astCalc.da(zl)*np.pi/180./3600 * 1e3
    
sig_crit = lenslib.sig_crit(zl,zs) # solar masses per Mpc^2
sig_crit /= (1e3)**2. # solar masses per kpc^2

Rein = cat[name]['Lens 1 b']*0.05*scale

fits = kmodel.fits
r = np.logspace(-5,5,2501)

# build lens galaxy light profiles
if model.galno == 1.:
    gal1 = model.gals[0]
    gal1.re *= 0.05*scale
    sb = gal1.eval(r)

else:
    gal1,gal2 = model.gals
    gal1.re *= 0.05*scale
    gal2.re *= 0.05*scale
    sb = fits[0]*gal1.eval(r) + fits[1]*gal2.eval(r)

# stellar mass profile
lr,light = deproject(r,sb)

sigma_dm = np.load('/data/ljo31b/EELs/phys_models/models/interpolators/gNFW_aperture_mass_measure_'+name+'.npy')[()]

GAMMAS = np.random.uniform(0.01,2.9,100)
RSS = np.random.uniform(1.5,550.,100)
BS = np.random.uniform(1.,50,100)

s1,s2 = RSS*0.,BS*0.
for ii in range(len(GAMMAS)):
    gg,rr,bb = GAMMAS[ii],RSS[ii],BS[ii]
    vd1 = run(kresult,result,name,gg,rr,bb)
    
    arr = np.column_stack((rr,gg))
   
    B,G = bb/rr,gg
    R0 = rr
    N=3.0001
    drb = 2*R0/B*(beta(((N-3.)/2.),((3.-G)/2.)) - beta(((N-3.)/2.),(3./2.)) * (1+B**2)**((3-N)/2.) * hyp((N-3.)/2.,G/2.,N/2., 1./(1+B**2))) 
    B = 1./rr
    dr1 = 2*R0/B*(beta(((N-3.)/2.),((3.-G)/2.)) - beta(((N-3.)/2.),(3./2.)) * (1+B**2)**((3-N)/2.) * hyp((N-3.)/2.,G/2.,N/2., 1./(1+B**2))) 
    kappa = bb * dr1/drb

    vd2 = kappa * sigma_dm.eval(arr)

    print '%f,%f'%(vd1**0.5, vd2**0.5)
    s1[ii],s2[ii] = vd1**0.5, vd2**0.5

pl.figure()
pl.scatter(s1,s2,c='SteelBlue',s=40)
x = np.linspace(0,500,100)
pl.plot(x,x,'k',ls='--')
pl.show()
