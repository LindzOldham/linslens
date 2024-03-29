from tools import solarmag
import numpy as np, pylab as pl, pyfits as py
import cPickle
from jeans.makemodel import *
from astLib import astCalc
from imageSim import SBObjects
from itertools import product
from linslens import EELsModels_huge as L

def run(result,name):

    model = L.EELs(result,name)
    model.Initialise()
    RE,_ = model.GetSourceSize(kpc=True)
    fits = model.fits
    scale = model.scale
    z = model.z
    r = np.logspace(-5,5,1501)

    if name in ['J0901','J1218','J1323']:
        gal1 = model.srcs[0]
        gal1.re *= 0.05*scale
        sb = fits[0][-2]*gal1.eval(r)

    elif name == 'J0837':
        gal1 = model.srcs[0]
        gal1.re *= 0.05*scale
        sb = fits[0][-3]*gal1.eval(r) # the other one is the dust lane!

    else:
        gal1,gal2 = model.srcs
        gal1.re *= 0.05*scale
        gal2.re *= 0.05*scale
        sb = fits[0][-3]*gal1.eval(r) + fits[0][-2]*gal2.eval(r)# in image units, but is normalised by the total mass
    
    
    # stellar mass profile
    lr,light = deproject(r,sb)
    #Mlum = light2mass(lr,light,1.)
    #fac = Mlum[-1]/(10**10)
    #Mlum /= fac
    #Mlum_model = splrep(lr,Mlum)
    light_model = splrep(lr,light)

    sr = lr[:-100]
    Sigma = sr*0.
    for i in range(sr.size):
        #r = lr[i:]
        R = lr[i]
        rr = np.logspace(-5,0.5*np.log10(lr[-1]**2-R**2),1001)
        y = (rr**2+R**2)**0.5
        f = splev(y,light_model)
        model = splrep(rr,f)
        #model = splrep(r,Mlum[i:]*r/(r**2-R**2)**0.5)
        Sigma[i] = 2.*splint(rr[0],rr[-1],model)

    model = splrep(sr,Sigma*2.*np.pi*sr)
    Sig_einstein = splint(0,rein,model)
    #Sig_einstein = splev(rein,model)
  
    print Sig_einstein
    return Sig_einstein


name = 'J0901'

dir = '/data/ljo31/Lens/LensModels/twoband/'

lz = np.load('/data/ljo31/Lens/LensParams/LensRedshiftsUpdated.npy')[()]
sz = np.load('/data/ljo31/Lens/LensParams/SourceRedshiftsUpdated.npy')[()]
zl,zs = lz[name][0],sz[name][0]
scale = astCalc.da(zl)*np.pi/180./3600 * 1e3

cat = np.load('/data/ljo31/Lens/LensParams/Structure_lensgals_2src.npy')[0]
rein = cat[name]['Lens 1 b']*0.05*scale

if name in ['J0913','J1125','J1144','J1347','J1446','J1605','J1619','J2228']:
    result = np.load(dir+name+'_212')
elif name in ['J0837','J0901','J1218','J1323']:
    result = np.load(dir+name+'_211')
else:
    print 'missing eel!'

sigma_star = run(result,name)

np.savetxt('/data/ljo31b/EELs/phys_models/models/sigma_star_sigma_einstein_'+name+'.dat',np.array([sigma_star]))


