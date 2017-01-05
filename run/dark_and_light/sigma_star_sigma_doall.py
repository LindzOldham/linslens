from tools import solarmag
import numpy as np, pylab as pl, pyfits as py
import cPickle
from jeans.makemodel import *
from astLib import astCalc
from imageSim import SBObjects
from itertools import product
from linslens import EELsModels_huge as L, EELsKeckModels as K

def run(kresult,result,name):

    zl = lz[name][0]
    scale = astCalc.da(zl)*np.pi/180./3600 * 1e3
    rein = cat[name]['Lens 1 b']*0.05*scale

    model = L.EELs(result,name)
    kmodel = K.EELs(kresult,result,name)

    model.Initialise()
    kmodel.Initialise()

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
    light_model = splrep(lr,light)

    sr = lr[:-100]
    Sigma = sr*0.
    for i in range(sr.size):
        R = lr[i]
        rr = np.logspace(-5,0.5*np.log10(lr[-1]**2-R**2),1001)
        y = (rr**2+R**2)**0.5
        f = splev(y,light_model)
        model = splrep(rr,f)
        Sigma[i] = 2.*splint(rr[0],rr[-1],model)

    model = splrep(sr,Sigma*2.*np.pi*sr)
    Sig_einstein = splint(0,rein,model)
  
    print Sig_einstein
    return Sig_einstein

names = ['J0837','J0901','J0913','J1125','J1144','J1218','J1323','J1347','J1446','J1605','J1606','J1619','J2228']
dir = '/data/ljo31/Lens/LensModels/twoband/'
lz = np.load('/data/ljo31/Lens/LensParams/LensRedshiftsUpdated.npy')[()]
cat = np.load('/data/ljo31/Lens/LensParams/Structure_lensgals_2src.npy')[0]

for name in names:
    if name == 'J1248':
        continue
    try:
        result = np.load(dir+name+'_212')
        kresult = np.load(dir+name+'_Kp_212')
    except:
        if name == 'J1347':
            result = np.load(dir+name+'_112')
            kresult = np.load('/data/ljo31/Lens/J1347/twoband_Kp_112_2')
        elif name == 'J1619':
            #result = np.load(dir+name+'_212')
            #kresult = np.load(dir+name+'_Kp_212_lensandgalon')
            #print 'J1619'                                                                                                                                                                                                       
            continue
        else:
            result = np.load(dir+name+'_211')
            kresult = np.load(dir+name+'_Kp_211')
            print 'here'

    sigma_star = run(kresult,result,name)

    np.savetxt('/data/ljo31b/EELs/phys_models/models/sigma_star_sigma_einstein_'+name+'.dat',np.array([sigma_star]))



