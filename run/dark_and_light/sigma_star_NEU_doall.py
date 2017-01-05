from tools import solarmag
import numpy as np, pylab as pl, pyfits as py
import cPickle
from jeans.makemodel import *
from astLib import astCalc
from imageSim import SBObjects
from itertools import product
from linslens import EELsModels_huge as L, EELsKeckModels as K

def run(kresult,result,name):

    model = L.EELs(result,name)
    kmodel = K.EELs(kresult,result,name)

    model.Initialise()
    kmodel.Initialise()
    zl = lz[name][0]
    scale = astCalc.da(zl)*np.pi/180./3600 * 1e3

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
    Mlum = light2mass(lr,light,1.)
    fac = Mlum[-1]/(10**10)
    Mlum /= fac
    
    # take sigma within the effective radius of the galaxy
    sigma_star = veldisp(r,sb,Mlum,ap=scale)#[-0.5*scale,0.5*scale,-0.5*scale,0.5*scale])
    #print sigma_star
    return sigma_star**0.5


#name = 'J0901'
names = ['J0837','J0901','J0913','J1125','J1144','J1218','J1323','J1347','J1446','J1605','J1606','J1619','J2228']
dir = '/data/ljo31/Lens/LensModels/twoband/'
lz = np.load('/data/ljo31/Lens/LensParams/LensRedshiftsUpdated.npy')[()]

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
            result = np.load(dir+name+'_212')
            kresult = np.load(dir+name+'_Kp_212_lensandgalon')
            print 'J1619'                                                                                                                                                                                                       
        else:
            result = np.load(dir+name+'_211')
            kresult = np.load(dir+name+'_Kp_211')
            print 'here'
 
    sigma_star = run(kresult,result,name)
    print name, sigma_star
    np.savetxt('/data/ljo31b/EELs/phys_models/models/sigma_star_'+name+'.dat',np.array([sigma_star]))


