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
    Mlum = light2mass(lr,light,1.)
    fac = Mlum[-1]/(10**10)
    #Mlum /= fac
    
    # take sigma within the effective radius of the galaxy
    sigma_star = veldisp(r,sb,Mlum,ap=RE)
    print sigma_star
    return sigma_star


name = 'J0901'

dir = '/data/ljo31/Lens/LensModels/twoband/'



if name in ['J0913','J1125','J1144','J1347','J1446','J1605','J1619','J2228']:
    result = np.load(dir+name+'_212')
elif name in ['J0837','J0901','J1218','J1323']:
    result = np.load(dir+name+'_211')
else:
    print 'missing eel!'

sigma_star = run(result,name)

np.savetxt('/data/ljo31b/EELs/phys_models/models/sigma_star_'+name+'.dat',np.array([sigma_star]))


