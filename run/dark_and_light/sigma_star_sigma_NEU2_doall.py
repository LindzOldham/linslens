from tools import solarmag
import numpy as np, pylab as pl, pyfits as py
import cPickle
from jeans.makemodel import *
from astLib import astCalc
from imageSim import SBObjects
from itertools import product
from linslens import EELsModels_huge as L, EELsKeckModels as K
import lenslib

def run(kresult,result,name):

    model = L.EELs(result,name)
    model.Initialise()
    
    if kresult is not None:
        kmodel = K.EELs(kresult,result,name)
        kmodel.Initialise()
        fits = kmodel.fits

    zl,zs = lz[name][0],sz[name][0]
    scale = astCalc.da(zl)*np.pi/180./3600 * 1e3
    Rein_tot = cat[name]['Lens 1 b']*0.05*scale

    # get relative contributions from the K band
    fracs = [1.]
    Fracs = [1.,0.]
    if model.galno == 2.:
        galk1,galk2 = kmodel.gals
        ii = np.where(kmodel.fits[0:2]!=0.)
        if len(ii[0])<len(kmodel.gals):
            print 'one of the components has no flux!'
            print 'component ', '%.2f'%np.where(kmodel.fits[0:2]==0), 'has no flux so we are removing it'
            model.galno = 1.
            model.gals = [model.gals[ii[0]]]
            # preserve fractions
            Fracs = [0,0]
            Fracs[ii[0]] = 1.
        else:
            #print galk1.re,galk2.re
            mk1,mk2 = galk1.getMag(kmodel.fits[0],0), galk2.getMag(kmodel.fits[1],0)
            Fk1,Fk2 = 10**(-0.4*mk1), 10**(-0.4*mk2) 
            f1,f2 = Fk1/(Fk1+Fk2), Fk2/(Fk1+Fk2)
            print 'there are two components with relative fluxes ' '%.2f'%f1, 'and ', '%.2f'%f2
            fracs = [f1,f2]
            Fracs = fracs

    print name, Fracs
    # build lenses
    if model.galno == 1.:
        gal = model.gals[0]
        lens1 = PROFILES.Sersic(x=0,y=0,n=gal.n,re=gal.re*0.05*scale,q=gal.q,pa=gal.pa,b=50.,zl=zl,zs=zs)
        lens1.setbFromMass(1e12)
        lenses = [lens1]
    else:
        gal1,gal2 = model.gals
        lens1 = PROFILES.Sersic(x=0,y=0,n=gal1.n,re=gal1.re*0.05*scale,q=gal1.q,pa=gal1.pa,b=50.,zl=zl,zs=zs)
        lens2 = PROFILES.Sersic(x=0,y=0,n=gal2.n,re=gal2.re*0.05*scale,q=gal2.q,pa=gal2.pa,b=50.,zl=zl,zs=zs)
        # can set b from total mass
        lens1.setbFromMass(1e12)
        lens2.setbFromMass(1e12)
        if name == 'J2228':
            lens1.setbFromMass(5e12)
            lens2.setbFromMass(5e12)
        lenses = [lens1, lens2]

    r = np.logspace(-5,5,2501)
    lr = r[:-100]
    mass = lr*0.
    sigma = lr*0.

    for ii in range(len(lenses)):
        lens = lenses[ii]
        sigma += fracs[ii] * lens.sigma(lr)
        #print fracs[ii], lens.sigma(lr)
    
    sigmod = splrep(lr,sigma*2.*np.pi*lr)
    M_ein = splint(0,Rein_tot,sigmod)
    print M_ein
    return M_ein



#name = 'J0901'
names = ['J0837','J0901','J0913','J1125','J1144','J1218','J1323','J1347','J1446','J1605','J1606','J1619','J2228']
dir = '/data/ljo31/Lens/LensModels/twoband/'
lz = np.load('/data/ljo31/Lens/LensParams/LensRedshiftsUpdated.npy')[()]
sz = np.load('/data/ljo31/Lens/LensParams/SourceRedshiftsUpdated.npy')[()]
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
            result = np.load(dir+name+'_212')
            kresult = np.load(dir+name+'_Kp_212_lensandgalon')
            print 'J1619'                                                                                                                                                                                                       
        else:
            result = np.load(dir+name+'_211')
            kresult = np.load(dir+name+'_Kp_211')
            print 'here'
    if name == 'J1144':
        result = np.load('/data/ljo31/Lens/J1144/twoband_darkandlightprep')
        kresult = None

    sigma_star = run(kresult,result,name)
    print name, sigma_star
    np.savetxt('/data/ljo31b/EELs/phys_models/models/sigma_star_sigma'+name+'.dat',np.array([sigma_star]))


