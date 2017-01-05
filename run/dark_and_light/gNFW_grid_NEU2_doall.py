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

def run(kresult,result,name):

    model = L.EELs(result,name)
    model.Initialise()
    
    if kresult is not None:
        kmodel = K.EELs(kresult,result,name)
        kmodel.Initialise()
        fits = kmodel.fits


    zl,zs = lz[name][0],sz[name][0]
    scale = astCalc.da(zl)*np.pi/180./3600 * 1e3
    
    sig_crit = lenslib.sig_crit(zl,zs) # solar masses per Mpc^2
    sig_crit /= (1e3)**2. # solar masses per kpc^2

    Rein = cat[name]['Lens 1 b']*0.05*scale

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
    arr = [[r,r0grid[m],gammagrid[n],zl,zs] for m,n in product(range(len(r0grid)) ,range(len(gammagrid)))]
    print 'one'
    Mdm = np.zeros((lr.size,r0grid.size,gammagrid.size))
    sigma_star = np.zeros((r0grid.size,gammagrid.size))
    out = p.map(gridgNFW,arr)
    for i in range(len(arr)):
        Mdm[:,idx[i][0],idx[i][1]] = out[i]
    print 'two'
    # also multiprocess sigma star!
    arr = [[r,sb,Mdm[:,idx[i][0],idx[i][1]],scale] for i in range(len(arr))] #not stellar component
    out = p.map(gridveldisp,arr)
    print 'three'
    for i in range(len(arr)):
        print i
        sigma_star[idx[i][0],idx[i][1]] = out[i]

    vd = sigma_star

    # SB
    print name
    return vd


phot = py.open('/data/ljo31/Lens/LensParams/Phot_1src_huge_new_new.fits')[1].data
names = phot['name']
lz = np.load('/data/ljo31/Lens/LensParams/LensRedshiftsUpdated.npy')[()]
sz = np.load('/data/ljo31/Lens/LensParams/SourceRedshiftsUpdated.npy')[()]
cat = np.load('/data/ljo31/Lens/LensParams/Structure_lensgals_2src.npy')[0]

#cat_nugg = []

dir = '/data/ljo31/Lens/LensModels/twoband/'

r0grid = np.arange(1,600.,2.5)
gammagrid = np.arange(0,2.9,0.05)
idx = list(product(range(len(r0grid)) ,range(len(gammagrid))))
p = Pool(8)

#name = 'J0901'
names = ['J0837','J0901','J0913','J1125','J1144','J1218','J1323','J1347','J1446','J1605','J1606','J1619','J2228']

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
    if name != 'J1144':
        continue

    vd = run(kresult,result,name)
    
    print name
    # build interplator
    ax = {}
    ax[0] = splrep(r0grid,np.arange(r0grid.size),k=1,s=0)
    ax[1] = splrep(gammagrid,np.arange(gammagrid.size),k=1,s=0)
    obj = ndinterp.ndInterp(ax,vd,order=3)
    np.save('/data/ljo31b/EELs/phys_models/models/interpolators/gNFW_aperture_mass_measure_'+name,obj)


