from tools import solarmag
import numpy as np, pylab as pl, pyfits as py
import cPickle
from jeans.makemodel import *
from astLib import astCalc
from imageSim import SBObjects
from itertools import product
from linslens import EELsModels_huge as L
import ndinterp
from multiprocessing import Pool

def run(result,name):

    # make model and extract useful properties
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
    # just get mass and project it

    lr,light = deproject(r,sb)

    arr = [[lr,r0grid[m],gammagrid[n],rein] for m,n in product(range(len(r0grid)) ,range(len(gammagrid)))]

    Mdm = np.zeros((r0grid.size,gammagrid.size))
    out = p.map(gridgNFW3_sigma,arr)
    for i in range(len(arr)):
        Mdm[idx[i][0],idx[i][1]] = out[i]
    
    
    print name
    return Mdm


phot = py.open('/data/ljo31/Lens/LensParams/Phot_1src_huge_new_new.fits')[1].data
names = phot['name']

cat_nugg = []

dir = '/data/ljo31/Lens/LensModels/twoband/'

r0grid = np.arange(10.,600.,5)
gammagrid = np.arange(0,2.9,0.05)
idx = list(product(range(len(r0grid)) ,range(len(gammagrid))))
p = Pool(8)

name = 'J0901'

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

vd = run(result,name)
    
# build interplator
ax = {}
ax[0] = splrep(r0grid,np.arange(r0grid.size),k=1,s=0)
ax[1] = splrep(gammagrid,np.arange(gammagrid.size),k=1,s=0)
obj = ndinterp.ndInterp(ax,vd,order=3)
np.save('/data/ljo31b/EELs/phys_models/models/interpolators/gNFW_aperture_mass_measure_Sig_einstein_'+name,obj)


