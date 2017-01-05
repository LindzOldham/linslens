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

def run(name):

    zl,zs = lz[name][0],sz[name][0]
    scale = astCalc.da(zl)*np.pi/180./3600 * 1e3
    Rein_tot = cat[name]['Lens 1 b']*0.05*scale

    r = np.logspace(-5,5,2501)
    lr = r[:-100]

    arr = [[r,r0grid[m],gammagrid[n],zl,zs,Rein_tot] for m,n in product(range(len(r0grid)) ,range(len(gammagrid)))]

    Mdm = np.zeros((r0grid.size,gammagrid.size))
    out = p.map(gridgNFW_sigma,arr)
    for i in range(len(arr)):
        Mdm[idx[i][0],idx[i][1]] = out[i]
    
    
    print name, Mdm
    return Mdm


phot = py.open('/data/ljo31/Lens/LensParams/Phot_1src_huge_new_new.fits')[1].data
names = phot['name']


dir = '/data/ljo31/Lens/LensModels/twoband/'

r0grid = np.arange(2.,600.,5)
gammagrid = np.arange(0,2.9,0.05)
idx = list(product(range(len(r0grid)) ,range(len(gammagrid))))
p = Pool(8)

lz = np.load('/data/ljo31/Lens/LensParams/LensRedshiftsUpdated.npy')[()]
sz = np.load('/data/ljo31/Lens/LensParams/SourceRedshiftsUpdated.npy')[()]

cat = np.load('/data/ljo31/Lens/LensParams/Structure_lensgals_2src.npy')[0]


for name in names:
    
    '''if name == 'J1248':
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
            print 'here' '''

    vd = run(name)
    
    # build interplator
    ax = {}
    ax[0] = splrep(r0grid,np.arange(r0grid.size),k=1,s=0)
    ax[1] = splrep(gammagrid,np.arange(gammagrid.size),k=1,s=0)
    obj = ndinterp.ndInterp(ax,vd,order=3)
    pl.figure()
    pl.imshow((obj.z-obj.spline)/obj.z,interpolation='nearest',origin='lower',aspect='auto')
    pl.colorbar()
    pl.show()

    np.save('/data/ljo31b/EELs/phys_models/models/interpolators/gNFW_aperture_mass_measure_Sig_einstein_'+name,obj)


