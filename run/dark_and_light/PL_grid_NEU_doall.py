from tools import solarmag
import numpy as np, pylab as pl, pyfits as py
import cPickle
from jeans.makemodel import *
from astLib import astCalc
from imageSim import SBObjects
from itertools import product
from linslens import EELsModels_huge as L, EELsKeckModels as K
import ndinterp
from multiprocessing import Pool
import lenslib

def run(kresult,result,name):

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
    arr = [[lr,etagrid[m],Rein,sig_crit] for m in (range(len(etagrid)))]

    Mdm = np.zeros((lr.size,etagrid.size))
    sigma_star = np.zeros((etagrid.size))
    out = p.map(grid_powerlaw,arr)
    for i in range(len(arr)):
        Mdm[:,idx[i]] = out[i]
    
    # also multiprocess sigma star!
    arr = [[r,sb,Mdm[:,idx[i]],scale] for i in range(len(arr))] # scale = 1 arcsec in kpc
    out = p.map(gridveldisp,arr)
    for i in range(len(arr)):
        sigma_star[idx[i]] = out[i]

    vd = sigma_star**0.5

    # SB
    print name, vd
    return vd


phot = py.open('/data/ljo31/Lens/LensParams/Phot_1src_huge_new_new.fits')[1].data
names = phot['name']
lz = np.load('/data/ljo31/Lens/LensParams/LensRedshiftsUpdated.npy')[()]
sz = np.load('/data/ljo31/Lens/LensParams/SourceRedshiftsUpdated.npy')[()]
cat = np.load('/data/ljo31/Lens/LensParams/Structure_lensgals_2src.npy')[0]

dir = '/data/ljo31/Lens/LensModels/twoband/'

etagrid = np.arange(0.1,2.,0.025)
idx = list(range(len(etagrid)))
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
 
    vd = run(kresult,result,name)
    print name
    # build interplator
    ax = {}
    ax[0] = splrep(etagrid,np.arange(etagrid.size),k=1,s=0)
    obj = ndinterp.ndInterp(ax,vd,order=3)
    np.save('/data/ljo31b/EELs/phys_models/models/interpolators/PL_aperture_mass_measure_'+name,obj)


