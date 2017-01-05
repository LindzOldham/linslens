from linslens import EELsModels_huge as L
import numpy as np
import pylab as pl
from linslens import EELsKeckModels as K
from linslens.Plotter import *
from linslens.GrabImages_huge import *
from scipy import ndimage
import glob

def clip(arr,nsig=3.5):
    a = arr.flatten()
    m,s,l = a.mean(),a.std(),a.size
    while 1:
        a = a[abs(a-m)<s*nsig]
        if a.size==l:
            return m,s
        m,s,l = a.mean(),a.std(),a.size

cuts = [['J0901',[25,25,25,25]],['J1125',[50,50,50,50]],['J2228',[45,45]]]
cuts = dict(cuts)

def MakeMap(name):
    # get background noise from old sigma map
    # V BAND
    sigfile = glob.glob('/data/ljo31/Lens/'+name+'/'+bands[name]+'_noisemap_*.fits')
    if name in ['J1323','J1347']:
        sigfile = glob.glob('/data/ljo31/Lens/'+name+'/*'+bands[name]+'_noisemap_huge.fits')
    elif name == 'J2228':
        sigfile = glob.glob('/data/ljo31/Lens/'+name+'/F606W_noisemap.fits')
    sigma = py.open(sigfile[0])[0].data[0,0]
    poisson = sigma**2.

    im = py.open('/data/ljo31b/EELs/galsub/images/'+name+'_maxlnL.fits')[1].data
    if name in ['J0913','J1125']:
        whtfile = glob.glob('/data/ljo31/Lens/'+name+'/*'+bands[name]+'_wht_cutout_huge2.fits')[0]
    elif name == 'J2228':
        whtfile = glob.glob('/data/ljo31/Lens/'+name+'/*'+bands[name]+'_wht_cutout.fits')[0]
    else:
        whtfile = glob.glob('/data/ljo31/Lens/'+name+'/*'+bands[name]+'_wht_cutout_huge.fits')[0]
    wht = py.open(whtfile)[0].data
    if name in cuts.keys():
        cut = cuts[name]
        print len(cut)
        if len(cut)==4:
            print wht.shape
            wht = wht[cut[0]:-cut[1],cut[2]:-cut[3]]
            print wht.shape
        elif len(cut)==2:
            wht = wht[cut[0]:-cut[1]]
        else:
            print 'huh?'

    print im.shape, wht.shape
    smooth = ndimage.gaussian_filter(im,0.7)
    print smooth.shape
    print sigma, poisson
    noisemap = np.where((smooth>0.7*sigma)&(im>0),im/wht+poisson, poisson)**0.5
    
    ## get rid of nans
    ii = np.where(np.isnan(noisemap)==True)
    noisemap[ii] = np.amax(noisemap[np.isnan(noisemap)==False])
    noisemap_v = noisemap.copy()

    ### I BAND
    sigfile = glob.glob('/data/ljo31/Lens/'+name+'/F814W_noisemap_*.fits')
    if name in ['J1323','J1347']:
        sigfile = glob.glob('/data/ljo31/Lens/'+name+'/*F814W_noisemap_huge.fits')
    elif name == 'J2228':
        sigfile = glob.glob('/data/ljo31/Lens/'+name+'/F814W_noisemap.fits')
    sigma = py.open(sigfile[0])[0].data[0,0]
    poisson = sigma**2.

    im = py.open('/data/ljo31b/EELs/galsub/images/'+name+'_maxlnL.fits')[2].data
    if name in ['J0913','J1125']:
        whtfile = glob.glob('/data/ljo31/Lens/'+name+'/*'+bands[name]+'_wht_cutout_huge2.fits')[0]
    elif name == 'J2228':
        whtfile = glob.glob('/data/ljo31/Lens/'+name+'/*'+'F814W_wht_cutout.fits')[0]
    else:
        whtfile = glob.glob('/data/ljo31/Lens/'+name+'/*'+'F814W_wht_cutout_huge.fits')[0]
    wht = py.open(whtfile)[0].data
    if name in cuts.keys():
        cut = cuts[name]
        print len(cut)
        if len(cut)==4:
            print wht.shape
            wht = wht[cut[0]:-cut[1],cut[2]:-cut[3]]
            print wht.shape
        elif len(cut)==2:
            wht = wht[cut[0]:-cut[1]]
        else:
            print 'huh?'

    smooth = ndimage.gaussian_filter(im,0.7)
    noisemap = np.where((smooth>0.7*sigma)&(im>0),im/wht+poisson, poisson)**0.5
    
    ## get rid of nans
    ii = np.where(np.isnan(noisemap)==True)
    noisemap[ii] = np.amax(noisemap[np.isnan(noisemap)==False])
    noisemap_i = noisemap.copy()


    if name == 'J2228':
        pl.figure(figsize=(10,7))
        pl.subplot(121)
        pl.imshow(noisemap_v,origin='lower',interpolation='nearest',vmax=0.1)
        pl.colorbar()
        pl.subplot(122)
        pl.imshow(noisemap_i,origin='lower',interpolation='nearest',vmax=0.1)
        pl.colorbar()
        pl.show()


   
    outname = '/data/ljo31b/EELs/galsub/images/'+name+'_maxlnL_sig.fits'
    hdu = py.HDUList()
    phdu = py.PrimaryHDU()
    hdr = phdu.header
    hdr['object'] = name
    vband   = py.ImageHDU(noisemap_v,name=bands[name])
    iband = py.ImageHDU(noisemap_i,name='F814W')
    hdu.append(phdu)
    hdu.append(vband)
    hdu.append(iband)
    hdu.writeto(outname,clobber=True)
    
    

dir = '/data/ljo31/Lens/LensModels/twoband/'
sz = np.load('/data/ljo31/Lens/LensParams/SourceRedshiftsUpdated.npy')[()]
bands = np.load('/data/ljo31/Lens/LensParams/HSTBands.npy')[()]

names = sz.keys()
names.sort()
        
for name in names:
    if name == 'J1248':
        continue
    MakeMap(name)

