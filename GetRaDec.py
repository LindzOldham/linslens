import glob
import pyfits as py, pylab as pl, numpy as np

sz = np.load('/data/ljo31/Lens/LensParams/SourceRedshiftsUpdated.npy')[()]
lz = np.load('/data/ljo31/Lens/LensParams/LensRedshiftsUpdated.npy')[()]

names = sz.keys()
names.sort()

dir = '/data/ljo31b/EELs/esi/raw/'
for name in names:
    #if name == 'J1248':
    #    continue
    try:
        file = glob.glob('/data/ljo31b/EELs/esi/*/EEL_'+name+'*.fits')[0]
    except:
        file = glob.glob(dir+'*/run/EEL_'+name+'*_bgsub.fits')[0]
    hdr = py.open(file)[0].header
    ra,dec = hdr['RA'],hdr['DEC']
    print name, '& ', ra, '& ', dec, '& ', '%.4f'%lz[name][0], '& ', '%.4f'%sz[name][0], r'\\'



for name in names:
    #if name == 'J1248':
    #    continue
    try:
        file = glob.glob('/data/ljo31b/EELs/esi/*/EEL_'+name+'*.fits')[0]
    except:
        file = glob.glob(dir+'*/run/EEL_'+name+'*_bgsub.fits')[0]
    hdr = py.open(file)[0].header
    print hdr['DATE-OBS']
