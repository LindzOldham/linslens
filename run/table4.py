import numpy as np, pyfits as py
from astLib import astCalc


phot = py.open('/data/ljo31/Lens/LensParams/Phot_1src_huge_new.fits')[1].data
names = phot['name']
masses = np.load('/data/ljo31b/EELs/inference/new/huge/masses_212.npy')
ml,ml_l, ml_u, ms, ms_l, ms_u = masses
dml = np.median((ml_l,ml_u),axis=0)
dms = np.median((ms_l,ms_u),axis=0)

for ii in range(len(names)):
    name = names[ii]
    print name, '& $', '%.2f'%ml[ii], '\pm', '%.2f'%dml[ii], '$ & $', '%.2f'%ms[ii], '\pm','%.2f'% dms[ii], r'$\\'

