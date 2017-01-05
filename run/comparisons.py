import numpy as np, pyfits as py, pylab as pl
from astLib import astCalc


photnew = py.open('/data/ljo31/Lens/LensParams/Phot_2src_huge_new.fits')[1].data
photold = py.open('/data/ljo31/Lens/LensParams/Phot_2src_new.fits')[1].data

pl.figure()
pl.scatter(photnew['Re v'],photnew['Re v']-photold['Re v'],c='SteelBlue',s=40)
pl.axhline(0)
pl.xlabel('$r_e$ (pixels, new)')
pl.ylabel('new-old')
pl.title('212 models')
pl.savefig('/data/ljo31/public_html/Lens/big_boxes/re_212.pdf')

pl.figure()
pl.scatter(photnew['mag v'],photnew['mag v']-photold['mag v'],c='SteelBlue',s=40)
pl.axhline(0)
pl.xlabel('$V$ (new)')
pl.ylabel('new-old')
pl.title('212 models')
pl.savefig('/data/ljo31/public_html/Lens/big_boxes/magv_212.pdf')


pl.figure()
pl.scatter(photnew['mag i'],photnew['mag i']-photold['mag i'],c='SteelBlue',s=40)
pl.axhline(0)
pl.xlabel('$I$ (new)')
pl.ylabel('new-old')
pl.title('212 models')
pl.savefig('/data/ljo31/public_html/Lens/big_boxes/magi_212.pdf')

pl.figure()
pl.scatter(photnew['mu i'],photnew['mu i']-photold['mu i'],c='SteelBlue',s=40)
pl.axhline(0)
pl.xlabel('$\mu_I$ (new)')
pl.ylabel('new-old')
pl.title('212 models')
pl.savefig('/data/ljo31/public_html/Lens/big_boxes/mui_212.pdf')

pl.figure()
pl.scatter(photnew['mu v'],photnew['mu v']-photold['mu v'],c='SteelBlue',s=40)
pl.axhline(0)
pl.xlabel('$\mu_V$ (new)')
pl.ylabel('new-old')
pl.title('212 models')
pl.savefig('/data/ljo31/public_html/Lens/big_boxes/muv_212.pdf')

# conclusion -- changes ARE NOT systematic

# also n of lens gals
