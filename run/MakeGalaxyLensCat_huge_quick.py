from linslens import EELsLensModels as L
import numpy as np

def MakeTab(file,name):
    result = np.load(file)
    model = L.EELs(result,name)
    model.Initialise()
    mag = model.GalaxyObsMag()
    Re = model.GalaxySize(kpc=True)
    
    return Re,mag


files = ['/data/ljo31/Lens/J0837/twoband_0', '/data/ljo31/Lens/J0913/twoband_1', '/data/ljo31/Lens/J1144/twoband_1', '/data/ljo31/Lens/J1218/twoband_0', '/data/ljo31/Lens/J1323/twoband_0', '/data/ljo31/Lens/J1347/twoband_0', '/data/ljo31/Lens/J1446/twoband_0', '/data/ljo31/Lens/J1605/twoband_0', '/data/ljo31/Lens/J1606/twoband_0', '/data/ljo31/Lens/J1619/twoband_0', '/data/ljo31/Lens/J2228/twoband_0']

names = ['J0837','J0913','J1144','J1218','J1323','J1347','J1446','J1605','J1606','J1619','J2228']

Res,mags = [],[]
for ii in range(len(names)):
    Re,mag  = MakeTab(files[ii],names[ii])
    Res.append(Re)
    mags.append(mag)

names = np.array(names)
Res = np.array(Res)
mags = np.array(mags)

from astropy.io.fits import *

#print mags, mags.shape, Res, Res.shape
c1 = Column(name='name', format='A5',array=names)
c2 = Column(name='mag v', format='D',array=mags[:,0])
c3 = Column(name='mag i', format='D',array=mags[:,1])
c4 = Column(name='re v', format='D',array=Res[:,0])
c5 = Column(name='re i', format='D',array=Res[:,1])




coldefs = ColDefs([c1,c2,c3,c4,c5])
tbhdu = BinTableHDU.from_columns(coldefs)
tbhdu.writeto('/data/ljo31/Lens/LensParams/Phot_1src_lensgals_huge_quick.fits',clobber=True)

print 'FINISHED 1-SRC MODELS'

