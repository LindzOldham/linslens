import pyfits as py,numpy as np,pylab as pl
import indexTricks as iT
import colorImage
from scipy import ndimage
from linslens import EELsModels_huge as L

dir = '/data/ljo31/Lens/LensModels/twoband/'
sz = np.load('/data/ljo31/Lens/LensParams/SourceRedshifts.npy')[()]
names = sz.keys()
names.sort()

CI = colorImage.ColorImage()
CI.nonlin = 40.
for name in names:
    if name in ['J0837','J0901','J1218']:#,'J1323']:
        result = np.load(dir+name+'_211')
    elif name == 'J1248':
        continue
    else:
        try:
            result = np.load(dir+name+'_212')
        except:
            result = np.load(dir+name+'_112')
    model = L.EELs(result,name)
    model.Initialise()
    V,I = model.models
    sy,sx = V.shape
    cx,cy = (sx-100.)/2., (sy-100.)/2.
    if name == 'J1606':
        cy-=20
    mid = 0.5*(V+I)
    img = CI.createModel(V[cy:-cy,cx:-cx],mid[cy:-cy,cx:-cx],I[cy:-cy,cx:-cx])
    pl.imshow(img,origin='lower',interpolation='nearest')
    pl.gca().xaxis.set_ticks([])
    pl.gca().yaxis.set_ticks([])
    pl.savefig('/data/ljo31/public_html/Lens/model_colour_images/'+name+'_colour.png')
    pl.show()

