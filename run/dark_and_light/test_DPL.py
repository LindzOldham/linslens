from MWApython.pylens import MassModels,pylens
import numpy as np, pylab as pl
import indexTricks as iT
from imageSim import SBObjects

B,Q,ETA = 10.,1.,1.
X,Y = 0.1,0.1
ETA1,ETA2,RS = 2.,0.,6.
l1 = MassModels.PowerLaw('l1',{'x':X,'y':Y,'b':B,'eta':ETA,'q':Q,'pa':0.})
l2 = MassModels.DPL('l2',{'x':X,'y':Y,'b':B,'eta1':ETA1,'eta2':ETA2,'q':Q,'pa':0.,'rs':RS})
l3 = MassModels.sGNFW('l3',{'x':X,'y':Y,'b':B,'rs':1000000,'eta':2})

s = SBObjects.Sersic('s',{'x':X+1,'y':Y-1,'re':25,'n':4,'q':1,'pa':0})
y,x = iT.coords((81,81))-40

import time
st = time.time()
x1,y1 = pylens.getDeflections(l1,[x,y])
st1 = time.time()
print 'time 1', st1-st
x2,y2 = pylens.getDeflections(l2,[x,y])
st2 = time.time()
print 'time 2', st2-st1
x3,y3 = pylens.getDeflections(l3,[x,y])
st3 = time.time()
print 'time 3', st3-st2

s1 = s.pixeval(x1,y1,csub=31)
s2 = s.pixeval(x2,y2,csub=31)
s3 = s.pixeval(x3,y3,csub=31)


pl.figure(figsize=(16,8))
pl.subplot(231)
pl.imshow(s1,interpolation='nearest',origin='lower')
pl.colorbar()
pl.subplot(232)
pl.imshow(s2,interpolation='nearest',origin='lower')
pl.colorbar()
pl.subplot(233)
pl.imshow(s3,interpolation='nearest',origin='lower')
pl.colorbar()
## residuals
pl.subplot(234)
pl.imshow((s1-s1)/s1,interpolation='nearest',origin='lower')
pl.colorbar()
pl.subplot(235)
pl.imshow((s1-s2)/s1,interpolation='nearest',origin='lower')
pl.colorbar()
pl.subplot(236)
pl.imshow((s1-s3)/s1,interpolation='nearest',origin='lower')
pl.colorbar()


pl.show()

# I'm not sure I like this because it is very different from the elliptical power law where they are meant to be the same. This potential thing...

# Matt's sGNFW is 5 times slower than Tom's DPL. And seems to get a different result? Has it not been debugged?
