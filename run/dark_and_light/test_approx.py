import numpy as np, pylab as pl, pyfits as py
from MWApython.pylens import MassModels
import indexTricks as iT
import time

m1 = MassModels.DPL('TC',{'x':0,'y':0,'b':10,'rs':100.1,'eta1':2.1,'eta2':0.9,'q':0.8,'pa':43.})
m2 = MassModels.sGNFW('MA',{'x':0,'y':0,'b':10,'rs':10,'eta':2.})
m3 = MassModels.DPL_2('TC',{'x':0,'y':0,'b':10,'rs':100.1,'eta':2.1,'q':0.8,'pa':43.})

y,x = iT.coords((100,100))
t1=time.time()
x1,y1 = m1.deflections(x,y)
t2=time.time()
x2,y2 = m2.deflections(x,y)
t3=time.time()
x3,y3 = m3.deflections(x,y)
t4=time.time()
##
x1,y1 = m1.deflections(x,y)
t5=time.time()
x2,y2 = m2.deflections(x,y)
t6=time.time()
x3,y3 = m3.deflections(x,y)
t7=time.time()

print t7-t6,t4-t3
print t6-t5,t3-t2
print t5-t4,t2-t1

'''pl.figure(figsize=(23,8))
pl.subplot(131)
pl.imshow(x1,interpolation='nearest',origin='lower')
pl.colorbar()
pl.subplot(132)
pl.imshow(x3,interpolation='nearest',origin='lower')
pl.colorbar()
pl.subplot(133)
pl.imshow(((x3-x1)/x3),interpolation='nearest',origin='lower',vmin=0.001,vmax=0.0011)
pl.colorbar()
pl.show()'''
