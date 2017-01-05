import numpy as np, pylab as pl, pyfits as py
from tools import gus_plotting as g

res1 = np.load('/data/ljo31b/EELs/aperture_mass_measure_inference')
res2 = np.load('/data/ljo31b/EELs/galsub/emceeruns/J0901_parametric_DPL_0')

lp1,trace1,dic1,_ = res1
lp2,trace2,dic2,_ = res2

ch1 = g.changechain(trace1[2000:,0])
ch2 = g.changechain(trace2[:,0])

#g.triangle_plot(ch1)
#g.triangle_plot(ch2)

#pl.show()

# save as files to use Matt's corner plotter

ch1 = np.column_stack((dic1['ML'][:,0].ravel(),dic1['gamma'][:,0].ravel()))
ch2 = np.column_stack((dic2['stellar mass'][:,0].ravel(),dic2['Lens 1 gamma'][:,0].ravel()))

#np.savetxt('/data/ljo31b/EELs/job_applications/inference_apertures.dat',ch1)
#np.savetxt('/data/ljo31b/EELs/job_applications/inference_pixels.dat',ch2)

from scipy import ndimage 

ch2[:,0] *= 10.
# y is gamma
# x is ML
xbins = np.linspace(0,0.2,20)
ybins = np.linspace(2,2.8,20)
smooth=1.  # experiment

pl.xlim([xbins[0],xbins[-1]])
pl.ylim([ybins[0],ybins[-1]])

H,x,y = pl.histogram2d(ch2[:,0],ch2[:,1],bins=[xbins,ybins])
H = ndimage.gaussian_filter(H,smooth)
sortH = np.sort(H.flatten())
cumH = sortH.cumsum()

lvl00 = 2*sortH.max()
lvl68 = sortH[cumH>cumH.max()*0.32].min()
lvl95 = sortH[cumH>cumH.max()*0.05].min()
lvl99 = sortH[cumH>cumH.max()*0.003].min()

pl.contourf(H.T,[lvl95,lvl68],colors='blue',alpha=0.4,extent=(xbins[0],xbins[-1],ybins[0],ybins[-1]))
pl.contourf(H.T,[lvl68,lvl00],colors='blue',alpha=0.7,extent=(xbins[0],xbins[-1],ybins[0],ybins[-1]))
pl.contour(H.T,[lvl68,lvl95],colors='blue',extent=(xbins[0],xbins[-1],ybins[0],ybins[-1]))

pl.show()

#import iCornerPlotter as i
#dir = '/data/ljo31b/EELs/job_applications/'
#i.CornerPlotter([dir+'inference_pixels.dat,blue',dir+'inference_apertures.dat,red'])
