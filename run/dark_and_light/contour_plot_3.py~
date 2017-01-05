import numpy as np, pylab as pl, pyfits as py
from scipy import ndimage 

ch2 = np.loadtxt('/data/ljo31b/EELs/job_applications/inference_pixels.dat')
ch1 = np.loadtxt('/data/ljo31b/EELs/job_applications/inference_apertures.dat')

ch2[:,0] *= 10.
ch1[:,0] /= 1.2

xbins = np.linspace(0,2.2,20)
ybins = np.linspace(1.5,2.8,20)
smooth=1.  # experiment

pl.xlim([xbins[0],xbins[-1]])
pl.ylim([ybins[0],ybins[-1]])

H,x,y = pl.histogram2d(ch1[:,0],ch1[:,1],bins=[xbins,ybins])
H = ndimage.gaussian_filter(H,smooth)
sortH = np.sort(H.flatten())
cumH = sortH.cumsum()

lvl00 = 2*sortH.max()
lvl68 = sortH[cumH>cumH.max()*0.32].min()
lvl95 = sortH[cumH>cumH.max()*0.05].min()
lvl99 = sortH[cumH>cumH.max()*0.003].min()

#pl.contourf(H.T,[lvl95,lvl68],colors='red',alpha=0.4,extent=(xbins[0],xbins[-1],ybins[0],ybins[-1]))
#pl.contourf(H.T,[lvl68,lvl00],colors='red',alpha=0.7,extent=(xbins[0],xbins[-1],ybins[0],ybins[-1]))
pl.contour(H.T,[lvl68,lvl95],colors='Crimson',extent=(xbins[0],xbins[-1],ybins[0],ybins[-1]),linestyles='dashed')


smooth=0.85
H,x,y = pl.histogram2d(ch2[:,0],ch2[:,1],bins=[xbins,ybins])
H = ndimage.gaussian_filter(H,smooth)
sortH = np.sort(H.flatten())
cumH = sortH.cumsum()

lvl00 = 2*sortH.max()
lvl68 = sortH[cumH>cumH.max()*0.32].min()
lvl95 = sortH[cumH>cumH.max()*0.05].min()
lvl99 = sortH[cumH>cumH.max()*0.003].min()

#pl.contourf(H.T,[lvl95,lvl68],colors='blue',alpha=0.4,extent=(xbins[0],xbins[-1],ybins[0],ybins[-1]))
#pl.contourf(H.T,[lvl68,lvl00],colors='blue',alpha=0.7,extent=(xbins[0],xbins[-1],ybins[0],ybins[-1]))
pl.contour(H.T,[lvl68,lvl95],colors='CornflowerBlue',extent=(xbins[0],xbins[-1],ybins[0],ybins[-1]))


pl.plot(-1,-1,color='CornflowerBlue',label='pixel modelling')
pl.plot(-1,-1,color='Crimson',label='aperture modelling')

pl.ylabel('dark halo inner slope',fontsize=30)
pl.xlabel('stellar mass ($10^{11}$M$_{\odot}$)',fontsize=30)
pl.legend(loc='lower right')

pl.show()
