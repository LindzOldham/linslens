import numpy as np, pylab as pl, pyfits as py
from scipy import ndimage 
import glob
from astLib import astCalc
import sys

name = sys.argv[1] 
print name
#name = 'J0901'

outFile = '/data/ljo31b/EELs/aperture_mass_measure_inference_gNFW_'+name
result = np.load(outFile)
lp,trace,dic,_ = result
ml, gamma = dic['ML'][500:,0].ravel()*10, dic['gamma'][500:,0].ravel()

dir = '/data/ljo31b/EELs/galsub/emceeruns/'
file = glob.glob(dir+name+'_parametric_DPL__VI_*')
file.sort()
f = file[-1]
print f 
result_DPL = np.load(f)

lz = np.load('/data/ljo31/Lens/LensParams/LensRedshiftsUpdated.npy')[()]
sz = np.load('/data/ljo31/Lens/LensParams/SourceRedshiftsUpdated.npy')[()]
zl,zs = lz[name][0],sz[name][0]
scale = astCalc.da(zl)*np.pi/180./3600 * 1e3

lp_DPL,trace_DPL,dic_DPL,_ = result_DPL
ml_dpl, gamma_dpl = dic_DPL['stellar mass'][:,0].ravel()*10, dic_DPL['Lens 1 gamma'][:,0].ravel()

xbins = np.linspace(0,4,20)
ybins = np.linspace(0,4,20)
smooth=1.  # experiment

pl.xlim([xbins[0],xbins[-1]])
pl.ylim([ybins[0],ybins[-1]])

H,x,y = pl.histogram2d(gamma,ml,bins=[xbins,ybins])
H = ndimage.gaussian_filter(H,1)
sortH = np.sort(H.flatten())
cumH = sortH.cumsum()

lvl00 = 2*sortH.max()
lvl68 = sortH[cumH>cumH.max()*0.32].min()
lvl95 = sortH[cumH>cumH.max()*0.05].min()
lvl99 = sortH[cumH>cumH.max()*0.003].min()

#pl.contourf(H.T,[lvl95,lvl68],colors='red',alpha=0.4,extent=(xbins[0],xbins[-1],ybins[0],ybins[-1]))
#pl.contourf(H.T,[lvl68,lvl00],colors='red',alpha=0.7,extent=(xbins[0],xbins[-1],ybins[0],ybins[-1]))
pl.contour(H.T,[lvl95,lvl68],colors='Crimson',extent=(xbins[0],xbins[-1],ybins[0],ybins[-1]),linestyles='dashed')

smooth=0.85
smooth = 0.75
smoot=0.5
H,x,y = pl.histogram2d(gamma_dpl,ml_dpl,bins=[xbins,ybins])
H = ndimage.gaussian_filter(H,smooth)
sortH = np.sort(H.flatten())
cumH = sortH.cumsum()

lvl00 = 2*sortH.max()
lvl68 = sortH[cumH>cumH.max()*0.32].min()
lvl95 = sortH[cumH>cumH.max()*0.05].min()
lvl99 = sortH[cumH>cumH.max()*0.003].min()

#pl.contourf(H.T,[lvl95,lvl68],colors='blue',alpha=0.4,extent=(xbins[0],xbins[-1],ybins[0],ybins[-1]))
#pl.contourf(H.T,[lvl68,lvl00],colors='blue',alpha=0.7,extent=(xbins[0],xbins[-1],ybins[0],ybins[-1]))
pl.contour(H.T,[lvl95,lvl68],colors='CornflowerBlue',extent=(xbins[0],xbins[-1],ybins[0],ybins[-1]))


pl.plot(-1,-1,color='CornflowerBlue',label='pixel modelling')
pl.plot(-1,-1,color='Crimson',label='aperture modelling')

pl.xlabel('dark halo inner slope',fontsize=30)
pl.ylabel('stellar mass ($10^{11}$M$_{\odot}$)',fontsize=30)
pl.legend(loc='upper left')
pl.title(name+': gNFW')
pl.savefig('/data/ljo31/public_html/Lens/dark_and_light/compare_profiles_new/contour_gNFW_VI_%s.pdf'%name)
#pl.show()

# need to make corner plot as we now have r_s as well!
g1,r1,ml1 = dic['gamma'][500:,0].ravel(), dic['RS'][500:,0].ravel(),dic['ML'][500:,0].ravel()*10.
g2,r2,ml2 = dic_DPL['Lens 1 gamma'][:,0].ravel(), dic_DPL['Lens 1 rs'][:,0].ravel()*0.05*scale,dic_DPL['stellar mass'][:,0].ravel()*10.

ch1 = np.column_stack((g1,r1,ml1))
ch2 = np.column_stack((g2,r2,ml2))

from corner_plot import corner_plot
corner_plot(ch1,axis_labels=['$\gamma$', '$r_s$', 'M$_{\star}$'])
pl.savefig('/data/ljo31/public_html/Lens/dark_and_light/compare_profiles_new/%s_gNFW_aperture.pdf'%name)
corner_plot(ch2,axis_labels=['$\gamma$', '$r_s$', 'M$_{\star}$'],cmap='Purples')
pl.savefig('/data/ljo31/public_html/Lens/dark_and_light/compare_profiles_new/%s_gNFW_pixels.pdf'%name)
pl.show()

np.savetxt('APERTURE.dat',ch1)
np.savetxt('PIXEL.dat',ch2)
