import numpy as np, pylab as pl, pyfits as py
from scipy import ndimage 
import glob
from astLib import astCalc
import sys


def plot(name):    
    dir = '/data/ljo31b/EELs/galsub/emceeruns/'
    file = glob.glob(dir+name+'_parametric_DPL__VI_*')
    file.sort()
    f = file[-1]
    print f 
    if name == 'J0837':
        f = dir+'darwin/J0837_parametric_DPL_darwin_2_Iband_ctd'

    result_DPL = np.load(f)

    zl,zs = lz[name][0],sz[name][0]
    scale = astCalc.da(zl)*np.pi/180./3600 * 1e3

    lp_DPL,trace_DPL,dic_DPL,_ = result_DPL
    ml_dpl, gamma_dpl = dic_DPL['stellar mass'][:,0].ravel()*10, dic_DPL['Lens 1 gamma'][:,0].ravel()
    
    Mstar = chab_m[name==names]
    Mstar = 10**Mstar / (1e11)

    DML = (ml_dpl-Mstar)/Mstar
    
    xbins = np.linspace(0,3.,20)
    ybins = np.linspace(-1,1.75,30)
    smooth=1.  # experiment

    pl.xlim([xbins[0],xbins[-1]])
    pl.ylim([ybins[0],ybins[-1]])

    smooth=0.85
    smooth = 0.75
    #smooth=0.5
    H,x,y = pl.histogram2d(gamma_dpl,DML,bins=[xbins,ybins])
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

    Mstar = chab_m[name==names]
    print Mstar
    Mstar = 10**Mstar / (1e11)

    pl.xlabel('dark halo inner slope',fontsize=30)
    pl.ylabel('stellar mass \n relative to a Milky-Way-like IMF',fontsize=30)
    #pl.axvline(1,color='k')
    #pl.axhline(0,color='k')
    #pl.text(0.81,1,'NFW',rotation=90,fontsize=30)
    #pl.text(0.5,-0.19,'MILKY \n WAY',fontsize=30,horizontalalignment='center')

    pl.savefig('/data/ljo31/public_html/Lens/dark_and_light/compare_profiles_new/contour_plot__%s.pdf'%name)
    pl.show()



lz = np.load('/data/ljo31/Lens/LensParams/LensRedshiftsUpdated.npy')[()]
sz = np.load('/data/ljo31/Lens/LensParams/SourceRedshiftsUpdated.npy')[()]

chab = np.load('/data/ljo31b/EELs/inference/new/huge/chabrier_masses_212.npy')
salp = np.load('/data/ljo31b/EELs/inference/new/huge/salpeter_masses_212.npy')
chab_m,salp_m = chab[0], salp[0]

names = np.array(['J0837','J0901','J0913','J1125','J1144','J1218','J1323','J1347','J1446','J1605','J1619','J2228'])

for name in ['J0837','J0901','J0913','J1144']:
    plot(name)
