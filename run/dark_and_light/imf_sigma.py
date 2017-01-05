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

    ml_dpl = np.median(ml_dpl)
    dml = np.percentile(dic_DPL['stellar mass'][:,0].ravel()*10,84)-ml_dpl

    Mstar = salp_m[name==names]
    Mstar = 10**Mstar / (1e11)
    dMstar = Mstar * dsalp

    DML = ml_dpl/Mstar

    sigma = sigmas[names==name]
    dsigma = dsigmas[names==name]

    dDML = (dml/Mstar)**2. + (dMstar * ml_dpl/Mstar)**2.
    dDML = dDML**0.5

    pl.errorbar(np.log10(sigma),np.log10(DML),xerr=dsigma/sigma/np.log(10),yerr=dDML/DML/np.log(10.),color='k')
    
    




lz = np.load('/data/ljo31/Lens/LensParams/LensRedshiftsUpdated.npy')[()]
sz = np.load('/data/ljo31/Lens/LensParams/SourceRedshiftsUpdated.npy')[()]

chab = np.load('/data/ljo31b/EELs/inference/new/huge/chabrier_masses_212.npy')
salp = np.load('/data/ljo31b/EELs/inference/new/huge/salpeter_masses_212.npy')
chab_m,salp_m = chab[0], salp[0]
dchab = np.mean((chab[1],chab[2]))
dsalp = np.mean((salp[1],salp[2]))


names = np.array(['J0837','J0901','J0913','J1125','J1144','J1218','J1323','J1347','J1446','J1605','J1619','J2228'])
colours = np.array(['SteelBlue','Crimson','SeaGreen','DarkOrange'])

# DATA
fp = np.load('/data/ljo31b/EELs/esi/kinematics/inference/results_0.30_source_indous_vdfit_jul2016_J2228.npy')
l,m,u = fp
d = np.mean((l,u),axis=0)
dvl,dvs,dsigmal,dsigmas = d.T
vl,vs,sigmal,sigmas = m.T
dsigmas = sigmas*0.05


ii=0
for name in ['J0837','J0901','J0913','J1144']:
    plot(name)
    ii+=1

ls = np.linspace(10**2.15,10**2.50,50)
tt = 1.31 * np.log10(ls) - 3.14
tt_l = (1.31-0.16) * np.log10(ls) - 3.15
tt_u = (1.31+0.16) * np.log10(ls) - 3.13
pl.plot(np.log10(ls),tt,c='k',ls='--')
pl.fill_between(np.log10(ls),tt_l,tt_u,color='LightGray',alpha=0.5)
pl.xlabel(r'$\log \sigma$')
pl.ylabel(r'$\log \alpha_{salp}$')
pl.show()
