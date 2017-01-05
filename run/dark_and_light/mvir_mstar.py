import numpy as np, pylab as pl, pyfits as py
from MWApython.pylens import MassModels
import indexTricks as iT
import glob
from CalculateConc import *
from jeans.makemodel import buildhalo,virialRadius
from astLib import astCalc

names = ['J0837','J0901','J0913','J1125','J1144']
lz = np.load('/data/ljo31/Lens/LensParams/LensRedshiftsUpdated.npy')[()]
sz = np.load('/data/ljo31/Lens/LensParams/SourceRedshiftsUpdated.npy')[()]

dir = '/data/ljo31b/EELs/galsub/emceeruns/'

pl.figure()
# add Behroozi's relations
mstars = np.logspace(9.5,11.8,150)
mvirs = mstars*0.
for ii in range(mstars.size):
    mvirs[ii] = buildhalo(mstars[ii])

pl.plot(np.log10(mstars),mvirs,c='DarkGray',label='Behroozi+2010')
pl.fill_between(np.log10(mstars),mvirs-0.5,mvirs+0.5,color='LightGray',alpha=0.5)
pl.xlabel('$\log M_{\star}/M_{\odot}$')
pl.ylabel('$\log M_{vir}/M_{\odot}$')


for name in names:
    file = glob.glob(dir+name+'_parametric_*')
    file.sort()
    f = file[-1]
    while 1:
        if 'DPL' in f:
            file = file[:-1]
            f = file[-1]
        else:
            break
    result_PL = np.load(f)

    file = glob.glob(dir+name+'_parametric_DPL*')
    file.sort()
    f = file[-1]
    if name == 'J0901':
        f = file[-2]
    result_DPL = np.load(f)

    lp_PL,trace_PL,dic_PL,_ = result_PL
    lp_DPL,trace_DPL,dic_DPL,_ = result_DPL
    
    a1_PL,a3_PL = np.unravel_index(lp_PL[:,0].argmax(),lp_PL[:,0].shape)
    a1_DPL,a3_DPL = np.unravel_index(lp_DPL[:,0].argmax(),lp_DPL[:,0].shape)

    q_pl = np.median(dic_PL['Lens 1 q'][:,0])
    mstar_pl = np.median(dic_PL['stellar mass'][:,0])*1e12 
    b_pl = np.median(dic_PL['Lens 1 b'][:,0])
    eta_pl = np.median(dic_PL['Lens 1 eta'][:,0]) 
    
    q_dpl = np.median(dic_DPL['Lens 1 q'][:,0]) 
    mstar_dpl = np.median(dic_DPL['stellar mass'][:,0])*1e12 
    b_dpl = np.median(dic_DPL['Lens 1 b'][:,0]) 
    gamma_dpl = np.median(dic_DPL['Lens 1 gamma'][:,0]) 
    rs_dpl = np.median(dic_DPL['Lens 1 rs'][:,0])

    '''q_pl = dic_PL['Lens 1 q'][a1_PL,0,a3_PL] 
    mstar_pl = dic_PL['stellar mass'][a1_PL,0,a3_PL]*1e12 
    b_pl = dic_PL['Lens 1 b'][a1_PL,0,a3_PL] 
    eta_pl = dic_PL['Lens 1 eta'][a1_PL,0,a3_PL] 
    
    q_dpl = dic_DPL['Lens 1 q'][a1_DPL,0,a3_DPL] 
    mstar_dpl = dic_DPL['stellar mass'][a1_DPL,0,a3_DPL]*1e12 
    b_dpl = dic_DPL['Lens 1 b'][a1_DPL,0,a3_DPL] 
    gamma_dpl = dic_DPL['Lens 1 gamma'][a1_DPL,0,a3_DPL] 
    rs_dpl = dic_DPL['Lens 1 rs'][a1_DPL,0,a3_DPL] '''

    zl,zs = lz[name][0],sz[name][0]
    scale = astCalc.da(zl)*np.pi/180./3600 * 1e3

    # get virial masses
    mvir_dpl, c_dpl, rvir_dpl = CalcLensDPL(zl,zs,gamma_dpl,rs_dpl*0.05*scale,b_dpl*0.05*scale)
    mvir_pl, rvir_pl = CalcLensPL(zl,zs,eta_pl,b_dpl*0.05*scale)
    print name, '%.2f'%rvir_pl, '%.2f'%rvir_dpl,'%.2f'%np.log10(mvir_pl), '%.2f'%np.log10(mvir_dpl)

    # figure
    if name == 'J0837':
        pl.scatter(np.log10(mstar_dpl),np.log10(mvir_dpl),c='SteelBlue',s=40,label='gNFW')
        pl.scatter(np.log10(mstar_pl),np.log10(mvir_pl),c='Crimson',s=40,label='power law')
    else:
        pl.scatter(np.log10(mstar_dpl),np.log10(mvir_dpl),c='SteelBlue',s=40)
        pl.scatter(np.log10(mstar_pl),np.log10(mvir_pl),c='Crimson',s=40)
    
pl.legend(loc='upper left',fontsize=20,scatterpoints=1)
pl.show()
