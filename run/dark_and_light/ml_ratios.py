import numpy as np, pylab as pl, pyfits as py
from MWApython.pylens import MassModels
import indexTricks as iT
import glob
from CalculateConc import *
from jeans.makemodel import buildhalo,virialRadius
from astLib import astCalc

names = np.array(['J0837','J0901','J0913','J1125','J1144'])
lz = np.load('/data/ljo31/Lens/LensParams/LensRedshiftsUpdated.npy')[()]
sz = np.load('/data/ljo31/Lens/LensParams/SourceRedshiftsUpdated.npy')[()]

dir = '/data/ljo31b/EELs/galsub/emceeruns/'

pl.figure()

#pl.xlabel('$\log M_{\star}/M_{\odot}$')
pl.xlabel('$\log L_{V}/L_{V,\odot}$')

pl.ylabel('$M_{\star}/L_{V}$')
pl.scatter(0,0,c='Crimson',label='power law')
pl.scatter(0,0,c='SteelBlue',label='gNFW')
x=np.linspace(8,12,10)
y=np.ones(10)
pl.plot(x,y*2.1,color='DarkGray',label=r'Chabrier ($6 \times 10^9$ Gy)')
pl.plot(x,y*3.8,color='k',label=r'Salpeter ($6 \times 10^9$ Gyr)')


Lb = np.load('/data/ljo31/Lens/LensParams/B_redshift0_model_ABSMAG_LENSES.npy')
Lv = np.load('/data/ljo31/Lens/LensParams/V_redshift0_model_ABSMAG_LENSES.npy')

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

    q_pl = dic_PL['Lens 1 q'][a1_PL,0,a3_PL] 
    mstar_pl = dic_PL['stellar mass'][a1_PL,0,a3_PL]*1e12 
    b_pl = dic_PL['Lens 1 b'][a1_PL,0,a3_PL] 
    eta_pl = dic_PL['Lens 1 eta'][a1_PL,0,a3_PL] 
    
    q_dpl = dic_DPL['Lens 1 q'][a1_DPL,0,a3_DPL] 
    mstar_dpl = dic_DPL['stellar mass'][a1_DPL,0,a3_DPL]*1e12 
    b_dpl = dic_DPL['Lens 1 b'][a1_DPL,0,a3_DPL] 
    gamma_dpl = dic_DPL['Lens 1 gamma'][a1_DPL,0,a3_DPL] 
    rs_dpl = dic_DPL['Lens 1 rs'][a1_DPL,0,a3_DPL] 

    Lv_star = Lv[names==name][0]
    ML_PL = mstar_pl/Lv_star
    ML_DPL = mstar_dpl/Lv_star

    #pl.scatter(np.log10(mstar_pl),ML_PL,c='Crimson',s=40)
    #pl.scatter(np.log10(mstar_dpl),ML_DPL,c='SteelBlue',s=40)

    pl.scatter(np.log10(Lv_star),ML_PL,c='Crimson',s=40)
    pl.scatter(np.log10(Lv_star),ML_DPL,c='SteelBlue',s=40)
    print name, '%.2f'%ML_PL,'%.2f'%ML_DPL

pl.legend(loc='upper left',scatterpoints=1,fontsize=15)
pl.xlim(9.75,11.5)
pl.show()
    
