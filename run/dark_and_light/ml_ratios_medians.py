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
        f = file[-3]
    result_DPL = np.load(f)

    lp_PL,trace_PL,dic_PL,_ = result_PL
    lp_DPL,trace_DPL,dic_DPL,_ = result_DPL
    
    q_pl = np.median(dic_PL['Lens 1 q'][:,0])
    mstar_pl = np.median(dic_PL['stellar mass'][:,0])*1e12 
    b_pl = np.median(dic_PL['Lens 1 b'][:,0])
    eta_pl = np.median(dic_PL['Lens 1 eta'][:,0]) 
    
    q_dpl = np.median(dic_DPL['Lens 1 q'][:,0]) 
    mstar_dpl = np.median(dic_DPL['stellar mass'][:,0])*1e12 
    b_dpl = np.median(dic_DPL['Lens 1 b'][:,0]) 
    gamma_dpl = np.median(dic_DPL['Lens 1 gamma'][:,0]) 
    rs_dpl = np.median(dic_DPL['Lens 1 rs'][:,0])

    Lv_star = Lv[names==name][0]
    ML_PL = mstar_pl/Lv_star
    ML_DPL = mstar_dpl/Lv_star

    #pl.scatter(np.log10(mstar_pl),ML_PL,c='Crimson',s=40)
    #pl.scatter(np.log10(mstar_dpl),ML_DPL,c='SteelBlue',s=40)

    pl.scatter(np.log10(Lv_star),ML_PL,c='Crimson',s=40)
    pl.scatter(np.log10(Lv_star),ML_DPL,c='SteelBlue',s=40)
    #print name, '%.2f'%ML_PL,'%.2f'%ML_DPL
    #print name, '%.2f'%np.log10(mstar_pl), '%.2f'%np.log10(mstar_dpl)

    zl,zs = lz[name][0],sz[name][0]
    scale = astCalc.da(zl)*np.pi/180./3600 * 1e3
    #print name, '%.2f,%.2f,%.2f,%.2f,%.2f'%(eta_pl, gamma_dpl, rs_dpl*0.05*scale,q_pl, q_dpl)
    print name, '%.2f,%.2f'%(b_pl*0.05*scale, b_dpl*scale*0.05)


pl.legend(loc='upper left',scatterpoints=1,fontsize=15)
pl.xlim(9.75,11.5)
#pl.show()
    
