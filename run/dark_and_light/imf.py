import numpy as np, pylab as pl, pyfits as py, glob

chab = np.load('/data/ljo31b/EELs/inference/new/huge/chabrier_masses_212.npy')
salp = np.load('/data/ljo31b/EELs/inference/new/huge/salpeter_masses_212.npy')

chab_m,salp_m = chab[0], salp[0]

names = np.array(['J0837','J0901','J0913','J1125','J1144'])#,'J1218','J1323','J1347','J1446','J1605','J1606','J1619','J2228'])

dir = '/data/ljo31b/EELs/galsub/emceeruns/'
dir3 = '/data/ljo31b/EELs/esi/kinematics/inference/vdfit/NEW/'

for name in names:
    n = np.where(names==name)
    m_chab, m_salp = chab_m[n], salp_m[n]

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

    m_PL, m_DPL = np.median(dic_PL['stellar mass'][:,0])*1e12, np.median(dic_DPL['stellar mass'][:,0])*1e12
    m_PL, m_DPL = np.log10(m_PL), np.log10(m_DPL)

    
    alpha_PL_chab = 10**(m_PL-m_chab)
    alpha_DPL_chab = 10**(m_DPL-m_chab)
    alpha_PL_salp = 10**(m_PL-m_salp)
    alpha_DPL_salp = 10**(m_DPL-m_salp)
    #print name, alpha_PL_chab, alpha_DPL_chab

    result_K = np.load(dir3+name+'_1.00_lens_esi_indous_vdfit_LENS')
    lp_K, trace_K, dic_K, _ = result_K
    s_lens = np.median(dic_K['lens dispersion'])
    pl.scatter(s_lens,alpha_PL_chab,c='Crimson',s=40)
    pl.scatter(s_lens,alpha_DPL_chab,c='SteelBlue',s=40)
    print name, '%.f'%s_lens
    
x=np.linspace(210,270,10)
y=np.ones(10)
pl.plot(x,y,label='Chabrier',color='LightGray')
pl.plot(x,y*1.7,label='Salpeter',color='k')
pl.scatter(0,0,color='Crimson',label='power law',s=40)
pl.scatter(0,0,color='SteelBlue',label='gNFW',s=40)
pl.legend(loc='upper left',scatterpoints=1,fontsize=20)
pl.xlim(215,270)
pl.xlabel('$\sigma$/kms$^{-1}$')
pl.ylabel('$M_{\star,PL}/M_{\star,chab}$')
pl.show()
