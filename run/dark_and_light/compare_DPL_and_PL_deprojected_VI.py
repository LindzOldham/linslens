from linslens.Profiles import *
from astLib import astCalc
import glob, numpy as np, pylab as pl
import lenslib
from jeans.makemodel import deproject
import GetStellarMass
import cPickle
from scipy.interpolate import splev,splrep,splint
from linslens import PROFILES
import sys

name = sys.argv[1]

names = ['J0837','J0901','J0913','J1125','J1144','J1218','J1323','J1347','J1446','J1605','J1606','J1619','J2228']
Mvir, dMvir, r_eff, dr_eff = np.load('/data/ljo31/Lens/LensParams/Mvir.npy')
lz = np.load('/data/ljo31/Lens/LensParams/LensRedshiftsUpdated.npy')[()]
sz = np.load('/data/ljo31/Lens/LensParams/SourceRedshiftsUpdated.npy')[()]

#name = 'J1144'

dir = '/data/ljo31b/EELs/galsub/emceeruns/'
file = glob.glob(dir+name+'_parametric_VI_*')
file.sort()
f = file[-1]
result_PL = np.load(f)

file = glob.glob(dir+name+'_parametric_DPL__VI_*')
file.sort()
f = file[-1]
print f 
result_DPL = np.load(f)


lp_PL,trace_PL,dic_PL,_ = result_PL
lp_DPL,trace_DPL,dic_DPL,_ = result_DPL
                 
### POWER LAW ####
## construct power-law sigma and mass
zl,zs = lz[name][0], sz[name][0]
sig_crit = lenslib.sig_crit(zl,zs) # solar masses per Mpc^2
sig_crit /= (1e3)**2. # solar masses per kpc^2
scale = astCalc.da(zl)*np.pi/180./3600 * 1e3
rein = np.median(dic_PL['Lens 1 b'][:,0])*0.05*scale # rein in kpc
eta = np.median(dic_PL['Lens 1 eta'][:,0])
q = np.median(dic_PL['Lens 1 q'][:,0])

# make mass objects
r = np.logspace(-7,4,4500)
lr = r[:-100]
m_PL = PROFILES.PowerLaw(x=0.,y=0.,eta=eta,pa=0.,q=q,b=rein,zl=zl,zs=zs)
mass_PL = m_PL.mass(r)
rho_PL = m_PL.rho(r)

#### GNFW ####
### gNFW sigma and mass
rs, rein = np.median(dic_DPL['Lens 1 rs'][:,0])*0.05*scale,np.median(dic_DPL['Lens 1 b'][:,0])*0.05*scale 
gamma, q = np.median(dic_DPL['Lens 1 gamma'][:,0]), np.median(dic_DPL['Lens 1 q'][:,0])
m_DPL = PROFILES.gNFW_TC(x=0.,y=0.,eta=gamma,rs=rs,pa=0.,q=q,b=rein,zl=zl,zs=zs)
mass_DPL = m_DPL.mass(r)
rho_DPL = m_DPL.rho(r)


### also old result -- Einstein radius
cat = np.load('/data/ljo31/Lens/LensParams/Structure_lensgals_2src.npy')[0]
rein = cat[name]['Lens 1 b']*0.05*scale

# add in LM!
dir = '/data/ljo31/Lens/LensModels/twoband/'
try:
    oresult = np.load(dir+name+'_212')
except:
    if name == 'J1347':
        oresult = np.load(dir+name+'_112')
    else:
        oresult = np.load(dir+name+'_211')

Mstar_PL = np.median(dic_PL['stellar mass'][:,0])
Mstar_DPL = np.median(dic_DPL['stellar mass'][:,0])
_,sb_M,cum_M = GetStellarMass.GetStellarMass(oresult,1.,r,scale,name,deprojected=True)
sb_M_PL, sb_M_DPL = sb_M*Mstar_PL, sb_M*Mstar_DPL
cum_M_PL, cum_M_DPL = cum_M*Mstar_PL, cum_M*Mstar_DPL


## finally, add original inference for total mass
rein, eta = cat[name]['Lens 1 b']*0.05*scale, cat[name]['Lens 1 eta']
q = cat[name]['Lens 1 q']
m_TPL = PROFILES.PowerLaw(x=0.,y=0.,eta=eta,pa=0.,q=q,b=rein,zl=zl,zs=zs)
mass_TPL = m_TPL.mass(r)
rho_TPL = m_TPL.rho(r)




## non-cum, projected
pl.figure(figsize=(18,7))
pl.subplot(121)
pl.plot(lr,sb_M_DPL,color='SteelBlue',ls='--')#,label='LM gNFW')
pl.plot(lr,sb_M_PL,color='Crimson',ls='--')#,label='LM PL')
pl.plot(lr,rho_DPL[:-100],color='SteelBlue',ls=':')#,label='DM gNFW')
pl.plot(lr,rho_PL[:-100],color='Crimson',ls=':')#,label='DM PL')
pl.plot(lr,sb_M_DPL+rho_DPL[:-100],color='SteelBlue',ls='-')
pl.plot(lr,sb_M_PL+rho_PL[:-100],color='Crimson',ls='-')
pl.plot(lr,rho_TPL[:-100],color='k',ls='-')

pl.yscale('log')
print ((sb_M_PL-sb_M_DPL)/sb_M_PL)[-1]
pl.axis([0,20,5e6,1e12])
pl.ylabel('3D mass density')
pl.axvline(rein,color='k',ls='--') # Einstein radius
names = np.array(names)
pl.axvline(r_eff[names==name],color='k',ls='--') # effective radius
pl.xlabel('r / kpc')
#pl.figtext(0.32,0.425,'$f_{DM}^{PL}(<0.5R_e) = $'+'%.2f'%f_DM_PL)
#pl.figtext(0.31,0.35,'$f_{DM}^{gNFW}(<0.5R_e) = $'+'%.2f'%f_DM_DPL)


# cum projected
pl.subplot(122)
pl.plot(lr,cum_M_DPL,color='SteelBlue',ls='--')#,label='LM gNFW')
pl.plot(lr,cum_M_PL,color='Crimson',ls='--')#,label='LM PL')
pl.plot(lr,mass_DPL,color='SteelBlue',ls=':')#,label='DM gNFW')
pl.plot(lr,mass_PL[:-100],color='Crimson',ls=':')#,label='DM PL')
pl.plot(lr,cum_M_DPL+mass_DPL,color='SteelBlue',ls='-')
pl.plot(lr,cum_M_PL+mass_PL[:-100],color='Crimson',ls='-')
pl.plot(lr,mass_TPL[:-100],color='k',ls='-')
pl.yscale('log')
pl.axis([0,10,1e8,6e11])
pl.ylabel('cumulative deprojected mass')
pl.xlabel('R')
pl.axvline(rein,color='k',ls='--') # Einstein radius
names = np.array(names)
pl.axvline(r_eff[names==name],color='k',ls='--') # effective radius
##
pl.plot(0,0,color='LightGray',ls='-',label='DM+LM')
pl.plot(0,0,color='LightGray',ls=':',label='DM')
pl.plot(0,0,color='LightGray',ls='--',label='LM')
pl.plot(0,0,color='SteelBlue',ls='-',label='gNFW')
pl.plot(0,0,color='Crimson',ls='-',label='PL')
pl.plot(0,0,color='k',ls='-',label='TPL')
pl.legend(loc='lower right',fontsize=20)
pl.suptitle(name)
pl.savefig('/data/ljo31/public_html/Lens/dark_and_light/compare_profiles_new/deprojected_VI_%s.pdf'%name)
pl.show()

#pl.figure()
#pl.plot(lr,cumsig_DPL[:-100]/(cumsig_DPL[:-100]+cum_M_DPL),color='SteelBlue')
#pl.plot(lr,cumsig_PL[:-100]/(cumsig_PL[:-100]+cum_M_PL),color='Crimson')
#pl.xlim([0,20])
#pl.show()

apl,adpl = np.min(dic_PL['stellar mass'][:,0].ravel()*10), np.min(dic_DPL['stellar mass'][:,0].ravel()*10)
bpl,bdpl = np.max(dic_PL['stellar mass'][:,0].ravel()*10), np.max(dic_DPL['stellar mass'][:,0].ravel()*10)


pl.figure()
pl.hist(dic_PL['stellar mass'][:,0].ravel()*10,bins=np.arange(apl-0.05,bpl+0.05,(bpl-apl+0.1)/30.),alpha=0.5,label='power law',histtype='stepfilled',normed=True)
pl.hist(dic_DPL['stellar mass'][:,0].ravel()*10,bins=np.arange(adpl-0.05,bdpl+0.05,(bdpl-adpl+0.1)/30.),alpha=0.5,label='gNFW',histtype='stepfilled',normed=True)
pl.legend(loc='upper left',fontsize=25)
pl.xlabel('M$_{\star}$ ( 10$^{11}$ M$_{\odot}$)')
pl.ylabel('probability density')
pl.savefig('/data/ljo31/public_html/Lens/dark_and_light/compare_profiles_new/stellarmass_VI_%s.pdf'%name)
pl.show()

pl.figure(figsize=(15,7))
pl.subplot(121)
pl.plot(lp_PL[:,0])
pl.subplot(122)
pl.plot(lp_DPL[:,0])
pl.show()

apl,adpl = np.min(dic_PL['Lens 1 eta'][:,0].ravel()+1.), np.min(dic_DPL['Lens 1 gamma'][:,0].ravel())
bpl,bdpl = np.max(dic_PL['Lens 1 eta'][:,0].ravel()+1.), np.max(dic_DPL['Lens 1 gamma'][:,0].ravel())


pl.figure()
pl.hist(dic_PL['Lens 1 eta'][:,0].ravel()+1.,bins=np.arange(apl-0.05,bpl+0.05,(bpl-apl+0.1)/30.),alpha=0.5,label='power law',histtype='stepfilled',normed=True)

pl.hist(dic_DPL['Lens 1 gamma'][:,0].ravel(),bins=np.arange(adpl-0.05,bdpl+0.05,(bdpl-adpl+0.1)/30.),alpha=0.5,label='gNFW',histtype='stepfilled',normed=True)
pl.legend(loc='upper left',fontsize=25)
pl.xlabel('3D density slope')
pl.ylabel('probability density')
pl.savefig('/data/ljo31/public_html/Lens/dark_and_light/compare_profiles_new/slope_VI_%s.pdf'%name)
pl.show()


