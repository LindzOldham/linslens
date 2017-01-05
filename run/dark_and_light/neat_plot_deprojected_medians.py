from linslens.Profiles import *
from astLib import astCalc
import glob, numpy as np, pylab as pl
import lenslib
from jeans.makemodel import deproject
import GetStellarMass
import cPickle
from scipy.interpolate import splev,splrep,splint

names = ['J0837','J0901','J0913','J1125','J1144','J1218','J1323','J1347','J1446','J1605','J1606','J1619','J2228']
Mvir, dMvir, r_eff, dr_eff = np.load('/data/ljo31/Lens/LensParams/Mvir.npy')
lz = np.load('/data/ljo31/Lens/LensParams/LensRedshiftsUpdated.npy')[()]
sz = np.load('/data/ljo31/Lens/LensParams/SourceRedshiftsUpdated.npy')[()]

name = 'J0901'

dir = '/data/ljo31b/EELs/galsub/emceeruns/'
file = glob.glob(dir+name+'_parametric_DPL*')
file.sort()
f = file[-1]
print f
result_DPL = np.load(f)

lp_DPL,trace_DPL,dic_DPL,_ = result_DPL
a1_DPL,a3_DPL = np.unravel_index(lp_DPL[:,0].argmax(),lp_DPL[:,0].shape)
                 

## construct power-law sigma and mass
zl,zs = lz[name][0], sz[name][0]
sig_crit = lenslib.sig_crit(zl,zs) # solar masses per Mpc^2
sig_crit /= (1e3)**2. # solar masses per kpc^2
scale = astCalc.da(zl)*np.pi/180./3600 * 1e3

r = np.logspace(-7,4,4500)
lr = r[:-100]

### gNFW sigma and mass
rs, rein = np.median(dic_DPL['Lens 1 rs'])*0.05*scale,np.median(dic_DPL['Lens 1 b'])*0.05*scale 
gamma = np.median(dic_DPL['Lens 1 gamma'])
cumsig_DPL,sigma_DPL = gNFW(r,sig_crit,rein,gamma,rs,projected=False)

### also old result -- Einstein radius
cat = np.load('/data/ljo31/Lens/LensParams/Structure_lensgals_2src.npy')[0]

# add in LM!
dir = '/data/ljo31/Lens/LensModels/twoband/'
try:
    oresult = np.load(dir+name+'_212')
except:
    if name == 'J1347':
        oresult = np.load(dir+name+'_112')
    else:
        oresult = np.load(dir+name+'_211')

Mstar_DPL = np.median(dic_DPL['stellar mass'])
_,sb_M,cum_M = GetStellarMass.GetStellarMass(oresult,1.,r,scale,name,deprojected=True)
sb_M_DPL = sb_M*Mstar_DPL
cum_M_DPL = cum_M*Mstar_DPL

## finally, add original inference for total mass
rein, eta = cat[name]['Lens 1 b']*0.05*scale, cat[name]['Lens 1 eta']
reff = r_eff[names==name]

## non-cum, projected
pl.figure()
pl.plot(lr,sb_M_DPL,color='Crimson',label='stellar mass')#,label='LM gNFW')
pl.plot(lr,sigma_DPL[:-100],color='SteelBlue',label='dark matter')#,label='DM gNFW')
pl.yscale('log')
pl.xscale('log')
pl.axis([0.14,140,6000,9e9])
pl.ylabel(r'$\rho$ / M$_{\star}$kpc$^{-3}$')
pl.axvline(rein,color='k',ls='--') # Einstein radius
names = np.array(names)
pl.axvline(r_eff[names==name],color='k',ls='--') # effective radius
pl.xlabel('r / kpc')
pl.figtext(0.461,0.45,'Einstein radius',rotation=90)
pl.figtext(0.576,0.45,'effective radius',rotation=90)

# uncertaintiezzz
rein, gamma, rs, Mstar = dic_DPL['Lens 1 b'][:,0].ravel()[::5000]*0.05*scale, dic_DPL['Lens 1 gamma'][:,0].ravel()[::5000], dic_DPL['Lens 1 rs'][:,0].ravel()[::5000]*0.05*scale, dic_DPL['stellar mass'][:,0].ravel()[::5000]
sigma = np.zeros((rein.size, r.size))
for ii in range(rein.size):
    REIN, GAMMA, RS, MSTAR = rein[ii], gamma[ii], rs[ii], Mstar[ii]
    CUMSIG,SIGMA = gNFW(r,sig_crit,REIN,GAMMA,RS,projected=False)
    sigma[ii] = SIGMA
    print ii

l, u, m = r*0.,r*0.,r*0.
for ii in range(r.size):
    s = sigma[:,ii]
    m[ii] = np.median(s)
    l[ii],u[ii] = np.percentile(s,16),np.percentile(s,84)

mm,ll,uu = np.median(Mstar), np.percentile(Mstar,1), np.percentile(Mstar,99)
#pl.plot(r,m,color='SteelBlue')
pl.fill_between(r,l,u,color='LightBlue',alpha=0.8)
pl.fill_between(lr,ll*sb_M, uu*sb_M,color='LightPink',alpha=0.5)
pl.legend(loc='upper right',fontsize=25)

pl.show()
