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

name = 'J0837'

dir = '/data/ljo31b/EELs/galsub/emceeruns/'
file = glob.glob(dir+name+'_parametric_*')
file.sort()
f = file[-1]
while 1:
    if 'DPL' in f:
        file = file[:-1]
        f = file[-1]
    else:
        break
print f
result_PL = np.load(f)

file = glob.glob(dir+name+'_parametric_DPL*')
file.sort()
f = file[-1]
print f
result_DPL = np.load(f)

lp_PL,trace_PL,dic_PL,_ = result_PL
lp_DPL,trace_DPL,dic_DPL,_ = result_DPL
a1_PL,a3_PL = np.unravel_index(lp_PL[:,0].argmax(),lp_PL[:,0].shape)
a1_DPL,a3_DPL = np.unravel_index(lp_DPL[:,0].argmax(),lp_DPL[:,0].shape)

## HACK
if name == 'J1446':
    nn = 0
    while nn<10:
        row = lp_PL[-1,0,:]
        lp_PL = lp_PL[:,:,row>row.min()]
        trace_PL = trace_PL[:,:,row>row.min(),:]
        for key in dic_PL.keys():
            dic_PL[key] = dic_PL[key][:,:,row>row.min()]
        nn+=1
    nn = 0
    while nn<15:
        row = lp_PL[0,0,:]
        lp_PL = lp_PL[:,:,row>row.min()]
        trace_PL = trace_PL[:,:,row>row.min(),:]
        for key in dic_PL.keys():
            dic_PL[key] = dic_PL[key][:,:,row>row.min()]
        nn+=1
    nn = 0
    while nn<10:
        row = lp_PL[250,0,:]
        lp_PL = lp_PL[:,:,row>row.min()]
        trace_PL = trace_PL[:,:,row>row.min(),:]
        for key in dic_PL.keys():
            dic_PL[key] = dic_PL[key][:,:,row>row.min()]
        nn+=1
    

    pl.figure()
    pl.plot(lp_PL[:,0])
    pl.show()
                                

## construct power-law sigma and mass
zl,zs = lz[name][0], sz[name][0]
sig_crit = lenslib.sig_crit(zl,zs) # solar masses per Mpc^2
sig_crit /= (1e3)**2. # solar masses per kpc^2
scale = astCalc.da(zl)*np.pi/180./3600 * 1e3
rein = np.median(dic_PL['Lens 1 b'])*0.05*scale # rein in kpc
eta = np.median(dic_PL['Lens 1 eta'])

r = np.logspace(-7,4,4500)
lr = r[:-100]
cumsig_PL,sigma_PL = PowerLaw(r,sig_crit,rein,eta,projected=True)

### gNFW sigma and mass
rs, rein = np.median(dic_DPL['Lens 1 rs'])*0.05*scale,np.median(dic_DPL['Lens 1 b'])*0.05*scale 
gamma = np.median(dic_DPL['Lens 1 gamma'])
cumsig_DPL,sigma_DPL = gNFW(r,sig_crit,rein,gamma,rs,projected=True)

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

Mstar_PL = np.median(dic_PL['stellar mass'])
Mstar_DPL = np.median(dic_DPL['stellar mass'])
_,sb_M,cum_M = GetStellarMass.GetStellarMass(oresult,1.,r,scale,name,deprojected=False)
sb_M_PL, sb_M_DPL = sb_M*Mstar_PL, sb_M*Mstar_DPL
cum_M_PL, cum_M_DPL = cum_M*Mstar_PL, cum_M*Mstar_DPL

## finally, add original inference for total mass
rein, eta = cat[name]['Lens 1 b']*0.05*scale, cat[name]['Lens 1 eta']
cumsig_TPL,sigma_TPL = PowerLaw(r,sig_crit,rein,eta,projected=True)

### get projected DM fractions
f_DM_PL =  cumsig_PL[:-100]/(cumsig_PL[:-100]+cum_M_PL)
f_DM_DPL = cumsig_DPL[:-100]/(cumsig_DPL[:-100]+cum_M_DPL)

reff = r_eff[names==name]

mod_F_PL = splrep(lr[100:],f_DM_PL[100:])
print splev(scale,mod_F_PL)
f_DM_PL = splev(0.5*reff,mod_F_PL)
mod_F_DPL = splrep(lr[100:],f_DM_DPL[100:])
f_DM_DPL = splev(0.5*reff,mod_F_DPL)

print f_DM_PL, f_DM_DPL



## non-cum, projected
pl.figure(figsize=(18,7))
pl.subplot(121)
pl.plot(lr,sb_M_DPL,color='SteelBlue',ls='--')#,label='LM gNFW')
pl.plot(lr,sb_M_PL,color='Crimson',ls='--')#,label='LM PL')
pl.plot(lr,sigma_DPL[:-100],color='SteelBlue',ls=':')#,label='DM gNFW')
pl.plot(lr,sigma_PL[:-100],color='Crimson',ls=':')#,label='DM PL')
pl.yscale('log')
print ((sb_M_PL-sb_M_DPL)/sb_M_PL)[-1]
pl.axis([0,20,5e6,1e12])
pl.ylabel('projected mass surface density')
pl.axvline(rein,color='k',ls='--') # Einstein radius
names = np.array(names)
pl.axvline(r_eff[names==name],color='k',ls='--') # effective radius
pl.xlabel('R / kpc')
pl.figtext(0.32,0.425,'$f_{DM}^{PL}(<0.5R_e) = $'+'%.2f'%f_DM_PL)
pl.figtext(0.31,0.35,'$f_{DM}^{gNFW}(<0.5R_e) = $'+'%.2f'%f_DM_DPL)


# cum projected
pl.subplot(122)
pl.plot(lr,cum_M_DPL,color='SteelBlue',ls='--')#,label='LM gNFW')
pl.plot(lr,cum_M_PL,color='Crimson',ls='--')#,label='LM PL')
pl.plot(lr,cumsig_DPL[:-100],color='SteelBlue',ls=':')#,label='DM gNFW')
pl.plot(lr,cumsig_PL[:-100],color='Crimson',ls=':')#,label='DM PL')
pl.yscale('log')
pl.axis([0,10,1e8,6e11])
pl.ylabel('cumulative projected mass')
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
pl.legend(loc='lower right',fontsize=20)
pl.suptitle(name)
pl.show()

#pl.figure()
#pl.plot(lr,cumsig_DPL[:-100]/(cumsig_DPL[:-100]+cum_M_DPL),color='SteelBlue')
#pl.plot(lr,cumsig_PL[:-100]/(cumsig_PL[:-100]+cum_M_PL),color='Crimson')
#pl.xlim([0,20])
#pl.show()

