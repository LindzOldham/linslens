from linslens.Profiles import *
from astLib import astCalc
import glob, numpy as np, pylab as pl
import lenslib
from jeans.makemodel import deproject
import GetStellarMass
import cPickle

names = ['J0837','J0901','J0913','J1125','J1144','J1218','J1323','J1347','J1446','J1605','J1606','J1619','J2228']
Mvir, dMvir, r_eff, dr_eff = np.load('/data/ljo31/Lens/LensParams/Mvir.npy')
lz = np.load('/data/ljo31/Lens/LensParams/LensRedshiftsUpdated.npy')[()]
sz = np.load('/data/ljo31/Lens/LensParams/SourceRedshiftsUpdated.npy')[()]

name = 'J0913'

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
result_param = np.load(f)

file = glob.glob(dir+name+'_pixellated_*')
file.sort()
f = file[-1]
print f
result_pix = np.load(f)

lp_param,trace_param,dic_param,_ = result_param
lp_pix,trace_pix,dic_pix,_ = result_pix
a1_param,a3_param = np.unravel_index(lp_param[:,0].argmax(),lp_param[:,0].shape)
a1_pix,a3_pix = np.unravel_index(lp_pix[:,0].argmax(),lp_pix[:,0].shape)

## construct power-law sigma and mass
zl,zs = lz[name][0], sz[name][0]
sig_crit = lenslib.sig_crit(zl,zs) # solar masses per Mpc^2
sig_crit /= (1e3)**2. # solar masses per kpc^2
scale = astCalc.da(zl)*np.pi/180./3600 * 1e3
r = np.logspace(-7,3,3500)
lr = r[:-100]

# param
rein_param = dic_param['Lens 1 b'][a1_param,0,a3_param]*0.05*scale # rein in kpc
eta_param = dic_param['Lens 1 eta'][a1_param,0,a3_param]
DM_param = PowerLaw(r,sig_crit,rein_param,eta_param)[0]

# pix
rein_pix = dic_pix['Lens b'][a1_pix,0,a3_pix]*0.05*scale # rein in kpc
eta_pix = dic_pix['Lens eta'][a1_pix,0,a3_pix]
DM_pix = PowerLaw(r,sig_crit,rein_pix,eta_pix)[0]

# add in LM!
dir = '/data/ljo31/Lens/LensModels/twoband/'
try:
    oresult = np.load(dir+name+'_212')
except:
    if name == 'J1347':
        oresult = np.load(dir+name+'_112')
    else:
        oresult = np.load(dir+name+'_211')

Mstar_param = dic_param['stellar mass'][a1_param,0,a3_param]
Mstar_pix = dic_pix['stellar mass'][a1_pix,0,a3_pix]
LM = GetStellarMass.GetStellarMass(oresult,1.,r,scale,name)
LM_param, LM_pix = LM*Mstar_param, LM*Mstar_pix

## finally, add original inference for total mass
cat = np.load('/data/ljo31/Lens/LensParams/Structure_lensgals_2src.npy')[0]
rein_TPL, eta_TPL = cat[name]['Lens 1 b']*0.05*scale, cat[name]['Lens 1 eta']
M_TPL = PowerLaw(r,sig_crit,rein_TPL,eta_TPL)[0]

# now make mass plots
pl.figure(figsize=(18,7))
pl.subplot(121)
pl.plot(0,0,color='LightGray',ls='-',label='DM+LM')
pl.plot(0,0,color='LightGray',ls=':',label='DM')
pl.plot(0,0,color='LightGray',ls='--',label='LM')
###
pl.plot(lr,DM_pix[:-100],color='SteelBlue',ls=':')
pl.plot(lr,LM_pix,color='SteelBlue',ls='--')
pl.plot(lr,LM_pix+DM_pix[:-100],color='SteelBlue',label='pixellated')
###
pl.plot(lr,DM_param[:-100],color='Crimson',ls=':')
pl.plot(lr,LM_param,color='Crimson',ls='--')
pl.plot(lr,LM_param+DM_param[:-100],color='Crimson',label='parametric')
###
pl.plot(lr,M_TPL[:-100],color='k',ls='-',label='orig param (DM+LM)')
###
pl.yscale('log')
pl.legend(loc='lower right',ncol=2,fontsize=15)
pl.ylabel(r'M (M$_{\odot}$)')
pl.xlabel('r (kpc)')
pl.axvline(rein_TPL,color='k',ls='-.',label='Einstein radius',lw=4) # Einstein radius
pl.axvline(rein_TPL*cat[name]['Lens 1 q'],color='k',ls='-.',label='Einstein radius',lw=4) # Einstein radius

pl.axvline(r_eff[names==name],color='k',ls='-.',label='effective radius') # effective radius
pl.errorbar(r_eff[names==name], Mvir[names==name],xerr = dr_eff[names==name],yerr=dMvir[names==name],color='k',marker='o',lw=4)
pl.axis([0,15,1e9,9e11])
pl.title(name)

### finally, put some info on the plot
pl.figtext(0.40,0.47,'$\eta_{param} = $'+'%.2f'%(eta_param),fontsize=15)
pl.figtext(0.40,0.42,'$\eta_{pix} = $'+'%.2f'%(eta_pix),fontsize=15)
##
pl.figtext(0.39,0.37,'$\eta_{orig} = $'+'%.2f'%eta_TPL,fontsize=15)
##
apar,apix = np.min(dic_param['stellar mass'][:,0].ravel()*10), np.min(dic_pix['stellar mass'][:,0].ravel()*10)
bpar,bpix = np.max(dic_param['stellar mass'][:,0].ravel()*10), np.max(dic_pix['stellar mass'][:,0].ravel()*10)

pl.subplot(122)
pl.hist(dic_param['stellar mass'][:,0].ravel()*10,bins=np.linspace(apar-0.05,bpar+0.05,30),alpha=0.5,label='parametric',histtype='stepfilled',normed=True)
pl.hist(dic_pix['stellar mass'][:,0].ravel()*10,bins=np.linspace(apix-0.05,bpix+0.05,30),alpha=0.5,label='pixellated',histtype='stepfilled',normed=True)
pl.legend(loc='upper left',fontsize=25)
pl.xlabel('M$_{\star}$ ( 10$^{11}$ M$_{\odot}$)')
pl.ylabel('probability density')
pl.show()
'''
## now try plotting with uncertainties!
# go through chains, which could take a while
olp,otrace,odic,_ = oresult
otrace = otrace[:,0]
for key in odic.keys():
    odic[key] = odic[key][:,0].ravel()
masses_TPL = np.zeros((odic[key].size,r.size))
for n in range(odic[key].size):
    eta, rein = odic['Lens 1 eta'][n], odic['Lens 1 b'][n]*0.05*scale
    masses_TPL[n] = PowerLaw(r,sig_crit,rein,eta)[0]
masses_TPL = np.row_stack((r,masses_TPL))
outname = '/data/ljo31b/EELs/galsub/params_from_chains/'+name+'_masses_TPL'
f = open(outname,'wb')
cPickle.dump(masses_TPL,f,2)
f.close()
masses_TPL = np.load('/data/ljo31b/EELs/galsub/params_from_chains/'+name+'_masses_TPL')

# at each radius, find median and uncertainties
masses_TPL = masses_TPL[1:]
med_TPL = np.median(masses_TPL,axis=0)
lo_TPL, hi_TPL = np.percentile(masses_TPL,16,axis=0), np.percentile(masses_TPL,84,axis=0)



# try the same for other runs!
# PL (DM+LM)
# could do both at once
masses_PL,masses_DPL = np.zeros((dic_PL[key].size,r.size)), np.zeros((dic_DPL[key].size,r.size))


for key in dic_PL.keys():
    dic_PL[key] = dic_PL[key][:,0].ravel()
for key in dic_DPL.keys():
    dic_DPL[key] = dic_DPL[key][:,0].ravel()

for n in range(dic_PL[key].size):
    eta, rein = dic_PL['Lens 1 eta'][n], dic_PL['Lens 1 b'][n]*0.05*scale
    masses_PL[n] = PowerLaw(r,sig_crit,rein,eta)[0]
    gamma, rein, rs = dic_DPL['Lens 1 gamma'][n], dic_DPL['Lens 1 b'][n]*0.05*scale, dic_DPL['Lens 1 rs'][n]*0.05*scale
    masses_DPL[n] = gNFW(r,sig_crit,rein,gamma,rs)[0]

masses_PL = np.row_stack((r,masses_PL))
masses_DPL = np.row_stack((r,masses_DPL))

outname_PL = '/data/ljo31b/EELs/galsub/params_from_chains/'+name+'_masses_PL'
outname_DPL = '/data/ljo31b/EELs/galsub/params_from_chains/'+name+'_masses_DPL'
f = open(outname_PL,'wb')
cPickle.dump(masses_PL,f,2)
f.close()
f = open(outname_DPL,'wb')
cPickle.dump(masses_DPL,f,2)
f.close()

masses_PL = np.load('/data/ljo31b/EELs/galsub/params_from_chains/'+name+'_masses_PL')
masses_DPL = np.load('/data/ljo31b/EELs/galsub/params_from_chains/'+name+'_masses_DPL')
# at each radius, find median and uncertainties
masses_PL = masses_PL[1:]
med_PL = np.median(masses_PL,axis=0)
lo_PL, hi_PL = np.percentile(masses_PL,16,axis=0), np.percentile(masses_PL,84,axis=0)
                 
masses_DPL = masses_DPL[1:]
med_DPL = np.median(masses_DPL,axis=0)
lo_DPL, hi_DPL = np.percentile(masses_DPL,16,axis=0), np.percentile(masses_DPL,84,axis=0)

Mstar_med_PL, Mstar_med_DPL = np.median(dic_PL['stellar mass']), np.median(dic_DPL['stellar mass'])
Mstar_lo_PL, Mstar_lo_DPL = np.percentile(dic_PL['stellar mass'],16), np.percentile(dic_DPL['stellar mass'],16)
Mstar_hi_PL, Mstar_hi_DPL = np.percentile(dic_PL['stellar mass'],84), np.percentile(dic_DPL['stellar mass'],84)




pl.figure()
# legend
pl.plot(0,0,color='LightGray',ls='-',label='DM+LM')
pl.plot(0,0,color='LightGray',ls=':',label='DM')
pl.plot(0,0,color='LightGray',ls='--',label='LM')
# TPL
pl.plot(r,med_TPL,color='k')
pl.fill_between(r,lo_TPL,hi_TPL,color='DarkGray')#,alpha=0.5)
# PL
# DM
pl.plot(lr,med_DPL[:-100],color='SteelBlue',ls=':')
pl.fill_between(lr,lo_DPL[:-100],hi_DPL[:-100],color='LightBlue',alpha=0.5)
# LM
pl.plot(lr,LM*Mstar_med_DPL,color='SteelBlue',ls='--')
pl.fill_between(lr,LM*Mstar_lo_DPL,LM*Mstar_hi_DPL,color='LightBlue',alpha=0.5)
# total
pl.plot(lr,med_DPL[:-100]+LM*Mstar_med_DPL,color='SteelBlue',label='gNFW')
pl.fill_between(lr, lo_DPL[:-100]+LM*Mstar_lo_DPL,hi_DPL[:-100]+LM*Mstar_hi_DPL,color='LightBlue',alpha=0.5)

# add in PL if this works!
# DM
pl.plot(lr,med_PL[:-100],color='Crimson',ls=':')
pl.fill_between(lr,lo_PL[:-100],hi_PL[:-100],color='LightPink',alpha=0.5)
# LM
pl.plot(lr,LM*Mstar_med_PL,color='Crimson',ls='--')
pl.fill_between(lr,LM*Mstar_lo_PL,LM*Mstar_hi_PL,color='LightPink',alpha=0.5)
# total
pl.plot(lr,med_PL[:-100]+LM*Mstar_med_PL,color='Crimson',label='power law')
pl.fill_between(lr, lo_PL[:-100]+LM*Mstar_lo_PL,hi_PL[:-100]+LM*Mstar_hi_PL,color='LightPink',alpha=0.5)



pl.yscale('log')
pl.legend(loc='lower right',ncol=2,fontsize=15)
pl.ylabel(r'M (M$_{\odot}$)')
pl.xlabel('r (kpc)')
pl.axvline(rein,color='k',ls='-.',label='Einstein radius',lw=4) # Einstein radius
pl.axvline(rein*cat[name]['Lens 1 q'],color='k',ls='-.',label='Einstein radius',lw=4) # Einstein radius

pl.axvline(r_eff[names==name],color='k',ls='-.',label='effective radius') # effective radius
pl.errorbar(r_eff[names==name], Mvir[names==name],xerr = dr_eff[names==name],yerr=dMvir[names==name],color='k',marker='o',lw=4)
pl.axis([0,15,1e9,9e11])

pl.show()
'''
