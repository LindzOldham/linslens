from linslens.Profiles import *
from astLib import astCalc
import glob, numpy as np, pylab as pl
import lenslib
from jeans.makemodel import deproject
import GetStellarMass
from scipy.interpolate import splrep, splev, splint

names = ['J0837','J0901','J0913','J1125','J1144','J1218','J1323','J1347','J1446','J1605','J1606','J1619','J2228']
lz = np.load('/data/ljo31/Lens/LensParams/LensRedshiftsUpdated.npy')[()]
sz = np.load('/data/ljo31/Lens/LensParams/SourceRedshiftsUpdated.npy')[()]

name = 'J0901'

dir = '/data/ljo31b/EELs/galsub/emceeruns/'
file = glob.glob(dir+name+'_parametric_*')
file.sort()
f = file[-1]
if 'DPL' in f:
    f = file[-2]

result_PL = np.load(f)

file = glob.glob(dir+name+'_parametric_DPL*')
file.sort()
f = file[-1]
result_DPL = np.load(f)

lp_PL,trace_PL,dic_PL,_ = result_PL
lp_DPL,trace_DPL,dic_DPL,_ = result_DPL
a1_PL,a3_PL = np.unravel_index(lp_PL[:,0].argmax(),lp_PL[:,0].shape)
a1_DPL,a3_DPL = np.unravel_index(lp_DPL[:,0].argmax(),lp_DPL[:,0].shape)

## quantities
zl,zs = lz[name][0], sz[name][0]
scale = astCalc.da(zl)*np.pi/180./3600 * 1e3
eta = dic_PL['Lens 1 eta'][a1_PL,0,a3_PL]
eta_l = np.percentile(dic_PL['Lens 1 eta'][:,0].ravel(),16)
eta_u = np.percentile(dic_PL['Lens 1 eta'][:,0].ravel(),84)


gamma = dic_DPL['Lens 1 gamma'][a1_DPL,0,a3_DPL]
rs = dic_DPL['Lens 1 rs'][a1_DPL,0,a3_DPL]*0.05*scale
gamma_l = np.percentile(dic_DPL['Lens 1 gamma'][:,0].ravel(),16)
gamma_u = np.percentile(dic_DPL['Lens 1 gamma'][:,0].ravel(),84)
rs_l = np.percentile(dic_DPL['Lens 1 rs'][:,0].ravel(),16)
rs_u = np.percentile(dic_DPL['Lens 1 rs'][:,0].ravel(),84)

r = np.logspace(-7,3,3500)
R = np.linspace(1,15,1000)

dpdr_PL = dlogrho_dlogr_PL(R,eta)
dpdr_PL_l = dlogrho_dlogr_PL(R,eta_l)
dpdr_PL_u = dlogrho_dlogr_PL(R,eta_u)


dpdr_DPL = dlogrho_dlogr_gNFW(R,gamma,rs)
dpdr_DPL_l = dlogrho_dlogr_gNFW(R,gamma_l,rs)
dpdr_DPL_u = dlogrho_dlogr_gNFW(R,gamma_u,rs)


print np.mean(dpdr_DPL[R<3.])

# mass-weighted mean slope within the einstein radius?
sig_crit = lenslib.sig_crit(zl,zs) # solar masses per Mpc^2
sig_crit /= (1e3)**2. # solar masses per kpc^2

rein = dic_PL['Lens 1 b'][a1_PL,0,a3_PL]*0.05*scale # rein in kpc
eta = dic_PL['Lens 1 eta'][a1_PL,0,a3_PL]
DM_PL = np.array(PowerLaw(r,sig_crit,rein,eta)[0])

rein = dic_DPL['Lens 1 b'][a1_DPL,0,a3_DPL]*0.05*scale 
DM_DPL = np.array(gNFW(r,sig_crit,rein,gamma,rs)[0])

# try in log scale?
model_PL = splrep(r,DM_PL)
DM_PL = splev(R,model_PL)
model_DPL = splrep(r,DM_DPL)
DM_DPL = splev(R,model_DPL)


print DM_DPL.sum(), DM_PL.sum()

#m1 = [DM_DPL[ii] * dpdr_DPL[ii] / DM_DPL[:ii].sum() for ii in range(1,len(R))]
#m2 = [DM_PL[ii] * dpdr_PL[ii] / DM_PL[:ii].sum() for ii in range(1,len(R))]
cat = np.load('/data/ljo31/Lens/LensParams/Structure_lensgals_2src.npy')[0]
rein = cat[name]['Lens 1 b']*0.05*scale
q = cat[name]['Lens 1 q']

m1 = np.sum(DM_DPL[R<rein*q] * dpdr_DPL[R<rein*q]) / DM_DPL[R<rein*q].sum()
m2 = np.sum(DM_PL[R<rein*q] * dpdr_PL[R<rein*q]) / DM_PL[R<rein*q].sum()

print m1, m2

'''
pl.figure()
pl.plot(R[1:],m1,color='Crimson',label='power law')
pl.plot(R[1:],m2,color='SteelBlue',label='gNFW')
pl.xlim([0,15])

pl.show()
'''
