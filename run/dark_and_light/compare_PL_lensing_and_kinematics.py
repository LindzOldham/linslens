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
result_PL = np.load(f)

lp_PL,trace_PL,dic_PL,_ = result_PL
a1_PL,a3_PL = np.unravel_index(lp_PL[:,0].argmax(),lp_PL[:,0].shape)

   

## construct power-law sigma and mass
zl,zs = lz[name][0], sz[name][0]
sig_crit = lenslib.sig_crit(zl,zs) # solar masses per Mpc^2
sig_crit /= (1e3)**2. # solar masses per kpc^2
scale = astCalc.da(zl)*np.pi/180./3600 * 1e3
rein = np.median(dic_PL['Lens 1 b'])*0.05*scale # rein in kpc
eta = np.median(dic_PL['Lens 1 eta'])
eta_PL = np.median(dic_PL['Lens 1 eta'])
rein_PL = np.median(dic_PL['Lens 1 b'])*0.05*scale # rein in kpc

r = np.logspace(-7,4,4500)
#lr = r[:-100]

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

Mstar_PL = np.median(dic_PL['stellar mass'])#*100

print 'stellar masses:', Mstar_PL

## finally, add original inference for total mass
rein, eta = cat[name]['Lens 1 b']*0.05*scale, cat[name]['Lens 1 eta']

## get kinematics
sigma_star = np.loadtxt('/data/ljo31b/EELs/phys_models/models/sigma_star_'+name+'.dat')[()]
sigma_dm_PL = np.load('/data/ljo31b/EELs/phys_models/models/interpolators/PL_aperture_mass_measure_'+name+'.npy')[()]

# stellar mass normalised to 10^10 solar masses
sigma_PL = (Mstar_PL*sigma_star**2. + sigma_dm_PL.eval(np.array([eta_PL]))*rein_PL**(2.-eta_PL))**0.5
sigma_TPL = (sigma_dm_PL.eval(np.array([eta]))*rein**(2.-eta))**0.5

print '####'
print sigma_PL
#print (Mstar_PL*sigma_star**2.)**0.5
#print (sigma_dm_PL.eval(np.array([eta_PL]))*rein_PL**(2.-eta_PL))**0.5
print sigma_TPL

# now read in inference and compare!
dir3 = '/data/ljo31b/EELs/esi/kinematics/inference/vdfit/NEW/'
result_K = np.load(dir3+name+'_1.00_lens_esi_indous_vdfit_LENS')
lp_K, trace_K, dic_K, _ = result_K
s_lens = np.median(dic_K['lens dispersion'])
print s_lens
