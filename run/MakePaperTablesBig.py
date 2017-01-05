import numpy as np, pyfits as py
from astLib import astCalc

table = np.load('/data/ljo31/Lens/LensParams/Structure_1src.npy')[()]
m,l,u = table
phot = py.open('/data/ljo31/Lens/LensParams/Phot_1src_new.fits')[1].data
kphot = np.load('/data/ljo31/Lens/LensParams/KeckPhot_1src_new_dict.npy')[()]
names = phot['name']
sz = np.load('/data/ljo31/Lens/LensParams/SourceRedshifts.npy')[()]
lz = np.load('/data/ljo31/Lens/LensParams/LensRedshifts.npy')[()]

print r'lens & $z_l$ & $z_s$ & $r_{Ein}$ & $\eta$ & $q_{lens}$ & $\phi_{lens}$ & $\gamma_{ext}$ & $\phi_{ext}$ & $V$ (mag) & $I$ (mag) & $K$ (mag) & $n$ & $R_e$ (kpc) & $\phi$ (deg) & $q$ & $\mu$ \\'


for ii in range(len(names)):
    name = names[ii]
    print name, '&', lz[name], '&', sz[name], '& $', '%.2f'%m[name]['Lens 1 b'], '\pm', '%.2f'%l[name]['Lens 1 b'], '$ & $', '%.2f'%m[name]['Lens 1 eta'], '\pm', '%.2f'%l[name]['Lens 1 eta'], '$ & $', '%.2f'%m[name]['Lens 1 q'], '\pm', '%.2f'%l[name]['Lens 1 q'], '$ & $','%.2f'%m[name]['Lens 1 pa'], '\pm', '%.2f'%l[name]['Lens 1 pa'], '$ & $', '%.2f'%m[name]['extShear'], '\pm', '%.2f'%l[name]['extShear'], '$ & $', '%.2f'%m[name]['extShear PA'], '\pm', '%.2f'%l[name]['extShear PA'], '$ & $','%.2f'%phot['mag v'][ii], '\pm', '%.2f'%np.min((phot['mag v lo'],phot['mag v hi'])),'$ & $', '%.2f'%phot['mag i'][ii], '\pm', '%.2f'%np.min((phot['mag i lo'],phot['mag i hi'])),'$ & $', '%.2f'%kphot['mag k'][ii], '\pm', '%.2f'%np.min((kphot['mag k lo'],kphot['mag k hi'])),'$ & $', '%.2f'%m[name]['Source 1 n'], '\pm', '%.2f'%l[name]['Source 1 n'], '$ & $', '%.2f'%m[name]['Source 1 re'], '\pm', '%.2f'%l[name]['Source 1 re'], '$ & $', '%.2f'%m[name]['Source 1 pa'], '\pm', '%.2f'%l[name]['Source 1 pa'], '$ & $', '%.2f'%m[name]['Source 1 q'], '\pm', '%.2f'%l[name]['Source 1 q'], '$ & $',r'\\'

# make two tables:
#lens models and source models
table = np.load('/data/ljo31/Lens/LensParams/Structure_2src.npy')[()]
m,l,u = table


print r'lens & $z_l$ & $z_s$ & $r_{Ein}$ & $\eta$ & $q_{lens}$ & $\phi_{lens}$ & $\gamma_{ext}$ & $\phi_{ext}$ \\'
for ii in range(len(names)):
    name = names[ii]
    Da = astCalc.da(lz[name])
    scale = Da*1e3*np.pi/180./3600.
    print name, '&', lz[name], '&', sz[name], '& $', '%.2f'%(scale*0.05*m[name]['Lens 1 b']), '\pm', '%.2f'%(scale*0.05*l[name]['Lens 1 b']), '$ & $', '%.2f'%m[name]['Lens 1 eta'], '\pm', '%.2f'%l[name]['Lens 1 eta'], '$ & $', '%.2f'%m[name]['Lens 1 q'], '\pm', '%.2f'%l[name]['Lens 1 q'], '$ & $','%.2f'%m[name]['Lens 1 pa'], '\pm', '%.2f'%l[name]['Lens 1 pa'], '$ & $', '%.2f'%m[name]['extShear'], '\pm', '%.2f'%l[name]['extShear'], '$ & $', '%.2f'%m[name]['extShear PA'], '\pm', '%.2f'%l[name]['extShear PA'], r'$ \\'


phot = py.open('/data/ljo31/Lens/LensParams/Phot_1src_new.fits')[1].data
kphot = np.load('/data/ljo31/Lens/LensParams/KeckPhot_1src_new_dict.npy')[()]
phot2 = py.open('/data/ljo31/Lens/LensParams/Phot_2src_new.fits')[1].data
kphot2 = np.load('/data/ljo31/Lens/LensParams/KeckPhot_2src_new_dict.npy')[()]
table = np.load('/data/ljo31/Lens/LensParams/Structure_1src.npy')[()]
m,l,u = table
table2 = np.load('/data/ljo31/Lens/LensParams/Structure_2src.npy')[()]
m2,l2,u2 = table2
for i in range(3):
    print '%%%'

print r'lens & $V$ (mag) & $I$ (mag) & $K$ (mag) & $n$ (X model)& $R_e$ (X model; kpc) &$R_e$ (Y model; kpc) & $\phi$ (deg) & $q$ \\'
for ii in range(len(names)):
    name = names[ii]
    Da = astCalc.da(sz[name])
    scale = Da*1e3*np.pi/180./3600.
    print name, '& $','%.2f'%phot2['mag v'][ii], '\pm', '%.2f'%np.min((phot2['mag v lo'],phot['mag v hi'])),'$ & $', '%.2f'%phot2['mag i'][ii], '\pm', '%.2f'%np.min((phot2['mag i lo'],phot2['mag i hi'])),'$ & $', '%.2f'%kphot2['mag k'][ii], '\pm', '%.2f'%np.min((kphot2['mag k lo'],kphot2['mag k hi'])),'$ & $', '%.2f'%m[name]['Source 1 n'], '\pm', '%.2f'%l[name]['Source 1 n'], '$ & $', '%.2f'%(scale*0.05*m[name]['Source 1 re']), '\pm', '%.2f'%(scale*0.05*l[name]['Source 1 re']), '$ & $', '%.2f'%phot2['Re v'][ii], '\pm', '%.2f'%(np.min((phot2['Re v lo'][ii],phot2['Re v hi'][ii]))), '$ & $','%.2f'%m[name]['Source 1 pa'], '\pm', '%.2f'%l[name]['Source 1 pa'], '$ & $', '%.2f'%m[name]['Source 1 q'], '\pm', '%.2f'%l[name]['Source 1 q'], '$',r'\\'
