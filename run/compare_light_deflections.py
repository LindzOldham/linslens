import numpy as np, pylab as pl

fracs_I = np.load('/data/ljo31b/EELs/galsub/fracs_Iband.npy')[()]
fracs_V =  np.load('/data/ljo31b/EELs/galsub/fracs_Vband.npy')[()]
fracs_K =  np.load('/data/ljo31b/EELs/galsub/fracs.npy')[()]

keys = fracs_K.keys()
keys.sort()

for name in keys:
    print name
    print 'V', '%.3f,%.3f'%(fracs_V[name][0],fracs_V[name][1])
    print 'I', '%.3f,%.3f'%(fracs_I[name][0],fracs_I[name][1])
    print 'K', '%.3f,%.3f'%(fracs_K[name][0],fracs_K[name][1])

