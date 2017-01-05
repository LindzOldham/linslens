import numpy as np, pylab as pl, pyfits as py
import glob

name = 'J2228'

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
if name == 'J1218':
    f = '/data/ljo31b/EELs/galsub/emceeruns/J1218_parametric_0'
result_PL = np.load(f)

file = glob.glob(dir+name+'_parametric_DPL*')
file.sort()
f = file[-1]
if name == 'J1218':
    f = '/data/ljo31b/EELs/galsub/emceeruns/J1218_parametric_DPL_0'
print f
result_DPL = np.load(f)

lp_PL,trace_PL,dic_PL,_ = result_PL
lp_DPL,trace_DPL,dic_DPL,_ = result_DPL
a1_PL,a3_PL = np.unravel_index(lp_PL[:,0].argmax(),lp_PL[:,0].shape)
a1_DPL,a3_DPL = np.unravel_index(lp_DPL[:,0].argmax(),lp_DPL[:,0].shape)

for key in dic_PL.keys():
    if key == 'Lens 1 eta':
        continue
    pl.figure()
    pl.title(key)
    pl.hist(dic_PL[key][:,0].ravel(),30,normed=True,alpha=0.5,label='power law',histtype='stepfilled')
    pl.hist(dic_DPL[key][:,0].ravel(),30,normed=True,alpha=0.5,label='gNFW',histtype='stepfilled')
    pl.show()

# J1218 has very different solutions? And offset DM haloes?
