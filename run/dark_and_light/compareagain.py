import numpy as np, pylab as pl

result=np.load('/data/ljo31b/EELs/galsub/emceeruns/J0837_parametric_DPL_4')
lpI,traceI,dicI,_ = result

result=np.load('/data/ljo31b/EELs/galsub/emceeruns/J0837_parametric_DPL_1')
lpV,traceV,dicV,_ = result

for key in ['Lens 1 b','Lens 1 gamma','stellar mass','Lens 1 q']:
    pl.figure()
    pl.title(key)
    pl.hist(dicV[key][:,0].ravel(),30,alpha=0.5,histtype='stepfilled',label='V',normed=True)
    pl.hist(dicI[key][:,0].ravel(),30,alpha=0.5,histtype='stepfilled',label='I',normed=True)
    pl.legend(loc='upper right')
    pl.show()
