import PROFILES
import numpy as np, pylab as pl
from scipy.special import hyp2f1

g = 1.5

r = np.logspace(-5,5,4000)
mass = PROFILES.gNFW(x=0.,y=0.,eta=g,pa=0.,q=1.,b=1.,zl=0.3,zs=0.7,rs=1.)
s1 = mass.sigma(r) # keeton's catalogue
s2 = 2. * (1.+r**2.)**-1 * hyp2f1(1.,g/2.,1.5,(1.+r**2)**-1) # keeton again - double-checking
s3 = 2. * r**-2 * hyp2f1(1.,1.5-0.5*g,1.5,-1./r**2.) # mathematica

scale = (s1/s2)[0]

pl.figure()
pl.loglog(r,s1/scale)
pl.loglog(r,s2)
pl.loglog(r,s3)
