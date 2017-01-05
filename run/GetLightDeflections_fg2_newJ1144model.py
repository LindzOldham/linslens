from linslens import EELsModels_huge as L
import numpy as np
import pylab as pl
from linslens import EELsKeckModels as K
from linslens.Plotter import *
import lenslib
from pylens import MassModels,pylens
import indexTricks as iT
import cPickle
from astLib import astCalc

''' assuming we have been able to make one-component models for everything! 
Either: assume a b_ein, do both, fixing ratio to be that of the photometry.  Or: do one-component models for everything. Actually, isn't the former just as good as the latter?! Yes it is.
And!!! We don't even need a grid. Because it just scales in amplitude. That's the only place b comes into it. Set b equal to unity, then scale. But we need to keep the deflection angles for each COMPONENT separate, because their amplitudes also depend on b.'''

def clip(arr,nsig=3.5):
    a = arr.flatten()
    m,s,l = a.mean(),a.std(),a.size
    while 1:
        a = a[abs(a-m)<s*nsig]
        if a.size==l:
            return m,s
        m,s,l = a.mean(),a.std(),a.size


def Deflect(model,name):
    model.Initialise()
    model.GetFits(plotresid=True)
    pl.show()
    zl,zs = lz[name][0],sz[name][0]
    sig_crit = lenslib.sig_crit(zl,zs) # this comes back in solar masses per Mpc^2 

    # recast in M_sun / pixel^2
    scale = astCalc.da(lz[name][0])*np.pi/180./3600
    sig_crit *= (0.05*scale)**2.
                
    print name, 'has ', int(model.galno), 'galaxy components'

    # define coords for the THREE bands
    yc,xc = iT.coords(model.imgs[0].shape)
    xv,yv = xc + model.Dx, yc+model.Dy
    xi,yi = xv+model.Ddic['xoffset'], yv+model.Ddic['yoffset']

    # get relative contributions from the K band
    fracs = [1.]
   
    # build lenses
    gal = model.gals[0]
    lens1 = MassModels.Sersic('lens 1',{'x':gal.x,'y':gal.y,'n':gal.n,'re':gal.re,'q':gal.q,'pa':gal.pa,'b':50.})
    lens1.setbFromMass(1e12,sig_crit)
    lenses = [lens1]


    # sum deflection angles
    V,I = [],[]
    n = 0
    x0v,y0v = xv*0.,yv*0.
    x0i,y0i = xi*0.,yi*0.
    for lens in lenses:
        lens.setPars()
        XV,YV = lens.deflections(xv,yv)
        XI,YI = lens.deflections(xi,yi)
        x0v += XV*fracs[n]
        x0i += XI*fracs[n]
        y0v += YV*fracs[n]
        y0i += YI*fracs[n]
        print n, fracs[n]
        n += 1
    V, I = [XV,YV],[XI,YI]

    pl.figure(figsize=(12,7))
    pl.subplot(121)
    pl.imshow(model.img1,interpolation='nearest',origin='lower',vmin=0,vmax=0.5)
    pl.colorbar()
    pl.subplot(122)
    pl.imshow(V[0],interpolation='nearest',origin='lower')
    pl.colorbar()
    pl.show()

    outFile = '/data/ljo31b/EELs/galsub/light_deflections_J1144_newmodel.npy'
    deflections = [V,I]
    f = open(outFile,'wb')
    cPickle.dump(deflections,f,2)
    f.close()
    #print fracs
    #print 'FRACS',FRACS
    #print 'FRACS',FRACS
    #np.save('/data/ljo31b/EELs/galsub/fracs.npy',dict(FRACS))


dir = '/data/ljo31/Lens/LensModels/twoband/'
sz = np.load('/data/ljo31/Lens/LensParams/SourceRedshiftsUpdated.npy')[()]
lz = np.load('/data/ljo31/Lens/LensParams/LensRedshiftsUpdated.npy')[()]

bands = np.load('/data/ljo31/Lens/LensParams/HSTBands.npy')[()]

names = sz.keys()
names.sort()
     
#ZPs = np.load('/data/ljo31/Lens/LensParams/Keck_zeropoints.npy')[()]
FRACS = []

result = np.load('/data/ljo31/Lens/J1144/twoband_darkandlightprep')
model = L.EELs(result, 'J1144')

Deflect(model,'J1144')


