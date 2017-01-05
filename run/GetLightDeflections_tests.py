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
from imageSim import SBModels

def Deflect(model,kmodel,name):
    model.Initialise()
    kmodel.Initialise()
    print model.gals[0].re, kmodel.gals[0].re

    zl,zs = lz[name][0],sz[name][0]
    sig_crit = lenslib.sig_crit(zl,zs) # this comes back in solar masses per Mpc^2 
    #print np.log10(sig_crit)
    # recast in M_sun / pixel^2
    scale = astCalc.da(lz[name][0])*np.pi/180./3600
    sig_crit *= (0.05*scale)**2.
    #print 'sigma crit', np.log10(sig_crit)
                
    print name, 'has ', int(model.galno), 'galaxy components'

    # define coords for the THREE bands
    yc,xc = iT.coords(model.imgs[0].shape)
    xv,yv = xc + model.Dx, yc+model.Dx
    xi,yi = xv+model.Ddic['xoffset'], yv+model.Ddic['yoffset']
    yk,xk = iT.coords(kmodel.img.shape)
    xk,yk = xk*kmodel.pix + kmodel.Ddic['xoffset'] + kmodel.Dx, yk*kmodel.pix + kmodel.Ddic['yoffset'] + kmodel.Dy

    # get relative contributions from the K band
    fracs = [1.,1.]

    # build lenses
    gal = model.gals[0]
    lens1 = MassModels.Sersic('lens 1',{'x':gal.x,'y':gal.y,'n':1,'re':gal.re,'q':1,'pa':gal.pa,'b':50.})
    lens1.setbFromMass(1e12,sig_crit)
    lens2 = MassModels.ExtShear('lens 1',{'x':gal.x,'y':gal.y,'b':0.3,'pa':0.})
    lenses = [lens1,lens2]
    
    # sum deflection angles
    V,I,K = [],[],[]
    n = 0
    x0v,y0v = xv*0.,yv*0.
    x0i,y0i = xi*0.,yi*0.
    x0k,y0k = xk*0.,yk*0.
    for lens in lenses:
        lens.setPars()
        XV,YV = lens.deflections(xv,yv)
        XI,YI = lens.deflections(xi,yi)
        XK,YK = lens.deflections(xk,yk)
        x0v += XV*fracs[n]
        x0i += XI*fracs[n]
        x0k += XK*fracs[n]
        y0v += YV*fracs[n]
        y0i += YI*fracs[n]
        y0k += YK*fracs[n]
        n += 1
    V, I, K = [XV,YV],[XI,YI],[XK,YK]

    # plot to see what they look like again
    '''pl.figure()
    pl.imshow(XV,interpolation='nearest',origin='lower')
    pl.colorbar()
    pl.figure()
    pl.imshow(YV,interpolation='nearest',origin='lower')
    pl.colorbar()
    pl.show()'''

    # what happens when I put a source right behind the lens?
    src = SBModels.Sersic('source',{'x':gal.x,'y':gal.y,'n':1,'re':10.,'q':1,'pa':gal.pa})
    im = src.pixeval(xv-x0v,yv-y0v,1.,csub=23)
    pl.figure(figsize=(25,7))
    pl.subplot(131)
    pl.imshow(im,interpolation='nearest',origin='lower')
    pl.colorbar()
    pl.subplot(132)
    pl.imshow(src.pixeval(xv,yv,1.,csub=23.),interpolation='nearest',origin='lower')
    pl.colorbar()
    pl.subplot(133)
    pl.imshow(src.pixeval(x0v,y0v,1.,csub=23.),interpolation='nearest',origin='lower')
    pl.colorbar()
    pl.show()
              
    #outFile = '/data/ljo31b/EELs/galsub/light_deflections_'+name+'.npy'
    #deflections = [V,I,K]
    #f = open(outFile,'wb')
    #cPickle.dump(deflections,f,2)
    #f.close()

dir = '/data/ljo31/Lens/LensModels/twoband/'
sz = np.load('/data/ljo31/Lens/LensParams/SourceRedshiftsUpdated_0.30_source_indous_vdfit_jul2016_J2228.npy')[()]
lz = np.load('/data/ljo31/Lens/LensParams/LensRedshiftsUpdated.npy')[()]

bands = np.load('/data/ljo31/Lens/LensParams/HSTBands.npy')[()]

names = sz.keys()
names.sort()
     
#ZPs = np.load('/data/ljo31/Lens/LensParams/Keck_zeropoints.npy')[()]
 
name = 'J1347'
result = np.load(dir+name+'_112')
kresult = np.load('/data/ljo31/Lens/J1347/twoband_Kp_112_2')

model = L.EELs(result, name)
kmodel = K.EELs(kresult,result,name)
print name
Deflect(model,kmodel,name)

