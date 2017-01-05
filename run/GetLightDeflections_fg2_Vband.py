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


def Deflect(model,kmodel,name):
    model.Initialise()
    kmodel.Initialise()
    
    zl,zs = lz[name][0],sz[name][0]
    sig_crit = lenslib.sig_crit(zl,zs) # this comes back in solar masses per Mpc^2 

    # recast in M_sun / pixel^2
    scale = astCalc.da(lz[name][0])*np.pi/180./3600
    sig_crit *= (0.05*scale)**2.
                
    #print name, 'has ', int(model.galno), 'galaxy components'

    # define coords for the THREE bands
    yc,xc = iT.coords(model.imgs[0].shape)
    xv,yv = xc + model.Dx, yc+model.Dy
    xi,yi = xv+model.Ddic['xoffset'], yv+model.Ddic['yoffset']
    yk,xk = iT.coords(kmodel.img.shape)
    xk,yk = xk*kmodel.pix + kmodel.Ddic['xoffset'] + kmodel.Dx, yk*kmodel.pix + kmodel.Ddic['yoffset'] + kmodel.Dy

    # get relative contributions from the I band
    fracs = [1.]
    Fracs = [1.,0.]
    if model.galno == 2.:
        gali1,gali2 = model.gals
        #print 'fits',model.fits
        ii = np.where(model.fits[0][0:2]!=0.)
        if len(ii[0])<len(model.gals):
            print 'one of the components has no flux!'
            print 'component ', '%.2f'%np.where(model.fits[0][0:2]==0), 'has no flux so we are removing it'
            model.galno = 1.
            model.gals = [model.gals[ii[0]]]
            # preserve fractions
            Fracs = [0,0]
            Fracs[ii[0]] = 1.
            print Fracs
        else:
            #print galk1.re,galk2.re
            mk1,mk2 = gali1.getMag(model.fits[0][0],0), gali2.getMag(model.fits[0][1],0)
            Fk1,Fk2 = 10**(-0.4*mk1), 10**(-0.4*mk2) 
            f1,f2 = Fk1/(Fk1+Fk2), Fk2/(Fk1+Fk2)
            print 'there are two components with relative fluxes ' '%.2f'%f1, 'and ', '%.2f'%f2
            fracs = [f1,f2]
            Fracs = fracs

    # build lenses
    if model.galno == 1.:
        gal = model.gals[0]
        lens1 = MassModels.Sersic('lens 1',{'x':gal.x,'y':gal.y,'n':gal.n,'re':gal.re,'q':gal.q,'pa':gal.pa,'b':50.})
        lens1.setbFromMass(1e12,sig_crit)
        lenses = [lens1]
    else:
        gal1,gal2 = model.gals
        lens1 = MassModels.Sersic('lens 1',{'x':gal1.x,'y':gal1.y,'n':gal1.n,'re':gal1.re,'q':gal1.q,'pa':gal1.pa,'b':50.})
        lens2 = MassModels.Sersic('lens 2',{'x':gal2.x,'y':gal2.y,'n':gal2.n,'re':gal2.re,'q':gal2.q,'pa':gal2.pa,'b':50.})
        # can set b from total mass
        lens1.setbFromMass(1e12,sig_crit)
        lens2.setbFromMass(1e12,sig_crit)
        if name == 'J2228':
            lens1.setbFromMass(5e12,sig_crit)
            lens2.setbFromMass(5e12,sig_crit)
        lenses = [lens1, lens2]

    # sum deflection angles
    '''V,I,K = [],[],[]
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
        #print n, fracs[n]
        n += 1
    V, I, K = [XV,YV],[XI,YI],[XK,YK]

    outFile = '/data/ljo31b/EELs/galsub/light_deflections_'+name+'_IBAND.npy'
    deflections = [V,I,K]
    f = open(outFile,'wb')
    cPickle.dump(deflections,f,2)
    f.close()
    #print fracs'''
    print 'FRACS',FRACS
    FRACS.append([name,Fracs])
    np.save('/data/ljo31b/EELs/galsub/fracs_Vband.npy',dict(FRACS))



dir = '/data/ljo31/Lens/LensModels/twoband/'
sz = np.load('/data/ljo31/Lens/LensParams/SourceRedshiftsUpdated.npy')[()]
lz = np.load('/data/ljo31/Lens/LensParams/LensRedshiftsUpdated.npy')[()]

bands = np.load('/data/ljo31/Lens/LensParams/HSTBands.npy')[()]

names = sz.keys()
names.sort()
     
FRACS = []
for ii in range(len(names)):
    name = names[ii]
    if name == 'J1248':
        continue
    try:
        result = np.load(dir+name+'_212')
        kresult = np.load(dir+name+'_Kp_212')
    except:
        if name == 'J1347':
            result = np.load(dir+name+'_112')
            kresult = np.load('/data/ljo31/Lens/J1347/twoband_Kp_112_2')
        elif name == 'J1619':
            result = np.load(dir+name+'_212')
            kresult = np.load(dir+name+'_Kp_212_lensandgalon')
            print 'J1619'                                                                                                                                                                                                       
        else:
            result = np.load(dir+name+'_211')
            kresult = np.load(dir+name+'_Kp_211')
            print 'here'
    model = L.EELs(result, name)
    kmodel = K.EELs(kresult,result,name)
    print name
    Deflect(model,kmodel,name)



