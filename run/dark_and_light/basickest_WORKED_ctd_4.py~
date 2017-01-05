import numpy,pyfits,pylab
import indexTricks as iT
from pylens import MassModels,pylens,adaptTools as aT,pixellatedTools as pT
from imageSim import SBModels,convolve
from scipy.sparse import diags
import pymc,cPickle
from scipy import optimize
import myEmcee_blobs as myEmcee #updateEmcee as myEmcee
import numpy as np, pylab as pl, pyfits as py
from pylens import lensModel
from scipy.interpolate import RectBivariateSpline
import adaptToolsBug as BB
from linslens import EELsKeckLensModels as K, EELsImages_huge as Image, EELsModels_huge as L
import indexTricks as iT

from matplotlib import rcParams
rcParams['xtick.direction'] = 'in'
rcParams['ytick.direction'] = 'onut'

name = 'J0837'
file = py.open('/data/ljo31b/EELs/galsub/images/'+name+'.fits')
img1,img2 = file[1].data, file[2].data

# make images 120 pixels across
my,mx = img1.shape[0]/2., img1.shape[1]/2.
img1,img2 = img1[my-45:my+45,mx-45:mx+45], img2[my-45:my+45,mx-45:mx+45]

'''SIG1 = py.open('/data/ljo31/Lens/J0837/F606W_noisemap_huge.fits')[0].data.copy()[my-45:my+45,mx-45:mx+45]
PSF1 = py.open('/data/ljo31/Lens/J0837/F606W_psf1.fits')[0].data.copy()
PSF1 = PSF1/np.sum(PSF1)

MASK = py.open('/data/ljo31/Lens/J0837/mask22.fits')[0].data.copy()[my-45:my+45,mx-45:mx+45]
MASK = MASK==1'''

mask = py.open('/data/ljo31b/EELs/galsub/masks/'+name+'.fits')[0].data.copy()[my-45:my+45,mx-45:mx+45]
mask = mask==1
_,sig1,psf1,_,sig2,psf2,DX,DY,OVRS,_ = Image.J0837()  
sig1,sig2 = sig1[my-45:my+45,mx-45:mx+45], sig2[my-45:my+45,mx-45:mx+45]
#DX,DY = -100.,-100.
y,x = iT.coords(img1.shape)
x,y = x+DX+(mx-45.), y+DY+(my-45.) 

'''pl.figure()
pl.imshow(sig1-SIG1,interpolation='nearest',origin='lower')
pl.colorbar()
pl.figure()
pl.imshow(psf1-PSF1,interpolation='nearest',origin='lower')
pl.colorbar()
pl.figure()
pl.imshow(mask-MASK,interpolation='nearest',origin='lower')
pl.colorbar()
pl.show()'''

result = np.load('/data/ljo31/Lens/LensModels/twoband/'+name+'_211')
lp,trace,dic,_=result
a1,a3 = numpy.unravel_index(lp[:,0].argmax(),lp[:,0].shape)
a2=0.

Npnts = 1  # Defines `fineness' of source reconstruction (bigger is coarser)

# Function to make a `nice' plot
def showRes(x,y,src,psf,img,sig,mask,iflt,vflt,cmat,reg,niter,npix):
    oy,ox = iT.coords((npix,npix))
    oy -= oy.mean()
    ox -= ox.mean()
    span = max(x.max()-x.min(),y.max()-y.min())
    oy *= span/npix
    ox *= span/npix
    ox += x.mean()
    oy += y.mean()
    lmat = psf*src.lmat
    rmat = src.rmat
    print reg
    res,fit,model,rhs,regg = aT.getModelG(iflt,vflt,lmat,cmat,rmat,reg,niter=niter)
    print regg
    osrc = src.eval(ox.ravel(),oy.ravel(),fit).reshape(ox.shape)

    oimg = img*numpy.nan
    oimg[mask] = (lmat*fit)

    ext = [0,img.shape[1],0,img.shape[0]]
    ext2 = [x.mean()-span/2.,x.mean()+span/2.,y.mean()-span/2.,y.mean()+span/2.]
    pylab.figure()
    pylab.subplot(221)
    img[~mask] = numpy.nan
    pylab.imshow(img,origin='lower',interpolation='nearest',extent=ext,vmin=0,vmax=1,cmap='jet',aspect='auto')
    pylab.colorbar()
    pylab.subplot(222)
    pylab.imshow(oimg,origin='lower',interpolation='nearest',extent=ext,vmin=0,vmax=1,cmap='jet',aspect='auto')
    pylab.colorbar()
    pylab.subplot(223)
    pylab.imshow((img-oimg)/sig,origin='lower',interpolation='nearest',extent=ext,vmin=-3,vmax=3,cmap='jet',aspect='auto')
    pylab.colorbar()
    pylab.subplot(224)
    pylab.imshow(osrc,origin='lower',interpolation='nearest',extent=ext2,vmin=0,vmax=1,cmap='jet',aspect='auto')
    pylab.colorbar()
    return osrc

img,sig,psf = img1, sig1,psf1

cpsf = convolve.convolve(img,psf)[1]
ifltm = img[mask]
sfltm = sig[mask]
vfltm = sfltm**2
cmatm = diags(1./sfltm,0)
xm = x[mask]
ym = y[mask]
coords = [xm,ym]

PSF = pT.getPSFMatrix(psf,img.shape)
PSFm = pT.maskPSFMatrix(PSF,mask)

iflt = img.flatten()
sflt = sig.flatten()
vflt = sflt**2

xflt = x.flatten()
yflt = y.flatten()
 
src = aT.AdaptiveSource(ifltm/sfltm,ifltm.size/Npnts)

vx,vy = dic['Lens 1 x'][a1,0,a3], dic['Lens 1 y'][a1,0,a3]
X = pymc.Uniform('Lens 1 x',vx-5,vx+5,value=vx)
Y = pymc.Uniform('Lens 1 y',vy-5,vy+5,vy)
B = pymc.Uniform('Lens 1 b',0.5,100.,value=dic['Lens 1 b'][a1,0,a3])
Q = pymc.Uniform('Lens 1 q',0.1,1.0,value=dic['Lens 1 q'][a1,0,a3])
ETA = pymc.Uniform('Lens 1 eta',0.5,1.5,value=dic['Lens 1 eta'][a1,0,a3])
PA = pymc.Uniform('Lens 1 pa',-180,180.,value=dic['Lens 1 pa'][a1,0,a3])

SH = pymc.Uniform('extShear',-0.3,0.3,value=dic['extShear'][a1,0,a3])
SHPA = pymc.Uniform('extShear PA',-180.,180,value=dic['extShear PA'][a1,0,a3])

lens1 = MassModels.PowerLaw('Lens 1',{'x':X,'y':Y,'b':B,'eta':ETA,'q':Q,'pa':PA})
shear = MassModels.ExtShear('shear',{'x':X,'y':Y,'b':SH,'pa':SHPA})
lenses = [lens1,shear]
pars = [X,Y,B,Q,ETA,PA,SH,SHPA]
cov = [1.,1.,0.5,0.1,0.1,5.,0.05,5.]
#cov = [0.05,0.05,0.2,0.05,0.1,1.,0.05,0.5]

xl,yl = pylens.getDeflections(lenses,coords)
src.update(xl,yl)

reg=4.
previousResult = None

import time
def doFit(p=None,doReg=True,updateReg=True,checkImgs=True,levMar=False):
    global reg
    # Check if using levMar-style parameters
    if p is not None:
        print 'p is not none'
        for i in range(len(p)):
            pars[i].value = p[i]
            # If the parameter is out-of-bounds return a bad fit
            try:
                a = pars[i].logp
            except:
                return iflt/sflt

    for l in lenses:
        l.setPars()
    xl,yl = pylens.getDeflections(lenses,coords)

    src.update(xl,yl,doReg=doReg)
    lmat = PSFm*src.lmat
    if doReg==True:
        rmat = src.rmat
    else:
        rmat = None
    nupdate = 0
    if doReg==True and updateReg==True:
        nupdate = 10
    res,fit,model,_,regg = aT.getModelG(ifltm,vfltm,lmat,cmatm,rmat,reg,nupdate)
    reg = regg[0]
    
    if checkImgs is False:
        if levMar:
            res = res**0.5+ifltm*0.
        return -0.5*res
    # This checks is images are formed outside of the masked region
    xl,yl = pylens.getDeflections(lenses,[xflt,yflt])
    oimg,pix = src.eval(xl,yl,fit,domask=False)
    oimg = PSF*oimg
    res = (iflt-oimg)/sflt
    if levMar:
        return res
    return -0.5*(res**2).sum()


@pymc.observed
def likelihood(value=0.,tmp=pars):
    return doFit(None,True,False,True,False)

cov = numpy.array(cov)

print 'about to do doFit - i.e. get the regularisation for the current model'

doFit(None,True,False,False)
doFit(None,True,False,False)
print 'done doFit'

xl,yl = pylens.getDeflections(lenses,coords)
src.update(xl,yl)
print 'reg',reg
osrc = showRes(xl,yl,src,PSFm,img,sig,mask,ifltm,vfltm,cmatm,reg,0,400)
pylab.show()
print 'reg',reg

S = myEmcee.PTEmcee(pars+[likelihood],cov=cov,nthreads=34,nwalkers=100,ntemps=3)
S.sample(500)

print 'done emcee'
outFile = '/data/ljo31/Lens/J0837/basickest_NEWIMS_BIGIMS_3'

f = open(outFile,'wb')
cPickle.dump(S.result(),f,2)
f.close()

print 'cPickled!'

result = S.result()
lp,trace,dic,_ = result
a1,a3 = numpy.unravel_index(lp[:,0].argmax(),lp[:,0].shape)
a2=0
for i in range(len(pars)):
    pars[i].value = trace[a1,0,a3,i]
    print "%18s  %8.3f"%(pars[i].__name__,pars[i].value)

print 'reg',reg
doFit(None,True,False,False)
doFit(None,True,False,False)
print 'reg',reg
xl,yl = pylens.getDeflections(lenses,coords)
src.update(xl,yl)
osrc = showRes(xl,yl,src,PSFm,img,sig,mask,ifltm,vfltm,cmatm,reg,0,400)
reg=4.

jj=0
for jj in range(12):
    S.p0 = trace[-1]
    print 'sampling'
    S.sample(500)

    f = open(outFile,'wb')
    cPickle.dump(S.result(),f,2)
    f.close()

    result = S.result()
    lp = result[0]

    trace = numpy.array(result[1])
    a1,a3 = numpy.unravel_index(lp[:,0].argmax(),lp[:,0].shape)
    for i in range(len(pars)):
        pars[i].value = trace[a1,0,a3,i]
    print jj
    jj+=1
    doFit(None,True,False,False)
    doFit(None,True,False,False)
    print reg
    xl,yl = pylens.getDeflections(lenses,coords)
    src.update(xl,yl)
    reg = 4.

print 'das Ende'

print reg
doFit(None,True,False,False)
doFit(None,True,False,False)
print reg
xl,yl = pylens.getDeflections(lenses,coords)
src.update(xl,yl)
osrc = showRes(xl,yl,src,PSFm,img,sig,mask,ifltm,vfltm,cmatm,reg,0,400)
pylab.show()

