import numpy,pyfits,pylab
import indexTricks as iT
from pylens import MassModels,pylens
from imageSim import SBModels,convolve
import pymc,cPickle
from scipy import optimize
import myEmcee_blobs as myEmcee #updateEmcee as myEmcee
import numpy as np, pylab as pl, pyfits as py
from linslens.GrabImages_huge import *

from matplotlib import rcParams
rcParams['xtick.direction'] = 'in'
rcParams['ytick.direction'] = 'onut'

name = 'J0901'
file = py.open('/data/ljo31b/EELs/galsub/images/'+name+'_maxlnL.fits')
img1,img2 = file[1].data, file[2].data

# make images 120 pixels across
XX=60.
my,mx = img1.shape[0]/2., img1.shape[1]/2.
img1,img2 = img1[my-XX:my+XX,mx-XX:mx+XX], img2[my-XX:my+XX,mx-XX:mx+XX]
#sig1 = py.open('/data/ljo31b/EELs/galsub/images/'+name+'_maxlnL_sig.fits')[1].data[my-XX:my+XX,mx-XX:mx+XX]

_,sig1,psf1,_,sig2,psf2,DX,DY,_,mask = EasyAddImages(name)
psf1 = psf1/np.sum(psf1)
psf2 = psf2/np.sum(psf2)
sig1,sig2,mask = sig1[my-XX:my+XX,mx-XX:mx+XX], sig2[my-XX:my+XX,mx-XX:mx+XX], mask[my-XX:my+XX,mx-XX:mx+XX]
mask = mask==0.
OVRS=1

yc,xc = iT.coords(img1.shape)
xc,yc = xc+DX+(mx-XX), yc+DY+(my-XX) 

V,I,_ = np.load('/data/ljo31b/EELs/galsub/deflections/light_deflections_'+name+'.npy')
V = [V[ii][my-XX:my+XX,mx-XX:mx+XX] for ii in range(len(V))]
I = [I[ii][my-XX:my+XX,mx-XX:mx+XX] for ii in range(len(I))]
xstars, ystars = [V[0].flatten(), I[0].flatten()], [V[1].flatten(), I[1].flatten()]
xstarms, ystarms = [V[0][mask],I[0][mask]], [V[1][mask], I[1][mask]]

result = np.load('/data/ljo31/Lens/LensModels/twoband/'+name+'_211')
#result = np.load('/data/ljo31/Lens/J0901/basickest_6temps_200walkers_120pixels_smallcov')
lp,trace,dic,_=result
a1,a3 = numpy.unravel_index(lp[:,0].argmax(),lp[:,0].shape)
a2=0.

# now - make a parametric model as the pix model is being so weird!!!
imgs = [img1,img2]
sigs = [sig1,sig2]
psfs = [psf1,psf2]
PSFs = []
for i in range(len(psfs)):
    psf = psfs[i]
    image = imgs[i]
    psf /= psf.sum()
    psf = convolve.convolve(image,psf)[1]
    PSFs.append(psf)

dxs = [0.,dic['xoffset'][a1,0,a3]]
dys = [0.,dic['yoffset'][a1,0,a3]]

# lens
vx,vy = dic['Lens 1 x'][a1,0,a3], dic['Lens 1 y'][a1,0,a3]
X = pymc.Uniform('Lens 1 x',vx-5,vx+5,value=vx)
Y = pymc.Uniform('Lens 1 y',vy-5,vy+5,vy)
B = pymc.Uniform('Lens 1 b',0.0,100.,value=dic['Lens 1 b'][a1,0,a3]/1.5)
Q = pymc.Uniform('Lens 1 q',0.1,1.0,value=dic['Galaxy 2 q'][a1,0,a3])
ETA = pymc.Uniform('Lens 1 eta',0.5,1.5,value=1.)
PA = pymc.Uniform('Lens 1 pa',-180,180.,value=dic['Lens 1 pa'][a1,0,a3])

SH = pymc.Uniform('extShear',-0.3,0.3,value=dic['extShear'][a1,0,a3])
SHPA = pymc.Uniform('extShear PA',-180.,180,value=dic['extShear PA'][a1,0,a3])

lens1 = MassModels.PowerLaw('Lens 1',{'x':X,'y':Y,'b':B,'eta':ETA,'q':Q,'pa':PA})
shear = MassModels.ExtShear('shear',{'x':X,'y':Y,'b':SH,'pa':SHPA})
lenses = [lens1,shear]
pars = [X,Y,B,Q,ETA,PA,SH,SHPA]
cov = [0.5,0.5,0.5,0.1,0.1,5.,0.05,5.]

# source
sx,sy = dic['Source 1 x'][a1,0,a3]+vx, dic['Source 1 y'][a1,0,a3]+vy
SX = pymc.Uniform('Source 1 x',sx-5,sx+5,sx)
SY = pymc.Uniform('Source 1 y',sy-5,sy+5,sy)
SR = pymc.Uniform('Source 1 re',0.5,100,dic['Source 1 re'][a1,0,a3])
SN = pymc.Uniform('Source 1 n',0.5,8,dic['Source 1 n'][a1,0,a3])
SPA = pymc.Uniform('Source 1 pa',-90,90,dic['Source 1 pa'][a1,0,a3])
SQ = pymc.Uniform('Source 1 q',0.2,1.,dic['Source 1 q'][a1,0,a3])
src = SBModels.Sersic('Source 1',{'x':SX,'y':SY,'pa':SPA,'q':SQ,'re':SR,'n':SN})
pars += [SX, SY, SR, SN, SPA, SQ]
cov += [0.5,0.5,0.5,0.2,1.,0.05]

# if needed, could simplify this further by fixing source re, n, pa and q!
#SR, SN = dic['Source 1 re'][a1,0,a3], dic['Source 1 n'][a1,0,a3]
#SPA, SQ = dic['Source 1 q'][a1,0,a3], dic['Source 1 q'][a1,0,a3]

Mstar = pymc.Uniform('stellar mass',0.0,10.,0.1)
pars += [Mstar]
cov += [0.5]

# with M/L deflections
@pymc.deterministic
def logP(value=0.,p=pars):
    lp = 0.
    for i in range(len(imgs)):
        dx,dy = dxs[i],dys[i]
        xp,yp = xc+dx,yc+dy
        image,sigma,psf = imgs[i],sigs[i],PSFs[i]
        imin,sigin,xin,yin = image[mask], sigma[mask],xp[mask],yp[mask]
        # model
        model = np.empty((2,imin.size))
        for lens in lenses:
            lens.setPars()
        x0,y0 = pylens.getDeflections(lenses,[xin,yin])
        x0 -= Mstar.value*xstarms[i]
        y0 -= Mstar.value*ystarms[i]  # check I accounted for offsets?
        src.setPars()
        tmp = xc*0.
        tmp[mask] = src.pixeval(x0,y0,1.,csub=21)
        tmp = iT.resamp(tmp,OVRS,True)
        tmp = convolve.convolve(tmp,psf,False)[0]
        model[0] = tmp[mask].ravel()
        model[1] = np.ones(model[1].size)
        # nnls
        rhs = (imin/sigin) # data
        op = (model/sigin).T # model matrix
        fit, chi = optimize.nnls(op,rhs)
        model = (model.T*fit).sum(1)
        resid = (model-imin)/sigin
        lp += -0.5*(resid**2.).sum()
    return lp

@pymc.observed
def likelihood(value=0.,lp=logP):
    return lp

cov = numpy.array(cov)
print len(pars), len(cov)

'''
S = myEmcee.PTEmcee(pars+[likelihood],cov=cov,nthreads=32,nwalkers=60,ntemps=3)
S.sample(500)

print 'done emcee'
outFile = '/data/ljo31/Lens/J0901/basic_parametric_0'

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
   

print 'das Ende'
'''
pl.figure()
pl.plot(lp[:,0])

colours = ['F555W', 'F814W']
models = []
fits = []
for i in range(len(imgs)):
    dx,dy = dxs[i],dys[i]
    xp,yp = xc+dx,yc+dy
    image,sigma,psf = imgs[i],sigs[i],PSFs[i]
    imin,sigin,xin,yin = image.flatten(), sigma.flatten(),xp.flatten(),yp.flatten()
    model = np.empty((2,imin.size))
    for lens in lenses:
        lens.setPars()
    x0,y0 = pylens.getDeflections(lenses,[xin,yin])
    x0 -= Mstar.value*xstars[i]
    y0 -= Mstar.value*ystars[i] 

    src.setPars()
    tmp = xc*0.
    tmp = src.pixeval(x0,y0,1.,csub=21).reshape(xc.shape)
    tmp = iT.resamp(tmp,1,True)
    tmp = convolve.convolve(tmp,psf,False)[0]
    model[0] = tmp.ravel()
    model[1] = np.ones(model[1].shape)

    rhs = image[mask]/sigma[mask]
    mmodel = model.reshape((2,image.shape[0],image.shape[1]))
    mmmodel = np.empty((2,image[mask].size))
    for m in range(mmodel.shape[0]):
        mmmodel[m] = mmodel[m][mask]
    op = (mmmodel/sigma[mask]).T
    rhs = image[mask]/sigma[mask]
    fit, chi = optimize.nnls(op,rhs)
    components = (model.T*fit).T.reshape((2,image.shape[0],image.shape[1]))
    model = components.sum(0)
    models.append(model)
    NotPlicely(image,model,sigma,colours[i])
    pl.show()
       

