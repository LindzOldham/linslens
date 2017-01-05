import numpy as np, pylab as pl
from linslens import EELsKeckLensModels as K, EELsImages_huge as Image, EELsModels_huge as L
from linslens.Plotter import *
from linslens.pixplot import *
import indexTricks as iT
from pylens import MassModels,pylens,adaptTools as aT,pixellatedTools as pT
from imageSim import SBModels,convolve
from scipy.sparse import diags
import pymc,cPickle
import myEmcee_blobs as myEmcee

''' lets do these one at a time for now, to avoid errors '''
''' start off with power laws '''
''' pixellated sources '''
''' band by band for now? pixellated sources on many bands? '''

def MakeModel(name,result,Npnts=2):

    # read model so far, then put it back in
    lp,trace,dic,_ = result
    a2=0
    a1,a3 = np.unravel_index(lp[:,0].argmax(),lp[:,0].shape)

    x_start, y_start = dic['Lens x'][a1,a2,a3], dic['Lens y'][a1,a2,a3]

    # now set up priors
    LX = pymc.Uniform('Lens x',x_start-10,x_start+10,x_start)
    LY = pymc.Uniform('Lens y',y_start-10,y_start+10,y_start)
    LB = pymc.Uniform('Lens b',0.,100.,dic['Lens b'][a1,a2,a3])
    LETA = pymc.Uniform('Lens eta',0.5,2.5,dic['Lens eta'][a1,a2,a3])
    LQ = pymc.Uniform('Lens q',0.1,1.,dic['Lens q'][a1,a2,a3])
    LPA = pymc.Uniform('Lens pa',-180,180,dic['Lens pa'][a1,a2,a3])
    SH = pymc.Uniform('shear',-0.3,0.3,dic['shear'][a1,a2,a3])
    SHPA = pymc.Uniform('shear pa',-180,180,dic['shear pa'][a1,a2,a3])
    Mstar = pymc.Uniform('stellar mass',0.01,10.,dic['stellar mass'][a1,a2,a3])

    lens = MassModels.PowerLaw('lens',{'x':LX,'y':LY,'b':LB,'eta':LETA,'q':LQ,'pa':LPA})
    shear = MassModels.ExtShear('shear',{'x':LX,'y':LY,'b':SH,'pa':SHPA})
    lenses = [lens,shear]

    pars = [LX,LY,LB,LETA,LQ,LPA,SH,SHPA,Mstar]
    cov = np.array([0.5,0.5,2.,0.2,0.1,10.,0.01,10.,2.])

    # now set up the inference
    imgs = [img1,img2]
    sigs = [sig1,sig2]
    ifltms = [img[pix_mask] for img in imgs]
    sfltms = [sig[pix_mask] for sig in sigs]
    vfltms = [sfltm**2 for sfltm in sfltms]
    cmatms = [diags(1./sfltm,0) for sfltm in sfltms]
    xm,ym = x[pix_mask],y[pix_mask]
    coords = [[xm,ym],[xm+vmodel.Ddic['xoffset'],ym+vmodel.Ddic['yoffset']]]

    # stellar mass deflection angles
    x_stars, y_stars = [V[0].flatten(), I[0].flatten()], [V[1].flatten(), I[1].flatten()]
    x_starms, y_starms = [V[0][pix_mask],I[0][pix_mask]], [V[1][pix_mask], I[1][pix_mask]]

    PSFs = [pT.getPSFMatrix(psf, img1.shape) for psf in [psf1,psf2]]
    PSFms = [pT.maskPSFMatrix(PSF,pix_mask) for PSF in PSFs]
    
    iflts = [img1.flatten(),img2.flatten()]
    sflts = [sig1.flatten(),sig2.flatten()]
    vflts = [sflt**2. for sflt in sflts]
    xflt,yflt = x.flatten(), y.flatten()

    srcs = []
    for ii in range(len(iflts)):
        srcs.append(aT.AdaptiveSource(ifltms[ii]/sfltms[ii],ifltms[ii].size/Npnts))
        xl,yl = pylens.getDeflections(lenses,coords[ii])
        xl,yl = xl - x_starms[ii]*Mstar.value, yl - y_starms[ii]*Mstar.value
        srcs[ii].update(xl,yl)

    import time
    reg=1.

    def doFit(p=None,updateReg=False,checkImgs=False,doReg=True):
        global reg
        reg=1.
        lp = 0.
        
        for l in lenses:
            l.setPars()
        for ii in range(len(ifltms)):
            src,ifltm,sfltm,vfltm,PSFm,cmatm = srcs[ii],ifltms[ii],sfltms[ii],vfltms[ii],PSFms[ii],cmatms[ii]
            PSF,coord,iflt,sflt = PSFs[ii],coords[ii],iflts[ii],sflts[ii]
            x_star,y_star,x_starm,y_starm =x_stars[ii],y_stars[ii],x_starms[ii],y_starms[ii]
            
            xl,yl = pylens.getDeflections(lenses,coord)
            xl,yl = xl - x_starm*Mstar.value, yl - y_starm*Mstar.value

            src.update(xl,yl,doReg=doReg)

            lmat = PSFm*src.lmat
            rmat = src.rmat

            nupdate = 0
            if updateReg==True:
                nupdate = 10

            res,fit,model,_,regg = aT.getModelG(ifltm,vfltm,lmat,cmatm,rmat,reg,nupdate)   
            reg = regg[0]
            if checkImgs is False:
                lp += -0.5*res

            else:
                xl,yl = pylens.getDeflections(lenses,[xflt,yflt])
                xl,yl = xl - x_star*Mstar.value, yl - y_star*Mstar.value
                oimg,pix = src.eval(xl,yl,fit,domask=False)
                oimg = PSF*oimg
                res = (iflt-oimg)/sflt
                lp+= -0.5*(res**2).sum()
        return lp
            

    @pymc.observed
    def likelihood(value=0.,tmp=pars):
        return doFit(None,False,True,True)

    # check initial model
    doFit(False,True)
    doFit(False,True)
    print 'current regularisation (set by hand): ', '%.1f'%reg
    
    vms = [2.,6.]
    for ii in range(len(imgs)):
        src,ifltm,sfltm,vfltm,PSFm,cmatm = srcs[ii],ifltms[ii],sfltms[ii],vfltms[ii],PSFms[ii],cmatms[ii]
        x_star,y_star,x_starm,y_starm,PSF =x_stars[ii],y_stars[ii],x_starms[ii],y_starms[ii],PSFs[ii]
        img,sig,coord = imgs[ii],sigs[ii],coords[ii]
            
        xl,yl = pylens.getDeflections(lenses,coord)
        xl,yl = xl - x_starm*Mstar.value, yl - y_starm*Mstar.value
        print xl.shape
        osrc = showRes(xl,yl,src,PSFm,img,sig,pix_mask,ifltm,vfltm,cmatm,reg,0,400,vmax_src=vms[ii])
        pl.show()

        # just check deflections all look sensible
        check=False
        if check:
            pl.figure(figsize=(15,5))
            pl.subplot(131)
            pl.imshow(lenses[0].deflections(x,y)[0],origin='lower',interpolation='nearest')
            pl.title('power law')
            pl.colorbar()
            pl.subplot(132)
            pl.imshow(V[0]*0.1,origin='lower',interpolation='nearest')
            pl.title('sersic')
            pl.colorbar()
            pl.subplot(133)
            pl.imshow(lenses[1].deflections(x,y)[0],origin='lower',interpolation='nearest')
            pl.colorbar()
            pl.title('shear')
            pl.suptitle('x deflections')
            ####
            pl.figure(figsize=(15,5))
            pl.subplot(131)
            pl.imshow(lenses[0].deflections(x,y)[1],origin='lower',interpolation='nearest')
            pl.title('power law')
            pl.colorbar()
            pl.subplot(132)
            pl.imshow(V[1]*0.1,origin='lower',interpolation='nearest')
            pl.title('sersic')
            pl.colorbar()
            pl.subplot(133)
            pl.imshow(lenses[1].deflections(x,y)[1],origin='lower',interpolation='nearest')
            pl.colorbar()
            pl.title('shear')
            pl.suptitle('y deflections')
        
            pl.show()
        
    pl.figure()
    pl.plot(lp[:,0])
    pl.show()
        
                  
    #return
    # now do some sampling. This will be slow -- make mask tighter?
    S = myEmcee.PTEmcee(pars+[likelihood],cov=cov,nthreads=24,nwalkers=300,ntemps=3)
    S.sample(500)

    outFile = '/data/ljo31b/EELs/galsub/emceeruns/'+name+'_'+str(X)
    f = open(outFile,'wb')
    cPickle.dump(S.result(),f,2)
    f.close()

    result = S.result()
    lp,trace,dic,_ = result
    a1,a3 = numpy.unravel_index(lp[:,0].argmax(),lp[:,0].shape)
    a2=0
    for i in range(len(pars)):
        pars[i].value = trace[a1,0,a3,i]
        print "%18s  %8.3f"%(pars[i].__name__,pars[i].value)

    # check inferred model
    '''doFit(False,True)
    doFit(False,True)
    print 'current regularisation (set by hand): ', '%.1f'%reg
    
    for ii in range(len(imgs)):
        src,ifltm,sfltm,vfltm,PSFm,cmatm = srcs[ii],ifltms[ii],sfltms[ii],vfltms[ii],PSFms[ii],cmatms[ii]
        x_star,y_star,x_starm,y_starm,PSF =x_stars[ii],y_stars[ii],x_starms[ii],y_starms[ii],PSFs[ii]
        img,sig,coord = imgs[ii],sigs[ii],coords[ii]
            
        xl,yl = pylens.getDeflections(lenses,coord)
        xl,yl = xl - x_starm*Mstar.value, yl - y_starm*Mstar.value

        osrc = showRes(xl,yl,src,PSFm,img,sig,pix_mask,ifltm,vfltm,cmatm,reg,0,400)
        pl.show()'''

    # re-sample a load of times
    kk=0
    for kk in range(10):
        S.p0 = trace[-1]
        S.sample(500)
        
        f = open(outFile,'wb')
        cPickle.dump(S.result(),f,2)
        f.close()

        result = S.result()
        lp,trace,dic,_=result
        
        a1,a3 = numpy.unravel_index(lp[:,0].argmax(),lp[:,0].shape)
        for i in range(len(pars)):
            pars[i].value = trace[a1,0,a3,i]
        print kk
        kk+=1
        
        # check inferred model each time!
        '''doFit(False,True)
        doFit(False,True)
        print 'current regularisation (set by hand): ', '%.1f'%reg
    
        for ii in range(len(imgs)):
            src,ifltm,sfltm,vfltm,PSFm,cmatm = srcs[ii],ifltms[ii],sfltms[ii],vfltms[ii],PSFms[ii],cmatms[ii]
            x_star,y_star,x_starm,y_starm,PSF =x_stars[ii],y_stars[ii],x_starms[ii],y_starms[ii],PSFs[ii]
            img,sig,coord = imgs[ii],sigs[ii],coords[ii]
            
            xl,yl = pylens.getDeflections(lenses,coord)
            xl,yl = xl - x_starm*Mstar.value, yl - y_starm*Mstar.value

            osrc = showRes(xl,yl,src,PSFm,img,sig,pix_mask,ifltm,vfltm,cmatm,reg,0,400)
            pl.show()'''


# make image smaller
name = 'J0901'
X = 2
file = py.open('/data/ljo31b/EELs/galsub/images/'+name+'.fits')
img1,img2 = file[1].data, file[2].data

# make images 160 pixels across
my,mx = img1.shape[0]/2., img1.shape[1]/2.
img1,img2 = img1[my-80:my+80,mx-80:mx+80], img2[my-80:my+80,mx-80:mx+80] 
pix_mask = py.open('/data/ljo31b/EELs/galsub/masks/'+name+'.fits')[0].data.copy()[my-80:my+80,mx-80:mx+80]
pix_mask = pix_mask==1
_,sig1,psf1,_,sig2,psf2,DX,DY,OVRS,_ = Image.J0901()  
sig1,sig2 = sig1[my-80:my+80,mx-80:mx+80], sig2[my-80:my+80,mx-80:mx+80]
y,x = iT.coords(img1.shape)
x,y = x+DX+(mx-80.), y+DY+(my-80.) 


# Start off DM with same q and pa as light. Try to make it as automated as possible
# stellar mass deflection angles
V,I,_ = np.load('/data/ljo31b/EELs/galsub/deflections/light_deflections_'+name+'.npy')
# make these the right shapes!
V = [V[ii][my-80:my+80,mx-80:mx+80] for ii in range(len(V))]
I = [I[ii][my-80:my+80,mx-80:mx+80] for ii in range(len(I))]

masses = np.load('/data/ljo31b/EELs/inference/new/huge/masses_212.npy')
names = ['J0837','J0901','J0913','J1125','J1144','J1218','J1323','J1347','J1446','J1605','J1606','J1619','J2228']
logM = masses[0]
# load results for making model
dir = '/data/ljo31/Lens/LensModels/twoband/'
try:
    result = np.load(dir+name+'_212')
except:
    if name == 'J1347':
        result = np.load(dir+name+'_112')
    elif name == 'J1619':
        result = np.load(dir+name+'_212')
    else:
        result = np.load(dir+name+'_211')

vmodel = L.EELs(result,name)
vmodel.Initialise()

# load results for making model
result = np.load('/data/ljo31b/EELs/galsub/emceeruns/'+name+'_1')
lp,trace,dic,_ = result

MakeModel(name,result)
