import numpy as np, pylab as pl
from linslens.Plotter import *
import indexTricks as iT
from pylens import pylens
from imageSim import SBObjects,convolve
from scipy.sparse import diags
import pymc,cPickle
import myEmcee_blobs as myEmcee
from linslens.GrabImages_huge import *
from MWApython.pylens import MassModels
import time


def MakeModel(name,result,oldResult=None):
    print 'starting a dark+light parametric model for ', name, '!'

    lp,trace,dic,_=result
    if len(lp.shape)==3.:
        lp = lp[:,0]
        for key in dic.keys():
            dic[key] = dic[key][:,0]
    a1,a2 = np.unravel_index(lp.argmax(),lp.shape)
    
    # set up data
    imgs = [img1]#,img2]
    sigs = [sig1]#,sig2]
    psfs = [psf1]#,psf2]
    PSFs = []
    for i in range(len(psfs)):
        psf = psfs[i]
        image = imgs[i]
        psf /= psf.sum()
        psf = convolve.convolve(image,psf)[1]
        PSFs.append(psf)

    # set up parameters
    dxs = [0.,dic['xoffset'][a1,a2]]
    dys = [0.,dic['yoffset'][a1,a2]]
    #print dic['Lens 1 b'][a1,a2]*0.05*5.
    # lens - eGNFW with n = 3
    vx,vy = dic['Lens 1 x'][a1,a2], dic['Lens 1 y'][a1,a2]
    X = pymc.Uniform('Lens 1 x',vx-10,vx+10,value=vx)
    Y = pymc.Uniform('Lens 1 y',vy-10,vy+10,vy)
    B = pymc.Uniform('Lens 1 b',0.05,200.,value=dic['Lens 1 b'][a1,a2]/2.0)
    Q = pymc.Uniform('Lens 1 q',0.1,1.0,value=dic['Lens 1 q'][a1,a2])
    GAMMA = pymc.Uniform('Lens 1 gamma',0.05,2.75,value=1.)#dic['Lens 1 eta'][a1,a2])
    RS = pymc.Uniform('Lens 1 rs',0.05,200.,value=10.)#dic['Lens 1 b'][a1,a2])
    PA = pymc.Uniform('Lens 1 pa',90.,450.,value=dic['Lens 1 pa'][a1,a2]+360.)
    SH = pymc.Uniform('extShear',-0.3,0.3,value=dic['extShear'][a1,a2])
    SHPA = pymc.Uniform('extShear PA',-180.,180,value=dic['extShear PA'][a1,a2])
    #print Q.value
    lens1 = MassModels.eGNFW('Lens 1',{'x':X,'y':Y,'b':B,'gammain':GAMMA,'rs':RS,'q':Q,'pa':PA,'trunc':1e-5})
    shear = MassModels.ExtShear('shear',{'x':X,'y':Y,'b':SH,'pa':SHPA})
    lenses = [lens1,shear]
    pars = [X,Y,B,Q,GAMMA,RS,PA,SH,SHPA]
    cov = [0.5,0.5,5.,0.1,0.5,20.,5.,0.05,5.]

    # source 1
    if 'Source 1 x' in dic.keys():
        sx,sy = dic['Source 1 x'][a1,a2]+vx, dic['Source 1 y'][a1,a2]+vy
    elif name in ['J1323','J1347']:
        sx,sy = dic['Source 2 x'][a1,a2], dic['Source 2 y'][a1,a2]
    else:
        sx,sy = dic['Source 2 x'][a1,a2]+vx, dic['Source 2 y'][a1,a2]+vy
    print sx,sy
    SX = pymc.Uniform('Source 1 x',sx-20,sx+20,sx)
    SY = pymc.Uniform('Source 1 y',sy-20,sy+20,sy)
    SR = pymc.Uniform('Source 1 re',0.5,100,dic['Source 1 re'][a1,a2])
    SN = pymc.Uniform('Source 1 n',0.5,9,dic['Source 1 n'][a1,a2])
    try:
        spa = dic['Source 1 pa'][a1,a2]
        SPA = pymc.Uniform('Source 1 pa',spa-90,spa+90,spa)
    except:
        spa = dic['Source 2 pa'][a1,a2]
        SPA = pymc.Uniform('Source 1 pa',spa-90,spa+90,spa)
        print 'there are two source components with PAs fixed to be the same'
    print 
    SQ = pymc.Uniform('Source 1 q',0.1,1.,dic['Source 1 q'][a1,a2])
    src = SBObjects.Sersic('Source 1',{'x':SX,'y':SY,'pa':SPA,'q':SQ,'re':SR,'n':SN,'c':2.})
    pars += [SX, SY, SR, SN, SPA, SQ]
    cov += [0.5,0.5,0.5,0.2,1.,0.05]
    srcs = [src]

    # source 2
    if 'Source 2 re' in dic.keys():
        spa = dic['Source 2 pa'][a1,a2]
        if name in ['J0837','J1323','J1347']:
            sx,sy = dic['Source 2 x'][a1,a2], dic['Source 2 y'][a1,a2]
        else:
            sx,sy = dic['Source 2 x'][a1,a2]+vx, dic['Source 2 y'][a1,a2]+vy
        print sx,sy
        SX2 = pymc.Uniform('Source 2 x',sx-20,sx+20,sx)
        SY2 = pymc.Uniform('Source 2 y',sy-20,sy+20,sy)
        SR2 = pymc.Uniform('Source 2 re',0.5,100,dic['Source 2 re'][a1,a2])
        SN2 = pymc.Uniform('Source 2 n',0.5,8,dic['Source 2 n'][a1,a2])
        SPA2 = pymc.Uniform('Source 2 pa',spa-90,spa+90,spa)
        SQ2 = pymc.Uniform('Source 2 q',0.2,1.,dic['Source 2 q'][a1,a2])
        src2 = SBObjects.Sersic('Source 2',{'x':SX2,'y':SY2,'pa':SPA2,'q':SQ2,'re':SR2,'n':SN2,'c':2.})
        pars += [SX2, SY2, SR2, SN2, SPA2, SQ2]
        cov += [0.5,0.5,0.5,0.2,1.,0.05]
        srcs.append(src2)

    print '####'
    if name == 'J1606':
        BOX = pymc.Uniform('boxiness',1,4,value=dic['boxiness'][a1,a2])
        pars += [BOX]
        cov += [0.1]

    log_Mstar = logM[names==name]
    Mstar_start = (10**log_Mstar) * 1e-12  *0.7
    Mstar = pymc.Uniform('stellar mass',0.0,10.,Mstar_start)
    pars += [Mstar]
    cov += [0.1]

    if oldResult is not None:
        print 'initialising based on previous chain'
        lp,trace,dic,_ = oldResult
        a1,a3 = np.unravel_index(lp[:,0].argmax(),lp[:,0].shape)
        for i in range(len(pars)):
            pars[i].value = trace[a1,0,a3,i]
        
    # with M/L deflections
    @pymc.deterministic
    def logP(value=0.,p=pars):
        lp = 0.
        for i in range(len(imgs)):
            dx,dy = dxs[i],dys[i]
            xp,yp = xc+dx,yc+dy
            image,sigma,psf = imgs[i],sigs[i],PSFs[i]
            imin,sigin,xin,yin = image[mask], sigma[mask],xp[mask2],yp[mask2]
            # model
            model = np.empty((len(srcs)+1,imin.size))
            for lens in lenses:
                lens.setPars()
            x0,y0 = pylens.getDeflections(lenses,[xin,yin])
            x0 -= Mstar.value*xstarms[i]
            y0 -= Mstar.value*ystarms[i]  # check I accounted for offsets?
            n=0
            for src in srcs:
                src.setPars()
                tmp = xc*0.
                if name == 'J1606' and src.name == 'Source 2':
                    tmp[mask2] = src.boxypixeval(x0,y0,1./OVRS,csub=31,c=BOX.value)
                else:
                    tmp[mask2] = src.pixeval(x0,y0,1./OVRS,csub=31)
                tmp = iT.resamp(tmp,OVRS,True)
                tmp = convolve.convolve(tmp,psf,False)[0]
                model[n] = tmp[mask].ravel()
                if name == 'J0837' and src.name == 'Source 2':
                    model[n] *= -1.
                n+=1
            model[n] = np.ones(model[n].size)
            # nnls
            rhs = (imin/sigin) # data
            op = (model/sigin).T # model matrix
            fit, chi = optimize.nnls(op,rhs)
            model = (model.T*fit).sum(1)
            resid = (model-imin)/sigin
            lp += -0.5*(resid**2.).sum()
        print lp
        return lp

    @pymc.observed
    def likelihood(value=0.,lp=logP):
        return lp

    cov = np.array(cov)
    print len(pars), len(cov)
    
    for i in range(len(pars)):
        print "%18s  %8.3f"%(pars[i].__name__,pars[i].value)

    '''S = myEmcee.PTEmcee(pars+[likelihood],cov=cov,nthreads=2,nwalkers=100,ntemps=6)#3)
    #if oldResult is not None:
    #    S.p0 = trace[-1]
    S.sample(500)

    print 'done emcee'
    outFile = '/data/ljo31b/EELs/galsub/emceeruns/'+name+'_parametric_eGNFW_0'

    f = open(outFile,'wb')
    cPickle.dump(S.result(),f,2)
    f.close()

    print 'cPickled!'

    result = S.result()
    lp,trace,dic,_ = result
    a1,a3 = np.unravel_index(lp[:,0].argmax(),lp[:,0].shape)
    a2=0
    for i in range(len(pars)):
        pars[i].value = trace[a1,0,a3,i]
        print "%18s  %8.3f"%(pars[i].__name__,pars[i].value)


    jj=0
    for jj in range(10):
        S.p0 = trace[-1]
        print 'sampling'
        S.sample(1000)

        f = open(outFile,'wb')
        cPickle.dump(S.result(),f,2)
        f.close()

        result = S.result()
        lp = result[0]

        trace = np.array(result[1])
        a1,a3 = np.unravel_index(lp[:,0].argmax(),lp[:,0].shape)
        for i in range(len(pars)):
            pars[i].value = trace[a1,0,a3,i]
            print "%18s  %8.3f"%(pars[i].__name__,pars[i].value)
        print jj
        jj+=1'''
    
    
    colours = ['V', 'I']
    models = []
    fits = []
    for i in range(len(imgs)):
        dx,dy = dxs[i],dys[i]
        xp,yp = xc+dx,yc+dy
        image,sigma,psf = imgs[i],sigs[i],PSFs[i]
        imin,sigin,xin,yin = image.flatten(), sigma.flatten(),xp.flatten(),yp.flatten()
        model = np.empty((len(srcs)+1,imin.size))
        for lens in lenses:
            lens.setPars()
        x0,y0 = pylens.getDeflections(lenses,[xin,yin])
        x0 -= Mstar.value*xstars[i]
        y0 -= Mstar.value*ystars[i] 
        
        n = 0
        for src in srcs:
            src.setPars()
            tmp = xc*0.
            if name == 'J1606' and src.name == 'Source 2':
                tmp = src.boxypixeval(x0,y0,1./OVRS,csub=31,c=BOX.value).reshape(xc.shape)
            else:
                tmp = src.pixeval(x0,y0,1./OVRS,csub=31).reshape(xc.shape)
            tmp = iT.resamp(tmp,OVRS,True)
            tmp = convolve.convolve(tmp,psf,False)[0]
            model[n] = tmp.ravel()
            if name == 'J0837' and src.name == 'Source 2':
                model[n] *= -1.
                print 'it is J0837'
            n+=1
        model[n] = np.ones(model[n].shape)
        n+=1
        rhs = image[mask]/sigma[mask]
        mmodel = model.reshape((n,image.shape[0],image.shape[1]))
        mmmodel = np.empty((n,image[mask].size))
        for m in range(mmodel.shape[0]):
            mmmodel[m] = mmodel[m][mask]
        op = (mmmodel/sigma[mask]).T
        rhs = image[mask]/sigma[mask]
        fit, chi = optimize.nnls(op,rhs)
        components = (model.T*fit).T.reshape((n,image.shape[0],image.shape[1]))
        model = components.sum(0)
        models.append(model)
        print fit
        NotPlicely(image,model,sigma,colours[i])
        pl.show()
        #for ii in range(3):
        #    pl.figure()
        #    pl.imshow(components[ii],interpolation='nearest',origin='lower')
        #    pl.colorbar()
        #pl.show()
                  
  


# make image smaller
name = 'J0913'
file = py.open('/data/ljo31b/EELs/galsub/images/'+name+'_maxlnL.fits')
img1,img2 = file[1].data, file[2].data

# make images 120 pixels across
XX=60.
my,mx = img1.shape[0]/2., img1.shape[1]/2.
_,sig1,psf1,_,sig2,psf2,DX,DY,_,mask = EasyAddImages(name)

img1,img2 = img1[my-XX:my+XX,mx-XX:mx+XX], img2[my-XX:my+XX,mx-XX:mx+XX]
sig1,sig2,mask = sig1[my-XX:my+XX,mx-XX:mx+XX],sig2[my-XX:my+XX,mx-XX:mx+XX],mask[my-XX:my+XX,mx-XX:mx+XX]

psf1 = psf1/np.sum(psf1)
psf2 = psf2/np.sum(psf2)

OVRS=1
if name in ['J1125','J1347']:
    OVRS=2
yc,xc = iT.overSample(img1.shape,OVRS)
yo,xo = iT.coords(img1.shape)
xc,yc = xc+DX+(mx-XX), yc+DY+(my-XX) 
xo,yo = xo+DX+(mx-XX), yo+DY+(my-XX) 
if OVRS >1:
    tck = RectBivariateSpline(yo[:,0],xo[0],mask)
    mask2 = tck.ev(xc,yc)
    mask2[mask2<0.5] = 0
    mask2[mask2>0.5] = 1
    mask2 = mask2==0
mask = mask==0
if OVRS == 1:
    mask2 = mask

# Start off DM with same q and pa as light. Try to make it as automated as possible
# stellar mass deflection angles
V,I,_ = np.load('/data/ljo31b/EELs/galsub/deflections/light_deflections_'+name+'.npy')
# make these the right shapes!
V = [V[ii][my-XX:my+XX,mx-XX:mx+XX] for ii in range(len(V))]
I = [I[ii][my-XX:my+XX,mx-XX:mx+XX] for ii in range(len(I))]
if OVRS>1:
    for i in range(len(V)):
        tck = RectBivariateSpline(yo[:,0],xo[0],V[i])
        V[i] = tck.ev(yc,xc)
        tck = RectBivariateSpline(yo[:,0],xo[0],I[i])
        I[i] = tck.ev(yc,xc)

xstars, ystars = [V[0].flatten(), I[0].flatten()], [V[1].flatten(), I[1].flatten()]
xstarms, ystarms = [V[0][mask2],I[0][mask2]], [V[1][mask2], I[1][mask2]]

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
    else:
        result = np.load(dir+name+'_211')

#if name == 'J0913':
#    print 'IMPORTANT: doing a 211 model for J0913 to see if this helps'
#    result = np.load('/data/ljo31/Lens/LensModels/J0913_212') #dir+name+'_211')

#oldResult = np.load('/data/ljo31b/EELs/galsub/emceeruns/'+name+'_parametric_DPL_1')
MakeModel(name,result)#,oldResult)

