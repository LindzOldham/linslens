import numpy as np, pylab as pl
from linslens.Plotter import *
from linslens.pixplot import *
import indexTricks as iT
from pylens import MassModels,pylens,adaptTools as aT,pixellatedTools as pT
from imageSim import SBModels,convolve
from scipy.sparse import diags
import pymc,cPickle
from linslens.GrabImages_huge import *
import glob
from scipy.ndimage.morphology import binary_erosion as erode

''' lets do these one at a time for now, to avoid errors '''
''' start off with power laws '''
''' pixellated sources '''
''' band by band for now? pixellated sources on many bands? '''

def MakeModel(name,result,oresult,Npnts=1):
    
    # read model
    lp,trace,dic,_ = result
    a2=0
    a1,a3 = np.unravel_index(lp[:,0].argmax(),lp[:,0].shape)

    olp,_,odic,_= oresult
    if len(olp.shape)==3.:
        olp = olp[:,0]
        for key in odic.keys():
            odic[key] = odic[key][:,0]
    oa1,oa2 = np.unravel_index(olp.argmax(),olp.shape)

    # now set up lenses
    LX = dic['Lens 1 x'][a1,a2,a3]
    LY = dic['Lens 1 y'][a1,a2,a3]
    LB = dic['Lens 1 b'][a1,a2,a3]
    LETA = dic['Lens 1 eta'][a1,a2,a3]
    LQ = dic['Lens 1 q'][a1,a2,a3]
    LPA = dic['Lens 1 pa'][a1,a2,a3]
    SH = dic['shear'][a1,a2,a3]
    SHPA = dic['shear pa'][a1,a2,a3]
    Mstar = dic['stellar mass'][a1,a2,a3]

    lens = MassModels.PowerLaw('lens',{'x':LX,'y':LY,'b':LB,'eta':LETA,'q':LQ,'pa':LPA})
    shear = MassModels.ExtShear('shear',{'x':LX,'y':LY,'b':SH,'pa':SHPA})
    lenses = [lens,shear]
    for l in lenses:
        l.setPars()

    if 'DELTAX' in dic.keys():
        print 'adding stellar mass offsets...'
        DELTAX = dic['DELTAX'][a1,0,a3]
        DELTAY = dic['DELTAY'][a1,0,a3]
    else:
        DELTAX,DELTAY = 0.,0.

    # now set up the images etc
    imgs = [img1,img2]
    sigs = [sig1,sig2]
    ifltms = [img[pix_mask] for img in imgs]
    sfltms = [sig[pix_mask] for sig in sigs]
    vfltms = [sfltm**2 for sfltm in sfltms]
    cmatms = [diags(1./sfltm,0) for sfltm in sfltms]
    xm,ym = x[pix_mask],y[pix_mask]
    coords = [[xm,ym],[xm+odic['xoffset'][oa1,oa2],ym+odic['yoffset'][oa1,oa2]]]

    # stellar mass deflection angles
    x_stars, y_stars = [V[0].flatten(), I[0].flatten()], [V[1].flatten(), I[1].flatten()]
    x_starms, y_starms = [V[0][pix_mask],I[0][pix_mask]], [V[1][pix_mask], I[1][pix_mask]]

    PSFs = [pT.getPSFMatrix(psf, img1.shape) for psf in [psf1,psf2]]
    PSFms = [pT.maskPSFMatrix(PSF,pix_mask) for PSF in PSFs]
    
    iflts = [img1.flatten(),img2.flatten()]
    sflts = [sig1.flatten(),sig2.flatten()]
    vflts = [sflt**2. for sflt in sflts]
    xflts = [x.flatten(), (x+odic['xoffset'][oa1,oa2]).flatten()]
    yflts = [y.flatten(), (y+odic['yoffset'][oa1,oa2]).flatten()]

    srcs = []
    for ii in range(len(iflts)-1):
        srcs.append(aT.AdaptiveSource(ifltms[ii]/sfltms[ii],ifltms[ii].size/Npnts))
        xl,yl = pylens.getDeflections(lenses,coords[ii])
        xl,yl = xl - x_starms[ii]*Mstar, yl - y_starms[ii]*Mstar
        srcs[ii].update(xl,yl)

    import time
    #reg=5.

    def doFit(p=None):
        global reg
        #reg=5.
        lp = 0.
        
        for l in lenses:
            l.setPars()
        for ii in range(len(ifltms)-1):
            src,ifltm,sfltm,vfltm,PSFm,cmatm = srcs[ii],ifltms[ii],sfltms[ii],vfltms[ii],PSFms[ii],cmatms[ii]
            PSF,coord,iflt,sflt = PSFs[ii],coords[ii],iflts[ii],sflts[ii]
            coord[0] += DELTAX
            coord[1]+=DELTAY

            x_star,y_star,x_starm,y_starm =x_stars[ii],y_stars[ii],x_starms[ii],y_starms[ii]
            xflt,yflt = xflts[ii],yflts[ii]
            xflt+=DELTAX
            yflt+=DELTAY

            xl,yl = pylens.getDeflections(lenses,coord)
            xl,yl = xl - x_starm*Mstar, yl - y_starm*Mstar

            src.update(xl,yl,doReg=True)
            lmat = PSFm*src.lmat
            rmat = src.rmat

            nupdate = 0
            res,fit,model,_,regg = aT.getModelG(ifltm,vfltm,lmat,cmatm,rmat,reg,nupdate)   
            reg = regg[0]
            print 'reg',reg
            xl,yl = pylens.getDeflections(lenses,[xflt,yflt])
            xl,yl = xl - x_star*Mstar, yl - y_star*Mstar
            oimg,pix = src.eval(xl,yl,fit,domask=False)
            oimg = PSF*oimg
            res = (iflt-oimg)/sflt
            lp+= -0.5*(res**2).sum()
        return lp
            
    # check initial model
    doFit()
    doFit()
    print 'current regularisation (set by hand): ', '%.1f'%reg
    
    for ii in range(len(imgs)-1):
        src,ifltm,sfltm,vfltm,PSFm,cmatm = srcs[ii],ifltms[ii],sfltms[ii],vfltms[ii],PSFms[ii],cmatms[ii]
        x_star,y_star,x_starm,y_starm,PSF =x_stars[ii],y_stars[ii],x_starms[ii],y_starms[ii],PSFs[ii]
        img,sig,coord = imgs[ii],sigs[ii],coords[ii]
        coord[0]+=DELTAX
        coord[1]+=DELTAY
        
        xl,yl = pylens.getDeflections(lenses,coord)
        xl,yl = xl - x_starm*Mstar, yl - y_starm*Mstar
        
        print xl.shape
        srcs[ii].update(xl,yl)
        osrc = showRes(xl,yl,srcs[ii],PSFm,img,sig,pix_mask,ifltm,vfltm,cmatm,reg,0,400)
        pl.show()
       

# make image smaller
name = 'J0901'
X = 0
file = py.open('/data/ljo31b/EELs/galsub/images/'+name+'_maxlnL.fits')
img1,img2 = file[1].data, file[2].data

# make images 120 pixels across
XX=30.
my,mx = img1.shape[0]/2., img1.shape[1]/2.
_,sig1,psf1,_,sig2,psf2,DX,DY,_,_ = EasyAddImages(name)

img1,img2 = img1[my-XX:my+XX,mx-XX:mx+XX], img2[my-XX:my+XX,mx-XX:mx+XX]
sig1,sig2 = sig1[my-XX:my+XX,mx-XX:mx+XX],sig2[my-XX:my+XX,mx-XX:mx+XX]

psf1 = psf1/np.sum(psf1)
psf2 = psf2/np.sum(psf2)

y,x = iT.coords(img1.shape)
x,y = x+DX+(mx-XX), y+DY+(my-XX) 

pix_mask = py.open('/data/ljo31b/EELs/galsub/masks/'+name+'.fits')[0].data.copy()[my-XX:my+XX,mx-XX:mx+XX]
pix_mask = pix_mask==1
pix_mask=erode(pix_mask,iterations=3)

# stellar mass deflection angles
V,I,_ = np.load('/data/ljo31b/EELs/galsub/deflections/light_deflections_'+name+'.npy')
# make these the right shapes!
V = [V[ii][my-XX:my+XX,mx-XX:mx+XX] for ii in range(len(V))]
I = [I[ii][my-XX:my+XX,mx-XX:mx+XX] for ii in range(len(I))]



# just need to grab a couple of things
dir = '/data/ljo31/Lens/LensModels/twoband/'
try:
    oresult = np.load(dir+name+'_212')
except:
    if name == 'J1347':
        oresult = np.load(dir+name+'_112')
    elif name == 'J1619':
        oresult = np.load(dir+name+'_212')
    else:
        oresult = np.load(dir+name+'_211')

reg=0.5
# load results for making model
files = glob.glob('/data/ljo31b/EELs/galsub/emceeruns/'+name+'_pixellated_*')
files.sort()
result = np.load(files[-3])
print 'reading in file ', files[-3]
lp,trace,dic,_ = result

MakeModel(name,result,oresult)

pl.figure()
pl.plot(lp[:,0])
#for key in dic.keys():
#    pl.figure()
#    pl.hist(dic[key][:,0].ravel(),30)
#    pl.title(key)

pl.show()
