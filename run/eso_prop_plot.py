from linslens import EELsModels_huge as L
import numpy as np
import pylab as pl
from linslens import EELsKeckModels as K
from linslens.Plotter import *
import colorImage
from scipy import ndimage
from jeans.makemodel import deproject
import lenslib
import glob
from linslens.Profiles import *
from astLib import astCalc
from linslens.GetStellarMass import *
from scipy.interpolate import splrep, splev
import matplotlib

matplotlib.rc('axes',edgecolor='white')

def clip(arr,nsig=3.5):
    a = arr.flatten()
    m,s,l = a.mean(),a.std(),a.size
    while 1:
        a = a[abs(a-m)<s*nsig]
        if a.size==l:
            return m,s
        m,s,l = a.mean(),a.std(),a.size


def MakeTab(model,kmodel,name):
    model.Initialise()
    imgs,sigs,mods = model.imgs, model.sigs,model.models
    xm,ym = int(imgs[0].shape[1]/2.),int(imgs[0].shape[0]/2.)
    if name == 'J1606':
        ym -= 10
    elif name == 'J0837':
        ym -= 5
    for i in range(len(imgs)):
        imgs[i] = imgs[i][ym-30:ym+30,xm-30:xm+30]
        sigs[i] = sigs[i][ym-30:ym+30,xm-30:xm+30]
        mods[i] = mods[i][ym-30:ym+30,xm-30:xm+30]

    xx,yy = np.arange(5,26,1), np.ones(21)*5.
    if name == 'J1125':
        pl.subplot(231)
    elif name == 'J1347':
        pl.subplot(232)
    elif name == 'J1144':
        pl.subplot(233)
    elif name == 'J1619':
        pl.subplot(234)
    elif name == 'J1606':
        pl.subplot(235)
    elif name == 'J1446':
        pl.subplot(236)
    
    ###
    kmodel.Initialise()
    img,sig,mod = kmodel.img,kmodel.sig,kmodel.model
    if name in ['J1218','J1347','J1606','J2228']:
        scale = 5./3
    elif name in ['J0837','J0901','J0913','J1125','J1144','J1323','J1446','J1605','J1619']:
        scale = 5.
    xm,ym = int(img.shape[1]/2.), int(img.shape[0]/2.)
    if name == 'J1606':
        ym -= 10*scale
    elif name == 'J0837':
        ym -= 5.*scale
    img = img[ym-30*scale:ym+30*scale,xm-30*scale:xm+30*scale]
    sig = sig[ym-30*scale:ym+30*scale,xm-30*scale:xm+30*scale]
    mod = mod[ym-30*scale:ym+30*scale,xm-30*scale:xm+30*scale]
    ############pl.subplot(144)
    ###
    
    CI = colorImage.ColorImage()
    CI.nonlin = 10.
    V,I,K = imgs[0],imgs[1],img
    
    dxi,dyi = model.Ddic['xoffset'],model.Ddic['yoffset']
    dxk,dyk = kmodel.Ddic['xoffset']+kmodel.Dx,kmodel.Ddic['yoffset']+kmodel.Dy
    I = ndimage.shift(I,[dyi,dxi])
    
    K = ndimage.zoom(K,1./scale)
    K = ndimage.shift(K,[-shifts[name][1],-shifts[name][0]])

    vki = CI.createModel(V,I,K)
    pl.imshow(vki,interpolation='nearest',origin='lower')
    if name == 'J1619':
        pl.plot(xx,yy,color='White',lw=4)
        pl.imshow(vki,interpolation='nearest',origin='lower')  
    #pl.figtext(0.14,0.1,"$1''$",color='White',fontsize=45,weight=1000,family='sans-serif',stretch='ultra-expanded') 
    pl.figtext(0.185,0.15,"$1''$",color='White',fontsize=45,weight=1000,family='sans-serif',stretch='ultra-expanded') 

    pl.gca().xaxis.set_ticks([])
    pl.gca().yaxis.set_ticks([])
    #pl.spines['all'].set_color('white')

 

dir = '/data/ljo31/Lens/LensModels/twoband/'
sz = np.load('/data/ljo31/Lens/LensParams/SourceRedshiftsUpdated.npy')[()]
names = sz.keys()
names.sort()
shifts = [['J0837',[1.3,2.9]],['J0901',[0,3]],['J0913',[0,0]],['J1125',[0.9,2.1]],['J1144',[-0.8,4.55]],['J1218',[-3.2,-6.8]],['J1323',[1.0,1.1]],['J1347',[3.7,0]],['J1446',[0,0]],['J1605',[-3.5,0.1]],['J1606',[-4.8,7.6]],['J1619',[3.8,0]],['J2228',[-2.4,-3.6]]]
        
shifts = dict(shifts)
lz = np.load('/data/ljo31/Lens/LensParams/LensRedshiftsUpdated.npy')[()]

fig = pl.figure(tight_layout=False,figsize=(12.4,8))

'''ax = [fig.add_subplot(2,2,i+1) for i in range(4)]

for a in ax:
    a.set_xticklabels([])
    a.set_yticklabels([])
    a.set_aspect('equal')

fig.subplots_adjust(wspace=0, hspace=0)
i=0'''
#fig.subplots_adjust(wspace=0.1, hspace=0)

for name in ['J1125','J1619','J1606','J1446','J1347','J1144']:
    #axis = ax[i]
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
    
    model = L.EELs(result, name)
    kmodel = K.EELs(kresult,result,name)
    MakeTab(model,kmodel,name)
    #i+=1

fig.subplots_adjust(hspace=0.0,wspace=0.0)
pl.savefig('/data/ljo31b/proposals/eels/EELS.pdf')
pl.show()
