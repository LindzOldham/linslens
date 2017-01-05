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
import VDfit
from math import sqrt,log,log10
import scipy
nfit = 6
bg = 'polynomial'
import special_functions as sf
bias=100.
import numpy
light = 299792.458
from scipy import optimize

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
    if name == 'J0837':
        ax = ax1
    elif name == 'J0913':
        ax = ax2
    elif name == 'J1446':
        ax = ax3
    elif name == 'J1606':
        ax = ax4
    
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
    ax.imshow(vki,interpolation='nearest',origin='lower')
    if name == 'J1446':
        ax.plot(xx,yy,color='White',lw=4)
        ax.imshow(vki,interpolation='nearest',origin='lower')  
    pl.figtext(0.05,0.15,"$1''$",color='White',fontsize=45,weight=1000,family='sans-serif',stretch='ultra-expanded')   

    #ax.gca().xaxis.set_ticks([])
    #ax.gca().yaxis.set_ticks([])
    
    ax.set_xticklabels([])
    ax.set_yticklabels([])
 

dir = '/data/ljo31/Lens/LensModels/twoband/'
sz = np.load('/data/ljo31/Lens/LensParams/SourceRedshiftsUpdated.npy')[()]
names = sz.keys()
names.sort()
shifts = [['J0837',[1.3,2.9]],['J0901',[0,3]],['J0913',[0,0]],['J1125',[0.9,2.1]],['J1144',[-0.8,4.55]],['J1218',[-3.2,-6.8]],['J1323',[1.0,1.1]],['J1347',[3.7,0]],['J1446',[0,0]],['J1605',[-3.5,0.1]],['J1606',[-4.8,7.6]],['J1619',[3.8,0]],['J2228',[-2.4,-3.6]]]
        
shifts = dict(shifts)
lz = np.load('/data/ljo31/Lens/LensParams/LensRedshiftsUpdated.npy')[()]

fig = pl.figure(tight_layout=True,figsize=(24,10))

ax1 = pl.subplot2grid((2,5),(0,0))
ax2 = pl.subplot2grid((2,5),(0,1))
ax3 = pl.subplot2grid((2,5),(1,0))
ax4 = pl.subplot2grid((2,5),(1,1))
ax5 = pl.subplot2grid((2,5),(0,2),colspan=3,rowspan=2) # spectrum goes here


'''ax = [fig.add_subplot(2,2,i+1) for i in range(4)]

for a in ax:
    a.set_xticklabels([])
    a.set_yticklabels([])
    a.set_aspect('equal')

fig.subplots_adjust(wspace=0, hspace=0)
i=0'''

for name in ['J0837','J0913','J1606','J1446']:
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

#fig.subplots_adjust(left=0.15,right=0.95)

# now make spectrum

ax = ax5
ax.set_xlabel(r'wavelength ($\mathrm{\AA}$)')
ax.set_ylabel('flux')

result = np.load('/data/ljo31b/EELs/esi/kinematics/inference/vdfit/NEW/J0837_0.31_source_esi_indous_vdfit')
mask = np.array([[7580,7700],[5570,5585],[6860,6900]])
sigsci = lambda wave: 20.40

t1 = VDfit.INDOUS(sigsci)
t2 = VDfit.PICKLES(sigsci)

sz = np.load('/data/ljo31/Lens/LensParams/SourceRedshiftsUpdated.npy')[()]
lz = np.load('/data/ljo31/Lens/LensParams/LensRedshiftsUpdated.npy')[()]
name = 'J0837'
wid = '0.31'

scispec = py.open('/data/ljo31b/EELs/esi/kinematics/apertures/final/'+name+'_ap_'+str(wid)+'_spec_source.fits')[0].data
varspec = py.open('/data/ljo31b/EELs/esi/kinematics/apertures/final/'+name+'_ap_'+str(wid)+'_var_source.fits')[0].data
sciwave = py.open('/data/ljo31b/EELs/esi/kinematics/apertures/final/'+name+'_ap_'+str(wid)+'_wl_source.fits')[0].data

zl,zs = lz[name][0], sz[name][0]
lim = 3800.*(1.+zl)
maxlim = 4800.*(1.+zs)
srclim = 3650.*(1.+zs)+40.

if srclim<=lim:
    srclim = lim+10.

edges = np.where(np.isnan(scispec)==False)[0]
start = edges[0]
end = np.where(sciwave>np.log10(maxlim))[0][0]
scispec = scispec[start:end]
varspec = varspec[start:end]
datascale = sciwave[1]-sciwave[0] 
sciwave = 10**sciwave[start:end]

# convert esi air wavelengths to vacuum wavelengths
# sdss website, cf, morton 1991
vac = np.linspace(3800,11000,(11000-3800+1)*5)
air = vac / (1. + 2.735182e-4 + (131.4182/vac**2.) + (2.76249e8/vac**8))
model = splrep(air,vac)
sciwave = splev(sciwave,model)
outwave = np.log10(sciwave)

zp = scispec.mean()
scispec /= zp
varspec /= zp**2

### reading routine
outwave,scispec,varspec = outwave[outwave>log10(lim)], scispec[outwave>log10(lim)], varspec[outwave>log10(lim)]
origwave,origsci,origvar = outwave.copy(),scispec.copy(),varspec.copy()
ma = np.ones(outwave.size)
for M in mask:
    cond = np.where((outwave>np.log10(M[0]))&(outwave<np.log10(M[1])))
    ma[cond]=0
ma=ma==1
outwave,scispec,varspec=outwave[ma],scispec[ma],varspec[ma]
isig = 1./varspec**0.5
ntemps1,ntemps2 = t1.nspex, t2.nspex

# Create the polynomial fit components
BIAS = scispec*0.
operator = scipy.zeros((scispec.size,2*ntemps1+ntemps2+nfit))
X = np.arange(outwave.size)
X = 2.*X/X.max() - 1.

for i in range(nfit):
    p = scipy.zeros((nfit,1))
    p[i] = 1.
    coeff = {'coeff':p,'type':bg}
    poly = sf.genfunc(X,0.,coeff)
    operator[:,i+2*ntemps1+ntemps2] = poly
    BIAS += bias*poly

oper = operator.T 
cond = np.where(outwave<=np.log10(srclim),True,False)


lp,trace,dic,_=result
a1,a2 = numpy.unravel_index(lp.argmax(),lp.shape)
velL,sigL,velS,sigS = trace[a1,a2]

zL, zS = zl+velL/light, zs+velS/light
linwave = 10**outwave
oper[:ntemps1] = t1.getSpectra(outwave,zL,sigL)
oper[ntemps1:2.*ntemps1,~cond] = t1.getSpectra(outwave[~cond],zS,sigS)
oper[2*ntemps1:2*ntemps1+ntemps2,cond] = t2.getSpectra(linwave[cond],zS,sigS)
    
op = (oper*isig).T
rhs = (scispec+BIAS)*isig
fit,chi = optimize.nnls(op,rhs)
for i in range(nfit):
    fit[ntemps1*2+ntemps2+i] -= bias
maskmodel = scipy.dot(oper.T,fit)
# unmasked
if mask is not None or restmask is not None:
    operator = scipy.zeros((origsci.size,2*ntemps1+ntemps2+nfit))
    X = np.arange(origwave.size)
    X = 2.*X/X.max() - 1.

    for i in range(nfit):
        p = scipy.zeros((nfit,1))
        p[i] = 1.
        coeff = {'coeff':p,'type':bg}
        poly = sf.genfunc(X,0.,coeff)
        operator[:,i+2*ntemps1+ntemps2] = poly

    oper = operator.T 
    origcond = np.where(origwave<=np.log10(srclim),True,False)
    linorigwave = 10**origwave

    oper[:ntemps1] = t1.getSpectra(origwave,zL,sigL)
    oper[ntemps1:2.*ntemps1,~origcond] = t1.getSpectra(origwave[~origcond],zS,sigS)
    oper[2*ntemps1:2*ntemps1+ntemps2,origcond] = t2.getSpectra(linorigwave[origcond],zS,sigS)

      
outmodel = scipy.dot(oper.T,fit)
lens = scipy.dot(oper[:ntemps1].T,fit[:ntemps1])
source = scipy.dot(oper[2*ntemps1:2*ntemps1+ntemps2].T,fit[2*ntemps1:2*ntemps1+ntemps2]) + scipy.dot(oper[ntemps1:2*ntemps1].T,fit[ntemps1:2*ntemps1])
cont = scipy.dot(oper[2*ntemps1+ntemps2:].T,fit[2*ntemps1+ntemps2:])

for M in mask:
    ax.axvspan(M[0], M[1], color='DarkGray')
ax.plot(10**origwave,origsci,'LightGray',)
ax.plot(10**origwave,outmodel,'k',)
ax.plot(10**origwave,lens,'SteelBlue',label='lens')
ax.plot(10**origwave,source,'Crimson',label='source')
ax.plot(10**origwave,cont,'Navy')
ax.legend(loc='upper left',frameon=False,bbox_to_anchor=(0.,1),fontsize=60)
ax.set_xlim(lim,maxlim)
ax.set_ylim(-0.5,4)

pl.savefig('/data/ljo31b/proposals/eels/EELS.pdf')
pl.show()


'''

13 47 04.96 
-01 01 03.57



'''
