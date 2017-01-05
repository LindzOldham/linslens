import cPickle,numpy,pyfits as py
import pymc
from pylens import *
from imageSim import SBModels,convolve,SBObjects
import indexTricks as iT
from SampleOpt import AMAOpt
import pylab as pl, numpy as np
import myEmcee_blobs as myEmcee
from scipy import optimize
from scipy.interpolate import RectBivariateSpline
from linslens import EELsKeckLensModels as K, EELsImages_huge as Image, EELsModels_huge as L

# grab image
name = 'J2228'
#X = 3
file = py.open('/data/ljo31b/EELs/galsub/images/'+name+'.fits')
img1,img2 = file[1].data, file[2].data

# make images 50 pixels across
my,mx = img1.shape[0]/2., img1.shape[1]/2.
img1,img2 = img1[my-60:my+60,mx-60:mx+60], img2[my-60:my+60,mx-60:mx+60]
_,sig1,psf1,_,sig2,psf2,DX,DY,OVRS,_ = Image.J2228()  
y,x = iT.coords(img1.shape)
x,y = x+DX+(mx-60.), y+DY+(my-60.) 

# Start off DM with same q and pa as light. Try to make it as automated as possible
# stellar mass deflection angles
V,_,_ = np.load('/data/ljo31b/EELs/galsub/deflections/light_deflections_'+name+'.npy')
# make these the right shapes!
V = [V[ii][my-60:my+60,mx-60:mx+60] for ii in range(len(V))]
xv,yv = V

masses = np.load('/data/ljo31b/EELs/inference/new/huge/masses_212.npy')
names = ['J0837','J0901','J0913','J1125','J1144','J1218','J1323','J1347','J1446','J1605','J1606','J1619','J2228']
logM = masses[0]

#pl.figure()
#pl.imshow(img1,interpolation='nearest',origin='lower')
#pl.show()

# write down the brightest pixels

# these are the indices
points = np.array([[75,60.],[51,69]])


xpoints,ypoints = points.T
xvn = [xv[ypoints[ii],xpoints[ii]] for ii in range(len(points))]
yvn = [yv[ypoints[ii],xpoints[ii]] for ii in range(len(points))]

#print xvn, xvn
#pl.figure()
#pl.imshow(xv,interpolation='nearest',origin='lower')
#pl.scatter(xpoints,ypoints,color='k',s=40)
#pl.show()


# write down initial model
dir = '/data/ljo31/Lens/LensModels/twoband/'

result = np.load(dir+name+'_211')
kresult = np.load(dir+name+'_Kp_211')
model = K.EELs(kresult,result,name)
model.Initialise()
b_start, sh_start, sh_pa_start = model.lenses[0].b, model.lenses[1].b, model.lenses[1].pa
gal = model.gals[0]
q_start, pa_start, x_start, y_start = gal.q, gal.pa, gal.x, gal.y
log_Mstar = logM[names==name]
Mstar_start = (10**log_Mstar) * 1e-12 * 0.75
print 'mstar', Mstar_start

print x_start,y_start,b_start,pa_start,q_start
LX = x_start#pymc.Uniform('Lens x',x_start-10,x_start+10,x_start)
LY = y_start#pymc.Uniform('Lens y',y_start-10,y_start+10,y_start)
LB = pymc.Uniform('Lens b',0.,b_start*2,b_start) # assuming half is light, half dark for now
#LETA = pymc.Uniform('Lens eta',0.5,1.5,1.)
LQ = q_start#pymc.Uniform('Lens q',0.1,1.,q_start)
LPA = pa_start# pymc.Uniform('Lens pa',-180,180,pa_start)
SH = sh_start#pymc.Uniform('shear',-0.3,0.3,sh_start)
SHPA = sh_pa_start#pymc.Uniform('shear pa',-180,180,sh_pa_start)
Mstar = pymc.Uniform('stellar mass',0.01,100.,Mstar_start)

lens = MassModels.PowerLaw('lens',{'x':LX,'y':LY,'b':LB,'eta':1,'q':LQ,'pa':LPA})
shear = MassModels.ExtShear('shear',{'x':LX,'y':LY,'b':SH,'pa':SHPA})
lenses = [lens,shear]

pars = [LB,Mstar]
cov = np.array([1.,2.])

@pymc.deterministic
def logP(value=0.,p=pars):
    lp = 0.
    for lens in lenses:
        lens.setPars()
    xl,yl = pylens.getDeflections(lenses,[xpoints,ypoints])
    xl,yl = xl-Mstar.value*xvn, yl-Mstar.value*yvn
    mx,my = np.mean(xl),np.mean(yl)
    rad = (xl-mx)**2. + (yl-my)**2.
    lp = -1*rad.sum()
    return lp

@pymc.observed
def likelihood(value=0.,lp=logP):
    return lp


SS = AMAOpt(pars,[likelihood],[logP],cov=cov)
SS.sample(20000)
lp,trace,det = SS.result()
pl.figure()
pl.plot(lp)
pl.show() 
print 'results from optimisation:'
for i in range(len(pars)):
    pars[i].value = trace[-1,i]
    print "%18s  %8.3f"%(pars[i].__name__,pars[i].value)
    

#for key in det.keys():
#    pl.figure()
#    pl.plot(det[key])
#    pl.title(key)

pl.show()

for lens in lenses:
    lens.setPars()
xl,yl = pylens.getDeflections(lenses,[xpoints,ypoints])
xl,yl = xl-Mstar.value*xvn, yl-Mstar.value*yvn
pl.figure()
pl.scatter(xl,yl)
pl.show()

mx,my = np.mean(xl),np.mean(yl)
rad = (xl-mx)**2. + (yl-my)**2.
print rad**0.5

# this seems to be a stable solution! Now let's put it in the gui/emcee!
