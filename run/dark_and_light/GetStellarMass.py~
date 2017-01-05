from linslens.Profiles import *
from astLib import astCalc
import glob, numpy as np, pylab as pl
import lenslib
from jeans.makemodel import deproject
from imageSim import SBModels

def GetStellarMass(oresult,Mstar,r,scale,name,deprojected=False):
    fracs= np.load('/data/ljo31b/EELs/galsub/fracs.npy')[()]
    olp,otrace,odic,_ = oresult
    if len(olp.shape)==3:
        olp = olp[:,0]
        for key in odic.keys():
            odic[key] = odic[key][:,0]
    oa1,oa2 = np.unravel_index(olp.argmax(),olp.shape)

    if 'Galaxy 1 x' in odic.keys():
        gx1,gy1 = odic['Galaxy 1 x'][oa1,oa2]*0.05*scale, odic['Galaxy 1 y'][oa1,oa2]*0.05*scale
    else:
        gx1,gy1 = odic['Galaxy 2 x'][oa1,oa2]*0.05*scale, odic['Galaxy 2 y'][oa1,oa2]*0.05*scale
    gr1,gn1 = odic['Galaxy 1 re'][oa1,oa2]*0.05*scale, odic['Galaxy 1 n'][oa1,oa2]
    gpa1,gq1 = odic['Galaxy 1 pa'][oa1,oa2], odic['Galaxy 1 q'][oa1,oa2]
    frac = [1.,0.]
    gal1 = SBModels.Sersic('galaxy 1',{'x':0,'y':0,'re':gr1,'n':gn1,'pa':gpa1,'q':gq1})
    gals = [gal1]
    if 'Galaxy 2 re' in odic.keys():
        if 'Galaxy 2 x' not in odic.keys():
            gx2,gy2 = gx1, gy1
        else:
            gx2,gy2 = odic['Galaxy 2 x'][oa1,oa2]*0.05*scale, odic['Galaxy 2 y'][oa1,oa2]*0.05*scale
        gr2,gn2 = odic['Galaxy 2 re'][oa1,oa2]*0.05*scale, odic['Galaxy 2 n'][oa1,oa2]
        try:
            gpa2,gq2 = odic['Galaxy 2 pa'][oa1,oa2], odic['Galaxy 2 q'][oa1,oa2]
        except:
            gpa2,gq2 = odic['Galaxy 1 pa'][oa1,oa2], odic['Galaxy 2 q'][oa1,oa2]
        frac = fracs[name]
        gal2 = SBModels.Sersic('galaxy 2',{'x':gx2-gx1,'y':gy2-gy1,'re':gr2,'n':gn2,'pa':gpa2,'q':gq2})
        gals.append(gal2) 

    if frac[1]>frac[0]:
        # centroid on the brightest component
        gals[0].x += (gx1-gx2)
        gals[0].y += (gy1-gy2)
        gals[1].x,gals[1].y = 0.,0.
    sb = np.zeros(r.size)

    for ii in range(len(gals)):
        sb += frac[ii]*gals[ii].eval(r)

    #deproject
    lr,light = deproject(r,sb)
    if deprojected:
        return lr,light
    # cumulatively sum
    lightmod = splrep(lr,light*4.*np.pi*lr**2.)
    cumlight = np.array([splint(0,lr[ii],lightmod) for ii in range(lr.size)])
    cumlight /= cumlight[-1]
    LM = cumlight*Mstar*1e12
    if name == 'J2228':
        LM *= 5. # different deflection angle normalisation!
    LM[0] = LM[1]
    return LM
