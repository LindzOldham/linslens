import pyfits as py, numpy as np
from imageSim import SBModels,convolve
from pylens import *
import indexTricks as iT
from scipy import optimize
from scipy.interpolate import splrep, splev, splint,RectBivariateSpline
from astLib import astCalc
import EELsImages_huge as EELsImages, Plotter
from GetPhotometry import *
import SBBModels, SBBProfiles
import itertools
from Plotter import *

def EasyAddImages(name):
    if name == 'J0837':
        img1,sig1,psf1,img2,sig2,psf2,Dx,Dy,OVRS,mask = EELsImages.J0837()
    elif name == 'J0901':
        img1,sig1,psf1,img2,sig2,psf2,Dx,Dy,OVRS,mask  = EELsImages.J0901()
    elif name == 'J0913':
        img1,sig1,psf1,img2,sig2,psf2,Dx,Dy,OVRS,mask  = EELsImages.J0913(2)
    elif name == 'J1125':
        img1,sig1,psf1,img2,sig2,psf2,Dx,Dy,OVRS,mask  = EELsImages.J1125()
    elif name == 'J1144':
        img1,sig1,psf1,img2,sig2,psf2,Dx,Dy,OVRS,mask   = EELsImages.J1144()
    elif name == 'J1218':
        img1,sig1,psf1,img2,sig2,psf2,Dx,Dy,OVRS,mask  = EELsImages.J1218()
    elif name == 'J1248':
        img1,sig1,psf1,Dx,Dy,OVRS,mask  = EELsImages.J1248()
    elif name == 'J1323':
        img1,sig1,psf1,img2,sig2,psf2,Dx,Dy,OVRS,mask  = EELsImages.J1323(2)
    elif name == 'J1347':
        img1,sig1,psf1,img2,sig2,psf2,Dx,Dy,OVRS,mask = EELsImages.J1347()
    elif name == 'J1446':
        img1,sig1,psf1,img2,sig2,psf2,Dx,Dy,OVRS,mask = EELsImages.J1446()
    elif name == 'J1605':
        img1,sig1,psf1,img2,sig2,psf2,Dx,Dy,OVRS,mask = EELsImages.J1605(2)
    elif name == 'J1606':
        img1,sig1,psf1,img2,sig2,psf2,Dx,Dy,OVRS,mask = EELsImages.J1606()
    elif name == 'J1619':
        img1,sig1,psf1,img2,sig2,psf2,Dx,Dy,OVRS,mask = EELsImages.J1619()
    elif name == 'J2228':
        img1,sig1,psf1,img2,sig2,psf2,Dx,Dy,OVRS,mask = EELsImages.J2228()
    else:
        print 'are you sure this is an EEL?'
        return
    return img1,sig1,psf1,img2,sig2,psf2,Dx,Dy,OVRS,mask
