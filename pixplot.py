import numpy,pyfits,pylab
import indexTricks as iT
from pylens import MassModels,pylens,adaptTools as aT,pixellatedTools as pT
from pylens import lensModel


# Function to make a `nice' plot
def showRes(x,y,src,psf,img,sig,mask,iflt,vflt,cmat,reg,niter,npix,vmax_src=2):
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
    pylab.imshow(img,origin='lower',interpolation='nearest',extent=ext,vmin=0,vmax=3,cmap='jet',aspect='auto')
    pylab.colorbar()
    pylab.subplot(222)
    pylab.imshow(oimg,origin='lower',interpolation='nearest',extent=ext,vmin=0,vmax=3,cmap='jet',aspect='auto')
    pylab.colorbar()
    pylab.subplot(223)
    pylab.imshow((img-oimg)/sig,origin='lower',interpolation='nearest',extent=ext,vmin=-5,vmax=5,cmap='jet',aspect='auto')
    pylab.colorbar()
    pylab.subplot(224)
    pylab.imshow(osrc,origin='lower',interpolation='nearest',extent=ext2,vmin=0,vmax=vmax_src,cmap='jet',aspect='auto')
    pylab.colorbar()
    return osrc
