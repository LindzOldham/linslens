import pylab as pl, pyfits as py, numpy as np

def NotPlicely(image,im,sigma,colour,cmap='afmhot'):
    ext = [0,image.shape[0],0,image.shape[1]]
    #vmin,vmax = numpy.amin(image), numpy.amax(image)
    pl.figure()
    pl.subplot(221)
    pl.imshow(image,origin='lower',interpolation='nearest',extent=ext,cmap=cmap,aspect='auto',vmin=0,vmax=np.amax(image)*0.5) #,vmin=vmin,vmax=vmax)
    pl.colorbar()
    pl.title('data')
    pl.subplot(222)
    pl.imshow(im,origin='lower',interpolation='nearest',extent=ext,cmap=cmap,aspect='auto',vmin=0,vmax=np.amax(image)*0.5) #,vmin=vmin,vmax=vmax)
    pl.colorbar()
    pl.title('model')
    pl.subplot(223)
    pl.imshow(image-im,origin='lower',interpolation='nearest',extent=ext,vmin=-0.25,vmax=0.25,cmap=cmap,aspect='auto')
    pl.colorbar()
    pl.title('data-model')
    pl.subplot(224)
    pl.imshow((image-im)/sigma,origin='lower',interpolation='nearest',extent=ext,vmin=-5,vmax=5,cmap=cmap,aspect='auto')
    pl.title('signal-to-noise residuals')
    pl.colorbar()
    pl.suptitle(colour)

def SotPleparately(image,im,sigma,col,vmin=None,vmax=None,cmap='afmhot'):
    ext = [0,image.shape[0],0,image.shape[1]]
    if vmin is None:
        vmin = -5
    if vmax is None:
        vmax = 5
    pl.figure()
    pl.imshow(image,origin='lower',interpolation='nearest',extent=ext,cmap=cmap,vmin=0,aspect='auto') #,vmin=vmin,vmax=vmax)
    #pl.colorbar()
    #pl.title('data - '+str(col))
    pl.gca().xaxis.set_major_locator(pl.NullLocator())
    pl.gca().yaxis.set_major_locator(pl.NullLocator())

    pl.figure()
    pl.imshow(im,origin='lower',interpolation='nearest',extent=ext,cmap=cmap,vmin=0,aspect='auto') #,vmin=vmin,vmax=vmax)
    #pl.colorbar()
    #pl.title('model - '+str(col))
    pl.gca().xaxis.set_major_locator(pl.NullLocator())
    pl.gca().yaxis.set_major_locator(pl.NullLocator())
    pl.figure()
    pl.imshow((image-im)/sigma,origin='lower',interpolation='nearest',extent=ext,vmin=vmin,vmax=vmax,cmap=cmap,aspect='auto')
    #pl.title('signal-to-noise residuals - '+str(col))
    #pl.colorbar()
    pl.gca().xaxis.set_major_locator(pl.NullLocator())
    pl.gca().yaxis.set_major_locator(pl.NullLocator())

# this has to be edited depending on number of components...
def CotSomponents(components,col):
    pl.figure()
    pl.subplot(221)
    pl.imshow(components[0],interpolation='nearest',origin='lower',cmap='afmhot',aspect='auto',vmax=np.amax(components[0])*0.7)
    pl.colorbar()
    pl.title('galaxy 1 ')
    pl.subplot(222)
    pl.imshow(components[1],interpolation='nearest',origin='lower',cmap='afmhot',aspect='auto',vmax=np.amax(components[1])*0.7)
    pl.colorbar()
    pl.title('galaxy 2 ')
    pl.subplot(223)
    pl.imshow(components[2],interpolation='nearest',origin='lower',cmap='afmhot',aspect='auto',vmax=np.amax(components[2])*0.7)
    pl.colorbar()
    pl.title('source 1 ')
    pl.suptitle(col)
    pl.subplot(224)
    pl.imshow(components[3],interpolation='nearest',origin='lower',cmap='afmhot',aspect='auto',vmax=np.amax(components.sum(0))*0.7)
    pl.colorbar()
    pl.title('model ')
    pl.suptitle(col)

def SotPleparatelyPaper(image,im,sigma,vmin=None,vmax=None):
    ext = [0,image.shape[0],0,image.shape[1]]
    if vmin is None:
        vmin = -5
    if vmax is None:
        vmax = 5
    pl.figure()
    pl.imshow(image,origin='lower',interpolation='nearest',extent=ext,cmap='afmhot_r',vmin=0,vmax=0.75*np.amax(image))
    #pl.colorbar()
    pl.figure()
    pl.imshow(im,origin='lower',interpolation='nearest',extent=ext,cmap='afmhot_r',vmin=0,vmax=0.75*np.amax(image))
    #pl.colorbar()
    pl.figure()
    pl.imshow((image-im)/sigma,origin='lower',interpolation='nearest',extent=ext,vmin=vmin,vmax=vmax,cmap='afmhot')
    #pl.colorbar()
