from linslens import EELsModels_huge as L
import numpy as np, pylab as pl, pyfits as py, cPickle

''' this time, doing every tenth entry rather than every 400th as I think this is fairer. It will just take AN AGE '''

def MakeTab(file,name):
    result = np.load(file)
    model = L.EELs(result,name)
    model.Initialise()
    mags = model.GetIntrinsicMags()
    Res = model.GetSourceSize(kpc=True)
    restmags, Ls = model.GetPhotometry()
    mus = model.GetSB()
    restmus = model.GetRestSB()
    model.MakePDFDict()
    model.GetPDFs(kpc=True)
    med,lo,hi = model.UncertaintiesFromPDF(makecat=True)
    np.save('/data/ljo31/Lens/PDFs/211_huge_new_new_'+str(model.name), [np.array(model.muPDF),np.array(model.magPDF),np.array(model.RePDF)])
    return med,lo,hi

dir = '/data/ljo31/Lens/LensModels/twoband/'
files = [dir+'J0837_211',dir+'J0901_211',dir+'J0913_212',dir+'J1125_212',dir+'J1144_212',dir+'J1218_211',dir+'J1323_212',dir+'J1347_212',dir+'J1446_212',dir+'J1605_212',dir+'J1606_212',dir+'J1619_212',dir+'J2228_212']

names = ['J0837','J0901','J0913','J1125','J1144','J1218','J1323','J1347','J1446','J1605','J1606','J1619','J2228']

'''for ii in range(len(names)):
    result = np.load(files[ii])
    model = L.EELs(result,names[ii])
    model.Initialise()
    model.GetFits(plotresid=True)
    pl.figure()
    lp = result[0]
    if lp.shape[1]<10.:
        pl.plot(lp[:,0])
    else:
        pl.plot(lp)
    pl.show()'''

cats_m, cats_l, cats_h = [], [], []

for ii in range(len(names)):
    med,lo,hi = MakeTab(files[ii],names[ii])
    cats_m.append(med)
    cats_l.append(lo)
    cats_h.append(hi)



m, l, h = np.array(cats_m), np.array(cats_l), np.array(cats_h)


from astropy.io.fits import *
names = np.array(names)

c1 = Column(name='name', format='A5',array=names)
c2 = Column(name='mag v', format='D',array=m[:,0])
c3 = Column(name='mag i', format='D',array=m[:,1])
c4 = Column(name='mu v', format='D',array=m[:,2])
c5 = Column(name='mu i', format='D',array=m[:,3])
c6 = Column(name='Re v', format='D',array=m[:,4])
c7 = Column(name='Re i', format='D',array=m[:,5])
c8 = Column(name='v-i', format='D',array=m[:,6])

### and uncertainties - lower bounds
c2l = Column(name='mag v lo', format='D',array=l[:,0])
c3l = Column(name='mag i lo', format='D',array=l[:,1])
c4l = Column(name='mu v lo', format='D',array=l[:,2])
c5l = Column(name='mu i lo', format='D',array=l[:,3])
c6l = Column(name='Re v lo', format='D',array=l[:,4])
c7l = Column(name='Re i lo', format='D',array=l[:,5])
c8l = Column(name='v-i lo', format='D',array=l[:,6])


### and upper bounds
c2h = Column(name='mag v hi', format='D',array=h[:,0])
c3h = Column(name='mag i hi', format='D',array=h[:,1])
c4h = Column(name='mu v hi' , format='D',array=h[:,2])
c5h = Column(name='mu i hi', format='D',array=h[:,3])
c6h = Column(name='Re v hi', format='D',array=h[:,4])
c7h = Column(name='Re i hi', format='D',array=h[:,5])
c8h = Column(name='v-i hi', format='D',array=h[:,6])



coldefs = ColDefs([c1,c2,c3,c4,c5,c6,c7,c8,c2l,c3l,c4l,c5l,c6l,c7l,c8l,c2h,c3h,c4h,c5h,c6h,c7h,c8h])
tbhdu = BinTableHDU.from_columns(coldefs)
tbhdu.writeto('/data/ljo31/Lens/LensParams/Phot_2src_huge_new_new.fits',clobber=True)

for i in range(5):
    print r'%%%'

dir = '/data/ljo31/Lens/LensModels/twoband/'
files = [dir+'J0837_211',dir+'J0901_211',dir+'J0913_211',dir+'J1125_211',dir+'J1144_211',dir+'J1218_211',dir+'J1323_211','/data/ljo31/Lens/J1347/twoband_1_ctd_9',dir+'J1446_211',dir+'J1605_211',dir+'J1606_211',dir+'J1619_211',dir+'J2228_211']
names = ['J0837','J0901','J0913','J1125','J1144','J1218','J1323','J1347','J1446','J1605','J1606','J1619','J2228']

cats_m, cats_l, cats_h = [], [], []

for ii in range(len(names)):
    med,lo,hi = MakeTab(files[ii],names[ii])
    cats_m.append(med)
    cats_l.append(lo)
    cats_h.append(hi)



m, l, h = np.array(cats_m), np.array(cats_l), np.array(cats_h)


from astropy.io.fits import *
names = np.array(names)

c1 = Column(name='name', format='A5',array=names)
c2 = Column(name='mag v', format='D',array=m[:,0])
c3 = Column(name='mag i', format='D',array=m[:,1])
c4 = Column(name='mu v', format='D',array=m[:,2])
c5 = Column(name='mu i', format='D',array=m[:,3])
c6 = Column(name='Re v', format='D',array=m[:,4])
c7 = Column(name='Re i', format='D',array=m[:,5])
c8 = Column(name='v-i', format='D',array=m[:,6])

### and uncertainties - lower bounds
c2l = Column(name='mag v lo', format='D',array=l[:,0])
c3l = Column(name='mag i lo', format='D',array=l[:,1])
c4l = Column(name='mu v lo', format='D',array=l[:,2])
c5l = Column(name='mu i lo', format='D',array=l[:,3])
c6l = Column(name='Re v lo', format='D',array=l[:,4])
c7l = Column(name='Re i lo', format='D',array=l[:,5])
c8l = Column(name='v-i lo', format='D',array=l[:,6])


### and upper bounds
c2h = Column(name='mag v hi', format='D',array=h[:,0])
c3h = Column(name='mag i hi', format='D',array=h[:,1])
c4h = Column(name='mu v hi' , format='D',array=h[:,2])
c5h = Column(name='mu i hi', format='D',array=h[:,3])
c6h = Column(name='Re v hi', format='D',array=h[:,4])
c7h = Column(name='Re i hi', format='D',array=h[:,5])
c8h = Column(name='v-i hi', format='D',array=h[:,6])



coldefs = ColDefs([c1,c2,c3,c4,c5,c6,c7,c8,c2l,c3l,c4l,c5l,c6l,c7l,c8l,c2h,c3h,c4h,c5h,c6h,c7h,c8h])
tbhdu = BinTableHDU.from_columns(coldefs)
tbhdu.writeto('/data/ljo31/Lens/LensParams/Phot_1src_huge_new_new.fits',clobber=True)
