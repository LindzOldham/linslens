import numpy as np, pylab as pl, pyfits as py
from linslens import EELsModels_huge as L

def MakeTab(file,name):
    result = np.load(file)
    model = L.EELs(result,name)
    model.Initialise()
    model.GetFits(plotresid=True)
    pl.show()

def MakeTab2(file,name):
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
    np.save('/data/ljo31/Lens/PDFs/211_huge__'+str(model.name), [np.array(model.muPDF),np.array(model.magPDF),np.array(model.RePDF)])
    return med,lo,hi

def MakeTab3(file,name):
    result = np.load(file)
    model = L.EELs(result,name)
    model.Initialise()
    lo,med,hi = model.Ldic, model.Ddic, model.Udic
    #model.GetFits(plotresid=True)
    return lo,med,hi

files = ['/data/ljo31/Lens/J0837/twoband_0', '/data/ljo31/Lens/J0913/twoband_1', '/data/ljo31/Lens/J1144/twoband_1', '/data/ljo31/Lens/J1218/twoband_0', '/data/ljo31/Lens/J1323/twoband_0', '/data/ljo31/Lens/J1347/twoband_0', '/data/ljo31/Lens/J1446/twoband_0', '/data/ljo31/Lens/J1605/twoband_0', '/data/ljo31/Lens/J1606/twoband_0', '/data/ljo31/Lens/J1619/twoband_0', '/data/ljo31/Lens/J2228/twoband_0']

names = ['J0837','J0913','J1144','J1218','J1323','J1347','J1446','J1605','J1606','J1619','J2228']

cats_m, cats_l, cats_h = [], [], []

for ii in range(len(names)):
    lo,med,hi = MakeTab3(files[ii],names[ii])
    cats_m.append([names[ii],med])
    cats_l.append([names[ii],lo])
    cats_h.append([names[ii],hi])
    

m, l, h = dict(cats_m), dict(cats_l), dict(cats_h)

np.save('/data/ljo31/Lens/LensParams/Structure_1src_huge',[m,l,h])



'''
cats_m, cats_l, cats_h = [], [], []

for ii in range(len(names)):
    med,lo,hi = MakeTab2(files[ii],names[ii])
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
tbhdu.writeto('/data/ljo31/Lens/LensParams/Phot_1src_huge.fits',clobber=True)
'''
