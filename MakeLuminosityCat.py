import EELsModels as L
import numpy as np

def MakeTab(model):
    Lvs,Lbs,Lis = [],[],[]
    model.MakeDict()
    model.BuildLenses()
    model.BuildGalaxies()
    model.BuildSources()
    model.EasyAddImages()
    model.GetFits(plotresid=False)
    mags = model.GetIntrinsicMags()
    Res = model.GetSourceSize(kpc=True)
    for i in range(len(ages)):
        lv,li,lb = model.GetPhotometryAtAge(K=False,age=ages[i])
        Lvs.append(lv)
        Lbs.append(lb)
        Lis.append(li)
    return np.array(Lvs),np.array(Lis),np.array(Lbs)

ages = ['0.010','0.125','0.250','0.375','0.500','0.625','0.750','0.875','1.000','1.250','1.500','1.700','1.750','2.000','2.200','2.250','2.500','2.750','3.000','3.250','3.500','3.750','4.000','4.250','4.500','4.750','5.000','5.250','5.500','5.750','6.000','7.000','8.000','9.000','10.00','12.00','15.00','20.00']

lumvs,lumis,lumbs = [],[],[]
result = np.load('/data/ljo31/Lens/LensModels/J0837_211')
model = L.EELs(result,name='J0837')
print 'J0837'
Lvs,Lis,Lbs=MakeTab(model)
lumvs.append(Lvs.flatten())
lumis.append(Lis.flatten())
lumbs.append(Lbs.flatten())


result = np.load('/data/ljo31/Lens/LensModels/J0901_211')
model = L.EELs(result,name='J0901')
print 'J0901'
Lvs,Lis,Lbs=MakeTab(model)
lumvs.append(Lvs.flatten())
lumis.append(Lis.flatten())
lumbs.append(Lbs.flatten())

result = np.load('/data/ljo31/Lens/LensModels/J0913_211')
model = L.EELs(result,name='J0913')
print 'J0913'
Lvs,Lis,Lbs=MakeTab(model)
lumvs.append(Lvs.flatten())
lumis.append(Lis.flatten())
lumbs.append(Lbs.flatten())

result = np.load('/data/ljo31/Lens/LensModels/J1125_211')
model = L.EELs(result,name='J1125')
print 'J1125'
Lvs,Lis,Lbs=MakeTab(model)
lumvs.append(Lvs.flatten())
lumis.append(Lis.flatten())
lumbs.append(Lbs.flatten())

result = np.load('/data/ljo31/Lens/LensModels/J1144_211')
model = L.EELs(result,name='J1144')
print 'J1144'
Lvs,Lis,Lbs=MakeTab(model)
lumvs.append(Lvs.flatten())
lumis.append(Lis.flatten())
lumbs.append(Lbs.flatten())


result = np.load('/data/ljo31/Lens/LensModels/J1218_211')
model = L.EELs(result,name='J1218')
print 'J1218'
Lvs,Lis,Lbs=MakeTab(model)
lumvs.append(Lvs.flatten())
lumis.append(Lis.flatten())
lumbs.append(Lbs.flatten())

result = np.load('/data/ljo31/Lens/LensModels/J1323_211') 
model = L.EELs(result,name='J1323')
print 'J1323'
Lvs,Lis,Lbs=MakeTab(model)
lumvs.append(Lvs.flatten())
lumis.append(Lis.flatten())
lumbs.append(Lbs.flatten())

result = np.load('/data/ljo31/Lens/LensModels/J1347_211')
model = L.EELs(result,name='J1347')
print 'J1347'
Lvs,Lis,Lbs=MakeTab(model)
lumvs.append(Lvs.flatten())
lumis.append(Lis.flatten())
lumbs.append(Lbs.flatten())

result = np.load('/data/ljo31/Lens/LensModels/J1446_211')
model = L.EELs(result,name='J1446')
print 'J1446'
Lvs,Lis,Lbs=MakeTab(model)
lumvs.append(Lvs.flatten())
lumis.append(Lis.flatten())
lumbs.append(Lbs.flatten())

result = np.load('/data/ljo31/Lens/LensModels/J1605_211') 
model = L.EELs(result,name='J1605')
print 'J1605'
Lvs,Lis,Lbs=MakeTab(model)
lumvs.append(Lvs.flatten())
lumis.append(Lis.flatten())
lumbs.append(Lbs.flatten())

result = np.load('/data/ljo31/Lens/LensModels/J1606_211')
model = L.EELs(result,name='J1606')
print 'J1606'
Lvs,Lis,Lbs=MakeTab(model)
lumvs.append(Lvs.flatten())
lumis.append(Lis.flatten())
lumbs.append(Lbs.flatten())

result = np.load('/data/ljo31/Lens/LensModels/J2228_211')
model = L.EELs(result,name='J2228')
print 'J2228'
Lvs,Lis,Lbs=MakeTab(model)
lumvs.append(Lvs.flatten())
lumis.append(Lis.flatten())
lumbs.append(Lbs.flatten())


lumbs = np.array(lumbs)
lumvs = np.array(lumvs)
lumis = np.array(lumis)

from astropy.io.fits import *
names = np.array(['J0837','J0901','J0913','J1125','J1144','J1218','J1323','J1347','J1446','J1605','J1606','J2228'])
c1 = Column(name='name', format='A5',array=names)
c2 = Column(name='0.010', format='D',array=lumbs[:,0])
c3 = Column(name='0.125', format='D',array=lumbs[:,1])
c4 = Column(name='0.250', format='D',array=lumbs[:,2])
c5 = Column(name='0.375', format='D',array=lumbs[:,3])
c6 = Column(name='0.500', format='D',array=lumbs[:,4])
c7 = Column(name='0.625', format='D',array=lumbs[:,5])
c8 = Column(name='0.750', format='D',array=lumbs[:,6])
c9 = Column(name='0.875', format='D',array=lumbs[:,7])
c10 = Column(name='1.000', format='D',array=lumbs[:,8])
c11 = Column(name='1.250', format='D',array=lumbs[:,9])
c12 = Column(name='1.500', format='D',array=lumbs[:,10])
c13 = Column(name='1.700', format='D',array=lumbs[:,11])
c14 = Column(name='1.750', format='D',array=lumbs[:,12])
c15 = Column(name='2.000', format='D',array=lumbs[:,13])
c16 = Column(name='2.200', format='D',array=lumbs[:,14])
c17 = Column(name='2.250', format='D',array=lumbs[:,15])
c18 = Column(name='2.500', format='D',array=lumbs[:,16])
c19 = Column(name='2.750', format='D',array=lumbs[:,17])
c20 = Column(name='3.000', format='D',array=lumbs[:,18])
c21 = Column(name='3.250', format='D',array=lumbs[:,19])
c22 = Column(name='3.500', format='D',array=lumbs[:,20])
c23 = Column(name='3.750', format='D',array=lumbs[:,21])
c24 = Column(name='4.000', format='D',array=lumbs[:,22])
c25 = Column(name='4.250', format='D',array=lumbs[:,23])
c26 = Column(name='4.500', format='D',array=lumbs[:,24])
c27 = Column(name='4.750', format='D',array=lumbs[:,25])
c28 = Column(name='5.000', format='D',array=lumbs[:,26])
c29 = Column(name='5.250', format='D',array=lumbs[:,27])
c30 = Column(name='5.500', format='D',array=lumbs[:,28])
c31 = Column(name='5.750', format='D',array=lumbs[:,29])
c32 = Column(name='6.000', format='D',array=lumbs[:,30])
c33 = Column(name='7.000', format='D',array=lumbs[:,31])
c34 = Column(name='8.000', format='D',array=lumbs[:,32])
c35 = Column(name='9.000', format='D',array=lumbs[:,33])
c36 = Column(name='10.00', format='D',array=lumbs[:,34])
c37 = Column(name='12.00', format='D',array=lumbs[:,35])
c38 = Column(name='15.00', format='D',array=lumbs[:,36])
c39 = Column(name='20.00', format='D',array=lumbs[:,37])


coldefs = ColDefs([c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,c21,c22,c23,c24,c25,c26,c27,c28,c29,c30,c31,c32,c33,c34,c35,c36,c37,c38,c39])
tbhdu = BinTableHDU.from_columns(coldefs)
tbhdu.writeto('/data/ljo31/Lens/LensParams/Lumb_ages_1src.fits',clobber=True)


## 2 src models
lumvs,lumis,lumbs = [],[],[]


#result = np.load('/data/ljo31/Lens/LensModels/J0913_212_concentric')
#model = L.EELs(result,name='J0913')
#print 'J0913 concentric'
#MakeTab(model)

result = np.load('/data/ljo31/Lens/LensModels/J0913_212_nonconcentric')
model = L.EELs(result,name='J0913')
print 'J0913 nonconcentric'
Lvs,Lis,Lbs=MakeTab(model)
lumvs.append(Lvs.flatten())
lumis.append(Lis.flatten())
lumbs.append(Lbs.flatten())



#result = np.load('/data/ljo31/Lens/LensModels/J1125_212_concentric')
#model = L.EELs(result,name='J1125')
#print 'J1125 concentric'
#MakeTab(model)

result = np.load('/data/ljo31/Lens/LensModels/J1125_212_nonconcentric')
model = L.EELs(result,name='J1125')
print 'J1125 nonconcentric'
Lvs,Lis,Lbs=MakeTab(model)
lumvs.append(Lvs.flatten())
lumis.append(Lis.flatten())
lumbs.append(Lbs.flatten())



result = np.load('/data/ljo31/Lens/LensModels/J1144_212_allparams')
model = L.EELs(result,name='J1144')
print 'J1144'
Lvs,Lis,Lbs=MakeTab(model)
lumvs.append(Lvs.flatten())
lumis.append(Lis.flatten())
lumbs.append(Lbs.flatten())


result = np.load('/data/ljo31/Lens/LensModels/J1323_212') 
model = L.EELs(result,name='J1323')
print 'J1323'
Lvs,Lis,Lbs=MakeTab(model)
lumvs.append(Lvs.flatten())
lumis.append(Lis.flatten())
lumbs.append(Lbs.flatten())


result = np.load('/data/ljo31/Lens/LensModels/J1347_112')
model = L.EELs(result,name='J1347')
print 'J1347'
Lvs,Lis,Lbs=MakeTab(model)
lumvs.append(Lvs.flatten())
lumis.append(Lis.flatten())
lumbs.append(Lbs.flatten())



result = np.load('/data/ljo31/Lens/LensModels/J1446_212')
model = L.EELs(result,name='J1446')
print 'J1446'
Lvs,Lis,Lbs=MakeTab(model)
lumvs.append(Lvs.flatten())
lumis.append(Lis.flatten())
lumbs.append(Lbs.flatten())



result = np.load('/data/ljo31/Lens/LensModels/J1605_212_final') 
model = L.EELs(result,name='J1605')
Lvs,Lis,Lbs=MakeTab(model)
lumvs.append(Lvs.flatten())
lumis.append(Lis.flatten())
lumbs.append(Lbs.flatten())




result = np.load('/data/ljo31/Lens/LensModels/J1606_112')
model = L.EELs(result,name='J1606')
Lvs,Lis,Lbs=MakeTab(model)
lumvs.append(Lvs.flatten())
lumis.append(Lis.flatten())
lumbs.append(Lbs.flatten())



result = np.load('/data/ljo31/Lens/LensModels/J2228_212')
model = L.EELs(result,name='J2228')
Lvs,Lis,Lbs=MakeTab(model)
lumvs.append(Lvs.flatten())
lumis.append(Lis.flatten())
lumbs.append(Lbs.flatten())




lumbs = np.array(lumbs)
lumvs = np.array(lumvs)
lumis = np.array(lumis)

from astropy.io.fits import *
names = np.array(['J0913','J1125','J1144','J1323','J1347','J1446','J1605','J1606','J2228'])
c1 = Column(name='name', format='A5',array=names)
c2 = Column(name='0.010', format='D',array=lumbs[:,0])
c3 = Column(name='0.125', format='D',array=lumbs[:,1])
c4 = Column(name='0.250', format='D',array=lumbs[:,2])
c5 = Column(name='0.375', format='D',array=lumbs[:,3])
c6 = Column(name='0.500', format='D',array=lumbs[:,4])
c7 = Column(name='0.625', format='D',array=lumbs[:,5])
c8 = Column(name='0.750', format='D',array=lumbs[:,6])
c9 = Column(name='0.875', format='D',array=lumbs[:,7])
c10 = Column(name='1.000', format='D',array=lumbs[:,8])
c11 = Column(name='1.250', format='D',array=lumbs[:,9])
c12 = Column(name='1.500', format='D',array=lumbs[:,10])
c13 = Column(name='1.700', format='D',array=lumbs[:,11])
c14 = Column(name='1.750', format='D',array=lumbs[:,12])
c15 = Column(name='2.000', format='D',array=lumbs[:,13])
c16 = Column(name='2.200', format='D',array=lumbs[:,14])
c17 = Column(name='2.250', format='D',array=lumbs[:,15])
c18 = Column(name='2.500', format='D',array=lumbs[:,16])
c19 = Column(name='2.750', format='D',array=lumbs[:,17])
c20 = Column(name='3.000', format='D',array=lumbs[:,18])
c21 = Column(name='3.250', format='D',array=lumbs[:,19])
c22 = Column(name='3.500', format='D',array=lumbs[:,20])
c23 = Column(name='3.750', format='D',array=lumbs[:,21])
c24 = Column(name='4.000', format='D',array=lumbs[:,22])
c25 = Column(name='4.250', format='D',array=lumbs[:,23])
c26 = Column(name='4.500', format='D',array=lumbs[:,24])
c27 = Column(name='4.750', format='D',array=lumbs[:,25])
c28 = Column(name='5.000', format='D',array=lumbs[:,26])
c29 = Column(name='5.250', format='D',array=lumbs[:,27])
c30 = Column(name='5.500', format='D',array=lumbs[:,28])
c31 = Column(name='5.750', format='D',array=lumbs[:,29])
c32 = Column(name='6.000', format='D',array=lumbs[:,30])
c33 = Column(name='7.000', format='D',array=lumbs[:,31])
c34 = Column(name='8.000', format='D',array=lumbs[:,32])
c35 = Column(name='9.000', format='D',array=lumbs[:,33])
c36 = Column(name='10.00', format='D',array=lumbs[:,34])
c37 = Column(name='12.00', format='D',array=lumbs[:,35])
c38 = Column(name='15.00', format='D',array=lumbs[:,36])
c39 = Column(name='20.00', format='D',array=lumbs[:,37])



coldefs = ColDefs([c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,c21,c22,c23,c24,c25,c26,c27,c28,c29,c30,c31,c32,c33,c34,c35,c36,c37,c38,c39])
tbhdu = BinTableHDU.from_columns(coldefs)
tbhdu.writeto('/data/ljo31/Lens/LensParams/Lumb_ages_2src.fits',clobber=True)
