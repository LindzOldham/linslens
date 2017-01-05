from linslens import EELsModels_huge as L
import numpy as np

def MakeTab(model):
    model.Initialise()
    lo,med,hi = model.Ldic, model.Ddic, model.Udic
    return lo,med,hi
    

cats_m, cats_l, cats_h = [], [], []


dir = '/data/ljo31/Lens/LensModels/twoband/'
files = [dir+'J0837_211',dir+'J0901_211',dir+'J0913_211',dir+'J1125_211',dir+'J1144_211',dir+'J1218_211',dir+'J1323_211','/data/ljo31/Lens/J1347/twoband_1_ctd_9',dir+'J1446_211',dir+'J1605_211',dir+'J1606_211',dir+'J1619_211',dir+'J2228_211']
names = ['J0837','J0901','J0913','J1125','J1144','J1218','J1323','J1347','J1446','J1605','J1606','J1619','J2228']

for ii in range(len(names)):
    result = np.load(files[ii])
    model = L.EELs(result,name=names[ii])
    lo,med,hi = MakeTab(model)
    cats_m.append([model.name,med])
    cats_l.append([model.name,lo])
    cats_h.append([model.name,hi])



m, l, h = dict(cats_m), dict(cats_l), dict(cats_h)

np.save('/data/ljo31/Lens/LensParams/Structure_1src_huge_new',[m,l,h])



## 2 src models
cats_m, cats_l, cats_h = [], [], []

dir = '/data/ljo31/Lens/LensModels/twoband/'
files = [dir+'J0837_211',dir+'J0901_211',dir+'J0913_212',dir+'J1125_212',dir+'J1144_212',dir+'J1218_211',dir+'J1323_212',dir+'J1347_212',dir+'J1446_212',dir+'J1605_212',dir+'J1606_212',dir+'J1619_212',dir+'J2228_212']
names = ['J0837','J0901','J0913','J1125','J1144','J1218','J1323','J1347','J1446','J1605','J1606','J1619','J2228']

for ii in range(len(names)):
    result = np.load(files[ii])
    model = L.EELs(result,name=names[ii])
    lo,med,hi = MakeTab(model)
    cats_m.append([model.name,med])
    cats_l.append([model.name,lo])
    cats_h.append([model.name,hi])


m, l, h = dict(cats_m), dict(cats_l), dict(cats_h)

np.save('/data/ljo31/Lens/LensParams/Structure_2src_huge_new',[m,l,h])

