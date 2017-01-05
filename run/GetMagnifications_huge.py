from linslens import EELsModels_huge as L
import numpy as np

def MakeTab(model):
    model.Initialise()
    mags = model.GetIntrinsicMags()
    Res = model.GetSourceSize(kpc=True)
    Mv,Mi = model.magnification()
    return Mv,Mi

dir = '/data/ljo31/Lens/LensModels/twoband/'
files = [dir+'J0837_211',dir+'J0901_211',dir+'J0913_211',dir+'J1125_211',dir+'J1144_211',dir+'J1218_211',dir+'J1323_211','/data/ljo31/Lens/J1347/twoband_1_ctd_9',dir+'J1446_211',dir+'J1605_211',dir+'J1606_211',dir+'J1619_211',dir+'J2228_211']
names = ['J0837','J0901','J0913','J1125','J1144','J1218','J1323','J1347','J1446','J1605','J1606','J1619','J2228']

cat = []

for ii in range(len(names)):
    result = np.load(files[ii])
    model = L.EELs(result,name=names[ii])
    Mv,Mi = MakeTab(model)
    cat.append((model.name,[Mv,Mi]))

cat = dict(cat)

np.save('/data/ljo31/Lens/LensParams/magnifications_211_huge',cat)


cat = []

dir = '/data/ljo31/Lens/LensModels/twoband/'
files = [dir+'J0837_211',dir+'J0901_211',dir+'J0913_212',dir+'J1125_212',dir+'J1144_212',dir+'J1218_211',dir+'J1323_212',dir+'J1347_212',dir+'J1446_212',dir+'J1605_212',dir+'J1606_212',dir+'J1619_212',dir+'J2228_212']
names = ['J0837','J0901','J0913','J1125','J1144','J1218','J1323','J1347','J1446','J1605','J1606','J1619','J2228']

for ii in range(len(names)):
    result = np.load(files[ii])
    model = L.EELs(result,name=names[ii])
    Mv,Mi = MakeTab(model)
    cat.append((model.name,[Mv,Mi]))


cat = dict(cat)

np.save('/data/ljo31/Lens/LensParams/magnifications_212_huge',cat)




