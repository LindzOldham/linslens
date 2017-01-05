import numpy as np
from linslens import EELsModels as L

def BT(model):
    model.Initialise()
    model.GetIntrinsicMags()
    BTv,BTi = model.BT()
    print '%.2f'%BTv,'%.2f'%BTi

result = np.load('/data/ljo31/Lens/LensModels/J0913_212_nonconcentric')
model = L.EELs(result,name='J0913')
BT(model)

result = np.load('/data/ljo31/Lens/LensModels/J1125_212')
model = L.EELs(result,name='J1125')
BT(model)

result = np.load('/data/ljo31/Lens/LensModels/J1144_212_allparams')
model = L.EELs(result,name='J1144')
BT(model)

result = np.load('/data/ljo31/Lens/LensModels/J1323_212') 
model = L.EELs(result,name='J1323')
BT(model)

result = np.load('/data/ljo31/Lens/LensModels/J1347_112')
model = L.EELs(result,name='J1347')
BT(model)

result = np.load('/data/ljo31/Lens/LensModels/J1446_212')
model = L.EELs(result,name='J1446')
BT(model)

result = np.load('/data/ljo31/Lens/LensModels/J1605_212_final') 
model = L.EELs(result,name='J1605')
BT(model)

result = np.load('/data/ljo31/Lens/LensModels/J1606_112')
model = L.EELs(result,name='J1606')
BT(model)

result = np.load('/data/ljo31/Lens/LensModels/J1619_212')
model = L.EELs(result,name='J1619')
BT(model)

result = np.load('/data/ljo31/Lens/LensModels/J2228_212')
model = L.EELs(result,name='J2228')
BT(model)


