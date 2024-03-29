import numpy as np
from linslens import EELsModels_huge as L

dir = '/data/ljo31/Lens/LensModels/twoband/'
files = ['J0913','J1125','J1144','J1323','J1347','J1446','J1605','J1606','J1619','J2228']
for file in files:
    result = np.load(dir+file+'_212')
    model = L.EELs(result,file)
    model.Initialise()
    model.GetIntrinsicMags()
    BTv,BTi = model.BT()
    print model.name,'%.2f'%BTv,'%.2f'%BTi

for file in files:
    result = np.load(dir+file+'_212')
    model = L.EELs(result,file)
    model.Initialise()
    print model.name, '$', '%.2f'%model.Ddic['Source 1 n'], r'\pm', '%.2f'%np.mean((model.Ldic['Source 1 n'],model.Udic['Source 1 n'])), '$ & $', '%.2f'%model.Ddic['Source 2 n'], r'\pm', '%.2f'%np.mean((model.Ldic['Source 2 n'],model.Udic['Source 2 n']))

files = ['J0837','J0901','J0913','J1125','J1144','J1218','J1323','J1347','J1446','J1605','J1606','J1619','J2228']
for file in files:
    result = np.load(dir+file+'_211')
    model = L.EELs(result,file)
    model.Initialise()
    print model.name, '$', '%.2f'%model.Ddic['Source 1 n'], r'\pm', '%.2f'%np.mean((model.Ldic['Source 1 n'],model.Udic['Source 1 n'])), '$'
