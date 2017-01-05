from linslens import EELsLensModels_huge as L
import numpy as np

def MakeTab(model):
    model.Initialise()
    model.GalaxyObsMag()
    model.GalaxySize(kpc=True)
    model.MakePDFDict()
    model.GetPDFs()
    med,lo,hi = model.UncertaintiesFromPDF(makecat=True)
    np.save('/data/ljo31/Lens/PDFs/211_lensgal_huge_'+str(model.name), [np.array([model.muPDF]),np.array([model.magPDF]),np.array([ model.fitPDF]), model.viPDF,np.array([model.RePDF])])
    #model.EasyAddPDFs()
    magPDF = model.magPDF
    RePDF = model.RePDF
    #print 'ICI',magPDF, RePDF
    return magPDF, RePDF


result = np.load('/data/ljo31/Lens/LensModels/twoband/J0837_211')
model = L.EELs(result,name='J0837')
name = 'J0837'
magPDF, RePDF = MakeTab(model)
np.save('/data/ljo31/Lens/Analysis/ReMagPDFs_twoband_lensgals_'+str(name),[RePDF,magPDF])

result = np.load('/data/ljo31/Lens/LensModels/twoband/J0901_211')
model = L.EELs(result,name='J0901')
name = 'J0901'
magPDF, RePDF = MakeTab(model)
np.save('/data/ljo31/Lens/Analysis/ReMagPDFs_twoband_lensgals_'+str(name),[RePDF,magPDF])


result = np.load('/data/ljo31/Lens/LensModels/twoband/J0913_211')
model = L.EELs(result,name='J0913')
name = 'J0913'
magPDF, RePDF = MakeTab(model)
np.save('/data/ljo31/Lens/Analysis/ReMagPDFs_twoband_lensgals_'+str(name),[RePDF,magPDF])

result = np.load('/data/ljo31/Lens/LensModels/twoband/J1125_211')
model = L.EELs(result,name='J1125')
name = 'J1125'
magPDF, RePDF = MakeTab(model)
np.save('/data/ljo31/Lens/Analysis/ReMagPDFs_twoband_lensgals_'+str(name),[RePDF,magPDF])

result = np.load('/data/ljo31/Lens/LensModels/twoband/J1144_211')
model = L.EELs(result,name='J1144')
name = 'J1144'
magPDF, RePDF = MakeTab(model)
np.save('/data/ljo31/Lens/Analysis/ReMagPDFs_twoband_lensgals_'+str(name),[RePDF,magPDF])

result = np.load('/data/ljo31/Lens/LensModels/twoband/J1218_211')
model = L.EELs(result,name='J1218')
name = 'J1218'
magPDF, RePDF = MakeTab(model)
np.save('/data/ljo31/Lens/Analysis/ReMagPDFs_twoband_lensgals_'+str(name),[RePDF,magPDF])

result = np.load('/data/ljo31/Lens/LensModels/twoband/J1323_211') 
model = L.EELs(result,name='J1323')
name = 'J1323'
magPDF, RePDF = MakeTab(model)
np.save('/data/ljo31/Lens/Analysis/ReMagPDFs_twoband_lensgals_'+str(name),[RePDF,magPDF])

result = np.load('/data/ljo31/Lens/J1347/twoband_0')
model = L.EELs(result,name='J1347')
name = 'J1347'
magPDF, RePDF = MakeTab(model)
np.save('/data/ljo31/Lens/Analysis/ReMagPDFs_twoband_lensgals_'+str(name),[RePDF,magPDF])

result = np.load('/data/ljo31/Lens/LensModels/twoband/J1446_211')
model = L.EELs(result,name='J1446')
name = 'J1446'
magPDF, RePDF = MakeTab(model)
np.save('/data/ljo31/Lens/Analysis/ReMagPDFs_twoband_lensgals_'+str(name),[RePDF,magPDF])

result = np.load('/data/ljo31/Lens/LensModels/twoband/J1605_211') 
model = L.EELs(result,name='J1605')
name = 'J1605'
magPDF, RePDF = MakeTab(model)
np.save('/data/ljo31/Lens/Analysis/ReMagPDFs_twoband_lensgals_'+str(name),[RePDF,magPDF])

result = np.load('/data/ljo31/Lens/LensModels/twoband/J1606_211')
model = L.EELs(result,name='J1606')
name = 'J1606'
magPDF, RePDF = MakeTab(model)
np.save('/data/ljo31/Lens/Analysis/ReMagPDFs_twoband_lensgals_'+str(name),[RePDF,magPDF])

result = np.load('/data/ljo31/Lens/LensModels/twoband/J1619_211')
model = L.EELs(result,name='J1619')
name = 'J1619'
magPDF, RePDF = MakeTab(model)
np.save('/data/ljo31/Lens/Analysis/ReMagPDFs_twoband_lensgals_'+str(name),[RePDF,magPDF])


result = np.load('/data/ljo31/Lens/LensModels/twoband/J2228_211')
model = L.EELs(result,name='J2228')
name = 'J2228'
magPDF, RePDF = MakeTab(model)
np.save('/data/ljo31/Lens/Analysis/ReMagPDFs_twoband_lensgals_'+str(name),[RePDF,magPDF])



'''

result = np.load('/data/ljo31/Lens/LensModels/J0913_212_concentric')
model = L.EELs(result,name='J0913')
print 'J0913 concentric'
MakeTab(model)

result = np.load('/data/ljo31/Lens/LensModels/J0913_212_nonconcentric')
model = L.EELs(result,name='J0913')
print 'J0913 nonconcentric'
MakeTab(model)

result = np.load('/data/ljo31/Lens/LensModels/J1125_212_concentric')
model = L.EELs(result,name='J1125')
print 'J1125 concentric'
MakeTab(model)

result = np.load('/data/ljo31/Lens/LensModels/J1125_212_nonconcentric')
model = L.EELs(result,name='J1125')
print 'J1125 nonconcentric'
MakeTab(model)


result = np.load('/data/ljo31/Lens/LensModels/J1144_212_allparams')
model = L.EELs(result,name='J1144')
print 'J1144'
MakeTab(model)


result = np.load('/data/ljo31/Lens/LensModels/J1323_212') 
model = L.EELs(result,name='J1323')
print 'J1323'
MakeTab(model)

result = np.load('/data/ljo31/Lens/LensModels/J1347_112')
model = L.EELs(result,name='J1347')
print 'J1347'
MakeTab(model)

result = np.load('/data/ljo31/Lens/LensModels/J1446_211')
model = L.EELs(result,name='J1446')
print 'J1446'
MakeTab(model)

result = np.load('/data/ljo31/Lens/LensModels/J1605_212_final') 
model = L.EELs(result,name='J1605')
MakeTab(model)


result = np.load('/data/ljo31/Lens/LensModels/J2228_212')
model = L.EELs(result,name='J2228')
print 'J2228'
MakeTab(model)


result = np.load('/data/ljo31/Lens/LensModels/J1606_211')
model = L.EELs(result,name='J1606')
print 'J1606'
MakeTab(model)

'''
