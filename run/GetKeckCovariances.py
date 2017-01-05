import EELsKeckModels as L
import numpy as np

def MakeTab(model):
    model.MakeDict()
    model.BuildLenses()
    model.BuildGalaxies()
    model.BuildSources()
    model.EasyAddImages()
    model.GetFits(plotresid=False)
    model.GetIntrinsicMags()
    model.GetSourceSize()
    model.GetPhotometry()
    model.GetSB()
    model.GetRestSB()
    model.MakePDFDict()
    print len(model.dictionaries)
    model.GetPDFs(kpc=True)
    magPDF = model.magPDF
    rePDF = model.RePDF
    return magPDF, rePDF



result = np.load('/data/ljo31/Lens/J0837/Kp_0')
hstresult = np.load('/data/ljo31/Lens/LensModels/J0837_211')
model = L.EELs(result,hstresult,name='J0837')
name = 'J0837'
magPDF, RePDF = MakeTab(model)
np.save('/data/ljo31/Lens/Analysis/ReMagKeckPDFs_'+str(name),[RePDF,magPDF])


hstresult = np.load('/data/ljo31/Lens/LensModels/J0901_211')
result = np.load('/data/ljo31/Lens/J0901/Kp_0')
model = L.EELs(result,hstresult,name='J0901')
name = 'J0901'
magPDF, RePDF = MakeTab(model)
np.save('/data/ljo31/Lens/Analysis/ReMagKeckPDFs_'+str(name),[RePDF,magPDF])



hstresult = np.load('/data/ljo31/Lens/LensModels/J0913_211')
result = np.load('/data/ljo31/Lens/J0913/Kp_211_0')
model = L.EELs(result,hstresult,name='J0913')
name = 'J0913'
magPDF, RePDF = MakeTab(model)
np.save('/data/ljo31/Lens/Analysis/ReMagKeckPDFs_'+str(name),[RePDF,magPDF])


result = np.load('/data/ljo31/Lens/J1125/Kp_211_0')
hstresult = np.load('/data/ljo31/Lens/LensModels/J1125_211')
model = L.EELs(result,hstresult,name='J1125')
name = 'J1125'
magPDF, RePDF = MakeTab(model)
np.save('/data/ljo31/Lens/Analysis/ReMagKeckPDFs_'+str(name),[RePDF,magPDF])


result = np.load('/data/ljo31/Lens/J1144/Kp_211_0')
hstresult = np.load('/data/ljo31/Lens/LensModels/J1144_211')
model = L.EELs(result,hstresult,name='J1144')
name = 'J1144'
magPDF, RePDF = MakeTab(model)
np.save('/data/ljo31/Lens/Analysis/ReMagKeckPDFs_'+str(name),[RePDF,magPDF])



result = np.load('/data/ljo31/Lens/J1218/Kp_0')
hstresult = np.load('/data/ljo31/Lens/LensModels/J1218_211')
model = L.EELs(result,hstresult,name='J1218')
name = 'J1218'
magPDF, RePDF = MakeTab(model)
np.save('/data/ljo31/Lens/Analysis/ReMagKeckPDFs_'+str(name),[RePDF,magPDF])


result = np.load('/data/ljo31/Lens/J1323/Kp_211_7')
hstresult = np.load('/data/ljo31/Lens/LensModels/J1323_211') 
model = L.EELs(result,hstresult,name='J1323')
name = 'J1323'
magPDF, RePDF = MakeTab(model)
np.save('/data/ljo31/Lens/Analysis/ReMagKeckPDFs_'+str(name),[RePDF,magPDF])


result = np.load('/data/ljo31/Lens/J1347/Kp_211_0')
hstresult = np.load('/data/ljo31/Lens/LensModels/J1347_211')
model = L.EELs(result,hstresult,name='J1347')
name = 'J1347'
magPDF, RePDF = MakeTab(model)
np.save('/data/ljo31/Lens/Analysis/ReMagKeckPDFs_'+str(name),[RePDF,magPDF])


result = np.load('/data/ljo31/Lens/J1446/Kp_211_0')
hstresult = np.load('/data/ljo31/Lens/LensModels/J1446_211')
model = L.EELs(result,hstresult,name='J1446')
name = 'J1446'
magPDF, RePDF = MakeTab(model)
np.save('/data/ljo31/Lens/Analysis/ReMagKeckPDFs_'+str(name),[RePDF,magPDF])

result = np.load('/data/ljo31/Lens/J1605/Kp_211_0')
hstresult = np.load('/data/ljo31/Lens/LensModels/J1605_211') 
model = L.EELs(result,hstresult,name='J1605')
name = 'J1605'
magPDF, RePDF = MakeTab(model)
np.save('/data/ljo31/Lens/Analysis/ReMagKeckPDFs_'+str(name),[RePDF,magPDF])


result = np.load('/data/ljo31/Lens/J1606/Kp_211_0')
hstresult = np.load('/data/ljo31/Lens/LensModels/J1606_211')
model = L.EELs(result,hstresult,name='J1606')
name = 'J1606'
magPDF, RePDF = MakeTab(model)
np.save('/data/ljo31/Lens/Analysis/ReMagKeckPDFs_'+str(name),[RePDF,magPDF])

result = np.load('/data/ljo31/Lens/J2228/Kp_211_0')
hstresult = np.load('/data/ljo31/Lens/LensModels/J2228_211')
model = L.EELs(result,hstresult,name='J2228')
name = 'J2228'
magPDF, RePDF = MakeTab(model)
np.save('/data/ljo31/Lens/Analysis/ReMagKeckPDFs_'+str(name),[RePDF,magPDF])


'''
result = np.load('/data/ljo31/Lens/J0913/Kp_0')
hstresult = np.load('/data/ljo31/Lens/LensModels/J0913_212_concentric')
model = L.EELs(result,hstresult,name='J0913')
print 'J0913 concentric'
MakeTab(model)


result = np.load('/data/ljo31/Lens/J0913/Kp_0')
hstresult = np.load('/data/ljo31/Lens/LensModels/J0913_212_nonconcentric')
model = L.EELs(result,hstresult,name='J0913')
print 'J0913 nonconcentric'
MakeTab(model)


result = np.load('/data/ljo31/Lens/J1125/Kp_0')
hstresult = np.load('/data/ljo31/Lens/LensModels/J1125_212_concentric')
model = L.EELs(result,hstresult,name='J1125')
print 'J1125 concentric'
MakeTab(model)


result = np.load('/data/ljo31/Lens/J1125/Kp_0')
hstresult = np.load('/data/ljo31/Lens/LensModels/J1125_212_nonconcentric')
model = L.EELs(result,hstresult,name='J1125')
print 'J1125 nonconcentric'
MakeTab(model)

result = np.load('/data/ljo31/Lens/J1144/Kp_212_0')
hstresult = np.load('/data/ljo31/Lens/LensModels/J1144_212_allparams')
model = L.EELs(result,hstresult,name='J1144')
print 'J1144'
MakeTab(model)

result = np.load('/data/ljo31/Lens/J1323/Kp_0')
hstresult = np.load('/data/ljo31/Lens/LensModels/J1323_212') 
model = L.EELs(result,hstresult,name='J1323')
print 'J1323'
MakeTab(model)

result = np.load('/data/ljo31/Lens/J1347/Kp_0')
hstresult = np.load('/data/ljo31/Lens/LensModels/J1347_112')
model = L.EELs(result,hstresult,name='J1347')
print 'J1347'
MakeTab(model)

result = np.load('/data/ljo31/Lens/J1446/Kp_0')
hstresult = np.load('/data/ljo31/Lens/LensModels/J1446_212')
model = L.EELs(result,hstresult,name='J1446')
print 'J1446'
MakeTab(model)

result = np.load('/data/ljo31/Lens/J1605/Kp_0')
hstresult = np.load('/data/ljo31/Lens/LensModels/J1605_212_final') 
model = L.EELs(result,hstresult,name='J1605')
MakeTab(model)

result = np.load('/data/ljo31/Lens/J1606/Kp_1')
hstresult = np.load('/data/ljo31/Lens/LensModels/J1606_112')
model = L.EELs(result,hstresult,name='J1606')
print 'J1606'
MakeTab(model)

result = np.load('/data/ljo31/Lens/J2228/Kp_0')
hstresult = np.load('/data/ljo31/Lens/LensModels/J2228_212')
model = L.EELs(result,hstresult,name='J2228')
print 'J2228'
MakeTab(model)
'''
