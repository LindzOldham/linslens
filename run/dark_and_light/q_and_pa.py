import numpy as np, pylab as pl
import indexTricks as iT
import pymc,cPickle
import glob

def compares(name):
    ## inference on dark halo STRUCTURE
    file = glob.glob('/data/ljo31b/EELs/galsub/emceeruns/'+name+'_parametric_*')
    file.sort()
    f = file[-1]
    while 1:
        if 'DPL' in f:
            file = file[:-1]
            f = file[-1]
        else:
            break
    haloResult = np.load(f)

    # stellar mass STRUCTURE result
    dir = '/data/ljo31/Lens/LensModels/twoband/'
    if name == 'J1144':
        result = np.load('/data/ljo31/Lens/J1144/twoband_darkandlightprep_ctd')
    else:
        try:
            result = np.load(dir+name+'_212')
        except:
            if name == 'J1347':
                result = np.load(dir+name+'_112')
            else:
                result = np.load(dir+name+'_211')

    lp,trace,dic,_=result
    if len(lp.shape)==3.:
        lp = lp[:,0]
        for key in dic.keys():
            dic[key] = dic[key][:,0]
    a1,a2 = np.unravel_index(lp.argmax(),lp.shape)
    
    lph,traceh,dich,_= haloResult
    a1h,a2h = np.unravel_index(lph[:,0].argmax(),lph[:,0].shape)

    q_h, pa_h = dich['Lens 1 q'][a1h,0,a2h], dich['Lens 1 pa'][a1h,0,a2h]
    q_s1, pa_s1 = dic['Galaxy 1 q'][a1,a2], dic['Galaxy 1 pa'][a1,a2]
    if 'Galaxy 2 q' in dic.keys():
        q_s2, pa_s2 = dic['Galaxy 2 q'][a1,a2], dic['Galaxy 2 pa'][a1,a2]
    else:
        q_s2,pa_s2 = 0.,0.
    q_l, pa_l = dic['Lens 1 q'][a1,a2], dic['Lens 1 pa'][a1,a2]

    print '%s'%name
    print '%s,%.3f,%.3f'%('halo',q_h,pa_h)
    print '%s,%.3f,%.3f'%('stellar 1',q_s1,pa_s1)
    print '%s,%.3f,%.3f'%('stellar 2',q_s2,pa_s2)
    print '%s,%.3f,%.3f'%('lens',q_l,pa_l)

names = np.array(['J0837','J0901','J0913','J1144'])
for name in names:
    compares(name)

