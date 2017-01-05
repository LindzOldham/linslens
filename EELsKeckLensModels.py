import pyfits as py, numpy as np
from imageSim import SBObjects,convolve
from pylens import *
import indexTricks as iT
from scipy import optimize
from scipy.interpolate import splrep, splev, splint
from astLib import astCalc
import EELsImages_huge_forpaper as EELsImages, Plotter
from GetPhotometry import *
import SBBModels, SBBProfiles
import itertools
from Plotter import *

source_redshifts = dict([('J0837',0.6411),('J0901',0.586),('J0913',0.539),('J1125',0.689),('J1144',0.706),('J1218',0.6009),('J1248',0.528),('J1323',0.4641),('J1347',0.63),('J1446',0.585),('J1605',0.542),('J1606',0.6549),('J1619',0.6137),('J2228',0.4370)])

lens_redshifts = dict([('J0837',0.4256),('J0901',0.311),('J0913',0.395),('J1125',0.442),('J1144',0.372),('J1218',0.3182),('J1248',0.304),('J1323',0.3194),('J1347',0.39),('J1446',0.317),('J1605',0.306),('J1606',0.3816),('J1619',0.3638),('J2228',0.2391)])

ZPs = np.load('/data/ljo31/Lens/LensParams/Keck_zeropoints.npy')[()]

class EELs:
    def __init__(self,result,hstresult,name,fits=None):
        self.result = result
        if self.result is not None:
            self.lp, self.trace, self.dic,_ = self.result
        self.lp,self.trace = self.lp[:,0], self.trace[:,0]
        for key in self.dic.keys():
            self.dic[key] = self.dic[key][:,0]
        hlp,htrace,hdic1,_ = hstresult
        hdic = []
        if hlp.shape[1]>10:
            a1,a2 = np.unravel_index(hlp.argmax(), hlp.shape)
            for key in hdic1.keys():
                if key not in self.dic.keys():
                    hdic.append((key,hdic1[key][a1,a2]))
        else:
            a1,a3 = np.unravel_index(hlp[:,0].argmax(), hlp[:,0].shape)
            for key in hdic1.keys():
                if key not in self.dic.keys():
                   hdic.append((key,hdic1[key][a1,0,a3])) 
        hdic = dict(hdic)
        self.hdic = hdic
        self.name = name
        self.fits = fits
        self.ZP = ZPs[name]

                          
    def MakeDict(self):
        if 'Source 2 re' not in self.dic.keys():
            self.srcno = 1
        else:
            self.srcno = 2
        if 'Galaxy 2 re' not in self.hdic.keys() and 'Galaxy 3 re' not in self.hdic.keys() and 'Galaxy 2 re' not in self.dic.keys() and 'Galaxy 3 re' not in self.dic.keys():
            self.galno = 1
        elif 'Galaxy 2 re' in self.hdic.keys() and 'Galaxy 3 re' not in self.hdic.keys():
            self.galno = 2
        elif 'Galaxy 2 re' in self.dic.keys() and 'Galaxy 3 re' not in self.dic.keys():
            self.galno = 2
        else:
            self.galno = 3
        if self.srcno==2 and (self.name == 'J1323' or self.name == 'J1347' or self.name == 'J0837') :
            print "repositioning J1323's/J1347's/J0837's source"
            self.dic['Source 2 x'] -= self.hdic['Lens 1 x']
            self.dic['Source 2 y'] -= self.hdic['Lens 1 y']
        if 'Source 1 pa' not in self.dic.keys():
            self.dic['Source 1 pa'] = self.dic['Source 2 pa']
        if self.galno == 2 and 'Galaxy 2 pa' not in self.hdic.keys() and 'Galaxy 2 pa' not in self.dic.keys():
            self.hdic['Galaxy 2 pa'] = self.hdic['Galaxy 1 pa']
            
        l,u,d = [], [], []
        for key in self.dic.keys():
            f = self.dic[key].reshape((self.trace.shape[0]*self.trace.shape[1]))
            lo,med,up = np.percentile(f,50)-np.percentile(f,16), np.percentile(f,50), np.percentile(f,84)-np.percentile(f,50) 
            d.append((key,med))
            l.append((key,lo))
            u.append((key,up))
        
        if 'Source 1 x' not in self.dic.keys():
            for key in 'x', 'y':
                f = self.dic['Source 2 '+key].reshape((self.trace.shape[0]*self.trace.shape[1]))
                lo,med,up = np.percentile(f,50)-np.percentile(f,16), np.percentile(f,50), np.percentile(f,84)-np.percentile(f,50) 
                d.append(('Source 1 '+key,med))
                l.append(('Source 1 '+key,lo))
                u.append(('Source 1 '+key,up))
                self.dic['Source 1 '+key] = self.dic['Source 2 '+key]
        if 'Galaxy 2 x' not in self.hdic.keys() and 'Galaxy 1 x' in self.hdic.keys() and self.galno == 2:
            for key in 'x', 'y':
                self.hdic['Galaxy 2 '+key] = self.hdic['Galaxy 1 '+key]
        if 'Galaxy 2 x' not in self.dic.keys() and 'Galaxy 1 x' in self.dic.keys() and self.galno == 2:
            for key in 'x', 'y':
                f = self.dic['Galaxy 1 '+key].reshape((self.trace.shape[0]*self.trace.shape[1]))
                lo,med,up = np.percentile(f,50)-np.percentile(f,16), np.percentile(f,50), np.percentile(f,84)-np.percentile(f,50) 
                d.append(('Galaxy 2 '+key,med))
                l.append(('Galaxy 2 '+key,lo))
                u.append(('Galaxy 2 '+key,up))
                self.dic['Galaxy 2 '+key] = self.dic['Galaxy 1 '+key]

        self.Ddic = dict(d)                    
        self.Ldic = dict(l)
        self.Udic = dict(u)

        # if name == j1619, add lens params to hdic for now
        if self.name == 'J1619':
            for key in 'x','y','b','eta','pa','q':
                self.hdic['Lens 1 '+key] = self.Ddic['Lens 1 '+key]
        #return self.Ddic,self.Ldic,self.Udic

    def PrintTable(self):
        ### print source model
        print r'\begin{table}[H]'
        print r'\centering'
        print r'\begin{tabular}{|c|cccccc|}\hline'
        print r' object & x & y & re & n & pa & q \\\hline'
        print 'source 1 & $', '%.2f'%(self.Ddic['Source 1 x']+self.hdic['Lens 1 x']), '_{-', '%.2f'%self.Ldic['Source 1 x'],'}^{+','%.2f'%self.Udic['Source 1 x'], '}$ & $', '%.2f'%(self.Ddic['Source 1 y']+self.hdic['Lens 1 y']),'_{-', '%.2f'%self.Ldic['Source 1 y'],'}^{+', '%.2f'%self.Udic['Source 1 y'], '}$ & $', '%.2f'%self.Ddic['Source 1 re'],'_{-', '%.2f'%self.Ldic['Source 1 re'],'}^{+', '%.2f'%self.Udic['Source 1 re'], '}$ & $', '%.2f'%self.Ddic['Source 1 n'],'_{-', '%.2f'%self.Ldic['Source 1 n'],'}^{+', '%.2f'%self.Udic['Source 1 n'], '}$ & $','%.2f'%self.Ddic['Source 1 pa'],'_{-', '%.2f'%self.Ldic['Source 1 pa'],'}^{+', '%.2f'%self.Udic['Source 1 pa'], '}$ & $','%.2f'%self.Ddic['Source 1 q'],'_{-', '%.2f'%self.Ldic['Source 1 q'],'}^{+', '%.2f'%self.Udic['Source 1 q'], '}$',r'\\'
        ###
        if self.srcno ==2:
            print 'source 2 & $', '%.2f'%(self.Ddic['Source 2 x']+self.hdic['Lens 1 x']), '_{-', '%.2f'%self.Ldic['Source 2 x'],'}^{+','%.2f'%self.Udic['Source 2 x'], '}$ & $', '%.2f'%(self.Ddic['Source 2 y']+self.hdic['Lens 1 x']),'_{-', '%.2f'%self.Ldic['Source 2 y'],'}^{+', '%.2f'%self.Udic['Source 2 y'], '}$ & $', '%.2f'%self.Ddic['Source 2 re'],'_{-', '%.2f'%self.Ldic['Source 2 re'],'}^{+', '%.2f'%self.Udic['Source 2 re'], '}$ & $', '%.2f'%self.Ddic['Source 2 n'],'_{-', '%.2f'%self.Ldic['Source 2 n'],'}^{+', '%.2f'%self.Udic['Source 2 n'], '}$ & $','%.2f'%self.Ddic['Source 2 pa'],'_{-', '%.2f'%self.Ldic['Source 2 pa'],'}^{+', '%.2f'%self.Udic['Source 2 pa'], '}$ & $','%.2f'%self.Ddic['Source 2 q'],'_{-', '%.2f'%self.Ldic['Source 2 q'],'}^{+', '%.2f'%self.Udic['Source 2 q'], '}$',r'\\'
        ###
        if 'Galaxy 1 re' not in self.hdic.keys():
            print 'galaxy 1 & $', '%.2f'%self.Ddic['Galaxy 1 x'], '_{-', '%.2f'%self.Ldic['Galaxy 1 x'],'}^{+','%.2f'%self.Udic['Galaxy 1 x'], '}$ & $', '%.2f'%self.Ddic['Galaxy 1 y'],'_{-', '%.2f'%self.Ldic['Galaxy 1 y'],'}^{+', '%.2f'%self.Udic['Galaxy 1 y'], '}$ & $', '%.2f'%self.Ddic['Galaxy 1 re'],'_{-', '%.2f'%self.Ldic['Galaxy 1 re'],'}^{+', '%.2f'%self.Udic['Galaxy 1 re'], '}$ & $', '%.2f'%self.Ddic['Galaxy 1 n'],'_{-', '%.2f'%self.Ldic['Galaxy 1 n'],'}^{+', '%.2f'%self.Udic['Galaxy 1 n'], '}$ & $','%.2f'%self.Ddic['Galaxy 1 pa'],'_{-', '%.2f'%self.Ldic['Galaxy 1 pa'],'}^{+', '%.2f'%self.Udic['Galaxy 1 pa'], '}$ & $','%.2f'%self.Ddic['Galaxy 1 q'],'_{-', '%.2f'%self.Ldic['Galaxy 1 q'],'}^{+', '%.2f'%self.Udic['Galaxy 1 q'], '}$',r'\\'
            ###
            print 'galaxy 2 & $', '%.2f'%self.Ddic['Galaxy 2 x'], '_{-', '%.2f'%self.Ldic['Galaxy 2 x'],'}^{+','%.2f'%self.Udic['Galaxy 2 x'], '}$ & $', '%.2f'%self.Ddic['Galaxy 2 y'],'_{-', '%.2f'%self.Ldic['Galaxy 2 y'],'}^{+', '%.2f'%self.Udic['Galaxy 2 y'], '}$ & $', '%.2f'%self.Ddic['Galaxy 2 re'],'_{-', '%.2f'%self.Ldic['Galaxy 2 re'],'}^{+', '%.2f'%self.Udic['Galaxy 2 re'], '}$ & $', '%.2f'%self.Ddic['Galaxy 2 n'],'_{-', '%.2f'%self.Ldic['Galaxy 2 n'],'}^{+', '%.2f'%self.Udic['Galaxy 2 n'], '}$ & $','%.2f'%self.Ddic['Galaxy 2 pa'],'_{-', '%.2f'%self.Ldic['Galaxy 2 pa'],'}^{+', '%.2f'%self.Udic['Galaxy 2 pa'], '}$ & $','%.2f'%self.Ddic['Galaxy 2 q'],'_{-', '%.2f'%self.Ldic['Galaxy 2 q'],'}^{+', '%.2f'%self.Udic['Galaxy 2 q'], '}$',r'\\'
            ###
        print r'\end{tabular}'
        print r"\caption{'Source properties in the K band'}$}"
        print r'\end{table}'
        ### print PSF model
        print r'\begin{table}[H]'
        print r'\centering'
        print r'\begin{tabular}{|cccc|}\hline'
        print r'$\sigma$ & $q$ & $pa$ & amp \\\hline'
        print '%.2f'%self.Ddic['sigma 1'], '&', '%.2f'%self.Ddic['q 1'],'&','%.2f'%self.Ddic['pa 1'],'&','%.2f'%self.Ddic['amp 1'], r'\\'
        print '%.2f'%self.Ddic['sigma 2'], '&', '%.2f'%self.Ddic['q 2'],'&','%.2f'%self.Ddic['pa 2'],'&','%.2f'%self.Ddic['amp 2'], r'\\'
        print '%.2f'%self.Ddic['sigma 3'], '&', '%.2f'%self.Ddic['q 3'],'&','%.2f'%self.Ddic['pa 3'],'&','%.2f'%self.Ddic['amp 3'], r'\\\hline'
        print r'\end{tabular}'
        print r'\caption{PSF components}'
        print r'\end{table}'

    def BuildSources(self):
        self.srcs = []
        for number in range(1,1+self.srcno):
            name = 'Source '+str(number)
            p = {}
            for key in 'q','re','n','pa':
                p[key] = self.Ddic[name+' '+key]
            for key in 'x','y':
                p[key] = self.Ddic[name+' '+key]+self.hdic['Lens 1 '+key]
            self.srcs.append(SBBModels.Sersic(name,p))

    def BuildGalaxies(self):
        self.gals = []
        for number in range(1,1+self.galno):
            name = 'Galaxy '+str(number)
            p = {}
            if 'Galaxy 1 re' in self.hdic.keys():
                for key in 'x','y','q','re','n','pa':
                    p[key] = self.hdic[name+' '+key]
            else: 
                for key in 'x','y','q','re','n','pa':
                    p[key] = self.Ddic[name+' '+key]
            self.gals.append(SBBModels.Sersic(name,p))

    def BuildLenses(self):
        self.lenses = []
        p = {}
        for key in 'x','y','q','pa','b','eta':
            p[key] = self.hdic['Lens 1 '+key]
        self.lenses.append(MassModels.PowerLaw('Lens 1',p))
        p = {}
        p['x'] = self.lenses[0].pars['x']
        p['y'] = self.lenses[0].pars['y']
        p['b'] = self.hdic['extShear']
        p['pa'] = self.hdic['extShear PA']
        self.lenses.append(MassModels.ExtShear('shear',p))

    def AddImages(self,img,sig,psf=None,Dx=None,Dy=None):
        self.img=img
        self.sig=sig
        self.psf=psf
        self.Dx=Dx
        self.Dy=Dy

    def EasyAddImages(self):
        if self.name == 'J0837':
            self.img,self.Dx,self.Dy,self.OVRS,self.pix = EELsImages.J0837_K()
        elif self.name == 'J0901':
            self.img,self.Dx,self.Dy,self.OVRS,self.pix = EELsImages.J0901_K()
        elif self.name == 'J0913':
            self.img,self.Dx,self.Dy,self.OVRS,self.pix = EELsImages.J0913_K(self.srcno)
        elif self.name == 'J1125':
            self.img,self.Dx,self.Dy,self.OVRS,self.pix = EELsImages.J1125_K()
        elif self.name == 'J1144':
            self.img,self.Dx,self.Dy,self.OVRS,self.pix  = EELsImages.J1144_K()
        elif self.name == 'J1218':
            self.img,self.Dx,self.Dy,self.OVRS,self.pix = EELsImages.J1218_K()
        elif self.name == 'J1248':
            self.img,self.Dx,self.Dy,self.OVRS,self.pix = EELsImages.J1248_K()
        elif self.name == 'J1323':
            self.img,self.Dx,self.Dy,self.OVRS,self.pix = EELsImages.J1323_K()
        elif self.name == 'J1347':
            self.img,self.Dx,self.Dy,self.OVRS,self.pix = EELsImages.J1347_K()
        elif self.name == 'J1446':
            self.img,self.Dx,self.Dy,self.OVRS,self.pix = EELsImages.J1446_K()
        elif self.name == 'J1605':
            self.img,self.Dx,self.Dy,self.OVRS,self.pix = EELsImages.J1605_K()
        elif self.name == 'J1606':
            self.img,self.Dx,self.Dy,self.OVRS,self.pix = EELsImages.J1606_K()
        elif self.name == 'J1619':
            self.img,self.Dx,self.Dy,self.OVRS,self.pix = EELsImages.J1619_K()
        elif self.name == 'J2228':
            self.img,self.Dx,self.Dy,self.OVRS,self.pix = EELsImages.J2228_K()
        else:
            print 'you discovered a new EEL?!?'
        self.sig = np.ones(self.img.shape)
        

    def GetFits(self,plotsep=False,plotresid=False,plotcomps=False):
        OVRS=self.OVRS
        yo,xo = iT.coords(self.img.shape)
        yc,xc=iT.overSample(self.img.shape,OVRS)
        dx = self.Ddic['xoffset']
        dy = self.Ddic['yoffset']
        xp,yp = xc*self.pix+dx+self.Dx,yc*self.pix+dy+self.Dy
        xop,yop = xo*self.pix+dy+self.Dx,yo*self.pix+dy+self.Dy
        image = self.img
        sigma = self.sig
        xpsf,ypsf = iT.coords((81,81))-40
        psfObj1 = SBObjects.Gauss('psf 1',{'x':0,'y':0,'sigma':self.Ddic['sigma 1'],'q':self.Ddic['q 1'],'pa':self.Ddic['pa 1'],'amp':self.Ddic['amp 1']})
        psfObj2 = SBObjects.Gauss('psf 2',{'x':0,'y':0,'sigma':self.Ddic['sigma 2'],'q':self.Ddic['q 2'],'pa':self.Ddic['pa 2'],'amp':self.Ddic['amp 2']})
        psfObj3 = SBObjects.Gauss('psf 3',{'x':0,'y':0,'sigma':self.Ddic['sigma 3'],'q':self.Ddic['q 3'],'pa':self.Ddic['pa 3'],'amp':self.Ddic['amp 3']})
        psf1 = psfObj1.pixeval(xpsf,ypsf)  / (np.pi*2.*self.Ddic['sigma 1']**2.)
        psf2 = psfObj2.pixeval(xpsf,ypsf) / (np.pi*2.*self.Ddic['sigma 2']**2.)
        psf3 = psfObj3.pixeval(xpsf,ypsf)  / (np.pi*2.*self.Ddic['sigma 3']**2.)
        psf = psf1 + psf2 + psf3 
        psf /= psf.sum()
        psf = convolve.convolve(image,psf)[1]
        imin,sigin,xin,yin = image.flatten(), sigma.flatten(),xp.flatten(),yp.flatten()
        n = 0
        model = np.empty(((len(self.gals) + len(self.srcs)+1),imin.size))
        for gal in self.gals:
            gal.setPars()
            tmp = xc*0.
            tmp = gal.pixeval(xp,yp,self.pix/OVRS,csub=23) # evaulate on the oversampled grid. OVRS = number of new pixels per old pixel.
            tmp = iT.resamp(tmp,OVRS,True) # convert it back to original size
            tmp = convolve.convolve(tmp,psf,False)[0]
            model[n] = tmp.ravel()
            n +=1
        for lens in self.lenses:
            lens.setPars()
        x0,y0 = pylens.lens_images(self.lenses,self.srcs,[xp,yp],1./OVRS,getPix=True)
        for src in self.srcs:
            src.setPars()
            tmp = xc*0.
            if 'boxiness' in self.dic.keys() and src.pars['q']<0.2:
                print 'this source has a boxy/disky component with c = ', '%.2f'%self.Ddic['boxiness'], '. I hope this is J1606!'
                tmp = src.boxypixeval(x0,y0,self.pix/OVRS,csub=23,c=self.Ddic['boxiness'])
            else:
                tmp = src.pixeval(x0,y0,self.pix/OVRS,csub=23)
            tmp = iT.resamp(tmp,OVRS,True)
            tmp = convolve.convolve(tmp,psf,False)[0]
            model[n] = tmp.ravel()
            if self.name == 'J0837':
                if src.name == 'Source 2':
                    print "just making J0837's second source component into a dust lane..."
                    model[n] *= -1
            n +=1
        model[n] = np.ones(model[n].shape)
        n +=1
        rhs = (imin/sigin) # data
        op = (model/sigin).T # model matrix
        fit, chi = optimize.nnls(op,rhs)
        components = (model.T*fit).T.reshape((n,image.shape[0],image.shape[1]))
        model = components.sum(0)
        if plotsep:
            SotPleparately(image,model,sigma,'K')
        if plotresid:
            NotPlicely(image,model,sigma,'K')
        if plotcomps:
            CotSomponents(components,'K')
        self.fits = fit
        self.model = model
        self.components = components
        print self.components

    def Initialise(self):
        self.MakeDict()
        self.BuildSources()
        self.BuildGalaxies()
        self.BuildLenses()
        self.EasyAddImages()
        self.GetFits(plotresid=False)


    def GalaxySize(self,kpc=False):
        self.z = lens_redshifts[self.name]
        self.Da = astCalc.da(self.z)
        self.scale = self.Da*1e3*np.pi/180./3600.
        if self.galno == 1:
            self.Re_k = self.hdic['Galaxy 1 re']*0.05
        elif self.galno == 2:
            Xgrid = np.logspace(-4,5,1501)
            Res = []
            galaxy = self.fits[0]*self.gals[0].eval(Xgrid) + self.fits[1]*self.gals[1].eval(Xgrid)
            R = Xgrid.copy()
            light = galaxy*2.*np.pi*R
            mod = splrep(R,light,t=np.logspace(-3.8,4.8,1301))
            intlight = np.zeros(len(R))
            for i in range(len(R)):
                intlight[i] = splint(0,R[i],mod)
            model = splrep(intlight[:-300],R[:-300])
            if len(model[1][np.where(np.isnan(model[1])==True)]>0):
                print "arrays need to be increasing monotonically! But don't worry about it"
                model = splrep(intlight[:-450],R[:-450])
            if len(model[1][np.where(np.isnan(model[1])==True)]>0):
                print "take two: arrays need to be increasing monotonically! But don't worry about it"
                model = splrep(intlight[:-600],R[:-600])
            if len(model[1][np.where(np.isnan(model[1])==True)]>0):
                print "take three: arrays need to be increasing monotonically! But don't worry about it"
                model = splrep(intlight[:-700],R[:-700])
            reff = splev(0.5*intlight[-1],model)
            self.Re_k = reff*0.05
        else:
            print 'do you have more than two galaxy components? and do you NEED them really?'
        if kpc:
            return self.Re_k*self.scale
        return self.Re_k

    def GalaxyObsMag(self):
        if self.galno == 1:
            self.mag_k = self.gals[0].getMag(self.fits[0],self.ZP)
        elif self.galno == 2:
            if np.any(self.fits[0:2]==0):
                print ' something is zero!'
                ii = np.where(self.fits[0:2] !=0)[0]
                self.mag_k = self.gals[ii].getMag(self.fits[ii],self.ZP)
            else:
                mk1,mk2 = self.gals[0].getMag(self.fits[0],self.ZP), self.gals[1].getMag(self.fits[1],self.ZP)
                Fk = 10**(0.4*(self.ZP-mk1)) + 10**(0.4*(self.ZP-mk2))
                self.mag_k = -2.5*np.log10(Fk) + self.ZP
        else:
            print 'you have more than two galaxy components?!'
        return self.mag_k


    def GetSB(self):
        self.mu_k = self.mag_k + 2.5*np.log10(2.*np.pi*self.Re_k**2.)
        return self.mu_k

    def GetRestSB(self):
        self.mu_k_rest = self.mag_k_rest + 2.5*np.log10(2.*np.pi*self.Re_k**2.)
        return self.mu_k_rest

    def GetPhotometry(self):
        self.Lk,self.mag_k_rest = MassK(self.mag_k, self.z)
        return [self.mag_k_rest, self.Lk]

    def GetPhotometryAtAge(self,age=4.):
       self.Lk,self.mag_rest = MassK(self.mag_k, self.z, age=age)
       return self.Lk
        
    def Arcsec2Kpc(self,z=None):
        self.z=z
        self.Da = astCalc.da(self.z)
        self.scale = self.Da*1e3*np.pi/180./3600.
        self.Re_k_kpc = self.Re_k*scale
        return self.Re_k_kpc


    def MakePDFDict(self):
        self.dictionaries = []
        for b1,b2 in itertools.product(range(self.trace.shape[0]), range(self.trace.shape[1])):
            go = []
            for key in self.dic.keys():
                go.append((key,self.dic[key][b1,b2]))
            go = dict(go)
            self.dictionaries.append(go)
        print len(self.dictionaries)

    def GetPDFs(self,kpc=False):
        self.muPDF = []
        self.murestPDF = []
        self.RePDF = []
        self.magrestPDF = []
        self.LumPDF = []
        self.magPDF = []
        self.fitPDF = []
        OVRS=self.OVRS
        yo,xo = iT.coords(self.img.shape)
        yc,xc=iT.overSample(self.img.shape,OVRS)
        image, sigma = self.img, self.sig
        xpsf,ypsf = iT.coords((81,81))-40
        for b1 in range(0,len(self.dictionaries),50):
            if 'Galaxy 1 re' not in self.hdic.keys():
                gals = []
                for number in range(1,1+self.galno):
                    name = 'Galaxy '+str(number)
                    p = {}
                    for key in 'x','y','q','re','n','pa':
                        p[key] = self.dictionaries[b1][name+' '+key]
                    gals.append(SBBModels.Sersic(name,p))
            else:
                gals = self.gals
            srcs = []
            for number in range(1,1+self.srcno):
                name = 'Source '+str(number)
                p = {}
                for key in 'q','re','n','pa':
                    p[key] = self.dictionaries[b1][name+' '+key]
                for key in 'x','y':
                    p[key] = self.dictionaries[b1][name+' '+key]+self.hdic['Lens 1 '+key]
                srcs.append(SBBModels.Sersic(name,p))
            # fits
            dx = self.dictionaries[b1]['xoffset']
            dy = self.dictionaries[b1]['yoffset']
            xp,yp = xc*self.pix+dx+self.Dx,yc*self.pix+dy+self.Dy
            xop,yop = xo*self.pix+dy+self.Dx,yo*self.pix+dy+self.Dy
            psfObj1 = SBObjects.Gauss('psf 1',{'x':0,'y':0,'sigma':self.dictionaries[b1]['sigma 1'],'q':self.dictionaries[b1]['q 1'],'pa':self.dictionaries[b1]['pa 1'],'amp':self.dictionaries[b1]['amp 1']})
            psfObj2 = SBObjects.Gauss('psf 2',{'x':0,'y':0,'sigma':self.dictionaries[b1]['sigma 2'],'q':self.dictionaries[b1]['q 2'],'pa':self.dictionaries[b1]['pa 2'],'amp':self.dictionaries[b1]['amp 2']})
            psfObj3 = SBObjects.Gauss('psf 3',{'x':0,'y':0,'sigma':self.dictionaries[b1]['sigma 3'],'q':self.dictionaries[b1]['q 3'],'pa':self.dictionaries[b1]['pa 3'],'amp':self.dictionaries[b1]['amp 3']})
            psf1 = psfObj1.pixeval(xpsf,ypsf)  / (np.pi*2.*self.dictionaries[b1]['sigma 1']**2.)
            psf2 = psfObj2.pixeval(xpsf,ypsf) / (np.pi*2.*self.dictionaries[b1]['sigma 2']**2.)
            psf3 = psfObj3.pixeval(xpsf,ypsf)  / (np.pi*2.*self.dictionaries[b1]['sigma 3']**2.)
            psf = psf1 + psf2 + psf3 
            psf /= psf.sum()
            psf = convolve.convolve(image,psf)[1]
            imin,sigin,xin,yin = image.flatten(), sigma.flatten(),xp.flatten(),yp.flatten()
            n = 0
            model = np.empty(((len(gals) + len(srcs)+1),imin.size))
            for gal in gals:
                gal.setPars()
                tmp = xc*0.
                tmp = gal.pixeval(xp,yp,self.pix/OVRS,csub=11) 
                tmp = iT.resamp(tmp,OVRS,True) 
                tmp = convolve.convolve(tmp,psf,False)[0]
                model[n] = tmp.ravel()
                n +=1
            for lens in self.lenses:
                lens.setPars()
            x0,y0 = pylens.lens_images(self.lenses,srcs,[xp,yp],1./OVRS,getPix=True)
            for src in srcs:
                src.setPars()
                tmp = xc*0.
                if 'boxiness' in self.dic.keys() and src.name == 'Source 2':
                    tmp = src.boxypixeval(x0,y0,self.pix/OVRS,csub=11,c=self.dictionaries[b1]['boxiness']) 
                else:
                    tmp = src.pixeval(x0,y0,self.pix/OVRS,csub=11)
                tmp = iT.resamp(tmp,OVRS,True)
                tmp = convolve.convolve(tmp,psf,False)[0]
                model[n] = tmp.ravel()
                if self.name == 'J0837' and src.name == 'Source 2':
                   model[n] *= -1
                n +=1
            model[n] = np.ones(model[n].shape)
            n +=1
            rhs = (imin/sigin) # data
            op = (model/sigin).T # model matrix
            fit, chi = optimize.nnls(op,rhs)
            # source plane (observed) magnitudes
            if self.galno == 1:
                mag = gals[0].getMag(fit[0],self.ZP)
            elif self.galno == 2:
                if np.any(fit[0:2]==0):
                    ii = np.where(fit[0:2] !=0)[0]
                    mag = gals[ii].getMag(fit[ii],self.ZP)
                else:
                    mk1,mk2 = gals[0].getMag(fit[0],self.ZP), gals[1].getMag(fit[1],self.ZP)
                    Fk = 10**(0.4*(self.ZP-mk1)) + 10**(0.4*(self.ZP-mk2))
                    mag = -2.5*np.log10(Fk) + self.ZP
            # intrinsic magnitudes
            Lk,k_rest = MassK(mag,self.z)
            # sizes
            if self.galno == 1:
                Re = gals[0].pars['re']*0.05
            elif self.galno == 2:
                Xgrid = np.logspace(-4,5,1501)
                Ygrid = np.logspace(-4,5,1501)
                galaxy = fit[0]*gals[0].eval(Xgrid) + fit[1]*gals[1].eval(Xgrid)
                R = Xgrid.copy()
                light = galaxy*2.*np.pi*R
                mod = splrep(R,light,t=np.logspace(-3.8,4.8,1301))
                intlight = np.zeros(len(R))
                for i in range(len(R)):
                    intlight[i] = splint(0,R[i],mod)
                model = splrep(intlight[:-300],R[:-300])
                if len(model[1][np.where(np.isnan(model[1])==True)]>0):
                    model = splrep(intlight[:-600],R[:-600])
                if len(model[1][np.where(np.isnan(model[1])==True)]>0):
                    model = splrep(intlight[:-700],R[:-700])
                reff = splev(0.5*intlight[-1],model)
                Re = reff*0.05
            # surface brightnesses
            mu_rest = k_rest +  2.5*np.log10(2.*np.pi*Re**2.)
            mu = mag +  2.5*np.log10(2.*np.pi*Re**2.)
            #print mag
            self.muPDF.append(mu)
            self.murestPDF.append(mu_rest)
	    if kpc:
		self.RePDF.append(Re*self.scale)
	    else:
		self.RePDF.append(Re)
            self.magrestPDF.append(k_rest)
            self.LumPDF.append(Lk)
            self.magPDF.append(mag)
            self.fitPDF.append(fit)

        def EasyAddPDFs(self,mod):
            try:
                self.muPDF, self.magPDF,self.RePDF = np.load('/data/ljo31/Lens/PDFs/'+str(mod)+'_lensgals_Kp_'+str(self.name+'.npy'))
            except:
                print 'must be 211,112,212 or some such combination'

    def PrintPDFTable(self):
        magPDF,magrestPDF,muPDF,RePDF,LumPDF = np.array(self.magPDF),np.array(self.magrestPDF),np.array(self.muPDF),np.array(self.RePDF),np.array(self.LumPDF)
        print r'\begin{table}[H]'
        print r'\centering'
        print r'\begin{tabular}{|c|cccccc|}\hline'
        print r' band & $m_{int}$ & $m_{rest}$ & $\mu_{e}$ & $R_e$ / ', r"$''$ & $R_e$ / kpc & $\log(L)$ ",r'\\\hline'
        print 'K', '& $', '%.2f'%(np.percentile(magPDF,50)), '_{-','%.2f'%(np.percentile(magPDF,50)-np.percentile(magPDF,16)), '}^{+', '%.2f'%(np.percentile(magPDF,84)-np.percentile(magPDF,50)), '}$ & $', '%.2f'%(np.percentile(magrestPDF,50)), '_{-','%.2f'%(np.percentile(magrestPDF,50)-np.percentile(magrestPDF,16)), '}^{+', '%.2f'%(np.percentile(magrestPDF,84)-np.percentile(magrestPDF,50)), '}$ & $', '%.2f'%(np.percentile(muPDF,50)), '_{-','%.2f'%(np.percentile(muPDF,50)-np.percentile(muPDF,16)), '}^{+', '%.2f'%(np.percentile(muPDF,84)-np.percentile(muPDF,50)), '}$ & $', '%.2f'%(np.percentile(RePDF,50)), '_{-','%.2f'%(np.percentile(RePDF,50)-np.percentile(RePDF,16)), '}^{+', '%.2f'%(np.percentile(RePDF,84)-np.percentile(RePDF,50)), '}$ & $', '%.2f'%(np.percentile(RePDF,50)*self.scale), '_{-','%.2f'%(self.scale*(np.percentile(RePDF,50)-np.percentile(RePDF,16))), '}^{+', '%.2f'%(self.scale*(np.percentile(RePDF,84)-np.percentile(RePDF,50))), '}$ & $', '%.2f'%(np.percentile(LumPDF,50)), '_{-','%.2f'%(np.percentile(LumPDF,50)-np.percentile(LumPDF,16)), '}^{+', '%.2f'%(np.percentile(LumPDF,84)-np.percentile(LumPDF,50)), r'}$ \\\hline'
        print r'\end{tabular}'
        print r'\end{table}'


    def UncertaintiesFromPDF(self,makecat=False):
        magPDF,muPDF,RePDF,= np.array(self.magPDF),np.array(self.muPDF),np.array(self.RePDF)
        keys = ['magPDF','muPDF','RePDF']
        vals = [magPDF,muPDF,RePDF]
        l,d,u=[],[],[]
        f1,f2,f3 = np.zeros(len(keys)),np.zeros(len(keys)),np.zeros(len(keys))
        for i in range(len(vals)):
            key = keys[i]
            f = vals[i]
            K1,K2,K3 = np.percentile(f,50)-np.percentile(f,16), np.percentile(f,50), np.percentile(f,84) - np.percentile(f,50)
            d.append((key,K2))
            l.append((key,K1))
            u.append((key,K3))
            f1[i],f2[i],f3[i] = [K2,K1,K3]
        self.photometry = [l,d,u]
        if makecat:
            return f1.flatten().tolist(),f2.flatten().tolist(),f3.flatten().tolist()
        return d,l,u


