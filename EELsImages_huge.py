import numpy as np, pylab as pl, pyfits as py

def J0837():
    img1 = py.open('/data/ljo31/Lens/J0837/F606W_sci_cutout_huge.fits')[0].data.copy()#[30:-30,30:-30] 
    sig1 = py.open('/data/ljo31/Lens/J0837/F606W_noisemap_huge.fits')[0].data.copy()#[30:-30,30:-30] 
    psf1 = py.open('/data/ljo31/Lens/J0837/F606W_psf1.fits')[0].data.copy()
    psf1 = psf1/np.sum(psf1)
    img2 = py.open('/data/ljo31/Lens/J0837/F814W_sci_cutout_huge.fits')[0].data.copy()#[30:-30,30:-30]
    sig2 = py.open('/data/ljo31/Lens/J0837/F814W_noisemap_huge.fits')[0].data.copy()#[30:-30,30:-30] 
    psf2 = py.open('/data/ljo31/Lens/J0837/F814W_psf3.fits')[0].data.copy()
    psf2 = psf2/np.sum(psf2)
    Dx,Dy = -100,-100
    OVRS=1
    mask = py.open('/data/ljo31/Lens/J0837/mask_huge.fits')[0].data
    return img1,sig1,psf1,img2,sig2,psf2,Dx,Dy,OVRS,mask

def J0901():
    img1 = py.open('/data/ljo31/Lens/J0901/F606W_sci_cutout_huge.fits')[0].data.copy()[25:-25,25:-25]
    sig1 = py.open('/data/ljo31/Lens/J0901/F606W_noisemap_huge.fits')[0].data.copy()[25:-25,25:-25]
    psf1 = py.open('/data/ljo31/Lens/J0901/F606W_psf2.fits')[0].data.copy()
    psf1 = psf1[5:-6,5:-6]
    psf1 = psf1/np.sum(psf1)
    img2 = py.open('/data/ljo31/Lens/J0901/F814W_sci_cutout_huge.fits')[0].data.copy()[25:-25,25:-25]
    sig2 = py.open('/data/ljo31/Lens/J0901/F814W_noisemap_huge.fits')[0].data.copy()[25:-25,25:-25]
    psf2 = py.open('/data/ljo31/Lens/J0901/F814W_psf2.fits')[0].data.copy()
    psf2 = psf2[4:-6,3:-6]
    psf2 = psf2/np.sum(psf2)
    Dx,Dy = -85.,-85.
    OVRS=2
    mask =  py.open('/data/ljo31/Lens/J0901/mask_huge.fits')[0].data[25:-25,25:-25]
    return img1,sig1,psf1,img2,sig2,psf2,Dx,Dy,OVRS,mask

def J0913(srcno):
    img1 = py.open('/data/ljo31/Lens/J0913/F555W_sci_cutout_huge2.fits')[0].data.copy()
    sig1 = py.open('/data/ljo31/Lens/J0913/F555W_noisemap_huge2.fits')[0].data.copy()
    psf1 = py.open('/data/ljo31/Lens/J0913/F555W_psf4.fits')[0].data.copy()
    psf1 = psf1/np.sum(psf1)
    img2 = py.open('/data/ljo31/Lens/J0913/F814W_sci_cutout_huge2.fits')[0].data.copy()
    sig2 = py.open('/data/ljo31/Lens/J0913/F814W_noisemap_huge2.fits')[0].data.copy()
    psf2 = py.open('/data/ljo31/Lens/J0913/F814W_psf3.fits')[0].data.copy()
    psf2 = psf2/np.sum(psf2)
    Dx,Dy = -80,-80
    if srcno==2:
        Dx,Dy = -95.,-100.
    OVRS=2
    mask = py.open('/data/ljo31/Lens/J0913/mask_huge2.fits')[0].data
    return img1,sig1,psf1,img2,sig2,psf2,Dx,Dy,OVRS,mask

def J1125():
    img1 = py.open('/data/ljo31/Lens/J1125/F606W_sci_cutout_huge2.fits')[0].data.copy()[50:-50,50:-50]
    sig1 = py.open('/data/ljo31/Lens/J1125/F606W_noisemap_huge2.fits')[0].data.copy()[50:-50,50:-50]
    psf1 = py.open('/data/ljo31/Lens/J1125/F606W_psf3_filledin.fits')[0].data.copy()
    psf1 = psf1[5:-7,5:-6]
    psf1 = psf1/np.sum(psf1)
    img2 = py.open('/data/ljo31/Lens/J1125/F814W_sci_cutout_huge2.fits')[0].data.copy()[50:-50,50:-50]
    sig2 = py.open('/data/ljo31/Lens/J1125/F814W_noisemap_huge2.fits')[0].data.copy()[50:-50,50:-50]
    psf2 = py.open('/data/ljo31/Lens/J1125/F814W_psf3_filledin.fits')[0].data.copy()
    psf2 = psf2[5:-8,5:-6]
    psf2 = psf2/np.sum(psf2)
    Dx,Dy = -85,-85
    OVRS=2
    mask =  py.open('/data/ljo31/Lens/J1125/mask_huge2.fits')[0].data[50:-50,50:-50]
    return img1,sig1,psf1,img2,sig2,psf2,Dx,Dy,OVRS,mask

def J1144():
    img1 = py.open('/data/ljo31/Lens/J1144/F606W_sci_cutout_huge.fits')[0].data.copy()
    sig1 = py.open('/data/ljo31/Lens/J1144/F606W_noisemap_huge.fits')[0].data.copy()
    psf1 = py.open('/data/ljo31/Lens/J1144/F606W_psf1.fits')[0].data.copy()
    psf1 = psf1/np.sum(psf1)
    img2 = py.open('/data/ljo31/Lens/J1144/F814W_sci_cutout_huge.fits')[0].data.copy()
    sig2 = py.open('/data/ljo31/Lens/J1144/F814W_noisemap_huge.fits')[0].data.copy()
    psf2 = py.open('/data/ljo31/Lens/J1144/F814W_psf1.fits')[0].data.copy()
    psf2 = psf2/np.sum(psf2)
    Dx,Dy = -95,-90
    OVRS=1
    mask = py.open('/data/ljo31/Lens/J1144/mask_huge.fits')[0].data
    return img1,sig1,psf1,img2,sig2,psf2,Dx,Dy,OVRS,mask

def J1218():
    img1 = py.open('/data/ljo31/Lens/J1218/F606W_sci_cutout_huge.fits')[0].data.copy()#[10:-10,10:-10]
    sig1 = py.open('/data/ljo31/Lens/J1218/F606W_noisemap_huge.fits')[0].data.copy()#[1#0:-10,10:-10]
    psf1 = py.open('/data/ljo31/Lens/J1218/F606W_psf1.fits')[0].data.copy()
    psf1 = psf1/np.sum(psf1)
    img2 = py.open('/data/ljo31/Lens/J1218/F814W_sci_cutout_huge.fits')[0].data.copy()#[10:-10,10:-10]
    sig2 = py.open('/data/ljo31/Lens/J1218/F814W_noisemap_huge.fits')[0].data.copy()#[#10:-10,10:-10]
    psf2 = py.open('/data/ljo31/Lens/J1218/F814W_psf1.fits')[0].data.copy()
    psf2 = psf2/np.sum(psf2)
    Dx,Dy = -160,-160.
    OVRS=1
    mask = py.open('/data/ljo31/Lens/J1218/mask_huge.fits')[0].data
    return img1,sig1,psf1,img2,sig2,psf2,Dx,Dy,OVRS,mask

def J1248():
    img1 = py.open('/data/ljo31/Lens/J1248/F555W_sci_cutout.fits')[0].data.copy()[10:-10,20:-25]
    sig1 = py.open('/data/ljo31/Lens/J1248/F555W_noisemap.fits')[0].data.copy()[10:-10,20:-25]
    psf1 = py.open('/data/ljo31/Lens/J1248/F555W_psf1.fits')[0].data.copy()
    psf1 = psf1/np.sum(psf1)
    img2 = py.open('/data/ljo31/Lens/J1248/F814W_sci_cutout.fits')[0].data.copy()[10:-10,20:-25]
    sig2 = py.open('/data/ljo31/Lens/J1248/F814W_noisemap.fits')[0].data.copy()[10:-10,20:-25]
    psf2 = py.open('/data/ljo31/Lens/J1248/psf1_nopedestal.fits')[0].data.copy()[8:-8,8:-8]
    psf2 = psf2/np.sum(psf2)
    Dx,Dy=10.,0
    OVRS=1
    return img1,sig1,psf1,Dx,Dy,OVRS

def J1323(srcno):
    img1 = py.open('/data/ljo31/Lens/J1323/SDSSJ1323+3946_F555W_sci_cutout_huge.fits')[0].data.copy()
    sig1 = py.open('/data/ljo31/Lens/J1323/SDSSJ1323+3946_F555W_noisemap_huge.fits')[0].data.copy()
    psf1 = py.open('/data/ljo31/Lens/J1323/SDSSJ1323+3946_F555W_psf2.fits')[0].data.copy()
    psf1 = psf1[10:-10,11:-10] 
    psf1 = psf1/np.sum(psf1)
    img2 = py.open('/data/ljo31/Lens/J1323/SDSSJ1323+3946_F814W_sci_cutout_huge.fits')[0].data.copy()
    sig2 = py.open('/data/ljo31/Lens/J1323/SDSSJ1323+3946_F814W_noisemap_huge.fits')[0].data.copy()
    psf2 = py.open('/data/ljo31/Lens/J1323/SDSSJ1323+3946_F814W_psf3.fits')[0].data.copy()
    if srcno ==2:
        psf2 = psf2[8:-8,9:-8]
    psf2 = psf2/np.sum(psf2)
    Dx,Dy = -100.,-100.
    OVRS=2
    mask = py.open('/data/ljo31/Lens/J1323/mask_huge.fits')[0].data
    return img1,sig1,psf1,img2,sig2,psf2,Dx,Dy,OVRS,mask

def J1347():
    img1 = py.open('/data/ljo31/Lens/J1347/SDSSJ1347-0101_F606W_sci_cutout_huge.fits')[0].data.copy()
    sig1 = py.open('/data/ljo31/Lens/J1347/SDSSJ1347-0101_F606W_noisemap_huge.fits')[0].data.copy()
    psf1 = py.open('/data/ljo31/Lens/J1347/SDSSJ1347-0101_F606W_psf.fits')[0].data.copy()
    psf1 = psf1[15:-15,15:-15]
    psf1 /= psf1.sum()
    img2 = py.open('/data/ljo31/Lens/J1347/SDSSJ1347-0101_F814W_sci_cutout_huge.fits')[0].data.copy()
    sig2 = py.open('/data/ljo31/Lens/J1347/SDSSJ1347-0101_F814W_noisemap_huge.fits')[0].data.copy()
    psf2 = py.open('/data/ljo31/Lens/J1347/SDSSJ1347-0101_F814W_psf_#2.fits')[0].data.copy()
    psf2 = psf2[15:-15,15:-16]
    psf2 /= psf2.sum()
    Dx,Dy = -89,-87
    OVRS=2
    mask = py.open('/data/ljo31/Lens/J1347/mask_huge.fits')[0].data
    return img1,sig1,psf1,img2,sig2,psf2,Dx,Dy,OVRS,mask

def J1446():
    img1 = py.open('/data/ljo31/Lens/J1446/F606W_sci_cutout_huge.fits')[0].data.copy()
    sig1 = py.open('/data/ljo31/Lens/J1446/F606W_noisemap_huge.fits')[0].data.copy()
    psf1 = py.open('/data/ljo31/Lens/J1446/F606W_psf1.fits')[0].data.copy()
    psf1 = psf1[5:-5,5:-5]
    psf1 = psf1/np.sum(psf1)
    img2 = py.open('/data/ljo31/Lens/J1446/F814W_sci_cutout_huge.fits')[0].data.copy()
    sig2 = py.open('/data/ljo31/Lens/J1446/F814W_noisemap_huge.fits')[0].data.copy()
    psf2 = py.open('/data/ljo31/Lens/J1446/F814W_psf1.fits')[0].data.copy()
    psf2 = psf2[6:-6,7:-7]
    psf2 = psf2/np.sum(psf2)
    Dx,Dy = -95.,-110.
    OVRS=2
    mask = py.open('/data/ljo31/Lens/J1446/mask_huge.fits')[0].data
    return img1,sig1,psf1,img2,sig2,psf2,Dx,Dy,OVRS,mask

def J1605(srcno):
    img1 = py.open('/data/ljo31/Lens/J1605/F555W_sci_cutout_huge.fits')[0].data.copy()
    sig1 = py.open('/data/ljo31/Lens/J1605/F555W_noisemap_huge.fits')[0].data.copy() 
    psf1 = py.open('/data/ljo31/Lens/J1605/SDSSJ1605+3811_F555W_psf.fits')[0].data.copy()
    psf1 = psf1[10:-10,10:-10]
    psf1 = psf1/np.sum(psf1)
    img2 = py.open('/data/ljo31/Lens/J1605/F814W_sci_cutout_huge.fits')[0].data.copy()
    sig2 = py.open('/data/ljo31/Lens/J1605/F814W_noisemap_huge.fits')[0].data.copy()
    if srcno ==1:
        psf2 = py.open('/data/ljo31/Lens/J1605/F814W_psf_#3.fits')[0].data.copy()  
        psf2= psf2[20:-20,19:-20]
    elif srcno ==2:
        psf2 = py.open('/data/ljo31/Lens/J1605/F814W_psf1new.fits')[0].data.copy()  
    psf2 /= psf2.sum()
    Dx,Dy = -115.,-115.
    OVRS=1
    mask = py.open('/data/ljo31/Lens/J1605/mask_huge.fits')[0].data
    return img1,sig1,psf1,img2,sig2,psf2,Dx,Dy,OVRS,mask

def J1606():
    img1 = py.open('/data/ljo31/Lens/J1606/F606W_sci_cutout_huge.fits')[0].data.copy()
    sig1 = py.open('/data/ljo31/Lens/J1606/F606W_noisemap_huge.fits')[0].data.copy()
    psf1 = py.open('/data/ljo31/Lens/J1606/F606W_psf.fits')[0].data.copy()
    psf1 = psf1/np.sum(psf1)
    img2 = py.open('/data/ljo31/Lens/J1606/F814W_sci_cutout_huge.fits')[0].data.copy()
    sig2 = py.open('/data/ljo31/Lens/J1606/F814W_noisemap_huge.fits')[0].data.copy()
    psf2 = py.open('/data/ljo31/Lens/J1606/F814W_psf.fits')[0].data.copy()
    psf2 = psf2/np.sum(psf2)
    Dx,Dy = -145,-135
    OVRS=1
    mask = py.open('/data/ljo31/Lens/J1606/mask_huge.fits')[0].data
    return img1,sig1,psf1,img2,sig2,psf2,Dx,Dy,OVRS,mask

def J1619():
    img1 = py.open('/data/ljo31/Lens/J1619/F606W_sci_cutout_huge.fits')[0].data.copy()
    sig1 = py.open('/data/ljo31/Lens/J1619/F606W_noisemap_huge.fits')[0].data.copy()*2.
    psf1 = py.open('/data/ljo31/Lens/J1619/F606W_psf2new.fits')[0].data.copy()
    psf1 = psf1/np.sum(psf1)
    img2 = py.open('/data/ljo31/Lens/J1619/F814W_sci_cutout_huge.fits')[0].data.copy()
    sig2 = py.open('/data/ljo31/Lens/J1619/F814W_noisemap_huge.fits')[0].data.copy()
    sig2[sig2>0.07] = 0.07
    sig2 *= 2.
    psf2 = py.open('/data/ljo31/Lens/J1619/F814W_psf1neu.fits')[0].data.copy()
    psf2 = psf2/np.sum(psf2)
    psf3 = psf2.copy()
    psf3[5:-5,5:-5] = np.nan
    cond = np.isfinite(psf3)
    m = psf3[cond].mean()
    psf2 = psf2 - m
    psf2 = psf2/np.sum(psf2)
    Dx,Dy = -100,-100
    OVRS=2
    mask = py.open('/data/ljo31/Lens/J1619/mask_huge.fits')[0].data
    return img1,sig1,psf1,img2,sig2,psf2,Dx,Dy,OVRS,mask

def J2228():
    img1 = py.open('/data/ljo31/Lens/J2228/F555W_sci_cutout.fits')[0].data.copy()[45:-45]
    sig1 = py.open('/data/ljo31/Lens/J2228/F606W_noisemap.fits')[0].data.copy()[45:-45]    
    psf1 = py.open('/data/ljo31/Lens/J2228/F606W_psf1.fits')[0].data.copy()
    psf1 = psf1/np.sum(psf1)
    img2 = py.open('/data/ljo31/Lens/J2228/F814W_sci_cutout.fits')[0].data.copy()[45:-45] 
    sig2 = py.open('/data/ljo31/Lens/J2228/F814W_noisemap.fits')[0].data.copy()[45:-45]      
    psf2 = py.open('/data/ljo31/Lens/J2228/F814W_psf1.fits')[0].data.copy()
    psf2 = psf2/np.sum(psf2)
    Dx,Dy = -110,-135
    OVRS=1
    mask = py.open('/data/ljo31/Lens/J2228/mask_huge.fits')[0].data[45:-45]
    return img1,sig1,psf1,img2,sig2,psf2,Dx,Dy,OVRS,mask


### K band data

def J0837_K():
    image = py.open('/data/ljo31/Lens/J0837/J0837_Kp_narrow_med.fits')[0].data.copy()[810:1100,790:1105]
    Dx,Dy = 16.,23.
    OVRS,pix=1,0.2
    return image,Dx,Dy,OVRS,pix

def J0901_K():
    image = py.open('/data/ljo31/Lens/J0901/J0901_Kp_narrow.fits')[0].data.copy()[500:720,510:760] # update this maybe?
    Dx,Dy = 9.,9.
    OVRS,pix=1,0.2
    return image,Dx,Dy,OVRS,pix

def J0913_K(srcno):
    image = py.open('/data/ljo31/Lens/J0913/J0913_nirc2_n_Kp_6x6.fits')[0].data.copy()[141:475,155:460]
    if srcno == 1:
        Dx,Dy = 5.,5.
    elif srcno ==2:
        Dx,Dy = -10.,-15.
    OVRS,pix=1,0.2
    return image,Dx,Dy,OVRS,pix

def J1125_K():
    image = py.open('/data/ljo31/Lens/J1125/Kp_J1125_nirc2_n.fits')[0].data.copy()[650:905,640:915]
    Dx,Dy = 10.,12.
    OVRS,pix=1,0.2
    return image,Dx,Dy,OVRS,pix

def J1144_K():
    image = py.open('/data/ljo31/Lens/J1144/J1144_nirc2_n_Kp_6x6.fits')[0].data.copy()[92:455,130:501]
    Dx,Dy = 8.,7.
    OVRS,pix=1,0.2
    return image,Dx,Dy,OVRS,pix

def J1218_K():
    image = py.open('/data/ljo31/Lens/J1218/J1218_med.fits')[0].data.copy()[640:820,460:640]
    Dx,Dy = 10,0
    OVRS,pix = 1,0.6
    return image,Dx,Dy,OVRS,pix

def J1248_K():
    ()


def J1323_K():
    image = py.open('/data/ljo31/Lens/J1323/J1323_nirc2.fits')[0].data.copy()[615:835,855:1090]
    Dx,Dy = 20.,22.
    OVRS,pix=1,0.2
    return image,Dx,Dy,OVRS,pix

def J1347_K():
    image = py.open('/data/ljo31/Lens/J1347/J1347_med.fits')[0].data.copy()[900:1030,910:1025]
    Dx,Dy = -6.,-5.
    OVRS,pix = 4,0.6
    return image,Dx,Dy,OVRS,pix

def J1446_K():
    image = py.open('/data/ljo31/Lens/J1446/EEL1446_med.fits')[0].data.copy()[730:895,630:920]
    Dx,Dy = 21.,14.
    OVRS,pix=1,0.2
    return image,Dx,Dy,OVRS,pix

def J1605_K():
    image = py.open('/data/ljo31/Lens/J1605/J1605_Kp_narrow_med.fits')[0].data.copy()[535:740,590:835]
    Dx,Dy = 4.,5.
    OVRS,pix=1,0.2
    return image,Dx,Dy,OVRS,pix


def J1606_K():
    image = py.open('/data/ljo31/Lens/J1606/J1606_med.fits')[0].data.copy()[630:830,880:1050]
    Dx,Dy = 3.,-7.
    OVRS,pix=1,0.6
    return image,Dx,Dy,OVRS,pix

def J1619_K():
    image = py.open('/data/ljo31/Lens/J1619/J1619_nirc2_n_Kp_6x6.fits')[0].data.copy()[200:400,200:400]
    Dx,Dy = 26.,29.
    OVRS,pix=1,0.2
    return image,Dx,Dy,OVRS,pix


def J2228_K():
    image = py.open('/data/ljo31/Lens/J2228/J2228_med.fits')[0].data.copy()[550:730,550:730]
    Dx,Dy = 8.,8.
    OVRS,pix = 1,0.6
    return image,Dx,Dy,OVRS,pix


## weight images for unlensed modelling
def J0837_w():
    V,I =  py.open('/data/ljo31/Lens/J0837/F606W_wht_cutout.fits')[0].data[30:-30,30:-30] ,py.open('/data/ljo31/Lens/J0837/F814W_wht_cutout.fits')[0].data[30:-30,30:-30] 
    return V,I

def J0901_w():
    V,I =  py.open('/data/ljo31/Lens/J0901/F606W_wht_cutout.fits')[0].data,py.open('/data/ljo31/Lens/J0901/F814W_wht_cutout.fits')[0].data
    return V,I

def J0913_w():
    V,I =  py.open('/data/ljo31/Lens/J0913/F555W_wht_cutout_big.fits')[0].data,py.open('/data/ljo31/Lens/J0913/F814W_wht_cutout_big.fits')[0].data
    return V,I

def J1125_w():
    V,I =  py.open('/data/ljo31/Lens/J1125/F606W_wht_cutout.fits')[0].data,py.open('/data/ljo31/Lens/J1125/F814W_wht_cutout.fits')[0].data
    return V,I


def J1144_w():
    V,I =  py.open('/data/ljo31/Lens/J1144/F606W_wht_cutout_biggerigger.fits')[0].data,py.open('/data/ljo31/Lens/J1144/F814W_wht_cutout_biggerigger.fits')[0].data
    return V,I

def J1218_w():
    V,I =  py.open('/data/ljo31/Lens/J1218/F606W_wht_cutout.fits')[0].data[10:-10,10:-10],py.open('/data/ljo31/Lens/J1218/F814W_wht_cutout.fits')[0].data[10:-10,10:-10]
    return V,I

def J1248_w():
    ()

def J1323_w():
    V,I =  py.open('/data/ljo31/Lens/J1323/SDSSJ1323+3946_F555W_sci_cutout_big.fits')[0].data,py.open('/data/ljo31/Lens/J1323/SDSSJ1323+3946_F814W_sci_cutout_big.fits')[0].data
    return V,I

def J1347_w():
    V,I =  py.open('/data/ljo31/Lens/J1347/SDSSJ1347-0101_F606W_wht_cutout.fits')[0].data,py.open('/data/ljo31/Lens/J1347/SDSSJ1347-0101_F814W_wht_cutout.fits')[0].data
    return V,I

def J1446_w():
    V,I =  py.open('/data/ljo31/Lens/J1446/F606W_wht_cutout.fits')[0].data,py.open('/data/ljo31/Lens/J1446/F814W_wht_cutout.fits')[0].data
    return V,I

def J1605_w():
    V,I =  py.open('/data/ljo31/Lens/J1605/SDSSJ1605+3811_F555W_wht_cutout2.fits')[0].data,py.open('/data/ljo31/Lens/J1605/SDSSJ1605+3811_F814W_wht_cutout2.fits')[0].data
    return V,I

def J1606_w():
    V,I =  py.open('/data/ljo31/Lens/J1606/F606W_wht_cutout.fits')[0].data[20:-20,20:-40],py.open('/data/ljo31/Lens/J1606/F814W_wht_cutout.fits')[0].data[20:-20,20:-40]
    return V,I

def J1619_w():
    V,I =  py.open('/data/ljo31/Lens/J1619/F606W_wht_cutout.fits')[0].data[30:-38,25:-30],py.open('/data/ljo31/Lens/J1619/F814W_wht_cutout.fits')[0].data[30:-38,25:-30]
    return V,I

def J2228_w():
    V,I =  py.open('/data/ljo31/Lens/J2228/F555W_wht_cutout.fits')[0].data[120:-110,25:-25] ,py.open('/data/ljo31/Lens/J2228/F814W_wht_cutout.fits')[0].data[120:-110,25:-25] 
    return V,I

