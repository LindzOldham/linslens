import corner_plot
import numpy as np, pylab as pl

result = np.load('/data/ljo31b/EELs/galsub/emceeruns/J0901_parametric_sDPL_3')
lp,trace,dic,_ = result
result2 = np.load('/data/ljo31b/EELs/galsub/emceeruns/J0901_parametric_DPL_3')
lp2,trace2,dic2,_ = result2

'''arr = []
for key in dic.keys():
    dic[key] = dic[key][500:,0].ravel()
for key in dic.keys():
    arr.append(dic[key])
arr=np.array(arr)

labels_1 = ['$R_{Ein}$', '$y_l$', '$n_s$', '$\log r_s$','$M_{\star}$','$y_s$',r'$\theta_s$','$q_l$','$q_s$','$x_l$','SH',r'$\gamma$','$r_{e,s}$','$x_s$','SHPA','LPA']

labels_2 = ['$R_{Ein}$', '$y_l$', '$n_s$', '$\log r_s$','$M_{\star}$','$y_s$',r'$\theta_s$','$q_s$','$x_l$','SH',r'$\gamma$','$r_{e,s}$','$x_s$','SHPA']

#arr[3] = np.log10(arr[3])
corner_plot.corner_plot(arr.T,axis_labels=labels_1)
pl.show()'''

for key in ['Lens 1 b','Lens 1 rs','Lens 1 gamma','stellar mass']:
    pl.figure()
    pl.hist(dic[key][:,0].ravel(),30,alpha=0.5,histtype='stepfilled')
    pl.hist(dic2[key][:,0].ravel(),30,alpha=0.5,histtype='stepfilled')
    pl.title(key)
    
pl.show()

# spherical: rs tends to zero. But it also started around there. 
# Elliptical: rs is still large, but it started there too. It hasn't explored well yet!
