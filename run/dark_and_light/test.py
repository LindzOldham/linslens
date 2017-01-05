import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as tic

fig = plt.figure(tight_layout=False)

x = np.arange(100)
y = 3.*np.sin(x*2.*np.pi/100.)

for i in range(2):
    temp = 211 + i
    print temp
    ax = plt.subplot(temp)
    plt.plot(x,y)
    plt.subplots_adjust(hspace = 0)
    temp = tic.MaxNLocator(3)
    ax.yaxis.set_major_locator(temp)
    ax.set_xticklabels(())
    ax.title.set_visible(False)

plt.show()

"""
To create plots that share a common axes (visually) you can set the
hspace between the subplots close to zero (do not use zero itself).
Normally you'll want to turn off the tick labels on all but one of the
axes.

In this example the plots share a common xaxis but you can follow the
same logic to supply a common y axis.
"""
import matplotlib.pyplot as plt
import numpy as np

t = np.arange(0.0, 2.0, 0.01)

s1 = np.sin(2*np.pi*t)
s2 = np.exp(-t)
s3 = s1*s2

# axes rect in relative 0,1 coords left, bottom, width, height.  Turn
# off xtick labels on all but the lower plot


f = plt.figure()
plt.subplots_adjust(hspace=0.001)


ax1 = plt.subplot(311)
ax1.plot(t, s1)
plt.yticks(np.arange(-0.9, 1.0, 0.4))
plt.ylim(-1, 1)

ax2 = plt.subplot(312, sharex=ax1)
ax2.plot(t, s2)
plt.yticks(np.arange(0.1, 1.0, 0.2))
plt.ylim(0, 1)

ax3 = plt.subplot(313, sharex=ax1)
ax3.plot(t, s3)
plt.yticks(np.arange(-0.9, 1.0, 0.4))
plt.ylim(-1, 1)

xticklabels = ax1.get_xticklabels() + ax2.get_xticklabels()
plt.setp(xticklabels, visible=False)

plt.show()
