#!/usr/bin/python3
# -*- coding: UTF-8 -*-
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
def expDtofloat(s):
    return float(s.decode().replace('D','E'))

data = np.loadtxt('axs-ext.dat', skiprows=15,
                  converters={0:expDtofloat,1:expDtofloat,
                              2:expDtofloat,3:expDtofloat},
                  )
data_size=len(set(data[:,0]))
data = data[:,3].reshape((-1,data_size))

negval = sum(1 for i in data.flatten() if i < 0)
print('Negative', negval, 'of',len(data.flatten()))

data = np.abs(data)

plt.imshow(data,
           origin='lower',
           cmap='jet',
           aspect='auto',
           norm=LogNorm(
               vmin = np.min(data)*1.0,
               vmax = np.max(data)*1.0
               )
           )
plt.show()
