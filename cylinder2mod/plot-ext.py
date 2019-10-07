#!/usr/bin/python3
# -*- coding: UTF-8 -*-
from matplotlib import pyplot as plt
import numpy as np
def expDtofloat(s):
    return float(s.decode().replace('D','E'))

data = np.loadtxt('axs-ext.dat', skiprows=15,
                  converters={0:expDtofloat,1:expDtofloat,
                              2:expDtofloat,3:expDtofloat},
                  )
data_size=len(set(data[:,0]))
data = data[:,3].reshape((-1,data_size))

plt.imshow(data)
plt.show()
