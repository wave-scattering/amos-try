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

import configparser
config = configparser.ConfigParser()
config.read('default.ini')
x_min = float(config['beam']['x_min'])
x_max = float(config['beam']['x_max'])

rl_min = float(config['cylinder']['rl_min'])
rl_max = float(config['cylinder']['rl_max'])
plt.figure()
# spectra = data[0,:]
spectra = data[7,:]
spectra_x = np.linspace(x_min, x_max, len(spectra))
plt.plot(spectra_x, spectra)
spectra = data[-1,:]
plt.plot(spectra_x, spectra)
# plt.title(f'$r/L = ${rl_max}')
plt.tight_layout()
# plt.show()
plt.savefig('spectra.pdf')

plt.figure()
plt.imshow(data,
           origin='lower',
           cmap='hot',
           aspect='auto',
           # vmin = np.mean(data)*0.1, vmax = np.mean(data)*4,
           # extent=(x_min/2/np.pi/50*300, x_max/2/np.pi/50*300, rl_min, rl_max),
           extent=(x_min, x_max, rl_min, rl_max),
           norm=LogNorm(
               vmin = np.min(data)*1.0,
               vmax = np.max(data)*1.0
               )
           )
plt.xlabel(r'$kr$')
plt.ylabel(r'$r/L$')
plt.tight_layout()
# plt.colorbar()
# plt.show()
plt.savefig('map.pdf')
