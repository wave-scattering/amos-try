#!/usr/bin/python3
# -*- coding: UTF-8 -*-
from matplotlib import pyplot as plt
import numpy as np
data = np.loadtxt('index.txt')
data_int = np.loadtxt('index_int.txt')
plt.scatter(data[:,0], data[:,1], label="Re",marker='o', lw=5)
plt.scatter(data[:,0], data[:,2], label="Im", lw=5)
plt.scatter(data_int[:,0], data_int[:,2], label="int Im", lw=1)
plt.scatter(data_int[:,0], data_int[:,1], label="int Re", lw=1)
plt.legend()
plt.show()
