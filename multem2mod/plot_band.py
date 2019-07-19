#!/usr/bin/python3
# -*- coding: UTF-8 -*-
import matplotlib.pyplot as plt
import numpy as np


freq = []
ck_all = []
ck = []
with open('fort.9') as f: 
    for line in f: 
        # print(line[:-1])
        split = line.split()
        if '   ' in line and line[:3] != '   ':
            # print('new line')
            if len(ck) != 0:
                ck_all.append(ck)
                ck = []
            freq.append(float(split[0]))
            split = split[1:]
        for i in range(int(len(split)/2)):
            ck.append(float(split[2*i])+1j*float(split[2*i+1]))
        if len(split) % 2 != 0:
            ck.append(float(split[-1]))
        # print(freq[-1],ck)
    ck_all.append(ck)

print("TOTAL READ: Freq len -> ", len(freq), "  ck_len ->",len(ck_all))
for i in range(len(freq)):
    # print(freq[i], ck_all[i])
    band = np.real(np.array(ck_all[i]))
    # band = np.abs(np.array(ck_all[i]))
    plt.scatter(band,
                [freq[i]#/2/np.pi
                 ]*len(ck_all[i]),
                c='blue',
                s=10
                )
    plt.scatter(-band,
                [freq[i]#/2/np.pi
                 ]*len(ck_all[i]),
                c='blue',
                s=10
                )
    # band2 = 1-band[band>0.5]
    # plt.scatter(band2,
    #             [freq[i]/2/np.pi]*len(band2),
    #             c='red',
    #             s=10
    #             )
    # band2 = 1+band[band<-0.5]
    # plt.scatter(band2,
    #             [freq[i]/2/np.pi]*len(band2),
    #             c='red',
    #             s=10
    #             )

#plt.xlim(0.5,0)
# plt.xlim(0,1)
# plt.hlines([0.8, 0.9, 1.0], -1, 1, lw = 0.5)
plt.savefig('band.png')
plt.show()
