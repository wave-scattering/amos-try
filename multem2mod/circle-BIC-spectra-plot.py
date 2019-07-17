#!/usr/bin/python3
# -*- coding: UTF-8 -*-
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
import os
import sys
import subprocess
import multiprocessing
import ctypes

def create_input(npts,ak1, ak2, zinf, zsup, polar, lmax, r_ratio):
    str = ('           ********************************************\n'
           '           ********INPUT FILE FOR TRANSMISSION*********\n'
           '           ********************************************\n'
           '   KTYPE = 2   KSCAN = 1   KEMB  = 0    LMAX ='+'%2i'%(lmax)+'   NCOMP = 1   NUNIT = 1\n'
           ' ALPHA =    1.000000  BETA =    1.000000   FAB =   60.000000  RMAX =   16.000000\n'
           '  NP ='+'%4i'%(npts)+'  ZINF =  '+
           '%11.8f'%(zinf)
           +'  ZSUP =  '+'%12.9f'%(zsup)+'\n'
           '  THETA/AK(1) =  '+'%11.8f'%(ak1)+'     FI/AK(2) =  '+'%11.8f'%(ak2)+'   POLAR ='+polar+'     FEIN =   0.00\n'
           '\n'
           'Give information for the "NCOMP" components \n'
           '\n'
           '     IT  = 2\n'
           '     MUMED =   1.00000000   0.00000000     EPSMED=   1.00000000   0.00000000\n'
           '   NPLAN = 1  NLAYER = 1\n'
           '       S =   '+'%10.8f'%(r_ratio)+'     MUSPH =   1.00000000   0.00000000     EPSSPH=  15.00000000   0.00000000\n'
           'xyzDL 0.0  0.0  0.0\n'
           'xyzDR 0.0  0.0  1.0\n'
           '     MUEMBL=   1.00000000   0.00000000    EPSEMBL=   1.00000000   0.00000000\n'
           '     MUEMBR=   1.00000000   0.00000000    EPSEMBR=   1.00000000   0.00000000\n')
    with open('fort.10','w') as f:
        print(str, file=f)


def get_shared_array(length, width, depth):
    shared_array_base = multiprocessing.Array(ctypes.c_double, length*width*depth)
    shared_array = np.ctypeslib.as_array(shared_array_base.get_obj())
    shared_array = shared_array.reshape(length, width, depth)
    return shared_array


def eval(i):
    create_input(npts, ak1[i], ak2[i],zinf, zsup, polar,lmax, r_ratio)
    if os.path.isfile('multem2'):
        subprocess.run(['./multem2'], stdout=subprocess.DEVNULL)
    #FREQUENCY   TRANSMITTANCE  Reflectance   Absorbance
    d = np.loadtxt('fort.8').T
    C_all[i,:,:] = d



if __name__ == "__main__":
    plt.figure(figsize=(11,6))
    npts = 1404
    lmax=10
    r_ratio = 0.47050000
    # from_k = 3.7128
    # to_k = 3.7129
    plot_factor = 1#/3.666*1.75
    #from_k = 1.6/plot_factor
    #to_k = 2.25/plot_factor

    # As in article
    from_k = 3.73
    to_k = 3.741

    # from_k = 3.3
    # to_k = 4.6
    # kk = 0.01
    #
    zinf = from_k
    zsup = to_k
    ak_step = 0.501
    # abs_k = np.arange(0.42,0.421, ak_step)/2/np.pi
    # kk = 0.41
    kk = 0.01
    abs_k = np.arange(kk,kk+ak_step/2, ak_step)/2/np.pi
    # abs_k = np.arange(0.3,0.301, ak_step)/2/np.pi
    phi_all = np.arange(0,65,100)
    polar='S' # S or P
    # ak1 = np.arange(0.1,0.15, ak_step)/2/np.pi
    # ak2 = np.arange(0.1,0.15, ak_step)/2/np.pi
    for phi_grad in phi_all:
        phi = np.pi*phi_grad/180
        zinf = from_k
        zsup = to_k
        sign = ('npts%i'%(npts)+'_lmax%i'%(lmax)+'_abs_k_%g_%g'%(abs_k[0],abs_k[-1])+
                '_phi_%g'%(phi_grad)+'_z_%g_%g'%(zinf, zsup)+'_pol'+polar+'_radius%g'%(r_ratio))
        ak1 = np.cos(phi)*abs_k
        ak2 = np.sin(phi)*abs_k
        kpts = len(ak1)
        print(sign)
        if os.path.isfile(sign+'.npz'):# and False:
            data = np.load(sign+'.npz')
            C_all = data['C_all']
        else:
            C_all = get_shared_array(kpts,4,npts)
            for i in range(kpts):
                print(i+1, 'of', kpts)
                eval(i)
            np.savez_compressed(sign, C_all=C_all)
            # sys.exit(0)
        for i in range(kpts):
            d = C_all[i,:,:]
            z = np.linspace(from_k, to_k, npts)
            plt.title('Reflectance New\n'+sign)
            plt.plot(z#*plot_factor
                     , d[1], label='trans |k|=%g'%(abs_k[i]*2*np.pi), lw=0.4)
            plt.plot(z#*plot_factor
                     , d[2], label='refl |k|=%g'%(abs_k[i]*2*np.pi), lw=0.4)
            plt.plot(z#*plot_factor
                     , (d[1]+d[2]-1)*1e5, label='(trans+refl-1)*1e5 '%(abs_k[i]*2*np.pi), lw=0.4)
            # plt.plot(z#*plot_factor
            #          , d, label='source angle=%g'%(abs_k[i]), lw=0.4)
            # plt.xlim(from_k#*plot_factor
            #          ,to_k#*plot_factor
            # )
            plt.ylim(-0.01,1.01)
            plt.legend()
        plt.tight_layout()
        #plt.show()
        # plt.show()
        # plt.clf(); plt.close()
        # plt.figure(figsize=(16,9))
        # plt.imshow(np.flipud(C_all.T), aspect='auto',
        #            # interpolation='quadric',
        #            interpolation='none',
        #            norm=LogNorm(vmin=np.min(C_all)*1e3, vmax=np.max(C_all)*2e0),
        #            # norm=LogNorm(vmin=np.min(test)*7e3, vmax=np.max(test)/3
        #            # ),
        #            cmap=plt.cm.plasma,
        #            extent =(ak1[0],ak1[-1], zinf/2/np.pi, zsup/2/np.pi)
        #            )
        # # plt.legend()
        # plt.xlabel(r'$\theta$, deg')
        # plt.ylabel(r'$\omega d / 2 \pi c$')
        # plt.title('Multem')
        plt.savefig(sign+'.png', dpi=150)
        plt.clf(); plt.close()

        # plt.show()
