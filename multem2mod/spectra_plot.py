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

def create_input(npts,ak1, zinf, zsup, polar, lmax):
    str = ('           ********************************************\n'
           '           ********INPUT FILE FOR TRANSMISSION*********\n'
           '           ********************************************\n'
           '   KTYPE = 1   KSCAN = 1   KEMB  = 0    LMAX ='+'%2i'%(lmax)+'   NCOMP = 1   NUNIT = 1\n'
           ' ALPHA =    1.000000  BETA =    1.000000   FAB =   90.000000  RMAX =   16.000000\n'
           '  NP ='+'%4i'%(npts)+'  ZINF =  '+
           '%11.8f'%(zinf)
           +'  ZSUP =  '+'%12.9f'%(zsup)+'\n'
           '  THETA/AK(1) =  '+'%11.8f'%(ak1)+'     FI/AK(2) =   0.00000000   POLAR ='+polar+'     FEIN =   0.00\n'
           '\n'
           'Give information for the "NCOMP" components \n'
           '\n'
           '     IT  = 2\n'
           '     MUMED =   1.00000000   0.00000000     EPSMED=   1.00000000   0.00000000\n'
           '   NPLAN = 1  NLAYER = 1\n'
           '       S =   0.16666667     MUSPH =   1.00000000   0.00000000     EPSSPH=  12.00000000   0.00000000\n'
           'xyzDL 0.0  0.0  0.0\n'
           'xyzDR 0.0  0.0  1.0\n'
           '     MUEMBL=   1.00000000   0.00000000    EPSEMBL=   1.00000000   0.00000000\n'
           '     MUEMBR=   1.00000000   0.00000000    EPSEMBR=   1.00000000   0.00000000\n')
    with open('fort.10','w') as f:
        print(str, file=f)


def get_shared_array(length, width):
    shared_array_base = multiprocessing.Array(ctypes.c_double, length*width)
    shared_array = np.ctypeslib.as_array(shared_array_base.get_obj())
    shared_array = shared_array.reshape(length, width)
    return shared_array


def eval(i):
    create_input(npts, ak1[i], zinf, zsup, polar,lmax)
    if os.path.isfile('multem2'):
        subprocess.run(['./multem2'], stdout=subprocess.DEVNULL)
    #FREQUENCY   TRANSMITTANCE  Reflectance   Absorbance
    d = np.loadtxt('fort.8').T
    C_all[i,:] = d[2]


if __name__ == "__main__":

    # npts = 1500
    # ak1 = 0.15
    #ak1 = 0.15
    # ak1 = np.arange(0.,0.17, 0.01)

    # npts = 500
    # from_k = 0.6
    # to_k = 0.75
    # ak_step = 0.5
    # ak1 = np.arange(25,50, ak_step)
    # polar='S'

    npts = 501
    lmax=7
    from_k = 0.7
    to_k = 0.8
    # to_k = 0.64
    ak_step = 1
    # BIC at 19.471220634
    # ak1 = np.arange(35,37, ak_step)
    ak1 = np.arange(15,25, ak_step)
    polar='P' # S or P

    # npts = 101
    # from_k = 0.6
    # to_k = 0.75
    # ak_step = 0.5
    # ak1 = np.arange(25,50, ak_step)
    # polar='P'

    zinf = from_k*2*np.pi
    zsup = to_k*2*np.pi
    kpts = len(ak1)
    sign = 'npts%i'%(npts)+'_lmax%i'%(lmax)+'_ak1_%g_%g'%(ak1[0],ak1[-1])+'_z_%g_%g'%(zinf, zsup)+'_pol'+polar
    print(sign)
    if os.path.isfile(sign+'.npz'):
        data = np.load(sign+'.npz')
        C_all = data['C_all']
    else:
        C_all = get_shared_array(kpts,npts)
        for i in range(kpts):
            print(i, 'of', kpts)
            eval(i)
        np.savez_compressed(sign, C_all=C_all)
        sys.exit(0)
    for i in range(kpts):
        d = C_all[i,:]
        z = np.linspace(from_k, to_k, npts)
        plt.plot(z, d, label='Reflectance ak1=%4f'%(ak1[i]))
        plt.xlim(from_k, to_k)
    plt.show()
    plt.clf(); plt.close()
    plt.figure(figsize=(16,9))
    plt.imshow(np.flipud(C_all.T), aspect='auto',
               # interpolation='quadric',
               interpolation='none',
               norm=LogNorm(vmin=np.min(C_all)*1e3, vmax=np.max(C_all)*2e0),
               # norm=LogNorm(vmin=np.min(test)*7e3, vmax=np.max(test)/3
               # ),
               cmap=plt.cm.plasma,
               extent =(ak1[0],ak1[-1], zinf/2/np.pi, zsup/2/np.pi)
               )
    # plt.legend()
    plt.xlabel(r'$\theta$, deg')
    plt.ylabel(r'$\omega d / 2 \pi c$')
    plt.title('Multem')
    plt.savefig(sign+'.png', dpi=600)
    plt.show()
