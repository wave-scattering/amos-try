#!/usr/bin/python3
import matplotlib.pyplot as plt
import numpy as np
import os
import subprocess
import multiprocessing
import ctypes


def create_input(npts, theta, fi, zinf, zsup, polar, lmax, r_ratio):

    #-------------------------------- input file for Sylvia Swiecicki(2017) ---------------------------------

    str = ('           ********************************************\n'
           '           ********INPUT FILE FOR TRANSMISSION*********\n'
           '           ********************************************\n'
           '   KTYPE = 1   KSCAN = 2   KEMB  = 0    LMAX ='+'%2i'%(lmax)+'   NCOMP = 1   NUNIT = 1\n'
           ' ALPHA =    1.000000  BETA =    1.000000   FAB =   60.000000  RMAX =   5.000000\n'
           '  NP ='+'%4i'%(npts)+'  ZINF =  '+
           '%11.8f'%(zinf)
           +'  ZSUP =  '+'%12.9f'%(zsup)+'\n'
            '  THETA/AK(1) =  '+'%11.8f'%(theta)+'     FI/AK(2) =  '+'%11.8f'%(fi)+'   POLAR ='+polar+'     FEIN =   0.00\n'
            '\n'
            'Give information for the "NCOMP" components \n'
            '\n'
            '     IT  = 2\n'
            '     MUMED =   1.00000000   0.00000000     EPSMED=   2.10250000   0.00000000\n'
            '   NPLAN = 1  NLAYER = 1\n')
                                                                                #AU at 0.75 um eps = -20.1480000  1.24700000
                                                                                #AU at 0.9 um eps = -32.7190000    1.99550000
                                                                                #AU at 0.65 um eps = -12.9530000   1.12090000
    if lambda_incident == 650:
        str += '       S =   '+'%10.8f'%(r_ratio)+'     MUSPH =   1.00000000   0.00000000     EPSSPH=  -12.9530000   1.12090000\n'
    if lambda_incident == 750:
        str += '       S =   '+'%10.8f'%(r_ratio)+'     MUSPH =   1.00000000   0.00000000     EPSSPH=  -20.1480000   1.24700000\n'
    if lambda_incident == 900:
        str += '       S =   '+'%10.8f'%(r_ratio)+'     MUSPH =   1.00000000   0.00000000     EPSSPH=  -32.7190000   1.99550000\n'
    str += ('xyzDL 0.0  0.0  0.0\n'
            'xyzDR 0.0  0.0  1.0\n'
            '     MUEMBL=   1.00000000   0.00000000    EPSEMBL=   2.10250000   0.00000000\n'
            '     MUEMBR=   1.00000000   0.00000000    EPSEMBR=   2.10250000   0.00000000\n')

    #-------------------------------- input file for Zarina Sadrieva(2019) ---------------------------------

    # str = ('           ********************************************\n'
    #        '           ********INPUT FILE FOR TRANSMISSION*********\n'
    #        '           ********************************************\n'
    #        '   KTYPE = 1   KSCAN = 2   KEMB  = 0    LMAX ='+'%2i'%(lmax)+'   NCOMP = 1   NUNIT = 1\n'
    #        ' ALPHA =    1.000000  BETA =    1.000000   FAB =   90.000000  RMAX =   20.000000\n'
    #        '  NP ='+'%4i'%(npts)+'  ZINF =  '+
    #        '%11.8f'%(zinf)
    #        +'  ZSUP =  '+'%12.9f'%(zsup)+'\n'
    #         '  THETA/AK(1) =  '+'%11.8f'%(theta)+'     FI/AK(2) =  '+'%11.8f'%(fi)+'   POLAR ='+polar+'     FEIN =   0.00\n'
    #         '\n'
    #         'Give information for the "NCOMP" components \n'
    #         '\n'
    #         '     IT  = 2\n'
    #         '     MUMED =   1.00000000   0.00000000     EPSMED=   1.00000000   0.00000000\n'
    #         '   NPLAN = 1  NLAYER = 1\n'
    #         '       S =   '+'%10.8f'%(r_ratio)+'     MUSPH =   1.00000000   0.00000000     EPSSPH=  12.0000000   0.00000000\n'
    #         'xyzDL 0.0  0.0  0.0\n'
    #         'xyzDR 0.0  0.0  1.0\n'
    #         '     MUEMBL=   1.00000000   0.00000000    EPSEMBL=   1.00000000   0.00000000\n'
    #         '     MUEMBR=   1.00000000   0.00000000    EPSEMBR=   1.00000000   0.00000000\n')


    with open('fort.10','w') as f:
        print(str, file=f)


def eval(i):
    create_input(npts, np.arcsin(theta[i])/np.pi*180, fi, zinf, zsup, polar, lmax, r_ratio)
    if os.path.isfile('multem2'):
        my_env = os.environ.copy()
        my_env["OMP_NUM_THREADS"] = "1"
        # print("Running multem...")
        subprocess.run(['./multem2'],
                       stdout=subprocess.DEVNULL,
                       env=my_env)
        # print("Done.")
    #FREQUENCY   TRANSMITTANCE  Reflectance   Absorbance
    d = np.loadtxt('fort.8').T
    data_arr[i,:] = d[:,1]
    # np.savetxt(data_file, d[:,1])

if __name__ == "__main__":
    plt.figure(figsize=(11,6))


    #-------------------------------- input file for Zarina Sadrieva(2019) ---------------------------------
    fi = 0
    # from_sin_theta = 0.4
    # to_sin_theta = 0.8
    # to_sin_theta = 0.35
    n_theta = 200

    # lambda_incident = 650
    lambda_incident = 750
    from_sin_theta = 0.0
    to_sin_theta = 0.4
    # lambda_incident = 900
    # from_sin_theta = 0.4
    # to_sin_theta = 0.8
    a = 475.0
    s = 100.0
    r_ratio = s/a
    zinf = lambda_incident/a
    zsup = (lambda_incident+0.01)/a
    npts = 2
    lmax= 2
    polar='S' # S or P

    theta = np.linspace(from_sin_theta, to_sin_theta, n_theta)
    kpts = len(theta)

    data_arr = np.empty((kpts, 4))

    for i in range(kpts):
        print(i+1, 'of', kpts)
        eval(i)

    data_arr = data_arr.transpose()
    # from_sin_theta = np.sin(np.pi * from_theta_deg/180)
    # to_sin_theta = np.sin(np.pi * to_theta_deg/180)
    x = np.linspace(from_sin_theta, to_sin_theta, kpts)
    # plt.plot(x, data_arr[1], label='trans lmax=%i'%lmax, lw=1.0)
    plt.plot(x, data_arr[2], label='refl', lw=2.0)
    # plt.plot(x, (data_arr[1]+data_arr[2]-1)*1e1, label='(trans+refl-1)*1e1', lw=0.4)
    plt.ylim(-0.01,1.01)
    plt.legend()
    plt.tight_layout()
    plt.show()
