# ------------------------------------------------------------------------------------------------------------------------
#  this script draws figures from PHYSICAL REVIEW B 95, 195406 (2017)
#     Surface-lattice resonances in two-dimensional arrays of spheres: Multipolar interactions and a mode analysis
#     Sylvia D. Swiecicki and J. E. Sipe
#  using MULTEM and multipole decomposition modification
# ------------------------------------------------------------------------------------------------------------------------

import matplotlib.pyplot as plt
import numpy as np
import os
import subprocess


def create_input(npts, ak1, ak2, zinf, zsup, polar, lmax, r_ratio, rmax, epssph_re, epssph_im, type, order,
                 is_multipole_type_selected, is_multipole_order_selected, is_m_projection_selected):


    str_fort10 = ('           ********************************************\n'
                  '           ********INPUT FILE FOR TRANSMISSION*********\n'
                  '           ********************************************\n'
                  '   KTYPE ='+'%2i'%(ktype)+'   KSCAN = 2   KEMB  = 0    LMAX ='+'%2i'%(lmax)+'   NCOMP = 1   NUNIT = 1\n'
                  ' ALPHA =    1.000000  BETA =    1.000000   FAB =   60.000000  RMAX ='+'%11.6f'%(rmax)+'\n'
                  '  NP ='+'%4i'%(npts)+'  ZINF ='+
                  '%19.15f'%(zinf)+'  ZSUP ='+'%19.15f'%(zsup)+'\n'
                  '  THETA/AK(1) ='+'%19.15f'%(ak1)+'     FI/AK(2) ='+'%19.15f'%(ak2)+'   POLAR ='+polar+'     FEIN =   0.00\n'
                  '\n'
                  'Give information for the "NCOMP" components \n'
                  '\n'
                  '     IT  = 2\n'
                  '     MUMED =   1.00000000   0.00000000     EPSMED=   2.10250000   0.00000000\n'
                  '   NPLAN = 1  NLAYER = 1\n'
                  '       S =   '+'%10.8f'%(r_ratio)+'     MUSPH =   1.00000000   0.00000000     EPSSPH=  '+'%11.6f'%(epssph_re)+'   '+'%11.8f'%(epssph_im)+'\n'
                  'xyzDL 0.0  0.0  0.0\n'
                  'xyzDR 0.0  0.0  1.0\n'
                  '     MUEMBL=   1.00000000   0.00000000    EPSEMBL=   1.00000000   0.00000000\n'
                  '     MUEMBR=   1.00000000   0.00000000    EPSEMBR=   1.00000000   0.00000000\n')

    with open('fort.10','w') as f:
        print(str_fort10, file=f)

    str_ini = ('[selectors]\n'
    'is_multipole_type_selected = '+'%s'%(is_multipole_type_selected)+'\n'
    'is_multipole_order_selected = '+'%s'%(is_multipole_order_selected)+'\n'
    'is_m_projection_selected = '+'%s'%(is_m_projection_selected)+'\n'
    '\n'
    '[regime]\n'
    'multipole_type = '+'%s'%(type)+'\n'
    'multipole_order = '+'%s'%(order)+'\n'
    'm_projection = '+'%s'%(m)+'\n')

    with open('multipole_regime_parameters.ini','w') as f:
        print(str_ini, file=f)


def eval(i, zinf, zsup, epssph_re, epssph_im, rmax):
    create_input(npts, theta[i], fi, zinf, zsup, polar, lmax, r_ratio, rmax, epssph_re, epssph_im,
                 type, order, is_multipole_type_selected, is_multipole_order_selected, is_m_projection_selected)
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
    print(data_arr[i,:])


if __name__ == "__main__":

#plot figure 7 from PhysRevB

    fig, ax1= plt.subplots(1, 1, figsize = (8, 10))
    ktype = 1
    lmax = 3
    rmax_arr = np.linspace(5.0, 20.0, 1)

    for rmax in rmax_arr:
        print('RMAX = ', rmax,'is calculating...')
        fig.suptitle(f'RMAX = {rmax}')

#  figure 7 a
        fi = 0
        from_sin_theta = 0.0
        to_sin_theta = 0.4
        n_theta = 100
        sin_theta = np.linspace(from_sin_theta, to_sin_theta, n_theta)
        x = sin_theta
        theta = np.arcsin(sin_theta)*180/np.pi
        kpts = n_theta
        data_arr = np.empty((kpts, 4))
        a = 475.0
        s = 100.0
        r_ratio = s/a
        lambda_incident = 750
        epssph_re = -20.1480000
        epssph_im = 1.24700000
        zinf = lambda_incident/a
        zsup = (lambda_incident+0.01)/a
        npts = 2
        polar='S' # S or P

# ps multipole
        is_multipole_type_selected = '1'
        is_multipole_order_selected = '1'
        is_m_projection_selected = '1'
        type = '0 0'
        order = '1 1'
        m = '-1 1'

        for i in range(kpts):
            eval(i, zinf, zsup, epssph_re, epssph_im, rmax)

        ax1.plot(x, data_arr[:, 2], 'r--', label='ps', lw = 2.0)

# ps + mz + qks
        is_multipole_type_selected = '1'
        is_multipole_order_selected = '1'
        is_m_projection_selected = '1'
        type =  '0 0 1 0 0'
        order = '1 1 1 2 2'
        m =     '-1 1 0 -2 2'

        for i in range(kpts):
            eval(i, zinf, zsup, epssph_re, epssph_im, rmax)

        ax1.plot(x, data_arr[:, 2], 'r', label='ps+mz+qks', lw = 2.0)

# plot config
        ax1.set_title(f'{lambda_incident} nm')
        ax1.set_ylabel('R')
        ax1.set_xlabel(r'sin$\theta$')
        ax1.set_ylim(-0.01, 1.01)
        ax1.legend()

        plt.show()








