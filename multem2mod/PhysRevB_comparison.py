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


def create_input(npts, theta, fi, zinf, zsup, polar, lmax, r_ratio, rmax, epssph_re, epssph_im, type, order,
                 is_multipole_type_selected, is_multipole_order_selected, is_m_projection_selected):

    #-------------------------------- input file for Sylvia Swiecicki(2017) ---------------------------------

    str_fort10 = ('           ********************************************\n'
                  '           ********INPUT FILE FOR TRANSMISSION*********\n'
                  '           ********************************************\n'
                  '   KTYPE = 1   KSCAN = 2   KEMB  = 0    LMAX ='+'%2i'%(lmax)+'   NCOMP = 1   NUNIT = 1\n'
                  ' ALPHA =    1.000000  BETA =    1.000000   FAB =   60.000000  RMAX ='+'%11.6f'%(rmax)+'\n'
                  '  NP ='+'%4i'%(npts)+'  ZINF =  '+
                  '%11.8f'%(zinf)+'  ZSUP =  '+'%12.9f'%(zsup)+'\n'
                  '  THETA/AK(1) =  '+'%11.8f'%(theta)+'     FI/AK(2) =  '+'%11.8f'%(fi)+'   POLAR ='+polar+'     FEIN =   0.00\n'
                  '\n'
                  'Give information for the "NCOMP" components \n'
                  '\n'
                  '     IT  = 2\n'
                  '     MUMED =   1.00000000   0.00000000     EPSMED=   2.10250000   0.00000000\n'
                  '   NPLAN = 1  NLAYER = 1\n'
                  '       S =   '+'%10.8f'%(r_ratio)+'     MUSPH =   1.00000000   0.00000000     EPSSPH=  '+'%11.8f'%(epssph_re)+'   '+'%11.8f'%(epssph_im)+'\n'
                  'xyzDL 0.0  0.0  0.0\n'
                  'xyzDR 0.0  0.0  1.0\n'
                  '     MUEMBL=   1.00000000   0.00000000    EPSEMBL=   2.10250000   0.00000000\n'
                  '     MUEMBR=   1.00000000   0.00000000    EPSEMBR=   2.10250000   0.00000000\n')


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

#         [regime]
# ;;                  ps   mz
# ;;                 / \
#     ;multipole_type =  0 0   1
# ;multipole_order = 1 1   1
# ;m_projections =  -1 1   0
#
#
# ;;                  ps
# ;;                 / \
#     ;multipole_type =  0 0
# ;multipole_order = 1 1
# ;m_projections =  -1 1


def eval(i, zinf, zsup, epssph_re, epssph_im, rmax):
    create_input(npts, np.arcsin(theta[i])/np.pi*180, fi, zinf, zsup, polar, lmax, r_ratio, rmax, epssph_re, epssph_im,
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


def plot_figure_4(lambda_incident_arr, epssph_re_arr, epssph_im_arr, rmax_arr):
    plt.figure(figsize=(11,6))
    data_arr_for_several_calc = np.empty((len(lambda_incident_arr), kpts, 4))
    x = np.linspace(from_sin_theta, to_sin_theta, kpts)
    folder_name = 'figure4_PhysRevB'
    if not (os.path.exists(folder_name)):
        os.makedirs(folder_name)

    for rmax in rmax_arr:
        print('RMAX = ', rmax,'is calculating...')
        for n in range(len(lambda_incident_arr)):
            lambda_incident = lambda_incident_arr[n]
            zinf = lambda_incident/a
            zsup = (lambda_incident+0.01)/a
            epssph_re = epssph_re_arr[n]
            epssph_im = epssph_im_arr[n]

            for i in range(kpts):
                print('sin(theta)=',theta[i])
                eval(i, zinf, zsup, epssph_re, epssph_im, rmax)

            data_arr_for_several_calc[n, :, :] = data_arr
            plt.plot(x, data_arr_for_several_calc[n, :, 2] + n, label='refl lambda=%i'%lambda_incident, lw = 2.0)
            plt.title('Rmax=%f,'%rmax)
            plt.ylim(-0.01, 3.01)
            plt.legend()
            plt.tight_layout()

        file_name = f'rmax{rmax}.png'
        plt.savefig(folder_name+'/'+file_name)
        plt.clf()
        # plt.show()

    print('Calculation has finished')


if __name__ == "__main__":

#-------------------------------- plot figure 4 from PhysRevB ----------------------------------------------------------
    # plt.figure(figsize=(11,6))
    #
    # fi = 0
    # from_sin_theta = 0.0
    # to_sin_theta = 0.99
    # n_theta = 1000
    # theta = np.linspace(from_sin_theta, to_sin_theta, n_theta)
    # kpts = len(theta)
    #
    # a = 475.0
    # s = 100.0
    # r_ratio = s/a
    # npts = 2
    # lmax= 2
    # polar='S' # S or P
    #
    # type = -1
    # order = -1
    # s_m = 0
    # m = 1
    #
    # data_arr = np.empty((kpts, 4))
    #
    # rmax_arr = np.linspace(1.0, 1.0, 1)
    # lambda_incident_arr = [650, 750, 900]
    # epssph_re_arr = [-12.9530000, -20.1480000, -32.7190000]
    # epssph_im_arr = [1.12090000, 1.24700000, 1.99550000]
    #
    # plot_figure_4(lambda_incident_arr, epssph_re_arr, epssph_im_arr, rmax_arr)

#------------------------------------- plot figure 7 from PhysRevB ----------------------------------------------------------------------

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize = (10, 8))

    lmax= 2

    rmax_arr = np.linspace(5.0, 19.0, 1)

    for rmax in rmax_arr:
        print('RMAX = ', rmax,'is calculating...')
        fig.suptitle(f'RMAX = {rmax}')

# -------- first plot --------------------------------------------------------------------------------------------------

        fi = 0
        from_sin_theta = 0.0
        to_sin_theta = 0.4
        n_theta = 100
        theta = np.linspace(from_sin_theta, to_sin_theta, n_theta)
        kpts = len(theta)
        data_arr = np.empty((kpts, 4))
        x = np.linspace(from_sin_theta, to_sin_theta, kpts)

        a = 475.0
        s = 100.0
        r_ratio = s/a
        lambda_incident = 750
        zinf = lambda_incident/a
        zsup = (lambda_incident+0.01)/a
        epssph_re = -20.1480000
        epssph_im = 1.24700000
        npts = 2
        polar='S' # S or P

        lmax = 2


        is_multipole_type_selected = '0'
        is_multipole_order_selected = '0'
        is_m_projection_selected = '0'
        type = '-1'
        order = '1'
        m = '0'

        for i in range(kpts):
            eval(i, zinf, zsup, epssph_re, epssph_im, rmax)

        ax1.plot(x, data_arr[:, 2], 'b--', label='lmax = 1', lw = 2.0)

        lmax = 2

        is_multipole_type_selected = '1'
        is_multipole_order_selected = '1'
        is_m_projection_selected = '1'
        type = '0 0'
        order = '1 1'
        m = '-1 1'

        # TODO check if is_m_proj_selected = 0 and rest are 1

        for i in range(kpts):
            eval(i, zinf, zsup, epssph_re, epssph_im, rmax)

        ax1.plot(x, data_arr[:, 2], 'r', label='Ps', lw = 2.0)

#------------ plot config ----------------------------------------------------------------------------------------------
        ax1.set_title(f'{lambda_incident} nm')
        ax1.set_ylabel('R')
        ax1.set_xlabel(r'sin$\theta$')
        ax1.set_ylim(-0.01, 1.01)
        ax1.legend()


# -------- second plot --------------------------------------------------------------------------------------------------

        # from_sin_theta = 0.4
        # to_sin_theta = 0.8
        # n_theta = 1000
        # theta = np.linspace(from_sin_theta, to_sin_theta, n_theta)
        # kpts = len(theta)
        # x = np.linspace(from_sin_theta, to_sin_theta, kpts)
        #
        # lambda_incident = 900
        # zinf = lambda_incident/a
        # zsup = (lambda_incident+0.01)/a
        # epssph_re = -32.7190000
        # epssph_im = 1.99550000
        # a = 475.0
        # s = 100.0
        # r_ratio = s/a
        # npts = 2
        # polar='S' # S or P
        #
        # type = 0
        # order = 1
        # s_m = 1
        # m = 0
        #
        # for i in range(kpts):
        #     eval(i, zinf, zsup, epssph_re, epssph_im, rmax)
        #
        # ax2.plot(x, data_arr[:, 2], 'r--', label='Ps', lw = 2.0)
        #
        # type = -1
        # order = -1
        # s_m = 0
        # m = 0
        #
        # for i in range(kpts):
        #     eval(i, zinf, zsup, epssph_re, epssph_im, rmax)
        #
        # ax2.plot(x, data_arr[:, 2], 'r', label='all multipoles', lw = 2.0)

#------------ plot config ----------------------------------------------------------------------------------------------
        # ax2.set_title(f'{lambda_incident} nm')
        # ax2.set_xlabel(r'sin$\theta$')
        # ax2.set_ylim(-0.01, 1.01)
        # ax2.legend()
        # plt.tight_layout()

# ----------- plot saving -----------------------------------------------------------------------------------------------
#         folder_name = 'figure7_PhysRevB'
#         if not (os.path.exists(folder_name)):
#             os.makedirs(folder_name)
#         file_name = f'rmax{rmax}.png'
#         plt.savefig(folder_name+'/'+file_name)


        plt.show()
        ax1.clear()
        ax2.clear()






