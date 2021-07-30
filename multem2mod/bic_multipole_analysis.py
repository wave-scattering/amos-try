import matplotlib.pyplot as plt
import numpy as np
import os
import subprocess
import matplotlib.cm as cm
from matplotlib.colors import LogNorm
import math as m

def create_input(npts, ak1, ak2, zinf, zsup, polar, lmax, r_ratio, rmax, epssph_re, epssph_im, type, order,
                 is_multipole_type_selected, is_multipole_order_selected, is_m_projection_selected):

    #-------------------------------- input file for Sylvia Swiecicki(2017) ---------------------------------

    str_fort10 = ('           ********************************************\n'
                  '           ********INPUT FILE FOR TRANSMISSION*********\n'
                  '           ********************************************\n'
                  '   KTYPE ='+'%2i'%(ktype)+'   KSCAN = 1   KEMB  = 0    LMAX ='+'%2i'%(lmax)+'   NCOMP = 1   NUNIT = 1\n'
                  ' ALPHA =    1.000000  BETA =    1.000000   FAB =   90.000000  RMAX ='+'%11.6f'%(rmax)+'\n'
                  '  NP ='+'%4i'%(npts)+'  ZINF =  '+
                  '%11.8f'%(zinf)+'  ZSUP =  '+'%12.9f'%(zsup)+'\n'
                  '  THETA/AK(1) =  '+'%11.8f'%(ak1)+'     FI/AK(2) =  '+'%11.8f'%(ak2)+'   POLAR ='+polar+'     FEIN =   0.00\n'
                  '\n'
                  'Give information for the "NCOMP" components \n'
                  '\n'
                  '     IT  = 2\n'
                  '     MUMED =   1.00000000   0.00000000     EPSMED=   1.00000000   0.00000000\n'
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


def eval(i, npts, zinf):
    create_input(npts, angle_param1[i], angle_param2, zinf, zsup, polar, lmax, r_ratio, rmax, epssph_re, epssph_im,
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
    return d


def save_result(regime, folder_name, filename, data):
    if not (os.path.exists(folder_name)):
        os.makedirs(folder_name)

    if regime == 'jpg':
        plt.savefig(folder_name+'/'+filename+'.png')
        plt.close()

    if regime == 'txt':
        # np.savetxt('M103spectra/103,k='+sign+'.txt', spectra_w, delimiter = ',')
        np.savetxt(folder_name+'/'+filename+'.txt', data, delimiter = ',')
    print('results were saved')


is_fano_fit_data_needed = 1

lmax = 5
rmax = 5
a = 350
s = 100
r_ratio = s/a
polar='S' # S or P
# from_y = 0.5
# to_y = 0.55
# zinf = from_y*2*np.pi
# zsup = to_y*2*np.pi
npts = 1000
epssph_re = 50.0
epssph_im = 0.000000
is_multipole_type_selected = '1'
is_multipole_order_selected = '1'
is_m_projection_selected = '1'
type = '1'
order = '5'
m = '5'
angle_param2 = 0
kpts = 5
data = np.empty((kpts, 4, npts))
R = np.empty((kpts, npts))
ktype = 2

k_values = np.linspace(3.0e-1, 3.1e-1, 2)
show = 0

k = 0
r = 0

for k_value in k_values:
    k += 1
    print(k, 'is calculating from', len(k_values))
    is_save_needed = 1

    #M1X3
    # from_y = 0.45
    # to_y = 0.46

    #M1X5
    # from_y = 0.643124949
    # to_y =   0.64312495559
    from_y = 0.6431249633
    to_y =   0.643124985
    delta_w = to_y - from_y
    # from_y = 0.6431249725660459
    # to_y =   0.6431249725660465
    zinf = from_y*2*np.pi
    zsup = to_y*2*np.pi
    #TODO DRY
    if ktype == 1:
        from_angle_param1 = 0.01
        to_angle_param1 = 10.0
        angle_param1 = np.linspace(from_angle_param1, to_angle_param1, kpts)
        x = angle_param1

        for i in range(kpts):
            print(i+1, 'of', kpts)
            R[i,:] = eval(i, npts, zinf)[2,:]

    if ktype == 2:
        # from_angle_param1 = 1.5e-1/2
        # to_angle_param1 = (2.2e-1 - 0.000)/2
        delta_ap = k_value/kpts
        from_angle_param1 = (k_value - delta_ap)/2
        to_angle_param1 = (k_value + delta_ap)/2
        angle_param1 = np.linspace(from_angle_param1, to_angle_param1, kpts)
        x = angle_param1*2


    #TODO this works bad, needed to fix
        for i in range(kpts):
            print(i+1, 'of', kpts)
            delta_z = (zsup - zinf)/npts
            kx_max = zinf/a
            if (angle_param1[i]*2*np.pi/a - kx_max >= 1e-7):
                current_zinf = angle_param1[i]*2*np.pi + 0.001
                print('new y =', current_zinf/2/np.pi)
                current_npts = int((zsup - current_zinf)/delta_z)
                nan_till = npts-current_npts
                R[i, 0:nan_till] = np.nan
                d = eval(i, current_npts, current_zinf)
                R[i, nan_till:] = d[2, :]

            else:
                R[i,:] = eval(i, npts, zinf)[2,:]

    R[R<1e-10] = 1e-10

    fig = plt.figure(figsize = (10,10))
    plt.rcParams['font.size'] = '14'
    im = plt.imshow(R.T, extent = (np.amin(x), np.amax(x), to_y, from_y), cmap=cm.hot, norm=LogNorm(), aspect='auto')#, interpolation = 'nearest')
    cb = plt.colorbar(im)
    cb.set_label('reflectance')
    #------------------------
    plt.ylabel(r'${\omega d / 2\pi c }$')
    plt.xlabel(r'${k_x d/\pi}$')
    # plt.xlabel(r'${\theta}$')
    # ax1.set_xticks(np.arange(min(x)-0.001, max(x)+0.01, 0.1))
    sign_ax1 = ('d=%i'%(a)+'npts%i'%(npts)+'__pol_'+polar+'_epssph%f'%(epssph_re))
    plt.title(sign_ax1)
    plt.gca().invert_yaxis()
    if show:
        plt.show()


    if is_fano_fit_data_needed:
        th = 0.009
        dots_needed = int(npts*th)
        j = 1
        spectra_optimization_counter = 0
        while True:
            # spectra_optimization_counter += 1
            print('counter', spectra_optimization_counter)
            # if (spectra_optimization_counter > 20 or delta_w < 1e-16):
            #     is_save_needed = 0
            #     print('cant optimize')
            #     plt.close()
            #     break
            for i in range(kpts):
                print(i+1, 'of', kpts)
                d = eval(i, npts, zinf)
                R[i, :] = d[2, :]

            const_x = np.median(x)
            index_const_theta = np.abs(x - const_x).argmin()
            R_slice = R[index_const_theta, :]
            y = np.linspace(from_y, to_y, npts)
            index_Rmax = np.where(R_slice == np.amax(R_slice))[0]
            if (len(y[index_Rmax]) == 1):
                print(R_slice[index_Rmax])
            else:
                print('overflow')
        #TODO DRY and j counter and w_c is not defined yet and len
            if (len(y[index_Rmax]) > 1):
                # from_y = w_c - 1.9**j*delta_w
                # print(from_y)
                # to_y = w_c + 1.9**j*delta_w
                # print(to_y)
                # zinf = from_y*2*np.pi
                # zsup = to_y*2*np.pi
                j += 1
                # from_y += 2*spectra_optimization_counter*delta_w
                # to_y -= 2*spectra_optimization_counter*delta_w
                # from_y += 0.01*(20-spectra_optimization_counter)*delta_w
                # to_y -= 0.01*(20-spectra_optimization_counter)*delta_w
                # from_y += 4*factor**spectra_optimization_counter*delta_w
                # to_y -= 4*factor**spectra_optimization_counter*delta_w
                from_y = w_c - 1.6*j*factor*delta_w
                to_y = w_c + 1.6*j*factor*delta_w
                print(from_y)
                print(to_y)
                zinf = from_y*2*np.pi
                zsup = to_y*2*np.pi
                if (spectra_optimization_counter > 20 or delta_w < 1e-17):
                    is_save_needed = 0
                    print('cant optimize')
                    plt.close()
                    break
                continue
            j = 1
            spectra_optimization_counter += 1
            num_of_dots = len(np.where(R_slice >= 0.2)[0])
            if (num_of_dots/npts >= th):
                #TODO DRY and low and high correction
                if (index_Rmax <= int(len(y)*0.05)):
                    adj_val = 1.5*factor*delta_w
                    print('too low')
                    from_y -= adj_val
                    # print(from_y)
                    to_y -= adj_val
                    # print(to_y)
                    zinf = from_y*2*np.pi
                    zsup = to_y*2*np.pi
                    continue
                elif (index_Rmax >= int(len(y)*0.95)):
                    adj_val = 1.5*factor*delta_w
                    print('too high')
                    from_y += adj_val
                    to_y += adj_val
                    zinf = from_y*2*np.pi
                    zsup = to_y*2*np.pi
                    continue
                else:
                    print('spectra data for fano fit are ready')
                    r += 1
                    break

            print(num_of_dots, 'dots from', dots_needed, 'needed')
            w_c = float(y[index_Rmax])

            if delta_w > 1e-9:
                delta_w = 0.7**spectra_optimization_counter*abs(y[-1]-y[0])
            else:
                # delta_w = 0.9**spectra_optimization_counter*abs(y[-1]-y[0])
                delta_w = (3.1-0.05*spectra_optimization_counter)*1e-10
            #TODO binary search base on overflow
            #0.8 for 103 105
            # print('range', y[-1]-y[0])
        #TODO figure out the space of break
            if (spectra_optimization_counter > 20 or delta_w < 1e-19):
                is_save_needed = 0
                print('cant optimize')
                plt.close()
                break
            #-------------------
            print('delta_w = ', delta_w)
            factor = 1
            from_y = w_c - factor*delta_w
            print(from_y)
            to_y = w_c + factor*delta_w
            print(to_y)
            zinf = from_y*2*np.pi
            zsup = to_y*2*np.pi


            fig, (ax1, ax2) = plt.subplots(2, 1, figsize = (10,10))
            plt.rcParams['font.size'] = '14'
            im = ax1.imshow(R.T, extent = (np.amin(x), np.amax(x), to_y, from_y), cmap=cm.hot, norm=LogNorm(), aspect='auto')#, interpolation = 'nearest')
            # TODO fix it
            # cb = plt.colorbar(im)
            # cb.set_label('reflectance')
            #------------------------
            ax1.set_ylabel(r'${\omega d / 2\pi c }$')
            # ax1.set_xlabel(r'${k_x d/\pi}$')
            ax1.set_xlabel(r'${\theta}$')
            # ax1.set_xticks(np.arange(min(x)-0.001, max(x)+0.01, 0.1))
            sign_ax1 = ('d=%i'%(a)+'npts%i'%(npts)+'__pol_'+polar+'_epssph%f'%(epssph_re))
            ax1.set_title(sign_ax1)
            ax1.invert_yaxis()
            ax2.plot(y, R_slice)
            sign_ax2 = ('theta=%f'%(const_x))
            ax2.set_title(sign_ax2)
            if show:
                plt.show()
            # plt.close()

        fig, (ax1, ax2) = plt.subplots(2, 1, figsize = (10,10))
        plt.rcParams['font.size'] = '14'
        im = ax1.imshow(R.T, extent = (np.amin(x), np.amax(x), to_y, from_y), cmap=cm.hot, norm=LogNorm(), aspect='auto')#, interpolation = 'nearest')
        # TODO fix it
        # cb = plt.colorbar(im)
        # cb.set_label('reflectance')
        #------------------------
        ax1.set_ylabel(r'${\omega d / 2\pi c }$')
        # ax1.set_xlabel(r'${k_x d/\pi}$')
        ax1.set_xlabel(r'${\theta}$')
        # ax1.set_xticks(np.arange(min(x)-0.001, max(x)+0.01, 0.1))
        sign_ax1 = ('d=%i'%(a)+'npts%i'%(npts)+'__pol_'+polar+'_epssph%f'%(epssph_re))
        ax1.set_title(sign_ax1)
        ax1.invert_yaxis()
        ax2.plot(y, R_slice)
        sign_ax2 = ('theta=%f'%(const_x))
        ax2.set_title(sign_ax2)
        # plt.show()

        # sign_jpg = sign_ax1 + sign_ax2
        # sign_txt = const_x

        txt_filename = '155,k='+str(const_x)

        print('fitted:', r)
        if is_save_needed:
            spectra_w_param = np.array(list(zip(y, R_slice)))
            save_result('jpg', 'M155figures', str(const_x), data=None)
            save_result('txt', 'M155spectra', txt_filename, spectra_w_param)
