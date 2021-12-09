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
                                                                                                                                                                                      '  NP ='+'%4i'%(npts)+'  ZINF ='+
                  '%19.15f'%(zinf)+'  ZSUP ='+'%19.15f'%(zsup)+'\n'
                                                               '  THETA/AK(1) ='+'%19.15f'%(ak1)+'     FI/AK(2) ='+'%19.15f'%(ak2)+'   POLAR ='+polar+'     FEIN =   0.00\n'
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


def eval(i,j):
    create_input(npts, angle_param1[i], angle_param2[j], from_y*2*np.pi, to_y*2*np.pi, polar, lmax, r_ratio, rmax, epssph_re, epssph_im,
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




fig = plt.figure(figsize = (10,10))
plt.rcParams['font.size'] = '14'
show = 1
is_save_needed = 1
lmax = 5
rmax = 16
a = 400
s = 100
r_ratio = s/a
polar='S' # S or P
from_y = 0.502
to_y = from_y + 1e-7*from_y
npts = 2
epssph_re = 220.0
epssph_im = 0.0000
is_multipole_type_selected = '1'
is_multipole_order_selected = '1'
is_m_projection_selected = '1'
type = '1'
order = '5'
m = '0'
# angle_param2 = 0
kxpts = 100
kypts = 100
# kpts = 5
# data = np.empty((4, npts))
R = np.empty((kxpts, kypts))
# R = np.empty((2*kpts, 2*kpts))
ktype = 2


if ktype == 2:
    plt.xlabel(r'${k_x d/\pi}$')
    from_angle_param1 = (-1.0+0.0001)/2
    to_angle_param1 = (1.0-0.0001)/2
    angle_param1 = np.linspace(from_angle_param1, to_angle_param1, kxpts)
    angle_param2 = np.linspace(from_angle_param1, to_angle_param1, kypts)
    counter = 1
    for i in range(kxpts):
        for j in range(kypts):
            print(counter, 'of', kxpts*kypts)
            R[i,j] = eval(i,j)[2,0]
            counter += 1



num = 1e-15
R[R<num] = num

im = plt.imshow(R.T, extent = (-1, 1, 1, -1), cmap=cm.hot, norm=LogNorm(), aspect='auto', interpolation = 'none')
cb = plt.colorbar(im)
cb.set_label('reflectance')
#------------------------
plt.ylabel(r'${\omega d / 2\pi c }$')
sign_ax1 = ('d=%i'%(a)+'npts%i'%(npts)+'__pol_'+polar+'_epssph%f'%(epssph_re))
plt.title(sign_ax1)
plt.gca().invert_yaxis()
if show:
    plt.show()

# if is_save_needed:
#     save_result('txt', 'kx_ky', 'x_y', angle_param1*2)
#     save_result('txt', 'kx_ky', str(a)+ '_eps_'+ str(epssph_re), R.T)