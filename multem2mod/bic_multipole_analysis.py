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


lmax = 3
rmax = 5
a = 250
s = 100
r_ratio = s/a
polar='S' # S or P
from_y = 0.1
to_y = 1.5
zinf = from_y*2*np.pi
zsup = to_y*2*np.pi
npts = 1000
epssph_re = 20.0
epssph_im = 0.0
is_multipole_type_selected = '1'
is_multipole_order_selected = '1'
is_m_projection_selected = '1'
type = '1'
order = '3'
m = '0'
angle_param2 = 0
kpts = 50
data = np.empty((kpts, 4, npts))
R = np.empty((kpts, npts))
ktype = 2
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
    from_angle_param1 = 0.001/2
    to_angle_param1 = (1.0 - 0.001)/2
    angle_param1 = np.linspace(from_angle_param1, to_angle_param1, kpts)
    # /2/np.pi*a
    x = angle_param1*a/np.pi

    for i in range(kpts):
        print(i+1, 'of', kpts)
        delta_z = (zsup - zinf)/npts
        # ------ triangle draw ------
        # kx = ak1[i]*a/np.pi
        # print('kx=', ak1[i])



        kx_max = zinf/a
        # print('kx_max=', kx_max)
        if (angle_param1[i]*2*np.pi/a - kx_max >= 1e-7):
            # print('test')
            current_zinf = angle_param1[i]*2*np.pi + 0.001
            print('new y =', current_zinf/2/np.pi)
            # current_npts = m.ceil((zsup - current_zinf)/delta_z)
            current_npts = int((zsup - current_zinf)/delta_z)
            nan_till = npts-current_npts
            R[i, 0:nan_till] = np.nan
            d = eval(i, current_npts, current_zinf)
            R[i, nan_till:] = d[2, :]
            # R[i,nan_till:] = eval(i, current_npts, current_zinf)[2,:]

        else:
            R[i,:] = eval(i, npts, zinf)[2,:]



R[R<1e-5] = 1e-5

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
plt.show()