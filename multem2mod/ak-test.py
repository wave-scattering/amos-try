import matplotlib.pyplot as plt
import numpy as np
import os
import subprocess
import matplotlib.cm as cm
from matplotlib.colors import LogNorm


def create_input(npts, ak1, ak2, zinf, zsup, polar, lmax, r_ratio, rmax, epssph_re, epssph_im, type, order,
                 is_multipole_type_selected, is_multipole_order_selected, is_m_projection_selected):

    #-------------------------------- input file for Sylvia Swiecicki(2017) ---------------------------------

    str_fort10 = ('           ********************************************\n'
                  '           ********INPUT FILE FOR TRANSMISSION*********\n'
                  '           ********************************************\n'
                  '   KTYPE = 2   KSCAN = 1   KEMB  = 0    LMAX ='+'%2i'%(lmax)+'   NCOMP = 1   NUNIT = 1\n'
                  ' ALPHA =    1.000000  BETA =    1.000000   FAB =   60.000000  RMAX ='+'%11.6f'%(rmax)+'\n'
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


def eval(i):
    create_input(npts, ak1[i], ak2, zinf, zsup, polar, lmax, r_ratio, rmax, epssph_re, epssph_im,
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
    data[i,:, :] = d



plt.figure(figsize=(11,6))

lmax = 4
rmax = 5
# a = 300
# s = 180
# a = 600
# s = 300
r_ratio = 0.47050000
npts = 500
polar='S' # S or P

zinf = 3.73
zsup = 3.741
# zinf = 1
# zsup = 4
epssph_re = 15.0
epssph_im = 0.0

ak2 = 0
# from_ak1 = 0.0
# to_ak1 = 0.5
# n_ak1 = 100
# ak1 = np.linspace(from_ak1, to_ak1, n_ak1)
ak1 = [0.05/2/np.pi]
kpts = len(ak1)
data= np.empty((kpts, 4, npts))
T = np.empty((kpts, npts))

# is_multipole_type_selected = '1'
# is_multipole_order_selected = '1'
# is_m_projection_selected = '1'
# type = '0 0 0 1 1 1'
# order = '1 1 1 1 1 1'
# m = '-1 1 0 -1 0 1'

is_multipole_type_selected = '0'
is_multipole_order_selected = '0'
is_m_projection_selected = '0'
type = '0 0 0'
order = '1 1 '
m = '1'

for i in range(kpts):
    print(i+1, 'of', kpts)
    eval(i)
    T[i, :] = data[i, 1, :]

x = np.linspace(zinf, zsup, npts)
plt.plot(x, T[0, :])
plt.show()