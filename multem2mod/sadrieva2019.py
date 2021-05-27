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


def eval(i):
    create_input(npts, ak1[i]/2/np.pi*a, ak2, zinf, zsup, polar, lmax, r_ratio, rmax, epssph_re, epssph_im,
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
    # R = d[2, :]
    # lambda_arr = d[0, :]
    data[i,:,:] = d


plt.rcParams['font.size'] = '14'
fig, ax1 = plt.subplots(1, 1, figsize = (6,20))

lmax = 2
rmax = 5
a = 600
s = 100
r_ratio = s/a
# r_ratio = 0.4705
# s = 100
# a = s/r_ratio
#gamma tm3
polar='P' # S or P
from_y = 0.91
to_y = 0.92
#----------------------
#gamma te1
# polar='S' # S or P
# from_y = 0.814
# to_y = 0.816

# polar = 'P'
# from_y = 0.5
# to_y = 1.0

# polar='P' # S or P
# from_y = 0.5
# to_y = 1.0
# zinf = from_y*2*np.pi
# zsup = to_y
zinf = from_y*2*np.pi
zsup = to_y*2*np.pi
# zinf = 3.73
# zsup = 3.75
npts = 1000
# x = np.linspace(zinf, zsup, npts)
epssph_re = 12.0
epssph_im = 0.0
# epssph_re = 15.0
# epssph_im = 0.0

ak2 = 0
from_ak1 = 0.0*np.pi/a
to_ak1 = 0.1*np.pi/a
n_ak1 = 300
ak1 = np.linspace(from_ak1, to_ak1, n_ak1)
# ak1 = [0.01/a, 0.05/a]
kpts = len(ak1)
data= np.empty((kpts, 4, npts))
R = np.empty((kpts, npts))

#------------- all --------------------------
# is_multipole_type_selected = '0'
# is_multipole_order_selected = '0'
# is_m_projection_selected = '0'
# type = '1'
# order = '1'
# m = '-1'

is_multipole_type_selected = '1'
is_multipole_order_selected = '1'
is_m_projection_selected = '1'
type = '0'
order = '1'
m = '0'

if (is_multipole_type_selected == '0'):
    multipole = 'all multipoles'
else:
    if (type == '0' and m == '0'):
        multipole = 'pz'
    if (type == '1' and m == '0'):
        multipole = 'mz'

for i in range(kpts):
    print(i+1, 'of', kpts)
    eval(i)
    R[i, :] = data[i,2,:]
# y = 1/data[0,0,:]
# y = a/lambda_arr
# y = a/lambda_arr
R[R<1e-5] = 1e-5
kx = ak1*a/np.pi

im = ax1.imshow(R.T, extent = (np.amin(kx), np.amax(kx), to_y, from_y), cmap=cm.hot, norm=LogNorm(), aspect='auto', interpolation = 'nearest')
cb = plt.colorbar(im)
cb.set_label('reflectance')
ax1.set_ylabel(r'${\omega d / 2\pi c }$')
ax1.set_xlabel(r'${k_x d/\pi}$')
# sign = ('npts%i'%(npts)+'_lmax%i'%(lmax)+'_pol_'+polar+'_multipole_'+multipole)
# plt.title(sign)
plt.gca().invert_yaxis()
# plt.gca().invert_xaxis()
plt.show()
