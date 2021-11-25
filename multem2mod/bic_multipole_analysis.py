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


def eval(ap1, ap2, from_y, to_y, npts):
    create_input(npts, ap1, ap2, from_y*2*np.pi, to_y*2*np.pi, polar, lmax, r_ratio, rmax, epssph_re, epssph_im,
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

#TODO use kpts_ap2 for r = f(kx,ky)
def calc_R_map(ap1_start, ap1_end, kpts_ap1, ap2_start, ap2_end, kpts_ap2, from_y, to_y, npts, R_min, regime):
    #regime: 'theta_scan' or 'phi_scan'
    #TODO autocheck for ap1 and ap2 size in case of regime
    #TODO dry for if regime
        ap1 = np.linspace(ap1_start, ap1_end, kpts_ap1)
        ap2 = np.linspace(ap2_start, ap2_end, kpts_ap2)

        if regime == 'theta_scan':
            x = ap1
            y = np.linspace(from_y, to_y, npts)
            R = np.empty((kpts_ap1, npts))
        if regime == 'phi_scan':
            x = ap2
            y = np.linspace(from_y, to_y, npts)
            R = np.empty((kpts_ap2, npts))

        if regime == 'theta_scan':
            for i in range(kpts_ap1):
                print(i+1, 'of', kpts_ap1)
                R[i,:] = eval(ap1[i], ap2, from_y, to_y, npts)[2,:]
                R[R<R_min] = R_min

        if regime == 'phi_scan':
            for i in range(kpts_ap2):
                print(i+1, 'of', kpts_ap2)
                R[i,:] = eval(ap1, ap2[i], from_y, to_y, npts)[2,:]
                R[R<R_min] = R_min


        if regime == 'freq_scan':
            x = np.linspace(from_y, to_y, npts)
            R = eval(ap1, ap2, from_y, to_y, npts)[2,:]
            R[R<R_min] = R_min
            y = R


        return x, y, R


def get_half_width_half_maxima_and_x0(x, y):
    #TODO corner cases
    idx = [np.argmax(y), np.argmin(y)]
    ymax = y[idx[0]]
    ymin = y[idx[1]]
    amp = ymax - ymin
    idx_hwhm = np.abs(y - 0.5*amp).argmin()
    x0 = x[idx[0]]
    hwhm = np.abs(x0 - x[idx_hwhm])

    return hwhm, x0

def q_factor_estimate(hwhm, x0, q_factor_limit):
    is_enough_spectra = 0
    q = x0/(2*hwhm)
    if (q >= q_factor_limit):
        is_enough_spectra = 1

    return is_enough_spectra


#TODO DRY
def find_spectrum(k_value, x, th):
    from_x = x[0]
    to_x = x[-1]
    while True:
        y = eval(k_value, 0, from_x, to_x, npts)[2,:]
        x = np.linspace(from_x, to_x, npts)
        x_range = to_x - from_x
        show_spectrum((10,10), x, y, 'theta=20')
        hwhm, x0 = get_half_width_half_maxima_and_x0(x, y)
        hwhm_factor = 2
        from_x = x0 - hwhm_factor*hwhm
        to_x = x0 + hwhm_factor*hwhm
        if (hwhm/x_range*100 >= th): break
    y = eval(k_value, 0, from_x, to_x, npts)[2,:]
    x = np.linspace(from_x, to_x, npts)
    hwhm, x0 = get_half_width_half_maxima_and_x0(x, y)
    is_enough_spectra = q_factor_estimate(hwhm, x0, 1e11)
    # show_spectrum((10,10), x, y, k_value)
    Q = x0/hwhm

    return x, y, is_enough_spectra, Q


def show_spectrum(figsize, x, y, sign):
    plt.figure(figsize=figsize)
    plt.plot(x, y, 'r', lw=0.5)
    plt.scatter(x, y)
    plt.title(sign)
    #only for example
    plt.xlabel(r'${\omega d / 2\pi c }$')
    plt.ylabel('reflectance')
    plt.show()

    # save_result('jpg', 'M144figures', str(sign), data=None)
    plt.clf(); plt.close()


def show_R_map(figsize, x, y, R, ktype):
    plt.figure(figsize = figsize)
    im = plt.imshow(R.T, extent = (np.amin(x), np.amax(x), np.amax(y), np.amin(y)), cmap=cm.hot, norm=LogNorm(), aspect='auto', interpolation = 'none')
    cb = plt.colorbar(im)
    cb.set_label('reflectance')
    plt.ylabel(r'${\omega d / 2\pi c }$')
    if ktype == 1: plt.xlabel(r'${\theta}$')
    if ktype == 2: plt.xlabel(r'${k_x d/\pi}$')
    # ax1.set_xticks(np.arange(min(x)-0.001, max(x)+0.01, 0.1))
    sign_ax1 = ('d=%i'%(a)+'npts%i'%(npts)+'__pol_'+polar+'_epssph%f'%(epssph_re)+'_l_'+order+'_m_'+m)
    plt.title(sign_ax1)
    plt.gca().invert_yaxis()
    plt.show()
    # plt.clf(); plt.close()


def save_result(regime, folder_name, filename, data):
    if not (os.path.exists(folder_name)):
        os.makedirs(folder_name)

    if regime == 'jpg':
        plt.savefig(folder_name+'/'+filename+'.png')
        plt.close()

    if regime == 'txt':
        np.savetxt(folder_name+'/'+filename+'.txt', data, delimiter = ',')
    print('results were saved')


def multi_loadtxt(dir, filelist):
    output = ()
    for fname in filelist:
        out = np.loadtxt(dir+"/"+fname)
        output += (out,)
    return output

# figsize = (10,10)
plt.rcParams['font.size'] = '14'
is_fit_data_needed = 1
is_save_needed = 0
lmax = 5
# a = 400
a = 400
rmax = 7
s = 100
r_ratio = s/a
polar='S' # S or P
# epssph_re = 250.0
epssph_re = 50.0
epssph_im = 0.0000

is_multipole_type_selected = '1'
is_multipole_order_selected = '1'
is_m_projection_selected = '1'
type =   '1'
order =  '3'
m =      '0'

ktype = 1
R_min = 1e-20

# k_values = np.logspace(40, -9, num = 1, endpoint=True, base=2)
k_values = np.linspace(40, 1, 1)
# k_values = k_values/2
show = 1

k = 0

# for k_value in k_values:
    # k += 1
    # print(k, 'is calculating from', len(k_values))
    # print(k, 'is calculating from', len(theta_values))

    #M133
    # from_y = 0.452095418
    # to_y =   0.452095463

    # from_y = 0.5
    # to_y =   0.51

    #M103
    # from_y = 0.4520912
    # to_y =   0.4520984

    # from_y = 0.45
    # to_y = 0.46

    #M1X4 lossless
    # from_y = 0.5488589
    # to_y =   0.548858993

    #M1X4 lossy
    # from_y = 0.548
    # to_y = 0.549

    #M1X5
    # a_values = np.linspace(400, 400, 1)
    # for a in a_values:
    #     r_ratio = s/a
    #     # delta_ap = k_value/(100*kpts_ap1)
    #     # from_angle_param1 = (k_value - delta_ap)/2; to_angle_param1 = (k_value + delta_ap)/2
    #     from_y = 0.516; to_y = 0.518; npts = 100
    #     from_angle_param1 = 30; to_angle_param1 = 90.1; kpts_ap1 = 1
    #     from_angle_param2 = 0; to_angle_param2 = 90.1; kpts_ap2 = 1
    #     regime = 'freq_scan'
    #     x, y, R = calc_R_map(from_angle_param1, to_angle_param1, kpts_ap1, from_angle_param2, to_angle_param2, kpts_ap2, from_y, to_y, npts, R_min, regime)
    # #     if show:
    # #         show_R_map(figsize, x, y, R, ktype)
    # # plt.show()
    #     show_spectrum(figsize, x, y, 'test')
    #     # save_result('txt', 'rest maps', 'x', x)
    #     # save_result('txt', 'rest maps', 'y', y)
    #     # save_result('txt', 'rest maps', 'R', R.T)


# q_factors = []
# d_val = []
# from_angle_param1 = 40
# a_values = np.arange(250, 370, 1)
# from_y = 0.3; to_y = 0.5; npts = 300
# for a in a_values:
#     r_ratio = s/a
#     show = 0
#     # if is_fit_data_needed:
#     #         th = 25
#     #         x = np.linspace(from_y, to_y, npts)
#     #         x, y, is_enough_spectra, Q = find_spectrum(from_angle_param1, x, th)
#     #         q_factors.append(Q)
#     #         d_val.append(a)
#     #         if is_enough_spectra:
#     #             print('q is more than expected')
#     #             # break
#     #         if show:
#     #             show_spectrum(figsize, x, y, 'd_'+str(a)+'theta_'+str(k_value)+'_Q_'+str(Q))
#
# i = 1
# # a_values = np.arange(371, 450, 1)
# a_values = np.arange(465, 480, 0.05)
# from_y = 0.47; to_y = 0.6; npts = 300
# for a in a_values:
#     print(i, 'is calculating from', len(a_values))
#     i = i + 1
#     r_ratio = s/a
#     if is_fit_data_needed:
#         th = 25
#         x = np.linspace(from_y, to_y, npts)
#         x, y, is_enough_spectra, Q = find_spectrum(from_angle_param1, x, th)
#         q_factors.append(Q)
#         d_val.append(a)
#         if is_enough_spectra:
#             print('q is more than expected')
#             # break
#         # if show:
#             # show_spectrum(figsize, x, y, 'd_'+str(a)+'theta_'+str(k_value)+'_Q_'+str(Q))


# a_values = np.arange(451, 500, 1)
# from_y = 0.55; to_y = 0.7; npts = 300
# for a in a_values:
#     r_ratio = s/a
#     if is_fit_data_needed:
#         th = 25
#         x = np.linspace(from_y, to_y, npts)
#         x, y, is_enough_spectra, Q = find_spectrum(from_angle_param1, x, th)
#         q_factors.append(Q)
#         d_val.append(a)
#         if is_enough_spectra:
#             print('q is more than expected')
#             # break
#         if show:
#             show_spectrum(figsize, x, y, 'd_'+str(a)+'theta_'+str(k_value)+'_Q_'+str(Q))

# save_result('txt', '40_detailed', 'q', q_factors)
# save_result('txt', '40_detailed', 'd', d_val)
# plt.plot(d_val, q_factors)
# plt.show()


# plt.rcParams.update({
#     "font.family": "sans-serif",
#     "font.size": 14})

# fig, axs = plt.subplots(2, 2, figsize=(18,12))
# fig.tight_layout(pad=2, w_pad=0.5, h_pad=5)
# dirs = ['15', '30', '40', '50', '65']
# for dir in dirs:
#     x, y = multi_loadtxt(dir, ("d.txt", "q.txt"))
#     axs[0,0].plot(x, y, label=dir+' deg.', lw=2, alpha=1)
#     axs[0,0].set_yscale('log')
#     axs[0,0].legend()
#     axs[0,0].set_ylabel('Q factor')
#     axs[0,0].set_xlabel('d, nm')
#
#
# x, y = multi_loadtxt('40_detailed', ("d.txt", "q.txt"))
# axs[1,0].plot(x, y, color='red', lw=1, alpha=1)
# axs[1,0].set_yscale('log')
#
# x, y = multi_loadtxt('50_detailed', ("d.txt", "q.txt"))
# axs[0,1].plot(x, y, color='red', lw=1, alpha=1)
# axs[0,1].set_yscale('log')
#
# x, y = multi_loadtxt('65_detailed', ("d.txt", "q.txt"))
# axs[1,1].plot(x, y, color='red', lw=1, alpha=1)
# axs[1,1].set_yscale('log')
#
#
# plt.show()


    # if is_save_needed:
    #     c = 3e8
    #     d = a*1e-9
    #     x = x*2*np.pi*c/d
    #     spectra_data = np.array(list(zip(x, y)))
    #     # save_result('jpg', 'M105figures_2', str(k_value), data=None)
    #     save_result('txt', 'M105spectra', '103,k='+str(2*k_value), spectra_data)



#exmples for progress report
R_min = 1e-10
from_y = 0.5; to_y = 0.6; npts = 300
from_angle_param1 = 0; to_angle_param1 = 90.1; kpts_ap1 = 100
from_angle_param2 = 0; to_angle_param2 = 0; kpts_ap2 = 1
regime = 'theta_scan'
# x, y, R = calc_R_map(from_angle_param1, to_angle_param1, kpts_ap1, from_angle_param2, to_angle_param2, kpts_ap2, from_y, to_y, npts, R_min, regime)
# show_R_map((10,10), x, y, R, ktype)
if is_fit_data_needed:
    theta = 20
    th = 25
    x = np.linspace(from_y, to_y, npts)
    x, y, is_enough_spectra, Q = find_spectrum(theta, x, th)
    # if is_enough_spectra:
        # print('q is more than expected')
        # break
    if show:
        show_spectrum((10,10), x, y, 'theta='+str(theta))