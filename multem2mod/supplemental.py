import matplotlib.pyplot as plt
import numpy as np
import os
import subprocess
import matplotlib.cm as cm
from matplotlib.colors import LogNorm
import math as m
from matplotlib import rc
from matplotlib.ticker import ScalarFormatter, NullFormatter
from matplotlib.ticker import FormatStrFormatter
import matplotlib.ticker as mticker


def create_input(npts, ak1, ak2, zinf, zsup, polar, lmax, r_ratio, rmax, epssph_re, epssph_im, type, order,
                 is_multipole_type_selected, is_multipole_order_selected, is_m_projection_selected):


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
def find_spectrum(ap1, ap2, x, th):
    from_x = x[0]
    to_x = x[-1]
    npts = len(x)
    while True:
        y = eval(ap1, ap2, from_x, to_x, npts)[2,:]
        x = np.linspace(from_x, to_x, npts)
        # show_spectrum((10,10), x, y, ap1)
        x_range = to_x - from_x
        hwhm, x0 = get_half_width_half_maxima_and_x0(x, y)
        hwhm_factor = 2
        from_x = x0 - hwhm_factor*hwhm
        to_x = x0 + hwhm_factor*hwhm
        if (hwhm/x_range*100 >= th): break
    y = eval(ap1, ap2, from_x, to_x, npts)[2,:]
    x = np.linspace(from_x, to_x, npts)
    hwhm, x0 = get_half_width_half_maxima_and_x0(x, y)
    is_enough_spectra = q_factor_estimate(hwhm, x0, 1e11)
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


#TODO use kpts_ap2 for r = f(kx,ky)
def calc_R_map(ap1_start, ap1_end, kpts_ap1, ap2_start, ap2_end, kpts_ap2, from_y, to_y, npts, R_min):
    #regime: 'theta_scan' or 'phi_scan'
    #TODO autocheck for ap1 and ap2 size in case of regime
    #TODO dry for if regime
    ap1 = np.linspace(ap1_start, ap1_end, kpts_ap1)
    ap2 = np.linspace(ap2_start, ap2_end, kpts_ap2)

    x = ap1
    y = np.linspace(from_y, to_y, npts)
    R = np.empty((kpts_ap1, npts))
    for i in range(kpts_ap1):
        print(i+1, 'of', kpts_ap1)
        R[i,:] = eval(ap1[i], ap2, from_y, to_y, npts)[2,:]
        R[R<R_min] = R_min

    return x, y, R


def show_R_map(figsize, x, y, R, ktype):
    plt.figure(figsize = figsize)
    if ktype == 1: plt.xlabel(r'${\theta}$')
    if ktype == 2:
        x = x*2
        plt.xlabel(r'${k_x d/\pi}$')
    im = plt.imshow(R.T, extent = (np.amin(x), np.amax(x), np.amax(y), np.amin(y)), cmap=cm.hot, norm=LogNorm(), aspect='auto', interpolation = 'none')
    cb = plt.colorbar(im)
    cb.set_label('reflectance')
    plt.ylabel(r'${\omega d / 2\pi c }$')
    # ax1.set_xticks(np.arange(min(x)-0.001, max(x)+0.01, 0.1))
    # sign_ax1 = ('d=%i'%(a)+'npts%i'%(npts)+'__pol_'+polar+'_epssph%f'%(epssph_re)+'_l_'+order+'_m_'+m)
    # plt.title(sign_ax1)
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


def binary_data_processing(data, criterion, min_value, max_value):
    data = np.array(data)
    for i in range(np.size(data)):
        if (data[i] <= criterion):
            data[i] = min_value
        else:
            data[i] = max_value

    return data


#q factor at off-gamma point for M103
#collecting spectra

# lmax = 3
# a = 350
# rmax = 7
# s = 100
# r_ratio = s/a
# polar='S' # S or P
# epssph_re = 50.0
# epssph_im = 0.0000
# is_multipole_type_selected = '1'
# is_multipole_order_selected = '1'
# is_m_projection_selected = '1'
# type =   '1'
# order =  '3'
# m =      '0'
# ktype = 2
#
#
# ap1_start = 0.0; ap1_end = 0.42 + 0.001; npts_ap1 = 30
# ap2_start = 0; ap2_end = 0; npts_ap2 = 1
# from_y = 0.45; to_y = 0.455; npts_y = 300
# R_min = 1e-15
# #calculation of map of reflectance
# figsize = (10,10)
# # x, y, R = calc_R_map(ap1_start, ap1_end, npts_ap1, ap2_start, ap2_end, npts_ap2, from_y, to_y, npts_y, R_min)
# # show_R_map(figsize, x, y, R, ktype)
#
# ap1_values = np.linspace(0.05, 0.8, 30)
# # ap1_values = [0.808722]
# for ap1 in ap1_values:
#     if ktype == 2: ap1 = ap1/2
#     th = 25
#     freq = np.linspace(from_y, to_y, npts_y)
#     x, y, is_enough_spectra, Q = find_spectrum(ap1, 0, freq, th)
#     if is_enough_spectra:
#         print('q is more than expected')
#         print('ap1', ap1*2)
#         continue
#     if ktype == 2: ap1 = ap1*2
#     # show_spectrum(figsize, x, y, str(ap1))
#     c = 3e8
#     d = a*1e-9
#     x = x*2*np.pi*c/d
#     spectra_data = np.array(list(zip(x, y)))
#     save_result('txt', 'M103spectra', '103,k='+str(ap1), spectra_data)
#     plt.plot(x, y, 'r', lw=0.5)
#     plt.scatter(x, y)
#     save_result('jpg', 'off_gamma_pictures_test', '103,k='+str(ap1), '')

#-----------------------------------------------------------------------------------------------------------------------


# electric octopole N103 reflectance map
# lmax = 3
# a = 350
# rmax = 7
# s = 100
# r_ratio = s/a
# polar='P' # S or P
# epssph_re = 80.0
# epssph_im = 0.0000
# is_multipole_type_selected = '1'
# is_multipole_order_selected = '1'
# is_m_projection_selected = '1'
# type =   '0'
# order =  '3'
# m =      '0'
# ktype = 1
#
# ap1_start = 0; ap1_end = 89.99; npts_ap1 = 30
# ap2_start = 0; ap2_end = 0; npts_ap2 = 1
# from_y = 0.3; to_y = 0.6; npts_y = 300
# R_min = 1e-15
# figsize = (10,10)
# x, y, R = calc_R_map(ap1_start, ap1_end, npts_ap1, ap2_start, ap2_end, npts_ap2, from_y, to_y, npts_y, R_min)
# show_R_map(figsize, x, y, R, ktype)
#
# save_result('txt', 'N103', 'x', x)
# save_result('txt', 'N103', 'y', y)
# save_result('txt', 'N103', 'R', R.T)

#----------------------------------------------------------------------------------------------------------------------

#q factor at gamma point for M144 at gamma-m valley
#collecting spectra
# lmax = 4
# a = 400
# rmax = 7
# s = 100
# r_ratio = s/a
# polar='S' # S or P
# epssph_re = 50.0
# epssph_im = 0.0000
# is_multipole_type_selected = '1'
# is_multipole_order_selected = '1'
# is_m_projection_selected = '1'
# type =   '1'
# order =  '4'
# m =      '4'
#
# # ktype = 1
# # ap1_start = 0; ap1_end = 89.99; npts_ap1 = 10
# # ap2_start = 45; ap2_end = 45; npts_ap2 = 1
# from_y = 0.5; to_y = 0.7; npts_y = 500
# # R_min = 1e-15
# # figsize = (10,10)
# # x, y, R = calc_R_map(ap1_start, ap1_end, npts_ap1, ap2_start, ap2_end, npts_ap2, from_y, to_y, npts_y, R_min)
# # show_R_map(figsize, x, y, R, ktype)
# ktype = 2
# ap_values = np.logspace(-1.5, -4, base=2, num=15, endpoint=True)
# for ap1 in ap_values:
#     if ktype == 2: ap1 = ap1/2
#     #for gamma-m valley kx and ky should be equal
#     ap2 = ap1
#     th = 25
#     freq = np.linspace(from_y, to_y, npts_y)
#     x, y, is_enough_spectra, Q = find_spectrum(ap1, ap2, freq, th)
#     if is_enough_spectra:
#         print('q is more than expected')
#         break
#     if ktype == 2: ap1 = ap1*2
#     # show_spectrum(figsize, x, y, str(ap1))
#     c = 3e8
#     d = a*1e-9
#     x = x*2*np.pi*c/d
#     spectra_data = np.array(list(zip(x, y)))
#     save_result('txt', 'M144spectra_gamma', '103,k='+str(ap1), spectra_data)
#     plt.plot(x, y, 'r', lw=0.5)
#     plt.scatter(x, y)
#     save_result('jpg', 'gamma_pictures', '103,k='+str(ap1), '')


#-----------------------------------------------------------------------------------------------------------------------

# q factor for N103 collecting spectra

# lmax = 3
# a = 350
# rmax = 7
# s = 100
# r_ratio = s/a
# polar='P' # S or P
# epssph_re = 80.0
# epssph_im = 0.0000
# is_multipole_type_selected = '1'
# is_multipole_order_selected = '1'
# is_m_projection_selected = '1'
# type =   '0'
# order =  '3'
# m =      '0'
# ktype = 2
# #
# ap1_start = 0; ap1_end = 0.5 - 0.001; npts_ap1 = 50
# ap2_start = 0; ap2_end = 0; npts_ap2 = 1
# from_y = 0.433; to_y = 0.434; npts_y = 300
# R_min = 1e-15
# figsize = (10,10)
# # x, y, R = calc_R_map(ap1_start, ap1_end, npts_ap1, ap2_start, ap2_end, npts_ap2, from_y, to_y, npts_y, R_min)
# # show_R_map(figsize, x, y, R, ktype)
#
# # ap1_values = np.logspace(-4, -12, num = 20, endpoint=True, base=2) #gamma
# # ap1_values = np.linspace(0.75, 0.8, 50) #off-gamma
# ap1_values = [0.775]
# for ap1 in ap1_values:
#     if ktype == 2: ap1 = ap1/2
#     th = 25
#     freq = np.linspace(from_y, to_y, npts_y)
#     x, y, is_enough_spectra, Q = find_spectrum(ap1, 0, freq, th)
#     if is_enough_spectra:
#         print('q is more than expected')
#         break
#     if ktype == 2: ap1 = ap1*2
#     # show_spectrum(figsize, x, y, str(ap1))
#     c = 3e8
#     d = a*1e-9
#     x = x*2*np.pi*c/d
#     spectra_data = np.array(list(zip(x, y)))
#     # save_result('txt', 'N103spectra_gamma', '103,k='+str(ap1), spectra_data)
#     save_result('txt', 'N103spectra', '103,k='+str(ap1), spectra_data)
#     plt.plot(x, y, 'r', lw=0.5)
#     plt.scatter(x, y)
#     save_result('jpg', 'gamma_pictures', '103,k='+str(ap1), '')

#-----------------------------------------------------------------------------------------------------------------------
#q factor vs kx visualisation
# plt.rcParams.update({
#     "font.family": "sans-serif",
#     "font.size": 14})
# plt.rcParams['ps.useafm'] = True
# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# plt.rcParams['pdf.fonttype'] = 42
# #
# main_dir = 'supplemental_material'
# dir = 'M144_gamma_x_valley'
# fig,ax = plt.subplots(figsize=(10,10))
# x, y = multi_loadtxt(main_dir+'/'+dir, ("RLns-fano.txt", "Q0s-fano.txt"))
# log_x, log_y = np.log(x), np.log(y)
# curve_fit = np.polyfit(log_x, log_y, 1)
# q_fitted = curve_fit[0]*log_x + curve_fit[1]
# q_fitted = np.exp(q_fitted)
# # plt.title(curve_fit[0])
# ax.plot(x, q_fitted, color='red', lw=1, alpha=0.5)
# ax.plot(x, y, marker='o', linestyle='None', mfc='None', mec='black', mew=2, ms=8)
# ax.set_xscale('log')
# ax.set_yscale('log')
# ax.tick_params(which='both', direction='in')
# ax.legend(['approximation', 'simulation'], loc=3)
# ax.set_ylim([1e8, 1e11])
# # ax.set_xlim([0.12, 0.4])
# ax.xaxis.set_major_formatter(ScalarFormatter())
# ax.xaxis.set_minor_formatter(NullFormatter())
# ax.xaxis.set_ticklabels([])
# ax.yaxis.set_ticklabels([])
# # plt.show()
# plt.savefig('/home/ashalev/Pictures/'+dir+'.eps')
# plt.close()
#-----------------------------------------------------------------------------------------------------------------------


# # supplm.material figure #TODO fill
# # M103 all band
# plt.rcParams.update({
#         "font.family": "sans-serif",
#         "font.size": 14})
# plt.rcParams['ps.useafm'] = True
# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# plt.rcParams['pdf.fonttype'] = 42
# fig,ax = plt.subplots(figsize=(10,10))
# dir = 'supplemental_material/M103_off_gamma'
# x1, y1 = multi_loadtxt(dir, ("RLns-fano.txt", "Q0s-fano.txt"))
# # x2, y2 = multi_loadtxt(main_folder2, ("RLns-fano.txt", "Q0s-fano.txt"))
# # x = np.concatenate((x1, x2), axis=0)
# # y = np.concatenate((y1, y2), axis=0)
# # idx = x.argsort()
# # x = x[idx]
# # y = y[idx]
# # save_result('txt', 'test', 'RLns-fano', x)
# # save_result('txt', 'test', 'Q0s-fano', y)
# # ax.plot(x1, y1, color='black', linestyle='--', lw=2, alpha=1)
# ax.plot(x1, y1, marker='o', linestyle='None', mfc='None', mec='black', mew=2, ms=8)
#
#
# # t = np.linspace(0.0001169,1,1000)
# # ax.plot(t*(0.452*2),(19e1*(np.abs(2/(3*t - 15*t* (1 - t**2)))))**2, color='red', lw=2)
# ax.tick_params(which='both', direction='in')
# #off gamma
# t = np.linspace(0.00001,1,10000)
# ax.plot(t*(0.4521*2),(13e-7*np.abs(1e8*2/(3*t - 15*t* (1 - t**2))))**2,  color='red', lw=1, alpha=0.5)
# ax.set_xlim([0.8, 0.82])
# ax.legend(['numerical data', 'formula'], loc=3)
# ax.xaxis.set_major_formatter(ScalarFormatter())
# ax.xaxis.set_minor_formatter(NullFormatter())
# ax.yaxis.set_major_formatter(ScalarFormatter())
# ax.yaxis.set_minor_formatter(NullFormatter())
# ax.xaxis.set_ticklabels([])
# ax.yaxis.set_ticklabels([])
# ax.set_yscale('log')
# # plt.show()
# plt.savefig('/home/ashalev/Desktop/M103_off_gamma.eps', format='eps')
# plt.close()


#-----------------------------------------------------------------------------------------------------------------------


# supplm.material figure 3
# fig,ax = plt.subplots(figsize=(10,10))
# N103 all band ---------------------------------------------------------------------------------------------
# dir = 'supplemental_material/N103_qfactor_all_band'
# x, y = multi_loadtxt(dir, ("RLns-fano.txt", "Q0s-fano.txt"))
# ax.plot(x, y, color='black', linestyle='--', lw=2, alpha=1)
# t = np.linspace(0.0001169,1,1000)
# ax.plot(t*(0.433*2),(90e1*(np.abs(2/(3*t - 15*t* (1 - t**2)))))**2, color='red', lw=2)
# ax.tick_params(which='both', direction='in')
# N103 off gamma ----------------------------------------------------------------------------------------------
# dir = 'supplemental_material/N103_off_gamma'
# x, y = multi_loadtxt(dir, ("RLns-fano.txt", "Q0s-fano.txt"))
# t = np.linspace(0.00001,1,10000)
# ax.plot(x, y, marker='o', linestyle='None', mfc='None', mec='black', mew=2, ms=8)
# ax.plot(t*(0.4333*2),(65e-7*np.abs(1e8*2/(3*t - 15*t* (1 - t**2))))**2,  color='red', lw=1, alpha=0.5)
# ax.set_xlim([0.75, 0.8])
# ax.legend(['numerical data', 'formula'], loc=3)
#--------------------------------------------------------------------------------------------------------------
# ax.xaxis.set_major_formatter(ScalarFormatter())
# ax.xaxis.set_minor_formatter(NullFormatter())
# ax.yaxis.set_major_formatter(ScalarFormatter())
# ax.yaxis.set_minor_formatter(NullFormatter())
# ax.xaxis.set_ticklabels([])
# ax.yaxis.set_ticklabels([])
# ax.set_yscale('log')
# # plt.show()
# plt.savefig('/home/ashalev/Desktop/N103_off_gamma.eps', format='eps')
# plt.close()


#--------------------------------------------------------------------------------------------------------------
# polarization maps. what figure?
# main_dir = 'supplemental_material/polarization_maps'
# dir = 'N103'
# kx, ky= multi_loadtxt(main_dir+'/'+dir, ("kx.txt", "ky.txt"))
# fig,ax = plt.subplots(figsize=(10,10))
# if dir == 'N103':
#     F2, F3 = multi_loadtxt(main_dir+'/'+dir, ("F_2.csv", "F_3.csv"))
#     Cx, Cy = multi_loadtxt(main_dir+'/'+dir, ("Cx.csv", "Cy.csv"))
# if dir == 'M103':
#     theta = np.arcsin(np.sqrt(kx**2 + ky**2))
#     phi = np.angle(kx+1j*ky)
#     F2, F3 = multi_loadtxt(main_dir+'/'+dir, ("F_2.csv", "F_3.csv"))
#     Cx = F2*np.cos(theta)*np.cos(phi) - F3*np.sin(phi)
#     Cy = F2*np.cos(theta)*np.sin(phi) + F3*np.cos(phi)
#
# step = 1
# amp = np.sqrt(Cx**2 + Cy**2)
# Cx, Cy = Cx/amp, Cy/amp
# amp = binary_data_processing(F3, 0, -1, 1)
# ax.quiver(kx[::step], ky[::step], Cx[::step], Cy[::step], amp[::step], cmap='rainbow', width=0.004, pivot="mid", headwidth=3.0, headlength=3.0, headaxislength=2)
#
# ax.xaxis.set_major_formatter(ScalarFormatter())
# ax.xaxis.set_minor_formatter(NullFormatter())
# ax.yaxis.set_major_formatter(ScalarFormatter())
# ax.yaxis.set_minor_formatter(NullFormatter())
# ax.xaxis.set_ticklabels([])
# ax.yaxis.set_ticklabels([])
#
# # plt.show()
#
# plt.savefig('/home/ashalev/Pictures/'+dir+'.eps', format='eps')
# plt.close()

#--------------------------------------------------------------------------------------------------------------


#figure 2


# def plot_qfactor_all_range(dir, figsize, is_analytics_needed, is_subplot):
#     fig, ax = plt.subplots(figsize=figsize)
#     x, y = multi_loadtxt(dir, ("RLns-fano.txt", "Q0s-fano.txt"))
#     ax.plot(x, y, color='black', linestyle='--', lw=3, alpha=1)
#     if is_analytics_needed:
#         t = np.linspace(0.00001169,1,10000)
#         #M103
#         # ax.plot(t*(0.452*2),(19e1*(np.abs(2/(3*t - 15*t*(1 - t**2)))))**2, color='red', lw=2)
#         #N103
#         ax.plot(t*(0.433*2),(90e1*(np.abs(2/(3*t - 15*t* (1 - t**2)))))**2, color='red', lw=2)
#     ax.set_yscale('log')
#     ax.tick_params(which='both', direction='in', length=5, width=1)
#     # ax.legend(['simulation', 'analytics'], loc=9)
#     ax.set_ylim([1e5, 2e11])
#     # ax.set_ylabel(r'${Q}$')
#     # ax.set_xlabel(r'${k_x d/\pi}$')
#     ax.yaxis.set_ticks(np.logspace(5, 11, 4))
#     ax.xaxis.set_ticks(np.arange(0.0, 1.0+0.1, 0.2))
#
#     ax.xaxis.set_major_formatter(ScalarFormatter())
#     ax.xaxis.set_minor_formatter(NullFormatter())
#     ax.yaxis.set_major_formatter(ScalarFormatter())
#     ax.yaxis.set_minor_formatter(NullFormatter())
#     ax.xaxis.set_ticklabels([])
#     ax.yaxis.set_ticklabels([])
#
#     # plt.show()
#
# def plot_qfactor_gamma_point(dir, figsize, is_subplot):
#     fig, ax = plt.subplots(figsize=figsize)
#     x, y = multi_loadtxt(dir, ("RLns-fano.txt", "Q0s-fano.txt"))
#     log_x, log_y = np.log(x), np.log(y)
#     curve_fit = np.polyfit(log_x, log_y, 1)
#     xx = np.linspace(np.min(x)-0.1*(np.max(x) - np.min(x)), np.max(x)+0.1*(np.max(x) - np.min(x)), 1000)
#     log_xx = np.log(xx)
#     q_fitted = curve_fit[0]*log_xx + curve_fit[1]
#     q_fitted = np.exp(q_fitted)
#     ax.plot(xx, q_fitted, color='red', lw=2, alpha=1)
#     ax.plot(x, y, marker='o', linestyle='None', mfc='none', mec='black', mew=2, ms=8)
#     ax.set_ylim([1e7, 2e11])
#     ax.set_xlim([3e-4, 4e-2])
#     # ax.set_ylabel('Q')
#     # ax.set_xlabel(r'${k_x d/\pi}$')
#     # ax.legend(['simulation', 'analytics'], loc=3)
#     ax.set_xscale('log')
#     ax.set_yscale('log')
#     ax.tick_params(which='both', direction='in', length=5, width=1)
#     ax.xaxis.set_major_formatter(ScalarFormatter())
#     ax.xaxis.set_minor_formatter(NullFormatter())
#     ax.yaxis.set_major_formatter(ScalarFormatter())
#     ax.yaxis.set_minor_formatter(NullFormatter())
#     ax.xaxis.set_ticklabels([])
#     ax.yaxis.set_ticklabels([])
#     # plt.show()
#
#
# def plot_qfactor_offgamma_point(dir, figsize, is_subplot):
#     fig, ax = plt.subplots(figsize=figsize)
#     x, y = multi_loadtxt(dir, ("RLns-fano.txt", "Q0s-fano.txt"))
#
#     t = np.linspace(0.01,1,30000)
#     #M103
#     # ax.plot(t*(0.4521*2),(13e-7*np.abs(1e8*2/(3*t - 15*t* (1 - t**2))))**2,  color='red', lw=2, alpha=1)
#     #N103
#     ax.plot(t*(0.4333*2),(65e-7*np.abs(1e8*2/(3*t - 15*t* (1 - t**2))))**2,  color='red', lw=2, alpha=1)
#     ax.plot(x, y, marker='o', linestyle='None', mfc='none', mec='black', mew=2, ms=8)
#     # ax.set_xlim([0.8, 0.82])
#     ax.set_xlim([0.75, 0.8])
#     ax.set_xscale('log')
#     ax.set_yscale('log')
#     # ax.legend(['simulation', 'analytics'], loc=1)
#
#     ax.set_ylim([5e5, 2e12])
#     # ax.yaxis.set_ticks(np.logspace(3, 11, 5))
#     ax.tick_params(which='both', direction='in', length=5, width=1)
#     # plt.show()
#     ax.yaxis.set_ticks(np.logspace(6, 12, 4))
#
#     ax.xaxis.set_major_locator(mticker.FixedLocator([0.75, 0.76, 0.77, 0.78, 0.79, 0.8]))
#     ax.xaxis.set_major_formatter(ScalarFormatter())
#     ax.xaxis.set_minor_formatter(NullFormatter())
#     ax.xaxis.set_ticklabels([])
#     ax.yaxis.set_ticklabels([])
#
#
#
# def multipole_custom_plot(x1, y1, x2, y2, ax=None, first_plt_kwargs={}, second_plt_kwargs={}):
#     if ax == None:
#         ax = plt.gca()
#     ax.plot(x1, y1, **first_plt_kwargs)   #plot numerical data
#     ax.plot(x2, y2, **second_plt_kwargs)  #plot analytics
#     return (ax)
#
# plt.rcParams.update({
#     "font.family": "sans-serif",
#     "font.size": 20})

#------------------------------------------------------------------------------------
# #supl Fig. 2a
# fig, axs = plt.subplots(1, 3, figsize=(15.65,8.65))
# #numerical data
# dir = 'supplemental_material/M103_qfactor_all_band'
# x1, y1 = multi_loadtxt(dir, ("RLns-fano.txt", "Q0s-fano.txt"))
# num_data_plt_params = {'color':'black', 'linestyle':'--', 'linewidth': 3}
# #analytics
# t = np.linspace(0.00001169,1,10000)
# x2, y2 = t*(0.452*2), (19e1*(np.abs(2/(3*t - 15*t*(1 - t**2)))))**2
# analytics_plt_params = {'color':'red', 'lw':2}
#
# multipole_custom_plot(x1, y1, x2, y2, ax=axs[0], first_plt_kwargs=num_data_plt_params, second_plt_kwargs=analytics_plt_params)
# #TODO hide all ax settings in function
# axs[0].set_yscale('log')
# axs[0].tick_params(which='both', direction='in', length=5, width=1)
# axs[0].legend(['simulation', 'analytics'], loc=3)
# axs[0].set_ylim([1e3, 2e11])
# axs[0].set_ylabel(r'${Q}$')
# axs[0].set_xlabel(r'${k_x d/\pi}$')
# axs[0].yaxis.set_ticks(np.logspace(3, 11, 5))
# axs[0].xaxis.set_ticks(np.arange(0.0, 1.0+0.1, 0.2))
# plt.show()
#------------------------------------------------------------------------------------

#---- fig2 suppl -----------------------------------------------------------
# dir = 'supplemental_material/M103_qfactor_all_band'
# figsize = (10,10)
# is_analytics_needed = True
# is_subplot = False
# # plot_qfactor_all_range(dir, figsize, is_analytics_needed, is_subplot)
# dir = 'supplemental_material/M103_qfactor_gamma'
# # plot_qfactor_gamma_point(dir, figsize, is_subplot)
# dir = 'supplemental_material/M103_qfactor_offgamma'
# plot_qfactor_offgamma_point(dir, figsize, is_subplot)
#
#
#
# fn = 'M103qfactor_offgamma_wo_legend'
# plt.savefig('/home/ashalev/Pictures/'+fn+'.eps', format='eps')
# # plt.savefig('/home/ashalev/Pictures/'+fn+'.png', format='png')
# plt.close()
# plt.show()
# plt.close()
#------------------------------------------------------------------------------------

#----------- N103 ------------------------
# dir = 'supplemental_material/N103_qfactor_all_band'
# figsize = (10,10)
# is_analytics_needed = True
# is_subplot = False
# plot_qfactor_all_range(dir, figsize, is_analytics_needed, is_subplot)


# dir = 'supplemental_material/N103_qfactor_gamma'
# figsize = (10,10)
# is_subplot = False
# # plot_qfactor_gamma_point(dir, figsize, is_subplot)
# #
# dir = 'supplemental_material/N103_qfactor_offgamma'
# figsize = (10,10)
# is_subplot = False
# plot_qfactor_offgamma_point(dir, figsize, is_subplot)
#
# fn = 'N103qfactor_offgamma_wo_legend'
# plt.savefig('/home/ashalev/Pictures/'+fn+'.eps', format='eps')
# # plt.savefig('/home/ashalev/Pictures/'+fn+'.png', format='png')
# plt.close()
# # plt.show()




# plot_qfactor_offgamma_point(dir)
# construct_plots(nrows, ncols)



#fig 3 main paper
lmax = 5
a = 400
rmax = 16
s = 100
r_ratio = s/a
polar='S' # S or P
epssph_re = 220.0
epssph_im = 0.0000
is_multipole_type_selected = '1'
is_multipole_order_selected = '1'
is_m_projection_selected = '1'
type =   '1'
order =  '5'
m =      '0'
ktype = 1

ap1_start = 0; ap1_end = 90.1; npts_ap1 = 100
# ap1_start = 0; ap1_end = (1. - 0.001)/2; npts_ap1 = 100
ap2_start = 0.0; ap2_end = ap2_start; npts_ap2 = 1
from_y = 0.495; to_y = 0.51; npts_y = 100
R_min = 1e-15
figsize = (10,10)
x, y, R = calc_R_map(ap1_start, ap1_end, npts_ap1, ap2_start, ap2_end, npts_ap2, from_y, to_y, npts_y, R_min)
show_R_map(figsize, x, y, R, ktype)
