import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import scipy.optimize
from tqdm import tqdm
from scipy import signal
import peakutils


def gaussian(x, x0, gamma):
    return np.exp(-(x-x0)**2/(2*gamma**2))


def lorentzian(x, I, x0, gamma, bg_level):
    return I*(gamma**2/ ((x - x0)**2 + gamma**2)) + bg_level
    # return I* (gamma**2
    #         /( (G/2.)**2 + (x-fl)**2 )
    #         +bg_level)
# return ((I**2+0.1)**0.5 /(1+q**2)
#         *(q*G/2. + x -fl  )**2
#         /( (G/2.)**2 + (x-fl)**2 )
#         +bg_level#+x*slope# + slope2*x**2
#         )


def multi_savetxt(output_dir, in_var):
    for pair in in_var:
        np.savetxt(output_dir+"/"+pair[1], pair[0])


def multi_loadtxt(dir, filelist):
    output = ()
    for fname in filelist:
        out = np.loadtxt(dir+"/"+fname)
        output += (out,)
    return output

def set_pbar(total, desc):
    return tqdm(total=total,
                desc=desc,
                file=sys.stdout,
                bar_format='{l_bar}{bar}| elapsed: {elapsed} remaining: {remaining}')

def read_data_from_files(dir_name, is_all_spec_needed):
    data_names = []
    data_h=[]
    for filename in os.listdir(dir_name):
        path = dir_name+"/"+filename
        if path.count('/') > 1: continue
        try:
            tt=float(filename[6:-4])
            # print(tt)
        except:
            continue
        data_names.append(filename)
        data_h.append(tt)
    # data_names.sort()
    data_val_names = np.array([[x, y] for x,y in sorted(zip(data_h,data_names))]).T
    data_h = [float(x) for x in data_val_names[0,:]]
    data_names = data_val_names[1,:]
    # print(data_val_names)
    # data_h = sorted(data_h)

    # print(data_names)
    # print(data_h)
    if not is_all_spec_needed:
        data_names = data_names[from_n:to_n]
        data_h = data_h[from_n:to_n]
    print(data_names, data_h)
    data_list = []
    pbar = set_pbar(total=len(data_names), desc='Reading data        ')
    for name in data_names:
        print(name)
        data = np.loadtxt(dir_name+"/"+name, comments=['#','!'], delimiter=',')
        # data[:,2] = 1. - data[:,2]
        data_list.append(data)
        pbar.update()
    pbar.close()

    return np.array(data_list), np.array(data_names), np.array(data_h)

def get_freq_estimate(data_list, data_h):
    #find the x, y, h range of all input data
    plt.figure('estimate')
    ax = plt.gca()
    x1 = data_list[0][:,0];        y1 = data_list[0][:,1];     h1 = data_h[0]
    x2 = data_list[-1][:,0];       y2 = data_list[-1][:,1];    h2 = data_h[-1]
    idx1, peak_widths, idx_widths = get_estimated_peaks_and_widths(x1, y1, ax)
    idx2, peak_widths, idx_widths = get_estimated_peaks_and_widths(x2, y2, ax)
    # plt.show()
    plt.clf(); plt.close()
    return idx1, idx2

def get_estimated_peaks_and_widths(x, y, plt_ax):
    # span = int(x.size/20) # Moving average size (to substract as background)
    span = int(x.size/10) # Moving average size (to substract as background)
    # make span to be even number for correct padding and convolution
    span = span+1 if span%2 == 0 else span
    #peak_thres = 0.24
    peak_thres = 0.05
    #peak_thres = 0.253

    # span = 181
    # peak_thres = 0.25
    #scipy filter
    b, a = signal.butter(1, 0.06)
    y_butter = signal.filtfilt(b, a, y)
    # Analog of Matalb smooth - moving average filter.
    y_pad = np.pad(y, int(span/2), 'edge')
    y_aver = np.convolve(y_pad, np.ones((span,))/span, mode='valid')
    # print(len(y), len(y_aver))
    y_filtered=y_butter-y_aver # filter some noise and background

    idxs = peakutils.indexes(y_filtered,
                             thres=peak_thres,
                             #min_dist=51
                             )
    max_try = 10000
    test_thres = peak_thres
    for i in range(max_try):
        if len(idxs) < 2: break
        test_thres *= 1.002
        idxs = peakutils.indexes(y_filtered,
                                 thres=test_thres,
                                 #min_dist=51
                                 )
        # print(i, test_thres, len(idxs))

    #get estimates of peak widths
    sign_switch = np.pad((np.diff(np.sign(y_filtered)) != 0)*1, (0,1), 'constant')
    peak_width_mask = np.zeros(sign_switch.shape)
    peak_widths = []
    idx_widths = []
    for i in idxs:
        j=0
        width = 0.
        idx_width = 0
        factor = 3
        while (j*factor+i+1 < len(x)):
            j+=1
            if sign_switch[i+j] == 0: continue
            if i-j*factor < 0: continue
            idx_width = int(j*factor*1)
            peak_width_mask[i-idx_width] = i
            peak_width_mask[i+idx_width] = i
            width = x[i+idx_width]-x[i-idx_width]
            break
        peak_widths.append(width)
        idx_widths.append(idx_width)
    # print(peak_widths)
    x_switch = x[peak_width_mask>0]
    y_switch = y_filtered[peak_width_mask>0]

    plt_ax.plot(x, y_filtered, lw=0.3, color = "black")
    plt_ax.plot(x, y, lw=0.3, marker='o', markersize=0.8, alpha=0.5, color = "blue")
    # plt_ax.plot(x, y_aver, lw=0.5, color = "cyan")
    plt_ax.scatter(x[idxs],y_filtered[idxs], s=30, alpha=0.8, c='r')
    plt_ax.scatter(x[idxs],y[idxs], s=30, alpha=0.5, c='r')
    plt_ax.scatter(x_switch,y_switch, s=30, alpha=0.5, c='b')
    #plt.show()
    return idxs[0], peak_widths, idx_widths

def lorentz_curve_estimate(x, y):
    idx = [np.argmin(y),np.argmax(y)]
    ymin, ymax = y[idx[0]], y[idx[1]]
    # print('y_max', ymax)
    # print('y_min', ymin)
    bg_level = ymin
    I = ymax-ymin
    # print('I_estimated = ', I)
    x0 = x[idx[1]]
    # print('x0',x0)
    idx_hwhm = np.abs(y - 0.5*I).argmin()
    # print('test', idx_hwhm)
    # print('y_hwhm', y[idx_hwhm])
    # print('x_hwhm', x[idx_hwhm])
    gamma = np.abs(x0 - x[idx_hwhm])
    # print('gamma', gamma)

    return I, x0, gamma, bg_level


def make_lorentzian_fit(data):
    x_data = data[:,0]
    # y_range = np.max(data[:,1]-np.min(data[:,1]))
    # y_min = np.min(data[:,1])
    # y_data = (data[:,1]-y_min)/y_range
    y_data = data[:,1]
    # x_data = data[:,0]-np.abs(data[0,0]+data[-1,0])/2
    # y_range = np.max(data[:,1]-np.min(data[:,1]))
    # y_min = np.min(data[:,1])
    # y_data = (data[:,1]-y_min)/y_range
    xfit = np.linspace(x_data[0], x_data[-1],300)
    # xfit = xfit + np.abs(data[0,0]+data[-1,0])/2

    try:
        # I = 1.0
        # I, x0, gamma = lorentz_curve_estimate(x_data, y_data, I)
        I, x0, gamma, bg_level = lorentz_curve_estimate(x_data, y_data)
        print('I_est', I)
        popt, perr = scipy.optimize.curve_fit(lorentzian, x_data, y_data, p0=[I, x0, gamma, bg_level])
        I, x0, gamma, bg_level  = tuple(popt)

        sigma_data = 1/lorentzian(x_data, I, x0, 2*gamma, bg_level)
        # sigma_data = 1/gaussian(x_data, x0, 2*gamma)
        plt.plot(x_data, 1/sigma_data, lw=1, alpha = 0.4, color='r', ls='--')
        popt, perr = scipy.optimize.curve_fit(lorentzian, x_data, y_data, p0=[I, x0, gamma, bg_level], sigma=sigma_data)
        I, x0, gamma, bg_level  = tuple(popt)

        isCF = True
        print("I:",I, "x0:",x0, "gamma:", gamma, "bg_level:", bg_level)
    except:
        I, x0, gamma, bg_level = 0.,0.,0.,0.
        isCF = False
    yfit = lorentzian(xfit, I, x0, gamma, bg_level)


    # print("Fit     :  ",I,x0, gamma)
    return x_data, y_data, xfit, yfit, I, x0, gamma, bg_level, isCF, sigma_data


q_factor = []; angle_param = []
d = 350e-9
c = 3e8
is_all_spec_needed = 1
from_n = 0; to_n = 15
spectra_dir = 'M133spectra'
output_dir = spectra_dir+"_output"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
data_list, data_names, data_h = read_data_from_files(spectra_dir, is_all_spec_needed)

# idx1, idx2 = get_freq_estimate(data_list, data_h)
#TODO what is idx1, idx2??
# print(idx1)
# print(idx2)

for i in range(len(data_list)):
    print('-----------------------------------------------')
    print(i, 'of', len(data_list))
    print('-----------------------------------------------')

    plt.figure('fit',figsize=(6,1*(len(data_list))))
    data = data_list[i]
    x_data = data[:,0]
    y_data = data[:,1]
    # xdd, ydd, xfit, yfit, I, x0, gamma, isCF = make_lorentzian_fit(data)

    # I_est, x0_est, gamma_est, bg_level = lorentz_curve_estimate(x_data, y_data)
    # plt.plot(x_data, lorentzian(x_data, I_est, x0_est, gamma_est, bg_level), 'k')


    xdd, ydd, xfit, yfit, I, x0, gamma, bg_level, isCF, sigma_data = make_lorentzian_fit(data)
    # plt.scatter(x_data, y_data)
    # plt.plot(x_data, lorentzian(x_data, I, x0, gamma, bg_level), '-r', linewidth = 2)
    # plt.plot(xfit, yfit, '-g')
    # plt.plot(x_data, lorentzian(x_data, (I+I_est)/2, (x0+x0_est)/2, (gamma+gamma_est)/2), '-r', linewidth = 2)
    # plt.show()

    if isCF:
        yfit = lorentzian(xfit, I, x0, gamma, bg_level)

    step = -1.7
    plt.plot(xfit, yfit+step*i, lw=1.5, alpha = 0.7)
    plt.plot(x_data, 1/sigma_data+step*i, lw=0.5, alpha = 0.4, color='black', ls='--')
    plt.plot(x_data, y_data+step*i, lw=0.5, alpha = 0.5, marker=',', markersize=0.1)

    angle_param.append(data_h[i])
    FWHM = 2*gamma
    q = x0/FWHM
    q_factor.append(q)
    # RLs.append(#cyl_R/
    #     (data_h[i]#+0.15
    #               )); QQs.append(q); QQs2stdev.append(q_stdev); Q0s.append(np.abs((x_mid+fl)/G))
    # RLns.append(#cyl_R/
    #     (data_h[i]));

multi_savetxt(output_dir,((angle_param, "angle_param.txt"),(q_factor, "q_factor.txt")))

plt.figure('fit')
# plt.title(r"$\frac{I_0}{1+q^2} \frac{(qG/2 + x -f_l)^2 }{ (G/2)^2 + (x-f_l)^2 } +I_{bg}$", fontsize=10)
# plt.title(model, mode="plain", fontsize=8)
# plt.ylim(-22,2)
plt.savefig(output_dir+'/fano_fit.pdf')
plt.clf(); plt.close()
ap,q = multi_loadtxt(output_dir, ("angle_param.txt", "q_factor.txt"))
plt.scatter(ap,q)
plt.plot(ap, 1e6/ap**4, 'r')
# plt.legend()
plt.xscale('log')
plt.yscale('log')
plt.xlabel('kxd/pi')
plt.ylabel('Q')
plt.savefig('nq.pdf')
plt.clf(); plt.close()