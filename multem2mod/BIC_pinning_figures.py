import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter, NullFormatter
from scipy.interpolate import interp2d
from PIL import Image
from matplotlib import cm
from matplotlib.colors import ListedColormap
from matplotlib import rc



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


plt.rcParams.update({
    "font.family": "sans-serif",
    "font.size": 14})

plt.rcParams['ps.useafm'] = True
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
plt.rcParams['pdf.fonttype'] = 42

#TODO new coloprmap

# rainbow = cm.get_cmap('rainbow', 256)
# newcolors = rainbow(np.linspace(0, 1, 256))
# red = np.array([256/256, 0/256, 0/256, 1])
# blue = np.array([0/256, 0/256, 256/256, 1])
# newcolors[:100, :] = blue
# newcolors[0:100, :] = red
# newcmp = ListedColormap(newcolors)

# vsh_105 = Image.open('/home/ashalev/Pictures/M105_mod.png')
# vsh_105 = np.array(vsh_105).astype(np.float) / 255
# vsh_103 = Image.open('/home/ashalev/Pictures/M103.png')
# vsh_103 = np.array(vsh_103).astype(np.float) / 255
# vsh_133 = Image.open('/home/ashalev/Pictures/M133_mod.png')
# vsh_133 = np.array(vsh_133).astype(np.float) / 255
# vsh_144 = Image.open('/home/ashalev/Pictures/M144_mod.png')
# vsh_144 = np.array(vsh_144).astype(np.float) / 255

fig, axs = plt.subplots(2, 4, figsize=(15.65,8.65), sharey='row')
fig.tight_layout(pad=2, w_pad=0.5, h_pad=5)
# main_folder = 'q_factor_results/13_09_results(fano)_double_precision/eps50_d350_a100/'
main_folder = 'q_factor_results/13_09_results(fano)_double_precision/new_try/'

dirs = ['M103spectra_output', 'M105spectra_output', 'M133spectra_output', 'M144spectra_output']
i = 0
j = 0
for dir in dirs:
    x, y = multi_loadtxt(main_folder+dir, ("RLs-fano.txt", "Q0s-fano.txt"))
    log_x, log_y = np.log(x), np.log(y)
    curve_fit = np.polyfit(log_x, log_y, 1)
    q_fitted = curve_fit[0]*log_x + curve_fit[1]
    q_fitted = np.exp(q_fitted)
    axs[i,j].plot(x, q_fitted, color='red', lw=1, alpha=0.5)
    axs[i,j].plot(x, y, marker='o', linestyle='None', mfc='None', mec='black', mew=2, ms=8)
    # axs[i,j].set_title(r'$'+dir[:4]+'$')
    axs[i,j].set_xscale('log')
    axs[i,j].set_yscale('log')
    # power = str(round(curve_fit[0]))
    # axs[i,j].legend([r'$\sim 1/k_x^{'+power+'}$', r'$multem$'], loc=3)
    # axs[i,j].set_xlabel(r'${k_x d/\pi}$')
    axs[i,j].tick_params(which='both', direction='in')
    axs[i,j].xaxis.set_ticklabels([])
    axs[i,j].yaxis.set_ticklabels([])
    j += 1

axs[0,2].xaxis.set_major_formatter(ScalarFormatter())
axs[0,2].xaxis.set_minor_formatter(NullFormatter())
axs[0,3].xaxis.set_major_formatter(ScalarFormatter())
axs[0,3].xaxis.set_minor_formatter(NullFormatter())
axs[0,2].xaxis.set_ticklabels([])
axs[0,3].xaxis.set_ticklabels([])

# axs[0,0].figure.figimage(vsh_103, 410, 680)
# axs[0,1].figure.figimage(vsh_105, 850, 680)
# axs[0,2].figure.figimage(vsh_133, 1205, 680)
# axs[0,3].figure.figimage(vsh_144, 1635, 680)
axs[0,0].legend(['approximation', 'simulation'], loc=3)
# axs[0,0].set_ylabel(r'$Q$')
# axs[0,2].xaxis.set_major_formatter(ScalarFormatter())
# axs[0,2].xaxis.set_minor_formatter(NullFormatter())
# axs[0,2].set_xticks([2e-2, 8e-2])
# axs[0,3].xaxis.set_major_formatter(ScalarFormatter())
# axs[0,3].xaxis.set_minor_formatter(NullFormatter())
# axs[0,3].set_xticks([0.12, 0.5])
# for j in range(4):
#     axs[j].tick_params(which='both', direction='in')

#creating a circle
# rad = 1.0
# circle1 = plt.Circle((0, 0), rad, color='black', fill=False)


main_folder = 'new_4fig'
kx, ky= multi_loadtxt(main_folder, ("kx.txt", "ky.txt"))
dirs = ['M103', 'M105', 'M133', 'M144']
step = 1
ticks = [-1, -0.5, 0, 0.5, 1]
i = 1
j = 0
for dir in dirs:
    theta = np.arcsin(np.sqrt(kx**2 + ky**2))
    phi = np.angle(kx+1j*ky)
    F2, F3 = multi_loadtxt(main_folder+'/'+dir, ("F_2.csv", "F_3.csv"))
    Cx = F2*np.cos(theta)*np.cos(phi) - F3*np.sin(phi)
    Cy = F2*np.cos(theta)*np.sin(phi) + F3*np.cos(phi)
    # kxi = np.linspace(kx.min(), kx.max(), kx.size)
    # kyi = np.linspace(ky.min(), ky.max(), ky.size)
    # Cx = interp2d(kx, ky, Cx)(kxi, kyi)
    # Cy = interp2d(kx, ky, Cy)(kxi, kyi)


    amp = np.sqrt(Cx**2 + Cy**2)
    Cx, Cy = Cx/amp, Cy/amp
    amp = binary_data_processing(F3, 0, -1, 1)
    # x = np.linspace(0, 1, np.size(kx))
    # y = np.linspace(0, 1, np.size(ky))
    # X, Y = np.meshgrid(x, y)
    # CX, CY = np.meshgrid(Cx, Cy)
    # print(np.size(Cx), axis=0)
    # amp_norm = amp/np.max(amp)
    # axs[i,j].quiver(kx[::step], ky[::step], Cx[::step], Cy[::step], amp[::step], cmap='rainbow',
    #               width=0.02, pivot="mid", headwidth=2, headaxislength=5, scale=30)
    # axs[i,j].streamplot(x, y, CX, CY, density=0.5)
    # axs[i,j].quiver(kx[::step], ky[::step], Cx[::step], Cy[::step], amp[::step], cmap='rainbow', width=0.01, scale=25)
    # axs[i,j].add_patch(circle1)
    axs[i,j].quiver(kx[::step], ky[::step], Cx[::step], Cy[::step], amp[::step], cmap='rainbow', width=0.004, pivot="mid", headwidth=3.0, headlength=3.0, headaxislength=2)
    # axs[i,j].quiver(kx[::step], ky[::step], Cx[::step], Cy[::step], amp[::step], cmap='rainbow', width=0.01, pivot="mid", headwidth=3, headlength=5, minshaft=1, minlength=2)


    # axs[i,j].set_xlabel(r'${k_x/k_0}$')
    axs[i,j].tick_params(which='both', direction='in')
    axs[i,j].set_xticks(ticks)
    axs[i,j].xaxis.set_ticklabels([])
    axs[i,j].yaxis.set_ticklabels([])
    j += 1
# axs[1,0].add_patch(circle1)
# axs[1,1].add_patch(circle1)
# axs[1,2].add_patch(circle1)
# axs[1,3].add_patch(circle1)
# axs[1,0].set_ylabel(r'${k_y/k_0}$')
# axs[1,0].set_yticks(ticks)
plt.show()
# plt.savefig('/home/ashalev/Desktop/4fig.png')
# plt.savefig('/home/ashalev/Desktop/4fig.eps', format='eps')
# plt.close()

# dir = 'M133spectra_output'
# x, y = multi_loadtxt(main_folder+dir, ("RLs-fano.txt", "Q0s-fano.txt"))
# plt.figure(figsize=(10,10))
# plt.xticks([2e-2, 8e-2])
#
# plt.plot(x, y, 'r', lw=0.5)
# plt.show()

