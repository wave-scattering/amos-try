import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import LogNorm
#figure 3b
#------------------
#blank
#------------------
#figure 3d -------------------------------------------------------------------------------------------------------------
# dir = 'kx_ky'
# z_data = np.loadtxt(dir+'/'+'400_eps_220.0.txt', dtype=float, delimiter=',')
# fig, axs = plt.subplots(1,1)
# im = plt.imshow(z_data.T, extent = (-1, 1, 1, -1), cmap=cm.hot, norm=LogNorm(), aspect='auto', interpolation = 'none')
# cb = plt.colorbar(im)
# cb.set_label('reflectance')
# plt.ylabel(r'${\omega d / 2\pi c }$')
# plt.gca().invert_yaxis()
# plt.show()
#-----------------------------------------------------------------------------------------------------------------------
#figure 3d new-------------------------------------------------------------------------------------------------------------
dir = 'theta_phi'
R = np.loadtxt(dir+'/'+'400_eps_220.0.txt', dtype=float, delimiter=',')
theta = np.loadtxt(dir+'/'+'x_y.txt', dtype=float)
phi = np.linspace(0, 360.0, 300)
theta = np.linspace(0, 90.0, 300)
print(R.shape)
a = 400e-9
k0 = 0.502*2*np.pi/a
kx = k0*np.sin(theta)*np.cos(phi)
x = kx*a/np.pi
print(np.max(x))
print(np.min(x))
ky = k0*np.sin(theta)*np.sin(phi)
y = ky*a/np.pi
print(np.max(y))
print(np.min(y))
fig, axs = plt.subplots(1,1)
im = plt.imshow(R, extent = (-1, 1, 1, -1), cmap=cm.hot, norm=LogNorm(), aspect='auto', interpolation = 'none')
cb = plt.colorbar(im)
cb.set_label('reflectance')
plt.gca().invert_yaxis()
plt.show()
#---------------