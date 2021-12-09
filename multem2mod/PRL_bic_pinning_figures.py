import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import LogNorm
#figure 3b
#------------------
#blank
#------------------
#figure 3d -------------------------------------------------------------------------------------------------------------
dir = 'kx_ky'
z_data = np.loadtxt(dir+'/'+'400_eps_220.0.txt', dtype=float, delimiter=',')
fig, axs = plt.subplots(1,1)
im = plt.imshow(z_data.T, extent = (-1, 1, 1, -1), cmap=cm.hot, norm=LogNorm(), aspect='auto', interpolation = 'none')
cb = plt.colorbar(im)
cb.set_label('reflectance')
plt.ylabel(r'${\omega d / 2\pi c }$')
plt.gca().invert_yaxis()
plt.show()
#-----------------------------------------------------------------------------------------------------------------------
