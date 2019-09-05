# wolfram_alpha
# exp(-(4+9i)^2) * (1.0 - erf(-i * (4+9i)))
# 0.05225352990488347758787552890909121671242855292674236050... +
# 0.02298883403846503952236217596768748188473022278897852839... i
#
# z_old   (5.22535299048835084E-002,2.29888340384650433E-002)
# z_lib  (-2.27750539897804546E-003,8.67894276156025435E-003)
#
# python (-0.002277505398978044+0.008678942761560251j)
from scipy.special import erf,erfc
import numpy as np

# print(erf(4+9j))

# z = 4+9j
# z = 250+250j
z = -13.905919999999979-13.752619999999947j # (-3.4862337365824946e-06-4.430872355966719e-05j)
cerf_py = np.exp(-z**2) * (1.0 - erf(-1j * z))
print(cerf_py)
cerf_py = np.exp(-z**2) * (erfc(-1j * z))
print(cerf_py)


from mpmath import *
mp.dps = 150; mp.pretty = True
cerf_py = np.exp(-z**2) * (1.0 - erf(-1j * z))
print(cerf_py)
cerf_py = np.exp(-z**2) * (erfc(-1j * z))
print(cerf_py)


