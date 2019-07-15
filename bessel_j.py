#!/usr/bin/env python3
from mpmath import *
from scipy.special import jv, jve
import numpy as np
mp.dps = 35; mp.pretty = True
v_arg2 = 50
z_arg = (1+1j)*500
for v_arg in range(15):
    res_mp = ((besselj(v_arg, z_arg)))
    res = jv(v_arg, z_arg)
    res2 = jve(v_arg, z_arg)/exp(-abs(z_arg.imag))
    print("v =",v_arg, "bess diff:",complex(res/res_mp-1), "\t",complex(res2/res_mp-1))
res = besselj(v_arg2, z_arg)
print(res)
print("   %20.17g  %20.17g"%(np.real(res), np.imag(res)))

