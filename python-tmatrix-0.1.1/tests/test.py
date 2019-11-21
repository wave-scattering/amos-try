from tmatrix import tmfixed, PhaseMatrix
from scipy import arccos

(S, walb, errcode) = tmfixed(re = 10.0, lmbda = arccos(-1.0)*2.0, m = 1.5 + 0.02j, eps = 0.5, np = -1, ddelt = 0.001, ndgs = 2, alpha = 145.0, beta = 52.0, theta0 = 56.0, theta = 65.0, phi0 = 114.0, phi = 128.0, rat = 0.1)

print S, walb, errcode
print PhaseMatrix(S) 

