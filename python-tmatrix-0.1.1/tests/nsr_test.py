from tmatrix import RandomMultiSphere
from scipy import pi, zeros, ones

Ri = ones(7, dtype = float)
SPi = zeros((7, 3), dtype = float)
SPi[1, 2] = SPi[2, 1] = SPi[3, 0] = 2.0
SPi[4, 2] = SPi[5, 1] = SPi[6, 0] = -2.0
mi = ones(7)*2.0 + 0.1j

print RandomMultiSphere(Lambda = 2.0*pi, Ri = Ri, SPi = SPi, mi = mi, xScale = 0.25)


