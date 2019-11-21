(SPH, CYL, GCH, NANOROD) = (-1, -2, -3, -9)
(OBLATE, PROLATE) = range(2)
(CGM, OOSM) = range(2)
MAXSPHERES = 150

def getErrorDescFixed(ErrCode):
    if ErrCode == 1:
        return 'Convergence is not obtained for the current value of NPN1'
    elif ErrCode == 2:
        return 'NGAUSS is greater than NPNG1'
    elif ErrCode == 3:
        return 'An angular parameter is outside its allowable range'
    elif ErrCode == 4:
        return 'NMAX is greater than NPNG1'
    else:
        return 'Unknown error code'

def getErrorDescRandom(ErrCode):
    if ErrCode == 1:
        return 'NK is greater than 1000'
    elif ErrCode == 2:
        return 'Convergence is not obtained for the current value of NPN1'
    elif ErrCode == 3:
        return 'NGAUSS is greater than NPNG1'
    elif ErrCode == 4:
        return 'NMAX1 is greater than NPNG1'
    elif ErrCode == 5:
        return 'NMAX is greater than NPNG1'
    elif ErrCode == 6:
        return 'Error in subroutine CCG'
    elif ErrCode == 7:
        return 'Error in subroutine CCGIN'
    else:
        return 'Unknown error code'

def getErrorDescBiSphere(ErrCode):
    if ErrCode == 1:
        return 'R12 is smaller than R1 + R2'
    elif ErrCode == 2:
        return 'R2 is greater than R1'
    elif ErrCode == 3:
        return 'NOD dimension exceeded'
    elif ErrCode == 4:
        return 'NOTD dimension exceeded'
    elif ErrCode == 5:
        return 'Error in subroutine CCG'
    elif ErrCode == 6:
        return 'Error in subroutine CCGIN'
    else:
        return 'Unknown error code'

def getErrorDescFixedSphereCluster(ErrCode):
#    if ErrCode == 1:
#        return 'Maximum number of iterations was exceeded'
#    else:
        return 'Unknown error code'

def getErrorDescRandomSphereCluster(ErrCode):
    if ErrCode == 1:
        return 'Maximum number of iterations was exceeded'
    else:
        return 'Unknown error code'


def tmfixed(lmbda, np, re, rat, eps,  m, epse=1, error = 1.0e-4, ndgs = 4, alpha = 0.0, beta = 90.0, theta0 = 0.0, theta = 0.0, phi0 = 0.0, Phi = 0.0):
    import sys
    if sys.version_info.major >= 3:
        from tmatrix._tmfixed import tmfixed as tmf
    else:
        from _tmfixed import tmfixed as tmf

    from scipy import zeros

    S = zeros((2, 2), dtype = complex)
    Qext = 0.0
    Qsca = 0.0
    errcode = 0
    (Qext, Qsca, S[0, 0], S[0, 1], S[1, 0], S[1, 1], errcode) = tmf(re, rat, lmbda, m.real, m.imag, eps, epse, np, error, ndgs, alpha, beta, theta0, theta, phi0, Phi)

    if errcode > 0:
        raise ValueError(getErrorDescFixed(errcode))

    Qabs = Qext - Qsca
    Albedo = Qsca/Qext

    return (Qext, Qsca, Qabs, Albedo), S

def tmrandom(Lambda, NP, Re, Rat, Eps, m, Error = 1.0e-4, NDGS = 4, epse = 1):
    from _tmrandom import tmrandom as tmr
    from scipy import pi

    Cext = 0.0
    Csca = 0.0
    g = 0.0
    errcode = 0
    epse = 1.
    (Cext, Csca, g, errcode) = tmr(Re, Rat, Lambda, m.real, m.imag, Eps, epse, NP, Error, NDGS)

    if errcode > 0:
        raise ValueError(getErrorDescRandom(errcode))

    Qext = Cext/(pi*Re*Re)
    Qsca = Csca/(pi*Re*Re)
    Qabs = Qext - Qsca
    Albedo = Qsca/Qext

    return (Qext, Qsca, Qabs, Albedo, g)

def tmbisphere(Lambda, R1, m1, R2, m2, d12):
    from _tmbisphere import tmbisphere as tmbs
    from scipy import pi

    Cext = 0.0
    Csca = 0.0
    g = 0.0
    errcode = 0
    (Cext, Csca, g, errcode) = tmbs(Lambda, R1, m1.real, m1.imag, R2, m2.real, m2.imag, R1 + d12 + R2)

    if errcode > 0:
        raise ValueError(getErrorDescBiSphere(errcode))

    Qext = Cext/(pi*R1*R1)
    Qsca = Csca/(pi*R1*R1)
    Qabs = Qext - Qsca
    Albedo = Qsca/Qext

    return (Qext, Qsca, Qabs, Albedo, g)

def tmnsfixed(Lambda, Ri, SPi, mi, xScale = 1.0, rScale = 1.0, mrScale = 1.0, miScale = 1.0, MaxIter = 100, Meth = CGM, Error = 1.0e-6, SphError = 1.0e-3, CtrError = 1.0e-9, Phi = 0.0, Theta = 0.0):
    from _tmnsfixed import tmnsfixed as tmnsf
    from scipy import array, pi, zeros, hstack, vstack

    if len(Ri) > MAXSPHERES:
        raise ValueError("The number of spheres (%i) is greater than the maximum allowed (%i)" % (len(Ri), MAXSPHERES))

    NPart = len(Ri)
    Qext = 0.0
    Qabs = 0.0
    Qsca = 0.0
    g = 0.0
    errcode = 0
    Qei = zeros(len(Ri), dtype = float)
    Qai = zeros(len(Ri), dtype = float)
    SA = zeros(4, dtype = complex)
    SM = zeros((4, 4), dtype = float)

    Xi = 2.0*pi*Ri/Lambda
    XPi = 2.0*pi*SPi/Lambda
    Xi = hstack((Xi, zeros(MAXSPHERES - len(Xi), dtype = float)))
    XPi = vstack((XPi, zeros((MAXSPHERES - XPi.shape[0], XPi.shape[1]), dtype = float)))
    mii = hstack((mi, zeros(MAXSPHERES - len(mi), dtype = complex)))

    (Qext, Qabs, Qsca, g, Qei, Qai, SA, SM, errcode) = tmnsf(NPart, Xi, XPi[:, 0], XPi[:, 1], XPi[:, 2], mii.real, mii.imag, xScale, rScale, mrScale, miScale, MaxIter, Meth, Error, SphError, CtrError, Phi, Theta)
    S = -1.0j*array([[SA[1], -SA[2]], [-SA[3], SA[0]]])

    if errcode > 0:
        raise ValueError(getErrorDescFixedSphereCluster(errcode))

    Albedo = Qsca/Qext

    return (Qext, Qsca, Qabs, Albedo, g), S, SM, Qei[:NPart], Qai[:NPart]

def tmnsrandom(Lambda, Ri, SPi, mi, xScale = 1.0, rScale = 1.0, mrScale = 1.0, miScale = 1.0, MaxIter = 100, Rlx = 0.0, Error = 1.0e-6, SphError = 1.0e-3, CtrError = 1.0e-9, Theta = 0.0):
    from _tmnsrandom import tmnsrandom as tmnsr
    from scipy import pi, zeros, hstack, vstack

    if len(Ri) > MAXSPHERES:
        raise ValueError("The number of spheres (%i) is greater than the maximum allowed (%i)" % (len(Ri), MAXSPHERES))

    NPart = len(Ri)
    Qext = 0.0
    Qabs = 0.0
    Qsca = 0.0
    g = 0.0
    errcode = 0
    Qei = zeros(len(Ri), dtype = float)
    Qai = zeros(len(Ri), dtype = float)
    S = zeros((2, 2), dtype = complex)

    Xi = 2.0*pi*Ri/Lambda
    XPi = 2.0*pi*SPi/Lambda
    Xi = hstack((Xi, zeros(MAXSPHERES - len(Xi), dtype = float)))
    XPi = vstack((XPi, zeros((MAXSPHERES - XPi.shape[0], XPi.shape[1]), dtype = float)))
    mii = hstack((mi, zeros(MAXSPHERES - len(mi), dtype = complex)))

    (Qext, Qabs, Qsca, g, Qei, Qai, S, errcode) = tmnsr(NPart, Xi, XPi[:, 0], XPi[:, 1], XPi[:, 2], mii.real, mii.imag, xScale, rScale, mrScale, miScale, MaxIter, Rlx, Error, SphError, CtrError, Theta)

    if errcode > 0:
        raise ValueError(getErrorDescRandomSphereCluster(errcode))

    Albedo = Qsca/Qext

    return (Qext, Qsca, Qabs, Albedo, g), S, Qei[:NPart], Qai[:NPart]


# T-Matrix codes for particles in fixed position
def FixedSpheroid(Lambda, Rv, Eps, m, Type = OBLATE, Error = 1.0e-4, NDGS = 4, Alpha = 0.0, Beta = 90.0, Theta = 0.0, Phi = 0.0):
    Theta0 = 0.0
    Phi0 = 0.0
    if Type == OBLATE:
        return tmfixed(Lambda, SPH, Rv, 1.0, Eps, m, Error, NDGS, Alpha, Beta, Theta0, Theta, Phi0, Phi)
    else:
        return tmfixed(Lambda, SPH, Rv, 1.0, 1.0/Eps, m, Error, NDGS, Alpha, Beta, Theta0, Theta, Phi0, Phi)

def FixedCylinder(Lambda, Rv, Eps, m, Type = OBLATE, Error = 1.0e-4, NDGS = 4, Alpha = 0.0, Beta = 90.0, Theta = 0.0, Phi = 0.0):
    Theta0 = 0.0
    Phi0 = 0.0
    if Type == OBLATE:
        return tmfixed(Lambda, CYL, Rv, 1.0, Eps, m, Error, NDGS, Alpha, Beta, Theta0, Theta, Phi0, Phi)
    else:
        return tmfixed(Lambda, CYL, Rv, 1.0, 1.0/Eps, m, Error, NDGS, Alpha, Beta, Theta0, Theta, Phi0, Phi)

def FixedChebyshev(Lambda, Rv, Eps, m, Degree, Type = OBLATE, Error = 1.0e-4, NDGS = 4, Alpha = 0.0, Beta = 90.0, Theta = 0.0, Phi = 0.0):
    Theta0 = 0.0
    Phi0 = 0.0
    if Type == OBLATE:
        return tmfixed(Lambda, Degree, Rv, 1.0, Eps, m, Error, NDGS, Alpha, Beta, Theta0, Theta, Phi0, Phi)
    else:
        return tmfixed(Lambda, Degree, Rv, 1.0, 1.0/Eps, m, Error, NDGS, Alpha, Beta, Theta0, Theta, Phi0, Phi)

def FixedSphereCluster(Lambda, Ri, SPi, mi, MaxIter = 100, Meth = CGM, Error = 1.0e-6, SphError = 1.0e-3, CtrError = 1.0e-9, Phi = 0.0, Theta = 0.0):
    xScale = 1.0
    rScale = 1.0
    mrScale = 1.0
    miScale = 1.0
    return tmnsfixed(Lambda, Ri, SPi, mi, xScale, rScale, mrScale, miScale, MaxIter, Meth, Error, SphError, CtrError, Phi, Theta)


# T-Matrix codes for particles in random orientations
def RandomSpheroid(Lambda, Rv, Eps, m, Type = OBLATE, Error = 1.0e-4, NDGS = 4):
    if Type == OBLATE:
        return tmrandom(Lambda, SPH, Rv, 1.0, Eps, m, Error, NDGS)
    else:
        return tmrandom(Lambda, SPH, Rv, 1.0, 1.0/Eps, m, Error, NDGS)

def RandomCylinder(Lambda, Rv, Eps, m, Type = OBLATE, Error = 1.0e-4, NDGS = 4):
    if Type == OBLATE:
        return tmrandom(Lambda, CYL, Rv, 1.0, Eps, m, Error, NDGS)
    else:
        return tmrandom(Lambda, CYL, Rv, 1.0, 1.0/Eps, m, Error, NDGS)


def RandomNanorod(Lambda, Rv, Eps, m, Type = OBLATE, Error = 1.0e-4, NDGS = 4, epse=1):
    return tmrandom(Lambda, NANOROD, Rv, 1.0, Eps, m, Error, NDGS, epse=epse)


def RandomChebyshev(Lambda, Rv, Eps, m, Degree, Type = OBLATE, Error = 1.0e-4, NDGS = 4):
    if Type == OBLATE:
        return tmrandom(Lambda, Degree, Rv, 1.0, Eps, m, Error, NDGS)
    else:
        return tmrandom(Lambda, Degree, Rv, 1.0, 1.0/Eps, m, Error, NDGS)

def RandomBiSphere(Lambda, R1, m1, R2, m2, d12 = 0):
    return tmbisphere(Lambda, R1, m1, R2, m2, d12)

def RandomSphereCluster(Lambda, Ri, SPi, mi, MaxIter = 100, Rlx = 0.0, Error = 1.0e-6, SphError = 1.0e-3, CtrError = 1.0e-9, Theta = 0.0):
    xScale = 1.0
    rScale = 1.0
    mrScale = 1.0
    miScale = 1.0
    return tmnsrandom(Lambda, Ri, SPi, mi, xScale, rScale, mrScale, miScale, MaxIter, Rlx, Error, SphError, CtrError, Theta)

