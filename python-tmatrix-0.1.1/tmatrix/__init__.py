import sys
if sys.version_info.major >= 3:
    from tmatrix.tmatrix import OBLATE, PROLATE
    from tmatrix.tmatrix import tmfixed, tmrandom, tmbisphere, tmnsfixed, tmnsrandom
    from tmatrix.tmatrix import FixedSpheroid, RandomSpheroid, FixedCylinder, RandomCylinder, FixedChebyshev, RandomChebyshev, RandomBiSphere, FixedSphereCluster, RandomSphereCluster, RandomNanorod
else:
    from tmatrix import OBLATE, PROLATE
    from tmatrix import tmfixed, tmrandom, tmbisphere, tmnsfixed, tmnsrandom
    from tmatrix import FixedSpheroid, RandomSpheroid, FixedCylinder, RandomCylinder, FixedChebyshev, RandomChebyshev, RandomBiSphere, FixedSphereCluster, RandomSphereCluster, RandomNanorod
