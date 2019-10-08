       program axspartcl1

* Warning in module AXSPARTCL in file axspartcl.f: Variables set but never used:
*   RMF set at line 257 file axspartcl.f
*    XS set at line 682 file axspartcl.f
*
*  If changing LMX, one has to adjust NPN1 in AMPLDR to the same value
*                   NPNG1 is only within the AMPLDR
c Extinction for a single homogeneous sphere of radius 300nm  (300.001/300)
c host dielectric constant= (1.00000000000000,0.000000000000000E+000)
c sphere diel. constant= (2.10250000000000,0.000000000000000E+000)
c for lambda_0=1000nm is sg_tot=1.21460719430552
c QSCA=   2.15779683977996
c QEXT=  -2.15779683977996
C         FAC=LAM**2/(2.d0*P**2*REV**2)     !=2/xs**2
C relates effciences to normalized cross sections. For the setting above:
C         FAC=0,56289546467965428579933035116515
c S11=0.32957D+03 + i*0.17171D+03
c S12=0.00000D+00 + i*0.00000D+00
c S21=0.00000D+00 + i*0.00000D+00
c S22=0.32957D+03 + i*0.17171D+03
c For unpolarized light, extinction cross section is
c [(2.140) and (2.159) of {MTL}]:
c
c         C_{ext}= \fr{2\pi}{k_1} \mb{Im } (S_{11}+_{22})
c
c Homogeneous particle
c Particle parameters:
c radius =   300.000333332963
c particle diel. constant= (2.10250000000000,0.000000000000000E+000)
c background dielectric constant ZEPS0= (1.00000000000000,0.000000000000000E+000)
C
C==========================
C      Room for improvement:
C            1) it would be more resonable to replace rev by the
C               size parameter k*rev
C            2)
C            3)
C            4)
C--------/---------/---------/---------/---------/---------/---------/--
C  This routines calculates the single particle scattering properties
C  (including coated particles)
C
C                    make -f mkaxspsc
C
C k_l length in units (2*PI/A=) PI:    xkl= 0.8660254037844386d0
C
C Outputs the total elastic scattering cross section TCS
C
C Parameters:
C
C Partial wave expansion is used, which is badly convergent
C for large size parameters $x> 100$. In numerical applications, the
C series is to be cut off after
C          LMAX (LMX parameter here) \approx x+4x^{1/3}+2$.
C In the case of LMX=50 this means that x <= 35.
C If one wants to observe ripples, the cutoff for a given x has to
C be even larger
C
C  ALPHA and BETA - Euler angles (in degrees) specifying the orientation
C            of the scattering particle relative to the laboratory reference
C            frame (Refs. 6 and 7).
C  THET0 - zenith angle of the incident beam in degrees
C  THET - zenith angle of the scattered beam in degrees    !with respect to
C  PHI0 - azimuth angle of the incident beam in degrees    !the laboratory frame!!!
C  PHI - azimuth angle of the scattered beam in degrees
C
C  ICHOICE=1 if NAG library is available, otherwise ICHOICE=2
C
C  NCHECK  -  .EQ.0  THEN  NGSS=2*NGAUSS, FACTOR=1D0
C             .EQ.1  THEN  NGSS = NGAUSS, FACTOR=2D0: theta=pi/2 is mirror
C                          symmetry plane as in the case of Chebysh. particle,
C                          ellipsoid, and cylinder
C
C  NAXSM   -  .EQ.0 : Gauss abscissas do not have +/- theta symmetry
C             .EQ.1 : Gauss abscissas have +/- theta symmetry
C
C  NDGS - controlling the number ND=NDGS*NMAX of division points in
C  computing integrals over the particle surface (Ref. 5).
C  For compact particles, the
C  recommended value is 2. For highly aspherical particles larger
C  values (3, 4,...) may be necessary to obtain convergence.
C  The code does not check convergence over this parameter.
C  Therefore, control comparisons of results obtained with
C  different NDGS-values are recommended.
C
C
C  Computation can be speed up if one sets YNCHECK=.FALSE..
C  However, then the check of Gauss integrations
C  convergence is not performed.
C---------------------------------------------------------------------
      implicit none
      integer LMX,LCS,ILCS,ikl,ieps,istep,ide,ndefp,itter
      integer NOUT,NOUTI,NSTEP,NFIN,NMAT,NP,NPP,NDGS,NDGSP
      real*8 TOL,DEFP,DEFPP,DDELT,DDELTP,x_max,x_min
      real*8 hlength_max,hlength_min,rl_min,rl_max
      real*8 lambda_min,lambda_max,omega_max
      COMPLEX*16 ZEPS0,CCEPS,CSEPS           !,ZARTAN
      character*1 ync,yncv
      logical ynperfcon,ynperfconv,ynintens,ynoptth,ynbrug,yncheck
cc      external ZARTAN

c Parameters:
C ::: number of the output unit for cross sections and scattering matrix
      PARAMETER (NOUT=35)
C ::: number of the output unit for the field intensity
      PARAMETER (NOUTI=60)
c Maximal number of spherical harmonics used. The floating number is
c specified below by the value of LMAX parameter
      PARAMETER (lmx=100)
*
* If convergence test in the calculation of the scattering cross sections
* is to be performed, yncheck=.true., otherwise yncheck=.false.
      parameter (yncheck=.false.)
*
* If particle is coated, ync=y, otherwise ync=n
      parameter (ync='n')
*
* ynperfcon=.true. if core is a perfect conductor, otherwise
* ynperfcon=.false.
*
      PARAMETER (ynperfcon=.false.)
*
* ynintens=.true. if the field intensity is to be calculated; otherwise
* ynintens=.false.
*
      PARAMETER (ynintens=.false.)
*
c number of coatings
      parameter (lcs=1)
c The coating layer to which material data are read in
      parameter (ilcs=1)
c if coated, the ratio 'core radius/particle radius'
c   (If lcs.ne.1, program is singular for rff=0. and 1.) - use homogeneous
c    particle instead!!!
c       PARAMETER (rff=0.95d0)
c background dielectric constant
      PARAMETER (ZEPS0=1.D0)
c
c Temporarily, before the coated part is finished,
c read in particle core dielectric function to CSEPS
c particle (core) dielectric constant (depending whether lcs=1 or lcs=2)
c
      PARAMETER (CCEPS=(80.0d0,0.0000000d0))
C      PARAMETER (CCEPS=(1.45D0,0.0d0)**2)    !SiO2
c      PARAMETER (CCEPS=(1.2D0,0.01d0)**2)    !to test lisac
c      PARAMETER (CCEPS=(-70.720839d0,7.05596d0))   !-70.720839,7.05596   Au for 1319nm
c      PARAMETER (CCEPS=(-10.84D0,0.762d0))   !JAP89_5774   ellipsoid for ld=633
c      PARAMETER (CCEPS=(-2.03D0,0.602d0))    !JAP89_5774   sphere for ld=354
C >>>     SPHERE (OUTER SHELL SCATTERER) PERMITTIVITY                  <<<
*  n(silica)=1.45  <--->    EPS(1)=2.1025D0
*  n(ZnS)=2.       <--->    EPS(1)=4.D0
      PARAMETER (CSEPS=(1.005d0,0.d0)**2)
c      PARAMETER (CSEPS=(1.05d0,0.d0)**2)    !to test lisac
c      PARAMETER (CSEPS=(-10.84D0,0.762d0))    !JAP89_5774

* material code number
c   NMAT=0             dispersionless dielectric
c   NMAT=1             Drude metal
c   NMAT=2             Ag
c   NMAT=3             Au
c   NMAT=4             ZnS
c   NMAT=5             Cu
c   NMAT=6             Al
c   NMAT=7             Pt
c   NMAT=8             Si
C  NMAT = 9           water
*
      PARAMETER(NMAT=0)
*
c Temporarily option for reading of the real data for the dielectric constant
c The number of the entries in a material data file to be read below
* Data files for Au,Cu,Al,Pt should be ordered with the decreased wavelength
* (omega increases in the loop and is oriented along the data file)
*
c          AGC.DAT                NFIN=73       ! from Palik
c          Audat.dat              NFIN=66       ! from Palik
c          Au_2dat.dat            NFIN=76       ! from JAW
c          Au*new.dat             NFIN=142
c          Cudat.dat              NFIN=47       ! from Palik
c          Aldat.dat              NFIN=80       ! from Palik
c          Ptdat.dat              NFIN=53       ! from Palik
c          Nidat.dat              NFIN=68       ! from Palik
c          Sidat.dat              NFIN=291
c          sieps.dat              NFIN=117
c          measured_Water_dispersion_T=24.txt       NFIN = 3600
*
      PARAMETER (NFIN=3600)
*
C ::: relative error allowed for the TCS. If the convergence
*     within TOL is not reached, program issues warning
      PARAMETER (TOL=1.d-3)
*
c If ynbrug=.true., performs Bruggeman approximation for ZEPS1. Otherwise
c ynbrug=false.
       parameter (ynbrug=.false.)
******************************************************************
c Declarations:

      integer ICHOICE,LMAX,NCHECK,NAXSM,NCHECKP,NAXSMP  !common block variables

      REAL*8 RMF(lcs),rff(lcs),RMUF,ff,filfrac   !ff for Bruggemann; filfrac for ZnS
      real*8 pi,xs,lambda,rap,rev_beg
      real*8 enw,xstep,revf,revin,revinl,revl
      real*8 omf(NFIN),omxf,reepsz,plasma,omxp
      real*8 delo,omega0,omega,rsnm,hlength
      real*8 RAT,RATP,AXI,REV,REVP,ALPHA,BETA,ALPHAE,BETAE      !common block variables
      real*8 THET0,THET,THETV,PHI0,PHI     !common block variables
      real*16 ceps1real(NFIN),  ceps1imag(nfin)


      complex*16 ceps1(NFIN),ZEPS1,ZEPS0V,Z1,Z2
      COMPLEX*16 ci,zeps(lcs+1)
      real*16 global_eps_r, global_eps_i
*
      COMMON /TOAMPLD/RAT,REV,ALPHA,BETA,DDELT
*
* transfers real*8 RAT,REV,ALPHA,BETA,DDELT  from the main to AMPLDR
*
      COMMON /TOTAMPLD/THET0,THET,PHI0,PHI
*
* transfers real*8 THET0,THET,PHI0,PHI from the main to AMPLDR
*
      COMMON /TOIAMPLD/NCHECK,NAXSM,NDGS
*
* transfers integers NCHECK,NAXSM,NDGS from the main
* to AMPLDR
*
      COMMON /TOITMT/ICHOICE,NPP,NCHECKP,NAXSMP,NDGSP

* transfers integers ICHOICE,NP,NCHECK,NAXSM,NDGS from the main to TMTAXSPV
*
      COMMON /TOTMT/DEFPP,RATP,REVP,ALPHAE,BETAE,DDELTP
*
* transfers real*8 DEFP,RAT,REV,ALPHAE,BETAE,DDELT from the main to TMTAXSPV
*
      COMMON /TOLTMT/ ynoptth
*
* transfers logical ynoptth from the main to TMTAXSP
*
      COMMON /DIELF/ zeps0v
        COMMON /REVF/ revf
*
* transfers ZEPS0,REV from the main to AMPL
*
      COMMON /TOSPHERECR/ rmuf
      COMMON /TOSPHERECH/ yncv
      COMMON /TOSPHERECL/ ynperfconv
      common /CYLPAR/ rsnm, hlength
      common /EPSGG/ global_eps_r, global_eps_i
*
* From here to spherec
*---------------------------------------------------------------
c Data:
      DATA PI/3.141592653589793d0/
      DATA ci/(0.d0,1.d0)/
*
* Convergence variable (has to be at least equal to 2):
*
      lmax=25  !lmx
*
*  FCC parameters:

c f=0.05  (--)
c       RMUF=0.2879411911484862d0
c f=0.06
c       RMUF=0.3059831741945870d0
c f=0.07
c       RMUF=0.3221166265075572d0
c f=0.08
c       RMUF=0.3367780601921260d0
c f=0.09
c       RMUF=0.3502632974822208d0
c f=0.1  (0f)
c       RMUF=0.3627831678597810d0
c f=0.12
c       RMUF=0.3855146420814098d0
c f=0.15  (--)
c       RMUF=0.4152830592077074d0
c f=0.2  (f)
c       RMUF=0.4570781497340833d0
c f=0.25  (--)
c       RMUF=0.4923725109213483d0
c f=0.3  (--)
c       RMUF=0.5232238679605294d0
c f=0.35  (--)
c       RMUF=0.5508116833525640d0
c f=0.36  (vATL)
c       RMUF=0.5560083268891277d0
c f=0.4  (ff)
c       RMUF=0.5758823822969722d0
c f=0.45  (--)
c       RMUF=0.5989418136982620d0
c f=0.5  (--)
c       RMUF=0.6203504908994001d0
c f=0.55  (--)
c       RMUF=0.6403754763690468d0
c f=0.58  (2w)
c       RMUF=0.6518131631368212d0
c f=0.6  (3f)
c       RMUF=0.6592207650508868d0
c f=0.62  (3w)
c       RMUF=0.6664655293794289d0
c f=0.64 (4f)
c       RMUF=0.6735561203842518d0
c f=0.65 (--)
c       RMUF=0.6770461107710686d0
c f=0.66 (4w)
c        RMUF=0.6805004874579641d0
c f=0.68 (5f)
c       RMUF=0.6873059437985093d0
c f=0.7  (6f)
c       RMUF=0.6939792343839249d0
c f=0.72 (7f)
c       RMUF=0.7005265949644416d0
c close packed :
c       RMUF=1.d0/DSQRT(2.D0)
      RMUF=1.d0
      rmf(lcs)=rmuf

CCCCCCCCCCCCCCCCCC    Assignement of common variables CCCCCCCCCCC
*
* ynoptth=.true. if you want to check optical theorem
      ynoptth=.false.
*
      ynperfconv=ynperfcon
      yncv=ync
        zeps0v=zeps0

* Preinitialization necessary since below rsnm,rev, and defp (all
* feeded in evanesc) are not always specified:
*
      rsnm=1.d0
      rev=1.d0
      defp=1.d0
*
      write(6,*)'Chose particle shape'
      write(6,*)'Only axially symmetric particle shapes allowed'
      write(6,*)'(The axis of axial symmetry is Z-axis)'
      write(6,*)
*
C--------/---------/---------/---------/---------/---------/---------/--
      write(6,*)'Chebyshev particles: type a positive number equal to
     1 the order of the Chebyshev polynomial=the number of wrinkles
     1 on the sphere surface'
      write(6,*)'oblate/prolate spheroids: type -1'
      write(6,*)'oblate/prolate cylinders: type -2'
      write(6,*)'generalized Chebyshev particle/droplet: type -3'
      write(6,*)'sphere cut by a plane on its top type -4'
      write(6,*)'sphere cut by a plane on its bottom: type -5'
      write(6,*)'upwardly oriented cone: type -6'
      write(6,*)'conus on a finite cylinder: type -7'
      write(6,*)'intensity around homogeneous/coated sphere: type -50'

      Open(unit = 90,file = 'epsWater.txt', status = 'unknown')
*
!     read(5,*) NP
      NP=-2
      write(6,*)'Auto-select: -2 (cylinder)'
cz      NP=-1                      !temporarily
*
      NPP=NP
C
C     NP
C
C     positive number     Chebyshev particles
C              r(\theta)=r_0[1+\eps*\cos(NP*\theta)]
C     -1                  oblate/prolate spheroids
C              r(\theta)=a\left[\sin^2\theta + (a^2/b^2)\cos^2\theta]^{-1/2}
C     -2                  oblate/prolate cylinders
C     -3                  generalized Chebyshev particle/droplet
C     -4                  sphere cut by a plane on its top
C     -5                  sphere cut by a plane on its bottom
C     -6                  cone
C     -7                  conus on a finite cylinder (in preparation)
*
C      PARAMETER (NP=-1)
*
* specify the shape of particles within a given NP class:
C     NP.gt.0 - DEFP = deformation parameter of a Chebyshev particle
C     NP=-1 - DEFP = the ratio of the horizontal to rotational axes. DEFP is
C             larger than 1 for oblate spheroids and smaller than 1 for
C             prolate spheroids.
C     NP=-2 - DEFP = the ratio of the cylinder diameter to its length.
C     NP=-3 - no DEFP is specified
C     NP=-4 - DEFP is the height (along the axial symmetry axis)
C             of the resulting cut sphere
C     NP=-5 - DEFP is the height (along the axial symmetry axis)
C             of the resulting cut sphere
C                Note that always DEFP.LT.2*REV specified
C     NP=-6 - DEFP is the height (along the axial symmetry axis) divided
C             by the width of a base
C     NP=-7 - DEFP  is the height (along the axial symmetry axis) divided
C             by the width of a base
C
C Warning:
C   In computations for spheres, use DEFP=1.000001 instead of DEFP=1.
C   DEFP=1 can cause overflows in some rare cases.
*
c      PARAMETER (DEFP=1.000001D0)
c      DEFP=1.000001
*
*
C--------/---------/---------/---------/---------/---------/---------/--
       Open(UNIT=10, FILE='test.txt')
*
      if (NP.gt.0) then
*
      write(6,*)'Radius of the undeformed Chebyshev particle in
     1 your units (in nm if dispersive data used)'
      read(5,*) rsnm
*
      write(6,*)'The amplitude of wrinkles on the sphere surface'
      read(5,*) defp
*
      write(6,*)'The equal-volume-sphere radius'
      read(5,*) rev
*
      RAT=1. D0
*
      else if (NP.eq.-50) then

      if(.not.ynintens) then
      write(6,*)'NP.eq.-50 option is only for intensity calculation!!!'
      stop
      end if

      write(6,*)'The (coated) sphere radius in your units'
      read(5,*) rsnm
      rev=rsnm
*
      else if (NP.eq.-1) then

      write(6,*)'The half-length of the spheroid along
     1   the ROTATIONAL AXIS z-axis in your units
     2  (in nm if dispersive data used)'
      read(5,*) hlength
cz      hlength=63.3d0

      write(6,*)'The half-length of the spheroid along the
     1  HORIZONTAL AXIS (in theta=pi/2 plane) in your units
     2 (in nm if dispersive data used)'
      read(5,*) rsnm
cz      rsnm=21.1d0
*
C     NP=-1 - DEFP = the ratio of the horizontal to rotational axes. DEFP is
C             larger than 1 for oblate spheroids and smaller than 1 for
C             prolate spheroids.
*
      DEFP=rsnm/hlength              !always revolution axis length
*                                    !in the denominator
*
      rev=rsnm/DEFP**(1.D0/3.D0)     !=equal-volume-sphere radius
                             !Room for improvement here - it would be
                                   !more resonable to replace rev
                               !by the size parameter k*rev
*
      if (lcs.gt.1) then

      write(6,*)'The ratio of the inner to outer spheroid half-length
     1  along the HORIZONTAL AXIS'
      read(5,*) revin

      if (revin.ge.1) then
           write(6,*)'The ratio cannot be greater or equal 1'
      stop
      end if

      revin=rev*revin
      end if               !lcs.gt.1

*
      RAT=1.D0
*
C--------/---------/---------/---------/---------/---------/---------/--
*
      else if (NP.eq.-2) then
*
      write(6,*)'Read cylinder maximal r/l'
!     read(5,*) rl_max
      rl_max = 0.75D0
      write(6,*)'Auto-set cylinder maximal r/l',rl_max

      write(6,*)'Read cylinder minimal r/l'
!     read(5,*) rl_min
      rl_min = 0.45D0
      write(6,*)'Auto-set cylinder minimal r/l', rl_min

      write(6,*)'Read amount of steps in length'
!     read(5,*) ndefp
      ndefp = 15
      write(6,*)'Auto-set amount of steps in length',ndefp

      write(6,*)'Read cylinder radius'
!     read(5,*) rsnm
      rsnm = 1.0D0
      write(6,*)'Auto-set cylinder radius',rsnm

       rsnm =  rsnm*2

* specify the shape:
C NP=-2 - DEFP = the ratio of the cylinder diameter to its length.
*
      hlength_max = rsnm/rl_min/2.d0
      hlength_min = rsnm/rl_max/2.d0
      DEFP=rsnm/hlength_max
      rsnm=rsnm/2.d0                                  !cylinder radius
      hlength_max=hlength_max/2.d0
      hlength_min=hlength_min/2.d0        !cylinder half-length
      rev=hlength_max*(3.D0*DEFP*DEFP/2.D0)**(1.D0/3.D0)  !=equal-volume-sphere radius

      RAT=1. D0
*
      else if (NP.eq.-3) then
C     NP=-3 - no DEFP is specified
      write(6,*)'The length of in your units'
      read(5,*) rsnm
      rev=rsnm
C      write(6,*)'rev(rsnm) not yet determined for NP=-3'
C      pause
*
      RAT=1. D0
*
      else if ((NP.eq.-4).or.(NP.eq.-5))  then
C     NP=-4,-5 - DEFP is the height (along the axial symmetry axis)
C                 of the resulting cut sphere
      write(6,*)'The radius of the original uncut sphere in your units'
C--------/---------/---------/---------/---------/---------/---------/--
      read(5,*) rsnm
      rev=rsnm
*
      write(6,*)'The height of the cut sphere in your units'
      read(5,*) defp
      defp=defp/rsnm
*
      RAT=1. D0
*
      else if (NP.eq.-6) then
*
C     NP=-6 -
C
      write(6,*)'The cone base diameter in your units'
      read(5,*) rsnm

      write(6,*)'The cone heigth of in your units'
      read(5,*) hlength
*
      rsnm=rsnm/2.d0
      rev=(hlength*rsnm**2/4.d0)**(1.d0/3.d0)   !=equal-volume-sphere radius
      RAT=1. D0
*
      else if (NP.eq.-7) then
*
C     NP=-7 -
C
      write(6,*)'Not ready yet!'
      pause
cc      write(6,*)'The length of in your units'
cc      read(5,*) rsnm
*
      end if                         ! end NP if
****************************************
*
      defpp=defp
*
*
      if (RAT.eq.1.) then
      write(6,*)'Particle size specified in terms of
     1    the equal-volume-sphere radius'
      else if (RAT.ne.1.) then
      write(6,*)'Particle size specified in terms of
     1    the equal-surface-area-sphere radius'
      end if
*
* equivalent-(volume/surface-area)-sphere radius
*
cc      write(6,*)'Read equal-volume-sphere radius in nm'
cc      read(5,*) rev
*
cc      rev=300.d0                         !feeded as REV to RSP* routines
*
      AXI=rev
      revf=rev
*
C  Equivalent equal-(volume/surface-area)-sphere radius
*
cc      REV=RAT*AXI                      !feeded as REV to RSP* routines
*

C  NCHECK  -  .EQ.0  THEN  NGSS=2*NGAUSS, FACTOR=1D0
C             .EQ.1  THEN  NGSS = NGAUSS, FACTOR=2D0: theta=pi/2 is mirror
C                          symmetry plane as in the case of Chebysh. particle,
C                          ellipsoid, and cylinder
*
      NCHECK=0
*
      IF (NP.EQ.-1.OR.NP.EQ.-2) NCHECK=1         !ellipsoid(sphere) and cylinder
      IF (NP.GT.0.AND.(-1)**NP.EQ.1) NCHECK=1    !Chebysh. particle
*
C If theta=pi/2 is not a scatterer mirror symmetry plane:
C  NAXSM   -  .EQ.0 : Gauss abscissas do not have +/- theta symmetry
C             .EQ.1 : Gauss abscissas have +/- theta symmetry

      NAXSM=1

      IF (NP.LE.-4) NAXSM=0

C--------/---------/---------/---------/---------/---------/---------/--
C
C  ALPHA and BETA - Euler angles (in degrees) specifying the orientation
C    of the scattering particle relative to the laboratory reference
C    frame (Refs. 6 and 7).
*
C       ALPHA=90.D0           !laboratory frame coincides with particle frame
C      BETA=0.D0

      Write(6,*)'ALPHA and BETA - Euler angles (in degrees) specifying
     &    the orientation  of the scattering particle relative to the
     & laboratory reference  frame'
!     READ(5,*)  ALPHA, BETA
      ALPHA = 0.D0
      BETA = 0.D0
      write(6,*)'Auto-set ALPHA and BETA',ALPHA,BETA



      if((ALPHA.eq.0).and.(BETA.eq.0))
     & write(6,*)'Laboratory frame coincides with particle frame'

* DDELT - the desired absolute accuracy of computing the
* expansion coefficients of a normalized scattering matrix.
* (This accuracy is usually worse by a factor of 10 than
* the accuracy of computing the optical cross sections.)
* Since convergence test is only performed for the accuracy
* of computing the optical cross sections, DDELT is reset
* later on to DDELT=0.1D0*DDELT
*
      DDELT=TOL
*
C      THET0 - zenith angle of the incident beam in degrees
C      THET - zenith angle of the scattered beam in degrees
C      PHI0 - azimuth angle of the incident beam in degrees
C      PHI - azimuth angle of the scattered beam in degrees
*
      if (.not.ynintens) then
      write(6,*)'Specify (theta,phi) angles of the incident beam
     1             (in degrees)'
!     read(5,*) THET0,PHI0
      THET0 =90.D0
      PHI0 = 0.D0
      write(6,*)'Auto-set...'
*
      write(6,*)'Specify (theta,phi) angles of the scattered beam
     1             (in degrees)'
!     read(5,*) THET,PHI
      THET = 90.D0
      PHI = 0.D0
      write(6,*)'Auto-set...'
*
      else if (ynintens) then
      write(6,*)'PHI angle of incidence of the incident plane wave
     &  is for intensity calculations set to zero'
      write(6,*)'Specify Theta angle of incidence (in degrees)'
      read(5,*) THET0
cz      THET0=90.d0
      PHI0=0.d0
      THETV=THET0
      end if
*
C--------/---------/---------/---------/---------/---------/---------/--
cc      THET0=56.D0
cc      THET=65.D0
cc      PHI0=114.D0
cc      PHI=128.D0
*
* test setup:

      if ((thet0.gt.180.).or.(thet.gt.180.)) then
      write(6,*)'Theta angles has to be smaller than 180'
      stop
      end if

      if ((thet0.lt.0.).or.(thet.lt.0.)) then
      write(6,*)'Theta angles has to be positive'
      stop
      end if

      if ((phi0.gt.360.).or.(phi.gt.360.)) then
      write(6,*)'Phi angles has to be smaller than 360'
      stop
      end if

        IF ((NP.EQ.-4).or.(NP.EQ.-5)) THEN
      if ((defp.le.0).or.(defp.ge.2)) then
      write(6,*)'DEFP has to be >0 and <2'
      stop
      end if
      END IF

C--------/---------/---------/---------/---------/---------/---------/--
*
* If NAG library is available, set ICHOICE=1, otherwise ICHOICE=2

      ICHOICE=2

*  controlling the number ND=NDGS*NMAX of division points in
C  computing integrals over the particle surface (Ref. 5).
C  For compact particles, the
C  recommended value is 2. For highly aspherical particles larger
C  values (3, 4,...) may be necessary to obtain convergence.
C  The code does not check convergence over this parameter.
C  Therefore, control comparisons of results obtained with
C  different NDGS-values are recommended.
*
*  Check that NDGS*LAMXD does not exceed NPNG1 value in subroutines
*  For a current values of LMAXD=50 and NPNG1=800 then NDGS<= 16!!!

      IF ((NP.EQ.-4).or.(NP.EQ.-5)) THEN
         NDGS=16
      ELSE IF (NP.EQ.-6) THEN
         NDGS=16
      ELSE IF ((NP.EQ.-1).and.(max(defp,1.d0/defp).gt.1.5d0)) THEN !spheroids
         if (ynintens) then
            NDGS=min(40.d0,14*max(defp,1.d0/defp))
         else
            NDGS=min(16.d0,4*max(defp,1.d0/defp))
         end if
      ELSE
         NDGS=4
      END IF
*
      WRITE(6,*) 'NDGS=',NDGS
*
      IF (YNINTENS) THEN

          RATP=RAT
          NDGSP=NDGS
          NAXSMP=NAXSM
          NCHECKP=NCHECK
          revp=rev
          DDELTP=DDELT
          ALPHAE=ALPHA
          BETAE=BETA

      END IF
*
      IF (ICHOICE.EQ.1) THEN
      WRITE(6,*) 'NAG ROUTINES USED FOR THE MATRIX INVERSION'
      ELSE IF (ICHOICE.EQ.0) THEN
      WRITE(6,*) 'NAG ROUTINES (FOR THE MATRIX INVERSION) ARE NOT USED'
      END IF
      WRITE(6,*)
*
      IF (NCHECK.EQ.0) THEN
      WRITE(6,*) 'Particle without theta=pi/2 mirror symmetry'
      ELSE IF (NCHECK.EQ.1) THEN
      WRITE(6,*) 'Particle has  theta=pi/2 mirror symmetry'
      END IF
      WRITE(6,*)
*
      IF (NAXSM.EQ.0) THEN
      WRITE(6,*) 'Gauss abscissas not +/- theta symmetric'
      ELSE IF (NAXSM.EQ.1) THEN
      WRITE(6,*) 'Gauss abscissas +/- theta symmetric'
      END IF
      WRITE(6,*)
*
      WRITE(NOUT,5454) ICHOICE,NCHECK
 5454 FORMAT ('ICHOICE=',I1,'  NCHECK=',I1)
      WRITE(NOUT,*)'NAXSM=', NAXSM

      IF(NP.EQ.-1.AND.DEFP.GE.1D0) PRINT 7000,DEFP
      IF(NP.EQ.-1.AND.DEFP.LT.1D0) PRINT 7001,DEFP
      IF(NP.GE.0) PRINT 7100,NP,DEFP
      IF(NP.EQ.-2.AND.DEFP.GE.1D0) PRINT 7150,DEFP
      IF(NP.EQ.-2.AND.DEFP.LT.1D0) PRINT 7151,DEFP
      IF(NP.EQ.-3) PRINT 7160
      IF(NP.EQ.-4) PRINT 7170,DEFP
      PRINT 7200,DDELT
      IF (DABS(RAT-1D0).LE.1D-6) PRINT 8003, AXI
      IF (DABS(RAT-1D0).GT.1D-6) PRINT 8004, AXI

 7000 FORMAT('OBLATE SPHEROIDS, A/B=',F11.7)
 7001 FORMAT('PROLATE SPHEROIDS, A/B=',F11.7)
 7100 FORMAT('CHEBYSHEV PARTICLES, T',
     &       I1,'(',F5.2,')')
 7150 FORMAT('OBLATE CYLINDERS, D/L=',F11.7)
 7151 FORMAT('PROLATE CYLINDERS, D/L=',F11.7)
 7160 FORMAT('GENERALIZED CHEBYSHEV PARTICLES')
 7170 FORMAT('SHERE CUT BY A PLANE, DEFP=H=',F11.7)
 7200 FORMAT ('ACCURACY OF COMPUTATIONS DDELT = ',D8.2)
 8003 FORMAT('EQUAL-VOLUME-SPHERE RADIUS=',F8.4)
 8004 FORMAT('EQUAL-SURFACE-AREA-SPHERE RADIUS=',F8.4)
C--------/---------/---------/---------/---------/---------/---------/--
* Checking set up:

*
      if((np.gt.0).and.(dabs(defp).ge.1.d0)) then
      write(6,*)'Absolute value of defp has to be less than 1.!!!'
      stop
      end if
*
      if((np.eq.-4).and.(defp.ge.2.d0*REV)) then
      WRITE(6,*)'Invalid parameters for a cut sphere!'
      WRITE(6,*)'Execution stopped!'
      write(6,*)'The defp has to be less than 2*(sphere radius) !!!'
      stop
      end if
*
      if((np.eq.-1).and.(defp.eq.1)) then
       write(6,*)'Use DEFP=1.000001 instead of DEFP=1'
      end if
*
      if (nmat.gt.1) then
        write(6,*)'Real material data are to be provided'
      if (ynbrug) write(6,*)'Bruggeman approx. used!'
      if (ynbrug) write(nout,*)'#Bruggeman approximation performed'
      end if
*
      if ((ync.eq.'y'.and.lcs.eq.1).or.(ync.eq.'n'.and.lcs.ne.1)) then
      write(6,*)'Check compatibility of YNC and LCS'
      stop
      end if

C--------------------------------------------------------------------
* Reading in the input data:
c      write(6,*)'Read the particle (core) dielectric constant'
c      read(5,*) zeps(1)

*  n(silica)=1.45  <--->    ZEPS(1)=2.1025D0
*  n(ZnS)=2.       <--->    ZEPS(1)=4.D0
      ZEPS(1)=cceps
      if(lcs.gt.1) zeps(lcs)=cseps
      zeps(lcs+1)=zeps0

******************************
*
      if (lcs.ge.2) then                  !coated particle
*
cc      write(6,*)'Read equal-volume-core radius in nm'
cc      read(5,*) rff(1)
cc      rff(1)=204.9d0
cc      rff(1)=rff(1)/rs
*
      if (np.eq.-50) then
*
      write(6,*)'Coated sphere core radii r(l) labelled from 1 for the
     & inner core till LCS for the outer shell
     &   ===> r(lcs) is the sphere radius'
C--------/---------/---------/---------/---------/---------/---------/--
      do ikl=1,lcs-1
*
      write(6,*)'Read in r(l) for l=',ikl
      read (5,*) rff(ikl)
      rff(ikl)=rff(ikl)/rsnm

c      rff(1)=0.75d0
      rmf(ikl)=rff(ikl)*rmuf

      if ((lcs.gt.2).and.(ikl.ge.2)) then
      write(6,*)'Read in the lth-sphere layer diel. const. for l=',ikl
      write(6,*)'(In case of dispersive component, give 1. )'
      read (5,*) zeps(ikl)
*
      end if         !lcs.ikl

      enddo          !ikl
*
      end if         !np.eq.-50
*
      end if         !lcs.ge.2
*
*oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
* Scanning over frequency interval:
*-------------------------------------------------
*
      if (ynintens) then
      write(6,*)'ELECTRIC FIELD INTENSITY CALCULATION'
      write(6,*)'Read lambda (in the vacuum) in nm'
      else if (.not.ynintens) then
      write(6,*)'READ INITIAL (MINIMAL) x-parameter'
      end if
C--------/---------/---------/---------/---------/---------/---------/--

!     read(5,*) x_min
      x_min = 0.41D0
      write(6,*)'Auto-set x_min', x_min

      lambda = 2 *pi * rsnm/x_min
cv      lambda=633.d0                      !JAP89_5776 for gold ellipsod
c      lambda=354.d0                      !JAP89_5776 for silver sphere
cc      lambda=500d0

*
* size parameter is customarily defined as the ratio of
* circumference of particle to the wavelength in the host medium
* in which the particle is embedded
*                    x=kr=\sg a=2*pi*r/lambda
*      xs=2.d0*pi*rs*dble(sqrt(zeps0))/lambda
* convert lambda to the lambda in vacuum:
c      lambda=lambda*dble(sqrt(zeps0))
c  omega=2.d0*pi*rsnm/(lambda*rmuf)=xs/(rmuf*dble(sqrt(zeps0))),
c where rs is the particle radius (in nm) and  lambda
c is the wavelengths (in nm)
c in the vacuum:

         omega=2.d0*pi*rev/(lambda*rmuf)
         omega0=omega

c      write(6,*)'Read omega'
c      read(5,*) omega
c      WRITE(6,*)'omega=', omega
c      if (1.eq.1) nstep=0             ! temporarily
c      delo=0.d0
c      omega0=omega
c      go to 11
*
* Option for omega input:
c      write(6,*)'Read omega ='
c      read(5,*) omega
*
c      xs=RMUF*omega*dble(sqrt(zeps0))
*
c       write(6,*)'Equiv. size parameter x=2*pi*rs*n_0/lambda=',xs
*
      if (.not.ynintens) then
      write(6,*)'Scan up to x-parameter (in nm)'
!     read(5,*) x_max
      x_max = 0.55D0
      write(6,*)'Auto-set x_max', x_max

      enw = 2*pi*rsnm/x_max
c      enw=500
      write(6,*)'Amount of scanning steps'
!     read(5,*) nstep
      nstep = 20
      write(6,*)'Auto-set nstep',nstep

c      xstep=5
      end if
*
C ::: number of steps on frequency interval:
      if ((lambda.eq.enw).or.(ynintens)) then
       XSTEP=0
       delo=0.d0
      else
       XSTEP=(lambda-enw)/nstep
C ::: width of the searched frequency interval
       ENW=2.d0*pi*rev/(enw*rmuf)
       enw=enw-omega0
       delo=enw/dble(nstep-1)
      end if
*
C--------/---------/---------/---------/---------/---------/---------/--
*                  --------------------------------
* output initial statements

      OPEN(UNIT=NOUT,FILE='axs-scs.dat')
      rewind(NOUT)
      WRITE(NOUT,*)'#Orientationally averaged scattering cs for a single
     & particle'
      WRITE(NOUT,*)'#(cross sections normalized per
     & equal-volume-sphere surface S=pi*rev**2)'
      write(nout,*)
      write(nout,*)'#Equiv.-volume-sphere radius in nm=', rev
      write(nout,*)
      write(NOUT,*)'#Angular momentum cut-off LMAX=',LMAX
        write(NOUT,*)'#In the number ND=NDGS*LMAX of GIP points,
     & NDGS=',NDGS
        write(NOUT,*)'#DEFP=',DEFP
      write(nout,*)'#Material number =',NMAT
        if (ynbrug) write(nout,*)'#Bruggeman approximation performed'
      if (ync .eq.'n') then
        write(nout,*)'#Homogeneous particle'
      else if (ync.eq.'y') then
       write(nout,*)'#coated particle'
      end if
      write(nout,*)
      write(nout,*)'#host dielectric constant=', zeps(lcs+1)
      if (nmat.ge.1) write(nout,*)'#Dispersive layer number=', ilcs
      if ((lcs.eq.1).and.(nmat.eq.0))
     &     write(nout,*)'#sphere diel. const.=', zeps(1)
      if (lcs.ge.2) then
      write(nout,*)'#core radius/sphere radius =',rff(1)
      if ((ilcs.ne.1).or.(nmat.eq.0))
     &  write(nout,*)'#sphere core diel. const.=', zeps(1)
      if ((ilcs.ne.lcs).or.(nmat.eq.0))
     &  write(nout,*)'#coating diel. const.=',zeps(lcs)
C--------/---------/---------/---------/---------/---------/---------/--
      end if
      WRITE(NOUT,*)'#Vacuum lambda_0 and sg_sc in columns'
      write(nout,*)

      OPEN(UNIT=NOUT+1,FILE='axs-oavext.dat')
      rewind(NOUT+1)
      WRITE(NOUT+1,*)'#Orientationally averaged extinction cs for a
     & single (coated) particle'
      WRITE(NOUT+1,*)'#(cross sections normalized per
     & equal-volume-sphere surface S=pi*rev**2)'
      write(nout+1,*)
      write(nout+1,*)'#Equiv.-volume-sphere radius in nm=', rev
      write(nout+1,*)
        write(NOUT+1,*)'#Angular momentum cut-off LMAX=',LMAX
        write(NOUT+1,*)'#In the number ND=NDGS*LMAX of GIP points,
     & NDGS=',NDGS
        write(NOUT+1,*)'#DEFP=',DEFP
      write(nout+1,*)'#Material number =',NMAT
        if (ynbrug) write(nout+1,*)'#Bruggeman approximation performed'
      if (ync .eq.'n') then
        write(nout+1,*)'#Homogeneous particle'
      else if (ync.eq.'y') then
       write(nout+1,*)'#coated particle'
      end if
      write(nout+1,*)
      write(nout+1,*)'#host dielectric constant=', zeps(lcs+1)
      if (nmat.ge.1) write(nout+1,*)'#Dispersive layer number=', ilcs
      if ((lcs.eq.1).and.(nmat.eq.0))
     &     write(nout+1,*)'#sphere diel. const.=', zeps(1)
      if (lcs.ge.2) then
      write(nout+1,*)'#core radius/sphere radius =',rff(1)
      if ((ilcs.ne.1).or.(nmat.eq.0))
     &  write(nout+1,*)'#sphere core diel. const.=', zeps(1)
      if ((ilcs.ne.lcs).or.(nmat.eq.0))
     &  write(nout+1,*)'#coating diel. const.=',zeps(lcs)
C--------/---------/---------/---------/---------/---------/---------/--
      end if
      WRITE(NOUT+1,*)'#Vacuum lambda_0 and sg_tot in columns'
      write(nout+1,*)


      OPEN(UNIT=NOUT+2,FILE='axs-abs.dat')
      rewind(NOUT+2)
      WRITE(NOUT+2,*)'#Orientationally averaged absorption cs for a
     & single (coated) particle'
      WRITE(NOUT+2,*)'#(cross sections normalized per
     & equal-volume-sphere surface S=pi*rev**2)'
      write(nout+2,*)
      write(nout+2,*)'#Equiv.-volume-sphere radius in nm=', rev
      write(nout+2,*)
        write(NOUT+2,*)'#Angular momentum cut-off LMAX=',LMAX
        write(NOUT+2,*)'#In the number ND=NDGS*LMAX of GIP points,
     & NDGS=',NDGS
        write(NOUT+2,*)'#DEFP=',DEFP
      write(nout+2,*)'#Material number =',NMAT
        if (ynbrug) write(nout+2,*)'#Bruggeman approximation performed'

      if (ync .eq.'n') then
        write(nout+2,*)'#Homogeneous particle'
      else if (ync.eq.'y') then
       write(nout+2,*)'#coated particle'
      end if
      write(nout+2,*)
      write(nout+2,*)'# host dielectric constant=', zeps(lcs+1)
      if (nmat.ge.1) write(nout+2,*)'#Dispersive layer number=', ilcs
      if ((lcs.eq.1).and.(nmat.eq.0))
     &     write(nout+2,*)'#sphere diel. const.=', zeps(1)
      if (lcs.ge.2) then
      write(nout+2,*)'#core radius/sphere radius =',rff(1)
      if ((ilcs.ne.1).or.(nmat.eq.0))
     &  write(nout+2,*)'#sphere core diel. const.=', zeps(1)
      if ((ilcs.ne.lcs).or.(nmat.eq.0))
     &  write(nout+2,*)'#coating diel. const.=',zeps(lcs)
C--------/---------/---------/---------/---------/---------/---------/--
      end if
      WRITE(NOUT+2,*)'#Vacuum lambda_0 and sg_abs in columns'
      write(nout+2,*)

C--------/---------/---------/---------/---------/---------/---------/--
      OPEN(UNIT=NOUT+3,FILE='axs-ext.dat')
      rewind(NOUT+3)
      WRITE(NOUT+3,*)'#Extinction cs for a
     & single (coated) particle in a fixed orientation'
      WRITE(NOUT+3,*)'#(cross sections normalized per
     & equal-volume-sphere surface S=pi*rev**2)'
      write(NOUT+3,*)
      write(NOUT+3,*)'#Equiv.-volume-sphere radius in nm=', rev
      write(NOUT+3,*)
        write(NOUT+3,*)'#Angular momentum cut-off LMAX=',LMAX
        write(NOUT+3,*)'#In the number ND=NDGS*LMAX of GIP points,
     & NDGS=',NDGS
        write(NOUT+3,*)'#DEFP=',DEFP
      write(NOUT+3,*)'#Material number =',NMAT
        if (ynbrug) write(NOUT+3,*)'#Bruggeman approximation performed'
      if (ync .eq.'n') then
        write(NOUT+3,*)'#Homogeneous particle'
      else if (ync.eq.'y') then
       write(NOUT+3,*)'#coated particle'
      end if
      write(NOUT+3,*)
      write(NOUT+3,*)'#host dielectric constant=', zeps(lcs+1)
      if (nmat.ge.1) write(NOUT+3,*)'#Dispersive layer number=', ilcs
      if ((lcs.eq.1).and.(nmat.eq.0))
     &     write(NOUT+3,*)'#sphere diel. const.=', zeps(1)
      if (lcs.ge.2) then
      write(NOUT+3,*)'#core radius/sphere radius =',rff(1)
      if ((ilcs.ne.1).or.(nmat.eq.0))
     &  write(NOUT+3,*)'#sphere core diel. const.=', zeps(1)
      if ((ilcs.ne.lcs).or.(nmat.eq.0))
     &  write(NOUT+3,*)'#coating diel. const.=',zeps(lcs)
C--------/---------/---------/---------/---------/---------/---------/--
      end if
      WRITE(NOUT+3,*)'#In columns: vacuum lambda_0 and sg_tot
     & (normalized and unnormalized)'
      write(NOUT+3,*)



      OPEN(UNIT=NOUT+5,FILE='axs-albedo.dat')
      rewind(NOUT+5)
      WRITE(NOUT+5,*)'#Orientationally averaged albedo for a single
     & (coated) particle'
      write(nout+5,*)'#Equiv.-volume-sphere radius in nm=', rev
      write(nout+5,*)
        write(NOUT+5,*)'#Angular momentum cut-off LMAX=',LMAX
        write(NOUT+5,*)'#In the number ND=NDGS*LMAX of GIP points,
     & NDGS=',NDGS
        write(NOUT+5,*)'#DEFP=',DEFP
      write(nout+5,*)'#Material number =',NMAT
        if (ynbrug) write(nout+5,*)'#Bruggeman approximation performed'
      if (ync .eq.'n') then
        write(nout+5,*)'#Homogeneous particle'
      else if (ync.eq.'y') then
       write(nout+5,*)'#coated particle'
      end if
      write(nout+5,*)
      write(nout+5,*)'#host dielectric constant=', zeps(lcs+1)
      if (nmat.ge.1) write(nout+5,*)'#Dispersive layer number=', ilcs
      if ((lcs.eq.1).and.(nmat.eq.0))
     &     write(nout+5,*)'#sphere diel. const.=', zeps(1)
      if (lcs.ge.2) then
      write(nout+5,*)'#core radius/sphere radius =',rff(1)
      if ((ilcs.ne.1).or.(nmat.eq.0))
     &  write(nout+5,*)'#sphere core diel. const.=', zeps(1)
      if ((ilcs.ne.lcs).or.(nmat.eq.0))
     &  write(nout+5,*)'#coating diel. const.=',zeps(lcs)
C--------/---------/---------/---------/---------/---------/---------/--
      end if
        WRITE(NOUT+5,*)'#In columns:'
      WRITE(NOUT+5,*)'#Vacuum lambda_0, albedo, and tcs-(acs+tsc)'
      write(nout+5,*)

      OPEN(UNIT=NOUT+6,FILE='axs-dipolext.dat')
      rewind(NOUT+6)
      WRITE(NOUT+6,*)'#Orientationally averaged dipole ext. for a
     & single (coated) particle'
      write(nout+6,*)'#Equiv.-volume-sphere radius in nm=', rev
      write(nout+6,*)
        write(NOUT+6,*)'#Angular momentum cut-off LMAX=',LMAX
        write(NOUT+6,*)'#In the number ND=NDGS*LMAX of GIP points,
     & NDGS=',NDGS
        write(NOUT+6,*)'#DEFP=',DEFP
      write(nout+6,*)'#Material number =',NMAT
      if (ync .eq.'n') then
        write(nout+6,*)'#Homogeneous particle'
      else if (ync.eq.'y') then
       write(nout+6,*)'#coated particle'
      end if
      write(nout+6,*)
      write(nout+6,*)'#host dielectric constant=', zeps(lcs+1)
      if (nmat.ge.1) write(nout+6,*)'#Dispersive layer number=', ilcs
      if ((lcs.eq.1).and.(nmat.eq.0))
     &     write(nout+6,*)'#sphere diel. const.=', zeps(1)
      if (lcs.ge.2) then
      write(nout+6,*)'#core radius/sphere radius =',rff(1)
      if ((ilcs.ne.1).or.(nmat.eq.0))
     &  write(nout+6,*)'#sphere core diel. const.=', zeps(1)
      if ((ilcs.ne.lcs).or.(nmat.eq.0))
     &  write(nout+6,*)'#coating diel. const.=',zeps(lcs)
C--------/---------/---------/---------/---------/---------/---------/--
      end if
      WRITE(NOUT+6,*)'#Vacuum lambda_0 and dipole extinction in columns'
      write(nout+6,*)
C--------/---------/---------/---------/---------/---------/---------/--
      OPEN(UNIT=NOUT+7,FILE='axs-quadrext.dat')
      rewind(NOUT+7)
      WRITE(NOUT+7,*)'#Orientationally averaged quadrupole ext. for a
     & single (coated) particle'
      write(nout+7,*)'#Equiv.-volume-sphere radius in nm=', rev
      write(nout+7,*)
        write(NOUT+7,*)'#Angular momentum cut-off LMAX=',LMAX
        write(NOUT+7,*)'#In the number ND=NDGS*LMAX of GIP points,
     & NDGS=',NDGS
        write(NOUT+7,*)'#DEFP=',DEFP
      write(nout+7,*)'#Material number =',NMAT
        if (ynbrug) write(nout+7,*)'#Bruggeman approximation performed'
      if (ync .eq.'n') then
        write(nout+7,*)'#Homogeneous particle'
      else if (ync.eq.'y') then
       write(nout+7,*)'#coated particle'
      end if
      write(nout+7,*)
      write(nout+7,*)'#host dielectric constant=', zeps(lcs+1)
      if (nmat.ge.1) write(nout+7,*)'#Dispersive layer number=', ilcs
      if ((lcs.eq.1).and.(nmat.eq.0))
     &     write(nout+7,*)'#sphere diel. const.=', zeps(1)
      if (lcs.ge.2) then
      write(nout+7,*)'#core radius/sphere radius =',rff(1)
      if ((ilcs.ne.1).or.(nmat.eq.0))
     &  write(nout+7,*)'#sphere core diel. const.=', zeps(1)
      if ((ilcs.ne.lcs).or.(nmat.eq.0))
     &  write(nout+7,*)'#coating diel. const.=',zeps(lcs)
C--------/---------/---------/---------/---------/---------/---------/--
      end if
        WRITE(NOUT+7,*)'#In columns:'
      WRITE(NOUT+7,*)'#Vacuum lambda_0 and quadrupole extinction'
      write(nout+7,*)
C--------/---------/---------/---------/---------/---------/---------/--
      OPEN(UNIT=NOUT+10,FILE='axs-phmat.dat')
      rewind(NOUT+10)
      WRITE(NOUT+10,5000)
      write(nout+10,*)
      write(nout+10,*)'#Equiv.-volume-sphere radius in nm=', rev
      write(nout+10,*)
        write(NOUT+10,*)'#Angular momentum cut-off LMAX=',LMAX
        write(NOUT+10,*)'#In the number ND=NDGS*LMAX of GIP points,
     & NDGS=',NDGS
        write(NOUT+10,*)'#DEFP=',DEFP
      write(nout+10,*)'#Material number =',NMAT
        if (ynbrug) write(nout+10,*)'#Bruggeman approximation performed'
      WRITE(NOUT+10,1005) THET0,THET,PHI0,PHI,ALPHA,BETA
      write(nout+10,*)

C--------/---------/---------/---------/---------/---------/---------/--
      OPEN(UNIT=NOUT+12,FILE='axs-ampmat.dat')
      rewind(NOUT+12)
      WRITE(NOUT+12,1006)
      write(nout+12,*)
      write(nout+12,*)'#Equiv.-volume-sphere radius in nm=', rev
      write(nout+12,*)
        write(NOUT+12,*)'#Angular momentum cut-off LMAX=',LMAX
        write(NOUT+12,*)'#In the number ND=NDGS*LMAX of GIP points,
     & NDGS=',NDGS
        write(NOUT+12,*)'#DEFP=',DEFP
      write(nout+12,*)'#Material number =',NMAT
        if (ynbrug) write(nout+12,*)'#Bruggeman approximation performed'
      WRITE(NOUT+12,1005) THET0,THET,PHI0,PHI,ALPHA,BETA
      WRITE(NOUT+12,*)'#In columns:'
      WRITE(NOUT+12,*)'#   Vacuum lambda_0, ReVV, ImVV, ReVH, ImVH'
        WRITE(NOUT+12,*)'#   Vacuum lambda_0, ReHV, ImHV, ReHH, ImHH'
      write(nout+12,*)
C--------/---------/---------/---------/---------/---------/---------/--
      OPEN(UNIT=NOUT+13,FILE='axs-intv.dat')
      rewind(NOUT+13)
      write(nout+13,*)'Intensity in theta for
     &      (un)correlated light source'
      write(nout+13,*)
      write(nout+13,*)'#Equiv.-volume-sphere radius in nm=', rev
      write(nout+13,*)
        write(NOUT+13,*)'#Angular momentum cut-off LMAX=',LMAX
        write(NOUT+13,*)'#In the number ND=NDGS*LMAX of GIP points,
     & NDGS=',NDGS
        write(NOUT+13,*)'#DEFP=',DEFP
      write(nout+13,*)'#Material number =',NMAT
        if (ynbrug) write(nout+13,*)'#Bruggeman approximation performed'
      WRITE(NOUT+13,1005) THET0,THET,PHI0,PHI,ALPHA,BETA
        write(nout+13,*)'#host dielectric constant=', zeps(lcs+1)
      if (nmat.ge.1) write(nout,*)'#Dispersive layer number=', ilcs
      if ((lcs.eq.1).and.(nmat.eq.0))
     & write(nout+13,*)'#sphere diel. const.=', zeps(1)
      if (lcs.ge.2) then
      write(nout+13,*)'#core radius/sphere radius =',rff(1)
      if ((ilcs.ne.1).or.(nmat.eq.0))
     &  write(nout+13,*)'#sphere core diel. const.=', zeps(1)
      if ((ilcs.ne.lcs).or.(nmat.eq.0))
     &  write(nout+13,*)'#coating diel. const.=',zeps(lcs)
C--------/---------/---------/---------/---------/---------/---------/--
      end if
      WRITE(NOUT+13,*)'#In columns:'
        WRITE(NOUT+13,*)'#Vacuum lambda_0,|VV|**2+|VH|**2,|VV+VH|**2'
      write(nout+13,*)

C--------/---------/---------/---------/---------/---------/---------/--
      OPEN(UNIT=NOUT+14,FILE='axs-inth.dat')
      rewind(NOUT+14)
      write(nout+14,*)'Intensity in phi for (un)correlated light source'
      write(nout+14,*)
      write(nout+14,*)
      write(nout+14,*)'#Equiv.-volume-sphere radius in nm=', rev
      write(nout+14,*)
        write(NOUT+14,*)'#Angular momentum cut-off LMAX=',LMAX
        write(NOUT+14,*)'#In the number ND=NDGS*LMAX of GIP points,
     & NDGS=',NDGS
        write(NOUT+14,*)'#DEFP=',DEFP
      write(nout+14,*)'#Material number =',NMAT
        if (ynbrug) write(nout+14,*)'#Bruggeman approximation performed'
      WRITE(NOUT+14,1005) THET0,THET,PHI0,PHI,ALPHA,BETA
        write(nout+14,*)'#host dielectric constant=', zeps(lcs+1)
      if (nmat.ge.1) write(nout,*)'#Dispersive layer number=', ilcs
      if ((lcs.eq.1).and.(nmat.eq.0))
     & write(nout+14,*)'#sphere diel. const.=', zeps(1)
      if (lcs.ge.2) then
      write(nout+14,*)'#core radius/sphere radius =',rff(1)
      if ((ilcs.ne.1).or.(nmat.eq.0))
     &  write(nout+14,*)'#sphere core diel. const.=', zeps(1)
      if ((ilcs.ne.lcs).or.(nmat.eq.0))
     &  write(nout+14,*)'#coating diel. const.=',zeps(lcs)
C--------/---------/---------/---------/---------/---------/---------/--
      end if
      WRITE(NOUT+14,*)'#In columns:'
      WRITE(NOUT+14,*)'#Vacuum lambda_0, |HV|**2+|HH|**2, |HV+HH|**2'
      write(nout+14,*)

      OPEN(UNIT=NOUT+15,FILE='tr1diag.dat')
      rewind(NOUT+15)


      Open(10, File='maps.txt',status = 'unknown')

*****************************   ZEPS1  ***********************************
*                  --------------------------------
* Sphere optical constants in the case of a dispersion
* READING IN MATERIAL  DATA:
* Reading real material data, e.g., according to Palik's  book
* requires reading data files OMF and CEPS1 of dimension NFIN
* OMF is reepsz/omega and CEPS1 contains the sphere EPS
*                       material constant reading:
*
      if (nmat.le.1) then

        go to 2         !no reading of material data

      else if (nmat.eq.2) then            ! silver data

      OPEN(UNIT=30,FILE='agc.dat')
      rewind(30)
        do ieps=1,nfin
          read(30,*) omf(ieps),ceps1(ieps)
        enddo
       close(30)

      else if (nmat.eq.3) then        ! Gold data

c      OPEN(UNIT=30,FILE='Au293Knew.dat')       !Gold data for different T
      OPEN(UNIT=30,FILE='Audat.dat')          !Gold data in nm
      write(6,*)'Gold particles'
      rewind(30)
        do ieps=1, nfin
          read(30,*) omf(ieps),ceps1(ieps)
c          omf(ieps)=2.d0*pi*rev*omf(ieps)/(1240.d0*rmuf)
          omf(ieps)=2.d0*pi*rev/(omf(ieps)*rmuf)
        enddo
       close(30)

cc      else if (nmat.eq.4) then

      else if (nmat.eq.5) then        ! Copper data

      OPEN(UNIT=30,FILE='Cudat.dat')          !Copper data in nm
      write(6,*)'Copper particles'
      rewind(30)
        do ieps=1, nfin
          read(30,*) omf(ieps),ceps1(ieps)
          omf(ieps)=2.d0*pi*rev/(omf(ieps)*rmuf)
        enddo
      close(30)

      else if (nmat.eq.6) then        ! Aluminium data

      OPEN(UNIT=30,FILE='Aldat.dat')          !Aluminium data in nm
      write(6,*)'Aluminum particles'
      rewind(30)
        do ieps=1, nfin
          read(30,*) omf(ieps),ceps1(ieps)
          omf(ieps)=2.d0*pi*rev/(omf(ieps)*rmuf)
        enddo
      close(30)

      else if (nmat.eq.7) then        ! Platinum data

      OPEN(UNIT=30,FILE='Ptdat.dat')          !Platinum data in nm
      write(6,*)'Platinum particles'
      rewind(30)
        do ieps=1, nfin
          read(30,*) omf(ieps),ceps1(ieps)
          omf(ieps)=2.d0*pi*rev/(omf(ieps)*rmuf)
        enddo
      close(30)

      else if (nmat.eq.8) then        ! Silicon data

c     OPEN(UNIT=30,FILE='sieps.dat')  !Silicon data in nm
        OPEN(UNIT=30,FILE='Sidat.dat')   !Silicon data in nm for larger interval
      write(6,*)'Silicon particles'
      rewind(30)
        do ieps=1, nfin
          read(30,*) omf(ieps),ceps1(ieps)
C          ceps1imag(ieps)= Aimag(ceps1(ieps))
c         ceps1real(ieps) = REAL(ceps1(ieps))
c          ceps1imag(ieps)=0
c          ceps1(ieps) = CMPLX(ceps1real(ieps), ceps1imag(ieps))
          omf(ieps)=2.d0*pi*rev/(omf(ieps)*rmuf)
        enddo
      close(30)

       else if (nmat.eq.9) then        ! Water data


         OPEN(UNIT=30,FILE='measured_Water_dispersion_T=24.txt')   !Water in GHz
      write(6,*)'Water particles'
      rewind(30)
        do ieps=1, nfin
          read(30,*) omf(ieps),ceps1real(ieps),ceps1imag(ieps)
          ceps1(ieps) = CMPLX(ceps1real(ieps), ceps1imag(ieps))
          omf(ieps)=2.d0*pi*rev/((3.d8/omf(ieps))*rmuf)
        enddo
      close(30)

      end if                      ! material constant reading

*********************

*                     --------------------------------
* begin main scanning loop:
      write(6,*) 'begining of main loop'
   2   rev_beg = rev



       do 199 itter = 1,ndefp

      DEFP=rsnm/hlength_max + (rsnm/hlength_min - rsnm/hlength_max)*
     & DBLE(itter-1)/DBLE(ndefp-1)
      DEFPP = DEFP
      hlength = rsnm/DEFP
      rev=hlength*(3.D0*DEFP*DEFP/2.D0)**(1.D0/3.D0)
       do 200 istep=1,nstep




      omega=omega0 + dble(istep-1)*delo


      write(6,*) itter , istep, DEFP, hlength
C     omega_max = omega0 + dble(nstep)*delo
c     lambda_min=2.d0*pi*REV/(omega_max*RMUF)
c      lambda_max=2.d0*pi*REV/(omega0*RMUF)
c      lambda = lambda_min + (lambda_max - lambda_min)*DBLE(istep)
c     & /dble(NSTEP)
      lambda=2.d0*pi*rev_beg/(omega*RMUF)

cc      xs=RMUF*omega*dble(sqrt(zeps0))  !Equiv. size parameter

      if ((nmat.eq.0).or.(ynperfcon)) go to 8      !dispersionless dielectric
                                                   !or ideal metal
* In case of a dispersion, EPSSPH is modified.
* For ideal Drude metal
*     plasma=2.d0*pi*sphere radius in nm/(lambda_z in nm*rmuf)
* where lambda_z is the wavelength for which Re eps_s=0.

       reepsz=2.d0*pi*rev/(323.83d0*rmuf)

      IF (NMAT.EQ.1) THEN              !Material decision IF - Drude metal

      plasma=reepsz
        omxp=plasma/omega
        zeps1=1.d0-omxp**2/(1.d0+ci*plasma/(144.d0*omega))
      go to 5
*
        ELSE IF (nmat.eq.4) then             !Material decision IF - ZnS
*
       filfrac=0.62d0         ! filfrac of ZnS in ZnS core
       call  znsrefind(LAMBDA,FILFRAC,zeps1)
       go to 5
*
      ELSE IF (NMAT.EQ.2) THEN         !Material decision IF - Ag

c >>> real material data:           !silver
*                         lambda_z=323.83d0
*                         lambda_p=164.d0
* When real material data are used,
* reepsz differs from plasma!!! The plasma wavelength is
* calculated below:

       plasma=reepsz*7.2d0/3.8291d0

* security trap - remainder (not optimized!)
      omxf=omega/reepsz
      if (omxf.gt.omf(1)) then
       write(6,*)'Calculation of has to stop with'
       write(6,*)' OMF(1)'
       write(6,*)' OMXF=', omxf
       stop
      end if

      if (omxf.lt.omf(nfin)) then
        omxp=plasma/omega
        zeps1=1.d0-omxp**2/(1.d0+ci*plasma/(144.d0*omega))
* damping coefficient for silver is plasma/144 where plasma is different from
* the Re eps zero crossing at 3.8291 eV according to Palik!!!
       go to 5
      else if (omxf.eq.omf(1)) then
       zeps1=ceps1(1)
       go to 5
      else
      do ieps=2,nfin
* data file ordered with the increased wavelength
* omxf increases in the loop and is oriented opposite to the data file
       if (omxf.gt.omf(ieps)) then     ! linear interpolation
       zeps1=ceps1(ieps)+(omxf-omf(ieps))*(ceps1(ieps-1)-ceps1(ieps))
     1 /(omf(ieps-1)-omf(ieps))
       go to 5
       end if
      enddo
       end if   !end Ag

      ELSE IF ((NMAT.GE.3).or.((nmat.ge.5).and.(nmat.le.7))) then   !Material decision IF
                                                                    !Au,Cu,Al,Pt
c >>>
* data file ordered with the decreased wavelength
* omega increases in the loop and is oriented along the data file
*
      if ( (omega.lt.omf(1)).or.(omega.gt.omf(nfin)) ) then
cc       write(6,*)'Material data not available for this wavelength'
cc       stop
*
      call sordalc(NMAT,lambda,ZEPS1)
      go to 5
*
      end if
*
      if (omega.eq.omf(nfin)) then
       zeps1=ceps1(nfin)
       go to 5
      else
      do ieps=1,nfin-1
       if (omega.lt.omf(ieps+1)) then     ! linear interpolation
       zeps1=ceps1(ieps)+(omega-omf(ieps))*(ceps1(ieps+1)-ceps1(ieps))
     1 /(omf(ieps+1)-omf(ieps))
       go to 5
       end if
      enddo
      end if

      ELSE IF (NMAT.EQ.8) then           !Material decision IF - Silicon
c >>>
* data file ordered with the decreased wavelength
* omega increases in the loop and is oriented along the data file
*
      if ( (omega.lt.omf(1)).or.(omega.gt.omf(nfin)) ) then
       write(6,*)'Material data not available for this wavelength'
       stop
*
      end if
*
      if (omega.eq.omf(nfin)) then
       zeps1=ceps1(nfin)
       go to 5
      else
      do ieps=1,nfin-1
       if (omega.lt.omf(ieps+1)) then     ! linear interpolation
       zeps1=ceps1(ieps)+(omega-omf(ieps))*(ceps1(ieps+1)-ceps1(ieps))
     1 /(omf(ieps+1)-omf(ieps))
       go to 5
       end if
      enddo
      end if

      ELSE IF (NMAT.EQ.9) then           !Material decision IF - water
c >>>
* data file ordered with the decreased wavelength
* omega increases in the loop and is oriented along the data file
*
      if ( (omega.lt.omf(1)).or.(omega.gt.omf(nfin)) ) then
       write(6,*)'Material data not available for this wavelength'
       stop
*
      end if
*
      if (omega.eq.omf(nfin)) then
       zeps1=ceps1(nfin)
       go to 5
      else
      do ieps=1,nfin-1
       if (omega.lt.omf(ieps+1)) then     ! linear interpolation
       zeps1=ceps1(ieps)+(omega-omf(ieps))*(ceps1(ieps+1)-ceps1(ieps))
     1 /(omf(ieps+1)-omf(ieps))
       go to 5
       end if
      enddo
      end if


      END IF                  ! END of Material decision IF

* The end of reading real data according to Palik's  book
*_____________________________________
* activate Bruggeman:




  5   if (ynbrug) then
      ff=0.8d0
      z1 = (3.d0*ff-1.d0)*zeps1+(2.d0 - 3.d0*ff)*zeps0
      z2 =  sqrt(z1*z1 + 8.d0*zeps1*zeps0)
*
       if (IMAG(z2).GE.0.0) then
         zeps1= (z1 + z2)/4.d0
       else
         zeps1= (z1 - z2)/4.d0
       end if
       end if

      zeps(ilcs)=zeps1
*______________________________________

  8   continue



cc      write(6,*)'LAMBDA in AXSPARTCL=', LAMBDA

      IF (YNINTENS) goto 500
*
      if ((np.eq.-1).and.(lcs.gt.1)) then

        revinl=revin*2.d0*pi*sqrt(zeps0)/lambda

        revl=rev*2.d0*pi*sqrt(zeps0)/lambda

ctest
cc      defpp=0.5
cc      revinl=0.5
cc      revl=10.d0
ctest
      if (yncheck) then
        ide=2
        else
        ide=4
        end if

      call lisac(ide,lcs,lmax,lmax,ndgs*lmax,lambda,defpp,
     & revinl/revl,revl,zeps)

      else
      global_eps_r =REAL( zeps(1))
      global_eps_i = AIMAG(zeps(1))


C      write(90,*) lambda,  global_eps_r,
C     & global_eps_i
      call ampldr(yncheck,lmax,ichoice,npp,defpp,rsnm,hlength,lambda,
     &  zeps(1),zeps0)
      end if
*
C--------/---------/---------/---------/---------/---------/---------/--
cc      if(nstep.gt.10) go to 200
c      write(6,*) 'istep=', istep
cc      write(6,*) 'lambda=',lambda
cc      write(6,*)
cc      write(6,*)'Scattering coefficient='
cc      write(6,*)'tsc=', tsc

 200  continue
 199  continue

* <<<
      close(nout)
      close(nout+1)
        close(nout+2)
        close(nout+3)
      close(nout+5)
        close(nout+6)
        close(nout+7)
      close(nout+10)
      close(nout+12)
      close(nout+13)
      close(nout+14)
      close(nout+15)
      close(10)
      close(90)
* <<<
      if (ync .eq.'n') then
        write(6,*)'Homogeneous particle'
      else if (ync.eq.'y') then
       write(6,*)'coated particle'
      end if
*
      write(6,*)'Particle parameters:'
      write(6,*)
      write(6,*)'Equivalent sphere radius =', rev
      if (ync.eq.'n') write(6,*)'particle diel. constant=', zeps(1)
      write(6,*)'background dielectric constant ZEPS0=', zeps0
      if (ync.eq.'y') write(6,*)'core diel. constant=', zeps(1)
      if (ync.eq.'y') write(6,*)'coating diel. constant=', zeps(lcs)
      if (ync.eq.'y') write(6,*)'core radius/particle radius =',
     & rff(1)
      write(6,*)
        write(6,*)'OA scattering cs versus wavelength in axs-scs.dat'
      write(6,*)'OA Extinction versus wavelength in axs-oavext.dat'
      write(6,*)'OA Absorption versus wavelength in axs-abs.dat'
        write(6,*)'OA Albedo versus wavelength in axs-albedo.dat'
        write(6,*)'  [3rd column displays qext-(qsca +qabs)]  '
        write(6,*)'Extinction versus wavelength in axs-ext.dat'
      write(6,*)'Phase matrix vs wavelength in axsphmat.dat'
      write(6,*)'Amplitude matrix vs wavelength in axsampmat.dat'
      write(6,*)'OA Dipole extinction in axs-dipolext.dat'
      write(6,*)'OA Quadrupole extinction in axs-quadrext.dat'
*--------/---------/---------/---------/---------/---------/---------/--

      IF (.NOT.YNINTENS) go to 1000

C--------/---------/---------/---------/---------/---------/---------/--
*
* output initial statements

 500  OPEN(UNIT=NOUTI,FILE='intnsty.dat')
      rewind(NOUTI)
      WRITE(NOUTI,*)'#Field intensity profile around a tip'
      write(NOUTI,*)
      write(NOUTI,*)'#Equiv.-volume-sphere radius in nm=', rev
      write(NOUTI,*)
        write(NOUTI,*)'#Angular momentum cut-off LMAX=',LMAX
        write(NOUTI,*)'#In the number ND=NDGS*LMAX of GIP points,
     & NDGS=',NDGS
        write(NOUTI,*)'#DEFP=',DEFP
      write(NOUTI,*)'#Material number =',NMAT
      if (ync .eq.'n') then
        write(NOUTI,*)'#Homogeneous particle'
      else if (ync.eq.'y') then
       write(NOUTI,*)'#coated particle'
      end if
      write(NOUTI,*)
      write(NOUTI,*)'#host dielectric constant=', zeps(lcs+1)
      if (ync.eq.'n')
     &   write(NOUTI,*)'#particle diel. constant=', zeps(1)
      if (ync.eq.'y')
     &   write(NOUTI,*)'#coating diel. constant=',zeps(lcs)
C--------/---------/---------/---------/---------/---------/---------/--
      if (ync.eq.'y')
     & write(NOUTI,*)'#core radius/particle radius =',rff(1)
      WRITE(NOUTI,*)'#phi,cos(theta),r(theta), and |E| in columns'
      write(NOUTI,*)

        OPEN(UNIT=NOUTI+1,FILE='elfcomp.dat')
      rewind(NOUTI+1)
      WRITE(NOUTI+1,*)'#Total electric field components around a tip'
      write(NOUTI+1,*)
      write(NOUTI+1,*)'#Equiv.-volume-sphere radius in nm=', rev
      write(NOUTI+1,*)
      write(NOUTI+1,*)'#Angular momentum cut-off LMAX=',LMAX
        write(NOUTI+1,*)'#In the number ND=NDGS*LMAX of GIP points,
     & NDGS=',NDGS
        write(NOUTI+1,*)'#DEFP=',DEFP
      write(NOUTI+1,*)'#Material number =',NMAT
      if (ync .eq.'n') then
        write(NOUTI+1,*)'#Homogeneous particle'
      else if (ync.eq.'y') then
       write(NOUTI+1,*)'#coated particle'
      end if
      write(NOUTI+1,*)
      write(NOUTI+1,*)'#host dielectric constant=', zeps(lcs+1)
      if (ync.eq.'n')
     &   write(NOUTI+1,*)'#particle diel. constant=', zeps(1)
      if (ync.eq.'y')
     &   write(NOUTI+1,*)'#coating diel. constant=',zeps(lcs)
C--------/---------/---------/---------/---------/---------/---------/--
      if (ync.eq.'y')
     & write(NOUTI+1,*)'#core radius/particle radius =',rff(1)
      WRITE(NOUTI+1,*)'#phi,cos(theta),r(theta), E_r,E_theta,E_phi
     & in columns'
      write(NOUTI+1,*)
*
        OPEN(UNIT=NOUTI+2,FILE='elfscat.dat')
      rewind(NOUTI+2)
      WRITE(NOUTI+2,*)'#Scattered E-field components around a tip'
      write(NOUTI+2,*)
      write(NOUTI+2,*)'#Equiv.-volume-sphere radius in nm=', rev
      write(NOUTI+2,*)
      write(NOUTI+2,*)'#Angular momentum cut-off LMAX=',LMAX
        write(NOUTI+2,*)'#In the number ND=NDGS*LMAX of GIP points,
     & NDGS=',NDGS
        write(NOUTI+2,*)'#DEFP=',DEFP
      write(NOUTI+2,*)'#Material number =',NMAT
      if (ync .eq.'n') then
        write(NOUTI+2,*)'#Homogeneous particle'
      else if (ync.eq.'y') then
       write(NOUTI+2,*)'#coated particle'
      end if
      write(NOUTI+2,*)
      write(NOUTI+2,*)'#host dielectric constant=', zeps(lcs+1)
      if (ync.eq.'n')
     &   write(NOUTI+2,*)'#particle diel. constant=', zeps(1)
      if (ync.eq.'y')
     &   write(NOUTI+2,*)'#coating diel. constant=',zeps(lcs)
C--------/---------/---------/---------/---------/---------/---------/--
      if (ync.eq.'y')
     & write(NOUTI+2,*)'#core radius/particle radius =',rff(1)
      WRITE(NOUTI+2,*)'#phi,cos(theta),r(theta), E_r,E_theta,E_phi
     & in columns'
      write(NOUTI+2,*)
*
      call evanesc(lmax,npp,lcs,revp,defpp,rsnm,hlength,lambda,thetv,
     & rmf,zeps)
*
      close(nouti)
      close(nouti+1)
      close(nouti+2)

      write(6,*)'Particle parameters:'
      write(6,*)
      write(6,*)'Equivalent sphere radius =', rev
      if (ync.eq.'n') write(6,*)'particle diel. constant=', zeps(1)
      write(6,*)'background dielectric constant ZEPS0=', zeps0
      if (ync.eq.'y') write(6,*)'core diel. constant=', zeps(1)
      if (ync.eq.'y') write(6,*)'coating diel. constant=', zeps(lcs)
      if (ync.eq.'y') write(6,*)'core radius/particle radius =',
     & rff(1)
      write(6,*)
      write(6,*)'Total |E|**2 in intnsty.dat'
      write(6,*)'Total electric field components around a tip
     & in elfcomp.dat'
      write(6,*)'Scattered electric field components around
     &  a tip in elfscat.dat'

 1000 CONTINUE


 1005 FORMAT ('thet0=',F6.2,'  thet=',F6.2,'  phi0=',F6.2,
     &        '  phi=',F6.2,'  alpha=',F6.2,'  beta=',F6.2)
 1006 FORMAT ('AMPLITUDE or S-MATRIX')
 5000 FORMAT ('4X4 PHASE MATRIX')

      end



      SUBROUTINE TMTAXSPV(nmax,lambda,rsnm,ht,zeps1,zeps0,TMT)
C--------/---------/---------/---------/---------/---------/---------/--
c Warning in module TMTAXSP in file tmtaxsp.f: Variables set but never used:
c    NGGG set at line 182 file tmtaxsp.f
c Warning in module TMTAXSP in file tmtaxsp.f: Variables may be used before set:
c    QEXT1 used at line 215 file tmtaxsp.f
c    QEXT1 set at line 220 file tmtaxsp.f
c    QSCA1 used at line 214 file tmtaxsp.f
c    QSCA1 set at line 221 file tmtaxsp.f
C--------/---------/---------/---------/---------/---------/---------/--
C NMAX - angular momentum cut off
C LAMBDA - vacuum wavelength
C RAP=S(1,1)*KAPPA0/2.D0/PI     !=rmuf*ALPHA/LAMBDA =rsnm/LAMBDA
C
C    RETURNS the T matrix of a general axially symmetric scatterer
C    The latter has also non-zero of-diagonal (mixed EH and HE terms):
C
C            |  TMT(2,*) |  TMT(3,*)   |   |  TMT(MM) |  TMT(ME)   |
C    TMT  =  | ----------+-------------| = |----------+------------|
C            |  TMT(4,*) |  TMT(1,*)   |   |  TMT(EM) |  TMT(EE)   |
C
C    TMT(1,*) terms corresponds to TEE scattering matrices
C    TMT(2,*) terms corresponds to TMM scattering matrices
C    TMT(3,*) terms corresponds to TME scattering matrices
C    TMT(4,*) terms corresponds to TEM scattering matrices
C    TMT(4,*)=-TMT(3,*)^t where t denotes transposed TMT(3,*) submatrix
C
C    TMT's equal to i*sin(eta)*exp(i*eta), where eta is a phase-shift
C             ====>    S=1+2*T=exp(2*i*eta)
C    and the unitarity of S-matrix implies
C            T^\dagger*T=T*T^\dagger=-(1/2)*(T^\dagger+T)
C______________________________
C
C LOCAL VARIABLES:
C ===============
C ICHOICE=1 if NAG library is available, otherwise ICHOICE=2
C
C NP,EPS: specifies the shape of particles within a given NP class:
C     NP.gt.0 - EPS = deformation parameter of a Chebyshev particle
C     NP=-1 - EPS = the ratio of the horizontal to rotational axes. EPS is
C             larger than 1 for oblate spheroids and smaller than 1 for
C             prolate spheroids.
C     NP=-2 - EPS = the ratio of the cylinder diameter to its length.
C     NP=-3 - no EPS is specified
C     NP=-4 - EPS is the height (along the axial symmetry axis)
C             of the resulting cut sphere
C                Note that always EPS.LT.2*REV specified
C
C Warning:
C  In computations for spheres, use EPS=1.000001 instead of EPS=1.
C  EPS=1 can cause overflows in some rare cases.
C
C  LAM - the wavelength of incident light in the ambient.
C
C  RAT = 1 - particle size is specified in terms of the
C                equal-volume-sphere radius
C  RAT.ne.1 - particle size is specified in terms of the
C                equal-surface-area-sphere radius
C  AXI ... equivalent-(volume/surface-area)-sphere radius
C  REV=A=RAT*AXI ... equal-volume-sphere radius
C                  (feeded as REV to RSP* routines)
C  DDELT - required precision
C
C  ALPHA and BETA - Euler angles (in degrees) specifying the
C          orientation of the scattering particle relative to
C          the laboratory reference frame (Refs. 6 and 7).
C
C  For axially symmetric scatterers, when the T matrix is computed in
C  natural coordinate system with the $z$ axis along the axis of particle
C  axial symmetry, one can show that the T matrix is {\em diagonal} with
C  respect to the azimuthal indices $m$ and $m'$ \cite{Wat},
C
C              T_{lm,l'm'}^{ij}=\delta_{mm'} T_{lm,l'm},
C
C  and that it satisfies reciprocity relation \cite{GuS,Mis36},
C
C               T_{lm,l'm}^{ij}=(-1)^{i+j} T_{l'm,lm}^{ji}.
C
C  \cite{Mis91} also claims the relation:
C
C                T_{lm,l'm}^{ij}= (-1)^{i+j} T_{l-m,l'-m}^{ij}
C
C--------/---------/---------/---------/---------/---------/---------/--
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER LMAXD,LMAX1D,LMTD
      INTEGER NAXSM,ICHOICEV,ICHOICE

      PARAMETER (LMAXD=50,LMAX1D=LMAXD+1,LMTD=LMAX1D*LMAX1D-1)

*      INCLUDE 'ampld.par.f'
      Parameter (NPN1=100, NPNG1=2600, NPNG2=2*NPNG1, NPN2=2*NPN1,
     &           NPL=NPN2+1, NPN3=NPN1+1,
     &           NPN4=NPN1, NPN5=2*NPN4, NPN6=NPN4+1)
* number of the output unit
*
      REAL*8  LAM,LAMBDA,MRR,MRI,RSNM,HT,X(NPNG2),W(NPNG2),
     *        S(NPNG2),SS(NPNG2),AN(NPN1),R(NPNG2),DR(NPNG2),
     *        DDR(NPNG2),DRR(NPNG2),DRI(NPNG2),ANN(NPN1,NPN1)
      REAL*8 TR1(NPN2,NPN2),TI1(NPN2,NPN2)
c      REAL*8 XALPHA(300),XBETA(300),WALPHA(300),WBETA(300)

      COMPLEX*16 CZERO,COPTE,COPTM
      COMPLEX*16 zeps1,zeps0
      COMPLEX*16 TMT(4,LMTD,LMTD)
*
      logical ynoptth
*
      COMMON /CT/ TR1,TI1
* transfers the real and imaginary part of the T matrix (2*NMAX,2*NMAX)
* array for a given value of M from TT via (TMATR0 and TMATR) to the main
*
cc      COMMON /TMAT/ RT11,RT12,RT21,RT22,IT11,IT12,IT21,IT22
* transfers T matrix arrays obtained from TR1,TI1 in the main to the
* AMPL routine --->  NOT USED HERE
*

      COMMON /CHOICE/ ICHOICE
* transfers the choice of inversion from here to TT
*
      COMMON /TOITMT/ICHOICEV,NP,NCHECK,NAXSM,NDGS
*
* transfers integers ICHOICEV,NP,NCHECK,NAXSM,NDGS from the main here
*
      COMMON /TOTMT/EPS,RAT,REV,ALPHA,BETA,DDELT
*
* transfers real*8 EPS(DEFP),RAT,A(REV),ALPHA,BETA,DDELT
* from the main here
*
      COMMON /TOLTMT/ ynoptth
*
* transfers logical ynoptth from the main here
*****************************************************************
      DATA CZERO/(0.D0,0.D0)/
*
      P=DACOS(-1D0)                 !local PI constant
*
      ICHOICE=ICHOICEV
      A=REV
      LAM=LAMBDA/SQRT(ZEPS0)       !vacuum wavelength devided by SQRT(ZEPS0)

      write(6,*)'LAM,LAMBDA in TMTAXSP=', LAM, LAMBDA
*
* the real part of the refractive index contrast
*
      MRR=DBLE(SQRT(ZEPS1/ZEPS0))
*
* the imaginary  part of the refractive index contrast
*
      MRI=DIMAG(SQRT(ZEPS1/ZEPS0))
*
      DDELT=0.1D0*DDELT               !conv. test is switched off now!!!
*
* DDELT is used to test the accuracy of computing the
* optical cross sections. This accuracy is usually better
* than the absolute accuracy of computing the expansion coefficients
* of a normalized scattering matrix by a factor of 10. Therefore,
* the desired accuracy of computing the expansion coefficients
* is rescaled by a factor 0.1 before entering the test of the
* accuracy of computing the optical cross sections.
*
* Other local constants:
*
      LMTOT=(NMAX+1)**2-1

      IF (DABS(RAT-1D0).GT.1D-8.AND.NP.EQ.-1) CALL SAREA (EPS,RAT)
      IF (DABS(RAT-1D0).GT.1D-8.AND.NP.GE.0) CALL SURFCH(NP,EPS,RAT)
      IF (DABS(RAT-1D0).GT.1D-8.AND.NP.EQ.-2) CALL SAREAC (EPS,RAT)
      IF (NP.EQ.-3) CALL DROP (RAT)

*___________________________________________________
* Determination of the Wiscombe value of the floating
C angular momentum cutoff NMAX:

      XEV=2D0*P*A/LAM                     !size parameter
      IXXX=XEV+4.05D0*XEV**0.333333D0     !Wiscombe conv. criterion for NMAX
      IF (XEV.GT.1.) THEN
      INM1=MAX0(3,IXXX)
      ELSE
      INM1=2                 !The Bessel package routine RYB
                             !requires NMAX to be at least 2
      END IF
*
      IF (INM1.GE.NPN1) PRINT 7333, NPN1
      IF (INM1.GE.NPN1) STOP
 7333 FORMAT('CONVERGENCE IS NOT OBTAINED FOR NPN1=',I3,
     &       '.  EXECUTION TERMINATED')

      NGAUSS=NMAX*NDGS
cc      NNNGGG=NGAUSS+1

      IF (NGAUSS.EQ.NPNG1) PRINT 7336
 7336    FORMAT('WARNING: NGAUSS=NPNG1')
*
* GIF division points and weights + other numerical constants
*
        CALL CONST(NGAUSS,NMAX,X,W,AN,ANN,S,SS,NP,EPS,RSNM,HT)    !In TMTAXSPV
*
* specify particle shape:
*
        CALL VARY(LAM,MRR,MRI,A,EPS,RSNM,HT,NP,NGAUSS,X,P,
     &              PPI,PIR,PII,R,DR,DDR,DRR,DRI,NMAX)
*
* determine m=m'=0 elements of the T matrix
*
         CALL TMATR0 (NGAUSS,X,W,AN,ANN,PPI,PIR,PII,R,DR,
     &                 DDR,DRR,DRI,NMAX,NCHECK,NAXSM)
*
cc      WRITE(NOUT,*)'NMAX=',NMAX
cc      WRITE(NOUT,*)'NGAUSS=',NGAUSS
*<<<
*
* TMT initialization:

cc         DO JA=1,LMTOT
cc         DO JB=1,LMTOT
*
cc            TMT(1,JB,JA)=CZERO
cc            TMT(2,JB,JA)=CZERO
cc            TMT(3,JB,JA)=CZERO
cc            TMT(4,JB,JA)=CZERO
cc
cc         ENDDO
cc         ENDDO
*
C            |  TMT(2,*) |  TMT(3,*)   |   |  TMT(MM) |  TMT(ME)   |
C    TMT  =  | ----------+-------------| = |----------+------------|
C            |  TMT(4,*) |  TMT(1,*)   |   |  TMT(EM) |  TMT(EE)   |
C
C    TMT(2,*) terms corresponds to TMM scattering sub-matrix
C    TMT(3,*) terms corresponds to TME scattering sub-matrix
C    TMT(4,*) terms corresponds to TEM scattering sub-matrix
C    TMT(1,*) terms corresponds to TEE scattering sub-matrix
C    TMT(4,*)=-TMT(3,*)^t where t denotes transposed TMT(3,*) sub-matrix

*****************     Assign  m=m=0 elements of TMT matrix   ***********

         DO L1=1,NMAX
         DO L2=1,NMAX
          N1=L1+NMAX
          N2=L2+NMAX

            JA=L1*(L1+1)       ! (l,m) index with (1-1)=1
            JB=L2*(L2+1)

* see (5.39) of {MTL}:  !!!Iff plane of symmetry perpendicular to the
                        !!!axis of rotation

          if ((NAXSM.eq.1).and.((-1)**(L1+L2).ne.1)) then
            TMT(2,JA,JB)=CZERO
            TMT(1,JA,JB)=CZERO
          else
            TMT(2,JA,JB)=DCMPLX(TR1(L1,L2),TI1(L1,L2))
            TMT(1,JA,JB)=DCMPLX(TR1(N1,N2),TI1(N1,N2))
          end if
* see (5.37) of {MTL}:
            TMT(4,JA,JB)=CZERO         !DCMPLX(TR1(N1,L2),TI1(N1,L2))
            TMT(3,JB,JA)=CZERO         !-TMT(4,JA,JB)

cd      if (ja.eq.2) then
cd         write(6,*) 'jb, tmt(2,2,jb)=',jb, tmt(2,ja,jb),jb
cd         write(6,*) 'jb, tmt(3,jb,2)=',jb, tmt(3,jb,ja),jb
cd      end if

         ENDDO
         ENDDO
*

*****************    Assign  m=m'>0 elements of the T matrix   ***********

      DO 220 M=1,NMAX
*
c         CALL TMATR(M,NGAUSS,X,W,AN,ANN,S,SS,PPI,PIR,PII,R,DR,
c     &               DDR,DRR,DRI,NMAX,NCHECK,NAXSM)

         CALL TMTR(M,NGAUSS,X,W,AN,ANN,PPI,PIR,PII,R,DR,DDR,DRR,
     &               DRI,NMAX,NCHECK,NAXSM)
*
* <<< returns  m=m'>0 elements of the T matrix
* For a given M, TMATR returns Re T and Im T matrices ordered from
* 1 to 2*(NMAX-M+1)
*
         NM=NMAX-M+1             !size of a T block returned by TT

         DO L1=M,NMAX
         DO L2=M,NMAX

          K1=L1-M+1               !K1,K2,KK1,KK2 label the entries of
          K2=L2-M+1               !a TT returned T block
          KK1=L1-M+1+NM
          KK2=L2-M+1+NM


            JA =L1*(L1+1)+M       ! (l,m) index with (1-1)=1
            JB =L2*(L2+1)+M
            JAM=L1*(L1+1)-M       ! (l,m) index with (1-1)=1
            JBM=L2*(L2+1)-M
*
* see (5.39) of {MTL}: !!!Iff plane of symmetry perpendicular to the
                       !!!axis of rotation

          if ((NAXSM.eq.1).and.((-1)**(L1+L2).ne.1)) then
            TMT(2,JA,JB)  =CZERO
            TMT(2,JAM,JBM)=CZERO
            TMT(1,JA,JB)  =CZERO
            TMT(1,JAM,JBM)=CZERO
          else
            TMT(2,JA,JB)   = DCMPLX(TR1(K1,K2),TI1(K1,K2))
            TMT(2,JAM,JBM) = TMT(2,JA,JB)
            TMT(1,JA,JB)   = DCMPLX(TR1(KK1,KK2),TI1(KK1,KK2))
            TMT(1,JAM,JBM) = TMT(1,JA,JB)
          end if

          if ((NAXSM.eq.1).and.((-1)**(L1+L2).ne.-1)) then
            TMT(4,JA,JB)  =CZERO
            TMT(4,JAM,JBM)=CZERO
            TMT(3,JA,JB)  =CZERO
            TMT(3,JAM,JBM)=CZERO
          else
            TMT(4,JA,JB)   = DCMPLX(TR1(KK1,K2),TI1(KK1,K2))
            TMT(4,JAM,JBM) =-TMT(4,JA,JB)
            TMT(3,JB,JA)   =-TMT(4,JA,JB)
            TMT(3,JBM,JAM) = TMT(4,JA,JB)    !=-TMT(3,JB,JA)
          end if
*
*  Using reciprocity (Eq. (15) of Ref. \ct{Mis97}):
*
*       T_{lm,l'm}^{ij}=(-1)^{i+j} T_{l'm,lm}^{ji}
*
*  and (see Eq. (36) \JQSRT{55}):
*
*       T_{lm,l'm'}^{ij}=(-1)^{m+m'} T_{l'-m',l-m}^{ji}
*
*  Moreover, for axially symmetric particles one has
*  (see Eq. (31),(36) \JQSRT{55}):
*
*  T_{lm,l'm}^{ij}=(-1)^{i+j} T_{l-m,l'-m}^{ij} = T_{l'-m',l-m}^{ji}
*
         ENDDO
         ENDDO

  220 CONTINUE    !end of loop over m's

* Activate the cz-part below if you want to see the diagonal elements
* of the T-matrix

cz      JB =(NMAX+1)**2-1
cz      do 230 JA=1,min(24,JB)
cz      write(16,*)'JA,TMT(1,JA,JA)=',JA,TMT(1,JA,JA)
cz      write(16,*)'JA,TMT(2,JA,JA)=',JA,TMT(2,JA,JA)

cz  230 CONTINUE    !end of loop over JA


      if (.not.ynoptth) go to 400
*  Check optical theorem
*
* \sum_{p_1l_1} \left|T_{l_1m,lm}^{p_1p}\right|^2
*              = -\mb{Re}\, T_{lm,lm}^{pp}.
*

      DO 300 M=0,NMAX
      do 300 L1=MAX(1,M),NMAX

      N1=M

      do 280 i=1,2

      JA =L1*(L1+1)+N1       ! (l,m) index with (1-1)=1

      copte=czero
      coptm=czero

      DO 250 L2=MAX(1,M),NMAX

      JB =L2*(L2+1)+N1

      copte=copte+ ABS(TMT(1,JB,JA))**2+ ABS(TMT(3,JB,JA))**2
      coptm=coptm+ ABS(TMT(2,JB,JA))**2+ ABS(TMT(4,JB,JA))**2

  250 CONTINUE

      copte=copte+dble(TMT(1,JA,JA))

      if (abs(copte).gt.1.4d-2) then
      write(6,*)'Optical theorem violated for EE-mode and
     & L1=',L1,'M=',N1
      write(6,*)'abs(copte)=',abs(copte)
      write(6,*)'(True violation only in the absorptionless case!)'
      pause
      end if

      coptm=coptm+dble(TMT(2,JA,JA))

      if (abs(coptm).gt.1.4d-2) then
      write(6,*)'Optical theorem violated for MM-mode and
     & L1=',L1,'M=',N1
      write(6,*)'abs(coptm)=',abs(coptm)
      write(6,*)'(True violation only in the absorptionless case!)'
      pause
      end if
*
      if (n1.ne.0) then
            N1=-N1
      else
            go to 300
      end if
*
  280 CONTINUE    !polarization loop
*
  300 CONTINUE    !optical theorem check

  400 RETURN
      END

C**********************************************************************

      SUBROUTINE AMPLDR(yncheck,nmax,ichoicev,np,eps,rsnm,ht,lambda,
     1                  zeps1,zeps0)

C Warning in module AMPLDR in file ampldr.f: Variables set but never used:
C    NGGG set at line 493 file ampldr.f
C--------/---------/---------/---------/---------/---------/---------/--
C YNCHECK=.TRUE. if you want to check Gauss integrations
C convergence; otherwise YNCHECK=.FALSE.
C NMAX - angular momentum cut off
C LAMBDA - the vacuum wavelength
C
C Outputs to common block T matrix
C
C                    |  TMT(M,M) |  TMT(M,E)   |
C            TMT  =  | ----------+-------------|
C                    |  TMT(E,M) |  TMT(E,E)   |
C
C    TMT(1,*) corresponds to TEE scattering matrix
C    TMT(2,*) corresponds to TMM scattering matrix
C    TMT(3,*) corresponds to TME scattering matrix
C    TMT(4,*) corresponds to TEM scattering matrix
C    TMT(4,*)=-TMT(3,*)^t where t denotes transposed TMT(3,*) submatrix
C
C ICHOICE=1 if NAG library is available, otherwise ICHOICE=2
C
C NP,EPS: specifies the shape of particles within a given NP class:
C     NP.gt.0 - EPS = deformation parameter of a Chebyshev particle
C     NP=-1 - EPS = the ratio of the horizontal to rotational axes. EPS is
C             larger than 1 for oblate spheroids and smaller than 1 for
C             prolate spheroids.
C     NP=-2 - EPS = the ratio of the cylinder diameter to its length.
C     NP=-3 - no EPS is specified
C     NP=-4 - EPS is the height (along the axial symmetry axis)
C             of the resulting cut sphere
C                Note that always EPS.LT.2*REV specified
C
C Warning:
C  In computations for spheres, use EPS=1.000001 instead of EPS=1.
C  EPS=1 can cause overflows in some rare cases.
C
C  LAM - the wavelength of incident light in the ambient.
C                  LAM=LAMBDA/SQRT(ZEPS0) here
C
C  RAT = 1 - particle size is specified in terms of the
C                equal-volume-sphere radius
C  RAT.ne.1 - particle size is specified in terms of the
C                equal-surface-area-sphere radius
C  AXI ... equivalent-(volume/surface-area)-sphere radius
C  REV=A=RAT*AXI ... equal-volume-sphere radius
C                  (feeded as REV to RSP* routines)
C  DDELT - required precision
C  XS  - Equiv. size parameter x=2*pi*rev*n_0/lambda
C
C  ALPHA and BETA - Euler angles (in degrees) specifying the
C          orientation of the scattering particle relative to
C          the laboratory reference frame (Refs. 6 and 7).
C
C  THET0 - zenith angle of the incident beam in degrees
C  THET - zenith angle of the scattered beam in degrees
C  PHI0 - azimuth angle of the incident beam in degrees
C  PHI - azimuth angle of the scattered beam in degrees
C--------/---------/---------/---------/---------/---------/---------/--
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER NOUT,NAXSM,ICHOICEV,ICHOICE
      LOGICAL YNCHECK

*      INCLUDE 'ampld.par.f'
      Parameter (NPN1=100, NPNG1=2600, NPNG2=2*NPNG1, NPN2=2*NPN1,
     &           NPL=NPN2+1, NPN3=NPN1+1,
     &           NPN4=NPN1, NPN5=2*NPN4, NPN6=NPN4+1)
* number of the output unit
      PARAMETER (NOUT=35)
*
      REAL*8  LAM,LAMBDA,MRR,MRI,RSNM,HT,DDELT,DDELTA,
     *        X(NPNG2),W(NPNG2),
     *        S(NPNG2),SS(NPNG2),AN(NPN1),R(NPNG2),DR(NPNG2),
     *        DDR(NPNG2),DRR(NPNG2),DRI(NPNG2),ANN(NPN1,NPN1)
      REAL*8 TR1(NPN2,NPN2),TI1(NPN2,NPN2)
c      REAL*8 XALPHA(300),XBETA(300),WALPHA(300),WBETA(300)
      REAL*4
     &     RT11(NPN6,NPN4,NPN4),RT12(NPN6,NPN4,NPN4),
     &     RT21(NPN6,NPN4,NPN4),RT22(NPN6,NPN4,NPN4),
     &     IT11(NPN6,NPN4,NPN4),IT12(NPN6,NPN4,NPN4),
     &     IT21(NPN6,NPN4,NPN4),IT22(NPN6,NPN4,NPN4)
      COMPLEX*16 S11,S12,S21,S22
      COMPLEX*16 zeps1,zeps0
*
      COMMON /CT/ TR1,TI1
* transfers the real and imaginary part of the T matrix (2*NMAX,2*NMAX)
* array for a given value of M from TMATR0 and TMATR to the AMPLDR
*
      COMMON /TMAT/ RT11,RT12,RT21,RT22,IT11,IT12,IT21,IT22
* transfers T matrix arrays obtained from TR1,TI1 in the AMPLDR
* to the AMPL routine
*
      COMMON /CHOICE/ ICHOICE
* transfers the choice of inversion to relevant matrix inversion
* routines
*
      COMMON /TOAMPLD/RAT,REV,ALPHA,BETA,DDELT
*
* transfers real*8 RAT,A(REV),ALPHA,BETA,DDELT from the main here
*
      COMMON /TOTAMPLD/THET0,THET,PHI0,PHI
*
* transfers real*8 THET0,THET,PHI0,PHI from the main here

      COMMON /TOIAMPLD/NCHECK,NAXSM,NDGS

* transfers integers NCHECK,NAXSM,NDGS from the main here
*
*****************************************************************
*
      P=DACOS(-1D0)                   !local PI constant
*
      ICHOICE=ICHOICEV
      A=REV
      LAM=LAMBDA/SQRT(ZEPS0)          !wavelength in the ambient

cc      write(6,*)'LAM,LAMBDA in AMPL=', LAM, LAMBDA
*
* the real part of the refractive index contrast
*
      MRR=DBLE(SQRT(ZEPS1/ZEPS0))
*
* the imaginary  part of the refractive index contrast
*
      MRI=DIMAG(SQRT(ZEPS1/ZEPS0))
*
      DDELTA=0.1D0*DDELT
*
* DDELT is used to test the accuracy of computing the
* optical cross sections. This accuracy is usually better
* than the absolute accuracy of computing the expansion coefficients
* of a normalized scattering matrix by a factor of 10. Therefore,
* the desired accuracy of computing the expansion coefficients
* is rescaled by a factor 0.1 before entering the test of the
* accuracy of computing the optical cross sections.

      IF (DABS(RAT-1D0).GT.1D-8.AND.NP.EQ.-1) CALL SAREA (EPS,RAT)
      IF (DABS(RAT-1D0).GT.1D-8.AND.NP.GE.0) CALL SURFCH(NP,EPS,RAT)
      IF (DABS(RAT-1D0).GT.1D-8.AND.NP.EQ.-2) CALL SAREAC (EPS,RAT)
      IF (NP.EQ.-3) CALL DROP (RAT)

      PRINT 7400, LAM,MRR,MRI

 7400 FORMAT('LAM=',F12.6,3X,'MRR=',D10.4,3X,'MRI=',D10.4)
*
*___________________________________________________
* Determination of the Wiscombe value of the floating
C angular momentum cutoff NMAX:

      XEV=2D0*P*A/LAM
      IXXX=XEV+4.05D0*XEV**0.333333D0     !Wiscombe conv. criterion for NMAX
      IF (XEV.GT.1.) THEN
      INM1=MAX0(3,IXXX)
      ELSE
      INM1=2                 !The Bessel package routine RYB
                             !requires NMAX to be at least 2
      END IF
*
      IF (INM1.GE.NPN1) PRINT 7333, NPN1
      IF (INM1.GE.NPN1) STOP
 7333 FORMAT('CONVERGENCE IS NOT OBTAINED FOR NPN1=',I3,
     &       '.  EXECUTION TERMINATED')

*_______________________________________________________________

      NGAUSS=NMAX*NDGS

      IF (YNCHECK) THEN

         write(6,*)
         write(6,*)'NMAX-convergence test'
         write(6,*)
         write(6,*)'(NGAUSS=NMAX*NDGS)'


c Internal determination of the floating angular momentum cutoff
c NMAX using convergence criterion of {Mis32}. It begins convergence
c convergence test with the Wiscombe value for the floating angular
c momentum cutoff NMAX with its subsequent increase by one, till
c the convergence criterion {Mis32} is satisfied
c
      QEXT1=0D0
      QSCA1=0D0

      DO 50 NMA=INM1,NPN1
         NMAX=NMA
         NGAUSS=NMAX*NDGS    !the number of the Gauss integration points

         IF (NGAUSS.GT.NPNG1) PRINT 7340, NGAUSS
         IF (NGAUSS.GT.NPNG1) STOP

 7340    FORMAT('NGAUSS =',I3,' I.E. IS GREATER THAN NPNG1.',
     &          '  EXECUTION TERMINATED')
c 7334    FORMAT(' NMAX =', I3,'  DC2=',D8.2,'   DC1=',D8.2)
*
         CALL CONST(NGAUSS,NMAX,X,W,AN,ANN,S,SS,NP,EPS,RSNM,HT)      !In AMPLDR
*
* specify particle shape:
         CALL VARY(LAM,MRR,MRI,A,EPS,RSNM,HT,NP,NGAUSS,X,P,
     &              PPI,PIR,PII,R,DR,DDR,DRR,DRI,NMAX)
*
* determine m=m'=0 elements of the T matrix
*
         CALL TMATR0 (NGAUSS,X,W,AN,ANN,PPI,PIR,PII,R,DR,
     &                 DDR,DRR,DRI,NMAX,NCHECK,NAXSM)
*
         QEXT=0D0
         QSCA=0D0
*
* make convergence test {Mis32} for a given NMAX:
*
         DO 4 N=1,NMAX
            N1=N+NMAX
            TR1NN=TR1(N,N)
            TI1NN=TI1(N,N)
            TR1NN1=TR1(N1,N1)
            TI1NN1=TI1(N1,N1)
            DN1=DFLOAT(2*N+1)
            QSCA=QSCA+DN1*(TR1NN*TR1NN+TI1NN*TI1NN
     &                    +TR1NN1*TR1NN1+TI1NN1*TI1NN1)
            QEXT=QEXT+(TR1NN+TR1NN1)*DN1
    4    CONTINUE

*>>> for debugging:
cc      OPEN(NOUT+1,FILE='tr1diag.dat')
cc      OPEN(NOUT+2,FILE='ti1diag.dat')
cc              DO N=1,2*NMAX
cc                  write(nout+1,*) TR1(N,N)
cc                  write(nout+2,*) TI1(N,N)
cc              enddo
cc      close(nout+1)
cc      close(nout+2)
*<<<
         write(6,*)'NMAX=',NMAX
         write(6,*)'NGAUSS=',NGAUSS
         write(6,*)'QSCA1=',QSCA1
         write(6,*)'QSCA=',QSCA
         write(6,*)'QEXT1=',QEXT1
         write(6,*)'QEXT=',QEXT
*<<<
         DSCA=DABS((QSCA1-QSCA)/QSCA)
         DEXT=DABS((QEXT1-QEXT)/QEXT)
         QEXT1=QEXT
         QSCA1=QSCA

c        PRINT 7334, NMAX,DSCA,DEXT

         IF(DSCA.LE.DDELTA.AND.DEXT.LE.DDELTA) GO TO 55
         IF (NMA.EQ.NPN1) PRINT 7333, NPN1
         IF (NMA.EQ.NPN1) STOP

   50 CONTINUE                   !Successful L-convergence test exit

* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   55 CONTINUE                   !Begin NGAUSS-convergence test

         write(6,*)
         write(6,*)'NGAUSS-convergence test'
         write(6,*)

      NNNGGG=NGAUSS+1

      IF (NGAUSS.EQ.NPNG1) PRINT 7336
 7336    FORMAT('WARNING: NGAUSS=NPNG1')

      IF (NGAUSS.EQ.NPNG1) GO TO 160

      DO 150 NGAUS=NNNGGG,NPNG1
*
         IF (NGAUS.EQ.NPNG1) PRINT 7336
*
         NGAUSS=NGAUS
cc         NGGG=2*NGAUSS
*
* GIF division points and weights + other numerical constants
*
         CALL CONST(NGAUSS,NMAX,X,W,AN,ANN,S,SS,NP,EPS,RSNM,HT)     !In AMPLDR
*
* specify particle shape:
*
         CALL VARY(LAM,MRR,MRI,A,EPS,RSNM,HT,NP,NGAUSS,X,P,
     &              PPI,PIR,PII,R,DR,DDR,DRR,DRI,NMAX)
*
* determine m=m'=0 elements of the T matrix
*
         CALL TMATR0 (NGAUSS,X,W,AN,ANN,PPI,PIR,PII,R,DR,
     &                 DDR,DRR,DRI,NMAX,NCHECK,NAXSM)
*
         QEXT=0D0
         QSCA=0D0

         DO 104 N=1,NMAX
            N1=N+NMAX
            TR1NN=TR1(N,N)
            TI1NN=TI1(N,N)
            TR1NN1=TR1(N1,N1)
            TI1NN1=TI1(N1,N1)

            DN1=DFLOAT(2*N+1)

            QSCA=QSCA+DN1*(TR1NN*TR1NN+TI1NN*TI1NN
     &                    +TR1NN1*TR1NN1+TI1NN1*TI1NN1)
            QEXT=QEXT+(TR1NN+TR1NN1)*DN1

  104    CONTINUE

         DSCA=DABS((QSCA1-QSCA)/QSCA)
         DEXT=DABS((QEXT1-QEXT)/QEXT)

c        PRINT 7337, NGGG,DSCA,DEXT
c 7337    FORMAT(' NG=',I3,'  DC2=',D8.2,'   DC1=',D8.2)
*<<<
         write(6,*)'NGAUSS=',NGAUSS
         write(6,*)'QSCA1=',QSCA1
         write(6,*)'QSCA=',QSCA
         write(6,*)'QEXT1=',QEXT1
         write(6,*)'QEXT=',QEXT

         IF(DSCA.LE.DDELTA.AND.DEXT.LE.DDELTA) GO TO 160
*<<<
         QEXT1=QEXT
         QSCA1=QSCA
*
  150 CONTINUE

* %%%%%%%%%%%%%%%%%% Successful NGAUSS-convergence test %%%%%%%%%%%%%%%

      ELSE  IF (.NOT.YNCHECK) THEN

* GIF division points and weights + other numerical constants
*
         CALL CONST(NGAUSS,NMAX,X,W,AN,ANN,S,SS,NP,EPS,RSNM,HT)     !In AMPLDR
*
* specify particle shape:
*
         CALL VARY(LAM,MRR,MRI,A,EPS,RSNM,HT,NP,NGAUSS,X,P,
     &              PPI,PIR,PII,R,DR,DDR,DRR,DRI,NMAX)
*
* determine m=m'=0 elements of the T matrix
*
         CALL TMATR0 (NGAUSS,X,W,AN,ANN,PPI,PIR,PII,R,DR,
     &                 DDR,DRR,DRI,NMAX,NCHECK,NAXSM)
*
*<<<
*
      END IF                     !YNCHECK

 160  CONTINUE
*<<<
      WRITE(6,*)
      WRITE(6,*)'NMAX=',NMAX
      WRITE(6,*)'NGAUSS=',NGAUSS
      WRITE(6,*)
cc      WRITE(NOUT,*)'NMAX=',NMAX
cc      WRITE(NOUT,*)'NGAUSS=',NGAUSS

*<<<

*************   Calculation of scattering cross sections   *********

*Initialization:

      QSCA=0D0
      QEXT=0D0
      NNM=2*NMAX

*   >>>  DETERMINATION OF QEXT AND QSCA CONTRIBUTIONS FOR M=0

      DO 204 N=1,NNM

         QEXT=QEXT+TR1(N,N)

cc         if ((n.le.5).or.(((n-nmax).le.5).and.((n-nmax).gt.0))) then
cc         xx=-dble(TR1(N,N))       ! sin^2\eta_l
cc         write(nout+15,*) 'n, sin^2\eta_l', n, xx
cc         end if

  204 CONTINUE


* Given RT1 and IT1 matrices from TMATR0 routine,
* assigning of RT^{ij} and IT^{ij} matrix entries to be
* used later by AMPL routine

      DO 213 N2=1,NMAX
         NN2=N2+NMAX
         DO 213 N1=1,NMAX
            NN1=N1+NMAX
            ZZ1=TR1(N1,N2)
            RT11(1,N1,N2)=ZZ1
            ZZ2=TI1(N1,N2)
            IT11(1,N1,N2)=ZZ2
            ZZ3=TR1(N1,NN2)
            RT12(1,N1,N2)=ZZ3
            ZZ4=TI1(N1,NN2)
            IT12(1,N1,N2)=ZZ4
            ZZ5=TR1(NN1,N2)
            RT21(1,N1,N2)=ZZ5
            ZZ6=TI1(NN1,N2)
            IT21(1,N1,N2)=ZZ6
            ZZ7=TR1(NN1,NN2)
            RT22(1,N1,N2)=ZZ7
            ZZ8=TI1(NN1,NN2)
            IT22(1,N1,N2)=ZZ8
*
            QSCA=QSCA+ZZ1*ZZ1+ZZ2*ZZ2+ZZ3*ZZ3+ZZ4*ZZ4
     &           +ZZ5*ZZ5+ZZ6*ZZ6+ZZ7*ZZ7+ZZ8*ZZ8
*
  213 CONTINUE   !end of the loop over orbital numbers
*________________


*<<<

      if (abs(qsca).gt.(1.0001d0*abs(qext))) then
         write(6,*)'M=',0
         write(6,*)'QSCA=',QSCA
         write(6,*)'QEXT=',QEXT
       write(6,*)
     & 'WARNING: abs(qsca).gt.abs(qext)!!!'
c     pause
      end if

      if (qext.gt.1.d-7) then
       write(6,*)
     & 'WARNING: Partial sum QEXT has to be negative!'
c      pause
      end if


cc         write(nout,*)'M=',M
cc         write(nout,*)'QSCA=',QSCA
cc         write(nout,*)'QSC=',QSC
cc         write(nout,*)'QEXT=',QEXT
cc         write(nout,*)'QXT=',QXT
*<<<

*   >>>  DETERMINATION OF QEXT AND QSCA CONTRIBUTIONS FOR M >< 0

      DO 220 M=1,NMAX
*
c         CALL TMATR(M,NGAUSS,X,W,AN,ANN,S,SS,PPI,PIR,PII,R,DR,
c     &               DDR,DRR,DRI,NMAX,NCHECK,NAXSM)

         CALL TMTR(M,NGAUSS,X,W,AN,ANN,PPI,PIR,PII,R,DR,DDR,DRR,
     &               DRI,NMAX,NCHECK,NAXSM)
*
* <<< returns  m=m'>0 elements of the T matrix
*
         NM=NMAX-M+1
         M1=M+1
         QSC=0D0

* Given RT1 and IT1 matrices from TMATR routine,
* assigning of RT^{ij} and IT^{ij} matrix entries to be
* used later by AMPL routine.
*
         DO 214 N2=1,NM              !summation over orbital numbers

* conversion of the N22 index of RT1 and IT1 matrices
* to the index NN2 of RT^{ij} and IT^{ij} matrices

            NN2=N2+M-1        !from M to NMAX
            N22=N2+NM         !from NMAX+1 to 2*NMAX-M+1

            DO 214 N1=1,NM           !summation over orbital numbers

* conversion of the N11 index of RT1 and IT1 matrices
* to the index NN1 of RT^{ij} and IT^{ij} matrices

               NN1=N1+M-1        !from M to NMAX
               N11=N1+NM         !from NMAX+1 to 2*NMAX-M+1

               ZZ1=TR1(N1,N2)
               RT11(M1,NN1,NN2)=ZZ1
               ZZ2=TI1(N1,N2)
               IT11(M1,NN1,NN2)=ZZ2
               ZZ3=TR1(N1,N22)
               RT12(M1,NN1,NN2)=ZZ3
               ZZ4=TI1(N1,N22)
               IT12(M1,NN1,NN2)=ZZ4
               ZZ5=TR1(N11,N2)
               RT21(M1,NN1,NN2)=ZZ5
               ZZ6=TI1(N11,N2)
               IT21(M1,NN1,NN2)=ZZ6
               ZZ7=TR1(N11,N22)
               RT22(M1,NN1,NN2)=ZZ7
               ZZ8=TI1(N11,N22)
               IT22(M1,NN1,NN2)=ZZ8
*
               QSC=QSC+(ZZ1*ZZ1+ZZ2*ZZ2+ZZ3*ZZ3+ZZ4*ZZ4
     &                 +ZZ5*ZZ5+ZZ6*ZZ6+ZZ7*ZZ7+ZZ8*ZZ8)*2D0
*
* multiplication by 2d0 here accounts for +/-M symmetry of resulting
* expressions

  214    CONTINUE     !end of the loop over orbital numbers

         NNM=2*NM
         QXT=0D0

         DO 215 N=1,NNM

            QXT=QXT+TR1(N,N)*2D0       !multiplication by 2d0 accounts
                                       !for +/-M symmetry of resulting
                                       !expressions
  215    CONTINUE

*<<<
* Summation over magnetic quantum number:

         QSCA=QSCA+QSC
         QEXT=QEXT+QXT

*<<<

      if (abs(qsc).gt.(1.0001d0*abs(qxt))) then
         write(6,*)'M=',M
         write(6,*)'QSCA=',QSCA
         write(6,*)'QSC=',QSC
         write(6,*)'QEXT=',QEXT
         write(6,*)'QXT=',QXT
       write(6,*)
     & 'WARNING: abs(qsc).gt.abs(qxt)!!!'
c     pause
      end if

      if (qxt.ge.1d-7) then
       write(6,*)
     & 'WARNING: Partial sum QXT has to be negative!'
C      pause
      end if

cc         write(nout,*)'M=',M
cc         write(nout,*)'QSCA=',QSCA
cc         write(nout,*)'QSC=',QSC
cc         write(nout,*)'QEXT=',QEXT
cc         write(nout,*)'QXT=',QXT
*
c        PRINT 7800,M,DABS(QXT),QSC,NMAX
c 7800    FORMAT(' m=',I3,'  qxt=',D12.6,'  qsc=',D12.6,
c     &          '  nmax=',I3)

  220 CONTINUE    !end of loop over m's

*<<<
         write(6,*)'QSCA=',QSCA
         write(6,*)'QEXT=',QEXT

*
* 'QSCA' and '-QEXT' are now 'efficiency factors' for scattering
* and extinction (=\sum_{AL} \sin^2\eta_{AL}).


      QABS=-QEXT-QSCA       !absorption
      WALB=-QSCA/QEXT       !albedo

      IF (ABS(WALB).GT.1D0+DDELTA) THEN
      PRINT 9111
 9111 FORMAT ('WARNING: THE ALBEDO WALB IS GREATER THAN 1')
      WRITE(6,*)'WALB=',WALB
      END IF

*<<<
C In order to convert the efficiencies 'QSCA' and '-QEXT' into
C normalized (per scatterer effective surface S=pi*rev**2)
C cross-sections
C         QEXT=(2/x**2) \sum_{AL} \sin^2\eta_{AL}  for eta real
C         QEXT=(2/x**2) \sum_{AL} Re (T)    for a general eta
C         QSCA=(2/x**2) \sum_{AL} |T|^2
C (cf Eq.(2.135-8) of Newton�s book)
C At the moment, the prefactor (2/xev**2) is still missing.
C (lambda here is the wavelength in the exterior ambient medium)
cc      write(6,*)'LAM in AMPL=', LAM
c         FAC=LAM**2/(2.d0*P**2*REV**2)     !=2/xs**2
         FAC=2.d0/XEV**2
         write(nout,*)    lambda, FAC*QSCA
         write(nout+1,*)  lambda,-FAC*QEXT
         write(nout+2,*)  lambda, FAC*QABS
         write(nout+5,*)  lambda, FAC*WALB
         write(nout+10,*)
         write(nout+10,*) lambda
         write(nout+12,*)
cc         write(nout+12,*) lambda
cc         write(nout+13,*) lambda
cc         write(nout+16,*) -qext

*<<<
*_________________________________________________________
C  COMPUTATION OF THE AMPLITUDE AND PHASE MATRICES
C  AMPLITUDE MATRIX [Eqs. (2)-(4) of Ref. 6]
*
      CALL AMPL (NMAX,LAM,THET0,THET,PHI0,PHI,ALPHA,BETA,
     &           S11,S12,S21,S22)
*
C  PHASE MATRIX [Eqs. (13)-(29) of Ref. 6]
      Z11=0.5D0*(S11*DCONJG(S11)+S12*DCONJG(S12)
     &          +S21*DCONJG(S21)+S22*DCONJG(S22))
      Z12=0.5D0*(S11*DCONJG(S11)-S12*DCONJG(S12)
     &          +S21*DCONJG(S21)-S22*DCONJG(S22))
      Z13=-S11*DCONJG(S12)-S22*DCONJG(S21)
      Z14=(0D0,1D0)*(S11*DCONJG(S12)-S22*DCONJG(S21))
      Z21=0.5D0*(S11*DCONJG(S11)+S12*DCONJG(S12)
     &          -S21*DCONJG(S21)-S22*DCONJG(S22))
      Z22=0.5D0*(S11*DCONJG(S11)-S12*DCONJG(S12)
     &          -S21*DCONJG(S21)+S22*DCONJG(S22))
      Z23=-S11*DCONJG(S12)+S22*DCONJG(S21)
      Z24=(0D0,1D0)*(S11*DCONJG(S12)+S22*DCONJG(S21))
      Z31=-S11*DCONJG(S21)-S22*DCONJG(S12)
      Z32=-S11*DCONJG(S21)+S22*DCONJG(S12)
      Z33=S11*DCONJG(S22)+S12*DCONJG(S21)
      Z34=(0D0,-1D0)*(S11*DCONJG(S22)+S21*DCONJG(S12))
      Z41=(0D0,1D0)*(S21*DCONJG(S11)+S22*DCONJG(S12))
      Z42=(0D0,1D0)*(S21*DCONJG(S11)-S22*DCONJG(S12))
      Z43=(0D0,-1D0)*(S22*DCONJG(S11)-S12*DCONJG(S21))
      Z44=S22*DCONJG(S11)-S12*DCONJG(S21)


      WRITE(NOUT+10,5001) Z11,Z12,Z13,Z14
      WRITE(NOUT+10,5001) Z21,Z22,Z23,Z24
      WRITE(NOUT+10,5001) Z31,Z32,Z33,Z34
      WRITE(NOUT+10,5001) Z41,Z42,Z43,Z44

 5001 FORMAT (4F10.4)

c      ITIME=MCLOCK()
c      TIME=DFLOAT(ITIME)/6000D0
c      PRINT 1001,TIME
c 1001 FORMAT (' time =',F8.2,' min')

      RETURN
      END

C********************************************************************

      SUBROUTINE AMPL (NMAX,DLAM,TL,TL1,PL,PL1,ALPHA,BETA,
     &                 VV,VH,HV,HH)
C--------/---------/---------/---------/---------/---------/---------/--
C >>> NMAX,DLAM,TL,TL1,PL,PL1,ALPHA,BETA
C <<< VV,VH,HV,HH
C=================
C
C    GIVEN T MATRIX IN COMMON BLOCK IT CALCULATES THE AMPLITUDE MATRIX
C
C    This routine closely follows exposition by
C       M. I. Mishchenko, Calculation of the amplitude matrix
C       for a nonspherical particle in a fixed orientation,
C       Appl. Opt. vol. 39, 1026-1031 (2000).
C
C   NMAX - angular momentum cutoff
C   DLAM=LAMBDA/SQRT(ZEPS0)  - wavelength of incident light in the ambient
C                     (vacuum wavelength divided by SQRT(ZEPS0))
C
C   LAMBDA - vacuum wavelength. Determined as DLAM*SQRT(ZEPS0) and
C            only used here for the write out purposes
C   TL,TL1,PL,PL1 ... angles in degrees
C                     determined w.r.t laboratory frame:
C   TL (THET0 IN MAIN) - zenith angle of the incident beam in degrees
C   TL1 (THET IN MAIN) - zenith angle of the scattered beam in degrees
C   PL (PHI0 IN MAIN) - azimuth angle of the incident beam in degrees
C   PL1 (PHI IN MAIN) - azimuth angle of the scattered beam in degrees
C
C   ALPHA and BETA - Euler angles (in degrees) specifying the
C         orientation of the scattering particle relative to the
C         laboratory reference frame (Refs. 6 and 7).
C   VV,VH,HV,HH ... amplitude scattering matrix elements S11,S12,S21,S22
C--------/---------/---------/---------/---------/---------/---------/--
      IMPLICIT REAL*8 (A-B,D-H,O-Z)
      IMPLICIT COMPLEX*16 (C)
      INTEGER NOUT

* number of the output unit
      PARAMETER (NOUT=35)
*     INCLUDE 'ampld.par.f'
      Parameter (NPN1=100, NPNG1=2600, NPNG2=2*NPNG1, NPN2=2*NPN1,
     &           NPL=NPN2+1, NPN3=NPN1+1,
     &           NPN4=NPN1, NPN5=2*NPN4, NPN6=NPN4+1)
      REAL*8 DLAM,LAMBDA,CEXT,CEXT1,CEXT2
      REAL*8 AL(3,2),AL1(3,2),AP(2,3),AP1(2,3),B(3,3),
     *       R(2,2),R1(2,2),C(3,2),CA,CB,CT,CP,CTP,CPP,CT1,CP1,
     *       CTP1,CPP1
      REAL*8 DV1(NPN6),DV2(NPN6),DV01(NPN6),DV02(NPN6)
      REAL*4
     &     TR11(NPN6,NPN4,NPN4),TR12(NPN6,NPN4,NPN4),
     &     TR21(NPN6,NPN4,NPN4),TR22(NPN6,NPN4,NPN4),
     &     TI11(NPN6,NPN4,NPN4),TI12(NPN6,NPN4,NPN4),
     &     TI21(NPN6,NPN4,NPN4),TI22(NPN6,NPN4,NPN4)
      COMPLEX*16 CAL(NPN4,NPN4),VV,VH,HV,HH,ZEPS0
*_____
      COMMON /TMAT/ TR11,TR12,TR21,TR22,TI11,TI12,TI21,TI22
*_
      COMMON /DIELF/ zeps0
      COMMON /REVF/ rev
      COMMON /CYLPAR/ rsnm, hlength

*
* transfers ZEPS0,REV here from the main
*____
*
C Checking the initial set of angles TL,TL1,PL,PL1,ALPHA,BETA
C for allowability

      IF (ALPHA.LT.0D0.OR.ALPHA.GT.360D0.OR.
     &    BETA.LT.0D0.OR.BETA.GT.180D0.OR.
     &    TL.LT.0D0.OR.TL.GT.180D0.OR.
     &    TL1.LT.0D0.OR.TL1.GT.180D0.OR.
     &    PL.LT.0D0.OR.PL.GT.360D0.OR.
     &    PL1.LT.0D0.OR.PL1.GT.360D0) THEN
          WRITE(NOUT,2000)
          STOP
      ELSE
          CONTINUE
      ENDIF
 2000 FORMAT ('AN ANGULAR PARAMETER IS OUTSIDE ITS',
     &        ' ALLOWABLE RANGE')

* SPECIFYING NUMERICAL CONSTANTS:

      PIN=DACOS(-1D0)         !=PI
      PIN2=PIN*0.5D0          !=PI/2
      PI=PIN/180D0            !=PI/180

* conversion from degrees to radians:
      ALPH=ALPHA*PI
      BET=BETA*PI
      THETL=TL*PI
      PHIL=PL*PI
      THETL1=TL1*PI
      PHIL1=PL1*PI

* initialization of the vacuum wavelength LAMBDA

      LAMBDA=DLAM*SQRT(ZEPS0)         !vacuum wavelength

      EPS=1D-7
      IF (THETL.LT.PIN2) THETL=THETL+EPS
      IF (THETL.GT.PIN2) THETL=THETL-EPS
      IF (THETL1.LT.PIN2) THETL1=THETL1+EPS
      IF (THETL1.GT.PIN2) THETL1=THETL1-EPS
      IF (PHIL.LT.PIN) PHIL=PHIL+EPS
      IF (PHIL.GT.PIN) PHIL=PHIL-EPS
      IF (PHIL1.LT.PIN) PHIL1=PHIL1+EPS
      IF (PHIL1.GT.PIN) PHIL1=PHIL1-EPS
      IF (BET.LE.PIN2.AND.PIN2-BET.LE.EPS) BET=BET-EPS
      IF (BET.GT.PIN2.AND.BET-PIN2.LE.EPS) BET=BET+EPS

C   Given TL,TL1,PL,PL1 in laboratory frame
C   COMPUTE THETP, PHIP, THETP1, AND PHIP1 in particle frame
C   (see EQS. (9), (20), AND (21))
C
C incident beam:

      CB=DCOS(BET)
      SB=DSIN(BET)
      CT=DCOS(THETL)
      ST=DSIN(THETL)
      CP=DCOS(PHIL-ALPH)
      SP=DSIN(PHIL-ALPH)

      CTP=CT*CB+ST*SB*CP             !Eq. (9)
      THETP=DACOS(CTP)
      CPP=CB*ST*CP-SB*CT             !Eq. (20)
      SPP=ST*SP                      !Eq. (21)
      PHIP=DATAN(SPP/CPP)

      IF (PHIP.GT.0D0.AND.SP.LT.0D0) PHIP=PHIP+PIN
      IF (PHIP.LT.0D0.AND.SP.GT.0D0) PHIP=PHIP+PIN
      IF (PHIP.LT.0D0) PHIP=PHIP+2D0*PIN

C scattered beam:

      CT1=DCOS(THETL1)
      ST1=DSIN(THETL1)
      CP1=DCOS(PHIL1-ALPH)
      SP1=DSIN(PHIL1-ALPH)

      CTP1=CT1*CB+ST1*SB*CP1          !Eq. (9)
      THETP1=DACOS(CTP1)
      CPP1=CB*ST1*CP1-SB*CT1          !Eq. (20)
      SPP1=ST1*SP1                    !Eq. (21)
      PHIP1=DATAN(SPP1/CPP1)

      IF (PHIP1.GT.0D0.AND.SP1.LT.0D0) PHIP1=PHIP1+PIN
      IF (PHIP1.LT.0D0.AND.SP1.GT.0D0) PHIP1=PHIP1+PIN
      IF (PHIP1.LT.0D0) PHIP1=PHIP1+2D0*PIN

C____________COMPUTE MATRIX BETA, EQ. (22) of {Mis39}

      CA=DCOS(ALPH)
      SA=DSIN(ALPH)
      B(1,1)=CA*CB
      B(1,2)=SA*CB
      B(1,3)=-SB
      B(2,1)=-SA
      B(2,2)=CA
      B(2,3)=0D0
      B(3,1)=CA*SB
      B(3,2)=SA*SB
      B(3,3)=CB

C____________COMPUTE 3x2 MATRICES AL AND AL1 for incident and
C              scattered beams in laboratory frame
C                   [Eq. (15)  of {Mis39}]

      CP=DCOS(PHIL)
      SP=DSIN(PHIL)
      CP1=DCOS(PHIL1)
      SP1=DSIN(PHIL1)

C incident beam:
      AL(1,1)=CT*CP
      AL(1,2)=-SP
      AL(2,1)=CT*SP
      AL(2,2)=CP
      AL(3,1)=-ST
      AL(3,2)=0D0

C scattered beam:
      AL1(1,1)=CT1*CP1
      AL1(1,2)=-SP1
      AL1(2,1)=CT1*SP1
      AL1(2,2)=CP1
      AL1(3,1)=-ST1
      AL1(3,2)=0D0

C____________COMPUTE 2X3 MATRICES AP^(-1) AND AP1^(-1) for incident
C             and scattered beams in particle frame
C                   [Eq. (16)  of {Mis39}]
      CT=CTP
      ST=DSIN(THETP)
      CP=DCOS(PHIP)
      SP=DSIN(PHIP)
      CT1=CTP1
      ST1=DSIN(THETP1)
      CP1=DCOS(PHIP1)
      SP1=DSIN(PHIP1)
C incident beam:
      AP(1,1)=CT*CP
      AP(1,2)=CT*SP
      AP(1,3)=-ST
      AP(2,1)=-SP
      AP(2,2)=CP
      AP(2,3)=0D0
C scattered beam:
      AP1(1,1)=CT1*CP1
      AP1(1,2)=CT1*SP1
      AP1(1,3)=-ST1
      AP1(2,1)=-SP1
      AP1(2,2)=CP1
      AP1(2,3)=0D0

C____________COMPUTE MATRICES R AND R^(-1), EQ. (14)

C Computation of R for incident beam
      DO I=1,3
         DO J=1,2
            X=0D0
            DO K=1,3
               X=X+B(I,K)*AL(K,J)
            ENDDO
            C(I,J)=X
         ENDDO
      ENDDO
      DO I=1,2
         DO J=1,2
            X=0D0
            DO K=1,3
               X=X+AP(I,K)*C(K,J)
            ENDDO
            R(I,J)=X
         ENDDO
      ENDDO

C Computation of R^(-1) for scattered beam using Cramers rule:

      DO I=1,3
         DO J=1,2
            X=0D0
            DO K=1,3
               X=X+B(I,K)*AL1(K,J)
            ENDDO
            C(I,J)=X
         ENDDO
      ENDDO
      DO I=1,2
         DO J=1,2
            X=0D0
            DO K=1,3
               X=X+AP1(I,K)*C(K,J)
            ENDDO
            R1(I,J)=X
         ENDDO
      ENDDO
C====
C R for scattered beam determined, now Cramers rule:

      D=1D0/(R1(1,1)*R1(2,2)-R1(1,2)*R1(2,1))
      X=R1(1,1)
      R1(1,1)=R1(2,2)*D
      R1(1,2)=-R1(1,2)*D
      R1(2,1)=-R1(2,1)*D
      R1(2,2)=X*D

C____________MATRICES R AND R^(-1) determined
C=========================================================
C      THE AMPLITUDE MATRIX IN PARTICLE FRAME
C
      CI=(0D0,1D0)

C >>> ALPHA numerical prefactors without phi-angles
C     (following Eq. (28))

      DO 5 NN=1,NMAX
         DO 5 N=1,NMAX
            CN=CI**(NN-N-1)
            DNN=DFLOAT((2*N+1)*(2*NN+1))
            DNN=DNN/DFLOAT( N*NN*(N+1)*(NN+1) )
            RN=DSQRT(DNN)
            CAL(N,NN)=CN*RN
    5 CONTINUE

      DCTH0=CTP             !\cos\vartheta_{inc}^P
      DCTH=CTP1             !\cos\vartheta_{sca}^P
      PH=PHIP1-PHIP         !(\varphi_{sca}^P-\varphi_{inc}^P)

* amplitude scattering matrix elements S11,S12,S21,S22 initialization

      VV=(0D0,0D0)
      VH=(0D0,0D0)
      HV=(0D0,0D0)
      HH=(0D0,0D0)
C______________________________________________________________
C Main summation loop:

      DO 500 M=0,NMAX
         M1=M+1
         NMIN=MAX(M,1)
*
* Specify pi- and tau- scattering functions:

         CALL VIGAMPL (DCTH, NMAX, M, DV1, DV2)
         CALL VIGAMPL (DCTH0, NMAX, M, DV01, DV02)
*
         FC=2D0*DCOS(M*PH)    !takes into account +/- m contribution
         FS=2D0*DSIN(M*PH)
*
         DO 400 NN=NMIN,NMAX

            DV1NN=DV01(NN)           !\pi-functions
            DV2NN=DV02(NN)           !\tau-functions

            DO 400 N=NMIN,NMAX
               DV1N=DV1(N)           !\pi-functions
               DV2N=DV2(N)           !\tau-functions

               CT11=DCMPLX(TR11(M1,N,NN),TI11(M1,N,NN))
               CT22=DCMPLX(TR22(M1,N,NN),TI22(M1,N,NN))

               IF (M.EQ.0) THEN     !T^{21}=T^{12}=0 in particle frame

                  CN=CAL(N,NN)*DV2N*DV2NN

                  VV=VV+CN*CT22
                  HH=HH+CN*CT11

                 ELSE   !T^{21}\neq T^{12}\neq 0

                  CT12=DCMPLX(TR12(M1,N,NN),TI12(M1,N,NN))
                  CT21=DCMPLX(TR21(M1,N,NN),TI21(M1,N,NN))

* complete \alpha-factors (Eq. (28)) taking
* into account w.r.t. summation over +/- m in particle frame:
*
*     T^{11}_{-mnn'} = T^{11}_{mnn'}; T^{22}_{-mnn'} = T^{22}_{mnn'}
*  T^{12}_{-mnn'} = - T^{12}_{mnn'}; T^{21}_{-mnn'} = - T^{21}_{mnn'}

                  CN1=CAL(N,NN)*FC
                  CN2=CAL(N,NN)*FS

                  D11=DV1N*DV1NN    !\pi-\pi
                  D12=DV1N*DV2NN    !\pi-\tau
                  D21=DV2N*DV1NN    !\tau-\pi
                  D22=DV2N*DV2NN    !\tau-\tau

                  VV=VV+(CT11*D11+CT21*D21
     &                  +CT12*D12+CT22*D22)*CN1

                  VH=VH+(CT11*D12+CT21*D22
     &                  +CT12*D11+CT22*D21)*CN2

                  HV=HV-(CT11*D21+CT21*D11
     &                  +CT12*D22+CT22*D12)*CN2

                  HH=HH+(CT11*D22+CT21*D12
     &                  +CT12*D21+CT22*D11)*CN1
               ENDIF

  400    CONTINUE      !(over n,n')
  500 CONTINUE         !end of main summation loop (over m)

C Final multiplication of S11,S12,S21,S22 by (1/k) in the
C original code:

      DK=2D0*PIN/DLAM     !wavevector in surrounding medium
      VV=VV/DK
      VH=VH/DK
      HV=HV/DK
      HH=HH/DK

C   amplitude scattering matrix elements S11,S12,S21,S22 determined
C==================================================================
C TRANSFORMATION OF THE AMPLITUDE MATRIX FROM PARTICLE TO
c LABORATORY FRAME:

      CVV=VV*R(1,1)+VH*R(2,1)
      CVH=VV*R(1,2)+VH*R(2,2)
      CHV=HV*R(1,1)+HH*R(2,1)
      CHH=HV*R(1,2)+HH*R(2,2)
      VV=R1(1,1)*CVV+R1(1,2)*CHV
      VH=R1(1,1)*CVH+R1(1,2)*CHH
      HV=R1(2,1)*CVV+R1(2,2)*CHV
      HH=R1(2,1)*CVH+R1(2,2)*CHH

      PRINT 1101, VV
      PRINT 1102, VH
      PRINT 1103, HV
      PRINT 1104, HH

* For particles with plane of symmetry:
C If THET0=THET=90, then the incidence is perpendicular to
C the axis of axial symmetry  ===> E_\theta is along the axis
C of axial symmetry, whereas E_\phi is perpendicular to
C the symmetry axis. Then
C      C_{ext} = (4\pi/k_1} \mb{Im} S_{11}    if E_\phi=0
C
C or
C
C      C_{ext} = (4\pi/k_1} \mb{Im} S_{22}    if E_\theta=0
C
C For particles with plane of symmetry:
C Extiction for E along the axis of axial symmetry:
      cext1=2.d0*dlam*dimag(VV)      !Eq. (2.159)

C Extiction for E perpendicular to the axis of axial symmetry:
      cext2=2.d0*dlam*dimag(HH)      !Eq. (2.159)

C Orientationally averaged extiction
      cext=dlam*dimag(VV+HH)         !Eq. (5.97)

      write(6,*)'C_{ext}=\fr{2\pi}{k_1} Im (S_{11}+S_{22})=',
     & cext               !=2.d0*PIN*dimag(VV+HH)/k_1

      write(10, *) lambda, cext1, cext2
cc      FAC=lambda**2/(2.d0*PIN**2*REV**2)     !=2/xs**2
      FAC=1.d0/(PIN*REV**2)   !an effective geom. cross section

      write(nout+3,1107)lambda,fac*cext,fac*cext1,fac*cext2
      write(nout+12,1105) lambda, VV,VH
      write(nout+12,1105) lambda, HV,HH
      write(nout+13,1106) lambda, dble(VV*dconjg(vv) + VH*dconjg(vh)),
     &                               dble((vv+vh)*dconjg(vv+vh))
      write(nout+14,1106) lambda, dble(hV*dconjg(hv) + hH*dconjg(hh)),
     &                               dble((hv+hh)*dconjg(hv+hh))


 1101 FORMAT ('S11=',D11.5,' + i*',D11.5)
 1102 FORMAT ('S12=',D11.5,' + i*',D11.5)
 1103 FORMAT ('S21=',D11.5,' + i*',D11.5)
 1104 FORMAT ('S22=',D11.5,' + i*',D11.5)
 1105 FORMAT (F8.2,5X,D11.5,2X,D11.5,5X,D11.5,2X,D11.5)
 1106 FORMAT (F8.2,5X,D11.5,5X,D11.5)
 1107 FORMAT (F10.4,3(5X,D16.8))

      RETURN
      END

      SUBROUTINE VIGAMPL (X, NMAX, M, DDV1, DV2)
C--------/---------/---------/---------/---------/---------/---------/--
C >>> X,NMAX,M
C <<< DDV1, DV2   (related to scattering \pi and \tau functions)
C =============
C     For a given azimuthal number M.GE.0 returns
C
C  \pi scattering function:
C     DDV1(N)=dvig(0,m,n,arccos x)/sin(arccos x) = m*d_{0m}^{(l)}/ sin\theta
C
C  \tau scattering function:
C     DV2(N)=[d/d(arccos x)] dvig(0,m,n,arccos x) = d d_{0m}^{(l)}/d\theta
C
C     for 1.LE.N.LE.LMAX and 0.LE.X.LE.1
C      DDV1 is calculated because (DV1/sin\theta) is singular for
C             either \beta=0 or \beta=\pi
C     (For a given M.NEQ.0, only the M.LE.N.LE.LMAX terms are determined!)
C
C     According to Eq. (4.1.24) of Ref. \ct{Ed}:
C
C             d_{00}^{(l)}(\theta)= P_l(\cos\theta)
C
C     (Rodrigues formula [Eq. (2.5.14) of Ref. \ct{Ed}] then yields
C                       P_1(x)=x; P_2=(3x^2-1)/2; etc.
C
C     Similar to routine VIG, which however returns
C     DV1(N)=dvig(0,m,n,arccos x)   ! = d_{0m}^{(l)}
C
C     In addition, VIGAMPL(V) has a block treating the case when
C     arccos x is very small option
C
C     (There is a missing $l$ factor in the 2nd term in the curly bracket
C     in recurrence (35) of Ref. \ct{Mis39} for DV2).
C
C     X=cos(theta), where theta is the polar angle
C     LMAXD ... maximal angular momentum cutoff
C     NMAX ... floating  angular momentum cutoff
C--------/---------/---------/---------/---------/---------/---------/--
*      INCLUDE 'ampld.par.f'
      Parameter (NPN1=100, NPNG1=2600, NPNG2=2*NPNG1, NPN2=2*NPN1,
     &           NPL=NPN2+1, NPN3=NPN1+1,
     &           NPN4=NPN1, NPN5=2*NPN4, NPN6=NPN4+1)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 DDV1(NPN6), DV2(NPN6)

      integer n,nmax,M,I,I2
      REAL*8 A,X,QS,QS1,DSI,D1,D2,D3,DER,DN,DX,QN,QN1,QN2,
     & QNM,QNM1,QMM

* DDV1 and DV2 initialization
      DO 1 N=1,NMAX
          DDV1(N)=0.D0
         DV2(N) =0.D0
    1 CONTINUE

      DX=DABS(X)
      A=1.D0
      QS=DSQRT(1D0-X*X)                        !sin\theta
*
*********************************************************
C   For theta=0 [see Eqs. above]:
C              d_{00}^{(0)}(0)=1
C              d_{01}^{(1)}(0)=0
C              d_{02}^{(2)}(\beta)=0
C     and
C         d d_{00}^{(0)}(\beta)/d\beta=0
C         d d_{01}^{(1)}(\beta)/d\beta=1/\sqrt{2}
C         d d_{02}^{(2)}(\beta)/d\beta=0
C
C  See Eqs. (4.1-4) of \ct{Mis91}:
C
C   (m/\sin\theta) d_{0m}^l(0)=(\delta_{m\pm 1}/2) \sqrt{l(l+1)}
C      d d_{0m}^l(0)/d\beta   =(m\delta_{m\pm 1}/2) \sqrt{l(l+1)}
C
*
*  (4.2.1) of \ct{Ed}:
*   d_{0m}^{(l)}(pi) = (-1)^{l+m} \dt_{0,m}
*
*  (4.2.3) of \ct{Ed}:
*   d_{0m}^{(l)}(0) = (-1)^{m} \dt_{0,m} = \dt_{0,m}
*=======================================
*
*  If X^l_m=(m/\sin\theta) d_{0m}^{(l)}, then, according to (3.29) of {TKS}:
*
*  X^{m+1}_{m+1}=\sin\theta \sqrt{\fr{2m+1}{2m+2}}
*                           \left(\fr{m+1}{m}\right)X^{m}_{m}
*
*  According to (3.30) of {TKS}:
*  X^{m+1}_{m}= -\sqrt{2m+1}\,\cos\theta X^{m}_{m}
*
* According to (3.31) of {TKS}:
*  X^{l}_{m}=\fr{1}{\sqrt{l^2-m^2}}\,\left[(2l-1)\cos\theta
*          X^{l-1}_{m} - \sqrt{(l-1)^2-m^2}}\,\X^{l-2}_{m} \right]
*
* Initial recurrence values are X^1_1=\sqrt{2}/2 and X^l_0=0
***********************************************************************
*                   NONZERO DDV1/DV2 INITIALIZATION
*
*                          M = 0

 100  IF (M.EQ.0) THEN     !all DDV1(N)=X^l_0=0; see (3.33) of {TKS}:

* According to (3.37) of {TKS}, DV2(0)=0.d0

      DV2(1)=-QS

      IF (NMAX.GE.2) DV2(2)=3*X*DV2(1)

      IF (NMAX.LT.3) RETURN
*
      DO N=3,NMAX           !recurrence (3.36) of {TKS},
      DV2(N)=(2*N-1)*X*DV2(N-1)/(N-1)-N*DV2(N-2)/(N-1)
      ENDDO
***********************************************************************
*                           M > 0

       ELSE IF (M.GT.0) THEN
*
* >>> Determine X^m_m according to Eq. (3.29) of {TKS}:

      A=1.d0/DSQRT(2.D0)               !X^1_1=A_1

      DO I=1,M-1
      A=QS*DBLE(I+1)*DSQRT(2*I+1.d0)*A/(I*DSQRT(2*I+2.d0))
      ENDDO

* <<< A is now X^m_m; see (3.29) of {TKS}

      DDV1(M)=A
      DV2(M)=X*A                        !see (3.34) of {TKS}

* >>> Determine X^{m+1}_m:

      IF (M.EQ.NMAX)  GO TO 120

      DER=X*DSQRT(2*M+1.d0)*A          ! DER=X^{m+1}_m; see (3.30) of {TKS}
      DDV1(M+1)=DER
      DV2(M+1)=((M+1)*X*DER-A*DSQRT(2*M+1.d0))/DBLE(M)  !(3.35) of {TKS}

* >>> Determine remaining X^{l}_m's

      IF ((M+2).EQ.NMAX)  GO TO 120

       DO N=M+2,NMAX
       D3=DSQRT(DBLE(N)**2-DBLE(M)**2)
       DDV1(N)=((2*N-1)*X*DDV1(N-1) -
     &                DSQRT(DBLE(N-1)**2-DBLE(M)**2)*DDV1(N-2))/D3
                                                      !see (3.31) of {TKS}
       DV2(N)=(N*X*DDV1(N)-DDV1(N-1)*D3)/DBLE(M)      !see (3.35) of {TKS}
       ENDDO

      END IF

  120 RETURN
      END
C**********************************************************************

      SUBROUTINE CONST (NGAUSS,NMAX,X,W,AN,ANN,S,SS,NP,EPS,RX,HT)
C--------/---------/---------/---------/---------/---------/---------/--
C >>> NGAUSS,NMAX,NP,EPS
C <<< X,W,AN,ANN,S,SS
C=====================
C
C  NGAUSS - the number of GIF division points
C  NMAX - angular momentum cutoff
C  P=PI=DACOS(-1d0)
C  NP - parameter specifying the particle shape
C  EPS - deformation parameter for a given particle shape
C  RX ... 1st char. dimension of a particle
C  HT ... 2nd char. dimension of a particle
C
C  X=\cos\theta  - GIF division points
C  W - GIF weights
C  AN(N)=N*(N+1)
C  ANN(l_1,l_2)=\sqrt{\fr{(2 l_1+1)}{l_1(l_1+1)} }
C                       \sqrt{\fr{(2 l_2+1)}{l_2(l_2+1)} }/2
C  S  ... 1/(|\sin\theta|)
C  SS ... 1/(\sin^2\theta)
C--------/---------/---------/---------/---------/---------/---------/--
      IMPLICIT REAL*8 (A-H,O-Z)
*      INCLUDE 'ampld.par.f'
      Parameter (NPN1=100, NPNG1=2600, NPNG2=2*NPNG1, NPN2=2*NPN1,
     &           NPL=NPN2+1, NPN3=NPN1+1,
     &           NPN4=NPN1, NPN5=2*NPN4, NPN6=NPN4+1)
      INTEGER NEPS,JG,I,J
      REAL*8 EE,EE1,CC,SI,XI1,XI2,XAV
      REAL*8 XTHETA,THETA0,PI,EPS,RX,HT
      REAL*8 X(NPNG2),W(NPNG2),X1(NPNG2),W1(NPNG2),
     *        X2(NPNG2),W2(NPNG2),
     *        S(NPNG2),SS(NPNG2),
     *        AN(NPN1),ANN(NPN1,NPN1),DD(NPN1)
*
      DATA PI/3.141592653589793d0/
*
      if (npng1.le.ngauss) then
      write(6,*)'In CONST increase the value of NPNG1 parameter!'
      write(6,*)'NPNG1=', NPNG1,'should be larger or equal to
     & NGAUSS=', NGAUSS
      stop
      end if
*
      DO 10 N=1,NMAX
           NN=N*(N+1)
           AN(N)=DFLOAT(NN)
           D=DSQRT(DFLOAT(2*N+1)/DFLOAT(NN))
           DD(N)=D
           DO 10 N1=1,N
                DDD=D*DD(N1)*0.5D0
                ANN(N,N1)=DDD
                ANN(N1,N)=DDD
   10 CONTINUE

      NG=2*NGAUSS
*
* GIF division points (ordered from -1 to 1) and weights
*
      NEPS=MAX(EPS,1.d0/EPS)      !pre-number of Gauss integration
                                  !intervals from EPS

      IF (NP.EQ.-1) THEN         ! spheroid

      IF(NEPS.EQ.1) THEN

      CALL GAUSS(NG,0,0,X,W)

      ELSE IF ((NEPS.GE.2).and.(NEPS.LE.2)) THEN

      CALL GAULEG(-1.d0,0.d0,X,W,NGAUSS)

      DO I=1,NGAUSS

         W(I+NGAUSS)=W(NGAUSS-I+1)
         X(I+NGAUSS)=-X(NGAUSS-I+1)

      ENDDO

      ELSE IF (NEPS.GT.2) THEN

c      NEPS=NEPS
      EE=EPS*EPS
      EE1=EE-1D0

      XAV=0.d0

      DO I=1,NEPS

      XI1=DBLE(I)/(NEPS+1)

          CC=XI1*XI1
          SI=1D0-CC
          X2(I)=ABS(XI1*DSQRT(SI)*EE1/DSQRT(SI+EE*CC))   ! ~ |dr(theta)/dtheta|
          XAV=XAV+1.d0/X2(I)

      ENDDO

      XAV=XAV           !sum over 1/|dr(theta)/dtheta|

c_____ estimate integration intervals:

       DO I=1,NEPS

      X2(I)=1.d0/(XAV*X2(I))   !the normalized length of
                               !the i-th integration interval
                               ! \sum x2(i) = 1
      ENDDO

      JG=NGAUSS

      DO I=1,NEPS

      IF(i.eq.1) then

      XI1=0.d0
      XI2=X2(1)

      else

      XI2=XI2+X2(I)

      end if

      NG1=DFLOAT(NGAUSS)*X2(I) !the number of Gauss points on
                               !the i-th integration interval

      IF(I.EQ.NEPS) NG1=NG-JG

      CALL GAULEG(XI1,XI2,X1,W1,NG1)

      XI1=XI2

      DO  J=1,NG1
         W(JG+J)=W1(J)
         X(JG+J)=X1(J)
      ENDDO         !J

      JG=JG+NG1

      ENDDO         !I

*
* Assuming mirror symmetry in the $\theta=\pi/2$ plane
*
      DO  I=1,NGAUSS
         W(I)=W(NG-I+1)
         X(I)=-X(NG-I+1)
      ENDDO

      ENDIF           !NEPS

      ELSE IF (NP.EQ.-2) THEN         ! cylinder

******************   Only involves cylinders  **********************

      NG1=DFLOAT(NGAUSS)/2D0
      NG2=NGAUSS-NG1
      XX=-DCOS(DATAN(EPS))        !-COS OF SEPARATION ANGLE BETWEEN
                                  !HORIZONTAL AND VERTICAL CYLINDER
                                  !FACES
*
* GIF division points and weights
*
      CALL GAUSS(NG1,0,0,X1,W1)         !for (0,NG1)
      CALL GAUSS(NG2,0,0,X2,W2)         !for (NG1+1,NGAUSS=NG1+NG2)
*
C In GAUSS (N,IND1,IND2,Z,W):
C IND1 = 0 - INTERVAL (-1,1),
C IND1 = 1 - (0,1)
C IND2 = 1 RESULTS ARE PRINTED.
*
*
      DO 12 I=1,NG1
         W(I)=0.5D0*(XX+1D0)*W1(I)
         X(I)=0.5D0*(XX+1D0)*X1(I)+0.5D0*(XX-1D0)
   12 CONTINUE
      DO 14 I=1,NG2
         W(I+NG1)=-0.5D0*XX*W2(I)
         X(I+NG1)=-0.5D0*XX*X2(I)+0.5D0*XX
   14 CONTINUE
*
* Assuming mirror symmetry in the $\theta=\pi/2$ plane
*
      DO 16 I=1,NGAUSS
         W(NG-I+1)=W(I)
         X(NG-I+1)=-X(I)
   16 CONTINUE
******************************************************************
      ELSE IF (NP.EQ.-4) THEN         ! cut sphere on top

      XTHETA=DACOS(EPS-1.d0)       !separation angle viewed from the sphere origin

* ratio of the integration lengths along the plane part to the total
* integration length:

      XX=DSIN(XTHETA)/(PI-XTHETA)
      NG2=XX*DBLE(NG)
      NG1=NG-NG2
      THETA0=1.D0/SQRT(8.D0/EPS-3.D0)  !cosine of the separation angle
*
      CALL GAULEG(-1.D0,THETA0,X1,W1,NG1)       !for (0,NG1) along the
                                                !spherical part
      CALL GAULEG(THETA0,1.D0,X2,W2,NG2)        !for (NG2+1,NG=NG1+NG2)
                                                !along the plane part
*
      DO  I=1,NG1
         W(I)=W1(I)
         X(I)=X1(I)
      ENDDO

      DO I=1,NG2

         W(I+NG1)=W2(I)
         X(I+NG1)=X2(I)

      ENDDO
*
******************************************************************
      ELSE IF (NP.EQ.-5) THEN            ! cut sphere on its bottom

      XTHETA=DACOS(EPS-1.d0)
      XX=DSIN(XTHETA)/(PI-XTHETA)
      NG1=XX*DBLE(NG)
      NG2=NG-NG1
      THETA0=-1.D0/SQRT(8.D0/EPS-3.D0)     !cosine of the separation angle
*
      CALL GAULEG(-1.D0,THETA0,X1,W1,NG1)       !for (0,NG1) along the plane part
      CALL GAULEG(THETA0,1.D0,X2,W2,NG2)        !for (NG2+1,NG=NG1+NG2) along
                                                !the spherical part
*
      DO  I=1,NG1
         W(I)=W1(I)
         X(I)=X1(I)
      ENDDO

      DO I=1,NG2

         W(I+NG1)=W2(I)
         X(I+NG1)=X2(I)

      ENDDO
******************************************************************
      ELSE IF (NP.EQ.-6) THEN        ! upwardly oriented cone

      XX=dsqrt(HT**2+RX**2)           !the length of the cone slant
      XTHETA=XX/(XX+RX)               !ratio of the integration lengths along
                                      !the cone slant to the total inteq. length

      NG2=XTHETA*DBLE(NG)
      NG1=NG-NG2

      XX=dsqrt(XX**2+RX**2/2.d0)/2.d0    !=the length of the median of the slant

      THETA0=-HT/(2.d0*XX)               !=cos of the separation angle
                                         !(always negative)
*
      CALL GAULEG(-1.D0,THETA0,X1,W1,NG1)       !for (0,NG1) | along base
      CALL GAULEG(THETA0,1.D0,X2,W2,NG2)        !for (NG1+1,NG=NG1+NG2) | along slant
*
      DO  I=1,NG1
         W(I)=W1(I)
         X(I)=X1(I)
      ENDDO

      DO I=1,NG2

         W(I+NG1)=W2(I)
         X(I+NG1)=X2(I)

      ENDDO

******************************************************************
*
      ELSE
*
      CALL GAUSS(NG,0,0,X,W)
C--------/---------/---------/---------/---------/---------/---------/--
C >>> N,IND1,IND2
C <<< X,W
C=================
C    CALCULATION OF POINTS AND WEIGHTS OF GAUSSIAN QUADRATURE
C    FORMULA. IF IND1 = 0 - ON INTERVAL (-1,1), IF IND1 = 1 - ON
C    INTERVAL  (0,1). IF  IND2 = 1 RESULTS ARE PRINTED.
C
C    N - NUMBER OF GIF DIVISION POINTS (mostly N=NGAUSS in main program)
C    X - DIVISION POINTS (FROM -1 to 1)
C    W - WEIGHTS
C--------/---------/---------/---------/---------/---------/---------/--
c
c      CALL GAULEG(-1.D0,1.D0,X,W,NG)
*
      END IF

      if (np.gt.-4) then           !mirror symmetry present

       DO 20 I=1,NGAUSS
           Y=X(I)
           Y=1D0/(1D0-Y*Y)
           SS(I)=Y                 !1/sin**2(theta)
           SS(NG-I+1)=Y
           Y=DSQRT(Y)              !1/|sin(theta)|
           S(I)=Y
           S(NG-I+1)=Y
   20 CONTINUE

      else                         !mirror symmetry absent

       DO 30 I=1,NG
           Y=X(I)
           Y=1D0/(1D0-Y*Y)         !1/sin**2(theta)
           SS(I)=Y
           Y=DSQRT(Y)              !1/|sin(theta)|
           S(I)=Y
   30 CONTINUE

      END IF

      RETURN
      END

C**********************************************************************

      SUBROUTINE VARY (LAM,MRR,MRI,A,EPS,RSNM,HT,NP,NGAUSS,X,
     *                 P,PPI,PIR,PII,R,DR,DDR,DRR,DRI,NMAX)
C--------/---------/---------/---------/---------/---------/---------/--
C >>> LAM,MRR,MRI,A,EPS,NP,NGAUSS,X,P,NMAX
C <<< PPI,PIR,PII,R,DR,DDR,DRR,DRI
C=========================
C  LAM - wavelength of incident light in the ambient
C  MRR - the real part of the refractive index
C  MRI - the imaginary  part of the refractive index
C  A=RAT*AXI, where RAT and AXI are the main program input parameters
C      RAT = 1 - particle size is specified in terms of the
C                equal-volume-sphere radius
C      RAT.ne.1 - particle size is specified in terms of the
C                equal-surface-area-sphere radius
C  AXI - equivalent-(volume/surface-area)-sphere radius
C  NP - particle shape class
C  EPS - shape deformation parameter within a given particle shape class
C  NGAUSS - the number of Gauss integration division points
C           in the integral over theta
C  NMAX - angular momentum cutoff
C  P=DACOS(-1D0)
C  PI=P*2D0/LAM - wave vector
C  PPI=PI*PI
C  PIR=PPI*MRR
C  PII=PPI*MRI
C  R=r^2(\theta)            for axially symmetric particles
C  DR=dr(\theta)/(d\theta)  for axially symmetric particles
C  DDR=\lambda/[2*\pi*r(\theta)]
C  DRR=(MRR/(MRR**2+MRI**2))*(\lambda/[2*\pi*r(\theta)])
C  DRI=-(MRI/(MRR**2+MRI**2))*(\lambda/[2*\pi*r(\theta)])
C--------/---------/---------/---------/---------/---------/---------/--
*      INCLUDE 'ampld.par.f'
      Parameter (NPN1=100, NPNG1=2600, NPNG2=2*NPNG1, NPN2=2*NPN1,
     &           NPL=NPN2+1, NPN3=NPN1+1,
     &           NPN4=NPN1, NPN5=2*NPN4, NPN6=NPN4+1)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8  X(NPNG2),R(NPNG2),DR(NPNG2),MRR,MRI,LAM,
     *        Z(NPNG2),ZR(NPNG2),ZI(NPNG2),
     *        DDR(NPNG2),DRR(NPNG2),DRI(NPNG2)
cc     *        J(NPNG2,NPN1),Y(NPNG2,NPN1),JR(NPNG2,NPN1),
cc     *        JI(NPNG2,NPN1),DJ(NPNG2,NPN1),DY(NPNG2,NPN1),
cc     *        DJR(NPNG2,NPN1),DJI(NPNG2,NPN1)
cc      COMMON /CBESS/ J,Y,JR,JI,DJ,DY,DJR,DJI

      NG=NGAUSS*2

* decision tree to specify particle shape:

      IF (NP.GT.0)  CALL RSP2(X,NG,A,EPS,NP,R,DR)       ! Chebyshev particle
      IF (NP.EQ.-1) CALL RSP1(X,NG,NGAUSS,A,EPS,R,DR)   ! oblate/prolate spheroids
      IF (NP.EQ.-2) CALL RSP3(X,NG,NGAUSS,A,EPS,R,DR)   ! oblate/prolate cylinder
      IF (NP.EQ.-3) CALL RSP4(X,NG,A,R,DR)              ! a distorted Chebyshev droplet
      IF (NP.EQ.-4) CALL RSP5(X,NG,RSNM,EPS,R,DR)       ! sphere cut by a plane on its top
      IF (NP.EQ.-5) CALL RSPI5(X,NG,RSNM,EPS,R,DR)      ! sphere cut by a plane on its bottom
      IF (NP.EQ.-6) CALL RSP6(X,NG,RSNM,HT,R,DR)        ! upwardly oriented cone
cc      IF (NP.EQ.-7) CALL RSP7(X,NG,RSNM,HT,R,DR)          ! cone cut on its top
cc      IF (NP.EQ.-8) CALL RSP8(X,NG,RSNM,HT,R,DR)          ! cone on a cylinder

*
      PI=P*2D0/LAM                 !wave vector
      PPI=PI*PI
      PIR=PPI*MRR
      PII=PPI*MRI
      V=1D0/(MRR*MRR+MRI*MRI)
      PRR=MRR*V
      PRI=-MRI*V
      TA=0D0
      DO 10 I=1,NG
           VV=DSQRT(R(I))
           V=VV*PI
           TA=MAX(TA,V)            !Max. size parameter
           VV=1D0/V
           DDR(I)=VV
           DRR(I)=PRR*VV
           DRI(I)=PRI*VV
           V1=V*MRR
           V2=V*MRI
           Z(I)=V              !=(2\pi/\lambda)*r
           ZR(I)=V1            !=(2\pi/\lambda)*r*MRR
           ZI(I)=V2            !=(2\pi/\lambda)*r*MRI
   10 CONTINUE
      IF (NMAX.GT.NPN1) PRINT 9000,NMAX,NPN1
      IF (NMAX.GT.NPN1) STOP
 9000 FORMAT(' NMAX = ',I2,', i.e., greater than ',I3)
*
* TA is the ``max. size parameter", MAX(2*PI*SQRT(RI)/LAMBDA)

      TB=TA*DSQRT(MRR*MRR+MRI*MRI)     !=TA*EPSIN
      TB=DMAX1(TB,DFLOAT(NMAX))
*
      NNMAX1=1.2D0*DSQRT(DMAX1(TA,DFLOAT(NMAX)))+3D0
      NNMAX2=(TB+4D0*(TB**0.33333D0)+1.2D0*DSQRT(TB))  !Wiscombe bound
      NNMAX2=NNMAX2-NMAX+5
*
* generate arrays of Bessel functions at NGAUSS GIF division
* points and store them in the common block /CBESS/
*
      CALL BESS(Z,ZR,ZI,NG,NMAX,NNMAX1,NNMAX2)
*
      RETURN
      END

C**********************************************************************

      SUBROUTINE RSP1(X,NG,NGAUSS,REV,EPS,R,DR)
C--------/---------/---------/---------/---------/---------/---------/--
C >>> X,NG,NGAUSS,REV,EPS
C <<< R,DR
C=========================
C   Activated for NP=-1
C
C   Calculation of the functions
C              R(I)=r(y)**2 and DR(I)=((d/dy)r(y))/r(y)
C   for an oblate/prolate spheroids droplet specified by the parameters
C   REV and EPS at NGAUSS Gauss integration formula (GIF) division points
C   in the integral over theta. Here Y=ACOS(X)=THETA.
C
C      r(\theta,\phi)=a\left[\sin^2\theta + (a^2/b^2)\cos^2\theta]^{-1/2}
C     dr(\theta,\phi)/d\theta =r(\theta,\phi)*\cos\theta*\sin\theta
C                        *[(a^2/b^2)-1]/[\sin^2\theta + (a^2/b^2)\cos^2\theta]
C
C   X - GIF division points \cos\theta_j -  Y = arccos X
C   REV ... equal-volume-sphere radius
C   EPS=a/b = the ratio of the horizontal to rotational axes.  EPS is
C             larger than 1 for oblate spheroids and smaller than 1 for
C             prolate spheroids.
C   NGAUSS ... the number of GIF division points
C   NG=2*NGAUSS
C
C   1.LE.I.LE.NGAUSS
C--------/---------/---------/---------/---------/---------/---------/--
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 X(NG),R(NG),DR(NG)

      A=REV*EPS**(1.D0/3.D0)
      AA=A*A
      EE=EPS*EPS
      EE1=EE-1D0

      DO 50 I=1,NGAUSS
          C=X(I)
          CC=C*C
          SS=1D0-CC
          S=DSQRT(SS)                !=\sin\theta
          RR=1D0/(SS+EE*CC)
          R(I)=AA*RR                 !=r(\theta)**2
          R(NG-I+1)=R(I)
          DR(I)=RR*C*S*EE1           !=[dr(theta)/d theta]/r(theta)
          DR(NG-I+1)=-DR(I)
   50 CONTINUE

      RETURN
      END

C**********************************************************************

      SUBROUTINE RSP2 (X,NG,REV,EPS,N,R,DR)
C--------/---------/---------/---------/---------/---------/---------/--
C >>> X,NG,REV,EPS,N
C <<< R,DR
C=========================
C   Activated for NP.gt.0
C
C   Calculation of the functions R(I)=r(y)**2 and
C   DR(I)=((d/dy)r(y))/r(y) for a Chebyshev particle
C   specified by the parameters REV, EPS, and N (Y=ACOS(X)=THETA).
C
C       r(\theta,\phi)=r_0[1+\eps T_n(\cos\theta)]    (*)
C
C   EPS ... deformation parameter of a Chebyshev particle; |EPS|<1
C   N   ... the degree of the Chebyshev polynomial
C   All Chebyshev particles with N.GE.2 become partially concave
C   as the absolute value of the deformation parameter EPS increases
C   and exhibit surface roughness in the form of waves running
C   completely around the particle.
C
C   X - GIF division points \cos\theta_j -  Y = arccos X
C   REV ... equal-volume-sphere radius r_ev
C   NGAUSS ... the number of GIF division points
C   NG=2*NGAUSS
C
C   1.LE.I.LE.NGAUSS
C
C--------/---------/---------/---------/---------/---------/---------/--
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 X(NG),R(NG),DR(NG)
      DNP=DFLOAT(N)
      DN=DNP*DNP
      DN4=DN*4D0                   !=4*N**2
      EP=EPS*EPS                   !=EPS**2
      A=1D0+1.5D0*EP*(DN4-2D0)/(DN4-1D0)
      I=(DNP+0.1D0)*0.5D0
      I=2*I
      IF (I.EQ.N) A=A-3D0*EPS*(1D0+0.25D0*EP)/
     *              (DN-1D0)-0.25D0*EP*EPS/(9D0*DN-1D0)
      R0=REV*A**(-1D0/3D0)
      DO 50 I=1,NG
         XI=DACOS(X(I))*DNP
         RI=R0*(1D0+EPS*DCOS(XI))    !the Chebyshev shape function (*)
         R(I)=RI*RI
         DR(I)=-R0*EPS*DNP*DSIN(XI)/RI
c        WRITE(NOUT,*) I,R(I),DR(I)
   50 CONTINUE

      RETURN
      END

C**********************************************************************

      SUBROUTINE RSP3 (X,NG,NGAUSS,REV,EPS,R,DR)
C--------/---------/---------/---------/---------/---------/---------/--
C >>> X,NG,NGAUSS,REV,EPS
C <<< R,DR
C=========================
C   Activated for NP=-2
C
C   Calculation of the functions R(I)=r(y)**2 and
C   DR(I)=((d/dy)r(y))/r(y) for an oblate/prolate cylinder
C   specified by the parameters REV and EPS  at NGAUSS  Gauss
C   integration points in the integral over theta.
C
C   X - GIF division points \cos\theta_j -  Y = arccos X
C   REV ... equal-volume-sphere radius r_ev
C   EPS ... the ratio of the cylinder diameter to its length
C   H   ... half-length of the cylinder
C   A=H*EPS  ... cylinder radius   ====>
C
C   4*PI*REV**3/3=2*H*PI*A**2=2*PI*H**3*EPS**2 <====>
C                H=REV*( (2D0/(3D0*EPS*EPS))**(1D0/3D0) )
C
C
C   NGAUSS ... the number of GIF division points
C   NG=2*NGAUSS
C
C   1.LE.I.LE.NGAUSS
C
C--------/---------/---------/---------/---------/---------/---------/--
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 X(NG),R(NG),DR(NG)

* Determine half-length of the cylinder
      H=REV*( (2D0/(3D0*EPS*EPS))**(1D0/3D0) )

* Determine cylinder radius:
      A=H*EPS

      DO 50 I=1,NGAUSS
         CO=-X(I)
         SI=DSQRT(1D0-CO*CO)

         IF (SI/CO.GT.A/H) GO TO 20

* Along the plane cuts:

         RAD=H/CO
         RTHET=H*SI/(CO*CO)
         GO TO 30

* Along the circular surface:
   20    CONTINUE
         RAD=A/SI
         RTHET=-A*CO/(SI*SI)
cc         RAD=1.D-10
cc         RTHET=0.D0

   30    R(I)=RAD*RAD
         R(NG-I+1)=R(I)          !using mirror symmetry

         DR(I)=-RTHET/RAD
         DR(NG-I+1)=-DR(I)       !using mirror symmetry

   50 CONTINUE

      RETURN
      END

C**********************************************************************

      SUBROUTINE RSP4 (X,NG,REV,R,DR)
C--------/---------/---------/---------/---------/---------/---------/--
C >>> X,NG
C <<< R,DR
C=========================
C   Activated for NP=-3
C
C   Calculation of the functions R(I)=r(y)**2 and
C   DR(I)=((d/dy)r(y))/r(y) for a distorted
C   droplet (generalized Chebyshev particle) specified by the
C   parameters REV and c_n (Chebyshev expansion coefficients).
C   (Y=ACOS(X)=THETA).
C   The coefficients of the Chebyshev polynomial expansion are
C   specified in the subroutine DROP.
C
C   X - GIF division points \cos\theta_j -  Y = arccos X
C   REV ... equal-volume-sphere radius  r_ev
C   NGAUSS ... the number of GIF division points
C   NG=2*NGAUSS
C
C   1.LE.I.LE.NGAUSS
C
C--------/---------/---------/---------/---------/---------/---------/--
      PARAMETER (NC=10)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 X(NG),R(NG),DR(NG),C(0:NC)
      COMMON /CDROP/ C,R0V
      R0=REV*R0V
      DO I=1,NG
         XI=DACOS(X(I))
         RI=1D0+C(0)
         DRI=0D0
         DO N=1,NC
            XIN=XI*N
            RI=RI+C(N)*DCOS(XIN)
            DRI=DRI-C(N)*N*DSIN(XIN)
         ENDDO
         RI=RI*R0
         DRI=DRI*R0
         R(I)=RI*RI
         DR(I)=DRI/RI
c        WRITE(NOUT,*) I,R(I),DR(I)
      ENDDO

      RETURN
      END


C**********************************************************************

      SUBROUTINE RSP5 (X,NG,REV,EPS,R,DR)
C--------/---------/---------/---------/---------/---------/---------/--
C >>> X,NG,REV,EPS
C <<< R,DR
C=========================
C   Activated for NP=-4
C
C   Calculation of the functions
C              R(I)=r(y)**2 and DR(I)=((d/dy)r(y))/r(y)
C   for a sphere cut by the plane specified by the parameters
C   REV and  EPS at NGAUSS Gauss integration formula (GIF) division points
C   in the integral over theta. Here Y=ACOS(X)=THETA and EPS=2*R_0/H,
C   where R_0 is the radius of the original uncut sphere, whereas H
C   is the height (along the axial symmetry axis) of the resulting
C   cut sphere.
C
C   The origin of coordinates is located along the axis of symmetry
C                  midway the plane and sphere top.
C
C             ===>    Note that always EPS.GT.1
C ===
C   X - GIF division points \cos\theta_j -  Y = arccos X
C   REV ... the radius of the original uncut sphere
C   EPS ...  H is the height (along the axial symmetry axis)
C            of the resulting cut sphere. Note that always EPS.LT.2*REV
C   THETA0 ... a COSINE of the separation angle between two different
C              functional dependences of r(\theta), that along the sphere
C              surface and that along the plane surface
C   NG=2*NGAUSS ... the number of GIF division points
C
C   1.LE.I.LE.NGAUSS
C--------/---------/---------/---------/---------/---------/---------/--
      IMPLICIT REAL*8 (A-H,O-Z)

      INTEGER I,NG
      REAL*8 A,CC,CO,SS,SI,RAD,REV,EPS,THETA0,RTHET
      REAL*8 X(NG),R(NG),DR(NG)

      IF (EPS.GE.2.d0*REV) THEN
      WRITE(6,*)'Invalid parameters for a cut sphere!'
      WRITE(6,*)'Execution stopped!'
      STOP
      END IF

      THETA0=1.D0/SQRT(8.D0*REV/EPS-3.D0)  !cosine of the separation angle
                                           !    (always positive)
      DO 50 I=1,NG

          CO=X(I)
          CC=CO*CO
          SS=1.D0-CC
          SI=DSQRT(SS)                    !=\sin\theta

          IF (CO.LT.THETA0) THEN        ! r(\theta) along the sphere surface

          A=REV-EPS/2.D0

          RAD=A*CO+SQRT(REV**2-(A*SI)**2)
          RTHET=-A*SI - CO*SI*A**2/SQRT(REV**2-(A*SI)**2)
cc          RAD=1.D-10
cc          RTHET=0.D0


          ELSE IF (CO.GE.THETA0) THEN   ! r(\theta) along the plane surface
                                        !    (CO positive)
          RAD=EPS/(2.D0*CO)
          RTHET=EPS*SI/(2.D0*CO**2)

          END IF

          DR(I)=RTHET/RAD
          R(I)=RAD*RAD

   50 CONTINUE

      RETURN
      END

C**********************************************************************

      SUBROUTINE RSPI5 (X,NG,REV,EPS,R,DR)
C--------/---------/---------/---------/---------/---------/---------/--
C >>> X,NG,REV,EPS
C <<< R,DR
C=========================
C   Activated for NP=-5
C
C   SIMILAR TO RSP5, EXCEPT FOR THAT THE PLANE CUT IS ON THE SPHERE "BOTTOM"
C   ===> cosine of the separation angle has the same magnitude as in RSP5,
C        but is always negative in the present case !!!
C
C   Calculation of the functions
C              R(I)=r(y)**2 and DR(I)=((d/dy)r(y))/r(y)
C   for a sphere cut by the plane specified by the parameters
C   REV and  EPS at NGAUSS Gauss integration formula (GIF) division points
C   in the integral over theta. Here Y=ACOS(X)=THETA and EPS=2*R_0/H,
C   where R_0 is the radius of the original uncut sphere, whereas H
C   is the height (along the axial symmetry axis) of the resulting
C   cut sphere.
C
C   The origin of coordinates is located along the axis of symmetry
C                  midway the plane and sphere top.
C
C                  ===>  Note that always EPS.GT.1
C   ===
C   X - GIF division points \cos\theta_j -  Y = arccos X
C   REV ... the radius of the original uncut sphere
C   EPS ...  H is the height (along the axial symmetry axis)
C            of the resulting cut sphere. Note that always EPS.LT.2*REV
C   THETA0 ... a COSINE of the separation angle between two different
C              functional dependences of r(\theta), that along the sphere
C              surface and that along the plane surface
C   NG=2*NGAUSS ... the number of GIF division points
C
C   1.LE.I.LE.NGAUSS
C--------/---------/---------/---------/---------/---------/---------/--
      IMPLICIT REAL*8 (A-H,O-Z)

      INTEGER I,NG
      REAL*8 A,CC,CO,SS,SI,RAD,REV,EPS,THETA0,RTHET
      REAL*8 X(NG),R(NG),DR(NG)

      IF (EPS.GE.2.d0*REV) THEN
      WRITE(6,*)'Invalid parameters for a cut sphere!'
      WRITE(6,*)'Execution stopped!'
      STOP
      END IF

      THETA0=-1.D0/SQRT(8.D0*REV/EPS-3.D0)  !cosine of the separation angle
                                            !is always negative in the present
                                            !case
      DO 50 I=1,NG

          CO=X(I)
          CC=CO*CO
          SS=1.D0-CC
          SI=DSQRT(SS)                  !=\sin\theta

          IF (CO.GT.THETA0) THEN        ! r(\theta) along the sphere surface

          A=REV-EPS/2.D0

          RAD=-A*CO+SQRT(REV**2-(A*SI)**2)
          RTHET=A*SI - CO*SI*A**2/SQRT(REV**2-(A*SI)**2)
cc          RAD=1.D-10
cc          RTHET=0.D0

          ELSE IF (CO.LE.THETA0) THEN   ! r(\theta) along the plane surface
                                        !        (CO negative)
          RAD=-EPS/(2.D0*CO)
          RTHET=-EPS*SI/(2.D0*CO**2)

          END IF

          DR(I)=RTHET/RAD
          R(I)=RAD*RAD

   50 CONTINUE

      RETURN
      END

C*********************************************************************

      SUBROUTINE RSP6(X,NG,REV,HT,R,DR)
C--------/---------/---------/---------/---------/---------/---------/--
C >>> X,NG,REV,HT
C <<< R,DR
C=========================
C   Activated for NP=-6
C
C   Calculation of the functions
C              R(I)=r(y)**2 and DR(I)=((d/dy)r(y))/r(y)
C   for an upwardly pointing singular cone specified by the parameters
C   REV and EPS at NGAUSS Gauss integration formula (GIF) division points
C   in the integral over theta.
C   REV ... the the half width of the cone base
C   HT ...  cone height (along the axial symmetry axis)
C
C   The origin of coordinates is placed on the axis of axial symmetry,
C         at the distance (H/3) from the base to the cone
C                                    <===>
C     The base and slant of the cone form a triangle, and the origin
C          of coordinates is placed at the triangle centroid.
C
C   ===>  the length of the cone slant A = DSQRT(H**2+REV**2)
C   ===>  the length of the median of the cone slant MA = dsqrt(ma**2+rev**2/2.d0)/2.d0
C
C   The cone base end points and the traingle centroid form a second triangle.
C   If we denote by 2*Theta1 the second triangle angle at the centroid,
C   the separation angle theta_s as viewed from the origin (between the slant and
C   base part of the cone surfaces) is given by
C
C               - cos (theta_s) =  cos (theta1) = (h/3)/(2*MA/3)  = h/(2*MA)
C
C   Then, for  theta>theta_s (i.e., cos (theta) < cos (theta_s) ):
C
C                      r(theta)=-h/(3*cos(theta))
C
C   For  theta<theta_s (i.e., cos (theta) > cos (theta_s) ), one first considers
C   a triangle formed by the cone apex, the origin of coordinates and r(theta).
C   Let Theta_v denote the angle of this triangle at the cone apex.
C   According to the law of cosines:
C
C         c^2 = r^2+4h^2/9 - [4rh\cos(\theta)]/3
C         r^2 = c^2+4h^2/9 - [4rh\cos(\theta_v)]/3
C
C   where c is the length of the traingle opposite to the origin of coordinates.
C   Upon combining the last two equations one arrives
C
C                 r\cos(\theta)+c\cos(\theta_v) = 2h/3    (*)
C
C   According to the law of sines:
C
C                      c/sin(\theta) = r(\theta)/sin(\theta_v)
C
C   When substituting back to (*), one finds
C
C               r = (2h/3)/[\cos(\theta)+\cot(\theta_v) \sin(\theta)]
C
C   \cot(\theta_v) can be easily determined from the very first triangle,
C
C              \cot(\theta_v) = h/(b/2) = 2h/b
C
C   Conus volume = (PI*REV**2*H)/3.d0
C   Here Y=ACOS(X)=THETA and EPS=2*R_0/H,
C   where R_0 is the radius of the original uncut sphere, whereas H
C   is the height (along the axial symmetry axis) of the resulting
C   cut sphere.
C   ===>
C
C   X - GIF division points \cos\theta_j -  Y = arccos X
C
C   NG=2*NGAUSS ... the number of GIF division points
C
C   1.LE.I.LE.NGAUSS
C--------/---------/---------/---------/---------/---------/---------/--
      IMPLICIT NONE

      INTEGER I,NG
      REAL*8 HT,MA,CO,SI,CC,SS,RAD,REV,THETA0,RTHET
      REAL*8 X(NG),R(NG),DR(NG)

      MA=DSQRT(HT**2+REV**2)              !=the length of the cone slant
      MA=dsqrt(ma**2+8.d0*rev**2)/2.d0    !=the length of the median of the slant

      THETA0=-HT/(2.d0*MA)                !=cos of the separation angle
                                          !  (always negative)

      DO 50 I=1,NG

          CO=X(I)
          CC=CO*CO
          SS=1.D0-CC
          SI=DSQRT(SS)                    !=\sin\theta

          IF (CO.GT.THETA0) THEN          ! theta < theta0,
                                          ! i.e. r(\theta) along the cone slant surface

          MA=CO+HT*SI/REV
          RAD=2.d0*HT/(3.d0*MA)
          RTHET=2.d0*HT*(SI-HT*CO/REV)/(3.d0*MA**2)

          ELSE IF (CO.LE.THETA0) THEN     ! theta > theta0,
                                          ! i.e. r(\theta) along the cone base surface

          RAD=-HT/(3.d0*CO)               !always positive (CO<0 on the base)
          RTHET=RAD*SI/CO                 !always negative

          END IF

          R(I) = RAD*RAD
          DR(I)= RTHET/RAD

   50 CONTINUE

      RETURN
      END


C*********************************************************************

      SUBROUTINE BESS(X,XR,XI,NG,NMAX,NNMAX1,NNMAX2)
C--------/---------/---------/---------/---------/---------/---------/--
C >>> X,XR,XI,NG,NMAX,NNMAX1,NNMAX2
C <<< Output J,Y,JR,JI,DJ,DY,DJR,DJI  to common block CBESS
C==========================================================
C  Generates Bessel functions for each Gauss integration point
C
C  X =(2\pi/\lambda)*r
C  XR=(2\pi/\lambda)*r*MRR, MRR ... real part of the rel. refractive index
C  XI=(2\pi/\lambda)*r*MRI, MRI ... imag. part of the rel. refractive index
C  NG=2*NGAUSS or 60 ... the number of Gauss integration points
C  J,Y,JR,JI ... arrays of Bessel functions
C  DJ,DY,DJR,DJI  ... arrays of Bessel functions derivatives of the form
C                           [xf(x)]'/x                   (A)
C                 where prime denotes derivative with respect to x.
C                 (Note that Bessel function derivatives enter Eqs. (39)
C                  \cite{TKS} only in the (A) combination!!!!)
C  NMAX   ... angular momentum cutoff
C  NNMAX1 ... angular momentum cutoff - DETERMINES NUMERICAL ACCURACY
C  NNMAX2 ... angular momentum cutoff - DETERMINES NUMERICAL ACCURACY
C--------/---------/---------/---------/---------/---------/---------/--
*      INCLUDE 'ampld.par.f'
      Parameter (NPN1=100, NPNG1=2600, NPNG2=2*NPNG1, NPN2=2*NPN1,
     &           NPL=NPN2+1, NPN3=NPN1+1,
     &           NPN4=NPN1, NPN5=2*NPN4, NPN6=NPN4+1)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 X(NG),XR(NG),XI(NG),
     *        J(NPNG2,NPN1),Y(NPNG2,NPN1),JR(NPNG2,NPN1),
     *        JI(NPNG2,NPN1),DJ(NPNG2,NPN1),DY(NPNG2,NPN1),
     *        DJR(NPNG2,NPN1),DJI(NPNG2,NPN1),
     *        AJ(NPN1),AY(NPN1),AJR(NPN1),AJI(NPN1),
     *        ADJ(NPN1),ADY(NPN1),ADJR(NPN1),
     *        ADJI(NPN1)
      COMMON /CBESS/ J,Y,JR,JI,DJ,DY,DJR,DJI    !arrays of generated Bessel functions
*
      DO 10 I=1,NG
           XX=X(I)
*
           CALL RJB(XX,AJ,ADJ,NMAX,NNMAX1)
           CALL RYB(XX,AY,ADY,NMAX)
*
           YR=XR(I)
           YI=XI(I)
*
           CALL CJB(YR,YI,AJR,AJI,ADJR,ADJI,NMAX,2)
*
           DO 10 N=1,NMAX
                J(I,N)=AJ(N)
                Y(I,N)=AY(N)
                JR(I,N)=AJR(N)
                JI(I,N)=AJI(N)
                DJ(I,N)=ADJ(N)
                DY(I,N)=ADY(N)
                DJR(I,N)=ADJR(N)
                DJI(I,N)=ADJI(N)
   10 CONTINUE

      RETURN
      END

C**********************************************************************

      SUBROUTINE RJB(X,Y,U,NMAX,NNMAX)
C--------/---------/---------/---------/---------/---------/---------/--
C >>> EPS
C <<< RAT
C=================
C  X =(2\pi/\lambda)*r
C  Y ...
C  NMAX - angular momentum cutoff
C  NNMAX - angular momentum cutoff - DETERMINES NUMERICAL ACCURACY
C--------/---------/---------/---------/---------/---------/---------/--
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 Y(NMAX),U(NMAX),Z(800)
*
      L=NMAX+NNMAX
      XX=1D0/X
      Z(L)=1D0/(DFLOAT(2*L+1)*XX)
      L1=L-1
      DO 5 I=1,L1
         I1=L-I
         Z(I1)=1D0/(DFLOAT(2*I1+1)*XX-Z(I1+1))
    5 CONTINUE
      Z0=1D0/(XX-Z(1))
      Y0=Z0*DCOS(X)*XX
      Y1=Y0*Z(1)
      U(1)=Y0-Y1*XX
      Y(1)=Y1
      DO 10 I=2,NMAX
         YI1=Y(I-1)
         YI=YI1*Z(I)
         U(I)=YI1-DFLOAT(I)*YI*XX
         Y(I)=YI
   10 CONTINUE

      RETURN
      END

C**********************************************************************

      SUBROUTINE RYB(X,Y,V,NMAX)
C--------/---------/---------/---------/---------/---------/---------/--
C >>> EPS
C <<< RAT
C=================
C  X =(2\pi/\lambda)*r
C  NMAX - angular momentum cutoff
C--------/---------/---------/---------/---------/---------/---------/--
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 Y(NMAX),V(NMAX)
*
      C=DCOS(X)
      S=DSIN(X)
      X1=1D0/X
      X2=X1*X1
      X3=X2*X1
      Y1=-C*X2-S*X1
      Y(1)=Y1
      Y(2)=(-3D0*X3+X1)*C-3D0*X2*S
      NMAX1=NMAX-1
      DO 5 I=2,NMAX1
    5     Y(I+1)=DFLOAT(2*I+1)*X1*Y(I)-Y(I-1)
      V(1)=-X1*(C+Y1)
      DO 10 I=2,NMAX
  10       V(I)=Y(I-1)-DFLOAT(I)*X1*Y(I)
      RETURN
      END

C**********************************************************************

      SUBROUTINE CJB(XR,XI,YR,YI,UR,UI,NMAX,NNMAX)
C--------/---------/---------/---------/---------/---------/---------/--
C
C   CALCULATION OF SPHERICAL BESSEL FUNCTIONS OF THE FIRST KIND
C   J=JR+I*JI OF COMPLEX ARGUMENT X=XR+I*XI OF ORDERS FROM 1 TO NMAX
C   BY USING BACKWARD RECURSION. PARAMETER NNMAX DETERMINES NUMERICAL
C   ACCURACY. U=UR+I*UI - FUNCTION (1/X)(D/DX)(X*J(X))=J(X)/X + J'(X)
C
C  XR=(2\pi/\lambda)*r*MRR, MRR ... real part of the rel. refractive index
C  XI=(2\pi/\lambda)*r*MRI, MRI ... imag. part of the rel. refractive index
C
C   NMAX  - angular momentum cutoff
C   NNMAX - angular momentum cutoff - DETERMINES NUMERICAL ACCURACY
C--------/---------/---------/---------/---------/---------/---------/--
*      INCLUDE 'ampld.par.f'
      Parameter (NPN1=100, NPNG1=2600, NPNG2=2*NPNG1, NPN2=2*NPN1,
     &           NPL=NPN2+1, NPN3=NPN1+1,
     &           NPN4=NPN1, NPN5=2*NPN4, NPN6=NPN4+1)
      IMPLICIT REAL*8 (A-H,O-Z)

      REAL*8 YR(NMAX),YI(NMAX),UR(NMAX),UI(NMAX)
      REAL*8 CYR(NPN1),CYI(NPN1),CZR(1200),CZI(1200)
c     *       CUR(NPN1),CUI(NPN1)
*
      L=NMAX+NNMAX
      XRXI=1D0/(XR*XR+XI*XI)
      CXXR=XR*XRXI             !Re [1/(XR+i*XI)]
      CXXI=-XI*XRXI            !Im [1/(XR+i*XI)]
      QF=1D0/DFLOAT(2*L+1)
      CZR(L)=XR*QF
      CZI(L)=XI*QF
      L1=L-1
      DO I=1,L1
         I1=L-I
         QF=DFLOAT(2*I1+1)
         AR=QF*CXXR-CZR(I1+1)
         AI=QF*CXXI-CZI(I1+1)
         ARI=1D0/(AR*AR+AI*AI)
         CZR(I1)=AR*ARI
         CZI(I1)=-AI*ARI
      ENDDO

      AR=CXXR-CZR(1)
      AI=CXXI-CZI(1)
      ARI=1D0/(AR*AR+AI*AI)
      CZ0R=AR*ARI
      CZ0I=-AI*ARI
      CR=DCOS(XR)*DCOSH(XI)
      CI=-DSIN(XR)*DSINH(XI)
      AR=CZ0R*CR-CZ0I*CI
      AI=CZ0I*CR+CZ0R*CI
      CY0R=AR*CXXR-AI*CXXI
      CY0I=AI*CXXR+AR*CXXI
      CY1R=CY0R*CZR(1)-CY0I*CZI(1)
      CY1I=CY0I*CZR(1)+CY0R*CZI(1)
      AR=CY1R*CXXR-CY1I*CXXI
      AI=CY1I*CXXR+CY1R*CXXI
      CU1R=CY0R-AR
      CU1I=CY0I-AI
      CYR(1)=CY1R
      CYI(1)=CY1I
c      CUR(1)=CU1R
c      CUI(1)=CU1I
      YR(1)=CY1R
      YI(1)=CY1I
      UR(1)=CU1R
      UI(1)=CU1I

      DO I=2,NMAX
         QI=DFLOAT(I)
         CYI1R=CYR(I-1)
         CYI1I=CYI(I-1)
         CYIR=CYI1R*CZR(I)-CYI1I*CZI(I)
         CYII=CYI1I*CZR(I)+CYI1R*CZI(I)
         AR=CYIR*CXXR-CYII*CXXI            !Re [J/(XR+i*XI)]
         AI=CYII*CXXR+CYIR*CXXI            !Im [J/(XR+i*XI)]
         CUIR=CYI1R-QI*AR
         CUII=CYI1I-QI*AI
         CYR(I)=CYIR
         CYI(I)=CYII
c         CUR(I)=CUIR
c         CUI(I)=CUII
         YR(I)=CYIR
         YI(I)=CYII
         UR(I)=CUIR
         UI(I)=CUII
      ENDDO
*
      RETURN
      END

C**********************************************************************

      SUBROUTINE TMATR0(NGAUSS,X,W,AN,ANN,PPI,PIR,PII,R,DR,DDR,
     *                  DRR,DRI,NMAX,NCHECK,NAXSM)
C--------/---------/---------/---------/---------/---------/---------/--
C >>> NGAUSS,X,W,AN,ANN,PPI,PIR,PII,R,DR,DDR,DRR,DRI,NMAX,NCHECK
C <<< common blocks /TMAT99/, /CT/ (for main),  and /CTT/ (for TT)
C=====================
C
C  Determines the T-matrix of an axially symmetric scatterer
C                           for M=0
C
C  NGAUSS - the number of GIF division points
C  X=\cos\theta  - GIF division points
C  W - GIF weights
C  AN(N)=N*(N+1)
C  ANN(l_1,l_2)=\sqrt{\fr{(2 l_1+1)}{l_1(l_1+1)} }
C                       \sqrt{\fr{(2 l_2+1)}{l_2(l_2+1)} }/2
C  NMAX - angular momentum cutoff
C  NCHECK  -  .EQ.0  THEN  NGSS=2*NGAUSS, FACTOR=1D0
C             .EQ.1  THEN  NGSS = NGAUSS, FACTOR=2D0
C  NAXSM   -  .EQ.0 : Gauss abscissas not +/- theta symmetry
C             .EQ.1 : Gauss abscissas  +/- theta symmetric
C  P=DACOS(-1D0)
C  PI=P*2D0/LAM - wave vector
C  PPI=PI*PI
C  PIR=PPI*MRR
C  PII=PPI*MRI
C  R=r^2(\theta)                        for axially symmetric particles
C  DR=[dr(\theta)/(d\theta)]/r(\theta)  for axially symmetric particles
C  DDR=\lambda/[2*\pi*r(\theta)]=1/(k_out*r)
C  DRR=(MRR/(MRR**2+MRI**2))*(\lambda/[2*\pi*r(\theta)])
C                  = Re 1/(k_in*r)
C  DRI=-(MRI/(MRR**2+MRI**2))*(\lambda/[2*\pi*r(\theta)])
C                  = Im 1/(k_in*r)
C  Refractive index outside is assumed to real, whereas inside
C  a scatterer, refractive index is allowed to be complex in general.
C  Consequently, the Bessel function j_l(k_in*r) will in general
C  be complex. The routine below performs Waterman surface integral
C  separately for the real and imaginary parts of the integrand.
C
C--------/---------/---------/---------/---------/---------/---------/--
      IMPLICIT REAL*8 (A-H,O-Z)
*      INCLUDE 'ampld.par.f'
      Parameter (NPN1=100, NPNG1=2600, NPNG2=2*NPNG1, NPN2=2*NPN1,
     &           NPL=NPN2+1, NPN3=NPN1+1,
     &           NPN4=NPN1, NPN5=2*NPN4, NPN6=NPN4+1)
      INTEGER NOUT
* number of the output unit
      PARAMETER (NOUT=35)

      REAL*8  X(NPNG2),W(NPNG2),AN(NPN1),
     *        R(NPNG2),DR(NPNG2),SIG(NPN2),
     *        J(NPNG2,NPN1),Y(NPNG2,NPN1),
     *        JR(NPNG2,NPN1),JI(NPNG2,NPN1),DJ(NPNG2,NPN1),
     *        DY(NPNG2,NPN1),DJR(NPNG2,NPN1),
     *        DJI(NPNG2,NPN1),DDR(NPNG2),DRR(NPNG2),
     *        D1(NPNG2,NPN1),D2(NPNG2,NPN1),
     *        DRI(NPNG2),RR(NPNG2),
     *        DV1(NPN1),DV2(NPN1)

      REAL*8  R11(NPN1,NPN1),R12(NPN1,NPN1),
     *        R21(NPN1,NPN1),R22(NPN1,NPN1),
     *        I11(NPN1,NPN1),I12(NPN1,NPN1),
     *        I21(NPN1,NPN1),I22(NPN1,NPN1),
     *        RG11(NPN1,NPN1),RG12(NPN1,NPN1),
     *        RG21(NPN1,NPN1),RG22(NPN1,NPN1),
     *        IG11(NPN1,NPN1),IG12(NPN1,NPN1),
     *        IG21(NPN1,NPN1),IG22(NPN1,NPN1),
     *        ANN(NPN1,NPN1),
     *        QR(NPN2,NPN2),QI(NPN2,NPN2),
     *        RGQR(NPN2,NPN2),RGQI(NPN2,NPN2),
     *        TQR(NPN2,NPN2),TQI(NPN2,NPN2),
     *        TRGQR(NPN2,NPN2),TRGQI(NPN2,NPN2)
cc      REAL*8 TR1(NPN2,NPN2),TI1(NPN2,NPN2)
*
      COMMON /TMAT99/
     &            R11,R12,R21,R22,I11,I12,I21,I22,RG11,RG12,RG21,RG22,
     &            IG11,IG12,IG21,IG22           !only between TMATR routines
      COMMON /CBESS/ J,Y,JR,JI,DJ,DY,DJR,DJI
cc      COMMON /CT/ TR1,TI1                       !output from TT routine
      COMMON /CTT/ QR,QI,RGQR,RGQI              !input for TT routine
*
      MM1=1
      NNMAX=NMAX+NMAX
      NG=2*NGAUSS
      FACTOR=1D0
*
      IF (NCHECK.EQ.1) THEN          !Theta=pi/2 is scatterer mirror symmetry plane
            NGSS=NGAUSS
            FACTOR=2D0
      ELSE IF (NCHECK.EQ.0) THEN     !Theta=pi/2 is not a scatterer mirror symmetry plane
            NGSS=NG
      ENDIF
*
      SI=1D0
      DO 5 N=1,NNMAX
           SI=-SI
           SIG(N)=SI              !=(-1)**N
    5 CONTINUE
*
* Assigning Wigner d-matrices - assuming Gauss abscissas having
* mirror symmetry in the \theta=\pi/2 plane:

      DO 25 I=1,NGAUSS

         I1=NGAUSS-I+1
         I2=NGAUSS+I
*
         CALL VIG( X(I1), NMAX, 0, DV1, DV2)
*
         DO N=1,NMAX

            DD1=DV1(N)
            DD2=DV2(N)
            D1(I1,N)=DD1
            D2(I1,N)=DD2

         IF (NAXSM.EQ.1) THEN         !Gauss abscissas chosen +/- symmetric
*
* using (4.2.4) and (4.2.6) of {Ed},
*           d_{0m}^{(l)}(\pi-\theta) = (-1)^{l+m} d_{0m}^{(l)}(\theta)

            SI=SIG(N)                  !=(-1)**N
            D1(I2,N)=DD1*SI
            D2(I2,N)=-DD2*SI

         END IF
         ENDDO

         IF (NAXSM.EQ.0) THEN        !Gauss abscissas not chosen +/- symmetric
*
         CALL VIG( X(I2), NMAX, 0, DV1, DV2)
*
          DO N=1,NMAX
            DD1=DV1(N)
            DD2=DV2(N)
            D1(I2,N)=DD1
            D2(I2,N)=DD2
          ENDDO

          END IF

   25 CONTINUE
*
*  Assigning r^2(\theta)*weight product:

      DO 40 I=1,NGSS
           RR(I)=W(I)*R(I)

cc           if (dr(i).eq.0.d0) RR(I)=0.d0   !temporarily only

   40 CONTINUE
*
      DO 300 N1=MM1,NMAX
           AN1=AN(N1)
           DO 300 N2=MM1,NMAX
                AN2=AN(N2)

                AR12=0D0
                AR21=0D0
                AI12=0D0
                AI21=0D0
                GR12=0D0
                GR21=0D0
                GI12=0D0
                GI21=0D0

c        OPEN(NOUT+20,FILE='surfint.dat')   !Gauss convergence check

                IF (NCHECK.EQ.1.AND.SIG(N1+N2).LT.0D0) GO TO 205
*
* Gauss integration loop:
*
                DO 200 I=1,NGSS    !=NGAUSS   if NCHECK.EQ.1
                                   !=2*NGAUSS if NCHECK.EQ.0

                    D1N1=D1(I,N1)
                    D2N1=D2(I,N1)
                    D1N2=D1(I,N2)
                    D2N2=D2(I,N2)
                    A12=D1N1*D2N2
                    A21=D2N1*D1N2
                    A22=D2N1*D2N2
c                    AA1=A12+A21

* Vector spherical harmonics:
C  Since refractive index is allowed to be complex in general,
C  the Bessel function j_l(k_in*r) is complex. The code below
C  performs a separation of the complex integrand in Waterman's
C  surface integral into its respective real and imaginary
C  parts:

* Bessel functions of the exterior argument:

                    QJ1=J(I,N1)
                    QY1=Y(I,N1)
                    QDJ1=DJ(I,N1)
                    QDY1=DY(I,N1)

* Bessel functions of the interior argument:

                    QJR2=JR(I,N2)
                    QJI2=JI(I,N2)
                    QDJR2=DJR(I,N2)
                    QDJI2=DJI(I,N2)
*_____________________

* Re and Im of j_{n2}(k_{in}r) j_{n1}(k_{out}r):

                    C1R=QJR2*QJ1
                    C1I=QJI2*QJ1

* Re and Im of j_{n2}(k_{in}r) h_{n1}(k_{out}r):

                    B1R=C1R-QJI2*QY1
                    B1I=C1I+QJR2*QY1

* Re and Im of j_{n2}(k_{in}r) [k_{out}r j_{n1}(k_{out}r)]'/(k_{out}r):

                    C2R=QJR2*QDJ1
                    C2I=QJI2*QDJ1

* Re and Im of j_{n2}(k_{in}r) [k_{out}r h_{n1}(k_{out}r)]'/(k_{out}r):

                    B2R=C2R-QJI2*QDY1
                    B2I=C2I+QJR2*QDY1

                    DDRI=DDR(I)               !1/(k_{out}r)

* Re and Im of [1/(k_{out}r)]*j_{n2}(k_{in}r) j_{n1}(k_{out}r)

                    C3R=DDRI*C1R
                    C3I=DDRI*C1I

* Re and Im of [1/(k_{out}r)]*j_{n2}(k_{in}r) h_{n1}(k_{out}r):

                    B3R=DDRI*B1R
                    B3I=DDRI*B1I

* Re and Im of  [k_{in}r j_{n2}(k_{in}r)]'/(k_{in}r)
*                          * j_{n1}(k_{out}r):

                    C4R=QDJR2*QJ1
                    C4I=QDJI2*QJ1

* Re and Im of [k_{in}r j_{n2}(k_{in}r)]'/(k_{in}r)
*                          *  h_{n1}(k_{out}r):

                    B4R=C4R-QDJI2*QY1
                    B4I=C4I+QDJR2*QY1

                    DRRI=DRR(I)               !Re[1/(k_{in}r)]
                    DRII=DRI(I)               !Im[1/(k_{in}r)]

* Re and Im of [1/(k_{in}r)] j_{n2}(k_{in}r) j_{n1}(k_{out}r):

                    C5R=C1R*DRRI-C1I*DRII
                    C5I=C1I*DRRI+C1R*DRII


* Re and Im of [1/(k_{in}r)] j_{n2}(k_{in}r) h_{n1}(k_{out}r):

                    B5R=B1R*DRRI-B1I*DRII
                    B5I=B1I*DRRI+B1R*DRII

*%%%%%%%  Forming integrands of J-matrices (J^{11}=J^{22}=0 for m=0): %%%%%%%%

                    URI=DR(I)        !dr/(d\theta)
                    RRI=RR(I)        !w(i)*r^2(\theta)

* w(i)*r^2(\theta)*D2N1*D2N2:
                    F1=RRI*A22      !prefactor containing r^2(\theta)<->hat{r} part


* N1*(N1+1)*w(i)*r(\theta)*[dr/(d\theta)]*D1N1*D2N2:
                    F2=RRI*URI*AN1*A12     !prefactor containing r(\theta)*[dr/(d\theta)]
                                           !hat{theta} part


                    AR12=AR12+F1*B2R+F2*B3R        !~Re J^{12}
                    AI12=AI12+F1*B2I+F2*B3I        !~Im J^{12}

                    GR12=GR12+F1*C2R+F2*C3R        !~Re Rg J^{12}
                    GI12=GI12+F1*C2I+F2*C3I        !~Im Rg J^{12}

* N2*(N2+1)*w(i)*r(\theta)*[dr/(d\theta)]*D2N1*D1N2:
                    F2=RRI*URI*AN2*A21     !prefactor containing r(\theta)*[dr/(d\theta)]
                                           !hat{theta} part

                    AR21=AR21+F1*B4R+F2*B5R        !~Re J^{21}
                    AI21=AI21+F1*B4I+F2*B5I        !~Im J^{21}

                    GR21=GR21+F1*C4R+F2*C5R        !~Re Rg J^{21}
                    GI21=GI21+F1*C4I+F2*C5I        !~Im Rg J^{21}

  200           CONTINUE               !end of Gauss integration

c                write(nout+20,*)'N1=',N1,'   N2=',N2
c                write(nout+20,*)'AR12=', AR12
c                write(nout+20,*)'AI12=', AI12
c                write(nout+20,*)'AR21=', AR21
c                write(nout+20,*)'AI21=', AI21
c                write(nout+20,*)'GR12=', GR12
c                write(nout+20,*)'GI12=', GI12
c                write(nout+20,*)'GR21=', GR21
c                write(nout+20,*)'GI21=', GI21

*%%%%%%%%%%%%%  Forming J-matrices (J^{11}=J^{22}=0 for m=0):

  205           AN12=ANN(N1,N2)*FACTOR

                R12(N1,N2)=AR12*AN12
                R21(N1,N2)=AR21*AN12
                I12(N1,N2)=AI12*AN12
                I21(N1,N2)=AI21*AN12

                RG12(N1,N2)=GR12*AN12
                RG21(N1,N2)=GR21*AN12
                IG12(N1,N2)=GI12*AN12
                IG21(N1,N2)=GI21*AN12

  300 CONTINUE            !end of the loop over angular momenta

c      close(nout+20)

*%%%%%%%%%%%%%%%%%%%%%%%  Forming Q and RgQ -matrices

      TPIR=PIR                 !Re [1/k_{in}^2]
      TPII=PII                 !Im [1/k_{in}^2]
      TPPI=PPI                 !1/k_{out}^2

      NM=NMAX

      DO 310 N1=MM1,NMAX
           K1=N1-MM1+1
           KK1=K1+NM
           DO 310 N2=MM1,NMAX
                K2=N2-MM1+1
                KK2=K2+NM

                TAR12= I12(N1,N2)
                TAI12=-R12(N1,N2)
                TGR12= IG12(N1,N2)
                TGI12=-RG12(N1,N2)

                TAR21=-I21(N1,N2)
                TAI21= R21(N1,N2)
                TGR21=-IG21(N1,N2)
                TGI21= RG21(N1,N2)

                TQR(K1,K2)=TPIR*TAR21-TPII*TAI21+TPPI*TAR12
                TQI(K1,K2)=TPIR*TAI21+TPII*TAR21+TPPI*TAI12
                TRGQR(K1,K2)=TPIR*TGR21-TPII*TGI21+TPPI*TGR12
                TRGQI(K1,K2)=TPIR*TGI21+TPII*TGR21+TPPI*TGI12

                TQR(K1,KK2)=0D0
                TQI(K1,KK2)=0D0
                TRGQR(K1,KK2)=0D0
                TRGQI(K1,KK2)=0D0

                TQR(KK1,K2)=0D0
                TQI(KK1,K2)=0D0
                TRGQR(KK1,K2)=0D0
                TRGQI(KK1,K2)=0D0

                TQR(KK1,KK2)=TPIR*TAR12-TPII*TAI12+TPPI*TAR21
                TQI(KK1,KK2)=TPIR*TAI12+TPII*TAR12+TPPI*TAI21
                TRGQR(KK1,KK2)=TPIR*TGR12-TPII*TGI12+TPPI*TGR21
                TRGQI(KK1,KK2)=TPIR*TGI12+TPII*TGR12+TPPI*TGI21
  310 CONTINUE

      NNMAX=2*NM
      DO 320 N1=1,NNMAX
           DO 320 N2=1,NNMAX
                QR(N1,N2)=TQR(N1,N2)
                QI(N1,N2)=TQI(N1,N2)
                RGQR(N1,N2)=TRGQR(N1,N2)
                RGQI(N1,N2)=TRGQI(N1,N2)
  320 CONTINUE

*%%%%%%%%%%%%%%%%%%%%%%%  Forming resulting T-matrix
*
* Calculate the product Q^{-1} Rg Q
*
      CALL TT(NMAX,NCHECK)
*
      RETURN
      END

C**********************************************************************

      SUBROUTINE TMATR(M,NGAUSS,X,W,AN,ANN,S,SS,PPI,PIR,PII,R,DR,DDR,
     *                  DRR,DRI,NMAX,NCHECK,NAXSM)
C--------/---------/---------/---------/---------/---------/---------/--
C >>> M,NGAUSS,X,W,AN,ANN,S,SS,PPI,PIR,PII,R,DR,DDR,DRR,DRI,NMAX,NCHECK
C <<< common blocks /TMAT99/, /CT/ (for main),  and /CTT/ (for TT)
C=====================
C
C  Determines the T-matrix of an axially symmetric scatterer
C                           for M.GT.0
C
C  M      - azimuthal number
C  NGAUSS - the number of GIF division points
C  X=\cos\theta  - GIF division points
C  W - GIF weights
C  AN(N)=N*(N+1)
C  ANN(l_1,l_2)=\sqrt{\fr{(2 l_1+1)}{l_1(l_1+1)} }
C                       \sqrt{\fr{(2 l_2+1)}{l_2(l_2+1)} }/2
C  NMAX - angular momentum cutoff
C  NCHECK  -  .EQ.0  THEN  NGSS=2*NGAUSS, FACTOR=1D0
C             .EQ.1  THEN  NGSS = NGAUSS, FACTOR=2D0
C  NAXSM   -  .EQ.0 : Gauss abscissas do not have +/- theta symmetry
C             .EQ.1 : Gauss abscissas have +/- theta symmetry
C  NCHECK - specifies whether NG=2*NGAUSS or otherwise
C  P=DACOS(-1D0)
C  PI=P*2D0/LAM - wave vector
C  PPI=PI*PI
C  PIR=PPI*MRR
C  PII=PPI*MRI
C  R=r^2(\theta)                       for axially symmetric particles
C  DR=[dr(\theta)/(d\theta)]/r(\theta) for axially symmetric particles
C  DDR=\lambda/[2*\pi*r(\theta)]=1/(k_out*r)
C  DRR=(MRR/(MRR**2+MRI**2))*(\lambda/[2*\pi*r(\theta)])
C                  = Re 1/(k_in*r)
C  DRI=-(MRI/(MRR**2+MRI**2))*(\lambda/[2*\pi*r(\theta)])
C                  = Im 1/(k_in*r)
C  S(I)=1/(|\sin\theta|)
C  SS(I)=1/(\sin^2\theta)   (supplied from CONSTE)
C
C  Refractive index outside is assumed to real, whereas inside
C  a scatterer, refractive index is allowed to be complex in general.
C  Consequently, the Bessel function j_l(k_in*r) will in general
C  be complex. The routine below performs Waterman surface integral
C  separately for the real and imaginary parts of the integrand.
C
C--------/---------/---------/---------/---------/---------/---------/--
*      INCLUDE 'ampld.par.f'
      Parameter (NPN1=100, NPNG1=2600, NPNG2=2*NPNG1, NPN2=2*NPN1,
     &           NPL=NPN2+1, NPN3=NPN1+1,
     &           NPN4=NPN1, NPN5=2*NPN4, NPN6=NPN4+1)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8  X(NPNG2),W(NPNG2),AN(NPN1),S(NPNG2),SS(NPNG2),
     *        R(NPNG2),DR(NPNG2),SIG(NPN2),
     *        J(NPNG2,NPN1),Y(NPNG2,NPN1),
     *        JR(NPNG2,NPN1),JI(NPNG2,NPN1),DJ(NPNG2,NPN1),
     *        DY(NPNG2,NPN1),DJR(NPNG2,NPN1),
     *        DJI(NPNG2,NPN1),DDR(NPNG2),DRR(NPNG2),
     *        D1(NPNG2,NPN1),D2(NPNG2,NPN1),
     *        DRI(NPNG2),DS(NPNG2),DSS(NPNG2),RR(NPNG2),
     *        DV1(NPN1),DV2(NPN1)

      REAL*8  R11(NPN1,NPN1),R12(NPN1,NPN1),
     *        R21(NPN1,NPN1),R22(NPN1,NPN1),
     *        I11(NPN1,NPN1),I12(NPN1,NPN1),
     *        I21(NPN1,NPN1),I22(NPN1,NPN1),
     *        RG11(NPN1,NPN1),RG12(NPN1,NPN1),
     *        RG21(NPN1,NPN1),RG22(NPN1,NPN1),
     *        IG11(NPN1,NPN1),IG12(NPN1,NPN1),
     *        IG21(NPN1,NPN1),IG22(NPN1,NPN1),
     *        ANN(NPN1,NPN1),
     *        QR(NPN2,NPN2),QI(NPN2,NPN2),
     *        RGQR(NPN2,NPN2),RGQI(NPN2,NPN2),
     *        TQR(NPN2,NPN2),TQI(NPN2,NPN2),
     *        TRGQR(NPN2,NPN2),TRGQI(NPN2,NPN2)
cc      REAL*8 TR1(NPN2,NPN2),TI1(NPN2,NPN2)
*________
      COMMON /TMAT99/
     &            R11,R12,R21,R22,I11,I12,I21,I22,RG11,RG12,RG21,RG22,
     &            IG11,IG12,IG21,IG22          !only between TMATR routines
      COMMON /CBESS/ J,Y,JR,JI,DJ,DY,DJR,DJI
cc      COMMON /CT/ TR1,TI1                      !output from TT routine
      COMMON /CTT/ QR,QI,RGQR,RGQI             !input for TT routine
*________
      MM1=M
      QM=DFLOAT(M)
      QMM=QM*QM
      NG=2*NGAUSS
      NM=NMAX+NMAX
      FACTOR=1D0
*
      IF (NCHECK.EQ.1) THEN          !Theta=pi/2 is scatterer mirror symmetry plane
            NGSS=NGAUSS
            FACTOR=2D0
      ELSE IF (NCHECK.EQ.0) THEN     !Theta=pi/2 is not a scatterer mirror symmetry plane
            NGSS=NG
      ENDIF
*
      SI=1D0
      DO 5 N=1,NM                 !NM=2*NMAX
           SI=-SI
           SIG(N)=SI              !=(-1)**N
    5 CONTINUE
*
* Assigning Wigner d-matrices:

      DO 25 I=1,NGAUSS

         I1=NGAUSS-I+1
         I2=NGAUSS+I
*
         CALL VIG( X(I1), NMAX, M, DV1, DV2)
C--------/---------/---------/---------/---------/---------/---------/--
C     DV1(N)=dvig(0,m,n,arccos x)/sin(arccos x)   != d_{0m}^{(l)}
C     DV2(N)=[d/d(arccos x)] dvig(0,m,n,arccos x) != d d_{0m}^{(l)}/d\theta
C--------/---------/---------/---------/---------/---------/---------/--
*
         DO N=1,NMAX

            DD1=DV1(N)
            DD2=DV2(N)
            D1(I1,N)=DD1
            D2(I1,N)=DD2

         IF (NAXSM.EQ.1) THEN         !Gauss abscissas chosen +/- symmetric
*
* using (4.2.4) and (4.2.6) of {Ed},
*           d_{0m}^{(l)}(\pi-\theta) = (-1)^{l+m} d_{0m}^{(l)}(\theta)

            SI=SIG(N+M)           !=(-1)**(N+M)
                                  !exactly what follows from {Ed}
            D1(I2,N)=DD1*SI
            D2(I2,N)=-DD2*SI

         END IF
         ENDDO

         IF (NAXSM.EQ.0) THEN        !Gauss abscissas not chosen +/- symmetric
*
         CALL VIG( X(I2), NMAX, M, DV1, DV2)
*
          DO N=1,NMAX
            DD1=DV1(N)
            DD2=DV2(N)
            D1(I2,N)=DD1
            D2(I2,N)=DD2
          ENDDO

          END IF

   25 CONTINUE
*
*  Assigning r^2(\theta)*weight product:

      DO 40 I=1,NGSS
           WR=W(I)*R(I)

cc          if (dr(i).eq.0.d0) WR=0.d0   !temporarily only

           DS(I)=S(I)*QM*WR       !=DFLOAT(M)*W(I)*r^2(\theta)/(|\sin\theta|)
           DSS(I)=SS(I)*QMM       !=DFLOAT(M)**2/(\sin^2\theta)
           RR(I)=WR

   40 CONTINUE
*
      DO 300  N1=MM1,NMAX         !MM1=M below
           AN1=AN(N1)

           DO 300 N2=MM1,NMAX
                AN2=AN(N2)
                AR11=0D0
                AR12=0D0
                AR21=0D0
                AR22=0D0
                AI11=0D0
                AI12=0D0
                AI21=0D0
                AI22=0D0
                GR11=0D0
                GR12=0D0
                GR21=0D0
                GR22=0D0
                GI11=0D0
                GI12=0D0
                GI21=0D0
                GI22=0D0
                SI=SIG(N1+N2)

                DO 200 I=1,NGSS    !=NGAUSS   if NCHECK.EQ.1
                                   !=2*NGAUSS if NCHECK.EQ.0
                    D1N1=D1(I,N1)
                    D2N1=D2(I,N1)
                    D1N2=D1(I,N2)
                    D2N2=D2(I,N2)
                    A11=D1N1*D1N2
                    A12=D1N1*D2N2
                    A21=D2N1*D1N2
                    A22=D2N1*D2N2
                    AA1=A12+A21            != D1N1*D2N2+D2N1*D1N2
                    AA2=A11*DSS(I)+A22     !=(D1N1*D1N2)*DFLOAT(M)**2/(\sin^2\theta)
                                           ! +D2N1*D2N2

* Vector spherical harmonics:
C  Since refractive index is allowed to be complex in general,
C  the Bessel function j_l(k_in*r) is complex. The code below
C  performs a separation of the complex integrand in Waterman's
C  surface integral into its respective real and imaginary
C  parts:

* Bessel functions of the exterior argument:

                    QJ1=J(I,N1)
                    QY1=Y(I,N1)
                    QDJ1=DJ(I,N1)
                    QDY1=DY(I,N1)

* Bessel functions of the interior argument:

                    QJR2=JR(I,N2)
                    QJI2=JI(I,N2)
                    QDJR2=DJR(I,N2)
                    QDJI2=DJI(I,N2)

* Re and Im of j_{n2}(k_{in}r) j_{n1}(k_{out}r):

                    C1R=QJR2*QJ1
                    C1I=QJI2*QJ1

* Re and Im of j_{n2}(k_{in}r) h_{n1}(k_{out}r):

                    B1R=C1R-QJI2*QY1
                    B1I=C1I+QJR2*QY1

* Re and Im of j_{n2}(k_{in}r) j_{n1}'(k_{out}r)/(k_{out}r):

                    C2R=QJR2*QDJ1
                    C2I=QJI2*QDJ1

* Re and Im of j_{n2}(k_{in}r) h_{n1}'(k_{out}r)/(k_{out}r):

                    B2R=C2R-QJI2*QDY1
                    B2I=C2I+QJR2*QDY1

                    DDRI=DDR(I)               !1/(k_{out}r)

* Re and Im of [1/(k_{out}r)]*j_{n2}(k_{in}r) j_{n1}(k_{out}r)

                    C3R=DDRI*C1R
                    C3I=DDRI*C1I

* Re and Im of [1/(k_{out}r)]*j_{n2}(k_{in}r) h_{n1}(k_{out}r):

                    B3R=DDRI*B1R
                    B3I=DDRI*B1I

* Re and Im of [k_{in}r j_{n2}(k_{in}r)]'/(k_{in}r)
*                          * j_{n1}(k_{out}r):

                    C4R=QDJR2*QJ1
                    C4I=QDJI2*QJ1

* Re and Im of [k_{in}r j_{n2}(k_{in}r)]'/(k_{in}r)
*                          *  h_{n1}(k_{out}r):

                    B4R=C4R-QDJI2*QY1
                    B4I=C4I+QDJR2*QY1

                    DRRI=DRR(I)               !Re[1/(k_{in}r)]
                    DRII=DRI(I)               !Im[1/(k_{in}r)]

* Re and Im of [1/(k_{in}r)] j_{n2}(k_{in}r) j_{n1}(k_{out}r):

                    C5R=C1R*DRRI-C1I*DRII
                    C5I=C1I*DRRI+C1R*DRII

* Re and Im of [1/(k_{in}r)] j_{n2}(k_{in}r) h_{n1}(k_{out}r):

                    B5R=B1R*DRRI-B1I*DRII
                    B5I=B1I*DRRI+B1R*DRII


* Re and Im of j_{n2}'(k_{in}r) j_{n1}'(k_{out}r):

                    C6R=QDJR2*QDJ1
                    C6I=QDJI2*QDJ1

* Re and Im of j_{n2}'(k_{in}r) h_{n1}'(k_{out}r):

                    B6R=C6R-QDJI2*QDY1
                    B6I=C6I+QDJR2*QDY1

* Re and Im of [1/(k_{out}r)] j_{n2}'(k_{in}r) j_{n1}(k_{out}r):

                    C7R=C4R*DDRI
                    C7I=C4I*DDRI

* Re and Im of [1/(k_{out}r)] j_{n2}'(k_{in}r) h_{n1}(k_{out}r):

                    B7R=B4R*DDRI
                    B7I=B4I*DDRI

* Re and Im of [1/(k_{in}r)] j_{n2}(k_{in}r) j_{n1}'(k_{out}r):

                    C8R=C2R*DRRI-C2I*DRII
                    C8I=C2I*DRRI+C2R*DRII

* Re and Im of [1/(k_{in}r)] j_{n2}(k_{in}r) h_{n1}'(k_{out}r):

                    B8R=B2R*DRRI-B2I*DRII
                    B8I=B2I*DRRI+B2R*DRII


* %%%%%%%%%  Forming integrands of J-matrices (J^{11}=J^{22}=0 for m.eq.0):

                    URI=DR(I)
                    DSI=DS(I)          !prop to  m/sin(theta)
                    RRI=RR(I)

                    IF (NCHECK.EQ.1.AND.SI.GT.0D0) GO TO 150

* [DFLOAT(M)*W(I)*r^2(I)/(|\sin\theta|)]*(D1N1*D2N2+D2N1*D1N2):

                    E1=DSI*AA1             ! <-- AA1

                    AR11=AR11+E1*B1R
                    AI11=AI11+E1*B1I
                    GR11=GR11+E1*C1R
                    GI11=GI11+E1*C1I

                    IF (NCHECK.EQ.1) GO TO 160

  150               CONTINUE


* w(i)*r^2(\theta)*[(D1N1*D1N2)*DFLOAT(M)**2/(\sin^2\theta)+D2N1*D2N2]
* (prefactor containing r^2(\theta)<->hat{r} part)

                    F1=RRI*AA2             ! <-- AA2

* N1*(N1+1)*w(i)*r(\theta)*[dr/(d\theta)]*D1N1*D2N2:
*  (prefactor containing r(\theta)*[dr/(d\theta)] - hat{theta} part)

                    F2=RRI*URI*AN1*A12             ! <-- A12

                    AR12=AR12+F1*B2R+F2*B3R        !~Re J^{12}
                    AI12=AI12+F1*B2I+F2*B3I        !~Im J^{12}

                    GR12=GR12+F1*C2R+F2*C3R        !~Re Rg J^{12}
                    GI12=GI12+F1*C2I+F2*C3I        !~Im Rg J^{12}

* N2*(N2+1)*w(i)*r(\theta)*[dr/(d\theta)]*D2N1*D1N2:
* (!prefactor containing r(\theta)*[dr/(d\theta)] - hat{theta} part)

                    F2=RRI*URI*AN2*A21             ! <-- A21

                    AR21=AR21+F1*B4R+F2*B5R
                    AI21=AI21+F1*B4I+F2*B5I

                    GR21=GR21+F1*C4R+F2*C5R
                    GI21=GI21+F1*C4I+F2*C5I

                    IF (NCHECK.EQ.1) GO TO 200

  160               E2=DSI*URI*A11
                    E3=E2*AN2
                    E2=E2*AN1

                    AR22=AR22+E1*B6R+E2*B7R+E3*B8R
                    AI22=AI22+E1*B6I+E2*B7I+E3*B8I

                    GR22=GR22+E1*C6R+E2*C7R+E3*C8R
                    GI22=GI22+E1*C6I+E2*C7I+E3*C8I

  200           CONTINUE           !Gauss integration

*%%%%%%%%%%%%%  Forming J-matrices (J^{11}=J^{22}=0 for m.eq.0):

                AN12=ANN(N1,N2)*FACTOR

                R11(N1,N2)=AR11*AN12       !Re J^{11}
                R12(N1,N2)=AR12*AN12       !Re J^{12}
                R21(N1,N2)=AR21*AN12       !Re J^{21}
                R22(N1,N2)=AR22*AN12       !Re J^{22}
                I11(N1,N2)=AI11*AN12       !Im J^{11}
                I12(N1,N2)=AI12*AN12       !Im J^{12}
                I21(N1,N2)=AI21*AN12       !Im J^{21}
                I22(N1,N2)=AI22*AN12       !Im J^{22}

                RG11(N1,N2)=GR11*AN12       !Re (Rg J^{11})
                RG12(N1,N2)=GR12*AN12       !Re (Rg J^{12})
                RG21(N1,N2)=GR21*AN12       !Re (Rg J^{21})
                RG22(N1,N2)=GR22*AN12       !Re (Rg J^{22})
                IG11(N1,N2)=GI11*AN12       !Im (Rg J^{11})
                IG12(N1,N2)=GI12*AN12       !Im (Rg J^{12})
                IG21(N1,N2)=GI21*AN12       !Im (Rg J^{21})
                IG22(N1,N2)=GI22*AN12       !Im (Rg J^{22})

  300 CONTINUE

*%%%%%%%%%%%%%%%%%%%%%%%  Forming Q and RgQ -matrices

      TPIR=PIR                 !Re [1/k_{in}^2]
      TPII=PII                 !Im [1/k_{in}^2]
      TPPI=PPI                 !1/k_{out}^2

      NM=NMAX-MM1+1
      DO 310 N1=MM1,NMAX
           K1=N1-MM1+1                !from 1 to NMAX-MM1+1
           KK1=K1+NM                  !from NMAX-MM1+2 to 2*(NMAX-MM1+1)

           DO 310 N2=MM1,NMAX
                K2=N2-MM1+1           !from 1 to NMAX-MM1+1
                KK2=K2+NM             !from NMAX-MM1+2 to 2*(NMAX-MM1+1)

                TAR11=-R11(N1,N2)
                TAI11=-I11(N1,N2)
                TGR11=-RG11(N1,N2)
                TGI11=-IG11(N1,N2)

                TAR12= I12(N1,N2)
                TAI12=-R12(N1,N2)
                TGR12= IG12(N1,N2)
                TGI12=-RG12(N1,N2)

                TAR21=-I21(N1,N2)
                TAI21= R21(N1,N2)
                TGR21=-IG21(N1,N2)
                TGI21= RG21(N1,N2)

                TAR22=-R22(N1,N2)
                TAI22=-I22(N1,N2)
                TGR22=-RG22(N1,N2)
                TGI22=-IG22(N1,N2)

                TQR(K1,K2)=TPIR*TAR21-TPII*TAI21+TPPI*TAR12
                TQI(K1,K2)=TPIR*TAI21+TPII*TAR21+TPPI*TAI12
                TRGQR(K1,K2)=TPIR*TGR21-TPII*TGI21+TPPI*TGR12
                TRGQI(K1,K2)=TPIR*TGI21+TPII*TGR21+TPPI*TGI12

                TQR(K1,KK2)=TPIR*TAR11-TPII*TAI11+TPPI*TAR22
                TQI(K1,KK2)=TPIR*TAI11+TPII*TAR11+TPPI*TAI22
                TRGQR(K1,KK2)=TPIR*TGR11-TPII*TGI11+TPPI*TGR22
                TRGQI(K1,KK2)=TPIR*TGI11+TPII*TGR11+TPPI*TGI22

                TQR(KK1,K2)=TPIR*TAR22-TPII*TAI22+TPPI*TAR11
                TQI(KK1,K2)=TPIR*TAI22+TPII*TAR22+TPPI*TAI11
                TRGQR(KK1,K2)=TPIR*TGR22-TPII*TGI22+TPPI*TGR11
                TRGQI(KK1,K2)=TPIR*TGI22+TPII*TGR22+TPPI*TGI11

                TQR(KK1,KK2)=TPIR*TAR12-TPII*TAI12+TPPI*TAR21
                TQI(KK1,KK2)=TPIR*TAI12+TPII*TAR12+TPPI*TAI21
                TRGQR(KK1,KK2)=TPIR*TGR12-TPII*TGI12+TPPI*TGR21
                TRGQI(KK1,KK2)=TPIR*TGI12+TPII*TGR12+TPPI*TGI21

  310 CONTINUE

      NNMAX=2*NM
      DO 320 N1=1,NNMAX
           DO 320 N2=1,NNMAX
                QR(N1,N2)=TQR(N1,N2)
                QI(N1,N2)=TQI(N1,N2)
                RGQR(N1,N2)=TRGQR(N1,N2)
                RGQI(N1,N2)=TRGQI(N1,N2)
  320 CONTINUE
*
      CALL TT(NM,NCHECK)
*
      RETURN
      END


C*****************************************************************

      SUBROUTINE VIG(X, NMAX, M, DV1, DV2)
C--------/---------/---------/---------/---------/---------/---------/--
C >>> X,NMAX,M
C <<< DV1, DV2
C =============
C     For a given azimuthal number M, calculation of the Wigner d-functions
C     DV1(N)=dvig(0,m,n,arccos x)/sin(arccos x)   != d_{0m}^{(l)}
C     and
C     DV2(N)=[d/d(arccos x)] dvig(0,m,n,arccos x) != d d_{0m}^{(l)}/d\theta
C     for 1.LE.N.LE.NMAX and 0.LE.X.LE.1.
C     For M.NEQ.0, only the  M.LE.N.LE.NMAX terms are determined
C
C     Made using recurrences of  Ref. \ct{Mis39}
C     (There is a missing $l$ factor in the 2nd term in the curly bracket
C     in recurrence (35) of Ref. \ct{Mis39} for DV2).
C
C     X=cos(theta), where theta is the polar angle
C     NMAX - angular momentum cutoff
C
C--------/---------/---------/---------/---------/---------/---------/--
*      INCLUDE 'ampld.par.f'
      Parameter (NPN1=100, NPNG1=2600, NPNG2=2*NPNG1, NPN2=2*NPN1,
     &           NPL=NPN2+1, NPN3=NPN1+1,
     &           NPN4=NPN1, NPN5=2*NPN4, NPN6=NPN4+1)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 DV1(NPN1),DV2(NPN1)

      A=1D0
      QS=DSQRT(1D0-X*X)
      QS1=1D0/QS
      DO N=1,NMAX
         DV1(N)=0D0
         DV2(N)=0D0
      ENDDO

      IF (M.NE.0) GO TO 20

      D1=1D0
      D2=X
      DO N=1,NMAX
         QN=DFLOAT(N)
         QN1=DFLOAT(N+1)
         QN2=DFLOAT(2*N+1)
         D3=(QN2*X*D2-QN*D1)/QN1          !recurrence (31) of Ref. {Mis39}
         DER=QS1*(QN1*QN/QN2)*(-D1+D3)    !recurrence (35) of Ref. {Mis39}
         DV1(N)=D2
         DV2(N)=DER
         D1=D2
         D2=D3
      ENDDO
      RETURN

   20 QMM=DFLOAT(M*M)

*A_m initialization - recurrence (34) of Ref. {Mis39}
      DO I=1,M
         I2=I*2
         A=A*DSQRT(DFLOAT(I2-1)/DFLOAT(I2))*QS
      ENDDO
*
      D1=0D0
      D2=A

      DO N=M,NMAX
         QN=DFLOAT(N)
         QN2=DFLOAT(2*N+1)
         QN1=DFLOAT(N+1)
         QNM=DSQRT(QN*QN-QMM)
         QNM1=DSQRT(QN1*QN1-QMM)
         D3=(QN2*X*D2-QNM*D1)/QNM1              !recurrence (31) of Ref. {Mis39}
         DER=QS1*(-QN1*QNM*D1+QN*QNM1*D3)/QN2   !recurrence (35) of Ref. {Mis39}
         DV1(N)=D2
         DV2(N)=DER
         D1=D2
         D2=D3
      ENDDO
*
      RETURN
      END

C**********************************************************************

      SUBROUTINE TT(NMAX,NCHECK)
C--------/---------/---------/---------/---------/---------/---------/--
C >>> NMAX,NCHECK
C <<< COMMON BLOCKS
C=================
C  NMAX=NMAX-M+1 here, where NMAX is the angular momentum cutoff in main
C  NCHECK -
C
C   CALCULATION OF THE MATRIX    T = - RG(Q) * (Q**(-1))
C
C   INPUT  IN COMMON /CTT/
C   OUTPUT IN COMMON /CT/
C
C--------/---------/---------/---------/---------/---------/---------/--
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER NOUT

* number of the output unit
      PARAMETER (NOUT=35)
*      INCLUDE 'ampld.par.f'
      Parameter (NPN1=100, NPNG1=2600, NPNG2=2*NPNG1, NPN2=2*NPN1,
     &           NPL=NPN2+1, NPN3=NPN1+1,
     &           NPN4=NPN1, NPN5=2*NPN4, NPN6=NPN4+1)
      REAL*8  QR(NPN2,NPN2),QI(NPN2,NPN2),EMACH,
     *       RGQR(NPN2,NPN2),RGQI(NPN2,NPN2)
cc      REAL*8 F(NPN2,NPN2),B(NPN2),WORK(NPN2),
cc     *       A(NPN2,NPN2),C(NPN2,NPN2),D(NPN2,NPN2),E(NPN2,NPN2)
      REAL*8 TR1(NPN2,NPN2),TI1(NPN2,NPN2)
      COMPLEX*16 ZQ(NPN2,NPN2),ZX(NPN2),ZW(NPN2)
      INTEGER IPIV(NPN2),IPVT(NPN2)
*
      COMMON /CHOICE/ ICHOICE
      COMMON /CT/ TR1,TI1
      COMMON /CTT/ QR,QI,RGQR,RGQI
*
      DATA EMACH/0.D0/                 !/1.D-17/
*
      NNMAX=2*NMAX

      DO I=1,NNMAX
       DO J=1,NNMAX
          ZQ(I,J)=DCMPLX(QR(I,J),QI(I,J))
       ENDDO
      ENDDO

      IF (ICHOICE.EQ.2) GOTO 5    ! NAG or not NAG decision tree

********************************************************************
*   Inversion from NAG-LIB or Waterman's method    !NAG library used
*
       INFO=0
*
cc           CALL F07ARF(NNMAX,NNMAX,ZQ,NPN2,IPIV,INFO)
cc           IF (INFO.NE.0) WRITE(NOUT,1100) INFO
cc           CALL F07AWF(NNMAX,ZQ,NPN2,IPIV,ZX,NPN2,INFO)
cc           IF (INFO.NE.0) WRITE(NOUT,1100) INFO
*
 1100      FORMAT ('WARNING:  info=', i2)

* Calculate T-matrix = - RG(Q) * (Q**(-1))
*
       DO I=1,NNMAX
          DO J=1,NNMAX
             TR=0D0
             TI=0D0
             DO K=1,NNMAX
                    ARR=RGQR(I,K)
                    ARI=RGQI(I,K)
                    AR=ZQ(K,J)
                    AI=DIMAG(ZQ(K,J))
                    TR=TR-ARR*AR+ARI*AI
                    TI=TI-ARR*AI-ARI*AR
                 ENDDO
             TR1(I,J)=TR
             TI1(I,J)=TI
          ENDDO
       ENDDO

       GOTO 70                        !Return

*********************************************************************

C  Gaussian elimination             !NAG library not used

  5   CALL ZGER(ZQ,IPIV,NNMAX,NPN2,EMACH)  !Gauss elimination of ZQ to
                                           !a lower diagonal matrix
      DO 6 I=1,NNMAX
              DO K=1,NNMAX    !Initialization of the right-hand side ZB
                              !(a row vector) of the matrix equation ZX*ZQ=ZB

              ZX(K)=DCMPLX(RGQR(I,K),RGQI(I,K))
              ENDDO

      CALL ZSUR(ZQ,IPIV,ZX,NNMAX,NPN2,EMACH)  !Solving ZX*ZQ=ZB by
                                               !backsubstition
                                               !(ZX overwritten on exit)
             DO K=1,NNMAX
*
* Assign T-matrix elements = - RG(Q) * (Q**(-1))
*
             TR1(I,K)=-DBLE(ZX(K))
             TI1(I,K)=-DIMAG(ZX(K))
             ENDDO
  6   CONTINUE
*
   70 RETURN
      END

C********************************************************************

      SUBROUTINE PROD(A,B,C,NDIM,N)
C--------/---------/---------/---------/---------/---------/---------/--
C >>> A,B,NDIM,N
C <<< C=A*B
C=================
C--------/---------/---------/---------/---------/---------/---------/--
      REAL*8 A(NDIM,N),B(NDIM,N),C(NDIM,N),cij
*
      DO 10 I=1,N
           DO 10 J=1,N
                CIJ=0d0
                DO 5 K=1,N
                     CIJ=CIJ+A(I,K)*B(K,J)
    5           CONTINUE
                C(I,J)=CIJ
   10 CONTINUE
*
      RETURN
      END

C**********************************************************************

      SUBROUTINE INV1 (NMAX,F,A)
C--------/---------/---------/---------/---------/---------/---------/--
C >>> EPS
C <<< RAT
C=================
C  NMAX - angular momentum cutoff
C--------/---------/---------/---------/---------/---------/---------/--
      IMPLICIT REAL*8 (A-H,O-Z)
*      INCLUDE 'ampld.par.f'
      Parameter (NPN1=100, NPNG1=2600, NPNG2=2*NPNG1, NPN2=2*NPN1,
     &           NPL=NPN2+1, NPN3=NPN1+1,
     &           NPN4=NPN1, NPN5=2*NPN4, NPN6=NPN4+1)
      REAL*8  A(NPN2,NPN2),F(NPN2,NPN2),B(NPN1),
     *        WORK(NPN1),Q1(NPN1,NPN1),Q2(NPN1,NPN1),
     &        P1(NPN1,NPN1),P2(NPN1,NPN1)
      INTEGER IPVT(NPN1),IND1(NPN1),IND2(NPN1)
*
      NDIM=NPN1
      NN1=(DFLOAT(NMAX)-0.1D0)*0.5D0+1D0
      NN2=NMAX-NN1
*
      DO 5 I=1,NMAX
         IND1(I)=2*I-1
         IF(I.GT.NN1) IND1(I)=NMAX+2*(I-NN1)
         IND2(I)=2*I
         IF(I.GT.NN2) IND2(I)=NMAX+2*(I-NN2)-1
    5 CONTINUE
      NNMAX=2*NMAX
*
      DO 15 I=1,NMAX
         I1=IND1(I)
         I2=IND2(I)
         DO 15 J=1,NMAX
            J1=IND1(J)
            J2=IND2(J)
            Q1(J,I)=F(J1,I1)
            Q2(J,I)=F(J2,I2)
   15 CONTINUE
*
      CALL INVERT(NDIM,NMAX,Q1,P1,COND,IPVT,WORK,B)
      CALL INVERT(NDIM,NMAX,Q2,P2,COND,IPVT,WORK,B)
*
      DO 30 I=1,NNMAX
         DO 30 J=1,NNMAX
            A(J,I)=0D0
   30 CONTINUE
      DO 40 I=1,NMAX
         I1=IND1(I)
         I2=IND2(I)
         DO 40 J=1,NMAX
            J1=IND1(J)
            J2=IND2(J)
            A(J1,I1)=P1(J,I)
            A(J2,I2)=P2(J,I)
   40 CONTINUE
*
      RETURN
      END

C*********************************************************************

      SUBROUTINE INVERT(NDIM,N,A,X,COND,IPVT,WORK,B)
C--------/---------/---------/---------/---------/---------/---------/--
C >>> EPS
C <<< RAT
C=================
C--------/---------/---------/---------/---------/---------/---------/--
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 A(NDIM,N),X(NDIM,N),WORK(N),B(N)
      INTEGER IPVT(N)
*
      CALL DECOMP(NDIM,N,A,COND,IPVT,WORK)
*
      IF (COND+1D0.EQ.COND) PRINT 5,COND
C     IF (COND+1D0.EQ.COND) STOP
    5 FORMAT(' THE MATRIX IS SINGULAR FOR THE GIVEN NUMERICAL ACCURACY '
     *      ,'COND = ',D12.6)

      DO 30 I=1,N
           DO 10 J=1,N
                B(J)=0D0
                IF (J.EQ.I) B(J)=1D0
  10       CONTINUE
*
           CALL SOLVE (NDIM,N,A,B,IPVT)
*
           DO 30 J=1,N
                X(J,I)=B(J)
   30 CONTINUE
*
      RETURN
      END

C********************************************************************

      SUBROUTINE DECOMP(NDIM,N,A,COND,IPVT,WORK)
C--------/---------/---------/---------/---------/---------/---------/--
C >>> EPS
C <<< RAT
C=================
C--------/---------/---------/---------/---------/---------/---------/--
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 A(NDIM,N),COND,WORK(N)
      INTEGER IPVT(N)
*
      IPVT(N)=1
      IF(N.EQ.1) GO TO 80
      NM1=N-1
      ANORM=0D0
      DO 10 J=1,N
          T=0D0
          DO 5 I=1,N
              T=T+DABS(A(I,J))
    5     CONTINUE
          IF (T.GT.ANORM) ANORM=T
   10 CONTINUE
      DO 35 K=1,NM1
          KP1=K+1
          M=K
          DO 15 I=KP1,N
              IF (DABS(A(I,K)).GT.DABS(A(M,K))) M=I
   15     CONTINUE
          IPVT(K)=M
          IF (M.NE.K) IPVT(N)=-IPVT(N)
          T=A(M,K)
          A(M,K)=A(K,K)
          A(K,K)=T
          IF (T.EQ.0d0) GO TO 35
          DO 20 I=KP1,N
              A(I,K)=-A(I,K)/T
   20     CONTINUE
          DO 30 J=KP1,N
              T=A(M,J)
              A(M,J)=A(K,J)
              A(K,J)=T
              IF (T.EQ.0D0) GO TO 30
              DO 25 I=KP1,N
                  A(I,J)=A(I,J)+A(I,K)*T
   25         CONTINUE
   30     CONTINUE
   35 CONTINUE
      DO 50 K=1,N
          T=0D0
          IF (K.EQ.1) GO TO 45
          KM1=K-1
          DO 40 I=1,KM1
              T=T+A(I,K)*WORK(I)
   40     CONTINUE
   45     EK=1D0
          IF (T.LT.0D0) EK=-1D0
          IF (A(K,K).EQ.0D0) GO TO 90
          WORK(K)=-(EK+T)/A(K,K)
   50 CONTINUE
      DO 60 KB=1,NM1
          K=N-KB
          T=0D0
          KP1=K+1
          DO 55 I=KP1,N
              T=T+A(I,K)*WORK(K)
   55     CONTINUE
          WORK(K)=T
          M=IPVT(K)
          IF (M.EQ.K) GO TO 60
          T=WORK(M)
          WORK(M)=WORK(K)
          WORK(K)=T
   60 CONTINUE
      YNORM=0D0
      DO 65 I=1,N
          YNORM=YNORM+DABS(WORK(I))
   65 CONTINUE
*
      CALL SOLVE (NDIM,N,A,WORK,IPVT)
*
      ZNORM=0D0
      DO 70 I=1,N
          ZNORM=ZNORM+DABS(WORK(I))
   70 CONTINUE
      COND=ANORM*ZNORM/YNORM
      IF (COND.LT.1d0) COND=1D0
      RETURN
   80 COND=1D0
      IF (A(1,1).NE.0D0) RETURN
   90 COND=1D52
*
      RETURN
      END

C**********************************************************************

      SUBROUTINE SOLVE (NDIM,N,A,B,IPVT)
C--------/---------/---------/---------/---------/---------/---------/--
C >>> EPS
C <<< RAT
C=================
C--------/---------/---------/---------/---------/---------/---------/--
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 A(NDIM,N),B(N)
      INTEGER IPVT(N)
*
      IF (N.EQ.1) GO TO 50
      NM1=N-1
      DO 20 K=1,NM1
          KP1=K+1
          M=IPVT(K)
          T=B(M)
          B(M)=B(K)
          B(K)=T
          DO 10 I=KP1,N
              B(I)=B(I)+A(I,K)*T
   10     CONTINUE
   20 CONTINUE
      DO 40 KB=1,NM1
          KM1=N-KB
          K=KM1+1
          B(K)=B(K)/A(K,K)
          T=-B(K)
          DO 30 I=1,KM1
              B(I)=B(I)+A(I,K)*T
   30     CONTINUE
   40 CONTINUE
   50 B(1)=B(1)/A(1,1)
*
      RETURN
      END

C*****************************************************************

      SUBROUTINE SAREA(D,RAT)
C--------/---------/---------/---------/---------/---------/---------/--
C >>> EPS
C <<< RAT
C=================
C--------/---------/---------/---------/---------/---------/---------/--
      IMPLICIT REAL*8 (A-H,O-Z)
      IF (D.GE.1) GO TO 10
      E=DSQRT(1D0-D*D)
      R=0.5D0*(D**(2D0/3D0) + D**(-1D0/3D0)*DASIN(E)/E)
      R=DSQRT(R)
      RAT=1D0/R
      RETURN
   10 E=DSQRT(1D0-1D0/(D*D))
      R=0.25D0*(2D0*D**(2D0/3D0) + D**(-4D0/3D0)*DLOG((1D0+E)/(1D0-E))
     &   /E)
      R=DSQRT(R)
      RAT=1D0/R
*
      return
      END

c****************************************************************

      SUBROUTINE SURFCH(N,E,RAT)
C--------/---------/---------/---------/---------/---------/---------/--
C >>> N,E,RAT
C <<< RAT
C=================
C--------/---------/---------/---------/---------/---------/---------/--
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 X(60),W(60)
*
      DN=DFLOAT(N)
      EN=E*DN
      NG=60
*
* GIF division points and weights
*
      CALL GAUSS(NG,0,0,X,W)
*
      S=0D0
      V=0D0
      DO 10 I=1,NG
         XI=X(I)
         DX=DACOS(XI)
         DXN=DN*DX
         DS=DSIN(DX)
         DSN=DSIN(DXN)
         DCN=DCOS(DXN)
         A=1D0+E*DCN
         A2=A*A
         ENS=EN*DSN
         S=S+W(I)*A*DSQRT(A2+ENS*ENS)
         V=V+W(I)*(DS*A+XI*ENS)*DS*A2
   10 CONTINUE
      RS=DSQRT(S*0.5D0)
      RV=(V*3D0/4D0)**(1D0/3D0)
      RAT=RV/RS
*
      RETURN
      END

C********************************************************************

      SUBROUTINE SAREAC(EPS,RAT)
C--------/---------/---------/---------/---------/---------/---------/--
C >>> EPS
C <<< RAT
C=================
C--------/---------/---------/---------/---------/---------/---------/--
      IMPLICIT REAL*8 (A-H,O-Z)
*
      RAT=(1.5D0/EPS)**(1D0/3D0)
      RAT=RAT/DSQRT( (EPS+2D0)/(2D0*EPS) )
*
      RETURN
      END

C**********************************************************************

      SUBROUTINE DROP(RAT)
C--------/---------/---------/---------/---------/---------/---------/--
C >>> EPS
C <<< RAT
C=================
C--------/---------/---------/---------/---------/---------/---------/--
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER NOUT
* number of the output unit
      PARAMETER (NOUT=35)
      PARAMETER (NC=10, NG=60)

      REAL*8 X(NG),W(NG),C(0:NC)
      COMMON /CDROP/ C,R0V
      C(0)=-0.0481 D0
      C(1)= 0.0359 D0
      C(2)=-0.1263 D0
      C(3)= 0.0244 D0
      C(4)= 0.0091 D0
      C(5)=-0.0099 D0
      C(6)= 0.0015 D0
      C(7)= 0.0025 D0
      C(8)=-0.0016 D0
      C(9)=-0.0002 D0
      C(10)= 0.0010 D0
*
* GIF division points and weights
*
      CALL GAUSS(NG,0,0,X,W)
*
      S=0D0
      V=0D0
      DO I=1,NG
         XI=DACOS(X(I))
         WI=W(I)
         RI=1D0+C(0)
         DRI=0D0
         DO N=1,NC
            XIN=XI*N
            RI=RI+C(N)*DCOS(XIN)
            DRI=DRI-C(N)*N*DSIN(XIN)
         ENDDO
         SI=DSIN(XI)
         CI=X(I)
         RISI=RI*SI
         S=S+WI*RI*DSQRT(RI*RI+DRI*DRI)
         V=V+WI*RI*RISI*(RISI-DRI*CI)
      ENDDO
      RS=DSQRT(S*0.5D0)
      RV=(V*3D0*0.25D0)**(1D0/3D0)
      IF (DABS(RAT-1D0).GT.1D-8) RAT=RV/RS
      R0V=1D0/RV
      WRITE(NOUT,1000) R0V
      DO N=0,NC
         WRITE(NOUT,1001) N,C(N)
      ENDDO
 1000 FORMAT ('r_0/r_ev=',F7.4)
 1001 FORMAT ('c_',I2,'=',F7.4)

      RETURN
      END


      SUBROUTINE GAUSS(N,IND1,IND2,Z,W)
C--------/---------/---------/---------/---------/---------/---------/--
C >>> N,IND1,IND2
C <<< Z,W
C=================
C    CALCULATION OF POINTS AND WEIGHTS OF GAUSSIAN QUADRATURE
C    FORMULA. IF IND1 = 0 - ON INTERVAL (-1,1), IF IND1 = 1 - ON
C    INTERVAL  (0,1). IF  IND2 = 1 RESULTS ARE PRINTED.
C
C    N - NUMBER OF GIF DIVISION POINTS (mostly N=NGAUSS in main program)
C    Z - DIVISION POINTS
C    W - WEIGHTS
C--------/---------/---------/---------/---------/---------/---------/--
      IMPLICIT REAL*8 (A-H,P-Z)
      REAL*8 Z(N),W(N)
      DATA A,B,C /1D0,2D0,3D0/
      IND=MOD(N,2)
      K=N/2+IND
      F=DFLOAT(N)
      DO 100 I=1,K
          M=N+1-I
          IF(I.EQ.1) X=A-B/((F+A)*F)
          IF(I.EQ.2) X=(Z(N)-A)*4D0+Z(N)
          IF(I.EQ.3) X=(Z(N-1)-Z(N))*1.6D0+Z(N-1)
          IF(I.GT.3) X=(Z(M+1)-Z(M+2))*C+Z(M+3)
          IF(I.EQ.K.AND.IND.EQ.1) X=0D0
          NITER=0
          CHECK=1D-16
   10     PB=1D0
          NITER=NITER+1
          IF (NITER.LE.100) GO TO 15
          CHECK=CHECK*10D0
   15     PC=X
          DJ=A
          DO 20 J=2,N
              DJ=DJ+A
              PA=PB
              PB=PC
   20         PC=X*PB+(X*PB-PA)*(DJ-A)/DJ
          PA=A/((PB-X*PC)*F)
          PB=PA*PC*(A-X*X)
          X=X-PB
          IF(DABS(PB).GT.CHECK*DABS(X)) GO TO 10
          Z(M)=X
          W(M)=PA*PA*(A-X*X)
          IF(IND1.EQ.0) W(M)=B*W(M)
          IF(I.EQ.K.AND.IND.EQ.1) GO TO 100
          Z(I)=-Z(M)
          W(I)=W(M)
  100 CONTINUE
      IF(IND2.NE.1) GO TO 110
      PRINT 1100,N
 1100 FORMAT(' ***  POINTS AND WEIGHTS OF GAUSSIAN QUADRATURE FORMULA',
     * ' OF ',I4,'-TH ORDER')
      DO 105 I=1,K
          ZZ=-Z(I)
  105     PRINT 1200,I,ZZ,I,W(I)
 1200 FORMAT(' ',4X,'X(',I4,') = ',F17.14,5X,'W(',I4,') = ',F17.14)
      GO TO 115
  110 CONTINUE
C     PRINT 1300,N
C 1300 FORMAT(' GAUSSIAN QUADRATURE FORMULA OF ',I4,'-TH ORDER IS USED')
  115 CONTINUE
      IF(IND1.EQ.0) GO TO 140
      DO 120 I=1,N
  120     Z(I)=(A+Z(I))/B
  140 CONTINUE

      RETURN
      END

      SUBROUTINE gauleg(x1,x2,x,w,n)
C--------/---------/---------/---------/---------/---------/---------/--
C  Given the lower and upper limits of integration x1 and x2, and given n
C  this routine returns arrays x(1:n) and w(1:n) of length n, containing
C  the abscissas and weights of the Gaussian-Legendre n-point quadrature
C  formula.
C--------/---------/---------/---------/---------/---------/---------/--
      INTEGER n
      REAL*8 x1,x2,x(n),w(n)
      DOUBLE PRECISION EPS
      PARAMETER (EPS=3.d-14)
      INTEGER i,j,m
      DOUBLE PRECISION p1,p2,p3,pp,xl,xm,z,z1

      m=(n+1)/2          !The roots are symmetric in the interval, so we only
      xm=0.5d0*(x2+x1)   !have to find half of them
      xl=0.5d0*(x2-x1)

* Loop over the desired roots:

      do 12 i=1,m
        z=cos(3.141592654d0*(i-.25d0)/(n+.5d0))
              !Starting with the above approximation to the ith root, we enter
              !the main loop of refinement by Newton's method.
 1      continue
          p1=1.d0
          p2=0.d0

          do 11 j=1,n         !Loop up the recurrence relation to get Legendre
            p3=p2             !polynomial evaluated at z.
            p2=p1
            p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
 11       continue

* p1 is now the desired  Legendre polynomial. We next compute pp, its derivative,
* by a standard relation involving also p2, the polynomial of one lower order:

          pp=n*(z*p1-p2)/(z*z-1.d0)
          z1=z
          z=z1-p1/pp                   !Newton's method

        if (abs(z-z1).gt.EPS) goto 1

* Scale the root to the desired interval, and put in its symmetric counterpart:
        x(i)=xm-xl*z
        x(n+1-i)=xm+xl*z

* Compute the weight and its symmetric counterpart:
        w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
        w(n+1-i)=w(i)

 12   continue

      return
      END
      SUBROUTINE AMPMAT(NMAX,DLAM,TP,TP1,PP,PP1,KEX,CT)
C--------/---------/---------/---------/---------/---------/---------/--
C >>> NMAX,DLAM,TP,TP1,PP,PP1,CT
C <<< VV,VH,HV,HH,CEXT  (written directly to output files)
C=================
C* transfers ZEPS0,REV here from the main via COMMON
C
C    GIVEN T MATRIX IN COMMON BLOCK IT CALCULATES
C     THE AMPLITUDE MATRIX IN PARTICLE FRAME
C
C    This routine closely follows exposition by
C       M. I. Mishchenko, Calculation of the amplitude matrix
C       for a nonspherical particle in a fixed orientation,
C       Appl. Opt. vol. 39, 1026-1031 (2000).
C >>>
C   NMAX - angular momentum cutoff
C   DLAM=LAMBDA/SQRT(ZEPS0)  - wavelength of incident light in the ambient
C                     (vacuum wavelength divided by SQRT(ZEPS0))
C
C   LAMBDA - vacuum wavelength. Determined as DLAM*SQRT(ZEPS0) and
C            only used here for the write out purposes
C   TP,TP1,PP,PP1 ... respective incident and scattering beam angles
C                  (in degrees) determined w.r.t laboratory frame:
C   TP (THET0 IN MAIN) - zenith angle of the incident beam in degrees
C   TP1 (THET IN MAIN) - zenith angle of the scattered beam in degrees
C   PP (PHI0 IN MAIN) - azimuth angle of the incident beam in degrees
C   PP1 (PHI IN MAIN) - azimuth angle of the scattered beam in degrees
C   KEX ... the size (equiv. volume sphere) parameter in the ambient
C <<<
C   VV,VH,HV,HH ... amplitude scattering matrix elements S11,S12,S21,S22
C   CEXT ... extinction cross section for a fixed particle orientation
C--------/---------/---------/---------/---------/---------/---------/--
C      IMPLICIT NONE
      INTEGER NOUT,NPN1,NPN4,NPN6
        INTEGER SU,SUP,ST,ST1,ST2
        logical ynmishch

* number of the output unit
      PARAMETER (NOUT=35)
* Either Mischenko:
        Parameter (NPN1=100, NPN4=NPN1, NPN6=NPN4+1)
* or Bohren:
        PARAMETER (SU=50,SUP=0,ST=SU+SUP,ST1=ST+1,ST2=ST+ST)
* If routine is to be used with Mischenko convention, ynnmishch=.true.,
*  otherwise ynnmishch=.false.
        PARAMETER (ynmishch=.false.)

        INTEGER NMAX,NMIN,M,M1,N,NN
      REAL*8 DLAM,DK,LAMBDA,KEX,FAC,CEXT
        REAL*8 TP,TP1,PP,PP1,FC,FS,REV,PIN,PIN2,PI,THETP,PHIP,THETP1,
     & PHIP1,EPS,DNN,RN,DCTH0,DCTH,PH,DV1NN,DV2NN,DV1N,DV2N,D11,D12,
     & D21,D22
      REAL*8 DV1(NPN6),DV2(NPN6),DV01(NPN6),DV02(NPN6)
* Either Mischenko:
      REAL*4
     &     TR11(NPN6,NPN4,NPN4),TR12(NPN6,NPN4,NPN4),
     &     TR21(NPN6,NPN4,NPN4),TR22(NPN6,NPN4,NPN4),
     &     TI11(NPN6,NPN4,NPN4),TI12(NPN6,NPN4,NPN4),
     &     TI21(NPN6,NPN4,NPN4),TI22(NPN6,NPN4,NPN4)
* or Bohren:
        COMPLEX*16 CT(ST2,ST2,ST1)
*
      COMPLEX*16 CN,CN1,CN2,CI,VV,VH,HV,HH,ZEPS0,CT11,CT12,CT21,CT22,
     & CAL(NPN4,NPN4)

*_____
      COMMON /TMAT/ TR11,TR12,TR21,TR22,TI11,TI12,TI21,TI22
*_
      COMMON /DIELF/ zeps0
      COMMON /REVF/ rev
*
* transfers ZEPS0,REV here from the main
*____
*
C Checking the initial set of angles TP,TP1,PP,PP1
C for allowability

      IF (TP.LT.0D0.OR.TP.GT.180D0.OR.
     &    TP1.LT.0D0.OR.TP1.GT.180D0.OR.
     &    PP.LT.0D0.OR.PP.GT.360D0.OR.
     &    PP1.LT.0D0.OR.PP1.GT.360D0) THEN
          WRITE(NOUT,2000)
          STOP
      ELSE
          CONTINUE
      ENDIF
 2000 FORMAT ('AN ANGULAR PARAMETER IS OUTSIDE ITS',
     &        ' ALLOWABLE RANGE')

* SPECIFYING NUMERICAL CONSTANTS:

      PIN=DACOS(-1D0)         !=PI
      PIN2=PIN*0.5D0          !=PI/2
      PI=PIN/180D0            !=PI/180
* conversion from degrees to radians:
      THETP=TP*PI
      PHIP=PP*PI
      THETP1=TP1*PI
      PHIP1=PP1*PI
* initialization of the vacuum wavelength LAMBDA

        LAMBDA=DLAM*SQRT(ZEPS0)         !vacuum wavelength

      EPS=1D-7
      IF (THETP.LT.PIN2) THETP=THETP+EPS
      IF (THETP.GT.PIN2) THETP=THETP-EPS
      IF (THETP1.LT.PIN2) THETP1=THETP1+EPS
      IF (THETP1.GT.PIN2) THETP1=THETP1-EPS
      IF (PHIP.LT.PIN) PHIP=PHIP+EPS
      IF (PHIP.GT.PIN) PHIP=PHIP-EPS
      IF (PHIP1.LT.PIN) PHIP1=PHIP1+EPS
      IF (PHIP1.GT.PIN) PHIP1=PHIP1-EPS
C=========================================================
C      THE AMPLITUDE MATRIX IN PARTICLE FRAME
C
      CI=(0D0,1D0)

C >>> ALPHA numerical prefactors without phi-angles
C     (following Eq. (28))

      DO 5 NN=1,NMAX
         DO 5 N=1,NMAX
            CN=CI**(NN-N-1)
            DNN=DFLOAT((2*N+1)*(2*NN+1))
            DNN=DNN/DFLOAT( N*NN*(N+1)*(NN+1) )
            RN=DSQRT(DNN)
            CAL(N,NN)=CN*RN
    5 CONTINUE

      DCTH0=COS(TP)             !\cos\vartheta_{inc}^P
      DCTH=COS(TP1)             !\cos\vartheta_{sca}^P
      PH=PHIP1-PHIP         !(\varphi_{sca}^P-\varphi_{inc}^P)

* amplitude scattering matrix elements S11,S12,S21,S22 initialization

      VV=(0D0,0D0)
      VH=(0D0,0D0)
      HV=(0D0,0D0)
      HH=(0D0,0D0)
C______________________________________________________________
C Main summation loop:

      DO 500 M=0,NMAX
         M1=M+1
         NMIN=MAX(M,1)              !Bohren MTOPE
*
* * Specify pi- and tau- scattering functions:

         CALL VIGAMPL (DCTH, NMAX, M, DV1, DV2)
         CALL VIGAMPL (DCTH0, NMAX, M, DV01, DV02)
*
         FC=2D0*DCOS(M*PH)    !takes into account +/- m contribution
         FS=2D0*DSIN(M*PH)
*
         DO 400 NN=NMIN,NMAX

            DV1NN=M*DV01(NN)        !\pi-functions
            DV2NN=DV02(NN)          !\tau-functions

            DO 400 N=NMIN,NMAX
               DV1N=M*DV1(N)        !\pi-functions
               DV2N=DV2(N)          !\tau-functions

         if (ynmishch) then
* In Mischenko's  notation:
               CT11=DCMPLX(TR11(M1,N,NN),TI11(M1,N,NN))
               CT22=DCMPLX(TR22(M1,N,NN),TI22(M1,N,NN))
         else
* In Barber's notation:
               CT11=CT(N,NN,M1)
               CT22=CT(N+NMAX,NN+NMAX,M1)
         end if

               IF (M.EQ.0) THEN   !T^{21}=T^{12}=0 in particle frame

                  CN=CAL(N,NN)*DV2N*DV2NN

                  VV=VV+CN*CT22
                  HH=HH+CN*CT11

                 ELSE   !T^{21}\neq T^{12}\neq 0

                  if (ynmishch) then
* In Mischenko's  notation:
                  CT12=DCMPLX(TR12(M1,N,NN),TI12(M1,N,NN))
                  CT21=DCMPLX(TR21(M1,N,NN),TI21(M1,N,NN))
                  else
* In Barber's notation:
                  CT12=CT(N,NN+NMAX,M1)
                  CT21=CT(N+NMAX,NN,M1)
                  end if

* complete \alpha-factors (Eq. (28)) taking
* into account w.r.t. summation over +/- m in particle frame:
*
*     T^{11}_{-mnn'} = T^{11}_{mnn'}; T^{22}_{-mnn'} = T^{22}_{mnn'}
*  T^{12}_{-mnn'} = - T^{12}_{mnn'}; T^{21}_{-mnn'} = - T^{21}_{mnn'}
*
                  CN1=CAL(N,NN)*FC
                  CN2=CAL(N,NN)*FS

                  D11=DV1N*DV1NN    !\pi-\pi
                  D12=DV1N*DV2NN    !\pi-\tau
                  D21=DV2N*DV1NN    !\tau-\pi
                  D22=DV2N*DV2NN    !\tau-\tau

                  VV=VV+(CT11*D11+CT21*D21
     &                  +CT12*D12+CT22*D22)*CN1

                  VH=VH+(CT11*D12+CT21*D22
     &                  +CT12*D11+CT22*D21)*CN2

                  HV=HV-(CT11*D21+CT21*D11
     &                  +CT12*D22+CT22*D12)*CN2

                  HH=HH+(CT11*D22+CT21*D12
     &                  +CT12*D21+CT22*D11)*CN1
               ENDIF

  400    CONTINUE      !(over n,n')
  500 CONTINUE         !end of main summation loop (over m)

C Final multiplication of S11,S12,S21,S22 by (1/k)

      DK=2D0*PIN/DLAM   !wavevector in surrounding medium
      VV=VV/DK
      VH=VH/DK
      HV=HV/DK
      HH=HH/DK

C   amplitude scattering matrix elements S11,S12,S21,S22 determined


      PRINT 1101, VV
      PRINT 1102, VH
      PRINT 1103, HV
      PRINT 1104, HH

* For particles with plane of symmetry:

      cext=2.d0*PIN*dimag(VV+HH)/dlam       !Eq. (5.97)

      write(6,*)'C_{ext}= \fr{2\pi}{lambda} \mb{Im } (S_{11}+S_{22})=',
     & cext               !=2.d0*PIN*dimag(VV+HH)/lambda
*
      if (ynmishch) then
      FAC=lambda**2/(2.d0*PIN**2*REV**2)     !=2/xs**2
      else
      FAC=2.d0/KEX**2
        end if
*
        write(nout+3,1105) lambda, fac*cext, cext
      write(nout+12,1105) lambda, VV,VH
      write(nout+12,1105) lambda, HV,HH
      write(nout+13,1106) lambda, dble(VV*dconjg(vv) + VH*dconjg(vh)),
     &                               dble((vv+vh)*dconjg(vv+vh))
      write(nout+14,1106) lambda, dble(hV*dconjg(hv) + hH*dconjg(hh)),
     &                               dble((hv+hh)*dconjg(hv+hh))


 1101 FORMAT ('S11=',D11.5,' + i*',D11.5)
 1102 FORMAT ('S12=',D11.5,' + i*',D11.5)
 1103 FORMAT ('S21=',D11.5,' + i*',D11.5)
 1104 FORMAT ('S22=',D11.5,' + i*',D11.5)
 1105 FORMAT (F8.2,5X,D11.5,2X,D11.5,5X,D11.5,2X,D11.5)
 1106 FORMAT (F8.2,5X,D11.5,5X,D11.5)


      RETURN
      END
        SUBROUTINE beschb(x,gam1,gam2,gampl,gammi)
c--------/---------/---------/---------/---------/---------/---------/--
c Modified to a double precision subroutine Aug. 17, 1998
c Evaluates
C          gam1=\fr{1}{2\nu} [\Gamma^{-1}(1-\nu)-\Gamma^{-1}(1+\nu)]
C          gam2=  \fr{1}{2}  [\Gamma^{-1}(1-\nu)+\Gamma^{-1}(1+\nu)]
c by Chebyshev expansion for |x|\leq 2
c by Chebyshev expansion for |x|\leq 1/2. Also returns \Gamma^{-1}(1-\nu)
c and \Gamma^{-1}(1+\nu)
c--------/---------/---------/---------/---------/---------/---------/--
      INTEGER NUSE1,NUSE2
      DOUBLE PRECISION gam1,gam2,gammi,gampl,x
      PARAMETER (NUSE1=7,NUSE2=8)
CU    USES chebev
      REAL*8 xx,c1(7),c2(8),chebev
      SAVE c1,c2
      DATA c1/-1.142022680371172d0,6.516511267076d-3,3.08709017308d-4,
     *-3.470626964d-6,6.943764d-9,3.6780d-11,-1.36d-13/
      DATA c2/1.843740587300906d0,-.076852840844786d0,1.271927136655d-3,
     *-4.971736704d-6,-3.3126120d-8,2.42310d-10,-1.70d-13,-1.d-15/
      xx=8.d0*x*x-1.d0
      gam1=chebev(-1.d0,1.d0,c1,NUSE1,xx)
      gam2=chebev(-1.d0,1.d0,c2,NUSE2,xx)
      gampl=gam2-x*gam1
      gammi=gam2+x*gam1
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software

      SUBROUTINE bessik(x,xnu,ri,rkmu,rip,rkmup)
C--------/---------/---------/---------/---------/---------/---------/--
c >>> x,xnu
c <<< ri,rk,rip,rkp
c                ==========================================
c           Calculates the modified Bessel functions of arbitrary order
c Returns ri=I_\nu, rip=I_\nu', rk=K_\mu, rkmup=K_\mu'
c                ==========================================
c Recurrence relations for the Bessel functions are stable
c     1) downwards for J_\nu and J_\nu':
c                 J_{\nu-1} =(\nu/x) J_\nu + J_\nu'
c                 J_{\nu-1}'=((\nu-1)/x) J_{\nu-1} - J_\nu'
c     2) upwards for Y_\nu and Y_\nu':
c                 Y_{\nu+1} =(2\nu/x) Y_\nu - Y_{\nu-1}
c                 Y_{\nu}'=(\nu/x) Y_{\nu} - Y_{\nu+1}
C-----------------------------------------------------------------------
      INTEGER MAXIT
      REAL*8 ri,rip,rk,rkp,x,xnu,XMIN
      DOUBLE PRECISION EPS,FPMIN,PI
      PARAMETER (EPS=1.d-16,FPMIN=1.d-30,MAXIT=10000,XMIN=2.d0,
     *PI=3.141592653589793d0)
CU    USES beschb
      INTEGER i,l,nl
      DOUBLE PRECISION a,a1,b,c,d,del,del1,delh,dels,e,f,fact,fact2,ff,
     *gam1,gam2,gammi,gampl,h,p,pimu,q,q1,q2,qnew,ril,ril1,rimu,rip1,
     *ripl,ritemp,rk1,rkmu,rkmup,rktemp,s,sum,sum1,x2,xi,xi2,xmu,xmu2
      if(x.le.0..or.xnu.lt.0.) pause 'bad arguments in bessik'
      nl=int(xnu+.5d0)
      xmu=xnu-nl
      xmu2=xmu*xmu
      xi=1.d0/x
      xi2=2.d0*xi
      h=xnu*xi
      if(h.lt.FPMIN)h=FPMIN
      b=xi2*xnu
      d=0.d0
      c=h
      do 11 i=1,MAXIT
        b=b+xi2
        d=1.d0/(b+d)
        c=b+1.d0/c
        del=c*d
        h=del*h
        if(abs(del-1.d0).lt.EPS)goto 1
11    continue
      pause 'x too large in bessik; try asymptotic expansion'
1     continue
      ril=FPMIN
      ripl=h*ril
      ril1=ril
      rip1=ripl
      fact=xnu*xi
      do 12 l=nl,1,-1
        ritemp=fact*ril+ripl
        fact=fact-xi
        ripl=fact*ritemp+ril
        ril=ritemp
12    continue
      f=ripl/ril
      if(x.lt.XMIN) then
        x2=.5d0*x
        pimu=PI*xmu
        if(abs(pimu).lt.EPS)then
          fact=1.d0
        else
          fact=pimu/sin(pimu)
        endif
        d=-log(x2)
        e=xmu*d
        if(abs(e).lt.EPS)then
          fact2=1.d0
        else
          fact2=sinh(e)/e
        endif
        call beschb(xmu,gam1,gam2,gampl,gammi)
        ff=fact*(gam1*cosh(e)+gam2*fact2*d)
        sum=ff
        e=exp(e)
        p=0.5d0*e/gampl
        q=0.5d0/(e*gammi)
        c=1.d0
        d=x2*x2
        sum1=p
        do 13 i=1,MAXIT
          ff=(i*ff+p+q)/(i*i-xmu2)
          c=c*d/i
          p=p/(i-xmu)
          q=q/(i+xmu)
          del=c*ff
          sum=sum+del
          del1=c*(p-i*ff)
          sum1=sum1+del1
          if(abs(del).lt.abs(sum)*EPS)goto 2
13      continue
        pause 'bessk series failed to converge'
2       continue
        rkmu=sum
        rk1=sum1*xi2
      else
        b=2.d0*(1.d0+x)
        d=1.d0/b
        delh=d
        h=delh
        q1=0.d0
        q2=1.d0
        a1=.25d0-xmu2
        c=a1
        q=c
        a=-a1
        s=1.d0+q*delh
        do 14 i=2,MAXIT
          a=a-2*(i-1)
          c=-a*c/i
          qnew=(q1-b*q2)/a
          q1=q2
          q2=qnew
          q=q+c*qnew
          b=b+2.d0
          d=1.d0/(b+a*d)
          delh=(b*d-1.d0)*delh
          h=h+delh
          dels=q*delh
          s=s+dels
          if(abs(dels/s).lt.EPS)goto 3
14      continue
        pause 'bessik: failure to converge in cf2'
3       continue
********************************************************
* determination of I_\mu, I_\mu', K_\mu, and K_\mu'
        h=a1*h
        rkmu=sqrt(PI/(2.d0*x))*exp(-x)/s
        rk1=rkmu*(xmu+x+.5d0-h)*xi
* rk1=K_{\mu+1} here
      endif
      rkmup=xmu*xi*rkmu-rk1
* rkmup=K_{\mu}'
* scaling to obtain the original  J_\nu, J_\nu'
      rimu=xi/(f*rkmu-rkmup)
      ri=(rimu*ril1)/ril
      rip=(rimu*rip1)/ril
c* the stable upwards recurrence for K_\mu
c      do 15 i=1,nl
c        rktemp=(xmu+i)*xi2*rk1+rkmu
c        rkmu=rk1
c        rk1=rktemp
c15    continue
c      rk=rkmu
c      rkp=xnu*xi*rkmu-rk1
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software v1]3"iA#.

      SUBROUTINE bessjy(x,xnu,rj,rymu,rjp,ry1)
c Variables declared but never referenced:
c        RY                RYP             RYTEMP
c--------/---------/---------/---------/---------/---------/---------/--
c >>> x,xnu
c <<< j,ry,rjp,ryp
c                ==========================================
c       Calculates the cylindrical Bessel functions of arbitrary order.
c Returns rj=j_\nu, rjp=J_\nu', rymu=Y_\mu , ry1=Y_{\mu+1}
c       Applicable in the interval \approx (0, 10.000)
C       The spherical Bessel functions can be found using
c               j_l=\sqrt{\fr{\pi}{2z}}\, J_{n+1/2}, etc.
c
c  TESTED UP TO XNU=19 and X=50 - WORKS WITH AN EXCELLENT PRECISION !!!
c                ==========================================
c Recurrence relations for the Bessel functions are stable
c     1) downwards for J_\nu and J_\nu':
c                 J_{\nu-1} =(\nu/x) J_\nu + J_\nu'
c                 J_{\nu-1}'=((\nu-1)/x) J_{\nu-1} - J_\nu'
c     2) upwards for Y_\nu and Y_\nu':
c                 Y_{\nu+1} =(2\nu/x) Y_\nu - Y_{\nu-1}
c                 Y_{\nu}'=(\nu/x) Y_{\nu} - Y_{\nu+1}
c                                                 [cf. (9.1.27) of \ct{AS}]
c
c Uses two recurrence relations,
c    - CF1 which determines the ratio J_\nu/J_\nu'
c    - CF2 which determines the ratio p+iq=(J_\nu'+iY_\nu')/(J_\nu+iY_\nu)
c and the expression for the Wronskian
c    - W= J_\nu Y_\nu' - J_\nu' Y_\nu = 2/(\pi x)
c
c The rate of convergence of CF1 is determined by the position of the
c turning point x_{tp}=\sqrt{\nu(\nu+1)}, beyond which the Bessel functions
c become oscillatory. If x\leq x_{tp}, convergence is very rapid.
c On the other hand, CF2 converges rapidly for x\geq x_{tp}, while the
c convergence fails for x\rar 0.
c
c  COMMENT: Using (9.1.2) of \ct{AS} one has
c         (\cos[(n+1/2)\pi]=0,       \sin[(n+1/2)\pi]=(-1)^n)
c
c               Y_{n+1/2} = (-1)^{n+1} J_{-(n+1/2)}
c  and hence
c                   y_{n} = (-1)^{n+1} j_{-n-1}
c
c  and one can continue with a stable recurrence downwards to generate
c  also the relevant Neumann (Weber) functions.
c
c                        E X E C U T I O N
c
c 1) Uses CF1 recurrence relation to calculate the ratio
c                     f_\nu=J_\nu/J_\nu'
c
c 2) For x<x_{min}=2 one uses downwards recurrence relations down to
c   |xmu|<1/2, while for x>x_{min}=2 one uses downwards recurrence
c   relations so that x\geq x_{tp} with the initial value of
c                       J_\nu=1.d-30
c
c 3) Exaluates CF2 for xmu. Then with W one has enough relations to solve
c    for all four functions and  determine them at xmu.
c    The sign of J_\xmu is chosen to be the same as the sign of the
c    initial J_\nu. For x<2, Temme's series are used instead of CF2.
c
c 4) Once all four functions are determined at xmu, J_\nu and J_\nu'
c    are determined by simply scaling the initial values by the ratio
c    of CF2+W result to that found after applying the recurrence 1).
c    The quantities Y_\nu and Y_\nu' are found using the stable upwards
c    recurrence 2).
c
c Modified to a double precision subroutine Aug. 17, 1998
c--------/---------/---------/---------/---------/---------/---------/--
      INTEGER MAXIT
      DOUBLE PRECISION rj,rjp,ry,ryp,x,xnu,XMIN
      DOUBLE PRECISION EPS,FPMIN,PI
      PARAMETER (EPS=1.d-16,FPMIN=1.d-30,MAXIT=10000,XMIN=2.d0,
     *PI=3.141592653589793d0)
CU    USES beschb
      INTEGER i,isign,l,nl
      DOUBLE PRECISION a,b,br,bi,c,cr,ci,d,del,del1,den,di,dlr,dli,dr,e,
     *f,fact,fact2,fact3,ff,gam,gam1,gam2,gammi,gampl,h,p,pimu,pimu2,q,
     *r,rjl,rjl1,rjmu,rjp1,rjpl,rjtemp,ry1,rymu,rymup,rytemp,sum,sum1,
     *temp,w,x2,xi,xi2,xmu,xmu2
      if(x.le.0..or.xnu.lt.0.) pause 'bad arguments in bessjy'
      if(x.lt.XMIN)then
        nl=int(xnu+.5d0)
      else
        nl=max(0,int(xnu-x+1.5d0))
      endif
      xmu=xnu-nl
      xmu2=xmu*xmu
      xi=1.d0/x
      xi2=2.d0*xi
      w=xi2/PI
      isign=1
      h=xnu*xi
      if(h.lt.FPMIN)h=FPMIN
      b=xi2*xnu
      d=0.d0
      c=h
      do 11 i=1,MAXIT
        b=b+xi2
        d=b-d
        if(abs(d).lt.FPMIN)d=FPMIN
        c=b-1.d0/c
        if(abs(c).lt.FPMIN)c=FPMIN
        d=1.d0/d
        del=c*d
        h=del*h
* h=J_\nu/J_\nu' is the result of CF1 relation
        if(d.lt.0.d0)isign=-isign
        if(abs(del-1.d0).lt.EPS)goto 1
11    continue
      pause 'x too large in bessjy; try asymptotic expansion'
1     continue
* Initialization of J_\nu and J_\nu' for the stable downwards recurrence 1)
      rjl=isign*FPMIN
      rjpl=h*rjl
      rjl1=rjl
      rjp1=rjpl
      fact=xnu*xi
      do 12 l=nl,1,-1
        rjtemp=fact*rjl+rjpl
        fact=fact-xi
        rjpl=fact*rjtemp-rjl
        rjl=rjtemp
12    continue
      if(rjl.eq.0.d0)rjl=EPS
      f=rjpl/rjl
* Now have unnormalized J_\mu and J_\mu' and f=J_\mu/J_\mu'
***********************************************************
* Dichotomy 0<x<2 or x\geq 2
      if(x.lt.XMIN) then
        x2=.5d0*x
        pimu=PI*xmu
        if(abs(pimu).lt.EPS)then
          fact=1.d0
        else
          fact=pimu/sin(pimu)
        endif
        d=-log(x2)
        e=xmu*d
        if(abs(e).lt.EPS)then
          fact2=1.d0
        else
          fact2=sinh(e)/e
        endif
        call beschb(xmu,gam1,gam2,gampl,gammi)
        ff=2.d0/PI*fact*(gam1*cosh(e)+gam2*fact2*d)
        e=exp(e)
        p=e/(gampl*PI)
        q=1.d0/(e*PI*gammi)
        pimu2=0.5d0*pimu
        if(abs(pimu2).lt.EPS)then
          fact3=1.d0
        else
          fact3=sin(pimu2)/pimu2
        endif
        r=PI*pimu2*fact3*fact3
        c=1.d0
        d=-x2*x2
        sum=ff+r*q
        sum1=p
        do 13 i=1,MAXIT
          ff=(i*ff+p+q)/(i*i-xmu2)
          c=c*d/i
          p=p/(i-xmu)
          q=q/(i+xmu)
          del=c*(ff+r*q)
          sum=sum+del
          del1=c*p-i*del
          sum1=sum1+del1
          if(abs(del).lt.(1.d0+abs(sum))*EPS)goto 2
13      continue
        pause 'bessy series failed to converge'
2       continue
        rymu=-sum
        ry1=-sum1*xi2
        rymup=xmu*xi*rymu-ry1
        rjmu=w/(rymup-f*rymu)
      else
        a=.25d0-xmu2
        p=-.5d0*xi
        q=1.d0
        br=2.d0*x
        bi=2.d0
        fact=a*xi/(p*p+q*q)
        cr=br+q*fact
        ci=bi+p*fact
        den=br*br+bi*bi
        dr=br/den
        di=-bi/den
        dlr=cr*dr-ci*di
        dli=cr*di+ci*dr
        temp=p*dlr-q*dli
        q=p*dli+q*dlr
        p=temp
        do 14 i=2,MAXIT
          a=a+2*(i-1)
          bi=bi+2.d0
          dr=a*dr+br
          di=a*di+bi
          if(abs(dr)+abs(di).lt.FPMIN)dr=FPMIN
          fact=a/(cr*cr+ci*ci)
          cr=br+cr*fact
          ci=bi-ci*fact
          if(abs(cr)+abs(ci).lt.FPMIN)cr=FPMIN
          den=dr*dr+di*di
          dr=dr/den
          di=-di/den
          dlr=cr*dr-ci*di
          dli=cr*di+ci*dr
          temp=p*dlr-q*dli
          q=p*dli+q*dlr
          p=temp
          if(abs(dlr-1.d0)+abs(dli).lt.EPS)goto 3
14      continue
        pause 'cf2 failed in bessjy'
3       continue
********************************************************
* determination of J_\mu, J_\mu', Y_\mu, and Y_\mu'
        gam=(p-f)/q
        rjmu=sqrt(w/((p-f)*gam+q))
        rjmu=sign(rjmu,rjl)
        rymu=rjmu*gam
        rymup=rymu*(p+q/gam)
        ry1=xmu*xi*rymu-rymup
* ry1=Y_{\mu+1} here
      endif
* scaling to obtain the original  J_\nu, J_\nu'
      fact=rjmu/rjl
      rj=rjl1*fact
      rjp=rjp1*fact
c* the stable upwards recurrence 2)  (xi2=2/x)
c      do 15 i=1,nl
c        rytemp=(xmu+i)*xi2*ry1-rymu
c        rymu=ry1
c        ry1=rytemp
c15    continue
c      ry=rymu
c      ryp=xnu*xi*rymu-ry1
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software v1]3"iA#.
      SUBROUTINE BIGA(CIOR,XX,NTRM,NOABS,YESANG,CBIGA)
ct      program biga
C--------/---------/---------/---------/---------/---------/---------/--
c    Calculate logarithmic derivatives of J-Bessel-function
c    Input :  CIOR, XX, NTRM, NOABS, YESANG
c    Output :  RBIGA or CBIGA
c    Routines called :  CONFRA
c              ERRMSG('ConFra--Iteration failed to converge',.TRUE.)
c  =====
c
c  COMPLEX*16   CBIGA( LMAXD )
c  call BIGA (CIOR, XX, LMAXD, .false., .false., RBIGA, CBIGA )
c
c According to:
c  [1] Lentz ???
c  [2] W. J. Wiscombe, Appl. Opt. 19, 1505-1509 (1980).
c                      Improved Mie scattering algorithms.
c     (http://www.opticsinfobase.org/ao/abstract.cfm?URI=ao-19-9-1505)
c  [R] W. J. Wiscombe, NCAR Technical Note
c
c =====
c  CIOR            Complex index of refraction. Should be
c                  now in scattering notation (Im part of n always nonnegative)
c  XX               size parameter
c  NTRM            No. of terms in Mie series (my LMAXV)
c  NOABS           TRUE, sphere non-absorbing (determined by MIMCUT=0.d0)
c  YESANG          TRUE if scattering amplitudes are to be calculated
c
c  RBIGA(N)        Bessel function ratio A_N=\psi_N'/\psi_N (Ref. 2, Eq. 2)
c                     (REAL version, for when imag refrac index = 0)
c  CBIGA(N)        Bessel function ratio A_N=\psi_N'/\psi_N (Ref. 2, Eq. 2)
c                     ( COMPLEX version )
c
c    INTERNAL VARIABLES :
c
c       CONFRA     Value of Lentz continued fraction for CBIGA(NTRM),
c                     used to initialize downward recurrence
c
c       DOWN       = True, use down-recurrence.  False, do not.
c
c       F1,F2,F3   Arithmetic statement functions used in determining
c                     whether to use up-  or down-recurrence
c                     ( Ref. 2, Eqs. 6-8 )
c
c       MRE        Real refractive index
c       MIM        Imaginary refractive index
c
c       REZINV     1 / ( MRE * XX ); temporary variable for recurrence
c       ZINV       1 / ( CIOR * XX ); temporary variable for recurrence
C--------/---------/---------/---------/---------/---------/---------/--

      IMPLICIT  NONE

c     .. Scalar Arguments ..

      INTEGER    NTRM,LMAXD
        PARAMETER (LMAXD=60)

      LOGICAL    NOABS, YESANG
      REAL*8     XX
      COMPLEX*16 CIOR
c     ..
c     .. Array Arguments ..

      REAL*8      RBIGA(LMAXD)
      COMPLEX*16   CBIGA(LMAXD)
c     ..
c     .. Local Scalars ..

      LOGICAL   DOWN
      INTEGER   N
      REAL*8      MIM, MRE, REZINV, RTMP
      COMPLEX*16   CTMP,ZINV
c     ..
c     .. External Functions ..

      COMPLEX*16   CONFRA
      EXTERNAL  CONFRA
c     ..
c     .. Intrinsic Functions ..
c
c      INTRINSIC ABS, AIMAG, COS, EXP, REAL, SIN
c     ..
c     .. Statement Functions ..
c
      REAL*8      F1,F2    !,F3  - not used here
c     ..
c     .. Statement Function definitions ..
c
c                                                   ** Eq. R47c
      F1(MRE) = -8.d0 + MRE**2*( 26.22d0 +
     &  MRE*( -0.4474d0 + MRE**3*( 0.00204d0 - 0.000175d0*MRE ) ) )

c                                                   ** Eq. R47b
      F2(MRE) = 3.9d0 + MRE*( -10.8d0 + 13.78d0*MRE )
c                                                   ** Eq. R47a
cc      F3(MRE) = -15.04d0 + MRE*( 8.42d0 + 16.35d0*MRE )
c     ..

ct      NOABS=.false.
ct      YESANG=.false.
ct      XX=7.31806289774349d0
ct      CIOR=(1.45d0,1.d0)

c                                  ** Decide whether BigA can be
c                                  ** calculated by up-recurrence
      MRE = DBLE(CIOR)
      MIM = ABS(IMAG(CIOR))

      IF( MRE.LT.1.0 .OR. MRE.GT.10.0 .OR. MIM.GT.10.0 ) THEN

         DOWN = .TRUE.

      ELSE IF (YESANG) THEN

         DOWN = .TRUE.
c                                                    ** Eq. R48
         IF( MIM*XX .LT. F2(MRE) ) DOWN = .FALSE.

      ELSE

         DOWN = .TRUE.
c                                                    ** Eq. R48
         IF( MIM*XX .LT. F1(MRE) ) DOWN = .FALSE.

      END IF


      ZINV   = 1.d0 /(CIOR*XX)
      REZINV = 1.d0 /(MRE*XX)


      IF(DOWN) THEN
c                          ** Compute initial high-order BigA using
c                          ** Lentz method ( Ref. 1, pp. 17-20 )

         CTMP = CONFRA(NTRM,ZINV)

c                                   *** Downward recurrence for BigA
         IF(NOABS) THEN
c                                        ** No-absorption case; Eq (R23)
            RBIGA(NTRM) = DBLE(CTMP)
            CBIGA(NTRM) = RBIGA(NTRM)

            DO 10 N = NTRM, 2, -1

               RBIGA(N - 1) = ( N*REZINV ) -
     &                        1.d0 / ( ( N*REZINV ) + RBIGA(N) )

               CBIGA(N-1) = RBIGA(N-1)

   10       CONTINUE

         ELSE
c                                         ** Absorptive case; Eq (R23)
            CBIGA(NTRM) = CTMP

            DO 20 N = NTRM, 2, -1
            CBIGA(N-1) = (N*ZINV) - 1.d0 / ( (N*ZINV) + CBIGA(N) )
   20       CONTINUE

         END IF


      ELSE
c                            *** Upward recurrence for BigA
         IF(NOABS) THEN
c                                  ** No-absorption case; Eq (R20,21)
            RTMP = SIN( MRE*XX )
            RBIGA(1) = - REZINV + RTMP /
     &                   ( RTMP*REZINV - COS(MRE*XX) )
            CBIGA(1) = RBIGA(1)

            DO 30 N = 2, NTRM
               RBIGA(N) = -( N*REZINV ) +
     &                      1.d0 / ( (N*REZINV) - RBIGA(N - 1) )
             CBIGA(N) = RBIGA(N)
   30       CONTINUE

         ELSE
c                                     ** Absorptive case; Eq (R20,22)

            CTMP = EXP((0.d0,2.d0)*CIOR*XX)
            CBIGA(1) = - ZINV + (1.d0-CTMP) /
     &              (ZINV * (1.d0-CTMP) + (0.d0,1.d0)*(1.d0+CTMP))

            DO 40 N = 2, NTRM
            CBIGA(N) = - (N*ZINV) + 1.d0/((N*ZINV) - CBIGA(N-1))
   40       CONTINUE

         END IF

      END IF

      RETURN
      END

      COMPLEX*16 FUNCTION CONFRA(N,ZINV)
C--------/---------/---------/---------/---------/---------/---------/--
c         Compute Bessel function ratio A-sub-N from its
c         continued fraction using Lentz method
c
c         ZINV = Reciprocal of argument of A
c
c    I N T E R N A L    V A R I A B L E S
c    ------------------------------------
c    CAK      Term in continued fraction expansion of A (Eq. R25)
c    CAPT     Factor used in Lentz iteration for A (Eq. R27)
c    CNUMER   Numerator   in capT  ( Eq. R28A )
c    CDENOM   Denominator in capT  ( Eq. R28B )
c    CDTD     Product of two successive denominators of capT factors
c                 ( Eq. R34C )
c    CNTN     Product of two successive numerators of capT factors
c                 ( Eq. R34B )
c    EPS1     Ill-conditioning criterion
c    EPS2     Convergence criterion
c    KK       Subscript k of cAk  ( Eq. R25B )
c    KOUNT    Iteration counter ( used to prevent infinite looping )
c    MAXIT    Max. allowed no. of iterations
c    MM       + 1  and - 1, alternately
C--------/---------/---------/---------/---------/---------/---------/--
      IMPLICIT  NONE

c     .. Scalar Arguments ..

      INTEGER   N
      COMPLEX*16   ZINV
c     ..
c     .. Local Scalars ..

      INTEGER   KK, KOUNT, MAXIT, MM
      REAL*8      EPS1, EPS2
      COMPLEX*16   CAK, CAPT, CDENOM, CDTD, CNTN, CNUMER
c     ..
c     .. External Subroutines ..

cc      EXTERNAL  ERRMSG
c     ..
c     .. Intrinsic Functions ..
c
c      INTRINSIC ABS, AIMAG, REAL
c     ..
      DATA      EPS1 / 1.d-2 / , EPS2 / 1.d-8 /
      DATA      MAXIT / 10000 /

c                                 ** Eq. R25a
      CONFRA = ( N + 1 ) * ZINV
      MM     = - 1
      KK     = 2*N + 3
c                                 ** Eq. R25b, k=2
      CAK    = ( MM*KK ) * ZINV
      CDENOM = CAK
      CNUMER = CDENOM + 1.d0 / CONFRA
      KOUNT  = 1

   10 CONTINUE
      KOUNT = KOUNT + 1

      IF (KOUNT.GT.MAXIT) then
cx     &    CALL ERRMSG('ConFra--Iteration failed to converge',.TRUE.)
          write(6,*)'ConFra--Iteration failed to converge'
          pause
      end if

      MM  = - MM
      KK  = KK + 2
c                                 ** Eq. R25b
      CAK = ( MM*KK ) * ZINV
c                                          ** Eq. R32
      IF( ABS( CNUMER / CAK ).LE.EPS1 .OR.
     &    ABS( CDENOM / CAK ).LE.EPS1 ) THEN

c                                  ** Ill-conditioned case -- stride
c                                  ** two terms instead of one
c                                       ** Eq. R34
         CNTN   = CAK * CNUMER + 1.d0
         CDTD   = CAK * CDENOM + 1.d0
c                                           ** Eq. R33
         CONFRA = ( CNTN / CDTD ) * CONFRA

         MM  = - MM
         KK  = KK + 2
c                                 ** Eq. R25b
         CAK = ( MM*KK ) * ZINV
c                                      ** Eq. R35
         CNUMER = CAK + CNUMER / CNTN
         CDENOM = CAK + CDENOM / CDTD
         KOUNT  = KOUNT + 1
         GO TO  10

      ELSE
c                           *** Well-conditioned case
c                                  ** Eq. R27
         CAPT   = CNUMER / CDENOM
c                                  ** Eq. R26
         CONFRA = CAPT * CONFRA
c                                  ** Check for convergence; Eq. R31

         IF (      ABS( DBLE (CAPT) - 1.d0 ).GE.EPS2
     &        .OR. ABS( IMAG(CAPT) )      .GE.EPS2 )  THEN

c                                        ** Eq. R30
            CNUMER = CAK + 1.d0 / CNUMER
            CDENOM = CAK + 1.d0 / CDENOM

            GO TO  10

         END IF

      END IF

      RETURN
      END

            DOUBLE PRECISION FUNCTION chebev(a,b,c,m,x)
      INTEGER m
      REAL*8 a,b,x,c(m)
      INTEGER j
      REAL*8 d,dd,sv,y,y2
      if ((x-a)*(x-b).gt.0.) pause 'x not in range in chebev'
      d=0.d0
      dd=0.d0
      y=(2.d0*x-a-b)/(b-a)
      y2=2.d0*y
      do 11 j=m,2,-1
        sv=d
        d=y2*d-dd+c(j)
        dd=sv
11    continue
      chebev=y*d-dd+0.5d0*c(1)
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software v1]3"iA#.
      subroutine GNRCBSH(ZX,LMX,zeta,dzeta)
C--------/---------/---------/---------/---------/---------/---------/--
C >>> CQEPS,XX,LMX
C <<< zeta,dzeta from l=0 up to l=LMX
C =====
C Calculates an array of Riccati-Hankel functions of the
C argument CQEPS*XX using Mackowski et al recurrences
C--------/---------/---------/---------/---------/---------/---------/--
      implicit none
      integer LMAXD

C  >>> ANGULAR MOMENTUM CUTOFF ON ARRAY DIMENSION
C      (The actual angular-momentum cutoff on summation is specified
C       by the value of variable LMX)

      PARAMETER (LMAXD=60)

      real*8 xt
      integer lmx,l
      complex*16 zx,z1,ci,cone,CQEPSt
      complex*16 dr1(0:LMAXD),dr3(0:LMAXD),zeta(0:lmx),dzeta(0:lmx)

      DATA ci/(0.d0,1.d0)/,cone/(1.d0,0.d0)/
C--------/---------/---------/---------/---------/---------/---------/--

      z1=dcmplx(dble(sqrt(zx*dconjg(zx))),0.d0)
      cqepst=zx/z1
        xt=dble(z1)

      if (imag(CQEPST).ne.0.d0) then
         call BIGA(CQEPST,xt,LMX,.false.,.false.,dr1(1))
      else if (imag(CQEPST).eq.0.d0) then
         call BIGA(CQEPST,xt,LMX,.true.,.false.,dr1(1))
      end if

        z1=-ci*exp(ci*zx)*sin(zx)
        DR3(0)=ci
        DR1(0)=cone/zx - cone/(cone/zx + dr1(1))
        zeta(0)=-ci*exp(ci*zx)

        do 15  l=1,lmx

         z1=z1*(-DR1(l-1)+dble(l)/zx)*(-DR3(l-1)+dble(l)/zx)
         DR3(l)=DR1(l)+ci/z1
           zeta(l)=zeta(l-1)*(-DR3(l-1)+dble(l)/zx)
           dzeta(l)=dr3(l)*zeta(l)

  15     continue

      return
      end

            subroutine GNRICBESSH(CQEPS,XX,LMX,zeta,dzeta)
C--------/---------/---------/---------/---------/---------/---------/--
C >>> CQEPS,XX,LMX
C <<< zeta,dzeta of argument CQEPS*XX from l=0 up to l=LMX
C =====
C Calculates an array of Riccati-Hankel functions of the
C argument CQEPS*XX using recurrences in Eqs. (63)-(64), (66)-(67) of
C
C [1] D. W. Mackowski, R. A. Altenkirch, and M. P. Menguc,
C "Internal absorption cross sections in a stratified sphere,"
C Appl. Opt. 29, 1551-1559 (1990)
C http://www.opticsinfobase.org/ao/abstract.cfm?URI=ao-29-10-1551
C
C z1     ... \psi_l\xi_l
C DR1(l) ... DR^{(1)}_l=\psi_l'/\psi_l
C DR3(l) ... DR^{(3)}_l= \xi_l'/\xi_l
C
C An array DR1(l) is supplied by biga.f, wherein it is
C calculated by the downward recurrence relation [e.g. (62) of [1])
C--------/---------/---------/---------/---------/---------/---------/--
      implicit none
      integer LMAXD

C  >>> ANGULAR MOMENTUM CUTOFF ON ARRAY DIMENSION
C      (The actual angular-momentum cutoff on summation is specified
C       by the value of variable LMX)

      PARAMETER (LMAXD=60)

      integer lmx,l
      real*8 xx
      complex*16 zx,z1,ci,cone,CQEPS
      complex*16 dr1(0:LMAXD),dr3(0:LMAXD),zeta(0:lmx),dzeta(0:lmx)

      DATA ci/(0.d0,1.d0)/,cone/(1.d0,0.d0)/
C--------/---------/---------/---------/---------/---------/---------/--

C Determine the Bessel function ratio A_1=\psi_1'/\psi_1

      if (imag(CQEPS).ne.0.d0) then
         call BIGA(CQEPS,xx,LMX,.false.,.false.,dr1(1))
      else if (imag(CQEPS).eq.0.d0) then
         call BIGA(CQEPS,xx,LMX,.true.,.false.,dr1(1))
      end if

        zx=xx*cqeps

* initial values of the recurrence (63),(64),(67) of [1]
* \psi_0=zx*j_0=sin(zx); \xi_0=zx*h_0=-ci*exp(ci*zx)

        z1=-ci*exp(ci*zx)*sin(zx)
        DR3(0)=ci
        DR1(0)=cone/zx - cone/(cone/zx + dr1(1))
        zeta(0)=-ci*exp(ci*zx)

*
        do 15  l=1,lmx   !Eqs. (63),(64),(67) of Mackowski et al

         z1=z1*(-DR1(l-1)+dble(l)/zx)*(-DR3(l-1)+dble(l)/zx)
         DR3(l)=DR1(l)+ci/z1
           zeta(l)=zeta(l-1)*(-DR3(l-1)+dble(l)/zx)
           dzeta(l)=dr3(l)*zeta(l)

  15     continue

      return
      end

            SUBROUTINE gnzbess(omega,lmax,jl,djl,nl,dnl)

*   Variables declared but never referenced:
*        CRY1              CRYMU
*   Variables used before set:
*       CRJ               CRJP
C--------/---------/---------/---------/---------/---------/---------/--
C   ROUTINE TO GENERATE ARRAYS OF THE SPHERICAL BESSEL FUNCTION
C       AND THEIR DERIVATIVES FOR A COMPLEX ARGUMENT AND REAL ORDER.
C           (RETURNS jl(i*|OM|), nl(i*|OM|), etc. IF OM IMAGINARY)
C                 ================================
C Limitation : bessjy and bessik require real argument on the input
C
C*****************************************************************
C                      >>> FOR OMEGA REAL  <<<
C  Calculated values are compared with Table 10.5 p. 465 of \ct{AS}
C      (The number (-n) in parenthesis there means 10^{-n})
C
C According to (9.1.60) of \ct{AS}:
C
C           |J_\nu(x)|\leq 1           for  \nu\geq 0
C           |J_\nu(x)|\leq 1/\sqrt{2}  for  \nu\geq 1
C One has
C            j_l=\sqrt{\fr{\pi}{2z}}\, J_{l+1/2}       (10.1.1)
C which when combined with (9.1.60) gives bounds on j_l.
C Calls routine bessjy from Numerical Recipies which returns
C the cylindrical Bessel functions and their derivatives of a given
C order. The routine then determines the rest using the stable downwards
C recurrence relations  for j_\nu and j_\nu' [cf. (10.1.21-22) of \ct{AS}]:
c                 j_{\nu-1} =((\nu+1)/x) j_\nu + j_\nu'
c                 j_{\nu-1}'=((\nu-1)/x) j_{\nu-1} - j_\nu/(2.d0*x).
C The spherical y_n Bessel functions are then recovered using the formulae
C of \ct{AS}
C                 y_{n} = (-1)^{n+1} j_{-n-1}           (9.1.2)
C*****************************************************************
C                >>> FOR PURELY IMAGINARY OMEGA   <<<
C
C Calls routine bessik from Numerical Recipies which returns
C the modified cylindrical Bessel functions and their derivatives of a given
C order. The routine then determines the rest using the stable downwards
C recurrence relations for j_\nu and j_\nu' [cf. (10.1.21-22) of \ct{AS}]:
c
C  Calculated values are compared with Table 10.10 p. 473 of \ct{AS}
C      (The number (-n) in parenthesis there means 10^{-n})
C--------------------------
C According to (9.6.27) of \ct{AS}, I_0'(z)=I_1(z)
C Recurrence for I_\nu [cf. (9.6.26) of \ct{AS}] is stored in bessel.bal:
c                 I_{\nu-1} =(\nu/x) I_\nu + I_\nu'
c                 I_{\nu-1}'=((\nu-1)/x) I_{\nu-1} + I_\nu
c According to (9.6.6) of \ct{AS}: I_{-\nu}(z)=I_{\nu}(z);
c K_{-\nu}(z)=K_{\nu}(z).
C The spherical Bessel functions for complex argument are found
C according to (10.2.2-3) of \ct{AS}:
c
C         j_l(ix)=e^{l\pi i/2} \sqrt{\fr{\pi}{2x}}\, I_{l+1/2}(x)
C      y_l(ix)=e^{-3(l+1)\pi i/2} \sqrt{\fr{\pi}{2x}}\, I_{-l-1/2}(x)
C*****************************************************************
C                >>> FOR GENERAL COMPLEX OMEGA   <<<
C AMOS SLAC package from http://www.netlib.org/amos/ is used
C -----------------------------------------------------------------
      implicit none
      integer LMAX,IKODE,i,l,LMX
      REAL*8 pi,pib2z,prod,fnu
      PARAMETER (PI=3.141592653589793d0)
* LMX is internal LMAX here (for AMOS):
      PARAMETER (LMX=82)
* IKODE=1 ... unscaled Bessel functions
* IKODE=2 ... scaled Bessel functions (for AMOS):
      PARAMETER (IKODE=1)
* FNU=ORDER OF INITIAL Y FUNCTION, FNU.GE.0.0D0(for AMOS):
      PARAMETER (FNU=0.5d0)
*
      integer nz,ierr
      REAL*8 xom,xnu,rj,rymu,rjp,ry1,xi
      complex*16 JL(0:lmax),NL(0:lmax),DJL(0:lmax),DNL(0:lmax)
      complex*16 ei,omega,xjl,cxi,czi,cpib2z
      complex*16 crj,crymu,crjp,cry1
      real*8 zr,zi
      real*8 cyr(lmx+1),cyi(lmx+1)
c      character*1 yesn

      ei=dcmplx(0.d0,1.d0)
      zi=dimag(omega)
      zr=dble(omega)

*-------------
*  remainders:

      if (lmx.lt.lmax) then
      write(6,*)'In gnzbess.f:'
      write(6,*)'Internal LMAX (LMX).lt.LMAX'
      stop
      end if

      if ((zi.ne.0.d0).and.(zr.ne.0.d0)) then
c      write(6,*)'OMEGA is neither real nor purely imaginary'
      go to 200
      end if
*--------------
c      yesn='n'
      xom=zabs(omega)
      pib2z=sqrt(pi/(2.d0*xom))
      xnu=dble(lmax+1/2.)

*
      if (dimag(omega).ne.0.) go to 100
*
*******************************************************************
C                       >>> OMEGA  REAL  <<<
*                  >>>  ASYMPTOTIC   PART  <<<
*
* security for omega exceedingly small:

      if (zabs(omega).gt.1.d-9) go to 10

* using asymptotic expansion [cf. (10.1.2) of \ct{AS}] to determine
* JL(LMAX) and DJL(LMAX):

      prod=1.d0
        do i=0,LMAX
         prod=dble(2*i+1)*prod
        enddo
         JL(LMAX)=(1.d0-omega**2/dble(2*(2*LMAX+3)) )*omega**LMAX
     &   /prod
         XJL=(1.d0-omega**2/dble(2*(2*LMAX+5)) )*omega**(LMAX+1)
     &   /(prod*dble(2*lmax+3))
         xnu=dble(lmax)
         xi=1.d0/xom
         DJL(LMAX)=JL(LMAX)*xnu*xi-XJL
      do i=1,lmax
      prod=prod/dble(2*(LMAX-i)+3)
        JL(LMAX-i)=(1.d0-omega**2/dble(2*(2*(LMAX-i)+3)))*
     &   omega**(LMAX-i)/prod
        xnu=dble(lmax-i)
        DJL(LMAX-i)=JL(LMAX-i)*xnu*xi-JL(LMAX+1-i)
      enddo

* continuation of the previous loop
      NL(0)=JL(0)*xi+DJL(0)
      DNL(0)=-NL(0)*xi-JL(0)
      NL(0)=-NL(0)
      DNL(0)=-DNL(0)
      do i=1,lmax
        xnu=dble(i-1)
        NL(i)=NL(i-1)*xnu*xi-DNL(i-1)
        xnu=xnu+2.d0
        DNL(i)=-NL(i)*xnu*xi+NL(i-1)
      enddo

      return
**************************************************************
C                      >>>  OMEGA  REAL  <<<
*                    >>>   NORMAL  PART   <<<
*
 10   call bessjy(xom,xnu,rj,rymu,rjp,ry1)

      JL(LMAX)=pib2z*rj
      DJL(LMAX)=pib2z*rjp-JL(LMAX)/(2.d0*xom)

* Stable downwards recurrence for J_\nu and
c Using (9.1.2) of \ct{AS} one has
c         (\cos[(n+1/2)\pi]=0,       \sin[(n+1/2)\pi]=(-1)^n)
c
c               Y_{n+1/2} = (-1)^{n+1} J_{-(n+1/2)}
c  and hence
c                   y_{n} = (-1)^{n+1} j_{-n-1}
c
c  and one can continue with a stable recurrence downwards to generate
c  also the relevant Neumann (Weber) functions.

      xi=1.d0/xom
      do i=1,lmax
        xnu=dble(lmax-i+2)
        JL(LMAX-i)=JL(LMAX+1-i)*xnu*xi+DJL(LMAX+1-i)
        xnu=xnu-2.d0
        DJL(LMAX-i)=JL(LMAX-i)*xnu*xi-JL(LMAX+1-i)
      enddo

* continuation of the previous loop
      NL(0)=JL(0)*xi+DJL(0)
      DNL(0)=-NL(0)*xi-JL(0)
      NL(0)=-NL(0)
      DNL(0)=-DNL(0)
      do i=1,lmax
        xnu=dble(i-1)
        NL(i)=NL(i-1)*xnu*xi-DNL(i-1)
        xnu=xnu+2.d0
        DNL(i)=-NL(i)*xnu*xi+NL(i-1)
      enddo
* security trap:
      do L=0,LMAX
        if (abs(JL(L)).gt.pib2z) then
        write(6,*)'GNCBESS gives incorrect Bessel functions'
        write(6,*)'L=', L, ', JL(L)=', JL(L),'.gt.',pib2z
        write(6,*)'Violation of the bound (9.1.60) of {AS}'
        stop
        end if
      enddo
      return
**************************************************************
*
 100  continue
*
C                   >>>  PURELY IMAGINARY OMEGA   <<<
*                     >>>   ASYMPTOTIC   PART  <<<
*
* security for omega exceedingly small:

      if (zabs(omega).gt.1.d-9) go to 110

* using asymptotic expansion [cf. (10.1.2) of \ct{AS}] to determine
* JL(LMAX) and DJL(LMAX):

      prod=1.d0
        do i=0,LMAX
         prod=dble(2*i+1)*prod
        enddo
         JL(LMAX)=(1.d0-omega**2/dble(2*(2*LMAX+3)) )*omega**LMAX
     &   /prod
         XJL=(1.d0-omega**2/dble(2*(2*LMAX+5)) )*omega**(LMAX+1)
     &   /(prod*dble(2*lmax+3))
         xnu=dble(lmax)
         cxi=-ei*1.d0/xom
         DJL(LMAX)=JL(LMAX)*xnu*cxi-XJL
      do i=1,lmax
      prod=prod/dble(2*(LMAX-i)+3)
      JL(LMAX-i)=(1.d0-omega**2/dble(2*(2*(LMAX-i)+3)))*
     &   omega**(LMAX-i)/prod
        xnu=dble(lmax-i)
        DJL(LMAX-i)=JL(LMAX-i)*xnu*cxi-JL(LMAX+1-i)
      enddo

* continuation of the previous loop
      NL(0)=JL(0)*cxi+DJL(0)
      DNL(0)=-NL(0)*cxi-JL(0)
      NL(0)=-NL(0)
      DNL(0)=-DNL(0)
      do i=1,lmax
        xnu=dble(i-1)
        NL(i)=NL(i-1)*xnu*cxi-DNL(i-1)
        xnu=xnu+2.d0
        DNL(i)=-NL(i)*xnu*cxi+NL(i-1)
      enddo
      return

*----------------------------------------------------------------------
C                   >>>  PURELY IMAGINARY OMEGA   <<<
*                     >>>   NORMAL   PART  <<<
*
 110  call bessik(xom,xnu,rj,rymu,rjp,ry1)

      JL(LMAX)=pib2z*rj*ei**lmax
      DJL(LMAX)=-ei*pib2z*rjp*ei**lmax+ei*JL(LMAX)/(2.d0*xom)

      cxi=-ei*1.d0/xom
      do i=1,lmax
        xnu=dble(lmax-i+2)
        JL(LMAX-i)=JL(LMAX+1-i)*xnu*cxi+DJL(LMAX+1-i)
        xnu=xnu-2.d0
        DJL(LMAX-i)=JL(LMAX-i)*xnu*cxi-JL(LMAX+1-i)
      enddo

* continuation of the previous loop
      NL(0)=JL(0)*cxi+DJL(0)
      DNL(0)=-NL(0)*cxi-JL(0)
      NL(0)=-NL(0)
      DNL(0)=-DNL(0)
      do i=1,lmax
        xnu=dble(i-1)
        NL(i)=NL(i-1)*xnu*cxi-DNL(i-1)
        xnu=xnu+2.d0
        DNL(i)=-NL(i)*xnu*cxi+NL(i-1)
      enddo
      return

*----------------------------------------------------------------------
C                   >>>  GENERAL COMPLEX OMEGA   <<<

 200  czi=dcmplx(1.d0,0.d0)/omega
      cpib2z=sqrt(czi)*sqrt(pi/2.d0)

*                  >>>  ASYMPTOTIC   PART  <<<
*
* security for omega exceedingly small:

      if (zabs(omega).gt.1.d-9) go to 275   !250 --> Bess,
                                            !275 --> Sandia, 300 --> NR

* using asymptotic expansion [cf. (10.1.2) of \ct{AS}] to determine
* JL(LMAX) and DJL(LMAX):

      prod=1.d0
        do i=0,LMAX
         prod=dble(2*i+1)*prod
        enddo
         JL(LMAX)=(1.d0-omega**2/dble(2*(2*LMAX+3)) )*omega**LMAX
     &   /prod
         XJL=(1.d0-omega**2/dble(2*(2*LMAX+5)) )*omega**(LMAX+1)
     &   /(prod*dble(2*lmax+3))
         xnu=dble(lmax)
         DJL(LMAX)=JL(LMAX)*xnu*czi-XJL
      do i=1,lmax
      prod=prod/dble(2*(LMAX-i)+3)
        JL(LMAX-i)=(1.d0-omega**2/dble(2*(2*(LMAX-i)+3)))*
     &   omega**(LMAX-i)/prod
        xnu=dble(lmax-i)
        DJL(LMAX-i)=JL(LMAX-i)*xnu*czi-JL(LMAX+1-i)
      enddo

* continuation of the previous loop
      NL(0)=JL(0)*czi+DJL(0)
      DNL(0)=-NL(0)*czi-JL(0)
      NL(0)=-NL(0)
      DNL(0)=-DNL(0)
      do i=1,lmax
        xnu=dble(i-1)
        NL(i)=NL(i-1)*xnu*czi-DNL(i-1)
        xnu=xnu+2.d0
        DNL(i)=-NL(i)*xnu*czi+NL(i-1)
      enddo

      return

*------------------------------
*                    >>>   NORMAL   PART  <<<

250   call bess(lmax,omega,jl,djl)

      do i=lmax,1,-1
      jl(i)=jl(i-1)
      djl(i)=djl(i-1)
      enddo

      do i=lmax,lmax
        xnu=dble(lmax-i+2)
* (10.1.21) of AS:
        JL(LMAX-i)=xnu*czi*JL(LMAX-i+1) + djl(LMAX-i+1)
        xnu=xnu-2.d0
* (10.1.22) of AS:
        DJL(LMAX-i)=JL(LMAX-i)*xnu*czi-JL(LMAX+1-i)
      enddo

* continuation of the previous loop
* (10.1.21) of AS calculating of JL(-1) and DJL(-1):
      NL(0)=JL(0)*czi+DJL(0)
      DNL(0)=-NL(0)*czi-JL(0)
*  (10.1.15) of AS: NL(0)=JL(-1) : sign correction
      NL(0)=-NL(0)
      DNL(0)=-DNL(0)
      do i=1,lmax
        xnu=dble(i-1)
* (10.1.22) of AS:
        NL(i)=NL(i-1)*xnu*czi-DNL(i-1)
        xnu=xnu+2.d0
* (10.1.21) of AS:
        DNL(i)=-NL(i)*xnu*czi+NL(i-1)
      enddo
      RETURN
*                    >>>   NORMAL   PART  SANDIA <<<

 275  call zbesj(zr,zi,fnu,ikode,lmax+1,cyr,cyi,nz,ierr)
      if(ierr.gt.0) then
      write(6,*)'ERROR in ZBESJ. IERR=', ierr
      write(6,*)'IERR=0, NORMAL RETURN - COMPUTATION COMPLETED'
C--------/---------/---------/---------/---------/---------/---------/--
       if(ierr.eq.1) then
      write(6,*)'IERR=1, INPUT ERROR-NO COMPUTATION'
      else if(ierr.eq.2) then
      write(6,*)'IERR=2, OVERFLOW - NO COMPUTATION, AIMAG(Z) TOO LARGE
     1 ON KODE=1'
        else if(ierr.eq.3) then
      write(6,*)'IERR=3, CABS(Z) OR FNU+N-1 LARGE - COMPUTATION DONE
     1 BUT LOSSES OF SIGNIFCANCE BY ARGUMENT REDUCTION PRODUCE LESS THAN
     2 HALF OF MACHINE ACCURACY'
      else if(ierr.eq.4) then
      write(6,*)'IERR=4, CABS(Z) OR FNU+N-1 TOO LARGE - NO COMPUTATION
     1  BECAUSE OF COMPLETE LOSSES OF SIGNIFIANCE BY ARGUMENT REDUCTION'
      else if(ierr.eq.5) then
      write(6,*)'IERR=5, ERROR - NO COMPUTATION
     1 ALGORITHM TERMINATION CONDITION NOT MET'
      end if
      stop
      end if

      JL(LMAX)=cpib2z*dcmplx(cyr(LMAX+1),cyi(LMAX+1))
      JL(LMAX-1)=cpib2z*dcmplx(cyr(LMAX),cyi(LMAX))
* (10.1.21) of AS:
      DJL(LMAX)=JL(LMAX-1)-dble(lmax+1)*JL(LMAX)*czi
*
      do i=1,lmax
        xnu=dble(lmax-i+2)
        JL(LMAX-i)=cpib2z*dcmplx(cyr(LMAX-i+1),cyi(LMAX-i+1))
        xnu=xnu-2.d0
* (10.1.22) of AS:
        DJL(LMAX-i)=JL(LMAX-i)*xnu*czi-JL(LMAX+1-i)
      enddo

* continuation of the previous loop
* (10.1.21) of AS calculating of JL(-1) and DJL(-1):
      NL(0)=JL(0)*czi+DJL(0)
      DNL(0)=-NL(0)*czi-JL(0)
*  (10.1.15) of AS: NL(0)=JL(-1) : sign correction
      NL(0)=-NL(0)
      DNL(0)=-DNL(0)
      do i=1,lmax
        xnu=dble(i-1)
* (10.1.22) of AS:
        NL(i)=NL(i-1)*xnu*czi-DNL(i-1)
        xnu=xnu+2.d0
* (10.1.21) of AS:
        DNL(i)=-NL(i)*xnu*czi+NL(i-1)
      enddo
      RETURN
*-------------------------  UNDER CONSTRUCTION  -------------------------
*                     >>>   NORMAL   PART  NR <<<
 300      xnu=dble(lmax+1/2.)
C      call zbessjy(omega,xnu,crj,crymu,crjp,cry1)

      JL(LMAX)=cpib2z*crj
      DJL(LMAX)=cpib2z*crjp-JL(LMAX)*czi/2.d0

* Stable downwards recurrence for J_\nu and
c Using (9.1.2) of \ct{AS} one has
c         (\cos[(n+1/2)\pi]=0,       \sin[(n+1/2)\pi]=(-1)^n)
c
c               Y_{n+1/2} = (-1)^{n+1} J_{-(n+1/2)}
c  and hence
c                   y_{n} = (-1)^{n+1} j_{-n-1}
c
c  and one can continue with a stable recurrence downwards to generate
c  also the relevant Neumann (Weber) functions.

      do i=1,lmax
        xnu=dble(lmax-i+2)
        JL(LMAX-i)=JL(LMAX+1-i)*xnu*czi+DJL(LMAX+1-i)
        xnu=xnu-2.d0
        DJL(LMAX-i)=JL(LMAX-i)*xnu*czi-JL(LMAX+1-i)
      enddo

* continuation of the previous loop
      NL(0)=JL(0)*czi+DJL(0)
      DNL(0)=-NL(0)*czi-JL(0)
      NL(0)=-NL(0)
      DNL(0)=-DNL(0)
      do i=1,lmax
        xnu=dble(i-1)
        NL(i)=NL(i-1)*xnu*czi-DNL(i-1)
        xnu=xnu+2.d0
        DNL(i)=-NL(i)*xnu*czi+NL(i-1)
      enddo
* security trap: here use (9.1.62)
c      do L=0,LMAX
c        if (abs(JL(L)).gt.pib2z) then
c        write(6,*)'GNZBESS gives incorrect Bessel functions'
c        write(6,*)'L=', L, ', JL(L)=', JL(L),'.gt.',pib2z
c        write(6,*)'Violation of the bound (9.1.62) of {AS}'
        write(6,*)'Bound (9.1.62) of {AS} not checked'
c        stop
c        end if
c      enddo
      return

* End of evaluation
      END
C (C) Copr. 7/2000  Alexander Moroz

      subroutine intens(lmax,lcs,rsnm,lambda,thetv,rmf,zeps)
C--------/---------/---------/---------/---------/---------/---------/--
C >>> lmax,lcs,rsnm,lambda,thetv,rmf,zeps
C <<<
C =============
C f77 -g -C intens.f -o rnintens
C =============
C
C    PROGRAM TO DETERMINE TOTAL ELECTRIC FIELD FOR A PLANE WAVE
C                    INCIDENCE ALONG THE Z-AXIS
C
C >>> INPUT VARIABLES:
C
C    rsnm ... nanoshell radius
C    rmf  ... array of shell radiii in units of rsnm
C    rsnm - characteristic particle dimensions
C    LAMBDA - vacuum wavelength
C    THETV - zenith (theta) angle of the incident beam in degrees
C===========================
C    OTHER VARIABLES&ARRAYS:
C
C    zeinc(3,2)  ... incident plane wave electric field
C    zescat(3,2) ... scattered wave electric field
C                First column labels 1,2,3 correspond to the respective
C                     r, theta, and phi vector components
C                of the respective incident and scattered electric fields.
C                Second column label 1 corresponds to the
C                  (\vx+i\vy) circular polarization component,
C                whereas second column label 2 corresponds to the
C                  (\vx-i\vy) circular polarization component,
C                of the respective incident and scattered electric fields.
C
C     xpole(NGD,3) ... various multipoles at NGD different positions on the
C                      unit sphere in the spherical coordinates.
C              xpole(*,1) is the \ve_\vr component
C              xpole(*,2) is the \ve_\theta component
C              xpole(*,3) is the \ve_\vphi component
C
C     LMAXD ... maximal angular momentum cutoff. If changing LMAXD, you
C               should also change it in GNRICBESSH, BIGA, SPHRINT,
C               VCTHARM, VIGAMPL
C     LMAX ... floating angular momentum cutoff for plane wave expansion
C     LMAXS ... angular momentum cutoff for T-matrix
C                 LMAXS <= LMAX <= LMAXD
C     NGAUSSD ... maximal number of allowed Gaussian integration points
C     NGAUSS  ... floating  number of different theta points
C     NGD    ... maximal number of Gaussian integration points
C     NG,NGV=NG+3 ... floating number of different theta points; NG goes
C                 to CONSTE, NGV then to RSPP routines
C     NRM   ... maximal number of shels around a particle
C     NRMAX ... maximal number of allowed different R points
C--------/---------/---------/---------/---------/---------/---------/--
      implicit none
      integer lmaxd,lmaxd1,LMTD,NRM,NPHI,NOUT,NX
      real*8 TOL,RINT,rsnm,lambda,thetv,thet

* number of the output unit
      PARAMETER (NOUT=60)
* angular-momentum cut-offs
      PARAMETER (LMAXD=60,LMAXD1=LMAXD+1,LMTD=LMAXD1*LMAXD1-1)
* cut-offs on the number of different theta, phi, and R-shells
      PARAMETER(NRM=3)

C ::: number of scan steps along radial direction - fill in odd number:
      PARAMETER (NX=201)

C ::: radial interval on which the decay rates are calculated in the units
C     of the sphere radius (RMUF=1)
      PARAMETER (RINT=2.d0)
*
C ::: relative error allowed. If the convergence
*     within TOL is not reached, program issues warning
      PARAMETER (TOL=1.d-2)
*
      integer i,ij,ip,il,j,l,ilm,ilmp,lcs1,nscan
      INTEGER LCS,lmax,lmaxs,lmaxt,m,ml,NGV
c      logical ynphic

      real*8 pi,omega,xelinv,delo,delt,rmff,rmfc
      REAL*8 phi,xint,aint,xx
      REAL*8 rmf(lcs),X(100),r(100)

      COMPLEX*16 RX(2),SG(2),CQEPS(2)                  !ZSIGMA
      complex*16  cdl,zrfac,zsfac,zint(3),zeps(lcs+1)
      complex*16  ci,czero,zeinc(3,2),zescat(3,2)
      complex*16  zmfac,zmpifac,zpifac,zelinc(2)
      complex*16  cpole(lmtd,3),bpole(lmtd,3),ppole(lmtd,3)
      complex*16  cpolinc(lmtd,3),bpolinc(lmtd,3),ppolinc(lmtd,3)
      complex*16  mrpole(lmtd,3),nrpole(lmtd,3),mspole(lmtd,3),
     & nspole(lmtd,3),mrpolinc(lmtd,3),nrpolinc(lmtd,3)
      COMPLEX*16 ZA(LMTD),ZB(LMTD)
      COMPLEX*16 AM(lmaxd,lcs+1),AE(lmaxd,lcs+1),BM(lmaxd,lcs+1),
     & BE(lmaxd,lcs+1)
C                      =====================

      COMPLEX*16 JL(0:lmaxd), NL(0:lmaxd)
      COMPLEX*16 DRJL(0:lmaxd), DRNL(0:lmaxd)
      COMPLEX*16 HL(0:lmaxd),DRHL(0:lmaxd)
*
C--------/---------/---------/---------/---------/---------/---------/--
C  INPUT DATA ********************************************************
*
      DATA PI/3.141592653589793d0/
      DATA ci/(0.d0,1.d0)/,czero/(0.D0,0.D0)/
        lcs1=lcs+1

* specify the number of different PHI points around the axis of symmetry
      NPHI=40
        NGV=40       !number of different theta points

        if(ngv.gt.100) write(6,*)
     & 'Increase dimension of R- and X-arrays to at least NGV+1=',
     & NGV+1
*
*
* test setup:

      if (LMAX.gt.LMAXD) then
      write(6,*)'In INTENS'
      write(6,*)'LMAX has to be smaller than LMAXD in INTENS'
      stop
      end if

C  TEMPORARY OPTIONS ******************************************

      lmaxs=lmax

ct      LMAXS=max(6.d0,1+2*pi*rsnm/lambda)   !cutoff on scattered field

        if(lmaxs.gt.lmax) then
        write(6,*)'In INTENS:'
        write(6,*)'Increase LMAX so that LMAXS.LE:LMAX'
        pause
        stop
        end if
*
        write(nout,*)  '#Cutoff for T-matrix LMAXS=',LMAXS
        write(nout+1,*)'#Cutoff for T-matrix LMAXS=',LMAXS
        write(nout+2,*)'#Cutoff for T-matrix LMAXS=',LMAXS
*

C  END OF TEMPORARY OPTIONS ******************************************
cc      write(6,*)'r-component of the incident field amplitude'
cc      (5,*) zelinc(1)
      write(6,*)
cd      write(6,*)'FILL IN THE INCIDENT ELECTRIC FIELD COMPONENTS'
cd      write(6,*)
cd      write(6,*)'THETA COMPONENT:'
cd      read(5,*) zelinc(1)
cd
cd      write(6,*)'PHI COMPONENT:'
cd      read(5,*) zelinc(2)

      write(6,*)'LINEARLY POLARIZED INCIDENT WAVE OF UNIT AMPLITUDE'
      write(6,*)'FILL IN POLARIZATION ANGLE ALPHA IN DEGREES'
        write(6,*)'[between (k,z)-plane and polarization vector]'
      write(6,*)'(For incidence along z-axis: '
      write(6,*)'alpha=0 for x-polarization'
        write(6,*)'alpha=90 for y-polarization)'
        read(5,*) xx
        write(nout,*)
     &'#polarization angle alpha w.r.t. meridional plane in deg=',xx
        write(nout+1,*)
     &'#polarization angle alpha w.r.t. meridional plane in deg=',xx
        write(nout+2,*)
     &'#polarization angle alpha w.r.t. meridional plane in deg=',xx

        xx=xx*pi/180.d0

* Theta component - per definition lies in the meridional plane
* i.e., the plane through k and the z-axis. Provided that k coincides
* with the z-axis, the meridional plane is taken to be the (z,x)-plane
        zelinc(1)=cos(xx)

* Phi component- per definition perpendicular to the meridional plane
        zelinc(2)=sin(xx)
*
      xelinv=sqrt(dble( zelinc(1)*dconjg(zelinc(1))+
     1      zelinc(2)*dconjg(zelinc(2)) ) )
*
      thetv=thetv*pi/180.d0       !Theta angle of incidence (in radians)
      xx=cos(thetv)
*
* Constants initialization
C RAP=S(1,1)*KAPPA0/2.D0/PI=rmuf*ALPHA/LAMBDA =rsnm/LAMBDA

  2   write(6,*)
        write(6,*)'For semi-circular theta-scan write 1'
        write(6,*)'For scan along the x-axis write 2'
      write(6,*)'For scan along the z-axis write 3'
      write(6,*)'(by default, incidence in the (x,z)-plane)'
      read(5,*) nscan
      write(6,*)

      IF ((nscan.lt.1).or.(nscan.gt.4)) then
        write(6,*)'Incorrect choice of nscan. Please again.'
        go to 2
        end if

        IF (nscan.eq.1) THEN

      write(6,*)
      write(6,*)'Chose a particular value of out-PHI in degrees'
      read(5,*) phi
      phi=phi*pi/180.d0     !conversion to radians
        NPHI=1

        ELSE IF (nscan.eq.2) THEN

cc      write(6,*)
cc      write(6,*)'Chose a particular value of out-THETA in degrees'
cc      read(5,*) thet
        thet=90.d0
      thet=thet*pi/180.d0                         !conversion to radians
        x(1)=cos(thet)
        phi=0.d0
      NGV=0

        ELSE IF (nscan.eq.3) THEN

        phi=0.d0     !can in principle be arbitrary value
                   !The result should not depend on the phi-choice
      NGV=0
        NPHI=1

        END IF
*
      zpifac=zexp(ci*phi)
        omega=2.d0*pi*rsnm/lambda   !size parameter
*
        write(nout,*)  '#In radians: Phi_inc=0, Phi_out=', phi
        write(nout+1,*)'#In radians: Phi_inc=0, Phi_out=', phi
        write(nout+2,*)'#In radians: Phi_inc=0, Phi_out=', phi
        WRITE(NOUT,*)'#cos(theta),r(theta), and |E| in columns'
        WRITE(NOUT+1,*)'#cos(theta),r(theta), E_r,E_theta,E_phi
     & in columns'
        WRITE(NOUT+2,*)'#cos(theta),r(theta), E_r,E_theta,E_phi
     & in columns'
        write(nout,*)
        write(nout+1,*)
        write(nout+2,*)
***************************************************************************
* Incident plane wave expansion into spherical harmonics:

      call vctharm(LMAX,xx,bpolinc,cpolinc,ppolinc)
*
* If bpolinc(*,1) or cpolinc(*,1) is nonzero then there is
* a mismatch of array dimensions between main and vctharm
*
* Scalar product for expansion coefficients - see Eqs. (3.13-14) of {Mis91}
* or Eqs. (C.56-59) of {MTL}
* Using that C^*_{mn} = (-1)^m C_{-m n}, B^*_{mn} = (-1)^m B_{-m n},
* P^*_{mn} = (-1)^m P_{-m n}
*
      do 10 l=1,lmax

      cdl=dsqrt( dble(2*l+1)/(4.d0*pi*dble(l*(l+1))) )

      do 10 m=-l,l

      ilm=l*(l+1)+m         ! (l,m)  index with (1-1)=1
      ilmp=l*(l+1)-m        ! (l,-m) index with (1-1)=1

      za(ilm)=czero
      zb(ilm)=czero

* Scalar product for expansion coefficients - see Eqs. (3.13-14) of {Mis91}
* or Eqs. (C.56-59) of {MTL}
* (using that C^*_{mn} = (-1)^m C_{-m n}, B^*_{mn} = (-1)^m B_{-m n})

      do 5 ij=2,3              !1st (r-)components of C and B poles are zero
           za(ilm)=za(ilm)+cpolinc(ilmp,ij)*zelinc(ij-1)
           zb(ilm)=zb(ilm)+bpolinc(ilmp,ij)*zelinc(ij-1)

  5   continue

      za(ilm)=4.d0*pi*cdl*za(ilm)*ci**l
      zb(ilm)=4.d0*pi*cdl*zb(ilm)*ci**(l-1)    !in-phi set to zero by default

* In Eqs. (3.13-14) of {Mis91} the factor (-1.d0)**m in the expansion
* coefficiens  cancels out against the same factor (-1.d0)**m when the
* complex conjugate C and B poles are rewritten in terms of the
* nonconjugate ones (see index ilmp above) using formulae
*       C^*_{mn} = (-1)^m C_{-m n}, B^*_{mn} = (-1)^m B_{-m n}

 10   continue
C--------/---------/---------/---------/---------/---------/---------/--
*
* Calculation of the intensities in respective shells:
*
      call sphrint(lmax,lcs,lambda,rsnm,rmf,zeps,am,ae,bm,be)
C--------/---------/---------/---------/---------/---------/---------/--
C  Given the incident unit plane wave expansion coefficients C, returns
C  the arrays A, B for each shell of a multilayered sphere
C--------/---------/---------/---------/---------/---------/---------/--

*>>> scanning radial path over the interval RINT:

c      rmfc=1.d0       !determines a lower cutoff from which
                       !incident wave amplitude is calculated
                               !for consistency check - should be at least 1.d0

      do 600 il=0,nx            !nx            !over the interval RINT

c      do 500 ipi=0,NPHI                      !over different phi angles
cc      phi=2.d0*pi*dble(ipi)/dble(NPHI)

        if (nscan.le.2) then

        delo=rint/dble(nx)
      rmff=1.d-4 + dble(il)*delo


        else if (nscan.eq.3) then

        delo=2.d0*rint/dble(nx)
      rmff=-rint+ 1.d-4 + dble(il)*delo

* For theta=\pi/2 incidence this should read
c     x(1)=xx
c      if (rmff.lt.0.d0)  phi = -pi
c      if (rmff.ge.0.d0)  phi = 0.d0
c
* For theta=0 incidence this should read
c     phi=0.d0
c
        if (rmff.lt.0.d0)  x(1)=-1.d0          !=cos(theta)
        if (rmff.ge.0.d0)  x(1)= 1.d0          !=cos(theta)
      rmff=abs(rmff)

      end if            !nscan

        if (nscan.eq.2) x(1)=0.d0

        if (rmff.lt.1.d-4) rmff=1.d-4

cd        rmff=(rsnm+1.d0)/rsnm
cx      rmff= 0.8771d0 + dble(il)*delo

* identifying the layer for a given rmff

      do 45 ml=1,lcs

      if (rmff.lt.rmf(ml)) then
      j=ml
      go to 48
      end if

 45   continue

      if (rmff.gt.rmf(lcs)) j=lcs1

 48   CQEPS(1)=SQRT(ZEPS(j))
      SG(1)=OMEGA*CQEPS(1)      !wave vector in the dipole medium
      RX(1)=SG(1)*rmff

* generate arrays of Bessel functions
*
      lmaxt=lmax
        if (rmff.lt.1.d-2) lmaxt=min(20,lmax,lmaxs)   !cut-off to control an
                                                    !under/overflow of Bessel functions
*
      CALL GNZBESS(RX(1),LMAXT,jl,drjl,nl,drnl)
*
      if (j.gt.1) CALL GNRICBESSH(CQEPS(1),OMEGA*rmff,LMAXS,hl,drhl)
C--------/---------/---------/---------/---------/---------/---------/--
C Returns Riccati-Bessel functions zeta,dzeta of the argument
C RX(i)=CQEPS(i)*OMEGA*rmff(j)
C--------/---------/---------/---------/---------/---------/---------/--

      DO L=1,LMAXT

      if ((j.gt.1).and.(l.le.lmaxs))  HL(L)=HL(L)/RX(1)      !HL now h_l(u)
        if ((j.gt.1).and.(l.le.lmaxs))  DRHL(L)=DRHL(L)/RX(1)  !DRHL now [uh_l(u)]'/u

        DRJL(L)=JL(L)/RX(1) + DRJL(L)           !DRJL now [uj_l(u)]'/u

        ENDDO
*
* Bessel function determined ===>
*          determine electric and magnetic multipoles:

      do 200 i=1,NGV+1                     !over different theta�s

      if (nscan.eq.1) then

      delt=pi/dble(ngv)
      x(i)=cos(dble(i-1)*delt)

      end if
*********************************************************
* Forming M and N poles as in Eqs. (3.3-4) of {Mis91}
* Initialization of the vector spherical harmonics:
*
      call vctharm(lmax,x(i),bpole,cpole,ppole)
*
* If bpole(*,1) or cpole(*,1) is nonzero then there is
* a mismatch of array dimensions between main and vctharm
*

*
      do 55 l=1,lmaxt

      cdl=dsqrt( dble(2*l+1)/(4.d0*pi*dble(l*(l+1))) )

      zrfac=cdl*jl(l)

      if ((j.gt.1).and.(l.le.lmaxs)) zsfac=cdl*hl(l)

      do 55 m=-l,l

      zmfac=(-1.d0)**m
*
      ilm=l*(l+1)+m         ! (l,m) index with (1-1)=1
*
      do 50 ij=1,3                   !over vector components

      if (rmff.gt.1.d0) then

* incident multipoles:  (PHI=0 for incident wave)
      mrpolinc(ilm,ij)=zmfac*zrfac*cpolinc(ilm,ij)
      nrpolinc(ilm,ij)=(zrfac*ppolinc(ilm,ij)*dble(l*(l+1))/rx(1)+
     & cdl*drjl(l)*bpolinc(ilm,ij))*zmfac
      end if

* scattered multipoles:

      mrpole(ilm,ij)=zmfac*zrfac*cpole(ilm,ij)*zpifac**m
*
      if ((j.gt.1).and.(l.le.lmaxs))
     &  mspole(ilm,ij)=zmfac*zsfac*cpole(ilm,ij)*zpifac**m
*
*
      nrpole(ilm,ij)=(zrfac*ppole(ilm,ij)*dble(l*(l+1))/rx(1)+
     & cdl*drjl(l)*bpole(ilm,ij))*zmfac*zpifac**m
*
      if ((j.gt.1).and.(l.le.lmaxs))
     & nspole(ilm,ij)=
     &     (zsfac*ppole(ilm,ij)*dble(l*(l+1))/rx(1)+
     &                cdl*drhl(l)*bpole(ilm,ij))*zmfac*zpifac**m
*
* M and N poles are as in Eqs. (3.3-4) of {Mis91}.

 50   continue             !over vector components
 55   continue             !over angular-momentum indices

*************************                  *********************
*
* Initialization of the incident & scattered fields:
*
      do 70 ip=1,1
        do 70 ij=1,3
         if (rmff.gt.1.d0) zeinc(ij,ip)=czero
          zescat(ij,ip)=czero
 70   continue

      if (rmff.gt.1.d0) then
* Recalculate incident field for consistency check:

      do 100 l=1,lmaxt
      do 100 m=-l,l

      ilm=l*(l+1)+m        ! (l,m) index with (1-1)=1

cc      zmpifac=zpifac**m     !zpifac=zexp(ci*phi) here

      do ij=1,3

* incident wave:

       zeinc(ij,1)=zeinc(ij,1)
     1  +za(ilm)*mrpolinc(ilm,ij)+zb(ilm)*nrpolinc(ilm,ij)

cc Plane-wave incidence along z-direction:
cc       zeinc(ij,1)=zeinc(ij,1)
cc     1        +zfac(l)*(mrpole(ilm+1,ij)+nrpole(ilm+1,ij))*zpifac
cc       zeinc(ij,2)=zeinc(ij,2)
cc     1        +zfac(l)*(mrpole(ilm-1,ij)-nrpole(ilm-1,ij))/zpifac
cc
      enddo
*
 100   continue                    !over angular momenta
      end if                       !rmff.gt.1.d0
*******************************************************************
* Local field intensity:

      do 130 l=1,lmaxs
      do 130 m=-l,l

      ilm=l*(l+1)+m            ! (l,m)  index with (1-1)=1

cc      zmpifac=zpifac**m        !zpifac=zexp(ci*phi) here

      do ij=1,3

      if (j.eq.lcs1) then
       zescat(ij,1)=zescat(ij,1)
     1   + za(ilm)*bm(l,j)*mspole(ilm,ij)+
     2              zb(ilm)*be(l,j)*nspole(ilm,ij)
      else if (j.eq.1) then
       zescat(ij,1)=zescat(ij,1)
     1   + za(ilm)*am(l,j)*mrpole(ilm,ij)+
     2              zb(ilm)*ae(l,j)*nrpole(ilm,ij)
      else
       zescat(ij,1)=zescat(ij,1)
     1   +za(ilm)*am(l,j)*mrpole(ilm,ij)+
     2              zb(ilm)*ae(l,j)*nrpole(ilm,ij)+
     2              za(ilm)*bm(l,j)*mspole(ilm,ij)+
     2              zb(ilm)*be(l,j)*nspole(ilm,ij)
      end if

      enddo

 130   continue                    !over angular momenta

************************************************************
* Final field intensity:

      do ip=1,1
      xint=0.d0
      aint=0.d0

         do ij=1,3

      if (rmff.gt.1.d0) aint=aint+zeinc(ij,ip)*dconjg(zeinc(ij,ip))

         if (j.eq.lcs1) then
           zint(ij)=zeinc(ij,ip)+zescat(ij,ip)
         else
           zint(ij)=zescat(ij,ip)
         end if

         xint=xint + zint(ij)*dconjg(zint(ij))

         enddo              !ij-loop

      if (rmff.gt.1.d0)  aint=sqrt(aint)
         xint=sqrt(xint)

      if (nscan.le.2) r(i)=rsnm*rmff
      if (nscan.eq.3) r(i)=x(i)*rsnm*rmff

      write(6,*)'r=', r(i)
      write(6,*)'phi=',phi
      write(6,*)'cos(theta)=',x(i)

      if (rmff.gt.1.d0) then
      if (abs((aint-xelinv)/xelinv).gt.5.d-2) then
      write(6,*)'INCREASE LMAX IN AXSPARTCL!'
      write(6,*)'LMAX is not enough to ensure the convergence
     1 of the incident field intensity!!!'
      pause
      end if
      end if

      if (rmff.gt.1.d0)
     &  write(6,*)'Incident el. field magnitude=',aint
      write(6,*)'Total el. field magnitude=',xint
      if (rmff.gt.1.d0)
     &  write(6,*)'El. field enhancement factor=', xint/aint

cc      if (ip.eq.1) write(6,*)'El. intensity for plus circ. pol.=',xint
cc      if (ip.eq.2) write(6,*)'El. intensity for minus circ. pol.=',xint
C--------/---------/---------/---------/---------/---------/---------/--

      write(nout,1000)   x(i),r(i),xint
      write(nout+1,1010) x(i),r(i),
     1 sqrt(dble(zint(1)*dconjg(zint(1)))),
     2 sqrt(dble(zint(2)*dconjg(zint(2)))),
     3 sqrt(dble(zint(3)*dconjg(zint(3))))
* sqrt(dble(zint*dconjg(zint)))
* (1),zint(2),zint(3)
      write(nout+2,1010) x(i),r(i),
     1 sqrt(dble(zescat(1,1)*dconjg(zescat(1,1)))),
     2 sqrt(dble(zescat(2,1)*dconjg(zescat(2,1)))),
     3 sqrt(dble(zescat(3,1)*dconjg(zescat(3,1))))
* (1),zescat(2),zescat(3)


cc      if (ip.eq.2) write(nout+1,*) phi,x(i),rsnm,xint
      enddo                !ip-loop

      if (rmff.lt.1.d-3) go to 500

 200   continue                    !over theta angle

 500   continue                    !over different phi angles

 600   continue                    !end of scanning path
*
 1000 FORMAT (2F10.4,F16.8)
 1010      FORMAT (3F10.4,3F14.8)

      return
      end

C**********************************************************************

      subroutine vctharm(lmax,x,bpole,cpole,ppole)
C--------/---------/---------/---------/---------/---------/---------/--
C >>> x=cos(theta),bpole,cpole,ppole
C                   (r,theta,phi)-components in the order
C =============
C
C    PROGRAM TO DETERMINE THE REDUCED VECTOR SPHERICAL HARMONICS B,C,P
C    that are only a function of theta as defined by
C    to Eqs. (3.5-8) of \ct{Mi91}, resp. to Eqs. (C.19-21) of \ct{MTL}
C
C     LMAXD ... maximal angular momentum cutoff
C     LMAX ... floating  angular momentum cutoff
C--------/---------/---------/---------/---------/---------/---------/--
      implicit none
      integer lmaxd,lmaxd1,Lmtd
      PARAMETER (LMAXD=60,LMAXD1=LMAXD+1,LMTD=LMAXD1**2-1)

      integer ilm,ilp,l,lmax,m,mp
      REAL*8 DDV1(LMAXD1),DV1(LMAXD1),DV2(LMAXD1),X
      complex*16  cpole(lmtd,3),bpole(lmtd,3),ppole(lmtd,3)
      complex*16  ci,czero

C  INPUT DATA ********************************************************
*
      DATA ci/(0.d0,1.d0)/,czero/(0.D0,0.D0)/

*
C  TEMPORARY OPTIONS ******************************************

cc      LMAX=LMAXD

*********************************************************************
*
* determine \pi and \tau scattering functions in terms of the
* Wigner d-functions:

cc      si=dsqrt(1.d0-x*x)                 !DX=DABS(X)

      do 20 m=0,lmax

cc Activate only for testing purposes
cc      if (DABS(1D0-DABS(X)).gt.5D-2) call VIGG (X, LMAX, M, DV1, DV2)

      call vigamplv(x,LMAX,M,DV1,DV2,DDV1)
C--------/---------/---------/---------/---------/---------/---------/--
C  ===> X=cos(theta),LMAX,M (only nonnegative)
C <===  DV1,DV2,DDV1
C============================
C     DV1(N)=dvig(0,m,n,arccos x)/sin(arccos x)    ! = d_{0m}^{(l)}
C
C     DV2(N)=[d/d(arccos x)] dvig(0,m,n,arccos x)  ! = d d_{0m}^{(l)}/d\theta
C
C     DDV1(N)=m*dvig(0,m,n,arccos x)/sin(arccos x) ! = m*d_{0m}^{(l)}/ sin\theta
C     and
C
C     One can show that $d_{00}^{(1)}(\theta)=\cos\theta$
C     One has (see Eq. (4.2.5) of \ct{Ed}):
C                       $d_{0m}^{(l)}=(-1)^m d_{0-m}^{(l)}$
C     and (see Eq. (35) of \ct{Mis91}):
C            $dd_{0m}^{(l)}/(d\theta)=(-1)^m dd_{0-m}^{(l)}/(d\theta)$
C--------/---------/---------/---------/---------/---------/---------/--
*
      mp=max(m,1)
      do 10 l=mp,lmax       ! m is only positive here

      ilp=L*(L+1)+M        ! (l,m) index with (1-1)=1

      bpole(ilp,1)= czero                      !Eq. (3.5) of \ct{Mis91}
      bpole(ilp,2)= dcmplx(dv2(l))
      bpole(ilp,3)= ci*ddv1(l)

      cpole(ilp,1)= czero                      !Eq. (3.6) of \ct{Mis91}
      cpole(ilp,2)= ci*ddv1(l)
      cpole(ilp,3)=-dcmplx(dv2(l))

      ppole(ilp,1)= dcmplx(dv1(l))             !Eq. (3.7) of \ct{Mis91}
      ppole(ilp,2)= czero
      ppole(ilp,3)= czero

      if (m.eq.0) go to 10

* Assigning harmonics for negative m using that
* X_{l-m} = (-1)^m X_{lm}^*, resp. that
* d^l_{0-m} = (-1)^m d^l_{0-m}:

      ilm=L*(L+1)-M        ! (l,-m) index with (1-1)=1

      bpole(ilm,1)=czero
      bpole(ilm,2)=dcmplx(dv2(l))*(-1.d0)**m
      bpole(ilm,3)=-ci*ddv1(l)*(-1.d0)**m

      cpole(ilm,1)=czero
      cpole(ilm,2)=-ci*ddv1(l)*(-1.d0)**m
      cpole(ilm,3)=-dcmplx(dv2(l))*(-1.d0)**m

      ppole(ilm,1)=dcmplx(dv1(l)*(-1.d0)**m)
      ppole(ilm,2)=czero
      ppole(ilm,3)=czero

 10   continue
 20   continue

      return
      end

C*********************************************************************

      SUBROUTINE VIGAMPLV (X,LMAX,M,DV1,DV2,DDV1)
C--------/---------/---------/---------/---------/---------/---------/--
C >>> X,LMAX,M (only nonnegative)
C <<< DV1,DV2,DDV1
C =============
C
C     X=cos(theta), where theta is the polar angle
C     LMAXD ... maximal angular momentum cutoff
C     LMAX ... floating  angular momentum cutoff
C
C Returns \pi and \tau scattering functions in terms of the
C Wigner d-functions. Algorithm as described in Eqs. (31-35)
C  of Ref. \cite{Mis39} used. (Note however a missing $n$
C  factor in the 2nd term in the curly bracket in
C   Eq. (35) of Ref. \cite{Mis39}.)
C
C     For a given azimuthal number M.GE.0 returns
C      the Wigner d-functions
C            DV1(N)=dvig(0,m,n,arccos x) = d_{0m}^{(l)}
C
C  \pi scattering function:
C     DDV1(N)=dvig(0,m,n,arccos x)/sin(arccos x)
C                              = m*d_{0m}^{(l)}/ sin\theta
C
C  \tau scattering function:
C     DV2(N)=[d/d(arccos x)] dvig(0,m,n,arccos x)
C                              = d d_{0m}^{(l)}/d\theta
C
C     for 1.LE.N.LE.LMAX and 0.LE.X.LE.1
C      DDV1 is calculated because (DV1/sin\theta) is singular for
C             either \beta=0 or \beta=\pi
C     (For a given M.NEQ.0, only the M.LE.N.LE.LMAX terms are determined!)
C =====
C
C     In the present case (Eq. (B.28) of Ref. \ct{MTL})
C                       (cf. (4.1.24) of \ct{Ed}):
C
C             d_{00}^{(l)}(\theta)= P_l(\cos\theta)
C      d_{0m}^{(l)}(\theta)= (-1)^m \sqrt{(l-m)!/(l+m)!}P_l^m(\cos\theta)
C     (assuming Jackson's P^m_l), where $d^{(l)}_{m0}=(-1)^m d^{(l)}_{0m}$)
C     (Edmonds P^m_l is without (-1)**m prefactor; cf. (4.1.24) therein)
C
C Special values:
C
C     (Rodrigues formula [Eq. (2.5.14) of Ref. \ct{Ed}] then yields
C                       P_1(x)=x; P_2=(3x^2-1)/2; etc.
C
C      P_l^m(x) = (-1)^m (1-x^2)^{m/2} \fr{d^m P_l(x)}{dx} ===>
C                       P_1^1(\cos\beta)=-\sin\beta
C     Therefore,
C              d_{00}^{(1)}(\beta)=\cos\beta
C              d_{01}^{(1)}(\beta)=\sin\beta/\sqrt{2}
C         d d_{00}^{(1)}(\beta)/d\beta=-\sin\beta
C         d d_{01}^{(1)}(\beta)/d\beta=\cos\beta/\sqrt{2}
C
C     Acc. Eq. (34) of {Mis39}:
C
C     A_0=1, A_1=1/\sqrt{2}, A_2=\sqrt{3}/(2*\sqrt{2})
C
C     Therefore [Eq. (32) of {Mis39}]:
C              d_{00}^{(0)}(\beta)=1
C              d_{01}^{(1)}(\beta)=\sin\beta/\sqrt{2}
C              d_{02}^{(2)}(\beta)=\sqrt{3}\sin^2\beta/(2*\sqrt{2})
C     and
C         d d_{00}^{(0)}(\beta)/d\beta=0
C         d d_{01}^{(1)}(\beta)/d\beta=\cos\beta/\sqrt{2}
C         d d_{02}^{(2)}(\beta)/d\beta=\sqrt{3}\sin\beta \cos\beta/\sqrt{2}
C                                = \sqrt{3}\sin (2\beta) /(2*\sqrt{2})
C =====
C     Similar to routine VIG, which however only returns
C            DV1(N)=dvig(0,m,n,arccos x) = d_{0m}^{(l)}
C
C     When arccos x is very small, a care has to be exercise to generate
C     nboth DDV1 and DV2. That part has been made using recurrences of
C     Ref. \ct{TKS}
C--------/---------/---------/---------/---------/---------/---------/--

      IMPLICIT none
      INTEGER LMAXD,LMAXD1

      PARAMETER (LMAXD=60,LMAXD1=LMAXD+1)

      integer n,LMAX,M,I,I2
      REAL*8 A,X,QS,D1,D2,D3,DER,DN,DX,QN,QN1,QN2,
     & QNM,QNM1,QMM
      REAL*8 DDV1(LMAXD1), DV1(LMAXD1), DV2(LMAXD1)

* DV1, DDV1, and DV2 initialization
      DO 1 N=1,LMAX
         DV1(N) =0.D0
         DDV1(N)=0.D0
         DV2(N) =0.D0
    1 CONTINUE

      DX=DABS(X)
      A=1.D0                       !initial A_0
      QS=DSQRT(1D0-X*X)            !sin\theta
***********************************************************************
*                      NONZERO DV1 INITIALIZATION
*
      IF (M.NE.0) GO TO 20
*
* DDV1(N)=0.d0 [see (3.33) of {TKS}]
* D1,D2, and D3 below are the three consequent terms
*       d_{0m}^{n-1}, d_{0m}^{n}, and d_{0m}^{n+1} beginning
*       with n=m
*=============
* Recurrence initialization following d^l_{00}=P_l
*         [see Eq. (B.27) of {MTL}]
*
      D1=1.D0                  !d^0_{00}=P_0=1   (see Sec. 6.8 of {NR})
      D2=X                     !d^1_{00}=P_1=x

      DO 5 N=1,LMAX
         QN=DFLOAT(N)
         QN1=DFLOAT(N+1)
         QN2=DFLOAT(2*N+1)
         D3=(QN2*X*D2-QN*D1)/QN1       !recurrence (31) of Ref. {Mis39} for d^{N+1}_{00}
         DV1(N)=D2                     !d^N_{00}
         D1=D2                         !becomes d^{N-1}_{00} in D3
         D2=D3                         !becomes d^{N}_{00} in D3
  5   CONTINUE
*
      go to 100

***********************************************************
*                           M\neq 0 part
*  d_{0m}^m \equiv A_m*(sin\theta)**m   initialization =
*  (33) and recurrence (34) of Ref. {Mis39}

   20 CONTINUE
      DO 25 I=1,M
         I2=I*2
         A=A*DSQRT(DFLOAT(I2-1)/DFLOAT(I2))*QS  !recurrence (33,34) of Ref. {Mis39} f
   25 CONTINUE

*
* Recurrence initialization following Eqs. (32,33) of Ref. {Mis39}

      D1=0.D0                 !=DV1(M-1); see Eq. (32) of Ref. {Mis39}
      D2=A                    !=DV1(M);   see Eq. (33) of Ref. {Mis39}
      QMM=DFLOAT(M*M)
*
      DO 30 N=M,LMAX
         QN=DFLOAT(N)
         QN1=DFLOAT(N+1)
         QN2=DFLOAT(2*N+1)
         QNM=DSQRT(QN*QN-QMM)
         QNM1=DSQRT(QN1*QN1-QMM)
         D3=(QN2*X*D2-QNM*D1)/QNM1    !recurrence (31) of Ref. {Mis39} for d^{N+1}_{0M}
         DV1(N)=D2                    !d^N_{0M}
         D1=D2                        !becomes d^{N-1}_{0M} in D3,DER
         D2=D3                        !becomes d^{N}_{0M} in D3
   30 CONTINUE

      go to 100

*                        DV1 INITIALIZED
*
*             It remains to determine DDV1 and DV2
*********************************************************
*  (1-cos\theta) is very small:
*
C   For theta=0 [see Eqs. above]:
C              d_{00}^{(0)}(0)=1
C              d_{01}^{(1)}(0)=0
C              d_{02}^{(2)}(\beta)=0
C     and
C         d d_{00}^{(0)}(\beta)/d\beta=0
C         d d_{01}^{(1)}(\beta)/d\beta=1/\sqrt{2}
C         d d_{02}^{(2)}(\beta)/d\beta=0
C
C  See Eqs. (4.1-4) of \ct{Mis91}:
C
C   (m/\sin\theta) d_{0m}^l(0)=(\delta_{m\pm 1}/2) \sqrt{l(l+1)}
C      d d_{0m}^l(0)/d\beta   =(m\delta_{m\pm 1}/2) \sqrt{l(l+1)}
C
*
*  (4.2.1) of \ct{Ed}:
*   d_{0m}^{(l)}(pi) = (-1)^{l+m} \dt_{0,m}
*
*  (4.2.3) of \ct{Ed}:
*   d_{0m}^{(l)}(0) = (-1)^{m} \dt_{0,m} = \dt_{0,m}
*=======================================
*
*  If X^l_m=(m/\sin\theta) d_{0m}^{(l)}, then, according to (3.29) of {TKS}:
*
*  X^{m+1}_{m+1}=\sin\theta \sqrt{\fr{2m+1}{2m+2}}
*                           \left(\fr{m+1}{m}\right)X^{m}_{m}
*
*  According to (3.30) of {TKS}:
*  X^{m+1}_{m}= -\sqrt{2m+1}\,\cos\theta X^{m}_{m}
*
* According to (3.31) of {TKS}:
*  X^{l}_{m}=\fr{1}{\sqrt{l^2-m^2}}\,\left[(2l-1)\cos\theta
*          X^{l-1}_{m} - \sqrt{(l-1)^2-m^2}}\,\X^{l-2}_{m} \right]
*
* Initial recurrence values are X^1_1=\sqrt{2}/2 and X^l_0=0
***********************************************************************
*                   NONZERO DDV1/DV2 INITIALIZATION
*
*                          M = 0

 100  IF (M.EQ.0) THEN     !all DDV1(N)=X^l_0=0; see (3.33) of {TKS}:

* According to (3.37) of {TKS}, DV2(0)=0.d0

      DV2(1)=-QS

      IF (LMAX.GE.2) DV2(2)=3*X*DV2(1)

      IF (LMAX.LT.3) RETURN
*
      DO N=3,LMAX           !recurrence (3.36) of {TKS},
      DV2(N)=(2*N-1)*X*DV2(N-1)/(N-1)-N*DV2(N-2)/(N-1)
      ENDDO
***********************************************************************
*                           M > 0

       ELSE IF (M.GT.0) THEN
*
* >>> Determine X^m_m according to Eq. (3.29) of {TKS}:

      A=1.d0/DSQRT(2.D0)               !X^1_1=A_1

      DO I=1,M-1
      A=QS*DBLE(I+1)*DSQRT(2*I+1.d0)*A/(I*DSQRT(2*I+2.d0))
      ENDDO

* <<< A is now X^m_m; see (3.29) of {TKS}

      DDV1(M)=A
      DV2(M)=X*A                        !see (3.34) of {TKS}

* >>> Determine X^{m+1}_m:

      IF (M.EQ.LMAX)  GO TO 120

      DER=X*DSQRT(2*M+1.d0)*A          ! DER=X^{m+1}_m; see (3.30) of {TKS}
      DDV1(M+1)=DER
      DV2(M+1)=((M+1)*X*DER-A*DSQRT(2*M+1.d0))/DBLE(M)  !(3.35) of {TKS}

* >>> Determine remaining X^{l}_m's

      IF ((M+2).EQ.LMAX)  GO TO 120

       DO N=M+2,LMAX
       D3=DSQRT(DBLE(N)**2-DBLE(M)**2)
       DDV1(N)=((2*N-1)*X*DDV1(N-1) -
     &                DSQRT(DBLE(N-1)**2-DBLE(M)**2)*DDV1(N-2))/D3
                                                      !see (3.31) of {TKS}
       DV2(N)=(N*X*DDV1(N)-DDV1(N-1)*D3)/DBLE(M)      !see (3.35) of {TKS}
       ENDDO

      END IF

cv  100 IF (M.NE.1) RETURN
cv
cv      DO 110 N=1,LMAX
cv         DN=DFLOAT(N*(N+1))
cv         DN=0.5D0*DSQRT(DN)
cv         IF (X.LT.0D0) DN=DN*(-1)**(N+1)
cv         DV1(N)=DN
cv         IF (X.LT.0D0) DN=-DN
cv         DV2(N)=DN
cv  110 CONTINUE

  120 RETURN
      END
C*********************************************************************

      SUBROUTINE VIGG (X, LMAX, M, DV1, DV2)
C--------/---------/---------/---------/---------/---------/---------/--
C Original VIGAMPL from Mishchenko code. Maintained here only for
C testing purposes.
C--------/---------/---------/---------/---------/---------/---------/--
      IMPLICIT none
      INTEGER LMAXD,LMAXD1

      PARAMETER (LMAXD=60,LMAXD1=LMAXD+1)

      integer n,LMAX,M,I,I2
      REAL*8 A,X,QS,QS1,D1,D2,D3,DER,DN,DSI,DX,QN,QN1,QN2,
     & QNM,QNM1,QMM
      REAL*8 DDV1(LMAXD1), DV1(LMAXD1), DV2(LMAXD1)

      DO 1 N=1,LMAX
         DV1(N)=0D0
         DV2(N)=0D0
    1 CONTINUE
      DX=DABS(X)
      IF (DABS(1D0-DX).LE.1D-10) GO TO 100
      A=1D0
      QS=DSQRT(1D0-X*X)
      QS1=1D0/QS
      DSI=QS1
      IF (M.NE.0) GO TO 20
      D1=1D0
      D2=X
      DO 5 N=1,LMAX
         QN=DFLOAT(N)
         QN1=DFLOAT(N+1)
         QN2=DFLOAT(2*N+1)
         D3=(QN2*X*D2-QN*D1)/QN1
         DER=QS1*(QN1*QN/QN2)*(-D1+D3)
         DV1(N)=D2*DSI
         DV2(N)=DER
         D1=D2
         D2=D3
    5 CONTINUE
      RETURN
   20 QMM=DFLOAT(M*M)
      DO 25 I=1,M
         I2=I*2
         A=A*DSQRT(DFLOAT(I2-1)/DFLOAT(I2))*QS
   25 CONTINUE
      D1=0D0
      D2=A
      DO 30 N=M,LMAX
         QN=DFLOAT(N)
         QN2=DFLOAT(2*N+1)
         QN1=DFLOAT(N+1)
         QNM=DSQRT(QN*QN-QMM)
         QNM1=DSQRT(QN1*QN1-QMM)
         D3=(QN2*X*D2-QNM*D1)/QNM1
         DER=QS1*(-QN1*QNM*D1+QN*QNM1*D3)/QN2
         DV1(N)=D2*DSI
         DV2(N)=DER
         D1=D2
         D2=D3
   30 CONTINUE
      RETURN
  100 IF (M.NE.1) RETURN
      DO 110 N=1,LMAX
         DN=DFLOAT(N*(N+1))
         DN=0.5D0*DSQRT(DN)
         IF (X.LT.0D0) DN=DN*(-1)**(N+1)
         DV1(N)=DN
         IF (X.LT.0D0) DN=-DN
         DV2(N)=DN
  110 CONTINUE
      RETURN
      END

C (C) Copr. 10/2005  Alexander Moroz
      subroutine lisac(de,lcs,nmx1,nmx2,ngauss,lambda,eps,ki,kex,
     & THET0,THET,PHI0,PHI,zeps)
C--------/---------/---------/---------/---------/---------/---------/--
C      DE ... control parameter. If DE.ne.4 than NGAUSS=N1MAX*NDGS, where
C             NDGS=4 (in PUJOL)
C      If DE.eq.4, one has to supply NMX1,NMX2,NGAUSS on the input
C      NMX1,NMX2  ... angular-momentum cutoffs. Their values are transferred
C                below to those of variables N1MAX,N2MAX. Note that
C                convergence routine PUJOL may modify on the output
C                                     N1MAX,N2MAX
C
C     MTOPE=MAX(M,1) renamed, in concordance with MIshchenko code, to NMIN
C
C      NGAUSS ... the number of Gauss abscissas (quadrature points)
C      NGAUD ... internal cut-off of the number of Gauss abscissas - relevant
C               for LISAC, TNATURAL and GAUSSA
C     lambda ... vacuum wavelength
C      EPS=b/a  where a is the spheroid revoluton axis
C
C      Following the OCP=1 option:
C      KI ... size parameter of the core w.r.t. lambda in the ambient
C      KEX ... size parameter of the shell w.r.t. lambda in the ambient
C
C      Following the OCP=2 option:
C      KI ... size parameter of the core w.r.t. lambda in the ambient
C      KEX ... size parameter ratio Rint/Rext of the shell w.r.t.
C             lambda in the ambient
C
C      mr1=sqrt(zeps(1)/zeps(lcs+1)) ... relative RI of the core w.r.t. ambient
C      mr2=sqrt(zeps(1)/zeps(lcs+1)) ... relative RI of the shell w.r.t. ambient
C      MR=MR1/MR2
C
C =====
C      ZKE ... corresponding size parameter:
C           KI*MR2 = size parameter of the core w.r.t. lambda in the shell (CA=1)
C      or
C           KEX  (CA=2)
C           ===> ZKE is a kind of Mishchenko's REV
C                      equivalent-volume-sphere size parameter
C      ZKR=ZKE*EPS**(1.0D0/3.0D0)
C      ZKR=ZKR/SQRT(SENTH*SENTH+(COSTH*EPS)**2)
C      MKR=MR*ZKR
C
C      Delta ... accuracy parameter
C
C      SA is the number of angles where calculations will be made
C      SU is the maximum size of the T-matrix (=2*nmax)
C      T is a matrix of dimension (2*n1max)*(2*n1max)
C      N2max is the "dimension" needed to calculate cross sections
C      Since usually N2max<N1MAX, this will speed calculations
C
C       Adapted from the source code LISA of
C      ARTURO QUIRANTES SIERRA
C      Department of Applied Physics, Faculty of Sciences
C      University of Granada, 18071 Granada (SPAIN)
C      http://www.ugr.es/local/aquiran/codigos.htm
C      aquiran@ugr.es
C ===========
C >>> OUTPUT
C      q.dat ... contains Qext, Qsca, Qabs, efficiencies
C                     (calculated in TNATURAL)
C      failures.dat ... about the correctness of several inequalities
C                       among the scattering matrices
C      mm1.dat ... the scattering matrix elements F11 and F22
C      mm2.dat ... the scattering matrix elements F33 and F44
C      mm3.dat ... the scattering matrix elements F12 and F34
C                  All Muller elements normalized according to
C         (1/2) \int_0^\pi F_{11}(\theta) \sin\theta d\theta = 1
C
C      If one wants Cross Sections, multiply CEXT, CSCA by 2*pi/(k*k)
C                where k = 2*pi/wavelength
C      Or, to get Qext, Qsca, multiply Cext, Csca by 2/(X*X)
C               where X is the dimensionless size parameter
C
C--------/---------/---------/---------/---------/---------/---------/--
      implicit none
*
      INTEGER NGAUD,NOUT,SU,SUP,ST,SA,SS,S1,SD,ST1,ST2,SE    !,SF
      INTEGER SAA,SSS,SUU,SDD,S11,SEE,LCS
*
C      The following parameters set the memory usage
C      SU is the maximum size of the T-matrix (=2*nmax)
C      SA is the number of angles where calculations will be made
C      If you have computer memory problems, I suggest you lower the value
C      of SU.  Remember, the higher SU, the larger the maximum particle size
C      But you have to change SU and/or SA in the entire code, not just here
C      A simple "replace all" in your text editor will do
*
      PARAMETER (SU=50,SUP=0,ST=SU+SUP,SA=1802,SS=SU+SU,S1=SU+1,
     & SD=SU+SU+1,ST1=ST+1,ST2=ST+ST,SE=SD+2)    !,SF=SU+SU+SU+SU+2
      PARAMETER(NGAUD=300)   !internal cut-off of the number of Gauss abscissas
C ::: number of the output unit for cross sections and scattering matrix
      PARAMETER (NOUT=35)
*
cc      PARAMETER (LCS=2)
C      Sometimes you only need to calculate cross sections, not intensities.
C      In that case, use the line PARAMETER(SAa=1...
C      Advantage?  You will have more memory left, so SU can be higher.
C      Otherwise, for Mueller matrix calculations, do PARAMETER(SAa=SA..

c      PARAMETER(SAa=7,SSs=1,SUu=1,SDd=1,S11=1,SEe=1)
      PARAMETER(SAa=SA,SSs=SS,SUu=SU,SDd=SD,S11=S1,SEe=SE)
*
      INTEGER NOCON1,NOCON2,NOCON,NCERO,STOPE,CFA(5),NC,CF
      INTEGER M,N,N1,NPR,NSO,M1,S,NGAUSS,N1MAX,N2MAX,N1MAX1,N1MAX2
      INTEGER NMX,NMX1,NMX2,IAM,IAM1,IAM2,IBM,NTOPE,DE,CA,OCP,ILCN

cc      INTEGER IAJ,IAY,IBJ,ICJ,NX
      REAL*8 LAMBDA,EPS,KI,KEX,THET0,THET,PHI0,PHI
      REAL*8 CEXT,CSCA,QEXT,QSCA,QABS,CGN1,CGN2,DLAM
      REAL*8 DELTA,SEN,COSS,XM,CGN(200)
      REAL*8 AP1(SDd)
      REAL*8 DPRIMA,PI,PG,H,DB(st2+1)
      real*8 PX1(SSs,SAa),PX2(SSs,SAa),PX3(SSs,SAa),PX4(SSs,SAa)
      real*8 AP2(SDd),AP3(SDd),AP4(SDd),BP1(SDd),BP2(SDd)
cc      real*8       ASIM,ALFA,CTE,CZ,DN,EM,EPS1,EPS2,EPS3,KE1,KE2,KE3,
cc    &  KI1,KI2,KI3,KM,KIM,sum,lam,NH,pc1,pc2,PC,PE,QBACK
cc      real*8 AE1(SAa),AE2(SAa),AE3(SAa),AE4(SAa),BE1(SAa),BE2(SAa)
      real*8 A1(SAa),A2(SAa),A3(SAa),A4(SAa),B1(SAa),B2(SAa),ANG(SAa)
      real*8 FALLOS(SAa)
      COMPLEX*16 D3(SUu,SUu,SDd),D4(SUu,SUu,SDd),D5(SUu,SUu,SDd)
      COMPLEX*16 AX1(SUu,SUu,SDd),AX2(SUu,SUu,SDd)
      COMPLEX*16 BX1(SUu,0:SDd,SDd),BX2(SUu,0:SDd,SDd)
      COMPLEX*16 G1(SDd),G2(SDd)
      COMPLEX*16 CT,Z1,Z2,Z3,Z4,Z5,TA(2),ZKE,MR,MR1,MR2
      COMPLEX*16 T1(ST2,ST2,ST1),T(ST2,ST2,ST1)
      complex*16 G3(SDd),G4(SDd),G5(SDd),ZEPS(lcs+1)

      EXTERNAL PUJOL,TNATURAL,CGORDANR

      COMMON/CAUNO/N1MAX1            !to TNATURAL
      COMMON/CADOS/N1MAX2            !to PUJOL

      PI=3.141592653589793238462643D0
      PG=PI/180.0D0
C--------/---------/---------/---------/---------/---------/---------/--
      if (ngauss.gt.ngaud) then
      write(6,*)'NGAUSS on the input greater than the internal
     & cutoff NGAUD=', NGAUD
      write(6,*)'Increase the value of NGAUD to at least',NGAUSS
        pause
      stop
      end if
*
C      T is a matrix of dimension (2*NMX1)*(2*NMX1)
C      NMX2 is the "dimension" needed to calculate cross sections
C      Since usually N2max<NMX1, this will speed calculations
C
      if ((nmx1.gt.su).or.(nmx2.gt.su)) then
      write(6,*)'NMX1 or NMX2 on the input greater than the internal
     & cutoff SU=', SU
      write(6,*)'Increase the value of SU to at least',max(NMX1,NMX2)
        pause
      stop
      end if
*

C      mr1= Relative RI of the core (complex) w.r.t. ambient
      mr1=sqrt(zeps(1)/zeps(lcs+1))        !(1.2d0,0.01d0)

C      mr2= Relative RI of the shell (complex) w.r.t. ambient
      mr2=sqrt(zeps(lcs)/zeps(lcs+1))      !1.05d0

c >>>      The ambient RI has earlier been renormalized to 1
cx
cx      WRITE (*,*)'Please enter Delta accuracy parameter'
cx      WRITE (*,*)'(see Eq. 21 in Appl.Opt. 32, 4652-4666, 1993)'
cx      WRITE (*,*)'(1): Delta = 0.01  (higher speed)'
cx      WRITE (*,*)'(2): Delta = 0.001 (higler accuracy)'
cx      WRITE (*,*)'(3): User-selected Delta'
cx      WRITE (*,*)'(4): No delta. NMAX and NGAUSS are entered'
cx      WRITE (*,*)'     (See reference above)'
cx 10    READ (*,*) DE
*
C      T is a matrix of dimension (2*n1max)*(2*n1max)
C      N2max is the "dimension" needed to calculate cross sections
C      Since usually N2max<N1MAX, this will speed calculations
C      And NGAUSS is the number of Gauss quadrature points needed
C      for the calculation of the T-matrix elements
C      If DE=4, you can set those elements yourself
C      Otherwise, the PUJOL convergence subroutine will do it for you
C      I personally use option 2, unless I need more accurady (DE=3)
*
c      IF(DE.NE.1.AND.DE.NE.2.AND.DE.NE.3.AND.DE.NE.4) GO TO 10
      IF(DE.EQ.1) DELTA=1.0D-2
      IF(DE.EQ.2) DELTA=1.0D-3
      IF(DE.EQ.3) THEN
20        WRITE (*,*) 'Please enter a value for Delta:'
          READ (*,*) DELTA
      IF(DELTA.LE.0) GO TO 20
      END IF

      IF(DE.EQ.4) THEN
cx            WRITE (*,*) 'Enter a value for N1MAX'
cx            READ (*,*) N1MAX
cx            WRITE (*,*) 'Enter a value for N2MAX'
cx            READ (*,*) N2MAX
      IF(NMX2.GT.NMX1) NMX2=NMX1
cx            IF(N2MAX.GT.N1MAX) N2MAX=N1MAX
cx            WRITE (*,*) 'Enter a value for NGAUSS'
cx            READ (*,*) NGAUSS
ct            DELTA=0
      END IF
*
C      If you want to calculate intensities, polarizations and the like,
C      you have to set the number of angles (NTOPE)
C      in this example, theta=0,10,20...180 degrees
      NTOPE=7
      NC=1
      DO N=1,NTOPE
            ANG(N)=DFLOAT(N-1)*180.0d0/dfloat(ntope-1)
            ANG(N)=ANG(N)*PG
      END DO
*
      OPEN(1,FILE='q.dat')
      OPEN(2,FILE='failures.dat')
      OPEN(3,FILE='mm1.dat')
      OPEN(4,FILE='mm2.dat')
      OPEN(11,FILE='mm3.dat')
*
cx      WRITE (*,*)
cx      WRITE (*,*)'You can input the particle`s dimensions in two ways:'
cx      WRITE (*,*)'1) Outer and Inner size parameters (K*Rext, K*Rint)'
cx      WRITE (*,*)'2) Outer size (K*Rext) and core/particle ratio
cx     & (Rint/Rext)'
cx      WRITE (*,*) 'What`s your choice? (1/2)'
cx50      READ (*,*) OCP
cx      IF(OCP.NE.1.AND.OCP.NE.2) GO TO 50
cx      WRITE (*,*)
cx      WRITE (*,*)'Now, enter size parameters'
cx      WRITE (*,*)'Remember, these are dimensionless size parameters'
cx      WRITE (*,*)'Please refer to subroutine TNATURAL
cx     & for more information'
cx      WRITE (*,*)
cx      WRITE (*,*) 'Enter value for K*Rext (outer surface)'
cx      READ (*,*) KEX
cx      IF (OCP.EQ.1)WRITE (*,*) 'Enter value for K*Rint (inner surface)'
cx      IF (OCP.EQ.2)WRITE (*,*) 'Enter value for Rint/Rext ratio'
cx      READ (*,*) KI
cx      WRITE (*,*) 'Entre the value for the nonsphericity parameter'
cx      WRITE (*,*) 'Please refer to subroutine TNATURAL
cx     &  for more information'
cx      READ (*,*) EPS
*
      OCP=2
      NOCON=0
*
C      Now the T-matrix core (TNATURAL) has to be run twice:
C      once for the inner surface (CA=1), and one for the outer surface (CA=2)
C      If you set NMAX, NGAUSS (DE=4), the code will go right to TNATURAL
C      Otherwise, the convergence subroutine PUJOL will be called
*
C============================================================================
      do 15 ilcn=1,lcs

            CA=ilcn
            MR=sqrt(zeps(ilcn)/zeps(ilcn+1))

      if (ilcn.eq.1) then

            ZKE=KI*MR2                        !KI=K*Rint    (OCP=1)
                                            !KI=Rint/Rext (OCP=2)
                                    !Rint,Rext ... equal-volume-sphere radii
            IF(OCP.EQ.2) ZKE=ZKE*KEX          !KEX=K*Rext*n_h of outer surface
                                            !ZKE=K*Rint*n_2 ... core size parameter
                                            !       in shell medium
            N2MAX=NMX1      !m-angular momentum cutoff

      else if (ilcn.eq.2) then

c            MR=sqrt(zeps(ilcn))

            ZKE=KEX         !ZKE=KEX=K*Rext*n_h ...  size param. of outer surface
                            ! shell in ambient

      end if    !ilcn

      IF (DE.EQ.4) THEN


C   CA=1 ... TMAT for a core assuming surrounding medium
C            having the same refractive index as coating
*
      CALL TNATURAL(CA,MR,ZKE,EPS,NMX1,N2MAX,NGAUSS,CEXT,CSCA,T,T1)
C--------/---------/---------/---------/---------/---------/---------/--
C      >>> CA,MR,ZKE,EPS,NMX1,N2MAX,NGAUSS
C      <<< CEXT,CSCA,T,T1
C      CA=1 ... TNATURAL runs for the inner surface
C                  Returns TMAT for a core assuming surrounding medium
C                 having the same refractive index as coating
C      CA=2 ... TNATURAL runs for the outer surface: returns TMAT
C                   for a particle with core
C      NMX1 ... angular-momentum cutoff (not changed on the output!)
C      N2MAX ...   azimuthal parameter (not changed on the output!)
C      N1MAX1=NMX1 ... transferred in via COMMON; used only for CA=2
C--------/---------/---------/---------/---------/---------/---------/--
      N2MAX=NMX2
      N1MAX1=NMX1     !Transferred in via COMMON; used only for CA=2

*
      ELSE
*
C      This is for an automatic convergence (nonzero zeta values)
C      NOCON1,NOCON2 are convergence data (see PUJOL for more details)
*
      if (ilcn.eq.1) then

      CALL PUJOL(CA,MR,ZKE,EPS,DELTA,CEXT,CSCA,T1,NOCON,N1MAX,N2MAX)
C--------/---------/---------/---------/---------/---------/---------/--
C      >>> CA,MR,ZKE,EPS,DELTA
C      <<< CEXT,CSCA,T1,NOCON,N1MAX,N2MAX
C         (may modify N1MAX,N2MAX on the output)
C--------/---------/---------/---------/---------/---------/---------/--
*
          N1MAX1=N1MAX
          NOCON1=NOCON

      else if (ilcn.eq.2) then
*
      CALL PUJOL(CA,MR,ZKE,EPS,DELTA,CEXT,CSCA,T,NOCON,N1MAX,N2MAX)
C--------/---------/---------/---------/---------/---------/---------/--
C      >>> CA,MR,ZKE,EPS,DELTA
C      <<< CEXT,CSCA,T,NOCON,N1MAX,N2MAX
C         (may modify N1MAX,N2MAX on the output)
C--------/---------/---------/---------/---------/---------/---------/--
*
       write(6,*)'LMAX should be at least equal or greater than',N1MAX

       NOCON2=NOCON

      end if      !ilcn
*
      END IF      !DE.NE.4
*
 15      CONTINUE    !ilcn
C============================================================================
*
C      The T-matrix calculations are finished.
C      Now let's calculate the extinction, scattering and
C      absorption efficiencies.

      if (de.eq.4) then
        NMX=NMX1
      else
        NMX=N1MAX
      end if

      DLAM=LAMBDA/SQRT(ZEPS(lcs+1))

      CALL AMPMAT (NMX,DLAM,THET0,THET,PHI0,PHI,KEX,T)
C--------/---------/---------/---------/---------/---------/---------/--
C >>> NMAX,DLAM,TP,TP1,PP,PP1,T
C <<< VV,VH,HV,HH,CEXT  (written directly to output files)
C=================
C* transfers ZEPS0,REV here from the main via COMMON
C
C    GIVEN T MATRIX IN COMMON BLOCK IT CALCULATES
C     THE AMPLITUDE MATRIX IN PARTICLE FRAME
C
C    This routine closely follows exposition by
C       M. I. Mishchenko, Calculation of the amplitude matrix
C       for a nonspherical particle in a fixed orientation,
C       Appl. Opt. vol. 39, 1026-1031 (2000).
C >>>
C   NMAX - angular momentum cutoff
C   DLAM=LAMBDA/SQRT(ZEPS0)  - wavelength of incident light in the ambient
C                     (vacuum wavelength divided by SQRT(ZEPS0))
C
C   LAMBDA - vacuum wavelength. Determined as DLAM*SQRT(ZEPS0) and
C            only used here for the write out purposes
C   TP,TP1,PP,PP1 ... respective incident and scattering beam angles
C                  (in degrees) determined w.r.t laboratory frame:
C   TP (THET0 IN MAIN) - zenith angle of the incident beam in degrees
C   TP1 (THET IN MAIN) - zenith angle of the scattered beam in degrees
C   PP (PHI0 IN MAIN) - azimuth angle of the incident beam in degrees
C   PP1 (PHI IN MAIN) - azimuth angle of the scattered beam in degrees
C   KEX ... the shell size (equiv. volume sphere) parameter in the ambient
C <<<
C   VV,VH,HV,HH ... amplitude scattering matrix elements S11,S12,S21,S22
C   CEXT ... extinction cross section for a fixed particle orientation
C--------/---------/---------/---------/---------/---------/---------/--

C      Please note that Cext, Csca are NOT cross sections ===>
C      Normalization w.r.t. (2/x**2), where KEX is the shell size
C       (equiv. volume sphere) parameter in the ambient
*
      QEXT=CEXT*2.0D0/(KEX*KEX)
      QSCA=CSCA*2.0D0/(KEX*KEX)
      QABS=QEXT-QSCA
*
C      If all you need is Qext, Qsca, Qabs please add the following line:

cc      go to 444
*
*======================================================================
C                       INTENSITY CALCULATION
C      That way, you'll skip the long and tedious intensity calculations
C      On the other hand, if you want Mueller matrix elements, do nothing more
C      Still here?  Great.
C      Now to the expansion coefficients and orientation averaging
C      I suggest JOSA A 8, 871-882 (from now on, JOSA) for a clearer picture
C      Calculation of Ann'n1 (JOSA, Eqs. 4.18)
*
      do n=0,n2max+n2max
          db(n+1)=dfloat(n+n+1)
      end do
*
      DO 303 N=1,N2MAX                !always >= 1
      DO 302 N1=0,N2MAX+N2MAX         !from 0 upto 2*N2MAX
      DO 301 NPR=MAX(1,ABS(N-N1)),MIN(N+N1,N2MAX)   !n' of Eq. (4.18); always >= 1
*
      IAM1=0
      IAM2=MIN(N,NPR)
      IBM=0
      CALL CGORDANR(N,N1,NPR,IAM1,IAM2,IBM,CGN)
            IF(MOD((N+N1+NPR),2).EQ.0) THEN
                  AX1(N,NPR,N1+1)=CGN(0-IAM1+1)*(T(N,NPR,1)
     *                    +T(N+N1MAX,NPR+N1MAX,1))*0.5D0
                  AX2(N,NPR,N1+1)=CGN(0-IAM1+1)*(T(N,NPR,1)
     *                    -T(N+N1MAX,NPR+N1MAX,1))*0.5D0
            ELSE
                  AX1(N,NPR,N1+1)=0.0D0
                  AX2(N,NPR,N1+1)=0.0D0
            END IF
      DO 300 M1=1,MIN(N,NPR)
      if(CGN(M1-IAM1+1).ne.0) then
            IF(MOD((N+N1+NPR),2).EQ.0) THEN
                  TA(1)=T(N,NPR,M1+1)+T(N+N1MAX,NPR+N1MAX,M1+1)
                  TA(2)=T(N,NPR,M1+1)-T(N+N1MAX,NPR+N1MAX,M1+1)
            ELSE
                  TA(1)=T(N,NPR+N1MAX,M1+1)+T(N+N1MAX,NPR,M1+1)
                  TA(2)=T(N,NPR+N1MAX,M1+1)-T(N+N1MAX,NPR,M1+1)
            END IF
            AX1(N,NPR,N1+1)=AX1(N,NPR,N1+1)+TA(1)*CGN(M1-IAM1+1)
            AX2(N,NPR,N1+1)=AX2(N,NPR,N1+1)+TA(2)*CGN(M1-IAM1+1)
      end if
 300  CONTINUE
            CT=2.0D0*(0.0D0,1.0D0)**(NPR-N)/SQRT(db(npr+1))
            AX1(N,NPR,N1+1)=AX1(N,NPR,N1+1)*CT
            AX2(N,NPR,N1+1)=AX2(N,NPR,N1+1)*CT
*
 301  CONTINUE
 302  CONTINUE
 303  CONTINUE
*
C      Calculation of Bmnn1 (JOSA, Eqs. 4.17)
      DO 330 N=1,N2MAX
      DO 320 N1=0,N2MAX+N2MAX
      DO M=-N2MAX,N2MAX       !BX1,BX2 initialization to zero
            if(BX1(N,N1,M+N2MAX+1).ne.0) BX1(N,N1,M+N2MAX+1)=0.0D0
            if(BX2(N,N1,M+N2MAX+1).ne.0) BX2(N,N1,M+N2MAX+1)=0.0D0
            END DO
      DO 310 NPR=MAX(1,ABS(N-N1)),MIN(N+N1,N2MAX)
            IAM1=-N2MAX
            IAM2=N2MAX
            IBM=-1
*
      CALL CGORDANR(N,NPR,N1,IAM1,IAM2,IBM,CGN)
*
            H=SQRT(DB(NPR+1)/DB(N1+1))
       DO IAM=IAM1,IAM2
         CGN(IAM-IAM1+1)=CGN(IAM-IAM1+1)*H
         IF(MOD((N-IAM),2).NE.0) CGN(IAM-IAM1+1)=-CGN(IAM-IAM1+1)
       END DO
*
      DO 309 M=-N2MAX,N2MAX
        if(CGN(M-IAM1+1).ne.0) then
        BX1(N,N1,M+N2MAX+1)=BX1(N,N1,M+N2MAX+1)
     &                         +CGN(M-IAM1+1)*AX1(N,NPR,N1+1)
        BX2(N,N1,M+N2MAX+1)=BX2(N,N1,M+N2MAX+1)
     &                         +CGN(M-IAM1+1)*AX2(N,NPR,N1+1)
            end if
*
 309  CONTINUE
 310  CONTINUE
 320  CONTINUE
 330  CONTINUE
*
C      Calculation of Dmnnso (JOSA, Eqs. 4.29-4.33)
*
      DO 349 N=1,N2MAX
      DO 348 NSO=1,N2MAX
*
        DO 342 M=-MIN(N,NSO),MIN(N,NSO)
        ax1(N,NSO,M+N2MAX+1)=0.0D0
        ax2(N,NSO,M+N2MAX+1)=0.0D0
        DO 341 N1=ABS(M-1),N2MAX+N2MAX
      if(BX1(N,N1,M+N2MAX+1).ne.0.and.BX1(NSO,N1,M+N2MAX+1).ne.0)
     & then
            ax1(N,NSO,M+N2MAX+1)=ax1(N,NSO,M+N2MAX+1)+DB(N1+1)*
     *            BX1(N,N1,M+N2MAX+1)*DCONJG(BX1(NSO,N1,M+N2MAX+1))
            end if
      if(BX2(N,N1,M+N2MAX+1).ne.0.and.BX2(NSO,N1,M+N2MAX+1).ne.0)
     &      then
            ax2(N,NSO,M+N2MAX+1)=ax2(N,NSO,M+N2MAX+1)+DB(N1+1)*
     *            BX2(N,N1,M+N2MAX+1)*DCONJG(BX2(NSO,N1,M+N2MAX+1))
            end if
*
 341  CONTINUE
 342  CONTINUE
*
         DO 344 M=MAX(-N,-NSO+2),MIN(N,NSO+2)
            D3(N,NSO,M+N2MAX+1)=0.0D0
            D4(N,NSO,M+N2MAX+1)=0.0D0
            D5(N,NSO,M+N2MAX+1)=0.0D0
         DO 343 N1=ABS(M-1),N2MAX+N2MAX
*
      if(BX1(N,N1,M+N2MAX+1).ne.0.and.BX1(NSO,N1,2-M+N2MAX+1).ne.0)
     &      then
            D3(N,NSO,M+N2MAX+1)=D3(N,NSO,M+N2MAX+1)+DB(N1+1)*
     *               BX1(N,N1,M+N2MAX+1)*DCONJG(BX1(NSO,N1,2-M+N2MAX+1))
            end if
            if(BX2(N,N1,M+N2MAX+1).ne.0) then
            if(BX2(NSO,N1,2-M+N2MAX+1).ne.0) then
      D4(N,NSO,M+N2MAX+1)=D4(N,NSO,M+N2MAX+1)+DB(N1+1)*
     *               BX2(N,N1,M+N2MAX+1)*DCONJG(BX2(NSO,N1,2-M+N2MAX+1))
      if(BX1(NSO,N1,2-M+N2MAX+1).ne.0) then
      D5(N,NSO,M+N2MAX+1)=D5(N,NSO,M+N2MAX+1)+DB(N1+1)*
     *         BX2(N,N1,M+N2MAX+1)*DCONJG(BX1(NSO,N1,2-M+N2MAX+1))
      end if
      end if
      end if
 343  CONTINUE
 344  CONTINUE
 348  CONTINUE
 349  CONTINUE

C      Calculation of the G1-G5 (JOSA Eqs. 4.23-4.27)
C      and As,Bs (Id. Eqs. 2.31-2.36)
*
      DO S=0,N2MAX+N2MAX
            G1(S+1)=0.0D0
            G2(S+1)=0.0D0
            G3(S+1)=0.0D0
            G4(S+1)=0.0D0
            G5(S+1)=0.0D0
            AP1(S+1)=0.0D0
            AP2(S+1)=0.0D0
            AP3(S+1)=0.0D0
            AP4(S+1)=0.0D0
            BP1(S+1)=0.0D0
            BP2(S+1)=0.0D0
      END DO
*
      DO 410 S=0,N2MAX+N2MAX
      DO 400 N=1,N2MAX
*
      IF(ABS(N-S).GT.N2MAX) GO TO 400
      DO 390 NSO=MAX(1,ABS(N-S)),MIN(N+S,N2MAX)
      Z1=0.0D0
      Z2=0.0D0
      Z3=0.0D0
      Z4=0.0D0
      Z5=0.0D0
c      To save memory, D1 and D2 are stored at AX1,AX2
      IAM1=0
      IAM2=MIN(N,NSO)
      IBM=0
*
      CALL CGORDANR(N,S,NSO,IAM1,IAM2,IBM,CGN)
*
      DO 360 M=-MIN(N,NSO),MIN(N,NSO)
*
      IF(M.LT.0) THEN
            IF(CGN(-M-IAM1+1).NE.0) THEN
            IF(MOD((N+S-NSO),2).EQ.0) THEN
                  Z1=Z1+CGN(-M-IAM1+1)*ax1(N,NSO,M+N2MAX+1)
                  Z2=Z2+CGN(-M-IAM1+1)*ax2(N,NSO,M+N2MAX+1)
            ELSE
                  Z1=Z1-CGN(-M-IAM1+1)*ax1(N,NSO,M+N2MAX+1)
                  Z2=Z2-CGN(-M-IAM1+1)*ax2(N,NSO,M+N2MAX+1)
            END IF
            END IF
      ELSE               !M.GE.0
            IF(CGN(M-IAM1+1).NE.0) THEN   !CGN1 may be used before set
            IF(CGN1.EQ.0.AND.M.EQ.1) CGN1=CGN(M-IAM1+1)       !=CGN(2)
            Z1=Z1+CGN(M-IAM1+1)*ax1(N,NSO,M+N2MAX+1)
            Z2=Z2+CGN(M-IAM1+1)*ax2(N,NSO,M+N2MAX+1)
      END IF
      END IF
*
360      CONTINUE
*
C      D4 and D5 are stored as D1 and D2, to save memory
      IF(S.GE.2) THEN
      IAM1=-MIN(N,NSO+2)
      IAM2=MIN(N,NSO-2)
      IBM=2
*
      CALL CGORDANR(N,S,NSO,IAM1,IAM2,IBM,CGN)
*
      DO 380 M=MAX(-N,-NSO+2),MIN(N,NSO+2)
*
            IF(CGN(-M-IAM1+1).NE.0) THEN    !CGN2 may be used before set
                  IF(CGN2.EQ.0.AND.M.EQ.1) CGN2=CGN(-M-IAM1+1)  !=CGN(0)
                  Z3=Z3+CGN(-M-IAM1+1)*D3(N,NSO,M+N2MAX+1)
                  Z4=Z4+CGN(-M-IAM1+1)*D4(N,NSO,M+N2MAX+1)
                  Z5=Z5+CGN(-M-IAM1+1)*D5(N,NSO,M+N2MAX+1)
            END IF
380      CONTINUE
*
      END IF
      H=SQRT(DB(N+1)/DB(NSO+1))*DB(S+1)/(CSCA+CSCA)
      IF(MOD((N+NSO-S),2).EQ.0) THEN
            G2(S+1)=G2(S+1)+H*Z2*CGN1
            IF(S.GE.2) G4(S+1)=G4(S+1)+H*Z4*CGN2
      ELSE
            G2(S+1)=G2(S+1)-H*Z2*CGN1
            IF(S.GE.2) G4(S+1)=G4(S+1)-H*Z4*CGN2
      END IF
      G1(S+1)=G1(S+1)+H*Z1*CGN1
      IF(S.GE.2) G3(S+1)=G3(S+1)+H*Z3*CGN2
      IF(S.GE.2) G5(S+1)=G5(S+1)-H*Z5*CGN1
      CGN1=0.0D0
      CGN2=0.0D0
*
390      CONTINUE
400      CONTINUE
*
      AP1(S+1)=DREAL(G1(S+1)+G2(S+1))
      AP2(S+1)=DREAL(G3(S+1)+G4(S+1))
      AP3(S+1)=DREAL(G3(S+1)-G4(S+1))
      AP4(S+1)=DREAL(G1(S+1)-G2(S+1))
      BP1(S+1)=2.0D0*DREAL(G5(S+1))
      BP2(S+1)=2.0D0*IMAG(G5(S+1))
      XM=ABS(AP1(S+1))

      IF(ABS(AP2(S+1)).GT.XM) XM=ABS(AP2(S+1))
      IF(ABS(AP3(S+1)).GT.XM) XM=ABS(AP3(S+1))
      IF(ABS(AP4(S+1)).GT.XM) XM=ABS(AP4(S+1))
      IF(ABS(BP1(S+1)).GT.XM) XM=ABS(BP1(S+1))
      IF(ABS(BP2(S+1)).GT.XM) XM=ABS(BP2(S+1))

      XM=XM*1000.0D0
C====================================================================
C                 MUELLER MATRIX CALCULATION
C====================================================================

      IF((DELTA.EQ.0).AND.(XM.LT.1.0D-3)) THEN
            STOPE=S
            GO TO 420
      ELSE IF((DELTA.NE.0).AND.(XM.LT.DELTA)) THEN
            STOPE=S
          GO TO 420
      END IF

410      CONTINUE
            STOPE=N2MAX+N2MAX
420      CONTINUE
*
C      Calculation of angle-dependent coefficients for a given angle
C      \Theta=ANG(N)
*
      DO N=1,NTOPE
*
      SEN=SIN(ANG(N))
      COSS=COS(ANG(N))
      PX1(1,N)=COSS
      PX1(2,N)=(3.0D0*COSS*COSS-1.0D0)/2.0D0
      PX2(1,N)=0.0D0
      PX2(2,N)=-SEN*SEN*DSQRT(0.375D0)
      PX3(1,N)=0.0D0
      PX3(2,N)=(1.0D0+COSS)*(1.0D0+COSS)/4.0D0
      PX4(1,N)=0.0D0
      PX4(2,N)=(1.0D0-COSS)*(1.0D0-COSS)/4.0D0
      A1(N)=AP1(1)+AP1(2)*PX1(1,N)+AP1(3)*PX1(2,N)
      A4(N)=AP4(1)+AP4(2)*PX1(1,N)+AP4(3)*PX1(2,N)
      A2(N)=(AP2(3)*(PX3(2,N)+PX4(2,N))
     &                +AP3(3)*(PX3(2,N)-PX4(2,N)))/2.0D0
      A3(N)=(AP2(3)*(PX3(2,N)-PX4(2,N))
     &                +AP3(3)*(PX3(2,N)+PX4(2,N)))/2.0D0
      B1(N)=BP1(3)*PX2(2,N)
      B2(N)=BP2(3)*PX2(2,N)
*
      DO S=3,STOPE
        PX1(S,N)=(2.0D0*S-1.0D0)*COSS*PX1(S-1,N)-(S-1.0D0)*PX1(S-2,N)
        PX1(S,N)=PX1(S,N)/S
        PX2(S,N)=(2.0D0*S-1.0D0)*COSS*PX2(S-1,N)
        PX2(S,N)=PX2(S,N)-DSQRT((S-1.0D0)*(S-1.0D0)-4.0D0)*PX2(S-2,N)
        PX2(S,N)=PX2(S,N)/DSQRT(S*S-4.0D0)
        PX3(S,N)=(2.0D0*S-1.0D0)*(S*(S-1.0D0)*COSS-4.0D0)*PX3(S-1,N)
        PX3(S,N)=PX3(S,N)-S*((S-1.0D0)*(S-1.0D0)-4.0D0)*PX3(S-2,N)
        PX3(S,N)=PX3(S,N)/((S-1.0D0)*(S*S-4.0D0))
        PX4(S,N)=(2.0D0*S-1.0D0)*(S*(S-1.0D0)*COSS+4.0D0)*PX4(S-1,N)
        PX4(S,N)=PX4(S,N)-S*((S-1.0D0)*(S-1.0D0)-4.0D0)*PX4(S-2,N)
        PX4(S,N)=PX4(S,N)/((S-1.0D0)*(S*S-4.0D0))
        A1(N)=A1(N)+AP1(S+1)*PX1(S,N)
        A2(N)=A2(N)+(AP2(S+1)*(PX3(S,N)+PX4(S,N))+AP3(S+1)*(PX3(S,N)-
     *               PX4(S,N)))/2.0D0
        A3(N)=A3(N)+(AP2(S+1)*(PX3(S,N)-PX4(S,N))+AP3(S+1)*(PX3(S,N)+
     *               PX4(S,N)))/2.0D0
        A4(N)=A4(N)+AP4(S+1)*PX1(S,N)
        B1(N)=B1(N)+BP1(S+1)*PX2(S,N)
        B2(N)=B2(N)+BP2(S+1)*PX2(S,N)
      END DO        !S=3,STOPE
*
C      Now for a few inequalities that have to be fulfilled
C      If any of these fail (to withing an accuracy of DELTA/100),
C      Then there's something wrong the calculation process
C      (Most likely, convergence failures)
C      CFA(1)..CFA(5) will tell you which inequality fails,
C      and for how many angles

      FALLOS(1)=(A1(N)+A2(N))*(A1(N)+A2(N))
      FALLOS(1)=FALLOS(1)-(A3(N)+A4(N))*(A3(N)+A4(N))
      FALLOS(1)=FALLOS(1)-4.0D0*(B1(N)*B1(N)+B2(N)*B2(N))
      FALLOS(1)=FALLOS(1)/(A1(N)*A1(N))
      FALLOS(2)=1.0D0-(A2(N)/A1(N))
      FALLOS(2)=FALLOS(2)-ABS((A3(N)/A1(N))-(A4(N)/A1(N)))
      FALLOS(3)=1.0D0-(B1(N)/A1(N))
      FALLOS(3)=FALLOS(3)-ABS((A2(N)/A1(N))-(B1(N)/A1(N)))
      FALLOS(4)=1.0D0+(B1(N)/A1(N))
      FALLOS(4)=FALLOS(4)-ABS((A2(N)/A1(N))+(B1(N)/A1(N)))
*
C      DPRIMA is the accuracy parameter allowed to the inequalities
C      If you want to change it, here's your chance
*
      DPRIMA=DELTA/100
      IF(DELTA.EQ.0) DPRIMA=1.0D-6
      IF(FALLOS(1).LT.(-DPRIMA)) CFA(1)=CFA(1)+1 !CFA may be used before set
      IF(FALLOS(2).LT.(-DPRIMA)) CFA(2)=CFA(2)+1
      IF(FALLOS(3).LT.(-DPRIMA)) CFA(3)=CFA(3)+1
      IF(FALLOS(4).LT.(-DPRIMA)) CFA(4)=CFA(4)+1
*
      END DO   !N=1,NTOPE
*
      CF=CFA(1)+CFA(2)+CFA(3)+CFA(4)+CFA(5)
*
C      If the code has come this far, congratulations, you've made it
C      The output parameters are A1(n)..A4(n),B1(n),B2(n)
C      That means that, for an angle ANG(n), the Mueller matrix is:
C
C            (A1 B1  0   0)                          (F11 F12  0   0 )
C            (B1 A2  0   0)      or, in different      (F12 F22  0   0 )
C            (0  0   A3 B2)      notation ....            ( 0   0  F33 F34)
C            (0  0  -B2 A4)                        ( 0   0 -F34 F44)
C
C      This is the case for axisymmetric, randomly-oriented particles
C      See e.g. Bohren and Huffman's book
C      Of course, if you chose NOT to calculate any Mueller element...
*
C                       END OF INTENSITY CALCULATION
*======================================================================
*
 444      continue
*
C
C      And here are the final results
C      First, the efficiencies
      WRITE (*,*)
      WRITE (*,*) 'LAMBDA=', lambda
      WRITE (*,*) 'The results are:'
      WRITE (*,*) 'Qext=',Qext
      WRITE (*,*) 'Qsca=',Qsca
      WRITE (*,*) 'Qabs=',qabs
*
      if (de.ne.4) then
*
      if (((nocon1.gt.10000000).or.(nocon2.gt.10000000))
     & .and.(NMX**2.lt.10000)) then
      WRITE (6,*) 'Convergence not achieved'
       else
      WRITE (6,*) 'Convergence achieved with:'
       end if
*
      WRITE (*,*) 'nocon1=',nocon1,', nocon2=',nocon2
      else              !de.eq.4
      WRITE (*,*) 'convergence check not performed'
      end if
*
      WRITE (1,*) qext,qsca,qabs
      WRITE (nout,*) lambda,qsca
      WRITE (nout+1,*) lambda,qext
      WRITE (nout+2,*) lambda,qabs
      WRITE (nout+5,*) lambda, qsca/qext, qext-(qsca +qabs)
      if (qsca/qext.gt.1)
     & write(6,*)'WARNING - albedo greater than 1 for lambda=',lambda

C      Then, the warning in case the inequalities have failed.
      IF(CF.NE.0) then
      WRITE (*,*)
      WRITE (*,*)'W A R N I N G'
        WRITE (*,*)'Some checking inequalities have failed'
      WRITE (*,*)'Please see the efficiencies.dat file for details'
      write (*,*)
      write (2,*) cfa(1),fallos(1)
      write (2,*) cfa(2),fallos(2)
      write (2,*) cfa(3),fallos(3)
      write (2,*) cfa(4),fallos(4)
       END IF
      CONTINUE
cx      WRITE (*,*) 'Want Mueller matrix data stored in file?
cx     & (1=yes, 0=no)'
cx 1010      READ (*,*) NCERO

      NCERO=1

cx      IF(NCERO.NE.1.AND.NCERO.NE.0) GO TO 1010
      IF (NCERO.EQ.1) THEN
      DO S=1,NTOPE
            ANG(S)=ANG(S)/PG
            WRITE (3,*) S,REAL(ANG(S)),A1(S),A2(S)
            WRITE (4,*) S,REAL(ANG(S)),A3(S),A4(S)
            WRITE (11,*) S,REAL(ANG(S)),B1(S),B2(S)
      END DO
      WRITE (*,*) 'Mueller matrix elements stored on files
     & mm1.dat, mm2.dat and and mm3.dat'
      END IF
cx      WRITE (*,*)
cx      WRITE (*,*) 'End of the program'

      RETURN
      END
C

      SUBROUTINE PUJOL(CA,MR,ZKE,EPS,DELTA,CEXT,CSCA,T,NOCON,
     *  N1MAX,N2MAX)
C--------/---------/---------/---------/---------/---------/---------/--
C      >>> CA,MR,ZKE,EPS,DELTA
C      <<< CEXT,CSCA,T,NOCON,N1MAX,N2MAX
C
C      SUBROUTINE FOR CONVERGENCE IN N1MAX,N2MAX,NGAUSS
C      Beginnig with an ad-hoc value, the routine sequntially
C      increases the value of the respective angular-momentum cutoffs
C      N1MAX,N2MAX for the core and exterior shell, till desired
C      converged of cross sections is reached.
C
C      NDGS ... the step size for Gauss integration
C               The number of integration points is set as NMAX*NDGS
C--------/---------/---------/---------/---------/---------/---------/--
      implicit none
      INTEGER SU,SUP,ST,ST1,ST2
      PARAMETER (SU=50,SUP=0,ST=SU+SUP,ST1=ST+1,ST2=ST+ST)
      REAL*8 CBX,CEXT,CSCA,K1,K2,K3,K4,K5,K6,DELTA,C1(2),C2(2),EPS
      COMPLEX*16 ZKE,MR,CCX,T(ST2,ST2,ST1),T1(ST2,ST2,ST1)
      INTEGER*4 ADHOC,I,J,NOCON
      INTEGER NMAX,N1MAX,N2MAX,NGAUSS,ML,CA,N1MAX2,NDGS
      EXTERNAL TNATURAL
      COMMON/CADOS/N1MAX2     !enters as N1MAX for CA=2
      NOCON=0

C       ========== Convergence for N1max ==========
C      We need a starting value for NMAX
C      The first logical step is Wiscombe's criterion x+4x^0.3333+2
C      However, this usually means a great degree of accuracy
C      So, it's better to lower a bit, by using ADHOC
C      Try your own values.  You can try to raise it, if needed.

      ADHOC=-10
      if(ca.eq.2) ADHOC=-15
      NMAX=INT(ZKE+4.0D0*(ZKE**(1.0D0/3.0D0))+2)  !Wiscombe's criterion
      NMAX=NMAX+ADHOC

C      Of course, we do want NMAX to be positive, hmm?
      IF(NMAX.LT.1) NMAX=1

C      And we also want to have NMAX(core)<NMAX(shell)
      IF(CA.EQ.2.AND.NMAX.LT.N1MAX2) NMAX=N1MAX2

C      NDGS=4 is a good value, but you might want to raise it for
C      highly nonspherical particles.

      NDGS=4
      NGAUSS=NMAX*NDGS

C      For convergence checking, we only need m=0
      ML=0
*
      CALL TNATURAL(CA,MR,ZKE,EPS,NMAX,ML,NGAUSS,CEXT,CSCA,T,T1)
*
      C1(2)=CEXT
      C2(2)=CSCA
      K3=1000.0D0
      K4=1000.0D0
      K1=100.0D0
      K2=100.0D0

30      NMAX=NMAX+1

C      If NMAX is larger than we can afford, the we must stop
      IF(NMAX.GT.ST) THEN
            WRITE (*,*) 'NMAX IS HIGHER THAT THE MAXIMUM POSSIBLE
     &                 VALUE OF',SU
                pause
            STOP 777
      END IF

      NGAUSS=NMAX*NDGS   !with the new value of NMAX increased by one
*
      CALL TNATURAL(CA,MR,ZKE,EPS,NMAX,ML,NGAUSS,CEXT,CSCA,T,T1)
*
      C1(1)=C1(2)
      C2(1)=C2(2)
      C1(2)=CEXT
      C2(2)=CSCA
      K1=ABS(1-C1(1)/C1(2))
      K2=ABS(1-C2(1)/C2(2))

C      K1, K2 are the differences between Cext, Csca for Nmax and for Nmax-1
C      (To be precise, Cext and Csca are only calculated using m=0)
C      If K1,K2<0.1*Delta, then our convergence goal has been achieved
C      Otherwise, we make NMAX+1, and continue
*
      IF(MAX(K1,K2).GT.DELTA/10) GO TO 30
      CONTINUE
      N1MAX=NMAX
*
C       ========== Convergence for N2max ==========
C      We have calculated the T matrix for m=0 to within the desired accuracy
C      However, many times we do not need too use all their elements
C      So we will get N2MAX, and will not use T-matrix elements with a larger index
C      If CA=1 (inner surface), it's better to use the complete matrix: N2max=N1max

      IF(CA.EQ.1) THEN
            N2MAX=N1MAX
            GO TO 60
      END IF

C      For the outer surface, let's see how many T-matrix elements are needed
C      We will calculate Cext, Csca for N2max and for N1max
C      When those sections converge to within 0.1*Delta, we will stop
C      A good starting point will be half N1max

ct      NMX=N1MAX
      N2MAX=MAX(1,N1MAX/2)    !N1MAX/2        !

C      Let's now calculate Cext,Csca
C      We will again suppose m=0
      CBX=0.0D0
      CCX=(0.0D0,0.0D0)
      DO I=1,N2MAX
        DO J=1,N2MAX
          CBX=CBX+ABS(T(I,J,1))*ABS(T(I,J,1))
          CBX=CBX+ABS(T(I+NMAX,J+NMAX,1))*ABS(T(I+NMAX,J+NMAX,1))
         END DO
         CCX=CCX+T(I,I,1)+T(I+NMAX,I+NMAX,1)
      END DO
      C1(2)=-DREAL(CCX)
      C2(2)=CBX
50      N2MAX=N2MAX+1
C      If N2MAX=N1MAX, it means we need all T-matrix elements; well...
      IF(N2MAX.EQ.N1MAX) GO TO 60
      CBX=0.0D0
      CCX=(0.0D0,0.0D0)
      DO I=1,N2MAX
        DO J=1,N2MAX
          CBX=CBX+ABS(T(I,J,1))*ABS(T(I,J,1))
          CBX=CBX+ABS(T(I+NMAX,J+NMAX,1))*ABS(T(I+NMAX,J+NMAX,1))
         END DO
         CCX=CCX+T(I,I,1)+T(I+NMAX,I+NMAX,1)
      END DO
      C1(2)=-DREAL(CCX)
      C2(2)=CBX
      K1=ABS(1-C1(2)/CEXT)
      K2=ABS(1-C2(2)/CSCA)
C      Just like before, if we get convergente better than 0.1*Delta, success!
      IF(MAX(K1,K2).GT.DELTA/10) GO TO 50
60      CONTINUE

C       ========== Convergence in NGAUSS ==========
C      And now, let's see how many Gauss quadrature points we need
C      We will start with NDGS*NMAX

      NMAX=N1MAX
      ML=0
      NGAUSS=NMAX*NDGS

      CALL TNATURAL(CA,MR,ZKE,EPS,NMAX,ML,NGAUSS,CEXT,CSCA,T,T1)

      K3=1000.0D0
      K4=1000.0D0
      K1=100.0D0
      K2=100.0D0
70      C1(1)=CEXT
      C2(1)=CSCA

C      And now, let's increase NGAUSS by NDGS
      NGAUSS=NGAUSS+NDGS
*
      CALL TNATURAL(CA,MR,ZKE,EPS,NMAX,ML,NGAUSS,CEXT,CSCA,T,T1)
*
C      And now for a little difference
C      Here I demand that the differences in Cext,Csca decreases
C      for two consecutive values of Ngauss:
C      That is, Cext(ngauss)<Cext(ngauss-1)<Cext(ngauss-2)
C      The reason is, K1 and K2 can increase before going down again
C      So I let them go up on the condition that they go down again
C      If they don't, it might mean the convergence cannot be achieved
C      And further increases in Ngauss might only worsen results
C      In the past, I did the same check in N1MAX, but since it cannot
C      go higher than SU, it would have to stop anyway.

      K5=K3
      K6=K4
      K3=K1
      K4=K2
      K1=ABS(1-C1(1)/CEXT)
      K2=ABS(1-C2(1)/CSCA)
      IF((MAX(K1,K2).GT.MAX(K3,K4)).AND.(MAX(K3,K4).GT.MAX(K5,K6)))
     & THEN
        NOCON=NOCON+10000000
        NGAUSS=NGAUSS-NDGS-NDGS
        GO TO 80
      END IF
      IF(MAX(K1,K2).GT.DELTA/10) GO TO 70
 80   CONTINUE
*
      CALL TNATURAL(CA,MR,ZKE,EPS,N1MAX,N2MAX,NGAUSS,CEXT,CSCA,T,T1)
*
      NOCON=INT(NGAUSS+1000*N2MAX+100000*N1MAX)+NOCON

C      Please take a look at NOCON.  It's a reminder of the T-matrix size
C      It has the form abbccddd, where bb=N1max, cc=N2max, ddd=Ngauss
C      a=1 if there is no convergence in Ngauss, 0 otherwise

      IF(CA.EQ.1) N1MAX2=N1MAX

      RETURN
      END
C

      SUBROUTINE TNATURAL(CA,MR,ZKE,EPS,NMAX,ML,NGAUSS,CEXT,CSCA,T,T1)
C--------/---------/---------/---------/---------/---------/---------/--
C      >>> CA,MR,ZKE,EPS,NMAX,ML,NGAUSS
C      <<< CEXT,CSCA,T
C=======
C      CA=1 ... TNATURAL runs for the inner surface
C                  Returns TMAT for a core assuming surrounding medium
C                   having the same refractive index as coating
C      CA=2 ... TNATURAL runs for the outer surface: returns TMAT
C                   for a particle with no core
C      MR(=MR1/MR2) ... refractive index ratio, where
C            mr1= Absolute RI of the core (complex)
C            mr2= Absolute RI of the shell (complex)
C      ZKE ... corresponding subshell size parameter in shell medium
C                     K*Rint*n_2 (CA=1) or K*R_ext*n_h (CA=2)
C             where K is the wave vector in the vacuum
C
C      EPS ... excentricity
C      NMAX ... angular-momentum cutoff (not changed on the output!)
C      ML ...   azimuthal parameter (not changed on the output!)
C      NGAUSS ... Number of abscissas and weights in Gaussian integration
C     NCHECK ... to distinguish the cases of plane-symmetric particles,
C                  NCHECK=0 ... particles are not +/- theta symmetric
C                  NCHECK=1 ... particles have +/- theta symmetry
C      N1MAX1=NMX1 ... transferred in via COMMON; used only for CA=2
C=======
C      SUBROUTINE FOR CALCULATION OF THE NATURAL T-MATRIX
C      (Natural=Z axis along the symmetry axis)
C
C      COMPLEX*16 T(ST2,ST2,ST1)
C
C     T_n= [Q_n^{11} - Q_n^{13} T_{n-1}] x
C                 [Q_n^{31} - Q_n^{33} T_{n-1}]^{-1}
C
C      This is the core of the apple, so to speak
C      Made along the lines of Barber and Hill: "Light Scattering
C      by Particles: Computational Methods" (Barber-Hill from now on)
C      With some changes e.g. the use of Wigner instead of Legendre functions
C      So don't be surprised if the equations slightly differ from that book's
C      In case of doubt, read my lips: it works!
C
C      COMPLEX*16 HA(ST1),BC(ST1),BC2(ST1) ... h,j,y spherical Bessel functions
C      Rx,Rw: abscissas and weights
C--------/---------/---------/---------/---------/---------/---------/--
      implicit none
      INTEGER NGAUD,SU,SUP,ST,ST1,ST2    !,ST3
      PARAMETER (SU=50,SUP=0,ST=SU+SUP,ST1=ST+1,ST2=ST+ST) !,ST3=ST2+1
      PARAMETER(NGAUD=300)   !internal cut-off of the number of Gauss abscissas

      INTEGER REL(ST1),IPIV(ST2)                  !,REL2(ST1)
      INTEGER M,I,J,K,IREL,JREL,TBUCLE,N1MAX1,NQ
      INTEGER N,NMIN,M1,NMAX,ML,NGAUSS,NCHECK,CA,NN,NK
      DOUBLE PRECISION SENTH,COSTH,SENT2,THETA,WTSEN
      DOUBLE PRECISION EPS,RX(NGAUD),RW(NGAUD),WIG(ST1)
      DOUBLE PRECISION GM(ST1)   !PI,
      DOUBLE PRECISION DEL(ST),D(ST1),DP(ST1),GM1
      DOUBLE PRECISION KAYUDA,CEXT,CSCA,DELT,EMACH
      COMPLEX*16 ZX1,ZAX,ZDX,ZIA(ST1),ZIZ(ST1),ZIB(ST1),
     & ZJA(ST1),ZJZ(ST1),ZKA(ST1),CB(ST1)
      COMPLEX*16 A(ST2,ST2),B(ST2,ST2),AA(ST2,ST2),BB(ST2,ST2)
      COMPLEX*16 HA(ST1),BC(ST1),BC2(ST1),ZWIG2(ST1)
      COMPLEX*16 ZDERI,ZDESEN,ZI1,ZA,ZB,ZC,ZD,ZE,ZF,ZFX,ZFY,ZPJ,ZR
      COMPLEX*16 ZGM2,ZKE,ZKR,MR,ZMD,MKR,CZERO
      COMPLEX*16 T(ST2,ST2,ST1),T1(ST2,ST2,ST1),ZX(ST2),CCX
      COMPLEX*16 JL(0:ST),NL(0:ST),DRJL(0:ST),DRNL(0:ST),
     & HL(0:ST),DRHL(0:ST)

      EXTERNAL GNZBESS,GNRCBSH,WIGNER,GAUSS     !BESSEL,BESSEL2,HANKEL,
      EXTERNAL ZGER,ZSUR,LULABY
*
      COMMON/CAUNO/N1MAX1                      !from main
*
      DATA EMACH/0.D0/                         ! /1.D-17/
      DATA CZERO/(0.D0,0.D0)/
*
C      PI=3.141592653589793238462643D0
C       CA=1: core.     CA=2: coating.
C      if m1=m2, we have a homgeneous scatterer, then T1=0 and over for CA=1
C      Same if the core has zero size

      if(ca.eq.1.and.(mr.eq.1.or.zke.eq.0)) then
            do m=0,ml
              NMIN=MAX(1,M)
              M1=M+1
              N=NMAX-NMIN+1
              DO I=1,N+N         !T1 initialization for CA=1
                DO J=1,N+N
                IF(T1(I,J,M1).NE.0) T1(I,J,M1)=CZERO
                END DO
                END DO
            end do
            cext=1.0d0
            csca=1.0d0
            go to 120
      end if

C     Gauss quadrature subroutine.  Rx,Rw: abscissas and weights
*
      CALL GAUSSA(NGAUSS,RX,RW)
*
      ZMD=MR*MR    !relative "subshell in a shell" diel. constant
C      If the particle has a plane of symmetry perpendicular to the symmetry
C      axis, calculations get simpler.  Half of the T-matrix elements need not be
C      computed (they equal zero), and the other half need surface integration
C      only from 0 to pi/2
C      For the cases of plane-symmetric particles, NCHECK=1, otherwise NCHECK=0
      NCHECK=1

C      Now, let's calculate n and n*(n+1) in hight-precision
      DO I=1,NMAX+1
            D(I)=DFLOAT(I)
            DP(I)=DFLOAT(I)*DFLOAT(I+1)
            if(I.NE.0) GM(I)=SQRT((2.0D0*D(I)+1.0D0)/DP(I))
      END DO

* Initialization
*
      DO J=1,2*NMAX
      DO I=1,2*NMAX
       A(I,J)=CZERO
       B(I,J)=CZERO
      if (ca.eq.2) AA(I,J)=CZERO
      if (ca.eq.2) BB(I,J)=CZERO
      ENDDO
      ENDDO
*
*======================================================================
C      Loop for M (azimuthal parameter)
C      That means, the T-matrix decomposes into m T-submatrices
C      That happens for axisymmetric particles.
*
      DO 40 M=0,ML            !MAIN LOOP
*
      NMIN=MAX(1,M)
      M1=M+1
      N=NMAX-NMIN+1
      DO I=NMIN,NMAX
            REL(I)=I-NMIN+1    !from 1 up to NMAX-MAX(1,M)+1
      END DO
*
C**********************************************************************
C               PRINCIPAL INTEGRATION AND TMAT ASSIGNING LOOP
C
C       Loop for TBUCLE, (angle)
      DO 30 TBUCLE=1,NGAUSS
C      As said before, if NCHECK=1 we need not integrate from pi/2 to pi
C      In that case, the value of the A,B matrices is halved
C      But we will calculate B*A^-1, so there's really no need to
C      multiply by 2 for each A,B, matrix element
C      Maybe there's a more elegant means to halve Ngauss prior to the
C      GAUSS subroutine, but this will also do
            IF(NCHECK.EQ.1.AND.TBUCLE.GT.(NGAUSS/2)) GO TO 30
            THETA=RX(TBUCLE)
            SENTH=SIN(THETA)
            SENT2=SENTH*SENTH          ! |\sin\theta|**2
            COSTH=COS(THETA)
            WTSEN=RW(TBUCLE)*SENTH     !gaussweight*\sin\theta
*
            CALL WIGNER (SENTH,COSTH,M,NMAX,WIG)
C--------/---------/---------/---------/---------/---------/---------/--
C      >>> SENTH,COSTH,M,NMA
C      <<< WIG
C
C      Subroutine to calculate Wigner functions
C      They don't explode for hight n,m values like Legendre functions
C      WIG(N+1)=dm0_n
C       Mishchenko also calculates Wigner functions - see routine VIGAMPL
C--------/---------/---------/---------/---------/---------/---------/--
C      Now we have to specify the particle surface: r(theta) and derivative
C      This is easy, as long as the particle is axisymmetric
C      I've used a spheroid, and of course you can put your own shape
C      You can even specify different surfaces for core and coating
C      e.g. an inner spheroid inside a cylinder.
C      All you have to add is an "IF CA=1 ... else ..."
C      The limit is your imagination.
C      Of course, I suggest you make sure that the core fits INSIDE the particle!
C
C     PARTICLE SHAPE:
C
C      ZKR=r(\theta)
C      ZDERI=dr(\theta)/d\theta
C
C      Spheroidal particle with axes a,b (a=axis of revolution)
C      ZKE=Equivalent-volume size parameter = (a*b*b)^0.33333, eps=b/a

      NCHECK=1

      ZKR=ZKE*EPS**(1.0D0/3.0D0)                 !cf. A=REV*EPS**(1D0/3D0) in Mi
      ZKR=ZKR/SQRT(SENTH*SENTH+(COSTH*EPS)**2)   !cf. RR=1D0/(SS+EE*CC) followed by
                                               !    R(I)=A*A*RR in Mi where
                                               !    R(I)=r(\theta)**2=ZKR**2

      ZDERI=SENTH*COSTH*(EPS*EPS-1.0D0)*ZKR*ZKR*ZKR
      ZDERI=ZDERI/(ZKE*ZKE*EPS**(2.0D0/3.0D0))  !cf. DR(I)=RR*C*S*(EE-1D0) in Mi
                              ! >>> hence ZDERI=ZKR*RR*C*S*(EE-1D0)=dr(\theta)/d\theta
                              !Note that in Mi: DR(I)=[dr(theta)/d theta]/r(theta)
                              !ZKE=ZKR/(SQRT(RR)*EPS**(1.0D0/3.0D0))
C >>> BESSEL FUNCTION PART
C      Calculation of Hankel and Bessel function
C      HA(n+1)=h_n(kr), BC(n+1)=j_n(mkr)
*
      MKR=MR*ZKR           !argument of Bessel functions along the inteface
*
      CALL GNRCBSH(ZKR,NMAX,HL,drhl)
*
      CALL GNZBESS(MKR,NMAX,jl,drjl,nl,drnl)
*
      DO I=1,NMAX+1
            BC(I)=JL(I-1)            !BC(1)=j_0 of {AS}
          HA(I)=HL(I-1)/ZKR        !HA(1)=h_0 of {AS}
      END DO
*
      CALL GNZBESS(ZKR,NMAX,jl,drjl,nl,drnl)

C      And in the outer surface, we also need Hankel functions for mkr
      IF(CA.EQ.2) THEN
*
      CALL GNRCBSH(MKR,NMAX,HL,drhl)
               DO I=1,NMAX+1
               BC2(I)=HL(I-1)/MKR        !BC2(1)=h_0 of {AS}
               END DO
      END IF               !CA.EQ.2
C <<<

            ZDESEN=ZDERI*SENTH
C      Now for another variation
C      In Barber-Hill (pp. 90-91), Legendre polynomials are used in pairs:
C      P(i)*P(j), P(i-1)*P(j), P(i)*P(j-1), P(i-1)*P(j-1)
C      First, they were substituted for Wigner functions
C      But there's a problem: when RX nears unity, P's are also close to one
C      So, ZE (see below) was calculated as a difference of nearly-unity terms
C      Therefore increasing the chance of roundoff errors
C      I've since introduced a Delta function as:
C      Del(i)=Cos(theta)*Del(i-1)-i*Sen(theta)*Sen(theta)*Wigner(i)
C      (that's for m=0; when m<>0, its slightly different
C      That way, calculations are simplified and can be done with a bit more
C      accuracy
      IF(M.EQ.0) THEN
            DEL(1)=-SENT2                       !SENT2=[sin(theta)]**2
            DEL(2)=-(SENT2+SENT2+SENT2)*COSTH
            DO I=3,NMAX
            DEL(I)=DEL(I-1)*COSTH-D(I)*SENT2*WIG(I)
            END DO
      ELSE               !M.NE.0
            DO I=NMIN,NMAX
      DEL(I)=D(I)*COSTH*WIG(I+1)-
     &                WIG(I)*SQRT((D(I)-D(M))*(D(I)+D(M)))
            END DO
      END IF
*
C      And now, let's calculate A, B matrix elements
C      I=row, J=column
C      See that A,B =0 when (i or j are less than m) and m>1
C      So the loop goes from NMIN to NMAX, NMIN being the minimum value
C      of i,j so that we only calculate nonzero matrix elements
*
      DO I=NMIN,NMAX                 !NMIN=MAX(1,M)
            ZIZ(I)=HA(I+1)*WTSEN
            ZIA(I)=HA(I)/HA(I+1)
            ZJZ(I)=JL(I)*WTSEN           !orig. DREAL(ZIZ(I))
            ZJA(I)=JL(I-1)/JL(I)         !orig. DREAL(HA(I))/DREAL(HA(I+1))
            ZIB(I)=MR*BC(I)/BC(I+1)
            ZKA(I)=ZJA(I)-ZIA(I)
            ZWIG2(I)=WIG(I+1)*DP(I)*ZDESEN            !real*8 WIG2 is complex*16 now
            IF(CA.EQ.2) CB(I)=MR*BC2(I)/BC2(I+1)
      END DO

      DO 20 I=NMIN,NMAX              !NMIN=MAX(1,M)
       DO 10 J=NMIN,NMAX
*
            ZPJ=D(I)*D(J)/ZKR                     !real*8 IPJ is complex*16 now
            ZI1=WIG(I+1)*WIG(J+1)*ZDESEN         !real*8 I1 is complex*16 now
            GM1=0.5D0*GM(I)*GM(J)/SENT2
            IF(MOD(M,2).NE.0) GM1=-GM1
            ZGM2=(0.0D0,1.0D0)*GM1
*
C      Let's now use the following sub-matrix description for A,B
C      First quadrant (I) is the (nmax*nmax) left upper side of the matrix
C      Second quadrant (II) is the (nmax*nmax) right upper side of the matrix
C      Third (III) and fourth (IV) are the same for the left and right lower side
C
C      II and III calculations.
C      For plane-symmetric particles and i+j=even, matrix elements are zero
      IF((NCHECK.EQ.1.AND.MOD((I+J),2).NE.0).OR.NCHECK.EQ.0) THEN
*
      IF(M.NE.0) THEN  !if M=0, II, III quadrant elements equal zero
*
        ZA=ZKR+ZIA(I)*(ZKR*ZIB(J)-D(J))-D(I)*ZIB(J)+ZPJ  !D(I)=dble(I)
        ZB=DEL(I)*WIG(J+1)+DEL(J)*WIG(I+1)
        ZB=ZB*ZKR
        ZC=ZI1*(DP(I)*ZIB(J)+DP(J)*ZIA(I)-ZPJ*(D(I)+D(J)+2.0D0))
        ZR=ZC+ZKA(I)*DP(J)*ZI1                         !DP(J)=dble(J*(J+1))
        ZX1=D(M)*GM1*BC(J+1)

C      ======== III quadrant calculations ============
* A,B may be used before set

      A(N+REL(I),REL(J))=A(N+REL(I),REL(J))
     &                         -ZIZ(I)*ZX1*(ZA*ZB+ZC)
       ZAX=ZA+ZKA(I)*(ZKR*ZIB(J)-D(J))
      B(N+REL(I),REL(J))=B(N+REL(I),REL(J))
     &                        -ZJZ(I)*ZX1*(ZAX*ZB+ZR)
C      ======== II quadrant calculations  ============
      ZA=ZA+ZKR*(ZMD-1.0D0)
      A(REL(I),N+REL(J))=A(REL(I),N+REL(J))
     &                       -ZIZ(I)*ZX1*(ZA*ZB+ZC)/MR
      ZAX=ZA+ZKA(I)*(ZKR*ZIB(J)-D(J))
      B(REL(I),N+REL(J))=B(REL(I),N+REL(J))
     &                      -ZJZ(I)*ZX1*(ZAX*ZB+ZR)/MR
C
C       =============================================================
C                                 CA=2
C       Additional matrices for outer surface (CA=2); first part
C        Forming    AA corresponds to Q_n^{13}
C                   BB corresponds to Q_n^{33}
C       =============================================================

      IF(CA.EQ.2) THEN

        ZA=ZKR*(1.0D0+ZIA(I)*CB(J))-D(J)*ZIA(I)-D(I)*CB(J)+ZPJ
        ZC=ZC+(CB(J)-ZIB(J))*DP(I)*ZI1
        ZR=ZC+ZKA(I)*DP(J)*ZI1
        ZX1=D(M)*GM1*BC2(J+1)

* AA,BB may be used before set

      AA(N+REL(I),REL(J))=AA(N+REL(I),REL(J))
     &                           -ZIZ(I)*ZX1*(ZA*ZB+ZC)
      ZAX=ZA+ZKA(I)*(ZKR*CB(J)-D(J))
      BB(N+REL(I),REL(J))=BB(N+REL(I),REL(J))
     &                       -ZJZ(I)*ZX1*(ZAX*ZB+ZR)
      ZA=ZKR*(ZMD+ZIA(I)*CB(J))-D(J)*ZIA(I)-D(I)*CB(J)+ZPJ
      AA(REL(I),N+REL(J))=AA(REL(I),N+REL(J))
     &                     -ZIZ(I)*ZX1*(ZA*ZB+ZC)/MR
      ZAX=ZA+ZKA(I)*(ZKR*CB(J)-D(J))
      BB(REL(I),N+REL(J))=BB(REL(I),N+REL(J))
     &                        -ZJZ(I)*ZX1*(ZAX*ZB+ZR)/MR

      END IF                   !CA.EQ.2
C
C       =============================================================
C       =============================================================
      END IF                !M.NE.0
      END IF                !NCHECK.EQ.1
C      I and IV quadrants calculations
C      For plane-symmetric particles and i+j=odd,
C      matrix elements are zero
*
      IF((NCHECK.EQ.1.AND.MOD((I+J),2).EQ.0).OR.NCHECK.NE.1) THEN
        ZD=ZKR*(ZIB(J)-ZMD*ZIA(I))+ZMD*D(I)-D(J)
        ZE=DEL(I)*DEL(J)
        IF(M.NE.0) ZE=ZE+(DP(M)-D(M))*WIG(I+1)*WIG(J+1)
        ZE=ZE*ZKR
        ZFX=ZWIG2(J)*DEL(I)
        ZFY=ZWIG2(I)*DEL(J)
        ZX1=ZGM2*BC(J+1)
C       ======== IV quadrant calculations ==============
      ZF=ZFX-ZMD*ZFY
      A(N+REL(I),N+REL(J))=A(N+REL(I),N+REL(J))
     &                        +ZIZ(I)*ZX1*(ZD*ZE+ZF)/MR
      ZDX=ZD-(ZJA(I)-ZIA(I))*ZKR*ZMD
      B(N+REL(I),N+REL(J))=B(N+REL(I),N+REL(J))
     &                        +ZJZ(I)*ZX1*(ZDX*ZE+ZF)/MR
C       ======== I quadrant calculations  ==============
      ZD=ZKR*(ZIB(J)-ZIA(I))+D(I)-D(J)
      ZF=ZFX-ZFY
      A(REL(I),REL(J))=A(REL(I),REL(J))+ZIZ(I)*ZX1*(ZD*ZE+ZF)
      ZDX=ZD-(ZJA(I)-ZIA(I))*ZKR
      B(REL(I),REL(J))=B(REL(I),REL(J))+ZJZ(I)*ZX1*(ZDX*ZE+ZF)
C
C       =============================================================
C                                     CA=2
C       Second part of forming
C            Q_n^{13} <===> AA
C            Q_n^{33} <===> BB
C       =============================================================
C
      IF(CA.EQ.2) THEN
        ZD=ZKR*(CB(J)-ZMD*ZIA(I))+ZMD*D(I)-D(J)
        ZF=ZFX-ZMD*ZFY
        ZX1=ZGM2*BC2(J+1)
      AA(N+REL(I),N+REL(J))=AA(N+REL(I),N+REL(J))
     &                       +ZIZ(I)*ZX1*(ZD*ZE+ZF)/MR
      ZDX=ZD-(ZJA(I)-ZIA(I))*ZKR*ZMD
      BB(N+REL(I),N+REL(J))=BB(N+REL(I),N+REL(J))
     &               +ZJZ(I)*ZX1*(ZDX*ZE+ZF)/MR
      ZD=ZKR*(CB(J)-ZIA(I))+D(I)-D(J)
      ZF=ZFX-ZFY
      AA(REL(I),REL(J))=AA(REL(I),REL(J))+ZIZ(I)*ZX1*(ZD*ZE+ZF)
      ZDX=ZD-(ZJA(I)-ZIA(I))*ZKR
      BB(REL(I),REL(J))=BB(REL(I),REL(J))+ZJZ(I)*ZX1*(ZDX*ZE+ZF)
      END IF                          !CA.EQ.2
C
C       =============================================================
C       =============================================================
      END IF
C      End of the J (column) loop
 10   CONTINUE
C           End of the I (row) loop
 20   CONTINUE
C      End of the TBUCLE (surface angle) loop
 30   CONTINUE
C       =============================================================
C                                     CA=2
C      Applying recurrence relation:
C
C            T_n= [Q_n^{11} - Q_n^{13} T_{n-1}] x
C                     [Q_n^{31} - Q_n^{33} T_{n-1}]^{-1}
C
C   ??? The values of T1 are only assigned for CA=1 at the end of
C       d040-loop. However, the values of T1 are only used here
C       below when CA=2 ??? Up to here, the T1 entries have been
C       merely initialized to zero, but only for CA=1???
C       =============================================================
      IF(CA.EQ.2) THEN
*
      NN=N1MAX1-NMIN+1
      NK=MIN(N,NN)
      IF(NK.EQ.0) GO TO 45
*
      DO 43 I=1,NN          !AA corresponds to Q_n^{13}
      DO 42 J=1,NN          !BB corresponds to Q_n^{33}
      DO 41 K=1,NK
*
      IF((NCHECK.EQ.1.AND.MOD((I+J),2).EQ.0).OR.NCHECK.NE.1) THEN
        A(I,J)=A(I,J)+AA(I,K)*T1(K,J,M1)+AA(I,K+N)*T1(K+NN,J,M1)
        B(I,J)=B(I,J)+BB(I,K)*T1(K,J,M1)+BB(I,K+N)*T1(K+NN,J,M1)
        A(I+N,J+N)=A(I+N,J+N)+AA(I+N,K)*T1(K,J+NN,M1)
     &                      +AA(I+N,K+N)*T1(K+NN,J+NN,M1)
        B(I+N,J+N)=B(I+N,J+N)+BB(I+N,K)*T1(K,J+NN,M1)
     &                      +BB(I+N,K+N)*T1(K+NN,J+NN,M1)
      END IF
*
      IF(M.NE.0) THEN
      IF((NCHECK.EQ.1.AND.MOD((I+J),2).NE.0).OR.NCHECK.NE.1) THEN
      A(I,J+N)=A(I,J+N)+AA(I,K)*T1(K,J+NN,M1)
     &             +AA(I,K+N)*T1(K+NN,J+NN,M1)
      B(I,J+N)=B(I,J+N)+BB(I,K)*T1(K,J+NN,M1)
     &                +BB(I,K+N)*T1(K+NN,J+NN,M1)
      A(I+N,J)=A(I+N,J)+AA(I+N,K)*T1(K,J,M1)
     &                +AA(I+N,K+N)*T1(K+NN,J,M1)
      B(I+N,J)=B(I+N,J)+BB(I+N,K)*T1(K,J,M1)
     &                  +BB(I+N,K+N)*T1(K+NN,J,M1)
      END IF
      END IF
*
41      CONTINUE
42      CONTINUE
43      CONTINUE
*
      END IF                !CA.EQ.2
*======================================================================
*======================================================================
*
45      CONTINUE
*
C      Enough with multiplications
C      Now we have to invert+multiply:  T = -B*A^(-1)
C      If CA=1, T is the T1 matrix; otherwise, if CA=2, T is the final T-matrix
C      The LULABY matrix does the product B*A^(-1), storing it in B
C      nq is the size of the A, B, matrix to invert
      nq=nmax-max(1,m)+1
      nq=nq+nq

cz      OPEN(26,FILE='amatdg.dat')
cz      rewind(26)
cz      OPEN(27,FILE='bmatdg.dat')
cz      rewind(27)
cz      do i=1,nq
cz      write (26,*)'i=', i, a(i,i)
cz      write (27,*)'i=', i, b(i,i)
cz      enddo
cz      close(26)
cz      close(27)
*
cq      CALL LULABY(A,B,NQ,ca,ml,m)
*
cq      go to 66
C  Gaussian elimination

      CALL ZGER(A,IPIV,NQ,ST2,EMACH)  !Gauss elimination of A to
                                      !a lower diagonal matrix
      DO 6 I=1,NQ
              DO K=1,NQ       !Initialization of the right-hand side B
                              !(a row vector) of the matrix equation ZX*A=B

              ZX(K)=B(I,K)
              ENDDO

      CALL ZSUR (A,IPIV,ZX,NQ,ST2,EMACH)      !Solving ZX*A=B by
                                              !backsubstitution
                                              !(ZX overwritten on exit)
             DO K=1,NQ
                B(I,K)=ZX(K)
             ENDDO
  6   CONTINUE
  66  CONTINUE

cz       OPEN(28,FILE='bmatdg1.dat')
cz      rewind(28)
cz      do i=1,nq
cz      write (28,*)'i=', i, b(i,i)
cz      enddo
cz      close(28)

      DO J=1,NMAX+NMAX
      DO I=1,NMAX+NMAX
        IF(CA.EQ.1.AND.T1(I,J,M1).NE.0) T1(I,J,M1)=CZERO
                                    !T1 initialization for CA=1
        IF(T(I,J,M1).NE.0) T(I,J,M1)=CZERO
                                    !T initialization for any CA
      END DO
      END DO

C  And now we "decompress" the T (or T1) matrix
*
*  Assign T-matrix elements = - RG(Q) * (Q**(-1)) = -B
*

      DO I=1,N              !N=NMAX-NMIN+1
      DO J=1,N
            IREL=I+NMIN-1    !From NMIN(=NMIN) up to NMAX in Mishchenko notation
            JREL=J+NMIN-1    !From NMIN(=NMIN) up to NMAX in Mishchenko notation
            T(IREL,JREL,M1)=-B(I,J)
            T(IREL+NMAX,JREL,M1)=-B(I+N,J)
            T(IREL,JREL+NMAX,M1)=-B(I,J+N)
            T(IREL+NMAX,JREL+NMAX,M1)=-B(I+N,J+N)
*
* The only place where the values of T1 are assigned.
* Note that this only happens for CA=1. The values of
* T1 are then used above in the routine when CA=2
*
            IF(CA.EQ.1) THEN
                  T1(I,J,M1)=-B(I,J)
                  T1(I+N,J,M1)=-B(I+N,J)
                  T1(I,J+N,M1)=-B(I,J+N)
                  T1(I+N,J+N,M1)=-B(I+N,J+N)
            END IF
*
            IF (A(I,J).NE.0) A(I,J)=CZERO
            IF (B(I,J).NE.0) B(I,J)=CZERO
            IF (A(I+N,J).NE.0) A(I+N,J)=CZERO
            IF (B(I+N,J).NE.0) B(I+N,J)=CZERO
            IF (A(I,J+N).NE.0) A(I,J+N)=CZERO
            IF (B(I,J+N).NE.0) B(I,J+N)=CZERO
            IF (A(I+N,J+N).NE.0) A(I+N,J+N)=CZERO
            IF (B(I+N,J+N).NE.0) B(I+N,J+N)=CZERO
*
            IF(CA.EQ.2) THEN
              IF (AA(I,J).NE.0) AA(I,J)=CZERO
              IF (BB(I,J).NE.0) BB(I,J)=CZERO
              IF (AA(I+N,J).NE.0) AA(I+N,J)=CZERO
              IF (BB(I+N,J).NE.0) BB(I+N,J)=CZERO
              IF (AA(I,J+N).NE.0) AA(I,J+N)=CZERO
              IF (BB(I,J+N).NE.0) BB(I,J+N)=CZERO
              IF (AA(I+N,J+N).NE.0) AA(I+N,J+N)=CZERO
              IF (BB(I+N,J+N).NE.0) BB(I+N,J+N)=CZERO
            END IF
*
      END DO
      END DO
*
40      CONTINUE          !Enf of the azimuth (m) loop
*
*=================================================================
C      Scattering cross-section averaged over the uniform
C     orientation distribution:
*
*      C_{sca}= \fr{2\pi}{k_1^2 } \mb{Tr}\,|T_{mnmn}^{jj} (P)|^2
*

      CSCA=0.0D0
      CEXT=0.0D0
      CCX=CZERO
      CONTINUE
*
      DO M=0,ML
cc        MTOPE=MAX(M,1)-1
*
        do i=1,nmax
        do j=1,nmax
*
            IF (M.EQ.0) THEN
                  DELT=1.0D0
            ELSE
                  DELT=2.0D0
            END IF
*
      KAYUDA=(ABS(T(I,J,M+1)))*(ABS(T(I,J,M+1)))
      KAYUDA=KAYUDA+(ABS(T(I+NMAX,J,M+1)))*(ABS(T(I+NMAX,J,M+1)))
      KAYUDA=KAYUDA+(ABS(T(I,J+NMAX,M+1)))*(ABS(T(I,J+NMAX,M+1)))
      KAYUDA=KAYUDA+
     &    (ABS(T(I+NMAX,J+NMAX,M+1)))*(ABS(T(I+NMAX,J+NMAX,M+1)))
      CSCA=CSCA+DELT*KAYUDA
      KAYUDA=0.0D0
*
      END DO
      END DO
      END DO
*=================================================================
C      Extinction cross-section averaged over the uniform
C     orientation distribution:
*
*       C_{ext}= - \fr{2\pi}{k_1^2 } \mb{Tr Re}\, T_{mnmn}^{jj}(P)
*
      ccx=CZERO
*
      DO I=1,NMAX
            CCX=CCX+(T(I,I,1)+T(I+NMAX,I+NMAX,1))
      END DO
*
      DO M=1,ML
        DO I=M,NMAX
          CCX=CCX+2.0D0*(T(I,I,M+1)+T(I+NMAX,I+NMAX,M+1))
        END DO
      END DO
*
      CONTINUE
*
      CEXT=-DREAL(CCX)
*
C      If you want Cross Sections, multiply CEXT, CSCA by 2*pi/(k*k)
C      where k = 2*pi/wavelength
C      Or, to get Qext, Qsca, multiply Cext, Csca by 2/(X*X)
C      where X is the dimensionless size parameter
C      You can do it in the main routine.  Here we're through
120      RETURN
      END
C

      SUBROUTINE WIGNER(SENTH,COSTH,M,NMAX,WIG)
C--------/---------/---------/---------/---------/---------/---------/--
C      >>> SENTH,COSTH,M,NMAX
C      <<< WIG
C======
C      Subroutine to calculate the Wigner functions
C      WIG(n+1)=dm0,n
C--------/---------/---------/---------/---------/---------/---------/--
      implicit none
      INTEGER SU,SUP,ST,ST1
      PARAMETER (SU=50,SUP=0,ST=SU+SUP,ST1=ST+1)
*
      INTEGER TI,M,NMAX
      DOUBLE PRECISION SENTH,COSTH,WIG(ST1)
C      M=0
      IF(M.EQ.0) THEN
        WIG(1)=1.0D0
        WIG(2)=COSTH
        DO TI=2,NMAX
            WIG(TI+1)=(DFLOAT(TI+TI)-1.0D0)*COSTH*WIG(TI)
            WIG(TI+1)=WIG(TI+1)-(DFLOAT(TI)-1.0D0)*WIG(TI-1)
            WIG(TI+1)=WIG(TI+1)/(DFLOAT(TI))
        END DO
        RETURN
C      M<>0
      ELSE
        WIG(M+1)=1.0D0
        IF(MOD(M,2).NE.0) WIG(M+1)=-WIG(M+1)
      DO TI=0,M-1
          WIG(TI+1)=0.0D0
          WIG(M+1)=WIG(M+1)*SENTH*DSQRT(DFLOAT(TI+1)*DFLOAT(M+TI+1))
        WIG(M+1)=WIG(M+1)/(2.0D0*DFLOAT(TI+1))
      END DO
      DO TI=M+1,NMAX
       WIG(TI+1)=DFLOAT(TI+TI-1)*COSTH*WIG(TI)
       WIG(TI+1)=WIG(TI+1)
     &      -DSQRT(DFLOAT((TI-1)*(TI-1)-M*M))*WIG(TI-1)
      WIG(TI+1)=WIG(TI+1)/DSQRT(DFLOAT(TI*TI-M*M))
      END DO
      RETURN
      END IF
      CONTINUE
      END
C

      SUBROUTINE GAUSSA(NGAUSS,RX,RW)
C--------/---------/---------/---------/---------/---------/---------/--
C      >>> NGAUSS
C      <<< RX,RW
C=======
C      Subroutine to calculate weights and abscissas in Gauss integration.
C      Reference: "Light Scattering by Particles: Computational Methods",
C      P.W. Barber y S.C. Hill, pp. 172-173
C--------/---------/---------/---------/---------/---------/---------/--
      implicit none
*
      INTEGER NGAUD
      PARAMETER(NGAUD=300)   !internal cut-off of the number of Gauss abscissas
*
      DOUBLE PRECISION RX(NGAUD),RW(NGAUD),RN1(NGAUD),RN2(NGAUD)
      DOUBLE PRECISION CN,CNN1,CON1,CON2,XI,X,PM1,PM2,RN,AUX,DER1P
      DOUBLE PRECISION DER2P,APPFCT,B,BISQ,BFROOT,RATIO,PROD
      DOUBLE PRECISION CONST,TOL,C1,C2,C3,C4,PI,P,PMP1,PP
      INTEGER K,NDIV2,NP1,NM1,NM2,IN,NGAUSS
*
      DATA PI/3.141592653589793238462643D0/
      DATA CONST/0.148678816357D0/
      DATA TOL/1.0D-15/
      DATA C1,C2/0.125D0,-0.0807291666D0/
      DATA C3,C4/0.2460286458D0,-1.824438767D0/
*
C      Warning: if TOL is lower than your computer's roundoff error,
C      your code might fall into an endless loop!
*
      IF(NGAUSS.EQ.1) THEN
            RX(1)=0.577350269189626D0
            RW(1)=1.0D0
            RETURN
      END IF
*
*
      IF(NGAUSS.GT.NGAUD) THEN
      write(6,*)'Raise NGAUD in LISAC, GAUSSA and TNATURAL to at least
     & NGAUSS=', NGAUSS
      END IF
*
      CN=DFLOAT(NGAUSS)
      NDIV2=NGAUSS/2
      NP1=NGAUSS+1
      CNN1=CN*(CN+1.0D0)
      APPFCT=1.0D0/SQRT((CN+0.5D0)**2+CONST)
      CON1=0.5D0*PI
      CON2=0.5D0*PI
      do in=2,ngauss
            rn=dfloat(in)
            rn2(in)=(rn-1.0d0)/rn     !=(in-1)/in
            rn1(in)=rn2(in)+1.0d0     !=(2*in-1)/in
      end do
*
*
*
      DO 1030 K=1,NDIV2
         B=(DFLOAT(K)-0.25D0)*PI
         BISQ=1.0D0/(B*B)          !Bisq is a first approximation to Rx
*
         BFROOT=B*(1.0D0+BISQ*(C1+BISQ*(C2+BISQ*(C3+C4*BISQ))))
         XI=COS(APPFCT*BFROOT)
1010         X=XI
         PM2=1.0D0
         PM1=X
*
            DO 1020 IN=2,NGAUSS
             RN=DFLOAT(IN)
               P=((2.0D0*RN-1.0D0)*X*PM1-(RN-1.0D0)*PM2)/RN
               PM2=PM1
               PM1=P
1020            CONTINUE
*
         PM1=PM2
         AUX=1.0D0/(1.0D0-X*X)
         DER1P=CN*(PM1-X*P)*AUX
         DER2P=((X+X)*DER1P-CNN1*P)*AUX
         RATIO=P/DER1P
         XI=X-RATIO*(1.0D0+RATIO*DER2P/(2.0D0*DER1P))
         IF(ABS(XI-X).GT.TOL) GO TO 1010
*
C      Now a bit of modification here to calculate Pn with higher accuracy
C      The above calculates PM1=Pn-1(X), not Pn-1(XI)
C      If X and XI differ, it might enlarge roundoff errors
C      So just to play in the safe side...
         X=XI
         PM2=1.0D0
         PM1=X
         PMP1=1.0D0
            DO IN=2,NGAUSS
                  RN=DFLOAT(IN)
                  P=RN1(IN)*X*PM1-RN2(IN)*PM2
                  PM2=PM1
                  PM1=P
                  PP=X*PMP1+RN*PM2
                  PMP1=PP
            END DO
*
C      End of the fine-tuning modification
C      The weights are usually calculated as wi=2*(1-xi^2)/(n*Pn-1)^2
C      However, for small xi values, the relative error of the weight is higher
C      Another way to calculate it is wi=2/[(1-xi^2)(Pn')^2]
C      This works better at the endpoints
*
         RX(K)=-XI
         RW(K)=2.0d0/((1.0D0-XI*XI)*PMP1*PMP1)
         RX(NP1-K)=-RX(K)
         RW(NP1-K)=RW(K)
*
*
*
 1030 CONTINUE          !K-loop
*
      IF(MOD(NGAUSS,2).NE.0) THEN
      RX(NDIV2+1)=0.0D0
      NM1=NGAUSS-1
      NM2=NGAUSS-2
      PROD=CN
      DO 1040 K=1,NM2,2
             PROD=PROD*DFLOAT(NM1-K)/DFLOAT(NGAUSS-K)
 1040 CONTINUE
      RW(NDIV2+1)=2.0D0/(PROD*PROD)
      END IF
*
      DO 1050 K=1,NGAUSS
            RX(K)=CON1*RX(K)+CON2
            RW(K)=CON1*RW(K)
 1050 CONTINUE
*
      RETURN
      END
C
*
      SUBROUTINE CGORDANR(IAJ,IBJ,ICJ,IAM1,IAM2,IBM,CGN)
C--------/---------/---------/---------/---------/---------/---------/--
C      >>> IAJ,IBJ,ICJ,IAM1,IAM2,IBM
C      <<< CGN ... C-G(aj,bj,cj;am,bm,cm)
C      =======
C      Subroutine for Clebsch-Gordan calculations
C      Since we need a huge lot of calculations, this sobroutine
C      does not calculate just one. Instead, a set of CG are
C      calculated by using recurrence relations.
C      Specifically, we will calculate C-G(aj,bj,cj;am,bm,cm)
C      For given values of aj,bj,cj,bm and all possible values of am
C      (possible means all am values so that C-G is nonzero, of course)
C      It has been checked to Messiah "Quantum mechanics" Appendix C, Eq. 14a
C      The sum is thumbs-up to within a 10^-12 factor when N2MAX=30
C      Reference: Schulten and Gordon, J. Math. Phys. 16, 1961-1970 (1975)
C--------/---------/---------/---------/---------/---------/---------/--
      implicit none
      real*8 cc1,cc2,dd1,dd2,cgn(200),CX(200),C1,SC,CMAX,CMAX2
      INTEGER IAJ,IBJ,ICJ,IBM,ICM,IAM1,IAM2,IAM,IAMEDIO
      INTEGER DN1,KP,KP2,IAMA,IAMB,DNTOPE,IAMIN,IAMAX,TIC
*
C      First, let's check that icj and ibm falls within limits
      if(icj.lt.abs(iaj-ibj).or.icj.gt.(iaj+ibj)) go to 100
      if(abs(ibm).gt.ibj) go to 200
*
C      IAMIN, IAMAX are the lowest and highest possible values for am
      IAMIN=-MIN(IAJ,ICJ+IBM)
      IAMAX=MIN(IAJ,ICJ-IBM)
C      IAMIN, IAMAX might be such that CG=0
C      If so happens, let's jump up till CG<>0
      KP=MAX(0,IAMIN-IAM1)
      KP2=MAX(0,IAM2-IAMAX)
      DNTOPE=10
      TIC=0
C      Calculation of every possible CG value
      IAMA=IAMIN
      IAMB=IAMAX
      DN1=IAMAX-IAMIN
      CX(1)=1.0D0
      CMAX=CX(1)
      SC=CX(1)*CX(1)
      if(dn1.eq.0) go      to 100
      dd2=dfloat(iaj*(iaj+1)-ibj*(ibj+1)+icj*(icj+1))
      cc2=dfloat((iaj-iama+1)*(iaj+iama))
      cc2=cc2*dfloat((icj-iama-ibm+1)*(icj+iama+ibm))
      cc2=-sqrt(cc2)
C      If there are not many CG to calculate, upward recurrence alone will do
*
      IF(DN1.GE.DNTOPE) IAMB=IAMIN+INT(DN1/2)-1
*
C      Upward recurrence
*
      iam=iama
5      iam=iam+1
            icm=iam+ibm
            if(abs(icm).gt.icj) go to 8
            cc1=cc2
            cc2=dfloat((iaj-iam+1)*(iaj+iam))
            cc2=cc2*dfloat((icj-icm+1)*(icj+icm))
            cc2=-sqrt(cc2)
c      C-G(iaj,iaj,icj/iam,iam,icm)=0 if icj=odd, so...
            if(iaj.eq.ibj.and.iam.eq.ibm.and.mod(icj,2).ne.0) then
                  CX(IAM-IAMIN+1)=0.0d0
                  go to 5
            end if
            dd1=dd2-2.0d0*dfloat((iam-1)*(icm-1))
            CX(IAM-IAMIN+1)=-CX(IAM-IAMIN)*DD1/CC2
            if(dn1.eq.1) then
                  SC=SC+CX(IAM-IAMIN+1)*CX(IAM-IAMIN+1)
                  go to 100
            end if
            if(iam.eq.(iama+1)) go to 7
            CX(IAM-IAMIN+1)=CX(IAM-IAMIN+1)-CX(IAM-IAMIN-1)*CC1/CC2
   7        SC=SC+CX(IAM-IAMIN+1)*CX(IAM-IAMIN+1)
            CMAX2=ABS(CX(IAM-IAMIN+1))
            IF(CMAX2.GT.CMAX) CMAX=CMAX2
*
C      Now let's see if we need go on with upward recurrence
*
            if(iam.lt.iamb) go to 5
            if(tic.eq.1) go to 8
*
c      CX(middle) is used for re-scaling.  If equals zero, then disaster!
c      So let us just check and, if zero, go for the next one.
c      C-G(iaj,iaj,icj/iam,iam,icm)=0 if icj=odd, so...
*
      if(iaj.eq.ibj.and.iam.eq.ibm.and.mod(icj,2).ne.0) go to 5
*
c      Due to round-off errors, we have to be a bit more tolerant: X<10^-14
c      (X being the smallest/largest C-G ratio in our series)
c      However, for high values of iaj,ibj... many C-G are very small, and it
c      becomes difficult to know whether C-G is zero because of roundoff errors
c      or whether that is its real value.
c      My solution: use the next C-G in the series by "jumping" just once
*
      if(dn1.ge.dntope.and.iam.ge.iamb.and.(cmax2/cmax).le.(1.0D-14))
     & then
            tic=1
            go to 5
      end if
  8   continue
      if(dn1.lt.dntope) go to 100
      iamedio=iam
      C1=CX(IAMEDIO-IAMIN+1)
      SC=SC-C1*C1
C      Downward recurrence
      IAMA=IAMEDIO
      iamb=iamax
      CX(IAMB-IAMIN+1)=1.0D0
      cc1=dfloat((iaj-iamb)*(iaj+iamb+1))
      cc1=cc1*dfloat((icj-iamb-ibm)*(icj+iamb+ibm+1))
      cc1=-sqrt(cc1)
      do iam=iamb-1,iama,-1
            icm=iam+ibm
            if(icm.gt.icj) go to 100
            cc2=cc1
            cc1=dfloat((iaj-iam)*(iaj+iam+1))
            cc1=cc1*dfloat((icj-icm)*(icj+icm+1))
            cc1=-sqrt(cc1)
c      C-G(iaj,iaj,icj/iam,iam,icm)=0 if icj=odd, so...
      if(iaj.eq.ibj.and.iam.eq.ibm.and.mod(icj,2).ne.0) then
            CX(IAM-IAMIN+1)=0.0d0
            go to 18
      end if
            if(-icm.gt.icj) go to 18
            dd1=dd2-2.0d0*dfloat((iam+1)*(icm+1))
            CX(IAM-IAMIN+1)=-CX(IAM-IAMIN+2)*DD1/CC1
            if(iam.eq.(iamb-1)) go to 18
            CX(IAM-IAMIN+1)=CX(IAM-IAMIN+1)-CX(IAM-IAMIN+3)*CC2/CC1
  18  end do
      C1=C1/CX(IAMEDIO-IAMIN+1)
C      Now, rescaling
      do iam=IAMEDIO,IAMAX
            cx(iam-iamin+1)=cx(iam-iamin+1)*c1
            sc=sc+cx(iam-iamin+1)*cx(iam-iamin+1)
      end do
C      Here's where the combined (upward+downward) recurrence ends
  100 CONTINUE
      SC=(dfloat(2*ICJ+1)/dfloat(2*IBJ+1))/SC
        SC=SQRT(SC)
C      Now, let's set the sign for SC
      if(cx(iamax-iamin+1).ne.abs(cx(iamax-iamin+1))) SC=-SC
      if(mod((iaj+iamax),2).ne.0) SC=-SC
      if(kp.ne.0) then
            do iam=0,kp-1
                  cgn(iam+1)=0.0d0
            end do
      end if
      if(kp2.ne.0) then
            do iam=iam2-kp2+1,iam2
                  cgn(iam-iam1+1)=0.0d0
            end do
      end if
      do iam=iam1+kp,iam2-kp2
            cgn(iam-iam1+1)=cx(iam-iamin+1)*SC
      end do
  200 CONTINUE
      RETURN
      END
C

      SUBROUTINE LULABY(A,B,NQ,ca,ml,m)           !ca,ml,m are dummy arguments
C--------/---------/---------/---------/---------/---------/---------/--
C      >>> A,NQ,ca,ml,m
C      <<< B
C======

C      Subroutine for LU-based inversion and multiplication
C      of a square nq*nq matrix.
C      The LULABY matrix does the product B*A^(-1), storing it in B
C      The final result (B) is equal to the product (A^-1)*B
C      Since we need to do B*(A^-1) we transpose the A, B matrices:
C      And then we will transpose the final matrix:
C      (At^-1)*Bt = (B*(A^-1))t = -Tt
C      This is a tailored adaptation of the ZGESV (=ZGETRF+ZGETRS)
C      subroutine from the LAPACK package
C--------/---------/---------/---------/---------/---------/---------/--
      implicit none
      INTEGER SU,SUP,ST2                                !,ST3
      PARAMETER (SU=50,SUP=0,ST2=SU+SU+SUP+SUP)
      DOUBLE COMPLEX a(ST2,ST2),b(ST2,ST2),uno,cero,ztemp
      DOUBLE PRECISION SMAX
      INTEGER IPIV(ST2),I,J,K,LL,IMAX,JP,INFO,IP,NT,NQ,ca,ml,m
      uno=(1.0d+0,0.0d+0)
      cero=(0.0d+0,0.0d+0)
      info=0
C      First of all, transpose the A, B, matrices
      do j=1,nq
      do i=1,j-1
            ztemp=a(i,j)
            a(i,j)=a(j,i)
            a(j,i)=ztemp
            ztemp=b(i,j)
            b(i,j)=b(j,i)
            b(j,i)=ztemp
      end do
      end do
C      ZGETRF routine to make an LU decomposition
      DO 100 J=1,NQ
C      Search for pivoting and checking for singularity
C      Sub-subroutine IZAMAX to get the value of i such that a(i,j)=max for a given j
            nt=nq-j+1
            imax=1
            if(nt.eq.1) go to 20
                smax=cdabs(a(j,j))
            do 10 i=2,nt
                  if(abs(a(j+i-1,j)).le.smax) go to 10
                  imax=i
                  smax=cdabs(a(j+i-1,j))
  10        continue
  20        continue
            nt=nq-j
C      End of sub-subroutine IZAMAX
            jp=j-1+imax
            ipiv(j)=jp
            if(a(jp,j).ne.cero) then
C      Sub-subroutine ZSWAP for row interchange
                  if(jp.ne.j) then
                        do i=1,nq
                              ztemp=a(j,i)
                              a(j,i)=a(jp,i)
                              a(jp,i)=ztemp
                        end do
C      End of sub-subroutine ZSWAP
                  end if
C      Sub-subroutine ZSCAL to rescale: divide by por a(j,j)
                  if(j.lt.nq) then
                        ztemp=uno/a(j,j)
                        do i=1,nt
                              a(j+i,j)=ztemp*a(j+i,j)
                        end do
C      End of sub-subroutine ZSCAL
                  end if
C      Now if info=1 that means that U(l,l)=0 and invertion of A is not feasible
            else if(info.eq.0) then
                  info=j
                  write (*,*) 'MATRIX IS SINGULAR - info=',J
                        pause
                  stop
            end if
            if(j.lt.nq) then
C      Sub-subroutine ZGERU para substitute A by A*x*y' (x,y vectors)
C      And therefore update the 'trailing submatrix'
      do k=1,nt
        if(a(j,j+k).ne.cero) then
          ztemp=-uno*a(j,j+k)
            do ll=1,nt
            a(j+ll,j+k)=a(j+ll,j+k)+a(j+ll,j)*ztemp
            end do
          end if
      end do
C      End of sub-subroutine ZGERU
            end if
C      End of the J loop J
 100  continue
C      End of subroutine ZGETRF
C
C      A is now an LU decomposition of a rowwise permutation of A
C      Diagonal unity elements belong tl T.  That is, Lii=1
C
C      Resolving to obtain X = A^-1 * B
c      The procedure is to solve L*X=B, and then U*X=B
C      Sub-subroutine ZLASWP for row interchange in B
      do i=1,nq                 !transposition
            ip=ipiv(i)
            if(ip.ne.i) then
                  do k=1,nq
                        ztemp=b(i,k)
                        b(i,k)=b(ip,k)
                        b(ip,k)=ztemp
                  end do
            end if
c      End of subroutine ZLASWP
      end do
c      Subroutine ZTRSM to solve L*X=B.  The result is overwritten in B
      do j=1,nq
        do k=1,nq
          if(b(k,j).ne.cero) then
            do i=k+1,nq
            b(i,j)=b(i,j)-b(k,j)*a(i,k)
            end do
          end if
        end do
      end do
c      Subroutine ZTRSM to solve U*X=B.  The final result is overwritten in B
      do j=1,nq
        do k=nq,1,-1
          if(b(k,j).ne.cero) then
            b(k,j)=b(k,j)/a(k,k)
            do i=1,k-1
                  b(i,j)=b(i,j)-b(k,j)*a(i,k)
            end do
          end if
        end do
      end do
C      End of subroutine ZTRSM
C      And now, let's transpose again
C      The final result, B, equals the product B*A^(-1)
      do j=1,nq
      do i=1,j-1                !transposition
            ztemp=b(i,j)
            b(i,j)=b(j,i)
            b(j,i)=ztemp
      end do
      end do
C      End of subroutine LULABY

      return
      end
      subroutine medium(ynbrug,nmat,nfin,omega,lambda,rmuf,
     &                  rsnm,omf,ceps1,zeps0,zeps1)
C--------/---------/---------/---------/---------/---------/---------/--
C >>> ynbrug,nmat,nfin,omega,lambda,rmuf,rsnm,omf,ceps1,zeps0
C <<< zeps1
C
C
C--------/---------/---------/---------/---------/---------/---------/--
      implicit none
      integer nmat,nfin,ieps
      real*8 ff,filfrac,pi,xs,lambda,reepsz,rsnm,rmuf,omxf,omega,
     * omxp,plasma
      real*8 omf(nfin)
      complex*16 ci,z1,z2,zeps0,ZEPS1
      complex*16 ceps1(NFIN)
      logical ynbrug

      DATA PI/3.141592653589793d0/
      DATA ci/(0.d0,1.d0)/

* For ideal Drude metal
*     plasma=2.d0*pi*sphere radius in nm/(lambda_z in nm*rmuf)
* where lambda_z is the wavelength for which Re eps_s=0.

      reepsz=2.d0*pi*rsnm/(323.83d0*rmuf)

      IF (NMAT.EQ.1) THEN              !Material decision IF - Drude metal

      plasma=reepsz
        omxp=plasma/omega
        zeps1=1.d0-omxp**2/(1.d0+ci*plasma/(144.d0*omega))
      go to 5
*
      ELSE IF (nmat.eq.4) then             !Material decision IF - ZnS
*
       filfrac=0.62d0         ! filfrac of ZnS in ZnS core

       write(6,*)'Fill. fraction of ZnS in ZnS core=', filfrac

       call  znsrefind(LAMBDA,FILFRAC,zeps1)
       go to 5
*
      ELSE IF (NMAT.EQ.2) THEN         !Material decision IF - Ag

c >>> real material data:           !silver
*                         lambda_z=323.83d0
*                         lambda_p=164.d0
* When real material data are used,
* reepsz differs from plasma!!! The plasma wavelength is
* calculated below:

       plasma=reepsz*7.2d0/3.8291d0

* security trap - remainder (not optimized!)
      omxf=omega/reepsz
      if (omxf.gt.omf(1)) then
       write(6,*)'Calculation of has to stop with'
       write(6,*)' OMF(1)'
       write(6,*)' OMXF=', omxf
       stop
      end if

      if (omxf.lt.omf(nfin)) then
        omxp=plasma/omega
        zeps1=1.d0-omxp**2/(1.d0+ci*plasma/(144.d0*omega))
* damping coefficient for silver is plasma/144 where plasma is different from
* the Re eps zero crossing at 3.8291 eV according to Palik!!!
       go to 5
      else if (omxf.eq.omf(1)) then
       zeps1=ceps1(1)
       go to 5
      else
      do ieps=2,nfin
* data file ordered with the increased wavelength
* omxf increases in the loop and is oriented opposite to the data file
       if (omxf.gt.omf(ieps)) then     ! linear interpolation
       zeps1=ceps1(ieps)+(omxf-omf(ieps))*(ceps1(ieps-1)-ceps1(ieps))
     1 /(omf(ieps-1)-omf(ieps))
       go to 5
       end if
      enddo
       end if   !end Ag

      ELSE IF ((NMAT.GE.3).or.((nmat.ge.5).and.(nmat.le.7))) then   !Material decision IF
                                                                    !Au,Cu,Al,Pt
c >>>
* data file ordered with the decreased wavelength
* omega increases in the loop and is oriented along the data file
*
      if ( (omega.lt.omf(1)).or.(omega.gt.omf(nfin)) ) then
cc       write(6,*)'Material data not available for this wavelength'
cc       stop
*
      call sordalc(NMAT,lambda,ZEPS1)
      go to 5
*
      end if
*
      if (omega.eq.omf(nfin)) then
       zeps1=ceps1(nfin)
       go to 5
      else
      do ieps=1,nfin-1
       if (omega.lt.omf(ieps+1)) then     ! linear interpolation
       zeps1=ceps1(ieps)+(omega-omf(ieps))*(ceps1(ieps+1)-ceps1(ieps))
     1 /(omf(ieps+1)-omf(ieps))
       go to 5
       end if
      enddo
      end if

      ELSE IF (NMAT.EQ.8) then           !Material decision IF - Silicon
c >>>
* data file ordered with the decreased wavelength
* omega increases in the loop and is oriented along the data file
*
      if ( (omega.lt.omf(1)).or.(omega.gt.omf(nfin)) ) then
       write(6,*)'Material data not available for this wavelength'
       stop
*
      end if
*
      if (omega.eq.omf(nfin)) then
       zeps1=ceps1(nfin)
       go to 5
      else
      do ieps=1,nfin-1
       if (omega.lt.omf(ieps+1)) then     ! linear interpolation
       zeps1=ceps1(ieps)+(omega-omf(ieps))*(ceps1(ieps+1)-ceps1(ieps))
     1 /(omf(ieps+1)-omf(ieps))
       go to 5
       end if
      enddo
      end if

      END IF                  ! END of Material decision IF

* The end of reading real data according to Palik's  book
*_____________________________________
* activate Bruggeman:

  5   if (ynbrug) then
      ff=0.8d0

      write(6,*)'Bruggeman with ff=', ff

      z1 = (3.d0*ff-1.d0)*zeps1+(2.d0 - 3.d0*ff)*zeps0
      z2 =  sqrt(z1*z1 + 8.d0*zeps1*zeps0)
*
       if (IMAG(z2).GE.0.0) then
         zeps1= (z1 + z2)/4.d0
       else
         zeps1= (z1 - z2)/4.d0
       end if
       end if
*______________________________________

      RETURN
      END
       subroutine readmat(nmat,nfin,rev,rmuf,omf,ceps1)
C--------/---------/---------/---------/---------/---------/---------/--
C >>> nmat,nfin,rev,rmuf
C <<< omf,ceps1
C         ROUTINE TO READ IN VARIOUS MATERIAL DATA
C--------/---------/---------/---------/---------/---------/---------/--
      implicit none
      integer nmat,nfin,ieps
      real*8 pi,rev,rmuf,omf(NFIN)
      complex*16 ceps1(NFIN),ZEPS1

      DATA PI/3.141592653589793d0/
*****************************   ZEPS1  ***********************************
* Reading real material data, e.g., according to Palik's  book
* requires reading data files OMF and CEPS1 of dimension NFIN
* OMF is reepsz/omega and CEPS1 contains the sphere EPS
*                       material constant reading:
*
      if (nmat.eq.2) then            ! silver data

      OPEN(UNIT=30,FILE='agc.dat')
      rewind(30)
        do ieps=1,nfin
          read(30,*) omf(ieps),ceps1(ieps)
        enddo
       close(30)

      else if (nmat.eq.3) then        ! Gold data

c      OPEN(UNIT=30,FILE='Au293Knew.dat')       !Gold data for different T
      OPEN(UNIT=30,FILE='Aumdat.dat')          !Gold data in nm
      write(6,*)'Gold particles'
      rewind(30)
        do ieps=1, nfin
          read(30,*) omf(ieps),ceps1(ieps)
c          omf(ieps)=2.d0*pi*rev*omf(ieps)/(1240.d0*rmuf)
          omf(ieps)=2.d0*pi*rev/(omf(ieps)*rmuf)
        enddo
       close(30)

cc      else if (nmat.eq.4) then

      else if (nmat.eq.5) then        ! Copper data

      OPEN(UNIT=30,FILE='Cudat.dat')          !Copper data in nm
      write(6,*)'Copper particles'
      rewind(30)
        do ieps=1, nfin
          read(30,*) omf(ieps),ceps1(ieps)
          omf(ieps)=2.d0*pi*rev/(omf(ieps)*rmuf)
        enddo
      close(30)

      else if (nmat.eq.6) then        ! Aluminium data

      OPEN(UNIT=30,FILE='Aldat.dat')          !Aluminium data in nm
      write(6,*)'Aluminum particles'
      rewind(30)
        do ieps=1, nfin
          read(30,*) omf(ieps),ceps1(ieps)
          omf(ieps)=2.d0*pi*rev/(omf(ieps)*rmuf)
        enddo
      close(30)

      else if (nmat.eq.7) then        ! Platinum data

      OPEN(UNIT=30,FILE='Ptdat.dat')          !Platinum data in nm
      write(6,*)'Platinum particles'
      rewind(30)
        do ieps=1, nfin
          read(30,*) omf(ieps),ceps1(ieps)
          omf(ieps)=2.d0*pi*rev/(omf(ieps)*rmuf)
        enddo
      close(30)

      else if (nmat.eq.8) then        ! Silicon data

c     OPEN(UNIT=30,FILE='sieps.dat')  !Silicon data in nm
        OPEN(UNIT=30,FILE='Sidat.dat')   !Silicon data in nm for larger interval
      write(6,*)'Silicon particles'
      rewind(30)
        do ieps=1, nfin
          read(30,*) omf(ieps),ceps1(ieps)
          omf(ieps)=2.d0*pi*rev/(omf(ieps)*rmuf)
        enddo
      close(30)

      end if                      ! material constant reading

*********************
  120 RETURN
      END
        subroutine sordalc(NMAT,LAMBDA,ZEPS)
C----------------------------------------------------------------------
C        SUBROUTINE TO CALCULATE THE DIELECTRIC CONSTANT OF METALS
C             ACCORDING TO AN ARTICLE BY M. A. Ordal et al,
C  Optical properties of fourteen metals in the infrared and far infrared:
C   Al, Co, Cu, Au, Fe, Pb, Mo, Ni, Pd, Pt, Ag, Ti, V, and W,
C   Appl. Opt. {\bf 22}, 1099 (1983); ibid. {\bf 24}, 4493 (1985)
C
C             f77 -g -check_bounds ordalc.f -o rnordalc
C
C   omega=1/\lambda [cm^{-1}]  in Ordal [spectroscopic convention]
C
C   Re eps = - \fr{\om_p^2}{\om^2+\om_\tau^2}
C   Im eps =   \fr{\om_p^2 \om_\tau}{\om^3+\om \om_\tau^2}
C   1eV=1.24 \mu m
C----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER NMAT
      REAL*8 lambda,plasma,tau
      COMPLEX*16 zeps
C                       -------------------------------
      REAL*8 pi,x,y,omega
      DATA PI/3.141592653589793d0/
C   ---------
C ::: speed of light in vacuum in nm/s
C      PARAMETER (c0=2.99792458d17)
C
C According to Table I of Ordal et al, Appl. Opt. {\bf 24}, 4493 (1985):
C                 (plasma=omega_plasma and tau=omega_tau below)
C
C          plasma[THz]/eV/cm-1           tau[THz]/meV/cm-1      LN
C   Al            3570/14.75/119000      19.4/81.8/660          79
C   Cu            1914/7.3890/59600      8.34/9.075/73.2        46
C   Au            2175/9.026/72800       6.5/26.7/215           65
C   Ag            2175/9.013/72700       4.35/18/145            73
C   Pt            1244/5.1450/41500      16.73/69.2/558
C
C ::: conversion factor between normal angular frequency and eV:
c      PARAMETER(XCN=4.13566727d-15)

      if (nmat.eq.3) then           !Au
C ::: Convert the plasma frequency from Table I of Ordal et al
C ::: from [cm-1] to [nm-1] ===> conversion factor 10^{-7}
         PLASMA=72800.d-7            !d12 in Hz/d-7 in [nm-1]
C ::: Convert the tau frequency from Table I of Ordal et al
C ::: from [cm-1] to [nm-1] ===> conversion factor 10^{-7}
          TAU=215.d-7
      else if (nmat.eq.5) then      !Cu
        PLASMA=59600d-7
        TAU=73.2d-7
      else if (nmat.eq.6) then      !Al
        PLASMA=119000d-7
        TAU=660.d-7
      else if (nmat.eq.7) then      !Pt
        PLASMA=41150d-7
        TAU=558.d-7
      end if
C                       -------------------------------

c      write(6,*)'Read in wavelength in nm'
c      read(5,*) lambda
       omega=1.d0/lambda
*
       X=-plasma**2/(omega**2+tau**2)
       Y=tau*plasma**2/(omega**3+omega*tau**2)
*
c      zeps=dcmplx(X,Y)
        zeps =CMPLX(80,0)
*
       END
*
C (C) Copr. 04/2002  Alexander Moroz
      subroutine spherec(lmax,lcs,lambda,rsnm,rmf,zeps,tmt)
C--------/---------/---------/---------/---------/---------/---------/--
C >>> lmax,lcs,lambda,rsnm,rmf,zeps
C <<< tmt
C
C Returns T-matrix
C
C
C                         |  TMT(1,*) |  TMT(4,*)   |
C                 TMT  =  | ----------+-------------|
C                         |  TMT(3,*) |  TMT(2,*)   |
C
C    TMT(1,*) terms corresponds to TEE scattering matrices
C    TMT(2,*) terms corresponds to TMM scattering matrices
C    TMT(3,*) terms corresponds to TME scattering matrices
C    TMT(4,*) terms corresponds to TEM scattering matrices
C    TMT(4,*)=-TMT(3,*)^t where t denotes transposed TMT(3,*) submatrix
C
C    TMT's equal to i*sin(eta)*exp(i*eta), where eta is a phase-shift
C ===============
C  This routines calculates the single sphere scattering properties
C          (including coated spheres)
C
C     xcore ... core size parameter
C     xs    ... sphere size parameter
C     xsa   ... array of dim lcs containing the core and sphere size parameters
C
C
C k_l length in units (2*PI/A=) PI:    xkl= 0.8660254037844386d0
C
C ALPHA, BETA ... arrays of the electric and magnetic phase shifts
C
C LCS ... number of layers of the coated sphere.
C         If lcs=1 - homogeneous sphere
C ILCS ... the coating layer to which material data are read in
C
C ZEPS ... array of dielectric constants of respective shells
C
C Partial wave expansion is used, which is badly convergent
C for large size parameters $x> 100$. In numerical applications, the
C series is to be cut off after
C          LMAX (LMX parameter here) \approx x+4x^{1/3}+2$.
C In the case of LMX=50 this means that x <= 35.
C If one wants to observe ripples, the cutoff for a given x has to
C be even larger
C---------------------------------------------------------------------
      implicit none
      integer LMAXD,LCS,LMAX,LMAXD1,LMTD
      integer NFIN
      real*8 TOL
      character*1 ync
      logical ynperfcon,yntest

c Parameters:
c number of spherical harmonics used
      PARAMETER (LMAXD=50,LMAXD1=LMAXD+1,LMTD=LMAXD1*LMAXD1-1)

c maximal number of sphere coatings
cc      parameter (lcs=3)

* material code number
c   NMAT=0             dispersionless dielectric
c   NMAT=1             Drude metal
c   NMAT=2             Ag
c   NMAT=3             Au
c   NMAT=4             ZnS
c   NMAT=5             Si
*
cc      PARAMETER(NMAT=0)
*
c Temporarily option for reading of the real data for the dielectric constant
c The number of the entries in a material data file to be read below
c          AGC.DAT                NFIN=73       ! from Palik
c          Audat.dat              NFIN=65       ! from Palik
c          Au_2dat.dat            NFIN=76       ! from JAW
c          Au*new.dat             NFIN=142
c          Sidat.dat              NFIN=291
c          sieps.dat              NFIN=66
*
cc      PARAMETER (NFIN=76)
c
C ::: relative error allowed for the TCS. If the convergence
*     within TOL is not reached, program issues warning
      PARAMETER (TOL=1.d-6)
*
c Declarations:
      integer i,j,l,ikl,ij,ij1,m
      real*8 ZEPS0
      real*8 rmf(lcs),rmuf,xs,rsnm,pi
      real*8 omega,lambda
*
      COMPLEX*16 ci,czero,zeps(lcs+1),cqeps(2)
      COMPLEX*16 RX(2),SG(2),zs,zz1,zz2
      COMPLEX*16 KMT(2,lmaxd),TMT(4,LMTD,LMTD),alpha(lmaxd),beta(lmaxd)
*
      COMPLEX*16 cm(lcs,lmaxd),dm(lcs,lmaxd),ce(lcs,lmaxd),
     & de(lcs,lmaxd)
      COMPLEX*16 tt1(2,2,lmaxd,2),tt2(2,2,lmaxd,2)
      COMPLEX*16 AM(lmaxd),AE(lmaxd),BM(lmaxd),BE(lmaxd)
*
      COMPLEX*16 JL(0:lmaxd),NL(0:lmaxd)
      COMPLEX*16 DRJL(0:lmaxd),DRNL(0:lmaxd)
      COMPLEX*16 UL(2,0:lmaxd),VL(2,0:lmaxd)
      COMPLEX*16 DRUL(2,0:lmaxd),DRVL(2,0:lmaxd)
      COMPLEX*16 ZARTAN
*
      external ZARTAN
*
*
cc      COMMON /TOSPHERECC/ zeps
      COMMON /TOSPHERECR/ rmuf
      COMMON /TOSPHERECH/ ync
      COMMON /TOSPHERECL/ ynperfcon
*
* From the main here
*---------------------------------------------------------------
c Data:
      DATA PI/3.141592653589793d0/
      DATA ci/(0.d0,1.d0)/,czero/(0.D0,0.D0)/
C--------/---------/---------/---------/---------/---------/---------/--
* Initialization:

      do ij=1,lmtd
       do ikl=1,lmtd
         tmt(3,ikl,ij)=czero
         tmt(4,ikl,ij)=czero

           if (ikl.ne.ij) then
             tmt(1,ikl,ij)=czero
             tmt(2,ikl,ij)=czero
           end if

       enddo
      enddo

C--------/---------/---------/---------/---------/---------/---------/--
* Checking set up:

      if ((ync.eq.'y'.and.lcs.eq.1).or.(ync.eq.'n'.and.lcs.ne.1)) then
      write(6,*)'Check compatibility of YNC and LCS'
      stop
      end if

C--------------------------------------------------------------------
* Reading in the input data:
c close packed :
c       RMUF=1.d0/DSQRT(2.D0)
*
* size parameter is customarily defined as the ratio of
* circumference of sphere to the wavelength in the host medium
* in which the sphere is embedded
*                    x=kr=\sg a=2*pi*r/lambda
*      xs=2.d0*pi*rs*dble(sqrt(zeps0))/lambda
* convert lambda to the lambda in vacuum:
c      lambda=lambda*dble(sqrt(zeps0))
c  omega=2.d0*pi*rs/(lambda*rmuf)=xs/(rmuf*dble(sqrt(zeps0))),
c where rs is the sphere radius (in nm) and  lambda is the wavelengths (in nm)
c in the vacuum:
         xs=2.d0*pi*rsnm/lambda
         zeps0=zeps(lcs+1)
         omega=xs/(sqrt(zeps0)*rmuf)
*
* Option for omega input:
c      write(6,*)'Read omega ='
c      read(5,*) omega
*
c      xs=RMUF*omega*dble(sqrt(zeps0))
*
c       write(6,*)'Size parameter x=2*pi*rs*n_0/lambda=', xs
*
*                  --------------------------------
*______________________________________
*

  7   if (.not.ynperfcon) then

      ij1=1

      do l=1,lmax
      AM(l)=dcmplx(1.d0,0.d0)
      AE(l)=dcmplx(1.d0,0.d0)
      BM(l)=dcmplx(0.d0,0.d0)
      BE(l)=dcmplx(0.d0,0.d0)
      enddo

      else if (ynperfcon) then

      CQEPS(2)=SQRT(ZEPS(2))
      SG(2)=omega*CQEPS(2)
      RX(1)=SG(2)*RMF(1)
*
*
      call gnzbess(RX(1),lmaxd,jl,drjl,nl,drnl)
*
      DO 10 L=1,lmax
C >>> (AS 10.1.22):
      UL(1,L)=RMF(1)*JL(L)
      VL(1,L)=RMF(1)*NL(L)
      DRJL(L)=SG(2)*DRJL(L)
      DRNL(L)=SG(2)*DRNL(L)
      DRUL(1,L)=JL(L)+RMF(1)*DRJL(L)
      DRVL(1,L)=NL(L)+RMF(1)*DRNL(L)
      AM(l)= NL(L)                ! cm(1,l)
      BM(l)=-JL(L)                ! dm(1,l)
      AE(l)= DRVL(1,L)            ! ce(1,l)
      BE(l)=-DRUL(1,L)            ! de(1,l)

* cf. Jackson 1962, p. 571, Eqs. (16.147);
*                        B/A should yield -tan(phase shift)


  10  continue

      if (lcs.eq.1) go to 30

      ij1=2

      end if
C********************************************************************
c Execution:
* Calculation of the phase shifts
*
      DO 28 j=ij1,lcs

      CQEPS(1)=SQRT(ZEPS(j))
      SG(1)=omega*CQEPS(1)
*
      CQEPS(2)=SQRT(ZEPS(j+1))
      SG(2)=omega*CQEPS(2)
*
      DO 25 I=1,2
*
      RX(I)=SG(I)*RMF(j)
c      WRITE(6,*)'i, rx(i)=', i, rx(i)
C >>>
*
      call gnzbess(RX(I),lmaxd,jl,drjl,nl,drnl)
*
c      write(6,*)'jl=', jl
      DO 15 L=1,lmax
C >>> (AS 10.1.22):
      UL(I,L)=RMF(j)*JL(L)
      VL(I,L)=RMF(j)*NL(L)
      DRJL(L)=SG(I)*DRJL(L)
      DRNL(L)=SG(I)*DRNL(L)
      DRUL(I,L)=JL(L)+RMF(j)*DRJL(L)
      DRVL(I,L)=NL(L)+RMF(j)*DRNL(L)

  15  continue

  25  CONTINUE
*
c      write(6,*)'ul=', ul
*
C >>>  END OF THE LOOP TO ASSIGN VALUES OF BESSEL FUNCTIONS
C      JL and NL start to oscillate after RX.GT. approx 2.5
C********************************************************************
C           Transfer matrix for a layered (coated) sphere
C********************************************************************
*
      do l=1,lmax
*
*   magnetic part
*
      tt1(1,1,l,1)= UL(1,L)
      tt1(1,2,l,1)= VL(1,L)
      tt1(2,1,l,1)= DRUL(1,L)
      tt1(2,2,l,1)= DRVL(1,L)
*
      tt2(1,1,l,1)= sg(2)*DRVL(2,L)
      tt2(1,2,l,1)= - sg(2)*VL(2,L)
      tt2(2,1,l,1)= - sg(2)*DRUL(2,L)
      tt2(2,2,l,1)= sg(2)*UL(2,L)
*
*   electric part
*
      tt1(1,1,l,2)=cqeps(1)*UL(1,L)
      tt1(1,2,l,2)=cqeps(1)*VL(1,L)
      tt1(2,1,l,2)=DRUL(1,L)/cqeps(1)
      tt1(2,2,l,2)= DRVL(1,L)/cqeps(1)
*
      tt2(1,1,l,2)= sg(2)*DRVL(2,L)/cqeps(2)
      tt2(1,2,l,2)= -sg(2)*cqeps(2)*VL(2,L)
      tt2(2,1,l,2)= -sg(2)*DRUL(2,L)/cqeps(2)
      tt2(2,2,l,2)= sg(2)*cqeps(2)*UL(2,L)
*
* m-part
*
      cm(j,l)=AM(l)*(tt2(1,1,l,1)*tt1(1,1,l,1)
     1 +tt2(1,2,l,1)*tt1(2,1,l,1))+BM(l)*(
     2 tt2(1,1,l,1)*tt1(1,2,l,1)+tt2(1,2,l,1)*tt1(2,2,l,1))
*
      dm(j,l)=AM(l)*(tt2(2,1,l,1)*tt1(1,1,l,1)
     1 +tt2(2,2,l,1)*tt1(2,1,l,1))+BM(l)*(
     2 tt2(2,1,l,1)*tt1(1,2,l,1)+tt2(2,2,l,1)*tt1(2,2,l,1))
*
* e-part
*
      ce(j,l)=AE(l)*(tt2(1,1,l,2)*tt1(1,1,l,2)
     1 +tt2(1,2,l,2)*tt1(2,1,l,2))+BE(l)*(
     2 tt2(1,1,l,2)*tt1(1,2,l,2)+tt2(1,2,l,2)*tt1(2,2,l,2))
*
      de(j,l)=AE(l)*(tt2(2,1,l,2)*tt1(1,1,l,2)
     1 +tt2(2,2,l,2)*tt1(2,1,l,2))+BE(l)*(
     2 tt2(2,1,l,2)*tt1(1,2,l,2)+tt2(2,2,l,2)*tt1(2,2,l,2))
*
      AM(l)=cm(j,l)
      BM(l)=dm(j,l)
      AE(l)=ce(j,l)
      BE(l)=de(j,l)
c      write(6,*) AM(l), BM(l)
c      write(6,*) AE(l), BE(l)
*
      enddo
*
  28  CONTINUE
c      write(6,*)'am=', am
c      write(6,*)'bm=', bm
C--------/---------/---------/---------/---------/---------/---------/--
C     ASSIGNING VALUES TO ELEMENTS OF THE K-MATRIX
C >>>
  30  CONTINUE
*
      DO 40 L=1,lmax

* In the following, one needs only phase shifts, so that division
* by SG(2) is omitted:
*
      KMT(1,L)=bm(l)/am(l)
      KMT(2,L)=be(l)/ae(l)
*
 40   CONTINUE
c      write(6,*)'kmt=', kmt
C********************************************************************
      DO 60 j=1,lmax

* >>> Extracting phase-shifts from KMT
*  >>> magnetic part:

         zs=-kmt(1,j)
         alpha(j)=zartan(zs)
*
* Under normal circumstances, Im. part of a phase shift \geq 0 !!!)
*
      if(dimag(alpha(j)).lt.0.d0) then
      write(6,*)'dimag(alpha(j))=',dimag(alpha(j)) ,' is negative'
      write(6,*)'omega, j=', omega, j
      pause
      end if
c      write(6,*)'j, alpha(j)=', j, alpha(j)
*>>> electric part:

         zs=-kmt(2,j)
         beta(j)=zartan(zs)
*
* Under normal circumstances, Im. part of a phase shift \geq 0 !!!)
*
      if(dimag(beta(j)).lt.0.d0) then
      write(6,*)'dimag(alpha(j))=',dimag(beta(j)) ,' is negative'
      write(6,*)'omega, j=', omega, j
      pause
      end if
c      write(6,*)'beta(j)=', beta(j)

  60  CONTINUE

      do 100 l=1,lmax

      zz1 = ci*sin(beta(l))*exp(ci*beta(l))
      zz2 = ci*sin(alpha(l))*exp(ci*alpha(l))

      do 100 m=-l,l

        ij=l*(l+1)+m

        tmt(1,ij,ij)= zz1
        tmt(2,ij,ij)= zz2

 100  continue
*
      return
      end

C (C) Copr. 6/2003  Alexander Moroz

      subroutine sphrd(lambda,xrot,xperp,rev,eps0,zeps)
C--------/---------/---------/---------/---------/---------/---------/--
C      >>>  lambda,xc,xb,eps0,zeps
C      <<<  sext(3)
C
C  xrot ... the half-length of the spheroid along the rotational z-axis
C  xperp ... the half-length of the spheroid along the perpendicular axis
C  rev ... equal-volume-sphere radius
C--------/---------/---------/---------/---------/---------/---------/--
      implicit none
      integer NOUT
C ::: number of the output unit for cross sections and scattering matrix
      PARAMETER (NOUT=35)

      integer ij,ns,npol
      real*8 lambda,xrot,xperp,rev,eps0,pi,sext(6),xk,xmn,xmj,xx,xe,
     &  xlz,xlx,xdz,xdzrn,xdx,xvol,xlp,xdp,xdprn,xcr,xfx,xapl,pf,qf
      complex*16 ci,cone,zalph(6),zeps

      DATA PI/3.141592653589793d0/
      data ci/(0.d0,1.d0)/,cone/(1.d0,0.d0)/

      xk=2.d0*pi*sqrt(eps0)/lambda
      xcr=pi*rev**2        !an effective geom. cross section
      xvol=xrot*xperp**2/3.d0   !V/(4.d0*pi)

      if (xrot.gt.xperp) then
      ns=1         !prolate spheroid
      else
      ns=2         !oblate spheroid
      end if

      xmj=max(xrot,xperp)   !major semiaxis
      xmn=min(xrot,xperp)   !minor semiaxis

      xe=(xmj**2-xmn**2)/xmj**2
      xe=sqrt(xe)            !eccentricity

      if (ns.eq.1) then ! prolate

      xlz=log((1.d0+xe)/(1.d0-xe))
      xlz=(1.d0-xe**2)*(-1.d0 + xlz/(2.d0*xe))/xe**2
*
      xdz=1.d0 + xlz*(1.d0+xe**2)/(1.d0-xe**2)
      xdz=3.d0*xdz/4.d0
      xdzrn=1.d0/xe**2 + (5.d0*xe**2 -3.d0)/((1.d0-xe**2)*xe**2)*xlz
      xdzrn=3.d0*xdzrn/4.d0
      xdx=(3.d0*log((1.d0+xe)/(1.d0-xe))/(2.d0*xe) - xdz)/2.d0

      else if (ns.eq.2) then ! oblate

      xlz=1.d0 - sqrt(1.d0-xe**2)*dasin(xe)/xe
      xlz=xlz/xe**2
*
      xdz=1.d0 + xlz*(1.d0-2.d0*xe**2)
      xdz=3.d0*xdz/4.d0
      xdzrn=(2.d0*xe**2 + 3.d0)*xlz - 1.d0
      xdzrn=3.d0*xdzrn*(1.d0-xe**2)/(4.d0*xe**2)
      xdx=(3.d0*sqrt(1.d0-xe**2)*dasin(xe)/xe - xdz)/2.d0

      end if

      xlx=(1.d0-xlz)/2.d0

      do 20 npol=1,2

      if (npol.eq.1) then        !polarization along the rotation axis
        xlp=xlz
        xdp=xdz
        xdprn=xdzrn
        xfx=1.d0-(xk*xperp)**2/10.d0
        xapl=xrot
      else if (npol.eq.2) then   !perpendicular polarization
        xlp=xlx
        xdp=xdx
        xfx=1.d0-(xk*xrot)**2/10.d0
        xapl=xperp
      end if

* Static Rayleigh polarizability:
      zalph(1)=xvol*(zeps-cone)/(cone+xlp*(zeps-cone))

* MLWA polarizability
      zalph(2)=zalph(1)/(cone-zalph(1)*xk**2/xmj
     & - ci*2.d0*xk**3*zalph(1)/3.d0)

* MLWA polarizability with a ddepol factor:
      zalph(3)=zalph(1)/(cone-xdp*zalph(1)*xk**2/xapl
     &  - ci*2.d0*xk**3*zalph(1)/3.d0)

      pf=0.37d0
      qf=1.d0-pf

* MLWA polarizability with averaged ddepol factor:
      zalph(4)=zalph(1)/(cone-(qf*xdp/xrot +pf/xapl)*zalph(1)*xk**2
     &  - ci*2.d0*xk**3*zalph(1)/3.d0)

* MLWA polarizability with a ddepol factor + averaging:
      zalph(5)=xfx*zalph(1)/(cone-xdp*xfx*zalph(1)*xk**2/xrot
     &  - ci*2.d0*xk**3*xfx*zalph(1)/3.d0)

* MLWA polarizability with a renormalized ddepol factor:
      zalph(6)=zalph(1)/(cone-xdprn*zalph(1)*xk**2/xapl
     &  - ci*2.d0*xk**3*zalph(1)/3.d0)

      do ij=1,6
        sext(ij)=4.d0*pi*xk*imag(zalph(ij))/xcr    ! extinction cross section
*
* xcr=pi*rev**2 - an effective geom. cross section

      if (ij.eq.6)   write(NOUT+18,1108) lambda, sext(6)

      enddo

      if (npol.eq.1) write(NOUT+16,1107) lambda, sext
      if (npol.eq.2) write(NOUT+17,1107) lambda, sext

 20      continue

 1107 FORMAT (F8.2,6(3X,D12.6))
 1108 FORMAT (F8.2,3X,D12.6)

      return
      end

      subroutine sphrint(lmax,lcs,lambda,rsnm,rmf,zeps,am,ae,bm,be)
C--------/---------/---------/---------/---------/---------/---------/--
C >>> lmax,lcs,lambda,rsnm,rmf,zeps
C <<< AM,AE,BM,BE
C
C  Given the incident plane wave expansion coefficients C, returns
C  the arrays A, B for each shell of a multilayered sphere
C
C   TT2 is the forward transfer matrix which translates the regular solution
c       from inside to outside the sphere. Their elements determine the ratio
C                         D/C=T1(21)/T1(11)
C
C   TT1 is the backward transfer matrix which translates the asymptotic solution
c       from outside to inside the sphere. Given the ratio (*), they translate
C       the coefficients A and B backwards into internal layers one by one.
C
C CM,CE... incident plane wave expansion coefficients
C AM,AE,BM,BE ... arrays of the expansion coefficients within the sphere
C
C LCS ... number of layers of the coated sphere.
C         If lcs=1 - homogeneous sphere
C ILCS ... the coating layer to which material data are read in
C
C ZEPS ... array of dielectric constants of respective shells
C
C
C Partial wave expansion is used, which is badly convergent
C for large size parameters $x> 100$. In numerical applications, the
C series is to be cut off after
C          LMAX (LMAXD parameter here) \approx x+4x^{1/3}+2$.
C In the case of LMAXD=50 this means that x <= 35.
C If one wants to observe ripples, the cutoff for a given x has to
C be even larger
C---------------------------------------------------------------------
      implicit none
      integer NOUT,LMAX,LCS,LMAXD
      real*8 TOL,OMEGA,pi,rsnm,rmuf
      logical ynperfcon
C--------/---------/---------/---------/---------/---------/---------/--

C     INPUT REQUIRED:
C--------------------------------------------------------------------
C     PARAMETER INPUT :
C--------------------------------------------------------------------
C ::: number of the output unit

      PARAMETER (NOUT=60)

C  >>> ANGULAR MOMENTUM CUTOFF ON ARRAY DIMENSION
C      (The actual angular-momentum cutoff on summation is specified
C       by the value of variable LMAXD below)

      PARAMETER (LMAXD=60)

C ::: relative error allowed. If the convergence
*     within TOL is not reached, program issues warning

      PARAMETER (TOL=1.d-3)
*
* ynperfcon=.true. if core is a perfect conductor, otherwise
* ynperfcon=.false.
*
      PARAMETER (ynperfcon=.false.)
*
C***********************************************************************
C     >>>               DECLARATIONS:               <<<
C***********************************************************************

      integer i,j,l
      real*8 xs,xd1,xd2,lambda
*
      complex*16 cp,zeps0
      COMPLEX*16 ci,czero
*
      COMPLEX*16 AM(LMAXD,lcs+1),AE(LMAXD,lcs+1),BM(LMAXD,lcs+1),
     & BE(LMAXD,lcs+1)

* Array declarations:

      REAL*8 RMF(lcs)
      COMPLEX*16 RX(2),SG(2),CQEPS(2),zeps(lcs+1)
*
* coated sphere declarations:
*                    moving lmax ===>
*
      COMPLEX*16 cmx11(1:lmaxd),cmx12(1:lmaxd),cmx21(1:lmaxd),
     & cmx22(1:lmaxd)
      COMPLEX*16 tt1(2,2,1:lmaxd,2,lcs),tt2(2,2,1:lmaxd,2,lcs)
*
C--------/---------/---------/---------/---------/---------/---------/--
C  KMT,NML,DML,NEL,DEL,EPS,RX,SG, and Bessel functions below must be
C  declared complex in the case of absorptiom and/or amplification
C--------------------------------------------------------------------
C                      ====================
C                           SUBR MODIF:
C                      =====================

      COMPLEX*16 JL(0:lmaxd), NL(0:lmaxd)
      COMPLEX*16 DRJL(0:lmaxd), DRNL(0:lmaxd)
      COMPLEX*16 UL(2,0:lmaxd),DRUL(2,0:lmaxd)
      COMPLEX*16 WL(2,0:lmaxd), DRWL(2,0:lmaxd)
      COMPLEX*16 HL(0:lmaxd),DRHL(0:lmaxd)
C
C********************************************************************
C     READING THE DATA :
C
      DATA PI/3.141592653589793d0/
      DATA ci/(0.d0,1.d0)/,czero/(0.d0,0.d0)/
      DATA RMUF/1.0D0/
*
C--------------------------------------------------------------------
* Reading in the input data:
c close packed :
c       RMUF=1.d0/DSQRT(2.D0)
*
* size parameter is customarily defined as the ratio of
* circumference of sphere to the wavelength in the host medium
* in which the sphere is embedded
*                    x=kr=\sg a=2*pi*r/lambda
*      xs=2.d0*pi*rs*dble(sqrt(zeps0))/lambda
* convert lambda to the lambda in vacuum:
c      lambda=lambda*dble(sqrt(zeps0))
c  omega=2.d0*pi*rs/(lambda*rmuf)=xs/(rmuf*dble(sqrt(zeps0))),
c where rs is the sphere radius (in nm) and  lambda is the wavelengths (in nm)
c in the vacuum:
         xs=2.d0*pi*rsnm/lambda
         zeps0=zeps(lcs+1)
         omega=xs/(sqrt(zeps0)*rmuf)
*
* Option for omega input:
c      write(6,*)'Read omega ='
c      read(5,*) omega
*
c      xs=RMUF*omega*dble(sqrt(zeps0))
*
c       write(6,*)'Size parameter x=2*pi*rs*n_0/lambda=', xs
*
C--------/---------/---------/---------/---------/---------/---------/--


      omega=2.d0*pi*rsnm/(lambda*rmuf)   !size parameter for default value of rmuf=1

* By activating line below you can feed in a desired sphere size parameter
* to check program results against earlier calculations (for instance those
* of Chew), or you can use an output of other programs to feed in here
* the size parameter corresponding to a Mie resonance

cc      omega=2.51022320078408d0          !temporarily check of Chew's results


C--------/---------/---------/---------/---------/---------/---------/--
*
C********************************************************************
C >>> ASSIGNING THE VALUES FOR BESSEL FUNCTIONS and their
C     derivatives inside and outside the sphere boundary.
C     THE DATA ARE PRODUCED BY THE SUBROUTINE BESSEL(Y,LMAX,PHI,PSI)
C     RECURSION RELATIONS ARE USED ACCORDING (AS 10.1.21-22)
C     NL(=YL) ARE CONSTRUCTED FROM JL BY USING (AS 10.1.15)
C     THE PREFIX DR MEANS THE DERIVATIVE WITH RESPECT TO RMUF AND
C     NOT RX(I)!!
C
*
      DO 28 j=1,lcs

      CQEPS(1)=SQRT(ZEPS(j))
      SG(1)=OMEGA*CQEPS(1)
      CQEPS(2)=SQRT(ZEPS(j+1))
      SG(2)=OMEGA*CQEPS(2)
*
      DO 25 I=1,2
*
      RX(I)=SG(I)*RMF(j)
*
      CALL GNZBESS(RX(I),LMAX,jl,drjl,nl,drnl)
*
      CALL GNRICBESSH(CQEPS(I),OMEGA*rmf(j),LMAX,hl,drhl)
C--------/---------/---------/---------/---------/---------/---------/--
C Returns Riccati-Bessel functions zeta,dzeta of the argument
C RX(i)=CQEPS(i)*OMEGA*rmff(j)
C--------/---------/---------/---------/---------/---------/---------/--
*
      DO 15 L=1,LMAX

      WL(I,L)=HL(L)/SG(I)
      DRWL(I,L)=DRHL(L)

C >>> (AS 10.1.22):
      UL(I,L)=RMF(j)*JL(L)
      DRJL(L)=SG(I)*DRJL(L)
      DRUL(I,L)=JL(L)+RMF(j)*DRJL(L)

  15  continue

  25  CONTINUE
*
C >>>  END OF THE LOOP TO ASSIGN VALUES OF BESSEL FUNCTIONS
C      JL and NL start to oscillate after RX.GT. approx 2.5
C********************************************************************
C Calculation of the regular and asymptotic solutions.
C   TT1 is the backward transfer matrix which translates the asymptotic solution
c       from outside to inside the sphere
C   TT2 is the forward transfer matrix which translates the regular solution
c       from inside to outside the sphere
*
      do 27 l=1,lmax
*
*   magnetic part
*
C--------/---------/---------/---------/---------/---------/---------/--

      tt1(1,1,l,1,j)= -ci*sg(1)*UL(2,L)*WL(1,L)*(DRWL(1,L)/WL(1,L)
     &     - DRUL(2,L)/UL(2,L))
      tt1(1,2,l,1,j)= -ci*sg(1)*WL(2,L)*WL(1,L)*(DRWL(1,L)/WL(1,L)
     &     - DRWL(2,L)/WL(2,L))
      tt1(2,1,l,1,j)= -ci*sg(1)*UL(2,L)*UL(1,L)*(-DRUL(1,L)/UL(1,L)
     &     +DRUL(2,L)/UL(2,L))
      tt1(2,2,l,1,j)= -ci*sg(1)*WL(2,L)*UL(1,L)*(-DRUL(1,L)/UL(1,L)
     &      + DRWL(2,L)/WL(2,L))
*
      tt2(1,1,l,1,j)=-ci*sg(2)*UL(1,L)*WL(2,L)*(DRWL(2,L)/WL(2,L)
     &    - DRUL(1,L)/UL(1,L))
      tt2(1,2,l,1,j)=-ci*sg(2)*WL(1,L)*WL(2,L)*(DRWL(2,L)/WL(2,L)
     &    -  DRWL(1,L)/WL(1,L))
      tt2(2,1,l,1,j)=-ci*sg(2)*UL(1,L)*UL(2,L)*(-DRUL(2,L)/UL(2,L)
     &      +DRUL(1,L)/UL(1,L))
      tt2(2,2,l,1,j)=-ci*sg(2)*WL(1,L)*UL(2,L)*(-DRUL(2,L)/UL(2,L)
     & +DRWL(1,L)/WL(1,L))
*
*   electric part
*

       cp=sqrt(zeps(j)/zeps(j+1))

C--------/---------/---------/---------/---------/---------/---------/--

      tt1(1,1,l,2,j)=-ci*sg(1)*UL(2,L)*WL(1,L)*
     1    (DRWL(1,L)/(cp*WL(1,L)) - cp*DRUL(2,L)/UL(2,L))
      tt1(1,2,l,2,j)=-ci*sg(1)*WL(2,L)*WL(1,L)*
     1    (DRWL(1,L)/(cp*WL(1,L)) - cp*DRWL(2,L)/WL(2,L))
      tt1(2,1,l,2,j)=-ci*sg(1)*UL(2,L)*UL(1,L)*
     1    (-DRUL(1,L)/(cp*UL(1,L)) + cp*DRUL(2,L)/UL(2,L))
      tt1(2,2,l,2,j)=-ci*sg(1)*WL(2,L)*UL(1,L)*
     1    (-DRUL(1,L)/(cp*UL(1,L)) + cp*DRWL(2,L)/WL(2,L))
*
      tt2(1,1,l,2,j)= -ci*sg(2)*UL(1,L)*WL(2,L)*
     1   (cp*DRWL(2,L)/WL(2,L) - DRUL(1,L)/(cp*UL(1,L)))
      tt2(1,2,l,2,j)= -ci*sg(2)*WL(1,L)*WL(2,L)*
     1   (cp*DRWL(2,L)/WL(2,L) - DRWL(1,L)/(cp*WL(1,L)))
      tt2(2,1,l,2,j)= -ci*sg(2)*UL(1,L)*UL(2,L)*
     1   (-cp*DRUL(2,L)/UL(2,L) +DRUL(1,L)/(cp*UL(1,L)))
      tt2(2,2,l,2,j)= -ci*sg(2)*WL(1,L)*UL(2,L)*
     1   (-cp*DRUL(2,L)/UL(2,L) +DRWL(1,L)/(cp*WL(1,L)))
C--------/---------/---------/---------/---------/---------/---------/--
*
  27  CONTINUE       ! over l
*
  28  CONTINUE      ! over shells
*
*    Forward and backward transfer matrices
*
      if (lcs.eq.1) go to 40
*
****
      DO 35 j=2,lcs

      DO 30 l=1,lmax

* =========  TT2 (FORWARD) PART - {cal T}-matrices (Eq. (30)) ===========
* {cal T}-matrices below are defined as {cal T}(n) = \prod_{j=1}^n T^+(j),
* a slight modification of Eq. (30)).
*                  Their elements determine the ratio
C                         D/C=T1(21)/T1(11)
*
* m-part
*
      CMX11(l)=tt2(1,1,l,1,j)*tt2(1,1,l,1,j-1) +
     &             tt2(1,2,l,1,j)*tt2(2,1,l,1,j-1)
      CMX12(l)=tt2(1,1,l,1,j)*tt2(1,2,l,1,j-1) +
     &             tt2(1,2,l,1,j)*tt2(2,2,l,1,j-1)
      CMX21(l)=tt2(2,1,l,1,j)*tt2(1,1,l,1,j-1) +
     &             tt2(2,2,l,1,j)*tt2(2,1,l,1,j-1)
      CMX22(l)=tt2(2,1,l,1,j)*tt2(1,2,l,1,j-1) +
     &             tt2(2,2,l,1,j)*tt2(2,2,l,1,j-1)
*
      tt2(1,1,l,1,j)= CMX11(l)
      tt2(1,2,l,1,j)= CMX12(l)
      tt2(2,1,l,1,j)= CMX21(l)
      tt2(2,2,l,1,j)= CMX22(l)
*
* e-part
*
      CMX11(l)=tt2(1,1,l,2,j)*tt2(1,1,l,2,j-1) +
     &             tt2(1,2,l,2,j)*tt2(2,1,l,2,j-1)
      CMX12(l)=tt2(1,1,l,2,j)*tt2(1,2,l,2,j-1) +
     &             tt2(1,2,l,2,j)*tt2(2,2,l,2,j-1)
      CMX21(l)=tt2(2,1,l,2,j)*tt2(1,1,l,2,j-1) +
     &             tt2(2,2,l,2,j)*tt2(2,1,l,2,j-1)
      CMX22(l)=tt2(2,1,l,2,j)*tt2(1,2,l,2,j-1) +
     &             tt2(2,2,l,2,j)*tt2(2,2,l,2,j-1)
*
      tt2(1,1,l,2,j)= CMX11(l)
      tt2(1,2,l,2,j)= CMX12(l)
      tt2(2,1,l,2,j)= CMX21(l)
      tt2(2,2,l,2,j)= CMX22(l)
*
  30  CONTINUE          ! over l
*
  35  CONTINUE          ! over shells

  40  CONTINUE

      DO l=1,lmax

      AM(l,lcs+1)=1.d0                    !CM(l)
      BM(l,lcs+1)=tt2(2,1,l,1,lcs)/tt2(1,1,l,1,lcs)

      AE(l,lcs+1)=1.d0                    !CE(l)
      BE(l,lcs+1)=tt2(2,1,l,2,lcs)/tt2(1,1,l,2,lcs)

      enddo

      xd1=0.d0          !convergence test parameter
      xd2=0.d0          !convergence test parameter

      DO 50 j=lcs,1,-1

      DO 45 l=1,lmax
*
* =========  TT1  (BACKWARD) PART - {cal M}-matrices  (Eq. (30)) ===========
*      Given the ratio D/C=T1(21)/T1(11), they translate the coefficients
*         A and B backwards into internal layers one by one.
C*
* m-part
*
      AM(l,j)=tt1(1,1,l,1,j)*AM(l,j+1) +
     &             tt1(1,2,l,1,j)*BM(l,j+1)
      BM(l,j)=tt1(2,1,l,1,j)*AM(l,j+1) + tt1(2,2,l,1,j)*BM(l,j+1)
*
* e-part
*
      AE(l,j)=tt1(1,1,l,2,j)*AE(l,j+1) + tt1(1,2,l,2,j)*BE(l,j+1)
      BE(l,j)=tt1(2,1,l,2,j)*AE(l,j+1) + tt1(2,2,l,2,j)*BE(l,j+1)
*
*  consistency check:
      if (j.eq.1) then
      if (xd1.lt.abs(BM(l,1))) xd1=abs(BM(l,1))
        BM(l,1) = czero
      if (xd2.lt.abs(BE(l,1))) xd2=abs(BE(l,1))
         BE(l,1) = czero
      end if

*
  45  CONTINUE          ! over l
*
  50  CONTINUE          ! over shells

      if (xd1.gt.1.d-4) write(6,*)'xd1=',xd1,'.gt.1.d-4'
      if (xd2.gt.1.d-4) write(6,*)'xd2=',xd2,'.gt.1.d-4'


      return
      END
*
      SUBROUTINE TMTR (M,NGAUSS,X,W,AN,ANN,PPI,PIR,PII,R,DR,DDR,
     *                  DRR,DRI,NMAX,NCHECK,NAXSM)
C--------/---------/---------/---------/---------/---------/---------/--
C >>> M,NGAUSS,X,W,AN,ANN,S,SS,PPI,PIR,PII,R,DR,DDR,DRR,DRI,NMAX,NCHECK
C <<< common blocks /TMAT99/, /CT/ (for main),  and /CTT/ (for TT)
C=====================
C
C  Determines the T-matrix of an axially symmetric scatterer
C                           for M.GT.0
C
C  M      - azimuthal number
C  NGAUSS - the number of GIF division points
C  X=\cos\theta  - GIF division points
C  W - GIF weights
C  AN(N)=N*(N+1)
C  ANN(l_1,l_2)=\sqrt{\fr{(2 l_1+1)}{l_1(l_1+1)} }
C                       \sqrt{\fr{(2 l_2+1)}{l_2(l_2+1)} }/2
C  NMAX - angular momentum cutoff
C  NCHECK  -  .EQ.0  THEN  NGSS=2*NGAUSS, FACTOR=1D0
C             .EQ.1  THEN  NGSS = NGAUSS, FACTOR=2D0
C  NAXSM   -  .EQ.0 : Gauss abscissas do not have +/- theta symmetry
C             .EQ.1 : Gauss abscissas have +/- theta symmetry
C  NCHECK - specifies whether NG=2*NGAUSS or otherwise
C  P=DACOS(-1D0)
C  PI=P*2D0/LAM - wave vector
C  PPI=PI*PI
C  PIR=PPI*MRR
C  PII=PPI*MRI
C  R=r^2(\theta)                       for axially symmetric particles
C  DR=[dr(\theta)/(d\theta)]/r(\theta) for axially symmetric particles
C  DDR=\lambda/[2*\pi*r(\theta)]=1/(k_out*r)
C  DRR=(MRR/(MRR**2+MRI**2))*(\lambda/[2*\pi*r(\theta)])
C                  = Re 1/(k_in*r)
C  DRI=-(MRI/(MRR**2+MRI**2))*(\lambda/[2*\pi*r(\theta)])
C                  = Im 1/(k_in*r)
C
C  Refractive index outside is assumed to real, whereas inside
C  a scatterer, refractive index is allowed to be complex in general.
C  Consequently, the Bessel function j_l(k_in*r) will in general
C  be complex. The routine below performs Waterman surface integral
C  separately for the real and imaginary parts of the integrand.
C
C--------/---------/---------/---------/---------/---------/---------/--
*      INCLUDE 'ampld.par.f'
      Parameter (NPN1=100, NPNG1=2600, NPNG2=2*NPNG1, NPN2=2*NPN1,
     &           NPL=NPN2+1, NPN3=NPN1+1,
     &           NPN4=NPN1, NPN5=2*NPN4, NPN6=NPN4+1)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8  X(NPNG2),W(NPNG2),AN(NPN1),
     *        R(NPNG2),DR(NPNG2),SIG(NPN2),
     *        J(NPNG2,NPN1),Y(NPNG2,NPN1),
     *        JR(NPNG2,NPN1),JI(NPNG2,NPN1),DJ(NPNG2,NPN1),
     *        DY(NPNG2,NPN1),DJR(NPNG2,NPN1),
     *        DJI(NPNG2,NPN1),DDR(NPNG2),DRR(NPNG2),
     *        D1(NPNG2,NPN1),D2(NPNG2,NPN1),D3(NPNG2,NPN1),
     *        DRI(NPNG2),RR(NPNG2),
     *        DV1(NPN3),DDV1(NPN3),DV2(NPN3),
     *        DD1,DD2,DD3

      REAL*8  R11(NPN1,NPN1),R12(NPN1,NPN1),
     *        R21(NPN1,NPN1),R22(NPN1,NPN1),
     *        I11(NPN1,NPN1),I12(NPN1,NPN1),
     *        I21(NPN1,NPN1),I22(NPN1,NPN1),
     *        RG11(NPN1,NPN1),RG12(NPN1,NPN1),
     *        RG21(NPN1,NPN1),RG22(NPN1,NPN1),
     *        IG11(NPN1,NPN1),IG12(NPN1,NPN1),
     *        IG21(NPN1,NPN1),IG22(NPN1,NPN1),
     *        ANN(NPN1,NPN1),
     *        QR(NPN2,NPN2),QI(NPN2,NPN2),
     *        RGQR(NPN2,NPN2),RGQI(NPN2,NPN2),
     *        TQR(NPN2,NPN2),TQI(NPN2,NPN2),
     *        TRGQR(NPN2,NPN2),TRGQI(NPN2,NPN2)
cc      REAL*8 TR1(NPN2,NPN2),TI1(NPN2,NPN2)
*________
      COMMON /TMAT99/
     &            R11,R12,R21,R22,I11,I12,I21,I22,RG11,RG12,RG21,RG22,
     &            IG11,IG12,IG21,IG22          !only between TMATR routines
      COMMON /CBESS/ J,Y,JR,JI,DJ,DY,DJR,DJI
cc      COMMON /CT/ TR1,TI1                      !output from TT routine
      COMMON /CTT/ QR,QI,RGQR,RGQI             !input for TT routine
*________
      MM1=M
      QM=DFLOAT(M)
      QMM=QM*QM
      NG=2*NGAUSS
      NM=NMAX+NMAX
      FACTOR=1D0
*
      IF (NCHECK.EQ.1) THEN          !Theta=pi/2 is scatterer mirror symmetry plane
            NGSS=NGAUSS
            FACTOR=2D0
      ELSE IF (NCHECK.EQ.0) THEN     !Theta=pi/2 is not a scatterer mirror symmetry plane
            NGSS=NG
      ENDIF
*
      SI=1D0
      DO 5 N=1,NM                 !NM=2*NMAX
           SI=-SI
           SIG(N)=SI              !=(-1)**N
    5 CONTINUE
*
* Assigning Wigner d-matrices:

      DO 25 I=1,NGAUSS

         I1=NGAUSS-I+1
         I2=NGAUSS+I
*
      CALL VIGF(X(I1),NMAX,M,DV1,DV2,DDV1)
C--------/---------/---------/---------/---------/---------/---------/--
C >>> X,NMAX,M (only nonnegative)
C <<< DV1,DV2,DDV1
C =============
C
C     X=cos(theta), where theta is the polar angle
C     NMAX ... floating  angular momentum cutoff
C
C Returns \pi and \tau scattering functions in terms of the
C Wigner d-functions. Algorithm as described in Eqs. (31-35)
C  of Ref. \cite{Mis39} used. (Note however a missing $n$
C  factor in the 2nd term in the curly bracket in
C   Eq. (35) of Ref. \cite{Mis39}.)
C
C  For a given azimuthal number M.GE.0 returns
C  the Wigner d-functions
C            DV1(N)=dvig(0,m,n,arccos x) = d_{0m}^{(l)}
C
C  \pi scattering function:
C     DDV1(N)=dvig(0,m,n,arccos x)/sin(arccos x)
C                              = m*d_{0m}^{(l)}/ sin\theta
C
C  \tau scattering function:
C     DV2(N)=[d/d(arccos x)] dvig(0,m,n,arccos x)
C                              = d d_{0m}^{(l)}/d\theta
C--------/---------/---------/---------/---------/---------/---------/--
*
         DO N=1,NMAX

            DD1=DDV1(N)
            DD2=DV2(N)
            DD3=DV1(N)
            D1(I1,N)=DD1
            D2(I1,N)=DD2
            D3(I1,N)=DD3

         IF (NAXSM.EQ.1) THEN         !Gauss abscissas chosen +/- symmetric
*
* using (4.2.4) and (4.2.6) of {Ed},
*           d_{0m}^{(l)}(\pi-\theta) = (-1)^{l+m} d_{0m}^{(l)}(\theta)

            SI=SIG(N+M)           !=(-1)**(N+M)
                                  !exactly what follows from {Ed}
            D1(I2,N)=DD1*SI
            D2(I2,N)=-DD2*SI
            D3(I2,N)=DD3*SI

         END IF
         ENDDO

         IF (NAXSM.EQ.0) THEN        !Gauss abscissas not chosen +/- symmetric
*
         CALL VIGF(X(I2),NMAX,M,DV1,DV2,DDV1)
*
          DO N=1,NMAX
            DD1=DDV1(N)
            DD2=DV2(N)
            DD3=DV1(N)
            D1(I2,N)=DD1
            D2(I2,N)=DD2
            D3(I2,N)=DD3
          ENDDO

          END IF

   25 CONTINUE
*
*  Assigning r^2(\theta)*weight product:

      DO 40 I=1,NGSS
           WR=W(I)*R(I)

cc          if (dr(i).eq.0.d0) WR=0.d0   !temporarily only

           RR(I)=WR            !W(I)*r^2(\theta)

   40 CONTINUE
*
      DO 300  N1=MM1,NMAX         !MM1=M below
           AN1=AN(N1)

           DO 300 N2=MM1,NMAX
                AN2=AN(N2)
                AR11=0D0
                AR12=0D0
                AR21=0D0
                AR22=0D0
                AI11=0D0
                AI12=0D0
                AI21=0D0
                AI22=0D0
                GR11=0D0
                GR12=0D0
                GR21=0D0
                GR22=0D0
                GI11=0D0
                GI12=0D0
                GI21=0D0
                GI22=0D0
                SI=SIG(N1+N2)

                DO 200 I=1,NGSS    !=NGAUSS   if NCHECK.EQ.1
                                   !=2*NGAUSS if NCHECK.EQ.0
                    D1N1=D1(I,N1)
                    D2N1=D2(I,N1)
                    D3N1=D3(I,N1)
                    D1N2=D1(I,N2)
                    D2N2=D2(I,N2)
                    D3N2=D3(I,N2)
                    A11=D1N1*D3N2            !pi(N1)*D(N2)
                    A12=D3N1*D2N2            !D(N1)*tau(N2)
                    A21=D2N1*D3N2            !tau(N1)*D(N2)
                    A22=D2N1*D2N2            !tau(N1)*tau(N2)
                    AA1=D1N1*D2N2+D2N1*D1N2  !pi(N1)*tau(N2)+tau(N1)*pi(N2)
                    AA2=D1N1*D1N2 +A22       !pi(N1)*pi(N2)+tau(N1)*tau(N2)

* Vector spherical harmonics:
C  Since refractive index is allowed to be complex in general,
C  the Bessel function j_l(k_in*r) is complex. The code below
C  performs a separation of the complex integrand in Waterman's
C  surface integral into its respective real and imaginary
C  parts:

* Bessel functions of the exterior argument:

                    QJ1=J(I,N1)
                    QY1=Y(I,N1)
                    QDJ1=DJ(I,N1)
                    QDY1=DY(I,N1)

* Bessel functions of the interior argument:

                    QJR2=JR(I,N2)
                    QJI2=JI(I,N2)
                    QDJR2=DJR(I,N2)
                    QDJI2=DJI(I,N2)

* Re and Im of j_{n2}(k_{in}r) j_{n1}(k_{out}r):

                    C1R=QJR2*QJ1
                    C1I=QJI2*QJ1

* Re and Im of j_{n2}(k_{in}r) h_{n1}(k_{out}r):

                    B1R=C1R-QJI2*QY1
                    B1I=C1I+QJR2*QY1

* Re and Im of j_{n2}(k_{in}r) j_{n1}'(k_{out}r)/(k_{out}r):

                    C2R=QJR2*QDJ1
                    C2I=QJI2*QDJ1

* Re and Im of j_{n2}(k_{in}r) h_{n1}'(k_{out}r)/(k_{out}r):

                    B2R=C2R-QJI2*QDY1
                    B2I=C2I+QJR2*QDY1

                    DDRI=DDR(I)               !1/(k_{out}r)

* Re and Im of [1/(k_{out}r)]*j_{n2}(k_{in}r) j_{n1}(k_{out}r)

                    C3R=DDRI*C1R
                    C3I=DDRI*C1I

* Re and Im of [1/(k_{out}r)]*j_{n2}(k_{in}r) h_{n1}(k_{out}r):

                    B3R=DDRI*B1R
                    B3I=DDRI*B1I

* Re and Im of [k_{in}r j_{n2}(k_{in}r)]'/(k_{in}r)
*                          * j_{n1}(k_{out}r):

                    C4R=QDJR2*QJ1
                    C4I=QDJI2*QJ1

* Re and Im of [k_{in}r j_{n2}(k_{in}r)]'/(k_{in}r)
*                          *  h_{n1}(k_{out}r):

                    B4R=C4R-QDJI2*QY1
                    B4I=C4I+QDJR2*QY1

                    DRRI=DRR(I)               !Re[1/(k_{in}r)]
                    DRII=DRI(I)               !Im[1/(k_{in}r)]

* Re and Im of [1/(k_{in}r)] j_{n2}(k_{in}r) j_{n1}(k_{out}r):

                    C5R=C1R*DRRI-C1I*DRII
                    C5I=C1I*DRRI+C1R*DRII

* Re and Im of [1/(k_{in}r)] j_{n2}(k_{in}r) h_{n1}(k_{out}r):

                    B5R=B1R*DRRI-B1I*DRII
                    B5I=B1I*DRRI+B1R*DRII


* Re and Im of j_{n2}'(k_{in}r) j_{n1}'(k_{out}r):

                    C6R=QDJR2*QDJ1
                    C6I=QDJI2*QDJ1

* Re and Im of j_{n2}'(k_{in}r) h_{n1}'(k_{out}r):

                    B6R=C6R-QDJI2*QDY1
                    B6I=C6I+QDJR2*QDY1

* Re and Im of [1/(k_{out}r)] j_{n2}'(k_{in}r) j_{n1}(k_{out}r):

                    C7R=C4R*DDRI
                    C7I=C4I*DDRI

* Re and Im of [1/(k_{out}r)] j_{n2}'(k_{in}r) h_{n1}(k_{out}r):

                    B7R=B4R*DDRI
                    B7I=B4I*DDRI

* Re and Im of [1/(k_{in}r)] j_{n2}(k_{in}r) j_{n1}'(k_{out}r):

                    C8R=C2R*DRRI-C2I*DRII
                    C8I=C2I*DRRI+C2R*DRII

* Re and Im of [1/(k_{in}r)] j_{n2}(k_{in}r) h_{n1}'(k_{out}r):

                    B8R=B2R*DRRI-B2I*DRII
                    B8I=B2I*DRRI+B2R*DRII


* %%%%%%%%%  Forming integrands of J-matrices (J^{11}=J^{22}=0 for m.eq.0):

                    URI=DR(I)
                    RRI=RR(I)

                    IF (NCHECK.EQ.1.AND.SI.GT.0D0) GO TO 150

* W(I)*r^2(I)*(pi(N1)*tau(N2)+tau(N1)*pi(N2):

                    E1=RR(I)*AA1             ! <-- AA1

                    AR11=AR11+E1*B1R
                    AI11=AI11+E1*B1I
                    GR11=GR11+E1*C1R
                    GI11=GI11+E1*C1I

                    IF (NCHECK.EQ.1) GO TO 160

  150               CONTINUE


* w(i)*r^2(\theta)*[pi(N1)*pi(N2)+tau(N1)*tau(N2)]
* (prefactor containing r^2(\theta)<->hat{r} part)

                    F1=RRI*AA2             ! <-- AA2

* N1*(N1+1)*w(i)*r(\theta)*[dr/(d\theta)]*D(N1)*tau(N2):
*  (prefactor containing r(\theta)*[dr/(d\theta)] - hat{theta} part)

                    F2=RRI*URI*AN1*A12             ! <-- A12

                    AR12=AR12+F1*B2R+F2*B3R        !~Re J^{12}
                    AI12=AI12+F1*B2I+F2*B3I        !~Im J^{12}

                    GR12=GR12+F1*C2R+F2*C3R        !~Re Rg J^{12}
                    GI12=GI12+F1*C2I+F2*C3I        !~Im Rg J^{12}

* N2*(N2+1)*w(i)*r(\theta)*[dr/(d\theta)]*tau(N1)*D(N2):
* (!prefactor containing r(\theta)*[dr/(d\theta)] - hat{theta} part)

                    F2=RRI*URI*AN2*A21             ! <-- A21

                    AR21=AR21+F1*B4R+F2*B5R
                    AI21=AI21+F1*B4I+F2*B5I

                    GR21=GR21+F1*C4R+F2*C5R
                    GI21=GI21+F1*C4I+F2*C5I

                    IF (NCHECK.EQ.1) GO TO 200

  160  CONTINUE

* w(i)*r^2(\theta)*[dr/(d\theta)]*pi(N1)*D(N2):
* (!prefactor containing r^2(\theta)*[dr/(d\theta)] - hat{theta} part)

                    E2=RRI*URI*A11
                    E3=E2*AN2
                    E2=E2*AN1

                    AR22=AR22+E1*B6R+E2*B7R+E3*B8R
                    AI22=AI22+E1*B6I+E2*B7I+E3*B8I

                    GR22=GR22+E1*C6R+E2*C7R+E3*C8R
                    GI22=GI22+E1*C6I+E2*C7I+E3*C8I

  200           CONTINUE           !Gauss integration

*%%%%%%%%%%%%%  Forming J-matrices (J^{11}=J^{22}=0 for m.eq.0):

                AN12=ANN(N1,N2)*FACTOR

                R11(N1,N2)=AR11*AN12       !Re J^{11}
                R12(N1,N2)=AR12*AN12       !Re J^{12}
                R21(N1,N2)=AR21*AN12       !Re J^{21}
                R22(N1,N2)=AR22*AN12       !Re J^{22}
                I11(N1,N2)=AI11*AN12       !Im J^{11}
                I12(N1,N2)=AI12*AN12       !Im J^{12}
                I21(N1,N2)=AI21*AN12       !Im J^{21}
                I22(N1,N2)=AI22*AN12       !Im J^{22}

                RG11(N1,N2)=GR11*AN12       !Re (Rg J^{11})
                RG12(N1,N2)=GR12*AN12       !Re (Rg J^{12})
                RG21(N1,N2)=GR21*AN12       !Re (Rg J^{21})
                RG22(N1,N2)=GR22*AN12       !Re (Rg J^{22})
                IG11(N1,N2)=GI11*AN12       !Im (Rg J^{11})
                IG12(N1,N2)=GI12*AN12       !Im (Rg J^{12})
                IG21(N1,N2)=GI21*AN12       !Im (Rg J^{21})
                IG22(N1,N2)=GI22*AN12       !Im (Rg J^{22})

  300 CONTINUE

*%%%%%%%%%%%%%%%%%%%%%%%  Forming Q and RgQ -matrices

      TPIR=PIR                 !Re [1/k_{in}^2]
      TPII=PII                 !Im [1/k_{in}^2]
      TPPI=PPI                 !1/k_{out}^2

      NM=NMAX-MM1+1
      DO 310 N1=MM1,NMAX
           K1=N1-MM1+1                !from 1 to NMAX-MM1+1
           KK1=K1+NM                  !from NMAX-MM1+2 to 2*(NMAX-MM1+1)

           DO 310 N2=MM1,NMAX
                K2=N2-MM1+1           !from 1 to NMAX-MM1+1
                KK2=K2+NM             !from NMAX-MM1+2 to 2*(NMAX-MM1+1)

                TAR11=-R11(N1,N2)
                TAI11=-I11(N1,N2)
                TGR11=-RG11(N1,N2)
                TGI11=-IG11(N1,N2)

                TAR12= I12(N1,N2)
                TAI12=-R12(N1,N2)
                TGR12= IG12(N1,N2)
                TGI12=-RG12(N1,N2)

                TAR21=-I21(N1,N2)
                TAI21= R21(N1,N2)
                TGR21=-IG21(N1,N2)
                TGI21= RG21(N1,N2)

                TAR22=-R22(N1,N2)
                TAI22=-I22(N1,N2)
                TGR22=-RG22(N1,N2)
                TGI22=-IG22(N1,N2)

                TQR(K1,K2)=TPIR*TAR21-TPII*TAI21+TPPI*TAR12
                TQI(K1,K2)=TPIR*TAI21+TPII*TAR21+TPPI*TAI12
                TRGQR(K1,K2)=TPIR*TGR21-TPII*TGI21+TPPI*TGR12
                TRGQI(K1,K2)=TPIR*TGI21+TPII*TGR21+TPPI*TGI12

                TQR(K1,KK2)=TPIR*TAR11-TPII*TAI11+TPPI*TAR22
                TQI(K1,KK2)=TPIR*TAI11+TPII*TAR11+TPPI*TAI22
                TRGQR(K1,KK2)=TPIR*TGR11-TPII*TGI11+TPPI*TGR22
                TRGQI(K1,KK2)=TPIR*TGI11+TPII*TGR11+TPPI*TGI22

                TQR(KK1,K2)=TPIR*TAR22-TPII*TAI22+TPPI*TAR11
                TQI(KK1,K2)=TPIR*TAI22+TPII*TAR22+TPPI*TAI11
                TRGQR(KK1,K2)=TPIR*TGR22-TPII*TGI22+TPPI*TGR11
                TRGQI(KK1,K2)=TPIR*TGI22+TPII*TGR22+TPPI*TGI11

                TQR(KK1,KK2)=TPIR*TAR12-TPII*TAI12+TPPI*TAR21
                TQI(KK1,KK2)=TPIR*TAI12+TPII*TAR12+TPPI*TAI21
                TRGQR(KK1,KK2)=TPIR*TGR12-TPII*TGI12+TPPI*TGR21
                TRGQI(KK1,KK2)=TPIR*TGI12+TPII*TGR12+TPPI*TGI21

  310 CONTINUE

      NNMAX=2*NM
      DO 320 N1=1,NNMAX
           DO 320 N2=1,NNMAX
                QR(N1,N2)=TQR(N1,N2)
                QI(N1,N2)=TQI(N1,N2)
                RGQR(N1,N2)=TRGQR(N1,N2)
                RGQI(N1,N2)=TRGQI(N1,N2)
  320 CONTINUE
*
      CALL TT(NM,NCHECK)
*
      RETURN
      END

       SUBROUTINE VIGF(X,LMAX,M,DV1,DV2,DDV1)
C--------/---------/---------/---------/---------/---------/---------/--
C >>> X,LMAX,M (only nonnegative)
C <<< DV1,DV2,DDV1
C =============
C
C     X=cos(theta), where theta is the polar angle
C     LMAXD ... maximal angular momentum cutoff
C     LMAX ... floating  angular momentum cutoff
C
C Returns \pi and \tau scattering functions in terms of the
C Wigner d-functions. Algorithm as described in Eqs. (31-35)
C  of Ref. \cite{Mis39} used. (Note however a missing $n$
C  factor in the 2nd term in the curly bracket in
C   Eq. (35) of Ref. \cite{Mis39}.)
C
C     For a given azimuthal number M.GE.0 returns
C      the Wigner d-functions
C            DV1(N)=dvig(0,m,n,arccos x) = d_{0m}^{(l)}
C
C  \pi scattering function:
C     DDV1(N)=dvig(0,m,n,arccos x)/sin(arccos x)
C                              = m*d_{0m}^{(l)}/ sin\theta
C
C  \tau scattering function:
C     DV2(N)=[d/d(arccos x)] dvig(0,m,n,arccos x)
C                              = d d_{0m}^{(l)}/d\theta
C
C     for 1.LE.N.LE.LMAX and 0.LE.X.LE.1
C      DDV1 is calculated because (DV1/sin\theta) is singular for
C             either \beta=0 or \beta=\pi
C     (For a given M.NEQ.0, only the M.LE.N.LE.LMAX terms are determined!)
C =====
C
C     In the present case (Eq. (B.28) of Ref. \ct{MTL})
C                       (cf. (4.1.24) of \ct{Ed}):
C
C             d_{00}^{(l)}(\theta)= P_l(\cos\theta)
C      d_{0m}^{(l)}(\theta)= (-1)^m \sqrt{(l-m)!/(l+m)!}P_l^m(\cos\theta)
C     (assuming Jackson's P^m_l), where $d^{(l)}_{m0}=(-1)^m d^{(l)}_{0m}$)
C     (Edmonds P^m_l is without (-1)**m prefactor; cf. (4.1.24) therein)
C
C Special values:
C
C     (Rodrigues formula [Eq. (2.5.14) of Ref. \ct{Ed}] then yields
C                       P_1(x)=x; P_2=(3x^2-1)/2; etc.
C
C      P_l^m(x) = (-1)^m (1-x^2)^{m/2} \fr{d^m P_l(x)}{dx} ===>
C                       P_1^1(\cos\beta)=-\sin\beta
C     Therefore,
C              d_{00}^{(1)}(\beta)=\cos\beta
C              d_{01}^{(1)}(\beta)=\sin\beta/\sqrt{2}
C         d d_{00}^{(1)}(\beta)/d\beta=-\sin\beta
C         d d_{01}^{(1)}(\beta)/d\beta=\cos\beta/\sqrt{2}
C
C     Acc. Eq. (34) of {Mis39}:
C
C     A_0=1, A_1=1/\sqrt{2}, A_2=\sqrt{3}/(2*\sqrt{2})
C
C     Therefore [Eq. (32) of {Mis39}]:
C              d_{00}^{(0)}(\beta)=1
C              d_{01}^{(1)}(\beta)=\sin\beta/\sqrt{2}
C              d_{02}^{(2)}(\beta)=\sqrt{3}\sin^2\beta/(2*\sqrt{2})
C     and
C         d d_{00}^{(0)}(\beta)/d\beta=0
C         d d_{01}^{(1)}(\beta)/d\beta=\cos\beta/\sqrt{2}
C         d d_{02}^{(2)}(\beta)/d\beta=\sqrt{3}\sin\beta \cos\beta/\sqrt{2}
C                                = \sqrt{3}\sin (2\beta) /(2*\sqrt{2})
C =====
C     Similar to routine VIG, which however only returns
C            DV1(N)=dvig(0,m,n,arccos x) = d_{0m}^{(l)}
C
C     When arccos x is very small, a care has to be exercise to generate
C     nboth DDV1 and DV2. That part has been made using recurrences of
C     Ref. \ct{TKS}
C--------/---------/---------/---------/---------/---------/---------/--

      IMPLICIT none
      INTEGER NPN1,NPNG1,NPNG2,NPN2,NPL,NPN3,NPN4,NPN5,NPN6
      Parameter (NPN1=100, NPNG1=2600, NPNG2=2*NPNG1, NPN2=2*NPN1,
     &           NPL=NPN2+1, NPN3=NPN1+1,
     &           NPN4=NPN1, NPN5=2*NPN4, NPN6=NPN4+1)

      integer n,LMAX,M,I,I2
      REAL*8 A,X,QS,D1,D2,D3,DER,DN,DX,QN,QN1,QN2,
     & QNM,QNM1,QMM
      REAL*8 DDV1(NPN6), DV1(NPN6), DV2(NPN6)

* DV1, DDV1, and DV2 initialization
      DO 1 N=1,LMAX
         DV1(N) =0.D0
         DDV1(N)=0.D0
         DV2(N) =0.D0
    1 CONTINUE

      DX=DABS(X)
      A=1.D0                       !initial A_0
      QS=DSQRT(1D0-X*X)            !sin\theta
***********************************************************************
*                      NONZERO DV1 INITIALIZATION
*
      IF (M.NE.0) GO TO 20
*
* DDV1(N)=0.d0 [see (3.33) of {TKS}]
* D1,D2, and D3 below are the three consequent terms
*       d_{0m}^{n-1}, d_{0m}^{n}, and d_{0m}^{n+1} beginning
*       with n=m
*=============
* Recurrence initialization following d^l_{00}=P_l
*         [see Eq. (B.27) of {MTL}]
*
      D1=1.D0                  !d^0_{00}=P_0=1   (see Sec. 6.8 of {NR})
      D2=X                     !d^1_{00}=P_1=x

      DO 5 N=1,LMAX
         QN=DFLOAT(N)
         QN1=DFLOAT(N+1)
         QN2=DFLOAT(2*N+1)
         D3=(QN2*X*D2-QN*D1)/QN1       !recurrence (31) of Ref. {Mis39} for d^{N+1}_{00}
         DV1(N)=D2                     !d^N_{00}
         D1=D2                         !becomes d^{N-1}_{00} in D3
         D2=D3                         !becomes d^{N}_{00} in D3
  5   CONTINUE
*
      go to 100

***********************************************************
*                           M\neq 0 part
*  d_{0m}^m \equiv A_m*(sin\theta)**m   initialization - (33) and recurrence (34) of Ref. {Mis39}

   20 CONTINUE
      DO 25 I=1,M
         I2=I*2
         A=A*DSQRT(DFLOAT(I2-1)/DFLOAT(I2))*QS  !recurrence (33,34) of Ref. {Mis39} f
   25 CONTINUE

*
* Recurrence initialization following Eqs. (32,33) of Ref. {Mis39}

      D1=0.D0                 !=DV1(M-1); see Eq. (32) of Ref. {Mis39}
      D2=A                    !=DV1(M);   see Eq. (33) of Ref. {Mis39}
      QMM=DFLOAT(M*M)
*
      DO 30 N=M,LMAX
         QN=DFLOAT(N)
         QN1=DFLOAT(N+1)
         QN2=DFLOAT(2*N+1)
         QNM=DSQRT(QN*QN-QMM)
         QNM1=DSQRT(QN1*QN1-QMM)
         D3=(QN2*X*D2-QNM*D1)/QNM1             !recurrence (31) of Ref. {Mis39} for d^{N+1}_{0M}
         DV1(N)=D2                             !d^N_{0M}
         D1=D2                                 !becomes d^{N-1}_{0M} in D3,DER
         D2=D3                                 !becomes d^{N}_{0M} in D3
   30 CONTINUE

      go to 100

*                        DV1 INITIALIZED
*             It remains to determine DDV1 and DV2
*********************************************************
*  (1-cos\theta) is very small:
*
C   For theta=0 [see Eqs. above]:
C              d_{00}^{(0)}(0)=1
C              d_{01}^{(1)}(0)=0
C              d_{02}^{(2)}(\beta)=0
C     and
C         d d_{00}^{(0)}(\beta)/d\beta=0
C         d d_{01}^{(1)}(\beta)/d\beta=1/\sqrt{2}
C         d d_{02}^{(2)}(\beta)/d\beta=0
C
C  See Eqs. (4.1-4) of \ct{Mis91}:
C
C   (m/\sin\theta) d_{0m}^l(0)=(\delta_{m\pm 1}/2) \sqrt{l(l+1)}
C      d d_{0m}^l(0)/d\beta   =(m\delta_{m\pm 1}/2) \sqrt{l(l+1)}
C
*
*  (4.2.1) of \ct{Ed}:
*   d_{0m}^{(l)}(pi) = (-1)^{l+m} \dt_{0,m}
*
*  (4.2.3) of \ct{Ed}:
*   d_{0m}^{(l)}(0) = (-1)^{m} \dt_{0,m} = \dt_{0,m}
*=======================================
*
*  If X^l_m=(m/\sin\theta) d_{0m}^{(l)}, then, according to (3.29) of {TKS}:
*
*  X^{m+1}_{m+1}=\sin\theta \sqrt{\fr{2m+1}{2m+2}}
*                           \left(\fr{m+1}{m}\right)X^{m}_{m}
*
*  According to (3.30) of {TKS}:
*  X^{m+1}_{m}= -\sqrt{2m+1}\,\cos\theta X^{m}_{m}
*
* According to (3.31) of {TKS}:
*  X^{l}_{m}=\fr{1}{\sqrt{l^2-m^2}}\,\left[(2l-1)\cos\theta
*          X^{l-1}_{m} - \sqrt{(l-1)^2-m^2}}\,\X^{l-2}_{m} \right]
*
* Initial recurrence values are X^1_1=\sqrt{2}/2 and X^l_0=0
***********************************************************************
*                   NONZERO DDV1/DV2 INITIALIZATION
*
*                          M = 0

 100  IF (M.EQ.0) THEN     !all DDV1(N)=X^l_0=0; see (3.33) of {TKS}:

* According to (3.37) of {TKS}, DV2(0)=0.d0

      DV2(1)=-QS

      IF (LMAX.GE.2) DV2(2)=3*X*DV2(1)

      IF (LMAX.LT.3) RETURN
*
      DO N=3,LMAX           !recurrence (3.36) of {TKS},
      DV2(N)=(2*N-1)*X*DV2(N-1)/(N-1)-N*DV2(N-2)/(N-1)
      ENDDO
***********************************************************************
*                           M > 0

       ELSE IF (M.GT.0) THEN
*
* >>> Determine X^m_m according to Eq. (3.29) of {TKS}:

      A=1.d0/DSQRT(2.D0)               !X^1_1=A_1

      DO I=1,M-1
      A=QS*DBLE(I+1)*DSQRT(2*I+1.d0)*A/(I*DSQRT(2*I+2.d0))
      ENDDO

* <<< A is now X^m_m; see (3.29) of {TKS}

      DDV1(M)=A
      DV2(M)=X*A                        !see (3.34) of {TKS}

* >>> Determine X^{m+1}_m:

      IF (M.EQ.LMAX)  GO TO 120

      DER=X*DSQRT(2*M+1.d0)*A          ! DER=X^{m+1}_m; see (3.30) of {TKS}
      DDV1(M+1)=DER
      DV2(M+1)=((M+1)*X*DER-A*DSQRT(2*M+1.d0))/DBLE(M)  !(3.35) of {TKS}

* >>> Determine remaining X^{l}_m's

      IF ((M+2).EQ.LMAX)  GO TO 120

       DO N=M+2,LMAX
       D3=DSQRT(DBLE(N)**2-DBLE(M)**2)
       DDV1(N)=((2*N-1)*X*DDV1(N-1) -
     &                DSQRT(DBLE(N-1)**2-DBLE(M)**2)*DDV1(N-2))/D3
                                                      !see (3.31) of {TKS}
       DV2(N)=(N*X*DDV1(N)-DDV1(N-1)*D3)/DBLE(M)      !see (3.35) of {TKS}
       ENDDO

      END IF

cv  100 IF (M.NE.1) RETURN
cv
cv      DO 110 N=1,LMAX
cv         DN=DFLOAT(N*(N+1))
cv         DN=0.5D0*DSQRT(DN)
cv         IF (X.LT.0D0) DN=DN*(-1)**(N+1)
cv         DV1(N)=DN
cv         IF (X.LT.0D0) DN=-DN
cv         DV2(N)=DN
cv  110 CONTINUE

  120 RETURN
      END


      DOUBLE PRECISION FUNCTION D1MACH(I)
C***BEGIN PROLOGUE  D1MACH
C***DATE WRITTEN   750101   (YYMMDD)
C***REVISION DATE  890213   (YYMMDD)
C***CATEGORY NO.  R1
C***KEYWORDS  LIBRARY=SLATEC,TYPE=DOUBLE PRECISION(R1MACH-S D1MACH-D),
C             MACHINE CONSTANTS
C***AUTHOR  FOX, P. A., (BELL LABS)
C           HALL, A. D., (BELL LABS)
C           SCHRYER, N. L., (BELL LABS)
C***PURPOSE  Returns double precision machine dependent constants
C***DESCRIPTION
C
C   D1MACH can be used to obtain machine-dependent parameters
C   for the local machine environment.  It is a function
C   subprogram with one (input) argument, and can be called
C   as follows, for example
C
C        D = D1MACH(I)
C
C   where I=1,...,5.  The (output) value of D above is
C   determined by the (input) value of I.  The results for
C   various values of I are discussed below.
C
C   D1MACH( 1) = B**(EMIN-1), the smallest positive magnitude.
C   D1MACH( 2) = B**EMAX*(1 - B**(-T)), the largest magnitude.
C   D1MACH( 3) = B**(-T), the smallest relative spacing.
C   D1MACH( 4) = B**(1-T), the largest relative spacing.
C   D1MACH( 5) = LOG10(B)
C
C   Assume double precision numbers are represented in the T-digit,
C   base-B form
C
C              sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
C
C   where 0 .LE. X(I) .LT. B for I=1,...,T, 0 .LT. X(1), and
C   EMIN .LE. E .LE. EMAX.
C
C   The values of B, T, EMIN and EMAX are provided in I1MACH as
C   follows:
C   I1MACH(10) = B, the base.
C   I1MACH(14) = T, the number of base-B digits.
C   I1MACH(15) = EMIN, the smallest exponent E.
C   I1MACH(16) = EMAX, the largest exponent E.
C
C   To alter this function for a particular environment,
C   the desired set of DATA statements should be activated by
C   removing the C from column 1.  Also, the values of
C   D1MACH(1) - D1MACH(4) should be checked for consistency
C   with the local operating system.
C
C***REFERENCES  FOX P.A., HALL A.D., SCHRYER N.L.,*FRAMEWORK FOR A
C                 PORTABLE LIBRARY*, ACM TRANSACTIONS ON MATHEMATICAL
C                 SOFTWARE, VOL. 4, NO. 2, JUNE 1978, PP. 177-188.
C***ROUTINES CALLED  XERROR
C***END PROLOGUE  D1MACH
C
      INTEGER SMALL(4)
      INTEGER LARGE(4)
      INTEGER RIGHT(4)
      INTEGER DIVER(4)
      INTEGER LOG10(4)
C
      DOUBLE PRECISION DMACH(5)
      SAVE DMACH
C
      EQUIVALENCE (DMACH(1),SMALL(1))
      EQUIVALENCE (DMACH(2),LARGE(1))
      EQUIVALENCE (DMACH(3),RIGHT(1))
      EQUIVALENCE (DMACH(4),DIVER(1))
      EQUIVALENCE (DMACH(5),LOG10(1))
C
C     MACHINE CONSTANTS FOR THE AMIGA
C     ABSOFT FORTRAN COMPILER USING THE 68020/68881 COMPILER OPTION
C
C     DATA SMALL(1), SMALL(2) / Z'00100000', Z'00000000' /
C     DATA LARGE(1), LARGE(2) / Z'7FEFFFFF', Z'FFFFFFFF' /
C     DATA RIGHT(1), RIGHT(2) / Z'3CA00000', Z'00000000' /
C     DATA DIVER(1), DIVER(2) / Z'3CB00000', Z'00000000' /
C     DATA LOG10(1), LOG10(2) / Z'3FD34413', Z'509F79FF' /
C
C     MACHINE CONSTANTS FOR THE AMIGA
C     ABSOFT FORTRAN COMPILER USING SOFTWARE FLOATING POINT
C
C     DATA SMALL(1), SMALL(2) / Z'00100000', Z'00000000' /
C     DATA LARGE(1), LARGE(2) / Z'7FDFFFFF', Z'FFFFFFFF' /
C     DATA RIGHT(1), RIGHT(2) / Z'3CA00000', Z'00000000' /
C     DATA DIVER(1), DIVER(2) / Z'3CB00000', Z'00000000' /
C     DATA LOG10(1), LOG10(2) / Z'3FD34413', Z'509F79FF' /
C
C     MACHINE CONSTANTS FOR THE APOLLO
C
C     DATA SMALL(1), SMALL(2) / 16#00100000, 16#00000000 /
C     DATA LARGE(1), LARGE(2) / 16#7FFFFFFF, 16#FFFFFFFF /
C     DATA RIGHT(1), RIGHT(2) / 16#3CA00000, 16#00000000 /
C     DATA DIVER(1), DIVER(2) / 16#3CB00000, 16#00000000 /
C     DATA LOG10(1), LOG10(2) / 16#3FD34413, 16#509F79FF /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM
C
C     DATA SMALL(1) / ZC00800000 /
C     DATA SMALL(2) / Z000000000 /
C     DATA LARGE(1) / ZDFFFFFFFF /
C     DATA LARGE(2) / ZFFFFFFFFF /
C     DATA RIGHT(1) / ZCC5800000 /
C     DATA RIGHT(2) / Z000000000 /
C     DATA DIVER(1) / ZCC6800000 /
C     DATA DIVER(2) / Z000000000 /
C     DATA LOG10(1) / ZD00E730E7 /
C     DATA LOG10(2) / ZC77800DC0 /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 5700 SYSTEM
C
C     DATA SMALL(1) / O1771000000000000 /
C     DATA SMALL(2) / O0000000000000000 /
C     DATA LARGE(1) / O0777777777777777 /
C     DATA LARGE(2) / O0007777777777777 /
C     DATA RIGHT(1) / O1461000000000000 /
C     DATA RIGHT(2) / O0000000000000000 /
C     DATA DIVER(1) / O1451000000000000 /
C     DATA DIVER(2) / O0000000000000000 /
C     DATA LOG10(1) / O1157163034761674 /
C     DATA LOG10(2) / O0006677466732724 /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 6700/7700 SYSTEMS
C
C     DATA SMALL(1) / O1771000000000000 /
C     DATA SMALL(2) / O7770000000000000 /
C     DATA LARGE(1) / O0777777777777777 /
C     DATA LARGE(2) / O7777777777777777 /
C     DATA RIGHT(1) / O1461000000000000 /
C     DATA RIGHT(2) / O0000000000000000 /
C     DATA DIVER(1) / O1451000000000000 /
C     DATA DIVER(2) / O0000000000000000 /
C     DATA LOG10(1) / O1157163034761674 /
C     DATA LOG10(2) / O0006677466732724 /
C
C     MACHINE CONSTANTS FOR THE CDC 170/180 SERIES USING NOS/VE
C
C     DATA SMALL(1) / Z"3001800000000000" /
C     DATA SMALL(2) / Z"3001000000000000" /
C     DATA LARGE(1) / Z"4FFEFFFFFFFFFFFE" /
C     DATA LARGE(2) / Z"4FFE000000000000" /
C     DATA RIGHT(1) / Z"3FD2800000000000" /
C     DATA RIGHT(2) / Z"3FD2000000000000" /
C     DATA DIVER(1) / Z"3FD3800000000000" /
C     DATA DIVER(2) / Z"3FD3000000000000" /
C     DATA LOG10(1) / Z"3FFF9A209A84FBCF" /
C     DATA LOG10(2) / Z"3FFFF7988F8959AC" /
C
C     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES
C
C     DATA SMALL(1) / 00564000000000000000B /
C     DATA SMALL(2) / 00000000000000000000B /
C     DATA LARGE(1) / 37757777777777777777B /
C     DATA LARGE(2) / 37157777777777777777B /
C     DATA RIGHT(1) / 15624000000000000000B /
C     DATA RIGHT(2) / 00000000000000000000B /
C     DATA DIVER(1) / 15634000000000000000B /
C     DATA DIVER(2) / 00000000000000000000B /
C     DATA LOG10(1) / 17164642023241175717B /
C     DATA LOG10(2) / 16367571421742254654B /
C
C     MACHINE CONSTANTS FOR THE CELERITY C1260
C
C     DATA SMALL(1), SMALL(2) / Z'00100000', Z'00000000' /
C     DATA LARGE(1), LARGE(2) / Z'7FEFFFFF', Z'FFFFFFFF' /
C     DATA RIGHT(1), RIGHT(2) / Z'3CA00000', Z'00000000' /
C     DATA DIVER(1), DIVER(2) / Z'3CB00000', Z'00000000' /
C     DATA LOG10(1), LOG10(2) / Z'3FD34413', Z'509F79FF' /
C
C     MACHINE CONSTANTS FOR THE CONVEX C-1
C
C     DATA SMALL(1), SMALL(2) / '00100000'X,'00000000'X /
C     DATA LARGE(1), LARGE(2) / '7FFFFFFF'X,'FFFFFFFF'X /
C     DATA RIGHT(1), RIGHT(2) / '3CC00000'X,'00000000'X /
C     DATA DIVER(1), DIVER(2) / '3CD00000'X,'00000000'X /
C     DATA LOG10(1), LOG10(2) / '3FF34413'X,'509F79FF'X /
C
C     MACHINE CONSTANTS FOR THE CRAY-1
C
C     DATA SMALL(1) / 201354000000000000000B /
C     DATA SMALL(2) / 000000000000000000000B /
C     DATA LARGE(1) / 577767777777777777777B /
C     DATA LARGE(2) / 000007777777777777774B /
C     DATA RIGHT(1) / 376434000000000000000B /
C     DATA RIGHT(2) / 000000000000000000000B /
C     DATA DIVER(1) / 376444000000000000000B /
C     DATA DIVER(2) / 000000000000000000000B /
C     DATA LOG10(1) / 377774642023241175717B /
C     DATA LOG10(2) / 000007571421742254654B /
C
C     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200
C
C     NOTE - IT MAY BE APPROPRIATE TO INCLUDE THE FOLLOWING CARD -
C     STATIC DMACH(5)
C
C     DATA SMALL /    20K, 3*0 /
C     DATA LARGE / 77777K, 3*177777K /
C     DATA RIGHT / 31420K, 3*0 /
C     DATA DIVER / 32020K, 3*0 /
C     DATA LOG10 / 40423K, 42023K, 50237K, 74776K /
C
C     MACHINE CONSTANTS FOR THE ELXSI 6400
C     (ASSUMING REAL*8 IS THE DEFAULT DOUBLE PRECISION)
C
C     DATA SMALL(1), SMALL(2) / '00100000'X,'00000000'X /
C     DATA LARGE(1), LARGE(2) / '7FEFFFFF'X,'FFFFFFFF'X /
C     DATA RIGHT(1), RIGHT(2) / '3CB00000'X,'00000000'X /
C     DATA DIVER(1), DIVER(2) / '3CC00000'X,'00000000'X /
C     DATA LOG10(1), LOG10(2) / '3FD34413'X,'509F79FF'X /
C
C     MACHINE CONSTANTS FOR THE HARRIS 220
C
C     DATA SMALL(1), SMALL(2) / '20000000, '00000201 /
C     DATA LARGE(1), LARGE(2) / '37777777, '37777577 /
C     DATA RIGHT(1), RIGHT(2) / '20000000, '00000333 /
C     DATA DIVER(1), DIVER(2) / '20000000, '00000334 /
C     DATA LOG10(1), LOG10(2) / '23210115, '10237777 /
C
C     MACHINE CONSTANTS FOR THE HONEYWELL 600/6000 SERIES
C
C     DATA SMALL(1), SMALL(2) / O402400000000, O000000000000 /
C     DATA LARGE(1), LARGE(2) / O376777777777, O777777777777 /
C     DATA RIGHT(1), RIGHT(2) / O604400000000, O000000000000 /
C     DATA DIVER(1), DIVER(2) / O606400000000, O000000000000 /
C     DATA LOG10(1), LOG10(2) / O776464202324, O117571775714 /
C
C     MACHINE CONSTANTS FOR THE HP 2100
C     THREE WORD DOUBLE PRECISION OPTION WITH FTN4
C
C     DATA SMALL(1), SMALL(2), SMALL(3) / 40000B,       0,       1 /
C     DATA LARGE(1), LARGE(2), LARGE(3) / 77777B, 177777B, 177776B /
C     DATA RIGHT(1), RIGHT(2), RIGHT(3) / 40000B,       0,    265B /
C     DATA DIVER(1), DIVER(2), DIVER(3) / 40000B,       0,    276B /
C     DATA LOG10(1), LOG10(2), LOG10(3) / 46420B,  46502B,  77777B /
C
C     MACHINE CONSTANTS FOR THE HP 2100
C     FOUR WORD DOUBLE PRECISION OPTION WITH FTN4
C
C     DATA SMALL(1), SMALL(2) /  40000B,       0 /
C     DATA SMALL(3), SMALL(4) /       0,       1 /
C     DATA LARGE(1), LARGE(2) /  77777B, 177777B /
C     DATA LARGE(3), LARGE(4) / 177777B, 177776B /
C     DATA RIGHT(1), RIGHT(2) /  40000B,       0 /
C     DATA RIGHT(3), RIGHT(4) /       0,    225B /
C     DATA DIVER(1), DIVER(2) /  40000B,       0 /
C     DATA DIVER(3), DIVER(4) /       0,    227B /
C     DATA LOG10(1), LOG10(2) /  46420B,  46502B /
C     DATA LOG10(3), LOG10(4) /  76747B, 176377B /
C
C     MACHINE CONSTANTS FOR THE HP 9000
C
C     DATA SMALL(1), SMALL(2) / 00040000000B, 00000000000B /
C     DATA LARGE(1), LARGE(2) / 17737777777B, 37777777777B /
C     DATA RIGHT(1), RIGHT(2) / 07454000000B, 00000000000B /
C     DATA DIVER(1), DIVER(2) / 07460000000B, 00000000000B /
C     DATA LOG10(1), LOG10(2) / 07764642023B, 12047674777B /
C
C     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
C     THE XEROX SIGMA 5/7/9, THE SEL SYSTEMS 85/86, AND
C     THE PERKIN ELMER (INTERDATA) 7/32.
C
C     DATA SMALL(1), SMALL(2) / Z00100000, Z00000000 /
C     DATA LARGE(1), LARGE(2) / Z7FFFFFFF, ZFFFFFFFF /
C     DATA RIGHT(1), RIGHT(2) / Z33100000, Z00000000 /
C     DATA DIVER(1), DIVER(2) / Z34100000, Z00000000 /
C     DATA LOG10(1), LOG10(2) / Z41134413, Z509F79FF /
C
C     MACHINE CONSTANTS FOR THE IBM PC
C     ASSUMES THAT ALL ARITHMETIC IS DONE IN DOUBLE PRECISION
C     ON 8088, I.E., NOT IN 80 BIT FORM FOR THE 8087.
C
C     DATA SMALL(1),SMALL(2) /  2002288515,    1050897 /
C     DATA LARGE(1),LARGE(2) /  1487780761, 2146426097 /
C     DATA RIGHT(1),RIGHT(2) / -1209488034, 1017118298 /
C     DATA DIVER(1),DIVER(2) / -1209488034, 1018166874 /
C     DATA LOG10(1),LOG10(2) /  1352628735, 1070810131 /
C
C     MACHINE CONSTANTS FOR THE IBM RS 6000         ! For SOLO
C
cc     DATA SMALL(1), SMALL(2) / Z'00100000', Z'00000000' /
cc     DATA LARGE(1), LARGE(2) / Z'7FEFFFFF', Z'FFFFFFFF' /
cc     DATA RIGHT(1), RIGHT(2) / Z'3CA00000', Z'00000000' /
cc     DATA DIVER(1), DIVER(2) / Z'3CB00000', Z'00000000' /
cc     DATA LOG10(1), LOG10(2) / Z'3FD34413', Z'509F79FF' /
C
C     MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR)
C
C     DATA SMALL(1), SMALL(2) / "033400000000, "000000000000 /
C     DATA LARGE(1), LARGE(2) / "377777777777, "344777777777 /
C     DATA RIGHT(1), RIGHT(2) / "113400000000, "000000000000 /
C     DATA DIVER(1), DIVER(2) / "114400000000, "000000000000 /
C     DATA LOG10(1), LOG10(2) / "177464202324, "144117571776 /
C
C     MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR)
C
C     DATA SMALL(1), SMALL(2) / "000400000000, "000000000000 /
C     DATA LARGE(1), LARGE(2) / "377777777777, "377777777777 /
C     DATA RIGHT(1), RIGHT(2) / "103400000000, "000000000000 /
C     DATA DIVER(1), DIVER(2) / "104400000000, "000000000000 /
C     DATA LOG10(1), LOG10(2) / "177464202324, "476747767461 /
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
C     32-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL).
C
C     DATA SMALL(1), SMALL(2) /    8388608,           0 /
C     DATA LARGE(1), LARGE(2) / 2147483647,          -1 /
C     DATA RIGHT(1), RIGHT(2) /  612368384,           0 /
C     DATA DIVER(1), DIVER(2) /  620756992,           0 /
C     DATA LOG10(1), LOG10(2) / 1067065498, -2063872008 /
C
C     DATA SMALL(1), SMALL(2) / O00040000000, O00000000000 /
C     DATA LARGE(1), LARGE(2) / O17777777777, O37777777777 /
C     DATA RIGHT(1), RIGHT(2) / O04440000000, O00000000000 /
C     DATA DIVER(1), DIVER(2) / O04500000000, O00000000000 /
C     DATA LOG10(1), LOG10(2) / O07746420232, O20476747770 /
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
C     16-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL).
C
C     DATA SMALL(1), SMALL(2) /    128,      0 /
C     DATA SMALL(3), SMALL(4) /      0,      0 /
C     DATA LARGE(1), LARGE(2) /  32767,     -1 /
C     DATA LARGE(3), LARGE(4) /     -1,     -1 /
C     DATA RIGHT(1), RIGHT(2) /   9344,      0 /
C     DATA RIGHT(3), RIGHT(4) /      0,      0 /
C     DATA DIVER(1), DIVER(2) /   9472,      0 /
C     DATA DIVER(3), DIVER(4) /      0,      0 /
C     DATA LOG10(1), LOG10(2) /  16282,   8346 /
C     DATA LOG10(3), LOG10(4) / -31493, -12296 /
C
C     DATA SMALL(1), SMALL(2) / O000200, O000000 /
C     DATA SMALL(3), SMALL(4) / O000000, O000000 /
C     DATA LARGE(1), LARGE(2) / O077777, O177777 /
C     DATA LARGE(3), LARGE(4) / O177777, O177777 /
C     DATA RIGHT(1), RIGHT(2) / O022200, O000000 /
C     DATA RIGHT(3), RIGHT(4) / O000000, O000000 /
C     DATA DIVER(1), DIVER(2) / O022400, O000000 /
C     DATA DIVER(3), DIVER(4) / O000000, O000000 /
C     DATA LOG10(1), LOG10(2) / O037632, O020232 /
C     DATA LOG10(3), LOG10(4) / O102373, O147770 /
C
C     MACHINE CONSTANTS FOR THE SILICON GRAPHICS IRIS
C
C >>> ruunat:
      data dmach(1) / 2.22507 38585 072013 d-308 /
      data dmach(2) / 1.79769 31348 623158 d+308 /
      data dmach(3) / 2.22044 60492 503131 d-16  /
      data dmach(4) / 4.44089 20985 006262 d-16  /
      data dmach(5) / 0.30102 99956 639812  d0   /
C >>> SARA:
*      data dmach(1) / 2.22507 38585 0720138 d-308 /
*      data dmach(2) / 1.79769 31348 6231571 d+308 /
*      data dmach(3) / 2.22044 60492 5031308 d-16  /
*      data dmach(4) / 4.44089 20985 0062616 d-16  /
*      data dmach(5) / 0.30102 99956 63981198   d0 /
C
C     DATA SMALL(1), SMALL(2) / Z'00100000',Z'00000000' /
C     DATA LARGE(1), LARGE(2) / Z'7FEFFFFF',Z'FFFFFFFF' /
C     DATA RIGHT(1), RIGHT(2) / Z'3CB00000',Z'00000000' /
C     DATA DIVER(1), DIVER(2) / Z'3CC00000',Z'00000000' /
C     DATA LOG10(1), LOG10(2) / Z'3FD34413',Z'509F79FF' /
C
C     MACHINE CONSTANTS FOR THE SUN
C
C     from SLATEC CML committee - work for Sun3, Sun4, and Sparc
C
C     DATA SMALL(1), SMALL(2) / Z'00100000', Z'00000000' /
C     DATA LARGE(1), LARGE(2) / Z'7FEFFFFF', Z'FFFFFFFF' /
C     DATA RIGHT(1), RIGHT(2) / Z'3CA00000', Z'00000000' /
C     DATA DIVER(1), DIVER(2) / Z'3CB00000', Z'00000000' /
C     DATA LOG10(1), LOG10(2) / Z'3FD34413', Z'509F79FF' /
C
C     from Sun Microsystems - work for Sun 386i
C
C     DATA SMALL(1), SMALL(2) / Z'00000000', Z'00100000' /
C     DATA LARGE(1), LARGE(2) / Z'FFFFFFFF', Z'7FEFFFFF' /
C     DATA RIGHT(1), RIGHT(2) / Z'00000000', Z'3CA00000' /
C     DATA DIVER(1), DIVER(2) / Z'00000000', Z'3CB00000' /
C     DATA LOG10(1), LOG10(2) / Z'509F79FF', Z'3FD34413' /
C
C     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES FTN COMPILER
C
C     DATA SMALL(1), SMALL(2) / O000040000000, O000000000000 /
C     DATA LARGE(1), LARGE(2) / O377777777777, O777777777777 /
C     DATA RIGHT(1), RIGHT(2) / O170540000000, O000000000000 /
C     DATA DIVER(1), DIVER(2) / O170640000000, O000000000000 /
C     DATA LOG10(1), LOG10(2) / O177746420232, O411757177572 /
C
C     MACHINE CONSTANTS FOR VAX 11/780
C     (EXPRESSED IN INTEGER AND HEXADECIMAL)
C     THE HEX FORMAT BELOW MAY NOT BE SUITABLE FOR UNIX SYSYEMS
C     THE INTEGER FORMAT SHOULD BE OK FOR UNIX SYSTEMS
C
C     DATA SMALL(1), SMALL(2) /        128,           0 /
C     DATA LARGE(1), LARGE(2) /     -32769,          -1 /
C     DATA RIGHT(1), RIGHT(2) /       9344,           0 /
C     DATA DIVER(1), DIVER(2) /       9472,           0 /
C     DATA LOG10(1), LOG10(2) /  546979738,  -805796613 /
C
C     DATA SMALL(1), SMALL(2) / Z00000080, Z00000000 /
C     DATA LARGE(1), LARGE(2) / ZFFFF7FFF, ZFFFFFFFF /
C     DATA RIGHT(1), RIGHT(2) / Z00002480, Z00000000 /
C     DATA DIVER(1), DIVER(2) / Z00002500, Z00000000 /
C     DATA LOG10(1), LOG10(2) / Z209A3F9A, ZCFF884FB /
C
C     MACHINE CONSTANTS FOR VAX 11/780 (G-FLOATING)
C     (EXPRESSED IN INTEGER AND HEXADECIMAL)
C     THE HEX FORMAT BELOW MAY NOT BE SUITABLE FOR UNIX SYSYEMS
C     THE INTEGER FORMAT SHOULD BE OK FOR UNIX SYSTEMS
C
C     DATA SMALL(1), SMALL(2) /         16,           0 /
C     DATA LARGE(1), LARGE(2) /     -32769,          -1 /
C     DATA RIGHT(1), RIGHT(2) /      15552,           0 /
C     DATA DIVER(1), DIVER(2) /      15568,           0 /
C     DATA LOG10(1), LOG10(2) /  1142112243, 2046775455 /
C
C     DATA SMALL(1), SMALL(2) / Z00000010, Z00000000 /
C     DATA LARGE(1), LARGE(2) / ZFFFF7FFF, ZFFFFFFFF /
C     DATA RIGHT(1), RIGHT(2) / Z00003CC0, Z00000000 /
C     DATA DIVER(1), DIVER(2) / Z00003CD0, Z00000000 /
C     DATA LOG10(1), LOG10(2) / Z44133FF3, Z79FF509F /
C
C
C***FIRST EXECUTABLE STATEMENT  D1MACH
      IF (I .LT. 1  .OR.  I .GT. 5)
     1   CALL XERROR ('D1MACH -- I OUT OF BOUNDS', 25, 1, 2)
C
      D1MACH = DMACH(I)
      RETURN
C
      END
      DOUBLE PRECISION FUNCTION DGAMLN(Z,IERR)
C***BEGIN PROLOGUE  DGAMLN
C***DATE WRITTEN   830501   (YYMMDD)
C***REVISION DATE  830501   (YYMMDD)
C***CATEGORY NO.  B5F
C***KEYWORDS  GAMMA FUNCTION,LOGARITHM OF GAMMA FUNCTION
C***AUTHOR  AMOS, DONALD E., SANDIA NATIONAL LABORATORIES
C***PURPOSE  TO COMPUTE THE LOGARITHM OF THE GAMMA FUNCTION
C***DESCRIPTION
C
C               **** A DOUBLE PRECISION ROUTINE ****
C         DGAMLN COMPUTES THE NATURAL LOG OF THE GAMMA FUNCTION FOR
C         Z.GT.0.  THE ASYMPTOTIC EXPANSION IS USED TO GENERATE VALUES
C         GREATER THAN ZMIN WHICH ARE ADJUSTED BY THE RECURSION
C         G(Z+1)=Z*G(Z) FOR Z.LE.ZMIN.  THE FUNCTION WAS MADE AS
C         PORTABLE AS POSSIBLE BY COMPUTIMG ZMIN FROM THE NUMBER OF BASE
C         10 DIGITS IN A WORD, RLN=AMAX1(-ALOG10(R1MACH(4)),0.5E-18)
C         LIMITED TO 18 DIGITS OF (RELATIVE) ACCURACY.
C
C         SINCE INTEGER ARGUMENTS ARE COMMON, A TABLE LOOK UP ON 100
C         VALUES IS USED FOR SPEED OF EXECUTION.
C
C     DESCRIPTION OF ARGUMENTS
C
C         INPUT      Z IS D0UBLE PRECISION
C           Z      - ARGUMENT, Z.GT.0.0D0
C
C         OUTPUT      DGAMLN IS DOUBLE PRECISION
C           DGAMLN  - NATURAL LOG OF THE GAMMA FUNCTION AT Z.NE.0.0D0
C           IERR    - ERROR FLAG
C                     IERR=0, NORMAL RETURN, COMPUTATION COMPLETED
C                     IERR=1, Z.LE.0.0D0,    NO COMPUTATION
C
C
C***REFERENCES  COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
C                 BY D. E. AMOS, SAND83-0083, MAY, 1983.
C***ROUTINES CALLED  I1MACH,D1MACH
C***END PROLOGUE  DGAMLN
      DOUBLE PRECISION CF, CON, FLN, FZ, GLN, RLN, S, TLG, TRM, TST,
     * T1, WDTOL, Z, ZDMY, ZINC, ZM, ZMIN, ZP, ZSQ, D1MACH
      INTEGER I, IERR, I1M, K, MZ, NZ, I1MACH
      DIMENSION CF(22), GLN(100)
C           LNGAMMA(N), N=1,100
      DATA GLN(1), GLN(2), GLN(3), GLN(4), GLN(5), GLN(6), GLN(7),
     1     GLN(8), GLN(9), GLN(10), GLN(11), GLN(12), GLN(13), GLN(14),
     2     GLN(15), GLN(16), GLN(17), GLN(18), GLN(19), GLN(20),
     3     GLN(21), GLN(22)/
     4     0.00000000000000000D+00,     0.00000000000000000D+00,
     5     6.93147180559945309D-01,     1.79175946922805500D+00,
     6     3.17805383034794562D+00,     4.78749174278204599D+00,
     7     6.57925121201010100D+00,     8.52516136106541430D+00,
     8     1.06046029027452502D+01,     1.28018274800814696D+01,
     9     1.51044125730755153D+01,     1.75023078458738858D+01,
     A     1.99872144956618861D+01,     2.25521638531234229D+01,
     B     2.51912211827386815D+01,     2.78992713838408916D+01,
     C     3.06718601060806728D+01,     3.35050734501368889D+01,
     D     3.63954452080330536D+01,     3.93398841871994940D+01,
     E     4.23356164607534850D+01,     4.53801388984769080D+01/
      DATA GLN(23), GLN(24), GLN(25), GLN(26), GLN(27), GLN(28),
     1     GLN(29), GLN(30), GLN(31), GLN(32), GLN(33), GLN(34),
     2     GLN(35), GLN(36), GLN(37), GLN(38), GLN(39), GLN(40),
     3     GLN(41), GLN(42), GLN(43), GLN(44)/
     4     4.84711813518352239D+01,     5.16066755677643736D+01,
     5     5.47847293981123192D+01,     5.80036052229805199D+01,
     6     6.12617017610020020D+01,     6.45575386270063311D+01,
     7     6.78897431371815350D+01,     7.12570389671680090D+01,
     8     7.46582363488301644D+01,     7.80922235533153106D+01,
     9     8.15579594561150372D+01,     8.50544670175815174D+01,
     A     8.85808275421976788D+01,     9.21361756036870925D+01,
     B     9.57196945421432025D+01,     9.93306124547874269D+01,
     C     1.02968198614513813D+02,     1.06631760260643459D+02,
     D     1.10320639714757395D+02,     1.14034211781461703D+02,
     E     1.17771881399745072D+02,     1.21533081515438634D+02/
      DATA GLN(45), GLN(46), GLN(47), GLN(48), GLN(49), GLN(50),
     1     GLN(51), GLN(52), GLN(53), GLN(54), GLN(55), GLN(56),
     2     GLN(57), GLN(58), GLN(59), GLN(60), GLN(61), GLN(62),
     3     GLN(63), GLN(64), GLN(65), GLN(66)/
     4     1.25317271149356895D+02,     1.29123933639127215D+02,
     5     1.32952575035616310D+02,     1.36802722637326368D+02,
     6     1.40673923648234259D+02,     1.44565743946344886D+02,
     7     1.48477766951773032D+02,     1.52409592584497358D+02,
     8     1.56360836303078785D+02,     1.60331128216630907D+02,
     9     1.64320112263195181D+02,     1.68327445448427652D+02,
     A     1.72352797139162802D+02,     1.76395848406997352D+02,
     B     1.80456291417543771D+02,     1.84533828861449491D+02,
     C     1.88628173423671591D+02,     1.92739047287844902D+02,
     D     1.96866181672889994D+02,     2.01009316399281527D+02,
     E     2.05168199482641199D+02,     2.09342586752536836D+02/
      DATA GLN(67), GLN(68), GLN(69), GLN(70), GLN(71), GLN(72),
     1     GLN(73), GLN(74), GLN(75), GLN(76), GLN(77), GLN(78),
     2     GLN(79), GLN(80), GLN(81), GLN(82), GLN(83), GLN(84),
     3     GLN(85), GLN(86), GLN(87), GLN(88)/
     4     2.13532241494563261D+02,     2.17736934113954227D+02,
     5     2.21956441819130334D+02,     2.26190548323727593D+02,
     6     2.30439043565776952D+02,     2.34701723442818268D+02,
     7     2.38978389561834323D+02,     2.43268849002982714D+02,
     8     2.47572914096186884D+02,     2.51890402209723194D+02,
     9     2.56221135550009525D+02,     2.60564940971863209D+02,
     A     2.64921649798552801D+02,     2.69291097651019823D+02,
     B     2.73673124285693704D+02,     2.78067573440366143D+02,
     C     2.82474292687630396D+02,     2.86893133295426994D+02,
     D     2.91323950094270308D+02,     2.95766601350760624D+02,
     E     3.00220948647014132D+02,     3.04686856765668715D+02/
      DATA GLN(89), GLN(90), GLN(91), GLN(92), GLN(93), GLN(94),
     1     GLN(95), GLN(96), GLN(97), GLN(98), GLN(99), GLN(100)/
     2     3.09164193580146922D+02,     3.13652829949879062D+02,
     3     3.18152639620209327D+02,     3.22663499126726177D+02,
     4     3.27185287703775217D+02,     3.31717887196928473D+02,
     5     3.36261181979198477D+02,     3.40815058870799018D+02,
     6     3.45379407062266854D+02,     3.49954118040770237D+02,
     7     3.54539085519440809D+02,     3.59134205369575399D+02/
C             COEFFICIENTS OF ASYMPTOTIC EXPANSION
      DATA CF(1), CF(2), CF(3), CF(4), CF(5), CF(6), CF(7), CF(8),
     1     CF(9), CF(10), CF(11), CF(12), CF(13), CF(14), CF(15),
     2     CF(16), CF(17), CF(18), CF(19), CF(20), CF(21), CF(22)/
     3     8.33333333333333333D-02,    -2.77777777777777778D-03,
     4     7.93650793650793651D-04,    -5.95238095238095238D-04,
     5     8.41750841750841751D-04,    -1.91752691752691753D-03,
     6     6.41025641025641026D-03,    -2.95506535947712418D-02,
     7     1.79644372368830573D-01,    -1.39243221690590112D+00,
     8     1.34028640441683920D+01,    -1.56848284626002017D+02,
     9     2.19310333333333333D+03,    -3.61087712537249894D+04,
     A     6.91472268851313067D+05,    -1.52382215394074162D+07,
     B     3.82900751391414141D+08,    -1.08822660357843911D+10,
     C     3.47320283765002252D+11,    -1.23696021422692745D+13,
     D     4.88788064793079335D+14,    -2.13203339609193739D+16/
C
C             LN(2*PI)
      DATA CON                    /     1.83787706640934548D+00/
C
C***FIRST EXECUTABLE STATEMENT  DGAMLN
      IERR=0
      IF (Z.LE.0.0D0) GO TO 70
      IF (Z.GT.101.0D0) GO TO 10
      NZ = INT(SNGL(Z))
      FZ = Z - FLOAT(NZ)
      IF (FZ.GT.0.0D0) GO TO 10
      IF (NZ.GT.100) GO TO 10
      DGAMLN = GLN(NZ)
      RETURN
   10 CONTINUE
      WDTOL = D1MACH(4)
      WDTOL = DMAX1(WDTOL,0.5D-18)
      I1M = I1MACH(14)
      RLN = D1MACH(5)*FLOAT(I1M)
      FLN = DMIN1(RLN,20.0D0)
      FLN = DMAX1(FLN,3.0D0)
      FLN = FLN - 3.0D0
      ZM = 1.8000D0 + 0.3875D0*FLN
      MZ = INT(SNGL(ZM)) + 1
      ZMIN = FLOAT(MZ)
      ZDMY = Z
      ZINC = 0.0D0
      IF (Z.GE.ZMIN) GO TO 20
      ZINC = ZMIN - FLOAT(NZ)
      ZDMY = Z + ZINC
   20 CONTINUE
      ZP = 1.0D0/ZDMY
      T1 = CF(1)*ZP
      S = T1
      IF (ZP.LT.WDTOL) GO TO 40
      ZSQ = ZP*ZP
      TST = T1*WDTOL
      DO 30 K=2,22
        ZP = ZP*ZSQ
        TRM = CF(K)*ZP
        IF (DABS(TRM).LT.TST) GO TO 40
        S = S + TRM
   30 CONTINUE
   40 CONTINUE
      IF (ZINC.NE.0.0D0) GO TO 50
      TLG = DLOG(Z)
      DGAMLN = Z*(TLG-1.0D0) + 0.5D0*(CON-TLG) + S
      RETURN
   50 CONTINUE
      ZP = 1.0D0
      NZ = INT(SNGL(ZINC))
      DO 60 I=1,NZ
        ZP = ZP*(Z+FLOAT(I-1))
   60 CONTINUE
      TLG = DLOG(ZDMY)
      DGAMLN = ZDMY*(TLG-1.0D0) - DLOG(ZP) + 0.5D0*(CON-TLG) + S
      RETURN
C
C
   70 CONTINUE
      IERR=1
      RETURN
      END
*DECK I1MACH
      INTEGER FUNCTION I1MACH(I)
C***BEGIN PROLOGUE  I1MACH
C***DATE WRITTEN   750101   (YYMMDD)
C***REVISION DATE  890213   (YYMMDD)
C***CATEGORY NO.  R1
C***KEYWORDS  LIBRARY=SLATEC,TYPE=INTEGER(I1MACH-I),MACHINE CONSTANTS
C***AUTHOR  FOX, P. A., (BELL LABS)
C           HALL, A. D., (BELL LABS)
C           SCHRYER, N. L., (BELL LABS)
C***PURPOSE  Returns integer machine dependent constants
C***DESCRIPTION
C
C     I1MACH can be used to obtain machine-dependent parameters
C     for the local machine environment.  It is a function
C     subroutine with one (input) argument, and can be called
C     as follows, for example
C
C          K = I1MACH(I)
C
C     where I=1,...,16.  The (output) value of K above is
C     determined by the (input) value of I.  The results for
C     various values of I are discussed below.
C
C  I/O unit numbers.
C    I1MACH( 1) = the standard input unit.
C    I1MACH( 2) = the standard output unit.
C    I1MACH( 3) = the standard punch unit.
C    I1MACH( 4) = the standard error message unit.
C
C  Words.
C    I1MACH( 5) = the number of bits per integer storage unit.
C    I1MACH( 6) = the number of characters per integer storage unit.
C
C  Integers.
C    assume integers are represented in the S-digit, base-A form
C
C               sign ( X(S-1)*A**(S-1) + ... + X(1)*A + X(0) )
C
C               where 0 .LE. X(I) .LT. A for I=0,...,S-1.
C    I1MACH( 7) = A, the base.
C    I1MACH( 8) = S, the number of base-A digits.
C    I1MACH( 9) = A**S - 1, the largest magnitude.
C
C  Floating-Point Numbers.
C    Assume floating-point numbers are represented in the T-digit,
C    base-B form
C               sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
C
C               where 0 .LE. X(I) .LT. B for I=1,...,T,
C               0 .LT. X(1), and EMIN .LE. E .LE. EMAX.
C    I1MACH(10) = B, the base.
C
C  Single-Precision
C    I1MACH(11) = T, the number of base-B digits.
C    I1MACH(12) = EMIN, the smallest exponent E.
C    I1MACH(13) = EMAX, the largest exponent E.
C
C  Double-Precision
C    I1MACH(14) = T, the number of base-B digits.
C    I1MACH(15) = EMIN, the smallest exponent E.
C    I1MACH(16) = EMAX, the largest exponent E.
C
C  To alter this function for a particular environment,
C  the desired set of DATA statements should be activated by
C  removing the C from column 1.  Also, the values of
C  I1MACH(1) - I1MACH(4) should be checked for consistency
C  with the local operating system.
C
C***REFERENCES  FOX P.A., HALL A.D., SCHRYER N.L.,*FRAMEWORK FOR A
C                 PORTABLE LIBRARY*, ACM TRANSACTIONS ON MATHEMATICAL
C                 SOFTWARE, VOL. 4, NO. 2, JUNE 1978, PP. 177-188.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  I1MACH
C
      INTEGER IMACH(16),OUTPUT
      SAVE IMACH
      EQUIVALENCE (IMACH(4),OUTPUT)
C
C     MACHINE CONSTANTS FOR THE AMIGA
C     ABSOFT COMPILER
C
C     DATA IMACH(1) /    5 /
C     DATA IMACH(2) /    6 /
C     DATA IMACH(3) /    5 /
C     DATA IMACH(4) /    6 /
C     DATA IMACH(5) /   32 /
C     DATA IMACH(6) /    4 /
C     DATA IMACH(7) /    2 /
C     DATA IMACH(8) /   31 /
C     DATA IMACH(9) / 2147483647 /
C     DATA IMACH(10)/    2 /
C     DATA IMACH(11)/   24 /
C     DATA IMACH(12)/ -126 /
C     DATA IMACH(13)/  127 /
C     DATA IMACH(14)/   53 /
C     DATA IMACH(15)/ -1022 /
C     DATA IMACH(16)/  1023 /
C
C     MACHINE CONSTANTS FOR THE APOLLO
C
C     DATA IMACH(1) /    5 /
C     DATA IMACH(2) /    6 /
C     DATA IMACH(3) /    6 /
C     DATA IMACH(4) /    6 /
C     DATA IMACH(5) /   32 /
C     DATA IMACH(6) /    4 /
C     DATA IMACH(7) /    2 /
C     DATA IMACH(8) /   31 /
C     DATA IMACH(9) / 2147483647 /
C     DATA IMACH(10)/    2 /
C     DATA IMACH(11)/   24 /
C     DATA IMACH(12)/ -125 /
C     DATA IMACH(13)/  129 /
C     DATA IMACH(14)/   53 /
C     DATA IMACH(15)/ -1021 /
C     DATA IMACH(16)/  1025 /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM
C
C     DATA IMACH( 1) /    7 /
C     DATA IMACH( 2) /    2 /
C     DATA IMACH( 3) /    2 /
C     DATA IMACH( 4) /    2 /
C     DATA IMACH( 5) /   36 /
C     DATA IMACH( 6) /    4 /
C     DATA IMACH( 7) /    2 /
C     DATA IMACH( 8) /   33 /
C     DATA IMACH( 9) / Z1FFFFFFFF /
C     DATA IMACH(10) /    2 /
C     DATA IMACH(11) /   24 /
C     DATA IMACH(12) / -256 /
C     DATA IMACH(13) /  255 /
C     DATA IMACH(14) /   60 /
C     DATA IMACH(15) / -256 /
C     DATA IMACH(16) /  255 /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 5700 SYSTEM
C
C     DATA IMACH( 1) /   5 /
C     DATA IMACH( 2) /   6 /
C     DATA IMACH( 3) /   7 /
C     DATA IMACH( 4) /   6 /
C     DATA IMACH( 5) /  48 /
C     DATA IMACH( 6) /   6 /
C     DATA IMACH( 7) /   2 /
C     DATA IMACH( 8) /  39 /
C     DATA IMACH( 9) / O0007777777777777 /
C     DATA IMACH(10) /   8 /
C     DATA IMACH(11) /  13 /
C     DATA IMACH(12) / -50 /
C     DATA IMACH(13) /  76 /
C     DATA IMACH(14) /  26 /
C     DATA IMACH(15) / -50 /
C     DATA IMACH(16) /  76 /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 6700/7700 SYSTEMS
C
C     DATA IMACH( 1) /   5 /
C     DATA IMACH( 2) /   6 /
C     DATA IMACH( 3) /   7 /
C     DATA IMACH( 4) /   6 /
C     DATA IMACH( 5) /  48 /
C     DATA IMACH( 6) /   6 /
C     DATA IMACH( 7) /   2 /
C     DATA IMACH( 8) /  39 /
C     DATA IMACH( 9) / O0007777777777777 /
C     DATA IMACH(10) /   8 /
C     DATA IMACH(11) /  13 /
C     DATA IMACH(12) / -50 /
C     DATA IMACH(13) /  76 /
C     DATA IMACH(14) /  26 /
C     DATA IMACH(15) / -32754 /
C     DATA IMACH(16) /  32780 /
C
C     MACHINE CONSTANTS FOR THE CDC 170/180 SERIES USING NOS/VE
C
C     DATA IMACH( 1) /     5 /
C     DATA IMACH( 2) /     6 /
C     DATA IMACH( 3) /     7 /
C     DATA IMACH( 4) /     6 /
C     DATA IMACH( 5) /    64 /
C     DATA IMACH( 6) /     8 /
C     DATA IMACH( 7) /     2 /
C     DATA IMACH( 8) /    63 /
C     DATA IMACH( 9) / 9223372036854775807 /
C     DATA IMACH(10) /     2 /
C     DATA IMACH(11) /    47 /
C     DATA IMACH(12) / -4095 /
C     DATA IMACH(13) /  4094 /
C     DATA IMACH(14) /    94 /
C     DATA IMACH(15) / -4095 /
C     DATA IMACH(16) /  4094 /
C
C     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES
C
C     DATA IMACH( 1) /    5 /
C     DATA IMACH( 2) /    6 /
C     DATA IMACH( 3) /    7 /
C     DATA IMACH( 4) /6LOUTPUT/
C     DATA IMACH( 5) /   60 /
C     DATA IMACH( 6) /   10 /
C     DATA IMACH( 7) /    2 /
C     DATA IMACH( 8) /   48 /
C     DATA IMACH( 9) / 00007777777777777777B /
C     DATA IMACH(10) /    2 /
C     DATA IMACH(11) /   47 /
C     DATA IMACH(12) / -929 /
C     DATA IMACH(13) / 1070 /
C     DATA IMACH(14) /   94 /
C     DATA IMACH(15) / -929 /
C     DATA IMACH(16) / 1069 /
C
C     MACHINE CONSTANTS FOR THE CELERITY C1260
C
C     DATA IMACH(1) /    5 /
C     DATA IMACH(2) /    6 /
C     DATA IMACH(3) /    6 /
C     DATA IMACH(4) /    0 /
C     DATA IMACH(5) /   32 /
C     DATA IMACH(6) /    4 /
C     DATA IMACH(7) /    2 /
C     DATA IMACH(8) /   31 /
C     DATA IMACH(9) / Z'7FFFFFFF' /
C     DATA IMACH(10)/    2 /
C     DATA IMACH(11)/   24 /
C     DATA IMACH(12)/ -126 /
C     DATA IMACH(13)/  127 /
C     DATA IMACH(14)/   53 /
C     DATA IMACH(15)/ -1022 /
C     DATA IMACH(16)/  1023 /
C
C     MACHINE CONSTANTS FOR THE CONVEX C-1
C
C     DATA IMACH( 1) /     5/
C     DATA IMACH( 2) /     6/
C     DATA IMACH( 3) /     7/
C     DATA IMACH( 4) /     6/
C     DATA IMACH( 5) /    32/
C     DATA IMACH( 6) /     4/
C     DATA IMACH( 7) /     2/
C     DATA IMACH( 8) /    31/
C     DATA IMACH( 9) /2147483647/
C     DATA IMACH(10) /     2/
C     DATA IMACH(11) /    24/
C     DATA IMACH(12) /  -128/
C     DATA IMACH(13) /   127/
C     DATA IMACH(14) /    53/
C     DATA IMACH(15) / -1024/
C     DATA IMACH(16) /  1023/
C
C     MACHINE CONSTANTS FOR THE CRAY-1
C
C     DATA IMACH( 1) /   100 /
C     DATA IMACH( 2) /   101 /
C     DATA IMACH( 3) /   102 /
C     DATA IMACH( 4) /   101 /
C     DATA IMACH( 5) /    64 /
C     DATA IMACH( 6) /     8 /
C     DATA IMACH( 7) /     2 /
C     DATA IMACH( 8) /    63 /
C     DATA IMACH( 9) /  777777777777777777777B /
C     DATA IMACH(10) /     2 /
C     DATA IMACH(11) /    47 /
C     DATA IMACH(12) / -8189 /
C     DATA IMACH(13) /  8190 /
C     DATA IMACH(14) /    94 /
C     DATA IMACH(15) / -8099 /
C     DATA IMACH(16) /  8190 /
C
C     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200
C
C     DATA IMACH( 1) /   11 /
C     DATA IMACH( 2) /   12 /
C     DATA IMACH( 3) /    8 /
C     DATA IMACH( 4) /   10 /
C     DATA IMACH( 5) /   16 /
C     DATA IMACH( 6) /    2 /
C     DATA IMACH( 7) /    2 /
C     DATA IMACH( 8) /   15 /
C     DATA IMACH( 9) /32767 /
C     DATA IMACH(10) /   16 /
C     DATA IMACH(11) /    6 /
C     DATA IMACH(12) /  -64 /
C     DATA IMACH(13) /   63 /
C     DATA IMACH(14) /   14 /
C     DATA IMACH(15) /  -64 /
C     DATA IMACH(16) /   63 /
C
C     MACHINE CONSTANTS FOR THE ELXSI 6400
C
C     DATA IMACH( 1) /     5/
C     DATA IMACH( 2) /     6/
C     DATA IMACH( 3) /     6/
C     DATA IMACH( 4) /     6/
C     DATA IMACH( 5) /    32/
C     DATA IMACH( 6) /     4/
C     DATA IMACH( 7) /     2/
C     DATA IMACH( 8) /    32/
C     DATA IMACH( 9) /2147483647/
C     DATA IMACH(10) /     2/
C     DATA IMACH(11) /    24/
C     DATA IMACH(12) /  -126/
C     DATA IMACH(13) /   127/
C     DATA IMACH(14) /    53/
C     DATA IMACH(15) / -1022/
C     DATA IMACH(16) /  1023/
C
C     MACHINE CONSTANTS FOR THE HARRIS 220
C
C     DATA IMACH( 1) /       5 /
C     DATA IMACH( 2) /       6 /
C     DATA IMACH( 3) /       0 /
C     DATA IMACH( 4) /       6 /
C     DATA IMACH( 5) /      24 /
C     DATA IMACH( 6) /       3 /
C     DATA IMACH( 7) /       2 /
C     DATA IMACH( 8) /      23 /
C     DATA IMACH( 9) / 8388607 /
C     DATA IMACH(10) /       2 /
C     DATA IMACH(11) /      23 /
C     DATA IMACH(12) /    -127 /
C     DATA IMACH(13) /     127 /
C     DATA IMACH(14) /      38 /
C     DATA IMACH(15) /    -127 /
C     DATA IMACH(16) /     127 /
C
C     MACHINE CONSTANTS FOR THE HONEYWELL 600/6000 SERIES
C
C     DATA IMACH( 1) /    5 /
C     DATA IMACH( 2) /    6 /
C     DATA IMACH( 3) /   43 /
C     DATA IMACH( 4) /    6 /
C     DATA IMACH( 5) /   36 /
C     DATA IMACH( 6) /    6 /
C     DATA IMACH( 7) /    2 /
C     DATA IMACH( 8) /   35 /
C     DATA IMACH( 9) / O377777777777 /
C     DATA IMACH(10) /    2 /
C     DATA IMACH(11) /   27 /
C     DATA IMACH(12) / -127 /
C     DATA IMACH(13) /  127 /
C     DATA IMACH(14) /   63 /
C     DATA IMACH(15) / -127 /
C     DATA IMACH(16) /  127 /
C
C     MACHINE CONSTANTS FOR THE HP 2100
C     3 WORD DOUBLE PRECISION OPTION WITH FTN4
C
C     DATA IMACH(1) /      5/
C     DATA IMACH(2) /      6 /
C     DATA IMACH(3) /      4 /
C     DATA IMACH(4) /      1 /
C     DATA IMACH(5) /     16 /
C     DATA IMACH(6) /      2 /
C     DATA IMACH(7) /      2 /
C     DATA IMACH(8) /     15 /
C     DATA IMACH(9) /  32767 /
C     DATA IMACH(10)/      2 /
C     DATA IMACH(11)/     23 /
C     DATA IMACH(12)/   -128 /
C     DATA IMACH(13)/    127 /
C     DATA IMACH(14)/     39 /
C     DATA IMACH(15)/   -128 /
C     DATA IMACH(16)/    127 /
C
C     MACHINE CONSTANTS FOR THE HP 2100
C     4 WORD DOUBLE PRECISION OPTION WITH FTN4
C
C     DATA IMACH(1) /      5 /
C     DATA IMACH(2) /      6 /
C     DATA IMACH(3) /      4 /
C     DATA IMACH(4) /      1 /
C     DATA IMACH(5) /     16 /
C     DATA IMACH(6) /      2 /
C     DATA IMACH(7) /      2 /
C     DATA IMACH(8) /     15 /
C     DATA IMACH(9) /  32767 /
C     DATA IMACH(10)/      2 /
C     DATA IMACH(11)/     23 /
C     DATA IMACH(12)/   -128 /
C     DATA IMACH(13)/    127 /
C     DATA IMACH(14)/     55 /
C     DATA IMACH(15)/   -128 /
C     DATA IMACH(16)/    127 /
C
C     MACHINE CONSTANTS FOR THE HP 9000
C
C     DATA IMACH(1)  /    5 /
C     DATA IMACH(2)  /    6 /
C     DATA IMACH(3)  /    6 /
C     DATA IMACH(3)  /    7 /
C     DATA IMACH(5)  /   32 /
C     DATA IMACH(6)  /    4 /
C     DATA IMACH(7)  /    2 /
C     DATA IMACH(8)  /   32 /
C     DATA IMACH(9)  /2147483647 /
C     DATA IMACH(10) /    2 /
C     DATA IMACH(11) /   24 /
C     DATA IMACH(12) / -126 /
C     DATA IMACH(13) /  127 /
C     DATA IMACH(14) /   53 /
C     DATA IMACH(15) /-1015 /
C     DATA IMACH(16) / 1017 /
C
C     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
C     THE XEROX SIGMA 5/7/9, THE SEL SYSTEMS 85/86, AND
C     THE PERKIN ELMER (INTERDATA) 7/32.
C
C     DATA IMACH( 1) /   5 /
C     DATA IMACH( 2) /   6 /
C     DATA IMACH( 3) /   7 /
C     DATA IMACH( 4) /   6 /
C     DATA IMACH( 5) /  32 /
C     DATA IMACH( 6) /   4 /
C     DATA IMACH( 7) /  16 /
C     DATA IMACH( 8) /  31 /
C     DATA IMACH( 9) / Z7FFFFFFF /
C     DATA IMACH(10) /  16 /
C     DATA IMACH(11) /   6 /
C     DATA IMACH(12) / -64 /
C     DATA IMACH(13) /  63 /
C     DATA IMACH(14) /  14 /
C     DATA IMACH(15) / -64 /
C     DATA IMACH(16) /  63 /
C
C     MACHINE CONSTANTS FOR THE IBM PC
C
C     DATA IMACH( 1) /     5 /
C     DATA IMACH( 2) /     6 /
C     DATA IMACH( 3) /     0 /
C     DATA IMACH( 4) /     0 /
C     DATA IMACH( 5) /    32 /
C     DATA IMACH( 6) /     4 /
C     DATA IMACH( 7) /     2 /
C     DATA IMACH( 8) /    31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /     2 /
C     DATA IMACH(11) /    24 /
C     DATA IMACH(12) /  -125 /
C     DATA IMACH(13) /   127 /
C     DATA IMACH(14) /    53 /
C     DATA IMACH(15) / -1021 /
C     DATA IMACH(16) /  1023 /
C
C     MACHINE CONSTANTS FOR THE IBM RS 6000
C
cc     DATA IMACH( 1) /          5 /
cc     DATA IMACH( 2) /          6 /
cc     DATA IMACH( 3) /          6 /
cc     DATA IMACH( 4) /          0 /
cc     DATA IMACH( 5) /         32 /
cc     DATA IMACH( 6) /          4 /
cc     DATA IMACH( 7) /          2 /
cc     DATA IMACH( 8) /         31 /
cc     DATA IMACH( 9) / 2147483647 /
cc     DATA IMACH(10) /          2 /
cc     DATA IMACH(11) /         24 /
cc     DATA IMACH(12) /       -125 /
cc     DATA IMACH(13) /        128 /
cc     DATA IMACH(14) /         53 /
cc     DATA IMACH(15) /      -1021 /
cc     DATA IMACH(16) /       1024 /
C
C     MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR)
C
C     DATA IMACH( 1) /    5 /
C     DATA IMACH( 2) /    6 /
C     DATA IMACH( 3) /    5 /
C     DATA IMACH( 4) /    6 /
C     DATA IMACH( 5) /   36 /
C     DATA IMACH( 6) /    5 /
C     DATA IMACH( 7) /    2 /
C     DATA IMACH( 8) /   35 /
C     DATA IMACH( 9) / "377777777777 /
C     DATA IMACH(10) /    2 /
C     DATA IMACH(11) /   27 /
C     DATA IMACH(12) / -128 /
C     DATA IMACH(13) /  127 /
C     DATA IMACH(14) /   54 /
C     DATA IMACH(15) / -101 /
C     DATA IMACH(16) /  127 /
C
C     MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR)
C
C     DATA IMACH( 1) /    5 /
C     DATA IMACH( 2) /    6 /
C     DATA IMACH( 3) /    5 /
C     DATA IMACH( 4) /    6 /
C     DATA IMACH( 5) /   36 /
C     DATA IMACH( 6) /    5 /
C     DATA IMACH( 7) /    2 /
C     DATA IMACH( 8) /   35 /
C     DATA IMACH( 9) / "377777777777 /
C     DATA IMACH(10) /    2 /
C     DATA IMACH(11) /   27 /
C     DATA IMACH(12) / -128 /
C     DATA IMACH(13) /  127 /
C     DATA IMACH(14) /   62 /
C     DATA IMACH(15) / -128 /
C     DATA IMACH(16) /  127 /
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
C     32-BIT INTEGER ARITHMETIC.
C
C     DATA IMACH( 1) /    5 /
C     DATA IMACH( 2) /    6 /
C     DATA IMACH( 3) /    5 /
C     DATA IMACH( 4) /    6 /
C     DATA IMACH( 5) /   32 /
C     DATA IMACH( 6) /    4 /
C     DATA IMACH( 7) /    2 /
C     DATA IMACH( 8) /   31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /    2 /
C     DATA IMACH(11) /   24 /
C     DATA IMACH(12) / -127 /
C     DATA IMACH(13) /  127 /
C     DATA IMACH(14) /   56 /
C     DATA IMACH(15) / -127 /
C     DATA IMACH(16) /  127 /
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
C     16-BIT INTEGER ARITHMETIC.
C
C     DATA IMACH( 1) /    5 /
C     DATA IMACH( 2) /    6 /
C     DATA IMACH( 3) /    5 /
C     DATA IMACH( 4) /    6 /
C     DATA IMACH( 5) /   16 /
C     DATA IMACH( 6) /    2 /
C     DATA IMACH( 7) /    2 /
C     DATA IMACH( 8) /   15 /
C     DATA IMACH( 9) / 32767 /
C     DATA IMACH(10) /    2 /
C     DATA IMACH(11) /   24 /
C     DATA IMACH(12) / -127 /
C     DATA IMACH(13) /  127 /
C     DATA IMACH(14) /   56 /
C     DATA IMACH(15) / -127 /
C     DATA IMACH(16) /  127 /
C
C     MACHINE CONSTANTS FOR THE SILICON GRAPHICS IRIS
C
      DATA IMACH( 1) /     5 /
      DATA IMACH( 2) /     6 /
      DATA IMACH( 3) /     6 /
      DATA IMACH( 4) /     0 /
      DATA IMACH( 5) /    32 /
      DATA IMACH( 6) /     4 /
      DATA IMACH( 7) /     2 /
      DATA IMACH( 8) /    31 /
      DATA IMACH( 9) / 2147483647 /
      DATA IMACH(10) /     2 /
      DATA IMACH(11) /    23 /
      DATA IMACH(12) /  -126 /
      DATA IMACH(13) /   127 /
      DATA IMACH(14) /    52 /
      DATA IMACH(15) / -1022 /
      DATA IMACH(16) /  1023 /
C
C     MACHINE CONSTANTS FOR THE SUN
C
C     DATA IMACH(1) /    5 /
C     DATA IMACH(2) /    6 /
C     DATA IMACH(3) /    6 /
C     DATA IMACH(4) /    6 /
C     DATA IMACH(5) /   32 /
C     DATA IMACH(6) /    4 /
C     DATA IMACH(7) /    2 /
C     DATA IMACH(8) /   31 /
C     DATA IMACH(9) /2147483647 /
C     DATA IMACH(10)/    2 /
C     DATA IMACH(11)/   24 /
C     DATA IMACH(12)/ -125 /
C     DATA IMACH(13)/  128 /
C     DATA IMACH(14)/   53 /
C     DATA IMACH(15)/ -1021 /
C     DATA IMACH(16)/  1024 /
C
C     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES FTN COMPILER
C
C
C     DATA IMACH( 1) /    5 /
C     DATA IMACH( 2) /    6 /
C     DATA IMACH( 3) /    1 /
C     DATA IMACH( 4) /    6 /
C     DATA IMACH( 5) /   36 /
C     DATA IMACH( 6) /    4 /
C     DATA IMACH( 7) /    2 /
C     DATA IMACH( 8) /   35 /
C     DATA IMACH( 9) / O377777777777 /
C     DATA IMACH(10) /    2 /
C     DATA IMACH(11) /   27 /
C     DATA IMACH(12) / -128 /
C     DATA IMACH(13) /  127 /
C     DATA IMACH(14) /   60 /
C     DATA IMACH(15) /-1024 /
C     DATA IMACH(16) / 1023 /
C
C     MACHINE CONSTANTS FOR THE VAX 11/780
C
C     DATA IMACH(1) /    5 /
C     DATA IMACH(2) /    6 /
C     DATA IMACH(3) /    5 /
C     DATA IMACH(4) /    6 /
C     DATA IMACH(5) /   32 /
C     DATA IMACH(6) /    4 /
C     DATA IMACH(7) /    2 /
C     DATA IMACH(8) /   31 /
C     DATA IMACH(9) /2147483647 /
C     DATA IMACH(10)/    2 /
C     DATA IMACH(11)/   24 /
C     DATA IMACH(12)/ -127 /
C     DATA IMACH(13)/  127 /
C     DATA IMACH(14)/   56 /
C     DATA IMACH(15)/ -127 /
C     DATA IMACH(16)/  127 /
C
C     MACHINE CONSTANTS FOR THE VAX 11/780, G-FLOAT OPTION
C
C     DATA IMACH(1) /    5 /
C     DATA IMACH(2) /    6 /
C     DATA IMACH(3) /    5 /
C     DATA IMACH(4) /    6 /
C     DATA IMACH(5) /   32 /
C     DATA IMACH(6) /    4 /
C     DATA IMACH(7) /    2 /
C     DATA IMACH(8) /   31 /
C     DATA IMACH(9) /2147483647 /
C     DATA IMACH(10)/    2 /
C     DATA IMACH(11)/   24 /
C     DATA IMACH(12)/ -127 /
C     DATA IMACH(13)/  127 /
C     DATA IMACH(14)/   53 /
C     DATA IMACH(15)/ -1022 /
C     DATA IMACH(16)/  1023 /
C
C     MACHINE CONSTANTS FOR THE Z80 MICROPROCESSOR
C
C     DATA IMACH( 1) /     1/
C     DATA IMACH( 2) /     1/
C     DATA IMACH( 3) /     0/
C     DATA IMACH( 4) /     1/
C     DATA IMACH( 5) /    16/
C     DATA IMACH( 6) /     2/
C     DATA IMACH( 7) /     2/
C     DATA IMACH( 8) /    15/
C     DATA IMACH( 9) / 32767/
C     DATA IMACH(10) /     2/
C     DATA IMACH(11) /    24/
C     DATA IMACH(12) /  -127/
C     DATA IMACH(13) /   127/
C     DATA IMACH(14) /    56/
C     DATA IMACH(15) /  -127/
C     DATA IMACH(16) /   127/
C
C
C***FIRST EXECUTABLE STATEMENT  I1MACH
      IF (I .LT. 1  .OR.  I .GT. 16) GO TO 10
C
      I1MACH = IMACH(I)
      RETURN
C
   10 CONTINUE
      WRITE (UNIT = OUTPUT, FMT = 9000)
 9000 FORMAT ('1ERROR    1 IN I1MACH - I OUT OF BOUNDS')
C
C     CALL FDUMP
C
C
      STOP
      END
      SUBROUTINE XERROR(MESS,NMESS,L1,L2)
C
C     THIS IS A DUMMY XERROR ROUTINE TO PRINT ERROR MESSAGES WITH NMESS
C     CHARACTERS. L1 AND L2 ARE DUMMY PARAMETERS TO MAKE THIS CALL
C     COMPATIBLE WITH THE SLATEC XERROR ROUTINE. THIS IS A FORTRAN 77
C     ROUTINE.
C
      CHARACTER*(*) MESS
      NN=NMESS/70
      NR=NMESS-70*NN
      IF(NR.NE.0) NN=NN+1
      K=1
      PRINT 900
  900 FORMAT(/)
      DO 10 I=1,NN
        KMIN=MIN0(K+69,NMESS)
        PRINT *, MESS(K:KMIN)
        K=K+70
   10 CONTINUE
      PRINT 900
      RETURN
      END
      DOUBLE PRECISION FUNCTION ZSABS(ZR, ZI)
C***BEGIN PROLOGUE  ZSABS
C***REFER TO  ZBESH,ZBESI,ZBESJ,ZBESK,ZBESY,ZAIRY,ZBIRY
C
C     ZSABS COMPUTES THE ABSOLUTE VALUE OR MAGNITUDE OF A DOUBLE
C     PRECISION COMPLEX VARIABLE CMPLX(ZR,ZI)
C
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  ZSABS
      DOUBLE PRECISION ZR, ZI, U, V, Q, S
      U = DABS(ZR)
      V = DABS(ZI)
      S = U + V
C-----------------------------------------------------------------------
C     S*1.0D0 MAKES AN UNNORMALIZED UNDERFLOW ON CDC MACHINES INTO A
C     TRUE FLOATING ZERO
C-----------------------------------------------------------------------
      S = S*1.0D+0
      IF (S.EQ.0.0D+0) GO TO 20
      IF (U.GT.V) GO TO 10
      Q = U/V
      ZSABS = V*DSQRT(1.D+0+Q*Q)
      RETURN
   10 Q = V/U
      ZSABS = U*DSQRT(1.D+0+Q*Q)
      RETURN
   20 ZSABS = 0.0D+0
      RETURN
      END
      SUBROUTINE ZACAI(ZR, ZI, FNU, KODE, MR, N, YR, YI, NZ, RL, TOL,
     * ELIM, ALIM)
C***BEGIN PROLOGUE  ZACAI
C***REFER TO  ZAIRY
C
C     ZACAI APPLIES THE ANALYTIC CONTINUATION FORMULA
C
C         K(FNU,ZN*EXP(MP))=K(FNU,ZN)*EXP(-MP*FNU) - MP*I(FNU,ZN)
C                 MP=PI*MR*CMPLX(0.0,1.0)
C
C     TO CONTINUE THE K FUNCTION FROM THE RIGHT HALF TO THE LEFT
C     HALF Z PLANE FOR USE WITH ZAIRY WHERE FNU=1/3 OR 2/3 AND N=1.
C     ZACAI IS THE SAME AS ZACON WITH THE PARTS FOR LARGER ORDERS AND
C     RECURRENCE REMOVED. A RECURSIVE CALL TO ZACON CAN RESULT IF ZACON
C     IS CALLED FROM ZAIRY.
C
C***ROUTINES CALLED  ZASYI,ZBKNU,ZMLRI,ZSERI,ZS1S2,D1MACH,ZSABS
C***END PROLOGUE  ZACAI
C     COMPLEX CSGN,CSPN,C1,C2,Y,Z,ZN,CY
      DOUBLE PRECISION ALIM, ARG, ASCLE, AZ, CSGNR, CSGNI, CSPNR,
     * CSPNI, C1R, C1I, C2R, C2I, CYR, CYI, DFNU, ELIM, FMR, FNU, PI,
     * RL, SGN, TOL, YY, YR, YI, ZR, ZI, ZNR, ZNI, D1MACH, ZSABS
      INTEGER INU, IUF, KODE, MR, N, NN, NW, NZ
      DIMENSION YR(N), YI(N), CYR(2), CYI(2)
      DATA PI / 3.14159265358979324D0 /
      NZ = 0
      ZNR = -ZR
      ZNI = -ZI
      AZ = ZSABS(ZR,ZI)
      NN = N
      DFNU = FNU + DBLE(FLOAT(N-1))
      IF (AZ.LE.2.0D0) GO TO 10
      IF (AZ*AZ*0.25D0.GT.DFNU+1.0D0) GO TO 20
   10 CONTINUE
C-----------------------------------------------------------------------
C     POWER SERIES FOR THE I FUNCTION
C-----------------------------------------------------------------------
      CALL ZSERI(ZNR, ZNI, FNU, KODE, NN, YR, YI, NW, TOL, ELIM, ALIM)
      GO TO 40
   20 CONTINUE
      IF (AZ.LT.RL) GO TO 30
C-----------------------------------------------------------------------
C     ASYMPTOTIC EXPANSION FOR LARGE Z FOR THE I FUNCTION
C-----------------------------------------------------------------------
      CALL ZASYI(ZNR, ZNI, FNU, KODE, NN, YR, YI, NW, RL, TOL, ELIM,
     * ALIM)
      IF (NW.LT.0) GO TO 80
      GO TO 40
   30 CONTINUE
C-----------------------------------------------------------------------
C     MILLER ALGORITHM NORMALIZED BY THE SERIES FOR THE I FUNCTION
C-----------------------------------------------------------------------
      CALL ZMLRI(ZNR, ZNI, FNU, KODE, NN, YR, YI, NW, TOL)
      IF(NW.LT.0) GO TO 80
   40 CONTINUE
C-----------------------------------------------------------------------
C     ANALYTIC CONTINUATION TO THE LEFT HALF PLANE FOR THE K FUNCTION
C-----------------------------------------------------------------------
      CALL ZBKNU(ZNR, ZNI, FNU, KODE, 1, CYR, CYI, NW, TOL, ELIM, ALIM)
      IF (NW.NE.0) GO TO 80
      FMR = DBLE(FLOAT(MR))
      SGN = -DSIGN(PI,FMR)
      CSGNR = 0.0D0
      CSGNI = SGN
      IF (KODE.EQ.1) GO TO 50
      YY = -ZNI
      CSGNR = -CSGNI*DSIN(YY)
      CSGNI = CSGNI*DCOS(YY)
   50 CONTINUE
C-----------------------------------------------------------------------
C     CALCULATE CSPN=EXP(FNU*PI*I) TO MINIMIZE LOSSES OF SIGNIFICANCE
C     WHEN FNU IS LARGE
C-----------------------------------------------------------------------
      INU = INT(SNGL(FNU))
      ARG = (FNU-DBLE(FLOAT(INU)))*SGN
      CSPNR = DCOS(ARG)
      CSPNI = DSIN(ARG)
      IF (MOD(INU,2).EQ.0) GO TO 60
      CSPNR = -CSPNR
      CSPNI = -CSPNI
   60 CONTINUE
      C1R = CYR(1)
      C1I = CYI(1)
      C2R = YR(1)
      C2I = YI(1)
      IF (KODE.EQ.1) GO TO 70
      IUF = 0
      ASCLE = 1.0D+3*D1MACH(1)/TOL
      CALL ZS1S2(ZNR, ZNI, C1R, C1I, C2R, C2I, NW, ASCLE, ALIM, IUF)
      NZ = NZ + NW
   70 CONTINUE
      YR(1) = CSPNR*C1R - CSPNI*C1I + CSGNR*C2R - CSGNI*C2I
      YI(1) = CSPNR*C1I + CSPNI*C1R + CSGNR*C2I + CSGNI*C2R
      RETURN
   80 CONTINUE
      NZ = -1
      IF(NW.EQ.(-2)) NZ=-2
      RETURN
      END


      SUBROUTINE ZDIV(AR, AI, BR, BI, CR, CI)
C***BEGIN PROLOGUE  ZDIV
C***REFER TO  ZBESH,ZBESI,ZBESJ,ZBESK,ZBESY,ZAIRY,ZBIRY
C
C     DOUBLE PRECISION COMPLEX DIVIDE C=A/B.
C
C***ROUTINES CALLED  ZSABS
C***END PROLOGUE  ZDIV
      DOUBLE PRECISION AR, AI, BR, BI, CR, CI, BM, CA, CB, CC, CD
      DOUBLE PRECISION ZSABS
      BM = 1.0D0/ZSABS(BR,BI)
      CC = BR*BM
      CD = BI*BM
      CA = (AR*CC+AI*CD)*BM
      CB = (AI*CC-AR*CD)*BM
      CR = CA
      CI = CB
      RETURN
      END
      SUBROUTINE ZSEXP(AR, AI, BR, BI)
C***BEGIN PROLOGUE  ZSEXP
C***REFER TO  ZBESH,ZBESI,ZBESJ,ZBESK,ZBESY,ZAIRY,ZBIRY
C
C     DOUBLE PRECISION COMPLEX EXPONENTIAL FUNCTION B=EXP(A)
C
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  ZSEXP
      DOUBLE PRECISION AR, AI, BR, BI, ZM, CA, CB
      ZM = DEXP(AR)
      CA = ZM*DCOS(AI)
      CB = ZM*DSIN(AI)
      BR = CA
      BI = CB
      RETURN
      END
      SUBROUTINE ZKSCL(ZRR,ZRI,FNU,N,YR,YI,NZ,RZR,RZI,ASCLE,TOL,ELIM)
C***BEGIN PROLOGUE  ZKSCL
C***REFER TO  ZBESK
C
C     SET K FUNCTIONS TO ZERO ON UNDERFLOW, CONTINUE RECURRENCE
C     ON SCALED FUNCTIONS UNTIL TWO MEMBERS COME ON SCALE, THEN
C     RETURN WITH MIN(NZ+2,N) VALUES SCALED BY 1/TOL.
C
C***ROUTINES CALLED  ZUCHK,ZSABS,ZSLOG
C***END PROLOGUE  ZKSCL
C     COMPLEX CK,CS,CY,CZERO,RZ,S1,S2,Y,ZR,ZD,CELM
      DOUBLE PRECISION ACS, AS, ASCLE, CKI, CKR, CSI, CSR, CYI,
     * CYR, ELIM, FN, FNU, RZI, RZR, STR, S1I, S1R, S2I,
     * S2R, TOL, YI, YR, ZEROI, ZEROR, ZRI, ZRR, ZSABS,
     * ZDR, ZDI, CELMR, ELM, HELIM, ALAS
      INTEGER I, IC, IDUM, KK, N, NN, NW, NZ
      DIMENSION YR(N), YI(N), CYR(2), CYI(2)
      DATA ZEROR,ZEROI / 0.0D0 , 0.0D0 /
C
      NZ = 0
      IC = 0
      NN = MIN0(2,N)
      DO 10 I=1,NN
        S1R = YR(I)
        S1I = YI(I)
        CYR(I) = S1R
        CYI(I) = S1I
        AS = ZSABS(S1R,S1I)
        ACS = -ZRR + DLOG(AS)
        NZ = NZ + 1
        YR(I) = ZEROR
        YI(I) = ZEROI
        IF (ACS.LT.(-ELIM)) GO TO 10
        CALL ZSLOG(S1R, S1I, CSR, CSI, IDUM)
        CSR = CSR - ZRR
        CSI = CSI - ZRI
        STR = DEXP(CSR)/TOL
        CSR = STR*DCOS(CSI)
        CSI = STR*DSIN(CSI)
        CALL ZUCHK(CSR, CSI, NW, ASCLE, TOL)
        IF (NW.NE.0) GO TO 10
        YR(I) = CSR
        YI(I) = CSI
        IC = I
        NZ = NZ - 1
   10 CONTINUE
      IF (N.EQ.1) RETURN
      IF (IC.GT.1) GO TO 20
      YR(1) = ZEROR
      YI(1) = ZEROI
      NZ = 2
   20 CONTINUE
      IF (N.EQ.2) RETURN
      IF (NZ.EQ.0) RETURN
      FN = FNU + 1.0D0
      CKR = FN*RZR
      CKI = FN*RZI
      S1R = CYR(1)
      S1I = CYI(1)
      S2R = CYR(2)
      S2I = CYI(2)
      HELIM = 0.5D0*ELIM
      ELM = DEXP(-ELIM)
      CELMR = ELM
      ZDR = ZRR
      ZDI = ZRI
C
C     FIND TWO CONSECUTIVE Y VALUES ON SCALE. SCALE RECURRENCE IF
C     S2 GETS LARGER THAN EXP(ELIM/2)
C
      DO 30 I=3,N
        KK = I
        CSR = S2R
        CSI = S2I
        S2R = CKR*CSR - CKI*CSI + S1R
        S2I = CKI*CSR + CKR*CSI + S1I
        S1R = CSR
        S1I = CSI
        CKR = CKR + RZR
        CKI = CKI + RZI
        AS = ZSABS(S2R,S2I)
        ALAS = DLOG(AS)
        ACS = -ZDR + ALAS
        NZ = NZ + 1
        YR(I) = ZEROR
        YI(I) = ZEROI
        IF (ACS.LT.(-ELIM)) GO TO 25
        CALL ZSLOG(S2R, S2I, CSR, CSI, IDUM)
        CSR = CSR - ZDR
        CSI = CSI - ZDI
        STR = DEXP(CSR)/TOL
        CSR = STR*DCOS(CSI)
        CSI = STR*DSIN(CSI)
        CALL ZUCHK(CSR, CSI, NW, ASCLE, TOL)
        IF (NW.NE.0) GO TO 25
        YR(I) = CSR
        YI(I) = CSI
        NZ = NZ - 1
        IF (IC.EQ.KK-1) GO TO 40
        IC = KK
        GO TO 30
   25   CONTINUE
        IF(ALAS.LT.HELIM) GO TO 30
        ZDR = ZDR - ELIM
        S1R = S1R*CELMR
        S1I = S1I*CELMR
        S2R = S2R*CELMR
        S2I = S2I*CELMR
   30 CONTINUE
      NZ = N
      IF(IC.EQ.N) NZ=N-1
      GO TO 45
   40 CONTINUE
      NZ = KK - 2
   45 CONTINUE
      DO 50 I=1,NZ
        YR(I) = ZEROR
        YI(I) = ZEROI
   50 CONTINUE
      RETURN
      END
      SUBROUTINE ZSLOG(AR, AI, BR, BI, IERR)
C***BEGIN PROLOGUE  ZSLOG
C***REFER TO  ZBESH,ZBESI,ZBESJ,ZBESK,ZBESY,ZAIRY,ZBIRY
C
C     DOUBLE PRECISION COMPLEX LOGARITHM B=CLOG(A)
C     IERR=0,NORMAL RETURN      IERR=1, Z=CMPLX(0.0,0.0)
C***ROUTINES CALLED  ZSABS
C***END PROLOGUE  ZSLOG
      DOUBLE PRECISION AR, AI, BR, BI, ZM, DTHETA, DPI, DHPI
      DOUBLE PRECISION ZSABS
      DATA DPI , DHPI  / 3.141592653589793238462643383D+0,
     1                   1.570796326794896619231321696D+0/
C
      IERR=0
      IF (AR.EQ.0.0D+0) GO TO 10
      IF (AI.EQ.0.0D+0) GO TO 20
      DTHETA = DATAN(AI/AR)
      IF (DTHETA.LE.0.0D+0) GO TO 40
      IF (AR.LT.0.0D+0) DTHETA = DTHETA - DPI
      GO TO 50
   10 IF (AI.EQ.0.0D+0) GO TO 60
      BI = DHPI
      BR = DLOG(DABS(AI))
      IF (AI.LT.0.0D+0) BI = -BI
      RETURN
   20 IF (AR.GT.0.0D+0) GO TO 30
      BR = DLOG(DABS(AR))
      BI = DPI
      RETURN
   30 BR = DLOG(AR)
      BI = 0.0D+0
      RETURN
   40 IF (AR.LT.0.0D+0) DTHETA = DTHETA + DPI
   50 ZM = ZSABS(AR,AI)
      BR = DLOG(ZM)
      BI = DTHETA
      RETURN
   60 CONTINUE
      IERR=1
      RETURN
      END
      SUBROUTINE ZMLRI(ZR, ZI, FNU, KODE, N, YR, YI, NZ, TOL)
C***BEGIN PROLOGUE  ZMLRI
C***REFER TO  ZBESI,ZBESK
C
C     ZMLRI COMPUTES THE I BESSEL FUNCTION FOR RE(Z).GE.0.0 BY THE
C     MILLER ALGORITHM NORMALIZED BY A NEUMANN SERIES.
C
C***ROUTINES CALLED  DGAMLN,D1MACH,ZSABS,ZSEXP,ZSLOG,ZMLT
C***END PROLOGUE  ZMLRI
C     COMPLEX CK,CNORM,CONE,CTWO,CZERO,PT,P1,P2,RZ,SUM,Y,Z
      DOUBLE PRECISION ACK, AK, AP, AT, AZ, BK, CKI, CKR, CNORMI,
     * CNORMR, CONEI, CONER, FKAP, FKK, FLAM, FNF, FNU, PTI, PTR, P1I,
     * P1R, P2I, P2R, RAZ, RHO, RHO2, RZI, RZR, SCLE, STI, STR, SUMI,
     * SUMR, TFNF, TOL, TST, YI, YR, ZEROI, ZEROR, ZI, ZR, DGAMLN,
     * D1MACH, ZSABS
      INTEGER I, IAZ, IDUM, IFNU, INU, ITIME, K, KK, KM, KODE, M, N, NZ
      DIMENSION YR(N), YI(N)
      DATA ZEROR,ZEROI,CONER,CONEI / 0.0D0, 0.0D0, 1.0D0, 0.0D0 /
      SCLE = D1MACH(1)/TOL
      NZ=0
      AZ = ZSABS(ZR,ZI)
      IAZ = INT(SNGL(AZ))
      IFNU = INT(SNGL(FNU))
      INU = IFNU + N - 1
      AT = DBLE(FLOAT(IAZ)) + 1.0D0
      RAZ = 1.0D0/AZ
      STR = ZR*RAZ
      STI = -ZI*RAZ
      CKR = STR*AT*RAZ
      CKI = STI*AT*RAZ
      RZR = (STR+STR)*RAZ
      RZI = (STI+STI)*RAZ
      P1R = ZEROR
      P1I = ZEROI
      P2R = CONER
      P2I = CONEI
      ACK = (AT+1.0D0)*RAZ
      RHO = ACK + DSQRT(ACK*ACK-1.0D0)
      RHO2 = RHO*RHO
      TST = (RHO2+RHO2)/((RHO2-1.0D0)*(RHO-1.0D0))
      TST = TST/TOL
C-----------------------------------------------------------------------
C     COMPUTE RELATIVE TRUNCATION ERROR INDEX FOR SERIES
C-----------------------------------------------------------------------
      AK = AT
      DO 10 I=1,80
        PTR = P2R
        PTI = P2I
        P2R = P1R - (CKR*PTR-CKI*PTI)
        P2I = P1I - (CKI*PTR+CKR*PTI)
        P1R = PTR
        P1I = PTI
        CKR = CKR + RZR
        CKI = CKI + RZI
        AP = ZSABS(P2R,P2I)
        IF (AP.GT.TST*AK*AK) GO TO 20
        AK = AK + 1.0D0
   10 CONTINUE
      GO TO 110
   20 CONTINUE
      I = I + 1
      K = 0
      IF (INU.LT.IAZ) GO TO 40
C-----------------------------------------------------------------------
C     COMPUTE RELATIVE TRUNCATION ERROR FOR RATIOS
C-----------------------------------------------------------------------
      P1R = ZEROR
      P1I = ZEROI
      P2R = CONER
      P2I = CONEI
      AT = DBLE(FLOAT(INU)) + 1.0D0
      STR = ZR*RAZ
      STI = -ZI*RAZ
      CKR = STR*AT*RAZ
      CKI = STI*AT*RAZ
      ACK = AT*RAZ
      TST = DSQRT(ACK/TOL)
      ITIME = 1
      DO 30 K=1,80
        PTR = P2R
        PTI = P2I
        P2R = P1R - (CKR*PTR-CKI*PTI)
        P2I = P1I - (CKR*PTI+CKI*PTR)
        P1R = PTR
        P1I = PTI
        CKR = CKR + RZR
        CKI = CKI + RZI
        AP = ZSABS(P2R,P2I)
        IF (AP.LT.TST) GO TO 30
        IF (ITIME.EQ.2) GO TO 40
        ACK = ZSABS(CKR,CKI)
        FLAM = ACK + DSQRT(ACK*ACK-1.0D0)
        FKAP = AP/ZSABS(P1R,P1I)
        RHO = DMIN1(FLAM,FKAP)
        TST = TST*DSQRT(RHO/(RHO*RHO-1.0D0))
        ITIME = 2
   30 CONTINUE
      GO TO 110
   40 CONTINUE
C-----------------------------------------------------------------------
C     BACKWARD RECURRENCE AND SUM NORMALIZING RELATION
C-----------------------------------------------------------------------
      K = K + 1
      KK = MAX0(I+IAZ,K+INU)
      FKK = DBLE(FLOAT(KK))
      P1R = ZEROR
      P1I = ZEROI
C-----------------------------------------------------------------------
C     SCALE P2 AND SUM BY SCLE
C-----------------------------------------------------------------------
      P2R = SCLE
      P2I = ZEROI
      FNF = FNU - DBLE(FLOAT(IFNU))
      TFNF = FNF + FNF
      BK = DGAMLN(FKK+TFNF+1.0D0,IDUM) - DGAMLN(FKK+1.0D0,IDUM) -
     * DGAMLN(TFNF+1.0D0,IDUM)
      BK = DEXP(BK)
      SUMR = ZEROR
      SUMI = ZEROI
      KM = KK - INU
      DO 50 I=1,KM
        PTR = P2R
        PTI = P2I
        P2R = P1R + (FKK+FNF)*(RZR*PTR-RZI*PTI)
        P2I = P1I + (FKK+FNF)*(RZI*PTR+RZR*PTI)
        P1R = PTR
        P1I = PTI
        AK = 1.0D0 - TFNF/(FKK+TFNF)
        ACK = BK*AK
        SUMR = SUMR + (ACK+BK)*P1R
        SUMI = SUMI + (ACK+BK)*P1I
        BK = ACK
        FKK = FKK - 1.0D0
   50 CONTINUE
      YR(N) = P2R
      YI(N) = P2I
      IF (N.EQ.1) GO TO 70
      DO 60 I=2,N
        PTR = P2R
        PTI = P2I
        P2R = P1R + (FKK+FNF)*(RZR*PTR-RZI*PTI)
        P2I = P1I + (FKK+FNF)*(RZI*PTR+RZR*PTI)
        P1R = PTR
        P1I = PTI
        AK = 1.0D0 - TFNF/(FKK+TFNF)
        ACK = BK*AK
        SUMR = SUMR + (ACK+BK)*P1R
        SUMI = SUMI + (ACK+BK)*P1I
        BK = ACK
        FKK = FKK - 1.0D0
        M = N - I + 1
        YR(M) = P2R
        YI(M) = P2I
   60 CONTINUE
   70 CONTINUE
      IF (IFNU.LE.0) GO TO 90
      DO 80 I=1,IFNU
        PTR = P2R
        PTI = P2I
        P2R = P1R + (FKK+FNF)*(RZR*PTR-RZI*PTI)
        P2I = P1I + (FKK+FNF)*(RZR*PTI+RZI*PTR)
        P1R = PTR
        P1I = PTI
        AK = 1.0D0 - TFNF/(FKK+TFNF)
        ACK = BK*AK
        SUMR = SUMR + (ACK+BK)*P1R
        SUMI = SUMI + (ACK+BK)*P1I
        BK = ACK
        FKK = FKK - 1.0D0
   80 CONTINUE
   90 CONTINUE
      PTR = ZR
      PTI = ZI
      IF (KODE.EQ.2) PTR = ZEROR
      CALL ZSLOG(RZR, RZI, STR, STI, IDUM)
      P1R = -FNF*STR + PTR
      P1I = -FNF*STI + PTI
      AP = DGAMLN(1.0D0+FNF,IDUM)
      PTR = P1R - AP
      PTI = P1I
C-----------------------------------------------------------------------
C     THE DIVISION CEXP(PT)/(SUM+P2) IS ALTERED TO AVOID OVERFLOW
C     IN THE DENOMINATOR BY SQUARING LARGE QUANTITIES
C-----------------------------------------------------------------------
      P2R = P2R + SUMR
      P2I = P2I + SUMI
      AP = ZSABS(P2R,P2I)
      P1R = 1.0D0/AP
      CALL ZSEXP(PTR, PTI, STR, STI)
      CKR = STR*P1R
      CKI = STI*P1R
      PTR = P2R*P1R
      PTI = -P2I*P1R
      CALL ZMLT(CKR, CKI, PTR, PTI, CNORMR, CNORMI)
      DO 100 I=1,N
        STR = YR(I)*CNORMR - YI(I)*CNORMI
        YI(I) = YR(I)*CNORMI + YI(I)*CNORMR
        YR(I) = STR
  100 CONTINUE
      RETURN
  110 CONTINUE
      NZ=-2
      RETURN
      END
      SUBROUTINE ZMLT(AR, AI, BR, BI, CR, CI)
C***BEGIN PROLOGUE  ZMLT
C***REFER TO  ZBESH,ZBESI,ZBESJ,ZBESK,ZBESY,ZAIRY,ZBIRY
C
C     DOUBLE PRECISION COMPLEX MULTIPLY, C=A*B.
C
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  ZMLT
      DOUBLE PRECISION AR, AI, BR, BI, CR, CI, CA, CB
      CA = AR*BR - AI*BI
      CB = AR*BI + AI*BR
      CR = CA
      CI = CB
      RETURN
      END
      SUBROUTINE ZRATI(ZR, ZI, FNU, N, CYR, CYI, TOL)
C***BEGIN PROLOGUE  ZRATI
C***REFER TO  ZBESI,ZBESK,ZBESH
C
C     ZRATI COMPUTES RATIOS OF I BESSEL FUNCTIONS BY BACKWARD
C     RECURRENCE.  THE STARTING INDEX IS DETERMINED BY FORWARD
C     RECURRENCE AS DESCRIBED IN J. RES. OF NAT. BUR. OF STANDARDS-B,
C     MATHEMATICAL SCIENCES, VOL 77B, P111-114, SEPTEMBER, 1973,
C     BESSEL FUNCTIONS I AND J OF COMPLEX ARGUMENT AND INTEGER ORDER,
C     BY D. J. SOOKNE.
C
C***ROUTINES CALLED  ZSABS,ZDIV
C***END PROLOGUE  ZRATI
C     COMPLEX Z,CY(1),CONE,CZERO,P1,P2,T1,RZ,PT,CDFNU
      DOUBLE PRECISION AK, AMAGZ, AP1, AP2, ARG, AZ, CDFNUI, CDFNUR,
     * CONEI, CONER, CYI, CYR, CZEROI, CZEROR, DFNU, FDNU, FLAM, FNU,
     * FNUP, PTI, PTR, P1I, P1R, P2I, P2R, RAK, RAP1, RHO, RT2, RZI,
     * RZR, TEST, TEST1, TOL, TTI, TTR, T1I, T1R, ZI, ZR, ZSABS
      INTEGER I, ID, IDNU, INU, ITIME, K, KK, MAGZ, N
      DIMENSION CYR(N), CYI(N)
      DATA CZEROR,CZEROI,CONER,CONEI,RT2/
     1 0.0D0, 0.0D0, 1.0D0, 0.0D0, 1.41421356237309505D0 /
      AZ = ZSABS(ZR,ZI)
      INU = INT(SNGL(FNU))
      IDNU = INU + N - 1
      MAGZ = INT(SNGL(AZ))
      AMAGZ = DBLE(FLOAT(MAGZ+1))
      FDNU = DBLE(FLOAT(IDNU))
      FNUP = DMAX1(AMAGZ,FDNU)
      ID = IDNU - MAGZ - 1
      ITIME = 1
      K = 1
      PTR = 1.0D0/AZ
      RZR = PTR*(ZR+ZR)*PTR
      RZI = -PTR*(ZI+ZI)*PTR
      T1R = RZR*FNUP
      T1I = RZI*FNUP
      P2R = -T1R
      P2I = -T1I
      P1R = CONER
      P1I = CONEI
      T1R = T1R + RZR
      T1I = T1I + RZI
      IF (ID.GT.0) ID = 0
      AP2 = ZSABS(P2R,P2I)
      AP1 = ZSABS(P1R,P1I)
C-----------------------------------------------------------------------
C     THE OVERFLOW TEST ON K(FNU+I-1,Z) BEFORE THE CALL TO CBKNU
C     GUARANTEES THAT P2 IS ON SCALE. SCALE TEST1 AND ALL SUBSEQUENT
C     P2 VALUES BY AP1 TO ENSURE THAT AN OVERFLOW DOES NOT OCCUR
C     PREMATURELY.
C-----------------------------------------------------------------------
      ARG = (AP2+AP2)/(AP1*TOL)
      TEST1 = DSQRT(ARG)
      TEST = TEST1
      RAP1 = 1.0D0/AP1
      P1R = P1R*RAP1
      P1I = P1I*RAP1
      P2R = P2R*RAP1
      P2I = P2I*RAP1
      AP2 = AP2*RAP1
   10 CONTINUE
      K = K + 1
      AP1 = AP2
      PTR = P2R
      PTI = P2I
      P2R = P1R - (T1R*PTR-T1I*PTI)
      P2I = P1I - (T1R*PTI+T1I*PTR)
      P1R = PTR
      P1I = PTI
      T1R = T1R + RZR
      T1I = T1I + RZI
      AP2 = ZSABS(P2R,P2I)
      IF (AP1.LE.TEST) GO TO 10
      IF (ITIME.EQ.2) GO TO 20
      AK = ZSABS(T1R,T1I)*0.5D0
      FLAM = AK + DSQRT(AK*AK-1.0D0)
      RHO = DMIN1(AP2/AP1,FLAM)
      TEST = TEST1*DSQRT(RHO/(RHO*RHO-1.0D0))
      ITIME = 2
      GO TO 10
   20 CONTINUE
      KK = K + 1 - ID
      AK = DBLE(FLOAT(KK))
      T1R = AK
      T1I = CZEROI
      DFNU = FNU + DBLE(FLOAT(N-1))
      P1R = 1.0D0/AP2
      P1I = CZEROI
      P2R = CZEROR
      P2I = CZEROI
      DO 30 I=1,KK
        PTR = P1R
        PTI = P1I
        RAP1 = DFNU + T1R
        TTR = RZR*RAP1
        TTI = RZI*RAP1
        P1R = (PTR*TTR-PTI*TTI) + P2R
        P1I = (PTR*TTI+PTI*TTR) + P2I
        P2R = PTR
        P2I = PTI
        T1R = T1R - CONER
   30 CONTINUE
      IF (P1R.NE.CZEROR .OR. P1I.NE.CZEROI) GO TO 40
      P1R = TOL
      P1I = TOL
   40 CONTINUE
      CALL ZDIV(P2R, P2I, P1R, P1I, CYR(N), CYI(N))
      IF (N.EQ.1) RETURN
      K = N - 1
      AK = DBLE(FLOAT(K))
      T1R = AK
      T1I = CZEROI
      CDFNUR = FNU*RZR
      CDFNUI = FNU*RZI
      DO 60 I=2,N
        PTR = CDFNUR + (T1R*RZR-T1I*RZI) + CYR(K+1)
        PTI = CDFNUI + (T1R*RZI+T1I*RZR) + CYI(K+1)
        AK = ZSABS(PTR,PTI)
        IF (AK.NE.CZEROR) GO TO 50
        PTR = TOL
        PTI = TOL
        AK = TOL*RT2
   50   CONTINUE
        RAK = CONER/AK
        CYR(K) = RAK*PTR*RAK
        CYI(K) = -RAK*PTI*RAK
        T1R = T1R - CONER
        K = K - 1
   60 CONTINUE
      RETURN
      END
      SUBROUTINE ZS1S2(ZRR, ZRI, S1R, S1I, S2R, S2I, NZ, ASCLE, ALIM,
     * IUF)
C***BEGIN PROLOGUE  ZS1S2
C***REFER TO  ZBESK,ZAIRY
C
C     ZS1S2 TESTS FOR A POSSIBLE UNDERFLOW RESULTING FROM THE
C     ADDITION OF THE I AND K FUNCTIONS IN THE ANALYTIC CON-
C     TINUATION FORMULA WHERE S1=K FUNCTION AND S2=I FUNCTION.
C     ON KODE=1 THE I AND K FUNCTIONS ARE DIFFERENT ORDERS OF
C     MAGNITUDE, BUT FOR KODE=2 THEY CAN BE OF THE SAME ORDER
C     OF MAGNITUDE AND THE MAXIMUM MUST BE AT LEAST ONE
C     PRECISION ABOVE THE UNDERFLOW LIMIT.
C
C***ROUTINES CALLED  ZSABS,ZSEXP,ZSLOG
C***END PROLOGUE  ZS1S2
C     COMPLEX CZERO,C1,S1,S1D,S2,ZR
      DOUBLE PRECISION AA, ALIM, ALN, ASCLE, AS1, AS2, C1I, C1R, S1DI,
     * S1DR, S1I, S1R, S2I, S2R, ZEROI, ZEROR, ZRI, ZRR, ZSABS
      INTEGER IUF, IDUM, NZ
      DATA ZEROR,ZEROI  / 0.0D0 , 0.0D0 /
      NZ = 0
      AS1 = ZSABS(S1R,S1I)
      AS2 = ZSABS(S2R,S2I)
      IF (S1R.EQ.0.0D0 .AND. S1I.EQ.0.0D0) GO TO 10
      IF (AS1.EQ.0.0D0) GO TO 10
      ALN = -ZRR - ZRR + DLOG(AS1)
      S1DR = S1R
      S1DI = S1I
      S1R = ZEROR
      S1I = ZEROI
      AS1 = ZEROR
      IF (ALN.LT.(-ALIM)) GO TO 10
      CALL ZSLOG(S1DR, S1DI, C1R, C1I, IDUM)
      C1R = C1R - ZRR - ZRR
      C1I = C1I - ZRI - ZRI
      CALL ZSEXP(C1R, C1I, S1R, S1I)
      AS1 = ZSABS(S1R,S1I)
      IUF = IUF + 1
   10 CONTINUE
      AA = DMAX1(AS1,AS2)
      IF (AA.GT.ASCLE) RETURN
      S1R = ZEROR
      S1I = ZEROI
      S2R = ZEROR
      S2I = ZEROI
      NZ = 1
      IUF = 0
      RETURN
      END
      SUBROUTINE ZSERI(ZR, ZI, FNU, KODE, N, YR, YI, NZ, TOL, ELIM,
     * ALIM)
C***BEGIN PROLOGUE  ZSERI
C***REFER TO  ZBESI,ZBESK
C
C     ZSERI COMPUTES THE I BESSEL FUNCTION FOR REAL(Z).GE.0.0 BY
C     MEANS OF THE POWER SERIES FOR LARGE CABS(Z) IN THE
C     REGION CABS(Z).LE.2*SQRT(FNU+1). NZ=0 IS A NORMAL RETURN.
C     NZ.GT.0 MEANS THAT THE LAST NZ COMPONENTS WERE SET TO ZERO
C     DUE TO UNDERFLOW. NZ.LT.0 MEANS UNDERFLOW OCCURRED, BUT THE
C     CONDITION CABS(Z).LE.2*SQRT(FNU+1) WAS VIOLATED AND THE
C     COMPUTATION MUST BE COMPLETED IN ANOTHER ROUTINE WITH N=N-ABS(NZ).
C
C***ROUTINES CALLED  DGAMLN,D1MACH,ZUCHK,ZSABS,ZDIV,ZSLOG,ZMLT
C***END PROLOGUE  ZSERI
C     COMPLEX AK1,CK,COEF,CONE,CRSC,CSCL,CZ,CZERO,HZ,RZ,S1,S2,Y,Z
      DOUBLE PRECISION AA, ACZ, AK, AK1I, AK1R, ALIM, ARM, ASCLE, ATOL,
     * AZ, CKI, CKR, COEFI, COEFR, CONEI, CONER, CRSCR, CZI, CZR, DFNU,
     * ELIM, FNU, FNUP, HZI, HZR, RAZ, RS, RTR1, RZI, RZR, S, SS, STI,
     * STR, S1I, S1R, S2I, S2R, TOL, YI, YR, WI, WR, ZEROI, ZEROR, ZI,
     * ZR, DGAMLN, D1MACH, ZSABS
      INTEGER I, IB, IDUM, IFLAG, IL, K, KODE, L, M, N, NN, NZ, NW
      DIMENSION YR(N), YI(N), WR(2), WI(2)
      DATA ZEROR,ZEROI,CONER,CONEI / 0.0D0, 0.0D0, 1.0D0, 0.0D0 /
C
      NZ = 0
      AZ = ZSABS(ZR,ZI)
      IF (AZ.EQ.0.0D0) GO TO 160
      ARM = 1.0D+3*D1MACH(1)
      RTR1 = DSQRT(ARM)
      CRSCR = 1.0D0
      IFLAG = 0
      IF (AZ.LT.ARM) GO TO 150
      HZR = 0.5D0*ZR
      HZI = 0.5D0*ZI
      CZR = ZEROR
      CZI = ZEROI
      IF (AZ.LE.RTR1) GO TO 10
      CALL ZMLT(HZR, HZI, HZR, HZI, CZR, CZI)
   10 CONTINUE
      ACZ = ZSABS(CZR,CZI)
      NN = N
      CALL ZSLOG(HZR, HZI, CKR, CKI, IDUM)
   20 CONTINUE
      DFNU = FNU + DBLE(FLOAT(NN-1))
      FNUP = DFNU + 1.0D0
C-----------------------------------------------------------------------
C     UNDERFLOW TEST
C-----------------------------------------------------------------------
      AK1R = CKR*DFNU
      AK1I = CKI*DFNU
      AK = DGAMLN(FNUP,IDUM)
      AK1R = AK1R - AK
      IF (KODE.EQ.2) AK1R = AK1R - ZR
      IF (AK1R.GT.(-ELIM)) GO TO 40
   30 CONTINUE
      NZ = NZ + 1
      YR(NN) = ZEROR
      YI(NN) = ZEROI
      IF (ACZ.GT.DFNU) GO TO 190
      NN = NN - 1
      IF (NN.EQ.0) RETURN
      GO TO 20
   40 CONTINUE
      IF (AK1R.GT.(-ALIM)) GO TO 50
      IFLAG = 1
      SS = 1.0D0/TOL
      CRSCR = TOL
      ASCLE = ARM*SS
   50 CONTINUE
      AA = DEXP(AK1R)
      IF (IFLAG.EQ.1) AA = AA*SS
      COEFR = AA*DCOS(AK1I)
      COEFI = AA*DSIN(AK1I)
      ATOL = TOL*ACZ/FNUP
      IL = MIN0(2,NN)
      DO 90 I=1,IL
        DFNU = FNU + DBLE(FLOAT(NN-I))
        FNUP = DFNU + 1.0D0
        S1R = CONER
        S1I = CONEI
        IF (ACZ.LT.TOL*FNUP) GO TO 70
        AK1R = CONER
        AK1I = CONEI
        AK = FNUP + 2.0D0
        S = FNUP
        AA = 2.0D0
   60   CONTINUE
        RS = 1.0D0/S
        STR = AK1R*CZR - AK1I*CZI
        STI = AK1R*CZI + AK1I*CZR
        AK1R = STR*RS
        AK1I = STI*RS
        S1R = S1R + AK1R
        S1I = S1I + AK1I
        S = S + AK
        AK = AK + 2.0D0
        AA = AA*ACZ*RS
        IF (AA.GT.ATOL) GO TO 60
   70   CONTINUE
        S2R = S1R*COEFR - S1I*COEFI
        S2I = S1R*COEFI + S1I*COEFR
        WR(I) = S2R
        WI(I) = S2I
        IF (IFLAG.EQ.0) GO TO 80
        CALL ZUCHK(S2R, S2I, NW, ASCLE, TOL)
        IF (NW.NE.0) GO TO 30
   80   CONTINUE
        M = NN - I + 1
        YR(M) = S2R*CRSCR
        YI(M) = S2I*CRSCR
        IF (I.EQ.IL) GO TO 90
        CALL ZDIV(COEFR, COEFI, HZR, HZI, STR, STI)
        COEFR = STR*DFNU
        COEFI = STI*DFNU
   90 CONTINUE
      IF (NN.LE.2) RETURN
      K = NN - 2
      AK = DBLE(FLOAT(K))
      RAZ = 1.0D0/AZ
      STR = ZR*RAZ
      STI = -ZI*RAZ
      RZR = (STR+STR)*RAZ
      RZI = (STI+STI)*RAZ
      IF (IFLAG.EQ.1) GO TO 120
      IB = 3
  100 CONTINUE
      DO 110 I=IB,NN
        YR(K) = (AK+FNU)*(RZR*YR(K+1)-RZI*YI(K+1)) + YR(K+2)
        YI(K) = (AK+FNU)*(RZR*YI(K+1)+RZI*YR(K+1)) + YI(K+2)
        AK = AK - 1.0D0
        K = K - 1
  110 CONTINUE
      RETURN
C-----------------------------------------------------------------------
C     RECUR BACKWARD WITH SCALED VALUES
C-----------------------------------------------------------------------
  120 CONTINUE
C-----------------------------------------------------------------------
C     EXP(-ALIM)=EXP(-ELIM)/TOL=APPROX. ONE PRECISION ABOVE THE
C     UNDERFLOW LIMIT = ASCLE = D1MACH(1)*SS*1.0D+3
C-----------------------------------------------------------------------
      S1R = WR(1)
      S1I = WI(1)
      S2R = WR(2)
      S2I = WI(2)
      DO 130 L=3,NN
        CKR = S2R
        CKI = S2I
        S2R = S1R + (AK+FNU)*(RZR*CKR-RZI*CKI)
        S2I = S1I + (AK+FNU)*(RZR*CKI+RZI*CKR)
        S1R = CKR
        S1I = CKI
        CKR = S2R*CRSCR
        CKI = S2I*CRSCR
        YR(K) = CKR
        YI(K) = CKI
        AK = AK - 1.0D0
        K = K - 1
        IF (ZSABS(CKR,CKI).GT.ASCLE) GO TO 140
  130 CONTINUE
      RETURN
  140 CONTINUE
      IB = L + 1
      IF (IB.GT.NN) RETURN
      GO TO 100
  150 CONTINUE
      NZ = N
      IF (FNU.EQ.0.0D0) NZ = NZ - 1
  160 CONTINUE
      YR(1) = ZEROR
      YI(1) = ZEROI
      IF (FNU.NE.0.0D0) GO TO 170
      YR(1) = CONER
      YI(1) = CONEI
  170 CONTINUE
      IF (N.EQ.1) RETURN
      DO 180 I=2,N
        YR(I) = ZEROR
        YI(I) = ZEROI
  180 CONTINUE
      RETURN
C-----------------------------------------------------------------------
C     RETURN WITH NZ.LT.0 IF CABS(Z*Z/4).GT.FNU+N-NZ-1 COMPLETE
C     THE CALCULATION IN CBINU WITH N=N-IABS(NZ)
C-----------------------------------------------------------------------
  190 CONTINUE
      NZ = -NZ
      RETURN
      END
      SUBROUTINE ZSHCH(ZR, ZI, CSHR, CSHI, CCHR, CCHI)
C***BEGIN PROLOGUE  ZSHCH
C***REFER TO  ZBESK,ZBESH
C
C     ZSHCH COMPUTES THE COMPLEX HYPERBOLIC FUNCTIONS CSH=SINH(X+I*Y)
C     AND CCH=COSH(X+I*Y), WHERE I**2=-1.
C
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  ZSHCH
C
      DOUBLE PRECISION CCHI, CCHR, CH, CN, CSHI, CSHR, SH, SN, ZI, ZR,
     * DCOSH, DSINH
      SH = DSINH(ZR)
      CH = DCOSH(ZR)
      SN = DSIN(ZI)
      CN = DCOS(ZI)
      CSHR = SH*CN
      CSHI = CH*SN
      CCHR = CH*CN
      CCHI = SH*SN
      RETURN
      END
      SUBROUTINE ZSSQRT(AR, AI, BR, BI)
C***BEGIN PROLOGUE  ZSSQRT
C***REFER TO  ZBESH,ZBESI,ZBESJ,ZBESK,ZBESY,ZAIRY,ZBIRY
C
C     DOUBLE PRECISION COMPLEX SQUARE ROOT, B=CSQRT(A)
C
C***ROUTINES CALLED  ZSABS
C***END PROLOGUE  ZSSQRT
      DOUBLE PRECISION AR, AI, BR, BI, ZM, DTHETA, DPI, DRT
      DOUBLE PRECISION ZSABS
      DATA DRT , DPI / 7.071067811865475244008443621D-1,
     1                 3.141592653589793238462643383D+0/
      ZM = ZSABS(AR,AI)
      ZM = DSQRT(ZM)
      IF (AR.EQ.0.0D+0) GO TO 10
      IF (AI.EQ.0.0D+0) GO TO 20
      DTHETA = DATAN(AI/AR)
      IF (DTHETA.LE.0.0D+0) GO TO 40
      IF (AR.LT.0.0D+0) DTHETA = DTHETA - DPI
      GO TO 50
   10 IF (AI.GT.0.0D+0) GO TO 60
      IF (AI.LT.0.0D+0) GO TO 70
      BR = 0.0D+0
      BI = 0.0D+0
      RETURN
   20 IF (AR.GT.0.0D+0) GO TO 30
      BR = 0.0D+0
      BI = DSQRT(DABS(AR))
      RETURN
   30 BR = DSQRT(AR)
      BI = 0.0D+0
      RETURN
   40 IF (AR.LT.0.0D+0) DTHETA = DTHETA + DPI
   50 DTHETA = DTHETA*0.5D+0
      BR = ZM*DCOS(DTHETA)
      BI = ZM*DSIN(DTHETA)
      RETURN
   60 BR = ZM*DRT
      BI = ZM*DRT
      RETURN
   70 BR = ZM*DRT
      BI = -ZM*DRT
      RETURN
      END
      SUBROUTINE ZUCHK(YR, YI, NZ, ASCLE, TOL)
C***BEGIN PROLOGUE  ZUCHK
C***REFER TO ZSERI,ZUOIK,ZUNK1,ZUNK2,ZUNI1,ZUNI2,ZKSCL
C
C      Y ENTERS AS A SCALED QUANTITY WHOSE MAGNITUDE IS GREATER THAN
C      EXP(-ALIM)=ASCLE=1.0E+3*D1MACH(1)/TOL. THE TEST IS MADE TO SEE
C      IF THE MAGNITUDE OF THE REAL OR IMAGINARY PART WOULD UNDERFLOW
C      WHEN Y IS SCALED (BY TOL) TO ITS PROPER VALUE. Y IS ACCEPTED
C      IF THE UNDERFLOW IS AT LEAST ONE PRECISION BELOW THE MAGNITUDE
C      OF THE LARGEST COMPONENT; OTHERWISE THE PHASE ANGLE DOES NOT HAVE
C      ABSOLUTE ACCURACY AND AN UNDERFLOW IS ASSUMED.
C
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  ZUCHK
C
C     COMPLEX Y
      DOUBLE PRECISION ASCLE, SS, ST, TOL, WR, WI, YR, YI
      INTEGER NZ
      NZ = 0
      WR = DABS(YR)
      WI = DABS(YI)
      ST = DMIN1(WR,WI)
      IF (ST.GT.ASCLE) RETURN
      SS = DMAX1(WR,WI)
      ST = ST/TOL
      IF (SS.LT.ST) NZ = 1
      RETURN
      END
      SUBROUTINE ZUNHJ(ZR, ZI, FNU, IPMTR, TOL, PHIR, PHII, ARGR, ARGI,
     * ZETA1R, ZETA1I, ZETA2R, ZETA2I, ASUMR, ASUMI, BSUMR, BSUMI)
C***BEGIN PROLOGUE  ZUNHJ
C***REFER TO  ZBESI,ZBESK
C
C     REFERENCES
C         HANDBOOK OF MATHEMATICAL FUNCTIONS BY M. ABRAMOWITZ AND I.A.
C         STEGUN, AMS55, NATIONAL BUREAU OF STANDARDS, 1965, CHAPTER 9.
C
C         ASYMPTOTICS AND SPECIAL FUNCTIONS BY F.W.J. OLVER, ACADEMIC
C         PRESS, N.Y., 1974, PAGE 420
C
C     ABSTRACT
C         ZUNHJ COMPUTES PARAMETERS FOR BESSEL FUNCTIONS C(FNU,Z) =
C         J(FNU,Z), Y(FNU,Z) OR H(I,FNU,Z) I=1,2 FOR LARGE ORDERS FNU
C         BY MEANS OF THE UNIFORM ASYMPTOTIC EXPANSION
C
C         C(FNU,Z)=C1*PHI*( ASUM*AIRY(ARG) + C2*BSUM*DAIRY(ARG) )
C
C         FOR PROPER CHOICES OF C1, C2, AIRY AND DAIRY WHERE AIRY IS
C         AN AIRY FUNCTION AND DAIRY IS ITS DERIVATIVE.
C
C               (2/3)*FNU*ZETA**1.5 = ZETA1-ZETA2,
C
C         ZETA1=0.5*FNU*CLOG((1+W)/(1-W)), ZETA2=FNU*W FOR SCALING
C         PURPOSES IN AIRY FUNCTIONS FROM CAIRY OR CBIRY.
C
C         MCONJ=SIGN OF AIMAG(Z), BUT IS AMBIGUOUS WHEN Z IS REAL AND
C         MUST BE SPECIFIED. IPMTR=0 RETURNS ALL PARAMETERS. IPMTR=
C         1 COMPUTES ALL EXCEPT ASUM AND BSUM.
C
C***ROUTINES CALLED  ZSABS,ZDIV,ZSLOG,ZSSQRT,D1MACH
C***END PROLOGUE  ZUNHJ
C     COMPLEX ARG,ASUM,BSUM,CFNU,CONE,CR,CZERO,DR,P,PHI,PRZTH,PTFN,
C    *RFN13,RTZTA,RZTH,SUMA,SUMB,TFN,T2,UP,W,W2,Z,ZA,ZB,ZC,ZETA,ZETA1,
C    *ZETA2,ZTH
      DOUBLE PRECISION ALFA, ANG, AP, AR, ARGI, ARGR, ASUMI, ASUMR,
     * ATOL, AW2, AZTH, BETA, BR, BSUMI, BSUMR, BTOL, C, CONEI, CONER,
     * CRI, CRR, DRI, DRR, EX1, EX2, FNU, FN13, FN23, GAMA, GPI, HPI,
     * PHII, PHIR, PI, PP, PR, PRZTHI, PRZTHR, PTFNI, PTFNR, RAW, RAW2,
     * RAZTH, RFNU, RFNU2, RFN13, RTZTI, RTZTR, RZTHI, RZTHR, STI, STR,
     * SUMAI, SUMAR, SUMBI, SUMBR, TEST, TFNI, TFNR, THPI, TOL, TZAI,
     * TZAR, T2I, T2R, UPI, UPR, WI, WR, W2I, W2R, ZAI, ZAR, ZBI, ZBR,
     * ZCI, ZCR, ZEROI, ZEROR, ZETAI, ZETAR, ZETA1I, ZETA1R, ZETA2I,
     * ZETA2R, ZI, ZR, ZTHI, ZTHR, ZSABS, AC, D1MACH
      INTEGER IAS, IBS, IPMTR, IS, J, JR, JU, K, KMAX, KP1, KS, L, LR,
     * LRP1, L1, L2, M, IDUM
      DIMENSION AR(14), BR(14), C(105), ALFA(180), BETA(210), GAMA(30),
     * AP(30), PR(30), PI(30), UPR(14), UPI(14), CRR(14), CRI(14),
     * DRR(14), DRI(14)
      DATA AR(1), AR(2), AR(3), AR(4), AR(5), AR(6), AR(7), AR(8),
     1     AR(9), AR(10), AR(11), AR(12), AR(13), AR(14)/
     2     1.00000000000000000D+00,     1.04166666666666667D-01,
     3     8.35503472222222222D-02,     1.28226574556327160D-01,
     4     2.91849026464140464D-01,     8.81627267443757652D-01,
     5     3.32140828186276754D+00,     1.49957629868625547D+01,
     6     7.89230130115865181D+01,     4.74451538868264323D+02,
     7     3.20749009089066193D+03,     2.40865496408740049D+04,
     8     1.98923119169509794D+05,     1.79190200777534383D+06/
      DATA BR(1), BR(2), BR(3), BR(4), BR(5), BR(6), BR(7), BR(8),
     1     BR(9), BR(10), BR(11), BR(12), BR(13), BR(14)/
     2     1.00000000000000000D+00,    -1.45833333333333333D-01,
     3    -9.87413194444444444D-02,    -1.43312053915895062D-01,
     4    -3.17227202678413548D-01,    -9.42429147957120249D-01,
     5    -3.51120304082635426D+00,    -1.57272636203680451D+01,
     6    -8.22814390971859444D+01,    -4.92355370523670524D+02,
     7    -3.31621856854797251D+03,    -2.48276742452085896D+04,
     8    -2.04526587315129788D+05,    -1.83844491706820990D+06/
      DATA C(1), C(2), C(3), C(4), C(5), C(6), C(7), C(8), C(9), C(10),
     1     C(11), C(12), C(13), C(14), C(15), C(16), C(17), C(18),
     2     C(19), C(20), C(21), C(22), C(23), C(24)/
     3     1.00000000000000000D+00,    -2.08333333333333333D-01,
     4     1.25000000000000000D-01,     3.34201388888888889D-01,
     5    -4.01041666666666667D-01,     7.03125000000000000D-02,
     6    -1.02581259645061728D+00,     1.84646267361111111D+00,
     7    -8.91210937500000000D-01,     7.32421875000000000D-02,
     8     4.66958442342624743D+00,    -1.12070026162229938D+01,
     9     8.78912353515625000D+00,    -2.36408691406250000D+00,
     A     1.12152099609375000D-01,    -2.82120725582002449D+01,
     B     8.46362176746007346D+01,    -9.18182415432400174D+01,
     C     4.25349987453884549D+01,    -7.36879435947963170D+00,
     D     2.27108001708984375D-01,     2.12570130039217123D+02,
     E    -7.65252468141181642D+02,     1.05999045252799988D+03/
      DATA C(25), C(26), C(27), C(28), C(29), C(30), C(31), C(32),
     1     C(33), C(34), C(35), C(36), C(37), C(38), C(39), C(40),
     2     C(41), C(42), C(43), C(44), C(45), C(46), C(47), C(48)/
     3    -6.99579627376132541D+02,     2.18190511744211590D+02,
     4    -2.64914304869515555D+01,     5.72501420974731445D-01,
     5    -1.91945766231840700D+03,     8.06172218173730938D+03,
     6    -1.35865500064341374D+04,     1.16553933368645332D+04,
     7    -5.30564697861340311D+03,     1.20090291321635246D+03,
     8    -1.08090919788394656D+02,     1.72772750258445740D+00,
     9     2.02042913309661486D+04,    -9.69805983886375135D+04,
     A     1.92547001232531532D+05,    -2.03400177280415534D+05,
     B     1.22200464983017460D+05,    -4.11926549688975513D+04,
     C     7.10951430248936372D+03,    -4.93915304773088012D+02,
     D     6.07404200127348304D+00,    -2.42919187900551333D+05,
     E     1.31176361466297720D+06,    -2.99801591853810675D+06/
      DATA C(49), C(50), C(51), C(52), C(53), C(54), C(55), C(56),
     1     C(57), C(58), C(59), C(60), C(61), C(62), C(63), C(64),
     2     C(65), C(66), C(67), C(68), C(69), C(70), C(71), C(72)/
     3     3.76327129765640400D+06,    -2.81356322658653411D+06,
     4     1.26836527332162478D+06,    -3.31645172484563578D+05,
     5     4.52187689813627263D+04,    -2.49983048181120962D+03,
     6     2.43805296995560639D+01,     3.28446985307203782D+06,
     7    -1.97068191184322269D+07,     5.09526024926646422D+07,
     8    -7.41051482115326577D+07,     6.63445122747290267D+07,
     9    -3.75671766607633513D+07,     1.32887671664218183D+07,
     A    -2.78561812808645469D+06,     3.08186404612662398D+05,
     B    -1.38860897537170405D+04,     1.10017140269246738D+02,
     C    -4.93292536645099620D+07,     3.25573074185765749D+08,
     D    -9.39462359681578403D+08,     1.55359689957058006D+09,
     E    -1.62108055210833708D+09,     1.10684281682301447D+09/
      DATA C(73), C(74), C(75), C(76), C(77), C(78), C(79), C(80),
     1     C(81), C(82), C(83), C(84), C(85), C(86), C(87), C(88),
     2     C(89), C(90), C(91), C(92), C(93), C(94), C(95), C(96)/
     3    -4.95889784275030309D+08,     1.42062907797533095D+08,
     4    -2.44740627257387285D+07,     2.24376817792244943D+06,
     5    -8.40054336030240853D+04,     5.51335896122020586D+02,
     6     8.14789096118312115D+08,    -5.86648149205184723D+09,
     7     1.86882075092958249D+10,    -3.46320433881587779D+10,
     8     4.12801855797539740D+10,    -3.30265997498007231D+10,
     9     1.79542137311556001D+10,    -6.56329379261928433D+09,
     A     1.55927986487925751D+09,    -2.25105661889415278D+08,
     B     1.73951075539781645D+07,    -5.49842327572288687D+05,
     C     3.03809051092238427D+03,    -1.46792612476956167D+10,
     D     1.14498237732025810D+11,    -3.99096175224466498D+11,
     E     8.19218669548577329D+11,    -1.09837515608122331D+12/
      DATA C(97), C(98), C(99), C(100), C(101), C(102), C(103), C(104),
     1     C(105)/
     2     1.00815810686538209D+12,    -6.45364869245376503D+11,
     3     2.87900649906150589D+11,    -8.78670721780232657D+10,
     4     1.76347306068349694D+10,    -2.16716498322379509D+09,
     5     1.43157876718888981D+08,    -3.87183344257261262D+06,
     6     1.82577554742931747D+04/
      DATA ALFA(1), ALFA(2), ALFA(3), ALFA(4), ALFA(5), ALFA(6),
     1     ALFA(7), ALFA(8), ALFA(9), ALFA(10), ALFA(11), ALFA(12),
     2     ALFA(13), ALFA(14), ALFA(15), ALFA(16), ALFA(17), ALFA(18),
     3     ALFA(19), ALFA(20), ALFA(21), ALFA(22)/
     4    -4.44444444444444444D-03,    -9.22077922077922078D-04,
     5    -8.84892884892884893D-05,     1.65927687832449737D-04,
     6     2.46691372741792910D-04,     2.65995589346254780D-04,
     7     2.61824297061500945D-04,     2.48730437344655609D-04,
     8     2.32721040083232098D-04,     2.16362485712365082D-04,
     9     2.00738858762752355D-04,     1.86267636637545172D-04,
     A     1.73060775917876493D-04,     1.61091705929015752D-04,
     B     1.50274774160908134D-04,     1.40503497391269794D-04,
     C     1.31668816545922806D-04,     1.23667445598253261D-04,
     D     1.16405271474737902D-04,     1.09798298372713369D-04,
     E     1.03772410422992823D-04,     9.82626078369363448D-05/
      DATA ALFA(23), ALFA(24), ALFA(25), ALFA(26), ALFA(27), ALFA(28),
     1     ALFA(29), ALFA(30), ALFA(31), ALFA(32), ALFA(33), ALFA(34),
     2     ALFA(35), ALFA(36), ALFA(37), ALFA(38), ALFA(39), ALFA(40),
     3     ALFA(41), ALFA(42), ALFA(43), ALFA(44)/
     4     9.32120517249503256D-05,     8.85710852478711718D-05,
     5     8.42963105715700223D-05,     8.03497548407791151D-05,
     6     7.66981345359207388D-05,     7.33122157481777809D-05,
     7     7.01662625163141333D-05,     6.72375633790160292D-05,
     8     6.93735541354588974D-04,     2.32241745182921654D-04,
     9    -1.41986273556691197D-05,    -1.16444931672048640D-04,
     A    -1.50803558053048762D-04,    -1.55121924918096223D-04,
     B    -1.46809756646465549D-04,    -1.33815503867491367D-04,
     C    -1.19744975684254051D-04,    -1.06184319207974020D-04,
     D    -9.37699549891194492D-05,    -8.26923045588193274D-05,
     E    -7.29374348155221211D-05,    -6.44042357721016283D-05/
      DATA ALFA(45), ALFA(46), ALFA(47), ALFA(48), ALFA(49), ALFA(50),
     1     ALFA(51), ALFA(52), ALFA(53), ALFA(54), ALFA(55), ALFA(56),
     2     ALFA(57), ALFA(58), ALFA(59), ALFA(60), ALFA(61), ALFA(62),
     3     ALFA(63), ALFA(64), ALFA(65), ALFA(66)/
     4    -5.69611566009369048D-05,    -5.04731044303561628D-05,
     5    -4.48134868008882786D-05,    -3.98688727717598864D-05,
     6    -3.55400532972042498D-05,    -3.17414256609022480D-05,
     7    -2.83996793904174811D-05,    -2.54522720634870566D-05,
     8    -2.28459297164724555D-05,    -2.05352753106480604D-05,
     9    -1.84816217627666085D-05,    -1.66519330021393806D-05,
     A    -1.50179412980119482D-05,    -1.35554031379040526D-05,
     B    -1.22434746473858131D-05,    -1.10641884811308169D-05,
     C    -3.54211971457743841D-04,    -1.56161263945159416D-04,
     D     3.04465503594936410D-05,     1.30198655773242693D-04,
     E     1.67471106699712269D-04,     1.70222587683592569D-04/
      DATA ALFA(67), ALFA(68), ALFA(69), ALFA(70), ALFA(71), ALFA(72),
     1     ALFA(73), ALFA(74), ALFA(75), ALFA(76), ALFA(77), ALFA(78),
     2     ALFA(79), ALFA(80), ALFA(81), ALFA(82), ALFA(83), ALFA(84),
     3     ALFA(85), ALFA(86), ALFA(87), ALFA(88)/
     4     1.56501427608594704D-04,     1.36339170977445120D-04,
     5     1.14886692029825128D-04,     9.45869093034688111D-05,
     6     7.64498419250898258D-05,     6.07570334965197354D-05,
     7     4.74394299290508799D-05,     3.62757512005344297D-05,
     8     2.69939714979224901D-05,     1.93210938247939253D-05,
     9     1.30056674793963203D-05,     7.82620866744496661D-06,
     A     3.59257485819351583D-06,     1.44040049814251817D-07,
     B    -2.65396769697939116D-06,    -4.91346867098485910D-06,
     C    -6.72739296091248287D-06,    -8.17269379678657923D-06,
     D    -9.31304715093561232D-06,    -1.02011418798016441D-05,
     E    -1.08805962510592880D-05,    -1.13875481509603555D-05/
      DATA ALFA(89), ALFA(90), ALFA(91), ALFA(92), ALFA(93), ALFA(94),
     1     ALFA(95), ALFA(96), ALFA(97), ALFA(98), ALFA(99), ALFA(100),
     2     ALFA(101), ALFA(102), ALFA(103), ALFA(104), ALFA(105),
     3     ALFA(106), ALFA(107), ALFA(108), ALFA(109), ALFA(110)/
     4    -1.17519675674556414D-05,    -1.19987364870944141D-05,
     5     3.78194199201772914D-04,     2.02471952761816167D-04,
     6    -6.37938506318862408D-05,    -2.38598230603005903D-04,
     7    -3.10916256027361568D-04,    -3.13680115247576316D-04,
     8    -2.78950273791323387D-04,    -2.28564082619141374D-04,
     9    -1.75245280340846749D-04,    -1.25544063060690348D-04,
     A    -8.22982872820208365D-05,    -4.62860730588116458D-05,
     B    -1.72334302366962267D-05,     5.60690482304602267D-06,
     C     2.31395443148286800D-05,     3.62642745856793957D-05,
     D     4.58006124490188752D-05,     5.24595294959114050D-05,
     E     5.68396208545815266D-05,     5.94349820393104052D-05/
      DATA ALFA(111), ALFA(112), ALFA(113), ALFA(114), ALFA(115),
     1     ALFA(116), ALFA(117), ALFA(118), ALFA(119), ALFA(120),
     2     ALFA(121), ALFA(122), ALFA(123), ALFA(124), ALFA(125),
     3     ALFA(126), ALFA(127), ALFA(128), ALFA(129), ALFA(130)/
     4     6.06478527578421742D-05,     6.08023907788436497D-05,
     5     6.01577894539460388D-05,     5.89199657344698500D-05,
     6     5.72515823777593053D-05,     5.52804375585852577D-05,
     7     5.31063773802880170D-05,     5.08069302012325706D-05,
     8     4.84418647620094842D-05,     4.60568581607475370D-05,
     9    -6.91141397288294174D-04,    -4.29976633058871912D-04,
     A     1.83067735980039018D-04,     6.60088147542014144D-04,
     B     8.75964969951185931D-04,     8.77335235958235514D-04,
     C     7.49369585378990637D-04,     5.63832329756980918D-04,
     D     3.68059319971443156D-04,     1.88464535514455599D-04/
      DATA ALFA(131), ALFA(132), ALFA(133), ALFA(134), ALFA(135),
     1     ALFA(136), ALFA(137), ALFA(138), ALFA(139), ALFA(140),
     2     ALFA(141), ALFA(142), ALFA(143), ALFA(144), ALFA(145),
     3     ALFA(146), ALFA(147), ALFA(148), ALFA(149), ALFA(150)/
     4     3.70663057664904149D-05,    -8.28520220232137023D-05,
     5    -1.72751952869172998D-04,    -2.36314873605872983D-04,
     6    -2.77966150694906658D-04,    -3.02079514155456919D-04,
     7    -3.12594712643820127D-04,    -3.12872558758067163D-04,
     8    -3.05678038466324377D-04,    -2.93226470614557331D-04,
     9    -2.77255655582934777D-04,    -2.59103928467031709D-04,
     A    -2.39784014396480342D-04,    -2.20048260045422848D-04,
     B    -2.00443911094971498D-04,    -1.81358692210970687D-04,
     C    -1.63057674478657464D-04,    -1.45712672175205844D-04,
     D    -1.29425421983924587D-04,    -1.14245691942445952D-04/
      DATA ALFA(151), ALFA(152), ALFA(153), ALFA(154), ALFA(155),
     1     ALFA(156), ALFA(157), ALFA(158), ALFA(159), ALFA(160),
     2     ALFA(161), ALFA(162), ALFA(163), ALFA(164), ALFA(165),
     3     ALFA(166), ALFA(167), ALFA(168), ALFA(169), ALFA(170)/
     4     1.92821964248775885D-03,     1.35592576302022234D-03,
     5    -7.17858090421302995D-04,    -2.58084802575270346D-03,
     6    -3.49271130826168475D-03,    -3.46986299340960628D-03,
     7    -2.82285233351310182D-03,    -1.88103076404891354D-03,
     8    -8.89531718383947600D-04,     3.87912102631035228D-06,
     9     7.28688540119691412D-04,     1.26566373053457758D-03,
     A     1.62518158372674427D-03,     1.83203153216373172D-03,
     B     1.91588388990527909D-03,     1.90588846755546138D-03,
     C     1.82798982421825727D-03,     1.70389506421121530D-03,
     D     1.55097127171097686D-03,     1.38261421852276159D-03/
      DATA ALFA(171), ALFA(172), ALFA(173), ALFA(174), ALFA(175),
     1     ALFA(176), ALFA(177), ALFA(178), ALFA(179), ALFA(180)/
     2     1.20881424230064774D-03,     1.03676532638344962D-03,
     3     8.71437918068619115D-04,     7.16080155297701002D-04,
     4     5.72637002558129372D-04,     4.42089819465802277D-04,
     5     3.24724948503090564D-04,     2.20342042730246599D-04,
     6     1.28412898401353882D-04,     4.82005924552095464D-05/
      DATA BETA(1), BETA(2), BETA(3), BETA(4), BETA(5), BETA(6),
     1     BETA(7), BETA(8), BETA(9), BETA(10), BETA(11), BETA(12),
     2     BETA(13), BETA(14), BETA(15), BETA(16), BETA(17), BETA(18),
     3     BETA(19), BETA(20), BETA(21), BETA(22)/
     4     1.79988721413553309D-02,     5.59964911064388073D-03,
     5     2.88501402231132779D-03,     1.80096606761053941D-03,
     6     1.24753110589199202D-03,     9.22878876572938311D-04,
     7     7.14430421727287357D-04,     5.71787281789704872D-04,
     8     4.69431007606481533D-04,     3.93232835462916638D-04,
     9     3.34818889318297664D-04,     2.88952148495751517D-04,
     A     2.52211615549573284D-04,     2.22280580798883327D-04,
     B     1.97541838033062524D-04,     1.76836855019718004D-04,
     C     1.59316899661821081D-04,     1.44347930197333986D-04,
     D     1.31448068119965379D-04,     1.20245444949302884D-04,
     E     1.10449144504599392D-04,     1.01828770740567258D-04/
      DATA BETA(23), BETA(24), BETA(25), BETA(26), BETA(27), BETA(28),
     1     BETA(29), BETA(30), BETA(31), BETA(32), BETA(33), BETA(34),
     2     BETA(35), BETA(36), BETA(37), BETA(38), BETA(39), BETA(40),
     3     BETA(41), BETA(42), BETA(43), BETA(44)/
     4     9.41998224204237509D-05,     8.74130545753834437D-05,
     5     8.13466262162801467D-05,     7.59002269646219339D-05,
     6     7.09906300634153481D-05,     6.65482874842468183D-05,
     7     6.25146958969275078D-05,     5.88403394426251749D-05,
     8    -1.49282953213429172D-03,    -8.78204709546389328D-04,
     9    -5.02916549572034614D-04,    -2.94822138512746025D-04,
     A    -1.75463996970782828D-04,    -1.04008550460816434D-04,
     B    -5.96141953046457895D-05,    -3.12038929076098340D-05,
     C    -1.26089735980230047D-05,    -2.42892608575730389D-07,
     D     8.05996165414273571D-06,     1.36507009262147391D-05,
     E     1.73964125472926261D-05,     1.98672978842133780D-05/
      DATA BETA(45), BETA(46), BETA(47), BETA(48), BETA(49), BETA(50),
     1     BETA(51), BETA(52), BETA(53), BETA(54), BETA(55), BETA(56),
     2     BETA(57), BETA(58), BETA(59), BETA(60), BETA(61), BETA(62),
     3     BETA(63), BETA(64), BETA(65), BETA(66)/
     4     2.14463263790822639D-05,     2.23954659232456514D-05,
     5     2.28967783814712629D-05,     2.30785389811177817D-05,
     6     2.30321976080909144D-05,     2.28236073720348722D-05,
     7     2.25005881105292418D-05,     2.20981015361991429D-05,
     8     2.16418427448103905D-05,     2.11507649256220843D-05,
     9     2.06388749782170737D-05,     2.01165241997081666D-05,
     A     1.95913450141179244D-05,     1.90689367910436740D-05,
     B     1.85533719641636667D-05,     1.80475722259674218D-05,
     C     5.52213076721292790D-04,     4.47932581552384646D-04,
     D     2.79520653992020589D-04,     1.52468156198446602D-04,
     E     6.93271105657043598D-05,     1.76258683069991397D-05/
      DATA BETA(67), BETA(68), BETA(69), BETA(70), BETA(71), BETA(72),
     1     BETA(73), BETA(74), BETA(75), BETA(76), BETA(77), BETA(78),
     2     BETA(79), BETA(80), BETA(81), BETA(82), BETA(83), BETA(84),
     3     BETA(85), BETA(86), BETA(87), BETA(88)/
     4    -1.35744996343269136D-05,    -3.17972413350427135D-05,
     5    -4.18861861696693365D-05,    -4.69004889379141029D-05,
     6    -4.87665447413787352D-05,    -4.87010031186735069D-05,
     7    -4.74755620890086638D-05,    -4.55813058138628452D-05,
     8    -4.33309644511266036D-05,    -4.09230193157750364D-05,
     9    -3.84822638603221274D-05,    -3.60857167535410501D-05,
     A    -3.37793306123367417D-05,    -3.15888560772109621D-05,
     B    -2.95269561750807315D-05,    -2.75978914828335759D-05,
     C    -2.58006174666883713D-05,    -2.41308356761280200D-05,
     D    -2.25823509518346033D-05,    -2.11479656768912971D-05,
     E    -1.98200638885294927D-05,    -1.85909870801065077D-05/
      DATA BETA(89), BETA(90), BETA(91), BETA(92), BETA(93), BETA(94),
     1     BETA(95), BETA(96), BETA(97), BETA(98), BETA(99), BETA(100),
     2     BETA(101), BETA(102), BETA(103), BETA(104), BETA(105),
     3     BETA(106), BETA(107), BETA(108), BETA(109), BETA(110)/
     4    -1.74532699844210224D-05,    -1.63997823854497997D-05,
     5    -4.74617796559959808D-04,    -4.77864567147321487D-04,
     6    -3.20390228067037603D-04,    -1.61105016119962282D-04,
     7    -4.25778101285435204D-05,     3.44571294294967503D-05,
     8     7.97092684075674924D-05,     1.03138236708272200D-04,
     9     1.12466775262204158D-04,     1.13103642108481389D-04,
     A     1.08651634848774268D-04,     1.01437951597661973D-04,
     B     9.29298396593363896D-05,     8.40293133016089978D-05,
     C     7.52727991349134062D-05,     6.69632521975730872D-05,
     D     5.92564547323194704D-05,     5.22169308826975567D-05,
     E     4.58539485165360646D-05,     4.01445513891486808D-05/
      DATA BETA(111), BETA(112), BETA(113), BETA(114), BETA(115),
     1     BETA(116), BETA(117), BETA(118), BETA(119), BETA(120),
     2     BETA(121), BETA(122), BETA(123), BETA(124), BETA(125),
     3     BETA(126), BETA(127), BETA(128), BETA(129), BETA(130)/
     4     3.50481730031328081D-05,     3.05157995034346659D-05,
     5     2.64956119950516039D-05,     2.29363633690998152D-05,
     6     1.97893056664021636D-05,     1.70091984636412623D-05,
     7     1.45547428261524004D-05,     1.23886640995878413D-05,
     8     1.04775876076583236D-05,     8.79179954978479373D-06,
     9     7.36465810572578444D-04,     8.72790805146193976D-04,
     A     6.22614862573135066D-04,     2.85998154194304147D-04,
     B     3.84737672879366102D-06,    -1.87906003636971558D-04,
     C    -2.97603646594554535D-04,    -3.45998126832656348D-04,
     D    -3.53382470916037712D-04,    -3.35715635775048757D-04/
      DATA BETA(131), BETA(132), BETA(133), BETA(134), BETA(135),
     1     BETA(136), BETA(137), BETA(138), BETA(139), BETA(140),
     2     BETA(141), BETA(142), BETA(143), BETA(144), BETA(145),
     3     BETA(146), BETA(147), BETA(148), BETA(149), BETA(150)/
     4    -3.04321124789039809D-04,    -2.66722723047612821D-04,
     5    -2.27654214122819527D-04,    -1.89922611854562356D-04,
     6    -1.55058918599093870D-04,    -1.23778240761873630D-04,
     7    -9.62926147717644187D-05,    -7.25178327714425337D-05,
     8    -5.22070028895633801D-05,    -3.50347750511900522D-05,
     9    -2.06489761035551757D-05,    -8.70106096849767054D-06,
     A     1.13698686675100290D-06,     9.16426474122778849D-06,
     B     1.56477785428872620D-05,     2.08223629482466847D-05,
     C     2.48923381004595156D-05,     2.80340509574146325D-05,
     D     3.03987774629861915D-05,     3.21156731406700616D-05/
      DATA BETA(151), BETA(152), BETA(153), BETA(154), BETA(155),
     1     BETA(156), BETA(157), BETA(158), BETA(159), BETA(160),
     2     BETA(161), BETA(162), BETA(163), BETA(164), BETA(165),
     3     BETA(166), BETA(167), BETA(168), BETA(169), BETA(170)/
     4    -1.80182191963885708D-03,    -2.43402962938042533D-03,
     5    -1.83422663549856802D-03,    -7.62204596354009765D-04,
     6     2.39079475256927218D-04,     9.49266117176881141D-04,
     7     1.34467449701540359D-03,     1.48457495259449178D-03,
     8     1.44732339830617591D-03,     1.30268261285657186D-03,
     9     1.10351597375642682D-03,     8.86047440419791759D-04,
     A     6.73073208165665473D-04,     4.77603872856582378D-04,
     B     3.05991926358789362D-04,     1.60315694594721630D-04,
     C     4.00749555270613286D-05,    -5.66607461635251611D-05,
     D    -1.32506186772982638D-04,    -1.90296187989614057D-04/
      DATA BETA(171), BETA(172), BETA(173), BETA(174), BETA(175),
     1     BETA(176), BETA(177), BETA(178), BETA(179), BETA(180),
     2     BETA(181), BETA(182), BETA(183), BETA(184), BETA(185),
     3     BETA(186), BETA(187), BETA(188), BETA(189), BETA(190)/
     4    -2.32811450376937408D-04,    -2.62628811464668841D-04,
     5    -2.82050469867598672D-04,    -2.93081563192861167D-04,
     6    -2.97435962176316616D-04,    -2.96557334239348078D-04,
     7    -2.91647363312090861D-04,    -2.83696203837734166D-04,
     8    -2.73512317095673346D-04,    -2.61750155806768580D-04,
     9     6.38585891212050914D-03,     9.62374215806377941D-03,
     A     7.61878061207001043D-03,     2.83219055545628054D-03,
     B    -2.09841352012720090D-03,    -5.73826764216626498D-03,
     C    -7.70804244495414620D-03,    -8.21011692264844401D-03,
     D    -7.65824520346905413D-03,    -6.47209729391045177D-03/
      DATA BETA(191), BETA(192), BETA(193), BETA(194), BETA(195),
     1     BETA(196), BETA(197), BETA(198), BETA(199), BETA(200),
     2     BETA(201), BETA(202), BETA(203), BETA(204), BETA(205),
     3     BETA(206), BETA(207), BETA(208), BETA(209), BETA(210)/
     4    -4.99132412004966473D-03,    -3.45612289713133280D-03,
     5    -2.01785580014170775D-03,    -7.59430686781961401D-04,
     6     2.84173631523859138D-04,     1.10891667586337403D-03,
     7     1.72901493872728771D-03,     2.16812590802684701D-03,
     8     2.45357710494539735D-03,     2.61281821058334862D-03,
     9     2.67141039656276912D-03,     2.65203073395980430D-03,
     A     2.57411652877287315D-03,     2.45389126236094427D-03,
     B     2.30460058071795494D-03,     2.13684837686712662D-03,
     C     1.95896528478870911D-03,     1.77737008679454412D-03,
     D     1.59690280765839059D-03,     1.42111975664438546D-03/
      DATA GAMA(1), GAMA(2), GAMA(3), GAMA(4), GAMA(5), GAMA(6),
     1     GAMA(7), GAMA(8), GAMA(9), GAMA(10), GAMA(11), GAMA(12),
     2     GAMA(13), GAMA(14), GAMA(15), GAMA(16), GAMA(17), GAMA(18),
     3     GAMA(19), GAMA(20), GAMA(21), GAMA(22)/
     4     6.29960524947436582D-01,     2.51984209978974633D-01,
     5     1.54790300415655846D-01,     1.10713062416159013D-01,
     6     8.57309395527394825D-02,     6.97161316958684292D-02,
     7     5.86085671893713576D-02,     5.04698873536310685D-02,
     8     4.42600580689154809D-02,     3.93720661543509966D-02,
     9     3.54283195924455368D-02,     3.21818857502098231D-02,
     A     2.94646240791157679D-02,     2.71581677112934479D-02,
     B     2.51768272973861779D-02,     2.34570755306078891D-02,
     C     2.19508390134907203D-02,     2.06210828235646240D-02,
     D     1.94388240897880846D-02,     1.83810633800683158D-02,
     E     1.74293213231963172D-02,     1.65685837786612353D-02/
      DATA GAMA(23), GAMA(24), GAMA(25), GAMA(26), GAMA(27), GAMA(28),
     1     GAMA(29), GAMA(30)/
     2     1.57865285987918445D-02,     1.50729501494095594D-02,
     3     1.44193250839954639D-02,     1.38184805735341786D-02,
     4     1.32643378994276568D-02,     1.27517121970498651D-02,
     5     1.22761545318762767D-02,     1.18338262398482403D-02/
      DATA EX1, EX2, HPI, GPI, THPI /
     1     3.33333333333333333D-01,     6.66666666666666667D-01,
     2     1.57079632679489662D+00,     3.14159265358979324D+00,
     3     4.71238898038468986D+00/
      DATA ZEROR,ZEROI,CONER,CONEI / 0.0D0, 0.0D0, 1.0D0, 0.0D0 /
C
      RFNU = 1.0D0/FNU
C-----------------------------------------------------------------------
C     OVERFLOW TEST (Z/FNU TOO SMALL)
C-----------------------------------------------------------------------
      TEST = D1MACH(1)*1.0D+3
      AC = FNU*TEST
      IF (DABS(ZR).GT.AC .OR. DABS(ZI).GT.AC) GO TO 15
      ZETA1R = 2.0D0*DABS(DLOG(TEST))+FNU
      ZETA1I = 0.0D0
      ZETA2R = FNU
      ZETA2I = 0.0D0
      PHIR = 1.0D0
      PHII = 0.0D0
      ARGR = 1.0D0
      ARGI = 0.0D0
      RETURN
   15 CONTINUE
      ZBR = ZR*RFNU
      ZBI = ZI*RFNU
      RFNU2 = RFNU*RFNU
C-----------------------------------------------------------------------
C     COMPUTE IN THE FOURTH QUADRANT
C-----------------------------------------------------------------------
      FN13 = FNU**EX1
      FN23 = FN13*FN13
      RFN13 = 1.0D0/FN13
      W2R = CONER - ZBR*ZBR + ZBI*ZBI
      W2I = CONEI - ZBR*ZBI - ZBR*ZBI
      AW2 = ZSABS(W2R,W2I)
      IF (AW2.GT.0.25D0) GO TO 130
C-----------------------------------------------------------------------
C     POWER SERIES FOR CABS(W2).LE.0.25D0
C-----------------------------------------------------------------------
      K = 1
      PR(1) = CONER
      PI(1) = CONEI
      SUMAR = GAMA(1)
      SUMAI = ZEROI
      AP(1) = 1.0D0
      IF (AW2.LT.TOL) GO TO 20
      DO 10 K=2,30
        PR(K) = PR(K-1)*W2R - PI(K-1)*W2I
        PI(K) = PR(K-1)*W2I + PI(K-1)*W2R
        SUMAR = SUMAR + PR(K)*GAMA(K)
        SUMAI = SUMAI + PI(K)*GAMA(K)
        AP(K) = AP(K-1)*AW2
        IF (AP(K).LT.TOL) GO TO 20
   10 CONTINUE
      K = 30
   20 CONTINUE
      KMAX = K
      ZETAR = W2R*SUMAR - W2I*SUMAI
      ZETAI = W2R*SUMAI + W2I*SUMAR
      ARGR = ZETAR*FN23
      ARGI = ZETAI*FN23
      CALL ZSSQRT(SUMAR, SUMAI, ZAR, ZAI)
      CALL ZSSQRT(W2R, W2I, STR, STI)
      ZETA2R = STR*FNU
      ZETA2I = STI*FNU
      STR = CONER + EX2*(ZETAR*ZAR-ZETAI*ZAI)
      STI = CONEI + EX2*(ZETAR*ZAI+ZETAI*ZAR)
      ZETA1R = STR*ZETA2R - STI*ZETA2I
      ZETA1I = STR*ZETA2I + STI*ZETA2R
      ZAR = ZAR + ZAR
      ZAI = ZAI + ZAI
      CALL ZSSQRT(ZAR, ZAI, STR, STI)
      PHIR = STR*RFN13
      PHII = STI*RFN13
      IF (IPMTR.EQ.1) GO TO 120
C-----------------------------------------------------------------------
C     SUM SERIES FOR ASUM AND BSUM
C-----------------------------------------------------------------------
      SUMBR = ZEROR
      SUMBI = ZEROI
      DO 30 K=1,KMAX
        SUMBR = SUMBR + PR(K)*BETA(K)
        SUMBI = SUMBI + PI(K)*BETA(K)
   30 CONTINUE
      ASUMR = ZEROR
      ASUMI = ZEROI
      BSUMR = SUMBR
      BSUMI = SUMBI
      L1 = 0
      L2 = 30
      BTOL = TOL*(DABS(BSUMR)+DABS(BSUMI))
      ATOL = TOL
      PP = 1.0D0
      IAS = 0
      IBS = 0
      IF (RFNU2.LT.TOL) GO TO 110
      DO 100 IS=2,7
        ATOL = ATOL/RFNU2
        PP = PP*RFNU2
        IF (IAS.EQ.1) GO TO 60
        SUMAR = ZEROR
        SUMAI = ZEROI
        DO 40 K=1,KMAX
          M = L1 + K
          SUMAR = SUMAR + PR(K)*ALFA(M)
          SUMAI = SUMAI + PI(K)*ALFA(M)
          IF (AP(K).LT.ATOL) GO TO 50
   40   CONTINUE
   50   CONTINUE
        ASUMR = ASUMR + SUMAR*PP
        ASUMI = ASUMI + SUMAI*PP
        IF (PP.LT.TOL) IAS = 1
   60   CONTINUE
        IF (IBS.EQ.1) GO TO 90
        SUMBR = ZEROR
        SUMBI = ZEROI
        DO 70 K=1,KMAX
          M = L2 + K
          SUMBR = SUMBR + PR(K)*BETA(M)
          SUMBI = SUMBI + PI(K)*BETA(M)
          IF (AP(K).LT.ATOL) GO TO 80
   70   CONTINUE
   80   CONTINUE
        BSUMR = BSUMR + SUMBR*PP
        BSUMI = BSUMI + SUMBI*PP
        IF (PP.LT.BTOL) IBS = 1
   90   CONTINUE
        IF (IAS.EQ.1 .AND. IBS.EQ.1) GO TO 110
        L1 = L1 + 30
        L2 = L2 + 30
  100 CONTINUE
  110 CONTINUE
      ASUMR = ASUMR + CONER
      PP = RFNU*RFN13
      BSUMR = BSUMR*PP
      BSUMI = BSUMI*PP
  120 CONTINUE
      RETURN
C-----------------------------------------------------------------------
C     CABS(W2).GT.0.25D0
C-----------------------------------------------------------------------
  130 CONTINUE
      CALL ZSSQRT(W2R, W2I, WR, WI)
      IF (WR.LT.0.0D0) WR = 0.0D0
      IF (WI.LT.0.0D0) WI = 0.0D0
      STR = CONER + WR
      STI = WI
      CALL ZDIV(STR, STI, ZBR, ZBI, ZAR, ZAI)
      CALL ZSLOG(ZAR, ZAI, ZCR, ZCI, IDUM)
      IF (ZCI.LT.0.0D0) ZCI = 0.0D0
      IF (ZCI.GT.HPI) ZCI = HPI
      IF (ZCR.LT.0.0D0) ZCR = 0.0D0
      ZTHR = (ZCR-WR)*1.5D0
      ZTHI = (ZCI-WI)*1.5D0
      ZETA1R = ZCR*FNU
      ZETA1I = ZCI*FNU
      ZETA2R = WR*FNU
      ZETA2I = WI*FNU
      AZTH = ZSABS(ZTHR,ZTHI)
      ANG = THPI
      IF (ZTHR.GE.0.0D0 .AND. ZTHI.LT.0.0D0) GO TO 140
      ANG = HPI
      IF (ZTHR.EQ.0.0D0) GO TO 140
      ANG = DATAN(ZTHI/ZTHR)
      IF (ZTHR.LT.0.0D0) ANG = ANG + GPI
  140 CONTINUE
      PP = AZTH**EX2
      ANG = ANG*EX2
      ZETAR = PP*DCOS(ANG)
      ZETAI = PP*DSIN(ANG)
      IF (ZETAI.LT.0.0D0) ZETAI = 0.0D0
      ARGR = ZETAR*FN23
      ARGI = ZETAI*FN23
      CALL ZDIV(ZTHR, ZTHI, ZETAR, ZETAI, RTZTR, RTZTI)
      CALL ZDIV(RTZTR, RTZTI, WR, WI, ZAR, ZAI)
      TZAR = ZAR + ZAR
      TZAI = ZAI + ZAI
      CALL ZSSQRT(TZAR, TZAI, STR, STI)
      PHIR = STR*RFN13
      PHII = STI*RFN13
      IF (IPMTR.EQ.1) GO TO 120
      RAW = 1.0D0/DSQRT(AW2)
      STR = WR*RAW
      STI = -WI*RAW
      TFNR = STR*RFNU*RAW
      TFNI = STI*RFNU*RAW
      RAZTH = 1.0D0/AZTH
      STR = ZTHR*RAZTH
      STI = -ZTHI*RAZTH
      RZTHR = STR*RAZTH*RFNU
      RZTHI = STI*RAZTH*RFNU
      ZCR = RZTHR*AR(2)
      ZCI = RZTHI*AR(2)
      RAW2 = 1.0D0/AW2
      STR = W2R*RAW2
      STI = -W2I*RAW2
      T2R = STR*RAW2
      T2I = STI*RAW2
      STR = T2R*C(2) + C(3)
      STI = T2I*C(2)
      UPR(2) = STR*TFNR - STI*TFNI
      UPI(2) = STR*TFNI + STI*TFNR
      BSUMR = UPR(2) + ZCR
      BSUMI = UPI(2) + ZCI
      ASUMR = ZEROR
      ASUMI = ZEROI
      IF (RFNU.LT.TOL) GO TO 220
      PRZTHR = RZTHR
      PRZTHI = RZTHI
      PTFNR = TFNR
      PTFNI = TFNI
      UPR(1) = CONER
      UPI(1) = CONEI
      PP = 1.0D0
      BTOL = TOL*(DABS(BSUMR)+DABS(BSUMI))
      KS = 0
      KP1 = 2
      L = 3
      IAS = 0
      IBS = 0
      DO 210 LR=2,12,2
        LRP1 = LR + 1
C-----------------------------------------------------------------------
C     COMPUTE TWO ADDITIONAL CR, DR, AND UP FOR TWO MORE TERMS IN
C     NEXT SUMA AND SUMB
C-----------------------------------------------------------------------
        DO 160 K=LR,LRP1
          KS = KS + 1
          KP1 = KP1 + 1
          L = L + 1
          ZAR = C(L)
          ZAI = ZEROI
          DO 150 J=2,KP1
            L = L + 1
            STR = ZAR*T2R - T2I*ZAI + C(L)
            ZAI = ZAR*T2I + ZAI*T2R
            ZAR = STR
  150     CONTINUE
          STR = PTFNR*TFNR - PTFNI*TFNI
          PTFNI = PTFNR*TFNI + PTFNI*TFNR
          PTFNR = STR
          UPR(KP1) = PTFNR*ZAR - PTFNI*ZAI
          UPI(KP1) = PTFNI*ZAR + PTFNR*ZAI
          CRR(KS) = PRZTHR*BR(KS+1)
          CRI(KS) = PRZTHI*BR(KS+1)
          STR = PRZTHR*RZTHR - PRZTHI*RZTHI
          PRZTHI = PRZTHR*RZTHI + PRZTHI*RZTHR
          PRZTHR = STR
          DRR(KS) = PRZTHR*AR(KS+2)
          DRI(KS) = PRZTHI*AR(KS+2)
  160   CONTINUE
        PP = PP*RFNU2
        IF (IAS.EQ.1) GO TO 180
        SUMAR = UPR(LRP1)
        SUMAI = UPI(LRP1)
        JU = LRP1
        DO 170 JR=1,LR
          JU = JU - 1
          SUMAR = SUMAR + CRR(JR)*UPR(JU) - CRI(JR)*UPI(JU)
          SUMAI = SUMAI + CRR(JR)*UPI(JU) + CRI(JR)*UPR(JU)
  170   CONTINUE
        ASUMR = ASUMR + SUMAR
        ASUMI = ASUMI + SUMAI
        TEST = DABS(SUMAR) + DABS(SUMAI)
        IF (PP.LT.TOL .AND. TEST.LT.TOL) IAS = 1
  180   CONTINUE
        IF (IBS.EQ.1) GO TO 200
        SUMBR = UPR(LR+2) + UPR(LRP1)*ZCR - UPI(LRP1)*ZCI
        SUMBI = UPI(LR+2) + UPR(LRP1)*ZCI + UPI(LRP1)*ZCR
        JU = LRP1
        DO 190 JR=1,LR
          JU = JU - 1
          SUMBR = SUMBR + DRR(JR)*UPR(JU) - DRI(JR)*UPI(JU)
          SUMBI = SUMBI + DRR(JR)*UPI(JU) + DRI(JR)*UPR(JU)
  190   CONTINUE
        BSUMR = BSUMR + SUMBR
        BSUMI = BSUMI + SUMBI
        TEST = DABS(SUMBR) + DABS(SUMBI)
        IF (PP.LT.BTOL .AND. TEST.LT.BTOL) IBS = 1
  200   CONTINUE
        IF (IAS.EQ.1 .AND. IBS.EQ.1) GO TO 220
  210 CONTINUE
  220 CONTINUE
      ASUMR = ASUMR + CONER
      STR = -BSUMR*RFN13
      STI = -BSUMI*RFN13
      CALL ZDIV(STR, STI, RTZTR, RTZTI, BSUMR, BSUMI)
      GO TO 120
      END
      
      SUBROUTINE ZGER(A,INT,N,NC,EMACH)
C     ------------------------------------------------------------------
C     ZGER IS A STANDARD SUBROUTINE TO PERFORM GAUSSIAN ELIMINATION ON
C     A NC*NC MATRIX 'A' PRIOR TO INVERSION FOR THE MATRIX PROBLEM
C     WHERE THE MATRIX A MULTIPLIES A VECTOR X ON THE RIGHT:
C
C                         X*A=B
C
C     ZGER MAKES AN LOWER DIAGONAL MATRIX AND HAS ONLY 23 PROGRAM LINES!
C     THIS ROUTINE DOES NOT BOTHER ABOUT ELEMENTS DIRECTLY ABOVE
C     THE MATRIX DIAGONAL AS THEY ARE NOT USED EXPLICITLY IN AN
C     ACCOMAPANYING ZSER ROUTINE.
C
C     INT   RECORDS PIVOTING DETAILS OF THE GAUSS ELIMINATION.
C     EMACH IS A CUTOFF ON THE MATRIX ELEMENTS
C     ------------------------------------------------------------------
      IMPLICIT NONE
C
C ..  SCALAR ARGUMENTS  ..
C
      INTEGER N,NC
      REAL*8 EMACH
C
C ..  ARRAY ARGUMENTS  ..
C
      INTEGER    INT(NC)
      COMPLEX*16 A(NC,NC)
C
C ..  LOCAL SCALARS  ..
C
      INTEGER    I,II,IN,J,K
      COMPLEX*16 YR,DUM
C
C ..  INTRINSIC FUNCTIONS  ..
C
*      INTRINSIC ABS
C     ------------------------------------------------------------------
C
      DO 10 II=2,N
      I=II-1          !from 1 to N-1
      YR=A(I,I)
      IN=I
*
* PIVOTING:
* Finding an element with the largest magnitude in the I-th row
* to the right of the matrix diagonal (I,I)-element
* (including the diag. element):

      DO 2 J=II,N    !J from I+1 to N
      IF(ABS(YR)-ABS(A(I,J)))1,2,2
   1  YR=A(I,J)
      IN=J
   2  CONTINUE
      INT(I)=IN      !The largest element in the I-th row to the right of
                     !the matrix diagonal (I,I)-element is
                           !in the IN-th column and is denoted by YR
*
* Executing pivoting:

      IF(IN-I)3,5,3  !If IN.NE.I exchange the I-th and IN-th columns on the
   3  DO 4 J=I,N     !right of and including the matrix diagonal (I,I)-element
      DUM=A(J,I)     !Only column parts below and including I-th component
      A(J,I)=A(J,IN) !are exchanged - the remaining are ignored.
   4  A(J,IN)=DUM

   5  IF(ABS(YR)-EMACH)10,10,6

* For J from I+1 to N, the A(I,J)/A(I,I) multiple of
* the Ith column is subtracted from the Jth column,
* as if it were the Gaussian elimination of matrix elements
* in the Ith row to the right of the diagonal (I,I)-element,
* but with the following differences:
* 1) only the elements of the Jth column below the Ith element
*    are affected by the subtraction,
* and
* 2) the subtraction is performed if and only if both
*    A(I,I) and A(I,J) are sufficiently large!!!
*
* 3) Additionally, A(I,I) is not rescaled to 1
*
   6  DO 9 J=II,N                   !J from I+1 to N

      IF(ABS(A(I,J))-EMACH)9,9,7

   7  A(I,J)=A(I,J)/YR

      DO 8 K=II,N
   8  A(K,J)=A(K,J)-A(I,J)*A(K,I)   !In Jth column from I+1 to N,
   9  CONTINUE                      !the elements above the (I+1)th
                                    !element are not affected
*
  10  CONTINUE                      !end of "column loop"
      RETURN
      END
      SUBROUTINE znsrefind(LAMBDA,FILFRAC,zeps)
C--------/---------/---------/---------/---------/---------/---------/--
C   FILFRAC ... filfrac of ZnS in ZnS core
C--------/---------/---------/---------/---------/---------/---------/--
      IMPLICIT NONE
      COMPLEX*16              zeps
      REAL*8          F,FILFRAC,LAMBDA,EPSHOST
      REAL*8          EPSZnS
      REAL*8          eMG1,EPSPAR
*
c if the host is different from the medium (silica n=1.45)

      EPSHOST = 1.45d0*1.45d0

c ZnS filling fraction f

      f = FILFRAC
c Particle material: ZnS bulk dielectric constant 350 nm - 900 nm

      EPSZnS = 5.164d0 + 1.208d+7/(lambda*lambda*100 - 0.732d+7)

* particle material

      EPSPAR = EPSZnS

c               Bruggeman (effective medium host)
c               wor = (3.d0*F-1.d0)*EPSPAR+(2.d0 - 3.d0*F)*EPSHOST
c               wor1 =  sqrt(wor*wor + 8.d0*EPSPAR*EPSHOST)
c       if (AIMAG(wor1).GE.0.0) then
c               eBr = (wor + wor1 )/4.d0
c       else
c               e_Br = (wor - wor1 )/4.d0
c       end if
c               Maxwell-Garnett MG1    (medium material host)

      eMG1=EPSHOST*(2.d0*EPSHOST+EPSPAR + 2.d0*F*(EPSPAR - EPSHOST))/
     &            (2.d0*EPSHOST + EPSPAR - F*(EPSPAR - EPSHOST))
C--------/---------/---------/---------/---------/---------/---------/--
c               Maxwell-Garnett MG2 (particle material host)
C       eMG2 = EPSPAR*(2.*EPSPAR + EPSHOST + 2.*(1.- F)*(EPSHOST - EPSPAR))/
C     &          (2.*EPSPAR + EPSHOST -(1.- F)*(EPSHOST - EPSPAR))
c      WRITE(*,*) LAMBDA, eMG1
*
      zeps = eMG1
*
      RETURN
      END

      SUBROUTINE ZSUR(A,INT,X,N,NC,EMACH)
C     ------------------------------------------------------------------
C     ZSUR IS A STANDARD BACK-SUBSTITUTION SUBROUTINE USING THE
C     OUTPUT OF ZGER TO CALCULATE X TIMES A-INVERSE, RETURNED IN X.
C     IT HAS ONLY 21 PROGRAM LINES.
C     FOR THE MATRIX PROBLEM WHERE THE MATRIX 'A' MULTIPLIES A
C     VECTOR X ON THE RIGHT:
C
C                         X*A=B
C
C     INT RECORDS PIVOTING DETAILS OF THE GAUSS ELIMINATION
C         PERFOMED IN ZGER.
C     EMACH IS A CUTOFF ON THE MATRIX ELEMENTS.
C     ------------------------------------------------------------------
      IMPLICIT NONE
C
C ..  SCALAR ARGUMENTS  ..
C
      INTEGER N,NC
      REAL*8 EMACH
C
C ..  ARRAY ARGUMENTS  ..
C
      INTEGER    INT(NC)
      COMPLEX*16 A(NC,NC),X(NC)
C
C ..  LOCAL SCALARS  ..
C
      INTEGER    I,II,IN,J,IJ
      COMPLEX*16 DUM
C
C ..  INTRINSIC FUNCTIONS  ..
C
*      INTRINSIC ABS
C     ------------------------------------------------------------------
C
C Rearanging the "right-hand side" using pivoting information from
C ZGER:

      DO 5 II=2,N
      I=II-1               !from 1 to N-1
      IF(INT(I)-I)1,2,1    !If INT(I).NE.I exchange the I-th and INT(I)-th
                           !elements of the vector X
   1  IN=INT(I)
      DUM=X(IN)
      X(IN)=X(I)
      X(I)=DUM
*
C Rearanging the "right-hand side" by mirroring the alegbaric
C operations performed on columns of 'A' matrix in ZGER:
* If A(I,J) is sufficiently large, X(J)=X(J)-X(I)*A(I,J)
*
   2  DO 4 J=II,N                 !from I+1 to N
      IF(ABS(A(I,J))-EMACH)4,4,3
   3  X(J)=X(J)-X(I)*A(I,J)
   4  CONTINUE

   5  CONTINUE               !the I-th row of A multiplied by X(I)
                             !subtracted from X
*
* BACKSUBSTITUTION:
*
      DO 10 II=1,N
      I=N-II+1         !from N   to 1
      IJ=I+1           !from N+1 to 2
      IF(I-N)6,8,6
   6  DO 7 J=IJ,N
   7  X(I)=X(I)-X(J)*A(J,I)

* Imposing a cutoff on the possible underflow.
* You can change the factor 1.0D-7 below to some
* other number, depending on the precision of your
* computer

   8  IF(ABS(A(I,I))-EMACH*1.0D-7)9,10,10
   9  A(I,I)=EMACH*1.0D-7*(1.D0,1.D0)
  10  X(I)=X(I)/A(I,I)
      RETURN
      END

C (C) Copr. 03/2003  Alexander Moroz

C (C) Copr. 09/1998  Alexander Moroz

C (C) Copr. 03/2003  Alexander Moroz

C (C) Copr. 10/2005  Alexander Moroz
C (C) Copr. 01/2003  Alexander Moroz
C (C) Copr. 01/2003  Alexander Moroz
      function zartan(zs)
*--------/---------/---------/---------/---------/---------/---------/--
* For a complex argument using that
*              atan (z)= -ln((1+iz)/(1-iz))*i/2.d0
* Hower, a direct application of this formula often
* results in a nonzero imaginary part of  $\arctan z$ even
* for a  purely real $z=x$.
* Therefore, in order to avoid this, in a numerical
* implementation, the complex logarithm is written
* explicitly as
*
*  ln\left(\fr{1+iz}{1-iz}\right)= (ar1,ar2)= log(
*  (1.d0 + abs(z)**2 -2.d0*imag(z))/1.d0+abs(z)**2 +2.d0*imag(z)))/2.d0
* +
* ci*atan(2.d0*dble(z)/(1-abs(z)**2).
*
* For a real z=x, |x|<1,  log here is purely imaginary and
* equals to i*atan(2x/(1-x**2))\equiv 2*i*atan x
*--------/---------/---------/---------/---------/---------/---------/--
      complex*16 zartan,zs
      real*8 xxs,xas,ar1,ar2,pi
      DATA PI/3.141592653589793d0/

      if (dimag(zs).eq.0.) then
        zartan=dcmplx(atan(dble(zs)),0.d0)
      return
      end if
         xas=abs(zs)**2
      ar1=log((1.d0+xas-2.d0*dimag(zs))/(1.d0+xas+2.d0*dimag(zs)))/2.d0
         xxs=dble(zs)
* special case:
         if(xas.eq.1.) then
           if(xxs.ge.0.) then
            ar2=pi/2.d0
           else if (xxs.lt.0.) then
            ar2=-pi/2.d0
           end if
         zartan =dcmplx(ar2,- ar1)/2.d0
         end if

* remaining cases:
         ar2=2.d0*xxs/(1.d0-xas)

         if(xas.lt.1.d0)  then     ! 1st and 4th quadrant
         zartan=dcmplx(atan(ar2),- ar1)/2.d0
         else if (xas.gt.1. .and. xxs.ge.0.) then       ! 2nd quadrant
         zartan=dcmplx(pi+atan(ar2),- ar1)/2.d0
         else if(xas.gt.1. .and. xxs.lt.0.) then        ! 3rd quadrant
         zartan=dcmplx(-pi+atan(ar2),- ar1)/2.d0
         end if
       return
       end
C (C) Copr. 1/1999  Alexander Moroz
