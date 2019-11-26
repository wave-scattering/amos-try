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

      use libcylinder
      implicit none
      integer LCS,ILCS,ikl,ieps,istep,ide,ndefp,itter
      integer NOUT,NOUTI,NSTEP,NFIN,NMAT,NP,NPP,NDGS,NDGSP
      real(dp) TOL,DEFP,DEFPP,DDELT,DDELTP,x_max,x_min
      real(dp) hlength_max,hlength_min,rl_min,rl_max, nanorod_cap_hr
      complex(dp) ZEPS0,CCEPS,CSEPS           !,ZARTAN
      character(1) ync,yncv
      logical ynperfcon,ynperfconv,ynintens,ynoptth,ynbrug,yncheck
cc      external ZARTAN

c Parameters:
C ::: number of the output unit for cross sections and scattering matrix
      PARAMETER (NOUT=35)
C ::: number of the output unit for the field intensity
      PARAMETER (NOUTI=60)
c Maximal number of spherical harmonics used. The floating number is
c specified below by the value of LMAX parameter
!     PARAMETER (lmx=100)
*
* If convergence test in the calculation of the scattering cross sections
* is to be performed, yncheck=.true., otherwise yncheck=.false.
!     parameter (yncheck=.false.)
      parameter (yncheck=.true.)
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

c Temporarily, before the coated part is finished,
c read in particle core dielectric function to CSEPS
c particle (core) dielectric constant (depending whether lcs=1 or lcs=2)
c
! Should be set in *.ini config file
C      PARAMETER (CCEPS=(80.0d0,0.0000000d0))
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

      real(dp) RMF(lcs),rff(lcs),RMUF,ff,filfrac   !ff for Bruggemann; filfrac for ZnS
      real(dp) lambda,rev_beg
      real(dp) enw,xstep,revf,revin,revinl,revl
      real(dp) omf(NFIN),omxf,reepsz,plasma,omxp
      real(dp) delo,omega0,omega,rsnm,hlength
      real(dp) RAT,RATP,AXI,REV,REVP,ALPHA,BETA,ALPHAE,BETAE      !common block variables
      real(dp) THET0,THET,THETV,PHI0,PHI     !common block variables
      real(dp) ceps1real(NFIN),  ceps1imag(nfin)


      complex(dp) ceps1(NFIN),ZEPS1,ZEPS0V,Z1,Z2
      complex(dp) zeps(lcs+1)
      real(dp) global_eps_r, global_eps_i
*
      COMMON /TOAMPLD/RAT,REV,ALPHA,BETA,DDELT
*
* transfers real(dp) RAT,REV,ALPHA,BETA,DDELT  from the main to AMPLDR
*
      COMMON /TOTAMPLD/THET0,THET,PHI0,PHI
*
* transfers real(dp) THET0,THET,PHI0,PHI from the main to AMPLDR
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
* transfers real(dp) DEFP,RAT,REV,ALPHAE,BETAE,DDELT from the main to TMTAXSPV
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
C--------/---------/---------/---------/---------/---------/---------/--
      !---------------------------------------------------------------
      !     Read config
      !---------------------------------------------------------------

      call cli_parse
      call ini_parse

c background dielectric constant
!     PARAMETER (ZEPS0=1.D0) !set in ini_parse
c

*
* From here to spherec
*---------------------------------------------------------------
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
!     rsnm=1.d0 ! Should be set in *.ini file
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
      write(6,*)'nanorod/finite cylinder capped by half-spheroids: -9'
      write(6,*)'intensity around homogeneous/coated sphere: type -50'

      Open(unit = 90,file = 'epsWater.txt', status = 'unknown')
*
      if (NP == 0) then
      write(6,*) 'Input the particle type:'
      read(5,*) NP
      else
      write(6,*)'Auto-select from *.ini file: ', NP
      end if
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
!     -9                  nanorod (cylinder capped with half-spheroids)
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
C     NP=-9 - DEFP = the ratio of the cylinder diameter to its length.
C             Spheroid caps of the nanorod are defined by additional
C             parameter.
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
      else if ((NP .eq. -2) .or. (NP .eq. -9)) then
*
      if (rl_max < 0_dp) then
        write(6,*)'Enter cylinder maximal r/l:'
        read(5,*) rl_max
      else
        write(6,*)'Auto-set cylinder maximal r/l from *.ini', rl_max
      end if

      if (rl_min < 0_dp) then
        write(6,*)'Enter cylinder minimal r/l:'
        read(5,*) rl_min
      else
        write(6,*)'Auto-set cylinder minimal r/l from *.ini', rl_min
      end if

      if (ndefp <= 0) then
        write(6,*)'Enter amount of steps in length:'
        read(5,*) ndefp
      else
        write(6,*)'Auto-set amount of steps in length from *.ini',ndefp
      end if

      if (rsnm <= 0_dp) then
        write(6,*)'Enter cylinder radius:'
        read(5,*) rsnm
      else
        write(6,*)'Auto-set cylinder radius from *.ini',rsnm
      end if

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
C      stop
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
      stop
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
      !ellipsoid(sphere), cylinder, and nanorod
      IF (NP.EQ.-1.OR.NP.EQ.-2.OR.NP.EQ.-9) NCHECK=1
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
      write(6,*)'Auto-set from *.ini:', thet0, phi0
*
      write(6,*)'Specify (theta,phi) angles of the scattered beam
     1             (in degrees)'
!     read(5,*) THET,PHI
      write(6,*)'Auto-set from *.ini', thet, phi
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

      if (x_min < 0) then
        write(6,*)'Enter minimum value of x-parameter:'
        read(5,*) x_min
      else
        write(6,*)'Auto-set x_min from *.ini ', x_min
      end if

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
      if (x_max < 0) then
        write(6,*)'Enter maximum value of x-parameter:'
        read(5,*) x_max
      else
        write(6,*)'Auto-set x_max from *.ini ', x_max
      end if

      enw = 2*pi*rsnm/x_max
c      enw=500
      if (nstep <= 0) then
        write(6,*)'Enter amount of scanning steps:'
        read(5,*) nstep
      else
        write(6,*)'Auto-set nstep from *.ini ',nstep
      end if
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
       if (aIMAG(z2).GE.0.0) then
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

! TODO: commented out lisac just to compile, call has not enough arguments.
!subroutine lisac(de,lcs,nmx1,nmx2,ngauss,lambda,eps,ki,kex,
!     & THET0,THET,PHI0,PHI,zeps)
!     call lisac(ide,lcs,lmax,lmax,ndgs*lmax,lambda,defpp,
!    & revinl/revl,revl,zeps)

      else
      global_eps_r =REAL( zeps(1))
      global_eps_i = AIMAG(zeps(1))


C      write(90,*) lambda,  global_eps_r,
C     & global_eps_i
      call ampldr(yncheck,lmax,ichoice,npp,defpp,nanorod_cap_hr,
     & rsnm,hlength,lambda,zeps(1),zeps0)
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

      contains

C--------/---------/---------/---------/---------/---------/---------/--
      subroutine ini_parse()
      character(len=:), allocatable :: items(:,:) !< Items pairs.
      integer                       :: error
      character(len=:), allocatable :: string         !< String option.
      real(dp)                      :: double, d_re, d_im
      integer                       :: num
      write(6,*) 'Reading config from file: ', ini_config_file_name
      call fini%load(filename=ini_config_file_name)

      allocate(character(999) :: string)
      string = repeat(' ', 999)

      call fini%get(section_name='general', option_name='particle_type',
     &  val=string, error=error)
      NP = 0
      if ((trim(string) .eq. 'cylinder').or.(trim(string)=='-2')) NP=-2
      if ((trim(string) .eq. 'nanorod').or.(trim(string)=='-9')) NP=-9

      call fini%get(section_name='nanorod',
     &  option_name='nanorod_cap_hr', val=double, error=error)
      nanorod_cap_hr = 1_dp ! default is round cap
      if (error==0) nanorod_cap_hr = double

      call fini%get(section_name='cylinder', option_name='rl_min',
     &  val=double, error=error)
      rl_min = -1_dp
      if (error==0) rl_min = double

      call fini%get(section_name='cylinder', option_name='rl_max',
     &  val=double, error=error)
      rl_max = -1_dp
      if (error==0) rl_max = double

      call fini%get(section_name='cylinder', option_name='rl_steps',
     &  val=num, error=error)
      ndefp = 0
      if (error==0) ndefp = num

      call fini%get(section_name='cylinder', option_name='radius',
     &  val=double, error=error)
      rsnm = -1_dp
      if (error==0) rsnm = double

      call fini%get(section_name='cylinder', option_name='alpha',
     &  val=double, error=error)
      alpha = 0_dp
      if (error==0) alpha = double

      call fini%get(section_name='cylinder', option_name='beta',
     &  val=double, error=error)
      beta = 0_dp
      if (error==0) beta = double

      call fini%get(section_name='general',
     &  option_name='background_epsilon',
     &  val=double, error=error)
      zeps0 = 1_dp
      if (error==0) zeps0 = double

      d_re = 1_dp
      call fini%get(section_name='cylinder', option_name='eps_real',
     &  val=d_re, error=error)
      d_im = 0_dp
      call fini%get(section_name='cylinder', option_name='eps_imag',
     &  val=d_im, error=error)
      cceps = d_re + ci*d_im

      call fini%get(section_name='beam', option_name='theta0',
     &  val=double, error=error)
      thet0 = 0_dp
      if (error==0) thet0 = double

      call fini%get(section_name='beam', option_name='phi0',
     &  val=double, error=error)
      phi0 = 0_dp
      if (error==0) phi0 = double

      call fini%get(section_name='beam', option_name='theta',
     &  val=double, error=error)
      thet = 0_dp
      if (error==0) thet = double

      call fini%get(section_name='beam', option_name='phi',
     &  val=double, error=error)
      phi = 0_dp
      if (error==0) phi = double

      call fini%get(section_name='beam', option_name='x_min',
     &  val=double, error=error)
      x_min = -1_dp
      if (error==0) x_min = double

      call fini%get(section_name='beam', option_name='x_max',
     &  val=double, error=error)
      x_max = -1_dp
      if (error==0) x_max = double

      call fini%get(section_name='beam', option_name='x_steps',
     &  val=num, error=error)
      nstep = 0
      if (error==0) nstep = num

      end subroutine ini_parse



C**********************************************************************

      SUBROUTINE AMPLDR(yncheck,nmax,ichoicev,np,eps,nanorod_cap_hr,
     &                  rsnm,ht,lambda,zeps1,zeps0)

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
!     NP=-9 - EPS = the ratio of the cylinder diameter to its length,
!             nanorod caps are set independantly.
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
      IMPLICIT real(dp) (A-H,O-Z)
      INTEGER NOUT,NAXSM,ICHOICEV,ICHOICE
      LOGICAL YNCHECK
      integer nmax, np, inm1, ixxx, m, m1, n, n1, n11, n2 ,n22, ncheck,
     & ndgs, ngaus, ngauss, nm, nma, nn1, nn2, nnm, nnnggg
      real(dp) nanorod_cap_hr

       INCLUDE 'ampld.par.f'
* number of the output unit
      PARAMETER (NOUT=35)
*
      real(dp)  LAM,LAMBDA,MRR,MRI,RSNM,HT,DDELT,DDELTA,
     *        X(NPNG2),W(NPNG2),
     *        S(NPNG2),SS(NPNG2),AN(NPN1),R(NPNG2),DR(NPNG2),
     *        DDR(NPNG2),DRR(NPNG2),DRI(NPNG2),ANN(NPN1,NPN1)
      real(dp) TR1(NPN2,NPN2),TI1(NPN2,NPN2)
c      real(dp) XALPHA(300),XBETA(300),WALPHA(300),WBETA(300)
      real(dp)
     &     RT11(NPN6,NPN4,NPN4),RT12(NPN6,NPN4,NPN4),
     &     RT21(NPN6,NPN4,NPN4),RT22(NPN6,NPN4,NPN4),
     &     IT11(NPN6,NPN4,NPN4),IT12(NPN6,NPN4,NPN4),
     &     IT21(NPN6,NPN4,NPN4),IT22(NPN6,NPN4,NPN4)
      complex(dp) S11,S12,S21,S22
      complex(dp) zeps1,zeps0
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
* transfers real(dp) RAT,A(REV),ALPHA,BETA,DDELT from the main here
*
      COMMON /TOTAMPLD/THET0,THET,PHI0,PHI
*
* transfers real(dp) THET0,THET,PHI0,PHI from the main here

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
      MRI=aimag(SQRT(ZEPS1/ZEPS0))
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
      IF (DABS(RAT-1D0).GT.1D-8.AND.NP.EQ.-9)
     &  CALL SAREAnanorod (EPS,RAT,nanorod_cap_hr)
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
         CALL VARY(LAM,MRR,MRI,A,EPS,nanorod_cap_hr,
     &              RSNM,HT,NP,NGAUSS,X,P,
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
            DN1=dble(2*N+1)
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
         IF (NGAUS.gt.300.and.NGAUS.lt.2595) cycle
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
         CALL VARY(LAM,MRR,MRI,A,EPS,nanorod_cap_hr,
     &              RSNM,HT,NP,NGAUSS,X,P,
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

            DN1=dble(2*N+1)

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
         CALL VARY(LAM,MRR,MRI,A,EPS,nanorod_cap_hr,
     &             RSNM,HT,NP,NGAUSS,X,P,
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

      DO N2=1,NMAX
         NN2=N2+NMAX
         DO N1=1,NMAX
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
         end do
      end do !end of the loop over orbital numbers
*________________


*<<<

      if (abs(qsca).gt.(1.0001d0*abs(qext))) then
         write(6,*)'M=',0
         write(6,*)'QSCA=',QSCA
         write(6,*)'QEXT=',QEXT
       write(6,*)
     & 'WARNING: abs(qsca).gt.abs(qext)!!!'
c     stop
      end if

      if (qext.gt.1.d-7) then
       write(6,*)
     & 'WARNING: Partial sum QEXT has to be negative!'
c      stop
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
         DO N2=1,NM              !summation over orbital numbers

* conversion of the N22 index of RT1 and IT1 matrices
* to the index NN2 of RT^{ij} and IT^{ij} matrices

            NN2=N2+M-1        !from M to NMAX
            N22=N2+NM         !from NMAX+1 to 2*NMAX-M+1

            DO N1=1,NM           !summation over orbital numbers

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
            end do
         end do     !end of the loop over orbital numbers

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
c     stop
      end if

      if (qxt.ge.1d-7) then
       write(6,*)
     & 'WARNING: Partial sum QXT has to be negative!'
C      stop
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
      Z11=0.5D0*(S11*conjg(S11)+S12*conjg(S12)
     &          +S21*conjg(S21)+S22*conjg(S22))
      Z12=0.5D0*(S11*conjg(S11)-S12*conjg(S12)
     &          +S21*conjg(S21)-S22*conjg(S22))
      Z13=-S11*conjg(S12)-S22*conjg(S21)
      Z14=(0D0,1D0)*(S11*conjg(S12)-S22*conjg(S21))
      Z21=0.5D0*(S11*conjg(S11)+S12*conjg(S12)
     &          -S21*conjg(S21)-S22*conjg(S22))
      Z22=0.5D0*(S11*conjg(S11)-S12*conjg(S12)
     &          -S21*conjg(S21)+S22*conjg(S22))
      Z23=-S11*conjg(S12)+S22*conjg(S21)
      Z24=(0D0,1D0)*(S11*conjg(S12)+S22*conjg(S21))
      Z31=-S11*conjg(S21)-S22*conjg(S12)
      Z32=-S11*conjg(S21)+S22*conjg(S12)
      Z33=S11*conjg(S22)+S12*conjg(S21)
      Z34=(0D0,-1D0)*(S11*conjg(S22)+S21*conjg(S12))
      Z41=(0D0,1D0)*(S21*conjg(S11)+S22*conjg(S12))
      Z42=(0D0,1D0)*(S21*conjg(S11)-S22*conjg(S12))
      Z43=(0D0,-1D0)*(S22*conjg(S11)-S12*conjg(S21))
      Z44=S22*conjg(S11)-S12*conjg(S21)


      WRITE(NOUT+10,5001) Z11,Z12,Z13,Z14
      WRITE(NOUT+10,5001) Z21,Z22,Z23,Z24
      WRITE(NOUT+10,5001) Z31,Z32,Z33,Z34
      WRITE(NOUT+10,5001) Z41,Z42,Z43,Z44

 5001 FORMAT (4F10.4)

c      ITIME=MCLOCK()
c      TIME=dble(ITIME)/6000D0
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
      IMPLICIT real(dp) (A-B,D-H,O-Z)
      IMPLICIT complex(dp) (C)
      INTEGER NOUT
      integer nmax, i, j, k, m, m1, n, nmin, nn

* number of the output unit
      PARAMETER (NOUT=35)
      INCLUDE 'ampld.par.f'
      real(dp) DLAM,LAMBDA,CEXT,CEXT1,CEXT2
      real(dp) AL(3,2),AL1(3,2),AP(2,3),AP1(2,3),B(3,3),
     *       R(2,2),R1(2,2),C(3,2),CA,CB,CT,CP,CTP,CPP,CT1,CP1,
     *       CTP1,CPP1
      real(dp) DV1(NPN6),DV2(NPN6),DV01(NPN6),DV02(NPN6)
      REAL(dp)
     &     TR11(NPN6,NPN4,NPN4),TR12(NPN6,NPN4,NPN4),
     &     TR21(NPN6,NPN4,NPN4),TR22(NPN6,NPN4,NPN4),
     &     TI11(NPN6,NPN4,NPN4),TI12(NPN6,NPN4,NPN4),
     &     TI21(NPN6,NPN4,NPN4),TI22(NPN6,NPN4,NPN4)
      complex(dp) CAL(NPN4,NPN4),VV,VH,HV,HH,ZEPS0
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
      PIgrad=PIN/180D0            !=PI/180

* conversion from degrees to radians:
      ALPH=ALPHA*PIgrad
      BET=BETA*PIgrad
      THETL=TL*PIgrad
      PHIL=PL*PIgrad
      THETL1=TL1*PIgrad
      PHIL1=PL1*PIgrad

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
!     CI=(0D0,1D0)

C >>> ALPHA numerical prefactors without phi-angles
C     (following Eq. (28))

      DO NN=1,NMAX
         DO N=1,NMAX
            CN=CI**(NN-N-1)
            DNN=dble((2*N+1)*(2*NN+1))
            DNN=DNN/dble( N*NN*(N+1)*(NN+1) )
            RN=DSQRT(DNN)
            CAL(N,NN)=CN*RN
         end do
      end do

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

      DO M=0,NMAX
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
         DO NN=NMIN,NMAX

            DV1NN=DV01(NN)           !\pi-functions
            DV2NN=DV02(NN)           !\tau-functions

            DO  N=NMIN,NMAX
               DV1N=DV1(N)           !\pi-functions
               DV2N=DV2(N)           !\tau-functions

               CT11=cmplx_dp(TR11(M1,N,NN),TI11(M1,N,NN))
               CT22=cmplx_dp(TR22(M1,N,NN),TI22(M1,N,NN))

               IF (M.EQ.0) THEN     !T^{21}=T^{12}=0 in particle frame

                  CN=CAL(N,NN)*DV2N*DV2NN

                  VV=VV+CN*CT22
                  HH=HH+CN*CT11

                 ELSE   !T^{21}\neq T^{12}\neq 0

                  CT12=cmplx_dp(TR12(M1,N,NN),TI12(M1,N,NN))
                  CT21=cmplx_dp(TR21(M1,N,NN),TI21(M1,N,NN))

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
            end do
         end do     !(over n,n')
      end do        !end of main summation loop (over m)

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
      cext1=2.d0*dlam*aimag(VV)      !Eq. (2.159)

C Extiction for E perpendicular to the axis of axial symmetry:
      cext2=2.d0*dlam*aimag(HH)      !Eq. (2.159)

C Orientationally averaged extiction
      cext=dlam*aimag(VV+HH)         !Eq. (5.97)

      write(6,*)'C_{ext}=\fr{2\pi}{k_1} Im (S_{11}+S_{22})=',
     & cext               !=2.d0*PIN*aimag(VV+HH)/k_1

      write(10, *) lambda, cext1, cext2
cc      FAC=lambda**2/(2.d0*PIN**2*REV**2)     !=2/xs**2
      FAC=1.d0/(PIN*REV**2)   !an effective geom. cross section

      write(nout+3,1107)lambda,fac*cext,fac*cext1,fac*cext2
      write(nout+12,1105) lambda, VV,VH
      write(nout+12,1105) lambda, HV,HH
      write(nout+13,1106) lambda, dble(VV*conjg(vv) + VH*conjg(vh)),
     &                               dble((vv+vh)*conjg(vv+vh))
      write(nout+14,1106) lambda, dble(hV*conjg(hv) + hH*conjg(hh)),
     &                               dble((hv+hh)*conjg(hv+hh))


 1101 FORMAT ('S11=',D11.5,' + i*',D11.5)
 1102 FORMAT ('S12=',D11.5,' + i*',D11.5)
 1103 FORMAT ('S21=',D11.5,' + i*',D11.5)
 1104 FORMAT ('S22=',D11.5,' + i*',D11.5)
 1105 FORMAT (F8.2,5X,D11.5,2X,D11.5,5X,D11.5,2X,D11.5)
 1106 FORMAT (F8.2,5X,D11.5,5X,D11.5)
 1107 FORMAT (F10.4,3(5X,D16.8))

      RETURN
      END

C********************************************************************
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
       INCLUDE 'ampld.par.f'
      IMPLICIT real(dp) (A-H,O-Z)
      integer m, ngauss, nmax, ncheck, naxsm, i, i1, i2, k1, k2, kk1,
     & kk2, mm1, n, n1, n2, ng, ngss, nm, nnmax
      real(dp)  X(NPNG2),W(NPNG2),AN(NPN1),
     *        R(NPNG2),DR(NPNG2),SIG(NPN2),
     *        J(NPNG2,NPN1),Y(NPNG2,NPN1),
     *        JR(NPNG2,NPN1),JI(NPNG2,NPN1),DJ(NPNG2,NPN1),
     *        DY(NPNG2,NPN1),DJR(NPNG2,NPN1),
     *        DJI(NPNG2,NPN1),DDR(NPNG2),DRR(NPNG2),
     *        D1(NPNG2,NPN1),D2(NPNG2,NPN1),D3(NPNG2,NPN1),
     *        DRI(NPNG2),RR(NPNG2),
     *        DV1(NPN3),DDV1(NPN3),DV2(NPN3),
     *        DD1,DD2,DD3

      real(dp)  R11(NPN1,NPN1),R12(NPN1,NPN1),
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
cc      real(dp) TR1(NPN2,NPN2),TI1(NPN2,NPN2)
*________
      COMMON /TMAT99/
     &            R11,R12,R21,R22,I11,I12,I21,I22,RG11,RG12,RG21,RG22,
     &            IG11,IG12,IG21,IG22          !only between TMATR routines
      COMMON /CBESS/ J,Y,JR,JI,DJ,DY,DJR,DJI
cc      COMMON /CT/ TR1,TI1                      !output from TT routine
      COMMON /CTT/ QR,QI,RGQR,RGQI             !input for TT routine
*________
      MM1=M
      QM=dble(M)
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
      DO N=1,NM                 !NM=2*NMAX
           SI=-SI
           SIG(N)=SI              !=(-1)**N
      end do
*
* Assigning Wigner d-matrices:

      DO I=1,NGAUSS

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
      end do
*
*  Assigning r^2(\theta)*weight product:

      DO I=1,NGSS
           WR=W(I)*R(I)

cc          if (dr(i).eq.0.d0) WR=0.d0   !temporarily only

           RR(I)=WR            !W(I)*r^2(\theta)
      end do
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
      CALL TT(NM)
*
      RETURN
      END


C (C) Copr. 09/1998  Alexander Moroz
C (C) Copr. 10/2005  Alexander Moroz
      end program
