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
      integer LMX,LCS,ILCS,ikl,ieps,istep,ide,ndefp,itter
      integer NOUT,NOUTI,NSTEP,NFIN,NMAT,NP,NPP,NDGS,NDGSP
      real(dp) TOL,DEFP,DEFPP,DDELT,DDELTP,x_max,x_min
      real(dp) hlength_max,hlength_min,rl_min,rl_max
      real(dp) lambda_min,lambda_max,omega_max
      complex(dp) ZEPS0,CCEPS,CSEPS           !,ZARTAN
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

      real(dp) RMF(lcs),rff(lcs),RMUF,ff,filfrac   !ff for Bruggemann; filfrac for ZnS
      real(dp) xs,lambda,rap,rev_beg
      real(dp) enw,xstep,revf,revin,revinl,revl
      real(dp) omf(NFIN),omxf,reepsz,plasma,omxp
      real(dp) delo,omega0,omega,rsnm,hlength
      real(dp) RAT,RATP,AXI,REV,REVP,ALPHA,BETA,ALPHAE,BETAE      !common block variables
      real(dp) THET0,THET,THETV,PHI0,PHI     !common block variables
      real*16 ceps1real(NFIN),  ceps1imag(nfin)


      complex(dp) ceps1(NFIN),ZEPS1,ZEPS0V,Z1,Z2
      complex(dp) zeps(lcs+1)
      real*16 global_eps_r, global_eps_i
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
      rsnm = 100.0D0
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
      x_min = 0.45D0
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

      contains



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
      IMPLICIT real(dp) (A-H,O-Z)
      integer nmax, i, inm1, ixxx, ja, jam, jb, jbm, k1, k2, kk1, kk2,
     & l1, l2, lmtot, m, n1, n2, ncheck, ndgs, ngauss, nm, np
      INTEGER LMAXD,LMAX1D,LMTD
      INTEGER NAXSM,ICHOICEV,ICHOICE

      PARAMETER (LMAXD=50,LMAX1D=LMAXD+1,LMTD=LMAX1D*LMAX1D-1)

       INCLUDE 'ampld.par.f'
* number of the output unit
*
      real(dp)  LAM,LAMBDA,MRR,MRI,RSNM,HT,X(NPNG2),W(NPNG2),
     *        S(NPNG2),SS(NPNG2),AN(NPN1),R(NPNG2),DR(NPNG2),
     *        DDR(NPNG2),DRR(NPNG2),DRI(NPNG2),ANN(NPN1,NPN1)
      real(dp) TR1(NPN2,NPN2),TI1(NPN2,NPN2)
c      real(dp) XALPHA(300),XBETA(300),WALPHA(300),WBETA(300)

      complex(dp) CZERO,COPTE,COPTM
      complex(dp) zeps1,zeps0
      complex(dp) TMT(4,LMTD,LMTD)
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
* transfers real(dp) EPS(DEFP),RAT,A(REV),ALPHA,BETA,DDELT
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
      stop
      end if

      coptm=coptm+dble(TMT(2,JA,JA))

      if (abs(coptm).gt.1.4d-2) then
      write(6,*)'Optical theorem violated for MM-mode and
     & L1=',L1,'M=',N1
      write(6,*)'abs(coptm)=',abs(coptm)
      write(6,*)'(True violation only in the absorptionless case!)'
      stop
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
      IMPLICIT real(dp) (A-H,O-Z)
      INTEGER NOUT,NAXSM,ICHOICEV,ICHOICE
      LOGICAL YNCHECK
      integer nmax, np, inm1, ixxx, m, m1, n, n1, n11, n2 ,n22, ncheck,
     & ndgs, ngaus, ngauss, nm, nma, nn1, nn2, nnm, nnnggg

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
      REAL*4
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
      REAL*4
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

C********************************************************************

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
      real(dp) DLAM,DK,LAMBDA,KEX,FAC,CEXT
      real(dp) TP,TP1,PP,PP1,FC,FS,REV,PIN,PIN2,PIgrad,THETP,PHIP,
     & THETP1,
     & PHIP1,EPS,DNN,RN,DCTH0,DCTH,PH,DV1NN,DV2NN,DV1N,DV2N,D11,D12,
     & D21,D22
      real(dp) DV1(NPN6),DV2(NPN6),DV01(NPN6),DV02(NPN6)
* Either Mischenko:
      REAL*4
     &     TR11(NPN6,NPN4,NPN4),TR12(NPN6,NPN4,NPN4),
     &     TR21(NPN6,NPN4,NPN4),TR22(NPN6,NPN4,NPN4),
     &     TI11(NPN6,NPN4,NPN4),TI12(NPN6,NPN4,NPN4),
     &     TI21(NPN6,NPN4,NPN4),TI22(NPN6,NPN4,NPN4)
* or Bohren:
        complex(dp) CT(ST2,ST2,ST1)
*
      complex(dp) CN,CN1,CN2,CI,VV,VH,HV,HH,ZEPS0,CT11,CT12,CT21,CT22,
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
      PIgrad=PIN/180D0            !=PI/180
* conversion from degrees to radians:
      THETP=TP*PIgrad
      PHIP=PP*PIgrad
      THETP1=TP1*PIgrad
      PHIP1=PP1*PIgrad
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

      real(dp) xt
      integer lmx,l
      complex(dp) zx,z1,ci,cone,CQEPSt
      complex(dp) dr1(0:LMAXD),dr3(0:LMAXD),zeta(0:lmx),dzeta(0:lmx)

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
      real(dp) TOL,RINT,rsnm,lambda,thetv,thet

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

      real(dp) pi,omega,xelinv,delo,delt,rmff,rmfc
      real(dp) phi,xint,aint,xx
      real(dp) rmf(lcs),X(100),r(100)

      complex(dp) RX(2),SG(2),CQEPS(2)                  !ZSIGMA
      complex(dp)  cdl,zrfac,zsfac,zint(3),zeps(lcs+1)
      complex(dp)  ci,czero,zeinc(3,2),zescat(3,2)
      complex(dp)  zmfac,zmpifac,zpifac,zelinc(2)
      complex(dp)  cpole(lmtd,3),bpole(lmtd,3),ppole(lmtd,3)
      complex(dp)  cpolinc(lmtd,3),bpolinc(lmtd,3),ppolinc(lmtd,3)
      complex(dp)  mrpole(lmtd,3),nrpole(lmtd,3),mspole(lmtd,3),
     & nspole(lmtd,3),mrpolinc(lmtd,3),nrpolinc(lmtd,3)
      complex(dp) ZA(LMTD),ZB(LMTD)
      complex(dp) AM(lmaxd,lcs+1),AE(lmaxd,lcs+1),BM(lmaxd,lcs+1),
     & BE(lmaxd,lcs+1)
C                      =====================

      complex(dp) JL(0:lmaxd), NL(0:lmaxd)
      complex(dp) DRJL(0:lmaxd), DRNL(0:lmaxd)
      complex(dp) HL(0:lmaxd),DRHL(0:lmaxd)
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
        stop
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
      stop
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
      real(dp) DDV1(LMAXD1),DV1(LMAXD1),DV2(LMAXD1),X
      complex(dp)  cpole(lmtd,3),bpole(lmtd,3),ppole(lmtd,3)
      complex(dp)  ci,czero

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
      real(dp) A,X,QS,D1,D2,D3,DER,DN,DX,QN,QN1,QN2,
     & QNM,QNM1,QMM
      real(dp) DDV1(LMAXD1), DV1(LMAXD1), DV2(LMAXD1)

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
      real(dp) A,X,QS,QS1,D1,D2,D3,DER,DN,DSI,DX,QN,QN1,QN2,
     & QNM,QNM1,QMM
      real(dp) DDV1(LMAXD1), DV1(LMAXD1), DV2(LMAXD1)

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
      real(dp) LAMBDA,EPS,KI,KEX,THET0,THET,PHI0,PHI
      real(dp) CEXT,CSCA,QEXT,QSCA,QABS,CGN1,CGN2,DLAM
      real(dp) DELTA,SEN,COSS,XM,CGN(200)
      real(dp) AP1(SDd)
      real(dp) DPRIMA,PI,PG,H,DB(st2+1)
      real(dp) PX1(SSs,SAa),PX2(SSs,SAa),PX3(SSs,SAa),PX4(SSs,SAa)
      real(dp) AP2(SDd),AP3(SDd),AP4(SDd),BP1(SDd),BP2(SDd)
cc      real(dp)       ASIM,ALFA,CTE,CZ,DN,EM,EPS1,EPS2,EPS3,KE1,KE2,KE3,
cc    &  KI1,KI2,KI3,KM,KIM,sum,lam,NH,pc1,pc2,PC,PE,QBACK
cc      real(dp) AE1(SAa),AE2(SAa),AE3(SAa),AE4(SAa),BE1(SAa),BE2(SAa)
      real(dp) A1(SAa),A2(SAa),A3(SAa),A4(SAa),B1(SAa),B2(SAa),ANG(SAa)
      real(dp) FALLOS(SAa)
      complex(dp) D3(SUu,SUu,SDd),D4(SUu,SUu,SDd),D5(SUu,SUu,SDd)
      complex(dp) AX1(SUu,SUu,SDd),AX2(SUu,SUu,SDd)
      complex(dp) BX1(SUu,0:SDd,SDd),BX2(SUu,0:SDd,SDd)
      complex(dp) G1(SDd),G2(SDd)
      complex(dp) CT,Z1,Z2,Z3,Z4,Z5,TA(2),ZKE,MR,MR1,MR2
      complex(dp) T1(ST2,ST2,ST1),T(ST2,ST2,ST1)
      complex(dp) G3(SDd),G4(SDd),G5(SDd),ZEPS(lcs+1)

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
        stop
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
        stop
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
      real(dp) CBX,CEXT,CSCA,K1,K2,K3,K4,K5,K6,DELTA,C1(2),C2(2),EPS
      complex(dp) ZKE,MR,CCX,T(ST2,ST2,ST1),T1(ST2,ST2,ST1)
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
                stop
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
C      complex(dp) T(ST2,ST2,ST1)
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
C      complex(dp) HA(ST1),BC(ST1),BC2(ST1) ... h,j,y spherical Bessel functions
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
      DOUBLE PRECISION DEL(ST),D(ST1),drp(ST1),GM1
      DOUBLE PRECISION KAYUDA,CEXT,CSCA,DELT,EMACH
      complex(dp) ZX1,ZAX,ZDX,ZIA(ST1),ZIZ(ST1),ZIB(ST1),
     & ZJA(ST1),ZJZ(ST1),ZKA(ST1),CB(ST1)
      complex(dp) A(ST2,ST2),B(ST2,ST2),AA(ST2,ST2),BB(ST2,ST2)
      complex(dp) HA(ST1),BC(ST1),BC2(ST1),ZWIG2(ST1)
      complex(dp) ZDERI,ZDESEN,ZI1,ZA,ZB,ZC,ZD,ZE,ZF,ZFX,ZFY,ZPJ,ZR
      complex(dp) ZGM2,ZKE,ZKR,MR,ZMD,MKR,CZERO
      complex(dp) T(ST2,ST2,ST1),T1(ST2,ST2,ST1),ZX(ST2),CCX
      complex(dp) JL(0:ST),NL(0:ST),DRJL(0:ST),DRNL(0:ST),
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
            drp(I)=DFLOAT(I)*DFLOAT(I+1)
            if(I.NE.0) GM(I)=SQRT((2.0D0*D(I)+1.0D0)/drp(I))
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
            ZWIG2(I)=WIG(I+1)*drp(I)*ZDESEN            !real(dp) WIG2 is complex(dp) now
            IF(CA.EQ.2) CB(I)=MR*BC2(I)/BC2(I+1)
      END DO

      DO 20 I=NMIN,NMAX              !NMIN=MAX(1,M)
       DO 10 J=NMIN,NMAX
*
            ZPJ=D(I)*D(J)/ZKR                     !real(dp) IPJ is complex(dp) now
            ZI1=WIG(I+1)*WIG(J+1)*ZDESEN         !real(dp) I1 is complex(dp) now
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
        ZC=ZI1*(drp(I)*ZIB(J)+drp(J)*ZIA(I)-ZPJ*(D(I)+D(J)+2.0D0))
        ZR=ZC+ZKA(I)*drp(J)*ZI1                         !drp(J)=dble(J*(J+1))
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
        ZC=ZC+(CB(J)-ZIB(J))*drp(I)*ZI1
        ZR=ZC+ZKA(I)*drp(J)*ZI1
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
        IF(M.NE.0) ZE=ZE+(drp(M)-D(M))*WIG(I+1)*WIG(J+1)
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
      real(dp) cc1,cc2,dd1,dd2,cgn(200),CX(200),C1,SC,CMAX,CMAX2
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
                        stop
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
      real(dp) ff,filfrac,pi,xs,lambda,reepsz,rsnm,rmuf,omxf,omega,
     * omxp,plasma
      real(dp) omf(nfin)
      complex(dp) ci,z1,z2,zeps0,ZEPS1
      complex(dp) ceps1(NFIN)
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
      real(dp) pi,rev,rmuf,omf(NFIN)
      complex(dp) ceps1(NFIN),ZEPS1

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
      real(dp) TOL
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
      real(dp) ZEPS0
      real(dp) rmf(lcs),rmuf,xs,rsnm,pi
      real(dp) omega,lambda
*
      complex(dp) ci,czero,zeps(lcs+1),cqeps(2)
      complex(dp) RX(2),SG(2),zs,zz1,zz2
      complex(dp) KMT(2,lmaxd),TMT(4,LMTD,LMTD),alpha(lmaxd),beta(lmaxd)
*
      complex(dp) cm(lcs,lmaxd),dm(lcs,lmaxd),ce(lcs,lmaxd),
     & de(lcs,lmaxd)
      complex(dp) tt1(2,2,lmaxd,2),tt2(2,2,lmaxd,2)
      complex(dp) AM(lmaxd),AE(lmaxd),BM(lmaxd),BE(lmaxd)
*
      complex(dp) JL(0:lmaxd),NL(0:lmaxd)
      complex(dp) DRJL(0:lmaxd),DRNL(0:lmaxd)
      complex(dp) UL(2,0:lmaxd),VL(2,0:lmaxd)
      complex(dp) DRUL(2,0:lmaxd),DRVL(2,0:lmaxd)
      complex(dp) ZARTAN
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
      stop
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
      stop
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
      real(dp) lambda,xrot,xperp,rev,eps0,pi,sext(6),xk,xmn,xmj,xx,xe,
     &  xlz,xlx,xdz,xdzrn,xdx,xvol,xlp,xdp,xdprn,xcr,xfx,xapl,pf,qf
      complex(dp) ci,cone,zalph(6),zeps

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
      real(dp) TOL,OMEGA,pi,rsnm,rmuf
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
      real(dp) xs,xd1,xd2,lambda
*
      complex(dp) cp,zeps0
      complex(dp) ci,czero
*
      complex(dp) AM(LMAXD,lcs+1),AE(LMAXD,lcs+1),BM(LMAXD,lcs+1),
     & BE(LMAXD,lcs+1)

* Array declarations:

      real(dp) RMF(lcs)
      complex(dp) RX(2),SG(2),CQEPS(2),zeps(lcs+1)
*
* coated sphere declarations:
*                    moving lmax ===>
*
      complex(dp) cmx11(1:lmaxd),cmx12(1:lmaxd),cmx21(1:lmaxd),
     & cmx22(1:lmaxd)
      complex(dp) tt1(2,2,1:lmaxd,2,lcs),tt2(2,2,1:lmaxd,2,lcs)
*
C--------/---------/---------/---------/---------/---------/---------/--
C  KMT,NML,DML,NEL,DEL,EPS,RX,SG, and Bessel functions below must be
C  declared complex in the case of absorptiom and/or amplification
C--------------------------------------------------------------------
C                      ====================
C                           SUBR MODIF:
C                      =====================

      complex(dp) JL(0:lmaxd), NL(0:lmaxd)
      complex(dp) DRJL(0:lmaxd), DRNL(0:lmaxd)
      complex(dp) UL(2,0:lmaxd),DRUL(2,0:lmaxd)
      complex(dp) WL(2,0:lmaxd), DRWL(2,0:lmaxd)
      complex(dp) HL(0:lmaxd),DRHL(0:lmaxd)
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

      use libcylinder
      IMPLICIT none
!     INTEGER NPN1,NPNG1,NPNG2,NPN2,NPL,NPN3,NPN4,NPN5,NPN6
!     INCLUDE 'ampld.par.f'

      integer n,LMAX,M,I,I2
      real(dp) A,X,QS,D1,D2,D3,DER,DN,DX,QN,QN1,QN2,
     & QNM,QNM1,QMM
      real(dp) DDV1(NPN6), DV1(NPN6), DV2(NPN6)

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


      SUBROUTINE znsrefind(LAMBDA,FILFRAC,zeps)
C--------/---------/---------/---------/---------/---------/---------/--
C   FILFRAC ... filfrac of ZnS in ZnS core
C--------/---------/---------/---------/---------/---------/---------/--
      IMPLICIT NONE
      complex(dp)              zeps
      real(dp)          F,FILFRAC,LAMBDA,EPSHOST
      real(dp)          EPSZnS
      real(dp)          eMG1,EPSPAR
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


C (C) Copr. 09/1998  Alexander Moroz
C (C) Copr. 10/2005  Alexander Moroz
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
      complex(dp) zartan,zs
      real(dp) xxs,xas,ar1,ar2,pi
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
      end program
