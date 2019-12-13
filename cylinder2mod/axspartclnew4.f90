program axspartcl1
    ! warning in module axspartcl in file axspartcl.f: variables set but never used:
    !   rmf set at line 257 file axspartcl.f
    !    xs set at line 682 file axspartcl.f
    !
    !  if changing lmx, one has to adjust npn1 in ampldr to the same value
    !                   npng1 is only within the ampldr
    ! extinction for a single homogeneous sphere of radius 300nm  (300.001/300)
    ! host dielectric constant= (1.00000000000000,0.000000000000000e+000)
    ! sphere diel. constant= (2.10250000000000,0.000000000000000e+000)
    ! for lambda_0=1000nm is sg_tot=1.21460719430552
    ! qsca=   2.15779683977996
    ! qext=  -2.15779683977996
    !         fac=lam**2/(2.d0*p**2*rev**2)     !=2/xs**2
    ! relates effciences to normalized cross sections. for the setting above:
    !         fac=0,56289546467965428579933035116515
    ! s11=0.32957d+03 + i*0.17171d+03
    ! s12=0.00000d+00 + i*0.00000d+00
    ! s21=0.00000d+00 + i*0.00000d+00
    ! s22=0.32957d+03 + i*0.17171d+03
    ! for unpolarized light, extinction cross section is
    ! [(2.140) and (2.159) of {mtl}]:
    !
    !         c_{ext}= \fr{2\pi}{k_1} \mb{im } (s_{11}+_{22})
    !
    ! homogeneous particle
    ! particle parameters:
    ! radius =   300.000333332963
    ! particle diel. constant= (2.10250000000000,0.000000000000000e+000)
    ! background dielectric constant zeps0= (1.00000000000000,0.000000000000000e+000)
    !
    !==========================
    !      room for improvement:
    !            1) it would be more resonable to replace rev by the
    !               size parameter k*rev
    !            2)
    !            3)
    !            4)
    !--------/---------/---------/---------/---------/---------/---------/--
    !  this routines calculates the single particle scattering properties
    !  (including coated particles)
    !
    !                    make -f mkaxspsc
    !
    ! k_l length in units (2*pi/a=) pi:    xkl= 0.8660254037844386d0
    !
    ! outputs the total elastic scattering cross section tcs
    !
    ! parameters:
    !
    ! partial wave expansion is used, which is badly convergent
    ! for large size parameters $x> 100$. in numerical applications, the
    ! series is to be cut off after
    !          lmax (lmx parameter here) \approx x+4x^{1/3}+2$.
    ! in the case of lmx=50 this means that x <= 35.
    ! if one wants to observe ripples, the cutoff for a given x has to
    ! be even larger
    !
    !  alpha and beta - euler angles (in degrees) specifying the orientation
    !            of the scattering particle relative to the laboratory reference
    !            frame (refs. 6 and 7).
    !  thet0 - zenith angle of the incident beam in degrees
    !  thet - zenith angle of the scattered beam in degrees    !with respect to
    !  phi0 - azimuth angle of the incident beam in degrees    !the laboratory frame!!!
    !  phi - azimuth angle of the scattered beam in degrees
    !
    !  ichoice=1 if nag library is available, otherwise ichoice=2
    !
    !  ncheck  -  .eq.0  then  ngss=2*ngauss, factor=1d0
    !             .eq.1  then  ngss = ngauss, factor=2d0: theta=pi/2 is mirror
    !                          symmetry plane as in the case of chebysh. particle,
    !                          ellipsoid, and cylinder
    !
    !  naxsm   -  .eq.0 : gauss abscissas do not have +/- theta symmetry
    !             .eq.1 : gauss abscissas have +/- theta symmetry
    !
    !  ndgs - controlling the number nd=ndgs*nmax of division points in
    !  computing integrals over the particle surface (ref. 5).
    !  for compact particles, the
    !  recommended value is 2. for highly aspherical particles larger
    !  values (3, 4,...) may be necessary to obtain convergence.
    !  the code does not check convergence over this parameter.
    !  therefore, control comparisons of results obtained with
    !  different ndgs-values are recommended.
    !
    !
    !  computation can be speed up if one sets mpar%yncheck=0
    !  however, then the check of gauss integrations
    !  convergence is not performed.
    !---------------------------------------------------------------------

    use libcylinder
    implicit none
    integer lcs, ilcs, ikl, ieps, istep, ide, ndefp, itter
    integer nout, nouti, nstep, nfin, nmat, np, npp, ndgs, ndgsp
    real(dp) tol, defp, defpp, ddelt, ddeltp, x_max, x_min
    real(dp) hlength_max, hlength_min, rl_min, rl_max
    complex(dp) zeps0, cceps, cseps           !,zartan
    character(1) ync, yncv
    logical ynperfcon, ynperfconv, ynintens, ynoptth, ynbrug
    !c      external zartan

    ! parameters:
    ! ::: number of the output unit for cross sections and scattering matrix
    parameter (nout = 35)
    ! ::: number of the output unit for the field intensity
    parameter (nouti = 60)
    ! maximal number of spherical harmonics used. the floating number is
    ! specified below by the value of lmax parameter
    !     parameter (lmx=100)
    !
    ! if convergence test in the calculation of the scattering cross sections
    ! is to be performed, yncheck=.true., otherwise yncheck=.false. ! Moved to *.ini config
!   parameter (yncheck=.false.)
!   parameter (yncheck = .true.)

    !
    ! if particle is coated, ync=y, otherwise ync=n
    parameter (ync = 'n')
    !
    ! ynperfcon=.true. if core is a perfect conductor, otherwise
    ! ynperfcon=.false.
    !
    parameter (ynperfcon = .false.)
    !
    ! ynintens=.true. if the field intensity is to be calculated; otherwise
    ! ynintens=.false.
    !
    parameter (ynintens = .false.)
    !
    ! number of coatings
    parameter (lcs = 1)
    ! the coating layer to which material data are read in
    parameter (ilcs = 1)
    ! if coated, the ratio 'core radius/particle radius'
    !   (if lcs.ne.1, program is singular for rff=0. and 1.) - use homogeneous
    !    particle instead!!!
    !       parameter (rff=0.95d0)

    ! temporarily, before the coated part is finished,
    ! read in particle core dielectric function to cseps
    ! particle (core) dielectric constant (depending whether lcs=1 or lcs=2)
    !
    ! should be set in *.ini config file
    !      parameter (cceps=(80.0d0,0.0000000d0))
    !      parameter (cceps=(1.45d0,0.0d0)**2)    !sio2
    !      parameter (cceps=(1.2d0,0.01d0)**2)    !to test lisac
    !      parameter (cceps=(-70.720839d0,7.05596d0))   !-70.720839,7.05596   au for 1319nm
    !      parameter (cceps=(-10.84d0,0.762d0))   !jap89_5774   ellipsoid for ld=633
    !      parameter (cceps=(-2.03d0,0.602d0))    !jap89_5774   sphere for ld=354
    ! >>>     sphere (outer shell scatterer) permittivity                  <<<
    !  n(silica)=1.45  <--->    eps(1)=2.1025d0
    !  n(zns)=2.       <--->    eps(1)=4.d0
    parameter (cseps = (1.005d0, 0.d0)**2)
    !      parameter (cseps=(1.05d0,0.d0)**2)    !to test lisac
    !      parameter (cseps=(-10.84d0,0.762d0))    !jap89_5774
    ! material code number
    !   nmat=0             dispersionless dielectric
    !   nmat=1             drude metal
    !   nmat=2             ag
    !   nmat=3             au
    !   nmat=4             zns
    !   nmat=5             cu
    !   nmat=6             al
    !   nmat=7             pt
    !   nmat=8             si
    !  nmat = 9           water
    !
    parameter(nmat = 0)
    !
    ! temporarily option for reading of the real data for the dielectric constant
    ! the number of the entries in a material data file to be read below
    ! data files for au,cu,al,pt should be ordered with the decreased wavelength
    ! (omega increases in the loop and is oriented along the data file)
    !
    !          agc.dat                nfin=73       ! from palik
    !          audat.dat              nfin=66       ! from palik
    !          au_2dat.dat            nfin=76       ! from jaw
    !          au*new.dat             nfin=142
    !          cudat.dat              nfin=47       ! from palik
    !          aldat.dat              nfin=80       ! from palik
    !          ptdat.dat              nfin=53       ! from palik
    !          nidat.dat              nfin=68       ! from palik
    !          sidat.dat              nfin=291
    !          sieps.dat              nfin=117
    !          measured_water_dispersion_t=24.txt       nfin = 3600
    !
    parameter (nfin = 3600)
    !
    ! ::: relative error allowed for the tcs. if the convergence
    !     within tol is not reached, program issues warning
    parameter (tol = 1.d-3)
    !
    ! if ynbrug=.true., performs bruggeman approximation for zeps1. otherwise
    ! ynbrug=false.
    parameter (ynbrug = .false.)
    !*****************************************************************
    ! declarations:

    integer lmax, ncheck, naxsm, ncheckp, naxsmp  !common block variables

    real(dp) rmf(lcs), rff(lcs), rmuf, ff, filfrac   !ff for bruggemann; filfrac for zns
    real(dp) lambda, rev_beg
    real(dp) enw, xstep, revf, revin, revinl, revl
    real(dp) omf(nfin), omxf, reepsz, plasma, omxp
    real(dp) delo, omega0, omega, rsnm, hlength
    real(dp) rat, ratp, axi, rev, revp, alpha, beta, alphae, betae      !common block variables
    real(dp) thet0, thet, thetv, phi0, phi     !common block variables
    real(dp) ceps1real(nfin), ceps1imag(nfin)

    complex(dp) ceps1(nfin), zeps1, zeps0v, z1, z2
    complex(dp) zeps(lcs + 1)
    real(dp) global_eps_r, global_eps_i
    !
    common /toampld/rat, rev, alpha, beta, ddelt
    !
    ! transfers real(dp) rat,rev,alpha,beta,ddelt  from the main to ampldr
    !
    common /totampld/thet0, thet, phi0, phi
    !
    ! transfers real(dp) thet0,thet,phi0,phi from the main to ampldr
    !
    common /toiampld/ncheck, naxsm, ndgs
    !
    ! transfers integers ncheck,naxsm,ndgs from the main
    ! to ampldr
    !
    common /toitmt/npp, ncheckp, naxsmp, ndgsp
    ! transfers integers np,ncheck,naxsm,ndgs from the main to tmtaxspv
    !
    common /totmt/defpp, ratp, revp, alphae, betae, ddeltp
    !
    ! transfers real(dp) defp,rat,rev,alphae,betae,ddelt from the main to tmtaxspv
    !
    common /toltmt/ ynoptth
    !
    ! transfers logical ynoptth from the main to tmtaxsp
    !
    common /dielf/ zeps0v
    common /revf/ revf
    !
    ! transfers zeps0,rev from the main to ampl
    !
    common /tospherecr/ rmuf
    common /tospherech/ yncv
    common /tospherecl/ ynperfconv
    common /cylpar/ rsnm, hlength
    common /epsgg/ global_eps_r, global_eps_i
    !--------/---------/---------/---------/---------/---------/---------/--
    !     !---------------------------------------------------------------
    !     !     read config
    !     !---------------------------------------------------------------

    call cli_parse
    call ini_parse

    cceps = mpar%cceps
    zeps0 = mpar%zeps0

    ! background dielectric constant
    !     parameter (zeps0=1.d0) !set in ini_parse
    !
    !
    ! from here to spherec
    !---------------------------------------------------------------
    ! convergence variable (has to be at least equal to 2):
    !
    lmax = 25  !lmx
    !
    !  fcc parameters:

    ! f=0.05  (--)
    !       rmuf=0.2879411911484862d0
    ! f=0.06
    !       rmuf=0.3059831741945870d0
    ! f=0.07
    !       rmuf=0.3221166265075572d0
    ! f=0.08
    !       rmuf=0.3367780601921260d0
    ! f=0.09
    !       rmuf=0.3502632974822208d0
    ! f=0.1  (0f)
    !       rmuf=0.3627831678597810d0
    ! f=0.12
    !       rmuf=0.3855146420814098d0
    ! f=0.15  (--)
    !       rmuf=0.4152830592077074d0
    ! f=0.2  (f)
    !       rmuf=0.4570781497340833d0
    ! f=0.25  (--)
    !       rmuf=0.4923725109213483d0
    ! f=0.3  (--)
    !       rmuf=0.5232238679605294d0
    ! f=0.35  (--)
    !       rmuf=0.5508116833525640d0
    ! f=0.36  (vatl)
    !       rmuf=0.5560083268891277d0
    ! f=0.4  (ff)
    !       rmuf=0.5758823822969722d0
    ! f=0.45  (--)
    !       rmuf=0.5989418136982620d0
    ! f=0.5  (--)
    !       rmuf=0.6203504908994001d0
    ! f=0.55  (--)
    !       rmuf=0.6403754763690468d0
    ! f=0.58  (2w)
    !       rmuf=0.6518131631368212d0
    ! f=0.6  (3f)
    !       rmuf=0.6592207650508868d0
    ! f=0.62  (3w)
    !       rmuf=0.6664655293794289d0
    ! f=0.64 (4f)
    !       rmuf=0.6735561203842518d0
    ! f=0.65 (--)
    !       rmuf=0.6770461107710686d0
    ! f=0.66 (4w)
    !        rmuf=0.6805004874579641d0
    ! f=0.68 (5f)
    !       rmuf=0.6873059437985093d0
    ! f=0.7  (6f)
    !       rmuf=0.6939792343839249d0
    ! f=0.72 (7f)
    !       rmuf=0.7005265949644416d0
    ! close packed :
    !       rmuf=1.d0/dsqrt(2.d0)
    rmuf = 1.d0
    rmf(lcs) = rmuf

    !ccccccccccccccccc    assignement of common variables ccccccccccc
    !
    ! ynoptth=.true. if you want to check optical theorem
    ynoptth = .false.
    !
    ynperfconv = ynperfcon
    yncv = ync
    zeps0v = zeps0
    ! preinitialization necessary since below rsnm,rev, and defp (all
    ! feeded in evanesc) are not always specified:
    !
    !     rsnm=1.d0 ! should be set in *.ini file
    rev = 1.d0
    defp = 1.d0
    !
    write(6, *)'chose particle shape'
    write(6, *)'only axially symmetric particle shapes allowed'
    write(6, *)'(the axis of axial symmetry is z-axis)'
    write(6, *)
    !
    !--------/---------/---------/---------/---------/---------/---------/--
    write(6, *)'chebyshev particles: type a positive number equal to', &
            'the order of the chebyshev polynomial=the number of wrinkles', &
            'on the sphere surface'
    write(6, *)'oblate/prolate spheroids: type -1'
    write(6, *)'oblate/prolate cylinders: type -2'
    write(6, *)'generalized chebyshev particle/droplet: type -3'
    write(6, *)'sphere cut by a plane on its top type -4'
    write(6, *)'sphere cut by a plane on its bottom: type -5'
    write(6, *)'upwardly oriented cone: type -6'
    write(6, *)'conus on a finite cylinder: type -7'
    write(6, *)'nanorod/finite cylinder capped by half-spheroids: -9'
    write(6, *)'intensity around homogeneous/coated sphere: type -50'

    open(unit = 90, file = 'epswater.txt', status = 'unknown')
    !
    np = mpar%np
    if (np == 0) then
        write(6, *) 'input the particle type:'
        read(5, *) np
    else
        write(6, *)'auto-select from *.ini file: ', np
    end if
    !z      np=-1                      !temporarily
    !
    npp = np
    !
    !     np
    !
    !     positive number     chebyshev particles
    !              r(\theta)=r_0[1+\eps*\cos(np*\theta)]
    !     -1                  oblate/prolate spheroids
    !              r(\theta)=a\left[\sin^2\theta + (a^2/b^2)\cos^2\theta]^{-1/2}
    !     -2                  oblate/prolate cylinders
    !     -3                  generalized chebyshev particle/droplet
    !     -4                  sphere cut by a plane on its top
    !     -5                  sphere cut by a plane on its bottom
    !     -6                  cone
    !     -7                  conus on a finite cylinder (in preparation)
    !     -9                  nanorod (cylinder capped with half-spheroids)
    !
    !      parameter (np=-1)
    !
    ! specify the shape of particles within a given np class:
    !     np.gt.0 - defp = deformation parameter of a chebyshev particle
    !     np=-1 - defp = the ratio of the horizontal to rotational axes. defp is
    !             larger than 1 for oblate spheroids and smaller than 1 for
    !             prolate spheroids.
    !     np=-2 - defp = the ratio of the cylinder diameter to its length.
    !     np=-3 - no defp is specified
    !     np=-4 - defp is the height (along the axial symmetry axis)
    !             of the resulting cut sphere
    !     np=-5 - defp is the height (along the axial symmetry axis)
    !             of the resulting cut sphere
    !                note that always defp.lt.2*rev specified
    !     np=-6 - defp is the height (along the axial symmetry axis) divided
    !             by the width of a base
    !     np=-7 - defp  is the height (along the axial symmetry axis) divided
    !             by the width of a base
    !     np=-9 - defp = the ratio of the cylinder diameter to its length.
    !             spheroid caps of the nanorod are defined by additional
    !             parameter.
    !
    ! warning:
    !   in computations for spheres, use defp=1.000001 instead of defp=1.
    !   defp=1 can cause overflows in some rare cases.
    !
    !      parameter (defp=1.000001d0)
    !      defp=1.000001
    !
    !
    !--------/---------/---------/---------/---------/---------/---------/--
    open(unit = 10, file = 'test.txt')
    !
    if (np>0) then
        !
        write(6, *)'radius of the undeformed chebyshev particle in', &
                'your units (in nm if dispersive data used)'
        read(5, *) rsnm
        !
        write(6, *)'the amplitude of wrinkles on the sphere surface'
        read(5, *) defp
        !
        write(6, *)'the equal-volume-sphere radius'
        read(5, *) rev
        !
        rat = 1.d0
        !
    else if (np==-50) then

        if(.not.ynintens) then
            write(6, *)'np.eq.-50 option is only for intensity calculation!!!'
            stop
        end if

        write(6, *)'the (coated) sphere radius in your units'
        read(5, *) rsnm
        rev = rsnm
        !
    else if (np==-1) then

        write(6, *)'the half-length of the spheroid along', &
                'the rotational axis z-axis in your units', &
                '(in nm if dispersive data used)'
        read(5, *) hlength
        !z      hlength=63.3d0

        write(6, *)'the half-length of the spheroid along the', &
                'horizontal axis (in theta=pi/2 plane) in your units', &
                '(in nm if dispersive data used)'
        read(5, *) rsnm
        !z      rsnm=21.1d0
        !
        !     np=-1 - defp = the ratio of the horizontal to rotational axes. defp is
        !             larger than 1 for oblate spheroids and smaller than 1 for
        !             prolate spheroids.
        !
        defp = rsnm/hlength              !always revolution axis length
        !                                    !in the denominator
        !
        rev = rsnm/defp**(1.d0/3.d0)     !=equal-volume-sphere radius
        !                            !room for improvement here - it would be
        !                                  !more resonable to replace rev
        !                              !by the size parameter k*rev
        !
        if (lcs>1) then

            write(6, *)'the ratio of the inner to outer spheroid half-length along the horizontal axis'
            read(5, *) revin

            if (revin>=1) then
                write(6, *)'the ratio cannot be greater or equal 1'
                stop
            end if

            revin = rev*revin
        end if               !lcs.gt.1
        !
        rat = 1.d0
        !
        !--------/---------/---------/---------/---------/---------/---------/--
        !
    else if ((np == -2) .or. (np == -9)) then
        !
        rl_max = mpar%rl_max
        if (rl_max < 0_dp) then
            write(6, *)'enter cylinder maximal r/l:'
            read(5, *) rl_max
        else
            write(6, *)'auto-set cylinder maximal r/l from *.ini', rl_max
        end if

        rl_min = mpar%rl_min
        if (rl_min < 0_dp) then
            write(6, *)'enter cylinder minimal r/l:'
            read(5, *) rl_min
        else
            write(6, *)'auto-set cylinder minimal r/l from *.ini', rl_min
        end if

        ndefp = mpar%ndefp
        if (ndefp <= 0) then
            write(6, *)'enter amount of steps in length:'
            read(5, *) ndefp
        else
            write(6, *)'auto-set amount of steps in length from *.ini', ndefp
        end if

        rsnm = mpar%rsnm
        if (rsnm <= 0_dp) then
            write(6, *)'enter cylinder radius:'
            read(5, *) rsnm
        else
            write(6, *)'auto-set cylinder radius from *.ini', rsnm
        end if

        rsnm = rsnm*2
        ! specify the shape:
        ! np=-2 - defp = the ratio of the cylinder diameter to its length.
        !
        hlength_max = rsnm/rl_min/2.d0
        hlength_min = rsnm/rl_max/2.d0
        defp = rsnm/hlength_max
        rsnm = rsnm/2.d0                                  !cylinder radius
        hlength_max = hlength_max/2.d0
        hlength_min = hlength_min/2.d0        !cylinder half-length
        rev = hlength_max*(3.d0*defp*defp/2.d0)**(1.d0/3.d0)  !=equal-volume-sphere radius

        rat = 1.d0
        !
    else if (np==-3) then
        !     np=-3 - no defp is specified
        write(6, *)'the length of in your units'
        read(5, *) rsnm
        rev = rsnm
        !      write(6,*)'rev(rsnm) not yet determined for np=-3'
        !      stop
        !
        rat = 1.d0
        !
    else if ((np==-4).or.(np==-5))  then
        !     np=-4,-5 - defp is the height (along the axial symmetry axis)
        !                 of the resulting cut sphere
        write(6, *)'the radius of the original uncut sphere in your units'
        !--------/---------/---------/---------/---------/---------/---------/--
        read(5, *) rsnm
        rev = rsnm
        !
        write(6, *)'the height of the cut sphere in your units'
        read(5, *) defp
        defp = defp/rsnm
        !
        rat = 1.d0
        !
    else if (np==-6) then
        !
        !     np=-6 -
        !
        write(6, *)'the cone base diameter in your units'
        read(5, *) rsnm

        write(6, *)'the cone heigth of in your units'
        read(5, *) hlength
        !
        rsnm = rsnm/2.d0
        rev = (hlength*rsnm**2/4.d0)**(1.d0/3.d0)   !=equal-volume-sphere radius
        rat = 1.d0
        !
    else if (np==-7) then
        !
        !     np=-7 -
        !
        write(6, *)'not ready yet!'
        stop
        !c      write(6,*)'the length of in your units'
        !c      read(5,*) rsnm
        !
    end if                         ! end np if
    !***************************************
    !
    defpp = defp
    !
    !
    if (rat==1.) then
        write(6, *)'particle size specified in terms of the equal-volume-sphere radius'
    else if (rat/=1.) then
        write(6, *)'particle size specified in terms of the equal-surface-area-sphere radius'
    end if
    !
    ! equivalent-(volume/surface-area)-sphere radius
    !
    !c      write(6,*)'read equal-volume-sphere radius in nm'
    !c      read(5,*) rev
    !
    !c      rev=300.d0                         !feeded as rev to rsp* routines
    !
    axi = rev
    revf = rev
    !
    !  equivalent equal-(volume/surface-area)-sphere radius
    !
    !c      rev=rat*axi                      !feeded as rev to rsp* routines
    !

    !  ncheck  -  .eq.0  then  ngss=2*ngauss, factor=1d0
    !             .eq.1  then  ngss = ngauss, factor=2d0: theta=pi/2 is mirror
    !                          symmetry plane as in the case of chebysh. particle,
    !                          ellipsoid, and cylinder
    !
    ncheck = 0
    !
    !     !ellipsoid(sphere), cylinder, and nanorod
    if (np==-1.or.np==-2.or.np==-9) ncheck = 1
    if (np>0.and.(-1)**np==1) ncheck = 1    !chebysh. particle
    !
    ! if theta=pi/2 is not a scatterer mirror symmetry plane:
    !  naxsm   -  .eq.0 : gauss abscissas do not have +/- theta symmetry
    !             .eq.1 : gauss abscissas have +/- theta symmetry

    naxsm = 1

    if (np<=-4) naxsm = 0

    !--------/---------/---------/---------/---------/---------/---------/--
    !
    !  alpha and beta - euler angles (in degrees) specifying the orientation
    !    of the scattering particle relative to the laboratory reference
    !    frame (refs. 6 and 7).
    !
    !       alpha=90.d0           !laboratory frame coincides with particle frame
    !      beta=0.d0

    write(6, *)'alpha and beta - euler angles (in degrees) specifying', &
            'the orientation  of the scattering particle relative to the', &
            'laboratory reference  frame'
    !     read(5,*)  alpha, beta
    alpha = mpar%alpha
    beta = mpar%beta
    write(6, *)'auto-set alpha and beta', alpha, beta

    if((alpha==0).and.(beta==0))&
            write(6, *)'laboratory frame coincides with particle frame'
    ! ddelt - the desired absolute accuracy of computing the
    ! expansion coefficients of a normalized scattering matrix.
    ! (this accuracy is usually worse by a factor of 10 than
    ! the accuracy of computing the optical cross sections.)
    ! since convergence test is only performed for the accuracy
    ! of computing the optical cross sections, ddelt is reset
    ! later on to ddelt=0.1d0*ddelt
    !
    ddelt = tol
    !
    !      thet0 - zenith angle of the incident beam in degrees
    !      thet - zenith angle of the scattered beam in degrees
    !      phi0 - azimuth angle of the incident beam in degrees
    !      phi - azimuth angle of the scattered beam in degrees
    !
    if (.not.ynintens) then
        write(6, *)'specify (theta,phi) angles of the incident beam (in degrees)'
        !     read(5,*) thet0,phi0
        thet0 = mpar%thet0
        phi0 = mpar%phi0
        write(6, *)'auto-set from *.ini:', thet0, phi0
        !
        write(6, *)'specify (theta,phi) angles of the scattered beam (in degrees)'
        !     read(5,*) thet,phi
        thet = mpar%thet
        phi = mpar%phi
        write(6, *)'auto-set from *.ini', thet, phi
        !
    else if (ynintens) then
        write(6, *)'phi angle of incidence of the incident plane wave', &
                'is for intensity calculations set to zero'
        write(6, *)'specify theta angle of incidence (in degrees)'
        read(5, *) thet0
        !z      thet0=90.d0
        phi0 = 0.d0
        thetv = thet0
    end if
    !
    !--------/---------/---------/---------/---------/---------/---------/--
    !c      thet0=56.d0
    !c      thet=65.d0
    !c      phi0=114.d0
    !c      phi=128.d0
    !
    ! test setup:

    if ((thet0>180.).or.(thet>180.)) then
        write(6, *)'theta angles has to be smaller than 180'
        stop
    end if

    if ((thet0<0.).or.(thet<0.)) then
        write(6, *)'theta angles has to be positive'
        stop
    end if

    if ((phi0>360.).or.(phi>360.)) then
        write(6, *)'phi angles has to be smaller than 360'
        stop
    end if

    if ((np==-4).or.(np==-5)) then
        if ((defp<=0).or.(defp>=2)) then
            write(6, *)'defp has to be >0 and <2'
            stop
        end if
    end if

    !--------/---------/---------/---------/---------/---------/---------/--
    !
    ! if nag library is available, set mpar%ichoice=1, otherwise mpar%ichoice=2

    mpar%ichoice = 2
    !  controlling the number nd=ndgs*nmax of division points in
    !  computing integrals over the particle surface (ref. 5).
    !  for compact particles, the
    !  recommended value is 2. for highly aspherical particles larger
    !  values (3, 4,...) may be necessary to obtain convergence.
    !  the code does not check convergence over this parameter.
    !  therefore, control comparisons of results obtained with
    !  different ndgs-values are recommended.
    !
    !  check that ndgs*lamxd does not exceed npng1 value in subroutines
    !  for a current values of lmaxd=50 and npng1=800 then ndgs<= 16!!!

    if ((np==-4).or.(np==-5)) then
        ndgs = 16
    else if (np==-6) then
        ndgs = 16
    else if ((np==-1).and.(max(defp, 1.d0/defp)>1.5d0)) then !spheroids
        if (ynintens) then
            ndgs = min(40.d0, 14*max(defp, 1.d0/defp))
        else
            ndgs = min(16.d0, 4*max(defp, 1.d0/defp))
        end if
    else
        ndgs = 4
    end if
    !
    write(6, *) 'ndgs=', ndgs
    !
    if (ynintens) then

        ratp = rat
        ndgsp = ndgs
        naxsmp = naxsm
        ncheckp = ncheck
        revp = rev
        ddeltp = ddelt
        alphae = alpha
        betae = beta

    end if
    !
    if (mpar%ichoice==1) then
        write(6, *) 'lapack routines are used for the matrix inversion'
    else if (mpar%ichoice==2) then
        write(6, *) 'zger and zsur (for the matrix inversion) are used'
    end if
    write(6, *)
    !
    if (ncheck==0) then
        write(6, *) 'particle without theta=pi/2 mirror symmetry'
    else if (ncheck==1) then
        write(6, *) 'particle has  theta=pi/2 mirror symmetry'
    end if
    write(6, *)
    !
    if (naxsm==0) then
        write(6, *) 'gauss abscissas not +/- theta symmetric'
    else if (naxsm==1) then
        write(6, *) 'gauss abscissas +/- theta symmetric'
    end if
    write(6, *)
    !
    write(nout, 5454) mpar%ichoice, ncheck
    5454 format ('ichoice=', i1, '  ncheck=', i1)
    write(nout, *)'naxsm=', naxsm

    if(np==-1.and.defp>=1d0) print 7000, defp
    if(np==-1.and.defp<1d0) print 7001, defp
    if(np>=0) print 7100, np, defp
    if(np==-2.and.defp>=1d0) print 7150, defp
    if(np==-2.and.defp<1d0) print 7151, defp
    if(np==-3) print 7160
    if(np==-4) print 7170, defp
    print 7200, ddelt
    if (dabs(rat - 1d0)<=1d-6) print 8003, axi
    if (dabs(rat - 1d0)>1d-6) print 8004, axi

    7000 format('oblate spheroids, a/b=', f11.7)
    7001 format('prolate spheroids, a/b=', f11.7)
    7100 format('chebyshev particles, t', &
            i1, '(', f5.2, ')')
    7150 format('oblate cylinders, d/l=', f11.7)
    7151 format('prolate cylinders, d/l=', f11.7)
    7160 format('generalized chebyshev particles')
    7170 format('shere cut by a plane, defp=h=', f11.7)
    7200 format ('accuracy of computations ddelt = ', d8.2)
    8003 format('equal-volume-sphere radius=', f8.4)
    8004 format('equal-surface-area-sphere radius=', f8.4)
    !--------/---------/---------/---------/---------/---------/---------/--
    ! checking set up:
    !
    if((np>0).and.(dabs(defp)>=1.d0)) then
        write(6, *)'absolute value of defp has to be less than 1.!!!'
        stop
    end if
    !
    if((np==-4).and.(defp>=2.d0*rev)) then
        write(6, *)'invalid parameters for a cut sphere!'
        write(6, *)'execution stopped!'
        write(6, *)'the defp has to be less than 2*(sphere radius) !!!'
        stop
    end if
    !
    if((np==-1).and.(defp==1)) then
        write(6, *)'use defp=1.000001 instead of defp=1'
    end if
    !
    if (nmat>1) then
        write(6, *)'real material data are to be provided'
        if (ynbrug) write(6, *)'bruggeman approx. used!'
        if (ynbrug) write(nout, *)'#bruggeman approximation performed'
    end if
    !
    if ((ync=='y'.and.lcs==1).or.(ync=='n'.and.lcs/=1)) then
        write(6, *)'check compatibility of ync and lcs'
        stop
    end if

    !--------------------------------------------------------------------
    ! reading in the input data:
    !      write(6,*)'read the particle (core) dielectric constant'
    !      read(5,*) zeps(1)
    !  n(silica)=1.45  <--->    zeps(1)=2.1025d0
    !  n(zns)=2.       <--->    zeps(1)=4.d0
    zeps(1) = cceps
    if(lcs>1) zeps(lcs) = cseps
    zeps(lcs + 1) = zeps0
    !*****************************
    !
    if (lcs>=2) then                  !coated particle
        !
        !c      write(6,*)'read equal-volume-core radius in nm'
        !c      read(5,*) rff(1)
        !c      rff(1)=204.9d0
        !c      rff(1)=rff(1)/rs
        !
        if (np==-50) then
            !
            write(6, *)'coated sphere core radii r(l) labelled from 1 for the', &
                    'inner core till lcs for the outer shell', &
                    '===> r(lcs) is the sphere radius'
            !--------/---------/---------/---------/---------/---------/---------/--
            do ikl = 1, lcs - 1
                !
                write(6, *)'read in r(l) for l=', ikl
                read (5, *) rff(ikl)
                rff(ikl) = rff(ikl)/rsnm

                !      rff(1)=0.75d0
                rmf(ikl) = rff(ikl)*rmuf

                if ((lcs>2).and.(ikl>=2)) then
                    write(6, *)'read in the lth-sphere layer diel. const. for l=', ikl
                    write(6, *)'(in case of dispersive component, give 1. )'
                    read (5, *) zeps(ikl)
                    !
                end if         !lcs.ikl

            enddo          !ikl
            !
        end if         !np.eq.-50
        !
    end if         !lcs.ge.2
    !
    !oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
    ! scanning over frequency interval:
    !-------------------------------------------------
    !
    if (ynintens) then
        write(6, *)'electric field intensity calculation'
        write(6, *)'read lambda (in the vacuum) in nm'
    else if (.not.ynintens) then
        write(6, *)'read initial (minimal) x-parameter'
    end if
    !--------/---------/---------/---------/---------/---------/---------/--
    x_min = mpar%x_min
    if (x_min < 0) then
        write(6, *)'enter minimum value of x-parameter:'
        read(5, *) x_min
    else
        write(6, *)'auto-set x_min from *.ini ', x_min
    end if

    lambda = 2*pi*rsnm/x_min
    !v      lambda=633.d0                      !jap89_5776 for gold ellipsod
    !      lambda=354.d0                      !jap89_5776 for silver sphere
    !c      lambda=500d0
    !
    ! size parameter is customarily defined as the ratio of
    ! circumference of particle to the wavelength in the host medium
    ! in which the particle is embedded
    !                    x=kr=\sg a=2*pi*r/lambda
    !      xs=2.d0*pi*rs*dble(sqrt(zeps0))/lambda
    ! convert lambda to the lambda in vacuum:
    !      lambda=lambda*dble(sqrt(zeps0))
    !  omega=2.d0*pi*rsnm/(lambda*rmuf)=xs/(rmuf*dble(sqrt(zeps0))),
    ! where rs is the particle radius (in nm) and  lambda
    ! is the wavelengths (in nm)
    ! in the vacuum:

    omega = 2.d0*pi*rev/(lambda*rmuf)
    omega0 = omega

    !      write(6,*)'read omega'
    !      read(5,*) omega
    !      write(6,*)'omega=', omega
    !      if (1.eq.1) nstep=0             ! temporarily
    !      delo=0.d0
    !      omega0=omega
    !      goto 11
    !
    ! option for omega input:
    !      write(6,*)'read omega ='
    !      read(5,*) omega
    !
    !      xs=rmuf*omega*dble(sqrt(zeps0))
    !
    !       write(6,*)'equiv. size parameter x=2*pi*rs*n_0/lambda=',xs
    !
    if (.not.ynintens) then
        write(6, *)'scan up to x-parameter (in nm)'
        x_max = mpar%x_max
        if (x_max < 0) then
            write(6, *)'enter maximum value of x-parameter:'
            read(5, *) x_max
        else
            write(6, *)'auto-set x_max from *.ini ', x_max
        end if

        enw = 2*pi*rsnm/x_max
        !      enw=500
        nstep = mpar%nstep
        if (nstep <= 0) then
            write(6, *)'enter amount of scanning steps:'
            read(5, *) nstep
        else
            write(6, *)'auto-set nstep from *.ini ', nstep
        end if
        !      xstep=5
    end if
    !
    ! ::: number of steps on frequency interval:
    if ((lambda==enw).or.(ynintens)) then
        xstep = 0
        delo = 0.d0
    else
        xstep = (lambda - enw)/nstep
        ! ::: width of the searched frequency interval
        enw = 2.d0*pi*rev/(enw*rmuf)
        enw = enw - omega0
        delo = enw/dble(nstep - 1)
    end if
    !
    !--------/---------/---------/---------/---------/---------/---------/--
    !                  --------------------------------
    ! output initial statements

    open(unit = nout, file = 'axs-scs.dat')
    rewind(nout)
    write(nout, *)'#orientationally averaged scattering cs for a single particle'
    write(nout, *)'#(cross sections normalized per equal-volume-sphere surface s=pi*rev**2)'
    write(nout, *)
    write(nout, *)'#equiv.-volume-sphere radius in nm=', rev
    write(nout, *)
    write(nout, *)'#angular momentum cut-off lmax=', lmax
    write(nout, *)'#in the number nd=ndgs*lmax of gip points, ndgs=', ndgs
    write(nout, *)'#defp=', defp
    write(nout, *)'#material number =', nmat
    if (ynbrug) write(nout, *)'#bruggeman approximation performed'
    if (ync =='n') then
        write(nout, *)'#homogeneous particle'
    else if (ync=='y') then
        write(nout, *)'#coated particle'
    end if
    write(nout, *)
    write(nout, *)'#host dielectric constant=', zeps(lcs + 1)
    if (nmat>=1) write(nout, *)'#dispersive layer number=', ilcs
    if ((lcs==1).and.(nmat==0))&
            write(nout, *)'#sphere diel. const.=', zeps(1)
    if (lcs>=2) then
        write(nout, *)'#core radius/sphere radius =', rff(1)
        if ((ilcs/=1).or.(nmat==0))&
                write(nout, *)'#sphere core diel. const.=', zeps(1)
        if ((ilcs/=lcs).or.(nmat==0))&
                write(nout, *)'#coating diel. const.=', zeps(lcs)
        !--------/---------/---------/---------/---------/---------/---------/--
    end if
    write(nout, *)'#vacuum lambda_0 and sg_sc in columns'
    write(nout, *)

    open(unit = nout + 1, file = 'axs-oavext.dat')
    rewind(nout + 1)
    write(nout + 1, *)'#orientationally averaged extinction cs for a single (coated) particle'
    write(nout + 1, *)'#(cross sections normalized per equal-volume-sphere surface s=pi*rev**2)'
    write(nout + 1, *)
    write(nout + 1, *)'#equiv.-volume-sphere radius in nm=', rev
    write(nout + 1, *)
    write(nout + 1, *)'#angular momentum cut-off lmax=', lmax
    write(nout + 1, *)'#in the number nd=ndgs*lmax of gip points, ndgs=', ndgs
    write(nout + 1, *)'#defp=', defp
    write(nout + 1, *)'#material number =', nmat
    if (ynbrug) write(nout + 1, *)'#bruggeman approximation performed'
    if (ync =='n') then
        write(nout + 1, *)'#homogeneous particle'
    else if (ync=='y') then
        write(nout + 1, *)'#coated particle'
    end if
    write(nout + 1, *)
    write(nout + 1, *)'#host dielectric constant=', zeps(lcs + 1)
    if (nmat>=1) write(nout + 1, *)'#dispersive layer number=', ilcs
    if ((lcs==1).and.(nmat==0))&
            write(nout + 1, *)'#sphere diel. const.=', zeps(1)
    if (lcs>=2) then
        write(nout + 1, *)'#core radius/sphere radius =', rff(1)
        if ((ilcs/=1).or.(nmat==0))&
                write(nout + 1, *)'#sphere core diel. const.=', zeps(1)
        if ((ilcs/=lcs).or.(nmat==0))&
                write(nout + 1, *)'#coating diel. const.=', zeps(lcs)
        !--------/---------/---------/---------/---------/---------/---------/--
    end if
    write(nout + 1, *)'#vacuum lambda_0 and sg_tot in columns'
    write(nout + 1, *)

    open(unit = nout + 2, file = 'axs-abs.dat')
    rewind(nout + 2)
    write(nout + 2, *)'#orientationally averaged absorption cs for a single (coated) particle'
    write(nout + 2, *)'#(cross sections normalized per equal-volume-sphere surface s=pi*rev**2)'
    write(nout + 2, *)
    write(nout + 2, *)'#equiv.-volume-sphere radius in nm=', rev
    write(nout + 2, *)
    write(nout + 2, *)'#angular momentum cut-off lmax=', lmax
    write(nout + 2, *)'#in the number nd=ndgs*lmax of gip points, ndgs=', ndgs
    write(nout + 2, *)'#defp=', defp
    write(nout + 2, *)'#material number =', nmat
    if (ynbrug) write(nout + 2, *)'#bruggeman approximation performed'

    if (ync =='n') then
        write(nout + 2, *)'#homogeneous particle'
    else if (ync=='y') then
        write(nout + 2, *)'#coated particle'
    end if
    write(nout + 2, *)
    write(nout + 2, *)'# host dielectric constant=', zeps(lcs + 1)
    if (nmat>=1) write(nout + 2, *)'#dispersive layer number=', ilcs
    if ((lcs==1).and.(nmat==0))&
            write(nout + 2, *)'#sphere diel. const.=', zeps(1)
    if (lcs>=2) then
        write(nout + 2, *)'#core radius/sphere radius =', rff(1)
        if ((ilcs/=1).or.(nmat==0))&
                write(nout + 2, *)'#sphere core diel. const.=', zeps(1)
        if ((ilcs/=lcs).or.(nmat==0))&
                write(nout + 2, *)'#coating diel. const.=', zeps(lcs)
        !--------/---------/---------/---------/---------/---------/---------/--
    end if
    write(nout + 2, *)'#vacuum lambda_0 and sg_abs in columns'
    write(nout + 2, *)

    !--------/---------/---------/---------/---------/---------/---------/--
    open(unit = nout + 3, file = 'axs-ext.dat')
    rewind(nout + 3)
    write(nout + 3, *)'#extinction cs for a single (coated) particle in a fixed orientation'
    write(nout + 3, *)'#(cross sections normalized per equal-volume-sphere surface s=pi*rev**2)'
    write(nout + 3, *)
    write(nout + 3, *)'#equiv.-volume-sphere radius in nm=', rev
    write(nout + 3, *)
    write(nout + 3, *)'#angular momentum cut-off lmax=', lmax
    write(nout + 3, *)'#in the number nd=ndgs*lmax of gip points, ndgs=', ndgs
    write(nout + 3, *)'#defp=', defp
    write(nout + 3, *)'#material number =', nmat
    if (ynbrug) write(nout + 3, *)'#bruggeman approximation performed'
    if (ync =='n') then
        write(nout + 3, *)'#homogeneous particle'
    else if (ync=='y') then
        write(nout + 3, *)'#coated particle'
    end if
    write(nout + 3, *)
    write(nout + 3, *)'#host dielectric constant=', zeps(lcs + 1)
    if (nmat>=1) write(nout + 3, *)'#dispersive layer number=', ilcs
    if ((lcs==1).and.(nmat==0))&
            write(nout + 3, *)'#sphere diel. const.=', zeps(1)
    if (lcs>=2) then
        write(nout + 3, *)'#core radius/sphere radius =', rff(1)
        if ((ilcs/=1).or.(nmat==0))&
                write(nout + 3, *)'#sphere core diel. const.=', zeps(1)
        if ((ilcs/=lcs).or.(nmat==0))&
                write(nout + 3, *)'#coating diel. const.=', zeps(lcs)
        !--------/---------/---------/---------/---------/---------/---------/--
    end if
    write(nout + 3, *)'#in columns: vacuum lambda_0 and sg_tot (normalized and unnormalized)'
    write(nout + 3, *)

    open(unit = nout + 5, file = 'axs-albedo.dat')
    rewind(nout + 5)
    write(nout + 5, *)'#orientationally averaged albedo for a single (coated) particle'
    write(nout + 5, *)'#equiv.-volume-sphere radius in nm=', rev
    write(nout + 5, *)
    write(nout + 5, *)'#angular momentum cut-off lmax=', lmax
    write(nout + 5, *)'#in the number nd=ndgs*lmax of gip points, ndgs=', ndgs
    write(nout + 5, *)'#defp=', defp
    write(nout + 5, *)'#material number =', nmat
    if (ynbrug) write(nout + 5, *)'#bruggeman approximation performed'
    if (ync =='n') then
        write(nout + 5, *)'#homogeneous particle'
    else if (ync=='y') then
        write(nout + 5, *)'#coated particle'
    end if
    write(nout + 5, *)
    write(nout + 5, *)'#host dielectric constant=', zeps(lcs + 1)
    if (nmat>=1) write(nout + 5, *)'#dispersive layer number=', ilcs
    if ((lcs==1).and.(nmat==0))&
            write(nout + 5, *)'#sphere diel. const.=', zeps(1)
    if (lcs>=2) then
        write(nout + 5, *)'#core radius/sphere radius =', rff(1)
        if ((ilcs/=1).or.(nmat==0))&
                write(nout + 5, *)'#sphere core diel. const.=', zeps(1)
        if ((ilcs/=lcs).or.(nmat==0))&
                write(nout + 5, *)'#coating diel. const.=', zeps(lcs)
        !--------/---------/---------/---------/---------/---------/---------/--
    end if
    write(nout + 5, *)'#in columns:'
    write(nout + 5, *)'#vacuum lambda_0, albedo, and tcs-(acs+tsc)'
    write(nout + 5, *)

    open(unit = nout + 6, file = 'axs-dipolext.dat')
    rewind(nout + 6)
    write(nout + 6, *)'#orientationally averaged dipole ext. for a single (coated) particle'
    write(nout + 6, *)'#equiv.-volume-sphere radius in nm=', rev
    write(nout + 6, *)
    write(nout + 6, *)'#angular momentum cut-off lmax=', lmax
    write(nout + 6, *)'#in the number nd=ndgs*lmax of gip points, ndgs=', ndgs
    write(nout + 6, *)'#defp=', defp
    write(nout + 6, *)'#material number =', nmat
    if (ync =='n') then
        write(nout + 6, *)'#homogeneous particle'
    else if (ync=='y') then
        write(nout + 6, *)'#coated particle'
    end if
    write(nout + 6, *)
    write(nout + 6, *)'#host dielectric constant=', zeps(lcs + 1)
    if (nmat>=1) write(nout + 6, *)'#dispersive layer number=', ilcs
    if ((lcs==1).and.(nmat==0))&
            write(nout + 6, *)'#sphere diel. const.=', zeps(1)
    if (lcs>=2) then
        write(nout + 6, *)'#core radius/sphere radius =', rff(1)
        if ((ilcs/=1).or.(nmat==0))&
                write(nout + 6, *)'#sphere core diel. const.=', zeps(1)
        if ((ilcs/=lcs).or.(nmat==0))&
                write(nout + 6, *)'#coating diel. const.=', zeps(lcs)
        !--------/---------/---------/---------/---------/---------/---------/--
    end if
    write(nout + 6, *)'#vacuum lambda_0 and dipole extinction in columns'
    write(nout + 6, *)
    !--------/---------/---------/---------/---------/---------/---------/--
    open(unit = nout + 7, file = 'axs-quadrext.dat')
    rewind(nout + 7)
    write(nout + 7, *)'#orientationally averaged quadrupole ext. for a single (coated) particle'
    write(nout + 7, *)'#equiv.-volume-sphere radius in nm=', rev
    write(nout + 7, *)
    write(nout + 7, *)'#angular momentum cut-off lmax=', lmax
    write(nout + 7, *)'#in the number nd=ndgs*lmax of gip points, ndgs=', ndgs
    write(nout + 7, *)'#defp=', defp
    write(nout + 7, *)'#material number =', nmat
    if (ynbrug) write(nout + 7, *)'#bruggeman approximation performed'
    if (ync =='n') then
        write(nout + 7, *)'#homogeneous particle'
    else if (ync=='y') then
        write(nout + 7, *)'#coated particle'
    end if
    write(nout + 7, *)
    write(nout + 7, *)'#host dielectric constant=', zeps(lcs + 1)
    if (nmat>=1) write(nout + 7, *)'#dispersive layer number=', ilcs
    if ((lcs==1).and.(nmat==0))&
            write(nout + 7, *)'#sphere diel. const.=', zeps(1)
    if (lcs>=2) then
        write(nout + 7, *)'#core radius/sphere radius =', rff(1)
        if ((ilcs/=1).or.(nmat==0))&
                write(nout + 7, *)'#sphere core diel. const.=', zeps(1)
        if ((ilcs/=lcs).or.(nmat==0))&
                write(nout + 7, *)'#coating diel. const.=', zeps(lcs)
        !--------/---------/---------/---------/---------/---------/---------/--
    end if
    write(nout + 7, *)'#in columns:'
    write(nout + 7, *)'#vacuum lambda_0 and quadrupole extinction'
    write(nout + 7, *)
    !--------/---------/---------/---------/---------/---------/---------/--
    open(unit = nout + 10, file = 'axs-phmat.dat')
    rewind(nout + 10)
    write(nout + 10, 5000)
    write(nout + 10, *)
    write(nout + 10, *)'#equiv.-volume-sphere radius in nm=', rev
    write(nout + 10, *)
    write(nout + 10, *)'#angular momentum cut-off lmax=', lmax
    write(nout + 10, *)'#in the number nd=ndgs*lmax of gip points, ndgs=', ndgs
    write(nout + 10, *)'#defp=', defp
    write(nout + 10, *)'#material number =', nmat
    if (ynbrug) write(nout + 10, *)'#bruggeman approximation performed'
    write(nout + 10, 1005) thet0, thet, phi0, phi, alpha, beta
    write(nout + 10, *)

    !--------/---------/---------/---------/---------/---------/---------/--
    open(unit = nout + 12, file = 'axs-ampmat.dat')
    rewind(nout + 12)
    write(nout + 12, 1006)
    write(nout + 12, *)
    write(nout + 12, *)'#equiv.-volume-sphere radius in nm=', rev
    write(nout + 12, *)
    write(nout + 12, *)'#angular momentum cut-off lmax=', lmax
    write(nout + 12, *)'#in the number nd=ndgs*lmax of gip points, ndgs=', ndgs
    write(nout + 12, *)'#defp=', defp
    write(nout + 12, *)'#material number =', nmat
    if (ynbrug) write(nout + 12, *)'#bruggeman approximation performed'
    write(nout + 12, 1005) thet0, thet, phi0, phi, alpha, beta
    write(nout + 12, *)'#in columns:'
    write(nout + 12, *)'#   vacuum lambda_0, revv, imvv, revh, imvh'
    write(nout + 12, *)'#   vacuum lambda_0, rehv, imhv, rehh, imhh'
    write(nout + 12, *)
    !--------/---------/---------/---------/---------/---------/---------/--
    open(unit = nout + 13, file = 'axs-intv.dat')
    rewind(nout + 13)
    write(nout + 13, *)'intensity in theta for (un)correlated light source'
    write(nout + 13, *)
    write(nout + 13, *)'#equiv.-volume-sphere radius in nm=', rev
    write(nout + 13, *)
    write(nout + 13, *)'#angular momentum cut-off lmax=', lmax
    write(nout + 13, *)'#in the number nd=ndgs*lmax of gip points, ndgs=', ndgs
    write(nout + 13, *)'#defp=', defp
    write(nout + 13, *)'#material number =', nmat
    if (ynbrug) write(nout + 13, *)'#bruggeman approximation performed'
    write(nout + 13, 1005) thet0, thet, phi0, phi, alpha, beta
    write(nout + 13, *)'#host dielectric constant=', zeps(lcs + 1)
    if (nmat>=1) write(nout, *)'#dispersive layer number=', ilcs
    if ((lcs==1).and.(nmat==0))&
            write(nout + 13, *)'#sphere diel. const.=', zeps(1)
    if (lcs>=2) then
        write(nout + 13, *)'#core radius/sphere radius =', rff(1)
        if ((ilcs/=1).or.(nmat==0))&
                write(nout + 13, *)'#sphere core diel. const.=', zeps(1)
        if ((ilcs/=lcs).or.(nmat==0))&
                write(nout + 13, *)'#coating diel. const.=', zeps(lcs)
        !--------/---------/---------/---------/---------/---------/---------/--
    end if
    write(nout + 13, *)'#in columns:'
    write(nout + 13, *)'#vacuum lambda_0,|vv|**2+|vh|**2,|vv+vh|**2'
    write(nout + 13, *)

    !--------/---------/---------/---------/---------/---------/---------/--
    open(unit = nout + 14, file = 'axs-inth.dat')
    rewind(nout + 14)
    write(nout + 14, *)'intensity in phi for (un)correlated light source'
    write(nout + 14, *)
    write(nout + 14, *)
    write(nout + 14, *)'#equiv.-volume-sphere radius in nm=', rev
    write(nout + 14, *)
    write(nout + 14, *)'#angular momentum cut-off lmax=', lmax
    write(nout + 14, *)'#in the number nd=ndgs*lmax of gip points, ndgs=', ndgs
    write(nout + 14, *)'#defp=', defp
    write(nout + 14, *)'#material number =', nmat
    if (ynbrug) write(nout + 14, *)'#bruggeman approximation performed'
    write(nout + 14, 1005) thet0, thet, phi0, phi, alpha, beta
    write(nout + 14, *)'#host dielectric constant=', zeps(lcs + 1)
    if (nmat>=1) write(nout, *)'#dispersive layer number=', ilcs
    if ((lcs==1).and.(nmat==0))&
            write(nout + 14, *)'#sphere diel. const.=', zeps(1)
    if (lcs>=2) then
        write(nout + 14, *)'#core radius/sphere radius =', rff(1)
        if ((ilcs/=1).or.(nmat==0))&
                write(nout + 14, *)'#sphere core diel. const.=', zeps(1)
        if ((ilcs/=lcs).or.(nmat==0))&
                write(nout + 14, *)'#coating diel. const.=', zeps(lcs)
        !--------/---------/---------/---------/---------/---------/---------/--
    end if
    write(nout + 14, *)'#in columns:'
    write(nout + 14, *)'#vacuum lambda_0, |hv|**2+|hh|**2, |hv+hh|**2'
    write(nout + 14, *)

    open(unit = nout + 15, file = 'tr1diag.dat')
    rewind(nout + 15)

    open(10, file = 'maps.txt', status = 'unknown')
    !****************************   zeps1  ***********************************
    !                  --------------------------------
    ! sphere optical constants in the case of a dispersion
    ! reading in material  data:
    ! reading real material data, e.g., according to palik's  book
    ! requires reading data files omf and ceps1 of dimension nfin
    ! omf is reepsz/omega and ceps1 contains the sphere eps
    !                       material constant reading:
    !
    select case (nmat)
    case (1) !no reading of material data
    case (2) ! silver data

        open(unit = 30, file = 'agc.dat')
        rewind(30)
        do ieps = 1, nfin
            read(30, *) omf(ieps), ceps1(ieps)
        enddo
        close(30)

    case (3) ! gold data

        !      open(unit=30,file='au293knew.dat')       !gold data for different t
        open(unit = 30, file = 'audat.dat')          !gold data in nm
        write(6, *)'gold particles'
        rewind(30)
        do ieps = 1, nfin
            read(30, *) omf(ieps), ceps1(ieps)
            !          omf(ieps)=2.d0*pi*rev*omf(ieps)/(1240.d0*rmuf)
            omf(ieps) = 2.d0*pi*rev/(omf(ieps)*rmuf)
        enddo
        close(30)

        !c      else if (nmat.eq.4) then

    case (5) ! copper data

        open(unit = 30, file = 'cudat.dat')          !copper data in nm
        write(6, *)'copper particles'
        rewind(30)
        do ieps = 1, nfin
            read(30, *) omf(ieps), ceps1(ieps)
            omf(ieps) = 2.d0*pi*rev/(omf(ieps)*rmuf)
        enddo
        close(30)

    case (6) ! aluminium data

        open(unit = 30, file = 'aldat.dat')          !aluminium data in nm
        write(6, *)'aluminum particles'
        rewind(30)
        do ieps = 1, nfin
            read(30, *) omf(ieps), ceps1(ieps)
            omf(ieps) = 2.d0*pi*rev/(omf(ieps)*rmuf)
        enddo
        close(30)

    case (7) ! platinum data

        open(unit = 30, file = 'ptdat.dat')          !platinum data in nm
        write(6, *)'platinum particles'
        rewind(30)
        do ieps = 1, nfin
            read(30, *) omf(ieps), ceps1(ieps)
            omf(ieps) = 2.d0*pi*rev/(omf(ieps)*rmuf)
        enddo
        close(30)

    case (8) ! silicon data

        !     open(unit=30,file='sieps.dat')  !silicon data in nm
        open(unit = 30, file = 'sidat.dat')   !silicon data in nm for larger interval
        write(6, *)'silicon particles'
        rewind(30)
        do ieps = 1, nfin
            read(30, *) omf(ieps), ceps1(ieps)
            !          ceps1imag(ieps)= aimag(ceps1(ieps))
            !         ceps1real(ieps) = real(ceps1(ieps))
            !          ceps1imag(ieps)=0
            !          ceps1(ieps) = cmplx(ceps1real(ieps), ceps1imag(ieps))
            omf(ieps) = 2.d0*pi*rev/(omf(ieps)*rmuf)
        enddo
        close(30)

    case (9) ! water data

        open(unit = 30, file = 'measured_water_dispersion_t=24.txt')   !water in ghz
        write(6, *)'water particles'
        rewind(30)
        do ieps = 1, nfin
            read(30, *) omf(ieps), ceps1real(ieps), ceps1imag(ieps)
            ceps1(ieps) = cmplx(ceps1real(ieps), ceps1imag(ieps))
            omf(ieps) = 2.d0*pi*rev/((3.d8/omf(ieps))*rmuf)
        enddo
        close(30)

    end select ! material constant reading
    !********************
    !                     --------------------------------
    ! begin main scanning loop:
    write(6, *) 'begining of main loop'
    rev_beg = rev

    do itter = 1, ndefp

        defp = rsnm/hlength_max + (rsnm/hlength_min - rsnm/hlength_max)*&
                dble(itter - 1)/dble(ndefp - 1)
        defpp = defp
        hlength = rsnm/defp
        rev = hlength*(3.d0*defp*defp/2.d0)**(1.d0/3.d0)
        do istep = 1, nstep

            omega = omega0 + dble(istep - 1)*delo

            write(6, *) itter, istep, defp, hlength
            !omega_max = omega0 + dble(nstep)*delo
            !lambda_min=2.d0*pi*rev/(omega_max*rmuf)
            !lambda_max=2.d0*pi*rev/(omega0*rmuf)
            !lambda = lambda_min + (lambda_max - lambda_min)*dble(istep) &
            !         /dble(nstep)
            lambda = 2.d0*pi*rev_beg/(omega*rmuf)

            !xs=rmuf*omega*dble(sqrt(zeps0))  !equiv. size parameter

            ! it is neither dispersionless dielectric nor ideal metal
            if ((nmat/=0).and.(.not.ynperfcon)) then
                ! in case of a dispersion, epssph is modified.
                ! for ideal drude metal
                !     plasma=2.d0*pi*sphere radius in nm/(lambda_z in nm*rmuf)
                ! where lambda_z is the wavelength for which re eps_s=0.

                reepsz = 2.d0*pi*rev/(323.83d0*rmuf)

                select case (nmat)
                case (1) !material decision if - drude metal

                    plasma = reepsz
                    omxp = plasma/omega
                    zeps1 = 1.d0 - omxp**2/(1.d0 + ci*plasma/(144.d0*omega))

                case (4) !material decision if - zns
                    !
                    filfrac = 0.62d0         ! filfrac of zns in zns core
                    call  znsrefind(lambda, filfrac, zeps1)

                case (2) !material decision if - ag

                    ! >>> real material data:           !silver
                    !                         lambda_z=323.83d0
                    !                         lambda_p=164.d0
                    ! when real material data are used,
                    ! reepsz differs from plasma!!! the plasma wavelength is
                    ! calculated below:

                    plasma = reepsz*7.2d0/3.8291d0
                    ! security trap - remainder (not optimized!)
                    omxf = omega/reepsz
                    if (omxf>omf(1)) then
                        write(6, *)'calculation of has to stop with'
                        write(6, *)' omf(1)'
                        write(6, *)' omxf=', omxf
                        stop
                    end if

                    if (omxf<omf(nfin)) then
                        omxp = plasma/omega
                        zeps1 = 1.d0 - omxp**2/(1.d0 + ci*plasma/(144.d0*omega))
                        ! damping coefficient for silver is plasma/144 where plasma is different from
                        ! the re eps zero crossing at 3.8291 ev according to palik!!!
                    else if (omxf==omf(1)) then
                        zeps1 = ceps1(1)
                    else
                        do ieps = 2, nfin
                            ! data file ordered with the increased wavelength
                            ! omxf increases in the loop and is oriented opposite to the data file
                            if (omxf>omf(ieps)) then     ! linear interpolation
                                zeps1 = ceps1(ieps) + (omxf - omf(ieps))*(ceps1(ieps - 1) - ceps1(ieps))&
                                       /(omf(ieps - 1) - omf(ieps))
                                exit
                            end if
                        enddo
                    end if   !end ag

                case (3, 5, 7)
                    !material decision if&
                    !au,cu,al,pt
                    ! >>>
                    ! data file ordered with the decreased wavelength
                    ! omega increases in the loop and is oriented along the data file
                    !
                    if ((omega<omf(1)).or.(omega>omf(nfin))) then
                        !c       write(6,*)'material data not available for this wavelength'
                        !c       stop

                        call sordalc(nmat, lambda, zeps1)

                    end if
                    !
                    if (omega==omf(nfin)) then
                        zeps1 = ceps1(nfin)
                    else
                        do ieps = 1, nfin - 1
                            if (omega<omf(ieps + 1)) then     ! linear interpolation
                                zeps1 = ceps1(ieps) + (omega - omf(ieps))*(ceps1(ieps + 1) - ceps1(ieps))&
                                       /(omf(ieps + 1) - omf(ieps))
                                exit
                            end if
                        enddo
                    end if

                case (8) !material decision if - silicon
                    ! >>>
                    ! data file ordered with the decreased wavelength
                    ! omega increases in the loop and is oriented along the data file
                    !
                    if ((omega<omf(1)).or.(omega>omf(nfin))) then
                        write(6, *)'material data not available for this wavelength'
                        stop
                        !
                    end if
                    !
                    if (omega==omf(nfin)) then
                        zeps1 = ceps1(nfin)
                    else
                        do ieps = 1, nfin - 1
                            if (omega<omf(ieps + 1)) then     ! linear interpolation
                                zeps1 = ceps1(ieps) + (omega - omf(ieps))*(ceps1(ieps + 1) - ceps1(ieps))&
                                       /(omf(ieps + 1) - omf(ieps))
                                exit
                            end if
                        enddo
                    end if

                case (9) !material decision if - water
                    ! >>>
                    ! data file ordered with the decreased wavelength
                    ! omega increases in the loop and is oriented along the data file
                    !
                    if ((omega<omf(1)).or.(omega>omf(nfin))) then
                        write(6, *)'material data not available for this wavelength'
                        stop
                        !
                    end if
                    !
                    if (omega==omf(nfin)) then
                        zeps1 = ceps1(nfin)
                    else
                        do ieps = 1, nfin - 1
                            if (omega<omf(ieps + 1)) then     ! linear interpolation
                                zeps1 = ceps1(ieps) + (omega - omf(ieps))*(ceps1(ieps + 1) - ceps1(ieps))&
                                       /(omf(ieps + 1) - omf(ieps))
                                exit
                            end if
                        enddo
                    end if

                end select ! end of material decision if
                ! the end of reading real data according to palik's  book
                !_____________________________________
                ! activate bruggeman:

                if (ynbrug) then
                    ff = 0.8d0
                    z1 = (3.d0*ff - 1.d0)*zeps1 + (2.d0 - 3.d0*ff)*zeps0
                    z2 = sqrt(z1*z1 + 8.d0*zeps1*zeps0)
                    !
                    if (aimag(z2)>=0.0) then
                        zeps1 = (z1 + z2)/4.d0
                    else
                        zeps1 = (z1 - z2)/4.d0
                    end if
                end if

                zeps(ilcs) = zeps1
                !______________________________________
            end if

            !c      write(6,*)'lambda in axspartcl=', lambda

            if (ynintens) goto 500
            !
            if ((np==-1).and.(lcs>1)) then

                revinl = revin*2.d0*pi*sqrt(zeps0)/lambda

                revl = rev*2.d0*pi*sqrt(zeps0)/lambda

                !test
                !c      defpp=0.5
                !c      revinl=0.5
                !c      revl=10.d0
                !test
                if (mpar%yncheck /= 0) then
                    ide = 2
                else
                    ide = 4
                end if

                ! todo: commented out lisac just to compile, call has not enough arguments.
                !subroutine lisac(de,lcs,nmx1,nmx2,ngauss,lambda,eps,ki,kex,
                !     & thet0,thet,phi0,phi,zeps)
                !     call lisac(ide,lcs,lmax,lmax,ndgs*lmax,lambda,defpp,
                !    & revinl/revl,revl,zeps)

            else
                global_eps_r = real(zeps(1))
                global_eps_i = aimag(zeps(1))


                !      write(90,*) lambda,  global_eps_r,
                !     & global_eps_i
                call ampldr(lmax, npp, defpp, &
                        rsnm, hlength, lambda, zeps(1), zeps0)
            end if
            !
            !--------/---------/---------/---------/---------/---------/---------/--
            !c      if(nstep.gt.10) goto 200
            !      write(6,*) 'istep=', istep
            !c      write(6,*) 'lambda=',lambda
            !c      write(6,*)
            !c      write(6,*)'scattering coefficient='
            !c      write(6,*)'tsc=', tsc

        end do
    end do
    ! <<<
    close(nout)
    close(nout + 1)
    close(nout + 2)
    close(nout + 3)
    close(nout + 5)
    close(nout + 6)
    close(nout + 7)
    close(nout + 10)
    close(nout + 12)
    close(nout + 13)
    close(nout + 14)
    close(nout + 15)
    close(10)
    close(90)
    ! <<<
    if (ync =='n') then
        write(6, *)'homogeneous particle'
    else if (ync=='y') then
        write(6, *)'coated particle'
    end if
    !
    write(6, *)'particle parameters:'
    write(6, *)
    write(6, *)'equivalent sphere radius =', rev
    if (ync=='n') write(6, *)'particle diel. constant=', zeps(1)
    write(6, *)'background dielectric constant zeps0=', zeps0
    if (ync=='y') write(6, *)'core diel. constant=', zeps(1)
    if (ync=='y') write(6, *)'coating diel. constant=', zeps(lcs)
    if (ync=='y') write(6, *)'core radius/particle radius =', rff(1)
    write(6, *)
    write(6, *)'oa scattering cs versus wavelength in axs-scs.dat'
    write(6, *)'oa extinction versus wavelength in axs-oavext.dat'
    write(6, *)'oa absorption versus wavelength in axs-abs.dat'
    write(6, *)'oa albedo versus wavelength in axs-albedo.dat'
    write(6, *)'  [3rd column displays qext-(qsca +qabs)]  '
    write(6, *)'extinction versus wavelength in axs-ext.dat'
    write(6, *)'phase matrix vs wavelength in axsphmat.dat'
    write(6, *)'amplitude matrix vs wavelength in axsampmat.dat'
    write(6, *)'oa dipole extinction in axs-dipolext.dat'
    write(6, *)'oa quadrupole extinction in axs-quadrext.dat'
    !--------/---------/---------/---------/---------/---------/---------/--

    if (.not.ynintens) goto 1000

    !--------/---------/---------/---------/---------/---------/---------/--
    !
    ! output initial statements

    500  open(unit = nouti, file = 'intnsty.dat')
    rewind(nouti)
    write(nouti, *)'#field intensity profile around a tip'
    write(nouti, *)
    write(nouti, *)'#equiv.-volume-sphere radius in nm=', rev
    write(nouti, *)
    write(nouti, *)'#angular momentum cut-off lmax=', lmax
    write(nouti, *)'#in the number nd=ndgs*lmax of gip points, ndgs=', ndgs
    write(nouti, *)'#defp=', defp
    write(nouti, *)'#material number =', nmat
    if (ync =='n') then
        write(nouti, *)'#homogeneous particle'
    else if (ync=='y') then
        write(nouti, *)'#coated particle'
    end if
    write(nouti, *)
    write(nouti, *)'#host dielectric constant=', zeps(lcs + 1)
    if (ync=='n')&
            write(nouti, *)'#particle diel. constant=', zeps(1)
    if (ync=='y')&
            write(nouti, *)'#coating diel. constant=', zeps(lcs)
    !--------/---------/---------/---------/---------/---------/---------/--
    if (ync=='y')&
            write(nouti, *)'#core radius/particle radius =', rff(1)
    write(nouti, *)'#phi,cos(theta),r(theta), and |e| in columns'
    write(nouti, *)

    open(unit = nouti + 1, file = 'elfcomp.dat')
    rewind(nouti + 1)
    write(nouti + 1, *)'#total electric field components around a tip'
    write(nouti + 1, *)
    write(nouti + 1, *)'#equiv.-volume-sphere radius in nm=', rev
    write(nouti + 1, *)
    write(nouti + 1, *)'#angular momentum cut-off lmax=', lmax
    write(nouti + 1, *)'#in the number nd=ndgs*lmax of gip points, ndgs=', ndgs
    write(nouti + 1, *)'#defp=', defp
    write(nouti + 1, *)'#material number =', nmat
    if (ync =='n') then
        write(nouti + 1, *)'#homogeneous particle'
    else if (ync=='y') then
        write(nouti + 1, *)'#coated particle'
    end if
    write(nouti + 1, *)
    write(nouti + 1, *)'#host dielectric constant=', zeps(lcs + 1)
    if (ync=='n')&
            write(nouti + 1, *)'#particle diel. constant=', zeps(1)
    if (ync=='y')&
            write(nouti + 1, *)'#coating diel. constant=', zeps(lcs)
    !--------/---------/---------/---------/---------/---------/---------/--
    if (ync=='y')&
            write(nouti + 1, *)'#core radius/particle radius =', rff(1)
    write(nouti + 1, *)'#phi,cos(theta),r(theta), e_r,e_theta,e_phi in columns'
    write(nouti + 1, *)
    !
    open(unit = nouti + 2, file = 'elfscat.dat')
    rewind(nouti + 2)
    write(nouti + 2, *)'#scattered e-field components around a tip'
    write(nouti + 2, *)
    write(nouti + 2, *)'#equiv.-volume-sphere radius in nm=', rev
    write(nouti + 2, *)
    write(nouti + 2, *)'#angular momentum cut-off lmax=', lmax
    write(nouti + 2, *)'#in the number nd=ndgs*lmax of gip points, ndgs=', ndgs
    write(nouti + 2, *)'#defp=', defp
    write(nouti + 2, *)'#material number =', nmat
    if (ync =='n') then
        write(nouti + 2, *)'#homogeneous particle'
    else if (ync=='y') then
        write(nouti + 2, *)'#coated particle'
    end if
    write(nouti + 2, *)
    write(nouti + 2, *)'#host dielectric constant=', zeps(lcs + 1)
    if (ync=='n')&
            write(nouti + 2, *)'#particle diel. constant=', zeps(1)
    if (ync=='y')&
            write(nouti + 2, *)'#coating diel. constant=', zeps(lcs)
    !--------/---------/---------/---------/---------/---------/---------/--
    if (ync=='y')&
            write(nouti + 2, *)'#core radius/particle radius =', rff(1)
    write(nouti + 2, *)'#phi,cos(theta),r(theta), e_r,e_theta,e_phi in columns'
    write(nouti + 2, *)
    !
    call evanesc(lmax, npp, lcs, revp, defpp, rsnm, hlength, lambda, thetv, &
            rmf, zeps)
    !
    close(nouti)
    close(nouti + 1)
    close(nouti + 2)

    write(6, *)'particle parameters:'
    write(6, *)
    write(6, *)'equivalent sphere radius =', rev
    if (ync=='n') write(6, *)'particle diel. constant=', zeps(1)
    write(6, *)'background dielectric constant zeps0=', zeps0
    if (ync=='y') write(6, *)'core diel. constant=', zeps(1)
    if (ync=='y') write(6, *)'coating diel. constant=', zeps(lcs)
    if (ync=='y') write(6, *)'core radius/particle radius =', &
            rff(1)
    write(6, *)
    write(6, *)'total |e|**2 in intnsty.dat'
    write(6, *)'total electric field components around a tip in elfcomp.dat'
    write(6, *)'scattered electric field components around a tip in elfscat.dat'

    1000 continue

    1005 format ('thet0=', f6.2, '  thet=', f6.2, '  phi0=', f6.2, &
            '  phi=', f6.2, '  alpha=', f6.2, '  beta=', f6.2)
    1006 format ('amplitude or s-matrix')
    5000 format ('4x4 phase matrix')

contains

    !**********************************************************************

    subroutine ampldr(nmax, np, eps, &
            rsnm, ht, lambda, zeps1, zeps0)

        ! warning in module ampldr in file ampldr.f: variables set but never used:
        !    nggg set at line 493 file ampldr.f
        !--------/---------/---------/---------/---------/---------/---------/--
        ! mpar%yncheck== 1 if you want to check gauss integrations
        ! convergence; otherwise mpar%yncheck=0
        ! nmax - angular momentum cut off
        ! lambda - the vacuum wavelength
        !
        ! outputs to common block t matrix
        !
        !                    |  tmt(m,m) |  tmt(m,e)   |
        !            tmt  =  | ----------+-------------|
        !                    |  tmt(e,m) |  tmt(e,e)   |
        !
        !    tmt(1,*) corresponds to tee scattering matrix
        !    tmt(2,*) corresponds to tmm scattering matrix
        !    tmt(3,*) corresponds to tme scattering matrix
        !    tmt(4,*) corresponds to tem scattering matrix
        !    tmt(4,*)=-tmt(3,*)^t where t denotes transposed tmt(3,*) submatrix
        !
        ! ichoice=1 if nag library is available, otherwise ichoice=2
        !
        ! np,eps: specifies the shape of particles within a given np class:
        !     np.gt.0 - eps = deformation parameter of a chebyshev particle
        !     np=-1 - eps = the ratio of the horizontal to rotational axes. eps is
        !             larger than 1 for oblate spheroids and smaller than 1 for
        !             prolate spheroids.
        !     np=-2 - eps = the ratio of the cylinder diameter to its length.
        !     np=-3 - no eps is specified
        !     np=-4 - eps is the height (along the axial symmetry axis)
        !             of the resulting cut sphere
        !                note that always eps.lt.2*rev specified
        !     np=-9 - eps = the ratio of the cylinder diameter to its length,
        !             nanorod caps are set independantly.
        !
        ! warning:
        !  in computations for spheres, use eps=1.000001 instead of eps=1.
        !  eps=1 can cause overflows in some rare cases.
        !
        !  lam - the wavelength of incident light in the ambient.
        !                  lam=lambda/sqrt(zeps0) here
        !
        !  rat = 1 - particle size is specified in terms of the
        !                equal-volume-sphere radius
        !  rat.ne.1 - particle size is specified in terms of the
        !                equal-surface-area-sphere radius
        !  axi ... equivalent-(volume/surface-area)-sphere radius
        !  rev=a=rat*axi ... equal-volume-sphere radius
        !                  (feeded as rev to rsp* routines)
        !  ddelt - required precision
        !  xs  - equiv. size parameter x=2*pi*rev*n_0/lambda
        !
        !  alpha and beta - euler angles (in degrees) specifying the
        !          orientation of the scattering particle relative to
        !          the laboratory reference frame (refs. 6 and 7).
        !
        !  thet0 - zenith angle of the incident beam in degrees
        !  thet - zenith angle of the scattered beam in degrees
        !  phi0 - azimuth angle of the incident beam in degrees
        !  phi - azimuth angle of the scattered beam in degrees
        !--------/---------/---------/---------/---------/---------/---------/--
        implicit none
        integer nout, naxsm
        integer nmax, np, inm1, ixxx, m, m1, n, n1, n11, n2, n22, ncheck, &
                ndgs, ngaus, ngauss, nm, nma, nn1, nn2, nnm, nnnggg
        ! number of the output unit
        parameter (nout = 35)

        real(dp) eps, a, alpha, beta, dext, dsca, dn1, &
                fac, p, phi, phi0, pir, pii, ppi, &
                qabs, qext, qext1, qsc, qsca, qsca1, &
                qxt, rat, rev, thet, thet0, &
                tr1nn, tr1nn1, ti1nn, ti1nn1, &
                walb, xev, z11, z12, z13, z14, &
                z21, z22, z23, z24, &
                z31, z32, z33, z34, &
                z41, z42, z43, z44, &
                zz1, zz2, zz3, zz4, &
                zz5, zz6, zz7, zz8
        real(dp)  lam, lambda, mrr, mri, rsnm, ht, ddelt, ddelta, &
                x(npng2), w(npng2), &
                s(npng2), ss(npng2), an(npn1), r(npng2), dr(npng2), &
                ddr(npng2), drr(npng2), dri(npng2), ann(npn1, npn1)
        real(dp) tr1(npn2, npn2), ti1(npn2, npn2)
        !      real(dp) xalpha(300),xbeta(300),walpha(300),wbeta(300)
        real(dp) rt11(npn6, npn4, npn4), rt12(npn6, npn4, npn4), &
                rt21(npn6, npn4, npn4), rt22(npn6, npn4, npn4), &
                it11(npn6, npn4, npn4), it12(npn6, npn4, npn4), &
                it21(npn6, npn4, npn4), it22(npn6, npn4, npn4)
        complex(dp) s11, s12, s21, s22
        complex(dp) zeps1, zeps0
        !
        common /ct/ tr1, ti1
        ! transfers the real and imaginary part of the t matrix (2*nmax,2*nmax)
        ! array for a given value of m from tmatr0 and tmatr to the ampldr
        !
        common /tmat/ rt11, rt12, rt21, rt22, it11, it12, it21, it22
        ! transfers t matrix arrays obtained from tr1,ti1 in the ampldr
        ! to the ampl routine
        !
        !     common /choice/ ichoice
        ! transfers the choice of inversion to relevant matrix inversion
        ! routines
        !
        common /toampld/rat, rev, alpha, beta, ddelt
        !
        ! transfers real(dp) rat,a(rev),alpha,beta,ddelt from the main here
        !
        common /totampld/thet0, thet, phi0, phi
        !
        ! transfers real(dp) thet0,thet,phi0,phi from the main here

        common /toiampld/ncheck, naxsm, ndgs
        ! transfers integers ncheck,naxsm,ndgs from the main here
        !
        !****************************************************************
        !
        mpar%eps = eps
        mpar%rev = rev
        p = dacos(-1d0)                   !local pi constant
        !
        a = rev
        lam = lambda/sqrt(zeps0)          !wavelength in the ambient

        !c      write(6,*)'lam,lambda in ampl=', lam, lambda
        !
        ! the real part of the refractive index contrast
        !
        mrr = dble(sqrt(zeps1/zeps0))
        !
        ! the imaginary  part of the refractive index contrast
        !
        mri = aimag(sqrt(zeps1/zeps0))
        !
        ddelta = 0.1d0*ddelt
        !
        ! ddelt is used to test the accuracy of computing the
        ! optical cross sections. this accuracy is usually better
        ! than the absolute accuracy of computing the expansion coefficients
        ! of a normalized scattering matrix by a factor of 10. therefore,
        ! the desired accuracy of computing the expansion coefficients
        ! is rescaled by a factor 0.1 before entering the test of the
        ! accuracy of computing the optical cross sections.

        if (dabs(rat - 1d0) > 1d-8) then
            select case(np)
            case(0:)
                call radii_ratio_chebyshev(np, eps, rat)
            case(-1)
                call radii_ratio_spheroid(eps, rat)
            case(-2)
                call radii_ratio_cylinder(eps, rat)
            case(-3)
                call radii_ratio_droplet (rat)
            case(-9)
                call radii_ratio_nanorod(eps, mpar%nanorod_cap_hr, rat)
            end select
        end if

        print 7400, lam, mrr, mri

        7400 format('lam=', f12.6, 3x, 'mrr=', d10.4, 3x, 'mri=', d10.4)
        !
        !___________________________________________________
        ! determination of the wiscombe value of the floating
        ! angular momentum cutoff nmax:

        xev = 2d0*p*a/lam
        ixxx = xev + 4.05d0*xev**0.333333d0     !wiscombe conv. criterion for nmax
        if (xev>1.) then
            inm1 = max0(3, ixxx)
        else
            inm1 = 2                 !the bessel package routine ryb
            !                        !requires nmax to be at least 2
        end if
        !
        if (inm1>=npn1) print 7333, npn1
        if (inm1>=npn1) stop
        7333 format('convergence is not obtained for npn1=', i3, &
                '.  execution terminated')
        !_______________________________________________________________

        ngauss = nmax*ndgs

        if (mpar%yncheck /= 0) then

            write(6, *)
            write(6, *)'nmax-convergence test'
            write(6, *)
            write(6, *)'(ngauss=nmax*ndgs)'


            ! internal determination of the floating angular momentum cutoff
            ! nmax using convergence criterion of {mis32}. it begins convergence
            ! convergence test with the wiscombe value for the floating angular
            ! momentum cutoff nmax with its subsequent increase by one, till
            ! the convergence criterion {mis32} is satisfied
            !
            qext1 = 0d0
            qsca1 = 0d0

            do nma = inm1, npn1
                nmax = nma
                ngauss = nmax*ndgs    !the number of the gauss integration points

                if (ngauss>npng1) print 7340, ngauss
                if (ngauss>npng1) stop

                7340      format('ngauss =', i3, ' i.e. is greater than npng1.', &
                        '  execution terminated')
                ! 7334    format(' nmax =', i3,'  dc2=',d8.2,'   dc1=',d8.2)
                !
                call const(ngauss, nmax, x, w, an, ann, s, ss, np, eps, rsnm, ht)      !in ampldr
                !
                ! specify particle shape:
                call vary(lam, mrr, mri, a, eps, &
                        rsnm, ht, np, ngauss, x, p, &
                        ppi, pir, pii, r, dr, ddr, drr, dri, nmax)
                !
                ! determine m=m'=0 elements of the t matrix
                !
                if (mpar%yn_adaptive == 0) then
                    call tmatr0 (ngauss, x, w, an, ann, ppi, pir, pii, r, dr, &
                            ddr, drr, dri, nmax, ncheck, naxsm)
                else
                    call tmatr0_adapt (ngauss, x, w, an, ann, ppi, pir, pii, r, dr, &
                            ddr, drr, dri, nmax, ncheck, naxsm)
                endif
                !
                qext = 0d0
                qsca = 0d0
                !
                ! make convergence test {mis32} for a given nmax:
                !
                do n = 1, nmax
                    n1 = n + nmax
                    tr1nn = tr1(n, n)
                    ti1nn = ti1(n, n)
                    tr1nn1 = tr1(n1, n1)
                    ti1nn1 = ti1(n1, n1)
                    dn1 = dble(2*n + 1)
                    qsca = qsca + dn1*(tr1nn*tr1nn + ti1nn*ti1nn &
                            + tr1nn1*tr1nn1 + ti1nn1*ti1nn1)
                    qext = qext + (tr1nn + tr1nn1)*dn1
                end do
                !>>> for debugging:
                !c      open(nout+1,file='tr1diag.dat')
                !c      open(nout+2,file='ti1diag.dat')
                !c              do n=1,2*nmax
                !c                  write(nout+1,*) tr1(n,n)
                !c                  write(nout+2,*) ti1(n,n)
                !c              enddo
                !c      close(nout+1)
                !c      close(nout+2)
                !<<<
                write(6, *)'nmax=', nmax
                write(6, *)'ngauss=', ngauss
                write(6, *)'qsca1=', qsca1
                write(6, *)'qsca=', qsca
                write(6, *)'qext1=', qext1
                write(6, *)'qext=', qext
                !<<<
                dsca = dabs((qsca1 - qsca)/qsca)
                dext = dabs((qext1 - qext)/qext)
                qext1 = qext
                qsca1 = qsca

                !        print 7334, nmax,dsca,dext

                if(dsca<=ddelta.and.dext<=ddelta) exit
                if (nma==npn1) print 7333, npn1
                if (nma==npn1) stop

                !successful l-convergence test exit&
            end do
            ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            !begin ngauss-convergence test

            write(6, *)
            write(6, *)'ngauss-convergence test'
            write(6, *)

            nnnggg = ngauss + 1

            if (ngauss==npng1) print 7336
            7336    format('warning: ngauss=npng1')

            if (ngauss/=npng1) then
                do ngaus = nnnggg, npng1
                    !
                    if (ngaus>300.and.ngaus<2595) cycle
                    if (ngaus==npng1) print 7336
                    !
                    ngauss = ngaus
                    !c         nggg=2*ngauss
                    !
                    ! gif division points and weights + other numerical constants
                    !
                    call const(ngauss, nmax, x, w, an, ann, s, ss, np, eps, rsnm, ht)     !in ampldr
                    !
                    ! specify particle shape:
                    !
                    call vary(lam, mrr, mri, a, eps, &
                            rsnm, ht, np, ngauss, x, p, &
                            ppi, pir, pii, r, dr, ddr, drr, dri, nmax)
                    !
                    ! determine m=m'=0 elements of the t matrix
                    !
                    if (mpar%yn_adaptive == 0) then
                        call tmatr0 (ngauss, x, w, an, ann, ppi, pir, pii, r, dr, &
                                ddr, drr, dri, nmax, ncheck, naxsm)
                    else
                        call tmatr0_adapt (ngauss, x, w, an, ann, ppi, pir, pii, r, dr, &
                                ddr, drr, dri, nmax, ncheck, naxsm)
                    endif
                    !
                    qext = 0d0
                    qsca = 0d0

                    do n = 1, nmax
                        n1 = n + nmax
                        tr1nn = tr1(n, n)
                        ti1nn = ti1(n, n)
                        tr1nn1 = tr1(n1, n1)
                        ti1nn1 = ti1(n1, n1)

                        dn1 = dble(2*n + 1)

                        qsca = qsca + dn1*(tr1nn*tr1nn + ti1nn*ti1nn&
                                + tr1nn1*tr1nn1 + ti1nn1*ti1nn1)
                        qext = qext + (tr1nn + tr1nn1)*dn1

                    end do

                    dsca = dabs((qsca1 - qsca)/qsca)
                    dext = dabs((qext1 - qext)/qext)

                    !        print 7337, nggg,dsca,dext
                    ! 7337    format(' ng=',i3,'  dc2=',d8.2,'   dc1=',d8.2)
                    !<<<
                    write(6, *)'ngauss=', ngauss
                    write(6, *)'qsca1=', qsca1
                    write(6, *)'qsca=', qsca
                    write(6, *)'qext1=', qext1
                    write(6, *)'qext=', qext

                    if(dsca<=ddelta.and.dext<=ddelta) exit
                    !<<<
                    qext1 = qext
                    qsca1 = qsca

                end do
                ! %%%%%%%%%%%%%%%%%% successful ngauss-convergence test %%%%%%%%%%%%%%%
            end if
        else if (mpar%yncheck == 0) then
                ! gif division points and weights + other numerical constants
                !
                call const(ngauss, nmax, x, w, an, ann, s, ss, np, eps, rsnm, ht)     !in ampldr
                !
                ! specify particle shape:
                !
                call vary(lam, mrr, mri, a, eps, &
                        rsnm, ht, np, ngauss, x, p, &
                        ppi, pir, pii, r, dr, ddr, drr, dri, nmax)
                !
                ! determine m=m'=0 elements of the t matrix
                !
                if (mpar%yn_adaptive == 0) then
                    call tmatr0 (ngauss, x, w, an, ann, ppi, pir, pii, r, dr, &
                            ddr, drr, dri, nmax, ncheck, naxsm)
                else
                    call tmatr0_adapt (ngauss, x, w, an, ann, ppi, pir, pii, r, dr, &
                            ddr, drr, dri, nmax, ncheck, naxsm)
                endif
                !
                !<<<
                !
        end if !yncheck

        write(6, *)
        write(6, *)'nmax=', nmax
        write(6, *)'ngauss=', ngauss
        write(6, *)
        !c      write(nout,*)'nmax=',nmax
        !c      write(nout,*)'ngauss=',ngauss
        !<<<
        !************   calculation of scattering cross sections   *********
        !initialization:

        qsca = 0d0
        qext = 0d0
        nnm = 2*nmax
        !   >>>  determination of qext and qsca contributions for m=0

        do n = 1, nnm

            qext = qext + tr1(n, n)

            !c         if ((n.le.5).or.(((n-nmax).le.5).and.((n-nmax).gt.0))) then
            !c         xx=-dble(tr1(n,n))       ! sin^2\eta_l
            !c         write(nout+15,*) 'n, sin^2\eta_l', n, xx
            !c         end if

        end do
        ! given rt1 and it1 matrices from tmatr0 routine,
        ! assigning of rt^{ij} and it^{ij} matrix entries to be
        ! used later by ampl routine

        do n2 = 1, nmax
            nn2 = n2 + nmax
            do n1 = 1, nmax
                nn1 = n1 + nmax
                zz1 = tr1(n1, n2)
                rt11(1, n1, n2) = zz1
                zz2 = ti1(n1, n2)
                it11(1, n1, n2) = zz2
                zz3 = tr1(n1, nn2)
                rt12(1, n1, n2) = zz3
                zz4 = ti1(n1, nn2)
                it12(1, n1, n2) = zz4
                zz5 = tr1(nn1, n2)
                rt21(1, n1, n2) = zz5
                zz6 = ti1(nn1, n2)
                it21(1, n1, n2) = zz6
                zz7 = tr1(nn1, nn2)
                rt22(1, n1, n2) = zz7
                zz8 = ti1(nn1, nn2)
                it22(1, n1, n2) = zz8
                !
                qsca = qsca + zz1*zz1 + zz2*zz2 + zz3*zz3 + zz4*zz4&
                        + zz5*zz5 + zz6*zz6 + zz7*zz7 + zz8*zz8
            end do
        end do !end of the loop over orbital numbers
        !________________
        !<<<

        if (abs(qsca)>(1.0001d0*abs(qext))) then
            write(6, *)'m=', 0
            write(6, *)'qsca=', qsca
            write(6, *)'qext=', qext
            write(6, *)&
                    'warning: abs(qsca).gt.abs(qext)!!!'
            !     stop
        end if

        if (qext>1.d-7) then
            write(6, *)&
                    'warning: partial sum qext has to be negative!'
            !      stop
        end if


        !c         write(nout,*)'m=',m
        !c         write(nout,*)'qsca=',qsca
        !c         write(nout,*)'qsc=',qsc
        !c         write(nout,*)'qext=',qext
        !c         write(nout,*)'qxt=',qxt
        !<<<
        !   >>>  determination of qext and qsca contributions for m >< 0

        do m = 1, nmax
            !
            !         call tmatr(m,ngauss,x,w,an,ann,s,ss,ppi,pir,pii,r,dr,
            !     &               ddr,drr,dri,nmax,ncheck,naxsm)
            if (mpar%yn_adaptive == 0) then
                call tmatr(m, ngauss, x, w, an, ann, s, ss, ppi, pir, pii, r, dr, &
                        ddr, drr, dri, nmax, ncheck, naxsm)
            else
                call tmatr_adapt(m, ngauss, x, w, an, ann, s, ss, ppi, pir, pii, r, dr, &
                        ddr, drr, dri, nmax, ncheck, naxsm)
            endif

!            call tmtr(m, ngauss, x, w, an, ann, ppi, pir, pii, r, dr, ddr, drr, &
!                    dri, nmax, ncheck, naxsm)
            !
            ! <<< returns  m=m'>0 elements of the t matrix
            !
            nm = nmax - m + 1
            m1 = m + 1
            qsc = 0d0
            ! given rt1 and it1 matrices from tmatr routine,
            ! assigning of rt^{ij} and it^{ij} matrix entries to be
            ! used later by ampl routine.
            !
            do n2 = 1, nm              !summation over orbital numbers
                ! conversion of the n22 index of rt1 and it1 matrices
                ! to the index nn2 of rt^{ij} and it^{ij} matrices

                nn2 = n2 + m - 1        !from m to nmax
                n22 = n2 + nm         !from nmax+1 to 2*nmax-m+1

                do n1 = 1, nm           !summation over orbital numbers
                    ! conversion of the n11 index of rt1 and it1 matrices
                    ! to the index nn1 of rt^{ij} and it^{ij} matrices

                    nn1 = n1 + m - 1        !from m to nmax
                    n11 = n1 + nm         !from nmax+1 to 2*nmax-m+1

                    zz1 = tr1(n1, n2)
                    rt11(m1, nn1, nn2) = zz1
                    zz2 = ti1(n1, n2)
                    it11(m1, nn1, nn2) = zz2
                    zz3 = tr1(n1, n22)
                    rt12(m1, nn1, nn2) = zz3
                    zz4 = ti1(n1, n22)
                    it12(m1, nn1, nn2) = zz4
                    zz5 = tr1(n11, n2)
                    rt21(m1, nn1, nn2) = zz5
                    zz6 = ti1(n11, n2)
                    it21(m1, nn1, nn2) = zz6
                    zz7 = tr1(n11, n22)
                    rt22(m1, nn1, nn2) = zz7
                    zz8 = ti1(n11, n22)
                    it22(m1, nn1, nn2) = zz8
                    !
                    qsc = qsc + (zz1*zz1 + zz2*zz2 + zz3*zz3 + zz4*zz4&
                            + zz5*zz5 + zz6*zz6 + zz7*zz7 + zz8*zz8)*2d0
                    !
                    ! multiplication by 2d0 here accounts for +/-m symmetry of resulting
                    ! expressions
                end do
            end do     !end of the loop over orbital numbers

            nnm = 2*nm
            qxt = 0d0

            do n = 1, nnm
                !multiplication by 2d0 accounts
                !for +/-m symmetry of resulting
                !expressions
                qxt = qxt + tr1(n, n)*2d0
            end do
            !<<<
            ! summation over magnetic quantum number:

            qsca = qsca + qsc
            qext = qext + qxt
            !<<<

            if (abs(qsc)>(1.0001d0*abs(qxt))) then
                write(6, *)'m=', m
                write(6, *)'qsca=', qsca
                write(6, *)'qsc=', qsc
                write(6, *)'qext=', qext
                write(6, *)'qxt=', qxt
                write(6, *)&
                        'warning: abs(qsc).gt.abs(qxt)!!!'
                !     stop
            end if

            if (qxt>=1d-7) then
                write(6, *)&
                        'warning: partial sum qxt has to be negative!'
                !      stop
            end if

            !write(nout,*)'m=',m
            !write(nout,*)'qsca=',qsca
            !write(nout,*)'qsc=',qsc
            !write(nout,*)'qext=',qext
            !write(nout,*)'qxt=',qxt
            !
            !print 7800,m,dabs(qxt),qsc,nmax
            !7800 format(' m=',i3,'  qxt=',d12.6,'  qsc=',d12.6, '  nmax=',i3)

        end do !end of loop over m's
        !<<<
        write(6, *)'qsca=', qsca
        write(6, *)'qext=', qext

        ! 'qsca' and '-qext' are now 'efficiency factors' for scattering
        ! and extinction (=\sum_{al} \sin^2\eta_{al}).

        qabs = -qext - qsca       !absorption
        walb = -qsca/qext       !albedo

        if (abs(walb)>1d0 + ddelta) then
            print 9111
            9111 format ('warning: the albedo walb is greater than 1')
            write(6, *)'walb=', walb
        end if
        !<<<
        ! in order to convert the efficiencies 'qsca' and '-qext' into
        ! normalized (per scatterer effective surface s=pi*rev**2)
        ! cross-sections
        !         qext=(2/x**2) \sum_{al} \sin^2\eta_{al}  for eta real
        !         qext=(2/x**2) \sum_{al} re (t)    for a general eta
        !         qsca=(2/x**2) \sum_{al} |t|^2
        ! (cf eq.(2.135-8) of newton�s book)
        ! at the moment, the prefactor (2/xev**2) is still missing.
        ! (lambda here is the wavelength in the exterior ambient medium)
        !c      write(6,*)'lam in ampl=', lam
        !         fac=lam**2/(2.d0*p**2*rev**2)     !=2/xs**2
        fac = 2.d0/xev**2
        write(nout, *)    lambda, fac*qsca
        write(nout + 1, *)  lambda, -fac*qext
        write(nout + 2, *)  lambda, fac*qabs
        write(nout + 5, *)  lambda, fac*walb
        write(nout + 10, *)
        write(nout + 10, *) lambda
        write(nout + 12, *)
        !write(nout+12,*) lambda
        !write(nout+13,*) lambda
        !write(nout+16,*) -qext
        !<<<
        !_________________________________________________________
        !  computation of the amplitude and phase matrices
        !  amplitude matrix [eqs. (2)-(4) of ref. 6]
        !
        call ampl (nmax, lam, thet0, thet, phi0, phi, alpha, beta, &
                s11, s12, s21, s22)
        !
        !  phase matrix [eqs. (13)-(29) of ref. 6]
        z11 = 0.5d0*(s11*conjg(s11) + s12*conjg(s12)&
                + s21*conjg(s21) + s22*conjg(s22))
        z12 = 0.5d0*(s11*conjg(s11) - s12*conjg(s12)&
                + s21*conjg(s21) - s22*conjg(s22))
        z13 = -s11*conjg(s12) - s22*conjg(s21)
        z14 = (0d0, 1d0)*(s11*conjg(s12) - s22*conjg(s21))
        z21 = 0.5d0*(s11*conjg(s11) + s12*conjg(s12)&
                - s21*conjg(s21) - s22*conjg(s22))
        z22 = 0.5d0*(s11*conjg(s11) - s12*conjg(s12)&
                - s21*conjg(s21) + s22*conjg(s22))
        z23 = -s11*conjg(s12) + s22*conjg(s21)
        z24 = (0d0, 1d0)*(s11*conjg(s12) + s22*conjg(s21))
        z31 = -s11*conjg(s21) - s22*conjg(s12)
        z32 = -s11*conjg(s21) + s22*conjg(s12)
        z33 = s11*conjg(s22) + s12*conjg(s21)
        z34 = (0d0, -1d0)*(s11*conjg(s22) + s21*conjg(s12))
        z41 = (0d0, 1d0)*(s21*conjg(s11) + s22*conjg(s12))
        z42 = (0d0, 1d0)*(s21*conjg(s11) - s22*conjg(s12))
        z43 = (0d0, -1d0)*(s22*conjg(s11) - s12*conjg(s21))
        z44 = s22*conjg(s11) - s12*conjg(s21)

        write(nout + 10, 5001) z11, z12, z13, z14
        write(nout + 10, 5001) z21, z22, z23, z24
        write(nout + 10, 5001) z31, z32, z33, z34
        write(nout + 10, 5001) z41, z42, z43, z44

        5001 format (4f10.4)

        !itime=mclock()
        !time=dble(itime)/6000d0
        !print 1001,time
        !1001 format (' time =',f8.2,' min')

        return
    end

    !********************************************************************

    subroutine ampl (nmax, dlam, tl, tl1, pl, pl1, alpha, beta, &
            vv, vh, hv, hh)
        !--------/---------/---------/---------/---------/---------/---------/--
        ! >>> nmax,dlam,tl,tl1,pl,pl1,alpha,beta
        ! <<< vv,vh,hv,hh
        !=================
        !
        !    given t matrix in common block it calculates the amplitude matrix
        !
        !    this routine closely follows exposition by
        !       m. i. mishchenko, calculation of the amplitude matrix
        !       for a nonspherical particle in a fixed orientation,
        !       appl. opt. vol. 39, 1026-1031 (2000).
        !
        !   nmax - angular momentum cutoff
        !   dlam=lambda/sqrt(zeps0)  - wavelength of incident light in the ambient
        !                     (vacuum wavelength divided by sqrt(zeps0))
        !
        !   lambda - vacuum wavelength. determined as dlam*sqrt(zeps0) and
        !            only used here for the write out purposes
        !   tl,tl1,pl,pl1 ... angles in degrees
        !                     determined w.r.t laboratory frame:
        !   tl (thet0 in main) - zenith angle of the incident beam in degrees
        !   tl1 (thet in main) - zenith angle of the scattered beam in degrees
        !   pl (phi0 in main) - azimuth angle of the incident beam in degrees
        !   pl1 (phi in main) - azimuth angle of the scattered beam in degrees
        !
        !   alpha and beta - euler angles (in degrees) specifying the
        !         orientation of the scattering particle relative to the
        !         laboratory reference frame (refs. 6 and 7).
        !   vv,vh,hv,hh ... amplitude scattering matrix elements s11,s12,s21,s22
        !--------/---------/---------/---------/---------/---------/---------/--
        !implicit real(dp) (a-b, d-h, o-z)
        !implicit complex(dp) (c)
        implicit none
        integer nout
        integer nmax, i, j, k, m, m1, n, nmin, nn
        ! number of the output unit
        parameter (nout = 35)
        !include 'ampld.par.f'
        real(dp) rev, eps, hlength, x, y
        real(dp) tl, tl1, pl, pl1, alpha, beta, &
                alph, bet, d, d11, d12, d21, d22, &
                dcth, dcth0, dk, dnn, dv1n, dv1nn, &
                dv2n, dv2nn, fac, fc, fs, ph, phil, &
                phil1, phip, phip1, pigrad, pin, pin2, &
                rn, rsnm, sa, sb, sp, sp1, spp, spp1, &
                st, st1, thetl, thetl1, thetp, thetp1
        real(dp) dlam, lambda, cext, cext1, cext2
        real(dp) al(3, 2), al1(3, 2), ap(2, 3), ap1(2, 3), b(3, 3), &
                r(2, 2), r1(2, 2), c(3, 2), ca, cb, ct, cp, ctp, cpp, ct1, cp1, &
                ctp1, cpp1
        real(dp) dv1(npn6), dv2(npn6), dv01(npn6), dv02(npn6)
        real(dp) tr11(npn6, npn4, npn4), tr12(npn6, npn4, npn4), &
                tr21(npn6, npn4, npn4), tr22(npn6, npn4, npn4), &
                ti11(npn6, npn4, npn4), ti12(npn6, npn4, npn4), &
                ti21(npn6, npn4, npn4), ti22(npn6, npn4, npn4)

        complex(dp) chh, chv, cn, cn1, cn2, &
                ct11, ct12, ct21, ct22, &
                cvh, cvv
        complex(dp) cal(npn4, npn4), vv, vh, hv, hh, zeps0
        !_____
        common /tmat/ tr11, tr12, tr21, tr22, ti11, ti12, ti21, ti22
        !_
        common /dielf/ zeps0
        common /revf/ rev
        common /cylpar/ rsnm, hlength

        ! transfers zeps0,rev here from the main

        ! checking the initial set of angles tl,tl1,pl,pl1,alpha,beta
        ! for allowability

        if (alpha<0d0.or.alpha>360d0.or.&
                beta<0d0.or.beta>180d0.or.&
                tl<0d0.or.tl>180d0.or.&
                tl1<0d0.or.tl1>180d0.or.&
                pl<0d0.or.pl>360d0.or.&
                pl1<0d0.or.pl1>360d0) then
            write(nout, 2000)
            stop
        else

        endif
        2000 format ('an angular parameter is outside its', &
                ' allowable range')
        ! specifying numerical constants:

        pin = dacos(-1d0)         !=pi
        pin2 = pin*0.5d0        !=pi/2
        pigrad = pin/180d0      !=pi/180
        ! conversion from degrees to radians:
        alph = alpha*pigrad
        bet = beta*pigrad
        thetl = tl*pigrad
        phil = pl*pigrad
        thetl1 = tl1*pigrad
        phil1 = pl1*pigrad
        ! initialization of the vacuum wavelength lambda

        lambda = dlam*sqrt(zeps0)         !vacuum wavelength

        eps = 1d-7
        if (thetl<pin2) thetl = thetl + eps
        if (thetl>pin2) thetl = thetl - eps
        if (thetl1<pin2) thetl1 = thetl1 + eps
        if (thetl1>pin2) thetl1 = thetl1 - eps
        if (phil<pin) phil = phil + eps
        if (phil>pin) phil = phil - eps
        if (phil1<pin) phil1 = phil1 + eps
        if (phil1>pin) phil1 = phil1 - eps
        if (bet<=pin2.and.pin2 - bet<=eps) bet = bet - eps
        if (bet>pin2.and.bet - pin2<=eps) bet = bet + eps

        !   given tl,tl1,pl,pl1 in laboratory frame
        !   compute thetp, phip, thetp1, and phip1 in particle frame
        !   (see eqs. (9), (20), and (21))
        !
        ! incident beam:

        cb = dcos(bet)
        sb = dsin(bet)
        ct = dcos(thetl)
        st = dsin(thetl)
        cp = dcos(phil - alph)
        sp = dsin(phil - alph)

        ctp = ct*cb + st*sb*cp       !eq. (9)
        thetp = dacos(ctp)
        cpp = cb*st*cp - sb*ct       !eq. (20)
        spp = st*sp                      !eq. (21)
        phip = datan(spp/cpp)

        if (phip>0d0.and.sp<0d0) phip = phip + pin
        if (phip<0d0.and.sp>0d0) phip = phip + pin
        if (phip<0d0) phip = phip + 2d0*pin

        ! scattered beam:

        ct1 = dcos(thetl1)
        st1 = dsin(thetl1)
        cp1 = dcos(phil1 - alph)
        sp1 = dsin(phil1 - alph)

        ctp1 = ct1*cb + st1*sb*cp1    !eq. (9)
        thetp1 = dacos(ctp1)
        cpp1 = cb*st1*cp1 - sb*ct1    !eq. (20)
        spp1 = st1*sp1                    !eq. (21)
        phip1 = datan(spp1/cpp1)

        if (phip1>0d0.and.sp1<0d0) phip1 = phip1 + pin
        if (phip1<0d0.and.sp1>0d0) phip1 = phip1 + pin
        if (phip1<0d0) phip1 = phip1 + 2d0*pin

        !____________compute matrix beta, eq. (22) of {mis39}

        ca = dcos(alph)
        sa = dsin(alph)
        b(1, 1) = ca*cb
        b(1, 2) = sa*cb
        b(1, 3) = -sb
        b(2, 1) = -sa
        b(2, 2) = ca
        b(2, 3) = 0d0
        b(3, 1) = ca*sb
        b(3, 2) = sa*sb
        b(3, 3) = cb

        !____________compute 3x2 matrices al and al1 for incident and
        !              scattered beams in laboratory frame
        !                   [eq. (15)  of {mis39}]

        cp = dcos(phil)
        sp = dsin(phil)
        cp1 = dcos(phil1)
        sp1 = dsin(phil1)

        ! incident beam:
        al(1, 1) = ct*cp
        al(1, 2) = -sp
        al(2, 1) = ct*sp
        al(2, 2) = cp
        al(3, 1) = -st
        al(3, 2) = 0d0

        ! scattered beam:
        al1(1, 1) = ct1*cp1
        al1(1, 2) = -sp1
        al1(2, 1) = ct1*sp1
        al1(2, 2) = cp1
        al1(3, 1) = -st1
        al1(3, 2) = 0d0

        !____________compute 2x3 matrices ap^(-1) and ap1^(-1) for incident
        !             and scattered beams in particle frame
        !                   [eq. (16)  of {mis39}]
        ct = ctp
        st = dsin(thetp)
        cp = dcos(phip)
        sp = dsin(phip)
        ct1 = ctp1
        st1 = dsin(thetp1)
        cp1 = dcos(phip1)
        sp1 = dsin(phip1)
        ! incident beam:
        ap(1, 1) = ct*cp
        ap(1, 2) = ct*sp
        ap(1, 3) = -st
        ap(2, 1) = -sp
        ap(2, 2) = cp
        ap(2, 3) = 0d0
        ! scattered beam:
        ap1(1, 1) = ct1*cp1
        ap1(1, 2) = ct1*sp1
        ap1(1, 3) = -st1
        ap1(2, 1) = -sp1
        ap1(2, 2) = cp1
        ap1(2, 3) = 0d0

        !____________compute matrices r and r^(-1), eq. (14)

        ! computation of r for incident beam
        do i = 1, 3
            do j = 1, 2
                x = 0d0
                do k = 1, 3
                    x = x + b(i, k)*al(k, j)
                enddo
                c(i, j) = x
            enddo
        enddo
        do i = 1, 2
            do j = 1, 2
                x = 0d0
                do k = 1, 3
                    x = x + ap(i, k)*c(k, j)
                enddo
                r(i, j) = x
            enddo
        enddo

        ! computation of r^(-1) for scattered beam using cramers rule:

        do i = 1, 3
            do j = 1, 2
                x = 0d0
                do k = 1, 3
                    x = x + b(i, k)*al1(k, j)
                enddo
                c(i, j) = x
            enddo
        enddo
        do i = 1, 2
            do j = 1, 2
                x = 0d0
                do k = 1, 3
                    x = x + ap1(i, k)*c(k, j)
                enddo
                r1(i, j) = x
            enddo
        enddo
        !====
        ! r for scattered beam determined, now cramers rule:

        d = 1d0/(r1(1, 1)*r1(2, 2) - r1(1, 2)*r1(2, 1))
        x = r1(1, 1)
        r1(1, 1) = r1(2, 2)*d
        r1(1, 2) = -r1(1, 2)*d
        r1(2, 1) = -r1(2, 1)*d
        r1(2, 2) = x*d

        !____________matrices r and r^(-1) determined
        !=========================================================
        !      the amplitude matrix in particle frame
        !
        !     ci=(0d0,1d0)

        ! >>> alpha numerical prefactors without phi-angles
        !     (following eq. (28))

        do nn = 1, nmax
            do n = 1, nmax
                cn = ci**(nn - n - 1)
                dnn = dble((2*n + 1)*(2*nn + 1))
                dnn = dnn/dble(n*nn*(n + 1)*(nn + 1))
                rn = dsqrt(dnn)
                cal(n, nn) = cn*rn
            end do
        end do

        dcth0 = ctp             !\cos\vartheta_{inc}^p
        dcth = ctp1             !\cos\vartheta_{sca}^p
        ph = phip1 - phip       !(\varphi_{sca}^p-\varphi_{inc}^p)
        ! amplitude scattering matrix elements s11,s12,s21,s22 initialization

        vv = (0d0, 0d0)
        vh = (0d0, 0d0)
        hv = (0d0, 0d0)
        hh = (0d0, 0d0)
        !______________________________________________________________
        ! main summation loop:

        do m = 0, nmax
            m1 = m + 1
            nmin = max(m, 1)
            !
            ! specify pi- and tau- scattering functions:

            call vigampl (dcth, nmax, m, dv1, dv2)
            call vigampl (dcth0, nmax, m, dv01, dv02)
            !
            fc = 2d0*dcos(m*ph)    !takes into account +/- m contribution
            fs = 2d0*dsin(m*ph)
            !
            do nn = nmin, nmax

                dv1nn = dv01(nn)           !\pi-functions
                dv2nn = dv02(nn)           !\tau-functions

                do  n = nmin, nmax
                    dv1n = dv1(n)           !\pi-functions
                    dv2n = dv2(n)           !\tau-functions

                    ct11 = cmplx_dp(tr11(m1, n, nn), ti11(m1, n, nn))
                    ct22 = cmplx_dp(tr22(m1, n, nn), ti22(m1, n, nn))

                    if (m==0) then     !t^{21}=t^{12}=0 in particle frame

                        cn = cal(n, nn)*dv2n*dv2nn

                        vv = vv + cn*ct22
                        hh = hh + cn*ct11

                    else   !t^{21}\neq t^{12}\neq 0

                        ct12 = cmplx_dp(tr12(m1, n, nn), ti12(m1, n, nn))
                        ct21 = cmplx_dp(tr21(m1, n, nn), ti21(m1, n, nn))
                        ! complete \alpha-factors (eq. (28)) taking
                        ! into account w.r.t. summation over +/- m in particle frame:
                        !
                        !     t^{11}_{-mnn'} = t^{11}_{mnn'}; t^{22}_{-mnn'} = t^{22}_{mnn'}
                        !  t^{12}_{-mnn'} = - t^{12}_{mnn'}; t^{21}_{-mnn'} = - t^{21}_{mnn'}

                        cn1 = cal(n, nn)*fc
                        cn2 = cal(n, nn)*fs

                        d11 = dv1n*dv1nn    !\pi-\pi
                        d12 = dv1n*dv2nn    !\pi-\tau
                        d21 = dv2n*dv1nn    !\tau-\pi
                        d22 = dv2n*dv2nn    !\tau-\tau

                        vv = vv + (ct11*d11 + ct21*d21&
                                + ct12*d12 + ct22*d22)*cn1

                        vh = vh + (ct11*d12 + ct21*d22&
                                + ct12*d11 + ct22*d21)*cn2

                        hv = hv - (ct11*d21 + ct21*d11&
                                + ct12*d22 + ct22*d12)*cn2

                        hh = hh + (ct11*d22 + ct21*d12&
                                + ct12*d21 + ct22*d11)*cn1
                    endif
                end do
            end do     !(over n,n')
        end do        !end of main summation loop (over m)

        ! final multiplication of s11,s12,s21,s22 by (1/k) in the
        ! original code:

        dk = 2d0*pin/dlam     !wavevector in surrounding medium
        vv = vv/dk
        vh = vh/dk
        hv = hv/dk
        hh = hh/dk

        !   amplitude scattering matrix elements s11,s12,s21,s22 determined
        !==================================================================
        ! transformation of the amplitude matrix from particle to
        ! laboratory frame:

        cvv = vv*r(1, 1) + vh*r(2, 1)
        cvh = vv*r(1, 2) + vh*r(2, 2)
        chv = hv*r(1, 1) + hh*r(2, 1)
        chh = hv*r(1, 2) + hh*r(2, 2)
        vv = r1(1, 1)*cvv + r1(1, 2)*chv
        vh = r1(1, 1)*cvh + r1(1, 2)*chh
        hv = r1(2, 1)*cvv + r1(2, 2)*chv
        hh = r1(2, 1)*cvh + r1(2, 2)*chh

        print 1101, vv
        print 1102, vh
        print 1103, hv
        print 1104, hh
        ! for particles with plane of symmetry:
        ! if thet0=thet=90, then the incidence is perpendicular to
        ! the axis of axial symmetry  ===> e_\theta is along the axis
        ! of axial symmetry, whereas e_\phi is perpendicular to
        ! the symmetry axis. then
        !      c_{ext} = (4\pi/k_1} \mb{im} s_{11}    if e_\phi=0
        !
        ! or
        !
        !      c_{ext} = (4\pi/k_1} \mb{im} s_{22}    if e_\theta=0
        !
        ! for particles with plane of symmetry:
        ! extiction for e along the axis of axial symmetry:
        cext1 = 2.d0*dlam*aimag(vv)      !eq. (2.159)

        ! extiction for e perpendicular to the axis of axial symmetry:
        cext2 = 2.d0*dlam*aimag(hh)      !eq. (2.159)

        ! orientationally averaged extiction
        cext = dlam*aimag(vv + hh)         !eq. (5.97)

        write(6, *)'c_{ext}=\fr{2\pi}{k_1} im (s_{11}+s_{22})=', &
                cext               !=2.d0*pin*aimag(vv+hh)/k_1

        write(10, *) lambda, cext1, cext2
        !c      fac=lambda**2/(2.d0*pin**2*rev**2)     !=2/xs**2
        fac = 1.d0/(pin*rev**2)   !an effective geom. cross section

        write(nout + 3, 1107)lambda, fac*cext, fac*cext1, fac*cext2
        write(nout + 12, 1105) lambda, vv, vh
        write(nout + 12, 1105) lambda, hv, hh
        write(nout + 13, 1106) lambda, dble(vv*conjg(vv) + vh*conjg(vh)), &
                dble((vv + vh)*conjg(vv + vh))
        write(nout + 14, 1106) lambda, dble(hv*conjg(hv) + hh*conjg(hh)), &
                dble((hv + hh)*conjg(hv + hh))

        1101 format ('s11=', d11.5, ' + i*', d11.5)
        1102 format ('s12=', d11.5, ' + i*', d11.5)
        1103 format ('s21=', d11.5, ' + i*', d11.5)
        1104 format ('s22=', d11.5, ' + i*', d11.5)
        1105 format (f8.2, 5x, d11.5, 2x, d11.5, 5x, d11.5, 2x, d11.5)
        1106 format (f8.2, 5x, d11.5, 5x, d11.5)
        1107 format (f10.4, 3(5x, d16.8))

        return
    end

    !********************************************************************
    subroutine tmtr (m, ngauss, x, w, an, ann, ppi, pir, pii, r, dr, ddr, &
            drr, dri, nmax, ncheck, naxsm)
        !--------/---------/---------/---------/---------/---------/---------/--
        ! >>> m,ngauss,x,w,an,ann,s,ss,ppi,pir,pii,r,dr,ddr,drr,dri,nmax,ncheck
        ! <<< common blocks /tmat99/, /ct/ (for main),  and /ctt/ (for tt)
        !=====================
        !
        !  determines the t-matrix of an axially symmetric scatterer
        !                           for m.gt.0
        !
        !  m      - azimuthal number
        !  ngauss - the number of gif division points
        !  x=\cos\theta  - gif division points
        !  w - gif weights
        !  an(n)=n*(n+1)
        !  ann(l_1,l_2)=\sqrt{\fr{(2 l_1+1)}{l_1(l_1+1)} }
        !                       \sqrt{\fr{(2 l_2+1)}{l_2(l_2+1)} }/2
        !  nmax - angular momentum cutoff
        !  ncheck  -  .eq.0  then  ngss=2*ngauss, factor=1d0
        !             .eq.1  then  ngss = ngauss, factor=2d0
        !  naxsm   -  .eq.0 : gauss abscissas do not have +/- theta symmetry
        !             .eq.1 : gauss abscissas have +/- theta symmetry
        !  ncheck - specifies whether ng=2*ngauss or otherwise
        !  p=dacos(-1d0)
        !  pi=p*2d0/lam - wave vector
        !  ppi=pi*pi
        !  pir=ppi*mrr
        !  pii=ppi*mri
        !  r=r^2(\theta)                       for axially symmetric particles
        !  dr=[dr(\theta)/(d\theta)]/r(\theta) for axially symmetric particles
        !  ddr=\lambda/[2*\pi*r(\theta)]=1/(k_out*r)
        !  drr=(mrr/(mrr**2+mri**2))*(\lambda/[2*\pi*r(\theta)])
        !                  = re 1/(k_in*r)
        !  dri=-(mri/(mrr**2+mri**2))*(\lambda/[2*\pi*r(\theta)])
        !                  = im 1/(k_in*r)
        !
        !  refractive index outside is assumed to real, whereas inside
        !  a scatterer, refractive index is allowed to be complex in general.
        !  consequently, the bessel function j_l(k_in*r) will in general
        !  be complex. the routine below performs waterman surface integral
        !  separately for the real and imaginary parts of the integrand.
        !
        !--------/---------/---------/---------/---------/---------/---------/--
        !include 'ampld.par.f'
        implicit none
        integer m, ngauss, nmax, ncheck, naxsm, i, i1, i2, k1, k2, kk1, &
                kk2, mm1, n, n1, n2, ng, ngss, nm, nnmax
        real(dp) ppi, pir, pii, a11, a12, a21, a22, &
                aa1, aa2, ar11, ar12, ar21, ar22, &
                ai11, ai12, ai21, ai22, an1, an2, &
                an12, b1r, b1i, b2r, b2i, b3r, b3i, &
                b4r, b4i, b5r, b5i, b6r, b6i, &
                b7r, b7i, b8r, b8i, &
                c1r, c1i, c2r, c2i, &
                c3r, c3i, c4r, c4i, &
                c5r, c5i, c6r, c6i, &
                c7r, c7i, c8r, c8i, &
                d1n1, d1n2, d2n1, d2n2, &
                d3n1, d3n2, ddri, drii, &
                drri, e1, e2, e3, f1, f2, &
                factor, gr11, gr12, gr21, gr22, &
                gi11, gi12, gi21, gi22, &
                qdj1, qdjr2, qdji2, &
                qdy1, qj1, qjr2, qji2, qm, qmm, &
                qy1, rri, si, &
                tar11, tar12, tar21, tar22, &
                tai11, tai12, tai21, tai22, &
                tgr11, tgr12, tgr21, tgr22, &
                tgi11, tgi12, tgi21, tgi22, &
                tpir, tpii, tppi, uri, wr
        real(dp)  x(npng2), w(npng2), an(npn1), &
                r(npng2), dr(npng2), sig(npn2), &
                ddr(npng2), drr(npng2), &
                d1(npng2, npn1), d2(npng2, npn1), d3(npng2, npn1), &
                dri(npng2), rr(npng2), &
                dv1(npn3), ddv1(npn3), dv2(npn3), &
                dd1, dd2, dd3

        real(dp)  r11(npn1, npn1), r12(npn1, npn1), &
                r21(npn1, npn1), r22(npn1, npn1), &
                i11(npn1, npn1), i12(npn1, npn1), &
                i21(npn1, npn1), i22(npn1, npn1), &
                rg11(npn1, npn1), rg12(npn1, npn1), &
                rg21(npn1, npn1), rg22(npn1, npn1), &
                ig11(npn1, npn1), ig12(npn1, npn1), &
                ig21(npn1, npn1), ig22(npn1, npn1), &
                ann(npn1, npn1), &
                qr(npn2, npn2), qi(npn2, npn2), &
                rgqr(npn2, npn2), rgqi(npn2, npn2), &
                tqr(npn2, npn2), tqi(npn2, npn2), &
                trgqr(npn2, npn2), trgqi(npn2, npn2)
        !real(dp) tr1(npn2,npn2),ti1(npn2,npn2)
        !________
        common /tmat99/&
                r11, r12, r21, r22, i11, i12, i21, i22, rg11, rg12, rg21, rg22, &
                ig11, ig12, ig21, ig22          !only between tmatr routines
        !common /ct/ tr1,ti1                    !output from tt routine
        common /ctt/ qr, qi, rgqr, rgqi         !input for tt routine
        !________
        mm1 = m
        qm = dble(m)
        qmm = qm*qm
        ng = 2*ngauss
        nm = nmax + nmax
        factor = 1d0
        !
        if (ncheck==1) then          !theta=pi/2 is scatterer mirror symmetry plane
            ngss = ngauss
            factor = 2d0
        else if (ncheck==0) then     !theta=pi/2 is not a scatterer mirror symmetry plane
            ngss = ng
        endif
        !
        si = 1d0
        do n = 1, nm                 !nm=2*nmax
            si = -si
            sig(n) = si              !=(-1)**n
        end do
        !
        ! assigning wigner d-matrices:

        do i = 1, ngauss

            i1 = ngauss - i + 1
            i2 = ngauss + i
            !
            call vigf(x(i1), nmax, m, dv1, dv2, ddv1)
            !--------/---------/---------/---------/---------/---------/---------/--
            ! >>> x,nmax,m (only nonnegative)
            ! <<< dv1,dv2,ddv1
            ! =============
            !
            !     x=cos(theta), where theta is the polar angle
            !     nmax ... floating  angular momentum cutoff
            !
            ! returns \pi and \tau scattering functions in terms of the
            ! wigner d-functions. algorithm as described in eqs. (31-35)
            !  of ref. \cite{mis39} used. (note however a missing $n$
            !  factor in the 2nd term in the curly bracket in
            !   eq. (35) of ref. \cite{mis39}.)
            !
            !  for a given azimuthal number m.ge.0 returns
            !  the wigner d-functions
            !            dv1(n)=dvig(0,m,n,arccos x) = d_{0m}^{(l)}
            !
            !  \pi scattering function:
            !     ddv1(n)=dvig(0,m,n,arccos x)/sin(arccos x)
            !                              = m*d_{0m}^{(l)}/ sin\theta
            !
            !  \tau scattering function:
            !     dv2(n)=[d/d(arccos x)] dvig(0,m,n,arccos x)
            !                              = d d_{0m}^{(l)}/d\theta
            !--------/---------/---------/---------/---------/---------/---------/--

            do n = 1, nmax

                dd1 = ddv1(n)
                dd2 = dv2(n)
                dd3 = dv1(n)
                d1(i1, n) = dd1
                d2(i1, n) = dd2
                d3(i1, n) = dd3

                if (naxsm==1) then         !gauss abscissas chosen +/- symmetric

                    ! using (4.2.4) and (4.2.6) of {ed},
                    !           d_{0m}^{(l)}(\pi-\theta) = (-1)^{l+m} d_{0m}^{(l)}(\theta)

                    si = sig(n + m)           !=(-1)**(n+m)
                    !                                 !exactly what follows from {ed}
                    d1(i2, n) = dd1*si
                    d2(i2, n) = -dd2*si
                    d3(i2, n) = dd3*si

                end if
            enddo

            if (naxsm==0) then        !gauss abscissas not chosen +/- symmetric

                call vigf(x(i2), nmax, m, dv1, dv2, ddv1)

                do n = 1, nmax
                    dd1 = ddv1(n)
                    dd2 = dv2(n)
                    dd3 = dv1(n)
                    d1(i2, n) = dd1
                    d2(i2, n) = dd2
                    d3(i2, n) = dd3
                enddo

            end if
        end do

        !  assigning r^2(\theta)*weight product:

        do i = 1, ngss
            wr = w(i)*r(i)

            !c          if (dr(i).eq.0.d0) wr=0.d0   !temporarily only

            rr(i) = wr            !w(i)*r^2(\theta)
        end do

        do n1 = mm1, nmax         !mm1=m below
            an1 = an(n1)

            do n2 = mm1, nmax
                an2 = an(n2)
                ar11 = 0d0
                ar12 = 0d0
                ar21 = 0d0
                ar22 = 0d0
                ai11 = 0d0
                ai12 = 0d0
                ai21 = 0d0
                ai22 = 0d0
                gr11 = 0d0
                gr12 = 0d0
                gr21 = 0d0
                gr22 = 0d0
                gi11 = 0d0
                gi12 = 0d0
                gi21 = 0d0
                gi22 = 0d0
                si = sig(n1 + n2)

                do i = 1, ngss    !=ngauss   if ncheck.eq.1
                    !                                  !=2*ngauss if ncheck.eq.0
                    d1n1 = d1(i, n1)
                    d2n1 = d2(i, n1)
                    d3n1 = d3(i, n1)
                    d1n2 = d1(i, n2)
                    d2n2 = d2(i, n2)
                    d3n2 = d3(i, n2)
                    a11 = d1n1*d3n2            !pi(n1)*d(n2)
                    a12 = d3n1*d2n2            !d(n1)*tau(n2)
                    a21 = d2n1*d3n2            !tau(n1)*d(n2)
                    a22 = d2n1*d2n2            !tau(n1)*tau(n2)
                    aa1 = d1n1*d2n2 + d2n1*d1n2  !pi(n1)*tau(n2)+tau(n1)*pi(n2)
                    aa2 = d1n1*d1n2 + a22       !pi(n1)*pi(n2)+tau(n1)*tau(n2)
                    ! vector spherical harmonics:
                    !  since refractive index is allowed to be complex in general,
                    !  the bessel function j_l(k_in*r) is complex. the code below
                    !  performs a separation of the complex integrand in waterman's
                    !  surface integral into its respective real and imaginary
                    !  parts:
                    ! bessel functions of the exterior argument:

                    qj1 = cbess%j(i, n1)
                    qy1 = cbess%y(i, n1)
                    qdj1 = cbess%dj(i, n1)
                    qdy1 = cbess%dy(i, n1)
                    ! bessel functions of the interior argument:

                    qjr2 = cbess%jr(i, n2)
                    qji2 = cbess%ji(i, n2)
                    qdjr2 = cbess%djr(i, n2)
                    qdji2 = cbess%dji(i, n2)
                    ! re and im of j_{n2}(k_{in}r) j_{n1}(k_{out}r):

                    c1r = qjr2*qj1
                    c1i = qji2*qj1
                    ! re and im of j_{n2}(k_{in}r) h_{n1}(k_{out}r):

                    b1r = c1r - qji2*qy1
                    b1i = c1i + qjr2*qy1
                    ! re and im of j_{n2}(k_{in}r) j_{n1}'(k_{out}r)/(k_{out}r):

                    c2r = qjr2*qdj1
                    c2i = qji2*qdj1
                    ! re and im of j_{n2}(k_{in}r) h_{n1}'(k_{out}r)/(k_{out}r):

                    b2r = c2r - qji2*qdy1
                    b2i = c2i + qjr2*qdy1

                    ddri = ddr(i)               !1/(k_{out}r)
                    ! re and im of [1/(k_{out}r)]*j_{n2}(k_{in}r) j_{n1}(k_{out}r)

                    c3r = ddri*c1r
                    c3i = ddri*c1i
                    ! re and im of [1/(k_{out}r)]*j_{n2}(k_{in}r) h_{n1}(k_{out}r):

                    b3r = ddri*b1r
                    b3i = ddri*b1i
                    ! re and im of [k_{in}r j_{n2}(k_{in}r)]'/(k_{in}r)
                    !                         *j_{n1}(k_{out}r):

                    c4r = qdjr2*qj1
                    c4i = qdji2*qj1
                    ! re and im of [k_{in}r j_{n2}(k_{in}r)]'/(k_{in}r)
                    !                         * h_{n1}(k_{out}r):

                    b4r = c4r - qdji2*qy1
                    b4i = c4i + qdjr2*qy1

                    drri = drr(i)               !re[1/(k_{in}r)]
                    drii = dri(i)               !im[1/(k_{in}r)]
                    ! re and im of [1/(k_{in}r)] j_{n2}(k_{in}r) j_{n1}(k_{out}r):

                    c5r = c1r*drri - c1i*drii
                    c5i = c1i*drri + c1r*drii
                    ! re and im of [1/(k_{in}r)] j_{n2}(k_{in}r) h_{n1}(k_{out}r):

                    b5r = b1r*drri - b1i*drii
                    b5i = b1i*drri + b1r*drii
                    ! re and im of j_{n2}'(k_{in}r) j_{n1}'(k_{out}r):

                    c6r = qdjr2*qdj1
                    c6i = qdji2*qdj1
                    ! re and im of j_{n2}'(k_{in}r) h_{n1}'(k_{out}r):

                    b6r = c6r - qdji2*qdy1
                    b6i = c6i + qdjr2*qdy1
                    ! re and im of [1/(k_{out}r)] j_{n2}'(k_{in}r) j_{n1}(k_{out}r):

                    c7r = c4r*ddri
                    c7i = c4i*ddri
                    ! re and im of [1/(k_{out}r)] j_{n2}'(k_{in}r) h_{n1}(k_{out}r):

                    b7r = b4r*ddri
                    b7i = b4i*ddri
                    ! re and im of [1/(k_{in}r)] j_{n2}(k_{in}r) j_{n1}'(k_{out}r):

                    c8r = c2r*drri - c2i*drii
                    c8i = c2i*drri + c2r*drii
                    ! re and im of [1/(k_{in}r)] j_{n2}(k_{in}r) h_{n1}'(k_{out}r):

                    b8r = b2r*drri - b2i*drii
                    b8i = b2i*drri + b2r*drii
                    ! %%%%%%%%%  forming integrands of j-matrices (j^{11}=j^{22}=0 for m.eq.0):

                    uri = dr(i)
                    rri = rr(i)

                    if (ncheck==1.and.si>0d0) goto 150
                    ! w(i)*r^2(i)*(pi(n1)*tau(n2)+tau(n1)*pi(n2):

                    e1 = rr(i)*aa1             ! <-- aa1

                    ar11 = ar11 + e1*b1r
                    ai11 = ai11 + e1*b1i
                    gr11 = gr11 + e1*c1r
                    gi11 = gi11 + e1*c1i

                    if (ncheck==1) goto 160

                    150 continue
                    ! w(i)*r^2(\theta)*[pi(n1)*pi(n2)+tau(n1)*tau(n2)]
                    ! (prefactor containing r^2(\theta)<->hat{r} part)

                    f1 = rri*aa2             ! <-- aa2
                    ! n1*(n1+1)*w(i)*r(\theta)*[dr/(d\theta)]*d(n1)*tau(n2):
                    !  (prefactor containing r(\theta)*[dr/(d\theta)] - hat{theta} part)

                    f2 = rri*uri*an1*a12             ! <-- a12

                    ar12 = ar12 + f1*b2r + f2*b3r        !~re j^{12}
                    ai12 = ai12 + f1*b2i + f2*b3i        !~im j^{12}

                    gr12 = gr12 + f1*c2r + f2*c3r        !~re rg j^{12}
                    gi12 = gi12 + f1*c2i + f2*c3i        !~im rg j^{12}
                    ! n2*(n2+1)*w(i)*r(\theta)*[dr/(d\theta)]*tau(n1)*d(n2):
                    ! (!prefactor containing r(\theta)*[dr/(d\theta)] - hat{theta} part)

                    f2 = rri*uri*an2*a21             ! <-- a21

                    ar21 = ar21 + f1*b4r + f2*b5r
                    ai21 = ai21 + f1*b4i + f2*b5i

                    gr21 = gr21 + f1*c4r + f2*c5r
                    gi21 = gi21 + f1*c4i + f2*c5i

                    if (ncheck==1) cycle

                    160  continue
                    ! w(i)*r^2(\theta)*[dr/(d\theta)]*pi(n1)*d(n2):
                    ! (!prefactor containing r^2(\theta)*[dr/(d\theta)] - hat{theta} part)

                    e2 = rri*uri*a11
                    e3 = e2*an2
                    e2 = e2*an1

                    ar22 = ar22 + e1*b6r + e2*b7r + e3*b8r
                    ai22 = ai22 + e1*b6i + e2*b7i + e3*b8i

                    gr22 = gr22 + e1*c6r + e2*c7r + e3*c8r
                    gi22 = gi22 + e1*c6i + e2*c7i + e3*c8i

                end do   !gauss integration
                !%%%%%%%%%%%%%  forming j-matrices (j^{11}=j^{22}=0 for m.eq.0):

                an12 = ann(n1, n2)*factor

                r11(n1, n2) = ar11*an12       !re j^{11}
                r12(n1, n2) = ar12*an12       !re j^{12}
                r21(n1, n2) = ar21*an12       !re j^{21}
                r22(n1, n2) = ar22*an12       !re j^{22}
                i11(n1, n2) = ai11*an12       !im j^{11}
                i12(n1, n2) = ai12*an12       !im j^{12}
                i21(n1, n2) = ai21*an12       !im j^{21}
                i22(n1, n2) = ai22*an12       !im j^{22}

                rg11(n1, n2) = gr11*an12       !re (rg j^{11})
                rg12(n1, n2) = gr12*an12       !re (rg j^{12})
                rg21(n1, n2) = gr21*an12       !re (rg j^{21})
                rg22(n1, n2) = gr22*an12       !re (rg j^{22})
                ig11(n1, n2) = gi11*an12       !im (rg j^{11})
                ig12(n1, n2) = gi12*an12       !im (rg j^{12})
                ig21(n1, n2) = gi21*an12       !im (rg j^{21})
                ig22(n1, n2) = gi22*an12       !im (rg j^{22})

            end do
        end do
        !%%%%%%%%%%%%%%%%%%%%%%%  forming q and rgq -matrices

        tpir = pir                 !re [1/k_{in}^2]
        tpii = pii                 !im [1/k_{in}^2]
        tppi = ppi                 !1/k_{out}^2

        nm = nmax - mm1 + 1
        do n1 = mm1, nmax
            k1 = n1 - mm1 + 1                !from 1 to nmax-mm1+1
            kk1 = k1 + nm                  !from nmax-mm1+2 to 2*(nmax-mm1+1)

            do n2 = mm1, nmax
                k2 = n2 - mm1 + 1           !from 1 to nmax-mm1+1
                kk2 = k2 + nm             !from nmax-mm1+2 to 2*(nmax-mm1+1)

                tar11 = -r11(n1, n2)
                tai11 = -i11(n1, n2)
                tgr11 = -rg11(n1, n2)
                tgi11 = -ig11(n1, n2)

                tar12 = i12(n1, n2)
                tai12 = -r12(n1, n2)
                tgr12 = ig12(n1, n2)
                tgi12 = -rg12(n1, n2)

                tar21 = -i21(n1, n2)
                tai21 = r21(n1, n2)
                tgr21 = -ig21(n1, n2)
                tgi21 = rg21(n1, n2)

                tar22 = -r22(n1, n2)
                tai22 = -i22(n1, n2)
                tgr22 = -rg22(n1, n2)
                tgi22 = -ig22(n1, n2)

                tqr(k1, k2) = tpir*tar21 - tpii*tai21 + tppi*tar12
                tqi(k1, k2) = tpir*tai21 + tpii*tar21 + tppi*tai12
                trgqr(k1, k2) = tpir*tgr21 - tpii*tgi21 + tppi*tgr12
                trgqi(k1, k2) = tpir*tgi21 + tpii*tgr21 + tppi*tgi12

                tqr(k1, kk2) = tpir*tar11 - tpii*tai11 + tppi*tar22
                tqi(k1, kk2) = tpir*tai11 + tpii*tar11 + tppi*tai22
                trgqr(k1, kk2) = tpir*tgr11 - tpii*tgi11 + tppi*tgr22
                trgqi(k1, kk2) = tpir*tgi11 + tpii*tgr11 + tppi*tgi22

                tqr(kk1, k2) = tpir*tar22 - tpii*tai22 + tppi*tar11
                tqi(kk1, k2) = tpir*tai22 + tpii*tar22 + tppi*tai11
                trgqr(kk1, k2) = tpir*tgr22 - tpii*tgi22 + tppi*tgr11
                trgqi(kk1, k2) = tpir*tgi22 + tpii*tgr22 + tppi*tgi11

                tqr(kk1, kk2) = tpir*tar12 - tpii*tai12 + tppi*tar21
                tqi(kk1, kk2) = tpir*tai12 + tpii*tar12 + tppi*tai21
                trgqr(kk1, kk2) = tpir*tgr12 - tpii*tgi12 + tppi*tgr21
                trgqi(kk1, kk2) = tpir*tgi12 + tpii*tgr12 + tppi*tgi21
            end do
        end do

        nnmax = 2*nm
        do n1 = 1, nnmax
            do n2 = 1, nnmax
                qr(n1, n2) = tqr(n1, n2)
                qi(n1, n2) = tqi(n1, n2)
                rgqr(n1, n2) = trgqr(n1, n2)
                rgqi(n1, n2) = trgqi(n1, n2)
            end do
        end do
        !
        call tt(nm)
        !
        return
    end

end program
! (c) copr. 09/1998  Alexander Moroz
! (c) copr. 10/2005  Alexander Moroz
