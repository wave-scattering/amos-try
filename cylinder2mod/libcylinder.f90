module libcylinder
    use constants
    use model_parameters
    use special_functions
    use cylinder_blas
    !use dense_solve
    !use multem_blas
    !use amos
    !use errfun, only : wpop

    implicit none

    type, private :: cdrop_values
        real(dp) :: c(0:10), r0v
    end type cdrop_values
    type(cdrop_values), public :: cdrop

    public cmplx_dp, rsp_spheroid, rsp_cylinder, rsp_nanorod, &
            rsp_chebyshev, rsp_droplet, radii_ratio_spheroid, &
            radii_ratio_cylinder, radii_ratio_nanorod, &
            radii_ratio_chebyshev, radii_ratio_droplet


contains
    !=======================================================================
    complex(dp) function cmplx_dp(re, im)
        real(dp), intent(in) :: re, im
        cmplx_dp = cmplx(re, im, kind = dp)
    end function cmplx_dp

    !=======================================================================
    subroutine rsp_spheroid(x, r, dr)
        !--------/---------/---------/---------/---------/---------/---------/--
        ! >>> x,ng,ngauss,rev,eps
        ! <<< r,dr
        !=========================
        !   activated for np=-1
        !
        !   calculation of the functions
        !              r(i)=r(y)**2 and dr(i)=((d/dy)r(y))/r(y)
        !   for an oblate/prolate spheroids droplet specified by the parameters
        !   rev and eps at ngauss gauss integration formula (gif) division points
        !   in the integral over theta. here y=acos(x)=theta.
        !
        !      r(\theta,\phi)=a\left[\sin^2\theta + (a^2/b^2)\cos^2\theta]^{-1/2}
        !
        !   x - gif division points \cos\theta_j -  y = arccos x
        !   rev ... equal-volume-sphere radius
        !   eps ... the ratio of the horizontal to rotational axes.  eps is
        !           larger than 1 for oblate spheroids and smaller than 1 for
        !           prolate spheroids.
        !   ngauss ... the number of gif division points
        !   ng=2*ngauss
        !
        !   1 <= i <= ngauss
        !--------/---------/---------/---------/---------/---------/---------/--

        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: r(:), dr(:)

        integer ng, ngauss, i
        real(dp) rev, eps
        real(dp) a, aa, c, cc, ee, ee1, rr
        real(dp) s, ss

        ! Assign model parameters
        rev = mpar%rev
        eps = mpar%eps

        ! Set ng and ngauss
        ng = size(x)
        ngauss = ng/2
        if (ng == 1) ngauss = 1

        a = rev*eps**(1d0/3d0)
        aa = a*a
        ee = eps*eps
        ee1 = ee - 1d0

        do i = 1, ngauss
            c = x(i)
            cc = c*c
            ss = 1d0 - cc
            s = dsqrt(ss)                 !=\sin\theta
            rr = 1d0/(ss + ee*cc)
            r(i) = aa*rr
            r(ng - i + 1) = r(i)
            dr(i) = rr*c*s*ee1
            dr(ng - i + 1) = -dr(i)
        end do

        return
    end

    !=======================================================================
    subroutine rsp_cylinder(x, r, dr)
        !--------/---------/---------/---------/---------/---------/---------/--
        ! >>> x,ng,ngauss,rev,eps
        ! <<< r,dr
        !=========================
        !   activated for np=-2
        !
        !   calculation of the functions r(i)=r(y)**2 and
        !   dr(i)=((d/dy)r(y))/r(y) for an oblate/prolate cylinder
        !   specified by the parameters rev and eps  at ngauss  gauss
        !   integration points in the integral over theta.
        !
        !   x - gif division points \cos\theta_j -  y = arccos x
        !   rev ... equal-volume-sphere radius r_ev
        !   eps ... the ratio of the cylinder diameter to its length
        !   h   ... half-length of the cylinder
        !   a=h*eps  ... cylinder radius   ====>
        !
        !   4*pi*rev**3/3=2*h*pi*a**2=2*pi*h**3*eps**2 <====>
        !                h=rev*( (2_dp/(3_dp*eps*eps))**(1_dp/3_dp) )
        !
        !
        !   ngauss ... the number of gif division points
        !   ng=2*ngauss
        !
        !   1 <= i <= ngauss
        !
        !--------/---------/---------/---------/---------/---------/---------/--

        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: r(:), dr(:)

        integer ng, ngauss, i
        real(dp) rev, eps
        real(dp) h, a, si, co, rthet, rad

        ! Assign model parameters
        rev = mpar%rev
        eps = mpar%eps

        ! Set ng and ngauss
        ng = size(x)
        ngauss = ng/2
        if (ng == 1) ngauss = 1

        ! determine half-length of the cylinder
        h = rev*((2.0_dp/(3.0_dp*eps*eps))**(1.0_dp/3.0_dp))
        ! determine cylinder radius:
        a = h*eps

        do i = 1, ngauss
            co = -x(i)
            si = dsqrt(1.0_dp - co*co)

            if ((h*si)>(a*co)) then
                ! along the circular surface:
                rad = a/si
                rthet = -a*co/(si*si)
                !rad=1.d-10
                !rthet=0.0_dp
                !write(*,*) 'cyl'
            else
                ! along the plane cuts:
                rad = h/co
                rthet = h*si/(co*co)
            end if

            r(i) = rad*rad
            r(ng - i + 1) = r(i)          !using mirror symmetry

            dr(i) = -rthet/rad
            dr(ng - i + 1) = -dr(i)       !using mirror symmetry

        end do

        return
    end

    !=======================================================================
    subroutine rsp_nanorod(x, r, dr)
        !--------/---------/---------/---------/---------/---------/---------/--
        ! >>> x,ng,ngauss,rev,eps,epse
        ! <<< r,dr
        !=========================
        !   activated for np=-9
        !
        !   calculation of the functions r(i)=r(y)**2 and
        !   dr(i)=((d/dy)r(y))/r(y) for a nanorod particle
        !   specified by the parameters rev, eps, and cap  at ngauss  gauss
        !   integration points in the integral over theta.
        !
        !   x - gif division points \cos\theta_j -  y = arccos x
        !   rev ... equal-volume-sphere radius r_ev
        !   eps ... the ratio of the nanorod diameter to its length
        !   epse ... the ratio of half-spheroid caps height to cylider radius
        !   epsc ... the ration of cylinder half-height to cylinder radius
        !   h   ... half-length of the nanorod
        !   he .. cap height
        !   hc .. cylinder height
        !   a  ... cylinder radius   ====>
        !
        !   ngauss ... the number of gif division points
        !   ng=2*ngauss
        !
        !   1 <= i <= ngauss
        !
        !--------/---------/---------/---------/---------/---------/---------/--

        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: r(:), dr(:)

        integer ng, ngauss, i
        real(dp) rev, eps, epse
        real(dp) h, a, co, si, rad, rthet
        real(dp) valdr, c2, s2, alpha, beta, he, hc, epsc

        ! Assign model parameters
        rev = mpar%rev
        eps = mpar%eps
        epse = mpar%nanorod_cap_hr

        ! Set ng and ngauss
        ng = size(x)
        ngauss = ng/2
        if (ng == 1) ngauss = 1

        !rev = 2._dp
        !eps = 0.5_dp
        !epse = 0.65_dp
        !xstep = pi/(ng-1)
        !do i = 1,ng, 1
        !    todo: remove testing setting of x
        !    x(i)=cos(pi-(i-1)*xstep)
        !end do

        !cylinder radius
        a = rev*(2d0*eps/(3d0 - eps*epse)) ** (1d0/3d0)
        h = a/eps      ! nanorod half-height
        he = a*epse  ! spheroid cap half-height
        hc = h - he     ! cylinder half-height
        epsc = hc/a    ! cylider aspect ratio

        do i = 1, ngauss, 1
            co = -x(i)
            si = dsqrt(1_dp - co*co)

            if ((hc*si)>(a*co)) then
                ! along the circular surface:
                rad = a/si
                rthet = -a*co/(si*si)
                valdr = -rthet/rad
            else
                !  along elliptic cap
                c2 = co**2
                s2 = si**2
                ! solution of square equation of ellipse move from the origin
                alpha = dsqrt((epse**2 - epsc**2)*s2 + c2)
                beta = epse**2*s2 + c2
                rad = (hc*co + he*alpha)/beta
                valdr = -((-alpha*hc*si&
                        + he*(epse**2 - epsc**2 - 1._dp)*si*co&
                        )/(he*alpha**2 + alpha*hc*co)&
                        - (2*(epse**2 - 1._dp)*si*co/beta))

            endif
            r(i) = rad*rad
            r(ng - i + 1) = r(i)          !using mirror symmetry

            dr(i) = valdr
            dr(ng - i + 1) = -dr(i)       !using mirror symmetry

        enddo

        return
    end

    !=======================================================================
    subroutine rsp_chebyshev(x, n, r, dr)
        !--------/---------/---------/---------/---------/---------/---------/--
        ! >>> x,ng,rev,eps,n
        ! <<< r,dr
        !=========================
        !   activated for np > 0
        !
        !   calculation of the functions r(i)=r(y)**2 and
        !   dr(i)=((d/dy)r(y))/r(y) for a chebyshev particle
        !   specified by the parameters rev, eps, and n,
        !
        !       r(\theta,\phi)=r_0[1+\eps t_n(\cos\theta)]    (*)
        !
        !   eps ... deformation parameter of a chebyshev particle; |eps|<1
        !   n   ... the degree of the chebyshev polynomial
        !   all chebyshev particles with n.ge.2 become partially concave
        !   as the absolute value of the deformation parameter eps increases
        !   and exhibit surface roughness in the form of waves running
        !   completely around the particle.
        !
        !   x - gif division points \cos\theta_j -  y = arccos x
        !   rev ... equal-volume-sphere radius r_ev
        !   ngauss ... the number of gif division points
        !   ng=2*ngauss
        !
        !   1 <= i <= ngauss
        !
        !--------/---------/---------/---------/---------/---------/---------/--

        integer, intent(in) :: n
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: r(:), dr(:)

        integer ng, i
        real(dp) rev, eps
        real(dp) ep, a, r0, ri, xi
        real(dp) dn, dnp, dn4

        ! Assign model parameters
        rev = mpar%rev
        eps = mpar%eps

        ! Set ng
        ng = size(x)

        dnp = dble(n)
        dn = dnp*dnp
        dn4 = dn*4d0
        ep = eps*eps
        a = 1d0 + 1.5d0*ep*(dn4 - 2d0)/(dn4 - 1d0)
        i = (dnp + 0.1d0)*0.5d0
        i = 2*i
        if (i ==n ) then
            a = a - 3d0*eps*(1d0 + 0.25d0*ep)/ &
                    (dn - 1d0) - 0.25d0*ep*eps/(9d0*dn - 1d0)
        end if
        r0 = rev*a**(-1d0/3d0)
        do i = 1, ng
            xi = dacos(x(i))*dnp
            ri = r0*(1d0 + eps*dcos(xi))    !the chebyshev shape function (*)
            r(i) = ri*ri
            dr(i) = -r0*eps*dnp*dsin(xi)/ri
            !write(nout,*) i,r(i),dr(i)
        end do

        return
    end

    !=======================================================================
    subroutine rsp_droplet(x, r, dr)
        !--------/---------/---------/---------/---------/---------/---------/--
        ! >>> x,ng,rev
        ! <<< r,dr
        !=========================
        !   activated for np=-3
        !
        !   calculation of the functions r(i)=r(y)**2 and
        !   dr(i)=((d/dy)r(y))/r(y) for a distorted
        !   droplet (generalized chebyshev particle) specified by the
        !   parameters rev and c_n (chebyshev expansion coefficients).
        !   the coefficients of the chebyshev polynomial expansion are
        !   specified in the subroutine drop.
        !
        !   x - gif division points \cos\theta_j -  y = arccos x
        !   rev ... equal-volume-sphere radius  r_ev
        !   ngauss ... the number of gif division points
        !   ng=2*ngauss
        !
        !   1 <= i <= ngauss
        !
        !--------/---------/---------/---------/---------/---------/---------/--

        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: r(:), dr(:)

        integer, parameter :: nc = 10

        integer ng, n, i
        real(dp) rev
        real(dp) dri, r0, ri, xi, xin

        ! Assign model parameters
        rev = mpar%rev

        ! Set ng
        ng = size(x)

        r0 = rev*cdrop%r0v
        do i = 1, ng
            xi = dacos(x(i))
            ri = 1d0 + cdrop%c(0)
            dri = 0d0
            do n = 1, nc
                xin = xi*n
                ri = ri + cdrop%c(n)*dcos(xin)
                dri = dri - cdrop%c(n)*n*dsin(xin)
            enddo
            ri = ri*r0
            dri = dri*r0
            r(i) = ri*ri
            dr(i) = dri/ri
            !write(nout,*) i,r(i),dr(i)
        enddo

        return
    end

    !=======================================================================
    function radii_ratio_spheroid() result(rat)
        !--------/---------/---------/---------/---------/---------/---------/--
        ! >>> eps
        ! <<< rat
        !=========================
        !   activated for np=-1
        !
        !   calculation of the ratio between the volume-equivalent and
        !   the surface-equivalent radii for the spheroid

        real(dp) :: eps, rat
        real(dp) :: e, r

        ! Assign model parameters
        eps = mpar%eps

        if (eps < 1d0) then ! Prolate spheroid
            e = dsqrt(1d0 - eps*eps)
            r = 0.5d0*(eps**(2d0/3d0) + eps**(-1d0/3d0)*dasin(e)/e)
        else ! Oblate spheroid
            e = dsqrt(1d0 - 1d0/(eps*eps))
            r = 0.25d0*(2d0*eps**(2d0/3d0) &
                    + eps**(-4d0/3d0)*dlog((1d0 + e)/(1d0 - e))/e)
        end if

        rat = 1d0/dsqrt(r)
    end function radii_ratio_spheroid

    !=======================================================================
    function radii_ratio_cylinder() result(rat)
        !--------/---------/---------/---------/---------/---------/---------/--
        ! >>> eps
        ! <<< rat
        !=========================
        !   activated for np=-2
        !
        !   calculation of the ratio between the volume-equivalent and
        !   the surface-equivalent radii for the cylinder (rv/rs)

        implicit none
        real(dp) :: eps, rat

        ! Assign model parameters
        eps = mpar%eps

        rat = (1.5d0/eps)**(1d0/3d0)/ &
                dsqrt((eps + 2d0)/(2d0*eps))
    end function radii_ratio_cylinder

    !=======================================================================
    function radii_ratio_nanorod() result(rat)
        !--------/---------/---------/---------/---------/---------/---------/--
        ! >>> eps, epse
        ! <<< rat
        !=========================
        !   activated for np=-9
        !
        !   calculation of the ratio between the volume-equivalent and
        !   the surface-equivalent radii for the nanorod

        real(dp) :: eps, epse, epsc, rat
        real(dp) :: e, rv, rs

        ! Assign model parameters
        eps = mpar%eps
        epse = mpar%nanorod_cap_hr

        ! TODO: replace with computation of the nanorod surface
        ! (2019-12-09) In principle it is done but needs some testing
        epsc = (1d0/eps) - epse

        rv = (1.5d0*epsc + epse)**(1d0/3d0)

        if (dabs(epse) <= 1d-6) then ! Plane caps (ie, cylinder)
            rs = dsqrt(epsc + 0.5d0)
        elseif (dabs(epse - 1d0) <= 1d-6) then ! Spherical caps
            rs = dsqrt(epsc + 1d0)
        elseif (epse > 1d0) then ! Prolate spheroid
            e = dsqrt(1d0 - 1d0/(epse*epse))
            rs = dsqrt(epsc + 0.5d0*(1d0 + epse*dasin(e)/e))
        else ! Oblate spheroid
            e = dsqrt(1d0 - eps*eps)
            rs = dsqrt(epsc + 0.5d0*(1d0 + epse*dlog((1d0 + e)/(1d0 - e))/e))
        end if

        rat = rv/rs
    end function radii_ratio_nanorod

    !=======================================================================
    function radii_ratio_chebyshev(n) result(rat)
        !--------/---------/---------/---------/---------/---------/---------/--
        ! >>> n, eps
        ! <<< rat
        !   activated for np>=0
        !
        !   calculation of the ratio between the volume-equivalent and
        !   the surface-equivalent radii for a chebyshev particle

        integer, parameter :: ng = 60

        integer, intent(in) :: n

        integer :: i
        real(dp) :: eps, rat
        real(dp) :: dn, s, v, xi, dx, dxn, ds, &
                dsn, dcn, a, a2, en, ens, rs, rv
        real(dp) :: x(ng), w(ng)

        ! Assign model parameters
        eps = mpar%eps

        dn = dble(n)
        en = eps*dn

        ! gif division points and weights
        call gauss(ng, 0, 0, x, w)

        s = 0d0
        v = 0d0
        do i = 1, ng
            xi = x(i)
            dx = dacos(xi)
            dxn = dn*dx
            ds = dsin(dx)
            dsn = dsin(dxn)
            dcn = dcos(dxn)
            a = 1d0 + eps*dcn
            a2 = a*a
            ens = en*dsn
            s = s + w(i)*a*dsqrt(a2 + ens*ens)
            v = v + w(i)*(ds*a + xi*ens)*ds*a2
        end do
        rs = dsqrt(s*0.5d0)
        rv = (v*3d0/4d0)**(1d0/3d0)

        rat = rv/rs
    end function radii_ratio_chebyshev

    !=======================================================================
    function radii_ratio_droplet() result(rat)
        !=========================
        !   activated for np>=0
        !
        !   calculation of the ratio between the volume-equivalent and
        !   the surface-equivalent radii for a distorted chebyshev droplet

        ! nout: number of the output unit
        integer, parameter :: nout = 35, ng = 60, nc = 10

        integer i, n
        real(dp) s, v, rat, cki, dri, rv, ri, si, risi, rs, wi, xi, xin
        real(dp) x(ng), w(ng)

        cdrop%c(0) = -0.0481_dp
        cdrop%c(1) = 0.0359_dp
        cdrop%c(2) = -0.1263_dp
        cdrop%c(3) = 0.0244_dp
        cdrop%c(4) = 0.0091_dp
        cdrop%c(5) = -0.0099_dp
        cdrop%c(6) = 0.0015_dp
        cdrop%c(7) = 0.0025_dp
        cdrop%c(8) = -0.0016_dp
        cdrop%c(9) = -0.0002_dp
        cdrop%c(10) = 0.0010_dp

        ! gif division points and weights
        call gauss(ng, 0, 0, x, w)

        s = 0d0
        v = 0d0
        do i = 1, ng
            xi = dacos(x(i))
            wi = w(i)
            ri = 1d0 + cdrop%c(0)
            dri = 0d0
            do n = 1, nc
                xin = xi*n
                ri = ri + cdrop%c(n)*dcos(xin)
                dri = dri - cdrop%c(n)*n*dsin(xin)
            enddo
            si = dsin(xi)
            cki = x(i)
            risi = ri*si
            s = s + wi*ri*dsqrt(ri*ri + dri*dri)
            v = v + wi*ri*risi*(risi - dri*cki)
        enddo
        rs = dsqrt(s*0.5d0)
        rv = (v*3d0*0.25d0)**(1d0/3d0)

        cdrop%r0v = 1d0/rv
        write(nout, 1000) cdrop%r0v
        do n = 0, nc
            write(nout, 1001) n, cdrop%c(n)
        enddo
        1000 format('r_0/r_ev=', f7.4)
        1001 format('c_', i2, '=', f7.4)

        if (dabs(rat - 1d0) > 1d-8) then
            rat = rv/rs
        else
            rat = 1d0
        end if
    end function radii_ratio_droplet

    !=======================================================================
    subroutine gauleg(x1, x2, x, w, n)
        !--------/---------/---------/---------/---------/---------/---------/--
        !  given the lower and upper limits of integration x1 and x2, and given n
        !  this routine returns arrays x(1:n) and w(1:n) of length n, containing
        !  the abscissas and weights of the Gaussian-Legendre n-point quadrature
        !  formula.
        !--------/---------/---------/---------/---------/---------/---------/--
        integer, intent(in) :: n
        real(dp), intent(in) :: x1, x2
        real(dp), intent(out) :: x(n), w(n)
        real(dp), parameter :: eps = 3.d-14
        integer i, j, m
        double precision p1, p2, p3, pp, xl, xm, z, z1

        m = (n + 1)/2    !the roots are symmetric in the interval, so we only
        xm = 0.5d0*(x2 + x1)   !have to find half of them
        xl = 0.5d0*(x2 - x1)

        ! loop over the desired roots:

        do i = 1, m
            z = cos(3.141592654d0*(i - .25d0)/(n + .5d0))
            ! starting with the above approximation to the ith root, we enter
            ! the main loop of refinement by newton's method.
            do
                p1 = 1.d0
                p2 = 0.d0

                do j = 1, n !loop up the recurrence relation to get legendre
                    p3 = p2     !polynomial evaluated at z.
                    p2 = p1
                    p1 = ((2.d0*j - 1.d0)*z*p2 - (j - 1.d0)*p3)/j
                end do

                ! p1 is now the desired  legendre polynomial. we next compute pp, its derivative,
                ! by a standard relation involving also p2, the polynomial of one lower order:

                pp = n*(z*p1 - p2)/(z*z - 1.d0)
                z1 = z
                z = z1 - p1/pp                   !newton's method

                if (abs(z - z1)<=eps) exit
            end do
            ! scale the root to the desired interval, and put in its symmetric counterpart:
            x(i) = xm - xl*z
            x(n + 1 - i) = xm + xl*z
            ! compute the weight and its symmetric counterpart:
            w(i) = 2.d0*xl/((1.d0 - z*z)*pp*pp)
            w(n + 1 - i) = w(i)

        end do

        return
    end

    !=======================================================================
    subroutine gauss(n, ind1, ind2, z, w)
        !--------/---------/---------/---------/---------/---------/---------/--
        ! >>> n,ind1,ind2
        ! <<< z,w
        !=================
        !    calculation of points and weights of gaussian quadrature
        !    formula. if ind1 = 0 - on interval (-1,1), if ind1 = 1 - on
        !    interval  (0,1). if  ind2 = 1 results are printed.
        !
        !    n - number of gif division points (mostly n=ngauss in main program)
        !    z - division points
        !    w - weights
        !--------/---------/---------/---------/---------/---------/---------/--

        integer, intent(in) :: n, ind1, ind2
        real(dp), intent(out) :: z(:), w(:)

        integer i, k, ind, j, m, niter
        real(dp) x, a, b, c, check, dj, f, pa, pb, pc, zz

        a = 1d0
        b = 2d0
        c = 3d0
        !data a, b, c /1d0, 2d0, 3d0/
        ind = mod(n, 2)
        k = n/2 + ind
        f = dble(n)
        do i = 1, k
            m = n + 1 - i
            select case (i)
            case(1)
                x = a - b/((f + a)*f)
            case(2)
                x = (z(n) - a)*4d0 + z(n)
            case(3)
                x = (z(n - 1) - z(n))*1.6d0 + z(n - 1)
            case default
                x = (z(m + 1) - z(m + 2))*c + z(m + 3)
            end select

            if(i == k .and. ind ==1 ) x = 0d0
            niter = 0
            check = 1d-16
            do
                pb = 1d0
                niter = niter + 1
                if (niter > 100) then
                    check = check*10d0
                end if
                pc = x
                dj = a
                do j = 2, n
                    dj = dj + a
                    pa = pb
                    pb = pc
                    pc = x*pb + (x*pb - pa)*(dj - a)/dj
                end do
                pa = a/((pb - x*pc)*f)
                pb = pa*pc*(a - x*x)
                x = x - pb
                if (dabs(pb)<=check*dabs(x)) exit
            end do
            z(m) = x
            w(m) = pa*pa*(a - x*x)
            if(ind1==0) w(m) = b*w(m)
            if(i/=k.or.ind/=1) then
                z(i) = -z(m)
                w(i) = w(m)
            end if
        end do
        if (ind2 == 1) then
            print 1100, n
            1100 format(' ***  points and weights of gaussian quadrature formula', &
                    ' of ', i4, '-th order')
            do i = 1, k
                zz = -z(i)
                print 1200, i, zz, i, w(i)
            end do
            1200 format(' ', 4x, 'x(', i4, ') = ', f17.14, 5x, 'w(', i4, ') = ', f17.14)
            !     print 1300,n
            ! 1300 format(' gaussian quadrature formula of ',i4,'-th order is used')
        else
            if(ind1 /= 0) then
                do i = 1, n
                    z(i) = (a + z(i))/b
                end do
            end if
        end if

        return
    end
    !=======================================================================

end module libcylinder
