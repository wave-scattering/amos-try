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

    public cmplx_dp


contains
    !=======================================================================
    complex(dp) function cmplx_dp(re, im)
        real(dp), intent(in) :: re, im
        cmplx_dp = cmplx(re, im, kind = dp)
    end function cmplx_dp
    !=======================================================================
    !=======================================================================
    !=======================================================================
    !=======================================================================
    !=======================================================================
    !=======================================================================
    subroutine rsp_cylinder (x, r, dr)
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
        !                h=rev*( (2d0/(3d0*eps*eps))**(1d0/3d0) )
        !
        !
        !   ngauss ... the number of gif division points
        !   ng=2*ngauss
        !
        !   1.le.i.le.ngauss
        !
        !--------/---------/---------/---------/---------/---------/---------/--
        !TODO: convert to f90
        implicit none
        real(dp) rev, eps, h, a, si, co, rthet, rad
        integer ng, ngauss, i
        real(dp) x(:), r(:), dr(:)
        rev = mpar%rev
        eps = mpar%eps
        ! determine half-length of the cylinder
        h = rev * ((2d0 / (3d0 * eps * eps))**(1d0 / 3d0))
        ! determine cylinder radius:
        a = h * eps
        ng = size(x)
        ngauss = ng / 2
        if (ng == 1) ngauss = 1
        do i = 1, ngauss
            co = -x(i)
            si = dsqrt(1d0 - co * co)

            if ((h * si)>(a * co)) then
                ! along the circular surface:
                rad = a / si
                rthet = -a * co / (si * si)
                !rad=1.d-10
                !rthet=0.d0
            else
                ! along the plane cuts:
                rad = h / co
                rthet = h * si / (co * co)
            end if

            r(i) = rad * rad
            r(ng - i + 1) = r(i)          !using mirror symmetry

            dr(i) = -rthet / rad
            dr(ng - i + 1) = -dr(i)       !using mirror symmetry

        end do

        return
    end
    !=======================================================================
    subroutine drop (rat)
        !=================
        implicit none
        integer nout, nc, ng, i, n
        real(dp) s, v, rat, cki, dri, rv, ri, si, risi, rs, wi, xi, xin
        ! number of the output unit
        parameter (nout = 35)
        parameter (ng = 60)

        real(dp) x(ng), w(ng)
        nc = 10
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
        !
        ! gif division points and weights
        !
        call gauss (ng, 0, 0, x, w)
        !
        s = 0d0
        v = 0d0
        do i = 1, ng
            xi = dacos(x(i))
            wi = w(i)
            ri = 1d0 + cdrop%c(0)
            dri = 0d0
            do n = 1, nc
                xin = xi * n
                ri = ri + cdrop%c(n) * dcos(xin)
                dri = dri - cdrop%c(n) * n * dsin(xin)
            enddo
            si = dsin(xi)
            cki = x(i)
            risi = ri * si
            s = s + wi * ri * dsqrt(ri * ri + dri * dri)
            v = v + wi * ri * risi * (risi - dri * cki)
        enddo
        rs = dsqrt(s * 0.5d0)
        rv = (v * 3d0 * 0.25d0)**(1d0 / 3d0)
        if (dabs(rat - 1d0)>1d-8) rat = rv / rs
        cdrop%r0v = 1d0 / rv
        write(nout, 1000) cdrop%r0v
        do n = 0, nc
            write(nout, 1001) n, cdrop%c(n)
        enddo
        1000 format ('r_0/r_ev=', f7.4)
        1001 format ('c_', i2, '=', f7.4)

        return
    end
    !=======================================================================
    subroutine sareac (eps, rat)
        !--------/---------/---------/---------/---------/---------/---------/--
        ! >>> eps
        ! <<< rat
        !=================
        implicit none
        real(dp) :: rat, eps
        !
        rat = (1.5d0 / eps)**(1d0 / 3d0)
        rat = rat / dsqrt((eps + 2d0) / (2d0 * eps))
        !
        return
    end
    !=======================================================================
    subroutine sareananorod (eps, rat)
        !--------/---------/---------/---------/---------/---------/---------/--
        ! >>> eps
        ! <<< rat
        !=================
        !--------/---------/---------/---------/---------/---------/---------/--
        real(dp) :: rat, eps
        ! TODO: replace with computation of the nanorod surface
        rat = (1.5d0 / eps)**(1d0 / 3d0)
        rat = rat / dsqrt((eps + 2d0) / (2d0 * eps))
        !
        return
    end
    !=======================================================================
    subroutine sarea (d, rat)
        !--------/---------/---------/---------/---------/---------/---------/--
        ! >>> eps
        ! <<< rat
        !=================
        !--------/---------/---------/---------/---------/---------/---------/--
        real(dp) :: e, d, rat, r
        if (d < 1) then
            e = dsqrt(1d0 - d * d)
            r = 0.5d0 * (d**(2d0 / 3d0) + d**(-1d0 / 3d0) * dasin(e) / e)
            r = dsqrt(r)
            rat = 1d0 / r
            return
        else
            e = dsqrt(1d0 - 1d0 / (d * d))
            r = 0.25d0 * (2d0 * d**(2d0 / 3d0) &
                    + d**(-4d0 / 3d0) * dlog((1d0 + e) / (1d0 - e)) / e)
            r = dsqrt(r)
            rat = 1d0 / r
            !
        end if
        return
    end
    !=======================================================================
    subroutine surfch (n, e, rat)
        !--------/---------/---------/---------/---------/---------/---------/--
        ! >>> n,e,rat
        ! <<< rat
        !=================
        !--------/---------/---------/---------/---------/---------/---------/--
        !        implicit real(dp) (a-h, o-z)
        integer n, ng, i
        real(dp) e, rat, dn, s, v, xi, dx, dxn, ds, dsn, dcn, a, a2, en, ens, rs, rv
        real(dp) x(60), w(60)
        !
        dn = dble(n)
        en = e * dn
        ng = 60
        !
        ! gif division points and weights
        !
        call gauss (ng, 0, 0, x, w)
        !
        s = 0d0
        v = 0d0
        do i = 1, ng
            xi = x(i)
            dx = dacos(xi)
            dxn = dn * dx
            ds = dsin(dx)
            dsn = dsin(dxn)
            dcn = dcos(dxn)
            a = 1d0 + e * dcn
            a2 = a * a
            ens = en * dsn
            s = s + w(i) * a * dsqrt(a2 + ens * ens)
            v = v + w(i) * (ds * a + xi * ens) * ds * a2
        end do
        rs = dsqrt(s * 0.5d0)
        rv = (v * 3d0 / 4d0)**(1d0 / 3d0)
        rat = rv / rs
        !
        return
    end
    !=======================================================================
    subroutine gauleg(x1, x2, x, w, n)
        !--------/---------/---------/---------/---------/---------/---------/--
        !  given the lower and upper limits of integration x1 and x2, and given n
        !  this routine returns arrays x(1:n) and w(1:n) of length n, containing
        !  the abscissas and weights of the gaussian-legendre n-point quadrature
        !  formula.
        !--------/---------/---------/---------/---------/---------/---------/--
        integer, intent(in) :: n
        real(dp), intent(in) :: x1, x2
        real(dp), intent(out) :: x(n), w(n)
        real(dp) eps
        parameter (eps = 3.d-14)
        integer i, j, m
        double precision p1, p2, p3, pp, xl, xm, z, z1

        m = (n + 1) / 2    !the roots are symmetric in the interval, so we only
        xm = 0.5d0 * (x2 + x1)   !have to find half of them
        xl = 0.5d0 * (x2 - x1)

        ! loop over the desired roots:

        do i = 1, m
            z = cos(3.141592654d0 * (i - .25d0) / (n + .5d0))
            ! starting with the above approximation to the ith root, we enter
            ! the main loop of refinement by newton's method.
            do
                p1 = 1.d0
                p2 = 0.d0

                do j = 1, n !loop up the recurrence relation to get legendre
                    p3 = p2     !polynomial evaluated at z.
                    p2 = p1
                    p1 = ((2.d0 * j - 1.d0) * z * p2 - (j - 1.d0) * p3) / j
                end do

                ! p1 is now the desired  legendre polynomial. we next compute pp, its derivative,
                ! by a standard relation involving also p2, the polynomial of one lower order:

                pp = n * (z * p1 - p2) / (z * z - 1.d0)
                z1 = z
                z = z1 - p1 / pp                   !newton's method

                if (abs(z - z1)<=eps) exit
            end do
            ! scale the root to the desired interval, and put in its symmetric counterpart:
            x(i) = xm - xl * z
            x(n + 1 - i) = xm + xl * z
            ! compute the weight and its symmetric counterpart:
            w(i) = 2.d0 * xl / ((1.d0 - z * z) * pp * pp)
            w(n + 1 - i) = w(i)

        end do

        return
    end

    !=======================================================================
    subroutine gauss (n, ind1, ind2, z, w)
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
        !        implicit real(dp) (a-h, p-z)
        integer, intent(in) :: n, ind1, ind2
        integer i, k, ind, j, m, niter
        real(dp) x, a, b, c, check, dj, f, pa, pb, pc, zz
        real(dp), intent(out) :: z(:), w(:)
        a = 1d0
        b = 2d0
        c = 3d0
        !        data a, b, c /1d0, 2d0, 3d0/
        ind = mod(n, 2)
        k = n / 2 + ind
        f = dble(n)
        do i = 1, k
            m = n + 1 - i
            select case (i)
            case(1)
                x = a - b / ((f + a) * f)
            case(2)
                x = (z(n) - a) * 4d0 + z(n)
            case(3)
                x = (z(n - 1) - z(n)) * 1.6d0 + z(n - 1)
            case default
                x = (z(m + 1) - z(m + 2)) * c + z(m + 3)
            end select

            if(i==k.and.ind==1) x = 0d0
            niter = 0
            check = 1d-16
            do
                pb = 1d0
                niter = niter + 1
                if (niter>100) then
                    check = check * 10d0
                end if
                pc = x
                dj = a
                do j = 2, n
                    dj = dj + a
                    pa = pb
                    pb = pc
                    pc = x * pb + (x * pb - pa) * (dj - a) / dj
                end do
                pa = a / ((pb - x * pc) * f)
                pb = pa * pc * (a - x * x)
                x = x - pb
                if (dabs(pb)<=check * dabs(x)) exit
            end do
            z(m) = x
            w(m) = pa * pa * (a - x * x)
            if(ind1==0) w(m) = b * w(m)
            if(i/=k.or.ind/=1) then
                z(i) = -z(m)
                w(i) = w(m)
            end if
        end do
        if (ind2==1) then
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
            if(ind1/=0) then
                do i = 1, n
                    z(i) = (a + z(i)) / b
                end do
            end if
        end if

        return
    end
    !=======================================================================

end module libcylinder
