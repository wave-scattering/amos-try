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

    type, private :: integrand_values
        integer :: n1, n2, nmax
    end type integrand_values
    type(integrand_values), public :: integrand

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
    function f03 ( x )

        !*****************************************************************************80
        !
        !! F03 is the integrand function LOG(X)/SQRT(X).
        !
        !  Modified:
        !
        !    11 September 2015
        !
        !  Author:
        !
        !    John Burkardt
        !
        implicit none

        real ( kind = 8 ) f03
        real ( kind = 8 ) x

        if ( x <= 0.0D+00 ) then
            f03 = 0.0D+00
        else
            f03 = log ( x ) / sqrt ( x )
        end if

        return
    end

    !=======================================================================
    subroutine integrate_p (fun, a, b, result)
        !    DQAGP is an adaptive integrator that can handle singularities
        !    of the integrand at user specified points,
        implicit none

        integer ( kind = 4 ), parameter :: limit = 1500
        integer ( kind = 4 ), parameter :: npts = 1
        integer ( kind = 4 ), parameter :: npts2 = 2 * npts

        integer ( kind = 4 ), parameter :: leniw = 2 * limit + npts2
        integer ( kind = 4 ), parameter :: lenw = leniw * 2 - npts2

        real ( kind = 8 ) :: a, b
        real ( kind = 8 ) abserr
        real ( kind = 8 ), parameter :: epsabs = 0.0D+00
        real ( kind = 8 ), parameter :: epsrel = 1.0D-5
        real ( kind = 8 ), external :: fun
        integer ( kind = 4 ) ier
        integer ( kind = 4 ) iwork(leniw)
        integer ( kind = 4 ) last
        integer ( kind = 4 ) neval
        real ( kind = 8 ) points(npts2)
        real ( kind = 8 ) result
        real ( kind = 8 ) work(lenw)
        !
        !  Singularity points:
        !
        points(1) = 1.0_dp/dsqrt(mpar%eps**2 +1.0_dp)
        write(*,*) 'Points = ', points

        call dqagp ( fun, a, b, npts2, points, epsabs, epsrel, result, abserr, &
                neval, ier, leniw, lenw, last, iwork, work )

!        write ( *, '(a,g14.6)' ) '  Integral left endpoint A =    ', a
!        write ( *, '(a,g14.6)' ) '  Integral right endpoint B =   ', b
        write ( *, '(a,g14.6)' ) '  Estimated integral is         ', result
        write ( *, '(a,g14.6)' ) '  Estimated integral error =    ', abserr
        write ( *, '(a,i8)' ) '  Number of function evaluations, NEVAL = ', neval
!        write ( *, '(a,i8)' ) '  Error return code IER = ', ier

        return
    end
    !=======================================================================
    subroutine integrate (fun, a, b, result)
        !    DQAGS is an adaptive integrator for endpoint singularities.
        !
        implicit none
        integer ( kind = 4 ), parameter :: limit = 500
        integer ( kind = 4 ), parameter :: lenw = limit * 4

        real ( kind = 8 ) :: a,b
        real ( kind = 8 ) abserr
        real ( kind = 8 ), parameter :: epsabs = 0.0D+00
        real ( kind = 8 ), parameter :: epsrel = 1.0D-2
        real ( kind = 8 ), external :: fun
        integer ( kind = 4 ) ier
        integer ( kind = 4 ) iwork(limit)
        integer ( kind = 4 ) last
        integer ( kind = 4 ) neval
        real ( kind = 8 ) result
        real ( kind = 8 ) work(lenw)

        call dqags ( fun, a, b, epsabs, epsrel, result, abserr, neval, ier, &
                limit, lenw, last, iwork, work )

!            write ( *, '(a,g14.6)' ) '  Integral left endpoint A =    ', a
!            write ( *, '(a,g14.6)' ) '  Integral right endpoint B =   ', b
!            write ( *, '(a,g14.6)' ) '  Estimated integral is         ', result
!            write ( *, '(a,g14.6)' ) '  Estimated integral error =    ', abserr
!            write ( *, '(a,i8)' ) '  Number of function evaluations, NEVAL = ', neval
!            write ( *, '(a,i8)' ) '  Error return code IER = ', ier
!        if (ier /= 0 .and. ier/=3) then
!            stop 'dqags integration error'
!        end if

        return
    end

    !=======================================================================
    function ar12_m0_integrand(xi)
        implicit none
        integer n1, n2
        real(dp) ar12_m0_integrand(2)
        real(dp) x_to_rsp(1), r_from_x(1), dr_from_x(1)
        real(dp) xi, d1n1, d2n1, d1n2, d2n2, a12, a22, qj1, qdj1, &
                qy1, qdy1, qjr2, qji2, qdjr2, qdji2, c1r, bir, c2r, &
                b2r, ddri, b3r, rr1, f1, f2, an1, b1r, uri
        n1 = integrand%n1
        n2 = integrand%n2
        an1 = dble(n1*(n1+1))
        if (mpar%np.ne.-2) stop 'adaptive integration is only implemented for cylinder'
        x_to_rsp(1) = xi
        call rsp_cylinder(x_to_rsp, r_from_x, dr_from_x)
        call vig_1v ( xi, n1, 0, d1n1, d2n1)
        call vig_1v ( xi, n2, 0, d1n2, d2n2)
        a12=d1n1*d2n2
        a22=d2n1*d2n2
        call cbessjdj(r_from_x(1),n1, qj1, qdj1)
        call cbessydy(r_from_x(1),n1, qy1, qdy1)
        ! bessel functions of the interior argument:
        call cbesscjcdj(r_from_x(1),n2, integrand%nmax, qjr2, qji2, qdjr2, qdji2)
        ! re and im of j_{n2}(k_{in}r) j_{n1}(k_{out}r):

        c1r=qjr2*qj1
        ! re and im of j_{n2}(k_{in}r) h_{n1}(k_{out}r):
        b1r=c1r-qji2*qy1
        ! re and im of j_{n2}(k_{in}r) [k_{out}r j_{n1}(k_{out}r)]'/(k_{out}r):
        c2r=qjr2*qdj1
        ! re and im of j_{n2}(k_{in}r) [k_{out}r h_{n1}(k_{out}r)]'/(k_{out}r):
        b2r=c2r-qji2*qdy1
        ddri=1.0_dp/(dsqrt(r_from_x(1))*cbess%wv)
        ! re and im of [1/(k_{out}r)]*j_{n2}(k_{in}r) h_{n1}(k_{out}r):
        b3r=ddri*b1r
        !%%%%%%%  forming integrands of j-matrices (j^{11}=j^{22}=0 for m=0): %%%%%%%%
        !                   uri=dr(i)        !dr/(d\theta)
        uri=-dr_from_x(1)  ! todo: why 'minus' sign was need to fit previous line?
        !        rri=rr(i)        !w(i)*r^2(\theta)
        ! w(i)*r^2(\theta)*d2n1*d2n2:
        f1=r_from_x(1)*a22      !prefactor containing r^2(\theta)<->hat{r} part
        ! n1*(n1+1)*w(i)*r(\theta)*[dr/(d\theta)]*d1n1*d2n2:
        f2=r_from_x(1)*uri*an1*a12     !prefactor containing r(\theta)*[dr/(d\theta)]
        !hat{theta} part
        ar12_m0_integrand = (/f1*b2r+f2*b3r, 0.0_dp/)        !~re j^{12}
        return
    end function ar12_m0_integrand
    !=======================================================================
    function m0_integrand(xi)
        implicit none
        integer n1, n2
        real(dp) m0_integrand(8)
        real(dp) x_to_rsp(1), r_from_x(1), dr_from_x(1)
        real(dp) xi, d1n1, d2n1, d1n2, d2n2, a12, a21, a22, qj1, qdj1, &
                qy1, qdy1, qjr2, qji2, qdjr2, qdji2, c1r, c1i, b1r, b1i, c2r, c2i, &
                b2r, b2i,ddri, c3r, c3i, b3r, b3i, c4r, c4i, b4r, b4i, &
                drri, drii, c5r, c5i, b5r, b5i, rr1, f1, f2, an1, an2, uri, &
                ar12, ai12, gr12, gi12, ar21, ai21, gr21, gi21, v

        n1 = integrand%n1
        n2 = integrand%n2
        an1 = dble(n1*(n1+1))
        an2 = dble(n2*(n2+1))
        if (mpar%np.ne.-2) stop 'adaptive integration is only implemented for cylinder'
        x_to_rsp(1) = xi
        call rsp_cylinder(x_to_rsp, r_from_x, dr_from_x)
        call vig_1v ( xi, n1, 0, d1n1, d2n1)
        call vig_1v ( xi, n2, 0, d1n2, d2n2)
        a12 = d1n1 * d2n2
        a21 = d2n1 * d1n2
        a22 = d2n1 * d2n2

        call cbessjdj(r_from_x(1),n1, qj1, qdj1)
        call cbessydy(r_from_x(1),n1, qy1, qdy1)
        ! bessel functions of the interior argument:
        call cbesscjcdj(r_from_x(1),n2, integrand%nmax, qjr2, qji2, qdjr2, qdji2)
        ! re and im of j_{n2}(k_{in}r) j_{n1}(k_{out}r):

        c1r = qjr2 * qj1
        c1i = qji2 * qj1
        ! re and im of j_{n2}(k_{in}r) h_{n1}(k_{out}r):
        b1r = c1r - qji2 * qy1
        b1i = c1i + qjr2 * qy1
        ! re and im of j_{n2}(k_{in}r) [k_{out}r j_{n1}(k_{out}r)]'/(k_{out}r):
        c2r = qjr2 * qdj1
        c2i = qji2 * qdj1
        ! re and im of j_{n2}(k_{in}r) [k_{out}r h_{n1}(k_{out}r)]'/(k_{out}r):
        b2r = c2r - qji2 * qdy1
        b2i = c2i + qjr2 * qdy1

        ddri=1.0_dp/(dsqrt(r_from_x(1))*cbess%wv) !1/(k_{out}r)
        ! re and im of [1/(k_{out}r)]*j_{n2}(k_{in}r) j_{n1}(k_{out}r)

        c3r = ddri * c1r
        c3i = ddri * c1i
        ! re and im of [1/(k_{out}r)]*j_{n2}(k_{in}r) h_{n1}(k_{out}r):

        b3r = ddri * b1r
        b3i = ddri * b1i
        ! re and im of  [k_{in}r j_{n2}(k_{in}r)]'/(k_{in}r)
        !                          * j_{n1}(k_{out}r):

        c4r = qdjr2 * qj1
        c4i = qdji2 * qj1
        ! re and im of [k_{in}r j_{n2}(k_{in}r)]'/(k_{in}r)
        !                          *  h_{n1}(k_{out}r):

        b4r = c4r - qdji2 * qy1
        b4i = c4i + qdjr2 * qy1

        v = 1.0_dp / (cbess%mrr**2 + cbess%mri**2)
        drri = cbess%mrr * v * ddri               !re[1/(k_{in}r)]
        drii = -cbess%mri * v * ddri               !im[1/(k_{in}r)]
        ! re and im of [1/(k_{in}r)] j_{n2}(k_{in}r) j_{n1}(k_{out}r):

        c5r = c1r * drri - c1i * drii
        c5i = c1i * drri + c1r * drii
        ! re and im of [1/(k_{in}r)] j_{n2}(k_{in}r) h_{n1}(k_{out}r):

        b5r = b1r * drri - b1i * drii
        b5i = b1i * drri + b1r * drii
        !%%%%%%%  forming integrands of j-matrices (j^{11}=j^{22}=0 for m=0): %%%%%%%%
        !                   uri=dr(i)        !dr/(d\theta)
        uri = -dr_from_x(1)  ! todo: why 'minus' sign was need to fit previous line?
!        rri=rr(i)        !w(i)*r^2(\theta)
        ! w(i)*r^2(\theta)*d2n1*d2n2:
        f1 = r_from_x(1) * a22      !prefactor containing r^2(\theta)<->hat{r} part
        ! n1*(n1+1)*w(i)*r(\theta)*[dr/(d\theta)]*d1n1*d2n2:
        f2 = r_from_x(1) * uri * an1 * a12     !prefactor containing r(\theta)*[dr/(d\theta)]
        !hat{theta} part
        ar12 = f1*b2r + f2*b3r        !~re j^{12}
        ai12 = f1 * b2i + f2 * b3i        !~im j^{12}

        gr12 = f1 * c2r + f2 * c3r        !~re rg j^{12}
        gi12 = f1 * c2i + f2 * c3i        !~im rg j^{12}

        !*  n2*(n2+1)*w(i)*r(\theta)*[dr/(d\theta)]*d2n1*d1n2:
        f2 = r_from_x(1) * uri * an2 * a21     !prefactor containing r(\theta)*[dr/(d\theta)]
        !                                          !hat{theta} part

        ar21 = f1 * b4r + f2 * b5r        !~re j^{21}
        ai21 = f1 * b4i + f2 * b5i        !~im j^{21}

        gr21 = f1 * c4r + f2 * c5r        !~re rg j^{21}
        gi21 = f1 * c4i + f2 * c5i        !~im rg j^{21}

        m0_integrand = (/ ar12, ai12, gr12, gi12, ar21, ai21, gr21, gi21/)        !~re j^{12}
        return
    end function m0_integrand
    !=======================================================================
    function m0_ar12(xi)
        real(dp) m0_ar12, xi, val(8)
        val = m0_integrand(xi)
        m0_ar12 = val(1)
        return
    end function
    !=======================================================================
    function m0_ai12(xi)
        real(dp) m0_ai12, xi, val(8)
        val = m0_integrand(xi)
        m0_ai12 = val(2)
        return
    end function
    !=======================================================================
    function m0_gr12(xi)
        real(dp) m0_gr12, xi, val(8)
        val = m0_integrand(xi)
        m0_gr12 = val(3)
        return
    end function
    !=======================================================================
    function m0_gi12(xi)
        real(dp) m0_gi12, xi, val(8)
        val = m0_integrand(xi)
        m0_gi12 = val(4)
        return
    end function
    !=======================================================================
    function m0_ar21(xi)
        real(dp) m0_ar21, xi, val(8)
        val = m0_integrand(xi)
        m0_ar21 = val(5)
        return
    end function
    !=======================================================================
    function m0_ai21(xi)
        real(dp) m0_ai21, xi, val(8)
        val = m0_integrand(xi)
        m0_ai21 = val(6)
        return
    end function
    !=======================================================================
    function m0_gr21(xi)
        real(dp) m0_gr21, xi, val(8)
        val = m0_integrand(xi)
        m0_gr21 = val(7)
        return
    end function
    !=======================================================================
    function m0_gi21(xi)
        real(dp) m0_gi21, xi, val(8)
        val = m0_integrand(xi)
        m0_gi21 = val(8)
        return
    end function
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
        !                h=rev*( (2_dp/(3_dp*eps*eps))**(1_dp/3_dp) )
        !
        !
        !   ngauss ... the number of gif division points
        !   ng=2*ngauss
        !
        !   1.le.i.le.ngauss
        !
        !--------/---------/---------/---------/---------/---------/---------/--
        implicit none
        real(dp) rev, eps, h, a, si, co, rthet, rad
        integer ng, ngauss, i
        real(dp) x(:), r(:), dr(:)
        rev = mpar%rev
        eps = mpar%eps
        ! determine half-length of the cylinder
        h = rev * ((2.0_dp / (3.0_dp * eps * eps))**(1.0_dp / 3.0_dp))
        ! determine cylinder radius:
        a = h * eps
        ng = size(x)
        ngauss = ng / 2
        if (ng == 1) ngauss = 1
        do i = 1, ngauss
            co = -x(i)
            si = dsqrt(1.0_dp - co * co)

            if ((h * si)>(a * co)) then
                ! along the circular surface:
                rad = a / si
                rthet = -a * co / (si * si)
                !rad=1.d-10
                !rthet=0.0_dp
!                write(*,*) 'cyl'
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
        !  the abscissas and weights of the Gaussian-Legendre n-point quadrature
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
        implicit none
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
