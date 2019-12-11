module lib_adaptive_integration
    use constants
    use special_functions
    use libcylinder

    implicit none

    type, private :: integrand_values
        integer :: n1, n2, nmax, m
    end type integrand_values
    type(integrand_values), public :: integrand

    public cmplx_dp


contains
    !=======================================================================
    subroutine integrate_p (fun, a, b, result)
        !    DQAGP is an adaptive integrator that can handle singularities
        !    of the integrand at user specified points,
        implicit none

        integer ( kind = 4 ), parameter :: limit = 1500
        integer ( kind = 4 ), parameter :: npts = 1
        integer ( kind = 4 ), parameter :: npts2 = 2*npts

        integer ( kind = 4 ), parameter :: leniw = 2*limit + npts2
        integer ( kind = 4 ), parameter :: lenw = leniw*2 - npts2

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
        integer ( kind = 4 ), parameter :: lenw = limit*4

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
    function mi_integrand(xi)
        implicit none
        integer n1, n2, m
        real(dp) mi_integrand(16)
        real(dp) x_to_rsp(1), r_from_x(1), dr_from_x(1)
        real(dp) xi, d1n1, d2n1, d1n2, d2n2, a11, a12, a21, a22, qj1, qdj1, &
                qy1, qdy1, qjr2, qji2, qdjr2, qdji2, c1r, c1i, b1r, b1i, c2r, c2i, &
                b2r, b2i,ddri, c3r, c3i, b3r, b3i, c4r, c4i, b4r, b4i, &
                drri, drii, c5r, c5i, b5r, b5i, rr1, f1, f2, an1, an2, uri, &
                ar12, ai12, gr12, gi12, ar21, ai21, gr21, gi21, v, &
                si, rri, ar11, ai11, gr11, gi11, ar22, ai22, gr22, gi22, &
                e1, e2, e3, dssi, dsi, aa1, aa2, b6i, b6r, b7i, b7r, b8i, b8r, &
                c6i, c6r, c7i, c7r, c8i, c8r


        n1 = integrand%n1
        n2 = integrand%n2
        m = integrand%m
        si = (-1.0_dp)**(n1 + n2 + 1)
        an1 = dble(n1*(n1 + 1))
        an2 = dble(n2*(n2 + 1))
        if (mpar%np.ne.-2) stop 'adaptive integration is only implemented for cylinder'
        x_to_rsp(1) = xi
        call rsp_cylinder(x_to_rsp, r_from_x, dr_from_x)
        call vig_1v (xi, n1, m, d1n1, d2n1)
        call vig_1v (xi, n2, m, d1n2, d2n2)
        a11 = d1n1*d1n2
        a12 = d1n1*d2n2
        a21 = d2n1*d1n2
        a22 = d2n1*d2n2

        dssi = dble(m**2)/ (1d0 - xi*xi)       !=dble(m)**2/(\sin^2\theta)
        aa1 = a12 + a21            != d1n1*d2n2+d2n1*d1n2
        aa2 = a11*dssi + a22   !=(d1n1*d1n2)*dble(m)**2/(\sin^2\theta)
        !                          ! +d2n1*d2n2
        ! vector spherical harmonics:
        !  since refractive index is allowed to be complex in general,
        !  the bessel function j_l(k_in*r) is complex. the code below
        !  performs a separation of the complex integrand in waterman's
        !  surface integral into its respective real and imaginary
        !  parts:

        call cbessjdj(r_from_x(1),n1, qj1, qdj1)
        call cbessydy(r_from_x(1),n1, qy1, qdy1)
        ! bessel functions of the interior argument:
        call cbesscjcdj(r_from_x(1),n2, integrand%nmax, qjr2, qji2, qdjr2, qdji2)
        ! re and im of j_{n2}(k_{in}r) j_{n1}(k_{out}r):

        c1r = qjr2*qj1
        c1i = qji2*qj1
        ! re and im of j_{n2}(k_{in}r) h_{n1}(k_{out}r):

        b1r = c1r - qji2*qy1
        b1i = c1i + qjr2*qy1
        ! re and im of j_{n2}(k_{in}r) j_{n1}'(k_{out}r):

        c2r = qjr2*qdj1
        c2i = qji2*qdj1
        ! re and im of j_{n2}(k_{in}r) h_{n1}'(k_{out}r):

        b2r = c2r - qji2*qdy1
        b2i = c2i + qjr2*qdy1

        ddri=1.0_dp/(dsqrt(r_from_x(1))*cbess%wv) !1/(k_{out}r)
        ! re and im of [1/(k_{out}r)]*j_{n2}(k_{in}r) j_{n1}(k_{out}r)

        c3r = ddri*c1r
        c3i = ddri*c1i
        ! re and im of [1/(k_{out}r)]*j_{n2}(k_{in}r) h_{n1}(k_{out}r):

        b3r = ddri*b1r
        b3i = ddri*b1i
        ! re and im of j_{n2}'(k_{in}r) j_{n1}(k_{out}r):

        c4r = qdjr2*qj1
        c4i = qdji2*qj1
        ! re and im of j_{n2}'(k_{in}r) h_{n1}(k_{out}r):

        b4r = c4r - qdji2*qy1
        b4i = c4i + qdjr2*qy1
        ! re and im of [1/(k_{in}r)] j_{n2}(k_{in}r) j_{n1}(k_{out}r):

        v = 1.0_dp/(cbess%mrr**2 + cbess%mri**2)
        drri = cbess%mrr*v*ddri               !re[1/(k_{in}r)]
        drii = -cbess%mri*v*ddri               !im[1/(k_{in}r)]

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
        ! %%%%%%%%%  forming integrands of j-matrices (j^{11}=j^{22}=0 for m=0):

        uri = -dr_from_x(1)  ! todo: why 'minus' sign was need to fit previous line?
        !                uri = dr(i)

        !                wr = w(i)*r_from_x(1)
        dsi = r_from_x(1)*dsqrt(1d0/(1d0 - xi*xi))*dble(m)      !=dble(m)*w(i)*r^2(\theta)/(|\sin\theta|)
        rri = r_from_x(1)

        ! w(i)*r^2(\theta)*[(d1n1*d1n2)*dble(m)**2/(\sin^2\theta)+d2n1*d2n2]:
        f1 = rri*aa2            !prefactor containing r^2(\theta)<->hat{r} part
        ! n1*(n1+1)*w(i)*r(\theta)*[dr/(d\theta)]*d1n1*d2n2:
        f2 = rri*uri*an1*a12     !prefactor containing r(\theta)*[dr/(d\theta)]
        !                                          !hat{theta} part

        ar12 = f1*b2r + f2*b3r        !~re j^{12}
        ai12 = f1*b2i + f2*b3i        !~im j^{12}

        gr12 = f1*c2r + f2*c3r        !~re rg j^{12}
        gi12 = f1*c2i + f2*c3i        !~im rg j^{12}
        ! n2*(n2+1)*w(i)*r(\theta)*[dr/(d\theta)]*d2n1*d1n2:
        f2 = rri*uri*an2*a21     !prefactor containing r(\theta)*[dr/(d\theta)]
        !                                          !hat{theta} part

        ar21 = f1*b4r + f2*b5r
        ai21 = f1*b4i + f2*b5i

        gr21 = f1*c4r + f2*c5r
        gi21 = f1*c4i + f2*c5i
        ! [dble(m)*w(i)*r^2(i)/(|\sin\theta|)]*(d1n1*d2n2+d2n1*d1n2):
        e1 = dsi*aa1

        ar11 = e1*b1r
        ai11 = e1*b1i
        gr11 = e1*c1r
        gi11 = e1*c1i

        e2 = dsi*uri*a11
        e3 = e2*an2
        e2 = e2*an1

        ar22 = e1*b6r + e2*b7r + e3*b8r
        ai22 = e1*b6i + e2*b7i + e3*b8i

        gr22 = e1*c6r + e2*c7r + e3*c8r
        gi22 = e1*c6i + e2*c7i + e3*c8i

        mi_integrand = (/ ar12, ai12, gr12, gi12, ar21, ai21, gr21, gi21, ar11, ai11, gr11, gi11, ar22, ai22, gr22, gi22/)        !~re j^{12}

    end function mi_integrand
    !=======================================================================
    function mi_ar12(xi)
        real(dp) mi_ar12, xi, val(16)
        val = mi_integrand(xi)
        mi_ar12 = val(1)
        return
    end function
    !=======================================================================
    function mi_ai12(xi)
        real(dp) mi_ai12, xi, val(16)
        val = mi_integrand(xi)
        mi_ai12 = val(2)
        return
    end function
    !=======================================================================
    function mi_gr12(xi)
        real(dp) mi_gr12, xi, val(16)
        val = mi_integrand(xi)
        mi_gr12 = val(3)
        return
    end function
    !=======================================================================
    function mi_gi12(xi)
        real(dp) mi_gi12, xi, val(16)
        val = mi_integrand(xi)
        mi_gi12 = val(4)
        return
    end function
    !=======================================================================
    function mi_ar21(xi)
        real(dp) mi_ar21, xi, val(16)
        val = mi_integrand(xi)
        mi_ar21 = val(5)
        return
    end function
    !=======================================================================
    function mi_ai21(xi)
        real(dp) mi_ai21, xi, val(16)
        val = mi_integrand(xi)
        mi_ai21 = val(6)
        return
    end function
    !=======================================================================
    function mi_gr21(xi)
        real(dp) mi_gr21, xi, val(16)
        val = mi_integrand(xi)
        mi_gr21 = val(7)
        return
    end function
    !=======================================================================
    function mi_gi21(xi)
        real(dp) mi_gi21, xi, val(16)
        val = mi_integrand(xi)
        mi_gi21 = val(8)
        return
    end function
    !=======================================================================
    
    function mi_ar11(xi)
        real(dp) mi_ar11, xi, val(16)
        val = mi_integrand(xi)
        mi_ar11 = val(9)
        return
    end function
    !=======================================================================
    function mi_ai11(xi)
        real(dp) mi_ai11, xi, val(16)
        val = mi_integrand(xi)
        mi_ai11 = val(10)
        return
    end function
    !=======================================================================
    function mi_gr11(xi)
        real(dp) mi_gr11, xi, val(16)
        val = mi_integrand(xi)
        mi_gr11 = val(11)
        return
    end function
    !=======================================================================
    function mi_gi11(xi)
        real(dp) mi_gi11, xi, val(16)
        val = mi_integrand(xi)
        mi_gi11 = val(12)
        return
    end function
    !=======================================================================
    function mi_ar22(xi)
        real(dp) mi_ar22, xi, val(16)
        val = mi_integrand(xi)
        mi_ar22 = val(13)
        return
    end function
    !=======================================================================
    function mi_ai22(xi)
        real(dp) mi_ai22, xi, val(16)
        val = mi_integrand(xi)
        mi_ai22 = val(14)
        return
    end function
    !=======================================================================
    function mi_gr22(xi)
        real(dp) mi_gr22, xi, val(16)
        val = mi_integrand(xi)
        mi_gr22 = val(15)
        return
    end function
    !=======================================================================
    function mi_gi22(xi)
        real(dp) mi_gi22, xi, val(16)
        val = mi_integrand(xi)
        mi_gi22 = val(16)
        return
    end function
    !=======================================================================

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
        a12 = d1n1*d2n2
        a21 = d2n1*d1n2
        a22 = d2n1*d2n2

        call cbessjdj(r_from_x(1),n1, qj1, qdj1)
        call cbessydy(r_from_x(1),n1, qy1, qdy1)
        ! bessel functions of the interior argument:
        call cbesscjcdj(r_from_x(1),n2, integrand%nmax, qjr2, qji2, qdjr2, qdji2)
        ! re and im of j_{n2}(k_{in}r) j_{n1}(k_{out}r):

        c1r = qjr2*qj1
        c1i = qji2*qj1
        ! re and im of j_{n2}(k_{in}r) h_{n1}(k_{out}r):
        b1r = c1r - qji2*qy1
        b1i = c1i + qjr2*qy1
        ! re and im of j_{n2}(k_{in}r) [k_{out}r j_{n1}(k_{out}r)]'/(k_{out}r):
        c2r = qjr2*qdj1
        c2i = qji2*qdj1
        ! re and im of j_{n2}(k_{in}r) [k_{out}r h_{n1}(k_{out}r)]'/(k_{out}r):
        b2r = c2r - qji2*qdy1
        b2i = c2i + qjr2*qdy1

        ddri=1.0_dp/(dsqrt(r_from_x(1))*cbess%wv) !1/(k_{out}r)
        ! re and im of [1/(k_{out}r)]*j_{n2}(k_{in}r) j_{n1}(k_{out}r)

        c3r = ddri*c1r
        c3i = ddri*c1i
        ! re and im of [1/(k_{out}r)]*j_{n2}(k_{in}r) h_{n1}(k_{out}r):

        b3r = ddri*b1r
        b3i = ddri*b1i
        ! re and im of  [k_{in}r j_{n2}(k_{in}r)]'/(k_{in}r)
        !                         *j_{n1}(k_{out}r):

        c4r = qdjr2*qj1
        c4i = qdji2*qj1
        ! re and im of [k_{in}r j_{n2}(k_{in}r)]'/(k_{in}r)
        !                         * h_{n1}(k_{out}r):

        b4r = c4r - qdji2*qy1
        b4i = c4i + qdjr2*qy1

        v = 1.0_dp/(cbess%mrr**2 + cbess%mri**2)
        drri = cbess%mrr*v*ddri               !re[1/(k_{in}r)]
        drii = -cbess%mri*v*ddri               !im[1/(k_{in}r)]
        ! re and im of [1/(k_{in}r)] j_{n2}(k_{in}r) j_{n1}(k_{out}r):

        c5r = c1r*drri - c1i*drii
        c5i = c1i*drri + c1r*drii
        ! re and im of [1/(k_{in}r)] j_{n2}(k_{in}r) h_{n1}(k_{out}r):

        b5r = b1r*drri - b1i*drii
        b5i = b1i*drri + b1r*drii
        !%%%%%%%  forming integrands of j-matrices (j^{11}=j^{22}=0 for m=0): %%%%%%%%
        !                   uri=dr(i)        !dr/(d\theta)
        uri = -dr_from_x(1)  ! todo: why 'minus' sign was need to fit previous line?
!        rri=rr(i)        !w(i)*r^2(\theta)
        ! w(i)*r^2(\theta)*d2n1*d2n2:
        f1 = r_from_x(1)*a22      !prefactor containing r^2(\theta)<->hat{r} part
        ! n1*(n1+1)*w(i)*r(\theta)*[dr/(d\theta)]*d1n1*d2n2:
        f2 = r_from_x(1)*uri*an1*a12     !prefactor containing r(\theta)*[dr/(d\theta)]
        !hat{theta} part
        ar12 = f1*b2r + f2*b3r        !~re j^{12}
        ai12 = f1*b2i + f2*b3i        !~im j^{12}

        gr12 = f1*c2r + f2*c3r        !~re rg j^{12}
        gi12 = f1*c2i + f2*c3i        !~im rg j^{12}

        !*  n2*(n2+1)*w(i)*r(\theta)*[dr/(d\theta)]*d2n1*d1n2:
        f2 = r_from_x(1)*uri*an2*a21     !prefactor containing r(\theta)*[dr/(d\theta)]
        !                                          !hat{theta} part

        ar21 = f1*b4r + f2*b5r        !~re j^{21}
        ai21 = f1*b4i + f2*b5i        !~im j^{21}

        gr21 = f1*c4r + f2*c5r        !~re rg j^{21}
        gi21 = f1*c4i + f2*c5i        !~im rg j^{21}

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
end module lib_adaptive_integration
