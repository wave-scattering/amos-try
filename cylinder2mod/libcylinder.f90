module libcylinder
!    use dense_solve
!    use multem_blas
!    use amos
!    use errfun, only : wpop

    implicit none
    integer, parameter, public :: dp = kind(0.0D0)
    complex(dp), parameter, public :: ci = (0.0_dp, 1.0_dp)
    complex(dp), parameter, public :: czero = (0.0_dp, 0.0_dp)
    complex(dp), parameter, public :: cone = (1.0_dp, 0.0_dp)
    complex(dp), parameter, public :: ctwo = (2.0_dp, 0.0_dp)
    real(dp), parameter, public :: pi = 4.0_dp * ATAN(1.0_dp)
    public cmplx_dp

    integer, parameter, public :: NPN1=200, NPNG1=2600, NPNG2=2*NPNG1, NPN2=2*NPN1, &
               NPL=NPN2+1, NPN3=NPN1+1,&
               NPN4=NPN1, NPN5=2*NPN4, NPN6=NPN4+1

contains
    complex(dp) function cmplx_dp(re, im)
        real(dp), intent(in) :: re, im
        cmplx_dp = cmplx(re,im, kind=dp)
    end function cmplx_dp


    function zartan(zs)
!        *--------/---------/---------/---------/---------/---------/---------/--
!        * For a complex argument using that
!        *              atan (z)= -ln((1+iz)/(1-iz))*i/2.d0
!        * Hower, a direct application of this formula often
!        * results in a nonzero imaginary part of  $\arctan z$ even
!        * for a  purely real $z=x$.
!        * Therefore, in order to avoid this, in a numerical
!        * implementation, the complex logarithm is written
!        * explicitly as
!        *
!        *  ln\left(\fr{1+iz}{1-iz}\right)= (ar1,ar2)= log(
!                        *  (1.d0 + abs(z)**2 -2.d0*imag(z))/1.d0+abs(z)**2 +2.d0*imag(z)))/2.d0
!        * +
!        * ci*atan(2.d0*dble(z)/(1-abs(z)**2).
!        *
!        * For a real z=x, |x|<1,  log here is purely imaginary and
!        * equals to i*atan(2x/(1-x**2))\equiv 2*i*atan x
!        *--------/---------/---------/---------/---------/---------/---------/--
        complex(dp) zartan
        complex(dp), intent(in) :: zs
        real(dp) xxs,xas,ar1,ar2

        if (aimag(zs).eq.0.) then
            zartan=cmplx_dp(atan(dble(zs)),0.d0)
            return
        end if
        xas=abs(zs)**2
        ar1=log((1.d0+xas-2.d0*aimag(zs))/(1.d0+xas+2.d0*aimag(zs)))/2.d0
        xxs=dble(zs)
!       * special case:
        if(xas.eq.1.) then
            if(xxs.ge.0.) then
                ar2=pi/2.d0
            else if (xxs.lt.0.) then
                ar2=-pi/2.d0
            end if
            zartan =cmplx_dp(ar2,- ar1)/2.d0
        end if
!       * remaining cases:
        ar2=2.d0*xxs/(1.d0-xas)
        if(xas.lt.1.d0)  then     ! 1st and 4th quadrant
            zartan=cmplx_dp(atan(ar2),- ar1)/2.d0
        else if (xas.gt.1. .and. xxs.ge.0.) then       ! 2nd quadrant
            zartan=cmplx_dp(pi+atan(ar2),- ar1)/2.d0
        else if(xas.gt.1. .and. xxs.lt.0.) then        ! 3rd quadrant
            zartan=cmplx_dp(-pi+atan(ar2),- ar1)/2.d0
        end if
        return
    end function zartan


    subroutine vigf(x, lmax, m, dv1, dv2, ddv1)
        !--------/---------/---------/---------/---------/---------/---------/--
        ! >>> x,lmax,m (only nonnegative)
        ! <<< dv1,dv2,ddv1
        ! =============
        !
        !     x=cos(theta), where theta is the polar angle
        !     lmaxd ... maximal angular momentum cutoff
        !     lmax ... floating  angular momentum cutoff
        !
        ! returns \pi and \tau scattering functions in terms of the
        ! wigner d-functions. algorithm as described in eqs. (31-35)
        !  of ref. \cite{mis39} used. (note however a missing $n$
        !  factor in the 2nd term in the curly bracket in
        !   eq. (35) of ref. \cite{mis39}.)
        !
        !     for a given azimuthal number m.ge.0 returns
        !      the wigner d-functions
        !            dv1(n)=dvig(0,m,n,arccos x) = d_{0m}^{(l)}
        !
        !  \pi scattering function:
        !     ddv1(n)=dvig(0,m,n,arccos x)/sin(arccos x)
        !                              = m*d_{0m}^{(l)}/ sin\theta
        !
        !  \tau scattering function:
        !     dv2(n)=[d/d(arccos x)] dvig(0,m,n,arccos x)
        !                              = d d_{0m}^{(l)}/d\theta
        !
        !     for 1.le.n.le.lmax and 0.le.x.le.1
        !      ddv1 is calculated because (dv1/sin\theta) is singular for
        !             either \beta=0 or \beta=\pi
        !     (for a given m.neq.0, only the m.le.n.le.lmax terms are determined!)
        ! =====
        !
        !     in the present case (eq. (b.28) of ref. \ct{mtl})
        !                       (cf. (4.1.24) of \ct{ed}):
        !
        !             d_{00}^{(l)}(\theta)= p_l(\cos\theta)
        !      d_{0m}^{(l)}(\theta)= (-1)^m \sqrt{(l-m)!/(l+m)!}p_l^m(\cos\theta)
        !     (assuming jackson's p^m_l), where $d^{(l)}_{m0}=(-1)^m d^{(l)}_{0m}$)
        !     (edmonds p^m_l is without (-1)**m prefactor; cf. (4.1.24) therein)
        !
        ! special values:
        !
        !     (rodrigues formula [eq. (2.5.14) of ref. \ct{ed}] then yields
        !                       p_1(x)=x; p_2=(3x^2-1)/2; etc.
        !
        !      p_l^m(x) = (-1)^m (1-x^2)^{m/2} \fr{d^m p_l(x)}{dx} ===>
        !                       p_1^1(\cos\beta)=-\sin\beta
        !     therefore,
        !              d_{00}^{(1)}(\beta)=\cos\beta
        !              d_{01}^{(1)}(\beta)=\sin\beta/\sqrt{2}
        !         d d_{00}^{(1)}(\beta)/d\beta=-\sin\beta
        !         d d_{01}^{(1)}(\beta)/d\beta=\cos\beta/\sqrt{2}
        !
        !     acc. eq. (34) of {mis39}:
        !
        !     a_0=1, a_1=1/\sqrt{2}, a_2=\sqrt{3}/(2*\sqrt{2})
        !
        !     therefore [eq. (32) of {mis39}]:
        !              d_{00}^{(0)}(\beta)=1
        !              d_{01}^{(1)}(\beta)=\sin\beta/\sqrt{2}
        !              d_{02}^{(2)}(\beta)=\sqrt{3}\sin^2\beta/(2*\sqrt{2})
        !     and
        !         d d_{00}^{(0)}(\beta)/d\beta=0
        !         d d_{01}^{(1)}(\beta)/d\beta=\cos\beta/\sqrt{2}
        !         d d_{02}^{(2)}(\beta)/d\beta=\sqrt{3}\sin\beta \cos\beta/\sqrt{2}
        !                                = \sqrt{3}\sin (2\beta) /(2*\sqrt{2})
        ! =====
        !     similar to routine vig, which however only returns
        !            dv1(n)=dvig(0,m,n,arccos x) = d_{0m}^{(l)}
        !
        !     when arccos x is very small, a care has to be exercise to generate
        !     nboth ddv1 and dv2. that part has been made using recurrences of
        !     ref. \ct{tks}
        !--------/---------/---------/---------/---------/---------/---------/--
        implicit none

        integer n, lmax, m, i, i2
        real(dp) a, x, qs, d1, d2, d3, der, dn, dx, qn, qn1, qn2, &
                qnm, qnm1, qmm
        real(dp) ddv1(:), dv1(:), dv2(:)

        !* dv1, ddv1, and dv2 initialization
        do n = 1, lmax
            dv1(n) = 0.d0
            ddv1(n) = 0.d0
            dv2(n) = 0.d0
        end do

        dx = dabs(x)
        a = 1.d0                       !initial a_0
        qs = dsqrt(1d0 - x * x)            !sin\theta
        !***********************************************************************
        !*                      nonzero dv1 initialization
        !*
        if (m.ne.0) go to 20
        !*
        !* ddv1(n)=0.d0 [see (3.33) of {tks}]
        !* d1,d2, and d3 below are the three consequent terms
        !*       d_{0m}^{n-1}, d_{0m}^{n}, and d_{0m}^{n+1} beginning
        !*       with n=m
        !*=============
        !* recurrence initialization following d^l_{00}=p_l
        !*         [see eq. (b.27) of {mtl}]
        !*
        d1 = 1.d0                  !d^0_{00}=p_0=1   (see sec. 6.8 of {nr})
        d2 = x                     !d^1_{00}=p_1=x

        do n = 1, lmax
            qn = dble(n)
            qn1 = dble(n + 1)
            qn2 = dble(2 * n + 1)
            !        * recurrence (31) of ref. {mis39} for d^{n+1}_{00}
            d3 = (qn2 * x * d2 - qn * d1) / qn1
            dv1(n) = d2                     !d^n_{00}
            d1 = d2                         !becomes d^{n-1}_{00} in d3
            d2 = d3                         !becomes d^{n}_{00} in d3
        end do
        !*
        go to 100
        !***********************************************************
        !*                           m\neq 0 part
        !*  d_{0m}^m \equiv a_m*(sin\theta)**m   initialization - (33) and recurrence (34) of ref. {mis39}

        20 continue
        do i = 1, m
            i2 = i * 2
            !          !recurrence (33,34) of ref. {mis39} f
            a = a * dsqrt(dble(i2 - 1) / dble(i2)) * qs
        end do
        !*
        !* recurrence initialization following eqs. (32,33) of ref. {mis39}

        d1 = 0.d0                 !=dv1(m-1); see eq. (32) of ref. {mis39}
        d2 = a                    !=dv1(m);   see eq. (33) of ref. {mis39}
        qmm = dble(m * m)
        !*
        do n = m, lmax
            qn = dble(n)
            qn1 = dble(n + 1)
            qn2 = dble(2 * n + 1)
            qnm = dsqrt(qn * qn - qmm)
            qnm1 = dsqrt(qn1 * qn1 - qmm)
            !          !recurrence (31) of ref. {mis39} for d^{n+1}_{0m}
            d3 = (qn2 * x * d2 - qnm * d1) / qnm1
            dv1(n) = d2                             !d^n_{0m}
            d1 = d2                           !becomes d^{n-1}_{0m} in d3,der
            d2 = d3                                 !becomes d^{n}_{0m} in d3
        end do

        go to 100
        !*                        dv1 initialized
        !*             it remains to determine ddv1 and dv2
        !*********************************************************
        !*  (1-cos\theta) is very small:
        !*
        !   for theta=0 [see eqs. above]:
        !              d_{00}^{(0)}(0)=1
        !              d_{01}^{(1)}(0)=0
        !              d_{02}^{(2)}(\beta)=0
        !     and
        !         d d_{00}^{(0)}(\beta)/d\beta=0
        !         d d_{01}^{(1)}(\beta)/d\beta=1/\sqrt{2}
        !         d d_{02}^{(2)}(\beta)/d\beta=0
        !
        !  see eqs. (4.1-4) of \ct{mis91}:
        !
        !   (m/\sin\theta) d_{0m}^l(0)=(\delta_{m\pm 1}/2) \sqrt{l(l+1)}
        !      d d_{0m}^l(0)/d\beta   =(m\delta_{m\pm 1}/2) \sqrt{l(l+1)}
        !
        !*
        !*  (4.2.1) of \ct{ed}:
        !*   d_{0m}^{(l)}(pi) = (-1)^{l+m} \dt_{0,m}
        !*
        !*  (4.2.3) of \ct{ed}:
        !*   d_{0m}^{(l)}(0) = (-1)^{m} \dt_{0,m} = \dt_{0,m}
        !*=======================================
        !*
        !*  if x^l_m=(m/\sin\theta) d_{0m}^{(l)}, then, according to (3.29) of {tks}:
        !*
        !*  x^{m+1}_{m+1}=\sin\theta \sqrt{\fr{2m+1}{2m+2}}
        !*                           \left(\fr{m+1}{m}\right)x^{m}_{m}
        !*
        !*  according to (3.30) of {tks}:
        !*  x^{m+1}_{m}= -\sqrt{2m+1}\,\cos\theta x^{m}_{m}
        !*
        !* according to (3.31) of {tks}:
        !*  x^{l}_{m}=\fr{1}{\sqrt{l^2-m^2}}\,\left[(2l-1)\cos\theta
        !*          x^{l-1}_{m} - \sqrt{(l-1)^2-m^2}}\,\x^{l-2}_{m} \right]
        !*
        !* initial recurrence values are x^1_1=\sqrt{2}/2 and x^l_0=0
        !***********************************************************************
        !*                   nonzero ddv1/dv2 initialization
        !*
        !*                          m = 0

        100  if (m.eq.0) then     !all ddv1(n)=x^l_0=0; see (3.33) of {tks}:
            !* according to (3.37) of {tks}, dv2(0)=0.d0

            dv2(1) = -qs

            if (lmax.ge.2) dv2(2) = 3 * x * dv2(1)

            if (lmax.lt.3) return
            !*
            do n = 3, lmax           !recurrence (3.36) of {tks},
                dv2(n) = (2 * n - 1) * x * dv2(n - 1) / (n - 1) - n * dv2(n - 2) / (n - 1)
            enddo
            !***********************************************************************
            !*                           m > 0

        else if (m.gt.0) then
            !*
            !* >>> determine x^m_m according to eq. (3.29) of {tks}:

            a = 1.d0 / dsqrt(2.d0)               !x^1_1=a_1

            do i = 1, m - 1
                a = qs * dble(i + 1) * dsqrt(2 * i + 1.d0) * a / (i * dsqrt(2 * i + 2.d0))
            enddo
            !* <<< a is now x^m_m; see (3.29) of {tks}

            ddv1(m) = a
            dv2(m) = x * a                        !see (3.34) of {tks}
            !* >>> determine x^{m+1}_m:

            if (m.eq.lmax)  go to 120
            !     ! der=x^{m+1}_m; see (3.30) of {tks}
            der = x * dsqrt(2 * m + 1.d0) * a
            ddv1(m + 1) = der
            dv2(m + 1) = ((m + 1) * x * der - a * dsqrt(2 * m + 1.d0)) / dble(m)  !(3.35) of {tks}
            !* >>> determine remaining x^{l}_m's

            if ((m + 2).eq.lmax)  go to 120

            do n = m + 2, lmax
                d3 = dsqrt(dble(n)**2 - dble(m)**2)
                ddv1(n) = ((2 * n - 1) * x * ddv1(n - 1) - &
                        dsqrt(dble(n - 1)**2 - dble(m)**2) * ddv1(n - 2)) / d3
                !                                                 !see (3.31) of {tks}
                dv2(n) = (n * x * ddv1(n) - ddv1(n - 1) * d3) / dble(m)  !see (3.35) of {tks}
            enddo

        end if

        120 return
    end
end module libcylinder
