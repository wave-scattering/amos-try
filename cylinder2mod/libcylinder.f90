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

end module libcylinder
