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

    integer, parameter, public :: NPN1=100, NPNG1=2600, NPNG2=2*NPNG1, NPN2=2*NPN1, &
               NPL=NPN2+1, NPN3=NPN1+1,&
               NPN4=NPN1, NPN5=2*NPN4, NPN6=NPN4+1

contains
    complex(dp) function cmplx_dp(re, im)
        real(dp), intent(in) :: re, im
        cmplx_dp = cmplx(re,im, kind=dp)
    end function cmplx_dp

end module libcylinder
