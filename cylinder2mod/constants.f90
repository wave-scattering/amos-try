module constants

    implicit none
    integer, parameter, public :: dp = kind(0.0D0)
    complex(dp), parameter, public :: ci = (0.0_dp, 1.0_dp)
    complex(dp), parameter, public :: czero = (0.0_dp, 0.0_dp)
    complex(dp), parameter, public :: cone = (1.0_dp, 0.0_dp)
    complex(dp), parameter, public :: ctwo = (2.0_dp, 0.0_dp)
    real(dp), parameter, public :: pi = 4.0_dp * ATAN(1.0_dp)

    ! Moroz
    integer, parameter, public :: NPN1=200, NPNG1=2600, NPNG2=2*NPNG1, NPN2=2*NPN1, &
               NPL=NPN2+1, NPN3=NPN1+1,&
               NPN4=NPN1, NPN5=2*NPN4, NPN6=NPN4+1

!    ! Mishenko
!    integer, parameter, public :: NPN1=100, NPNG1=500, NPNG2=2*NPNG1, NPN2=2*NPN1, &
!               NPL=NPN2+1, NPN3=NPN1+1, &
!               NPN4=NPN1, NPN5=2*NPN4, NPN6=NPN4+1

end module constants
