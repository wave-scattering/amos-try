module constants

    implicit none
    integer, parameter, public :: dp = kind(0.0D0)
    complex(dp), parameter, public :: ci = (0.0_dp, 1.0_dp)
    complex(dp), parameter, public :: czero = (0.0_dp, 0.0_dp)
    complex(dp), parameter, public :: cone = (1.0_dp, 0.0_dp)
    complex(dp), parameter, public :: ctwo = (2.0_dp, 0.0_dp)
    real(dp), parameter, public :: pi = 4.0_dp * ATAN(1.0_dp)

    ! Moroz
    integer, parameter, public :: npn1 = 200, npn2 = 2 * npn1, npng1 = 2600, npng2 = 2 * npng1, &
            npl = npn2 + 1, npn3 = npn1 + 1, npn4 = npn1, npn5 = 2 * npn4, npn6 = npn4 + 1

    ! Mishenko
    ! integer, parameter, public :: npn1 = 100, npng1 = 500, npng2 = 2 * npng1, npn2 = 2 * npn1, &
    !        npl = npn2 + 1, npn3 = npn1 + 1, npn4 = npn1, npn5 = 2 * npn4, npn6 = npn4 + 1

end module constants