module model_parameters
    use constants
!    use dense_solve
!    use multem_blas
!    use amos
!    use errfun, only : wpop

    implicit none

    type, private :: particle_type_values
        integer :: spheroid, cylinder, droplet, sphere_cut_top, sphere_cut_bottom, &
                cone, cone_cut_top, cone_on_cylinder, nanorod
    end type particle_type_values
    type(particle_type_values), public, parameter :: particle_type = particle_type_values( &
            -1, -2, -3, -4, -5, -6, -7, -8, -9)

!    TYPE Point
!        REAL :: x, y
!    END TYPE Point

    !    TYPE Circle
!        TYPE (point) :: Center
!        REAL :: Radius
!    END TYPE Circle

!    TYPE Circle :: round

!    round%Radius = 10.
!    round%Center%x = 0
!    round%Center%y = 10
    contains

end module model_parameters
