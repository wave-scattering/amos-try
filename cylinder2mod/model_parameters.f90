module model_parameters
    use constants
    use flap, only : command_line_interface
    use penf, only : I4P
    use finer, only :  file_ini

    !    use dense_solve
!    use multem_blas
!    use amos
!    use errfun, only : wpop

    implicit none

    character(999), public        :: ini_config_file_name  !< Name of INI file.
    type(file_ini), public        :: fini       !< INI file handler.

    type, private :: particle_type_values
        integer :: spheroid, cylinder, droplet, sphere_cut_top, sphere_cut_bottom, &
                cone, cone_cut_top, cone_on_cylinder, nanorod
    end type particle_type_values
    type(particle_type_values), public, parameter :: particle_type = particle_type_values( &
            -1, -2, -3, -4, -5, -6, -7, -8, -9)

    public cli_parse
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

    subroutine cli_parse()
        !       !< build and parse test cli.
        type(command_line_interface) :: cli  !< command line interface.
        integer(i4p) :: error !< error trapping flag.

        call cli%init(progname = 'cylinder2mod', &
                authors = '', &
                help = 'usage: ', &
                examples = ["cylinder2mod -i config.ini"], &
                epilog = new_line('a') // "all done")

        call cli%add(switch = '--ini', &
                switch_ab = '-i', &
                help = 'name of ini file', &
                required = .false., &
                def = 'default.ini', &
                act = 'store')

        call cli%parse(error = error) ; if (error/=0) stop

        call cli%get(switch = '--ini', val = ini_config_file_name)
    endsubroutine cli_parse

end module model_parameters
