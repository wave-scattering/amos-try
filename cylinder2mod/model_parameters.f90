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

    type, private :: model_parameters_type
        real(dp) :: rsnm, x_min, x_max, rl_min, rl_max, nanorod_cap_hr, &
                alpha, beta, thet, thet0, phi, phi0
        integer :: np, nstep, ndefp
        complex(dp) :: cceps, zeps0
    end type model_parameters_type
    type(model_parameters_type), public :: mpar

    character(999), private        :: ini_config_file_name  !< Name of INI file.
    type(file_ini), private        :: fini       !< INI file handler.

    type, private :: particle_type_values
        integer :: spheroid, cylinder, droplet, sphere_cut_top, sphere_cut_bottom, &
                cone, cone_cut_top, cone_on_cylinder, nanorod
    end type particle_type_values
    type(particle_type_values), public, parameter :: particle_type = particle_type_values( &
            -1, -2, -3, -4, -5, -6, -7, -8, -9)

    public cli_parse, ini_parse
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

    !--------/---------/---------/---------/---------/---------/---------/--
    subroutine ini_parse()
!        character(len = :), allocatable :: items(:, :) !< items pairs.
        integer :: error
        character(len = :), allocatable :: string         !< string option.
        real(dp) :: double, d_re, d_im
        integer :: num
        write(6, *) 'reading config from file: ', ini_config_file_name
        call fini%load(filename = ini_config_file_name)

        allocate(character(999) :: string)
        string = repeat(' ', 999)

        call fini%get(section_name = 'general', option_name = 'particle_type', &
                val = string, error = error)
        mpar%np = 0
        if ((trim(string) .eq. 'cylinder').or.(trim(string)=='-2')) mpar%np = -2
        if ((trim(string) .eq. 'nanorod').or.(trim(string)=='-9')) mpar%np = -9

        call fini%get(section_name = 'nanorod', &
                option_name = 'nanorod_cap_hr', val = double, error = error)
        mpar%nanorod_cap_hr = 1_dp ! default is round cap
        if (error==0) mpar%nanorod_cap_hr = double

        call fini%get(section_name = 'cylinder', option_name = 'rl_min', &
                val = double, error = error)
        mpar%rl_min = -1_dp
        if (error==0) mpar%rl_min = double

        call fini%get(section_name = 'cylinder', option_name = 'rl_max', &
                val = double, error = error)
        mpar%rl_max = -1_dp
        if (error==0) mpar%rl_max = double

        call fini%get(section_name = 'cylinder', option_name = 'rl_steps', &
                val = num, error = error)
        mpar%ndefp = 0
        if (error==0) mpar%ndefp = num

        call fini%get(section_name = 'cylinder', option_name = 'radius', &
                val = double, error = error)
        mpar%rsnm = -1_dp
        if (error==0) mpar%rsnm = double

        call fini%get(section_name = 'cylinder', option_name = 'alpha', &
                val = double, error = error)
        mpar%alpha = 0_dp
        if (error==0) mpar%alpha = double

        call fini%get(section_name = 'cylinder', option_name = 'beta', &
                val = double, error = error)
        mpar%beta = 0_dp
        if (error==0) mpar%beta = double

        call fini%get(section_name = 'general', &
                option_name = 'background_epsilon', &
                val = double, error = error)
        mpar%zeps0 = 1_dp
        if (error==0) mpar%zeps0 = double

        d_re = 1_dp
        call fini%get(section_name = 'cylinder', option_name = 'eps_real', &
                val = d_re, error = error)
        d_im = 0_dp
        call fini%get(section_name = 'cylinder', option_name = 'eps_imag', &
                val = d_im, error = error)
        mpar%cceps = d_re + ci * d_im

        call fini%get(section_name = 'beam', option_name = 'theta0', &
                val = double, error = error)
        mpar%thet0 = 0_dp
        if (error==0) mpar%thet0 = double

        call fini%get(section_name = 'beam', option_name = 'phi0', &
                val = double, error = error)
        mpar%phi0 = 0_dp
        if (error==0) mpar%phi0 = double

        call fini%get(section_name = 'beam', option_name = 'theta', &
                val = double, error = error)
        mpar%thet = 0_dp
        if (error==0) mpar%thet = double

        call fini%get(section_name = 'beam', option_name = 'phi', &
                val = double, error = error)
        mpar%phi = 0_dp
        if (error==0) mpar%phi = double

        call fini%get(section_name = 'beam', option_name = 'x_min', &
                val = double, error = error)
        mpar%x_min = -1_dp
        if (error==0) mpar%x_min = double

        call fini%get(section_name = 'beam', option_name = 'x_max', &
                val = double, error = error)
        mpar%x_max = -1_dp
        if (error==0) mpar%x_max = double

        call fini%get(section_name = 'beam', option_name = 'x_steps', &
                val = num, error = error)
        mpar%nstep = 0
        if (error==0) mpar%nstep = num

    end subroutine ini_parse
    !--------/---------/---------/---------/---------/---------/---------/--


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