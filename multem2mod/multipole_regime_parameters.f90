module multipole_regime_parameters
    use penf, only : I4P
    use flap, only : command_line_interface
    use finer, only : file_ini
    use finer, only : file_ini

    implicit none

    type, private :: multipole_parameters
        integer :: s_type, s_ord, s_m, type, order, m_numb
    end type multipole_parameters

    type(multipole_parameters), public :: mrp

    character(999), private :: ini_config_file_name  !< Name of INI file.
    type(file_ini), private :: fini       !< INI file handler.

    public cli_parse, ini_parse


    contains

    subroutine ini_parse()

        integer :: error, num
        write(6, *) 'reading config from file: ', ini_config_file_name
        call fini%load(filename = ini_config_file_name)


        call fini%get(section_name = 'selectors', option_name = 'select_multipole_type', val = num, error = error)
        mrp%s_type = 0
        if (error==0 .and. (num == 0 .or. num == 1)) mrp%s_type = num

        call fini%get(section_name = 'selectors', option_name = 'select_multipole_order', val = num, error = error)
        mrp%s_ord = 0
        if (error==0 .and. (num == 0 .or. num == 1)) mrp%s_ord = num

        call fini%get(section_name = 'selectors', option_name = 'select_m_projections', val = num, error = error)
        mrp%s_m = 0
        if (error==0 .and. (num == 0 .or. num == 1)) mrp%s_m = num

        call fini%get(section_name = 'regime', option_name = 'type', val = num, error = error)
        mrp%type = 0
        if (error==0 .and. (num == 0 .or. num == 1)) mrp%type = num

        call fini%get(section_name = 'regime', option_name = 'order', val = num, error = error)
        mrp%order = 0
        if (error==0) mrp%order = num

        call fini%get(section_name = 'regime', option_name = 'm_projection', val = num, error = error)
        mrp%m_numb = 0
        if (error==0) mrp%m_numb = num

    end subroutine ini_parse

        subroutine cli_parse()
            !       !< build and parse test cli.
            type(command_line_interface) :: cli  !< command line interface.
            integer(i4p) :: error !< error trapping flag.

            call cli%init(progname = 'multem2', &
                    authors = '', &
                    help = 'usage: ', &
                    examples = ["multem2 -i config.ini"], &
                    epilog = new_line('a') // "all done")

            call cli%add(switch = '--ini', &
                    switch_ab = '-i', &
                    help = 'name of ini file', &
                    required = .false., &
                    def = 'multipole_regime_parameters.ini', &
                    act = 'store')

            call cli%parse(error = error) ; if (error/=0) stop

            call cli%get(switch = '--ini', val = ini_config_file_name)
        endsubroutine cli_parse



end module multipole_regime_parameters