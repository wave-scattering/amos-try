module libreadnk
    implicit none
    integer, parameter, public :: dp = kind(0.0D0)
    public get_refractive_index, read_refractive_index_from_yaml
contains
    !--------/---------/---------/---------/---------/---------/---------/--
    !--------/---------/---------/---------/---------/---------/---------/--
    subroutine get_refractive_index(filename, total_points, from_WL, to_WL, spectra_int)
        character(*), intent(in) :: filename
        real(dp), intent(in) :: from_WL, to_WL
        integer :: ios, i
        integer, intent(in) :: total_points
        real(dp), allocatable :: spectra(:,:)
        real(dp), allocatable, intent(out) :: spectra_int(:,:)
        real(dp), allocatable :: WLs(:), y(:)
        !============================================================
        call read_refractive_index_from_yaml(filename, ios, spectra)

        allocate(spectra_int(total_points,3))
        allocate(WLs(total_points))
        allocate(y(total_points))

        call linspace(from=from_WL, to=to_WL, array=WLs)
        spectra_int(:,1) = WLs
        y = interpolate(spectra(:,1), spectra(:,2), WLs)
        spectra_int(:,2) = y
        y = interpolate(spectra(:,1), spectra(:,3), WLs)
        spectra_int(:,3) = y

    end subroutine
    !--------/---------/---------/---------/---------/---------/---------/--
    !--------/---------/---------/---------/---------/---------/---------/--
    subroutine read_refractive_index_from_yaml(filename, iostat, spectra)
        integer, intent(out)          :: iostat
        character(*), intent(in) :: filename
        real(dp), allocatable, intent(out) :: spectra(:,:)
        character(len=:), allocatable :: fname
        real(dp) :: wl_point(3)
        real(dp), allocatable :: temp(:,:)
        integer, parameter :: read_unit = 99
        character(len=:), allocatable :: line
        logical :: is_not_in_data_block, is_not_needed_type
        integer :: n, i, key, dsize
        !============================================================

        fname = trim(adjustl(filename))
        !        write(*,*) '==========================='
        !        write(*,*) 'Reading from file', fname
        open(unit=read_unit, file=fname, iostat=iostat)
        if ( iostat /= 0 ) then
            write(*,*) "Error opening data *.yml data file"//fname
            stop 1
        end if

        n = 0
        key = 0
        is_not_in_data_block = .true.
        is_not_needed_type = .true.
        dsize = 0; allocate(spectra(0,3))

        do
            call get_line(read_unit, line, iostat)
            if (iostat /= 0) exit
            if (line == 'DATA:') then
                is_not_in_data_block = .false.
                cycle
            end if
            ! if a new block, which is not a DATA:
            if (len(trim(adjustl(line))) == len(trim(line))) then
                is_not_in_data_block = .true.
                is_not_needed_type = .true.
            end if
            if (is_not_in_data_block ) cycle  ! skip lines in non-data blocks
            if (trim(line) == '    data: |') cycle  ! ignore this line in yaml file
            key = index(line, ':')
            if (key > 0) then
                is_not_needed_type = .true.  ! any data of unexpected type is ignored
                if (trim(line) == '  - type: tabulated nk') is_not_needed_type = .false.
                cycle
            end if
            if (is_not_needed_type) cycle
            ! Now it should be OK to parse data from line ...
            do i = 1,3,1
                line = trim(adjustl(line))//' '
                key = index(line,' ')
                read(line(1:key),*) wl_point(i)
                line(1:key) = ' '
            end do
            ! Save spectra point
            n = n + 1
            allocate(temp(1:n,3))
            if ( allocated (spectra) ) temp (1:n-1,:) = spectra
            call move_alloc (temp, spectra)
            spectra(n,:) = wl_point
        end do
        close(read_unit)
    end subroutine read_refractive_index_from_yaml

    !--------/---------/---------/---------/---------/---------/---------/--
    !--------/---------/---------/---------/---------/---------/---------/--
    subroutine get_line(lun, line, iostat &!, iomsg
            )
        !       From https://stackoverflow.com/a/34752340/4280547
        integer, intent(in)           :: lun
        character(len=:), intent(out), allocatable :: line
        integer, intent(out)          :: iostat
        !        character(*), intent(inout)   :: iomsg
        integer, parameter            :: buffer_len = 80
        character(len=buffer_len)     :: buffer
        integer                       :: size_read
        !============================================================
        line = ''
        do
            read ( lun, '(A)',  &
                    iostat = iostat,  &
                    !                    iomsg = iomsg,  &
                    advance = 'no',  &
                    size = size_read ) buffer
            if (is_iostat_eor(iostat)) then
                line = line // buffer(:size_read)
                iostat = 0
                exit
            else if (iostat == 0) then
                line = line // buffer
            else
                exit
            end if
        end do
    end subroutine get_line

    !--------/---------/---------/---------/---------/---------/---------/--
    !--------/---------/---------/---------/---------/---------/---------/--

    function interpolate(x, y, p) result(r)
        ! From http://fortranwiki.org/fortran/show/interpolation
        !! This function constructs a piecewise cubic Hermitian interpolation of an array y(x) based on discrete numerical data,
        !! and evaluates the interpolation at points p. Note that the mesh spacing of x does not necessarily have to be uniform.
        real(dp), intent(in)  :: x(:)         !! Variable x
        real(dp), intent(in)  :: y(size(x))   !! Function y(x)
        real(dp), intent(in)  :: p(:)         !! Interpolation domain p
        real(dp)              :: r(size(p))   !! Interpolation result y(p)
        real(dp)              :: di(size(p))  !! Interpolation result d(y(p))/dp

        real(dp)              :: d(size(x))
        integer               :: err
        !============================================================
        if (x(1)>p(1) .or. x(size(x))<p(size(p))) stop 'Interpolation domain exceeds available data wavelengths'
        ! Create a PCHIP interpolation of the input data
        !       call dpchez ( nd, xd, yd, dd, spline, wk, nwk, ierr )
        call dpchez(size(x), x, y, d, .false., 0, 0, err)

        ! Extract the interpolated data at provided points
        call dpchev(size(x),  x,  y,  d, size(p),   p,  r, di, err )
        !The original
        !        call dpchfe(size(x),  x,  y,  d, size(p),   p,  r, err)
    end function

    !--------/---------/---------/---------/---------/---------/---------/--
    !--------/---------/---------/---------/---------/---------/---------/--

    subroutine linspace(from, to, array)
        ! From https://stackoverflow.com/a/57211848/4280547
        ! Generates evenly spaced numbers from `from` to `to` (inclusive).
        !
        ! Inputs:
        ! -------
        !
        ! from, to : the lower and upper boundaries of the numbers to generate
        !
        ! Outputs:
        ! -------
        !
        ! array : Array of evenly spaced numbers
        !
        real(dp), intent(in) :: from, to
        real(dp), intent(out) :: array(:)
        real(dp) :: range
        integer :: n, i
        !============================================================
        n = size(array)
        range = to - from

        if (n == 0) return

        if (n == 1) then
            array(1) = from
            return
        end if


        do i=1, n
            array(i) = from + range * (i - 1) / (n - 1)
        end do
    end subroutine

end module libreadnk
