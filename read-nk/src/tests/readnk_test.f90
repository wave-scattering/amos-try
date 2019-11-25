! Execute readnk_test and check the result with $python3 plot-inidex.py
program readnk_test
    use libreadnk
    implicit none
    integer :: ios, i, total_points
    real(dp) :: from_WL, to_WL
    real(dp), allocatable :: spectra(:,:), spectra_int(:,:)
    real(dp), allocatable :: WLs(:), y(:)
    character(len=:), allocatable :: filename, prefix
    !=======================================================================
    total_points = 101
    from_WL = 0.3_dp
    to_WL = 0.7_dp
    prefix = '../../../'
!    filename = prefix//'nk-data/Au-Johnson.yml'
!    filename = prefix//'nk-data/Au-Olmon-ev.yml'
!    filename = prefix//'nk-data/Ag-Johnson.yml'
!    filename = prefix//'nk-data/Au-Rakic-BB.yml'
!    filename = prefix//'nk-data/Cr-Johnson.yml'
    filename = prefix//'nk-data/Si-Green-2008.yml'
!    filename = prefix//'nk-data/TiO2-Bodurov.yml'
!    filename = prefix//'nk-data/TiO2-Devore-o.yml'
    call read_refractive_index_from_yaml(filename, ios, spectra)

    open(unit=10, file="index.txt", iostat=ios)
    do i = 1, size(spectra,1),1
        write(10,*) spectra(i,:)
    end do
    close(10)

    call get_refractive_index(filename, total_points, from_WL, to_WL, spectra_int)

    open(unit=10, file="index_int.txt", iostat=ios)
    do i = 1, size(spectra_int,1),1
        write(10,*) spectra_int(i,:)
    end do
    close(10)

    write(*,*) "end"

end program readnk_test
