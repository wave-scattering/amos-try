module dense_solve
    implicit none
    integer, parameter:: dp=kind(0.d0)
    private
    public zgetrf_wrap
contains
    subroutine zgetrf_wrap ( mat_array, permutation_indices )
    implicit none
    integer, parameter:: dp=kind(0.d0)
    complex(dp) ::mat_array(:,:)
    integer :: permutation_indices(:)
    integer :: info, lda, m, n!, array_rank
    integer :: array_shape(2)

    !allocate(array_shape(1:rank(mat_array)))
    array_shape = shape(mat_array)
    m = array_shape(1)
    n = array_shape(2)
    lda = size(permutation_indices)
    info = 0
    if( m < 0 ) then
        info = -1
    else if( n < 0 ) then
        info = -2
    else if( m /= n ) then
        info = -3
    else if( lda < max( 1, m ) ) then
        info = -4
    end if
    if( info /= 0 ) then
        write(*,*) "Erorr in zgetrf_wrap", -info
        stop 1
        return
    end if

    call zgetrf( m, n, mat_array, lda, permutation_indices, info )

    if( info /= 0 ) then
        write(*,*) "Erorr in zgetrf_wrap", -info
        stop 1
    end if
    !stop
    return
end subroutine zgetrf_wrap
end module



