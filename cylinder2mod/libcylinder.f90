module libcylinder
    use constants
    use model_parameters
    use special_functions
!    use dense_solve
!    use multem_blas
!    use amos
!    use errfun, only : wpop

    implicit none
    public cmplx_dp, zge, zsu


contains
    !=======================================================================
    complex(dp) function cmplx_dp(re, im)
        real(dp), intent(in) :: re, im
        cmplx_dp = cmplx(re,im, kind=dp)
    end function cmplx_dp
    !=======================================================================
    !=======================================================================
    !=======================================================================
    !=======================================================================
    !=======================================================================
    !=======================================================================
    !=======================================================================
    !=======================================================================
    !=======================================================================
    !=======================================================================
    !=======================================================================
    !=======================================================================
    subroutine solve (a, b, ipvt)
        !--------/---------/---------/---------/---------/---------/---------/--
        ! >>> eps
        ! <<< rat
        !=================
        !--------/---------/---------/---------/---------/---------/---------/--
        integer ndim, n, nm1, k, kp1, i, kb, km1, m
        real(dp) t
        real(dp), intent(in) :: a(:, :)
        real(dp), intent(out) :: b(:)
        integer :: array_shape(2)
        integer, intent(in) :: ipvt(:)
        !
        array_shape = shape(a)
        ndim = array_shape(1)
        n = array_shape(2)
        if (n.eq.1) then
            b(1) = b(1) / a(1, 1)
            return
        end if
        nm1 = n - 1
        do k = 1, nm1
            kp1 = k + 1
            m = ipvt(k)
            t = b(m)
            b(m) = b(k)
            b(k) = t
            do i = kp1, n
                b(i) = b(i) + a(i, k) * t
            end do
        end do
        do kb = 1, nm1
            km1 = n - kb
            k = km1 + 1
            b(k) = b(k) / a(k, k)
            t = -b(k)
            do i = 1, km1
                b(i) = b(i) + a(i, k) * t
            end do
        end do
        b(1) = b(1) / a(1, 1)
        !
        return
    end
    !=======================================================================
    subroutine zge(a, int, n, nc, emach)

        !     ------------------------------------------------------------------
        !     zge is a standard subroutine to perform gaussian elimination on
        !     a nc*nc matrix 'a' prior  to inversion, details stored in 'int'
        !     ------------------------------------------------------------------
        !
        ! ..  scalar arguments  ..
        !
        integer n, nc
        real(dp) emach
        !
        ! ..  array arguments  ..
        !
        integer    int(nc)
        complex(dp) a(nc, nc)
        !
        ! ..  local scalars  ..
        !
        integer    i, ii, in, j, k
        complex(dp) yr, dum
        !     ------------------------------------------------------------------
        !
        do ii = 2, n
            i = ii - 1
            yr = a(i, i)
            in = i
            do j = ii, n
                if(abs(yr) - abs(a(j, i)) >= 0) cycle
                yr = a(j, i)
                in = j
            end do
            int(i) = in
            if(in - i /= 0) then
                do j = i, n
                    dum = a(i, j)
                    a(i, j) = a(in, j)
                    a(in, j) = dum
                end do
            end if
            if(abs(yr) - emach > 0) then
                do j = ii, n
                    if(abs(a(j, i)) - emach<=0) cycle
                    a(j, i) = a(j, i) / yr
                    do k = ii, n
                        a(j, k) = a(j, k) - a(i, k) * a(j, i)
                    end do
                end do
            end if
        end do
        return
    end subroutine
    !=======================================================================
    subroutine zsu(a, int, x, n, nc, emach)

        !     ------------------------------------------------------------------
        !     zsu  is  a standard back-substitution  subroutine  using the
        !     output of zge to calculate  a-inverse times x, returned in x
        !     ------------------------------------------------------------------
        !
        ! ..  scalar arguments  ..
        !
        integer n, nc
        real(dp) emach
        !
        ! ..  array arguments  ..
        !
        integer    int(nc)
        complex(dp) a(nc, nc), x(nc)
        !
        ! ..  local scalars  ..
        !
        integer    i, ii, in, j, ij
        complex(dp) dum
        !     ------------------------------------------------------------------
        !
        do ii = 2, n
            i = ii - 1
            if(int(i) - i /= 0) then
                in = int(i)
                dum = x(in)
                x(in) = x(i)
                x(i) = dum
            end if
            do j = ii, n
                if(abs(a(j, i)) - emach >0) x(j) = x(j) - a(j, i) * x(i)
            end do
        end do
        do ii = 1, n
            i = n - ii + 1
            ij = i + 1
            if(i - n /= 0) then
                do j = ij, n
                    x(i) = x(i) - a(i, j) * x(j)
                end do
            end if
            if(abs(a(i, i)) - emach * 1.0d-7 < 0) then
                a(i, i) = emach * 1.0d-7 * (1.0_dp, 1.0_dp)
            else
                x(i) = x(i) / a(i, i)
            endif

        end do
        return
    end subroutine

    !=======================================================================

end module libcylinder
