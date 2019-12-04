module cylinder_blas
    use constants
    use model_parameters

    !use dense_solve
    !use multem_blas
    !use amos
    !use errfun, only : wpop

    implicit none
    ! This subroutines are not used in the main code at the moment,
    ! so they are private in the module and commented out
    ! private solve, zge, zsu, inv1, invert, decomp, prod


contains
    !=======================================================================
    !=======================================================================
    !=======================================================================
    !=======================================================================
    !=======================================================================
    subroutine zsur(a, int, x, emach) !TODO add intent specs
        !     ------------------------------------------------------------------
        !     zsur is  a standard back-substitution  subroutine  using the
        !     output of zge to calculate x times a-inverse, returned in x
        !     ------------------------------------------------------------------
        implicit none
        integer n
        real(dp) emach
        integer    int(:)
        complex(dp) a(:, :), x(:)
        integer    i, ii, in, j, ij
        complex(dp) dum
        integer :: array_shape(2)
        !     ------------------------------------------------------------------
        array_shape = shape(a)
        n = array_shape(1)
        do ii = 2, n
            i = ii - 1
            if((int(i) - i)/=0) then    !if int(i).ne.i switch the i-th and
                                        !int(i)-th elements of the vector x
                in = int(i)
                dum = x(in)
                x(in) = x(i)
                x(i) = dum
            end if
            !
            ! forming a matrix product
            !
            do j = ii, n
                if((abs(a(i, j)) - emach)>0) x(j) = x(j) - x(i) * a(i, j)
            end do
        end do
        !                 !the i-th row of a multiplied by x(i)
        !                            !subtracted from x
        !
        do ii = 1, n
            i = n - ii + 1
            ij = i + 1
            if((i - n)/=0) then
                do j = ij, n
                    x(i) = x(i) - x(j) * a(j, i)
                end do
            end if
            if((abs(a(i, i)) - emach * 1.0d-7)<0)  then
                a(i, i) = emach * (1.0d-7) * (1.d0, 1.d0)
            end if
            x(i) = x(i) / a(i, i)
        end do
        return
    end

    !=======================================================================
    subroutine zger(a, int, emach)
        !     ------------------------------------------------------------------
        !     zge is a standard subroutine to perform gaussian elimination on
        !     a nc*nc matrix 'a' prior  to inversion, details stored in 'int'
        !                   makes an lower diagonal matrix
        !     this routine does not bother about elements directly above
        !     the matrix diagonal as they are not used explicitly in an
        !     accomapanying zse routine
        !     ------------------------------------------------------------------
        integer n
        real(dp), intent(in) :: emach
        integer, intent(out) :: int(:)
        complex(dp), intent(inout) :: a(:, :)
        integer    i, ii, in, j, k
        complex(dp) yr, dum
        integer :: array_shape(2)

        !     ------------------------------------------------------------------
        array_shape = shape(a)
        n = array_shape(1)
        do ii = 2, n
            i = ii - 1
            yr = a(i, i)
            in = i
            !
            ! finding an element with the largest magnitude in the i-th column
            ! below the matrix diagonal (including the diag. element):

            do j = ii, n
                if((abs(yr) - abs(a(i, j)))<0) then
                    yr = a(i, j)
                else
                    exit
                end if
                in = j
            end do

            int(i) = in      !the largest element in the i-th row above the matrix
                             !diagonal is in the in-th row and is denoted by yr

            if ((in - i)/=0) then  !if in.ne.i switch the i-th and in-th columns to the
                do j = i, n     !left of the matrix diagonal (including the diag. element)
                    dum = a(j, i)
                    a(j, i) = a(j, in)
                    a(j, in) = dum
                end do
            end if

            if((abs(yr) - emach)<=0) exit
            !     ! gaussian elemination of matrix elements above the matrix diagonal.
            !     ! subtraction of (a(i,j)/a(i,i)) multiple of the ith column from the jth column
            !     ! in the sub-matrix beginning with the (i+1,i+1) diag. element
            do j = ii, n
                if((abs(a(i, j)) - emach)>0) then
                    a(i, j) = a(i, j) / yr
                    do k = ii, n
                        a(k, j) = a(k, j) - a(i, j) * a(k, i)   !k-t element of j-th column
                    end do
                end if
            end do                   !the elements in the ith row
            !                              !above diagonal are not set to zero
            !
        end do                   !end of "column loop"
        return
    end
    !=======================================================================
    !    subroutine prod(a, b, c, ndim, n)
    !        !--------/---------/---------/---------/---------/---------/---------/--
    !        ! >>> a,b,ndim,n
    !        ! <<< c=a*b
    !        !=================
    !        !--------/---------/---------/---------/---------/---------/---------/--
    !        integer ndim, n, i,j,k
    !        real(dp) a(ndim, n), b(ndim, n), c(ndim, n), cij
    !        !
    !        do i = 1, n
    !            do j = 1, n
    !                cij = 0d0
    !                do k = 1, n
    !                    cij = cij + a(i, k) * b(k, j)
    !                end do
    !                c(i, j) = cij
    !            end do
    !        end do
    !        !
    !        return
    !    end
    !=======================================================================
    !    subroutine inv1 (nmax, f, a)
    !        !--------/---------/---------/---------/---------/---------/---------/--
    !        ! >>> eps
    !        ! <<< rat
    !        !=================
    !        !  nmax - angular momentum cutoff
    !        !--------/---------/---------/---------/---------/---------/---------/--
    !        real(dp), intent(out):: a(:, :) !npn2
    !        real(dp), intent(in) :: f(:, :) !npn2
    !        integer, intent(in) :: nmax
    !        integer npn1, npn2, i, i1, i2, j, j1, j2, ndim, nn1, nn2, nnmax
    !        real(dp) cond
    !        !npn1 = npn2/2
    !        real(dp), allocatable :: b(:), &
    !                work(:), q1(:, :), q2(:, :), &
    !                p1(:, :), p2(:, :)
    !        integer, allocatable :: ipvt(:), ind1(:), ind2(:)
    !        integer :: array_shape(2)
    !
    !        array_shape = shape(f)
    !        npn2 = array_shape(1)
    !        npn1 = npn2/2
    !        allocate(b(npn1), work(npn1), q1(npn1,npn1), q2(npn1,npn1), &
    !                p1(npn1,npn1), p2(npn1,npn1), ipvt(npn1), ind1(npn1), ind2(npn1))
    !
    !        ndim = npn1
    !        nn1 = (dble(nmax) - 0.1d0) * 0.5d0 + 1d0
    !        nn2 = nmax - nn1
    !        !
    !        do i = 1, nmax
    !            ind1(i) = 2 * i - 1
    !            if(i.gt.nn1) ind1(i) = nmax + 2 * (i - nn1)
    !            ind2(i) = 2 * i
    !            if(i.gt.nn2) ind2(i) = nmax + 2 * (i - nn2) - 1
    !        end do
    !        nnmax = 2 * nmax
    !        !
    !        do i = 1, nmax
    !            i1 = ind1(i)
    !            i2 = ind2(i)
    !            do j = 1, nmax
    !                j1 = ind1(j)
    !                j2 = ind2(j)
    !                q1(j, i) = f(j1, i1)
    !                q2(j, i) = f(j2, i2)
    !            end do
    !        end do
    !        !
    !        call invert(ndim, nmax, q1, p1, cond, ipvt, work, b)
    !        call invert(ndim, nmax, q2, p2, cond, ipvt, work, b)
    !        !
    !        a = 0d0
    !
    !        do i = 1, nmax
    !            i1 = ind1(i)
    !            i2 = ind2(i)
    !            do j = 1, nmax
    !                j1 = ind1(j)
    !                j2 = ind2(j)
    !                a(j1, i1) = p1(j, i)
    !                a(j2, i2) = p2(j, i)
    !            end do
    !        end do
    !        !
    !        deallocate(b, work, q1, q2, p1, p2, ipvt, ind1, ind2)
    !
    !        return
    !    end
    !    !=======================================================================
    !    subroutine invert (ndim, n, a, x, cond, ipvt, work, b)
    !        !--------/---------/---------/---------/---------/---------/---------/--
    !        ! >>> eps
    !        ! <<< rat
    !        !=================
    !        !--------/---------/---------/---------/---------/---------/---------/--
    !        integer n, ndim, i, j
    !        real(dp) cond
    !        real(dp) a(ndim, n), x(ndim, n), work(n), b(n)
    !        integer ipvt(n)
    !        !
    !        call decomp (ndim, n, a, cond, ipvt, work)
    !        !
    !        if (cond + 1d0.eq.cond) print 5, cond
    !        !     if (cond+1d0.eq.cond) stop
    !        5  format(' the matrix is singular for the given numerical accuracy '&
    !                , 'cond = ', d12.6)
    !
    !        do i = 1, n
    !            do j = 1, n
    !                b(j) = 0d0
    !                if (j.eq.i) b(j) = 1d0
    !            end do
    !            !
    !            call solve (a, b, ipvt)
    !            !
    !            do j = 1, n
    !                x(j, i) = b(j)
    !            end do
    !        end do
    !        !
    !        return
    !    end
    !    !=======================================================================
    !    subroutine decomp (ndim, n, a, cond, ipvt, work)
    !        integer n, ndim, i,j,k,kb,km1,kp1, m, nm1
    !        real anorm,t,ek, ynorm, znorm
    !        real(dp) a(ndim, n), cond, work(n)
    !        integer ipvt(n)
    !        !
    !        ipvt(n) = 1
    !        if(n.eq.1) go to 80
    !        nm1 = n - 1
    !        anorm = 0d0
    !        do j = 1, n
    !            t = 0d0
    !            do i = 1, n
    !                t = t + dabs(a(i, j))
    !            end do
    !            if (t.gt.anorm) anorm = t
    !        end do
    !        do k = 1, nm1
    !            kp1 = k + 1
    !            m = k
    !            do i = kp1, n
    !                if (dabs(a(i, k)).gt.dabs(a(m, k))) m = i
    !            end do
    !            ipvt(k) = m
    !            if (m.ne.k) ipvt(n) = -ipvt(n)
    !            t = a(m, k)
    !            a(m, k) = a(k, k)
    !            a(k, k) = t
    !            if (t.eq.0d0) go to 35
    !            do i = kp1, n
    !                a(i, k) = -a(i, k) / t
    !            end do
    !            do j = kp1, n
    !                t = a(m, j)
    !                a(m, j) = a(k, j)
    !                a(k, j) = t
    !                if (t.eq.0d0) go to 30
    !                do i = kp1, n
    !                    a(i, j) = a(i, j) + a(i, k) * t
    !                end do
    !                30     continue
    !            end do
    !            35 continue
    !        end do
    !        do k = 1, n
    !            t = 0d0
    !            if (k.eq.1) go to 45
    !            km1 = k - 1
    !            do i = 1, km1
    !                t = t + a(i, k) * work(i)
    !            end do
    !            45       ek = 1d0
    !            if (t.lt.0d0) ek = -1d0
    !            if (a(k, k).eq.0d0) go to 90
    !            work(k) = -(ek + t) / a(k, k)
    !        end do
    !        do kb = 1, nm1
    !            k = n - kb
    !            t = 0d0
    !            kp1 = k + 1
    !            do i = kp1, n
    !                t = t + a(i, k) * work(k)
    !            end do
    !            work(k) = t
    !            m = ipvt(k)
    !            if (m.eq.k) go to 60
    !            t = work(m)
    !            work(m) = work(k)
    !            work(k) = t
    !            60 continue
    !        end do
    !        ynorm = 0d0
    !        do i = 1, n
    !            ynorm = ynorm + dabs(work(i))
    !        end do
    !        !
    !        call solve (a, work, ipvt)
    !        !
    !        znorm = 0d0
    !        do i = 1, n
    !            znorm = znorm + dabs(work(i))
    !        end do
    !        cond = anorm * znorm / ynorm
    !        if (cond.lt.1d0) cond = 1d0
    !        return
    !        80 cond = 1d0
    !        if (a(1, 1).ne.0d0) return
    !        90 cond = 1d52
    !        !
    !        return
    !    end
    !    !=======================================================================
    !    subroutine solve (a, b, ipvt)
    !        !--------/---------/---------/---------/---------/---------/---------/--
    !        ! >>> eps
    !        ! <<< rat
    !        !=================
    !        !--------/---------/---------/---------/---------/---------/---------/--
    !        integer ndim, n, nm1, k, kp1, i, kb, km1, m
    !        real(dp) t
    !        real(dp), intent(in) :: a(:, :)
    !        real(dp), intent(out) :: b(:)
    !        integer :: array_shape(2)
    !        integer, intent(in) :: ipvt(:)
    !        !
    !        array_shape = shape(a)
    !        ndim = array_shape(1)
    !        n = array_shape(2)
    !        if (n.eq.1) then
    !            b(1) = b(1) / a(1, 1)
    !            return
    !        end if
    !        nm1 = n - 1
    !        do k = 1, nm1
    !            kp1 = k + 1
    !            m = ipvt(k)
    !            t = b(m)
    !            b(m) = b(k)
    !            b(k) = t
    !            do i = kp1, n
    !                b(i) = b(i) + a(i, k) * t
    !            end do
    !        end do
    !        do kb = 1, nm1
    !            km1 = n - kb
    !            k = km1 + 1
    !            b(k) = b(k) / a(k, k)
    !            t = -b(k)
    !            do i = 1, km1
    !                b(i) = b(i) + a(i, k) * t
    !            end do
    !        end do
    !        b(1) = b(1) / a(1, 1)
    !        !
    !        return
    !    end
    !    !=======================================================================
    !    subroutine zge(a, int, n, nc, emach)
    !
    !        !     ------------------------------------------------------------------
    !        !     zge is a standard subroutine to perform gaussian elimination on
    !        !     a nc*nc matrix 'a' prior  to inversion, details stored in 'int'
    !        !     ------------------------------------------------------------------
    !        !
    !        ! ..  scalar arguments  ..
    !        !
    !        integer n, nc
    !        real(dp) emach
    !        !
    !        ! ..  array arguments  ..
    !        !
    !        integer    int(nc)
    !        complex(dp) a(nc, nc)
    !        !
    !        ! ..  local scalars  ..
    !        !
    !        integer    i, ii, in, j, k
    !        complex(dp) yr, dum
    !        !     ------------------------------------------------------------------
    !        !
    !        do ii = 2, n
    !            i = ii - 1
    !            yr = a(i, i)
    !            in = i
    !            do j = ii, n
    !                if(abs(yr) - abs(a(j, i)) >= 0) cycle
    !                yr = a(j, i)
    !                in = j
    !            end do
    !            int(i) = in
    !            if(in - i /= 0) then
    !                do j = i, n
    !                    dum = a(i, j)
    !                    a(i, j) = a(in, j)
    !                    a(in, j) = dum
    !                end do
    !            end if
    !            if(abs(yr) - emach > 0) then
    !                do j = ii, n
    !                    if(abs(a(j, i)) - emach<=0) cycle
    !                    a(j, i) = a(j, i) / yr
    !                    do k = ii, n
    !                        a(j, k) = a(j, k) - a(i, k) * a(j, i)
    !                    end do
    !                end do
    !            end if
    !        end do
    !        return
    !    end subroutine
    !    !=======================================================================
    !    subroutine zsu(a, int, x, n, nc, emach)
    !
    !        !     ------------------------------------------------------------------
    !        !     zsu  is  a standard back-substitution  subroutine  using the
    !        !     output of zge to calculate  a-inverse times x, returned in x
    !        !     ------------------------------------------------------------------
    !        !
    !        ! ..  scalar arguments  ..
    !        !
    !        integer n, nc
    !        real(dp) emach
    !        !
    !        ! ..  array arguments  ..
    !        !
    !        integer    int(nc)
    !        complex(dp) a(nc, nc), x(nc)
    !        !
    !        ! ..  local scalars  ..
    !        !
    !        integer    i, ii, in, j, ij
    !        complex(dp) dum
    !        !     ------------------------------------------------------------------
    !        !
    !        do ii = 2, n
    !            i = ii - 1
    !            if(int(i) - i /= 0) then
    !                in = int(i)
    !                dum = x(in)
    !                x(in) = x(i)
    !                x(i) = dum
    !            end if
    !            do j = ii, n
    !                if(abs(a(j, i)) - emach >0) x(j) = x(j) - a(j, i) * x(i)
    !            end do
    !        end do
    !        do ii = 1, n
    !            i = n - ii + 1
    !            ij = i + 1
    !            if(i - n /= 0) then
    !                do j = ij, n
    !                    x(i) = x(i) - a(i, j) * x(j)
    !                end do
    !            end if
    !            if(abs(a(i, i)) - emach * 1.0d-7 < 0) then
    !                a(i, i) = emach * 1.0d-7 * (1.0_dp, 1.0_dp)
    !            else
    !                x(i) = x(i) / a(i, i)
    !            endif
    !
    !        end do
    !        return
    !    end subroutine
    !
    !    !=======================================================================

end module cylinder_blas
