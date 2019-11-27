module libcylinder
    use constants
    use model_parameters
    use special_functions
    use cylinder_blas
!    use dense_solve
!    use multem_blas
!    use amos
!    use errfun, only : wpop

    implicit none
    public cmplx_dp


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
    subroutine gauss (n, ind1, ind2, z, w)
        !--------/---------/---------/---------/---------/---------/---------/--
        ! >>> n,ind1,ind2
        ! <<< z,w
        !=================
        !    calculation of points and weights of gaussian quadrature
        !    formula. if ind1 = 0 - on interval (-1,1), if ind1 = 1 - on
        !    interval  (0,1). if  ind2 = 1 results are printed.
        !
        !    n - number of gif division points (mostly n=ngauss in main program)
        !    z - division points
        !    w - weights
        !--------/---------/---------/---------/---------/---------/---------/--
!        implicit real(dp) (a-h, p-z)
        integer, intent(in) :: n, ind1, ind2
        integer i,k, ind,j,m,niter
        real(dp) x, a,b,c, check, dj,f, pa,pb,pc,zz
        real(dp), intent(out) :: z(:), w(:)
        data a, b, c /1d0, 2d0, 3d0/
        ind = mod(n, 2)
        k = n / 2 + ind
        f = dble(n)
        do 100 i = 1, k
            m = n + 1 - i
            if(i.eq.1) x = a - b / ((f + a) * f)
            if(i.eq.2) x = (z(n) - a) * 4d0 + z(n)
            if(i.eq.3) x = (z(n - 1) - z(n)) * 1.6d0 + z(n - 1)
            if(i.gt.3) x = (z(m + 1) - z(m + 2)) * c + z(m + 3)
            if(i.eq.k.and.ind.eq.1) x = 0d0
            niter = 0
            check = 1d-16
            10     pb = 1d0
            niter = niter + 1
            if (niter.le.100) go to 15
            check = check * 10d0
            15     pc = x
            dj = a
            do j = 2, n
                dj = dj + a
                pa = pb
                pb = pc
                pc = x * pb + (x * pb - pa) * (dj - a) / dj
            end do
            pa = a / ((pb - x * pc) * f)
            pb = pa * pc * (a - x * x)
            x = x - pb
            if(dabs(pb).gt.check * dabs(x)) go to 10
            z(m) = x
            w(m) = pa * pa * (a - x * x)
            if(ind1.eq.0) w(m) = b * w(m)
            if(i.eq.k.and.ind.eq.1) go to 100
            z(i) = -z(m)
            w(i) = w(m)
        100 continue
        if(ind2.ne.1) go to 110
        print 1100, n
        1100 format(' ***  points and weights of gaussian quadrature formula', &
                ' of ', i4, '-th order')
        do i = 1, k
            zz = -z(i)
            print 1200, i, zz, i, w(i)
        end do
        1200 format(' ', 4x, 'x(', i4, ') = ', f17.14, 5x, 'w(', i4, ') = ', f17.14)
        go to 115
        110 continue
        !     print 1300,n
        ! 1300 format(' gaussian quadrature formula of ',i4,'-th order is used')
        115 continue
        if(ind1.eq.0) go to 140
        do i = 1, n
            z(i) = (a + z(i)) / b
        end do
        140 continue

        return
    end
    !=======================================================================

end module libcylinder
