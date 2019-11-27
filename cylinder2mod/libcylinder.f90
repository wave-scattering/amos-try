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
    subroutine gauleg(x1, x2, x, w, n)
        !--------/---------/---------/---------/---------/---------/---------/--
        !  given the lower and upper limits of integration x1 and x2, and given n
        !  this routine returns arrays x(1:n) and w(1:n) of length n, containing
        !  the abscissas and weights of the gaussian-legendre n-point quadrature
        !  formula.
        !--------/---------/---------/---------/---------/---------/---------/--
        integer, intent(in) :: n
        real(dp), intent(in) :: x1, x2
        real(dp), intent(out) :: x(n), w(n)
        real(dp) eps
        parameter (eps = 3.d-14)
        integer i, j, m
        double precision p1, p2, p3, pp, xl, xm, z, z1

        m = (n + 1) / 2    !the roots are symmetric in the interval, so we only
        xm = 0.5d0 * (x2 + x1)   !have to find half of them
        xl = 0.5d0 * (x2 - x1)

        ! loop over the desired roots:

        do i = 1, m
            z = cos(3.141592654d0 * (i - .25d0) / (n + .5d0))
            ! starting with the above approximation to the ith root, we enter
            ! the main loop of refinement by newton's method.
            1        continue
            p1 = 1.d0
            p2 = 0.d0

            do j = 1, n !loop up the recurrence relation to get legendre
                p3 = p2     !polynomial evaluated at z.
                p2 = p1
                p1 = ((2.d0 * j - 1.d0) * z * p2 - (j - 1.d0) * p3) / j
            end do

            ! p1 is now the desired  legendre polynomial. we next compute pp, its derivative,
            ! by a standard relation involving also p2, the polynomial of one lower order:

            pp = n * (z * p1 - p2) / (z * z - 1.d0)
            z1 = z
            z = z1 - p1 / pp                   !newton's method

            if (abs(z - z1).gt.eps) goto 1
            ! scale the root to the desired interval, and put in its symmetric counterpart:
            x(i) = xm - xl * z
            x(n + 1 - i) = xm + xl * z
            ! compute the weight and its symmetric counterpart:
            w(i) = 2.d0 * xl / ((1.d0 - z * z) * pp * pp)
            w(n + 1 - i) = w(i)

        end do

        return
    end

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
        a = 1d0
        b = 2d0
        c = 3d0
!        data a, b, c /1d0, 2d0, 3d0/
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
