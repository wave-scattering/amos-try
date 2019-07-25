module multem_blas
    use dense_solve
    use libmultem2b, only: dp, cmplx_dp
    implicit none
    private
    public comlr2
contains
    !=======================================================================
    !=======================================================================
    subroutine comlr2(nm, n, low, igh, int, hr, hi, wr, wi, zr, zi, ierr)

        !     ------------------------------------------------------------------
        !     this subroutine finds  the  eigenvalues and  eigenvectors  of  a
        !     complex upper hessenberg  matrix by the modified  lr method. the
        !     eigenvectors  of a complex  general matrix  can also be found if
        !     comhes has been used to reduce this general matrix to hessenberg
        !     form.
        !
        !     on input--->
        !        nm      must  be set to the row dimension  of two-dimensional
        !                array  parameters as  declared in the calling program
        !                dimension statement
        !        n       is the order of the matrix
        !        low,igh are integers determined by the  balancing  subroutine
        !                cbal.  if  cbal  has not been used,  set low=1, igh=n
        !        int     contains information on the rows and  columns  inter-
        !                changed in the reduction by comhes,if performed. only
        !                elements low through igh are used.if the eigenvectors
        !                of the hessenberg matrix are desired,set int(j)=j for
        !                these elements
        !        hr,hi   contain the real and imaginary parts, respectively,of
        !                the complex upper hessenberg matrix. their lower tri-
        !                angles  below the subdiagonal contain the multipliers
        !                which   were  used  in  the  reduction by  comhes, if
        !                performed.  if  the  eigenvectors  of  the hessenberg
        !                matrix are desired,these elements must be set to zero
        !
        !      on output--->
        !                the   upper hessenberg portions of hr and hi have been
        !                destroyed, but  the location hr(1,1) contains the norm
        !                of the triangularized matrix,
        !        wr,wi   contain the real and imaginary parts, respectively, of
        !                the   eigenvalues.  if  an  error  exit  is  made, the
        !                eigenvalues should be correct for indices ierr+1,...,n
        !        zr,zi   contain the real and imaginary parts, respectively, of
        !                the eigenvectors.the eigenvectors are unnormalized. if
        !                an error exit is  made, none of the  eigenvectors  has
        !                been found
        !        ierr    is set to  zero for normal return,
        !          j     if the j-th  eigenvalue has not been  determined after
        !                30 iterations.
        !
        !     arithmetic  is  real  except  for the  replacement  of  the algol
        !     procedure cdiv by  complex division and  use of  the  subroutines
        !     csqrt and cmplx in computing complex square roots.
        !     ------------------------------------------------------------------
        !
        ! ..  scalar arguments  ..
        !
        integer nm, n, low, igh, ierr
        !
        ! ..  array arguments  ..
        !
        integer int(igh)
        real(dp) hr(nm, n), hi(nm, n), wr(n), wi(n), zr(nm, n), zi(nm, n)
        !
        ! ..  local scalars  ..
        !
        integer    i, j, k, l, m, en, ii, jj, ll, mm, nn, im1, ip1, its, mp1, enm1, iend
        real(dp)   si, sr, ti, tr, xi, xr, yi, yr, zzi, zzr, norm, machep
        complex(dp) z3
        !     ------------------------------------------------------------------
        !
        !     ********** machep is a machine dependent parameter specifying
        !                the relative precision of floating point arithmetic.
        !
        machep = 2.0_dp**(-47)
        !
        ierr = 0
        !     ********** initialize eigenvector matrix **********
        do i = 1, n
            !
            do j = 1, n
                zr(i, j) = 0.0_dp
                zi(i, j) = 0.0_dp
                if (i == j) zr(i, j) = 1.0_dp
            end do
        end do
        !     ********** form the matrix of accumulated transformations
        !                from the information left by comhes **********
        iend = igh - low - 1
        if (iend <= 0) go to 180
        !     ********** for i=igh-1 step -1 until low+1 do -- **********
        do ii = 1, iend
            i = igh - ii
            ip1 = i + 1
            !
            do k = ip1, igh
                zr(k, i) = hr(k, i - 1)
                zi(k, i) = hi(k, i - 1)
            end do
            !
            j = int(i)
            if (i == j) go to 160
            !
            do k = i, igh
                zr(i, k) = zr(j, k)
                zi(i, k) = zi(j, k)
                zr(j, k) = 0.0_dp
                zi(j, k) = 0.0_dp
            end do
            !
            zr(j, i) = 1.0_dp
            160 continue
        end do
        !     ********** store roots isolated by cbal **********
        180 do i = 1, n
            if (i >= low .and. i <= igh) go to 200
            wr(i) = hr(i, i)
            wi(i) = hi(i, i)
            200 continue
        end do
        !
        en = igh
        tr = 0.0_dp
        ti = 0.0_dp
        !     ********** search for next eigenvalue **********
        220 if (en < low) go to 680
        its = 0
        enm1 = en - 1
        !     ********** look for single small sub-diagonal element
        !                for l=en step -1 until low do -- **********
        240 do ll = low, en
            l = en + low - ll
            if (l == low) go to 300
            if (abs(hr(l, l - 1)) + abs(hi(l, l - 1)) <=&
                    machep * (abs(hr(l - 1, l - 1)) + abs(hi(l - 1, l - 1))&
                            + abs(hr(l, l)) + abs(hi(l, l)))) go to 300
        end do
        !     ********** form shift **********
        300 if (l == en) go to 660
        if (its == 30) go to 1000
        if (its == 10 .or. its == 20) go to 320
        sr = hr(en, en)
        si = hi(en, en)
        xr = hr(enm1, en) * hr(en, enm1) - hi(enm1, en) * hi(en, enm1)
        xi = hr(enm1, en) * hi(en, enm1) + hi(enm1, en) * hr(en, enm1)
        if (xr == 0.0_dp .and. xi == 0.0_dp) go to 340
        yr = (hr(enm1, enm1) - sr) / 2.0_dp
        yi = (hi(enm1, enm1) - si) / 2.0_dp
        z3 = sqrt(cmplx_dp(yr**2 - yi**2 + xr, 2.0_dp * yr * yi + xi))
        zzr = dble(z3)
        zzi = aimag(z3)
        if (yr * zzr + yi * zzi >= 0.0_dp) go to 310
        zzr = -zzr
        zzi = -zzi
        310 z3 = cmplx_dp(xr, xi) / cmplx_dp(yr + zzr, yi + zzi)
        sr = sr - dble(z3)
        si = si - aimag(z3)
        go to 340
        !     ********** form exceptional shift **********
        320 sr = abs(hr(en, enm1)) + abs(hr(enm1, en - 2))
        si = abs(hi(en, enm1)) + abs(hi(enm1, en - 2))
        !
        340 do i = low, en
            hr(i, i) = hr(i, i) - sr
            hi(i, i) = hi(i, i) - si
        end do
        !
        tr = tr + sr
        ti = ti + si
        its = its + 1
        !     ********** look for two consecutive small
        !                sub-diagonal elements **********
        xr = abs(hr(enm1, enm1)) + abs(hi(enm1, enm1))
        yr = abs(hr(en, enm1)) + abs(hi(en, enm1))
        zzr = abs(hr(en, en)) + abs(hi(en, en))
        !     ********** for m=en-1 step -1 until l do -- **********
        do mm = l, enm1
            m = enm1 + l - mm
            if (m == l) go to 420
            yi = yr
            yr = abs(hr(m, m - 1)) + abs(hi(m, m - 1))
            xi = zzr
            zzr = xr
            xr = abs(hr(m - 1, m - 1)) + abs(hi(m - 1, m - 1))
            if (yr <= machep * zzr / yi * (zzr + xr + xi)) go to 420
        end do
        !     ********** triangular decomposition h=l*r **********
        420 mp1 = m + 1
        !
        do i = mp1, en
            im1 = i - 1
            xr = hr(im1, im1)
            xi = hi(im1, im1)
            yr = hr(i, im1)
            yi = hi(i, im1)
            if (abs(xr) + abs(xi) >= abs(yr) + abs(yi)) go to 460
            !     ********** interchange rows of hr and hi **********
            do j = im1, n
                zzr = hr(im1, j)
                hr(im1, j) = hr(i, j)
                hr(i, j) = zzr
                zzi = hi(im1, j)
                hi(im1, j) = hi(i, j)
                hi(i, j) = zzi
            end do
            !
            z3 = cmplx_dp(xr, xi) / cmplx_dp(yr, yi)
            wr(i) = 1.0_dp
            go to 480
            460      z3 = cmplx_dp(yr, yi) / cmplx_dp(xr, xi)
            wr(i) = -1.0_dp
            480      zzr = dble(z3)
            zzi = aimag(z3)
            hr(i, im1) = zzr
            hi(i, im1) = zzi
            !
            do j = i, n
                hr(i, j) = hr(i, j) - zzr * hr(im1, j) + zzi * hi(im1, j)
                hi(i, j) = hi(i, j) - zzr * hi(im1, j) - zzi * hr(im1, j)
            end do
            !
        end do
        !     ********** composition r*l=h **********
        do j = mp1, en
            xr = hr(j, j - 1)
            xi = hi(j, j - 1)
            hr(j, j - 1) = 0.0_dp
            hi(j, j - 1) = 0.0_dp
            !     ********** interchange columns of hr, hi, zr, and zi,
            !                if necessary **********
            if (wr(j) <= 0.0_dp) go to 580
            !
            do i = 1, j
                zzr = hr(i, j - 1)
                hr(i, j - 1) = hr(i, j)
                hr(i, j) = zzr
                zzi = hi(i, j - 1)
                hi(i, j - 1) = hi(i, j)
                hi(i, j) = zzi
            end do
            !
            do i = low, igh
                zzr = zr(i, j - 1)
                zr(i, j - 1) = zr(i, j)
                zr(i, j) = zzr
                zzi = zi(i, j - 1)
                zi(i, j - 1) = zi(i, j)
                zi(i, j) = zzi
            end do
            !
            580   do i = 1, j
                hr(i, j - 1) = hr(i, j - 1) + xr * hr(i, j) - xi * hi(i, j)
                hi(i, j - 1) = hi(i, j - 1) + xr * hi(i, j) + xi * hr(i, j)
            end do
            !     ********** accumulate transformations **********
            do i = low, igh
                zr(i, j - 1) = zr(i, j - 1) + xr * zr(i, j) - xi * zi(i, j)
                zi(i, j - 1) = zi(i, j - 1) + xr * zi(i, j) + xi * zr(i, j)
            end do
            !
        end do
        !
        go to 240
        !     ********** a root found **********
        660 hr(en, en) = hr(en, en) + tr
        wr(en) = hr(en, en)
        hi(en, en) = hi(en, en) + ti
        wi(en) = hi(en, en)
        en = enm1
        go to 220
        !     ********** all roots found.  backsubstitute to find
        !                vectors of upper triangular form **********
        680 norm = 0.0_dp
        !
        do i = 1, n
            !
            do j = i, n
                norm = norm + abs(hr(i, j)) + abs(hi(i, j))
            end do
        end do
        !
        hr(1, 1) = norm
        if (n == 1 .or. norm == 0.0_dp) go to 1001
        !     ********** for en=n step -1 until 2 do -- **********
        do nn = 2, n
            en = n + 2 - nn
            xr = wr(en)
            xi = wi(en)
            enm1 = en - 1
            !     ********** for i=en-1 step -1 until 1 do -- **********
            do ii = 1, enm1
                i = en - ii
                zzr = hr(i, en)
                zzi = hi(i, en)
                if (i == enm1) go to 760
                ip1 = i + 1
                !
                do j = ip1, enm1
                    zzr = zzr + hr(i, j) * hr(j, en) - hi(i, j) * hi(j, en)
                    zzi = zzi + hr(i, j) * hi(j, en) + hi(i, j) * hr(j, en)
                end do
                !
                760           yr = xr - wr(i)
                yi = xi - wi(i)
                if (yr == 0.0_dp .and. yi == 0.0_dp) yr = machep * norm
                z3 = cmplx_dp(zzr, zzi) / cmplx_dp(yr, yi)
                hr(i, en) = dble(z3)
                hi(i, en) = aimag(z3)
            end do
            !
        end do
        !     ********** end backsubstitution **********
        enm1 = n - 1
        !     ********** vectors of isolated roots **********
        do i = 1, enm1
            if (i >= low .and. i <= igh) go to 840
            ip1 = i + 1
            !
            do j = ip1, n
                zr(i, j) = hr(i, j)
                zi(i, j) = hi(i, j)
            end do
            !
            840 continue
        end do
        !     ********** multiply by transformation matrix to give
        !                vectors of original full matrix.
        !                for j=n step -1 until low+1 do -- **********
        do jj = low, enm1
            j = n + low - jj
            m = min0(j - 1, igh)
            !
            do i = low, igh
                zzr = zr(i, j)
                zzi = zi(i, j)
                !
                do k = low, m
                    zzr = zzr + zr(i, k) * hr(k, j) - zi(i, k) * hi(k, j)
                    zzi = zzi + zr(i, k) * hi(k, j) + zi(i, k) * hr(k, j)
                end do
                !
                zr(i, j) = zzr
                zi(i, j) = zzi
            end do
        end do
        !
        go to 1001
        !     ********** set error -- no convergence to an
        !                eigenvalue after 30 iterations **********
        1000 ierr = en
        1001 return
    end subroutine
end module