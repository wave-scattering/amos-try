module libmultem2b
    use dense_solve
    use amos
    use errfun, only: erf_pop

    implicit none
    private
    integer, parameter, public:: dp=kind(0.0D0)
    complex(dp), parameter, public :: czero = (0.0_dp, 0.0_dp)
    complex(dp), parameter, public :: ci    = (0.0_dp, 1.0_dp)
    complex(dp), parameter, public :: cone  = (1.0_dp, 0.0_dp)
    complex(dp), parameter, public :: ctwo  = (2.0_dp, 0.0_dp)
    real(dp), parameter, public :: pi=4.0_dp*ATAN(1.0_dp)
    public bessel, tmtrx, sphrm4, ceven, codd, scat, hoslab, blm, elmgen, &
           pair, cmplx_dp, cerf, lat2d, reduce, dlmkg, plw, setup, xmat
contains
    !=======================================================================
    !=======================================================================
    !=======================================================================
    !=======================================================================
    subroutine xmat(xodd, xeven, lmax, kappa, ak, elm, emach, ar1, ar2)
        !     ------------------------------------------------------------------
        !     xmat calculates the matrix describing multiple scatering  within
        !     a  layer, returning  it as :  xodd,  corresponding  to  odd  l+m,
        !     with lm=(10),(2-1),(21),... and xeven, corresponding to even l+m,
        !     with lm=(00),(1-1),(11),(2-2),(20),(22),...
        !     the  program  assumes  that  the  layer is a bravais lattice. the
        !     summation over the lattice follows the ewald method  suggested by
        !     kambe. emach is the machine accuracy.
        !     ------------------------------------------------------------------
        ! ..  parameter statements  ..
        integer   ndend
        parameter (ndend = 1240)
        ! ..  scalar arguments  ..
        integer    lmax
        real(dp)   emach
        complex(dp) kappa
        ! ..  array arguments  ..
        real(dp)  ar1(2), ar2(2)
        real(dp)   ak(2), elm(:)
        complex(dp) xodd(:, :), xeven(:, :)
        ! ..  local scalars  ..
        integer    l2max, ll2, ii, i, nndlm, k, kk, l, mm, nn, m, j1, j2, i1, i2, i3, n1
        integer    na, lll, n, il, nm, in, l2, il2, m2, il3, l3, m3, la1, lb1, la11, lb11
        integer    ll, j, l1
        real(dp)   ab1, ab2, ac, acsq, ad, al, an, an1, an2, ap, ap1, ap2, ar, b
        real(dp)   dnorm, rtpi, rtv, test, test1, test2, tv
        complex(dp) alpha, rta, rtai, kapsq, kant, knsq, xpk, xpa, cf, cp, cx, cz
        complex(dp) z, zz, w, ww, a, acc, gpsq, gp, bt, aa, ab, u, u1, u2, gam
        complex(dp) gk, gkk, sd, alm
        !
        ! ..  local arrays  ..
        !
        real(dp)   denom(ndend), r(2), b1(2), b2(2), akpt(2), fac(4 * lmax + 1)
        complex(dp) gkn(lmax + 1), agk(2 * lmax + 1), xpm(2 * lmax + 1), &
                pref((lmax + 1)**2)
        complex(dp) dlm((lmax + 1) * (2 * lmax + 1))
        !
        ! ..  arrays in common  ..
        !
        !      common/x1/ar1,ar2
        !----------------------------------------------------------------------

        !
        !     ak(1)  and  ak(2)  are the x  and y components of the
        !     momentum parallel to the surface, modulo a reciprocal
        !     lattice vector
        !
        rtpi = sqrt(pi)
        kapsq = kappa * kappa
        !
        !     the factorial  function  is tabulated  in fac . the array
        !     dlm will contain non-zero,i.e. l+m even,values as defined
        !     by kambe.dlm=dlm1+dlm2+dlm3.with lm=(00),(1-1),(11),(2-2)...
        !
        l2max = lmax + lmax
        ll2 = l2max + 1
        fac(1) = 1.0_dp
        ii = l2max + l2max
        do i = 1, ii
            fac(i + 1) = dble(i) * fac(i)
        end do
        nndlm = l2max * (l2max + 3) / 2 + 1
        do i = 1, nndlm
            dlm(i) = czero
        end do
        !
        !     the formula of kambe for the separation constant,alpha,is
        !     used,subject to a restriction which is imposed to control
        !     later rounding errors
        !
        tv = abs(ar1(1) * ar2(2) - ar1(2) * ar2(1))
        alpha = tv / (4.0_dp * pi) * kapsq
        al = abs(alpha)
        if(exp(al) * emach - 5.0d-5 > 0) al = log(5.0d-5 / emach)
        alpha = cmplx_dp(al, 0.0_dp)
        rta = sqrt(alpha)
        !
        !     dlm1 , the  sum  over  reciprocal   lattice  vectors  , is
        !     calculated first. the prefactor p1 is  tabulated  for even
        !     values of l+|m|,thus lm=(00),(11),(2 0),(22),(l2max,l2max)
        !     the  factorial  factor  f1  is simultaneously tabulated in
        !     denom,for all values of n=0,(l-|m|)/2
        !
        k = 1
        kk = 1
        ap1 = -2.0_dp / tv
        ap2 = -1.0_dp
        cf = ci / kappa
        do l = 1, ll2
            ap1 = ap1 / 2.0_dp
            ap2 = ap2 + 2.0_dp
            cp = cf
            mm = 1
            if(mod  (l, 2) == 0) then
                mm = 2
                cp = ci * cp
            end if
            nn = (l - mm) / 2 + 2
            do m = mm, l, 2
                j1 = l + m - 1
                j2 = l - m + 1
                ap = ap1 * sqrt(ap2 * fac(j1) * fac(j2))
                pref(kk) = ap * cp
                cp = -cp
                kk = kk + 1
                nn = nn - 1
                do i = 1, nn
                    i1 = i
                    i2 = nn - i + 1
                    i3 = nn + m - i
                    denom(k) = 1.0_dp / (fac(i1) * fac(i2) * fac(i3))
                    k = k + 1
                end do
            end do
        end do
        !
        !     the  reciprocal  lattice is  defined by  b1,b2 . the  summation
        !     begins with the origin point of the lattice , and  continues in
        !     steps of 8*n1 points , each  step involving the  perimeter of a
        !     parallelogram of lattice points about the origin,of side 2*n1+1
        !     each step begins at label 9.
        !     akpt=the current lattice vector in the sum
        !
        rtv = 2.0_dp * pi / tv
        b1(1) = -ar1(2) * rtv
        b1(2) = ar1(1) * rtv
        b2(1) = -ar2(2) * rtv
        b2(2) = ar2(1) * rtv
        test1 = 1.0d6
        ii = 1
        n1 = -1
        do
            n1 = n1 + 1
            na = n1 + n1 + ii
            an1 = dble(n1)
            an2 = -an1 - 1.0_dp
            do i1 = 1, na
                an2 = an2 + 1.0_dp
                do i2 = 1, 4
                    !     write(16,307) i1,i2
                    ! 307 format(33x,'i1=',i2,' , i2=',i2/33x,12('='))
                    an = an1
                    an1 = -an2
                    an2 = an
                    ab1 = an1 * b1(1) + an2 * b2(1)
                    ab2 = an1 * b1(2) + an2 * b2(2)
                    akpt(1) = ak(1) + ab1
                    akpt(2) = ak(2) + ab2
                    !
                    !     for  every lattice vector of the sum, three short arrays are
                    !     initialised as below. and used as tables:
                    !     xpm(m) contains values of xpk**|m|
                    !     agk(i) contains values of (ac/kappa)**i
                    !     gkn(n) contains values of (gp/kappa)**(2*n-1)*gam(n,z)
                    !     where l=0,l2max;m=-l,l;n=0,(l-|m|)/2;i=l-2*n
                    !     gam is the incomplete gamma function, which is calculated by
                    !     recurrence  from  the value  for n=0, which  in turn can  be
                    !     expressed in terms of the complex error function cerf
                    !     ac=mod(akpt). note special action if ac=0
                    !
                    acsq = akpt(1) * akpt(1) + akpt(2) * akpt(2)
                    gpsq = kapsq - acsq
                    if(abs(gpsq)<emach * emach)   then
                        write(7, 100)
                        100     format(13x, 'fatal error from xmat:'/3x, 'gpsq is too small.'&
                                /3x, 'give a small but nonzero value for "epsilon"'/3x, &
                                'in the data statement of the main program.'&
                                /3x, 'this defines a small imaginary part'&
                                /3x, 'in the frequency or wavelength value.')
                        stop
                    endif
                    ac = sqrt(acsq)
                    gp = sqrt(gpsq)
                    xpk = czero
                    gk = czero
                    gkk = cmplx_dp(1.0_dp, 0.0_dp)
                    if(ac - emach > 0) then
                        xpk = cmplx_dp(akpt(1) / ac, akpt(2) / ac)
                        gk = ac / kappa
                        gkk = gpsq / kapsq
                    end if
                    xpm(1) = cmplx_dp(1.0_dp, 0.0_dp)
                    agk(1) = cmplx_dp(1.0_dp, 0.0_dp)
                    do i = 2, ll2
                        xpm(i) = xpm(i - 1) * xpk
                        agk(i) = agk(i - 1) * gk
                    end do
                    cf = kappa / gp
                    zz = -alpha * gkk
                    cz = sqrt(-zz)
                    z = -ci * cz
                    cx = exp(-zz)
                    gam = rtpi * cerf(cz)
                    gkn(1) = cf * cx * gam
                    bt = z
                    b = 0.5_dp
                    lll = l2max / 2 + 1
                    do i = 2, lll
                        bt = bt / zz
                        b = b - 1.0_dp
                        gam = (gam - bt) / b
                        cf = cf * gkk
                        gkn(i) = cf * cx * gam
                    end do
                    !
                    !     the contribution to the sum dlm1 for a particular
                    !     reciprocal lattice vector is now accumulated into
                    !     the  elements of dlm,note special action if  ac=0
                    !
                    k = 1
                    kk = 1
                    do l = 1, ll2
                        mm = 1
                        if(mod  (l, 2) == 0) mm = 2
                        n = (l * l + mm) / 2
                        nn = (l - mm) / 2 + 2
                        do m = mm, l, 2
                            acc = czero
                            nn = nn - 1
                            il = l
                            do i = 1, nn
                                acc = acc + denom(k) * agk(il) * gkn(i)
                                il = il - 2
                                k = k + 1
                            end do
                            acc = pref(kk) * acc
                            if(ac - 1.0d-6 > 0) dlm(n) = dlm(n) + acc / xpm(m)
                            if((ac - 1.0d-6 <= 0) .or. (m - 1 /= 0)) then
                                nm = n - m + 1
                                dlm(nm) = dlm(nm) + acc * xpm(m)
                            end if
                            kk = kk + 1
                            n = n + 1
                        end do
                    end do
                    if(ii > 0) exit
                end do
                ii = 0
            end do
            !
            !     after each step of the summation a test on the
            !     convergence  of the  elements of  dlm is  made
            test2 = 0.0_dp
            do i = 1, nndlm
                dnorm = abs(dlm(i))
                test2 = test2 + dnorm * dnorm
            end do
            test = abs((test2 - test1) / test1)
            test1 = test2
            if(test - 0.001d0 <= 0) exit ! TODO: convergence constant
            if(n1 - 10 >= 0) exit
        end do
        if(test - 0.001d0 > 0) then ! TODO: convergence constant
            write(16, 26)n1
            26  format(//13x, 'dlm1,s not converged by n1=', i2)
        else
            write(16, 28)n1
            28    format(//13x, 'dlm1,s converged by n1=', i2)
        end if
        !     write(16,250)dlm
        !250  format(5h0dlm1,//,45(2e13.5,/))
        !
        !     dlm2, the sum over real space lattice vectors, begins with
        !     the adjustment of the array pref, to contain values of the
        !     prefactor  'p2' for lm=(00),(11),(20),(22),...
        !
        kk = 1
        ap1 = tv / (4.0_dp * pi)
        cf = kapsq / ci
        do l = 1, ll2
            cp = cf
            mm = 1
            if(mod  (l, 2) == 0) then
                mm = 2
                cp = -ci * cp
            end if
            j1 = (l - mm) / 2 + 1
            j2 = j1 + mm - 1
            in = j1 + l - 2
            ap2 = ((-1.0_dp)**in) * ap1
            do m = mm, l, 2
                ap = ap2 / (fac(j1) * fac(j2))
                pref(kk) = ap * cp * pref(kk)
                j1 = j1 - 1
                j2 = j2 + 1
                ap2 = -ap2
                cp = -cp
                kk = kk + 1
            end do
        end do
        !
        !     the summation proceeds in steps of 8*n1 lattice points
        !     as before, but this  time excluding  the origin  point
        !     r=the current lattice vector in the sum
        !     ar=mod(r)
        !
        n1 = 0
        do
            n1 = n1 + 1
            na = n1 + n1
            an1 = dble(n1)
            an2 = -an1 - 1.0_dp
            do i1 = 1, na
                an2 = an2 + 1.0_dp
                do i2 = 1, 4
                    an = an1
                    an1 = -an2
                    an2 = an
                    r(1) = an1 * ar1(1) + an2 * ar2(1)
                    r(2) = an1 * ar1(2) + an2 * ar2(2)
                    ar = sqrt(r(1) * r(1) + r(2) * r(2))
                    xpk = cmplx_dp(r(1) / ar, r(2) / ar)
                    xpm(1) = cmplx_dp(1.0_dp, 0.0_dp)
                    do i = 2, ll2
                        xpm(i) = xpm(i - 1) * xpk
                    end do
                    ad = ak(1) * r(1) + ak(2) * r(2)
                    sd = exp(-ad * ci)
                    !
                    !     for each lattice vector the integral 'u' is obtained
                    !     from the recurrence relation in l suggested by kambe
                    !     u1 and u2 are the  initial terms of this recurrence,
                    !     for l#-1 and l=0, and they are evaluated in terms of
                    !     the complex error function cerf
                    !
                    kant = 0.5_dp * ar * kappa
                    knsq = kant * kant
                    z = ci * kant / rta
                    zz = rta - z
                    z = rta + z
                    ww = cerf(-zz)
                    w = cerf(z)
                    aa = 0.5_dp * rtpi * (w - ww) / ci
                    ab = 0.5_dp * rtpi * (w + ww)
                    a = alpha - knsq / alpha
                    xpa = exp(a)
                    u1 = aa * xpa
                    u2 = ab * xpa / kant
                    !
                    !     the contribution to dlm2 from a particular lattice
                    !     vector  is  accumulated into  the elements of  dlm
                    !     this procedure includes the term (kant**l) and the
                    !     recurrence for the integral 'u'
                    !
                    kk = 1
                    al = -0.5_dp
                    cp = rta
                    cf = cmplx_dp(1.0_dp, 0.0_dp)
                    do l = 1, ll2
                        mm = 1
                        if(mod  (l, 2) == 0) mm = 2
                        n = (l * l + mm) / 2
                        do m = mm, l, 2
                            acc = pref(kk) * u2 * cf * sd
                            dlm(n) = dlm(n) + acc / xpm(m)
                            if(m - 1 /= 0) then
                                nm = n - m + 1
                                dlm(nm) = dlm(nm) + acc * xpm(m)
                            end if
                            kk = kk + 1
                            n = n + 1
                        end do
                        al = al + 1.0_dp
                        cp = cp / alpha
                        u = (al * u2 - u1 + cp * xpa) / knsq
                        u1 = u2
                        u2 = u
                        cf = kant * cf
                    end do
                end do
            end do
            !
            !     after each step of the summation a test on the
            !     convergence of the elements of dlm is made
            !
            test2 = 0.0_dp
            do i = 1, nndlm
                dnorm = abs(dlm(i))
                test2 = test2 + dnorm * dnorm
            end do
            test = abs((test2 - test1) / test1)
            test1 = test2
            if(test - 0.001d0 <= 0) exit ! TODO: convergence constant
            if(n1 - 10>=0) exit
        end do
        if(test - 0.001d0 > 0) then ! TODO: convergence constant
            write(16, 44)n1
            44    format(//3x, 'dlm2,s not converged by n1=', i2)
        else
            write(16, 46)n1
            46    format(//3x, 'dlm2,s converged by n1=', i2)
        end if
        !
        !     the term dlm3 has a non-zero contribution  only
        !     when l=m=0.it is evaluated here in terms of the
        !     complex error function cerf
        !
        xpa = exp(-alpha)
        rtai = 1.0_dp / (rtpi * rta)
        acc = kappa * (ci * (xpa - cerf(rta)) - rtai) / xpa
        ap = -0.5_dp / rtpi
        dlm(1) = dlm(1) + ap * acc
        !
        !     finally the elements of dlm are multiplied by the
        !     factor (-1.0_dp)**((m+|m|)/2)
        !
        do l = 2, ll2, 2
            n = l * l / 2 + 1
            do m = 2, l, 2
                dlm(n) = -dlm(n)
                n = n + 1
            end do
        end do
        !     write(16,251) dlm
        ! 251 format(15h0dlm1+dlm2+dlm3,//45(2e13.5,/))
        !
        !     summation over the clebsch-gordon type coefficients
        !     elm proceeds, first for  xodd, and then  for xeven.
        !     this gives the kambe elements  a(l2,m2;l3,m3) which
        !     give the elements  x(l3,m3;l2,m2) of xodd and xeven
        !
        k = 1
        ii = 0
        48  ll = lmax + ii
        i = 1
        do il2 = 1, ll
            l2 = il2 - ii
            m2 = -l2 + 1 - ii
            do i2 = 1, il2
                j = 1
                do il3 = 1, ll
                    l3 = il3 - ii
                    m3 = -l3 + 1 - ii
                    do i3 = 1, il3
                        alm = czero
                        la1 = max0(iabs(l2 - l3), iabs(m2 - m3))
                        lb1 = l2 + l3
                        n = (la1 * (la1 + 2) + m2 - m3 + 2) / 2
                        nn = 2 * la1 + 4
                        lb11 = lb1 + 1
                        la11 = la1 + 1
                        do l1 = la11, lb11, 2
                            alm = alm + elm(k) * dlm(n)
                            n = n + nn
                            nn = nn + 4
                            k = k + 1
                        end do
                        alm = alm / kappa
                        if(i - j == 0) alm = alm + ci
                        if(ii <= 0) then
                            xodd(j, i) = ci * alm
                        else
                            xeven(j, i) = ci * alm
                        end if
                        m3 = m3 + 2
                        j = j + 1
                    end do
                end do
                m2 = m2 + 2
                i = i + 1
            end do
        end do
        if(ii<=0) then
            ii = 1
            goto 48
        endif
        return
    end subroutine
    !=======================================================================
    subroutine setup(lmax, xeven, xodd, te, th, xxmat1, xxmat2)
        !     ------------------------------------------------------------------
        !     this subroutine constructs the secular matrix
        !     ------------------------------------------------------------------
        ! ..  scalar arguments ..
        integer lmax
        ! ..  array arguments ..
        complex(dp) xeven(:, :), xodd(:, :)
        complex(dp) xxmat2(:, :)
        complex(dp) te(:), th(:), xxmat1(:, :)
        ! ..  local scalars ..
        integer ia, la, ma, lmtot, ltt, lmax1, ib, lb, mb, i, lmxod, iaod, iaev, ibod
        integer ibev
        real(dp)  c0, signus, up, c, b1, b2, b3, u1, u2, a, down
        real(dp)  alpha1, alpha2, beta1, beta2
        complex(dp) omega1, omega2, z1, z2, z3
        !     ------------------------------------------------------------------
        lmax1 = lmax + 1
        lmtot = lmax1 * lmax1 - 1
        lmxod = (lmax * lmax1) / 2
        c0 = sqrt(8.0_dp * pi / 3.0_dp)
        signus = 1.0_dp
        iaod = 0
        iaev = lmxod
        do la = 1, lmax
            do ma = -la, la
                if(mod((la + ma), 2)==0) then
                    iaev = iaev + 1
                    ia = iaev
                else
                    iaod = iaod + 1
                    ia = iaod
                end if
                up = dble(2 * la + 1)
                signus = -signus
                c = signus * c0
                b1 = 0.0_dp
                if(abs(ma + 1)<=(la - 1)) b1 = blm(la - 1, ma + 1, 1, -1, la, -ma, lmax)
                b2 = 0.0_dp
                if(abs(ma - 1)<=(la - 1)) b2 = blm(la - 1, ma - 1, 1, 1, la, -ma, lmax)
                u1 = dble((la + ma) * (la - ma))
                u2 = dble((2 * la - 1) * (2 * la + 1))
                b3 = sqrt(u1 / u2)
                alpha1 = sqrt(dble((la - ma) * (la + ma + 1))) / 2.0_dp
                beta1 = sqrt(dble((la + ma) * (la - ma + 1))) / 2.0_dp
                ibod = 0
                ibev = lmxod
                do lb = 1, lmax
                    do mb = -lb, lb
                        if(mod((lb + mb), 2)==0) then
                            ibev = ibev + 1
                            ib = ibev
                        else
                            ibod = ibod + 1
                            ib = ibod
                        end if
                        a = dble(lb * (lb + 1) * la * (la + 1))
                        down = sqrt(a)
                        alpha2 = sqrt(dble((lb - mb) * (lb + mb + 1))) / 2.0_dp
                        beta2 = sqrt(dble((lb + mb) * (lb - mb + 1))) / 2.0_dp
                        ltt = la + ma + lb + mb
                        if(mod(ltt, 2)/=0)           then
                            if(mod((la + ma), 2)==0)       then
                                z1 = ceven(lb, mb + 1, la - 1, ma + 1, xeven)
                                z2 = ceven(lb, mb - 1, la - 1, ma - 1, xeven)
                                z3 = codd (lb, mb, la - 1, ma, xodd)
                                z1 = c * alpha2 * b1 * z1
                                z2 = -c * beta2 * b2 * z2
                                z3 = dble(mb) * b3 * z3
                                omega2 = up * (z1 + z2 + z3) / down
                                xxmat1(ia, ib) = -th(la + 1) * omega2
                                xxmat2(ia, ib) = te(la + 1) * omega2
                            else
                                z1 = codd (lb, mb + 1, la - 1, ma + 1, xodd)
                                z2 = codd (lb, mb - 1, la - 1, ma - 1, xodd)
                                z3 = ceven(lb, mb, la - 1, ma, xeven)
                                z1 = c * alpha2 * b1 * z1
                                z2 = -c * beta2 * b2 * z2
                                z3 = dble(mb) * b3 * z3
                                omega2 = up * (z1 + z2 + z3) / down
                                xxmat1(ia, ib) = te(la + 1) * omega2
                                xxmat2(ia, ib) = -th(la + 1) * omega2
                            end if
                        else
                            if(mod((la + ma), 2)==0)       then
                                z1 = codd (lb, mb - 1, la, ma - 1, xodd)
                                z2 = codd (lb, mb + 1, la, ma + 1, xodd)
                                z3 = ceven(lb, mb, la, ma, xeven)
                                z1 = 2.0_dp * beta1 * beta2 * z1
                                z2 = 2.0_dp * alpha1 * alpha2 * z2
                                z3 = dble(ma) * dble(mb) * z3
                                omega1 = (z1 + z2 + z3) / down
                                xxmat1(ia, ib) = -th(la + 1) * omega1
                                xxmat2(ia, ib) = -te(la + 1) * omega1
                            else
                                z1 = ceven(lb, mb - 1, la, ma - 1, xeven)
                                z2 = ceven(lb, mb + 1, la, ma + 1, xeven)
                                z3 = codd (lb, mb, la, ma, xodd)
                                z1 = 2.0_dp * beta1 * beta2 * z1
                                z2 = 2.0_dp * alpha1 * alpha2 * z2
                                z3 = dble(ma) * dble(mb) * z3
                                omega1 = (z1 + z2 + z3) / down
                                xxmat1(ia, ib) = -te(la + 1) * omega1
                                xxmat2(ia, ib) = -th(la + 1) * omega1
                            end if
                        end if
                    end do
                end do
            end do
        end do
        do i = 1, lmtot
            xxmat1(i, i) = cone + xxmat1(i, i)
            xxmat2(i, i) = cone + xxmat2(i, i)
        end do
        return
    end subroutine
    !=======================================================================
    subroutine plw(kappa, gk, lmax, ae, ah)

        !     ------------------------------------------------------------------
        !     this routine calculates the expansion coefficients 'ae,ah' of an
        !     incident plane electromagnetic wave of wave vector  'kappa' with
        !     components parallel to the surface equal to   '(gk(1),gk(2))'.
        !     ------------------------------------------------------------------
        ! ..  scalar arguments  ..
        !
        integer    lmax
        complex(dp) kappa
        !
        ! ..  array arguments  ..
        !
        complex(dp) ae(:, :), ah(:, :), gk(:)
        !
        ! ..  local scalars  ..
        !
        integer    m, ii, l, i, k
        real(dp)   akpar, fpi, a, signus, akg1, akg2
        complex(dp) ct, st, cf, n1, n2, n3, cc, cc1, z1, z2, z3
        !
        ! ..  local arrays  ..
        !
        complex(dp) ylm((lmax + 1)**2)
        !-----------------------------------------------------------------------
        !
        akg1 = dble(gk(1))
        akg2 = dble(gk(2))
        do k = 1, 2
            ae(k, 1) = czero
            ah(k, 1) = czero
        end do
        fpi = 4.0_dp * pi
        akpar = sqrt(akg1 * akg1 + akg2 * akg2)
        ct = gk(3) / kappa
        st = akpar / kappa
        cf = cone
        if(akpar>1.d-8) cf = cmplx(akg1 / akpar, akg2 / akpar, kind = dp)
        n1 = akg1 / kappa
        n2 = akg2 / kappa
        n3 = gk(3) / kappa
        call sphrm4(ylm, ct, st, cf, lmax)
        ii = 1
        cc = cmplx_dp(fpi, 0.0_dp)
        signus = -1.0_dp
        do l = 1, lmax
            cc = cc * ci
            a = dble(l * (l + 1))
            cc1 = cc / sqrt(a)
            do m = -l, l
                signus = -signus
                ii = ii + 1
                if(abs(m + 1)<=l)  then
                    i = l * l + l - m
                    z1 = cc1 * sqrt(dble((l - m) * (l + m + 1))) * ylm(i) / 2.0_dp
                else
                    z1 = czero
                end if
                if(abs(m - 1)<=l)  then
                    i = l * l + l - m + 2
                    z2 = cc1 * sqrt(dble((l + m) * (l - m + 1))) * ylm(i) / 2.0_dp
                else
                    z2 = czero
                end if
                i = l * l + l - m + 1
                z3 = cc1 * dble(m) * ylm(i)
                ae(1, ii) = signus * ci * (cf * z1 - conjg(cf) * z2)
                ae(2, ii) = -signus * (ct * cf * z1 + st * z3 + ct * conjg(cf) * z2)
                ah(1, ii) = signus * (ct * cf * z1 + st * z3 + ct * conjg(cf) * z2)
                ah(2, ii) = signus * ci * (cf * z1 - conjg(cf) * z2)
            end do
        end do
        return
    end subroutine
    !=======================================================================
    subroutine dlmkg(lmax, a0, gk, signus, kappa, dlme, dlmh, emach)
        !     ------------------------------------------------------------------
        !     this subroutine calculates the coefficients dlm(kg)
        !     ------------------------------------------------------------------
        ! ..  arguments  ..
        integer    lmax
        real(dp)   a0, signus, emach
        complex(dp) kappa
        complex(dp) dlme(2, (lmax + 1)**2), dlmh(2, (lmax + 1)**2), gk(3)
        !  .. local
        integer    k, ii, l, m, i
        real(dp)   akpar, alpha, beta, akg1, akg2
        complex(dp) c0, cc, coef, z1, z2, z3
        complex(dp) ct, st, cf
        complex(dp) ylm((lmax + 1)**2)
        !     ------------------------------------------------------------------
        akg1 = dble(gk(1))
        akg2 = dble(gk(2))
        do k = 1, 2
            dlme(k, 1) = czero
            dlmh(k, 1) = czero
        end do
        if(abs(gk(3))<emach)   then
            write(7, 101)
            stop
        endif
        c0 = 2.0_dp * pi / (kappa * a0 * gk(3) * signus)
        akpar = sqrt(akg1 * akg1 + akg2 * akg2)
        ct = gk(3) / kappa
        st = akpar / kappa
        cf = cone
        if(akpar>1.d-8) cf = cmplx_dp(akg1 / akpar, akg2 / akpar)
        call sphrm4(ylm, ct, st, cf, lmax)
        ii = 1
        cc = cone
        do l = 1, lmax
            cc = cc / ci
            coef = c0 * cc / sqrt(dble(l * (l + 1)))
            do m = -l, l
                ii = ii + 1
                alpha = sqrt(dble((l - m) * (l + m + 1))) / 2.0_dp
                beta = sqrt(dble((l + m) * (l - m + 1))) / 2.0_dp
                if(abs(m + 1)<=l)  then
                    i = l * l + l + m + 2
                    z1 = ylm(i)
                else
                    z1 = czero
                end if
                if(abs(m - 1)<=l)  then
                    i = l * l + l + m
                    z2 = ylm(i)
                else
                    z2 = czero
                end if
                i = l * l + l + m + 1
                z3 = ylm(i)
                dlmh(1, ii) = coef * (beta * ct * cf * z2 - dble(m) * st * z3&
                        + alpha * ct * conjg(cf) * z1)
                dlmh(2, ii) = coef * ci * (beta * cf * z2 - alpha * conjg(cf) * z1)
                dlme(1, ii) = coef * ci * (beta * cf * z2 - alpha * conjg(cf) * z1)
                dlme(2, ii) = -coef * (beta * ct * cf * z2 - dble(m) * st * z3&
                        + alpha * ct * conjg(cf) * z1)
            end do
        end do
        return
        101 format(13x, 'fatal error from dlmkg:'/3x, 'gk(3) is too small.'&
                /3x, 'give a small but nonzero value for "epsilon"'/3x, &
                'in the data statement of the main program.'&
                /3x, 'this defines a small imaginary part'&
                /3x, 'in the frequency or wavelength value.')
    end subroutine
    !=======================================================================
    subroutine reduce(ar1, ar2, ak, igmax, g, ig0, emach)

        !----------------------------------------------------------------------
        !     given the primitive vectors ar1,ar2 of a 2d lattice (in units of
        !     alpha), this  subroutine  reduces a wavector "ak" (in  units  of
        !     2*pi/alpha) within the sbz by adding an appropriate  reciprocal-
        !     lattice vector g(ig0)
        !----------------------------------------------------------------------
        ! ..  scalar arguments  ..
        !
        integer igmax, ig0
        real(dp) emach
        !
        ! ..  array arguments  ..
        !
        real(dp) ar1(2), ar2(2), ak(2), g(:, :)
        !
        ! ..  local scalars  ..
        !
        integer i, j, n, i1, i2
        real(dp) d, b, p, afi, ax, ay, akx, aky, fi0, am, bm, alpha, ra
        !
        ! ..  local arrays ..
        !
        real(dp) vx(6), vy(6), fi(6), x(6), y(6)
        !----------------------------------------------------------------------
        alpha = ar1(1)
        ra = 2.0_dp * pi / alpha
        d = ar2(1)
        b = ar2(2)
        if((dabs(d) - 0.5_dp)>emach) stop 'improper lattice vectors'
        if((dabs(dabs(d) - 0.5_dp)<emach).and.(dabs(dabs(b) - &
                dsqrt(3.0_dp) / 2.0_dp)>emach)) then
            b = 2.0_dp * b           ! centred rectangular lattice
            if((dabs(b) - 1.0_dp)<0.0_dp) then
                vx(1) = 1.0_dp
                vy(1) = 0.5_dp * (1.0_dp / b - b)
                vx(2) = 0.0_dp
                vy(2) = 0.5_dp * (1.0_dp / b + b)
                vx(3) = -1.0_dp
                vy(3) = vy(1)
                vx(4) = -1.0_dp
                vy(4) = -vy(1)
                vx(5) = 0.0_dp
                vy(5) = -vy(2)
                vx(6) = 1.0_dp
                vy(6) = -vy(1)
            else
                vx(1) = 0.5_dp + 0.5_dp / b / b
                vy(1) = 0.0_dp
                vx(2) = 0.5_dp - 0.5_dp / b / b
                vy(2) = 1.0_dp / b
                vx(3) = -vx(2)
                vy(3) = vy(2)
                vx(4) = -vx(1)
                vy(4) = 0.0_dp
                vx(5) = -vx(2)
                vy(5) = -vy(2)
                vx(6) = vx(2)
                vy(6) = -vy(2)
            endif
        else             !oblique or hexagonal lattice
            if(d>0.0_dp) then
                p = 0.5_dp * d * (d - 1) / b / b
                vx(1) = 0.5_dp - p
                vy(1) = 0.5_dp * (1.0_dp - 2.0_dp * d) / b
                vx(2) = 0.5_dp + p
                vy(2) = 0.5_dp / b
            else
                p = 0.5_dp * d * (d + 1) / b / b
                vx(1) = 0.5_dp - p
                vy(1) = -0.5_dp * (1.0_dp + 2.0_dp * d) / b
                vx(2) = 0.5_dp + p
                vy(2) = -0.5_dp / b
            endif
            vx(3) = -vx(2)
            vy(3) = vy(2)
            vx(4) = -vx(1)
            vy(4) = -vy(1)
            vx(5) = -vx(2)
            vy(5) = -vy(2)
            vx(6) = vx(2)
            vy(6) = -vy(2)
        endif
        n = 6
        do i = 1, n
            x(i) = vx(i) * ra
            y(i) = vy(i) * ra
        end do
        if(dabs(d)<emach) then
            n = 4              !rectangular or square lattice
            if(b>0.0_dp) then
                x(1) = vx(6) * ra
                y(1) = vy(6) * ra
                x(2) = vx(4) * ra
                y(2) = vy(4) * ra
                x(4) = vx(1) * ra
                y(4) = vy(1) * ra
            else
                x(2) = vx(3) * ra
                y(2) = vy(3) * ra
                x(3) = vx(4) * ra
                y(3) = vy(4) * ra
                x(4) = vx(6) * ra
                y(4) = vy(6) * ra
            endif
        endif
        !*****vertices are arranged in ascending order of the polar angle fi
        do i = 1, n
            fi(i) = datan2(y(i), x(i))
            if(fi(i)<0.0_dp) fi(i) = fi(i) + 2.0_dp * pi
        end do
        do j = 2, n
            afi = fi(j)
            ax = x(j)
            ay = y(j)
            do i = j - 1, 1, -1
                if(fi(i)<=afi) goto 5
                fi(i + 1) = fi(i)
                x(i + 1) = x(i)
                y(i + 1) = y(i)
            end do
            i = 0
            5    fi(i + 1) = afi
            x(i + 1) = ax
            y(i + 1) = ay
        end do
        !*****"ak" is reduced within the sbz
        ig0 = 1
        6  continue
        akx = ak(1) - g(1, ig0)
        aky = ak(2) - g(2, ig0)
        if((abs(akx)<emach).and.(abs(aky)<emach)) return
        fi0 = datan2(aky, akx)   ! find polar angles of the wavevector
        if(fi0<0.0_dp) fi0 = fi0 + 2.0_dp * pi
        i1 = n
        i = 1
        7       continue
        i2 = i
        if(fi0<fi(i))  go to 8
        i = i + 1
        i1 = i2
        if(i<=n) go to 7
        i1 = n
        i2 = 1
        8       continue
        am = abs(y(i2) * x(i1) - x(i2) * y(i1))
        bm = abs((x(i1) - x(i2)) * aky + (y(i2) - y(i1)) * akx)
        if(am>=bm) then
            ak(1) = akx
            ak(2) = aky
            return
        endif
        ig0 = ig0 + 1
        if(ig0>igmax) stop   'error from reduce:  insufficient nr. of&
                reciprocal lattice vectors '
        goto 6
    end subroutine
    !=======================================================================
    subroutine lat2d(a, b, rmax, imax, id, nta, ntb, vecmod)
        !     --------------------------------------------------------------
        !     given a two dimensional bravais lattice with primitive vectors
        !     (a(1),a(2)) , (b(1),b(2)) , defined so that 'b' is longer than
        !     'a' and their scalar product is positive,this routine calcula-
        !     tes the 'imax' lattice vectors: nta(i) * a + ntb(i) * b,having
        !     length 'vecmod(i)' less than 'rmax'.
        !     --------------------------------------------------------------
        ! ..  arguments ..
        integer imax, id
        real(dp) rmax
        integer nta(id), ntb(id)
        real(dp)a(2), b(2), vecmod(id)
        ! ..  local scalars ..
        integer i, na, nb, na0, j, nma, nmb, iord
        real(dp)rmax2, sp, amod2, bmod2, dum, vmod2, vm
        !     ------------------------------------------------------------------
        rmax2 = rmax * rmax
        !***  check if primitive vectors have positive scalar product
        sp = a(1) * b(1) + a(2) * b(2)
        if(sp<-1.d-06)  then
            b(1) = -b(1)
            b(2) = -b(2)
            sp = -sp
            write(6, 100) a(1), a(2), b(1), b(2)
        end if
        !***  check if 'b' is longer than 'a'
        amod2 = a(1) * a(1) + a(2) * a(2)
        bmod2 = b(1) * b(1) + b(2) * b(2)
        if(bmod2<amod2) then
            write(6, 101)
            do j = 1, 2
                dum = a(j)
                a(j) = b(j)
                b(j) = dum
            end do
            dum = amod2
            amod2 = bmod2
            bmod2 = dum
        endif
        !
        i = 0
        nb = 0
        do while ((nb * nb * bmod2)<=rmax2)
            na = 0
            do
                vmod2 = na * na * amod2 + nb * nb * bmod2 + 2 * na * nb * sp
                if(vmod2>rmax2) exit
                i = i + 1
                if(i>id)  go to 13
                nta(i) = na
                ntb(i) = nb
                vecmod(i) = sqrt(vmod2)
                if(na==0.and.nb==0) then
                    na = na + 1
                    cycle
                end if
                i = i + 1
                if(i>id) go to 13
                nta(i) = -na
                ntb(i) = -nb
                vecmod(i) = sqrt(vmod2)
                na = na + 1
            end do
            nb = nb + 1
        end do
        !
        na0 = sp / amod2 + 1
        nb = 1
        do
            if((nb * nb * (bmod2 - sp * sp / amod2))>rmax2) exit
            na = na0
            do
                vmod2 = na * na * amod2 + nb * nb * bmod2 - 2 * na * nb * sp
                if(vmod2>rmax2) exit
                i = i + 1
                if(i>id)  go to 13
                nta(i) = na
                ntb(i) = -nb
                vecmod(i) = sqrt(vmod2)
                i = i + 1
                if(i>id)  go to 13
                nta(i) = -na
                ntb(i) = nb
                vecmod(i) = sqrt(vmod2)
                na = na + 1
            end do
            na = na0 - 1

            do
                vmod2 = na * na * amod2 + nb * nb * bmod2 - 2 * na * nb * sp
                if(vmod2>rmax2.or.na<=0) exit
                i = i + 1
                if(i>id)  go to 13
                nta(i) = na
                ntb(i) = -nb
                vecmod(i) = sqrt(vmod2)
                i = i + 1
                if(i>id) go to 13
                nta(i) = -na
                ntb(i) = nb
                vecmod(i) = sqrt(vmod2)
                na = na - 1
            end do
            nb = nb + 1
        end do
        imax = i
        !
        do iord = 1, imax
            vm = vecmod(iord)
            do i = imax, iord, -1
                if(vecmod(i)>vm)  cycle
                vm = vecmod(i)
                vecmod(i) = vecmod(iord)
                vecmod(iord) = vm
                nma = nta(i)
                nta(i) = nta(iord)
                nta(iord) = nma
                nmb = ntb(i)
                ntb(i) = ntb(iord)
                ntb(iord) = nmb
            end do
        end do
        !
        return
        13 imax = i - 1
        write(6, 102) imax
        do i = 1, imax
            write(6, 103) i, nta(i), a(1), a(2), ntb(i), b(1), b(2), vecmod(i)
        end do
        stop
        !
        100 format(/13x, 'new primitive vectors defined to have positive scalar&
                product'/13x, 'a=(', 2e14.6, ')'/13x, 'b=(', 2e14.6, ')')
        101 format(/13x, 'w a r n i n g ! !'/'interchange primitive vectors in&
                call lat2d'/)
        102 format(//33x, 'from lat2d: maximum number of neighbours=', i4, &
                '  exceeded'//6x, 'lattice points found (non ordered)')
        103 format(i3, 3x, i5, '*(', 2e14.6, ') +', i5, '*(', 2e14.6, ')', 8x, e14.6)
        !
    end subroutine
    !=======================================================================
    complex(dp) function cerf(z)
        !     cerf,given complex argument z,provides the complex error function:
        !     w(z)=exp(-z**2)*(1.0-erf(-i*z))
        complex(dp), intent(in):: z
        cerf=exp(-z**2)*(1.0-erf_pop(-ci*z))
        return
    end function
    !=======================================================================
    complex(dp) function cmplx_dp(re, im)
        real(dp), intent(in) :: re, im
        cmplx_dp = cmplx(re,im, kind=dp)
    end function cmplx_dp
    !=======================================================================
    subroutine pair(igkmax, igkd, qil, qiil, qiiil, qivl, qir, qiir, qiiir, qivr)

        !     ------------------------------------------------------------------
        !     this subroutine calculates scattering q-matrices for a  double
        !     layer, from the corresponding matrices of the individual, left
        !     (l) and right (r), layers. the results are stored in q*l.
        !     -----------------------------------------------------------------
        ! ..  scalar arguments  ..
        integer igkmax
        ! ..  array arguments  ..
        integer   igkd
        complex(dp) qil (:, :), qiil(:, :), qiiil(:, :)
        complex(dp) qivl(:, :)
        complex(dp) qir (:, :), qiir(:, :), qiiir(:, :)
        complex(dp) qivr(:, :)
        !
        ! ..  local
        integer    igk1, igk2, igk3
        real(dp)   emach
        integer    int(igkd), jnt(igkd)
        complex(dp) qinv1(igkd, igkd), qinv2(igkd, igkd), w1(igkd, igkd)
        complex(dp) w2(igkd, igkd), w3(igkd, igkd), w4(igkd, igkd)
        !
        data emach/1.d-8/
        !-----------------------------------------------------------------------
        !
        do igk1 = 1, igkmax
            do igk2 = 1, igkmax
                qinv1(igk1, igk2) = qil (igk1, igk2)
                qinv2(igk1, igk2) = qivr(igk1, igk2)
                w2(igk1, igk2) = czero
                w3(igk1, igk2) = czero
            end do
        end do

        do igk1 = 1, igkmax
            w2(igk1, igk1) = cone
            w3(igk1, igk1) = cone
            do igk2 = 1, igkmax
                do igk3 = 1, igkmax
                    w2(igk1, igk2) = w2(igk1, igk2)&
                            - qiil (igk1, igk3) * qiiir(igk3, igk2)
                    w3(igk1, igk2) = w3(igk1, igk2)&
                            - qiiir(igk1, igk3) * qiil (igk3, igk2)
                end do
            end do
        end do

        call zgetrf_wrap(w2, int)
        call zgetrf_wrap(w3, jnt)
        do igk2 = 1, igkmax
            call zgetrs_wrap(w2, qinv1(:, igk2), int)
            call zgetrs_wrap(w3, qinv2(:, igk2), jnt)
        end do

        do igk1 = 1, igkmax
            do igk2 = 1, igkmax
                w1(igk1, igk2) = czero
                w2(igk1, igk2) = czero
                w3(igk1, igk2) = czero
                w4(igk1, igk2) = czero
                do igk3 = 1, igkmax
                    w1(igk1, igk2) = w1(igk1, igk2)&
                            + qir  (igk1, igk3) * qinv1(igk3, igk2)
                    w2(igk1, igk2) = w2(igk1, igk2)&
                            + qiil (igk1, igk3) * qinv2(igk3, igk2)
                    w3(igk1, igk2) = w3(igk1, igk2)&
                            + qiiir(igk1, igk3) * qinv1(igk3, igk2)
                    w4(igk1, igk2) = w4(igk1, igk2)&
                            + qivl (igk1, igk3) * qinv2(igk3, igk2)
                end do
            end do
        end do

        do igk1 = 1, igkmax
            do igk2 = 1, igkmax
                qinv1(igk1, igk2) = qiir (igk1, igk2)
                qinv2(igk1, igk2) = qiiil(igk1, igk2)
                do igk3 = 1, igkmax
                    qinv1(igk1, igk2) = qinv1(igk1, igk2)&
                            + qir (igk1, igk3) * w2(igk3, igk2)
                    qinv2(igk1, igk2) = qinv2(igk1, igk2)&
                            + qivl(igk1, igk3) * w3(igk3, igk2)
                end do
            end do
        end do

        do igk1 = 1, igkmax
            do igk2 = 1, igkmax
                qil  (igk1, igk2) = w1   (igk1, igk2)
                qiil (igk1, igk2) = qinv1(igk1, igk2)
                qiiil(igk1, igk2) = qinv2(igk1, igk2)
                qivl (igk1, igk2) = w4   (igk1, igk2)
            end do
        end do

        return
    end subroutine
    !=======================================================================
    subroutine elmgen(elm, nelmd, lmax)

        !     ------------------------------------------------------------------
        !     routine to tabulate the clebsch-gordon type coefficients elm,  for
        !     use with the subroutine xmat. the non-zero elm are tabulated first
        !     for  l2,m2; and l3,m3; odd. then for l2,m2; and l3,m3; even, using
        !     the same scheme as that by which they are accessed in xmat.
        !     ------------------------------------------------------------------
        !
        ! ..  scalar arguments  ..
        !
        integer nelmd, lmax
        !
        ! ..  array arguments  ..
        !
        real(dp) elm(nelmd)
        !
        ! ..  local scalars  ..
        !
        integer k, ii, ll, il2, l2, m2, i2, il3, l3, m3, i3, la1, lb1, la11, lb11, m1
        integer l11, l1, l
        real(dp) fourpi
        !     ------------------------------------------------------------------
        fourpi = 4.0_dp * pi
        k = 1
        ii = 0
        do ii = 0, 1
            ll = lmax + ii
            do il2 = 1, ll
                l2 = il2 - ii
                m2 = -l2 + 1 - ii
                do i2 = 1, il2
                    do il3 = 1, ll
                        l3 = il3 - ii
                        m3 = -l3 + 1 - ii
                        do i3 = 1, il3
                            la1 = max0(iabs(l2 - l3), iabs(m2 - m3))
                            lb1 = l2 + l3
                            la11 = la1 + 1
                            lb11 = lb1 + 1
                            m1 = m2 - m3
                            do l11 = la11, lb11, 2
                                l1 = l11 - 1
                                l = (l2 - l3 - l1) / 2 + m2
                                elm(k) = ((-1.0_dp)**l) * fourpi * blm(l1, m1, l3, m3, l2, -m2, lmax&
                                        )
                                k = k + 1
                            end do
                            m3 = m3 + 2
                        end do
                    end do
                    m2 = m2 + 2
                end do
            end do
        end do
        return
    end subroutine
    !=======================================================================
    real(dp) function blm(l1, m1, l2, m2, l3, m3, lmax)

        !-----------------------------------------------------------------------
        !     function blm  provides  the  integral  of  the  product  of three
        !     spherical harmonics,each of which can be expressed as a prefactor
        !     times  a  legendre  function. the  three  prefactors  are  lumped
        !     together as  factor 'c'; and   the integral of the three legendre
        !     functions follows gaunt summation scheme set out by slater(atomic
        !     structure, vol1, 309,310
        !-----------------------------------------------------------------------
        ! ..  scalar arguments  ..
        !
        integer l1, m1, l2, m2, l3, m3, lmax
        !
        ! ..  local scalars  ..
        !
        integer i, ia1, ia2, ia3, ia4, ia5, ia6, ia7, ia8, ia9, ib1, ib2, ib3, ib4
        integer ib5, ic, ic1, ic2, ic3, ic4, ic5, ic6, is, it, it1, it2, nl1, nl2
        integer nl3, nm1, nm2, nm3, ntemp, nn
        real(dp) sign, a, ad, an, b, bd, bn, c, cd, cn
        !
        ! ..  local arrays  ..
        !
        real(dp)fac(4 * lmax + 2)
        !-----------------------------------------------------------------------
        fac(1) = 1.0_dp
        nn = 4 * lmax + 1
        do i = 1, nn
            fac(i + 1) = dble(i) * fac(i)
        end do

        if(m1 + m2 + m3 /= 0) then
           blm = 0.0_dp
           return
        end if
        if ((l1 - lmax - lmax > 0).or.(l2 - lmax >0).or.(l3 - lmax >0) &
           .or.(l1 - iabs(m1)<0).or.(l2 - iabs(m2)<0).or.(l3 - iabs(m3)<0)) then
            stop 'invalid arguments for blm'
        end if
        if(mod  (l1 + l2 + l3, 2) /= 0) then
            blm = 0.0_dp
            return
        end if
        nl1 = l1
        nl2 = l2
        nl3 = l3
        nm1 = iabs(m1)
        nm2 = iabs(m2)
        nm3 = iabs(m3)
        ic = (nm1 + nm2 + nm3) / 2
        if(max0(nm1, nm2, nm3) - nm1 > 0) then
            if(max0(nm2, nm3) - nm2 > 0) then
                nl1 = l3
                nl3 = l1
                nm1 = nm3
                nm3 = iabs(m1)
            else
                nl1 = l2
                nl2 = l1
                nm1 = nm2
                nm2 = iabs(m1)
            end if
        end if
        if(nl2 - nl3 <0) then
            ntemp = nl2
            nl2 = nl3
            nl3 = ntemp
            ntemp = nm2
            nm2 = nm3
            nm3 = ntemp
        end if
        if(nl3 - iabs(nl2 - nl1) < 0) then
                blm = 0.0_dp
                return
        end if
        !
        !     calculation of factor  'a'.
        !
        is = (nl1 + nl2 + nl3) / 2
        ia1 = is - nl2 - nm3
        ia2 = nl2 + nm2
        ia3 = nl2 - nm2
        ia4 = nl3 + nm3
        ia5 = nl1 + nl2 - nl3
        ia6 = is - nl1
        ia7 = is - nl2
        ia8 = is - nl3
        ia9 = nl1 + nl2 + nl3 + 1
        an = ((-1.0_dp)**ia1) * fac(ia2 + 1) * fac(ia4 + 1) * fac(ia5 + 1) * fac(is + 1)
        ad = fac(ia3 + 1) * fac(ia6 + 1) * fac(ia7 + 1) * fac(ia8 + 1) * fac(ia9 + 1)
        a = an / ad
        !
        !     calculation of sum 'b'
        !
        ib1 = nl1 + nm1
        ib2 = nl2 + nl3 - nm1
        ib3 = nl1 - nm1
        ib4 = nl2 - nl3 + nm1
        ib5 = nl3 - nm3
        it1 = max0(0, -ib4) + 1
        it2 = min0(ib2, ib3, ib5) + 1
        b = 0.0_dp
        sign = (-1.0_dp)**(it1)
        ib1 = ib1 + it1 - 2
        ib2 = ib2 - it1 + 2
        ib3 = ib3 - it1 + 2
        ib4 = ib4 + it1 - 2
        ib5 = ib5 - it1 + 2
        do it = it1, it2
            sign = -sign
            ib1 = ib1 + 1
            ib2 = ib2 - 1
            ib3 = ib3 - 1
            ib4 = ib4 + 1
            ib5 = ib5 - 1
            bn = sign * fac(ib1 + 1) * fac(ib2 + 1)
            bd = fac(it) * fac(ib3 + 1) * fac(ib4 + 1) * fac(ib5 + 1)
            b = b + (bn / bd)
        end do
        !
        !       calculation of factor 'c'
        !
        ic1 = nl1 - nm1
        ic2 = nl1 + nm1
        ic3 = nl2 - nm2
        ic4 = nl2 + nm2
        ic5 = nl3 - nm3
        ic6 = nl3 + nm3
        cn = dble((2 * nl1 + 1) * (2 * nl2 + 1) * (2 * nl3 + 1)) * fac(ic1 + 1) * fac(ic3 + 1) * &
                fac(ic5 + 1)
        cd = fac(ic2 + 1) * fac(ic4 + 1) * fac(ic6 + 1)
        c = cn / (pi * cd)
        c = (sqrt(c)) / 2.0_dp
        blm = ((-1.0_dp)**ic) * a * b * c
        return
    end function
    !=======================================================================
    subroutine hoslab(igmax, kappa1, kappa2, kappa3, ak, g, dl, dr, d, &
            qi, qii, qiii, qiv, emach)
        !-----------------------------------------------------------------------
        !     this subroutine calculates the  q-matrices for a homogeneous
        !     plate  '2' of thickness 'd', having the semi-infinite medium
        !     '1' on its left and the semi-infinite medium '3' on its right
        !     ------------------------------------------------------------------
        !  .. arguments ..
        integer    igmax
        real(dp)   emach, d
        complex(dp) kappa1, kappa2, kappa3
        real(dp)   ak(:), g(:, :), dl(3), dr(3)
        complex(dp) qi(:, :), qii(:, :), qiii(:, :)
        complex(dp) qiv(:, :)
        !  .. local
        integer    i, j, ia, ib, ja, ig1, igkmax
        real(dp)   gkkpar
        complex(dp) gkkz1, gkkz2, gkkz3, z1, z2, z3, cqi, cqii
        complex(dp) cqiii, cqiv, denoma, denomb, gkkdum
        complex(dp) t(4, 2), r(4, 2), x(4), p(4, 2)
        !     -----------------------------------------------------------------
        igkmax = 2 * igmax
        do ia = 1, igkmax
            do ib = 1, igkmax
                qi  (ia, ib) = czero
                qii (ia, ib) = czero
                qiii(ia, ib) = czero
                qiv (ia, ib) = czero
            end do
        end do
        x(1) = kappa1 / kappa2
        x(2) = cone / x(1)
        x(3) = kappa2 / kappa3
        x(4) = cone / x(3)
        do ig1 = 1, igmax
            gkkpar = sqrt((ak(1) + g(1, ig1)) * (ak(1) + g(1, ig1)) + &
                    (ak(2) + g(2, ig1)) * (ak(2) + g(2, ig1)))
            gkkz1 = sqrt(kappa1 * kappa1 - gkkpar * gkkpar)
            gkkz2 = sqrt(kappa2 * kappa2 - gkkpar * gkkpar)
            gkkz3 = sqrt(kappa3 * kappa3 - gkkpar * gkkpar)
            do j = 1, 2
                denoma = x(j) * x(j) * gkkz2 + gkkz1
                denomb = gkkz2 + gkkz1
                if(abs(denoma)<emach.or.abs(denomb)<emach) goto 20
                r(j, 1) = (gkkz1 - x(j) * x(j) * gkkz2) / denoma
                r(j, 2) = (gkkz1 - gkkz2) / denomb
                t(j, 1) = ctwo * x(j) * gkkz1 / denoma
                t(j, 2) = ctwo * gkkz1 / denomb
                gkkdum = gkkz1
                gkkz1 = gkkz2
                gkkz2 = gkkdum
            end do
            ! asdfa
            do j = 3, 4
                denoma = x(j) * x(j) * gkkz3 + gkkz2
                denomb = gkkz3 + gkkz2
                if(abs(denoma)<emach.or.abs(denomb)<emach) goto 20
                r(j, 1) = (gkkz2 - x(j) * x(j) * gkkz3) / denoma
                r(j, 2) = (gkkz2 - gkkz3) / denomb
                t(j, 1) = ctwo * x(j) * gkkz2 / denoma
                t(j, 2) = ctwo * gkkz2 / denomb
                gkkdum = gkkz2
                gkkz2 = gkkz3
                gkkz3 = gkkdum
            end do
            z1 = exp(ci * gkkz2 * d)
            z2 = z1 * z1
            do i = 1, 2
                z3 = cone / (cone - z2 * r(2, i) * r(3, i))
                p(1, i) = t(3, i) * z3 * z1 * t(1, i)
                p(2, i) = r(4, i) + t(4, i) * r(2, i) * t(3, i) * z2 * z3
                p(3, i) = r(1, i) + t(2, i) * r(3, i) * t(1, i) * z2 * z3
                p(4, i) = t(2, i) * z3 * z1 * t(4, i)
            end do
            cqi = exp(ci * ((ak(1) + g(1, ig1)) * (dl(1) + dr(1)) + &
                    (ak(2) + g(2, ig1)) * (dl(2) + dr(2)) + &
                    gkkz1 * dl(3) + gkkz3 * dr(3)))
            cqii = exp(ctwo * ci * gkkz3 * dr(3))
            cqiii = exp(ctwo * ci * gkkz1 * dl(3))
            cqiv = exp(-ci * ((ak(1) + g(1, ig1)) * (dl(1) + dr(1)) + &
                    (ak(2) + g(2, ig1)) * (dl(2) + dr(2)) - &
                    gkkz1 * dl(3) - gkkz3 * dr(3)))
            do ja = 1, 2
                ia = 2 * ig1 - 2 + ja
                qi  (ia, ia) = cqi * p(1, ja)
                qii (ia, ia) = cqii * p(2, ja)
                qiii(ia, ia) = cqiii * p(3, ja)
                qiv (ia, ia) = cqiv * p(4, ja)
            end do
        end do
        return
        20 stop 'fatal error in hoslab'
    end subroutine
    !=======================================================================
    subroutine scat(igmax, zval, ak, g, kapin, kapout, eincid, qi, qiii)
        !     ------------------------------------------------------------------
        !     this subroutine calculates the reflectivity, transmittance and
        !     absorbance of a finite slab, characterized by transmission and
        !     reflection matrices qi and qiii, respectively.
        !     ------------------------------------------------------------------
        ! ..  arguments  ..
        integer    igmax
        real(dp)   zval, kapin, kapout
        real(dp)   ak(:), g(:, :)
        complex(dp) qi(:, :), qiii(:, :), eincid(:)
        ! ..  local
        integer   igkd
        integer    igk1, ig1, k1, igk2, igkmax
        real(dp)   down, refle, trans, absor, gkzin, gkzout, tes1
        complex(dp), allocatable :: etrans(:), erefle(:)
        !     ------------------------------------------------------------------
        igkd = size(eincid)
        allocate(etrans(1:igkd)); allocate(erefle(1:igkd))
        down = 0.0_dp
        refle = 0.0_dp
        trans = 0.0_dp
        igkmax = 2 * igmax
        igk1 = 0
        do ig1 = 1, igmax
            tes1 = (ak(1) + g(1, ig1)) * (ak(1) + g(1, ig1)) + (ak(2) + g(2, ig1)) * (ak(2) + &
                    g(2, ig1))
            gkzin = 0.0_dp
            gkzout = 0.0_dp
            if((kapin * kapin - tes1)>0.0_dp) gkzin = sqrt(kapin * kapin - tes1)
            if((kapout * kapout - tes1)>0.0_dp) gkzout = sqrt(kapout * kapout - tes1)
            do k1 = 1, 2
                igk1 = igk1 + 1
                etrans(igk1) = czero
                erefle(igk1) = czero
                do igk2 = 1, igkmax
                    etrans(igk1) = etrans(igk1) + qi  (igk1, igk2) * eincid(igk2)
                    erefle(igk1) = erefle(igk1) + qiii(igk1, igk2) * eincid(igk2)
                end do
                !
                down = down + eincid(igk1) * conjg(eincid(igk1)) * gkzin
                trans = trans + etrans(igk1) * conjg(etrans(igk1)) * gkzout
                refle = refle + erefle(igk1) * conjg(erefle(igk1)) * gkzin
            end do
        end do
        trans = trans / down
        refle = refle / down
        absor = 1.d0 - trans - refle
        write(8, 101)  zval, trans, refle, absor
        write(6, 101)  zval, trans, refle, absor
        return
        !
        101 format(5e14.6)
    end subroutine
    !=======================================================================
    complex(dp) function codd(l,m,l1,m1,xodd)
        integer, intent(in) ::l,m,l1,m1
        complex(dp), intent(in) :: xodd(:,:)
        integer i,j
        !     ------------------------------------------------------------------
        if(abs(m)<=l.and.abs(m1)<=l1) then
            i=(l*l+m+1)/2
            j=(l1*l1+m1+1)/2
            codd=xodd(i,j)
        else
            codd=czero
        end if
        return
    end function
    !=======================================================================
    complex(dp) function ceven(l,m,l1,m1,xeven)
        integer l,m,l1,m1
        complex(dp) xeven(:,:)
        integer i,j
        !     ------------------------------------------------------------------
        if(abs(m)<=l.and.abs(m1)<=l1) then
            i=(l*l+2*l+m+2)/2
            j=(l1*l1+2*l1+m1+2)/2
            ceven=xeven(i,j)
        else
            ceven=czero
        end if
        return
    end function
    !=======================================================================
    subroutine sphrm4(ylm,ct,st,cf,lmax)
        !     -----------------------------------------------------------------
        !     given  ct=cos(theta),  st=sin(theta),  and cf=exp(i*fi), this
        !     subroutine  calculates  all the  ylm(theta,fi) up to  l=lmax.
        !     subscripts are ordered thus:(l,m)=(0,0),(1,-1),(1,0),(1,1)...
        !     -----------------------------------------------------------------
        integer    lmax
        complex(dp) ct,st,cf
        complex(dp) ylm(:)
        ! ..  local
        integer    l,ll,lm,lm2,lm3,ln,lo,lp,lq,m
        real(dp)   a,asg,b,cl,cm
        complex(dp) sf,sa
        real(dp)   fac1(lmax+1),fac3(lmax+1),fac2((lmax+1)**2)
        !-----------------------------------------------------------------------
        lm=0
        cl=0.0_dp
        a=1.0_dp
        b=1.0_dp
        asg=1.0_dp
        ll=lmax+1
        !****** multiplicative factors required ******
        do l=1,ll
            fac1(l)=asg*sqrt((2.0_dp*cl+1.0_dp)*a/(4.0_dp*pi*b*b))
            fac3(l)=sqrt(2.0_dp*cl)
            cm=-cl
            ln=l+l-1
            do m=1,ln
                lo=lm+m
                fac2(lo)=sqrt((cl+1.0_dp+cm)*(cl+1.0_dp-cm) &
                        & /((2.0_dp*cl+3.0_dp)*(2.0_dp*cl+1.0_dp)))
                cm=cm+1.0_dp
            end do
            cl=cl+1.0_dp
            a=a*2.0_dp*cl*(2.0_dp*cl-1.0_dp)/4.0_dp
            b=b*cl
            asg=-asg
            lm=lm+ln
        end do
        !****** first all the ylm for m=+-l and m=+-(l-1) are ******
        !****** calculated by explicit formulae               ******
        lm=1
        cl=1.0_dp
        asg=-1.0_dp
        sf=cf
        sa= cmplx(1.0_dp,0.0_dp, kind=dp)
        ylm(1)=cmplx(fac1(1),0.0_dp, kind=dp)
        do l=1,lmax
            ln=lm+l+l+1
            ylm(ln)=fac1(l+1)*sa*sf*st
            ylm(lm+1)=asg*fac1(l+1)*sa*st/sf
            ylm(ln-1)=-fac3(l+1)*fac1(l+1)*sa*sf*ct/cf
            ylm(lm+2)=asg*fac3(l+1)*fac1(l+1)*sa*ct*cf/sf
            sa=st*sa
            sf=sf*cf
            cl=cl+1.0_dp
            asg=-asg
            lm=ln
        end do
        !****** using ylm and yl(m-1) in a recurence relation ******
        !****** yl(m+1) is calculated                         ******
        lm=1
        ll=lmax-1
        do l=1,ll
            ln=l+l-1
            lm2=lm+ln+4
            lm3=lm-ln
            do m=1,ln
                lo=lm2+m
                lp=lm3+m
                lq=lm+m+1
                ylm(lo)=-(fac2(lp)*ylm(lp)-ct*ylm(lq))/fac2(lq)
            end do
            lm=lm+l+l+1
        end do
        return
    end subroutine
    !=======================================================================
    subroutine bessel(BJ,Y,H,arg)
        !     ------------------------------------------------------------------
        !     THIS  SUBROUTINE COMPUTES THE  SPHERICAL BESSEL FUNCTIONS OF
        !     FIRST, SECOND  AND  THIRD  KIND  using Amos lib
        !
        !     2019.07.17 Change to use Amos lib
        !
        !     ON INPUT--->
        !     ARG    ARGUMENT OF THE BESSEL FUNCTIONS
        !     ON OUTPUT--->
        !     BJ     AN ARRAY CONTAINING THE BESSEL FUNCTIONS OF
        !            THE FIRST KIND UP TO LMAX1 IF LJ IS TRUE.
        !            REMEMBER, THAT BJ(1) CONTAINS THE FUNCTION OF
        !            L=0 AND SO ON.
        !     Y      AN ARRAY CONTAINING THE BESSEL FUNCTIONS OF
        !            THE SECOND KIND UP TO LMAX1 IF LY IS TRUE.
        !            REMEMBER,THAT  Y(1) CONTAINS THE FUNCTION OF L=0 AND SO ON.
        !     H      AN ARRAY CONTAINING THE BESSEL FUNCTIONS OF
        !            THE THIRD KIND UP TO LMAX1 IF LH IS TRUE.
        !            REMEMBER,THAT H (1) CONTAINS THE FUNCTION OF L=0 AND SO ON.
        !
        !     THE BESSEL FUNCTIONS OF 3RD KIND ARE DEFINED AS: H(L)=BJ(L)+I*Y(L)
        !     ------------------------------------------------------------------
        complex(dp), intent(in) :: arg
        complex(dp), intent(out) :: BJ(:),H(:),Y(:)
        ! local
        integer      :: lmax1
        INTEGER KODE, N, NZ, IERR
        real(dp)     :: zr, zi, FNU
        real(dp), allocatable :: cyr(:), cyi(:), cwrkr(:), cwrki(:)
        complex(dp),allocatable :: cy(:)
        !-----------------------------------------------------------------------
        lmax1 = size(BJ) ! to store from l=0 to l=lmax
        if (size(BJ)/=size(Y) .or. size(BJ)/=size(H)) stop 1
        allocate(cy(1:lmax1)); allocate(cyr(1:lmax1)); allocate(cyi(1:lmax1))
        allocate(cwrki(1:lmax1));  allocate(cwrkr(1:lmax1))
        zr = real(arg); zi = aimag(arg)
        FNU = 0.5_dp;   KODE=1;   N=lmax1
        call ZBESJ(zr, zi, FNU, KODE, N, CYR, CYI, NZ, IERR)
        if (IERR /= 0) stop 1
        ! Convert to spherical function
        cy = (cyr+ ci*cyi)*sqrt(pi/2.0_dp/arg)
        BJ = cy
        cwrkr=0.0_dp; cwrki=0.0_dp
        call ZBESY(zr, zi, FNU, KODE, N, CYR, CYI, NZ, cwrkr, cwrki, IERR)
        if (IERR /= 0) stop 1
        cy = (cyr+ ci*cyi)*sqrt(pi/2.0_dp/arg)
        Y = cy
        H=BJ+ci*Y
    end subroutine
    !=======================================================================
    subroutine tmtrx(rap,epssph,epsmed,mumed,musph,TE,TH)
        !     ------------------------------------------------------------------
        !     THIS SUBROUTINE  CALCULATES  THE  T-MATRIX FOR THE SCATTERING
        !     OF ELECTROMAGNETIC  FIELD  OF  WAVE-LENGHT LAMDA  BY A SINGLE
        !     SPHERE OF RADIUS S.  (RAP=S/LAMDA).
        !     EPSSPH : COMPLEX RELATIVE DIELECTRIC CONSTANT OF THE SPHERE.
        !     EPSMED : COMPLEX RELATIVE DIELECTRIC CONSTANT OF THE MEDIUM.
        !     LMAX   : MAXIMUM ANGULAR MOMENTUM from TE(0..LMAX) and TH
        !     ------------------------------------------------------------------
        complex(dp), intent(in) :: EPSSPH,EPSMED,MUSPH,MUMED,RAP
        complex(dp), intent(out) :: TE(:),TH(:)
        ! local
        INTEGER  ::  l1, lmax, lmax1, b_size
        complex(dp) :: C1,C2,C3,C4,C5,C6,AN,AJ,BN,BJ,ARG,ARGM,XISQ,XISQM,AR
        complex(dp), allocatable:: J(:),Y(:),H(:),JM(:),YM(:),HM(:)
        !-----------------------------------------------------------------------
        lmax1 = size(TE)
        ! to evaluate TE(0..lmax) we need one more oder in Bessel functions
        b_size = lmax1+1
        allocate(J (1:b_size)); allocate(Y (1:b_size)); allocate(H (1:b_size))
        allocate(JM(1:b_size)); allocate(YM(1:b_size)); allocate(HM(1:b_size))
        lmax = lmax1-1
        ar=2.0_dp*pi*rap
        xisq =sqrt(epsmed*mumed);  arg =xisq *ar
        xisqm=sqrt(epssph*musph);  argm=xisqm*ar
        call bessel(J,Y,H,arg);  call bessel(JM,YM,HM,argm)
        c1=epssph-epsmed;   c2=epsmed*argm;   c3=-epssph*arg
        c4= musph -mumed;   c5= mumed*argm;   c6= -musph*arg
        do  L1=1,LMAX1
            an=C1*L1*JM(L1)*Y(L1)+C2*JM(L1+1)*Y(L1)+C3*JM(L1)*Y(L1+1)
            aj=C1*L1*JM(L1)*J(L1)+C2*JM(L1+1)*J(L1)+C3*JM(L1)*J(L1+1)
            bn=C4*L1*JM(L1)*Y(L1)+C5*JM(L1+1)*Y(L1)+C6*JM(L1)*Y(L1+1)
            bj=C4*L1*JM(L1)*J(L1)+C5*JM(L1+1)*J(L1)+C6*JM(L1)*J(L1+1)
            TE(L1)=-aj/(aj+ci*an)
            TH(L1)=-bj/(bj+ci*bn)
        end do
        return
    end subroutine

end module