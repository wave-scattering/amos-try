!********************************************************************

!   calculation of the amplitude matrix

subroutine ampl (nmax, dlam, tl, tl1, pl, pl1, alpha, beta, &
        vv, vh, hv, hh, errcode)
    include 'amplq.par.f'
    implicit real(kind = 8) (a-b, d-h, o-z), complex(kind = 8) (c)
    real(kind = 8) al(3, 2), al1(3, 2), ap(2, 3), ap1(2, 3), b(3, 3), &
            r(2, 2), r1(2, 2), c(3, 2), ca, cb, ct, cp, ctp, cpp, ct1, cp1, &
            ctp1, cpp1
    real(kind = 8) dv1(npn6), dv2(npn6), dv01(npn6), dv02(npn6)
    real(kind = 4)&
            tr11(npn6, npn4, npn4), tr12(npn6, npn4, npn4), &
            tr21(npn6, npn4, npn4), tr22(npn6, npn4, npn4), &
            ti11(npn6, npn4, npn4), ti12(npn6, npn4, npn4), &
            ti21(npn6, npn4, npn4), ti22(npn6, npn4, npn4)
    complex(kind = 8) cal(npn4, npn4), vv, vh, hv, hh
    integer errcode

    common /tmat/ tr11, tr12, tr21, tr22, ti11, ti12, ti21, ti22

    if (alpha<0d0.or.alpha>360d0.or.&
            beta<0d0.or.beta>180d0.or.&
            tl<0d0.or.tl>180d0.or.&
            tl1<0d0.or.tl1>180d0.or.&
            pl<0d0.or.pl>360d0.or.&
            pl1<0d0.or.pl1>360d0) then
        !          write (6,2000)
        !          stop
        errcode = 3
        return
    else

    endif
    ! 2000 format ('an angular parameter is outside its',
    !     &        ' allowable range')
    pin = dacos(-1d0)
    pin2 = pin * 0.5d0
    pi = pin / 180d0
    alph = alpha * pi
    bet = beta * pi
    thetl = tl * pi
    phil = pl * pi
    thetl1 = tl1 * pi
    phil1 = pl1 * pi

    eps = 1d-8
    if (thetl<pin2) thetl = thetl + eps
    if (thetl>pin2) thetl = thetl - eps
    if (thetl1<pin2) thetl1 = thetl1 + eps
    if (thetl1>pin2) thetl1 = thetl1 - eps
    if (phil<pin) phil = phil + eps
    if (phil>pin) phil = phil - eps
    if (phil1<pin) phil1 = phil1 + eps
    if (phil1>pin) phil1 = phil1 - eps
    if (bet<=pin2.and.pin2 - bet<=eps) bet = bet - eps
    if (bet>pin2.and.bet - pin2<=eps) bet = bet + eps

    !_____________compute thetp, phip, thetp1, and phip1, eqs. (8), (19), and (20)

    cb = dcos(bet)
    sb = dsin(bet)
    ct = dcos(thetl)
    st = dsin(thetl)
    cp = dcos(phil - alph)
    sp = dsin(phil - alph)
    ctp = ct * cb + st * sb * cp
    thetp = dacos(ctp)
    cpp = cb * st * cp - sb * ct
    spp = st * sp
    phip = datan(spp / cpp)
    if (phip>0d0.and.sp<0d0) phip = phip + pin
    if (phip<0d0.and.sp>0d0) phip = phip + pin
    if (phip<0d0) phip = phip + 2d0 * pin

    ct1 = dcos(thetl1)
    st1 = dsin(thetl1)
    cp1 = dcos(phil1 - alph)
    sp1 = dsin(phil1 - alph)
    ctp1 = ct1 * cb + st1 * sb * cp1
    thetp1 = dacos(ctp1)
    cpp1 = cb * st1 * cp1 - sb * ct1
    spp1 = st1 * sp1
    phip1 = datan(spp1 / cpp1)
    if (phip1>0d0.and.sp1<0d0) phip1 = phip1 + pin
    if (phip1<0d0.and.sp1>0d0) phip1 = phip1 + pin
    if (phip1<0d0) phip1 = phip1 + 2d0 * pin

    !____________compute matrix beta, eq. (21)

    ca = dcos(alph)
    sa = dsin(alph)
    b(1, 1) = ca * cb
    b(1, 2) = sa * cb
    b(1, 3) = -sb
    b(2, 1) = -sa
    b(2, 2) = ca
    b(2, 3) = 0d0
    b(3, 1) = ca * sb
    b(3, 2) = sa * sb
    b(3, 3) = cb

    !____________compute matrices al and al1, eq. (14)

    cp = dcos(phil)
    sp = dsin(phil)
    cp1 = dcos(phil1)
    sp1 = dsin(phil1)
    al(1, 1) = ct * cp
    al(1, 2) = -sp
    al(2, 1) = ct * sp
    al(2, 2) = cp
    al(3, 1) = -st
    al(3, 2) = 0d0
    al1(1, 1) = ct1 * cp1
    al1(1, 2) = -sp1
    al1(2, 1) = ct1 * sp1
    al1(2, 2) = cp1
    al1(3, 1) = -st1
    al1(3, 2) = 0d0

    !____________compute matrices ap^(-1) and ap1^(-1), eq. (15)

    ct = ctp
    st = dsin(thetp)
    cp = dcos(phip)
    sp = dsin(phip)
    ct1 = ctp1
    st1 = dsin(thetp1)
    cp1 = dcos(phip1)
    sp1 = dsin(phip1)
    ap(1, 1) = ct * cp
    ap(1, 2) = ct * sp
    ap(1, 3) = -st
    ap(2, 1) = -sp
    ap(2, 2) = cp
    ap(2, 3) = 0d0
    ap1(1, 1) = ct1 * cp1
    ap1(1, 2) = ct1 * sp1
    ap1(1, 3) = -st1
    ap1(2, 1) = -sp1
    ap1(2, 2) = cp1
    ap1(2, 3) = 0d0

    !____________compute matrices r and r^(-1), eq. (13)
    do i = 1, 3
        do j = 1, 2
            x = 0d0
            do k = 1, 3
                x = x + b(i, k) * al(k, j)
            enddo
            c(i, j) = x
        enddo
    enddo
    do i = 1, 2
        do j = 1, 2
            x = 0d0
            do k = 1, 3
                x = x + ap(i, k) * c(k, j)
            enddo
            r(i, j) = x
        enddo
    enddo
    do i = 1, 3
        do j = 1, 2
            x = 0d0
            do k = 1, 3
                x = x + b(i, k) * al1(k, j)
            enddo
            c(i, j) = x
        enddo
    enddo
    do i = 1, 2
        do j = 1, 2
            x = 0d0
            do k = 1, 3
                x = x + ap1(i, k) * c(k, j)
            enddo
            r1(i, j) = x
        enddo
    enddo
    d = 1d0 / (r1(1, 1) * r1(2, 2) - r1(1, 2) * r1(2, 1))
    x = r1(1, 1)
    r1(1, 1) = r1(2, 2) * d
    r1(1, 2) = -r1(1, 2) * d
    r1(2, 1) = -r1(2, 1) * d
    r1(2, 2) = x * d

    ci = (0d0, 1d0)
    do nn = 1, nmax
        do n = 1, nmax
            cn = ci**(nn - n - 1)
            dnn = dfloat((2 * n + 1) * (2 * nn + 1))
            dnn = dnn / dfloat(n * nn * (n + 1) * (nn + 1))
            rn = dsqrt(dnn)
            cal(n, nn) = cn * rn
        end do
    end do
    dcth0 = ctp
    dcth = ctp1
    ph = phip1 - phip
    vv = (0d0, 0d0)
    vh = (0d0, 0d0)
    hv = (0d0, 0d0)
    hh = (0d0, 0d0)
    do m = 0, nmax
        m1 = m + 1
        nmin = max(m, 1)
        call vigampl (dcth, nmax, m, dv1, dv2)
        call vigampl (dcth0, nmax, m, dv01, dv02)
        fc = 2d0 * dcos(m * ph)
        fs = 2d0 * dsin(m * ph)
        do nn = nmin, nmax
            dv1nn = m * dv01(nn)
            dv2nn = dv02(nn)
            do n = nmin, nmax
                dv1n = m * dv1(n)
                dv2n = dv2(n)

                ct11 = dcmplx(tr11(m1, n, nn), ti11(m1, n, nn))
                ct22 = dcmplx(tr22(m1, n, nn), ti22(m1, n, nn))

                if (m==0) then

                    cn = cal(n, nn) * dv2n * dv2nn

                    vv = vv + cn * ct22
                    hh = hh + cn * ct11

                else

                    ct12 = dcmplx(tr12(m1, n, nn), ti12(m1, n, nn))
                    ct21 = dcmplx(tr21(m1, n, nn), ti21(m1, n, nn))

                    cn1 = cal(n, nn) * fc
                    cn2 = cal(n, nn) * fs

                    d11 = dv1n * dv1nn
                    d12 = dv1n * dv2nn
                    d21 = dv2n * dv1nn
                    d22 = dv2n * dv2nn

                    vv = vv + (ct11 * d11 + ct21 * d21&
                            + ct12 * d12 + ct22 * d22) * cn1

                    vh = vh + (ct11 * d12 + ct21 * d22&
                            + ct12 * d11 + ct22 * d21) * cn2

                    hv = hv - (ct11 * d21 + ct21 * d11&
                            + ct12 * d22 + ct22 * d12) * cn2

                    hh = hh + (ct11 * d22 + ct21 * d12&
                            + ct12 * d21 + ct22 * d11) * cn1
                endif
            end do
        end do
    end do
    dk = 2d0 * pin / dlam
    vv = vv / dk
    vh = vh / dk
    hv = hv / dk
    hh = hh / dk
    cvv = vv * r(1, 1) + vh * r(2, 1)
    cvh = vv * r(1, 2) + vh * r(2, 2)
    chv = hv * r(1, 1) + hh * r(2, 1)
    chh = hv * r(1, 2) + hh * r(2, 2)
    vv = r1(1, 1) * cvv + r1(1, 2) * chv
    vh = r1(1, 1) * cvh + r1(1, 2) * chh
    hv = r1(2, 1) * cvv + r1(2, 2) * chv
    hh = r1(2, 1) * cvh + r1(2, 2) * chh

    !      write (6,1005) tl,tl1,pl,pl1,alpha,beta
    !      write (6,1006)
    !      print 1101, vv
    !      print 1102, vh
    !      print 1103, hv
    !      print 1104, hh
    ! 1101 format ('s11=',d11.5,' + i*',d11.5)
    ! 1102 format ('s12=',d11.5,' + i*',d11.5)
    ! 1103 format ('s21=',d11.5,' + i*',d11.5)
    ! 1104 format ('s22=',d11.5,' + i*',d11.5)
    ! 1005 format ('thet0=',f6.2,'  thet=',f6.2,'  phi0=',f6.2,
    !     &        '  phi=',f6.2,'  alpha=',f6.2,'  beta=',f6.2)
    ! 1006 format ('amplitude matrix')
    return
end


!*****************************************************************
!
!     calculation of the functions
!     dv1(n)=dvig(0,m,n,arccos x)/sin(arccos x)
!     and
!     dv2(n)=[d/d(arccos x)] dvig(0,m,n,arccos x)
!     1<=n<=nmax
!     0<=x<=1

subroutine vigampl (x, nmax, m, dv1, dv2)
    include 'amplq.par.f'
    implicit real(kind = 8) (a-h, o-z)
    real(kind = 8) dv1(npn6), dv2(npn6)
    do n = 1, nmax
        dv1(n) = 0d0
        dv2(n) = 0d0
    end do
    dx = dabs(x)
    if (dabs(1d0 - dx)>1d-10) then
        a = 1d0
        qs = dsqrt(1d0 - x * x)
        qs1 = 1d0 / qs
        dsi = qs1
        if (m==0) then
            d1 = 1d0
            d2 = x
            do n = 1, nmax
                qn = dfloat(n)
                qn1 = dfloat(n + 1)
                qn2 = dfloat(2 * n + 1)
                d3 = (qn2 * x * d2 - qn * d1) / qn1
                der = qs1 * (qn1 * qn / qn2) * (-d1 + d3)
                dv1(n) = d2 * dsi
                dv2(n) = der
                d1 = d2
                d2 = d3
            end do
            return
        end if
        qmm = dfloat(m * m)
        do i = 1, m
            i2 = i * 2
            a = a * dsqrt(dfloat(i2 - 1) / dfloat(i2)) * qs
        end do
        d1 = 0d0
        d2 = a
        do n = m, nmax
            qn = dfloat(n)
            qn2 = dfloat(2 * n + 1)
            qn1 = dfloat(n + 1)
            qnm = dsqrt(qn * qn - qmm)
            qnm1 = dsqrt(qn1 * qn1 - qmm)
            d3 = (qn2 * x * d2 - qnm * d1) / qnm1
            der = qs1 * (-qn1 * qnm * d1 + qn * qnm1 * d3) / qn2
            dv1(n) = d2 * dsi
            dv2(n) = der
            d1 = d2
            d2 = d3
        end do
        return
    end if
    if (m/=1) return
    do n = 1, nmax
        dn = dfloat(n * (n + 1))
        dn = 0.5d0 * dsqrt(dn)
        if (x<0d0) dn = dn * (-1)**(n + 1)
        dv1(n) = dn
        if (x<0d0) dn = -dn
        dv2(n) = dn
    end do
    return
end


!**********************************************************************
!                                                                     *
!   input parameters:                                                 *
!                                                                     *
!   ng = 2*ngauss - number of quadrature points on the                *
!                   interval  (-1,1). ngauss.le.npng1                 *
!   nmax,mmax - maximum dimensions of the arrays.  nmax.le.npn1       *
!               mmax.le.npn1                                          *
!   p - pi                                                            *
!                                                                     *
!   output parameters:                                                *
!                                                                     *
!   x,w - points and weights of the quadrature formula                *
!   an(n) = n*(n+1)                                                   *
!   ann(n1,n2) = (1/2)*sqrt((2*n1+1)*(2*n2+1)/(n1*(n1+1)*n2*(n2+1)))  *
!   s(i)=1/sin(arccos(x(i)))                                          *
!   ss(i)=s(i)**2                                                     *
!                                                                     *
!**********************************************************************

subroutine constant (ngauss, nmax, mmax, p, x, w, an, ann, s, ss, np, eps)
    implicit real(kind = 16) (a-h, o-z)
    include 'amplq.par.f'
    real(kind = 16) x(npng2), w(npng2), x1(npng1), w1(npng1), &
            x2(npng1), w2(npng1), &
            s(npng2), ss(npng2), &
            an(npn1), ann(npn1, npn1), dd(npn1)

    do n = 1, nmax
        nn = n * (n + 1)
        an(n) = qfloat(nn)
        d = qsqrt(qfloat(2 * n + 1) / qfloat(nn))
        dd(n) = d
        do n1 = 1, n
            ddd = d * dd(n1) * 0.5q0
            ann(n, n1) = ddd
            ann(n1, n) = ddd
        end do
    end do
    ng = 2 * ngauss
    if (np==-2) then
        ng1 = dfloat(ngauss) / 2d0
        ng2 = ngauss - ng1
        xx = -qcos(qatan(eps))
        call qgauss(ng1, 0, 0, x1, w1)
        call qgauss(ng2, 0, 0, x2, w2)
        do i = 1, ng1
            w(i) = 0.5q0 * (xx + 1q0) * w1(i)
            x(i) = 0.5q0 * (xx + 1q0) * x1(i) + 0.5q0 * (xx - 1q0)
        end do
        do i = 1, ng2
            w(i + ng1) = -0.5q0 * xx * w2(i)
            x(i + ng1) = -0.5q0 * xx * x2(i) + 0.5q0 * xx
        end do
        do i = 1, ngauss
            w(ng - i + 1) = w(i)
            x(ng - i + 1) = -x(i)
        end do
    else
        call qgauss(ng, 0, 0, x, w)
        do i = 1, ngauss
            y = x(i)
            y = 1q0 / (1q0 - y * y)
            ss(i) = y
            ss(ng - i + 1) = y
            y = qsqrt(y)
            s(i) = y
            s(ng - i + 1) = y
        end do
    end if
    return
end

!***************************************************************

subroutine qgauss (n, ind1, ind2, z, w)
    implicit real(kind = 16) (a-h, p-z)
    real(kind = 16) z(n), w(n)
    a = 1q0
    b = 2q0
    c = 3q0
    ind = mod(n, 2)
    k = n / 2 + ind
    f = qfloat(n)
    do i = 1, k
        m = n + 1 - i
        select case (i)
        case(1)
            x = a - b / ((f + a) * f)
        case(2)
            x = (z(n) - a) * 4q0 + z(n)
        case(3)
            x = (z(n - 1) - z(n)) * 1.6q0 + z(n - 1)
        case default
            x = (z(m + 1) - z(m + 2)) * c + z(m + 3)
        end select

        if(i==k.and.ind==1) x = 0q0
        niter = 0
        check = 1q-32
        do
            pb = 1q0
            niter = niter + 1
            if (niter>100) then
                !print 5000, check
                check = check * 10q0
            end if
            pc = x
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
            if (qabs(pb)<=check * qabs(x)) exit
        end do
        z(m) = x
        w(m) = pa * pa * (a - x * x)
        if(ind1==0) w(m) = b * w(m)
        if(i/=k.or.ind/=1) then
            z(i) = -z(m)
            w(i) = w(m)
        end if
    end do
    ! 5000 format ('qgauss does not converge, check=',q10.3)
    if (ind2==1) then
        !      print 1100,n
        ! 1100 format(' ***  points and weights of gaussian quadrature formula',
        !     * ' of ',i4,'-th order')
        do i = 1, k
            zz = -z(i)
            !  105     print 1200,i,zz,i,w(i)
        end do
        ! 1200 format(' ',4x,'x(',i4,') = ',f17.14,5x,'w(',i4,') = ',f17.14)
    else
        !     print 1300,n
        ! 1300 format(' gaussian quadrature formula of ',i4,'-th order is used')
        if (ind1/=0) then
            do i = 1, n
                z(i) = (a + z(i)) / b
            end do
        end if
    end if
    return
end