module special_functions
    use constants
    !use dense_solve
    !use multem_blas
    use amos
    !use errfun, only : wpop
    implicit none

    type, private :: abess_values
        integer :: nnmax1
        real(dp), allocatable :: AY(:), ADY(:), AJ(:), ADJ(:), &
                AJR(:), ADJR(:), AJI(:), ADJI(:)
    end type abess_values
    type(abess_values), private :: abess

    type, private :: cbess_values
        real(dp) :: wv ! wave vector = 2*pi/lambda
        real(dp) :: mrr, mri ! real and imag parts of refractive index
        real(dp), allocatable :: J(:, :), Y(:, :), JR(:, :), JI(:, :), &
                DJ(:, :), DY(:, :), DJR(:, :), DJI(:, :)
    end type cbess_values
    type(cbess_values), public :: cbess

    public vig, vigf !,zartan


contains
    !=======================================================================
    !=======================================================================
    subroutine cbesscjcdj(x, n, nmax, jr, ji, djr, dji)
        real(dp) x, xx, xr, xi, jr, ji, djr, dji
        integer nmax, n
        !
        call reallocate_abess(nmax)
        xx = dsqrt(x)*cbess%wv
        xr = xx*cbess%mrr
        xi = xx*cbess%mri
        call cjb(xr, xi, abess%ajr, abess%aji, &
                abess%adjr, abess%adji, nmax, 2)
        jr = abess%ajr(n)
        ji = abess%aji(n)
        djr = abess%adjr(n)
        dji = abess%adji(n)
        return
    end subroutine cbesscjcdj
    !=======================================================================
    subroutine cbessydy(x, nmax, y, dy)
        real(dp) x, xx, y, dy
        integer nmax
        !
        call reallocate_abess(nmax)
        xx = dsqrt(x)*cbess%wv
        call ryb(xx, abess%ay, abess%ady, nmax)
        y = abess%ay(nmax)
        dy = abess%ady(nmax)
        return
    end subroutine cbessydy
    !=======================================================================
    subroutine cbessjdj(x, nmax, j, dj)
        real(dp) x, xx, j, dj
        integer nmax
        !
        call reallocate_abess(nmax)
        xx = dsqrt(x)*cbess%wv
        call rjb(xx, abess%aj, abess%adj, nmax, abess%nnmax1)
        j = abess%aj(nmax)
        dj = abess%adj(nmax)
        return
    end subroutine cbessjdj
    !=======================================================================
    subroutine bess_amos(bj, y, h, arg)
        !     ------------------------------------------------------------------
        !     This  subroutine computes the  spherical Bessel functions of
        !     first, second  and  third  kind  using Amos lib
        !
        !     2019.07.17 change to use Amos lib
        !
        !     on input--->
        !     arg    argument of the Bessel functions
        !     on output--->
        !     bj     an array containing the Bessel functions of
        !            the first kind up to lmax1 if lj is true.
        !            remember, that bj(1) contains the function of
        !            l=0 and so on.
        !     y      an array containing the Bessel functions of
        !            the second kind up to lmax1 if ly is true.
        !            remember,that  y(1) contains the function of l=0 and so on.
        !     h      an array containing the Bessel functions of
        !            the third kind up to lmax1 if lh is true.
        !            remember,that h (1) contains the function of l=0 and so on.
        !
        !     the Bessel functions of 3rd kind are defined as: h(l)=bj(l)+i*y(l)
        !     ------------------------------------------------------------------
        complex(dp), intent(in) :: arg
        complex(dp), intent(out) :: bj(:), h(:), y(:)
        ! local
        integer :: lmax1
        integer kode, n, nz, ierr
        real(dp) :: zr, zi, fnu
        real(dp), allocatable :: cyr(:), cyi(:), cwrkr(:), cwrki(:)
        complex(dp), allocatable :: cy(:)
        !-----------------------------------------------------------------------
        lmax1 = size(bj) ! to store from l=0 to l=lmax
        if (size(bj)/=size(y) .or. size(bj)/=size(h)) stop 1
        allocate(cy(1:lmax1)); allocate(cyr(1:lmax1)); allocate(cyi(1:lmax1))
        allocate(cwrki(1:lmax1));  allocate(cwrkr(1:lmax1))
        zr = real(arg); zi = aimag(arg)
        fnu = 0.5_dp;   kode = 1;   n = lmax1
        call zbesj(zr, zi, fnu, kode, n, cyr, cyi, nz, ierr)
        if (ierr /= 0) stop 1
        ! convert to spherical function
        cy = (cyr + ci * cyi) * sqrt(pi / 2.0_dp / arg)
        bj = cy
        cwrkr = 0.0_dp; cwrki = 0.0_dp
        call zbesy(zr, zi, fnu, kode, n, cyr, cyi, nz, cwrkr, cwrki, ierr)
        if (ierr /= 0) stop 1
        cy = (cyr + ci * cyi) * sqrt(pi / 2.0_dp / arg)
        y = cy
        h = bj + ci * y
    end subroutine
    !=======================================================================
    subroutine bess (x, xr, xi, ng, nmax, nnmax1)
        !--------/---------/---------/---------/---------/---------/---------/--
        ! >>> x,xr,xi,ng,nmax,nnmax1,nnmax2
        ! <<< output j,y,jr,ji,dj,dy,djr,dji  to common block cbess
        !==========================================================
        !  generates bessel functions for each gauss integration point
        !
        !  x =(2\pi/\lambda)*r
        !  xr=(2\pi/\lambda)*r*mrr, mrr ... real part of the rel. refractive index
        !  xi=(2\pi/\lambda)*r*mri, mri ... imag. part of the rel. refractive index
        !  ng=2*ngauss or 60 ... the number of gauss integration points
        !  j,y,jr,ji ... arrays of bessel functions
        !  dj,dy,djr,dji  ... arrays of bessel functions derivatives of the form
        !                           [xf(x)]'/x                   (a)
        !                 where prime denotes derivative with respect to x.
        !                 (note that bessel function derivatives enter eqs. (39)
        !                  \cite{tks} only in the (a) combination!!!!)
        !  nmax   ... angular momentum cutoff
        !  nnmax1 ... angular momentum cutoff - determines numerical accuracy
        !  nnmax2 ... angular momentum cutoff - determines numerical accuracy
        !--------/---------/---------/---------/---------/---------/---------/--
        integer ng, nnmax1, i, n
        real(dp) xx, yi, yr
        integer, intent(in) :: nmax
        real(dp) x(ng), xr(ng), xi(ng) ! TODO make them as x(:), etc.
        !
        call reallocate_abess(nmax)
        call reallocate_cbess(ng, nmax)
        abess%nnmax1 = nnmax1
        ng = size(x)
        do i = 1, ng
            xx = x(i)
            yr = xr(i)
            yi = xi(i)
            !
            call rjb(xx, abess%aj, abess%adj, nmax, nnmax1)
            call ryb(xx, abess%ay, abess%ady, nmax)
            call cjb(yr, yi, abess%ajr, abess%aji, abess%adjr, abess%adji, nmax, 2)
!            call cjb_amos(yr, yi, abess%ajr, abess%aji, abess%adjr, abess%adji, nmax)
            !
            do n = 1, nmax
                cbess%j(i, n) = abess%aj(n)
                cbess%y(i, n) = abess%ay(n)
                cbess%jr(i, n) = abess%ajr(n)
                cbess%ji(i, n) = abess%aji(n)
                cbess%dj(i, n) = abess%adj(n)
                cbess%dy(i, n) = abess%ady(n)
                cbess%djr(i, n) = abess%adjr(n)
                cbess%dji(i, n) = abess%adji(n)
            end do
        end do
        return
    end
    !=======================================================================
    subroutine reallocate_cbess(ng, nmax)
        integer nmax, nmax_old, ng, ng_old
        integer :: array_shape(2)
        integer :: i, k
        !     ------------------------------------------------------------------
        ng_old = ng;             nmax_old = nmax
        i = 1;                   k = 1
        if (allocated(cbess%J)) then
            ! reallocate only if needed, so to increase the computational speed
            array_shape = shape(cbess%J)
            ng_old = array_shape(1)
            nmax_old = array_shape(2)
            if (ng_old < ng) i = 2
            if (nmax_old < nmax) k = 2
            if (nmax_old < nmax .or. ng_old < ng) then
                deallocate(cbess%J, cbess%Y, cbess%JR, cbess%JI, &
                        cbess%DJ, cbess%DY, cbess%DJR, cbess%DJI)
            endif
        endif
        if (.not.allocated(cbess%J)) then
            ! Overallocate to reduce number of allocations
            allocate(cbess%J (ng*i, nmax*k), cbess%Y  (ng*i, nmax*k), &
                    cbess%JR (ng*i, nmax*k), cbess%JI (ng*i, nmax*k), &
                    cbess%DJ (ng*i, nmax*k), cbess%DY (ng*i, nmax*k), &
                    cbess%DJR(ng*i, nmax*k), cbess%DJI(ng*i, nmax*k))
        endif
    end

    !=======================================================================
    subroutine reallocate_abess(nmax)
        integer nmax, nmax_old
        integer, parameter :: g = 2
        if (allocated(abess%AY)) then
            ! reallocate only if needed, so to increase the computational speed
            nmax_old = size(abess%AY)
            if (nmax_old < nmax) then
                deallocate(abess%AY, abess%ADY, abess%AJ, abess%ADJ)
                deallocate(abess%AJR, abess%ADJR, abess%AJI, abess%ADJI)
            endif
        endif
        if (.not.allocated(abess%AY)) then
            ! Overallocate to reduce number of allocations
            allocate(abess%AY (nmax*g), abess%ADY (nmax*g), abess%AJ (nmax*g), abess%ADJ (nmax*g))
            allocate(abess%AJR(nmax*g), abess%ADJR(nmax*g), abess%AJI(nmax*g), abess%ADJI(nmax*g))
        endif
    end
    !=======================================================================
    subroutine cjb_amos(zr, zi, yr, yi, ur, ui, nmax)
        !     ------------------------------------------------------------------
        !     This subroutine computes the spherical Bessel functions of
        !     first kind and its derivative using Amos lib
        !     from 1 to nmax (without zero order)
        !     ------------------------------------------------------------------
        real(dp), intent(in) :: zr, zi
        real(dp), intent(out) :: yr(:), yi(:), ur(:), ui(:)
        integer, intent(in) :: nmax
        ! local
        integer kode, n, nz, ierr, nmax1, k
        real(dp) :: fnu
        real(dp), allocatable :: cyr(:), cyi(:), cwrkr(:), cwrki(:)
        complex(dp), allocatable :: cy(:), cdy(:)
        complex(dp) :: z
        !-----------------------------------------------------------------------
        nmax1 = nmax+1
        allocate(cy(1:nmax1));
        allocate(cdy(1:nmax));
        allocate(cyr(1:nmax1)); allocate(cyi(1:nmax1))
        allocate(cwrki(1:nmax1));  allocate(cwrkr(1:nmax1))
        fnu = 0.5_dp;   kode = 1;   n = nmax1
        call zbesj(zr, zi, fnu, kode, n, cyr, cyi, nz, ierr)
        if (ierr /= 0) stop "Error in zbesj() call"
        ! convert to spherical function
        z = zr + ci*zi
        cy = (cyr + ci*cyi) * sqrt(pi / (2.0_dp*z))
        yr(1:nmax) = real(cy(2:nmax+1))
        yi(1:nmax) = aimag(cy(2:nmax+1))
        do k=2, nmax
            n = k
            ! cdy(n) stores n-th derivative, n=1,2,...
            ! cy(n+1) stores n th value, n=0,1,2,...
            cdy(n-1) = cy(n-1) - (n/z)*cy(n) ! DLFM, eq. 10.51.2
        end do
        ! u = j(x)/x + j'(x)
        ur(1:nmax) = real(cy(2:nmax+1)/z + cdy(1:nmax))
        ui(1:nmax) = aimag(cy(2:nmax+1)/z + cdy(1:nmax))
    end subroutine

    !=======================================================================
    subroutine cjb (xr, xi, yr, yi, ur, ui, nmax, nnmax)
        !--------/---------/---------/---------/---------/---------/---------/--
        !
        !   calculation of spherical bessel functions of the first kind
        !   j=jr+i*ji of complex argument x=xr+i*xi of orders from 1 to nmax
        !   by using backward recursion. parameter nnmax determines numerical
        !   accuracy. u=ur+i*ui - function (1/x)(d/dx)(x*j(x))=j(x)/x + j'(x)
        !
        !  xr=(2\pi/\lambda)*r*mrr, mrr ... real part of the rel. refractive index
        !  xi=(2\pi/\lambda)*r*mri, mri ... imag. part of the rel. refractive index
        !
        !   nmax  - angular momentum cutoff
        !   nnmax - angular momentum cutoff - determines numerical accuracy

        !--------/---------/---------/---------/---------/---------/---------/--
        integer nmax, nnmax, l, l1, i, i1
        real(dp) xr, xi, xrxi, cxxr, qf, ar, ai, ari, cz0r, cz0i, cr, cci, &
                cy0r, cy0i, cy1r, cy1i, cu1r, cu1i, qi, cuii, cuir, cxxi, cyi1i, cyi1r, cyii, &
                cyir
        real(dp) yr(nmax), yi(nmax), ur(nmax), ui(nmax)
        real(dp) cyr(npn1), cyi(npn1), czr(1200), czi(1200)
        !    *      cur(npn1),cui(npn1)
        !
        l = nmax + nnmax
        xrxi = 1d0/(xr*xr + xi*xi)
        cxxr = xr*xrxi             !re [1/(xr+i*xi)]
        cxxi = -xi*xrxi            !im [1/(xr+i*xi)]
        qf = 1d0/dble(2*l + 1)
        czr(l) = xr*qf
        czi(l) = xi*qf
        l1 = l - 1
        do i = 1, l1
            i1 = l - i
            qf = dble(2*i1 + 1)
            ar = qf*cxxr - czr(i1 + 1)
            ai = qf*cxxi - czi(i1 + 1)
            ari = 1d0/(ar*ar + ai*ai)
            czr(i1) = ar*ari
            czi(i1) = -ai*ari
        enddo

        ar = cxxr - czr(1)
        ai = cxxi - czi(1)
        ari = 1d0/(ar*ar + ai*ai)
        cz0r = ar*ari
        cz0i = -ai*ari
        cr = dcos(xr)*dcosh(xi)
        cci = -dsin(xr)*dsinh(xi)
        ar = cz0r*cr - cz0i*cci
        ai = cz0i*cr + cz0r*cci
        cy0r = ar*cxxr - ai*cxxi
        cy0i = ai*cxxr + ar*cxxi
        cy1r = cy0r*czr(1) - cy0i*czi(1)
        cy1i = cy0i*czr(1) + cy0r*czi(1)
        ar = cy1r*cxxr - cy1i*cxxi
        ai = cy1i*cxxr + cy1r*cxxi
        cu1r = cy0r - ar
        cu1i = cy0i - ai
        cyr(1) = cy1r
        cyi(1) = cy1i
        !      cur(1)=cu1r
        !      cui(1)=cu1i
        yr(1) = cy1r
        yi(1) = cy1i
        ur(1) = cu1r
        ui(1) = cu1i

        do i = 2, nmax
            qi = dble(i)
            cyi1r = cyr(i - 1)
            cyi1i = cyi(i - 1)
            cyir = cyi1r*czr(i) - cyi1i*czi(i)
            cyii = cyi1i*czr(i) + cyi1r*czi(i)
            ar = cyir*cxxr - cyii*cxxi            !re [j/(xr+i*xi)]
            ai = cyii*cxxr + cyir*cxxi            !im [j/(xr+i*xi)]
            cuir = cyi1r - qi*ar
            cuii = cyi1i - qi*ai
            cyr(i) = cyir
            cyi(i) = cyii
            !         cur(i)=cuir
            !         cui(i)=cuii
            yr(i) = cyir
            yi(i) = cyii
            ur(i) = cuir
            ui(i) = cuii
        enddo
        !
        return
    end
    !=======================================================================
    subroutine rjb(x, y, u, nmax, nnmax)
        !=================
        !  x =(2\pi/\lambda)*r
        !  y ...
        !  nmax - angular momentum cutoff
        !  nnmax - angular momentum cutoff - determines numerical accuracy
        !--------/---------/---------/---------/---------/---------/---------/--
        integer nmax, nnmax, l, l1, i, i1
        real(dp) x, xx, z0, y0, y1, yi, yi1
        real(dp), intent(out) :: y(:), u(:)
        real(dp), allocatable :: z(:)

        l = nmax + nnmax
        allocate(z(l))
        xx = 1d0/x
        z(l) = 1d0/(dble(2*l + 1)*xx)
        l1 = l - 1
        do i = 1, l1
            i1 = l - i
            z(i1) = 1d0/(dble(2*i1 + 1)*xx - z(i1 + 1))
        end do
        z0 = 1d0/(xx - z(1))
        y0 = z0*dcos(x)*xx
        y1 = y0*z(1)
        u(1) = y0 - y1*xx
        y(1) = y1
        do i = 2, nmax
            yi1 = y(i - 1)
            yi = yi1*z(i)
            u(i) = yi1 - dble(i)*yi*xx
            y(i) = yi
        end do

        deallocate(z)
        return

    end

    !=======================================================================
    subroutine ryb(x, y, v, nmax)
        !=================
        !  x =(2\pi/\lambda)*r
        !  nmax - angular momentum cutoff
        !--------/---------/---------/---------/---------/---------/---------/--
        integer i, nmax, nmax1
        real(dp) :: c, s, x, x1, x2, x3, y1
        real(dp) y(:), v(:)
        !
        c = dcos(x)
        s = dsin(x)
        x1 = 1d0/x
        x2 = x1*x1
        x3 = x2*x1
        y1 = -c*x2 - s*x1
        y(1) = y1
        y(2) = (-3d0*x3 + x1)*c - 3d0*x2*s
        nmax1 = nmax - 1
        do i = 2, nmax1
            y(i + 1) = dble(2*i + 1)*x1*y(i) - y(i - 1)
        end do
        v(1) = -x1*(c + y1)
        do  i = 2, nmax
            v(i) = y(i - 1) - dble(i)*x1*y(i)
        end do
        return
    end
    !=======================================================================

    !    function zartan(zs)
    !!        *--------/---------/---------/---------/---------/---------/---------/--
    !!       *For a complex argument using that
    !!       *             atan (z)= -ln((1+iz)/(1-iz))*i/2.d0
    !!       *Hower, a direct application of this formula often
    !!       *results in a nonzero imaginary part of  $\arctan z$ even
    !!       *for a  purely real $z=x$.
    !!       *Therefore, in order to avoid this, in a numerical
    !!       *implementation, the complex logarithm is written
    !!       *explicitly as
    !!        *
    !!       * ln\left(\fr{1+iz}{1-iz}\right)= (ar1,ar2)= log(
    !!                       * (1.d0 + abs(z)**2 -2.d0*imag(z))/1.d0+abs(z)**2 +2.d0*imag(z)))/2.d0
    !!       *+
    !!       *ci*atan(2.d0*dble(z)/(1-abs(z)**2).
    !!        *
    !!       *For a real z=x, |x|<1,  log here is purely imaginary and
    !!       *equals to i*atan(2x/(1-x**2))\equiv 2*i*atan x
    !!        *--------/---------/---------/---------/---------/---------/---------/--
    !        complex(dp) zartan
    !        complex(dp), intent(in) :: zs
    !        real(dp) xxs,xas,ar1,ar2
    !
    !        if (aimag(zs).eq.0.) then
    !            zartan=cmplx_dp(atan(dble(zs)),0.d0)
    !            return
    !        end if
    !        xas=abs(zs)**2
    !        ar1=log((1.d0+xas-2.d0*aimag(zs))/(1.d0+xas+2.d0*aimag(zs)))/2.d0
    !        xxs=dble(zs)
    !!      *special case:
    !        if(xas.eq.1.) then
    !            if(xxs.ge.0.) then
    !                ar2=pi/2.d0
    !            else if (xxs.lt.0.) then
    !                ar2=-pi/2.d0
    !            end if
    !            zartan =cmplx_dp(ar2,- ar1)/2.d0
    !        end if
    !!      *remaining cases:
    !        ar2=2.d0*xxs/(1.d0-xas)
    !        if(xas.lt.1.d0)  then     ! 1st and 4th quadrant
    !            zartan=cmplx_dp(atan(ar2),- ar1)/2.d0
    !        else if (xas.gt.1. .and. xxs.ge.0.) then       ! 2nd quadrant
    !            zartan=cmplx_dp(pi+atan(ar2),- ar1)/2.d0
    !        else if(xas.gt.1. .and. xxs.lt.0.) then        ! 3rd quadrant
    !            zartan=cmplx_dp(-pi+atan(ar2),- ar1)/2.d0
    !        end if
    !        return
    !    end function zartan

    !=======================================================================

    subroutine vigf(x, lmax, m, dv1, dv2, ddv1)
        !--------/---------/---------/---------/---------/---------/---------/--
        ! >>> x,lmax,m (only nonnegative)
        ! <<< dv1,dv2,ddv1
        ! =============
        !
        !     x=cos(theta), where theta is the polar angle
        !     lmaxd ... maximal angular momentum cutoff
        !     lmax ... floating  angular momentum cutoff
        !
        ! returns \pi and \tau scattering functions in terms of the
        ! wigner d-functions. algorithm as described in eqs. (31-35)
        !  of ref. \cite{mis39} used. (note however a missing $n$
        !  factor in the 2nd term in the curly bracket in
        !   eq. (35) of ref. \cite{mis39}.)
        !
        !     for a given azimuthal number m.ge.0 returns
        !      the wigner d-functions
        !            dv1(n)=dvig(0,m,n,arccos x) = d_{0m}^{(l)}
        !
        !  \pi scattering function:
        !     ddv1(n)=dvig(0,m,n,arccos x)/sin(arccos x)
        !                              = m*d_{0m}^{(l)}/ sin\theta
        !
        !  \tau scattering function:
        !     dv2(n)=[d/d(arccos x)] dvig(0,m,n,arccos x)
        !                              = d d_{0m}^{(l)}/d\theta
        !
        !     for 1.le.n.le.lmax and 0.le.x.le.1
        !      ddv1 is calculated because (dv1/sin\theta) is singular for
        !             either \beta=0 or \beta=\pi
        !     (for a given m.neq.0, only the m.le.n.le.lmax terms are determined!)
        ! =====
        !
        !     in the present case (eq. (b.28) of ref. \ct{mtl})
        !                       (cf. (4.1.24) of \ct{ed}):
        !
        !             d_{00}^{(l)}(\theta)= p_l(\cos\theta)
        !      d_{0m}^{(l)}(\theta)= (-1)^m \sqrt{(l-m)!/(l+m)!}p_l^m(\cos\theta)
        !     (assuming jackson's p^m_l), where $d^{(l)}_{m0}=(-1)^m d^{(l)}_{0m}$)
        !     (edmonds p^m_l is without (-1)**m prefactor; cf. (4.1.24) therein)
        !
        ! special values:
        !
        !     (rodrigues formula [eq. (2.5.14) of ref. \ct{ed}] then yields
        !                       p_1(x)=x; p_2=(3x^2-1)/2; etc.
        !
        !      p_l^m(x) = (-1)^m (1-x^2)^{m/2} \fr{d^m p_l(x)}{dx} ===>
        !                       p_1^1(\cos\beta)=-\sin\beta
        !     therefore,
        !              d_{00}^{(1)}(\beta)=\cos\beta
        !              d_{01}^{(1)}(\beta)=\sin\beta/\sqrt{2}
        !         d d_{00}^{(1)}(\beta)/d\beta=-\sin\beta
        !         d d_{01}^{(1)}(\beta)/d\beta=\cos\beta/\sqrt{2}
        !
        !     acc. eq. (34) of {mis39}:
        !
        !     a_0=1, a_1=1/\sqrt{2}, a_2=\sqrt{3}/(2*\sqrt{2})
        !
        !     therefore [eq. (32) of {mis39}]:
        !              d_{00}^{(0)}(\beta)=1
        !              d_{01}^{(1)}(\beta)=\sin\beta/\sqrt{2}
        !              d_{02}^{(2)}(\beta)=\sqrt{3}\sin^2\beta/(2*\sqrt{2})
        !     and
        !         d d_{00}^{(0)}(\beta)/d\beta=0
        !         d d_{01}^{(1)}(\beta)/d\beta=\cos\beta/\sqrt{2}
        !         d d_{02}^{(2)}(\beta)/d\beta=\sqrt{3}\sin\beta \cos\beta/\sqrt{2}
        !                                = \sqrt{3}\sin (2\beta) /(2*\sqrt{2})
        ! =====
        !     similar to routine vig, which however only returns
        !            dv1(n)=dvig(0,m,n,arccos x) = d_{0m}^{(l)}
        !
        !     when arccos x is very small, a care has to be exercise to generate
        !     nboth ddv1 and dv2. that part has been made using recurrences of
        !     ref. \ct{tks}
        !--------/---------/---------/---------/---------/---------/---------/--
        implicit none

        integer n, lmax, m, i, i2
        real(dp) a, x, qs, d1, d2, d3, der, dx, qn, qn1, qn2, &
                qnm, qnm1, qmm
        real(dp) ddv1(:), dv1(:), dv2(:)

        !* dv1, ddv1, and dv2 initialization
        do n = 1, lmax
            dv1(n) = 0.d0
            ddv1(n) = 0.d0
            dv2(n) = 0.d0
        end do

        dx = dabs(x)
        a = 1.d0                       !initial a_0
        qs = dsqrt(1d0 - x*x)            !sin\theta
        !***********************************************************************
        !*                      nonzero dv1 initialization
        !*
        if (m==0) then
            !*
            !* ddv1(n)=0.d0 [see (3.33) of {tks}]
            !* d1,d2, and d3 below are the three consequent terms
            !*       d_{0m}^{n-1}, d_{0m}^{n}, and d_{0m}^{n+1} beginning
            !*       with n=m
            !*=============
            !* recurrence initialization following d^l_{00}=p_l
            !*         [see eq. (b.27) of {mtl}]
            !*
            d1 = 1.d0                  !d^0_{00}=p_0=1   (see sec. 6.8 of {nr})
            d2 = x                     !d^1_{00}=p_1=x

            do n = 1, lmax
                qn = dble(n)
                qn1 = dble(n + 1)
                qn2 = dble(2*n + 1)
                !       *recurrence (31) of ref. {mis39} for d^{n+1}_{00}
                d3 = (qn2*x*d2 - qn*d1)/qn1
                dv1(n) = d2                     !d^n_{00}
                d1 = d2                         !becomes d^{n-1}_{00} in d3
                d2 = d3                         !becomes d^{n}_{00} in d3
            end do
        else
            !***********************************************************
            !*                           m\neq 0 part
            !*  d_{0m}^m \equiv a_m*(sin\theta)**m   initialization - (33) and recurrence (34) of ref. {mis39}
            do i = 1, m
                i2 = i*2
                !          !recurrence (33,34) of ref. {mis39} f
                a = a*dsqrt(dble(i2 - 1)/dble(i2))*qs
            end do
            !* recurrence initialization following eqs. (32,33) of ref. {mis39}

            d1 = 0.d0                 !=dv1(m-1); see eq. (32) of ref. {mis39}
            d2 = a                    !=dv1(m);   see eq. (33) of ref. {mis39}
            qmm = dble(m*m)
            !*
            do n = m, lmax
                qn = dble(n)
                qn1 = dble(n + 1)
                qn2 = dble(2*n + 1)
                qnm = dsqrt(qn*qn - qmm)
                qnm1 = dsqrt(qn1*qn1 - qmm)
                !          !recurrence (31) of ref. {mis39} for d^{n+1}_{0m}
                d3 = (qn2*x*d2 - qnm*d1)/qnm1
                dv1(n) = d2                             !d^n_{0m}
                d1 = d2                           !becomes d^{n-1}_{0m} in d3,der
                d2 = d3                                 !becomes d^{n}_{0m} in d3
            end do
        end if

        !*                        dv1 initialized
        !*             it remains to determine ddv1 and dv2
        !*********************************************************
        !*  (1-cos\theta) is very small:
        !*
        !   for theta=0 [see eqs. above]:
        !              d_{00}^{(0)}(0)=1
        !              d_{01}^{(1)}(0)=0
        !              d_{02}^{(2)}(\beta)=0
        !     and
        !         d d_{00}^{(0)}(\beta)/d\beta=0
        !         d d_{01}^{(1)}(\beta)/d\beta=1/\sqrt{2}
        !         d d_{02}^{(2)}(\beta)/d\beta=0
        !
        !  see eqs. (4.1-4) of \ct{mis91}:
        !
        !   (m/\sin\theta) d_{0m}^l(0)=(\delta_{m\pm 1}/2) \sqrt{l(l+1)}
        !      d d_{0m}^l(0)/d\beta   =(m\delta_{m\pm 1}/2) \sqrt{l(l+1)}
        !
        !*
        !*  (4.2.1) of \ct{ed}:
        !*   d_{0m}^{(l)}(pi) = (-1)^{l+m} \dt_{0,m}
        !*
        !*  (4.2.3) of \ct{ed}:
        !*   d_{0m}^{(l)}(0) = (-1)^{m} \dt_{0,m} = \dt_{0,m}
        !*=======================================
        !*
        !*  if x^l_m=(m/\sin\theta) d_{0m}^{(l)}, then, according to (3.29) of {tks}:
        !*
        !*  x^{m+1}_{m+1}=\sin\theta \sqrt{\fr{2m+1}{2m+2}}
        !*                           \left(\fr{m+1}{m}\right)x^{m}_{m}
        !*
        !*  according to (3.30) of {tks}:
        !*  x^{m+1}_{m}= -\sqrt{2m+1}\,\cos\theta x^{m}_{m}
        !*
        !* according to (3.31) of {tks}:
        !*  x^{l}_{m}=\fr{1}{\sqrt{l^2-m^2}}\,\left[(2l-1)\cos\theta
        !*          x^{l-1}_{m} - \sqrt{(l-1)^2-m^2}}\,\x^{l-2}_{m} \right]
        !*
        !* initial recurrence values are x^1_1=\sqrt{2}/2 and x^l_0=0
        !***********************************************************************
        !*                   nonzero ddv1/dv2 initialization
        !*
        !*                          m = 0

        if (m==0) then     !all ddv1(n)=x^l_0=0; see (3.33) of {tks}:
            !* according to (3.37) of {tks}, dv2(0)=0.d0

            dv2(1) = -qs

            if (lmax>=2) dv2(2) = 3*x*dv2(1)

            if (lmax<3) return
            !*
            do n = 3, lmax           !recurrence (3.36) of {tks},
                dv2(n) = (2*n - 1)*x*dv2(n - 1)/(n - 1) - n*dv2(n - 2)/(n - 1)
            enddo
            !***********************************************************************
            !*                           m > 0

        else if (m>0) then
            !*
            !* >>> determine x^m_m according to eq. (3.29) of {tks}:

            a = 1.d0/dsqrt(2.d0)               !x^1_1=a_1

            do i = 1, m - 1
                a = qs*dble(i + 1)*dsqrt(2*i + 1.d0)*a/(i*dsqrt(2*i + 2.d0))
            enddo
            !* <<< a is now x^m_m; see (3.29) of {tks}

            ddv1(m) = a
            dv2(m) = x*a                        !see (3.34) of {tks}
            !* >>> determine x^{m+1}_m:

            if (m==lmax)  return
            !     ! der=x^{m+1}_m; see (3.30) of {tks}
            der = x*dsqrt(2*m + 1.d0)*a
            ddv1(m + 1) = der
            dv2(m + 1) = ((m + 1)*x*der - a*dsqrt(2*m + 1.d0))/dble(m)  !(3.35) of {tks}
            !* >>> determine remaining x^{l}_m's

            if ((m + 2)==lmax) return

            do n = m + 2, lmax
                d3 = dsqrt(dble(n)**2 - dble(m)**2)
                ddv1(n) = ((2*n - 1)*x*ddv1(n - 1) - &
                        dsqrt(dble(n - 1)**2 - dble(m)**2)*ddv1(n - 2))/d3
                !                                                 !see (3.31) of {tks}
                dv2(n) = (n*x*ddv1(n) - ddv1(n - 1)*d3)/dble(m)  !see (3.35) of {tks}
            enddo

        end if

        return
    end

    !=======================================================================

    subroutine vig (x, nmax, m, dv1, dv2)
        !--------/---------/---------/---------/---------/---------/---------/--
        ! >>> x,nmax,m (only nonnegative)
        ! <<< dv1, dv2
        ! =============
        !     for a given azimuthal number m.ge.0 returns
        !      the wigner d-functions , i.e.,
        !
        !     dv1(n)=dvig(0,m,n,arccos x) = d_{0m}^{(l)}
        !     and
        !     dv2(n)=[d/d(arccos x)] dvig(0,m,n,arccos x)  ! = d d_{0m}^{(l)}/d\theta
        !
        !     for 1.le.n.le.nmax and 0.le.x.le.1
        !     (for a given m.neq.0, only the m.le.n.le.nmax terms are determined!)
        !     according to eq. (4.1.24) of ref. \ct{ed}:
        !
        !             d_{00}^{(l)}(\theta)= p_l(\cos\theta) ===>
        !
        !     (rodrigues formula [eq. (2.5.14) of ref. \ct{ed}] then yields
        !                       p_1(x)=x; p_2=(3x^2-1)/2; etc.
        !     one can show that $d_{00}^{(1)}(\theta)=\cos\theta$
        !
        !     similar to routine vigampl, which however returns the wigner d-functions
        !     divided by sin\theta, i.e.,
        !     dv1(n)=dvig(0,m,n,arccos x)/sin(arccos x) = d_{0m}^{(l)}/sin\theta
        !
        !     made using recurrences of  ref. \ct{mis39}
        !     (there is a missing $l$ factor in the 2nd term in the curly bracket
        !     in recurrence (35) of ref. \ct{mis39} for dv2).
        !
        !     one has (see eq. (4.2.5) of \ct{ed}):
        !                       $d_{0m}^{(l)}=(-1)^m d_{0-m}^{(l)}$
        !     and (see eq. (35) of \ct{mis91}):
        !            $dd_{0m}^{(l)}/(d\theta)=(-1)^m dd_{0-m}^{(l)}/(d\theta)$
        !
        !     x=cos(theta), where theta is the polar angle
        !     nmax - angular momentum cutoff
        !
        !     called by tmatr and tmatr0 routines
        !--------/---------/---------/---------/---------/---------/---------/--
        implicit none
        real(dp) dv1(:), dv2(:)
        integer n, nmax, m, i, i2
        real(dp) a, x, qs, qs1, d1, d2, d3, der, qn, qn1, qn2, &
                qnm, qnm1, qmm

        a = 1d0
        qs = dsqrt(1d0 - x*x)
        qs1 = 1d0/qs
        do n = 1, nmax
            dv1(n) = 0d0
            dv2(n) = 0d0
        enddo

        if (m==0) then
            d1 = 1d0
            d2 = x
            do n = 1, nmax
                qn = dble(n)
                qn1 = dble(n + 1)
                qn2 = dble(2*n + 1)
                !        !recurrence (31) of ref. {mis39}
                d3 = (qn2*x*d2 - qn*d1)/qn1
                !        !recurrence (35) of ref. {mis39}
                der = qs1*(qn1*qn/qn2)*(-d1 + d3)
                dv1(n) = d2
                dv2(n) = der
                d1 = d2
                d2 = d3
            enddo
        else
            qmm = dble(m*m)
            !*a_m initialization - recurrence (34) of ref. {mis39}
            do i = 1, m
                i2 = i*2
                a = a*dsqrt(dble(i2 - 1)/dble(i2))*qs
            enddo
            !*
            d1 = 0d0
            d2 = a

            do n = m, nmax
                qn = dble(n)
                qn2 = dble(2*n + 1)
                qn1 = dble(n + 1)
                qnm = dsqrt(qn*qn - qmm)
                qnm1 = dsqrt(qn1*qn1 - qmm)
                !        !recurrence (31) of ref. {mis39}
                d3 = (qn2*x*d2 - qnm*d1)/qnm1
                !        !recurrence (35) of ref. {mis39}
                der = qs1*(-qn1*qnm*d1 + qn*qnm1*d3)/qn2
                dv1(n) = d2
                dv2(n) = der
                d1 = d2
                d2 = d3
            enddo
        end if
        return
    end

    !=======================================================================

    subroutine vig_1v (x, nmax, m, dv1, dv2)
        !--------/---------/---------/---------/---------/---------/---------/--
        ! Similar to vig(), but return values just for nmax
        !--------/---------/---------/---------/---------/---------/---------/--
        implicit none
        real(dp) dv1, dv2
        integer n, nmax, m, i, i2
        real(dp) a, x, qs, qs1, d1, d2, d3, der, qn, qn1, qn2, &
                qnm, qnm1, qmm

        a = 1d0
        qs = dsqrt(1d0 - x*x)
        qs1 = 1d0/qs
        dv1 = 0d0
        dv2 = 0d0

        if (m==0) then

            d1 = 1d0
            d2 = x
            do n = 1, nmax
                qn = dble(n)
                qn1 = dble(n + 1)
                qn2 = dble(2*n + 1)
                !        !recurrence (31) of ref. {mis39}
                d3 = (qn2*x*d2 - qn*d1)/qn1
                !        !recurrence (35) of ref. {mis39}
                der = qs1*(qn1*qn/qn2)*(-d1 + d3)
                dv1 = d2
                dv2 = der
                d1 = d2
                d2 = d3
            enddo
        else
            qmm = dble(m*m)
            !*a_m initialization - recurrence (34) of ref. {mis39}
            do i = 1, m
                i2 = i*2
                a = a*dsqrt(dble(i2 - 1)/dble(i2))*qs
            enddo
            !*
            d1 = 0d0
            d2 = a

            do n = m, nmax
                qn = dble(n)
                qn2 = dble(2*n + 1)
                qn1 = dble(n + 1)
                qnm = dsqrt(qn*qn - qmm)
                qnm1 = dsqrt(qn1*qn1 - qmm)
                !        !recurrence (31) of ref. {mis39}
                d3 = (qn2*x*d2 - qnm*d1)/qnm1
                !        !recurrence (35) of ref. {mis39}
                der = qs1*(-qn1*qnm*d1 + qn*qnm1*d3)/qn2
                dv1 = d2
                dv2 = der
                d1 = d2
                d2 = d3
            enddo
            !*
        end if
        return
    end

end module special_functions
