subroutine tmtaxsp(nmax, rap, zeps1, tmt)

    ! warning in module tmtaxsp in file tmtaxsp.f: variables set but never used:
    !    nggg set at line 182 file tmtaxsp.f

    ! warning in module tmtaxsp in file tmtaxsp.f: variables may be used before set:
    !    qext1 used at line 215 file tmtaxsp.f
    !    qext1 set at line 220 file tmtaxsp.f
    !    qsca1 used at line 214 file tmtaxsp.f
    !    qsca1 set at line 221 file tmtaxsp.f

    !--------/---------/---------/---------/---------/---------/---------/--
    ! nmax - angular momentum cut off
    ! rap=s(1,1)*kappa0/2.d0/pi     !=rmuf*alpha/lambda =rsnm/lambda
    !
    !    returns the t matrix of a general axially symmetric scatterer
    !    the latter has also non-zero of-diagonal (mixed eh and he terms):
    !
    !                         |  tmt(1,*) |  tmt(4,*)   |
    !                 tmt  =  | ----------+-------------|
    !                         |  tmt(3,*) |  tmt(2,*)   |
    !
    !    tmt(1,*) terms corresponds to tee scattering matrix,
    !    tmt(2,*) terms corresponds to tmm scattering matrix
    !    tmt(3,*) terms corresponds to tme scattering matrix
    !    tmt(4,*) terms corresponds to tem scattering matrix
    !    tmt(4,*)=-tmt(3,*)^t where t denotes transposed tmt(3,*) submatrix
    !
    !    tmt's equal to i*sin(eta)*exp(i*eta), where eta is a phase-shift
    !             ====>    s=1+2*t=exp(2*i*eta)
    !    and the unitarity of s-matrix implies
    !            t^\dagger*t=t*t^\dagger=-(1/2)*(t^\dagger+t)
    !______________________________
    !
    ! local variables:
    ! ===============
    ! mpar%ichoice=1 if nag library is available, otherwise mpar%ichoice=2
    !
    ! np,eps: specifies the shape of particles within a given np class:
    !     np>0 - eps = deformation parameter of a chebyshev particle
    !     np=-1 - eps = the ratio of the horizontal to rotational axes. eps is
    !             larger than 1 for oblate spheroids and smaller than 1 for
    !             prolate spheroids.
    !     np=-2 - eps = the ratio of the cylinder diameter to its length.
    !     np=-3 - no eps is specified
    !     np=-4 - eps is the height (along the axial symmetry axis)
    !             of the resulting cut sphere
    !                note that always eps.lt.2*rev specified
    !     np=-9 - eps = the ratio of the cylinder diameter to its length,
    !             nanorod caps are set independantly.
    !
    ! warning:
    !  in computations for spheres, use eps=1.000001 instead of eps=1.
    !  eps=1 can cause overflows in some rare cases.
    !
    !  lam - the (vacuum) wavelength of incident light. changed to
    !                       lam=lam*sqrt(mpar%zeps0) here
    !
    !  rat = 1 - particle size is specified in terms of the
    !                equal-volume-sphere radius
    !  rat.ne.1 - particle size is specified in terms of the
    !                equal-surface-area-sphere radius
    !  axi ... equivalent-(volume/surface-area)-sphere radius
    !  rev=a=rat*axi ... equal-volume-sphere radius
    !                  (feeded as rev to rsp* routines)
    !  ddelt - required precision
    !
    !  alpha and beta - euler angles (in degrees) specifying the
    !          orientation of the scattering particle relative to
    !          the laboratory reference frame (refs. 6 and 7).
    !
    !  for axially symmetric scatterers, when the t matrix is computed in
    !  natural coordinate system with the $z$ axis along the axis of particle
    !  axial symmetry, one can show that the t matrix is {\em diagonal} with
    !  respect to the azimuthal indices $m$ and $m'$ \cite{wat},
    !
    !              t_{lm,l'm'}^{ij}=\delta_{mm'} t_{lm,l'm},
    !
    !  and that it satisfies reciprocity relation \cite{gus,mis36},
    !
    !               t_{lm,l'm}^{ij}=(-1)^{i+j} t_{l'm,lm}^{ji}.
    !
    !  \cite{mis91} also claims the relation:
    !
    !                t_{lm,l'm}^{ij}= (-1)^{i+j} t_{l-m,l'-m}^{ij}
    !
    !--------/---------/---------/---------/---------/---------/---------/--
    use libcylinder
    implicit none

    integer, parameter :: lmaxd = 8, lmax1d = lmaxd + 1, lmtd = lmax1d*lmax1d - 1

    integer ndgs, ngauss, n, nm, np, naxsm
    integer n1, n2, m, l1, l2, lmtot, k1, k2, kk1, kk2
    integer nmax, ncheck, inm1, ixxx, ja, jam, jb, jbm
    real(dp) rev, eps, ht, a, alpha, beta, ddelt, dn1
    real(dp) p, pir, pii, ppi, qext, qsca, rat, rap, rsnm
    real(dp) tr1nn, ti1nn, tr1nn1, ti1nn1, xev
    real(dp)  lam, mrr, mri, x(npng2), w(npng2), s(npng2), ss(npng2), &
            an(npn1), r(npng2), dr(npng2), &
            ddr(npng2), drr(npng2), dri(npng2), ann(npn1, npn1)
    real(dp) tr1(npn2, npn2), ti1(npn2, npn2)
    !real(dp) xalpha(300),xbeta(300),walpha(300),wbeta(300)

    !complex(dp) czero
    complex(dp) zeps1
    complex(dp) tmt(4, lmtd, lmtd)

    common /ct/ tr1, ti1
    ! transfers the real and imaginary part of the t matrix (2*nmax,2*nmax)
    ! array for a given value of m from tt via (tmatr0 and tmatr) to the main

    !common /tmat/ rt11,rt12,rt21,rt22,it11,it12,it21,it22
    ! transfers t matrix arrays obtained from tr1,ti1 in the main to the
    ! ampl routine --->  not used here

    !common /choice/ ichoice
    ! transfers the choice of inversion from here to tt

    common/revval/ a
    ! transfers rev to const routine
    common /toitmt/np, ncheck, naxsm, ndgs
    ! transfers integers np,ncheck,naxsm,ndgs from the main here

    common /totmt/eps, rat, rev, alpha, beta, ddelt
    ! transfers real(dp) rat,a(rev),alpha,beta,ddelt from the main here

    !****************************************************************
    !data czero/(0.d0,0.d0)/

    p = dacos(-1d0)                 !local pi constant

    a = rev
    lam = rev*sqrt(mpar%zeps0)/rap       !vacuum wavelength times sqrt(mpar%zeps0)/

    ! the real part of the refractive index contrast
    mrr = dble(sqrt(zeps1/mpar%zeps0))

    ! the imaginary  part of the refractive index contrast
    mri = aimag(sqrt(zeps1/mpar%zeps0))

    ddelt = 0.1d0*ddelt               !conv. test is switched off now!!!

    ! ddelt is used to test the accuracy of computing the
    ! optical cross sections. this accuracy is usually better
    ! than the absolute accuracy of computing the expansion coefficients
    ! of a normalized scattering matrix by a factor of 10. therefore,
    ! the desired accuracy of computing the expansion coefficients
    ! is rescaled by a factor 0.1 before entering the test of the
    ! accuracy of computing the optical cross sections.
    !
    ! other local constants:
    !
    lmtot = (nmax + 1)**2 - 1

    ! Check if the radius is area-equivalent instead of volume-equivalent
    if (dabs(rat - 1d0) > 1d-8) then
        select case (np)
        case(0:)
!            rat = radii_ratio_chebyshev(np)
        case (-1) ! oblate/prolate spheroids
            rat = radii_ratio_spheroid()
        case (-2) ! oblate/prolate cylinder
            rat = radii_ratio_cylinder()
        case (-3) ! distorted chebyshev droplet
            rat = radii_ratio_droplet()
        case (-9) ! nanorod
            ! todo make a correct evaluation of surface-equivalent radius for nanorod
            ! (2019-12-09) In principle it is done but needs some testing
            rat = radii_ratio_nanorod()
        end select
    end if

    !___________________________________________________
    ! determination of the wiscombe value of the floating
    ! angular momentum cutoff nmax:

    xev = 2d0*p*a/lam
    ixxx = xev + 4.05d0*xev**0.333333d0     !wiscombe conv. criterion for max
    inm1 = max0(4, ixxx)
    !
    if (inm1 >= npn1) then
        print 7333, npn1
        stop
    end if
    7333 format('convergence is not obtained for npn1=', i3, &
            '.  execution terminated')
    !_______________________________________________________________

    ngauss = nmax*ndgs
    !c      nnnggg=ngauss+1

    if (ngauss == npng1) print 7336
    7336    format('warning: ngauss=npng1')
    !
    ! gif division points and weights + other numerical constants
    !
    call const(ngauss, nmax, x, w, an, ann, s, ss, np, eps)
    !
    ! specify particle shape:
    !
    call vary(lam, mrr, mri, a, eps, &
            rsnm, ht, np, ngauss, x, p, ppi, pir, pii, r, &
            dr, ddr, drr, dri, nmax)

    !         subroutine vary (lam,mrr,mri,a,eps,mpar%nanorod_cap_hr,
    !    &                 rsnm,ht,
    !                   np,ngauss,x,p,ppi,pir,pii,r,
    !                   dr,ddr,drr,dri,nmax)
    !
    !
    ! determine m=m'=0 elements of the t matrix
    !
    if (mpar%integration_type == 0) then
        call tmatr0 (ngauss, x, w, an, ann, ppi, pir, pii, r, dr, &
                ddr, drr, dri, nmax, ncheck, naxsm)
    else if (mpar%integration_type == 1) then
        call tmatr0_adapt (ngauss, x, w, an, ann, ppi, pir, pii, r, dr, &
                ddr, drr, dri, nmax, ncheck, naxsm)
    else
        call tmatr0_leru (ngauss, x, w, an, ann, ppi, pir, pii, r, dr, &
                ddr, drr, dri, nmax, ncheck, naxsm)
    endif

    qext = 0d0
    qsca = 0d0

    do n = 1, nmax
        n1 = n + nmax
        tr1nn = tr1(n, n)
        ti1nn = ti1(n, n)
        tr1nn1 = tr1(n1, n1)
        ti1nn1 = ti1(n1, n1)

        dn1 = dble(2*n + 1)

        qsca = qsca + dn1*(tr1nn*tr1nn + ti1nn*ti1nn&
                + tr1nn1*tr1nn1 + ti1nn1*ti1nn1)
        qext = qext + (tr1nn + tr1nn1)*dn1
    end do
    !<<<
    !c      write(nout,*)'nmax=',nmax
    !c      write(nout,*)'ngauss=',ngauss
    !<<<
    !
    ! tmt initialization:

    !c         do ja=1,lmtot
    !c         do jb=1,lmtot
    !
    !c            tmt(1,jb,ja)=czero
    !c            tmt(2,jb,ja)=czero
    !c            tmt(3,jb,ja)=czero
    !c            tmt(4,jb,ja)=czero

    !c         enddo
    !c         enddo
    !
    !
    !                         |  tmt(1,*) |  tmt(4,*)   |
    !                 tmt  =  | ----------+-------------|
    !                         |  tmt(3,*) |  tmt(2,*)   |
    !
    !    tmt(1,*) terms corresponds to tee scattering sub-matrix
    !    tmt(4,*) terms corresponds to tem scattering sub-matrix
    !    tmt(2,*) terms corresponds to tmm scattering sub-matrix
    !    tmt(3,*) terms corresponds to tme scattering sub-matrix
    !    tmt(4,*)=-tmt(3,*)^t where t denotes transposed tmt(3,*) submatrix
    !****************     assign  m=m=0 elements of tmt matrix   ***********

    do l1 = 1, nmax
        do l2 = 1, nmax
            n1 = l1 + nmax
            n2 = l2 + nmax

            ja = l1*(l1 + 1)       ! (l,m) index with (1-1)=1
            jb = l2*(l2 + 1)
            !
            ! see (5.39) of {mtl}:  !!!iff plane of symmetry perpendicular to the
            !                       !!!axis of rotation

            if ((naxsm == 1) .and. ((-1)**(l1 + l2) /= 1)) then
                tmt(2, ja, jb) = czero
                tmt(1, ja, jb) = czero
            else
                tmt(2, ja, jb) = cmplx_dp(tr1(l1, l2), ti1(l1, l2))
                tmt(1, ja, jb) = cmplx_dp(tr1(n1, n2), ti1(n1, n2))
            end if
            ! see (5.37) of {mtl}:
            tmt(4, ja, jb) = czero         !cmplx_dp(tr1(n1,l2),ti1(n1,l2))
            tmt(3, jb, ja) = czero         !-tmt(4,ja,jb)

        enddo
    enddo
    !
    !****************    assign  m=m'>0 elements of the t matrix   ***********

    do m = 1, nmax
        !
        if (mpar%integration_type == 0) then
            call tmatr(m, ngauss, x, w, an, ann, s, ss, ppi, pir, pii, r, dr, &
                    ddr, drr, dri, nmax, ncheck, naxsm)
        else if (mpar%integration_type == 1) then
            call tmatr_adapt(m, ngauss, x, w, an, ann, s, ss, ppi, pir, pii, r, dr, &
                    ddr, drr, dri, nmax, ncheck, naxsm)
        else
            call tmatr_leru(m, ngauss, x, w, an, ann, s, ss, ppi, pir, pii, r, dr, &
                    ddr, drr, dri, nmax, ncheck, naxsm)
        endif
        !
        ! <<< returns  m=m'>0 elements of the t matrix
        !
        nm = nmax - m + 1             !size of a t block returned by tt

        do l1 = m, nmax
            do l2 = m, nmax

                k1 = l1 - m + 1               !k1,k2,kk1,kk2 label the entries of
                k2 = l2 - m + 1               !a tt returned t block
                kk1 = l1 - m + 1 + nm
                kk2 = l2 - m + 1 + nm

                jam = l1*(l1 + 1) - m       ! (l,m) index with (1-1)=1
                jbm = l2*(l2 + 1) - m
                ja = l1*(l1 + 1) + m       ! (l,m) index with (1-1)=1
                jb = l2*(l2 + 1) + m
                !
                !
                ! see (5.39) of {mtl}: !!!iff plane of symmetry perpendicular to the
                !                        !!!axis of rotation

                if ((naxsm == 1) .and. ((-1)**(l1 + l2) /= 1)) then
                    tmt(2, ja, jb) = czero
                    tmt(2, jam, jbm) = czero
                    tmt(1, ja, jb) = czero
                    tmt(1, jam, jbm) = czero
                else
                    tmt(2, ja, jb) = cmplx_dp(tr1(k1, k2), ti1(k1, k2))
                    tmt(2, jam, jbm) = tmt(2, ja, jb)
                    tmt(1, ja, jb) = cmplx_dp(tr1(kk1, kk2), ti1(kk1, kk2))
                    tmt(1, jam, jbm) = tmt(1, ja, jb)
                end if

                if ((naxsm == 1) .and. ((-1)**(l1 + l2) /= -1)) then
                    tmt(4, ja, jb) = czero
                    tmt(4, jam, jbm) = czero
                    tmt(3, ja, jb) = czero
                    tmt(3, jam, jbm) = czero
                else
                    tmt(4, ja, jb) = cmplx_dp(tr1(kk1, k2), ti1(kk1, k2))
                    tmt(4, jam, jbm) = -tmt(4, ja, jb)
                    tmt(3, jb, ja) = -tmt(4, ja, jb)
                    tmt(3, jbm, jam) = tmt(4, ja, jb)    !=-tmt(3,jb,ja)
                end if
                !
                !  using reciprocity (eq. (15) of ref. \ct{mis97}):
                !
                !       t_{lm,l'm}^{ij}=(-1)^{i+j} t_{l'm,lm}^{ji}
                !
                !  and (see eq. (36) \jqsrt{55}):
                !
                !       t_{lm,l'm'}^{ij}=(-1)^{m+m'} t_{l'-m',l-m}^{ji}
                !
                !  moreover, for axially symmetric particles one has
                !  (see eq. (31),(36) \jqsrt{55}):
                !
                !  t_{lm,l'm}^{ij}=(-1)^{i+j} t_{l-m,l'-m}^{ij} = t_{l'-m',l-m}^{ji}
                !
            enddo
        enddo
        !
    end do    !end of loop over m's

    return
end

!********************************************************************

subroutine vigampl(x, nmax, m, ddv1, dv2)
    !--------/---------/---------/---------/---------/---------/---------/--
    ! >>> x,nmax,m (only nonnegative)
    ! <<< dv1, dv2
    ! =============
    !     for a given azimuthal number m.ge.0 returns
    !      the wigner d-functions divided by sin\theta, i.e.,
    !
    !     ddv1(n)=dvig(0,m,n,arccos x)/sin(arccos x)   ! = m*d_{0m}^{(l)}/ sin\theta
    !     and
    !     dv2(n)=[d/d(arccos x)] dvig(0,m,n,arccos x)  ! = d d_{0m}^{(l)}/d\theta
    !
    !     for 1 <= n <= nmax and 0 <= x <= 1
    !     (for a given m.neq.0, only the m <= n <= nmax terms are determined!)
    !     according to eq. (4.1.24) of ref. \ct{ed}:
    !
    !             d_{00}^{(l)}(\theta)= p_l(\cos\theta)
    !
    !     (rodrigues formula [eq. (2.5.14) of ref. \ct{ed}] then yields
    !                       p_1(x)=x; p_2=(3x^2-1)/2; etc.
    !     one can show that $d_{00}^{(1)}(\theta)=\cos\theta$
    !
    !     similar to routine vig, which however returns
    !     dv1(n)=dvig(0,m,n,arccos x)   ! = d_{0m}^{(l)}
    !
    !     in addition, vigampl has a block treating the case when
    !     arccos x is very small option
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
    !     lmaxd ... maximal angular momentum cutoff
    !     nmax ... floating  angular momentum cutoff
    !     called only by the ampl routine!!!
    !--------/---------/---------/---------/---------/---------/---------/--
    use libcylinder
    implicit none
    real(dp) ddv1(npn1), dv2(npn1)

    !      implicit none
    !      integer lmaxd,lmaxd1
    !      parameter (lmaxd=50,lmaxd1=lmaxd+1)
    !      real(dp) ddv1(lmaxd1),dv1(lmaxd1),dv2(lmaxd1)

    integer n, nmax, m, i, i2
    real(dp) a, x, qs, d3, der, dx, qmm
    ! ddv1 and dv2 initialization
    do n = 1, nmax
        ddv1(n) = 0.d0
        dv2(n) = 0.d0
    end do

    dx = dabs(x)
    a = 1.d0
    qs = dsqrt(1d0 - x*x)                        !sin\theta

    if (m /= 0) then          ! m\neq 0 part

        !    a_m*(sin\theta)**m   initialization - (33) and recurrence (34) of ref. {mis39}

        qmm = dble(m*m)

        do i = 1, m
            i2 = i*2
            a = a*dsqrt(dble(i2 - 1)/dble(i2))*qs  !recurrence (33,34) of ref. {mis39} f
        end do

    end if

    !********************************************************
    !  (1-cos\theta) is very small:
    !
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
    !
    !  (4.2.1) of \ct{ed}:
    !   d_{0m}^{(l)}(pi) = (-1)^{l+m} \dt_{0,m}
    !
    !  (4.2.3) of \ct{ed}:
    !   d_{0m}^{(l)}(0) = (-1)^{m} \dt_{0,m} = \dt_{0,m}
    !=======================================
    !
    !  if x^l_m=(m/\sin\theta) d_{0m}^{(l)}, then, according to (3.29) of {tks}:
    !
    !  x^{m+1}_{m+1}=\sin\theta \sqrt{\fr{2m+1}{2m+2}}
    !                           \left(\fr{m+1}{m}\right)x^{m}_{m}
    !
    !  according to (3.30) of {tks}:
    !  x^{m+1}_{m}= -\sqrt{2m+1}\,\cos\theta x^{m}_{m}
    !
    ! according to (3.31) of {tks}:
    !  x^{l}_{m}=\fr{1}{\sqrt{l^2-m^2}}\,\left[(2l-1)\cos\theta
    !          x^{l-1}_{m} - \sqrt{(l-1)^2-m^2}}\,\x^{l-2}_{m} \right]
    !
    ! initial recurrence values are x^1_1=\sqrt{2}/2 and x^l_0=0
    !**********************************************************************
    !                   nonzero ddv1/dv2 initialization
    !                          m = 0

    if (m==0) then     !all ddv1(n)=x^l_0=0; see (3.33) of {tks}:
        ! according to (3.37) of {tks}, dv2(0)=0.d0

        dv2(1) = qs

        if (nmax>=2) dv2(2) = 3*x*dv2(1)

        if (nmax<3) return
        !
        do n = 3, nmax           !recurrence (3.36) of {tks},
            dv2(n) = (2*n - 1)*x*dv2(n - 1)/(n - 1) - n*dv2(n - 2)/(n - 1)
        enddo
        !**********************************************************************
        !                           m > 0

    else if (m>0) then

        ! >>> determine x^m_m according to eq. (3.29) of {tks}:

        a = 1.d0/dsqrt(2.d0)               !x^1_1=a_1

        do i = 1, m - 1
            a = qs*dble(i + 1)*dsqrt(2*i + 1.d0)*a/(i*dsqrt(2*i + 2.d0))
        enddo
        ! <<< a is now x^m_m; see (3.29) of {tks}

        ddv1(m) = a
        dv2(m) = x*a                        !see (3.34) of {tks}
        ! >>> determine x^{m+1}_m:

        if (m/=nmax)  then
            !see (3.30) of {tks}
            der = x*dsqrt(2*m + 1.d0)*a          ! der=x^{m+1}_m;
            ddv1(m + 1) = der
            !(3.35) of {tks}
            dv2(m + 1) = ((m + 1)*x*der - a*dsqrt(2*m + 1.d0))/dble(m)
            ! >>> determine remaining x^{l}_m's
        end if

        if ((m + 2)/=nmax)  then
            do n = m + 2, nmax
                !see (3.31) of {tks}
                d3 = dsqrt(dble(n)**2 - dble(m)**2)
                ddv1(n) = ((2*n - 1)*x*ddv1(n - 1) - &
                        dsqrt(dble(n - 1)**2 - dble(m)**2)*ddv1(n - 2))/d3

                !see (3.35) of {tks}
                dv2(n) = (n*x*ddv1(n) - ddv1(n - 1)*d3)/dble(m)
            enddo
        end if

    end if

    return
end
!**********************************************************************

subroutine const(ngauss, nmax, x, w, an, ann, s, ss, np, eps)
    !--------/---------/---------/---------/---------/---------/---------/--
    ! >>> ngauss,nmax,np,eps
    ! <<< x,w,an,ann,s,ss
    !=====================
    !
    !  ngauss - the number of gif division points
    !  nmax - angular momentum cutoff
    !  p=pi=dacos(-1d0)
    !  np - parameter specifying the particle shape
    !  eps - deformation parameter for a given particle shape
    !
    !  x=\cos\theta  - gif division points
    !  w - gif weights
    !  an(n)=n*(n+1)
    !  ann(l_1,l_2)=\sqrt{\fr{(2 l_1+1)}{l_1(l_1+1)} }
    !                       \sqrt{\fr{(2 l_2+1)}{l_2(l_2+1)} }/2
    !  s  ... 1/(|\sin\theta|)
    !  ss ... 1/(\sin^2\theta)
    !--------/---------/---------/---------/---------/---------/---------/--
    use libcylinder
    implicit none
    integer, intent(in):: ngauss, nmax, np
    integer ng, ng1, ng2, nn
    integer neps, jg, i, j, n, n1
    real(dp), intent(in):: eps
    real(dp) d, ddd, xx, y
    real(dp) ee, ee1, cc, si, xi1, xi2, xav
    real(dp) xtheta, theta0, rx
    real(dp), intent(out):: x(npng2), w(npng2),&
            an(npn1), ann(npn1, npn1),s(npng2), ss(npng2)
    real(dp) x1(npng2), w1(npng2), &
            x2(npng2), w2(npng2), dd(npn1)
    common/revval/ rx

    !data pi/3.141592653589793d0/

    do n = 1, nmax
        nn = n*(n + 1)
        an(n) = dble(nn)
        d = dsqrt(dble(2*n + 1)/dble(nn))
        dd(n) = d
        do n1 = 1, n
            ddd = d*dd(n1)*0.5d0
            ann(n, n1) = ddd
            ann(n1, n) = ddd
        end do
    end do

    ng = 2*ngauss

    ! gif division points and weights

    !number of gauss integration
    !intervals from eps
    neps = max(eps, 1.d0/eps)

    select case (np)
    case(-1) ! spheroid

        if(neps==1) then

            call gauss(ng, 0, 0, x, w)

        else if (neps==2) then

            call gauleg(-1.d0, 0.d0, x, w, ngauss)

            do i = 1, ngauss

                w(i + ngauss) = w(ngauss - i + 1)
                x(i + ngauss) = -x(ngauss - i + 1)

            enddo

        else if (neps>2) then

            neps = 3*neps
            ng1 = dble(ngauss)/neps

            ee = eps*eps
            ee1 = ee - 1d0

            xav = 0.d0

            do i = 1, neps

                xi1 = dble(i)/(neps + 1)

                cc = xi1*xi1
                si = 1d0 - cc
                x2(i) = abs(xi1*si*ee1/(si + ee*cc))       !|dr(theta)/dtheta|
                xav = xav + 1.d0/x2(i)

            enddo

            xav = xav           !averaged 1/|dr(theta)/dtheta|

            !_____ estimate integration intervals:

            do i = 1, neps

                x2(i) = 1.d0/(xav*x2(i))

            enddo

            do i = 1, neps

                if(i==1) then

                    xi1 = 0.d0
                    xi2 = x2(1)

                else

                    xi2 = xi2 + x2(i)

                end if

                jg = ngauss + (i - 1)*ng1

                if(i==neps) ng1 = ngauss - (i - 1)*ng1

                call gauleg(xi1, xi2, x1, w1, ng1)

                xi1 = xi2

                do  j = 1, ng1
                    w(jg + j) = w1(j)
                    x(jg + j) = x1(j)
                enddo !j

            enddo !i

            ! assuming mirror symmetry in the $\theta=\pi/2$ plane

            do  i = 1, ngauss
                w(i) = w(ng - i + 1)
                x(i) = -x(ng - i + 1)
            enddo

        endif           !neps

    case(-2, -9) ! cylinder or nanorod
        !*****************   only involves cylinders  **********************

        ng1 = dble(ngauss)/2d0
        ng2 = ngauss - ng1
        !-cos of separation angle between
        !horizontal and vertical cylinder
        !faces
        xx = -dcos(datan(eps))

        ! gif division points and weights

        call gauss(ng1, 0, 0, x1, w1)         !for (0,ng1)
        call gauss(ng2, 0, 0, x2, w2)         !for (ng1+1,ngauss=ng1+ng2)

        ! in gauss (n,ind1,ind2,z,w):
        ! ind1 = 0 - interval (-1,1),
        ! ind1 = 1 - (0,1)
        ! ind2 = 1 results are printed.

        do i = 1, ng1
            w(i) = 0.5d0*(xx + 1d0)*w1(i)
            x(i) = 0.5d0*(xx + 1d0)*x1(i) + 0.5d0*(xx - 1d0)
        end do
        do i = 1, ng2
            w(i + ng1) = -0.5d0*xx*w2(i)
            x(i + ng1) = -0.5d0*xx*x2(i) + 0.5d0*xx
        end do

        ! assuming mirror symmetry in the $\theta=\pi/2$ plane

        do i = 1, ngauss
            w(ng - i + 1) = w(i)
            x(ng - i + 1) = -x(i)
        end do
        !*****************************************************************
    case(-4) ! cut sphere on top

        xtheta = dacos((eps - rx)/rx)
        xx = dsin(xtheta)/(pi - xtheta)
        ng2 = xx*dble(ng)
        ng1 = ng - ng2
        theta0 = 1.d0/sqrt(8.d0*rx/eps - 3.d0)  !cosine of the separation angle
        xx = theta0

        call gauleg(-1.d0, theta0, x1, w1, ng1)       !for (0,ng1)
        call gauleg(theta0, 1.d0, x2, w2, ng2)        !for (ng2+1,ng=ng1+ng2)

        do  i = 1, ng1
            w(i) = w1(i)
            x(i) = x1(i)
        enddo

        do i = 1, ng2

            w(i + ng1) = w2(i)
            x(i + ng1) = x2(i)

        enddo

        !*****************************************************************
    case(-5) ! cut sphere on its bottom

        xtheta = dacos((eps - rx)/rx)
        xx = dsin(xtheta)/(pi - xtheta)
        ng1 = xx*dble(ng)
        ng2 = ng - ng1
        theta0 = -1.d0/sqrt(8.d0*rx/eps - 3.d0)  !cosine of the separation angle

        call gauleg(-1.d0, theta0, x1, w1, ng1)       !for (0,ng1)
        call gauleg(theta0, 1.d0, x2, w2, ng2)        !for (ng2+1,ng=ng1+ng2)

        do  i = 1, ng1
            w(i) = w1(i)
            x(i) = x1(i)
        enddo

        do i = 1, ng2

            w(i + ng1) = w2(i)
            x(i + ng1) = x2(i)

        enddo
        !*****************************************************************

    case default

        call gauss(ng, 0, 0, x, w)
        !call gauleg(-1.d0,1.d0,x,w,ng)

    end select

    if (np>-4) then !mirror symmetry present

        do i = 1, ngauss
            y = x(i)
            y = 1d0/(1d0 - y*y)
            ss(i) = y
            ss(ng - i + 1) = y
            y = dsqrt(y)
            s(i) = y
            s(ng - i + 1) = y
        end do

    else !mirror symmetry absent

        do i = 1, ng
            y = x(i)
            y = 1d0/(1d0 - y*y)
            ss(i) = y
            y = dsqrt(y)
            s(i) = y
        end do

    end if

    return
end

!**********************************************************************

subroutine vary(lam, mrr, mri, a, eps, &
        rsnm, ht, np, ngauss, x, &
        p, ppi, pir, pii, r, dr, ddr, drr, dri, nmax)
    !--------/---------/---------/---------/---------/---------/---------/--
    ! >>> lam,mrr,mri,a,eps,np,ngauss,x,p,nmax
    ! <<< ppi,pir,pii,r,dr,ddr,drr,dri
    !=========================
    !  lam - wavelength of incident light in the ambient
    !  mrr - the real part of the refractive index
    !  mri - the imaginary  part of the refractive index
    !  a=rat*axi, where rat and axi are the main program input parameters
    !      rat = 1 - particle size is specified in terms of the
    !                equal-volume-sphere radius
    !      rat.ne.1 - particle size is specified in terms of the
    !                equal-surface-area-sphere radius
    !  axi - equivalent-(volume/surface-area)-sphere radius
    !  np - particle shape class
    !  eps - shape deformation parameter within a given particle shape class
    !  ngauss - the number of gauss integration division points
    !           in the integral over theta
    !  nmax - angular momentum cutoff
    !  p=dacos(-1d0)
    !  wv=p*2d0/lam - wave vector
    !  ppi=wv*wv
    !  pir=ppi*mrr
    !  pii=ppi*mri
    !  r=r^2(\theta)            for axially symmetric particles
    !  dr=dr(\theta)/(d\theta)  for axially symmetric particles
    !  ddr=\lambda/[2*\pi*r(\theta)]=1/(k_out*r(\theta))
    !  drr=(mrr/(mrr**2+mri**2))*(\lambda/[2*\pi*r(\theta)])
    !  dri=-(mri/(mrr**2+mri**2))*(\lambda/[2*\pi*r(\theta)])
!--------/---------/---------/---------/---------/---------/---------/--
    use libcylinder
    implicit none
    integer np, ng, ngauss, nmax, nnmax1, nnmax2, i
    real(dp), intent(in)::  a, eps, mrr, mri, lam
    real(dp), intent(out):: ppi, pir, pii
    real(dp) p, pri, prr, ta, tb,rsnm, ht, wv, v, v1, v2, vv
    real(dp), intent(in)::  x(npng2)
    real(dp), intent(out):: r(npng2), dr(npng2), ddr(npng2), &
            drr(npng2), dri(npng2)
    real(dp) z(npng2), zr(npng2), zi(npng2)

    ng = ngauss*2
    ht = 0d0

    ! decision tree to specify particle shape:
    select case (np)
    case(0:) ! chebyshev particle
        call rsp_chebyshev(x, np, r, dr)
    case (-1) ! oblate/prolate spheroids
        call rsp_spheroid(x, r, dr)
    case (-2) ! oblate/prolate cylinder
        call rsp_cylinder(x, r, dr)
    case (-3) ! distorted chebyshev droplet
        call rsp_droplet(x, r, dr)
    case (-4) ! sphere cut by a plane on its top
        call rsp_sphere_cut_top(x, ng, rsnm, eps, r, dr)
    case (-5) ! sphere cut by a plane on its bottom
        call rsp_sphere_cut_bottom(x, ng, rsnm, eps, r, dr)
    case (-6) ! upwardly oriented cone
        call rsp_cone_up(x, ng, rsnm, ht, r, dr)
    !case (-7) ! cone cut on its top
    !    call rsp_cone_cut_top(x, ng, rsnm, ht, r, dr)
    !case (-8) ! cone on a cylinder
    !    call rsp_cone_on_cylinder(x, ng, rsnm, ht, r, dr)
    case (-9) ! nanorod
        call rsp_nanorod (x, r, dr)
    end select

    wv = p*2d0/lam                 !wave vector k_out
    ppi = wv*wv
    pir = ppi*mrr
    pii = ppi*mri
    v = 1d0/(mrr*mrr + mri*mri)
    prr = mrr*v
    pri = -mri*v
    ta = 0d0

    do i = 1, ng
        vv = dsqrt(r(i))           !r(\theta)
        v = vv*wv                  !k_out*r(\theta)
        ta = max(ta, v)            !max. size parameter
        vv = 1d0/v
        ddr(i) = vv                !1/(k_out*r(\theta))
        drr(i) = prr*vv
        dri(i) = pri*vv

!Assingning arguments of the spherical Bessel functions:
        v1 = v*mrr
        v2 = v*mri
        z(i) = v              !k_out*r(\theta)=(2\pi/\lambda)*r
        zr(i) = v1            !=(2\pi/\lambda)*r*mrr
        zi(i) = v2            !=(2\pi/\lambda)*r*mri
    end do

    if (nmax>npn1) then
        print 90, nmax, npn1
        stop
    end if
    90 format(' nmax = ', i2, ', i.e., greater than ', i3)

    ! ta is the ``max. size parameter", max(2*pi*sqrt(ri)/lambda)

    tb = ta*dsqrt(mrr*mrr + mri*mri)     !=ta*epsin
    tb = dmax1(tb, dble(nmax))
    !
    nnmax1 = 1.2d0*dsqrt(dmax1(ta, dble(nmax))) + 3d0
    nnmax2 = (tb + 4d0*(tb**0.33333d0) + 1.2d0*dsqrt(tb))  !wiscombe bound
    nnmax2 = nnmax2 - nmax + 5

    ! generate arrays of bessel functions at ngauss gif division
    ! points and store them in the common block /cbess/
    !
    cbess%wv = WV
    cbess%mrr = mrr
    cbess%mri = mri
    call bess(z, zr, zi, ng, nmax, nnmax1)

    return
end

!**********************************************************************

subroutine rsp_sphere_cut_top(x, ng, rev, eps, r, dr)
    !--------/---------/---------/---------/---------/---------/---------/--
    ! >>> x,ng,rev,eps
    ! <<< r,dr
    !=========================
    !   activated for np=-4
    !
    !   similar to rsp_sphere_cut_bottom, except for that the plane cut is on the sphere "top"
    !   ===> cosine of the separation angle has the same magnitude as in rsp_sphere_cut_bottom,
    !        but is always positive in the present case !!!
    !
    !   calculation of the functions
    !              r(i)=r(y)**2 and dr(i)=((d/dy)r(y))/r(y)
    !   for a sphere cut by the plane specified by the parameters
    !   rev and  eps at ngauss gauss integration formula (gif) division points
    !   in the integral over theta. here y=acos(x)=theta and eps=2*r_0/h,
    !   where r_0 is the radius of the original uncut sphere, whereas h
    !   is the height (along the axial symmetry axis) of the resulting
    !   cut sphere.
    !
    !   the origin of coordinates is located along the axis of symmetry
    !                  midway the plane and sphere bottom.
    !
    !           ===> note that always eps > 1
    !   ===
    !   x - gif division points \cos\theta_j -  y = arccos x
    !   rev ... the radius of the original uncut sphere
    !   eps ...  h is the height (along the axial symmetry axis)
    !            of the resulting cut sphere. note that always eps.lt.2*rev
    !   theta0 ... a cosine of the separation angle between two different
    !              functional dependences of r(\theta), that along the sphere
    !              surface and that along the plane surface
    !   ng=2*ngauss ... the number of gif division points
    !
    !   1 <= i <= ngauss
    !--------/---------/---------/---------/---------/---------/---------/--
    use libcylinder
    implicit none
    integer, intent(in) :: ng
    real(dp), intent(in) :: rev, eps, x(ng)
    real(dp), intent(out) :: r(ng), dr(ng)

    integer i
    real(dp) a, rad, co, cc, rthet, si, ss, theta0

    if (eps>=2.d0*rev) then
        write(6, *)'invalid parameters for a cut sphere!'
        write(6, *)'execution stopped!'
        stop
    end if

    theta0 = 1.d0/sqrt(8.d0*rev/eps - 3.d0)  !cosine of the separation angle

    do i = 1, ng

        co = x(i)
        cc = co*co
        ss = 1.d0 - cc
        si = dsqrt(ss)                    !=\sin\theta

        if (co<theta0) then        ! r(\theta) along the sphere surface

            a = rev - eps/2.d0

            rad = a*co + sqrt(rev**2 - (a*si)**2)
            rthet = -a*si - co*si*a**2/sqrt(rev**2 - (a*si)**2)
            !c          rad=1.d-10
            !c          rthet=0.d0

        else if (co>=theta0) then   ! r(\theta) along the plane surface
            !                                         !    (co positive)
            rad = eps/(2.d0*co)
            rthet = eps*si/(2.d0*co**2)

        end if

        dr(i) = rthet/rad
        r(i) = rad*rad

    end do

    return
end

!**********************************************************************

subroutine rsp_sphere_cut_bottom(x, ng, rev, eps, r, dr)
    !--------/---------/---------/---------/---------/---------/---------/--
    ! >>> x,ng,rev,eps
    ! <<< r,dr
    !=========================
    !   activated for np=-5
    !
    !   similar to rsp_sphere_cut_top, except for that the plane cut is on the sphere "bottom"
    !   ===> cosine of the separation angle has the same magnitude as in rsp_sphere_cut_top,
    !        but is always negative in the present case !!!
    !
    !   calculation of the functions
    !              r(i)=r(y)**2 and dr(i)=((d/dy)r(y))/r(y)
    !   for a sphere cut by the plane specified by the parameters
    !   rev and  eps at ngauss gauss integration formula (gif) division points
    !   in the integral over theta. here y=acos(x)=theta and eps=2*r_0/h,
    !   where r_0 is the radius of the original uncut sphere, whereas h
    !   is the height (along the axial symmetry axis) of the resulting
    !   cut sphere.
    !
    !   the origin of coordinates is located along the axis of symmetry
    !                  midway the plane and sphere top.
    !
    !                  ===>  note that always eps > 1
    !   ===
    !   x - gif division points \cos\theta_j -  y = arccos x
    !   rev ... the radius of the original uncut sphere
    !   eps ...  h is the height (along the axial symmetry axis)
    !            of the resulting cut sphere. note that always eps.lt.2*rev
    !   theta0 ... a cosine of the separation angle between two different
    !              functional dependences of r(\theta), that along the sphere
    !              surface and that along the plane surface
    !   ng=2*ngauss ... the number of gif division points
    !
    !   1 <= i <= ngauss
    !--------/---------/---------/---------/---------/---------/---------/--
    use libcylinder
    implicit none
    integer, intent(in) :: ng
    real(dp), intent(in) :: rev, eps, x(ng)
    real(dp), intent(out) :: r(ng), dr(ng)

    integer i
    real(dp) rad, rthet, a, co, cc, si, ss, theta0

    if (eps>=2.d0*rev) then
        write(6, *)'invalid parameters for a cut sphere!'
        write(6, *)'execution stopped!'
        stop
    end if

    theta0 = -1.d0/sqrt(8.d0*rev/eps - 3.d0)  !cosine of the separation angle resent
    !                                           !case
    do i = 1, ng

        co = x(i)
        cc = co*co
        ss = 1.d0 - cc
        si = dsqrt(ss)                  !=\sin\theta

        if (co>theta0) then        !r(\theta) along the sphere surface

            a = rev - eps/2.d0

            rad = -a*co + sqrt(rev**2 - (a*si)**2)
            rthet = a*si - co*si*a**2/sqrt(rev**2 - (a*si)**2)
            !c          rad=1.d-10
            !c          rthet=0.d0

        else if (co<=theta0) then   ! r(\theta) along the plane surface
            !                                         !        (co negative)
            rad = -eps/(2.d0*co)
            rthet = -eps*si/(2.d0*co**2)

        end if

        dr(i) = rthet/rad
        r(i) = rad*rad

    end do

    return
end

!**********************************************************************

subroutine rsp_cone_up(x, ng, rev, ht, r, dr)
    !--------/---------/---------/---------/---------/---------/---------/--
    ! >>> x,ng,rev,ht
    ! <<< r,dr
    !=========================
    !   activated for np=-6
    !
    !   calculation of the functions
    !              r(i)=r(y)**2 and dr(i)=((d/dy)r(y))/r(y)
    !   for an upwardly pointing singular cone specified by the parameters
    !   rev and eps at ngauss gauss integration formula (gif) division points
    !   in the integral over theta.
    !   rev ... the the half width of the cone base
    !   ht ...  cone height (along the axial symmetry axis)
    !
    !   the origin of coordinates is placed on the axis of axial symmetry,
    !         at the distance (h/3) from the base to the cone
    !                                    <===>
    !     the base and slant of the cone form a triangle, and the origin
    !          of coordinates is placed at the triangle centroid.
    !
    !   ===>  the length of the cone slant a = dsqrt(h**2+rev**2)
    !   ===>  the length of the median of the cone slant ma = dsqrt(ma**2+rev**2/2.d0)/2.d0
    !
    !   the cone base end points and the traingle centroid form a second triangle.
    !   if we denote by 2*theta1 the second triangle angle at the centroid,
    !   the separation angle theta_s as viewed from the origin (between the slant and
    !   base part of the cone surfaces) is given by
    !
    !               - cos (theta_s) =  cos (theta1) = (h/3)/(2*ma/3)  = h/(2*ma)
    !
    !   then, for  theta>theta_s (i.e., cos (theta) < cos (theta_s) ):
    !
    !                      r(theta)=-h/(3*cos(theta))
    !
    !   for  theta<theta_s (i.e., cos (theta) > cos (theta_s) ), one first considers
    !   a triangle formed by the cone apex, the origin of coordinates and r(theta).
    !   let theta_v denote the angle of this triangle at the cone apex.
    !   according to the law of cosines:
    !
    !         c^2 = r^2+4h^2/9 - [4rh\cos(\theta)]/3
    !         r^2 = c^2+4h^2/9 - [4rh\cos(\theta_v)]/3
    !
    !   where c is the length of the traingle opposite to the origin of coordinates.
    !   upon combining the last two equations one arrives
    !
    !                 r\cos(\theta)+c\cos(\theta_v) = 2h/3    (*)
    !
    !   according to the law of sines:
    !
    !                      c/sin(\theta) = r(\theta)/sin(\theta_v)
    !
    !   when substituting back to (*), one finds
    !
    !               r = (2h/3)/[\cos(\theta)+\cot(\theta_v) \sin(\theta)]
    !
    !   \cot(\theta_v) can be easily determined from the very first triangle,
    !
    !              \cot(\theta_v) = h/(b/2) = 2h/b
    !
    !   conus volume = (pi*rev**2*h)/3.d0
    !   here y=acos(x)=theta and eps=2*r_0/h,
    !   where r_0 is the radius of the original uncut sphere, whereas h
    !   is the height (along the axial symmetry axis) of the resulting
    !   cut sphere.
    !   ===>
    !
    !   x - gif division points \cos\theta_j -  y = arccos x
    !
    !   ng=2*ngauss ... the number of gif division points
    !
    !   1 <= i <= ngauss
    !--------/---------/---------/---------/---------/---------/---------/--
    use libcylinder
    implicit none
    integer, intent(in) :: ng
    real(dp), intent(in) :: rev, ht, x(ng)
    real(dp), intent(out) :: r(ng), dr(ng)

    integer :: i
    real(dp) :: ma, co, si, cc, ss, rad, theta0, rthet

    if (ht<1e-8) stop 'ht should be more than zero'
    ma = dsqrt(ht**2 + rev**2)              !=the length of the cone slant
    ma = dsqrt(ma**2 + 8.d0*rev**2)/2.d0    !=the length of the median of the slant

    theta0 = -ht/(2.d0*ma)                !=cos of the separation angle
    !                                         !  (always negative)

    do i = 1, ng

        co = x(i)
        cc = co*co
        ss = 1.d0 - cc
        si = dsqrt(ss)                    !=\sin\theta

        if (co>theta0) then          ! theta < theta0,
            !                                           ! i.e. r(\theta) along the cone slant surface

            ma = co + ht*si/rev
            rad = 2.d0*ht/(3.d0*ma)
            rthet = 2.d0*ht*(si - ht*co/rev)/(3.d0*ma**2)

        else if (co<=theta0) then     ! theta > theta0, i.e. r(\theta) along the cone base surface

            rad = -ht/(3.d0*co)       ! always positive (co<0 on the base)
            rthet = rad*si/co         ! always negative

        end if

        r(i) = rad*rad
        dr(i) = rthet/rad

    end do

    return
end


!*********************************************************************
subroutine tmatr0(ngauss, x, w, an, ann, ppi, pir, pii, r, dr, ddr, &
        drr, dri, nmax, ncheck, naxsm)
    !--------/---------/---------/---------/---------/---------/---------/--
    ! >>> ngauss,x,w,an,ann,ppi,pir,pii,r,dr,ddr,drr,dri,nmax,ncheck
    ! <<< common blocks /tmat99/, /ct/ (for main),  and /ctt/ (for tt)
    !=====================
    !
    !  determines the t-matrix of an axially symmetric scatterer
    !                           for m=0
    !
    !  ngauss - the number of gif division points
    !  x=\cos\theta  - gif division points
    !  w - gif weights
    !  an(n)=n*(n+1)
    !  ann(l_1,l_2)=\sqrt{\fr{(2 l_1+1)}{l_1(l_1+1)} }
    !                       \sqrt{\fr{(2 l_2+1)}{l_2(l_2+1)} }/2
    !  nmax - angular momentum cutoff
    !  ncheck  -  .eq.0  then  ngss=2*ngauss, factor=1d0
    !             .eq.1  then  ngss = ngauss, factor=2d0
    !  naxsm   -  .eq.0 : gauss abscissas do not have +/- theta symmetry
    !             .eq.1 : gauss abscissas have +/- theta symmetry
    !  p=dacos(-1d0)
    !  pi=p*2d0/lam - wave vector
    !  ppi=pi*pi
    !  pir=ppi*mrr
    !  pii=ppi*mri
    !  r=r^2(\theta)                        for axially symmetric particles
    !  dr=[dr(\theta)/(d\theta)]/r(\theta)  for axially symmetric particles
    !  ddr=\lambda/[2*\pi*r(\theta)]=1/(k_out*r)
    !  drr=(mrr/(mrr**2+mri**2))*(\lambda/[2*\pi*r(\theta)])
    !                  = re 1/(k_in*r)
    !  dri=-(mri/(mrr**2+mri**2))*(\lambda/[2*\pi*r(\theta)])
    !                  = im 1/(k_in*r)
    !  refractive index outside is assumed to real, whereas inside
    !  a scatterer, refractive index is allowed to be complex in general.
    !  consequently, the bessel function j_l(k_in*r) will in general
    !  be complex. the routine below performs waterman surface integral
    !  separately for the real and imaginary parts of the integrand.
    !
    !--------/---------/---------/---------/---------/---------/---------/--
    use libcylinder
    implicit none
    integer ngauss, nmax, ncheck, naxsm, i, i1, i2, k1, k2, kk1, kk2, &
            mm1, n, n1, n2, ng, ngss, nm, nnmax
    real(dp) :: ppi, pir, pii, d1n1, a12, a21, a22, ai12, ai21, an1, &
            an12, an2, ar12, ar21, b1i, b1r, b2i, b2r, b3i, b3r, b4i, b4r, &
            b5r, b5i, c1r, c1i, c2r, c2i, c3r, c3i, c4r, c4i, c5r, c5i, &
            d1n2, d2n1, d2n2, dd1, dd2, ddri, drii, drri, f1, f2, &
            factor, gi12, gi21, gr12, gr21, &
            qdj1, qdji2, qdjr2, qdy1, qj1, &
            qji2, qjr2, qy1, rri, si, tai12, tai21, &
            tar12, tar21, tgi12, tgi21, &
            tgr12, tgr21, tppi, tpii, tpir, uri
    real(dp)  x(npng2), w(npng2), an(npn1), &
            r(npng2), dr(npng2), sig(npn2), &
            ddr(npng2), drr(npng2), &
            d1(npng2, npn1), d2(npng2, npn1), &
            dri(npng2), rr(npng2), &
            dv1(npn1), dv2(npn1)
    real(dp)  r11(npn1, npn1), r12(npn1, npn1), &
            r21(npn1, npn1), r22(npn1, npn1), &
            i11(npn1, npn1), i12(npn1, npn1), &
            i21(npn1, npn1), i22(npn1, npn1), &
            rg11(npn1, npn1), rg12(npn1, npn1), &
            rg21(npn1, npn1), rg22(npn1, npn1), &
            ig11(npn1, npn1), ig12(npn1, npn1), &
            ig21(npn1, npn1), ig22(npn1, npn1), &
            ann(npn1, npn1), &
            qr(npn2, npn2), qi(npn2, npn2), &
            rgqr(npn2, npn2), rgqi(npn2, npn2), &
            tqr(npn2, npn2), tqi(npn2, npn2), &
            trgqr(npn2, npn2), trgqi(npn2, npn2)
    !c      real(dp) tr1(npn2,npn2),ti1(npn2,npn2)
    !
    common /tmat99/&
            r11, r12, r21, r22, i11, i12, i21, i22, rg11, rg12, rg21, rg22, &
            ig11, ig12, ig21, ig22           !only between tmatr routines
    !c      common /ct/ tr1,ti1                       !output from tt routine
    common /ctt/ qr, qi, rgqr, rgqi              !input for tt routine
    !
    mm1 = 1
    nnmax = nmax + nmax
    ng = 2*ngauss
    ngss = ng
    factor = 1d0
    !
    if (ncheck==1) then         !theta=pi/2 is scatterer mirror symmetry plane
        ngss = ngauss
        factor = 2d0
    else                       !theta=pi/2 is not a scatterer mirror symmetry plane
    endif
    !
    si = 1d0
    do n = 1, nnmax
        si = -si
        sig(n) = si              !=(-1)**n
    end do
    !
    ! assigning Wigner d-matrices - assuming mirror symmetry
    ! in the \theta=\pi/2 plane:
    do i = 1, ngauss
        i1 = ngauss - i + 1
        !c         i2=ngauss+i
        if (ncheck==0) i2 = ngauss + i
        !
        call vig (x(i1), nmax, 0, dv1, dv2)
        !
        do n = 1, nmax
            dd1 = dv1(n)
            dd2 = dv2(n)
            d1(i1, n) = dd1
            d2(i1, n) = dd2
            ! if theta=pi/2 is not a scatterer mirror symmetry plane but
            ! gauss abscissas are still chosen symmetrically:
            if ((ncheck==0).and.(naxsm==1)) then
                ! using (4.2.4) and (4.2.6) of {ed},
                !           d_{0m}^{(l)}(\pi-\theta) = (-1)^{l+m} d_{0m}^{(l)}(\theta)
                ! (4.2.5) of {ed}:                   = (-1)^{l} d_{0 -m}^{(l)}(\theta)
                si = sig(n)                  !=(-1)**n
                d1(i2, n) = dd1*si
                d2(i2, n) = -dd2*si
            end if
        enddo
        ! if neither scatterer nor gauss abscissas have theta=pi/2
        ! as a mirror symmetry plane:
        if ((ncheck==0).and.(naxsm==0)) then
            !
            call vig (x(i2), nmax, 0, dv1, dv2)
            !
            do n = 1, nmax
                dd1 = dv1(n)
                dd2 = dv2(n)
                d1(i2, n) = dd1
                d2(i2, n) = dd2
            enddo
        end if
    end do
    !
    !  assigning r^2(\theta)*weight product:
    do i = 1, ngss
        rr(i) = w(i)*r(i)
        !c           if (dr(i).eq.0.d0) rr(i)=0.d0   !temporarily only
    end do
    !
    do n1 = mm1, nmax
        an1 = an(n1)
        do n2 = mm1, nmax
            an2 = an(n2)
            ar12 = 0d0
            ar21 = 0d0
            ai12 = 0d0
            ai21 = 0d0
            gr12 = 0d0
            gr21 = 0d0
            gi12 = 0d0
            gi21 = 0d0
            !        open(nout+3,file='surfint.dat')   !gauss convergence check
            if (ncheck/=1.or.sig(n1 + n2)>=0d0) then
                !
                ! gauss integration loop:
                !
                do i = 1, ngss    !=ngauss   if ncheck.eq.1
                    !                                  !=2*ngauss if ncheck.eq.0
                    d1n1 = d1(i, n1)
                    d2n1 = d2(i, n1)
                    d1n2 = d1(i, n2)
                    d2n2 = d2(i, n2)
                    a12 = d1n1*d2n2
                    a21 = d2n1*d1n2
                    a22 = d2n1*d2n2
                    !                    aa1=a12+a21
                    ! vector spherical harmonics:
                    !  since refractive index is allowed to be complex in general,
                    !  the bessel function j_l(k_in*r) is complex. the code below
                    !  performs a separation of the complex integrand in waterman's
                    !  surface integral into its respective real and imaginary
                    !  parts:
                    ! bessel functions of the exterior argument:
                    qj1 = cbess%j(i, n1)
                    qy1 = cbess%y(i, n1)
                    qdj1 = cbess%dj(i, n1)
                    qdy1 = cbess%dy(i, n1)
                    ! bessel functions of the interior argument:
                    qjr2 = cbess%jr(i, n2)
                    qji2 = cbess%ji(i, n2)
                    qdjr2 = cbess%djr(i, n2)
                    qdji2 = cbess%dji(i, n2)
                    !_____________________
                    ! re and im of j_{n2}(k_{in}r) j_{n1}(k_{out}r):
                    c1r = qjr2*qj1
                    c1i = qji2*qj1
                    ! re and im of j_{n2}(k_{in}r) h_{n1}(k_{out}r):
                    b1r = c1r - qji2*qy1
                    b1i = c1i + qjr2*qy1
                    ! re and im of j_{n2}(k_{in}r) [k_{out}r j_{n1}(k_{out}r)]'/(k_{out}r):
                    c2r = qjr2*qdj1
                    c2i = qji2*qdj1
                    ! re and im of j_{n2}(k_{in}r) [k_{out}r h_{n1}(k_{out}r)]'/(k_{out}r):
                    b2r = c2r - qji2*qdy1
                    b2i = c2i + qjr2*qdy1
                    ddri = ddr(i)               !1/(k_{out}r)
                    ! re and im of [1/(k_{out}r)]*j_{n2}(k_{in}r) j_{n1}(k_{out}r)
                    c3r = ddri*c1r
                    c3i = ddri*c1i
                    ! re and im of [1/(k_{out}r)]*j_{n2}(k_{in}r) h_{n1}(k_{out}r):
                    b3r = ddri*b1r
                    b3i = ddri*b1i
                    ! re and im of  [k_{in}r j_{n2}(k_{in}r)]'/(k_{in}r)
                    !                         *j_{n1}(k_{out}r):
                    c4r = qdjr2*qj1
                    c4i = qdji2*qj1
                    ! re and im of [k_{in}r j_{n2}(k_{in}r)]'/(k_{in}r)
                    !                         * h_{n1}(k_{out}r):
                    b4r = c4r - qdji2*qy1
                    b4i = c4i + qdjr2*qy1
                    drri = drr(i)               !re[1/(k_{in}r)]
                    drii = dri(i)               !im[1/(k_{in}r)]
                    ! re and im of [1/(k_{in}r)] j_{n2}(k_{in}r) j_{n1}(k_{out}r):
                    c5r = c1r*drri - c1i*drii
                    c5i = c1i*drri + c1r*drii
                    ! re and im of [1/(k_{in}r)] j_{n2}(k_{in}r) h_{n1}(k_{out}r):
                    b5r = b1r*drri - b1i*drii
                    b5i = b1i*drri + b1r*drii
                    !%%%%%%%  forming integrands of j-matrices (j^{11}=j^{22}=0 for m=0): %%%%%%%%
                    uri = dr(i)        !dr/(d\theta)
                    rri = rr(i)        !w(i)*r^2(\theta)
                    ! w(i)*r^2(\theta)*d2n1*d2n2:
                    f1 = rri*a22      !prefactor containing r^2(\theta)<->hat{r} part
                    ! n1*(n1+1)*w(i)*r(\theta)*[dr/(d\theta)]*d1n1*d2n2:
                    f2 = rri*uri*an1*a12     !prefactor containing r(\theta)*[dr/(d\theta)]
                    !                                          !hat{theta} part
                    ar12 = ar12 + f1*b2r + f2*b3r        !~re j^{12}
                    ai12 = ai12 + f1*b2i + f2*b3i        !~im j^{12}
                    gr12 = gr12 + f1*c2r + f2*c3r        !~re rg j^{12}
                    gi12 = gi12 + f1*c2i + f2*c3i        !~im rg j^{12}
                    !*  n2*(n2+1)*w(i)*r(\theta)*[dr/(d\theta)]*d2n1*d1n2:
                    f2 = rri*uri*an2*a21     !prefactor containing r(\theta)*[dr/(d\theta)]
                    !                                          !hat{theta} part
                    ar21 = ar21 + f1*b4r + f2*b5r        !~re j^{21}
                    ai21 = ai21 + f1*b4i + f2*b5i        !~im j^{21}
                    gr21 = gr21 + f1*c4r + f2*c5r        !~re rg j^{21}
                    gi21 = gi21 + f1*c4i + f2*c5i        !~im rg j^{21}
                end do               !end of gauss integration
                !                write(nout+3,*)'n1=',n1,'   n2=',n2
                !                write(nout+3,*)'ar12=', ar12
                !                write(nout+3,*)'ai12=', ai12
                !                write(nout+3,*)'ar21=', ar21
                !                write(nout+3,*)'ai21=', ai21
                !                write(nout+3,*)'gr12=', gr12
                !                write(nout+3,*)'gi12=', gi12
                !                write(nout+3,*)'gr21=', gr21
                !                write(nout+3,*)'gi21=', gi21
                !%%%%%%%%%%%%%  forming j-matrices (j^{11}=j^{22}=0 for m=0):
            end if
            an12 = ann(n1, n2)*factor
            r12(n1, n2) = ar12*an12
            r21(n1, n2) = ar21*an12
            i12(n1, n2) = ai12*an12
            i21(n1, n2) = ai21*an12
            rg12(n1, n2) = gr12*an12
            rg21(n1, n2) = gr21*an12
            ig12(n1, n2) = gi12*an12
            ig21(n1, n2) = gi21*an12
        end do
    end do            !end of the loop over angular momenta
    !      close(nout+3)
    !%%%%%%%%%%%%%%%%%%%%%%%  forming q and rgq -matrices
    tpir = pir                 !re [1/k_{in}^2]
    tpii = pii                 !im [1/k_{in}^2]
    tppi = ppi                 !1/k_{out}^2
    nm = nmax
    do n1 = mm1, nmax
        k1 = n1 - mm1 + 1
        kk1 = k1 + nm
        do n2 = mm1, nmax
            k2 = n2 - mm1 + 1
            kk2 = k2 + nm
            tar12 = i12(n1, n2)
            tai12 = -r12(n1, n2)
            tgr12 = ig12(n1, n2)
            tgi12 = -rg12(n1, n2)
            tar21 = -i21(n1, n2)
            tai21 = r21(n1, n2)
            tgr21 = -ig21(n1, n2)
            tgi21 = rg21(n1, n2)
            tqr(k1, k2) = tpir*tar21 - tpii*tai21 + tppi*tar12
            tqi(k1, k2) = tpir*tai21 + tpii*tar21 + tppi*tai12
            trgqr(k1, k2) = tpir*tgr21 - tpii*tgi21 + tppi*tgr12
            trgqi(k1, k2) = tpir*tgi21 + tpii*tgr21 + tppi*tgi12
            tqr(k1, kk2) = 0d0
            tqi(k1, kk2) = 0d0
            trgqr(k1, kk2) = 0d0
            trgqi(k1, kk2) = 0d0
            tqr(kk1, k2) = 0d0
            tqi(kk1, k2) = 0d0
            trgqr(kk1, k2) = 0d0
            trgqi(kk1, k2) = 0d0
            tqr(kk1, kk2) = tpir*tar12 - tpii*tai12 + tppi*tar21
            tqi(kk1, kk2) = tpir*tai12 + tpii*tar12 + tppi*tai21
            trgqr(kk1, kk2) = tpir*tgr12 - tpii*tgi12 + tppi*tgr21
            trgqi(kk1, kk2) = tpir*tgi12 + tpii*tgr12 + tppi*tgi21
        end do
    end do
    nnmax = 2*nm
    do n1 = 1, nnmax
        do n2 = 1, nnmax
            qr(n1, n2) = tqr(n1, n2)
            qi(n1, n2) = tqi(n1, n2)
            rgqr(n1, n2) = trgqr(n1, n2)
            rgqi(n1, n2) = trgqi(n1, n2)
        end do
    end do
    !%%%%%%%%%%%%%%%%%%%%%%%  forming resulting t-matrix
    !
    ! calculate the product q^{-1} rg q
    !
    call tt(nmax)
    return
end

!*********************************************************************

subroutine tmatr0_adapt(ngauss, x, w, an, ann, ppi, pir, pii, r, dr, ddr, &
        drr, dri, nmax, ncheck, naxsm)
    !--------/---------/---------/---------/---------/---------/---------/--
    ! >>> ngauss,x,w,an,ann,ppi,pir,pii,r,dr,ddr,drr,dri,nmax,ncheck
    ! <<< common blocks /tmat99/, /ct/ (for main),  and /ctt/ (for tt)
    !=====================
    !
    !  determines the t-matrix of an axially symmetric scatterer
    !                           for m=0
    !
    !  ngauss - the number of gif division points
    !  x=\cos\theta  - gif division points
    !  w - gif weights
    !  an(n)=n*(n+1)
    !  ann(l_1,l_2)=\sqrt{\fr{(2 l_1+1)}{l_1(l_1+1)} }
    !                       \sqrt{\fr{(2 l_2+1)}{l_2(l_2+1)} }/2
    !  nmax - angular momentum cutoff
    !  ncheck  -  .eq.0  then  ngss=2*ngauss, factor=1d0
    !             .eq.1  then  ngss = ngauss, factor=2d0
    !  naxsm   -  .eq.0 : gauss abscissas do not have +/- theta symmetry
    !             .eq.1 : gauss abscissas have +/- theta symmetry
    !  p=dacos(-1d0)
    !  pi=p*2d0/lam - wave vector
    !  ppi=pi*pi
    !  pir=ppi*mrr
    !  pii=ppi*mri
    !  r=r^2(\theta)                        for axially symmetric particles
    !  dr=[dr(\theta)/(d\theta)]/r(\theta)  for axially symmetric particles
    !  ddr=\lambda/[2*\pi*r(\theta)]=1/(k_out*r)
    !  drr=(mrr/(mrr**2+mri**2))*(\lambda/[2*\pi*r(\theta)])
    !                  = re 1/(k_in*r)
    !  dri=-(mri/(mrr**2+mri**2))*(\lambda/[2*\pi*r(\theta)])
    !                  = im 1/(k_in*r)
    !  refractive index outside is assumed to real, whereas inside
    !  a scatterer, refractive index is allowed to be complex in general.
    !  consequently, the bessel function j_l(k_in*r) will in general
    !  be complex. the routine below performs waterman surface integral
    !  separately for the real and imaginary parts of the integrand.
    !
    !--------/---------/---------/---------/---------/---------/---------/--
    use libcylinder
    use lib_adaptive_integration
    implicit none
    integer ngauss, nmax, ncheck, naxsm, i, i1, i2, k1, k2, kk1, kk2, &
            mm1, n, n1, n2, ng, ngss, nm, nnmax
    real(dp) :: ppi, pir, pii, d1n1, a12, a21, a22, ai12, ai21, an1, &
            an12, an2, ar12, ar21, b1i, b1r, b2i, b2r, b3i, b3r, b4i, b4r, &
            b5r, b5i, c1r, c1i, c2r, c2i, c3r, c3i, c4r, c4i, c5r, c5i, &
            d1n2, d2n1, d2n2, dd1, dd2, ddri, drii, drri, f1, f2, &
            factor, gi12, gi21, gr12, gr21, &
            qdj1, qdji2, qdjr2, qdy1, qj1, &
            qji2, qjr2, qy1, si, tai12, tai21, &
            tar12, tar21, tgi12, tgi21, &
            tgr12, tgr21, tppi, tpii, tpir, uri, &
            x_to_rsp(1), r_from_x(1), dr_from_x(1), xi, &
            v
    real(dp)  x(npng2), w(npng2), an(npn1), &
            r(npng2), dr(npng2), sig(npn2), &
            ddr(npng2), drr(npng2), &
            d1(npng2, npn1), d2(npng2, npn1), &
            dri(npng2), rr(npng2), &
            dv1(npn1), dv2(npn1)

    real(dp)  r11(npn1, npn1), r12(npn1, npn1), &
            r21(npn1, npn1), r22(npn1, npn1), &
            i11(npn1, npn1), i12(npn1, npn1), &
            i21(npn1, npn1), i22(npn1, npn1), &
            rg11(npn1, npn1), rg12(npn1, npn1), &
            rg21(npn1, npn1), rg22(npn1, npn1), &
            ig11(npn1, npn1), ig12(npn1, npn1), &
            ig21(npn1, npn1), ig22(npn1, npn1), &
            ann(npn1, npn1), &
            qr(npn2, npn2), qi(npn2, npn2), &
            rgqr(npn2, npn2), rgqi(npn2, npn2), &
            tqr(npn2, npn2), tqi(npn2, npn2), &
            trgqr(npn2, npn2), trgqi(npn2, npn2)
    !real(dp) tr1(npn2,npn2),ti1(npn2,npn2)

    common /tmat99/&
            r11, r12, r21, r22, i11, i12, i21, i22, rg11, rg12, rg21, rg22, &
            ig11, ig12, ig21, ig22           !only between tmatr routines
    !common /ct/ tr1,ti1                       !output from tt routine
    common /ctt/ qr, qi, rgqr, rgqi              !input for tt routine

    mm1 = 1
    nnmax = nmax + nmax
    ng = 2*ngauss
    ngss = ng
    factor = 1d0

    if (ncheck==1) then         !theta=pi/2 is scatterer mirror symmetry plane
        ngss = ngauss
        factor = 2d0
    else                        !theta=pi/2 is not a scatterer mirror symmetry plane

    endif

    si = 1d0
    do n = 1, nnmax
        si = -si
        sig(n) = si              !=(-1)**n
    end do

    ! Assigning Wigner d-matrices - assuming mirror symmetry
    ! in the \theta=\pi/2 plane:

    do i = 1, ngauss

        i1 = ngauss - i + 1
        !i2=ngauss+i
        if (ncheck==0) i2 = ngauss + i

        call vig (x(i1), nmax, 0, dv1, dv2)

        do n = 1, nmax

            dd1 = dv1(n)
            dd2 = dv2(n)
            d1(i1, n) = dd1
            d2(i1, n) = dd2

            ! if theta=pi/2 is not a scatterer mirror symmetry plane but
            ! gauss abscissas are still chosen symmetrically:

            if ((ncheck==0).and.(naxsm==1)) then

                ! using (4.2.4) and (4.2.6) of {ed},
                !           d_{0m}^{(l)}(\pi-\theta) = (-1)^{l+m} d_{0m}^{(l)}(\theta)
                ! (4.2.5) of {ed}:                   = (-1)^{l} d_{0 -m}^{(l)}(\theta)

                si = sig(n)                  !=(-1)**n
                d1(i2, n) = dd1*si
                d2(i2, n) = -dd2*si

            end if

        enddo

        ! if neither scatterer nor gauss abscissas have theta=pi/2
        ! as a mirror symmetry plane:

        if ((ncheck==0).and.(naxsm==0)) then
            !
            call vig (x(i2), nmax, 0, dv1, dv2)
            !
            do n = 1, nmax
                dd1 = dv1(n)
                dd2 = dv2(n)
                d1(i2, n) = dd1
                d2(i2, n) = dd2
            enddo

        end if
    end do

    !  Assigning r^2(\theta)*weight product and Bessel function independent factors:
    !  In Ru et al formulas:
    !  s=n_2/n_1$ is the refractive index constrast and $x(\theta)=k_1 r(\theta)$ is the
    !  {\em polar} angle dependent size parameter of an axially symmetric scatterer.

    do i = 1, ngss
        rr(i) = w(i)*r(i)          !r^2(\cos\theta)+w(i)
        !if (dr(i).eq.0.d0) rr(i)=0.d0   !temporarily only
    end do

    !  Assigning Bessel function dependent factors:
    !  n1,n2\ge |m|

    do n1 = mm1, nmax
        an1 = an(n1)
        do n2 = mm1, nmax
            an2 = an(n2)
            ! Init:
            ar12 = 0d0
            ar21 = 0d0
            ai12 = 0d0
            ai21 = 0d0
            gr12 = 0d0
            gr21 = 0d0
            gi12 = 0d0
            gi21 = 0d0

            !        open(nout+3,file='surfint.dat')   !gauss convergence check

            if (ncheck/=1.or.sig(n1 + n2)>=0d0) then
                !
                ! gauss integration loop for ar12 only:
                ! Gauss integration loop:
                ! Traditional integration transforms \int_0^\pi X* \sin\theta d\theta --> \int_{-1}^1 X* d(\cos\theta)
                ! and generates weights for the independent variable "\cos\theta"
                !
                integrand%n1 = n1
                integrand%n2 = n2
                integrand%nmax = nmax ! used only for bessel evaluation

                ! gauss integration loop (other vars):
                !
                !do i = 1, ngss    !=ngauss   if ncheck.eq.1
                !                  !=2*ngauss if ncheck.eq.0
                !    ar12 = ar12 + w(i)*m0_ar12(x(i))         !~re j^{12}
                !    ai12 = ai12 + w(i)*m0_ai12(x(i))        !~im j^{12}
                !
                !    gr12 = gr12 + w(i)*m0_gr12(x(i))        !~re rg j^{12}
                !    gi12 = gi12 + w(i)*m0_gi12(x(i))        !~im rg j^{12}
                !
                !    ar21 = ar21 + w(i)*m0_ar21(x(i))        !~re j^{21}
                !    ai21 = ai21 + w(i)*m0_ai21(x(i))        !~im j^{21}
                !
                !    gr21 = gr21 + w(i)*m0_gr21(x(i))       !~re rg j^{21}
                !    gi21 = gi21 + w(i)*m0_gi21(x(i))        !~im rg j^{21}
                !end do            !end of gauss integration
                call integrate(m0_ar12, -1.0_dp, 0.0_dp, ar12)
                call integrate(m0_ai12, -1.0_dp, 0.0_dp, ai12)
                call integrate(m0_gr12, -1.0_dp, 0.0_dp, gr12)
                call integrate(m0_gi12, -1.0_dp, 0.0_dp, gi12)
                call integrate(m0_ar21, -1.0_dp, 0.0_dp, ar21)
                call integrate(m0_ai21, -1.0_dp, 0.0_dp, ai21)
                call integrate(m0_gr21, -1.0_dp, 0.0_dp, gr21)
                call integrate(m0_gi21, -1.0_dp, 0.0_dp, gi21)

                !write(nout+3,*)'n1=',n1,'   n2=',n2
                !write(nout+3,*)'ar12=', ar12
                !write(nout+3,*)'ai12=', ai12
                !write(nout+3,*)'ar21=', ar21
                !write(nout+3,*)'ai21=', ai21
                !write(nout+3,*)'gr12=', gr12
                !write(nout+3,*)'gi12=', gi12
                !write(nout+3,*)'gr21=', gr21
                !write(nout+3,*)'gi21=', gi21
                !%%%%%%%%%%%%%  forming j-matrices (j^{11}=j^{22}=0 for m=0):
            end if

            an12 = ann(n1, n2)*factor

            r12(n1, n2) = ar12*an12
            r21(n1, n2) = ar21*an12
            i12(n1, n2) = ai12*an12
            i21(n1, n2) = ai21*an12

            rg12(n1, n2) = gr12*an12
            rg21(n1, n2) = gr21*an12
            ig12(n1, n2) = gi12*an12
            ig21(n1, n2) = gi21*an12
        end do
    end do            !end of the loop over angular momenta

    !      close(nout+3)
    !%%%%%%%%%%%%%%%%%%%%%%%  forming q and rgq -matrices

    tpir = pir                 !re [1/k_{in}^2]
    tpii = pii                 !im [1/k_{in}^2]
    tppi = ppi                 !1/k_{out}^2

    nm = nmax

    do n1 = mm1, nmax
        k1 = n1 - mm1 + 1
        kk1 = k1 + nm
        do n2 = mm1, nmax
            k2 = n2 - mm1 + 1
            kk2 = k2 + nm

            tar12 = i12(n1, n2)
            tai12 = -r12(n1, n2)
            tgr12 = ig12(n1, n2)
            tgi12 = -rg12(n1, n2)

            tar21 = -i21(n1, n2)
            tai21 = r21(n1, n2)
            tgr21 = -ig21(n1, n2)
            tgi21 = rg21(n1, n2)

            tqr(k1, k2) = tpir*tar21 - tpii*tai21 + tppi*tar12
            tqi(k1, k2) = tpir*tai21 + tpii*tar21 + tppi*tai12
            trgqr(k1, k2) = tpir*tgr21 - tpii*tgi21 + tppi*tgr12
            trgqi(k1, k2) = tpir*tgi21 + tpii*tgr21 + tppi*tgi12

            tqr(k1, kk2) = 0d0
            tqi(k1, kk2) = 0d0
            trgqr(k1, kk2) = 0d0
            trgqi(k1, kk2) = 0d0

            tqr(kk1, k2) = 0d0
            tqi(kk1, k2) = 0d0
            trgqr(kk1, k2) = 0d0
            trgqi(kk1, k2) = 0d0

            tqr(kk1, kk2) = tpir*tar12 - tpii*tai12 + tppi*tar21
            tqi(kk1, kk2) = tpir*tai12 + tpii*tar12 + tppi*tai21
            trgqr(kk1, kk2) = tpir*tgr12 - tpii*tgi12 + tppi*tgr21
            trgqi(kk1, kk2) = tpir*tgi12 + tpii*tgr12 + tppi*tgi21
        end do
    end do

    nnmax = 2*nm
    do n1 = 1, nnmax
        do n2 = 1, nnmax
            qr(n1, n2) = tqr(n1, n2)
            qi(n1, n2) = tqi(n1, n2)
            rgqr(n1, n2) = trgqr(n1, n2)
            rgqi(n1, n2) = trgqi(n1, n2)
        end do
    end do
    !%%%%%%%%%%%%%%%%%%%%%%%  forming resulting t-matrix
    !
    ! calculate the product q^{-1} rg q
    !
    call tt(nmax)

    return
end

!*********************************************************************

subroutine tmatr0_leru(ngauss, x, w, an, ann, ppi, pir, pii, r, dr, ddr, &
        drr, dri, nmax, ncheck, naxsm)
    !--------/---------/---------/---------/---------/---------/---------/--
    ! >>> ngauss,x,w,an,ann,ppi,pir,pii,r,dr,ddr,drr,dri,nmax,ncheck
    ! <<< common blocks /tmat99/, /ct/ (for main),  and /ctt/ (for tt)
    !=====================
    ! complex(dp), parameter:: czero=cmplx_dp(0.0_dp,0.0_dp), ci=cmplx_dp(0.0_dp, 1.0_dp)
    ! are supplied here by the "constants" routine
    !
    !  determines the t-matrix of an axially symmetric scatterer
    !                           for m=0
    !
    !  ngauss - the number of gif division points
    !  x=\cos\theta  - gif division points
    !  w - gif weights
    !  an(n)=n*(n+1)
    !  ann(l_1,l_2)=\sqrt{\fr{(2 l_1+1)}{l_1(l_1+1)} }
    !                       \sqrt{\fr{(2 l_2+1)}{l_2(l_2+1)} }/2
    !  nmax - angular momentum cutoff
    !  ncheck  -  .eq.0  then  ngss=2*ngauss, factor=1d0
    !             .eq.1  then  ngss = ngauss, factor=2d0
    !  naxsm   -  .eq.0 : gauss abscissas do not have +/- theta symmetry
    !             .eq.1 : gauss abscissas have +/- theta symmetry
    !  MRR - the real part of the rel.  refractive index
    !  MRI - the imaginary of the rel. the refractive index
    !  p=dacos(-1d0)
    !  pi=p*2d0/lam      !wave vector k=2\pi/\ld in the host medium
    !  ppi=pi*pi         !k^2
    !  pir=ppi*mrr       !k^2*mrr
    !  pii=ppi*mri       !k^2*mri
    !  r=r^2(\theta)                        for axially symmetric particles
    !  dr=[dr(\theta)/(d\theta)]/r(\theta)  for axially symmetric particles
    !  ddr=\lambda/[2*\pi*r(\theta)]=1/(k_out*r)
    !  drr=(mrr/(mrr**2+mri**2))*(\lambda/[2*\pi*r(\theta)])
    !                  = re 1/(k_in*r)
    !  dri=-(mri/(mrr**2+mri**2))*(\lambda/[2*\pi*r(\theta)])
    !                  = im 1/(k_in*r)
    !  Refractive index outside is assumed to real, whereas inside
    !  a scatterer, refractive index is allowed to be complex in general.
    !  Consequently, the bessel function j_l(k_in*r) will in general
    !  be complex. The routine below performs Waterman surface integral
    !  separately for the real and imaginary parts of the integrand.
    !--------/---------/---------/---------/---------/---------/---------/--
    use libcylinder
    implicit none
    integer, intent(in) :: ngauss, nmax, ncheck, naxsm
    integer i, i1, i2, k1, k2, kk1, kk2, &
            mm1, n, n1, n2, ng, ngss, nm, nnmax
    real(dp), intent(in) :: ppi, pir, pii
    real(dp) d1n1, a12, a21, a22, ai12, ai21, an1, &
            an12, an2, ar12, ar21, b1i, b1r, b2i, b2r, b3i, b3r, b4i, b4r, &
            b5r, b5i, c1r, c1i, c2r, c2i, c3r, c3i, c4r, c4i, c5r, c5i, &
            d1n2, d2n1, d2n2, dd1, dd2, ddri, drii, drri, f1, f2, &
            factor, gi12, gi21, gr12, gr21, &
            qdj1, qdji2, qdjr2, qdy1, qj1, &
            qji2, qjr2, qy1, rri, si, tai12, tai21, &
            tar12, tar21, tgi12, tgi21, &
            tgr12, tgr21, tppi, tpii, tpir, uri

    real(dp), intent(in) :: an(npn1), ann(npn1, npn1), r(npng2), dr(npng2), &
            ddr(npng2), drr(npng2), dri(npng2), x(npng2), w(npng2)
    real(dp) sig(npn2), d1(npng2, npn1), d2(npng2, npn1), &
            rr(npng2), dv1(npn1), dv2(npn1)

    real(dp) r11(npn1, npn1), r12(npn1, npn1), &
            r21(npn1, npn1), r22(npn1, npn1), &
            i11(npn1, npn1), i12(npn1, npn1), &
            i21(npn1, npn1), i22(npn1, npn1), &
            rg11(npn1, npn1), rg12(npn1, npn1), &
            rg21(npn1, npn1), rg22(npn1, npn1), &
            ig11(npn1, npn1), ig12(npn1, npn1), &
            ig21(npn1, npn1), ig22(npn1, npn1), &
            qr(npn2, npn2), qi(npn2, npn2), &
            rgqr(npn2, npn2), rgqi(npn2, npn2), &
            tqr(npn2, npn2), tqi(npn2, npn2), &
            trgqr(npn2, npn2), trgqi(npn2, npn2)
    !c      real(dp) tr1(npn2,npn2),ti1(npn2,npn2)
    !Ru_patch:
    !real(dp)complex(dp)
    !integer m    m\equiv 0 in tmatr0
    real(dp) :: a11,drd, qs, xkr, A(nmax)
    complex(dp) :: znf, zxipsi, zxidpsi, zdxipsi, zdxidpsi, &
            zl1, zl2, zl3, zl4, zl1i, zl2i, zl3i, zl4i
    complex(dp) zq11(npn1, npn1), zq12(npn1, npn1), zq21(npn1, npn1), &
            zq22(npn1, npn1)
    !
    common /tmat99/&
            r11, r12, r21, r22, i11, i12, i21, i22, rg11, rg12, rg21, rg22, &
            ig11, ig12, ig21, ig22           !only between tmatr routines
    !c      common /ct/ tr1,ti1                       !output from tt routine
    common /ctt/ qr, qi, rgqr, rgqi              !input for tt routine
    common /ztt/ zq11, zq12, zq21, zq22
    !
    mm1 = 1     !even if m=0 one has still mm1=1
    nnmax = nmax + nmax
    ng = 2*ngauss
    ngss = ng
    factor = 1d0
    !
    if (ncheck==1) then         !theta=pi/2 is scatterer mirror symmetry plane
        ngss = ngauss
        factor = 2d0
    else                       !theta=pi/2 is not a scatterer mirror symmetry plane

    endif
    !
    si = 1d0
    do n = 1, nnmax
        si = -si
        sig(n) = si              !=(-1)**n
    end do
!
! assigning wigner d-matrices - assuming mirror symmetry
! in the \theta=\pi/2 plane:

    do i = 1, ngauss

        i1 = ngauss - i + 1
        !c         i2=ngauss+i
        if (ncheck==0) i2 = ngauss + i
        !
        call vig (x(i1), nmax, 0, dv1, dv2)
        !
        do n = 1, nmax

            dd1 = dv1(n)
            dd2 = dv2(n)
            d1(i1, n) = dd1
            d2(i1, n) = dd2

! if theta=pi/2 is not a scatterer mirror symmetry plane but
! gauss abscissas are still chosen symmetrically:

            if ((ncheck==0).and.(naxsm==1)) then
! using (4.2.4) and (4.2.6) of {ed},
!           d_{0m}^{(l)}(\pi-\theta) = (-1)^{l+m} d_{0m}^{(l)}(\theta)
! (4.2.5) of {ed}:                   = (-1)^{l} d_{0 -m}^{(l)}(\theta)

                si = sig(n)                  !=(-1)**n
                d1(i2, n) = dd1*si
                d2(i2, n) = -dd2*si

            end if
        enddo
! if neither scatterer nor gauss abscissas have theta=pi/2
! as a mirror symmetry plane:

        if ((ncheck==0).and.(naxsm==0)) then
            !
            call vig (x(i2), nmax, 0, dv1, dv2)
            !
            do n = 1, nmax
                dd1 = dv1(n)
                dd2 = dv2(n)
                d1(i2, n) = dd1
                d2(i2, n) = dd2
            enddo

        end if
    end do
!  Assigning r^2(\theta)*weight product and Bessel function independent factors:
!  In Ru et al formulas:
!  s=n_2/n_1$ is the refractive index constrast and $x(\theta)=k_1 r(\theta)$ is the
!  {\em polar} angle dependent size parameter of an axially symmetric scatterer.
!  assigning r^2(\theta)*weight product:

    do i = 1, ngss
        rr(i) = w(i)*r(i)           !r^2(\theta)*weight

        !c           if (dr(i).eq.0.d0) rr(i)=0.d0   !temporarily only
    end do
    !
    ! Before the n1, n2 (or l,l') loop:

    znf = cmplx_dp(pir, pii)/ppi       !in general complex ref. index contrast n_2/n_1;
    !the parameter 's' in Ru et al formulas

    do n=mm1,nmax
        A(n) = sqrt(dble(2*n + 1)/(2*n*(n + 1)))
    enddo

    do n1 = mm1, nmax
        an1 = an(n1)
        do n2 = mm1, nmax
            an2 = an(n2)
! Init:
            ar12 = 0d0
            ar21 = 0d0
            ai12 = 0d0
            ai21 = 0d0
            gr12 = 0d0
            gr21 = 0d0
            gi12 = 0d0
            gi21 = 0d0

! Init of Ru et al integrals:
            zl1 = czero
            zl2 = czero
            zl3 = czero
            zl4 = czero

!        open(nout+3,file='surfint.dat')   !gauss convergence check

            if (ncheck/=1.or.sig(n1 + n2)>=0d0) then

! Gauss integration loop:
! Traditional integration transforms
! \int_0^\pi X* \sin\theta d\theta --> \int_{-1}^1 X* d(\cos\theta)
! and generates weights for the independent variable "\cos\theta"
        do i = 1, ngss    !=ngauss   if ncheck.eq.1
            !=2*ngauss if ncheck.eq.0

            d1n1 = d1(i, n1)
            d2n1 = d2(i, n1)
            d1n2 = d1(i, n2)
            d2n2 = d2(i, n2)

!aij store important products d_n.d_{n'}, d_n.\tau_{n'}, \tau_n.d_{n'}
!necessary for forming the integrands of Q-matrices.
!The first index of aij is 1 (2) if the factor with n1-index is d1 (d2)
!the second index of aij is 1 (2) if the factor with n2-index is d1 (d2)
            a11 = d1n1*d1n2     !for LeRu
            a12 = d1n1*d2n2
            a21 = d2n1*d1n2
            a22 = d2n1*d2n2
            !aa1=a12+a21

! Spherical Bessel functions:
! Since refractive index is allowed to be complex in general,
! the Bessel function j_l(k_in*r) is complex. The code below
! performs a separation of the complex integrand in Waterman's
! surface integral into its respective real and imaginary
! parts.

! Bessel functions of the exterior argument:

            qj1 = cbess%j(i, n1)
            qy1 = cbess%y(i, n1)
            qdj1 = cbess%dj(i, n1)
            qdy1 = cbess%dy(i, n1)

! Bessel functions of the interior argument:
            qjr2 = cbess%jr(i, n2)
            qji2 = cbess%ji(i, n2)
            qdjr2 = cbess%djr(i, n2)
            qdji2 = cbess%dji(i, n2)

!--------/---------/---------/---------/---------/---------/---------/--
!Assigning integrads by eqs 2-7 of Ru2012
!Bessel function part:
!Amend for the Riccati-Bessel functions of Ru et al.
!The Riccati-Bessel functions \xi_l(z)=zh_l(z) (\psi_l(z)=zj_l(z)) are the functions of
!the exterior (interior) argument there.
! \xi or d\xi carry always the subscript n (=n1 here) and
! \psi or d\psi carry always the subscript k (=n2 here)
!The latter allows one to initiate their complex valued products.
!  PI=P*2D0/LAM                 !wave vector k=2\pi/\ld in the host medium
!  PPI=PI*PI                    !k^2
!  PIR=PPI*MRR                  !k^2*mrr
!  PII=PPI*MRI                  !k^2*mri
!  R(I):=r^2(\theta)

        xkr = sqrt(ppi*r(i))        !k_out*r(\theta)

        zxipsi = cmplx_dp(pir, pii)*cmplx_dp(qj1, qy1)*cmplx_dp(qjr2, qji2)*r(i)

        zxidpsi = xkr*cmplx_dp(qj1, qy1)*cmplx_dp(qjr2, qji2) + &
                cmplx_dp(pir, pii)*cmplx_dp(qj1, qy1)*cmplx_dp(qdjr2, qdji2)*r(i)

        zdxipsi = znf*xkr*cmplx_dp(qj1, qy1)*cmplx_dp(qjr2, qji2) + &
                cmplx_dp(pir, pii)*cmplx_dp(qdj1, qdy1)*cmplx_dp(qjr2, qji2)*r(i)

        zdxidpsi = cmplx_dp(qj1, qy1)*cmplx_dp(qjr2, qji2) + &
                znf*xkr*cmplx_dp(qj1, qy1)*cmplx_dp(qdjr2, qdji2) + &
                xkr*cmplx_dp(qdj1, qdy1)*cmplx_dp(qjr2, qji2) + &
                cmplx_dp(pir, pii)*cmplx_dp(qdj1, qdy1)*cmplx_dp(qdjr2, qdji2)*r(i)
!--------/---------/---------/---------/---------/---------/---------/--
! Gauss integration weigths and variables are in the variable '\cos\theta'.
! At the present code of Ref. \ct{MTL}:
! r(i) contains $r^2(\theta)$, whereas we need $r(\theta)$
! dr(i) contains $r'(\theta)/r(\theta)$, whereas we need $r'(\theta)$
! ds(i) = s(i)*qm*wr       !=dble(m)*w(i)*r^2(\theta)/(|\sin\theta|)

!       qs = dsqrt(1.d0 - x(i)**2)          !\sin\theta
        drd = xkr*dr(i)                     !x_\theta=k_out*r'(\theta) of Ru et al

! Amend for the Riccati-Bessel functions of Ru et al
! on using that they only contain combinations of
! \xi or d\xi with the subscript n (=n1 here) and
! \psi or d\psi with the subscript k (=n2 here)
! Assigning integrads by eqs 2-7 of Ru2012:

            if (ncheck==1) then         !theta=pi/2 is scatterer mirror symmetry plane
                !eqs 18-19 of Ru2012
                if (mod(n1 + n2, 2)==0) then     !n1+n2 even

                    zl1i = a21*drd*zxipsi
                    zl2i = a12*drd*zxipsi
                    zl3i = a21*drd*zdxidpsi &
                            - dble(n1*(n1 + 1))*a11*zxidpsi
                    zl4i = a12*znf*drd*zdxidpsi &
                            - dble(n2*(n2 + 1))*a11*zdxipsi

                else                    !n1+n2 odd

                    zl1i = czero
                    zl2i = czero
                    zl3i = czero
                    zl4i = czero

                end if

            else if (ncheck==0) then        !theta=pi/2 is not a scatterer mirror symmetry plane

                zl1i = a21*drd*zxipsi
                zl2i = a12*drd*zxipsi
                zl3i = a21*drd*zdxidpsi &
                        - dble(n1*(n1 + 1))*a11*zxidpsi
                zl4i = a12*znf*drd*zdxidpsi &
                        - dble(n2*(n2 + 1))*a11*zdxipsi
!--------/---------/---------/---------/---------/---------/---------/--

            end if

! Integrate Ru et al integrals:
            zl1 = zl1 + zl1i*w(i)
            zl2 = zl2 + zl2i*w(i)
            zl3 = zl3 + zl3i*w(i)
            zl4 = zl4 + zl4i*w(i)

! todo: a switch can be implemented here to run the code
! either by following Ru et al integration or
! following the original MTL path
!__________________
! re and im of j_{n2}(k_{in}r) j_{n1}(k_{out}r):

                    c1r = qjr2*qj1
                    c1i = qji2*qj1
                    ! re and im of j_{n2}(k_{in}r) h_{n1}(k_{out}r):

                    b1r = c1r - qji2*qy1
                    b1i = c1i + qjr2*qy1
                    ! re and im of j_{n2}(k_{in}r) [k_{out}r j_{n1}(k_{out}r)]'/(k_{out}r):

                    c2r = qjr2*qdj1
                    c2i = qji2*qdj1
                    ! re and im of j_{n2}(k_{in}r) [k_{out}r h_{n1}(k_{out}r)]'/(k_{out}r):

                    b2r = c2r - qji2*qdy1
                    b2i = c2i + qjr2*qdy1

                    ddri = ddr(i)               !1/(k_{out}r)
                    ! re and im of [1/(k_{out}r)]*j_{n2}(k_{in}r) j_{n1}(k_{out}r)

                    c3r = ddri*c1r
                    c3i = ddri*c1i
                    ! re and im of [1/(k_{out}r)]*j_{n2}(k_{in}r) h_{n1}(k_{out}r):

                    b3r = ddri*b1r
                    b3i = ddri*b1i
                    ! re and im of  [k_{in}r j_{n2}(k_{in}r)]'/(k_{in}r)
                    !                         *j_{n1}(k_{out}r):

                    c4r = qdjr2*qj1
                    c4i = qdji2*qj1
                    ! re and im of [k_{in}r j_{n2}(k_{in}r)]'/(k_{in}r)
                    !                        *h_{n1}(k_{out}r):

                    b4r = c4r - qdji2*qy1
                    b4i = c4i + qdjr2*qy1

                    drri = drr(i)               !re[1/(k_{in}r)]
                    drii = dri(i)               !im[1/(k_{in}r)]
                    ! re and im of [1/(k_{in}r)] j_{n2}(k_{in}r) j_{n1}(k_{out}r):

                    c5r = c1r*drri - c1i*drii
                    c5i = c1i*drri + c1r*drii
                    ! re and im of [1/(k_{in}r)] j_{n2}(k_{in}r) h_{n1}(k_{out}r):

                    b5r = b1r*drri - b1i*drii
                    b5i = b1i*drri + b1r*drii
                    !%%%%%%%  forming integrands of j-matrices (j^{11}=j^{22}=0 for m=0): %%%%%%%%

                    uri = dr(i)        !dr/(d\theta)
                    rri = rr(i)        !w(i)*r^2(\theta)
                    ! w(i)*r^2(\theta)*d2n1*d2n2:
                    f1 = rri*a22      !prefactor containing r^2(\theta)<->hat{r} part
                    ! n1*(n1+1)*w(i)*r(\theta)*[dr/(d\theta)]*d1n1*d2n2:
                    f2 = rri*uri*an1*a12     !prefactor containing r(\theta)*[dr/(d\theta)]
                    !                                          !hat{theta} part

                    ar12 = ar12 + f1*b2r + f2*b3r        !~re j^{12}
                    ai12 = ai12 + f1*b2i + f2*b3i        !~im j^{12}

                    gr12 = gr12 + f1*c2r + f2*c3r        !~re rg j^{12}
                    gi12 = gi12 + f1*c2i + f2*c3i        !~im rg j^{12}

                    !*  n2*(n2+1)*w(i)*r(\theta)*[dr/(d\theta)]*d2n1*d1n2:
                    f2 = rri*uri*an2*a21     !prefactor containing r(\theta)*[dr/(d\theta)]
                    !                                          !hat{theta} part

                    ar21 = ar21 + f1*b4r + f2*b5r        !~re j^{21}
                    ai21 = ai21 + f1*b4i + f2*b5i        !~im j^{21}

                    gr21 = gr21 + f1*c4r + f2*c5r        !~re rg j^{21}
                    gi21 = gi21 + f1*c4i + f2*c5i        !~im rg j^{21}
                end do               !end of gauss integration

                ! Taking into account the "factor":
                ! factor=1 or 2 depending on ncheck==0 or ncheck==1

                if (abs(factor-2d0)<=1.d-10) then
                    zl1 = zl1*factor
                    zl2 = zl2*factor
                    zl3 = zl3*factor
                    zl4 = zl4*factor
                end if

!After performing surface integrals, forming the Q-matrices
!by eqs 12-17 of Ru2012

        if (n1/=n2) then
            zq11(n1, n2) = ci*A(n1)*A(n2)*(znf**2 - cone)*&
                    (n1*(n1 + 1)*zl2 - n2*(n2 + 1)*zl1)/((n1*(n1 + 1) - n2*(n2 + 1))*znf)
            zq22(n1, n2) = ci*A(n1)*A(n2)*(znf**2 - cone)*&
                    (zl3 + znf*n1*(n1 + 1)*(zl2 - zl1)/(n1*(n1 + 1) - n2*(n2 + 1)))/znf
        else                !diagonal elements
            zq11(n1, n2) = -ci*A(n1)*A(n2)*(-znf*zl1 + zl3 + (zl2 - zl4)/znf)
            zq22(n1, n2) = -ci*A(n1)*A(n2)*(-zl1 + zl2 - zl4 + zl3/znf)
        end if

        zq12(n1, n2) = czero
        zq21(n1, n2) = czero

! todo: a switch can be implemented here to run the code
! either by following Ru et al integration or
! following the original MTL path
!Original MTL part continues here:

                !write(nout+3,*)'n1=',n1,'   n2=',n2
                !write(nout+3,*)'ar12=', ar12
                !write(nout+3,*)'ai12=', ai12
                !write(nout+3,*)'ar21=', ar21
                !write(nout+3,*)'ai21=', ai21
                !write(nout+3,*)'gr12=', gr12
                !write(nout+3,*)'gi12=', gi12
                !write(nout+3,*)'gr21=', gr21
                !write(nout+3,*)'gi21=', gi21
                !%%%%%%%%%%%%%  forming j-matrices (j^{11}=j^{22}=0 for m=0):
            end if

            an12 = ann(n1, n2)*factor

            r12(n1, n2) = ar12*an12
            r21(n1, n2) = ar21*an12
            i12(n1, n2) = ai12*an12
            i21(n1, n2) = ai21*an12

            rg12(n1, n2) = gr12*an12
            rg21(n1, n2) = gr21*an12
            ig12(n1, n2) = gi12*an12
            ig21(n1, n2) = gi21*an12
        end do
    end do            !end of the loop over angular momenta

    !      close(nout+3)
    !%%%%%%%%%%%%%%%%%%%%%%%  forming q and rgq -matrices

    tpir = pir                 !re [1/k_{in}^2]
    tpii = pii                 !im [1/k_{in}^2]
    tppi = ppi                 !1/k_{out}^2

    nm = nmax

    do n1 = mm1, nmax
        k1 = n1 - mm1 + 1
        kk1 = k1 + nm
        do n2 = mm1, nmax
            k2 = n2 - mm1 + 1
            kk2 = k2 + nm

            tar12 = i12(n1, n2)
            tai12 = -r12(n1, n2)
            tgr12 = ig12(n1, n2)
            tgi12 = -rg12(n1, n2)

            tar21 = -i21(n1, n2)
            tai21 = r21(n1, n2)
            tgr21 = -ig21(n1, n2)
            tgi21 = rg21(n1, n2)

            tqr(k1, k2) = tpir*tar21 - tpii*tai21 + tppi*tar12
            tqi(k1, k2) = tpir*tai21 + tpii*tar21 + tppi*tai12
            trgqr(k1, k2) = tpir*tgr21 - tpii*tgi21 + tppi*tgr12
            trgqi(k1, k2) = tpir*tgi21 + tpii*tgr21 + tppi*tgi12

            tqr(k1, kk2) = 0d0
            tqi(k1, kk2) = 0d0
            trgqr(k1, kk2) = 0d0
            trgqi(k1, kk2) = 0d0

            tqr(kk1, k2) = 0d0
            tqi(kk1, k2) = 0d0
            trgqr(kk1, k2) = 0d0
            trgqi(kk1, k2) = 0d0

            tqr(kk1, kk2) = tpir*tar12 - tpii*tai12 + tppi*tar21
            tqi(kk1, kk2) = tpir*tai12 + tpii*tar12 + tppi*tai21
            trgqr(kk1, kk2) = tpir*tgr12 - tpii*tgi12 + tppi*tgr21
            trgqi(kk1, kk2) = tpir*tgi12 + tpii*tgr12 + tppi*tgi21
        end do
    end do

    nnmax = 2*nm
    do n1 = 1, nnmax
        do n2 = 1, nnmax
            qr(n1, n2) = tqr(n1, n2)
            qi(n1, n2) = tqi(n1, n2)
            rgqr(n1, n2) = trgqr(n1, n2)
            rgqi(n1, n2) = trgqi(n1, n2)
        end do
    end do
!%%%%%%%%%%%%%%%%%%%%%%%  forming resulting t-matrix
!
! calculate the product q^{-1} rg q
!
    call tt(nmax)

    return
end


!**********************************************************************
subroutine tmatr(m, ngauss, x, w, an, ann, s, ss, ppi, pir, pii, r, dr, ddr, &
        drr, dri, nmax, ncheck, naxsm)
    !--------/---------/---------/---------/---------/---------/---------/--
    ! >>> ngauss,x,w,an,ann,s,ss,ppi,pir,pii,r,dr,ddr,drr,dri,nmax,ncheck
    ! <<< common blocks /tmat99/, /ct/ (for main),  and /ctt/ (for tt)
    !=====================
    !
    !  determines the t-matrix of an axially symmetric scatterer
    !                           for m > 0
    !
    !  m      - azimuthal number
    !  ngauss - the number of gif division points
    !  x=\cos\theta  - gif division points
    !  w - gif weights
    !  an(n)=n*(n+1)
    !  ann(l_1,l_2)=\sqrt{\fr{(2 l_1+1)}{l_1(l_1+1)} }
    !                       \sqrt{\fr{(2 l_2+1)}{l_2(l_2+1)} }/2
    !  nmax - angular momentum cutoff
    !  ncheck  -  .eq.0  then  ngss=2*ngauss, factor=1d0
    !             .eq.1  then  ngss = ngauss, factor=2d0
    !  naxsm   -  .eq.0 : gauss abscissas do not have +/- theta symmetry
    !             .eq.1 : gauss abscissas have +/- theta symmetry
    !  ncheck - specifies whether ng=2*ngauss or otherwise
    !  p=dacos(-1d0)
    !  pi=p*2d0/lam - wave vector
    !  ppi=pi*pi
    !  pir=ppi*mrr
    !  pii=ppi*mri
    !  r=r^2(\theta)                       for axially symmetric particles
    !  dr=[dr(\theta)/(d\theta)]/r(\theta) for axially symmetric particles
    !  ddr=\lambda/[2*\pi*r(\theta)]=1/(k_out*r)
    !  drr=(mrr/(mrr**2+mri**2))*(\lambda/[2*\pi*r(\theta)])
    !                  = re 1/(k_in*r)
    !  dri=-(mri/(mrr**2+mri**2))*(\lambda/[2*\pi*r(\theta)])
    !                  = im 1/(k_in*r)
    !  refractive index outside is assumed to real, whereas inside
    !  a scatterer, refractive index is allowed to be complex in general.
    !  consequently, the bessel function j_l(k_in*r) will in general
    !  be complex. the routine below performs waterman surface integral
    !  separately for the real and imaginary parts of the integrand.
    !
    !--------/---------/---------/---------/---------/---------/---------/--
    use libcylinder
    !include 'ampld.par.f'
    implicit none
    integer m, ngauss, nmax, ncheck, naxsm
    integer n, ng, ngss, nm, nnmax, &
            n1, n2, i, i1, i2, &
            k1, k2, kk1, kk2, mm1
    real(dp) ppi, pir, pii
    real(dp) a11, a12, a21, a22, aa1, aa2, &
            ar11, ar12, ar21, ar22, &
            ai11, ai12, ai21, ai22, &
            an1, an2, an12, &
            b1r, b1i, b2r, b2i, b3r, b3i, &
            b4r, b4i, b5r, b5i, b6r, b6i, &
            b7r, b7i, b8r, b8i, &
            c1r, c1i, c2r, c2i, c3r, c3i, &
            c4r, c4i, c5r, c5i, c6r, c6i, &
            c7r, c7i, c8r, c8i, &
            d1n1, d1n2, d2n1, d2n2, &
            dd1, dd2, ddri, drii, drri, &
            dsi, e1, e2, e3, f1, f2, &
            factor, &
            gr11, gr12, gr21, gr22, &
            gi11, gi12, gi21, gi22, &
            qdj1, qdjr2, qdji2, qdy1, &
            qj1, qjr2, qji2, qm, qmm, &
            qy1, rri, si, &
            tar11, tar12, tar21, tar22, &
            tai11, tai12, tai21, tai22, &
            tgr11, tgr12, tgr21, tgr22, &
            tgi11, tgi12, tgi21, tgi22, &
            tpir, tpii, tppi, uri, wr
    real(dp)  x(npng2), w(npng2), an(npn1), s(npng2), ss(npng2), &
            r(npng2), dr(npng2), sig(npn2), &
            ddr(npng2), drr(npng2), &
            d1(npng2, npn1), d2(npng2, npn1), &
            dri(npng2), ds(npng2), dss(npng2), rr(npng2), &
            dv1(npn1), dv2(npn1)
    real(dp)  r11(npn1, npn1), r12(npn1, npn1), &
            r21(npn1, npn1), r22(npn1, npn1), &
            i11(npn1, npn1), i12(npn1, npn1), &
            i21(npn1, npn1), i22(npn1, npn1), &
            rg11(npn1, npn1), rg12(npn1, npn1), &
            rg21(npn1, npn1), rg22(npn1, npn1), &
            ig11(npn1, npn1), ig12(npn1, npn1), &
            ig21(npn1, npn1), ig22(npn1, npn1), &
            ann(npn1, npn1), &
            qr(npn2, npn2), qi(npn2, npn2), &
            rgqr(npn2, npn2), rgqi(npn2, npn2), &
            tqr(npn2, npn2), tqi(npn2, npn2), &
            trgqr(npn2, npn2), trgqi(npn2, npn2)
    !c      real(dp) tr1(npn2,npn2),ti1(npn2,npn2)
    !________
    common /tmat99/&
            r11, r12, r21, r22, i11, i12, i21, i22, rg11, rg12, rg21, rg22, &
            ig11, ig12, ig21, ig22          !only between tmatr routines
    !c      common /ct/ tr1,ti1                      !output from tt routine
    common /ctt/ qr, qi, rgqr, rgqi             !input for tt routine
    !________
    mm1 = m
    qm = dble(m)
    qmm = qm*qm
    ng = 2*ngauss
    ngss = ng
    factor = 1d0
    !
    if (ncheck==1) then          !theta=pi/2 is mirror symmetry plane
        ngss = ngauss
        factor = 2d0
    endif
    !
    si = 1d0
    nm = nmax + nmax
    do n = 1, nm
        si = -si
        sig(n) = si              !=(-1)**n
    end do
    !
    ! assigning wigner d-matrices - assuming mirror symmetry
    ! in the \theta=\pi/2 plane:
    do i = 1, ngauss
        i1 = ngauss - i + 1
        i2 = ngauss + i
        !
        call vig (x(i1), nmax, m, dv1, dv2)
        !
        do n = 1, nmax
            dd1 = dv1(n)
            dd2 = dv2(n)
            d1(i1, n) = dd1
            d2(i1, n) = dd2
            if (naxsm==1) then         !gauss abscissas chosen +/- symmetric
                ! using (4.2.4) and (4.2.6) of {ed},
                !           d_{0m}^{(l)}(\pi-\theta) = (-1)^{l+m} d_{0m}^{(l)}(\theta)
                si = sig(n + m)                  !=(-1)**(n+m)
                !                                          !exactly what follows from {ed}
                d1(i2, n) = dd1*si
                d2(i2, n) = -dd2*si
            end if
        enddo
        !
        if (naxsm==0) then        !gauss abscissas not chosen +/- symmetric
            !
            call vig (x(i2), nmax, m, dv1, dv2)
            !
            do n = 1, nmax
                dd1 = dv1(n)
                dd2 = dv2(n)
                d1(i2, n) = dd1
                d2(i2, n) = dd2
            enddo
        end if
    end do
    !
    !  assigning r^2(\theta)*weight product:
    do i = 1, ngss
        wr = w(i)*r(i)
        !c           if (dr(i).eq.0.d0) wr=0.d0   !temporarily only
        ds(i) = s(i)*qm*wr       !=dble(m)*w(i)*r^2(\theta)/(|\sin\theta|)
        dss(i) = ss(i)*qmm       !=dble(m)**2/(\sin^2\theta)
        rr(i) = wr
    end do
    !
    do n1 = mm1, nmax
        an1 = an(n1)
        do n2 = mm1, nmax
            an2 = an(n2)
            ar11 = 0d0
            ar12 = 0d0
            ar21 = 0d0
            ar22 = 0d0
            ai11 = 0d0
            ai12 = 0d0
            ai21 = 0d0
            ai22 = 0d0
            gr11 = 0d0
            gr12 = 0d0
            gr21 = 0d0
            gr22 = 0d0
            gi11 = 0d0
            gi12 = 0d0
            gi21 = 0d0
            gi22 = 0d0
            si = sig(n1 + n2)
            do i = 1, ngss
                d1n1 = d1(i, n1)
                d2n1 = d2(i, n1)
                d1n2 = d1(i, n2)
                d2n2 = d2(i, n2)
                a11 = d1n1*d1n2
                a12 = d1n1*d2n2
                a21 = d2n1*d1n2
                a22 = d2n1*d2n2
                aa1 = a12 + a21            != d1n1*d2n2+d2n1*d1n2
                aa2 = a11*dss(i) + a22     !=(d1n1*d1n2)*dble(m)**2/(\sin^2\theta)
                !                          ! +d2n1*d2n2
                !  Bessel functions:
                !  since refractive index is allowed to be complex in general,
                !  the bessel function j_l(k_in*r) is complex. the code below
                !  performs a separation of the complex integrand in waterman's
                !  surface integral into its respective real and imaginary
                !  parts:
                qj1 = cbess%j(i, n1)
                qy1 = cbess%y(i, n1)
                qjr2 = cbess%jr(i, n2)
                qji2 = cbess%ji(i, n2)
                qdjr2 = cbess%djr(i, n2)
                qdji2 = cbess%dji(i, n2)
                qdj1 = cbess%dj(i, n1)
                qdy1 = cbess%dy(i, n1)
                ! re and im of j_{n2}(k_{in}r) j_{n1}(k_{out}r):
                c1r = qjr2*qj1
                c1i = qji2*qj1
                ! re and im of j_{n2}(k_{in}r) h_{n1}(k_{out}r):
                b1r = c1r - qji2*qy1
                b1i = c1i + qjr2*qy1
                ! re and im of j_{n2}(k_{in}r) j_{n1}'(k_{out}r):
                c2r = qjr2*qdj1
                c2i = qji2*qdj1
                ! re and im of j_{n2}(k_{in}r) h_{n1}'(k_{out}r):
                b2r = c2r - qji2*qdy1
                b2i = c2i + qjr2*qdy1
                ddri = ddr(i)               !1/(k_{out}r)
                ! re and im of [1/(k_{out}r)]*j_{n2}(k_{in}r) j_{n1}(k_{out}r)
                c3r = ddri*c1r
                c3i = ddri*c1i
                ! re and im of [1/(k_{out}r)]*j_{n2}(k_{in}r) h_{n1}(k_{out}r):
                b3r = ddri*b1r
                b3i = ddri*b1i
                ! re and im of j_{n2}'(k_{in}r) j_{n1}(k_{out}r):
                c4r = qdjr2*qj1
                c4i = qdji2*qj1
                ! re and im of j_{n2}'(k_{in}r) h_{n1}(k_{out}r):
                b4r = c4r - qdji2*qy1
                b4i = c4i + qdjr2*qy1
                drri = drr(i)               !re[1/(k_{in}r)]
                drii = dri(i)               !im[1/(k_{in}r)]
                ! re and im of [1/(k_{in}r)] j_{n2}(k_{in}r) j_{n1}(k_{out}r):
                c5r = c1r*drri - c1i*drii
                c5i = c1i*drri + c1r*drii
                ! re and im of [1/(k_{in}r)] j_{n2}(k_{in}r) h_{n1}(k_{out}r):
                b5r = b1r*drri - b1i*drii
                b5i = b1i*drri + b1r*drii
                ! re and im of j_{n2}'(k_{in}r) j_{n1}'(k_{out}r):
                c6r = qdjr2*qdj1
                c6i = qdji2*qdj1
                ! re and im of j_{n2}'(k_{in}r) h_{n1}'(k_{out}r):
                b6r = c6r - qdji2*qdy1
                b6i = c6i + qdjr2*qdy1
                ! re and im of [1/(k_{out}r)] j_{n2}'(k_{in}r) j_{n1}(k_{out}r):
                c7r = c4r*ddri
                c7i = c4i*ddri
                ! re and im of [1/(k_{out}r)] j_{n2}'(k_{in}r) h_{n1}(k_{out}r):
                b7r = b4r*ddri
                b7i = b4i*ddri
                ! re and im of [1/(k_{in}r)] j_{n2}(k_{in}r) j_{n1}'(k_{out}r):
                c8r = c2r*drri - c2i*drii
                c8i = c2i*drri + c2r*drii
                ! re and im of [1/(k_{in}r)] j_{n2}(k_{in}r) h_{n1}'(k_{out}r):
                b8r = b2r*drri - b2i*drii
                b8i = b2i*drri + b2r*drii
                ! %%%%%%%%%  forming integrands of j-matrices (j^{11}=j^{22}=0 for m=0):
                uri = dr(i)
                dsi = ds(i)
                rri = rr(i)
                if (ncheck==1.and.si>0d0) then
                    ! w(i)*r^2(\theta)*[(d1n1*d1n2)*dble(m)**2/(\sin^2\theta)+d2n1*d2n2]:
                    f1 = rri*aa2            !prefactor containing r^2(\theta)<->hat{r} part
                    ! n1*(n1+1)*w(i)*r(\theta)*[dr/(d\theta)]*d1n1*d2n2:
                    f2 = rri*uri*an1*a12     !prefactor containing r(\theta)*[dr/(d\theta)]
                    !                                          !hat{theta} part
                    ar12 = ar12 + f1*b2r + f2*b3r        !~re j^{12}
                    ai12 = ai12 + f1*b2i + f2*b3i        !~im j^{12}
                    gr12 = gr12 + f1*c2r + f2*c3r        !~re rg j^{12}
                    gi12 = gi12 + f1*c2i + f2*c3i        !~im rg j^{12}
                    ! n2*(n2+1)*w(i)*r(\theta)*[dr/(d\theta)]*d2n1*d1n2:
                    f2 = rri*uri*an2*a21     !prefactor containing r(\theta)*[dr/(d\theta)]
                    !                                          !hat{theta} part
                    ar21 = ar21 + f1*b4r + f2*b5r
                    ai21 = ai21 + f1*b4i + f2*b5i
                    gr21 = gr21 + f1*c4r + f2*c5r
                    gi21 = gi21 + f1*c4i + f2*c5i
                    cycle
                end if
                ! [dble(m)*w(i)*r^2(i)/(|\sin\theta|)]*(d1n1*d2n2+d2n1*d1n2):
                e1 = dsi*aa1
                ar11 = ar11 + e1*b1r
                ai11 = ai11 + e1*b1i
                gr11 = gr11 + e1*c1r
                gi11 = gi11 + e1*c1i
                e2 = dsi*uri*a11
                e3 = e2*an2
                e2 = e2*an1
                ar22 = ar22 + e1*b6r + e2*b7r + e3*b8r
                ai22 = ai22 + e1*b6i + e2*b7i + e3*b8i
                gr22 = gr22 + e1*c6r + e2*c7r + e3*c8r
                gi22 = gi22 + e1*c6i + e2*c7i + e3*c8i
            end do  ! end of Gauss integration
            an12 = ann(n1, n2)*factor
            r11(n1, n2) = ar11*an12
            r12(n1, n2) = ar12*an12
            r21(n1, n2) = ar21*an12
            r22(n1, n2) = ar22*an12
            i11(n1, n2) = ai11*an12
            i12(n1, n2) = ai12*an12
            i21(n1, n2) = ai21*an12
            i22(n1, n2) = ai22*an12
            rg11(n1, n2) = gr11*an12
            rg12(n1, n2) = gr12*an12
            rg21(n1, n2) = gr21*an12
            rg22(n1, n2) = gr22*an12
            ig11(n1, n2) = gi11*an12
            ig12(n1, n2) = gi12*an12
            ig21(n1, n2) = gi21*an12
            ig22(n1, n2) = gi22*an12
        end do
    end do
    tpir = pir                 !re [1/k_{in}^2]
    tpii = pii                 !im [1/k_{in}^2]
    tppi = ppi                 !1/k_{out}^2
    nm = nmax - mm1 + 1
    do n1 = mm1, nmax
        k1 = n1 - mm1 + 1
        kk1 = k1 + nm
        do n2 = mm1, nmax
            k2 = n2 - mm1 + 1
            kk2 = k2 + nm
            tar11 = -r11(n1, n2)
            tai11 = -i11(n1, n2)
            tgr11 = -rg11(n1, n2)
            tgi11 = -ig11(n1, n2)
            tar12 = i12(n1, n2)
            tai12 = -r12(n1, n2)
            tgr12 = ig12(n1, n2)
            tgi12 = -rg12(n1, n2)
            tar21 = -i21(n1, n2)
            tai21 = r21(n1, n2)
            tgr21 = -ig21(n1, n2)
            tgi21 = rg21(n1, n2)
            tar22 = -r22(n1, n2)
            tai22 = -i22(n1, n2)
            tgr22 = -rg22(n1, n2)
            tgi22 = -ig22(n1, n2)
            tqr(k1, k2) = tpir*tar21 - tpii*tai21 + tppi*tar12
            tqi(k1, k2) = tpir*tai21 + tpii*tar21 + tppi*tai12
            trgqr(k1, k2) = tpir*tgr21 - tpii*tgi21 + tppi*tgr12
            trgqi(k1, k2) = tpir*tgi21 + tpii*tgr21 + tppi*tgi12
            tqr(k1, kk2) = tpir*tar11 - tpii*tai11 + tppi*tar22
            tqi(k1, kk2) = tpir*tai11 + tpii*tar11 + tppi*tai22
            trgqr(k1, kk2) = tpir*tgr11 - tpii*tgi11 + tppi*tgr22
            trgqi(k1, kk2) = tpir*tgi11 + tpii*tgr11 + tppi*tgi22
            tqr(kk1, k2) = tpir*tar22 - tpii*tai22 + tppi*tar11
            tqi(kk1, k2) = tpir*tai22 + tpii*tar22 + tppi*tai11
            trgqr(kk1, k2) = tpir*tgr22 - tpii*tgi22 + tppi*tgr11
            trgqi(kk1, k2) = tpir*tgi22 + tpii*tgr22 + tppi*tgi11
            tqr(kk1, kk2) = tpir*tar12 - tpii*tai12 + tppi*tar21
            tqi(kk1, kk2) = tpir*tai12 + tpii*tar12 + tppi*tai21
            trgqr(kk1, kk2) = tpir*tgr12 - tpii*tgi12 + tppi*tgr21
            trgqi(kk1, kk2) = tpir*tgi12 + tpii*tgr12 + tppi*tgi21
        end do
    end do
    nnmax = 2*nm
    do n1 = 1, nnmax
        do n2 = 1, nnmax
            qr(n1, n2) = tqr(n1, n2)
            qi(n1, n2) = tqi(n1, n2)
            rgqr(n1, n2) = trgqr(n1, n2)
            rgqi(n1, n2) = trgqi(n1, n2)
        end do
    end do
    call tt(nm)
    !
    return
end

!**********************************************************************

subroutine tmatr_adapt(m, ngauss, x, w, an, ann, s, ss, ppi, pir, pii, r, dr, ddr, &
        drr, dri, nmax, ncheck, naxsm)
    !--------/---------/---------/---------/---------/---------/---------/--
    ! >>> ngauss,x,w,an,ann,s,ss,ppi,pir,pii,r,dr,ddr,drr,dri,nmax,ncheck
    ! <<< common blocks /tmat99/, /ct/ (for main),  and /ctt/ (for tt)
    !=====================
    !
    !  determines the t-matrix of an axially symmetric scatterer
    !                           for m > 0
    !
    !  m      - azimuthal number
    !  ngauss - the number of gif division points
    !  x=\cos\theta  - gif division points
    !  w - gif weights
    !  an(n)=n*(n+1)
    !  ann(l_1,l_2)=\sqrt{\fr{(2 l_1+1)}{l_1(l_1+1)} }
    !                       \sqrt{\fr{(2 l_2+1)}{l_2(l_2+1)} }/2
    !  nmax - angular momentum cutoff
    !  ncheck  -  .eq.0  then  ngss=2*ngauss, factor=1d0
    !             .eq.1  then  ngss = ngauss, factor=2d0
    !  naxsm   -  .eq.0 : gauss abscissas do not have +/- theta symmetry
    !             .eq.1 : gauss abscissas have +/- theta symmetry
    !  ncheck - specifies whether ng=2*ngauss or otherwise
    !  p=dacos(-1d0)
    !  pi=p*2d0/lam - wave vector
    !  ppi=pi*pi
    !  pir=ppi*mrr
    !  pii=ppi*mri
    !  r=r^2(\theta)                       for axially symmetric particles
    !  dr=[dr(\theta)/(d\theta)]/r(\theta) for axially symmetric particles
    !  ddr=\lambda/[2*\pi*r(\theta)]=1/(k_out*r)
    !  drr=(mrr/(mrr**2+mri**2))*(\lambda/[2*\pi*r(\theta)])
    !                  = re 1/(k_in*r)
    !  dri=-(mri/(mrr**2+mri**2))*(\lambda/[2*\pi*r(\theta)])
    !                  = im 1/(k_in*r)
    !  refractive index outside is assumed to real, whereas inside
    !  a scatterer, refractive index is allowed to be complex in general.
    !  consequently, the bessel function j_l(k_in*r) will in general
    !  be complex. the routine below performs waterman surface integral
    !  separately for the real and imaginary parts of the integrand.
    !
    !--------/---------/---------/---------/---------/---------/---------/--
    use libcylinder
    use lib_adaptive_integration
    !include 'ampld.par.f'
    implicit none
    integer m, ngauss, nmax, ncheck, naxsm
    integer n, ng, ngss, nm, nnmax, &
            n1, n2, i, i1, i2, &
            k1, k2, kk1, kk2, mm1
    real(dp) ppi, pir, pii
    real(dp) a11, a12, a21, a22, aa1, aa2, &
            ar11, ar12, ar21, ar22, &
            ai11, ai12, ai21, ai22, &
            an1, an2, an12, &
            b1r, b1i, b2r, b2i, b3r, b3i, &
            b4r, b4i, b5r, b5i, b6r, b6i, &
            b7r, b7i, b8r, b8i, &
            c1r, c1i, c2r, c2i, c3r, c3i, &
            c4r, c4i, c5r, c5i, c6r, c6i, &
            c7r, c7i, c8r, c8i, &
            d1n1, d1n2, d2n1, d2n2, &
            dd1, dd2, ddri, drii, drri, &
            dsi, e1, e2, e3, f1, f2, &
            factor, &
            gr11, gr12, gr21, gr22, &
            gi11, gi12, gi21, gi22, &
            qdj1, qdjr2, qdji2, qdy1, &
            qj1, qjr2, qji2, qm, qmm, &
            qy1, rri, si, &
            tar11, tar12, tar21, tar22, &
            tai11, tai12, tai21, tai22, &
            tgr11, tgr12, tgr21, tgr22, &
            tgi11, tgi12, tgi21, tgi22, &
            tpir, tpii, tppi, uri, wr, &
            x_to_rsp(1), r_from_x(1), dr_from_x(1), xi, &
            v, dssi

    real(dp)  x(npng2), w(npng2), an(npn1), s(npng2), ss(npng2), &
            r(npng2), dr(npng2), sig(npn2), &
            ddr(npng2), drr(npng2), &
            d1(npng2, npn1), d2(npng2, npn1), &
            dri(npng2), ds(npng2), dss(npng2), rr(npng2), &
            dv1(npn1), dv2(npn1)

    real(dp)  r11(npn1, npn1), r12(npn1, npn1), &
            r21(npn1, npn1), r22(npn1, npn1), &
            i11(npn1, npn1), i12(npn1, npn1), &
            i21(npn1, npn1), i22(npn1, npn1), &
            rg11(npn1, npn1), rg12(npn1, npn1), &
            rg21(npn1, npn1), rg22(npn1, npn1), &
            ig11(npn1, npn1), ig12(npn1, npn1), &
            ig21(npn1, npn1), ig22(npn1, npn1), &
            ann(npn1, npn1), &
            qr(npn2, npn2), qi(npn2, npn2), &
            rgqr(npn2, npn2), rgqi(npn2, npn2), &
            tqr(npn2, npn2), tqi(npn2, npn2), &
            trgqr(npn2, npn2), trgqi(npn2, npn2)
    !c      real(dp) tr1(npn2,npn2),ti1(npn2,npn2)
    !________
    common /tmat99/&
            r11, r12, r21, r22, i11, i12, i21, i22, rg11, rg12, rg21, rg22, &
            ig11, ig12, ig21, ig22          !only between tmatr routines
    !c      common /ct/ tr1,ti1                      !output from tt routine
    common /ctt/ qr, qi, rgqr, rgqi             !input for tt routine
    !________
    mm1 = m
    qm = dble(m)
    qmm = qm*qm
    ng = 2*ngauss
    ngss = ng
    factor = 1d0
    !
    if (ncheck==1) then          !theta=pi/2 is mirror symmetry plane
        ngss = ngauss
        factor = 2d0
    endif
    !
    si = 1d0
    nm = nmax + nmax

    do n = 1, nm
        si = -si
        sig(n) = si              !=(-1)**n
    end do
    !
    ! assigning wigner d-matrices - assuming mirror symmetry
    ! in the \theta=\pi/2 plane:

    do i = 1, ngauss

        i1 = ngauss - i + 1
        i2 = ngauss + i
        !
        call vig (x(i1), nmax, m, dv1, dv2)
        !
        do n = 1, nmax

            dd1 = dv1(n)
            dd2 = dv2(n)
            d1(i1, n) = dd1
            d2(i1, n) = dd2

            if (naxsm==1) then         !gauss abscissas chosen +/- symmetric
                ! using (4.2.4) and (4.2.6) of {ed},
                !           d_{0m}^{(l)}(\pi-\theta) = (-1)^{l+m} d_{0m}^{(l)}(\theta)

                si = sig(n + m)                  !=(-1)**(n+m)
                !                                !exactly what follows from {ed}
                d1(i2, n) = dd1*si
                d2(i2, n) = -dd2*si

            end if
        enddo
        !
        if (naxsm==0) then        !gauss abscissas not chosen +/- symmetric
            !
            call vig (x(i2), nmax, m, dv1, dv2)
            !
            do n = 1, nmax
                dd1 = dv1(n)
                dd2 = dv2(n)
                d1(i2, n) = dd1
                d2(i2, n) = dd2
            enddo

        end if

    end do
    !
    !  assigning r^2(\theta)*weight product:

    do i = 1, ngss
        wr = w(i)*r(i)

        !c           if (dr(i).eq.0.d0) wr=0.d0   !temporarily only

        ds(i) = s(i)*qm*wr       !=dble(m)*w(i)*r^2(\theta)/(|\sin\theta|)
        dss(i) = ss(i)*qmm       !=dble(m)**2/(\sin^2\theta)
        rr(i) = wr
    end do
    !
    do n1 = mm1, nmax
        an1 = an(n1)

        do n2 = mm1, nmax
            an2 = an(n2)
            ar11 = 0d0
            ar12 = 0d0
            ar21 = 0d0
            ar22 = 0d0
            ai11 = 0d0
            ai12 = 0d0
            ai21 = 0d0
            ai22 = 0d0
            gr11 = 0d0
            gr12 = 0d0
            gr21 = 0d0
            gr22 = 0d0
            gi11 = 0d0
            gi12 = 0d0
            gi21 = 0d0
            gi22 = 0d0
            si = sig(n1 + n2)

            integrand%n1 = n1
            integrand%n2 = n2
            integrand%nmax = nmax ! used only for bessel evaluation
            integrand%m = m

            if (ncheck==1.and.si>0d0) then
                !do i = 1, ngss
                !    ar12 = ar12 + w(i)*mi_ar12(x(i))       !~re j^{12}
                !    ai12 = ai12 + w(i)*mi_ai12(x(i))       !~im j^{12}

                !    gr12 = gr12 + w(i)*mi_gr12(x(i))       !~re rg j^{12}
                !    gi12 = gi12 + w(i)*mi_gi12(x(i))        !~im rg j^{12}

                !    ar21 = ar21 + w(i)*mi_ar21(x(i))
                !    ai21 = ai21 + w(i)*mi_ai21(x(i))

                !    gr21 = gr21 + w(i)*mi_gr21(x(i))
                !    gi21 = gi21 + w(i)*mi_gi21(x(i))
                !end do ! end of gauss integration loop
                call integrate(mi_ar12, -1.0_dp, 0.0_dp, ar12)
                call integrate(mi_ai12, -1.0_dp, 0.0_dp, ai12)
                call integrate(mi_gr12, -1.0_dp, 0.0_dp, gr12)
                call integrate(mi_gi12, -1.0_dp, 0.0_dp, gi12)
                call integrate(mi_ar21, -1.0_dp, 0.0_dp, ar21)
                call integrate(mi_ai21, -1.0_dp, 0.0_dp, ai21)
                call integrate(mi_gr21, -1.0_dp, 0.0_dp, gr21)
                call integrate(mi_gi21, -1.0_dp, 0.0_dp, gi21)

            else
                !do i = 1, ngss
                !    ar11 = ar11 + w(i)*mi_ar11(x(i))
                !    ai11 = ai11 + w(i)*mi_ai11(x(i))
                !    gr11 = gr11 + w(i)*mi_gr11(x(i))
                !    gi11 = gi11 + w(i)*mi_gi11(x(i))

                !    ar22 = ar22 + w(i)*mi_ar22(x(i))
                !    ai22 = ai22 + w(i)*mi_ai22(x(i))

                !    gr22 = gr22 + w(i)*mi_gr22(x(i))
                !    gi22 = gi22 + w(i)*mi_gi22(x(i))
                !end do ! end of gauss integration loop
                call integrate(mi_ar11, -1.0_dp, 0.0_dp, ar11)
                call integrate(mi_ai11, -1.0_dp, 0.0_dp, ai11)
                call integrate(mi_gr11, -1.0_dp, 0.0_dp, gr11)
                call integrate(mi_gi11, -1.0_dp, 0.0_dp, gi11)
                call integrate(mi_ar22, -1.0_dp, 0.0_dp, ar22)
                call integrate(mi_ai22, -1.0_dp, 0.0_dp, ai22)
                call integrate(mi_gr22, -1.0_dp, 0.0_dp, gr22)
                call integrate(mi_gi22, -1.0_dp, 0.0_dp, gi22)
            end if

            an12 = ann(n1, n2)*factor

            r11(n1, n2) = ar11*an12
            r12(n1, n2) = ar12*an12
            r21(n1, n2) = ar21*an12
            r22(n1, n2) = ar22*an12
            i11(n1, n2) = ai11*an12
            i12(n1, n2) = ai12*an12
            i21(n1, n2) = ai21*an12
            i22(n1, n2) = ai22*an12

            rg11(n1, n2) = gr11*an12
            rg12(n1, n2) = gr12*an12
            rg21(n1, n2) = gr21*an12
            rg22(n1, n2) = gr22*an12
            ig11(n1, n2) = gi11*an12
            ig12(n1, n2) = gi12*an12
            ig21(n1, n2) = gi21*an12
            ig22(n1, n2) = gi22*an12
        end do
    end do

    tpir = pir                 !re [1/k_{in}^2]
    tpii = pii                 !im [1/k_{in}^2]
    tppi = ppi                 !1/k_{out}^2

    nm = nmax - mm1 + 1
    do n1 = mm1, nmax
        k1 = n1 - mm1 + 1
        kk1 = k1 + nm

        do n2 = mm1, nmax
            k2 = n2 - mm1 + 1
            kk2 = k2 + nm

            tar11 = -r11(n1, n2)
            tai11 = -i11(n1, n2)
            tgr11 = -rg11(n1, n2)
            tgi11 = -ig11(n1, n2)

            tar12 = i12(n1, n2)
            tai12 = -r12(n1, n2)
            tgr12 = ig12(n1, n2)
            tgi12 = -rg12(n1, n2)

            tar21 = -i21(n1, n2)
            tai21 = r21(n1, n2)
            tgr21 = -ig21(n1, n2)
            tgi21 = rg21(n1, n2)

            tar22 = -r22(n1, n2)
            tai22 = -i22(n1, n2)
            tgr22 = -rg22(n1, n2)
            tgi22 = -ig22(n1, n2)

            tqr(k1, k2) = tpir*tar21 - tpii*tai21 + tppi*tar12
            tqi(k1, k2) = tpir*tai21 + tpii*tar21 + tppi*tai12
            trgqr(k1, k2) = tpir*tgr21 - tpii*tgi21 + tppi*tgr12
            trgqi(k1, k2) = tpir*tgi21 + tpii*tgr21 + tppi*tgi12

            tqr(k1, kk2) = tpir*tar11 - tpii*tai11 + tppi*tar22
            tqi(k1, kk2) = tpir*tai11 + tpii*tar11 + tppi*tai22
            trgqr(k1, kk2) = tpir*tgr11 - tpii*tgi11 + tppi*tgr22
            trgqi(k1, kk2) = tpir*tgi11 + tpii*tgr11 + tppi*tgi22

            tqr(kk1, k2) = tpir*tar22 - tpii*tai22 + tppi*tar11
            tqi(kk1, k2) = tpir*tai22 + tpii*tar22 + tppi*tai11
            trgqr(kk1, k2) = tpir*tgr22 - tpii*tgi22 + tppi*tgr11
            trgqi(kk1, k2) = tpir*tgi22 + tpii*tgr22 + tppi*tgi11

            tqr(kk1, kk2) = tpir*tar12 - tpii*tai12 + tppi*tar21
            tqi(kk1, kk2) = tpir*tai12 + tpii*tar12 + tppi*tai21
            trgqr(kk1, kk2) = tpir*tgr12 - tpii*tgi12 + tppi*tgr21
            trgqi(kk1, kk2) = tpir*tgi12 + tpii*tgr12 + tppi*tgi21
        end do
    end do

    nnmax = 2*nm
    do n1 = 1, nnmax
        do n2 = 1, nnmax
            qr(n1, n2) = tqr(n1, n2)
            qi(n1, n2) = tqi(n1, n2)
            rgqr(n1, n2) = trgqr(n1, n2)
            rgqi(n1, n2) = trgqi(n1, n2)
        end do
    end do

    call tt(nm)
    !
    return
end

!**********************************************************************

subroutine tmatr_leru(m, ngauss, x, w, an, ann, s, ss, ppi, pir, pii,  &
        r, dr, ddr, drr, dri, nmax, ncheck, naxsm)
    !--------/---------/---------/---------/---------/---------/---------/--
    ! >>> ngauss,x,w,an,ann,s,ss,ppi,pir,pii,r,dr,ddr,drr,dri,nmax,ncheck
    ! <<< common blocks /tmat99/, /ct/ (for main),  and /ctt/ (for tt)
    !=====================
    ! complex(dp), parameter:: czero=cmplx_dp(0.0_dp,0.0_dp), ci=cmplx_dp(0.0_dp, 1.0_dp)
    ! are supplied here by the "constants" routine
    !
    !  determines the t-matrix of an axially symmetric scatterer
    !                           for m > 0
    !
    !  m      - azimuthal number
    !  ngauss - the number of gif division points
    !  x=\cos\theta  - gif division points
    !  w - gif weights
    !  an(n)=n*(n+1)
    !  ann(l_1,l_2)=\sqrt{\fr{(2 l_1+1)}{l_1(l_1+1)} }
    !                       \sqrt{\fr{(2 l_2+1)}{l_2(l_2+1)} }/2
    !  nmax - angular momentum cutoff
    !  ncheck  -  .eq.0  then  ngss=2*ngauss, factor=1d0
    !             .eq.1  then  ngss = ngauss, factor=2d0
    !  naxsm   -  .eq.0 : gauss abscissas do not have +/- theta symmetry
    !             .eq.1 : gauss abscissas have +/- theta symmetry
    !  ncheck - specifies whether ng=2*ngauss or otherwise
    !  MRR - the real part of the rel.  refractive index
    !  MRI - the imaginary of the rel. the refractive index
    !  p=dacos(-1d0)
    !  pi=p*2d0/lam      !wave vector k=2\pi/\ld in the host medium
    !  ppi=pi*pi         !k^2
    !  pir=ppi*mrr       !k^2*mrr
    !  pii=ppi*mri       !k^2*mri
    !  r=r^2(\theta)                       for axially symmetric particles
    !  dr=[dr(\theta)/(d\theta)]/r(\theta) for axially symmetric particles
    !  ddr=\lambda/[2*\pi*r(\theta)]=1/(k_out*r)
    !  drr=(mrr/(mrr**2+mri**2))*(\lambda/[2*\pi*r(\theta)])
    !                  = re 1/(k_in*r)
    !  dri=-(mri/(mrr**2+mri**2))*(\lambda/[2*\pi*r(\theta)])
    !                  = im 1/(k_in*r)
    !  Refractive index outside is assumed to real, whereas inside
    !  a scatterer, refractive index is allowed to be complex in general.
    !  Consequently, the bessel function j_l(k_in*r) will in general
    !  be complex. The routine below performs Waterman surface integral
    !  separately for the real and imaginary parts of the integrand.
    !--------/---------/---------/---------/---------/---------/---------/--
    use libcylinder
    !include 'ampld.par.f'
    implicit none
    integer, intent(in) :: m, ngauss, nmax, ncheck, naxsm
    integer i, i1, i2, k1, k2, kk1, kk2, &
            mm1, n, n1, n2, ng, ngss, nm, nnmax
! npn1, npn2,npng2 have been declared as public integers

    real(dp), intent(in) :: ppi, pir, pii
    real(dp) a11, a12, a21, a22, aa1, aa2, &
            ar11, ar12, ar21, ar22, &
            ai11, ai12, ai21, ai22, &
            an1, an2, an12, &
            b1r, b1i, b2r, b2i, b3r, b3i, &
            b4r, b4i, b5r, b5i, b6r, b6i, &
            b7r, b7i, b8r, b8i, &
            c1r, c1i, c2r, c2i, c3r, c3i, &
            c4r, c4i, c5r, c5i, c6r, c6i, &
            c7r, c7i, c8r, c8i, &
            d1n1, d1n2, d2n1, d2n2, &
            dd1, dd2, ddri, drii, drri, &
            dsi, e1, e2, e3, f1, f2, &
            factor, &
            gr11, gr12, gr21, gr22, &
            gi11, gi12, gi21, gi22, &
            qdj1, qdjr2, qdji2, qdy1, &
            qj1, qjr2, qji2, qm, qmm, &
            qy1, rri, si, &
            tar11, tar12, tar21, tar22, &
            tai11, tai12, tai21, tai22, &
            tgr11, tgr12, tgr21, tgr22, &
            tgi11, tgi12, tgi21, tgi22, &
            tpir, tpii, tppi, uri, wr

    real(dp), intent(in) :: an(npn1), ann(npn1, npn1), s(npng2), &
            ss(npng2), r(npng2), dr(npng2), ddr(npng2), drr(npng2), &
            dri(npng2), x(npng2), w(npng2)

    real(dp)  sig(npn2), &
            d1(npng2, npn1), d2(npng2, npn1), &
            ds(npng2), dss(npng2), rr(npng2), &
            dv1(npn1), dv2(npn1)

    real(dp)  r11(npn1, npn1), r12(npn1, npn1), &
            r21(npn1, npn1), r22(npn1, npn1), &
            i11(npn1, npn1), i12(npn1, npn1), &
            i21(npn1, npn1), i22(npn1, npn1), &
            rg11(npn1, npn1), rg12(npn1, npn1), &
            rg21(npn1, npn1), rg22(npn1, npn1), &
            ig11(npn1, npn1), ig12(npn1, npn1), &
            ig21(npn1, npn1), ig22(npn1, npn1), &
            qr(npn2, npn2), qi(npn2, npn2), &
            rgqr(npn2, npn2), rgqi(npn2, npn2), &
            tqr(npn2, npn2), tqi(npn2, npn2), &
            trgqr(npn2, npn2), trgqi(npn2, npn2)
    !c      real(dp) tr1(npn2,npn2),ti1(npn2,npn2)
    !Ru_patch:
    real(dp) :: drd, qs, xkr, a(nmax)
    complex(dp) :: znf, zxipsi, zxidpsi, zdxipsi, zdxidpsi, &
            zk1, zk2, zl1, zl2, zl3, zl4, zk1i, zk2i, zl1i, zl2i, zl3i, zl4i
    complex(dp) zq11(npn1, npn1), zq12(npn1, npn1), zq21(npn1, npn1), &
            zq22(npn1, npn1)
    !________
    common /tmat99/&
            r11, r12, r21, r22, i11, i12, i21, i22, rg11, rg12, rg21, rg22, &
            ig11, ig12, ig21, ig22          !only between tmatr routines
    !c      common /ct/ tr1,ti1                      !output from tt routine
    common /ctt/ qr, qi, rgqr, rgqi             !input for tt routine
    common /ztt/ zq11, zq12, zq21, zq22
    !________
    mm1 = m                !in tmatr0 routine mm1=1 for m=0
    qm = dble(m)
    qmm = qm*qm
    ng = 2*ngauss
    ngss = ng
    factor = 1d0
    !
    if (ncheck==1) then          !theta=pi/2 is mirror symmetry plane
        ngss = ngauss
        factor = 2d0
    endif
    !
    si = 1d0
    nm = nmax + nmax

    do n = 1, nm
        si = -si
        sig(n) = si              !=(-1)**n
    end do

! Assigning Wigner d-matrices for in total 2*ngauss points

    do i = 1, ngauss

        i1 = ngauss - i + 1
        i2 = ngauss + i
        !
        call vig (x(i1), nmax, m, dv1, dv2)
        !
        do n = 1, nmax

            dd1 = dv1(n)
            dd2 = dv2(n)
            d1(i1, n) = dd1
            d2(i1, n) = dd2

            if (naxsm==1) then         !gauss abscissas chosen +/- symmetric
                ! using (4.2.4) and (4.2.6) of {ed},
                !           d_{0m}^{(l)}(\pi-\theta) = (-1)^{l+m} d_{0m}^{(l)}(\theta)

                si = sig(n + m)                  !=(-1)**(n+m)
                                                 !exactly what follows from {ed}
                d1(i2, n) = dd1*si
                d2(i2, n) = -dd2*si

            end if
        enddo
        !
        if (naxsm==0) then                !gauss abscissas not chosen +/- symmetric
            !
            call vig (x(i2), nmax, m, dv1, dv2)
            !
            do n = 1, nmax
                dd1 = dv1(n)
                dd2 = dv2(n)
                d1(i2, n) = dd1
                d2(i2, n) = dd2
            enddo

        end if

    end do   !over to ngauss

!  Assigning r^2(\theta)*weight product and Bessel function independent factors:
!  The arrays ss, s have been determined in CONST as 1/(\sin^2\theta) and
!  1/(|\sin\theta|), respectively:

    do i = 1, ngss
        wr = w(i)*r(i)           !r^2(\theta)*weight
        !c           if (dr(i).eq.0.d0) wr=0.d0   !temporarily only
        ds(i) = s(i)*qm*wr       !=dble(m)*w(i)*r^2(\theta)/(|\sin\theta|)
        dss(i) = ss(i)*qmm       !=dble(m)**2/(\sin^2\theta)
        rr(i) = wr               !r^2(\theta)*weight
    end do

!  Before the n1, n2 (or l,l') loop:
!  In Ru et al formulas:
!  s=n_2/n_1$ is the refractive index constrast and $x(\theta)=k_1 r(\theta)$ is the
!  polar angle dependent size parameter of an axially symmetric scatterer.

    znf = cmplx_dp(pir, pii)/ppi       !in general complex ref. index contrast n_2/n_1;
                                       !the parameter 's' in Ru et al formulas

    do n=mm1,nmax
        A(n) = sqrt(dble(2*n + 1)/(2*n*(n + 1)))
    enddo

    do n1 = mm1, nmax          !mm1 = m
        an1 = an(n1)

        do n2 = mm1, nmax      !mm1 = m
            an2 = an(n2)
!Init:
            ar11 = 0d0
            ar12 = 0d0
            ar21 = 0d0
            ar22 = 0d0
            ai11 = 0d0
            ai12 = 0d0
            ai21 = 0d0
            ai22 = 0d0
            gr11 = 0d0
            gr12 = 0d0
            gr21 = 0d0
            gr22 = 0d0
            gi11 = 0d0
            gi12 = 0d0
            gi21 = 0d0
            gi22 = 0d0
            si = sig(n1 + n2)

! Init of Ru et al integrals:
            zk1 = czero
            zk2 = czero
            zl1 = czero
            zl2 = czero
            zl3 = czero
            zl4 = czero

! Gauss integration loop:
! Traditional integration transform:
! \int_0^\pi X* \sin\theta d\theta --> \int_{-1}^1 X* d(\cos\theta)
! and generates weights for the independent variable "\cos\theta"

            do i = 1, ngss       !=ngauss   if ncheck.eq.1
                                 !=2*ngauss if ncheck.eq.0
                d1n1 = d1(i, n1)
                d2n1 = d2(i, n1)
                d1n2 = d1(i, n2)
                d2n2 = d2(i, n2)

!aij store important products d_n.d_{n'}, d_n.\tau_{n'}, \tau_n.d_{n'},...
!necessary for forming the integrands of Q-matrices. In Le Ru patch
!a22 is no longer needed.
!The first index of aij is 1 (2) if the factor with n1-index is d1 (d2)
!the second index of aij is 1 (2) if the factor with n2-index is d1 (d2)

                a11 = d1n1*d1n2
                a12 = d1n1*d2n2
                a21 = d2n1*d1n2
                a22 = d2n1*d2n2
                aa1 = a12 + a21            != d1n1*d2n2+d2n1*d1n2
                aa2 = a11*dss(i) + a22     !=(d1n1*d1n2)*dble(m)**2/(\sin^2\theta)
                                           ! +d2n1*d2n2

! Spherical Bessel functions: outside (1); inside (2)
! Since refractive index is allowed to be complex in general,
! the Bessel function j_l(k_in*r) is complex. The code below
! performs a separation of the complex integrand in Waterman's
! surface integral into its respective real and imaginary
! parts.

! Bessel functions of the real exterior argument:

                qj1 = cbess%j(i, n1)
                qy1 = cbess%y(i, n1)
                qdj1 = cbess%dj(i, n1)
                qdy1 = cbess%dy(i, n1)

! Bessel functions of the, in principle complex, interior argument:
                qjr2 = cbess%jr(i, n2)
                qji2 = cbess%ji(i, n2)
                qdjr2 = cbess%djr(i, n2)
                qdji2 = cbess%dji(i, n2)

!--------/---------/---------/---------/---------/---------/---------/--
!Assigning integrads by eqs 2-7 of Ru2012
!Bessel function part:
!Amend for the Riccati-Bessel functions of Ru et al.
!The Riccati-Bessel functions \xi_l(z)=zh_l(z) (\psi_l(z)=zj_l(z)) are the functions of
!the exterior (interior) argument there.
! \xi or d\xi carry always the subscript n (=n1 here) and
! \psi or d\psi carry always the subscript k (=n2 here)
!The latter allows one to initiate their complex valued products.
!  PI=P*2D0/LAM                 !wave vector k=2\pi/\ld in the host medium
!  PPI=PI*PI                    !k^2
!  PIR=PPI*MRR                  !k^2*mrr
!  PII=PPI*MRI                  !k^2*mri
!  R(I):=r^2(\theta)

        xkr = sqrt(ppi*r(i))        !k_out*r(\theta)

        zxipsi = cmplx_dp(pir, pii)*cmplx_dp(qj1, qy1)*cmplx_dp(qjr2, qji2)*r(i)

        zxidpsi = xkr*cmplx_dp(qj1, qy1)*cmplx_dp(qjr2, qji2) + &
                cmplx_dp(pir, pii)*cmplx_dp(qj1, qy1)*cmplx_dp(qdjr2, qdji2)*r(i)

        zdxipsi = znf*xkr*cmplx_dp(qj1, qy1)*cmplx_dp(qjr2, qji2) + &
                cmplx_dp(pir, pii)*cmplx_dp(qdj1, qdy1)*cmplx_dp(qjr2, qji2)*r(i)

        zdxidpsi = cmplx_dp(qj1, qy1)*cmplx_dp(qjr2, qji2) + &
                znf*xkr*cmplx_dp(qj1, qy1)*cmplx_dp(qdjr2, qdji2) + &
                xkr*cmplx_dp(qdj1, qdy1)*cmplx_dp(qjr2, qji2) + &
                cmplx_dp(pir, pii)*cmplx_dp(qdj1, qdy1)*cmplx_dp(qdjr2, qdji2)*r(i)
!--------/---------/---------/---------/---------/---------/---------/--
! Gauss integration weigths and variables are in the variable '\cos\theta'.
! At the present code of Ref. \ct{MTL}:
! r(i) contains $r^2(\theta)$, whereas we need $r(\theta)$
! dr(i) contains $r'(\theta)/r(\theta)$, whereas we need $r'(\theta)$
! ds(i) = s(i)*qm*wr       !=dble(m)*w(i)*r^2(\theta)/(|\sin\theta|)

                qs = dsqrt(1.d0 - x(i)**2)          !\sin\theta             
                drd = xkr*dr(i)                     !x_\theta=k_out*r'(\theta) of Ru et al

! Amend for the Riccati-Bessel functions of Ru et al
! on using that they only contain combinations of
! \xi or d\xi with the subscript n (=n1 here) and
! \psi or d\psi with the subscript k (=n2 here)
! Assigning integrads by eqs 2-7 of Ru2012:

                if (ncheck==1) then         !theta=pi/2 is scatterer mirror symmetry plane
                                            !eqs 18-19 of Ru2012
                    if (mod(n1 + n2, 2)==0) then     !n1+n2 even

                        zk1i = czero
                        zk2i = czero

                        zl1i = a21*drd*zxipsi
                        zl2i = a12*drd*zxipsi
                        zl3i = a21*drd*zdxidpsi &
                                - dble(n1*(n1 + 1))*a11*zxidpsi
                        zl4i = a12*znf*drd*zdxidpsi &
                                - dble(n2*(n2 + 1))*a11*zdxipsi

                    else                    !n1+n2 odd
                        zk1i = qm*a11*drd*zxidpsi/qs
                        zk2i = qm*a11*drd*zdxipsi/qs
! qs=dsqrt(1.d0-x(i)**2) is \sin\theta;
! By dividing by qs, the integrand will, as before, be integrated with the measure
! d(\cos\theta) for which the GIF weights are supplied

                        zl1i = czero
                        zl2i = czero
                        zl3i = czero
                        zl4i = czero

                    end if

                else if (ncheck==0) then        !theta=pi/2 is not a scatterer mirror symmetry plane

                    zk1i = qm*a11*drd*zxidpsi/qs
                    zk2i = qm*a11*drd*zdxipsi/qs
                    zl1i = a21*drd*zxipsi
                    zl2i = a12*drd*zxipsi
                    zl3i = a21*drd*zdxidpsi &
                            - dble(n1*(n1 + 1))*a11*zxidpsi
                    zl4i = a12*znf*drd*zdxidpsi &
                            - dble(n2*(n2 + 1))*a11*zdxipsi
    !--------/---------/---------/---------/---------/---------/---------/--

                end if

    ! Integrate Ru et al integrals:
                zk1 = zk1 + zk1i*w(i)
                zk2 = zk2 + zk2i*w(i)
                zl1 = zl1 + zl1i*w(i)
                zl2 = zl2 + zl2i*w(i)
                zl3 = zl3 + zl3i*w(i)
                zl4 = zl4 + zl4i*w(i)

    ! todo: a switch can be implemented here to run the code
    ! either by following Ru et al integration or
    ! following the original MTL path

        ! re and im of j_{n2}(k_{in}r) j_{n1}(k_{out}r):
                c1r = qjr2*qj1
                c1i = qji2*qj1

        ! re and im of j_{n2}(k_{in}r) h_{n1}(k_{out}r):
                b1r = c1r - qji2*qy1
                b1i = c1i + qjr2*qy1

        ! re and im of j_{n2}(k_{in}r) j_{n1}'(k_{out}r):
                c2r = qjr2*qdj1
                c2i = qji2*qdj1

        ! re and im of j_{n2}(k_{in}r) h_{n1}'(k_{out}r):
                b2r = c2r - qji2*qdy1
                b2i = c2i + qjr2*qdy1

                ddri = ddr(i)               !1/(k_{out}r)

        !re and im of [1/(k_{out}r)]*j_{n2}(k_{in}r) j_{n1}(k_{out}r)
                c3r = ddri*c1r
                c3i = ddri*c1i

        ! re and im of [1/(k_{out}r)]*j_{n2}(k_{in}r) h_{n1}(k_{out}r):
                b3r = ddri*b1r
                b3i = ddri*b1i

        ! re and im of j_{n2}'(k_{in}r) j_{n1}(k_{out}r):
                c4r = qdjr2*qj1
                c4i = qdji2*qj1

        ! re and im of j_{n2}'(k_{in}r) h_{n1}(k_{out}r):
                b4r = c4r - qdji2*qy1
                b4i = c4i + qdjr2*qy1

                drri = drr(i)               !re[1/(k_{in}r)]
                drii = dri(i)               !im[1/(k_{in}r)]
        ! re and im of [1/(k_{in}r)] j_{n2}(k_{in}r) j_{n1}(k_{out}r):

                c5r = c1r*drri - c1i*drii
                c5i = c1i*drri + c1r*drii
        ! re and im of [1/(k_{in}r)] j_{n2}(k_{in}r) h_{n1}(k_{out}r):

                b5r = b1r*drri - b1i*drii
                b5i = b1i*drri + b1r*drii
        ! re and im of j_{n2}'(k_{in}r) j_{n1}'(k_{out}r):

                c6r = qdjr2*qdj1
                c6i = qdji2*qdj1
        ! re and im of j_{n2}'(k_{in}r) h_{n1}'(k_{out}r):

                b6r = c6r - qdji2*qdy1
                b6i = c6i + qdjr2*qdy1
        ! re and im of [1/(k_{out}r)] j_{n2}'(k_{in}r) j_{n1}(k_{out}r):

                c7r = c4r*ddri
                c7i = c4i*ddri
        ! re and im of [1/(k_{out}r)] j_{n2}'(k_{in}r) h_{n1}(k_{out}r):

                b7r = b4r*ddri
                b7i = b4i*ddri
        ! re and im of [1/(k_{in}r)] j_{n2}(k_{in}r) j_{n1}'(k_{out}r):

                c8r = c2r*drri - c2i*drii
                c8i = c2i*drri + c2r*drii
        ! re and im of [1/(k_{in}r)] j_{n2}(k_{in}r) h_{n1}'(k_{out}r):

                b8r = b2r*drri - b2i*drii
                b8i = b2i*drri + b2r*drii

! %%%%%%%%%  forming integrands of j-matrices (j^{11}=j^{22}=0 for m=0):

                uri = dr(i)
                dsi = ds(i)
                rri = rr(i)

                if (ncheck==1.and.si>0d0) then
                    ! w(i)*r^2(\theta)*[(d1n1*d1n2)*dble(m)**2/(\sin^2\theta)+d2n1*d2n2]:
                    f1 = rri*aa2            !prefactor containing r^2(\theta)<->hat{r} part
                    ! n1*(n1+1)*w(i)*r(\theta)*[dr/(d\theta)]*d1n1*d2n2:
                    f2 = rri*uri*an1*a12     !prefactor containing r(\theta)*[dr/(d\theta)]
                    !                                          !hat{theta} part

                    ar12 = ar12 + f1*b2r + f2*b3r        !~re j^{12}
                    ai12 = ai12 + f1*b2i + f2*b3i        !~im j^{12}

                    gr12 = gr12 + f1*c2r + f2*c3r        !~re rg j^{12}
                    gi12 = gi12 + f1*c2i + f2*c3i        !~im rg j^{12}
                    ! n2*(n2+1)*w(i)*r(\theta)*[dr/(d\theta)]*d2n1*d1n2:
                    f2 = rri*uri*an2*a21     !prefactor containing r(\theta)*[dr/(d\theta)]
                    !                                          !hat{theta} part

                    ar21 = ar21 + f1*b4r + f2*b5r
                    ai21 = ai21 + f1*b4i + f2*b5i

                    gr21 = gr21 + f1*c4r + f2*c5r
                    gi21 = gi21 + f1*c4i + f2*c5i

                    cycle         !causes the loop to skip the remainder of its body
                end if
        ! [dble(m)*w(i)*r^2(i)/(|\sin\theta|)]*(d1n1*d2n2+d2n1*d1n2):
                e1 = dsi*aa1

                ar11 = ar11 + e1*b1r
                ai11 = ai11 + e1*b1i
                gr11 = gr11 + e1*c1r
                gi11 = gi11 + e1*c1i

                e2 = dsi*uri*a11
                e3 = e2*an2
                e2 = e2*an1

                ar22 = ar22 + e1*b6r + e2*b7r + e3*b8r
                ai22 = ai22 + e1*b6i + e2*b7i + e3*b8i

                gr22 = gr22 + e1*c6r + e2*c7r + e3*c8r
                gi22 = gi22 + e1*c6i + e2*c7i + e3*c8i

            end do  ! end of Gauss integration

        ! Taking into account the "factor":
        ! factor=1 or 2 depending on ncheck==0 or ncheck==1

            if (abs(factor-2d0)<=1.d-10) then
                zk1 = zk1*factor
                zk2 = zk2*factor
                zl1 = zl1*factor
                zl2 = zl2*factor
                zl3 = zl3*factor
                zl4 = zl4*factor
            end if

!After performing surface integrals, forming the Q-matrices
!by eqs 12-17 of Ru2012

            if (n1/=n2) then
                zq11(n1, n2) = ci*A(n1)*A(n2)*(znf**2 - cone)*&
                        (n1*(n1 + 1)*zl2 - n2*(n2 + 1)*zl1)/((n1*(n1 + 1) - n2*(n2 + 1))*znf)
                zq22(n1, n2) = ci*A(n1)*A(n2)*(znf**2 - cone)*&
                        (zl3 + znf*n1*(n1 + 1)*(zl2 - zl1)/(n1*(n1 + 1) - n2*(n2 + 1)))/znf
            else                !diagonal elements
                zq11(n1, n2) = -ci*A(n1)*A(n2)*(-znf*zl1 + zl3 + (zl2 - zl4)/znf)
                zq22(n1, n2) = -ci*A(n1)*A(n2)*(-zl1 + zl2 - zl4 + zl3/znf)
            end if

            zq12(n1, n2) =  A(n1)*A(n2)*(znf**2 - cone)*zk1/znf
            zq21(n1, n2) = -A(n1)*A(n2)*(znf**2 - cone)*zk2/znf

! todo: a switch can be implemented here to run the code
! either by following Ru et al integration or
! following the original MTL path
! Original part continues here

            an12 = ann(n1, n2)*factor

            r11(n1, n2) = ar11*an12
            r12(n1, n2) = ar12*an12
            r21(n1, n2) = ar21*an12
            r22(n1, n2) = ar22*an12
            i11(n1, n2) = ai11*an12
            i12(n1, n2) = ai12*an12
            i21(n1, n2) = ai21*an12
            i22(n1, n2) = ai22*an12

            rg11(n1, n2) = gr11*an12
            rg12(n1, n2) = gr12*an12
            rg21(n1, n2) = gr21*an12
            rg22(n1, n2) = gr22*an12
            ig11(n1, n2) = gi11*an12
            ig12(n1, n2) = gi12*an12
            ig21(n1, n2) = gi21*an12
            ig22(n1, n2) = gi22*an12
        end do              !over n2
    end do                  !over n1

    tpir = pir                 !re [1/k_{in}^2]
    tpii = pii                 !im [1/k_{in}^2]
    tppi = ppi                 !1/k_{out}^2

    nm = nmax - mm1 + 1         !mm1 = m
    do n1 = mm1, nmax
        k1 = n1 - mm1 + 1       !k1 set to run from k1=1 for any mm1=m
        kk1 = k1 + nm           !kk1 set to run from kk1=nmax-mm1+2 for any mm1=m

        do n2 = mm1, nmax
            k2 = n2 - mm1 + 1
            kk2 = k2 + nm

            tar11 = -r11(n1, n2)
            tai11 = -i11(n1, n2)
            tgr11 = -rg11(n1, n2)
            tgi11 = -ig11(n1, n2)

            tar12 = i12(n1, n2)
            tai12 = -r12(n1, n2)
            tgr12 = ig12(n1, n2)
            tgi12 = -rg12(n1, n2)

            tar21 = -i21(n1, n2)
            tai21 = r21(n1, n2)
            tgr21 = -ig21(n1, n2)
            tgi21 = rg21(n1, n2)

            tar22 = -r22(n1, n2)
            tai22 = -i22(n1, n2)
            tgr22 = -rg22(n1, n2)
            tgi22 = -ig22(n1, n2)

            tqr(k1, k2) = tpir*tar21 - tpii*tai21 + tppi*tar12
            tqi(k1, k2) = tpir*tai21 + tpii*tar21 + tppi*tai12
            trgqr(k1, k2) = tpir*tgr21 - tpii*tgi21 + tppi*tgr12
            trgqi(k1, k2) = tpir*tgi21 + tpii*tgr21 + tppi*tgi12

            tqr(k1, kk2) = tpir*tar11 - tpii*tai11 + tppi*tar22
            tqi(k1, kk2) = tpir*tai11 + tpii*tar11 + tppi*tai22
            trgqr(k1, kk2) = tpir*tgr11 - tpii*tgi11 + tppi*tgr22
            trgqi(k1, kk2) = tpir*tgi11 + tpii*tgr11 + tppi*tgi22

            tqr(kk1, k2) = tpir*tar22 - tpii*tai22 + tppi*tar11
            tqi(kk1, k2) = tpir*tai22 + tpii*tar22 + tppi*tai11
            trgqr(kk1, k2) = tpir*tgr22 - tpii*tgi22 + tppi*tgr11
            trgqi(kk1, k2) = tpir*tgi22 + tpii*tgr22 + tppi*tgi11

            tqr(kk1, kk2) = tpir*tar12 - tpii*tai12 + tppi*tar21
            tqi(kk1, kk2) = tpir*tai12 + tpii*tar12 + tppi*tai21
            trgqr(kk1, kk2) = tpir*tgr12 - tpii*tgi12 + tppi*tgr21
            trgqi(kk1, kk2) = tpir*tgi12 + tpii*tgr12 + tppi*tgi21
        end do
    end do

    nnmax = 2*nm
    do n1 = 1, nnmax
        do n2 = 1, nnmax
            qr(n1, n2) = tqr(n1, n2)
            qi(n1, n2) = tqi(n1, n2)
            rgqr(n1, n2) = trgqr(n1, n2)
            rgqi(n1, n2) = trgqi(n1, n2)
        end do
    end do

    call tt(nm)
    !
    return
end

!**********************************************************************

subroutine tt(nmax)
!=================
!  nmax=nmax-m+1 here, where nmax is the angular momentum cutoff in main
!  ncheck
!
!   calculation of the matrix    t = - rg(q)*(q**(-1))
!
!   input in commons /ctt/ and  /ztt/
!   output in common /ct/
!
!--------/---------/---------/---------/---------/---------/---------/--
    use libcylinder
    implicit none
! npn1, npn2 have been declared as public integers
    integer nnmax, nmax, nm, i, info, j, k
    integer, allocatable :: ipiv(:)     !an allocatable array is assigned memory storage
                                        !during execution
    real(dp) ai, ar, ari, arr, ti, tr
    real(dp)  qr(npn2, npn2), qi(npn2, npn2), emach, &
            rgqr(npn2, npn2), rgqi(npn2, npn2)
    real(dp) tr1(npn2, npn2), ti1(npn2, npn2)    !
    complex(dp) zq11(npn1, npn1), zq12(npn1, npn1), zq21(npn1, npn1), &
            zq22(npn1, npn1)
    complex(dp), allocatable :: zq(:, :), zx(:), zw(:)
    !
    common /ct/ tr1, ti1
    common /ctt/ qr, qi, rgqr, rgqi
    common /ztt/ zq11, zq12, zq21, zq22
    !
    data emach/0.d0/                 !/1.d-17/
    !
    nnmax = 2*nmax
    allocate(zq(1:nnmax, 1:nnmax), ipiv(1:nnmax), &
            zx(1:nnmax), zw(1:nnmax))       !the actual bounds for the allocatable arrays
                                            !determined
    k=1
!>>>
    if (k==1) then
        do i = 1, nmax
            do j = 1, nmax
                zq(i, j) = zq11(i,j)
                zq(i, nmax+j) = zq12(i,j)
                zq(nmax+i, j) = zq21(i,j)
                zq(nmax+i, nmax+j) = zq22(i,j)
            end do
        end do

    else

        do i = 1, nnmax
            do j = 1, nnmax
                zq(i, j) = cmplx_dp(qr(i, j), qi(i, j))
            enddo
        enddo

    end if
!<<<
!
!Select Gaussian elimination:
!mpar%ichoice == 1 .... lapack routines for the matrix inversion
!mpar%ichoice == 2 ...  Pendry's zger and zsur are used

    if (mpar%ichoice==2) then       ! gaussian elimination; nag library not used

        call zger(zq, ipiv, emach)  !gauss elimination of zq to a
                                    !lower diagonal matrix
        do i = 1, nnmax
            do k = 1, nnmax    !initialization of the right-hand side zb
                               !(a row vector) of the matrix equation zx*zq=zb
                zx(k) = cmplx_dp(rgqr(i, k), rgqi(i, k))
            enddo

!solving zx*zq=zb by backsubstition (zx overwritten on exit)
            call zsur(zq, ipiv, zx, emach)

            do k = 1, nnmax

                ! assign t-matrix elements = - rg(q)*(q**(-1))
                tr1(i, k) = -dble(zx(k))
                ti1(i, k) = -aimag(zx(k))
            enddo
        end do

        deallocate(zq, ipiv, zx, zw)    !release arrays from memory

    else                 !matrix inversion from lapack

!        do i = 1, nnmax
!            do j = 1, nnmax
!                zq(i, j) = cmplx_dp(qr(i, j), qi(i, j))
!            enddo
!        enddo

        info = 0
        !     call zgetrf_wrap(zq, ipiv)
        call zgetrf(nnmax, nnmax, zq, nnmax, ipiv, info)
        if (info/=0) write (6, 1100) info
        call zgetri(nnmax, zq, nnmax, ipiv, zw, nnmax, info)
        if (info/=0) write (6, 1100) info

        1100      format ('warning:  info=', i2)

! calculate t-matrix = - rg(q)*(q**(-1))
        do i = 1, nnmax
            do j = 1, nnmax
                tr = 0d0
                ti = 0d0
                do k = 1, nnmax
                    arr = rgqr(i, k)
                    ari = rgqi(i, k)
                    ar = zq(k, j)
                    ai = aimag(zq(k, j))
                    tr = tr - arr*ar + ari*ai
                    ti = ti - arr*ai - ari*ar
                enddo
                tr1(i, j) = tr
                ti1(i, j) = ti
            enddo
        enddo
    end if

    return
end