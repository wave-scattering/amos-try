!DMQMULTEM, VERSION 2.  MULTEM 2: A NEW VERSION OF THE PROGRAM FOR
!   TRANSMISSION AND BAND-STRUCTURE CALCULATIONS OF PHOTONIC
!   CRYSTALS.  N. STEFANOU, V. YANNOPAPAS, A. MODINOS.
!EF. IN COMP. PHYS. COMMUN. 132 (2000) 189
!EADME
!!The program is in one file which contains:
!1.   The FORTRAN source code (MULTEM)
!2.   Two input data files  (one  for the transmission test run and
!     one for the band-structure test run)
!3.   The corresponding two files with  the  test run outputs which
!     will be produced by running the  program using the above sets
!     of input data.
!The comment lines starting with  "CCCCCCCCC----------->"  identify
!the different constituent parts. These lines must  be removed when
!you will extract the input data files.
!!Input data are read from unit 10.
!No special  job control  statements are  needed.  One has just to
!compile the FORTRAN  source code  and  execute  the  so  produced
!executable file, having  in the same  directory the desired input
!data in the file in unit 10 (named "ftn10" in HP-UX for instance)
!CCCCCCCC-----------> END OF README
!CCCCCCCC-----------> THIS IS THE FIRST LINE OF THE FILE
!CCCCCCCC-----------> HERE STARTS THE FORTRAN SOURCE CODE
!=======================================================================

program multem
    use libmultem2b
    implicit none
    !     ------------------------------------------------------------------
    !     A B S T R A C T
    !     this program calculates either the absorbance, reflectivity  and
    !     transmittance  of   light  by  a   finite  slab   consisting  of
    !     homogeneous   plates and   multilayers  of  spherical  particles
    !     arranged in  a two-dimensional  bravais lattice, or the  complex
    !     photonic  band structure of such an infinite periodic structure.
    !
    !     D E S C R I P T I O N    O F    I N P U T    D A T A
    !     ktype=     1: the direction of an incident  em wave is specified
    !                   by the polar angles of incidence "theta" and "fi".
    !                   the program calculates the transmission,reflection
    !                   and  absorption   coefficients of  a  finite  slab
    !                2: the direction of  an incident em wave is specified
    !                   by the components  of the wavevector  parallel  to
    !                   the  interfaces of the structure:
    !                   aq(1) and aq(2) (and the  frequency). the
    !                   program  calculates  the transmission, reflection,
    !                   absorption coefficients of a finite slab
    !                3: the program calculates  the photonic  complex band
    !                   structure of such  an infinite periodic  structure
    !                   for a  wavevector with components parallel to  the
    !                   interfaces of the structure: aq(1) and aq(2)
    !     kscan=     1: scanning over frequencies
    !                2: scanning over wavelengths
    !     kemb        : indicates the presence (=1) or absence (=0) of a
    !                   different embedding medium
    !     lmax        : cutoff in spherical waves expansions
    !     ncomp       : number of different components in the unit slice.
    !                   their type is specified  by the integer array
    !                   it(icomp)
    !     it=        1: homogeneous plate of thickness "d"
    !                2: multilayer  of spherical  particles arranged in  a
    !                   2d  bravais lattice.each layer consists of "nplan"
    !                   non-primitive  planes of spheres with the same 2-d
    !                   periodicity. the number of unit layers is equal to
    !                   2**(nlayer-1).
    !     dl, dr      : position vectors indicating the origin on the left
    !                   and on the right of the  unit,  respectively. both
    !                   are directed from left to right.
    !     al          : primitive  translation  vector  of the  unit slice
    !                   (effective only for band structure calculation).it
    !                   is given in program units.
    !     nunit       : specifies the number of unit slices (2**(nunit-1))
    !                   of the sample
    !     alpha,alphap: length of primitive vectors of the two-dimensional
    !                   lattice. in program units the size of alpha serves
    !                   as  the unit length.  thus  alpha must be equal to
    !                   1.d0
    !     fab         : angle (in deg) between alpha and alphap
    !     rmax        : upper limit for the length of  reciprocal  lattice
    !                   vectors (in units of 1/alpha) which  must be taken
    !                   into account
    !     zinf,zsup   : minimum  and  maximum  values  of  frequency   (in
    !                   program units: omega*alpha/c), or  wavelength  (in
    !                   program units: lamda/alpha  ),  according  to  the
    !                   value of kscan. c and lamda refer to vacuum
    !     np          : number of equally spaced points between zinf, zsup
    !     polar       : polarization ('s ' or  'p ') of the incident light
    !     aq(1,2)     : wavevector components parallel  to the  interfaces
    !                   of the structure (xy-plane) in units of 2*pi/alpha
    !     theta,fi    : polar angles of incidence (in deg) of the incident
    !                   light
    !     fein        : angle  (in deg) specifying  the direction  of  the
    !                   polarization  vector  for  normal  incidence.  not
    !                   effective otherwise
    !     eps*,mu*    : relative dielectric functions and magnetic permea-
    !                   bilities of the various media
    !     ------------------------------------------------------------------
    !
    ! ..  parameter statements ..
    !
    integer   lmaxd, igd, igkd, nelmd, ncompd, npland
    parameter (lmaxd = 14, igd = 21, igkd = 2 * igd, &
            nelmd = 165152, ncompd = 8, npland = 4)
    !
    ! ..  scalar variables ..
    !
    integer      lmax, i, igkmax, igk1, igk2, igmax, ktype, kscan, ncomp, ig1
    integer      n, np, ig0, nunit, icomp, kemb, iu, ipl, ilayer
    real(dp)     alpha, emach, epsilon
    real(dp)     a0, ra0, rmax, akxy
    real(dp)     zval, zstep, zinf, zsup, fab, alphap, theta, fi, fein
    complex(dp)   kappa, kappa0, akzin, muembl, epsembl
    complex(dp)   muembr, epsembr, d2, kapout
    complex(dp)   kappal, kappar, kappasl, d1, kapin, kapl, kapr
    complex(dp)   mlast, elast, mfirst, efirst, rap
    character(2)  polar
    character(17) text1(2)
    character(5)  dummy
    !
    ! ..  array variables ..
    !
    integer    nt1(igd), nt2(igd), it(ncompd)
    integer    nlayer(ncompd), nplan(ncompd)
    real(dp)   elm(nelmd), ak(2), vecmod(igd), dl(3, ncompd, npland)
    real(dp)   dr(3, ncompd, npland), g(2, igd), b1(2), b2(2)
    real(dp)   s(ncompd, npland), al(3), d(ncompd), vec0(3), aq(2)
    complex(dp) qil  (igkd, igkd), qiil(igkd, igkd), qiiil(igkd, igkd)
    complex(dp) qivl (igkd, igkd), qir (igkd, igkd), qiir (igkd, igkd)
    complex(dp) qiiir(igkd, igkd), qivr(igkd, igkd), wivl (igkd, igkd)
    complex(dp) wil  (igkd, igkd), wiil(igkd, igkd), wiiil(igkd, igkd)
    complex(dp) eincid(igkd), ein(2), eps2(ncompd), eps3(ncompd)
    complex(dp) mu1(ncompd), mu2(ncompd), mu3(ncompd), eps1(ncompd)
    complex(dp) musph(ncompd, npland), epssph(ncompd, npland)
    !
    ! ..  common blocks ..
    !
    real(dp)   ar1(2), ar2(2)
    !      common/x1/ar1,ar2
    !
    ! ..  data statements ..
    !
    data emach/1.d-8/, epsilon/0.d0/
    data eincid/igkd*(0.d0, 0.d0)/, vec0/3*0.d0/
    data text1/'homogeneous plate', 'photonic crystal'/
    !     ------------------------------------------------------------------
    !
    read(10, 200) ktype, kscan, kemb, lmax, ncomp, nunit
    if(ktype<=0.or.ktype>=4) stop 'illegal input value of ktype'
    if(kscan<=0.or.kscan>=3) stop 'illegal input value of kscan'
    if(kemb<0.or.kemb>=2)   stop 'illegal input value of kemb '
    if(lmax<=0.or.lmax>lmaxd.or.lmax>14)&
            stop 'lmax<=0.or.lmax>min0(14,lmaxd)'
    if(ncomp<=0.or.ncomp>ncompd)&
            stop 'illegal input value of ncomp'
    if(nunit<=0)           stop 'illegal input value of nunit'
    read(10, 202) alpha, alphap, fab, rmax
    fab = fab * pi / 180.0_dp
    read(10, 203) np, zinf, zsup
    if(np<=1)                  stop 'illegal input value of  np '
    if(ktype>=2) then
        read(10, 204) aq(1), aq(2), polar, fein
        fein = fein * pi / 180.d0
        aq(1) = 2.d0 * pi * aq(1)
        aq(2) = 2.d0 * pi * aq(2)
        if(ktype<3) then
            write(6, 222)
        else
            write(6, 223)
        endif
        if(ktype==2) write(6, 207) aq(1), aq(2), polar
        if(ktype==3) write(6, 225) aq(1), aq(2)
    else
        read(10, 204) theta, fi, polar, fein
        write(6, 208) theta, fi, polar
        fein = fein * pi / 180.d0
        theta = theta * pi / 180.d0
        fi = fi * pi / 180.d0
    endif
    do icomp = 1, ncomp
        read(10, 201) it(icomp)
        if(it(icomp)<=0.or.it(icomp)>2)&
                stop 'illegal component type'
        write(6, 209) icomp, text1(it(icomp))
        if(it(icomp)==1) then
            read(10, 204) d(icomp)
            read(10, 205) mu1(icomp), eps1(icomp), mu2(icomp), eps2(icomp), &
                    mu3(icomp), eps3(icomp)
            write(6, 210) mu1(icomp), mu2(icomp), mu3(icomp), eps1(icomp), &
                    eps2(icomp), eps3(icomp)
            read(10, *) dummy, (dl(i, icomp, 1), i = 1, 3)
            read(10, *) dummy, (dr(i, icomp, 1), i = 1, 3)
        else
            read(10, 205) mu1(icomp), eps1(icomp)
            if(dble(mu1(icomp))<=0.d0.or.dble(eps1(icomp))<=0.d0)&
                    then
                write(6, 226)
                stop
            endif
            read(10, 201) nplan(icomp), nlayer(icomp)
            do ipl = 1, nplan(icomp)
                read(10, 206) s(icomp, ipl), musph(icomp, ipl), epssph(icomp, ipl)
                read(10, *) dummy, (dl(i, icomp, ipl), i = 1, 3)
                read(10, *) dummy, (dr(i, icomp, ipl), i = 1, 3)
            end do
            write(6, 211)  mu1(icomp), (musph(icomp, ipl), ipl = 1, nplan(icomp))
            write(6, 220) eps1(icomp), (epssph(icomp, ipl), ipl = 1, nplan(icomp))
            write(6, 224) (s(icomp, ipl), ipl = 1, nplan(icomp))
            write(6, 212) 2**(nlayer(icomp) - 1)
        endif
    end do

    d1 = sqrt(mu1(1) * eps1(1))
    d2 = sqrt(mu1(ncomp) * eps1(ncomp))
    if(it(ncomp)==1) d2 = sqrt(mu3(ncomp) * eps3(ncomp))
    if(aimag(d1)/=0.d0) then
        write(6, 227)
        stop
    endif
    if(aimag(d2)/=0.d0) then
        write(6, 228)
        stop
    endif
    if(ktype/=3) then
        write(6, 221) 2**(nunit - 1)
        if(kemb==1) then
            read(10, 205) muembl, epsembl
            read(10, 205) muembr, epsembr
            d1 = sqrt(muembl * epsembl)
            d2 = sqrt(muembr * epsembr)
            if(aimag(d1)/=0.d0) then
                write(6, 227)
                stop
            endif
            if(aimag(d2)/=0.d0) then
                write(6, 228)
                stop
            endif
        endif
    else
        read(10, *) dummy, (al(i), i = 1, 3)
    endif
    call elmgen(elm, nelmd, lmax)
    !
    !****** define the 2d direct and reciprocal-lattice vectors ******
    !
    ar1(1) = alpha
    ar1(2) = 0.d0
    ar2(1) = alphap * cos(fab)
    ar2(2) = alphap * sin(fab)
    write(6, 213) ar1(1), ar1(2), ar2(1), ar2(2)
    a0 = abs(ar1(1) * ar2(2) - ar1(2) * ar2(1))
    ra0 = 2.d0 * pi / a0
    b1(1) = -ar1(2) * ra0
    b1(2) = ar1(1) * ra0
    b2(1) = -ar2(2) * ra0
    b2(2) = ar2(1) * ra0
    call lat2d(b1, b2, rmax, igmax, igd, nt1, nt2, vecmod)
    write(6, 214) b1(1), b1(2), b2(1), b2(2)
    igkmax = 2 * igmax
    do ig1 = 1, igmax
        g(1, ig1) = nt1(ig1) * b1(1) + nt2(ig1) * b2(1)
        g(2, ig1) = nt1(ig1) * b1(2) + nt2(ig1) * b2(2)
        write(6, 215) ig1, nt1(ig1), nt2(ig1), vecmod(ig1)
    end do
    zstep = (zsup - zinf) / dble(np - 1)
    zval = zinf - zstep
    if(ktype<3) then
        if(kscan==1) write(6, 216)
        if(kscan==2) write(6, 217)
        if(polar/='s '.and.polar/='p ') stop 'illegal polarization'
    else
        if(kscan==1) write(6, 218)
        if(kscan==2) write(6, 219)
    endif
    if(polar=='p ') then
        ein(1) = cone
        ein(2) = czero
    else
        ein(1) = czero
        ein(2) = cone
    end if
    do n = 1, np   !****** scanning over frequencies/wavelengths ******
        zval = zval + zstep
        if(kscan==1) kappa0 = cmplx_dp(zval, epsilon)
        if(kscan==2) kappa0 = cmplx_dp(2.d0 * pi / zval, epsilon)
        kapin = kappa0 * d1
        kapout = kappa0 * d2
        if(ktype==1) then
            ak(1) = dble(kapin) * sin(theta) * cos(fi)
            ak(2) = dble(kapin) * sin(theta) * sin(fi)
            do i = 1, igkmax
                eincid(i) = czero
            end do
        else
            ak(1) = aq(1)
            ak(2) = aq(2)
        endif
        if(ktype/=3) then !define the polarization vector from "ak"*****
            akxy = ak(1) * ak(1) + ak(2) * ak(2)
            akzin = sqrt(kapin * kapin - akxy)
            if(dble(akzin)<emach)      stop 'improper incident wave'
            akxy = sqrt(akxy)
            if(akxy<emach) then
                ein(1) = cmplx_dp(cos(fein), 0.d0)
                ein(2) = cmplx_dp(sin(fein), 0.d0)
            end if
            call reduce(ar1, ar2, ak, igmax, g, ig0, emach)   !"ak" in sbz*******
            do i = 1, 2
                eincid(2 * ig0 - 2 + i) = ein(i)
            end do
        else
            call reduce(ar1, ar2, ak, igmax, g, ig0, emach)   !"ak" in sbz*******
        endif
        !
        !****** construct the transfer matrix of the unit slice ******
        !
        if(it(1)==1) then
            kappal = sqrt(mu1(1) * eps1(1)) * kappa0
            kappasl = sqrt(mu2(1) * eps2(1)) * kappa0
            kappar = sqrt(mu3(1) * eps3(1)) * kappa0
            kapl = kappal
            kapr = kappar
            mlast = mu3(1)
            elast = eps3(1)
            mfirst = mu1(1)
            efirst = eps1(1)
            call hoslab(igmax, kappal, kappasl, kappar, ak, g, dl(1, 1, 1), &
                    dr(1, 1, 1), d(1), qil, qiil, qiiil, qivl, emach)
        else
            kappa = sqrt(mu1(1) * eps1(1)) * kappa0
            kapl = kappa
            kapr = kappa
            mlast = mu1(1)
            elast = eps1(1)
            mfirst = mu1(1)
            efirst = eps1(1)
            rap = s(1, 1) * kappa0 / 2.d0 / pi
            call pcslab(lmax, igmax, rap, eps1(1), epssph(1, 1), mu1(1), musph(1, 1), &
                    kappa, ak, dl(1, 1, 1), dr(1, 1, 1), g, elm, a0, emach, &
                    qil, qiil, qiiil, qivl, ar1, ar2)
            if(nplan(1)>=2) then
                do ipl = 2, nplan(1)
                    rap = s(1, ipl) * kappa0 / 2.d0 / pi
                    call pcslab(lmax, igmax, rap, eps1(1), epssph(1, ipl), mu1(1), &
                            musph(1, ipl), kappa, ak, dl(1, 1, ipl), dr(1, 1, ipl), &
                            g, elm, a0, emach, qir, qiir, qiiir, qivr, ar1, ar2)
                    call pair(igkmax, igkd, qil, qiil, qiiil, qivl, qir, qiir, qiiir, qivr)
                end do
            endif
            if(nlayer(1)>=2) then
                do ilayer = 1, nlayer(1) - 1
                    call pair(igkmax, igkd, qil, qiil, qiiil, qivl, qil, qiil, qiiil, qivl)
                end do
            endif
        endif
        if(ncomp>=2) then
            do icomp = 2, ncomp
                if(it(icomp)==1) then
                    kappal = sqrt(mu1(icomp) * eps1(icomp)) * kappa0
                    kappasl = sqrt(mu2(icomp) * eps2(icomp)) * kappa0
                    kappar = sqrt(mu3(icomp) * eps3(icomp)) * kappa0
                    kapr = kappar
                    if(abs(mu1(icomp) - mlast)/=0.d0.or.abs(eps1(icomp) - elast)/=&
                            0.d0) stop 'improper matching of successive host media'
                    mlast = mu3(icomp)
                    elast = eps3(icomp)
                    call hoslab(igmax, kappal, kappasl, kappar, ak, g, dl(1, icomp, 1), &
                            dr(1, icomp, 1), d(icomp), qir, qiir, qiiir, qivr, emach)
                else
                    kappa = sqrt(mu1(icomp) * eps1(icomp)) * kappa0
                    kapr = kappa
                    if(abs(mu1(icomp) - mlast)/=0.d0.or.abs(eps1(icomp) - elast)/=&
                            0.d0) stop 'improper matching of successive host media'
                    mlast = mu1(icomp)
                    elast = eps1(icomp)
                    rap = s(icomp, 1) * kappa0 / 2.d0 / pi
                    call pcslab(lmax, igmax, rap, eps1(icomp), epssph(icomp, 1), mu1(icomp)&
                            , musph(icomp, 1), kappa, ak, dl(1, icomp, 1), &
                            dr(1, icomp, 1), g, elm, a0, emach, qir, qiir, qiiir, qivr, &
                            ar1, ar2)
                    if(nplan(icomp)>=2) then
                        do igk1 = 1, igkmax
                            do igk2 = 1, igkmax
                                wil  (igk1, igk2) = qir  (igk1, igk2)
                                wiil (igk1, igk2) = qiir (igk1, igk2)
                                wiiil(igk1, igk2) = qiiir(igk1, igk2)
                                wivl (igk1, igk2) = qivr (igk1, igk2)
                            end do
                        end do
                        do ipl = 2, nplan(icomp)
                            rap = s(icomp, ipl) * kappa0 / 2.d0 / pi
                            call pcslab(lmax, igmax, rap, eps1(icomp), epssph(icomp, ipl), &
                                    mu1(icomp), musph(icomp, ipl), kappa, ak, &
                                    dl(1, icomp, ipl), dr(1, icomp, ipl), g, elm, a0, emach, &
                                    qir, qiir, qiiir, qivr, ar1, ar2)
                            call pair(igkmax, igkd, wil, wiil, wiiil, wivl, qir, qiir, qiiir, qivr)
                        end do
                        do igk1 = 1, igkmax
                            do igk2 = 1, igkmax
                                qir  (igk1, igk2) = wil  (igk1, igk2)
                                qiir (igk1, igk2) = wiil (igk1, igk2)
                                qiiir(igk1, igk2) = wiiil(igk1, igk2)
                                qivr (igk1, igk2) = wivl (igk1, igk2)
                            end do
                        end do
                    endif
                    if(nlayer(icomp)>=2) then
                        do ilayer = 1, nlayer(icomp) - 1
                            call pair(igkmax, igkd, qir, qiir, qiiir, qivr, qir, qiir, qiiir, qivr)
                        end do
                    endif
                endif
                call pair(igkmax, igkd, qil, qiil, qiiil, qivl, qir, qiir, qiiir, qivr)
            end do
        endif
        if(ktype<3) then
            !
            !****** the unit slice is defined. this can be repeated by the ******
            !****** doubling-layer  technique, interfaces can be added and ******
            !****** reflectivity/transmittance/absorbance are calculated.  ******
            !
            if(nunit/=1) then
                if(abs(mlast - mfirst)/=0.d0.or.abs(elast - efirst)/=0.d0)&
                        stop 'improper matching of successive host media'
                do iu = 1, nunit - 1
                    call pair(igkmax, igkd, qil, qiil, qiiil, qivl, &
                            qil, qiil, qiiil, qivl)
                end do
            end if
            if(kemb==1) then
                call hoslab(igmax, kapr, (kapr + kapout) / 2.d0, kapout, ak, g, vec0, &
                        vec0, 0.d0, qir, qiir, qiiir, qivr, emach)
                call pair(igkmax, igkd, qil, qiil, qiiil, qivl, qir, qiir, qiiir, qivr)
                do igk1 = 1, igkmax
                    do igk2 = 1, igkmax
                        qir  (igk1, igk2) = qil  (igk1, igk2)
                        qiir (igk1, igk2) = qiil (igk1, igk2)
                        qiiir(igk1, igk2) = qiiil(igk1, igk2)
                        qivr (igk1, igk2) = qivl (igk1, igk2)
                    end do
                end do
                call hoslab(igmax, kapin, (kapl + kapin) / 2.d0, kapl, ak, g, vec0, &
                        vec0, 0.d0, qil, qiil, qiiil, qivl, emach)
                call pair(igkmax, igkd, qil, qiil, qiiil, qivl, qir, qiir, qiiir, qivr)
            endif
            call scat(igmax, zval, ak, g, dble(kapin), dble(kapout), &
                    eincid, qil, qiiil)
        else
            !
            !****** alternatively, calculate complex photonic band structure ******
            !
            if(abs(mlast - mfirst)/=0.d0.or.abs(elast - efirst)/=0.d0)&
                    stop 'improper matching of successive host media'
            call band(igmax, zval, emach, ak, al, qil, qiil, qiiil, qivl)
        endif
    end do
    stop
    200 format(///, 6(10x, i2))
    201 format(6(10x, i2))
    202 format(4(8x, f12.6))
    203 format(6x, i4, 2(8x, f13.8))
    204 format(2(15x, f13.8), 10x, a2, 10x, f7.2///)
    205 format(2(12x, 2f13.8))
    206 format(10x, f13.8, 2(12x, 2f13.8))
    207 format(3x, 'k_parallel=', 2f12.6, 5x, a2, 'polarization')
    208 format(3x, 'angles of incidence (in rad):  theta=', f7.2, 3x, 'fi=', &
            f7.2, 5x, a2, 'polarization')
    209 format(3x, 'component nr.', i2, 3x, 'type:', 2x, a17)
    210 format(3x, 'mu :', 2f10.5, ' | ', 2f10.5, ' | ', 2f10.5/&
            3x, 'eps:', 2f10.5, ' | ', 2f10.5, ' | ', 2f10.5)
    211 format(3x, 'mu :', 2f10.5, ' | ', 4(3x, 2f10.5))
    212 format(29x, i6, ' unit layers')
    213 format(/13x, 'primitive lattice vectors'/13x, 'ar1 = (', 2f12.4, ')'/&
            13x, 'ar2 = (', 2f12.4, ')')
    214 format(13x, 'unit vectors in reciprocal space:'/13x, 'b1  = (', &
            2f12.4, ')'/13x, 'b2  = (', 2f12.4, ')'//3x, 'reciprocal vectors', 5x, &
            'length')
    215 format(i3, 4x, 2i5, 5x, e14.6)
    216 format(//4x, 'frequency   transmittance  reflectance   absorbance'&
            /60('-'))
    217 format(//4x, 'wavelength  transmittance  reflectance   absorbance'&
            /60('-'))
    218 format(//1x, 'frequency', 7x, 'normalized k_z'/&
            1x, 9('-'), 7x, 14('-'))
    219 format(//3x, 'wavelength versus normalized k_z'/3x, 32('-'))
    220 format(3x, 'eps:', 2f10.5, ' | ', 4(3x, 2f10.5))
    221 format(3x, 'the sample consists of ', i6, ' unit slices')
    222 format(5x, '****************************************************'/&
            5x, '*** output: transmittance/reflectance/absorbance ***'/&
            5x, '****************************************************')
    223 format(5x, '****************************************************'/&
            5x, '************** output: band structure **************'/&
            5x, '****************************************************')
    224 format(3x, '  s:', 23x, 4(3x, f10.5, 10x))
    225 format(3x, 'k_parallel=', 2f12.6)
    226 format(5x, '----------------------------------'/&
            5x, 'illegal input: spheres embedded in'/&
            5x, 'a  medium  of negative  dielectric'/&
            5x, 'constant. the ewald  summation  in'/&
            5x, 'subroutine xmat does not converge.'/&
            5x, 'direct - space summation is needed'/&
            5x, 'instead.'/&
            5x, '----------------------------------')
    227 format(5x, '----------------------------------'/&
            5x, 'illegal input:semi-infinite medium'/&
            5x, 'of complex refractive index on the'/&
            5x, 'left side of the slab.'/&
            5x, '----------------------------------')
    228 format(5x, '----------------------------------'/&
            5x, 'illegal input:semi-infinite medium'/&
            5x, 'of complex refractive index on the'/&
            5x, 'right side of the slab.'/&
            5x, '----------------------------------')
end program