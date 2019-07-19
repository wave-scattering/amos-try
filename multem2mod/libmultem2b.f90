module libmultem2b
    use dense_solve
    use amos
    implicit none
    private
    integer, parameter:: dp=kind(0.d0)
    complex(dp), parameter :: czero = (0.0_dp, 0.0_dp)
    complex(dp), parameter :: ci    = (0.0_dp, 1.0_dp)
    complex(dp), parameter :: cone  = (1.0_dp, 0.0_dp)
    complex(dp), parameter :: ctwo  = (2.0_dp, 0.0_dp)
    real(dp), parameter :: pi=4.0_dp*ATAN(1.0_dp)
    public bessel, tmtrx, sphrm4
contains
    !=======================================================================
    !=======================================================================
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