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
CCCCCCCCC-----------> END OF README
CCCCCCCCC-----------> THIS IS THE FIRST LINE OF THE FILE
CCCCCCCCC-----------> HERE STARTS THE FORTRAN SOURCE CODE
C=======================================================================
      module libmultem2a
      use dense_solve
      use libmultem2b
      use multem_blas
      implicit none
      private
      public cerf, ceven, codd, band, scat, pcslab,
     & reduce, lat2d
      contains
C=======================================================================
C======================================================================
C=======================================================================
C=======================================================================
      subroutine xmat(xodd,xeven,lmax,kappa,ak,elm,emach,ar1, ar2)
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
      parameter (ndend=1240)
! ..  scalar arguments  ..
      integer    lmax
      real(dp)   emach
      complex(dp) kappa
! ..  array arguments  ..
      real(dp)  ar1(2),ar2(2)
      real(dp)   ak(2),elm(:)
      complex(dp) xodd(:,:),xeven(:,:)
! ..  local scalars  ..
      integer    l2max,ll2,ii,i,nndlm,k,kk,l,mm,nn,m,j1,j2,i1,i2,i3,n1
      integer    na,lll,n,il,nm,in,l2,il2,m2,il3,l3,m3,la1,lb1,la11,lb11
      integer    ll,j,l1
      real(dp)   ab1,ab2,ac,acsq,ad,al,an,an1,an2,ap,ap1,ap2,ar,b
      real(dp)   dnorm,rtpi,rtv,test,test1,test2,tv
      complex(dp) alpha,rta,rtai,kapsq,kant,knsq,xpk,xpa,cf,cp,cx,cz
      complex(dp) z,zz,w,ww,a,acc,gpsq,gp,bt,aa,ab,u,u1,u2,gam
      complex(dp) gk,gkk,sd,alm
!
! ..  local arrays  ..
!
      real(dp)   denom(ndend),r(2),b1(2),b2(2),akpt(2),fac(4*lmax+1)
      complex(dp) gkn(lmax+1),agk(2*lmax+1),xpm(2*lmax+1),
     &  pref((lmax+1)**2)
      complex(dp) dlm((lmax+1)*(2*lmax+1))
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
      rtpi=sqrt(pi)
      kapsq=kappa*kappa
!
!     the factorial  function  is tabulated  in fac . the array
!     dlm will contain non-zero,i.e. l+m even,values as defined
!     by kambe.dlm=dlm1+dlm2+dlm3.with lm=(00),(1-1),(11),(2-2)...
!
      l2max=lmax+lmax
      ll2=l2max+1
      fac(1)=1.0_dp
      ii=l2max+l2max
      do i=1,ii
        fac(i+1)=dble(i)*fac(i)
      end do
      nndlm=l2max*(l2max+3)/2+1
      do i=1,nndlm
        dlm(i)=czero
      end do
!
!     the formula of kambe for the separation constant,alpha,is
!     used,subject to a restriction which is imposed to control
!     later rounding errors
!
      tv=abs(ar1(1)*ar2(2)-ar1(2)*ar2(1))
      alpha=tv/(4.0_dp*pi)*kapsq
      al=abs(alpha)
      if(exp(al)*emach-5.0d-5 > 0) al=log(5.0d-5/emach)
      alpha=cmplx_dp(al,0.0_dp)
      rta=sqrt(alpha)
!
!     dlm1 , the  sum  over  reciprocal   lattice  vectors  , is
!     calculated first. the prefactor p1 is  tabulated  for even
!     values of l+|m|,thus lm=(00),(11),(2 0),(22),(l2max,l2max)
!     the  factorial  factor  f1  is simultaneously tabulated in
!     denom,for all values of n=0,(l-|m|)/2
!
      k=1
      kk=1
      ap1=-2.0_dp/tv
      ap2=-1.0_dp
      cf=ci/kappa
      do l=1,ll2
        ap1=ap1/2.0_dp
        ap2=ap2+2.0_dp
        cp=cf
        mm=1
        if(mod  (l,2) == 0) then
        mm=2
        cp=ci*cp
        end if
        nn=(l-mm)/2+2
        do m=mm,l,2
          j1=l+m-1
          j2=l-m+1
          ap=ap1*sqrt(ap2*fac(j1)*fac(j2))
          pref(kk)=ap*cp
          cp=-cp
          kk=kk+1
          nn=nn-1
          do i=1,nn
            i1=i
            i2=nn-i+1
            i3=nn+m-i
            denom(k)=1.0_dp/(fac(i1)*fac(i2)*fac(i3))
            k=k+1
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
      rtv=2.0_dp*pi/tv
      b1(1)=-ar1(2)*rtv
      b1(2)=ar1(1)*rtv
      b2(1)=-ar2(2)*rtv
      b2(2)=ar2(1)*rtv
      test1=1.0d6
      ii=1
      n1=-1
   9  n1=n1+1
      na=n1+n1+ii
      an1=dble(n1)
      an2=-an1-1.0_dp
      do i1=1,na
        an2=an2+1.0_dp
        do i2=1,4
!     write(16,307) i1,i2
! 307 format(33x,'i1=',i2,' , i2=',i2/33x,12('='))
          an=an1
          an1=-an2
          an2=an
          ab1=an1*b1(1)+an2*b2(1)
          ab2=an1*b1(2)+an2*b2(2)
          akpt(1)=ak(1)+ab1
          akpt(2)=ak(2)+ab2
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
          acsq=akpt(1)*akpt(1)+akpt(2)*akpt(2)
          gpsq=kapsq-acsq
          if(abs(gpsq)<emach*emach)   then
          write(7,100)
  100     format(13x,'fatal error from xmat:'/3x,'gpsq is too small.'
     &     /3x,'give a small but nonzero value for "epsilon"'/3x,
     &     'in the data statement of the main program.'
     &     /3x,'this defines a small imaginary part'
     &     /3x,'in the frequency or wavelength value.')
          stop
          endif
          ac=sqrt(acsq)
          gp=sqrt(gpsq)
          xpk=czero
          gk=czero
          gkk=cmplx_dp(1.0_dp,0.0_dp)
          if(ac-emach > 0) then
            xpk=cmplx_dp(akpt(1)/ac,akpt(2)/ac)
            gk=ac/kappa
            gkk=gpsq/kapsq
          end if
          xpm(1)=cmplx_dp(1.0_dp,0.0_dp)
          agk(1)=cmplx_dp(1.0_dp,0.0_dp)
          do i=2,ll2
            xpm(i)=xpm(i-1)*xpk
            agk(i)=agk(i-1)*gk
          end do
          cf=kappa/gp
          zz=-alpha*gkk
          cz=sqrt(-zz)
          z=-ci*cz
          cx=exp(-zz)
          gam=rtpi*cerf(cz)
          gkn(1)=cf*cx*gam
          bt=z
          b=0.5_dp
          lll=l2max/2+1
          do i=2,lll
            bt=bt/zz
            b=b-1.0_dp
            gam=(gam-bt)/b
            cf=cf*gkk
            gkn(i)=cf*cx*gam
          end do
!
!     the contribution to the sum dlm1 for a particular
!     reciprocal lattice vector is now accumulated into
!     the  elements of dlm,note special action if  ac=0
!
          k=1
          kk=1
          do l=1,ll2
            mm=1
            if(mod  (l,2) == 0) mm=2
            n=(l*l+mm)/2
            nn=(l-mm)/2+2
            do m=mm,l,2
              acc=czero
              nn=nn-1
              il=l
              do i=1,nn
                acc=acc+denom(k)*agk(il)*gkn(i)
                il=il-2
                k=k+1
              end do
              acc=pref(kk)*acc
              if(ac-1.0d-6)17,17,165
 165          dlm(n)=dlm(n)+acc/xpm(m)
              if(m-1)17,18,17
  17          nm=n-m+1
              dlm(nm)=dlm(nm)+acc*xpm(m)
  18          kk=kk+1
              n=n+1
            end do
          end do
          if(ii > 0) exit
        end do
      ii=0
      end do
!
!     after each step of the summation a test on the
!     convergence  of the  elements of  dlm is  made
      test2=0.0_dp
      do i=1,nndlm
        dnorm=abs(dlm(i))
        test2=test2+dnorm*dnorm
      end do
      test=abs((test2-test1)/test1)
      test1=test2
      if(test-0.001d0 > 0) then
      if(n1-10)9,25,25
  25  write(16,26)n1
  26  format(//13x, 'dlm1,s not converged by n1=', i2)
      else
        write(16,28)n1
  28    format(//13x, 'dlm1,s converged by n1=', i2)
      end if
!     write(16,250)dlm
!250  format(5h0dlm1,//,45(2e13.5,/))
!
!     dlm2, the sum over real space lattice vectors, begins with
!     the adjustment of the array pref, to contain values of the
!     prefactor  'p2' for lm=(00),(11),(20),(22),...
!
      kk=1
      ap1=tv/(4.0_dp*pi)
      cf=kapsq/ci
      do l=1,ll2
        cp=cf
        mm=1
        if(mod  (l,2) == 0) then
          mm=2
          cp=-ci*cp
        end if
        j1=(l-mm)/2+1
        j2=j1+mm-1
        in=j1+l-2
        ap2=((-1.0_dp)**in)*ap1
        do m=mm,l,2
          ap=ap2/(fac(j1)*fac(j2))
          pref(kk)=ap*cp*pref(kk)
          j1=j1-1
          j2=j2+1
          ap2=-ap2
          cp=-cp
          kk=kk+1
        end do
      end do
!
!     the summation proceeds in steps of 8*n1 lattice points
!     as before, but this  time excluding  the origin  point
!     r=the current lattice vector in the sum
!     ar=mod(r)
!
      n1=0
  32  n1=n1+1
      na=n1+n1
      an1=dble(n1)
      an2=-an1-1.0_dp
      do i1=1,na
        an2=an2+1.0_dp
        do i2=1,4
          an=an1
          an1=-an2
          an2=an
          r(1)=an1*ar1(1)+an2*ar2(1)
          r(2)=an1*ar1(2)+an2*ar2(2)
          ar=sqrt(r(1)*r(1)+r(2)*r(2))
          xpk=cmplx_dp(r(1)/ar,r(2)/ar)
          xpm(1)=cmplx_dp(1.0_dp,0.0_dp)
          do i=2,ll2
            xpm(i)=xpm(i-1)*xpk
          end do
          ad=ak(1)*r(1)+ak(2)*r(2)
          sd=exp(-ad*ci)
!
!     for each lattice vector the integral 'u' is obtained
!     from the recurrence relation in l suggested by kambe
!     u1 and u2 are the  initial terms of this recurrence,
!     for l#-1 and l=0, and they are evaluated in terms of
!     the complex error function cerf
!
          kant=0.5_dp*ar*kappa
          knsq=kant*kant
          z=ci*kant/rta
          zz=rta-z
          z=rta+z
          ww=cerf(-zz)
          w=cerf(z)
          aa=0.5_dp*rtpi*(w-ww)/ci
          ab=0.5_dp*rtpi*(w+ww)
          a=alpha-knsq/alpha
          xpa=exp(a)
          u1=aa*xpa
          u2=ab*xpa/kant
!
!     the contribution to dlm2 from a particular lattice
!     vector  is  accumulated into  the elements of  dlm
!     this procedure includes the term (kant**l) and the
!     recurrence for the integral 'u'
!
          kk=1
          al=-0.5_dp
          cp=rta
          cf=cmplx_dp(1.0_dp,0.0_dp)
          do l=1,ll2
            mm=1
            if(mod  (l,2) == 0) mm=2
            n=(l*l+mm)/2
            do m=mm,l,2
              acc=pref(kk)*u2*cf*sd
              dlm(n)=dlm(n)+acc/xpm(m)
              if(m-1 /= 0) then
                nm=n-m+1
                dlm(nm)=dlm(nm)+acc*xpm(m)
              end if
              kk=kk+1
              n=n+1
            end do
            al=al+1.0_dp
            cp=cp/alpha
            u=(al*u2-u1+cp*xpa)/knsq
            u1=u2
            u2=u
            cf=kant*cf
          end do
        end do
      end do
!
!     after each step of the summation a test on the
!     convergence of the elements of dlm is made
!
      test2=0.0_dp
      do i=1,nndlm
        dnorm=abs(dlm(i))
        test2=test2+dnorm*dnorm
      end do
      test=abs((test2-test1)/test1)
      test1=test2
      if(test-0.001d0 > 0) then
      if(n1-10) 32,43,43
  43    write(16,44)n1
  44    format(//3x, 'dlm2,s not converged by n1=', i2)
      else
        write(16,46)n1
  46    format(//3x, 'dlm2,s converged by n1=', i2)
      end if
!
!     the term dlm3 has a non-zero contribution  only
!     when l=m=0.it is evaluated here in terms of the
!     complex error function cerf
!
      xpa=exp(-alpha)
      rtai=1.0_dp/(rtpi*rta)
      acc=kappa*(ci*(xpa-cerf(rta))-rtai)/xpa
      ap=-0.5_dp/rtpi
      dlm(1)=dlm(1)+ap*acc
!
!     finally the elements of dlm are multiplied by the
!     factor (-1.0_dp)**((m+|m|)/2)
!
      do l=2,ll2,2
        n=l*l/2+1
        do m=2,l,2
          dlm(n)=-dlm(n)
          n=n+1
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
      k=1
      ii=0
  48  ll=lmax+ii
      i=1
      do il2=1,ll
        l2=il2-ii
        m2=-l2+1-ii
        do i2=1,il2
          j=1
          do il3=1,ll
            l3=il3-ii
            m3=-l3+1-ii
            do i3=1,il3
              alm=czero
              la1=max0(iabs(l2-l3),iabs(m2-m3))
              lb1=l2+l3
              n=(la1*(la1+2)+m2-m3+2)/2
              nn=2*la1+4
              lb11=lb1+1
              la11=la1+1
              do l1=la11,lb11,2
                alm=alm+elm(k)*dlm(n)
                n=n+nn
                nn=nn+4
                k=k+1
              end do
              alm=alm/kappa
              if(i-j == 0) alm=alm+ci
              if(ii <= 0) then
                xodd(j,i)=ci*alm
              else
                xeven(j,i)=ci*alm
              end if
              m3=m3+2
              j=j+1
            end do
          end do
          m2=m2+2
          i=i+1
        end do
      end do
      if(ii<=0) then
        ii=1
        goto 48
      endif
      return
      end subroutine
C=======================================================================
C=======================================================================
C=======================================================================
      SUBROUTINE PCSLAB(LMAX,IGMAX,RAP,EPSMED,EPSSPH,MUMED,MUSPH,KAPPA,
     &                  AK,DL,DR,G,ELM,A0,EMACH,QI,QII,QIII,QIV)

C     ------------------------------------------------------------------
C     THIS SUBROUTINE COMPUTES THE TRANSMISSION/REFLECTION MATRICES FOR
C     A PLANE OF SPHERES EMBEDDED IN A HOMOGENEOUS HOST MEDIUM.
C     ------------------------------------------------------------------
C
C ..  PARAMETER STATEMENTS ..
C
      INTEGER   IGD,NELMD
      PARAMETER (IGD=21,NELMD=165152)
C
C ..  SCALAR ARGUMENTS ..
C
      INTEGER    LMAX,IGMAX
      real(dp)     A0,EMACH
      complex(dp) EPSMED,EPSSPH,MUMED,MUSPH,KAPPA,RAP

C
C ..  ARRAY ARGUMENTS ..
C
      real(dp)     AK(2),DL(3),DR(3),G(2,IGD),ELM(NELMD)
      complex(dp) QI(:,:),QII(:,:),QIII(:,:)
      complex(dp) QIV(:,:)
C
C ..  LOCAL SCALARS  ..
C

      INTEGER    LMTOT,L,M,II,IGK1,IGK2,LMAX1,IGKMAX
      INTEGER    IG1,IG2,ISIGN1,ISIGN2,K1,K2
      INTEGER    LMXOD,IEV,IOD
      real(dp)     SIGN1,SIGN2,SIGNUS
      complex(dp) CQI,CQII,CQIII,CQIV
C
C ..  LOCAL ARRAYS  ..
C
      integer     :: int1((lmax+1)**2-1), int2((lmax+1)**2-1)
      complex(dp) AE(2,(lmax+1)**2),AH(2,(lmax+1)**2),GKK(3,IGD)
      complex(dp) GK(3),LAME(2),LAMH(2)
      complex(dp) XEVEN(((lmax+1)*(lmax+2))/2,((lmax+1)*(lmax+2))/2),
     & XODD((lmax*(lmax+1))/2,(lmax*(lmax+1))/2)
      complex(dp) :: TE(lmax+1),TH(lmax+1)
      complex(dp) BMEL1((LMAX+1)**2-1),BMEL2((LMAX+1)**2-1)
      complex(dp) XXMAT1((LMAX+1)**2-1,(LMAX+1)**2-1),
     & XXMAT2((LMAX+1)**2-1,(LMAX+1)**2-1)
      complex(dp) DLME(2,(lmax+1)**2),DLMH(2,(lmax+1)**2)
C
C ..  COMMON BLOCKS ..
C
      real(dp)     AR1(2),AR2(2)
      COMMON/X1/AR1,AR2
C     ------------------------------------------------------------------
      IGKMAX=2*IGMAX
      LMAX1=LMAX+1
      LMTOT=LMAX1*LMAX1-1
      LMXOD=(LMAX*LMAX1)/2
      DO IG1=1,IGMAX
        GKK(1,IG1)=cmplx((AK(1)+G(1,IG1)),0.0_dp, kind=dp)
        GKK(2,IG1)=cmplx((AK(2)+G(2,IG1)),0.0_dp, kind=dp)
        GKK(3,IG1)=SQRT(KAPPA*KAPPA-GKK(1,IG1)*GKK(1,IG1)-
     &                            GKK(2,IG1)*GKK(2,IG1))
      end do
      CALL TMTRX(RAP,EPSSPH,EPSMED,MUMED,MUSPH,TE,TH)
      CALL XMAT(XODD,XEVEN,LMAX,KAPPA,AK,ELM,EMACH,ar1, ar2)
      CALL SETUP(LMAX,XEVEN,XODD,TE,TH,XXMAT1,XXMAT2)
      call zgetrf_wrap ( XXMAT1, int1 )
      call zgetrf_wrap ( XXMAT2, int2 )
      ISIGN2=1
      SIGN2=3.0_dp-2.0_dp*ISIGN2
      IGK2=0
      DO IG2=1,IGMAX
        GK(1)=      GKK(1,IG2)
        GK(2)=      GKK(2,IG2)
        GK(3)=SIGN2*GKK(3,IG2)
        CALL PLW(KAPPA,GK,LMAX,AE,AH)
        DO K2=1,2
          IGK2=IGK2+1
          II=0
          IEV=LMXOD
          IOD=0
          DO L=1,LMAX
            DO M=-L,L
              II=II+1
              IF(MOD((L+M),2)==0)  THEN
                IEV=IEV+1
                BMEL1(IEV)=TH(L+1)*AH(K2,II+1)
                BMEL2(IEV)=TE(L+1)*AE(K2,II+1)
              ELSE
                IOD=IOD+1
                BMEL1(IOD)=TE(L+1)*AE(K2,II+1)
                BMEL2(IOD)=TH(L+1)*AH(K2,II+1)
              END IF
            end do
          end do

          call zgetrs_wrap(XXMAT1, BMEL1, INT1)
          call zgetrs_wrap(XXMAT2, BMEL2, INT2)
          DO ISIGN1=1,2
            SIGN1=3.0_dp-2.0_dp*ISIGN1
            IGK1=0
            do IG1=1,IGMAX
              GK(1)=      GKK(1,IG1)
              GK(2)=      GKK(2,IG1)
              GK(3)=SIGN1*GKK(3,IG1)
              CALL DLMKG(LMAX,A0,GK,SIGN1,KAPPA,DLME,DLMH,EMACH)
              DO K1=1,2
                LAME(K1)=CZERO;      LAMH(K1)=CZERO
                II=0; IOD=0
                IEV=LMXOD
                do L=1,lmax
                  do M=-L,L
                    II=II+1
                    if(MOD((L+M),2) == 0)  then
                      IEV=IEV+1
                      LAME(K1)=LAME(K1)+DLME(K1,II+1)*BMEL2(IEV)
                      LAMH(K1)=LAMH(K1)+DLMH(K1,II+1)*BMEL1(IEV)
                    else
                      IOD=IOD+1
                      LAME(K1)=LAME(K1)+DLME(K1,II+1)*BMEL1(IOD)
                      LAMH(K1)=LAMH(K1)+DLMH(K1,II+1)*BMEL2(IOD)
                    end if
                  end do
                end do
              end do
              do K1=1,2
                IGK1=IGK1+1
                IF(ISIGN1==1) QI  (IGK1,IGK2)=LAMH(K1)+LAME(K1)
                IF(ISIGN1==2) QIII(IGK1,IGK2)=LAMH(K1)+LAME(K1)
              end do
            end do
          end do
        end do
      end do

      IGK2=0
      DO IG2=1,IGMAX
        DO K2=1,2
          IGK2=IGK2+1
          IGK1=0
          DO IG1=1,IGMAX
            DO K1=1,2
              IGK1=IGK1+1
              SIGNUS=1.0_dp
              IF(K2/=K1) SIGNUS=-1.0_dp
              QII(IGK1,IGK2)=SIGNUS*QIII(IGK1,IGK2)
              QIV(IGK1,IGK2)=SIGNUS*QI  (IGK1,IGK2)
            end do
          end do
        end do
      end do

      DO IGK1=1,IGKMAX
        QI (IGK1,IGK1)=CONE + QI (IGK1,IGK1)
        QIV(IGK1,IGK1)=CONE + QIV(IGK1,IGK1)
      end do

      IGK2=0
      DO IG2=1,IGMAX
        DO IG1=1,IGMAX
      CQI  =EXP(CI*(GKK(1,IG1)*DR(1)+GKK(2,IG1)*DR(2)+GKK(3,IG1)*DR(3)+
     &              GKK(1,IG2)*DL(1)+GKK(2,IG2)*DL(2)+GKK(3,IG2)*DL(3)))
      CQII =EXP(CI*((GKK(1,IG1)-GKK(1,IG2))*DR(1)+(GKK(2,IG1)
     &     -GKK(2,IG2))*DR(2)+(GKK(3,IG1)+GKK(3,IG2))*DR(3)))
      CQIII=EXP(-CI*((GKK(1,IG1)-GKK(1,IG2))*DL(1)+(GKK(2,IG1)
     &     -GKK(2,IG2))*DL(2)-(GKK(3,IG1)+GKK(3,IG2))*DL(3)))
      CQIV =EXP(-CI*(GKK(1,IG1)*DL(1)+GKK(2,IG1)*DL(2)-GKK(3,IG1)*DL(3)+
     &              GKK(1,IG2)*DR(1)+GKK(2,IG2)*DR(2)-GKK(3,IG2)*DR(3)))
          DO K2=1,2
            IGK2=(IG2-1)*2+K2
            DO K1=1,2
              IGK1=(IG1-1)*2+K1
              QI  (IGK1,IGK2)=CQI  *QI  (IGK1,IGK2)
              QII (IGK1,IGK2)=CQII *QII (IGK1,IGK2)
              QIII(IGK1,IGK2)=CQIII*QIII(IGK1,IGK2)
              QIV (IGK1,IGK2)=CQIV *QIV (IGK1,IGK2)
            end do
          end do
        end do
      end do
      RETURN
      END subroutine
C=======================================================================
      subroutine band(igmax,zval,emach,ak,g,al,kapl,kapr,
     &                qi,qii,qiii,qiv)
!!     ------------------------------------------------------------------
!     this subroutine calculates the complex photonic band structure of
!     an infinite crystal. it  provides the  propagating and evanescent
!     eigenmodes of the em field in the given crystal, corresponding to
!     "ak" and a given frequency.
!     ------------------------------------------------------------------
!
! ..  parameter statements ..
!
      integer   igd,igkd,igk2d
      parameter(igd=21,igkd=2*igd,igk2d=2*igkd)

! ..  scalar arguments  ..

      integer    igmax
      real(dp)   zval,emach
      complex(dp) kapl,kapr

! ..  array arguments ..

      real(dp)   ak(2),g(2,igd),al(3)
      complex(dp) qi(igkd,igkd),qii(igkd,igkd)
      complex(dp) qiii(igkd,igkd),qiv(igkd,igkd)

! ..  scalar variables ..

      integer    ii,i,igk1,igk2,igkmax,j
      integer    kd,lib1,lib2,lu,lp,ln,ifail,igk3,igk2m
      real(dp)   aka,bkzre,bkzim
      complex(dp) eaka
!
! ..  array variables ..

      integer    int(igkd)
      real(dp)   ar(igk2d,igk2d),ai(igk2d,igk2d)
      complex(dp) :: a(igk2d,igk2d)
      real(dp)   rr(igk2d),ri(igk2d),vr(igk2d,igk2d),vi(igk2d,igk2d)
      complex(dp) :: r2(igk2d)
      real(dp)   akzap(igk2d),akzip(igk2d)
      real(dp)   akzrep(igk2d),akzimp(igk2d),akzren(igk2d),akzimn(igk2d)
      complex(dp) qh1(igkd,igkd),qh2(igkd,igkd),akz(igk2d)
      complex(dp) comvec(igk2d,igk2d)
      complex(dp) comvec2(igk2d,igk2d)
!     ------------------------------------------------------------------
      igkmax=2*igmax
      igk2m=2*igkmax
      aka=ak(1)*al(1)+ak(2)*al(2)
      eaka=exp(ci*aka)
      do igk1=1,igkmax
        do igk2=1,igkmax
          qh1(igk1,igk2)=czero
          qh2(igk1,igk2)=czero
        end do
      end do
      do igk1=1,igkmax
        qh2(igk1,igk1)=cone
        do igk2=1,igkmax
          do igk3=1,igkmax
            qh1(igk1,igk2)=qh1(igk1,igk2)-qiii(igk1,igk3)*qi (igk3,igk2)
            qh2(igk1,igk2)=qh2(igk1,igk2)-qiii(igk1,igk3)*qii(igk3,igk2)
          end do
        end do
      end do

!     call zgetrf_wrap(qiv, int)
!     do igk2=1,igkmax
!     call zgetrs_wrap(qiv, qh1(:,igk2), int)
!     call zgetrs_wrap(qiv, qh2(:,igk2), int)
!     end do

      call zge(qiv,int,igkmax,igkd,emach)
      do igk2=1,igkmax
        call zsu(qiv,int,qh1(1,igk2),igkmax,igkd,emach)
        call zsu(qiv,int,qh2(1,igk2),igkmax,igkd,emach)
      end do

      do igk1=1,igkmax
        do igk2=1,igkmax
          ar(igk1,igk2)=dble(qi(igk1,igk2))
          ai(igk1,igk2)=aimag(qi(igk1,igk2))
          ar(igk1,igkmax+igk2)=dble(qii(igk1,igk2))
          ai(igk1,igkmax+igk2)=aimag(qii(igk1,igk2))
          ar(igkmax+igk1,igk2)=dble(qh1(igk1,igk2))
          ai(igkmax+igk1,igk2)=aimag(qh1(igk1,igk2))
          ar(igkmax+igk1,igkmax+igk2)=dble(qh2(igk1,igk2))
          ai(igkmax+igk1,igkmax+igk2)=aimag(qh2(igk1,igk2))
        end do
      end do
      do igk1=1,igkmax
        do igk2=1,igkmax
          a(igk1,igk2) = qi(igk1,igk2)
          a(igk1,igkmax+igk2) = qii(igk1,igk2)
          a(igkmax+igk1,igk2) = qh1(igk1,igk2)
          ar(igkmax+igk1,igkmax+igk2)=qh2(igk1,igk2)
        end do
      end do
!     if (.true.) then
      if (.false.) then
        call zgeevx_wrap (a, r2, comvec)
        do ii=1,igk2m
!*****  the if-structure  which follows  can be  omitted  if the accuracy
!*****  'machep' of the subroutine comlr2 is chosen greater than 2**(-47)
!         if((rr(ii)==0.0_dp).and.(ri(ii)==0.0_dp)) then
!           rr(ii)=1.d-20
!           ri(ii)=1.d-20
!         endif
!         ! normalized k_z
          akz(ii)=(-ci/pi)*log(r2(ii)/eaka)
        end do
      else
        call cnaa(igk2d,igk2m,ar,ai,rr,ri,vr,vi,ifail)

        if(ifail/=0) then
          write(6,102) ifail
          stop
        endif
        do ii=1,igk2m
!*****  the if-structure  which follows  can be  omitted  if the accuracy
!*****  'machep' of the subroutine comlr2 is chosen greater than 2**(-47)
          if((rr(ii)==0.0_dp).and.(ri(ii)==0.0_dp)) then
            rr(ii)=1.d-20
            ri(ii)=1.d-20
          endif
          ! normalized k_z
          akz(ii)=(-ci/pi)*log(cmplx(rr(ii),ri(ii),kind=dp)/eaka)
        end do
      endif
      do lib2=1,igk2m
        do lib1=1,igk2m
          comvec(lib1,lib2)=vr(lib1,lib2)+ci*vi(lib1,lib2)
        end do
      end do
      lu=1
      lp=1
      ln=1
      do kd=1,igk2m
!*****warning!! the appropriate limits for aimag(akz(kd))
!*****depend strongly on igmax.
        if(aimag(akz(kd))>0.0_dp) then
          akzrep(lp)=dble(akz(kd))
          akzimp(lp)=aimag(akz(kd))
          lp=lp+1
        else
          akzren(ln)=dble(akz(kd))
          akzimn(ln)=aimag(akz(kd))
          ln=ln+1
        endif
        if(abs(aimag(akz(kd)))>1.0d-2) cycle
        akzap(lu)=dble(akz(kd))
        akzip(lu)=aimag(akz(kd))
        lu=lu+1
      end do

      if (lu<1.1d0) then
        do j=2,lp-1
          bkzim=akzimp(j)
          bkzre=akzrep(j)
          do i=j-1,1,-1
            if(akzimp(i)<=bkzim) go to 15
            akzimp(i+1)=akzimp(i)
            akzrep(i+1)=akzrep(i)
          end do
          i=0
   15     akzimp(i+1)=bkzim
          akzrep(i+1)=bkzre
        end do
        do j=2,ln-1
          bkzim=akzimn(j)
          bkzre=akzren(j)
          do i=j-1,1,-1
            if(akzimn(i)<=bkzim) go to 18
            akzimn(i+1)=akzimn(i)
            akzren(i+1)=akzren(i)
          end do
          i=0
   18     akzimn(i+1)=bkzim
          akzren(i+1)=bkzre
        end do
        write(6,101)  zval,akzrep(1),akzren(ln-1)
        write(6,103)  akzimp(1),akzimn(ln-1)
        write(9,101)  zval,akzrep(1),akzren(ln-1)
        write(9,103)  akzimp(1),akzimn(ln-1)
      else
        write(6,101)  zval,(akzap(i),i=1,lu-1)
        write(9,101)  zval,(akzap(i),i=1,lu-1)
      end if
      return
  101 format(e10.4,3x,10(e10.4,1x))
  102 format(//13x,'error in cnaa   ifail = ',i2)
  103 format(13x,10(e10.4,1x))
      end subroutine
C======================================================================
      END module

      PROGRAM MULTEM
      use libmultem2a
      use libmultem2b
      IMPLICIT NONE
C     ------------------------------------------------------------------
C     A B S T R A C T
C     THIS PROGRAM CALCULATES EITHER THE ABSORBANCE, REFLECTIVITY  AND
C     TRANSMITTANCE  OF   LIGHT  BY  A   FINITE  SLAB   CONSISTING  OF
C     HOMOGENEOUS   PLATES AND   MULTILAYERS  OF  SPHERICAL  PARTICLES
C     ARRANGED IN  A TWO-DIMENSIONAL  BRAVAIS LATTICE, OR THE  COMPLEX
C     PHOTONIC  BAND STRUCTURE OF SUCH AN INFINITE PERIODIC STRUCTURE.
C
C     D E S C R I P T I O N    O F    I N P U T    D A T A
C     KTYPE=     1: THE DIRECTION OF AN INCIDENT  EM WAVE IS SPECIFIED
C                   BY THE POLAR ANGLES OF INCIDENCE "THETA" AND "FI".
C                   THE PROGRAM CALCULATES THE TRANSMISSION,REFLECTION
C                   AND  ABSORPTION   COEFFICIENTS OF  A  FINITE  SLAB
C                2: THE DIRECTION OF  AN INCIDENT EM WAVE IS SPECIFIED
C                   BY THE COMPONENTS  OF THE WAVEVECTOR  PARALLEL  TO
C                   THE  INTERFACES OF THE STRUCTURE:
C                   AQ(1) AND AQ(2) (AND THE  FREQUENCY). THE
C                   PROGRAM  CALCULATES  THE TRANSMISSION, REFLECTION,
C                   ABSORPTION COEFFICIENTS OF A FINITE SLAB
C                3: THE PROGRAM CALCULATES  THE PHOTONIC  COMPLEX BAND
C                   STRUCTURE OF SUCH  AN INFINITE PERIODIC  STRUCTURE
C                   FOR A  WAVEVECTOR WITH COMPONENTS PARALLEL TO  THE
C                   INTERFACES OF THE STRUCTURE: AQ(1) AND AQ(2)
C     KSCAN=     1: SCANNING OVER FREQUENCIES
C                2: SCANNING OVER WAVELENGTHS
C     KEMB        : INDICATES THE PRESENCE (=1) OR ABSENCE (=0) OF A
C                   DIFFERENT EMBEDDING MEDIUM
C     LMAX        : CUTOFF IN SPHERICAL WAVES EXPANSIONS
C     NCOMP       : NUMBER OF DIFFERENT COMPONENTS IN THE UNIT SLICE.
C                   THEIR TYPE IS SPECIFIED  BY THE INTEGER ARRAY
C                   IT(ICOMP)
C     IT=        1: HOMOGENEOUS PLATE OF THICKNESS "D"
C                2: MULTILAYER  OF SPHERICAL  PARTICLES ARRANGED IN  A
C                   2D  BRAVAIS LATTICE.EACH LAYER CONSISTS OF "NPLAN"
C                   NON-PRIMITIVE  PLANES OF SPHERES WITH THE SAME 2-D
C                   PERIODICITY. THE NUMBER OF UNIT LAYERS IS EQUAL TO
C                   2**(NLAYER-1).
C     DL, DR      : POSITION VECTORS INDICATING THE ORIGIN ON THE LEFT
C                   AND ON THE RIGHT OF THE  UNIT,  RESPECTIVELY. BOTH
C                   ARE DIRECTED FROM LEFT TO RIGHT.
C     AL          : PRIMITIVE  TRANSLATION  VECTOR  OF THE  UNIT SLICE
C                   (EFFECTIVE ONLY FOR BAND STRUCTURE CALCULATION).IT
C                   IS GIVEN IN PROGRAM UNITS.
C     NUNIT       : SPECIFIES THE NUMBER OF UNIT SLICES (2**(NUNIT-1))
C                   OF THE SAMPLE
C     ALPHA,ALPHAP: LENGTH OF PRIMITIVE VECTORS OF THE TWO-DIMENSIONAL
C                   LATTICE. IN PROGRAM UNITS THE SIZE OF ALPHA SERVES
C                   AS  THE UNIT LENGTH.  THUS  ALPHA MUST BE EQUAL TO
C                   1.D0
C     FAB         : ANGLE (IN DEG) BETWEEN ALPHA AND ALPHAP
C     RMAX        : UPPER LIMIT FOR THE LENGTH OF  RECIPROCAL  LATTICE
C                   VECTORS (IN UNITS OF 1/ALPHA) WHICH  MUST BE TAKEN
C                   INTO ACCOUNT
C     ZINF,ZSUP   : MINIMUM  AND  MAXIMUM  VALUES  OF  FREQUENCY   (IN
C                   PROGRAM UNITS: OMEGA*ALPHA/C), OR  WAVELENGTH  (IN
C                   PROGRAM UNITS: LAMDA/ALPHA  ),  ACCORDING  TO  THE
C                   VALUE OF KSCAN. C AND LAMDA REFER TO VACUUM
C     NP          : NUMBER OF EQUALLY SPACED POINTS BETWEEN ZINF, ZSUP
C     POLAR       : POLARIZATION ('S ' OR  'P ') OF THE INCIDENT LIGHT
C     AQ(1,2)     : WAVEVECTOR COMPONENTS PARALLEL  TO THE  INTERFACES
C                   OF THE STRUCTURE (XY-PLANE) IN UNITS OF 2*PI/ALPHA
C     THETA,FI    : POLAR ANGLES OF INCIDENCE (IN DEG) OF THE INCIDENT
C                   LIGHT
C     FEIN        : ANGLE  (IN DEG) SPECIFYING  THE DIRECTION  OF  THE
C                   POLARIZATION  VECTOR  FOR  NORMAL  INCIDENCE.  NOT
C                   EFFECTIVE OTHERWISE
C     EPS*,MU*    : RELATIVE DIELECTRIC FUNCTIONS AND MAGNETIC PERMEA-
C                   BILITIES OF THE VARIOUS MEDIA
C     ------------------------------------------------------------------
C
C ..  PARAMETER STATEMENTS ..
C
      INTEGER   LMAXD,IGD,IGKD,NELMD,NCOMPD,NPLAND
      PARAMETER (LMAXD=14,IGD=21,IGKD=2*IGD,
     & NELMD=165152,NCOMPD=8,NPLAND=4)
C
C ..  SCALAR VARIABLES ..
C
      INTEGER      LMAX,I,IGKMAX,IGK1,IGK2,IGMAX,KTYPE,KSCAN,NCOMP,IG1
      INTEGER      N,NP,IG0,NUNIT,ICOMP,KEMB,IU,IPL,ILAYER
      REAL(dp)     ALPHA,EMACH,EPSILON
      REAL(dp)     A0,RA0,RMAX,AKXY
      REAL(dp)     ZVAL,ZSTEP,ZINF,ZSUP,FAB,ALPHAP,THETA,FI,FEIN
      COMPLEX(dp)   KAPPA,KAPPA0,AKZIN,MUEMBL,EPSEMBL
      COMPLEX(dp)   MUEMBR,EPSEMBR,D2,KAPOUT
      COMPLEX(dp)   KAPPAL,KAPPAR,KAPPASL,D1,KAPIN,KAPL,KAPR
      COMPLEX(dp)   MLAST,ELAST,MFIRST,EFIRST,RAP
      CHARACTER(2)  POLAR
      CHARACTER(17) TEXT1(2)
      CHARACTER(5)  DUMMY
C
C ..  ARRAY VARIABLES ..
C
      INTEGER    NT1(IGD),NT2(IGD),IT(NCOMPD)
      INTEGER    NLAYER(NCOMPD),NPLAN(NCOMPD)
      REAL(dp)   ELM(NELMD),AK(2),VECMOD(IGD),DL(3,NCOMPD,NPLAND)
      REAL(dp)   DR(3,NCOMPD,NPLAND),G(2,IGD),B1(2),B2(2)
      REAL(dp)   S(NCOMPD,NPLAND),AL(3),D(NCOMPD),VEC0(3),AQ(2)
      COMPLEX(dp) QIL  (IGKD,IGKD),QIIL(IGKD,IGKD),QIIIL(IGKD,IGKD)
      COMPLEX(dp) QIVL (IGKD,IGKD),QIR (IGKD,IGKD),QIIR (IGKD,IGKD)
      COMPLEX(dp) QIIIR(IGKD,IGKD),QIVR(IGKD,IGKD),WIVL (IGKD,IGKD)
      COMPLEX(dp) WIL  (IGKD,IGKD),WIIL(IGKD,IGKD),WIIIL(IGKD,IGKD)
      COMPLEX(dp) EINCID(IGKD),EIN(2),EPS2(NCOMPD),EPS3(NCOMPD)
      COMPLEX(dp) MU1(NCOMPD),MU2(NCOMPD),MU3(NCOMPD),EPS1(NCOMPD)
      COMPLEX(dp) MUSPH(NCOMPD,NPLAND), EPSSPH(NCOMPD,NPLAND)
C
C ..  COMMON BLOCKS ..
C
      REAL(dp)   AR1(2),AR2(2)
      COMMON/X1/AR1,AR2
C
C ..  DATA STATEMENTS ..
C
      DATA EMACH/1.D-8/,EPSILON/0.D0/
      DATA EINCID/IGKD*(0.D0,0.D0)/,VEC0/3*0.D0/
      DATA TEXT1/'HOMOGENEOUS PLATE','PHOTONIC CRYSTAL'/
C     ------------------------------------------------------------------
C
      READ(10,200) KTYPE,KSCAN,KEMB,LMAX,NCOMP,NUNIT
      IF(KTYPE<=0.OR.KTYPE>=4) STOP 'ILLEGAL INPUT VALUE OF KTYPE'
      IF(KSCAN<=0.OR.KSCAN>=3) STOP 'ILLEGAL INPUT VALUE OF KSCAN'
      IF(KEMB<0.OR.KEMB>=2)   STOP 'ILLEGAL INPUT VALUE OF KEMB '
      IF(LMAX<=0.OR.LMAX>LMAXD.OR.LMAX>14)
     &          STOP 'LMAX<=0.OR.LMAX>MIN0(14,LMAXD)'
      IF(NCOMP<=0.OR.NCOMP>NCOMPD)
     &                   STOP 'ILLEGAL INPUT VALUE OF NCOMP'
      IF(NUNIT<=0)           STOP 'ILLEGAL INPUT VALUE OF NUNIT'
      READ(10,202) ALPHA,ALPHAP,FAB,RMAX
      FAB=FAB*PI/180.D0
      READ(10,203) NP,ZINF,ZSUP
      IF(NP<=1)                  STOP 'ILLEGAL INPUT VALUE OF  NP '
              IF(KTYPE>=2) THEN
      READ(10,204) AQ(1),AQ(2),POLAR,FEIN
      FEIN=FEIN*PI/180.D0
      AQ(1)=2.D0*PI*AQ(1)
      AQ(2)=2.D0*PI*AQ(2)
      IF(KTYPE<3) THEN
      WRITE(6,222)
      ELSE
      WRITE(6,223)
      ENDIF
      IF(KTYPE==2) WRITE(6,207) AQ(1),AQ(2),POLAR
      IF(KTYPE==3) WRITE(6,225) AQ(1),AQ(2)
                             ELSE
      READ(10,204) THETA,FI,POLAR,FEIN
      WRITE(6,208) THETA,FI,POLAR
      FEIN=FEIN*PI/180.D0
      THETA=THETA*PI/180.D0
      FI=FI*PI/180.D0
                             ENDIF
      DO 3 ICOMP=1,NCOMP
      READ(10,201) IT(ICOMP)
      IF(IT(ICOMP)<=0.OR.IT(ICOMP)>2)
     &                    STOP 'ILLEGAL COMPONENT TYPE'
      WRITE(6,209) ICOMP,TEXT1(IT(ICOMP))
      IF(IT(ICOMP)==1) THEN
      READ(10,204) D(ICOMP)
      READ(10,205) MU1(ICOMP),EPS1(ICOMP),MU2(ICOMP),EPS2(ICOMP),
     &             MU3(ICOMP),EPS3(ICOMP)
      WRITE(6,210) MU1(ICOMP),MU2(ICOMP),MU3(ICOMP),EPS1(ICOMP),
     &             EPS2(ICOMP),EPS3(ICOMP)
      READ(10,*) DUMMY,(DL(I,ICOMP,1),I=1,3)
      READ(10,*) DUMMY,(DR(I,ICOMP,1),I=1,3)
                      ELSE
      READ(10,205) MU1(ICOMP),EPS1(ICOMP)
      IF(dble(MU1(ICOMP))<=0.D0.OR.dble(EPS1(ICOMP))<=0.D0)
     &THEN
      WRITE(6,226)
      STOP
      ENDIF
      READ(10,201) NPLAN(ICOMP),NLAYER(ICOMP)
      DO IPL=1,NPLAN(ICOMP)
        READ(10,206) S(ICOMP,IPL),MUSPH(ICOMP,IPL),EPSSPH(ICOMP,IPL)
        READ(10,*) DUMMY,(DL(I,ICOMP,IPL),I=1,3)
        READ(10,*) DUMMY,(DR(I,ICOMP,IPL),I=1,3)
      end do
      WRITE(6,211)  MU1(ICOMP),( MUSPH(ICOMP,IPL),IPL=1,NPLAN(ICOMP))
      WRITE(6,220) EPS1(ICOMP),(EPSSPH(ICOMP,IPL),IPL=1,NPLAN(ICOMP))
      WRITE(6,224) (S(ICOMP,IPL),IPL=1,NPLAN(ICOMP))
      WRITE(6,212) 2**(NLAYER(ICOMP)-1)
              ENDIF
    3 CONTINUE
                         D1=SQRT(MU1(1)    *EPS1(1))
                         D2=SQRT(MU1(NCOMP)*EPS1(NCOMP))
      IF(IT(NCOMP)==1) D2=SQRT(MU3(NCOMP)*EPS3(NCOMP))
      IF(aimag(D1)/=0.D0) THEN
      WRITE(6,227)
      STOP
      ENDIF
      IF(aimag(D2)/=0.D0) THEN
      WRITE(6,228)
      STOP
      ENDIF
      IF(KTYPE/=3) THEN
      WRITE(6,221) 2**(NUNIT-1)
      IF(KEMB==1) THEN
      READ(10,205) MUEMBL,EPSEMBL
      READ(10,205) MUEMBR,EPSEMBR
      D1=SQRT(MUEMBL*EPSEMBL)
      D2=SQRT(MUEMBR*EPSEMBR)
      IF(aimag(D1)/=0.D0) THEN
      WRITE(6,227)
      STOP
      ENDIF
      IF(aimag(D2)/=0.D0) THEN
      WRITE(6,228)
      STOP
      ENDIF
            ENDIF
             ELSE
      READ(10,*) DUMMY,(AL(I),I=1,3)
             ENDIF
      CALL ELMGEN(ELM,NELMD,LMAX)
C
C****** DEFINE THE 2D DIRECT AND RECIPROCAL-LATTICE VECTORS ******
C
      AR1(1)= ALPHA
      AR1(2)= 0.D0
      AR2(1)= ALPHAP*COS(FAB)
      AR2(2)= ALPHAP*SIN(FAB)
      WRITE(6,213) AR1(1),AR1(2),AR2(1),AR2(2)
      A0=ABS(AR1(1)*AR2(2)-AR1(2)*AR2(1))
      RA0=2.D0*PI/A0
      B1(1)=-AR1(2)*RA0
      B1(2)= AR1(1)*RA0
      B2(1)=-AR2(2)*RA0
      B2(2)= AR2(1)*RA0
      CALL LAT2D(B1,B2,RMAX,IGMAX,IGD,NT1,NT2,VECMOD)
      WRITE(6,214) B1(1),B1(2),B2(1),B2(2)
      IGKMAX=2*IGMAX
      DO IG1=1,IGMAX
        G(1,IG1)=NT1(IG1)*B1(1)+NT2(IG1)*B2(1)
        G(2,IG1)=NT1(IG1)*B1(2)+NT2(IG1)*B2(2)
        WRITE(6,215) IG1,NT1(IG1),NT2(IG1),VECMOD(IG1)
      end do
      ZSTEP=(ZSUP-ZINF)/dble(NP-1)
      ZVAL=ZINF-ZSTEP
      IF(KTYPE<3) THEN
      IF(KSCAN==1) WRITE(6,216)
      IF(KSCAN==2) WRITE(6,217)
      IF(POLAR/='S '.AND.POLAR/='P ') STOP 'ILLEGAL POLARIZATION'
                     ELSE
      IF(KSCAN==1) WRITE(6,218)
      IF(KSCAN==2) WRITE(6,219)
                     ENDIF
      IF(POLAR=='P ') THEN
      EIN(1)=CONE
      EIN(2)=CZERO
                        ELSE
      EIN(1)=CZERO
      EIN(2)=CONE
                        END IF
      DO 1 N=1,NP   !****** SCANNING OVER FREQUENCIES/WAVELENGTHS ******
      ZVAL=ZVAL+ZSTEP
      IF(KSCAN==1) KAPPA0=cmplx_dp(ZVAL,EPSILON)
      IF(KSCAN==2) KAPPA0=cmplx_dp(2.D0*PI/ZVAL,EPSILON)
      KAPIN =KAPPA0*D1
      KAPOUT=KAPPA0*D2
                                             IF(KTYPE==1) THEN
                           AK(1)=dble(KAPIN)*SIN(THETA)*COS(FI)
                           AK(2)=dble(KAPIN)*SIN(THETA)*SIN(FI)
               DO 50 I=1,IGKMAX
               EINCID(I)=CZERO
  50                       CONTINUE
                                ELSE
                           AK(1)=AQ(1)
               AK(2)=AQ(2)
                                                           ENDIF
      IF(KTYPE/=3) THEN !DEFINE THE POLARIZATION VECTOR FROM "AK"*****
      AKXY=AK(1)*AK(1)+AK(2)*AK(2)
      AKZIN=SQRT(KAPIN*KAPIN-AKXY)
      IF(dble(AKZIN)<EMACH)      STOP 'IMPROPER INCIDENT WAVE'
      AKXY=SQRT(AKXY)
      IF(AKXY<EMACH) THEN
      EIN(1)=cmplx_dp(COS(FEIN),0.D0)
      EIN(2)=cmplx_dp(SIN(FEIN),0.D0)
                        END IF
      CALL REDUCE(AR1,AR2,AK,IGMAX,G,IG0,EMACH)   !"AK" IN SBZ*******
      DO 2 I=1,2
      EINCID(2*IG0-2+I)=EIN(I)
    2 CONTINUE
                     ELSE
      CALL REDUCE(AR1,AR2,AK,IGMAX,G,IG0,EMACH)   !"AK" IN SBZ*******
             ENDIF
C
C****** CONSTRUCT THE TRANSFER MATRIX OF THE UNIT SLICE ******
C
      IF(IT(1)==1) THEN
      KAPPAL =SQRT(MU1(1)*EPS1(1))*KAPPA0
      KAPPASL=SQRT(MU2(1)*EPS2(1))*KAPPA0
      KAPPAR =SQRT(MU3(1)*EPS3(1))*KAPPA0
      KAPL=KAPPAL
      KAPR=KAPPAR
      MLAST=MU3(1)
      ELAST=EPS3(1)
      MFIRST=MU1(1)
      EFIRST=EPS1(1)
      CALL HOSLAB(IGMAX,KAPPAL,KAPPASL,KAPPAR,AK,G,DL(1,1,1),
     &           DR(1,1,1),D(1),QIL,QIIL,QIIIL,QIVL,EMACH)
                     ELSE
      KAPPA=SQRT(MU1(1)*EPS1(1))*KAPPA0
      KAPL=KAPPA
      KAPR=KAPPA
      MLAST=MU1(1)
      ELAST=EPS1(1)
      MFIRST=MU1(1)
      EFIRST=EPS1(1)
      RAP=S(1,1)*KAPPA0/2.D0/PI
      CALL PCSLAB(LMAX,IGMAX,RAP,EPS1(1),EPSSPH(1,1),MU1(1),MUSPH(1,1),
     &           KAPPA,AK,DL(1,1,1),DR(1,1,1),G,ELM,A0,EMACH,
     &         QIL,QIIL,QIIIL,QIVL)
      IF(NPLAN(1)>=2) THEN
      DO 13 IPL=2,NPLAN(1)
      RAP=S(1,IPL)*KAPPA0/2.D0/PI
      CALL PCSLAB(LMAX,IGMAX,RAP,EPS1(1),EPSSPH(1,IPL),MU1(1),
     &           MUSPH(1,IPL),KAPPA,AK,DL(1,1,IPL),DR(1,1,IPL),
     &           G,ELM,A0,EMACH, QIR,QIIR,QIIIR,QIVR)
      CALL PAIR(IGKMAX,igkd,QIL,QIIL,QIIIL,QIVL,QIR,QIIR,QIIIR,QIVR)
   13 CONTINUE
            ENDIF
      IF(NLAYER(1)>=2) THEN
      DO 14 ILAYER=1,NLAYER(1)-1
      CALL PAIR(IGKMAX,igkd,QIL,QIIL,QIIIL,QIVL,QIL,QIIL,QIIIL,QIVL)
   14 CONTINUE
             ENDIF
                      ENDIF
                                                    IF(NCOMP>=2) THEN
      DO 4 ICOMP=2,NCOMP
      IF(IT(ICOMP)==1) THEN
      KAPPAL =SQRT(MU1(ICOMP)*EPS1(ICOMP))*KAPPA0
      KAPPASL=SQRT(MU2(ICOMP)*EPS2(ICOMP))*KAPPA0
      KAPPAR =SQRT(MU3(ICOMP)*EPS3(ICOMP))*KAPPA0
      KAPR=KAPPAR
      IF(ABS(MU1(ICOMP)-MLAST)/=0.D0.OR.ABS(EPS1(ICOMP)-ELAST)/=
     &        0.D0) STOP 'IMPROPER MATCHING OF SUCCESSIVE HOST MEDIA'
      MLAST=MU3(ICOMP)
      ELAST=EPS3(ICOMP)
      CALL HOSLAB(IGMAX,KAPPAL,KAPPASL,KAPPAR,AK,G,DL(1,ICOMP,1),
     &            DR(1,ICOMP,1),D(ICOMP),QIR,QIIR,QIIIR,QIVR,EMACH)
                      ELSE
      KAPPA=SQRT(MU1(ICOMP)*EPS1(ICOMP))*KAPPA0
      KAPR=KAPPA
      IF(ABS(MU1(ICOMP)-MLAST)/=0.D0.OR.ABS(EPS1(ICOMP)-ELAST)/=
     &        0.D0) STOP 'IMPROPER MATCHING OF SUCCESSIVE HOST MEDIA'
      MLAST=MU1(ICOMP)
      ELAST=EPS1(ICOMP)
      RAP=S(ICOMP,1)*KAPPA0/2.D0/PI
      CALL PCSLAB(LMAX,IGMAX,RAP,EPS1(ICOMP),EPSSPH(ICOMP,1),MU1(ICOMP)
     &           ,MUSPH(ICOMP,1),KAPPA,AK,DL(1,ICOMP,1),
     &           DR(1,ICOMP,1),G,ELM,A0,EMACH,QIR,QIIR,QIIIR,QIVR)
      IF(NPLAN(ICOMP)>=2) THEN
         DO IGK1=1,IGKMAX
           DO IGK2=1,IGKMAX
             WIL  (IGK1,IGK2)=QIR  (IGK1,IGK2)
             WIIL (IGK1,IGK2)=QIIR (IGK1,IGK2)
             WIIIL(IGK1,IGK2)=QIIIR(IGK1,IGK2)
             WIVL (IGK1,IGK2)=QIVR (IGK1,IGK2)
           end do
         end do
      DO IPL=2,NPLAN(ICOMP)
        RAP=S(ICOMP,IPL)*KAPPA0/2.D0/PI
        CALL PCSLAB(LMAX,IGMAX,RAP,EPS1(ICOMP),EPSSPH(ICOMP,IPL),
     &             MU1(ICOMP),MUSPH(ICOMP,IPL),KAPPA,AK,
     &             DL(1,ICOMP,IPL),DR(1,ICOMP,IPL),G,ELM,A0,EMACH,
     &             QIR,QIIR,QIIIR,QIVR)
        CALL PAIR(IGKMAX,igkd,WIL,WIIL,WIIIL,WIVL,QIR,QIIR,QIIIR,QIVR)
      end do
         DO IGK1=1,IGKMAX
           DO IGK2=1,IGKMAX
             QIR  (IGK1,IGK2)=WIL  (IGK1,IGK2)
             QIIR (IGK1,IGK2)=WIIL (IGK1,IGK2)
             QIIIR(IGK1,IGK2)=WIIIL(IGK1,IGK2)
             QIVR (IGK1,IGK2)=WIVL (IGK1,IGK2)
           end do
         end do
            ENDIF
      IF(NLAYER(ICOMP)>=2) THEN
      DO ILAYER=1,NLAYER(ICOMP)-1
        CALL PAIR(IGKMAX,igkd,QIR,QIIR,QIIIR,QIVR,QIR,QIIR,QIIIR,QIVR)
      end do
             ENDIF
                      ENDIF
      CALL PAIR(IGKMAX,igkd,QIL,QIIL,QIIIL,QIVL,QIR,QIIR,QIIIR,QIVR)
    4 CONTINUE
                                                                   ENDIF
      IF(KTYPE<3) THEN
C
C****** THE UNIT SLICE IS DEFINED. THIS CAN BE REPEATED BY THE ******
C****** DOUBLING-LAYER  TECHNIQUE, INTERFACES CAN BE ADDED AND ******
C****** REFLECTIVITY/TRANSMITTANCE/ABSORBANCE ARE CALCULATED.  ******
C
             IF(NUNIT==1) GO TO 30
         IF(ABS(MLAST-MFIRST)/=0.D0.OR.ABS(ELAST-EFIRST)/=0.D0)
     &       STOP 'IMPROPER MATCHING OF SUCCESSIVE HOST MEDIA'
         DO IU=1,NUNIT-1
             CALL PAIR(IGKMAX,igkd,QIL,QIIL,QIIIL,QIVL,
     &                 QIL,QIIL,QIIIL,QIVL)
         end do
   30        CONTINUE
             IF(KEMB==1) THEN
             CALL HOSLAB(IGMAX,KAPR,(KAPR+KAPOUT)/2.D0,KAPOUT,AK,G,VEC0,
     &                   VEC0,0.D0,QIR,QIIR,QIIIR,QIVR,EMACH)
         CALL PAIR(IGKMAX,igkd,QIL,QIIL,QIIIL,QIVL,QIR,QIIR,QIIIR,QIVR)
         DO IGK1=1,IGKMAX
           DO IGK2=1,IGKMAX
             QIR  (IGK1,IGK2)=QIL  (IGK1,IGK2)
             QIIR (IGK1,IGK2)=QIIL (IGK1,IGK2)
             QIIIR(IGK1,IGK2)=QIIIL(IGK1,IGK2)
             QIVR (IGK1,IGK2)=QIVL (IGK1,IGK2)
           end do
         end do
             CALL HOSLAB(IGMAX,KAPIN,(KAPL+KAPIN)/2.D0,KAPL,AK,G,VEC0,
     &                   VEC0,0.D0,QIL,QIIL,QIIIL,QIVL,EMACH)
         CALL PAIR(IGKMAX,igkd,QIL,QIIL,QIIIL,QIVL,QIR,QIIR,QIIIR,QIVR)
                      ENDIF
             CALL SCAT(IGMAX,ZVAL,AK,G,dble(KAPIN),dble(KAPOUT),
     &                 EINCID,QIL,QIIIL)
                     ELSE
C
C****** ALTERNATIVELY, CALCULATE COMPLEX PHOTONIC BAND STRUCTURE ******
C
      IF(ABS(MLAST-MFIRST)/=0.D0.OR.ABS(ELAST-EFIRST)/=0.D0)
     &    STOP 'IMPROPER MATCHING OF SUCCESSIVE HOST MEDIA'
      CALL BAND(IGMAX,ZVAL,EMACH,AK,G,AL,KAPL,KAPR,QIL,QIIL,QIIIL,QIVL)
                     ENDIF
    1 CONTINUE
      STOP
  200 FORMAT(///,6(10X,I2))
  201 FORMAT(6(10X,I2))
  202 FORMAT(4(8X,F12.6))
  203 FORMAT(6X,I4,2(8X,F13.8))
  204 FORMAT(2(15X,F13.8),10X,A2,10X,F7.2///)
  205 FORMAT(2(12X,2F13.8))
  206 FORMAT(10X,F13.8,2(12X,2F13.8))
  207 FORMAT(3X,'K_PARALLEL=',2F12.6,5X,A2,'POLARIZATION')
  208 FORMAT(3X,'ANGLES OF INCIDENCE (IN RAD):  THETA=',F7.2,3X,'FI=',
     &       F7.2,5X,A2,'POLARIZATION')
  209 FORMAT(3X,'COMPONENT NR.',I2,3X,'TYPE:',2X,A17)
  210 FORMAT(3X,'MU :',2F10.5,' | ',2F10.5,' | ',2F10.5/
     &       3X,'EPS:',2F10.5,' | ',2F10.5,' | ',2F10.5)
  211 FORMAT(3X,'MU :',2F10.5,' | ',4(3X,2F10.5))
  212 FORMAT(29X,I6,' UNIT LAYERS')
  213 FORMAT(/13X,'PRIMITIVE LATTICE VECTORS'/13X,'AR1 = (',2F12.4,')'/
     &        13X,'AR2 = (',2F12.4,')')
  214 FORMAT(13X,'UNIT VECTORS IN RECIPROCAL SPACE:'/13X,'B1  = (',
     &2F12.4, ')'/13X,'B2  = (',2F12.4,')'//3X,'RECIPROCAL VECTORS',5X,
     &       'LENGTH')
  215 FORMAT(I3,4X,2I5,5X,E14.6)
  216 FORMAT(//4X,'FREQUENCY   TRANSMITTANCE  REFLECTANCE   ABSORBANCE'
     &        /60('-'))
  217 FORMAT(//4X,'WAVELENGTH  TRANSMITTANCE  REFLECTANCE   ABSORBANCE'
     &        /60('-'))
  218 FORMAT(//1X,'FREQUENCY',7X,'NORMALIZED K_Z'/
     &         1X,9('-'),7X,14('-'))
  219 FORMAT(//3X,'WAVELENGTH VERSUS NORMALIZED K_Z'/3X,32('-'))
  220 FORMAT(3X,'EPS:',2F10.5,' | ',4(3X,2F10.5))
  221 FORMAT(3X,'THE SAMPLE CONSISTS OF ',I6,' UNIT SLICES')
  222 FORMAT(5X,'****************************************************'/
     &       5X,'*** OUTPUT: TRANSMITTANCE/REFLECTANCE/ABSORBANCE ***'/
     &       5X,'****************************************************')
  223 FORMAT(5X,'****************************************************'/
     &       5X,'************** OUTPUT: BAND STRUCTURE **************'/
     &       5X,'****************************************************')
  224 FORMAT(3X,'  S:',23X,4(3X,F10.5,10X))
  225 FORMAT(3X,'K_PARALLEL=',2F12.6)
  226 FORMAT(5X,'----------------------------------'/
     &       5X,'ILLEGAL INPUT: SPHERES EMBEDDED IN'/
     &       5X,'A  MEDIUM  OF NEGATIVE  DIELECTRIC'/
     &       5X,'CONSTANT. THE EWALD  SUMMATION  IN'/
     &       5X,'SUBROUTINE XMAT DOES NOT CONVERGE.'/
     &       5X,'DIRECT - SPACE SUMMATION IS NEEDED'/
     &       5X,'INSTEAD.'/
     &       5X,'----------------------------------')
  227 FORMAT(5X,'----------------------------------'/
     &       5X,'ILLEGAL INPUT:SEMI-INFINITE MEDIUM'/
     &       5X,'OF COMPLEX REFRACTIVE INDEX ON THE'/
     &       5X,'LEFT SIDE OF THE SLAB.'/
     &       5X,'----------------------------------')
  228 FORMAT(5X,'----------------------------------'/
     &       5X,'ILLEGAL INPUT:SEMI-INFINITE MEDIUM'/
     &       5X,'OF COMPLEX REFRACTIVE INDEX ON THE'/
     &       5X,'RIGHT SIDE OF THE SLAB.'/
     &       5X,'----------------------------------')
      end program
