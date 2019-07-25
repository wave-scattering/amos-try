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

C=======================================================================
C======================================================================
C=======================================================================
      SUBROUTINE DLMKG(LMAX,A0,GK,SIGNUS,KAPPA,DLME,DLMH,EMACH)
C     ------------------------------------------------------------------
C     THIS SUBROUTINE CALCULATES THE COEFFICIENTS DLM(KG)
C     ------------------------------------------------------------------

C ..  ARGUMENTS  ..
      INTEGER    LMAX
      REAL(dp)   A0,SIGNUS,EMACH
      COMPLEX(dp) KAPPA
      COMPLEX(dp) DLME(2,(lmax+1)**2),DLMH(2,(lmax+1)**2),GK(3)
C  .. LOCAL
      INTEGER    K,II,L,M,I
      REAL(dp)   AKPAR,ALPHA,BETA,AKG1,AKG2
      COMPLEX(dp) C0,CC,COEF,Z1,Z2,Z3
      COMPLEX(dp) CT,ST,CF
      COMPLEX(dp) YLM((lmax+1)**2)
C     ------------------------------------------------------------------
      AKG1=dble(GK(1))
      AKG2=dble(GK(2))
      DO K=1,2
        DLME(K,1)=CZERO
        DLMH(K,1)=CZERO
      end do
      IF(ABS(GK(3))<EMACH)   THEN
      WRITE(7,101)
      STOP
      ENDIF
      C0=2.0_dp*PI/(KAPPA*A0*GK(3)*SIGNUS)
      AKPAR=SQRT(AKG1*AKG1+AKG2*AKG2)
      CT=GK(3)/KAPPA
      ST=AKPAR/KAPPA
      CF=CONE
      IF(AKPAR>1.D-8) CF=cmplx_dp(AKG1/AKPAR,AKG2/AKPAR)
      CALL SPHRM4(YLM,CT,ST,CF,LMAX)
      II=1
      CC=CONE
      DO L=1,LMAX
        CC=CC/CI
        COEF=C0*CC/SQRT(dble(L*(L+1)))
        DO M=-L,L
          II=II+1
          ALPHA=SQRT(dble((L-M)*(L+M+1)))/2.0_dp
          BETA =SQRT(dble((L+M)*(L-M+1)))/2.0_dp
          if(ABS(M+1)<=L)  then
            I=L*L+L+M+2
            Z1=YLM(I)
          else
            Z1=CZERO
          end if
          if(ABS(M-1)<=L)  then
            I=L*L+L+M
            Z2=YLM(I)
          else
            Z2=CZERO
          end if
          I=L*L+L+M+1
          Z3=YLM(I)
          DLMH(1,II)=COEF*(BETA*CT*CF*Z2-dble(M)*ST*Z3
     &             +ALPHA*CT*CONJG(CF)*Z1)
          DLMH(2,II)=COEF*CI*(BETA*CF*Z2-ALPHA*CONJG(CF)*Z1)
          DLME(1,II)=COEF*CI*(BETA*CF*Z2-ALPHA*CONJG(CF)*Z1)
          DLME(2,II)=-COEF*(BETA*CT*CF*Z2-dble(M)*ST*Z3
     &             +ALPHA*CT*CONJG(CF)*Z1)
        end do
      end do
      RETURN
  101 FORMAT(13X,'FATAL ERROR FROM DLMKG:'/3X,'GK(3) IS TOO SMALL.'
     & /3X,'GIVE A SMALL BUT NONZERO VALUE FOR "EPSILON"'/3X,
     & 'IN THE DATA STATEMENT OF THE MAIN PROGRAM.'
     & /3X,'THIS DEFINES A SMALL IMAGINARY PART'
     & /3X,'IN THE FREQUENCY OR WAVELENGTH VALUE.')
      END subroutine
C=======================================================================
C=======================================================================
C=======================================================================
      SUBROUTINE PLW(KAPPA,GK,LMAX,AE,AH)

C     ------------------------------------------------------------------
C     THIS ROUTINE CALCULATES THE EXPANSION COEFFICIENTS 'AE,AH' OF AN
C     INCIDENT PLANE ELECTROMAGNETIC WAVE OF WAVE VECTOR  'KAPPA' WITH
C     COMPONENTS PARALLEL TO THE SURFACE EQUAL TO   '(GK(1),GK(2))'.
C     ------------------------------------------------------------------
C ..  SCALAR ARGUMENTS  ..
C
      INTEGER    LMAX
      COMPLEX(dp) KAPPA
C
C ..  ARRAY ARGUMENTS  ..
C
      COMPLEX(dp) AE(:,:),AH(:,:),GK(:)
C
C ..  LOCAL SCALARS  ..
C
      INTEGER    M,II,L,I,K
      REAL(dp)   AKPAR,FPI,A,SIGNUS,AKG1,AKG2
      COMPLEX(dp) CT,ST,CF,N1,N2,N3,CC,CC1,Z1,Z2,Z3
C
C ..  LOCAL ARRAYS  ..
C
      COMPLEX(dp) YLM((lmax+1)**2)
C-----------------------------------------------------------------------
C
      AKG1=dble(GK(1))
      AKG2=dble(GK(2))
      DO 3 K=1,2
      AE(K,1)=CZERO
    3 AH(K,1)=CZERO
      FPI=4.0_dp*PI
      AKPAR=SQRT(AKG1*AKG1+AKG2*AKG2)
      CT=GK(3)/KAPPA
      ST=AKPAR/KAPPA
      CF=CONE
      IF(AKPAR>1.D-8) CF=CMPLX(AKG1/AKPAR,AKG2/AKPAR, kind=dp)
      N1=AKG1/KAPPA
      N2=AKG2/KAPPA
      N3=GK(3)/KAPPA
      CALL SPHRM4(YLM,CT,ST,CF,LMAX)
      II=1
      CC=cmplx_dp(FPI,0.0_dp)
      SIGNUS=-1.0_dp
      DO 1 L=1,LMAX
      CC=CC*CI
      A=dble(L*(L+1))
      CC1=CC/SQRT(A)
      DO 2 M=-L,L
      SIGNUS=-SIGNUS
      II=II+1
      IF(ABS(M+1)<=L)  then
      I=L*L+L-M
      Z1=CC1*SQRT(dble((L-M)*(L+M+1)))*YLM(I)/2.0_dp
                         else
      Z1=CZERO
                         end if
      IF(ABS(M-1)<=L)  then
      I=L*L+L-M+2
      Z2=CC1*SQRT(dble((L+M)*(L-M+1)))*YLM(I)/2.0_dp
                         else
      Z2=CZERO
                         end if
      I=L*L+L-M+1
      Z3=CC1*dble(M)*YLM(I)
      AE(1,II)= SIGNUS*CI*(CF*Z1-CONJG(CF)*Z2)
      AE(2,II)=-SIGNUS*(CT*CF*Z1+ST*Z3+CT*CONJG(CF)*Z2)
      AH(1,II)= SIGNUS*(CT*CF*Z1+ST*Z3+CT*CONJG(CF)*Z2)
      AH(2,II)= SIGNUS*CI*(CF*Z1-CONJG(CF)*Z2)
    2 CONTINUE
    1 CONTINUE
      RETURN
      END subroutine
C=======================================================================
      SUBROUTINE SETUP(LMAX,XEVEN,XODD,TE,TH,XXMAT1,XXMAT2)

C     ------------------------------------------------------------------
C     THIS SUBROUTINE CONSTRUCTS THE SECULAR MATRIX
C     ------------------------------------------------------------------
C
C ..  SCALAR ARGUMENTS ..
C
      INTEGER LMAX
C
C ..  ARRAY ARGUMENTS ..
C
      COMPLEX(dp) XEVEN(:,:),XODD(:,:)
      COMPLEX(dp) XXMAT2(:,:)
      COMPLEX(dp) TE(:),TH(:),XXMAT1(:,:)
C
C ..  LOCAL SCALARS ..
C
      INTEGER IA,LA,MA,LMTOT,LTT,LMAX1,IB,LB,MB,I,LMXOD,IAOD,IAEV,IBOD
      INTEGER IBEV
      REAL(dp)  C0,SIGNUS,UP,C,B1,B2,B3,U1,U2,A,DOWN
      REAL(dp)  ALPHA1,ALPHA2,BETA1,BETA2
      COMPLEX(dp) OMEGA1,OMEGA2,Z1,Z2,Z3
C     ------------------------------------------------------------------
      LMAX1=LMAX+1
      LMTOT=LMAX1*LMAX1-1
      LMXOD=(LMAX*LMAX1)/2
      C0=SQRT(8.0_dp*PI/3.0_dp)
      SIGNUS=1.0_dp
      IAOD=0
      IAEV=LMXOD
      DO 1 LA=1,LMAX
      DO 1 MA=-LA,LA
      IF(MOD((LA+MA),2)==0) THEN
      IAEV=IAEV+1
      IA=IAEV
                              ELSE
      IAOD=IAOD+1
      IA=IAOD
                              END IF
      UP=dble(2*LA+1)
      SIGNUS=-SIGNUS
      C=SIGNUS*C0
      B1=0.0_dp
      IF(ABS(MA+1)<=(LA-1)) B1=BLM(LA-1,MA+1,1,-1,LA,-MA,LMAX)
      B2=0.0_dp
      IF(ABS(MA-1)<=(LA-1)) B2=BLM(LA-1,MA-1,1, 1,LA,-MA,LMAX)
      U1=dble((LA+MA)*(LA-MA))
      U2=dble((2*LA-1)*(2*LA+1))
      B3=SQRT(U1/U2)
      ALPHA1=SQRT(dble((LA-MA)*(LA+MA+1)))/2.0_dp
      BETA1 =SQRT(dble((LA+MA)*(LA-MA+1)))/2.0_dp
      IBOD=0
      IBEV=LMXOD
      DO 2 LB=1,LMAX
      DO 2 MB=-LB,LB
      IF(MOD((LB+MB),2)==0) THEN
      IBEV=IBEV+1
      IB=IBEV
                              ELSE
      IBOD=IBOD+1
      IB=IBOD
                              END IF
      A=dble(LB*(LB+1)*LA*(LA+1))
      DOWN=SQRT(A)
      ALPHA2=SQRT(dble((LB-MB)*(LB+MB+1)))/2.0_dp
      BETA2 =SQRT(dble((LB+MB)*(LB-MB+1)))/2.0_dp
      LTT=LA+MA+LB+MB
          IF(MOD(LTT,2)/=0)           then
             IF(MOD((LA+MA),2)==0)       then
             Z1=CEVEN(LB,MB+1,LA-1,MA+1,XEVEN)
             Z2=CEVEN(LB,MB-1,LA-1,MA-1,XEVEN)
             Z3=CODD (LB,MB  ,LA-1,MA  ,XODD )
             Z1= C*ALPHA2*B1*Z1
             Z2=-C*BETA2* B2*Z2
             Z3=dble(MB)*B3*Z3
             OMEGA2=UP*(Z1+Z2+Z3)/DOWN
             XXMAT1(IA,IB)=-TH(LA+1)*OMEGA2
             XXMAT2(IA,IB)= TE(LA+1)*OMEGA2
                                           else
             Z1=CODD (LB,MB+1,LA-1,MA+1,XODD )
             Z2=CODD (LB,MB-1,LA-1,MA-1,XODD )
             Z3=CEVEN(LB,MB  ,LA-1,MA  ,XEVEN)
             Z1= C*ALPHA2*B1*Z1
             Z2=-C*BETA2* B2*Z2
             Z3=dble(MB)*B3*Z3
             OMEGA2=UP*(Z1+Z2+Z3)/DOWN
             XXMAT1(IA,IB)= TE(LA+1)*OMEGA2
             XXMAT2(IA,IB)=-TH(LA+1)*OMEGA2
                                           end if
                                        else
             IF(MOD((LA+MA),2)==0)       then
             Z1=CODD (LB,MB-1,LA,MA-1,XODD )
             Z2=CODD (LB,MB+1,LA,MA+1,XODD )
             Z3=CEVEN(LB,MB  ,LA,MA  ,XEVEN)
             Z1=2.0_dp*BETA1 *BETA2 *Z1
             Z2=2.0_dp*ALPHA1*ALPHA2*Z2
             Z3=dble(MA)*dble(MB)*Z3
             OMEGA1=(Z1+Z2+Z3)/DOWN
             XXMAT1(IA,IB)=-TH(LA+1)*OMEGA1
             XXMAT2(IA,IB)=-TE(LA+1)*OMEGA1
                                           else
             Z1=CEVEN(LB,MB-1,LA,MA-1,XEVEN)
             Z2=CEVEN(LB,MB+1,LA,MA+1,XEVEN)
             Z3=CODD (LB,MB  ,LA,MA  ,XODD )
             Z1=2.0_dp*BETA1 *BETA2 *Z1
             Z2=2.0_dp*ALPHA1*ALPHA2*Z2
             Z3=dble(MA)*dble(MB)*Z3
             OMEGA1=(Z1+Z2+Z3)/DOWN
             XXMAT1(IA,IB)=-TE(LA+1)*OMEGA1
             XXMAT2(IA,IB)=-TH(LA+1)*OMEGA1
                                           end if
                                        end if
    2 CONTINUE
    1 CONTINUE
      DO 3 I=1,LMTOT
      XXMAT1(I,I)=CONE+XXMAT1(I,I)
    3 XXMAT2(I,I)=CONE+XXMAT2(I,I)
      RETURN
      END subroutine
C=======================================================================
C=======================================================================
C=======================================================================
      SUBROUTINE XMAT(XODD,XEVEN,LMAX,KAPPA,AK,ELM,EMACH)

C     ------------------------------------------------------------------
C     XMAT CALCULATES THE MATRIX DESCRIBING MULTIPLE SCATERING  WITHIN
C     A  LAYER, RETURNING  IT AS :  XODD,  CORRESPONDING  TO  ODD  L+M,
C     WITH LM=(10),(2-1),(21),... AND XEVEN, CORRESPONDING TO EVEN L+M,
C     WITH LM=(00),(1-1),(11),(2-2),(20),(22),...
C     THE  PROGRAM  ASSUMES  THAT  THE  LAYER IS A BRAVAIS LATTICE. THE
C     SUMMATION OVER THE LATTICE FOLLOWS THE EWALD METHOD  SUGGESTED BY
C     KAMBE. EMACH IS THE MACHINE ACCURACY.
C     ------------------------------------------------------------------
C
C ..  PARAMETER STATEMENTS  ..
C
      INTEGER   NDEND
      PARAMETER (NDEND=1240)
C
C ..  SCALAR ARGUMENTS  ..
C
      INTEGER    LMAX
      REAL(dp)   EMACH
      COMPLEX(dp) KAPPA
C
C ..  ARRAY ARGUMENTS  ..
C
      REAL(dp)   AK(2),ELM(:)
      COMPLEX(dp) XODD(:,:),XEVEN(:,:)
C
C ..  LOCAL SCALARS  ..
C
      INTEGER    L2MAX,LL2,II,I,NNDLM,K,KK,L,MM,NN,M,J1,J2,I1,I2,I3,N1
      INTEGER    NA,LLL,N,IL,NM,IN,L2,IL2,M2,IL3,L3,M3,LA1,LB1,LA11,LB11
      INTEGER    LL,J,L1
      REAL(dp)   AB1,AB2,AC,ACSQ,AD,AL,AN,AN1,AN2,AP,AP1,AP2,AR,B
      REAL(dp)   DNORM,RTPI,RTV,TEST,TEST1,TEST2,TV
      COMPLEX(dp) ALPHA,RTA,RTAI,KAPSQ,KANT,KNSQ,XPK,XPA,CF,CP,CX,CZ
      COMPLEX(dp) Z,ZZ,W,WW,A,ACC,GPSQ,GP,BT,AA,AB,U,U1,U2,GAM
      COMPLEX(dp) GK,GKK,SD,ALM
C
C ..  LOCAL ARRAYS  ..
C
      REAL(dp)   DENOM(NDEND),R(2),B1(2),B2(2),AKPT(2),FAC(4*lmax+1)
      COMPLEX(dp) GKN(lmax+1),AGK(2*lmax+1),XPM(2*lmax+1),
     &  PREF((lmax+1)**2)
      complex(dp) DLM((lmax+1)*(2*lmax+1))
c
C ..  ARRAYS IN COMMON  ..
C
      REAL(dp)  AR1(2),AR2(2)
      COMMON/X1/AR1,AR2
C----------------------------------------------------------------------

C
C     AK(1)  AND  AK(2)  ARE THE X  AND Y COMPONENTS OF THE
C     MOMENTUM PARALLEL TO THE SURFACE, MODULO A RECIPROCAL
C     LATTICE VECTOR
C
      RTPI=SQRT(PI)
      KAPSQ=KAPPA*KAPPA
C
C     THE FACTORIAL  FUNCTION  IS TABULATED  IN FAC . THE ARRAY
C     DLM WILL CONTAIN NON-ZERO,I.E. L+M EVEN,VALUES AS DEFINED
C     BY KAMBE.DLM=DLM1+DLM2+DLM3.WITH LM=(00),(1-1),(11),(2-2)...
C
      L2MAX=LMAX+LMAX
      LL2=L2MAX+1
      FAC(1)=1.0_dp
      II=L2MAX+L2MAX
      DO 2 I=1,II
   2  FAC(I+1)=dble(I)*FAC(I)
      NNDLM=L2MAX*(L2MAX+3)/2+1
      DO 3 I=1,NNDLM
   3  DLM(I)=CZERO
C
C     THE FORMULA OF KAMBE FOR THE SEPARATION CONSTANT,ALPHA,IS
C     USED,SUBJECT TO A RESTRICTION WHICH IS IMPOSED TO CONTROL
C     LATER ROUNDING ERRORS
C
      TV=ABS(AR1(1)*AR2(2)-AR1(2)*AR2(1))
      ALPHA=TV/(4.0_dp*PI)*KAPSQ
      AL=ABS(ALPHA)
      IF(EXP(AL)*EMACH-5.0D-5)5,5,4
   4  AL=LOG(5.0D-5/EMACH)
   5  ALPHA=cmplx_dp(AL,0.0_dp)
      RTA=SQRT(ALPHA)
C
C     DLM1 , THE  SUM  OVER  RECIPROCAL   LATTICE  VECTORS  , IS
C     CALCULATED FIRST. THE PREFACTOR P1 IS  TABULATED  FOR EVEN
C     VALUES OF L+|M|,THUS LM=(00),(11),(2 0),(22),(L2MAX,L2MAX)
C     THE  FACTORIAL  FACTOR  F1  IS SIMULTANEOUSLY TABULATED IN
C     DENOM,FOR ALL VALUES OF N=0,(L-|M|)/2
C
      K=1
      KK=1
      AP1=-2.0_dp/TV
      AP2=-1.0_dp
      CF=CI/KAPPA
      DO 8 L=1,LL2
      AP1=AP1/2.0_dp
      AP2=AP2+2.0_dp
      CP=CF
      MM=1
      IF(MOD  (L,2))7,6,7
   6  MM=2
      CP=CI*CP
   7  NN=(L-MM)/2+2
      DO 8 M=MM,L,2
      J1=L+M-1
      J2=L-M+1
      AP=AP1*SQRT(AP2*FAC(J1)*FAC(J2))
      PREF(KK)=AP*CP
      CP=-CP
      KK=KK+1
      NN=NN-1
      DO 8 I=1,NN
      I1=I
      I2=NN-I+1
      I3=NN+M-I
      DENOM(K)=1.0_dp/(FAC(I1)*FAC(I2)*FAC(I3))
   8  K=K+1
C
C     THE  RECIPROCAL  LATTICE IS  DEFINED BY  B1,B2 . THE  SUMMATION
C     BEGINS WITH THE ORIGIN POINT OF THE LATTICE , AND  CONTINUES IN
C     STEPS OF 8*N1 POINTS , EACH  STEP INVOLVING THE  PERIMETER OF A
C     PARALLELOGRAM OF LATTICE POINTS ABOUT THE ORIGIN,OF SIDE 2*N1+1
C     EACH STEP BEGINS AT LABEL 9.
C     AKPT=THE CURRENT LATTICE VECTOR IN THE SUM
C
      RTV=2.0_dp*PI/TV
      B1(1)=-AR1(2)*RTV
      B1(2)=AR1(1)*RTV
      B2(1)=-AR2(2)*RTV
      B2(2)=AR2(1)*RTV
      TEST1=1.0D6
      II=1
      N1=-1
   9  N1=N1+1
      NA=N1+N1+II
      AN1=dble(N1)
      AN2=-AN1-1.0_dp
      DO 22 I1=1,NA
      AN2=AN2+1.0_dp
      DO 21 I2=1,4
C     WRITE(16,307) I1,I2
C 307 FORMAT(33X,'I1=',I2,' , I2=',I2/33X,12('='))
      AN=AN1
      AN1=-AN2
      AN2=AN
      AB1=AN1*B1(1)+AN2*B2(1)
      AB2=AN1*B1(2)+AN2*B2(2)
      AKPT(1)=AK(1)+AB1
      AKPT(2)=AK(2)+AB2
C
C     FOR  EVERY LATTICE VECTOR OF THE SUM, THREE SHORT ARRAYS ARE
C     INITIALISED AS BELOW. AND USED AS TABLES:
C     XPM(M) CONTAINS VALUES OF XPK**|M|
C     AGK(I) CONTAINS VALUES OF (AC/KAPPA)**I
C     GKN(N) CONTAINS VALUES OF (GP/KAPPA)**(2*N-1)*GAM(N,Z)
C     WHERE L=0,L2MAX;M=-L,L;N=0,(L-|M|)/2;I=L-2*N
C     GAM IS THE INCOMPLETE GAMMA FUNCTION, WHICH IS CALCULATED BY
C     RECURRENCE  FROM  THE VALUE  FOR N=0, WHICH  IN TURN CAN  BE
C     EXPRESSED IN TERMS OF THE COMPLEX ERROR FUNCTION CERF
C     AC=MOD(AKPT). NOTE SPECIAL ACTION IF AC=0
C
      ACSQ=AKPT(1)*AKPT(1)+AKPT(2)*AKPT(2)
      GPSQ=KAPSQ-ACSQ
      IF(ABS(GPSQ)<EMACH*EMACH)   THEN
      WRITE(7,100)
  100 FORMAT(13X,'FATAL ERROR FROM XMAT:'/3X,'GPSQ IS TOO SMALL.'
     & /3X,'GIVE A SMALL BUT NONZERO VALUE FOR "EPSILON"'/3X,
     & 'IN THE DATA STATEMENT OF THE MAIN PROGRAM.'
     & /3X,'THIS DEFINES A SMALL IMAGINARY PART'
     & /3X,'IN THE FREQUENCY OR WAVELENGTH VALUE.')
      STOP
      ENDIF
      AC=SQRT(ACSQ)
      GP=SQRT(GPSQ)
      XPK=CZERO
      GK=CZERO
      GKK=cmplx_dp(1.0_dp,0.0_dp)
      IF(AC-EMACH)11,11,10
  10  XPK=cmplx_dp(AKPT(1)/AC,AKPT(2)/AC)
      GK=AC/KAPPA
      GKK=GPSQ/KAPSQ
  11  XPM(1)=cmplx_dp(1.0_dp,0.0_dp)
      AGK(1)=cmplx_dp(1.0_dp,0.0_dp)
      DO 12 I=2,LL2
      XPM(I)=XPM(I-1)*XPK
  12  AGK(I)=AGK(I-1)*GK
      CF=KAPPA/GP
      ZZ=-ALPHA*GKK
      CZ=SQRT(-ZZ)
      Z=-CI*CZ
      CX=EXP(-ZZ)
      GAM=RTPI*CERF(CZ)
      GKN(1)=CF*CX*GAM
      BT=Z
      B=0.5_dp
      LLL=L2MAX/2+1
      DO 13 I=2,LLL
      BT=BT/ZZ
      B=B-1.0_dp
      GAM=(GAM-BT)/B
      CF=CF*GKK
  13  GKN(I)=CF*CX*GAM
C
C     THE CONTRIBUTION TO THE SUM DLM1 FOR A PARTICULAR
C     RECIPROCAL LATTICE VECTOR IS NOW ACCUMULATED INTO
C     THE  ELEMENTS OF DLM,NOTE SPECIAL ACTION IF  AC=0
C
      K=1
      KK=1
      DO 19 L=1,LL2
      MM=1
      IF(MOD  (L,2))15,14,15
  14  MM=2
  15  N=(L*L+MM)/2
      NN=(L-MM)/2+2
      DO 19 M=MM,L,2
      ACC=CZERO
      NN=NN-1
      IL=L
      DO 16 I=1,NN
      ACC=ACC+DENOM(K)*AGK(IL)*GKN(I)
      IL=IL-2
  16  K=K+1
      ACC=PREF(KK)*ACC
      IF(AC-1.0D-6)17,17,165
 165  DLM(N)=DLM(N)+ACC/XPM(M)
      IF(M-1)17,18,17
  17  NM=N-M+1
      DLM(NM)=DLM(NM)+ACC*XPM(M)
  18  KK=KK+1
  19  N=N+1
      IF(II)21,21,22
  21  CONTINUE
  22  II=0
C
C     AFTER EACH STEP OF THE SUMMATION A TEST ON THE
C     CONVERGENCE  OF THE  ELEMENTS OF  DLM IS  MADE
C
      TEST2=0.0_dp
      DO 23 I=1,NNDLM
      DNORM=ABS(DLM(I))
  23  TEST2=TEST2+DNORM*DNORM
      TEST=ABS((TEST2-TEST1)/TEST1)
      TEST1=TEST2
      IF(TEST-0.001D0)27,27,24
  24  IF(N1-10)9,25,25
  25  WRITE(16,26)N1
  26  FORMAT(29H**DLM1,S NOT CONVERGED BY N1=,I2)
      GOTO 285
  27  WRITE(16,28)N1
  28  FORMAT(25H DLM1,S CONVERGED BY N1 =,I2)
C     WRITE(16,250)DLM
C250  FORMAT(5H0DLM1,//,45(2E13.5,/))
C
C     DLM2, THE SUM OVER REAL SPACE LATTICE VECTORS, BEGINS WITH
C     THE ADJUSTMENT OF THE ARRAY PREF, TO CONTAIN VALUES OF THE
C     PREFACTOR  'P2' FOR LM=(00),(11),(20),(22),...
C
 285  KK=1
      AP1=TV/(4.0_dp*PI)
      CF=KAPSQ/CI
      DO 31 L=1,LL2
      CP=CF
      MM=1
      IF(MOD  (L,2))30,29,30
  29  MM=2
      CP=-CI*CP
  30  J1=(L-MM)/2+1
      J2=J1+MM-1
      IN=J1+L-2
      AP2=((-1.0_dp)**IN)*AP1
      DO 31 M=MM,L,2
      AP=AP2/(FAC(J1)*FAC(J2))
      PREF(KK)=AP*CP*PREF(KK)
      J1=J1-1
      J2=J2+1
      AP2=-AP2
      CP=-CP
  31  KK=KK+1
C
C     THE SUMMATION PROCEEDS IN STEPS OF 8*N1 LATTICE POINTS
C     AS BEFORE, BUT THIS  TIME EXCLUDING  THE ORIGIN  POINT
C     R=THE CURRENT LATTICE VECTOR IN THE SUM
C     AR=MOD(R)
C
      N1=0
  32  N1=N1+1
      NA=N1+N1
      AN1=dble(N1)
      AN2=-AN1-1.0_dp
      DO 40 I1=1,NA
      AN2=AN2+1.0_dp
      DO 40 I2=1,4
      AN=AN1
      AN1=-AN2
      AN2=AN
      R(1)=AN1*AR1(1)+AN2*AR2(1)
      R(2)=AN1*AR1(2)+AN2*AR2(2)
      AR=SQRT(R(1)*R(1)+R(2)*R(2))
      XPK=cmplx_dp(R(1)/AR,R(2)/AR)
      XPM(1)=cmplx_dp(1.0_dp,0.0_dp)
      DO 33 I=2,LL2
  33  XPM(I)=XPM(I-1)*XPK
      AD=AK(1)*R(1)+AK(2)*R(2)
      SD=EXP(-AD*CI)
C
C     FOR EACH LATTICE VECTOR THE INTEGRAL 'U' IS OBTAINED
C     FROM THE RECURRENCE RELATION IN L SUGGESTED BY KAMBE
C     U1 AND U2 ARE THE  INITIAL TERMS OF THIS RECURRENCE,
C     FOR L#-1 AND L=0, AND THEY ARE EVALUATED IN TERMS OF
C     THE COMPLEX ERROR FUNCTION CERF
C
      KANT=0.5_dp*AR*KAPPA
      KNSQ=KANT*KANT
      Z=CI*KANT/RTA
      ZZ=RTA-Z
      Z=RTA+Z
      WW=CERF(-ZZ)
      W=CERF(Z)
      AA=0.5_dp*RTPI*(W-WW)/CI
      AB=0.5_dp*RTPI*(W+WW)
      A=ALPHA-KNSQ/ALPHA
      XPA=EXP(A)
      U1=AA*XPA
      U2=AB*XPA/KANT
C
C     THE CONTRIBUTION TO DLM2 FROM A PARTICULAR LATTICE
C     VECTOR  IS  ACCUMULATED INTO  THE ELEMENTS OF  DLM
C     THIS PROCEDURE INCLUDES THE TERM (KANT**L) AND THE
C     RECURRENCE FOR THE INTEGRAL 'U'
C
      KK=1
      AL=-0.5_dp
      CP=RTA
      CF=cmplx_dp(1.0_dp,0.0_dp)
      DO 39 L=1,LL2
      MM=1
      IF(MOD  (L,2))35,34,35
  34  MM=2
  35  N=(L*L+MM)/2
      DO 38 M=MM,L,2
      ACC=PREF(KK)*U2*CF*SD
      DLM(N)=DLM(N)+ACC/XPM(M)
      IF(M-1)36,37,36
  36  NM=N-M+1
      DLM(NM)=DLM(NM)+ACC*XPM(M)
  37  KK=KK+1
  38  N=N+1
      AL=AL+1.0_dp
      CP=CP/ALPHA
      U=(AL*U2-U1+CP*XPA)/KNSQ
      U1=U2
      U2=U
  39  CF=KANT*CF
  40  CONTINUE
C
C     AFTER EACH STEP OF THE SUMMATION A TEST ON THE
C     CONVERGENCE OF THE ELEMENTS OF DLM IS MADE
C
      TEST2=0.0_dp
      DO 41 I=1,NNDLM
      DNORM=ABS(DLM(I))
  41  TEST2=TEST2+DNORM*DNORM
      TEST=ABS((TEST2-TEST1)/TEST1)
      TEST1=TEST2
      IF(TEST-0.001D0)45,45,42
  42  IF(N1-10)32,43,43
  43  WRITE(16,44)N1
  44  FORMAT(31H0**DLM2,S NOT CONVERGED BY N1 =,I2)
      GOTO 465
  45  WRITE(16,46)N1
  46  FORMAT(24H DLM2,S CONVERGED BY N1=,I2)
C
C     THE TERM DLM3 HAS A NON-ZERO CONTRIBUTION  ONLY
C     WHEN L=M=0.IT IS EVALUATED HERE IN TERMS OF THE
C     COMPLEX ERROR FUNCTION CERF
C
 465  XPA=EXP(-ALPHA)
      RTAI=1.0_dp/(RTPI*RTA)
      ACC=KAPPA*(CI*(XPA-CERF(RTA))-RTAI)/XPA
      AP=-0.5_dp/RTPI
      DLM(1)=DLM(1)+AP*ACC
C
C     FINALLY THE ELEMENTS OF DLM ARE MULTIPLIED BY THE
C     FACTOR (-1.0_dp)**((M+|M|)/2)
C
      DO 47 L=2,LL2,2
      N=L*L/2+1
      DO 47 M=2,L,2
      DLM(N)=-DLM(N)
  47  N=N+1
C     WRITE(16,251) DLM
C 251 FORMAT(15H0DLM1+DLM2+DLM3,//45(2E13.5,/))
C
C     SUMMATION OVER THE CLEBSCH-GORDON TYPE COEFFICIENTS
C     ELM PROCEEDS, FIRST FOR  XODD, AND THEN  FOR XEVEN.
C     THIS GIVES THE KAMBE ELEMENTS  A(L2,M2;L3,M3) WHICH
C     GIVE THE ELEMENTS  X(L3,M3;L2,M2) OF XODD AND XEVEN
C
      K=1
      II=0
  48  LL=LMAX+II
      I=1
      DO 56 IL2=1,LL
      L2=IL2-II
      M2=-L2+1-II
      DO 56 I2=1,IL2
      J=1
      DO 55 IL3=1,LL
      L3=IL3-II
      M3=-L3+1-II
      DO 55 I3=1,IL3
      ALM=CZERO
      LA1=MAX0(IABS(L2-L3),IABS(M2-M3))
      LB1=L2+L3
      N=(LA1*(LA1+2)+M2-M3+2)/2
      NN=2*LA1+4
      LB11=LB1+1
      LA11=LA1+1
      DO 49 L1=LA11,LB11,2
      ALM=ALM+ELM(K)*DLM(N)
      N=N+NN
      NN=NN+4
  49  K=K+1
      ALM=ALM/KAPPA
      IF(I-J)51,50,51
  50  ALM=ALM+CI
  51  IF(II)52,52,53
  52  XODD(J,I)=CI*ALM
      GOTO 54
  53  XEVEN(J,I)=CI*ALM
  54  M3=M3+2
  55  J=J+1
      M2=M2+2
  56  I=I+1
      IF(II)57,57,58
  57  II=1
      GOTO 48
  58  CONTINUE
      RETURN
      END subroutine
C=======================================================================
      SUBROUTINE ZGE(A,INT,N,NC,EMACH)

C     ------------------------------------------------------------------
C     ZGE IS A STANDARD SUBROUTINE TO PERFORM GAUSSIAN ELIMINATION ON
C     A NC*NC MATRIX 'A' PRIOR  TO INVERSION, DETAILS STORED IN 'INT'
C     ------------------------------------------------------------------
C
C ..  SCALAR ARGUMENTS  ..
C
      INTEGER N,NC
      REAL(dp) EMACH
C
C ..  ARRAY ARGUMENTS  ..
C
      INTEGER    INT(NC)
      COMPLEX(dp) A(NC,NC)
C
C ..  LOCAL SCALARS  ..
C
      INTEGER    I,II,IN,J,K
      COMPLEX(dp) YR,DUM
C     ------------------------------------------------------------------
C
      DO 10 II=2,N
      I=II-1
      YR=A(I,I)
      IN=I
      DO 2 J=II,N
      IF(ABS(YR)-ABS(A(J,I)))1,2,2
   1  YR=A(J,I)
      IN=J
   2  CONTINUE
      INT(I)=IN
      IF(IN-I)3,5,3
   3  DO 4 J=I,N
      DUM=A(I,J)
      A(I,J)=A(IN,J)
   4  A(IN,J)=DUM
   5  IF(ABS(YR)-EMACH)10,10,6
   6  DO 9 J=II,N
      IF(ABS(A(J,I))-EMACH)9,9,7
   7  A(J,I)=A(J,I)/YR
      DO 8 K=II,N
   8  A(J,K)=A(J,K)-A(I,K)*A(J,I)
   9  CONTINUE
  10  CONTINUE
      RETURN
      END subroutine
C=======================================================================
      SUBROUTINE ZSU(A,INT,X,N,NC,EMACH)

C     ------------------------------------------------------------------
C     ZSU  IS  A STANDARD BACK-SUBSTITUTION  SUBROUTINE  USING THE
C     OUTPUT OF ZGE TO CALCULATE  A-INVERSE TIMES X, RETURNED IN X
C     ------------------------------------------------------------------
C
C ..  SCALAR ARGUMENTS  ..
C
      INTEGER N,NC
      REAL(dp) EMACH
C
C ..  ARRAY ARGUMENTS  ..
C
      INTEGER    INT(NC)
      COMPLEX(dp) A(NC,NC),X(NC)
C
C ..  LOCAL SCALARS  ..
C
      INTEGER    I,II,IN,J,IJ
      COMPLEX(dp) DUM
C     ------------------------------------------------------------------
C
      DO 5 II=2,N
      I=II-1
      IF(INT(I)-I)1,2,1
   1  IN=INT(I)
      DUM=X(IN)
      X(IN)=X(I)
      X(I)=DUM
   2  DO 4 J=II,N
      IF(ABS(A(J,I))-EMACH)4,4,3
   3  X(J)=X(J)-A(J,I)*X(I)
   4  CONTINUE
   5  CONTINUE
      DO 10 II=1,N
      I=N-II+1
      IJ=I+1
      IF(I-N)6,8,6
   6  DO 7 J=IJ,N
   7  X(I)=X(I)-A(I,J)*X(J)
   8  IF(ABS(A(I,I))-EMACH*1.0D-7)9,10,10
   9  A(I,I)=EMACH*1.0D-7*(1.0_dp,1.0_dp)
  10  X(I)=X(I)/A(I,I)
      RETURN
      END subroutine
C=======================================================================
C=======================================================================
C=======================================================================
      SUBROUTINE CNAA(NDIM,N,AR,AI,EVR,EVI,VECR,VECI,IERR)

C     ------------------------------------------------------------------
C     'EISPACK'  IS A  COLLECTION  OF CODES FOR  SOLVING  THE ALGEBRAIC
C     EIGENVALUE  PROBLEM.  THE ORIGINAL  ALGOL  CODES WERE  WRITTEN BY
C     J. H. WILKINSON, ET.AL., AND SUBSEQUENTLY  TRANSLATED TO  FORTRAN
C     AND TESTED AT ARGONNE NATIONAL LABORATORY.
C
C     THIS   SUBROUTINE  COMPUTES  ALL  EIGENVALUES  AND  CORRESPONDING
C     EIGENVECTORS  OF  AN  ARBITRARY   COMPLEX  MATRIX.  THE MATRIX IS
C     BALANCED BY EXACT NORM  REDUCING  SIMILARITY  TRANSFORMATIONS AND
C     THEN  IS  REDUCED  TO  COMPLEX  HESSENBERG   FORM  BY  STABILIZED
C     ELEMENTARY SIMILARITY TRANSFORMATIONS. A MODIFIED LR ALGORITHM IS
C     USED TO COMPUTE THE EIGENVALUES OF THE HESSENBERG MATRIX.
C
C       ON INPUT--->
C          NDIM     MUST BE THE ROW DIMENSION OF THE ARRAYS AR,AI,VECR,
C                   AND VECI IN THE CALLING PROGRAM DIMENSION STATEMENT
C          N        IS THE ORDER OF THE MATRIX. N MUST NOT EXCEED NDIM.
C                   N*NDIM  MUST NOT EXCEED 22500=150*150=53744(OCTAL).
C                   N MUST NOT EXCEED 150.  N MAY BE 1.
C          AR,AI    ARRAYS WITH  EXACTLY  NDIM  ROWS  AND  AT  LEAST  N
C                   COLUMNS.  THE LEADING N BY N SUBARRAYS MUST CONTAIN
C                   THE REAL AND  IMAGINARY  PARTS  RESPECTIVELY OF THE
C                   ARBITRARY COMPLEX MATRIX WHOSE EIGENSYSTEM IS TO BE
C                   COMPUTED.
C
C        ON OUTPUT--->
C          EVR,EVI    CONTAIN THE REAL AND IMAGINARY PARTS RESPECTIVELY
C                     OF THE COMPUTED EIGENVALUES.  THE EIGENVALUES ARE
C                     NOT ORDERED IN ANY WAY.
C          VECR,VECI  CONTAIN IN THE LEADING N BY N  SUBARRAYS THE REAL
C                     AND IMAGINARY PARTS RESPECTIVELY  OF THE COMPUTED
C                     EIGENVECTORS.  THE J-TH COLUMNS  OF VECR AND VECI
C                     CONTAIN THE  EIGENVECTOR  ASSOCIATED  WITH EVR(J)
C                     AND  EVI(J).  THE EIGENVECTORS ARE NOT NORMALIZED
C                     IN ANY WAY.
C          IERR       IS A STATUS CODE.
C                   --NORMAL CODE.
C                     0 MEANS THE LR ITERATIONS CONVERGED.
C                   --ABNORMAL CODES.
C                     J MEANS THE J-TH EIGENVALUE HAS NOT BEEN FOUND IN
C                     30 ITERATIONS. THE FIRST J-1 ELEMENTS OF EVR  AND
C                     EVI CONTAIN THOSE EIGENVALUES  ALREADY  FOUND. NO
C                     EIGENVECTORS ARE COMPUTED.
C                    -1 MEANS THE INPUT VALUES OF N, NDIM ARE TOO LARGE
C                     OR INCONSISTENT.
C          AR,AI      ARE DESTROYED.
C     ------------------------------------------------------------------
C
C ..  SCALAR ARGUMENTS  ..
C
      INTEGER IERR,N,NDIM
C
C ..  ARRAY ARGUMENTS  ..
C
      REAL(dp)AR(NDIM,N),AI(NDIM,N),EVR(N),EVI(N)
      complex(dp)::a(ndim, n), w(n), vr(ndim,n)
      REAL(dp)VECR(NDIM,N),VECI(NDIM,N)
C
C ..  LOCAL SCALARS  ..
C
      INTEGER   I,IERRPI,IGH,LOW,NMIERR
C
C ..  LOCAL ARRAYS  ..
C
      INTEGER INT(270)
      REAL(dp)SCALE(270)
C     ------------------------------------------------------------------
C
      IF(NDIM<N .OR. N<1) GO TO 10
      IF(N*NDIM > 72900) GO TO 10
      a = ar + ci*ai
      call zgebal_wrap(a,w,vr, scale, low, igh)
      ar = dble(a)
      ai = aimag(a)
!       call zgehrd_wrap(a, w, vr, scale, low, igh)
!     ierr = 0
!     ar = dble(a)
!     ai = aimag(a)
!     evr = dble(w)
!     evi = aimag(w)
!     vecr = dble(vr)
!     veci = dble(vr)

      CALL COMHES(NDIM,N,LOW,IGH,AR,AI,INT)
      CALL COMLR2(NDIM,N,LOW,IGH,INT,AR,AI,EVR,EVI,VECR,VECI,IERR)
      IF(IERR==0) GO TO 2
!     CALL ERRCHK(54,54HIN CNAA  , SOME EIGENVALUE NOT FOUND IN 30 ITERA
!    1TIONS.)
      IF(IERR==N) GO TO 20
      NMIERR = N - IERR
      DO 1 I=1,NMIERR
      IERRPI = IERR + I
      EVR(I) = EVR(IERRPI)
 1    EVI(I) = EVI(IERRPI)
      GO TO 20
 2    CALL CBABK2(NDIM,N,LOW,IGH,SCALE,N,VECR,VECI)
      GO TO 20
10    write(*,*)"IN CNAA INPUT DIM IN ERROR OR MATRIX IS TOO BIG."
!CALL ERRCHK(58,58HIN CNAA  , INPUT DIMENSIONS IN ERROR OR MATRIX I
!    1S TOO BIG.)
      IERR=-1
20    IF(IERR > 0) IERR = N-IERR+1
      RETURN
      END subroutine
C=======================================================================
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
      CALL XMAT(XODD,XEVEN,LMAX,KAPPA,AK,ELM,EMACH)
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
      SUBROUTINE REDUCE(AR1,AR2,AK,IGMAX,G,IG0,EMACH)

C----------------------------------------------------------------------
C     GIVEN THE PRIMITIVE VECTORS AR1,AR2 OF A 2D LATTICE (IN UNITS OF
C     ALPHA), THIS  SUBROUTINE  REDUCES A WAVECTOR "AK" (IN  UNITS  OF
C     2*PI/ALPHA) WITHIN THE SBZ BY ADDING AN APPROPRIATE  RECIPROCAL-
C     LATTICE VECTOR G(IG0)
C----------------------------------------------------------------------
C ..  SCALAR ARGUMENTS  ..
C
      INTEGER IGMAX,IG0
      REAL(dp) EMACH
C
C ..  ARRAY ARGUMENTS  ..
C
      REAL(dp) AR1(2),AR2(2),AK(2),G(:,:)
C
C ..  LOCAL SCALARS  ..
C
      INTEGER I,J,N,I1,I2
      REAL(dp) D,B,P,AFI,AX,AY,AKX,AKY,FI0,AM,BM,ALPHA,RA
C
C ..  LOCAL ARRAYS ..
C
      REAL(dp) VX(6),VY(6),FI(6),X(6),Y(6)
C----------------------------------------------------------------------
      ALPHA=AR1(1)
      RA=2.0_dp*PI/ALPHA
      D=AR2(1)
      B=AR2(2)
      IF( (DABS(D)-0.5_dp)>EMACH) STOP 'IMPROPER LATTICE VECTORS'
      IF((DABS(DABS(D)-0.5_dp)<EMACH).AND.(DABS(DABS(B)-
     &    DSQRT(3.0_dp)/2.0_dp)>EMACH)) THEN
      B=2.0_dp*B           ! CENTRED RECTANGULAR LATTICE
        IF((DABS(B)-1.0_dp)<0.0_dp) THEN
      VX(1)= 1.0_dp
      VY(1)=0.5_dp*(1.0_dp/B - B)
      VX(2)= 0.0_dp
      VY(2)=0.5_dp*(1.0_dp/B + B)
      VX(3)=-1.0_dp
      VY(3)= VY(1)
      VX(4)=-1.0_dp
      VY(4)=-VY(1)
      VX(5)= 0.0_dp
      VY(5)=-VY(2)
      VX(6)= 1.0_dp
      VY(6)=-VY(1)
        ELSE
      VX(1)=0.5_dp+0.5_dp/B/B
      VY(1)=0.0_dp
      VX(2)=0.5_dp-0.5_dp/B/B
      VY(2)=1.0_dp/B
      VX(3)=-VX(2)
      VY(3)=VY(2)
      VX(4)=-VX(1)
      VY(4)=0.0_dp
      VX(5)=-VX(2)
      VY(5)=-VY(2)
      VX(6)=VX(2)
      VY(6)=-VY(2)
        ENDIF
      ELSE             !OBLIQUE OR HEXAGONAL LATTICE
        IF(D>0.0_dp) THEN
      P     =0.5_dp*D*(D-1)/B/B
      VX(1)=0.5_dp-P
      VY(1)=0.5_dp*(1.0_dp-2.0_dp*D)/B
      VX(2)=0.5_dp+P
      VY(2)=0.5_dp/B
        ELSE
      P     = 0.5_dp*D*(D+1)/B/B
      VX(1)= 0.5_dp-P
      VY(1)=-0.5_dp*(1.0_dp+2.0_dp*D)/B
      VX(2)= 0.5_dp+P
      VY(2)=-0.5_dp/B
        ENDIF
      VX(3)=-VX(2)
      VY(3)= VY(2)
      VX(4)=-VX(1)
      VY(4)=-VY(1)
      VX(5)=-VX(2)
      VY(5)=-VY(2)
      VX(6)= VX(2)
      VY(6)=-VY(2)
      ENDIF
      N=6
      DO 1 I=1,N
      X(I)=VX(I)*RA
      Y(I)=VY(I)*RA
  1   CONTINUE
      IF(DABS(D)<EMACH) THEN
      N=4              !RECTANGULAR OR SQUARE LATTICE
      IF(  B>0.0_dp) THEN
      X(1)=VX(6)*RA
      Y(1)=VY(6)*RA
      X(2)=VX(4)*RA
      Y(2)=VY(4)*RA
      X(4)=VX(1)*RA
      Y(4)=VY(1)*RA
      ELSE
      X(2)=VX(3)*RA
      Y(2)=VY(3)*RA
      X(3)=VX(4)*RA
      Y(3)=VY(4)*RA
      X(4)=VX(6)*RA
      Y(4)=VY(6)*RA
      ENDIF
      ENDIF
C*****VERTICES ARE ARRANGED IN ASCENDING ORDER OF THE POLAR ANGLE FI
      DO 2 I=1,N
      FI(I)=DATAN2(Y(I),X(I))
      IF(FI(I)<0.0_dp) FI(I)=FI(I) + 2.0_dp*PI
   2  CONTINUE
      DO 3  J=2,N
      AFI=FI(J)
      AX = X(J)
      AY = Y(J)
      DO 4  I=J-1,1,-1
      IF(FI(I)<=AFI) GOTO 5
      FI(I+1) =FI(I)
       X(I+1) =X(I)
       Y(I+1) =Y(I)
   4  CONTINUE
      I=0
   5  FI(I+1)=AFI
       X(I+1)=AX
       Y(I+1)=AY
   3  CONTINUE
C*****"AK" IS REDUCED WITHIN THE SBZ
      IG0=1
   6  CONTINUE
      AKX=AK(1)-G(1,IG0)
      AKY=AK(2)-G(2,IG0)
      IF((ABS(AKX)<EMACH).AND.(ABS(AKY)<EMACH)) RETURN
      FI0=DATAN2(AKY,AKX)   ! FIND POLAR ANGLES OF THE WAVEVECTOR
      IF(FI0<0.0_dp) FI0=FI0+2.0_dp*PI
        I1=N
        I=1
    7       CONTINUE
        I2=I
        IF(FI0<FI(I))  GO TO 8
        I=I+1
        I1=I2
        IF(I<=N) GO TO 7
        I1=N
        I2=1
    8       CONTINUE
      AM=ABS(Y(I2)*X(I1)-X(I2)*Y(I1))
      BM=ABS((X(I1)-X(I2))*AKY+(Y(I2)-Y(I1))*AKX)
      IF(AM>=BM) THEN
      AK(1)=AKX
      AK(2)=AKY
      RETURN
      ENDIF
      IG0=IG0+1
      IF(IG0>IGMAX) STOP   'ERROR FROM REDUCE:  INSUFFICIENT NR. OF
     &RECIPROCAL LATTICE VECTORS '
      GOTO 6
      END subroutine
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
      DO 8 IPL=1,NPLAN(ICOMP)
      READ(10,206) S(ICOMP,IPL),MUSPH(ICOMP,IPL),EPSSPH(ICOMP,IPL)
      READ(10,*) DUMMY,(DL(I,ICOMP,IPL),I=1,3)
      READ(10,*) DUMMY,(DR(I,ICOMP,IPL),I=1,3)
    8 CONTINUE
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
      DO 5 IG1=1,IGMAX
      G(1,IG1)=NT1(IG1)*B1(1)+NT2(IG1)*B2(1)
      G(2,IG1)=NT1(IG1)*B1(2)+NT2(IG1)*B2(2)
      WRITE(6,215) IG1,NT1(IG1),NT2(IG1),VECMOD(IG1)
    5 CONTINUE
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
         DO 17 IGK1=1,IGKMAX
         DO 17 IGK2=1,IGKMAX
         WIL  (IGK1,IGK2)=QIR  (IGK1,IGK2)
         WIIL (IGK1,IGK2)=QIIR (IGK1,IGK2)
         WIIIL(IGK1,IGK2)=QIIIR(IGK1,IGK2)
         WIVL (IGK1,IGK2)=QIVR (IGK1,IGK2)
   17        CONTINUE
      DO 15 IPL=2,NPLAN(ICOMP)
      RAP=S(ICOMP,IPL)*KAPPA0/2.D0/PI
      CALL PCSLAB(LMAX,IGMAX,RAP,EPS1(ICOMP),EPSSPH(ICOMP,IPL),
     &           MU1(ICOMP),MUSPH(ICOMP,IPL),KAPPA,AK,
     &           DL(1,ICOMP,IPL),DR(1,ICOMP,IPL),G,ELM,A0,EMACH,
     &           QIR,QIIR,QIIIR,QIVR)
      CALL PAIR(IGKMAX,igkd,WIL,WIIL,WIIIL,WIVL,QIR,QIIR,QIIIR,QIVR)
   15 CONTINUE
         DO 18 IGK1=1,IGKMAX
         DO 18 IGK2=1,IGKMAX
         QIR  (IGK1,IGK2)=WIL  (IGK1,IGK2)
         QIIR (IGK1,IGK2)=WIIL (IGK1,IGK2)
         QIIIR(IGK1,IGK2)=WIIIL(IGK1,IGK2)
         QIVR (IGK1,IGK2)=WIVL (IGK1,IGK2)
   18        CONTINUE
            ENDIF
      IF(NLAYER(ICOMP)>=2) THEN
      DO 16 ILAYER=1,NLAYER(ICOMP)-1
      CALL PAIR(IGKMAX,igkd,QIR,QIIR,QIIIR,QIVR,QIR,QIIR,QIIIR,QIVR)
   16 CONTINUE
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
         DO 9 IU=1,NUNIT-1
             CALL PAIR(IGKMAX,igkd,QIL,QIIL,QIIIL,QIVL,
     &                 QIL,QIIL,QIIIL,QIVL)
    9        CONTINUE
   30        CONTINUE
             IF(KEMB==1) THEN
             CALL HOSLAB(IGMAX,KAPR,(KAPR+KAPOUT)/2.D0,KAPOUT,AK,G,VEC0,
     &                   VEC0,0.D0,QIR,QIIR,QIIIR,QIVR,EMACH)
         CALL PAIR(IGKMAX,igkd,QIL,QIIL,QIIIL,QIVL,QIR,QIIR,QIIIR,QIVR)
         DO 11 IGK1=1,IGKMAX
         DO 11 IGK2=1,IGKMAX
         QIR  (IGK1,IGK2)=QIL  (IGK1,IGK2)
         QIIR (IGK1,IGK2)=QIIL (IGK1,IGK2)
         QIIIR(IGK1,IGK2)=QIIIL(IGK1,IGK2)
         QIVR (IGK1,IGK2)=QIVL (IGK1,IGK2)
   11        CONTINUE
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
