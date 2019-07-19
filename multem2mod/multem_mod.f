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
!      use errfun, only: erf_pop
      use libmultem2b
      implicit none
      integer, parameter:: dp=kind(0.d0)
      real(dp), parameter :: pi=4.0_dp*ATAN(1.0_dp)
      complex(dp), parameter :: czero = (0.0_dp, 0.0_dp)
      complex(dp), parameter :: ci    = (0.0_dp, 1.0_dp)
      complex(dp), parameter :: cone  = (1.0_dp, 0.0_dp)
      complex(dp), parameter :: ctwo  = (2.0_dp, 0.0_dp)
      private
      public cerf, blm, ceven, codd, band, scat, pair, hoslab, pcslab,
     & reduce, lat2d, elmgen
      contains
C=======================================================================

C=======================================================================
      SUBROUTINE SCAT(IGMAX,ZVAL,AK,G,KAPIN,KAPOUT,EINCID,QI,QIII)

C     ------------------------------------------------------------------
C     THIS SUBROUTINE CALCULATES THE REFLECTIVITY, TRANSMITTANCE AND
C     ABSORBANCE OF A FINITE SLAB, CHARACTERIZED BY TRANSMISSION AND
C     REFLECTION MATRICES QI AND QIII, RESPECTIVELY.
C     ------------------------------------------------------------------
C
C ..  PARAMETER STATEMENTS ..
C
      INTEGER   IGD,IGKD
      PARAMETER (IGD=21,IGKD=2*IGD)
C
C ..  SCALAR ARGUMENTS  ..
C
      INTEGER    IGMAX
      REAL(dp)   ZVAL,KAPIN,KAPOUT
C
C ..  ARRAY ARGUMENTS ..
C
      REAL(dp)   AK(2),G(2,IGD)
      COMPLEX(dp) QI(IGKD,IGKD),QIII(IGKD,IGKD),EINCID(IGKD)
C
C ..  LOCAL SCALARS  ..
C
      INTEGER    IGK1,IG1,K1,IGK2,IGKMAX
      REAL(dp)   DOWN,REFLE,TRANS,ABSOR,GKZIN,GKZOUT,TES1
      COMPLEX(dp) ETRANS(IGKD),EREFLE(IGKD)
C     ------------------------------------------------------------------
      DOWN=0.D0
      REFLE=0.D0
      TRANS=0.D0
      IGKMAX=2*IGMAX
      IGK1=0
      DO 1 IG1=1,IGMAX
      TES1=(AK(1)+G(1,IG1))*(AK(1)+G(1,IG1))+(AK(2)+G(2,IG1))*(AK(2)+
     &            G(2,IG1))
      GKZIN =0.D0
      GKZOUT=0.D0
      IF(( KAPIN*KAPIN -TES1)>0.D0) GKZIN =SQRT( KAPIN*KAPIN -TES1)
      IF((KAPOUT*KAPOUT-TES1)>0.D0) GKZOUT=SQRT(KAPOUT*KAPOUT-TES1)
      DO 1 K1=1,2
      IGK1=IGK1+1
      ETRANS(IGK1)=CZERO
      EREFLE(IGK1)=CZERO
      DO 2 IGK2=1,IGKMAX
      ETRANS(IGK1)=ETRANS(IGK1)+QI  (IGK1,IGK2)*EINCID(IGK2)
      EREFLE(IGK1)=EREFLE(IGK1)+QIII(IGK1,IGK2)*EINCID(IGK2)
    2 CONTINUE
C
      DOWN =DOWN +EINCID(IGK1)*CONJG(EINCID(IGK1))*GKZIN
      TRANS=TRANS+ETRANS(IGK1)*CONJG(ETRANS(IGK1))*GKZOUT
      REFLE=REFLE+EREFLE(IGK1)*CONJG(EREFLE(IGK1))*GKZIN
    1 CONTINUE
      TRANS=TRANS/DOWN
      REFLE=REFLE/DOWN
      ABSOR=1.D0-TRANS-REFLE
      WRITE(8,101)  ZVAL,TRANS,REFLE,ABSOR
      WRITE(6,101)  ZVAL,TRANS,REFLE,ABSOR
      RETURN
C
  101 FORMAT(5E14.6)
      END subroutine
C======================================================================
      SUBROUTINE HOSLAB(IGMAX,KAPPA1,KAPPA2,KAPPA3,AK,G,DL,DR,D,
     &                  QI,QII,QIII,QIV,EMACH)

C-----------------------------------------------------------------------
C     THIS SUBROUTINE CALCULATES THE  Q-MATRICES FOR A HOMOGENEOUS
C     PLATE  '2' OF THICKNESS 'D', HAVING THE SEMI-INFINITE MEDIUM
C     '1' ON ITS LEFT AND THE SEMI-INFINITE MEDIUM '3' ON ITS RIGHT
C     ------------------------------------------------------------------
C
C  .. PARAMETER STATEMENTS ..
C
      INTEGER IGD,IGKD
      PARAMETER (IGD=21,IGKD=2*IGD)
C
C  .. SCALAR ARGUMENTS ..
C
      INTEGER    IGMAX
      REAL(dp)   EMACH,D
      COMPLEX(dp) KAPPA1,KAPPA2,KAPPA3
C
C  .. ARRAY AGUMENTS ..
C
      REAL(dp)   AK(2),G(2,IGD),DL(3),DR(3)
      COMPLEX(dp) QI(IGKD,IGKD),QII(IGKD,IGKD),QIII(IGKD,IGKD)
      COMPLEX(dp) QIV(IGKD,IGKD)

C
C  .. LOCAL SCALARS ..
C
      INTEGER    I,J,IA,IB,JA,IG1,IGKMAX
      REAL(dp)   GKKPAR
      COMPLEX(dp) GKKZ1,GKKZ2,GKKZ3,Z1,Z2,Z3,CQI,CQII
      COMPLEX(dp) CQIII,CQIV,DENOMA,DENOMB,GKKDUM
C
C  .. LOCAL ARRAYS ..
C
      COMPLEX(dp) T(4,2),R(4,2),X(4),P(4,2)
C     -----------------------------------------------------------------
      IGKMAX=2*IGMAX
      DO 1 IA=1,IGKMAX
      DO 1 IB=1,IGKMAX
      QI  (IA,IB)=CZERO
      QII (IA,IB)=CZERO
      QIII(IA,IB)=CZERO
      QIV (IA,IB)=CZERO
    1 CONTINUE
      X(1)=KAPPA1/KAPPA2
      X(2)=CONE/X(1)
      X(3)=KAPPA2/KAPPA3
      X(4)=CONE/X(3)
      DO 3 IG1=1,IGMAX
      GKKPAR=SQRT((AK(1)+G(1,IG1))*(AK(1)+G(1,IG1))+
     &            (AK(2)+G(2,IG1))*(AK(2)+G(2,IG1)))
      GKKZ1=SQRT(KAPPA1*KAPPA1-GKKPAR*GKKPAR)
      GKKZ2=SQRT(KAPPA2*KAPPA2-GKKPAR*GKKPAR)
      GKKZ3=SQRT(KAPPA3*KAPPA3-GKKPAR*GKKPAR)
      DO 9 J=1,2
      DENOMA=X(J)*X(J)*GKKZ2+GKKZ1
      DENOMB=     GKKZ2+GKKZ1
      IF(ABS(DENOMA)<EMACH.OR.ABS(DENOMB)<EMACH) GO TO 20
      R(J,1)=(GKKZ1-X(J)*X(J)*GKKZ2)/DENOMA
      R(J,2)=           (GKKZ1-GKKZ2)/DENOMB
      T(J,1)=CTWO*X(J)*GKKZ1/DENOMA
      T(J,2)=CTWO*GKKZ1/DENOMB
      GKKDUM=GKKZ1
      GKKZ1 =GKKZ2
      GKKZ2 =GKKDUM
 9    CONTINUE
      DO 10 J=3,4
      DENOMA=X(J)*X(J)*GKKZ3+GKKZ2
      DENOMB=          GKKZ3+GKKZ2
      IF(ABS(DENOMA)<EMACH.OR.ABS(DENOMB)<EMACH) GO TO 20
      R(J,1)=(GKKZ2-X(J)*X(J)*GKKZ3)/DENOMA
      R(J,2)=           (GKKZ2-GKKZ3)/DENOMB
      T(J,1)=CTWO*X(J)*GKKZ2/DENOMA
      T(J,2)=CTWO*GKKZ2/DENOMB
      GKKDUM=GKKZ2
      GKKZ2 =GKKZ3
      GKKZ3 =GKKDUM
 10   CONTINUE
      Z1=EXP(CI*GKKZ2*D)
      Z2=Z1*Z1
      DO 5 I=1,2
      Z3=CONE/(CONE-Z2*R(2,I)*R(3,I))
      P(1,I)=T(3,I)*Z3*Z1*T(1,I)
      P(2,I)=R(4,I)+T(4,I)*R(2,I)*T(3,I)*Z2*Z3
      P(3,I)=R(1,I)+T(2,I)*R(3,I)*T(1,I)*Z2*Z3
      P(4,I)=T(2,I)*Z3*Z1*T(4,I)
    5 CONTINUE
      CQI  =EXP(CI*((AK(1)+G(1,IG1))*(DL(1)+DR(1))+
     &              (AK(2)+G(2,IG1))*(DL(2)+DR(2))+
     &                 GKKZ1*DL(3)+GKKZ3*DR(3)))
      CQII =EXP(CTWO*CI*GKKZ3*DR(3))
      CQIII=EXP(CTWO*CI*GKKZ1*DL(3))
      CQIV =EXP(-CI*((AK(1)+G(1,IG1))*(DL(1)+DR(1))+
     &               (AK(2)+G(2,IG1))*(DL(2)+DR(2))-
     &                 GKKZ1*DL(3)-GKKZ3*DR(3)))
      DO 7 JA=1,2
      IA=2*IG1-2+JA
      QI  (IA,IA)=CQI  *P(1,JA)
      QII (IA,IA)=CQII *P(2,JA)
      QIII(IA,IA)=CQIII*P(3,JA)
      QIV (IA,IA)=CQIV *P(4,JA)
    7 CONTINUE
    3 CONTINUE
      RETURN
   20 STOP 'FATAL ERROR IN HOSLAB'
      END subroutine
C=======================================================================
      SUBROUTINE DLMKG(LMAX,A0,GK,SIGNUS,KAPPA,DLME,DLMH,EMACH)
C     ------------------------------------------------------------------
C     THIS SUBROUTINE CALCULATES THE COEFFICIENTS DLM(KG)
C     ------------------------------------------------------------------
C ..  PARAMETER STATEMENTS  ..
      INTEGER LMAXD,LMAX1D,LM1SQD
      PARAMETER(LMAXD=14,LMAX1D=LMAXD+1,LM1SQD=LMAX1D*LMAX1D)
C ..  ARGUMENTS  ..
      INTEGER    LMAX
      REAL(dp)   A0,SIGNUS,EMACH
      COMPLEX(dp) KAPPA
      COMPLEX(dp) DLME(2,LM1SQD),DLMH(2,LM1SQD),GK(3)
C  .. LOCAL
      INTEGER    K,II,L,M,I
      REAL(dp)   AKPAR,ALPHA,BETA,AKG1,AKG2
      COMPLEX(dp) C0,CC,COEF,Z1,Z2,Z3
      COMPLEX(dp) CT,ST,CF
      COMPLEX(dp) YLM((lmax+1)**2)
C     ------------------------------------------------------------------
      IF(LMAX>LMAXD)  GO TO 10
      AKG1=DREAL(GK(1))
      AKG2=DREAL(GK(2))
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
      IF(AKPAR>1.D-8) CF=DCMPLX(AKG1/AKPAR,AKG2/AKPAR)
      CALL SPHRM4(YLM,CT,ST,CF,LMAX)
      II=1
      CC=CONE
      DO L=1,LMAX
        CC=CC/CI
        COEF=C0*CC/SQRT(DFLOAT(L*(L+1)))
        DO M=-L,L
          II=II+1
          ALPHA=SQRT(DFLOAT((L-M)*(L+M+1)))/2.0_dp
          BETA =SQRT(DFLOAT((L+M)*(L-M+1)))/2.0_dp
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
          DLMH(1,II)=COEF*(BETA*CT*CF*Z2-DFLOAT(M)*ST*Z3
     &             +ALPHA*CT*CONJG(CF)*Z1)
          DLMH(2,II)=COEF*CI*(BETA*CF*Z2-ALPHA*CONJG(CF)*Z1)
          DLME(1,II)=COEF*CI*(BETA*CF*Z2-ALPHA*CONJG(CF)*Z1)
          DLME(2,II)=-COEF*(BETA*CT*CF*Z2-DFLOAT(M)*ST*Z3
     &             +ALPHA*CT*CONJG(CF)*Z1)
        end do
      end do
      RETURN
   10 WRITE(6,100) LMAX,LMAXD
      STOP
  100 FORMAT(//13X,'FROM DLMKG  :LMAX=',I5,' IS GREATER THAN DIMENSIONED
     * LMAXD=',I5)
  101 FORMAT(13X,'FATAL ERROR FROM DLMKG:'/3X,'GK(3) IS TOO SMALL.'
     & /3X,'GIVE A SMALL BUT NONZERO VALUE FOR "EPSILON"'/3X,
     & 'IN THE DATA STATEMENT OF THE MAIN PROGRAM.'
     & /3X,'THIS DEFINES A SMALL IMAGINARY PART'
     & /3X,'IN THE FREQUENCY OR WAVELENGTH VALUE.')
      END subroutine
C=======================================================================
C=======================================================================
C=======================================================================
      SUBROUTINE ELMGEN(ELM,NELMD,LMAX)

C     ------------------------------------------------------------------
C     ROUTINE TO TABULATE THE CLEBSCH-GORDON TYPE COEFFICIENTS ELM,  FOR
C     USE WITH THE SUBROUTINE XMAT. THE NON-ZERO ELM ARE TABULATED FIRST
C     FOR  L2,M2; AND L3,M3; ODD. THEN FOR L2,M2; AND L3,M3; EVEN, USING
C     THE SAME SCHEME AS THAT BY WHICH THEY ARE ACCESSED IN XMAT.
C     ------------------------------------------------------------------
C
C ..  SCALAR ARGUMENTS  ..
C
      INTEGER NELMD,LMAX
C
C ..  ARRAY ARGUMENTS  ..
C
      REAL(dp) ELM(NELMD)
C
C ..  LOCAL SCALARS  ..
C
      INTEGER K,II,LL,IL2,L2,M2,I2,IL3,L3,M3,I3,LA1,LB1,LA11,LB11,M1
      INTEGER L11,L1,L
      REAL(dp) FOURPI
C     ------------------------------------------------------------------
      FOURPI=4.0_dp*PI
      K=1
      II=0
   1  LL=LMAX+II
      DO 6 IL2=1,LL
      L2=IL2-II
      M2=-L2+1-II
      DO 6 I2=1,IL2
      DO 5 IL3=1,LL
      L3=IL3-II
      M3=-L3+1-II
      DO 5 I3=1,IL3
      LA1=MAX0(IABS(L2-L3),IABS(M2-M3))
      LB1=L2+L3
      LA11=LA1+1
      LB11=LB1+1
      M1=M2-M3
      DO 3 L11=LA11,LB11,2
      L1=L11-1
      L=(L2-L3-L1)/2+M2
      ELM(K)=((-1.0_dp)**L)*FOURPI*BLM(L1,M1,L3,M3,L2,-M2,LMAX)
   3  K=K+1
   5  M3=M3+2
   6  M2=M2+2
      IF(II)7,7,8
   7  II=1
      GOTO 1
   8  CONTINUE
      RETURN
      END subroutine
C=======================================================================
      SUBROUTINE LAT2D(A,B,RMAX,IMAX,ID,NTA,NTB,VECMOD)

C     --------------------------------------------------------------
C     GIVEN A TWO DIMENSIONAL BRAVAIS LATTICE WITH PRIMITIVE VECTORS
C     (A(1),A(2)) , (B(1),B(2)) , DEFINED SO THAT 'B' IS LONGER THAN
C     'A' AND THEIR SCALAR PRODUCT IS POSITIVE,THIS ROUTINE CALCULA-
C     TES THE 'IMAX' LATTICE VECTORS: NTA(I) * A + NTB(I) * B,HAVING
C     LENGTH 'VECMOD(I)' LESS THAN 'RMAX'.
C     --------------------------------------------------------------
C
C ..  SCALAR ARGUMENTS ..
C
      INTEGER IMAX,ID
      REAL(dp)RMAX
C
C ..  ARRAY ARGUMENTS ..
C
      INTEGER NTA(ID),NTB(ID)
      REAL(dp)A(2),B(2),VECMOD(ID)
C
C ..  LOCAL SCALARS ..
C
      INTEGER I,NA,NB,NA0,J,NMA,NMB,IORD
      REAL(dp)RMAX2,SP,AMOD2,BMOD2,DUM,VMOD2,VM
C
C ..  INTRINSIC FUNCTIONS ..
C
      INTRINSIC SQRT
C     ------------------------------------------------------------------
C
      RMAX2=RMAX*RMAX
C
C***  CHECK IF PRIMITIVE VECTORS HAVE POSITIVE SCALAR PRODUCT
C
      SP=A(1)*B(1)+A(2)*B(2)
      IF(SP<-1.D-06)  THEN
      B(1)=-B(1)
      B(2)=-B(2)
      SP=-SP
      WRITE(6,100) A(1),A(2),B(1),B(2)
                        END     IF
C
C***  CHECK IF 'B' IS LONGER THAN 'A'
C
      AMOD2=A(1)*A(1)+A(2)*A(2)
      BMOD2=B(1)*B(1)+B(2)*B(2)
      IF(BMOD2<AMOD2) THEN
      WRITE(6,101)
      DO 10 J=1,2
      DUM=A(J)
      A(J)=B(J)
   10 B(J)=DUM
      DUM=AMOD2
      AMOD2=BMOD2
      BMOD2=DUM
                      ENDIF
C
      I=0
      NB=0
    9 CONTINUE
      IF((NB*NB*BMOD2)>RMAX2)  GO TO 8
      NA=0
    7 CONTINUE
      VMOD2=NA*NA*AMOD2+NB*NB*BMOD2+2*NA*NB*SP
      IF(VMOD2>RMAX2)  GO TO 6
      I=I+1
      IF(I>ID)  GO TO 13
      NTA(I)=NA
      NTB(I)=NB
      VECMOD(I)=SQRT(VMOD2)
      IF(NA==0.AND.NB==0) GO TO 11
      I=I+1
      IF(I>ID) GO TO 13
      NTA(I)=-NA
      NTB(I)=-NB
      VECMOD(I)=SQRT(VMOD2)
   11 NA=NA+1
      GO TO 7
    6 CONTINUE
      NB=NB+1
      GO TO 9
    8 CONTINUE
C
      NA0=SP/AMOD2 + 1
      NB=1
    5 CONTINUE
      IF((NB*NB*(BMOD2-SP*SP/AMOD2))>RMAX2) GO TO 4
      NA=NA0
    3 CONTINUE
      VMOD2=NA*NA*AMOD2+NB*NB*BMOD2-2*NA*NB*SP
      IF(VMOD2>RMAX2) GO TO 2
      I=I+1
      IF(I>ID)  GO TO 13
      NTA(I)=NA
      NTB(I)=-NB
      VECMOD(I)=SQRT(VMOD2)
      I=I+1
      IF(I>ID)  GO TO 13
      NTA(I)=-NA
      NTB(I)=NB
      VECMOD(I)=SQRT(VMOD2)
      NA=NA+1
      GO TO 3
    2 CONTINUE
      NA=NA0-1
    1 CONTINUE
      VMOD2=NA*NA*AMOD2+NB*NB*BMOD2-2*NA*NB*SP
      IF(VMOD2>RMAX2.OR.NA<=0)  GO TO 12
      I=I+1
      IF(I>ID)  GO TO 13
      NTA(I)=NA
      NTB(I)=-NB
      VECMOD(I)=SQRT(VMOD2)
      I=I+1
      IF(I>ID) GO TO 13
      NTA(I)=-NA
      NTB(I)=NB
      VECMOD(I)=SQRT(VMOD2)
      NA=NA-1
      GO TO 1
   12 CONTINUE
      NB=NB+1
      GO TO 5
    4 CONTINUE
      IMAX=I
C
      DO 15 IORD=1,IMAX
      VM=VECMOD(IORD)
      DO 16 I=IMAX,IORD,-1
      IF(VECMOD(I)>VM)  GO TO 16
      VM=VECMOD(I)
      VECMOD(I)=VECMOD(IORD)
      VECMOD(IORD)=VM
      NMA=NTA(I)
      NTA(I)=NTA(IORD)
      NTA(IORD)=NMA
      NMB=NTB(I)
      NTB(I)=NTB(IORD)
      NTB(IORD)=NMB
   16 CONTINUE
   15 CONTINUE
C
      RETURN
   13 IMAX=I-1
      WRITE(6,102) IMAX
      DO 14 I=1,IMAX
      WRITE(6,103) I,NTA(I),A(1),A(2),NTB(I),B(1),B(2),VECMOD(I)
   14 CONTINUE
      STOP
C
  100 FORMAT(/13X,'NEW PRIMITIVE VECTORS DEFINED TO HAVE POSITIVE SCALAR
     & PRODUCT'/13X,'A=(',2E14.6,')'/13X,'B=(',2E14.6,')')
  101 FORMAT(/13X,'W A R N I N G ! !'/'INTERCHANGE PRIMITIVE VECTORS IN
     &CALL LAT2D'/)
  102 FORMAT(//33X,'FROM LAT2D: MAXIMUM NUMBER OF NEIGHBOURS=',I4,
     &'  EXCEEDED'//6X,'LATTICE POINTS FOUND (NON ORDERED)')
  103 FORMAT(I3,3X,I5,'*(',2E14.6,') +',I5,'*(',2E14.6,')',8X,E14.6)
C
      END subroutine
C=======================================================================
      SUBROUTINE PLW(KAPPA,GK,LMAX,AE,AH)

C     ------------------------------------------------------------------
C     THIS ROUTINE CALCULATES THE EXPANSION COEFFICIENTS 'AE,AH' OF AN
C     INCIDENT PLANE ELECTROMAGNETIC WAVE OF WAVE VECTOR  'KAPPA' WITH
C     COMPONENTS PARALLEL TO THE SURFACE EQUAL TO   '(GK(1),GK(2))'.
C     ------------------------------------------------------------------
C
C ..  PARAMETER STATEMENTS  ..
C
      INTEGER LMAXD,LMAX1D,LM1SQD
      PARAMETER(LMAXD=14,LMAX1D=LMAXD+1,LM1SQD=LMAX1D*LMAX1D)
C
C ..  SCALAR ARGUMENTS  ..
C
      INTEGER    LMAX
      COMPLEX(dp) KAPPA
C
C ..  ARRAY ARGUMENTS  ..
C
      COMPLEX(dp) AE(2,LM1SQD),AH(2,LM1SQD),GK(3)
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
      IF(LMAX>LMAXD)  GO TO 10
      AKG1=DREAL(GK(1))
      AKG2=DREAL(GK(2))
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
      CC=DCMPLX(FPI,0.0_dp)
      SIGNUS=-1.0_dp
      DO 1 L=1,LMAX
      CC=CC*CI
      A=DFLOAT(L*(L+1))
      CC1=CC/SQRT(A)
      DO 2 M=-L,L
      SIGNUS=-SIGNUS
      II=II+1
      IF(ABS(M+1)<=L)  then
      I=L*L+L-M
      Z1=CC1*SQRT(DFLOAT((L-M)*(L+M+1)))*YLM(I)/2.0_dp
                         else
      Z1=CZERO
                         end if
      IF(ABS(M-1)<=L)  then
      I=L*L+L-M+2
      Z2=CC1*SQRT(DFLOAT((L+M)*(L-M+1)))*YLM(I)/2.0_dp
                         else
      Z2=CZERO
                         end if
      I=L*L+L-M+1
      Z3=CC1*DFLOAT(M)*YLM(I)
      AE(1,II)= SIGNUS*CI*(CF*Z1-CONJG(CF)*Z2)
      AE(2,II)=-SIGNUS*(CT*CF*Z1+ST*Z3+CT*CONJG(CF)*Z2)
      AH(1,II)= SIGNUS*(CT*CF*Z1+ST*Z3+CT*CONJG(CF)*Z2)
      AH(2,II)= SIGNUS*CI*(CF*Z1-CONJG(CF)*Z2)
    2 CONTINUE
    1 CONTINUE
      RETURN
   10 WRITE(6,100) LMAX,LMAXD
      STOP
  100 FORMAT(//13X,'FROM PLW:  LMAX=',I5,'  IS GREATER THAN DIMENSIONED
     * LMAXD=',I5)
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
      UP=DFLOAT(2*LA+1)
      SIGNUS=-SIGNUS
      C=SIGNUS*C0
      B1=0.0_dp
      IF(ABS(MA+1)<=(LA-1)) B1=BLM(LA-1,MA+1,1,-1,LA,-MA,LMAX)
      B2=0.0_dp
      IF(ABS(MA-1)<=(LA-1)) B2=BLM(LA-1,MA-1,1, 1,LA,-MA,LMAX)
      U1=DFLOAT((LA+MA)*(LA-MA))
      U2=DFLOAT((2*LA-1)*(2*LA+1))
      B3=SQRT(U1/U2)
      ALPHA1=SQRT(DFLOAT((LA-MA)*(LA+MA+1)))/2.0_dp
      BETA1 =SQRT(DFLOAT((LA+MA)*(LA-MA+1)))/2.0_dp
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
      A=DFLOAT(LB*(LB+1)*LA*(LA+1))
      DOWN=SQRT(A)
      ALPHA2=SQRT(DFLOAT((LB-MB)*(LB+MB+1)))/2.0_dp
      BETA2 =SQRT(DFLOAT((LB+MB)*(LB-MB+1)))/2.0_dp
      LTT=LA+MA+LB+MB
          IF(MOD(LTT,2)/=0)           then
             IF(MOD((LA+MA),2)==0)       then
             Z1=CEVEN(LB,MB+1,LA-1,MA+1,XEVEN)
             Z2=CEVEN(LB,MB-1,LA-1,MA-1,XEVEN)
             Z3=CODD (LB,MB  ,LA-1,MA  ,XODD )
             Z1= C*ALPHA2*B1*Z1
             Z2=-C*BETA2* B2*Z2
             Z3=DFLOAT(MB)*B3*Z3
             OMEGA2=UP*(Z1+Z2+Z3)/DOWN
             XXMAT1(IA,IB)=-TH(LA+1)*OMEGA2
             XXMAT2(IA,IB)= TE(LA+1)*OMEGA2
                                           else
             Z1=CODD (LB,MB+1,LA-1,MA+1,XODD )
             Z2=CODD (LB,MB-1,LA-1,MA-1,XODD )
             Z3=CEVEN(LB,MB  ,LA-1,MA  ,XEVEN)
             Z1= C*ALPHA2*B1*Z1
             Z2=-C*BETA2* B2*Z2
             Z3=DFLOAT(MB)*B3*Z3
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
             Z3=DFLOAT(MA)*DFLOAT(MB)*Z3
             OMEGA1=(Z1+Z2+Z3)/DOWN
             XXMAT1(IA,IB)=-TH(LA+1)*OMEGA1
             XXMAT2(IA,IB)=-TE(LA+1)*OMEGA1
                                           else
             Z1=CEVEN(LB,MB-1,LA,MA-1,XEVEN)
             Z2=CEVEN(LB,MB+1,LA,MA+1,XEVEN)
             Z3=CODD (LB,MB  ,LA,MA  ,XODD )
             Z1=2.0_dp*BETA1 *BETA2 *Z1
             Z2=2.0_dp*ALPHA1*ALPHA2*Z2
             Z3=DFLOAT(MA)*DFLOAT(MB)*Z3
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
   2  FAC(I+1)=DFLOAT(I)*FAC(I)
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
   5  ALPHA=DCMPLX(AL,0.0_dp)
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
      AN1=DFLOAT(N1)
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
      GKK=DCMPLX(1.0_dp,0.0_dp)
      IF(AC-EMACH)11,11,10
  10  XPK=DCMPLX(AKPT(1)/AC,AKPT(2)/AC)
      GK=AC/KAPPA
      GKK=GPSQ/KAPSQ
  11  XPM(1)=DCMPLX(1.0_dp,0.0_dp)
      AGK(1)=DCMPLX(1.0_dp,0.0_dp)
      DO 12 I=2,LL2
      XPM(I)=XPM(I-1)*XPK
  12  AGK(I)=AGK(I-1)*GK
      CF=KAPPA/GP
      ZZ=-ALPHA*GKK
      CZ=SQRT(-ZZ)
      Z=-CI*CZ
      CX=EXP(-ZZ)
      GAM=RTPI*CERF(CZ,EMACH)
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
      AN1=DFLOAT(N1)
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
      XPK=DCMPLX(R(1)/AR,R(2)/AR)
      XPM(1)=DCMPLX(1.0_dp,0.0_dp)
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
      WW=CERF(-ZZ,EMACH)
      W=CERF(Z,EMACH)
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
      CF=DCMPLX(1.0_dp,0.0_dp)
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
      ACC=KAPPA*(CI*(XPA-CERF(RTA,EMACH))-RTAI)/XPA
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
C
C ..  INTRINSIC FUNCTIONS  ..
C
      INTRINSIC ABS
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
C
C ..  INTRINSIC FUNCTIONS  ..
C
      INTRINSIC ABS
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
      SUBROUTINE CBABK2(NM,N,LOW,IGH,SCALE,M,ZR,ZI)

C     ------------------------------------------------------------------
C     THIS SUBROUTINE FORMS THE EIGENVECTORS OF A COMPLEX GENERAL
C     MATRIX BY BACK TRANSFORMING THOSE OF THE CORRESPONDING
C     BALANCED MATRIX DETERMINED BY  CBAL.
C
C     ON INPUT--->
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT,
C
C        N IS THE ORDER OF THE MATRIX,
C
C        LOW AND IGH ARE INTEGERS DETERMINED BY  CBAL,
C
C        SCALE CONTAINS INFORMATION DETERMINING THE PERMUTATIONS
C          AND SCALING FACTORS USED BY  CBAL,
C
C        M IS THE NUMBER OF EIGENVECTORS TO BE BACK TRANSFORMED,
C
C        ZR AND ZI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE EIGENVECTORS TO BE
C          BACK TRANSFORMED IN THEIR FIRST M COLUMNS.
C
C     ON OUTPUT--->
C        ZR AND ZI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE TRANSFORMED EIGENVECTORS
C          IN THEIR FIRST M COLUMNS.
C     ------------------------------------------------------------------
C
C ..  SCALAR ARGUMENTS  ..
C
      INTEGER NM,N,LOW,IGH,M
C
C ..  ARRAY ARGUMENTS  ..
C
      REAL(dp) SCALE(N),ZR(NM,M),ZI(NM,M)
C
C ..  LOCAL SCALARS  ..
C
      INTEGER I,J,K,II
      REAL(dp)S
C     ------------------------------------------------------------------
C
      IF (M==0) GO TO 200
      IF (IGH==LOW) GO TO 120
C
      DO 110 I = LOW, IGH
         S = SCALE(I)
C     ********** LEFT HAND EIGENVECTORS ARE BACK TRANSFORMED
C                IF THE FOREGOING STATEMENT IS REPLACED BY
C                S=1.0/SCALE(I). **********
         DO 100 J = 1, M
            ZR(I,J) = ZR(I,J) * S
            ZI(I,J) = ZI(I,J) * S
  100    CONTINUE
C
  110 CONTINUE
C     ********** FOR I=LOW-1 STEP -1 UNTIL 1,
C                IGH+1 STEP 1 UNTIL N DO -- **********
  120 DO 140 II = 1, N
         I = II
         IF (I >= LOW .AND. I <= IGH) GO TO 140
         IF (I < LOW) I = LOW - II
         K = SCALE(I)
         IF (K == I) GO TO 140
C
         DO 130 J = 1, M
            S = ZR(I,J)
            ZR(I,J) = ZR(K,J)
            ZR(K,J) = S
            S = ZI(I,J)
            ZI(I,J) = ZI(K,J)
            ZI(K,J) = S
  130    CONTINUE
C
  140 CONTINUE
C
  200 RETURN
      END subroutine
C=======================================================================
      SUBROUTINE CBAL(NM,N,AR,AI,LOW,IGH,SCALE)

C     ------------------------------------------------------------------
C     THIS SUBROUTINE BALANCES A COMPLEX MATRIX AND ISOLATES
C     EIGENVALUES WHENEVER POSSIBLE.
C
C     ON INPUT--->
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT,
C
C        N IS THE ORDER OF THE MATRIX,
C
C        AR AND AI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE COMPLEX MATRIX TO BE BALANCED.
C
C     ON OUTPUT--->
C        AR AND AI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE BALANCED MATRIX,
C
C        LOW AND IGH ARE TWO INTEGERS SUCH THAT AR(I,J) AND AI(I,J)
C          ARE EQUAL TO ZERO IF
C           (1) I IS GREATER THAN J AND
C           (2) J=1,...,LOW-1 OR I=IGH+1,...,N,
C
C        SCALE CONTAINS INFORMATION DETERMINING THE
C           PERMUTATIONS AND SCALING FACTORS USED.
C
C     SUPPOSE THAT THE PRINCIPAL SUBMATRIX IN ROWS LOW THROUGH IGH
C     HAS BEEN BALANCED, THAT P(J) DENOTES THE INDEX INTERCHANGED
C     WITH J DURING THE PERMUTATION STEP, AND THAT THE ELEMENTS
C     OF THE DIAGONAL MATRIX USED ARE DENOTED BY D(I,J).  THEN
C        SCALE(J) = P(J),    FOR J = 1,...,LOW-1
C                 = D(J,J)       J = LOW,...,IGH
C                 = P(J)         J = IGH+1,...,N.
C     THE ORDER IN WHICH THE INTERCHANGES ARE MADE IS N TO IGH+1,
C     THEN 1 TO LOW-1.
C
C     NOTE THAT 1 IS RETURNED FOR IGH IF IGH IS ZERO FORMALLY.
C
C     ARITHMETIC IS REAL THROUGHOUT.
C     ------------------------------------------------------------------
C
C ..  SCALAR ARGUMENTS  ..
C
      INTEGER NM,N,LOW,IGH
C
C ..  ARRAY ARGUMENTS  ..
C
      REAL(dp) AR(NM,N),AI(NM,N),SCALE(N)
C
C ..  LOCAL SCALARS  ..
C
      INTEGER I,J,K,L,M,JJ,IEXC
      REAL(dp)C,F,G,R,S,B2,RADIX
      LOGICAL NOCONV
C     ------------------------------------------------------------------
C
C     ********** RADIX IS A MACHINE DEPENDENT PARAMETER SPECIFYING
C                THE BASE OF THE MACHINE FLOATING POINT REPRESENTATION.
C
      RADIX = 2.0_dp
C
      B2 = RADIX * RADIX
      K = 1
      L = N
      GO TO 100
C     ******** IN-LINE PROCEDURE FOR ROW AND COLUMN EXCHANGE ********
   20 SCALE(M) = J
      IF (J == M) GO TO 50
C
      DO 30 I = 1, L
         F = AR(I,J)
         AR(I,J) = AR(I,M)
         AR(I,M) = F
         F = AI(I,J)
         AI(I,J) = AI(I,M)
         AI(I,M) = F
   30 CONTINUE
C
      DO 40 I = K, N
         F = AR(J,I)
         AR(J,I) = AR(M,I)
         AR(M,I) = F
         F = AI(J,I)
         AI(J,I) = AI(M,I)
         AI(M,I) = F
   40 CONTINUE
C
   50 GO TO (80,130), IEXC
C     ********** SEARCH FOR ROWS ISOLATING AN EIGENVALUE
C                AND PUSH THEM DOWN **********
   80 IF (L == 1) GO TO 280
      L = L - 1
C     ********** FOR J=L STEP -1 UNTIL 1 DO -- **********
  100 DO 120 JJ = 1, L
         J = L + 1 - JJ
C
         DO 110 I = 1, L
            IF (I == J) GO TO 110
            IF (AR(J,I) /= 0.0_dp.OR. AI(J,I) /=0.0_dp) GO TO 120
  110    CONTINUE
C
         M = L
         IEXC = 1
         GO TO 20
  120 CONTINUE
C
      GO TO 140
C     ********** SEARCH FOR COLUMNS ISOLATING AN EIGENVALUE
C                AND PUSH THEM LEFT **********
  130 K = K + 1
C
  140 DO 170 J = K, L
C
         DO 150 I = K, L
            IF (I == J) GO TO 150
            IF (AR(I,J) /= 0.0_dp .OR. AI(I,J) /= 0.0_dp) GO TO 170
  150    CONTINUE
C
         M = K
         IEXC = 2
         GO TO 20
  170 CONTINUE
C     ********** NOW BALANCE THE SUBMATRIX IN ROWS K TO L **********
      DO 180 I = K, L
  180 SCALE(I) = 1.0_dp
C     ********** ITERATIVE LOOP FOR NORM REDUCTION **********
  190 NOCONV = .FALSE.
C
      DO 270 I = K, L
         C = 0.0_dp
         R = 0.0_dp
C
         DO 200 J = K, L
            IF (J == I) GO TO 200
            C = C + ABS(AR(J,I)) + ABS(AI(J,I))
            R = R + ABS(AR(I,J)) + ABS(AI(I,J))
  200    CONTINUE
C     ********** GUARD AGAINST ZERO C OR R DUE TO UNDERFLOW **********
         IF (C == 0.0_dp .OR. R == 0.0_dp) GO TO 270
         G = R / RADIX
         F = 1.0_dp
         S = C + R
  210    IF (C >= G) GO TO 220
         F = F * RADIX
         C = C * B2
         GO TO 210
  220    G = R * RADIX
  230    IF (C < G) GO TO 240
         F = F / RADIX
         C = C / B2
         GO TO 230
C     ********** NOW BALANCE **********
  240    IF ((C + R) / F >= 0.95D0 * S) GO TO 270
         G = 1.0_dp / F
         SCALE(I) = SCALE(I) * F
         NOCONV = .TRUE.
C
         DO 250 J = K, N
            AR(I,J) = AR(I,J) * G
            AI(I,J) = AI(I,J) * G
  250    CONTINUE
C
         DO 260 J = 1, L
            AR(J,I) = AR(J,I) * F
            AI(J,I) = AI(J,I) * F
  260    CONTINUE
C
  270 CONTINUE
C
      IF (NOCONV) GO TO 190
C
  280 LOW = K
      IGH = L
      RETURN
      END subroutine
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
      CALL CBAL(NDIM,N,AR,AI,LOW,IGH,SCALE)
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
      SUBROUTINE COMHES(NM,N,LOW,IGH,AR,AI,INT)

C     ------------------------------------------------------------------
C     GIVEN A  COMPLEX  GENERAL  MATRIX, THIS  SUBROUTINE  REDUCES  A
C     SUBMATRIX SITUATED IN ROWS AND COLUMNS LOW THROUGH IGH TO UPPER
C     HESSENBERG FORM BY STABILIZED ELEMENTARY SIMILARITY TRANSFORMS.
C
C     ON INPUT--->
C        NM       MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C                 ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C                 DIMENSION STATEMENT
C        N        IS THE ORDER OF THE MATRIX
C        LOW,IGH  ARE INTEGERS DETERMINED BY THE BALANCING SUBROUTINE
C                 CBAL. IF  CBAL  HAS NOT BEEN USED, SET LOW=1, IGH=N
C        AR,AI    CONTAIN THE REAL AND IMAGINARY PARTS, RESPECTIVELY,
C                 OF THE COMPLEX INPUT MATRIX.
C
C     ON OUTPUT--->
C        AR,AI    CONTAIN THE REAL AND IMAGINARY PARTS, RESPECTIVELY,
C                 OF THE HESSENBERG MATRIX.THE MULTIPLIERS WHICH WERE
C                 USED IN THE  REDUCTION  ARE STORED IN THE REMAINING
C                 TRIANGLES UNDER THE HESSENBERG MATRIX,
C        INT      CONTAINS INFORMATION ON THE ROWS AND COLUMNS INTER-
C                 CHANGED IN THE REDUCTION. ONLY ELEMENTS LOW THROUGH
C                 IGH ARE USED.
C
C     ARITHMETIC IS REAL EXCEPT FOR THE REPLACEMENT OF THE ALGOL
C     PROCEDURE CDIV BY COMPLEX DIVISION USING SUBROUTINE CMPLX.
C     ------------------------------------------------------------------
C
C ..  SCALAR ARGUMENTS  ..
C
      INTEGER NM,N,LOW,IGH
C
C ..  ARRAY ARGUMENTS  ..
C
      INTEGER INT(IGH)
      REAL(dp)AR(NM,N),AI(NM,N)
C
C ..  LOCAL SCALARS  ..
C
      INTEGER    I,J,M,LA,KP1,MM1,MP1
      REAL(dp)   XR,XI,YR,YI
      COMPLEX(dp) Z3
C     ------------------------------------------------------------------
C
      LA = IGH - 1
      KP1 = LOW + 1
      IF (LA < KP1) GO TO 200
C
      DO 180 M = KP1, LA
         MM1 = M - 1
         XR = 0.0_dp
         XI = 0.0_dp
         I = M
C
         DO 100 J = M, IGH
            IF (ABS(AR(J,MM1)) + ABS(AI(J,MM1))
     X         <= ABS(XR) + ABS(XI)) GO TO 100
            XR = AR(J,MM1)
            XI = AI(J,MM1)
            I = J
  100    CONTINUE
C
         INT(M) = I
         IF (I == M) GO TO 130
C     ********** INTERCHANGE ROWS AND COLUMNS OF AR AND AI **********
         DO 110 J = MM1, N
            YR = AR(I,J)
            AR(I,J) = AR(M,J)
            AR(M,J) = YR
            YI = AI(I,J)
            AI(I,J) = AI(M,J)
            AI(M,J) = YI
  110    CONTINUE
C
         DO 120 J = 1, IGH
            YR = AR(J,I)
            AR(J,I) = AR(J,M)
            AR(J,M) = YR
            YI = AI(J,I)
            AI(J,I) = AI(J,M)
            AI(J,M) = YI
  120    CONTINUE
C     ********** END INTERCHANGE **********
  130    IF (XR == 0.0_dp .AND. XI == 0.0_dp) GO TO 180
         MP1 = M + 1
C
         DO 160 I = MP1, IGH
            YR = AR(I,MM1)
            YI = AI(I,MM1)
            IF (YR == 0.0_dp .AND. YI == 0.0_dp) GO TO 160
            Z3 = DCMPLX(YR,YI) / DCMPLX(XR,XI)
            YR = DREAL(Z3)
            YI = DIMAG (Z3)
            AR(I,MM1) = YR
            AI(I,MM1) = YI
C
            DO 140 J = M, N
               AR(I,J) = AR(I,J) - YR * AR(M,J) + YI * AI(M,J)
               AI(I,J) = AI(I,J) - YR * AI(M,J) - YI * AR(M,J)
  140       CONTINUE
C
            DO 150 J = 1, IGH
               AR(J,M) = AR(J,M) + YR * AR(J,I) - YI * AI(J,I)
               AI(J,M) = AI(J,M) + YR * AI(J,I) + YI * AR(J,I)
  150       CONTINUE
C
  160    CONTINUE
C
  180 CONTINUE
C
  200 RETURN
      END subroutine
C=======================================================================
      SUBROUTINE COMLR2(NM,N,LOW,IGH,INT,HR,HI,WR,WI,ZR,ZI,IERR)

C     ------------------------------------------------------------------
C     THIS SUBROUTINE FINDS  THE  EIGENVALUES AND  EIGENVECTORS  OF  A
C     COMPLEX UPPER HESSENBERG  MATRIX BY THE MODIFIED  LR METHOD. THE
C     EIGENVECTORS  OF A COMPLEX  GENERAL MATRIX  CAN ALSO BE FOUND IF
C     COMHES HAS BEEN USED TO REDUCE THIS GENERAL MATRIX TO HESSENBERG
C     FORM.
C
C     ON INPUT--->
C        NM      MUST  BE SET TO THE ROW DIMENSION  OF TWO-DIMENSIONAL
C                ARRAY  PARAMETERS AS  DECLARED IN THE CALLING PROGRAM
C                DIMENSION STATEMENT
C        N       IS THE ORDER OF THE MATRIX
C        LOW,IGH ARE INTEGERS DETERMINED BY THE  BALANCING  SUBROUTINE
C                CBAL.  IF  CBAL  HAS NOT BEEN USED,  SET LOW=1, IGH=N
C        INT     CONTAINS INFORMATION ON THE ROWS AND  COLUMNS  INTER-
C                CHANGED IN THE REDUCTION BY COMHES,IF PERFORMED. ONLY
C                ELEMENTS LOW THROUGH IGH ARE USED.IF THE EIGENVECTORS
C                OF THE HESSENBERG MATRIX ARE DESIRED,SET INT(J)=J FOR
C                THESE ELEMENTS
C        HR,HI   CONTAIN THE REAL AND IMAGINARY PARTS, RESPECTIVELY,OF
C                THE COMPLEX UPPER HESSENBERG MATRIX. THEIR LOWER TRI-
C                ANGLES  BELOW THE SUBDIAGONAL CONTAIN THE MULTIPLIERS
C                WHICH   WERE  USED  IN  THE  REDUCTION BY  COMHES, IF
C                PERFORMED.  IF  THE  EIGENVECTORS  OF  THE HESSENBERG
C                MATRIX ARE DESIRED,THESE ELEMENTS MUST BE SET TO ZERO
C
C      ON OUTPUT--->
C                THE   UPPER HESSENBERG PORTIONS OF HR AND HI HAVE BEEN
C                DESTROYED, BUT  THE LOCATION HR(1,1) CONTAINS THE NORM
C                OF THE TRIANGULARIZED MATRIX,
C        WR,WI   CONTAIN THE REAL AND IMAGINARY PARTS, RESPECTIVELY, OF
C                THE   EIGENVALUES.  IF  AN  ERROR  EXIT  IS  MADE, THE
C                EIGENVALUES SHOULD BE CORRECT FOR INDICES IERR+1,...,N
C        ZR,ZI   CONTAIN THE REAL AND IMAGINARY PARTS, RESPECTIVELY, OF
C                THE EIGENVECTORS.THE EIGENVECTORS ARE UNNORMALIZED. IF
C                AN ERROR EXIT IS  MADE, NONE OF THE  EIGENVECTORS  HAS
C                BEEN FOUND
C        IERR    IS SET TO  ZERO FOR NORMAL RETURN,
C          J     IF THE J-TH  EIGENVALUE HAS NOT BEEN  DETERMINED AFTER
C                30 ITERATIONS.
C
C     ARITHMETIC  IS  REAL  EXCEPT  FOR THE  REPLACEMENT  OF  THE ALGOL
C     PROCEDURE CDIV BY  COMPLEX DIVISION AND  USE OF  THE  SUBROUTINES
C     CSQRT AND CMPLX IN COMPUTING COMPLEX SQUARE ROOTS.
C     ------------------------------------------------------------------
C
C ..  SCALAR ARGUMENTS  ..
C
      INTEGER NM,N,LOW,IGH,IERR
C
C ..  ARRAY ARGUMENTS  ..
C
      INTEGER INT(IGH)
      REAL(dp) HR(NM,N),HI(NM,N),WR(N),WI(N),ZR(NM,N),ZI(NM,N)
C
C ..  LOCAL SCALARS  ..
C
      INTEGER    I,J,K,L,M,EN,II,JJ,LL,MM,NN,IM1,IP1,ITS,MP1,ENM1,IEND
      REAL(dp)   SI,SR,TI,TR,XI,XR,YI,YR,ZZI,ZZR,NORM,MACHEP
      COMPLEX(dp) Z3
C     ------------------------------------------------------------------
C
C     ********** MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING
C                THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.
C
      MACHEP = 2.0_dp**(-47)
C
      IERR = 0
C     ********** INITIALIZE EIGENVECTOR MATRIX **********
      DO 100 I = 1, N
C
         DO 100 J = 1, N
            ZR(I,J) = 0.0_dp
            ZI(I,J) = 0.0_dp
            IF (I == J) ZR(I,J) = 1.0_dp
  100 CONTINUE
C     ********** FORM THE MATRIX OF ACCUMULATED TRANSFORMATIONS
C                FROM THE INFORMATION LEFT BY COMHES **********
      IEND = IGH - LOW - 1
      IF (IEND <= 0) GO TO 180
C     ********** FOR I=IGH-1 STEP -1 UNTIL LOW+1 DO -- **********
      DO 160 II = 1, IEND
         I = IGH - II
         IP1 = I + 1
C
         DO 120 K = IP1, IGH
            ZR(K,I) = HR(K,I-1)
            ZI(K,I) = HI(K,I-1)
  120    CONTINUE
C
         J = INT(I)
         IF (I == J) GO TO 160
C
         DO 140 K = I, IGH
            ZR(I,K) = ZR(J,K)
            ZI(I,K) = ZI(J,K)
            ZR(J,K) = 0.0_dp
            ZI(J,K) = 0.0_dp
  140    CONTINUE
C
         ZR(J,I) = 1.0_dp
  160 CONTINUE
C     ********** STORE ROOTS ISOLATED BY CBAL **********
  180 DO 200 I = 1, N
         IF (I >= LOW .AND. I <= IGH) GO TO 200
         WR(I) = HR(I,I)
         WI(I) = HI(I,I)
  200 CONTINUE
C
      EN = IGH
      TR = 0.0_dp
      TI = 0.0_dp
C     ********** SEARCH FOR NEXT EIGENVALUE **********
  220 IF (EN < LOW) GO TO 680
      ITS = 0
      ENM1 = EN - 1
C     ********** LOOK FOR SINGLE SMALL SUB-DIAGONAL ELEMENT
C                FOR L=EN STEP -1 UNTIL LOW DO -- **********
  240 DO 260 LL = LOW, EN
         L = EN + LOW - LL
         IF (L == LOW) GO TO 300
         IF (ABS(HR(L,L-1)) + ABS(HI(L,L-1)) <=
     X      MACHEP * (ABS(HR(L-1,L-1)) + ABS(HI(L-1,L-1))
     X             + ABS(HR(L,L)) + ABS(HI(L,L)))) GO TO 300
  260 CONTINUE
C     ********** FORM SHIFT **********
  300 IF (L == EN) GO TO 660
      IF (ITS == 30) GO TO 1000
      IF (ITS == 10 .OR. ITS == 20) GO TO 320
      SR = HR(EN,EN)
      SI = HI(EN,EN)
      XR = HR(ENM1,EN) * HR(EN,ENM1) - HI(ENM1,EN) * HI(EN,ENM1)
      XI = HR(ENM1,EN) * HI(EN,ENM1) + HI(ENM1,EN) * HR(EN,ENM1)
      IF (XR == 0.0_dp .AND. XI == 0.0_dp) GO TO 340
      YR = (HR(ENM1,ENM1) - SR) / 2.0_dp
      YI = (HI(ENM1,ENM1) - SI) / 2.0_dp
      Z3 = SQRT(DCMPLX(YR**2-YI**2+XR,2.0_dp*YR*YI+XI))
      ZZR = DREAL(Z3)
      ZZI = DIMAG(Z3)
      IF (YR * ZZR + YI * ZZI >= 0.0_dp) GO TO 310
      ZZR = -ZZR
      ZZI = -ZZI
  310 Z3 = DCMPLX(XR,XI) / DCMPLX(YR+ZZR,YI+ZZI)
      SR = SR - DREAL(Z3)
      SI = SI - DIMAG(Z3)
      GO TO 340
C     ********** FORM EXCEPTIONAL SHIFT **********
  320 SR = ABS(HR(EN,ENM1)) + ABS(HR(ENM1,EN-2))
      SI = ABS(HI(EN,ENM1)) + ABS(HI(ENM1,EN-2))
C
  340 DO 360 I = LOW, EN
         HR(I,I) = HR(I,I) - SR
         HI(I,I) = HI(I,I) - SI
  360 CONTINUE
C
      TR = TR + SR
      TI = TI + SI
      ITS = ITS + 1
C     ********** LOOK FOR TWO CONSECUTIVE SMALL
C                SUB-DIAGONAL ELEMENTS **********
      XR = ABS(HR(ENM1,ENM1)) + ABS(HI(ENM1,ENM1))
      YR = ABS(HR(EN,ENM1)) + ABS(HI(EN,ENM1))
      ZZR = ABS(HR(EN,EN)) + ABS(HI(EN,EN))
C     ********** FOR M=EN-1 STEP -1 UNTIL L DO -- **********
      DO 380 MM = L, ENM1
         M = ENM1 + L - MM
         IF (M == L) GO TO 420
         YI = YR
         YR = ABS(HR(M,M-1)) + ABS(HI(M,M-1))
         XI = ZZR
         ZZR = XR
         XR = ABS(HR(M-1,M-1)) + ABS(HI(M-1,M-1))
         IF (YR <= MACHEP * ZZR / YI * (ZZR + XR + XI)) GO TO 420
  380 CONTINUE
C     ********** TRIANGULAR DECOMPOSITION H=L*R **********
  420 MP1 = M + 1
C
      DO 520 I = MP1, EN
         IM1 = I - 1
         XR = HR(IM1,IM1)
         XI = HI(IM1,IM1)
         YR = HR(I,IM1)
         YI = HI(I,IM1)
         IF (ABS(XR) + ABS(XI) >= ABS(YR) + ABS(YI)) GO TO 460
C     ********** INTERCHANGE ROWS OF HR AND HI **********
         DO 440 J = IM1, N
            ZZR = HR(IM1,J)
            HR(IM1,J) = HR(I,J)
            HR(I,J) = ZZR
            ZZI = HI(IM1,J)
            HI(IM1,J) = HI(I,J)
            HI(I,J) = ZZI
  440    CONTINUE
C
         Z3 = DCMPLX(XR,XI) / DCMPLX(YR,YI)
         WR(I) = 1.0_dp
         GO TO 480
  460    Z3 = DCMPLX(YR,YI) / DCMPLX(XR,XI)
         WR(I) = -1.0_dp
  480    ZZR = DREAL(Z3)
         ZZI = DIMAG(Z3)
         HR(I,IM1) = ZZR
         HI(I,IM1) = ZZI
C
         DO 500 J = I, N
            HR(I,J) = HR(I,J) - ZZR * HR(IM1,J) + ZZI * HI(IM1,J)
            HI(I,J) = HI(I,J) - ZZR * HI(IM1,J) - ZZI * HR(IM1,J)
  500    CONTINUE
C
  520 CONTINUE
C     ********** COMPOSITION R*L=H **********
      DO 640 J = MP1, EN
         XR = HR(J,J-1)
         XI = HI(J,J-1)
         HR(J,J-1) = 0.0_dp
         HI(J,J-1) = 0.0_dp
C     ********** INTERCHANGE COLUMNS OF HR, HI, ZR, AND ZI,
C                IF NECESSARY **********
         IF (WR(J) <= 0.0_dp) GO TO 580
C
         DO 540 I = 1, J
            ZZR = HR(I,J-1)
            HR(I,J-1) = HR(I,J)
            HR(I,J) = ZZR
            ZZI = HI(I,J-1)
            HI(I,J-1) = HI(I,J)
            HI(I,J) = ZZI
  540    CONTINUE
C
         DO 560 I = LOW, IGH
            ZZR = ZR(I,J-1)
            ZR(I,J-1) = ZR(I,J)
            ZR(I,J) = ZZR
            ZZI = ZI(I,J-1)
            ZI(I,J-1) = ZI(I,J)
            ZI(I,J) = ZZI
  560    CONTINUE
C
  580    DO 600 I = 1, J
            HR(I,J-1) = HR(I,J-1) + XR * HR(I,J) - XI * HI(I,J)
            HI(I,J-1) = HI(I,J-1) + XR * HI(I,J) + XI * HR(I,J)
  600    CONTINUE
C     ********** ACCUMULATE TRANSFORMATIONS **********
         DO 620 I = LOW, IGH
            ZR(I,J-1) = ZR(I,J-1) + XR * ZR(I,J) - XI * ZI(I,J)
            ZI(I,J-1) = ZI(I,J-1) + XR * ZI(I,J) + XI * ZR(I,J)
  620    CONTINUE
C
  640 CONTINUE
C
      GO TO 240
C     ********** A ROOT FOUND **********
  660 HR(EN,EN) = HR(EN,EN) + TR
      WR(EN) = HR(EN,EN)
      HI(EN,EN) = HI(EN,EN) + TI
      WI(EN) = HI(EN,EN)
      EN = ENM1
      GO TO 220
C     ********** ALL ROOTS FOUND.  BACKSUBSTITUTE TO FIND
C                VECTORS OF UPPER TRIANGULAR FORM **********
  680 NORM = 0.0_dp
C
      DO 720 I = 1, N
C
         DO 720 J = I, N
            NORM = NORM + ABS(HR(I,J)) + ABS(HI(I,J))
  720 CONTINUE
C
      HR(1,1) = NORM
      IF (N == 1 .OR. NORM == 0.0_dp) GO TO 1001
C     ********** FOR EN=N STEP -1 UNTIL 2 DO -- **********
      DO 800 NN = 2, N
         EN = N + 2 - NN
         XR = WR(EN)
         XI = WI(EN)
         ENM1 = EN - 1
C     ********** FOR I=EN-1 STEP -1 UNTIL 1 DO -- **********
         DO 780 II = 1, ENM1
            I = EN - II
            ZZR = HR(I,EN)
            ZZI = HI(I,EN)
            IF (I == ENM1) GO TO 760
            IP1 = I + 1
C
            DO 740 J = IP1, ENM1
               ZZR = ZZR + HR(I,J) * HR(J,EN) - HI(I,J) * HI(J,EN)
               ZZI = ZZI + HR(I,J) * HI(J,EN) + HI(I,J) * HR(J,EN)
  740       CONTINUE
C
  760       YR = XR - WR(I)
            YI = XI - WI(I)
            IF (YR == 0.0_dp .AND. YI == 0.0_dp) YR = MACHEP * NORM
            Z3 = DCMPLX(ZZR,ZZI) / DCMPLX(YR,YI)
            HR(I,EN) = DREAL(Z3)
            HI(I,EN) = DIMAG(Z3)
  780    CONTINUE
C
  800 CONTINUE
C     ********** END BACKSUBSTITUTION **********
      ENM1 = N - 1
C     ********** VECTORS OF ISOLATED ROOTS **********
      DO 840 I = 1, ENM1
         IF (I >= LOW .AND. I <= IGH) GO TO 840
         IP1 = I + 1
C
         DO 820 J = IP1, N
            ZR(I,J) = HR(I,J)
            ZI(I,J) = HI(I,J)
  820    CONTINUE
C
  840 CONTINUE
C     ********** MULTIPLY BY TRANSFORMATION MATRIX TO GIVE
C                VECTORS OF ORIGINAL FULL MATRIX.
C                FOR J=N STEP -1 UNTIL LOW+1 DO -- **********
      DO 880 JJ = LOW, ENM1
         J = N + LOW - JJ
         M = MIN0(J-1,IGH)
C
         DO 880 I = LOW, IGH
            ZZR = ZR(I,J)
            ZZI = ZI(I,J)
C
            DO 860 K = LOW, M
               ZZR = ZZR + ZR(I,K) * HR(K,J) - ZI(I,K) * HI(K,J)
               ZZI = ZZI + ZR(I,K) * HI(K,J) + ZI(I,K) * HR(K,J)
  860       CONTINUE
C
            ZR(I,J) = ZZR
            ZI(I,J) = ZZI
  880 CONTINUE
C
      GO TO 1001
C     ********** SET ERROR -- NO CONVERGENCE TO AN
C                EIGENVALUE AFTER 30 ITERATIONS **********
 1000 IERR = EN
 1001 RETURN
      END subroutine
C=======================================================================
C=======================================================================
      COMPLEX(dp) FUNCTION CERF(Z,EMACH)

C     ------------------------------------------------------------------
C     cerf,GIVEN COMPLEX ARGUMENT Z,PROVIDES THE COMPLEX ERROR FUNCTION:
C     W(Z)=EXP(-Z**2)*(1.0-ERF(-I*Z))
C     THE  EVALUATION  ALWAYS  TAKES   PLACE  IN  THE  FIRST   QUADRANT.
C     ONE  OF  THREE METHODS  IS  EXPLOYED  DEPENDING ON THE SIZE OF THE
C     ARGUMENT (A POWER SERIES,A RECURRENCE BASED ON CONTINUED FRACTIONS
C     THEORY, OR AN ASYMPTOTIC SERIES). EMACH IS THE MACHINE ACCURACY
C     ------------------------------------------------------------------
C
C ..  SCALAR ARGUMENTS  ..
C
      REAL(dp), intent(in)::     EMACH
      COMPLEX(dp), intent(in):: Z
C
C ..  LOCAL SCALARS  ..
C
      INTEGER    NN,N
      REAL(dp)   ABSZ,ABTERM,API,EPS,FACT,FACTD,FACTN
      REAL(dp)   Q,RTPI,TEST,X,Y,YY
      COMPLEX(dp) ZZ,SUM,ZZS,XZZS,CER!, erf
      COMPLEX(dp) H1,H2,H3,U1,U2,U3,TERM1,TERM2
C     ------------------------------------------------------------------
C
!     erf = erf_pop(z)
!     return
      EPS=5.0_dp*EMACH
      API=1.0_dp/PI
      IF(ABS(Z))2,1,2
   1  cerf=CONE
      GOTO 29
C
C     THE ARGUMENT IS TRANSLATED TO THE FIRST QUADRANT FROM
C     THE NN_TH QUADRANT, BEFORE THE METHOD FOR THE FUNCTION
C     EVALUATION IS CHOSEN
C
   2  X=DREAL(Z)
      Y=DIMAG(Z)
      YY=Y
      IF(Y)6,3,3
   3  IF(X)5,4,4
   4  ZZ=Z
      NN=1
      GOTO 9
   5  ZZ=DCMPLX(-X,Y)
      NN=2
      GOTO 9
   6  YY=-Y
      IF(X)7,8,8
   7  ZZ=-Z
      NN=3
      GOTO 9
   8  ZZ=DCMPLX(X,-Y)
      NN=4
   9  ZZS=ZZ*ZZ
      XZZS=EXP(-ZZS)
      ABSZ=ABS(ZZ)
      IF(ABSZ-10.0_dp)10,10,23
  10  IF(YY-1.0_dp)11,12,12
  11  IF(ABSZ-4.0_dp)13,18,18
  12  IF(ABSZ-1.0_dp)13,18,18
C
C     POWER SERIES(SEE ABRAMOWITZ AND STEGUN HANDBOOK OF
C     MATHEMATICAL FUNCTIONS, P297)
C
  13  Q=1.0_dp
      FACTN=-1.0_dp
      FACTD=1.0_dp
      TERM1=ZZ
      SUM=ZZ
  14  DO 15 N=1,5
      FACTN=FACTN+2.0_dp
      FACTD=FACTD+2.0_dp
      FACT=FACTN/(Q*FACTD)
      TERM1=FACT*ZZS*TERM1
      SUM=SUM+TERM1
  15  Q=Q+1.0_dp
      ABTERM=ABS(TERM1)
      IF(ABTERM-EPS)17,16,16
  16  IF(Q-100.0_dp)14,17,17
  17  FACT=2.0_dp*SQRT(API)
      SUM=FACT*CI*SUM
      CER=XZZS+XZZS*SUM
      GOTO 24
C
C     CONTINUED FRACTION THEORY(W(Z) IS RELATED TO THE LIMITING
C     VALUE OF U(N,Z)/H(N,Z), WHERE U AND H OBEY THE SAME
C     RECURRENCE RELATION IN N. SEE FADDEEVA AND TERENTIEV
C     (TABLES OF VALUES OF W(Z) FOR COMPLEX ARGUMENTS,PERGAMON
C       N.Y. 1961)
C
  18  TERM2=DCMPLX(1.D6,0.0_dp)
      Q=1.0_dp
      H1=CONE
      H2=2.0_dp*ZZ
      U1=CZERO
      RTPI=2.0_dp*SQRT(PI)
      U2=DCMPLX(RTPI,0.0_dp)
  19  TERM1=TERM2
      DO 20 N=1,5
      H3=H2*ZZ-Q*H1
      U3=U2*ZZ-Q*U1
      H1=H2
      H2=2.0_dp*H3
      U1=U2
      U2=2.0_dp*U3
  20  Q=Q+1.0_dp
      TERM2=U3/H3
      TEST=ABS((TERM2-TERM1)/TERM1)
      IF(TEST-EPS)22,21,21
  21  IF(Q-60.0_dp)19,19,13
  22  CER=API*CI*TERM2
      GOTO 24
C
C     ASYMPTOTIC SERIES: SEE ABRAMOWITZ AND STEGUN, P328
C
  23  CER=0.5124242D0/(ZZS-0.2752551D0)+0.05176536D0/(ZZS-2.724745D0)
      CER=CI*ZZ*CER
C
C     SYMMETRY RELATIONS ARE NOW USED TO TRANSFORM THE FUNCTION
C     BACK TO QUADRANT NN
C
  24  GOTO(28,26,27,25),NN
  25  CER=2.0_dp*XZZS-CER
  26  cerf=CONJG(CER)
      GOTO 29
  27  cerf=2.0_dp*XZZS-CER
      GOTO 29
  28  cerf=CER
  29  RETURN

      END function
C=======================================================================
      REAL(dp) FUNCTION BLM(L1,M1,L2,M2,L3,M3,LMAX)

C-----------------------------------------------------------------------
C     FUNCTION BLM  PROVIDES  THE  INTEGRAL  OF  THE  PRODUCT  OF THREE
C     SPHERICAL HARMONICS,EACH OF WHICH CAN BE EXPRESSED AS A PREFACTOR
C     TIMES  A  LEGENDRE  FUNCTION. THE  THREE  PREFACTORS  ARE  LUMPED
C     TOGETHER AS  FACTOR 'C'; AND   THE INTEGRAL OF THE THREE LEGENDRE
C     FUNCTIONS FOLLOWS GAUNT SUMMATION SCHEME SET OUT BY SLATER(ATOMIC
C     STRUCTURE, VOL1, 309,310
C-----------------------------------------------------------------------
C
C ..  PARAMETER STATEMENTS  ..
C
      INTEGER LMAXD,LMAX4D
      PARAMETER (LMAXD=14,LMAX4D=4*LMAXD+2)
C
C ..  SCALAR ARGUMENTS  ..
C
      INTEGER L1,M1,L2,M2,L3,M3,LMAX
C
C ..  LOCAL SCALARS  ..
C
      INTEGER I,IA1,IA2,IA3,IA4,IA5,IA6,IA7,IA8,IA9,IB1,IB2,IB3,IB4
      INTEGER IB5,IC,IC1,IC2,IC3,IC4,IC5,IC6,IS,IT,IT1,IT2,NL1,NL2
      INTEGER NL3,NM1,NM2,NM3,NTEMP,NN
      REAL(dp) SIGN,A,AD,AN,B,BD,BN,C,CD,CN
C
C ..  LOCAL ARRAYS  ..
C
      REAL(dp)FAC(LMAX4D)
C-----------------------------------------------------------------------
      FAC(1)=1.0_dp
      NN=4*LMAX+1
      DO 1 I=1,NN
   1  FAC(I+1)=DFLOAT(I)*FAC(I)
      IF(M1+M2+M3)8,21,8
  21  IF(L1-LMAX-LMAX)2,2,19
   2  IF(L2-LMAX)3,3,19
   3  IF(L3-LMAX)4,4,19
   4  IF(L1-IABS(M1))19,5,5
   5  IF(L2-IABS(M2))19,6,6
   6  IF(L3-IABS(M3))19,7,7
   7  IF(MOD  (L1+L2+L3,2))8,9,8
   8  BLM=0.0_dp
      RETURN
   9  NL1=L1
      NL2=L2
      NL3=L3
      NM1=IABS(M1)
      NM2=IABS(M2)
      NM3=IABS(M3)
      IC=(NM1+NM2+NM3)/2
      IF(MAX0(NM1,NM2,NM3)-NM1)13,13,10
  10  IF(MAX0(NM2,NM3)-NM2)11,11,12
  11  NL1=L2
      NL2=L1
      NM1=NM2
      NM2=IABS(M1)
      GOTO 13
  12  NL1=L3
      NL3=L1
      NM1=NM3
      NM3=IABS(M1)
  13  IF(NL2-NL3)14,15,15
  14  NTEMP=NL2
      NL2=NL3
      NL3=NTEMP
      NTEMP=NM2
      NM2=NM3
      NM3=NTEMP
  15  IF(NL3-IABS(NL2-NL1))16,17,17
  16  BLM=0.0_dp
      RETURN
C
C     CALCULATION OF FACTOR  'A'.
C
  17  IS=(NL1+NL2+NL3)/2
      IA1=IS-NL2-NM3
      IA2=NL2+NM2
      IA3=NL2-NM2
      IA4=NL3+NM3
      IA5=NL1+NL2-NL3
      IA6=IS-NL1
      IA7=IS-NL2
      IA8=IS-NL3
      IA9=NL1+NL2+NL3+1
      AN=((-1.0_dp)**IA1)*FAC(IA2+1)*FAC(IA4+1)*FAC(IA5+1)*FAC(IS+1)
      AD=FAC(IA3+1)*FAC(IA6+1)*FAC(IA7+1)*FAC(IA8+1)*FAC(IA9+1)
      A=AN/AD
C
C     CALCULATION OF SUM 'B'
C
      IB1=NL1+NM1
      IB2=NL2+NL3-NM1
      IB3=NL1-NM1
      IB4=NL2-NL3+NM1
      IB5=NL3-NM3
      IT1=MAX0(0,-IB4)+1
      IT2=MIN0(IB2,IB3,IB5)+1
      B=0.0_dp
      SIGN=(-1.0_dp)**(IT1)
      IB1=IB1+IT1-2
      IB2=IB2-IT1+2
      IB3=IB3-IT1+2
      IB4=IB4+IT1-2
      IB5=IB5-IT1+2
      DO 18 IT=IT1,IT2
      SIGN=-SIGN
      IB1=IB1+1
      IB2=IB2-1
      IB3=IB3-1
      IB4=IB4+1
      IB5=IB5-1
      BN=SIGN*FAC(IB1+1)*FAC(IB2+1)
      BD=FAC(IT)*FAC(IB3+1)*FAC(IB4+1)*FAC(IB5+1)
  18  B=B+(BN/BD)
C
C       CALCULATION OF FACTOR 'C'
C
      IC1=NL1-NM1
      IC2=NL1+NM1
      IC3=NL2-NM2
      IC4=NL2+NM2
      IC5=NL3-NM3
      IC6=NL3+NM3
      CN=DFLOAT((2*NL1+1)*(2*NL2+1)*(2*NL3+1))*FAC(IC1+1)*FAC(IC3+1)*
     1FAC(IC5+1)
      CD=FAC(IC2+1)*FAC(IC4+1)*FAC(IC6+1)
      C=CN/(PI*CD)
      C=(SQRT(C))/2.0_dp
      BLM=((-1.0_dp)**IC)*A*B*C
      RETURN
  19  WRITE(6,20)L1,L2,M2,L3,M3
  20  FORMAT(28H INVALID ARGUMENTS FOR BLM. ,5(I2,1H,))
      RETURN
      END function
C=======================================================================
C=======================================================================
      SUBROUTINE PAIR(IGKMAX,QIL,QIIL,QIIIL,QIVL,QIR,QIIR,QIIIR,QIVR)

C     ------------------------------------------------------------------
C     THIS SUBROUTINE CALCULATES SCATTERING Q-MATRICES FOR A  DOUBLE
C     LAYER, FROM THE CORRESPONDING MATRICES OF THE INDIVIDUAL, LEFT
C     (L) AND RIGHT (R), LAYERS. THE RESULTS ARE STORED IN Q*L.
C     -----------------------------------------------------------------
C
C ..  PARAMETER STATEMENTS  ..
C
      INTEGER   IGD,IGKD
      PARAMETER (IGD=21,IGKD=2*IGD)
C
C ..  SCALAR ARGUMENTS  ..
C
      INTEGER IGKMAX
C
C ..  ARRAY ARGUMENTS  ..
C
      COMPLEX(dp) QIL (IGKD,IGKD),QIIL(IGKD,IGKD),QIIIL(IGKD,IGKD)
      COMPLEX(dp) QIVL(IGKD,IGKD)
      COMPLEX(dp) QIR (IGKD,IGKD),QIIR(IGKD,IGKD),QIIIR(IGKD,IGKD)
      COMPLEX(dp) QIVR(IGKD,IGKD)
C
C ..  LOCAL
      INTEGER    IGK1,IGK2,IGK3
      REAL(dp)   EMACH
      INTEGER    INT(IGKD),JNT(IGKD)
      COMPLEX(dp) QINV1(IGKD,IGKD),QINV2(IGKD,IGKD),W1(IGKD,IGKD)
      COMPLEX(dp) W2(IGKD,IGKD),W3(IGKD,IGKD),W4(IGKD,IGKD)
C
      DATA EMACH/1.D-8/
C-----------------------------------------------------------------------
C
      DO IGK1=1,IGKMAX
        DO IGK2=1,IGKMAX
          QINV1(IGK1,IGK2)=QIL (IGK1,IGK2)
          QINV2(IGK1,IGK2)=QIVR(IGK1,IGK2)
          W2(IGK1,IGK2)=CZERO
          W3(IGK1,IGK2)=CZERO
        end do
      end do

      DO IGK1=1,IGKMAX
        W2(IGK1,IGK1)=CONE
        W3(IGK1,IGK1)=CONE
        DO IGK2=1,IGKMAX
          DO IGK3=1,IGKMAX
            W2(IGK1,IGK2)=W2(IGK1,IGK2)
     &         -QIIL (IGK1,IGK3)*QIIIR(IGK3,IGK2)
            W3(IGK1,IGK2)=W3(IGK1,IGK2)
     &         -QIIIR(IGK1,IGK3)*QIIL (IGK3,IGK2)
          end do
        end do
      end do

!     call zgetrf_wrap(w2, int)
!     call zgetrf_wrap(w3, jnt)
!     call zgetrs_wrap(w2, QINV1(1,:), int)
!     call zgetrs_wrap(w3, QINV2(1,:), jnt)

       CALL ZGE(W2,INT,IGKMAX,IGKD,EMACH)
       CALL ZGE(W3,JNT,IGKMAX,IGKD,EMACH)
       DO IGK2=1,IGKMAX
         CALL ZSU(W2,INT,QINV1(1,IGK2),IGKMAX,IGKD,EMACH)
         CALL ZSU(W3,JNT,QINV2(1,IGK2),IGKMAX,IGKD,EMACH)
       end do

      DO IGK1=1,IGKMAX
        DO IGK2=1,IGKMAX
          W1(IGK1,IGK2)=CZERO
          W2(IGK1,IGK2)=CZERO
          W3(IGK1,IGK2)=CZERO
          W4(IGK1,IGK2)=CZERO
          DO IGK3=1,IGKMAX
            W1(IGK1,IGK2)=W1(IGK1,IGK2)
     &                    +QIR  (IGK1,IGK3)*QINV1(IGK3,IGK2)
            W2(IGK1,IGK2)=W2(IGK1,IGK2)
     &                    +QIIL (IGK1,IGK3)*QINV2(IGK3,IGK2)
            W3(IGK1,IGK2)=W3(IGK1,IGK2)
     &                    +QIIIR(IGK1,IGK3)*QINV1(IGK3,IGK2)
            W4(IGK1,IGK2)=W4(IGK1,IGK2)
     &                    +QIVL (IGK1,IGK3)*QINV2(IGK3,IGK2)
          end do
        end do
      end do

      DO IGK1=1,IGKMAX
        DO IGK2=1,IGKMAX
          QINV1(IGK1,IGK2)=QIIR (IGK1,IGK2)
          QINV2(IGK1,IGK2)=QIIIL(IGK1,IGK2)
          DO IGK3=1,IGKMAX
            QINV1(IGK1,IGK2)=QINV1(IGK1,IGK2)
     &                       +QIR (IGK1,IGK3)*W2(IGK3,IGK2)
            QINV2(IGK1,IGK2)=QINV2(IGK1,IGK2)
     &                       +QIVL(IGK1,IGK3)*W3(IGK3,IGK2)
          end do
        end do
      end do

      DO IGK1=1,IGKMAX
        DO IGK2=1,IGKMAX
          QIL  (IGK1,IGK2)=W1   (IGK1,IGK2)
          QIIL (IGK1,IGK2)=QINV1(IGK1,IGK2)
          QIIIL(IGK1,IGK2)=QINV2(IGK1,IGK2)
          QIVL (IGK1,IGK2)=W4   (IGK1,IGK2)
        end do
      end do

      RETURN
      END subroutine
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
      INTEGER   LMAXD,LMAX1D,LM1SQD,IGD,IGKD,NELMD
      PARAMETER (LMAXD=14,LMAX1D=LMAXD+1)
      PARAMETER (LM1SQD=LMAX1D*LMAX1D,IGD=21,IGKD=2*IGD,NELMD=165152)
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
      complex(dp) QI(IGKD,IGKD),QII(IGKD,IGKD),QIII(IGKD,IGKD)
      complex(dp) QIV(IGKD,IGKD)
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
      complex(dp) AE(2,LM1SQD),AH(2,LM1SQD),GKK(3,IGD)
      complex(dp) GK(3),LAME(2),LAMH(2)
      complex(dp) XEVEN(((lmax+1)*(lmax+2))/2,((lmax+1)*(lmax+2))/2),
     & XODD((lmax*(lmax+1))/2,(lmax*(lmax+1))/2)
      complex(dp) :: TE(lmax+1),TH(lmax+1)
      complex(dp) BMEL1((LMAX+1)**2-1),BMEL2((LMAX+1)**2-1)
      complex(dp) XXMAT1((LMAX+1)**2-1,(LMAX+1)**2-1),
     & XXMAT2((LMAX+1)**2-1,(LMAX+1)**2-1)
      complex(dp) DLME(2,LM1SQD),DLMH(2,LM1SQD)
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
      SUBROUTINE BAND(IGMAX,ZVAL,EMACH,AK,G,AL,KAPL,KAPR,
     &                QI,QII,QIII,QIV)

C     ------------------------------------------------------------------
C     THIS SUBROUTINE CALCULATES THE COMPLEX PHOTONIC BAND STRUCTURE OF
C     AN INFINITE CRYSTAL. IT  PROVIDES THE  PROPAGATING AND EVANESCENT
C     EIGENMODES OF THE EM FIELD IN THE GIVEN CRYSTAL, CORRESPONDING TO
C     "AK" AND A GIVEN FREQUENCY.
C     ------------------------------------------------------------------
C
C ..  PARAMETER STATEMENTS ..
C
      INTEGER   IGD,IGKD,IGK2D
      PARAMETER(IGD=21,IGKD=2*IGD,IGK2D=2*IGKD)
C
C ..  SCALAR ARGUMENTS  ..
C
      INTEGER    IGMAX
      REAL(dp)   ZVAL,EMACH
      COMPLEX(dp) KAPL,KAPR
C
C ..  ARRAY ARGUMENTS ..
C
      REAL(dp)   AK(2),G(2,IGD),AL(3)
      COMPLEX(dp) QI(IGKD,IGKD),QII(IGKD,IGKD)
      COMPLEX(dp) QIII(IGKD,IGKD),QIV(IGKD,IGKD)
C
C ..  SCALAR VARIABLES ..
C
      INTEGER    II,I,IGK1,IGK2,IGKMAX,J
      INTEGER    KD,LIB1,LIB2,LU,LP,LN,IFAIL,IGK3,IGK2M
      REAL(dp)   AKA,BKZRE,BKZIM
      COMPLEX(dp) EAKA
C
C ..  ARRAY VARIABLES ..
C
      INTEGER    INT(IGKD)
      REAL(dp)   AR(IGK2D,IGK2D),AI(IGK2D,IGK2D)
      REAL(dp)   RR(IGK2D),RI(IGK2D),VR(IGK2D,IGK2D),VI(IGK2D,IGK2D)
      REAL(dp)   AKZAP(IGK2D),AKZIP(IGK2D)
      REAL(dp)   AKZREP(IGK2D),AKZIMP(IGK2D),AKZREN(IGK2D),AKZIMN(IGK2D)
      COMPLEX(dp) QH1(IGKD,IGKD),QH2(IGKD,IGKD),AKZ(IGK2D)
      COMPLEX(dp) COMVEC(IGK2D,IGK2D)
C     ------------------------------------------------------------------
      IGKMAX=2*IGMAX
      IGK2M=2*IGKMAX
      AKA=AK(1)*AL(1)+AK(2)*AL(2)
      EAKA=EXP(CI*AKA)
      DO 1 IGK1=1,IGKMAX
      DO 1 IGK2=1,IGKMAX
      QH1(IGK1,IGK2)=CZERO
      QH2(IGK1,IGK2)=CZERO
    1 CONTINUE
      DO 5 IGK1=1,IGKMAX
      QH2(IGK1,IGK1)=CONE
      DO 5 IGK2=1,IGKMAX
      DO 6 IGK3=1,IGKMAX
      QH1(IGK1,IGK2)=QH1(IGK1,IGK2)-QIII(IGK1,IGK3)*QI (IGK3,IGK2)
    6 QH2(IGK1,IGK2)=QH2(IGK1,IGK2)-QIII(IGK1,IGK3)*QII(IGK3,IGK2)
    5 CONTINUE

!     call zgetrf_wrap(QIV, int)
!     call zgetrs_wrap(QIV, QH1(1,:), int)
!     call zgetrs_wrap(QIV, QH2(1,:), int)

      CALL ZGE(QIV,INT,IGKMAX,IGKD,EMACH)
      do IGK2=1,IGKMAX
        CALL ZSU(QIV,INT,QH1(1,IGK2),IGKMAX,IGKD,EMACH)
        CALL ZSU(QIV,INT,QH2(1,IGK2),IGKMAX,IGKD,EMACH)
      end do

      DO 9 IGK1=1,IGKMAX
      DO 9 IGK2=1,IGKMAX
      AR(IGK1,IGK2)=DREAL(QI(IGK1,IGK2))
      AI(IGK1,IGK2)=DIMAG(QI(IGK1,IGK2))
      AR(IGK1,IGKMAX+IGK2)=DREAL(QII(IGK1,IGK2))
      AI(IGK1,IGKMAX+IGK2)=DIMAG(QII(IGK1,IGK2))
      AR(IGKMAX+IGK1,IGK2)=DREAL(QH1(IGK1,IGK2))
      AI(IGKMAX+IGK1,IGK2)=DIMAG(QH1(IGK1,IGK2))
      AR(IGKMAX+IGK1,IGKMAX+IGK2)=DREAL(QH2(IGK1,IGK2))
      AI(IGKMAX+IGK1,IGKMAX+IGK2)=DIMAG(QH2(IGK1,IGK2))
    9 CONTINUE
      CALL CNAA(IGK2D,IGK2M,AR,AI,RR,RI,VR,VI,IFAIL)
      IF(IFAIL/=0) THEN
      WRITE(6,102) IFAIL
                         STOP
                     ENDIF
      DO 8 II=1,IGK2M
C*****THE IF-STRUCTURE  WHICH FOLLOWS  CAN BE  OMITTED  IF THE ACCURACY
C*****'MACHEP' OF THE SUBROUTINE COMLR2 IS CHOSEN GREATER THAN 2**(-47)
      IF((RR(II)==0.0_dp).AND.(RI(II)==0.0_dp)) THEN
      RR(II)=1.D-20
      RI(II)=1.D-20
      ENDIF
      AKZ(II)=(-CI/PI)*LOG(DCMPLX(RR(II),RI(II))/EAKA) ! NORMALIZED K_Z
    8 CONTINUE
      DO 2 LIB2=1,IGK2M
      DO 2 LIB1=1,IGK2M
      COMVEC(LIB1,LIB2)=VR(LIB1,LIB2)+CI*VI(LIB1,LIB2)
    2 CONTINUE
      LU=1
      LP=1
      LN=1
      DO 10 KD=1,IGK2M
C*****WARNING!! THE APPROPRIATE LIMITS FOR DIMAG(AKZ(KD))
C*****DEPEND STRONGLY ON IGMAX.
      IF(DIMAG(AKZ(KD))>0.0_dp) THEN
      AKZREP(LP)=DREAL(AKZ(KD))
      AKZIMP(LP)=DIMAG(AKZ(KD))
      LP=LP+1
      ELSE
      AKZREN(LN)=DREAL(AKZ(KD))
      AKZIMN(LN)=DIMAG(AKZ(KD))
      LN=LN+1
      ENDIF
      IF(ABS(DIMAG(AKZ(KD)))>1.0D-2) GO TO 10
      AKZAP(LU)=DREAL(AKZ(KD))
      AKZIP(LU)=DIMAG(AKZ(KD))
      LU=LU+1
   10 CONTINUE
      IF (LU<1.1D0) THEN
      DO 13 J=2,LP-1
      BKZIM=AKZIMP(J)
      BKZRE=AKZREP(J)
      DO 14 I=J-1,1,-1
      IF(AKZIMP(I)<=BKZIM) GO TO 15
      AKZIMP(I+1)=AKZIMP(I)
      AKZREP(I+1)=AKZREP(I)
   14 CONTINUE
      I=0
   15 AKZIMP(I+1)=BKZIM
      AKZREP(I+1)=BKZRE
   13 CONTINUE
      DO 16 J=2,LN-1
      BKZIM=AKZIMN(J)
      BKZRE=AKZREN(J)
      DO 17 I=J-1,1,-1
      IF(AKZIMN(I)<=BKZIM) GO TO 18
      AKZIMN(I+1)=AKZIMN(I)
      AKZREN(I+1)=AKZREN(I)
   17 CONTINUE
      I=0
   18 AKZIMN(I+1)=BKZIM
      AKZREN(I+1)=BKZRE
   16 CONTINUE
      WRITE(6,101)  ZVAL,AKZREP(1),AKZREN(LN-1)
      WRITE(6,103)  AKZIMP(1),AKZIMN(LN-1)
      WRITE(9,101)  ZVAL,AKZREP(1),AKZREN(LN-1)
      WRITE(9,103)  AKZIMP(1),AKZIMN(LN-1)
      ELSE
      WRITE(6,101)  ZVAL,(AKZAP(I),I=1,LU-1)
      WRITE(9,101)  ZVAL,(AKZAP(I),I=1,LU-1)
      END IF
      RETURN
  101 FORMAT(E10.4,3X,10(E10.4,1X))
  102 FORMAT(//13X,'ERROR IN CNAA   IFAIL = ',I2)
  103 FORMAT(13X,10(E10.4,1X))
      END subroutine
C======================================================================
      SUBROUTINE REDUCE(AR1,AR2,AK,IGMAX,G,IG0,EMACH)

C----------------------------------------------------------------------
C     GIVEN THE PRIMITIVE VECTORS AR1,AR2 OF A 2D LATTICE (IN UNITS OF
C     ALPHA), THIS  SUBROUTINE  REDUCES A WAVECTOR "AK" (IN  UNITS  OF
C     2*PI/ALPHA) WITHIN THE SBZ BY ADDING AN APPROPRIATE  RECIPROCAL-
C     LATTICE VECTOR G(IG0)
C----------------------------------------------------------------------
C
C ..  PARAMETER STATEMENTS  ..
C
      INTEGER IGD
      PARAMETER (IGD=21)
C
C ..  SCALAR ARGUMENTS  ..
C
      INTEGER IGMAX,IG0
      REAL(dp) EMACH
C
C ..  ARRAY ARGUMENTS  ..
C
      REAL(dp) AR1(2),AR2(2),AK(2),G(2,IGD)
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
      integer, parameter:: dp=kind(0.d0)
      real(dp), parameter :: pi=4.0_dp*ATAN(1.0_dp)
      complex(dp), parameter :: czero = (0.0_dp, 0.0_dp)
!     complex(dp), parameter :: ci    = (0.0_dp, 1.0_dp)
      complex(dp), parameter :: cone  = (1.0_dp, 0.0_dp)
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
      CHARACTER*2  POLAR
      CHARACTER*17 TEXT1(2)
      CHARACTER*5  DUMMY
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
C ..  INTRINSIC FUNCTIONS ..
C
      INTRINSIC DCMPLX,SQRT,DREAL,SIN,COS,DFLOAT,ABS,DIMAG
C
C ..  EXTERNAL ROUTINES ..
C
C     EXTERNAL ELMGEN,LAT2D,PCSLAB,HOSLAB,PAIR,SCAT,BAND,REDUCE
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
      IF(DREAL(MU1(ICOMP))<=0.D0.OR.DREAL(EPS1(ICOMP))<=0.D0)
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
      IF(DIMAG(D1)/=0.D0) THEN
      WRITE(6,227)
      STOP
      ENDIF
      IF(DIMAG(D2)/=0.D0) THEN
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
      IF(DIMAG(D1)/=0.D0) THEN
      WRITE(6,227)
      STOP
      ENDIF
      IF(DIMAG(D2)/=0.D0) THEN
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
      ZSTEP=(ZSUP-ZINF)/DFLOAT(NP-1)
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
      IF(KSCAN==1) KAPPA0=DCMPLX(ZVAL,EPSILON)
      IF(KSCAN==2) KAPPA0=DCMPLX(2.D0*PI/ZVAL,EPSILON)
      KAPIN =KAPPA0*D1
      KAPOUT=KAPPA0*D2
                                             IF(KTYPE==1) THEN
                           AK(1)=DREAL(KAPIN)*SIN(THETA)*COS(FI)
                           AK(2)=DREAL(KAPIN)*SIN(THETA)*SIN(FI)
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
      IF(DREAL(AKZIN)<EMACH)      STOP 'IMPROPER INCIDENT WAVE'
      AKXY=SQRT(AKXY)
      IF(AKXY<EMACH) THEN
      EIN(1)=DCMPLX(COS(FEIN),0.D0)
      EIN(2)=DCMPLX(SIN(FEIN),0.D0)
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
      CALL PAIR(IGKMAX,QIL,QIIL,QIIIL,QIVL,QIR,QIIR,QIIIR,QIVR)
   13 CONTINUE
            ENDIF
      IF(NLAYER(1)>=2) THEN
      DO 14 ILAYER=1,NLAYER(1)-1
      CALL PAIR(IGKMAX,QIL,QIIL,QIIIL,QIVL,QIL,QIIL,QIIIL,QIVL)
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
      CALL PAIR(IGKMAX,WIL,WIIL,WIIIL,WIVL,QIR,QIIR,QIIIR,QIVR)
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
      CALL PAIR(IGKMAX,QIR,QIIR,QIIIR,QIVR,QIR,QIIR,QIIIR,QIVR)
   16 CONTINUE
             ENDIF
                      ENDIF
      CALL PAIR(IGKMAX,QIL,QIIL,QIIIL,QIVL,QIR,QIIR,QIIIR,QIVR)
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
             CALL PAIR(IGKMAX,QIL,QIIL,QIIIL,QIVL,QIL,QIIL,QIIIL,QIVL)
    9        CONTINUE
   30        CONTINUE
             IF(KEMB==1) THEN
             CALL HOSLAB(IGMAX,KAPR,(KAPR+KAPOUT)/2.D0,KAPOUT,AK,G,VEC0,
     &                   VEC0,0.D0,QIR,QIIR,QIIIR,QIVR,EMACH)
         CALL PAIR(IGKMAX,QIL,QIIL,QIIIL,QIVL,QIR,QIIR,QIIIR,QIVR)
         DO 11 IGK1=1,IGKMAX
         DO 11 IGK2=1,IGKMAX
         QIR  (IGK1,IGK2)=QIL  (IGK1,IGK2)
         QIIR (IGK1,IGK2)=QIIL (IGK1,IGK2)
         QIIIR(IGK1,IGK2)=QIIIL(IGK1,IGK2)
         QIVR (IGK1,IGK2)=QIVL (IGK1,IGK2)
   11        CONTINUE
             CALL HOSLAB(IGMAX,KAPIN,(KAPL+KAPIN)/2.D0,KAPL,AK,G,VEC0,
     &                   VEC0,0.D0,QIL,QIIL,QIIIL,QIVL,EMACH)
         CALL PAIR(IGKMAX,QIL,QIIL,QIIIL,QIVL,QIR,QIIR,QIIIR,QIVR)
                      ENDIF
             CALL SCAT(IGMAX,ZVAL,AK,G,DREAL(KAPIN),DREAL(KAPOUT),
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
