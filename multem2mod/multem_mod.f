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

      PROGRAM MULTEM
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
C      COMMON/X1/AR1,AR2
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
      FAB = FAB*pi/180.0_dp
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
      DO ICOMP=1,NCOMP
        READ(10,201) IT(ICOMP)
        IF(IT(ICOMP)<=0.OR.IT(ICOMP)>2)
     &                      STOP 'ILLEGAL COMPONENT TYPE'
        WRITE(6,209) ICOMP,TEXT1(IT(ICOMP))
        IF(IT(ICOMP)==1) THEN
        READ(10,204) D(ICOMP)
        READ(10,205) MU1(ICOMP),EPS1(ICOMP),MU2(ICOMP),EPS2(ICOMP),
     &               MU3(ICOMP),EPS3(ICOMP)
        WRITE(6,210) MU1(ICOMP),MU2(ICOMP),MU3(ICOMP),EPS1(ICOMP),
     &               EPS2(ICOMP),EPS3(ICOMP)
        READ(10,*) DUMMY,(DL(I,ICOMP,1),I=1,3)
        READ(10,*) DUMMY,(DR(I,ICOMP,1),I=1,3)
                        ELSE
        READ(10,205) MU1(ICOMP),EPS1(ICOMP)
        IF(dble(MU1(ICOMP))<=0.D0.OR.dble(EPS1(ICOMP))<=0.D0)
     &  THEN
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
      end do
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
               DO I=1,IGKMAX
               EINCID(I)=CZERO
               end do
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
      DO I=1,2
        EINCID(2*IG0-2+I)=EIN(I)
      end do
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
     &         QIL,QIIL,QIIIL,QIVL, ar1, ar2)
      IF(NPLAN(1)>=2) THEN
      DO IPL=2,NPLAN(1)
        RAP=S(1,IPL)*KAPPA0/2.D0/PI
        CALL PCSLAB(LMAX,IGMAX,RAP,EPS1(1),EPSSPH(1,IPL),MU1(1),
     &             MUSPH(1,IPL),KAPPA,AK,DL(1,1,IPL),DR(1,1,IPL),
     &             G,ELM,A0,EMACH, QIR,QIIR,QIIIR,QIVR, ar1, ar2)
        CALL PAIR(IGKMAX,igkd,QIL,QIIL,QIIIL,QIVL,QIR,QIIR,QIIIR,QIVR)
      end do
            ENDIF
      IF(NLAYER(1)>=2) THEN
      DO ILAYER=1,NLAYER(1)-1
        CALL PAIR(IGKMAX,igkd,QIL,QIIL,QIIIL,QIVL,QIL,QIIL,QIIIL,QIVL)
      end do
             ENDIF
                      ENDIF
                                                    IF(NCOMP>=2) THEN
      DO ICOMP=2,NCOMP
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
     &           DR(1,ICOMP,1),G,ELM,A0,EMACH,QIR,QIIR,QIIIR,QIVR,
     &           ar1, ar2)
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
     &             QIR,QIIR,QIIIR,QIVR, ar1, ar2)
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
      end do
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
