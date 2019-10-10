      SUBROUTINE BIGA (CIOR,XX,NTRM,NOABS,YESANG,CBIGA)
ct      program biga
C--------/---------/---------/---------/---------/---------/---------/--
c    Calculate logarithmic derivatives of J-Bessel-function
c    Input :  CIOR, XX, NTRM, NOABS, YESANG  
c    Output :  RBIGA or CBIGA  
c    Routines called :  CONFRA
c              ERRMSG('ConFra--Iteration failed to converge',.TRUE.)
c  =====
c
c  COMPLEX*16   CBIGA( LMAXD )
c  call BIGA (CIOR, XX, LMAXD, .false., .false., RBIGA, CBIGA )
c
c According to:
c  [1] Lentz ???
c  [2] W. J. Wiscombe, Appl. Opt. 19, 1505-1509 (1980).
c                      Improved Mie scattering algorithms.
c     (http://www.opticsinfobase.org/ao/abstract.cfm?URI=ao-19-9-1505)
c  [R] W. J. Wiscombe, NCAR Technical Note
c
c =====
c  CIOR            Complex index of refraction. Should be
c                  now in scattering notation (Im part of n always nonnegative)
c  XX               size parameter
c  NTRM            No. of terms in Mie series (my LMAXV)
c  NOABS           TRUE, sphere non-absorbing (determined by MIMCUT=0.d0)
c  YESANG          TRUE if scattering amplitudes are to be calculated
c
c  RBIGA(N)        Bessel function ratio A_N=\psi_N'/\psi_N (Ref. 2, Eq. 2)
c                     (REAL version, for when imag refrac index = 0)
c  CBIGA(N)        Bessel function ratio A_N=\psi_N'/\psi_N (Ref. 2, Eq. 2)
c                     ( COMPLEX version )
c
c    INTERNAL VARIABLES :
c
c       CONFRA     Value of Lentz continued fraction for CBIGA(NTRM),
c                     used to initialize downward recurrence
c
c       DOWN       = True, use down-recurrence.  False, do not.
c
c       F1,F2,F3   Arithmetic statement functions used in determining
c                     whether to use up-  or down-recurrence
c                     ( Ref. 2, Eqs. 6-8 )
c
c       MRE        Real refractive index
c       MIM        Imaginary refractive index
c
c       REZINV     1 / ( MRE * XX ); temporary variable for recurrence
c       ZINV       1 / ( CIOR * XX ); temporary variable for recurrence
C--------/---------/---------/---------/---------/---------/---------/--

      IMPLICIT  NONE

c     .. Scalar Arguments ..

      INTEGER    NTRM,LMAXD
      PARAMETER (LMAXD=60)

      LOGICAL    NOABS, YESANG
      REAL*8     XX
      COMPLEX*16 CIOR
c     ..
c     .. Array Arguments ..

      REAL*8      RBIGA(LMAXD)
      COMPLEX*16   CBIGA(LMAXD)
c     ..
c     .. Local Scalars ..

      LOGICAL   DOWN
      INTEGER   N
      REAL*8      MIM, MRE, REZINV, RTMP
      COMPLEX*16   CTMP,ZINV
c     ..
c     .. External Functions ..

      COMPLEX*16   CONFRA
      EXTERNAL  CONFRA
c     ..
c     .. Intrinsic Functions ..
c
c      INTRINSIC ABS, AIMAG, COS, EXP, REAL, SIN
c     ..
c     .. Statement Functions ..
c
      REAL*8      F1,F2
c     ,F3  - not used here
c     ..
c     .. Statement Function definitions ..
c
c                                                   ** Eq. R47c
      F1(MRE) = -8.d0 + MRE**2*( 26.22d0 +
     &  MRE*( -0.4474d0 + MRE**3*( 0.00204d0 - 0.000175d0*MRE ) ) )

c                                                   ** Eq. R47b
      F2(MRE) = 3.9d0 + MRE*( -10.8d0 + 13.78d0*MRE )
c                                                   ** Eq. R47a
cc      F3(MRE) = -15.04d0 + MRE*( 8.42d0 + 16.35d0*MRE )
c     ..

ct      NOABS=.false. 
ct      YESANG=.false. 
ct      XX=7.31806289774349d0
ct      CIOR=(1.45d0,1.d0)

c                                  ** Decide whether BigA can be
c                                  ** calculated by up-recurrence
      MRE = DBLE(CIOR)
      MIM = ABS(IMAG(CIOR))

      IF( MRE.LT.1.0 .OR. MRE.GT.10.0 .OR. MIM.GT.10.0 ) THEN

         DOWN = .TRUE.

      ELSE IF (YESANG) THEN

         DOWN = .TRUE.
c                                                    ** Eq. R48
         IF( MIM*XX .LT. F2(MRE) ) DOWN = .FALSE.

      ELSE

         DOWN = .TRUE.
c                                                    ** Eq. R48
         IF( MIM*XX .LT. F1(MRE) ) DOWN = .FALSE.

      END IF


      ZINV   = 1.d0 /(CIOR*XX)
      REZINV = 1.d0 /(MRE*XX)


      IF(DOWN) THEN
c                          ** Compute initial high-order BigA using
c                          ** Lentz method ( Ref. 1, pp. 17-20 )

         CTMP = CONFRA(NTRM,ZINV)

c                                   *** Downward recurrence for BigA
         IF(NOABS) THEN
c                                        ** No-absorption case; Eq (R23)
            RBIGA(NTRM) = DBLE(CTMP)
            CBIGA(NTRM) = RBIGA(NTRM)

            DO 10 N = NTRM, 2, -1

               RBIGA(N - 1) = ( N*REZINV ) -
     &                        1.d0 / ( ( N*REZINV ) + RBIGA(N) )

               CBIGA(N-1) = RBIGA(N-1)

   10       CONTINUE

         ELSE
c                                         ** Absorptive case; Eq (R23)
            CBIGA(NTRM) = CTMP

            DO 20 N = NTRM, 2, -1
            CBIGA(N-1) = (N*ZINV) - 1.d0 / ( (N*ZINV) + CBIGA(N) )
   20       CONTINUE

         END IF


      ELSE
c                            *** Upward recurrence for BigA
         IF(NOABS) THEN
c                                  ** No-absorption case; Eq (R20,21)
            RTMP = SIN( MRE*XX )
            RBIGA(1) = - REZINV + RTMP /
     &                   ( RTMP*REZINV - COS(MRE*XX) )
            CBIGA(1) = RBIGA(1)

            DO 30 N = 2, NTRM
               RBIGA(N) = -( N*REZINV ) +
     &                      1.d0 / ( (N*REZINV) - RBIGA(N - 1) )
             CBIGA(N) = RBIGA(N)
   30       CONTINUE

         ELSE
c                                     ** Absorptive case; Eq (R20,22)

            CTMP = EXP((0.d0,2.d0)*CIOR*XX)
            CBIGA(1) = - ZINV + (1.d0-CTMP) /
     &              (ZINV * (1.d0-CTMP) + (0.d0,1.d0)*(1.d0+CTMP))

            DO 40 N = 2, NTRM
            CBIGA(N) = - (N*ZINV) + 1.d0/((N*ZINV) - CBIGA(N-1))
   40       CONTINUE

         END IF

      END IF

      RETURN
      END

      COMPLEX*16 FUNCTION CONFRA(N,ZINV)
C--------/---------/---------/---------/---------/---------/---------/--
c         Compute Bessel function ratio A-sub-N from its
c         continued fraction using Lentz method
c
c         ZINV = Reciprocal of argument of A
c
c    I N T E R N A L    V A R I A B L E S
c    ------------------------------------
c    CAK      Term in continued fraction expansion of A (Eq. R25)
c    CAPT     Factor used in Lentz iteration for A (Eq. R27)
c    CNUMER   Numerator   in capT  ( Eq. R28A )
c    CDENOM   Denominator in capT  ( Eq. R28B )
c    CDTD     Product of two successive denominators of capT factors
c                 ( Eq. R34C )
c    CNTN     Product of two successive numerators of capT factors
c                 ( Eq. R34B )
c    EPS1     Ill-conditioning criterion
c    EPS2     Convergence criterion
c    KK       Subscript k of cAk  ( Eq. R25B )
c    KOUNT    Iteration counter ( used to prevent infinite looping )
c    MAXIT    Max. allowed no. of iterations
c    MM       + 1  and - 1, alternately
C--------/---------/---------/---------/---------/---------/---------/--
      IMPLICIT  NONE

c     .. Scalar Arguments ..

      INTEGER   N
      COMPLEX*16   ZINV
c     ..
c     .. Local Scalars ..

      INTEGER   KK, KOUNT, MAXIT, MM
      REAL*8      EPS1, EPS2
      COMPLEX*16   CAK, CAPT, CDENOM, CDTD, CNTN, CNUMER
c     ..
c     .. External Subroutines ..

cc      EXTERNAL  ERRMSG
c     ..
c     .. Intrinsic Functions ..
c
c      INTRINSIC ABS, AIMAG, REAL
c     ..
      DATA      EPS1 / 1.d-2 / , EPS2 / 1.d-8 /
      DATA      MAXIT / 10000 /

c                                 ** Eq. R25a
      CONFRA = ( N + 1 ) * ZINV
      MM     = - 1
      KK     = 2*N + 3
c                                 ** Eq. R25b, k=2
      CAK    = ( MM*KK ) * ZINV
      CDENOM = CAK
      CNUMER = CDENOM + 1.d0 / CONFRA
      KOUNT  = 1

   10 CONTINUE
      KOUNT = KOUNT + 1

      IF (KOUNT.GT.MAXIT) then
cx     &    CALL ERRMSG('ConFra--Iteration failed to converge',.TRUE.)
          write(6,*)'ConFra--Iteration failed to converge'
c    pause
      end if

      MM  = - MM
      KK  = KK + 2
c                                 ** Eq. R25b
      CAK = ( MM*KK ) * ZINV
c                                          ** Eq. R32
      IF( ABS( CNUMER / CAK ).LE.EPS1 .OR.
     &    ABS( CDENOM / CAK ).LE.EPS1 ) THEN

c                                  ** Ill-conditioned case -- stride
c                                  ** two terms instead of one
c                                       ** Eq. R34
         CNTN   = CAK * CNUMER + 1.d0
         CDTD   = CAK * CDENOM + 1.d0
c                                           ** Eq. R33
         CONFRA = ( CNTN / CDTD ) * CONFRA

         MM  = - MM
         KK  = KK + 2
c                                 ** Eq. R25b
         CAK = ( MM*KK ) * ZINV
c                                      ** Eq. R35
         CNUMER = CAK + CNUMER / CNTN
         CDENOM = CAK + CDENOM / CDTD
         KOUNT  = KOUNT + 1
         GO TO  10

      ELSE
c                           *** Well-conditioned case
c                                  ** Eq. R27
         CAPT   = CNUMER / CDENOM
c                                  ** Eq. R26
         CONFRA = CAPT * CONFRA
c                                  ** Check for convergence; Eq. R31

         IF (      ABS( DBLE (CAPT) - 1.d0 ).GE.EPS2
     &        .OR. ABS( IMAG(CAPT) )      .GE.EPS2 )  THEN

c                                        ** Eq. R30
            CNUMER = CAK + 1.d0 / CNUMER
            CDENOM = CAK + 1.d0 / CDENOM

            GO TO  10

         END IF

      END IF

      RETURN
      END
