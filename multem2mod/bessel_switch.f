      SUBROUTINE BESSEL(BJ,Y,H,ARG,LMX,lmax,LJ,LY,LH,LCALL)
      IMPLICIT NONE
!     ------------------------------------------------------------------
!     THIS  SUBROUTINE COMPUTES THE  SPHERICAL BESSEL FUNCTIONS OF
!     FIRST, SECOND  AND  THIRD  KIND  using Amos lib
!
!     2019.07.17 Change to use Amos lib
!
!     ON INPUT--->
!     ARG    ARGUMENT OF THE BESSEL FUNCTIONS
!     LMAX   MAX. ORDER OF THE BESSEL FUNCTIONS
!            (LIMITED UP TO 25 IN THE VERSION)
!     LJ     LOGICAL : IF LJ IS TRUE THE SPHERICAL BESSEL
!            FUNCTIONS OF THE FIRST KIND ARE CALCULATED UP TO LMAX
!     LY     LOGICAL : IF LY IS TRUE THE SPHERICAL BESSEL
!            FUNCTIONS OF THE SECOND KIND ARE CALCULATED UP TO LMAX
!     LH     LOGICAL : IF LH IS TRUE THE SPHERICAL BESSEL
!            FUNCTIONS OF THE THIRD KIND ARE CALCULATED UP TO LMAX
!     LCALL  LOGICAL : IF LCALL IS FALSE THE CHEBYCHEV
!            COEFFICIENTS ARE CALCULATED -THIS PART HAS TO
!            BE CALLED ONCE
!
!     ON OUTPUT--->
!     BJ     AN ARRAY CONTAINING THE BESSEL FUNCTIONS OF
!            THE FIRST KIND UP TO LMAX IF LJ IS TRUE.
!            REMEMBER, THAT BJ(1) CONTAINS THE FUNCTION OF
!            L=0 AND SO ON.
!     Y      AN ARRAY CONTAINING THE BESSEL FUNCTIONS OF
!            THE SECOND KIND UP TO LMAX IF LY IS TRUE.
!            REMEMBER,THAT  Y(1) CONTAINS THE FUNCTION OF L=0 AND SO ON.
!     H      AN ARRAY CONTAINING THE BESSEL FUNCTIONS OF
!            THE THIRD KIND UP TO LMAX IF LH IS TRUE.
!            REMEMBER,THAT H (1) CONTAINS THE FUNCTION OF L=0 AND SO ON.
!
!     THE BESSEL FUNCTIONS OF 3RD KIND ARE DEFINED AS: H(L)=BJ(L)+I*Y(L)
!     ------------------------------------------------------------------
      LOGICAL    LCALL,LH,LJ,LY
      INTEGER    lmax,LMX
      COMPLEX*16 ARG
      COMPLEX*16 BJ(LMX),H(LMX),Y(LMX)
      COMPLEX*16 Z, CY(lmx), i
      real(8) pi, ZR, ZI, FNU, CYR(lmx), CYI(lmx),
     & CWRKR(lmx), CWRKI(lmx)

      INTEGER KODE, N, NMAXD, NZ, IERR
      pi=4.D0*ATAN(1.D0)
      i = (0.0d0, 1.0d0)
      ZR = real(ARG)
      ZI = aimag(ARG)
      FNU = 0.5D0
      KODE=1
      N=lmax+1
!     call BESSEL_OLD(BJ,Y,H,ARG,LMX,LMAX,LJ,LY,LH,LCALL)
      CALL ZBESJ(ZR, ZI, FNU, KODE, N, CYR, CYI, NZ, IERR)
      ! Convert to spherical function
      CY = (CYR+ i*CYI)*sqrt(pi/2.0d0/arg)
      BJ = CY
      CWRKR=0.0d0
      CWRKI=0.0d0
      call ZBESY(ZR, ZI, FNU, KODE, N, CYR, CYI, NZ, CWRKR, CWRKI,
     *                 IERR)
      CY = (CYR+ i*CYI)*sqrt(pi/2.0d0/arg)
      Y = CY
      H=BJ+i*Y
      END

