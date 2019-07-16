      PROGRAM TESTBESJ
      COMPLEX Z, CY
      real(8) ZR, ZI, FNU, CYR, CYI, Y
      INTEGER KODE, N, NZ, IERR
      parameter(N=1)
      dimension CYR(N), CYI(N), CY(N)
      ZR = 500.d0
      ZI = 500.d0
      Z = CMPLX(ZR,ZI)
      KODE =1
      FNU = 50.d0
      CALL ZBESJ(ZR, ZI, FNU, KODE, N, CYR, CYI, NZ, IERR)
!     CALL CBESJ(Z, FNU, KODE, N, CY, NZ, IERR)
      Y= AIMAG(Z)
      write(*,*)CYR(1), CYI(1), IERR
!     write(*,*)CYR(1)/exp(-abs(Y)),CYI(1)/exp(-abs(Y)), IERR

!     write(*,*)CY(1), IERR
!     write(*,*)CY(1)/exp(-abs(Y)), IERR

      END
