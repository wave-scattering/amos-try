program test_besselj
  use amos
  implicit none
  integer, parameter :: dp=kind(0.d0)
!  integer, parameter :: dp = selected_real_kind(15, 307)
  integer :: i
  complex(dp) :: res
  real(dp) :: ZR, ZI, FNU, Y
  integer :: KODE, NZ, IERR
  integer, parameter :: N=1
  real(dp),dimension(N) :: CYR, CYI
  ZR = 500.0
  ZI = 500.0
!  Z = CMPLX(ZR,ZI)
  KODE =1
  FNU = 50.0
  call ZBESJ(ZR, ZI, FNU, KODE, N, CYR, CYI, NZ, IERR)

  write (*,*) cyr, cyi

end program test_besselj
