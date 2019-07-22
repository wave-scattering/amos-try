!$Id: errtest.f90 10 2016-08-15 19:43:56Z mexas $
! Copyright (c) 2016 Anton Shterenlikht, The University of Bristol, UK
!
! Simple tests of error functions in module errfun.
!*********************************************************************72

program errtest
use errfun
implicit none

! Use np points along x and y.
integer, parameter :: npx = 79, npy = 61 
real( kind=rkind ), parameter ::                                       &
  zero = 0.0_rkind,   one = 1.0_rkind,                                 &
  xmin = -3.9_rkind, xmax = 3.9_rkind,                                 &
  ymin = -3.0_rkind, ymax = 3.0_rkind,                                 &
  xdelta = 0.1_rkind, ydelta = 0.1_rkind
!  xdelta = (xmax-xmin)/(npx-1), ydelta = (ymax-ymin)/(npy-1)
complex( kind=rkind ) :: z(npx,npy), res1(npx,npy),                    &
  res2(npx,npy), z1
real( kind=rkind ) :: x(npx), y(npy), err = 0.06447_rkind * eps0,      &
  re1(npx,npy), im1(npx,npy), re2(npx,npy), im2(npx,npy)
real :: time2, time1

integer :: i,j,funit

! Generate the coordinate arrays
do concurrent( i=1:npx )
  x(i) = xmin + xdelta * (i-1)
end do
!x(40) = zero

do concurrent( i=1:npy )
  y(i) = ymin + ydelta * (i-1)
end do

write (*,*) "xdelta, ydelta:", xdelta, ydelta
write (*,*) "xmax, ymax:", x(npx), y(npy)

! Generate the array of z
do concurrent (i=1:npx, j=1:npy)
  z(i,j) = cmplx( x(i), y(j), kind=rkind ) 
end do

! Time the calculation
call cpu_time( time1 )
res1 = erf_zag( z, err )
call cpu_time( time2 )
write (*,*) "CPU time, s:", time2-time1

call cpu_time( time1 )
res2 = erf_pop( z )
call cpu_time( time2 )
write (*,*) "CPU time, s:", time2-time1


z1 = cmplx( zero, -3.0_rkind, kind=rkind )
write (*,*) erf_zag( z1, err )
write (*,*) erf_pop( z1 )
z1 = cmplx( zero, 3.0_rkind, kind=rkind )
write (*,*) erf_zag( z1, err )
write (*,*) erf_pop( z1 )

writeout: if ( .true. ) then
  ! Write the data out
  re1 = real( res1 )
  im1 = aimag( res1 )
  re2 = real( res2 )
  im2 = aimag( res2 )
  
  open(newunit=funit, file="out-fortran", form="formatted",              &
       status="replace")
  do i=1,npx
  do j=1,npy
    write (funit, "(6(es25.16,tr1))") x(i), y(j), re1(i,j), im1(i,j), &
         re2(i,j), im2(i,j)
  end do
  end do
  close(funit)

end if writeout

end program errtest
