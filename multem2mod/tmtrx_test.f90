program tmtrx_test
    use libmultem2b
    implicit none
    complex(dp) epsmed, epssph, mumed, musph, kappa0, rap
    real(dp) :: zinf, zsup, zstep, zval
    integer :: np, i
    complex(dp) :: TE(1), TH(1)

    rap = cmplx_dp(0.21052632_dp, 0.0_dp)
    epssph = cmplx_dp(-20.148_dp, 1.247_dp)
    epsmed= cmplx_dp(2.1025_dp, 0.0_dp)
    mumed = cmplx_dp(1.0_dp, 0.0_dp)
    musph = cmplx_dp(1.0_dp, 0.0_dp)

    call tmtrx(rap, epssph, epsmed, mumed, musph, TE, TH)

    print *, 'TE'
!    print *, (TE(i),i=1,3)
!    print *, (TE(i),i=4,6)
!    print *, (TE(i),i=7)
    print *, TE
    print *, count(TE /= 0)
    print *, '-------------------------------------------------------------'
    print *, 'TH'
!    print *, (TH(i),i=1,3)
!    print *, (TH(i),i=4,6)
!    print *, (TH(i),i=7)
    print *, TH
    print *, '-------------------------------------------------------------'
    print *, count(TH /= 0)

!    multipole_order = (/ 1, 1, 2/)
!    multipole_type = (/ 0, 1, 0/)

    !if count(multipole_order == -1) then
        !if count(TE /= 0) /= max(multipole_order)
            !stop

end program tmtrx_test