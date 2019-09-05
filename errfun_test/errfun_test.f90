program errfun_test
    use errfun
    implicit none
    integer, parameter:: dp=kind(0.d0)
    complex(dp) z_arg, z_res
    complex( kind=dp ), parameter ::                                    &
            cone = cmplx( 1.0_dp, 0.0_dp, kind=dp ),                            &
            ci = cmplx( 0.0_dp, 1.0_dp, kind=dp ),                            &
            czero = cmplx( 0.0_dp, 0.0_dp, kind=dp )

    z_arg = 3*cone + 3*czero
    z_res = cerf_lib(z_arg)
    write(*,*) "result: ", z_res
contains
    !=======================================================================
    complex(dp) function cerf_lib(z)
        !     cerf,given complex argument z,provides the complex error function:
        !     w(z)=exp(-z**2)*(1.0-erf(-i*z))
        complex(dp), intent(in) :: z
        cerf_lib = exp(-z**2) * (1.0 - erf_pop(-ci * z))
        return
    end function

end program errfun_test
