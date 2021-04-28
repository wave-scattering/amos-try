module errfun

    use iso_c_binding

    interface
        complex (c_double_complex) function faddeeva_w(z, relerr) bind(c, name = 'Faddeeva_w')
            use iso_c_binding
            complex(c_double_complex), value :: z
            real(c_double), value :: relerr
        end function faddeeva_w
    end interface

end module errfun
