program errfun_test
    use errfun
    implicit none
    integer, parameter:: dp=kind(0.d0)
    complex(dp) z_arg, z_res_lib, z_res_multem2
    complex( kind=dp ), parameter ::                                    &
            cone = cmplx( 1.0_dp, 0.0_dp, kind=dp ),                            &
            ci = cmplx( 0.0_dp, 1.0_dp, kind=dp ),                            &
            czero = cmplx( 0.0_dp, 0.0_dp, kind=dp )
    real(dp), parameter:: emach = 1.0d-8, step = 0.3
    integer i, j
    integer, parameter:: total =15
    real(dp), parameter:: start = -step*total, stop = step*total
    real(dp) x,y, z_abs
    !=======================================================================

    do x = start, stop, step
        do y = start, stop, step
            z_arg = cmplx(x,y, kind=dp)
            z_res_lib = cerf_lib(z_arg, emach)
            z_res_multem2 = cerf_multem2(z_arg,emach)
            z_abs = dble(abs(z_res_lib/z_res_multem2 -1.0_dp))
            if (z_abs > emach*10) then
                write(*,*) "arg", z_arg,"result: z_lib", z_res_lib, "z_old",z_res_multem2
                write(*,*) "z_lib/z_old - 1.0 = ", z_res_lib/z_res_multem2 -1.0_dp, "  arg:", x, y, "abs:", z_abs

            end if
        end do
    end do

contains
    !=======================================================================
    complex(dp) function cerf_lib(z, emach)
        !     cerf,given complex argument z,provides the complex error function:
        !     w(z)=exp(-z**2)*(1.0-erf(-i*z))
        complex(dp), intent(in) :: z
        real(dp) emach
!        cerf_lib = exp(-z**2) * (1.0 - erf_zag(-ci * z, emach))
        cerf_lib = exp(-z**2) * (1.0 - erf_pop(-ci * z))
!        cerf_lib = exp(-z**2) * (1.0 - erf(-ci * z))
        return
    end function
    !=======================================================================
    function cerf_multem2(z,emach)
        implicit none
        !     ------------------------------------------------------------------
        !     cerf,given complex argument z,provides the complex error function:
        !     w(z)=exp(-z**2)*(1.0-erf(-i*z))
        !     the  evaluation  always  takes   place  in  the  first   quadrant.
        !     one  of  three methods  is  exployed  depending on the size of the
        !     argument (a power series,a recurrence based on continued fractions
        !     theory, or an asymptotic series). emach is the machine accuracy
        !     ------------------------------------------------------------------
        !
        ! ..  scalar arguments  ..
        !
        real*8     emach
        complex*16 z
        !
        ! ..  local scalars  ..
        !
        integer    nn,n
        real*8     absz,abterm,api,eps,fact,factd,factn,pi
        real*8     q,rtpi,test,x,y,yy
        complex*16 zz,cone,ci,czero,sum,zzs,xzzs,cer,cerf_multem2
        complex*16 h1,h2,h3,u1,u2,u3,term1,term2
        !
        ! ..  intrinsic functions  ..
        !
        intrinsic dcmplx,exp,conjg
        !
        ! ..  data statements  ..
        !
        data pi/3.14159265358979d0/
        data cone/(1.d0,0.d0)/,ci/(0.d0,1.d0)/,czero/(0.d0,0.d0)/
        !     ------------------------------------------------------------------
        !
        eps=5.0d0*emach
        api=1.0d0/pi
        if(abs(z))2,1,2
        1  cerf_multem2=cone
        goto 29
        !
        !     the argument is translated to the first quadrant from
        !     the nn_th quadrant, before the method for the function
        !     evaluation is chosen
        !
        2  x=dreal(z)
        y=dimag(z)
        yy=y
        if(y)6,3,3
        3  if(x)5,4,4
        4  zz=z
        nn=1
        goto 9
        5  zz=dcmplx(-x,y)
        nn=2
        goto 9
        6  yy=-y
        if(x)7,8,8
        7  zz=-z
        nn=3
        goto 9
        8  zz=dcmplx(x,-y)
        nn=4
        9  zzs=zz*zz
        xzzs=exp(-zzs)
        absz=abs(zz)
        if(absz-10.0d0)10,10,23
        10  if(yy-1.0d0)11,12,12
        11  if(absz-4.0d0)13,18,18
        12  if(absz-1.0d0)13,18,18
        !
        !     power series(see abramowitz and stegun handbook of
        !     mathematical functions, p297)
        !
        13  q=1.0d0
        factn=-1.0d0
        factd=1.0d0
        term1=zz
        sum=zz
        14 do n=1,5
            factn=factn+2.0d0
            factd=factd+2.0d0
            fact=factn/(q*factd)
            term1=fact*zzs*term1
            sum=sum+term1
            q=q+1.0d0
        end do
        abterm=abs(term1)
        if(abterm-eps)17,16,16
        16  if(q-100.0d0)14,17,17
        17  fact=2.0d0*sqrt(api)
        sum=fact*ci*sum
        cer=xzzs+xzzs*sum
        goto 24
        !
        !     continued fraction theory(w(z) is related to the limiting
        !     value of u(n,z)/h(n,z), where u and h obey the same
        !     recurrence relation in n. see faddeeva and terentiev
        !     (tables of values of w(z) for complex arguments,pergamon
        !       n.y. 1961)
        !
        18  term2=dcmplx(1.d6,0.0d0)
        q=1.0d0
        h1=cone
        h2=2.0d0*zz
        u1=czero
        rtpi=2.0d0*sqrt(pi)
        u2=dcmplx(rtpi,0.0d0)
        19  term1=term2
        do n=1,5
            h3=h2*zz-q*h1
            u3=u2*zz-q*u1
            h1=h2
            h2=2.0d0*h3
            u1=u2
            u2=2.0d0*u3
            q=q+1.0d0
        end do
        term2=u3/h3
        test=abs((term2-term1)/term1)
        if(test-eps)22,21,21
        21  if(q-60.0d0)19,19,13
        22  cer=api*ci*term2
        goto 24
        !
        !     asymptotic series: see abramowitz and stegun, p328
        !
        23  cer=0.5124242d0/(zzs-0.2752551d0)+0.05176536d0/(zzs-2.724745d0)
        cer=ci*zz*cer
        !
        !     symmetry relations are now used to transform the function
        !     back to quadrant nn
        !
        24  goto(28,26,27,25),nn
        25  cer=2.0d0*xzzs-cer
        26  cerf_multem2=conjg(cer)
        goto 29
        27  cerf_multem2=2.0d0*xzzs-cer
        goto 29
        28  cerf_multem2=cer
        29  return
    end
    !=======================================================================


end program errfun_test
