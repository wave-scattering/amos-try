      SUBROUTINE bessjy(x,xnu,rj,rymu,rjp,ry1)
c Variables declared but never referenced:
c        RY                RYP             RYTEMP
c--------/---------/---------/---------/---------/---------/---------/--
c >>> x,xnu
c <<< j,ry,rjp,ryp
c                ==========================================
c       Calculates the cylindrical Bessel functions of arbitrary order.
c Returns rj=j_\nu, rjp=J_\nu', rymu=Y_\mu , ry1=Y_{\mu+1}
c       Applicable in the interval \approx (0, 10.000)
C       The spherical Bessel functions can be found using
c               j_l=\sqrt{\fr{\pi}{2z}}\, J_{n+1/2}, etc.
c
c  TESTED UP TO XNU=19 and X=50 - WORKS WITH AN EXCELLENT PRECISION !!!
c                ==========================================
c Recurrence relations for the Bessel functions are stable 
c     1) downwards for J_\nu and J_\nu':
c                 J_{\nu-1} =(\nu/x) J_\nu + J_\nu'
c                 J_{\nu-1}'=((\nu-1)/x) J_{\nu-1} - J_\nu'
c     2) upwards for Y_\nu and Y_\nu':
c                 Y_{\nu+1} =(2\nu/x) Y_\nu - Y_{\nu-1}
c                 Y_{\nu}'=(\nu/x) Y_{\nu} - Y_{\nu+1} 
c                                                 [cf. (9.1.27) of \ct{AS}]
c
c Uses two recurrence relations, 
c    - CF1 which determines the ratio J_\nu/J_\nu'
c    - CF2 which determines the ratio p+iq=(J_\nu'+iY_\nu')/(J_\nu+iY_\nu)
c and the expression for the Wronskian 
c    - W= J_\nu Y_\nu' - J_\nu' Y_\nu = 2/(\pi x)
c
c The rate of convergence of CF1 is determined by the position of the
c turning point x_{tp}=\sqrt{\nu(\nu+1)}, beyond which the Bessel functions
c become oscillatory. If x\leq x_{tp}, convergence is very rapid.
c On the other hand, CF2 converges rapidly for x\geq x_{tp}, while the
c convergence fails for x\rar 0.
c
c  COMMENT: Using (9.1.2) of \ct{AS} one has
c         (\cos[(n+1/2)\pi]=0,       \sin[(n+1/2)\pi]=(-1)^n) 
c                  
c               Y_{n+1/2} = (-1)^{n+1} J_{-(n+1/2)}
c  and hence
c                   y_{n} = (-1)^{n+1} j_{-n-1}
c 
c  and one can continue with a stable recurrence downwards to generate
c  also the relevant Neumann (Weber) functions.
c
c                        E X E C U T I O N 
c
c 1) Uses CF1 recurrence relation to calculate the ratio 
c                     f_\nu=J_\nu/J_\nu'
c
c 2) For x<x_{min}=2 one uses downwards recurrence relations down to
c   |xmu|<1/2, while for x>x_{min}=2 one uses downwards recurrence 
c   relations so that x\geq x_{tp} with the initial value of 
c                       J_\nu=1.d-30
c
c 3) Exaluates CF2 for xmu. Then with W one has enough relations to solve
c    for all four functions and  determine them at xmu.
c    The sign of J_\xmu is chosen to be the same as the sign of the 
c    initial J_\nu. For x<2, Temme's series are used instead of CF2.
c
c 4) Once all four functions are determined at xmu, J_\nu and J_\nu'
c    are determined by simply scaling the initial values by the ratio
c    of CF2+W result to that found after applying the recurrence 1).
c    The quantities Y_\nu and Y_\nu' are found using the stable upwards
c    recurrence 2).
c
c Modified to a double precision subroutine Aug. 17, 1998
c--------/---------/---------/---------/---------/---------/---------/--
      INTEGER MAXIT
      DOUBLE PRECISION rj,rjp,ry,ryp,x,xnu,XMIN
      DOUBLE PRECISION EPS,FPMIN,PI
      PARAMETER (EPS=1.d-16,FPMIN=1.d-30,MAXIT=10000,XMIN=2.d0,
     *PI=3.141592653589793d0)
CU    USES beschb
      INTEGER i,isign,l,nl
      DOUBLE PRECISION a,b,br,bi,c,cr,ci,d,del,del1,den,di,dlr,dli,dr,e,
     *f,fact,fact2,fact3,ff,gam,gam1,gam2,gammi,gampl,h,p,pimu,pimu2,q,
     *r,rjl,rjl1,rjmu,rjp1,rjpl,rjtemp,ry1,rymu,rymup,rytemp,sum,sum1,
     *temp,w,x2,xi,xi2,xmu,xmu2
      if(x.le.0..or.xnu.lt.0.) pause 'bad arguments in bessjy'
      if(x.lt.XMIN)then
        nl=int(xnu+.5d0)
      else
        nl=max(0,int(xnu-x+1.5d0))
      endif
      xmu=xnu-nl
      xmu2=xmu*xmu
      xi=1.d0/x
      xi2=2.d0*xi
      w=xi2/PI
      isign=1
      h=xnu*xi
      if(h.lt.FPMIN)h=FPMIN
      b=xi2*xnu
      d=0.d0
      c=h
      do 11 i=1,MAXIT
        b=b+xi2
        d=b-d
        if(abs(d).lt.FPMIN)d=FPMIN
        c=b-1.d0/c
        if(abs(c).lt.FPMIN)c=FPMIN
        d=1.d0/d
        del=c*d
        h=del*h
* h=J_\nu/J_\nu' is the result of CF1 relation  
        if(d.lt.0.d0)isign=-isign
        if(abs(del-1.d0).lt.EPS)goto 1
11    continue
      pause 'x too large in bessjy; try asymptotic expansion'
1     continue
* Initialization of J_\nu and J_\nu' for the stable downwards recurrence 1)
      rjl=isign*FPMIN
      rjpl=h*rjl
      rjl1=rjl
      rjp1=rjpl
      fact=xnu*xi
      do 12 l=nl,1,-1
        rjtemp=fact*rjl+rjpl
        fact=fact-xi
        rjpl=fact*rjtemp-rjl
        rjl=rjtemp
12    continue
      if(rjl.eq.0.d0)rjl=EPS
      f=rjpl/rjl
* Now have unnormalized J_\mu and J_\mu' and f=J_\mu/J_\mu' 
***********************************************************
* Dichotomy 0<x<2 or x\geq 2
      if(x.lt.XMIN) then
        x2=.5d0*x
        pimu=PI*xmu
        if(abs(pimu).lt.EPS)then
          fact=1.d0
        else
          fact=pimu/sin(pimu)
        endif
        d=-log(x2)
        e=xmu*d
        if(abs(e).lt.EPS)then
          fact2=1.d0
        else
          fact2=sinh(e)/e
        endif
        call beschb(xmu,gam1,gam2,gampl,gammi)
        ff=2.d0/PI*fact*(gam1*cosh(e)+gam2*fact2*d)
        e=exp(e)
        p=e/(gampl*PI)
        q=1.d0/(e*PI*gammi)
        pimu2=0.5d0*pimu
        if(abs(pimu2).lt.EPS)then
          fact3=1.d0
        else
          fact3=sin(pimu2)/pimu2
        endif
        r=PI*pimu2*fact3*fact3
        c=1.d0
        d=-x2*x2
        sum=ff+r*q
        sum1=p
        do 13 i=1,MAXIT
          ff=(i*ff+p+q)/(i*i-xmu2)
          c=c*d/i
          p=p/(i-xmu)
          q=q/(i+xmu)
          del=c*(ff+r*q)
          sum=sum+del
          del1=c*p-i*del
          sum1=sum1+del1
          if(abs(del).lt.(1.d0+abs(sum))*EPS)goto 2
13      continue
        pause 'bessy series failed to converge'
2       continue
        rymu=-sum
        ry1=-sum1*xi2
        rymup=xmu*xi*rymu-ry1
        rjmu=w/(rymup-f*rymu)
      else
        a=.25d0-xmu2
        p=-.5d0*xi
        q=1.d0
        br=2.d0*x
        bi=2.d0
        fact=a*xi/(p*p+q*q)
        cr=br+q*fact
        ci=bi+p*fact
        den=br*br+bi*bi
        dr=br/den
        di=-bi/den
        dlr=cr*dr-ci*di
        dli=cr*di+ci*dr
        temp=p*dlr-q*dli
        q=p*dli+q*dlr
        p=temp
        do 14 i=2,MAXIT
          a=a+2*(i-1)
          bi=bi+2.d0
          dr=a*dr+br
          di=a*di+bi
          if(abs(dr)+abs(di).lt.FPMIN)dr=FPMIN
          fact=a/(cr*cr+ci*ci)
          cr=br+cr*fact
          ci=bi-ci*fact
          if(abs(cr)+abs(ci).lt.FPMIN)cr=FPMIN
          den=dr*dr+di*di
          dr=dr/den
          di=-di/den
          dlr=cr*dr-ci*di
          dli=cr*di+ci*dr
          temp=p*dlr-q*dli
          q=p*dli+q*dlr
          p=temp
          if(abs(dlr-1.d0)+abs(dli).lt.EPS)goto 3
14      continue
        pause 'cf2 failed in bessjy'
3       continue
********************************************************
* determination of J_\mu, J_\mu', Y_\mu, and Y_\mu'
        gam=(p-f)/q
        rjmu=sqrt(w/((p-f)*gam+q))
        rjmu=sign(rjmu,rjl)
        rymu=rjmu*gam
        rymup=rymu*(p+q/gam)
        ry1=xmu*xi*rymu-rymup
* ry1=Y_{\mu+1} here
      endif
* scaling to obtain the original  J_\nu, J_\nu'
      fact=rjmu/rjl
      rj=rjl1*fact
      rjp=rjp1*fact
c* the stable upwards recurrence 2)  (xi2=2/x)
c      do 15 i=1,nl
c        rytemp=(xmu+i)*xi2*ry1-rymu
c        rymu=ry1
c        ry1=rytemp
c15    continue
c      ry=rymu
c      ryp=xnu*xi*rymu-ry1
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software v1]3"iA#.
