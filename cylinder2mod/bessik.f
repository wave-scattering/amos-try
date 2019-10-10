      SUBROUTINE bessik(x,xnu,ri,rkmu,rip,rkmup)
C--------/---------/---------/---------/---------/---------/---------/--
c >>> x,xnu
c <<< ri,rk,rip,rkp
c                ==========================================
c           Calculates the modified Bessel functions of arbitrary order
c Returns ri=I_\nu, rip=I_\nu', rk=K_\mu, rkmup=K_\mu'
c                ==========================================
c Recurrence relations for the Bessel functions are stable 
c     1) downwards for J_\nu and J_\nu':
c                 J_{\nu-1} =(\nu/x) J_\nu + J_\nu'
c                 J_{\nu-1}'=((\nu-1)/x) J_{\nu-1} - J_\nu'
c     2) upwards for Y_\nu and Y_\nu':
c                 Y_{\nu+1} =(2\nu/x) Y_\nu - Y_{\nu-1}
c                 Y_{\nu}'=(\nu/x) Y_{\nu} - Y_{\nu+1} 
C-----------------------------------------------------------------------
      INTEGER MAXIT
      REAL*8 ri,rip,rk,rkp,x,xnu,XMIN
      DOUBLE PRECISION EPS,FPMIN,PI
      PARAMETER (EPS=1.d-16,FPMIN=1.d-30,MAXIT=10000,XMIN=2.d0,
     *PI=3.141592653589793d0)
CU    USES beschb
      INTEGER i,l,nl
      DOUBLE PRECISION a,a1,b,c,d,del,del1,delh,dels,e,f,fact,fact2,ff,
     *gam1,gam2,gammi,gampl,h,p,pimu,q,q1,q2,qnew,ril,ril1,rimu,rip1,
     *ripl,ritemp,rk1,rkmu,rkmup,rktemp,s,sum,sum1,x2,xi,xi2,xmu,xmu2
      if(x.le.0..or.xnu.lt.0.) pause 'bad arguments in bessik'
      nl=int(xnu+.5d0)
      xmu=xnu-nl
      xmu2=xmu*xmu
      xi=1.d0/x
      xi2=2.d0*xi
      h=xnu*xi
      if(h.lt.FPMIN)h=FPMIN
      b=xi2*xnu
      d=0.d0
      c=h
      do 11 i=1,MAXIT
        b=b+xi2
        d=1.d0/(b+d)
        c=b+1.d0/c
        del=c*d
        h=del*h
        if(abs(del-1.d0).lt.EPS)goto 1
11    continue
      pause 'x too large in bessik; try asymptotic expansion'
1     continue
      ril=FPMIN
      ripl=h*ril
      ril1=ril
      rip1=ripl
      fact=xnu*xi
      do 12 l=nl,1,-1
        ritemp=fact*ril+ripl
        fact=fact-xi
        ripl=fact*ritemp+ril
        ril=ritemp
12    continue
      f=ripl/ril
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
        ff=fact*(gam1*cosh(e)+gam2*fact2*d)
        sum=ff
        e=exp(e)
        p=0.5d0*e/gampl
        q=0.5d0/(e*gammi)
        c=1.d0
        d=x2*x2
        sum1=p
        do 13 i=1,MAXIT
          ff=(i*ff+p+q)/(i*i-xmu2)
          c=c*d/i
          p=p/(i-xmu)
          q=q/(i+xmu)
          del=c*ff
          sum=sum+del
          del1=c*(p-i*ff)
          sum1=sum1+del1
          if(abs(del).lt.abs(sum)*EPS)goto 2
13      continue
        pause 'bessk series failed to converge'
2       continue
        rkmu=sum
        rk1=sum1*xi2
      else
        b=2.d0*(1.d0+x)
        d=1.d0/b
        delh=d
        h=delh
        q1=0.d0
        q2=1.d0
        a1=.25d0-xmu2
        c=a1
        q=c
        a=-a1
        s=1.d0+q*delh
        do 14 i=2,MAXIT
          a=a-2*(i-1)
          c=-a*c/i
          qnew=(q1-b*q2)/a
          q1=q2
          q2=qnew
          q=q+c*qnew
          b=b+2.d0
          d=1.d0/(b+a*d)
          delh=(b*d-1.d0)*delh
          h=h+delh
          dels=q*delh
          s=s+dels
          if(abs(dels/s).lt.EPS)goto 3
14      continue
        pause 'bessik: failure to converge in cf2'
3       continue
********************************************************
* determination of I_\mu, I_\mu', K_\mu, and K_\mu'
        h=a1*h
        rkmu=sqrt(PI/(2.d0*x))*exp(-x)/s
        rk1=rkmu*(xmu+x+.5d0-h)*xi
* rk1=K_{\mu+1} here
      endif
      rkmup=xmu*xi*rkmu-rk1
* rkmup=K_{\mu}' 
* scaling to obtain the original  J_\nu, J_\nu'
      rimu=xi/(f*rkmu-rkmup)
      ri=(rimu*ril1)/ril
      rip=(rimu*rip1)/ril
c* the stable upwards recurrence for K_\mu
c      do 15 i=1,nl
c        rktemp=(xmu+i)*xi2*rk1+rkmu
c        rkmu=rk1
c        rk1=rktemp
c15    continue
c      rk=rkmu
c      rkp=xnu*xi*rkmu-rk1
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software v1]3"iA#.
