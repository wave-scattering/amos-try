      SUBROUTINE gnzbess(omega,lmax,jl,djl,nl,dnl)

*   Variables declared but never referenced:
*        CRY1              CRYMU             
*   Variables used before set:
*       CRJ               CRJP                             
C--------/---------/---------/---------/---------/---------/---------/--
C   ROUTINE TO GENERATE ARRAYS OF THE SPHERICAL BESSEL FUNCTION 
C       AND THEIR DERIVATIVES FOR A COMPLEX ARGUMENT AND REAL ORDER.
C           (RETURNS jl(i*|OM|), nl(i*|OM|), etc. IF OM IMAGINARY)
C                 ================================
C Limitation : bessjy and bessik require real argument on the input
C
C*****************************************************************
C                      >>> FOR OMEGA REAL  <<<
C  Calculated values are compared with Table 10.5 p. 465 of \ct{AS}
C      (The number (-n) in parenthesis there means 10^{-n})
C
C According to (9.1.60) of \ct{AS}:
C
C           |J_\nu(x)|\leq 1           for  \nu\geq 0
C           |J_\nu(x)|\leq 1/\sqrt{2}  for  \nu\geq 1
C One has
C            j_l=\sqrt{\fr{\pi}{2z}}\, J_{l+1/2}       (10.1.1)
C which when combined with (9.1.60) gives bounds on j_l.
C Calls routine bessjy from Numerical Recipies which returns
C the cylindrical Bessel functions and their derivatives of a given
C order. The routine then determines the rest using the stable downwards 
C recurrence relations  for j_\nu and j_\nu' [cf. (10.1.21-22) of \ct{AS}]:
c                 j_{\nu-1} =((\nu+1)/x) j_\nu + j_\nu'
c                 j_{\nu-1}'=((\nu-1)/x) j_{\nu-1} - j_\nu/(2.d0*x).
C The spherical y_n Bessel functions are then recovered using the formulae 
C of \ct{AS}
C                 y_{n} = (-1)^{n+1} j_{-n-1}           (9.1.2)
C*****************************************************************
C                >>> FOR PURELY IMAGINARY OMEGA   <<<
C
C Calls routine bessik from Numerical Recipies which returns
C the modified cylindrical Bessel functions and their derivatives of a given
C order. The routine then determines the rest using the stable downwards 
C recurrence relations for j_\nu and j_\nu' [cf. (10.1.21-22) of \ct{AS}]:
c
C  Calculated values are compared with Table 10.10 p. 473 of \ct{AS}
C      (The number (-n) in parenthesis there means 10^{-n})
C--------------------------
C According to (9.6.27) of \ct{AS}, I_0'(z)=I_1(z)
C Recurrence for I_\nu [cf. (9.6.26) of \ct{AS}] is stored in bessel.bal:
c                 I_{\nu-1} =(\nu/x) I_\nu + I_\nu'
c                 I_{\nu-1}'=((\nu-1)/x) I_{\nu-1} + I_\nu
c According to (9.6.6) of \ct{AS}: I_{-\nu}(z)=I_{\nu}(z);
c K_{-\nu}(z)=K_{\nu}(z).
C The spherical Bessel functions for complex argument are found 
C according to (10.2.2-3) of \ct{AS}:
c
C         j_l(ix)=e^{l\pi i/2} \sqrt{\fr{\pi}{2x}}\, I_{l+1/2}(x)
C      y_l(ix)=e^{-3(l+1)\pi i/2} \sqrt{\fr{\pi}{2x}}\, I_{-l-1/2}(x)
C*****************************************************************
C                >>> FOR GENERAL COMPLEX OMEGA   <<<
C AMOS SLAC package from http://www.netlib.org/amos/ is used
C -----------------------------------------------------------------
      implicit none
      integer LMAX,IKODE,i,l,LMX
      REAL*8 pi,pib2z,prod,fnu
      PARAMETER (PI=3.141592653589793d0)
* LMX is internal LMAX here (for AMOS):
c      PARAMETER (LMX=82)
      PARAMETER (LMX=82)
* IKODE=1 ... unscaled Bessel functions
* IKODE=2 ... scaled Bessel functions (for AMOS):
      PARAMETER (IKODE=1)
* FNU=ORDER OF INITIAL Y FUNCTION, FNU.GE.0.0D0(for AMOS):
      PARAMETER (FNU=0.5d0)
*
      integer nz,ierr
      REAL*8 xom,xnu,rj,rymu,rjp,ry1,xi
      complex*16 JL(0:lmax),NL(0:lmax),DJL(0:lmax),DNL(0:lmax)
      complex*16 ei,omega,xjl,cxi,czi,cpib2z
      complex*16 crj,crymu,crjp,cry1
      real*8 zr,zi
      real*8 cyr(lmx+1),cyi(lmx+1)
c      character*1 yesn

      ei=dcmplx(0.d0,1.d0)
      zi=dimag(omega)
      zr=dble(omega)

*-------------
*  remainders:

      if (lmx.lt.lmax) then
      write(6,*)'In gnzbess.f:'
      write(6,*)'Internal LMAX (LMX).lt.LMAX'
      stop
      end if 

      if ((zi.ne.0.d0).and.(zr.ne.0.d0)) then
c      write(6,*)'OMEGA is neither real nor purely imaginary'
      go to 200
      end if 
*--------------
c      yesn='n' 
      xom=zabs(omega) 
      pib2z=sqrt(pi/(2.d0*xom))
      xnu=dble(lmax+1/2.)
      
*
      if (dimag(omega).ne.0.) go to 100
*
*******************************************************************
C                       >>> OMEGA  REAL  <<<
*                  >>>  ASYMPTOTIC   PART  <<<
*
* security for omega exceedingly small:

      if (zabs(omega).gt.1.d-9) go to 10

* using asymptotic expansion [cf. (10.1.2) of \ct{AS}] to determine
* JL(LMAX) and DJL(LMAX):

      prod=1.d0
        do i=0,LMAX 
         prod=dble(2*i+1)*prod    
        enddo
         JL(LMAX)=(1.d0-omega**2/dble(2*(2*LMAX+3)) )*omega**LMAX
     &   /prod
         XJL=(1.d0-omega**2/dble(2*(2*LMAX+5)) )*omega**(LMAX+1)
     &   /(prod*dble(2*lmax+3))
         xnu=dble(lmax)
         xi=1.d0/xom
         DJL(LMAX)=JL(LMAX)*xnu*xi-XJL
      do i=1,lmax
      prod=prod/dble(2*(LMAX-i)+3)
        JL(LMAX-i)=(1.d0-omega**2/dble(2*(2*(LMAX-i)+3)))*
     &   omega**(LMAX-i)/prod
        xnu=dble(lmax-i)
        DJL(LMAX-i)=JL(LMAX-i)*xnu*xi-JL(LMAX+1-i)       
      enddo

* continuation of the previous loop
      NL(0)=JL(0)*xi+DJL(0)
      DNL(0)=-NL(0)*xi-JL(0)
      NL(0)=-NL(0)
      DNL(0)=-DNL(0)
      do i=1,lmax
        xnu=dble(i-1)
        NL(i)=NL(i-1)*xnu*xi-DNL(i-1)
        xnu=xnu+2.d0
        DNL(i)=-NL(i)*xnu*xi+NL(i-1)
      enddo

      return
**************************************************************
C                      >>>  OMEGA  REAL  <<<
*                    >>>   NORMAL  PART   <<<
*
 10   call bessjy(xom,xnu,rj,rymu,rjp,ry1)

      JL(LMAX)=pib2z*rj
      DJL(LMAX)=pib2z*rjp-JL(LMAX)/(2.d0*xom)

* Stable downwards recurrence for J_\nu and 
c Using (9.1.2) of \ct{AS} one has
c         (\cos[(n+1/2)\pi]=0,       \sin[(n+1/2)\pi]=(-1)^n) 
c                  
c               Y_{n+1/2} = (-1)^{n+1} J_{-(n+1/2)}
c  and hence
c                   y_{n} = (-1)^{n+1} j_{-n-1}
c 
c  and one can continue with a stable recurrence downwards to generate
c  also the relevant Neumann (Weber) functions.
      
      xi=1.d0/xom
      do i=1,lmax
        xnu=dble(lmax-i+2)
        JL(LMAX-i)=JL(LMAX+1-i)*xnu*xi+DJL(LMAX+1-i)
        xnu=xnu-2.d0
        DJL(LMAX-i)=JL(LMAX-i)*xnu*xi-JL(LMAX+1-i)       
      enddo

* continuation of the previous loop
      NL(0)=JL(0)*xi+DJL(0)
      DNL(0)=-NL(0)*xi-JL(0)
      NL(0)=-NL(0)
      DNL(0)=-DNL(0)
      do i=1,lmax
        xnu=dble(i-1)
        NL(i)=NL(i-1)*xnu*xi-DNL(i-1)
        xnu=xnu+2.d0
        DNL(i)=-NL(i)*xnu*xi+NL(i-1)
      enddo
* security trap:
      do L=0,LMAX
        if (abs(JL(L)).gt.pib2z) then
        write(6,*)'GNCBESS gives incorrect Bessel functions'
        write(6,*)'L=', L, ', JL(L)=', JL(L),'.gt.',pib2z
        write(6,*)'Violation of the bound (9.1.60) of {AS}'
        stop
        end if
      enddo
      return
**************************************************************
*
 100  continue
*
C                   >>>  PURELY IMAGINARY OMEGA   <<<
*                     >>>   ASYMPTOTIC   PART  <<<
*
* security for omega exceedingly small:

      if (zabs(omega).gt.1.d-9) go to 110

* using asymptotic expansion [cf. (10.1.2) of \ct{AS}] to determine
* JL(LMAX) and DJL(LMAX):

      prod=1.d0
        do i=0,LMAX 
         prod=dble(2*i+1)*prod    
        enddo
         JL(LMAX)=(1.d0-omega**2/dble(2*(2*LMAX+3)) )*omega**LMAX
     &   /prod
         XJL=(1.d0-omega**2/dble(2*(2*LMAX+5)) )*omega**(LMAX+1)
     &   /(prod*dble(2*lmax+3))
         xnu=dble(lmax)
         cxi=-ei*1.d0/xom
         DJL(LMAX)=JL(LMAX)*xnu*cxi-XJL
      do i=1,lmax
      prod=prod/dble(2*(LMAX-i)+3)
      JL(LMAX-i)=(1.d0-omega**2/dble(2*(2*(LMAX-i)+3)))*
     &   omega**(LMAX-i)/prod
        xnu=dble(lmax-i)
        DJL(LMAX-i)=JL(LMAX-i)*xnu*cxi-JL(LMAX+1-i)       
      enddo

* continuation of the previous loop
      NL(0)=JL(0)*cxi+DJL(0)
      DNL(0)=-NL(0)*cxi-JL(0)
      NL(0)=-NL(0)
      DNL(0)=-DNL(0)
      do i=1,lmax
        xnu=dble(i-1)
        NL(i)=NL(i-1)*xnu*cxi-DNL(i-1)
        xnu=xnu+2.d0
        DNL(i)=-NL(i)*xnu*cxi+NL(i-1)
      enddo
      return

*----------------------------------------------------------------------
C                   >>>  PURELY IMAGINARY OMEGA   <<<
*                     >>>   NORMAL   PART  <<<
*    
 110  call bessik(xom,xnu,rj,rymu,rjp,ry1)

      JL(LMAX)=pib2z*rj*ei**lmax
      DJL(LMAX)=-ei*pib2z*rjp*ei**lmax+ei*JL(LMAX)/(2.d0*xom)
      
      cxi=-ei*1.d0/xom
      do i=1,lmax
        xnu=dble(lmax-i+2)
        JL(LMAX-i)=JL(LMAX+1-i)*xnu*cxi+DJL(LMAX+1-i)
        xnu=xnu-2.d0
        DJL(LMAX-i)=JL(LMAX-i)*xnu*cxi-JL(LMAX+1-i)       
      enddo

* continuation of the previous loop
      NL(0)=JL(0)*cxi+DJL(0)
      DNL(0)=-NL(0)*cxi-JL(0)
      NL(0)=-NL(0)
      DNL(0)=-DNL(0)
      do i=1,lmax
        xnu=dble(i-1)
        NL(i)=NL(i-1)*xnu*cxi-DNL(i-1)
        xnu=xnu+2.d0
        DNL(i)=-NL(i)*xnu*cxi+NL(i-1)
      enddo
      return

*----------------------------------------------------------------------
C                   >>>  GENERAL COMPLEX OMEGA   <<<

 200  czi=dcmplx(1.d0,0.d0)/omega
      cpib2z=sqrt(czi)*sqrt(pi/2.d0)

*                  >>>  ASYMPTOTIC   PART  <<<
*
* security for omega exceedingly small:

      if (zabs(omega).gt.1.d-9) go to 275   !250 --> Bess, 
                                            !275 --> Sandia, 300 --> NR

* using asymptotic expansion [cf. (10.1.2) of \ct{AS}] to determine
* JL(LMAX) and DJL(LMAX):

      prod=1.d0
        do i=0,LMAX 
         prod=dble(2*i+1)*prod    
        enddo
         JL(LMAX)=(1.d0-omega**2/dble(2*(2*LMAX+3)) )*omega**LMAX
     &   /prod
         XJL=(1.d0-omega**2/dble(2*(2*LMAX+5)) )*omega**(LMAX+1)
     &   /(prod*dble(2*lmax+3))
         xnu=dble(lmax)
         DJL(LMAX)=JL(LMAX)*xnu*czi-XJL
      do i=1,lmax
      prod=prod/dble(2*(LMAX-i)+3)
        JL(LMAX-i)=(1.d0-omega**2/dble(2*(2*(LMAX-i)+3)))*
     &   omega**(LMAX-i)/prod
        xnu=dble(lmax-i)
        DJL(LMAX-i)=JL(LMAX-i)*xnu*czi-JL(LMAX+1-i)       
      enddo

* continuation of the previous loop
      NL(0)=JL(0)*czi+DJL(0)
      DNL(0)=-NL(0)*czi-JL(0)
      NL(0)=-NL(0)
      DNL(0)=-DNL(0)
      do i=1,lmax
        xnu=dble(i-1)
        NL(i)=NL(i-1)*xnu*czi-DNL(i-1)
        xnu=xnu+2.d0
        DNL(i)=-NL(i)*xnu*czi+NL(i-1)
      enddo

      return

*------------------------------
*                    >>>   NORMAL   PART  <<<

 250  call bess(lmax,omega,jl,djl)

	do i=lmax,1,-1
	jl(i)=jl(i-1)
	djl(i)=djl(i-1)
	enddo

      do i=lmax,lmax
        xnu=dble(lmax-i+2)
* (10.1.21) of AS:
        JL(LMAX-i)=xnu*czi*JL(LMAX-i+1) + djl(LMAX-i+1)
        xnu=xnu-2.d0
* (10.1.22) of AS:
        DJL(LMAX-i)=JL(LMAX-i)*xnu*czi-JL(LMAX+1-i)       
      enddo

* continuation of the previous loop
* (10.1.21) of AS calculating of JL(-1) and DJL(-1):
      NL(0)=JL(0)*czi+DJL(0)
      DNL(0)=-NL(0)*czi-JL(0)
*  (10.1.15) of AS: NL(0)=JL(-1) : sign correction
      NL(0)=-NL(0)
      DNL(0)=-DNL(0)
      do i=1,lmax
        xnu=dble(i-1)
* (10.1.22) of AS:
        NL(i)=NL(i-1)*xnu*czi-DNL(i-1)
        xnu=xnu+2.d0
* (10.1.21) of AS:
        DNL(i)=-NL(i)*xnu*czi+NL(i-1)
      enddo
      RETURN
*                    >>>   NORMAL   PART  SANDIA <<<

 275  call zbesj(zr,zi,fnu,ikode,lmax+1,cyr,cyi,nz,ierr)
      if(ierr.gt.0) then
      write(6,*)'ERROR in ZBESJ. IERR=', ierr   
      write(6,*)'IERR=0, NORMAL RETURN - COMPUTATION COMPLETED'
C--------/---------/---------/---------/---------/---------/---------/--
       if(ierr.eq.1) then  
      write(6,*)'IERR=1, INPUT ERROR-NO COMPUTATION'
      else if(ierr.eq.2) then  
      write(6,*)'IERR=2, OVERFLOW - NO COMPUTATION, AIMAG(Z) TOO LARGE 
     1 ON KODE=1'      
        else if(ierr.eq.3) then  
      write(6,*)'IERR=3, CABS(Z) OR FNU+N-1 LARGE - COMPUTATION DONE       
     1 BUT LOSSES OF SIGNIFCANCE BY ARGUMENT REDUCTION PRODUCE LESS THAN 
     2 HALF OF MACHINE ACCURACY' 
      else if(ierr.eq.4) then  
      write(6,*)'IERR=4, CABS(Z) OR FNU+N-1 TOO LARGE - NO COMPUTATION
     1  BECAUSE OF COMPLETE LOSSES OF SIGNIFIANCE BY ARGUMENT REDUCTION' 
      else if(ierr.eq.5) then  
      write(6,*)'IERR=5, ERROR - NO COMPUTATION
     1 ALGORITHM TERMINATION CONDITION NOT MET'
      end if 
      stop
      end if 

      JL(LMAX)=cpib2z*dcmplx(cyr(LMAX+1),cyi(LMAX+1))
      JL(LMAX-1)=cpib2z*dcmplx(cyr(LMAX),cyi(LMAX))
* (10.1.21) of AS:
      DJL(LMAX)=JL(LMAX-1)-dble(lmax+1)*JL(LMAX)*czi
*  
      do i=1,lmax
        xnu=dble(lmax-i+2)
        JL(LMAX-i)=cpib2z*dcmplx(cyr(LMAX-i+1),cyi(LMAX-i+1))
        xnu=xnu-2.d0
* (10.1.22) of AS:
        DJL(LMAX-i)=JL(LMAX-i)*xnu*czi-JL(LMAX+1-i)       
      enddo

* continuation of the previous loop
* (10.1.21) of AS calculating of JL(-1) and DJL(-1):
      NL(0)=JL(0)*czi+DJL(0)
      DNL(0)=-NL(0)*czi-JL(0)
*  (10.1.15) of AS: NL(0)=JL(-1) : sign correction
      NL(0)=-NL(0)
      DNL(0)=-DNL(0)
      do i=1,lmax
        xnu=dble(i-1)
* (10.1.22) of AS:
        NL(i)=NL(i-1)*xnu*czi-DNL(i-1)
        xnu=xnu+2.d0
* (10.1.21) of AS:
        DNL(i)=-NL(i)*xnu*czi+NL(i-1)
      enddo
      RETURN
*-------------------------  UNDER CONSTRUCTION  -------------------------
*                     >>>   NORMAL   PART  NR <<<
 300      xnu=dble(lmax+1/2.)
C      call zbessjy(omega,xnu,crj,crymu,crjp,cry1)

      JL(LMAX)=cpib2z*crj
      DJL(LMAX)=cpib2z*crjp-JL(LMAX)*czi/2.d0

* Stable downwards recurrence for J_\nu and 
c Using (9.1.2) of \ct{AS} one has
c         (\cos[(n+1/2)\pi]=0,       \sin[(n+1/2)\pi]=(-1)^n) 
c                  
c               Y_{n+1/2} = (-1)^{n+1} J_{-(n+1/2)}
c  and hence
c                   y_{n} = (-1)^{n+1} j_{-n-1}
c 
c  and one can continue with a stable recurrence downwards to generate
c  also the relevant Neumann (Weber) functions.
      
      do i=1,lmax
        xnu=dble(lmax-i+2)
        JL(LMAX-i)=JL(LMAX+1-i)*xnu*czi+DJL(LMAX+1-i)
        xnu=xnu-2.d0
        DJL(LMAX-i)=JL(LMAX-i)*xnu*czi-JL(LMAX+1-i)       
      enddo

* continuation of the previous loop
      NL(0)=JL(0)*czi+DJL(0)
      DNL(0)=-NL(0)*czi-JL(0)
      NL(0)=-NL(0)
      DNL(0)=-DNL(0)
      do i=1,lmax
        xnu=dble(i-1)
        NL(i)=NL(i-1)*xnu*czi-DNL(i-1)
        xnu=xnu+2.d0
        DNL(i)=-NL(i)*xnu*czi+NL(i-1)
      enddo
* security trap: here use (9.1.62)
c      do L=0,LMAX
c        if (abs(JL(L)).gt.pib2z) then
c        write(6,*)'GNZBESS gives incorrect Bessel functions'
c        write(6,*)'L=', L, ', JL(L)=', JL(L),'.gt.',pib2z
c        write(6,*)'Violation of the bound (9.1.62) of {AS}'
        write(6,*)'Bound (9.1.62) of {AS} not checked'
c        stop
c        end if
c      enddo
      return

* End of evaluation 
      END
C (C) Copr. 7/2000  Alexander Moroz
