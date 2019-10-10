      subroutine sordalc(NMAT,LAMBDA,ZEPS)
C----------------------------------------------------------------------
C        SUBROUTINE TO CALCULATE THE DIELECTRIC CONSTANT OF METALS
C             ACCORDING TO AN ARTICLE BY ORDAL et al
C   [Appl. Opt. {\bf 22}, 1099 (1983); ibid. {\bf 24}, 4493 (1983)]
C                   IN THE INFRARED AND FAR INFRARED
C
C             f77 -g -check_bounds ordalc.f -o rnordalc
C
C  !!!  ALWAYS ADJUST LN PARAMETER FROM THE DATA WHICH ARE READ IN !!!
C          
C   omega=1/\lambda [cm^{-1}]  in Ordal [spectroscopic convention]
C
C   Re eps = - \fr{\om_p^2}{\om^2+\om_\tau^2}
C   Im eps =   \fr{\om_p^2 \om_\tau}{\om^3+\om \om_\tau^2}
C   1eV=1.24 \mu m
C----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER NMAT
      REAL*8 lambda,plasma,tau
      COMPLEX*16 zeps
C                       -------------------------------
      REAL*8 pi,x,y,omega
      DATA PI/3.141592653589793d0/
C   ---------      
C ::: speed of light in vacuum in nm/s
C      PARAMETER (c0=2.9927925d17)
C
C          plasma[THz]/plasma[cm-1]  tau[THz]/tau[cm-1]    
C   Al            3570/119000            19.4/660          
C   Cu            1914/59600             8.34/73.2          
C   Au            2175/72800              6.5/215           
C   Ag            2175/72700             4.35/145            
C
C ::: conversion factor between normal angular frequency and eV:
c      PARAMETER(XCN=4.13566727d-15)    

      if (nmat.eq.3) then
C ::: plasma frequency in THz:
         PLASMA=72800.d-7            !d12 in Hz/d-7 in [nm-1]
C ::: tau frequency in THz:
          TAU=215.d-7     
      else if (nmat.eq.5) then 
        PLASMA=59600d-7
        TAU=73.2d-7
      else if (nmat.eq.6) then
        PLASMA=119000d-7
        TAU=660.d-7
      end if
  
C                       -------------------------------

c      write(6,*)'Read in wavelength in nm'
c      read(5,*) lambda
       omega=1.d0/lambda
*
       X=-plasma**2/(omega**2+tau**2)
       Y=tau*plasma**2/(omega**3+omega*tau**2)       
*
       zeps=dcmplx(X,Y)  
*
       END
*
C (C) Copr. 04/2002  Alexander Moroz
