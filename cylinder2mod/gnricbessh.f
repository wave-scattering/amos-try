      subroutine GNRICBESSH(CQEPS,XX,LMX,zeta,dzeta)
C--------/---------/---------/---------/---------/---------/---------/--
C >>> CQEPS,XX,LMX
C <<< zeta,dzeta of argument CQEPS*XX from l=0 up to l=LMX
C =====
C Calculates an array of Riccati-Hankel functions of the 
C argument CQEPS*XX using recurrences in Eqs. (63)-(64), (66)-(67) of
C 
C [1] D. W. Mackowski, R. A. Altenkirch, and M. P. Menguc, 
C "Internal absorption cross sections in a stratified sphere," 
C Appl. Opt. 29, 1551-1559 (1990)
C http://www.opticsinfobase.org/ao/abstract.cfm?URI=ao-29-10-1551 
C
C z1     ... \psi_l\xi_l
C DR1(l) ... DR^{(1)}_l=\psi_l'/\psi_l
C DR3(l) ... DR^{(3)}_l= \xi_l'/\xi_l
C
C An array DR1(l) is supplied by biga.f, wherein it is 
C calculated by the downward recurrence relation [e.g. (62) of [1])
C--------/---------/---------/---------/---------/---------/---------/--
      implicit none
      integer LMAXD

C  >>> ANGULAR MOMENTUM CUTOFF ON ARRAY DIMENSION
C      (The actual angular-momentum cutoff on summation is specified
C       by the value of variable LMX) 

      PARAMETER (LMAXD=60)

      integer lmx,l
      real*8 xx
      complex*16 zx,z1,ci,cone,CQEPS
      complex*16 dr1(0:LMAXD),dr3(0:LMAXD),zeta(0:lmx),dzeta(0:lmx)

      DATA ci/(0.d0,1.d0)/,cone/(1.d0,0.d0)/ 
C--------/---------/---------/---------/---------/---------/---------/--

C Determine the Bessel function ratio A_1=\psi_1'/\psi_1

      if (imag(CQEPS).ne.0.d0) then 
         call BIGA(CQEPS,xx,LMX,.false.,.false.,dr1(1))
      else if (imag(CQEPS).eq.0.d0) then 
         call BIGA(CQEPS,xx,LMX,.true.,.false.,dr1(1))
      end if

        zx=xx*cqeps

* initial values of the recurrence (63),(64),(67) of [1]
* \psi_0=zx*j_0=sin(zx); \xi_0=zx*h_0=-ci*exp(ci*zx)

	  z1=-ci*exp(ci*zx)*sin(zx)   
	  DR3(0)=ci
	  DR1(0)=cone/zx - cone/(cone/zx + dr1(1))
	  zeta(0)=-ci*exp(ci*zx) 

*
	  do 15  l=1,lmx   !Eqs. (63),(64),(67) of Mackowski et al

	   z1=z1*(-DR1(l-1)+dble(l)/zx)*(-DR3(l-1)+dble(l)/zx)
	   DR3(l)=DR1(l)+ci/z1
           zeta(l)=zeta(l-1)*(-DR3(l-1)+dble(l)/zx)
           dzeta(l)=dr3(l)*zeta(l)

  15     continue

      return
      end