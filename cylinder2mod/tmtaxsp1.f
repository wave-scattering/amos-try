      SUBROUTINE TMTAXSP(nmax,RAP,zeps1,zeps0,TMT) 

c Warning in module TMTAXSP in file tmtaxsp.f: Variables set but never used:
c    NGGG set at line 182 file tmtaxsp.f

c Warning in module TMTAXSP in file tmtaxsp.f: Variables may be used before set:
c    QEXT1 used at line 215 file tmtaxsp.f
c    QEXT1 set at line 220 file tmtaxsp.f
c    QSCA1 used at line 214 file tmtaxsp.f
c    QSCA1 set at line 221 file tmtaxsp.f    

C--------/---------/---------/---------/---------/---------/---------/--
C NMAX - angular momentum cut off
C RAP=S(1,1)*KAPPA0/2.D0/PI     !=rmuf*ALPHA/LAMBDA =rsnm/LAMBDA 
C
C    RETURNS the T matrix of a general axially symmetric scatterer 
C    The latter has also non-zero of-diagonal (mixed EH and HE terms):
C
C                         |  TMT(1,*) |  TMT(4,*)   |
C                 TMT  =  | ----------+-------------|
C                         |  TMT(3,*) |  TMT(2,*)   |
C
C    TMT(1,*) terms corresponds to TEE scattering matrix,
C    TMT(2,*) terms corresponds to TMM scattering matrix
C    TMT(3,*) terms corresponds to TME scattering matrix
C    TMT(4,*) terms corresponds to TEM scattering matrix
C    TMT(4,*)=-TMT(3,*)^t where t denotes transposed TMT(3,*) submatrix
C
C    TMT's equal to i*sin(eta)*exp(i*eta), where eta is a phase-shift
C             ====>    S=1+2*T=exp(2*i*eta)
C    and the unitarity of S-matrix implies
C            T^\dagger*T=T*T^\dagger=-(1/2)*(T^\dagger+T) 
C______________________________
C
C LOCAL VARIABLES:
C ===============
C ICHOICE=1 if NAG library is available, otherwise ICHOICE=2 
C
C NP,EPS: specifies the shape of particles within a given NP class:
C     NP.gt.0 - EPS = deformation parameter of a Chebyshev particle
C     NP=-1 - EPS = the ratio of the horizontal to rotational axes. EPS is 
C             larger than 1 for oblate spheroids and smaller than 1 for        
C             prolate spheroids.
C     NP=-2 - EPS = the ratio of the cylinder diameter to its length.
C     NP=-3 - no EPS is specified
C     NP=-4 - EPS is the height (along the axial symmetry axis) 
C             of the resulting cut sphere
C                Note that always EPS.LT.2*REV specified
C
C Warning:
C  In computations for spheres, use EPS=1.000001 instead of EPS=1.    
C  EPS=1 can cause overflows in some rare cases.  
C
C  LAM - the (vacuum) wavelength of incident light. Changed to
C                       LAM=LAM*SQRT(ZEPS0) here
C
C  RAT = 1 - particle size is specified in terms of the            
C                equal-volume-sphere radius                             
C  RAT.ne.1 - particle size is specified in terms of the           
C                equal-surface-area-sphere radius
C  AXI ... equivalent-(volume/surface-area)-sphere radius
C  REV=A=RAT*AXI ... equal-volume-sphere radius 
C                  (feeded as REV to RSP* routines)
C  DDELT - required precision
C     
C  ALPHA and BETA - Euler angles (in degrees) specifying the 
C          orientation of the scattering particle relative to 
C          the laboratory reference frame (Refs. 6 and 7).
C
C  For axially symmetric scatterers, when the T matrix is computed in 
C  natural coordinate system with the $z$ axis along the axis of particle
C  axial symmetry, one can show that the T matrix is {\em diagonal} with 
C  respect to the azimuthal indices $m$ and $m'$ \cite{Wat},
C
C              T_{lm,l'm'}^{ij}=\delta_{mm'} T_{lm,l'm},
C
C  and that it satisfies reciprocity relation \cite{GuS,Mis36},
C
C               T_{lm,l'm}^{ij}=(-1)^{i+j} T_{l'm,lm}^{ji}.
C
C  \cite{Mis91} also claims the relation:
C
C                T_{lm,l'm}^{ij}= (-1)^{i+j} T_{l-m,l'-m}^{ij}
C
C--------/---------/---------/---------/---------/---------/---------/--
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER LMAXD,LMAX1D,LMTD
      INTEGER NAXSM,ICHOICEV,ICHOICE

      PARAMETER (LMAXD=8,LMAX1D=LMAXD+1,LMTD=LMAX1D*LMAX1D-1)
      INCLUDE 'ampld.par.f'
* number of the output unit
*
      REAL*8  LAM,MRR,MRI,X(NPNG2),W(NPNG2),S(NPNG2),SS(NPNG2),
     *        AN(NPN1),R(NPNG2),DR(NPNG2),
     *        DDR(NPNG2),DRR(NPNG2),DRI(NPNG2),ANN(NPN1,NPN1)
      REAL*8 TR1(NPN2,NPN2),TI1(NPN2,NPN2)            
c      REAL*8 XALPHA(300),XBETA(300),WALPHA(300),WBETA(300)

      COMPLEX*16 CZERO
      COMPLEX*16 zeps1,zeps0
      COMPLEX*16 TMT(4,LMTD,LMTD)
* 
      COMMON /CT/ TR1,TI1
* transfers the real and imaginary part of the T matrix (2*NMAX,2*NMAX) 
* array for a given value of M from TT via (TMATR0 and TMATR) to the main 
*
cc      COMMON /TMAT/ RT11,RT12,RT21,RT22,IT11,IT12,IT21,IT22
* transfers T matrix arrays obtained from TR1,TI1 in the main to the 
* AMPL routine --->  NOT USED HERE
*

      COMMON /CHOICE/ ICHOICE
* transfers the choice of inversion from here to TT
*
      COMMON/REVVAL/ A
* 
* transfers REV to CONST routine
*
      COMMON /TOITMT/ICHOICEV,NP,NCHECK,NAXSM,NDGS   
*
* transfers integers ICHOICEV,NP,NCHECK,NAXSM,NDGS from the main here
*
      COMMON /TOTMT/EPS,RAT,REV,ALPHA,BETA,DDELT   
* 
* transfers real*8 RAT,A(REV),ALPHA,BETA,DDELT from the main here
*     
*****************************************************************
      DATA CZERO/(0.D0,0.D0)/  
*
      P=DACOS(-1D0)                 !local PI constant
*
      ICHOICE=ICHOICEV
      A=REV
      LAM=REV*SQRT(ZEPS0)/RAP       !vacuum wavelength times SQRT(ZEPS0)/
*
* the real part of the refractive index contrast 
*      
      MRR=DBLE(SQRT(ZEPS1/ZEPS0))
*
* the imaginary  part of the refractive index contrast 
*
      MRI=DIMAG(SQRT(ZEPS1/ZEPS0)) 
* 
      DDELT=0.1D0*DDELT               !conv. test is switched off now!!!
*
* DDELT is used to test the accuracy of computing the 
* optical cross sections. This accuracy is usually better  
* than the absolute accuracy of computing the expansion coefficients 
* of a normalized scattering matrix by a factor of 10. Therefore,
* the desired accuracy of computing the expansion coefficients 
* is rescaled by a factor 0.1 before entering the test of the 
* accuracy of computing the optical cross sections.


*
* Other local constants:
*  
      LMTOT=(NMAX+1)**2-1

      IF (DABS(RAT-1D0).GT.1D-8.AND.NP.EQ.-1) CALL SAREA (EPS,RAT)
      IF (DABS(RAT-1D0).GT.1D-8.AND.NP.GE.0) CALL SURFCH(NP,EPS,RAT)
      IF (DABS(RAT-1D0).GT.1D-8.AND.NP.EQ.-2) CALL SAREAC (EPS,RAT)
      IF (NP.EQ.-3) CALL DROP (RAT)

*___________________________________________________
* Determination of the Wiscombe value of the floating 
C angular momentum cutoff NMAX:

      XEV=2D0*P*A/LAM
      IXXX=XEV+4.05D0*XEV**0.333333D0     !Wiscombe conv. criterion for NMAX
      INM1=MAX0(4,IXXX)
*
      IF (INM1.GE.NPN1) PRINT 7333, NPN1
      IF (INM1.GE.NPN1) STOP
 7333 FORMAT('CONVERGENCE IS NOT OBTAINED FOR NPN1=',I3,  
     &       '.  EXECUTION TERMINATED')

*_______________________________________________________________ 

      NGAUSS=NMAX*NDGS 
cc      NNNGGG=NGAUSS+1

      IF (NGAUSS.EQ.NPNG1) PRINT 7336
 7336    FORMAT('WARNING: NGAUSS=NPNG1')
*
* GIF division points and weights + other numerical constants
*
         CALL CONST(NGAUSS,NMAX,X,W,AN,ANN,S,SS,NP,EPS)
*
* specify particle shape:
*
         CALL VARY(LAM,MRR,MRI,A,EPS,NP,NGAUSS,X,P,PPI,PIR,PII,R,
     &              DR,DDR,DRR,DRI,NMAX)
*
* determine m=m'=0 elements of the T matrix
*
         CALL TMATR0 (NGAUSS,X,W,AN,ANN,PPI,PIR,PII,R,DR,
     &                 DDR,DRR,DRI,NMAX,NCHECK,NAXSM)
*
         QEXT=0D0
         QSCA=0D0
         
         DO 104 N=1,NMAX
            N1=N+NMAX
            TR1NN=TR1(N,N)
            TI1NN=TI1(N,N)
            TR1NN1=TR1(N1,N1)
            TI1NN1=TI1(N1,N1)
            
            DN1=DFLOAT(2*N+1) 
            
            QSCA=QSCA+DN1*(TR1NN*TR1NN+TI1NN*TI1NN
     &                    +TR1NN1*TR1NN1+TI1NN1*TI1NN1)
            QEXT=QEXT+(TR1NN+TR1NN1)*DN1
  104    CONTINUE

*<<<  
cc      WRITE(NOUT,*)'NMAX=',NMAX
cc      WRITE(NOUT,*)'NGAUSS=',NGAUSS
*<<<
*
* TMT initialization:

cc         DO JA=1,LMTOT     
cc         DO JB=1,LMTOT
*                                         
cc            TMT(1,JB,JA)=CZERO
cc            TMT(2,JB,JA)=CZERO 
cc            TMT(3,JB,JA)=CZERO
cc            TMT(4,JB,JA)=CZERO                       

cc         ENDDO
cc         ENDDO
*
C
C                         |  TMT(1,*) |  TMT(4,*)   |
C                 TMT  =  | ----------+-------------|
C                         |  TMT(3,*) |  TMT(2,*)   |
C
C    TMT(1,*) terms corresponds to TEE scattering sub-matrix
C    TMT(4,*) terms corresponds to TEM scattering sub-matrix
C    TMT(2,*) terms corresponds to TMM scattering sub-matrix
C    TMT(3,*) terms corresponds to TME scattering sub-matrix
C    TMT(4,*)=-TMT(3,*)^t where t denotes transposed TMT(3,*) submatrix

*****************     Assign  m=m=0 elements of TMT matrix   ***********

         DO L1=1,NMAX   
         DO L2=1,NMAX  
          N1=L1+NMAX
          N2=L2+NMAX
 
            JA=L1*(L1+1)       ! (l,m) index with (1-1)=1
            JB=L2*(L2+1)       
* 
* see (5.39) of {MTL}:  !!!Iff plane of symmetry perpendicular to the
                        !!!axis of rotation

	    if ((NAXSM.eq.1).and.((-1)**(L1+L2).ne.1)) then
	      TMT(2,JA,JB)=CZERO 
            TMT(1,JA,JB)=CZERO  
		else                                         
            TMT(2,JA,JB)=DCMPLX(TR1(L1,L2),TI1(L1,L2))
            TMT(1,JA,JB)=DCMPLX(TR1(N1,N2),TI1(N1,N2))
	    end if
* see (5.37) of {MTL}:
            TMT(4,JA,JB)=CZERO         !DCMPLX(TR1(N1,L2),TI1(N1,L2))
            TMT(3,JB,JA)=CZERO         !-TMT(4,JA,JB)  
                      
         ENDDO
         ENDDO
*
*****************    Assign  m=m'>0 elements of the T matrix   ***********

      DO 220 M=1,NMAX
*
         CALL TMATR(M,NGAUSS,X,W,AN,ANN,S,SS,PPI,PIR,PII,R,DR,
     &               DDR,DRR,DRI,NMAX,NCHECK,NAXSM)
*
* <<< returns  m=m'>0 elements of the T matrix
*
        NM=NMAX-M+1             !size of a T block returned by TT

         DO L1=M,NMAX   
         DO L2=M,NMAX 
 
          K1=L1-M+1               !K1,K2,KK1,KK2 label the entries of
	    K2=L2-M+1               !a TT returned T block
          KK1=L1-M+1+NM
          KK2=L2-M+1+NM

            JAM=L1*(L1+1)-M       ! (l,m) index with (1-1)=1
            JBM=L2*(L2+1)-M  
            JA =L1*(L1+1)+M       ! (l,m) index with (1-1)=1
            JB =L2*(L2+1)+M        
*
* 
* see (5.39) of {MTL}: !!!Iff plane of symmetry perpendicular to the
                       !!!axis of rotation

	    if ((NAXSM.eq.1).and.((-1)**(L1+L2).ne.1)) then
	      TMT(2,JA,JB)  =CZERO 
	      TMT(2,JAM,JBM)=CZERO
            TMT(1,JA,JB)  =CZERO  
	      TMT(1,JAM,JBM)=CZERO  
		else                                         
            TMT(2,JA,JB)   = DCMPLX(TR1(K1,K2),TI1(K1,K2))
            TMT(2,JAM,JBM) = TMT(2,JA,JB)            
            TMT(1,JA,JB)   = DCMPLX(TR1(KK1,KK2),TI1(KK1,KK2))
            TMT(1,JAM,JBM) = TMT(1,JA,JB)
	    end if

	    if ((NAXSM.eq.1).and.((-1)**(L1+L2).ne.-1)) then
	      TMT(4,JA,JB)  =CZERO 
	      TMT(4,JAM,JBM)=CZERO
            TMT(3,JA,JB)  =CZERO  
	      TMT(3,JAM,JBM)=CZERO  
		else  
            TMT(4,JA,JB)   = DCMPLX(TR1(KK1,K2),TI1(KK1,K2))
            TMT(4,JAM,JBM) =-TMT(4,JA,JB)
            TMT(3,JB,JA)   =-TMT(4,JA,JB) 
            TMT(3,JBM,JAM) = TMT(4,JA,JB)    !=-TMT(3,JB,JA) 
	    end if
*    
*  Using reciprocity (Eq. (15) of Ref. \ct{Mis97}):
*
*       T_{lm,l'm}^{ij}=(-1)^{i+j} T_{l'm,lm}^{ji}
*
*  and (see Eq. (36) \JQSRT{55}):
*
*       T_{lm,l'm'}^{ij}=(-1)^{m+m'} T_{l'-m',l-m}^{ji}
*                      
*  Moreover, for axially symmetric particles one has 
*  (see Eq. (31),(36) \JQSRT{55}):
*
*  T_{lm,l'm}^{ij}=(-1)^{i+j} T_{l-m,l'-m}^{ij} = T_{l'-m',l-m}^{ji}
*
         ENDDO
         ENDDO
*
  220 CONTINUE    !end of loop over m's
  
      RETURN
      END
 
C********************************************************************
                                                
      SUBROUTINE AMPL (NMAX,DLAM,TL,TL1,PL,PL1,ALPHA,BETA,
     &                 VV,VH,HV,HH)  
C--------/---------/---------/---------/---------/---------/---------/--
C >>> NMAX,DLAM,TL,TL1,PL,PL1,ALPHA,BETA
C <<< RAT
C=================
C
C    GIVEN T MATRIX IN COMMON BLOCK IT CALCULATES THE AMPLITUDE MATRIX
C
C   NMAX - angular momentum cutoff
C   DLAM - wavelength of incident light
C   TL (THET0 IN MAIN) - zenith angle of the incident beam in degrees
C   TL1 (THET IN MAIN) - zenith angle of the scattered beam in degrees    
C   PL (PHI0 IN MAIN) - azimuth angle of the incident beam in degrees    
C   PL1 (PHI IN MAIN) - azimuth angle of the scattered beam in degrees 
C   ALPHA and BETA - Euler angles (in degrees) specifying the 
C         orientation of the scattering particle relative to the 
C         laboratory reference frame (Refs. 6 and 7).  
C 
C--------/---------/---------/---------/---------/---------/---------/--
      IMPLICIT REAL*8 (A-B,D-H,O-Z)
      IMPLICIT COMPLEX*16 (C)
      INTEGER NOUT
      
* number of the output unit
      PARAMETER (NOUT=35)
      INCLUDE 'ampld.par.f'
      
      REAL*8 AL(3,2),AL1(3,2),AP(2,3),AP1(2,3),B(3,3),
     *       R(2,2),R1(2,2),C(3,2),CA,CB,CT,CP,CTP,CPP,CT1,CP1,
     *       CTP1,CPP1
      REAL*8 DV1(NPN6),DV2(NPN6),DV01(NPN6),DV02(NPN6)
      REAL*4
     &     TR11(NPN6,NPN4,NPN4),TR12(NPN6,NPN4,NPN4),
     &     TR21(NPN6,NPN4,NPN4),TR22(NPN6,NPN4,NPN4),
     &     TI11(NPN6,NPN4,NPN4),TI12(NPN6,NPN4,NPN4),
     &     TI21(NPN6,NPN4,NPN4),TI22(NPN6,NPN4,NPN4)
      COMPLEX*16 CAL(NPN4,NPN4),VV,VH,HV,HH
*_____
      COMMON /TMAT/ TR11,TR12,TR21,TR22,TI11,TI12,TI21,TI22
*_____

      IF (ALPHA.LT.0D0.OR.ALPHA.GT.360D0.OR.
     &    BETA.LT.0D0.OR.BETA.GT.180D0.OR.
     &    TL.LT.0D0.OR.TL.GT.180D0.OR.
     &    TL1.LT.0D0.OR.TL1.GT.180D0.OR.
     &    PL.LT.0D0.OR.PL.GT.360D0.OR.
     &    PL1.LT.0D0.OR.PL1.GT.360D0) THEN 
          WRITE(NOUT,2000)
          STOP
      ELSE
          CONTINUE
      ENDIF  
 2000 FORMAT ('AN ANGULAR PARAMETER IS OUTSIDE ITS',
     &        ' ALLOWABLE RANGE')

* SPECIFYING NUMERICAL CONSTANTS:

      PIN=DACOS(-1D0)         !=PI
      PIN2=PIN*0.5D0          !=PI/2
      PI=PIN/180D0            !=PI/180
      
* conversion from degrees to radians:
      ALPH=ALPHA*PI
      BET=BETA*PI
      THETL=TL*PI
      PHIL=PL*PI
      THETL1=TL1*PI
      PHIL1=PL1*PI

      EPS=1D-7
      IF (THETL.LT.PIN2) THETL=THETL+EPS
      IF (THETL.GT.PIN2) THETL=THETL-EPS
      IF (THETL1.LT.PIN2) THETL1=THETL1+EPS
      IF (THETL1.GT.PIN2) THETL1=THETL1-EPS
      IF (PHIL.LT.PIN) PHIL=PHIL+EPS
      IF (PHIL.GT.PIN) PHIL=PHIL-EPS
      IF (PHIL1.LT.PIN) PHIL1=PHIL1+EPS
      IF (PHIL1.GT.PIN) PHIL1=PHIL1-EPS
      IF (BET.LE.PIN2.AND.PIN2-BET.LE.EPS) BET=BET-EPS
      IF (BET.GT.PIN2.AND.BET-PIN2.LE.EPS) BET=BET+EPS
      
C_____________COMPUTE THETP, PHIP, THETP1, AND PHIP1, EQS. (8), (19), AND (20)

      CB=DCOS(BET)
      SB=DSIN(BET)
      CT=DCOS(THETL)
      ST=DSIN(THETL)
      CP=DCOS(PHIL-ALPH)
      SP=DSIN(PHIL-ALPH)
      
      CTP=CT*CB+ST*SB*CP
      THETP=DACOS(CTP)
      CPP=CB*ST*CP-SB*CT
      SPP=ST*SP
      PHIP=DATAN(SPP/CPP)
      
      IF (PHIP.GT.0D0.AND.SP.LT.0D0) PHIP=PHIP+PIN
      IF (PHIP.LT.0D0.AND.SP.GT.0D0) PHIP=PHIP+PIN
      IF (PHIP.LT.0D0) PHIP=PHIP+2D0*PIN

      CT1=DCOS(THETL1)
      ST1=DSIN(THETL1)
      CP1=DCOS(PHIL1-ALPH)
      SP1=DSIN(PHIL1-ALPH)
      
      CTP1=CT1*CB+ST1*SB*CP1
      THETP1=DACOS(CTP1)
      CPP1=CB*ST1*CP1-SB*CT1
      SPP1=ST1*SP1
      PHIP1=DATAN(SPP1/CPP1)
      
      IF (PHIP1.GT.0D0.AND.SP1.LT.0D0) PHIP1=PHIP1+PIN
      IF (PHIP1.LT.0D0.AND.SP1.GT.0D0) PHIP1=PHIP1+PIN
      IF (PHIP1.LT.0D0) PHIP1=PHIP1+2D0*PIN

C____________COMPUTE MATRIX BETA, EQ. (22) of {Mis39}

      CA=DCOS(ALPH)
      SA=DSIN(ALPH)
      B(1,1)=CA*CB
      B(1,2)=SA*CB
      B(1,3)=-SB
      B(2,1)=-SA
      B(2,2)=CA
      B(2,3)=0D0
      B(3,1)=CA*SB
      B(3,2)=SA*SB
      B(3,3)=CB

C____________COMPUTE MATRICES AL AND AL1, EQ. (14)  of {Mis39}

      CP=DCOS(PHIL)
      SP=DSIN(PHIL)
      CP1=DCOS(PHIL1)
      SP1=DSIN(PHIL1)

* Eq. (15) of {Mis39}:     
      AL(1,1)=CT*CP
      AL(1,2)=-SP
      AL(2,1)=CT*SP
      AL(2,2)=CP
      AL(3,1)=-ST
      AL(3,2)=0D0
      
* Eq. (16) of {Mis39}:       
      AL1(1,1)=CT1*CP1
      AL1(1,2)=-SP1
      AL1(2,1)=CT1*SP1
      AL1(2,2)=CP1
      AL1(3,1)=-ST1
      AL1(3,2)=0D0

C____________COMPUTE MATRICES AP^(-1) AND AP1^(-1), EQ. (15) 

      CT=CTP
      ST=DSIN(THETP) 
      CP=DCOS(PHIP)
      SP=DSIN(PHIP)
      CT1=CTP1
      ST1=DSIN(THETP1)
      CP1=DCOS(PHIP1)
      SP1=DSIN(PHIP1)
      AP(1,1)=CT*CP
      AP(1,2)=CT*SP
      AP(1,3)=-ST  
      AP(2,1)=-SP
      AP(2,2)=CP 
      AP(2,3)=0D0
      AP1(1,1)=CT1*CP1
      AP1(1,2)=CT1*SP1
      AP1(1,3)=-ST1   
      AP1(2,1)=-SP1
      AP1(2,2)=CP1 
      AP1(2,3)=0D0

C____________COMPUTE MATRICES R AND R^(-1), EQ. (13)
      DO I=1,3
         DO J=1,2
            X=0D0
            DO K=1,3
               X=X+B(I,K)*AL(K,J)
            ENDDO
            C(I,J)=X
         ENDDO
      ENDDO
      DO I=1,2
         DO J=1,2
            X=0D0
            DO K=1,3
               X=X+AP(I,K)*C(K,J)
            ENDDO
            R(I,J)=X
         ENDDO
      ENDDO
      DO I=1,3
         DO J=1,2
            X=0D0
            DO K=1,3
               X=X+B(I,K)*AL1(K,J)
            ENDDO
            C(I,J)=X
         ENDDO
      ENDDO
      DO I=1,2
         DO J=1,2
            X=0D0
            DO K=1,3
               X=X+AP1(I,K)*C(K,J)
            ENDDO
            R1(I,J)=X
         ENDDO
      ENDDO
      D=1D0/(R1(1,1)*R1(2,2)-R1(1,2)*R1(2,1))
      X=R1(1,1)
      R1(1,1)=R1(2,2)*D
      R1(1,2)=-R1(1,2)*D
      R1(2,1)=-R1(2,1)*D
      R1(2,2)=X*D

      CI=(0D0,1D0)
      DO 5 NN=1,NMAX
         DO 5 N=1,NMAX
            CN=CI**(NN-N-1)
            DNN=DFLOAT((2*N+1)*(2*NN+1)) 
            DNN=DNN/DFLOAT( N*NN*(N+1)*(NN+1) ) 
            RN=DSQRT(DNN)
            CAL(N,NN)=CN*RN
    5 CONTINUE
      DCTH0=CTP
      DCTH=CTP1 
      PH=PHIP1-PHIP
      VV=(0D0,0D0)
      VH=(0D0,0D0)
      HV=(0D0,0D0)
      HH=(0D0,0D0)
      DO 500 M=0,NMAX
         M1=M+1
         NMIN=MAX(M,1)
*
* Specify Wigner d-matrices:

         CALL VIGAMPL (DCTH, NMAX, M, DV1, DV2)
         CALL VIGAMPL (DCTH0, NMAX, M, DV01, DV02)
*
         FC=2D0*DCOS(M*PH)
         FS=2D0*DSIN(M*PH)
         DO 400 NN=NMIN,NMAX
            DV1NN=DV01(NN)         !changed from  M*DV01(NN)
            DV2NN=DV02(NN)
            DO 400 N=NMIN,NMAX
               DV1N=DV1(N)         !changed from  M*DV1(NN)
               DV2N=DV2(N)

               CT11=DCMPLX(TR11(M1,N,NN),TI11(M1,N,NN))
               CT22=DCMPLX(TR22(M1,N,NN),TI22(M1,N,NN))

               IF (M.EQ.0) THEN

                  CN=CAL(N,NN)*DV2N*DV2NN

                  VV=VV+CN*CT22  
                  HH=HH+CN*CT11

                 ELSE

                  CT12=DCMPLX(TR12(M1,N,NN),TI12(M1,N,NN))
                  CT21=DCMPLX(TR21(M1,N,NN),TI21(M1,N,NN))

                  CN1=CAL(N,NN)*FC
                  CN2=CAL(N,NN)*FS

                  D11=DV1N*DV1NN
                  D12=DV1N*DV2NN
                  D21=DV2N*DV1NN
                  D22=DV2N*DV2NN

                  VV=VV+(CT11*D11+CT21*D21   
     &                  +CT12*D12+CT22*D22)*CN1   

                  VH=VH+(CT11*D12+CT21*D22   
     &                  +CT12*D11+CT22*D21)*CN2

                  HV=HV-(CT11*D21+CT21*D11
     &                  +CT12*D22+CT22*D12)*CN2   

                  HH=HH+(CT11*D22+CT21*D12
     &                  +CT12*D21+CT22*D11)*CN1      
               ENDIF
  400    CONTINUE
  500 CONTINUE
      DK=2D0*PIN/DLAM
      VV=VV/DK
      VH=VH/DK
      HV=HV/DK
      HH=HH/DK
      CVV=VV*R(1,1)+VH*R(2,1)
      CVH=VV*R(1,2)+VH*R(2,2)
      CHV=HV*R(1,1)+HH*R(2,1)
      CHH=HV*R(1,2)+HH*R(2,2)
      VV=R1(1,1)*CVV+R1(1,2)*CHV
      VH=R1(1,1)*CVH+R1(1,2)*CHH
      HV=R1(2,1)*CVV+R1(2,2)*CHV
      HH=R1(2,1)*CVH+R1(2,2)*CHH
      
      
      PRINT 1101, VV
      PRINT 1102, VH
      PRINT 1103, HV
      PRINT 1104, HH
 1101 FORMAT ('S11=',D11.5,' + i*',D11.5)
 1102 FORMAT ('S12=',D11.5,' + i*',D11.5)
 1103 FORMAT ('S21=',D11.5,' + i*',D11.5)
 1104 FORMAT ('S22=',D11.5,' + i*',D11.5)

      RETURN
      END
     
C*********************************************************************

      SUBROUTINE VIGAMPL (X, NMAX, M, DDV1, DV2)
C--------/---------/---------/---------/---------/---------/---------/--
C >>> X,NMAX,M (only nonnegative)
C <<< DV1, DV2
C =============
C     For a given azimuthal number M.GE.0 returns 
C      the Wigner d-functions divided by sin\theta, i.e.,
C
C     DDV1(N)=dvig(0,m,n,arccos x)/sin(arccos x)   ! = m*d_{0m}^{(l)}/ sin\theta
C     and
C     DV2(N)=[d/d(arccos x)] dvig(0,m,n,arccos x)  ! = d d_{0m}^{(l)}/d\theta
C
C     for 1.LE.N.LE.NMAX and 0.LE.X.LE.1
C     (For a given M.NEQ.0, only the M.LE.N.LE.NMAX terms are determined!)
C     According to Eq. (4.1.24) of Ref. \ct{Ed}:
C
C             d_{00}^{(l)}(\theta)= P_l(\cos\theta)
C
C     (Rodrigues formula [Eq. (2.5.14) of Ref. \ct{Ed}] then yields 
C                       P_1(x)=x; P_2=(3x^2-1)/2; etc.
C     One can show that $d_{00}^{(1)}(\theta)=\cos\theta$
C
C     Similar to routine VIG, which however returns 
C     DV1(N)=dvig(0,m,n,arccos x)   ! = d_{0m}^{(l)}
C
C     In addition, VIGAMPL has a block treating the case when 
C     arccos x is very small option
C
C     Made using recurrences of  Ref. \ct{Mis39}  
C     (There is a missing $l$ factor in the 2nd term in the curly bracket 
C     in recurrence (35) of Ref. \ct{Mis39} for DV2).         
C
C     One has (see Eq. (4.2.5) of \ct{Ed}):
C                       $d_{0m}^{(l)}=(-1)^m d_{0-m}^{(l)}$
C     and (see Eq. (35) of \ct{Mis91}):
C            $dd_{0m}^{(l)}/(d\theta)=(-1)^m dd_{0-m}^{(l)}/(d\theta)$
C
C     X=cos(theta), where theta is the polar angle
C     LMAXD ... maximal angular momentum cutoff
C     NMAX ... floating  angular momentum cutoff
C     CALLED ONLY BY the AMPL routine!!!
C--------/---------/---------/---------/---------/---------/---------/--
      INCLUDE 'ampld.par.f'
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 DDV1(NPN1),DV2(NPN1)

c      IMPLICIT none
c      INTEGER LMAXD,LMAXD1     
c      PARAMETER (LMAXD=50,LMAXD1=LMAXD+1)
c      REAL*8 DDV1(LMAXD1),DV1(LMAXD1),DV2(LMAXD1)

      integer n,nmax,M,I,I2
      REAL*8 A,X,QS,QS1,DSI,D1,D2,D3,DER,DN,DX,QN,QN1,QN2,
     & QNM,QNM1,QMM

* DDV1 and DV2 initialization
      DO 1 N=1,NMAX
 	   DDV1(N)=0.D0
         DV2(N) =0.D0
    1 CONTINUE

      DX=DABS(X)
      A=1.D0
      QS=DSQRT(1D0-X*X)                        !sin\theta
*
      IF (M.NE.0) then          ! M\neq 0 part
*                          
*    A_m*(sin\theta)**m   initialization - (33) and recurrence (34) of Ref. {Mis39}

   20 QMM=DFLOAT(M*M)

      DO 25 I=1,M
         I2=I*2
         A=A*DSQRT(DFLOAT(I2-1)/DFLOAT(I2))*QS  !recurrence (33,34) of Ref. {Mis39} f
   25 CONTINUE

      end if
*
      go to 100

*********************************************************     
*  (1-cos\theta) is very small:
*
C   For theta=0 [see Eqs. above]:
C              d_{00}^{(0)}(0)=1 
C              d_{01}^{(1)}(0)=0
C              d_{02}^{(2)}(\beta)=0
C     and
C         d d_{00}^{(0)}(\beta)/d\beta=0
C         d d_{01}^{(1)}(\beta)/d\beta=1/\sqrt{2}  
C         d d_{02}^{(2)}(\beta)/d\beta=0
C
C  See Eqs. (4.1-4) of \ct{Mis91}:
C
C   (m/\sin\theta) d_{0m}^l(0)=(\delta_{m\pm 1}/2) \sqrt{l(l+1)}
C      d d_{0m}^l(0)/d\beta   =(m\delta_{m\pm 1}/2) \sqrt{l(l+1)}
C
*
*  (4.2.1) of \ct{Ed}:
*   d_{0m}^{(l)}(pi) = (-1)^{l+m} \dt_{0,m}
*  
*  (4.2.3) of \ct{Ed}:
*   d_{0m}^{(l)}(0) = (-1)^{m} \dt_{0,m} = \dt_{0,m}
*=======================================
*
*  If X^l_m=(m/\sin\theta) d_{0m}^{(l)}, then, according to (3.29) of {TKS}:
*
*  X^{m+1}_{m+1}=\sin\theta \sqrt{\fr{2m+1}{2m+2}}
*                           \left(\fr{m+1}{m}\right)X^{m}_{m}
*
*  According to (3.30) of {TKS}:
*  X^{m+1}_{m}= -\sqrt{2m+1}\,\cos\theta X^{m}_{m}
*
* According to (3.31) of {TKS}:
*  X^{l}_{m}=\fr{1}{\sqrt{l^2-m^2}}\,\left[(2l-1)\cos\theta
*          X^{l-1}_{m} - \sqrt{(l-1)^2-m^2}}\,\X^{l-2}_{m} \right]
* 
* Initial recurrence values are X^1_1=\sqrt{2}/2 and X^l_0=0
***********************************************************************
*                   NONZERO DDV1/DV2 INITIALIZATION
*                          M = 0

 100  IF (M.EQ.0) THEN     !all DDV1(N)=X^l_0=0; see (3.33) of {TKS}: 

* According to (3.37) of {TKS}, DV2(0)=0.d0

      DV2(1)=QS

      IF (NMAX.GE.2) DV2(2)=3*X*DV2(1)

      IF (NMAX.LT.3) RETURN
*
      DO N=3,NMAX           !recurrence (3.36) of {TKS}, 
      DV2(N)=(2*N-1)*X*DV2(N-1)/(N-1)-N*DV2(N-2)/(N-1)
	ENDDO
***********************************************************************
*                           M > 0

 	ELSE IF (M.GT.0) THEN       
*
* >>> Determine X^m_m according to Eq. (3.29) of {TKS}:

	A=1.d0/DSQRT(2.D0)               !X^1_1=A_1
	
	DO I=1,M-1
      A=QS*DBLE(I+1)*DSQRT(2*I+1.d0)*A/(I*DSQRT(2*I+2.d0))
	ENDDO              

* <<< A is now X^m_m; see (3.29) of {TKS}

      DDV1(M)=A
      DV2(M)=X*A                        !see (3.34) of {TKS}

* >>> Determine X^{m+1}_m:

	IF (M.EQ.NMAX)  GO TO 120        

      DER=X*DSQRT(2*M+1.d0)*A          ! DER=X^{m+1}_m; see (3.30) of {TKS}
      DDV1(M+1)=DER
      DV2(M+1)=((M+1)*X*DER-A*DSQRT(2*M+1.d0))/DBLE(M)  !(3.35) of {TKS}

* >>> Determine remaining X^{l}_m's

	IF ((M+2).EQ.NMAX)  GO TO 120 

       DO N=M+2,NMAX
       D3=DSQRT(DBLE(N)**2-DBLE(M)**2)
       DDV1(N)=((2*N-1)*X*DDV1(N-1) - 
     &                DSQRT(DBLE(N-1)**2-DBLE(M)**2)*DDV1(N-2))/D3      
                                                      !see (3.31) of {TKS}
	 DV2(N)=(N*X*DDV1(N)-DDV1(N-1)*D3)/DBLE(M)      !see (3.35) of {TKS}
       ENDDO

	END IF

  120 RETURN
      END 
C**********************************************************************

      SUBROUTINE CONST (NGAUSS,NMAX,X,W,AN,ANN,S,SS,NP,EPS)
C--------/---------/---------/---------/---------/---------/---------/--
C >>> NGAUSS,NMAX,NP,EPS
C <<< X,W,AN,ANN,S,SS
C=====================
C
C  NGAUSS - the number of GIF division points
C  NMAX - angular momentum cutoff
C  P=PI=DACOS(-1d0)
C  NP - parameter specifying the particle shape 
C  EPS - deformation parameter for a given particle shape
C
C  X=\cos\theta  - GIF division points
C  W - GIF weights
C  AN(N)=N*(N+1)
C  ANN(l_1,l_2)=\sqrt{\fr{(2 l_1+1)}{l_1(l_1+1)} }
C                       \sqrt{\fr{(2 l_2+1)}{l_2(l_2+1)} }/2
C  S  ... 1/(|\sin\theta|)
C  SS ... 1/(\sin^2\theta)
C--------/---------/---------/---------/---------/---------/---------/--
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'ampld.par.f'
      INTEGER NEPS,JG,I,J
	REAL*8 EE,EE1,CC,SI,XI1,XI2,XAV
      REAL*8 XTHETA,THETA0,RX,PI
      REAL*8 X(NPNG2),W(NPNG2),X1(NPNG2),W1(NPNG2),
     *        X2(NPNG2),W2(NPNG2),
     *        S(NPNG2),SS(NPNG2),
     *        AN(NPN1),ANN(NPN1,NPN1),DD(NPN1)
*
      COMMON/REVVAL/ RX
*
      DATA PI/3.141592653589793d0/ 
*      
      DO 10 N=1,NMAX
           NN=N*(N+1)
           AN(N)=DFLOAT(NN)
           D=DSQRT(DFLOAT(2*N+1)/DFLOAT(NN))
           DD(N)=D
           DO 10 N1=1,N
                DDD=D*DD(N1)*0.5D0
                ANN(N,N1)=DDD
                ANN(N1,N)=DDD
   10 CONTINUE
   
      NG=2*NGAUSS      
*
* GIF division points and weights
* 
	NEPS=MAX(EPS,1.d0/EPS)      !number of Gauss integration  
	                            !intervals from EPS

      IF (NP.EQ.-1) THEN         ! spheroid

	IF(NEPS.EQ.1) THEN
     
      CALL GAUSS(NG,0,0,X,W)

	ELSE IF (NEPS.EQ.2) THEN

	CALL GAULEG(-1.d0,0.d0,X,W,NGAUSS)

      DO I=1,NGAUSS

         W(I+NGAUSS)=W(NGAUSS-I+1)
         X(I+NGAUSS)=-X(NGAUSS-I+1)

      ENDDO

	ELSE IF (NEPS.GT.2) THEN
	
	NEPS=3*NEPS
      NG1=DFLOAT(NGAUSS)/NEPS
	
      EE=EPS*EPS
      EE1=EE-1D0

	XAV=0.d0

	DO I=1,NEPS

	XI1=DBLE(I)/(NEPS+1)
         
          CC=XI1*XI1
          SI=1D0-CC
          X2(I)=ABS(XI1*SI*EE1/(SI+EE*CC))       !|dr(theta)/dtheta|
	    XAV=XAV+1.d0/X2(I)

	ENDDO

	XAV=XAV           !averaged 1/|dr(theta)/dtheta|

c_____ estimate integration intervals:
   
 	DO I=1,NEPS
	
	X2(I)=1.d0/(XAV*X2(I))

	ENDDO
	      
	DO I=1,NEPS

	IF(i.eq.1) then

	XI1=0.d0
	XI2=X2(1)

	else 

	XI2=XI2+X2(I)

	end if
	 
	JG=NGAUSS+(I-1)*NG1

	IF(I.EQ.NEPS) NG1=NGAUSS-(I-1)*NG1
	
	CALL GAULEG(XI1,XI2,X1,W1,NG1)

	XI1=XI2

      DO  J=1,NG1
         W(JG+J)=W1(J)
         X(JG+J)=X1(J)
	ENDDO         !J

	ENDDO         !I

*
* Assuming mirror symmetry in the $\theta=\pi/2$ plane
*
      DO  I=1,NGAUSS
         W(I)=W(NG-I+1)
         X(I)=-X(NG-I+1)
	ENDDO

	ENDIF           !NEPS     

      ELSE IF (NP.EQ.-2) THEN         ! cylinder

******************   Only involves cylinders  ********************** 
     
      NG1=DFLOAT(NGAUSS)/2D0
      NG2=NGAUSS-NG1
      XX=-DCOS(DATAN(EPS))        !-COS OF SEPARATION ANGLE BETWEEN
                                  !HORIZONTAL AND VERTICAL CYLINDER
                                  !FACES
*
* GIF division points and weights
*
      CALL GAUSS(NG1,0,0,X1,W1)         !for (0,NG1)
      CALL GAUSS(NG2,0,0,X2,W2)         !for (NG1+1,NGAUSS=NG1+NG2)
*
C In GAUSS (N,IND1,IND2,Z,W):  
C IND1 = 0 - INTERVAL (-1,1), 
C IND1 = 1 - (0,1)
C IND2 = 1 RESULTS ARE PRINTED.
*
*
      DO 12 I=1,NG1
         W(I)=0.5D0*(XX+1D0)*W1(I)
         X(I)=0.5D0*(XX+1D0)*X1(I)+0.5D0*(XX-1D0)
   12 CONTINUE
      DO 14 I=1,NG2
         W(I+NG1)=-0.5D0*XX*W2(I)
         X(I+NG1)=-0.5D0*XX*X2(I)+0.5D0*XX
   14 CONTINUE
*
* Assuming mirror symmetry in the $\theta=\pi/2$ plane
*
      DO 16 I=1,NGAUSS
         W(NG-I+1)=W(I)
         X(NG-I+1)=-X(I)
   16 CONTINUE

******************************************************************  
      ELSE IF (NP.EQ.-4) THEN         ! cut sphere on top
           
      XTHETA=DACOS((EPS-RX)/RX)
      XX=DSIN(XTHETA)/(PI-XTHETA)
      NG2=XX*DBLE(NG)
      NG1=NG-NG2
      THETA0=1.D0/SQRT(8.D0*RX/EPS-3.D0)  !cosine of the separation angle 
      XX=THETA0
*
      CALL GAULEG(-1.D0,THETA0,X1,W1,NG1)       !for (0,NG1)
      CALL GAULEG(THETA0,1.D0,X2,W2,NG2)        !for (NG2+1,NG=NG1+NG2)
*
      DO  I=1,NG1
         W(I)=W1(I)
         X(I)=X1(I)
      ENDDO
      
      DO I=1,NG2

         W(I+NG1)=W2(I)
         X(I+NG1)=X2(I)

      ENDDO
*
******************************************************************  
      ELSE IF (NP.EQ.-5) THEN            ! cut sphere on its bottom
           
      XTHETA=DACOS((EPS-RX)/RX)
      XX=DSIN(XTHETA)/(PI-XTHETA)
      NG1=XX*DBLE(NG)
      NG2=NG-NG1
      THETA0=-1.D0/SQRT(8.D0*RX/EPS-3.D0)  !cosine of the separation angle 
*
      CALL GAULEG(-1.D0,THETA0,X1,W1,NG1)       !for (0,NG1)
      CALL GAULEG(THETA0,1.D0,X2,W2,NG2)        !for (NG2+1,NG=NG1+NG2)
*
      DO  I=1,NG1
         W(I)=W1(I)
         X(I)=X1(I)
      ENDDO
      
      DO I=1,NG2

         W(I+NG1)=W2(I)
         X(I+NG1)=X2(I)

      ENDDO

****************************************************************** 
*   
      ELSE 
*      
      CALL GAUSS(NG,0,0,X,W)
c      CALL GAULEG(-1.D0,1.D0,X,W,NG)
*     
      END IF
*
       if (np.gt.-4) then           !mirror symmetry present

       DO 20 I=1,NGAUSS
           Y=X(I)
           Y=1D0/(1D0-Y*Y)
           SS(I)=Y
           SS(NG-I+1)=Y
           Y=DSQRT(Y)
           S(I)=Y
           S(NG-I+1)=Y
   20 CONTINUE

      else                         !mirror symmetry absent
   
       DO 30 I=1,NG
           Y=X(I)
           Y=1D0/(1D0-Y*Y)
           SS(I)=Y
           Y=DSQRT(Y)
           S(I)=Y 
   30 CONTINUE

      END IF
*
      RETURN
      END
 
C**********************************************************************
 
      SUBROUTINE VARY (LAM,MRR,MRI,A,EPS,NP,NGAUSS,X,P,PPI,PIR,PII,
     *                 R,DR,DDR,DRR,DRI,NMAX)
C--------/---------/---------/---------/---------/---------/---------/--
C >>> LAM,MRR,MRI,A,EPS,NP,NGAUSS,X,P,NMAX
C <<< PPI,PIR,PII,R,DR,DDR,DRR,DRI
C=========================
C  LAM - wavelength of incident light
C  MRR - the real part of the refractive index 
C  MRI - the imaginary  part of the refractive index
C  A=RAT*AXI, where RAT and AXI are the main program input parameters                            
    
C      RAT = 1 - particle size is specified in terms of the            
C                equal-volume-sphere radius                             
C      RAT.ne.1 - particle size is specified in terms of the           
C                equal-surface-area-sphere radius
C  AXI - equivalent-(volume/surface-area)-sphere radius  
C  NP - particle shape class 
C  EPS - shape deformation parameter within a given particle shape class
C  NGAUSS - the number of Gauss integration division points 
C           in the integral over theta 
C  NMAX - angular momentum cutoff 
C  P=DACOS(-1D0)
C  PI=P*2D0/LAM - wave vector
C  PPI=PI*PI
C  PIR=PPI*MRR
C  PII=PPI*MRI
C  R=r^2(\theta)            for axially symmetric particles
C  DR=dr(\theta)/(d\theta)  for axially symmetric particles
C  DDR=\lambda/[2*\pi*r(\theta)]
C  DRR=(MRR/(MRR**2+MRI**2))*(\lambda/[2*\pi*r(\theta)])
C  DRI=-(MRI/(MRR**2+MRI**2))*(\lambda/[2*\pi*r(\theta)])
C--------/---------/---------/---------/---------/---------/---------/--
      INCLUDE 'ampld.par.f'
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8  X(NPNG2),R(NPNG2),DR(NPNG2),MRR,MRI,LAM,
     *        Z(NPNG2),ZR(NPNG2),ZI(NPNG2),
     *        DDR(NPNG2),DRR(NPNG2),DRI(NPNG2)
cc     *        J(NPNG2,NPN1),Y(NPNG2,NPN1),JR(NPNG2,NPN1),
cc     *        JI(NPNG2,NPN1),DJ(NPNG2,NPN1),DY(NPNG2,NPN1),
cc     *        DJR(NPNG2,NPN1),DJI(NPNG2,NPN1)
cc      COMMON /CBESS/ J,Y,JR,JI,DJ,DY,DJR,DJI
      NG=NGAUSS*2

* decision tree to specify particle shape:

      IF (NP.GT.0)  CALL RSP2(X,NG,A,EPS,NP,R,DR)       ! Chebyshev particle
      IF (NP.EQ.-1) CALL RSP1(X,NG,NGAUSS,A,EPS,R,DR)   ! oblate/prolate spheroids 
      IF (NP.EQ.-2) CALL RSP3(X,NG,NGAUSS,A,EPS,R,DR)   ! oblate/prolate cylinder
      IF (NP.EQ.-3) CALL RSP4(X,NG,A,R,DR)              ! a distorted Chebyshev droplet
      IF (NP.EQ.-4) CALL RSP5(X,NG,A,EPS,R,DR)          ! sphere cut by a plane on its top   
      IF (NP.EQ.-5) CALL RSPI5(X,NG,A,EPS,R,DR)         ! sphere cut by a plane on its bottom          
*
      PI=P*2D0/LAM                 !wave vector
      PPI=PI*PI
      PIR=PPI*MRR
      PII=PPI*MRI
      V=1D0/(MRR*MRR+MRI*MRI)
      PRR=MRR*V
      PRI=-MRI*V
      TA=0D0
      DO 10 I=1,NG
           VV=DSQRT(R(I))
           V=VV*PI
           TA=MAX(TA,V)            !Max. size parameter
           VV=1D0/V
           DDR(I)=VV
           DRR(I)=PRR*VV
           DRI(I)=PRI*VV
           V1=V*MRR
           V2=V*MRI
           Z(I)=V              !=(2\pi/\lambda)*r
           ZR(I)=V1            !=(2\pi/\lambda)*r*MRR
           ZI(I)=V2            !=(2\pi/\lambda)*r*MRI
   10 CONTINUE
      IF (NMAX.GT.NPN1) PRINT 9000,NMAX,NPN1
      IF (NMAX.GT.NPN1) STOP
 9000 FORMAT(' NMAX = ',I2,', i.e., greater than ',I3)
* 
* TA is the ``max. size parameter", MAX(2*PI*SQRT(RI)/LAMBDA)

      TB=TA*DSQRT(MRR*MRR+MRI*MRI)     !=TA*EPSIN
      TB=DMAX1(TB,DFLOAT(NMAX))
*      
      NNMAX1=1.2D0*DSQRT(DMAX1(TA,DFLOAT(NMAX)))+3D0
      NNMAX2=(TB+4D0*(TB**0.33333D0)+1.2D0*DSQRT(TB))  !Wiscombe bound
cc      NNMAX2=NNMAX2-NMAX+5
*
* generate arrays of Bessel functions at NGAUSS GIF division
* points and store them in the common block /CBESS/
*
      CALL BESS(Z,ZR,ZI,NG,NMAX,NNMAX1,NNMAX2)
*
      RETURN
      END
 
C**********************************************************************
 
      SUBROUTINE RSP1 (X,NG,NGAUSS,REV,EPS,R,DR)
C--------/---------/---------/---------/---------/---------/---------/-- 
C >>> X,NG,NGAUSS,REV,EPS
C <<< R,DR
C=========================
C   Activated for NP=-1
C
C   Calculation of the functions 
C              R(I)=r(y)**2 and DR(I)=((d/dy)r(y))/r(y) 
C   for an oblate/prolate spheroids droplet specified by the parameters 
C   REV and  EPS at NGAUSS Gauss integration formula (GIF) division points 
C   in the integral over theta. Here Y=ACOS(X)=THETA.
C
C      r(\theta,\phi)=a\left[\sin^2\theta + (a^2/b^2)\cos^2\theta]^{-1/2}
C
C   X - GIF division points \cos\theta_j -  Y = arccos X 
C   REV ... equal-volume-sphere radius
C   EPS ... the ratio of the horizontal to rotational axes.  EPS is 
C           larger than 1 for oblate spheroids and smaller than 1 for        
C           prolate spheroids. 
C   NGAUSS ... the number of GIF division points
C   NG=2*NGAUSS       
C                                                
C   1.LE.I.LE.NGAUSS                                                    
C--------/---------/---------/---------/---------/---------/---------/--
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 X(NG),R(NG),DR(NG)
      
      A=REV*EPS**(1D0/3D0)
      AA=A*A
      EE=EPS*EPS
      EE1=EE-1D0
      
      DO 50 I=1,NGAUSS
          C=X(I)
          CC=C*C
          SS=1D0-CC
          S=DSQRT(SS)                 !=\sin\theta
          RR=1D0/(SS+EE*CC)
          R(I)=AA*RR
          R(NG-I+1)=R(I)
          DR(I)=RR*C*S*EE1
          DR(NG-I+1)=-DR(I)
   50 CONTINUE

      RETURN
      END


C**********************************************************************
 
      SUBROUTINE RSP2 (X,NG,REV,EPS,N,R,DR)
C--------/---------/---------/---------/---------/---------/---------/-- 
C >>> X,NG,REV,EPS,N
C <<< R,DR
C=========================
C   Activated for NP.gt.0
C
C   Calculation of the functions R(I)=r(y)**2 and                     
C   DR(I)=((d/dy)r(y))/r(y) for a Chebyshev particle          
C   specified by the parameters REV, EPS, and N,
C
C       r(\theta,\phi)=r_0[1+\eps T_n(\cos\theta)]    (*)
C
C   EPS ... deformation parameter of a Chebyshev particle; |EPS|<1  
C   N   ... the degree of the Chebyshev polynomial 
C   All Chebyshev particles with N.GE.2 become partially concave
C   as the absolute value of the deformation parameter EPS increases
C   and exhibit surface roughness in the form of waves running
C   completely around the particle.
C
C   X - GIF division points \cos\theta_j -  Y = arccos X 
C   REV ... equal-volume-sphere radius r_ev
C   NGAUSS ... the number of GIF division points
C   NG=2*NGAUSS  
C                                                       
C   1.LE.I.LE.NGAUSS                                                  
C 
C--------/---------/---------/---------/---------/---------/---------/--
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 X(NG),R(NG),DR(NG)
      DNP=DFLOAT(N)
      DN=DNP*DNP
      DN4=DN*4D0
      EP=EPS*EPS
      A=1D0+1.5D0*EP*(DN4-2D0)/(DN4-1D0)
      I=(DNP+0.1D0)*0.5D0
      I=2*I
      IF (I.EQ.N) A=A-3D0*EPS*(1D0+0.25D0*EP)/
     *              (DN-1D0)-0.25D0*EP*EPS/(9D0*DN-1D0)
      R0=REV*A**(-1D0/3D0)
      DO 50 I=1,NG
         XI=DACOS(X(I))*DNP
         RI=R0*(1D0+EPS*DCOS(XI))    !the Chebyshev shape function (*)
         R(I)=RI*RI
         DR(I)=-R0*EPS*DNP*DSIN(XI)/RI
c        WRITE(NOUT,*) I,R(I),DR(I)
   50 CONTINUE

      RETURN
      END
 
C**********************************************************************
 
      SUBROUTINE RSP3 (X,NG,NGAUSS,REV,EPS,R,DR)
C--------/---------/---------/---------/---------/---------/---------/-- 
C >>> X,NG,NGAUSS,REV,EPS
C <<< R,DR
C=========================
C   Activated for NP=-2
C
C   Calculation of the functions R(I)=r(y)**2 and                     
C   DR(I)=((d/dy)r(y))/r(y) for an oblate/prolate cylinder  
C   specified by the parameters REV and EPS  at NGAUSS  Gauss 
C   integration points in the integral over theta.  
C
C   X - GIF division points \cos\theta_j -  Y = arccos X 
C   REV ... equal-volume-sphere radius r_ev
C   EPS ... the ratio of the cylinder diameter to its length
C   H   ... half-length of the cylinder
C   A=H*EPS  ... cylinder radius   ====>
C
C   4*PI*REV**3/3=2*H*PI*A**2=2*PI*H**3*EPS**2 <====> 
C                H=REV*( (2D0/(3D0*EPS*EPS))**(1D0/3D0) )  
C
C
C   NGAUSS ... the number of GIF division points
C   NG=2*NGAUSS  
C                                                       
C   1.LE.I.LE.NGAUSS                                                  
C  
C--------/---------/---------/---------/---------/---------/---------/--
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 X(NG),R(NG),DR(NG)

* Determine half-length of the cylinder
      H=REV*( (2D0/(3D0*EPS*EPS))**(1D0/3D0) )
      
* Determine cylinder radius:
      A=H*EPS
      
      DO 50 I=1,NGAUSS
         CO=-X(I)
         SI=DSQRT(1D0-CO*CO)
         
         IF (SI/CO.GT.A/H) GO TO 20 

* Along the plane cuts:  
  
         RAD=H/CO
         RTHET=H*SI/(CO*CO)
         GO TO 30
         
* Along the circular surface:         
   20    CONTINUE
         RAD=A/SI
         RTHET=-A*CO/(SI*SI)
cc         RAD=1.D-10
cc         RTHET=0.D0
         
   30    R(I)=RAD*RAD
         R(NG-I+1)=R(I)          !using mirror symmetry
                  
         DR(I)=-RTHET/RAD
         DR(NG-I+1)=-DR(I)       !using mirror symmetry
         
   50 CONTINUE

      RETURN
      END

      SUBROUTINE RSP4 (X,NG,REV,R,DR)
C--------/---------/---------/---------/---------/---------/---------/-- 
C >>> X,NG
C <<< R,DR
C=========================
C   Activated for NP=-3
C
C   Calculation of the functions R(I)=r(y)**2 and                   
C   DR(I)=((d/dy)r(y))/r(y) for a distorted                        
C   droplet (generalized Chebyshev particle) specified by the 
C   parameters REV and c_n (Chebyshev expansion coefficients).
C   The coefficients of the Chebyshev polynomial expansion are
C   specified in the subroutine DROP.
C 
C   X - GIF division points \cos\theta_j -  Y = arccos X 
C   REV ... equal-volume-sphere radius  r_ev
C   NGAUSS ... the number of GIF division points
C   NG=2*NGAUSS  
C                                              
C   1.LE.I.LE.NGAUSS                                                  
C
C--------/---------/---------/---------/---------/---------/---------/--
      PARAMETER (NC=10)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 X(NG),R(NG),DR(NG),C(0:NC)
      COMMON /CDROP/ C,R0V
      R0=REV*R0V
      DO I=1,NG
         XI=DACOS(X(I))
         RI=1D0+C(0)
         DRI=0D0
         DO N=1,NC
            XIN=XI*N
            RI=RI+C(N)*DCOS(XIN)
            DRI=DRI-C(N)*N*DSIN(XIN)
         ENDDO
         RI=RI*R0
         DRI=DRI*R0
         R(I)=RI*RI
         DR(I)=DRI/RI
c        WRITE(NOUT,*) I,R(I),DR(I)
      ENDDO

      RETURN
      END
      

C**********************************************************************
 
      SUBROUTINE RSP5 (X,NG,REV,EPS,R,DR)
C--------/---------/---------/---------/---------/---------/---------/-- 
C >>> X,NG,REV,EPS
C <<< R,DR
C=========================
C   Activated for NP=-4
C
C   SIMILAR TO RSPI5, EXCEPT FOR THAT THE PLANE CUT IS ON THE SPHERE "TOP"
C   ===> cosine of the separation angle has the same magnitude as in RSPI5,
C        but is always positive in the present case !!!
C
C   Calculation of the functions 
C              R(I)=r(y)**2 and DR(I)=((d/dy)r(y))/r(y) 
C   for a sphere cut by the plane specified by the parameters 
C   REV and  EPS at NGAUSS Gauss integration formula (GIF) division points 
C   in the integral over theta. Here Y=ACOS(X)=THETA and EPS=2*R_0/H,
C   where R_0 is the radius of the original uncut sphere, whereas H
C   is the height (along the axial symmetry axis) of the resulting
C   cut sphere. 
C
C   The origin of coordinates is located along the axis of symmetry
C                  midway the plane and sphere bottom.
C
C           ===> Note that always EPS.GT.1
C   ===
C   X - GIF division points \cos\theta_j -  Y = arccos X 
C   REV ... the radius of the original uncut sphere 
C   EPS ...  H is the height (along the axial symmetry axis) 
C            of the resulting cut sphere. Note that always EPS.LT.2*REV
C   THETA0 ... a COSINE of the separation angle between two different 
C              functional dependences of r(\theta), that along the sphere 
C              surface and that along the plane surface
C   NG=2*NGAUSS ... the number of GIF division points     
C                                                
C   1.LE.I.LE.NGAUSS                                                    
C--------/---------/---------/---------/---------/---------/---------/--
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 X(NG),R(NG),DR(NG)
      
      IF (EPS.GE.2.d0*REV) THEN
      WRITE(6,*)'Invalid parameters for a cut sphere!'
      WRITE(6,*)'Execution stopped!'
      STOP
      END IF
 
      THETA0=1.D0/SQRT(8.D0*REV/EPS-3.D0)  !cosine of the separation angle
    
      DO 50 I=1,NG
          
          CO=X(I)
          CC=CO*CO
          SS=1.D0-CC
          SI=DSQRT(SS)                    !=\sin\theta      

          IF (CO.LT.THETA0) THEN        ! r(\theta) along the sphere surface
          
          A=REV-EPS/2.D0
          
          RAD=A*CO+SQRT(REV**2-(A*SI)**2)
          RTHET=-A*SI - CO*SI*A**2/SQRT(REV**2-(A*SI)**2)
cc          RAD=1.D-10
cc          RTHET=0.D0                       

          
          ELSE IF (CO.GE.THETA0) THEN   ! r(\theta) along the plane surface
                                        !    (CO positive)             
          RAD=EPS/(2.D0*CO)
          RTHET=EPS*SI/(2.D0*CO**2)       
                   
          END IF 
          
          DR(I)=RTHET/RAD
          R(I)=RAD*RAD
                               
   50 CONTINUE

      RETURN
      END

C**********************************************************************
 
      SUBROUTINE RSPI5 (X,NG,REV,EPS,R,DR)
C--------/---------/---------/---------/---------/---------/---------/-- 
C >>> X,NG,REV,EPS
C <<< R,DR
C=========================
C   Activated for NP=-5
C
C   SIMILAR TO RSP5, EXCEPT FOR THAT THE PLANE CUT IS ON THE SPHERE "BOTTOM"
C   ===> cosine of the separation angle has the same magnitude as in RSP5,
C        but is always negative in the present case !!!
C
C   Calculation of the functions 
C              R(I)=r(y)**2 and DR(I)=((d/dy)r(y))/r(y) 
C   for a sphere cut by the plane specified by the parameters 
C   REV and  EPS at NGAUSS Gauss integration formula (GIF) division points 
C   in the integral over theta. Here Y=ACOS(X)=THETA and EPS=2*R_0/H,
C   where R_0 is the radius of the original uncut sphere, whereas H
C   is the height (along the axial symmetry axis) of the resulting
C   cut sphere. 
C
C   The origin of coordinates is located along the axis of symmetry
C                  midway the plane and sphere top.
C
C                  ===>  Note that always EPS.GT.1
C   ===
C   X - GIF division points \cos\theta_j -  Y = arccos X 
C   REV ... the radius of the original uncut sphere 
C   EPS ...  H is the height (along the axial symmetry axis) 
C            of the resulting cut sphere. Note that always EPS.LT.2*REV
C   THETA0 ... a COSINE of the separation angle between two different 
C              functional dependences of r(\theta), that along the sphere 
C              surface and that along the plane surface
C   NG=2*NGAUSS ... the number of GIF division points     
C                                                
C   1.LE.I.LE.NGAUSS                                                    
C--------/---------/---------/---------/---------/---------/---------/--
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 X(NG),R(NG),DR(NG)
      
      IF (EPS.GE.2.d0*REV) THEN
      WRITE(6,*)'Invalid parameters for a cut sphere!'
      WRITE(6,*)'Execution stopped!'
      STOP
      END IF
 
      THETA0=-1.D0/SQRT(8.D0*REV/EPS-3.D0)  !cosine of the separation angle
                                            !is always negative in the present
                                            !case  
      DO 50 I=1,NG
          
          CO=X(I)
          CC=CO*CO
          SS=1.D0-CC
          SI=DSQRT(SS)                  !=\sin\theta      

          IF (CO.GT.THETA0) THEN        !r(\theta) along the sphere surface
          
          A=REV-EPS/2.D0
          
          RAD=-A*CO+SQRT(REV**2-(A*SI)**2)
          RTHET=A*SI - CO*SI*A**2/SQRT(REV**2-(A*SI)**2)
cc          RAD=1.D-10
cc          RTHET=0.D0                       
          
          ELSE IF (CO.LE.THETA0) THEN   ! r(\theta) along the plane surface
                                        !        (CO negative)          
          RAD=-EPS/(2.D0*CO)
          RTHET=-EPS*SI/(2.D0*CO**2)       
                   
          END IF 
          
          DR(I)=RTHET/RAD
          R(I)=RAD*RAD
                                    
   50 CONTINUE

      RETURN
      END 

C*********************************************************************
 
      SUBROUTINE BESS (X,XR,XI,NG,NMAX,NNMAX1,NNMAX2)
C--------/---------/---------/---------/---------/---------/---------/--
C >>> X,XR,XI,NG,NMAX,NNMAX1,NNMAX2
C <<< Output J,Y,JR,JI,DJ,DY,DJR,DJI  to common block CBESS
C==========================================================
C  Generates Bessel functions for each Gauss integration point
C
C  X =(2\pi/\lambda)*r
C  XR=(2\pi/\lambda)*r*MRR, MRR ... real part of the rel. refractive index
C  XI=(2\pi/\lambda)*r*MRI, MRI ... imag. part of the rel. refractive index
C  NG=2*NGAUSS or 60 ... the number of Gauss integration points
C  J,Y,JR,JI ... arrays of Bessel functions
C  DJ,DY,DJR,DJI  ... arrays of Bessel functions derivatives of the form
C                           [xf(x)]'/x                   (A)
C                 where prime denotes derivative with respect to x.
C                 (Note that Bessel function derivatives enter Eqs. (39)
C                  \cite{TKS} only in the (A) combination!!!!)
C  NMAX   ... angular momentum cutoff
C  NNMAX1 ... angular momentum cutoff - DETERMINES NUMERICAL ACCURACY 
C  NNMAX2 ... angular momentum cutoff - DETERMINES NUMERICAL ACCURACY 
C--------/---------/---------/---------/---------/---------/---------/--
      INCLUDE 'ampld.par.f'
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 X(NG),XR(NG),XI(NG),
     *        J(NPNG2,NPN1),Y(NPNG2,NPN1),JR(NPNG2,NPN1),
     *        JI(NPNG2,NPN1),DJ(NPNG2,NPN1),DY(NPNG2,NPN1),
     *        DJR(NPNG2,NPN1),DJI(NPNG2,NPN1),
     *        AJ(NPN1),AY(NPN1),AJR(NPN1),AJI(NPN1),
     *        ADJ(NPN1),ADY(NPN1),ADJR(NPN1),
     *        ADJI(NPN1)
      COMMON /CBESS/ J,Y,JR,JI,DJ,DY,DJR,DJI    !arrays of generated Bessel functions
* 
      DO 10 I=1,NG
           XX=X(I)
*
           CALL RJB(XX,AJ,ADJ,NMAX,NNMAX1)
           CALL RYB(XX,AY,ADY,NMAX)
*
           YR=XR(I)
           YI=XI(I)
*
           CALL CJB(YR,YI,AJR,AJI,ADJR,ADJI,NMAX,2)
*
           DO 10 N=1,NMAX
                J(I,N)=AJ(N)
                Y(I,N)=AY(N)
                JR(I,N)=AJR(N)
                JI(I,N)=AJI(N)
                DJ(I,N)=ADJ(N)
                DY(I,N)=ADY(N)
                DJR(I,N)=ADJR(N)
                DJI(I,N)=ADJI(N)
   10 CONTINUE

      RETURN
      END
 
C**********************************************************************
 
      SUBROUTINE RJB(X,Y,U,NMAX,NNMAX)
C--------/---------/---------/---------/---------/---------/---------/--
C >>> EPS
C <<< RAT
C=================
C  X =(2\pi/\lambda)*r
C  Y ... 
C  NMAX - angular momentum cutoff
C  NNMAX - angular momentum cutoff - DETERMINES NUMERICAL ACCURACY  
C--------/---------/---------/---------/---------/---------/---------/--
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 Y(NMAX),U(NMAX),Z(800)
*
      L=NMAX+NNMAX
      XX=1D0/X
      Z(L)=1D0/(DFLOAT(2*L+1)*XX)
      L1=L-1
      DO 5 I=1,L1
         I1=L-I
         Z(I1)=1D0/(DFLOAT(2*I1+1)*XX-Z(I1+1))
    5 CONTINUE
      Z0=1D0/(XX-Z(1))
      Y0=Z0*DCOS(X)*XX
      Y1=Y0*Z(1)
      U(1)=Y0-Y1*XX
      Y(1)=Y1
      DO 10 I=2,NMAX
         YI1=Y(I-1)
         YI=YI1*Z(I)
         U(I)=YI1-DFLOAT(I)*YI*XX
         Y(I)=YI
   10 CONTINUE

      RETURN
      END
 
C**********************************************************************
 
      SUBROUTINE RYB(X,Y,V,NMAX)
C--------/---------/---------/---------/---------/---------/---------/--
C >>> EPS
C <<< RAT
C=================
C  X =(2\pi/\lambda)*r
C  NMAX - angular momentum cutoff
C--------/---------/---------/---------/---------/---------/---------/--
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 Y(NMAX),V(NMAX)
*
      C=DCOS(X)
      S=DSIN(X)
      X1=1D0/X
      X2=X1*X1
      X3=X2*X1
      Y1=-C*X2-S*X1
      Y(1)=Y1
      Y(2)=(-3D0*X3+X1)*C-3D0*X2*S
      NMAX1=NMAX-1
      DO 5 I=2,NMAX1
    5     Y(I+1)=DFLOAT(2*I+1)*X1*Y(I)-Y(I-1)
      V(1)=-X1*(C+Y1)
      DO 10 I=2,NMAX
  10       V(I)=Y(I-1)-DFLOAT(I)*X1*Y(I)
      RETURN
      END

C**********************************************************************
 
      SUBROUTINE CJB (XR,XI,YR,YI,UR,UI,NMAX,NNMAX)
C--------/---------/---------/---------/---------/---------/---------/--
C                                                                     
C   CALCULATION OF SPHERICAL BESSEL FUNCTIONS OF THE FIRST KIND       
C   J=JR+I*JI OF COMPLEX ARGUMENT X=XR+I*XI OF ORDERS FROM 1 TO NMAX  
C   BY USING BACKWARD RECURSION. PARAMETER NNMAX DETERMINES NUMERICAL  
C   ACCURACY. U=UR+I*UI - FUNCTION (1/X)(D/DX)(X*J(X))=J(X)/X + J'(X)                
C
C  XR=(2\pi/\lambda)*r*MRR, MRR ... real part of the rel. refractive index
C  XI=(2\pi/\lambda)*r*MRI, MRI ... imag. part of the rel. refractive index
C
C   NMAX  - angular momentum cutoff 
C   NNMAX - angular momentum cutoff - DETERMINES NUMERICAL ACCURACY                
                                                 
C--------/---------/---------/---------/---------/---------/---------/--
      INCLUDE 'ampld.par.f'
      IMPLICIT REAL*8 (A-H,O-Z)
      
      REAL*8 YR(NMAX),YI(NMAX),UR(NMAX),UI(NMAX)
      REAL*8 CYR(NPN1),CYI(NPN1),CZR(1200),CZI(1200)
c     *       CUR(NPN1),CUI(NPN1)
*
      L=NMAX+NNMAX
      XRXI=1D0/(XR*XR+XI*XI)
      CXXR=XR*XRXI             !Re [1/(XR+i*XI)]
      CXXI=-XI*XRXI            !Im [1/(XR+i*XI)] 
      QF=1D0/DFLOAT(2*L+1)
      CZR(L)=XR*QF
      CZI(L)=XI*QF
      L1=L-1
      DO I=1,L1
         I1=L-I
         QF=DFLOAT(2*I1+1)
         AR=QF*CXXR-CZR(I1+1)
         AI=QF*CXXI-CZI(I1+1)
         ARI=1D0/(AR*AR+AI*AI)
         CZR(I1)=AR*ARI
         CZI(I1)=-AI*ARI
      ENDDO   
      
      AR=CXXR-CZR(1)
      AI=CXXI-CZI(1)
      ARI=1D0/(AR*AR+AI*AI)
      CZ0R=AR*ARI
      CZ0I=-AI*ARI
      CR=DCOS(XR)*DCOSH(XI)
      CI=-DSIN(XR)*DSINH(XI)
      AR=CZ0R*CR-CZ0I*CI
      AI=CZ0I*CR+CZ0R*CI
      CY0R=AR*CXXR-AI*CXXI
      CY0I=AI*CXXR+AR*CXXI
      CY1R=CY0R*CZR(1)-CY0I*CZI(1)
      CY1I=CY0I*CZR(1)+CY0R*CZI(1)
      AR=CY1R*CXXR-CY1I*CXXI
      AI=CY1I*CXXR+CY1R*CXXI
      CU1R=CY0R-AR
      CU1I=CY0I-AI
      CYR(1)=CY1R
      CYI(1)=CY1I
c      CUR(1)=CU1R
c      CUI(1)=CU1I
      YR(1)=CY1R
      YI(1)=CY1I
      UR(1)=CU1R
      UI(1)=CU1I
      
      DO I=2,NMAX
         QI=DFLOAT(I)
         CYI1R=CYR(I-1)
         CYI1I=CYI(I-1)
         CYIR=CYI1R*CZR(I)-CYI1I*CZI(I)
         CYII=CYI1I*CZR(I)+CYI1R*CZI(I)
         AR=CYIR*CXXR-CYII*CXXI            !Re [J/(XR+i*XI)]
         AI=CYII*CXXR+CYIR*CXXI            !Im [J/(XR+i*XI)]
         CUIR=CYI1R-QI*AR
         CUII=CYI1I-QI*AI
         CYR(I)=CYIR
         CYI(I)=CYII
c         CUR(I)=CUIR
c         CUI(I)=CUII
         YR(I)=CYIR
         YI(I)=CYII
         UR(I)=CUIR
         UI(I)=CUII
      ENDDO 
*  
      RETURN
      END
 
C**********************************************************************
 
      SUBROUTINE TMATR0 (NGAUSS,X,W,AN,ANN,PPI,PIR,PII,R,DR,DDR,
     *                  DRR,DRI,NMAX,NCHECK,NAXSM)
C--------/---------/---------/---------/---------/---------/---------/--
C >>> NGAUSS,X,W,AN,ANN,PPI,PIR,PII,R,DR,DDR,DRR,DRI,NMAX,NCHECK
C <<< common blocks /TMAT99/, /CT/ (for main),  and /CTT/ (for TT)
C=====================
C
C  Determines the T-matrix of an axially symmetric scatterer
C                           for M=0
C
C  NGAUSS - the number of GIF division points
C  X=\cos\theta  - GIF division points
C  W - GIF weights
C  AN(N)=N*(N+1)
C  ANN(l_1,l_2)=\sqrt{\fr{(2 l_1+1)}{l_1(l_1+1)} }
C                       \sqrt{\fr{(2 l_2+1)}{l_2(l_2+1)} }/2
C  NMAX - angular momentum cutoff
C  NCHECK  -  .EQ.0  THEN  NGSS=2*NGAUSS, FACTOR=1D0
C             .EQ.1  THEN  NGSS = NGAUSS, FACTOR=2D0
C  NAXSM   -  .EQ.0 : Gauss abscissas do not have +/- theta symmetry
C             .EQ.1 : Gauss abscissas have +/- theta symmetry 
C  P=DACOS(-1D0)
C  PI=P*2D0/LAM - wave vector
C  PPI=PI*PI
C  PIR=PPI*MRR
C  PII=PPI*MRI
C  R=r^2(\theta)                        for axially symmetric particles
C  DR=[dr(\theta)/(d\theta)]/r(\theta)  for axially symmetric particles
C  DDR=\lambda/[2*\pi*r(\theta)]=1/(k_out*r)
C  DRR=(MRR/(MRR**2+MRI**2))*(\lambda/[2*\pi*r(\theta)])
C                  = Re 1/(k_in*r)
C  DRI=-(MRI/(MRR**2+MRI**2))*(\lambda/[2*\pi*r(\theta)])
C                  = Im 1/(k_in*r)
C  Refractive index outside is assumed to real, whereas inside
C  a scatterer, refractive index is allowed to be complex in general.
C  Consequently, the Bessel function j_l(k_in*r) will in general
C  be complex. The routine below performs Waterman surface integral
C  separately for the real and imaginary parts of the integrand.
C
C--------/---------/---------/---------/---------/---------/---------/--
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'ampld.par.f'
      INTEGER NOUT     
* number of the output unit
      PARAMETER (NOUT=35)
      
      REAL*8  X(NPNG2),W(NPNG2),AN(NPN1),
     *        R(NPNG2),DR(NPNG2),SIG(NPN2),
     *        J(NPNG2,NPN1),Y(NPNG2,NPN1),
     *        JR(NPNG2,NPN1),JI(NPNG2,NPN1),DJ(NPNG2,NPN1),
     *        DY(NPNG2,NPN1),DJR(NPNG2,NPN1),
     *        DJI(NPNG2,NPN1),DDR(NPNG2),DRR(NPNG2),
     *        D1(NPNG2,NPN1),D2(NPNG2,NPN1),
     *        DRI(NPNG2),RR(NPNG2),
     *        DV1(NPN1),DV2(NPN1)
 
      REAL*8  R11(NPN1,NPN1),R12(NPN1,NPN1),
     *        R21(NPN1,NPN1),R22(NPN1,NPN1),
     *        I11(NPN1,NPN1),I12(NPN1,NPN1),
     *        I21(NPN1,NPN1),I22(NPN1,NPN1),
     *        RG11(NPN1,NPN1),RG12(NPN1,NPN1),
     *        RG21(NPN1,NPN1),RG22(NPN1,NPN1),
     *        IG11(NPN1,NPN1),IG12(NPN1,NPN1),
     *        IG21(NPN1,NPN1),IG22(NPN1,NPN1),
     *        ANN(NPN1,NPN1),
     *        QR(NPN2,NPN2),QI(NPN2,NPN2),
     *        RGQR(NPN2,NPN2),RGQI(NPN2,NPN2),
     *        TQR(NPN2,NPN2),TQI(NPN2,NPN2),
     *        TRGQR(NPN2,NPN2),TRGQI(NPN2,NPN2)
cc      REAL*8 TR1(NPN2,NPN2),TI1(NPN2,NPN2)
*      
      COMMON /TMAT99/ 
     &            R11,R12,R21,R22,I11,I12,I21,I22,RG11,RG12,RG21,RG22,
     &            IG11,IG12,IG21,IG22           !only between TMATR routines
      COMMON /CBESS/ J,Y,JR,JI,DJ,DY,DJR,DJI
cc      COMMON /CT/ TR1,TI1                       !output from TT routine
      COMMON /CTT/ QR,QI,RGQR,RGQI              !input for TT routine
*      
      MM1=1
      NNMAX=NMAX+NMAX
      NG=2*NGAUSS
      NGSS=NG
      FACTOR=1D0
*
      IF (NCHECK.EQ.1) THEN         !Theta=pi/2 is scatterer mirror symmetry plane
            NGSS=NGAUSS
            FACTOR=2D0
         ELSE                       !Theta=pi/2 is not a scatterer mirror symmetry plane
            CONTINUE
      ENDIF
*
      SI=1D0
      DO 5 N=1,NNMAX
           SI=-SI
           SIG(N)=SI              !=(-1)**N
    5 CONTINUE
*
* Assigning Wigner d-matrices - assuming mirror symmetry 
* in the \theta=\pi/2 plane:

      DO 25 I=1,NGAUSS

         I1=NGAUSS-I+1 
cc         I2=NGAUSS+I 
         IF (NCHECK.EQ.0) I2=NGAUSS+I 
*
         CALL VIG ( X(I1), NMAX, 0, DV1, DV2)
*
         DO N=1,NMAX

            DD1=DV1(N)
            DD2=DV2(N)
            D1(I1,N)=DD1
            D2(I1,N)=DD2


* If theta=pi/2 is not a scatterer mirror symmetry plane but 
* Gauss abscissas are still chosen symmetrically:

         IF ((NCHECK.EQ.0).AND.(NAXSM.EQ.1)) THEN      
        
* using (4.2.4) and (4.2.6) of {Ed},  
*           d_{0m}^{(l)}(\pi-\theta) = (-1)^{l+m} d_{0m}^{(l)}(\theta)
* (4.2.5) of {Ed}:                   = (-1)^{l} d_{0 -m}^{(l)}(\theta)

            SI=SIG(N)                  !=(-1)**N           
            D1(I2,N)=DD1*SI       
            D2(I2,N)=-DD2*SI
           
         END IF
            
         ENDDO 
         
* If neither scatterer nor Gauss abscissas have theta=pi/2 
* as a mirror symmetry plane:   
       
         IF ((NCHECK.EQ.0).AND.(NAXSM.EQ.0)) THEN                                 
*
         CALL VIG ( X(I2), NMAX, 0, DV1, DV2)
*
          DO N=1,NMAX           
            DD1=DV1(N)
            DD2=DV2(N)
            D1(I2,N)=DD1
            D2(I2,N)=DD2
          ENDDO     
           
         END IF
                                                        
   25 CONTINUE
*
*  Assigning r^2(\theta)*weight product:

      DO 40 I=1,NGSS
           RR(I)=W(I)*R(I)
           
cc           if (dr(i).eq.0.d0) RR(I)=0.d0   !temporarily only
           
   40 CONTINUE
* 
      DO 300 N1=MM1,NMAX
           AN1=AN(N1)
           DO 300 N2=MM1,NMAX
                AN2=AN(N2)
                
                AR12=0D0
                AR21=0D0
                AI12=0D0
                AI21=0D0
                GR12=0D0
                GR21=0D0
                GI12=0D0
                GI21=0D0
                
c        OPEN(NOUT+3,FILE='surfint.dat')   !Gauss convergence check

                IF (NCHECK.EQ.1.AND.SIG(N1+N2).LT.0D0) GO TO 205
*
* Gauss integration loop:
*
                DO 200 I=1,NGSS    !=NGAUSS   if NCHECK.EQ.1
                                   !=2*NGAUSS if NCHECK.EQ.0 
                
                    D1N1=D1(I,N1)
                    D2N1=D2(I,N1)
                    D1N2=D1(I,N2)
                    D2N2=D2(I,N2)
                    A12=D1N1*D2N2
                    A21=D2N1*D1N2
                    A22=D2N1*D2N2
c                    AA1=A12+A21
 
* Vector spherical harmonics:
C  Since refractive index is allowed to be complex in general,
C  the Bessel function j_l(k_in*r) is complex. The code below 
C  performs a separation of the complex integrand in Waterman's
C  surface integral into its respective real and imaginary 
C  parts:

* Bessel functions of the exterior argument:

                    QJ1=J(I,N1)
                    QY1=Y(I,N1)
                    QDJ1=DJ(I,N1)
                    QDY1=DY(I,N1)
                    
* Bessel functions of the interior argument:                                        
                    
                    QJR2=JR(I,N2)
                    QJI2=JI(I,N2)
                    QDJR2=DJR(I,N2)
                    QDJI2=DJI(I,N2)                    
*_____________________    
                
* Re and Im of j_{n2}(k_{in}r) j_{n1}(k_{out}r): 

                    C1R=QJR2*QJ1
                    C1I=QJI2*QJ1
                    
* Re and Im of j_{n2}(k_{in}r) h_{n1}(k_{out}r):
                   
                    B1R=C1R-QJI2*QY1
                    B1I=C1I+QJR2*QY1

* Re and Im of j_{n2}(k_{in}r) [k_{out}r j_{n1}(k_{out}r)]'/(k_{out}r): 
 
                    C2R=QJR2*QDJ1
                    C2I=QJI2*QDJ1

* Re and Im of j_{n2}(k_{in}r) [k_{out}r h_{n1}(k_{out}r)]'/(k_{out}r):
                   
                    B2R=C2R-QJI2*QDY1
                    B2I=C2I+QJR2*QDY1
 
                    DDRI=DDR(I)               !1/(k_{out}r)

* Re and Im of [1/(k_{out}r)]*j_{n2}(k_{in}r) j_{n1}(k_{out}r)

                    C3R=DDRI*C1R
                    C3I=DDRI*C1I
                    
* Re and Im of [1/(k_{out}r)]*j_{n2}(k_{in}r) h_{n1}(k_{out}r):
                    
                    B3R=DDRI*B1R
                    B3I=DDRI*B1I

* Re and Im of  [k_{in}r j_{n2}(k_{in}r)]'/(k_{in}r) 
*                          * j_{n1}(k_{out}r): 
 
                    C4R=QDJR2*QJ1
                    C4I=QDJI2*QJ1
                    
* Re and Im of [k_{in}r j_{n2}(k_{in}r)]'/(k_{in}r) 
*                          *  h_{n1}(k_{out}r): 
                   
                    B4R=C4R-QDJI2*QY1
                    B4I=C4I+QDJR2*QY1
 
                    DRRI=DRR(I)               !Re[1/(k_{in}r)]
                    DRII=DRI(I)               !Im[1/(k_{in}r)]
                    
* Re and Im of [1/(k_{in}r)] j_{n2}(k_{in}r) j_{n1}(k_{out}r):
                    
                    C5R=C1R*DRRI-C1I*DRII
                    C5I=C1I*DRRI+C1R*DRII
                    
                    
* Re and Im of [1/(k_{in}r)] j_{n2}(k_{in}r) h_{n1}(k_{out}r):

                    B5R=B1R*DRRI-B1I*DRII
                    B5I=B1I*DRRI+B1R*DRII

*%%%%%%%  Forming integrands of J-matrices (J^{11}=J^{22}=0 for m=0): %%%%%%%%
 
                    URI=DR(I)        !dr/(d\theta)
                    RRI=RR(I)        !w(i)*r^2(\theta)
                    
* w(i)*r^2(\theta)*D2N1*D2N2:
                    F1=RRI*A22      !prefactor containing r^2(\theta)<->hat{r} part

                    
* N1*(N1+1)*w(i)*r(\theta)*[dr/(d\theta)]*D1N1*D2N2:                    
                    F2=RRI*URI*AN1*A12     !prefactor containing r(\theta)*[dr/(d\theta)]
                                           !hat{theta} part 
                       
                    
                    AR12=AR12+F1*B2R+F2*B3R        !~Re J^{12}
                    AI12=AI12+F1*B2I+F2*B3I        !~Im J^{12}
                    
                    GR12=GR12+F1*C2R+F2*C3R        !~Re Rg J^{12}
                    GI12=GI12+F1*C2I+F2*C3I        !~Im Rg J^{12}

* N2*(N2+1)*w(i)*r(\theta)*[dr/(d\theta)]*D2N1*D1N2: 
                    F2=RRI*URI*AN2*A21     !prefactor containing r(\theta)*[dr/(d\theta)]
                                           !hat{theta} part
                    
                    AR21=AR21+F1*B4R+F2*B5R        !~Re J^{21}
                    AI21=AI21+F1*B4I+F2*B5I        !~Im J^{21}
                    
                    GR21=GR21+F1*C4R+F2*C5R        !~Re Rg J^{21}
                    GI21=GI21+F1*C4I+F2*C5I        !~Im Rg J^{21}
                    
  200           CONTINUE               !end of Gauss integration
  
c                write(nout+3,*)'N1=',N1,'   N2=',N2 
c                write(nout+3,*)'AR12=', AR12 
c                write(nout+3,*)'AI12=', AI12 
c                write(nout+3,*)'AR21=', AR21
c                write(nout+3,*)'AI21=', AI21
c                write(nout+3,*)'GR12=', GR12 
c                write(nout+3,*)'GI12=', GI12 
c                write(nout+3,*)'GR21=', GR21
c                write(nout+3,*)'GI21=', GI21                

*%%%%%%%%%%%%%  Forming J-matrices (J^{11}=J^{22}=0 for m=0):
 
  205           AN12=ANN(N1,N2)*FACTOR
  
                R12(N1,N2)=AR12*AN12
                R21(N1,N2)=AR21*AN12
                I12(N1,N2)=AI12*AN12
                I21(N1,N2)=AI21*AN12
                
                RG12(N1,N2)=GR12*AN12
                RG21(N1,N2)=GR21*AN12
                IG12(N1,N2)=GI12*AN12
                IG21(N1,N2)=GI21*AN12
                
  300 CONTINUE            !end of the loop over angular momenta
 
c      close(nout+3)

*%%%%%%%%%%%%%%%%%%%%%%%  Forming Q and RgQ -matrices
 
      TPIR=PIR                 !Re [1/k_{in}^2]
      TPII=PII                 !Im [1/k_{in}^2]
      TPPI=PPI                 !1/k_{out}^2
 
      NM=NMAX
      
      DO 310 N1=MM1,NMAX
           K1=N1-MM1+1
           KK1=K1+NM
           DO 310 N2=MM1,NMAX
                K2=N2-MM1+1
                KK2=K2+NM
 
                TAR12= I12(N1,N2)
                TAI12=-R12(N1,N2)
                TGR12= IG12(N1,N2)
                TGI12=-RG12(N1,N2)
 
                TAR21=-I21(N1,N2)
                TAI21= R21(N1,N2)
                TGR21=-IG21(N1,N2)
                TGI21= RG21(N1,N2)
 
                TQR(K1,K2)=TPIR*TAR21-TPII*TAI21+TPPI*TAR12
                TQI(K1,K2)=TPIR*TAI21+TPII*TAR21+TPPI*TAI12
                TRGQR(K1,K2)=TPIR*TGR21-TPII*TGI21+TPPI*TGR12
                TRGQI(K1,K2)=TPIR*TGI21+TPII*TGR21+TPPI*TGI12
 
                TQR(K1,KK2)=0D0
                TQI(K1,KK2)=0D0
                TRGQR(K1,KK2)=0D0
                TRGQI(K1,KK2)=0D0
 
                TQR(KK1,K2)=0D0
                TQI(KK1,K2)=0D0
                TRGQR(KK1,K2)=0D0
                TRGQI(KK1,K2)=0D0
 
                TQR(KK1,KK2)=TPIR*TAR12-TPII*TAI12+TPPI*TAR21
                TQI(KK1,KK2)=TPIR*TAI12+TPII*TAR12+TPPI*TAI21
                TRGQR(KK1,KK2)=TPIR*TGR12-TPII*TGI12+TPPI*TGR21
                TRGQI(KK1,KK2)=TPIR*TGI12+TPII*TGR12+TPPI*TGI21
  310 CONTINUE
 
      NNMAX=2*NM
      DO 320 N1=1,NNMAX
           DO 320 N2=1,NNMAX
                QR(N1,N2)=TQR(N1,N2)
                QI(N1,N2)=TQI(N1,N2)
                RGQR(N1,N2)=TRGQR(N1,N2)
                RGQI(N1,N2)=TRGQI(N1,N2)
  320 CONTINUE
  
*%%%%%%%%%%%%%%%%%%%%%%%  Forming resulting T-matrix 
*
* Calculate the product Q^{-1} Rg Q
*
      CALL TT(NMAX,NCHECK)
*
      RETURN
      END
 
C**********************************************************************
 
      SUBROUTINE TMATR (M,NGAUSS,X,W,AN,ANN,S,SS,PPI,PIR,PII,R,DR,DDR,
     *                  DRR,DRI,NMAX,NCHECK,NAXSM)
C--------/---------/---------/---------/---------/---------/---------/--
C >>> NGAUSS,X,W,AN,ANN,S,SS,PPI,PIR,PII,R,DR,DDR,DRR,DRI,NMAX,NCHECK
C <<< common blocks /TMAT99/, /CT/ (for main),  and /CTT/ (for TT) 
C=====================
C
C  Determines the T-matrix of an axially symmetric scatterer
C                           for M.GT.0
C
C  M      - azimuthal number
C  NGAUSS - the number of GIF division points
C  X=\cos\theta  - GIF division points
C  W - GIF weights
C  AN(N)=N*(N+1)
C  ANN(l_1,l_2)=\sqrt{\fr{(2 l_1+1)}{l_1(l_1+1)} }
C                       \sqrt{\fr{(2 l_2+1)}{l_2(l_2+1)} }/2
C  NMAX - angular momentum cutoff
C  NCHECK  -  .EQ.0  THEN  NGSS=2*NGAUSS, FACTOR=1D0
C             .EQ.1  THEN  NGSS = NGAUSS, FACTOR=2D0
C  NAXSM   -  .EQ.0 : Gauss abscissas do not have +/- theta symmetry
C             .EQ.1 : Gauss abscissas have +/- theta symmetry  
C  NCHECK - specifies whether NG=2*NGAUSS or otherwise
C  P=DACOS(-1D0)
C  PI=P*2D0/LAM - wave vector
C  PPI=PI*PI
C  PIR=PPI*MRR
C  PII=PPI*MRI
C  R=r^2(\theta)                       for axially symmetric particles
C  DR=[dr(\theta)/(d\theta)]/r(\theta) for axially symmetric particles
C  DDR=\lambda/[2*\pi*r(\theta)]=1/(k_out*r)
C  DRR=(MRR/(MRR**2+MRI**2))*(\lambda/[2*\pi*r(\theta)])
C                  = Re 1/(k_in*r)
C  DRI=-(MRI/(MRR**2+MRI**2))*(\lambda/[2*\pi*r(\theta)])
C                  = Im 1/(k_in*r)
C  Refractive index outside is assumed to real, whereas inside
C  a scatterer, refractive index is allowed to be complex in general.
C  Consequently, the Bessel function j_l(k_in*r) will in general
C  be complex. The routine below performs Waterman surface integral
C  separately for the real and imaginary parts of the integrand.
C
C--------/---------/---------/---------/---------/---------/---------/--
      INCLUDE 'ampld.par.f'
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8  X(NPNG2),W(NPNG2),AN(NPN1),S(NPNG2),SS(NPNG2),
     *        R(NPNG2),DR(NPNG2),SIG(NPN2),
     *        J(NPNG2,NPN1),Y(NPNG2,NPN1),
     *        JR(NPNG2,NPN1),JI(NPNG2,NPN1),DJ(NPNG2,NPN1),
     *        DY(NPNG2,NPN1),DJR(NPNG2,NPN1),
     *        DJI(NPNG2,NPN1),DDR(NPNG2),DRR(NPNG2),
     *        D1(NPNG2,NPN1),D2(NPNG2,NPN1),
     *        DRI(NPNG2),DS(NPNG2),DSS(NPNG2),RR(NPNG2),
     *        DV1(NPN1),DV2(NPN1)
 
      REAL*8  R11(NPN1,NPN1),R12(NPN1,NPN1),
     *        R21(NPN1,NPN1),R22(NPN1,NPN1),
     *        I11(NPN1,NPN1),I12(NPN1,NPN1),
     *        I21(NPN1,NPN1),I22(NPN1,NPN1),
     *        RG11(NPN1,NPN1),RG12(NPN1,NPN1),
     *        RG21(NPN1,NPN1),RG22(NPN1,NPN1),
     *        IG11(NPN1,NPN1),IG12(NPN1,NPN1),
     *        IG21(NPN1,NPN1),IG22(NPN1,NPN1),
     *        ANN(NPN1,NPN1),
     *        QR(NPN2,NPN2),QI(NPN2,NPN2),
     *        RGQR(NPN2,NPN2),RGQI(NPN2,NPN2),
     *        TQR(NPN2,NPN2),TQI(NPN2,NPN2),
     *        TRGQR(NPN2,NPN2),TRGQI(NPN2,NPN2)
cc      REAL*8 TR1(NPN2,NPN2),TI1(NPN2,NPN2)
*________
      COMMON /TMAT99/ 
     &            R11,R12,R21,R22,I11,I12,I21,I22,RG11,RG12,RG21,RG22,
     &            IG11,IG12,IG21,IG22          !only between TMATR routines
      COMMON /CBESS/ J,Y,JR,JI,DJ,DY,DJR,DJI
cc      COMMON /CT/ TR1,TI1                      !output from TT routine
      COMMON /CTT/ QR,QI,RGQR,RGQI             !input for TT routine
*________
      MM1=M
      QM=DFLOAT(M)
      QMM=QM*QM
      NG=2*NGAUSS
      NGSS=NG
      FACTOR=1D0
*      
      IF (NCHECK.EQ.1) THEN          !THETA=PI/2 is mirror symmetry plane
            NGSS=NGAUSS
            FACTOR=2D0
         ELSE                        !THETA=PI/2 is not a mirror symmetry plane
            CONTINUE
      ENDIF
*
      SI=1D0
      NM=NMAX+NMAX
      
      DO 5 N=1,NM
           SI=-SI
           SIG(N)=SI              !=(-1)**N
    5 CONTINUE
*
* Assigning Wigner d-matrices - assuming mirror symmetry 
* in the \theta=\pi/2 plane:
    
      DO 25 I=1,NGAUSS
      
         I1=NGAUSS-I+1
         I2=NGAUSS+I
*
         CALL VIG (X(I1),NMAX,M,DV1,DV2)
*
         DO N=1,NMAX
         
            DD1=DV1(N)
            DD2=DV2(N)
            D1(I1,N)=DD1
            D2(I1,N)=DD2

         IF (NAXSM.EQ.1) THEN         !Gauss abscissas chosen +/- symmetric

* using (4.2.4) and (4.2.6) of {Ed},  
*           d_{0m}^{(l)}(\pi-\theta) = (-1)^{l+m} d_{0m}^{(l)}(\theta)

            SI=SIG(N+M)                  !=(-1)**(N+M)
                                         !exactly what follows from {Ed}           
            D1(I2,N)=DD1*SI       
            D2(I2,N)=-DD2*SI
           
         END IF            
         ENDDO 
*         
         IF (NAXSM.EQ.0) THEN        !Gauss abscissas not chosen +/- symmetric
*
         CALL VIG ( X(I2), NMAX, M, DV1, DV2)
*
          DO N=1,NMAX               
            DD1=DV1(N)
            DD2=DV2(N)
            D1(I2,N)=DD1
            D2(I2,N)=DD2
          ENDDO 
           
         END IF           

   25 CONTINUE
*
*  Assigning r^2(\theta)*weight product:
   
      DO 40 I=1,NGSS
           WR=W(I)*R(I)
           
cc           if (dr(i).eq.0.d0) WR=0.d0   !temporarily only
           
           DS(I)=S(I)*QM*WR       !=DFLOAT(M)*W(I)*r^2(\theta)/(|\sin\theta|)
           DSS(I)=SS(I)*QMM       !=DFLOAT(M)**2/(\sin^2\theta)
           RR(I)=WR
                      
   40 CONTINUE
* 
      DO 300  N1=MM1,NMAX
           AN1=AN(N1)
           
           DO 300 N2=MM1,NMAX
                AN2=AN(N2)
                AR11=0D0
                AR12=0D0
                AR21=0D0
                AR22=0D0
                AI11=0D0
                AI12=0D0
                AI21=0D0
                AI22=0D0
                GR11=0D0
                GR12=0D0
                GR21=0D0
                GR22=0D0
                GI11=0D0
                GI12=0D0
                GI21=0D0
                GI22=0D0
                SI=SIG(N1+N2)
 
                DO 200 I=1,NGSS
                    D1N1=D1(I,N1)
                    D2N1=D2(I,N1)
                    D1N2=D1(I,N2)
                    D2N2=D2(I,N2)
                    A11=D1N1*D1N2
                    A12=D1N1*D2N2
                    A21=D2N1*D1N2
                    A22=D2N1*D2N2
                    AA1=A12+A21            != D1N1*D2N2+D2N1*D1N2
                    AA2=A11*DSS(I)+A22     !=(D1N1*D1N2)*DFLOAT(M)**2/(\sin^2\theta)
                                           ! +D2N1*D2N2
                    
* Vector spherical harmonics:
C  Since refractive index is allowed to be complex in general,
C  the Bessel function j_l(k_in*r) is complex. The code below 
C  performs a separation of the complex integrand in Waterman's
C  surface integral into its respective real and imaginary 
C  parts:
                     
                    QJ1=J(I,N1)
                    QY1=Y(I,N1)
                    QJR2=JR(I,N2)
                    QJI2=JI(I,N2)
                    QDJR2=DJR(I,N2)
                    QDJI2=DJI(I,N2)
                    QDJ1=DJ(I,N1)
                    QDY1=DY(I,N1)

* Re and Im of j_{n2}(k_{in}r) j_{n1}(k_{out}r): 

                    C1R=QJR2*QJ1
                    C1I=QJI2*QJ1
                    
* Re and Im of j_{n2}(k_{in}r) h_{n1}(k_{out}r):
                    
                    B1R=C1R-QJI2*QY1
                    B1I=C1I+QJR2*QY1
 
* Re and Im of j_{n2}(k_{in}r) j_{n1}'(k_{out}r): 

                    C2R=QJR2*QDJ1
                    C2I=QJI2*QDJ1

* Re and Im of j_{n2}(k_{in}r) h_{n1}'(k_{out}r):                    
                    
                    B2R=C2R-QJI2*QDY1
                    B2I=C2I+QJR2*QDY1                    
 
                    DDRI=DDR(I)               !1/(k_{out}r)

* Re and Im of [1/(k_{out}r)]*j_{n2}(k_{in}r) j_{n1}(k_{out}r) 
    
                    C3R=DDRI*C1R
                    C3I=DDRI*C1I
                    
* Re and Im of [1/(k_{out}r)]*j_{n2}(k_{in}r) h_{n1}(k_{out}r):
                    
                    B3R=DDRI*B1R
                    B3I=DDRI*B1I

* Re and Im of j_{n2}'(k_{in}r) j_{n1}(k_{out}r): 
 
                    C4R=QDJR2*QJ1
                    C4I=QDJI2*QJ1
                                        
* Re and Im of j_{n2}'(k_{in}r) h_{n1}(k_{out}r): 
                    
                    B4R=C4R-QDJI2*QY1
                    B4I=C4I+QDJR2*QY1
                    
 
                    DRRI=DRR(I)               !Re[1/(k_{in}r)]
                    DRII=DRI(I)               !Im[1/(k_{in}r)]
                    
* Re and Im of [1/(k_{in}r)] j_{n2}(k_{in}r) j_{n1}(k_{out}r):  
                  
                    C5R=C1R*DRRI-C1I*DRII
                    C5I=C1I*DRRI+C1R*DRII
                    
* Re and Im of [1/(k_{in}r)] j_{n2}(k_{in}r) h_{n1}(k_{out}r):
                    
                    B5R=B1R*DRRI-B1I*DRII
                    B5I=B1I*DRRI+B1R*DRII
                    
                    
* Re and Im of j_{n2}'(k_{in}r) j_{n1}'(k_{out}r):  

                    C6R=QDJR2*QDJ1
                    C6I=QDJI2*QDJ1
                    
* Re and Im of j_{n2}'(k_{in}r) h_{n1}'(k_{out}r):

                    B6R=C6R-QDJI2*QDY1
                    B6I=C6I+QDJR2*QDY1
                    
* Re and Im of [1/(k_{out}r)] j_{n2}'(k_{in}r) j_{n1}(k_{out}r): 
 
                    C7R=C4R*DDRI
                    C7I=C4I*DDRI
                    
* Re and Im of [1/(k_{out}r)] j_{n2}'(k_{in}r) h_{n1}(k_{out}r): 
                    
                    B7R=B4R*DDRI
                    B7I=B4I*DDRI

* Re and Im of [1/(k_{in}r)] j_{n2}(k_{in}r) j_{n1}'(k_{out}r):  

                    C8R=C2R*DRRI-C2I*DRII
                    C8I=C2I*DRRI+C2R*DRII

* Re and Im of [1/(k_{in}r)] j_{n2}(k_{in}r) h_{n1}'(k_{out}r): 
                    
                    B8R=B2R*DRRI-B2I*DRII
                    B8I=B2I*DRRI+B2R*DRII


* %%%%%%%%%  Forming integrands of J-matrices (J^{11}=J^{22}=0 for m=0):
 
                    URI=DR(I)
                    DSI=DS(I)
                    RRI=RR(I)
 
                    IF (NCHECK.EQ.1.AND.SI.GT.0D0) GO TO 150

* [DFLOAT(M)*W(I)*r^2(I)/(|\sin\theta|)]*(D1N1*D2N2+D2N1*D1N2):
                    E1=DSI*AA1

                    AR11=AR11+E1*B1R
                    AI11=AI11+E1*B1I
                    GR11=GR11+E1*C1R
                    GI11=GI11+E1*C1I
                    
                    IF (NCHECK.EQ.1) GO TO 160
 
  150               CONTINUE

                    
* w(i)*r^2(\theta)*[(D1N1*D1N2)*DFLOAT(M)**2/(\sin^2\theta)+D2N1*D2N2]:
                    F1=RRI*AA2            !prefactor containing r^2(\theta)<->hat{r} part
                    
* N1*(N1+1)*w(i)*r(\theta)*[dr/(d\theta)]*D1N1*D2N2:                     
                    F2=RRI*URI*AN1*A12     !prefactor containing r(\theta)*[dr/(d\theta)]
                                           !hat{theta} part
                    
                    AR12=AR12+F1*B2R+F2*B3R        !~Re J^{12}
                    AI12=AI12+F1*B2I+F2*B3I        !~Im J^{12}
                    
                    GR12=GR12+F1*C2R+F2*C3R        !~Re Rg J^{12}
                    GI12=GI12+F1*C2I+F2*C3I        !~Im Rg J^{12}

* N2*(N2+1)*w(i)*r(\theta)*[dr/(d\theta)]*D2N1*D1N2:   
                    F2=RRI*URI*AN2*A21     !prefactor containing r(\theta)*[dr/(d\theta)]
                                           !hat{theta} part
                    
                    AR21=AR21+F1*B4R+F2*B5R
                    AI21=AI21+F1*B4I+F2*B5I
                    
                    GR21=GR21+F1*C4R+F2*C5R
                    GI21=GI21+F1*C4I+F2*C5I
                    
                    IF (NCHECK.EQ.1) GO TO 200
 
  160               E2=DSI*URI*A11
                    E3=E2*AN2
                    E2=E2*AN1
                    
                    AR22=AR22+E1*B6R+E2*B7R+E3*B8R
                    AI22=AI22+E1*B6I+E2*B7I+E3*B8I
                    
                    GR22=GR22+E1*C6R+E2*C7R+E3*C8R
                    GI22=GI22+E1*C6I+E2*C7I+E3*C8I
                    
  200           CONTINUE
  
                AN12=ANN(N1,N2)*FACTOR
                
                R11(N1,N2)=AR11*AN12
                R12(N1,N2)=AR12*AN12
                R21(N1,N2)=AR21*AN12
                R22(N1,N2)=AR22*AN12
                I11(N1,N2)=AI11*AN12
                I12(N1,N2)=AI12*AN12
                I21(N1,N2)=AI21*AN12
                I22(N1,N2)=AI22*AN12
                
                RG11(N1,N2)=GR11*AN12
                RG12(N1,N2)=GR12*AN12
                RG21(N1,N2)=GR21*AN12
                RG22(N1,N2)=GR22*AN12
                IG11(N1,N2)=GI11*AN12
                IG12(N1,N2)=GI12*AN12
                IG21(N1,N2)=GI21*AN12
                IG22(N1,N2)=GI22*AN12
 
  300 CONTINUE
  
      TPIR=PIR                 !Re [1/k_{in}^2]
      TPII=PII                 !Im [1/k_{in}^2]
      TPPI=PPI                 !1/k_{out}^2   
      
      NM=NMAX-MM1+1
      DO 310 N1=MM1,NMAX
           K1=N1-MM1+1
           KK1=K1+NM
           
           DO 310 N2=MM1,NMAX
                K2=N2-MM1+1
                KK2=K2+NM
 
                TAR11=-R11(N1,N2)
                TAI11=-I11(N1,N2)
                TGR11=-RG11(N1,N2)
                TGI11=-IG11(N1,N2)
 
                TAR12= I12(N1,N2)
                TAI12=-R12(N1,N2)
                TGR12= IG12(N1,N2)
                TGI12=-RG12(N1,N2)
 
                TAR21=-I21(N1,N2)
                TAI21= R21(N1,N2)
                TGR21=-IG21(N1,N2)
                TGI21= RG21(N1,N2)
 
                TAR22=-R22(N1,N2)
                TAI22=-I22(N1,N2)
                TGR22=-RG22(N1,N2)
                TGI22=-IG22(N1,N2)
 
                TQR(K1,K2)=TPIR*TAR21-TPII*TAI21+TPPI*TAR12
                TQI(K1,K2)=TPIR*TAI21+TPII*TAR21+TPPI*TAI12
                TRGQR(K1,K2)=TPIR*TGR21-TPII*TGI21+TPPI*TGR12
                TRGQI(K1,K2)=TPIR*TGI21+TPII*TGR21+TPPI*TGI12
 
                TQR(K1,KK2)=TPIR*TAR11-TPII*TAI11+TPPI*TAR22
                TQI(K1,KK2)=TPIR*TAI11+TPII*TAR11+TPPI*TAI22
                TRGQR(K1,KK2)=TPIR*TGR11-TPII*TGI11+TPPI*TGR22
                TRGQI(K1,KK2)=TPIR*TGI11+TPII*TGR11+TPPI*TGI22
 
                TQR(KK1,K2)=TPIR*TAR22-TPII*TAI22+TPPI*TAR11
                TQI(KK1,K2)=TPIR*TAI22+TPII*TAR22+TPPI*TAI11
                TRGQR(KK1,K2)=TPIR*TGR22-TPII*TGI22+TPPI*TGR11
                TRGQI(KK1,K2)=TPIR*TGI22+TPII*TGR22+TPPI*TGI11
 
                TQR(KK1,KK2)=TPIR*TAR12-TPII*TAI12+TPPI*TAR21
                TQI(KK1,KK2)=TPIR*TAI12+TPII*TAR12+TPPI*TAI21
                TRGQR(KK1,KK2)=TPIR*TGR12-TPII*TGI12+TPPI*TGR21
                TRGQI(KK1,KK2)=TPIR*TGI12+TPII*TGR12+TPPI*TGI21
                
  310 CONTINUE
 
      NNMAX=2*NM
      DO 320 N1=1,NNMAX
           DO 320 N2=1,NNMAX
                QR(N1,N2)=TQR(N1,N2)
                QI(N1,N2)=TQI(N1,N2)
                RGQR(N1,N2)=TRGQR(N1,N2)
                RGQI(N1,N2)=TRGQI(N1,N2)
  320 CONTINUE
* 
      CALL TT(NM,NCHECK)
* 
      RETURN
      END
 
C*****************************************************************
 
      SUBROUTINE VIG (X, NMAX, M, DV1, DV2)
C--------/---------/---------/---------/---------/---------/---------/--
C >>> X,NMAX,M (only nonnegative)
C <<< DV1, DV2
C =============
C     For a given azimuthal number M.GE.0 returns 
C      the Wigner d-functions , i.e.,
C
C     DV1(N)=dvig(0,m,n,arccos x) = d_{0m}^{(l)} 
C     and
C     DV2(N)=[d/d(arccos x)] dvig(0,m,n,arccos x)  ! = d d_{0m}^{(l)}/d\theta
C
C     for 1.LE.N.LE.NMAX and 0.LE.X.LE.1
C     (For a given M.NEQ.0, only the M.LE.N.LE.NMAX terms are determined!)
C     According to Eq. (4.1.24) of Ref. \ct{Ed}:
C
C             d_{00}^{(l)}(\theta)= P_l(\cos\theta) ===>
C
C     (Rodrigues formula [Eq. (2.5.14) of Ref. \ct{Ed}] then yields 
C                       P_1(x)=x; P_2=(3x^2-1)/2; etc.
C     One can show that $d_{00}^{(1)}(\theta)=\cos\theta$
C
C     Similar to routine VIGAMPL, which however returns the Wigner d-functions 
C     divided by sin\theta, i.e.,
C     DV1(N)=dvig(0,m,n,arccos x)/sin(arccos x) = d_{0m}^{(l)}/sin\theta
C
C     Made using recurrences of  Ref. \ct{Mis39}  
C     (There is a missing $l$ factor in the 2nd term in the curly bracket 
C     in recurrence (35) of Ref. \ct{Mis39} for DV2).  
C
C     One has (see Eq. (4.2.5) of \ct{Ed}):
C                       $d_{0m}^{(l)}=(-1)^m d_{0-m}^{(l)}$
C     and (see Eq. (35) of \ct{Mis91}):
C            $dd_{0m}^{(l)}/(d\theta)=(-1)^m dd_{0-m}^{(l)}/(d\theta)$
C
C     X=cos(theta), where theta is the polar angle
C     NMAX - angular momentum cutoff
C
C     CALLED BY TMATR AND TMATR0 routines
C--------/---------/---------/---------/---------/---------/---------/--
      INCLUDE 'ampld.par.f'
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 DV1(NPN1),DV2(NPN1)
 
      A=1D0
      QS=DSQRT(1D0-X*X)
      QS1=1D0/QS
      DO N=1,NMAX
         DV1(N)=0D0
         DV2(N)=0D0
      ENDDO   
      
      IF (M.NE.0) GO TO 20
      
      D1=1D0
      D2=X  
      DO N=1,NMAX
         QN=DFLOAT(N)
         QN1=DFLOAT(N+1)
         QN2=DFLOAT(2*N+1)
         D3=(QN2*X*D2-QN*D1)/QN1          !recurrence (31) of Ref. {Mis39}
         DER=QS1*(QN1*QN/QN2)*(-D1+D3)    !recurrence (35) of Ref. {Mis39}
         DV1(N)=D2
         DV2(N)=DER
         D1=D2
         D2=D3
      ENDDO   
      RETURN
      
   20 QMM=DFLOAT(M*M)
   
*A_m initialization - recurrence (34) of Ref. {Mis39}
      DO I=1,M
         I2=I*2
         A=A*DSQRT(DFLOAT(I2-1)/DFLOAT(I2))*QS  
      ENDDO 
*  
      D1=0D0
      D2=A 

      DO N=M,NMAX
         QN=DFLOAT(N)
         QN2=DFLOAT(2*N+1)
         QN1=DFLOAT(N+1)
         QNM=DSQRT(QN*QN-QMM)
         QNM1=DSQRT(QN1*QN1-QMM)
         D3=(QN2*X*D2-QNM*D1)/QNM1              !recurrence (31) of Ref. {Mis39}
         DER=QS1*(-QN1*QNM*D1+QN*QNM1*D3)/QN2   !recurrence (35) of Ref. {Mis39}
         DV1(N)=D2
         DV2(N)=DER
         D1=D2
         D2=D3
      ENDDO 
*  
      RETURN
      END 
 
C**********************************************************************
 
      SUBROUTINE TT(NMAX,NCHECK)
C--------/---------/---------/---------/---------/---------/---------/--
C >>> NMAX,NCHECK
C <<< COMMON BLOCKS
C=================
C  NMAX - angular momentum cutoff
C  NCHECK -
C                                                                     
C   CALCULATION OF THE MATRIX    T = - RG(Q) * (Q**(-1))              
C                                                                     
C   INPUT INFORMATION IS IN COMMON /CTT/                             
C   OUTPUT INFORMATION IS IN COMMON /CT/                              
C                                                                     
C--------/---------/---------/---------/---------/---------/---------/--
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER NOUT

* number of the output unit
      PARAMETER (NOUT=35)
      INCLUDE 'ampld.par.f'

      REAL*8  QR(NPN2,NPN2),QI(NPN2,NPN2),EMACH,
     *       RGQR(NPN2,NPN2),RGQI(NPN2,NPN2)
cc      REAL*8 F(NPN2,NPN2),B(NPN2),WORK(NPN2),
cc     *       A(NPN2,NPN2),C(NPN2,NPN2),D(NPN2,NPN2),E(NPN2,NPN2)
      REAL*8 TR1(NPN2,NPN2),TI1(NPN2,NPN2)
      COMPLEX*16 ZQ(NPN2,NPN2),ZX(NPN2),ZW(NPN2)
      INTEGER IPIV(NPN2),IPVT(NPN2)
*
      COMMON /CHOICE/ ICHOICE
      COMMON /CT/ TR1,TI1
      COMMON /CTT/ QR,QI,RGQR,RGQI
*
      DATA EMACH/1.D-10/
*
      NNMAX=2*NMAX
  
      DO I=1,NNMAX
       DO J=1,NNMAX
          ZQ(I,J)=DCMPLX(QR(I,J),QI(I,J))
       ENDDO
      ENDDO

      IF (ICHOICE.EQ.2) GOTO 5    ! NAG or not NAG decision tree

********************************************************************
*   Inversion from NAG-LIB or Waterman's method    !NAG library used
*
       INFO=0
*
      write(6,*)'NAG in TMTAXSP has been deactivated'
	write(6,*)   '(see lines 2805-2808 therein)'
	stop
ct           CALL F07ARF(NNMAX,NNMAX,ZQ,NPN2,IPIV,INFO)
ct           IF (INFO.NE.0) WRITE(NOUT,1100) INFO
ct           CALL F07AWF(NNMAX,ZQ,NPN2,IPIV,ZW,NPN2,INFO)
ct           IF (INFO.NE.0) WRITE(NOUT,1100) INFO
*
 1100      FORMAT ('WARNING:  info=', i2)

* Calculate T-matrix = - RG(Q) * (Q**(-1))
*
       DO I=1,NNMAX
          DO J=1,NNMAX
             TR=0D0
             TI=0D0
             DO K=1,NNMAX
                    ARR=RGQR(I,K)
                    ARI=RGQI(I,K)
                    AR=ZQ(K,J)
                    AI=DIMAG(ZQ(K,J))
                    TR=TR-ARR*AR+ARI*AI
                    TI=TI-ARR*AI-ARI*AR
                 ENDDO
             TR1(I,J)=TR
             TI1(I,J)=TI
          ENDDO
       ENDDO
       GOTO 70                        !Return

*********************************************************************
 
C  Gaussian elimination             !NAG library not used

  5   CALL ZGER(ZQ,IPIV,NNMAX,NPN2,EMACH)  !Gauss elimination of ZQ to
                                           !a lower diagonal matrix
      DO 6 I=1,NNMAX
              DO K=1,NNMAX    !Initialization of the right-hand side ZB
                              !(a row vector) of the matrix equation ZX*ZQ=ZB
 
              ZX(K)=DCMPLX(RGQR(I,K),RGQI(I,K))
              ENDDO 

      CALL ZSUR (ZQ,IPIV,ZX,NNMAX,NPN2,EMACH)  !Solving ZX*ZQ=ZB by 
                                               !backsubstition
                                               !(ZX overwritten on exit)
             DO K=1,NNMAX
*
* Assign T-matrix elements = - RG(Q) * (Q**(-1))
*
             TR1(I,K)=-DBLE(ZX(K))
             TI1(I,K)=-DIMAG(ZX(K))
             ENDDO
  6   CONTINUE
* 
   70 RETURN
      END
 
C********************************************************************
 
      SUBROUTINE PROD(A,B,C,NDIM,N)
C--------/---------/---------/---------/---------/---------/---------/--
C >>> A,B,NDIM,N
C <<< C=A*B
C=================
C--------/---------/---------/---------/---------/---------/---------/--
      REAL*8 A(NDIM,N),B(NDIM,N),C(NDIM,N),cij
*
      DO 10 I=1,N
           DO 10 J=1,N
                CIJ=0d0
                DO 5 K=1,N
                     CIJ=CIJ+A(I,K)*B(K,J)
    5           CONTINUE
                C(I,J)=CIJ
   10 CONTINUE
*
      RETURN
      END
 
C**********************************************************************
 
      SUBROUTINE INV1 (NMAX,F,A)
C--------/---------/---------/---------/---------/---------/---------/--
C >>> EPS
C <<< RAT
C=================
C  NMAX - angular momentum cutoff
C--------/---------/---------/---------/---------/---------/---------/--
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'ampld.par.f'
      REAL*8  A(NPN2,NPN2),F(NPN2,NPN2),B(NPN1),
     *        WORK(NPN1),Q1(NPN1,NPN1),Q2(NPN1,NPN1),
     &        P1(NPN1,NPN1),P2(NPN1,NPN1)
      INTEGER IPVT(NPN1),IND1(NPN1),IND2(NPN1)
*
      NDIM=NPN1
      NN1=(DFLOAT(NMAX)-0.1D0)*0.5D0+1D0 
      NN2=NMAX-NN1
*
      DO 5 I=1,NMAX
         IND1(I)=2*I-1
         IF(I.GT.NN1) IND1(I)=NMAX+2*(I-NN1)
         IND2(I)=2*I
         IF(I.GT.NN2) IND2(I)=NMAX+2*(I-NN2)-1
    5 CONTINUE
      NNMAX=2*NMAX
*
      DO 15 I=1,NMAX
         I1=IND1(I)
         I2=IND2(I)
         DO 15 J=1,NMAX
            J1=IND1(J)
            J2=IND2(J)
            Q1(J,I)=F(J1,I1)
            Q2(J,I)=F(J2,I2)
  15  CONTINUE
*
      CALL INVERT(NDIM,NMAX,Q1,P1,COND,IPVT,WORK,B)
      CALL INVERT(NDIM,NMAX,Q2,P2,COND,IPVT,WORK,B)
*
      DO 30 I=1,NNMAX
         DO 30 J=1,NNMAX
            A(J,I)=0D0
  30  CONTINUE
      DO 40 I=1,NMAX
         I1=IND1(I)
         I2=IND2(I)
         DO 40 J=1,NMAX
            J1=IND1(J)
            J2=IND2(J)
            A(J1,I1)=P1(J,I)
            A(J2,I2)=P2(J,I)
  40  CONTINUE
*
      RETURN
      END
 
C*********************************************************************
 
      SUBROUTINE INVERT (NDIM,N,A,X,COND,IPVT,WORK,B)
C--------/---------/---------/---------/---------/---------/---------/--
C >>> EPS
C <<< RAT
C=================
C--------/---------/---------/---------/---------/---------/---------/--
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 A(NDIM,N),X(NDIM,N),WORK(N),B(N)
      INTEGER IPVT(N)
*
      CALL DECOMP (NDIM,N,A,COND,IPVT,WORK)
*
      IF (COND+1D0.EQ.COND) PRINT 5,COND
C     IF (COND+1D0.EQ.COND) STOP
   5  FORMAT(' THE MATRIX IS SINGULAR FOR THE GIVEN NUMERICAL ACCURACY '
     *      ,'COND = ',D12.6)

      DO 30 I=1,N
           DO 10 J=1,N
                B(J)=0D0
                IF (J.EQ.I) B(J)=1D0
  10       CONTINUE
*
           CALL SOLVE (NDIM,N,A,B,IPVT)
*
           DO 30 J=1,N
                X(J,I)=B(J)
  30  CONTINUE
*
      RETURN
      END
 
C********************************************************************
 
      SUBROUTINE DECOMP (NDIM,N,A,COND,IPVT,WORK)
C--------/---------/---------/---------/---------/---------/---------/--
C >>> EPS
C <<< RAT
C=================
C--------/---------/---------/---------/---------/---------/---------/--
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 A(NDIM,N),COND,WORK(N)
      INTEGER IPVT(N)
*
      IPVT(N)=1
      IF(N.EQ.1) GO TO 80
      NM1=N-1
      ANORM=0D0
      DO 10 J=1,N
          T=0D0
          DO 5 I=1,N
              T=T+DABS(A(I,J))
    5     CONTINUE
          IF (T.GT.ANORM) ANORM=T
   10 CONTINUE
      DO 35 K=1,NM1
          KP1=K+1
          M=K
          DO 15 I=KP1,N
              IF (DABS(A(I,K)).GT.DABS(A(M,K))) M=I
   15     CONTINUE
          IPVT(K)=M
          IF (M.NE.K) IPVT(N)=-IPVT(N)
          T=A(M,K)
          A(M,K)=A(K,K)
          A(K,K)=T
          IF (T.EQ.0d0) GO TO 35
          DO 20 I=KP1,N
              A(I,K)=-A(I,K)/T
   20     CONTINUE
          DO 30 J=KP1,N
              T=A(M,J)
              A(M,J)=A(K,J)
              A(K,J)=T
              IF (T.EQ.0D0) GO TO 30
              DO 25 I=KP1,N
                  A(I,J)=A(I,J)+A(I,K)*T
   25         CONTINUE
   30     CONTINUE
   35 CONTINUE
      DO 50 K=1,N
          T=0D0
          IF (K.EQ.1) GO TO 45
          KM1=K-1
          DO 40 I=1,KM1
              T=T+A(I,K)*WORK(I)
   40     CONTINUE
   45     EK=1D0
          IF (T.LT.0D0) EK=-1D0
          IF (A(K,K).EQ.0D0) GO TO 90
          WORK(K)=-(EK+T)/A(K,K)
   50 CONTINUE
      DO 60 KB=1,NM1
          K=N-KB
          T=0D0
          KP1=K+1
          DO 55 I=KP1,N
              T=T+A(I,K)*WORK(K)
   55     CONTINUE
          WORK(K)=T
          M=IPVT(K)
          IF (M.EQ.K) GO TO 60
          T=WORK(M)
          WORK(M)=WORK(K)
          WORK(K)=T
   60 CONTINUE
      YNORM=0D0
      DO 65 I=1,N
          YNORM=YNORM+DABS(WORK(I))
   65 CONTINUE
*
      CALL SOLVE (NDIM,N,A,WORK,IPVT)
*
      ZNORM=0D0
      DO 70 I=1,N
          ZNORM=ZNORM+DABS(WORK(I))
   70 CONTINUE
      COND=ANORM*ZNORM/YNORM
      IF (COND.LT.1d0) COND=1D0
      RETURN
   80 COND=1D0
      IF (A(1,1).NE.0D0) RETURN
   90 COND=1D52
*
      RETURN
      END
 
C**********************************************************************
 
      SUBROUTINE SOLVE (NDIM,N,A,B,IPVT)
C--------/---------/---------/---------/---------/---------/---------/--
C >>> EPS
C <<< RAT
C=================
C--------/---------/---------/---------/---------/---------/---------/--
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 A(NDIM,N),B(N)
      INTEGER IPVT(N)
*
      IF (N.EQ.1) GO TO 50
      NM1=N-1
      DO 20 K=1,NM1
          KP1=K+1
          M=IPVT(K)
          T=B(M)
          B(M)=B(K)
          B(K)=T
          DO 10 I=KP1,N
              B(I)=B(I)+A(I,K)*T
   10     CONTINUE
   20 CONTINUE
      DO 40 KB=1,NM1
          KM1=N-KB
          K=KM1+1
          B(K)=B(K)/A(K,K)
          T=-B(K)
          DO 30 I=1,KM1
              B(I)=B(I)+A(I,K)*T
   30     CONTINUE
   40 CONTINUE
   50 B(1)=B(1)/A(1,1)
*
      RETURN
      END
 
C*****************************************************************
 
      SUBROUTINE SAREA (D,RAT)
C--------/---------/---------/---------/---------/---------/---------/--
C >>> EPS
C <<< RAT
C=================
C--------/---------/---------/---------/---------/---------/---------/--
      IMPLICIT REAL*8 (A-H,O-Z)
      IF (D.GE.1) GO TO 10
      E=DSQRT(1D0-D*D)
      R=0.5D0*(D**(2D0/3D0) + D**(-1D0/3D0)*DASIN(E)/E)
      R=DSQRT(R)
      RAT=1D0/R
      RETURN
   10 E=DSQRT(1D0-1D0/(D*D))
      R=0.25D0*(2D0*D**(2D0/3D0) + D**(-4D0/3D0)*DLOG((1D0+E)/(1D0-E))
     &   /E)
      R=DSQRT(R)
      RAT=1D0/R
*
      return
      END
 
c****************************************************************
 
      SUBROUTINE SURFCH (N,E,RAT)
C--------/---------/---------/---------/---------/---------/---------/--
C >>> N,E,RAT
C <<< RAT
C=================
C--------/---------/---------/---------/---------/---------/---------/--
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 X(60),W(60)
*
      DN=DFLOAT(N)
      EN=E*DN
      NG=60
*
* GIF division points and weights
*
      CALL GAUSS (NG,0,0,X,W)
*
      S=0D0
      V=0D0
      DO 10 I=1,NG
         XI=X(I)
         DX=DACOS(XI)
         DXN=DN*DX
         DS=DSIN(DX)
         DSN=DSIN(DXN)
         DCN=DCOS(DXN)
         A=1D0+E*DCN
         A2=A*A
         ENS=EN*DSN
         S=S+W(I)*A*DSQRT(A2+ENS*ENS)
         V=V+W(I)*(DS*A+XI*ENS)*DS*A2
   10 CONTINUE
      RS=DSQRT(S*0.5D0)
      RV=(V*3D0/4D0)**(1D0/3D0)
      RAT=RV/RS
*
      RETURN
      END
 
C********************************************************************
 
      SUBROUTINE SAREAC (EPS,RAT)
C--------/---------/---------/---------/---------/---------/---------/--
C >>> EPS
C <<< RAT
C=================
C--------/---------/---------/---------/---------/---------/---------/--
      IMPLICIT REAL*8 (A-H,O-Z)
*
      RAT=(1.5D0/EPS)**(1D0/3D0)
      RAT=RAT/DSQRT( (EPS+2D0)/(2D0*EPS) )
*
      RETURN
      END

C**********************************************************************

      SUBROUTINE DROP (RAT)
C--------/---------/---------/---------/---------/---------/---------/--
C >>> EPS
C <<< RAT
C=================
C--------/---------/---------/---------/---------/---------/---------/--      
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER NOUT
* number of the output unit
      PARAMETER (NOUT=35)  
      PARAMETER (NC=10, NG=60) 
                 
      REAL*8 X(NG),W(NG),C(0:NC)
      COMMON /CDROP/ C,R0V
      C(0)=-0.0481 D0
      C(1)= 0.0359 D0
      C(2)=-0.1263 D0
      C(3)= 0.0244 D0
      C(4)= 0.0091 D0
      C(5)=-0.0099 D0
      C(6)= 0.0015 D0
      C(7)= 0.0025 D0
      C(8)=-0.0016 D0
      C(9)=-0.0002 D0
      C(10)= 0.0010 D0
*
* GIF division points and weights
*
      CALL GAUSS (NG,0,0,X,W)
*
      S=0D0
      V=0D0
      DO I=1,NG
         XI=DACOS(X(I))
         WI=W(I)
         RI=1D0+C(0)
         DRI=0D0
         DO N=1,NC
            XIN=XI*N
            RI=RI+C(N)*DCOS(XIN)
            DRI=DRI-C(N)*N*DSIN(XIN)
         ENDDO
         SI=DSIN(XI)
         CI=X(I)
         RISI=RI*SI
         S=S+WI*RI*DSQRT(RI*RI+DRI*DRI)
         V=V+WI*RI*RISI*(RISI-DRI*CI)
      ENDDO
      RS=DSQRT(S*0.5D0)
      RV=(V*3D0*0.25D0)**(1D0/3D0)
      IF (DABS(RAT-1D0).GT.1D-8) RAT=RV/RS
      R0V=1D0/RV
      WRITE(NOUT,1000) R0V
      DO N=0,NC
         WRITE(NOUT,1001) N,C(N)
      ENDDO
 1000 FORMAT ('r_0/r_ev=',F7.4)
 1001 FORMAT ('c_',I2,'=',F7.4)
 
      RETURN
      END

C********************************************************************
 
      SUBROUTINE GAUSS (N,IND1,IND2,Z,W)
C--------/---------/---------/---------/---------/---------/---------/--
C >>> N,IND1,IND2
C <<< Z,W
C=================
C    CALCULATION OF POINTS AND WEIGHTS OF GAUSSIAN QUADRATURE         
C    FORMULA. IF IND1 = 0 - ON INTERVAL (-1,1), IF IND1 = 1 - ON      
C    INTERVAL  (0,1). IF  IND2 = 1 RESULTS ARE PRINTED. 
C              
C    N - NUMBER OF GIF DIVISION POINTS (mostly N=NGAUSS in main program)                                        
C    Z - DIVISION POINTS                                              
C    W - WEIGHTS                                                      
C--------/---------/---------/---------/---------/---------/---------/--
      IMPLICIT REAL*8 (A-H,P-Z)
      REAL*8 Z(N),W(N)
      DATA A,B,C /1D0,2D0,3D0/
      IND=MOD(N,2)
      K=N/2+IND
      F=DFLOAT(N)
      DO 100 I=1,K
          M=N+1-I
          IF(I.EQ.1) X=A-B/((F+A)*F)
          IF(I.EQ.2) X=(Z(N)-A)*4D0+Z(N)
          IF(I.EQ.3) X=(Z(N-1)-Z(N))*1.6D0+Z(N-1)
          IF(I.GT.3) X=(Z(M+1)-Z(M+2))*C+Z(M+3)
          IF(I.EQ.K.AND.IND.EQ.1) X=0D0
          NITER=0
          CHECK=1D-16
   10     PB=1D0
          NITER=NITER+1
          IF (NITER.LE.100) GO TO 15
          CHECK=CHECK*10D0
   15     PC=X
          DJ=A
          DO 20 J=2,N
              DJ=DJ+A
              PA=PB
              PB=PC
   20         PC=X*PB+(X*PB-PA)*(DJ-A)/DJ
          PA=A/((PB-X*PC)*F)
          PB=PA*PC*(A-X*X)
          X=X-PB
          IF(DABS(PB).GT.CHECK*DABS(X)) GO TO 10
          Z(M)=X
          W(M)=PA*PA*(A-X*X)
          IF(IND1.EQ.0) W(M)=B*W(M)
          IF(I.EQ.K.AND.IND.EQ.1) GO TO 100
          Z(I)=-Z(M)
          W(I)=W(M)
  100 CONTINUE
      IF(IND2.NE.1) GO TO 110
      PRINT 1100,N
 1100 FORMAT(' ***  POINTS AND WEIGHTS OF GAUSSIAN QUADRATURE FORMULA',
     * ' OF ',I4,'-TH ORDER')
      DO 105 I=1,K
          ZZ=-Z(I)
  105     PRINT 1200,I,ZZ,I,W(I)
 1200 FORMAT(' ',4X,'X(',I4,') = ',F17.14,5X,'W(',I4,') = ',F17.14)
      GO TO 115
  110 CONTINUE
C     PRINT 1300,N
C 1300 FORMAT(' GAUSSIAN QUADRATURE FORMULA OF ',I4,'-TH ORDER IS USED')
  115 CONTINUE
      IF(IND1.EQ.0) GO TO 140
      DO 120 I=1,N
  120     Z(I)=(A+Z(I))/B
  140 CONTINUE
  
      RETURN
      END

C********************************************************************

      SUBROUTINE gauleg(x1,x2,x,w,n)
C--------/---------/---------/---------/---------/---------/---------/--
C  Given the lower and upper limits of integration x1 and x2, and given n
C  this routine returns arrays x(1:n) and w(1:n) of length n, containing
C  the abscissas and weights of the Gaussian-Legendre n-point quadrature 
C  formula.
C--------/---------/---------/---------/---------/---------/---------/--
      INTEGER n
      REAL*8 x1,x2,x(n),w(n)
      DOUBLE PRECISION EPS
      PARAMETER (EPS=3.d-14)
      INTEGER i,j,m
      DOUBLE PRECISION p1,p2,p3,pp,xl,xm,z,z1

      m=(n+1)/2          !The roots are symmetric in the interval, so we only
      xm=0.5d0*(x2+x1)   !have to find half of them
      xl=0.5d0*(x2-x1)

* Loop over the desired roots:

      do 12 i=1,m
        z=cos(3.141592654d0*(i-.25d0)/(n+.5d0))
              !Starting with the above approximation to the ith root, we enter 
              !the main loop of refinement by Newton's method.
 1      continue
          p1=1.d0
          p2=0.d0

          do 11 j=1,n         !Loop up the recurrence relation to get Legendre
            p3=p2             !polynomial evaluated at z.
            p2=p1
            p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
 11       continue

* p1 is now the desired  Legendre polynomial. We next compute pp, its derivative,
* by a standard relation involving also p2, the polynomial of one lower order:

          pp=n*(z*p1-p2)/(z*z-1.d0)
          z1=z
          z=z1-p1/pp                   !Newton's method

        if (abs(z-z1).gt.EPS) goto 1

* Scale the root to the desired interval, and put in its symmetric counterpart:
        x(i)=xm-xl*z                   
        x(n+1-i)=xm+xl*z               

* Compute the weight and its symmetric counterpart:
        w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)  
        w(n+1-i)=w(i)    
               
 12   continue

      return
      END

C********************************************************************

      SUBROUTINE ZSUR(A,INT,X,N,NC,EMACH)  
C     ------------------------------------------------------------------  
C     ZSUR IS  A STANDARD BACK-SUBSTITUTION  SUBROUTINE  USING THE   
C     OUTPUT OF ZGE TO CALCULATE X TIMES A-INVERSE, RETURNED IN X  
C     ------------------------------------------------------------------  
      IMPLICIT NONE 
C  
C ..  SCALAR ARGUMENTS  ..  
C  
      INTEGER N,NC  
      REAL*8 EMACH  
C  
C ..  ARRAY ARGUMENTS  ..  
C  
      INTEGER    INT(NC)  
      COMPLEX*16 A(NC,NC),X(NC)  
C  
C ..  LOCAL SCALARS  ..  
C  
      INTEGER    I,II,IN,J,IJ  
      COMPLEX*16 DUM  
C  
C ..  INTRINSIC FUNCTIONS  ..  
C  
*      INTRINSIC ABS  
C     ------------------------------------------------------------------  
C  
      DO 5 II=2,N  
      I=II-1  
      IF(INT(I)-I)1,2,1    !If INT(I).NE.I switch the I-th and INT(I)-th 
                           !elements of the vector X
   1  IN=INT(I)  
      DUM=X(IN)  
      X(IN)=X(I)  
      X(I)=DUM  
*
* Forming a matrix product
*
   2  DO 4 J=II,N  
      IF(ABS(A(I,J))-EMACH)4,4,3  
   3  X(J)=X(J)-X(I)*A(I,J)
   4  CONTINUE  

   5  CONTINUE               !the I-th row of A multiplied by X(I) 
                             !subtracted from X 
*
      DO 10 II=1,N  
      I=N-II+1  
      IJ=I+1  
      IF(I-N)6,8,6  
   6  DO 7 J=IJ,N  
   7  X(I)=X(I)-X(J)*A(J,I)
   8  IF(ABS(A(I,I))-EMACH*1.0D-7)9,10,10  
   9  A(I,I)=EMACH*1.0D-7*(1.D0,1.D0)  
  10  X(I)=X(I)/A(I,I)  
      RETURN  
      END  

C********************************************************************

      SUBROUTINE ZGER(A,INT,N,NC,EMACH)  
C     ------------------------------------------------------------------  
C     ZGE IS A STANDARD SUBROUTINE TO PERFORM GAUSSIAN ELIMINATION ON  
C     A NC*NC MATRIX 'A' PRIOR  TO INVERSION, DETAILS STORED IN 'INT'  
C                   MAKES AN LOWER DIAGONAL MATRIX
C     THIS ROUTINE DOES NOT BOTHER ABOUT ELEMENTS DIRECTLY ABOVE
C     THE MATRIX DIAGONAL AS THEY ARE NOT USED EXPLICITLY IN AN
C     ACCOMAPANYING ZSE ROUTINE
C     ------------------------------------------------------------------  
      IMPLICIT NONE 
C  
C ..  SCALAR ARGUMENTS  ..  
C  
      INTEGER N,NC  
      REAL*8 EMACH  
C  
C ..  ARRAY ARGUMENTS  ..  
C  
      INTEGER    INT(NC)  
      COMPLEX*16 A(NC,NC)  
C  
C ..  LOCAL SCALARS  ..  
C  
      INTEGER    I,II,IN,J,K  
      COMPLEX*16 YR,DUM  
C  
C ..  INTRINSIC FUNCTIONS  ..  
C  
*      INTRINSIC ABS  
C     ------------------------------------------------------------------  
C  
      DO 10 II=2,N  
      I=II-1  
      YR=A(I,I)  
      IN=I  
*
* Finding an element with the largest magnitude in the I-th column
* below the matrix diagonal (including the diag. element):

      DO 2 J=II,N  
      IF(ABS(YR)-ABS(A(I,J)))1,2,2  
   1  YR=A(I,J)  
      IN=J  
   2  CONTINUE  
      INT(I)=IN      !The largest element in the I-th row above the matrix 
                     !diagonal is in the IN-th row and is denoted by YR
* 
      IF(IN-I)3,5,3  !If IN.NE.I switch the I-th and IN-th columns to the
   3  DO 4 J=I,N     !left of the matrix diagonal (including the diag. element)
      DUM=A(J,I)  
      A(J,I)=A(J,IN)  
   4  A(J,IN)=DUM  

   5  IF(ABS(YR)-EMACH)10,10,6  

*
* Gaussian elemination of matrix elements above the matrix diagonal.
* Subtraction of (A(I,J)/A(I,I)) multiple of the Ith column from the Jth column 
* in the sub-matrix beginning with the (I+1,I+1) diag. element
*
   6  DO 9 J=II,N  
      IF(ABS(A(I,J))-EMACH)9,9,7  
   7  A(I,J)=A(I,J)/YR  
      DO 8 K=II,N  
   8  A(K,J)=A(K,J)-A(I,J)*A(K,I)   !k-t element of j-th column
   9  CONTINUE                 !The elements in the Ith row
                               !above diagonal are not set to zero
*
  10  CONTINUE                 !end of "column loop"
      RETURN  
      END  

C (C) Copr. 03/2003  Alexander Moroz