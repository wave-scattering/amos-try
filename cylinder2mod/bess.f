      SUBROUTINE BESS (NMAX,zr,zj,zdj)
C--------/---------/---------/---------/---------/---------/---------/--
C >>> ZR,NMAX
C <<< Output ZJ,ZDJ
C==========================================================
C  Generates Bessel functions for each Gauss integration point
C
C  X =(2\pi/\lambda)*r
C  ZR=(XR,XI)
C  XR=(2\pi/\lambda)*r*MRR, MRR ... real part of the rel. refractive index
C  XI=(2\pi/\lambda)*r*MRI, MRI ... imag. part of the rel. refractive index
C  ZJ... arrays of Bessel functions
C  ZDJ... arrays of Bessel functions derivatives of the form
C                           [xf(x)]'/x                   (A)
C                 where prime denotes derivative with respect to x.
C                 (Note that Bessel function derivatives enter Eqs. (39)
C                  \cite{TKS} only in the (A) combination!!!!)
C  NMAX   ... angular momentum cutoff
C  ADJR+I*ADJI - FUNCTION (1/X)(D/DX)(X*J(X))=J(X)/X + J'(X)                
C--------/---------/---------/---------/---------/---------/---------/--
      IMPLICIT NONE
      Integer N,NMAX,NMAXD
c  NMAXD here should be by 2 larger than NMAXD in BESS      
      Parameter (NMAXD=62)   
      REAL*8 YR,YI
      REAL*8 AJR(NMAXD),AJI(NMAXD),ADJR(NMAXD),ADJI(NMAXD)
      COMPLEX*16 ZR,ZJ(NMAX),ZDJ(NMAX)
*
c           CALL RJB(XX,AJ,ADJ,NMAX,NNMAX1)
c           CALL RYB(XX,AY,ADY,NMAX)
*
           YR=DBLE(ZR)
           YI=IMAG(ZR)
*
           CALL CJB(YR,YI,AJR,AJI,ADJR,ADJI,NMAX,2)
*
           DO 10 N=1,NMAX
                ZJ(N)=dcmplx(AJR(N),AJI(N))
                ZDJ(N)=dcmplx(ADJR(N),ADJI(N))-ZJ(N)/ZR
   10 CONTINUE

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
C   NNMAX=2 - angular momentum cutoff - DETERMINES NUMERICAL ACCURACY                                                                    
C--------/---------/---------/---------/---------/---------/---------/--
      IMPLICIT NONE
      Integer NMAXD
c  NMAXD here should be by 2 larger than NMAXD in BESS      
      Parameter (NMAXD=62)

      Integer NMAX,NNMAX,I,I1,L,L1
      REAL*8 XRXI,XR,XI,QF,QI,AR,AI,ARI,CXXR,CXXI,CZ0R,CZ0I,CR,CI,
     & CY0R,CY0I,CY1R,CY1I,CU1R,CU1I,CYIR,CYII,CUIR,CUII,CYI1R,CYI1I

      REAL*8 YR(NMAXD),YI(NMAXD),UR(NMAXD),UI(NMAXD)
      REAL*8 CYR(NMAXD),CYI(NMAXD),CZR(1200),CZI(1200)
c     *       CUR(NMAXD),CUI(NMAXD)
*
      L=NMAX+NNMAX
      XRXI=1D0/(XR*XR+XI*XI)
      CXXR=XR*XRXI             !Re [1/(XR+i*XI)]
      CXXI=-XI*XRXI            !Im [1/(XR+i*XI)] 
      QF=1D0/DBLE(2*L+1)
      CZR(L)=XR*QF
      CZI(L)=XI*QF
      L1=L-1
      DO I=1,L1
         I1=L-I
         QF=DBLE(2*I1+1)
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
         QI=DBLE(I)
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
