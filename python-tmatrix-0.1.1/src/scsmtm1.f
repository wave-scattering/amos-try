C      program scsmtm
c Calculation of the T matrix and random-orientation scattering
c matrix for an ensemble of spheres.                                   c
c                                                                      c
c Usage:                                                               c
c                                                                      c
c   a.out inputfile                                                    c
c                                                                      c
c where                                                                c
c                                                                      c
c   a.out: executable file of this code                                c
c   inputfile: file containing the input parameters for calculation.   c
c                                                                      c
c The input parameters are                                             c
c                                                                      c
c  fpos: sphere data input file.                                       c
c  npart: number of spheres to be used in calculations                 c
c  xscale, rirscale, riiscale: scaling n,k: scaling factors for size   c
c    parameter and real and imaginary parts of refractive index for    c
c    sphere data listed in fpos.                                       c
c       spheres.  (k is positive)                                      c
c  fout: output file for results.                                      c
c  ftmtx: T matrix output file.                                        c
c  itermax: max number of iterations to be used in                     c
c           iterative solution of cluster T matrix (default:100)       c
c  eps: error tolerance in iterative solution (default: 1e-8)          c
c  rlx: relaxation parameter in iterative solution. rlx=1 gives        c
c       standard unrelaxed interation.  rlx<1 is underrelaxation.      c
c       rlx=0 switches to a conjugate gradient method.                 c
c  qeps1: error tolerance for determination of the number of           c
c         harmonic orders used to calculate the field from             c
c         each sphere (default: 1e-4).                                 c
c  qeps2: error tolerance for determination of the number of           c
c         harmonic orders in the cluster T matrix (1e-4)               c
c  nodri: if qeps1=0, all sphere expansions are truncated at n=nodri   c
c  thetamin, thetamax, ntheta: minimum, maximum, and number of angles  c
c        over which the random orientation scattering matrix elements  c
c        are calculated, in degrees.                                   c
c  x,xp,yp,zp,re(m),im(m)  (i=1,npart): size parameter, cartesian      c
c        coordinates in the cluster, and the real and imaginary parts  c
c        of the sphere refractive index, listed sequentially.          c
c        The units of size parameter and position are arbitrary -      c
c        touching identical spheres will then be separated by two size c
c        parameters.  Size parameters and positions are  multplied by  c
c        xscale, Re(m) by rirscale, and Im(m) by riiscale to obtain    c
c        the actual value used in the computations.                    c
c                                                                      c
c  This is an example of the inputfile format                          c
c                                                                      c
c pos.dat           'position file (blank if data below is used)'      c
c 7                 'npart (number of spheres)'                        c
c .25,2.,.1         'x, re(m), im(m) scaling factors'                  c
c test.dat          'output file'                                      c
c tmtx.dat          't matrix file'                                    c
c 100,1.d-8,0       'max iterations, solution eps, rlx parameter'      c
c 1.d-4,1.d-5       'qeps1,qeps2'                                      c
c 0,180,91          'thetamin,thetamax,ntheta'                         c
c 1.,0.,0.,0.,1.,1.                                                    c
c 1.,0.,0.,2.,1.,1.                                                    c
c 1.,0.,2.,0.,1.,1.                                                    c
c 1.,2.,0.,0.,1.,1.                                                    c
c 1.,0.,0.,-2.,1.,1.                                                   c
c 1.,0.,-2.,0.,1.,1.                                                   c
c 1.,-2.,0.,0.,1.,1.                                                   c
c                                                                      c
c Note that all positions and size parameters would be multipled by    c
c xscale=0.25 in the above example.   The sphere configuration         c
c represents 7 spheres located symmetrically on a lattice, and the     c
c spheres are in contact. Refractive index is 2+i, and is the same     c
c for all spheres.                                                     c
c                                                                      c
c Keyboard input option: if inputfile is left blank in the command     c
c line, the code prompts the user to input parameters.                 c
c                                                                      c
c If your compiler does not recognize the getarg(1,fdat) subroutine    c
c (which returns the first command line argument as the string fdat)   c
c substitute the corresponding command for your compliler or `hard     c
c wire' the name of the inputfile into the code.                       c
c                                                                      c
c                                                                      c
c     NOTES:                                                           c
c                                                                      c
c     position  file:  file  name  for  a file containing the size     c
c     parameters,  positions  and  refractive index values for the     c
c     npart  spheres.  The format of this file is the same as that     c
c     given  in  the  above  example,  i.e.,  a sequential list of     c
c     x,xp,yp,zp,  re(m),  im(m) for the NPART spheres. Scaling by     c
c     xscale,  rirscale,  and  riiscale  is  applied  as described     c
c     above.  If  NPART>  the  number of data entries in the file,     c
c     NPART  is  changed to the number of entries. If fpos is left     c
c     blank  in the input file, the data is assumed to be appended     c
c     to  the  end  of  the  input  file  (as shown in the example     c
c     above).                                                          c
c                                                                      c
c   output file: file to which the results are written                 c
c                                                                      c
c  T matrix file: file to which the T matrix is written.  See output   c
c        lines in code for format.  Leave blank if no output is        c
c        desired:  Warning:  this file can be huge (10's of MB)        c
c                                                                      c
c     Solution   parameter  rlx:  typically  rlx=0.8-0.9  works        c
c     fastest.  If  solution does not converge (i.e., iter>itermax),   c
c     change  to rlx=0 (conjugate gradient). Non-convergence usually   c
c     occurs   only   for   highly   conducting  spheres.  rlx=0  is   c
c     recommended if high accuracy in off-diagonal scattering matrix   c
c     elements is required.                                            c
c                                                                      c
c     Setting  itermax=0  will  give the non-interacting result,       c
c     i.e.,  the  result  corresponding  to  far-field wave addition   c
c     among the NPART spheres with no multiple scattering (sometimes   c
c     known  as  the  RGD approximation). The formulation used here,   c
c     however,  is  not  the  most  efficient  method to get the RGD   c
c     random-orientation approximation from a cluster of spheres.      c
c                                                                      c
c     qeps1  and  qeps2:  making these smaller will increase the       c
c     accuracy   of   the   calculations,  but  at  the  expense  of   c
c     significantly increased computational overhead. qeps1 controls   c
c     N(i)  (number  of  harmonic  orders for the sphere), and qeps2   c
c     controls  N  (number of orders of T matrix). The computational   c
c     overhead scales as NPART*N(i)^2 and N^2. N(i) is calculated as   c
c     the  number of orders required to give a Lorenz/Mie extinction   c
c     coefficient  within  qeps1  to the N(i)->infinity result, on a   c
c     relative  basis.  N  is  more complicated. An adequately small   c
c     value of qeps1 can be determined only by inspection of results   c
c     for  progressively  smaller  values  of  qeps1. qeps2 is small   c
c     enough  (usually!)  if  the  cluster  extinction coefficients,   c
c     calculated  from  the  sphere-centered  and cluster-centered T   c
c     matrices, are in close agreement.                                c
c                                                                      c
c     Setting qeps1=-N will allow the user to override the test and    c
c     fix the truncation limit at N for all spheres.                   c
c                                                                      c
c     The  code  first  calculates the cluster T matrix, per the       c
c     interative procedure described in "Calculation of the T matrix   c
c     and  scattering  matrix  for  ensembles  of  spheres,"  D.  W.   c
c     Mackowski  and  M.  I.  Mishchenko, J. Opt. Soc. Amer. A., 13,   c
c     2266-2278  (1996). It then calculates the orientation averaged   c
c     scattering   matrix   elements   for   the  cluster,  using  a   c
c     generalized spherical function formulation.                      c
c                                                                      c
c     The formulation here uses a different normalization than that    c
c     given in the paper.   Specifically, the addition coefficients    c
c     are scaled according to                                          c
c                                                                      c
c      A_{mnkl}(norm) = (E_mn/E_kl)^(1/2) A_{mnkl)                     c
c                                                                      c
c     and likewise for other quantities appearing in the analysis.     c
c     This eliminates the factorial functions throughout the code,     c
c     and greatly reduces the risk of overflow problems.               c
c                                                                      c
c     Efficiency  factors  for the cluster are defined with            c
c     respect  to  the  volume-mean  radius of the cluster. The        c
c     matrix  elements  are  normalized so that S11, integrated        c
c     over  theta,  is equal to Qsca. Extinction and absorption        c
c     factors  for  the  individual  spheres  are  defined with        c
c     respect to the sphere radius                                     c
c                                                                      c
c                                                                      c
c      The  expansion  coefficients  for  the  random  orientation     c
c      scattering  matrix are  written  to  the  output  file. The     c
c      formulas for the matrix elements are                            c
c                                                                      c
c     s11=sum_{n=0}^{2N} a11(n) d(0,0,n)                               c
c     s22=sum_{n=2}^{2N}(a22m(n) d(2,-2,n) + a22p(n) d(2,2,n))         c
c     s33=sum_{n=2}^{2N}(-a22m(n) d(2,-2,n) + a22p(n) d(2,2,n))        c
c     s44=sum_{n=0}^{2N} a44(n) d(0,0,n)                               c
c     s12=sum_{n=2}^{2N} a12(n) d(2,0,n)                               c
c     s23=sum_{n=2}^{2N}(a23m(n) d(2,-2,n) + a23p(n) d(2,2,n))         c
c     s34=sum_{n=2}^{2N} a34(n) d(2,0,n)                               c
c     s13=sum_{n=2}^{2N} a34(n) d(2,0,n)                               c
c     s24=sum_{n=2}^{2N} a24(n) d(2,0,n)                               c
c     s14=sum_{n=0}^{2N} a14(n) d(0,0,n)                               c
c                                                                      c
c     where d(m,k,n) is returned by subroutine ROTCOEF, as a function  c
c     of scattering angle, as DC(m,n*(n+1)+k), and N is the truncation c
c     limit of the T matrix.                                           c
c                                                                      c
c Important parameter definitions:                                     c
c                                                                      c
c  npd: maximum number of spheres in calculations                      c
c  nod: maximum order of sphere expansions                             c
c  notd: maximum order of T matrix expansion.                          c
c                                                                      c
c These are defined in an included file sctmdim.f.   They should     c
c be set to reflect the expected size of your problem and the          c
c memory capacity of your machine.   An example of sctmdim.f is      c
c given below:                                                         c
c                                                                      c
c    parameter(npd=50,nod=5,notd=30)                                   c
c                                                                      c
c It is difficult to estimate how much memory the code will use.       c
c You will obviously know if it uses too much.                         c
c                                                                      c
c Be sure to compile with the -fast option (or whatever option, on     c
c your compiler, that generates the fastest executable).               c
c                                                                      c
c The code was revised on 30 October 1998.                             c
c                                                                      c
c Further modifications:                                               c
c  16 January 1999: fixed T matrix convergence calculation to force    c
c        calculation of at least NODRMAX orders.  This insures that    c
c        the T matrix for a single sphere will be calculated           c
c        correctly                                                     c
c                                                                      c
c Questions/comments: contact                                          c
c    Daniel Mackowski                                                  c
c    dmckwski@eng.auburn.edu                                           c


      subroutine tmnsrandom(npart,xi,xp,yp,zp,sni,ski,
     &                      xscale,rscale,rirscale,riiscale,
     &                      itermax,rlx,eps,qeps1,qeps2,theta,
     &                      qet,qat,qst,g,qei,qai,sm,errcode)

      implicit real*8 (a-h,o-z)
      include 'tmnsdim.f'
c
      parameter(nbd=nod*(nod+2),nbd2=nbd+nbd,
     1          nbtd=notd*(notd+2),nbt1=notd*notd+4*notd+3,
     1          ntd=npd*nbd,notd2=notd+notd,nfd=notd2*(notd2+2),
     1          nbc=3*notd2+6)
      integer nodr(npd)
      real*8 xi(npd),sni(npd),ski(npd),
     1       xp(npd),yp(npd),zp(npd),qai(npd),qei(npd)
      complex*16 ci,tc(2,nbtd,2,nbtd)
      common/consts/bcof(0:nbc,0:nbc),fnr(0:2*nbc)
C      character fpos*30,fposo*30,fout*30,fdat*30,ftmtx*30
      data ci/(0.d0,1.d0)/
      real*8 sm(4,4)
      integer itermax
      real*8 rlx,eps,qeps1,qeps2,theta
C      data itermax,eps,rlx,qeps1,qeps2/100,1.d-8,.9,1.d-4,1.d-4/

      errcode=0

c
c calculation of constants
c
      pi=4.d0*datan(1.d0)
      do n=1,2*nbc
         fnr(n)=dsqrt(dble(n))
      enddo
      bcof(0,0)=1.d0
      do n=0,nbc-1
         do l=n+1,nbc
            bcof(n,l)=fnr(n+l)*bcof(n,l-1)/fnr(l)
            bcof(l,n)=bcof(n,l)
         enddo
         bcof(n+1,n+1)=fnr(n+n+2)*fnr(n+n+1)*bcof(n,n)/fnr(n+1)/fnr(n+1)
      enddo

C     I (Ovidio) added the next three lines (27/02/2010)
      tmin=theta
      tmax=theta
      nt=1
c
c data input
c
C      call getarg(1,fdat)
C   10 if(fdat.eq.' ') then
C         write(*,'('' particle position file:'',$)')
C         read(*,'(a)') fpos
C         if(fpos.eq.'/') fpos=fposo
C         fposo=fpos
C         if(fpos.ne.' ') open(1,file=fpos)
C         write(*,'('' npart:'',$)')
C         read(*,*) npart
C         write(*,'('' xscale, rirscale,riiscale:'',$)')
C         read(*,*) xscale,rirscale,riiscale
C         write(*,'('' output file (return for none):'',$)')
C         read(*,'(a)') fout
C         if(fout.eq.'/') fout=' '
C         write(*,'('' T matrix output file (return for none):'',$)')
C         read(*,'(a)') ftmtx
C         if(ftmtx.eq.'/') ftmtx=' '

C         write(*,'('' itermax, eps, rlx, qeps1, qeps2:'',$)')
C         read(*,*) itermax,eps,rlx,qeps1,qeps2
C         if(qeps1.lt.0) nodri=-qeps1
C         write(*,'('' thetamin, thetamax, ntheta:'',$)')
C         read(*,*) tmin,tmax,nt
C      else
C         open(1,file=fdat)
C         read(1,'(a)') fpos
C         fpos=fpos(:index(fpos,' '))
C         read(1,*) npart
C         read(1,*) xscale,rirscale,riiscale
C         read(1,'(a)') fout
C         fout=fout(:index(fout,' '))
C         read(1,'(a)') ftmtx
C         ftmtx=ftmtx(:index(ftmtx,' '))
C         read(1,*) itermax,eps,rlx
C         read(1,*) qeps1,qeps2
C         read(1,*) tmin,tmax,nt
C         if(fpos.ne.' ') then
C            close(1)
C            open(1,file=fpos)
C         endif
C      endif
      xv=0.
      xm=0.
      ym=0.
      zm=0.
      do i=1,npart
C         if(fdat.eq.' '.and.fpos.eq.' ') then
C            write(*,'('' x, xp, yp, zp, re(m), im(m)#'',i2,'':'',$)') i
C            read(*,*) xii,xpi,ypi,zpi,rri,rii
C         else
C            read(1,*,end=25,err=25) xii,xpi,ypi,zpi,rri,rii
C         endif
         xi(i)=xscale*xi(i)
         xp(i)=xscale*xp(i)
         yp(i)=xscale*yp(i)
         zp(i)=xscale*zp(i)
         sni(i)=sni(i)*rirscale
         ski(i)=ski(i)*riiscale
         xm=xm+xp(i)
         ym=ym+yp(i)
         zm=zm+zp(i)
         xv=xv+xi(i)*xi(i)*xi(i)
         if(qeps1.lt.0.) nodr(i)=-qeps1
      enddo
   25 npart=i-1
C      close(1)
c
c xv: volume mean size parameter
c
c this sorts the data according to distance from the cluster
c origin -- which is useful in the solution procedure
c
      xv=xv**(1./3.)
      xm=xm/dble(npart)
      ym=ym/dble(npart)
      zm=zm/dble(npart)
      do i=npart-1,1,-1
         do j=1,i
            j1=j+1
            xj=xp(j)-xm
            yj=yp(j)-ym
            zj=zp(j)-zm
            xj1=xp(j1)-xm
            yj1=yp(j1)-ym
            zj1=zp(j1)-zm
            rj=xj*xj+yj*yj+zj*zj
            rj1=xj1*xj1+yj1*yj1+zj1*zj1
            if(rj.gt.rj1) then
               call rswap(xp(j),xp(j1))
               call rswap(yp(j),yp(j1))
               call rswap(zp(j),zp(j1))
               call rswap(xi(j),xi(j1))
               call rswap(sni(j),sni(j1))
               call rswap(ski(j),ski(j1))
            endif
         enddo
      enddo
c
c calculate the T matrix
c
30    niter=itermax
      call tmtrx(npart,xp,yp,zp,sni,ski,xi,nodr,
     1            nodrt,niter,eps,qeps1,qeps2,rlx,tc,qei,qai)
      if(niter.ge.itermax) then
         errcode=1
         return
      endif
C      if(niter.ge.itermax) write(*,'('' maximum iterations '',
C     1 ''was exceeded'')')
C      write(*,'('' tmatrix order:'',i5)') nodrt
c
      nodrt2=2*nodrt
      nblkt=nodrt*(nodrt+2)
c
c This computes the total extinction and scattering efficiencies.
c efficiencies are defined with respect to xv.   Q's are calculated
c either directly from the T matrix or from the individual sphere
c results.   The two should agree.  If not, qeps2 is too large.
c
      qet=0.
      qst=0.
      do n=1,nblkt
         do ip=1,2
            qet=qet+tc(ip,n,ip,n)
            do l=1,nblkt
               do iq=1,2
                  qst=qst+tc(ip,n,iq,l)*conjg(tc(ip,n,iq,l))
               enddo
            enddo
         enddo
      enddo
      qet=2.*qet/xv/xv
      qst=2.*qst/xv/xv
      qat=qet-qst
C      qe=0.
C      qa=0.
C      do i=1,npart
C         qe=qe+qei(i)*xi(i)*xi(i)
C         qa=qa+qai(i)*xi(i)*xi(i)
C      enddo
C      qe=qe/xv/xv
C      qa=qa/xv/xv
C      qs=qe-qa
C      write(*,'('' qe, qa, qs (tmatrix):'',3e12.5)') qet,qat,qst
C      write(*,'('' qe, qa, qs (tij mtx):'',3e12.5)') qe,qa,qs
C      if(fdat.eq.' ') then
C         icalc=0
C         write(*,'('' recalculate T matrix (0/1):'',$)')
C         read(*,*) icalc
C         if(icalc.eq.1) then
C            write(*,'('' itermax, eps, rlx, qeps1, qeps2:'',$)')
C            read(*,*) itermax,eps,rlx,qeps1,qeps2
C            if(qeps1.lt.0) then
C               do i=1,npart
C                  nodr(i)=-qeps1
C               enddo
C            endif
C            goto 30
C         endif
C      endif

C      if(ftmtx.ne.' ') then
C         open(1,file=ftmtx)
C         write(1,'(2i5,e13.5)') nodrt,npart,xv
C         do n=1,nodrt
C            do m=-n,n
C               mn=n*(n+1)+m
C               do ip=1,2
C                  do l=n,nodrt
C                     do k=-l,l
C                        kl=l*(l+1)+k
C                        do iq=1,2
C                           if(cdabs(tc(ip,mn,iq,kl))/xv.gt..01*eps)
C     1                        then
C                              write(1,'(2e13.5)') tc(ip,mn,iq,kl)
C                           else
C                              write(1,*) 0,0
C                           endif
C                        enddo
C                     enddo
C                  enddo
C               enddo
C            enddo
C         enddo
C         close(1)
C      endif
c
c calculate the scattering matrix coefficients.
c nodrexp is the number of harmonic orders used in the generalized
c spherical function representation of the scattering matrix.
c Complete calculation of the scattering matrix requires nodrexp=
c 2*nodrt.   If only the asymmetry factor is required, set nodrexp=1
c
      nodrexp=2*nodrt
      call ranprops(tc,npart,nodrt,nodrexp,xscale,rirscale,riiscale,
     1    xv,nt,tmin,tmax,qet,qat,qst,a111,sm)
C     1    xv,fout,fdat,fpos,nt,tmin,tmax,qet,qat,qst,a111,sm)

C  200 if(fdat.eq.' ') then
C         write(*,'('' more (0/1):'',$)')
C         read(*,*) more
C         if(more.eq.1) goto 10
C      endif
C      stop
      return
      end
c
      subroutine ranprops(tc,npart,nodrt,nodrset,xscale,rirscale,
     *    riiscale,xv,nt,tmin,tmax,
C     *    riiscale,xv,fout,fdat,fpos,nt,tmin,tmax,
     *    qet,qat,qst,a111,sm)
      implicit real*8(a-h,o-z)
      include 'tmnsdim.f'
c
      parameter(nbd=nod*(nod+2),nbd2=nbd+nbd,
     1          nbtd=notd*(notd+2),nbt1=notd*notd+4*notd+3,
     1          ntd=npd*nbd,notd2=notd+notd,nfd=notd2*(notd2+2),
     1          nbc=3*notd2+6)
      real*8 dc(-2:2,0:nfd)
      real*8 sm(4,4)
      complex*16 aw(0:2,-1:1,0:notd2),bw(0:2,-1:1,0:notd2),
     1           cw(0:notd2),dw(0:notd2)
      complex*16 ci,tc(2,nbtd,2,nbtd)
      common/consts/bcof(0:nbc,0:nbc),fnr(0:2*nbc)
C      character fout*30,fdat*30,fpos*30
      data ci/(0.d0,1.d0)/

      pi=4.d0*datan(1.d0)
      nblkt=nodrt*(nodrt+2)
      qet=0.
      qst=0.
      do n=1,nblkt
         do ip=1,2
            qet=qet+tc(ip,n,ip,n)
            do l=1,nblkt
               do iq=1,2
                  qst=qst+tc(ip,n,iq,l)*conjg(tc(ip,n,iq,l))
               enddo
            enddo
         enddo
      enddo
      qet=2.*qet/xv/xv
      qst=2.*qst/xv/xv
      qat=qet-qst
      if(nodrset.eq.-1) then
         nodrexp=2*nodrt
      else
         nodrexp=nodrset
      endif
      call tmsm(tc,nodrt,nodrexp,nbtd,aw,bw,cw,dw)
c
c coefficients are scaled with respect to 2/xv^2, so
c a110 is the scattering efficiency, and a111 is the asymmetry
c factor
c
      a110=(aw(0,-1,0)+aw(0,1,0))*2./xv/xv
      do iw=0,nodrexp
         do k=-1,1
            do i=0,2
               aw(i,k,iw)=aw(i,k,iw)*2./xv/xv
               bw(i,k,iw)=bw(i,k,iw)*2./xv/xv
            enddo
         enddo
         cw(iw)=cw(iw)*2./xv/xv
         dw(iw)=dw(iw)*2./xv/xv
         a11=aw(0,-1,iw)+aw(0,1,iw)
         snr=a11/a110/dble(iw+iw+1)
         if (iw.eq.1) a111=snr
      enddo
c
c output file operations
c
C      if(fout.ne.' ') then
C         open(1,file=fout)
C         write(1,'('' scsmtm results'')')
C         if(fdat.eq.' ') then
C            write(1,'('' position file:'',a)') fpos
C         else
C            write(1,'('' input file:'',a)') fdat
C         endif
C         write(1,'('' number of spheres:'',i6)') npart
C         write(1,'('' xscale, rirscale, riiscale, xv:'',4e12.4)')
C     *         xscale,rirscale,riiscale,xv
C         write(1,'('' qext,qabs,qsca,<cos theta>:'',4e12.4)')
C     1      qet,qat,qst,a111
C         write(1,'('' theta    s11         s22         s33'',
C     1 ''         s44'',
C     1 ''         s21         s32         s43         s31''
C     2 ''         s42         s41'')')
C      endif
c
c scattering matrix calculation
c
100   do i=1,nt
         if(nt.eq.1) then
            th=tmin
         else
            th=tmin+(tmax-tmin)*dble(i-1)/dble(nt-1)
         endif
         ct=cos(th*pi/180.d0)
c
c  dc is the normalized generalized spherical function
c  dc(k,n*(n+1)+m) = ((n-k)!(n+m)!/(n+k)!/(n-m)!)^(1/2) D^k_{mn},
c  where D^k_{mn} is defined in M&M JOSA 96
c
         call rotcoef(ct,2,nodrexp,dc,2)
         s11t=0.
         s12t=0.
         s13t=0.
         s14t=0.
         s21t=0.
         s22t=0.
         s23t=0.
         s24t=0.
         s31t=0.
         s32t=0.
         s33t=0.
         s34t=0.
         s41t=0.
         s42t=0.
         s43t=0.
         s44t=0.
         do n=0,nodrexp
            nn0=n*(n+1)
            nnp2=nn0+2
            nnm2=nn0-2
            s11t=s11t+dc(0,nn0)*(aw(0,-1,n)+aw(0,1,n))
            s14t=s14t+dc(0,nn0)*(aw(0,1,n)-aw(0,-1,n))
            s44t=s44t+dc(0,nn0)*(bw(0,1,n)-bw(0,-1,n))
            s41t=s41t+dc(0,nn0)*(bw(0,-1,n)+bw(0,1,n))
            if(n.ge.2) then
               s12t=s12t+dc(2,nn0)*(aw(2,-1,n)+aw(2,1,n))
               s24t=s24t+dc(2,nn0)*(aw(2,1,n)-aw(2,-1,n))
               s34t=s34t+dc(2,nn0)*dimag(bw(2,-1,n)-bw(2,1,n))
               s31t=s31t-dc(2,nn0)*dimag(bw(2,-1,n)+bw(2,1,n))
               s13t=s13t+2.*dc(2,nn0)*dimag(aw(2,0,n))
               s42t=s42t+dc(2,nn0)*2.*bw(2,0,n)
               s22t=s22t+dc(2,nnp2)*cw(n)+dc(2,nnm2)*dw(n)
               s23t=s23t+dimag(dc(2,nnp2)*cw(n)+dc(2,nnm2)*dw(n))
               s33t=s33t+dc(2,nnp2)*cw(n)-dc(2,nnm2)*dw(n)
               s32t=s32t-dimag(dc(2,nnp2)*cw(n)-dc(2,nnm2)*dw(n))
            endif
         enddo
c
c here are the VV and HH differential cross sections
c
c         gvv=.25*(s11t+s22t-2.*s12t)
c         ghh=.25*(s11t+s22t+2.*s12t)
c
c only the diagonal and upper triangular elements of the
c scattering matrix are written
c
C         if(fout.ne.' ')
C     1    write(1,'(f8.2,10e12.4)') th,s11t,s22t,s33t,s44t,s12t,s32t,
C     1    s34t,s31t,s42t,s41t

C        I (Ovidio) added the following 16 lines (27/02/2010)
         sm(1,1)=s11t
         sm(1,2)=s12t
         sm(1,3)=s13t
         sm(1,4)=s14t
         sm(2,1)=s21t
         sm(2,2)=s22t
         sm(2,3)=s23t
         sm(2,4)=s24t
         sm(3,1)=s31t
         sm(3,2)=s32t
         sm(3,3)=s33t
         sm(3,4)=s34t
         sm(4,1)=s41t
         sm(4,2)=s42t
         sm(4,3)=s43t
         sm(4,4)=s44t
      enddo
C      if(fout.ne.' ') then
C         write(1,'('' scattering matrix expansion coefficients'')')
C         write(1,'(''    w  a11         a22m        a22p        '',
C     1    ''a23m        a23p        a44         a12         '',
C     1    ''a34         a13         a24         a14'')')
C         do iw=0,nodrexp
C            a11=dble(aw(0,-1,iw)+aw(0,1,iw))
C            a22mr=dw(iw)
C            a22pr=cw(iw)
C            a23mi=dimag(dw(iw))
C            a23pi=dimag(cw(iw))
C            a44=-dble(bw(0,-1,iw)-bw(0,1,iw))
C            a12=dble(aw(2,-1,iw)+aw(2,1,iw))
C            a34=dimag(bw(2,-1,iw)-bw(2,1,iw))
C            a13=2.*dimag(aw(2,0,iw))
C            a24=-dble(aw(2,-1,iw)-aw(2,1,iw))
C            a14=-dble(aw(0,-1,iw)-aw(0,1,iw))
C            write(1,'(i5,11e12.4)') iw,a11,a22mr,a22pr,a23mi,a23pi,
C     1           a44,a12,a34,a13,a24,a14
C         enddo
C         close(1)
C      endif
      return
      end


      subroutine rswap(a,b)
      real*8 a,b,c
      c=a
      a=b
      b=c
      return
      end
c
      subroutine iswap(a,b)
      integer a,b,c
      c=a
      a=b
      b=c
      return
      end
c
c random-orientation scattering matrix expansion coefficient
c calculation
c
      subroutine tmsm(tc,nodr,nodrw,ntd,aw,bw,cw,dw)
      implicit real*8(a-h,o-z)
      include 'tmnsdim.f'
      parameter(nbd=notd*(notd+2),nod2=notd+notd,
     1           nbd2=nod2*(nod2+2),nbc=3*nod2+6)
      complex*16 aw(0:2,-1:1,0:nod2),bw(0:2,-1:1,0:nod2),
     1           cw(0:nod2),dw(0:nod2)
      complex*16 ci,cin,bm(2,nbd,2),pp(notd,2,2),
     1           dm(-notd-1:notd+1,3,notd,2,notd,2),a,am(2,notd+1,2)
      complex*16 fm(3,notd,2,notd,2)
      complex*16 tc(2,ntd,2,ntd)
      real*8 vc(0:2*nod2+2)
      common/consts/bcof(0:nbc,0:nbc),fnr(0:2*nbc)
      data ci/(0.d0,1.d0)/
c
      nodr2=nodr+nodr
      nblk=nodr*(nodr+2)
      do n=1,nodr
         do ip=1,2
            do l=1,nodr
               do iq=1,2
                  do k=1,3
                     do iu=-nodr-1,nodr+1
                        dm(iu,k,n,ip,l,iq)=0.
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
c
      do n=1,nodr
         cin=ci**(n+1)
         pp(n,1,1) =-.5d0*cin*fnr(n+n+1)
         pp(n,2,1) =-pp(n,1,1)
         pp(n,1,2)=-pp(n,1,1)
         pp(n,2,2)=pp(n,2,1)
      enddo
      nblk=nodr*(nodr+2)
      do n=1,nblk
         do ip=1,2
            do l=1,nodr
               do k=-l,l
                  kl=l*(l+1)+k
                  a=tc(ip,n,1,kl)
                  tc(ip,n,1,kl)=tc(ip,n,1,kl)*pp(l,1,1)
     1               +tc(ip,n,2,kl)*pp(l,1,2)
                  tc(ip,n,2,kl)=a*pp(l,2,1)
     1               +tc(ip,n,2,kl)*pp(l,2,2)
               enddo
            enddo
         enddo
      enddo
c
C      write(*,'('' calculating D matrix'')')
      do iw=0,nodr2
C         write(*,'('' w:'',i2,''/'',i2)') iw,nodr2
         do iv=-iw,iw
            do n=1,nblk
               do ip=1,2
                  do k=1,2
                     bm(k,n,ip)=0.
                  enddo
               enddo
            enddo
            do n=1,nodr
               nn1=n*(n+1)
               do l=max(1,abs(iw-n)),min(nodr,iw+n)
                  am(1,l,1)=0.
                  am(1,l,2)=0.
                  am(2,l,1)=0.
                  am(2,l,2)=0.
               enddo
               do it=-n,n
                  itn=nn1+it
                  lmax=min(nodr,iw+n)
                  call vcfunc(iv,iw,-it,n,lmax,vc)
                  do l=max(1,abs(iv-it),abs(n-iw)),lmax
                     ll1=l*(l+1)
                     itvl=ll1+it-iv
                     do ik=1,2
                        do ip=1,2
                           am(ik,l,ip)=am(ik,l,ip)
     1                      +vc(l)*tc(ip,itn,ik,itvl)
                        enddo
                     enddo
                  enddo
               enddo
c
               do m=-n,n
                  mn=nn1+m
                  do ik=1,2
                     k=-3+2*ik
                     iu=m-k
                     if(abs(iu).le.iw) then
                        lmax=min(nodr,iw+n)
                        call vcfunc(-iu,iw,m,n,lmax,vc)
                        do l=max(1,abs(iw-n)),lmax
                           fl=-(-1)**l*vc(l)/dble(l+l+1)
                           do ip=1,2
                              bm(ik,mn,ip)=bm(ik,mn,ip)
     1                         +am(ik,l,ip)*fl
                           enddo
                        enddo
                     endif
                  enddo
               enddo
            enddo
            do iu=-min(iw,nodr+1),min(iw,nodr+1)
               do iku=1,3
                  if(iku.eq.1) then
                     k=-1
                     k1=-1
                  elseif(iku.eq.2) then
                     k=1
                     k1=1
                  else
                     k=1
                     k1=-1
                  endif
                  m=iu+k
                  ns=max(1,abs(m))
                  ik=(k+1)/2+1
                  ik1=(k1+1)/2+1
                  m1=iu+k1
                  do n=ns,nodr
                     nu=n*(n+1)+m
                     n1s=max(1,abs(m1),n-nodrw)
                     n1e=min(nodr,n+nodrw)
                     do n1=n1s,n1e
                        cin=ci**(n-n1)
                        nu1=n1*(n1+1)+m1
                        fnn1=-fnr(n+n+1)*fnr(n1+n1+1)*dble(iw+iw+1)
                        do ip=1,2
                           do ip1=1,2
                              a=bm(ik,nu,ip)*cin*fnn1
     1                             *conjg(bm(ik1,nu1,ip1))
                              dm(iu,iku,n,ip,n1,ip1)
     1                            = dm(iu,iku,n,ip,n1,ip1)+a
                           enddo
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
c
C      write(*,'('' calculating scattering matrix expansion'')')
      do iw=0,nodrw
         do k=-1,1
            do i=0,2
               aw(i,k,iw)=0.
               bw(i,k,iw)=0.
            enddo
         enddo
         cw(iw)=0.
         dw(iw)=0.
      enddo
      do iw=0,nodrw
         iu2=1
         do n=1,nodr
            n1s=max(1,abs(n-iw))
            n1e=min(nodr,n+iw)
            do n1=n1s,n1e
               do ik=1,3
                  do ip=1,2
                     do ip1=1,2
                        fm(ik,n,ip,n1,ip1)=0.
                     enddo
                  enddo
               enddo
            enddo
         enddo
         do iu=-nodr-1,nodr+1
            do k=-1,1,2
               m=iu+k
               ik=(k+1)/2+1
               ns=max(1,abs(m))
               do n=ns,nodr
                  n1max=min(iw+n,nodr)
                  call vcfunc(m,n,0,iw,n1max,vc)
                  do n1=ns,nodr
                     if((n+n1.lt.iw).or.(abs(n-n1).gt.iw)) goto 60
                     fc1=-(-1)**n*vc(n1)*fnr(iw+iw+1)/fnr(n1+n1+1)
                     do ip=1,2
                        do ip1=1,2
                           fm(ik,n,ip,n1,ip1)=fm(ik,n,ip,n1,ip1)
     1                     +dm(iu,ik,n,ip,n1,ip1)*fc1
                        enddo
                     enddo
60                enddo
               enddo
            enddo
            if(iw.lt.2) goto 75
            m=iu+1
            m1=iu-1
            ns=max(1,abs(m))
            n1s=max(1,abs(m1))
            do n=ns,nodr
               n1max=min(iw+n,nodr)
               call vcfunc(m,n,-2,iw,n1max,vc)
               do n1=n1s,nodr
                  if((n+n1.lt.iw).or.(abs(n-n1).gt.iw)) goto 70
                  fc1=-(-1)**n*vc(n1)*fnr(iw+iw+1)/fnr(n1+n1+1)
                  do ip=1,2
                     do ip1=1,2
                        fm(3,n,ip,n1,ip1)=fm(3,n,ip,n1,ip1)
     1                   +dm(iu,3,n,ip,n1,ip1)*fc1
                     enddo
                  enddo
70             enddo
            enddo
            iu2=1
75       enddo
63       do n=1,nodr
            n1s=max(1,abs(n-iw))
            n1e=min(nodr,n+iw)
            in=(-1)**n
            n1max=min(iw+n,nodr)
            call vcfunc(1,n,0,iw,n1max,vc)
            do n1=n1s,n1e
               fc1=2.d0*in*vc(n1)*fnr(iw+iw+1)/fnr(n1+n1+1)
               i=mod(n+n1-iw,2)+1
               do ip=1,2
                  ip1=(2-i)*ip+(i-1)*(3-ip)
                  do k=-1,1,2
                     ik=(k+1)/2+1
                     aw(0,k,iw)=aw(0,k,iw)+fm(ik,n,ip,n1,ip1)*fc1
                     bw(0,k,iw)=bw(0,k,iw)+fm(ik,n,ip,n1,3-ip1)*fc1
                  enddo
                  bw(2,0,iw)=bw(2,0,iw)+fm(3,n,ip,n1,3-ip1)*fc1
                  aw(2,0,iw)=aw(2,0,iw)+fm(3,n,ip,n1,ip1)*fc1
               enddo
            enddo
            if(iw.lt.2) goto 80
            call vcfunc(1,n,-2,iw,n1max,vc)
            do n1=n1s,n1e
               fc2=2.d0*in*vc(n1)*fnr(iw+iw+1)/fnr(n1+n1+1)
               i=mod(n+n1-iw,2)+1
               do ip=1,2
                  ip1=(2-i)*ip+(i-1)*(3-ip)
                  do k=-1,1,2
                     ik=(k+1)/2+1
                     aw(2,k,iw)=aw(2,k,iw)+fm(ik,n,ip,n1,ip1)*fc2
     1                          *(-1)**ip1
                     bw(2,k,iw)=bw(2,k,iw)+fm(ik,n,ip,n1,3-ip1)*fc2
     1                          *(-1)**(3-ip1)
                  enddo
               enddo
               fc3=2.*(-1)**(n1+iw)*vc(n1)*fnr(iw+iw+1)/fnr(n1+n1+1)
               do ip=1,2
                  do ip1=1,2
                     cw(iw)=cw(iw)+fm(3,n,ip,n1,ip1)*fc2*(-1)**ip1
                     dw(iw)=dw(iw)+fm(3,n,ip,n1,ip1)*fc3*(-1)**ip
                  enddo
               enddo
            enddo
80       enddo
      enddo
      return
      end
c
c calculation of cluster T matrix via iteration scheme
c
      subroutine tmtrx(npart,xp,yp,zp,sni,ski,xi,nodr,
     1            nodrtmax,niter,eps,qeps1,qeps2,rlx,tc,qei,qai)
      implicit real*8(a-h,o-z)
      include 'tmnsdim.f'
      parameter(nbd=nod*(nod+2),nbc=6*notd+6,
     1          nbtd=notd*(notd+2),nrd=.5*(npd-1)*(npd-2)+npd-1)
      integer nodr(npd),nblk(npd)
      logical icon(npd)
      real*8 xi(*),sni(*),ski(*),rp(npd),qe1(npd),
     1       xp(*),yp(*),zp(*),qei(*),qai(*),qsq(npd),qsq0(npd)
      complex*16 ci,an1(2,nod,npd),cn1(2,nod,npd),tc(2,nbtd,2,nbtd)
      real*8 drot(-nod:nod,0:nbd,nrd)
      real*8 drott(-nod:nod,0:nbtd,npd)
      complex*16 ephi,anpt(2,nbtd),
     1        a,amnl(2,nod,nbd,nrd),
     1        amnlt(2,notd,nbd,npd),pnp(2,nbd,npd),
     1        ek(-nod:nod,nrd),anp(2,nbd,npd),
     1        ekt(-notd:notd,npd)
      common/consts/bcof(0:nbc,0:nbc),fnr(0:2*nbc)
      data ci/(0.d0,1.d0)/
      itermax=0
      xv=0.
      xm=0.
      ym=0.
      zm=0.
      nodrmax=0
c
c mie coefficients
c
      do i=1,npart
         call mie1(xi(i),sni(i),ski(i),nodr(i),qeps1,qe1(i),
     1             qs1,an1(1,1,i),cn1(1,1,i),0)
         nblk(i)=nodr(i)*(nodr(i)+2)
         xm=xm+xp(i)
         ym=ym+yp(i)
         zm=zm+zp(i)
         xv=xv+xi(i)*xi(i)*xi(i)
         qei(i)=0.
         qai(i)=0.
         qsq(i)=0.
         qsq0(i)=0.
         icon(i)=.false.
         nodrmax=max(nodr(i),nodrmax)
      enddo
C      write(*,'('' single sphere max. nodr:'',i5)') nodrmax
      if(nodrmax.eq.nod) then
         write(*,'('' Warning: single--sphere error tolerance'',
     1   '' may not be attained'')')
         write(*,'('' decrease qeps1 and/or increase nod'')')
      endif
      xm=xm/dble(npart)
      ym=ym/dble(npart)
      zm=zm/dble(npart)
      xv=xv**(1./3.)
      nodrtmax=0
c
c estimate 0-i translation orders: max radius of cluster
c
      do i=1,npart
         x=xm-xp(i)
         y=ym-yp(i)
         z=zm-zp(i)
         rp(i)=sqrt(x*x+y*y+z*z)
         xc=rp(i)+xi(i)
         nodrt=max(nodr(i),nint(xc+4.*xc**(1./3.))+5)
         nodrtmax=max(nodrt,nodrtmax)
      enddo
      if(nodrtmax.gt.notd) then
         write(*,'('' Warning: notd dimension is too small.'',
     1       '' increase to '',i5)') nodrtmax
      endif
      nodrtmax=min(nodrtmax,notd)
C      write(*,'('' estimated T matrix order:'',i5)') nodrtmax
c
c 0-i translation coefficients
c
      do i=1,npart
         x=xm-xp(i)
         y=ym-yp(i)
         z=zm-zp(i)
         if(rp(i).ne.0.) then
            call trancoef(1,rp(i),nodr(i),nodrtmax,amnlt(1,1,1,i),notd)
            if((x.eq.0.).and.(y.eq.0.)) then
               ephi=1.
            else
               ephi=cdexp(ci*datan2(y,x))
            endif
            ekt(0,i)=1.d0
            do m=1,nodrtmax
               ekt(m,i)=ephi**m
               ekt(-m,i)=dconjg(ekt(m,i))
            enddo
            ct=z/rp(i)
            call rotcoef(ct,nodr(i),nodrtmax,drott(-nod,0,i),nod)
         endif
      enddo
      nblktmax=nodrtmax*(nodrtmax+2)
c
      do i=1,npart
         do j=i+1,npart
            ij=.5*(j-1)*(j-2)+j-i
            x=xp(i)-xp(j)
            y=yp(i)-yp(j)
            z=zp(i)-zp(j)
            if((x.eq.0.).and.(y.eq.0.)) then
               ephi=1.
            else
               ephi=cdexp(ci*datan2(y,x))
            endif
            ek(0,ij)=1.d0
            do m=1,max(nodr(i),nodr(j))
               ek(m,ij)=ephi**m
               ek(-m,ij)=dconjg(ek(m,ij))
            enddo
            r=sqrt(x*x+y*y+z*z)
            ct=z/r
            call rotcoef(ct,nodr(i),nodr(j),drot(-nod,0,ij),nod)
            call trancoef(3,r,nodr(i),nodr(j),amnl(1,1,1,ij),nod)
         enddo
      enddo
      do n=1,nblktmax
         do ip=1,2
            do l=1,nblktmax
               do iq=1,2
                  tc(ip,n,iq,l)=0.
               enddo
            enddo
         enddo
      enddo
c
      qe=0.
C      write(*,'(''  n   # its  del qe        qe '',
C     1  ''      qs(tij mtrx) qs(t mtrx)   qs error   # conv sphrs '')')
c
c loop over columns
c
      do np=1,nodrtmax
         npblk=np*(np+2)
         np1=np*(np+1)
         itermaxn=0
         dqe=0.
         itermaxn=0
         do mp=-np,np
            mnp=np1+mp
            mnpm=np1-mp
            do ipp=1,2
c
c set up rhs
c
               enorm=0.
               do i=1,npart
                  if(icon(i)) goto 20
                  do n=1,nodr(i)
                     nn1=n*(n+1)
                     nmax=min(np,n)
                     do m=-n,n
                        mn=nn1+m
                        do ip=1,2
                           ia=mod(ip+ipp,2)+1
                           if(rp(i).eq.0.) then
                              if((n.eq.np).and.(m.eq.mp)
     1                         .and.(ip.eq.ipp)) then
                                 pnp(ip,mn,i)=an1(ip,n,i)
                              else
                                 pnp(ip,mn,i)=0.
                              endif
                           else
                              a=0.
                              do k=-nmax,nmax
                                 kn=nn1+k
                                 a=a+(-1)**k*drott(k,mnpm,i)
     1                           *amnlt(ia,np,kn,i)*drott(-k,mn,i)
                              enddo
                              pnp(ip,mn,i)=a*an1(ip,n,i)
     1                          *(-1)**m*ekt(mp,i)*ekt(-m,i)
                           endif
                           anp(ip,mn,i)=pnp(ip,mn,i)
                           enorm=enorm+pnp(ip,mn,i)
     1                          *conjg(pnp(ip,mn,i))
                        enddo
                     enddo
                  enddo
20             enddo
               if(enorm.eq.0) goto 80
c
c solver
c
               if(niter.ne.0) then
                  call itersoln(npart,nodr,nblk,icon,eps,niter,
     1              rlx,ek,drot,amnl,an1,pnp,anp,iter,err)
                  itermaxn=max(itermaxn,iter)
               endif
c
c translate to 0, compute single-sphere efficiencies
c
               do i=1,npart
                  qaii=0.
                  qsqi=0.
                  if(icon(i)) goto 75
                  do n=1,nodr(i)
                     nn1=n*(n+1)
                     do m=-n,n
                        mn=nn1+m
                        anpt(1,mn)=anp(1,mn,i)
                        anpt(2,mn)=anp(2,mn,i)
                        anp12=anpt(1,mn)*conjg(anpt(1,mn))
                        anp22=anpt(2,mn)*conjg(anpt(2,mn))
                        qsqi=qsqi+anp12+anp22
                        qaii=qaii+dble(cn1(1,n,i))*anp12
     1                           +dble(cn1(2,n,i))*anp22
                     enddo
                  enddo
                  qai(i)=qai(i)+2.*qaii/xi(i)/xi(i)
                  qsq(i)=qsq(i)+2.*qsqi/xi(i)/xi(i)
                  if(rp(i).ne.0.) then
                     call vctran(anpt,1,nodr(i),np,ekt(-notd,i),
     1                       drott(-nod,0,i),amnlt(1,1,1,i),nod,notd)
                     nptrn=npblk
                  else
                     nptrn=min(nblk(i),npblk)
                  endif
c
c add to the T matrix
c
                  do n=1,nptrn
                     do ip=1,2
                        tc(ip,n,ipp,mnp)=tc(ip,n,ipp,mnp)+anpt(ip,n)
                     enddo
                  enddo
c
                  if(nptrn.le.npblk)
     1                qei(i)=qei(i)+anpt(ipp,mnp)*2./xi(i)/xi(i)
75             enddo
               dqe=dqe+tc(ipp,mnp,ipp,mnp)*2./xv/xv
80          enddo
         enddo
c
c calculation of convergence
c
         qstij=0.
         do i=1,npart
            qstij=qstij+(qei(i)-qai(i))*xi(i)*xi(i)/xv/xv
         enddo
         qst=0.
         do n=1,np
            nn1=n*(n+1)
            do l=1,n
               ll1=l*(l+1)
               if(n.eq.l) then
                  is=1
               else
                  is=2
               endif
               do m=-n,n
                  mn=nn1+m
                  do ip=1,2
                     do k=-l,l
                        kl=ll1+k
                        do iq=1,2
                           qst=qst+is*tc(iq,kl,ip,mn)
     1                       *conjg(tc(iq,kl,ip,mn))
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
         qst=2.*qst/xv/xv
         npartcon=0
         errimx=0.
         errimn=100.
         do i=1,npart
            if(icon(i)) then
               npartcon=npartcon+1
            else
               erri=abs((qsq(i)-qsq0(i))/qsq(i))
               if(erri.lt.qeps2.and.np.ge.nodrmax) then
                  icon(i)=.true.
                  npartcon=npartcon+1
               endif
               errimx=max(errimx,erri)
               errimn=min(errimn,erri)
               qsq0(i)=qsq(i)
            endif
         enddo
         qserr=abs(qstij-qst)/abs(qstij)
         if((npartcon.eq.npart.or.
     1       qserr.lt.qeps2).and.np.ge.nodrmax) then
            nodrtmax=np
C            write(*,'('' converged'')')
            goto 90
         endif
         itermax=max(itermaxn,itermax)
         qe=qe+dqe
C         write(*,'(2i4,5e13.5,i5)')
C     1            np,itermaxn,dqe,qe,qstij,qst,qserr,npartcon
c
c end of column order loop
c
      enddo
c
c calculate symmetric components of T matrix
c
90    niter=itermax
      do l=1,nodrtmax
         ll1=l*(l+1)
         do k=-l,l
            klp=ll1+k
            klm=ll1-k
            do n=l+1,nodrtmax
               nn1=n*(n+1)
               do m=-n,n
                  mnp=nn1+m
                  mnm=nn1-m
                  ikm=(-1)**(m+k)
                  do iq=1,2
                     do ip=1,2
                        tc(ip,mnp,iq,klp)=ikm*tc(iq,klm,ip,mnm)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
c
      return
      end
c
c iteration solver
c rlx=0: conjugate gradient
c 0<rlx<1: underrelaxation
c rlx=1: standard gauss-seidel iteration
c rlx>1: overrelaxation (does not work)
c
c usually rlx=.9 works fine.  If convergence is not obtained,
c use rlx=0.
c Thanks to Piotr Flatau
c
      subroutine itersoln(npart,nodr,nblk,icon,eps,niter,rlx,ek,drot,
     1                    amnl,an1,pnp,anp,iter,err)
      implicit real*8(a-h,o-z)
      include 'tmnsdim.f'

      parameter(nbd=nod*(nod+2),
     1          nbtd=notd*(notd+2),nrd=.5*(npd-1)*(npd-2)+npd-1)
      integer nodr(npd),nblk(npd)
      logical icon(npd)
      real*8 drot(-nod:nod,0:nbd,nrd)
      complex*16 anpt(2,nbtd),a,amnl(2,nod,nbd,nrd),an1(2,nod,npd),
     1        pnp(2,nbd,npd),ek(-nod:nod,nrd),anp(2,nbd,npd)
      complex*16 anptc(2,nbd),
     1        cr(2,nbd,npd),cp(2,nbd,npd),cw(2,nbd,npd),cq(2,nbd,npd),
     1        cap(2,nbd,npd),caw(2,nbd,npd),cak,csk,cbk,csk2
c
      err=0.
      iter=0
      enorm=0.
c
      do i=1,npart
         do n=1,nblk(i)
            enorm=enorm+pnp(1,n,i)*conjg(pnp(1,n,i))
     1                 +pnp(2,n,i)*conjg(pnp(2,n,i))
         enddo
      enddo
c
      if(enorm.eq.0.) return
      if(rlx.ne.0) goto 200
      csk=(0.,0.)
c
      do i=1,npart
         if(icon(i)) goto 30
         do n=1,nblk(i)
            cr(1,n,i)=0.
            cr(2,n,i)=0.
         enddo
         do j=1,npart
            if((j.ne.i).and.(.not.icon(j))) then
               if(i.lt.j) then
                  ij=.5*(j-1)*(j-2)+j-i
                  idir=1
               else
                  ij=.5*(i-1)*(i-2)+i-j
                  idir=2
               endif
               do n=1,nblk(j)
                  anpt(1,n)=anp(1,n,j)
                  anpt(2,n)=anp(2,n,j)
               enddo
               call vctran(anpt,idir,nodr(j),nodr(i),ek(-nod,ij),
     1              drot(-nod,0,ij),amnl(1,1,1,ij),nod,nod)
               do n=1,nblk(i)
                  cr(1,n,i)=cr(1,n,i)+anpt(1,n)
                  cr(2,n,i)=cr(2,n,i)+anpt(2,n)
               enddo
            endif
         enddo
         do n=1,nodr(i)
            nn1=n*(n+1)
            do m=-n,n
               mn=nn1+m
               do ip=1,2
                  cr(ip,mn,i)=an1(ip,n,i)*cr(ip,mn,i)+anp(ip,mn,i)
               enddo
            enddo
         enddo
         do n=1,nblk(i)
            do ip=1,2
               cr(ip,n,i)=pnp(ip,n,i)-cr(ip,n,i)
               cq(ip,n,i)=conjg(cr(ip,n,i))
               cw(ip,n,i)=cq(ip,n,i)
               cp(ip,n,i)=cr(ip,n,i)
               csk=csk+cr(ip,n,i)*cr(ip,n,i)
            enddo
         enddo
  30  enddo
      if(cdabs(csk).eq.0.) return
  40  continue
c     enorm=0.
      cak=(0.,0.)
      do i=1,npart
         if(icon(i)) goto 50
         do n=1,nblk(i)
            do ip=1,2
               caw(ip,n,i)=0.
               cap(ip,n,i)=0.
            enddo
         enddo
         do j=1,npart
            if((j.ne.i).and.(.not.icon(j))) then
               if(i.lt.j) then
                  ij=.5*(j-1)*(j-2)+j-i
                  idir=1
               else
                  ij=.5*(i-1)*(i-2)+i-j
                  idir=2
               endif
               in=1
               do n=1,nodr(j)
                  in=-in
                  ik=-in
                  nn1=n*(n+1)
                  do m=-n,n
                     ik=-ik
                     mn=nn1+m
                     mn1=nn1-m
                     do ip=1,2
                        anptc(ip,mn)=an1(ip,n,j)
     1                    *conjg(cw(ip,mn,j))
                        anpt(ip,mn)=cp(ip,mn,j)
                     enddo
                  enddo
               enddo
               call vctran(anpt,idir,nodr(j),nodr(i),ek(-nod,ij),
     1              drot(-nod,0,ij),amnl(1,1,1,ij),nod,nod)
               call vctran(anptc,5-idir,nodr(j),nodr(i),ek(-nod,ij),
     1              drot(-nod,0,ij),amnl(1,1,1,ij),nod,nod)
               in=1
               do n=1,nodr(j)
                  in=-in
                  ik=-in
                  nn1=n*(n+1)
                  do m=-n,n
                     ik=-ik
                     mn=nn1+m
                     mn1=nn1-m
                     cap(1,mn,i)=cap(1,mn,i)+anpt(1,mn)
                     cap(2,mn,i)=cap(2,mn,i)+anpt(2,mn)
                     caw(1,mn,i)=caw(1,mn,i)+anptc(1,mn)
                     caw(2,mn,i)=caw(2,mn,i)+anptc(2,mn)
                  enddo
               enddo
            endif
         enddo
         do n=1,nodr(i)
            nn1=n*(n+1)
            do m=-n,n
               mn=nn1+m
               do ip=1,2
                  caw(ip,mn,i)=conjg(caw(ip,mn,i))+cw(ip,mn,i)
                  cap(ip,mn,i)=an1(ip,n,i)*cap(ip,mn,i)+cp(ip,mn,i)
               enddo
            enddo
         enddo
         do n=1,nblk(i)
            do ip=1,2
               cak=cak+cap(ip,n,i)*conjg(cw(ip,n,i))
            enddo
         enddo
50    enddo
      cak=csk/cak
      csk2=(0.,0.)
      err=0.
      do i=1,npart
         if(.not.icon(i)) then
         do n=1,nblk(i)
            do ip=1,2
               anp(ip,n,i)=anp(ip,n,i)+cak*cp(ip,n,i)
c              enorm=enorm+anp(ip,n,i)*conjg(anp(ip,n,i))
               cr(ip,n,i)=cr(ip,n,i)-cak*cap(ip,n,i)
               cq(ip,n,i)=cq(ip,n,i)-conjg(cak)*caw(ip,n,i)
               csk2=csk2+cr(ip,n,i)*conjg(cq(ip,n,i))
               err=err+cr(ip,n,i)*conjg(cr(ip,n,i))
            enddo
         enddo
         endif
      enddo
      err=err/enorm
      if(err.lt. eps) return
      cbk=csk2/csk
      do i=1,npart
         if(.not.icon(i)) then
         do n=1,nblk(i)
            do ip=1,2
               cp(ip,n,i)=cr(ip,n,i)+cbk*cp(ip,n,i)
               cw(ip,n,i)=cq(ip,n,i)+conjg(cbk)*cw(ip,n,i)
            enddo
         enddo
         endif
      enddo
      csk=csk2
      iter=iter+1
      if(iter.le.niter) goto 40
      return
c
200   err=0.
      do i=1,npart
         if(icon(i)) goto 250
         do n=1,nblk(i)
            do ip=1,2
               cr(ip,n,i)=0.
            enddo
         enddo
         do j=1,npart
            if(i.eq.j.or.icon(j)) goto 240
            if(i.lt.j) then
               ij=.5*(j-1)*(j-2)+j-i
               idir=1
            else
               ij=.5*(i-1)*(i-2)+i-j
               idir=2
            endif
            do n=1,nblk(j)
               do ip=1,2
                  anpt(ip,n)=anp(ip,n,j)
               enddo
            enddo
            call vctran(anpt,idir,nodr(j),nodr(i),ek(-nod,ij),
     1           drot(-nod,0,ij),amnl(1,1,1,ij),nod,nod)
            do n=1,nblk(i)
               cr(1,n,i)=cr(1,n,i)+anpt(1,n)
               cr(2,n,i)=cr(2,n,i)+anpt(2,n)
            enddo
240      enddo
         do n=1,nodr(i)
            nn1=n*(n+1)
            do m=-n,n
               mn=nn1+m
               do ip=1,2
                  a=pnp(ip,mn,i)-cr(ip,mn,i)*an1(ip,n,i)
     1             -anp(ip,mn,i)
                  err=err+a*conjg(a)
                  anp(ip,mn,i)=anp(ip,mn,i)+a*rlx
c                 enorm=enorm+anp(ip,mn,i)*conjg(anp(ip,mn,i))
               enddo
            enddo
         enddo
250   enddo
      err=err/enorm
      iter=iter+1
      if((err.gt.eps).and.(iter.lt.niter)) goto 200
      return
      end
c
c
c this performs a vector harmonic expansion coefficient translation
c from one origin to the next.   Uses a rotation-translation-rotation
c formulation.
c
c idir=1: a2=h21 a1
c idir=2: a1=h12 a2
c idir=3: a2 h21=a1
c idir=4: a1 h12=a2
c
      subroutine vctran(anpt,idir,nodrj,nodri,ek,drot,amnl,ndd,nda)
      implicit real*8(a-h,o-z)
      include 'tmnsdim.f'
      parameter(nbtd=notd*(notd+2))
      real*8 drot(-ndd:ndd,0:*)
      complex*16 anpt(2,*),ant(2,notd),amt(2,-notd:notd),a,b,
     1           amnl(2,nda,*),ek(-nda:nda)
c
      if(idir.eq.3.or.idir.eq.4) then
         do n=1,nodrj
            nn1=n*(n+1)
            im=1
            do m=1,n
               im=-im
               mn=nn1+m
               mnm=nn1-m
               a=anpt(1,mn)
               anpt(1,mn)=im*anpt(1,mnm)
               anpt(1,mnm)=im*a
               a=anpt(2,mn)
               anpt(2,mn)=im*anpt(2,mnm)
               anpt(2,mnm)=im*a
            enddo
         enddo
      endif

      do n=1,nodrj
         nn1=n*(n+1)
         mmax=min(n,nodri)
         do m=-mmax,mmax
            amt(1,m)=0.
            amt(2,m)=0.
         enddo
         do k=-n,n
            knm=nn1-k
            kn=nn1+k
            a=ek(k)*anpt(1,kn)
            b=ek(k)*anpt(2,kn)
            do m=-mmax,mmax
               amt(1,m)=amt(1,m)+a*drot(-m,knm)
               amt(2,m)=amt(2,m)+b*drot(-m,knm)
            enddo
         enddo
         do m=-mmax,mmax
            mn=nn1+m
            anpt(1,mn)=amt(1,m)
            anpt(2,mn)=amt(2,m)
         enddo
      enddo
c
      mmax=min(nodrj,nodri)
      if(idir.eq.1.or.idir.eq.4) then
         do m=-mmax,mmax
            n1=max(1,abs(m))
            do n=n1,nodrj
               mn=n*(n+1)+m
               do ip=1,2
                  ant(ip,n)=anpt(ip,mn)
               enddo
            enddo
            do n=n1,nodri
               mn=n*(n+1)+m
               a=0.
               b=0.
               do l=n1,nodrj
                  ml=l*(l+1)+m
                  a=a+amnl(1,n,ml)*ant(1,l)
     1               +amnl(2,n,ml)*ant(2,l)
                  b=b+amnl(1,n,ml)*ant(2,l)
     1               +amnl(2,n,ml)*ant(1,l)
               enddo
               anpt(1,mn) = a
               anpt(2,mn) = b
            enddo
         enddo
      else
         do m=-mmax,mmax
            n1=max(1,abs(m))
            do n=n1,nodrj
               mn=n*(n+1)+m
               do ip=1,2
                  ant(ip,n)=anpt(ip,mn)
               enddo
            enddo
            do n=n1,nodri
               mn=n*(n+1)+m
               mnm=n*(n+1)-m
               a=0.
               b=0.
               do l=n1,nodrj
                  a=a+amnl(1,l,mnm)*ant(1,l)
     1               +amnl(2,l,mnm)*ant(2,l)
                  b=b+amnl(1,l,mnm)*ant(2,l)
     1               +amnl(2,l,mnm)*ant(1,l)
               enddo
               anpt(1,mn) = a
               anpt(2,mn) = b
            enddo
         enddo
      endif
c
      in=1
      do n=1,nodri
         in=-in
         nn1=n*(n+1)
         kmax=min(n,nodrj)
         do m=-n,n
            amt(1,m)=0.
            amt(2,m)=0.
         enddo
         ik=-(-1)**kmax
         do k=-kmax,kmax
            ik=-ik
            a=ik*anpt(1,nn1+k)
            b=ik*anpt(2,nn1+k)
            do m=-n,n
               amt(1,m)=amt(1,m)+a*drot(k,nn1+m)
               amt(2,m)=amt(2,m)+b*drot(k,nn1+m)
            enddo
         enddo
         ik=-in
         do m=-n,n
            ik=-ik
            mn=nn1+m
            anpt(1,mn)=amt(1,m)*ek(-m)*ik
            anpt(2,mn)=amt(2,m)*ek(-m)*ik
         enddo
      enddo
      if(idir.eq.3.or.idir.eq.4) then
         do n=1,nodri
            nn1=n*(n+1)
            im=1
            do m=1,n
               im=-im
               mn=nn1+m
               mnm=nn1-m
               a=anpt(1,mn)
               anpt(1,mn)=im*anpt(1,mnm)
               anpt(1,mnm)=im*a
               a=anpt(2,mn)
               anpt(2,mn)=im*anpt(2,mnm)
               anpt(2,mnm)=im*a
            enddo
         enddo
      endif
      return
      end
c
c  dc is the normalized generalized spherical function
c  dc(k,n*(n+1)+m) = ((n-k)!(n+m)!/(n+k)!/(n-m)!)^(1/2) D^k_{mn},
c  where D^k_{mn} is defined in M&M JOSA 96.
c  This routine calculates dc(k,n*(n+1)+m) for k=-nmax..nmax and
c  n=0...lmax.
c
      subroutine rotcoef(cbe,nmax,lmax,dc,ndim)
      include 'tmnsdim.f'
      parameter(nbtd=notd*(notd+2),notd2=notd+notd,nbc=3*notd2+6)
      implicit real*8(a-h,o-z)
      real*8 dc(-ndim:ndim,0:*)
      real*8 dk0(-notd2:notd2),dk01(-notd2:notd2)
      common/consts/bcof(0:nbc,0:nbc),fnr(0:2*nbc)
      data ci/(0.d0,1.d0)/
c
c calculation of rotation functions.
c D^k_{mn}(cos beta) = dc(k,n*(n+1)+m)
c
      sbe=dsqrt((1.d0+cbe)*(1.d0-cbe))
      cbe2=.5d0*(1.d0+cbe)
      sbe2=.5d0*(1.d0-cbe)
      in=1
      dk0(0)=1.d0
      sben=1.d0
      dc(0,0)=1.d0
      dk01(0)=0.
      do n=1,lmax
         kmax=min(n,nmax)
         nn1=n*(n+1)
         in=-in
         sben=sben*sbe/2.d0
         dk0(n)=in*sben*bcof(n,n)
         dk0(-n)=in*dk0(n)
         dk01(n)=0.
         dk01(-n)=0.
         dc(0,nn1+n)=dk0(n)
         dc(0,nn1-n)=dk0(-n)
         do k=-n+1,n-1
            kn=nn1+k
            dkt=dk01(k)
            dk01(k)=dk0(k)
            dk0(k)=(cbe*dble(n+n-1)*dk01(k)-fnr(n-k-1)*fnr(n+k-1)*dkt)
     1             /(fnr(n+k)*fnr(n-k))
            dc(0,kn)=dk0(k)
         enddo
         im=1
         do m=1,kmax
            im=-im
            fmn=1.d0/fnr(n-m+1)/fnr(n+m)
            m1=m-1
            dkm0=0.
            do k=-n,n
               kn=nn1+k
               dkm1=dkm0
               dkm0=dc(m1,kn)
               if(k.eq.n) then
                  dkn1=0.
               else
                  dkn1=dc(m1,kn+1)
               endif
               dc(m,kn)=(fnr(n+k)*fnr(n-k+1)*cbe2*dkm1
     1           -fnr(n-k)*fnr(n+k+1)*sbe2*dkn1
     1              -dble(k)*sbe*dc(m1,kn))*fmn
               dc(-m,nn1-k)=dc(m,kn)*(-1)**(k)*im
            enddo
         enddo
      enddo
      return
      end
c
c this computes the normalized translation coefficients for an
c axial translation of distance r.  For itype=1 or 3, the translation
c uses the spherical Bessel or Hankel functions as a basis function,
c respectively.    They are related to the coefficients appearing in
c M&M JOSA 96 by
c
c J^{ij}_{mnp mlq} = (E_{ml}/E_{mn})^(1/2) ac(s,n,l*(l+1)+m)
c
c where
c
c   E_{mn} = n(n+1)(n+m)!/((2n+1)(n-m)!)
c   s=mod(p+q,2)+1 (i.e., s=1 for the A coefficient, =2 for the B
c   coefficient)
c
c  ac(2,nd,*) is a complex*16 array.  ac is calculated for n=1,nmax
c  and l=1,lmax.
c
c  The calculation procedure is based on the recent derivation
c  of the addition theorem for vector harmonics, appearing in
c  Fuller and Mackowski, proc. Light Scattering by Nonspherical
c  Particles, NASA/GISS Sept. 1998.
c
      subroutine trancoef(itype,r,lmax,nmax,ac,nd)
      implicit real*8 (a-h,o-z)
      include 'tmnsdim.f'
      parameter(nbtd=notd*(notd+2),notd2=notd+notd+1,
     1         nfd=notd2*(notd2+2),nbc=6*notd+6)
      real*8 vc1(0:notd2+1),vc2(0:notd2+1)
      complex*16 ci,xi(0:notd2)
      complex*16 a,b,c,ac(2,nd,*)
      common/consts/bcof(0:nbc,0:nbc),fnr(0:2*nbc)
      data ci/(0.d0,1.d0)/
      if(r.eq.0.) then
         lblk=lmax*(lmax+2)
         do n=1,nmax
            do l=1,lblk
               do ip=1,2
                  ac(ip,n,l)=0.
               enddo
            enddo
         enddo
         if(itype.ne.1) return
         do n=1,min(nmax,lmax)
            do l=n*(n+1)-n,n*(n+1)+n
               ac(1,n,l)=1.
            enddo
         enddo
         return
      endif
      call bessel(nmax+lmax+1,r,xi)
5     nlmax=max(nmax,lmax)
      do n=0,nmax+lmax+1
         if(itype.eq.1) then
            xi(n)=dble(xi(n))/r*ci**n
         else
            xi(n)=xi(n)/r*ci**n
         endif
      enddo
      do n=1,nmax
         n21=n+n+1
         do l=1,lmax
            c=fnr(n21)*fnr(l+l+1)*ci**(n-l)
            ll1=l*(l+1)
            call vcfunc(-1,n,1,l,n+l,vc2)
            iwmn=abs(n-l)
            iwmx=n+l
            nlmin=min(l,n)
            do m=-nlmin,nlmin
               a=0.
               b=0.
               call vcfunc(-m,n,m,l,n+l,vc1)
               do iw=iwmn,iwmx
                  alnw=vc1(iw)*vc2(iw)
                  if(mod(n+l+iw,2).eq.0) then
                     a=a+alnw*xi(iw)
                  else
                     b=b+alnw*xi(iw)
                  endif
               enddo
               ac(1,n,ll1+m)=-c*a*(-1)**m
               ac(2,n,ll1+m)=-c*b*(-1)**m
            enddo
         enddo
      enddo
      return
      end
c
c vector coupling coefficients vc(iw) = C(m,n|k,l|m+k,iw)
c uses an upwards recurrence
c
      subroutine vcfunc(m,n,k,l,wmax,vcn)
      include 'tmnsdim.f'
      parameter(nbtd=notd*(notd+2),nbc=6*notd+6)
      implicit real*8(a-h,o-z)
      real*8 vcn(0:*)
      integer w,wmax,w1,w2
      common/consts/bcof(0:nbc,0:nbc),fnr(0:2*nbc)
      mk=abs(m+k)
      nl=abs(n-l)
      if(nl.ge.mk) then
         w=nl
         if(n.ge.l) then
            m1=m
            n1=n
            l1=l
            k1=k
         else
            m1=k
            n1=l
            k1=m
            l1=n
         endif
         vc1=(-1)**(k1+l1)*bcof(l1+k1,w-m1-k1)
     1     *bcof(l1-k1,w+m1+k1)/bcof(l1+l1,w+w+1)
      else
         w=mk
         if(m+k.ge.0) then
            vc1=(-1)**(n+m)*bcof(n-l+w,l-k)*bcof(l-n+w,n-m)
     1          /bcof(w+w+1,n+l-w)
         else
            vc1=(-1)**(l+k)*bcof(n-l+w,l+k)*bcof(l-n+w,n+m)
     1          /bcof(w+w+1,n+l-w)
         endif
      endif
      w1=w
      vcn(w)=vc1
      w=w1+1
      mk=m+k
      w2=min(wmax,n+l)
      if(w2.gt.w1) then
         t1=2*w*fnr(w+w+1)*fnr(w+w-1)/(fnr(w+mk)*fnr(w-mk)
     1     *fnr(n-l+w)*fnr(l-n+w)*fnr(n+l-w+1)*fnr(n+l+w+1))
         if(w1.eq.0) then
            t2=.5*dble(m-k)
         else
            t2=dble((m-k)*w*(w-1)-mk*n*(n+1)+mk*l*(l+1))
     1       /dble(2*w*(w-1))
         endif
         vcn(w)=t1*t2*vcn(w1)
      endif
      do w=w1+2,w2
         t1=2*w*fnr(w+w+1)*fnr(w+w-1)/(fnr(w+mk)*fnr(w-mk)
     1     *fnr(n-l+w)*fnr(l-n+w)*fnr(n+l-w+1)*fnr(n+l+w+1))
         t2=dble((m-k)*w*(w-1)-mk*n*(n+1)+mk*l*(l+1))
     1    /dble(2*w*(w-1))
         t3=fnr(w-mk-1)*fnr(w+mk-1)*fnr(l-n+w-1)*fnr(n-l+w-1)
     1     *fnr(n+l-w+2)*fnr(n+l+w)/(dble(2*(w-1))*fnr(2*w-3)
     1     *fnr(2*w-1))
         vcn(w)=t1*(t2*vcn(w-1)-t3*vcn(w-2))
      enddo
      return
      end
      subroutine bessel(n,ds,xi)
      implicit real*8(a-h,o-z)
      complex*16 xi(0:*)
c
c computes riccatti-bessel function xi(0:n)=psi(n)+i*chi(n)
c for real argument ds
c
c check if recurrence is forward or backwards
c
      if(int(ds).ge.n) goto 40
c
c backwards recurrence for log. derivative dn1(n)=psi'(n)/psi(n)
c
      ns=nint(ds+4.*(ds**.3333)+17)
      dns=0.d0
      do 20 i=ns-1,n,-1
         sn=dble(i+1)/ds
         dns=sn-1.d0/(dns+sn)
   20 continue
      xi(n)=dns
      xi(n-1)=dble(n)/ds-1.d0/(dns+dble(n)/ds)
      do 25 i=n-2,1,-1
         sn=dble(i+1)/ds
         xi(i)=sn-1.d0/(xi(i+1)+sn)
   25 continue
c
c forward recurrence: psi(n) computed from dn1(n)
c
      chi0=-dcos(ds)
      psi=dsin(ds)
      chi1=chi0/ds-psi
      xi(0)=dcmplx(psi,chi0)
      do 30 i=1,n
         chi2=dble(i+i+1)/ds*chi1-chi0
         psi=psi/(dble(i)/ds+xi(i))
         xi(i)=dcmplx(psi,chi1)
         chi0=chi1
         chi1=chi2
   30 continue
      return
c
c forward recurrence: applicable only if ds > n
c
   40 chi0=-dcos(ds)
      psi0=dsin(ds)
      chi1=chi0/ds-psi0
      psi1=psi0/ds+chi0
      xi(0)=dcmplx(psi0,chi0)
      do 50 i=1,n
         sn=dble(i+i+1)/ds
         chi2=sn*chi1-chi0
         psi2=sn*psi1-psi0
         xi(i)=dcmplx(psi1,chi1)
         chi0=chi1
         chi1=chi2
         psi0=psi1
         psi1=psi2
   50 continue
      return
      end
c
c single-sphere lorenz/mie coefficients
c
      subroutine mie1(x,sn,sk,nstop,qeps,qext,qsca,an,cn,iray)
      implicit real*8 (a-h,o-z)
      include 'tmnsdim.f'
      parameter(nomax2=2*nod)
      complex*16 y,ri,xip,pcp,da,db,pc(0:nomax2),xi(0:nomax2),
     1 an(2,nod),na,nb,cn(2,nod)
c
      if(qeps.gt.0.) nstop=nint(x+4.*x**(1./3.))+5.
      nstop=min(nstop,nod)
      ri=cmplx(sn,sk)
      y=x*ri
      call psic(nstop,y,pc)
      call bessel(nstop,x,xi)
      qsca=0.0
      qext=0.0
      do 300 n=1,nstop
         prn=dble(xi(n))
         pcp=pc(n-1)-n*pc(n)/y
         xip=xi(n-1)-n*xi(n)/x
         prp=dble(xip)
         da=ri*xip*pc(n)-xi(n)*pcp
         db=ri*xi(n)*pcp-xip*pc(n)
         na=ri*prp*pc(n)-prn*pcp
         nb=ri*prn*pcp-prp*pc(n)
         an(1,n)=na/da
         an(2,n)=nb/db
         if(na.eq.0.) then
            dn1=0.
         else
            dn1=1./na/dconjg(na)
         endif
         if(nb.eq.0.) then
            cn1=0.
         else
            cn1=1./nb/dconjg(nb)
         endif
         na=dcmplx(0.,1.)*pcp*dconjg(pc(n))
         cn(1,n)=na*dn1*dconjg(ri)
         cn(2,n)=na*cn1*ri
c        cn(1,n)=dcmplx(0.d0,1.d0)*ri/na
c        cn(2,n)=-dcmplx(0.d0,1.d0)*ri/nb
         qsca=qsca+(n+n+1)*(cdabs(an(1,n))*cdabs(an(1,n))
     1        +cdabs(an(2,n))*cdabs(an(2,n)))
         qext1=(n+n+1)*dble(an(1,n)+an(2,n))
         qext=qext+qext1
         err=abs(qext1)/abs(qext)
         if(err.lt.qeps.or.n.eq.nstop) goto 310
  300 continue
  310 nstop=min(n,nod)
      qsca=2./x/x*qsca
      qext=2./x/x*qext
      if((nstop.eq.1).and.(iray.gt.0)) then
         an(1,1)=-dcmplx(0.,1.)*2.*x*x*x*(ri-1.)*(ri+1.)/(ri*ri+2.)/3.
         if(iray.eq.2) an(1,1)=an(1,1)/(1.+an(1,1))
      endif
      return
      end
c
c spherical bessel function of a complex argument
c
      subroutine psic(n,ds,psi)
      complex*16 ds,psi(0:*),sn,psins
c
      ns=nint(cdabs(ds)+4.*(cdabs(ds)**.3333)+17)
      psins=(0.d0,0.d0)
      do 20 i=ns-1,n,-1
         sn=dble(i+1)/ds
         psins=sn-1.d0/(psins+sn)
   20 continue
      psi(n)=psins
      psi(n-1)=n/ds-1.d0/(psins+n/ds)
      do 25 i=n-2,1,-1
         sn=dble(i+1)/ds
         psi(i)=sn-1.d0/(psi(i+1)+sn)
   25 continue
      psins=cdsin(ds)
      psi(0)=psins
      do 30 i=1,n
         psins=psins/(dble(i)/ds+psi(i))
         psi(i)=psins
   30 continue
      return
      end
