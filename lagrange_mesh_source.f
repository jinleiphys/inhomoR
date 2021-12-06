      module lagrange_mesh_source
      real*8,allocatable::  g(:),tc(:,:)
      real*8,allocatable ::  rg(:),wg(:)
      real*8,dimension(:),allocatable :: fc,dfc,gc,dgc

      contains
      subroutine rmat_ini(ng,rmax)
c ng=number of Lagrange functions
c rmax=channel radius
c rg(1:ng)= array with the zeros of Legendre polynomilas (in [0,1])
c wg(1:ng)= weigths
      implicit real*8(a-h,o-z)
      if(allocated(g))deallocate(g,tc)
      if(allocated(rg)) deallocate(rg,wg)
      allocate(g(ng),tc(ng,ng))
      allocate(rg(ng),wg(ng))

      call legzo(ng,rg,wg) ! I think this can be replaced by gauleg, need to check

      !compute the  fn(a) in sofia's notes
      g(1:ng)=1/sqrt(rmax*rg(1:ng)*(1-rg(1:ng)))
      g(1:ng:2)=-g(1:ng:2) ! used for the (-1)^n factor

      do 1 i=1,ng
      ui=rg(i)
      do 1 j=1,ng
      uj=rg(j)
      if(i.eq.j)then
       t=4*ng*(ng+1)+3+(1-6*ui)/ui/(1-ui)
       t=t/(3*ui*(1-ui))   ! compute <fn|T+Lc(0)|fn> without hbarc^2 / (2 mu a^2 )
      else
       t=ng*(ng+1)+1+(ui+uj-2*ui*uj)/(ui-uj)**2
       t=t-1/(1-ui)-1/(1-uj)
       t=t/sqrt(ui*uj*(1-ui)*(1-uj))
       if(mod(i+j,2).ne.0)t=-t   ! compute <fn |T+Lc(0)|fn' > without hbarc^2 / (2 mu a^2 )
      end if
      tc(i,j)=t/rmax**2   ! <fn |T+Lc(0)| fn'> without hbarc^2 / (2 mu ) factor
    1 continue
      return
      end SUBROUTINE

      subroutine rmat_inho(ng,rmax,cpot,csou,ecm,eta,hm,l,cut,cft)
C     subroutine rmat_inho(ng,rmax,cpot,csou,ecm,eta,hm,l,cu,cut,cf,cft)
c cpot(1:ng)=complex arrray with the potentials at the mesh points
c csou(1:ng)=complex array with the source term
c ecm,eta,l: obvious
c hm=hbar**2/(2*mu)
cu,cut=scattering matrices without and with the source term
cf(1:ng),cft(1:ng)=complex arrays with the wave function at the mesh points
      implicit real*8(a,b,d-h,o-z)
      implicit complex*16(c)
      dimension c(ng,ng),crho(ng), cmat(ng,ng)
      dimension cf(ng),cft(ng),cpot(ng),csou(ng) !,fc(100),dfc(100),
!     1 gc(100),dgc(100)                         ! fc dfc, gc, dgc part can be changed
      dimension cmat_rho(ng), cmat_phi(ng) ! use for calculate C^{-1} \rho and C^{-1} \phi without invert matrix
      integer :: INFO,NRHS
      integer,dimension(1:ng) :: IPIV
c      allocate (c(ng,ng),crho(ng))
      c=hm*tc  ! <fn |T+Lc(0)| fn'>
      fac=hm*l*(l+1) ! l(l+1) hbarc^2 / (2 mu)  ! used for centrifugal_barrier
      xk=sqrt(ecm/hm)  ! k
      ak=xk*rmax    ! \rho=k*r  use for boundary condition

      do 1 i=1,ng
      r=rg(i)*rmax
      c(i,i)=c(i,i)+cpot(i)+fac/r**2-ecm   ! C-EI matrix in Sofia's notes
      crho(i)=csou(i)*sqrt(wg(i)*rmax) !  souce term function in Lagrange mesh basis < \phi_i | \rho> ! See details in pierre's notes
    1 continue

c      call cminv(c,ng,ng) ! invert the C matrix ! this can be replaced !!!!!!!!!!!!!!!!!!!


      NRHS=1
      cmat=c
      cmat_phi=g
      call ZGESV( ng, NRHS, cmat, ng, IPIV, cmat_phi, ng, INFO )  ! used for computed C^{-1} \phi
      If(INFO/=0) stop "error in calling ZGESV"

      NRHS=1
      cmat=c
      cmat_rho=crho
      call ZGESV( ng, NRHS, cmat, ng, IPIV, cmat_rho, ng, INFO )  ! used for computed C^{-1} \rho
      If(INFO/=0) stop "error in calling ZGESV"



      cr=0
      cr2=0
c      do 2 i1=1,ng
c      do 2 i2=1,ng
c      cr=cr+g(i1)*g(i2)*c(i1,i2)   ! R-marix without hbarc^2/ (2 mu a) factor
c      cr2=cr2+crho(i1)*g(i2)*c(i1,i2)  ! for the souce term in lagrange mesh
c    2 continue


      do i1=1, ng
        cr=cr+ g(i1)*cmat_phi(i1)
        cr2=cr2 + g(i1)* cmat_rho(i1)
      end do
      cr=cr*hm/rmax ! R-matrix

C      call coufra(ak,eta,l,l,fc,dfc,gc,dgc) ! call the coulomb function;  this can be replaced !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ci=dcmplx(gc(l),-fc(l))
      cip=dcmplx(dgc(l),-dfc(l))
      co=conjg(ci)
      cop=conjg(cip)

c      cu=(ci-ak*cip*cr)/(co-ak*cop*cr) ! S-matrix for homogeneous equation
      cut=-cr2/(co-ak*cop*cr)             ! S-matrix for inhomogeneous equation

c      crho(1:ng)=crho(1:ng)-cut*xk*cop*hm*g(1:ng)
c      cft(1:ng)=matmul(c(1:ng,1:ng),crho(1:ng))
      cft(1:ng) =  cmat_rho(1:ng) - cut*xk*cop*hm*cmat_phi(1:ng)
      cft(1:ng)=cft(1:ng)/sqrt(rmax*wg(1:ng))  ! compute the wave function of inhomogeneous equation

c      cf(1:ng)=matmul(c(1:ng,1:ng),g(1:ng))
c      cf(1:ng)=cf(1:ng)*hm*xk*(ci-cu*co)            ! compute the wave function of homogeneous equation
      return
      end SUBROUTINE

        SUBROUTINE LEGZO(N,X,W)
C
C       =========================================================
C       Purpose : Compute the zeros of Legendre polynomial Pn(x)
C                 in the interval [0,1], and the corresponding
C                 weighting coefficients for Gauss-Legendre
C                 integration
C       Input :   n    --- Order of the Legendre polynomial
C       Output:   X(n) --- Zeros of the Legendre polynomial
C                 W(n) --- Corresponding weighting coefficients
C       =========================================================
C
C Author: J. M. Jin
C Downloaded from http://jin.ece.illinois.edu/routines/routines.html
        implicit real*8 (a-h,o-z)
        dimension x(n),w(n)
        data pi/3.1415926535898d0/,one/1/
        n0=(n+1)/2
        do 45 nr=1,n0
           z=cos(pi*(nr-0.25d0)/(n+0.5d0))
10         z0=z
           p=1
           do 15 i=1,nr-1
15            p=p*(z-x(i))
           f0=1
           if (nr.eq.n0.and.n.ne.2*int(n/2)) z=0
           f1=z
           do 20 k=2,n
              pf=(2-one/k)*z*f1-(1-one/k)*f0
              pd=k*(f1-z*pf)/(1-z*z)
              f0=f1
20            f1=pf
           if (z.eq.0) go to 40
           fd=pf/p
           q=0
           do 35 i=1,nr-1
              wp=1
              do 30 j=1,nr-1
                 if (j.ne.i) wp=wp*(z-x(j))
30            continue
35            q=q+wp
           gd=(pd-q*fd)/p
           z=z-fd/gd
           if (abs(z-z0).gt.abs(z)*1.0d-15) go to 10
40         x(nr)=z
           x(n+1-nr)=-z
           w(nr)=2/((1-z*z)*pd*pd)
45         w(n+1-nr)=w(nr)
        x(1:n)=(1+x(n:1:-1))/2
        w(1:n)=w(n:1:-1)/2
        return
        end SUBROUTINE




       subroutine cminv(c,m1mod,mdim)
       implicit complex*16(c)
       allocatable cwork(:),inde(:)
       dimension c(mdim,mdim)
       lwork=64*mdim
       allocate(cwork(lwork),inde(mdim))
       call zsytrf('U',m1mod,c,mdim,inde,cwork,lwork,info)
       call zsytri('U',m1mod,c,mdim,inde,cwork,info2)
       do i=1,m1mod
       c(i,1:i)=c(1:i,i)
       enddo
       return
       end subroutine


       end module


*COUFRA
      SUBROUTINE COUFRA(RHO,ETA,MINL,MAXL,FC,FCP,GC,GCP)
C*** FONCTIONS COULOMBIENNES CALCULEES EN R = RHO PAR LA METHODE DES FRA
C*** CONTINUES DE STEED. MINL ET MAXL CORRESPONDENT AUX VRAIES VALEURS D
C*** VOIR BARNETT, FENG, STEED ET GOLDFARB, CPC 1974 *******************
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 K,K1,K2,K3,K4,M1,M2,M3,M4
      DIMENSION FC(MAXL),FCP(MAXL),GC(MAXL),GCP(MAXL)
      SAVE
      DATA ACCUR,STEP/1.D-7,100.0D0/
      PACE = STEP
      ACC = ACCUR
      R = RHO
      KTR = 1
      LMAX = MAXL
      LMIN1 = MINL+1
      XLL1 = MINL*LMIN1
      ETA2 = ETA**2
      TURN = ETA+SQRT(ETA2+XLL1)
      IF(R.LT.TURN.AND.ABS(ETA).GE.1.D-6) KTR = -1
      KTRP = KTR
      GO TO 2
    1 R = TURN
      TF = F
      TFP = FP
      LMAX = MINL
      KTRP = 1
    2 ETAR = ETA*R
   21 RHO2=R*R
      PL = LMAX+1
      PMX = PL+0.5D0
C** FRACTION CONTINUE POUR FP(MAXL)/F(MAXL) ; XL=F ; XLPRIME=FP ********
      FP = ETA/PL+PL/R
      DK = ETAR+ETAR
      DEL = 0
      D = 0
      F = 1
      K = (PL*PL-PL+ETAR)*(PL+PL-1)
      IF(PL*PL+PL+ETAR.NE.0.) GO TO 3
      R = R*1.0000001D0
      GO TO 2
    3 H = (PL*PL+ETA2)*(1-PL*PL)*RHO2
      K = K+DK+PL*PL*6
      D = 1/(D*H+K)
      DEL = DEL*(D*K-1)
      IF(PL.LT.PMX) DEL = -R*(PL*PL+ETA2)*(PL+1)*D/PL
      PL = PL+1
      FP = FP+DEL
      IF(D.LT.0) F = -F
      IF(PL.GT.20000.0D0) GO TO 11
      IF(ABS(DEL/FP).GE.ACC) GO TO 3
      FP = F*FP
      IF(LMAX.EQ.MINL) GO TO 5
      FC(LMAX+1) = F
      FCP(LMAX+1) = FP
C*** RECURRENCE ARRIERE POUR F ET FP ; GC,GCP UTILISES POUR STOCKAGE ***
      L = LMAX
      DO 4 LP=LMIN1,LMAX
      PL = L
      GC(L+1) = ETA/PL+PL/R
      GCP(L+1) = SQRT(ETA2+PL*PL)/PL
      FC(L) =(GC(L+1)*FC(L+1)+FCP(L+1))/GCP(L+1)
      FCP(L) = GC(L+1)*FC(L)-GCP(L+1)*FC(L+1)
    4 L = L-1
      F = FC(LMIN1)
      FP = FCP(LMIN1)
    5 IF(KTRP.EQ.-1) GO TO 1
C*** MEME CALCUL POUR R = TURN SI RHO.LT.TURN
C*** P + I.Q CALCULE EN MINL , EQUATION (32)
      P = 0
      Q = R-ETA
      PL = 0
      AR = -(ETA2+XLL1)
      AI = ETA
      BR = Q+Q
      BI = 2
      WI = ETA+ETA
      DR = BR/(BR*BR+BI*BI)
      DI = -BI/(BR*BR+BI*BI)
      DP = -(AR*DI+AI*DR)
      DQ = AR*DR-AI*DI
    6 P = P+DP
      Q = Q+DQ
      PL = PL+2
      AR = AR+PL
      AI = AI+WI
      BI = BI+2
      D = AR*DR-AI*DI+BR
      DI = AI*DR+AR*DI+BI
      T = 1/(D*D+DI*DI)
      DR = T*D
      DI = -T*DI
      H = BR*DR-BI*DI-1
      K = BI*DR+BR*DI
      T = DP*H-DQ*K
      DQ = DP*K+DQ*H
      DP = T
      IF(PL.GT.46000.0D0) GO TO 11
      IF(ABS(DP)+ABS(DQ).GE.(ABS(P)+ABS(Q))*ACC) GO TO 6
      P = P/R
      Q = Q/R
C*** CALCUL DE FP,G,GP, ET NORMALISATION DE F EN L = MINL **************
      G = (FP-P*F)/Q
      GP = P*G-Q*F
      W = 1/SQRT(FP*G-F*GP)
      G = W*G
      GP = W*GP
      IF(KTR.EQ.1) GO TO 8
      F = TF
      FP = TFP
      LMAX = MAXL
C*** CALCUL DE G(MINL) ET GP(MINL) PAR INTEGRATION RUNGE-KUTTA A PARTIR
C***         VOIR FOX ET MAYERS(1968) PG 202
      IF(RHO.LT.0.2D0*TURN) PACE = 999.0D0
      R3=1.0D0/3.0D0
      H = (RHO-TURN)/(PACE+1)
      H2 = H/2
      I2 = INT(PACE+0.001D0)
      ETAH = ETA*H
      H2LL = H2*XLL1
      S = (ETAH+H2LL/R)/R-H2
    7 RH2 = R+H2
      T = (ETAH+H2LL/RH2)/RH2-H2
      K1 = H2*GP
      M1 = S*G
      K2 = H2*(GP+M1)
      M2 = T*(G+K1)
      K3 = H*(GP+M2)
      M3 = T*(G+K2)
      M3 = M3+M3
      K4 = H2*(GP+M3)
      RH = R+H
      S = (ETAH+H2LL/RH)/RH-H2
      M4 = S*(G+K3)
      G = G+(K1+K2+K2+K3+K4)*R3
      GP = GP+(M1+M2+M2+M3+M4)*R3
      R = RH
      I2 = I2-1
      IF(ABS(GP).GT.1.D300) GO TO 11
      IF(I2.GE.0) GO TO 7
      W = 1/(FP*G-F*GP)
C*** RECURRENCE AVANT A PARTIR DE GC(MINL) ET GCP(MINL)
C*** RENORMALISATION DE FC ET FCP POUR CHAQUE VALEUR DE L **************
    8 GC(LMIN1) = G
      GCP(LMIN1) = GP
      IF(LMAX.EQ.MINL) GO TO 10
      DO 9 L=LMIN1,LMAX
      T = GC(L+1)
      GC(L+1) = (GC(L)*GC(L+1)-GCP(L))/GCP(L+1)
      GCP(L+1) = GC(L)*GCP(L+1)-GC(L+1)*T
      FC(L+1) = W*FC(L+1)
    9 FCP(L+1) = W*FCP(L+1)
      FC(LMIN1) = W*FC(LMIN1)
      FCP(LMIN1) = W*FCP(LMIN1)
      RETURN
   10 FC(LMIN1) = W*F
      FCP(LMIN1) = W*FP
      RETURN
   11 W = 0
      G = 0
      GP = 0
      GO TO 8
      END SUBROUTINE
