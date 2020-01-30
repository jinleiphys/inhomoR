C>     this version of gauleg is embedded in a module 
C>     and includes the standard transformation trns 
C>     tol is set to 1E-12, which should work with 
C>     REAL = double precision 
C>     Version vom 21.08.2008 (AN)

      MODULE gauss
       USE precision
       PRIVATE

       PUBLIC gauleg,trns,impulse,equidist,chebyshev,trnsmin,expdistr,simpson

       CONTAINS 
       
!    this routine only print out info on the version number 
       
       
       
       ! subroutine to generate the grid points and width for Simpson integration
       subroutine simpson(nnu,xival,rr,rw)
       implicit none 
       integer :: nxmx,nnu
       real*8 :: xival
       real*8,dimension(1:nnu) :: rr,rw
       real*8,dimension(1:nnu+1) :: sxx,wxx
       real*8 :: dx, d43,d23
       integer :: nxmx1, nxmx2,i 
       
        nxmx=nnu+1
        dx=xival/float(nxmx-1)
 
        d43=4.d0/3.d0
        d23=2.d0/3.d0
   
        wxx(1)=dx/3.d0
        wxx(nxmx)=dx/3.d0
        sxx(1)=0.d0
        sxx(nxmx)=float(nxmx-1)*dx

 
        nxmx1=nxmx-1
        nxmx2=nxmx-2
 
        do 50 i=2,nxmx1,2
        sxx(i)=float(i-1)*dx
  50    wxx(i)=d43*dx
 
        do 55 i=3,nxmx2,2
        sxx(i)=float(i-1)*dx
  55    wxx(i)=d23*dx
  
  
        do i=1,nnu
          rr(i)=sxx(i+1)
          rw(i)=wxx(i+1)
        end do 
       end subroutine  

        
C>     this subroutine calculates standard gauss Legendre points 
C>     between x1 and x2 (usually -1.0_dpreal and  1.0_dpreal) 
C>     N is the number of mesh points required. 
C>     The grid and the weights are stored in the arrays X and W 
C>     @param[in] x1 lower boundary 
C>     @param[in] x2 upper boundary 
C>     @param[in] N number of grid points 
C>     @param[out] X grid points 
C>     @param[out] W integration weights
       SUBROUTINE gauleg(N,x1,x2,X,W)
        IMPLICIT NONE
        INTEGER N        
        REAL(dpreal) x1,x2,X(N),W(N)
        REAL(dpreal) z1,z,xm,xl,pp,p3,p2,p1,pi,tol
        INTEGER m,i,j
        
        pi=acos(-1.0)
        tol=1.E-12
        
        m=(n+1)/2
        xm=0.5*(x2+x1)
        xl=0.5*(x2-x1)
        
        DO 10 i=1,m
         z=cos(pi*(i-0.25)/(N+0.5))
         
 20      CONTINUE
         p1=1.0E0
         p2=0.0E0
         DO 30 j=1,N
          p3=p2
          p2=p1
          p1=((2*j-1)*z*p2-(j-1)*p3)/j
 30      CONTINUE
         pp=N*(z*p1-p2)/(z*z-1.0E0)
         z1=z
         z=z1-p1/pp
         IF( abs(z1-z) .GT. tol) GOTO 20 ! Scheifenende
         
         X(i) = xm - xl*z
         X(n+1-i) = xm + xl*z
         W(i) = 2.E0*xl/((1.0-z*z)*pp*pp)
         W(n+1-i) = W(i)
 10     CONTINUE
       END SUBROUTINE gauleg


C     Impulsroutine von Dirk Hueber
C     for compatibility with V3NF codes 
C     no integration weights will be returned 

      SUBROUTINE IMPULSE(NN,PP,P,W,PBAR)
C
C     Diese Subroutine berechnet die Impulspunkte 
C
      PARAMETER (NMAX=30)
      INTEGER NN(4)
      REAL(dpreal) PP(5),P(NMAX),W(NMAX)
      REAL(dpreal) PBAR

      PBAR=PP(5)
      NUM=0
      DO 10 I=1,4
       IF (NN(I).EQ.0) GOTO 10
       AL=(PP(I+1)-PP(I))/NN(I)
       DO 20 J=1,NN(I)
        IF (I.EQ.2) THEN
         P(J+NUM)=PP(I)+(PP(I+1)-PP(I))*(REAL(J,dpreal)
     $         /REAL(NN(I),dpreal))**2
        ELSE
         P(J+NUM)=PP(I)+J*AL
        ENDIF
        W(J+NUM)=0.
20     CONTINUE
       NUM=NUM+NN(I)
10    CONTINUE

      END SUBROUTINE

C     the subsroutine trns uses gauleg to obtain Gauss-Legendre grid points 
C     and performs a hyperbolic transformation on the points and weights 
C     resulting new points and weights that are more even spread in the intervals
C     The grid has NP1/2 points between 0 and P1, NP1/2 points between P1 and P2 
C     and NP2 points between P2 and P3 
C     The complete grid runs from 0 to P3 and has NP=NP1+NP2 grid points
C     It is your responsibility to make sure that NP=NP1+NP2
C     P1,P2,P3 defines the intervals 
C     NP1,NP2,NP=NP1+NP2 the number of mesh points. 
C     grid and weights are stored in XP and AP on exit 


       SUBROUTINE TRNS(NP1,NP2,NP,P1,P2,P3,XP,AP)
       IMPLICIT NONE  
C     ===============
C
C     TRNS BELEGT DIE FELDER XP UND AP MIT TRANSFORMIERTEN
C     GAUSS-LEGENDRE-PUNKTEN UND GEWICHTEN
C
C     NP1 PUNKTE WERDEN UEBER DIE HYPERBOLISCHE TRANSFORMATION
C
C     X --> (1.+X) / (1./P1-(1./P1-2./P2)*X)
C
C     AUF DAS INTERVALL (0.;P2) ABGEBILDET, WOBEI
C     NP1/2 PUNKTE IN (0.;P1) UND
C     NP1/2 PUNKTE IN (P1;P2) LIEGEN
C
C     NP2 PUNKTE WERDEN UEBER DIE LINEARE TRANSFORMATION
C
C     X --> (P3+P2)/2. + (P3-P2)/2.*X
C
C     AUF DAS INTERVALL (P2;P3) ABGEBILDET
C
C     NP = NP1 + NP2
C
       REAL(dpreal) XP1(NP1),AP1(NP1),XP2(NP2),AP2(NP2)
       REAL(dpreal) XP(NP),AP(NP)
       REAL(dpreal) P1,P2,P3
       REAL(dpreal),PARAMETER :: eins=1.0_dpreal
       INTEGER NP1,NP2,NP
       REAL(dpreal) XX,X,A,DELPH
       INTEGER I

       CALL gauleg(NP1,-eins,eins,XP1,AP1)
       
       DO 1 I=1,NP1
        X=XP1(I)
        A=AP1(I)
        XX=1.0_dpreal/P1-(1.0_dpreal/P1-2.0_dpreal/P2)*X
        XP1(I)=(1.0_dpreal+X) / XX
 1      AP1(I)=(2.0_dpreal/P1-2.0_dpreal/P2)*A / XX**2
C
       IF(NP2 .NE. 0) THEN

        CALL gauleg(NP2,-eins,eins,XP2,AP2)


        DO 2 I=1,NP2
         X=XP2(I)
         A=AP2(I)
         DELPH=(P3-P2)/2.0_dpreal
         XP2(I)=(P3+P2)/2.0_dpreal + DELPH*X
 2       AP2(I)=DELPH*A
       ENDIF
C
       DO 3 I=1,NP1
        XP(I)=XP1(I)
 3      AP(I)=AP1(I)
C
       IF(NP2 .NE. 0) THEN
        DO 4 I=1,NP2
         XP(I+NP1)=XP2(I)
 4       AP(I+NP1)=AP2(I)
       ENDIF
C
       RETURN
      END SUBROUTINE trns
      
C  same as TRNS except that the first interval starts not at zero but at PMIN       

      SUBROUTINE TRNSMIN(NP1,NP2,NP,PMIN,P1,P2,P3,XP,AP)
       IMPLICIT NONE
       REAL(dpreal) XP1(NP1),AP1(NP1),XP2(NP2),AP2(NP2)
       REAL(dpreal) XP(NP),AP(NP)
       REAL(dpreal) P1,P2,P3,PMIN
       REAL(dpreal),PARAMETER :: eins=1.0_dpreal
       INTEGER NP1,NP2,NP,I
       REAL(dpreal) AA,BB,CC,X,A,XX1,XX2,DELPH
       
C     =============
C
C     TRNS BELEGT DIE FELDER XP UND AP MIT TRANSFORMIERTEN
C     GAUSS-LEGENDRE-PUNKTEN UND GEWICHTEN
C
C     NP1 PUNKTE WERDEN UEBER DIE HYPERBOLISCHE TRANSFORMATION
C
C     X --> (A+B*X) / (1.+C*X)
C     AA=P1 ; BB=( P1*(PMIN+P2)-2.*PMIN*P2 )/( P2-PMIN ) ;
C     CC=( 2.*P1-PMIN-P2 )/( P2-PMIN )
C
C     AUF DAS INTERVALL (PMIN;P2) ABGEBILDET, WOBEI
C     NP1/2 PUNKTE IN (PMIN;P1) UND
C     NP1/2 PUNKTE IN (P1;P2) LIEGEN
C
C     NP2 PUNKTE WERDEN UEBER DIE LINEARE TRANSFORMATION
C
C     X --> (P3+P2)/2. + (P3-P2)/2.*X
C
C     AUF DAS INTERVALL (P2;P3) ABGEBILDET
C
C     NP = NP1 + NP2
C
C
      AA=P1
      BB=( P1*(PMIN+P2)-2.*PMIN*P2 )/( P2-PMIN )
      CC=( 2.*P1-PMIN-P2 )/( P2-PMIN )
      CALL gauleg(NP1,-eins,eins,XP1,AP1)
      DO 1 I=1,NP1
      X=XP1(I)
      A=AP1(I)
      XX1=AA+BB*X
      XX2=1.+CC*X
      XP1(I)=XX1/XX2
    1 AP1(I)=( BB*XX2-XX1*CC )*A/(XX2*XX2)
      IF(NP2 .NE. 0) THEN
      CALL gauleg(NP2,-eins,eins,XP2,AP2)
      DO 2 I=1,NP2
      X=XP2(I)
      A=AP2(I)
      DELPH=(P3-P2)/2.
      XP2(I)=(P3+P2)/2. + DELPH*X
    2 AP2(I)=DELPH*A
      ENDIF
C
      DO 3 I=1,NP1
      XP(I)=XP1(I)
    3 AP(I)=AP1(I)
C                QQ
      IF(NP2 .NE. 0) THEN
      DO 4 I=1,NP2
      XP(I+NP1)=XP2(I)
    4 AP(I+NP1)=AP2(I)
      ENDIF
C
      RETURN
      END SUBROUTINE TRNSMIN

!! open formular of Numerical Recipes in FORTRAN 2nd ed. 
!!  based on equation 4.1.18
 
      SUBROUTINE equidist(NP1,NP2,NP,P1,P2,P3,XP,AP)
       IMPLICIT NONE

       INTEGER NP1,NP2,NP 
       REAL(dpreal) XP(NP),AP(NP),P1,P2,P3

       INTEGER ip
       REAL(dpreal) delta1,delta1_24th
       REAL(dpreal) delta2,delta2_24th


       IF(NP1+NP2.NE.NP) STOP 'NP not consistent with NP1,NP2'

       IF(NP1.GE.6) THEN  
         
         delta1=(P2-P1)/(NP1+1)
         delta1_24th=delta1/24.0
         

         XP(1)=P1+delta1
         XP(2)=P1+2.0*delta1
         XP(3)=P1+3.0*delta1
         AP(1)=delta1_24th*55.0
         AP(2)=-delta1_24th*4.0
         AP(3)=delta1_24th*33.0


         XP(NP1-2)=P1+(NP1-2)*delta1
         XP(NP1-1)=P1+(NP1-1)*delta1
         XP(NP1)=P1+NP1*delta1
         AP(NP1-2)=delta1_24th*33.0
         AP(NP1-1)=-delta1_24th*4.0
         AP(NP1)=delta1_24th*55.0


         DO ip=4,NP1-3
          XP(ip)=P1+ip*delta1
          AP(ip)=delta1
         END DO

       ELSE
         IF(NP1.EQ.0) THEN 
           CONTINUE
         ELSE IF(NP1.EQ.1) THEN
           AP(1)=P2-P1
           XP(1)=P1
         ELSE
           STOP 'NP1 LE 6 not implemented'
         END IF
       END IF
      
       IF(NP2.GE.6) THEN

         delta2=(P3-P2)/(NP2+1)
         delta2_24th=delta2/24.0
         

         XP(NP1+1)=P2+delta2
         XP(NP1+2)=P2+2.0*delta2
         XP(NP1+3)=P2+3.0*delta2
         AP(NP1+1)=delta2_24th*55.0
         AP(NP1+2)=-delta2_24th*4.0
         AP(NP1+3)=delta2_24th*33.0


         XP(NP1+NP2-2)=P2+(NP2-2)*delta2
         XP(NP1+NP2-1)=P2+(NP2-1)*delta2
         XP(NP1+NP2)=P2+NP2*delta2
         AP(NP1+NP2-2)=delta2_24th*33.0
         AP(NP1+NP2-1)=-delta2_24th*4.0
         AP(NP1+NP2)=delta2_24th*55.0


         DO ip=4,NP2-3
          XP(NP1+ip)=P2+ip*delta2
          AP(NP1+ip)=delta2
         END DO

       ELSE
         IF(NP2.EQ.0) THEN 
           CONTINUE
         ELSE IF(NP2.EQ.1) THEN
           AP(NP1+1)=P3-P2
           XP(NP1+1)=P2
         ELSE
           STOP 'NP2 LE 6 not implemented'
         END IF
       END IF


      END SUBROUTINE

      SUBROUTINE chebyshev(NP,xa,xb,XP,AP,gp,gpp)
       IMPLICIT NONE
       INTEGER NP
       REAL(dpreal) xa,xb,XP(NP),AP(NP),gp(NP,NP),gpp(NP,NP)
       INTEGER ip,ik,ij
       REAL(dpreal) ccoef(NP+1,NP),pikonst
       REAL(dpreal) sum,gfunc,gprime,x

       INTEGER k,i


       gfunc(x)=(1.0+x)/(1.0/xa-(1.0/xa-2.0/xb)*x)

       gprime(x)=((1.0/xa-(1.0/xa-2.0/xb)*x)
     X                +(1.0/xa-2.0/xb)*(1.0+x))
     X                 /(1.0/xa-(1.0/xa-2.0/xb)*x)**2


       AP(1:NP)=0.0
       pikonst=ACOS(-1.0)

       DO k=1,NP
        XP(NP+1-k)=gfunc(cos(pikonst*(k-0.5)/NP))

        DO i=2,NP-1,2
         AP(NP+1-k)=AP(NP+1-k)
     X        +2.0/(NP*(i-1))
     X         *(cos(pikonst*(i-2)*(k-0.5)/NP)
     X                 -cos(pikonst*i*(k-0.5)/NP))
     X            *gprime(cos(pikonst*(k-0.5)/NP))
        END DO
       END DO

       sum=0.0
       DO k=1,NP
        sum=sum+AP(k)
       END DO

       ccoef(1:NP+1,1:NP)=0.0

       DO ik=1,NP

        DO ij=NP,2,-1
         ccoef(ij-1,ik)=ccoef(ij+1,ik)
     X      +4.0*(ij-1.0)/NP*cos(pikonst*(ij-1.0)*(ik-0.5)/NP)
        END DO
       END DO

       
       DO ik=1,NP
        DO ip=1,NP
         gp(NP+1-ip,NP+1-ik)=-0.5*ccoef(1,ik)
         DO ij=1,NP-1
          gp(NP+1-ip,NP+1-ik)=gp(NP+1-ip,NP+1-ik)
     X         +ccoef(ij,ik)*cos(pikonst*(ij-1.0)*(ip-0.5)/NP)
         END DO
        END DO
       END DO
       

       DO ik=1,NP
        DO ip=1,NP
         gp(NP+1-ip,NP+1-ik)=gp(NP+1-ip,NP+1-ik)
     X      /gprime(cos(pikonst*(ip-0.5)/NP))
        END DO
       END DO

       DO ik=1,NP
        DO ip=1,NP
         gpp(ip,ik)=0.0
         DO ij=1,NP
          gpp(ip,ik)=gpp(ip,ik)+gp(ip,ij)*gp(ij,ik)
         END DO         
        END DO
       END DO       

       
      END SUBROUTINE
      
      SUBROUTINE expdistr(NP,P1,P2,XP,AP)
       IMPLICIT NONE
       INTEGER NP
       REAL(dpreal) P1,P2,XP(NP),AP(NP)
       REAL(dpreal) a,alpha,GP(NP),GW(NP)
       INTEGER ip

       CALL gauleg(NP,-1.0_dpreal,1.0_dpreal,GP,GW)
       
       alpha=0.5*log(P2/P1)
       a=P1

       DO ip=1,NP
        XP(ip)=a*exp(alpha*(GP(ip)+1.0))
        AP(ip)=alpha*XP(ip)*GW(ip)
       END DO

      END SUBROUTINE

      END MODULE gauss
