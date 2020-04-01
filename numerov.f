      module numerov
      implicit none
       real*8,private,dimension(:),allocatable :: nfc,ngc,nfcp,ngcp
       real*8,private,dimension(:),allocatable :: nfc1,ngc1,nfcp1,ngcp1
      contains


      subroutine Inhomo_numerov_method()
      use channels
      use mesh
      use coulfunc
      use input
      use constants
      use precision
      use pot
      implicit none
      real*8 :: eta, kr
      integer :: ifail
      integer :: nch, l
      real*8 :: s,ls,j
      complex*16,dimension(1:nr) :: wf
      complex*16 :: smat
      integer :: ir
      integer :: ie 



      allocate(nfc(0:lmax),ngc(0:lmax),nfcp(0:lmax),ngcp(0:lmax))
      allocate(nfc1(0:lmax),ngc1(0:lmax),nfcp1(0:lmax),ngcp1(0:lmax))
       
       
      ! add loop 
C     do ie=1, 1000
C     ecm=ecm+0.001
C     k=sqrt(2.0_dpreal*mu*ecm)/hbarc       
      ! compute the boundary conditions
      eta=za*zb*e2*mu/hbarc/hbarc/k
      kr=k*rmax
      call coul90(kr,eta,zero,lmax,nfc,ngc,nfcp,ngcp,0,ifail)   ! call for boundary of  inhomogeneous method
      kr=hcm*(nr-5)*k
      call coul90(kr,eta,zero,lmax,nfc1,ngc1,nfcp1,ngcp1,0,ifail)


      ! allocate source term
        deallocate(rho)
        allocate(rho(1:nr))

      if (readrho) call read_rho()





      do nch=1, alpha2b%nchmax
        l=alpha2b%l(nch)
        s=alpha2b%s(nch)
        j=alpha2b%j(nch)
        ls=0.5_dpreal*(j*(j+1)-l*(l+1)-s*(s+1))
        call potr(za*zb,ls)
        if (.not. readrho) then 
          do ir=1, nr 
            rho(ir)=v(ir)*sin(rr(ir))
          end do 
        end if 
        call numerov_inho(l,wf,smat)

        write(*,100) l,s,j,real(smat),aimag(smat)
        write(34,101)l,s,j
        write(44,101)l,s,j
        write(54,101)l,s,j

        do ir=1,nr
           write(34,*) rr(ir), real(wf(ir)),aimag(wf(ir))
           write(44,*) rr(ir), abs(wf(ir))
           write(54,*)  rr(ir), real(wf(ir))
        end do


      end do
      
      
C     end do  ! add loop 

100   format('l=',I3,' s=',f3.1, ' j=', f3.1, ' s-mat=(',2f10.6,')')
101   format('&l=',I3,' s=',f3.1, ' j=', f3.1)



      end subroutine


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine numerov_inho(l,wf,sl)
c      this subroutine is used to calculate the R_func in the notes
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        use mesh
        use constants
        use input
        use pot
        implicit none
        complex*16,dimension(1:nr) :: wf
        complex*16 :: smat
        ! for inhomogeneous equations
        complex*16,dimension(0:nr) :: inhomo,homo
        complex*16,dimension(1:nr) :: Tx
        complex*16,dimension(1:nr) :: Wx
        integer :: r0,ir,l
        complex*16,dimension(1:nr) :: kl
        real*8 :: const,r
        complex*16:: hc,hc1, sl,nl


        wf=0.0_dpreal



c      Numerov method to solve the differential radial equation
c      r0 the first point to do the integration
       inhomo=0.0_dpreal
       homo=0.0_dpreal
       r0=2*l
       rho=-rho*2.0_dpreal*mu/hbarc/hbarc                        !!!!!!!!
       inhomo(r0)=0.0d0    ! boundary condition
       inhomo(r0+1)=hcm**(l+1)  ! arbitrary value
       ir=r0+1; r=ir*hcm
       kl(ir)=2.*mu*ecm/hbarc**2-l*(l+1)/r**2-2.*mu*v(ir)/hbarc**2
       Tx(ir)=-hcm**2/12.0d0*kl(ir)
       Wx(ir)=(1-Tx(ir))*inhomo(ir)
       kl(ir+1)=2.*mu*ecm/hbarc**2-l*(l+1.)/((ir+1.)*hcm)**2-2.*mu*v(ir+1)/hbarc**2
       Tx(ir+1)=-hcm**2/12.0d0*kl(ir+1)
       Wx(ir+1)=(2+12.*Tx(ir)+12.*Tx(ir)**2)*Wx(ir)+hcm**2/11.*(rho(ir+1)+10*rho(ir))
       inhomo(ir+1)=Wx(ir+1)/(1.-Tx(ir+1))


       const=hcm**2/12.
       do ir=r0+2 ,nr-1
        kl(ir+1)=2.*mu*ecm/hbarc**2-l*(l+1.)/((ir+1.)*hcm)**2-
     &                  2.*mu*v(ir+1)/hbarc**2

        Tx(ir+1)=-hcm**2/12.0d0*kl(ir+1)
C       S=12.*Tx(ir)
        Wx(ir+1)=(2+12.*Tx(ir)+12.*Tx(ir)**2)*Wx(ir)-Wx(ir-1)
     &    +const*(rho(ir+1)+10*rho(ir)+rho(ir-1))
        inhomo(ir+1)=Wx(ir+1)/(1.-Tx(ir+1))
       end do




c************
c************homogeneous part
      homo(r0)=0.0d0    ! boundary condition
      homo(r0+1)=hcm**(l+1)  ! arbitrary value
c      homo(r0+1)=0.1
      ir=r0+1; r=ir*hcm
      kl(ir)=2.*mu*ecm/hbarc**2-l*(l+1)/r**2-2.*mu*v(ir)
     &                            /hbarc**2
      Tx(ir)=-hcm**2/12.0d0*kl(ir)
      Wx(ir)=(1-Tx(ir))*homo(ir)
      homo(r0+2)=2.*homo(r0+1)-hcm**2*kl(r0+1)*homo(r0+1)                             !!!!!!approximation

      ir=r0+2; r=ir*hcm

      kl(ir)=2.*mu*ecm/hbarc**2-l*(l+1)/r**2-2.*mu*v(ir)
     &                        /hbarc**2
      Tx(ir)=-hcm**2/12.0d0*kl(ir)
      Wx(ir)=(1-Tx(ir))*homo(ir)
      do ir=r0+2 ,nr-1
        kl(ir+1)=2.*mu*ecm/hbarc**2-l*(l+1.)/((ir+1.)*hcm)**2-
     &                  2.*mu*v(ir+1)/hbarc**2
        Tx(ir+1)=-hcm**2/12.0d0*kl(ir+1)
        Wx(ir+1)=(2+12.*Tx(ir)+12.*Tx(ir)**2)*Wx(ir)-Wx(ir-1)
        homo(ir+1)=Wx(ir+1)/(1.-Tx(ir+1))
       end do


c*************matching boundary condition**
      hc=cmplx(ngc(l),nfc(l),kind=8)
      hc1=cmplx(ngc1(l),nfc1(l),kind=8)
      sl=-(homo(nr)*inhomo(nr-5)-inhomo(nr)*
     &      homo(nr-5))/(homo(nr)*hc1-homo(nr-5)*hc)
      nl=(inhomo(nr)*hc1-inhomo(nr-5)*hc)/
     &                     (hc*(homo(nr-5))-hc1*homo(nr))



c--------
       wf=nl*homo(1:nr)+inhomo(1:nr)


       end subroutine
c-----------------------------------------------------------------------






      end module
