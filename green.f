      module green
      type greenf
         complex*16,allocatable,dimension(:,:) :: re,ir
      end type
      type(greenf) :: Gx
      real*8,private,dimension(:),allocatable :: nfc,ngc,nfcp,ngcp
      contains


      subroutine Inhomo_greenf_method()
      use mesh
      use channels
      use input
      use precision
      use pot
      use constants
      implicit none
      integer :: nch
      complex*16,dimension(1:nr) :: wf
      integer :: l
      real*8 :: s,j,ls
      integer :: ir, irp
      real*8 :: N



      call greenfunc()

      N=2.0_dpreal*mu/hbarc/hbarc/k
      do nch=1, alpha2b%nchmax
        write(*,*)"nch=",nch
        l=alpha2b%l(nch)
        s=alpha2b%s(nch)
        j=alpha2b%j(nch)
        ls=0.5_dpreal*(j*(j+1)-l*(l+1)-s*(s+1))
        call potr(za*zb,ls)
        if (.not. readrho) rho(:)=v(:)

        wf=0.0_dpreal
        do ir=1, nr
          do irp=1,nr
             wf(ir)=wf(ir) + N*Gx%re(min(ir,irp),nch)*Gx%ir(max(ir,irp),nch)*rho(irp)*rrw(irp)
          end do
        end do

        write(33,101)l,s,j
        write(43,101)l,s,j
        do ir=1,nr
           write(33,*) rr(ir), real(wf(ir)),aimag(wf(ir))
           write(43,*) rr(ir), abs(wf(ir))
        end do

      end do
101   format('&l=',I3,' s=',f3.1, ' j=', f3.1)



      end subroutine






c-----------------------------------------------------------------------
      subroutine greenfunc()
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       use channels
       use mesh
       use constants
       use input
       use pot
       use precision
       use coulfunc
       use interpolation
       use gauss
       implicit none
       integer :: l,nch
       real*8 :: s,j,ls
       integer :: ifail ! for subroutine coul90
       integer :: ir
       real*8 :: r, kr
       real*8 :: eta !Sommerfeld parameter
       complex*16 :: sl,nl !  s-matrix and normalization parameter
       integer :: r0
       complex*16  :: hlp,flp
       complex*16,dimension(0:irmax,alpha2b%nchmax) :: flx,hlx
       real*8,dimension(0:lmax) :: cph
       integer :: nrg

       if (allocated(nfc)) deallocate(nfc)
       if (allocated(ngc)) deallocate(ngc)
       if (allocated(nfcp)) deallocate(nfcp)
       if (allocated(ngcp)) deallocate(ngcp)
       if (allocated(Gx%re)) deallocate(Gx%re)
       if (allocated(Gx%ir)) deallocate(Gx%ir)


       allocate(Gx%re(1:nr,1:alpha2b%nchmax))
       allocate(Gx%ir(1:nr,1:alpha2b%nchmax))
       allocate(nfc(0:lmax),ngc(0:lmax),nfcp(0:lmax),ngcp(0:lmax))
       flx=0.0_dpreal
       hlx=0.0_dpreal



        

        kr=hcm*(irmax-2)*k
        eta=za*zb*e2*mu/hbarc/hbarc/k

        call coulph(eta,cph,lmax)
        call coul90(kr,eta,zero,lmax,nfc,ngc,nfcp,ngcp,0,ifail)
        if (ifail/=0) then
           write(*,*) 'coul90: ifail=',ifail; stop
        endif

        


       ! change grid points from gauleg to simpson
       nrg=nr
       deallocate(rr,rrw)
       nr=irmax
       allocate(rr(1:nr), rrw(1:nr))
       call simpson(nr,rmax,rr,rrw)
        write(*,*)"alpha2b%nchmax=",alpha2b%nchmax
       do nch=1,alpha2b%nchmax

         l=alpha2b%l(nch)
         s=alpha2b%s(nch)
         j=alpha2b%j(nch)
         ls=0.5_dpreal*(j*(j+1)-l*(l+1)-s*(s+1))
    
         write(*,*)"nch=",nch
         call potr(za*zb,ls)
         r0=2*l
       	 call sch(r0,v,l,flx(:,nch))
         call matching(l,k,flx(irmax-4:irmax,nch),sl,nl)
         

       	 flx(:,nch)=flx(:,nch)*nl

         call hlxkr(l,v,cph(l),hlx(:,nch))

     
           do ir=0,irmax
              write (200,*) hcm*ir, real(flx(ir,nch)), aimag(flx(ir,nch))
           end do
       end do

               
       ! change grid points from simpson  to gauleg
       deallocate(rr,rrw)
       nr=nrg
       allocate(rr(1:nr), rrw(1:nr))
       call gauleg(nr,0.0_dpreal,rmax,rr,rrw)

       do nch=1,alpha2b%nchmax
         do ir=1,nr
           r=rr(ir)
           Gx%re(ir,nch)=FFC(r/hcm,flx(:,nch),irmax+1)
           Gx%ir(ir,nch)=FFC(r/hcm,hlx(:,nch),irmax+1)
         end do
       end do





              

      end subroutine
c-----------------------------------------------------------------------





ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine sch(r0,vpot,l,rwfl)
c     mu! reduce mass
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       use mesh
       use constants, only:hbarc,e2
       use precision
       use input
       implicit none
       integer,intent(in) :: r0
       complex*16,dimension(1:irmax),intent(in) :: vpot
       integer,intent(in) :: l
       integer :: ir
       real*8 :: r
       complex*16,dimension(0:irmax),intent(out) :: rwfl   ! partial radial wave function
       complex*16,dimension(0:irmax) :: kl
       complex*16,dimension(0:irmax) :: Tx
       complex*16,dimension(0:irmax) :: Wx


        kl=0.0d0
        rwfl=0.0d0
        Tx=0.0d0


c      Numerov method to solve the differential radial equation
       rwfl(r0)=0      ! boundary condition


       ir=r0+1; r=ir*hcm
       rwfl(ir)=hcm**(l+1)  ! arbitrary value


       kl(ir)=2.*mu*ecm/hbarc**2-l*(l+1)/r**2-2.*mu*Vpot(ir)/hbarc**2
       Tx(ir)=-hcm**2/12.0d0*kl(ir)
       Wx(ir)=(1-Tx(ir))*rwfl(ir)

       rwfl(r0+2)=2.*rwfl(r0+1)-hcm**2*kl(r0+1)*rwfl(r0+1)
       ir=r0+2; r=ir*hcm
       kl(ir)=2.*mu*ecm/hbarc**2-l*(l+1.)/r**2-2.*mu*Vpot(ir)/hbarc**2
       Tx(ir)=-hcm**2/12.0d0*kl(ir)
       Wx(ir)=(1-Tx(ir))*rwfl(ir)


       do ir=r0+2 ,irmax-1

        kl(ir+1)=2.*mu*ecm/hbarc**2-l*(l+1.)/((ir+1.)*hcm)**2-2.*mu*Vpot(ir+1)/hbarc**2
        Tx(ir+1)=-hcm**2/12.0d0*kl(ir+1)
        Wx(ir+1)=(2+12.*Tx(ir)+12.*Tx(ir)**2)*Wx(ir)-Wx(ir-1)
        rwfl(ir+1)=Wx(ir+1)/(1.-Tx(ir+1))

       end do


      end subroutine sch


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine matching(l,k,wf,sl,nl) !wf has dimension (5)
c     nl*wf=0.5*i*(H(-)-sl*H(+))
c     nl*wfp=0.5*i*k*(H'(-)-sl*H'(+))
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       use precision, only:pi,iu
       use mesh,only:hcm
       implicit none

       integer,intent(in) :: l
       real*8,intent(in) :: k
       complex*16,intent(in),dimension(1:5) :: wf ! wavefunction at rmatch
       complex*16 :: wfp !derivative of wf
       complex*16,intent(out) :: nl !Normalization parameter
       complex*16 :: hc,hc1  !H(+),H(-)
       complex*16 ::hcp,hcp1 ! derivatives of H(+),H(-)
       complex*16,intent(out) :: sl ! S-matrix

       hc=cmplx(ngc(l),nfc(l),kind=8)
       hc1=cmplx(ngc(l),-nfc(l),kind=8)
       hcp=cmplx(ngcp(l),nfcp(l),kind=8)
       hcp1=cmplx(ngcp(l),-nfcp(l),kind=8)

       wfp=(-wf(5)+8.*wf(4)-8.*wf(2)+wf(1))/12./hcm
       nl=(hc1*hcp*iu*k-hc*hcp1*iu*k)/(2.*(hcp*wf(3)*k-hc*wfp))
       sl=(hc1*wfp-hcp1*wf(3)*k)/(hc*wfp-hcp*wf(3)*k)


      end subroutine matching
c----------------------------------------------------------------------
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine hlxkr(l,vpot,cphase,hlx)
c     This subroutine is used to calculate the hlx(r)
c              for each energy point
c           hlx(l,irmax)=ngc(l)+i*nfc(l)
c           use Numerov method from outside to inside
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       use mesh,only:irmax,hcm
       use channels,only:lmax
       use constants, only:hbarc,e2,zero
       use precision
       use coulfunc
       use input
       implicit none
       integer :: l
       complex*16,dimension(0:irmax) :: hlx
       complex*16,dimension(1:irmax) :: Vpot
       real*8 :: cphase
       real*8 :: r,eta,kr
       real*8 :: z12 ! For coulomb potential and eta zt*zx
       real*8 :: const
       integer :: ir ,ifail
       complex*16,dimension(1:irmax) :: kl
       real*8,dimension(0:lmax) :: nfc1,ngc1,nfcp1,ngcp1

       z12=za*zb
       eta=z12*e2*mu/hbarc/hbarc/k
       hlx=0.0d0

c----boundary condition
       kr=hcm*(irmax)*k
       call coul90(kr,eta,zero,lmax,nfc1,ngc1,nfcp1,ngcp1,0,ifail)
       hlx(irmax)=cmplx(ngc1(l),nfc1(l),kind=8)

       kr=hcm*(irmax-1)*k
       call coul90(kr,eta,zero,lmax,nfc1,ngc1,nfcp1,ngcp1,0,ifail)
       hlx(irmax-1)=cmplx(ngc1(l),nfc1(l),kind=8)



c      Numerov method to solve the differential radial equation

        ir=irmax;r=ir*hcm
       	kl(ir)=2.*mu*ecm/hbarc**2-l*(l+1.)/r**2-2.*mu*Vpot(ir)/hbarc**2

       	ir=irmax-1;r=ir*hcm
       	kl(ir)=2.*mu*ecm/hbarc**2-l*(l+1.)/r**2-2.*mu*Vpot(ir)/hbarc**2

       	const=hcm**2/12.


       	do ir=irmax-1,2,-1
           kl(ir-1)=2.*mu*ecm/hbarc**2-l*(l+1.)/((ir-1.)*hcm)**2-2.*mu*Vpot(ir-1)/hbarc**2
           hlx(ir-1)=((2.-10.*const*kl(ir))*hlx(ir)-(1.+const*kl(ir+1))*hlx(ir+1))/(1.+const*kl(ir-1))
       	end do


        write(250,*)'& l=',l
       	do ir=1,irmax
       	  write (250,*) hcm*ir, real(hlx(ir)),aimag(hlx(ir))
       	end do


      end subroutine hlxkr
c-----------------------------------------------------------------------




      end module
