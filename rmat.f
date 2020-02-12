      module rmat
      implicit none
      contains


      subroutine Inhomo_Rmatrix_method()
      use channels
      use lagrange_mesh_source
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
      real*8 :: hm
      complex*16,dimension(1:nr) :: wf
      complex*16 :: smat
      integer :: ir


      !initial Lagrange mesh
      call rmat_ini(nr,rmax)
      allocate(fc(0:lmax),gc(0:lmax),dfc(0:lmax),dgc(0:lmax))

      ! compute the boundary conditions
      eta=za*zb*e2*mu/hbarc/hbarc/k
      kr=k*rmax
      call coul90(kr,eta,zero,lmax,fc,gc,dfc,dgc,0,ifail)

      ! allocate source term
      if (.not. allocated(rho)) then
        allocate(rho(1:nr))
      end if
      if (readrho) call read_rho()

      hm=hbarc**2/(2.*mu)



      do nch=1, alpha2b%nchmax
        l=alpha2b%l(nch)
        s=alpha2b%s(nch)
        j=alpha2b%j(nch)
        ls=0.5_dpreal*(j*(j+1)-l*(l+1)-s*(s+1))
        call potr(za*zb,ls)
        if (.not. readrho) rho(:)=v(:)
        call rmat_inho(nr,rmax,v,rho,ecm,eta,hm,l,smat,wf)

        write(*,100) l,s,j,real(smat),aimag(smat)
        write(32,101)l,s,j
        write(42,101)l,s,j
        do ir=1,nr
           write(32,*) rr(ir), real(wf(ir)),aimag(wf(ir))
           write(42,*)rr(ir), abs(wf(ir))
        end do



      end do

100   format('l=',I3,' s=',f3.1, ' j=', f3.1, ' s-mat=(',2f10.6,')')
101   format('&l=',I3,' s=',f3.1, ' j=', f3.1)



      end subroutine




      end module
