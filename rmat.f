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
      implicit none
      real*8 :: eta, kr
      integer :: ifail
      integer :: nch, l
      real*8 :: s,ls,j


      !initial Lagrange mesh
      call rmat_ini(nr,rmax)
      allocate(fc(0:lmax),gc(0:lmax),dfc(0:lmax),dgc(0:lmax))

      ! compute the boundary conditions
      eta=za*zb*e2*mu/hbarc/hbarc/k
      kr=k*rmax
      call coul90(kr,eta,zero,lmax,fc,gc,dfc,dgc,0,ifail)


      do nch=1, alpha2b%nchmax


      end do 









      end subroutine




      end module
