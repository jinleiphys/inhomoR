      program inhomoR
      use input
      use rmat
      use numerov
      use mesh
      use gauss
      use precision
      use green
      implicit none

      call readin()

      ! use Lagrange mesh method to solve the inhomogeneous equation
      allocate(rr(1:nr),rrw(1:nr))
      call gauleg(nr,0.0_dpreal,rmax,rr,rrw)
      call Inhomo_Rmatrix_method()
      call Inhomo_greenf_method()
      deallocate(rr,rrw)

      ! use Numerov method to solve the inhomogeneous equation
      nr=irmax
      allocate(rr(1:nr), rrw(1:nr))
      call simpson(nr,rmax,rr,rrw)
      call Inhomo_numerov_method()


      call fkind()


      end program
