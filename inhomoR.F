      program inhomoR
      use input
      use rmat
      use numerov
      implicit none

      call readin()

      if(method==1) then
      ! use Lagrange mesh method to solve the inhomogeneous equation
        call Inhomo_Rmatrix_method()
      else
      ! use Numerov method to solve the inhomogeneous equation
        call Inhomo_numerov_method()

      end if




      end program
