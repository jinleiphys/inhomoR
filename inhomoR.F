      program inhomoR
      use input
      use rmat
      implicit none

      call readin()

      if(method==1) then
      ! use Lagrange mesh method to solve the inhomogeneous equation
        call Inhomo_Rmatrix_method()
      else
      ! use Numerov method to solve the inhomogeneous equation

      end if




      end program
