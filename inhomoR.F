      program inhomoR
      use input
      use rmat
      use numerov
      use mesh
      use gauss
      use precision
      use green
      implicit none
      real*8 ::  t1, t2 , tr, tg, tn

      call readin()

      ! use Lagrange mesh method to solve the inhomogeneous equation
      allocate(rr(1:nr),rrw(1:nr))
      call gauleg(nr,0.0_dpreal,rmax,rr,rrw)
      
      call cpu_time(t1)
      call Inhomo_Rmatrix_method()
      call cpu_time(t2)
      tr=t2-t1
      
      
      call cpu_time(t1)
      call Inhomo_greenf_method()
      call cpu_time(t2)
      tg=t2-t1
      deallocate(rr,rrw)

      ! use Numerov method to solve the inhomogeneous equation
      nr=irmax
      allocate(rr(1:nr), rrw(1:nr))
      call simpson(nr,rmax,rr,rrw)
      call cpu_time(t1)
      call Inhomo_numerov_method()
      call cpu_time(t2)
      tn=t2-t1
      


      call fkind()
      
      write(*,*) "R-matrix takes time: ", tr
      write(*,*) "Green's function takes time: ", tg
      write(*,*) "Numerov takes time: ", tn
       
      end program
