      module mesh
      integer :: nr
      real*8 ::  rmax
c      wave functions calculated at intervals of HCM up to abs(RMATCH).
      real*8 :: rmatch
      real*8 :: hcm
      integer :: irmax
      real*8,allocatable,dimension(:) :: rr, rrw
      end module
