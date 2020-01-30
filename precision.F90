! this module provides parameters that choose the correct kind of REAL,COMPLEX and INTEGER variables 


module precision
  integer, parameter   :: sint=selected_int_kind(9)    ! integers up to 10^9 usually 4 byte (MPI_INTEGER4)
  integer, parameter   :: lint=selected_int_kind(19)   ! integers up to 10^19 usually 8 byte (MPI_INTEGER8) 

  integer, parameter   :: spreal=selected_real_kind(6)   ! REAL/COMPLEX with 7 digits usually 4 byte (MPI_REAL4,MPI_COMPLEX8)
  integer, parameter   :: dpreal=selected_real_kind(12)  ! REAL/COMPLEX with 12 digits usually 8 byte (MPI_REAL8,MPI_COMPLEX16) 

  integer, parameter   :: maxinteger=2147483647        ! maximum single precision integer 2**(31)-1
  
  integer, parameter   :: nprint=0   ! controls the details printed out, production run = 0
  logical,public       :: print_infor = .true.
  real(dpreal),public :: eps_sp,eps_dp
  real(dpreal),parameter :: pi = acos(-1.0_dpreal)           ! pi and imaginary unit 
  complex(dpreal),parameter :: iu = (0.0_dpreal,1.0_dpreal)

   public init_eps,print_eps, print_data

 ! this subroutine initializes the machine precision for sp and dp

  contains

   subroutine init_eps
    implicit none

    real(spreal),external :: slamch
    real(dpreal),external :: dlamch
 
    eps_sp=slamch('E')
    eps_dp=dlamch('E')
    
   
   end subroutine init_eps

   subroutine print_eps
    implicit none
    
      ! first printout version running
    WRITE(*,*) 'Precision version: ',VERREV
    WRITE(*,*) 'Date             : ',VERDATE

    
    WRITE(*,*) 'Init precision: ',spreal,dpreal
    WRITE(*,*) 'eps:            ',eps_sp,eps_dp

    WRITE(*,*) 'COMPDATE: ',COMPDATE    
    WRITE(*,*) 'F90: ',MAKEF90
    WRITE(*,*) 'F77: ',MAKEF77
    WRITE(*,*) 'MAKEFOPT1:  ',MAKEFOPT1
    WRITE(*,*) 'MAKEFOPT2:  ',MAKEFOPT2
    WRITE(*,*) 'MAKEFOPT3:  ',MAKEFOPT3
    WRITE(*,*) 'MAKEFOPT4:  ',MAKEFOPT4
    WRITE(*,*) 'MAKEFOPT5:  ',MAKEFOPT5
    WRITE(*,*) 'MAKEFOPT6:  ',MAKEFOPT6
    WRITE(*,*) 'MAKEFOPT7:  ',MAKEFOPT7
    WRITE(*,*) 'MAKEFOPT8:  ',MAKEFOPT8
    WRITE(*,*) 'CC:  ',MAKECC
    WRITE(*,*) 'CC:  ',MAKECOPT
    WRITE(*,*) 'LIBSTD1:  ',MAKELIBSTD1
    WRITE(*,*) 'LIBSTD2:  ',MAKELIBSTD2
    WRITE(*,*) 'LIBSTD3:  ',MAKELIBSTD3
    WRITE(*,*) 'LIBSTD4:  ',MAKELIBSTD4
    WRITE(*,*) 'LIBSTD5:  ',MAKELIBSTD5
    WRITE(*,*) 'LIBSTD6:  ',MAKELIBSTD6
    WRITE(*,*) 'LIBSTD7:  ',MAKELIBSTD7
    WRITE(*,*) 'LIBSTD8:  ',MAKELIBSTD8
    
   end subroutine print_eps

   !myid rank of  each processes
   subroutine print_data(myid)
    implicit none
    integer:: myid

    print_infor = (myid .eq. 0)

   end subroutine

end module precision
