      module input
      implicit none
      real*8 :: Ecm, mu ! energy in CM frame and reduce mass : unit in MeV
      real*8 ::  k !  momentum
C     integer :: method ! method=1 lagrange mesh method, /=1 numerov method
      real*8 :: massa, massb !mass munber of interaction partciles
      real*8 :: za, zb ! charge of interaction partciles
      character*999 :: rho_file
      complex*16,allocatable,dimension(:) :: rho
      logical :: readrho
      logical,dimension(9999) :: written

      contains

      subroutine readin()
      use mesh
      use precision
      use constants
      use channels
      use pot
      use gauss
      implicit none



      namelist /global/ Ecm,hcm,rmax,lmin,lmax,nr
      namelist /systems/ massa, massb, za,zb,ja,jb
      namelist /source/ rho_file, readrho
      namelist /potential/ ptype,a1,a2,rc,uv,av,rv,uw,aw,rw,
     +                     vsov,rsov,asov,vsow,rsow,asow,vd,avd,rvd,wd,awd,rwd


      hcm=0.05;rmax=40;lmin=0;lmax=0;nr=60
      read(5,nml=global)
      irmax=nint(rmax/hcm)



      massa=1.0; massb=1.0; za=0.0;zb=0.0;ja=0.0;jb=0.0
      read(5,nml=systems)
      mu=massa*massb*amu/(massa+massb)
      write(*,*) "mu=",mu
      k=sqrt(2.0_dpreal*mu*ecm)/hbarc
      call alpha_2b()


      readrho=.False.
      read(5,nml=source)



       ptype=1;a1=0.0d0;a2=massb;rc=1.3d0
       uv=0.0d0;av=0.1d0;rv=0.0d0
       uw=0.0d0;aw=0.1d0;rw=0.0d0
       vsov=0.0d0;rsov=0.0d0;asov=0.1d0
       vsow=0.0d0;rsow=0.0d0;asow=0.1d0
       vd=0.0d0;avd=0.1d0;rvd=0.0d0
       wd=0.0d0;awd=0.1d0;rwd=0.0d0
       read(5,nml=potential)




      end subroutine

      subroutine read_rho()
      use interpolation
      use precision
      use mesh
      implicit none
      logical uu
      real*8 :: x
      complex*16 :: y
      real*8, allocatable,dimension(:) :: xv(:)
      complex*16,allocatable,dimension(:) ::  faux(:)
      integer :: n_rho, n
      real*8,parameter:: alpha=0
      real*8 :: r

      uu=.false.
      inquire(file=trim(rho_file), exist=uu)
      if (.not.uu) then
        write(0,*) 'source term file:',trim(rho_file),' not found!'
        stop
      endif
      write(*,'(8x, "source term file:",a20)') trim(rho_file)


      open(unit=99, file=trim(rho_file),status="old") ! dependent on the format of rho  decide later
      rewind(99)
      n_rho=0
20    read(99,*,end=200) x,y
      n_rho=n_rho+1
      goto 20

200   if(n_rho < 1) write(*,*) "error reading source term ! "

      ! readin the desity file
      rewind(99)
      allocate(xv(n_rho), faux(n_rho))
      do n=1, n_rho
        read(99,*)x, y
        xv(n)=x
        faux(n)=y
      end do

      do n=1, nr
        r=rr(n)
        rho(n) = cfival(r, xv, faux, n_rho, alpha)
        write(31,*) r, real(rho(n)), aimag(rho(n))
      end do


      close(99)

      end subroutine



      subroutine fkind()
       character*40 flkind(9999)
       integer writf(9999),nwrit,i
       flkind(:) = ' '

       flkind(8)='channels coupling index'
       written(8)=.TRUE.


       flkind(30)='potential used in the calcuation'
       written(30)=.TRUE.

       flkind(31)='local copy of source term'
       written(31)=.TRUE.

       flkind(32)='solutions of r-matrix'
       written(32)=.TRUE.
       flkind(33)='solutions of Green function'
       written(33)=.TRUE.
       flkind(34)='solutions of numerov'
       written(34)=.TRUE.

       flkind(42)='module solutions of r-matrix'
       written(42)=.TRUE.
       flkind(43)='module solutions of Green function'
       written(43)=.TRUE.
       flkind(44)='module solutions of numerov'
       written(44)=.TRUE.

       flkind(200)='regular solution of green function'
       written(200)=.TRUE.
       flkind(250)='irregular solution of green function'
       written(250)=.TRUE.


       nwrit = 0
       do i=1,9999
        if(written(i)) then
        flkind(i) = trim(flkind(i))//'.'
        nwrit = nwrit+1
        writf(nwrit) = i
        endif
       enddo
       write(*,990) (writf(i),flkind(writf(i)),i=1,nwrit)
990    format(/'  The following files have been created:',
     X  /(2x,2(i3,':',1x,a40)))
       return
       end subroutine



      end module
