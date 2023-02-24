module chorus_IO
  implicit none
  integer, parameter :: singlep = kind(1.0), doublep = kind(1.0d0)
  integer, parameter :: fp = doublep
  real(fp), parameter :: pi =  3.1415926535897932384626433832795028841971693993751_fp
  complex(fp), parameter :: mathj = (0.0_fp, 1.0_fp)
  real(fp), parameter :: small = 1.0e-10_fp, TINY_NUM = 1.0e-16_fp
  
  !!code units : [T]:1/omegap, [L]:lightc/omegap, [A]:m_e*lightc**2/e, 
  !!             [J,pcor,pcr]:m_e*lightc**2/omegap, [Density]:n0
  !!physical units: CGS
  !! J=mu-pcor-pcr,  pcor=pll/kmode-pcr, pcr=(omega-gyro)/kmode**2

  ! equilibrium magnetic field
  real(fp), parameter :: lightc = 2.99792458e10_fp ! speed of lights  (cm/s)
  real(fp), parameter :: Lshell = 5.0_fp           ! L-shell
  real(fp), parameter :: gyro0 = 0.2_fp            ! gyro-frequency at the equator 

  ! cold plasma frequency at the equator (plasma frequency omega_ce0/gyro0)
  real(fp), parameter :: omegap = 1.76e7_fp*3.12e-1_fp/(Lshell**3*gyro0) 
  ! earth radius in dimensionless unit
  real(fp), parameter :: RE = 6.375e8_fp*omegap/lightc  
  ! field line dependent of cold plasma density distribution
  ! wp2 = 1.0+cold2*lambda**2
  real(fp), parameter :: cold2 = 1.0_fp      

  ! Bi-Maxwellian distributed energetic particles
  real(fp), parameter :: loss_cone = 1e-15  ! depth of the loss cone
  real(fp), parameter :: vll = 0.15_fp           ! parallel thermal velocity = sqrt(Tll/me)
  real(fp), parameter :: vperp = 0.30_fp         ! perpendicular thermal velocity = sqrt(Tperp/me)
  real(fp), parameter :: aniso = (vperp/vll)**2  ! anisotropy Tperp/Tll
  real(fp), parameter :: nh = 0.002_fp           ! hot particle density 
  real(fp), parameter :: colla = 0.0_fp,  colld3 = 0.0_fp   ! collison rate

  ! background noise level
  real(fp) :: noise = 1.0e-2_fp  
 
  ! simulation parameters
  integer, parameter  :: NT0 = 0, NT = 800
  integer, parameter  :: NZ = 500, NJ = 1
  !integer, parameter  :: NZ = 500, NJ = 9
  !* simulation latitude region [-15deg, 15deg]
  real(fp), parameter :: LZ = 15.0_fp/180.0_fp*pi*Lshell*RE
  real(fp), parameter :: Lq = 2.0_fp*pi
  real(fp), parameter :: Lp = 0.20_fp
  integer, parameter  :: Nq = 30
  integer, parameter  :: Np = 400
  real(fp), parameter :: dq = Lq/Nq
  real(fp), parameter :: dp = Lp/Np
 
  ! Runge-Kutta adapative parameters
  integer, parameter :: upwind = 3 ! 2 : 2nd order upwind scheme
                                   ! 3 : 3nd order upwind scheme
  !* smallest time step dtminw for wave; smallest time step dtminp for particle
  real(fp), parameter :: dtminw = 1.0_fp,  dtminp = 1.0e-6_fp
  !* tolerance rtolw for wave; tolerance rtolp for particle
  real(fp), parameter :: rtolw = 1.0e-6_fp, rtolp = 1.0e-2_fp
  real(fp), parameter :: safety = 0.9_fp
  real(fp), parameter :: pgrow = -1.0_fp/3.0_fp, pshrink = -0.5_fp,               &
                         errcon = (5.0_fp/safety)**(1.0_fp/pgrow)
  ! numerical diffusive coefficient
  real(fp), parameter :: diff_num = 1.0_fp

  ! code control
  logical, parameter :: restart = .FALSE.
  integer, parameter :: nsave = 100, ncheck = 1000
  integer, parameter :: stdout = 1           ! stdout=0: output to screen; 1: write to .log
  integer, parameter :: nthread = 1         ! edison:48,  cori: 64
  character(len=10)  :: logfile = "chorus.log"
  character(len=9)   :: real_fmt = "ES50.30E5" 
  character(len=24)  :: cmpl_fmt = "ES50.30E5, 1X, ES50.30E5"
  !zhengjs added

  integer, parameter :: use_adv = 1
  real(fp), parameter :: dZ = LZ/NZ
  real(fp) :: dT = 48.898864632961590

  interface save_output
     module procedure save_output_real
     module procedure save_output_cmpl
     module procedure save_output_real_1d
     module procedure save_output_real_2d
     module procedure save_output_cmpl_1d
     module procedure save_output_cmpl_2d
  end interface save_output
   
  interface write_log
     module procedure write_log_real
     module procedure write_log_int
  end interface write_log

contains
   subroutine write_log_real(message, data)
     implicit none
     character(len=*), intent(in) :: message
     real(fp), optional, intent(in) :: data 
     logical :: lopen
     if (stdout == 0) then
        if (present(data)) then
           write(*,*) message, data
        else
           write(*,*) message
        end if
     else
        inquire(file=trim(logfile), opened=lopen)
        if (.not. lopen) then
           if (NT0 == 0) then
              open(unit=20, file=trim(logfile), action='write', status='replace', position='append')
           else 
              open(unit=20, file=trim(logfile), action='write', position='append')
           end if
        end if
        if (present(data)) then
           write(20, *) message, data
        else
           write(20, *) message
        end if
     end if
     return
   end subroutine write_log_real

   
   subroutine write_log_int(message, data)
     implicit none
     character(len=*), intent(in) :: message
     integer, intent(in) :: data 
     logical :: lopen
     if (stdout == 0) then
        write(*,*) message, data
     else
        inquire(file=trim(logfile), opened=lopen)
        if (.not. lopen) then
           if (NT0 == 0) then
              open(unit=20, file=trim(logfile), action='write', status='replace', position='append')
           else
              open(unit=20, file=trim(logfile), action='write', position='append')
           end if
        end if
        write(20, *) message, data
     end if
     return
   end subroutine write_log_int

   
   function format_pattern(dtype, col)
     implicit none
     character(len=*), intent(in) :: dtype
     integer, intent(in) :: col
     character(len=256) :: format_pattern
     character(len=100) :: pattern
     character(len=10) :: num_col
     if (dtype=="real") then
        pattern = real_fmt
     else if (dtype=="complex") then
        pattern = cmpl_fmt
     else
        call write_log('Specify I/O datatype! (real or complex)')
        stop
     end if
     if (col == 1) then
        format_pattern = "("//trim(pattern)//")"
     else if (col > 1) then
        write(num_col, '(I10)') col
        format_pattern = "("//trim(adjustl(num_col))//"("//trim(pattern)//", ','))"
     else
        call write_log('Number of I/O data should be greater than 0! (col>0)')
        stop
     end if
     return 
   end function format_pattern

   
   subroutine save_output_real(file, data, is_append)
     implicit none
     character(len=*), intent(in) :: file
     real(fp), intent(in) :: data
     logical, intent(in) :: is_append
     character(len=256) :: output_fmt
     logical :: lopen
     output_fmt = format_pattern("real", 1)
     if (trim(file)==trim(logfile)) then
        inquire(file=trim(logfile), opened=lopen)
        if (.not. lopen) then
           open(unit=20, file=trim(logfile), action='write', position='append')
        end if
        write(20, *) data
        return
     end if
     if (is_append) then
        open(unit=30, file=trim(file), action='write', position='append')
        write(30, trim(output_fmt)) data
        close(unit=30)
     else
        open(unit=30, file=trim(file), action='write', status='unknown')
        write(30, trim(output_fmt)) data
        close(unit=30)
     end if
     return
   end subroutine save_output_real


   subroutine save_output_cmpl(file, data, is_append)
     implicit none
     character(len=*), intent(in) :: file
     complex(fp), intent(in) :: data
     logical, intent(in) :: is_append
     character(len=256) :: output_fmt
     logical :: lopen
     output_fmt = format_pattern("complex", 1)
     if (trim(file)==trim(logfile)) then
        inquire(file=trim(logfile), opened=lopen)
        if (.not. lopen) then
           open(unit=20, file=trim(logfile), action='write', position='append')
        end if
        write(20, *) data
        return
     end if
     if (is_append) then
        open(unit=30, file=trim(file), action='write', position='append')
        write(30, trim(output_fmt)) data
        close(unit=30)
     else
        open(unit=30, file=trim(file), action='write', status='unknown')
        write(30, trim(output_fmt)) data
        close(unit=30)
     end if
     return
   end subroutine save_output_cmpl


   subroutine save_output_real_1d(file, array, nrow, is_append)
     implicit none
     character(len=*), intent(in) :: file
     integer, intent(in) :: nrow
     real(fp), intent(in) :: array(nrow)
     logical, intent(in) :: is_append
     character(len=256) :: output_fmt
     logical :: lopen
     integer :: i
     output_fmt = format_pattern("real", 1)
     if (trim(file)==trim(logfile)) then
        inquire(file=trim(logfile), opened=lopen)
        if (.not. lopen) then
           open(unit=20, file=trim(logfile), action='write', position='append')
        end if
        do i =1, nrow
           write(20, *) array(i)
        end do
        return
     end if
     if (is_append) then
        open(unit=30, file=trim(file), action='write', position='append')
        do i = 1, nrow
           write(30, trim(output_fmt)) array(i)
        end do
        close(unit=30)
     else
        open(unit=30, file=trim(file), action='write', status='unknown')
        do i = 1, nrow
           write(30, trim(output_fmt)) array(i)
        end do
        close(unit=30)
     end if
     return
   end subroutine save_output_real_1d


   subroutine save_output_cmpl_1d(file, array, nrow, is_append)
     implicit none
     character(len=*), intent(in) :: file
     integer, intent(in) :: nrow
     complex(fp), intent(in) :: array(nrow)
     logical, intent(in) :: is_append
     character(len=256) :: output_fmt
     logical :: lopen
     integer :: i
     output_fmt = format_pattern("complex", 1)
     if (trim(file)==trim(logfile)) then
        inquire(file=trim(logfile), opened=lopen)
        if (.not. lopen) then
           open(unit=20, file=trim(logfile), action='write', position='append')
        end if
        do i =1, nrow
           write(20, *) array(i)
        end do
        return
     end if
     if (is_append) then
        open(unit=30, file=trim(file), action='write', position='append')
        do i = 1, nrow
           write(30, trim(output_fmt)) array(i)
        end do
        close(unit=30)
     else
        open(unit=30, file=trim(file), action='write', status='unknown')
        do i = 1, nrow
           write(30, trim(output_fmt)) array(i)
        end do
        close(unit=30)
     end if
     return
   end subroutine save_output_cmpl_1d
  

   subroutine save_output_real_2d(file, array, nrow, ncol, is_append)
     implicit none
     character(len=*), intent(in) :: file
     integer, intent(in) :: nrow, ncol
     real(fp), intent(in) :: array(nrow, ncol)
     logical, intent(in) :: is_append
     character(len=256) :: output_fmt
     logical :: lopen
     integer :: i, j
     output_fmt = format_pattern("real", ncol)
     if (trim(file)==trim(logfile)) then
        inquire(file=trim(logfile), opened=lopen)
        if (.not. lopen) then
           open(unit=20, file=trim(logfile), action='write', position='append')
        end if
        do i =1, nrow
           do j = 1, ncol
              write(20, *) array(i,j)
           end do
        end do
        return
     end if
     if (is_append) then
        open(unit=30, file=trim(file), action='write', position='append')
        do i = 1, nrow
           write(30, trim(output_fmt)) (array(i,j), j=1,ncol)
        end do
        close(unit=30)
     else
        open(unit=30, file=trim(file), action='write', status='unknown')
        do i = 1, nrow
           write(30, trim(output_fmt)) (array(i,j), j=1,ncol)
        end do
        close(unit=30)
     end if
     return   
   end subroutine save_output_real_2d


   subroutine save_output_cmpl_2d(file, array, nrow, ncol, is_append)
    implicit none
     character(len=*), intent(in) :: file
     integer, intent(in) :: nrow, ncol
     complex(fp), intent(in) :: array(nrow, ncol)
     logical, intent(in) :: is_append
     character(len=256) :: output_fmt
     logical :: lopen
     integer :: i, j
     output_fmt = format_pattern("complex", ncol)
     if (trim(file)==trim(logfile)) then
        inquire(file=trim(logfile), opened=lopen)
        if (.not. lopen) then
           open(unit=20, file=trim(logfile), action='write', position='append')
        end if
        do i =1, nrow
           do j = 1, ncol
              write(20, *) array(i,j)
           end do
        end do
        return
     end if
     if (is_append) then
        open(unit=30, file=trim(file), action='write', position='append')
        do i = 1, nrow
           write(30, trim(output_fmt)) (array(i,j), j=1,ncol)
        end do
        close(unit=30)
     else
        open(unit=30, file=trim(file), action='write', status='unknown')
        do i = 1, nrow
           write(30, trim(output_fmt)) (array(i,j), j=1,ncol)
        end do
        close(unit=30)
     end if
     return   
   end subroutine save_output_cmpl_2d

   
   subroutine create_directory(newdirpath)
     implicit none
     character(len=*), intent(in) :: newdirpath
     logical :: dir_exist
     inquire(file=trim(newdirpath)//'/.', exist=dir_exist)
     if (dir_exist)  then
        call execute_command_line('rm -r '//trim(newdirpath))
     end if
     call execute_command_line('mkdir -p '//trim(newdirpath))
     return
   end subroutine create_directory

end module chorus_IO


!!$program test_chorus_IO
!!$  use chorus_IO
!!$  implicit none
!!$  call  write_log("parameter upwind : out of range (upwind=2,3)")
!!$  stop
!!$end program test_chorus_IO
