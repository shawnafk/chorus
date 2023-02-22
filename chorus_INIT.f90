module chorus_INIT
  use chorus_IO
  use chorus_TYPE
  use chorus_LIB
  use deriv_class
  implicit none
  private

  public :: set_init, get_most_unstable

  real(fp) :: omega0    ! most unstable mode frequency
  real(fp) :: ratio

contains
  !subroutine set_init()   
  subroutine set_init(eon, chorus)   
    implicit none
    type(energetic_eon), intent(out) :: eon(-NZ:NZ,NJ)
    type(chorus_mode), intent(out) :: chorus(-NZ:NZ+1)
    real(fp) :: omega, gamma
    real(fp) :: zpos(-NZ:NZ), Jact(-NZ:NZ, NJ)
    real(fp) :: gyro(-NZ:NZ), wp2(-NZ:NZ), kmode(-NZ:NZ)                                
    real(fp) :: pcr(-NZ:NZ), dgyro(-NZ:NZ), dkmode(-NZ:NZ), Jconst
    integer :: i, j, k, l
    
    call get_most_unstable(omega, gamma)
    !zhengjs
    !call set_zpos(zpos)
    call set_zpos_uni(zpos)
    !call set_Jact(Jact, zpos, omega)
    call set_Jact(Jact, zpos, omega,Jconst)
    if ( RESTART )  then
       open(unit=50, file='./restart/eon.out', action='read',                            &
            status='old', form='unformatted')
       read(50) eon
       close(unit=50)
       open(unit=60, file='./restart/chorus.out', action='read',                         &
            status='old', form='unformatted')
       read(60)  chorus
       close(unit=60)   
    else
       gyro = gyro0*(1.0_fp+4.5_fp*(zpos/(Lshell*RE))**2)
       wp2 = 1.0_fp+cold2*(zpos/(Lshell*RE))**2
       kmode = sqrt(omega**2+wp2*omega/(gyro-omega))
       pcr = (omega-gyro)/kmode**2
       dgyro = 9.0_fp*gyro0*zpos/(Lshell*RE)**2 
       dkmode = (2.0_fp*cold2*omega*(gyro-omega)-9.0_fp*omega*gyro0*wp2)                 &
              / (2.0_fp*kmode*(gyro-omega)**2)*zpos/(Lshell*RE)**2
       ! simple initial setup for eon, further initialization will be done in chorus_EON
       do l = 1, NJ
          do k = -NZ, NZ
             eon(k,l)%zpos = zpos(k)
             eon(k,l)%Jact = Jact(k,l)
             do i = 0, Nq
                eon(k,l)%qcor(i) = -0.5_fp*Lq + i*dq
             end do
             do j = 0, Np
                eon(k,l)%pcor(j) = -0.5_fp*Lp + j*dp
             end do
             eon(k,l)%const = Jconst 
             eon(k,l)%omega = omega
             eon(k,l)%gyro = gyro(k)
             eon(k,l)%wp2 = wp2(k)
             eon(k,l)%kmode = kmode(k)
             eon(k,l)%pcr = pcr(k)
             eon(k,l)%drg = Jact(k,l)/kmode(k)*dgyro(k)-(omega-gyro(k))**2              &
                          / (kmode(k)*kmode(k)*kmode(k)*kmode(k))*dkmode(k)
             if (Jact(k,l)+pcr(k) > 0.0_fp)  then
                eon(k,l)%is_resonant = .TRUE.
             else
                eon(k,l)%is_resonant = .FALSE.
             end if
          end do
       end do
       ! drag force at mode location
       do l = 1, NJ
          do k = NZ, -NZ+1, -1
             eon(k,l)%drg = 0.5_fp*(eon(k,l)%drg+eon(k-1,l)%drg)
             eon(k,l)%Jpos = 0.5_fp*(eon(k,l)%Jact+eon(k-1,l)%Jact)
          end do
          eon(-NZ,l)%drg = 2.0_fp*eon(-NZ,l)%drg-eon(-NZ+1,l)%drg
          eon(-NZ,l)%Jpos = 2.0_fp*eon(-NZ,l)%Jact-eon(-NZ+1,l)%Jpos
       end do
       ! simple initial setup for chorus, further initialization will be done in chorus_WAVE
       forall(k=-NZ+1:NZ)  chorus(k)%zpos = 0.5_fp*(zpos(k)+zpos(k-1))
       ! two ghost points, the left one sets as a noise background level, the right one sets as
       ! an out-going boundary condition
       chorus(-NZ)%zpos  = 2.0_fp*zpos(-NZ)-chorus(-NZ+1)%zpos
       chorus(NZ+1)%zpos = 2.0_fp*zpos(NZ)-chorus(NZ)%zpos
       forall(k=-NZ:NZ+1)  chorus(k)%omega  = omega
       forall(k=-NZ:NZ+1)  chorus(k)%gyro = gyro0*(1.0_fp+4.5_fp*(chorus(k)%zpos/(Lshell*RE))**2)
       forall(k=-NZ:NZ+1)  chorus(k)%wp2  = 1.0_fp+cold2*(chorus(k)%zpos/(Lshell*RE))**2
       forall(k=-NZ:NZ+1)  chorus(k)%kmode = sqrt(omega**2+chorus(k)%wp2*omega      &
                                           / (chorus(k)%gyro-omega))
       forall(k=-NZ:NZ+1)  chorus(k)%pcr = (omega-chorus(k)%gyro)/chorus(k)%kmode**2   
    end if
    write(*,*) 'most unstable frequency (unit of gyro-frequency): ', omega/gyro0
    write(*,*) 'most unstable growth rate (unit of linear frequency): ', gamma/omega
    write(*,*) 'time step used in the simulation: ', dT
    return
  end subroutine set_init


  subroutine set_zpos(zpos, omega)
    ! set the zpos grid determined by the equilibrium
    implicit none
    real(fp), intent(out) :: zpos(-NZ:NZ)
    real(fp), intent(in) :: omega
    real(fp) :: LT
    integer ::  k
    real(fp) :: gyro, wp2, kmode
    real(fp) :: rk(4)    
    integer, parameter :: limit = 500
    integer, parameter :: lenw = 4*limit, key = 6
    real(fp) :: abserr
    real(fp), parameter :: epsabs = 0.0_fp, epsrel = 1.0e-6_fp
    integer :: ier, iwork(limit), last, neval
    real(fp) :: work(lenw)
    omega0 = omega
    ! transit time over the simulation range
    call dqag(transit_time, 0, -LZ, epsabs, epsrel, key, LT, abserr, neval, ier,     &
              limit, lenw, last, iwork, work)
    ! time step
    dT = abs(LT) / NZ 
    zpos(0) = 0.0_fp
    ! zpos grid 
    do k = 0, -NZ+1, -1
       ! RK1
       gyro = gyro0*(1.0_fp+4.5_fp*(zpos(k)/(Lshell*RE))**2)
       wp2 = 1.0_fp+cold2*(zpos(k)/(Lshell*RE))**2
       kmode = sqrt(omega**2+wp2*omega/(gyro-omega))
       rk(1) = (omega-gyro)/kmode
       zpos(k-1) = zpos(k) + 0.5_fp*dT*rk(1)
       ! RK2
       gyro = gyro0*(1.0_fp+4.5_fp*(zpos(k-1)/(Lshell*RE))**2)
       wp2 = 1.0_fp+cold2*(zpos(k-1)/(Lshell*RE))**2
       kmode = sqrt(omega**2+wp2*omega/(gyro-omega))
       rk(2) = (omega-gyro)/kmode 
       zpos(k-1) = zpos(k) + 0.5_fp*dT*rk(2)
       ! RK3
       gyro = gyro0*(1.0_fp+4.5_fp*(zpos(k-1)/(Lshell*RE))**2)
       wp2 = 1.0_fp+cold2*(zpos(k-1)/(Lshell*RE))**2
       kmode = sqrt(omega**2+wp2*omega/(gyro-omega))
       rk(3) = (omega-gyro)/kmode
       zpos(k-1) = zpos(k) + dT*rk(3)
       ! RK4
       gyro = gyro0*(1.0_fp+4.5_fp*(zpos(k-1)/(Lshell*RE))**2)
       wp2 = 1.0_fp+cold2*(zpos(k-1)/(Lshell*RE))**2
       kmode = sqrt(omega**2+wp2*omega/(gyro-omega))
       rk(4) = (omega-gyro)/kmode
       zpos(k-1) = zpos(k) + dT*(rk(1)+2.0_fp*(rk(2)+rk(3))+rk(4))/6.0_fp
    end do
    forall(k=1:NZ)  zpos(k) = -zpos(-k)
    return
  end subroutine set_zpos

  subroutine set_zpos_uni(zpos)
   ! set the zpos grid determined by the equilibrium
   implicit none
   real(fp), intent(out) :: zpos(-NZ:NZ)
   integer :: k
   zpos(-Nz) = -LZ
   ! zpos grid 
   do k = -NZ, NZ-1
      zpos(k+1) = zpos(k) +dZ
   end do
   return
 end subroutine set_zpos_uni


  subroutine set_Jact(Jact, zpos, omega,Jconst)
    ! set the Jact grid determined by the equilibrium
    ! J between [-pcr, -pcr+LJ]
    implicit none
    real(fp), intent(out) :: Jact(-NZ:NZ,NJ), Jconst
    real(fp), intent(in) :: zpos(-NZ:NZ)
    real(fp), intent(in) :: omega
    integer :: NSOL, INFO, LWA
    real(fp) :: TOL
    real(fp) :: SOL(1), FVEC(1), WA(180)
    real(fp) :: kmode0, pcr0
    real(fp) :: const, gyro, wp2, kmode
    integer :: k, l
    NSOL = 1
    LWA = 180
    TOL = small
    ! sample the Jact grid and establish nonunifrom Jact grid at the equator
    if (NJ == 1)  then
       Jact(0,1) = vperp*vperp/gyro0
    else
       SOL(1) = 0.1_fp*vperp*vperp/gyro0
       do l = 1, NJ
          ratio = (real(l,fp)-0.1_fp) / NJ
          call hybrd1(sample_dist, NSOL, SOL, FVEC, TOL, INFO, WA, LWA)
          Jact(0,l) = SOL(1)
       end do
    end if
    kmode0 = sqrt(omega**2+omega/(gyro0-omega))
    pcr0 = (omega-gyro0)/(kmode0*kmode0)
    Jact(0,:) = Jact(0,:) - pcr0
    ! global Jact grid
    do l = 1, NJ
       const=(omega-gyro0)*Jact(0,l)+0.5_fp*(omega-gyro0)**2/kmode0**2
       Jconst = const
       do k = -NZ, -1
          gyro = gyro0*(1.0_fp+4.5_fp*(zpos(k)/(Lshell*RE))**2)
          wp2 = 1.0_fp+cold2*(zpos(k)/(Lshell*RE))**2
          kmode = sqrt(omega**2+wp2*omega/(gyro-omega))         
          Jact(k,l) = (const-0.5_fp*(omega-gyro)**2/kmode**2)/(omega-gyro)        
       end do
       do k = 1, NZ
          gyro = gyro0*(1.0_fp+4.5_fp*(zpos(k)/(Lshell*RE))**2)
          wp2 = 1.0_fp+cold2*(zpos(k)/(Lshell*RE))**2
          kmode = sqrt(omega**2+wp2*omega/(gyro-omega))         
          Jact(k,l) = (const-0.5_fp*(omega-gyro)**2/kmode**2)/(omega-gyro)        
       end do
    end do
    return
  end subroutine set_Jact


  subroutine get_most_unstable(omega, gamma)
  ! determine the most unstable mode triggered at the equator  
  ! the real frequency of most unstable mode: omega
  ! the linear growth rate of most unstable mode: gamma
    implicit none
    real(fp), intent(out) :: omega, gamma
    integer :: NSOL, INFO, LWA
    real(fp) :: TOL
    real(fp) :: SOL(1), FVEC(1), WA(180)
    real(fp) :: dgamma(1)
    NSOL = 1
    LWA = 180
    TOL = small
    SOL(1) = 0.35_fp*gyro0
    call hybrd1(get_dgamma, NSOL, SOL, FVEC, TOL, INFO, WA, LWA)
    omega = SOL(1)
    call get_growth_rate(omega, gamma, dgamma, 1)
    if (gamma .le. 0.0_fp) then
       call write_log('The growth rate has to be positive in the source region.', gamma)
       stop
    end if   
    return
  end subroutine get_most_unstable  


  subroutine get_dgamma(NSOL, SOL, FVEC, IFLAG)
    ! calculate dgamma/domega as a function of omega at the equator
    implicit none
    integer, intent(in) :: NSOL
    integer, intent(inout) :: IFLAG
    real(fp), intent(in) :: SOL(NSOL)
    real(fp), intent(out) :: FVEC(NSOL)
    real(fp) :: omega, gamma, dgamma(1)
    omega = SOL(1)
    call get_growth_rate(omega, gamma, dgamma, 1)
    FVEC(1) = dgamma(1)
    return
  end subroutine get_dgamma

  
  subroutine get_growth_rate(omega, gamma, dgamma, deriv)
    ! growth rate is the largest at the equator
    ! auto_deriv subroutine is used to automatically calculate dgamma/domega
    implicit none
    real(fp), intent(in) :: omega
    real(fp), intent(out) :: gamma
    real(fp), intent(out) :: dgamma(1)
    integer,  intent(in) :: deriv
    type(func) :: growth_rate, wvar, kvar, vg_var 
    call derivative(deriv)
    call independent(1, wvar, omega)
    kvar  = sqrt(wvar**2+wvar/(gyro0-wvar))
    vg_var = 2.0_fp*kvar/(2.0_fp*wvar+gyro0/(gyro0-wvar)**2)
    growth_rate = sqrt(2.0_fp*pi)*gyro0*vg_var*nh/(4.0_fp*kvar**2*vll)          &
                * exp (-0.5_fp*(wvar-gyro0)**2/(kvar*vll)**2)                   &
                * ((1.0_fp+loss_cone)*aniso*(gyro0-wvar)/gyro0-1.0_fp)   
    call extract(growth_rate, gamma, dgamma)
    return
  end subroutine get_growth_rate


  function transit_time(zpos)
    implicit none
    real(fp) :: transit_time
    real(fp), intent(in) :: zpos
    real(fp) :: gyro, wp2, kmode
    gyro = gyro0*(1.0_fp+4.5_fp*(zpos/(Lshell*RE))**2)
    wp2 = 1.0_fp+cold2*(zpos/(Lshell*RE))**2
    kmode = sqrt(omega0**2+wp2*omega0/(gyro-omega0))
    transit_time = kmode/(omega0-gyro)
    return
  end function transit_time

  
  subroutine sample_dist(NSOL, SOL, FVEC, IFLAG)
    implicit none
    integer, intent(in) :: NSOL
    integer, intent(inout) :: IFLAG
    real(fp), intent(in) :: SOL(NSOL)
    real(fp), intent(out) :: FVEC(NSOL)
    real(fp) :: mu
    mu = SOL(1)
    if(loss_cone /= 0 ) then
    FVEC(1) = 1.0_fp-ratio-(exp(-mu*gyro0/(vperp*vperp))-loss_cone              &
            * exp(-mu*gyro0/(loss_cone*vperp*vperp)))/(1.0_fp-loss_cone)  
    else
    FVEC(1) = 1.0_fp-ratio-(exp(-mu*gyro0/(vperp*vperp)))
    end if
    return
  end subroutine sample_dist

end module chorus_INIT


!!$program test_chorus_INIT
!!$  use chorus_IO
!!$  use chorus_TYPE
!!$  use chorus_LIB
!!$  use chorus_INIT
!!$  implicit none
!!$  type(energetic_eon) :: eon(-NZ:NZ, NJ)
!!$  type(chorus_mode) :: chorus(-NZ:NZ+1)
!!$  real(fp) :: omega, gamma
!!$  real(fp) :: zpos(-NZ:NZ), Jact(-NZ:NZ, NJ)
!!$  call get_most_unstable(omega, gamma)
!!$  call set_init(eon, chorus)
!!$  call save_output('zpos.out', chorus%zpos, 2*NZ+2, .FALSE.)
!!$  call save_output('gyro.out', chorus%gyro, 2*NZ+2,  .FALSE.)
!!$  call save_output('wp2.out', chorus%wp2, 2*NZ+2, .FALSE.)
!!$  call save_output('kmode.out', chorus%kmode, 2*NZ+2, .FALSE.)
!!$  call save_output('pcr.out', chorus%pcr, 2*NZ+2, .FALSE.)
!!$  call save_output('drg.out', eon%drg, 2*NZ+1, NJ, .FALSE.)
!!$  call save_output('Jact.out', eon%Jact, 2*NZ+1, NJ, .FALSE.)
!!$  stop 
!!$end program test_chorus_INIT
