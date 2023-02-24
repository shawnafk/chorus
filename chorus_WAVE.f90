module chorus_WAVE
   use chorus_IO
   use chorus_TYPE
   use chorus_LIB
   use chorus_INIT
   implicit none
   private

   public :: init_chorus, update_chorus, update_chorus_adv, intgl_source 

   complex(fp) :: lft_bnd

contains
   subroutine init_chorus(chorus)
      ! chorus initialization
      implicit none
      type(chorus_mode), intent(out) :: chorus(-NZ:NZ+1)
      real(fp) :: omega, gamma
      integer :: k
      call get_most_unstable(omega, gamma)
      lft_bnd = noise*gamma*gamma*(1.0_fp+mathj)
      if (.not. RESTART)  then
         ! initial amplitude
         ! ampl(1) = a
         ! ampl(2) = da/dt
         ! ampl(3) = int_0^t dtau da/dtau*exp(-mathj*(omega-gyro)*(t-tau))
         ! noise background: normal distributed noise
         do k = -NZ, NZ+1
            chorus(k)%ampl(1) = lft_bnd*(1.0_fp+noise*norm_rand(0.0_fp,1.0_fp))
            chorus(k)%ampl(2) = (0.0_fp, 0.0_fp)
            chorus(k)%ampl(3) = (0.0_fp, 0.0_fp)
         end do
      end if
      return
   end subroutine init_chorus


   subroutine update_chorus(chorus)
      ! solve the wave amplitude complex equations
      implicit none
      type(chorus_mode), intent(inout) :: chorus(-NZ:NZ+1)
      complex(fp) :: ampl(-NZ:NZ+1,3)
      real(fp) :: tstp(0:1), dtstp, err
      integer :: k
      ! update chorus amplitude
      tstp(0) = 0.0_fp
      tstp(1) = 0.1_fp*dT
      dtstp = tstp(1) - tstp(0)
      do while (.TRUE.)
         ! adaptive loop
         do while (.TRUE.)
            forall(k=-NZ:NZ+1)  ampl(k,1:3) = chorus(k)%ampl(1:3)
            call rk23_adstep(ampl, chorus, dtstp, err)
            err = err / rtolw
            ! tstp is too large
            if (err .gt. 1.0_fp)  then
               if (abs(dtstp) < dtminw) exit
               dtstp = max(abs(safety*(err**pshrink)*dtstp), 0.1_fp*dtstp)
               tstp(1) = tstp(0) + dtstp
               ! tstp is too small
            else
               if (err > errcon) then
                  dtstp = safety*(err**pgrow)*dtstp
               else
                  dtstp = 5.0_fp*dtstp
               end if
               exit
            end if
         end do
         forall(k=-NZ:NZ+1)  chorus(k)%ampl(1:3) = ampl(k,1:3)
         tstp(0) = tstp(1)
         tstp(1) = tstp(1) + dtstp
         ! reach the final position and exit the loop
         if (tstp(0) .ge. dT) then
            exit
         elseif (tstp(1) > dT) then
            tstp(1) = dT
            dtstp = tstp(1) - tstp(0)
            call rk23_adstep(ampl, chorus, dtstp, err)
            forall(k=-NZ:NZ+1)  chorus(k)%ampl(1:3) = ampl(k,1:3)
            exit
         end if
      end do
      return
   end subroutine update_chorus

   subroutine rk23_adstep(ampl, chorus, dtstp, err)
      ! update the chorus amplitude a each step
      implicit none
      complex(fp), intent(inout) :: ampl(-NZ:NZ+1,3)
      type(chorus_mode), intent(in) :: chorus(-NZ:NZ+1)
      real(fp), intent(in) :: dtstp
      real(fp), intent(out) :: err
      complex(fp) :: y(-NZ:NZ+1,3), dyodt(-NZ+1:NZ,3), rk(-NZ+1:NZ,3), enorm(-NZ+1:NZ,3)
      complex(fp) :: fluxb, rkb
      real(fp) :: vg
      integer :: k
      vg = 2.0_fp*chorus(NZ+1)%kmode/(2.0_fp*chorus(NZ+1)%omega                          &
         + chorus(NZ+1)%wp2*chorus(NZ+1)%gyro/(chorus(NZ+1)%gyro-chorus(NZ+1)%omega)**2)
      !!!!!!!!!!!!!!!!!!!!!!!! Runge-Kutta substep 1 !!!!!!!!!!!!!!!!!!!!!!!!!!!
      y = ampl
      dyodt = (0.0_fp, 0.0_fp)
      fluxb = -vg*(3.0_fp*y(NZ+1,1)-4.0_fp*y(NZ,1)+y(NZ-1,1))                            &
         / (chorus(NZ+1)%zpos-chorus(NZ-1)%zpos)
      do k = -NZ+1, NZ
         dyodt(k,1) = y(k,2) + diff_num                                                  &
            * (2.0_fp*y(k-1,1)/((chorus(k-1)%zpos-chorus(k)%zpos)                &
            * (chorus(k-1)%zpos-chorus(k+1)%zpos))+2.0_fp*y(k,1)                 &
            / ((chorus(k)%zpos-chorus(k+1)%zpos)                                 &
            * (chorus(k)%zpos-chorus(k-1)%zpos))                                 &
            + 2.0_fp*y(k+1,1)/((chorus(k+1)%zpos-chorus(k)%zpos)                 &
            * (chorus(k+1)%zpos-chorus(k-1)%zpos)))
         dyodt(k,2) = -2.0_fp*mathj*chorus(k)%omega*y(k,2)                               &
            + 2.0_fp*y(k-1,1)/((chorus(k-1)%zpos-chorus(k)%zpos)                 &
            * (chorus(k-1)%zpos-chorus(k+1)%zpos))+2.0_fp*y(k,1)                 &
            / ((chorus(k)%zpos-chorus(k+1)%zpos)                                 &
            * (chorus(k)%zpos-chorus(k-1)%zpos))                                 &
            + 2.0_fp*y(k+1,1)/((chorus(k+1)%zpos-chorus(k)%zpos)                 &
            * (chorus(k+1)%zpos-chorus(k-1)%zpos))-2.0_fp*mathj*chorus(k)%kmode  &
            * ((y(k,1)-y(k-1,1))/(chorus(k)%zpos-chorus(k-1)%zpos)               &
            + (y(k,1)-y(k+1,1))/(chorus(k)%zpos-chorus(k+1)%zpos)                &
            - (y(k+1,1)-y(k-1,1))/(chorus(k+1)%zpos-chorus(k-1)%zpos))           &
            - chorus(k)%wp2*chorus(k)%gyro/(chorus(k)%gyro-chorus(k)%omega)      &
            * y(k,3)-chorus(k)%Sr
         dyodt(k,3) = y(k,2)-mathj*(chorus(k)%omega-chorus(k)%gyro)*y(k,3)
         rk(k,1) = dyodt(k,1)
         rk(k,2) = dyodt(k,2)
         rk(k,3) = dyodt(k,3)
         y(k,1) = ampl(k,1) + rk(k,1)*dtstp
         y(k,2) = ampl(k,2) + rk(k,2)*dtstp
         y(k,3) = ampl(k,3) + rk(k,3)*dtstp
      end do
      rkb = fluxb
      y(NZ+1,1) = ampl(NZ+1,1) + rkb*dtstp
      y(-NZ,1) = ampl(-NZ,1)
      !!!!!!!!!!!!!!!!!!!!!!!! Runge-Kutta substep 2 !!!!!!!!!!!!!!!!!!!!!!!!!!!
      fluxb = -vg*(3.0_fp*y(NZ+1,1)-4.0_fp*y(NZ,1)+y(NZ-1,1))                            &
         / (chorus(NZ+1)%zpos-chorus(NZ-1)%zpos)
      dyodt = (0.0_fp, 0.0_fp)
      do k = -NZ+1, NZ
         dyodt(k,1) = y(k,2) + diff_num                                                  &
            * (2.0_fp*y(k-1,1)/((chorus(k-1)%zpos-chorus(k)%zpos)                &
            * (chorus(k-1)%zpos-chorus(k+1)%zpos))+2.0_fp*y(k,1)                 &
            / ((chorus(k)%zpos-chorus(k+1)%zpos)                                 &
            * (chorus(k)%zpos-chorus(k-1)%zpos))                                 &
            + 2.0_fp*y(k+1,1)/((chorus(k+1)%zpos-chorus(k)%zpos)                 &
            * (chorus(k+1)%zpos-chorus(k-1)%zpos)))
         dyodt(k,2) = -2.0_fp*mathj*chorus(k)%omega*y(k,2)                               &
            + 2.0_fp*y(k-1,1)/((chorus(k-1)%zpos-chorus(k)%zpos)                 &
            * (chorus(k-1)%zpos-chorus(k+1)%zpos))+2.0_fp*y(k,1)                 &
            / ((chorus(k)%zpos-chorus(k+1)%zpos)                                 &
            * (chorus(k)%zpos-chorus(k-1)%zpos))                                 &
            + 2.0_fp*y(k+1,1)/((chorus(k+1)%zpos-chorus(k)%zpos)                 &
            * (chorus(k+1)%zpos-chorus(k-1)%zpos))-2.0_fp*mathj*chorus(k)%kmode  &
            * ((y(k,1)-y(k-1,1))/(chorus(k)%zpos-chorus(k-1)%zpos)               &
            + (y(k,1)-y(k+1,1))/(chorus(k)%zpos-chorus(k+1)%zpos)                &
            - (y(k+1,1)-y(k-1,1))/(chorus(k+1)%zpos-chorus(k-1)%zpos))           &
            - chorus(k)%wp2*chorus(k)%gyro/(chorus(k)%gyro-chorus(k)%omega)      &
            * y(k,3)-chorus(k)%Sr
         dyodt(k,3) = y(k,2)-mathj*(chorus(k)%omega-chorus(k)%gyro)*y(k,3)
         rk(k,1) = rk(k,1) + dyodt(k,1)
         rk(k,2) = rk(k,2) + dyodt(k,2)
         rk(k,3) = rk(k,3) + dyodt(k,3)
         y(k,1) = ampl(k,1) + 0.25_fp*rk(k,1)*dtstp
         y(k,2) = ampl(k,2) + 0.25_fp*rk(k,2)*dtstp
         y(k,3) = ampl(k,3) + 0.25_fp*rk(k,3)*dtstp
      end do
      rkb = rkb + fluxb
      y(NZ+1,1) = ampl(NZ+1,1) + 0.25_fp*rkb*dtstp
      y(-NZ,1) = ampl(-NZ,1)
      enorm = rk
      !!!!!!!!!!!!!!!!!!!!!!!! Runge-Kutta substep 3 !!!!!!!!!!!!!!!!!!!!!!!!!!!
      fluxb = -vg*(3.0_fp*y(NZ+1,1)-4.0_fp*y(NZ,1)+y(NZ-1,1))                            &
         / (chorus(NZ+1)%zpos-chorus(NZ-1)%zpos)
      dyodt = (0.0_fp, 0.0_fp)
      do k = -NZ+1, NZ
         dyodt(k,1) = y(k,2) + diff_num                                                  &
            * (2.0_fp*y(k-1,1)/((chorus(k-1)%zpos-chorus(k)%zpos)                &
            * (chorus(k-1)%zpos-chorus(k+1)%zpos))+2.0_fp*y(k,1)                 &
            / ((chorus(k)%zpos-chorus(k+1)%zpos)                                 &
            * (chorus(k)%zpos-chorus(k-1)%zpos))                                 &
            + 2.0_fp*y(k+1,1)/((chorus(k+1)%zpos-chorus(k)%zpos)                 &
            * (chorus(k+1)%zpos-chorus(k-1)%zpos)))
         dyodt(k,2) = -2.0_fp*mathj*chorus(k)%omega*y(k,2)                               &
            + 2.0_fp*y(k-1,1)/((chorus(k-1)%zpos-chorus(k)%zpos)                 &
            * (chorus(k-1)%zpos-chorus(k+1)%zpos))+2.0_fp*y(k,1)                 &
            / ((chorus(k)%zpos-chorus(k+1)%zpos)                                 &
            * (chorus(k)%zpos-chorus(k-1)%zpos))                                 &
            + 2.0_fp*y(k+1,1)/((chorus(k+1)%zpos-chorus(k)%zpos)                 &
            * (chorus(k+1)%zpos-chorus(k-1)%zpos))-2.0_fp*mathj*chorus(k)%kmode  &
            * ((y(k,1)-y(k-1,1))/(chorus(k)%zpos-chorus(k-1)%zpos)               &
            + (y(k,1)-y(k+1,1))/(chorus(k)%zpos-chorus(k+1)%zpos)                &
            - (y(k+1,1)-y(k-1,1))/(chorus(k+1)%zpos-chorus(k-1)%zpos))           &
            - chorus(k)%wp2*chorus(k)%gyro/(chorus(k)%gyro-chorus(k)%omega)      &
            * y(k,3)-chorus(k)%Sr
         dyodt(k,3) = y(k,2)-mathj*(chorus(k)%omega-chorus(k)%gyro)*y(k,3)
         rk(k,1) = rk(k,1) + 4.0_fp*dyodt(k,1)
         rk(k,2) = rk(k,2) + 4.0_fp*dyodt(k,2)
         rk(k,3) = rk(k,3) + 4.0_fp*dyodt(k,3)
         y(k,1) = ampl(k,1) + rk(k,1)*dtstp/6.0_fp
         y(k,2) = ampl(k,2) + rk(k,2)*dtstp/6.0_fp
         y(k,3) = ampl(k,3) + rk(k,3)*dtstp/6.0_fp
      end do
      rkb = rkb + 4.0_fp*fluxb
      y(NZ+1,1) = ampl(NZ+1,1) + rkb*dtstp/6.0_fp
      ! noisy left boundary
      y(-NZ,1) = lft_bnd*(1.0_fp+noise*norm_rand(0.0_fp,1.0_fp))
      forall(k=-NZ:NZ+1)  ampl(k,1:3) = y(k,1:3)
      enorm = abs(1.5_fp*enorm-0.5_fp*rk) * dtstp/3.0_fp
      err = maxval(abs(enorm)/(abs(y(-NZ+1:NZ,1:3))+TINY_NUM))
      return
   end subroutine rk23_adstep

   subroutine update_chorus_adv(chorus)
      ! solve the wave amplitude complex equations
      implicit none
      type(chorus_mode), intent(inout) :: chorus(-NZ:NZ+1)
      !update chorus amplitude
      !call lax_wendroff(chorus, dT)
      call implicit_upwind(chorus,dT)
      !call upwind_1st(chorus, dT)
      return
   end subroutine update_chorus_adv


   subroutine implicit_upwind(chorus,dtstp)
      implicit none
      real(fp), intent(in) :: dtstp
      type(chorus_mode), intent(inout) :: chorus(-NZ:NZ+1)
      real(fp), parameter :: implicit = 1.0_fp
      real(fp) :: deltaZ, vg(-NZ:NZ+1)
      complex(fp) :: ampl(-NZ:NZ+1)
      integer :: k
      forall(k = -NZ:NZ+1)  vg(k) = 2.0_fp*chorus(k)%kmode                         &
         / (2.0_fp*chorus(k)%omega+chorus(k)%gyro         &
         / (chorus(k)%gyro-chorus(k)%omega)**2)
      forall(k = -NZ:NZ+1)  ampl(k) = chorus(k)%ampl(1)
      do k = -NZ+1, NZ+1
         deltaZ = chorus(k)%zpos-chorus(k-1)%zpos
         chorus(k)%ampl(1) = (0.5_fp*mathj*vg(k)/chorus(k)%kmode*chorus(k)%Sr    &
            + implicit*vg(k)/deltaZ*chorus(k-1)%ampl(1)+(1.0_fp/dtstp    &
            - (1.0_fp-implicit)*vg(k)/deltaZ)*ampl(k)                 &
            + (1.0_fp-implicit)*vg(k)/deltaZ*ampl(k-1))               &
            / (1.0_fp/dtstp+implicit*vg(k)/deltaZ)
      end do
      chorus(-NZ)%ampl(1)= lft_bnd*(1.0_fp+noise*norm_rand(0.0_fp,1.0_fp))
      return
   end subroutine implicit_upwind

   subroutine upwind_1st(chorus, dtstp)
      ! update the chorus amplitude a each step
      implicit none
      type(chorus_mode), intent(inout) :: chorus(-NZ:NZ+1)
      real(fp), intent(in) :: dtstp
      complex(fp) :: Sr
      real(fp) :: vg, deltaZ
      integer :: k
      complex(fp) :: ampl_pre(-NZ:NZ+1)
      forall(k=-NZ:NZ+1)  ampl_pre(k) = chorus(k)%ampl(1)
      !chorus from -NZ NZ+1, here we have left -NZ and NZ+1
      do k = -NZ+1, NZ
         deltaZ = chorus(k-1)%zpos - chorus(k)%zpos
         vg = 2.0_fp*chorus(k-1)%kmode/(2.0_fp*chorus(k-1)%omega + chorus(k-1)%wp2*chorus(k-1)%gyro/(chorus(k-1)%gyro-chorus(k-1)%omega)**2)
         Sr  = 0.5_fp*mathj*vg/chorus(k)%kmode*chorus(k-1)%Sr
         chorus(k)%ampl(1) = ampl_pre(k) - vg*(ampl_pre(k-1) - ampl_pre(k))/deltaZ*dtstp + Sr*dtstp
      end do
      chorus(-NZ)%ampl(1)= lft_bnd*(1.0_fp+noise*norm_rand(0.0_fp,1.0_fp))
   end subroutine upwind_1st

   subroutine lax_wendroff(chorus, dtstp)
      ! update the chorus amplitude a each step
      implicit none
      type(chorus_mode), intent(inout) :: chorus(-NZ:NZ+1)
      real(fp), intent(in) :: dtstp
      complex(fp) :: chorus_left_half, chorus_right_half,Sr_right_half,Sr_left_half,Sr_ave
      complex(fp) :: convert_coef
      real(fp) :: vgk(-NZ:NZ+1),vg_l, vg_r, vg_ave,vgb
      real(fp) :: dz_l, dz_r, dz_ave
      integer :: k
      complex(fp) :: ampl_pre(-NZ:NZ+1), fluxb
      do k = -NZ,NZ+1
         ampl_pre(k) = chorus(k)%ampl(1)
         vgk(k) = 2.0_fp*chorus(k)%kmode/(2.0_fp*chorus(k)%omega + chorus(k)%wp2*chorus(k)%gyro/(chorus(k)%gyro-chorus(k)%omega)**2)
      end do
      !chorus from -NZ NZ+1, here we have left -NZ and NZ+1
      do k = -NZ+1, NZ
         dz_l = chorus(k)%zpos - chorus(k-1)%zpos
         vg_l = (vgk(k-1)+vgk(k))/2
         Sr_left_half = (0.5_fp*mathj*vgk(k)/chorus(k)%kmode*chorus(k)%Sr + 0.5_fp*mathj*vgk(k-1)/chorus(k-1)%kmode*chorus(k-1)%Sr)/2
         chorus_left_half = (ampl_pre(k) + ampl_pre(k-1))/2 -vg_l*(ampl_pre(k)-ampl_pre(k-1))/dz_l *dtstp/2 + Sr_left_half*dtstp/2
         dz_r = chorus(k+1)%zpos - chorus(k)%zpos
         vg_r = (vgk(k+1)+vgk(k))/2
         Sr_right_half = (0.5_fp*mathj*vgk(k)/chorus(k)%kmode*chorus(k)%Sr + 0.5_fp*mathj*vgk(k+1)/chorus(k+1)%kmode*chorus(k+1)%Sr)/2
         chorus_right_half = (ampl_pre(k) + ampl_pre(k+1))/2 -vg_r*(ampl_pre(k+1)-ampl_pre(k))/dz_r *dtstp/2 + Sr_right_half*dtstp/2
         dz_ave = (dz_l + dz_r)/2
         vg_ave = (vg_l + vg_r)/2
         Sr_ave  = (Sr_right_half+Sr_left_half)/2
         chorus(k)%ampl(1) = ampl_pre(k) - vg_ave*(chorus_right_half - chorus_left_half)/dz_ave*dtstp + Sr_ave*dtstp
      end do
      vgb = 2.0_fp*chorus(NZ+1)%kmode/(2.0_fp*chorus(NZ+1)%omega                          &
         + chorus(NZ+1)%wp2*chorus(NZ+1)%gyro/(chorus(NZ+1)%gyro-chorus(NZ+1)%omega)**2)
      fluxb = -vgb*(3.0_fp*ampl_pre(NZ+1)-4.0_fp*ampl_pre(NZ)+ampl_pre(NZ-1))                            &
         / (chorus(NZ+1)%zpos-chorus(NZ-1)%zpos)
      !minus at fluxb
      chorus(NZ+1)%ampl(1) = chorus(NZ+1)%ampl(1) + fluxb*dtstp

      chorus(-NZ)%ampl(1)= lft_bnd*(1.0_fp+noise*norm_rand(0.0_fp,1.0_fp))
   end subroutine lax_wendroff

   subroutine intgl_source(chorus, eon)
      ! chorus%Sr = kmode*nh*\iint dqcor*dpcor*dJact
      !           * sqrt(2*gyro*(Jact+pcr+pcor))*feon*exp(mathj*qcor)
      implicit none
      type(chorus_mode), intent(inout) :: chorus(-NZ:NZ+1)
      type(energetic_eon), intent(in) :: eon(-NZ:NZ,NJ)
      complex(fp) :: jr(-NZ:NZ,NJ), Sr(-NZ:NZ)
      real(fp) :: fdp(0:Nq), fdpdq(0:Nq-1), fvp(0:Np-1,2)
      complex(fp) :: fqcor(0:Nq-1)
      integer :: i, j, k, l
      do l = 1, NJ
         do k = -NZ, NZ
            if ( eon(k,l)%is_resonant )   then
               ! fdp = \int dpcor sqrt(2*gyro*(Jact+pcr+pcor))*feon
               do i = 0, Nq
                 do j = 0, Np-1
                     fvp(j,1) = line_intg(sqrt(2.0_fp*eon(k,l)%gyro*max(eon(k,l)%Jact        &
                              + eon(k,l)%pcr+eon(k,l)%pcor(j), 0.0_fp)),                     &
                              sqrt(2.0_fp*eon(k,l)%gyro*max(eon(k,l)%Jact+eon(k,l)%pcr       &
                              + 0.5_fp*(eon(k,l)%pcor(j)+eon(k,l)%pcor(j+1)), 0.0_fp)),      &
                              sqrt(2.0_fp*eon(k,l)%gyro*max(eon(k,l)%Jact+eon(k,l)%pcr       &
                              + eon(k,l)%pcor(j+1), 0.0_fp)), eon(k,l)%feon(i,j,1),          &
                              eon(k,l)%feon(i,j,3), eon(k,l)%feon(i,j+1,1), dp)
                  end do
                  fdp(i) = sum(fvp(0:Np-1,1))
               end do
               ! fdpdq = \int_q^{q+dq} dqcor \int dpcor sqrt(2*gyro*(Jact+pcr+pcor))*feon
               do i = 0, Nq-1
                  do j = 0, Np-1
                     fvp(j,2) = line_intg(sqrt(2.0_fp*eon(k,l)%gyro*max(eon(k,l)%Jact        &
                              + eon(k,l)%pcr+eon(k,l)%pcor(j), 0.0_fp)),                     &
                              sqrt(2.0_fp*eon(k,l)%gyro*max(eon(k,l)%Jact+eon(k,l)%pcr       &
                              + 0.5_fp*(eon(k,l)%pcor(j)+eon(k,l)%pcor(j+1)), 0.0_fp)),      &
                              sqrt(2.0_fp*eon(k,l)%gyro*max(eon(k,l)%Jact+eon(k,l)%pcr       &
                              + eon(k,l)%pcor(j+1), 0.0_fp)), eon(k,l)%feon(i,j,2),          &
                              eon(k,l)%feon(i,j,4), eon(k,l)%feon(i,j+1,2), dp)
                  end do
                  fdpdq(i) = sum(fvp(0:Np-1,2))
               end do
               ! jr = kmode * \int dqcor*dpcor sqrt(2*gyro*(Jact+pcr+pcor))*feon*exp(mathj*\xi)
               do i = 0, Nq-1
                  fqcor(i) = line_intg(cos(eon(k,l)%qcor(i)),                                &
                             cos(0.5_fp*(eon(k,l)%qcor(i)+eon(k,l)%qcor(i+1))),              &
                             cos(eon(k,l)%qcor(i+1)), fdp(i), fdpdq(i), fdp(i+1), dq)        &
                           + mathj*line_intg(sin(eon(k,l)%qcor(i)),                          &
                             sin(0.5_fp*(eon(k,l)%qcor(i)+eon(k,l)%qcor(i+1))),              &
                             sin(eon(k,l)%qcor(i+1)), fdp(i), fdpdq(i), fdp(i+1), dq)
               end do
               jr(k,l) = eon(k,l)%kmode*sum(fqcor)
            else
               jr(k,l) = (0.0_fp, 0.0_fp)
            end if
         end do
      end do
      if (NJ == 1)  then
         forall(k=-NZ:NZ) Sr(k) = vperp*vperp/gyro0*exp(1.0_fp)*jr(k,1)
      else
         forall(k=-NZ:NZ) Sr(k) = trapez(eon(k,1:NJ)%Jact, jr(k,1:NJ), NJ)                   &
                                + 0.5_fp*(eon(k,1)%Jact+eon(k,1)%pcr)*jr(k,1)
      end if
      forall(k = -NZ+1:NZ)   chorus(k)%Sr = 0.5_fp*nh*(Sr(k)+Sr(k-1))
      chorus(-NZ)%Sr = (0.0_fp, 0.0_fp)
      chorus(NZ+1)%Sr = (0.0_fp, 0.0_fp)
      return
  end subroutine



   end module chorus_WAVE


   !!$program test_chorus_WAVE
   !!$  use chorus_IO
   !!$  use chorus_TYPE
   !!$  use chorus_LIB
   !!$  use chorus_INIT
   !!$  use chorus_WAVE
   !!$  implicit none
   !!$  type(chorus_mode) :: chorus(-NZ:NZ+1)
   !!$  type(energetic_eon) :: eon(-NZ:NZ,NJ)
   !!$  integer :: tstep, i, j, k, l
   !!$  call set_init(eon, chorus)
   !!$  call init_chorus(chorus)
   !!$  call save_output('chorus.out', chorus%ampl(1), 2*NZ+2, .FALSE.)
   !!$  do tstep = NT0+1, NT
   !!$     call update_chorus(chorus, eon)
   !!$     write(*,*) tstep
   !!$     if (mod(tstep,5) == 0) then
   !!$        call save_output('chorus.out', chorus%ampl(1), 2*NZ+2, .TRUE.)
   !!$     end if
   !!$  end do
   !!$  stop
!!$end program test_chorus_WAVE
