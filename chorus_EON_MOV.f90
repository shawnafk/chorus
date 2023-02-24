module chorus_EON_MOV
   use chorus_IO
   use chorus_TYPE
   use chorus_LIB
   use chorus_INIT
   use chorus_WAVE
   implicit none
   private

   public :: init_eon_mov, trap_eon_mov, shift_eon_mov, interp_ampl, interp_intf, intgl_source_mov

contains
   subroutine init_eon_mov(eon)
      implicit none
      type(energetic_eon), intent(inout) :: eon(-NZ:NZ,NJ)
      integer :: k, l
      if ( .not. RESTART)  then
         do l = 1, NJ
            do k = -NZ, NZ
               call set_feq(eon(k,l))
            end do
         end do
      end if
      return
   end subroutine init_eon_mov

   subroutine trap_eon_mov(eon)
      ! push one macro-particle
      implicit none
      type(energetic_eon), intent(inout) :: eon
      real(fp) :: tstp(0:1), dtstp, err
      real(fp) :: fstp(0:Nq,0:Np,4)
      if ( .not. eon%is_resonant )  then
         eon%feon = 0.0_fp
         return
      end if
      call get_Hamiltonian(eon)
      call build_eon(eon)
      tstp(0) = 0.0_fp
      tstp(1) = min(dT, 4.0_fp*pi/(Lp*Nq))
      dtstp = tstp(1) - tstp(0)
      do while (.TRUE.)
         ! adaptive loop
         do while (.TRUE.)
            fstp = eon%feon
            call rk23_adstep(fstp, dtstp, eon%dH_q, eon%dH_p, eon%Sr, err)
            err = err / rtolp
            ! tstp is too large
            if (err .gt. 1.0_fp) then
               if (abs(dtstp) < dtminp) exit
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
         eon%feon = fstp
         tstp(0) = tstp(1)
         tstp(1) = tstp(1) + dtstp
         ! reach the final position and exit the loop
         if (tstp(0) .ge. dT) then
            exit
         elseif (tstp(1) > dT)  then
            tstp(1) = dT
            dtstp = tstp(1) - tstp(0)
            call rk23_adstep(eon%feon, dtstp, eon%dH_q, eon%dH_p, eon%Sr, err)
            exit
         end if
      end do
      return
   end subroutine trap_eon_mov

   subroutine set_feq(eon)
      ! energetic electron equilibrium distribution for one macro-particle
      implicit none
      type(energetic_eon), intent(inout) :: eon
      real(fp) :: mu, Jdist1, Jdist2
      integer :: j
      do j = 0, Np
         mu = max(eon%Jpos+eon%pcr+eon%pcor(j), small)
         Jdist1 = exp(-mu*gyro0/(vperp*vperp))
         if (loss_cone /= 0.0_fp ) then
            Jdist2 = exp(-mu*gyro0/(loss_cone*vperp*vperp))
            eon%feq(j) = gyro0/((2.0_fp*pi)**1.5_fp*vperp*vperp*vll)/(1.0_fp-loss_cone)        &
               * exp(-0.5_fp*(eon%kmode*(eon%pcor(j)+eon%pcr)/vll)**2)           &
               * exp(-mu*(eon%gyro-gyro0)/(vll*vll))*(Jdist1-Jdist2)
            eon%dfeq(j) = -eon%feq(j)*((eon%kmode**2*eon%pcor(j)+eon%omega-gyro0)/(vll*vll)+gyro0/(vperp*vperp*loss_cone)*(loss_cone*Jdist1-Jdist2)/(Jdist1-Jdist2))
         else
            Jdist2 = 0.0_fp
            eon%feq(j) = gyro0/((2.0_fp*pi)**1.5_fp*vperp*vperp*vll)/(1.0_fp)        &
               * exp(-0.5_fp*(eon%kmode*(eon%pcor(j)+eon%pcr)/vll)**2)           &
               * exp(-mu*(eon%gyro-gyro0)/(vll*vll))*(Jdist1-Jdist2)
            eon%dfeq(j) = -eon%feq(j)*((eon%kmode**2*eon%pcor(j)+eon%omega-gyro0)        &
               / (vll*vll)+gyro0/(vperp*vperp)*(Jdist1)/(Jdist1-Jdist2))

         end if

      end do
      ! delta f
      eon%feon = 0.0_fp
      return
   end subroutine set_feq


   subroutine rk23_adstep(fstp, dtstp, dH_q, dH_p, Sr_eon, err)
      ! update the distribution at each step
      ! fstp(0:Nq,0:Np,1)   : f(qi,pj;t)
      ! fstp(0:Nq-1,0:Np,2) : \int_{qi}^{qi+dq} f(q,pj;t) dq
      !                       fstp(Nq,:,2) has no definition
      ! fstp(0:Nq,0:Np-1,3) : \int_{pj}^{pj+dp} f(qi,p;t) dp
      !                       fstp(:,Np,3) has no definition
      ! fstp(0:Nq-1,0:Np-1,4) : \int_{qi}^{qi+dq} \int_{pj}^{pj+dp} f(q,p;t) dqdp
      !                         fstp(Nq,:,4) and fstp(:,Np,4) has no definition
      implicit none
      real(fp), intent(inout) :: fstp(0:Nq,0:Np,4)
      real(fp), intent(in) :: dtstp
      real(fp), intent(in) :: dH_q(0:Nq,0:Np), dH_p(0:Nq,0:Np)
      real(fp), intent(in) :: Sr_eon(0:Nq,0:Np,4)
      real(fp), intent(out) :: err
      real(fp) :: frk(0:Nq,0:Np,4), dfrk(0:Nq,0:Np,4), dfodt(0:Nq,0:Np,4)
      real(fp) :: collision(0:Nq,0:Np,4)
      real(fp) :: rk(0:Nq,0:Np,4), enorm(0:Nq,0:Np,4)
      integer  :: i, j
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Runge-Kutta substep 1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! dfrk(0:Nq,0:Np,1) : df/dq(qi,pj)
      ! dfrk(0:Nq,0:Np,2) : df/dp(qi,pj)
      ! dfrk(0:Nq,0:Np-1,3) : d/dq \int_{pj}^{pj+dp} f(qi,p) dp
      ! dfrk(0:Nq-1,0:Np,4) : d/dp \int_{qi}^{qi+dq} f(q,pj) dq
      frk = fstp
      dfrk = 0.0_fp
      do j = 0, Np
         do i = 1, Nq-1
            ! calculate the spatial derivative of the distribution
            ! boundary value : dfrk(0,0:Np,1), dfrk(Nq,0:Np,1)
            dfrk(i,j,1) = upwind_deriv(frk(i-1,j,1), frk(i,j,1), frk(i+1,j,1),       &
               frk(i-1,j,2), frk(i,j,2), dH_p(i,j), dq)
         end do
      end do
      do j = 1, Np-1
         do i = 0, Nq
            ! calculate the velocity derivative of the distribution
            ! boundary value :: dfrk(0:Nq,0,2), dfrk(0:Nq,Np,2)
            dfrk(i,j,2) = upwind_deriv(frk(i,j-1,1), frk(i,j,1), frk(i,j+1,1),       &
               frk(i,j-1,3), frk(i,j,3), -dH_q(i,j), dp)
         end do
      end do
      do j = 0, Np-1
         do i = 1, Nq-1
            ! calculate the spatial derivatives of the line integrated distribution
! boundary value :: dfrk(0,0:Np-1,3), dfrk(Nq,0:Np-1,3)
            dfrk(i,j,3) = upwind_deriv(frk(i-1,j,3), frk(i,j,3), frk(i+1,j,3),       &
               frk(i-1,j,4), frk(i,j,4), dH_p(i,j)+dH_p(i,j+1), dq)
         end do
      end do
      do j = 1, Np-1
         do i = 0, Nq-1
            ! calculate the velocity derivative of the line integrated distribution
! boundary value :: dfrk(0:Nq-1,0,4), dfrk(0:Nq-1,Np,4)
            dfrk(i,j,4) = upwind_deriv(frk(i,j-1,2), frk(i,j,2), frk(i,j+1,2),       &
               frk(i,j-1,4), frk(i,j,4), -dH_q(i,j)-dH_q(i+1,j), dp)
         end do
      end do
      ! boundary condition
      call set_boundary(frk, dfrk, dH_q, dH_p)
      call collide_eon(collision, frk)
      ! calculate the time derivative of the distribution
      ! dfodt(0:Nq,0:Np,1)
      ! dfodt(0:Nq-1,0:Np,2)
      ! dfodt(0:Nq,0:Np-1,3)
! dfodt(0:Nq-1,0:Np-1,4)
      dfodt = 0.0_fp
      do j = 0, Np
         do i = 0, Nq
            dfodt(i,j,1) = -dH_p(i,j)*dfrk(i,j,1)+dH_q(i,j)*dfrk(i,j,2)             &
               + Sr_eon(i,j,1)+collision(i,j,1)
         end do
      end do
      do j = 0, Np
         do i = 0, Nq-1
            dfodt(i,j,2) = -line_intg(dH_p(i,j), 0.5_fp*(dH_p(i,j)+dH_p(i+1,j)), dH_p(i+1,j),     &
               dfrk(i,j,1), frk(i+1,j,1)-frk(i,j,1), dfrk(i+1,j,1), dq)    &
               + line_intg(dH_q(i,j), 0.5_fp*(dH_q(i,j)+dH_q(i+1,j)), dH_q(i+1,j),      &
               dfrk(i,j,2), dfrk(i,j,4), dfrk(i+1,j,2), dq)                 &
               + Sr_eon(i,j,2)+collision(i,j,2)
         end do
      end do
      do j = 0, Np-1
         do i = 0, Nq
            dfodt(i,j,3) = -line_intg(dH_p(i,j), 0.5_fp*(dH_p(i,j)+dH_p(i,j+1)), dH_p(i,j+1),     &
               dfrk(i,j,1), dfrk(i,j,3), dfrk(i,j+1,1), dp)                &
               + line_intg(dH_q(i,j), 0.5_fp*(dH_q(i,j)+dH_q(i,j+1)), dH_q(i,j+1),      &
               dfrk(i,j,2), frk(i,j+1,1)-frk(i,j,1), dfrk(i,j+1,2), dp)     &
               + Sr_eon(i,j,3)+collision(i,j,3)
         end do
      end do
      do j = 0, Np-1
         do i = 0, Nq-1
            dfodt(i,j,4) = -line_intg(dH_p(i+1,j), 0.5_fp*(dH_p(i+1,j)+dH_p(i+1,j+1)),            &
               dH_p(i+1,j+1),frk(i+1,j,1),frk(i+1,j,3),frk(i+1,j+1,1),dp)  &
               + line_intg(dH_p(i,j), 0.5_fp*(dH_p(i,j)+dH_p(i,j+1)), dH_p(i,j+1),      &
               frk(i,j,1), frk(i,j,3), frk(i,j+1,1), dp)                    &
               + line_intg(dH_q(i,j+1), 0.5_fp*(dH_q(i,j+1)+dH_q(i+1,j+1)),             &
               dH_q(i+1,j+1),frk(i,j+1,1),frk(i,j+1,2),frk(i+1,j+1,1),dq)   &
               - line_intg(dH_q(i,j), 0.5_fp*(dH_q(i,j)+dH_q(i+1,j)), dH_q(i+1,j),      &
               frk(i,j,1), frk(i,j,2), frk(i+1,j,1), dq)                    &
               + Sr_eon(i,j,4)+collision(i,j,4)
         end do
      end do
      rk = dfodt
      frk = fstp + rk * dtstp
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Runge-Kutta substep 2 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      dfrk = 0.0_fp
      do j = 0, Np
         do i = 1, Nq-1
            dfrk(i,j,1) = upwind_deriv(frk(i-1,j,1), frk(i,j,1), frk(i+1,j,1),      &
               frk(i-1,j,2), frk(i,j,2), dH_p(i,j), dq)
         end do
      end do
      do j = 1, Np-1
         do i = 0, Nq
            dfrk(i,j,2) = upwind_deriv(frk(i,j-1,1), frk(i,j,1), frk(i,j+1,1),      &
               frk(i,j-1,3), frk(i,j,3), -dH_q(i,j), dp)
         end do
      end do
      do j = 0, Np-1
         do i = 1, Nq-1
            dfrk(i,j,3) = upwind_deriv(frk(i-1,j,3), frk(i,j,3), frk(i+1,j,3),      &
               frk(i-1,j,4), frk(i,j,4), dH_p(i,j)+dH_p(i,j+1), dq)
         end do
      end do
      do j = 1, Np-1
         do i = 0, Nq-1
            dfrk(i,j,4) = upwind_deriv(frk(i,j-1,2), frk(i,j,2), frk(i,j+1,2),      &
               frk(i,j-1,4), frk(i,j,4), -dH_q(i,j)-dH_q(i+1,j), dp)
         end do
      end do
      call set_boundary(frk, dfrk, dH_q, dH_p)
      call collide_eon(collision, frk)
      dfodt = 0.0_fp
      do j = 0, Np
         do i = 0, Nq
            dfodt(i,j,1) = -dH_p(i,j)*dfrk(i,j,1)+dH_q(i,j)*dfrk(i,j,2)             &
               + Sr_eon(i,j,1)+collision(i,j,1)
         end do
      end do
      do j = 0, Np
         do i = 0, Nq-1
            dfodt(i,j,2) = -line_intg(dH_p(i,j), 0.5_fp*(dH_p(i,j)+dH_p(i+1,j)), dH_p(i+1,j),     &
               dfrk(i,j,1), frk(i+1,j,1)-frk(i,j,1), dfrk(i+1,j,1), dq)    &
               + line_intg(dH_q(i,j), 0.5_fp*(dH_q(i,j)+dH_q(i+1,j)), dH_q(i+1,j),      &
               dfrk(i,j,2), dfrk(i,j,4), dfrk(i+1,j,2), dq)                 &
               + Sr_eon(i,j,2)+collision(i,j,2)
         end do
      end do
      do j = 0, Np-1
         do i = 0, Nq
            dfodt(i,j,3) = -line_intg(dH_p(i,j), 0.5_fp*(dH_p(i,j)+dH_p(i,j+1)), dH_p(i,j+1),     &
               dfrk(i,j,1), dfrk(i,j,3), dfrk(i,j+1,1), dp)                &
               + line_intg(dH_q(i,j), 0.5_fp*(dH_q(i,j)+dH_q(i,j+1)), dH_q(i,j+1),      &
               dfrk(i,j,2), frk(i,j+1,1)-frk(i,j,1), dfrk(i,j+1,2), dp)     &
               + Sr_eon(i,j,3)+collision(i,j,3)
         end do
      end do
      do j = 0, Np-1
         do i = 0, Nq-1
            dfodt(i,j,4) = -line_intg(dH_p(i+1,j), 0.5_fp*(dH_p(i+1,j)+dH_p(i+1,j+1)),            &
               dH_p(i+1,j+1),frk(i+1,j,1),frk(i+1,j,3),frk(i+1,j+1,1),dp)  &
               + line_intg(dH_p(i,j), 0.5_fp*(dH_p(i,j)+dH_p(i,j+1)), dH_p(i,j+1),      &
               frk(i,j,1), frk(i,j,3), frk(i,j+1,1), dp)                    &
               + line_intg(dH_q(i,j+1), 0.5_fp*(dH_q(i,j+1)+dH_q(i+1,j+1)),             &
               dH_q(i+1,j+1),frk(i,j+1,1),frk(i,j+1,2),frk(i+1,j+1,1),dq)   &
               - line_intg(dH_q(i,j), 0.5_fp*(dH_q(i,j)+dH_q(i+1,j)), dH_q(i+1,j),      &
               frk(i,j,1), frk(i,j,2), frk(i+1,j,1), dq)                    &
               + Sr_eon(i,j,4)+collision(i,j,4)
         end do
      end do
      rk = rk + dfodt
      frk = fstp + 0.25_fp*rk * dtstp
      enorm = rk
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Runge-Kutta substep 3 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      dfrk = 0.0_fp
      do j = 0, Np
         do i = 1, Nq-1
            dfrk(i,j,1) = upwind_deriv(frk(i-1,j,1), frk(i,j,1), frk(i+1,j,1),      &
               frk(i-1,j,2), frk(i,j,2), dH_p(i,j), dq)
         end do
      end do
      do j = 1, Np-1
         do i = 0, Nq
            dfrk(i,j,2) = upwind_deriv(frk(i,j-1,1), frk(i,j,1), frk(i,j+1,1),      &
               frk(i,j-1,3), frk(i,j,3), -dH_q(i,j), dp)
         end do
      end do
      do j = 0, Np-1
         do i = 1, Nq-1
            dfrk(i,j,3) = upwind_deriv(frk(i-1,j,3), frk(i,j,3), frk(i+1,j,3),      &
               frk(i-1,j,4), frk(i,j,4), dH_p(i,j)+dH_p(i,j+1), dq)
         end do
      end do
      do j = 1, Np-1
         do i = 0, Nq-1
            dfrk(i,j,4) = upwind_deriv(frk(i,j-1,2), frk(i,j,2), frk(i,j+1,2),      &
               frk(i,j-1,4), frk(i,j,4), -dH_q(i,j)-dH_q(i+1,j), dp)
         end do
      end do
      call set_boundary(frk, dfrk, dH_q, dH_p)
      call collide_eon(collision, frk)
      dfodt = 0.0_fp
      do j = 0, Np
         do i = 0, Nq
            dfodt(i,j,1) = -dH_p(i,j)*dfrk(i,j,1)+dH_q(i,j)*dfrk(i,j,2)             &
               + Sr_eon(i,j,1)+collision(i,j,1)
         end do
      end do
      do j = 0, Np
         do i = 0, Nq-1
            dfodt(i,j,2) = -line_intg(dH_p(i,j), 0.5_fp*(dH_p(i,j)+dH_p(i+1,j)), dH_p(i+1,j),     &
               dfrk(i,j,1), frk(i+1,j,1)-frk(i,j,1), dfrk(i+1,j,1), dq)    &
               + line_intg(dH_q(i,j), 0.5_fp*(dH_q(i,j)+dH_q(i+1,j)), dH_q(i+1,j),      &
               dfrk(i,j,2), dfrk(i,j,4), dfrk(i+1,j,2), dq)                 &
               + Sr_eon(i,j,2)+collision(i,j,2)
         end do
      end do
      do j = 0, Np-1
         do i = 0, Nq
            dfodt(i,j,3) = -line_intg(dH_p(i,j), 0.5_fp*(dH_p(i,j)+dH_p(i,j+1)), dH_p(i,j+1),     &
               dfrk(i,j,1), dfrk(i,j,3), dfrk(i,j+1,1), dp)                &
               + line_intg(dH_q(i,j), 0.5_fp*(dH_q(i,j)+dH_q(i,j+1)), dH_q(i,j+1),      &
               dfrk(i,j,2), frk(i,j+1,1)-frk(i,j,1), dfrk(i,j+1,2), dp)     &
               + Sr_eon(i,j,3)+collision(i,j,3)
         end do
      end do
      do j = 0, Np-1
         do i = 0, Nq-1
            dfodt(i,j,4) = -line_intg(dH_p(i+1,j), 0.5_fp*(dH_p(i+1,j)+dH_p(i+1,j+1)),            &
               dH_p(i+1,j+1),frk(i+1,j,1),frk(i+1,j,3),frk(i+1,j+1,1),dp)  &
               + line_intg(dH_p(i,j), 0.5_fp*(dH_p(i,j)+dH_p(i,j+1)), dH_p(i,j+1),      &
               frk(i,j,1), frk(i,j,3), frk(i,j+1,1), dp)                    &
               + line_intg(dH_q(i,j+1), 0.5_fp*(dH_q(i,j+1)+dH_q(i+1,j+1)),             &
               dH_q(i+1,j+1),frk(i,j+1,1),frk(i,j+1,2),frk(i+1,j+1,1),dq)   &
               - line_intg(dH_q(i,j), 0.5_fp*(dH_q(i,j)+dH_q(i+1,j)), dH_q(i+1,j),      &
               frk(i,j,1), frk(i,j,2), frk(i+1,j,1), dq)                    &
               + Sr_eon(i,j,4)+collision(i,j,4)
         end do
      end do
      rk = rk + 4.0_fp * dfodt
      frk = fstp + rk * dtstp/6.0_fp
      enorm = abs(1.5_fp*enorm-0.5_fp*rk) * dtstp/3.0_fp
      dfrk = abs(fstp+rk*dtstp/6.0_fp)+TINY_NUM
      err = maxval(abs(enorm/dfrk))
      fstp = frk
      return
   end subroutine rk23_adstep


   subroutine set_boundary(fstp, dfstp, dH_q, dH_p)
      ! fstp(0:Nq,0:Np,1)   : f(qi,pj;t)
      ! fstp(0:Nq-1,0:Np,2) : \int_{qi}^{qi+dq} f(q,pj;t) dq
      !                       fstp(Nq,:,2) has no definition
      ! fstp(0:Nq,0:Np-1,3) : \int_{pj}^{pj+dp} f(qi,p;t) dp
      !                       fstp(:,Np,3) has no definition
      ! fstp(0:Nq-1,0:Np-1,4) : \int_{qi}^{qi+dq} \int_{pj}^{pj+dp} f(q,p;t) dqdp
      !                         fstp(Nq,:,4) and fstp(:,Np,4) has no definition
      ! dfstp(0:Nq,0:Np,1) : df/dq(qi,pj)
      ! dfstp(0:Nq,0:Np,2) : df/dp(qi,pj)
      ! dfstp(0:Nq,0:Np-1,3) : d/dq \int_{pj}^{pj+dp} f(qi,p) dp
      ! dfstp(0:Nq-1,0:Np,4) : d/dp \int_{qi}^{qi+dq} f(q,pj) dq
      !
      ! set the boundary conditions for,
      ! dfstp(0,0:Np,1),   dfstp(Nq,0:Np,1)
      ! dfstp(0:Nq,0,2),   dfstp(0:Nq,Np,2)
      ! dfstp(0,0:Np-1,3), dfstp(Nq,0:Np-1,3)
! dfstp(0:Nq-1,0,4), dfstp(0:Nq-1,Np,4)
      implicit none
      real(fp), intent(in) :: fstp(0:Nq,0:Np,4)
      real(fp), intent(out) :: dfstp(0:Nq,0:Np,4)
      real(fp), intent(in) :: dH_q(0:Nq,0:Np), dH_p(0:Nq,0:Np)
      real(fp) :: lft(0:Np,4), rht(0:Np,4)
      real(fp) :: bot(0:Nq,4), top(0:Nq,4)
      integer :: i, j
      lft = 0.0_fp
      rht = 0.0_fp
      ! periodic boundary condtions on both qcor ends
      forall(j = 0:Np)   lft(j,1) = fstp(Nq-1,j,1)
      forall(j = 0:Np)   lft(j,2) = fstp(Nq-1,j,2)
      forall(j = 0:Np-1) lft(j,3) = fstp(Nq-1,j,3)
      forall(j = 0:Np-1) lft(j,4) = fstp(Nq-1,j,4)
      forall(j = 0:Np)   rht(j,1) = fstp(1,j,1)
      forall(j = 0:Np)   rht(j,2) = fstp(0,j,2)
      forall(j = 0:Np-1) rht(j,3) = fstp(1,j,3)
      forall(j = 0:Np-1) rht(j,4) = fstp(0,j,4)
      forall(j = 0:Np) dfstp(0,j,1) = upwind_deriv(lft(j,1), fstp(0,j,1), fstp(1,j,1),         &
         lft(j,2), fstp(0,j,2), dH_p(0,j), dq)
      forall(j = 0:Np) dfstp(Nq,j,1) = upwind_deriv(fstp(Nq-1,j,1), fstp(Nq,j,1), rht(j,1),    &
         fstp(Nq-1,j,2), rht(j,2), dH_p(Nq,j), dq)
      forall(j = 0:Np-1) dfstp(0,j,3) = upwind_deriv(lft(j,3), fstp(0,j,3), fstp(1,j,3),       &
         lft(j,4), fstp(0,j,4),                    &
         dH_p(0,j)+dH_p(0,j+1), dq)
      forall(j = 0:Np-1) dfstp(Nq,j,3) = upwind_deriv(fstp(Nq-1,j,3), fstp(Nq,j,3), rht(j,3),  &
         fstp(Nq-1,j,4), rht(j,4),                &
         dH_p(Nq,j)+dH_p(Nq,j+1), dq)
      bot = 0.0_fp
      top = 0.0_fp
      forall(i = 0:Nq)  dfstp(i,0,2) = upwind_deriv(bot(i,1), fstp(i,0,1), fstp(i,1,1),        &
         bot(i,3), fstp(i,0,3), -dH_q(i,0), dp)
      forall(i = 0:Nq)  dfstp(i,Np,2) = upwind_deriv(fstp(i,Np-1,1), fstp(i,Np,1), top(i,1),   &
         fstp(i,Np-1,3), top(i,3), -dH_q(i,Np), dp)
      forall(i = 0:Nq-1)  dfstp(i,0,4) = upwind_deriv(bot(i,2), fstp(i,0,2), fstp(i,1,2),      &
         bot(i,4), fstp(i,0,4), -dH_q(i,0)-dH_q(i+1,0), dp)
      forall(i = 0:Nq-1)  dfstp(i,Np,4) = upwind_deriv(fstp(i,Np-1,2), fstp(i,Np,2), top(i,2), &
         fstp(i,Np-1,4), top(i,4), -dH_q(i,Np)-dH_q(i+1,Np), dp)
      return
   end subroutine set_boundary


   subroutine get_Hamiltonian(eon)
      ! Hamiltonian for one macro-particle
      implicit none
      type(energetic_eon), intent(inout) :: eon
      integer :: i, j
      real(fp) :: mu
      do j = 0, Np
         do i = 0, Nq
            mu = max(eon%Jpos+eon%pcr+eon%pcor(j), small)
            eon%H0(i,j) = 0.5_fp*(eon%kmode*eon%pcor(j))**2 + sqrt(2.0_fp*eon%gyro*mu) * real(eon%ampl*exp(-mathj*eon%qcor(i)))
            eon%dH0_q(i,j) = sqrt(2.0_fp*eon%gyro*mu) * aimag(eon%ampl*exp(-mathj*eon%qcor(i)))
            eon%dH_p(i,j) = eon%kmode*eon%kmode*eon%pcor(j) + eon%gyro/sqrt(2.0_fp*eon%gyro*mu) * real(eon%ampl*exp(-mathj*eon%qcor(i)))
            eon%dH_q(i,j) = eon%dH0_q(i,j) + eon%drg
         end do
      end do
      return
   end subroutine get_Hamiltonian


   subroutine collide_eon(collision, fstp)
      implicit none
      real(fp), intent(out) :: collision(0:Nq,0:Np,4)
      real(fp), intent(in) :: fstp(0:Nq,0:Np,4)
      real(fp) :: d1f(0:Nq,0:Np), d2f(0:Nq,0:Np)
      real(fp) :: bot(0:Nq,4), top(0:Nq,4)
      integer :: i, j
      ! no collision
      if (colla == 0.0_fp .and. colld3 == 0.0_fp)  then
         collision = 0.0_fp
      else
         bot = 0.0_fp
         top = 0.0_fp
         do j = 1, Np-1
            do i = 0, Nq
               call central_diff(d1f(i,j), d2f(i,j), fstp(i,j-1,1), fstp(i,j,1),           &
                  fstp(i,j+1,1), fstp(i,j-1,3), fstp(i,j,3), dp)
               collision(i,j,1) = -colla*fstp(i,j,1)+colld3*d2f(i,j)
            end do
         end do
         do i = 0, Nq
            call central_diff(d1f(i,0), d2f(i,0), bot(i,1), fstp(i,0,1),                   &
               fstp(i,1,1), bot(i,3), fstp(i,0,3), dp)
            collision(i,0,1) = -colla*fstp(i,0,1)+colld3*d2f(i,0)
         end do
         do i = 0, Nq
            call central_diff(d1f(i,Np), d2f(i,Np), fstp(i,Np-1,1), fstp(i,Np,1),          &
               top(i,1), fstp(i,Np-1,3), top(i,3), dp)
            collision(i,Np,1) = -colla*fstp(i,Np,1)+colld3*d2f(i,Np)
         end do
         do j = 0, Np-1
            do i = 0, Nq
               collision(i,j,3) = -colla*fstp(i,j,3)+colld3*(d1f(i,j+1)-d1f(i,j))
            end do
         end do
         do j = 1, Np-1
            do i = 0, Nq-1
               call central_diff(d1f(i,j), d2f(i,j), fstp(i,j-1,2), fstp(i,j,2),           &
                  fstp(i,j+1,2), fstp(i,j-1,4), fstp(i,j,4), dp)
               collision(i,j,2) = -colla*fstp(i,j,2)+colld3*d2f(i,j)
            end do
         end do
         do i = 0, Nq-1
            call central_diff(d1f(i,0), d2f(i,0), bot(i,2), fstp(i,0,2),                    &
               fstp(i,1,2), bot(i,4), fstp(i,0,4), dp)
            collision(i,0,2) = -colla*fstp(i,0,2)+colld3*d2f(i,0)
         end do
         do i = 0, Nq-1
            call central_diff(d1f(i,Np), d2f(i,Np), fstp(i,Np-1,2), fstp(i,Np,2),           &
               top(i,2), fstp(i,Np-1,4), top(i,4), dp)
            collision(i,Np,2) = -colla*fstp(i,Np,2)+colld3*d2f(i,Np)
         end do
         do j = 0, Np-1
            do i = 0, Nq-1
               collision(i,j,4) = -colla*fstp(i,j,4)+colld3*(d1f(i,j+1)-d1f(i,j))
            end do
         end do
      end if
      return
   end subroutine collide_eon

   subroutine build_eon(eon)
      ! equilibrium source term for one macro-particle
      implicit none
      type(energetic_eon), intent(inout) :: eon
      integer :: i, j
      eon%Sr = 0.0_fp
      do j = 0, Np
         do i = 0, Nq
            eon%Sr(i,j,1) = eon%dfeq(j)*eon%dH0_q(i,j)
         end do
      end do
      do j = 0, Np
         do i = 0, Nq-1
            eon%Sr(i,j,2) = eon%dfeq(j)*(eon%H0(i+1,j)-eon%H0(i,j))
         end do
      end do
      do j = 0, Np-1
         do i = 0, Nq
            eon%Sr(i,j,3) = line_intg(eon%dH0_q(i,j), 0.5_fp*(eon%dH0_q(i,j)             &
               + eon%dH0_q(i,j+1)), eon%dH0_q(i,j+1), eon%dfeq(j),            &
               eon%feq(j+1)-eon%feq(j), eon%dfeq(j+1), dp)
         end do
      end do
      do j = 0, Np-1
         do i = 0, Nq-1
            eon%Sr(i,j,4) = line_intg(eon%H0(i+1,j)-eon%H0(i,j), 0.5_fp*(eon%H0(i+1,j)   &
               - eon%H0(i,j)+eon%H0(i+1,j+1)-eon%H0(i,j+1)),                  &
               eon%H0(i+1,j+1)-eon%H0(i,j+1), eon%dfeq(j),                  &
               eon%feq(j+1)-eon%feq(j), eon%dfeq(j+1), dp)
         end do
      end do
      return
   end subroutine build_eon

   !eons push eon along z and interp feon to adject grids
   subroutine shift_eon_mov(eon)
      implicit none
      type(energetic_eon), intent(inout) :: eon(-NZ:NZ,NJ)
      real(fp) :: omega, gyro, wp2, kmode, pcr, zpos
      real(fp) :: rk(4)
      real(fp) :: newz,newj
      integer :: k, l
      do l = 1, NJ
         do k = -NZ, NZ
            zpos=eon(k,l)%zpos
            omega=eon(k,l)%omega
            ! update zpos
            ! RK1
            gyro = gyro0*(1.0_fp+4.5_fp*(zpos/(Lshell*RE))**2)
            wp2 = 1.0_fp+cold2*(zpos/(Lshell*RE))**2
            kmode = sqrt(omega**2+wp2*omega/(gyro-omega))
            rk(1) = (omega-gyro)/kmode
            newz = eon(k,l)%zpos + 0.5_fp*dT*rk(1)
            ! RK2
            gyro = gyro0*(1.0_fp+4.5_fp*(newz/(Lshell*RE))**2)
            wp2 = 1.0_fp+cold2*(newz/(Lshell*RE))**2
            kmode = sqrt(omega**2+wp2*omega/(gyro-omega))
            rk(2) = (omega-gyro)/kmode
            newz = eon(k,l)%zpos + 0.5_fp*dT*rk(2)
            ! RK3
            gyro = gyro0*(1.0_fp+4.5_fp*(newz/(Lshell*RE))**2)
            wp2 = 1.0_fp+cold2*(newz/(Lshell*RE))**2
            kmode = sqrt(omega**2+wp2*omega/(gyro-omega))
            rk(3) = (omega-gyro)/kmode
            newz = eon(k,l)%zpos + dT*rk(3)
            ! RK4
            gyro = gyro0*(1.0_fp+4.5_fp*(newz/(Lshell*RE))**2)
            wp2 = 1.0_fp+cold2*(newz/(Lshell*RE))**2
            kmode = sqrt(omega**2+wp2*omega/(gyro-omega))
            rk(4) = (omega-gyro)/kmode
            newz = eon(k,l)%zpos + dT*(rk(1)+2.0_fp*(rk(2)+rk(3))+rk(4))/6.0_fp
            !update Jpos and parameters
            gyro = gyro0*(1.0_fp+4.5_fp*(newz/(Lshell*RE))**2)
            wp2 = 1.0_fp+cold2*(newz/(Lshell*RE))**2
            kmode = sqrt(omega**2+wp2*omega/(gyro-omega))
            pcr = (omega-gyro)/kmode**2
            newj = (eon(k,l)%const-0.5_fp*(omega-gyro)**2/kmode**2)/(omega-gyro)
            eon(k,l)%zpos = newz
            eon(k,l)%Jpos = newj
            eon(k,l)%gyro = gyro
            eon(k,l)%wp2 = wp2
            eon(k,l)%kmode = kmode
            eon(k,l)%pcr = pcr
         end do
         eon(NZ,l)%feon = 0.0_fp
      end do
      return
   end subroutine shift_eon_mov


   subroutine intgl_source_mov(eon)
      ! chorus%Sr = kmode*nh*\iint dqcor*dpcor*dJact
      !           * sqrt(2*gyro*(Jact+pcr+pcor))*feon*exp(mathj*qcor)
      implicit none
      type(energetic_eon), intent(inout) :: eon(-NZ:NZ,NJ)
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
      do l = 1, NJ
         do k = -NZ+1, NZ
            eon(k,l)%intf = 0.5_fp*nh*(Sr(k)+Sr(k-1))
         end do
      end do
      return
   end subroutine intgl_source_mov

   subroutine cubic_abcd(a,b,c,d,df0,f0,f1,df1)
      implicit none
      complex(fp), intent(in) :: df0,f0,f1,df1
      complex(fp), intent(out) :: a,b,c,d
      d = f0
      c = df0
      b = -3*d/dZ**2 - 2 * c/dZ + 3*f1/dZ**2 - df1/dZ
      a = df1/dZ**2 - 2*f1/dZ**3 + c/dZ**2 + 2*d/dZ**3
   end subroutine
   subroutine icubic(y,x,xl,ps)
      complex(fp), intent(in) :: ps(4)
      real(fp), intent(in) :: x, xl
      complex(fp), intent(inout) :: y
      complex(fp) :: a,b,c,d
      call cubic_abcd(a,b,c,d,(ps(2)-ps(1))/dZ,ps(2),ps(3),(ps(4)-ps(3))/dZ)
      y = a*(x-xl)**3 + b*(x-xl)**2 + c*(x-xl) + d 
   end subroutine

   subroutine interp_ampl(eon,chorus,flag)
      implicit none
      type(energetic_eon), intent(out) :: eon(-NZ:NZ,NJ)
      type(chorus_mode), intent(in) :: chorus(-NZ:NZ+1)
      complex(fp) :: ampint(-NZ:NZ), ampgrid(-NZ:NZ+1)
      integer, intent(in) :: flag
      real(fp) :: zint(-NZ:NZ), zgrid(-NZ:NZ+1)
      integer :: k, l, i
      if (flag .eq. 0) then
      do l = 1, NJ
         do k = -NZ, NZ-1
            !locate eon
            i = INT(eon(k,l)%zpos/dZ)
            if (i>-NZ) then
               eon(k,l)%ampl = chorus(i)%ampl(1) + (chorus(i+1)%ampl(1) - chorus(i)%ampl(1))/(dZ)*(eon(k,l)%zpos-chorus(i)%zpos)
            else
               eon(k,l)%ampl = (0,0)
            end if
         end do
      end do
      else if (flag .eq. 1) then
         do l = 1, NJ
            do k = -NZ, NZ-1
               !locate eon
               i = INT(eon(k,l)%zpos/dZ)
               if (i<-NZ) then
                  eon(k,l)%ampl = (0,0)
               else if (i .eq. -NZ) then
                  call icubic( eon(k,l)%ampl, eon(k,l)%zpos,chorus(k)%zpos,(/(0.0_fp,0.0_fp), chorus(i)%ampl(1), chorus(i+1)%ampl(1) , chorus(i+2)%ampl(1) /))
               else
                  call icubic( eon(k,l)%ampl, eon(k,l)%zpos,chorus(k)%zpos,(/chorus(i-1)%ampl(1), chorus(i)%ampl(1), chorus(i+1)%ampl(1) , chorus(i+2)%ampl(1)/))
               end if
            end do
         end do
      else if (flag .eq. 2) then
         do k = -NZ, NZ+1
            ampgrid(k) = chorus(k)%ampl(1)
            zgrid(k) = chorus(k)%zpos
         end do
         do k = -NZ, NZ
            zint(k) = eon(k,1)%zpos
         end do
         call ispline(ampint, zint, 2*NZ+1, ampgrid, zgrid, 2*NZ+2, 1)
         do l = 1, NJ
            do k = -NZ, NZ-1
               eon(k,l)%ampl=ampint(k)
            end do
         end do
      end if
         return
      end subroutine interp_ampl

      subroutine interp_intf(eon,chorus)
         implicit none
         type(energetic_eon), intent(in) :: eon(-NZ:NZ,NJ)
         type(chorus_mode), intent(out) :: chorus(-NZ:NZ+1)
         complex(fp) :: intf
         integer :: k, i
         !reset current
         forall(k = -NZ+1:NZ)  chorus(k)%Sr = (0.0_fp, 0.0_fp)
         do k = -NZ, NZ-1
            i = INT(eon(k,1)%zpos/dZ)
            if (i>-NZ) then
               intf = eon(k,1)%intf
               chorus(i)%Sr = chorus(i)%Sr + (eon(k,1)%zpos-chorus(i)%zpos)/dZ*intf
               chorus(i+1)%Sr = chorus(i+1)%Sr - (eon(k,1)%zpos-chorus(i+1)%zpos)/dZ*intf
            end if
         end do
         chorus(-NZ)%Sr = (0.0_fp, 0.0_fp)
         chorus(NZ+1)%Sr = (0.0_fp, 0.0_fp)
         return
      end subroutine interp_intf
end module chorus_EON_MOV



 !!$program test_chorus_EON
 !!$  use chorus_IO
 !!$  use chorus_TYPE
 !!$  use chorus_LIB
 !!$  use chorus_EQ
 !!$  use chorus_WAVE
 !!$  use chorus_EON
 !!$  implicit none
 !!$  type(energetic_eon) :: eon(0:NZ)
 !!$  type(chorus_mode) :: chorus(0:NZ+1)
!!$  real(fp) :: dist(0:Nq,0:Np)
 !!$  integer :: tstep, i, j, k
 !!$  call init_chorus(chorus)
!!$  call init_eon(eon,chorus)
 !!$  call save_output('chorus.out', chorus%ampl(1), NZ+1, .FALSE.)
 !!$  do tstep = NT0+1, NT
 !!$     do k = 1, NZ
!!$        call trap_eon(eon(k), chorus(k))
 !!$        write(*,*) 'k : ', k
 !!$     end do
 !!$     call update_chorus(chorus, eon)
!!$     call shift_eon(eon)
 !!$     write(*,*) 'tstep : ', tstep,                                              &
 !!$          sqrt(sqrt(2.0_fp*chorus(NZ/2)%gyro*max(C0+chorus(NZ/2)%pres,0.0_fp))  &
 !!$        * abs(chorus(NZ/2)%ampl(1)))
 !!$     call save_output('chorus.out', chorus%ampl(1), NZ+1, .TRUE.)
 !!$     call save_output('dist.out', eon(0)%feon(0:Nq,0:Np,1), Nq+1, Np+1, .FALSE.)
 !!$     do k = 1, NZ-1
 !!$        call save_output('dist.out', eon(k)%feon(0:Nq,0:Np,1), Nq+1, Np+1, .TRUE.)
 !!$     end do
 !!$  end do
 !!$  stop
 !!$end program test_chorus_EON
