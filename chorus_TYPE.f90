module chorus_TYPE
  use chorus_IO
  implicit none
  
  type :: energetic_eon
     ! resonance cell location along the slowly varying coordinate
     real(fp) :: zpos, Jact
     ! fast varying phase space mesh grid
     real(fp) :: qcor(0:Nq), pcor(0:Np)
     ! Hamiltonian associated functions and derivatives
     real(fp) :: dH_q(0:Nq,0:Np), dH_p(0:Nq,0:Np)
     real(fp) :: H0(0:Nq,0:Np), dH0_q(0:Nq,0:Np)
     ! delta f distribution
     real(fp) :: feon(0:Nq,0:Np,4)
     ! equilibrium distribution and its derivative w.r.t Omega
     real(fp) :: feq(0:Np), dfeq(0:Np)
     ! equilibrium source terms
     real(fp) :: Sr(0:Nq,0:Np,4)
     ! slowly varying parameters:
     ! gyro-frequency at resonance location: gyro
     ! colde plasma frequency at resonance location: wp2
     ! central mode number at resonance location: kmode
     ! cyclotron resonance action at resonance location: pcr=(omega-gyro)/kmode**2
     real(fp) :: gyro, wp2, kmode, pcr
     ! drag force at mode location: drg
     real(fp) :: drg 
     ! Jact at mode location: Jpos
     real(fp) :: Jpos
     ! does the cell satisfy the resonant condition?
     logical :: is_resonant
     !zjs added
     real(fp) :: const, omega
     complex(fp) :: ampl
     complex(fp) :: intf
  end type energetic_eon


  type :: chorus_mode
     ! wave location
     real(fp) :: zpos
     ! wave variables: 
     ! ampl(1)=a
     ! ampl(2)=da/dt
     ! ampl(3)=int_0^t dtau da/dtau*exp(-mathj*(omega-gyro)*(t-tau))
     complex(fp) :: ampl(3)
     ! wave source : Sr = kmode*nh*\iiint dqcor*dpcor*dJact 
     !                  * sqrt(2*gyro*(Jact+pcor+pcr))*feon*exp(mathj*qcor)
     complex(fp) :: Sr
     ! slowly varying parameters: 
     ! most unstable mode frequency: omega
     ! gyro-frequency at mode location: gyro
     ! cold plasma frequency at mode location: wp2
     ! central mode number at mode location: kmode
     ! cycltron resonance action at mode location: pcr=(omega-gyro)/kmode**2
     real(fp) :: omega, gyro, wp2, kmode, pcr
  end type chorus_mode
   
end module chorus_TYPE


!!$program test_chorus_TYPE
!!$  use chorus_IO
!!$  use chorus_TYPE
!!$  implicit none
!!$  type(energetic_eon) :: eon(10)
!!$  type(chorus_mode) :: chorus(10)
!!$  stop
!!$end program test_chorus_TYPE
