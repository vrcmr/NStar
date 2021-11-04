!= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = !
!                                                                              !
! Copyright (C) 2010-2021  Victor Mourão Roque <victor.raphael@gmail.com>      !
!                                                                              !
!    This program is free software: you can redistribute it and/or modify      !
!    it under the terms of the GNU General Public License as published by      !
!    the Free Software Foundation, either version 3 of the License, or         !
!    (at your option) any later version.                                       !
!                                                                              !
!    This program is distributed in the hope that it will be useful,           !
!    but WITHOUT ANY WARRANTY; without even the implied warranty of            !
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             !
!    GNU General Public License for more details.                              !
!                                                                              !
!    You should have received a copy of the GNU General Public License         !
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.     !
!= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = !


MODULE EoS
  use constants , only : dp, PI_D
                         
  
  IMPLICIT NONE


!==============================================================================!  
! Description:                                                                 !
!   EoS for combustion in a NS                                                 !
!                                                                              !
! Author: Victor Mourão Roque                                                  !
! last modification:                                                           !
! date: 18/08/2021 ; by Victor Mourao Roque                                    !
!                                                                              !
! Refs:  O'Boyle PHYSICAL REVIEW D 102, 083027 (2020)                          !
! the equations of O'Boyle uses g/cm^3: 1 MeV/fm3 = 1.7827d12 g/cm^3           !
!                                                                              !
!                                                                              !
!==============================================================================!


  PRIVATE  ! Private scope by default


! Module variables
  ! real(dp), save, public :: gamma = 1.4d0
  real(dp), save :: ptr = 9.d0
  real(dp), save :: pstar = 24.d0
  real(dp), save :: gamma_th_red = 2.d0
  real(dp), save :: bag = 120.d0
  real(dp), save :: rho_0, rho_1, rho_2
  real(dp), save :: Gamma_1, Gamma_2, Gamma_3
  real(dp), save :: p_0, p_1, p_2, epsilon_0
  real(dp), save :: K_1, K_2, K_3
  real(dp), save :: a_1, a_2, a_3
  real(dp), save :: Lambda_1, Lambda_2, Lambda_3

!!!! LOW-density EOS (SLy) from Table II of O'Boyle 2020:
  real(dp), parameter :: c2 = 8.9875517874E20 ! squared speed of light
  real(dp), parameter :: rho_low1=0.d0, rho_low2=6.285d5, rho_low3=1.826d8,   &
                       rho_low4=3.350d11, rho_low5=5.317d11
  real(dp), parameter :: K_low1 =5.214d-9, K_low2=5.726d-8, K_low3 =1.662d-6, &
                       K_low4= -7.957d29, K_low5 =1.746d-8
  real(dp), parameter :: Gamma_low1=1.611d0, Gamma_low2=1.440d0,              &
                       Gamma_low3=1.269d0, Gamma_low4= -1.841d0,              &
                       Gamma_low5= 1.382d0  ! Gamma_low5= 1.393d0
  real(dp), parameter :: Lambda_low1=0.d0, Lambda_low2=-1.354d0,              &
                       Lambda_low3=-6.025d3, Lambda_low4=1.193d9,             &
                       Lambda_low5=7.077d8
  real(dp), parameter :: a_low1 =0.d0, a_low2=-1.861d-5, a_low3 =-5.278d-4,   &
                       a_low4= 1.035d-2, a_low5 =8.208d-3
  real(dp), parameter :: p_low1 = 0.d0, p_low2 = 11.452, p_low3 = 4.461d4,    &
                         p_low4 = 7.109d8, p_low5 = 9.875d8


! Public members
  PUBLIC :: initialize_EoS, eps, pressure, sound, k_tilde, density, ptr
  PUBLIC :: eps_q, eps_h, prs_q, prs_h, pstar, k_tilde_p


CONTAINS

!==============================================================================!
! Description:                                                                 !
!                                                                              !
! Author: Victor Mourão Roque                                                  !
! last modification:                                                           !
! date: dd/mm/yyyy ; by Victor Mourao Roque                                    !
!                                                                              !
! Refs:                                                                        !
!                                                                              !
!                                                                              !
!==============================================================================!

SUBROUTINE initialize_EoS()
use parameters, only : get_parameter

  IMPLICIT NONE
  character(len = 32), save     :: EoStype  = "stiff"

  call get_parameter("EoStype", EoStype)

  if (EoStype == "stiff" ) then
      rho_0 = 7.979946872680d13; rho_1 = 2.818382d14; rho_2 = 3.801893d14;
      K_1 = 6.02559d-28
      Gamma_1 = 2.764d0 ; Gamma_2 = 8.5d0; Gamma_3 = 3.2d0 ;
      gamma_th_red = 2.d0 - 1.d0 !
      ! ptr = 9.d0
  else if (EoStype == "soft" ) then
      rho_0 = 7.979946872680d13; rho_1 = 2.818382d14; rho_2 = 3.801893d14;
      K_1 = 6.02559d-28
      Gamma_1 = 2.752d0; Gamma_2 = 4.5d0; Gamma_3 = 3.5d0 ;
      ! gamma_th = 2.d0
      ! ptr = 9.d0 ! é preciso calcular o valor correto
  else if (EoStype == "intermediate" ) then
      rho_0 = 7.979946872680d13; rho_1 = 2.818382d14; rho_2 = 3.801893d14;
      K_1 = 6.02559d-28
      Gamma_1 = 2.764d0; Gamma_2 = 8.5d0; Gamma_3 = 3.2d0 ;
      ! gamma_th = 2.d0
      ! ptr = 9.d0 ! é preciso calcular o valor correto
  else
      write(*,*) "problem with EoS type"
  end if

  call get_parameter("Ptransition", ptr)
  call get_parameter("Pstar", pstar)
  write(*,*) "The transition occurs when the pressure reaches", ptr
  write(*,*) "The transition to three flavor quarks EoS occurs when the pressure reaches", pstar
  call get_parameter("bag", bag)

  !!!! HIGH-density EOS:

  K_2 = K_1 * Gamma_1 / Gamma_2 * rho_1**(Gamma_1 - Gamma_2)
  K_3 = K_2 * Gamma_2 / Gamma_3 * rho_2**(Gamma_2 - Gamma_3)

  p_0 = K_low5 * rho_0**Gamma_low5 + Lambda_low5

  Lambda_1 =  p_0 - K_1 * rho_0**Gamma_1      ! Eq. (4.10) of O'Boyle 2020
  Lambda_2 = Lambda_1 + (1.d0- Gamma_1/Gamma_2) *K_1 *rho_1**Gamma_1
  Lambda_3 = Lambda_2 + (1.d0- Gamma_2/Gamma_3) *K_2 *rho_2**Gamma_2


  p_1 = K_1 * rho_1**Gamma_1 + Lambda_1
  p_2 = K_2 * rho_2**Gamma_2 + Lambda_2

  epsilon_0 = K_low5 / (Gamma_low5 - 1.d0)  * rho_0**Gamma_low5 + &
              ( 1.d0 + a_low5) * rho_0  - Lambda_low5

  a_1 = epsilon_0 / rho_0 - K_1 /(Gamma_1 -1.d0) * rho_0**(Gamma_1 - 1.d0) + &
        Lambda_1/rho_0 - 1.d0
  a_2 = a_1 + Gamma_1 * (Gamma_2 -Gamma_1) / (Gamma_2 -1.d0) / (Gamma_1 -1.d0) &
        * K_1 * rho_1**(Gamma_1 -1.d0)
  a_3 = a_2 + Gamma_2 * (Gamma_3 -Gamma_2) / (Gamma_3 -1.d0) / (Gamma_2 -1.d0) &
        * K_2 * rho_2**(Gamma_2 -1.d0)

END SUBROUTINE initialize_EoS

elemental subroutine EoSparameters(rho, a, Gamma, K , Lambda)
  real(dp), intent(in)  :: rho
  real(dp), intent(out) :: a, Gamma, K , Lambda


  !!!! LOW-density EOS:
  if (rho .le. rho_low2)  then
  ! write(*,*) 'EOS low1'
    a = a_low1
    Gamma = Gamma_low1
    K = K_low1
    Lambda = Lambda_low1
  elseif ((rho .ge. rho_low2) .and. (rho .le. rho_low3)) then
  ! write(*,*) 'EOS low2'
    a = a_low2
    Gamma = Gamma_low2
    K = K_low2
    Lambda = Lambda_low2
  elseif ((rho .ge. rho_low3) .and. (rho .le. rho_low4)) then
  ! write(*,*) 'EOS low3'
    a = a_low3
    Gamma = Gamma_low3
    K = K_low3
    Lambda = Lambda_low3
  elseif ((rho .ge. rho_low4) .and. (rho .le. rho_low5)) then
  ! write(*,*) 'EOS low4'
    a = a_low4
    Gamma = Gamma_low4
    K = K_low4
    Lambda = Lambda_low4
  elseif ((rho .ge. rho_low5) .and. (rho .le. rho_0)) then
  ! write(*,*) 'EOS low5'
    a = a_low5
    Gamma = Gamma_low5
    K = K_low5
    Lambda = Lambda_low5
  elseif ((rho .ge. rho_0) .and. (rho .le. rho_1)) then
  ! write(*,*) 'EOS 1'
    a = a_1
    Gamma = Gamma_1
    K = K_1
    Lambda = Lambda_1
  elseif ((rho .ge. rho_1) .and. (rho .le. rho_2)) then
  ! write(*,*) 'EOS 2'
    a = a_2
    Gamma = Gamma_2
    K = K_2
    Lambda = Lambda_2
  else      !!!!!(p.ge.p_2)
  ! write(*,*) 'EOS 3'
    a = a_3
    Gamma = Gamma_3
    K = K_3
    Lambda = Lambda_3
  endif
end subroutine

!==============================================================================!
 real(dp) elemental FUNCTION eps_q(d,p)   result(internal_energy)
!==============================================================================!
! d = density,  p = pressure
!
 IMPLICIT NONE
 REAL(DP), INTENT(IN) :: d, p
 real(dp) :: E

  E = 3 * p + 4 * bag

  internal_energy = E / d - 1.d0

END FUNCTION eps_q


!==============================================================================!
 real(dp) elemental FUNCTION eps_h(d,p)   result(internal_energy)
!==============================================================================!
! d = density,  p = pressure
!
 IMPLICIT NONE
 REAL(DP), INTENT(IN) :: d, p
 real(dp) :: a, Gamma, K , Lambda, rho
 real(dp) :: eps_cold, eps_th, p_cold, p_th, E_cold, E

  rho = 1.7827d12 * d

  call EoSparameters(rho, a, Gamma, K , Lambda)
  ! (Eq. (4.5) of O'Boyle 2020)
  E_cold = K * rho**Gamma  / (Gamma - 1.d0) + (1.d0 + a) * rho - Lambda
  p_cold = K * rho**Gamma + Lambda

  E_cold = E_cold / 1.7827d12            ! energy density in MeV fm^{-3}
  p_cold = p_cold / 1.7827d12            ! pressure in MeV fm^{-3}

  eps_cold = E_cold/d - 1.d0

  p_th = p - p_cold

  eps_th = p_th / gamma_th_red / d

  internal_energy =  eps_th + eps_cold

END FUNCTION eps_h


!==============================================================================!
 real(dp) elemental FUNCTION eps(d,p)   result(internal_energy)
!==============================================================================!
! d = density,  p = pressure
!
 IMPLICIT NONE
 REAL(DP), INTENT(IN) :: d, p

if (p .ge. pstar)  then

  internal_energy = eps_q(d,p)

else

  internal_energy = eps_h(d,p)

endif

END FUNCTION eps



!==============================================================================!
 real(dp) elemental FUNCTION prs_q(d, epsilon)   result(prs)
!==============================================================================!
! d = density,  epsilon = internal energy
!
 IMPLICIT NONE
 REAL(DP),INTENT(IN) :: d, epsilon
 real(dp) :: E

 E = d  * (1.d0 + epsilon)

 prs = (E - 4*bag )/ 3.d0

END FUNCTION prs_q


!==============================================================================!
 real(dp) elemental FUNCTION prs_h(d, epsilon)   result(prs)
!==============================================================================!
! d = density,  epsilon = internal energy
!
 IMPLICIT NONE
 REAL(DP),INTENT(IN) :: d, epsilon
 real(dp) :: a, Gamma, K , Lambda
 real(dp) :: pth, E, gq, gh, rho
 real(dp) :: p_cold, eps_cold, p_th, E_cold


 rho = 1.7827d12 * d
 call EoSparameters(rho, a, Gamma, K , Lambda)
 ! (Eq. (4.4) of O'Boyle 2020)
 p_cold = K * rho**Gamma + Lambda
 E_cold = K * rho**Gamma  / (Gamma - 1.d0) + (1.d0 + a) * rho - Lambda

 E_cold = E_cold / 1.7827d12              ! energy density in MeV fm^{-3}
 p_cold = p_cold / 1.7827d12            ! pressure in MeV fm^{-3}

 eps_cold = E_cold/d - 1.d0
 p_th = gamma_th_red * d * (epsilon - eps_cold)
 prs = p_th + p_cold

END FUNCTION prs_h


!==============================================================================!
 real(dp) elemental FUNCTION pressure(d, epsilon)   result(prs)
!==============================================================================!
! d = density,  epsilon = internal energy
!
 IMPLICIT NONE
 REAL(DP),INTENT(IN) :: d, epsilon
 real(dp) :: a, Gamma, K , Lambda
 real(dp) :: pq, ph, pth, E, gq, gh, rho
 real(dp) :: p_cold, eps_cold, p_th, E_cold


 pq = prs_q(d, epsilon)

 ph = prs_h(d, epsilon)

  if (pq .ge. pstar) then
    ! write(*,*) "quarks"
    prs = pq

  else
    ! write(*,*) "hadrons"
    prs = ph

  endif


END FUNCTION pressure

!==============================================================================!
 real(dp) elemental FUNCTION sound( d, p)  result(cs)
!==============================================================================!
! d = density,  p = pressure, sound = velocity of sound
!
 IMPLICIT NONE
 REAL(DP), INTENT(IN) :: d,p
 real(dp) :: a, Gamma, K , Lambda, E, cs2, rho
 real(dp) :: eps_th, eps_cold, p_cold, p_th, E_cold, E_cgs, p_cgs, E_th

if ( p .ge. pstar) then

  cs = 0.5773502692d0 ! sqrt(1/3)

else

  rho = 1.7827d12 * d
  p_cgs = 1.7827d12 * p

  call EoSparameters(rho, a, Gamma, K , Lambda)

  p_cold = K * rho**Gamma + Lambda
  E_cgs = d * (1.d0 + eps(d,p)) * 1.7827d12 ! Etot

  ! E_cold = K * rho**Gamma  / (Gamma - 1.d0) + (1.d0 + a) * rho - Lambda

  ! p_th = p_cgs - p_cold
  ! E_th = p_th/ gamma_th_red
  ! E_cgs = E_cold + E_th


  cs2 = Gamma * (p_cold - Lambda) + (p_cgs - p_cold) &
        + gamma_th_red * (p_cgs + p_cold)
  cs2 = cs2 / (E_cgs + p_cgs )

  cs = sqrt(cs2)

endif

 END FUNCTION sound

!==============================================================================!
 real(dp) elemental FUNCTION k_tilde(d, epsilon)  result(ktil)
!==============================================================================!
! d = density,  p = pressure, sound = velocity of sound
!
! k_tilde is necessary to solve the hydro equations with a arbitrary EoS.
! In case of ideal gas, it takes this simple function.
! k_tilde = 1/d * \frac{\partial p}{\partial eps}
! Adicionar referência
!
 IMPLICIT NONE
 REAL(DP), INTENT(IN) :: d, epsilon
 real(dp) :: a, Gamma, K , Lambda, prs

 prs = pressure(d, epsilon)

if (prs .ge. pstar) then

  ! quarks phase
  ktil =  0.333333d0 ! 1.d0 / 3.d0

else

  ! hadrons phase
  ktil = gamma_th_red

endif


 END   FUNCTION k_tilde


!==============================================================================!
 real(dp) elemental FUNCTION k_tilde_p(d, p)  result(ktil)
!==============================================================================!
! d = density,  p = pressure, sound = velocity of sound
!
! k_tilde is necessary to solve the hydro equations with a arbitrary EoS.
! In case of ideal gas, it takes this simple function.
! k_tilde = 1/d * \frac{\partial p}{\partial eps}
! Adicionar referência
!
 IMPLICIT NONE
 REAL(DP), INTENT(IN) :: d, p

if (p .ge. pstar) then

  ! quarks phase
  ktil =  0.333333d0 ! 1.d0 / 3.d0

else

  ! hadrons phase
  ktil = gamma_th_red

endif


 END   FUNCTION k_tilde_p



!==============================================================================!
! Additional routines used in initial condition
!==============================================================================!



!==============================================================================!
 real(dp) FUNCTION density(prs)  result(rho)
!==============================================================================!
! d = density,  p = pressure, sound = velocity of sound
!
!
 IMPLICIT NONE
 REAL(DP), INTENT(IN) :: prs
 real(dp) :: a, Gamma, K , Lambda, p
 real(dp) :: E, mu, b, c, T0

 if ( prs .ge. pstar) then

  T0 = 40.d0 ! MeV

  a = 0.75d0 / PI_D**2
  b = 1.5d0 * T0**2
  c = 0.53 * PI_D**2 * T0**4 - (prs + bag) * 197.33**3

  mu = sqrt(0.5d0 * (-b + sqrt (b**2 - 4.d0*a*c)) / a )


  rho = 939.d0 * (T0**2 * mu + mu**3 / PI_D**2) / 197.33**3

    ! rho = 939.d0 * ( 1.333d0 * (prs + bag) )**0.75d0 / (1.772 * 197.33**0.75)

 else

  p = 1.7827d12 * prs !* 0.95d0 ! 5 % da pressão é referente a energia térmica

  call EoSparameters_p(p, a, Gamma, K , Lambda)
  rho = ((p - Lambda)/K)**(1/Gamma)
  ! E = K * rho**Gamma  / (Gamma - 1.d0) + (1.d0 + a) * rho - Lambda

  ! E = E/1.7827d12              ! energy density in MeV fm^{-3}
  ! rho = E/(1.01d0)

  ! rho = ((p - Lambda)/K)**(1/Gamma)

  rho =  rho/1.7827d12              ! energy density in MeV fm^{-3}
  ! essa densidade não está bem calculada porque a pressão que entra é a
  ! pressão total e não somente a p_cold

endif

END FUNCTION density


!==============================================================================!
 real(dp) elemental FUNCTION energy(prs)  result(E)
!==============================================================================!
! d = density,  p = pressure, sound = velocity of sound
!
! k_tilde is necessary to solve the hydro equations with a arbitrary EoS.
! In case of ideal gas, it takes this simple function.
! k_tilde = 1/d * \frac{\partial p}{\partial eps}
! Adicionar referência
!
 IMPLICIT NONE
 REAL(DP), INTENT(IN) :: prs
 real(dp) :: a, Gamma, K , Lambda, p, rho

 if (prs .ge. pstar) then

  E = 3.d0 * prs + 4.d0 * bag

 else

  p = 1.7827d12 * prs * 0.95d0 ! 5 % da pressão é referente a energia térmica

  call EoSparameters_p(p, a, Gamma, K , Lambda)
  rho = ((p - Lambda)/K)**(1/Gamma)

  E = K*rho**Gamma  / (Gamma - 1.d0) + (1.d0 + a)*rho  - Lambda    ! (Eq. (4.5) of O'Boyle 2020)

  E =  E/1.7827d12              ! energy density in MeV fm^{-3}
  ! essa energia não está bem calculada porque a pressão que entra é a
  ! pressão total e não somente a p_cold. Adicionei uma porcentagem, 95%, para
  ! tentar compensar isso

 endif

END FUNCTION energy

 ELEMENTAL subroutine EoSparameters_p(p, a, Gamma, K , Lambda)
  real(dp), intent(in)  :: p
  real(dp), intent(out) :: a, Gamma, K , Lambda

  !!!! LOW-density EOS:
  if (p .le. p_low2)  then
  ! write(*,*) 'EOS low1'
    a = a_low1
    Gamma = Gamma_low1
    K = K_low1
    Lambda = Lambda_low1
  elseif ((p .ge. p_low2) .and. (p .le. p_low3)) then
  ! write(*,*) 'EOS low2'
    a = a_low2
    Gamma = Gamma_low2
    K = K_low2
    Lambda = Lambda_low2
  elseif ((p .ge. p_low3) .and. (p .le. p_low4)) then
  ! write(*,*) 'EOS low3'
    a = a_low3
    Gamma = Gamma_low3
    K = K_low3
    Lambda = Lambda_low3
  elseif ((p .ge. p_low4) .and. (p .le. p_low5)) then
  ! write(*,*) 'EOS low4'
    a = a_low4
    Gamma = Gamma_low4
    K = K_low4
    Lambda = Lambda_low4
  elseif ((p .ge. p_low5) .and. (p .le. p_0)) then
  ! write(*,*) 'EOS low5'
    a = a_low5
    Gamma = Gamma_low5
    K = K_low5
    Lambda = Lambda_low5
  elseif ((p .ge. p_0) .and. (p .le. p_1)) then
  ! write(*,*) 'EOS 1'
    a = a_1
    Gamma = Gamma_1
    K = K_1
    Lambda = Lambda_1
  elseif ((p .ge. p_1) .and. (p .le. p_2)) then
  ! write(*,*) 'EOS 2'
    a = a_2
    Gamma = Gamma_2
    K = K_2
    Lambda = Lambda_2
  else      !!!!!(p.ge.p_2)
  ! write(*,*) 'EOS 3'
    a = a_3
    Gamma = Gamma_3
    K = K_3
    Lambda = Lambda_3
  endif
end subroutine




END MODULE EoS




