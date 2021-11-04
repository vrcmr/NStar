!= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = !
!                                                                              !
! Copyright (C) 2010-2017  Victor Mourão Roque <victor.raphael@gmail.com>      !
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
  use constants , only : dp
                         
  
  IMPLICIT NONE


!==============================================================================!  
! Description:                                                                 !
!   Newtonian Ideal equation of state used to complete the hydro equations     !
!                                                                              !
! Author: Victor Mourão Roque                                                  !
! last modification:                                                           !
! date: 08/07/2017 ; by Victor Mourao Roque                                    ! 
!                                                                              !
! Refs:                                                                        !
!                                                                              !
!                                                                              !
!==============================================================================!


  PRIVATE  ! Private scope by default


! Module variables
  real(dp), save, public :: gamma = 1.4d0


! Public members
  PUBLIC :: initialize_EoS, eps, pressure, sound, k_tilde, density


CONTAINS
  
!==============================================================================!  
! Description:                                                                 !
!   initialize the parameters to EoS, in this case the gamma                   !
!  Talvez fazer igual ao PLUTO que inicia gamma na condição inicial
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

  call get_parameter("gamma", gamma)


END SUBROUTINE initialize_EoS


!==============================================================================!
 real(dp) elemental FUNCTION eps(d,p)   result(internal_energy)
!==============================================================================!
! d = density,  p = pressure
!
 IMPLICIT NONE
 REAL(DP), INTENT(IN) :: d, p
 
 internal_energy = p / d / (gamma - 1.)  ! eq. 1.18 of Toro
 
END FUNCTION eps


!==============================================================================!
 real(dp) elemental FUNCTION pressure(d, epsilon)   result(prs)
!==============================================================================!
! d = density,  epsilon = internal energy
!
 IMPLICIT NONE
 REAL(DP),INTENT(IN) :: d, epsilon
 
  prs = epsilon * d * (gamma - 1.d0)    ! inverting eq. 1.18 of Toro

 END   FUNCTION pressure

!==============================================================================!
 real(dp) elemental FUNCTION sound( d, p)  result(cs) 
!==============================================================================!
! d = density,  p = pressure, sound = velocity of sound
!
 IMPLICIT NONE
 REAL(DP), INTENT(IN) :: d,p
 
  cs = dsqrt(dabs(gamma*p/d))  !Toro 1.35
 
 END FUNCTION sound

!==============================================================================!
 real(dp) elemental FUNCTION k_tilde(d, epsilon)  result(k)
!==============================================================================!
! d = density,  p = pressure, sound = velocity of sound
!
! k_tilde is necessary to solve the hydro equations with a arbitrary EoS. 
! In case of ideal gas, it takes this simple function. 
! k_tilde = \frac{\partial p}{\partial eps}
! Adicionar referência
!
 IMPLICIT NONE
 REAL(DP), INTENT(IN) :: d, epsilon
 
 
  k = gamma - 1.d0 

 END   FUNCTION k_tilde





!==============================================================================!
 real(dp) elemental FUNCTION density(prs)  result(rho)
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


 rho = prs**(1/gamma)


END FUNCTION density



END MODULE EoS
