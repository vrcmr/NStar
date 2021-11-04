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


MODULE init_cond
  USE constants,   only : i4b, dp, Pi_D
  USE EoS,         only : sound, eps, k_tilde
  USE mesh,        only : ibeg, iend, jbeg, jend, im, jm, xmax, ymax, xc, yc
  USE variables,   only : rho, vx1, vx2, prs, i_e, k_t, cs
  use variables,   only  : prim, cs
  
  IMPLICIT NONE

!==============================================================================!  
! Description:                                                                 !
!   The initial conditions are implemented here and the module will fullfil    !
! the block data
!                                                                              !
! Author: Victor Mourão Roque                                                  !
! last modification:                                                           !
! date: 08/08/2017 ; by Victor Mourao Roque                                    ! 
!                                                                              !
!==============================================================================!

! Private scope by default
!
  PRIVATE

  PUBLIC :: initia
  
  CONTAINS



!----------------------------------------------------------------------
 
 SUBROUTINE initia()  
!----------------------------------------------------------------------
! Purpose: to set initial conditions for the primtive variables d, u, p  
!and for the conserved variables cs
!teste hidrodinâmico: 
! KHI
! gamma = 1.3333333d0 length = (1.,2.) Nx2N => 256x512 ou 512x1024 
! Periodic Boundaries Conditions
! t = 3 (ou 6)         
!
  IMPLICIT NONE

! Input/Output Variables           


! Local Variables
  integer(i4b)          :: i,j
  REAL(DP), PARAMETER :: a = 0.01_dp , Vs = 0.5_dp, d0 = 0.505_dp , d1 = 0.495
  REAL(DP), PARAMETER :: A0 = 0.1_dp, sigma = 0.1_dp
  
  prim(:,:,4) = 1.d0


!==============================================================================
  
  do j = jbeg, jend
    do i = ibeg, iend

      IF(yc(j) .LE. 0.d0) THEN
        prim(i,j,1) = d0 - d1*tanh((yc(j)+.5_dp)/a)
        prim(i,j,2) = -Vs*tanh((yc(j)+.5_dp)/a)
        prim(i,j,3) = -A0*Vs*sin(2.d0*Pi_D*xc(i))*EXP(-(yc(j)+.5d0)**2/sigma)
      else
        prim(i,j,1) = d0 + d1*tanh((yc(j)-.5_dp)/a)
        prim(i,j,2) = Vs*tanh((yc(j)-.5_dp)/a)
        prim(i,j,3) = A0*Vs*sin(2*Pi_D*xc(i))*EXP(-(yc(j)-.5_dp)**2/sigma)
      ENDIF

      cs(i,j) = sound(prim(i,j,1),1.d0)
      prim(i,j, i_e) = eps(prim(i,j,rho),1.d0)
      prim(i,j, k_t) = k_tilde(prim(i,j,rho),prim(i,j,i_e))

    ENDDO       
  ENDDO

!==============================================================================



!==============================================================================

!==============================================================================
  
END SUBROUTINE



END MODULE init_cond

