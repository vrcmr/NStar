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
  USE EoS,         only : sound, initialize_EoS, k_tilde, eps
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
! RT
! gamma = 1.4d0 
! single-mode:
! Vy = 0.01[1 + cos(4πx)][1 + cos(3πy)]/4.
! length = {(-0.25,0.25),(-0.75,0.75)} Nx3N => 100x300 ou 200x600
! Boundaries Conditions: periodic on X and reflective on Y
! t = 12.75
! Ref: http://www.astro.princeton.edu/~jstone/Athena/tests/rt/rt.html
!
  IMPLICIT NONE

! Input/Output Variables           


! Local Variables
  integer(i4b)          :: i,j
  
  
  ! call initialize_EoS()

  


!==============================================================================
  
  do j = jbeg, jend
    do i = ibeg, iend

      IF(yc(j) .LE. 0.d0) THEN
        prim(i,j,rho) = 1.d0
      else
        prim(i,j,rho) = 2.d0
      ENDIF

      prim(i,j,vx1) = 0.d0

      prim(i,j,vx2) = 0.01d0*(1.d0 + cos(4.d0*PI_D*xc(i))) * &
                           (1.d0 + cos(3.d0*PI_D*yc(j)))/4.d0

      prim(i,j,prs) = 2.5d0 - prim(i,j,1)*0.1d0*yc(j)/.4d0

      cs(i,j) = sound(prim(i,j,rho),prim(i,j,prs))
      prim(i,j, i_e) = eps(prim(i,j,rho),prim(i,j,prs))
      prim(i,j, k_t) = k_tilde(prim(i,j,rho),prim(i,j,i_e))

    ENDDO       
  ENDDO

!==============================================================================



!==============================================================================

!==============================================================================
  
END SUBROUTINE



END MODULE init_cond

