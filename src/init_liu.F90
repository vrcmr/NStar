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
  USE constants,   only : i4b, dp
  USE EoS,         only : sound, k_tilde, eps
  USE mesh,        only : ibeg, iend, jbeg, jend, im, jm, xmax, ymax, xc, yc
  USE variables,   only : rho, vx1, vx2, prs, cs, i_e, k_t
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
 
 SUBROUTINE initia()  !(d, u, v, p, cs)
!----------------------------------------------------------------------
! Purpose: to set initial conditions for the primtive variables d, u, p  
!and for the conserved variables cs
!teste hidrodinâmico:  Liu 
! gamma = 1.4d0, t~0.23d0, xlength = 1.d0, ylength = 1.d0, cells = 100-800 
! CFL = 0.6 BC transmissive
! [1] G. A. Gerolymos, D. Sénéchal, and I. Vallet, 
!     Journal of Computational Physics 228, 8481 (2009).
!
  IMPLICIT NONE

! Input/Output Variables           
! REAL(DP), DIMENSION(-ng:cells(1)+ng,-ng:cells(2)+ng),INTENT(OUT) :: d,u,v,p,cs

! Local Variables
  integer(i4b)          :: i,j

!==============================================================================

  DO j = jbeg, jend
    DO i = ibeg, iend
     !set initial values in left section of domaim
      IF(xc(i) >= xmax*0.d0 .AND. yc(j) .LT. ymax*0.d0  ) THEN
        prim(i,j,rho) = 3.d0
        prim(i,j,vx1) = 0.75d0
        prim(i,j,vx2) = -0.5d0
        prim(i,j,prs) = 1.d0

        cs(i,j) = sound(3.d0,1.d0)
        prim(i,j, i_e) = eps(3.d0,1.d0)
        prim(i,j, k_t) = k_tilde(3.d0,1.d0)

      ELSEIF(xc(i) < xmax*0.d0 .AND. yc(j) >= ymax*0.d0  ) THEN
        prim(i,j,rho) = 2.d0
        prim(i,j,vx1) = -0.75d0
        prim(i,j,vx2) = 0.5d0
        prim(i,j,prs) = 1.d0

        cs(i,j) = sound(2.d0,1.d0)
        prim(i,j, i_e) = eps(2.d0,1.d0)
        prim(i,j, k_t) = k_tilde(2.d0,1.d0)

      ELSEIF(xc(i) >= xmax*0.d0 .AND. yc(j) >= ymax*0.d0  ) THEN
        prim(i,j,rho) = 1.d0
        prim(i,j,vx1) = -0.75d0
        prim(i,j,vx2) = -0.5d0
        prim(i,j,prs) = 1.d0

        cs(i,j) = sound(1.d0,1.d0)
        prim(i,j, i_e) = eps(1.d0,1.d0)
        prim(i,j, k_t) = k_tilde(1.d0,1.d0)

      ELSE
        prim(i,j,rho) = 1.d0
        prim(i,j,vx1) = 0.75d0
        prim(i,j,vx2) = 0.5d0
        prim(i,j,prs) = 1.d0

        cs(i,j) = sound(1.d0,1.d0)
        prim(i,j, i_e) = eps(1.d0,1.d0)
        prim(i,j, k_t) = k_tilde(1.d0,1.d0)

      ENDIF
    ENDDO       
  ENDDO

END SUBROUTINE



END MODULE init_cond

