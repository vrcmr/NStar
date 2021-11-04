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


MODULE sources
  use constants, only : dp, i4b
  use variables, only : prim 


  IMPLICIT NONE


!==============================================================================!  
! Description:                                                                 !
!   < Say what this module contains >                                          !
!                                                                              !
! Author: Victor Mourão Roque                                                  !
! last modification:                                                           !
! date: 02/09/2021 ; by Victor Mourao Roque                                    !
!                                                                              !
! Refs:                                                                        !
!                                                                              !
!                                                                              !
!==============================================================================!


  PRIVATE  ! Private scope by default


! Public members
  PUBLIC :: source, initialize_source

! Module variables
  real(dp), save, public :: grav  = -1.d-2

CONTAINS

SUBROUTINE initialize_source()
  use parameters, only : get_parameter

  IMPLICIT NONE

  call get_parameter("gravity", grav)

END SUBROUTINE initialize_source


! 	subroutine source(i,j,up,Sg)
subroutine source(i,j,up)
use mesh,       only : im, jm, ibeg, iend, jbeg, jend
use variables,  only : nc, rho,vx1,vx2,prs
USE EoS,        only : eps
implicit none
real(dp), dimension(1,1,nc), intent(inout) :: up
! 	real(dp), dimension(1,1,nc), intent(out) :: Sg
! Local variables
!
integer(i4b) :: i,j
  ! real(dp), parameter :: gravity = - 2.99792458d3! s^-1 ou 0.01 km^-1
  real(dp) :: rhoh, W2
  real(dp), dimension(1,1,nc) :: Sg


Sg(1,1,1) = 0
Sg(1,1,2) = 0

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j) 
!$OMP DO 
W2 = 1.d0/(1 - prim(i,j,vx2)**2 - prim(i,j,vx1)**2)
rhoh = prim(i,j,rho) * (1.d0 + eps(prim(i,j,rho), prim(i,j,prs))) &
     + prim(i,j,prs)
! DO j = jbeg,jend
!     DO i = ibeg, iend
Sg(1,1,3) = grav * rhoh * W2  ! grav*prim(i,j,rho)  !
Sg(1,1,4) = grav*prim(i,j,vx2) * rhoh * W2 !grav*prim(i,j,rho)*prim(i,j,vx2) !
! 		ENDDO
! ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

   up = up + Sg


  
  end subroutine source

END MODULE sources
