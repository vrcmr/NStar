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
  USE EoS,         only : sound, initialize_gamma
  USE mesh,        only : ibeg, iend, jbeg, jend, im, jm, xmax, ymax, xc, yc
  USE variables,   only : rho, vx1, vx2, prs, cs
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
!teste hidrodinâmico:  Sod's tube
! gamma = 1.4d0, t~0.2d0, length = 1.d0, tcells = 320, cells = 200
! [1] D. Radice and L. Rezzolla, arXiv.org astro-ph.IM, (2012).
!gamma = 1.4
!
  IMPLICIT NONE

! Input/Output Variables           
! REAL(DP), DIMENSION(-ng:cells(1)+ng,-ng:cells(2)+ng),INTENT(OUT) :: d,u,v,p,cs

! Local Variables
  integer(i4b)          :: i,j
  real(dp), dimension(:), allocatable :: d, u, v, p, css
  
  call initialize_gamma()

!==============================================================================
  allocate(d(jm))
  allocate(u(jm))
  allocate(v(jm))
  allocate(p(jm))
  allocate(css(jm))

  DO i = jbeg, jend

     !set initial values in left section of domaim
      IF(yc(i) .LE. ymax*.5) THEN
        d(i) = 1.d0
        u(i) = 0.d0
        v(i) = 0.d0
        p(i) = 1.d0
        css(i) = sound(d(i),p(i))
      else
    !set initial values in right section of domaim
        d(i) = 0.125
        u(i) = 0.d0
        v(i) = 0.d0
        p(i) = 0.1
        css(i) = sound(d(i),p(i))
      ENDIF
   ! ENDDO       
  ENDDO



do j = ibeg, iend

  prim(j,:,rho) = d(:)
  prim(j,:,vx1) = u(:)
  prim(j,:,vx2) = v(:)
  prim(j,:,prs) = p(:)
  cs(j,:) = css(:)
enddo

!==============================================================================

!  allocate(d(im))
!   allocate(u(im))
!   allocate(v(im))
!   allocate(p(im))
!   allocate(css(im))

!   DO i = ibeg,iend

!      !set initial values in left section of domaim
!       IF(xc(i) .LE. xmax*.5) THEN
!         d(i) = 1.d0
!         u(i) = 0.d0
!         v(i) = 0.d0
!         p(i) = 1.d0
!         css(i) = sound(d(i),p(i))
!       else
!     !set initial values in right section of domaim
!         d(i) = 0.125
!         u(i) = 0.d0
!         v(i) = 0.d0
!         p(i) = 0.1
!         css(i) = sound(d(i),p(i))
!       ENDIF
!    ! ENDDO       
!   ENDDO

! do j = jbeg, jend
!   prim(:,j,rho) = d(:)
!   prim(:,j,vx1) = u(:)
!   prim(:,j,vx2) = v(:)
!   prim(:,j,prs) = p(:)
!   cs(:,j) = css(:)
! enddo
!==============================================================================



  deallocate(p,v,u,d,css)
  
END SUBROUTINE



END MODULE init_cond

