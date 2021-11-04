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
!teste hidrodinâmico: Lax problem
! [1]   A. Suresh and H. T. huynh, 
! Journal of Computational Physics 136, 83 (1997).
! gamma = 1.4d0, t = .32, cells = 100, lenght = 2.d0 [-1,1]
! or
![1]    R. Borges, M. Carmona, B. Costa, and W. S. Don,
! Journal of Computational Physics 227, 3191 (2008).
! gamma = 1.4d0, t = 1.3, cells = 200, lenght = 10.d0
! transmissive BC
!
  IMPLICIT NONE

! Input/Output Variables           
! REAL(DP), DIMENSION(-ng:cells(1)+ng,-ng:cells(2)+ng),INTENT(OUT) :: d,u,v,p,cs

! Local Variables
  integer(i4b)          :: i,j
  real(dp), dimension(:), allocatable :: d, u, v, p, css
  
  call initialize_gamma()

!==============================================================================
!   allocate(d(jm))
!   allocate(u(jm))
!   allocate(v(jm))
!   allocate(p(jm))
!   allocate(css(jm))

!    DO j = jbeg,jend 
!      !set initial values in left section of domaim
!     IF(yc(j) .LE. 0.d0) THEN
!       d(j) = .445d0 
!       u(j) = 0.d0
!       v(j) = 0.698d0
!       p(j) = 3.528d0
!       css(j) = sound(d(j),p(j))
!     else
!     !set initial values in right section of domaim
!       d(j) = .5d0
!       u(j) = 0.d0
!       v(j) = 0.d0
!       p(j) = .571d0
!       css(j) = sound(d(j),p(j))
!     endif
! ENDDO  



! do i = ibeg, iend

!   prim(i,:,rho) = d(:)
!   prim(i,:,vx1) = u(:)
!   prim(i,:,vx2) = v(:)
!   prim(i,:,prs) = p(:)
!   cs(i,:) = css(:)
! enddo

!==============================================================================

  allocate(d(im))
  allocate(u(im))
  allocate(v(im))
  allocate(p(im))
  allocate(css(im))

  DO i = ibeg, iend
     !set initial values in left section of domaim
    IF(xc(i) .LE. 0.d0) THEN
      d(i) = .445d0 
      u(i) = 0.698d0
      v(i) = 0.d0
      p(i) = 3.528d0
      css(i) = sound(d(i),p(i))
    else
    !set initial values in right section of domaim
      d(i) = .5d0
      u(i) = 0.d0
      v(i) = 0.d0
      p(i) = .571d0
      css(i) = sound(d(i),p(i))
    endif
ENDDO  

do j = jbeg, jend
  prim(:,j,rho) = d(:)
  prim(:,j,vx1) = u(:)
  prim(:,j,vx2) = v(:)
  prim(:,j,prs) = p(:)
  cs(:,j) = css(:)
enddo
!==============================================================================



  deallocate(p,v,u,d,css)
  
END SUBROUTINE



END MODULE init_cond

