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
  USE EoS,         only : sound, density, pstar, k_tilde_p, eps
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
! Neutron Star combustion
! Condition made from hydrostatic equations using python
! p(x) = 2.64753683d1 - 2.66947193d0 * x - 6.90300373d-02 * x**2 +
!   9.82159699d-3 * x**3 - 7.08155943d-5 * x**4
!
  IMPLICIT NONE

! Input/Output Variables           


! Local Variables
  integer(i4b)          :: i,j
  ! real(dp), parameter :: rho_crit = 400.d0 ! 2.3e23 ! g/cm3
  

!==============================================================================
  
  do j = jbeg, jend
    do i = ibeg, iend

      prim(i,j,vx1) = 1.d-12
      prim(i,j,vx2) = 1.d-12

      prim(i,j,prs) = 2.64753683d1 - 2.66947193d0  * yc(j)     &
                                   - 6.90300373d-2 * yc(j)**2  &
                                   + 9.82159699d-3 * yc(j)**3  &
                                   - 7.08155943d-5 * yc(j)**4

      ! prim(i,j,prs) = 10.d0 * exp(- 0.01 *  rho_crit / 10 *yc(j))

      ! prim(i,j,prs) =20.5d0 - 12.5 *tanh(10.d0 * yc(j) - 10.d0)

      ! if ( yc(j) .le. 1.d0) then
      !   prim(i,j,prs) = 33.d0
      ! elseif ( yc(j) .gt. 12.d0  ) then
      !   prim(i,j,prs) = 1d-6
      ! else
      !   prim(i,j,prs) = 8.d0
      ! endif

      if (prim(i,j,prs) .lt. 1.d-6 .or. yc(j) .gt. 12.5d0) then
        prim(i,j,prs) = 1d-6
      endif

      if (xc(i)**2 + yc(j)**2 .le. 0.015625d0) then
        ! circle with 0.1 km radius
        prim(i,j,prs) = pstar
      endif


      prim(i,j,rho) = density(prim(i,j,prs))


      cs(i,j) = sound(prim(i,j,rho), prim(i,j,prs))

      prim(i,j, i_e) = eps(prim(i,j,rho),prim(i,j,prs))

      prim(i,j, k_t) = k_tilde_p(prim(i,j,rho),prim(i,j,prs))

        !k_tilde(prim(i,j,rho),prim(i,j,i_e))

    ENDDO       
  ENDDO

!==============================================================================



!==============================================================================

!==============================================================================
  
END SUBROUTINE



END MODULE init_cond

