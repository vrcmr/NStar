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


MODULE numflux
  use constants,   only  : dp, i4b

  IMPLICIT NONE


!==============================================================================!  
! Description:                                                                 !
!   This module uses the reconstruction and the flux split schemes and         !
! calculates the numerical fluxes to update the conservative variables         !
!                                                                              !
! Author: Victor Mourão Roque                                                  !
! last modification:                                                           !
! date: 19/07/2017 ; by Victor Mourao Roque                                    ! 
!                                                                              !
! Refs:                                                                        !
!       [1] D. S. Balsara and C.-W. Shu,                                       !
!           Journal of Computational Physics 160, 405–452 (2000).              !
!       [2] D. Radice and L. Rezzolla, arXiv astro-ph.IM, (2012).              !
!       [3] A. Mignone, P. Tzeferacos, and G. Bodo,                            !
!           Journal of Computational Physics 229, 5896–5920 (2010).            !
!                                                                              !
!==============================================================================!

  PRIVATE  ! Private scope by default


! Module constants
  !< Parameter declarations >


! Module variables
  !< Variable declarations >


! Public members
  PUBLIC :: num_flux, primitive_mean_values


CONTAINS

!==============================================================================!
!  Melhorias: Adicionar um initialize para escolher entre fazer a reconstrução !
! nos espaço das características ou no espaço dos fluxos                       !
!                                                                              !
! description of the Variables:                                                !
! pflux => physical flux                                                       !
! nflux => numerical flux                                                      !
! prim_mean => primitive variable mean                                         !
! lambda_i => eigenvalue in i+1/2 (interface)                                  !
! lambda_c => eigenvalue in central position of cell i                         !
! lev => left eigenvector 1 => plus / 2 => 0 / 3 => minus  in i+1/2            !
! rev => right eigenvector 1 => plus / 2 => 0 / 3 => minus  in i+1/2           !
!=============================================================================!
!==============================================================================!
  
  subroutine num_flux(coord, prim, cs, cons, nflux)
    use mesh,           only : r_s, cells
    use variables,      only : np, nc
    use physics,        only : phys_flux, eigenval, eigenvec, test_vectorsHD
    use fluxsplit,      only : flux_split, LLF
    use reconstruction, only : rec_cells
  
  !use reconstruction
  use physics, only : phys_flux, eigenval
  implicit none
  integer(i4b),             intent(in)     :: coord
  real(dp), dimension(cells(coord),np), intent(in)     :: prim
  real(dp), dimension(cells(coord),nc), intent(in)     :: cons
  real(dp), dimension(cells(coord)),    intent(in)     :: cs
  real(dp), dimension(cells(coord),nc), intent(inout)  :: nflux
! Local variables
  real(dp), dimension(:),     allocatable   :: cs_mean
  real(dp), dimension(:,:),   allocatable   :: prim_mean
  real(dp), dimension(:,:),   allocatable   :: lambda_i, lambda_c
  real(dp), dimension(:,:),   allocatable   :: pflux, nflux_r, nflux_l
  real(dp), dimension(:,:,:), allocatable   :: sflux_r, sflux_l
  real(dp), dimension(:,:,:), allocatable   :: lev, rev

 integer :: j


! Calculation of primitives simple mean eq.(16) of [2]
! alguns programas aplicam também a média de Roe, contudo segundo
! Zhang & MacFadyen (2006) o uso dessa não influência muito a qualidade 
! da resolução

  allocate( prim_mean(cells(coord),np) )
  allocate( cs_mean(cells(coord)) )

  call primitive_mean_values(coord, prim, cs, prim_mean, cs_mean)


!Eigenvalues and eigenvectors are computed in terms of prim_mean  
! Calculates de eigenvalues of interfaces i+1/2 values
!
  IF(.not. ALLOCATED(lambda_i)) allocate(lambda_i(cells(coord),nc))

  call eigenval(coord, prim_mean, cs_mean, lambda_i)

! Calculates de eigenvalues of central values
!
  IF(.not. ALLOCATED(lambda_c)) allocate(lambda_c(cells(coord),nc))

  call eigenval(coord, prim, cs, lambda_c)
 
 ! calculates left-eigenvector and right-eigenvector
  ALLOCATE(lev(cells(coord),nc,nc))
  ALLOCATE(rev(cells(coord),nc,nc))


  call eigenvec(coord, prim_mean, cs_mean, lambda_i, lev, rev)

!     call test_vectorsHD(coord, lev, rev)

  deallocate(cs_mean)
  deallocate(prim_mean)

! in order to compute nflux, we have to split pflux in a right-going, sflux_r 
!(f+), and a left-going, sflux_l (f−) look [2] page3
! split flux folowing eq (2.6a) of [1]
!
  IF(.not. ALLOCATED(pflux)) allocate(pflux(cells(coord),nc)) 

  call phys_flux(coord,prim,pflux)

  IF(.not. ALLOCATED(sflux_r)) allocate(sflux_r(cells(coord),nc,-r_s:r_s))
  IF(.not. ALLOCATED(sflux_l)) allocate(sflux_l(cells(coord),nc,-r_s:r_s))
  
  call flux_split(coord, pflux, cons, lambda_c, lambda_i, lev, sflux_r, sflux_l)

  deallocate(pflux)    ;  deallocate(lambda_c) 
  deallocate(lambda_i) ;  deallocate (lev)

  IF(.not. ALLOCATED(nflux_r)) allocate(nflux_r(cells(coord),nc))
  IF(.not. ALLOCATED(nflux_l)) allocate(nflux_l(cells(coord),nc))
  
  call rec_cells(coord,sflux_r,nflux_r)
  call rec_cells(coord,sflux_l,nflux_l)

  deallocate(sflux_r)  ;  deallocate(sflux_l)
  
  do j = 1, nc

    nflux(:,j) = ( nflux_r(:,1) + nflux_l(:,1) ) * rev(:,j,1) + &
                 ( nflux_r(:,2) + nflux_l(:,2) ) * rev(:,j,2) + &
                 ( nflux_r(:,3) + nflux_l(:,3) ) * rev(:,j,3) + &
                 ( nflux_r(:,4) + nflux_l(:,4) ) * rev(:,j,4) 

  enddo
                 
  deallocate(nflux_r)  ;  deallocate(nflux_l) ; deallocate(rev)

  end subroutine num_flux



!==============================================================================!
!  Subroutine needed to calculate the primitive variables on the interface of  !
! cells used in num_flux. The subroutine is called in the time integration.    !
! The primitive mean arrays are allocated and deallocated in the time          !
! integration also.                                                            !
!                                                                              !
!==============================================================================!

  subroutine primitive_mean_values(coord,prim,cs,p_mean,cs_mean)
    use mesh,       only  :  cells
    use variables,  only  :  np, i_e, k_t
    use EoS, only : eps
  implicit none
  integer(i4b),             intent(in)    :: coord
  real(dp), dimension(cells(coord)),   intent(in)    :: cs
  real(dp), dimension(cells(coord),np), intent(in)    :: prim
  real(dp), dimension(cells(coord)),   intent(out)   :: cs_mean
  real(dp), dimension(cells(coord),np), intent(out)   :: p_mean
! Local Variables
!  
  integer(i4b)  :: i, var

! $OMP parallel
! $OMP workshare

  !do i = cbeg(coord), cend(coord)
  do i = 1, cells(coord)-1
    do var = 1, np!-2
      p_mean(i,var) = ( prim(i+1,var) + prim(i,var) ) * 0.5d0
    enddo
    cs_mean(i) = DMAX1(cs(i+1),cs(i)) !(cs(i+1)+cs(i))/2.d0
    ! p_mean(i,i_e) = eps(p_mean(i,1), p_mean(i,4))
    ! p_mean(i,k_t) = ( prim(i+1, k_t) + prim(i,k_t) ) * 0.5d0
  enddo

  
! $OMP end workshare
! $OMP end parallel

  

  end subroutine primitive_mean_values
  

END MODULE numflux
