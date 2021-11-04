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


MODULE fluxsplit 
  use constants,  only : i4b, dp
 

  IMPLICIT NONE


!==============================================================================!  
! Description:                                                                 !
! Use the different methods to split the fluxes into a positive part [sflux_r] !
! (right going) and a negative part [sflux_l] (left going) for a arbitrary     !
! stencil size (-rs:rs)                                                        ! 
!                                                                              !
! j => i+1/2  related to interface cells                                       !
! last modification:                                                           !
! date: 27/06/2014 ; by Victor Mourao Roque                                    !
! based on: [1] D. S. Balsara and C.-W. Shu,                                   !
!               Journal of Computational Physics 160, 405–452 (2000).          !
!           [2] R. Donat, Journal of Computational Physics 146, 58–81 (1998).  !
!           [3] A. Mignone, P. Tzeferacos, and G. Bodo,                        !
!               Journal of Computational Physics 229, 5896–5920 (2010).        ! 
!                                                                              ! 
! lev => left eigenvector                                                      !  
! lambda_i => eigenvalue at cell center                                        ! 
! lambda_j => eigenvalue at cell interface calculed with the mean values       ! 
! r_s => size of the stencil                                                   ! 
! pflux => physical flux                                                       ! 
! conserved => conserved variables                                             ! 
!                                                                              !
!==============================================================================!



  PRIVATE  ! Private scope by default


! Module constants
!  < Parameter declarations >


! Module variables
 
  procedure(LLF), pointer, save :: flux_split => null()


! Public members
  PUBLIC :: initialize_split, finalize_split, flux_split, LLF


CONTAINS

  subroutine initialize_split()
    use parameters    , only : get_parameter
    implicit none
! local variables
!
    character(len=32) :: split_scheme = "LLF"

  ! select the integration method
 call get_parameter("split_scheme", split_scheme)

  select case(trim(split_scheme))

  case('LLF','llf')

    flux_split => LLF

  case('GLF','glf')

    !flux_split => GLF
    flux_split => LLF
    write (*,"(1x,a)") "The GLF flux split is not implemented:" //       &
                           " The used scheme is:" // trim(split_scheme)

  case('RF','rf')
    flux_split => LLF
    write (*,"(1x,a)") "The Roe flux split is not implemented:" //       &
                           " The used scheme is:" // trim(split_scheme)
    !flux_split => RF

  end select

  end subroutine initialize_split


subroutine finalize_split()
    implicit none
!
! release the procedure pointers
!
    nullify(flux_split)
  
  end subroutine finalize_split


!  

  subroutine LLF(coord, pflux, cons, lambda_c, lambda_i, lev, sfl_r, sfl_l)
    use mesh,        only : cend, cbeg, cells, r_s
    use variables,   only : nc

!==============================================================================! 
! Description:                                                                 !   
!    Split the physical flux into a stencil                                    !
! Size of stencil 2r+1 = order of the reconst. scheme                          !                         
!     [s] = i-r,...,i+r / [s'] = 2i - [s] +1                                   ! 
!                                                                              !
! Author: Victor Mourão Roque                                                  !
! last modification:                                                           !
! date: 13/07/2017 ; by Victor Mourao Roque                                    !
!                                                                              !
! melhorias:                                                                   !
! no gcc 7.1 o vetor sfl_r/l pode ser automático tb                            !  
!==============================================================================!
  implicit none
  integer(i4b), intent(in)  :: coord 
  real(dp), dimension(cells(coord),nc),      intent(in)  :: pflux, cons
  real(dp), dimension(cells(coord),nc),      intent(in)  :: lambda_c, lambda_i
  real(dp), dimension(cells(coord),nc,nc),   intent(in)  :: lev
  real(dp), dimension(cells(coord), nc, -r_s:r_s), intent(out) :: sfl_r, sfl_l

!  	real(dp), dimension(:,:),   intent(in)  :: pflux, cons
!  	real(dp), dimension(:,:),   intent(in)  :: lambda_c, lambda_i
!   	real(dp), dimension(:,:,:), intent(out) :: sfl_r, sfl_l
    
  
  ! Local Variables
  integer(i4b) :: i,k,s
  REAL(DP) :: alpha_k
  REAL(DP), PARAMETER :: qsi = 1.1d0

  DO i = cbeg(coord)-1, cend(coord)
    DO k = 1, nc

      alpha_k =  qsi*DMAX1(MAXVAL(DABS(lambda_c(i-r_s:i+r_s+1,k))),  &
                                  DABS(lambda_i(i,k)))  

      DO s = -r_s,r_s

        sfl_r(i,k,s) =  lev(i,1,k) * 0.5d0 * ( pflux(i+s,1)    &
                                   + alpha_k * cons(i+s,1) ) + & 
                        lev(i,2,k) * 0.5d0 * ( pflux(i+s,2)    &
                                   + alpha_k * cons(i+s,2) ) + & 
                        lev(i,3,k) * 0.5d0 * ( pflux(i+s,3)    &
                                   + alpha_k * cons(i+s,3) ) + & 
                        lev(i,4,k) * 0.5d0 * ( pflux(i+s,4)    &
                                   + alpha_k * cons(i+s,4) )
        
        sfl_l(i,k,s) =  lev(i,1,k) * 0.5d0 * ( pflux(i+1-s,1)   &
                                   - alpha_k * cons(i+1-s,1)) + & 
                        lev(i,2,k) * 0.5d0 * ( pflux(i+1-s,2)   &
                                   - alpha_k * cons(i+1-s,2)) + &
                        lev(i,3,k) * 0.5d0 * ( pflux(i+1-s,3)   &
                                   - alpha_k * cons(i+1-s,3)) + &
                        lev(i,4,k) * 0.5d0 * ( pflux(i+1-s,4)   &
                                   - alpha_k * cons(i+1-s,4))

      ENDDO 
    ENDDO
  ENDDO

end subroutine LLF
   

END MODULE fluxsplit
