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


MODULE variables
  use constants , only : i4b, dp

IMPLICIT NONE

!==============================================================================!  
! Description:                                                                 !
!   Module to inicialize and  store the variables block used in the code       !
!                                                                              !
! Author: Victor Mourão Roque                                                  !
! last modification:                                                           !
! date: 06/08/2017 ; by Victor Mourao Roque                                    ! 
!                                                                              !
!==============================================================================!


! By default everthing is private
!
PRIVATE


!type block_data
! 
!==============================================================================!  
! Description:                                                                 !
! por estar fazendo um primeiro teste em duas dimensões, essas matrizes terão  !
! somente 3 dimensões. Na primeira é definida a variável e nas outras duas a   !
! posição na rede                                                              !
!                                                                              !
! não sei qual seria a melhor abordagem, fazer dentro do type vários vetores ou!
! fazer um vetor type que contem em cada elemento todas as informações, ou     ! 
! seja, tenho um vetor com todas as células e cada elemento contem todas a     !
! informação da célula                                                         ! 
!                                                                              !
! Por hora vou trabalhar com dois vetores alocáveis, um para variáveis prim. e !
! outro para variáveis conservadas                                             !
!                                                                              !
! Author: Victor Mourão Roque                                                  !
! last modification:                                                           !
! date: 06/08/2017 ; by Victor Mourao Roque                                    ! 
!                                                                              !
!==============================================================================!  
! conservative variables  
! (1,:,:) => D 
! (2,:,:) => Sx ; (3,:,:) => Sy ;  (4,:,:) => \Tau ;
!  real(dp), dimension(:,:,:), allocatable       :: u
! primitive variables  
! (1,:,:) => \rho
! (2,:,:) => v_x ; (3,:,:) => v_y ; (4,:,:) => prs ; (5,:,:) => c_s
!  real(dp), dimension(:,:,:), allocatable       :: w
! Numerical flux ====>>> não tenho certeza que seja necessário. Melhor alocar um 
! vetor no RK. Mesma coisa acontece com o fluxo físico, que só é usado no split
! e assume o carater unidimensional
!real, dimension(:,:,:), allocatable       :: nflux
!
!end type block_data	

! type(block_data),save, public   :: var_block

! number of conservative variables
    integer, parameter, public  :: nc = 4
! number of primitive variables + internal energy + k_til
    integer, parameter, public  :: np = 6

! primitive variables and sound velocity
! (:,:,1) => \rho
! (:,:,2) => v_x ; (:,:,3) => v_y ; (:,:,4) => prs
  real(dp), dimension(:,:,:), allocatable  :: prim
  real(dp), dimension(:,:)  , allocatable  :: cs

! conservative variables
! (:,:,1) => D 
! (:,:,2) => Sx ; (:,:,3) => Sy ; (:,:,4) => \Tau ;
  real(dp), dimension(:,:,:), allocatable  :: cons 

! primitive variables index
!
integer(i4b), parameter, public  :: rho = 1, vx1 = 2, vx2 = 3, prs = 4
integer(i4b), parameter, public  :: i_e = 5, k_t = 6

! conservative variables index
!
integer(i4b), parameter, public  :: D = 1, Sx1 = 2, Sx2 = 3, tau = 4


public :: initialize_variables, finalize_variables
!public :: print_prim_var, print_cons_var, print_var_bc
public :: prim, cs, cons

 contains
!
!===============================================================================
!
! subroutine INITIALIZE_VARIABLES:
! -------------------------------
!
!   Subroutine allocates memory for all variable arrays.
!
!   Arguments:
!
!     im, jm, km - domain dimensions along the X, Y, and Z directions;
!
!===============================================================================
!


subroutine initialize_variables()
  use mesh , only : im, jm

  IMPLICIT NONE


! allocate conservative variables
!
 allocate(cons(im,jm,nc))
! allocate primitive variables and sound speed
!
 allocate(prim(im,jm,np))
 allocate(cs(im,jm))

end subroutine initialize_variables



!===============================================================================
!
! subroutine FINALIZE_VARIABLES:
! -----------------------------
!
!   Subroutine releases memory used by variable arrays.
!
!===============================================================================
!
 
  subroutine finalize_variables()
  IMPLICIT NONE
!
!-------------------------------------------------------------------------------
!
if (allocated(prim))  deallocate(prim)
if (allocated( cs ))  deallocate( cs ) 
if (allocated(cons))  deallocate(cons)

!-------------------------------------------------------------------------------
!
  end subroutine finalize_variables

  


! subroutine print_var_bc(idx)
!   use mesh, only : im,jm,ng,jcell
!   implicit none
!   integer(i4b), intent(in) :: idx 
!   integer(i4b) :: i,j
!   real(dp), dimension(im,jm) :: d
      
!     d(:,:) = prim(:,:,idx)

! !     do j = jm, 1, -1
! !       if((j == ng+1) ) then
! !         write(*,"(3f9.0,5x,a,8f9.0,5x,a,3f9.0)") d(1,j),d(2,j),d(3,j),'I', & 
! !                             d(4,j),d(5,j),d(6,j),d(7,j),d(8,j),d(9,j),d(10,j),d(11,j),'I', &
! !                             d(12,j),d(13,j),d(14,j)
! !         write(*,*) "                               --------------------------",&
! !         "-----------------------------------------------------               "
! !       elseif((j == jcell+ng) ) then
! !         write(*,*) "                               --------------------------",&
! !         "-----------------------------------------------------               "
! !         write(*,"(3f9.0,5x,a,8f9.0,5x,a,3f9.0)") d(1,j),d(2,j),d(3,j),'I', & 
! !                             d(4,j),d(5,j),d(6,j),d(7,j),d(8,j),d(9,j),d(10,j),d(11,j),'I', &
! !                             d(12,j),d(13,j),d(14,j)
        
! !       else
! !       write(*,"(3f9.0,5x,a,8f9.0,5x,a,3f9.0)") d(1,j),d(2,j),d(3,j),'I', & 
! !                             d(4,j),d(5,j),d(6,j),d(7,j),d(8,j),d(9,j),d(10,j),d(11,j),'I', &
! !                             d(12,j),d(13,j),d(14,j)
! !       endif
      
! !     enddo

!  do j = jm, 1, -1
!       if((j == ng+1) ) then
!         write(*,"(2f9.2,5x,a,6f9.2,5x,a,2f9.2)") d(1,j),d(2,j),'I',d(3,j), & 
!                             d(4,j),d(5,j),d(6,j),d(7,j),d(8,j),'I',d(9,j),d(10,j)
!         write(*,*) "                               --------------------------",&
!         "-----------------------------------------------------               "
!       elseif((j == jcell+ng) ) then
!         write(*,*) "                               --------------------------",&
!         "-----------------------------------------------------               "
!         write(*,"(2f9.2,5x,a,6f9.2,5x,a,2f9.2)") d(1,j),d(2,j),'I',d(3,j), & 
!                             d(4,j),d(5,j),d(6,j),d(7,j),d(8,j),'I',d(9,j),d(10,j)
!       else
!       write(*,"(2f9.2,5x,a,6f9.2,5x,a,2f9.2)") d(1,j),d(2,j),'I',d(3,j), & 
!                             d(4,j),d(5,j),d(6,j),d(7,j),d(8,j),'I',d(9,j),d(10,j)
!       endif
      
!     enddo
!     write(*,*)
! end subroutine print_var_bc

!   subroutine print_cons_var(idx)
!     use mesh, only : jm,im
!     implicit none
!     integer(i4b), intent(in) :: idx 
!     integer :: i,j

!     do j = 1,jm
!         write(*,*) (cons(i,j,idx), i = 1,im)
!     enddo
!     write(*,*)
  
!   end subroutine print_cons_var

!   subroutine print_prim_var(idx)
!     use mesh, only : jm,im
!     implicit none
!     integer(i4b), intent(in) :: idx 
!     integer :: i,j

!     do j = 1,jm
!         write(*,*) (prim(i,j,idx), i = 1,im)
!     enddo
!     write(*,*)

!   end subroutine


END MODULE variables
