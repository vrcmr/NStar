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


module boundaries
  use mesh,  only : ng
  use constants,  only :  i4b 
  use variables,  only :  prim, cs
  use variables,  only :  rho, vx1, vx2, prs, i_e, k_t

  IMPLICIT NONE


!==============================================================================!  
! Description:                                                                 !
!   Module responsable to choose the boundaries conditions                     !
!                                                                              !
! Author: Victor Mourão Roque                                                  !
! last modification:                                                           !
! date: 13/07/2017 ; by Victor Mourao Roque                                    ! 
!                                                                              !
! Melhorias:                                                                   !
! Implementar Inflow e Outflow. Transmissiva não me parece muito aplicável     !
! Outflow olhar e.g. userguide do PLUTO
!                                                                              !
! Refs:                                                                        !
!                                                                              !
!                                                                              !
!==============================================================================!
 
  PRIVATE 

! Local Variables
  procedure(periodic_x_beg), pointer, save :: bc_x_beg => null()
  procedure(periodic_x_end), pointer, save :: bc_x_end => null()
  procedure(periodic_y_beg), pointer, save :: bc_y_beg => null()
  procedure(periodic_y_end), pointer, save :: bc_y_end => null()
  
   
  PUBLIC :: boundary2d
  PUBLIC :: initialize_boundaries, finalize_boundaries


CONTAINS

  subroutine initialize_boundaries()
    use parameters,  only : get_parameter
  implicit none

 ! BC strings
 ! 
  character(len = 32), save     :: xbegbc  = "periodic"
  character(len = 32), save     :: xendbc  = "periodic"
  character(len = 32), save     :: ybegbc  = "periodic"
  character(len = 32), save     :: yendbc  = "periodic"
  
                          
! Obtain the type of BC
!  
  call get_parameter("xbegbc", xbegbc)
  call get_parameter("xendbc", xendbc)
  call get_parameter("ybegbc", ybegbc)
  call get_parameter("yendbc", yendbc)

! select the routine to boundary condition in x_beg
  select case(trim(xbegbc))
    case("periodic", "PERIODIC")
      
      bc_x_beg => periodic_x_beg

    case("transmissive", "TRANSMISSIVE")

      bc_x_beg => transm_x_beg

    case("reflective", "REFLECTIVE")

      bc_x_beg => refl_x_beg

     case default

      bc_x_beg => periodic_x_beg

  end select    

  ! select the routine to boundary condition in x_end
  select case(trim(xendbc))
    case("periodic", "PERIODIC")
    
      bc_x_end => periodic_x_end
    
    case("transmissive", "TRANSMISSIVE")

      bc_x_end => transm_x_end

    case("reflective", "REFLECTIVE")

      bc_x_end => refl_x_end

    case default

      bc_x_end => periodic_x_end

  end select

  ! select the routine to boundary condition in y_beg
  select case(trim(ybegbc))
    case("periodic", "PERIODIC")
      
      bc_y_beg => periodic_y_beg

    case("transmissive", "TRANSMISSIVE")

      bc_y_beg => transm_y_beg

    case("reflective", "REFLECTIVE")

      bc_y_beg => refl_y_beg

    case default

      bc_y_beg => periodic_y_beg

  end select    

  ! select the routine to boundary condition in x_end
  select case(trim(yendbc))
    case("periodic", "PERIODIC")

      bc_y_end => periodic_y_end

    case("transmissive", "TRANSMISSIVE")

      bc_y_end => transm_y_end

    case("reflective", "REFLECTIVE")

      bc_y_end => refl_y_end

    case default
      
      bc_y_end => periodic_y_end

  end select

  end subroutine initialize_boundaries

  
  subroutine finalize_boundaries()
  implicit none
!
! release the procedure pointers
!   
   nullify(bc_x_beg) ; nullify(bc_x_end)
   nullify(bc_y_beg) ; nullify(bc_y_end)

  
  end subroutine finalize_boundaries


  subroutine boundary2d()
    
    call bc_x_beg() ; call bc_x_end()
    call bc_y_beg() ; call bc_y_end()
  	
  end subroutine boundary2d

!==============================================================================!  
! Description:                                                                 !
!   Subroutines describe a PERIODIC BOUNDARY conditions. The ghost cells       !
!   repeat all variables of the oposite side.                                  !
!                                                                              !
! Author: Victor Mourão Roque                                                  !
! last modification:                                                           !
! date: 17/07/2017 ; by Victor Mourao Roque                                    ! 
!                                                                              !
!==============================================================================!


  subroutine periodic_x_beg() 
    use mesh,       only :  iend
  implicit none
    
    prim(1:ng,:,:) = prim(iend-ng+1:iend,:,:)
    cs(1:ng,:)     = cs(iend-ng+1:iend,:)

  end subroutine periodic_x_beg


  subroutine periodic_x_end()
    use mesh,       only :  ibeg, iend,im
  implicit none

    prim(iend+1:im,:,:) = prim(ibeg:ng+ng,:,:)
    cs(iend+1:im,:)     = cs(ibeg:ng+ng,:)

  end subroutine periodic_x_end

  subroutine periodic_y_beg() 
    use mesh,       only :  jend
    implicit none

     prim(:,1:ng,:) = prim(:,jend-ng+1:jend,:)
     cs(:,1:ng)     = cs(:,jend-ng+1:jend)

  end subroutine periodic_y_beg

  subroutine periodic_y_end() 
    use mesh,       only :  jbeg, jend, jm
    implicit none

    prim(:,jend+1:jm,:) = prim(:,jbeg:ng+ng,:)
    cs(:,jend+1:jm)     = cs(:,jbeg:ng+ng)
     
  end subroutine periodic_y_end


!==============================================================================!  
! Description:                                                                 !
!   Subroutines describe a REFLECTIVE BOUNDARY conditions. The NORMAL velocity !
! in ghost cells has a opposite sign while the other variables are symmetrized !
! across the boundary remains equal                                            !
!                                                                              !
! Author: Victor Mourão Roque                                                  !
! last modification:                                                           !
! date: 17/07/2017 ; by Victor Mourao Roque                                    ! 
!                                                                              !
!==============================================================================!

  subroutine refl_x_beg()
    use mesh,       only :  ibeg
  implicit none
! local variables
  integer(i4b) :: i
 
     do concurrent(i=1:ng)
      
      prim(i,:,rho) =  prim(ibeg+ng-i,:,rho)
      prim(i,:,vx1) = -prim(ibeg+ng-i,:,vx1)
      prim(i,:,vx2) =  prim(ibeg+ng-i,:,vx2)
      prim(i,:,prs) =  prim(ibeg+ng-i,:,prs)
      prim(i,:,i_e) =  prim(ibeg+ng-i,:,i_e)
      prim(i,:,k_t) =  prim(ibeg+ng-i,:,k_t)
      cs(i,:)       =  cs(ibeg+ng-i,:)

    end do

  end subroutine refl_x_beg


  subroutine refl_x_end()
    use mesh,       only :  iend

  implicit none
! local variables
  integer(i4b) :: i

    do concurrent(i=1:ng)
      prim(iend+i,:,rho) =  prim(iend-i+1,:,rho)
      prim(iend+i,:,vx1) = -prim(iend-i+1,:,vx1)
      prim(iend+i,:,vx2) =  prim(iend-i+1,:,vx2)
      prim(iend+i,:,prs) =  prim(iend-i+1,:,prs)
      prim(iend+i,:,i_e) =  prim(iend-i+1,:,i_e)
      prim(iend+i,:,k_t) =  prim(iend-i+1,:,k_t)

      cs(iend+i,:)       =  cs(iend-i+1,:)
    end do

  end subroutine refl_x_end


  subroutine refl_y_beg()
    use mesh,       only :  jbeg
  implicit none
! local variables
  integer(i4b) :: j
 
     do concurrent(j=1:ng)
      
      prim(:,j,rho) =  prim(:,jbeg+ng-j,rho)
      prim(:,j,vx1) =  prim(:,jbeg+ng-j,vx1)
      prim(:,j,vx2) = -prim(:,jbeg+ng-j,vx2)
      prim(:,j,prs) =  prim(:,jbeg+ng-j,prs)
      prim(:,j,i_e) =  prim(:,jbeg+ng-j,i_e)
      prim(:,j,k_t) =  prim(:,jbeg+ng-j,k_t)
      cs(:,j)       =  cs(:,jbeg+ng-j)

    end do

  end subroutine refl_y_beg


  subroutine refl_y_end()
    use mesh,       only :  jend
  implicit none
! local variables
  integer(i4b) :: j

    do concurrent(j=1:ng)
      prim(:,jend+j,rho) =  prim(:,jend-j+1,rho)
      prim(:,jend+j,vx1) =  prim(:,jend-j+1,vx1)
      prim(:,jend+j,vx2) = -prim(:,jend-j+1,vx2)
      prim(:,jend+j,prs) =  prim(:,jend-j+1,prs)
      prim(:,jend+j,i_e) =  prim(:,jend-j+1,i_e)
      prim(:,jend+j,k_t) =  prim(:,jend-j+1,k_t)
      cs(:,jend+j)       =  cs(:,jend-j+1)
    end do

  end subroutine refl_y_end

  

!==============================================================================!  
! Description:                                                                 !
!   Subroutines describe a TRANSMISSIVE BOUNDARY conditions. Just repeat the   !
! first/last cell in the boundaries                                            !
!                                                                              !
! Author: Victor Mourão Roque                                                  !
! last modification:                                                           !
! date: 17/07/2017 ; by Victor Mourao Roque                                    ! 
!                                                                              !
!==============================================================================!

    subroutine transm_x_beg()
      use mesh,       only :  ibeg
    implicit none
    integer :: i
    
    do concurrent(i=1:ng)
      prim(i,:,:) =  prim(ibeg,:,:)
      cs(i,:)       =  cs(ibeg,:)
    enddo

  end subroutine transm_x_beg


  subroutine transm_x_end()
    use mesh,       only :  iend
  implicit none
  integer  :: i

    do concurrent(i = 1:ng)
      prim(iend+i,:,:) =  prim(iend,:,:)
      cs(i,:)     =  cs(iend,:)   
    enddo 

  end subroutine transm_x_end

  subroutine transm_y_beg()
      use mesh,       only :  jbeg
    implicit none
    integer :: j
    
    do concurrent(j=1:ng)
      prim(:,j,:) =  prim(:,jbeg,:)
      cs(:,j)       =  cs(:,jbeg)
    enddo

  end subroutine transm_y_beg


  subroutine transm_y_end()
    use mesh,       only :  jend
  implicit none
  integer  :: j

    do concurrent(j = 1:ng)
      prim(:,jend+j,:) =  prim(:,jend,:)
      cs(j,:)     =  cs(:,jend)   
    enddo 

  end subroutine transm_y_end


end module boundaries
