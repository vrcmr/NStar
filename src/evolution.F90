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


MODULE evolution

  use constants,    only  :  i4b, dp

  IMPLICIT NONE


!==============================================================================!  
! Description:                                                                 !
!   Module responsable to choose the integration method, accomplish integration! 
! to evolve in time and determine the time step based on CFL coeficient        !
!                                                                              !
! Author: Victor Mourão Roque                                                  !
! last modification:                                                           !
! date: 13/07/2017 ; by Victor Mourao Roque                                    ! 
!                                                                              !
! Melhorias:                                                                   !
! atualizar o gcc para 7.1 e usar submodules para os esquemas para não o module!
! não ficar tão grande                                                         !
! Refs:                                                                        !
!                                                                              !
!                                                                              !
!==============================================================================!


  PRIVATE  ! Private scope by default


! Module constants
  


! Module variables
! Time evolution parameters
!
  real(dp), public, save :: cfl     = 2.5d-01

! Local Variables
  procedure(eulerian), pointer, save :: time_evolution => null()
  procedure(tstep_HD), pointer, save :: timestep => null()
  procedure(update_),   pointer, save :: update => null()


! Public members
  PUBLIC :: initialize_evolution,  finalize_evolution
  PUBLIC :: timestep, time_evolution


  interface
    module subroutine eulerian()
      use constants,  only  : i4b, dp
      use physics,    only  : prim2cons, cons2prim
      use mesh,       only  : im, iend, ibeg, jm, jend, jbeg, cells
      use variables,  only  : prim, np, cs
      use variables,  only  : cons, nc
      use numflux,    only  : num_flux
      use mesh,       only  : dt
      use boundaries, only  : boundary2d
      implicit none
      integer(i4b) :: i, j, coord
      real(dp), dimension(:,:,:), allocatable  :: cons_in
      real(dp), dimension(:,:,:), allocatable  :: updated
      real(dp), dimension(:,:,:), allocatable  :: nflux_x, nflux_y
      real(dp), dimension(:,:), allocatable :: nflux_1d, prim_1d, cons_1d
      real(dp), dimension(:),   allocatable :: cs_1d
    end subroutine eulerian
    
    module subroutine rgk2_2()
      use constants,  only  : i4b, dp
      use physics,    only  : prim2cons, cons2prim
      use mesh,       only  : im, iend, ibeg, jm, jend, jbeg, cells, dt
      use variables,  only  : prim, np, cs, cons, nc
      use numflux,    only  : num_flux
      use boundaries, only  : boundary2d
      implicit none
      integer(i4b) :: i, j, coord
      real(dp), dimension(:,:,:), allocatable  :: cons_in, cons_1
      real(dp), dimension(:,:,:), allocatable  :: updated
      real(dp), dimension(:,:,:), allocatable  :: nflux_x, nflux_y
      real(dp), dimension(:,:), allocatable :: nflux_1d, prim_1d, cons_1d
      real(dp), dimension(:),   allocatable :: cs_1d
    end subroutine rgk2_2

    module subroutine rgk3_3()
      use constants,  only  : i4b, dp
      use physics,    only  : prim2cons, cons2prim
      use mesh,       only  : im, iend, ibeg, jm, jend, jbeg, cells, dt
      use variables,  only  : prim, np, cs, cons, nc
      use numflux,    only  : num_flux
      use boundaries, only  : boundary2d
      implicit none
      integer(i4b) :: i, j, coord
      real(dp), dimension(:,:,:), allocatable  :: cons_in, cons_1
      real(dp), dimension(:,:,:), allocatable  :: updated
      real(dp), dimension(:,:,:), allocatable  :: nflux_x, nflux_y
      real(dp), dimension(:,:), allocatable :: nflux_1d, prim_1d, cons_1d
      real(dp), dimension(:),   allocatable :: cs_1d
    end subroutine rgk3_3

    module subroutine rgk4_5()
      use constants,  only  : i4b, dp
      use physics,    only  : prim2cons, cons2prim
      use mesh,       only  : im, iend, ibeg, jm, jend, jbeg, cells, dt
      use variables,  only  : prim, np, cs, cons, nc
      use numflux,    only  : num_flux
      use boundaries, only  : boundary2d
      implicit none
      integer(i4b) :: i, j, coord
      real(dp), dimension(:,:,:), allocatable  :: cons_in
      real(dp), dimension(:,:,:), allocatable  :: cons_1, cons_2, cons_3, cons_4
      real(dp), dimension(:,:,:), allocatable  :: updated1, updated2
      real(dp), dimension(:,:,:), allocatable  :: nflux_x, nflux_y
      real(dp), dimension(:,:), allocatable :: nflux_1d, prim_1d, cons_1d
      real(dp), dimension(:),   allocatable :: cs_1d
    end subroutine rgk4_5

  end interface

CONTAINS

  subroutine initialize_evolution()
    use parameters    , only : get_parameter
    IMPLICIT NONE
! local variables
!
    character(len=255) :: time_scheme = "rgk2_2"    
    character(len=32)  :: phys_eq = "HD"
    character(len=32)  :: source_term = "no"


  call get_parameter("phys_eq",phys_eq)

  select case(trim(phys_eq))
    case("HD","hd","hydro")

      timestep => tstep_HD

    case("RHD","rhd")

      timestep => tstep_RHD

    case default

      timestep => tstep_RHD

  end select      

  call get_parameter("source_term",source_term)

  select case(source_term)

  case("Y","Yes","on","ON")

    update => update_source

  case("N","No","off","OFF")

    update => update_

  case default
  
     update => update_

  endselect

    call get_parameter("cfl",cfl)

    if (cfl <= 0.0d+00) then
        write(*,*) "Incorrect value of the CFL coefficient!"
        write(*,*) "Setting to the default value cfl = 0.3 ..."
        write(*,*)

      cfl = 0.3d+00
  end if

! select the integration method
 call get_parameter("time_scheme", time_scheme)

  select case(trim(time_scheme))

  case('EULER','euler')

    time_evolution => eulerian

  case('rgk2_2','RGK2_2')

    time_evolution => rgk2_2

    case('rgk3_3','RGK3_3')

    time_evolution => rgk3_3

    case('rgk4_5','RGK4_5')

    time_evolution => rgk4_5

  case default

    write (*,"(1x,a)") "The selected time advance method is not " //       &
                           "implemented: " // trim(time_scheme)
!     time_evolution  => teste3

  end select

  end subroutine initialize_evolution

  subroutine finalize_evolution()
    implicit none
!
! release the procedure pointers
!
    nullify(time_evolution)
  
  end subroutine finalize_evolution


!==============================================================================! 
!                                                                              !
!        Purpose: to update the conserved and primitive variables              !
!                                                                              !
!==============================================================================! 

  subroutine update_(fi_x, fi_y, up) 
    use mesh,  only  :  ibeg, iend, dx, jbeg, jend, dy, im, jm
    use variables, only : nc
  implicit none
  real(dp), dimension(im,jm,nc), intent(in)   :: fi_x, fi_y
  real(dp), dimension(im,jm,nc), intent(out)  :: up 
! Local Variables
  integer(i4b)  :: i,j

!$OMP PARALLEL PRIVATE(i,j) SHARED(fi_x,fi_y,up)
!$OMP DO
!  DO CONCURRENT (j = 1:cells(2))
!     DO CONCURRENT (i = 1:cells(1))

  do j = jbeg, jend
    do i = ibeg, iend
      
    
      up(i,j,:) = - (-fi_x(i-1,j,:) + fi_x(i,j,:))/dx &
                  - (-fi_y(i,j-1,:) + fi_y(i,j,:))/dy

    enddo
  enddo

!$OMP END DO NOWAIT
!$OMP END PARALLEL

  end subroutine update_


  subroutine update_source(fi_x, fi_y, up) 
    use mesh,      only : ibeg, iend, dx, jbeg, jend, dy, im, jm
    use variables, only : nc
    use sources,   only : source
  implicit none
  real(dp), dimension(im,jm,nc), intent(in)   :: fi_x, fi_y
  real(dp), dimension(im,jm,nc), intent(out)  :: up 
! Local Variables
  integer(i4b)  :: i,j
!   REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: Sg

!$OMP PARALLEL PRIVATE(i,j) SHARED(fi_x,fi_y,up)
!$OMP DO
!  DO CONCURRENT (j = 1:cells(2))
!     DO CONCURRENT (i = 1:cells(1))

  do j = jbeg, jend
    do i = ibeg, iend
      
    
      up(i,j,:) = - (-fi_x(i-1,j,:) + fi_x(i,j,:))/dx &
                  - (-fi_y(i,j-1,:) + fi_y(i,j,:))/dy

      call source(i,j,up(i,j,:))
    enddo
  enddo
  
!$OMP END DO NOWAIT
!$OMP END PARALLEL

!   if(.not. allocated(Sg)) allocate(Sg(im,jm,nc))

!   call source(Sg)

!   up = up + Sg

  end subroutine update_source


  subroutine tstep_HD(timeout)
    use variables,    only  : prim, cs, vx1, vx2
    use mesh,         only  : NDIMS, ibeg, iend, jbeg, jend 
    use mesh,         only  : time, tend, step, dt, dx, dy
!==============================================================================! 
! Description:                                                                 !   
!    to apply the CFL condition to find a stable time step size DT             !
!    returns:  dt                                                              ! 
!                                                                              !
! Author: Victor Mourão Roque                                                  !
! last modification:                                                           !
! date: 13/07/2017 ; by Victor Mourao Roque                                    !                                                                               !
! Refs: A. Mignone, G. Bodo, S. Massaglia, T. Matsakos, O. Tesileanu,          !
!      C. Zanni, and A. Ferrari, arXiv.org astro-ph, 228 (2007).               !
!                                                                              !
!==============================================================================!

! Subroutine arguments
  real(dp),  intent(in)      :: timeout

!  local variables 
  real(dp), dimension(NDIMS) :: smax, ds
  real(dp)                   :: dt_l

      smax(1) = MAXVAL(abs(prim(ibeg:iend,jbeg:jend,vx1) & 
                   + cs(ibeg:iend,jbeg:jend)))
      smax(2) = MAXVAL(abs(prim(ibeg:iend,jbeg:jend,vx2) &   
                   + cs(ibeg:iend,jbeg:jend)))
   
   ds(1) = dx ; ds(2) = dy 

   dt_l = minval(ds/smax)

   !compute time step dt
   dt = cfl*dt_l
  
   
  ! Recalculates "dt" for early times (dt is reduced to compensate 
  !for approximate calculation of smax)

  if(step <= 5) then 
    dt = 2.d-1*dt
  endif
    
                
  ! Recalculates "dt" if output time "TIMEOUT" is exceeded
  if( (time + dt) > tend ) then
      dt = tend - time
  
  elseif( (time + dt) > timeout ) then
      dt = timeout - time
      
  endif  
  
  end subroutine tstep_HD


  subroutine tstep_RHD(timeout)
    use variables,    only  : prim, cs, vx1, vx2
    use mesh,         only  : dx, dy !,NDIMS
    use mesh,         only  : time, tend, step, dt
!==============================================================================! 
! Description:                                                                 !   
!    to apply the CFL condition to find a stable time step size DT             !
!    returns:  dt                                                              ! 
!                                                                              !
! Author: Victor Mourão Roque                                                  !
! last modification:                                                           !
! date: 13/07/2017 ; by Victor Mourao Roque                                    !                                                                               !
! Refs: A. Mignone, G. Bodo, S. Massaglia, T. Matsakos, O. Tesileanu,          !
!      C. Zanni, and A. Ferrari, arXiv.org astro-ph, 228 (2007).               !
!                                                                              !
!==============================================================================!

! Subroutine arguments
  real(dp),  intent(in)      :: timeout

!  local variables 
  real(dp) :: dt_l, dt_l_x, dt_l_y, smax1, smax2, smax3, smax_x, smax_y, ds

  smax1 = MAXVAL((prim(:,:,vx1) + cs(:,:))/(1.d0 + prim(:,:,vx1) * cs(:,:) ))
  smax2 = MAXVAL((prim(:,:,vx1) - cs(:,:))/(1.d0 - prim(:,:,vx1) * cs(:,:) ))
  smax3 = MAXVAL( prim(:,:,vx1) )
  smax_x = DMAX1(smax1,smax2,smax3)

  smax1 = MAXVAL((prim(:,:,vx2) + cs(:,:))/(1.d0 + prim(:,:,vx2) * cs(:,:) ))
  smax2 = MAXVAL((prim(:,:,vx2) - cs(:,:))/(1.d0 - prim(:,:,vx2) * cs(:,:) ))
  smax3 = MAXVAL( prim(:,:,vx2) )
  smax_y = DMAX1(smax1,smax2,smax3)
    
  ds = dx 

  dt_l_x = ds/smax_x

  ds = dy 

  dt_l_y = ds/smax_y

  dt_l = min(dt_l_x,dt_l_y)

  !compute time step dt
  dt = cfl*dt_l
  
   
  ! Recalculates "dt" for early times (dt is reduced to compensate 
  !for approximate calculation of smax)

  if(step <= 5) then 
    dt = 2.d-1*dt
  endif
    
                
  ! Recalculates "dt" if output time "TIMEOUT" is exceeded
  if( (time + dt) > tend ) then
      dt = tend - time
  
  elseif( (time + dt) > timeout ) then
      dt = timeout - time
      
  endif  
  
  end subroutine tstep_RHD


  
  

END MODULE evolution
