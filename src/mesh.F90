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

MODULE mesh
use constants, only : dp, i4b

IMPLICIT NONE

!==============================================================================!  
! Description:                                                                 !
!   Module responsable to generate the Mesh and coordenates                    !
!                                                                              !
! Author: Victor Mourão Roque                                                  !
! last modification:                                                           !
! date: 05/08/2017 ; by Victor Mourao Roque                                    ! 
!                                                                              !
!==============================================================================!


!=== MODULE VARIABLES ===
!
! time variables:
!
!   step  - integration step number
!   tstep - total number of time steps
!   time  - the evolution time in the code units
!   dts   - the current time step to snapshots
!   dti   - the initial time step
!   dt    - the new time step estimated from the CFL condition
!
  integer(i4b), save :: step    =  0
  integer(i4b), save :: tstep   =  0
  real(dp)    , save :: time    =  0.0d+00
  real(dp)    , save :: tini    =  0.0d+00
  real(dp)    , save :: tend    =  1.d+00
  real(dp)    , save :: dts     = -1.0d+00
  real(dp)    , save :: dti     =  1.0d-08
  real(dp)    , save :: dt      =  0.0d+00
  integer(i4b), save :: restart =  0

! domain dimensions:
!
!   icell - the resolution in the X direction
!   jcells - the resolution in the Y direction
!   kcells - the resolution in the Z direction
!   ng - the number of ghost cells
!   r_s - reconstruction order (2*r_s + 1) 
!
  integer(i4b), save :: NDIMS   = 2
  integer(i4b), save :: icell   = 1
  integer(i4b), save :: jcell   = 1
  integer(i4b), save :: kcell   = 1
  integer(i4b), save :: ng      = 1
  integer(i4b), save :: r_s     = 1

! domain bounds:
!
!   xmin, ymin, zmin - the lower coordinate bounds in the X, Y, and Z directions
!   xmax, ymax, zmax - the upper coordinate bounds in the X, Y, and Z directions
!   xlen, ylen, zlen - the domain sizes in the X, Y, and Z directions
!
  real(dp), save :: xmin  = 0.0d+00
  real(dp), save :: ymin  = 0.0d+00
  real(dp), save :: zmin  = 0.0d+00
  real(dp), save :: xmax  = 1.0d+00
  real(dp), save :: ymax  = 1.0d+00
  real(dp), save :: zmax  = 1.0d+00
  real(dp), save :: xlen  = 1.0d+00
  real(dp), save :: ylen  = 1.0d+00
  real(dp), save :: zlen  = 1.0d+00

! mesh dimensions:
!
!   im, jm, km          - the subdomain dimensions including the ghost zones
!   ibeg , jbeg , kbeg  - the first indices of domain (after the ghost zones)
!   iend , jend , kend  - the last indices of domain (before the ghost zones)
!
  integer(i4b), save :: im      = 1, jm      = 1, km      = 1
  integer(i4b), save :: ibeg    = 1, jbeg    = 1, kbeg    = 1
  integer(i4b), save :: iend    = 1, jend    = 1, kend    = 1

! coordinates increments
!
  real(dp), save :: dx    = 1.0d+00, dy    = 1.0d+00, dz    = 1.0d+00
  real(dp), save :: dxh   = 1.0d+00, dyh   = 1.0d+00, dzh   = 1.0d+00
  real(dp), save :: dxmin = 1.0d+00

! cell and face centered coordinate vectors
!
  real(dp), save, dimension(:), allocatable :: xc, yc, zc
  real(dp), save, dimension(:), allocatable :: xi, yi, zi

! the domain dimensions, its lengths, minima, maxima,
! and the bounds of the local subdomain
!
  integer(i4b), save, dimension(3) :: dims = 1
  integer(i4b), save, dimension(3) :: cells = 1
  integer(i4b), save, dimension(3) :: cbeg = 1
  integer(i4b), save, dimension(3) :: cend = 1
  real(dp), save, dimension(3) :: dlen = 1.0d+00
  real(dp), save, dimension(3) :: dmin = 0.0d+00
  real(dp), save, dimension(3) :: dmax = 1.0d+00
  !real(dp), save, dimension(3) :: ds   = 1.0d+00 ! size of cell in diff direct.


! Public scope by default
!
  PUBLIC

  private :: get_ghost_cell


CONTAINS

!==============================================================================!  
! Description:                                                                 !
!   Subroutine initializes mesh coordinates and other mesh parameters          !
!   none arguments needed                                                      !
!                                                                              !
! Author: Victor Mourão Roque                                                  !
! last modification:                                                           !
! date: dd/mm/yyyy ; by Victor Mourao Roque                                    ! 
!                                                                              !
!==============================================================================!	
!
  subroutine init_mesh()
    use parameters, only : get_parameter
  
  IMPLICIT NONE

! Local variables
  integer(i4b)       :: i, j, k 


!------------------------------------------------------------------------------!
! Time components
!
  !call get_parameter_real("dt", dt) => calculated by CFL condition
  call get_parameter("dti", dti) !=> calculeted in timestep subroutine
  call get_parameter("tend", tend)
  call get_parameter("dts", dts)

  call get_parameter("restart", restart)
  if (restart /= 0 ) then
    call get_parameter("tini", tini)
    write(*,*) "The simulation will restart from time:", tini
    time = tini
  else
    time = 0._dp
  endif
  
  if( dts <= 0 ) then
    write(*,*) "The value of the time step for output (dts) was not informed"
    write(*,*) "Using the number os time step to obtain dts:"
    call get_parameter("tstep", tstep)
    dts = (tend - time)/real(tstep,dp)
    write(*,*) "dts = ", dts
  endif


  


! grid components
!  
  call get_parameter("NDIMS", NDIMS)
  call get_parameter("icell", icell)
  call get_parameter("jcell", jcell)
#ifdef R3D
  call get_parameter("kcell", kcell)
#endif /* R3D */
  
  call get_parameter("xmin", xmin )
  call get_parameter("xmax", xmax )
  call get_parameter("ymin", ymin )
  call get_parameter("ymax", ymax )
#ifdef R3D
  call get_parameter("zmin", zmin )
  call get_parameter("zmax", zmax )
#endif /* R3D */

! initialize the domain dimensions array
!
  dims(1) = icell
  dims(2) = jcell
#ifdef R3D
  dims(3) = kcell
#endif /* R3D */

! initialize the ghost cell
!
  call get_ghost_cell()

! calculate dimensions and indices along the X direction
!
    if (icell > 1) then

      im  = icell + 2 * ng
      ibeg  =  1 + ng
      iend  = icell + ng
      cells(1) = im
      cbeg(1) = ibeg
      cend(1) = iend

    end if

! calculate dimensions and indices along the Y direction
!
    if (jcell > 1) then

      jm  = jcell + 2 * ng
      jbeg  =  1 + ng
      jend  = jcell + ng
      cells(2) = jm
      cbeg(2) = jbeg
      cend(2) = jend

    end if

! calculate dimensions and indices along the Z direction
!
    if (kcell > 1) then

      km  = kcell + 2 * ng
      kbeg  =  1 + ng
      kend  = kcell + ng
      cells(3) = km
      cbeg(3) = kbeg
      cend(3) = kend

    end if

! calculate the domain physical lengths
!
    xlen  = xmax - xmin
    ylen  = ymax - ymin
#ifdef R3D
    zlen  = zmax - zmin
#endif /* R3D */

! prepare the coordinate increments
!
    dx  = xlen / icell
 !   ds(1) = dx
    dxh = 0.5_dp*dx
    dy  = ylen / jcell
 !   ds(2) = dy
    dyh = 0.5_dp*dy
#ifdef R3D
    dz = zlen / kcell
    dzh = 0.5_dp*dz
!    ds(3) = dz
#endif /* R3D */

! initialize the domain lendths, minima and maxima
!
    dlen(1) = xlen
    dmin(1) = xmin
    dmax(1) = xmax

    dlen(2) = ylen
    dmin(2) = ymin
    dmax(2) = ymax
#ifdef R3D
    dlen(3) = zlen
    dmin(3) = zmin
    dmax(3) = zmax
#endif /* R3D */

! find the minimum coordinate interval
!
#ifdef R3D
    dxmin = min(dx, dy, dz)
#else /* R3D */
    dxmin = min(dx, dy)
#endif /* R3D */


! allocate the variables for coordinates
!
    allocate(xc(im))
    allocate(xi(im))
    allocate(yi(jm))
    allocate(yc(jm))
!#ifdef R3D    
    allocate(zc(km))
    allocate(zi(km))
!#endif /* R3D */

! generates coordinates
!
    if (icell > 1) then
      do concurrent (i = ibeg:im)
        xi(i) = xmin + REAL(i-ibeg)*dx
        xc(i) = xmin + REAL(i-ibeg)*dx + dxh 
      end do
! differents forms of assingn a array			
! 			do i = 1, im
! 				xi(i) = xmin + REAL(i-ibeg)*dx
! 				xc(i) = xmin + REAL(i-ibeg)*dx + dxh 
! 			enddo
!			 xi(:) = ((/(i, i=1,im)/)-ibeg)*dx + xmin
!			 xc(:) = xi(:) + dxh
    else
      xi(:) = 0._dp ; xc(:) = 0._dp
    endif

    if (jcell > 1) then
      do concurrent (j = jbeg:jm)
        yi(j) = ymin + REAL(j-jbeg)*dy
        yc(j) = ymin + REAL(j-jbeg)*dy + dyh 
      end do
    else
      yi(:) = 0._dp ; yc(:) = 0._dp
    endif

!#ifdef R3D    
    if (kcell > 1) then
      do concurrent (k = kbeg:km)
        zi(k) = zmin + REAL(k-kbeg)*dz
        zc(k) = zmin + REAL(k-kbeg)*dz + dzh 
      end do
    else
      zi(:) = 0._dp ; zc(:) = 0._dp
    endif
!#endif /* R3D */

END SUBROUTINE init_mesh


!==============================================================================!  
! Description:                                                                 !
!   Subroutine finalizes allocated mesh coordinates by deallocating them       !
!   none arguments needed                                                      !
!                                                                              !
! Author: Victor Mourão Roque                                                  !
! last modification:                                                           !
! date: dd/mm/yyyy ; by Victor Mourao Roque                                    ! 
!                                                                              !
!==============================================================================!
!
  SUBROUTINE finalize_mesh()

!-------------------------------------------------------------------------------

! deallocate coordinates
!
    if (allocated(xc)) deallocate(xc)
    if (allocated(yc)) deallocate(yc)
    if (allocated(zc)) deallocate(zc)
    if (allocated(xi)) deallocate(xi)
    if (allocated(yi)) deallocate(yi)
    if (allocated(zi)) deallocate(zi)


!-------------------------------------------------------------------------------
!
  END SUBROUTINE finalize_mesh


!==============================================================================!  
! Description:                                                                 !
!   Subroutine to determine the ghost cells based on the order of the          !
!   reconstruction method:                                                     !
! ng =                                                                         !
! 2 : to first and third order reconstruction WENO3                            !
! 3 : to five order reconstruction MP5, WENO5                                  !
! 4 : to five seventh reconstruction WENO7                                     !
! and the order parameter r_s needed to select the stencil used in             !
! reconstruction and split schemes                                             !
!   none arguments needed                                                      !
!                                                                              !
! Author: Victor Mourão Roque                                                  !
! last modification:                                                           !
! date: 19/07/2017 ; by Victor Mourao Roque                                    ! 
!                                                                              !
!==============================================================================!	

  subroutine get_ghost_cell()
    use parameters, only : get_parameter
  implicit none
    
    ! local variables
    !
    character(len=255) :: recons_method = "WENO5Z"

    call get_parameter("recons_method", recons_method)

    select case(trim(recons_method))
    case("NOREC")
      ng  = 2
      r_s = 0
    case("WENO3", "WENO3p")
      ng  = 2
      r_s = 1
    case("WENO5", "WENO5Z", "WENO5v", "WENO5Zv", "MP5")
      ng  = 3
      r_s = 2
    case("WENO7", "WENO7Zv")
      ng  = 4
      r_s = 3
    case default
    ! If the user did not choose appropriately, WENO5z is chosen by default
      ng  = 3
      r_s = 2
    end select

  END SUBROUTINE get_ghost_cell

!   subroutine print_mesh()
  
!     integer(i4b) :: i

!     print*,dx,dxh
!     print*,    

!     do i = 1,size(xi)
! 	    print*,xi(i)
!       print*,i,xc(i)
!     enddo
!   end subroutine print_mesh

end module mesh
