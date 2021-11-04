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


MODULE reconstruction
  use constants,   only  : dp, i4b

  IMPLICIT NONE

!==============================================================================!  
! Description:                                                                 !
!   < Say what this module contains >                                          !
!                                                                              !
! Author: Victor Mourão Roque                                                  !
! last modification:                                                           !
! date: 25/07/2017 ; by Victor Mourao Roque                                    ! 
!                                                                              !
! Refs:                                                                        !
!                                                                              !
!                                                                              !
!==============================================================================!


  PRIVATE  ! Private scope by default


! Module constants
!  < Parameter declarations >


! Local Variables
  procedure(rec_weno3), pointer, save :: rec_cells => null()

! Public members
  PUBLIC :: rec_cells, initialize_reconstr, finalize_reconstr


CONTAINS
  
subroutine initialize_reconstr()
  use parameters    , only : get_parameter
  implicit none
! local variables
!
    character(len=255) :: reconstr_method = "norec"    


! select the integration method
 call get_parameter("reconstr_method", reconstr_method)

  select case(trim(reconstr_method))

  case('norec','NOREC')

    rec_cells => norec

  case('weno3','WENO3')

    rec_cells => rec_weno3

  case('weno3p','WENO3p')

    rec_cells => rec_weno3p

  case('weno5','WENO5')

    rec_cells => rec_weno5

  case('weno5v','WENO5v')

    rec_cells => rec_weno5v

  case('weno5z','WENO5z')

    rec_cells => rec_weno5z

  case default

    write (*,"(1x,a)") "The selected time advance method is not " //       &
                           "implemented: " // trim(reconstr_method)
    rec_cells => rec_weno5

  end select

  end subroutine initialize_reconstr

  subroutine finalize_reconstr()
    implicit none
!
! release the procedure pointers
!
    nullify(rec_cells)
  
  end subroutine finalize_reconstr


  subroutine norec(coord,sflux,sf)
    use variables,  only : nc
    use mesh,       only : cells, r_s
    implicit none
! Input Variables
  integer(i4b), intent(in)  :: coord
  real(dp), dimension(cells(coord),nc,-r_s:r_s), intent(in)  :: sflux
! sf => splited flux reconstructed  	
!
  real(dp), dimension(cells(coord),nc), intent(out)          :: sf
 
  sf(:,:) =  sflux(:,:,0)  

 end subroutine norec


!==============================================================================!     
! A 3th order WENO                                                             !
! i => i+1/2                                                                   !
! last modification:                                                           !
! date: 27/06/2014 ; by Victor Mourao Roque                                    !
! based on: [1] A. Mignone, P. Tzeferacos, and G. Bodo,                        !
!               Journal of Computational Physics 229, 5896–5920 (2010).        ! 
! j = i+1/2 (interface)                                                        !
!==============================================================================! 

 subroutine rec_weno3(coord,sflux,sf)
  use variables,  only : nc
  use mesh,       only : cells, r_s, cend, cbeg
implicit none
! Input Variables
  integer(i4b), intent(in)  :: coord
  real(dp), dimension(cells(coord),nc,-r_s:r_s), intent(in)  :: sflux
! sf => splited flux reconstructed  	
!
  real(dp), dimension(cells(coord),nc), intent(out)          :: sf
!
! Local Variables
 real(dp) :: beta1, beta0
 real(dp) :: omega0, omega1 
 real(dp) :: a0, a1, sum_al

 real(dp), parameter :: eps_tol = 1.d-10
 real(dp), parameter :: d1 = 1.d0/3.d0
 real(dp), parameter :: d0 = 2.d0/3.d0
 integer(i4b)        :: i, j
 
!  do concurrent (i = 1:cells(coord))
do i = cbeg(coord)-1,cend(coord)   !cells(coord)
    DO j = 1, nc
    
! equation (25) of [1]
    beta0 = (sflux(i,j,1) - sflux(i,j,0))**2
    beta1 = (sflux(i,j,0) - sflux(i,j,-1))**2

! equation (28) of [1]  
    a0 = d0/(beta0 + eps_tol)**2
    a1 = d1/(beta1 + eps_tol)**2
    sum_al = a1+a0
  
    omega0 = a0/(sum_al)
    omega1 = a1/(sum_al)
 
! equation (27) of [1]
    sf(i,j) = .5d0*omega0*(sflux(i,j,0)+sflux(i,j,1)) + &
              .5d0*omega1*(-sflux(i,j,-1)+3.d0*sflux(i,j,0))
    ENDDO

  ENDDO

end subroutine rec_weno3



subroutine rec_weno3p(coord,sflux,sf)
  use variables,  only : nc
  use mesh,       only : cells, r_s, cend, cbeg
implicit none
! Input Variables
  integer(i4b), intent(in)  :: coord
  real(dp), dimension(cells(coord),nc,-r_s:r_s), intent(in)  :: sflux
! sf => splited flux reconstructed    
!
  real(dp), dimension(cells(coord),nc), intent(out)          :: sf
!
! Local Variables
  real(dp) :: beta1, beta0, difplus, difminus, dif
  real(dp) :: omega0, omega1, a0, a1, sum_al

  real(dp), parameter :: eps_tol = 1.d-10
  real(dp), parameter :: d1 = 1.d0/3.d0
  real(dp), parameter :: d0 = 2.d0/3.d0
  integer(i4b)        :: i, j
 
!  do concurrent (i = 1:cells(coord))
do i = cbeg(coord)-1,cend(coord)   !cells(coord)
  do j = 1, nc

  ! equation (25) of [1]
    difplus = sflux(i,j,1) - sflux(i,j,0) 
    difminus = sflux(i,j,0) - sflux(i,j,-1)
    beta0 = difplus*difplus
    beta1 = difminus*difminus

  ! equation (28) of [1]  
    dif = DABS(difplus - difminus)
    a0 = d0*(1.d0 + dif**2/(beta0 + eps_tol))
    a1 = d1*(1.d0 + dif**2/(beta1 + eps_tol))
    sum_al = a0 + a1

    omega0 = a0/(sum_al)
    omega1 = a1/(sum_al)

  ! equation (27) of [1]
    sf(i,j) = .5d0*omega0*( sflux(i,j,0)  + 1.d0*sflux(i,j,1)) + &
              .5d0*omega1*(-sflux(i,j,-1) + 3.d0*sflux(i,j,0))    

  enddo
enddo

end subroutine rec_weno3p
!==============================================================================!     
! A 5th order WENO                                                             !
! i => i+1/2                                                                   !
! last modification:                                                           !
! date: 27/06/2014 ; by Victor Mourao Roque                                    !
! based on: [1] A. Mignone, P. Tzeferacos, and G. Bodo,                        !
!               Journal of Computational Physics 229, 5896–5920 (2010).        ! 
!           [2] Z. He, X. Li, D. Fu, and Y. Ma,                                !
!               Sci. China Phys. Mech. Astron. 54, 511–522 (2011).             !
!           [3] A. Suresh and H. T. Huynh,                                     !
!               Journal of Computational Physics 136, 83–99 (1997).            !
! j = i+1/2 (interface)                                                        !
!                                                                              !
!            ALGUM FATOR ERRADO NESSA IMPLEMENTAÇÃO!!!!!!!!!!!                 !
!                                                                              !          
!==============================================================================!


subroutine rec_weno5(coord,sflux,sf)
  use variables,  only : nc
  use mesh,       only : r_s, cend, cells, cbeg
implicit none
! Input Variables
  integer(i4b), intent(in)  :: coord
  real(dp), dimension(cells(coord),nc,-r_s:r_s), intent(in)  :: sflux
! sf => splited flux reconstructed  	
!
  real(dp), dimension(cells(coord),nc), intent(out)          :: sf
!
! Local Variables
 real(dp) :: beta0, beta1, beta2
 real(dp) :: omega0, omega1, omega2 
 real(dp) :: a0, a1, a2, sum_al

 real(dp), parameter :: eps_tol = 1.d-10
 real(dp), parameter :: d2 = 0.3d0 !3.d0/10.d0
 real(dp), parameter :: d1 = 3.d0/5.d0
 real(dp), parameter :: d0 = 0.1d0 !1.d0/10.d0
 real(dp), parameter :: b1 = 13.d0/12.d0
 real(dp), parameter :: b2 = .25d0
 
 real(dp), parameter :: p = 2.d0
 integer(i4b)        :: i, j

 
  
  !DO i = 0, cells(coord)
  !DO CONCURRENT (i =1:cells(coord))
  DO i = cbeg(coord)-1, cend(coord) ! cells(coord)
    do j = 1, nc

    
! equation (35) of [1]
  beta0 = b1*(sflux(i,j,0)-2.d0*sflux(i,j,-1)+sflux(i,j,-2))**2 &
            + b2*(3.d0*sflux(i,j,0)-4.d0*sflux(i,j,-1)+sflux(i,j,-2))**2 

  beta1 = b1*(sflux(i,j,1)-2.d0*sflux(i,j,0)+sflux(i,j,-1))**2 &
            + b2*(sflux(i,j,1) - sflux(i,j,-1))**2 

  beta2 = b1*(sflux(i,j,2)-2.d0*sflux(i,j,1)+sflux(i,j,0))**2 &
            + b2*(3.d0*sflux(i,j,0)-4.d0*sflux(i,j,1)+sflux(i,j,2))**2 

  

! equation (34) of [1]  to WENO5
    a0 = d0/(beta0 + eps_tol)**p
    a1 = d1/(beta1 + eps_tol)**p
    a2 = d2/(beta2 + eps_tol)**p

  
    sum_al = 1._dp/(a0 + a1 + a2)
    omega0 = a0*sum_al
    omega1 = a1*sum_al
    omega2 = a2*sum_al
 

! equation (33) of [1]
    sf(i,j) = omega0 * &
        (2.d0*sflux(i,j,-2) - 7.d0*sflux(i,j,-1) + 11.d0*sflux(i,j,0)) &
              + omega1 * &
        (-1.d0*sflux(i,j,-1) + 5.d0*sflux(i,j,0) + 2.d0*sflux(i,j,1)) &
              + omega2 * &
        (2.d0*sflux(i,j,0) + 5.d0*sflux(i,j,1) - 1.d0*sflux(i,j,2))
     sf(i,j) = sf(i,j)/6.d0    
  enddo
  ENDDO


end subroutine rec_weno5


subroutine rec_weno5v(coord,sflux,sf)
  use variables,  only : nc
  use mesh,       only : r_s, cells
implicit none
! Input Variables
  integer(i4b), intent(in)  :: coord
  real(dp), dimension(cells(coord),nc,-r_s:r_s), intent(in)  :: sflux
! sf => splited flux reconstructed    
!
  real(dp), dimension(cells(coord),nc), intent(out)          :: sf
!
! Local Variables
 real(dp), dimension(cells(coord),nc) :: beta0, beta1, beta2
 real(dp), dimension(cells(coord),nc) :: omega0, omega1, omega2 
 real(dp), dimension(cells(coord),nc) :: a0, a1, a2, sum_al

 real(dp), parameter :: eps_tol = 1.d-10
 real(dp), parameter :: d2 = 0.3d0 !3.d0/10.d0
 real(dp), parameter :: d1 = 3.d0/5.d0
 real(dp), parameter :: d0 = 0.1d0 !1.d0/10.d0
 real(dp), parameter :: b1 = 13.d0/12.d0
 real(dp), parameter :: b2 = .25d0
 
 real(dp), parameter :: p = 2.d0
 

  ! equation (35) of [1]
  beta0(:,:) = b1*(sflux(:,:,0)-2.d0*sflux(:,:,-1)+sflux(:,:,-2))**2 &
            + b2*(3.d0*sflux(:,:,0)-4.d0*sflux(:,:,-1)+sflux(:,:,-2))**2 

  beta1(:,:) = b1*(sflux(:,:,1)-2.d0*sflux(:,:,0)+sflux(:,:,-1))**2 &
            + b2*(sflux(:,:,1) - sflux(:,:,-1))**2 

  beta2(:,:) = b1*(sflux(:,:,2)-2.d0*sflux(:,:,1)+sflux(:,:,0))**2 &
            + b2*(3.d0*sflux(:,:,0)-4.d0*sflux(:,:,1)+sflux(:,:,2))**2 

    

! equation (34) of [1]  to WENO5
    a0(:,:) = d0/(beta0(:,:) + eps_tol)**p
    a1(:,:) = d1/(beta1(:,:) + eps_tol)**p
    a2(:,:) = d2/(beta2(:,:) + eps_tol)**p

  
    sum_al = a0 + a1 + a2
    omega0 = a0/(sum_al)
    omega1 = a1/(sum_al)
    omega2 = a2/(sum_al)
 

! equation (33) of [1]
    sf(:,:) = omega0(:,:) * &
        (2.d0*sflux(:,:,-2) - 7.d0*sflux(:,:,-1) + 11.d0*sflux(:,:,0)) &
              + omega1(:,:) * &
        (-1.d0*sflux(:,:,-1) + 5.d0*sflux(:,:,0) + 2.d0*sflux(:,:,1)) &
              + omega2(:,:) * &
        (2.d0*sflux(:,:,0) + 5.d0*sflux(:,:,1) - 1.d0*sflux(:,:,2))



  sf = sf/6.d0


end subroutine rec_weno5v




!==============================================================================!     
! A 5th order WENO reconstruction for finite diference scheme                  !
! i => i+1/2                                                                   !
! last modification:                                                           !
! date: 01/10/2012 ; by Victor Mourao Roque                                    !
! based on: [1] A. Mignone, P. Tzeferacos, and G. Bodo,                        !
!               Journal of Computational Physics 229, 5896–5920 (2010).        ! 
!           [2] R. Borges, M. Carmona, B. Costa, and W. S. Don,                !
!               Journal of Computational Physics 227, 3191 (2008).             !
!           [3] M. Castro, B. Costa, and W. S. Don,                            !
!               Journal of Computational Physics 230, 1766 (2011).             !
!               Sci. China Phys. Mech. Astron. 54, 511–522 (2011).             !
!                                                                              !   
! j = i+1/2 (interface)                                                        !
! r => (2r-1) is related with the order of the method                          !
!==============================================================================!


subroutine rec_weno5z(coord,sflux,sf)
use variables,  only : nc
use mesh,       only : r_s, cells, cbeg
use constants,  only : i1b
implicit none
! Input Variables
  integer(i4b), intent(in)  :: coord
  real(dp), dimension(cells(coord),nc,-r_s:r_s), intent(in)  :: sflux
! sf => splited flux reconstructed    
!
  real(dp), dimension(cells(coord),nc), intent(out)          :: sf
!
! Local Variables
 integer(i1b), parameter :: r = 3 ! order (2r-1)
 real(dp), dimension(0:r-1) :: beta
 real(dp) :: omega0, omega1, omega2, tau_r 
 real(dp) :: a0, a1, a2, sum_al
 real(dp), parameter :: eps_tol = 1.d-40
 real(dp), parameter :: d2 = 0.3d0 !3.d0/10.d0
 real(dp), parameter :: d1 = 3.d0/5.d0
 real(dp), parameter :: d0 = 0.1d0 !1.d0/10.d0
 real(dp), parameter :: b1 = 13.d0/12.d0
 real(dp), parameter :: b2 = .25d0
 
 real(dp), parameter :: p = real(r,dp)
 integer(i4b)        :: i, j
 

!$OMP parallel do private(i,j,k, beta0, beta1, beta2, omega0, omega1, omega2, &
!$OMP&   a0, a1, a2, sum_al) shared(sf,sflux)
  
  !DO i = 0, cells(coord)
!   DO CONCURRENT (i = cbeg(coord):cells(coord))
DO i = cbeg(coord)-1,cells(coord)
    DO j=1,nc
        
! equation (35) of [1]
  !beta0
  beta(0) = b1*(sflux(i,j,0)-2.d0*sflux(i,j,-1)+sflux(i,j,-2))**2 &
            + b2*(3.d0*sflux(i,j,0)-4.d0*sflux(i,j,-1)+sflux(i,j,-2))**2 

  !beta1 
  beta(1)= b1*(sflux(i,j,1)-2.d0*sflux(i,j,0)+sflux(i,j,-1))**2 &
            + b2*(sflux(i,j,1) - sflux(i,j,-1))**2 

  !beta2
  beta(2) = b1*(sflux(i,j,2)-2.d0*sflux(i,j,1)+sflux(i,j,0))**2 &
            + b2*(3.d0*sflux(i,j,0)-4.d0*sflux(i,j,1)+sflux(i,j,2))**2 

  ! equation (28) of [3] WENOZ7
!     IF (MOD(r,2) .EQ. 1) THEN
!       tau_r = DABS(beta(0) - beta(r-1))
!     ELSE IF(MOD(r,2) .EQ. 0) THEN
      tau_r = DABS(beta(0) - beta(1) - beta(r-2) + beta(r-1))
!     ELSE
!       WRITE(*,*) "Problem in tau_r of WENOZ reconstruction"
!       !STOP "Problem in tau_r of WENOZ reconstruction"
!     END IF

! equation (34) of [1]  to WENO - Z
  a0 = d0*(1.d0 + (tau_r/(beta(0) + eps_tol))**p)
  a1 = d1*(1.d0 + (tau_r/(beta(1) + eps_tol))**p)
  a2 = d2*(1.d0 + (tau_r/(beta(2) + eps_tol))**p)
  
  
    sum_al = 1._dp/(a0 + a1 + a2)
    omega0 = a0*sum_al
    omega1 = a1*sum_al
    omega2 = a2*sum_al
 

! equation (33) of [1]
    sf(i,j) = omega0 * &
        (2.d0*sflux(i,j,-2) - 7.d0*sflux(i,j,-1) + 11.d0*sflux(i,j,0)) &
              + omega1 * &
        (-1.d0*sflux(i,j,-1) + 5.d0*sflux(i,j,0) + 2.d0*sflux(i,j,1)) &
              + omega2 * &
        (2.d0*sflux(i,j,0) + 5.d0*sflux(i,j,1) - 1.d0*sflux(i,j,2))

    sf(i,j) = sf(i,j)/6.d0 

    ENDDO

  ENDDO

!$OMP end do parallel

!   sf = sf/6.d0
  
 end subroutine  rec_weno5z

END MODULE reconstruction
