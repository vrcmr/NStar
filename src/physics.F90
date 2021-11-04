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

module physics
  use constants,  only : dp, i4b
  use mesh,       only  : im, jm, cells

IMPLICIT NONE

!==============================================================================!  
! Description:                                                                 !
!   This module determines the physical equations that will be used on the     !
! model, HD or RHD                                                             !
!                                                                              !
! Author: Victor Mourão Roque                                                  !
! last modification:                                                           !
! date: 18/07/2017 ; by Victor Mourao Roque                                    ! 
!                                                                              !
! Refs:                                                                        !
! [1] D. S. Balsara and C.-W. Shu,                                             !
!           Journal of Computational Physics 160, 405–452 (2000).              !
! [2] D. Radice and L. Rezzolla, arXiv astro-ph.IM, (2012).                    !
! [3] A. Mignone, P. Tzeferacos, and G. Bodo,                                  !
!           Journal of Computational Physics 229, 5896–5920 (2010).            !
!                                                                              !
! Melhorias:                                                                   !
!   utilizar submodules para descrever as equações                             !
!                                                                              !
!==============================================================================!

PRIVATE

! Local Variables
  procedure(p2c_RHD),       pointer, save :: prim2cons => null()
  procedure(c2p_RHD),       pointer, save :: cons2prim => null()
  procedure(phys_flux_RHD), pointer, save :: phys_flux => null()
  procedure(eigenval_RHD),  pointer, save :: eigenval  => null()
  procedure(eigenvec_RHD),  pointer, save :: eigenvec  => null()
 
PUBLIC :: initialize_physics, test_vectorsHD
PUBLIC :: prim2cons, cons2prim, phys_flux, eigenval, eigenvec

CONTAINS

subroutine initialize_physics()
  use parameters,  only  :  get_parameter


implicit none

  character(len=32)  ::  phys_eq = "RHD"


  call get_parameter("phys_eq",phys_eq)

  select case(trim(phys_eq))
    case("HD","hd","hydro")

      prim2cons => p2c_HD
      cons2prim => c2p_HD
      phys_flux => phys_flux_HD
      eigenval  => eigenval_HD
      eigenvec  => eigenvec_HD

    case("RHD","rhd")
      
      prim2cons => p2c_RHD
      cons2prim => c2p_RHD
      phys_flux => phys_flux_RHD
      eigenval  => eigenval_RHD
      eigenvec  => eigenvec_RHD

    case default
      
!       prim2cons => p2c_RHD
!       cons2prim => c2p_RHD
!       pflux     => pflux_RHD
!       eigenval  => eigenval_RHD

  end select

end subroutine initialize_physics


!==============================================================================!
!                                                                              !     
!                           NEWTONIAN HYDRODYNAMICS                            !
!                                                                              !
!==============================================================================! 


!==============================================================================!     
! calculates conserved variables from the primitive variables                  !
! last modification:                                                           !
! date: 11/06/2014 ; by Victor Mourao Roque                                    !
!==============================================================================! 

subroutine p2c_HD()
  use variables,  only  : prim, rho, vx1, vx2, i_e
  use variables,  only  : cons, D, Sx1, Sx2, tau
  implicit none
! Local Variables
  INTEGER(I4B) ::  i,j
  REAL(DP) :: V2
  

!$OMP PARALLEL PRIVATE(i,j,V2) &
!$OMP SHARED(cons)
!$OMP DO
  do j = 1, jm  
    do i= 1, im
      cons(i,j,D)   = prim(i,j,rho)
      cons(i,j,Sx1) = prim(i,j,rho)*prim(i,j,vx1)
      cons(i,j,Sx2) = prim(i,j,rho)*prim(i,j,vx2)
     
      V2 = prim(i,j,vx1)**2 + prim(i,j,vx2)**2

      cons(i,j,tau) = prim(i,j,rho) * &
                     ( prim(i,j,i_e)  + .5d0 * V2 )
      ! cons(i,j,tau) = prim(i,j,rho) * h - prim(i,j,prs)

    enddo
  enddo
!$OMP END DO NOWAIT
!$OMP END PARALLEL 


end subroutine p2c_HD

!==============================================================================!     
! calculates the primitive variables from conservative variables               !
! last modification:                                                           !
! date: 11/06/2014 ; by Victor Mourao Roque                                    !
!==============================================================================!


subroutine c2p_HD()
  use constants,  only  : i4b, dp
  use variables,  only  : prim, rho, vx1, vx2, prs, i_e, k_t, cs
  use variables,  only  : cons, D, Sx1, Sx2, tau
  use mesh,       only  : iend, ibeg, jbeg, jend
  use EoS,        only  : sound, pressure, k_tilde
  implicit none
! Local Variables
  integer(i4b)  :: i,j
  real(dp)      :: V2, int_energ

  do j= jbeg, jend
    do i= ibeg, iend

      prim(i,j,rho) = cons(i,j,D)
      prim(i,j,vx1) = cons(i,j,Sx1)/cons(i,j,D)
      prim(i,j,vx2) = cons(i,j,Sx2)/cons(i,j,D)
      

      V2 = prim(i,j,vx1)**2 + prim(i,j,vx2)**2
      int_energ = cons(i,j,tau)/cons(i,j,D) - .5_dp*V2
      
      prim(i,j,prs) = pressure( prim(i,j,rho), int_energ)

      cs(i,j) = sound( prim(i,j,rho), prim(i,j,prs) )

      prim(i,j,i_e) = int_energ
      prim(i,j,k_t) = k_tilde(prim(i,j,rho), int_energ)


    enddo
  enddo

end subroutine c2p_HD


!==============================================================================!     
! for given values of the primitive variables "d,u,p" it calculates the        !
! Physical Flux vector [mass, momentum and energy] for relativistic            !
! one-dimensional Euler equations                                              !
! last modification:                                                           !
! date: 11/06/2014 ; by Victor Mourao Roque                                    !
!==============================================================================!


subroutine phys_flux_HD(coord,prim,pflux)
  use variables,  only : rho, vx1, vx2, prs, i_e, nc, np
  implicit none
  integer(i4b),             intent(in)  :: coord
  real(dp), dimension(cells(coord),np), intent(in)  :: prim
  real(dp), dimension(cells(coord),nc), intent(out) :: pflux
! Local Variables
  real(dp) :: E, V2
  integer(i4b) :: i, j

  select case(coord)

  case(1)

    do i = 1, im
      ! eq. (23) of [2]
      pflux(i,1) = prim(i,rho) * prim(i,vx1)
      pflux(i,2) = prim(i,rho) * prim(i,vx1)*prim(i,vx1) + prim(i,prs)
      pflux(i,3) = prim(i,rho) * prim(i,vx1)*prim(i,vx2)
      

      V2 =  prim(i,vx1)**2 + prim(i,vx2)**2
      E  = prim(i,rho) * (.5d0*V2 + prim(i,i_e) )

      pflux(i,4) = prim(i,vx1)*(E + prim(i,prs))

    enddo

  case(2)

    do j = 1, jm
      ! eq. (23) of [2]
      pflux(j,1) = prim(j,rho) * prim(j,vx2)
      pflux(j,2) = prim(j,rho) * prim(j,vx1)*prim(j,vx2)
      pflux(j,3) = prim(j,rho) * prim(j,vx2)*prim(j,vx2) + prim(j,prs)

      V2 =  prim(j,vx1)**2 + prim(j,vx2)**2
      E = prim(j,rho) * (.5d0*V2 + prim(j,i_e) )

      pflux(j,4) = prim(j,vx2)*(E + prim(j,prs))

    enddo

  case default

    write(*,*) 'Somethig wrong with argument coord'
    stop

  end select

end subroutine phys_flux_HD


!==============================================================================!     
! Subroutine to calculate the eigenvalues in interfaces of the cells           ! 
! j => i+1/2                                                                   !
! last modification:                                                           !
! date: 11/06/2014 ; by Victor Mourao Roque                                    !
! based on: [1] R. Donat, Journal of Computational Physics 146, 58–81 (1998).  !
! Notes:                                                                       !
!==============================================================================!
! verificar se na hora de chamar essa subrotina pode somente chamar uma 
! velocidade vx1 ou vx2
subroutine eigenval_HD(coord,prim,cs,lambda)
  use variables,  only : vx1, vx2, np
  implicit none
  integer(i4b), intent(in)  :: coord
  real(dp)    , dimension(cells(coord),np), intent(in)  :: prim
  real(dp)    , dimension(cells(coord)),    intent(in)  :: cs
  real(dp)    , dimension(cells(coord),np), intent(out) :: lambda
! Local Variables
  integer(i4b) :: i,j

! eigenvalues eq.(14) and (15) of [1]
! acoustic waves (+/-) and material waves (0)

SELECT CASE(coord)

 CASE(1)

  DO CONCURRENT(i=1:im) 
    lambda(i,4) = prim(i,vx1) + cs(i)
    lambda(i,3) = prim(i,vx1)
    lambda(i,2) = prim(i,vx1)
    lambda(i,1) = prim(i,vx1) - cs(i)
  END DO

 CASE(2)

  DO CONCURRENT(j=1:jm)
    lambda(j,4) = prim(j,vx2) + cs(j)
    lambda(j,3) = prim(j,vx2)
    lambda(j,2) = prim(j,vx2)
    lambda(j,1) = prim(j,vx2) - cs(j)
  END DO

  case default

    write(*,*) 'Somethig wrong with argument coord'
    stop

 END SELECT

end subroutine eigenval_HD


!==============================================================================!     
! Subroutine to calculate the eigenvectors (+,0,-) in interfaces of the cells  ! 
! i => i+1/2                                                                   !
! last modification:                                                           !
! date: 12/06/2014 ; by Victor Mourao Roque                                    !
! based on: [1] A. Suresh and H. T. Huynh, JCP 136, 83 (1997).                 !
! rev(:,:,1) and lev(:,:,1) are the positive eigenvectors                      !
! rev(:,:,2) and lev(:,:,2) are the zero eigenvectors                          !
! rev(:,:,3) and lev(:,:,3) are the negative eigenvectors                      !
!                                                                              !
! subrotina com autovetores para uma equação de estado arbitrária              !
! na equação de estado deve-se descrever a função gamma como                   !
! gamma = 1 + 1/d{Del p/Del epsilon} e a velocidade do som como                !
! cs = SQRT[{Del p/Del d} + p/d^2{Del p/Del epsilon}] =                        !
!    = SQRT[{Del p/Del d} + p*(gamma -1)/d]                                    !
!                                                                              !
!==============================================================================!


subroutine eigenvec_HD(coord, prim, cs, lambda, lev, rev)
  use variables,  only : rho, vx1, vx2, prs, i_e, k_t, nc, np
  implicit none
  integer(i4b), intent(in)  :: coord
  real(dp), dimension(cells(coord),np),    intent(in)  :: prim
  real(dp), dimension(cells(coord)),       intent(in)  :: cs
  real(dp), dimension(cells(coord),nc),    intent(in)  :: lambda
  real(dp), dimension(cells(coord),nc,nc), intent(out) :: lev, rev
! Local Variables
  integer(i4b) :: i
  real(dp) :: Energy_dens, h,theta,b, fac_l, q2, cs2


 SELECT CASE(coord)

 CASE(1)

 DO i = 1,cells(coord)
    q2 = prim(i,vx1)**2 + prim(i,vx2)**2
    cs2 = cs(i)**2
    Energy_dens = prim(i,i_e) + .5_dp*q2
    h = Energy_dens + prim(i,prs)/prim(i,rho)
    b = prim(i,k_t)
    theta = q2 - h + cs2/b 
    fac_l = .5d0*b/cs2


   rev(i,:,1) = (/ 1.d0, lambda(i,1), prim(i,vx2), h - prim(i,vx1)*cs(i) /)

   rev(i,:,2) = (/ 0.d0, 0.d0, 1.d0, prim(i,vx2) /)

   rev(i,:,3) = (/ 1.d0, lambda(i,2), prim(i,vx2), h - cs2/b /)

   rev(i,:,4) = (/ 1.d0, lambda(i,4), prim(i,vx2), h + prim(i,vx1)*cs(i) /)


   lev(i,:,1) = fac_l* (/ theta + prim(i,vx1)*cs(i)/b, -prim(i,vx1) - cs(i)/b ,&
                        -prim(i,vx2), 1.d0/)

   lev(i,:,2) = fac_l*(/-2.d0*prim(i,vx2)*cs2/b , 0.d0, 2.d0*cs2/b, 0.d0/)

   lev(i,:,3) = fac_l* (/2.d0*h - 2.d0*q2, 2.d0*prim(i,vx1), &
                                           2.d0*prim(i,vx2), -2.d0/)

   lev(i,:,4) = fac_l*(/theta - prim(i,vx1)*cs(i)/b, - prim(i,vx1) + cs(i)/b , &
                        -prim(i,vx2), 1.d0/)


 ENDDO

 CASE(2)
  DO i = 1,cells(coord)

    q2 = prim(i,vx1)**2 + prim(i,vx2)**2
    cs2 = cs(i)**2
    Energy_dens = prim(i,i_e) + .5d0*q2
    h = Energy_dens + prim(i,prs)/prim(i,rho)
    b = prim(i,k_t)
    theta = q2 - h + cs2/b 
    fac_l = .5d0*b/cs2


   rev(i,:,1) = (/ 1.d0, prim(i,vx1), lambda(i,1), h - prim(i,vx2)*cs(i) /)

   rev(i,:,2) = (/ 0.d0, 1.d0, 0.d0, prim(i,vx1) /)

   rev(i,:,3) = (/ 1.d0, prim(i,vx1), lambda(i,2), h - cs2/b /)

   rev(i,:,4) = (/ 1.d0, prim(i,vx1), lambda(i,4), h + prim(i,vx2)*cs(i) /)


   lev(i,:,1) = fac_l* (/ theta + prim(i,vx2)*cs(i)/b, -prim(i,vx1), &
                         - prim(i,vx2) - cs(i)/b, 1.d0/)

   lev(i,:,2) = fac_l*(/-2.d0*prim(i,vx1)*cs2/b , 2.d0*cs2/b, 0.d0, 0.d0/)

   lev(i,:,3) = fac_l* (/2.d0*h - 2.d0*q2, 2.d0*prim(i,vx1), &
                                           2.d0*prim(i,vx2), -2.d0/)

   lev(i,:,4) = fac_l*(/theta - prim(i,vx2)*cs(i)/b, -prim(i,vx1), &
                        - prim(i,vx2) + cs(i)/b, 1.d0/)

  ENDDO

  case default

    write(*,*) 'Somethig wrong with argument coord'
    stop

 END SELECT


end subroutine eigenvec_HD








!==============================================================================!
!                                                                              !     
!                         RELATIVISTIC HYDRODYNAMICS                           !
!                                                                              !
!==============================================================================! 


!==============================================================================!     
! calculates conserved variables from the primitive variables                  !
! last modification:                                                           !
! date: 11/06/2014 ; by Victor Mourao Roque                                    !
!==============================================================================! 

subroutine p2c_RHD() 
  use variables,  only  : prim, rho, vx1, vx2, prs, i_e
  use variables,  only  : cons, D, Sx1, Sx2, tau
  use mesh,       only : cells, im, im
  implicit none
! Local Variables
  integer(i4b) ::  i, j
  real(dp) :: W, h, ss
  

!$OMP PARALLEL PRIVATE(i,j,ss,h,W) &
!$OMP SHARED(cons,prim)
!$OMP DO
  do j = 1, jm
    do i= 1, im

      W = 1.d0/SQRT(1._dp - prim(i,j,vx1)**2 - prim(i,j,vx2)**2)
      h = 1.d0 + prim(i,j,i_e) + prim(i,j,prs)/prim(i,j,rho)
      !
      cons(i,j,D) = prim(i,j,rho)*W
      !
      ss = cons(i,j,D)*h*W
      !
      cons(i,j,Sx1) = ss*prim(i,j,vx1)
      cons(i,j,Sx2) = ss*prim(i,j,vx2)
      cons(i,j,tau) = ss - prim(i,j,prs) - cons(i,j,D)
    enddo
  enddo
!$OMP END DO NOWAIT
!$OMP END PARALLEL 

end subroutine p2c_RHD

! !==============================================================================!     
! ! calculates the primitive variables from conservative variables               !
! ! last modification:                                                           !
! ! date: 11/06/2014 ; by Victor Mourao Roque                                    !
! !==============================================================================!


subroutine c2p_RHD()
  use minpackmod1
  use constants,  only  : i4b, dp
  use variables,  only  : prim, rho, vx1, vx2, prs, i_e, k_t, cs
  use variables,  only  : cons, D, Sx1, Sx2, tau
  use mesh,       only  : iend, ibeg, jend, jbeg, ng
  use EoS,        only  : sound, k_tilde_p, eps
  implicit none
! Local Variables
  integer(i4b)  :: i, j
 
  REAL(DP), DIMENSION(im,jm) :: p_in   
  REAL(DP) :: W

! Variables used in hybrd1
  INTEGER(I4B) :: iflag
  INTEGER(I4B), parameter :: n = 1
  REAL(DP), DIMENSION(n) :: x, fvec
  LOGICAL :: check 
  REAL(KIND=8) :: ddd, ssu, ssv, ttt
  REAL(DP), PARAMETER :: nls_tol = 1.d-12
  ! REAL(DP) :: int_e
  
! Uses pressure of previous timestep as initial value for Newton's method
  p_in = prim(:,:,prs)

  do j= jbeg, jend
    do i= ibeg, iend

      ddd = cons(i,j,D)
      ssu = cons(i,j,Sx1) 
      ssv = cons(i,j,Sx2) 
      ttt = cons(i,j,tau)

! Initial value for Newton's method 
      if(( p_in(i,j)/p_in(i+1,j+1) > 1.d+2 ) .or. &
         ( p_in(i,j)/p_in(i+1,j+1) < 1.d-2 ) .or. &
         ( p_in(i,j)/p_in(i-1,j-1) > 1.d+2 ) .or. &
         ( p_in(i,j)/p_in(i-1,j-1) < 1.d-2 )) then 
      
        x(1) = sum(p_in(i-1:i+1,j)) / 3.d0  ! for large gradients
      else
        x(1) = p_in(i,j)                    ! for small gradients
      endif  

      iflag = 1
   
      call hybrd1 ( func_c2p, n, x, fvec, nls_tol, iflag)

      if( iflag == 0 .or. iflag == 2) THEN
        print*,iflag
        write(*,*) "hybrd1 do not converge"
        stop
      endif

      prim(i,j,prs) = x(1)
      prim(i,j,vx1) = ssu/(ttt + ddd + prim(i,j,prs))
      prim(i,j,vx2) = ssv/(ttt + ddd + prim(i,j,prs))
      W = 1.d0/dsqrt(1.d0 - prim(i,j,vx1)**2 - prim(i,j,vx2)**2)
      prim(i,j,rho) = ddd/W


      cs(i,j) = sound(prim(i,j,rho),prim(i,j,prs))
      prim(i,j,k_t) = k_tilde_p(prim(i,j,rho),prim(i,j,prs))
      ! prim(i,j,i_e) = (ttt + ddd*(1.d0-W) + (1.d0-W**2)*prim(i,j,prs))/ddd/W
      prim(i,j,i_e) = eps(prim(i,j,rho),prim(i,j,prs))

      ! write(*,*) int_e, prim(i,j,i_e)

    enddo
  enddo

  CONTAINS

  subroutine func_c2p(n,x,fvec,iflag)
    use constants,  only : i4b, dp
    use EoS,        only : pressure
    ! use EoS,        only : prs_q, prs_h, eps_q, eps_h, pstar
  implicit none
  integer(i4b) :: n, iflag
  real(dp), dimension(n) :: x, fvec
  real(dp)     :: eps_, W_, d_, u_, v_, p_

    p_   = x(1)
    u_   = ssu/(ttt + ddd + p_)
    v_   = ssv/(ttt + ddd + p_)
    W_   = 1.d0/sqrt(1.d0 - u_**2 - v_**2)
    d_   = ddd/W_
    eps_ = (ttt + ddd*(1.d0 - W_) + (1.d0 - W_**2)*p_) / ddd / W_

    ! if(p_ .ge. pstar) then
    !   fvec(1) = prs_q(d_, eps_) - p_
    ! else
    !   fvec(1) = prs_h(d_, eps_) - p_
    ! endif

    ! if (abs(eps_ - eps_q(d_,p_))/eps_q(d_,p_) .le. &
    !     abs(eps_ - eps_h(d_,p_))/eps_h(d_,p_) ) then

    !    fvec(1) = prs_q(d_, eps_) - p_

    ! else

    !   fvec(1) = prs_h(d_, eps_) - p_

    ! endif

    fvec(1) = pressure(d_, eps_) - p_
 
  end subroutine func_c2p

end subroutine c2p_RHD


! !==============================================================================!     
! ! for given values of the primitive variables "d,u,p" it calculates the        !
! ! Physical Flux vector [mass, momentum and energy] for relativistic            !
! ! one-dimensional Euler equations                                              !
! ! last modification:                                                           !
! ! date: 11/06/2014 ; by Victor Mourao Roque                                    !
! ! W -> Lorentz factor                                                          !
! !==============================================================================!


subroutine phys_flux_RHD(coord,prim,pflux)
  use mesh,       only : cells, im, jm
  use variables,  only : rho, vx1, vx2, prs, i_e, nc, np
  implicit none
  integer(i4b),             intent(in)  :: coord
  real(dp), dimension(cells(coord),np), intent(in)  :: prim
  real(dp), dimension(cells(coord),nc), intent(out) :: pflux
! Local Variables
  real(dp) :: E, V2, W, H, s
  integer(i4b) :: i, j

  select case(coord)

  case(1)

    do i = 1,im

      W = 1.d0/SQRT(1.d0 - prim(i,vx1)**2 - prim(i,vx2)**2)
      h = 1.d0 + prim(i,i_e) + prim(i,prs)/prim(i,rho)

      ! eq. (31) of [2]
      pflux(i,1) = prim(i,rho) * W * prim(i,vx1)
      !
      s = pflux(i,1)*h*W
      !
      pflux(i,2) = s * prim(i,vx1) + prim(i,prs)    
      pflux(i,3) = s * prim(i,vx2)
      pflux(i,4) = s - pflux(i,1)

    enddo


  case(2)

    do j = 1,jm
      W = 1.d0/SQRT(1.d0 - prim(j,vx1)**2 - prim(j,vx2)**2)
      h = 1.d0 + prim(j,i_e) + prim(j,prs)/prim(j,rho)
      !
      ! eq. (31) of [2]
      pflux(j,1) = prim(j,rho) * W * prim(j,vx2)
      !
      s = pflux(j,1)*h*W
      !
      pflux(j,2) = s * prim(j,vx1) 
      pflux(j,3) = s * prim(j,vx2) + prim(j,prs)
      pflux(j,4) = s - pflux(j,1)

    enddo

  end select
end subroutine phys_flux_RHD


!==============================================================================!     
! Subroutine to calculate the eigenvalues in interfaces of the cells           ! 
! j => i+1/2                                                                   !
! last modification:                                                           !
! date: 11/06/2014 ; by Victor Mourao Roque                                    !
! based on: [1] R. Donat, Journal of Computational Physics 146, 58–81 (1998).  !
! Notes:                                                                       !
!==============================================================================!

subroutine eigenval_RHD(coord,prim,cs,lambda)
  use mesh,       only : cells
  use variables,  only : vx1, vx2, np
  implicit none
  integer(i4b), intent(in)  :: coord
  real(dp)    , dimension(cells(coord),np), intent(in)  :: prim
  real(dp)    , dimension(cells(coord)),    intent(in)  :: cs


  real(dp)    , dimension(cells(coord),np), intent(out) :: lambda
! Local Variables
  integer(i4b) :: i,j
  real(dp)     :: cs2, V2, uu, vv, u, v

! eigenvalues eq.(14) and (15) of [1]
! acoustic waves (+/-) and material waves (0)

select case(coord)

 case(1)

  do concurrent(i=1:im)
    cs2 = cs(i)**2
    u = prim(i,vx1) ; v = prim(i,vx2)
    uu = u**2 ; vv = v**2
    V2 = uu + vv
    lambda(i,4) = ( u*(1.d0 - cs2) + cs(i) * &
             DSQRT((1.d0 - V2)*(1.d0 - uu - (V2 - uu)*cs2)))/(1.d0 - cs2*v2)
    lambda(i,2) = u
    lambda(i,3) = u
    lambda(i,1) = ( u*(1.d0-cs2) - cs(i) * &
             DSQRT((1.d0 - V2)*(1.d0 - uu - (V2 - uu)*cs2)))/(1.d0 - cs2*v2)

  enddo


  case(2)

    do concurrent(j=1:jm)
      cs2 = cs(j)**2
      u = prim(j,vx1) ; v = prim(j,vx2)
      uu = u**2 ; vv = v**2
      V2 = uu + vv

      lambda(j,4) = ( v*(1.d0 - cs2) + cs(j) * &
               DSQRT((1.d0 - V2)*(1.d0 - vv - (V2 - vv)*cs2)))/(1.d0 - cs2*v2)
      lambda(j,2) = v
      lambda(j,3) = v
      lambda(j,1) = (v*(1.d0 - cs2) - cs(j) * &
               DSQRT((1.d0 - V2)*(1.d0 - vv - (V2 - vv)*cs2)))/(1.d0 - cs2*V2)

    enddo
  
 end select

end subroutine eigenval_RHD


!==============================================================================!     
! Subroutine to calculate the eigenvectors (+,0,-) in interfaces of the cells  ! 
! i => i+1/2                                                                   !
! last modification:                                                           !
! date: 12/06/2014 ; by Victor Mourao Roque                                    !
! based on: [1] A. Suresh and H. T. Huynh, JCP 136, 83 (1997).                 !
! rev(:,:,1) and lev(:,:,1) are the positive eigenvectors                      !
! rev(:,:,2) and lev(:,:,2) are the zero eigenvectors                          !
! rev(:,:,3) and lev(:,:,3) are the negative eigenvectors                      !
!                                                                              !
! subrotina com autovetores para uma equação de estado arbitrária              !
! na equação de estado deve-se descrever a função gamma como                   !
! gamma = 1 + 1/d{Del p/Del epsilon} e a velocidade do som como                !
! cs = SQRT[{Del p/Del d} + p/d^2{Del p/Del epsilon}] =                        !
!    = SQRT[{Del p/Del d} + p*(gamma -1)/d]                                    !
!                                                                              !
!==============================================================================!


subroutine eigenvec_RHD(coord, prim, cs, lambda, lev, rev)
  use mesh,       only : cells
  use variables,  only : rho, vx1, vx2, prs, i_e, k_t, nc, np
  implicit none
  integer(i4b), intent(in)  :: coord
  real(dp), dimension(cells(coord),np),    intent(in)  :: prim
  real(dp), dimension(cells(coord)),       intent(in)  :: cs
  real(dp), dimension(cells(coord),nc),    intent(in)  :: lambda
  real(dp), dimension(cells(coord),nc,nc), intent(out) :: lev, rev
! Local Variables
  integer(i4b) :: i, j
  REAL(dp) :: Delta, h, W, kappa, A_plus, A_minus, Deltah2
  REAL(DP) :: cs2, vv, uu, uv,W2,hW, v2, u, v
  real(dp) :: d, p, kt

  select case(coord)

  case(1)
    do i = 1, cells(coord)
      cs2 = cs(i)**2
      u = prim(i,vx1) ; v = prim(i,vx2)
      uu = u**2 ; vv = v**2 ; uv = u*v
      V2 = uu + vv
      d = prim(i,rho)  ;  p = prim(i,prs)
      kt = prim(i,k_t)

      W = 1.d0/dsqrt(1.d0 - v2)
      W2 = 1.d0/(1.d0 - v2)
      h  = 1.d0 + prim(i,i_e) + p/d
      hW = h*W

      kappa = kt / (kt - cs2)

      A_plus = (1.d0 - uu)/(1.d0 - u*lambda(i,4))
      A_minus = (1.d0 - uu)/(1.d0 - u*lambda(i,1))


! right eigenvectors given by eq. 17 of [1]
      rev(i,:,4) = (/ 1.d0, hW*A_plus * lambda(i,4),          &
                                                 hW*v, h*W*A_plus - 1.d0/) ! r_+

      rev(i,:,3) = (/ W*v, 2.d0*h*W2*uv, h*(1.d0+2.d0*W2*vv), &
                                                      2.d0*h*W2*v - W*v/) ! r_2

      rev(i,:,2) = (/ kappa/hW, u , v, 1.d0 - kappa/hW /)                 ! r_1

      rev(i,:,1) = (/ 1.d0, hW*A_minus * lambda(i,1),         &
                                                h*W*v, hW*A_minus - 1.d0/) ! r_-


! left eigenvectors given by eq. 3.3 of [1]

      Delta = h**3*W*(kappa - 1.d0) * (1.d0 - uu) * & 
              (A_plus*lambda(i,4) - A_minus*lambda(i,1))

      Deltah2 = h**2/Delta


      lev(i,:,1) = Deltah2 * (/hW*A_plus*(u - lambda(i,4)) &
           - u - W2*(v2 - uu)*(2.d0*kappa-1.d0)*(u - A_plus*lambda(i,4)) &
           + kappa * A_plus*lambda(i,4), &
            1.d0 + W2*(v2 - uu)*(2.d0*kappa-1.d0)*(1.d0-A_plus) - kappa*A_plus,&
            W2 * v *(2.d0*kappa-1.d0)*A_plus*(u - lambda(i,4)),&
              kappa*A_plus*lambda(i,4)- u - &
              W2*(v2 - uu)*(2.d0*kappa-1.d0)*(u - A_plus*lambda(i,4)) /) ! l_-

      lev(i,:,2) = (W/(kappa-1.d0))*(/h-W, W*u, W*v,-W/) ! l_1

      lev(i,:,3) = (1.d0/h/(1.d0-uu))*(/-v, uv, 1.d0-uu, -v/)  ! l_2

      lev(i,:,4) = -Deltah2 * (/hW*A_minus*(u - lambda(i,1)) &
           - u - W2*(v2 - uu)*(2.d0*kappa-1.d0)*(u - A_minus*lambda(i,1)) &
           + kappa*A_minus*lambda(i,1), &
            1.d0 + W2*(v2 - uu)*(2.d0*kappa-1.d0)*(1.d0-A_minus) - kappa*A_minus,&
            W2*v*(2.d0*kappa-1.d0)*A_minus*(u - lambda(i,1)),&
              kappa*A_minus*lambda(i,1) - u - &
              W2*(v2 - uu)*(2.d0*kappa-1.d0)*(u - A_minus*lambda(i,1)) /) !l_+

    enddo

  case(2)
    do i = 1, cells(coord)
      cs2 = cs(i)**2
      u = prim(i,vx1) ; v = prim(i,vx2)
      uu = u**2 ; vv = v**2 ; uv = u*v
      V2 = uu + vv
      d = prim(i,rho)  ;  p = prim(i,prs)
      kt = prim(i,k_t)

      W = 1.d0/dsqrt(1.d0 - v2)
      W2 = 1.d0/(1.d0 - v2)
      h  = 1.d0 + prim(i,i_e) + p/d
      hW = h*W

      kappa = kt / (kt - cs2)

      A_plus =  (1.d0 - vv)/(1.d0 - v*lambda(i,4))
      A_minus = (1.d0 - vv)/(1.d0 - v*lambda(i,1))

! right eigenvectors given by eq. 17 of [1]
      rev(i,:,4) = (/ 1.d0, hW*u, hW*A_plus * lambda(i,4),    &
                                                     h*W * A_plus - 1.d0/) ! r_+

      rev(i,:,3) = (/ W*u, h*(1.d0+2.d0*W2*uu), 2.d0*h*W2*uv, &
                                                        2.d0*h*W2*u-W*u/) ! r_2

      rev(i,:,2) = (/ kappa/hW, u, v, 1.d0 - kappa/hW /)                  ! r_1

      rev(i,:,1) = (/ 1.d0, hW*u, hW*A_minus*lambda(i,1),    &
                                                      hW*A_minus - 1.d0/) ! r_-

! left eigenvectors given by eq. 3.3 of [1]

      Delta = h**3*W*(kappa - 1.d0) * (1.d0 - vv) * & 
                       (A_plus*lambda(i,4) - A_minus*lambda(i,1))

      Deltah2 = h**2/Delta

      lev(i,:,1) = Deltah2 * (/hW*A_plus*(v - lambda(i,4)) &
           - v - W2*(v2 - vv)*(2.d0*kappa-1.d0)*(v - A_plus*lambda(i,4)) &
           + kappa*A_plus*lambda(i,4), &
            W2*u*(2.d0*kappa-1.d0)*A_plus*(v - lambda(i,4)),&
            1.d0 + W2*(v2 - vv)*(2.d0*kappa-1.d0)*(1.d0-A_plus) - kappa*A_plus,&
              kappa*A_plus*lambda(i,4)-v - &
              W2*(v2 - vv)*(2.d0*kappa-1.d0)*(v - A_plus*lambda(i,4)) /) ! l_-

      lev(i,:,2) = (W/(kappa-1.d0))*(/h-W,W*u, W*v, -W/)                 ! l_1

      lev(i,:,3) = (1.d0/h/(1.d0-vv))*(/-u, 1.d0-vv, uv, -u /)           ! l_2

      lev(i,:,4) = -Deltah2 * (/hW*A_minus*(v - lambda(i,1)) &
           - v - W2*(v2 - vv)*(2.d0*kappa-1.d0)*(v-A_minus*lambda(i,1)) &
           + kappa*A_minus*lambda(i,1), &
          W2*u*(2.d0*kappa-1.d0)*A_minus*(v - lambda(i,1)),&
          1.d0 + W2*(v2 - vv)*(2.d0*kappa-1.d0)*(1.d0-A_minus) - kappa*A_minus,&
              kappa*A_minus*lambda(i,1)- v - &
          W2*(v2 - vv)*(2.d0*kappa-1.d0)*(v - A_minus*lambda(i,1)) /)    ! l_+

    enddo
  end select

end subroutine eigenvec_RHD



 SUBROUTINE test_vectorsHD(coord, lev, rev)
  use mesh,       only : cells, cbeg, cend
  use variables,  only : rho, vx1, vx2, prs, nc, np
 IMPLICIT NONE

! Input Variables
 INTEGER(I4B), INTENT(IN) ::coord
 
! Output Variables
 REAL(DP), DIMENSION(cells(coord),nc,nc), INTENT(IN) :: lev
 REAL(DP), DIMENSION(cells(coord),nc,nc), INTENT(IN) :: rev
 
! Local Variables
 INTEGER(I4B) :: i,j,k,l
 
!  REAL(DP), DIMENSION(cells(coord),nc) :: unit_test
 REAL(DP), parameter :: eps = 9.d-10
 REAL(DP), DIMENSION(4,4) :: test, right, left, unit!,inv_r,inv_l

  
  unit = 0.d0
  FORALL(i=1:4,j=1:4, i.EQ.j)
    unit(i,j) = 1.d0
  END FORALL 


! If the are correct normalized 
!  DO i = 1, cells(coord)
!    do j=1,nc

!      unit_test(i,j) = DOT_PRODUCT(rev(i,:,j),lev(i,:,j))
  
!      if(DABS(unit_test(i,j)-1.d0) .GE. eps ) then
!         write(*,*) unit_test(i,j) -1.d0, i
!        stop"The eigenvectors are not correctly normalized"
!      endif
!    ENDDO
!  ENDDO


  DO i = cbeg(coord), cend(coord)

    test = 0.d0

    right(:,:) = rev(i,:,:)


    DO j = 1,4
      DO k =1,4
        left(k,j) = lev(i,j,k)
      ENDDO
    ENDDO

      test = MATMUL(right,left)


!    DO j = 1,4
!      DO k =1,4
!        IF(DABS(unit(j,k)-test(j,k)) .GE. eps) THEN
!        write(*,*) unit(j,k), test(j,k)
!        STOP
!        ENDIF
!      ENDDO
!    ENDDO

  !  CALL inverse(right,inv_r,3)

    
   ! DO j = 1,4
    !  DO k =1,4
      !  IF(DABS(inv_r(j,k)-left(j,k)) .GE. eps) THEN
      !  write(*,*) inv_r(j,k),left(j,k)
      !  write(*,*) inv_r(j,k)-left(j,k)
     !   WRITE(*,*)
       ! STOP
    !    ENDIF
     ! ENDDO
    !ENDDO

      print*,i, coord
      do l = 1,4
      write(*,*) test(l,1),test(l,2),test(l,3),test(l,4)
      enddo
      write(*,*)

  ENDDO


 END SUBROUTINE


end module physics
