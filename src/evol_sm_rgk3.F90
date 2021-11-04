submodule (evolution) rgk3_scheme
  contains

!=============================================================================!     
!  Runge-Kutta Third Order with three steps                                   !
! evolve the primitive variables in time t_0 to t_0 + dt  in one dimension    !
! last modification:                                                          !
! date: 27/06/2014 ; by Victor Mourao Roque                                   !
! based on: S. Gottlieb, D. I. Ketcheson, and C.-W. Shu,                      !
! J Sci Comput 38, 251–289 (2008).                                            !
!=============================================================================!  


module subroutine rgk3_3()
    use constants,  only  : i4b, dp
    use physics,    only  : prim2cons, cons2prim
    use mesh,       only  : im, iend, ibeg, jm, jend, jbeg, cells, dt
    use variables,  only  : prim, np, cs, cons, nc
    use numflux,    only  : num_flux
    use boundaries, only  : boundary2d
  implicit none
!local variables
  integer(i4b) :: i, j, coord
  real(dp), dimension(:,:,:), allocatable  :: cons_in, cons_1
  real(dp), dimension(:,:,:), allocatable  :: updated
  real(dp), dimension(:,:,:), allocatable  :: nflux_x, nflux_y
! 1D temporary arrays  
  real(dp), dimension(:,:), allocatable :: nflux_1d, prim_1d, cons_1d
  real(dp), dimension(:),   allocatable :: cs_1d

! Preparation Step
  
  
  call boundary2d()
  call prim2cons() 

  IF(.not. ALLOCATED(cons_in)) ALLOCATE(cons_in(im,jm,nc))
  
  cons_in = cons
   
  if(.not. allocated(nflux_x)) allocate(nflux_x(im,jm,nc))

! First step - X 
! 
 coord = 1
  if(.not. allocated(nflux_1d)) allocate(nflux_1d(cells(coord),nc))
  if(.not. allocated(prim_1d) ) allocate(prim_1d( cells(coord),np))
  if(.not. allocated(cons_1d) ) allocate(cons_1d( cells(coord),nc))
  if(.not. allocated(cs_1d) ) allocate(cs_1d( cells(coord)))

!$OMP PARALLEL SHARED(d,u,v,p,cs,nflux_x) &
!$OMP PRIVATE(j,nflux_1d,d_1,u_1,v_1,p_1,cs_1)
!$OMP DO
 
   y_fix_1: DO j = jbeg-1, jend

    prim_1d(:,:) = prim(:,j,:)
    cons_1d(:,:) = cons(:,j,:)
    cs_1d(:)     = cs(:,j)

    call num_flux(coord, prim_1d, cs_1d, cons_1d, nflux_1d) 

    nflux_x(:,j,:)  = nflux_1d(:,:)

  enddo y_fix_1

!$OMP END DO NOWAIT
!$OMP END PARALLEL

  deallocate(prim_1d)  ; deallocate(cons_1d)
  deallocate(nflux_1d) ; deallocate(cs_1d)


  if(.not. allocated(nflux_y)) allocate(nflux_y(im,jm,nc))

! First step - Y / 
! 
 coord = 2
  if(.not. allocated(nflux_1d)) allocate(nflux_1d(cells(coord),nc))
  if(.not. allocated(cons_1d) ) allocate(cons_1d( cells(coord),nc))
  if(.not. allocated(prim_1d) ) allocate(prim_1d( cells(coord),np))
  if(.not. allocated(cs_1d) ) allocate(cs_1d( cells(coord)))


!$OMP PARALLEL SHARED(d,u,v,p,cs,nflux_x) &
!$OMP PRIVATE(j,nflux_1d,d_1,u_1,v_1,p_1,cs_1)
!$OMP DO
 
  x_fix_1: DO i = ibeg-1, iend


    prim_1d(:,:) = prim(i,:,:)
    cons_1d(:,:) = cons(i,:,:)
    cs_1d(:)     = cs(i,:)

    call num_flux(coord, prim_1d, cs_1d, cons_1d, nflux_1d) 

    nflux_y(i,:,:)  = nflux_1d(:,:)

  enddo x_fix_1

!$OMP END DO NOWAIT
!$OMP END PARALLEL

  deallocate(prim_1d)  ; deallocate(cons_1d) 
  deallocate(nflux_1d) ; deallocate(cs_1d)

  if(.not. allocated(updated) ) allocate(updated(im,jm,nc))

  CALL update(nflux_x, nflux_y, updated)

  deallocate(nflux_x) ; deallocate(nflux_y) 

!$OMP PARALLEL SHARED(cons, cons_in, dt, updated,jbeg,jend,ibeg,iend) &
!$OMP PRIVATE(i,j)
!$OMP DO
  do j = jbeg, jend
    do i = ibeg, iend
      cons(i,j,:) = cons_in(i,j,:) + dt*updated(i,j,:)
    enddo
  enddo
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  DEALLOCATE(updated)

  call cons2prim()
  call boundary2d()

  IF(.not. ALLOCATED(cons_1)) ALLOCATE(cons_1(im,jm,nc))
  
  cons_1 = cons


! Second Step - X
!
 coord = 1
  if(.not. allocated(nflux_x)) allocate(nflux_x(im,jm,nc))
  if(.not. allocated(nflux_1d)) allocate(nflux_1d(cells(coord),nc))
  if(.not. allocated(cons_1d) ) allocate(cons_1d( cells(coord),nc))
  if(.not. allocated(prim_1d) ) allocate(prim_1d( cells(coord),np))
  if(.not. allocated(cs_1d) ) allocate(cs_1d( cells(coord)))

!$OMP PARALLEL SHARED(d,u,v,p,cs,nflux_x) &
!$OMP PRIVATE(j,nflux_1d,d_1,u_1,v_1,p_1,cs_1)
!$OMP DO
 
  y_fix_2: DO j = jbeg-1, jend

    prim_1d(:,:) = prim(:,j,:)
    cons_1d(:,:) = cons(:,j,:)
    cs_1d(:)     = cs(:,j)

    call num_flux(coord, prim_1d, cs_1d, cons_1d, nflux_1d) 

    nflux_x(:,j,:)  = nflux_1d(:,:)

  enddo y_fix_2

!$OMP END DO NOWAIT
!$OMP END PARALLEL

  deallocate(prim_1d)  ; deallocate(cons_1d)
  deallocate(nflux_1d) ; deallocate(cs_1d)


! Second Step - Y
!
 coord = 2
  if(.not. allocated(nflux_y)) allocate(nflux_y(im,jm,nc))
  if(.not. allocated(nflux_1d)) allocate(nflux_1d(cells(coord),nc))
  if(.not. allocated(cons_1d) ) allocate(cons_1d( cells(coord),nc))
  if(.not. allocated(prim_1d) ) allocate(prim_1d( cells(coord),np))
  if(.not. allocated(cs_1d) ) allocate(cs_1d( cells(coord)))


!$OMP PARALLEL SHARED(d,u,v,p,cs,nflux_x) &
!$OMP PRIVATE(j,nflux_1d,d_1,u_1,v_1,p_1,cs_1)
!$OMP DO

  x_fix_2: DO i = ibeg-1, iend

    prim_1d(:,:) = prim(i,:,:)
    cons_1d(:,:) = cons(i,:,:)
    cs_1d(:)     = cs(i,:)

    call num_flux(coord, prim_1d, cs_1d, cons_1d, nflux_1d) 

    nflux_y(i,:,:)  = nflux_1d(:,:)

  enddo x_fix_2

!$OMP END DO NOWAIT
!$OMP END PARALLEL

  deallocate(prim_1d)  ; deallocate(cons_1d) 
  deallocate(nflux_1d) ; deallocate(cs_1d)

  if(.not. allocated(updated) ) allocate(updated(im,jm,nc))

  CALL update(nflux_x, nflux_y, updated)

  deallocate(nflux_x) ; deallocate(nflux_y) 

!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!!!!!!
!!!!!  aqui entra a parte de source  !!!!!!
!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!!!!!!


!$OMP PARALLEL SHARED(cons, cons_in, cons_1, dt, updated) &
!$OMP PRIVATE(i,j)
!$OMP DO
  do j = jbeg, jend
    do i = ibeg, iend
      
      cons(i,j,:) = 0.75d0*cons_in(i,j,:) + 0.25d0*cons_1(i,j,:) + & 
                    0.25d0*dt*updated(i,j,:)
    enddo
  enddo
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  cons_1 = cons


! Third Step - X
!
 coord = 1
  if(.not. allocated(nflux_x)) allocate(nflux_x(im,jm,nc))
  if(.not. allocated(nflux_1d)) allocate(nflux_1d(cells(coord),nc))
  if(.not. allocated(cons_1d) ) allocate(cons_1d( cells(coord),nc))
  if(.not. allocated(prim_1d) ) allocate(prim_1d( cells(coord),np))
  if(.not. allocated(cs_1d) ) allocate(cs_1d( cells(coord)))

!$OMP PARALLEL SHARED(d,u,v,p,cs,nflux_x) &
!$OMP PRIVATE(j,nflux_1d,d_1,u_1,v_1,p_1,cs_1)
!$OMP DO
 
  y_fix_3: DO j = jbeg-1, jend

    prim_1d(:,:) = prim(:,j,:)
    cons_1d(:,:) = cons(:,j,:)
    cs_1d(:)     = cs(:,j)

    call num_flux(coord, prim_1d, cs_1d, cons_1d, nflux_1d) 

    nflux_x(:,j,:)  = nflux_1d(:,:)

  enddo y_fix_3

!$OMP END DO NOWAIT
!$OMP END PARALLEL

  deallocate(prim_1d)  ; deallocate(cons_1d)
  deallocate(nflux_1d) ; deallocate(cs_1d)


! Third Step - Y
!
 coord = 2
  if(.not. allocated(nflux_y)) allocate(nflux_y(im,jm,nc))
  if(.not. allocated(nflux_1d)) allocate(nflux_1d(cells(coord),nc))
  if(.not. allocated(cons_1d) ) allocate(cons_1d( cells(coord),nc))
  if(.not. allocated(prim_1d) ) allocate(prim_1d( cells(coord),np))
  if(.not. allocated(cs_1d) ) allocate(cs_1d( cells(coord)))


!$OMP PARALLEL SHARED(d,u,v,p,cs,nflux_x) &
!$OMP PRIVATE(j,nflux_1d,d_1,u_1,v_1,p_1,cs_1)
!$OMP DO

  x_fix_3: DO i = ibeg-1, iend

    prim_1d(:,:) = prim(i,:,:)
    cons_1d(:,:) = cons(i,:,:)
    cs_1d(:)     = cs(i,:)

    call num_flux(coord, prim_1d, cs_1d, cons_1d, nflux_1d) 

    nflux_y(i,:,:)  = nflux_1d(:,:)

  enddo x_fix_3

!$OMP END DO NOWAIT
!$OMP END PARALLEL

  deallocate(prim_1d)  ; deallocate(cons_1d) 
  deallocate(nflux_1d) ; deallocate(cs_1d)

  if(.not. allocated(updated) ) allocate(updated(im,jm,nc))

  CALL update(nflux_x, nflux_y, updated)

  deallocate(nflux_x) ; deallocate(nflux_y) 


!$OMP PARALLEL SHARED(cons, cons_in, cons_1, dt, updated) &
!$OMP PRIVATE(i,j)
!$OMP DO
  do j = jbeg, jend
    do i = ibeg, iend
      
      cons(i,j,:) = (cons_in(i,j,:) + 2.d0*cons_1(i,j,:) + & 
                    2.d0*dt*updated(i,j,:))/3.d0
    enddo
  enddo
!$OMP END DO NOWAIT
!$OMP END PARALLEL

! Final step  
  deallocate(updated) ; deallocate(cons_in) ; deallocate(cons_1)

  call cons2prim()
  call boundary2d()
  
  end subroutine rgk3_3

end submodule rgk3_scheme