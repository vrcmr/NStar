program drive2d
  use constants,      only : i4b, dp, sp
  use mesh
  use parameters,     only : read_parameters, finalize_parameters
  use parameters,     only : get_parameter
  use variables,      only : initialize_variables, finalize_variables
  !use variables,      only : print_prim_var, print_var_bc, print_cons_var
  use evolution
  use init_cond,      only : initia
  use boundaries,     only : initialize_boundaries, finalize_boundaries
  use physics,        only : initialize_physics, prim2cons, cons2prim
  use physics,        only : phys_flux
!   use variables,      only : prim, cs
  use fluxsplit,      only : initialize_split, finalize_split
  use reconstruction, only : initialize_reconstr, finalize_reconstr
  use writedata
  use EoS,            only : initialize_EoS
  use sources,        only : initialize_source
  
  IMPLICIT NONE

! local variables  
  integer(i4b)            :: iterm, n_int
  integer(i4b), parameter :: sucess  = 0
  integer(i4b), parameter :: ntmax   = 100000
  real(dp),     parameter :: timetol = 1.d-15 
  logical :: master = .true.
  character(len=32)  :: source_term = "no"

! final time to produce a snapshot of the simulation in one particular step
!
  real(dp)                :: timeout 

  call read_parameters(master, iterm)

if(iterm /= sucess) then

  write(*,*) "Parametros não lidos!"
  
endif


  call init_mesh() ! ; call print_mesh()
  call initialize_variables()
  call initialize_evolution()
  call initialize_boundaries()
  call initialize_physics()
  call initialize_split()
  call initialize_reconstr()
  call initialize_write()
  call initialize_EoS()


  call get_parameter("source_term",source_term)

  select case(source_term)

  case("Y","Yes","on","ON")

    call initialize_source()

  case default

    write(*,*) "No source term"

  endselect


!= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = ! 
!               Beginning of  Hydrodynamic calculation                         !
!= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = !   
  
  write(*,*)'!= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = ', & 
            '= = = = = = = = !'  
  write(*,*)'!                                                              ', & 
            '                !'  
  write(*,*)'!         Applying the initial conditions to the primitive', &   
            ' variables           !' 
  write(*,*)'!                                                              ', & 
            '                !'  
  write(*,*)'!= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = ', & 
            '= = = = = = = = !'
  WRITE(*,*) 
  

  call initia()

!   WRITE(FNAME,100) 'sod.',step, '.dat'
!   OPEN(UNIT=10, FILE=trim(fname), ACTION="WRITE")
!   do i=ibeg,iend
!     write(10,200) xc(i),prim(i,1),prim(i,2),prim(i,3)
!   enddo
!   close(10)

  call write_data()

! Loop that determine when the code will gives the output evolving just dts 

  timeout = tini
  DOdts: do while ( real(time,sp) < real(tend,sp) .and. timeout < tend)

    timeout = time + dts
    step = step + 1
  
    n_int = 0

! All evolution happens inside of loop EVOLVTIME    
    EVOLVTIME: do while (real(time,sp) < real(timeout,sp) )
     
      n_int = n_int + 1
  
! Determines dt based on CFL condition     

      call timestep(timeout)

! find current time   
      time = time + dt
   write(*,*) "tend: ", real(tend,sp)
      print*,"tempo:   ", real(time,sp)
      print*,"timeout: ", real(timeout,sp)
      print*,"dt:      ", dt
      write(*,*)
      write(*,*)

      call time_evolution()

    enddo EVOLVTIME


  call write_data()
   write(*,*) "snapshot:",step
    if(n_int >= ntmax) then
      write(*,*)
      write(*,*) "ATTENTION!"
      write(*,*) "The simulation has exceeded the maximum number of steps and",&
                 " will be finalized!"
      STOP
      write(*,*)
      stop
    endif
  
  enddo DOdts

! finalize all arrays, modules and variables 
  
  call finalize_reconstr()
  call finalize_split()
  call finalize_boundaries()
  call finalize_evolution()
  call finalize_variables()
  call finalize_mesh()  
  call finalize_parameters()

end program
