! template based on Jules standards:

! http://jules-lsm.github.io/coding_standards/index.html


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


MODULE writedata

  IMPLICIT NONE


!==============================================================================!  
! Description:                                                                 !
!   Writes the primitive variables into a .dat file                            !
!                                                                              !
! Author: Victor Mourão Roque                                                  !
! last modification:                                                           !
! date: 27/07/2017 ; by Victor Mourao Roque                                    ! 
!                                                                              !
! Refs:                                                                        !
!                                                                              !
!                                                                              !
!==============================================================================!


  PRIVATE  ! Private scope by default


! Module constants
  


! Module variables
  character(len=32), save :: fdir  
  character(len=32), save :: froot 
  character(len=32), save :: fext  
  


! Public members
  PUBLIC :: initialize_write, write_data


CONTAINS
  
  subroutine initialize_write()
    use parameters,  only : get_parameter
    use mesh,       only : icell, xmin, xmax
    use mesh,       only : jcell, ymin, ymax
    use mesh,       only : tstep, tini, tend
  implicit none
  character(len=32) :: filename


  fdir  = './data/' ; froot = 'data' ; fext  = 'dat'
  
   
    call get_parameter("froot", froot)
    call get_parameter("fdir", fdir)
    call get_parameter("fext", fext)

  WRITE(filename,100) trim(fdir), trim('info'), '.', trim('txt')
  100 FORMAT(A,A,A,A)

   OPEN(UNIT=20, FILE= trim(filename), FORM="FORMATTED", &
    STATUS="REPLACE", ACTION="WRITE")


    write(UNIT=20,FMT=110) 'x: ', icell, xmin, xmax
    write(UNIT=20,FMT=110) 'y: ', jcell, ymin, ymax
    write(UNIT=20,FMT=110) 't: ', tstep, tini, tend
    write(UNIT=20,FMT=120) 'var: ', 'dens', 'vx', 'vy', 'prs', 'cs'

    110 FORMAT(A3, I5.2, F10.3, F10.3)
    120 FORMAT(A5, A5, A3, A3, A4, A3)

   close(20)
  
  end subroutine initialize_write

  subroutine write_data()
    use mesh,       only : step
    use constants,  only : i4b, dp
    use mesh,       only : ibeg,iend, xc
    use mesh,       only : jbeg,jend, yc
    use variables,  only : prim, rho, vx1, vx2, prs, i_e, cs

  implicit none
! Local Variables
  character(len=32) :: filename
  integer(i4b)      :: i, j
  real(dp) :: E, int_e

   WRITE(filename,100) trim(fdir), trim(froot), '.', step, '.', trim(fext)
   100 FORMAT(A,A,A,I4.4,A,A)

   OPEN(UNIT=10, FILE= trim(filename), FORM="FORMATTED", &
    STATUS="REPLACE", ACTION="WRITE")

   
    do j = jbeg, jend
      do i = ibeg, iend

        int_e = prim(i,j,i_e)
        E = prim(i,j,rho) * (1.d0 + int_e)
     
        write(UNIT=10,FMT=200) prim(i,j,rho), prim(i,j,vx1), &
                             prim(i,j,vx2), prim(i,j,prs), cs(i,j), int_e, E
  
      enddo
      write(UNIT=10,FMT=200)
     enddo
   
   close(10)

   200 format(7(1X,ES20.8))

  end subroutine write_data


END MODULE writedata
