module constants
    use, intrinsic :: iso_fortran_env, only : real32, real64,real128
    use, intrinsic :: iso_fortran_env, only : int8, int16, int32, int64

    implicit none

! everthing is public by default
    public

! constants single, double and quadruple precision reals:
    integer, parameter :: sp = real32,  & ! single precision [4 bytes ]
                          dp = real64,  & ! double precision [8 bytes ]
                          qp = real128    ! quad precision   [16 bytes]

! constants for 8, 4, 2, and 1 byte integers:
    integer, parameter:: i1b = int8 ,  & ! [1 byte]
                         i2b = int16,  & ! [2 byte]
                         i4b = int32,  & ! [4 byte]
                         i8b = int64     ! [8 byte]

    integer, parameter :: spc = kind((real32,real32))
    integer, parameter :: dpc = kind((real64,real64))
    integer, parameter :: lgt = kind(.true.)

    real(SP), parameter :: PI=3.141592653589793238462643383279502884197_sp
    real(SP), parameter :: PIO2=1.57079632679489661923132169163975144209858_sp
    real(SP), parameter :: TWOPI=6.283185307179586476925286766559005768394_sp
    real(SP), parameter :: SQRT2=1.41421356237309504880168872420969807856967_sp 
    real(SP), parameter :: EULER=0.5772156649015328606065120900824024310422_sp
    real(DP), parameter :: PI_D=3.141592653589793238462643383279502884197_dp
    real(DP), parameter :: PIO2_D=1.57079632679489661923132169163975144209858_dp
    real(DP), parameter :: TWOPI_D=6.283185307179586476925286766559005768394_dp

!physical constants
    real(DP), parameter :: G_grav = 6.673d-8       !cm^3 g^-1 s^-2  
    real(DP), parameter :: clight = 2.99792458d10  !cm/s
! Ref vel. da luz: http://physics.nist.gov/cgi-bin/cuu/Value?c          



end module constants

