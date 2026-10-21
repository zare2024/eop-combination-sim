! constants.f90
! Working precision, physical constants and unit conversions.

module constants
  use, intrinsic :: iso_fortran_env, only: real64
  implicit none
  private

  ! double precision everywhere, as in SOFA
  integer, parameter, public :: dp = real64

  real(dp), parameter, public :: pi       = acos(-1.0_dp)
  real(dp), parameter, public :: c_light  = 299792458.0_dp   ! speed of light [m/s]
  real(dp), parameter, public :: mjd_zero = 2400000.5_dp     ! Julian Date of MJD 0

  ! unit conversions
  real(dp), parameter, public :: deg2rad = pi / 180.0_dp
  real(dp), parameter, public :: as2rad  = pi / 648000.0_dp     ! arcsecond to radian
  real(dp), parameter, public :: mas2rad = as2rad * 1.0e-3_dp   ! milliarcsecond to radian
  real(dp), parameter, public :: uas2rad = as2rad * 1.0e-6_dp   ! microarcsecond to radian
  real(dp), parameter, public :: us2s    = 1.0e-6_dp            ! microsecond to second
  real(dp), parameter, public :: s2ps    = 1.0e12_dp            ! second to picosecond

end module constants
