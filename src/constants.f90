! Working precision, physical constants and unit conversions.
module constants
  use, intrinsic :: iso_fortran_env, only: real64
  implicit none
  private

  integer, parameter, public :: dp = real64   ! double precision throughout, as in SOFA

  real(dp), parameter, public :: pi       = acos(-1.0_dp)
  real(dp), parameter, public :: c_light  = 299792458.0_dp   ! speed of light [m/s]
  real(dp), parameter, public :: mjd_zero = 2400000.5_dp     ! Julian Date of MJD 0

  real(dp), parameter, public :: deg2rad = pi / 180.0_dp
  real(dp), parameter, public :: as2rad  = pi / 648000.0_dp     ! arcseconds to radians
  real(dp), parameter, public :: mas2rad = as2rad * 1.0e-3_dp
  real(dp), parameter, public :: uas2rad = as2rad * 1.0e-6_dp
  real(dp), parameter, public :: us2s    = 1.0e-6_dp            ! microseconds to seconds
  real(dp), parameter, public :: s2ps    = 1.0e12_dp            ! seconds to picoseconds
end module constants
