! vlbi_model.f90
! Geometric VLBI delay, source visibility, and the partial derivatives of the
! celestial-to-terrestrial rotation with respect to the EOP.
!
! For a baseline b = r2 - r1 in the terrestrial frame and a source in unit
! direction k, the wavefront reaches station 2 later than station 1 by
!
!   tau = -(b . k_TRS) / c        with   k_TRS = RC2T * k_CRS
!
! Not modelled, in the simulated observations and in the analysis alike:
! aberration, gravitational delay, the retarded baseline, troposphere,
! ionosphere, station displacements (tides, loading) and source structure.

module vlbi_model
  use constants
  use sofa_interface, only: iau_GC2GD
  use data_io,        only: fail
  use frames,         only: n_eop, c2t_matrix
  implicit none
  private
  public :: source_vector, local_vertical, geometric_delay, elevation, rotation_partials

  ! central-difference steps: 1 mas for xp, yp, dX, dY and 0.1 ms for UT1-UTC
  real(dp), parameter :: step(n_eop) = [mas2rad, mas2rad, 1.0e-4_dp, mas2rad, mas2rad]

contains

  ! Unit vector towards a source at right ascension ra and declination dec [rad]
  pure function source_vector(ra, dec) result(k)
    real(dp), intent(in) :: ra, dec
    real(dp) :: k(3)
    k = [cos(dec) * cos(ra), cos(dec) * sin(ra), sin(dec)]
  end function source_vector

  ! Unit vector of the local vertical (WGS84 ellipsoid normal) at position xyz [m]
  function local_vertical(xyz) result(up)
    real(dp), intent(in) :: xyz(3)
    real(dp) :: up(3), elong, phi, height
    integer  :: j
    call iau_GC2GD(1, xyz, elong, phi, height, j)
    if (j /= 0) call fail('iau_GC2GD could not convert a station position')
    up = [cos(phi) * cos(elong), cos(phi) * sin(elong), sin(phi)]
  end function local_vertical

  ! Geometric delay tau = t2 - t1 [s] for stations at r1, r2 and direction k_trs.
  ! tau is linear in k_trs, so the same function maps the derivative of k_trs
  ! to the derivative of tau.
  pure function geometric_delay(r1, r2, k_trs) result(tau)
    real(dp), intent(in) :: r1(3), r2(3), k_trs(3)
    real(dp) :: tau
    tau = -dot_product(r2 - r1, k_trs) / c_light
  end function geometric_delay

  ! Elevation [rad] of direction k_trs seen from a station with local vertical up
  pure function elevation(up, k_trs) result(el)
    real(dp), intent(in) :: up(3), k_trs(3)
    real(dp) :: el
    el = asin(max(-1.0_dp, min(1.0_dp, dot_product(up, k_trs))))   ! clamp rounding beyond +-1
  end function elevation

  ! Derivatives of the rotation matrix RC2T with respect to each EOP,
  ! by central differences:   dR/dp = ( R(p + h) - R(p - h) ) / 2h
  ! Units: per radian for xp, yp, dX, dY and per second for UT1-UTC.
  function rotation_partials(mjd, eop) result(dr)
    real(dp), intent(in) :: mjd, eop(n_eop)
    real(dp) :: dr(3, 3, n_eop)
    real(dp) :: eop_plus(n_eop), eop_minus(n_eop)
    integer  :: p
    do p = 1, n_eop
      eop_plus  = eop
      eop_minus = eop
      eop_plus(p)  = eop(p) + step(p)
      eop_minus(p) = eop(p) - step(p)
      dr(:, :, p) = (c2t_matrix(mjd, eop_plus) - c2t_matrix(mjd, eop_minus)) / (2.0_dp * step(p))
    end do
  end function rotation_partials

end module vlbi_model
