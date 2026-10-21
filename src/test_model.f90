! test_model.f90
! Checks of the building blocks, independent of the session simulation.
! Run from the project folder with:  make test
!
!   1. the CIO-based rotation in frames equals SOFA's one-call iau_C2T06A
!   2. the numerical UT1 partial equals the analytic Earth-rotation formula
!   3. a source placed in a station's zenith has an elevation of 90 deg
!   4. EOP interpolation returns the tabulated value at a table epoch
!   5. least squares: exact line fit, and v'Pv from the normal equations
!      equals the sum of the squared residuals

program test_model
  use constants
  use sofa_interface
  use data_io
  use frames
  use vlbi_model
  use lsq
  implicit none

  ! two IVS stations from data/position.cat
  real(dp), parameter :: wettzell(3) = [4075539.5053_dp,   931735.6625_dp, 4801629.6156_dp]
  real(dp), parameter :: kokee(3)    = [-5543837.8378_dp, -2054566.3664_dp, 2387852.7011_dp]

  type(eop_table_t) :: table
  real(dp) :: mjd, eop(n_eop)
  integer  :: n_failed = 0

  call read_eop_c04('data/eopc04_2026-05.txt', table)
  mjd = table%mjd(8) + 0.3_dp            ! 07:12 UTC on a day inside the table
  eop = interpolate_eop(table, mjd)

  call test_rotation()
  call test_ut1_partial()
  call test_zenith()
  call test_interpolation()
  call test_least_squares()

  if (n_failed == 0) then
    print '(/, a)', 'All tests passed.'
  else
    print '(/, i0, a)', n_failed, ' test(s) failed.'
    stop 1
  end if

contains

  subroutine check(label, value, expected, tolerance)
    character(*), intent(in) :: label
    real(dp),     intent(in) :: value, expected, tolerance
    logical :: ok
    ok = abs(value - expected) <= tolerance
    if (.not. ok) n_failed = n_failed + 1
    print '(a, t62, es10.2, 3x, a)', label, value, merge('PASS', 'FAIL', ok)
  end subroutine check

  ! 1. Our step-by-step assembly (XYS06A, C2IXYS, ERA00, SP00, POM00, C2TCIO)
  !    must reproduce iau_C2T06A, which takes another route through SOFA.
  subroutine test_rotation()
    real(dp) :: e(n_eop), r_ours(3, 3), r_sofa(3, 3)
    real(dp) :: utc1, utc2, tai1, tai2, tt1, tt2, ut11, ut12
    integer  :: j
    e = eop
    e(i_dx) = 0.0_dp                     ! iau_C2T06A has no celestial pole offsets
    e(i_dy) = 0.0_dp
    r_ours = c2t_matrix(mjd, e)
    utc1 = mjd_zero + aint(mjd)
    utc2 = mjd - aint(mjd)
    call iau_UTCTAI(utc1, utc2, tai1, tai2, j)
    call iau_TAITT(tai1, tai2, tt1, tt2, j)
    call iau_UTCUT1(utc1, utc2, e(i_dut1), ut11, ut12, j)
    call iau_C2T06A(tt1, tt2, ut11, ut12, e(i_xp), e(i_yp), r_sofa)
    call check('1. rotation matrix vs iau_C2T06A, max difference', &
               maxval(abs(r_ours - r_sofa)), 0.0_dp, 1.0e-14_dp)
  end subroutine test_rotation

  ! 2. A change of UT1 rotates the Earth about the CIP axis z_cip. In the
  !    terrestrial frame the source direction then turns the other way:
  !        dk/dUT1 = -omega * (z_cip x k)
  !    so  dtau/dUT1 = (b . (z_cip x k)) * omega / c
  !    with omega = dERA/dUT1 and z_cip the third column of the polar motion matrix.
  subroutine test_ut1_partial()
    real(dp), parameter :: omega = 2.0_dp * pi * 1.00273781191135448_dp / 86400.0_dp  ! rad per s of UT1
    real(dp) :: k_crs(3), k_trs(3), dr(3, 3, n_eop), rpom(3, 3), numeric, analytic
    k_crs = source_vector(1.2_dp, 0.4_dp)
    k_trs = matmul(c2t_matrix(mjd, eop), k_crs)
    dr = rotation_partials(mjd, eop)
    numeric = geometric_delay(wettzell, kokee, matmul(dr(:, :, i_dut1), k_crs))
    call iau_POM00(eop(i_xp), eop(i_yp), 0.0_dp, rpom)    ! s' does not change the third column
    analytic = dot_product(kokee - wettzell, cross(rpom(:, 3), k_trs)) * omega / c_light
    call check('2. UT1 partial, numeric vs analytic, relative difference', &
               abs(numeric - analytic) / abs(analytic), 0.0_dp, 1.0e-6_dp)
  end subroutine test_ut1_partial

  ! 3. Direction of the local zenith, taken into the celestial frame and back
  subroutine test_zenith()
    real(dp) :: up(3), r(3, 3), k_crs(3)
    up = local_vertical(wettzell)
    r = c2t_matrix(mjd, eop)
    k_crs = matmul(transpose(r), up)
    call check('3. elevation of a source in the zenith [deg]', &
               elevation(up, matmul(r, k_crs)) / deg2rad, 90.0_dp, 1.0e-5_dp)
  end subroutine test_zenith

  ! 4. At a table epoch the Lagrange weights are 0, 1, 0, 0
  subroutine test_interpolation()
    real(dp) :: e(n_eop)
    e = interpolate_eop(table, table%mjd(6))
    call check('4. interpolated minus tabulated EOP at a table epoch', &
               maxval(abs(e - [table%xp(6), table%yp(6), table%dut1(6), table%dx(6), table%dy(6)])), &
               0.0_dp, 1.0e-15_dp)
  end subroutine test_interpolation

  ! 5. Line fit l = a + b t (same data as check_setup)
  subroutine test_least_squares()
    real(dp), parameter :: t(6) = [0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]
    real(dp), parameter :: l(6) = [2.1_dp, 2.4_dp, 3.1_dp, 3.4_dp, 4.1_dp, 4.4_dp]
    type(normal_eq_t) :: neq
    real(dp) :: x(2), q(2, 2), v(6)
    integer  :: i, info
    call neq_init(neq, 2)
    do i = 1, 6
      call neq_add(neq, [1.0_dp, t(i)], l(i), 1.0_dp)
    end do
    call neq_solve(neq, x, q, info)
    call check('5a. line fit: intercept', x(1), 2.042857142857143_dp, 1.0e-12_dp)
    call check('5b. line fit: slope', x(2), 0.482857142857143_dp, 1.0e-12_dp)
    v = x(1) + x(2) * t - l
    call check('5c. vTPv from normal equations minus sum(v**2)', &
               neq_vtpv(neq, x) - sum(v**2), 0.0_dp, 1.0e-12_dp)
  end subroutine test_least_squares

  pure function cross(a, b) result(c)
    real(dp), intent(in) :: a(3), b(3)
    real(dp) :: c(3)
    c = [a(2) * b(3) - a(3) * b(2), a(3) * b(1) - a(1) * b(3), a(1) * b(2) - a(2) * b(1)]
  end function cross

end program test_model
