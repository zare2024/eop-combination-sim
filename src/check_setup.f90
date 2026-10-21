! check_setup.f90
! Toolchain check for eop-combination-sim
!   Part 1  LAPACK  least-squares line fit through the normal equations
!   Part 2  SOFA    Earth rotation angle and celestial-to-terrestrial matrix
! The EOP values in Part 2 are illustrative test values, not real ones.
! Reference values were computed with ERFA (the C version of SOFA) for the same inputs.

program check_setup
  implicit none
  integer,  parameter :: dp   = kind(1.0d0)
  real(dp), parameter :: pi   = acos(-1.0_dp)
  real(dp), parameter :: as2r = pi / (180.0_dp * 3600.0_dp)   ! arcsec to rad

  ! Part 1
  integer, parameter :: m = 6, n = 2
  real(dp) :: t(m), l(m), A(m, n), Nmat(n, n), rhs(n)
  integer  :: info

  ! Part 2
  double precision, external :: iau_ERA00
  real(dp), parameter :: dut1 = 0.05_dp                 ! UT1-UTC in s (test value)
  real(dp), parameter :: xp = 0.15_dp, yp = 0.35_dp     ! polar motion in arcsec (test values)
  real(dp), parameter :: era_ref = 209.605127656480_dp  ! deg
  real(dp) :: utc1, utc2, tai1, tai2, tt1, tt2, ut11, ut12, era_deg
  real(dp) :: rc2t(3, 3), rc2t_ref(3, 3)
  integer  :: j, i

  logical :: all_ok
  all_ok = .true.

  ! ---------------- Part 1: normal equations with LAPACK ----------------
  ! Fit l = a + b*t. Each row of A holds the partial derivatives of one observation.
  t = [0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]
  l = [2.1_dp, 2.4_dp, 3.1_dp, 3.4_dp, 4.1_dp, 4.4_dp]
  A(:, 1) = 1.0_dp     ! d l / d a
  A(:, 2) = t          ! d l / d b
  Nmat = matmul(transpose(A), A)   ! N = A^T A  (unit weights)
  rhs  = matmul(transpose(A), l)   ! n = A^T l
  call dposv('U', n, 1, Nmat, n, rhs, n, info)   ! Cholesky solve of N x = n, solution lands in rhs

  print '(a)', 'Part 1  LAPACK, normal equations'
  if (info /= 0) then
     print '(2x, a, i0)', 'dposv failed, info = ', info
     all_ok = .false.
  else
     call report('intercept a', rhs(1), 2.042857142857143_dp, 1.0e-12_dp)
     call report('slope b    ', rhs(2), 0.482857142857143_dp, 1.0e-12_dp)
  end if

  ! ---------------- Part 2: SOFA Earth rotation ----------------
  call iau_DTF2D('UTC', 2026, 10, 21, 12, 0, 0.0_dp, utc1, utc2, j)
  call check_status('iau_DTF2D', j)
  call iau_UTCTAI(utc1, utc2, tai1, tai2, j)
  call check_status('iau_UTCTAI', j)
  call iau_TAITT(tai1, tai2, tt1, tt2, j)
  call check_status('iau_TAITT', j)
  call iau_UTCUT1(utc1, utc2, dut1, ut11, ut12, j)
  call check_status('iau_UTCUT1', j)

  era_deg = iau_ERA00(ut11, ut12) * 180.0_dp / pi
  call iau_C2T06A(tt1, tt2, ut11, ut12, xp * as2r, yp * as2r, rc2t)

  rc2t_ref = reshape([ &
     -8.6944769732037719e-01_dp, -4.9401967431928301e-01_dp,  2.2941250539773124e-03_dp, &
      4.9401804507350094e-01_dp, -8.6945072344027530e-01_dp, -1.2691142337630632e-03_dp, &
      2.6215960882806297e-03_dp,  2.9910726138201410e-05_dp,  9.9999656316374330e-01_dp], &
     [3, 3], order=[2, 1])

  print '(/, a)', 'Part 2  SOFA, 2026-10-21 12:00 UTC with test EOP values'
  call report('Earth rotation angle (deg)', era_deg, era_ref, 1.0e-10_dp)
  print '(2x, a)', 'Celestial-to-terrestrial matrix (IAU 2006/2000A):'
  do i = 1, 3
     print '(4x, 3f16.12)', rc2t(i, :)
  end do
  call report('max deviation from reference', maxval(abs(rc2t - rc2t_ref)), 0.0_dp, 1.0e-12_dp)

  print '(a)', ''
  if (all_ok) then
     print '(a)', 'All checks passed. Toolchain ready.'
  else
     print '(a)', 'Some checks failed. Send me the full output.'
  end if

contains

  subroutine report(label, value, expected, tolerance)
    character(*), intent(in) :: label
    real(dp),     intent(in) :: value, expected, tolerance
    character(4) :: verdict
    if (abs(value - expected) <= tolerance) then
       verdict = 'PASS'
    else
       verdict = 'FAIL'
       all_ok = .false.
    end if
    print '(2x, a, " = ", es22.14, "   expected ", es22.14, "   ", a)', label, value, expected, verdict
  end subroutine report

  subroutine check_status(routine, status)
    character(*), intent(in) :: routine
    integer,      intent(in) :: status
    ! SOFA convention: 0 ok, positive = warning, negative = error
    if (status /= 0) then
       print '(2x, a, " returned status ", i0)', routine, status
       if (status < 0) all_ok = .false.
    end if
  end subroutine check_status

end program check_setup
