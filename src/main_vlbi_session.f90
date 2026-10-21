! main_vlbi_session.f90
! Milestone 1: simulate a 24-hour VLBI session and estimate the five EOP.
!
!   1. read the settings (Fortran namelist) and the input catalogs
!   2. schedule the scans
!   3. for every observation compute the partial derivatives and the
!      noise-free observed-minus-computed delay (O-C):
!        observed  geometric delay with the true EOP (a priori + injected offsets)
!                  plus the clock difference of the two stations
!        computed  geometric delay with the a priori EOP, clocks zero
!   4. noise-free check: the adjustment must return the injected truth
!   5. Monte Carlo: add white noise, estimate, repeat, and compare the scatter
!      of the estimates with their formal errors
!
! Parameter vector x:
!   x(1:5)  EOP corrections: xp, yp [uas], UT1-UTC [us], dX, dY [uas]
!   then, for every station except the clock reference: clock offset [ps],
!   rate [ps/h] and quadratic term [ps/h^2], time counted from mid-session
!
! Usage:  ./vlbi_session [settings file]          (default: config/session.nml)

program vlbi_session
  use constants
  use sofa_interface, only: iau_CAL2JD, iau_JD2CAL
  use data_io
  use frames
  use vlbi_model
  use scheduler
  use lsq
  implicit none

  ! ---------- settings: namelist /session/ ----------
  integer, parameter :: max_stations = 32
  integer        :: start_year = 2026, start_month = 5, start_day = 4
  integer        :: start_hour = 17, start_minute = 0
  real(dp)       :: duration_h = 24.0_dp
  real(dp)       :: scan_interval_min = 3.0_dp
  real(dp)       :: elevation_cutoff_deg = 5.0_dp
  character(8)   :: station_names(max_stations) = ''
  character(8)   :: clock_reference = ''
  real(dp)       :: injected_eop(n_eop) = 0.0_dp
  real(dp)       :: sigma_delay_ps = 25.0_dp
  integer        :: n_monte_carlo = 200
  integer        :: seed = 1
  character(256) :: station_file = 'data/position.cat'
  character(256) :: source_file  = 'data/source.cat.geodetic.good'
  character(256) :: eop_file     = 'data/eopc04_2026-05.txt'
  character(256) :: output_dir   = 'results'
  namelist /session/ start_year, start_month, start_day, start_hour, start_minute, &
                     duration_h, scan_interval_min, elevation_cutoff_deg,        &
                     station_names, clock_reference, injected_eop,                &
                     sigma_delay_ps, n_monte_carlo, seed,                         &
                     station_file, source_file, eop_file, output_dir

  ! EOP report units (uas, us) to SI units (rad, s)
  real(dp),     parameter :: unit_to_si(n_eop) = [uas2rad, uas2rad, us2s, uas2rad, uas2rad]
  character(5), parameter :: unit_label(n_eop) = ['[uas]', '[uas]', '[us] ', '[uas]', '[uas]']

  ! ---------- data ----------
  type(station_t), allocatable :: stations(:)
  type(source_t),  allocatable :: sources(:)
  type(eop_table_t)            :: eop_table
  type(scan_t),    allocatable :: scans(:)
  real(dp) :: mjd_start, mjd_mid

  ! ---------- parameters and observations ----------
  integer  :: n_sta, n_par, n_obs, i_ref
  integer,  allocatable :: clock_col(:)      ! first clock column of each station, 0 for the reference
  real(dp), allocatable :: clock_true(:, :)  ! simulated clock polynomial of each station
  real(dp), allocatable :: x_true(:)         ! simulated truth of all parameters
  real(dp), allocatable :: at(:, :)          ! design matrix A stored transposed: column k holds
                                             ! the partials of observation k (contiguous in memory)
  real(dp), allocatable :: oc_true(:)        ! noise-free O-C [ps]
  integer,  allocatable :: obs_scan(:), obs_sta(:, :)
  real(dp) :: weight                         ! 1 / sigma^2, the same for every delay

  call read_settings()
  call read_stations(station_file, station_names(1:n_sta), stations)
  call read_sources(source_file, sources)
  call read_eop_c04(eop_file, eop_table)

  call make_schedule(mjd_start, duration_h, scan_interval_min, elevation_cutoff_deg, &
                     stations, sources, eop_table, scans)
  call setup_parameters()
  call build_observations()

  call execute_command_line('mkdir -p ' // trim(output_dir))
  call write_schedule()
  call print_session()
  call noise_free_check()
  call monte_carlo()
  print '(/, a)', 'Files written to ' // trim(output_dir) // '/: schedule.txt, residuals_run1.txt, monte_carlo.txt'

contains

  subroutine read_settings()
    character(256) :: settings_file, message
    integer  :: u, ios, j
    real(dp) :: djm0, djm

    call get_command_argument(1, settings_file)
    if (len_trim(settings_file) == 0) settings_file = 'config/session.nml'
    open(newunit=u, file=settings_file, status='old', action='read', iostat=ios)
    if (ios /= 0) call fail('cannot open settings file ' // trim(settings_file))
    read(u, nml=session, iostat=ios, iomsg=message)
    if (ios /= 0) call fail('cannot read the settings: ' // trim(message))
    close(u)

    n_sta = count(station_names /= '')
    if (n_sta < 2) call fail('station_names must list at least two stations')
    if (len_trim(clock_reference) == 0) clock_reference = station_names(1)
    if (sigma_delay_ps <= 0.0_dp) call fail('sigma_delay_ps must be positive')
    if (n_monte_carlo < 1) call fail('n_monte_carlo must be at least 1')

    call iau_CAL2JD(start_year, start_month, start_day, djm0, djm, j)
    if (j /= 0) call fail('invalid session start date')
    mjd_start = djm + (start_hour + start_minute / 60.0_dp) / 24.0_dp
    mjd_mid   = mjd_start + duration_h / 48.0_dp
    weight    = 1.0_dp / sigma_delay_ps**2
  end subroutine read_settings

  ! Parameter layout and the simulated truth.
  subroutine setup_parameters()
    integer :: i

    i_ref = 0
    do i = 1, n_sta
      if (stations(i)%name == clock_reference) i_ref = i
    end do
    if (i_ref == 0) call fail('clock_reference ' // trim(clock_reference) // ' is not in station_names')

    ! VLBI observes only clock differences between stations, so a clock offset
    ! common to all stations cannot be estimated. One clock is held fixed (the
    ! reference); every other station gets an offset, a rate and a quadratic term.
    allocate(clock_col(n_sta))
    n_par = n_eop
    do i = 1, n_sta
      if (i == i_ref) then
        clock_col(i) = 0
      else
        clock_col(i) = n_par + 1
        n_par = n_par + 3
      end if
    end do

    ! simulated truth: the injected EOP offsets and random clock polynomials
    call seed_random(seed)
    allocate(clock_true(3, n_sta), x_true(n_par))
    clock_true = 0.0_dp
    x_true(1:n_eop) = injected_eop
    do i = 1, n_sta
      if (i == i_ref) cycle
      clock_true(1, i) = 1000.0_dp * gauss()    ! offset [ps]
      clock_true(2, i) = 30.0_dp * gauss()      ! rate [ps/h]
      clock_true(3, i) = 0.5_dp * gauss()       ! quadratic term [ps/h^2]
      x_true(clock_col(i):clock_col(i) + 2) = clock_true(:, i)
    end do
  end subroutine setup_parameters

  ! Partial derivatives and noise-free O-C for every baseline of every scan.
  subroutine build_observations()
    real(dp) :: eop_apriori(n_eop), eop_true(n_eop), r_apriori(3, 3), r_true(3, 3)
    real(dp) :: dr(3, 3, n_eop), k_crs(3), k_apriori(3), k_true(3), dk(3, n_eop)
    real(dp) :: r1(3), r2(3), th, tau_apriori, tau_true
    integer  :: s, m, i, j, i1, i2, p, k

    ! every pair of stations in a scan forms one baseline, i.e. one delay
    n_obs = 0
    do s = 1, size(scans)
      m = size(scans(s)%stations)
      n_obs = n_obs + m * (m - 1) / 2
    end do
    allocate(at(n_par, n_obs), oc_true(n_obs), obs_scan(n_obs), obs_sta(2, n_obs))
    at = 0.0_dp

    k = 0
    do s = 1, size(scans)
      ! Earth orientation at the scan epoch: a priori (C04) and truth (a priori + offsets)
      eop_apriori = interpolate_eop(eop_table, scans(s)%mjd)
      eop_true    = eop_apriori + injected_eop * unit_to_si
      r_apriori   = c2t_matrix(scans(s)%mjd, eop_apriori)
      r_true      = c2t_matrix(scans(s)%mjd, eop_true)
      dr          = rotation_partials(scans(s)%mjd, eop_apriori)

      ! source direction in the terrestrial frame, and its derivative w.r.t. each EOP
      k_crs     = source_vector(sources(scans(s)%source)%ra, sources(scans(s)%source)%dec)
      k_apriori = matmul(r_apriori, k_crs)
      k_true    = matmul(r_true, k_crs)
      do p = 1, n_eop
        dk(:, p) = matmul(dr(:, :, p), k_crs)
      end do
      th = (scans(s)%mjd - mjd_mid) * 24.0_dp      ! hours from mid-session

      do i = 1, size(scans(s)%stations) - 1
        do j = i + 1, size(scans(s)%stations)
          i1 = scans(s)%stations(i)
          i2 = scans(s)%stations(j)
          r1 = stations(i1)%xyz
          r2 = stations(i2)%xyz
          k = k + 1
          obs_scan(k)   = s
          obs_sta(:, k) = [i1, i2]

          tau_apriori = geometric_delay(r1, r2, k_apriori)
          tau_true    = geometric_delay(r1, r2, k_true)
          oc_true(k)  = (tau_true - tau_apriori) * s2ps + clock(i2, th) - clock(i1, th)

          ! EOP partials: the delay formula applied to the derivative of the direction
          do p = 1, n_eop
            at(p, k) = geometric_delay(r1, r2, dk(:, p)) * s2ps * unit_to_si(p)
          end do
          ! clock partials: station 2 enters with +, station 1 with -
          call add_clock_partials(at(:, k), i2, +1.0_dp, th)
          call add_clock_partials(at(:, k), i1, -1.0_dp, th)
        end do
      end do
    end do
  end subroutine build_observations

  ! Simulated clock error of station i [ps], th hours from mid-session
  real(dp) function clock(i, th)
    integer,  intent(in) :: i
    real(dp), intent(in) :: th
    clock = clock_true(1, i) + clock_true(2, i) * th + clock_true(3, i) * th**2
  end function clock

  subroutine add_clock_partials(column, i, sgn, th)
    real(dp), intent(inout) :: column(:)
    integer,  intent(in)    :: i
    real(dp), intent(in)    :: sgn, th
    integer :: c
    c = clock_col(i)
    if (c == 0) return                 ! reference clock: not estimated
    column(c:c + 2) = column(c:c + 2) + sgn * [1.0_dp, th, th**2]
  end subroutine add_clock_partials

  ! Least-squares solution for the reduced observations oc [ps].
  subroutine solve(oc, x, q, sigma0)
    real(dp), intent(in)  :: oc(:)
    real(dp), intent(out) :: x(:), q(:, :), sigma0
    type(normal_eq_t) :: neq
    integer :: k, info
    call neq_init(neq, n_par)
    do k = 1, n_obs
      call neq_add(neq, at(:, k), oc(k), weight)
    end do
    call neq_solve(neq, x, q, info)
    if (info /= 0) call fail('the normal matrix is singular: some parameter is not estimable')
    sigma0 = sqrt(neq_vtpv(neq, x) / (n_obs - n_par))
  end subroutine solve

  ! Without noise the adjustment must return the injected truth, up to the
  ! rounding of double precision. This tests the delay model, the partials and
  ! the solver together. Deviations are measured in units of the formal error
  ! sqrt(q_ii), i.e. against the statistical uncertainty of each parameter.
  subroutine noise_free_check()
    real(dp) :: x(n_par), q(n_par, n_par), sigma0, dev(n_par)
    integer  :: p
    logical  :: ok
    call solve(oc_true, x, q, sigma0)
    do p = 1, n_par
      dev(p) = abs(x(p) - x_true(p)) / sqrt(q(p, p))
    end do
    ok = maxval(dev) < 1.0e-3_dp
    print '(/, a)', 'Noise-free check: the estimates must equal the injected truth'
    print '(2x, a, es8.1, a, es8.1, 4x, a)', 'largest |estimate - truth| / formal error:   EOP', &
          maxval(dev(1:n_eop)), '   clocks', maxval(dev(n_eop + 1:)), merge('PASS', 'FAIL', ok)
  end subroutine noise_free_check

  subroutine monte_carlo()
    real(dp) :: x(n_par), q(n_par, n_par), sigma0, oc(n_obs), v(n_obs)
    real(dp) :: err(n_eop, n_monte_carlo), formal(n_eop, n_monte_carlo), sigma0_run(n_monte_carlo)
    integer  :: run, k, p

    do run = 1, n_monte_carlo
      do k = 1, n_obs
        oc(k) = oc_true(k) + sigma_delay_ps * gauss()
      end do
      call solve(oc, x, q, sigma0)
      err(:, run) = x(1:n_eop) - x_true(1:n_eop)
      do p = 1, n_eop
        formal(p, run) = sigma0 * sqrt(q(p, p))
      end do
      sigma0_run(run) = sigma0
      if (run == 1) then
        v = matmul(x, at) - oc           ! residuals v = A x - l
        call print_run1(x, q, sigma0, v)
        call write_residuals(v)
      end if
    end do
    call print_monte_carlo(err, formal, sigma0_run)
    call write_monte_carlo(err, sigma0_run)
  end subroutine monte_carlo

  ! ---------- report ----------

  subroutine print_session()
    character(512) :: line
    integer :: i, n_seen

    n_seen = count([(any(scans%source == i), i = 1, size(sources))])
    line = ''
    do i = 1, n_sta
      line = trim(line) // ' ' // trim(stations(i)%name)
      if (i == i_ref) line = trim(line) // '*'
    end do

    print '(a)', 'eop-combination-sim | Milestone 1: VLBI session simulation and EOP estimation'
    print '(a)', repeat('=', 78)
    print '(3a, f0.1, a, f0.1, a, f0.1, a)', 'Session       ', utc_string(mjd_start), ' UTC, ', &
          duration_h, ' h, one scan every ', scan_interval_min, ' min, cutoff ', elevation_cutoff_deg, ' deg'
    print '(a, i0, 2a)', 'Network       ', n_sta, ' stations:', trim(line) // '   (* clock reference)'
    print '(a, i0, a, i0, a)', 'Sources       ', size(sources), ' in the catalog, ', n_seen, ' observed'
    print '(a, i0, a, i0, a, i0, a, i0, a)', 'Observations  ', size(scans), ' scans, ', n_obs, &
          ' delays, ', n_par, ' parameters (5 EOP + ', n_par - n_eop, ' clock)'
    print '(a, f0.1, a)', 'Noise         white, ', sigma_delay_ps, ' ps per delay'
  end subroutine print_session

  subroutine print_run1(x, q, sigma0, v)
    real(dp), intent(in) :: x(:), q(:, :), sigma0, v(:)
    real(dp) :: formal(n_par), corr(n_eop, n_eop), wrms
    integer  :: p, r, i, c

    do p = 1, n_par
      formal(p) = sigma0 * sqrt(q(p, p))
    end do
    wrms = sqrt(sum(v**2) / n_obs)          ! all weights are equal

    print '(/, a, i0)', 'Run 1 of ', n_monte_carlo
    print '(2x, a)', 'Parameter      Injected   Estimated   Formal error   (Est - Inj) / formal'
    do p = 1, n_eop
      print '(2x, a4, 1x, a5, f12.2, f12.2, f15.2, f23.2)', eop_names(p), unit_label(p), &
            x_true(p), x(p), formal(p), (x(p) - x_true(p)) / formal(p)
    end do
    print '(2x, a, f0.1, a, f5.3)', 'Post-fit residuals: WRMS ', wrms, ' ps, sigma0 ', sigma0

    print '(/, 2x, 3a)', 'Station clocks relative to ', trim(stations(i_ref)%name), &
          ': estimate minus truth (formal error)'
    print '(2x, a)', 'Station         offset [ps]         rate [ps/h]     quadratic [ps/h^2]'
    do i = 1, n_sta
      c = clock_col(i)
      if (c == 0) cycle
      print '(2x, a8, f12.2, " (", f5.2, ")", f12.3, " (", f5.3, ")", f13.4, " (", f6.4, ")")', &
            stations(i)%name, x(c) - x_true(c), formal(c), x(c + 1) - x_true(c + 1), formal(c + 1), &
            x(c + 2) - x_true(c + 2), formal(c + 2)
    end do

    ! correlation coefficients of the EOP estimates, from the cofactor matrix
    do r = 1, n_eop
      do p = 1, n_eop
        corr(r, p) = q(r, p) / sqrt(q(r, r) * q(p, p))
      end do
    end do
    print '(/, 2x, a)', 'Correlations of the EOP estimates'
    print '(10x, 5a8)', (adjustr(eop_names(p)), p = 1, n_eop)
    do r = 1, n_eop
      print '(4x, a4, 2x, 5f8.3)', eop_names(r), corr(r, 1:r)
    end do
  end subroutine print_run1

  subroutine print_monte_carlo(err, formal, sigma0_run)
    real(dp), intent(in) :: err(:, :), formal(:, :), sigma0_run(:)
    real(dp) :: mean_err, rms_err, mean_formal
    integer  :: p, n

    n = size(err, 2)
    print '(/, a, i0, a)', 'Monte Carlo over ', n, ' noise realisations (error = estimate - truth)'
    print '(2x, a)', 'Parameter    Mean error   RMS error   Mean formal error   RMS / formal'
    do p = 1, n_eop
      mean_err    = sum(err(p, :)) / n
      rms_err     = sqrt(sum(err(p, :)**2) / n)
      mean_formal = sum(formal(p, :)) / n
      print '(2x, a4, 1x, a5, f11.2, f12.2, f20.2, f15.2)', eop_names(p), unit_label(p), &
            mean_err, rms_err, mean_formal, rms_err / mean_formal
    end do
    print '(2x, a, f5.3)', 'Mean sigma0 ', sum(sigma0_run) / n
  end subroutine print_monte_carlo

  ! ---------- output files ----------

  subroutine write_schedule()
    character(512) :: line
    integer :: u, s, i
    open(newunit=u, file=trim(output_dir) // '/schedule.txt', status='replace', action='write')
    write(u, '(a)') '#  scan  epoch (UTC)        source    stations'
    do s = 1, size(scans)
      line = ''
      do i = 1, size(scans(s)%stations)
        line = trim(line) // ' ' // stations(scans(s)%stations(i))%name
      end do
      write(u, '(i7, 2x, a16, 3x, a8, a)') s, utc_string(scans(s)%mjd), &
            sources(scans(s)%source)%name, trim(line)
    end do
    close(u)
  end subroutine write_schedule

  subroutine write_residuals(v)
    real(dp), intent(in) :: v(:)
    integer :: u, k
    open(newunit=u, file=trim(output_dir) // '/residuals_run1.txt', status='replace', action='write')
    write(u, '(a)') '#  obs   scan   hours   station1  station2  source    residual [ps]'
    do k = 1, n_obs
      write(u, '(i6, i7, f8.3, 3(2x, a8), f12.2)') k, obs_scan(k), &
            (scans(obs_scan(k))%mjd - mjd_start) * 24.0_dp, stations(obs_sta(1, k))%name, &
            stations(obs_sta(2, k))%name, sources(scans(obs_scan(k))%source)%name, v(k)
    end do
    close(u)
  end subroutine write_residuals

  subroutine write_monte_carlo(err, sigma0_run)
    real(dp), intent(in) :: err(:, :), sigma0_run(:)
    integer :: u, run
    open(newunit=u, file=trim(output_dir) // '/monte_carlo.txt', status='replace', action='write')
    write(u, '(a)') '#  run   estimate minus truth:  xp [uas]  yp [uas]  dUT1 [us]  dX [uas]  dY [uas]   sigma0'
    do run = 1, size(err, 2)
      write(u, '(i6, 23x, 5f10.3, f10.4)') run, err(:, run), sigma0_run(run)
    end do
    close(u)
  end subroutine write_monte_carlo

  ! 'YYYY-MM-DD hh:mm' for an epoch in MJD (UTC), rounded to the minute
  function utc_string(mjd) result(text)
    real(dp), intent(in) :: mjd
    character(16) :: text
    integer  :: minutes, day, iy, im, id, j
    real(dp) :: fd
    minutes = nint(mjd * 1440.0_dp)
    day     = minutes / 1440
    minutes = minutes - day * 1440
    call iau_JD2CAL(mjd_zero, real(day, dp), iy, im, id, fd, j)
    write(text, '(i4.4, "-", i2.2, "-", i2.2, 1x, i2.2, ":", i2.2)') iy, im, id, minutes / 60, mod(minutes, 60)
  end function utc_string

  ! ---------- random numbers ----------

  subroutine seed_random(seed_value)
    integer, intent(in) :: seed_value
    integer, allocatable :: s(:)
    integer :: n, i
    call random_seed(size=n)
    allocate(s(n))
    s = seed_value + 7919 * [(i, i = 1, n)]
    call random_seed(put=s)
  end subroutine seed_random

  ! Standard normal random number (Box-Muller transform)
  function gauss() result(g)
    real(dp) :: g, u1, u2
    call random_number(u1)
    call random_number(u2)
    g = sqrt(-2.0_dp * log(1.0_dp - u1)) * cos(2.0_dp * pi * u2)
  end function gauss

end program vlbi_session
