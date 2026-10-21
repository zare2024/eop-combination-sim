! scheduler.f90
! A simple greedy scheduler. Every interval_min minutes it places one scan on
! the source with the highest score among those that at least two stations
! see above the elevation cutoff:
!
!   score = number of stations that see it + hours since it was last observed (at most 2)
!
! The first term favours many baselines per scan, the second spreads the
! scans over the sky. Slewing, source flux and SNR are not modelled.

module scheduler
  use constants
  use data_io,    only: station_t, source_t, eop_table_t
  use frames,     only: interpolate_eop, c2t_matrix
  use vlbi_model, only: source_vector, local_vertical, elevation
  implicit none
  private
  public :: scan_t, make_schedule

  type :: scan_t
    real(dp) :: mjd = 0.0_dp             ! scan epoch (UTC)
    integer  :: source = 0               ! index into the source list
    integer, allocatable :: stations(:)  ! indices of the observing stations
  end type scan_t

  real(dp), parameter :: max_bonus_h = 2.0_dp

contains

  subroutine make_schedule(mjd_start, duration_h, interval_min, cutoff_deg, &
                           stations, sources, eop_table, scans)
    real(dp),                  intent(in)  :: mjd_start, duration_h, interval_min, cutoff_deg
    type(station_t),           intent(in)  :: stations(:)
    type(source_t),            intent(in)  :: sources(:)
    type(eop_table_t),         intent(in)  :: eop_table
    type(scan_t), allocatable, intent(out) :: scans(:)

    type(scan_t), allocatable :: slots(:)
    real(dp) :: up(3, size(stations)), k_crs(3, size(sources)), last_seen(size(sources))
    logical  :: vis(size(stations)), best_vis(size(stations))
    real(dp) :: r(3, 3), t, score, best_score, cutoff
    integer  :: n_slots, n_scans, slot, i, j, best

    cutoff  = cutoff_deg * deg2rad
    n_slots = int(duration_h * 60.0_dp / interval_min + 1.0e-9_dp)
    allocate(slots(n_slots))

    do i = 1, size(stations)
      up(:, i) = local_vertical(stations(i)%xyz)
    end do
    do j = 1, size(sources)
      k_crs(:, j) = source_vector(sources(j)%ra, sources(j)%dec)
    end do
    last_seen = mjd_start - 1.0_dp     ! as if every source was last seen a day ago

    n_scans  = 0
    best_vis = .false.
    do slot = 0, n_slots - 1
      t = mjd_start + slot * interval_min / 1440.0_dp
      r = c2t_matrix(t, interpolate_eop(eop_table, t))
      best = 0
      best_score = -1.0_dp
      do j = 1, size(sources)
        vis = visible(up, matmul(r, k_crs(:, j)), cutoff)
        if (count(vis) < 2) cycle
        score = count(vis) + min((t - last_seen(j)) * 24.0_dp, max_bonus_h)
        if (score > best_score) then
          best       = j
          best_score = score
          best_vis   = vis
        end if
      end do
      if (best == 0) cycle               ! no source seen by two stations: leave the slot empty

      n_scans = n_scans + 1
      slots(n_scans)%mjd      = t
      slots(n_scans)%source   = best
      slots(n_scans)%stations = pack([(i, i = 1, size(stations))], best_vis)
      last_seen(best) = t
    end do

    scans = slots(1:n_scans)
  end subroutine make_schedule

  ! Which stations see direction k_trs above the elevation cutoff
  function visible(up, k_trs, cutoff) result(vis)
    real(dp), intent(in) :: up(:, :), k_trs(3), cutoff
    logical :: vis(size(up, 2))
    integer :: i
    do i = 1, size(up, 2)
      vis(i) = elevation(up(:, i), k_trs) >= cutoff
    end do
  end function visible

end module scheduler
