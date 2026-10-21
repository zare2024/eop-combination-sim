! frames.f90
! Earth orientation: a priori EOP from the C04 table, and the rotation from the
! celestial frame (GCRS, axes aligned with the ICRF) to the terrestrial frame
! (ITRS), CIO-based as in the IERS Conventions 2010, computed with SOFA:
!
!   r_ITRS = RPOM(xp, yp, s') * R3(ERA) * RC2I(X + dX, Y + dY, s) * r_GCRS
!
!   RC2I  precession-nutation: IAU 2006/2000A model plus the offsets dX, dY
!   ERA   Earth rotation angle, a linear function of UT1 = UTC + (UT1-UTC)
!   RPOM  polar motion xp, yp
!
! In the notation of the IERS Conventions, [GCRS] = Q R W [ITRS], so these
! SOFA matrices are the transposes: RC2I = Q^T, R3(ERA) = R^T, RPOM = W^T.
!
! An EOP vector holds the five parameters in SI units:
!   eop = (xp, yp, UT1-UTC, dX, dY)   in   (rad, rad, s, rad, rad)

module frames
  use constants
  use sofa_interface
  use data_io, only: eop_table_t, fail
  implicit none
  private
  public :: n_eop, i_xp, i_yp, i_dut1, i_dx, i_dy, eop_names
  public :: interpolate_eop, c2t_matrix

  integer, parameter :: n_eop = 5
  integer, parameter :: i_xp = 1, i_yp = 2, i_dut1 = 3, i_dx = 4, i_dy = 5
  character(4), parameter :: eop_names(n_eop) = ['xp  ', 'yp  ', 'dUT1', 'dX  ', 'dY  ']

contains

  ! A priori EOP at epoch mjd (UTC): 4-point Lagrange interpolation of the daily table.
  function interpolate_eop(table, mjd) result(eop)
    type(eop_table_t), intent(in) :: table
    real(dp),          intent(in) :: mjd
    real(dp) :: eop(n_eop)
    real(dp) :: w(4), tab(4)
    integer  :: k, i, j

    ! table interval [k, k+1] containing mjd, with one more point on each side
    k = 0
    do i = 2, size(table%mjd) - 2
      if (table%mjd(i) <= mjd .and. mjd < table%mjd(i + 1)) k = i
    end do
    if (k == 0) call fail('epoch outside the EOP table (it needs two days on each side)')

    ! Lagrange weights of the points k-1, k, k+1, k+2
    tab = table%mjd(k - 1:k + 2)
    do i = 1, 4
      w(i) = 1.0_dp
      do j = 1, 4
        if (j /= i) w(i) = w(i) * (mjd - tab(j)) / (tab(i) - tab(j))
      end do
    end do

    eop(i_xp)   = dot_product(w, table%xp(k - 1:k + 2))
    eop(i_yp)   = dot_product(w, table%yp(k - 1:k + 2))
    eop(i_dut1) = dot_product(w, table%dut1(k - 1:k + 2))
    eop(i_dx)   = dot_product(w, table%dx(k - 1:k + 2))
    eop(i_dy)   = dot_product(w, table%dy(k - 1:k + 2))
  end function interpolate_eop

  ! Rotation matrix from the celestial (GCRS) to the terrestrial (ITRS) frame
  ! at epoch mjd (UTC) for the given EOP.
  function c2t_matrix(mjd, eop) result(rc2t)
    real(dp), intent(in) :: mjd, eop(n_eop)
    real(dp) :: rc2t(3, 3)
    real(dp) :: utc1, utc2, tai1, tai2, tt1, tt2, ut11, ut12
    real(dp) :: x, y, s, era, sp, rc2i(3, 3), rpom(3, 3)
    integer  :: j

    ! Two-part date: Julian Date of 0h UTC plus the fraction of the day.
    ! This keeps the full time resolution, which the UT1 partials need.
    utc1 = mjd_zero + aint(mjd)
    utc2 = mjd - aint(mjd)

    ! time scales: TT for precession-nutation, UT1 for Earth rotation
    call iau_UTCTAI(utc1, utc2, tai1, tai2, j)
    call check_status(j, 'iau_UTCTAI')
    call iau_TAITT(tai1, tai2, tt1, tt2, j)
    call check_status(j, 'iau_TAITT')
    call iau_UTCUT1(utc1, utc2, eop(i_dut1), ut11, ut12, j)
    call check_status(j, 'iau_UTCUT1')

    ! precession-nutation: CIP coordinates from the model, corrected by dX, dY
    call iau_XYS06A(tt1, tt2, x, y, s)
    call iau_C2IXYS(x + eop(i_dx), y + eop(i_dy), s, rc2i)

    ! Earth rotation angle
    era = iau_ERA00(ut11, ut12)

    ! polar motion, including the TIO locator s'
    sp = iau_SP00(tt1, tt2)
    call iau_POM00(eop(i_xp), eop(i_yp), sp, rpom)

    ! rc2t = rpom * R3(era) * rc2i
    call iau_C2TCIO(rc2i, era, rpom, rc2t)
  end function c2t_matrix

  subroutine check_status(j, routine)
    integer,      intent(in) :: j
    character(*), intent(in) :: routine
    character(12) :: text
    if (j /= 0) then
      write(text, '(i0)') j
      call fail(routine // ' returned status ' // trim(text) // &
                ' (date outside the range this SOFA release accepts?)')
    end if
  end subroutine check_status

end module frames
