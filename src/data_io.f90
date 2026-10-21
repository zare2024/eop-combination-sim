! data_io.f90
! Readers for the three input files, plus the error handler used everywhere.
!   station positions   IVS scheduling catalog format (position.cat)
!   radio sources       IVS scheduling catalog format (source.cat.geodetic.good)
!   a priori EOP        IERS EOP 20 C04 format (eopc04.1962-now)

module data_io
  use, intrinsic :: iso_fortran_env, only: error_unit
  use constants
  implicit none
  private
  public :: station_t, source_t, eop_table_t
  public :: read_stations, read_sources, read_eop_c04, fail

  type :: station_t
    character(8) :: name = ''
    real(dp)     :: xyz(3) = 0.0_dp          ! terrestrial position [m]
  end type station_t

  type :: source_t
    character(8) :: name = ''                ! IVS name, e.g. 0552+398
    real(dp)     :: ra = 0.0_dp, dec = 0.0_dp  ! right ascension, declination [rad]
  end type source_t

  type :: eop_table_t
    ! daily values at 0h UTC: mjd, then xp, yp, dX, dY in rad and UT1-UTC in s
    real(dp), allocatable :: mjd(:), xp(:), yp(:), dut1(:), dx(:), dy(:)
  end type eop_table_t

  integer, parameter :: line_len = 512

contains

  ! Stations named in 'wanted', in that order, from a position.cat file.
  ! Data line:  ID  NAME  X  Y  Z  occupation-code  lon  lat  origin
  subroutine read_stations(filename, wanted, stations)
    character(*),                 intent(in)  :: filename
    character(*),                 intent(in)  :: wanted(:)
    type(station_t), allocatable, intent(out) :: stations(:)
    character(line_len) :: line
    character(8) :: id, name
    real(dp)     :: xyz(3)
    logical      :: found(size(wanted))
    integer      :: u, ios, i

    allocate(stations(size(wanted)))
    found = .false.
    u = open_file(filename)
    do
      read(u, '(a)', iostat=ios) line
      if (ios /= 0) exit
      if (.not. is_data_line(line, '*')) cycle
      read(line, *, iostat=ios) id, name, xyz
      if (ios /= 0) call fail('cannot read station line: ' // trim(line))
      do i = 1, size(wanted)
        if (name == wanted(i)) then
          stations(i) = station_t(name, xyz)
          found(i) = .true.
        end if
      end do
    end do
    close(u)

    do i = 1, size(wanted)
      if (.not. found(i)) call fail('station ' // trim(wanted(i)) // ' not found in ' // trim(filename))
    end do
  end subroutine read_stations

  ! All sources of a source.cat file.
  ! Data line:  IVS-name  common-name  hh mm ss.ssss  sdd mm ss.sssss  epoch ...
  subroutine read_sources(filename, sources)
    character(*),                intent(in)  :: filename
    type(source_t), allocatable, intent(out) :: sources(:)
    character(line_len) :: line
    character(8) :: name, common_name
    character(4) :: deg_text    ! read as text to keep the sign of declinations like -00 12 34
    integer      :: u, ios, n, rah, ram, decd, decm
    real(dp)     :: ras, decs, sgn

    allocate(sources(count_data_lines(filename, '*')))
    u = open_file(filename)
    n = 0
    do
      read(u, '(a)', iostat=ios) line
      if (ios /= 0) exit
      if (.not. is_data_line(line, '*')) cycle
      read(line, *, iostat=ios) name, common_name, rah, ram, ras, deg_text, decm, decs
      if (ios == 0) read(deg_text, *, iostat=ios) decd
      if (ios /= 0) call fail('cannot read source line: ' // trim(line))
      sgn = merge(-1.0_dp, 1.0_dp, index(deg_text, '-') > 0)
      n = n + 1
      sources(n)%name = name
      sources(n)%ra   = (rah + ram / 60.0_dp + ras / 3600.0_dp) * 15.0_dp * deg2rad
      sources(n)%dec  = sgn * (abs(decd) + decm / 60.0_dp + decs / 3600.0_dp) * deg2rad
    end do
    close(u)
  end subroutine read_sources

  ! Daily EOP from an IERS EOP 20 C04 file; comment lines start with '#'.
  ! Columns used: year month day hour MJD x["] y["] UT1-UTC[s] dX["] dY["]
  subroutine read_eop_c04(filename, eop)
    character(*),      intent(in)  :: filename
    type(eop_table_t), intent(out) :: eop
    character(line_len) :: line
    integer  :: u, ios, n, i, year, month, day, hour
    real(dp) :: mjd, x, y, ut1_utc, dx, dy

    n = count_data_lines(filename, '#')
    if (n < 4) call fail('need at least four days of EOP in ' // trim(filename))
    allocate(eop%mjd(n), eop%xp(n), eop%yp(n), eop%dut1(n), eop%dx(n), eop%dy(n))
    u = open_file(filename)
    n = 0
    do
      read(u, '(a)', iostat=ios) line
      if (ios /= 0) exit
      if (.not. is_data_line(line, '#')) cycle
      read(line, *, iostat=ios) year, month, day, hour, mjd, x, y, ut1_utc, dx, dy
      if (ios /= 0) call fail('cannot read EOP line: ' // trim(line))
      n = n + 1
      eop%mjd(n)  = mjd
      eop%xp(n)   = x * as2rad
      eop%yp(n)   = y * as2rad
      eop%dut1(n) = ut1_utc
      eop%dx(n)   = dx * as2rad
      eop%dy(n)   = dy * as2rad
    end do
    close(u)

    ! UT1-UTC jumps by 1 s at a leap second; interpolating across one would be wrong
    do i = 2, n
      if (abs(eop%dut1(i) - eop%dut1(i - 1)) > 0.5_dp) &
        call fail('leap second inside the EOP table is not handled: ' // trim(filename))
    end do
  end subroutine read_eop_c04

  ! Print an error message and stop with exit status 1.
  subroutine fail(message)
    character(*), intent(in) :: message
    write(error_unit, '(a)') 'Error: ' // trim(message)
    flush(error_unit)
    stop 1
  end subroutine fail

  integer function open_file(filename) result(u)
    character(*), intent(in) :: filename
    integer :: ios
    open(newunit=u, file=filename, status='old', action='read', iostat=ios)
    if (ios /= 0) call fail('cannot open ' // trim(filename))
  end function open_file

  ! A data line is neither blank nor a comment.
  logical function is_data_line(line, comment_char)
    character(*), intent(in) :: line
    character,    intent(in) :: comment_char
    character(len(line)) :: text
    text = adjustl(line)
    is_data_line = len_trim(text) > 0 .and. text(1:1) /= comment_char
  end function is_data_line

  integer function count_data_lines(filename, comment_char) result(n)
    character(*), intent(in) :: filename
    character,    intent(in) :: comment_char
    character(line_len) :: line
    integer :: u, ios
    n = 0
    u = open_file(filename)
    do
      read(u, '(a)', iostat=ios) line
      if (ios /= 0) exit
      if (is_data_line(line, comment_char)) n = n + 1
    end do
    close(u)
  end function count_data_lines

end module data_io
