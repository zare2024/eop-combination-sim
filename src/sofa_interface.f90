! Explicit interfaces for the SOFA Fortran 77 routines used here, so the compiler
! checks every call. SOFA itself is built separately (lib/libsofa.a).
! SOFA conventions: dates are two-part Julian Dates, angles are in radians.
module sofa_interface
  implicit none

  interface
    ! calendar date to MJD (djm0 + djm = JD)
    subroutine iau_CAL2JD(iy, im, id, djm0, djm, j)
      integer,          intent(in)  :: iy, im, id
      double precision, intent(out) :: djm0, djm
      integer,          intent(out) :: j
    end subroutine iau_CAL2JD

    ! Julian Date to calendar date and fraction of the day
    subroutine iau_JD2CAL(dj1, dj2, iy, im, id, fd, j)
      double precision, intent(in)  :: dj1, dj2
      integer,          intent(out) :: iy, im, id
      double precision, intent(out) :: fd
      integer,          intent(out) :: j
    end subroutine iau_JD2CAL

    ! UTC to TAI, applying the leap seconds
    subroutine iau_UTCTAI(utc1, utc2, tai1, tai2, j)
      double precision, intent(in)  :: utc1, utc2
      double precision, intent(out) :: tai1, tai2
      integer,          intent(out) :: j
    end subroutine iau_UTCTAI

    ! TAI to TT (TT = TAI + 32.184 s)
    subroutine iau_TAITT(tai1, tai2, tt1, tt2, j)
      double precision, intent(in)  :: tai1, tai2
      double precision, intent(out) :: tt1, tt2
      integer,          intent(out) :: j
    end subroutine iau_TAITT

    ! UTC to UT1, given UT1-UTC
    subroutine iau_UTCUT1(utc1, utc2, dut1, ut11, ut12, j)
      double precision, intent(in)  :: utc1, utc2, dut1
      double precision, intent(out) :: ut11, ut12
      integer,          intent(out) :: j
    end subroutine iau_UTCUT1

    ! CIP coordinates X, Y and CIO locator s, IAU 2006/2000A (TT)
    subroutine iau_XYS06A(date1, date2, x, y, s)
      double precision, intent(in)  :: date1, date2
      double precision, intent(out) :: x, y, s
    end subroutine iau_XYS06A

    ! celestial-to-intermediate matrix from X, Y and s
    subroutine iau_C2IXYS(x, y, s, rc2i)
      double precision, intent(in)  :: x, y, s
      double precision, intent(out) :: rc2i(3, 3)
    end subroutine iau_C2IXYS

    ! Earth rotation angle from UT1
    double precision function iau_ERA00(dj1, dj2)
      double precision, intent(in) :: dj1, dj2
    end function iau_ERA00

    ! TIO locator s' (TT)
    double precision function iau_SP00(date1, date2)
      double precision, intent(in) :: date1, date2
    end function iau_SP00

    ! polar motion matrix from xp, yp and s'
    subroutine iau_POM00(xp, yp, sp, rpom)
      double precision, intent(in)  :: xp, yp, sp
      double precision, intent(out) :: rpom(3, 3)
    end subroutine iau_POM00

    ! celestial-to-terrestrial matrix from its CIO-based parts
    subroutine iau_C2TCIO(rc2i, era, rpom, rc2t)
      double precision, intent(in)  :: rc2i(3, 3), era, rpom(3, 3)
      double precision, intent(out) :: rc2t(3, 3)
    end subroutine iau_C2TCIO

    ! the same matrix in one call, IAU 2006/2000A (only used in the tests)
    subroutine iau_C2T06A(tta, ttb, uta, utb, xp, yp, rc2t)
      double precision, intent(in)  :: tta, ttb, uta, utb, xp, yp
      double precision, intent(out) :: rc2t(3, 3)
    end subroutine iau_C2T06A

    ! geocentric to geodetic coordinates; n = 1 selects WGS84
    subroutine iau_GC2GD(n, xyz, elong, phi, height, j)
      integer,          intent(in)  :: n
      double precision, intent(in)  :: xyz(3)
      double precision, intent(out) :: elong, phi, height
      integer,          intent(out) :: j
    end subroutine iau_GC2GD
  end interface
end module sofa_interface
