# eop-combination-sim

A Fortran sandbox that simulates space geodetic observations and estimates Earth Orientation Parameters (EOP) by least squares. It is a self-study project, and the models are deliberately simplified.

**Milestone 1 (this version):** a 24-hour VLBI session with six IVS stations. The program schedules the scans, simulates the delays, and estimates the five EOP (polar motion xp, yp, UT1−UTC and the celestial pole offsets dX, dY) together with station clock polynomials through normal equations. A Monte Carlo run then checks that the formal errors match the actual scatter of the estimates.

## How it works

1. **Earth orientation.** The rotation from the celestial frame (GCRS, aligned with the ICRF) to the terrestrial frame (ITRS) follows the CIO-based procedure of the IERS Conventions 2010 and is computed with SOFA:

   `r_ITRS = RPOM(xp, yp) · R3(ERA) · RC2I(X + dX, Y + dY) · r_GCRS`

   These are SOFA's matrices, the transposes of W, R and Q in the IERS notation [GCRS] = Q·R·W·[ITRS]. Precession-nutation uses the IAU 2006/2000A model, and the Earth rotation angle ERA follows from UT1. The a priori EOP are IERS EOP 20 C04 values, interpolated with 4-point Lagrange polynomials.

2. **Scheduling.** Every 3 minutes a greedy scheduler picks one source from the IVS geodetic catalog. Among the sources that at least two stations see above 5° elevation, it takes the one with the highest score: the number of stations that see it plus the hours since it was last observed, capped at 2.

3. **Simulated delays.** For a baseline b = r2 − r1 and a source in direction k, the geometric delay is `tau = -(b · k_ITRS) / c` with `k_ITRS = RC2T · k_GCRS`. The simulated observation uses the true EOP (a priori plus injected offsets), adds quadratic clock errors for every station except the reference, and 25 ps of white noise. The computed delay uses the a priori EOP.

4. **Estimation.** Gauss-Markov model with normal equations `N = AᵀPA`, `b = AᵀPl`, solved by Cholesky decomposition (LAPACK). The EOP partials come from central differences of the rotation matrix. VLBI only observes clock differences, so one station's clock is held fixed as the reference, and every other station gets an offset, a rate and a quadratic term.

5. **Checks.** Without noise, the adjustment must return the injected values. With noise, 200 Monte Carlo runs compare the RMS of the estimation errors with the mean formal error, and a ratio near 1 means the stochastic model is consistent. `make test` checks the building blocks separately: the rotation against SOFA's one-call `iau_C2T06A`, the numerical UT1 partial against the analytic Earth rotation formula, elevations, interpolation and the least-squares solver.

## Results with the default settings

```
Session       2026-05-04 17:00 UTC, 24.0 h, one scan every 3.0 min, cutoff 5.0 deg
Network       6 stations: WETTZELL* WESTFORD KOKEE ISHIOKA HOBART26 HART15M   (* clock reference)
Sources       342 in the catalog, 121 observed
Observations  480 scans, 2888 delays, 20 parameters (5 EOP + 15 clock)
Noise         white, 25.0 ps per delay

Noise-free check: the estimates must equal the injected truth
  largest |estimate - truth| / formal error:   EOP 1.4E-05   clocks 2.5E-05    PASS

Monte Carlo over 200 noise realisations (error = estimate - truth)
  Parameter    Mean error   RMS error   Mean formal error   RMS / formal
  xp   [uas]       0.92       11.66               11.81           0.99
  yp   [uas]       0.50       11.99               11.41           1.05
  dUT1 [us]       -0.02        0.58                0.59           0.99
  dX   [uas]       0.70        8.71                8.31           1.05
  dY   [uas]      -0.48        7.95                8.64           0.92
  Mean sigma0 0.999
```

The full report, including run 1 and the correlation matrix, is in `results/summary.txt`. The formal errors are optimistic compared with real sessions, which also contain troposphere and clock noise. With daily offsets, dX and dY are almost uncorrelated with polar motion, while UT1 correlates with xp and yp at about 0.36 to 0.38 in magnitude through the network geometry. The exact digits depend on the compiler's random number generator, but the statistics do not.

## Build and run

You need gfortran, LAPACK and the SOFA Fortran 77 library. On Ubuntu, including WSL:

```bash
sudo apt install gfortran make liblapack-dev libblas-dev

# SOFA Fortran 77, release 2023-10-11, from https://www.iausofa.org
wget -O sofa_f.tar.gz https://www.iausofa.org/s/sofa_f-20231011tar.gz
tar xzf sofa_f.tar.gz
cd sofa/20231011/f77/src && make && make test
mkdir -p /path/to/eop-combination-sim/lib && cp libsofa.a /path/to/eop-combination-sim/lib/
```

Then, in the project folder:

```bash
make check   # toolchain: LAPACK and SOFA
make test    # model building blocks
make run     # Milestone 1 session, report also saved to results/summary.txt
```

The settings (date, network, noise, injected offsets, number of Monte Carlo runs) are in `config/session.nml`. To use another settings file, pass it as an argument: `./vlbi_session my_settings.nml`.

## Simplifications

- Geometric delay only: no aberration, gravitational delay, retarded baseline, troposphere or ionosphere.
- Station and source positions are fixed. No tides, loading or plate motion.
- Clocks are quadratic polynomials and the noise is white. Real clocks and the troposphere add random-walk-like noise.
- Daily EOP are interpolated without the sub-daily ocean tide and libration terms.
- The scheduler ignores slewing, source flux and SNR.
- The C04 offsets dX, dY are added to SOFA's IAU 2006/2000A model. Any difference in reference model cancels here, because the same a priori enters both the simulated and the computed delays.

## Layout

```
config/session.nml          settings (Fortran namelist)
data/                       stations, sources and EOP, see data/README.md
src/constants.f90           precision, constants, unit conversions
src/sofa_interface.f90      interfaces to the SOFA routines
src/data_io.f90             readers for the catalogs and the C04 file
src/frames.f90              EOP interpolation and the celestial-to-terrestrial rotation
src/vlbi_model.f90          delay, elevation and rotation partials
src/scheduler.f90           greedy scheduler
src/lsq.f90                 normal equations and Cholesky solution
src/main_vlbi_session.f90   Milestone 1 program
src/test_model.f90          tests (make test)
src/check_setup.f90         toolchain check (make check)
results/                    output of make run
```

## Data and credits

- **SOFA**: IAU Standards of Fundamental Astronomy, https://www.iausofa.org. This project calls SOFA routines and is not itself SOFA software. The library is built locally and not included here.
- **IERS EOP 20 C04**: IERS Earth Orientation Center, Paris Observatory, https://hpiers.obspm.fr/iers/eop/eopc04/eopc04.1962-now. `data/eopc04_2026-05.txt` is an extract.
- **IVS scheduling catalogs**: https://github.com/nvi-inc/sked_catalogs, source of `data/position.cat` (extract) and `data/source.cat.geodetic.good` (unchanged copy).
- Developed with AI assistance (Claude, Anthropic).

## Roadmap

- Ring laser observation model (Sagnac frequency) and its EOP partials
- Lunar Laser Ranging range model
- Combination of VLBI, ring laser and LLR on the normal equation level, including data gaps and simulations with several ring lasers
