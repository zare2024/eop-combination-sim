# Input data

## position.cat

Extract from `position.cat` of the IVS scheduling catalogs (https://github.com/nvi-inc/sked_catalogs, catalog version 2026-03-06_iGSFC). The station lines are copied unchanged. Positions come from the GSFC VLBI solution 2020c at epoch 2021-01-01, and station velocities are not applied.

## source.cat.geodetic.good

Unchanged copy of the IVS geodetic source catalog from the same repository (version 2024May28_iGSFC). It lists 342 sources with J2000 positions from ICRF2, ICRF3 and GSFC solutions; the last column of each line names the origin.

## eopc04_2026-05.txt

Extract (2026-04-26 to 2026-05-13) of the IERS EOP 20 C04 series, file `eopc04.1962-now` from the IERS Earth Orientation Center at Paris Observatory (https://hpiers.obspm.fr/iers/eop/eopc04/eopc04.1962-now), with its header. Daily values at 0h UTC: pole coordinates x, y, UT1−UTC, celestial pole offsets dX, dY, and their errors. To simulate another date, replace it with the full file from the IERS and change the date in `config/session.nml`.
