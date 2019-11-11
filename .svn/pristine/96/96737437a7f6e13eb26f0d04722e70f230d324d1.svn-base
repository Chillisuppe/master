!*******************************************************************************
!
!> @file  DG_storm_drag_garrat.f90
!> @brief contains module DG_storm_drag
!
!*******************************************************************************

MODULE DG_storm_drag

  USE GRID_api
  PUBLIC :: wind_drag

  CONTAINS

!*******************************************************************************
! DESCRIPTION of [FUNCTION wind_drag]:
!> @brief Calculates the drag coefficient for wind stress.
!>        This approach was first presented in the review "Review of Drag Coefficients
!>        over Oceans and Continents" (1977)
!>        It is only valid for wind speeds w between [4,21] m/s and implements Charnock's
!>        relation by a linear form C_d * 10^3 = 0.75 + 0.067 * w.
!>
!> @param[in]     r_wind_speed    magnitude of the wind in the cell
!> @return                        coefficient of drag
!>
  FUNCTION wind_drag(r_wind_speed) RESULT (r_wind_drag)

    IMPLICIT NONE

    REAL (KIND=GRID_SR), INTENT(in)   :: r_wind_speed
    REAL (KIND=GRID_SR)               :: r_wind_drag

    REAL (KIND=GRID_SR)               :: r_min_drag

!-- Minimal wind drag
    r_min_drag = 0.002_GRID_SR

    r_wind_drag = MIN(r_min_drag, (0.75_GRID_SR + 0.067_GRID_SR * r_wind_speed) * 0.001_GRID_SR)
    r_wind_drag = r_wind_drag * 1000.0_GRID_SR

  END FUNCTION wind_drag

!*******************************************************************************
END MODULE DG_storm_drag
