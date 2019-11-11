!*******************************************************************************
!
!> @file  DG_storm_drag_weisberg.f90
!> @brief contains module DG_storm_drag
!
!*******************************************************************************

MODULE DG_storm_drag

  USE GRID_api
  PUBLIC :: wind_drag

  CONTAINS

!*******************************************************************************
! DESCRIPTION of [FUNCTION wind_drag]:
!> @brief Calculates the drag coefficient for wind stress given a wind speed.
!>        Based on the presentation of the Tampa Bay study in Weisberg and
!>        Zheng (2006); see also Large and Pond (1981)
!>
!> @param[in]     r_wind_speed    magnitude of the wind in the cell
!> @return                        coefficient of drag
!>
  FUNCTION wind_drag(r_wind_speed) RESULT (r_wind_drag)

    IMPLICIT NONE

    REAL (KIND=GRID_SR), INTENT(in)   :: r_wind_speed
    REAL (KIND=GRID_SR)               :: r_wind_drag

    IF (r_wind_speed <= 11._GRID_SR) THEN
      r_wind_drag = 1.2_GRID_SR
    ELSE IF ((r_wind_speed > 11._GRID_SR).AND.(r_wind_speed <= 25._GRID_SR)) THEN
      r_wind_drag = 0.49_GRID_SR + 0.065_GRID_SR * r_wind_speed
    ELSE
      r_wind_drag = 0.49_GRID_SR + 0.065_GRID_SR * 25._GRID_SR
    ENDIF

    r_wind_drag = r_wind_drag * 1000.0_GRID_SR

  END FUNCTION wind_drag
!*******************************************************************************
END MODULE DG_storm_drag
