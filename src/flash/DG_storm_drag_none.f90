!*******************************************************************************
!
!> @file  DG_storm_drag_none.f90
!> @brief contains module DG_storm_drag
!
!*******************************************************************************

MODULE DG_storm_drag

  USE GRID_api
  PUBLIC :: wind_drag

  CONTAINS

!*******************************************************************************
! DESCRIPTION of [FUNCTION wind_drag]:
!> @brief This is a dummy module. It sets the wind drag to zero (for testing purposes only)
!>
!> @param[in]     r_wind_speed    magnitude of the wind in the cell
!> @return                        coefficient of drag
!>
  FUNCTION wind_drag(r_wind_speed) RESULT (r_wind_drag)

    IMPLICIT NONE

    REAL (KIND=GRID_SR), INTENT(in)   :: r_wind_speed
    REAL (KIND=GRID_SR)               :: r_wind_drag

    r_wind_drag = 0.0_GRID_SR

  END FUNCTION wind_drag

!*******************************************************************************
END MODULE DG_storm_drag
