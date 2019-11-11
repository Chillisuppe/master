!*******************************************************************************
!
!> @file  DG_storm_drag_powell.f90
!> @brief contains module DG_storm_drag
!
!*******************************************************************************

MODULE DG_storm_drag

  USE GRID_api
  PUBLIC :: wind_drag

  CONTAINS

!*******************************************************************************
! DESCRIPTION of [FUNCTION wind_drag]:
!> @brief Calculates the drag coefficient for wind stess depending on the
!>        wind speed based on "New findings on hurricane intensity, wind field
!>        extent, and surface drag coefficient behavior" by Powell (2007)
!>
!> @param[in]     r_wind_speed    magnitude of the wind in the cell
!> @return                        coefficient of drag
!>
  FUNCTION wind_drag(r_wind_speed) RESULT (r_wind_drag)

    IMPLICIT NONE

    REAL (KIND=GRID_SR), INTENT(in)   :: r_wind_speed
    REAL (KIND=GRID_SR)               :: r_wind_drag

    REAL (KIND=GRID_SR)               :: r_drag_l, r_drag_r, r_drag_rear, &
                                         r_weight, r_theta
    REAL (KIND=GRID_SR), DIMENSION(2) :: r_drag

    ! Calculate sector drags
    IF (r_wind_speed <= 15.708) THEN
      r_drag_l = 7.5E-4 + 6.6845E-5 * r_wind_speed
      r_drag_r = r_drag_l
      r_drag_rear = r_drag_l
    ELSEIF (15.708 < r_wind_speed .AND. r_wind_speed <= 18.7) THEN
      r_drag_l = 1.8E-3
      r_drag_r = 7.5E-4 + 6.6845E-5 * r_wind_speed
      r_drag_rear = r_drag_r
    ELSEIF (18.7 < r_wind_speed .AND. r_wind_speed <= 25.0) THEN
      r_drag_l = 1.8E-3
      r_drag_r = 2.0E-3
      r_drag_rear = r_drag_r
    ELSEIF (25.0 < r_wind_speed .AND. r_wind_speed <= 30.0) THEN
      r_drag_l = 1.8E-3 + 5.4E-3 * (r_wind_speed - 25.0)
      r_drag_r = 2.0E-3
      r_drag_rear = r_drag_r
    ELSEIF (30.0 < r_wind_speed .AND. r_wind_speed <= 35.0) THEN
      r_drag_l = 4.5E-3 - 2.33333E-4 * (r_wind_speed - 30.0)
      r_drag_r = 2.0E-3
      r_drag_rear = r_drag_r
    ELSEIF (35.0 < r_wind_speed .AND. r_wind_speed <= 45.0) THEN
      r_drag_l = 4.5E-3 - 2.33333E-4 * (r_wind_speed - 30.0)
      r_drag_r = 2.E-3 + 1.E-4 * (r_wind_speed - 35.0)
      r_drag_rear = 2.E-3 - 1.E-4 * (r_wind_speed - 35.0)
    ELSE
      r_drag_l = 1.0E-3
      r_drag_r = 3.0E-3
      r_drag_rear = r_drag_l
    ENDIF

    ! Calculate weights
    ! Left sector = [240., 20.] - Center = 310
    ! Left Right sector = [310, 85] - Width = 145
    ! Right sector = [ 20., 150.] - Center = 85
    ! Right Rear sector = [85, 195] - Width =
    ! Rear sector = [150., 240.] - 195
    ! Rear-Left sector = [85, 195] - Width =
    ! Left - Right sector
    IF (310.0 < r_theta .AND. r_theta <= 360.0) THEN
      r_weight = (r_theta - 310.0) / 135.0
      r_drag = [r_drag_l, r_drag_r]
      ! Left - Right sector
    ELSEIF (0.0 < r_theta .AND. r_theta <= 85.0) THEN
      r_weight = (r_theta + 50.0) / 135.0
      r_drag = [r_drag_l, r_drag_r]
      ! Right - Rear sector
    ELSEIF (85.0 < r_theta .AND. r_theta <= 195.0) THEN
      r_weight = (r_theta - 85.0) / 110.0
      r_drag = [r_drag_r, r_drag_rear]
      ! Rear - Left sector
    ELSEIF (195.0 < r_theta .AND. r_theta <= 310.0) THEN
      r_weight = (r_theta - 195.0) / 115.0
      r_drag = [r_drag_rear, r_drag_l]
    ENDIF

    r_wind_drag = r_drag(1) * (1.0 - r_weight) + r_drag(2) * r_weight
    r_wind_drag = r_wind_drag * 1000.0_GRID_SR

  END FUNCTION wind_drag

!*******************************************************************************
END MODULE DG_storm_drag
