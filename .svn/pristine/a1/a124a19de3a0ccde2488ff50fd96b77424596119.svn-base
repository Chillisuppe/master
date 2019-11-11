!*******************************************************************************
!
!> @file DG_storm_holland.f90
!> @brief contains module DG_storm_holland adapted from hurricane.f90 by S. Gopalakrisnan, 2011
!< based on Greg J. Holland "An Analytic Model of Wind and Pressure Profiles in Hurricanes"
!< Monthly Weather Review, Volume 108, pp 1212-1218, 1980
!
!*******************************************************************************

MODULE DG_storm_holland

  USE GRID_api
  USE DG_storm_drag

  PUBLIC  :: get_wind_holland, get_storm_position

  CONTAINS

!*******************************************************************************
    SUBROUTINE get_wind_holland(r_windfield, r_pos_storm, r_Pn, r_Pc, r_rho_air, r_A, r_B, &
                                r_x, r_y, r_time, r_tramp)

      IMPLICIT NONE

      REAL (KIND=GRID_SR), INTENT(in)                   :: r_Pn, r_Pc, r_rho_air, r_A, r_B
      REAL (KIND=GRID_SR), INTENT(in)                   :: r_x, r_y
      REAL (KIND=GRID_SR), INTENT(in)                   :: r_time, r_tramp
      REAL (KIND=GRID_SR), DIMENSION(2), INTENT(in)     :: r_pos_storm
      REAL (KIND=GRID_SR), DIMENSION(2), INTENT(inout)  :: r_windfield

      !-- Local declarations
      REAL (KIND=GRID_SR), DIMENSION(2)                 :: r_wind

      REAL (KIND=GRID_SR)                               :: r_C, r_Cd, r_f
      REAL (KIND=GRID_SR)                               :: r_w, r_r

!-- INITIALIZE CONSTANTS (TAKEN FROM HOLLAND'S PAPER/ FOR HURRICANE TRACY)
    r_Cd      =    0.025
    r_f       =    0.0 !CORIOLIS IS NEGLIGABLE

!-- Parameter constants
    r_C         = (r_A * r_B * (r_Pn-r_Pc) / r_rho_air) * 100.0_GRID_SR
    r_windfield = 0._GRID_SR

!-- Compute radius
    r_r = sqrt((r_x-r_pos_storm(1))**2 + &
               (r_y-r_pos_storm(2))**2 ) ! *10 !!*1000._GRID_SR

    !-- In the inner part, there are no velocities
    IF ((abs(r_r) < 10E-8)) THEN ! .OR. (abs(r_r) > 50000)) THEN !10E-1 for traveling storms
       r_wind(:)      = 0._GRID_SR
    ELSE

       !-- Formula (4) of Holland's paper for the gradient wind
       r_w = sqrt(r_C * exp(-r_A/r_r**r_B) / r_r**r_B + r_r**2 * r_f**2 / 4.0_GRID_SR) &
            - r_r * r_f / 2.0_GRID_SR
       !-- For earlier times - Initialization period 2 hours
        IF ( r_time <= r_tramp ) THEN
           r_w=r_w*exp(-(((r_time - r_tramp)/(r_tramp*0.45_GRID_SR))**2))
        END IF

        !-- Compute wind in (x,y)-/ normal- direction
        r_r = r_r !/10 !*1000._GRID_SR
        r_wind(1) = -r_w * r_y / r_r
        r_wind(2) =  r_w * r_x / r_r

        !-- Compute wind drag coefficient using private function
        r_Cd = wind_drag(r_w)

        !-- Determine wind field = drag*densitiy*velo*|velo|, formula (9) in Weisberg & Zheng Paper
        r_windfield(:) = r_Cd*r_rho_air*r_w*r_wind(:) / 1000.0_GRID_SR

     ENDIF

    END SUBROUTINE get_wind_holland

!*******************************************************************************
  SUBROUTINE get_storm_position(r_time, r_tstart, r_tfinal, r_tramp, r_pos_storm, r_pos_start, r_pos_final)

    IMPLICIT NONE

    REAL (KIND=GRID_SR), INTENT(in)                    :: r_time, r_tstart, r_tfinal
    REAL (KIND=GRID_SR), INTENT(in)                    :: r_tramp

    REAL (KIND=GRID_SR), DIMENSION(2), INTENT(inout)   :: r_pos_storm
    REAL (KIND=GRID_SR), DIMENSION(2), INTENT(in)      :: r_pos_start, r_pos_final

    !-- Local declarations
    REAL (KIND=GRID_SR)                                :: r_tnorm
    REAL (KIND=GRID_SR)                                :: r_alpha


!-- Keep hurricane in position 0 for the first two hours to ramp up speed smoothly
    IF (r_time <= r_tstart + r_tramp) THEN
       r_pos_storm(:)  = r_pos_start(:)
    ELSE

!-- set time counter back from 2h to 0s
      r_tnorm = r_time - (r_tstart + r_tramp)

      r_alpha     = r_tnorm / (r_tfinal - r_tstart)

      !write(*,*) r_alpha

      r_pos_storm(1) = r_pos_start(1) + r_alpha * (r_pos_final(1) - r_pos_start(1))
      r_pos_storm(2) = r_pos_start(2) + r_alpha * (r_pos_final(2) - r_pos_start(2))

    END IF !before 2 hours

    !write(*,*) r_pos_storm

  END SUBROUTINE get_storm_position

!*******************************************************************************
END MODULE DG_storm_holland
