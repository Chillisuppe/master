!*******************************************************************************
!
!> @file  TEST_eqsource_okada.f90
!> @brief program to test the eqsource_okada module
!
!*******************************************************************************
PROGRAM TEST_eqsource_okada

  USE GRID_api
  USE MISC_eqsource

  IMPLICIT NONE

  REAL (KIND = GRID_SR)                                 :: r_X, r_Y, r_length, r_width
  REAL (KIND = GRID_SR)                                 :: r_error, r_okadalift, r_uplift
  REAL (KIND = GRID_SR)                                 :: r_okadastrike = -2.747E-3_GRID_SR, &
                                                           r_okadadip    = -3.564E-2_GRID_SR, &
                                                           r_tol         =  1.0E-5_GRID_SR
  REAL (KIND = GRID_SR), PARAMETER                      :: r_earth       =  6.3675E6_GRID_SR
  REAL (KIND = GRID_SR)                                 :: DEG2RAD, LAT2METER

!--- compute value for pi since no grid_initialize is called
  GRID_PI = 4.0_GRID_SR * ATAN(1.0_GRID_SR)

!--- initialize some constants and calculate test for case 2 in the Okada paper 1985
  r_X      = 2.0_GRID_SR
  r_Y      = 3.0_GRID_SR
  r_length = 3.0_GRID_SR
  r_width  = 2.0_GRID_SR

  DEG2RAD   = GRID_PI / 180.0_GRID_SR
  LAT2METER = r_earth * DEG2RAD

!--- switch X to be at the center of fault
  r_X = r_X - r_length / 2.0_GRID_SR

!--- initialize source and compute uplift
  CALL initialize_source('Okada_Test.dat', &
                         [0._GRID_SR, r_X] / LAT2METER, &
                         [0._GRID_SR, r_Y] / LAT2METER, &
                         2_GRID_SI, 2_GRID_SI, 0_GRID_SI)

  r_uplift    = compute_uplift( [r_X, r_Y] / LAT2METER, 0.0_GRID_SR)
  r_okadalift = r_okadadip + r_okadastrike
  r_error     = ABS(r_uplift - r_okadalift)

  WRITE(*,*)
  WRITE(*,'(A,E12.5)') 'The calculated uplift is:          ', r_uplift
  WRITE(*,*)
  WRITE(*,'(A)')       'The values from the Okada paper (1985) are:'
  WRITE(*,'(A,E12.5)') 'strike-slip:                       ', r_okadastrike
  WRITE(*,'(A,E12.5)') 'dip-slip:                          ', r_okadadip
  WRITE(*,'(A,E12.5)') 'and the resulting uplift is:       ', r_okadalift
  WRITE(*,*)
  WRITE(*,'(A,E12.5)') 'The given test tolerance is:      ', r_tol

!--- check if test is passed (with tolerance r_tol)
  IF (ABS(r_error) < r_tol) THEN
    WRITE(*,'(A,E12.5)') 'The test passed with an error of: ', r_error
  ELSE
    WRITE(*,'(A,E12.5)') 'The test failed with an error of: ', r_error
  END IF

END PROGRAM TEST_eqsource_okada
