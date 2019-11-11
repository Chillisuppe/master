!*******************************************************************************
!
!> @file  IO_equation.f90
!> @brief contains module IO_equation
!
!*******************************************************************************
!
! VERSION(S):
!  1. original version                               s. vater     09/2018
!*******************************************************************************
! MODULE DESCRIPTION:
!> @brief equation specific input/output routines
!
MODULE IO_equation

  USE FLASH_parameters

  IMPLICIT NONE

  PRIVATE
  PUBLIC  :: p_equationparam, p_rteqninfo, ioeqn_initparams, &
             ioeqn_checkreadinput, ioeqn_putinputfile, ioeqn_runtimeinfo

!-------------------------------------------------------------------------------
! DESCRIPTION of num_param
!> structure for equation specific parameters
  TYPE eqn_param
    SEQUENCE
    REAL (KIND = GRID_SR)     :: r_gamma            !< specific heat ratio
    REAL (KIND = GRID_SR)     :: r_pref             !< reference pressure of gas
    REAL (KIND = GRID_SR)     :: r_rgas             !< spefic gas constant
    REAL (KIND = GRID_SR)     :: r_gravityswitch    !< handle for gravity source term
  END TYPE eqn_param
  TYPE (eqn_param)            :: p_equationparam

!-------------------------------------------------------------------------------
! DESCRIPTION of rt_eqninfo
!> structure for equation specific runtime information
  TYPE rt_eqninfo
    REAL (KIND = GRID_SR)     :: r_velmax
    REAL (KIND = GRID_SR)     :: r_densmax
    REAL (KIND = GRID_SR)     :: r_densmin
    REAL (KIND = GRID_SR)     :: r_pressmax
    REAL (KIND = GRID_SR)     :: r_pressmin
  END TYPE rt_eqninfo
  TYPE (rt_eqninfo)           :: p_rteqninfo

  CONTAINS
!*******************************************************************************
! DESCRIPTION of [SUBROUTINE io_initparams]:
!> @brief initializes structure p_equationparam with default values
!
  SUBROUTINE ioeqn_initparams

    IMPLICIT NONE

    p_equationparam%r_gamma         = 1.4_GRID_SR
    p_equationparam%r_pref          = 101325.0_GRID_SR
    p_equationparam%r_rgas          = 287.0_GRID_SR
    p_equationparam%r_gravityswitch = 1.0_GRID_SR

  END SUBROUTINE ioeqn_initparams

!*******************************************************************************
! DESCRIPTION of [FUNCTION ioeqn_checkreadinput]:
!> @brief reads user input from file
!>
!> @param[in]       c_param           parameter string to be checked
!> @param[in]       i_iounit          input unit of parameters file
!> @return                            true if parameter string is an equation parameter
!
  FUNCTION ioeqn_checkreadinput(c_param, i_iounit) RESULT (l_iseqnparam)

    IMPLICIT NONE

    CHARACTER (len=*), INTENT(in)                             :: c_param
    INTEGER (KIND = GRID_SI), INTENT(in)                      :: i_iounit
    LOGICAL                                                   :: l_iseqnparam

!--- local declarations

    l_iseqnparam = .TRUE.

!--- check if string is an equation specific parameter
    IF(c_param(1:14) == 'GAMMA         ') THEN
      read(i_iounit,*) p_equationparam%r_gamma
    ELSE IF(c_param(1:14) == 'REF_PRESSURE  ') THEN
      read(i_iounit,*) p_equationparam%r_pref
    ELSE IF(c_param(1:14) == 'GAS_CONSTANT  ') THEN
      read(i_iounit,*) p_equationparam%r_rgas
    ELSE IF(c_param(1:14) == 'GRAVITY_SWITCH') THEN
      read(i_iounit,*) p_equationparam%r_gravityswitch
    ELSE
      l_iseqnparam = .FALSE.
    END IF

  END FUNCTION ioeqn_checkreadinput

!*******************************************************************************
! DESCRIPTION of [SUBROUTINE ioeqn_putinputfile]:
!> @brief writes equation specific parameters to parameter file
!>
!> @param[in]       i_iounit          output unit for parameters file
!
  SUBROUTINE ioeqn_putinputfile(i_iounit)

    IMPLICIT NONE

    INTEGER (KIND = GRID_SI), INTENT(in)                      :: i_iounit

!--- local declarations

    WRITE(i_iounit,'(A,/,e12.4)') 'GAMMA      '      , p_equationparam%r_gamma
    WRITE(i_iounit,'(A,/,e12.4)') 'REF_PRESSURE'     , p_equationparam%r_pref
    WRITE(i_iounit,'(A,/,e12.4)') 'GAS_CONSTANT'     , p_equationparam%r_rgas
    WRITE(i_iounit,'(A,/,e12.4)') 'GRAVITY_SWITCH'   , p_equationparam%r_gravityswitch

  END SUBROUTINE ioeqn_putinputfile

!*******************************************************************************
! DESCRIPTION of [SUBROUTINE ioeqn_runtimeinfo]:
!> @brief prints some information on the current equation specific run time status
!>
!> @param[in]   i_iounit      IO unit where output is sent
!
  SUBROUTINE ioeqn_runtimeinfo(i_iounit)

    IMPLICIT NONE

    INTEGER (KIND = GRID_SI), INTENT(in)                      :: i_iounit

!--- local declarations

!--- output to i_iounit
    WRITE(i_iounit,'(A)')         ' +++++ ----- ----- ----- ----- ----- ----- ----- ----- +++++'
    WRITE(i_iounit,'(A,e12.4,A)') ' +++++ Max. velocity                      ', p_rteqninfo%r_velmax,   ' +++++'
    WRITE(i_iounit,'(A,e12.4,A)') ' +++++ Max. density                       ', p_rteqninfo%r_densmax,  ' +++++'
    WRITE(i_iounit,'(A,e12.4,A)') ' +++++ Min. density                       ', p_rteqninfo%r_densmin,  ' +++++'
    WRITE(i_iounit,'(A,e12.4,A)') ' +++++ Max. pressure                      ', p_rteqninfo%r_pressmax, ' +++++'
    WRITE(i_iounit,'(A,e12.4,A)') ' +++++ Min. pressure                      ', p_rteqninfo%r_pressmin, ' +++++'

  END SUBROUTINE ioeqn_runtimeinfo

END MODULE IO_equation
