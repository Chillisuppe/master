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
! DESCRIPTION of eqn_param
!> structure for equation specific parameters
  TYPE eqn_param
    SEQUENCE
    REAL (KIND = GRID_SR)     :: r_visc          !< viscosity
    REAL (KIND = GRID_SR)     :: r_gamma         !< bottom friction
    REAL (KIND = GRID_SR)     :: r_gammat        !< wind friction
    REAL (KIND = GRID_SR)     :: r_rho           !< density
    REAL (KIND = GRID_SR)     :: r_depth         !< background state of fluid depth
    REAL (KIND = GRID_SR)     :: r_wettol        !< cut off tolerance for wetting and drying
  END TYPE eqn_param
  TYPE (eqn_param)            :: p_equationparam

!-------------------------------------------------------------------------------
! DESCRIPTION of rt_eqninfo
!> structure for equation specific runtime information
  TYPE rt_eqninfo
    REAL (KIND = GRID_SR)     :: r_velmax
    REAL (KIND = GRID_SR)     :: r_sshmax
    REAL (KIND = GRID_SR)     :: r_sshmin
  END TYPE rt_eqninfo
  TYPE (rt_eqninfo)           :: p_rteqninfo

  CONTAINS
!*******************************************************************************
! DESCRIPTION of [SUBROUTINE io_initparams]:
!> @brief initializes structure p_equationparam with default values
!
  SUBROUTINE ioeqn_initparams

    IMPLICIT NONE

    p_equationparam%r_visc   =  0.0_GRID_SR
    p_equationparam%r_gamma  =  0.0_GRID_SR
    p_equationparam%r_gammat =  0.0_GRID_SR
    p_equationparam%r_rho    = 1000.0_GRID_SR
    p_equationparam%r_depth  =  0.0_GRID_SR
    p_equationparam%r_wettol = 0.00000001_GRID_SR

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
    IF(c_param(1:14) == 'VISCOSITY     ') THEN
      read(i_iounit,*) p_equationparam%r_visc
    ELSE IF(c_param(1:14) == 'BOTTOM_FRIC   ') THEN
      read(i_iounit,*) p_equationparam%r_gamma
    ELSE IF(c_param(1:14) == 'WIND_FRIC     ') THEN
      read(i_iounit,*) p_equationparam%r_gammat
    ELSE IF(c_param(1:14) == 'FLUID_DENSITY ') THEN
      read(i_iounit,*) p_equationparam%r_rho
    ELSE IF(c_param(1:14) == 'BATHY_PARAM_DE') THEN
      read(i_iounit,*) p_equationparam%r_depth
    ELSE IF(c_param(1:14) == 'WETDRY_TOL    ') THEN
      read(i_iounit,*) p_equationparam%r_wettol
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

    WRITE(i_iounit,'(A,/,e12.4)') 'VISCOSITY'        , p_equationparam%r_visc
    WRITE(i_iounit,'(A,/,e12.4)') 'BOTTOM_FRIC'      , p_equationparam%r_gamma
    WRITE(i_iounit,'(A,/,e12.4)') 'WIND_FRIC'        , p_equationparam%r_gammat
    WRITE(i_iounit,'(A,/,e12.4)') 'FLUID_DENSITY'    , p_equationparam%r_rho
    WRITE(i_iounit,'(A,/,e12.4)') 'BATHY_PARAM_DEPTH', p_equationparam%r_depth
    WRITE(i_iounit,'(A,/,e12.4)') 'WETDRY_TOL'       , p_equationparam%r_wettol

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
    WRITE(i_iounit,'(A,e12.4,A)') ' +++++ Max. velocity                      ', p_rteqninfo%r_velmax, ' +++++'
    WRITE(i_iounit,'(A,e12.4,A)') ' +++++ Max. ssh (dev. from r_depth)       ', p_rteqninfo%r_sshmax, ' +++++'
    WRITE(i_iounit,'(A,e12.4,A)') ' +++++ Min. ssh (dev. from r_depth)       ', p_rteqninfo%r_sshmin, ' +++++'

  END SUBROUTINE ioeqn_runtimeinfo

END MODULE IO_equation
