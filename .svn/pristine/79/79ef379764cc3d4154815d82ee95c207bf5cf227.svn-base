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
    REAL (KIND = GRID_SR), DIMENSION(GRID_DIMENSION)  :: r_velo !< background velocity
  END TYPE eqn_param
  TYPE (eqn_param)            :: p_equationparam

!-------------------------------------------------------------------------------
! DESCRIPTION of rt_eqninfo
!> structure for equation specific runtime information
  TYPE rt_eqninfo
    REAL (KIND = GRID_SR)     :: r_tracermax
    REAL (KIND = GRID_SR)     :: r_tracermin
  END TYPE rt_eqninfo
  TYPE (rt_eqninfo)           :: p_rteqninfo

  CONTAINS
!*******************************************************************************
! DESCRIPTION of [SUBROUTINE io_initparams]:
!> @brief initializes structure p_equationparam with default values
!
  SUBROUTINE ioeqn_initparams

    IMPLICIT NONE

    p_equationparam%r_velo = [0.1_GRID_SR, 0.1_GRID_SR]

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
    IF(c_param(1:14) == 'VELOCITY      ') THEN
      read(i_iounit,*) p_equationparam%r_velo
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

    WRITE(i_iounit,'(A,/,e12.4)') 'VELOCITY'         , p_equationparam%r_velo

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
    WRITE(i_iounit,'(A,e12.4,A)') ' +++++ Max. of tracer                     ', p_rteqninfo%r_tracermax, ' +++++'
    WRITE(i_iounit,'(A,e12.4,A)') ' +++++ Min. of tracer                     ', p_rteqninfo%r_tracermin, ' +++++'

  END SUBROUTINE ioeqn_runtimeinfo

END MODULE IO_equation
