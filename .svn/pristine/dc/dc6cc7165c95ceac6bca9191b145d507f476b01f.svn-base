!*******************************************************************************
!
! PROJECT:
!   StormFlash2d
!
! NAME:
!   TAM
! FUNCTION:
!   main program (driver routine for the adaptive scheme)
! SYNTAX:
!   TAM [options]
!
! LIBRARIES:
!   USEs fortran 90 modules:
!     FLASH_parameters
!     IO_utils
!     GRID_api
!     ADV_fvm
! REFERENCES:
!   this is based on TsunaFlash, but all different versions merged
! VERSION(S):
!   1. original version         j. behrens  11/1996
!   2. tidied up a little       j. behrens  07/1997
!   3. new control struct       j. behrens  12/1997
!   4. compliant to amatos 1.0  j. behrens  12/2000
!   5. compliant to amatos 1.2  j. behrens  03/2002
!   6. added visnetplot         f. klaschka 12/2003
!   7. first version of TAM     j. behrens  05/2009
!   8. removed visnetplot       s. vater    09/2012
!
!*******************************************************************************

PROGRAM TAM

!--- modules
  USE GRID_api
  USE FLASH_parameters
  USE IO_utils
  USE ADV_dg
  USE DG_equation, ONLY : equation_initialize, io_initialize

  IMPLICIT NONE

!--- variable declarations
  TYPE (control_struct)       :: p_contr

!--- set FLASH description in global datastruct
  GRID_parameters%program_name = 'StormFlash2d                                    '
  GRID_parameters%version      = 1
  GRID_parameters%subversion   = 0
  GRID_parameters%patchversion = 0
  GRID_parameters%datemonth    = 12
  GRID_parameters%dateyear     = 2011

!--- read command line options
  CALL io_getcmdline(p_contr)

!--- read user input
  CALL io_initparams(p_contr)

!--- read user input
  CALL io_getbatchinput(p_contr)

!--- register the 1 fem type
  IF(p_contr%io%i_numsigs == 1) THEN
    FEM_DG = grid_registerfemtype(p_contr%io%c_sigfiles(1))
  ELSE
    CALL grid_error(c_error='[StormFlash2d]: Inconsistent number of FEM types')
  END IF

!--- initialize grid generator
  IF(p_contr%cmd%l_output) THEN
    IF(p_contr%cmd%l_logging) THEN
      CALL grid_initialize(i_output=i_redirout, i_logging=i_redirlog, &
                           i_initvals=1_GRID_SI, l_preregvars=.FALSE.)
    ELSE
      CALL grid_initialize(i_output=i_redirout, &
                           i_initvals=1_GRID_SI, l_preregvars=.FALSE.)
    END IF
  ELSE
    IF(p_contr%cmd%l_logging) THEN
      CALL grid_initialize(i_logging=i_redirlog, &
                           i_initvals=1_GRID_SI, l_preregvars=.FALSE.)
    ELSE
      CALL grid_initialize(i_initvals=1_GRID_SI, l_preregvars=.FALSE.)
    END IF
  END IF

!--- register new variables
  IF(GRID_dimension /= 2) &
    CALL grid_error(c_error='[StormFlash2d] Inconsistent dimension in amatos, we require GRID_dimension=2')

  CALL equation_initialize(FEM_DG)

!--- initialize general IO
  CALL io_initialize()

!--- print global parameters
  CALL io_putparameters(p_contr)

!--- initialize grid, set up initial conditions etc.
  CALL fvm_initialize(p_grid, p_contr)

!--- call the (major) routine for timestepping
  CALL fvm_timestepping(p_grid, p_contr)

!--- terminate the simulation (gracefully free memory, terminate grid, etc.)
  CALL fvm_finish(p_grid, p_contr)

!--- terminate grid generator
  CALL grid_terminate

  STOP

!*******************************************************************************
END PROGRAM TAM
