!*******************************************************************************
!
!> @file  MISC_diag.F90
!> @brief contains module MISC_diag
!
!*******************************************************************************
!
! VERSION(S):
! 1. original version   j. behrens  8/2009
!
!*******************************************************************************
! MODULE DESCRIPTION:
!> @brief contains the diagnostics function
!
MODULE MISC_diag
  USE FLASH_parameters
  USE GRID_api
  USE MISC_utils
  USE MISC_timing

  INTEGER (KIND = GRID_SI)                              :: i_numdofs, i_iodiag

  PRIVATE
  PUBLIC  :: diag_initialize, diag_quit, diag_syncdisk, diag_diagnostics

  CONTAINS
!*******************************************************************************
! DESCRIPTION of [SUBROUTINE diag_initialize]:
!> @brief initializes diagnostics to the computed solution
!>
!> @param[in]       p_ghand         grid handling data structure
!> @param[in,out]   p_param         control structure for global parameters
!> @param[in]       i_faceunknowns  number of DOFs per element
!> @param[in]       r_einvvand      inverse Vandermonde matrix
!
  SUBROUTINE diag_initialize(p_ghand, p_param, i_faceunknowns, r_einvvand)

    IMPLICIT NONE

    TYPE (grid_handle), INTENT(in)                                :: p_ghand
    TYPE (control_struct), INTENT(inout)                          :: p_param
    INTEGER (KIND = GRID_SI), INTENT(in)                          :: i_faceunknowns
    REAL (KIND = GRID_SR), DIMENSION(:,:), INTENT(in)             :: r_einvvand

!--- local declarations
    CHARACTER (len=32)                                            :: c_file
    INTEGER (KIND = GRID_SI)                                      :: i_fst, i_tmp

    i_numdofs = i_faceunknowns

!--- open file for diagnostic output
    i_tmp = p_param%num%i_experiment
    WRITE(c_file,'(A, A6, I4.4)') TRIM(GRID_parameters%program_name), '_diag.', i_tmp
    OPEN(NEWUNIT=i_iodiag, file=c_file, action='write', form='formatted', iostat=i_fst)

    not_opened: IF(i_fst /= 0) THEN
      CALL grid_error(c_error='[diag_initialize]: could not open general diagnostics file')
    END IF not_opened

    IF(GRID_parameters%iolog > 0) &
      WRITE(GRID_parameters%iolog,*) 'INFO: Filename: ', c_file, ' opened on unit: ', i_iodiag

!--- print some version info and header
    WRITE(i_iodiag,1100) GRID_parameters%program_name, GRID_parameters%version, &
                         GRID_parameters%subversion, GRID_parameters%patchversion
    WRITE(i_iodiag,'(A14)') 'time'

 1100 FORMAT('************************************************', &
             '************************************************',/ &
             '***** PROGRAM: ',a15,61x,                   '*****',/ &
             '***** VERSION: ',i2.2,'.',i2.2,'.',i2.2,68x,'*****',/ &
             '***** Diagnostic output ',67x,              '*****',/ &
             '************************************************', &
             '************************************************',/ &
             '************************************************', &
             '************************************************')

  END SUBROUTINE diag_initialize

!*******************************************************************************
! DESCRIPTION of [SUBROUTINE diag_quit]:
!> @brief exits diagnostics to the computed solution
!>
!> @param[in]       p_ghand     grid handling data structure
!> @param[in,out]   p_param     control structure for global parameters
!
  SUBROUTINE diag_quit(p_ghand, p_param)

    IMPLICIT NONE

    TYPE (grid_handle), INTENT(in)                                :: p_ghand
    TYPE (control_struct), INTENT(inout)                          :: p_param

    WRITE(i_iodiag, 1200)

!--- close diagnostic output file
    CLOSE(i_iodiag)
    IF(GRID_parameters%iolog > 0) &
      write(GRID_parameters%iolog,*) 'INFO: Closed file on unit: ', i_iodiag

 1200 FORMAT('************************************************', &
             '************************************************')

  END SUBROUTINE diag_quit

!*******************************************************************************
! DESCRIPTION of [SUBROUTINE diag_syncdisk]:
!> @brief synchronizes diagnostics file to disk
!
  SUBROUTINE diag_syncdisk

    IMPLICIT NONE

!--- local declarations
    INTEGER (KIND = GRID_SI)                                      :: i_fst

!--- flush and sync
    FLUSH(i_iodiag)
    i_fst = fsync(FNUM(i_iodiag))

    IF (i_fst /= 0) CALL grid_error(c_error='[diag_syncdisk]: Error calling FSYNC')

  END SUBROUTINE diag_syncdisk

!*******************************************************************************
! DESCRIPTION of [SUBROUTINE diag_diagnostics]:
!> @brief applies diagnostics to the computed solution
!>
!> @param[in]       p_ghand         grid handling data structure
!> @param[in,out]   p_param         control structure for global parameters
!> @param[in]       p_tinfo         timestep info structure
!> @param[in]       i_elmtnodes     global element-node relation matrix
!> @param[in]       i_elmtdofs      global element-dof relation matrix
!> @param           r_coonod        coordinates for each grid node
!> @param           r_coodof        coordinates for each dof in the grid
!> @param           r_Q             discrete solution vector
!> @param           r_S             discrete fields used in source terms
!>
!> @note ON OUTPUT: plot output
!
  SUBROUTINE diag_diagnostics(p_ghand, p_param, p_tinfo, i_elmtnodes, i_elmtdofs, &
                              r_coonod, r_coodof, r_Q, r_S)

    IMPLICIT NONE

    TYPE (grid_handle), INTENT(in)                                :: p_ghand
    TYPE (control_struct), INTENT(inout)                          :: p_param
    TYPE (rt_info), INTENT(in)                                    :: p_tinfo
    INTEGER (KIND = GRID_SI), INTENT(in), DIMENSION(:,:)          :: i_elmtnodes, i_elmtdofs
    REAL (KIND = GRID_SR), INTENT(in), DIMENSION(:,:)             :: r_coonod, r_coodof, r_Q, r_S

!--- local declarations

    WRITE(i_iodiag,'(e14.6)') p_tinfo%r_modeltime

  END SUBROUTINE diag_diagnostics

!*******************************************************************************
END MODULE MISC_diag
