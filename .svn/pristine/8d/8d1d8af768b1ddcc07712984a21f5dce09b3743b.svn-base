!*******************************************************************************
!
!> @file  FLASH_parameters.f90
!> @brief contains module FLASH_parameters
!
!*******************************************************************************
!
! VERSION(S):
!  1. original version          j. behrens  12/2000
!  2. changed for TsunaFlash    j. behrens  01/2006
!  3. cleaned, uplift added     j. behrens  04/2009
!  4. added l_para              n.beisiegel 01/2011
!  5. added logical for linear
!     equations and limiter     n.beisiegel 08/2011
!  6. removed logical for linear
!     equations                 s. vater    07/2013
!  7. added structure for
!     generalized IO routines   m. baensch  07/2018
!*******************************************************************************
! MODULE DESCRIPTION:
!> @brief defines global control structure
!
MODULE FLASH_parameters
  USE GRID_api

  IMPLICIT NONE

  INTEGER (KIND = GRID_SI), PARAMETER :: io_fillen  = 48
  INTEGER (KIND = GRID_SI), PARAMETER :: i_redirout = 8
  INTEGER (KIND = GRID_SI), PARAMETER :: i_redirlog = 7
  INTEGER, PARAMETER                  :: i_comparlen = 14

  TYPE param_int
    SEQUENCE
    CHARACTER (len=i_comparlen)                     :: c_keyword
    INTEGER                                         :: i_size
    INTEGER (KIND = GRID_SI), POINTER, DIMENSION(:) :: p_value
  END TYPE param_int

  TYPE param_real
    SEQUENCE
    CHARACTER (len=i_comparlen)                   :: c_keyword
    INTEGER                                       :: i_size
    REAL (KIND = GRID_SR), POINTER, DIMENSION(:)  :: p_value
  END TYPE param_real

  TYPE param_char
    SEQUENCE
    CHARACTER (len=i_comparlen) :: c_keyword
    CHARACTER (len=io_fillen)   :: p_value
  END TYPE param_char

  TYPE param_log
    SEQUENCE
    CHARACTER (len=i_comparlen) :: c_keyword
    LOGICAL                     :: p_value
  END TYPE param_log

!-------------------------------------------------------------------------------
! DESCRIPTION of cmdline_param
!> structure for the command line
  TYPE cmdline_param
    SEQUENCE
    CHARACTER (len=io_fillen) :: c_batchfile    !< input file name
    LOGICAL                   :: l_output       !< redirect std output
    LOGICAL                   :: l_logging      !< enable logging (verbose)
  END TYPE cmdline_param

!-------------------------------------------------------------------------------
! DESCRIPTION of io_param
!> structure for the i/o behaviour
  TYPE io_param
    SEQUENCE
    INTEGER (KIND = GRID_SI)  :: i_plotoffset   !< timesteps between plots (constant timestep sizes)
    REAL (KIND = GRID_SR)     :: r_plotint      !< time intervals between plots (variable timestep sizes)
    INTEGER (KIND = GRID_SI)  :: i_saveoffset   !< timesteps between savesets
    INTEGER (KIND = GRID_SI)  :: i_diagoffset   !< timesteps between diagnostics
    INTEGER (KIND = GRID_SI)  :: i_savelast     !< indicator for last step saving
    INTEGER (KIND = GRID_SI)  :: i_numsigs      !< number of signature files
    CHARACTER (len=io_fillen) :: c_triangfile   !< file with initial triangulation
    CHARACTER (len=io_fillen), DIMENSION(:), POINTER :: c_sigfiles !< signature files
    LOGICAL                   :: l_netcdf       !< write NetCDF (offline) output
    LOGICAL                   :: l_para         !< write vtu (offline) output
    LOGICAL                   :: l_diagnostics  !< switch for diagnostics
    LOGICAL                   :: l_solplot      !< switch for analytic solution plot
    INTEGER                   :: i_subtriangpts !< number of equidistant subtriangulation points per edge (netcdf plot)
  END TYPE io_param

!-------------------------------------------------------------------------------
! DESCRIPTION of num_param
!> structure for global numeric and steering parameters
  TYPE num_param
    SEQUENCE
    REAL (KIND = GRID_SR)     :: r_deltatime     !< timestep length [s]
    REAL (KIND = GRID_SR)     :: r_reftolerance  !< tolerance for refinement
    REAL (KIND = GRID_SR)     :: r_crstolerance  !< tolerance for coarsening
    REAL (KIND = GRID_SR)     :: r_refwatermark  !< watermark for refinement
    REAL (KIND = GRID_SR)     :: r_crswatermark  !< watermark for coarsening
    REAL (KIND = GRID_SR)     :: r_cflmax        !< max. allowed cfl number
    INTEGER (KIND = GRID_SI)  :: i_experiment    !< current experiment identification
    INTEGER (KIND = GRID_SI)  :: i_crslevel      !< coarsest requested level
    INTEGER (KIND = GRID_SI)  :: i_reflevel      !< finest requested level
    INTEGER (KIND = GRID_SI)  :: i_frsttimestep  !< first timestep of experiment
    INTEGER (KIND = GRID_SI)  :: i_lasttimestep  !< last timestep of experiment
    REAL (KIND = GRID_SR)     :: r_starttime     !< initial time of experiment
    REAL (KIND = GRID_SR)     :: r_finaltime     !< final time of experiment
    INTEGER (KIND = GRID_SI)  :: i_refneigh      !< how many neighbours will be refined in the region around the fine region
    CHARACTER (len=io_fillen) :: c_timestepping  !< timestepping method
  END TYPE num_param

  TYPE test_param
    SEQUENCE
    TYPE (param_int),  DIMENSION(:), POINTER :: tint
    TYPE (param_real), DIMENSION(:), POINTER :: treal
    TYPE (param_char), DIMENSION(:), POINTER :: tchar
    TYPE (param_log),  DIMENSION(:), POINTER :: tlog
  END TYPE test_param

!-------------------------------------------------------------------------------
! DESCRIPTION of control_struct
!> global control structure
  TYPE control_struct
    SEQUENCE
    TYPE (num_param)          :: num
    TYPE (cmdline_param)      :: cmd
    TYPE (io_param)           :: io
  END TYPE control_struct

!-------------------------------------------------------------------------------
! DESCRIPTION of rt_info
!> structure for runtime information
  TYPE rt_info
    REAL (KIND = GRID_SR)     :: r_modeltime
    REAL (KIND = GRID_SR)     :: r_cputime
    REAL (KIND = GRID_SR)     :: r_timestep
    REAL (KIND = GRID_SR)     :: r_velmax
    REAL (KIND = GRID_SR)     :: r_sshmax
    REAL (KIND = GRID_SR)     :: r_sshmin
    REAL (KIND = GRID_SR)     :: r_cflnumber
    INTEGER (KIND = GRID_SI)  :: i_step
    INTEGER (KIND = GRID_SI)  :: i_adapit
    LOGICAL                   :: l_saved
    LOGICAL                   :: l_ploted
    LOGICAL                   :: l_diaged
  END TYPE rt_info
  TYPE (rt_info)              :: p_timestepinfo

!-------------------------------------------------------------------------------
! DESCRIPTION of io_vars
!> structure for i/o routines containing names of variables, units etc
  TYPE io_vars
    SEQUENCE
    REAL (KIND = GRID_SR)     :: r_factor         !< multiplication factor (e.g 1/GRID_GRAV)
    CHARACTER (len=io_fillen) :: c_varname        !< variable name
    CHARACTER (len=io_fillen) :: c_long_name      !< long name for variable
    CHARACTER (len=io_fillen) :: c_standard_name  !< standard name for variable
    CHARACTER (len=io_fillen) :: c_units          !< units for variable
  END TYPE io_vars

!--- constant for the FEM type
  INTEGER (KIND = GRID_SI)    :: FEM_DG      !< DG FEM type

!*******************************************************************************
END MODULE FLASH_parameters
