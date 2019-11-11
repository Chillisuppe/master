!*******************************************************************************
!
!> @file  MISC_diag.F90
!> @brief contains module MISC_diag
!
!*******************************************************************************
! MODULE DESCRIPTION:
!> @brief contains the diagnostics functionality to compute exact gauge timeseries
!>        at given points. These points are defined in the file gaugedata.dat in
!>        following format:
!>        - a line 'NUM_GAUGES' is followed by a line containing the number of
!>          gauges
!>        - a line 'GAUGE_COORDINATES' is followed by by NUM_GAUGES lines
!>          containing the (x,y) or (lon,lat) coordinates of the gauges (given
!>          by two numbers separated by spaces
!>        - a line 'GAUGE_NAMES' is followed by by NUM_GAUGES lines containing
!>          the gauge names used for the header in the diagnostics file
!
MODULE MISC_diag

  USE FLASH_parameters
  USE GRID_api
  USE MISC_timing
  USE MISC_utils
  USE DG_utils, ONLY: Vandermonde2D, bary

  INTEGER (KIND = GRID_SI)                              :: i_numdofs, i_iodiag, i_ngauges
  REAL (KIND = GRID_SR), DIMENSION(:,:), ALLOCATABLE    :: r_einvvander, r_coogauges
  CHARACTER (len=32), DIMENSION(:), ALLOCATABLE         :: c_gaugenames

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
    INTEGER (KIND = GRID_SI)                                      :: i_fst, i_tmp, i_alct, i_cnt

    i_numdofs = i_faceunknowns

    ALLOCATE(r_einvvander(i_numdofs, i_numdofs), stat=i_alct)
    IF(i_alct /= 0) &
      CALL grid_error(c_error='[diag_initialize]: Could not allocate variables')
    r_einvvander = r_einvvand

!--- read gauge data
    CALL read_gauges('gaugedata.dat')

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
    WRITE(i_iodiag,'(A14,1x)', advance='no') 'time                '
    DO i_cnt = 1, i_ngauges
      WRITE(i_iodiag,'(A14,1x)', advance='no') c_gaugenames(i_cnt)
    END DO

    WRITE(i_iodiag,*)

 1100 FORMAT('*********************************************', &
             '*********************************************',/ &
             '***** PROGRAM: ',a15,55x,                   '*****',/ &
             '***** VERSION: ',i2.2,'.',i2.2,'.',i2.2,62x,'*****',/ &
             '***** Diagnostic output ',61x,              '*****',/ &
             '*********************************************', &
             '*********************************************')

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

    IF(ALLOCATED(r_einvvander)) DEALLOCATE(r_einvvander)
    IF(ALLOCATED(r_coogauges))  DEALLOCATE(r_coogauges)
    IF(ALLOCATED(c_gaugenames)) DEALLOCATE(c_gaugenames)

    WRITE(i_iodiag, 1200)

!--- close diagnostic output file
    CLOSE(i_iodiag)
    IF(GRID_parameters%iolog > 0) &
      write(GRID_parameters%iolog,*) 'INFO: Closed file on unit: ', i_iodiag

 1200 FORMAT('*********************************************', &
             '*********************************************')

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
!> @param[in]       p_ghand     grid handling data structure
!> @param[in,out]   p_param     control structure for global parameters
!> @param[in]       p_tinfo     timestep info structure
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
    INTEGER (KIND = GRID_SI)                                      :: i_cnt
    INTEGER (KIND = GRID_SI), DIMENSION(i_ngauges)                :: i_gaugeelmts
    REAL (KIND = GRID_SR)                                         :: r_valgauge
    REAL (KIND = GRID_SR), DIMENSION(GRID_elementnodes,i_ngauges) :: r_bary
    REAL (KIND = GRID_SR), DIMENSION(i_ngauges,i_numdofs)         :: r_V, r_psi

    WRITE(i_iodiag,'(e14.7,1x)', advance='no') p_tinfo%r_modeltime

    DO i_cnt = 1, i_ngauges
      i_gaugeelmts(i_cnt) = grid_findelmt(r_coogauges(:, i_cnt), p_ghand, l_relative=.TRUE.)
      r_bary(:,i_cnt)     = bary(r_coogauges(:, i_cnt), &
                                 r_coonod(:,i_elmtnodes(:,i_gaugeelmts(i_cnt))))
    END DO

    CALL Vandermonde2D(GRID_femtypes%p_type(FEM_DG)%sig%i_degree, r_bary, r_V)
    r_psi = MATMUL(r_V, r_einvvander)

    DO i_cnt = 1, i_ngauges
      r_valgauge = DOT_PRODUCT(r_psi(i_cnt,:), &
        (r_Q(1, i_elmtdofs(:,i_gaugeelmts(i_cnt)))+r_S(1,i_elmtdofs(:,i_gaugeelmts(i_cnt))))/ GRID_GRAV)
      WRITE(i_iodiag,'(e14.7,1x)', advance='no') r_valgauge
    END DO

    WRITE(i_iodiag,*)

  END SUBROUTINE diag_diagnostics

!*******************************************************************************
! DESCRIPTION of [SUBROUTINE read_gauges]:
!> @brief short description
!
  SUBROUTINE read_gauges(c_filename)

    IMPLICIT NONE

    CHARACTER (len=*)                                     :: c_filename

!--- local declarations
    INTEGER (KIND = GRID_SI)                              :: i_iofil
    INTEGER (KIND = GRID_SI), PARAMETER                   :: i_ioerr = 0
    INTEGER (KIND = GRID_SI)                              :: i_iost, i_ioend, i_alct, i_cnt
    CHARACTER (len=80)                                    :: a_filrow

!--- open file for gauges information
    open(NEWUNIT=i_iofil, file=c_filename, status='OLD', action='READ', iostat=i_iost)
    file_notopen: IF(i_iost /= 0) THEN
      write(i_ioerr,*) 'ERROR: Filename: ', c_filename
      IF(GRID_parameters%iolog > 0) &
        write(GRID_parameters%iolog,*) 'ERROR: Filename: ', c_filename
      CALL grid_error(c_error='[read_gauges]: Could not open file for gauge coordinates')
    ELSE file_notopen
      IF(GRID_parameters%iolog > 0) THEN
        write(GRID_parameters%iolog,*) 'INFO: Filename: ', c_filename, ' opened on unit: ', i_iofil
      END IF
    END IF file_notopen

!--- read line by line

    read_loop: DO
      READ(i_iofil,'(a80)',iostat=i_ioend) a_filrow

!--- if file ended

      file_end: IF(i_ioend /= 0) THEN
        CLOSE(i_iofil)
        IF(GRID_parameters%iolog > 0) &
          WRITE(GRID_parameters%iolog,*) 'INFO: Closed file on unit: ', i_iofil
        EXIT
      ELSE file_end

!--- decide what to DO with line according to first character

        comment_line: IF(a_filrow(1:1) == '#' .or. a_filrow(1:1) == '!') THEN
          CYCLE read_loop
        ELSE IF (a_filrow(1:10) == 'NUM_GAUGES') THEN
          READ(i_iofil,*) i_ngauges
          ALLOCATE(r_coogauges(2, i_ngauges), c_gaugenames(i_ngauges), STAT=i_alct)
            IF(i_alct /= 0) CALL grid_error(c_error='[read_gauges]: could not allocate arrays')
        ELSE IF (a_filrow(1:10) == 'GAUGE_COOR') THEN
          DO i_cnt=1,i_ngauges
            READ(i_iofil,*) r_coogauges(:,i_cnt)
          END DO
        ELSE IF (a_filrow(1:10) == 'GAUGE_NAME') THEN
          DO i_cnt=1,i_ngauges
            READ(i_iofil,'(a)') c_gaugenames(i_cnt)
          END DO
        ELSE comment_line
          CALL grid_error(c_error='[read_gauges]: Dont know what to do with this parameter file')
        END IF comment_line
      END IF file_end
    END DO read_loop

  END SUBROUTINE read_gauges

!*******************************************************************************
END MODULE MISC_diag
