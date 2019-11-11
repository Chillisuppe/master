!*******************************************************************************
!
!> @file DG_initial_okushiri.f90
!> @brief contains module DG_initial
!
!*******************************************************************************
! MODULE DESCRIPTION:
!> @brief initializes fields with the Okushiri test case
!>        Note that the background fluid depth is set with BATHY_PARAM_DEPTH
!>        in the Parameters file.
!
MODULE DG_initial

  USE GRID_api
  USE FLASH_parameters
  USE DG_equation, ONLY : i_nprogvars, i_nsrcterms, i_valQ, i_valS
  USE IO_equation, ONLY : p_equationparam

  PRIVATE
  PUBLIC  i_ntstintparam, i_ntstrealparam, i_ntstcharparam, i_ntstlogparam, &
          initial_setparam, p_testparam, initialize_testcase, initialvalues_iter, &
          source_update, dg_elmt_solution, boundaryhgt

  TYPE (test_param)                                   :: p_testparam
  INTEGER, PARAMETER                                  :: i_ntstintparam  = 0
  INTEGER, PARAMETER                                  :: i_ntstrealparam = 0
  INTEGER, PARAMETER                                  :: i_ntstcharparam = 2
  INTEGER, PARAMETER                                  :: i_ntstlogparam  = 0

  REAL (KIND = GRID_SR), DIMENSION(393,244)               :: r_bathygrid
  REAL (KIND = GRID_SR), DIMENSION(393)                   :: r_xcoo
  REAL (KIND = GRID_SR), DIMENSION(244)                   :: r_ycoo
  REAL (KIND = GRID_SR), DIMENSION(2, 451)                :: r_bheight


  CONTAINS

!*******************************************************************************
! DESCRIPTION of [SUBROUTINE initial_setparam]:
!> @brief initializes test case specific parameter list
!
  SUBROUTINE initial_setparam

    IMPLICIT NONE

    ! file name for bottom topography
    p_testparam%tchar(1)%c_keyword = 'BATHY_FILE'

    ! file name for input wave
    p_testparam%tchar(2)%c_keyword = 'IWAVE_FILE'

  END SUBROUTINE initial_setparam

!*******************************************************************************
! DESCRIPTION of [SUBROUTINE initialize_testcase]:
!> @brief initializes testcase (e.g. setup constants, auxiliary fields)
!>
!> @param[in]       p_ghand       grid handling data structure
!> @param[in,out]   p_param       control structure for global parameters
!> @param[in]       r_time        current model time
!>
!> @note This is just a dummy in this case.
!
  SUBROUTINE initialize_testcase(p_ghand, p_param, r_time)

    IMPLICIT NONE

    TYPE (grid_handle), INTENT(in)                        :: p_ghand
    TYPE (control_struct), INTENT(inout)                  :: p_param
    REAL (KIND = GRID_SR), INTENT(in)                     :: r_time

    CALL read_bathy(p_testparam%tchar(1)%p_value)
    CALL read_bheight(p_testparam%tchar(2)%p_value)

  END SUBROUTINE initialize_testcase

!*******************************************************************************
! DESCRIPTION of [SUBROUTINE initialvalues_iter]:
!> @brief initializes a grid with values
!>
!> @param[in]       p_ghand       grid handling data structure
!> @param[in,out]   p_param       control structure for global parameters
!> @param[in]       r_time        current model time
!>
!> @note This routine initializes the flow field.
!
  SUBROUTINE initialvalues_iter(p_ghand, p_param, r_time)

    IMPLICIT NONE

    TYPE (grid_handle), INTENT(in)                        :: p_ghand
    TYPE (control_struct), INTENT(inout)                  :: p_param
    REAL (KIND = GRID_SR), INTENT(in)                     :: r_time

!--- local declarations
    INTEGER (KIND = GRID_SI)                              :: i_faceunknowns, i_numelmt, &
      i_num, i_alct
    INTEGER (KIND = GRID_SI), DIMENSION(:,:), ALLOCATABLE :: i_elmtdofs
    REAL (KIND = GRID_SR), DIMENSION(:,:), ALLOCATABLE    :: r_coodof, r_Q, r_S

!--- initialize some constants
    i_faceunknowns = GRID_femtypes%p_type(FEM_DG)%sig%i_unknowns
    i_numelmt      = p_ghand%i_enumfine
    i_num          = i_numelmt*i_faceunknowns

!--- allocate workspace
    ALLOCATE(i_elmtdofs(i_faceunknowns, i_numelmt), &
             r_coodof(GRID_dimension, i_num), &
             r_Q(i_nprogvars,i_num), r_S(i_nsrcterms, i_num), STAT=i_alct)
    IF(i_alct /= 0) CALL grid_error(c_error='[initialvalues_iter]: could not allocate grid field')

!--- get grid information
    CALL grid_getinfo(p_ghand, i_femtype=FEM_DG, l_relative=.TRUE., &
                      l_finelevel=.TRUE., i_elementdofs=i_elmtdofs, r_dofcoordinates=r_coodof)

!--- initialize fields
    r_Q = 0.0_GRID_SR
    r_S = 0.0_GRID_SR
    CALL source_update(p_param, r_time, i_numelmt, i_elmtdofs, r_coodof, .TRUE., &
                       r_Q, r_S)
    CALL field_init(p_param, i_numelmt, i_elmtdofs, r_coodof, r_Q, r_S)

!--- update grid information
    CALL grid_putinfo(p_ghand, l_finelevel=.TRUE., i_arraypoint=i_valQ, r_dofvalues=r_Q)
    CALL grid_putinfo(p_ghand, l_finelevel=.TRUE., i_arraypoint=i_valS, r_dofvalues=r_S)

!--- deallocate workspace
    DEALLOCATE(i_elmtdofs, r_coodof, r_Q, r_S)

  END SUBROUTINE initialvalues_iter

!*******************************************************************************
! DESCRIPTION of [SUBROUTINE source_update]:
!> @brief updates the source terms
!>
!> @param[in]       p_param         control structure for global parameters
!> @param[in]       r_time          model time at which source term is evaluated
!> @param[in]       i_numelmt       number of elements
!> @param[in]       i_elmtdofs      global element-dof relation matrix
!> @param[in]       r_coodof        coordinates of DOFs
!> @param[in]       l_gridchanged
!> @param[in,out]   r_Q             discrete solution vector
!> @param[in,out]   r_S             discrete fields used in source terms
!
  SUBROUTINE source_update(p_param, r_time, i_numelmt, i_elmtdofs, r_coodof, &
                           l_gridchanged, r_Q, r_S)

    IMPLICIT NONE

    TYPE (control_struct),                    INTENT(in)    :: p_param
    REAL (KIND = GRID_SR),                    INTENT(in)    :: r_time
    INTEGER (KIND = GRID_SI),                 INTENT(in)    :: i_numelmt
    INTEGER (KIND = GRID_SI), DIMENSION(:,:), INTENT(in)    :: i_elmtdofs
    REAL (KIND = GRID_SR), DIMENSION(:,:),    INTENT(in)    :: r_coodof
    LOGICAL,                                  INTENT(in)    :: l_gridchanged
    REAL (KIND = GRID_SR), DIMENSION(:,:),    INTENT(inout) :: r_Q, r_S

!--- local declarations
    INTEGER (KIND = GRID_SI)                                :: i_faceunknowns, &
      i_elmt, i_dof, i_icnt, i_jcnt, i_globaldof
    REAL (KIND = GRID_SR), DIMENSION(2)                     :: r_alpha
    REAL (KIND = GRID_SR)                                   :: r_aux

    IF (l_gridchanged) THEN
!--- initialize some constants
      i_faceunknowns = GRID_femtypes%p_type(FEM_DG)%sig%i_unknowns

!--- loop over the elements
      DO i_elmt=1, i_numelmt
        DO i_dof=1, i_faceunknowns
          i_globaldof = i_elmtdofs(i_dof,i_elmt)
          r_aux = r_Q(1,i_globaldof) + r_S(1,i_globaldof)
          ! search (i,j) indices of source_update data grid cell for linear interpolation
          DO i_icnt=2,393
            IF (r_coodof(1,i_globaldof) <= r_xcoo(i_icnt)) EXIT
          END DO
          r_alpha(1) = (r_coodof(1,i_globaldof) - r_xcoo(i_icnt-1))/(r_xcoo(i_icnt) - r_xcoo(i_icnt-1))
          DO i_jcnt=2,244
            IF (r_coodof(2,i_globaldof) <= r_ycoo(i_jcnt)) EXIT
          END DO
          r_alpha(2) = (r_coodof(2,i_globaldof) - r_ycoo(i_jcnt-1))/(r_ycoo(i_jcnt) - r_ycoo(i_jcnt-1))

!--- bilinear interpolation of the source_update data
          r_S(1,i_globaldof) = &
            -((1.0-r_alpha(1))*(1.0-r_alpha(2))*r_bathygrid(i_icnt-1,i_jcnt-1) + &
              (    r_alpha(1))*(1.0-r_alpha(2))*r_bathygrid(i_icnt  ,i_jcnt-1) + &
              (1.0-r_alpha(1))*(    r_alpha(2))*r_bathygrid(i_icnt-1,i_jcnt  ) + &
              (    r_alpha(1))*(    r_alpha(2))*r_bathygrid(i_icnt  ,i_jcnt  ) - p_equationparam%r_depth) * GRID_GRAV
          r_Q(1,i_globaldof) = MAX(0.0_GRID_SR, r_aux - r_S(1,i_globaldof))
        END DO
      END DO
    END IF

  END SUBROUTINE source_update

!*******************************************************************************
! DESCRIPTION of [SUBROUTINE dg_elmt_solution]:
!> @brief computes the exact solution approximated in the DG space to compare
!>        with in diagnostics
!>
!> @param[in]       r_coodof        coordinates where solution should be computed
!> @param[in]       r_time          model time
!> @param[in]       i_numelmt       number of elements
!> @param[in]       i_elmtdofs      DOFs of element
!> @param[in]       r_coonod        coordinates of element nodes
!> @param[in]       i_elmtnodes     nodes of element
!> @param           r_sol           computed exact solution
!
  SUBROUTINE dg_elmt_solution(r_coodof, r_time, i_numelmt, i_elmtdofs, r_coonod, i_elmtnodes, r_sol)

    IMPLICIT NONE

    REAL (KIND = GRID_SR),    DIMENSION(:,:), INTENT(IN)        :: r_coodof, r_coonod
    REAL (KIND = GRID_SR),                    INTENT(IN)        :: r_time
    INTEGER (KIND = GRID_SI),                 INTENT(IN)        :: i_numelmt
    INTEGER (KIND = GRID_SI), DIMENSION(:,:), INTENT(IN)        :: i_elmtdofs, i_elmtnodes
    REAL (KIND = GRID_SR),    DIMENSION(:,:), INTENT(INOUT)     :: r_sol

    r_sol = 0.0_GRID_SR

  END SUBROUTINE dg_elmt_solution

!*******************************************************************************
  FUNCTION boundaryhgt(r_time) RESULT(r_height)

    IMPLICIT NONE

    REAL(KIND = GRID_SR)                                  :: r_time, r_height

    INTEGER(KIND = GRID_SI)                               :: i_cnt, i_ts
    REAL (KIND = GRID_SR)                                 :: r_alpha

!--- search time interval for linear interpolation
    DO i_cnt=2,451
      IF (r_time <= r_bheight(1,i_cnt)) THEN
        i_ts = i_cnt
        EXIT
      END IF
      i_ts = i_cnt
    END DO

    IF (r_time > r_bheight(1,451)) THEN
      r_alpha = 1.0
      i_ts   = 451
    ELSE
      r_alpha = (r_time - r_bheight(1,i_ts-1))/(r_bheight(1,i_ts) - r_bheight(1,i_ts-1))
    END IF

    r_height = ((1.0-r_alpha)*r_bheight(2,i_ts-1) + r_alpha*r_bheight(2,i_ts)) * GRID_GRAV

  END FUNCTION boundaryhgt

!*******************************************************************************
! Only internal subroutines follow
!*******************************************************************************
! DESCRIPTION of [SUBROUTINE field_init]:
!> @brief initializes the height and momentum fields
!>
!> @param[in,out]   p_param       control structure for global parameters
!> @param           i_elmtdofs
!> @param           r_coodof      coordinates of DOFs
!
  SUBROUTINE field_init(p_param, i_numelmt, i_elmtdofs, r_coodof, r_Q, r_S)

    IMPLICIT NONE

    TYPE (control_struct),                    INTENT(inout) :: p_param
    INTEGER (KIND = GRID_SI),                 INTENT(in)    :: i_numelmt
    INTEGER (KIND = GRID_SI), DIMENSION(:,:), INTENT(in)    :: i_elmtdofs
    REAL (KIND = GRID_SR), DIMENSION(:,:),    INTENT(in)    :: r_coodof
    REAL (KIND = GRID_SR), DIMENSION(:,:),    INTENT(inout) :: r_Q
    REAL (KIND = GRID_SR), DIMENSION(:,:),    INTENT(in)    :: r_S

!--- local declarations
    INTEGER (KIND = GRID_SI)                                :: i_count

!--- loop over the DOFs
    dof_loop: DO i_count= 1, i_numelmt*GRID_femtypes%p_type(FEM_DG)%sig%i_unknowns
      r_Q(1,i_count) = MAX(0.0_GRID_SR, p_equationparam%r_depth * GRID_GRAV - r_S(1,i_count))
      r_Q(2,i_count) = r_Q(1,i_count) * 0.0_GRID_SR
      r_Q(3,i_count) = r_Q(1,i_count) * 0.0_GRID_SR
    END DO dof_loop

  END SUBROUTINE field_init

!*******************************************************************************
! DESCRIPTION of [SUBROUTINE read_bathy]:
!> @brief reads bathymetry from file
!>
!> @param[in]       c_bathyfile   filename for bathymetry file
!
  SUBROUTINE read_bathy(c_bathyfile)

    IMPLICIT NONE

    CHARACTER (len=*), INTENT(in)                           :: c_bathyfile

!--- local declarations
    INTEGER                                                 :: i_iounit, i_iost, &
      i_icnt, i_jcnt
    CHARACTER (len=80)                                      :: c_filrow
    REAL (KIND = GRID_SR), DIMENSION(3)                     :: r_rawdata

!--- open and read the source_update values (393x244 values)
      OPEN(NEWUNIT=i_iounit, file=p_testparam%tchar(1)%p_value, status='old', iostat=i_iost)
      IF(i_iost /= 0) &
        CALL grid_error(c_error='[source_update]: could not open source_update data file')
      IF(GRID_parameters%iolog > 0) &
        WRITE(GRID_parameters%iolog,*) 'INFO [source_update]: opened file', &
        p_testparam%tchar(1)%p_value,' on unit:',i_iounit

!--- read first comment line
      read(i_iounit,'(A80)',iostat=i_iost) c_filrow

!--- read data
      DO i_icnt = 1,393
        DO i_jcnt = 1,244
          READ(i_iounit,*) r_rawdata(1), r_rawdata(2), r_rawdata(3)
          r_bathygrid(i_icnt, i_jcnt) = r_rawdata(3)
          IF (i_icnt == 1) r_ycoo(i_jcnt) = r_rawdata(2)
        END DO
        r_xcoo(i_icnt) = r_rawdata(1)
      END DO

      CLOSE(i_iounit)
      IF(GRID_parameters%iolog > 0) &
        WRITE(GRID_parameters%iolog,*) 'INFO [read_bathy]: Closed file on unit: ', i_iounit

  END SUBROUTINE read_bathy

!*******************************************************************************
  SUBROUTINE read_bheight(c_iwavefile)

    IMPLICIT NONE

    CHARACTER (len=*), INTENT(in)                           :: c_iwavefile

!--- local declarations
    INTEGER                                                 :: i_iounit, i_iost, i_count
    CHARACTER (len=80)                                      :: a_filrow
    REAL (KIND = GRID_SR), DIMENSION(2)                     :: r_rawdata

!--- initialize consts
    i_iounit = 10

!--- open and read the boundary values
    OPEN(i_iounit, file=c_iwavefile, status='old', iostat=i_iost)
    IF(i_iost /= 0) &
      CALL grid_error(c_error='[read_bheight]: could not open boundary data file')
    IF(GRID_parameters%iolog > 0) &
      WRITE(GRID_parameters%iolog,*) 'INFO [read_bheight]: opened file', &
      c_iwavefile,' on unit:',i_iounit

!--- read first comment line
    read(i_iounit,'(A80)',iostat=i_iost) a_filrow

!--- read data
    DO i_count = 1,451
      READ(i_iounit,*) r_rawdata(1), r_rawdata(2)
      r_bheight(1,i_count) = r_rawdata(1)
      r_bheight(2,i_count) = r_rawdata(2)+0.13535_GRID_SR
    END DO

    close(i_iounit)
    IF(GRID_parameters%iolog > 0) &
      write(GRID_parameters%iolog,*) 'INFO [read_bheight]: Closed file on unit: ', i_iounit

  END SUBROUTINE read_bheight

!*******************************************************************************
END MODULE DG_initial
