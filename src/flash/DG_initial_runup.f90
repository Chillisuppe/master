!*******************************************************************************
!
!> @file  DG_initial_runup.f90
!> @brief contains module DG_initial
!
!*******************************************************************************
! MODULE DESCRIPTION:
!> @brief initializes fields with a quasi 1D tsunami runup onto a plane beach
!
MODULE DG_initial

  USE GRID_api
  USE FLASH_parameters
  USE DG_equation, ONLY : i_nprogvars, i_nsrcterms, i_valQ, i_valS

  PRIVATE
  PUBLIC  i_ntstintparam, i_ntstrealparam, i_ntstcharparam, i_ntstlogparam, &
          initial_setparam, p_testparam, initialize_testcase, initialvalues_iter, &
          source_update, dg_elmt_solution

  TYPE (test_param)                                   :: p_testparam
  INTEGER, PARAMETER                                  :: i_ntstintparam  = 0
  INTEGER, PARAMETER                                  :: i_ntstrealparam = 0
  INTEGER, PARAMETER                                  :: i_ntstcharparam = 1
  INTEGER, PARAMETER                                  :: i_ntstlogparam  = 0

  CONTAINS

!*******************************************************************************
! DESCRIPTION of [SUBROUTINE initial_setparam]:
!> @brief initializes test case specific parameter list
!
  SUBROUTINE initial_setparam

    IMPLICIT NONE

    p_testparam%tchar(1)%c_keyword = 'SSHINIT_FILE'

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
             r_Q(i_nprogvars, i_num), r_S(i_nsrcterms, i_num), STAT=i_alct)
    IF(i_alct /= 0) CALL grid_error(c_error='[initialvalues_iter]: could not allocate grid field')

!--- get grid information
    CALL grid_getinfo(p_ghand, i_femtype=FEM_DG, l_relative=.TRUE., &
                      l_finelevel=.TRUE., i_elementdofs=i_elmtdofs, r_dofcoordinates=r_coodof)

!--- initialize fields
    r_Q = 0.0_GRID_SR
    CALL source_update(p_param, r_time, i_numelmt, i_elmtdofs, r_coodof, .TRUE., &
                       r_Q, r_S)
    CALL field_init(p_param, i_numelmt, i_elmtdofs, r_coodof, r_Q)

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
      i_elmt, i_dof, i_globaldof
    REAL (KIND = GRID_SR)                                   :: r_aux

    IF (l_gridchanged) THEN
!--- initialize some constants
      i_faceunknowns = GRID_femtypes%p_type(FEM_DG)%sig%i_unknowns

!--- loop over the elements
      DO i_elmt=1, i_numelmt
        DO i_dof=1, i_faceunknowns
          i_globaldof = i_elmtdofs(i_dof,i_elmt)
          r_aux = r_Q(1,i_globaldof) + r_S(1,i_globaldof)
          r_S(1,i_globaldof) = &
            (5000.0_GRID_SR - 0.1_GRID_SR * r_coodof(1,i_globaldof)) * GRID_GRAV ! depth(x) = 1/10 * x
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
  SUBROUTINE dg_elmt_solution(r_coodof, r_time, i_numelmt, &
                              i_elmtdofs, r_coonod, i_elmtnodes, r_sol)

    IMPLICIT NONE

    REAL (KIND = GRID_SR),    DIMENSION(:,:), INTENT(IN)        :: r_coodof, r_coonod
    REAL (KIND = GRID_SR),                    INTENT(IN)        :: r_time
    INTEGER (KIND = GRID_SI),                 INTENT(IN)        :: i_numelmt
    INTEGER (KIND = GRID_SI), DIMENSION(:,:), INTENT(IN)        :: i_elmtdofs, i_elmtnodes
    REAL (KIND = GRID_SR),    DIMENSION(:,:), INTENT(INOUT)     :: r_sol

    r_sol(1,:) = 0.0_GRID_SR
    r_sol(2,:) = 0.0_GRID_SR
    r_sol(3,:) = 0.0_GRID_SR

  END SUBROUTINE dg_elmt_solution

!*******************************************************************************
! Only internal subroutines follow
!*******************************************************************************
! DESCRIPTION of [SUBROUTINE field_init]:
!> @brief initializes the height and momentum fields
!>
!> @param[in,out]   p_param       control structure for global parameters
!> @param           i_numelmt
!> @param[in]       r_coodof      coordinates of DOFs
!
  SUBROUTINE field_init(p_param, i_numelmt, i_elmtdofs, r_coodof, r_Q)

    IMPLICIT NONE

    TYPE (control_struct),                    INTENT(inout) :: p_param
    INTEGER (KIND = GRID_SI),                 INTENT(in)    :: i_numelmt
    INTEGER (KIND = GRID_SI), DIMENSION(:,:), INTENT(in)    :: i_elmtdofs
    REAL (KIND = GRID_SR), DIMENSION(:,:),    INTENT(in)    :: r_coodof
    REAL (KIND = GRID_SR), DIMENSION(:,:),    INTENT(inout) :: r_Q

!--- local declarations
    INTEGER                                               :: i_count, i_iounit, &
      i_filelen, i_alct, i_cnt, j_cnt, i_iost, i_ioend
    REAL (KIND = GRID_SR), DIMENSION(:,:), ALLOCATABLE    :: r_rawdata
    REAL (KIND = GRID_SR)                                 :: r_value, r_alpha
    LOGICAL                                               :: l_between
    CHARACTER (len=80)                                    :: a_filrow

!--- open and read the inital wave height values (2022 values)
    OPEN(NEWUNIT=i_iounit, file=p_testparam%tchar(1)%p_value, status='old', iostat=i_iost)
    IF(i_iost /= 0) &
      CALL grid_error(c_error='[field_init]: could not open initial data file')
    IF(GRID_parameters%iolog > 0) &
      WRITE(GRID_parameters%iolog,*) 'INFO [field_init]: opened file', &
      p_testparam%tchar(1)%p_value,' on unit:', i_iounit

!--- read line by line
    read_loop: DO
      read(i_iounit, '(a80)', iostat=i_ioend) a_filrow

!--- if file ended
      file_end: IF(i_ioend /= 0) THEN
        close(i_iounit)
        IF(GRID_parameters%iolog > 0) &
          write(GRID_parameters%iolog,*) 'INFO: Closed file on unit: ', i_iounit
        EXIT
      ELSE file_end
        comment_line: IF(a_filrow(1:1) == '#' .or. a_filrow(1:1) == '!') THEN
          CYCLE read_loop
        ELSE IF(a_filrow(1:8) == 'NUMDATAP') THEN comment_line
          read(i_iounit,*) i_filelen
          ALLOCATE(r_rawdata(2,i_filelen), stat=i_alct)
          IF(i_alct /= 0)  &
            CALL grid_error(c_error='[field_init]: could not allocate raw data array')
        ELSE IF(a_filrow(1:8) == 'DATAPOIN') THEN comment_line
          IF(.NOT. ALLOCATED(r_rawdata)) &
            CALL grid_error(c_error='[field_init]: raw data array not allocated')
          DO j_cnt=1,i_filelen
            READ(i_iounit,*) r_rawdata(1,j_cnt), r_rawdata(2,j_cnt)
          END DO
        END IF comment_line
      END IF file_end
    END DO read_loop

    IF(GRID_parameters%iolog > 0) &
      WRITE(GRID_parameters%iolog,*) 'INFO [field_init]: closed unit:',i_iounit

!--- determine the value of inital surface elevation for each dof
    dof_loop: DO i_count=1, i_numelmt*GRID_femtypes%p_type(FEM_DG)%sig%i_unknowns

!--- initialize aux variables
       i_cnt = 1
       l_between = .FALSE.

!--- find i such that the given x-coordinate is between two
!    given data points, i.e.
!    X(1,i) <= x-coordinate of the node <= X(1,i+1)
!    first make sure, X(1,1) <= x-coordinate!
      IF(r_rawdata(1,1) > r_coodof(1, i_count)) THEN !>
        l_between = .TRUE.
        i_cnt     = 2 !1
      ELSE IF(r_rawdata(1,i_filelen) < r_coodof(1,i_count)) THEN !<
        l_between = .TRUE.
        i_cnt     = i_filelen !- 1
      END IF

      search_loop: DO WHILE (.NOT. l_between)
        l_between = (r_rawdata(1,i_cnt) <= r_coodof(1,i_count)) .AND. &
                    (r_coodof(1,i_count) <= r_rawdata(1,i_cnt+1))
        i_cnt     = i_cnt + 1
      END DO search_loop

      i_cnt = i_cnt - 1

!--- interpolate the value at x linearly between X(1,i) and X(1,i+1)
      r_alpha = (r_coodof(1,i_count) - r_rawdata(1,i_cnt))/(r_rawdata(1,i_cnt+1) - r_rawdata(1,i_cnt))
      r_value = (r_rawdata(2,i_cnt+1) - r_rawdata(2, i_cnt))*r_alpha + r_rawdata(2,i_cnt)

      IF (r_coodof(1,i_count) < 0.0) THEN
        r_Q(1,i_count) = 0.0_GRID_SR * GRID_GRAV
      ELSE
        r_Q(1,i_count) = (r_value + MAX(0.0_GRID_SR, 0.1_GRID_SR * r_coodof(1,i_count))) * GRID_GRAV
      END IF
      r_Q(2,i_count) = 0.0_GRID_SR
      r_Q(3,i_count) = 0.0_GRID_SR
    END DO dof_loop

!--- deallocate workspace
    DEALLOCATE(r_rawdata)

  END SUBROUTINE field_init

!*******************************************************************************
END MODULE DG_initial
