!*******************************************************************************
!
!> @file  DG_initial_WB.f90
!> @brief contains module DG_initial
!
!*******************************************************************************
! MODULE DESCRIPTION:
!> @brief initializes fields for a simple wellbalancing testcase with a
!>        bathymetry bump and an initial flat sea surface and zero initial
!>        momentum. The general shape of the bump is given by BATHY_TYPE:
!>        1: parabolically shaped bump with surrounding flat bathymetry
!>        2: exponentially shaped bump
!>        Further parameters are given by ALPHA, BUMP_CENTRE, BUMP_HEIGHT.
!>        A small localized perturbation can be set with PERT_CENTRE and
!>        PERT_HEIGHT (note: in case PERT_HEIGHT > 0 only the initial condition
!>        is given by the exact solution).
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
          source_update, dg_elmt_solution

  TYPE (test_param)                                   :: p_testparam
  INTEGER, PARAMETER                                  :: i_ntstintparam  = 1
  INTEGER, PARAMETER                                  :: i_ntstrealparam = 5
  INTEGER, PARAMETER                                  :: i_ntstcharparam = 0
  INTEGER, PARAMETER                                  :: i_ntstlogparam  = 0

  CONTAINS

!*******************************************************************************
! DESCRIPTION of [SUBROUTINE initial_setparam]:
!> @brief initializes test case specific parameter list
!
  SUBROUTINE initial_setparam

    IMPLICIT NONE

    ! type of bathymetry bump
    p_testparam%tint(1)%c_keyword = 'BATHY_TYPE'
    p_testparam%tint(1)%i_size    = 1

    ! shape parameter for bathymetry bump
    p_testparam%treal(1)%c_keyword = 'ALPHA'
    p_testparam%treal(1)%i_size    = 1

    ! midpoint of bathymetry bump
    p_testparam%treal(2)%c_keyword = 'BUMP_CENTRE'
    p_testparam%treal(2)%i_size    = 2

    ! height of bathymetry bump
    p_testparam%treal(3)%c_keyword = 'BUMP_HEIGHT'
    p_testparam%treal(3)%i_size    = 1

    ! midpoint of perturbation in fluid surface
    p_testparam%treal(4)%c_keyword = 'PERT_CENTRE'
    p_testparam%treal(4)%i_size    = 2

    ! height of perturbation in fluid surface
    p_testparam%treal(5)%c_keyword = 'PERT_HEIGHT'
    p_testparam%treal(5)%i_size    = 1

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
    INTEGER (KIND = GRID_SI), DIMENSION(:,:), ALLOCATABLE :: i_elmtdofs, i_elmtnodes
    REAL (KIND = GRID_SR), DIMENSION(:,:), ALLOCATABLE    :: r_coodof, r_coonod, r_Q, r_S

!--- initialize some constants
    i_faceunknowns = GRID_femtypes%p_type(FEM_DG)%sig%i_unknowns
    i_numelmt      = p_ghand%i_enumfine
    i_num          = i_numelmt*i_faceunknowns

!--- allocate workspace
    ALLOCATE(i_elmtdofs(i_faceunknowns, i_numelmt), &
             i_elmtnodes(GRID_elementnodes, i_numelmt), &
             r_coodof(GRID_dimension, i_num), &
             r_coonod(GRID_dimension, i_numelmt*GRID_elementnodes), &
             r_Q(i_nprogvars,i_num), r_S(i_nsrcterms, i_num), STAT=i_alct)
    IF(i_alct /= 0) CALL grid_error(c_error='[initialvalues_iter]: could not allocate grid field')

!--- get grid information
    CALL grid_getinfo(p_ghand, i_femtype=FEM_DG, l_relative=.TRUE., &
                      l_finelevel=.TRUE., i_elementdofs=i_elmtdofs, &
                      i_elementnodes=i_elmtnodes, r_dofcoordinates=r_coodof, &
                      r_nodecoordinates=r_coonod)

!--- initialize fields
    r_Q = 0.0_GRID_SR
    CALL source_update(p_param, r_time, i_numelmt, i_elmtdofs, r_coodof, .TRUE., &
                       r_Q, r_S)
    CALL field_init(p_param, i_numelmt, i_elmtdofs, i_elmtnodes, r_coodof, r_coonod, r_Q)

!--- update grid information
    CALL grid_putinfo(p_ghand, l_finelevel=.TRUE., i_arraypoint=i_valQ, r_dofvalues=r_Q)
    CALL grid_putinfo(p_ghand, l_finelevel=.TRUE., i_arraypoint=i_valS, r_dofvalues=r_S)

!--- deallocate workspace
    DEALLOCATE(i_elmtdofs, i_elmtnodes, r_coodof, r_coonod, r_Q, r_S)

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
      i_elmt, i_dof, i_case

!--- initialize some constants
    i_case = p_testparam%tint(1)%p_value(1)

    IF (l_gridchanged) THEN
!--- initialize some constants
      i_faceunknowns = GRID_femtypes%p_type(FEM_DG)%sig%i_unknowns

!--- loop over the elements
      DO i_elmt=1, i_numelmt
        DO i_dof=1, i_faceunknowns
          r_S(1,i_elmtdofs(i_dof,i_elmt)) = &
            GRID_GRAV * bathy(r_coodof(:,i_elmtdofs(i_dof,i_elmt)), r_time, i_case)
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

!--- local declarations
    INTEGER (KIND = GRID_SI)                                    :: i_faceunknowns, &
                                                                   i_elmt, i_dof

!--- initialize some constants
    i_faceunknowns = GRID_femtypes%p_type(FEM_DG)%sig%i_unknowns

!--- loop over the elements and DOFs
    DO i_elmt=1, i_numelmt
      DO i_dof=1, i_faceunknowns
        r_sol(:,i_elmtdofs(i_dof,i_elmt)) = exact_solution(r_coodof(:,i_elmtdofs(i_dof,i_elmt)), r_time)
      END DO
    END DO

  END SUBROUTINE dg_elmt_solution

!*******************************************************************************
! Only internal subroutines follow
!*******************************************************************************
! DESCRIPTION of [FUNCTION exact_solution]:
!> @brief computes the exact solution to the problem at given time and coordinates
!>
!> @param[in]       r_coo         coordinates where solution should be computed
!> @param[in]       r_time        model time
!> @return                        exact solution
!
  FUNCTION exact_solution(r_coo, r_time) RESULT(r_sol)

    IMPLICIT NONE

    REAL (KIND = GRID_SR), DIMENSION(GRID_dimension), INTENT(IN)  :: r_coo
    REAL (KIND = GRID_SR),                            INTENT(IN)  :: r_time
    REAL (KIND = GRID_SR), DIMENSION(3)                           :: r_sol

!--- local declarations
    INTEGER (KIND = GRID_SI)                                      :: i_case
    REAL (KIND = GRID_SR)                                         :: r_pert, r_phgt
    REAL (KIND = GRID_SR), DIMENSION(2)                           :: r_pmid

!--- initialize some constants
    i_case = p_testparam%tint(1)%p_value(1)
    r_pmid  = p_testparam%treal(4)%p_value
    r_phgt  = p_testparam%treal(5)%p_value(1)

    r_pert = MAX(0.0_GRID_SR, r_phgt - 30.0_GRID_SR*((r_coo(1)-r_pmid(1))**2+(r_coo(2)-r_pmid(2))**2))

    IF ((r_time == 0.0_GRID_SR) .OR. (r_phgt == 0.0_GRID_SR)) THEN
      r_sol(1) = MAX(0.0_GRID_SR, p_equationparam%r_depth - bathy(r_coo, r_time, i_case)) + r_pert
    ELSE
      r_sol(1) = 0.0_GRID_SR
    END IF

    r_sol(2) = 0.0_GRID_SR
    r_sol(3) = 0.0_GRID_SR

  END FUNCTION exact_solution

!*******************************************************************************
! DESCRIPTION of [SUBROUTINE field_init]:
!> @brief initializes the height and momentum fields
!>
!> @param[in,out]   p_param       control structure for global parameters
!> @param[in]       i_elmtdofs
!> @param[in]       r_coodof      coordinates of DOFs
!
  SUBROUTINE field_init(p_param, i_numelmt, i_elmtdofs, i_elmtnodes, r_coodof, r_coonod, r_Q)

    IMPLICIT NONE

    TYPE (control_struct),                    INTENT(inout) :: p_param
    INTEGER (KIND = GRID_SI),                 INTENT(in)    :: i_numelmt
    INTEGER (KIND = GRID_SI), DIMENSION(:,:), INTENT(in)    :: i_elmtdofs, i_elmtnodes
    REAL (KIND = GRID_SR), DIMENSION(:,:),    INTENT(in)    :: r_coodof, r_coonod
    REAL (KIND = GRID_SR), DIMENSION(:,:),    INTENT(inout) :: r_Q

!--- local declarations
    INTEGER (KIND = GRID_SI)                                :: i_faceunknowns, i_alct
    REAL (KIND = GRID_SR), DIMENSION(:,:), ALLOCATABLE      :: r_sol

!--- initialize some constants
    i_faceunknowns = GRID_femtypes%p_type(FEM_DG)%sig%i_unknowns

!--- allocate workspace
    ALLOCATE(r_sol(i_nprogvars, i_numelmt*i_faceunknowns), stat=i_alct)
    IF(i_alct /= 0) THEN
      CALL grid_error(c_error='[field_init]: could not allocate workspace')
    END IF

    CALL dg_elmt_solution(r_coodof, 0.0_GRID_SR, i_numelmt, i_elmtdofs, r_coonod, &
                          i_elmtnodes, r_sol)
    r_Q = r_sol * GRID_GRAV

!--- deallocate workspace
    DEALLOCATE(r_sol)

  END SUBROUTINE field_init

!*******************************************************************************
! DESCRIPTION of [FUNCTION bathy]:
!> @brief computes the exact solution to the problem at given time and coordinates
!>
!> @param[in]       r_coo         coordinates where solution should be computed
!> @param[in]       r_time        model time
!> @param[in]       i_case        shape of bathymetry (1: parabolic bump with surrounding flat bathymetry)
!> @return                        bathymetry at r_coo
!
  FUNCTION bathy(r_coo, r_time, i_case) RESULT(r_bathy)

    IMPLICIT NONE

    REAL (KIND = GRID_SR),    DIMENSION(GRID_dimension), INTENT(IN)  :: r_coo
    REAL (KIND = GRID_SR),                               INTENT(IN)  :: r_time
    INTEGER (KIND = GRID_SI),                            INTENT(IN)  :: i_case
    REAL (KIND = GRID_SR)                                            :: r_bathy

!--- local declarations
    REAL (KIND = GRID_SR)                                            :: r_alpha, r_bhgt
    REAL (KIND = GRID_SR), DIMENSION(2)                              :: r_mid

!--- initialize some constants
    r_alpha = p_testparam%treal(1)%p_value(1)
    r_mid   = p_testparam%treal(2)%p_value
    r_bhgt  = p_testparam%treal(3)%p_value(1)

    SELECT CASE (i_case)
    CASE(1) ! parabolically shaped bump with surrounding flat bathymetry
      r_bathy = MAX(0.0_GRID_SR, r_bhgt - r_alpha*((r_coo(1)-r_mid(1))**2+(r_coo(2)-r_mid(2))**2))

    CASE(2) ! exponentially shaped bump
      r_bathy = r_bhgt * EXP(-((r_coo(1)-r_mid(1))**2+(r_coo(2)-r_mid(2))**2) / r_alpha)

    CASE DEFAULT
      CALL grid_error(c_error='[bathy]: bathymetry type not supported!')
    END SELECT

  END FUNCTION bathy

!*******************************************************************************
END MODULE DG_initial
