!*******************************************************************************
!
!> @file  DG_initial_standingwave.f90
!> @brief contains module DG_initial
!
!*******************************************************************************
! MODULE DESCRIPTION:
!> @brief Initial condition for a linear standing wave with constant source_update
!> and zero velocitites. The background depth is given by BATHY_PARAM_DEPTH
!> in the Parameters file. This test is adapted from Haiyang CUI 3.4.1 page 58.
!> 
!>

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
  INTEGER, PARAMETER                                  :: i_ntstintparam  = 0
  INTEGER, PARAMETER                                  :: i_ntstrealparam = 0
  INTEGER, PARAMETER                                  :: i_ntstcharparam = 0
  INTEGER, PARAMETER                                  :: i_ntstlogparam  = 0

  CONTAINS

!*******************************************************************************
! DESCRIPTION of [SUBROUTINE initial_setparam]:
!> @brief initializes test case specific parameter list
!
  SUBROUTINE initial_setparam

    IMPLICIT NONE

!     ! parameter description
!     p_testparam%tint(1)%c_keyword  = 'PARAM_INT'
!     p_testparam%tint(1)%i_size     = 1
!
     ! parameter description
     p_testparam%treal(1)%c_keyword = 'AMPLITUDE'
     p_testparam%treal(1)%i_size    = 1

     ! parameter description
     p_testparam%treal(2)%c_keyword = 'WAVELENGTH'
     p_testparam%treal(2)%i_size    = 1

     ! parameter description
     p_testparam%treal(3)%c_keyword = 'BATHY_PARAM_DEPTH'
     p_testparam%treal(3)%i_size    = 1
!
!     ! parameter description
!     p_testparam%tchar(1)%c_keyword = 'PARAM_CHAR'
!
!     ! parameter description
!     p_testparam%tlog(1)%c_keyword  = 'PARAM_LOG'

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
      i_elmt, i_dof

    IF (l_gridchanged) THEN
!--- initialize some constants
      i_faceunknowns = GRID_femtypes%p_type(FEM_DG)%sig%i_unknowns

!--- loop over the elements
      DO i_elmt=1, i_numelmt
        DO i_dof=1, i_faceunknowns
          r_S(1,i_elmtdofs(i_dof,i_elmt)) = 0.0_GRID_SR * GRID_GRAV
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

  FUNCTION exact_solution(r_coo, r_time) RESULT (r_sol)

    IMPLICIT NONE

    REAL (KIND = GRID_SR), DIMENSION(GRID_dimension), INTENT(in)  :: r_coo
    REAL (KIND = GRID_SR), INTENT(in)                             :: r_time
    REAL (KIND = GRID_SR), DIMENSION(3)                           :: r_sol

!--- local declarations
    REAL (KIND = GRID_SR)                                         :: r_x, r_w, r_eta, r_u, r_v

!--- initialize constants
!--- we neglect y direction. eta(x,y_0,t)=eta(x,y_1,t) for all y_0 and y_1
    r_x  = r_coo(1)
    r_w = 2.0_GRID_SR * GRID_PI/p_testparam%treal(2)%p_value(1) * &
      SQRT(0.5_GRID_SR * GRID_GRAV * p_testparam%treal(2)%p_value(1)/GRID_PI * &
      TANH(-p_testparam%treal(3)%p_value(1) * p_testparam%treal(2)%p_value(1)))       !Angular Velocity
    r_eta = p_testparam%treal(1)%p_value(1) * &
      SIN(2.0_GRID_SR*GRID_PI*(r_x-p_testparam%treal(2)%p_value(1)/4.0_GRID_SR)/p_testparam%treal(2)%p_value(1)) * &
      COS(r_w*r_time)
    r_u = 0.0_GRID_SR
    r_v = 0.0_GRID_SR


!---old code
    !r_c0 = SQRT(2.0_GRID_SR * p_equationparam%r_depth * GRID_GRAV)
    !r_x  = r_coo(1)
    !r_y  = r_coo(2)

    !r_eta = p_equationparam%r_depth + COS(GRID_PI*r_x) * COS(GRID_PI*r_y) * COS(r_c0*GRID_PI*r_time)
    !r_u   = 1.0_GRID_SR/r_c0 * SIN(GRID_PI*r_x) * COS(GRID_PI*r_y) * SIN(r_c0*GRID_PI*r_time)
    !r_v   = 1.0_GRID_SR/r_c0 * COS(GRID_PI*r_x) * SIN(GRID_PI*r_y) * SIN(r_c0*GRID_PI*r_time)

    r_sol(1) = r_eta
    r_sol(2) = r_eta*r_u
    r_sol(3) = r_eta*r_v

  END FUNCTION exact_solution

!*******************************************************************************
! DESCRIPTION of [SUBROUTINE field_init]:
!> @brief initializes the height and momentum fields
!>
!> @param[in,out]   p_param       control structure for global parameters
!> @param[in]       i_elmtdofs
!> @param           r_coodof      coordinates of DOFs
!
  SUBROUTINE field_init(p_param, i_numelmt, i_elmtdofs, i_elmtnodes, r_coodof, r_coonod, r_Q)

    IMPLICIT NONE

    TYPE (control_struct),                    INTENT(inout) :: p_param
    INTEGER (KIND = GRID_SI),                 INTENT(in)    :: i_numelmt
    INTEGER (KIND = GRID_SI), DIMENSION(:,:), INTENT(in)    :: i_elmtdofs, i_elmtnodes
    REAL (KIND = GRID_SR), DIMENSION(:,:),    INTENT(in)    :: r_coodof, r_coonod
    REAL (KIND = GRID_SR), DIMENSION(:,:),    INTENT(inout) :: r_Q

!--- local declarations
    INTEGER (KIND = GRID_SI)                                :: i_faceunknowns, i_alct, i_dof
    REAL (KIND = GRID_SR), DIMENSION(:,:), ALLOCATABLE      :: r_sol
    INTEGER (KIND = GRID_SI)                                :: i_elmt
    !REAL (KIND = GRID_SR), DIMENSION(1,:)                   :: r_tmp

!--- initialize some constants
    i_faceunknowns = GRID_femtypes%p_type(FEM_DG)%sig%i_unknowns

!--- allocate workspace
    ALLOCATE(r_sol(i_nprogvars, i_numelmt*i_faceunknowns), stat=i_alct)
    IF(i_alct /= 0) THEN
      CALL grid_error(c_error='[field_init]: could not allocate workspace')
    END IF

    CALL dg_elmt_solution(r_coodof, 0.0_GRID_SR, i_numelmt, i_elmtdofs, r_coonod, i_elmtnodes, r_sol)
    !r_Q(1,:) = (r_sol(1,:) - p_testparam%treal(3)%p_value(1)) * GRID_GRAV
    !r_Q(2,:) = 0.0_GRID_SR
    !r_Q(3,:) = 0.0_GRID_SR

    !DO i_elmt=1, i_numelmt
    !  r_Q(1,i_elmtdofs(:,i_elmt)) = (r_sol(1,i_elmtdofs(:,i_elmt)) - p_testparam%treal(3)%p_value(1)) * GRID_GRAV
    !  r_Q(2,i_elmtdofs(:,i_elmt)) = 0.0_GRID_SR
    !  r_Q(3,i_elmtdofs(:,i_elmt)) = 0.0_GRID_SR
    !END DO

!r_sol(:,i_elmtdofs(i_dof,i_elmt))

!Neuer Versuch
    !r_tmp = r_sol(1,:)
    !r_tmp = r_tmp * 0.0_GRID_SR
    !r_tmp = p_testparam%treal(3)%p_value(1)

!Versuch2
    !DO i_elmt=1, i_numelmt
     ! DO i_dof=1, i_faceunknowns
      !  r_Q(1,i_elmtdofs(i_dof,i_elmt)) = (r_sol(1,i_elmtdofs(i_dof,i_elmt)) - p_testparam%treal(3)%p_value(1)) * GRID_GRAV
       ! r_Q(2,i_elmtdofs(i_dof,i_elmt)) = 0.0_GRID_SR
        !r_Q(3,i_elmtdofs(i_dof,i_elmt)) = 0.0_GRID_SR
      !END DO
    !END DO
    r_Q = r_sol * GRID_GRAV



    

!--- deallocate workspace
    DEALLOCATE(r_sol)

  END SUBROUTINE field_init

!*******************************************************************************
END MODULE DG_initial
