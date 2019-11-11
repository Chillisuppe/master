!*******************************************************************************
!
!> @file  DG_initial_warmairbubble.f90
!> @brief contains module DG_initial
!
!*******************************************************************************
! MODULE DESCRIPTION:
!> @brief initialize fields with warm air bubble, see Giraldo and Restelli (2007)
!
MODULE DG_initial

  USE GRID_api
  USE FLASH_parameters
  USE DG_equation, ONLY : i_nprogvars, i_valQ, PrimitiveToProgvars
  USE IO_equation, ONLY : p_equationparam

  PRIVATE
  PUBLIC  i_ntstintparam, i_ntstrealparam, i_ntstcharparam, i_ntstlogparam, &
          initial_setparam, p_testparam, initialize_testcase, initialvalues_iter, &
          source_update, dg_elmt_solution

  TYPE (test_param)                                   :: p_testparam
  INTEGER, PARAMETER                                  :: i_ntstintparam  = 0
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

    ! background potential temperature
    p_testparam%treal(1)%c_keyword = 'PT_BACKGROUND'
    p_testparam%treal(1)%i_size    = 1

    ! potential temperature deviation
    p_testparam%treal(2)%c_keyword = 'PT_DEVIATION'
    p_testparam%treal(2)%i_size    = 1

    ! initial velocity (u0, v0)
    p_testparam%treal(3)%c_keyword = 'INITIAL_VELO'
    p_testparam%treal(3)%i_size    = 2

    ! bubble center position
    p_testparam%treal(4)%c_keyword = 'BUBBLE_CENTER'
    p_testparam%treal(4)%i_size    = 2

    ! bubble radius
    p_testparam%treal(5)%c_keyword = 'BUBBLE_RADIUS'
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
    REAL (KIND = GRID_SR), DIMENSION(:,:), ALLOCATABLE    :: r_coodof, r_coonod, r_Q

!--- initialize some constants
    i_faceunknowns = GRID_femtypes%p_type(FEM_DG)%sig%i_unknowns
    i_numelmt      = p_ghand%i_enumfine
    i_num          = i_numelmt*i_faceunknowns

!--- allocate workspace
    ALLOCATE(i_elmtdofs(i_faceunknowns, i_numelmt), &
             i_elmtnodes(GRID_elementnodes, i_numelmt), &
             r_coodof(GRID_dimension, i_num), &
             r_coonod(GRID_dimension, i_numelmt*GRID_elementnodes), &
             r_Q(i_nprogvars,i_num), STAT=i_alct)
    IF(i_alct /= 0) CALL grid_error(c_error='[initialvalues_iter]: could not allocate grid field')

!--- get grid information
    CALL grid_getinfo(p_ghand, i_femtype=FEM_DG, l_relative=.TRUE., &
                      l_finelevel=.TRUE., i_elementdofs=i_elmtdofs, &
                      i_elementnodes=i_elmtnodes, r_dofcoordinates=r_coodof, &
                      r_nodecoordinates=r_coonod)

!--- initialize fields
    CALL initial_values(p_param, i_numelmt, i_elmtdofs, i_elmtnodes, r_coodof, r_coonod, r_Q)

!--- update grid information
    CALL grid_putinfo(p_ghand, l_finelevel=.TRUE., i_arraypoint=i_valQ, r_dofvalues=r_Q)

!--- deallocate workspace
    DEALLOCATE(i_elmtdofs, i_elmtnodes, r_coodof, r_coonod, r_Q)

  END SUBROUTINE initialvalues_iter

!*******************************************************************************
! DESCRIPTION of [SUBROUTINE source_update]:
!> @brief updates the source terms
!>
!> @param[in]       p_param       control structure for global parameters
!> @param[in]       r_time        model time at which source term is evaluated
!> @param[in]       i_numelmt     number of elements
!> @param[in]       i_elmtdofs
!> @param           r_coodof      coordinates of DOFs
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

    r_sol = 0.0_GRID_SR

  END SUBROUTINE dg_elmt_solution

!*******************************************************************************
! Only internal subroutines follow
!*******************************************************************************
! DESCRIPTION of [SUBROUTINE initial_values]:
!> @brief calculates the warm air bubble as in Giraldo and Restelli (2007)
!>
!> @param[in,out]   p_param       control structure for global parameters
!> @param[in]       i_elmtdofs
!> @param[in]       r_coodof      coordinates of DOFs
!> @param[out]      r_Q           state vector
!>
!> @note the initial condition are taken from Giraldo and Restelli (2007)
!
  SUBROUTINE initial_values(p_param, i_numelmt, i_elmtdofs, i_elmtnodes, r_coodof, r_coonod, r_Q)

    IMPLICIT NONE

    TYPE (control_struct),                    INTENT(inout) :: p_param
    INTEGER (KIND = GRID_SI),                 INTENT(in)    :: i_numelmt
    INTEGER (KIND = GRID_SI), DIMENSION(:,:), INTENT(in)    :: i_elmtdofs, i_elmtnodes
    REAL (KIND = GRID_SR), DIMENSION(:,:),    INTENT(in)    :: r_coodof, r_coonod
    REAL (KIND = GRID_SR), DIMENSION(:,:),    INTENT(inout) :: r_Q

!--- local declarations
    INTEGER (KIND = GRID_SI)                                :: i_faceunknowns, i_dof
    REAL (KIND = GRID_SR)                                   :: r_bbleradius, &
      r_rho, r_thetabg, r_theta, &
      r_exner, r_radius, r_thetac, r_gammaaux, r_cp, r_pressure
    REAL (KIND = GRID_SR), DIMENSION(2)                     :: r_bblecntr, r_velo0

!--- initialize some constants
    i_faceunknowns = GRID_femtypes%p_type(FEM_DG)%sig%i_unknowns
    r_gammaaux = 1.0_GRID_SR / (p_equationparam%r_gamma - 1.0_GRID_SR) ! auxiliary variable
    r_thetabg  = p_testparam%treal(1)%p_value(1)                       ! background potential temperature
    r_thetac   = p_testparam%treal(2)%p_value(1)                       ! max. amplitude of potential temperature deviation
    r_velo0    = p_testparam%treal(3)%p_value                          ! initial velocity
    r_bblecntr = p_testparam%treal(4)%p_value                          ! center of the bubble
    r_bbleradius = p_testparam%treal(5)%p_value(1)                     ! radius around center of bubble that has constant deviation
    r_cp       = p_equationparam%r_gamma / &
      (p_equationparam%r_gamma-1.0_GRID_SR) * p_equationparam%r_rgas   ! specific heat capacity at constant pressure

!--- add perturbation to THETA and calculate variables
    DO i_dof = 1, i_faceunknowns*i_numelmt
      r_radius = SQRT( (r_coodof(1,i_dof)-r_bblecntr(1))**2 + &
                       (r_coodof(2,i_dof)-r_bblecntr(2))**2 )

!--- calculate potential temperature from deviation + background
      IF (r_radius <= r_bbleradius) THEN
        r_theta = r_thetabg + r_thetac / 2.0_GRID_SR * &
                  (1.0_GRID_SR + COS(GRID_PI * r_radius / r_bbleradius))
      ELSE
        r_theta = r_thetabg
      END IF
      r_exner    = 1.0_GRID_SR - GRID_GRAV/(r_cp*r_thetabg)*r_coodof(2,i_dof) ! this is the deviation of exner pressure
      r_rho      = p_equationparam%r_pref / &
        (p_equationparam%r_rgas*r_theta) * (r_exner)**(r_gammaaux)            ! this is the density
      r_pressure = p_equationparam%r_pref * &
        (p_equationparam%r_rgas * r_rho * r_theta/p_equationparam%r_pref)**p_equationparam%r_gamma !- GRID_GRAV * r_rho * r_coodof(2,i_dof)

!--- calculate total energy if total energy form is used
      CALL PrimitiveToProgvars([r_rho, r_velo0, r_pressure], r_Q(:,i_dof))

!--- r_Q are the conserved variables with
!  r_Q(1,i_dof) = r_rho           being rho
!  r_Q(2,i_dof) = r_rho*r_velo(1) being rho * u
!  r_Q(3,i_dof) = r_rho*r_velo(2) being rho * v
!  r_Q(4,i_dof) = r_rho*r_theta   being rho * theta

    END DO !i_dof

  END SUBROUTINE initial_values

!*******************************************************************************
END MODULE DG_initial
