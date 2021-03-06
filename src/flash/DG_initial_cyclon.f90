!*******************************************************************************
!
!> @file  DG_initial_cyclon.f90
!> @brief contains module DG_initial
!
!*******************************************************************************
! MODULE DESCRIPTION:
!> @brief initialize velocity field with linear advection and SUD vortex
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
  INTEGER, PARAMETER                                  :: i_ntstrealparam = 3
  INTEGER, PARAMETER                                  :: i_ntstcharparam = 0
  INTEGER, PARAMETER                                  :: i_ntstlogparam  = 0

  CONTAINS

!*******************************************************************************
! DESCRIPTION of [SUBROUTINE initial_setparam]:
!> @brief initializes test case specific parameter list
!
  SUBROUTINE initial_setparam

    IMPLICIT NONE

    ! initial position of vortex
    p_testparam%treal(1)%c_keyword = 'CENTER_POINT'
    p_testparam%treal(1)%i_size    = 2

    ! angle between direction of vortex and front
    p_testparam%treal(2)%c_keyword = 'ANGLE_TO_FRONT'
    p_testparam%treal(2)%i_size    = 1

    ! velocity of the backround drift
    p_testparam%treal(3)%c_keyword = 'DRIFT_VELOCITY'
    p_testparam%treal(3)%i_size    = 1

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
!> @note This routine initializes a SUD vortex on top of a background drift.
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
    CALL add_2suds(p_param, i_numelmt, i_elmtdofs, r_coodof, r_Q)

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
!> @param           r_S             discrete fields used in source terms
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
          r_S(1,i_elmtdofs(i_dof,i_elmt)) = 0.0_GRID_SR * GRID_GRAV          ! bathymetry
!           r_S(2,i_elmtdofs(i_dof,i_elmt)) = 2.0_GRID_SR * GRID_OMEGA * sin(r_coodof(2,i_cnt))
          r_S(2,i_elmtdofs(i_dof,i_elmt)) = 0.0_GRID_SR                      ! Coriolis term
          r_S(3,i_elmtdofs(i_dof,i_elmt)) = 0.0_GRID_SR                      ! wind in x-direction
          r_S(4,i_elmtdofs(i_dof,i_elmt)) = 0.0_GRID_SR                      ! wind in y-direction
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
    REAL (KIND = GRID_SR)                                         :: r_mid, r_L, r_alpha
    REAL (KIND = GRID_SR)                                         :: r_x, r_y

!--- initialize some constants
    r_mid   = 15.0_GRID_SR
    r_alpha = 0.25_GRID_SR
    r_L     = 10.0_GRID_SR

    r_x = r_coo(1)
    r_y = r_coo(2)

!    write(*,*) r_x, r_y, r_mid, r_alpha
    r_sol(1) = r_L*(1.0_GRID_SR-0.5_GRID_SR*EXP(-r_alpha*((r_x-r_mid)**2+(r_y-r_mid)**2)))
    r_sol(2) = 0.0_GRID_SR
    r_sol(3) = 0.0_GRID_SR

  END FUNCTION exact_solution

!*******************************************************************************
! DESCRIPTION of [SUBROUTINE field_init]:
!> @brief initializes the height and momentum fields
!>
!> @param[in,out]   p_param       control structure for global parameters
!> @param           r_coodof      coordinates of DOFs
!
  SUBROUTINE field_init(p_param, i_numelmt, i_elmtdofs, r_coodof, r_Q)

    IMPLICIT NONE

    TYPE (control_struct),                    INTENT(inout) :: p_param
    INTEGER (KIND = GRID_SI),                 INTENT(in)    :: i_numelmt
    INTEGER (KIND = GRID_SI), DIMENSION(:,:), INTENT(in)    :: i_elmtdofs
    REAL (KIND = GRID_SR), DIMENSION(:,:),    INTENT(in)    :: r_coodof
    REAL (KIND = GRID_SR), DIMENSION(:,:),    INTENT(inout) :: r_Q

!--- local declarations
    INTEGER  (KIND = GRID_SI)                               :: i_count, i_num, &
                                                               i_elmt, i_faceunknowns
    REAL (KIND = GRID_SR)                                   :: r_alpha, r_pi, r_velo

!--- initialize some constants
    r_pi    = 4.0*atan(1.0)
    r_alpha = p_testparam%treal(2)%p_value(1) * r_pi/180.0
    r_velo  = p_testparam%treal(3)%p_value(1)

    i_faceunknowns =  GRID_femtypes%p_type(FEM_DG)%sig%i_unknowns
    i_num          = i_numelmt*i_faceunknowns

!--- loop over the DOFs
    DO i_count= 1, i_num
      r_Q(1,i_count) = 10000.0_GRID_SR
      r_Q(2,i_count) = r_velo      !r_velo*sin(r_alpha)
      r_Q(3,i_count) = 0.0_GRID_SR !r_velo*cos(r_alpha)
    END DO

  END SUBROUTINE field_init

!*******************************************************************************
! DESCRIPTION of [SUBROUTINE add_sud]:
!> @brief initializes the velocity field
!>
!> @param           r_coodof      coordinates of DOFs
!>
!> @note A SUD vortex centered at the location specified in the parameters
!>       file is created.
!>
  SUBROUTINE add_sud(p_param, i_numelmt, i_elmtdofs, r_coodof, r_Q)

    IMPLICIT NONE

    TYPE (control_struct),                    INTENT(inout) :: p_param
    INTEGER (KIND = GRID_SI),                 INTENT(in)    :: i_numelmt
    INTEGER (KIND = GRID_SI), DIMENSION(:,:), INTENT(in)    :: i_elmtdofs
    REAL (KIND = GRID_SR), DIMENSION(:,:),    INTENT(in)    :: r_coodof
    REAL (KIND = GRID_SR), DIMENSION(:,:),    INTENT(inout) :: r_Q

    !--- local declarations
    INTEGER (KIND = GRID_SI)                                :: i_num, i_count
    REAL                                                    :: r_x, r_y    !< x,y-coordinates of a DOF [m]
    REAL                                                    :: r_radius    !< Fixed radius for given tangential velocity [m]
    REAL                                                    :: r_velocity  !< Given tangential velocity at r_radius [m/s]
    REAL                                                    :: r_x_init    !< Initial x coordinate of the cyclon center [m]
    REAL                                                    :: r_y_init    !< Initial y coordiante of the cyclon center [m]
    REAL                                                    :: r_r         !< Distance between DOF and cyclon center [m]
    REAL                                                    :: r_alpha     !< Angel between gripoint and cyclon center [radian]
    REAL                                                    :: r_s         !< Relation between r_radius and r_r []
    REAL                                                    :: r_v_tan     !< Tangential velocity [m/s]
    REAL                                                    :: r_a, r_b    !< Parameters for unsymmetric behaviour []
    REAL                                                    :: r_choke     !< Scaling parameter [0,1], s.t. velocity at r_rcut is zero []
    REAL                                                    :: r_rcut      !< Radius limit to cut the initial velocity function to zero [m]

!--- initialize some constants
    i_num = i_numelmt*GRID_femtypes%p_type(FEM_DG)%sig%i_unknowns

    r_x_init   = p_testparam%treal(1)%p_value(1) * (10.0**3)
    r_y_init   = p_testparam%treal(1)%p_value(2) * (10.0**3)
!--- Fixed SUD parameters
    r_radius   = 100000.    ! radius for given tangential velocity
    r_velocity = 71.521     ! tangential velocity at r_radius
    r_a        = 0.3398     ! SUD specific parameter
    r_b        = 0.0005377  ! SUD specific parameter
    r_rcut     = 7500000.   ! Radius limit to cut the initial velocity function to zero

    node_loop: DO i_count= 1, i_num
      r_x = r_coodof(1,i_count)
      r_y = r_coodof(2,i_count)
      r_r = sqrt((r_x-r_x_init)*(r_x-r_x_init) + (r_y-r_y_init)*(r_y-r_y_init))
      IF (r_r <= r_rcut) THEN
        r_choke = 1.-exp(-((r_r-r_rcut)**2.0)/(100000.**2.0))
      ELSE
        r_choke = 0.
      END IF
      r_s = r_r/r_radius
      r_alpha = atan2((r_y-r_y_init),(r_x-r_x_init))
      r_v_tan =  r_velocity * (r_s*(1+(6*r_b/(2*r_a))*(r_s**4))) / ((1+r_a*r_s*r_s+r_b*(r_s**6))**2)*r_choke
      r_Q(2,i_count) = -r_v_tan*sin(r_alpha) + r_Q(2,i_count)
      r_Q(3,i_count) =  r_v_tan*cos(r_alpha) + r_Q(3,i_count)
    END DO node_loop

  END SUBROUTINE add_sud

!*******************************************************************************
! DESCRIPTION of [SUBROUTINE add_2suds]:
!> @brief initializes the velocity field
!>
!> @param           r_coodof      coordinates of DOFs
!>
!> @note A SUD vortex centered at the location specified in the parameters
!>       file is created and a second vortex which is shifted is also inizialized.
!
  SUBROUTINE add_2suds(p_param, i_numelmt, i_elmtdofs, r_coodof, r_Q)

    IMPLICIT NONE

    TYPE (control_struct),                    INTENT(inout) :: p_param
    INTEGER (KIND = GRID_SI),                 INTENT(in)    :: i_numelmt
    INTEGER (KIND = GRID_SI), DIMENSION(:,:), INTENT(in)    :: i_elmtdofs
    REAL (KIND = GRID_SR), DIMENSION(:,:),    INTENT(in)    :: r_coodof
    REAL (KIND = GRID_SR), DIMENSION(:,:),    INTENT(inout) :: r_Q

    !--- local declarations
    INTEGER (KIND = GRID_SI)                                :: i_num, i_count
    INTEGER                                                 :: i_sud       !< Temporary counter for numer of suds to add
    REAL                                                    :: r_x, r_y    !< x,y-coordinates of a DOF [m]
    REAL                                                    :: r_radius    !< Fixed radius for given tangential velocity [m]
    REAL                                                    :: r_velocity  !< Given tangential velocity at r_radius [m/s]
    REAL                                                    :: r_x_init    !< Initial x coordinate of the cyclon center [m]
    REAL                                                    :: r_y_init    !< Initial y coordiante of the cyclon center [m]
    REAL                                                    :: r_x_dist    !< Distance in x-direction of sud centers [m]
    REAL                                                    :: r_y_dist    !< Distnace in y-direction of sud centers [m]
    REAL                                                    :: r_r         !< Distance between DOF and cyclon center [m]
    REAL                                                    :: r_alpha     !< Angel between gripoint and cyclon center [radian]
    REAL                                                    :: r_s         !< Relation between r_radius and r_r []
    REAL                                                    :: r_v_tan     !< Tangential velocity [m/s]
    REAL                                                    :: r_a, r_b    !< Parameters for unsymmetric behaviour []
    REAL                                                    :: r_choke     !< Scaling parameter [0,1], s.t. velocity at r_rcut is zero []
    REAL                                                    :: r_rcut      !< Radius limit to cut the initial velocity function to zero [m]

!--- allocate workspace
    i_num = i_numelmt*GRID_femtypes%p_type(FEM_DG)%sig%i_unknowns

    r_x_init   = p_testparam%treal(1)%p_value(1) * (10.0**3)
    r_y_init   = p_testparam%treal(1)%p_value(2) * (10.0**3)
!--- Fixed SUD parameters
    r_radius         = 100000.             ! radius for given tangential velocity
    r_velocity       = 15.0 !71.521        ! tangential velocity at r_radius
    r_a              = 0.3398              ! SUD specific parameter
    r_b              = 0.0005377           ! SUD specific parameter
    r_rcut           = 7500000.            ! Radius limit to cut the initial velocity function to zero
    r_x_dist         = 0.                  ! Distance of centers of vortexes in x direction
    r_y_dist         = 350000.             ! Distance of centers of vortexes in x direction

    sud_loop: DO i_sud=1,2
      node_loop: DO i_count= 1, i_num
        r_x = r_coodof(1,i_count)
        r_y = r_coodof(2,i_count)
        r_r = sqrt((r_x-r_x_init)*(r_x-r_x_init) + (r_y-r_y_init)*(r_y-r_y_init))
        IF (r_r <= r_rcut) THEN
          r_choke = 1.-exp(-((r_r-r_rcut)**2.0)/(100000.**2.0))
        ELSE
          r_choke = 0.
        END IF
        r_s = r_r/r_radius
        r_alpha = atan2((r_y-r_y_init),(r_x-r_x_init))
        r_v_tan =  r_velocity * (r_s*(1+(6*r_b/(2*r_a))*(r_s**4))) / ((1+r_a*r_s*r_s+r_b*(r_s**6))**2)*r_choke
        r_Q(2,i_count) = -r_v_tan*sin(r_alpha) + r_Q(2,i_count)
        r_Q(3,i_count) =  r_v_tan*cos(r_alpha) + r_Q(3,i_count)
      END DO node_loop

      r_x_init = r_x_init+r_x_dist
      r_y_init = r_y_init+r_y_dist
      r_radius = 100000.0
      r_velocity = 15.0 !71.521
    END DO sud_loop

  END SUBROUTINE add_2suds

!*******************************************************************************
END MODULE DG_initial
