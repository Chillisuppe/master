!*******************************************************************************
!
!> @file  DG_initial_vortex_benacchio.f90
!> @brief contains module DG_initial
!
!*******************************************************************************
! MODULE DESCRIPTION:
!> @brief initialize fields with a quasi-stationary vortex, see Benacchio (2014)
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
  INTEGER, PARAMETER                                  :: i_ntstrealparam = 4
  INTEGER, PARAMETER                                  :: i_ntstcharparam = 0
  INTEGER, PARAMETER                                  :: i_ntstlogparam  = 0

  CONTAINS

!*******************************************************************************
! DESCRIPTION of [SUBROUTINE initial_setparam]:
!> @brief initializes test case specific parameter list
!
  SUBROUTINE initial_setparam

    IMPLICIT NONE

    ! initial density
    p_testparam%treal(1)%c_keyword = 'INITIAL_DENSITY'
    p_testparam%treal(1)%i_size    = 1

    ! initial velocity (u0, v0)
    p_testparam%treal(2)%c_keyword = 'INITIAL_VELO'
    p_testparam%treal(2)%i_size    = 2

    ! vortex center position
    p_testparam%treal(3)%c_keyword = 'VORTEX_CENTER'
    p_testparam%treal(3)%i_size    = 2

    ! vortex radius
    p_testparam%treal(4)%c_keyword = 'VORTEX_RADIUS'
    p_testparam%treal(4)%i_size    = 1

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
!> @brief calculates the smooth vortex as used in Benacchio 2014 (dissertation)
!>
!> @param[in,out]   p_param       control structure for global parameters
!> @param[in]       i_elmtdofs
!> @param[in]       r_coodof      coordinates of DOFs
!> @param[out]      r_Q           state vector
!>
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
    REAL (KIND = GRID_SR)                                   :: r_theta, &
      r_rho0, r_radius, r_vortradius, r_rho, r_phase, r_pressure
    REAL (KIND = GRID_SR), DIMENSION(2)                     :: r_vortcntr, r_velo0, r_velo

!--- initialize some constants
    i_faceunknowns = GRID_femtypes%p_type(FEM_DG)%sig%i_unknowns
    r_rho0       = p_testparam%treal(1)%p_value(1)                ! initial gas density
    r_velo0      = p_testparam%treal(2)%p_value                   ! initial velocity
    r_vortcntr   = p_testparam%treal(3)%p_value                   ! center of the vortex
    r_vortradius = p_testparam%treal(4)%p_value(1)                ! vortex radius

!--- add perturbation to THETA and calculate variables
    DO i_dof = 1, i_faceunknowns*i_numelmt
      r_radius = SQRT( (r_coodof(1,i_dof)-r_vortcntr(1))**2 + &
                       (r_coodof(2,i_dof)-r_vortcntr(2))**2 ) / r_vortradius
      r_phase = ATAN2(r_coodof(2,i_dof)-r_vortcntr(2), r_coodof(1,i_dof)-r_vortcntr(1))

      IF (r_radius < 1.0_GRID_SR) THEN
        r_rho      = r_rho0 + 0.5_GRID_SR*(1.0_GRID_SR-r_radius**2)**6
        r_velo(1)  = r_velo0(1) - 1024*SIN(r_phase)*(1.0_GRID_SR-r_radius)**6*r_radius**6
        r_velo(2)  = r_velo0(2) + 1024*COS(r_phase)*(1.0_GRID_SR-r_radius)**6*r_radius**6
        r_pressure = init_pressure(r_radius,r_rho0)*r_rho0*(r_velo0(1)**2+r_velo0(2)**2) + &
                     p_equationparam%r_pref
      ELSE
        r_rho      = r_rho0
        r_velo(1)  = r_velo0(1)
        r_velo(2)  = r_velo0(2)
        r_pressure = p_equationparam%r_pref
      END IF

!--- calculate total energy if total energy form is used
      CALL PrimitiveToProgvars([r_rho, r_velo, r_pressure], r_Q(:,i_dof))

    END DO !i_dof

  END SUBROUTINE initial_values

!*******************************************************************************
! DESCRIPTION of [FUNCTION init_pressure]:
!> @brief pressure for the vortex
!>
!> @param[in]       r_r          radius r
!> @param[in]       r_d          center density
!> @return          r_pc         pressure coefficient
!
  FUNCTION init_pressure(r_r, r_d) RESULT(r_pc)

    IMPLICIT NONE

    REAL (KIND = GRID_SR), INTENT(IN)     :: r_r, r_d
    REAL (KIND = GRID_SR)                 :: r_pc

!--- local declarations
    REAL (KIND = GRID_SR)                 :: r_rdl
    REAL (KIND = GRID_SR), DIMENSION(25)  :: r_coeff, r_rpow, r_ones
    INTEGER (KIND = GRID_SI)              :: i_cnt

!--- these are powers of the radius r (polynomials)
    DO i_cnt = 1, 25
      r_rpow(i_cnt) = r_r**(i_cnt-1)
    END DO
    r_ones = 1.0_GRID_SR

!--- this is rho_dl
    r_rdl = (1.0_GRID_SR + 0.5_GRID_SR*(1-r_r**2)**6)/r_d

!--- these are the coefficients that get multiplied to the powers of r

    IF (r_d /= 0.5_GRID_SR) THEN
      r_coeff( 1) =   (   1.0_GRID_SR +   2.0_GRID_SR*r_rdl)                / 24.0_GRID_SR
      r_coeff( 2) = - (   1.0_GRID_SR +   2.0_GRID_SR*r_rdl) * 6.0_GRID_SR  / 13.0_GRID_SR
      r_coeff( 3) =   (   5.0_GRID_SR +  11.0_GRID_SR*r_rdl) * 3.0_GRID_SR  /  7.0_GRID_SR
      r_coeff( 4) = - (  37.0_GRID_SR + 110.0_GRID_SR*r_rdl) * 2.0_GRID_SR  / 15.0_GRID_SR
      r_coeff( 5) =   (  19.0_GRID_SR + 165.0_GRID_SR*r_rdl) * 3.0_GRID_SR  / 16.0_GRID_SR
      r_coeff( 6) =   (  29.0_GRID_SR - 132.0_GRID_SR*r_rdl) * 6.0_GRID_SR  / 17.0_GRID_SR
      r_coeff( 7) = - ( 269.0_GRID_SR - 462.0_GRID_SR*r_rdl)                /  9.0_GRID_SR
      r_coeff( 8) =   (  25.0_GRID_SR -  44.0_GRID_SR*r_rdl) * 18.0_GRID_SR / 19.0_GRID_SR
      r_coeff( 9) =   ( 119.0_GRID_SR + 110.0_GRID_SR*r_rdl) * 9.0_GRID_SR  / 40.0_GRID_SR
      r_coeff(10) = - ( 391.0_GRID_SR +  55.0_GRID_SR*r_rdl) * 4.0_GRID_SR  / 21.0_GRID_SR
      r_coeff(11) =   ( 510.0_GRID_SR/11.0_GRID_SR + 3.0_GRID_SR*r_rdl)
      r_coeff(12) =   (  85.0_GRID_SR -               r_rdl) * 12.0_GRID_SR / 23.0_GRID_SR
      r_coeff(13) = - (2210.0_GRID_SR -               r_rdl)                / 24.0_GRID_SR
    ELSE
      r_coeff( 1) =      1.0_GRID_SR / 12.0_GRID_SR
      r_coeff( 2) = -   12.0_GRID_SR / 13.0_GRID_SR
      r_coeff( 3) =      9.0_GRID_SR /  2.0_GRID_SR
      r_coeff( 4) = -  184.0_GRID_SR / 15.0_GRID_SR
      r_coeff( 5) =    609.0_GRID_SR / 32.0_GRID_SR
      r_coeff( 6) = -  222.0_GRID_SR / 17.0_GRID_SR
      r_coeff( 7) = -   38.0_GRID_SR /  9.0_GRID_SR
      r_coeff( 8) =     54.0_GRID_SR / 19.0_GRID_SR
      r_coeff( 9) =    783.0_GRID_SR / 20.0_GRID_SR
      r_coeff(10) = -  558.0_GRID_SR /  7.0_GRID_SR
      r_coeff(11) =   1053.0_GRID_SR / 22.0_GRID_SR
      r_coeff(12) =   1014.0_GRID_SR / 23.0_GRID_SR
      r_coeff(13) = - 1473.0_GRID_SR / 16.0_GRID_SR
    END IF

    r_coeff(14) =    204.0_GRID_SR /  5.0_GRID_SR
    r_coeff(15) =    510.0_GRID_SR / 13.0_GRID_SR
    r_coeff(16) = - 1564.0_GRID_SR / 27.0_GRID_SR
    r_coeff(17) =    153.0_GRID_SR /  8.0_GRID_SR
    r_coeff(18) =    450.0_GRID_SR / 29.0_GRID_SR
    r_coeff(19) = -  269.0_GRID_SR / 15.0_GRID_SR
    r_coeff(20) =    174.0_GRID_SR / 31.0_GRID_SR
    r_coeff(21) =     57.0_GRID_SR / 32.0_GRID_SR
    r_coeff(22) = -   74.0_GRID_SR / 33.0_GRID_SR
    r_coeff(23) =     15.0_GRID_SR / 17.0_GRID_SR
    r_coeff(24) = -    6.0_GRID_SR / 35.0_GRID_SR
    r_coeff(25) =      1.0_GRID_SR / 72.0_GRID_SR

    r_pc = 1024**2 * (r_r**12 * DOT_PRODUCT(r_rpow,r_coeff) - DOT_PRODUCT(r_ones,r_coeff))

  END FUNCTION init_pressure

!*******************************************************************************
END MODULE DG_initial
