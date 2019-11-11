!*******************************************************************************
!
!> @file  DG_initial_storm.f90
!> @brief contains module DG_initial
!
!*******************************************************************************
! MODULE DESCRIPTION:
!> @brief Initializes fields for a simple storm surge testcase (see Mandli 2011).
!>        A synthetic storm travels linearly from r_pos_start to r_pos_final
!>        within a time r_tfinal
!
MODULE DG_initial

  USE GRID_api
  USE FLASH_parameters
  USE DG_storm_holland
  USE DG_equation, ONLY : i_nprogvars, i_nsrcterms, i_valQ, i_valS

  PRIVATE
  PUBLIC  i_ntstintparam, i_ntstrealparam, i_ntstcharparam, i_ntstlogparam, &
          initial_setparam, p_testparam, initialize_testcase, initialvalues_iter, &
          source_update, dg_elmt_solution

  TYPE (test_param)                                   :: p_testparam
  INTEGER, PARAMETER                                  :: i_ntstintparam  = 0
  INTEGER, PARAMETER                                  :: i_ntstrealparam = 13
  INTEGER, PARAMETER                                  :: i_ntstcharparam = 0
  INTEGER, PARAMETER                                  :: i_ntstlogparam  = 0

  CONTAINS

!*******************************************************************************
! DESCRIPTION of [SUBROUTINE initial_setparam]:
!> @brief initializes test case specific parameter list
!
  SUBROUTINE initial_setparam

    IMPLICIT NONE

    p_testparam%treal( 1)%c_keyword = 'ALPHA1'
    p_testparam%treal( 1)%i_size    = 1

    p_testparam%treal( 2)%c_keyword = 'ALPHA2'
    p_testparam%treal( 2)%i_size    = 1

    p_testparam%treal( 3)%c_keyword = 'PRESSURE_CENT'
    p_testparam%treal( 3)%i_size    = 1

    p_testparam%treal( 4)%c_keyword = 'PRESSURE_N'
    p_testparam%treal( 4)%i_size    = 1

    p_testparam%treal( 5)%c_keyword = 'PRESSURE_A'
    p_testparam%treal( 5)%i_size    = 1

    p_testparam%treal( 6)%c_keyword = 'PARAM_A'
    p_testparam%treal( 6)%i_size    = 1

    p_testparam%treal( 7)%c_keyword = 'PARAM_B'
    p_testparam%treal( 7)%i_size    = 1

    p_testparam%treal( 8)%c_keyword = 'DENSITY_AIR'
    p_testparam%treal( 8)%i_size    = 1

    p_testparam%treal( 9)%c_keyword = 'START_POSITION'
    p_testparam%treal( 9)%i_size    = 2

    p_testparam%treal(10)%c_keyword = 'FINAL_POSITION'
    p_testparam%treal(10)%i_size    = 2

    p_testparam%treal(11)%c_keyword = 'START_TIME'
    p_testparam%treal(11)%i_size    = 1

    p_testparam%treal(12)%c_keyword = 'FINAL_TIME'
    p_testparam%treal(12)%i_size    = 1

    p_testparam%treal(13)%c_keyword = 'RAMP_TIME'
    p_testparam%treal(13)%i_size    = 1

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

!--- initialize fields, update source
    r_Q = 0._GRID_SR
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
    REAL(KIND = GRID_SR)                                    :: r_alpha1, r_alpha2, &
      r_Pc, r_Pn, r_PA, r_A, r_B, r_rho_air, r_tstart, r_tfinal, r_tramp
    REAL(KIND = GRID_SR), DIMENSION(2)                      :: r_pos_start, r_pos_final
    REAL (KIND = GRID_SR)                                   :: r_aux
    REAL (KIND = GRID_SR), DIMENSION(2)                     :: r_pos_storm
    REAL (KIND = GRID_SR), DIMENSION(2)                     :: r_windfield

    r_alpha1    = p_testparam%treal( 1)%p_value(1)
    r_alpha2    = p_testparam%treal( 2)%p_value(1)
    r_Pc        = p_testparam%treal( 3)%p_value(1)
    r_Pn        = p_testparam%treal( 4)%p_value(1)
    r_PA        = p_testparam%treal( 5)%p_value(1)
    r_A         = p_testparam%treal( 6)%p_value(1)
    r_B         = p_testparam%treal( 7)%p_value(1)
    r_rho_air   = p_testparam%treal( 8)%p_value(1)
    r_pos_start = p_testparam%treal( 9)%p_value
    r_pos_final = p_testparam%treal(10)%p_value
    r_tstart    = p_testparam%treal(11)%p_value(1)
    r_tfinal    = p_testparam%treal(12)%p_value(1)
    r_tramp     = p_testparam%treal(13)%p_value(1)
    r_pos_storm = 0._GRID_SR

    !IF (l_gridchanged) THEN
      !--- initialize some constants
      i_faceunknowns = GRID_femtypes%p_type(FEM_DG)%sig%i_unknowns

      !--- get the position of the hurricane; this includes a ramping up of 2h
      CALL get_storm_position(r_time, r_tstart, r_tfinal, r_tramp, &
                              r_pos_storm, r_pos_start, r_pos_final)

      !--- switch off all forcing after r_tfinal (plus the ramp up)
      IF (r_time .GT. r_tfinal + r_tramp) THEN
        r_S(3,:) = 0._GRID_SR
        r_S(4,:) = 0._GRID_SR
      ELSE
        !--- loop over the elements
        DO i_elmt=1, i_numelmt
          DO i_dof=1, i_faceunknowns
            i_globaldof = i_elmtdofs(i_dof,i_elmt)
            r_aux = r_Q(1,i_globaldof) + r_S(1,i_globaldof)
            r_windfield = 0._GRID_SR

            !--- compute velocity field using Holland's model
            CALL get_wind_holland(r_windfield, r_pos_storm, r_Pn, r_Pc, r_rho_air, r_A, r_B, &
                                  r_coodof(1,i_globaldof), r_coodof(2,i_globaldof), r_time, r_tramp)

            r_S(3,i_globaldof) = r_windfield(1)
            r_S(4,i_globaldof) = r_windfield(2)

            ! -- Only for t=0
            IF (r_coodof(1,i_globaldof) .LE. 350000.0) THEN
              r_S(1,i_globaldof) = 0._GRID_SR
            ELSE
              r_S(1,i_globaldof) = 0.0125_GRID_SR*(r_coodof(1,i_globaldof)-350000.0_GRID_SR) * GRID_GRAV
            END IF

            r_Q(1,i_globaldof) = MAX(0.0_GRID_SR, r_aux - r_S(1,i_globaldof))

          END DO
        END DO
      END IF
    !END IF !l_gridchanged

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
    INTEGER (KIND = GRID_SI)                                :: i_count, i_elmt, i_faceunknowns
    INTEGER (KIND = GRID_SI)                                :: i_globaldof
    REAL (KIND = GRID_SR)                                   :: r_x, r_bathy
    REAL (KIND = GRID_SR)                                   :: r_alpha1, r_alpha2

    r_alpha1    = p_testparam%treal(1)%p_value(1)
    r_alpha2    = p_testparam%treal(2)%p_value(1)

!--- allocate workspace
    i_faceunknowns = GRID_femtypes%p_type(FEM_DG)%sig%i_unknowns

    DO i_elmt=1, i_numelmt
      DO i_count=1,i_faceunknowns

        i_globaldof = i_elmtdofs(i_count,i_elmt)
        r_x = r_coodof(1,i_globaldof)

        IF (r_x .LE. 350000.0) THEN
          r_bathy = 0._GRID_SR
          !       ELSEIF ((r_x .GT. 350000.0) .AND. (r_x .LE. 450000.0)) THEN
          !         r_bathy = (r_x - 350000.0)*r_alpha1
          !       ELSEIF ((r_x .GT. 450000.0) .AND. (r_x .LE. 475000.0)) THEN
          !         r_bathy = 2700.0
          !       ELSEIF ( r_x .GT. 475000.0) THEN
          !         r_bathy= 2700.0+r_alpha2*(r_x - 475000.0)
        ELSE
          r_bathy= 0.0125_GRID_SR*(r_x-350000.0_GRID_SR)
        END IF

        r_Q(1,i_globaldof)   = MAX(3000.0_GRID_SR-r_bathy, 0.0_GRID_SR)*GRID_GRAV
        r_Q(2:3,i_globaldof) = 0.0_GRID_SR

      END DO !i_count
    END DO !i_elmt

  END SUBROUTINE field_init

!*******************************************************************************
END MODULE DG_initial
