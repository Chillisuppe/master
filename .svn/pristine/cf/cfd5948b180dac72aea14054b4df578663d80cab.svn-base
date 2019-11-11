!*******************************************************************************
!
!> @file  DG_initial_riemannproblem.f90
!> @brief contains module DG_initial
!
!*******************************************************************************
! MODULE DESCRIPTION:
!> @brief initializes fields with Riemann problem, and computes exact solution
!
MODULE DG_initial

  USE GRID_api
  USE FLASH_parameters
  USE DG_equation, ONLY : i_nprogvars, i_nsrcterms, i_valQ, i_valS
  USE MISC_quad

  PRIVATE swguess, swfunc1, swfunc2
  PUBLIC  i_ntstintparam, i_ntstrealparam, i_ntstcharparam, i_ntstlogparam, &
          initial_setparam, p_testparam, initialize_testcase, initialvalues_iter, &
          source_update, dg_elmt_solution

  TYPE (test_param)                                   :: p_testparam
  INTEGER, PARAMETER                                  :: i_ntstintparam  = 0
  INTEGER, PARAMETER                                  :: i_ntstrealparam = 3
  INTEGER, PARAMETER                                  :: i_ntstcharparam = 0
  INTEGER, PARAMETER                                  :: i_ntstlogparam  = 0

  TYPE (quadrature)                                             :: p_quad

  CONTAINS

!*******************************************************************************
! DESCRIPTION of [SUBROUTINE initial_setparam]:
!> @brief initializes test case specific parameter list
!
  SUBROUTINE initial_setparam

    IMPLICIT NONE

    p_testparam%treal(1)%c_keyword = 'XORIGIN'
    p_testparam%treal(1)%i_size    = 1

    p_testparam%treal(2)%c_keyword = 'STATE_LEFT'
    p_testparam%treal(2)%i_size    = 2

    p_testparam%treal(3)%c_keyword = 'STATE_RIGHT'
    p_testparam%treal(3)%i_size    = 2

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

!--- initialize quadrature (this creates p_quad with values)
    CALL p_quad%initialize(2*GRID_femtypes%p_type(FEM_DG)%sig%i_degree, FEM_DG)

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
             r_Q(i_nprogvars, i_num), r_S(i_nsrcterms, i_num), STAT=i_alct)
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
  SUBROUTINE dg_elmt_solution(r_coodof, r_time, i_numelmt, &
                              i_elmtdofs, r_coonod, i_elmtnodes, r_sol)

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

!--- compute L2 projection
    CALL p_quad%L2_projection(r_coonod, i_elmtnodes, i_elmtdofs, FUNCh , r_sol(1,:))
    CALL p_quad%L2_projection(r_coonod, i_elmtnodes, i_elmtdofs, FUNCmx, r_sol(2,:))
    CALL p_quad%L2_projection(r_coonod, i_elmtnodes, i_elmtdofs, FUNCmy, r_sol(3,:))

! !--- loop over the elements and DOFs
!     DO i_elmt=1, i_numelmt
!       DO i_dof=1, i_faceunknowns
!         r_sol(:,i_elmtdofs(i_dof,i_elmt)) = &
!           exact_solution(r_coodof(:,i_elmtdofs(i_dof,i_elmt)), r_time)
!       END DO
!     END DO

    CONTAINS

    FUNCTION FUNCh(x,y) RESULT(z)
      REAL (KIND = GRID_SR), INTENT(IN)   :: x, y
      REAL (KIND = GRID_SR)               :: z

      REAL (KIND = GRID_SR), DIMENSION(3) :: r_q

      r_q = exact_solution([x, y], r_time)
      z = r_q(1)

    END FUNCTION FUNCh

    FUNCTION FUNCmx(x,y) RESULT(z)
      REAL (KIND = GRID_SR), INTENT(IN)   :: x, y
      REAL (KIND = GRID_SR)               :: z

      REAL (KIND = GRID_SR), DIMENSION(3) :: r_q

      r_q = exact_solution([x, y], r_time)
      z = r_q(2)

    END FUNCTION FUNCmx

    FUNCTION FUNCmy(x,y) RESULT(z)
      REAL (KIND = GRID_SR), INTENT(IN)   :: x, y
      REAL (KIND = GRID_SR)               :: z

      REAL (KIND = GRID_SR), DIMENSION(3) :: r_q

      r_q = exact_solution([x, y], r_time)
      z = r_q(3)

    END FUNCTION FUNCmy

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
    INTEGER (KIND = GRID_SI)                                      :: i_it
    REAL (KIND = GRID_SR)                                         :: r_tol, r_x, &
      r_h_l, r_u_l, r_h_r, r_u_r, r_headl, r_taill, r_headr, r_tailr, r_hu_star, &
      r_a_l, r_a_r, r_Delta_u_crit, r_Delta_u, r_Delta_h, r_h_star, r_u_star, &
      r_a_star, r_f, r_df, r_xt, r_h0, r_u0, r_q_r, r_q_l, r_h_guess, S_r, S_l

    r_tol = 10E-6
    r_x   = r_coo(1) - p_testparam%treal(1)%p_value(1)
    r_h_l = p_testparam%treal(2)%p_value(1)
    r_u_l = p_testparam%treal(2)%p_value(2)
    r_h_r = p_testparam%treal(3)%p_value(1)
    r_u_r = p_testparam%treal(3)%p_value(2)

    IF (r_time < GRID_EPS) THEN
      IF (r_x < 0.0_GRID_SR) THEN
        r_h0 = r_h_l
        r_u0 = r_u_l
      ELSE
        r_h0 = r_h_r
        r_u0 = r_u_r
      END IF
    ELSE
      r_xt = r_x/r_time

!--- computational speed
      r_a_l = SQRT(GRID_GRAV*r_h_l)
      r_a_r = SQRT(GRID_GRAV*r_h_r)

!--- speed of the waves
      r_headl = r_u_l - r_a_l
      r_taill = r_u_l + 2.0_GRID_SR*r_a_l
      r_headr = r_u_r + r_a_r
      r_tailr = r_u_r - 2.0_GRID_SR*r_a_r

      IF ((r_h_l > 0.0) .AND. (r_h_r > 0.0)) THEN
        r_Delta_u      = r_u_r-r_u_l
        r_Delta_u_crit = 2.0_GRID_SR*(r_a_l+r_a_r)

        IF (r_Delta_u >= r_Delta_u_crit) THEN
        ! Data produce vacuum. Dry bed in the middle
          IF (r_xt < r_headl) THEN
            r_h0 = r_h_l
            r_u0 = r_u_l
          ELSEIF (r_xt < r_taill) THEN
            r_h0 = (r_u_l + 2.0_GRID_SR*r_a_l - r_xt)**2 / (9.0_GRID_SR*GRID_GRAV)
            r_u0 = (r_u_l + 2.0_GRID_SR*r_a_l + 2.0_GRID_SR*r_xt) / 3.0_GRID_SR
          ELSEIF (r_xt < r_tailr) THEN
            r_h0 = 0.0_GRID_SR
            r_u0 = 0.0_GRID_SR
          ELSEIF (r_xt < r_headr) THEN
            r_h0 = (-r_u_r + 2.0_GRID_SR*r_a_r + r_xt)**2 / (9.0_GRID_SR*GRID_GRAV)
            r_u0 = ( r_u_r - 2.0_GRID_SR*r_a_r + 2.0_GRID_SR*r_xt) / 3.0_GRID_SR
          ELSE
            r_h0 = r_h_r
            r_u0 = r_u_r
          END IF
        ELSE ! r_Delta_u >= r_Delta_u_crit
          ! Two rarefaction guess
          r_Delta_h = r_tol + 1.0_GRID_SR
          r_h_guess = swguess(r_h_l, r_u_l, r_a_l, r_h_r, r_u_r, r_a_r)
          i_it = 0
          DO WHILE ((r_Delta_h > r_tol) .AND. (i_it < 100))
            r_f  = swfunc1(r_h_guess,r_h_l) + swfunc1(r_h_guess,r_h_r) + r_Delta_u
            r_df = swfunc2(r_h_guess,r_h_l) + swfunc2(r_h_guess,r_h_r)
            r_h_star = r_h_guess - r_f/r_df
            r_Delta_h = 2.0_GRID_SR*ABS(r_h_star-r_h_guess)/(r_h_star+r_h_guess)
            i_it = i_it + 1
            r_h_guess = r_h_star
            r_u_star = 0.5_GRID_SR*(r_u_l+r_u_r) + &
                       0.5_GRID_SR*(swfunc1(r_h_star,r_h_r) - swfunc1(r_h_star,r_h_l))
            r_a_star = SQRT(GRID_GRAV*r_h_star)
          END DO

          r_a_star = SQRT(GRID_GRAV*r_h_star)
          IF (r_h_star <= r_h_l) THEN
            ! left rarefaction wave
            r_taill = r_u_star - r_a_star
          ELSE
            r_q_l = SQRT(0.5_GRID_SR*((r_h_star+r_h_l)*r_h_star/r_h_l**2))
            S_l   = r_u_l - r_a_l*r_q_l
          END IF

          IF (r_h_star <= r_h_r) THEN
            ! right rarefaction wave
            r_tailr = r_u_star + r_a_star
          ELSE
            r_q_r = SQRT(0.5_GRID_SR*((r_h_star+r_h_r)*r_h_star/r_h_r**2))
            S_r   = r_u_r + r_a_r*r_q_r
          END IF

          IF (r_xt < r_u_star) THEN
            IF (r_h_star <= r_h_l) THEN
              ! left rarefaction wave
              IF (r_xt < r_headl) THEN
                r_h0 = r_h_l
                r_u0 = r_u_l
              ELSEIF (r_xt <= r_taill) THEN
                r_h0 = (r_u_l + 2.0_GRID_SR*r_a_l - r_xt)**2 / (9.0_GRID_SR*GRID_GRAV)
                r_u0 = (r_u_l + 2.0_GRID_SR*r_a_l + 2.0_GRID_SR*r_xt) / 3.0_GRID_SR
              ELSE
                r_h0 = r_h_star
                r_u0 = r_u_star
              END IF
            ELSE
              ! left shock wave
              IF (r_xt < S_l) THEN
                r_h0 = r_h_l
                r_u0 = r_u_l
              ELSE
                r_h0 = r_h_star
                r_u0 = r_u_star
              END IF
            END IF
          ELSE
            IF (r_h_star <= r_h_r) THEN
              ! right rarefaction wave
              IF (r_xt > r_headr) THEN
                r_h0 = r_h_r
                r_u0 = r_u_r
              ELSEIF (r_xt >= r_tailr) THEN
                r_h0 = (-r_u_r + 2.0_GRID_SR*r_a_r + r_xt)**2 / (9.0_GRID_SR*GRID_GRAV)
                r_u0 = ( r_u_r - 2.0_GRID_SR*r_a_r + 2.0_GRID_SR*r_xt) / 3.0_GRID_SR
              ELSE
                r_h0 = r_h_star
                r_u0 = r_u_star
              END IF
            ELSE
              ! right shock wave
              IF (r_xt > S_r) THEN
                r_h0 = r_h_r
                r_u0 = r_u_r
              ELSE
                r_h0 = r_h_star
                r_u0 = r_u_star
              END IF
            END IF
          END IF
        END IF ! r_Delta_u >= r_Delta_u_crit
      ELSEIF ((r_h_l <= 0.0) .AND. (r_h_r > 0.0)) THEN
        ! dry bed region on the left
        IF (r_xt > r_headr) THEN
          r_h0 = r_h_r
          r_u0 = r_u_r
        ELSEIF (r_xt > r_tailr) THEN
          r_h0 = (-r_u_r + 2.0_GRID_SR*r_a_r + r_xt)**2 / (9.0_GRID_SR*GRID_GRAV)
          r_u0 = ( r_u_r - 2.0_GRID_SR*r_a_r + 2.0_GRID_SR*r_xt) / 3.0_GRID_SR
        ELSE
          r_h0 = 0.0_GRID_SR
          r_u0 = 0.0_GRID_SR
        END IF
      ELSEIF ((r_h_r <= 0.0) .AND. (r_h_l > 0.0)) THEN
        ! dry bed region on the right
        IF (r_xt < r_headl) THEN
          r_h0 = r_h_l
          r_u0 = r_u_l
        ELSEIF (r_xt < r_taill) THEN
          r_h0 = (r_u_l + 2.0_GRID_SR*r_a_l - r_xt)**2 / (9.0_GRID_SR*GRID_GRAV)
          r_u0 = (r_u_l + 2.0_GRID_SR*r_a_l + 2.0_GRID_SR*r_xt) / 3.0_GRID_SR
        ELSE
          r_h0 = 0.0_GRID_SR
          r_u0 = 0.0_GRID_SR
        END IF
      ELSE
        r_h0 = 0.0_GRID_SR
        r_u0 = 0.0_GRID_SR
      END IF
    END IF

    r_sol(1) = r_h0
    r_sol(2) = r_h0*r_u0
    r_sol(3) = 0.0_GRID_SR

  END FUNCTION exact_solution

!*******************************************************************************
! DESCRIPTION of [SUBROUTINE field_init]:
!> @brief initializes the height and momentum fields
!>
!> @param[in,out]   p_param       control structure for global parameters
!> @param           i_numelmt
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

    CALL dg_elmt_solution(r_coodof, 0.0_GRID_SR, i_numelmt, &
                          i_elmtdofs, r_coonod, i_elmtnodes, r_sol)
    r_Q = r_sol * GRID_GRAV

!--- deallocate workspace
    DEALLOCATE(r_sol)

  END SUBROUTINE field_init

!*******************************************************************************
! DESCRIPTION of [FUNCTION swguess]:
!> @brief
!
!> @param           r_hl
!> @param           r_ul
!> @param           r_al
!> @param           r_hr
!> @param           r_ur
!> @param           r_ar
!
  FUNCTION swguess(r_hl, r_ul, r_al, r_hr, r_ur, r_ar)

    IMPLICIT NONE

    REAL (KIND = GRID_SR)         :: r_hl, r_ul, r_al, r_hr, r_ur, r_ar
    REAL (KIND = GRID_SR)         :: swguess
    REAL (KIND = GRID_SR)         :: r_minh, r_h0, r_gr, r_gl, r_g

    r_g = GRID_GRAV

!--- minimum of h
    r_minh = MIN(r_hl,r_hr)

!--- two rarefaction guess
    swguess = (0.5*(r_al+r_ar) + 0.25*(r_ul-r_ur))**2 / r_g

    IF (swguess > r_minh) THEN
!--- use two shock guess
      r_h0 = swguess
      r_gl = SQRT(0.5*r_g * (r_h0+r_hl)/(r_h0*r_hl))
      r_gr = SQRT(0.5*r_g * (r_h0+r_hr)/(r_h0*r_hr))
      swguess = (r_gl*r_hl + r_gr*r_hr + r_ul - r_ur)/(r_gl + r_gr)
    END IF

  END FUNCTION swguess

!*******************************************************************************
! DESCRIPTION of [FUNCTION swfunc1]:
!> @brief
!
!> @param           r_h
!> @param           r_hk
!
  FUNCTION swfunc1(r_h, r_hk)

    IMPLICIT NONE

    REAL (KIND = GRID_SR)         :: r_h, r_hk
    REAL (KIND = GRID_SR)         :: r_a, r_ak, swfunc1

    r_a  = SQRT(GRID_GRAV*r_h)
    r_ak = SQRT(GRID_GRAV*r_hk)

    IF (r_h <= r_hk) THEN
      swfunc1 = 2.0*(r_a-r_ak)
    ELSE
      swfunc1 = (r_h-r_hk)*SQRT(0.5*GRID_GRAV*((r_h+r_hk)/(r_h*r_hk)))
    END IF

  END FUNCTION swfunc1

!*******************************************************************************
! DESCRIPTION of [FUNCTION swfunc2]:
!> @brief
!
!> @param           r_h
!> @param           r_hk
!
  FUNCTION swfunc2(r_h, r_hk)

    IMPLICIT NONE

    REAL (KIND = GRID_SR)         :: r_h, r_hk, r_a
    REAL (KIND = GRID_SR)         :: swfunc2, r_gk

    r_a = SQRT(GRID_GRAV*r_h)

    IF (r_h <= r_hk) THEN
      swfunc2 = GRID_GRAV/r_a
    ELSE
      r_gk    = SQRT(0.5*GRID_GRAV*(r_h+r_hk) / (r_h * r_hk))
      swfunc2 = r_gk - GRID_GRAV*(r_h-r_hk) / (4.0 * r_h**2 * r_gk)
    END IF

  END FUNCTION swfunc2

!*******************************************************************************
END MODULE DG_initial
