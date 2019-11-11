!*******************************************************************************
!
!> @file  DG_time_rosenbrock.F90
!> @brief contains module DG_timestepping
!
!*******************************************************************************
!
! VERSION(s):
! 1. version            michel baensch       02/2019
!
!*******************************************************************************
! MODULE DESCRIPTION:
!> @brief provides timestepping routine for DG method
!>
!> Implicit Rosenbrock-Wanner methods (combined with JFNK)
!
MODULE DG_timestepping

  USE GRID_api
  USE DG_limiter
  USE DG_flux
  USE DG_equation
  USE FLASH_parameters
  USE DG_time_utils

  PRIVATE
  PUBLIC   :: timestepping, read_coefficients, iotst_runtimeinfo

!-------------------------------------------------------------------------------
! DESCRIPTION of rt_tstinfo
!> structure for timestepping specific runtime information
  TYPE rt_tstinfo
    INTEGER (KIND = GRID_SI), DIMENSION(2)      :: i_niter
  END TYPE rt_tstinfo
  TYPE (rt_tstinfo)           :: p_rttstinfo

  CONTAINS

!*******************************************************************************
! DESCRIPTION of [SUBROUTINE timestepping]:
!> @brief Provides timestepping routine for DG method: - Implicit Rosenbrock-Wanner methods
!
!> @param[in]     p_ghand           grid handling data structure
!> @param         r_Q_curr          discrete solution vector from old time level
!> @param         r_Q_next          computed discrete solution vector from new time level
!> @param         r_S_curr          discrete fields used in source terms
!> @param         r_time            current model time
!> @param         r_dt              current time step length
!> @param         r_coeff           coefficient matrix of Rosenbrock time stepping
!> @param         i_numelmt         total number of elements
!> @param         i_numedge         total number of edges
!> @param         i_faceunknowns    number of DOFs per element
!> @param         i_degree          polynomial degree of element
!> @param         r_coonod          coordinates for each grid node
!> @param         r_coodof          coordinates for each dof in the grid
!> @param         r_metrics_inv     metric terms for each element; DIMENSION(2,2,i_numelmt)
!> @param         r_Dxi             derivative matrix in \f$\xi\f$ direction (from signature file)
!> @param         r_Deta            derivative matrix in \f$\eta\f$ direction (from signature file)
!> @param[in]     r_evander         Vandermonde matrix
!> @param[in]     r_einvvander      inverse Vandermonde matrix
!> @param         r_epsi            basis function evaluations at element quadrature points
!> @param         r_gpsi            basis function evaluations at edge quadrature points
!> @param         r_eMinvpsi
!> @param         r_eMinvdpsidxi
!> @param         r_eMinvdpsideta
!> @param         r_gMinvpsi
!> @param         i_equadpts        number of quadrature points per element
!> @param         i_gquadpts        number of quadrature points per edge
!> @param         r_eqwei           weights for the quadrature rule corresponding to the element quadrature points
!> @param         r_gqwei           weights for the quadrature rule corresponding to the edge quadrature points
!> @param         r_element_vol     volume of each element in grid
!> @param         r_edgelength      length of each edge in grid
!> @param         r_normals
!> @param         i_elementedges    edge indices of the edges of each element; DIMENSION(3,i_numelmt)
!> @param         i_elementnodes    node indices of each element; DIMENSION(3, i_numelmt)
!> @param         i_elementdofs     DOF indices corresponding to each element; DIMENSION(i_faceunknowns, i_numelmt)
!> @param         i_edgenodes       global node indices for each edge; DIMENSION(2,i_numedge)
!> @param         i_edgeinfo        edge-element relation for each edge: local edge index and global element index for adjacent elements; DIMENSION(4,i_numedge)
!> @param         i_edgeboundary    boundary flags for all edges
!
  SUBROUTINE timestepping(p_ghand, r_Q_curr, r_Q_next, r_S_curr, &
                          r_time, r_dt, r_coeff, i_numelmt, i_numedge, &
                          i_faceunknowns, i_degree, r_coonod, r_coodof, r_metrics_inv, &
                          r_Dxi, r_Deta, r_evander, r_einvvander, r_epsi, r_gpsi, &
                          r_eMinvpsi, r_eMinvdpsidxi, r_eMinvdpsideta, r_gMinvpsi, &
                          i_equadpts, i_gquadpts, r_eqwei, r_gqwei, r_elmtvolume, &
                          r_edgelength, r_normals, i_elmtedges, i_elmtnodes, &
                          i_elmtdofs, i_edgenodes, i_edgeinfo, i_edgeboundary)

    IMPLICIT NONE

    TYPE (grid_handle), INTENT(in)                          :: p_ghand
    REAL (KIND = GRID_SR), DIMENSION(:,:), INTENT(in)       :: r_Q_curr
    REAL (KIND = GRID_SR), DIMENSION(:,:), INTENT(out)      :: r_Q_next
    REAL (KIND = GRID_SR), DIMENSION(:,:,:), INTENT(in)     :: r_metrics_inv
    REAL (KIND = GRID_SR), DIMENSION(:,:), INTENT(in)       :: r_S_curr, r_coonod, r_coodof, &
                                                               r_normals, &
                                                               r_epsi, r_gpsi, r_eMinvpsi, &
                                                               r_gMinvpsi, r_eMinvdpsidxi, &
                                                               r_eMinvdpsideta, r_Dxi, r_Deta, &
                                                               r_evander, r_einvvander
    REAL (KIND = GRID_SR), DIMENSION(:), INTENT(in)         :: r_eqwei, r_gqwei, r_elmtvolume, &
                                                               r_edgelength
    REAL (KIND = GRID_SR), DIMENSION(:,:), INTENT(in)       :: r_coeff
    REAL (KIND = GRID_SR), INTENT(in)                       :: r_time, r_dt
    INTEGER (KIND = GRID_SI), INTENT(in)                    :: i_numedge, i_numelmt, &
                                                               i_faceunknowns, i_equadpts, i_gquadpts, &
                                                               i_degree
    INTEGER (KIND = GRID_SI), DIMENSION(:), INTENT(in)      :: i_edgeboundary
    INTEGER (KIND = GRID_SI), DIMENSION(:,:), INTENT(in)    :: i_elmtedges, i_elmtnodes, &
                                                               i_edgenodes, i_edgeinfo, i_elmtdofs

    INTEGER (KIND = GRID_SI)                                :: i_stages, i_stage, i_alct, i_cnt
    REAL (KIND = GRID_SR), DIMENSION(:,:), ALLOCATABLE      :: r_coeff_a, r_coeff_c, &
                                                               r_Q_laststage, r_Q_guess, r_flux, r_U
    REAL (KIND = GRID_SR), DIMENSION(:), ALLOCATABLE        :: r_coeff_m, &
                                                               r_rhs, r_lhs_guess, r_res
    REAL (KIND = GRID_SR)                                   :: r_coeff_gamma

!--- allocate workspace
    i_stages = SIZE(r_coeff,1)
    ALLOCATE(r_coeff_a(i_stages,i_stages), r_coeff_c(i_stages,i_stages), r_coeff_m(i_stages), &
             r_Q_laststage(i_nprogvars,i_numelmt*i_faceunknowns), &
             r_Q_guess(i_nprogvars,i_numelmt*i_faceunknowns), &
             r_flux(i_nprogvars,i_numelmt*i_faceunknowns), &
             r_U(i_nprogvars*i_numelmt*i_faceunknowns,i_stages), &
             r_rhs(i_nprogvars*i_numelmt*i_faceunknowns), &
             r_lhs_guess(i_nprogvars*i_numelmt*i_faceunknowns), &
             r_res(i_nprogvars*i_numelmt*i_faceunknowns), stat=i_alct)
    IF(i_alct /= 0) &
    CALL grid_error(c_error='[timestepping]: Could not allocate auxiliary arrays')

    r_coeff_a = MATMUL(r_coeff(:,2:2+i_stages), r_coeff(:,2+2*i_stages:1+3*i_stages))
    r_coeff_m = MATMUL(r_coeff(:,1), r_coeff(:,2+2*i_stages:1+3*i_stages))
    r_coeff_c = -r_coeff(:,2+2*i_stages:1+3*i_stages)
    r_coeff_gamma = r_coeff(1,2+i_stages)
    DO i_cnt = 1, i_stages
      r_coeff_c(i_cnt,i_cnt) = r_coeff_c(i_cnt,i_cnt) + 1.0_GRID_SR / r_coeff_gamma
    END DO

!--- initialize arrays for loop
    r_Q_laststage = r_Q_curr
    r_Q_guess     = r_Q_curr
    r_U           = 0.0_GRID_SR

!--- compute flux matrix
    CALL fvm_flux(r_flux, r_Q_laststage, r_S_curr, r_time, i_numelmt, i_numedge, &
                  i_faceunknowns, r_metrics_inv, r_Dxi, r_Deta, r_epsi, r_gpsi, &
                  r_eMinvpsi, r_eMinvdpsidxi, r_eMinvdpsideta, r_gMinvpsi, &
                  i_equadpts, i_gquadpts, r_eqwei, r_gqwei, r_elmtvolume, &
                  r_edgelength, r_normals, i_edgeinfo, i_edgeboundary, &
                  i_edgenodes, i_elmtdofs, r_coonod, r_coodof)

!     CALL preconditioner(r_flux, r_Q_curr, i_numelmt, i_faceunknowns, r_metrics_inv, &
!                         r_Dxi, r_Deta, i_elmtdofs)

!--- stage loop
    stage_loop: DO i_stage=1, i_stages
      r_flux = 0.0_GRID_SR

!--- calculate RHS
      DO i_cnt = 1, i_stage-1
        r_Q_laststage = r_Q_curr + r_coeff_a(i_stage,i_cnt) *  &
                                   RESHAPE(r_U(:,i_cnt), SHAPE(r_Q_curr))
      END DO

!--- compute flux matrix
      CALL fvm_flux(r_flux, r_Q_laststage, r_S_curr, r_time, i_numelmt, i_numedge, &
                    i_faceunknowns, r_metrics_inv, r_Dxi, r_Deta, r_epsi, r_gpsi, &
                    r_eMinvpsi, r_eMinvdpsidxi, r_eMinvdpsideta, r_gMinvpsi, &
                    i_equadpts, i_gquadpts, r_eqwei, r_gqwei, r_elmtvolume, &
                    r_edgelength, r_normals, i_edgeinfo, i_edgeboundary, &
                    i_edgenodes, i_elmtdofs, r_coonod, r_coodof)

      r_rhs = RESHAPE(r_flux, SHAPE(r_rhs)) * r_coeff_gamma * r_dt

      DO i_cnt = 1, i_stage-1
        r_rhs = r_rhs + r_coeff_c(i_stage,i_cnt) * r_U(:,i_cnt) * r_coeff_gamma
      END DO

      r_lhs_guess = JFNK_LHS(RESHAPE(r_Q_guess, SHAPE(r_rhs)), r_dt, r_coeff_gamma, &
                             r_flux, r_Q_curr, r_S_curr, r_time, i_numelmt, &
                             i_numedge, i_faceunknowns, &
                             r_metrics_inv, r_Dxi, r_Deta, r_epsi, r_gpsi, &
                             r_eMinvpsi, r_eMinvdpsidxi, r_eMinvdpsideta, r_gMinvpsi, &
                             i_equadpts, i_gquadpts, r_eqwei, r_gqwei, r_elmtvolume, &
                             r_edgelength, r_normals, i_edgeinfo, i_edgeboundary, &
                             i_edgenodes, i_elmtdofs, r_coonod, r_coodof)

!--- compute residual and then solve Au = RHS
      r_res = r_rhs - r_lhs_guess

      IF (NORM2(r_res) < GRID_EPS) EXIT

      CALL GMRES(r_U(:,i_stage), RESHAPE(r_Q_guess, SHAPE(r_rhs)), r_res, &
                 r_dt, r_coeff_gamma, &
                 r_flux, r_Q_curr, r_S_curr, r_time, i_numelmt, i_numedge, &
                 i_faceunknowns, r_metrics_inv, r_Dxi, r_Deta, r_epsi, r_gpsi, &
                 r_eMinvpsi, r_eMinvdpsidxi, r_eMinvdpsideta, r_gMinvpsi, &
                 i_equadpts, i_gquadpts, r_eqwei, r_gqwei, r_elmtvolume, &
                 r_edgelength, r_normals, i_edgeinfo, i_edgeboundary, &
                 i_edgenodes, i_elmtdofs, r_coonod, r_coodof, p_rttstinfo%i_niter(i_stage))!, &
!                  r_tolerance=1E-4_GRID_SR, i_maxiterations=200)

      r_Q_guess = r_Q_guess + r_coeff_m(i_stage) * RESHAPE(r_U(:,i_stage),SHAPE(r_Q_guess))

    END DO stage_loop

!--- compute next stage and apply limiter to it
    r_Q_next = r_Q_guess

    CALL limiter(p_ghand, r_Q_next, i_elmtedges, i_edgeinfo, &
                 r_edgelength, r_metrics_inv, i_numelmt, r_coonod, r_coodof, &
                 i_edgenodes, i_elmtnodes, i_elmtdofs, i_numedge, &
                 r_S_curr, i_faceunknowns, r_evander, r_einvvander)

!--- deallocate workspace
    DEALLOCATE (r_Q_laststage, r_flux, r_U, r_rhs, r_lhs_guess, r_res)

  END SUBROUTINE timestepping

!*******************************************************************************
! DESCRIPTION of [SUBROUTINE read_coefficients]:
!> @brief The coefficient structure is as follows:
!> @      r_coeff = beta, Alpha, Gamma. Gamma^(-1) so the coefficients are comprised of
!> @                sx1,  sxs,   sxs,   sxs arrays
!
!> @param[in]     p_param           control structure for global parameters
!> @param[out]    r_coeff
!
  SUBROUTINE read_coefficients(p_param, r_coeff)

    IMPLICIT NONE

    TYPE(control_struct), INTENT(in)                                :: p_param
    REAL (KIND = GRID_SR), DIMENSION(:,:), ALLOCATABLE, INTENT(out) :: r_coeff

    INTEGER                                                         :: i_alct

    SELECT CASE(p_param%num%c_timestepping)
      CASE('ros2')
        ! 2nd order Rosenbrock method
        ALLOCATE(r_coeff(2, 7), stat=i_alct)
        IF(i_alct /= 0) &
          CALL grid_error(c_error='[read_coefficients]: Could not allocate enough memory')
        r_coeff(1,:) = [0.5_GRID_SR, 0.0_GRID_SR, 0.0_GRID_SR,  0.292893218813453_GRID_SR, &
          0.000000000000000_GRID_SR, 3.414213562373095_GRID_SR, 0.000000000000000_GRID_SR]
        r_coeff(2,:) = [0.5_GRID_SR, 1.0_GRID_SR, 0.0_GRID_SR, -0.585786437626905_GRID_SR, &
          0.292893218813453_GRID_SR, 6.828427124746189_GRID_SR, 3.414213562373095_GRID_SR]

      CASE default
        CALL grid_error(c_error='[read_coefficients]: timestepping method not supported!')
    END SELECT

  END SUBROUTINE read_coefficients

!*******************************************************************************
! DESCRIPTION of [SUBROUTINE iotst_runtimeinfo]:
!> @brief prints some information on the current timestepping specific run time status
!>
!> @param[in]   i_iounit      IO unit where output is sent
!
  SUBROUTINE iotst_runtimeinfo(i_iounit)

    IMPLICIT NONE

    INTEGER (KIND = GRID_SI), INTENT(in)                      :: i_iounit

!--- local declarations

!--- output to i_iounit
    WRITE(i_iounit,'(A)')         ' +++++ ----- ----- ----- ----- ----- ----- ----- ----- +++++'
    WRITE(i_iounit,'(A,i8,A)') ' +++++ no. iterations (Rosenbrock), 1. stage  ', p_rttstinfo%i_niter(1), ' +++++'
    WRITE(i_iounit,'(A,i8,A)') ' +++++ no. iterations (Rosenbrock), 2. stage  ', p_rttstinfo%i_niter(2), ' +++++'

  END SUBROUTINE iotst_runtimeinfo

!*******************************************************************************
END MODULE DG_timestepping
