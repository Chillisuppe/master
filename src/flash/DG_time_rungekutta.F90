!*******************************************************************************
!
!> @file  DG_time_rungekutta.F90
!> @brief contains module DG_timestepping
!
!*******************************************************************************
!
! VERSION(s):
! 1. version            nicole beisiegel        12/2011
!
!*******************************************************************************
! MODULE DESCRIPTION:
!> @brief provides timestepping routine for DG method
!>
!> explicit SSP Runge Kutta methods
!
MODULE DG_timestepping

  USE GRID_api
  USE DG_limiter
  USE DG_flux
  USE DG_equation
  USE FLASH_parameters

  PRIVATE
  PUBLIC   :: timestepping, read_coefficients, iotst_runtimeinfo

  CONTAINS

!*******************************************************************************
! DESCRIPTION of [SUBROUTINE timestepping]:
!> @brief Provides timestepping routine for DG method: - Explicit SSP Runge Kutta methods
!
!> @param[in]     p_ghand           grid handling data structure
!> @param         r_Q_curr          discrete solution vector from old time level
!> @param         r_Q_next          computed discrete solution vector from new time level
!> @param         r_S_curr          discrete fields used in source terms
!> @param         r_time            current model time
!> @param         r_dt              current time step length
!> @param         r_coeff           coefficient matrix of Runge-Kutta time stepping
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
!> @param         r_elmtvolume      volume of each element in grid
!> @param         r_edgelength      length of each edge in grid
!> @param         r_normals
!> @param         i_elmtedges       edge indices of the edges of each element; DIMENSION(3,i_numelmt)
!> @param         i_elmtnodes       node indices of each element; DIMENSION(3, i_numelmt)
!> @param         i_elmtdofs        DOF indices corresponding to each element; DIMENSION(i_faceunknowns, i_numelmt)
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

    INTEGER (KIND = GRID_SI)                                :: i_stage, i_alct
    REAL (KIND = GRID_SR), DIMENSION(:,:), ALLOCATABLE      :: r_Q_secondstage, r_Q_laststage, r_flux

!--- allocate workspace
    ALLOCATE(r_Q_secondstage(i_nprogvars,i_numelmt*i_faceunknowns), &
             r_Q_laststage(i_nprogvars,i_numelmt*i_faceunknowns), &
             r_flux(i_nprogvars,i_numelmt*i_faceunknowns), stat=i_alct)
          IF(i_alct /= 0) &
            CALL grid_error(c_error='[timestepping]: Could not allocate auxiliary arrays')

!--- initialize arrays for Runge-Kutta loop
    r_Q_laststage   = r_Q_curr
    r_Q_secondstage = 0.0_GRID_SR

!--- Runge-Kutta loop
    runge_kutta_loop: DO i_stage=1, SIZE(r_coeff,1)
      r_flux = 0.0_GRID_SR

!--- compute flux matrix
      CALL fvm_flux(r_flux, r_Q_laststage, r_S_curr, r_time, i_numelmt, i_numedge, &
                    i_faceunknowns, r_metrics_inv, r_Dxi, r_Deta, r_epsi, r_gpsi, &
                    r_eMinvpsi, r_eMinvdpsidxi, r_eMinvdpsideta, r_gMinvpsi, &
                    i_equadpts, i_gquadpts, r_eqwei, r_gqwei, r_elmtvolume, &
                    r_edgelength, r_normals, i_edgeinfo, i_edgeboundary, &
                    i_edgenodes, i_elmtdofs, r_coonod, r_coodof)

!--- compute next stage and apply limiter to it
      r_Q_next =        r_coeff(i_stage,1) * r_Q_curr &
                 +      r_coeff(i_stage,2) * r_Q_laststage &
                 +      r_coeff(i_stage,3) * r_Q_secondstage &
                 + r_dt*r_coeff(i_stage,4) * r_flux

      CALL limiter(p_ghand, r_Q_next, i_elmtedges, i_edgeinfo, &
                   r_edgelength, r_metrics_inv, i_numelmt, r_coonod, r_coodof, &
                   i_edgenodes, i_elmtnodes, i_elmtdofs, i_numedge, &
                   r_S_curr, i_faceunknowns, r_evander, r_einvvander)

!---  update for next stage
      r_Q_laststage = r_Q_next

!--- store the 2nd stage for SSP(5,3)
      IF (i_stage == 2) r_Q_secondstage = r_Q_next

    END DO runge_kutta_loop

!--- deallocate workspace
    DEALLOCATE (r_Q_secondstage, r_Q_laststage, r_flux)

  END SUBROUTINE timestepping

!*******************************************************************************
! DESCRIPTION of [SUBROUTINE read_coefficients]:
!> @brief returns coefficient matrix for SSP Runge-Kutta schemes (SSP(k,s), where
!>        k is the order and s the number of stages). Implemented are: euler, rk22
!>        (Heun), rk33, rk34 and rk35
!>
!> The coefficient structure is as follows:
!>   1. column: /f$\alpha_{i,1}/f$
!>   2. column: /f$\alpha_{i,i}/f$, where the first entry is set to zero
!>   3. column: /f$\alpha_{i,2}/f$
!>   4. column: /f$\beta_{i,i}/f$
!>
!> See \cite Kubatko2007, \cite Ruuth2006, \cite Spiteri2002
!
!> @param[in]     p_param           control structure for global parameters
!> @param[out]    r_coeff
!
  SUBROUTINE read_coefficients(p_param, r_coeff)

    IMPLICIT NONE

    TYPE(control_struct), INTENT(in)                                :: p_param
    REAL (KIND = GRID_SR), DIMENSION(:,:), ALLOCATABLE, INTENT(out) :: r_coeff

    INTEGER                                                         :: i_alct, i_maxlen

!--- set constants
    i_maxlen = 4

    SELECT CASE(p_param%num%c_timestepping)
      CASE('euler')
        ! explicit Euler
        ALLOCATE(r_coeff(1, i_maxlen), stat=i_alct)
        IF(i_alct /= 0) &
          CALL grid_error(c_error='[read_coefficients]: Could not allocate enough memory')
        r_coeff(1,:) = [1.0_GRID_SR, 0.0_GRID_SR, 0.0_GRID_SR, 1.0_GRID_SR]
      CASE('rk22')
        ! SSPRK(2,2), Heun method
        ALLOCATE(r_coeff(2, i_maxlen), stat=i_alct)
        IF(i_alct /= 0) &
          CALL grid_error(c_error='[read_coefficients]: Could not allocate enough memory')
        r_coeff(1,:) = [1.0_GRID_SR, 0.0_GRID_SR, 0.0_GRID_SR, 1.0_GRID_SR]
        r_coeff(2,:) = [0.5_GRID_SR, 0.5_GRID_SR, 0.0_GRID_SR, 0.5_GRID_SR]
      CASE('rk33')
        ! SSPRK(3,3), third-order Runge-Kutta
        ALLOCATE(r_coeff(3, i_maxlen), stat=i_alct)
        IF(i_alct /= 0) &
          CALL grid_error(c_error='[read_coefficients]: Could not allocate enough memory')
        r_coeff(1,:) = [1.0000000000_GRID_SR, 0.0000000000_GRID_SR, 0.0_GRID_SR, 1.0000000000_GRID_SR]
        r_coeff(2,:) = [0.7500000000_GRID_SR, 0.2500000000_GRID_SR, 0.0_GRID_SR, 0.2500000000_GRID_SR]
        r_coeff(3,:) = [0.3333333333_GRID_SR, 0.6666666666_GRID_SR, 0.0_GRID_SR, 0.6666666666_GRID_SR]
      CASE('rk34')
        ! SSPRK(3,4)
        ALLOCATE(r_coeff(4, i_maxlen), stat=i_alct)
        IF(i_alct /= 0) &
          CALL grid_error(c_error='[read_coefficients]: Could not allocate enough memory')
        r_coeff(1,:) = [1.0000000000_GRID_SR, 0.0000000000_GRID_SR, 0.0_GRID_SR, 0.5000000000_GRID_SR]
        r_coeff(2,:) = [0.0000000000_GRID_SR, 1.0000000000_GRID_SR, 0.0_GRID_SR, 0.5000000000_GRID_SR]
        r_coeff(3,:) = [0.6666666666_GRID_SR, 0.3333333333_GRID_SR, 0.0_GRID_SR, 0.1666666666_GRID_SR]
        r_coeff(4,:) = [0.0000000000_GRID_SR, 1.0000000000_GRID_SR, 0.0_GRID_SR, 0.5000000000_GRID_SR]
      CASE('rk35')
        ! SSPRK(3,5)
        ALLOCATE(r_coeff(5, i_maxlen), stat=i_alct)
        IF(i_alct /= 0) &
          CALL grid_error(c_error='[read_coefficients]: Could not allocate enough memory')
        r_coeff(1,:) = [1.000000000000000_GRID_SR, 0.000000000000000_GRID_SR, 0.000000000000000_GRID_SR, 0.377268915331368_GRID_SR]
        r_coeff(2,:) = [0.000000000000000_GRID_SR, 1.000000000000000_GRID_SR, 0.000000000000000_GRID_SR, 0.377268915331368_GRID_SR]
        r_coeff(3,:) = [0.355909775063327_GRID_SR, 0.644090224936674_GRID_SR, 0.000000000000000_GRID_SR, 0.242995220537396_GRID_SR]
        r_coeff(4,:) = [0.367933791638137_GRID_SR, 0.632066208361863_GRID_SR, 0.000000000000000_GRID_SR, 0.238458932846290_GRID_SR]
        r_coeff(5,:) = [0.000000000000000_GRID_SR, 0.762406163401431_GRID_SR, 0.237593836598569_GRID_SR, 0.287632146308408_GRID_SR]
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
!     WRITE(i_iounit,'(A)')         ' +++++ ----- ----- ----- ----- ----- ----- ----- ----- +++++'

  END SUBROUTINE iotst_runtimeinfo

!*******************************************************************************
END MODULE DG_timestepping
