!*******************************************************************************
!
!> @file  DG_RS_rusanov.F90
!> @brief contains module DG_riemann_solver
!
!*******************************************************************************
! MODULE DESCRIPTION:
!> @brief computes exact or approximate solution of Riemann problem
!
MODULE DG_riemann_solver

  USE GRID_api
  USE DG_equation, ONLY : flux, velocity, i_nprogvars

  PRIVATE
  PUBLIC  :: riemannsolver, riemannsolverAdv

  CONTAINS

!*******************************************************************************
! DESCRIPTION of [FUNCTION riemannsolver]:
!> @brief Rusanov Riemann solver for shallow water equations
!>
!> @param[in]     r_Q_l           left state vector (fluid depth, momentum in x-, momentum in y-direction)
!> @param[in]     r_Q_r           right state vector (fluid depth, momentum in x-, momentum in y-direction)
!> @param[in]     r_normal        normal vector pointing from left to right element
!> @return                        solution of Riemann problem
!>
  FUNCTION riemannsolver(r_Q_l, r_Q_r, r_normal) RESULT(r_flux)

    IMPLICIT NONE

    REAL (KIND = GRID_SR), DIMENSION(i_nprogvars)                 :: r_flux
    REAL (KIND = GRID_SR), DIMENSION(i_nprogvars), INTENT(IN)     :: r_Q_l, r_Q_r
    REAL (KIND = GRID_SR), DIMENSION(GRID_DIMENSION), INTENT(IN)  :: r_normal

    REAL (KIND = GRID_SR)                                         :: r_u_l, r_u_r, &
      r_v_l, r_v_r, r_a_l, r_a_r, r_u_norm_l, r_u_norm_r, r_abs_l, r_abs_r, r_lambda
    REAL (KIND = GRID_SR), DIMENSION(i_nprogvars)                 :: r_dq
    REAL (KIND = GRID_SR), DIMENSION(i_nprogvars,GRID_DIMENSION)  :: r_flux_left, r_flux_right, &
                                                                     r_flux_help, r_fluxrus

!--- compute velocities and their normal components
    r_u_l = velocity(r_Q_l(1), r_Q_l(2))
    r_v_l = velocity(r_Q_l(1), r_Q_l(3))
    r_u_r = velocity(r_Q_r(1), r_Q_r(2))
    r_v_r = velocity(r_Q_r(1), r_Q_r(3))

    r_u_norm_l =  r_u_l*r_normal(1) + r_v_l*r_normal(2)
    r_u_norm_r =  r_u_r*r_normal(1) + r_v_r*r_normal(2)

!--- compute sound speeds
    r_a_l = sqrt(r_Q_l(1))
    r_a_r = sqrt(r_Q_r(1))

!--- compute maximal flux strength
    r_abs_l  = abs(r_u_norm_l) + r_a_l
    r_abs_r  = abs(r_u_norm_r) + r_a_r
    r_lambda = abs(max(r_abs_l, r_abs_r))

!--- compute left and right fluxes
    r_flux_left  = flux(r_Q_l)
    r_flux_right = flux(r_Q_r)

!--- compute q_r - q_l as three-dimensional vector
    r_dq = r_Q_r - r_Q_l

    r_flux_help = r_lambda*MATMUL(RESHAPE(r_dq, [3, 1]), RESHAPE(r_normal, [1, GRID_DIMENSION]))

    r_fluxrus = (r_flux_left + r_flux_right - r_flux_help)
    r_fluxrus = r_fluxrus/2.0_GRID_SR

    r_flux = MATMUL(r_fluxrus, r_normal)

  END FUNCTION riemannsolver

!*******************************************************************************
! DESCRIPTION of [FUNCTION riemannsolverAdv]:
!> @brief Rusanov Riemann Solver for Shallow Water Equations
!>
!> @param[in]     r_Q_l           left state vector (fluid depth, momentum in x-, momentum in y-direction)
!> @param[in]     r_Q_r           right state vector (fluid depth, momentum in x-, momentum in y-direction)
!> @param[in]     r_normal        normal vector pointing from left to right element
!> @return                        solution of Riemann problem
!
  FUNCTION riemannsolverAdv(r_Q_l, r_Q_r, r_normal) RESULT(r_flux)

    IMPLICIT NONE

    REAL (KIND = GRID_SR), DIMENSION(i_nprogvars)                 :: r_flux
    REAL (KIND = GRID_SR), DIMENSION(i_nprogvars), INTENT(IN)     :: r_Q_l, r_Q_r
    REAL (KIND = GRID_SR), DIMENSION(GRID_DIMENSION), INTENT(IN)  :: r_normal

    REAL (KIND = GRID_SR)                                         :: r_u_l, r_u_r, &
      r_v_l, r_v_r, r_u_norm_l, r_u_norm_r, r_lambda
    REAL (KIND = GRID_SR), DIMENSION(i_nprogvars)                 :: r_dq
    REAL (KIND = GRID_SR), DIMENSION(i_nprogvars,GRID_DIMENSION)  :: r_flux_left, r_flux_right, &
                                                                     r_flux_help, r_fluxrus

!--- compute velocities and their normal components
    r_u_l = velocity(r_Q_l(1), r_Q_l(2))
    r_v_l = velocity(r_Q_l(1), r_Q_l(3))
    r_u_r = velocity(r_Q_r(1), r_Q_r(2))
    r_v_r = velocity(r_Q_r(1), r_Q_r(3))

    r_u_norm_l =  r_u_l*r_normal(1) + r_v_l*r_normal(2)
    r_u_norm_r =  r_u_r*r_normal(1) + r_v_r*r_normal(2)

!--- compute maximal flux strength
    r_lambda = MAX(ABS(r_u_norm_l), ABS(r_u_norm_r))

!--- compute left and right fluxes
    r_flux_left  = flux(r_Q_l, r_iswet=0.0_GRID_SR)
    r_flux_right = flux(r_Q_r, r_iswet=0.0_GRID_SR)

!--- compute q_r - q_l as three-dimensional vector
    r_dq = r_Q_r - r_Q_l

    r_flux_help = r_lambda*MATMUL(RESHAPE(r_dq, [3, 1]), RESHAPE(r_normal, [1, GRID_DIMENSION]))

    r_fluxrus = (r_flux_left + r_flux_right - r_flux_help)
    r_fluxrus = r_fluxrus/2.0_GRID_SR

    r_flux = MATMUL(r_fluxrus, r_normal)

  END FUNCTION riemannsolverAdv

!*******************************************************************************
END MODULE DG_riemann_solver
