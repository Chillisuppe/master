!*******************************************************************************
!
!> @file  DG_RS_upwind.F90
!> @brief contains module DG_riemann_solver
!
!*******************************************************************************
! MODULE DESCRIPTION:
!> @brief computes exact or approximate solution of Riemann problem
!
MODULE DG_riemann_solver

  USE GRID_api
  USE DG_equation

  PRIVATE
  PUBLIC  :: riemannsolver

  CONTAINS

!*******************************************************************************
! DESCRIPTION of [FUNCTION riemannsolver]:
!> @brief Upwind Riemann solver for linear advection equation
!>
!> @param[in]     r_Q_l           left state vector (phi)
!> @param[in]     r_Q_r           right state vector (phi)
!> @param[in]     r_normal        normal vector pointing from left to right element
!> @return                        solution of Riemann problem
!>
  FUNCTION riemannsolver(r_Q_l, r_Q_r, r_normal) RESULT(r_flux)

    IMPLICIT NONE

    REAL (KIND = GRID_SR), DIMENSION(i_nprogvars)                :: r_flux
    REAL (KIND = GRID_SR), DIMENSION(i_nprogvars), INTENT(IN)    :: r_Q_l, r_Q_r
    REAL (KIND = GRID_SR), DIMENSION(GRID_DIMENSION), INTENT(in) :: r_normal

    REAL (KIND = GRID_SR), DIMENSION(GRID_DIMENSION)             :: r_velocity
    REAL (KIND = GRID_SR), DIMENSION(i_nprogvars,GRID_DIMENSION) :: r_flux_aux

!--- compute velocities and their normal components
    r_velocity(1) = velocity(1)
    r_velocity(2) = velocity(2)

!--- compute left and right fluxes
    IF (DOT_PRODUCT(r_velocity, r_normal) >= 0.0_GRID_SR) THEN
      r_flux_aux = flux(r_Q_l)
    ELSE
      r_flux_aux = flux(r_Q_r)
    END IF

    r_flux = MATMUL(r_flux_aux, r_normal)

  END FUNCTION riemannsolver

END MODULE DG_riemann_solver
