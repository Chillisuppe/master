!*******************************************************************************
!
!> @file  DG_friction_constant.f90
!> @brief contains module DG_friction
!
!*******************************************************************************
! MODULE DESCRIPTION:
!> @brief contains routines to take care of constant friction
!
MODULE DG_friction

  USE GRID_api
  USE DG_utils

  PUBLIC  :: friction

  CONTAINS
!*******************************************************************************
! DESCRIPTION of [FUNCTION friction]:
!> @brief
!
!> @param[in]     r_h
!> @param[in]     r_hu
!> @param[in]     r_hv
!> @param[in]     r_gamma       bottom friction parameter
!
!> @note if Manning n negative in param file - adjust variable friction
!
  FUNCTION friction(r_h, r_hu, r_hv, r_gamma)

    IMPLICIT NONE

    REAL(KIND=GRID_SR), INTENT(in)                 :: r_h, r_hu, r_hv, r_gamma
    REAL(KIND=GRID_SR)                             :: friction
    REAL(KIND=GRID_SR)                             :: r_l2velo, r_root, r_tol, r_heavi

!-- Initialize tolerance for wet/dry
    r_tol = 10E-10

    r_heavi = heaviside(r_h)
!-- Compute constant Manning friction
    r_l2velo = 0._GRID_SR + (1-r_heavi)*(sqrt((r_hu/(r_h+r_tol))**2+(r_hv/(r_h+r_tol))**2))
    r_root   = r_h**(4/3)
    friction = 0._GRID_SR + (1-r_heavi)*(r_l2velo/(r_root+r_tol)*r_gamma**2*GRID_GRAV)

     IF (abs(r_root) .LT. 10E-6)THEN
       friction = 0._GRID_SR
     END IF

  END FUNCTION friction

!*******************************************************************************
END MODULE DG_friction

