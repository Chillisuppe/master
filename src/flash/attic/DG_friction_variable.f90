!*******************************************************************************
!
!> @file  DG_friction_variable.f90
!> @brief contains module DG_friction
!
!*******************************************************************************
! MODULE DESCRIPTION:
!> @brief contains routines to take care of variable friction
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

!-- Local declarations

    INTEGER(KIND=GRID_SI)                          :: r_num_areas
    REAL(KIND=GRID_SR)                             :: r_n
    REAL(KIND=GRID_SR)                             :: r_max_fric, r_min_fric, &
                                                      r_max_height, r_min_height

    REAL(KIND=GRID_SR)                             :: r_l2velo, r_root, r_tol, r_heavi

    r_tol = 10E-10
    friction = 0._GRID_SR

    r_heavi = heaviside(r_h)

    r_num_areas = 2

    r_min_height = 0._GRID_SR
    r_max_height = 100.0_GRID_SR

    r_min_fric  = 0.02 !(clean earth channel)
    r_max_fric  = 0.05 !(0.075 or 0.15 if there is heavy brush or trees)

    r_n = r_h*(r_min_fric - r_max_fric)/(r_max_height)+ r_max_fric

!-- Compute Manning friction
    r_l2velo = 0._GRID_SR + (1-r_heavi)*(sqrt((r_hu/(r_h+r_tol))**2+(r_hv/(r_h+r_tol))**2))
    r_root   = r_h**(4/3)
    friction = 0._GRID_SR + (1-r_heavi)*(r_l2velo/(r_root+r_tol)*r_n**2*GRID_GRAV)

  END FUNCTION friction

!*******************************************************************************
END MODULE DG_friction

