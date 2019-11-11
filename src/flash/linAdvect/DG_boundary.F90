!*******************************************************************************
!
!> @file  DG_boundary.F90
!> @brief contains module DG_boundary
!
!*******************************************************************************
! DESCRIPTION:
!> @brief definition of boundary conditions on one edge
!
MODULE DG_boundary

  USE GRID_api
  USE FLASH_parameters

  PRIVATE
  PUBLIC :: compute_bc

  CONTAINS

!*******************************************************************************
! DESCRIPTION of [SUBROUTINE compute_bc]:
!> @brief compute boundary conditions on given edge (quadrature) points
!
!> @param[in]     r_time        current model time
!> @param[in]     r_coords      coordinates of boundary points
!> @param[in]     i_npoints     number of boundary points
!> @param[in]     i_edgebou     boundary flag of edge
!> @param[in]     r_normal      unit normal vector on boundary
!> @param[in]     r_Q_l         state vector at boundary points from left element (given)
!> @param[in]     r_S_l         discrete fields used in source terms from left element (given)
!> @param[out]    r_Q_r         state vector at boundary points from right element (to be computed)
!> @param[out]    r_S_r         discrete fields used in source terms from right element (to be computed)
!
  SUBROUTINE compute_bc(r_time, r_coords, i_npoints, i_edgebou, r_normal, &
                    r_Q_l, r_S_l, r_Q_r, r_S_r)

    IMPLICIT NONE

    REAL (KIND = GRID_SR),                 INTENT(IN)   :: r_time
    REAL (KIND = GRID_SR), DIMENSION(:,:), INTENT(IN)   :: r_coords
    INTEGER (KIND = GRID_SI),              INTENT(IN)   :: i_npoints, i_edgebou
    REAL (KIND = GRID_SR), DIMENSION(:),   INTENT(IN)   :: r_normal
    REAL (KIND = GRID_SR), DIMENSION(:,:), INTENT(IN)   :: r_Q_l, r_S_l
    REAL (KIND = GRID_SR), DIMENSION(:,:), INTENT(OUT)  :: r_Q_r, r_S_r

    SELECT CASE (i_edgebou)
    CASE(-4) ! extrapolation boundary conditions
      r_Q_r = r_Q_l
      r_S_r = r_S_l

    CASE DEFAULT
      CALL grid_error(c_error='[compute_bc]: boundary condition not supported!')

    END SELECT

  END SUBROUTINE compute_bc

!*******************************************************************************
END MODULE DG_boundary
