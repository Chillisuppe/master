!*******************************************************************************
!
!> @file  DG_boundary_atmatrest.F90
!> @brief contains module DG_boundary
!
!*******************************************************************************
! DESCRIPTION:
!> @brief definition of boundary conditions on one edge
!
MODULE DG_boundary

  USE GRID_api
  USE FLASH_parameters
  USE DG_equation, ONLY : velocity, PrimitiveToProgvars
  USE IO_equation, ONLY : p_equationparam
  USE DG_initial,  ONLY : p_testparam

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

    INTEGER (KIND = GRID_SI)                            :: i_point, i_alct
    REAL (KIND = GRID_SR)                               :: r_u_l, r_v_l, r_uv_norm, &
      r_uv_tang, r_theta, r_rhoref, r_gammaaux, r_pressure, r_rho

    SELECT CASE (i_edgebou)
    CASE(-1) ! Dirichlet boundary
      r_theta    = p_testparam%treal(1)%p_value(1)
      r_rhoref   = p_equationparam%r_pref/(p_equationparam%r_rgas*r_theta)
      r_gammaaux = (p_equationparam%r_gamma - 1.0_GRID_SR) / p_equationparam%r_gamma

      DO i_point = 1, i_npoints
        r_pressure = p_equationparam%r_pref * (1.0_GRID_SR - r_coords(2,i_point) * r_gammaaux * GRID_GRAV * &
                     r_rhoref / p_equationparam%r_pref)**(1.0_GRID_SR/r_gammaaux)
        r_rho      = r_rhoref * (r_pressure / p_equationparam%r_pref)**(1.0_GRID_SR/p_equationparam%r_gamma)

        r_u_l = velocity(r_Q_l(1,i_point), r_Q_l(2,i_point))
        r_v_l = velocity(r_Q_l(1,i_point), r_Q_l(3,i_point))

        CALL PrimitiveToProgvars([r_rho, p_testparam%treal(2)%p_value, r_pressure], r_Q_r(:,i_point))
        r_S_r = 0.0_GRID_SR
      END DO

    CASE(-3) ! reflecting boundary conditions
      r_Q_r(1,:) = r_Q_l(1,:)
      r_Q_r(4,:) = r_Q_l(4,:)
      r_S_r      = r_S_l

      DO i_point = 1, i_npoints
        r_u_l = velocity(r_Q_l(1,i_point), r_Q_l(2,i_point))
        r_v_l = velocity(r_Q_l(1,i_point), r_Q_l(3,i_point))

        r_uv_norm =  r_normal(1)*r_u_l + r_normal(2)*r_v_l ! this is (\vec u,\vec n)
        r_uv_tang = -r_normal(2)*r_u_l + r_normal(1)*r_v_l ! this is (\vec u,\vec t)

        r_Q_r(2,i_point) = r_Q_r(1,i_point)*(-r_uv_tang*r_normal(2) - r_uv_norm*r_normal(1))
        r_Q_r(3,i_point) = r_Q_r(1,i_point)*( r_uv_tang*r_normal(1) - r_uv_norm*r_normal(2))
      END DO

    CASE(-4) ! extrapolation boundary conditions
      r_Q_r = r_Q_l
      r_S_r = r_S_l

    CASE(-5) ! no-flux boundary conditions
      r_Q_r(1,:) = r_Q_l(1,:)
      r_Q_r(4,:) = r_Q_l(4,:)
      r_S_r      = r_S_l

      DO i_point = 1, i_npoints
        r_u_l = velocity(r_Q_l(1,i_point), r_Q_l(2,i_point))
        r_v_l = velocity(r_Q_l(1,i_point), r_Q_l(3,i_point))

        r_uv_norm =  r_normal(1)*r_u_l + r_normal(2)*r_v_l ! this is (\vec u,\vec n)
        r_uv_tang = -r_normal(2)*r_u_l + r_normal(1)*r_v_l ! this is (\vec u,\vec t)

        r_Q_r(2,i_point) = -r_Q_r(1,i_point)*r_uv_tang*r_normal(2)
        r_Q_r(3,i_point) =  r_Q_r(1,i_point)*r_uv_tang*r_normal(1)
      END DO

    CASE DEFAULT
      CALL grid_error(c_error='[compute_bc]: boundary condition not supported!')

    END SELECT

  END SUBROUTINE compute_bc

!*******************************************************************************
END MODULE DG_boundary
