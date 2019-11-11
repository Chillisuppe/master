!*******************************************************************************
!
!> @file  DG_boundary_conicalisland.F90
!> @brief contains module DG_boundary
!
!*******************************************************************************
! DESCRIPTION:
!> @brief definition of boundary conditions on one edge
!
MODULE DG_boundary

  USE GRID_api
  USE FLASH_parameters
  USE DG_equation, ONLY : velocity
  USE IO_equation, ONLY : p_equationparam
  USE DG_initial, ONLY : boundaryhgt

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

    INTEGER (KIND = GRID_SI)                            :: i_point
    REAL (KIND = GRID_SR)                               :: r_w1, r_w2, r_h0, r_u0
    REAL (KIND = GRID_SR), DIMENSION(2)                 :: r_r1, r_r2
    REAL (KIND = GRID_SR)                               :: r_u_l, r_v_l, r_uv_norm, r_uv_tang, r_hu_norm

    SELECT CASE (i_edgebou)
    CASE(-3) ! reflecting boundary conditions
      r_Q_r(1,:) = r_Q_l(1,:)
      r_S_r      = r_S_l

      DO i_point = 1, i_npoints
        r_u_l = velocity(r_Q_l(1,i_point), r_Q_l(2,i_point))
        r_v_l = velocity(r_Q_l(1,i_point), r_Q_l(3,i_point))

        r_uv_norm =  r_normal(1)*r_u_l + r_normal(2)*r_v_l ! this is (\vec u,\vec n)
        r_uv_tang = -r_normal(2)*r_u_l + r_normal(1)*r_v_l ! this is (\vec u,\vec t)

        r_Q_r(2,i_point) = r_Q_r(1,i_point)*(-r_uv_tang*r_normal(2) - r_uv_norm*r_normal(1))
        r_Q_r(3,i_point) = r_Q_r(1,i_point)*( r_uv_tang*r_normal(1) - r_uv_norm*r_normal(2))
      END DO

    CASE(-4) ! transparent boundary conditions
!       r_Q_r = r_Q_l

      r_h0 = p_equationparam%r_depth
      r_u0 = 0.0_GRID_SR

      DO i_point = 1, i_npoints
        IF (r_Q_l(1,i_point) < p_equationparam%r_wettol) THEN
          r_Q_r(:,i_point) = r_Q_l(:,i_point)
        ELSE
          r_u_l = velocity(r_Q_l(1,i_point), r_Q_l(2,i_point))
          r_v_l = velocity(r_Q_l(1,i_point), r_Q_l(3,i_point))

          r_uv_norm = -r_normal(1)*r_u_l - r_normal(2)*r_v_l ! this is (\vec u,\vec n)
          r_uv_tang = -r_normal(2)*r_u_l + r_normal(1)*r_v_l ! this is (\vec u,\vec t)

          r_r1 = [1.0_GRID_SR, r_u0-SQRT(GRID_GRAV*r_h0-r_S_l(1,i_point))]
          r_r2 = [1.0_GRID_SR, r_u0+SQRT(GRID_GRAV*r_h0-r_S_l(1,i_point))]

          r_w1 = 1.0_GRID_SR/(2.0_GRID_SR*SQRT(GRID_GRAV*r_h0-r_S_l(1,i_point))) * &
            ((r_u0+SQRT(GRID_GRAV*r_h0-r_S_l(1,i_point)))*((r_Q_l(1,i_point)+r_S_l(1,i_point))/GRID_GRAV-r_h0) - &
            (r_Q_l(1,i_point)/GRID_GRAV*r_uv_norm-(r_h0-r_S_l(1,i_point)/GRID_GRAV)*r_u0))
          r_w2 = 0.0_GRID_SR

          r_Q_r(1,i_point) = GRID_GRAV * (r_h0-r_S_l(1,i_point)/GRID_GRAV + r_w1*r_r1(1) + r_w2*r_r2(1))
          r_hu_norm        = GRID_GRAV * ((r_h0-r_S_l(1,i_point)/GRID_GRAV)*r_u0 + r_w1*r_r1(2) + r_w2*r_r2(2))

        r_Q_r(2,i_point) = -r_Q_r(1,i_point)*r_uv_tang*r_normal(2) - r_hu_norm*r_normal(1)
        r_Q_r(3,i_point) =  r_Q_r(1,i_point)*r_uv_tang*r_normal(1) - r_hu_norm*r_normal(2)

        END IF
      END DO

      r_S_r = r_S_l

    CASE(-6) ! inflow/outflow boundary

      r_h0 = p_equationparam%r_depth
      r_u0 = 0.0_GRID_SR

      DO i_point = 1, i_npoints
        IF (r_Q_l(1,i_point) < p_equationparam%r_wettol) THEN
          r_Q_r(:,i_point) = r_Q_l(:,i_point)
        ELSE
          r_u_l = velocity(r_Q_l(1,i_point), r_Q_l(2,i_point))
          r_v_l = velocity(r_Q_l(1,i_point), r_Q_l(3,i_point))

          r_uv_norm = -r_normal(1)*r_u_l - r_normal(2)*r_v_l ! this is (\vec u,\vec n)
          r_uv_tang = -r_normal(2)*r_u_l + r_normal(1)*r_v_l ! this is (\vec u,\vec t)

          r_r1 = [1.0_GRID_SR, r_u0-SQRT(GRID_GRAV*r_h0-r_S_l(1,i_point))]
          r_r2 = [1.0_GRID_SR, r_u0+SQRT(GRID_GRAV*r_h0-r_S_l(1,i_point))]

          r_w1 = 1.0_GRID_SR/(2.0_GRID_SR*SQRT(GRID_GRAV*r_h0-r_S_l(1,i_point))) * &
            ((r_u0+SQRT(GRID_GRAV*r_h0-r_S_l(1,i_point)))*((r_Q_l(1,i_point)+r_S_l(1,i_point))/GRID_GRAV-r_h0) - &
            (r_Q_l(1,i_point)/GRID_GRAV*r_uv_norm-(r_h0-r_S_l(1,i_point)/GRID_GRAV)*r_u0))
          r_w2 = 1.0_GRID_SR/(2.0_GRID_SR*SQRT(GRID_GRAV*r_h0-r_S_l(1,i_point))) * &
            ((-r_u0+SQRT(GRID_GRAV*r_h0-r_S_l(1,i_point)))*boundaryhgt(r_time) + &
            (boundaryhgt(r_time)+r_h0-r_S_l(1,i_point)/GRID_GRAV) * &
            2.0_GRID_SR*(SQRT((boundaryhgt(r_time     )+r_h0)*GRID_GRAV-r_S_l(1,i_point)) - &
                         SQRT((boundaryhgt(0.0_GRID_SR)+r_h0)*GRID_GRAV-r_S_l(1,i_point))))

          r_Q_r(1,i_point) = GRID_GRAV * (r_h0-r_S_l(1,i_point)/GRID_GRAV + r_w1*r_r1(1) + r_w2*r_r2(1))
          r_hu_norm        = GRID_GRAV * ((r_h0-r_S_l(1,i_point)/GRID_GRAV)*r_u0 + r_w1*r_r1(2) + r_w2*r_r2(2))

          r_Q_r(2,i_point) = -r_Q_r(1,i_point)*r_uv_tang*r_normal(2) - r_hu_norm*r_normal(1)
          r_Q_r(3,i_point) =  r_Q_r(1,i_point)*r_uv_tang*r_normal(1) - r_hu_norm*r_normal(2)

        END IF
      END DO

      r_S_r = r_S_l

    CASE DEFAULT
      CALL grid_error(c_error='[compute_bc]: boundary condition not supported!')

    END SELECT

  END SUBROUTINE compute_bc

!*******************************************************************************
END MODULE DG_boundary
