!*******************************************************************************
!
!> @file  DG_RS_hllc.F90
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
  PUBLIC  :: riemannsolver

  CONTAINS

!*******************************************************************************
! DESCRIPTION of [FUNCTION riemannsolver]:
!> @brief HLLC Riemann solver for shallow water equations
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

    REAL (KIND = GRID_SR)                                         :: r_u_norm_l, r_u_norm_r, &
                                                                     r_u_tang_r, r_u_tang_l, &
                                                                     r_a_l, r_a_r, r_phi_star, &
                                                                     r_qr, r_ql, &
                                                                     r_S_l, r_S_r, r_tol, &
                                                                     r_S_star
    REAL (KIND = GRID_SR), DIMENSION(i_nprogvars)                 :: r_Q_rot_l, r_Q_rot_r, &
                                                                     r_flux_hll, r_Q_star
    REAL (KIND = GRID_SR), DIMENSION(i_nprogvars,GRID_DIMENSION)  :: r_flux_l, r_flux_r

    r_tol = 10E-10

    IF ((r_Q_l(1)>r_tol) .OR. (r_Q_r(1)>r_tol)) THEN
!--- rotate state vectors
      r_Q_rot_l(1) =  r_Q_l(1)
      r_Q_rot_l(2) =  r_Q_l(2)*r_normal(1) + r_Q_l(3)*r_normal(2)
      r_Q_rot_l(3) = -r_Q_l(2)*r_normal(2) + r_Q_l(3)*r_normal(1)

      r_Q_rot_r(1) =  r_Q_r(1)
      r_Q_rot_r(2) =  r_Q_r(2)*r_normal(1) + r_Q_r(3)*r_normal(2)
      r_Q_rot_r(3) = -r_Q_r(2)*r_normal(2) + r_Q_r(3)*r_normal(1)

!--- compute velocity
      r_u_norm_l = velocity(r_Q_rot_l(1), r_Q_rot_l(2))
      r_u_norm_r = velocity(r_Q_rot_r(1), r_Q_rot_r(2))

      r_u_tang_l = velocity(r_Q_rot_l(1), r_Q_rot_l(3))
      r_u_tang_r = velocity(r_Q_rot_r(1), r_Q_rot_r(3))

!--- compute sound speeds
      r_a_l = SQRT(r_Q_rot_l(1))
      r_a_r = SQRT(r_Q_rot_r(1))

!--- compute approximation for $\phi$ in the star region as in Toro (2001), equation (10.17)
!       r_phi_star = 0.5_GRID_SR*(r_Q_rot_r(1)+r_Q_rot_l(1)) - &
!                    0.25_GRID_SR*(r_u_norm_r-r_u_norm_l)*(r_Q_rot_r(1)+r_Q_rot_l(1))/(r_a_l+r_a_r)

!--- compute approximation for $\phi$ in the star region as in Toro (2001), equation (10.18)
      r_phi_star = (0.5_GRID_SR*(r_a_l + r_a_r) + 0.25_GRID_SR*(r_u_norm_l - r_u_norm_r))**2

!--- compute qL/qR according to equation (10.23)
      IF (r_phi_star > r_Q_rot_l(1)) THEN
        r_ql = SQRT(0.5_GRID_SR * r_phi_star*(r_phi_star+r_Q_rot_l(1)) / (r_Q_rot_l(1)**2))
      ELSE
        r_ql = 1.0_GRID_SR
      END IF

      IF (r_phi_star > r_Q_rot_r(1)) THEN
        r_qr = SQRT(0.5_GRID_SR * r_phi_star*(r_phi_star+r_Q_rot_r(1)) / (r_Q_rot_r(1)**2))
      ELSE
        r_qr = 1.0_GRID_SR
      END IF

!--- compute wave speeds of fastest and slowest wave using the Fraccarollo and Toro fix for dry beds
      IF (r_Q_rot_l(1) < r_tol) THEN
        r_S_l = r_u_norm_r - 2.0_GRID_SR*r_a_r
      ELSE
        r_S_l = r_u_norm_l - r_a_l*r_ql
      END IF

      IF (r_Q_rot_r(1) < r_tol) THEN
        r_S_r = r_u_norm_l + 2.0_GRID_SR*r_a_l
      ELSE
        r_S_r = r_u_norm_r + r_a_r*r_qr
      END IF

!--- compute S_star according to equation (10.27)
      r_S_star = (r_S_l*r_Q_rot_r(1)*(r_u_norm_r - r_S_r) - &
                  r_S_r*r_Q_rot_l(1)*(r_u_norm_l - r_S_l)) / &
                 (r_Q_rot_r(1)*(r_u_norm_r - r_S_r) - &
                  r_Q_rot_l(1)*(r_u_norm_l - r_S_l))

!--- compute left and right fluxes
      r_flux_l = flux(r_Q_rot_l)
      r_flux_r = flux(r_Q_rot_r)

!--- compute HLLC flux
      IF (r_S_l >= 0.0_GRID_SR) THEN
        r_flux_hll = r_flux_l(:,1)
      ELSEIF ((r_S_l < 0.0_GRID_SR) .AND. (0.0_GRID_SR <= r_S_star)) THEN
        r_Q_star(1) = r_Q_rot_l(1)*(r_S_l-r_u_norm_l) / (r_S_l-r_S_star)
        r_Q_star(2) = r_Q_rot_l(1)*(r_S_l-r_u_norm_l) / (r_S_l-r_S_star)*r_S_star
        r_Q_star(3) = r_Q_rot_l(1)*(r_S_l-r_u_norm_l) / (r_S_l-r_S_star)*r_u_tang_l
        r_flux_hll  = r_flux_l(:,1) + r_S_l*(r_Q_star-r_Q_rot_l)
      ELSEIF ((r_S_star < 0.0_GRID_SR) .AND. (0.0_GRID_SR <= r_S_r)) THEN
        r_Q_star(1) = r_Q_rot_r(1)*(r_S_r-r_u_norm_r) / (r_S_r-r_S_star)
        r_Q_star(2) = r_Q_rot_r(1)*(r_S_r-r_u_norm_r) / (r_S_r-r_S_star)*r_S_star
        r_Q_star(3) = r_Q_rot_r(1)*(r_S_r-r_u_norm_r) / (r_S_r-r_S_star)*r_u_tang_r
        r_flux_hll  = r_flux_r(:,1) + r_S_r*(r_Q_star-r_Q_rot_r)
      ELSEIF (r_S_r < 0.0_GRID_SR) THEN
        r_flux_hll = r_flux_r(:,1)
      END IF

!--- rotation back (just for the 'velocity'-components)
      r_flux(1) = r_flux_hll(1)
      r_flux(2) = r_flux_hll(2)*r_normal(1) - r_flux_hll(3)*r_normal(2)
      r_flux(3) = r_flux_hll(2)*r_normal(2) + r_flux_hll(3)*r_normal(1)

    ELSE
      r_flux = 0.0_GRID_SR
    END IF

  END FUNCTION riemannsolver

!*******************************************************************************
END MODULE DG_riemann_solver
