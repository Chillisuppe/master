!*******************************************************************************
!
!> @file  DG_RS_roe.F90
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
!> @brief Roe Riemann solver for shallow water equations
!>
!> @note does not work for wet/ dry problem
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
                                                                     r_a_r, r_a_l, r_a_roe, &
                                                                     r_h_roe, r_u_roe, r_v_roe, &
                                                                     r_h_star_l, r_u_star_l, r_v_star_l, &
                                                                     r_h_star_r, r_u_star_r, r_v_star_r, &
                                                                     r_lambda1_l, r_lambda1_r, &
                                                                     r_lambda3_l, r_lambda3_r, &
                                                                     r_lambda_bar, r_tol
    REAL (KIND = GRID_SR), DIMENSION(i_nprogvars)                 :: r_Q_rot_l, r_Q_rot_r, &
                                                                     r_fluxtmp, r_lambda, &
                                                                     r_alpha, r_delta
    REAL (KIND = GRID_SR), DIMENSION(3,3)                         :: r_K
    REAL (KIND = GRID_SR), DIMENSION(i_nprogvars,GRID_DIMENSION)  :: r_flux_l, r_flux_r
    INTEGER (KIND = GRID_SI)                                      :: i_cnt

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

!--- compute Roe averages following the Roe-Pike approach (the passive scalar is v)
      r_h_roe = SQRT(r_Q_rot_l(1)*r_Q_rot_r(1))
      r_u_roe = (r_a_l*r_u_norm_l + r_a_r*r_u_norm_r) / (r_a_l + r_a_r)
      r_v_roe = (r_a_l*r_u_tang_l + r_a_r*r_u_tang_r) / (r_a_l + r_a_r)

!--- Roe sound speed
      r_a_roe = SQRT(0.5_GRID_SR*(r_a_l**2+r_a_r**2))

!--- compute jumps
      r_delta(1) = r_Q_rot_r(1)            - r_Q_rot_l(1)
      r_delta(2) = r_Q_rot_r(1)*r_u_norm_r - r_Q_rot_l(1)*r_u_norm_l
      r_delta(3) = r_Q_rot_r(1)*r_u_tang_r - r_Q_rot_l(1)*r_u_tang_l

!--- compute eigenvalues (u*n+a)
      r_lambda(1) = r_u_roe - r_a_roe
      r_lambda(2) = r_u_roe
      r_lambda(3) = r_u_roe + r_a_roe

!--- compute rotated eigenvectors (ith one = K(:,i))
      r_K(:,1) =  [1.0_GRID_SR, r_lambda(1), r_v_roe    ]
      r_K(:,2) =  [0.0_GRID_SR, 0.0_GRID_SR, 1.0_GRID_SR]
      r_K(:,3) =  [1.0_GRID_SR, r_lambda(3), r_v_roe    ]

!--- matrix of left eigenvectors
!          (  r_lambda(3)/(2*r_a_roe)  -1/(2*r_a_roe)   0  )
!       L= (               -r_v_roe                0    1  )
!          ( -r_lambda(1)/(2*r_a_roe)   1/(2*r_a_roe)   0  )

!--- compute wave strength alpha (computed with Maple K^{-1}*(ul*n-ur*n))
      r_alpha(1) = 0.5_GRID_SR*( r_lambda(3)*r_delta(1)-r_delta(2)) / r_a_roe
      r_alpha(2) = -r_v_roe*r_delta(1)+r_delta(3)
      r_alpha(3) = 0.5_GRID_SR*(-r_lambda(1)*r_delta(1)+r_delta(2)) / r_a_roe

!--- compute left and right fluxes
      r_flux_l = flux(r_Q_rot_l)
      r_flux_r = flux(r_Q_rot_r)

!--- compute Roe flux
      r_fluxtmp = r_flux_l(:,1)
      DO i_cnt = 1,3
        IF (r_lambda(i_cnt) < 0.0_GRID_SR) THEN
          r_fluxtmp = r_fluxtmp + r_alpha(i_cnt)*r_lambda(i_cnt)*r_K(:,i_cnt)
        END IF !r_lambda <=0
      END DO

!--- Harten-Hyman entropy fix?
!--- compute values in star region (U_star=U_L + alpha_1*K_1)
      r_h_star_l = r_Q_rot_l(1) + r_alpha(1)*r_K(1,1)
      r_u_star_l = r_u_norm_l   + r_alpha(1)*r_K(2,1)
      r_v_star_l = r_u_tang_l   + r_alpha(1)*r_K(3,1)

      r_h_star_r = r_Q_rot_r(1) - r_alpha(3)*r_K(1,3)
      r_u_star_r = r_u_norm_r   - r_alpha(3)*r_K(2,3)
      r_v_star_r = r_u_tang_r   - r_alpha(3)*r_K(3,3)

!--- again we need the normal direction
      r_lambda1_l = r_u_norm_l - r_a_l
      r_lambda1_r = r_u_star_l - SQRT(r_h_star_l)
      r_lambda3_l = r_u_norm_r + r_a_r
      r_lambda3_r = r_normal(1)+ SQRT(r_h_star_r)

!--- check for transonic rarefactions (star states with roe averages)
      !-- left transonic rarefaction
      IF ((r_lambda1_l < 0.0_GRID_SR) .AND. (r_lambda1_r > 0.0_GRID_SR)) THEN
        !-- compute intermediate speed (11.128 Toro)
        r_lambda_bar = r_lambda1_l * (r_lambda1_r-r_lambda(1)) / (r_lambda1_r-r_lambda1_l)
        r_fluxtmp    = r_flux_l(:,1)+r_lambda_bar*r_alpha(1)*r_K(:,1)

      !-- right transonic rarefaction
      ELSEIF ((r_lambda3_l < 0.0_GRID_SR) .AND. (r_lambda3_r > 0.0_GRID_SR)) THEN
        r_lambda_bar = r_lambda3_r * (r_lambda(3)-r_lambda3_l) / (r_lambda3_r-r_lambda3_l)
        r_fluxtmp    = r_flux_r(:,1)-r_lambda_bar*r_alpha(3)*r_K(:,3)

      END IF !transonic wave

!--- rotation back (just for the 'velocity'-components)
      r_flux(1) = r_fluxtmp(1)
      r_flux(2) = r_fluxtmp(2)*r_normal(1) - r_fluxtmp(3)*r_normal(2)
      r_flux(3) = r_fluxtmp(2)*r_normal(2) + r_fluxtmp(3)*r_normal(1)

    ELSE
      r_flux = 0.0_GRID_SR
    END IF

  END FUNCTION riemannsolver

!*******************************************************************************
END MODULE DG_riemann_solver
