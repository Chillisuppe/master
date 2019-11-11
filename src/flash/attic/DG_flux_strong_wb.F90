!*******************************************************************************
!
!> @file  DG_flux_strong_wb.F90
!> @brief contains module DG_flux
!
!*******************************************************************************
! MODULE DESCRIPTION
!> @brief provides routine for computing the right-hand side of the resulting ODE
!>        for shallow water equations in strong form.
!>
MODULE DG_flux

  USE GRID_api
  USE FLASH_parameters
  USE DG_boundary
  USE DG_riemann_solver
  USE DG_utils
  USE DG_friction
  USE DG_equation

  PRIVATE :: elmtflux, edgeflux
  PUBLIC  :: fvm_flux

  CONTAINS

!*******************************************************************************
! DESCRIPTION of [SUBROUTINE fvm_flux]:
!> @brief computes right-hand side of resulting ODE
!>
!> @param[out]    r_flux            DIMENSION(3,i_numelmt*i_faceunknowns);
!>                                  values of right-hand side of the
!>                                  equation /f$ U_t = H(U)/f$
!> @param[in]     r_Q               DIMENSION(3,i_numelmt*i_faceunknowns);
!>                                  state vector (height and momentum)
!> @param[in]     r_S               DIMENSION(4,i_numelmt*i_faceunknowns);
!>                                  bathymetry, Coriolis forcing, wind component in x and y direction
!> @param         p_phy             equation dependent parameter structure
!> @param         i_numelmt         total number of elements
!> @param         i_numedge         total number of edges
!> @param         i_faceunknowns    number of DOFs per element
!> @param         r_metrics_inv     metric terms for each element; DIMENSION(2,2,i_numelmt)
!> @param         r_Dxi             DIMENSION(i_faceunknowns, i_faceunknowns);
!>                                  differentiated basis functions in x direction at Lagrange points
!> @param         r_Deta            DIMENSION(i_faceunknowns, i_faceunknowns);
!>                                  differentiated basis functions in y direction at Lagrange points
!> @param         r_epsi            DIMENSION(i_equadpts, i_faceunknowns);
!>                                  matrix of basis function evaluation at element quadrature points
!> @param         r_gpsi            DIMENSION(i_gpsinonzero, i_gquadpts);
!>                                  matrix of basis function evaluation at edge quadrature points
!> @param         r_eMinvpsi        DIMENSION(i_faceunknowns, i_equadpts);
!>                                  inverse mass matrix interpolated on element quadrature points
!> @param         r_eMinvdpsidxi
!> @param         r_eMinvdpsideta
!> @param         r_gMinvpsi        DIMENSION(i_faceunknowns,i_gquadpts);
!>                                  basis functions evaluated at the edge quadrature points
!> @param         i_equadpts        number of quadrature points per element
!> @param         i_gquadpts        number of quadrature points per edge
!> @param         r_eqwei           weights for the quadrature rule corresponding to the element quadrature points; DIMENSION(i_equadpts)
!> @param         r_gqwei           weights for the quadrature rule corresponding to the edge quadrature points; DIMENSION(i_gquadpts)
!> @param         r_element_vol     volume of each element; DIMENSION(i_numelmt)
!> @param         r_edgelength      DIMENSION(i_numedge);
!>                                  length of each edge
!> @param         r_normals         DIMENSION(2, i_numedge);
!>                                  normals per edge
!> @param         i_edgeinfo        edge-element relation for each edge: local edge index and global element index for adjacent elements; DIMENSION(4,i_numedge)
!> @param         i_edgeboundary    boundary flags for all edges; DIMENSION(i_numedge)
!> @param         i_elementdofs     DOF indices corresponding to each element; DIMENSION(i_faceunknowns, i_numelmt)
!> @param         r_coodof
!
  SUBROUTINE fvm_flux(r_flux, r_Q, r_S, p_phy, &
                      i_numelmt, i_numedge, i_faceunknowns, &
                      r_metrics_inv, r_Dxi, r_Deta, r_epsi, r_gpsi, &
                      r_eMinvpsi, r_eMinvdpsidxi, r_eMinvdpsideta, r_gMinvpsi, &
                      i_equadpts, i_gquadpts, r_eqwei, r_gqwei, r_element_vol, &
                      r_edgelength, r_normals, &
                      i_edgeinfo, i_edgeboundary, i_elementdofs, r_coodof)

    IMPLICIT NONE

    REAL (KIND = GRID_SR), DIMENSION(:,:), INTENT(OUT)      :: r_flux
    REAL (KIND = GRID_SR), DIMENSION(:,:), INTENT(IN)       :: r_Q, r_S, r_Dxi, r_Deta, &
                                                               r_epsi, r_gpsi, r_eMinvpsi, &
                                                               r_eMinvdpsidxi, r_eMinvdpsideta, &
                                                               r_gMinvpsi, r_normals, r_coodof
    TYPE (phys_param), INTENT(IN)                           :: p_phy
    INTEGER (KIND = GRID_SI), INTENT(IN)                    :: i_numelmt, i_numedge, &
                                                               i_faceunknowns, &
                                                               i_equadpts, i_gquadpts
    REAL (KIND = GRID_SR), DIMENSION(:,:,:), INTENT(IN)     :: r_metrics_inv
    REAL (KIND = GRID_SR), DIMENSION(:), INTENT(IN)         :: r_eqwei, r_gqwei, r_element_vol, r_edgelength
    INTEGER (KIND = GRID_SI), DIMENSION(:,:), INTENT(IN)    :: i_edgeinfo, i_elementdofs
    INTEGER (KIND = GRID_SI), DIMENSION(:), INTENT(IN)      :: i_edgeboundary

    INTEGER (KIND = GRID_SI)                                :: i_elmt, i_elmt_l, i_elmt_r, &
                                                               i_eloc_l, i_eloc_r, i_cnt, &
                                                               i_edge, i_gpsinonzero
    REAL (KIND = GRID_SR)                                   :: r_edge_length, &
                                                               r_elmtvol_e, r_elmtvol_k, &
                                                               r_wet, r_minh_l, r_minh_r
    REAL (KIND = GRID_SR), DIMENSION(GRID_DIMENSION)        :: r_normal
    INTEGER (KIND = GRID_SI), DIMENSION(:), ALLOCATABLE     :: i_gpsiidx
    INTEGER (KIND = GRID_SI), DIMENSION(:,:), ALLOCATABLE   :: i_rotate, i_reflect, i_psiidx_l, i_psiidx_r
    REAL (KIND = GRID_SR), DIMENSION(:), ALLOCATABLE        :: r_h, r_b, r_h_e, r_h_k, &
                                                               r_hu_e, r_hu_k, &
                                                               r_hv_e, r_hv_k, &
                                                               r_b_e, r_b_k, &
                                                               r_h_q, r_hu_q, r_hv_q, r_cor_q, &
                                                               r_laphu_q, r_laphv_q, &
                                                               r_h_lq, r_hu_lq, r_hv_lq, r_b_lq, &
                                                               r_h_rq, r_hu_rq, r_hv_rq, r_b_rq
    REAL (KIND = GRID_SR), DIMENSION(:,:), ALLOCATABLE      :: r_rhs, r_rhs_l, r_rhs_r, &
                                                               r_gradb_q, r_tau_q, &
                                                               r_gradh_q, r_gradhu_q, r_gradhv_q, &
                                                               r_DxiDxi, r_DetaDeta, r_DxiDeta, &
                                                               r_ddx, r_ddy, r_lap, &
                                                               r_ddxq, r_ddyq, r_lapq
!-- Added for wellbalancing
    REAL (KIND=GRID_SR), DIMENSION(i_faceunknowns)          :: r_h_wb, r_b_taylor, &
                                                               r_cooq_x, r_cooq_y, r_x, r_y
    REAL (KIND=GRID_SR)                                     :: r_minh, r_maxh, r_bat_check, &
                                                               r_eps=10E-6

!--- allocate workspace

    i_gpsinonzero = GRID_femtypes%p_type(FEM_ZETA)%sig%i_gpsinonzero

    ALLOCATE(r_DxiDxi(i_faceunknowns, i_faceunknowns), &
             r_DxiDeta(i_faceunknowns, i_faceunknowns), &
             r_DetaDeta(i_faceunknowns, i_faceunknowns), &
             r_ddx(i_faceunknowns, i_faceunknowns), &
             r_ddy(i_faceunknowns, i_faceunknowns), &
             r_lap(i_faceunknowns, i_faceunknowns), &
             r_ddxq(i_faceunknowns, i_equadpts), &
             r_ddyq(i_faceunknowns, i_equadpts), &
             r_lapq(i_faceunknowns, i_equadpts), &
             i_gpsiidx(i_gpsinonzero), &
             i_psiidx_l(i_faceunknowns, GRID_elementnodes), &
             i_psiidx_r(i_faceunknowns, GRID_elementnodes), &
             i_rotate(i_faceunknowns, GRID_elementnodes), &
             i_reflect(i_faceunknowns, 2))

    r_DxiDxi   = MATMUL(r_Dxi , r_Dxi )
    r_DxiDeta  = MATMUL(r_Dxi , r_Deta)
    r_DetaDeta = MATMUL(r_Deta, r_Deta)

    i_gpsiidx  = GRID_femtypes%p_type(FEM_ZETA)%sig%i_gpsiidx
    i_rotate   = GRID_femtypes%p_type(FEM_ZETA)%sig%i_rotation
    i_reflect  = GRID_femtypes%p_type(FEM_ZETA)%sig%i_reflection
    i_psiidx_l = i_rotate
    i_psiidx_r(:,1) = i_reflect(i_rotate(:,1),2)
    i_psiidx_r(:,2) = i_reflect(i_rotate(:,3),2)
    i_psiidx_r(:,3) = i_reflect(i_rotate(:,2),2)

!--- initialize local quantities
    ALLOCATE(r_rhs(3,i_faceunknowns), r_h(i_faceunknowns), r_b(i_faceunknowns), &
             r_h_q(i_equadpts), r_hu_q(i_equadpts), r_hv_q(i_equadpts), &
             r_gradh_q(i_equadpts,2), r_gradhu_q(i_equadpts,2), r_gradhv_q(i_equadpts,2), &
             r_gradb_q(i_equadpts,2), r_cor_q(i_equadpts), r_tau_q(i_equadpts,2), &
             r_laphu_q(i_equadpts), r_laphv_q(i_equadpts))

    r_flux = 0._GRID_SR

!--- element loop for inner flux computation
    elmt_loop: DO i_elmt = 1,i_numelmt

!--- compute derivative operators
      r_ddx(:,:) = 0
      r_ddy(:,:) = 0
      r_lap(:,:) = 0
      DO i_cnt = 1, i_faceunknowns
        r_ddx(:,i_cnt) = r_Dxi(:,i_cnt)*r_metrics_inv(1,1,i_elmt) + r_Deta(:,i_cnt)*r_metrics_inv(1,2,i_elmt)
        r_ddy(:,i_cnt) = r_Dxi(:,i_cnt)*r_metrics_inv(2,1,i_elmt) + r_Deta(:,i_cnt)*r_metrics_inv(2,2,i_elmt)
        r_lap(:,i_cnt) = r_DxiDxi(:,i_cnt)  *r_metrics_inv(1,1,i_elmt)*r_metrics_inv(1,1,i_elmt) + &
                         r_DetaDeta(:,i_cnt)*r_metrics_inv(1,2,i_elmt)*r_metrics_inv(1,2,i_elmt) + &
                         2.0*r_DxiDeta(:,i_cnt) *r_metrics_inv(1,1,i_elmt)*r_metrics_inv(1,2,i_elmt) + &
                         r_DxiDxi(:,i_cnt)  *r_metrics_inv(2,1,i_elmt)*r_metrics_inv(2,1,i_elmt) + &
                         r_DetaDeta(:,i_cnt)*r_metrics_inv(2,2,i_elmt)*r_metrics_inv(2,2,i_elmt) + &
                         2.0*r_DxiDeta(:,i_cnt) *r_metrics_inv(2,1,i_elmt)*r_metrics_inv(2,2,i_elmt)
      END DO
      r_ddxq = MATMUL(r_ddx, r_epsi)
      r_ddyq = MATMUL(r_ddy, r_epsi)
      r_lapq = MATMUL(r_lap, r_epsi)

!--- compute unknowns and their derivatives at quadrature points
      r_h_q  = MATMUL(r_Q(1,i_elementdofs(:,i_elmt)), r_epsi)
      r_hu_q = MATMUL(r_Q(2,i_elementdofs(:,i_elmt)), r_epsi)
      r_hv_q = MATMUL(r_Q(3,i_elementdofs(:,i_elmt)), r_epsi)

      DO i_cnt = 1, i_faceunknowns
        IF(r_h_q(i_cnt) .LE. 10E-6) THEN
          r_h_q(i_cnt)  = 0._GRID_SR
          r_hu_q(i_cnt) = 0._GRID_SR
          r_hv_q(i_cnt) = 0._GRID_SR
        END IF
      END DO

      r_gradh_q(:,1)  = MATMUL(r_Q(1,i_elementdofs(:,i_elmt)), r_ddxq)
      r_gradh_q(:,2)  = MATMUL(r_Q(1,i_elementdofs(:,i_elmt)), r_ddyq)

      r_gradhu_q(:,1) = MATMUL(r_Q(2,i_elementdofs(:,i_elmt)), r_ddxq)
      r_gradhu_q(:,2) = MATMUL(r_Q(2,i_elementdofs(:,i_elmt)), r_ddyq)

      r_gradhv_q(:,1) = MATMUL(r_Q(3,i_elementdofs(:,i_elmt)), r_ddxq)
      r_gradhv_q(:,2) = MATMUL(r_Q(3,i_elementdofs(:,i_elmt)), r_ddyq)

      r_wet = 1.0_GRID_SR !- heaviside(MINVAL(r_Q(1,i_elementdofs(:,i_elmt)))-1E-10)

      r_gradb_q(:,1) = MATMUL(r_S(1,i_elementdofs(:,i_elmt)), r_ddxq)
      r_gradb_q(:,2) = MATMUL(r_S(1,i_elementdofs(:,i_elmt)), r_ddyq)
      r_cor_q        = MATMUL(r_S(2,i_elementdofs(:,i_elmt)), r_epsi)
      r_tau_q(:,1)   = MATMUL(r_S(3,i_elementdofs(:,i_elmt)), r_epsi)
      r_tau_q(:,2)   = MATMUL(r_S(4,i_elementdofs(:,i_elmt)), r_epsi)
      r_laphu_q      = MATMUL(r_Q(2,i_elementdofs(:,i_elmt)), r_lapq)
      r_laphv_q      = MATMUL(r_Q(3,i_elementdofs(:,i_elmt)), r_lapq)

!-- Fix for wellbalancing
      r_h = r_Q(1,i_elementdofs(:,i_elmt))
      r_b = r_S(1,i_elementdofs(:,i_elmt))
      r_minh = MINVAL(r_h)
      r_maxh = MAXVAL(r_h)

      IF ((r_minh .LE. 10E-10) .AND.(r_maxh .GT. 10E-3) .AND. &
          ((MAXVAL(abs(r_gradb_q(:,1))) .GE. 10E-10) .OR. &
           (MAXVAL(abs(r_gradb_q(:,2))) .GE. 10E-10))   .AND. &
          (MAXVAL(abs(r_hu_q)) .LE. 10E-13)             .AND. &
          (MAXVAL(abs(r_hv_q)) .LE. 10E-13)) THEN

!-- Get coordinates for DOF
        r_x = r_coodof(1,i_elementdofs(:,i_elmt))
        r_y = r_coodof(2,i_elementdofs(:,i_elmt))

!-- Get coordinates of quadrature point i_face
        r_cooq_x = MATMUL(r_x, r_epsi)
        r_cooq_y = MATMUL(r_y, r_epsi)

!-- Obtain modified fluid height
        r_b_taylor  = r_b + &
                      r_gradb_q(:,1)*(r_cooq_x - r_x) + &
                      r_gradb_q(:,2)*(r_cooq_y - r_y)
        r_h_wb      = MAX(0._GRID_SR, r_h-r_b_taylor)

        r_bat_check = MINVAL(r_h+r_b) - MAXVAL(r_b)

    !    r_Q(1,i_elementdofs(:,i_elmt)) = r_h_wb
    !    r_h_q = r_h_wb

        IF (MAXVAL(abs(r_h_wb)) .LE. 10E-8) THEN
          r_hu_q = 0._GRID_SR
          r_hv_q = 0._GRID_SR
        END IF

!-- Fix gradient erstmal nur, wenn die bathy "rausguckt"
        IF (r_bat_check .LT. r_eps) THEN
          r_wet = 0.0_GRID_SR
          r_gradb_q = r_gradh_q
        END IF
      END IF

!--- evaluate inner element flux
      CALL elmtflux(r_rhs, i_faceunknowns, i_equadpts, r_eqwei, r_eMinvpsi, &
                    r_h_q, r_hu_q, r_hv_q, r_gradh_q, r_gradhu_q, r_gradhv_q, r_wet, &
                    p_phy, r_gradb_q, r_cor_q, r_tau_q, r_laphu_q, r_laphv_q )

      r_flux(:,i_elementdofs(:,i_elmt)) = r_rhs

    END DO elmt_loop

    DEALLOCATE(r_rhs, r_h, r_b, r_h_q, r_hu_q, r_hv_q, &
               r_gradh_q, r_gradhu_q, r_gradhv_q, &
               r_gradb_q, r_cor_q, r_tau_q, r_laphu_q, r_laphv_q)

    ALLOCATE(r_rhs_l(3,i_faceunknowns), r_rhs_r(3,i_faceunknowns), &
             r_h_e(i_faceunknowns),     r_h_k(i_faceunknowns), &
             r_b_e(i_faceunknowns),     r_b_k(i_faceunknowns), &
             r_hu_e(i_faceunknowns),    r_hu_k(i_faceunknowns), &
             r_hv_e(i_faceunknowns),    r_hv_k(i_faceunknowns), &
             r_h_lq(i_gquadpts), r_hu_lq(i_gquadpts), r_hv_lq(i_gquadpts), r_b_lq(i_gquadpts), &
             r_h_rq(i_gquadpts), r_hu_rq(i_gquadpts), r_hv_rq(i_gquadpts), r_b_rq(i_gquadpts))

!-- edge loop for boundary flux computation
    edge_loop: DO i_edge= 1, i_numedge
    IF((i_edgeboundary(i_edge) <= 0) .OR. (i_edgeboundary(i_edge) == i_edge)) THEN
      i_eloc_l = i_edgeinfo(1,i_edge)
      i_eloc_r = i_edgeinfo(2,i_edge)
      i_elmt_l = i_edgeinfo(3,i_edge)
      i_elmt_r = i_edgeinfo(4,i_edge)
      IF(i_edgeboundary(i_edge) < 0) THEN
        i_elmt_r = i_elmt_l
        i_eloc_r = i_eloc_l
      ENDIF

      r_elmtvol_e = r_element_vol(i_elmt_l)
      r_h_e       = r_Q(1,i_elementdofs(i_psiidx_l(:,i_eloc_l),i_elmt_l))
      r_hu_e      = r_Q(2,i_elementdofs(i_psiidx_l(:,i_eloc_l),i_elmt_l))
      r_hv_e      = r_Q(3,i_elementdofs(i_psiidx_l(:,i_eloc_l),i_elmt_l))
      r_b_e       = r_S(1,i_elementdofs(i_psiidx_l(:,i_eloc_l),i_elmt_l))

!-- correct momentum for small water heights
      DO i_cnt = 1, i_faceunknowns
        IF(r_h_e(i_cnt) .LE. 10E-6) THEN
          r_h_e(i_cnt)  = 0._GRID_SR
          r_hu_e(i_cnt) = 0._GRID_SR
          r_hv_e(i_cnt) = 0._GRID_SR
        END IF
      END DO

      r_normal      = r_normals(:,i_edge)
      r_edge_length = r_edgelength(i_edge)

!-- compute values for the neighboring element
      CALL DG_bou(i_elmt_l, i_elmt_r, i_faceunknowns, i_elementdofs(i_psiidx_r(:,i_eloc_r),i_elmt_r), &
                  i_edgeboundary(i_edge), r_normal, r_Q, r_S(1,:), r_element_vol, &
                  r_h_e, r_hu_e, r_hv_e, r_b_e, r_h_k, r_hu_k, r_hv_k, r_b_k, r_elmtvol_k)

      r_minh_l = MINVAL(r_h_e)
      r_minh_r = MINVAL(r_h_k)

      r_h_lq  = MATMUL(r_h_e(i_gpsiidx) , r_gpsi)
      r_hu_lq = MATMUL(r_hu_e(i_gpsiidx), r_gpsi)
      r_hv_lq = MATMUL(r_hv_e(i_gpsiidx), r_gpsi)
      r_b_lq  = MATMUL(r_b_e(i_gpsiidx) , r_gpsi)

      r_h_rq  = MATMUL(r_h_k(i_gpsiidx) , r_gpsi)
      r_hu_rq = MATMUL(r_hu_k(i_gpsiidx), r_gpsi)
      r_hv_rq = MATMUL(r_hv_k(i_gpsiidx), r_gpsi)
      r_b_rq  = MATMUL(r_b_k(i_gpsiidx) , r_gpsi)

      CALL edgeflux(r_rhs_l, r_rhs_r, r_normal, r_minh_l, r_minh_r, &
                    i_faceunknowns, i_gquadpts, r_gqwei, r_gMinvpsi, &
                    r_h_lq, r_hu_lq, r_hv_lq, r_b_lq, r_h_rq, r_hu_rq, r_hv_rq, r_b_rq)

      r_flux(:,i_elementdofs(i_psiidx_l(:,i_eloc_l),i_elmt_l)) = &
        r_flux(:,i_elementdofs(i_psiidx_l(:,i_eloc_l),i_elmt_l)) + &
        r_edge_length/r_elmtvol_e*r_rhs_l

!--- if we have an inner or periodic edge, also update right element
      IF(i_edgeboundary(i_edge) >= 0) &
        r_flux(:,i_elementdofs(i_psiidx_r(:,i_eloc_r),i_elmt_r)) = &
          r_flux(:,i_elementdofs(i_psiidx_r(:,i_eloc_r),i_elmt_r)) - &
          r_edge_length/r_elmtvol_k*r_rhs_r

    END IF
    END DO edge_loop

!--- deallocate workspace
    DEALLOCATE(r_rhs_l, r_rhs_r, r_h_e, r_h_k, &
               r_hu_e, r_hu_k, r_hv_e, r_hv_k, r_b_e, r_b_k, &
               r_h_lq, r_hu_lq, r_hv_lq, r_b_lq, r_h_rq, r_hu_rq, r_hv_rq, r_b_rq)

    DEALLOCATE(r_DxiDxi, r_DxiDeta, r_DetaDeta, &
               r_ddx, r_ddy, r_lap, r_ddxq, r_ddyq, r_lapq, &
               i_gpsiidx, i_psiidx_l, i_psiidx_r, i_rotate, i_reflect)

  END SUBROUTINE fvm_flux

!*******************************************************************************
! DESCRIPTION of [SUBROUTINE elmtflux]:
!> @brief computes the divergence and source terms of the equation in consistent form
!
!> @param[out]    r_rhs           DIMENSION(3);
!>                                sum of divergence and source term per quadrature point
!> @param[in]     i_faceunknowns  number of DOFs per element
!> @param[in]     i_equadpts      number of quadrature points per element
!> @param[in]     r_eqwei         weights for the quadrature rule corresponding to the element quadrature points; DIMENSION(i_equadpts)
!> @param[in]     r_eMinvpsi      DIMENSION(i_faceunknowns, i_equadpts);
!>                                inverse mass matrix interpolated on element quadrature points
!> @param[in]     r_h             fluid depth at element quadrature points; DIMENSION(i_equadpts)
!> @param[in]     r_hu            momentum in x-direction at element quadrature points; DIMENSION(i_equadpts)
!> @param[in]     r_hv            momentum in y-direction at element quadrature points; DIMENSION(i_equadpts)
!> @param[in]     r_gradh         DIMENSION(i_equadpts, 2);
!>                                gradient of fluid height at quadrature points
!> @param[in]     r_gradhu        DIMENSION(i_equadpts, 2);
!>                                gradient of x-momentum at quadrature points
!> @param[in]     r_gradhv        DIMENSION(i_equadpts, 2);
!>                                gradient of y-momentum at quadrature points
!> @param[in]     r_wet
!> @param[in]     p_phy           equation dependent parameter structure
!> @param[in]     r_gradb         DIMENSION(i_equadpts, 2);
!>                                gradient of bottom topography at quadrature points
!> @param[in]     r_cor           DIMENSION(i_equadpts);
!>                                Coriolis forcing at quadrature points
!> @param[in]     r_tau           DIMENSION(i_equadpts, 2);
!>                                vector of wind forcing at quadrature points
!> @param[in]     r_laphu         DIMENSION(i_faceunknowns, i_equadpts);
!>                                Laplacian of hu at quadrature points
!> @param[in]     r_laphv         DIMENSION(i_faceunknowns, i_equadpts);
!>                                Laplacian of hv at quadrature points
!
  SUBROUTINE elmtflux(r_rhs, i_faceunknowns, i_equadpts, r_eqwei, r_eMinvpsi, &
                      r_h, r_hu, r_hv, r_gradh, r_gradhu, r_gradhv, r_wet, &
                      p_phy, r_gradb, r_cor, r_tau, r_laphu, r_laphv)

    IMPLICIT NONE

    REAL (KIND = GRID_SR), DIMENSION(:,:), INTENT(out)      :: r_rhs
    INTEGER (KIND = GRID_SI), INTENT(in)                    :: i_faceunknowns, i_equadpts
    REAL (KIND = GRID_SR), DIMENSION(:), INTENT(in)         :: r_eqwei, r_h, r_hu, r_hv
    REAL (KIND = GRID_SR), DIMENSION(:,:), INTENT(in)       :: r_gradh, r_gradhu, r_gradhv, r_eMinvpsi
    REAL (KIND = GRID_SR), INTENT(in)                       :: r_wet
    TYPE (phys_param), INTENT(in)                           :: p_phy
    REAL (KIND = GRID_SR), DIMENSION(:), INTENT(in)         :: r_cor, r_laphu, r_laphv
    REAL (KIND = GRID_SR), DIMENSION(:,:), INTENT(in)       :: r_gradb, r_tau

    INTEGER (KIND = GRID_SI)                                :: i_quad, i_dof
    REAL (KIND = GRID_SR)                                   :: r_heavi
    REAL (KIND = GRID_SR), DIMENSION(3)                     :: r_Div

    r_rhs(:,:) = 0.0_GRID_SR

!-- quadrature loop
    elmt_quad_loop: DO i_quad=1, i_equadpts

      r_heavi = heaviside(abs(r_h(i_quad)))

!--- divergence terms
      r_Div = divflux(r_h(i_quad), r_hu(i_quad), r_hv(i_quad), r_wet, r_gradh(i_quad,:), &
                      r_gradhu(i_quad,:), r_gradhv(i_quad,:))

      r_Div = r_Div+source(r_h(i_quad), r_hu(i_quad), r_hv(i_quad), r_heavi, r_cor(i_quad), &
                           r_wet, r_gradb(i_quad,:), r_tau(i_quad,:), 0.0_GRID_SR, &
                           r_laphu(i_quad), r_laphv(i_quad), p_phy)

!--- multiply with r_eqwei*Minv*psi for each dof
      elmt_dof_loop: DO i_dof=1,i_faceunknowns
        r_rhs(:,i_dof) = r_rhs(:,i_dof) - r_eqwei(i_quad)*r_eMinvpsi(i_dof,i_quad)*r_Div
      END DO elmt_dof_loop

    END DO elmt_quad_loop

  END SUBROUTINE elmtflux

!*******************************************************************************
! DESCRIPTION of [SUBROUTINE edgeflux]:
!> @brief computation of edge-based terms
!>
!> @param[out]    r_rhs_l         edge flux /f$ (F - F^{\ast})\cdot n /f$ for left element; DIMENSION(3,i_faceunknowns)
!> @param[out]    r_rhs_r         edge flux /f$ (F - F^{\ast})\cdot m /f$ for right element; DIMENSION(3,i_faceunknowns)
!> @param[in]     r_normal        DIMENSION(2);
!>                                local normal vector pointing from left to right element
!> @param[in]     r_minh_l
!> @param[in]     r_minh_r
!> @param[in]     i_faceunknowns  number of DOFs per element
!> @param[in]     i_gquadpts      number of quadrature points per edge
!> @param[in]     r_gqwei         weights for the quadrature rule corresponding to the edge quadrature points; DIMENSION(i_gquadpts)
!> @param[in]     r_gMinvpsi      DIMENSION(i_faceunknowns, i_equadpts);
!>                                inverse mass matrix interpolated on edge quadrature points
!> @param[in]     r_h_l           fluid depth at edge quadrature points of left element; DIMENSION(i_gquadpts)
!> @param[in]     r_hu_l          momentum in x-direction at edge quadrature points of left element; DIMENSION(i_gquadpts)
!> @param[in]     r_hv_l          momentum in y-direction at edge quadrature points of left element; DIMENSION(i_gquadpts)
!> @param[in]     r_h_r           fluid depth at edge quadrature points of right element; DIMENSION(i_gquadpts)
!> @param[in]     r_hu_r          momentum in x-direction at edge quadrature points of right element; DIMENSION(i_gquadpts)
!> @param[in]     r_hv_r          momentum in y-direction at edge quadrature points of right element; DIMENSION(i_gquadpts)
!
  SUBROUTINE edgeflux(r_rhs_l, r_rhs_r, r_normal, r_minh_l, r_minh_r, &
                      i_faceunknowns, i_gquadpts, r_gqwei, r_gMinvpsi, &
                      r_h_l, r_hu_l, r_hv_l, r_b_l, r_h_r, r_hu_r, r_hv_r, r_b_r)

    IMPLICIT NONE

    REAL (KIND = GRID_SR), DIMENSION(:,:), INTENT(out)      :: r_rhs_l, r_rhs_r
    REAL (KIND = GRID_SR), DIMENSION(2), INTENT(in)         :: r_normal
    REAL (KIND = GRID_SR), INTENT(in)                       :: r_minh_l, r_minh_r
    INTEGER (KIND = GRID_SI), INTENT(in)                    :: i_faceunknowns, i_gquadpts
    REAL (KIND = GRID_SR), DIMENSION(:), INTENT(in)         :: r_gqwei, &
                                                               r_h_l, r_hv_l, r_hu_l, r_b_l, &
                                                               r_h_r, r_hv_r, r_hu_r, r_b_r
    REAL (KIND = GRID_SR), DIMENSION(:,:), INTENT(in)       :: r_gMinvpsi

!--- local declarations
    INTEGER (KIND = GRID_SI)                                :: i_quad, i_dof
    REAL (KIND = GRID_SR), DIMENSION(3,2)                   :: r_F_l, r_F_r
    REAL (KIND = GRID_SR), DIMENSION(3)                     :: r_Fstar, r_Fleft, r_Frght

    r_rhs_l(:,:) = 0.0_GRID_SR
    r_rhs_r(:,:) = 0.0_GRID_SR

!--- perform edge quadrature
    edge_quad_loop: DO i_quad=1,i_gquadpts

      r_F_l = flux(r_h_l(i_quad), r_hu_l(i_quad), r_hv_l(i_quad))
      r_F_r = flux(r_h_r(i_quad), r_hu_r(i_quad), r_hv_r(i_quad))

      r_Fleft = MATMUL(r_F_l, r_normal)
      r_Frght = MATMUL(r_F_r, r_normal)

      r_Fstar  = riemannsolver(r_h_l(i_quad), r_h_r(i_quad), r_hu_l(i_quad), r_hu_r(i_quad), &
                               r_hv_l(i_quad), r_hv_r(i_quad), r_normal)

      edge_dof_loop: DO i_dof=1,i_faceunknowns
        r_rhs_l(:,i_dof) = r_rhs_l(:,i_dof) - r_gqwei(i_quad)*r_gMinvpsi(i_dof,i_quad)*(r_Fstar-r_Fleft)
        r_rhs_r(:,i_dof) = r_rhs_r(:,i_dof) - r_gqwei(i_quad)*r_gMinvpsi(i_dof,i_quad)*(r_Fstar-r_Frght)
      END DO edge_dof_loop
    END DO edge_quad_loop

  END SUBROUTINE edgeflux

!*******************************************************************************
END MODULE DG_flux
