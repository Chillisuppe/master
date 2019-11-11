!*******************************************************************************
!
!> @file  DG_flux_strong.F90
!> @brief contains module DG_flux
!
!*******************************************************************************
! MODULE DESCRIPTION:
!> @brief provides routine for computing the right-hand side of the resulting ODE
!>        for DG discretization in strong form.
!>
MODULE DG_flux

#ifdef _OPENMP
  USE OMP_lib
#endif
  USE GRID_api
  USE FLASH_parameters
  USE DG_boundary
  USE DG_riemann_solver
  USE DG_equation, ONLY : i_nprogvars, i_nsrcterms, source, flux, divflux

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
!> @param         r_time            current model time
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
!>                                  transformed basis functions at element quadrature points
!> @param         r_eMinvdpsidxi
!> @param         r_eMinvdpsideta
!> @param         r_gMinvpsi        DIMENSION(i_faceunknowns,i_gquadpts);
!>                                  transformed basis functions at edge quadrature points
!> @param         i_equadpts        number of quadrature points per element
!> @param         i_gquadpts        number of quadrature points per edge
!> @param         r_eqwei           weights for the quadrature rule corresponding to the element quadrature points; DIMENSION(i_equadpts)
!> @param         r_gqwei           weights for the quadrature rule corresponding to the edge quadrature points; DIMENSION(i_gquadpts)
!> @param         r_elmtvolume      volume of each element; DIMENSION(i_numelmt)
!> @param         r_edgelength      DIMENSION(i_numedge);
!>                                  length of each edge
!> @param         r_normals         DIMENSION(2, i_numedge);
!>                                  normals per edge
!> @param         i_edgeinfo        edge-element relation for each edge: local edge index and global element index for adjacent elements; DIMENSION(4,i_numedge)
!> @param         i_edgeboundary    boundary flags for all edges; DIMENSION(i_numedge)
!> @param         i_edgenodes       global node indices for each edge; DIMENSION(2,i_numedge)
!> @param         i_elmtdofs        DOF indices corresponding to each element; DIMENSION(i_faceunknowns, i_numelmt)
!> @param         r_coonod          coordinates for each grid node
!> @param         r_coodof          coordinates for each dof in the grid
!
  SUBROUTINE fvm_flux(r_flux, r_Q, r_S, r_time, i_numelmt, i_numedge, i_faceunknowns, &
                      r_metrics_inv, r_Dxi, r_Deta, r_epsi, r_gpsi, &
                      r_eMinvpsi, r_eMinvdpsidxi, r_eMinvdpsideta, r_gMinvpsi, &
                      i_equadpts, i_gquadpts, r_eqwei, r_gqwei, r_elmtvolume, &
                      r_edgelength, r_normals, i_edgeinfo, i_edgeboundary, &
                      i_edgenodes, i_elmtdofs, r_coonod, r_coodof)

    IMPLICIT NONE

    REAL (KIND = GRID_SR), DIMENSION(:,:), INTENT(OUT)      :: r_flux
    REAL (KIND = GRID_SR), DIMENSION(:,:), INTENT(IN)       :: r_Q, r_S, r_Dxi, r_Deta, &
                                                               r_epsi, r_gpsi, r_eMinvpsi, &
                                                               r_eMinvdpsidxi, r_eMinvdpsideta, &
                                                               r_gMinvpsi, r_normals, r_coonod, &
                                                               r_coodof
    REAL (KIND = GRID_SR), INTENT(in)                       :: r_time
    INTEGER (KIND = GRID_SI), INTENT(IN)                    :: i_numelmt, i_numedge, &
                                                               i_faceunknowns, &
                                                               i_equadpts, i_gquadpts
    REAL (KIND = GRID_SR), DIMENSION(:,:,:), INTENT(IN)     :: r_metrics_inv
    REAL (KIND = GRID_SR), DIMENSION(:), INTENT(IN)         :: r_eqwei, r_gqwei, r_elmtvolume, r_edgelength
    INTEGER (KIND = GRID_SI), DIMENSION(:,:), INTENT(IN)    :: i_edgeinfo, i_edgenodes, i_elmtdofs
    INTEGER (KIND = GRID_SI), DIMENSION(:), INTENT(IN)      :: i_edgeboundary

    INTEGER (KIND = GRID_SI)                                :: i_elmt, i_elmt_l, i_elmt_r, &
                                                               i_eloc_l, i_eloc_r, i_eqn, i_face, &
                                                               i_edge
    REAL (KIND = GRID_SR)                                   :: r_edgelen
    REAL (KIND = GRID_SR), DIMENSION(GRID_DIMENSION)        :: r_normal
    INTEGER (KIND = GRID_SI), DIMENSION(GRID_femtypes%p_type(FEM_DG)%sig%i_gpsinonzero) :: i_gpsiidx
    INTEGER (KIND = GRID_SI), DIMENSION(i_faceunknowns)     :: i_ledofs, i_ledofsl, i_ledofsr
    INTEGER (KIND = GRID_SI), DIMENSION(i_faceunknowns, 2)  :: i_reflect
    INTEGER (KIND = GRID_SI), DIMENSION(i_faceunknowns, GRID_elementnodes) :: i_rotate, i_psiidx_l, i_psiidx_r
    REAL (KIND = GRID_SR), DIMENSION(GRID_DIMENSION, i_gquadpts)  :: r_edgequadcoo
    REAL (KIND = GRID_SR), DIMENSION(i_faceunknowns, i_faceunknowns) :: r_DxiDxi, r_DetaDeta, r_DxiDeta, &
                                                                     r_ddx, r_ddy, r_lap
    REAL (KIND = GRID_SR), DIMENSION(i_nprogvars,i_equadpts)      :: r_Q_q, r_dQdx_q, r_dQdy_q
    REAL (KIND = GRID_SR), DIMENSION(i_nprogvars,i_faceunknowns)  :: r_rhs, r_Q_l, r_Q_r, &
                                                                     r_rhs_l, r_rhs_r
    REAL (KIND = GRID_SR), DIMENSION(i_nprogvars,i_gquadpts)      :: r_Q_lq, r_Q_rq
    REAL (KIND = GRID_SR), DIMENSION(i_nsrcterms,i_faceunknowns)  :: r_S_l, r_S_r
    REAL (KIND = GRID_SR), DIMENSION(i_nsrcterms,i_gquadpts)      :: r_S_lq, r_S_rq
    REAL (KIND = GRID_SR), DIMENSION(i_faceunknowns, i_equadpts)  :: r_ddxq, r_ddyq, r_lapq

!--- initialize local quantities
    r_DxiDxi   = MATMUL(r_Dxi , r_Dxi )
    r_DxiDeta  = MATMUL(r_Dxi , r_Deta)
    r_DetaDeta = MATMUL(r_Deta, r_Deta)

    i_gpsiidx  = GRID_femtypes%p_type(FEM_DG)%sig%i_gpsiidx
    i_rotate   = GRID_femtypes%p_type(FEM_DG)%sig%i_rotation
    i_reflect  = GRID_femtypes%p_type(FEM_DG)%sig%i_reflection
    i_psiidx_l = i_rotate
    i_psiidx_r(:,1) = i_reflect(i_rotate(:,1),2)
    i_psiidx_r(:,2) = i_reflect(i_rotate(:,3),2)
    i_psiidx_r(:,3) = i_reflect(i_rotate(:,2),2)

!--- element loop for inner flux computation
    !$OMP PARALLEL DO DEFAULT(SHARED) &
    !$OMP PRIVATE(r_ddx, r_ddy, r_lap, r_ddxq, r_ddyq, r_lapq, i_ledofs, &
    !$OMP r_Q_q, r_rhs, r_dQdx_q, r_dQdy_q) &
    !$OMP SCHEDULE(STATIC)
    elmt_loop: DO i_elmt = 1,i_numelmt

!--- compute derivative operators
      r_ddx = r_Dxi*r_metrics_inv(1,1,i_elmt) + r_Deta*r_metrics_inv(1,2,i_elmt)
      r_ddy = r_Dxi*r_metrics_inv(2,1,i_elmt) + r_Deta*r_metrics_inv(2,2,i_elmt)
      r_lap = r_DxiDxi     *r_metrics_inv(1,1,i_elmt)*r_metrics_inv(1,1,i_elmt) + &
              r_DetaDeta   *r_metrics_inv(1,2,i_elmt)*r_metrics_inv(1,2,i_elmt) + &
              2.0*r_DxiDeta*r_metrics_inv(1,1,i_elmt)*r_metrics_inv(1,2,i_elmt) + &
              r_DxiDxi     *r_metrics_inv(2,1,i_elmt)*r_metrics_inv(2,1,i_elmt) + &
              r_DetaDeta   *r_metrics_inv(2,2,i_elmt)*r_metrics_inv(2,2,i_elmt) + &
              2.0*r_DxiDeta*r_metrics_inv(2,1,i_elmt)*r_metrics_inv(2,2,i_elmt)

      r_ddxq = MATMUL(r_ddx, r_epsi)
      r_ddyq = MATMUL(r_ddy, r_epsi)
      r_lapq = MATMUL(r_lap, r_epsi)

!--- compute unknowns and their derivatives at quadrature points
      i_ledofs = i_elmtdofs(:,i_elmt)
      r_Q_q    = MATMUL(r_Q(:,i_ledofs), r_epsi)
      r_dQdx_q = MATMUL(r_Q(:,i_ledofs), r_ddxq)
      r_dQdy_q = MATMUL(r_Q(:,i_ledofs), r_ddyq)

!--- evaluate inner element flux and source term
      CALL elmtflux(r_rhs, i_faceunknowns, i_equadpts, r_eqwei, r_eMinvpsi, &
                    r_Q_q, r_dQdx_q, r_dQdy_q)
      CALL source(r_rhs, i_faceunknowns, i_equadpts, r_eqwei, r_eMinvpsi, r_epsi, &
                  r_ddxq, r_ddyq, r_lapq, r_Q(:,i_ledofs), r_S(:,i_ledofs), &
                  r_coodof(:,i_ledofs), r_time, 1.0_GRID_SR)

      r_flux(:,i_ledofs) = r_rhs(:,:)

    END DO elmt_loop
    !$OMP END PARALLEL DO

!--- edge loop for boundary flux computation
    !$OMP PARALLEL DO DEFAULT(SHARED) &
    !$OMP PRIVATE(i_eloc_l, i_eloc_r, i_elmt_l, i_elmt_r, i_ledofsl, i_ledofsr, &
    !$OMP r_Q_l, r_S_l, r_Q_r, r_S_r, r_Q_lq, r_S_lq, r_Q_rq, r_S_rq, &
    !$OMP r_normal, r_edgelen, r_edgequadcoo, r_rhs_l, r_rhs_r, i_eqn, i_face) &
    !$OMP SCHEDULE(STATIC)
    edge_loop: DO i_edge= 1, i_numedge
    ! do not visit "slave" edges in case of periodic b.c.
    IF((i_edgeboundary(i_edge) <= 0) .OR. (i_edgeboundary(i_edge) == i_edge)) THEN
      i_eloc_l = i_edgeinfo(1,i_edge)
      i_elmt_l = i_edgeinfo(3,i_edge)

      r_normal  = r_normals(:,i_edge)
      r_edgelen = r_edgelength(i_edge)

      i_ledofsl = i_elmtdofs(i_psiidx_l(:,i_eloc_l),i_elmt_l)
      r_Q_l     = r_Q(:,i_ledofsl)
      r_S_l     = r_S(:,i_ledofsl)
      r_Q_lq    = MATMUL(r_Q_l(:,i_gpsiidx), r_gpsi)
      r_S_lq    = MATMUL(r_S_l(:,i_gpsiidx), r_gpsi)

      ! check if we have an inner or periodic edge
      IF(i_edgeboundary(i_edge) >= 0) THEN
        ! compute values at quadrature points
        i_eloc_r  = i_edgeinfo(2,i_edge)
        i_elmt_r  = i_edgeinfo(4,i_edge)
        i_ledofsr = i_elmtdofs(i_psiidx_r(:,i_eloc_r),i_elmt_r)
        r_Q_r     = r_Q(:,i_ledofsr)
        r_S_r     = r_S(:,i_ledofsr)
        r_Q_rq    = MATMUL(r_Q_r(:,i_gpsiidx), r_gpsi)
        r_S_rq    = MATMUL(r_S_r(:,i_gpsiidx), r_gpsi)

      ELSE
        ! compute boundary conditions at edge quadrature points
        r_edgequadcoo = MATMUL(r_coonod(:,i_edgenodes(:,i_edge)), GRID_femtypes%p_type(FEM_DG)%sig%r_gquadcoo)
        CALL compute_bc(r_time, r_edgequadcoo, i_gquadpts, i_edgeboundary(i_edge), &
                        r_normal, r_Q_lq, r_S_lq, r_Q_rq, r_S_rq)
      ENDIF

      CALL edgeflux(r_rhs_l, r_rhs_r, r_normal, &
                    i_faceunknowns, i_gquadpts, r_gqwei, r_gMinvpsi, &
                    r_Q_lq, r_S_lq, r_Q_rq, r_S_rq)

      DO i_eqn = 1, i_nprogvars
        DO i_face = 1, i_faceunknowns
          !$OMP ATOMIC
          r_flux(i_eqn,i_ledofsl(i_face)) = r_flux(i_eqn,i_ledofsl(i_face)) + &
            r_edgelen/r_elmtvolume(i_elmt_l)*r_rhs_l(i_eqn,i_face)

!--- if we have an inner or periodic edge, also update right element
          IF(i_edgeboundary(i_edge) >= 0) THEN
            !$OMP ATOMIC
            r_flux(i_eqn,i_ledofsr(i_face)) = r_flux(i_eqn,i_ledofsr(i_face)) - &
              r_edgelen/r_elmtvolume(i_elmt_r)*r_rhs_r(i_eqn,i_face)
          END IF
        END DO
      END DO

    END IF
    END DO edge_loop
    !$OMP END PARALLEL DO

  END SUBROUTINE fvm_flux

!*******************************************************************************
! DESCRIPTION of [SUBROUTINE elmtflux]:
!> @brief computes the divergence and source terms of the equation in consistent form
!
!> @param[out]    r_rhs           DIMENSION(i_nprogvars);
!>                                sum of divergence and source term per quadrature point
!> @param[in]     i_faceunknowns  number of DOFs per element
!> @param[in]     i_equadpts      number of quadrature points per element
!> @param[in]     r_eqwei         weights for the quadrature rule corresponding to the element quadrature points; DIMENSION(i_equadpts)
!> @param[in]     r_eMinvpsi      DIMENSION(i_faceunknowns, i_equadpts);
!>                                inverse mass matrix interpolated on element quadrature points
!> @param[in]     r_Q             state vector (fluid depth, momentum in x-, momentum in y-direction) at element quadrature points; DIMENSION(i_nprogvars, i_equadpts)
!> @param[in]     r_dQdx          derivative of state vector in x-direction
!> @param[in]     r_dQdy          derivative of state vector in y-direction
!
  SUBROUTINE elmtflux(r_rhs, i_faceunknowns, i_equadpts, r_eqwei, r_eMinvpsi, &
                      r_Q, r_dQdx, r_dQdy)

    IMPLICIT NONE

    REAL (KIND = GRID_SR), DIMENSION(:,:), INTENT(out)      :: r_rhs
    INTEGER (KIND = GRID_SI), INTENT(in)                    :: i_faceunknowns, i_equadpts
    REAL (KIND = GRID_SR), DIMENSION(:), INTENT(in)         :: r_eqwei
    REAL (KIND = GRID_SR), DIMENSION(:,:), INTENT(in)       :: r_Q, r_dQdx, r_dQdy, r_eMinvpsi

    INTEGER (KIND = GRID_SI)                                :: i_quad, i_dof
    REAL (KIND = GRID_SR), DIMENSION(i_nprogvars)           :: r_Div

    r_rhs(:,:) = 0.0_GRID_SR

!--- quadrature loop
    elmt_quad_loop: DO i_quad=1, i_equadpts

!--- compute flux divergence
      r_Div = divflux(r_Q(:,i_quad), r_dQdx(:,i_quad), r_dQdy(:,i_quad))

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
!> @param[out]    r_rhs_l         edge flux \f$ (F - F^{\ast})\cdot n \f$ for left element; DIMENSION(3,i_faceunknowns)
!> @param[out]    r_rhs_r         edge flux \f$ (F - F^{\ast})\cdot m \f$ for right element; DIMENSION(3,i_faceunknowns)
!> @param[in]     r_normal        DIMENSION(2);
!>                                local normal vector pointing from left to right element
!> @param[in]     i_faceunknowns  number of DOFs per element
!> @param[in]     i_gquadpts      number of quadrature points per edge
!> @param[in]     r_gqwei         weights for the quadrature rule corresponding to the edge quadrature points; DIMENSION(i_gquadpts)
!> @param[in]     r_gMinvpsi      DIMENSION(i_faceunknowns, i_equadpts);
!>                                inverse mass matrix interpolated on edge quadrature points
!> @param[in]     r_Q_l           state vector (fluid depth, momentum in x-, momentum in y-direction) at edge quadrature points of left element; DIMENSION(i_nprogvars, i_gquadpts)
!> @param[in]     r_Q_r           state vector (fluid depth, momentum in x-, momentum in y-direction) at edge quadrature points of right element; DIMENSION(i_nprogvars, i_gquadpts)
!> @param[in]     r_S_l           source terms (bathymetry, Coriolis, wind) at edge quadrature points of left element; DIMENSION(i_nsrcterms, i_gquadpts)
!> @param[in]     r_S_r           source terms (bathymetry, Coriolis, wind) at edge quadrature points of right element; DIMENSION(i_nsrcterms, i_gquadpts)
!
  SUBROUTINE edgeflux(r_rhs_l, r_rhs_r, r_normal, &
                      i_faceunknowns, i_gquadpts, r_gqwei, r_gMinvpsi, &
                      r_Q_l, r_S_l, r_Q_r, r_S_r)

    IMPLICIT NONE

    REAL (KIND = GRID_SR), DIMENSION(:,:), INTENT(out)      :: r_rhs_l, r_rhs_r
    REAL (KIND = GRID_SR), DIMENSION(2), INTENT(in)         :: r_normal
    INTEGER (KIND = GRID_SI), INTENT(in)                    :: i_faceunknowns, i_gquadpts
    REAL (KIND = GRID_SR), DIMENSION(:), INTENT(in)         :: r_gqwei
    REAL (KIND = GRID_SR), DIMENSION(:,:), INTENT(in)       :: r_Q_l, r_Q_r, r_S_l, r_S_r, r_gMinvpsi

!--- local declarations
    INTEGER (KIND = GRID_SI)                                :: i_quad, i_dof
    REAL (KIND = GRID_SR), DIMENSION(i_nprogvars,2)         :: r_F_l, r_F_r
    REAL (KIND = GRID_SR), DIMENSION(i_nprogvars)           :: r_Fstar, r_Fleft, r_Frght

    r_rhs_l(:,:) = 0.0_GRID_SR
    r_rhs_r(:,:) = 0.0_GRID_SR

!--- perform edge quadrature
    edge_quad_loop: DO i_quad=1,i_gquadpts

      r_F_l = flux(r_Q_l(:,i_quad))
      r_F_r = flux(r_Q_r(:,i_quad))

      r_Fleft = MATMUL(r_F_l, r_normal)
      r_Frght = MATMUL(r_F_r, r_normal)
      r_Fstar = riemannsolver(r_Q_l(:,i_quad), r_Q_r(:,i_quad), r_normal)

      edge_dof_loop: DO i_dof=1,i_faceunknowns
        r_rhs_l(:,i_dof) = r_rhs_l(:,i_dof) - r_gqwei(i_quad)*r_gMinvpsi(i_dof,i_quad)*(r_Fstar-r_Fleft)
        r_rhs_r(:,i_dof) = r_rhs_r(:,i_dof) - r_gqwei(i_quad)*r_gMinvpsi(i_dof,i_quad)*(r_Fstar-r_Frght)
      END DO edge_dof_loop

    END DO edge_quad_loop

  END SUBROUTINE edgeflux

!*******************************************************************************
END MODULE DG_flux
