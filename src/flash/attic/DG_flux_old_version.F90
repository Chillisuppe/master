!*******************************************************************************
! MODULE NAME:
!       DG_flux
! FUNCTION:
!       Computes the strong form of flux for DG method.
!
! VERSION(s):
! 1. version            nicole beisiegel        02/2012
!*******************************************************************************
MODULE DG_flux

  USE GRID_api
  USE FLASH_parameters
  USE DG_boundary
  USE DG_riemann_solver
  USE DG_utils
  USE DG_fluxcorrection

  PRIVATE :: source
  PUBLIC  :: fvm_flux

  CONTAINS

!*******************************************************************************
  SUBROUTINE fvm_flux(r_phi_flux, r_phi, r_bathy, r_taux, r_tauy, r_coriol, r_gamma, r_gammat, r_visc, &
                      r_rho, i_numedge, i_numelmt, i_faceunknowns, i_equadpts, i_gquadpts, r_normals, &
                      r_metrics_inv, r_gpsi, r_epsi, r_eMinvpsi, r_Dxi, r_Deta, r_eqwei, &
                      r_element_vol, r_edgelength, i_edgeinfo, i_index, &
                      i_edgeboundary, r_gqwei, r_gMinvpsi)

    IMPLICIT NONE

    INTEGER (KIND = GRID_SI)                                ::  i_numedge, i_edge, &
                                                                i_elmt, i_elmt_k, i_numelmt, &
                                                                i_faceunknowns, jj, n, &
                                                                i_equadpts, i_gquadpts, &
                                                                i_count, i_quad, kj, i_cnt
    INTEGER (KIND = GRID_SI), DIMENSION(:)                  ::  i_edgeboundary
    INTEGER (KIND = GRID_SI), DIMENSION(:,:)                ::  i_edgeinfo, i_index
    REAL (KIND = GRID_SR)                                   ::  r_gamma, r_gammat, r_visc, r_rho, &
                                                                r_edge_length, r_absv, &
                                                                r_element_volume, r_element_vol_k
    REAL (KIND = GRID_SR), DIMENSION(GRID_DIMENSION)        ::  r_normal
    REAL (KIND = GRID_SR), DIMENSION(:), ALLOCATABLE        ::  r_Hq_l, r_Hq_l2, r_Hq_r, r_Hq_r2, r_Hq_e, r_source_e, &
                                                                r_ssh_k, r_ssh_e, &
                                                                r_height_k, r_height_e, &
                                                                r_bathy_e, r_bathy_k, &
                                                                r_u_k, r_u_e, r_hu_k, r_hu_e, &
                                                                r_v_e, r_v_k, r_hv_k, r_hv_e, &
                                                                r_taux_e, r_tauy_e, r_coriol_e
    REAL (KIND = GRID_SR), DIMENSION(3)                     ::  r_rhs, r_rhs_r, r_flux_e, r_flux_k, &
                                                                r_Fstar, r_Sq, &
                                                                r_force, r_Hq_l_sum !, r_Fstrong
#ifdef FLUXC
    REAL (KIND = GRID_SR), DIMENSION(i_faceunknowns)        ::  r_ssh_star_e, r_ssh_star_k
    REAL (KIND = GRID_SR), DIMENSION(3)                     ::  r_Fstrong_e, r_Fstrong_k
#endif
    REAL (KIND = GRID_SR), DIMENSION(2)                     ::  r_Div
    REAL (KIND = GRID_SR), DIMENSION(:)                     ::  r_edgelength, r_element_vol, r_eqwei, r_gqwei
    REAL (KIND = GRID_SR), DIMENSION(:,:)                   ::  r_taux, r_tauy, &
                                                                r_epsi, r_gpsi, r_eMinvpsi, r_gMinvpsi, r_Dxi, r_Deta
    REAL (KIND = GRID_SR), DIMENSION(3,2)                   ::  r_F_e, r_F_k !, r_flux_corr
    REAL (KIND = GRID_SR), DIMENSION(:,:), ALLOCATABLE      ::  r_Jacobi_e, r_Jacobi_k
    REAL (KIND = GRID_SR), DIMENSION(:,:), ALLOCATABLE      ::  r_N, r_N_r, r_Dxi_in, r_Deta_in
    REAL (KIND = GRID_SR), DIMENSION(:,:)                   ::  r_normals, r_bathy, r_coriol
    REAL (KIND = GRID_SR), DIMENSION(:,:,:), INTENT(inout)  ::  r_phi_flux
    REAL (KIND = GRID_SR), DIMENSION(:,:,:)                 ::  r_phi, r_metrics_inv
    INTEGER (KIND = GRID_SI), DIMENSION(:), ALLOCATABLE     ::  i_ind_r

!> Allocate workspace

    ALLOCATE( r_source_e(3), r_Hq_l(3), r_Hq_l2(3), r_Hq_r(3), r_Hq_r2(3), r_Hq_e(3), &
              r_Jacobi_e(2,2),r_Jacobi_k(2,2), r_ssh_k(i_faceunknowns), r_ssh_e(i_faceunknowns), &
              r_height_k(i_faceunknowns), r_height_e(i_faceunknowns), r_bathy_e(i_faceunknowns), r_bathy_k(i_faceunknowns), &
              r_u_k(i_faceunknowns), r_u_e(i_faceunknowns), r_v_e(i_faceunknowns), r_v_k(i_faceunknowns), &
              r_hu_k(i_faceunknowns), r_hu_e(i_faceunknowns), r_hv_e(i_faceunknowns), r_hv_k(i_faceunknowns), &
              r_taux_e(i_faceunknowns), r_tauy_e(i_faceunknowns), r_coriol_e(i_faceunknowns), &
              r_N(i_faceunknowns, i_faceunknowns), r_N_r(i_faceunknowns, i_faceunknowns), &
              i_ind_r(i_faceunknowns), r_Dxi_in(i_faceunknowns, i_faceunknowns), &
              r_Deta_in(i_faceunknowns, i_faceunknowns))

!> Initialize local quantities
    r_phi_flux = 0._GRID_SR
    r_Hq_l_sum = 0._GRID_SR

!> Compute interpolated D-Matrix
    r_Dxi_in  = matmul(r_Dxi, r_epsi)  !int dhde
    r_Deta_in = matmul(r_Deta,r_epsi)  !int bei frank dhdn 1/4

    r_Dxi_in  = matmul(r_Dxi_in, r_eMinvpsi)
    r_Deta_in = matmul(r_Deta_in, r_eMinvpsi)

    element_loop: DO i_elmt = 1,i_numelmt

      r_Jacobi_e = 0._GRID_SR
      r_Jacobi_k = 0._GRID_SR

! Store local variables for this element
      r_element_volume = r_element_vol(i_elmt)
      r_bathy_e        = r_bathy(i_elmt,:)
      r_ssh_e          = r_phi(i_elmt,:,1)
      r_height_e       = r_ssh_e  + r_bathy_e
      r_hu_e           = r_phi(i_elmt,:,2)
      r_hv_e           = r_phi(i_elmt,:,3)
      DO i_cnt = 1, i_faceunknowns
        r_absv = heaviside(r_height_e(i_cnt))
        r_u_e(i_cnt)     = r_absv*r_hu_e(i_cnt) + (1-r_absv)*(r_hu_e(i_cnt)/(r_height_e(i_cnt)+10E-20))
        r_v_e(i_cnt)     = r_absv*r_hv_e(i_cnt) + (1-r_absv)*(r_hv_e(i_cnt)/(r_height_e(i_cnt)+10E-20))
      END DO
      r_coriol_e       = r_coriol(i_elmt,:)
      r_taux_e         = r_taux(i_elmt,:)
      r_tauy_e         = r_tauy(i_elmt,:)
      r_Jacobi_e       = r_metrics_inv(:,:,i_elmt)

!> Face-Loop
      face_loop: DO i_count=1,i_faceunknowns
        r_Hq_r = 0._GRID_SR

!> Quadrature loop
        quad_loop: DO i_quad=1, i_faceunknowns

!> Initialize local vectors
          r_Hq_e = 0._GRID_SR
          r_F_e  = 0._GRID_SR
          r_Div  = 0._GRID_SR

!> Evaluate flux at quadrature points
          r_F_e  = flux(r_ssh_e(i_quad), r_bathy_e(i_quad), r_hu_e(i_quad), r_hv_e(i_quad))

!> Compute divergence of flux tensor
          r_Div  = matmul(r_Jacobi_e,[r_Dxi_in(i_count, i_quad), r_Deta_in(i_count,i_quad)] )
          r_Hq_e = matmul(r_F_e, r_Div)

          r_source_e = source(r_gamma, r_visc, r_gammat, r_rho, r_coriol_e, r_bathy_e, r_ssh_e, &
                              r_u_e, r_v_e, r_hu_e, r_hv_e, r_taux_e, r_tauy_e, r_Jacobi_e, &
                              i_quad, i_elmt, i_faceunknowns, r_Dxi, r_Deta, r_epsi, r_eMinvpsi)  !statt Jacobi_e

          r_Hq_r  = r_Hq_r  + (r_Hq_e+r_source_e)*r_eqwei(i_quad) !*r_eMinvpsi(i_count,i_quad)

        END DO quad_loop !i_quad

!> Add source and div F to flux vector
        r_phi_flux(i_elmt,i_count,:) = r_phi_flux(i_elmt,i_count,:) - r_Hq_r

      END DO face_loop
    END DO element_loop

!> Separate edge-loop

    edge_loop: DO i_edge= 1, i_numedge
      i_elmt   = i_edgeinfo(3,i_edge)
      i_elmt_k = i_edgeinfo(4,i_edge)

      r_element_volume = r_element_vol(i_elmt)
      r_bathy_e        = r_bathy(i_elmt,:)
      r_ssh_e          = r_phi(i_elmt,:,1)
      r_height_e       = r_ssh_e  + r_bathy_e
      r_hu_e           = r_phi(i_elmt,:,2)
      r_hv_e           = r_phi(i_elmt,:,3)
      DO i_cnt = 1, i_faceunknowns
        r_absv = heaviside(r_height_e(i_cnt))
        r_u_e(i_cnt)     = r_absv*r_hu_e(i_cnt) + (1-r_absv)*(r_hu_e(i_cnt)/(r_height_e(i_cnt)+10E-20))
        r_v_e(i_cnt)     = r_absv*r_hv_e(i_cnt) + (1-r_absv)*(r_hv_e(i_cnt)/(r_height_e(i_cnt)+10E-20))
      END DO
      r_coriol_e = r_coriol(i_elmt,:)
      r_taux_e   = r_taux(i_elmt,:)
      r_tauy_e   = r_tauy(i_elmt,:)
      r_Jacobi_e = r_metrics_inv(:,:,i_elmt)

      r_normal = r_normals(:,i_edge)
      n = 1
      r_edge_length = r_edgelength(i_edge)

! Further bookkeeping
      i_ind_r(:) = i_index(i_edge,:)

! Compute values for the second element
      CALL DG_bou(i_elmt, i_elmt_k, i_faceunknowns, i_edgeboundary(i_edge), r_element_vol_k, r_element_vol, &
                  r_ssh_e, r_ssh_k, r_u_e, r_v_e, r_u_k, r_v_k, r_hu_e, r_hv_e, r_hu_k, r_hv_k, r_bathy_e, r_bathy_k, &
                  r_height_k, r_bathy(i_elmt_k,:), r_normal, r_phi(i_elmt_k,:,:), r_Jacobi_k, r_metrics_inv(:,:,i_elmt_k))

! Depending on the edgenodes we choose the boundary matrix N^k
      SELECT CASE(i_edgeinfo(n,i_edge))
        CASE(1)
          r_N = GRID_femtypes%p_type(FEM_ZETA)%sig%r_bound_ma1
        CASE(2)
          r_N = GRID_femtypes%p_type(FEM_ZETA)%sig%r_bound_ma2
        CASE(3)
          r_N = GRID_femtypes%p_type(FEM_ZETA)%sig%r_bound_ma3
      END SELECT

      SELECT CASE(i_edgeinfo(2,i_edge))
        CASE(1)
          r_N_r = GRID_femtypes%p_type(FEM_ZETA)%sig%r_bound_ma1
        CASE(2)
          r_N_r = GRID_femtypes%p_type(FEM_ZETA)%sig%r_bound_ma2
        CASE(3)
          r_N_r = GRID_femtypes%p_type(FEM_ZETA)%sig%r_bound_ma3
      END SELECT

#ifdef FLUXC
!> Fix for well-balancedness acc. to Xing, Zhang & Shu, 2011 (just differs from ssh for discontinuous bathymetry)
      DO jj=1, i_faceunknowns
        r_ssh_star_e(jj) = DMAX1(0._GRID_SR, r_ssh_e(jj) + r_bathy_e(jj) - &
                                             max(r_bathy_e(jj), r_bathy_k(i_ind_r(jj))))
        r_ssh_star_k(jj) = DMAX1(0._GRID_SR, r_ssh_k(i_ind_r(jj)) + r_bathy_k(i_ind_r(jj)) - &
                                             max(r_bathy_e(jj), r_bathy_k(i_ind_r(jj))))
      END DO !jj
#endif

!> Perform edge quadrature
      faceloop2: DO i_count = 1, i_faceunknowns
        r_Hq_l  = 0._GRID_SR
        r_Hq_l2 = 0._GRID_SR

        innerloop4: DO kj=1,i_faceunknowns

          r_rhs   = 0._GRID_SR
          r_rhs_r = 0._GRID_SR
          r_Fstar = 0._GRID_SR

#ifdef FLUXC
        CALL flux_correction(r_ssh_star_e(i_quad), r_ssh_star_k(i_ind_r(i_quad)), &
                             r_bathy_e(i_quad), r_bathy_k(i_ind_r(i_quad)), &
                             r_u_e(i_quad), r_u_k(i_ind_r(i_quad)), r_v_e(i_quad), r_v_k(i_ind_r(i_quad)), &
                             r_hu_e(i_quad), r_hu_k(i_ind_r(i_quad)), r_hv_e(i_quad), r_hv_k(i_ind_r(i_quad)), &
                             r_height_e(i_quad), r_height_k(i_ind_r(i_quad)), r_ssh_e(i_quad), r_ssh_k(i_ind_r(i_quad)), &
                             r_normal, r_Fstrong_e, r_Fstrong_k)
#endif

          r_F_e    = flux(r_ssh_e(kj)         , r_bathy_e(kj)         , r_hu_e(kj)         , r_hv_e(kj)         )
          r_F_k    = flux(r_ssh_k(i_ind_r(kj)), r_bathy_k(i_ind_r(kj)), r_hu_k(i_ind_r(kj)), r_hv_k(i_ind_r(kj)))
          r_flux_e = matmul(r_F_e, r_normal)
          r_flux_k = matmul(r_F_k, r_normal)

          r_Fstar = riemannsolver(r_ssh_e(kj), r_ssh_k(i_ind_r(kj)), &
                                  r_bathy_e(kj), r_bathy_k(i_ind_r(kj)), &
                                  r_u_e(kj), r_u_k(i_ind_r(kj)), r_v_e(kj), r_v_k(i_ind_r(kj)), &
                                  r_hu_e(kj), r_hu_k(i_ind_r(kj)), r_hv_e(kj), r_hv_k(i_ind_r(kj)), &
                                  r_normal)

          r_rhs   = r_N(i_count, kj) * (r_flux_e - r_Fstar)
          r_rhs_r = r_N(i_count, kj) * (r_flux_k - r_Fstar)

          r_Hq_l  = r_Hq_l  + r_rhs  *r_edge_length
          r_Hq_l2 = r_Hq_l2 + r_rhs_r*r_edge_length

        END DO innerloop4   !kj

        r_phi_flux(i_elmt  ,i_count,:) = r_phi_flux(i_elmt,i_count,:) +r_Hq_l/r_element_volume

!> If edge is reflecting theres no update for the right element TO DO remove if
        IF(i_edgeboundary(i_edge) .EQ. 0) r_phi_flux(i_elmt_k,i_ind_r(i_count),:) = &
          r_phi_flux(i_elmt_k, i_ind_r(i_count),:)-r_Hq_l2/r_element_vol_k
      END DO faceloop2 !i_count

    END DO edge_loop

!------ Deallocate workspace
    DEALLOCATE(r_source_e, r_Hq_l, r_Hq_r, r_Hq_e, r_Jacobi_e, &
               r_Jacobi_k, r_ssh_k, r_ssh_e, r_height_k, r_height_e, r_bathy_e, &
               r_bathy_k, r_u_k, r_u_e, r_v_e, r_v_k, &
               r_taux_e, r_tauy_e, r_coriol_e, r_N, i_ind_r, r_hu_e, r_hu_k, r_hv_e, r_hv_k, &
               r_Dxi_in, r_Deta_in, r_N_r)

  END SUBROUTINE fvm_flux

!*******************************************************************************
  FUNCTION source(r_gamma, r_visc, r_gammat, r_rho, r_coriol, r_bathy, r_ssh, r_u, r_v, r_hu, r_hv, r_taux, r_tauy, &
                  r_metrics, i_face, i_elmt, i_faceunknowns, r_Dxi, r_Deta, r_epsi, r_eMinvpsi)

    IMPLICIT NONE

    INTEGER (KIND = GRID_SI)                            :: i_face, i_elmt
    REAL (KIND = GRID_SR)                               :: r_gamma, r_gammat, r_rho, &
                                                           r_manning, r_l2velo, r_root, &
                                                           r_absv, r_visc
    REAL (KIND = GRID_SR), DIMENSION(1,2)               :: r_grad_b, r_grad_ssh
    REAL (KIND = GRID_SR), DIMENSION(:), ALLOCATABLE    :: r_ssh, r_bathy, r_height, &
                                                           r_u, r_v, r_hv, r_hu, &
                                                           r_coriol, r_taux, r_tauy, &
                                                           source, r_dbdxi, r_dbdeta, &
                                                           r_dhdxi, r_dhdeta
    REAL (KIND=GRID_SR), DIMENSION(i_faceunknowns)      :: r_ssh_in
    REAL (KIND=GRID_SR), DIMENSION(2,2)                 :: r_laplace_hu, r_laplace_hv
    REAL(KIND=GRID_SR), DIMENSION(:,:), ALLOCATABLE     :: r_Dxi_in, r_Deta_in, r_PsiM
    REAL(KIND=GRID_SR), DIMENSION(:,:)                  :: r_epsi, r_Dxi, r_Deta, r_eMinvpsi
    REAL (KIND = GRID_SR), DIMENSION(2,2)               :: r_metrics, r_metrics_transp
    INTEGER (KIND = GRID_SI)                            :: i_faceunknowns

!------- Allocate workspace
    ALLOCATE(r_height(i_faceunknowns), source(3), r_dbdxi(i_faceunknowns), r_dbdeta(i_faceunknowns), &
             r_Dxi_in(i_faceunknowns, i_faceunknowns), &
             r_Deta_in(i_faceunknowns, i_faceunknowns), &
             r_PsiM(i_faceunknowns, i_faceunknowns))

    r_height = r_ssh + r_bathy

    r_Dxi_in  = matmul(r_Dxi, r_epsi)  !int dhde
    r_Deta_in = matmul(r_Deta,r_epsi)  !int bei frank dhdn

    r_PsiM    = matmul(r_eMinvpsi, r_epsi)

    r_Dxi_in  = matmul(r_eMinvpsi, r_Dxi_in)
    r_Deta_in = matmul(r_eMinvpsi, r_Deta_in) !TODO nochmal checken für grad bathy

    r_dbdxi   = matmul(r_Dxi_in, r_bathy)
    r_dbdeta  = matmul(r_Deta_in, r_bathy)

    r_dhdxi   = matmul(r_Dxi_in, r_ssh)
    r_dhdeta  = matmul(r_Deta_in, r_ssh)

    r_grad_b(1,1:2)=matmul(r_metrics,[r_dbdxi(i_face),r_dbdeta(i_face)])
    r_grad_ssh(1,1:2)=matmul(r_metrics,[r_dhdxi(i_face),r_dhdeta(i_face)])

    r_laplace_hu = laplace(r_hu, r_Dxi_in, r_Deta_in, r_metrics, i_faceunknowns, i_face)
    r_laplace_hv = laplace(r_hv, r_Dxi_in, r_Deta_in, r_metrics, i_faceunknowns, i_face)
    r_ssh_in = matmul(r_PsiM, r_ssh)

!> Compute Manning roughness coefficient - r_gamma is Manning coefficient
!> Assumption: Friction is zero is everything is dry.
     r_root     = r_height(i_face)**(4/3)
     r_l2velo   = sqrt(r_u(i_face)**2+r_v(i_face)**2)
     r_absv     = heaviside(r_height(i_face))
     r_manning  = 0._GRID_SR + (1-r_absv)*(r_l2velo/(r_root+10E-20)*r_gamma**2)

!------- compute the source term

    source(1) = 0._GRID_SR
    source(2) = 0._GRID_SR
    source(3) = 0._GRID_SR

#ifdef INHOMOGEN
    source(2) = -(-r_hv(i_face) * r_coriol(i_face) &
                + r_ssh_in(i_face)*r_grad_b(1,1) &
                + (1-r_absv)*(r_taux(i_face)*r_gammat /(r_rho*(r_height(i_face)+10E-20))) & !height(i_face)) &
                - r_manning*r_hu(i_face) &
                )! + r_visc*(r_laplace_hu(1,1)+r_laplace_hv(1,1)))

    source(3) = -(-r_hu(i_face) * r_coriol(i_face) &
                + r_ssh_in(i_face)*r_grad_b(1,2) &
                + (1-r_absv)*(r_tauy(i_face)*r_gammat /(r_rho*(r_height(i_face)+10E-20))) & !height(i_face)) &
                - r_manning*r_hv(i_face) &
                )!+ r_visc*(r_laplace_hu(2,2)+r_laplace_hv(2,2)))
#endif

!------- Deallocate workspace
    DEALLOCATE(r_height, source, r_dbdxi, r_dbdeta, r_Dxi_in, r_Deta_in, r_PsiM)

  END FUNCTION source

!*******************************************************************************
END MODULE DG_flux
