!*************************************************************
!
!> @file DG_limiter_nodal.F90
!> @brief contains module DG_limiter
!
!*************************************************************
MODULE DG_limiter

  USE GRID_api
  USE FLASH_parameters
  USE DG_equation
#ifdef BERNSTEIN
  USE DG_HO_utils
#endif
  PRIVATE  :: assemble_taylor, compute_alpha, &
              overbar
  PUBLIC   :: limiter

  CONTAINS

!*******************************************************************************
! DESCRIPTION of [SUBROUTINE limiter]:
!> @brief flux limiter for DG method
!>
!> @param[in]       p_ghand             grid handling data structure
!> @param[in,out]   r_Q                 discrete solution to be limited;
!>                                      DIMENSION(3,i_numelmt*i_faceunknowns)
!> @param[in]       i_elementedges      edge indices of the edges of each element; DIMENSION(3,i_numelmt)
!> @param[in]       i_edgeinfo          edge-element relation for each edge: local edge index and global element index for adjacent elements; DIMENSION(4,i_numedge)
!> @param[in]       r_edgelength        edge length; DIMENSION(i_numedge)
!> @param[in]       r_metrics_inv       metric terms for each element; DIMENSION(2,2,i_numelmt)
!> @param[in]       i_numelmt           total number of elements
!> @param[in]       r_coonod            node coordinates
!> @param[in]       r_coodof            DOF coordinates
!> @param[in]       i_edgenodes         global node indices for each edge; DIMENSION(2,i_numedge)
!> @param[in]       i_elementnodes      node indices of each element; DIMENSION(3, i_numelmt)
!> @param[in]       i_elementdofs       DOF indices corresponding to each element; DIMENSION(i_faceunknowns, i_numelmt)
!> @param[in]       i_numedge           total number of edges
!> @param[in]       r_S                 source terms at DOFs; DIMENSION(i_nsrcterms, i_numelmt*i_faceunknowns)
!> @param[in]       i_faceunknowns      number of DOFs per element
!> @param[in]       r_evander
!> @param[in]       r_einvvander
!

!**************************************************************
  SUBROUTINE limiter(p_ghand, r_QU, i_elementedges, i_edgeinfo, r_edgelength, r_metricsinv, &
                     i_numelmt, r_coonod, r_coodof, i_edgenodes, i_elementnodes, i_elementdofs, &
                     i_numedge, r_S, i_faceunknowns, r_evander, r_einvvander)

    IMPLICIT NONE

    TYPE(grid_handle)                                      :: p_ghand
    REAL(KIND=GRID_SR), DIMENSION(:,:), INTENT(inout)      :: r_QU
    REAL(KIND=GRID_SR), DIMENSION(:), INTENT(in)           :: r_edgelength
    REAL(KIND=GRID_SR), DIMENSION(:,:), INTENT(in)         :: r_coodof, r_coonod, r_evander, r_einvvander, r_S
    REAL(KIND=GRID_SR), DIMENSION(:,:,:), INTENT(in)       :: r_metricsinv
    INTEGER(KIND=GRID_SI), DIMENSION(:,:), INTENT(in)      :: i_edgeinfo, i_elementedges, i_edgenodes, &
                                                              i_elementnodes, i_elementdofs
    INTEGER(KIND=GRID_SI), INTENT(in)                      :: i_faceunknowns, i_numelmt, &
                                                              i_numedge

!-- Local declarations
    INTEGER(KIND=GRID_SI)  				   :: i_elmt, i_numnode, i_cnt, i_dof,  &
    							      i_iswet, i_node, i_elmt_p, i_patch

    INTEGER(KIND=GRID_SI), DIMENSION(:), ALLOCATABLE      :: i_ndmaster
    INTEGER(KIND=GRID_SI), DIMENSION(:,:), ALLOCATABLE    :: i_nodepatch

    REAL(KIND=GRID_SR)                                     :: r_alpha_x, r_alpha_y, r_alpha2_x, r_alpha2_y
    REAL(KIND=GRID_SR)                                     :: r_q_mean, r_hv_mean, r_hu_mean
    REAL(KIND=GRID_SR)                                     :: r_patch_mean, r_patch_mean2
    REAL(KIND=GRID_SR)                                     :: r_patch_mean_gradx, r_patch_mean_grady
    REAL(KIND=GRID_SR)                                     :: r_patch_mean2_gradx, r_patch_mean2_grady
    REAL(KIND=GRID_SR)                                     :: r_wettol = 10E-3


    REAL(KIND=GRID_SR), DIMENSION(2)                       :: r_b, r_alpha, r_alpha2
    REAL(KIND=GRID_SR), DIMENSION(2)                       :: r_q_meangrad, r_hu_meangrad, r_hv_meangrad

    REAL(KIND=GRID_SR), DIMENSION(3)                       :: r_patch_min, r_patch_max
    REAL(KIND=GRID_SR), DIMENSION(3)                       :: r_patch_min2, r_patch_max2

    REAL(KIND=GRID_SR), DIMENSION(i_faceunknowns)          :: r_bat, r_coox, r_cooy, r_omega
    REAL(KIND=GRID_SR), DIMENSION(i_faceunknowns)          :: r_q, r_qqc, r_q_corr, r_zero, r_qval, r_qval2, r_qq
    REAL(KIND=GRID_SR), DIMENSION(i_faceunknowns)          :: r_h, r_hu, r_hv, r_qhu, r_qhv
    REAL(KIND=GRID_SR), DIMENSION(i_faceunknowns)          :: r_val_corr
    REAL(KIND=GRID_SR), DIMENSION(i_faceunknowns)          :: r_q_grad_x, r_q_grad_y
    REAL(KIND=GRID_SR), DIMENSION(i_faceunknowns)          :: r_hu_grad_x, r_hu_grad_y
    REAL(KIND=GRID_SR), DIMENSION(i_faceunknowns)          :: r_hv_grad_x, r_hv_grad_y
    REAL(KIND=GRID_SR), DIMENSION(i_faceunknowns)          :: r_tmp_xx, r_tmp_xy, r_tmp_yy, r_tmp_yx
    REAL(KIND=GRID_SR), DIMENSION(i_faceunknowns)          :: r_tmp2_xx, r_tmp2_xy, r_tmp2_yy, r_tmp2_yx
    REAL(KIND=GRID_SR), DIMENSION(i_faceunknowns)          :: r_patch_val, r_patch_val2
    REAL(KIND=GRID_SR), DIMENSION(i_faceunknowns)          :: r_q_lim, r_hu_lim, r_hv_lim
    REAL(KIND=GRID_SR), DIMENSION(i_faceunknowns)          :: r_Psi_c

    REAL(KIND=GRID_SR), DIMENSION(2,2)			   :: r_metrics_inv
    REAL(KIND=GRID_SR), DIMENSION(2,2)                     :: r_q_meanlaplace, r_hu_meanlaplace, r_hv_meanlaplace

    REAL(KIND=GRID_SR), DIMENSION(i_faceunknowns, i_faceunknowns)	:: r_epsi, r_Dxi, r_Deta, r_ddx, r_ddy

    REAL(KIND=GRID_SR), DIMENSION(i_numelmt*i_faceunknowns)             :: r_QU_lim


!-- Allocate workspace
    i_numnode = p_ghand%i_nnumber

    ALLOCATE(i_ndmaster(i_numnode), &
             i_nodepatch(GRID_patchelements, i_numnode))

!-- Get nodepatch info via grid_getinfo
    CALL grid_getinfo(p_ghand, l_relative=.TRUE., &
                      l_finelevel=.TRUE., i_femtype=FEM_DG, i_nodepatch=i_nodepatch, i_nodemaster=i_ndmaster)

!-- Initialize matrices and arrays

    r_Dxi  = GRID_femtypes%p_type(FEM_DG)%sig%r_dpsidxi
    r_Deta = GRID_femtypes%p_type(FEM_DG)%sig%r_dpsideta
    r_epsi = GRID_femtypes%p_type(FEM_DG)%sig%r_epsiquad
    r_omega = GRID_femtypes%p_type(FEM_DG)%sig%r_equadwei

!-- Interpolation on center points.
 IF (i_faceunknowns .EQ. 3_GRID_SI) THEN
     r_Psi_c = (/ 0.333333333333333333_GRID_SR, &
          0.333333333333333333_GRID_SR, &
          0.333333333333333333_GRID_SR /)
  ELSE
     !-- only for quadratic functions
     r_Psi_c = (/ -0.11111111111111_GRID_SR, &
          -0.11111111111111_GRID_SR, &
          -0.11111111111111_GRID_SR, &
          0.44444444444444_GRID_SR, &
          0.44444444444444_GRID_SR, &
          0.4444444444444_GRID_SR /)

     !-- Take different values for Bernstein polynomials
#ifdef BERNSTEIN
     r_Psi_c = (/ 0.1111111111112_GRID_SR, &
          0.1111111111111_GRID_SR, &
          0.1111111111111_GRID_SR, &
          0.2222222222222_GRID_SR, &
          0.2222222222222_GRID_SR, &
          0.2222222222222_GRID_SR /)
#endif
  END IF

!-- STEP 1 APPLY LIMITER TO **FLUID HEIGHT/ TOTAL HEIGHT** --

!-- Do the actual limiting
    DO i_elmt=1, i_numelmt
       r_metrics_inv = r_metricsinv(:,:,i_elmt)

       !-- read in local fluid height
       r_h = r_QU(1,i_elementdofs(:,i_elmt))

      !-- Introduce index that is either 1 or 0 depending on if the elmt is wet or not
       i_iswet = 1
       IF (MINVAL(r_h) .LE. r_wettol) i_iswet = 0

!-- Limit on H for wet and h for semi-wet elements
       r_bat = r_S(1,i_elementdofs(:,i_elmt))
       r_q   = r_h + i_iswet*r_bat

       !-- local coordinates
       r_coox = r_coodof(1,i_elementdofs(:,i_elmt))
       r_cooy = r_coodof(2,i_elementdofs(:,i_elmt))
       !-- Compute barycentre
       r_b(1) = sum(r_coox(1:3))/3.0
       r_b(2) = sum(r_cooy(1:3))/3.0


!-- b) Compute the averages for q
       r_zero   = 0.0_GRID_SR
       r_qq     = MATMUL(r_q, r_epsi)
       !-- Corrected values for quadrature
       r_qqc    = MERGE(r_zero, r_qq, r_qq .LE. r_wettol)
       !-- Corrected values for gradient computation
       r_q_corr = MERGE(r_zero, r_q, r_q .LE. r_wettol)
       r_q_mean = 0.5*dot_product(r_qqc, r_omega)

       !-- Compute operators to perform differentiation up to second order
       DO i_cnt = 1, i_faceunknowns
          r_ddx(:,i_cnt)  = r_Dxi(:,i_cnt)*r_metrics_inv(1,1) + r_Deta(:,i_cnt)*r_metrics_inv(1,2)
          r_ddy(:,i_cnt)  = r_Dxi(:,i_cnt)*r_metrics_inv(2,1) + r_Deta(:,i_cnt)*r_metrics_inv(2,2)
       END DO !i_cnt

       !-- Compute all derivatives
       r_q_grad_x = MATMUL(r_q, r_ddx)
       r_q_grad_y = MATMUL(r_q, r_ddy) !instead of q_corr

       r_q_meangrad(1) = dot_product(r_Psi_c, r_q_grad_x)
       r_q_meangrad(2) = dot_product(r_Psi_c, r_q_grad_y)

       !-- Compute 2nd order derivatives inside the element (at centroid)
       r_tmp_xx = matmul(r_q_grad_x, r_ddx)
       r_tmp_xy = matmul(r_q_grad_x, r_ddy)
       r_tmp_yx = matmul(r_q_grad_y, r_ddx)
       r_tmp_yy = matmul(r_q_grad_y, r_ddy)

       ! Laplace is wrong
       r_q_meanlaplace(1,1) = 0.0 !0.5*dot_product(MATMUL(r_tmp_xx,r_epsi), r_omega)
       r_q_meanlaplace(1,2) = 0.0 !0.5*dot_product(MATMUL(r_tmp_xy,r_epsi), r_omega)
       r_q_meanlaplace(2,1) = 0.0 !0.5*dot_product(MATMUL(r_tmp_yx,r_epsi), r_omega)
       r_q_meanlaplace(2,2) = 0.0 !0.5*dot_product(MATMUL(r_tmp_yy,r_epsi), r_omega)

       IF (r_q_meanlaplace(1,1) /= 0.0 ) write(*,*) "Error Laplace computation xx", r_q_meanlaplace(1,1)
       IF (r_q_meanlaplace(1,2) /= 0.0 ) write(*,*) "Error Laplace computation xy", r_q_meanlaplace(1,2)
       IF (r_q_meanlaplace(2,1) /= 0.0 ) write(*,*) "Error Laplace computation yx", r_q_meanlaplace(2,1)
       IF (r_q_meanlaplace(2,2) /= 0.0 ) write(*,*) "Error Laplace computation yy", r_q_meanlaplace(2,2)


!-- c) Compute Patch Extrema

       !-- Initialize local min/max
       r_patch_min(1) = MIN( 10E+6, r_q_mean)
       r_patch_max(1) = MAX(-10E+6, r_q_mean)

       r_patch_min(2) = MIN( 10E+6, r_q_meangrad(1))
       r_patch_max(2) = MAX(-10E+6, r_q_meangrad(1))

       r_patch_min(3) = MIN( 10E+6, r_q_meangrad(2))
       r_patch_max(3) = MAX(-10E+6, r_q_meangrad(2))

       !-- Compute min/max for
       DO i_node=1,GRID_elementnodes
          DO i_elmt_p=1,GRID_patchelements
             i_patch = i_nodepatch(i_elmt_p, i_ndmaster(i_elementnodes(i_node, i_elmt)))
             IF (i_patch > 0) THEN

                !-- Compute mean value and mean grad of patch element
                r_zero = 0.0
                r_patch_val =  r_QU(1,i_elementdofs(:,i_patch))+i_iswet*r_S(1,i_elementdofs(:,i_patch))
                r_qval      = MATMUL(r_patch_val, r_epsi)
                !r_val_corr  = r_patch_val

                !-- correct for negative entries
                r_qval = MERGE(r_zero, r_qval, r_qval .LE. r_wettol)
                r_val_corr = MERGE(r_zero, r_patch_val, r_patch_val .LE. r_wettol) !val_corr, r_val_corr .LE. r_wettol)

                r_patch_mean = 0.5*dot_product(r_val_corr, r_omega) !qval --> corr

                r_patch_mean_gradx = 0.5*dot_product(MATMUL(MATMUL(r_qval, r_ddx), r_epsi), r_omega)
                r_patch_mean_grady = 0.5*dot_product(MATMUL(MATMUL(r_qval, r_ddy), r_epsi), r_omega) !qval --> corr

                r_patch_min(1) = MIN(r_patch_min(1), r_patch_mean)
                r_patch_max(1) = MAX(r_patch_max(1), r_patch_mean)
                r_patch_min(2) = MIN(r_patch_min(2), r_patch_mean_gradx)
                r_patch_max(2) = MAX(r_patch_max(2), r_patch_mean_gradx)
                r_patch_min(3) = MIN(r_patch_min(3), r_patch_mean_grady)
                r_patch_max(3) = MAX(r_patch_max(3), r_patch_mean_grady)
             END IF
          END DO !i_elmt_p
       END DO !i_node

!-- d) Compute alpha
       r_alpha= 0.0
       r_alpha(1) = compute_alpha(r_q,        r_q_mean,        r_patch_min(1), r_patch_max(1), i_faceunknowns)
       r_alpha_x  = compute_alpha(r_q_grad_x, r_q_meangrad(1), r_patch_min(2), r_patch_max(2), i_faceunknowns)
       r_alpha_y  = compute_alpha(r_q_grad_y, r_q_meangrad(2), r_patch_min(3), r_patch_max(3), i_faceunknowns)
       r_alpha(2) = MIN(r_alpha_x, r_alpha_y)

!-- e) Taylor expansion

5953       r_q_lim = assemble_taylor(r_q_mean, r_q_meangrad(1), r_q_meangrad(2), r_q_meanlaplace, &
                r_epsi, r_alpha, r_b, r_coox, r_cooy, i_faceunknowns)

       IF (MINVAL(r_q_lim) .LT. -10E-8) THEN
          IF((i_faceunknowns .EQ. 6 ) .AND. (r_alpha(2) .GT. 10E-6)) THEN
             r_alpha(2) = r_alpha(2)*0.5
          ELSEIF((i_faceunknowns .EQ. 6 ) .AND. (r_alpha(2) .LT. 10E-6)) THEN
             r_alpha(2) = 0.0
             r_alpha(1) = r_alpha(1)*0.5
          ELSEIF(r_alpha(1) .GT. 10E-6) THEN
             r_alpha(1) = r_alpha(1)*0.5
          ELSE
             r_alpha(1) = 0.0
          END IF
          GO TO 5953
       END IF

!-- f) Assign limited value to vector r_QU
    DO i_dof = 1, i_faceunknowns
       r_QU_lim(i_elementdofs(i_dof,i_elmt)) = MAX(r_q_lim(i_dof) - i_iswet*r_bat(i_dof),0.0)
    END DO !i_dof

    END DO !i_elmt

!-- Until here this should be testable.
!----------------------------------------------------------------!

!-- STEP 2 APPLY LIMITER TO **VELOCITY / MOMENTUM** --
! Do the actual limiting on the velocity
    DO i_elmt=1, i_numelmt
       r_metrics_inv = r_metricsinv(:,:,i_elmt)

       !-- Introduce index that is either 1 or 0 depending on if the elmt is wet or not
       !-- Limit on (hu,hv) for wet and (u,v) for semi-wet elements

       r_h  = r_QU(1, i_elementdofs(:,i_elmt))
       r_hu = r_QU(2, i_elementdofs(:,i_elmt))
       r_hv = r_QU(3, i_elementdofs(:,i_elmt))

       i_iswet = 1
       IF (MINVAL(r_h) .LE. r_wettol) i_iswet = 0

!-- compute to be limited quantities depending on i_iswet
!       DO i_dof=1, i_faceunknowns

          r_hu(:) = (i_iswet)*r_hu(:) + (1-i_iswet) * velocity(r_h(:), r_hu(:))
          r_hv(:) = (i_iswet)*r_hv(:) + (1-i_iswet) * velocity(r_h(:), r_hv(:))

          DO i_dof=1, i_faceunknowns
          IF (r_h(i_dof) .LE. r_wettol) THEN
             r_hu(i_dof) = 0.0
             r_hv(i_dof) = 0.0
          END IF

      END DO !i_dof

       !-- local coordinates
       r_coox = r_coodof(1,i_elementdofs(:,i_elmt))
       r_cooy = r_coodof(2,i_elementdofs(:,i_elmt))
       !-- Compute barycentre
       r_b(1) = sum(r_coox(1:3))/3.0
       r_b(2) = sum(r_cooy(1:3))/3.0

!-- b) Compute the averages for (hu,hv), (u,v) resp.
       r_qhu     = MATMUL(r_hu, r_epsi)
       r_qhv     = MATMUL(r_hv, r_epsi)

       r_hu_mean = 0.5*dot_product(r_qhu, r_omega)
       r_hv_mean = 0.5*dot_product(r_qhv, r_omega)

       !-- Compute operators to perform differentiation up to second order
       DO i_cnt = 1, i_faceunknowns
          r_ddx(:,i_cnt)  = r_Dxi(:,i_cnt)*r_metrics_inv(1,1) + r_Deta(:,i_cnt)*r_metrics_inv(1,2)
          r_ddy(:,i_cnt)  = r_Dxi(:,i_cnt)*r_metrics_inv(2,1) + r_Deta(:,i_cnt)*r_metrics_inv(2,2)
       END DO !i_cnt

       !-- Compute all derivatives
       r_hu_grad_x = MATMUL(r_hu, r_ddx)
       r_hu_grad_y = MATMUL(r_hu, r_ddy)

       r_hv_grad_x = MATMUL(r_hv, r_ddx)
       r_hv_grad_y = MATMUL(r_hv, r_ddy)

       r_hu_meangrad(1) = dot_product(r_Psi_c, r_hu_grad_x)
       r_hu_meangrad(2) = dot_product(r_Psi_c, r_hu_grad_y)

       r_hv_meangrad(1) = dot_product(r_Psi_c, r_hv_grad_x)
       r_hv_meangrad(2) = dot_product(r_Psi_c, r_hv_grad_y)

       !-- Compute 2nd order derivatives inside the element (at centroid)
       r_tmp_xx = matmul(r_hu_grad_x, r_ddx)
       r_tmp_xy = matmul(r_hu_grad_x, r_ddy)
       r_tmp_yx = matmul(r_hu_grad_y, r_ddx)
       r_tmp_yy = matmul(r_hu_grad_y, r_ddy)

       r_tmp2_xx = matmul(r_hv_grad_x, r_ddx)
       r_tmp2_xy = matmul(r_hv_grad_x, r_ddy)
       r_tmp2_yx = matmul(r_hv_grad_y, r_ddx)
       r_tmp2_yy = matmul(r_hv_grad_y, r_ddy)

       r_hu_meanlaplace(1,1) = 0.0 !0.5*dot_product(MATMUL(r_tmp_xx ,r_epsi), r_omega)
       r_hu_meanlaplace(1,2) = 0.0 !0.5*dot_product(MATMUL(r_tmp_xy ,r_epsi), r_omega)
       r_hu_meanlaplace(2,1) = 0.0 !0.5*dot_product(MATMUL(r_tmp_yx ,r_epsi), r_omega)
       r_hu_meanlaplace(2,2) = 0.0 !0.5*dot_product(MATMUL(r_tmp_yy ,r_epsi), r_omega)

       r_hv_meanlaplace(1,1) = 0.0 !0.5*dot_product(MATMUL(r_tmp2_xx,r_epsi), r_omega)
       r_hv_meanlaplace(1,2) = 0.0 !0.5*dot_product(MATMUL(r_tmp2_xy,r_epsi), r_omega)
       r_hv_meanlaplace(2,1) = 0.0 !0.5*dot_product(MATMUL(r_tmp2_yx,r_epsi), r_omega)
       r_hv_meanlaplace(2,2) = 0.0 !0.5*dot_product(MATMUL(r_tmp2_yy,r_epsi), r_omega)

!-- c) Compute Patch Extrema

       !-- Initialize local min/max
       r_patch_min(1) = MIN( 10E+4, r_hu_mean)
       r_patch_max(1) = MAX(-10E+4, r_hu_mean)

       r_patch_min(2) = MIN( 10E+4, r_hu_meangrad(1))
       r_patch_max(2) = MAX(-10E+4, r_hu_meangrad(1))

       r_patch_min(3) = MIN( 10E+4, r_hu_meangrad(2))
       r_patch_max(3) = MAX(-10E+4, r_hu_meangrad(2))


       r_patch_min2(1) = MIN( 10E+4, r_hv_mean)
       r_patch_max2(1) = MAX(-10E+4, r_hv_mean)

       r_patch_min2(2) = MIN( 10E+4, r_hv_meangrad(1))
       r_patch_max2(2) = MAX(-10E+4, r_hv_meangrad(1))

       r_patch_min2(3) = MIN( 10E+4, r_hv_meangrad(2))
       r_patch_max2(3) = MAX(-10E+4, r_hv_meangrad(2))

       !-- Compute min/max for
       DO i_node=1,GRID_elementnodes
          DO i_elmt_p=1,GRID_patchelements
             i_patch = i_nodepatch(i_elmt_p, i_ndmaster(i_elementnodes(i_node, i_elmt)))
             IF (i_patch > 0) THEN

                !-- Compute mean value and mean grad of patch element
                r_patch_val =  r_QU(2, i_elementdofs(:,i_patch))
                r_patch_val2 = r_QU(3, i_elementdofs(:,i_patch))

                IF (i_iswet .EQ. 0) THEN

                   !DO i_dof=1, i_faceunknowns


                      r_patch_val(:)  = velocity(r_QU(1,i_elementdofs(:,i_patch)), r_QU(2,i_elementdofs(:,i_patch)))
                      r_patch_val2(:) = velocity(r_QU(1,i_elementdofs(:,i_patch)), r_QU(3,i_elementdofs(:,i_patch)))

                 DO i_dof = 1, i_faceunknowns
                      IF (r_QU(1, i_elementdofs(i_dof,i_patch)) .LE. r_wettol) THEN

                         r_patch_val(i_dof)  = 0.0
                         r_patch_val2(i_dof) = 0.0


                     END IF

                END DO !i_dof
                   END IF

                r_qval  = MATMUL(r_patch_val, r_epsi)
                r_qval2 = MATMUL(r_patch_val2, r_epsi)

                r_patch_mean  = 0.5*dot_product(r_qval , r_omega)
                r_patch_mean2 = 0.5*dot_product(r_qval2, r_omega)

                r_patch_mean_gradx  = 0.5*dot_product(MATMUL(MATMUL(r_patch_val , r_ddx), r_epsi), r_omega)
                r_patch_mean_grady  = 0.5*dot_product(MATMUL(MATMUL(r_patch_val , r_ddy), r_epsi), r_omega)
                r_patch_mean2_gradx = 0.5*dot_product(MATMUL(MATMUL(r_patch_val2, r_ddx), r_epsi), r_omega)
                r_patch_mean2_grady = 0.5*dot_product(MATMUL(MATMUL(r_patch_val2, r_ddy), r_epsi), r_omega)

                r_patch_min(1) = MIN(r_patch_min(1), r_patch_mean)
                r_patch_max(1) = MAX(r_patch_max(1), r_patch_mean)
                r_patch_min(2) = MIN(r_patch_min(2), r_patch_mean_gradx)
                r_patch_max(2) = MAX(r_patch_max(2), r_patch_mean_gradx)
                r_patch_min(3) = MIN(r_patch_min(3), r_patch_mean_grady)
                r_patch_max(3) = MAX(r_patch_max(3), r_patch_mean_grady)

                r_patch_min2(1) = MIN(r_patch_min2(1), r_patch_mean2)
                r_patch_max2(1) = MAX(r_patch_max2(1), r_patch_mean2)
                r_patch_min2(2) = MIN(r_patch_min2(2), r_patch_mean2_gradx)
                r_patch_max2(2) = MAX(r_patch_max2(2), r_patch_mean2_gradx)
                r_patch_min2(3) = MIN(r_patch_min2(3), r_patch_mean2_grady)
                r_patch_max2(3) = MAX(r_patch_max2(3), r_patch_mean2_grady)
             END IF
          END DO !i_elmt_p
       END DO !i_node

!-- d) Compute alphas
       r_alpha = 0.0
       r_alpha2 = 0.0

       r_alpha(1) = compute_alpha(r_hu       , r_hu_mean       , r_patch_min(1), r_patch_max(1), i_faceunknowns)
       r_alpha_x  = compute_alpha(r_hu_grad_x, r_hu_meangrad(1), r_patch_min(2), r_patch_max(2), i_faceunknowns)
       r_alpha_y  = compute_alpha(r_hu_grad_y, r_hu_meangrad(2), r_patch_min(3), r_patch_max(3), i_faceunknowns)
       r_alpha(2) = MIN(r_alpha_x, r_alpha_y)

       r_alpha2(1) = compute_alpha(r_hv       , r_hv_mean       , r_patch_min2(1), r_patch_max2(1), i_faceunknowns)
       r_alpha2_x  = compute_alpha(r_hv_grad_x, r_hv_meangrad(1), r_patch_min2(2), r_patch_max2(2), i_faceunknowns)
       r_alpha2_y  = compute_alpha(r_hv_grad_y, r_hv_meangrad(2), r_patch_min2(3), r_patch_max2(3), i_faceunknowns)
       r_alpha2(2) = MIN(r_alpha2_x, r_alpha2_y)

!-- e) Taylor expansion

       r_hu_lim = assemble_taylor(r_hu_mean, r_hu_meangrad(1), r_hu_meangrad(2), r_hu_meanlaplace, &
       r_epsi, r_alpha,  r_b, r_coox, r_cooy, i_faceunknowns)

       r_hv_lim = assemble_taylor(r_hv_mean, r_hv_meangrad(1), r_hv_meangrad(2), r_hv_meanlaplace, &
       r_epsi, r_alpha2, r_b, r_coox, r_cooy, i_faceunknowns)

!-- f) Assign limited value to vector r_QU

       r_QU(1, i_elementdofs(:,i_elmt)) = r_QU_lim(i_elementdofs(:,i_elmt))

       DO i_dof=1, i_faceunknowns
           r_QU(2, i_elementdofs(i_dof,i_elmt)) = r_hu_lim(i_dof)*(i_iswet) + (1-i_iswet)*(r_hu_lim(i_dof) * &
                r_QU_lim(i_elementdofs(i_dof,i_elmt)))
           r_QU(3, i_elementdofs(i_dof,i_elmt)) = r_hv_lim(i_dof)*(i_iswet) + (1-i_iswet)*(r_hv_lim(i_dof) * &
                r_QU_lim(i_elementdofs(i_dof,i_elmt)))

          IF (ABS(r_QU_lim(i_elementdofs(i_dof,i_elmt))) .LE. r_wettol) THEN
             r_QU(2,i_elementdofs(i_dof,i_elmt)) = 0.0
             r_QU(3,i_elementdofs(i_dof,i_elmt)) = 0.0
          END IF
       END DO !i_dof

    END DO !i_elmt

!-- Deallocate workspace
    DEALLOCATE(i_ndmaster, i_nodepatch)

  END SUBROUTINE

!*********************************************************************************************************
  FUNCTION compute_alpha(r_var, r_mean, r_patch_min, r_patch_max, i_faceunknowns)

  IMPLICIT NONE

  INTEGER(KIND=GRID_SI), INTENT(in)			:: i_faceunknowns
  REAL(KIND=GRID_SR), DIMENSION(:), INTENT(in)		:: r_var
  REAL(KIND=GRID_SR), INTENT(in)    			:: r_mean, r_patch_min, r_patch_max
  INTEGER(KIND=GRID_SI)					:: i_cnt
  REAL(KIND=GRID_SR), DIMENSION(i_faceunknowns)		:: r_tmp, r_tmp2
  REAL(KIND=GRID_SR)					:: r_var_tmp, compute_alpha, r_tol, r_diff

  r_tol = 10E-12
  r_tmp = 0.0

#ifdef BERNSTEIN
     r_tmp2 = r_var
     r_tmp2 = bernstein_on_intpoints(r_tmp2, i_faceunknowns)
#endif

  DO i_cnt=1, 3 !i_faceunknowns
     r_var_tmp = r_var(i_cnt)

#ifdef BERNSTEIN
     r_var_tmp = r_tmp2(i_cnt)
#endif

     r_diff = r_var_tmp-r_mean
     IF (abs(r_diff) .LE. r_tol) THEN
        r_tmp(i_cnt) = 1._GRID_SR
     ELSEIF (r_diff .GT. r_tol) THEN
     	r_tmp(i_cnt) =  MIN(1._GRID_SR, (r_patch_max-r_mean)/r_diff)
     ELSEIF (r_diff .LT. r_tol) THEN
        r_tmp(i_cnt) =  MIN(1._GRID_SR, (r_patch_min-r_mean)/r_diff)
     END IF
  END DO !i_cnt

  compute_alpha = MINVAL(r_tmp(1:3))

  RETURN
  END FUNCTION compute_alpha
!*********************************************************************************************************
  FUNCTION assemble_taylor(r_mean, r_gradx, r_grady, r_laplace, r_epsi, &
  	   		   r_alpha, r_centr, r_coox, r_cooy, i_faceunknowns) RESULT(r_val)

  IMPLICIT NONE

  INTEGER(KIND=GRID_SI), INTENT(in)	   	    	:: i_faceunknowns
  REAL(KIND=GRID_SR), INTENT(in)			:: r_mean
  REAL(KIND=GRID_SR), DIMENSION(2), INTENT(in)		:: r_alpha
  REAL(KIND=GRID_SR), DIMENSION(2,2), INTENT(in)	:: r_laplace
  REAL(KIND=GRID_SR), INTENT(in)			:: r_gradx, r_grady
  REAL(KIND=GRID_SR), DIMENSION(2), INTENT(in)  	:: r_centr
  REAL(KIND=GRID_SR), DIMENSION(:), INTENT(in)		:: r_coox, r_cooy
  REAL(KIND=GRID_SR), DIMENSION(:,:), INTENT(in)        :: r_epsi

!-- Local declaration
  INTEGER(KIND=GRID_SI)		    			:: i_cnt
  REAL(KIND=GRID_SR)  					:: r_delta_x, r_delta_y
  REAL(KIND=GRID_SR), DIMENSION(i_faceunknowns)		:: r_val
  REAL(KIND=GRID_SR)                                    :: r_tmp
  REAL(KIND=GRID_SR)                                    :: r_B1, r_B2, r_B3, r_B4, r_B5, r_B6
  REAL(KIND=GRID_SR)                                    :: r_B4_tmp, r_B5_tmp, r_B6_tmp

  REAL(KIND=GRID_SR), DIMENSION(i_faceunknowns)                 :: r_weigh

!-- Quickfix get quad weights
  r_weigh = GRID_femtypes%p_type(FEM_DG)%sig%r_equadwei

  !-- Taylor expansion as in Kuzmins paper equation (16)
  IF ((r_alpha(1) .GT. 1.0 ) .OR. (r_alpha(1) .LT. 0.0)) write(*,*) 'wrong alpha 1 ', r_alpha(1)
  IF ((r_alpha(2) .GT. 1.0 ) .OR. (r_alpha(2) .LT. 0.0)) write(*,*) 'wrong alpha 2 ', r_alpha(2)

  !-- Compute correction term as in Luo et al. (2008) eq. (12)
  r_B4_tmp = overbar(r_coox, r_coox, r_centr(1), r_centr(1), i_faceunknowns, r_weigh, r_epsi)
  r_B5_tmp = overbar(r_cooy, r_cooy, r_centr(2), r_centr(2), i_faceunknowns, r_weigh, r_epsi)
  r_B6_tmp = overbar(r_coox, r_cooy, r_centr(1), r_centr(2), i_faceunknowns, r_weigh, r_epsi)

  DO i_cnt=1,i_faceunknowns

     !-- Assemble the Taylor basis
     r_delta_x = r_coox(i_cnt)-r_centr(1)
     r_delta_y = r_cooy(i_cnt)-r_centr(2)

     r_B1 = 1.0
     r_B2 = r_delta_x
     r_B3 = r_delta_y
     r_B4 = 0.5*(r_delta_x**2) - 0.5*r_B4_tmp
     r_B5 = 0.5*(r_delta_y**2) - 0.5*r_B5_tmp
     r_B6 = r_delta_y*r_delta_x -  r_B6_tmp

     r_val(i_cnt) = r_mean*r_B1 + &
          r_alpha(1)*r_gradx*r_B2 + &
          r_alpha(1)*r_grady*r_B3 + &
          r_alpha(2)*r_laplace(1,1)*r_B4 + &
          r_alpha(2)*r_laplace(2,2)*r_B5 + &
          r_alpha(2)*r_laplace(2,1)*r_B6

  END DO !i_cnt


  !--Map from function evaluations to coefficients
#ifdef BERNSTEIN
  r_val = bernstein_get_flambda(r_val, i_faceunknowns)
#endif

  RETURN
END FUNCTION assemble_taylor
!*********************************************************************************************************
  FUNCTION overbar(r_x, r_y, r_xc, r_yc, i_numquadpoints, r_weigh, r_epsi) RESULT(r_basis)

  IMPLICIT NONE

  INTEGER(KIND=GRID_SI), INTENT(in)		:: i_numquadpoints

  REAL(KIND=GRID_SR), INTENT(in)		:: r_xc, r_yc
  REAL(KIND=GRID_SR), DIMENSION(:), INTENT(in)	:: r_x, r_y
  REAL(KIND=GRID_SR), DIMENSION(:), INTENT(in)  :: r_weigh
  REAL(KIND=GRID_SR), DIMENSION(:,:), INTENT(in):: r_epsi

  REAL(KIND=GRID_SR)  		    		:: r_basis

!-- Local declarations
  INTEGER(KIND=GRID_SI)				:: i_cnt
  REAL(KIND=GRID_SR), DIMENSION(i_numquadpoints):: r_x_quad, r_y_quad

  r_basis = 0._GRID_SR

!-- Interpolate onto quadrature points
  r_x_quad = MATMUL(r_x, r_epsi)
  r_y_quad = MATMUL(r_y, r_epsi)

  DO i_cnt=1, i_numquadpoints
     r_basis = r_basis + 0.5_GRID_SR * r_weigh(i_cnt) * (r_x_quad(i_cnt) - r_xc) * (r_y_quad(i_cnt) - r_yc)
  END DO !i_cnt

  RETURN
  END FUNCTION overbar
!*********************************************************************************************************
END MODULE

