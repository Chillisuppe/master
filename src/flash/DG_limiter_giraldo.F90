!*******************************************************************************
!
!> @file  DG_limiter_giraldo.F90
!> @brief contains module DG_limiter
!
!*******************************************************************************
! MODULE DESCRIPTION:
!> @brief function to limit discrete DG solution
!
MODULE DG_limiter

  USE GRID_api
  USE FLASH_parameters

  PRIVATE
  PUBLIC  :: limiter

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
  SUBROUTINE limiter(p_ghand, r_Q, i_elementedges, i_edgeinfo, r_edgelength, r_metrics_inv, &
                     i_numelmt, r_coonod, r_coodof, i_edgenodes, i_elementnodes, i_elementdofs, &
                     i_numedge, r_S, i_faceunknowns, r_evander, r_einvvander)

    IMPLICIT NONE

    TYPE (grid_handle), INTENT(in)                            :: p_ghand
    REAL(KIND = GRID_SR), DIMENSION(:,:), INTENT(inout)       :: r_Q
    INTEGER(KIND = GRID_SI), DIMENSION(:,:), INTENT(in)       :: i_elementedges, i_edgeinfo, i_edgenodes, &
                                                                 i_elementnodes, i_elementdofs
    REAL(KIND = GRID_SR), DIMENSION(:), INTENT(in)            :: r_edgelength
    REAL(KIND = GRID_SR), DIMENSION(:,:,:), INTENT(in)        :: r_metrics_inv
    REAL(KIND = GRID_SR), DIMENSION(:,:), INTENT(in)          :: r_coonod, r_coodof, r_evander, r_einvvander, r_S
    INTEGER(KIND = GRID_SI), INTENT(in)                       :: i_numelmt, i_faceunknowns, &
                                                                 i_numedge

    INTEGER(KIND = GRID_SI)                                   :: i_alct, i_elmt, i_cnt, j_cnt, &
                                                                 i_edge, i_elmt_l, i_elmt_r
    REAL(KIND = GRID_SR)                                      :: r_elmt_min, r_elmt_max, &
                                                                 r_h, r_patch_min, r_patch_max
    REAL(KIND = GRID_SR), DIMENSION(i_faceunknowns)           :: r_modal, r_weight
    REAL(KIND = GRID_SR), DIMENSION(:), ALLOCATABLE           :: r_elmt_mean

!--- allocate workspace
    ALLOCATE(r_elmt_mean(i_numelmt), stat=i_alct)

    r_weight     = 0.0_GRID_SR
    r_weight(1)  = 1.0_GRID_SR

!--- element mean values
    DO i_elmt=1, i_numelmt
        r_elmt_mean(i_elmt) = DOT_PRODUCT(r_evander(:,1), r_Q(1, i_elementdofs(:, i_elmt))) / SUM(r_evander(:,1))
    END DO !i_elmt

!--- go through all elements
    DO i_elmt=1, i_numelmt

!--- compute minimum/maximum on every triangle
      r_elmt_min = r_elmt_mean(i_elmt)
      r_elmt_max = r_elmt_mean(i_elmt)

      DO i_cnt=1, i_faceunknowns
        r_h        = r_Q(1, i_elementdofs(i_cnt, i_elmt))
        r_elmt_min = min(r_h, r_elmt_min)
        r_elmt_max = max(r_h, r_elmt_max)
      END DO !i_cnt

!--- check for min/max of mean values in surrounding elements
      r_patch_min = r_elmt_mean(i_elmt)
      r_patch_max = r_elmt_mean(i_elmt)

      DO j_cnt=1, 3
        i_edge   = i_elementedges(j_cnt, i_elmt)
        i_elmt_l = i_edgeinfo(1, i_edge)
        i_elmt_r = i_edgeinfo(2, i_edge)

!--- swap elements if are we "right" from the edge
        IF(i_elmt .EQ. i_elmt_r) THEN
          i_elmt_r = i_elmt_l
        END IF

!--- compute max/min for the neighborhood
        IF (i_elmt_r .GT. 0) THEN
          r_patch_min = min(r_patch_min, r_elmt_mean(i_elmt_r))
          r_patch_max = max(r_patch_max, r_elmt_mean(i_elmt_r))
        END IF
      END DO !j_cnt

!--- are the extrema inside the element?
      IF((r_elmt_min < r_patch_min) .OR. (r_elmt_max > r_patch_max)) THEN

!--- map nodes to modes (first mode = mean)
        r_modal = MATMUL(r_Q(1, i_elementdofs(:, i_elmt)), r_evander)

!--- map modes back to nodes by weighting the modes
        DO i_cnt=1, i_faceunknowns
          r_modal(i_cnt) = r_modal(i_cnt)*r_weight(i_cnt)
        END DO

        r_Q(1, i_elementdofs(:, i_elmt)) = MATMUL(r_modal, r_einvvander)
      END IF

    END DO !i_elmt

    DEALLOCATE(r_elmt_mean)

  END SUBROUTINE limiter

!*******************************************************************************
END MODULE DG_limiter
