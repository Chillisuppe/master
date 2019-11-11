!*******************************************************************************
!
!> @file  DG_limiter_shu.F90
!> @brief contains module DG_limiter
!
!*******************************************************************************
! MODULE DESCRIPTION:
!> @brief function to limit discrete DG solution
!
MODULE DG_limiter

  USE GRID_api
  USE FLASH_parameters
  USE DG_limiter_utils

  PRIVATE
  PUBLIC  :: limiter

  CONTAINS

!*******************************************************************************
! DESCRIPTION of [SUBROUTINE limiter]:
!> @brief flux limiter for DG method
!>
!> Slope limiter of Cockburn and Shu (1998)
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

    REAL(KIND = GRID_SR)                                      :: r_nu, r_MM, r_det, &
                                                                 r_pos, r_neg, &
                                                                 r_theta_pl, r_theta_mn, r_deltax
    REAL(KIND = GRID_SR), DIMENSION(2)                        :: r_Deltab1, r_Deltab2, r_alpha
    REAL(KIND = GRID_SR), DIMENSION(3)                        :: r_Delta_hat, r_Deltaphibar, &
                                                                 r_phitilde, r_Deltaphi1, r_Deltaphi2
    REAL(KIND = GRID_SR), DIMENSION(2,2)                      :: r_matrix
    REAL(KIND = GRID_SR), DIMENSION(3,3)                      :: r_Delta
    REAL(KIND = GRID_SR), DIMENSION(:,:), ALLOCATABLE         :: r_vander_m, r_basis, &
                                                                 r_m, r_b, r_db, r_emean, r_dphi
    INTEGER(KIND = GRID_SI)                                   :: i_alct, i_loop, i_eqn, &
                                                                 i_degree, i_gloc, i_edge, i_edge2, &
                                                                 i_elmt

!--- initialize values, constants, allocate workspace
    r_nu = 1.5_GRID_SR    ! as in the paper by Cockburn & Shu
    r_MM = 50.0_GRID_SR

!--- allocate workspace
    ALLOCATE(r_vander_m(i_faceunknowns, 3), &
             r_basis(i_faceunknowns, 3), &
             r_b(2, i_numelmt), r_db(2, i_numedge), r_m(2, i_numedge), &
             r_emean(3, i_numelmt), r_dphi(3, i_numedge), stat=i_alct)
    IF (i_alct /= 0) CALL grid_error(c_error='[DG_limiter_shu]: Could not allocate enough memory')

!--- get info from signature and initialize
    i_degree  = GRID_femtypes%p_type(FEM_DG)%sig%i_degree

!---  matrix for converting basis coefficients to edge midpoint values
    SELECT CASE(i_degree)
      CASE(1)
        r_vander_m(:,1) = (/ 0.0, 0.5, 0.5 /)
        r_vander_m(:,2) = (/ 0.5, 0.0, 0.5 /)
        r_vander_m(:,3) = (/ 0.5, 0.5, 0.0 /)
      CASE(2)
        r_vander_m(:,1) = (/ 0.0, 0.0, 0.0, 1.0, 0.0, 0.0 /)
        r_vander_m(:,2) = (/ 0.0, 0.0, 0.0, 0.0, 1.0, 0.0 /)
        r_vander_m(:,3) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 /)
    END SELECT

!--- matrix for converting edge midpoint values to Lagrange points (basis coefficients)
    SELECT CASE (i_degree)
      CASE(1)
        r_basis(1,:) = (/ -1,  1,  1 /)
        r_basis(2,:) = (/  1, -1,  1 /)
        r_basis(3,:) = (/  1,  1, -1 /)
      CASE(2)
        r_basis(1,:) = (/ -1,  1,  1 /)
        r_basis(2,:) = (/  1, -1,  1 /)
        r_basis(3,:) = (/  1,  1, -1 /)
        r_basis(4,:) = (/  1,  0,  0 /)
        r_basis(5,:) = (/  0,  1,  0 /)
        r_basis(6,:) = (/  0,  0,  1 /)
    END SELECT

    r_b = 0._GRID_SR
    r_m = 0._GRID_SR

!--- compute mean values and barycenter for each element
    DO i_elmt=1, i_numelmt
      r_emean(1,i_elmt)= DOT_PRODUCT(r_evander(:,1), r_Q(1,i_elementdofs(:,i_elmt)) + &
                                     r_S(1,i_elementdofs(:,i_elmt))) / SUM(r_evander(:,1))
      r_emean(2,i_elmt)= DOT_PRODUCT(r_evander(:,1), r_Q(2,i_elementdofs(:,i_elmt))) / SUM(r_evander(:,1))
      r_emean(3,i_elmt)= DOT_PRODUCT(r_evander(:,1), r_Q(3,i_elementdofs(:,i_elmt))) / SUM(r_evander(:,1))
      r_b(:,i_elmt) = (r_coonod(:,i_elementnodes(1,i_elmt)) + r_coonod(:,i_elementnodes(2,i_elmt)) + &
                       r_coonod(:,i_elementnodes(3,i_elmt))) / 3.0_GRID_SR
    END DO !i_elmt

!--- compute coordinates of edge centres
    DO i_edge=1, i_numedge
      r_m(:,i_edge)  = 0.5 * (r_coonod(:,i_edgenodes(1,i_edge)) + r_coonod(:,i_edgenodes(2,i_edge)))
      IF ((i_edgeinfo(3,i_edge) /=0).AND.(i_edgeinfo(4,i_edge) /= 0)) THEN
        r_db(:,i_edge)   = r_b(:,i_edgeinfo(4,i_edge)) - r_b(:,i_edgeinfo(3,i_edge))
        r_dphi(:,i_edge) = r_emean(:,i_edgeinfo(4,i_edge))-r_emean(:,i_edgeinfo(3,i_edge))
      ELSE
        r_db(:,i_edge)   = r_m(:,i_edge) - r_b(:,i_edgeinfo(3,i_edge))
        i_gloc = i_edgeinfo(1,i_edge)
        i_elmt = i_edgeinfo(3,i_edge)
        r_dphi(1,i_edge) = DOT_PRODUCT(r_vander_m(:,i_gloc) , r_Q(1,i_elementdofs(:,i_elmt)) + &
                                       r_S(1,i_elementdofs(:,i_elmt))) - r_emean(1,i_elmt)
        r_dphi(2,i_edge) = DOT_PRODUCT(r_vander_m(:,i_gloc) , r_Q(2,i_elementdofs(:,i_elmt))) - r_emean(2,i_elmt)
        r_dphi(3,i_edge) = DOT_PRODUCT(r_vander_m(:,i_gloc) , r_Q(3,i_elementdofs(:,i_elmt))) - r_emean(3,i_elmt)
      END IF
    END DO

!--- loop over all elements
    elmt_loop: DO i_elmt=1, i_numelmt

      edge_loop: DO i_gloc=1,3

        i_edge = i_elementedges(i_gloc,i_elmt)

!--- get differences in the mean of edge-neighbours
        IF (i_edgeinfo(3,i_edge) == i_elmt) THEN
          r_Deltab1   = r_db(:,i_edge)
          r_Deltaphi1 = r_dphi(:,i_edge)
        ELSE
          r_Deltab1   = -r_db(:,i_edge)
          r_Deltaphi1 = -r_dphi(:,i_edge)
        END IF

!--- choose second neighbour for basis representation
        DO i_loop = 0,1

          i_edge2 = i_elementedges(MOD(i_gloc+i_loop,3)+1,i_elmt)

!--- get differences in the mean of edge-neighbours
          IF (i_edgeinfo(3,i_edge2) == i_elmt) THEN
            r_Deltab2   = r_db(:,i_edge2)
            r_Deltaphi2 = r_dphi(:,i_edge2)
          ELSE
            r_Deltab2   = -r_db(:,i_edge2)
            r_Deltaphi2 = -r_dphi(:,i_edge2)
          END IF

!--- determine coefficients alpha
          r_det = r_Deltab1(1) * r_Deltab2(2) - r_Deltab1(2) * r_Deltab2(1)

          r_matrix(1,1) =  r_Deltab2(2)/r_det
          r_matrix(1,2) = -r_Deltab2(1)/r_det
          r_matrix(2,1) = -r_Deltab1(2)/r_det
          r_matrix(2,2) =  r_Deltab1(1)/r_det

          r_alpha = MATMUL(r_matrix, r_m(:,i_edge)-r_b(:,i_elmt))

          IF ((r_alpha(1) >= 0) .AND. (r_alpha(2) >= 0)) THEN
            EXIT
          END IF
        END DO

!--- determine r_Deltaphibar
        r_Deltaphibar = r_alpha(1)*r_Deltaphi1 + r_alpha(2)*r_Deltaphi2

!--- compute r_phitilde for local edge
        r_phitilde(1) = DOT_PRODUCT(r_vander_m(:,i_gloc) , r_Q(1,i_elementdofs(:,i_elmt)) + &
                                    r_S(1,i_elementdofs(:,i_elmt))) - r_emean(1,i_elmt)
        r_phitilde(2) = DOT_PRODUCT(r_vander_m(:,i_gloc) , r_Q(2,i_elementdofs(:,i_elmt))) - r_emean(2,i_elmt)
        r_phitilde(3) = DOT_PRODUCT(r_vander_m(:,i_gloc) , r_Q(3,i_elementdofs(:,i_elmt))) - r_emean(3,i_elmt)

!--- compute Delta_i
        r_deltax = r_edgelength(i_edge)

        r_Delta(1,i_gloc) = TVB_minmod(r_phitilde(1), r_nu*r_Deltaphibar(1), r_MM, r_deltax)
        r_Delta(2,i_gloc) = TVB_minmod(r_phitilde(2), r_nu*r_Deltaphibar(2), r_MM, r_deltax)
        r_Delta(3,i_gloc) = TVB_minmod(r_phitilde(3), r_nu*r_Deltaphibar(3), r_MM, r_deltax)

      END DO edge_loop

!--- IF sum Delta_i = ...
      DO i_eqn=1,3
        IF (ABS(SUM(r_Delta(i_eqn,:))) <= 1E-10) THEN
          r_Q(i_eqn,i_elementdofs(:,i_elmt)) = r_emean(i_eqn,i_elmt) + MATMUL(r_basis, r_Delta(i_eqn,:))
        ELSE

          r_pos = MAX(0._GRID_SR, r_Delta(i_eqn,1)) + MAX(0._GRID_SR, r_Delta(i_eqn,2)) + MAX(0._GRID_SR, r_Delta(i_eqn,3))
          r_neg = MAX(0._GRID_SR,-r_Delta(i_eqn,1)) + MAX(0._GRID_SR,-r_Delta(i_eqn,2)) + MAX(0._GRID_SR,-r_Delta(i_eqn,3))

          ! Note (Stefan V.): Sometimes r_pos or r_neg can be zero, which must be handled to avoid
          !                   division by zero. This problem is not mentioned in the literature.
          IF(r_pos == 0.0_GRID_SR) THEN
            r_theta_pl = 1._GRID_SR
          ELSE
            r_theta_pl = MIN(1._GRID_SR, r_neg/r_pos)
          END IF

          IF(r_neg == 0.0_GRID_SR) THEN
            r_theta_mn = 1._GRID_SR
          ELSE
            r_theta_mn = MIN(1._GRID_SR, r_pos/r_neg)
          END IF

          DO i_gloc=1,3
            r_Delta_hat(i_gloc)= r_theta_pl*MAX(0._GRID_SR, r_Delta(i_eqn,i_gloc)) - &
                                 r_theta_mn*MAX(0._GRID_SR,-r_Delta(i_eqn,i_gloc))
          END DO

          r_Q(i_eqn,i_elementdofs(:,i_elmt)) = r_emean(i_eqn,i_elmt) + MATMUL(r_basis, r_Delta_hat)
        END IF !sum Delta=0
      END DO

      r_Q(1,i_elementdofs(:,i_elmt)) = r_Q(1,i_elementdofs(:,i_elmt)) - r_S(1,i_elementdofs(:,i_elmt))

    END DO elmt_loop

!--- deallocate workspace
    DEALLOCATE(r_vander_m, r_basis, r_emean, r_b, r_db, r_m, r_dphi)

  END SUBROUTINE limiter

!*******************************************************************************
END MODULE DG_limiter
