!*******************************************************************************
!
!> @file  DG_utils.F90
!> @brief contains module DG_utils
!
!*******************************************************************************
!
! VERSION(s):
! 1. version            nicole beisiegel        02/2012
!
!*******************************************************************************
! MODULE DESCRIPTION:
!> @brief provides some routines for DG method
!
MODULE DG_utils

  USE GRID_api

  PRIVATE
  PUBLIC  :: compute_edgeinfo, fvm_createmetrics, &
             heaviside, bary, Vandermonde2D

  CONTAINS

!*******************************************************************************
! DESCRIPTION of [SUBROUTINE compute_edgeinfo]:
!> @brief Bookkeeping routine for detecting dofs on edges
!> This is what we do..
!
!> @param[in]     i_numedge         total number of edges
!> @param[in]     i_elementedges    edge indices of the edges of each element; DIMENSION(3,i_numelmt)
!> @param[in]     i_edgeboundary    boundary flags for all edges; DIMENSION(i_numedge)
!> @param[in,out] i_edgeinfo        edge-element relation for each edge: local edge index and global element index for adjacent elements; DIMENSION(4,i_numedge)
!
  SUBROUTINE compute_edgeinfo(i_numedge, i_elementedges, i_edgeboundary, i_edgeinfo)

    IMPLICIT NONE

!---------- local declarations

    INTEGER (KIND = GRID_SI), INTENT(in)                    :: i_numedge
    INTEGER (KIND = GRID_SI), DIMENSION(:,:), INTENT(in)    :: i_elementedges
    INTEGER (KIND = GRID_SI), DIMENSION(:), INTENT(in)      :: i_edgeboundary
    INTEGER (KIND = GRID_SI), DIMENSION(:,:), INTENT(inout) :: i_edgeinfo

    INTEGER (KIND = GRID_SI)                                :: i_edge, i_elmt_l, i_elmt_r, &
                                                               i_eledge, i_cnt, i_alct
    INTEGER (KIND = GRID_SI), DIMENSION(:), ALLOCATABLE     :: i_edgemaster

!--- allocate arrays
    ALLOCATE(i_edgemaster(i_numedge), stat=i_alct)
    IF(i_alct /= 0) CALL grid_error(c_error='[compute_edgeinfo]: Could not allocate edgemaster array')

!--- Initialize global arrays

    i_edgemaster = (/ (i_cnt, i_cnt=1,i_numedge) /)
    i_edgemaster = MERGE(i_edgeboundary, i_edgemaster, i_edgeboundary > 0)

    DO i_edge=1, i_numedge

!--- Store the left and the right element of each edge.

      i_elmt_l = i_edgeinfo(3,i_edge)
      i_elmt_r = i_edgeinfo(4,i_edge)

!--- if this edge is a (true) boundary edge then store 0 as the right (left resp.) element
!--- store the local number of this edge for the left and right element

      i_edgeinfo(1:2,i_edge) = 0_GRID_SI

      DO i_eledge = 1, GRID_elementedges
        IF(i_edgemaster(i_elementedges(i_eledge,i_elmt_l)) == i_edgemaster(i_edge)) i_edgeinfo(1,i_edge) = i_eledge
        IF (i_elmt_r > 0) THEN
          IF(i_edgemaster(i_elementedges(i_eledge,i_elmt_r)) == i_edgemaster(i_edge)) i_edgeinfo(2,i_edge) = i_eledge
        END IF
      END DO !i_eledge
!       write (*,'(i4 i4 i6 i4)') i_edgeinfo(1,i_edge), i_edgeinfo(2,i_edge), i_edgeinfo(3,i_edge), i_edgeinfo(4,i_edge)

    END DO !i_edge

  END SUBROUTINE compute_edgeinfo

!*******************************************************************************
! DESCRIPTION of [SUBROUTINE fvm_createmetrics]:
!> @brief
!
!> @param[in]     p_ghand         grid handling data structure
!> @param[out]    r_metrics_inv   metric terms for each element; DIMENSION(2,2,i_numelmt)
!> @param[out]    r_elmtvolume    volume of each element; DIMENSION(i_numelmt)
!> @param[out]    r_edgelength
!> @param[out]    r_normals
!> @param[in]     i_edgeinfo      edge-element relation for each edge: local edge index and global element index for adjacent elements; DIMENSION(4,i_numedge)
!> @param[in]     i_elmtnodes     node indices of each element; DIMENSION(3, i_numelmt)
!> @param[in]     r_coonod
!
  SUBROUTINE fvm_createmetrics(p_ghand, r_metrics_inv, r_elmtvolume, r_edgelength, r_hmin, &
                               r_normals, i_edgeinfo, i_elmtnodes, i_elmtedges, r_coonod)

    IMPLICIT NONE

    TYPE (grid_handle), INTENT(in)                        :: p_ghand
    REAL (KIND = GRID_SR), DIMENSION(:,:,:), INTENT(out)  :: r_metrics_inv
    REAL (KIND = GRID_SR), DIMENSION(:), INTENT(out)      :: r_elmtvolume, r_edgelength, r_hmin
    REAL (KIND = GRID_SR), DIMENSION(:,:), INTENT(out)    :: r_normals
    INTEGER (KIND = GRID_SI), DIMENSION(:,:), INTENT(in)  :: i_edgeinfo, i_elmtnodes, i_elmtedges
    REAL (KIND = GRID_SR), DIMENSION(:,:), INTENT(in)     :: r_coonod

    INTEGER (KIND = GRID_SI)                              :: i_elmt, i_edge, i_numelmt, i_numedge, &
                                                             i_node1, i_node2
    REAL (KIND = GRID_SR), DIMENSION(3)                   :: r_x, r_y
    REAL (KIND = GRID_SR)                                 :: r_det_inv
    REAL (KIND = GRID_SR), DIMENSION(GRID_dimension,GRID_dimension) :: r_metrics

    i_numelmt = p_ghand%i_enumfine
    i_numedge = p_ghand%i_gnumfine

    r_hmin = HUGE(1.0_GRID_SR)

!--- calculate normals and edge length
    DO i_edge=1,i_numedge

      r_x = r_coonod(1, i_elmtnodes(:,i_edgeinfo(3,i_edge)))
      r_y = r_coonod(2, i_elmtnodes(:,i_edgeinfo(3,i_edge)))

      i_node1 = MOD(i_edgeinfo(1,i_edge)+0, GRID_elementedges) + 1
      i_node2 = MOD(i_edgeinfo(1,i_edge)+1, GRID_elementedges) + 1

      r_normals(1,i_edge) = r_y(i_node2) - r_y(i_node1)
      r_normals(2,i_edge) = r_x(i_node1) - r_x(i_node2)

      r_edgelength(i_edge) = SQRT(SUM(r_normals(:,i_edge)**2))

!--- set normals to unit length
      r_normals(:,i_edge) = r_normals(:,i_edge) / r_edgelength(i_edge)

    END DO ! i_edge

!--- calculate metric terms
    DO i_elmt=1,i_numelmt

      r_x = r_coonod(1, i_elmtnodes(:,i_elmt))
      r_y = r_coonod(2, i_elmtnodes(:,i_elmt))

!--- matrix from Axi+b = x i.e. A = dx/dxi used for transformation of the
!--- gradient A*grad_x = grad_xi
      r_metrics(1,1) = (r_x(2)-r_x(1)) / 2.0_GRID_SR
      r_metrics(2,1) = (r_y(2)-r_y(1)) / 2.0_GRID_SR
      r_metrics(1,2) = (r_x(3)-r_x(1)) / 2.0_GRID_SR
      r_metrics(2,2) = (r_y(3)-r_y(1)) / 2.0_GRID_SR

!--- this comes from Ax+b = xi \Rightarrow dxi/dx is this matrix
      r_metrics_inv(1,1,i_elmt) =  r_metrics(2,2)
      r_metrics_inv(2,1,i_elmt) = -r_metrics(1,2)
      r_metrics_inv(1,2,i_elmt) = -r_metrics(2,1)
      r_metrics_inv(2,2,i_elmt) =  r_metrics(1,1)

      r_det_inv = r_metrics(1,1)*r_metrics(2,2) - r_metrics(2,1)*r_metrics(1,2)
      r_metrics_inv(:,:,i_elmt) = r_metrics_inv(:,:,i_elmt) / r_det_inv
      r_elmtvolume(i_elmt) = ABS(r_det_inv) * 2.0_GRID_SR

!--- compute heights of triangles for cfl estimation from volume formula: A = a * h / 2
      r_hmin(i_elmtnodes(1,i_elmt)) = MIN(r_hmin(i_elmtnodes(1,i_elmt)), &
                                          2.0_GRID_SR * r_elmtvolume(i_elmt) / r_edgelength(i_elmtedges(1,i_elmt)))
      r_hmin(i_elmtnodes(2,i_elmt)) = MIN(r_hmin(i_elmtnodes(2,i_elmt)), &
                                          2.0_GRID_SR * r_elmtvolume(i_elmt) / r_edgelength(i_elmtedges(2,i_elmt)))
      r_hmin(i_elmtnodes(3,i_elmt)) = MIN(r_hmin(i_elmtnodes(3,i_elmt)), &
                                          2.0_GRID_SR * r_elmtvolume(i_elmt) / r_edgelength(i_elmtedges(3,i_elmt)))

    END DO ! i_elmt

  END SUBROUTINE fvm_createmetrics

!*******************************************************************************
! DESCRIPTION of [FUNCTION heaviside]:
!> @brief
!
!> @param         r_real
!
  FUNCTION heaviside(r_real)

    IMPLICIT NONE

    REAL (KIND = GRID_SR)           :: r_real
    REAL (KIND = GRID_SR)           :: heaviside

    heaviside = MERGE(1._GRID_SR, 0._GRID_SR, r_real <= 10E-10)

  END FUNCTION heaviside

!*******************************************************************************
! DESCRIPTION of [FUNCTION bary]:
!> @brief calculates barycentric coordinates for a given point and reference triangle
!>
!> @param[in]       r_coop      point coordinate
!> @param[in]       r_coov      coordinates of triangle vertices
!> @return          r_bary      resulting barycentric point coordinates
!
!> @note adapted from amatos
!
  FUNCTION bary(r_coop, r_coov) RESULT(r_bary)

    IMPLICIT NONE

    REAL (KIND = GRID_SR), DIMENSION(GRID_dimension), INTENT(in)    :: r_coop
    REAL (KIND = GRID_SR), DIMENSION(GRID_dimension, 3), INTENT(in) :: r_coov
    REAL (KIND = GRID_SR), DIMENSION(3)                             :: r_bary

    INTEGER (KIND = GRID_SI)                                        :: i_cnt
    REAL (KIND = GRID_SR)                                           :: r_dotg1g1, r_dotg3g3, r_dotg1g3, r_dotg1po, r_dotg3po
    REAL (KIND = GRID_SR)                                           :: r_areasq, r_sigma, r_tau
    REAL (KIND = GRID_SR), DIMENSION(GRID_dimension)                :: r_dcooedge1, r_dcooedge3, r_dcoopoint

!--- compute vectors of the two edges adjacent to the first vertex (bear order in mind!)
!--- (cf. Abb.4.1 in Diploma thesis of Thomas Heinze)
!--- and vector of point relative to the reference point in trinagle (first vertex)
    r_dcooedge1 = r_coov(:,2) - r_coov(:,1)
    r_dcooedge3 = r_coov(:,1) - r_coov(:,3)
    r_dcoopoint = r_coop      - r_coov(:,1)

!--- compute dot products and squared area of parallelogramm spanned by the two edges
    r_dotg1g1 = DOT_PRODUCT(r_dcooedge1, r_dcooedge1)
    r_dotg3g3 = DOT_PRODUCT(r_dcooedge3, r_dcooedge3)
    r_dotg1g3 = DOT_PRODUCT(r_dcooedge1, r_dcooedge3)
    r_dotg1po = DOT_PRODUCT(r_dcooedge1, r_dcoopoint)
    r_dotg3po = DOT_PRODUCT(r_dcooedge3, r_dcoopoint)
    r_areasq  = r_dotg1g3*r_dotg1g3 - r_dotg1g1*r_dotg3g3

!--- compute barycentric coordinates
    r_sigma = r_dotg3po*r_dotg1g3
    r_sigma = r_sigma - (r_dotg1po*r_dotg3g3)
    r_sigma = r_sigma / r_areasq

    r_tau = r_dotg3po*r_dotg1g1
    r_tau = r_tau - (r_dotg1po*r_dotg1g3)
    r_tau = r_tau / r_areasq

    r_bary(1) = 1.0_GRID_SR - r_sigma - r_tau
    r_bary(2) = r_sigma
    r_bary(3) = r_tau

    DO i_cnt = 1,3
      IF (ABS(r_bary(i_cnt)) <= GRID_EPS) r_bary(i_cnt) = 0.0_GRID_SR
    END DO

  END FUNCTION bary

!*******************************************************************************
! DESCRIPTION of [SUBROUTINE rstoab]:
!> @brief maps from (r,s) coordinates on right triangle to (a,b) coordinates on
!>        square [-1,1]^2 for evaluation of Jacobi polynomials
!
!> @param         r_r
!> @param         r_s
!> @param         r_a
!> @param         r_b
!
!> @note adapted from Hesthaven and Warburton, 2008
!
  SUBROUTINE rstoab(r_r, r_s, r_a, r_b)

    REAL (KIND = GRID_SR), DIMENSION(:), INTENT(in)       :: r_r, r_s
    REAL (KIND = GRID_SR), DIMENSION(:), INTENT(out)      :: r_a, r_b

    r_a = MERGE(2.0_GRID_SR*(1.0_GRID_SR+r_r)/(1.0_GRID_SR-r_s) - 1.0_GRID_SR, &
                -1.0_GRID_SR, ABS(r_s-1.0_GRID_SR) > GRID_EPS)
    r_b = r_s

  END SUBROUTINE rstoab

!*******************************************************************************
! DESCRIPTION of [SUBROUTINE JacobiP]:
!> @brief evaluates the Jacobi polynomial of type (alpha, beta), alpha, beta > -1,
!>        (alpha+beta) != -1, at points x for order i_degree and returns P(x). The
!>        polynomials are normalized to be orthonormal.
!
!> @param         r_x
!> @param         i_alpha
!> @param         i_beta
!> @param         i_degree
!
!> @note adapted from Hesthaven and Warburton, 2008
!
  SUBROUTINE JacobiP(r_x, i_alpha, i_beta, i_degree, r_P)

    REAL (KIND = GRID_SR), DIMENSION(:), INTENT(in)       :: r_x
    INTEGER (KIND = GRID_SI), INTENT(in)                  :: i_alpha, i_beta
    INTEGER (KIND = GRID_SI), INTENT(in)                  :: i_degree
    REAL (KIND = GRID_SR), DIMENSION(:), INTENT(out)      :: r_P

    INTEGER (KIND = GRID_SI)                              :: i_cnt
    REAL (KIND = GRID_SR)                                 :: r_alpha, r_beta, r_gamma0, r_gamma1, &
                                                             r_h1, r_aold, r_anew, r_bnew, r_cnt
    REAL (KIND = GRID_SR), DIMENSION(:,:), ALLOCATABLE    :: r_PL

    r_alpha = REAL(i_alpha, GRID_SR)
    r_beta  = REAL(i_beta , GRID_SR)

    ALLOCATE(r_PL(i_degree+1, SIZE(r_x)))

!--- initial values P_0(r_x) and P_1(r_x)
    r_gamma0  = 2.0_GRID_SR**(i_alpha+i_beta+1) / (r_alpha+r_beta+1.0_GRID_SR) * &
                GAMMA(r_alpha+1.0_GRID_SR) * GAMMA(r_beta+1.0_GRID_SR) / GAMMA(r_alpha+r_beta+1.0_GRID_SR)
    r_PL(1,:) = 1.0_GRID_SR / SQRT(r_gamma0)

    IF (i_degree == 0) THEN
      r_P = r_PL(1,:)
      DEALLOCATE(r_PL)
      RETURN
    END IF

    r_gamma1  = (r_alpha+1.0_GRID_SR) * (r_beta+1.0_GRID_SR) / (r_alpha+r_beta+3.0_GRID_SR) * r_gamma0
    r_PL(2,:) = ((r_alpha+r_beta+2.0_GRID_SR)*r_x/2.0_GRID_SR + (r_alpha-r_beta)/2.0_GRID_SR) / SQRT(r_gamma1)

    IF (i_degree == 1) THEN
      r_P = r_PL(2,:)
      DEALLOCATE(r_PL)
      RETURN
    END IF

!--- repeat value in recurrence
    r_aold = 2.0_GRID_SR / (r_alpha+r_beta+2.0_GRID_SR) * &
            SQRT((r_alpha+1.0_GRID_SR)*(r_beta+1.0_GRID_SR) / (r_alpha+r_beta+3.0_GRID_SR))

!--- forward recurrence using the symmetry of the recurrence
    DO i_cnt = 2,i_degree
      r_cnt = REAL(i_cnt-1, GRID_SR)
      r_h1   = 2.0_GRID_SR*r_cnt + r_alpha + r_beta
      r_anew = 2.0_GRID_SR/(r_h1+2.0_GRID_SR) * &
               SQRT((r_cnt+1.0_GRID_SR) * (r_cnt+1.0_GRID_SR+r_alpha+r_beta) * (r_cnt+1.0_GRID_SR+r_alpha) * &
                    (r_cnt+1.0_GRID_SR+r_beta) / ((r_h1+1.0_GRID_SR)*(r_h1 + 3.0_GRID_SR)))
      r_bnew          = - (r_alpha**2 - r_beta**2) / (r_h1*(r_h1+2.0_GRID_SR))
      r_PL(i_cnt+1,:) = 1.0_GRID_SR/r_anew * (-r_aold*r_PL(i_cnt-1,:) + (r_x-r_bnew)*r_PL(i_cnt,:))
      r_aold          = r_anew
    END DO

    r_P = r_PL(i_degree+1,:)

    DEALLOCATE(r_PL)

  END SUBROUTINE JacobiP

!*******************************************************************************
! DESCRIPTION of [SUBROUTINE Simplex2DP]:
!> @brief evaluates 2D orthonormal polynomial of order (i,j) on simplex at
!>        points (r,s).
!
!> @param         r_r
!> @param         r_s
!> @param         i
!> @param         j
!> @param         r_P
!
!> @note adapted from Hesthaven and Warburton, 2008
!
  SUBROUTINE Simplex2DP(r_r, r_s, i, j, r_P)

    REAL (KIND = GRID_SR), DIMENSION(:), INTENT(in)       :: r_r, r_s
    INTEGER (KIND = GRID_SI), INTENT(in)                  :: i, j
    REAL (KIND = GRID_SR), DIMENSION(:), INTENT(out)      :: r_P

    INTEGER (KIND = GRID_SI)                              :: i_size
    REAL (KIND = GRID_SR), DIMENSION(:), ALLOCATABLE      :: r_a, r_b, r_h1, r_h2

    i_size = SIZE(r_r)
    ALLOCATE(r_a(i_size), r_b(i_size), r_h1(i_size), r_h2(i_size))

    CALL rstoab(r_r, r_s, r_a, r_b)
    CALL JacobiP(r_a, 0    , 0, i, r_h1)
    CALL JacobiP(r_b, 2*i+1, 0, j, r_h2)
    r_P = SQRT(2.0_GRID_SR) * r_h1 * r_h2 * (1.0_GRID_SR-r_b)**i

    DEALLOCATE(r_a, r_b, r_h1, r_h2)

  END SUBROUTINE Simplex2DP

!*******************************************************************************
! DESCRIPTION of [SUBROUTINE Vandermonde2D]:
!> @brief computes the 2D Vandermonde matrix, V_{ij} = phi_j(r_i, s_i)
!
!> @param         i_degree
!> @param         r_bary
!> @param         r_V
!
!> @note adapted from Hesthaven and Warburton, 2008
!
  SUBROUTINE Vandermonde2D(i_degree, r_bary, r_V)

    IMPLICIT NONE

    INTEGER (KIND = GRID_SI), INTENT(in)                  :: i_degree
    REAL (KIND = GRID_SR), DIMENSION(:,:), INTENT(in)     :: r_bary
    REAL (KIND = GRID_SR), DIMENSION(:,:), INTENT(out)    :: r_V

    INTEGER (KIND = GRID_SI)                              :: i_cnt, j_cnt, i_sk
    REAL (KIND = GRID_SR), DIMENSION(SIZE(r_bary,2))      :: r_r, r_s

    r_V = 0.0_GRID_SR

!--- compute coordinates in straight sided reference triangle from barycentric ones
    DO i_cnt = 1, SIZE(r_r)
      r_r(i_cnt) = -r_bary(1,i_cnt) + r_bary(2,i_cnt) - r_bary(3,i_cnt)
      r_s(i_cnt) = -r_bary(1,i_cnt) - r_bary(2,i_cnt) + r_bary(3,i_cnt)
    END DO

!--- build the Vandermonde matrix
    i_sk = 1
    DO i_cnt=0, i_degree
      DO j_cnt=0, i_degree-i_cnt
        CALL Simplex2DP(r_r, r_s, i_cnt, j_cnt, r_V(:,i_sk))
        i_sk = i_sk+1
      END DO
    END DO

  END SUBROUTINE Vandermonde2D

!*******************************************************************************
END MODULE DG_utils
