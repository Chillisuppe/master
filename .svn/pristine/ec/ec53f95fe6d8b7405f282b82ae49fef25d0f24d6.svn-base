!*******************************************************************************
!
!> @file  DG_limiter_utils.F90
!> @brief contains module DG_limiter_utils
!
!> @author Nicole Beisiegel
!
!*******************************************************************************
! MODULE DESCRIPTION:
!> @brief
!
MODULE DG_limiter_utils

  USE GRID_api
  USE FLASH_parameters

  PRIVATE
  PUBLIC :: TVB_minmod, assemble_basis, P1part, LambdaPI, P1back

  CONTAINS

!*******************************************************************************
! DESCRIPTION of [FUNCTION TVB_minmod]:
!> @brief
!
!> @param         r_phi1
!> @param         r_phi2
!> @param         r_M
!> @param         r_delta
!
  FUNCTION TVB_minmod(r_phi1, r_phi2, r_M, r_delta)

    IMPLICIT NONE

    REAL(KIND = GRID_SR)        :: r_phi1, r_phi2, r_delta, r_M
    REAL(KIND = GRID_SR)        :: TVB_minmod

    IF (abs(r_phi1) < r_M*r_delta**2) THEN
      TVB_minmod = r_phi1
    ELSEIF (SIGN(1._GRID_SR, r_phi1) == SIGN(1._GRID_SR, r_phi2)) THEN
      TVB_minmod = SIGN(1._GRID_SR, r_phi1) * MIN(abs(r_phi1), abs(r_phi2))
    ELSE
      TVB_minmod = 0
    END IF
    RETURN

  END FUNCTION TVB_minmod

!*******************************************************************************
! DESCRIPTION of [SUBROUTINE LambdaPI]:
!> @brief
!>
!> @param[in]       i_elmt
!> @param[in]       i_faceunknowns    number of DOFs per element
!> @param[in]       i_degree          polynomial degree of element
!> @param[in]       i_edgeinfo        edge-element relation for each edge: local edge index and global element index for adjacent elements; DIMENSION(4,i_numedge)
!> @param[in]       i_elementedges    edge indices of the edges of each element; DIMENSION(3,i_numelmt)
!> @param[in]       r_edgelength
!> @param[in,out]   r_lambdapi_h
!> @param[in]       r_h
!> @param[in]       r_b0
!> @param[in]       r_b
!> @param[in]       r_m
!> @param[in]       r_elmt_mean
!> @param[in]       r_basis
!
  SUBROUTINE LambdaPI(i_elmt, i_faceunknowns, i_degree, &
                      i_edgeinfo, i_elementedges, r_edgelength, r_lambdapi_h, r_h, r_b0, r_b, r_m, &
                      r_elmt_mean, r_basis)

    IMPLICIT NONE

    INTEGER (KIND = GRID_SI), INTENT(in)                  :: i_faceunknowns, i_degree, i_elmt
    INTEGER (KIND = GRID_SI)                              :: k, i_edge, i_elmt_k, n, i_loop, i_elmt2,&
                                                             i_edge2, i_cnt

    INTEGER (KIND = GRID_SI), DIMENSION(i_faceunknowns)   :: i_ind_vec

    INTEGER (KIND = GRID_SI), DIMENSION(:,:), INTENT(in)  :: i_edgeinfo, i_elementedges

    REAL (KIND = GRID_SR), DIMENSION(:), INTENT(in)       :: r_edgelength
    REAL (KIND = GRID_SR)                                 :: r_det, r_neg, r_pos, r_theta_pl, r_theta_mn, &
                                                             r_deltax
    REAL (KIND = GRID_SR), DIMENSION(2,2)                 :: r_matrix
    REAL (KIND = GRID_SR), DIMENSION(3), INTENT(inout)    :: r_lambdapi_h
    REAL (KIND = GRID_SR), DIMENSION(3), INTENT(in)       :: r_h
    REAL (KIND = GRID_SR), DIMENSION(1,2), INTENT(in)     :: r_b0
    REAL (KIND = GRID_SR), DIMENSION(:,:), INTENT(in)     :: r_basis, r_b, r_m
    REAL (KIND = GRID_SR), DIMENSION(:), INTENT(in)       :: r_elmt_mean
    REAL (KIND = GRID_SR)                                 :: r_nu, r_MM
    REAL (KIND = GRID_SR), DIMENSION(2)                   :: r_b1, r_b2, r_x
    REAL (KIND = GRID_SR), DIMENSION(3)                   :: r_edge_mean, r_Delta, r_alpha1, r_alpha2, &
                                                             r_laplace_phi_bar, r_phi_tilde, r_Delta_hat
    REAL (KIND = GRID_SR), DIMENSION(3,1) :: r_vander_m


! Initialize local quantities
    r_nu = 1.50_GRID_SR         !As in the paper by Cockburn & Shu
    r_MM = 50.0_GRID_SR
    i_ind_vec = 0_GRID_SI

!------- Compute for every edge

    edge_loop: DO k= 1,3

      i_edge = i_elementedges(k,i_elmt)

      IF (i_edgeinfo(3,i_edge) .EQ. i_elmt) THEN
        i_elmt_k = i_edgeinfo(4,i_edge)
        n=1
      ELSE
        i_elmt_k = i_edgeinfo(3,i_edge)
        n=2
      END IF

      r_deltax = r_edgelength(i_edge)

!-------- Vandermondematrix Lagrangepoly(edgemidpoints)
      SELECT CASE(i_edgeinfo(n,i_edge))
        CASE(1)
          r_vander_m(:,1) = (/ 0.0, 0.5, 0.5 /)
        CASE(2)
          r_vander_m(:,1) = (/ 0.5, 0.0, 0.5 /)
        CASE(3)
          r_vander_m(:,1) = (/ 0.5, 0.5, 0.0 /)
      END SELECT

!-- This is always the linear part

      i_ind_vec(k)=(i_edgeinfo(n,i_edge))

!---- Compute r_phi(r_m)

      r_edge_mean(k) = dot_product(r_vander_m(:,1) , r_h)

!---- Compute barycentre of edge-neighbour

      IF((i_elmt_k .GT. 0).AND.(i_elmt_k .NE. i_elmt)) THEN
        r_b1(:) = r_b(i_elmt_k,:)
      ELSE
        r_b1(:) = r_m(i_edge,:)
      END IF

!------- Determine coefficients alpha

!--- Choose second neighbour for basis representation
!take that ones neighbour
      i_loop =0

 100  i_edge2 = i_elementedges(mod(k+i_loop,3)+1,i_elmt)

      IF ((i_edgeinfo(3,i_edge2) .EQ. i_elmt).AND.(i_edgeinfo(4,i_edge2).GT. 0)) THEN
        r_b2 = r_b(i_edgeinfo(4,i_edge2),:)
        i_elmt2 = i_edgeinfo(4,i_edge2)
      ELSEIF ((i_edgeinfo(4,i_edge2) .EQ. i_elmt).AND.(i_edgeinfo(3,i_edge2).GT. 0)) THEN !ist hier edgeinfo nicht immer gt 0?
        r_b2 = r_b(i_edgeinfo(3,i_edge2),:)
        i_elmt2 = i_edgeinfo(3,i_edge2)
      ELSE
        r_b2 = r_m(i_edge2, :)
        i_elmt2 = 0
      END IF

!--- Choose local e_1 and e_2 depending on edgeinfo

      ! m_i - b_0

      r_det=((r_b1(1)-r_b0(1,1))*(r_b2(2)-r_b0(1,2)))-((r_b1(2)-r_b0(1,2))*(r_b2(1)-r_b0(1,1)))

      r_matrix(1,1) = (r_b2(2)-r_b0(1,2))
      r_matrix(1,2) = -(r_b2(1)-r_b0(1,1))
      r_matrix(2,1) = -(r_b1(2)-r_b0(1,2))
      r_matrix(2,2) = (r_b1(1)-r_b0(1,1))

      r_matrix= r_matrix/r_det

!multipliziere mit m_i-b0

      r_x = matmul(r_matrix, r_m(i_edge,:)-r_b0(1,:))
      r_alpha1(k) = r_x(1)
      r_alpha2(k) = r_x(2)

      IF ((r_alpha1(k) .LT. 0 ) .OR. (r_alpha2(k) .LT. 0 )) THEN
        i_loop = i_loop +1
        IF (i_loop .LE. 1)THEN
          GO TO 100
        END IF
      END IF

!------- Determine r_phi_tilde, r_laplace_phi_bar

      r_phi_tilde(k)       = r_edge_mean(k)-r_elmt_mean(i_elmt)
      IF ((i_elmt2 .EQ. 0).AND.(i_elmt_k .EQ. 0)) THEN
        r_laplace_phi_bar(k) = r_alpha1(k)*(r_edge_mean(k)-r_elmt_mean(i_elmt))+ &
                               r_alpha2(k)*(r_edge_mean(mod(k,3)+1)-r_elmt_mean(i_elmt))
      ELSEIF ((i_elmt2 .GT. 0).AND.(i_elmt_k .EQ. 0)) THEN
        r_laplace_phi_bar(k) = r_alpha1(k)*(r_edge_mean(k)-r_elmt_mean(i_elmt))+ &
                               r_alpha2(k)*(r_elmt_mean(i_elmt2 )-r_elmt_mean(i_elmt))
      ELSEIF ((i_elmt2 .EQ. 0).AND.(i_elmt_k .GT. 0)) THEN
        r_laplace_phi_bar(k) = r_alpha1(k)*(r_elmt_mean(i_elmt_k)-r_elmt_mean(i_elmt))+ &
                               r_alpha2(k)*(r_edge_mean(mod(k,3)+1)-r_elmt_mean(i_elmt))
      ELSE
        r_laplace_phi_bar(k) = r_alpha1(k)*(r_elmt_mean(i_elmt_k)-r_elmt_mean(i_elmt))+ &
                               r_alpha2(k)*(r_elmt_mean(i_elmt2 )-r_elmt_mean(i_elmt))
      END IF

!------- Compute Delta_i

      r_Delta(k) = TVB_minmod(r_phi_tilde(k),r_nu*r_laplace_phi_bar(k), r_MM, r_deltax)

    END DO edge_loop

!------- IF sum Delta_i = ...
    IF (sum(r_Delta(:)) .LE. 10E-10) THEN
      r_lambdapi_h = r_elmt_mean(i_elmt) + matmul(r_basis,r_Delta)
    ELSE
      r_pos = 0._GRID_SR
      r_neg = 0._GRID_SR

      DO i_cnt=1,3
        r_pos = r_pos + DMAX1(0._GRID_SR, r_Delta(i_cnt))
        r_neg = r_neg + DMAX1(0._GRID_SR,-r_Delta(i_cnt))
      END DO !i_cnt

      r_theta_pl = DMIN1(1._GRID_SR, r_neg/r_pos)
      r_theta_mn = DMIN1(1._GRID_SR, r_pos/r_neg)

      DO i_cnt=1,3
        r_Delta_hat(i_cnt)= r_theta_pl*DMAX1(0._GRID_SR, r_Delta(i_cnt))-&
                            r_theta_mn*DMAX1(0._GRID_SR,-r_Delta(i_cnt))
      END DO !i_cnt

      r_lambdapi_h = r_elmt_mean(i_elmt) + matmul(r_basis,r_Delta_hat)
    END IF !sum Delta=0

  END SUBROUTINE LambdaPI

!*******************************************************************************
! DESCRIPTION of [SUBROUTINE assemble_basis]:
!> @brief
!>
!> @param[in,out]   r_basis
!> @param           i_degree
!
  SUBROUTINE assemble_basis(r_basis, i_degree)

    IMPLICIT NONE

    INTEGER (KIND = GRID_SI)                              :: i_degree
    REAL (KIND = GRID_SR), DIMENSION(:,:), INTENT(inout)  :: r_basis

!------ Lagrangelinearbasisatmidpnts(dofs)
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

  END SUBROUTINE assemble_basis

!*******************************************************************************
! DESCRIPTION of [SUBROUTINE P1part]:
!> @brief
!>
!> @param[in]       r_phi
!> @param[in,out]   r_p1
!> @param[in]       i_faceunknowns    number of DOFs per element
!
  SUBROUTINE P1part(r_phi, r_p1, i_faceunknowns, r_evander, r_einvvander)

    IMPLICIT NONE

    REAL (KIND = GRID_SR), DIMENSION(i_faceunknowns), INTENT(in)    :: r_phi
    REAL (KIND = GRID_SR), DIMENSION(:), INTENT(inout)              :: r_p1
    INTEGER (KIND = GRID_SI), INTENT(in)                            :: i_faceunknowns
    REAL (KIND = GRID_SR), DIMENSION(:,:), INTENT(in)               :: r_evander, r_einvvander

    INTEGER (KIND = GRID_SI)                                        :: i_cnt, j_cnt
    REAL (KIND = GRID_SR), DIMENSION(i_faceunknowns)                :: r_tmp

!--- initialize arrays
    r_p1 = 0._GRID_SR

    IF (i_faceunknowns == 3_GRID_SI) THEN
      r_p1 = r_phi(1:3)
    ELSE
!-- Compute coefficients for linear basis (Weierstrass) for the first three "modes"
      DO j_cnt=1, i_faceunknowns
        DO i_cnt=4,i_faceunknowns
!           r_evander(j_cnt,i_cnt) = 0.0
        END DO
      END DO
      r_tmp =  MATMUL(r_einvvander,(MATMUL(r_evander,r_phi)))
      r_p1 = r_tmp(1:3)
    END IF

!  r_p1 = r_phi(1:3)

  END SUBROUTINE P1part

!*******************************************************************************
! DESCRIPTION of [SUBROUTINE P1back]:
!> @brief flux limiter for DG method
!>
!> @param[in,out]   r_phi
!> @param[in]       r_p1
!> @param[in]       i_faceunknowns    number of DOFs per element
!> @param[in]       i_degree
!
  SUBROUTINE P1back(r_phi, r_p1, i_faceunknowns, i_degree)

    IMPLICIT NONE

    INTEGER (KIND = GRID_SI), INTENT(in)                            :: i_faceunknowns, &
                                                                       i_degree
    REAL (KIND = GRID_SR), DIMENSION(i_faceunknowns), INTENT(inout) :: r_phi
    REAL (KIND = GRID_SR), DIMENSION(3), INTENT(in)                 :: r_p1
    REAL (KIND = GRID_SR), DIMENSION(i_faceunknowns,3)              :: r_basis

!  CALL assemble_basis(r_basis, i_degree)
!-- Interpolation to basis that is used for computations

   r_phi=0.0
    IF (i_degree .EQ. 1) THEN
      r_phi=r_p1
    ELSE
      r_basis(1,:) = (/ 1.0, 0.0, 0.0 /)
      r_basis(2,:) = (/ 0.0, 1.0, 0.0 /)
      r_basis(3,:) = (/ 0.0, 0.0, 1.0 /)
      r_basis(4,:) = (/ 0.5, 0.5, 0.0 /)
      r_basis(5,:) = (/ 0.0, 0.5, 0.5 /)
      r_basis(6,:) = (/ 0.5, 0.0, 0.5 /)

      r_phi = matmul(r_basis, r_p1)
    END IF

  END SUBROUTINE P1back

!*******************************************************************************
END MODULE DG_limiter_utils
