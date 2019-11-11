!*******************************************************************************
!
!> @file  DG_bernstein_utils.F90
!> @brief contains module DG_bernstein_utils
!
!*******************************************************************************
!
! VERSION(s):
! 1. version            nicole beisiegel        02/2012
!
!*******************************************************************************
! MODULE DESCRIPTION:
!> @brief provides some routines for DG method with Bernstein polynomials as
!>        basis function
!
MODULE DG_bernstein_utils

  USE GRID_api
  USE FLASH_parameters, ONLY: FEM_DG

  PRIVATE
  PUBLIC  :: bernstein_quad_interpol_elmt

  CONTAINS

!*******************************************************************************
! DESCRIPTION of [FUNCTION bernstein_quad_interpol_elmt]:
!> @brief
!
!> @param[in,out] r_bernstein
!> @param         r_h
!> @param         i_faceunknowns    number of DOFs per element
!> @param         i_equadpts
!
  SUBROUTINE bernstein_quad_interpol_elmt(r_bernstein, r_h, i_faceunknowns, i_equadpts)

    IMPLICIT NONE

    INTEGER(KIND=GRID_SI)                                   :: i_approxorder, i_knt, &
                                                               i_equadpts, i_cnt, &
                                                               i_r, i_n, i_alct, &
                                                               i_faceunknowns, i_siz
    INTEGER(KIND=GRID_SI), DIMENSION(3)                     :: i_index, i_lambda
    INTEGER(KIND=GRID_SI), DIMENSION(3,3)                   :: i_maxindex
    INTEGER(KIND=GRID_SI), DIMENSION(:,:), ALLOCATABLE      :: i_lambda_vec

    REAL(KIND=GRID_SR), DIMENSION(:), INTENT(inout)         :: r_bernstein
    REAL(KIND=GRID_SR), DIMENSION(i_faceunknowns)           :: r_h
    REAL(KIND=GRID_SR), DIMENSION(:,:), ALLOCATABLE         :: r_bary_quadcoo, &
                                                               r_elmt_dofcoo
    REAL(KIND=GRID_SR), DIMENSION(3)                        :: r_tau, r_vec, r_e
    REAL(KIND=GRID_SR), DIMENSION(:,:,:,:), ALLOCATABLE     :: r_b

!-- Initialize constant
    i_approxorder = GRID_femtypes%p_type(FEM_DG)%sig%i_degree
    r_e=(/1,1,1/)

!-- Allocate workspace
    ALLOCATE(i_lambda_vec(i_faceunknowns ,3), &
             r_b(i_equadpts, i_approxorder +1, i_approxorder+1,i_approxorder+1), &
             r_bary_quadcoo(3,i_equadpts), &
             r_elmt_dofcoo(3, i_faceunknowns), stat=i_alct)
    IF(i_alct /= 0) CALL grid_error(c_error='[bernstein_quad_interpol_elmt]: could not allocate memory')

    r_bary_quadcoo= GRID_femtypes%p_type(FEM_DG)%sig%r_equadcoo
    r_elmt_dofcoo = GRID_femtypes%p_type(FEM_DG)%sig%r_edofcoo

!-- Initialize matrix for multi-index structure; use that first three dofs are the vertices
    i_maxindex = r_elmt_dofcoo(:,1:3)*i_approxorder

!-- Initialize ordinates r_b(0,...)
    DO i_cnt = 1, i_faceunknowns

!-- Barycentric coordinates of i_cnt-th dof
      r_vec = r_elmt_dofcoo(:,i_cnt)
      i_index = INT(MATMUL(i_maxindex, r_vec))+r_e
      r_b(1,i_index(1),i_index(2),i_index(3)) = r_h(i_cnt)

    END DO !i_cnt


!-- De Casteljau iteration for all quadpoints
    DO i_cnt = 1, i_equadpts

!-- Get barycentric coordinates of i_cnt-th quadrature point
      r_tau= r_bary_quadcoo(:,i_cnt)

      DO i_r = 2, i_approxorder+1 !just notation has to start with 2, because Fortran doesnt like zeros

!-- Find the multi-indices lambda depending on r with |\lambda|=i_approxorder-r
        i_n = i_approxorder+1-i_r

        i_lambda_vec(1:3,:)= (i_approxorder-1)*r_elmt_dofcoo(:,1:3)

!-- Number of multi-indices of level r_n
        i_siz = INT((i_n+1)*(i_n+2)*0.5)

        DO i_knt = 1, i_siz

!-- Set local i_lambda to corresponding multi-index
          i_lambda= i_lambda_vec(i_knt,:)+r_e

          r_b(i_r,i_lambda(1), i_lambda(2), i_lambda(3))= &
            r_tau(1)*r_b(i_r-1,i_lambda(1)+1, i_lambda(2), i_lambda(3)) + &
            r_tau(2)*r_b(i_r-1,i_lambda(1), i_lambda(2)+1, i_lambda(3)) + &
            r_tau(3)*r_b(i_r-1,i_lambda(1), i_lambda(2), i_lambda(3)+1)
        END DO !i_knt
      END DO !i_r

!-- Set value of quadpoint = b(n,0,0,0)
      r_bernstein(i_cnt) = r_b(i_approxorder+1,1,1,1)

    END DO !i_cnt

!-- Deallocate workspace
    DEALLOCATE(i_lambda_vec)

  END SUBROUTINE bernstein_quad_interpol_elmt

!*******************************************************************************
! DESCRIPTION of [FUNCTION bernstein_quad_interpol_edge]:
!> @brief
!
!> @param[in,out] r_bernstein
!> @param         r_h
!> @param         i_faceunknowns    number of DOFs per element
!> @param         i_gquadpts        number of quadrature points per edge
!
  SUBROUTINE bernstein_quad_interpol_edge(r_bernstein, r_h, i_faceunknowns, i_gquadpts)

    IMPLICIT NONE

    INTEGER(KIND=GRID_SI)                                   :: i_approxorder, i_knt, &
                                                               i_gquadpts, i_cnt, &
                                                               i_r, i_n, i_alct, &
                                                               i_faceunknowns, i_siz
    INTEGER(KIND=GRID_SI), DIMENSION(3)                     :: i_index, i_lambda
    INTEGER(KIND=GRID_SI), DIMENSION(3,3)                   :: i_maxindex
    INTEGER(KIND=GRID_SI), DIMENSION(:,:), ALLOCATABLE      :: i_lambda_vec

    REAL(KIND=GRID_SR), DIMENSION(:), INTENT(inout)         :: r_bernstein
    REAL(KIND=GRID_SR), DIMENSION(i_faceunknowns)           :: r_h
    REAL(KIND=GRID_SR), DIMENSION(:,:), ALLOCATABLE         :: r_bary_quadcoo, &
                                                               r_elmt_dofcoo
    REAL(KIND=GRID_SR), DIMENSION(3)                        :: r_tau, r_vec, r_e
    REAL(KIND=GRID_SR), DIMENSION(:,:,:,:), ALLOCATABLE     :: r_b

!-- Initialize constant
    i_approxorder = GRID_femtypes%p_type(FEM_DG)%sig%i_degree
    r_e=(/1,1,1/)

!-- Allocate workspace
    ALLOCATE(i_lambda_vec(i_faceunknowns ,3), &
             r_b(i_gquadpts, i_approxorder +1, i_approxorder+1,i_approxorder+1), &
             r_bary_quadcoo(3,i_gquadpts), &
             r_elmt_dofcoo(3, i_faceunknowns), stat=i_alct)
    IF(i_alct /= 0) CALL grid_error(c_error='[bernstein_quad_interpol_elmt]: could not allocate memory')

    r_bary_quadcoo= GRID_femtypes%p_type(FEM_DG)%sig%r_equadcoo
    r_elmt_dofcoo = GRID_femtypes%p_type(FEM_DG)%sig%r_edofcoo

!-- Initialize matrix for multi-index structure; use that first three dofs are the vertices
    i_maxindex = r_elmt_dofcoo(:,1:3)*i_approxorder

!-- Initialize ordinates r_b(0,...)
    DO i_cnt = 1, i_faceunknowns

!-- Barycentric coordinates of i_cnt-th dof
      r_vec = r_elmt_dofcoo(:,i_cnt)
      i_index = INT(MATMUL(i_maxindex, r_vec))+r_e
      r_b(1,i_index(1),i_index(2),i_index(3)) = r_h(i_cnt)

    END DO !i_cnt


!-- De Casteljau iteration for all quadpoints
    DO i_cnt = 1, i_gquadpts

!-- Get barycentric coordinates of i_cnt-th quadrature point
      r_tau= r_bary_quadcoo(:,i_cnt)

      DO i_r = 2, i_approxorder+1 !just notation has to start with 2, because Fortran doesnt like zeros

!-- Find the multi-indices lambda depending on r with |\lambda|=i_approxorder-r
        i_n = i_approxorder+1-i_r

        i_lambda_vec(1:3,:)= (i_approxorder-1)*r_elmt_dofcoo(:,1:3)

!-- Number of multi-indices of level r_n
        i_siz = INT((i_n+1)*(i_n+2)*0.5)

        DO i_knt = 1, i_siz

!-- Set local i_lambda to corresponding multi-index
          i_lambda= i_lambda_vec(i_knt,:)+r_e

          r_b(i_r,i_lambda(1), i_lambda(2), i_lambda(3))= &
            r_tau(1)*r_b(i_r-1,i_lambda(1)+1, i_lambda(2), i_lambda(3)) + &
            r_tau(2)*r_b(i_r-1,i_lambda(1), i_lambda(2)+1, i_lambda(3)) + &
            r_tau(3)*r_b(i_r-1,i_lambda(1), i_lambda(2), i_lambda(3)+1)
        END DO !i_knt
      END DO !i_r

!-- Set value of quadpoint = b(n,0,0,0)
      r_bernstein(i_cnt) = r_b(i_approxorder+1,1,1,1)

    END DO !i_cnt

!-- Deallocate workspace
    DEALLOCATE(i_lambda_vec)

  END SUBROUTINE bernstein_quad_interpol_edge

!*******************************************************************************
END MODULE DG_bernstein_utils
