!*******************************************************************************
!
!> @file  DG_time_utils.F90
!> @brief contains module DG_time_utils
!
!*******************************************************************************
!
! VERSION(S):
!  1. original version                               m. baensch     02/2019
!
!*******************************************************************************
! MODULE DESCRIPTION:
!> @brief routines for timestepping schemes (e.g. GMRES)
!
MODULE DG_time_utils

  USE GRID_api
  USE DG_equation, ONLY : i_nprogvars
  USE DG_flux,     ONLY : fvm_flux
  USE IO_equation, ONLY : p_equationparam
  USE FLASH_parameters

  PRIVATE

  PUBLIC :: GMRES, JFNK_LHS, preconditioner

  CONTAINS

!*******************************************************************************
! DESCRIPTION of [SUBROUTINE GMRES]:
!> @brief GMRES for solving linear system (see for example Habil. thesis of Birken 2012)
!>
!> @param[in]       r_solguess      solution guess (usually Q from last time step)
!> @param[in]       r_residual      residual for r_solguess
!> @param[in]       r_tolerance     (optional) termination tolerance; default = mach. prec.
!> @param[in]       i_maxiterations (optional) number of max. iterations; default = 100
!> @param[out]      r_solution      calculated GMRES solution
!
  SUBROUTINE GMRES(r_solution, r_solguess, r_residual, r_dt, r_gamma, &
                   r_flux_old, r_Q, r_Src, r_time, i_numelmt, i_numedge, &
                   i_faceunknowns, r_metrics_inv, r_Dxi, r_Deta, r_epsi, r_gpsi, &
                   r_eMinvpsi, r_eMinvdpsidxi, r_eMinvdpsideta, r_gMinvpsi, &
                   i_equadpts, i_gquadpts, r_eqwei, r_gqwei, r_elmtvolume, &
                   r_edgelength, r_normals, i_edgeinfo, i_edgeboundary, &
                   i_edgenodes, i_elmtdofs, r_coonod, r_coodof, &
                   i_niterations, r_tolerance, i_maxiterations)

    IMPLICIT NONE

    REAL (KIND = GRID_SR), DIMENSION(:),      INTENT(OUT)  :: r_solution
    REAL (KIND = GRID_SR), DIMENSION(:,:),    INTENT(IN)   :: r_flux_old, r_Q, r_Src, r_Dxi, r_Deta, &
                                                              r_epsi, r_gpsi, r_eMinvpsi, &
                                                              r_eMinvdpsidxi, r_eMinvdpsideta, &
                                                              r_gMinvpsi, r_normals, r_coonod, &
                                                              r_coodof
    REAL (KIND = GRID_SR),                    INTENT(IN)   :: r_time, r_dt, r_gamma
    INTEGER (KIND = GRID_SI),                 INTENT(IN)   :: i_numelmt, i_numedge, &
                                                              i_faceunknowns, &
                                                              i_equadpts, i_gquadpts
    REAL (KIND = GRID_SR), DIMENSION(:,:,:),  INTENT(IN)   :: r_metrics_inv
    REAL (KIND = GRID_SR), DIMENSION(:),      INTENT(IN)   :: r_solguess, r_residual, r_eqwei, r_gqwei, r_elmtvolume, r_edgelength
    INTEGER (KIND = GRID_SI), DIMENSION(:,:), INTENT(IN)   :: i_edgeinfo, i_edgenodes, i_elmtdofs
    INTEGER (KIND = GRID_SI), DIMENSION(:),   INTENT(IN)   :: i_edgeboundary
    INTEGER (KIND = GRID_SI),                 INTENT(OUT)  :: i_niterations

    REAL (KIND = GRID_SR),                    INTENT(IN), OPTIONAL  :: r_tolerance
    INTEGER (KIND = GRID_SI),                 INTENT(IN), OPTIONAL  :: i_maxiterations

!--- local declarations
    REAL (KIND = GRID_SR), DIMENSION(:,:), ALLOCATABLE  :: r_v, r_w, r_h
    REAL (KIND = GRID_SR), DIMENSION(:),   ALLOCATABLE  :: r_a, r_c, r_s, r_g
    REAL (KIND = GRID_SR)                               :: r_tol, r_tmp, r_b
    INTEGER (KIND = GRID_SI)                            :: i_maxiter, i_alct, i_iter, i_cnt, i_cnt2, i_numdofs
    LOGICAL                                             :: l_done = .FALSE.

!--- check for optional input parameters
    IF(PRESENT(r_tolerance)) THEN
      r_tol = r_tolerance
    ELSE
      r_tol = SQRT(GRID_EPS)
    END IF
    IF(PRESENT(i_maxiterations)) THEN
      i_maxiter = i_maxiterations
    ELSE
      i_maxiter = 50
    END IF

    i_numdofs = i_numelmt * i_faceunknowns

!--- allocate work space
!--- the variables correspond to the following in the typical pseudo code tables:
!    r_h = h,  r_a = alpha,  r_c = c,  r_s = s, r_g = gamma

    ALLOCATE(r_v(i_numdofs*i_nprogvars, i_maxiter), &
             r_w(i_numdofs*i_nprogvars, i_maxiter), &
             r_h(i_maxiter+1, i_maxiter), r_a(i_maxiter), &
             r_c(i_maxiter), r_s(i_maxiter), r_g(i_maxiter), stat=i_alct)
    IF(i_alct /= 0) CALL grid_error(c_error='[GMRES]: Could not allocate auxiliary variables')

!--- calculate first basis vector (r_v)
    r_g(1)    = NORM2(r_residual)
    r_v(:,1)  = r_residual / r_g(1)

!--- iteration
    iteration_loop: DO i_iter = 1, i_maxiter

!--- calculate r_w
      r_w(:,i_iter) = JFNK_LHS(r_v(:,i_iter), r_dt, r_gamma, &
                               r_flux_old, r_Q, r_Src, r_time, i_numelmt, i_numedge, i_faceunknowns, &
                               r_metrics_inv, r_Dxi, r_Deta, r_epsi, r_gpsi, &
                               r_eMinvpsi, r_eMinvdpsidxi, r_eMinvdpsideta, r_gMinvpsi, &
                               i_equadpts, i_gquadpts, r_eqwei, r_gqwei, r_elmtvolume, &
                               r_edgelength, r_normals, i_edgeinfo, i_edgeboundary, &
                               i_edgenodes, i_elmtdofs, r_coonod, r_coodof)

!--- Arnoldi iteration
      DO i_cnt = 1, i_iter
        r_h(i_cnt, i_iter) = DOT_PRODUCT(r_v(:,i_cnt), r_w(:,i_iter))
        r_w(:, i_iter)     = r_w(:, i_iter) - r_h(i_cnt, i_iter)*r_v(:,i_cnt)
      END DO

!--- exit if r_h for last iteration has been calculated
      IF (i_iter == i_maxiter) THEN
        r_v(:,1) = r_v(:,i_maxiter)
        r_g(1) = r_g(i_maxiter)
        EXIT
      END IF

!--- check for iteration stop and calculate next basis vector if not
      r_h(i_iter+1, i_iter) = NORM2(r_w(:,i_iter))
      IF (ABS(r_h(i_iter+1, i_iter)) < r_tol) THEN
        l_done = .TRUE.
      ELSE
        r_v(:,i_iter+1) = r_w(:,i_iter) / r_h(i_iter+1, i_iter)
      END IF

!--- Given's rotation (for next column of Hessenberg Matrix r_h)
      DO i_cnt = 1, i_iter-1
        r_tmp = r_h(i_cnt, i_iter)
        r_h(i_cnt,   i_iter) = +r_c(i_cnt)*r_tmp + r_s(i_cnt)*r_h(i_cnt+1, i_iter)
        r_h(i_cnt+1, i_iter) = -r_s(i_cnt)*r_tmp + r_c(i_cnt)*r_h(i_cnt+1, i_iter)
      END DO

      r_b                 = SQRT(r_h(i_iter, i_iter)**2 + r_h(i_iter+1, i_iter)**2)
      r_s(i_iter)         = r_h(i_iter+1, i_iter) / r_b
      r_c(i_iter)         = r_h(i_iter,   i_iter) / r_b
      r_h(i_iter, i_iter) = r_b
      r_g(i_iter+1)       = -r_s(i_iter) * r_g(i_iter)
      r_g(i_iter)         = r_c(i_iter) * r_g(i_iter)

      IF (l_done .OR. ABS(r_g(i_iter+1)) < r_tol) EXIT
    END DO iteration_loop

    i_niterations = i_iter

!--- calculate new solution
    r_a = 0.0_GRID_SR
    DO i_cnt = i_iter, 1, -1
      r_a(i_cnt) = (r_g(i_cnt) - DOT_PRODUCT(r_h(i_cnt,:),r_a)) / r_h(i_cnt, i_cnt)
    END DO

    r_solution = r_solguess + MATMUL(r_v, r_a)

    DEALLOCATE(r_v, r_w, r_h, r_a, r_c, r_s, r_g)

  END SUBROUTINE GMRES

!*******************************************************************************
! DESCRIPTION of [FUNCTION GMRES_LHS]:
!> @brief Calculate LHS for basis vector calculation
!
!
  FUNCTION JFNK_LHS(r_v, r_dt, r_gamma, &
                    r_flux_old, r_Q, r_S, r_time, i_numelmt, i_numedge, &
                    i_faceunknowns, r_metrics_inv, r_Dxi, r_Deta, r_epsi, r_gpsi, &
                    r_eMinvpsi, r_eMinvdpsidxi, r_eMinvdpsideta, r_gMinvpsi, &
                    i_equadpts, i_gquadpts, r_eqwei, r_gqwei, r_elmtvolume, &
                    r_edgelength, r_normals, i_edgeinfo, i_edgeboundary, &
                    i_edgenodes, i_elmtdofs, r_coonod, r_coodof) RESULT(r_w)

    IMPLICIT NONE

    REAL (KIND = GRID_SR), DIMENSION(:),      INTENT(IN)  :: r_v
    REAL (KIND = GRID_SR), DIMENSION(:,:),    INTENT(IN)  :: r_flux_old, r_Q, r_S, r_Dxi, r_Deta, &
                                                             r_epsi, r_gpsi, r_eMinvpsi, &
                                                             r_eMinvdpsidxi, r_eMinvdpsideta, &
                                                             r_gMinvpsi, r_normals, r_coonod, &
                                                             r_coodof
    REAL (KIND = GRID_SR),                    INTENT(IN)  :: r_time, r_dt, r_gamma
    INTEGER (KIND = GRID_SI),                 INTENT(IN)  :: i_numelmt, i_numedge, &
                                                             i_faceunknowns, &
                                                             i_equadpts, i_gquadpts
    REAL (KIND = GRID_SR), DIMENSION(:,:,:),  INTENT(IN)  :: r_metrics_inv
    REAL (KIND = GRID_SR), DIMENSION(:),      INTENT(IN)  :: r_eqwei, r_gqwei, r_elmtvolume, r_edgelength
    INTEGER (KIND = GRID_SI), DIMENSION(:,:), INTENT(IN)  :: i_edgeinfo, i_edgenodes, i_elmtdofs
    INTEGER (KIND = GRID_SI), DIMENSION(:),   INTENT(IN)  :: i_edgeboundary

!--- local declarations
    REAL (KIND = GRID_SR), DIMENSION(i_nprogvars,i_numelmt*i_faceunknowns)  :: r_flux_plus, r_flux_minus
    REAL (KIND = GRID_SR), DIMENSION(SIZE(r_v))                             :: r_w
    REAL (KIND = GRID_SR)                                                   :: r_eps

    IF (NORM2(r_v) < GRID_EPS) THEN
      r_w = 0.0_GRID_SR
    ELSE

!       r_eps = SQRT(GRID_EPS) / NORM2(r_v) * 2.0_GRID_SR
      r_eps = SQRT(1.0_GRID_SR + NORM2(RESHAPE(r_Q, SHAPE(r_v)))) / &
              NORM2(r_v) * (GRID_EPS / 2.0_GRID_SR)**(1.0_GRID_SR/3.0_GRID_SR)

      CALL fvm_flux(r_flux_plus, r_Q + r_eps * RESHAPE(r_v, SHAPE(r_Q)), &
                    r_S, r_time, i_numelmt, i_numedge, &
                    i_faceunknowns, r_metrics_inv, r_Dxi, r_Deta, r_epsi, r_gpsi, &
                    r_eMinvpsi, r_eMinvdpsidxi, r_eMinvdpsideta, r_gMinvpsi, &
                    i_equadpts, i_gquadpts, r_eqwei, r_gqwei, r_elmtvolume, &
                    r_edgelength, r_normals, i_edgeinfo, i_edgeboundary, &
                    i_edgenodes, i_elmtdofs, r_coonod, r_coodof)
      CALL fvm_flux(r_flux_minus, r_Q - r_eps * RESHAPE(r_v, SHAPE(r_Q)), &
                    r_S, r_time, i_numelmt, i_numedge, &
                    i_faceunknowns, r_metrics_inv, r_Dxi, r_Deta, r_epsi, r_gpsi, &
                    r_eMinvpsi, r_eMinvdpsidxi, r_eMinvdpsideta, r_gMinvpsi, &
                    i_equadpts, i_gquadpts, r_eqwei, r_gqwei, r_elmtvolume, &
                    r_edgelength, r_normals, i_edgeinfo, i_edgeboundary, &
                    i_edgenodes, i_elmtdofs, r_coonod, r_coodof)
      r_w = r_v - RESHAPE(r_flux_plus - r_flux_minus, SHAPE(r_v)) / r_eps / 2.0_GRID_SR * r_dt * r_gamma

    END IF

  END FUNCTION JFNK_LHS

!*******************************************************************************
! DESCRIPTION of [SUBROUTINE preconditioner]:
!> @brief performs a preconditioning of the linear system (using ILU(0))
!
!> @param[in]     r_flux            DIMENSION(i_nprogvars,i_numelmt*i_faceunknowns);
!> @param[in]     r_Q               DIMENSION(i_nprogvars,i_numelmt*i_faceunknowns);
!> @param         i_numelmt         total number of elements
!> @param         i_faceunknowns    number of DOFs per element
!> @param         r_metrics_inv     metric terms for each element; DIMENSION(2,2,i_numelmt)
!> @param         r_Dxi             DIMENSION(i_faceunknowns, i_faceunknowns);
!>                                  differentiated basis functions in x direction at Lagrange points
!> @param         r_Deta            DIMENSION(i_faceunknowns, i_faceunknowns);
!>                                  differentiated basis functions in y direction at Lagrange points
!> @param         i_elmtdofs        DOF indices corresponding to each element; DIMENSION(i_faceunknowns, i_numelmt)
!
!   SUBROUTINE preconditioner(r_flux, r_Q, i_numelmt, i_faceunknowns, r_metrics_inv, &
!                             r_Dxi, r_Deta, i_elmtdofs)
!
!     IMPLICIT NONE
!
!     REAL (KIND = GRID_SR), DIMENSION(:,:),    INTENT(IN)  :: r_flux, r_Q, r_Dxi, r_Deta
!     REAL (KIND = GRID_SR), DIMENSION(:,:,:),  INTENT(IN)  :: r_metrics_inv
!     INTEGER (KIND = GRID_SI),                 INTENT(IN)  :: i_numelmt, i_faceunknowns
!     INTEGER (KIND = GRID_SI), DIMENSION(:,:), INTENT(IN)  :: i_elmtdofs
!
! !--- local declarations
!     REAL (KIND = GRID_SR), DIMENSION(i_faceunknowns, i_faceunknowns)  :: r_ddx, r_ddy, r_Jac
!     REAL (KIND = GRID_SR), DIMENSION(i_nprogvars, i_faceunknowns)     :: r_dQdx, r_dQdy, r_ddQ
!
!     r_ddx = r_Dxi*r_metrics_inv(1,1,1) + r_Deta*r_metrics_inv(1,2,1)
!     r_ddy = r_Dxi*r_metrics_inv(2,1,1) + r_Deta*r_metrics_inv(2,2,1)
!
!     r_dQdx = MATMUL(r_Q(:,i_elmtdofs(:,1)), r_ddx)
!     r_dQdy = MATMUL(r_Q(:,i_elmtdofs(:,1)), r_ddy)
!
!     IF (NORM2(r_dQdx) < GRID_EPS) THEN
!       r_ddQ = 0.0_GRID_SR
!     ELSE
!       r_ddQ = MATMUL(1.0_GRID_SR/r_dQdx, r_ddx) + MATMUL(1.0_GRID_SR/r_dQdy, r_ddy)
!     END IF
!
!     r_Jac = MATMUL(TRANSPOSE(r_ddQ), r_flux(:,i_elmtdofs(:,1)))
!     write(*,*) r_Jac
!
!   END SUBROUTINE preconditioner

!*******************************************************************************
END MODULE DG_time_utils
