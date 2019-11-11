!*******************************************************************************
!
!> @file  DG_cfl_full.f90
!> @brief contains module DG_cfl
!
!*******************************************************************************
! MODULE DESCRIPTION:
!> @brief computation of CFL number
!
MODULE DG_cfl

  USE GRID_api
  USE FLASH_parameters
  USE DG_equation, ONLY : velocity
  USE IO_equation, ONLY : p_equationparam, p_rteqninfo

  PRIVATE
  PUBLIC  :: compute_cfl

  CONTAINS

!*******************************************************************************
! DESCRIPTION of [SUBROUTINE compute_cfl]:
!> @brief computation of CFL number according to local signal speeds
!
!> @param[in]     p_ghand         grid handling data structure
!> @param[in]     r_Q             state vector with all dofs
!> @param[in]     r_S             vector of source terms with all dofs
!> @param[in]     r_hmin          grid parameter for local cell size
!> @param[in]     i_faceunknowns  DOFs per element
!> @param[in]     i_elmtnodes     element-node relation
!> @param[in]     i_elmtdofs      element-DOF relation
!> @param[in,out] p_tsinfo        runtime information structure
!> @param[out]    r_cfllocal      local CFL numbers
!
  SUBROUTINE compute_cfl(p_ghand, r_Q, r_S, r_hmin, i_faceunknowns, i_elmtnodes, &
                         i_elmtdofs, p_tsinfo, r_cfllocal)

    IMPLICIT NONE

    TYPE(grid_handle), INTENT(in)                         :: p_ghand
    REAL (KIND = GRID_SR), DIMENSION(:,:), INTENT(in)     :: r_Q, r_S
    REAL (KIND = GRID_SR), DIMENSION(:), INTENT(in)       :: r_hmin
    INTEGER (KIND = GRID_SI), INTENT(in)                  :: i_faceunknowns
    INTEGER (KIND = GRID_SI), DIMENSION(:,:), INTENT(in)  :: i_elmtnodes, i_elmtdofs
    TYPE (rt_info), INTENT(inout)                         :: p_tsinfo
    REAL (KIND = GRID_SR), DIMENSION(:), INTENT(out)      :: r_cfllocal

    INTEGER (KIND = GRID_SI)                              :: i_elmt, i_alct
    REAL (KIND = GRID_SR), DIMENSION(i_faceunknowns)      :: r_uabs
    REAL (KIND = GRID_SR), DIMENSION(:), ALLOCATABLE      :: r_aux

!--- allocate workspace
    ALLOCATE(r_aux(SIZE(r_Q,2)), stat=i_alct)
    IF(i_alct /= 0) THEN
      CALL grid_error(c_error='[compute_cfl]: could not allocate workspace')
    END IF

    p_rteqninfo%r_velmax = 0.0_GRID_SR

    DO i_elmt = 1, p_ghand%i_enumfine
!--- compute velocity
      r_uabs = SQRT(velocity(r_Q(1,i_elmtdofs(:,i_elmt)), r_Q(2,i_elmtdofs(:,i_elmt)))**2 + &
                    velocity(r_Q(1,i_elmtdofs(:,i_elmt)), r_Q(3,i_elmtdofs(:,i_elmt)))**2)
      p_rteqninfo%r_velmax = MAX(p_rteqninfo%r_velmax, MAXVAL(r_uabs))

!--- estimate cfl = max(U+\sqrt(\phi))*\Delta t / \Delta x based on values at DOF points
      r_cfllocal(i_elmt) = MAXVAL(r_uabs + SQRT(r_Q(1,i_elmtdofs(:,i_elmt)))) * &
        p_tsinfo%r_timestep / MINVAL(r_hmin(i_elmtnodes(:,i_elmt)))

    END DO

!--- update timestep info structure
    p_tsinfo%r_cflnumber = MAXVAL(r_cfllocal)
    r_aux = MERGE((r_Q(1,:)+r_S(1,:)) / GRID_GRAV - p_equationparam%r_depth, &
                  r_Q(1,:) / GRID_GRAV, r_S(1,:) / GRID_GRAV < p_equationparam%r_depth)
    p_rteqninfo%r_sshmax = MAXVAL(r_aux)
    p_rteqninfo%r_sshmin = MINVAL(r_aux)


!--- deallocate workspace
    DEALLOCATE(r_aux)

  END SUBROUTINE compute_cfl

!*******************************************************************************
END MODULE DG_cfl
