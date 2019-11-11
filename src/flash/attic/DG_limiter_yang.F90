!*******************************************************************************
!
!> @file  DG_limiter_yang.F90
!> @brief contains module DG_limiter
!
! FUNCTION:
!   Limiter of Yang and Wang (2010)
!
!*******************************************************************************
MODULE DG_limiter
  USE GRID_api
  USE FLASH_parameters
  USE DG_utils

  PRIVATE
!   PUBLIC  :: limiter

  CONTAINS
!*******************************************************************************
! SUBROUTINE limiter(p_ghand, r_phi, i_elementedges, i_edgeinfo, r_edgelength, &
!                      i_numelmt, r_coo, i_edgenodes, i_elementnodes, i_numedge, &
!                      i_approxorder)

! IMPLICIT NONE

! TYPE (grid_handle), INTENT(in)                      :: p_ghand
! INTEGER (KIND=GRID_SI)                              :: i_numelmt, i_elmt, i_approxorder, &
!                                                        i_alct, i_cnt, i_face
! INTEGER (KIND=GRID_SI), DIMENSION(:,:)              ::
! REAL(KIND=GRID_SR), DIMENSION(:)                    :: r_edgelength
! REAL(KIND=GRID_SR), DIMENSION(:,:,:)                :: r_phi
! REAL(KIND=GRID_SR), DIMENSION(:,:,:), ALLOCATABLE   :: r_phi_lim
! REAL(KIND=GRID_SR), DIMENSION(:,:)                  :: r_coo
! REAL (KIND=GRID_SR), DIMENSION(:,:), ALLOCATABLE    :: r_mean_der, r_der
! REAL (KIND=GRID_SR)                                 :: r_beta = 1.5_GRID_SR


!--- Allocate workspace

! ALLOCATE(r_der(i_approxorder+1, i_faceunknowns), &
!          r_mean_der(i_approxorder+1, i_faceunknowns), &
!          r_phi_lim(i_elmt,1,i_faceunknowns), stat=i_alct)

! IF(i_alct .NE. 0) &
!   CALL grid_error(c_error:'[yang_wang]: Could not allocate enough memory.')

!--- Compute element means
!       DO i_elmt=1, i_numelmt
!          DO i_cnt=1, i_faceunknowns
!             r_elmt_mean(i_elmt)= dot_product(r_evander(1:i_faceunknowns,1),r_phi(i_elmt,1:i_faceunknowns,1))&
!                                   /sum(r_evander(:,1))
!          END DO !i_cnt
!        END DO !i_elmt

! 100 DO i_elmt=1, i_numelmt

!--- Compute higher order derivatives of phi in (xi,eta) and do the integration

!     DO i_cnt=1, i_approxorder
!        r_der(i_cnt+1, :) = matmul(r_phi, r_Dxi)
!        r_mean_der(i_cnt+1,:) =
!     END DO !i_cnt

!--- Check if limiting necessary after every step (reverse order)
!     DO i_cnt = 1, i_approxorder

!     r_forward  = r_mean(!x-x_i
!     r_backward =
!     r_minmod = TVB_minmod(r_mean_der(i_approxorder-i_cnt+1,:), r_beta*r_forward, &
!                          r_backward)

!     IF (r_minmod .EQ. r_mean_der(i_approxorder-i_cnt+1,:)) GO TO 100
!
!     END DO !i_cnt

! END DO !i_elmt

!--- Deallocate workspace

! DEALLOCATE(r_mean_der))

! END SUBROUTINE limiter

!*******************************************************************************
END MODULE DG_limiter
