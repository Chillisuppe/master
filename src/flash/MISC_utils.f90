!*******************************************************************************
!
!> @file  MISC_utils.F90
!> @brief contains module MISC_utils
!
!*******************************************************************************
!
! VERSION(S):
! 1. original version   S. Vater    06/2016
!
!*******************************************************************************
! MODULE DESCRIPTION:
!> @brief some general utility functions
!
MODULE MISC_utils

  USE GRID_api

  PRIVATE
  PUBLIC  :: fsync, compute_inverse

!--- Declare the interface for POSIX fsync function
  INTERFACE
    FUNCTION fsync(fd) bind(c,name="fsync")
    USE iso_c_binding, ONLY: c_int
      INTEGER(c_int), value             :: fd
      INTEGER(c_int)                    :: fsync
    END FUNCTION fsync
  END INTERFACE

  CONTAINS
!*******************************************************************************
! DESCRIPTION of [SUBROUTINE compute_inverse]:
!> @brief computes the inverse of a given matrix
!>
!> @param[in]       r_A         matrix A to be inverted
!> @param[out]      r_invA      inverse of matrix A
!>
  SUBROUTINE compute_inverse(r_A, r_invA)

    IMPLICIT NONE

    REAL (KIND = GRID_SR), DIMENSION(:,:), INTENT(in)         :: r_A
    REAL (KIND = GRID_SR), DIMENSION(:,:), INTENT(out)        :: r_invA

!--- local declarations

    INTEGER (KIND = GRID_SI)                                  :: i_sizeA, i_info, &
      i_alct, i_icnt, i_jcnt
    INTEGER, ALLOCATABLE, DIMENSION(:)                        :: i_IPIV
    REAL (KIND = GRID_SR), ALLOCATABLE, DIMENSION(:)          :: r_work
    REAL (KIND = GRID_SR), ALLOCATABLE, DIMENSION(:,:)        :: r_Atmp

    i_sizeA = SIZE(r_A, 1)
    ALLOCATE(i_IPIV(i_sizeA), r_work(i_sizeA), r_Atmp(i_sizeA,i_sizeA), stat=i_alct)
    IF(i_alct /= 0) CALL grid_error(c_error='[compute_inverse]: could not allocate data')

    DO i_jcnt = 1, i_sizeA
      DO i_icnt = 1, i_sizeA
        r_Atmp(i_icnt, i_jcnt) = r_A(i_icnt, i_jcnt)
      END DO
    END DO

    CALL DGETRF(i_sizeA, i_sizeA, r_Atmp, i_sizeA, i_IPIV, i_info)
    IF(i_info /= 0) CALL grid_error(c_error='[compute_inverse]: function ZGETRF failed')

    CALL DGETRI(i_sizeA, r_Atmp, i_sizeA, i_IPIV, r_work, i_sizeA, i_info)
    IF(i_info /= 0) CALL grid_error(c_error='[compute_inverse]: function ZGETRI failed')

    DO i_jcnt = 1, i_sizeA
      DO i_icnt = 1, i_sizeA
        r_invA(i_icnt, i_jcnt) = r_Atmp(i_icnt, i_jcnt)
      END DO
    END DO

    DEALLOCATE(i_IPIV, r_work, r_Atmp, stat=i_alct)
    IF(i_alct /= 0) CALL grid_error(c_error='[compute_inverse]: could not deallocate data')

  END SUBROUTINE compute_inverse

!*******************************************************************************
END MODULE MISC_utils
