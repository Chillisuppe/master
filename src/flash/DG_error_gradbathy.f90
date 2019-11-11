!*******************************************************************************
!
!> @file  DG_error_gradbathy.f90
!> @brief contains module DG_errorestimate
!
!*******************************************************************************
! MODULE DESCRIPTION:
!> @brief
!
MODULE DG_errorestimate

  USE GRID_api
  USE FLASH_parameters

  PRIVATE :: fvm_error
  PUBLIC  :: dg_errorest

  CONTAINS

!*******************************************************************************
! DESCRIPTION of [SUBROUTINE dg_errorest]:
!> @brief
!
!> @param[in]         p_ghand           grid handling data structure
!> @param[in,out]     r_errorloc        computed local error indicators
!> @param             i_faceunknowns    number of DOFs per element
!> @param             i_numdofs         total number of DOFs
!> @param             i_numelmt         total number of elements
!> @param             i_nnum            total number of nodes
!> @param             i_numedges        total number of edges
!> @param[in]         r_evander         Vandermonde matrix
!> @param[in]         r_elmtwidth
!> @param[in]         i_elmtdofs        global element-dof relation matrix
!> @param             i_edgeinfo        edge-element relation for each edge
!> @param             r_metrics_inv     metric terms for each element
!> @param[in]         r_Q               discrete solution vector
!> @param             r_S               discrete fields used in source terms
!
  SUBROUTINE dg_errorest(p_ghand, r_errorloc, i_faceunknowns, i_numdofs, i_numelmt, &
                         i_nnum, i_numedges, r_evander, r_elmtwidth, i_elmtdofs, &
                         i_edgeinfo, r_metrics_inv, r_Q, r_S)

    IMPLICIT NONE

    TYPE (grid_handle), INTENT(in)                        :: p_ghand
    REAL (KIND = GRID_SR), DIMENSION(:), INTENT(out)      :: r_errorloc
    INTEGER (KIND = GRID_SI), INTENT(in)                  :: i_faceunknowns, i_numdofs, &
      i_numelmt, i_nnum, i_numedges
    REAL (KIND = GRID_SR), DIMENSION(:,:), INTENT(in)     :: r_Q, r_S, r_evander
    REAL (KIND = GRID_SR), DIMENSION(:), INTENT(in)       :: r_elmtwidth
    INTEGER (KIND = GRID_SI), DIMENSION(:,:), INTENT(in)  :: i_elmtdofs, i_edgeinfo
    REAL (KIND = GRID_SR), DIMENSION(:,:,:), INTENT(in)   :: r_metrics_inv

!--- local declarations
    INTEGER (KIND = GRID_SI)                              :: i_elmt
    REAL (KIND = GRID_SR), DIMENSION(i_faceunknowns)      :: r_nod

!--- loop through all elements
    elmt_loop: DO i_elmt=1, i_numelmt
!--- elementwise gradient of bathymetry on dofs
      r_nod = r_S(1,i_elmtdofs(:,i_elmt)) / GRID_GRAV
      r_errorloc(i_elmt) = fvm_error(r_nod)
    END DO elmt_loop

    IF(GRID_parameters%iolog > 0) &
      write(GRID_parameters%iolog,*) 'ERROREST INFO: max. error: ', MAXVAL(r_errorloc), &
                                     ' min. error: ', MINVAL(r_errorloc)

  END SUBROUTINE dg_errorest

!*******************************************************************************
! DESCRIPTION of [FUNCTION fvm_error]:
!> @brief
!
!> @param             r_v
!> @return                      return value
!
  FUNCTION fvm_error(r_v) RESULT (r_esterr)

    IMPLICIT NONE

    REAL (KIND = GRID_SR), DIMENSION(3)                 :: r_v

!--- local declarations
    REAL (KIND = GRID_SR)                               :: r_esterr
    REAL (KIND = GRID_SR), PARAMETER                    :: r_1o3 = (1.0/3.0)
    REAL (KIND = GRID_SR)                               :: r_d1, r_d2, r_d3

!--- calculate differences
    r_d1  =  abs(r_v(1) - r_v(2))
    r_d2  =  abs(r_v(2) - r_v(3))
    r_d3  =  abs(r_v(3) - r_v(1))

!--- this is the estimated error
    r_esterr = (r_d1 + r_d2 + r_d3) * r_1o3

  END FUNCTION fvm_error
!*******************************************************************************
END MODULE DG_errorestimate
