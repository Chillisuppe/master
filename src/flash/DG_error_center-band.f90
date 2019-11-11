!*******************************************************************************
!
!> @file  DG_error_center-band.f90
!> @brief contains module DG_errorestimate
!
!*******************************************************************************
! MODULE DESCRIPTION:
!> @brief
!
MODULE DG_errorestimate

  USE GRID_api
  USE FLASH_parameters

  PRIVATE
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
    INTEGER (KIND = GRID_SI)                              :: i_alct, i_elmt
    REAL (KIND = GRID_SR)                                 :: r_dist
    REAL (KIND = GRID_SR), DIMENSION(2)                   :: r_centr
    REAL (KIND = GRID_SR), DIMENSION(GRID_dimension, i_faceunknowns) :: r_xy
    REAL (KIND = GRID_SR), DIMENSION(:,:), ALLOCATABLE    :: r_coo

!--- allocate auxiliary array
    ALLOCATE(r_coo(GRID_dimension, i_numdofs), stat=i_alct)
    IF(i_alct /= 0) THEN
      CALL grid_error(c_error='[dg_errorest]: Could not allocate auxiliary variables')
    END IF

!--- get coordinates
    CALL grid_getinfo(p_ghand, i_femtype=FEM_DG, l_finelevel=.TRUE., &
                      r_dofcoordinates=r_coo)

!--- position of the vortex
    r_centr = (/ 2.0, 1.0 /)

!--- loop through all elements
    elmt_loop: DO i_elmt=1, i_numelmt

!--- determine distance to position of the vortex
      r_xy   = r_coo(:,i_elmtdofs(:, i_elmt))
      r_dist = MAXVAL(ABS(r_xy(2,:)-r_centr(2)))
      r_errorloc(i_elmt) = 1.0_GRID_SR/r_dist
    END DO elmt_loop

    IF(GRID_parameters%iolog > 0) &
      write(GRID_parameters%iolog,*) 'ERROREST INFO: max. error: ', MAXVAL(r_errorloc), &
                                     ' min. error: ', MINVAL(r_errorloc)

!--- deallocate workspace
    DEALLOCATE(r_coo)

  END SUBROUTINE dg_errorest
!*******************************************************************************
END MODULE DG_errorestimate
