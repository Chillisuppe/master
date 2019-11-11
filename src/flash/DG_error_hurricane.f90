!*******************************************************************************
!
!> @file  DG_error_hurricane.f90
!> @brief contains module DG_errorestimate
!
!*******************************************************************************
! MODULE DESCRIPTION:
!> @brief
!
MODULE DG_errorestimate

  USE GRID_api
  USE FLASH_parameters
  USE DG_equation, ONLY: VAR1D_WINDX, VAR1D_WINDY

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
    REAL (KIND = GRID_SR), DIMENSION(:,:), INTENT(in)     :: r_Q, r_S,  r_evander
    REAL (KIND = GRID_SR), DIMENSION(:), INTENT(in)       :: r_elmtwidth
    INTEGER (KIND = GRID_SI), DIMENSION(:,:), INTENT(in)  :: i_elmtdofs, i_edgeinfo
    REAL (KIND = GRID_SR), DIMENSION(:,:,:), INTENT(in)   :: r_metrics_inv

!--- local declarations
    INTEGER (KIND = GRID_SI)                              :: i_alct, i_elmt
    INTEGER (KIND = GRID_SI), DIMENSION(2)                :: i_valind
    REAL (KIND = GRID_SR), DIMENSION(:,:), ALLOCATABLE    :: r_windfield

!--- allocate auxiliary array
    ALLOCATE(r_windfield(GRID_dimension,i_numdofs), stat=i_alct)
    IF(i_alct /= 0) THEN
      CALL grid_error(c_error='[dg_errorest]: Could not allocate auxiliary variables')
    END IF

!--- get values
    i_valind = (/ VAR1D_WINDX, VAR1D_WINDY /)
    CALL grid_getinfo(p_ghand, i_femtype=FEM_DG, l_finelevel=.TRUE., &
                      i_arraypoint=i_valind, r_dofvalues=r_windfield)

!--- loop through all elements
    elmt_loop: DO i_elmt=1, i_numelmt
      r_errorloc(i_elmt) = MAXVAL(r_windfield(1,i_elmtdofs(:,i_elmt))) + &
                           MAXVAL(r_windfield(2,i_elmtdofs(:,i_elmt)))
    END DO elmt_loop

    IF(GRID_parameters%iolog > 0) &
      write(GRID_parameters%iolog,*) 'ERROREST INFO: max. error: ', MAXVAL(r_errorloc), &
                                     ' min. error: ', MINVAL(r_errorloc)

!--- deallocate workspace
    DEALLOCATE(r_windfield)

  END SUBROUTINE dg_errorest
!*******************************************************************************
END MODULE DG_errorestimate
