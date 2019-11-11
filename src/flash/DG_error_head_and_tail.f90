!*******************************************************************************
!
!> @file  DG_error_head_and_tail.f90
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
    REAL (KIND = GRID_SR)                                 :: r_u_l, r_u_r, &
      r_h_l, r_h_r, r_ssh_star, r_hu_star, r_u_star, r_a_l, r_a_r, r_a_star, &
      r_q_r, r_head, r_tail, r_shock, r_bary
    REAL (KIND = GRID_SR), DIMENSION(:,:), ALLOCATABLE    :: r_coo

!--- define constants
    !--- Initial conditions repeated
    r_u_l = 0.0
    r_u_r = 0.0
    r_h_l = 3.0
    r_h_r = 1.0

    !--- Newton Raphson State (comes from python)
    r_ssh_star = 1.84857660309
    r_hu_star  = 4.31264029647
    r_u_star   = r_hu_star/r_ssh_star

    !--- Computational speed
    r_a_l    = SQRT(GRID_GRAV*r_h_l)
    r_a_r    = SQRT(GRID_GRAV*r_h_r)
    r_a_star = SQRT(GRID_GRAV*r_ssh_star)

    r_q_r = SQRT(0.5*((r_ssh_star+r_h_r)*r_ssh_star/r_h_r**2))

    !--- Speed of the waves
    r_head  = r_u_l-r_a_l
    r_tail  = r_u_star-r_a_star
    r_shock = r_u_r + r_a_r*r_q_r

!--- allocate auxiliary array
    ALLOCATE(r_coo(GRID_dimension, i_numdofs), stat=i_alct)
    IF(i_alct /= 0) THEN
      CALL grid_error(c_error='[dg_errorest]: Could not allocate auxiliary variables')
    END IF

!--- get values from grid
    CALL grid_getinfo(p_ghand, i_femtype=FEM_DG, l_finelevel=.TRUE., &
                      r_dofcoordinates=r_coo)

!--- loop through all elements
    elmt_loop: DO i_elmt=1, i_numelmt

      !-- Set refinement indicator to one at head, tail and shock and zero otherwise
      r_bary = SUM(r_coo(1,i_elmtdofs(:,i_elmt)))/3.0_GRID_SR

      r_errorloc(i_elmt) = (1.0/ABS(r_bary - 5.0-p_timestepinfo%r_modeltime*r_shock))/2.0_GRID_SR + &
                           (1.0/ABS(r_bary - 5.0-p_timestepinfo%r_modeltime*r_tail ))/3.0_GRID_SR

    END DO elmt_loop

    IF(GRID_parameters%iolog > 0) &
      write(GRID_parameters%iolog,*) 'ERROREST INFO: max. error: ', MAXVAL(r_errorloc), &
                                     ' min. error: ', MINVAL(r_errorloc)

!--- deallocate workspace
    DEALLOCATE(r_coo)

  END SUBROUTINE dg_errorest
!*******************************************************************************
END MODULE DG_errorestimate
