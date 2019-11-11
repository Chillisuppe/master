!*******************************************************************************
!
!> @file  DG_error_vort.f90
!> @brief contains module DG_errorestimate
!
!*******************************************************************************
! MODULE DESCRIPTION:
!> @brief
!
MODULE DG_errorestimate

  USE GRID_api
  USE FLASH_parameters
  USE DG_equation

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
    REAL (KIND = GRID_SR)                                 :: r_det_inv, r_vort
    REAL (KIND = GRID_SR), DIMENSION(2,2)                 :: r_Jacobi, r_metrics
    REAL (KIND = GRID_SR), DIMENSION(i_faceunknowns)      :: r_x, r_y, r_u, r_v, r_h
    REAL (KIND = GRID_SR), DIMENSION(:,:), ALLOCATABLE    :: r_coo

!--- allocate auxiliary array
    ALLOCATE(r_coo(GRID_dimension, i_numdofs), stat=i_alct)

    IF(i_alct /= 0) THEN
      CALL grid_error(c_error='[dg_errorest]: Could not allocate auxiliary variables')
    END IF

!--- get values
    CALL grid_getinfo(p_ghand, i_femtype=FEM_DG, l_finelevel=.TRUE., &
                      r_dofcoordinates=r_coo)

!--- loop through all elements
    elmt_loop: DO i_elmt=1, i_numelmt

      r_x = r_coo(1, i_elmtdofs(:,i_elmt))
      r_y = r_coo(2, i_elmtdofs(:,i_elmt))
      r_h = r_Q(1, i_elmtdofs(:,i_elmt))
      r_u = velocity(r_Q(1,i_elmtdofs(:,i_elmt)), r_Q(2,i_elmtdofs(:,i_elmt)))
      r_v = velocity(r_Q(1,i_elmtdofs(:,i_elmt)), r_Q(3,i_elmtdofs(:,i_elmt)))

!--- metric terms
      r_metrics(1,1) = (r_x(2)-r_x(1))/2
      r_metrics(1,2) = (r_x(3)-r_x(1))/2
      r_metrics(2,1) = (r_y(2)-r_y(1))/2
      r_metrics(2,2) = (r_y(3)-r_y(1))/2

!--- this comes from Ax+b=xi \Rightarrow dxi/dx is this matrix
      r_Jacobi(1,1) =  r_metrics(2,2)
      r_Jacobi(2,1) = -r_metrics(1,2)
      r_Jacobi(1,2) = -r_metrics(2,1)
      r_Jacobi(2,2) =  r_metrics(1,1)

      r_det_inv     = r_metrics(1,1)*r_metrics(2,2) - r_metrics(2,1)*r_metrics(1,2)
      r_Jacobi(:,:) = r_Jacobi(:,:) / r_det_inv

      CALL vort(r_u, r_v, r_Jacobi, i_faceunknowns, r_vort)

!--- elementwise gradient on dofs
      r_errorloc(i_elmt) = r_vort

    END DO elmt_loop

!--- deallocate workspace
    DEALLOCATE(r_coo)

    IF(GRID_parameters%iolog > 0) &
      write(GRID_parameters%iolog,*) 'ERROREST INFO: max. error: ', MAXVAL(r_errorloc), &
                                     ' min. error: ', MINVAL(r_errorloc)

  END SUBROUTINE dg_errorest

!*******************************************************************************
! DESCRIPTION of [FUNCTION vort]:
!> @brief
!
!> @param         r_u
!> @param         r_v
!> @param         r_Jacobi
!> @param         i_faceunknowns    number of DOFs per element
!> @param         r_vort
!
  SUBROUTINE vort(r_u, r_v, r_Jacobi, i_faceunknowns, r_vort)

    IMPLICIT NONE

    REAL (KIND = GRID_SR), DIMENSION(:), INTENT(in)                 :: r_u, r_v
    REAL (KIND = GRID_SR), DIMENSION(:,:), INTENT(in)               :: r_Jacobi
    INTEGER (KIND = GRID_SI), INTENT(in)                            :: i_faceunknowns
    REAL (KIND = GRID_SR), INTENT(out)                              :: r_vort

!--- local declarations
    REAL (KIND = GRID_SR), DIMENSION(i_faceunknowns,i_faceunknowns) :: r_epsi, r_Dxi, r_Deta
    REAL (KIND = GRID_SR), DIMENSION(i_faceunknowns,i_faceunknowns) :: r_Dxi_in, r_Deta_in
    REAL (KIND = GRID_SR), DIMENSION(i_faceunknowns,2)              :: r_grad_u, r_grad_v
    REAL (KIND = GRID_SR), DIMENSION(i_faceunknowns)                :: r_u_xi, r_v_xi, r_u_eta, &
                                                                       r_v_eta
!--- initialize constants
    r_Dxi   = GRID_femtypes%p_type(FEM_DG)%sig%r_dpsidxi
    r_Deta  = GRID_femtypes%p_type(FEM_DG)%sig%r_dpsideta
    r_epsi  = GRID_femtypes%p_type(FEM_DG)%sig%r_epsiquad

!--- intepolate Differential-Matrix
    r_Dxi_in  = MATMUL(r_Dxi, r_epsi)
    r_Deta_in = MATMUL(r_Deta, r_epsi)

!--- interpolate to inner dofs
    r_u_xi  = MATMUL(r_Dxi_in, r_u)
    r_u_eta = MATMUL(r_Deta_in, r_u)
    r_v_xi  = MATMUL(r_Dxi_in, r_v)
    r_v_eta = MATMUL(r_Deta_in, r_v)

!--- Construct gradient operator in physical variables
    r_grad_u(:,1) = r_Jacobi(1,1)*r_u_xi + r_Jacobi(1,2)*r_u_eta
    r_grad_u(:,2) = r_Jacobi(2,1)*r_u_xi + r_Jacobi(2,2)*r_u_eta
    r_grad_v(:,1) = r_Jacobi(1,1)*r_v_xi + r_Jacobi(1,2)*r_v_eta
    r_grad_v(:,2) = r_Jacobi(2,1)*r_v_xi + r_Jacobi(2,2)*r_v_eta

!-- Compute vorticity
    r_vort = r_grad_v(1,1)-r_grad_u(1,2)

  END SUBROUTINE vort

!*******************************************************************************
END MODULE DG_errorestimate
