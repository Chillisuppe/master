!*******************************************************************************
!
!> @file  DG_limiter_xingshu.F90
!> @brief contains module DG_limiter
!
!*******************************************************************************
! MODULE DESCRIPTION:
!> @brief function to limit discrete DG solution
!
MODULE DG_limiter

  USE GRID_api
  USE FLASH_parameters
  USE DG_utils
  USE DG_limiter_utils

  PRIVATE
  PUBLIC  :: limiter

  CONTAINS

!*******************************************************************************
! DESCRIPTION of [SUBROUTINE limiter]:
!> @brief flux limiter for DG method
!>
!> @param[in]       p_ghand             grid handling data structure
!> @param[in,out]   r_Q                 discrete solution to be limited;
!>                                      DIMENSION(3,i_numelmt*i_faceunknowns)
!> @param[in]       i_elementedges      edge indices of the edges of each element; DIMENSION(3,i_numelmt)
!> @param[in]       i_edgeinfo          edge-element relation for each edge: local edge index and global element index for adjacent elements; DIMENSION(4,i_numedge)
!> @param[in]       r_edgelength        edge length; DIMENSION(i_numedge)
!> @param[in]       r_metrics_inv       metric terms for each element; DIMENSION(2,2,i_numelmt)
!> @param[in]       i_numelmt           total number of elements
!> @param[in]       r_coonod            node coordinates
!> @param[in]       r_coodof            DOF coordinates
!> @param[in]       i_edgenodes         global node indices for each edge; DIMENSION(2,i_numedge)
!> @param[in]       i_elementnodes      node indices of each element; DIMENSION(3, i_numelmt)
!> @param[in]       i_elementdofs       DOF indices corresponding to each element; DIMENSION(i_faceunknowns, i_numelmt)
!> @param[in]       i_numedge           total number of edges
!> @param[in]       r_bathy             bathymetry of DOFs; DIMENSION(i_numelmt*i_faceunknowns)
!> @param[in]       i_faceunknowns      number of DOFs per element
!> @param[in]       r_evander
!> @param[in]       r_einvvander
!> @param[in]       i_limiter
!
  SUBROUTINE limiter(p_ghand, r_Q, i_elementedges, i_edgeinfo, r_edgelength, r_metrics_inv, &
                     i_numelmt, r_coonod, r_coodof, i_edgenodes, i_elementnodes, i_elementdofs, &
                     i_numedge, r_bathy, i_faceunknowns, r_evander, r_einvvander, i_limiter)

    IMPLICIT NONE

    TYPE (grid_handle), INTENT(in)                            :: p_ghand
    REAL(KIND = GRID_SR), DIMENSION(:,:), INTENT(inout)       :: r_Q
    INTEGER(KIND = GRID_SI), DIMENSION(:,:), INTENT(in)       :: i_elementedges, i_edgeinfo, i_edgenodes, &
                                                                 i_elementnodes, i_elementdofs
    REAL(KIND = GRID_SR), DIMENSION(:), INTENT(in)            :: r_edgelength, r_bathy
    REAL(KIND = GRID_SR), DIMENSION(:,:,:), INTENT(in)        :: r_metrics_inv
    REAL(KIND = GRID_SR), DIMENSION(:,:), INTENT(in)          :: r_coonod, r_coodof, r_evander, r_einvvander
    INTEGER(KIND = GRID_SI), INTENT(in)                       :: i_numelmt, i_faceunknowns, &
                                                                 i_numedge, i_limiter

    REAL(KIND = GRID_SR), DIMENSION(:), ALLOCATABLE           :: r_h, r_h_int, r_u, r_v, r_bat, &
                                                                 r_elmt_mean, r_elmt_mean2, r_mean
    REAL(KIND = GRID_SR), DIMENSION(2,2)                      :: r_matrix
    REAL(KIND = GRID_SR), DIMENSION(:,:), ALLOCATABLE         :: r_epsi, r_b, r_basis, r_m, r_vander_m
    REAL(KIND = GRID_SR)                                      :: r_min, r_heavi, r_emean, &
                                                                 r_emean2, r_emeank
    REAL(KIND = GRID_SR), DIMENSION(i_faceunknowns)           :: r_lambdapi_phi, r_vec
    REAL(KIND = GRID_SR), DIMENSION(3)                        :: r_p1, r_lambdapi_p1, r_diff

    INTEGER(KIND = GRID_SI)                                   :: i_elmt, i_cnt, i_elmt_k, &
                                                                 i_alct, k, n, &
                                                                 i_degree, i_loop, i_edge2, i_loc, i_loclmt, &
                                                                 i_edge, i_elmt2
    INTEGER(KIND = GRID_SI), DIMENSION(3)                     :: i_ind_vec
    REAL(KIND = GRID_SR), DIMENSION(i_faceunknowns)           :: r_dhdy, r_dhdx
    REAL(KIND = GRID_SR), DIMENSION(i_faceunknowns, i_faceunknowns) :: r_Dxi, r_Deta, r_ddx, r_ddy

!--- allocate workspace
    ALLOCATE(r_h(i_faceunknowns), r_h_int(i_faceunknowns), &
             r_epsi(i_faceunknowns, i_faceunknowns), &
             r_u(i_faceunknowns), r_v(i_faceunknowns), &
             r_bat(i_faceunknowns), &
             r_vander_m(i_faceunknowns, i_faceunknowns), &
             r_m(i_numedge,2), &
             r_b(i_numelmt,2), &
             r_basis(3, 3), &
             r_elmt_mean(i_numelmt), &
             r_mean(i_numelmt), &
             r_elmt_mean2(i_numelmt), stat=i_alct)
    IF (i_alct .NE. 0) CALL grid_error(55)

!  Get info from signature and initialize
    i_degree  = GRID_femtypes%p_type(FEM_DG)%sig%i_degree
    r_Dxi     = GRID_femtypes%p_type(FEM_DG)%sig%r_dpsidxi
    r_Deta    = GRID_femtypes%p_type(FEM_DG)%sig%r_dpsideta
    r_epsi    = GRID_femtypes%p_type(FEM_DG)%sig%r_epsiquad

    r_b  = 0._GRID_SR
    r_m  = 0._GRID_SR

! Preparation for further computations:
! Compute mean values for every element
    DO i_elmt=1, i_numelmt
      r_elmt_mean(i_elmt)  = dot_product(r_evander(1:i_faceunknowns,1),r_Q(i_elmt,1:i_faceunknowns,1))/sum(r_evander(:,1))
      r_elmt_mean2(i_elmt) = dot_product(r_evander(1:i_faceunknowns,1),r_Q(i_elmt,1:i_faceunknowns,1) + &
                             r_bathy(i_elmt,1:i_faceunknowns))/sum(r_evander(:,1))
      r_b(i_elmt,1) = (r_coonod(1,i_elementnodes(1,i_elmt))+r_coonod(1,i_elementnodes(2,i_elmt))+&
                       r_coonod(1,i_elementnodes(3,i_elmt)))/3
      r_b(i_elmt,2) = (r_coonod(2,i_elementnodes(1,i_elmt))+r_coonod(2,i_elementnodes(2,i_elmt))+&
                       r_coonod(2,i_elementnodes(3,i_elmt)))/3
    END DO !i_elmt

! Compute centres of edges
    DO i_edge=1, i_numedge
! Midpoint of current edge
      r_m(i_edge,1)= 0.5*(r_coonod(1,i_edgenodes(1,i_edge))+r_coonod(1,i_edgenodes(2,i_edge)))
      r_m(i_edge,2)= 0.5*(r_coonod(2,i_edgenodes(1,i_edge))+r_coonod(2,i_edgenodes(2,i_edge)))
    END DO !i_edge

!-- Initialize Lagrangebasis at edgemidpoints
    CALL assemble_basis(r_basis, 1_GRID_SI)

! The limiter element loop
    DO i_elmt=1, i_numelmt

! Initialize Parameter m used for the pp limiting
      r_min         = 0._GRID_SR
      r_lambdapi_p1 = 0._GRID_SR
      r_lambdapi_phi= 0._GRID_SR
      r_p1          = 0._GRID_SR

! Interpolate ssh on Gauss Lobatto points
      r_h     = r_Q(i_elmt,:,1)
      r_h_int = matmul(r_epsi,r_h)
      r_bat   = r_bathy(i_elmt,:)

! Compute point of min ssh for the first time
      r_min=MINVAL(r_h_int(:))

! Now apply the TVB limiter according to Cockburn and Shu; indicator based on h+b if m gt 0 or h if not.
      r_heavi = heaviside(r_min)       !-- Will return 1 for r_min = 0 and 0 otherwise
      r_vec   = r_h+(1-r_heavi)*r_bat
      r_mean  = r_heavi*r_elmt_mean+(1-r_heavi)*r_elmt_mean2

! (i) Compute P1 part of approximation
      CALL P1part(r_vec, r_p1, i_faceunknowns, r_evander, r_einvvander)

! (ii) Compute limited version of (i)
      CALL LambdaPI(i_elmt, i_faceunknowns, i_degree, i_edgeinfo, i_elementedges, &
                    r_edgelength, r_lambdapi_p1, r_p1, r_b(i_elmt,:), r_b, r_m, &
                    r_mean, r_basis)

! (iii) If (ii) .NE. (i) then (ii) is limited version of approximation
      r_diff = abs(r_lambdapi_p1-r_p1)

      IF (sum(r_diff(:)).LE. 10E-12) THEN
        r_lambdapi_phi = r_h
!write(*,*) 'OPT 1', r_lambdapi_phi
!-- Else: Limiting has to be done on phi
      ELSEIF(r_heavi .EQ. 1._GRID_SR)THEN
        CALL P1back(r_lambdapi_phi, r_lambdapi_p1, i_faceunknowns, i_degree)

!write(*,*) 'OPT 2', r_lambdapi_phi
      ELSE
        CALL P1part(r_h, r_p1, i_faceunknowns, r_evander, r_einvvander)
        CALL LambdaPI(i_elmt, i_faceunknowns, i_degree, i_edgeinfo, i_elementedges, &
                      r_edgelength, r_lambdapi_p1, r_p1, r_b(i_elmt,:), r_b, r_m, &
                      r_elmt_mean, r_basis)
        CALL P1back(r_lambdapi_phi, r_lambdapi_p1, i_faceunknowns, i_degree)

!write(*,*) 'OPT 3', r_lambdapi_phi
      END IF

      r_Q(i_elmt,:,1) = r_lambdapi_phi

    END DO !i_elmt

!PP LIMITER

    DO i_elmt=1, i_numelmt

! assemble derivative operators (Note: This might be taken into an external function in the future)
      r_ddx(:,:) = 0
      r_ddy(:,:) = 0
      DO i_cnt = 1, i_faceunknowns
        r_ddx(:,i_cnt) = r_Dxi(:,i_cnt)*r_metrics_inv(1,1,i_elmt) + r_Deta(:,i_cnt)*r_metrics_inv(1,2,i_elmt)
        r_ddy(:,i_cnt) = r_Dxi(:,i_cnt)*r_metrics_inv(2,1,i_elmt) + r_Deta(:,i_cnt)*r_metrics_inv(2,2,i_elmt)
      END DO
      r_ddx = MATMUL(r_ddx, r_epsi)
      r_ddy = MATMUL(r_ddy, r_epsi)

!> Interpol on lgl points
      r_h = r_Q(i_elmt,:,1)
      r_u = r_Q(i_elmt,:,2)
      r_v = r_Q(i_elmt,:,3)
      r_h_int = matmul(r_epsi,r_h)

!> Check if pp-limiting is needed; compute min again
      r_min=MINVAL(r_h_int(:))

      IF (r_min .GE. 0) THEN  !.GT.
        CALL pp_limiter(p_ghand, r_h, r_h_int, r_u, r_v, r_bat, r_min, i_faceunknowns, r_evander, r_epsi)

!> Put info back to phi (still interpolated on quadpoints)
        r_Q(i_elmt,:,1)=r_h
        r_Q(i_elmt,:,2)=r_u
        r_Q(i_elmt,:,3)=r_v
      END IF !m .gt. 0

      r_dhdx         = MATMUL(r_h, r_ddx)
      r_dhdy         = MATMUL(r_h, r_ddy)


    END DO !i_elmt

!> Deallocate workspace
    DEALLOCATE(r_h, r_h_int, r_epsi, r_u, r_v, r_bat, r_b, &
               r_vander_m, r_elmt_mean, r_mean, r_elmt_mean2, r_basis)

  END SUBROUTINE limiter

!*******************************************************************************
! DESCRIPTION of [SUBROUTINE pp_limiter]:
!> @brief positivity preserving limiter
!>
!> @param[in]     p_ghand         grid handling data structure
!> @param[in,out] r_h             DIMENSION(i_faceunknowns);
!>                                sea surface height at DOFs of element
!> @param[in,out] r_h_int         DIMENSION(i_faceunknowns);
!>                                sea surface height at quadrature points of element
!> @param[in,out] r_hu            DIMENSION(i_faceunknowns);
!>                                velocity in x-direction at DOFs of element
!> @param[in,out] r_hv            DIMENSION(i_faceunkonwns);
!>                                velocity in y-direction at DOFs of element
!> @param         r_bathy         DIMENSION(i_faceunknowns);
!>                                bathymetry at DOFs of element
!> @param         r_m             REAL; mitpoint of edge
!> @param         i_faceunknowns  number of DOFs per element
!> @param         r_evander       DIMENSION(i_faceunknowns, i_faceunknowns);
!>                                vandermonde matrix
!> @param         r_epsi          DIMENSION(i_faceunknowns, i_faceunknowns);
!>                                interpolation matrix onto quadrature points inside the element
!
  SUBROUTINE pp_limiter(p_ghand, r_h, r_h_int, r_hu, r_hv, r_bathy, r_m, i_faceunknowns, r_evander, r_epsi)

    IMPLICIT NONE

    TYPE(grid_handle), INTENT(in)                     :: p_ghand

    REAL(KIND=GRID_SR), DIMENSION(:), INTENT(inout)   :: r_h, r_hu, r_hv, r_h_int
    REAL(KIND=GRID_SR), DIMENSION(:)                  :: r_bathy
    REAL(KIND=GRID_SR), DIMENSION(:), ALLOCATABLE     :: r_height
    REAL(KIND=GRID_SR), DIMENSION(:,:), ALLOCATABLE   :: r_tilde
    REAL(KIND=GRID_SR), DIMENSION(:,:)                :: r_evander, r_epsi
    REAL(KIND=GRID_SR)                                :: r_m, r_theta, r_h_mean, r_height_mean, &
                                                         r_hu_mean, r_hv_mean, r_h_rec, &
                                                         r_elmt_mean, r_heavi

    INTEGER(KIND=GRID_SI)                               :: i_faceunknowns, i_cnt, i_alct

!> Allocate workspace
    ALLOCATE(r_height(i_faceunknowns), r_tilde(3,i_faceunknowns), stat=i_alct)
    IF (i_alct .NE. 0) CALL grid_error(55)

!> Initialize constants; compute height
    r_tilde       = 0._GRID_SR
    r_height      = r_h + r_bathy

!> Compute mean values for prognostic variables
    r_h_mean      = DOT_PRODUCT(r_evander(:,1)/sum(r_evander(:,1)),r_h)
    r_height_mean = DOT_PRODUCT(r_evander(:,1)/sum(r_evander(:,1)),r_height)
    r_hu_mean     = DOT_PRODUCT(r_evander(:,1)/sum(r_evander(:,1)),r_hu)
    r_hv_mean     = DOT_PRODUCT(r_evander(:,1)/sum(r_evander(:,1)),r_hv)

!> Work with either height or ssh mean depending on whether it's wet or dry
    r_heavi = heaviside(r_m)
    r_elmt_mean = r_h_mean !r_height_mean-r_heavi*(r_height_mean-r_h_mean)

!> Compute weighting parameter $/theta$:
    IF(abs(r_h_mean-r_m) .LE. 10E-30) THEN
      r_m=r_h_mean
    END IF

    r_theta=DMIN1(1._GRID_SR, r_h_mean/(r_h_mean-r_m))


!> If theta = 1 no limiting has to be done, otherwise a linear scaling is performed

    IF (r_theta .NE. 1) THEN !Do the actual limiting
      DO i_cnt=1, i_faceunknowns
        r_tilde(1,i_cnt)=r_theta*(r_h(i_cnt)-r_h_mean)+r_h_mean
        r_tilde(2,i_cnt)= r_theta*(r_hu(i_cnt)-r_hu_mean) + r_hu_mean
        r_tilde(3,i_cnt)= r_theta*(r_hv(i_cnt)-r_hv_mean) + r_hv_mean

!> Reconstruct ssh, u, v
        r_h(i_cnt)  = r_tilde(1,i_cnt)
        r_hu(i_cnt) = r_tilde(2,i_cnt)
        r_hv(i_cnt) = r_tilde(3,i_cnt)
      END DO !i_cnt

    END IF !theta not 1

    IF (r_m .LT. 0) write(*,*) 'Das sollte nicht passieren, wenn alles nass ist. :)'

!--- Deallocate workspace
    DEALLOCATE(r_height, r_tilde)

  END SUBROUTINE pp_limiter

!*******************************************************************************
END MODULE DG_limiter
