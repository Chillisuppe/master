!*******************************************************************************
!
!> @file  DG_limiter_BJSV.F90
!> @brief contains module DG_limiter
!
!*******************************************************************************
! MODULE DESCRIPTION:
!> @brief function to limit discrete DG solution
!
MODULE DG_limiter

  USE GRID_api
  USE FLASH_parameters
  USE DG_equation

  PRIVATE
  PUBLIC  :: limiter

  CONTAINS

!*******************************************************************************
! DESCRIPTION of [SUBROUTINE limiter]:
!> @brief flux limiter for DG method
!>
!> Slope limiter of Barth and Jesperson combined with Wetting and Drying fix
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
!> @param[in]       r_S                 source terms at DOFs; DIMENSION(i_nsrcterms, i_numelmt*i_faceunknowns)
!> @param[in]       i_faceunknowns      number of DOFs per element
!> @param[in]       r_evander
!> @param[in]       r_einvvander
!
  SUBROUTINE limiter(p_ghand, r_Q, i_elementedges, i_edgeinfo, r_edgelength, r_metrics_inv, &
                     i_numelmt, r_coonod, r_coodof, i_edgenodes, i_elementnodes, i_elementdofs, &
                     i_numedge, r_S, i_faceunknowns, r_evander, r_einvvander)

    IMPLICIT NONE

    TYPE (grid_handle), INTENT(in)                            :: p_ghand
    REAL(KIND = GRID_SR), DIMENSION(:,:), INTENT(inout)       :: r_Q
    INTEGER(KIND = GRID_SI), DIMENSION(:,:), INTENT(in)       :: i_elementedges, i_edgeinfo, i_edgenodes, &
                                                                 i_elementnodes, i_elementdofs
    REAL(KIND = GRID_SR), DIMENSION(:), INTENT(in)            :: r_edgelength
    REAL(KIND = GRID_SR), DIMENSION(:,:,:), INTENT(in)        :: r_metrics_inv
    REAL(KIND = GRID_SR), DIMENSION(:,:), INTENT(in)          :: r_coonod, r_coodof, r_evander, r_einvvander, r_S
    INTEGER(KIND = GRID_SI), INTENT(in)                       :: i_numelmt, i_faceunknowns, &
                                                                 i_numedge

    INTEGER(KIND = GRID_SI)                                   :: i_alct, i_edge, i_elmt
    REAL(KIND = GRID_SR)                                      :: r_alpha, r_hext, r_utmp
    INTEGER(KIND = GRID_SI), DIMENSION(:), ALLOCATABLE        :: i_dofssorted
    REAL(KIND = GRID_SR), DIMENSION(:), ALLOCATABLE           :: r_corr, r_htmp
    REAL(KIND = GRID_SR), DIMENSION(:,:), ALLOCATABLE         :: r_emean, r_emin, r_emax, &
                                                                 r_ulim, r_Qlim
    LOGICAL, DIMENSION(:), ALLOCATABLE                        :: l_loc

!--- allocate workspace
    ALLOCATE(i_dofssorted(i_faceunknowns), &
             r_corr(i_faceunknowns), &
             r_htmp(i_faceunknowns), &
             r_emean(6, i_numelmt), &
             r_emin(6, i_numelmt), &
             r_emax(6, i_numelmt), &
             r_ulim(i_faceunknowns, 3), &
             r_Qlim(3, i_numelmt*i_faceunknowns), &
             l_loc(i_faceunknowns), stat=i_alct)
    IF (i_alct /= 0) CALL grid_error(c_error='[DG_limiter_BJSV]: Could not allocate enough memory')

!--- compute mean, min and max values for each element
    DO i_elmt=1, i_numelmt
      r_emean(1,i_elmt) = DOT_PRODUCT(r_evander(:,1), r_Q(1,i_elementdofs(:,i_elmt))) / SUM(r_evander(:,1))
      r_emean(2,i_elmt) = DOT_PRODUCT(r_evander(:,1), r_Q(2,i_elementdofs(:,i_elmt))) / SUM(r_evander(:,1))
      r_emean(3,i_elmt) = DOT_PRODUCT(r_evander(:,1), r_Q(3,i_elementdofs(:,i_elmt))) / SUM(r_evander(:,1))
      r_emean(4,i_elmt) = DOT_PRODUCT(r_evander(:,1), r_Q(1,i_elementdofs(:,i_elmt)) + &
                                      r_S(1,i_elementdofs(:,i_elmt))) / SUM(r_evander(:,1))
      r_emean(5,i_elmt) = velocity(r_emean(1,i_elmt), r_emean(2,i_elmt))
      r_emean(6,i_elmt) = velocity(r_emean(1,i_elmt), r_emean(3,i_elmt))

      r_emin(1,i_elmt)  = r_emean(1,i_elmt)
      r_emin(2,i_elmt)  = r_emean(2,i_elmt)
      r_emin(3,i_elmt)  = r_emean(3,i_elmt)
      r_emin(4,i_elmt)  = r_emean(4,i_elmt)
      r_emin(5,i_elmt)  = r_emean(5,i_elmt)
      r_emin(6,i_elmt)  = r_emean(6,i_elmt)

      r_emax(1,i_elmt)  = r_emean(1,i_elmt)
      r_emax(2,i_elmt)  = r_emean(2,i_elmt)
      r_emax(3,i_elmt)  = r_emean(3,i_elmt)
      r_emax(4,i_elmt)  = r_emean(4,i_elmt)
      r_emax(5,i_elmt)  = r_emean(5,i_elmt)
      r_emax(6,i_elmt)  = r_emean(6,i_elmt)
    END DO


!--- loop over all edges to update min and max at inner edges
    DO i_edge=1, i_numedge
      IF ((i_edgeinfo(3,i_edge) > 0) .AND. (i_edgeinfo(4,i_edge) > 0)) THEN
        r_emin(1,i_edgeinfo(3,i_edge)) = MIN(r_emin(1,i_edgeinfo(3,i_edge)), r_emean(1,i_edgeinfo(4,i_edge)))
        r_emin(1,i_edgeinfo(4,i_edge)) = MIN(r_emin(1,i_edgeinfo(4,i_edge)), r_emean(1,i_edgeinfo(3,i_edge)))
        r_emin(2,i_edgeinfo(3,i_edge)) = MIN(r_emin(2,i_edgeinfo(3,i_edge)), r_emean(2,i_edgeinfo(4,i_edge)))
        r_emin(2,i_edgeinfo(4,i_edge)) = MIN(r_emin(2,i_edgeinfo(4,i_edge)), r_emean(2,i_edgeinfo(3,i_edge)))
        r_emin(3,i_edgeinfo(3,i_edge)) = MIN(r_emin(3,i_edgeinfo(3,i_edge)), r_emean(3,i_edgeinfo(4,i_edge)))
        r_emin(3,i_edgeinfo(4,i_edge)) = MIN(r_emin(3,i_edgeinfo(4,i_edge)), r_emean(3,i_edgeinfo(3,i_edge)))
        r_emin(4,i_edgeinfo(3,i_edge)) = MIN(r_emin(4,i_edgeinfo(3,i_edge)), r_emean(4,i_edgeinfo(4,i_edge)))
        r_emin(4,i_edgeinfo(4,i_edge)) = MIN(r_emin(4,i_edgeinfo(4,i_edge)), r_emean(4,i_edgeinfo(3,i_edge)))
        r_emin(5,i_edgeinfo(3,i_edge)) = MIN(r_emin(5,i_edgeinfo(3,i_edge)), r_emean(5,i_edgeinfo(4,i_edge)))
        r_emin(5,i_edgeinfo(4,i_edge)) = MIN(r_emin(5,i_edgeinfo(4,i_edge)), r_emean(5,i_edgeinfo(3,i_edge)))
        r_emin(6,i_edgeinfo(3,i_edge)) = MIN(r_emin(6,i_edgeinfo(3,i_edge)), r_emean(6,i_edgeinfo(4,i_edge)))
        r_emin(6,i_edgeinfo(4,i_edge)) = MIN(r_emin(6,i_edgeinfo(4,i_edge)), r_emean(6,i_edgeinfo(3,i_edge)))

        r_emax(1,i_edgeinfo(3,i_edge)) = MAX(r_emax(1,i_edgeinfo(3,i_edge)), r_emean(1,i_edgeinfo(4,i_edge)))
        r_emax(1,i_edgeinfo(4,i_edge)) = MAX(r_emax(1,i_edgeinfo(4,i_edge)), r_emean(1,i_edgeinfo(3,i_edge)))
        r_emax(2,i_edgeinfo(3,i_edge)) = MAX(r_emax(2,i_edgeinfo(3,i_edge)), r_emean(2,i_edgeinfo(4,i_edge)))
        r_emax(2,i_edgeinfo(4,i_edge)) = MAX(r_emax(2,i_edgeinfo(4,i_edge)), r_emean(2,i_edgeinfo(3,i_edge)))
        r_emax(3,i_edgeinfo(3,i_edge)) = MAX(r_emax(3,i_edgeinfo(3,i_edge)), r_emean(3,i_edgeinfo(4,i_edge)))
        r_emax(3,i_edgeinfo(4,i_edge)) = MAX(r_emax(3,i_edgeinfo(4,i_edge)), r_emean(3,i_edgeinfo(3,i_edge)))
        r_emax(4,i_edgeinfo(3,i_edge)) = MAX(r_emax(4,i_edgeinfo(3,i_edge)), r_emean(4,i_edgeinfo(4,i_edge)))
        r_emax(4,i_edgeinfo(4,i_edge)) = MAX(r_emax(4,i_edgeinfo(4,i_edge)), r_emean(4,i_edgeinfo(3,i_edge)))
        r_emax(5,i_edgeinfo(3,i_edge)) = MAX(r_emax(5,i_edgeinfo(3,i_edge)), r_emean(5,i_edgeinfo(4,i_edge)))
        r_emax(5,i_edgeinfo(4,i_edge)) = MAX(r_emax(5,i_edgeinfo(4,i_edge)), r_emean(5,i_edgeinfo(3,i_edge)))
        r_emax(6,i_edgeinfo(3,i_edge)) = MAX(r_emax(6,i_edgeinfo(3,i_edge)), r_emean(6,i_edgeinfo(4,i_edge)))
        r_emax(6,i_edgeinfo(4,i_edge)) = MAX(r_emax(6,i_edgeinfo(4,i_edge)), r_emean(6,i_edgeinfo(3,i_edge)))
      END IF
    END DO

!--- loop over all elements to do the limiting
    DO i_elmt=1, i_numelmt

!--- limit in surface elevation
      r_hext = MAXVAL(r_Q(1,i_elementdofs(:,i_elmt))+r_S(1,i_elementdofs(:,i_elmt)))
      r_alpha = 1.0_GRID_SR
      IF (r_hext > r_emean(4,i_elmt)+1e-16) &
        r_alpha = MIN(r_alpha, (r_emax(4,i_elmt) - r_emean(4,i_elmt)) / (r_hext - r_emean(4,i_elmt)))
      r_hext = MINVAL(r_Q(1,i_elementdofs(:,i_elmt))+r_S(1,i_elementdofs(:,i_elmt)))
      IF (r_hext < r_emean(4,i_elmt)-1e-16) &
        r_alpha = MIN(r_alpha, (r_emean(4,i_elmt) - r_emin(4,i_elmt)) / (r_emean(4,i_elmt) - r_hext))

      r_corr    = (r_alpha-1.0_GRID_SR) * (r_Q(1,i_elementdofs(:,i_elmt))+r_S(1,i_elementdofs(:,i_elmt)) - r_emean(4,i_elmt))
      r_corr(3) = -(r_corr(1) + r_corr(2))

      r_Qlim(1,i_elementdofs(:,i_elmt)) = r_Q(1,i_elementdofs(:,i_elmt)) + r_corr

!--- positivity preserving limiter ala Bunya et al. (2009)
      IF(MINVAL(r_Qlim(1,i_elementdofs(:,i_elmt))) < 0.0_GRID_SR) THEN
        l_loc = .TRUE.
        i_dofssorted(1) = MINLOC(r_Qlim(1,i_elementdofs(:,i_elmt)), 1, l_loc)
        l_loc(i_dofssorted(1)) = .FALSE.
        i_dofssorted(2) = MINLOC(r_Qlim(1,i_elementdofs(:,i_elmt)), 1, l_loc)
        l_loc(i_dofssorted(2)) = .FALSE.
        i_dofssorted(3) = MINLOC(r_Qlim(1,i_elementdofs(:,i_elmt)), 1, l_loc)

        r_htmp = r_Qlim(1,i_elementdofs(:,i_elmt))
        r_Qlim(1,i_elementdofs(i_dofssorted(1),i_elmt)) = 0.0_GRID_SR
        r_Qlim(1,i_elementdofs(i_dofssorted(2),i_elmt)) = &
          MAX(0.0_GRID_SR, r_htmp(i_dofssorted(2)) + r_htmp(i_dofssorted(1))/2.0_GRID_SR)
        r_Qlim(1,i_elementdofs(i_dofssorted(3),i_elmt)) = &
          MAX(0.0_GRID_SR, r_htmp(i_dofssorted(3)) + r_htmp(i_dofssorted(1)) - &     ! Note: the max-function is a fix to enforce positivity,
          (r_Qlim(1,i_elementdofs(i_dofssorted(2),i_elmt))-r_htmp(i_dofssorted(2)))) ! it might affect mass conservation
      END IF

!--- limit momentum based on velocity
!--- compute limited u-velocity
      r_utmp      = velocity(r_Q(1,i_elementdofs(1,i_elmt)), r_Q(2,i_elementdofs(1,i_elmt)))
      r_ulim(1,:) = MAX(MIN(r_utmp, r_emax(5, i_elmt)), r_emin(5, i_elmt))
      r_utmp      = velocity(r_Q(1,i_elementdofs(2,i_elmt)), r_Q(2,i_elementdofs(2,i_elmt)))
      r_ulim(2,:) = MAX(MIN(r_utmp, r_emax(5, i_elmt)), r_emin(5, i_elmt))
      r_utmp      = velocity(r_Q(1,i_elementdofs(3,i_elmt)), r_Q(2,i_elementdofs(3,i_elmt)))
      r_ulim(3,:) = MAX(MIN(r_utmp, r_emax(5, i_elmt)), r_emin(5, i_elmt))

!--- compute u-velocities at other nodes based on linear momentum assumption
      r_ulim(1,1) = velocity(r_Qlim(1,i_elementdofs(1,i_elmt)), &
                             3.0_GRID_SR*r_emean(2,i_elmt)-r_Qlim(1,i_elementdofs(2,i_elmt))*r_ulim(2,1) &
                                                          -r_Qlim(1,i_elementdofs(3,i_elmt))*r_ulim(3,1))
      r_ulim(2,2) = velocity(r_Qlim(1,i_elementdofs(2,i_elmt)), &
                             3.0_GRID_SR*r_emean(2,i_elmt)-r_Qlim(1,i_elementdofs(1,i_elmt))*r_ulim(1,2) &
                                                          -r_Qlim(1,i_elementdofs(3,i_elmt))*r_ulim(3,2))
      r_ulim(3,3) = velocity(r_Qlim(1,i_elementdofs(3,i_elmt)), &
                             3.0_GRID_SR*r_emean(2,i_elmt)-r_Qlim(1,i_elementdofs(1,i_elmt))*r_ulim(1,3) &
                                                          -r_Qlim(1,i_elementdofs(2,i_elmt))*r_ulim(2,3))

      IF((MAXVAL(r_ulim(:,1))-MINVAL(r_ulim(:,1)) <= MAXVAL(r_ulim(:,2))-MINVAL(r_ulim(:,2))) .AND. &
         (MAXVAL(r_ulim(:,1))-MINVAL(r_ulim(:,1)) <= MAXVAL(r_ulim(:,3))-MINVAL(r_ulim(:,3)))) THEN
        r_Qlim(2,i_elementdofs(:,i_elmt)) = r_Qlim(1,i_elementdofs(:,i_elmt))*r_ulim(:,1)
      ELSE
        IF(MAXVAL(r_ulim(:,2))-MINVAL(r_ulim(:,2)) <= MAXVAL(r_ulim(:,3))-MINVAL(r_ulim(:,3))) THEN
          r_Qlim(2,i_elementdofs(:,i_elmt)) = r_Qlim(1,i_elementdofs(:,i_elmt))*r_ulim(:,2)
        ELSE
          r_Qlim(2,i_elementdofs(:,i_elmt)) = r_Qlim(1,i_elementdofs(:,i_elmt))*r_ulim(:,3)
        END IF
      END IF

!--- compute limited v-velocity
      r_utmp      = velocity(r_Q(1,i_elementdofs(1,i_elmt)), r_Q(3,i_elementdofs(1,i_elmt)))
      r_ulim(1,:) = MAX(MIN(r_utmp, r_emax(6, i_elmt)), r_emin(6, i_elmt))
      r_utmp      = velocity(r_Q(1,i_elementdofs(2,i_elmt)), r_Q(3,i_elementdofs(2,i_elmt)))
      r_ulim(2,:) = MAX(MIN(r_utmp, r_emax(6, i_elmt)), r_emin(6, i_elmt))
      r_utmp      = velocity(r_Q(1,i_elementdofs(3,i_elmt)), r_Q(3,i_elementdofs(3,i_elmt)))
      r_ulim(3,:) = MAX(MIN(r_utmp, r_emax(6, i_elmt)), r_emin(6, i_elmt))

!--- compute v-velocities at other nodes based on linear momentum assumption
      r_ulim(1,1) = velocity(r_Qlim(1,i_elementdofs(1,i_elmt)), &
                             3.0_GRID_SR*r_emean(3,i_elmt)-r_Qlim(1,i_elementdofs(2,i_elmt))*r_ulim(2,1) &
                                                          -r_Qlim(1,i_elementdofs(3,i_elmt))*r_ulim(3,1))
      r_ulim(2,2) = velocity(r_Qlim(1,i_elementdofs(2,i_elmt)), &
                             3.0_GRID_SR*r_emean(3,i_elmt)-r_Qlim(1,i_elementdofs(1,i_elmt))*r_ulim(1,2) &
                                                          -r_Qlim(1,i_elementdofs(3,i_elmt))*r_ulim(3,2))
      r_ulim(3,3) = velocity(r_Qlim(1,i_elementdofs(3,i_elmt)), &
                             3.0_GRID_SR*r_emean(3,i_elmt)-r_Qlim(1,i_elementdofs(1,i_elmt))*r_ulim(1,3) &
                                                          -r_Qlim(1,i_elementdofs(2,i_elmt))*r_ulim(2,3))

      IF((MAXVAL(r_ulim(:,1))-MINVAL(r_ulim(:,1)) <= MAXVAL(r_ulim(:,2))-MINVAL(r_ulim(:,2))) .AND. &
         (MAXVAL(r_ulim(:,1))-MINVAL(r_ulim(:,1)) <= MAXVAL(r_ulim(:,3))-MINVAL(r_ulim(:,3)))) THEN
        r_Qlim(3,i_elementdofs(:,i_elmt)) = r_Qlim(1,i_elementdofs(:,i_elmt))*r_ulim(:,1)
      ELSE
        IF(MAXVAL(r_ulim(:,2))-MINVAL(r_ulim(:,2)) <= MAXVAL(r_ulim(:,3))-MINVAL(r_ulim(:,3))) THEN
          r_Qlim(3,i_elementdofs(:,i_elmt)) = r_Qlim(1,i_elementdofs(:,i_elmt))*r_ulim(:,2)
        ELSE
          r_Qlim(3,i_elementdofs(:,i_elmt)) = r_Qlim(1,i_elementdofs(:,i_elmt))*r_ulim(:,3)
        END IF
      END IF

    END DO

    r_Q(:,:) = r_Qlim(:,:)

!--- deallocate workspace
    DEALLOCATE(i_dofssorted, r_corr, r_htmp, r_emean, r_emin, r_emax, r_ulim, r_Qlim, l_loc)

  END SUBROUTINE limiter

!*******************************************************************************
END MODULE DG_limiter
