!*******************************************************************************
!
!> @file  DG_initial_coupled.f90
!> @brief contains module DG_initial
!
!*******************************************************************************
! MODULE DESCRIPTION:
!> @brief initialize fields for coupled test case (dynamic displacement)
!
MODULE DG_initial

  USE netcdf
  USE GRID_api
  USE FLASH_parameters
  USE DG_equation, ONLY : i_nprogvars, i_nsrcterms, i_valQ, i_valS
  USE IO_equation, ONLY : p_equationparam
  USE IO_ncugrid

  PRIVATE
  PUBLIC  i_ntstintparam, i_ntstrealparam, i_ntstcharparam, i_ntstlogparam, &
          initial_setparam, p_testparam, initialize_testcase, initialvalues_iter, &
          source_update, dg_elmt_solution

  TYPE (test_param)                                   :: p_testparam
  INTEGER, PARAMETER                                  :: i_ntstintparam  = 1
  INTEGER, PARAMETER                                  :: i_ntstrealparam = 2
  INTEGER, PARAMETER                                  :: i_ntstcharparam = 1
  INTEGER, PARAMETER                                  :: i_ntstlogparam  = 0

!--- variables for uplift information
    INTEGER (KIND = GRID_SI)                                    :: i_xdim, i_ydim, i_tdim
    REAL (KIND = GRID_SR), DIMENSION(:), ALLOCATABLE, SAVE      :: r_xvec, r_yvec, r_timevec
    REAL (KIND = GRID_SR), DIMENSION(:,:,:), ALLOCATABLE, SAVE  :: r_uplift


  CONTAINS

!*******************************************************************************
! DESCRIPTION of [SUBROUTINE initial_setparam]:
!> @brief initializes test case specific parameter list
!
  SUBROUTINE initial_setparam

    IMPLICIT NONE

    ! file name for uplift information
    p_testparam%tchar(1)%c_keyword = 'UPLIFT_FILE'

    ! type of source (1 = dynamic, 2 = uplift at end, 3 = max uplift)
    p_testparam%tint(1)%c_keyword  = 'SOURCE_TYPE'
    p_testparam%tint(1)%i_size     = 1

    ! beach slope
    p_testparam%treal(1)%c_keyword = 'BEACH_SLOPE'
    p_testparam%treal(1)%i_size    = 1

    ! position of beach toe
    p_testparam%treal(2)%c_keyword = 'POS_BEACHTOE'
    p_testparam%treal(2)%i_size    = 1

!     ! parameter description
!     p_testparam%tchar(1)%c_keyword = 'PARAM_CHAR'
!
!     ! parameter description
!     p_testparam%tlog(1)%c_keyword  = 'PARAM_LOG'

  END SUBROUTINE initial_setparam

!*******************************************************************************
! DESCRIPTION of [SUBROUTINE initialize_testcase]:
!> @brief initializes testcase (e.g. setup constants, auxiliary fields)
!>
!> @param[in]       p_ghand       grid handling data structure
!> @param[in,out]   p_param       control structure for global parameters
!> @param[in]       r_time        current model time
!>
!> @note This is just a dummy in this case.
!
  SUBROUTINE initialize_testcase(p_ghand, p_param, r_time)

    IMPLICIT NONE

    TYPE (grid_handle), INTENT(in)                        :: p_ghand
    TYPE (control_struct), INTENT(inout)                  :: p_param
    REAL (KIND = GRID_SR), INTENT(in)                     :: r_time

    CALL read_disp(p_testparam%tchar(1)%p_value)

  END SUBROUTINE initialize_testcase

!*******************************************************************************
! DESCRIPTION of [SUBROUTINE initialvalues_iter]:
!> @brief initializes a grid with values
!>
!> @param[in]       p_ghand       grid handling data structure
!> @param[in,out]   p_param       control structure for global parameters
!> @param[in]       r_time        current model time
!>
!> @note This routine initializes the flow field.
!
  SUBROUTINE initialvalues_iter(p_ghand, p_param, r_time)

    IMPLICIT NONE

    TYPE (grid_handle), INTENT(in)                        :: p_ghand
    TYPE (control_struct), INTENT(inout)                  :: p_param
    REAL (KIND = GRID_SR), INTENT(in)                     :: r_time

!--- local declarations
    INTEGER (KIND = GRID_SI)                              :: i_faceunknowns, i_numelmt, &
      i_num, i_alct
    INTEGER (KIND = GRID_SI), DIMENSION(:,:), ALLOCATABLE :: i_elmtdofs
    REAL (KIND = GRID_SR), DIMENSION(:,:), ALLOCATABLE    :: r_coodof, r_Q, r_S

!--- initialize some constants
    i_faceunknowns = GRID_femtypes%p_type(FEM_DG)%sig%i_unknowns
    i_numelmt      = p_ghand%i_enumfine
    i_num          = i_numelmt*i_faceunknowns

!--- allocate workspace
    ALLOCATE(i_elmtdofs(i_faceunknowns, i_numelmt), &
             r_coodof(GRID_dimension, i_num), &
             r_Q(i_nprogvars,i_num), r_S(i_nsrcterms, i_num), STAT=i_alct)
    IF(i_alct /= 0) CALL grid_error(c_error='[initialvalues_iter]: could not allocate grid field')

!--- get grid information
    CALL grid_getinfo(p_ghand, i_femtype=FEM_DG, l_relative=.TRUE., &
                      l_finelevel=.TRUE., i_elementdofs=i_elmtdofs, r_dofcoordinates=r_coodof)

!--- initialize fields
    CALL bathy_init(i_numelmt, r_coodof, r_S)
    CALL field_init(p_param, i_numelmt, i_elmtdofs, r_coodof, r_Q, r_S)
    CALL source_update(p_param, r_time, i_numelmt, i_elmtdofs, r_coodof, .TRUE., &
                       r_Q, r_S)

!--- update grid information
    CALL grid_putinfo(p_ghand, l_finelevel=.TRUE., i_arraypoint=i_valQ, r_dofvalues=r_Q)
    CALL grid_putinfo(p_ghand, l_finelevel=.TRUE., i_arraypoint=i_valS, r_dofvalues=r_S)

!--- deallocate workspace
    DEALLOCATE(i_elmtdofs, r_coodof, r_Q, r_S)

  END SUBROUTINE initialvalues_iter

!*******************************************************************************
! DESCRIPTION of [SUBROUTINE source_update]:
!> @brief updates the source terms
!>
!> @param[in]       p_param         control structure for global parameters
!> @param[in]       r_time          model time at which source term is evaluated
!> @param[in]       i_numelmt       number of elements
!> @param[in]       i_elmtdofs      global element-dof relation matrix
!> @param[in]       r_coodof        coordinates of DOFs
!> @param[in]       l_gridchanged
!> @param[in,out]   r_Q             discrete solution vector
!> @param[in,out]   r_S             discrete fields used in source terms
!
  SUBROUTINE source_update(p_param, r_time, i_numelmt, i_elmtdofs, r_coodof, &
                           l_gridchanged, r_Q, r_S)

    IMPLICIT NONE

    TYPE (control_struct),                    INTENT(in)    :: p_param
    REAL (KIND = GRID_SR),                    INTENT(in)    :: r_time
    INTEGER (KIND = GRID_SI),                 INTENT(in)    :: i_numelmt
    INTEGER (KIND = GRID_SI), DIMENSION(:,:), INTENT(in)    :: i_elmtdofs
    REAL (KIND = GRID_SR), DIMENSION(:,:),    INTENT(in)    :: r_coodof
    LOGICAL,                                  INTENT(in)    :: l_gridchanged
    REAL (KIND = GRID_SR), DIMENSION(:,:),    INTENT(inout) :: r_Q, r_S

    CALL bathy_init(i_numelmt, r_coodof, r_S)
    CALL bathy_lift(r_time, i_numelmt, r_coodof, r_S)

  END SUBROUTINE source_update

!*******************************************************************************
! DESCRIPTION of [SUBROUTINE dg_elmt_solution]:
!> @brief computes the exact solution approximated in the DG space to compare
!>        with in diagnostics
!>
!> @param[in]       r_coodof        coordinates where solution should be computed
!> @param[in]       r_time          model time
!> @param[in]       i_numelmt       number of elements
!> @param[in]       i_elmtdofs      DOFs of element
!> @param[in]       r_coonod        coordinates of element nodes
!> @param[in]       i_elmtnodes     nodes of element
!> @param           r_sol           computed exact solution
!
  SUBROUTINE dg_elmt_solution(r_coodof, r_time, i_numelmt, i_elmtdofs, r_coonod, i_elmtnodes, r_sol)

    IMPLICIT NONE

    REAL (KIND = GRID_SR),    DIMENSION(:,:), INTENT(IN)        :: r_coodof, r_coonod
    REAL (KIND = GRID_SR),                    INTENT(IN)        :: r_time
    INTEGER (KIND = GRID_SI),                 INTENT(IN)        :: i_numelmt
    INTEGER (KIND = GRID_SI), DIMENSION(:,:), INTENT(IN)        :: i_elmtdofs, i_elmtnodes
    REAL (KIND = GRID_SR),    DIMENSION(:,:), INTENT(INOUT)     :: r_sol

    r_sol = 0.0_GRID_SR

  END SUBROUTINE dg_elmt_solution

!*******************************************************************************
! Only internal subroutines follow
!*******************************************************************************
! DESCRIPTION of [SUBROUTINE bathy_init]:
!> @brief initializes the bathymetry field
!>
!> @param           i_numelmt
!> @param           r_coodof      coordinates of DOFs
!
  SUBROUTINE bathy_init(i_numelmt, r_coodof, r_S)

    IMPLICIT NONE

    INTEGER (KIND = GRID_SI),                 INTENT(in)    :: i_numelmt
    REAL (KIND = GRID_SR), DIMENSION(:,:),    INTENT(in)    :: r_coodof
    REAL (KIND = GRID_SR), DIMENSION(:,:),    INTENT(out)   :: r_S

!--- local declarations
    INTEGER (KIND = GRID_SI)                                :: i_count, i_num
    REAL (KIND = GRID_SR)                                   :: r_slope, r_x0

!--- initialize some constants
    i_num   = i_numelmt*GRID_femtypes%p_type(FEM_DG)%sig%i_unknowns
    r_slope = p_testparam%treal(1)%p_value(1)
    r_x0    = p_testparam%treal(2)%p_value(1)

!--- loop over the DOFs
    DO i_count=1,i_num
!       r_S(1,i_count) = MAX(0.0_GRID_SR, &
!         r_slope*(r_coodof(1,i_count)-COS(20.0/90.0*GRID_PI/2.0)*150000.0_GRID_SR) + &
!         p_equationparam%r_depth) * GRID_GRAV
      r_S(1,i_count) = MAX(0.0_GRID_SR, r_slope*(r_coodof(1,i_count)-r_x0)) * GRID_GRAV
    END DO

  END SUBROUTINE bathy_init

!*******************************************************************************
! DESCRIPTION of [SUBROUTINE bathy_lift:
!> @brief modifies the bathymetry field with uplift
!>
!> @param           i_numelmt
!> @param           r_coodof      coordinates of DOFs
!
  SUBROUTINE bathy_lift(r_time, i_numelmt, r_coodof, r_S)

    IMPLICIT NONE

    REAL (KIND = GRID_SR), INTENT(in)                     :: r_time
    INTEGER (KIND = GRID_SI),                 INTENT(in)  :: i_numelmt
    REAL (KIND = GRID_SR), DIMENSION(:,:),    INTENT(in)  :: r_coodof
    REAL (KIND = GRID_SR), DIMENSION(:,:), INTENT(inout)  :: r_S

!--- local declarations
    INTEGER                                               :: i_count, i_num, i_alct, &
                                                             i_dim, i_x, i_y, i_t, &
                                                             i_faceunknowns, i_numnode
    REAL (KIND = GRID_SR)                                 :: r_at, r_ax, r_ay, r_xx1, r_xx2, r_val, r_maxval, r_maxstep
!     INTEGER (KIND = GRID_SI), DIMENSION(:,:), ALLOCATABLE :: i_elmtnodes, i_elmtdofs
!     REAL (KIND = GRID_SR), DIMENSION(:), ALLOCATABLE      :: r_displ
!     REAL (KIND = GRID_SR), DIMENSION(:,:), ALLOCATABLE    :: r_coonod

 !--- allocate workspace
    i_faceunknowns = GRID_femtypes%p_type(FEM_DG)%sig%i_unknowns
!     i_numnode      = p_ghand%i_nnumber
    i_num          = i_numelmt*i_faceunknowns
!     ALLOCATE(i_elmtnodes(3,i_numelmt), i_elmtdofs(i_faceunknowns,i_numelmt), &
!              r_displ(i_numelmt), r_coonod(GRID_dimension,i_numnode), stat=i_alct)
    IF(i_alct /= 0) &
      CALL grid_error(c_error='[bathy_lift]: could not allocate workspace')

!--- determine time snapshot to be used
    SELECT CASE (p_testparam%tint(1)%p_value(1))
    CASE(1) ! dynamic displacement
      DO i_dim=1,i_tdim-1
        IF (r_timevec(i_dim+1)>=r_time) THEN
          i_t = i_dim
          EXIT
        END IF
      END DO
      IF (r_time > r_timevec(i_tdim)) THEN
        i_t  = i_tdim-1
        r_at = 1.0_GRID_SR
      ELSE
        r_at = (r_time - r_timevec(i_t)) / (r_timevec(i_t+1) - r_timevec(i_t))
      END IF

    CASE(2) ! final static displacement
      i_t  = i_tdim-1
      r_at = 1.0_GRID_SR

    CASE(3) ! maximum displacement
      r_maxval = MAXVAL(ABS(r_uplift(:,:,1)))
      i_t      = 1
      r_at     = 0.0_GRID_SR
      DO i_dim=1,i_tdim-1
        r_maxstep = MAXVAL(ABS(r_uplift(:,:,i_dim+1)))
        IF (r_maxstep > r_maxval) THEN
          i_t      = i_dim
          r_at     = 1.0_GRID_SR
          r_maxval = r_maxstep
        END IF
      END DO
      WRITE(*,*) "timestep with max value: ", i_t+1

    CASE DEFAULT
      CALL grid_error(c_error='[bathy_lift]: uplift type not supported!')
    END SELECT

!--- get dof coordinates
!     CALL grid_getinfo(p_ghand, i_femtype=FEM_DG, i_elementnodes=i_elmtnodes, &
!                       i_elementdofs=i_elmtdofs, r_nodecoordinates=r_coonod)

!--- loop over the DOFs
    dof_loop: DO i_count=1,i_num

!--- determine grid box surrounding coordinate value
      DO i_dim=1,i_xdim-1
        IF (r_xvec(i_dim+1)>=r_coodof(1,i_count)) THEN
          i_x = i_dim
          EXIT
        END IF
      END DO

      DO i_dim=1,i_ydim-1
        IF (r_yvec(i_dim+1)>=r_coodof(2,i_count)) THEN
          i_y = i_dim
          EXIT
        END IF
      END DO

!--- compute weights for linear interpolation
      r_ax = (r_coodof(1,i_count) - r_xvec(i_x)) / (r_xvec(i_x+1) - r_xvec(i_x))
      r_ay = (r_coodof(2,i_count) - r_yvec(i_y)) / (r_yvec(i_y+1) - r_yvec(i_y))

!--- interpolate linearly first along x direction, then along y direction
      r_xx1 = (1._GRID_SR - r_at) * ((1._GRID_SR - r_ay) * r_uplift(i_x  ,i_y,i_t  ) + r_ay * r_uplift(i_x  ,i_y+1,i_t  )) + &
               r_at               * ((1._GRID_SR - r_ay) * r_uplift(i_x  ,i_y,i_t+1) + r_ay * r_uplift(i_x  ,i_y+1,i_t+1))
      r_xx2 = (1._GRID_SR - r_at) * ((1._GRID_SR - r_ay) * r_uplift(i_x+1,i_y,i_t  ) + r_ay * r_uplift(i_x+1,i_y+1,i_t  )) + &
               r_at               * ((1._GRID_SR - r_ay) * r_uplift(i_x+1,i_y,i_t+1) + r_ay * r_uplift(i_x+1,i_y+1,i_t+1))

      r_val = (1._GRID_SR - r_ax) * r_xx1 + r_ax * r_xx2

      r_S(1,i_count) = r_S(1,i_count) + r_val * GRID_GRAV
    END DO dof_loop

! !--- save displacement in auxiliary data file
!   CALL ncugrid_createmesh('displ.nc', 'mesh2', i_numnode, i_numelmt, &
!                           i_elmtnodes, r_coonod, 'mesh', 'Universitaet Hamburg')
!
!   CALL ncugrid_putvariable('displ.nc', 'displacement', 'displacement', 'displacement', 'None', &
!                            'mesh2', 'face', r_displ, 0.0_GRID_SR)

!--- deallocate workspace
!     DEALLOCATE(r_displ, i_elmtnodes, i_elmtdofs, r_coonod)

  END SUBROUTINE bathy_lift

!*******************************************************************************
! DESCRIPTION of [SUBROUTINE field_init]:
!> @brief initializes the height and momentum fields
!>
!> @param[in,out]   p_param       control structure for global parameters
!> @param           i_elmtdofs
!> @param           r_coodof      coordinates of DOFs
!
  SUBROUTINE field_init(p_param, i_numelmt, i_elmtdofs, r_coodof, r_Q, r_S)

    IMPLICIT NONE

    TYPE (control_struct),                    INTENT(inout) :: p_param
    INTEGER (KIND = GRID_SI),                 INTENT(in)    :: i_numelmt
    INTEGER (KIND = GRID_SI), DIMENSION(:,:), INTENT(in)    :: i_elmtdofs
    REAL (KIND = GRID_SR), DIMENSION(:,:),    INTENT(in)    :: r_coodof
    REAL (KIND = GRID_SR), DIMENSION(:,:),    INTENT(inout) :: r_Q
    REAL (KIND = GRID_SR), DIMENSION(:,:),    INTENT(in)    :: r_S

!--- local declarations
    INTEGER (KIND = GRID_SI)                                :: i_count, i_num

!--- initialize some constants
    i_num = i_numelmt*GRID_femtypes%p_type(FEM_DG)%sig%i_unknowns

!--- loop over the DOFs
    DO i_count=1,i_num
      r_Q(1,i_count) = MAX(0.0_GRID_SR, p_equationparam%r_depth * GRID_GRAV - r_S(1,i_count))
      r_Q(2,i_count) = 0.0_GRID_SR
      r_Q(3,i_count) = 0.0_GRID_SR
    END DO

  END SUBROUTINE field_init

!*******************************************************************************
! DESCRIPTION of [SUBROUTINE read_disp]:
!> @brief reads dispacement a NetCDF file in COARDS format
!>
!
  SUBROUTINE read_disp(c_upliftfile)

    IMPLICIT NONE

!    INCLUDE "netcdf.inc"

    CHARACTER (len=*), INTENT(in)                           :: c_upliftfile

!---------- local declarations
    INTEGER (KIND = GRID_SI)                                :: i_alct, i_cnt, j_cnt, k_cnt
    INTEGER (KIND = GRID_SI)                                :: i_fileid
    INTEGER (KIND = GRID_SI)                                :: i_dimid, i_varid
    INTEGER (KIND = GRID_SI)                                :: i_ncstat         ! NetCDF library return status
    INTEGER (KIND = GRID_SI), DIMENSION(3)                  :: i_tmp

!---------- open bathymetry file
    i_ncstat = nf90_open(c_upliftfile, NF90_NOWRITE, i_fileid)
    IF(i_ncstat /= NF90_NOERR) &
      CALL grid_error(c_error='[read_disp]: could not open uplift data file')

!---------- determine array sizes.

    i_ncstat = nf90_inq_dimid(i_fileid, 'x', i_dimid)
    IF(i_ncstat /= NF90_NOERR) &
      CALL grid_error(c_error='[read_disp]: could not identify x dimension')
    i_ncstat = nf90_inquire_dimension(i_fileid, i_dimid, len=i_xdim)
    IF(i_ncstat /= NF90_NOERR) &
      CALL grid_error(c_error='[read_disp]: could not read x dimension')
    i_ncstat = nf90_inq_dimid(i_fileid, 'y', i_dimid)
    IF(i_ncstat /= NF90_NOERR) &
      CALL grid_error(c_error='[read_disp]: could not identify y dimension')
    i_ncstat = nf90_inquire_dimension(i_fileid, i_dimid, len=i_ydim)
    IF(i_ncstat /= NF90_NOERR) &
      CALL grid_error(c_error='[read_disp]: could not read y dimension')
    i_ncstat = nf90_inq_dimid(i_fileid, 'time', i_dimid)
    IF(i_ncstat /= NF90_NOERR) &
      CALL grid_error(c_error='[read_disp]: could not identify time dimension')
    i_ncstat = nf90_inquire_dimension(i_fileid, i_dimid, len=i_tdim)
    IF(i_ncstat /= NF90_NOERR) &
      CALL grid_error(c_error='[read_disp]: could not read time dimension')

!---------- allocate data array

    ALLOCATE(r_xvec(i_xdim), r_yvec(i_ydim), r_timevec(i_tdim), &
             r_uplift(i_xdim,i_ydim,i_tdim), stat=i_alct)
    IF(i_alct /= 0) &
      CALL grid_error(c_error='[read_etopo]: could not allocate raw data field')

!---------- read data arrays

    i_ncstat = nf90_inq_varid(i_fileid, 'x', i_varid)
    IF(i_ncstat /= NF90_NOERR) &
      CALL grid_error(c_error='[read_disp]: could not determine x varid')
    i_ncstat = nf90_get_var(i_fileid, i_varid, r_xvec)
    IF(i_ncstat /= NF90_NOERR) &
      CALL grid_error(c_error='[read_disp]: could not read x vector')

    i_ncstat = nf90_inq_varid(i_fileid, 'y', i_varid)
    IF(i_ncstat /= NF90_NOERR) &
      CALL grid_error(c_error='[read_disp]: could not determine y varid')
    i_ncstat = nf90_get_var(i_fileid, i_varid, r_yvec)
    IF(i_ncstat /= NF90_NOERR) &
      CALL grid_error(c_error='[read_disp]: could not read y vector')

    i_ncstat = nf90_inq_varid(i_fileid, 'time', i_varid)
    IF(i_ncstat /= NF90_NOERR) &
      CALL grid_error(c_error='[read_disp]: could not determine time varid')
    i_ncstat = nf90_get_var(i_fileid, i_varid, r_timevec)
    IF(i_ncstat /= NF90_NOERR) &
      CALL grid_error(c_error='[read_disp]: could not read time vector')

    i_ncstat = nf90_inq_varid(i_fileid, 'dz', i_varid)
    IF(i_ncstat /= NF90_NOERR) &
      CALL grid_error(c_error='[read_disp]: could not determine uplift varid')
    i_ncstat = nf90_get_var(i_fileid, i_varid, r_uplift)
    IF(i_ncstat /= NF90_NOERR) &
      CALL grid_error(c_error='[read_disp]: could not read uplift array')

!---------- close topography file

    i_ncstat = nf90_close(i_fileid)
    IF(i_ncstat /= NF90_NOERR) &
      CALL grid_error(c_error='[read_disp]: could not close bathy data file')
    IF(GRID_parameters%iolog > 0) &
      write(GRID_parameters%iolog,*) 'INFO: Closed uplift data file'

  END SUBROUTINE read_disp

!*******************************************************************************
END MODULE DG_initial
