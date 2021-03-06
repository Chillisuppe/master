!*******************************************************************************
!
!> @file  IO_netcdfplot_dg_remesh.f90
!> @brief contains module IO_netcdfplot
!
!*******************************************************************************
!
! VERSION(S):
!  1. original version                    j. behrens      01/2000
!  2. compliant to amatos 1.0             j. behrens      12/2000
!  3. compliant to amatos 1.2             j. behrens      03/2002
!  3. compliant to amatos 2.0             j. behrens      07/2003
!  4. rewritten to get more data          s. vater        05/2012
!  5. rewritten to use IO_ncugrid module  s. vater        03/2014
!
!*******************************************************************************
! MODULE DESCRIPTION:
!> @brief creates output in NetCDF file format
!
MODULE IO_netcdfplot
  USE IO_ncugrid
  USE FLASH_parameters
  USE GRID_api
  USE DG_equation, ONLY : p_iovars, i_nprogvars, i_nsrcterms
  USE DG_utils

  IMPLICIT NONE

  LOGICAL :: l_subtriangulation = .TRUE.

  PRIVATE
  PUBLIC :: plot_netcdf

  CONTAINS

!*******************************************************************************
! DESCRIPTION of [SUBROUTINE plot_netcdf]:
!> @brief creates output in NetCDF file format using the UGRID conventions.
!>        In this version the mesh is transformed to a mesh where each element
!>        is on its own and there is no connectivity to its neighbors. This
!>        preserves the discontiuous character of the data.
!>
!> @param[in]         p_ghand         grid handling data structure
!> @param[in]         i_time          time stamp for file naming (opt.)
!> @param[in]         r_time          current model time
!> @param[in]         i_numnode       total number of nodes
!> @param[in]         i_numelmt       total number of elements
!> @param[in]         i_numedge       total number of edges
!> @param[in]         i_elmtnodes     node indices of each element; DIMENSION(3, i_numelmt)
!> @param[in]         i_elmtdofs      DOF indices corresponding to each element; DIMENSION(i_faceunknowns, i_numelmt)
!> @param[in]         i_edgenodes     global node indices for each edge; DIMENSION(2,i_numedge)
!> @param[in]         r_coonod
!> @param[in]         r_coodof
!> @param[in]         r_Q             discrete solution vector
!> @param[in]         r_S
!> @param[in]         r_sol           discrete representation of exact solution vector
!> @param[in]         r_einvvander
!> @param[in]         i_nsubpts_edge
!> @param[in]         l_solplot       switch for output of exact solution
!
  SUBROUTINE plot_netcdf(p_ghand, i_time, r_time, i_numnode, i_numelmt, &
                         i_numedge, i_elmtnodes, i_elmtdofs, i_edgenodes, r_coonod, &
                         r_coodof, r_Q, r_S, r_sol, r_einvvander, i_nsubpts_edge, l_solplot)

    IMPLICIT NONE

    TYPE (grid_handle),                       INTENT(IN)        :: p_ghand
    INTEGER, OPTIONAL,                        INTENT(IN)        :: i_time
    REAL (KIND = GRID_SR),                    INTENT(IN)        :: r_time
    INTEGER,                                  INTENT(IN)        :: i_numnode, i_numelmt, i_numedge, i_nsubpts_edge
    INTEGER (KIND = GRID_SI), DIMENSION(:,:), INTENT(IN)        :: i_elmtnodes, i_elmtdofs, i_edgenodes
    REAL (KIND = GRID_SR),    DIMENSION(:,:), INTENT(IN)        :: r_coonod, r_coodof
    REAL (KIND = GRID_SR),    DIMENSION(:,:), INTENT(IN)        :: r_Q, r_S, r_sol, r_einvvander
    LOGICAL,                                  INTENT(IN)        :: l_solplot

!--- local declarations
    INTEGER (KIND = GRID_SI)                                    :: i_var, i_tcnt, &
      i_alct, i_elmt, i_eunknowns, i_patchelem, i_patchcnt, i_sub, i_shift, i_currelmt, i_cnt1, i_cnt2, i_niovars
    INTEGER (KIND = GRID_SI), SAVE                              :: i_nsubpts, i_nsubelmts
    INTEGER, SAVE                                               :: i_timecount = 0
    CHARACTER (len=32)                                          :: c_file, c_mesh
    CHARACTER (len=64)                                          :: c_title, c_tmp
    INTEGER, DIMENSION(GRID_elementnodes)                       :: i_ndindices
    TYPE (ncugrid_vardatatype), DIMENSION(:), ALLOCATABLE       :: p_vdataarr
    INTEGER (KIND = GRID_SI), DIMENSION(:,:), ALLOCATABLE       :: i_elmtnodesdg
    INTEGER (KIND = GRID_SI), DIMENSION(:,:), ALLOCATABLE, SAVE :: i_subtriang
    INTEGER (KIND = GRID_SI), DIMENSION(:), ALLOCATABLE, TARGET :: i_eltidx, i_nodeidx, i_eltlev, i_crseltlev
    REAL (KIND = GRID_SR), DIMENSION(:,:), ALLOCATABLE          :: r_coonoddg, r_tmpdat, r_bary, r_V
    REAL (KIND = GRID_SR), DIMENSION(:,:), ALLOCATABLE, SAVE    :: r_psi
    REAL (KIND = GRID_SR), DIMENSION(:,:), ALLOCATABLE, TARGET  :: r_vdata

!--- check input for optional parameter
    IF(PRESENT(i_time)) THEN
      i_tcnt = i_time
    ELSE
      i_tcnt = i_timecount
    END IF
    i_timecount = i_timecount + 1

    c_mesh = 'Mesh2'

    IF (l_solplot) THEN
      i_niovars = 2 * i_nprogvars + i_nsrcterms
    ELSE
      i_niovars = i_nprogvars + i_nsrcterms
    END IF

    ALLOCATE(p_vdataarr(i_niovars+3), stat=i_alct)
    IF (i_alct /= 0) CALL grid_error(c_error='[plot_netcdf]: could not allocate data arrays')

!--- create the title
    WRITE(c_title,*) 'netCDF output from ',TRIM(GRID_parameters%program_name)
    c_title = ADJUSTL(c_title)

!--- create generic file name
    WRITE(c_file, '(A, I8.8, A3)') TRIM(GRID_parameters%program_name), i_tcnt, '.nc'

!--- define data variables
    DO i_var = 1, i_niovars
      WRITE(c_tmp, '(3A)') TRIM(c_mesh), '_', TRIM(p_iovars(i_var)%c_varname)
      p_vdataarr(i_var)%c_varname       = TRIM(c_tmp)
      p_vdataarr(i_var)%c_long_name     = p_iovars(i_var)%c_long_name
      p_vdataarr(i_var)%c_standard_name = p_iovars(i_var)%c_standard_name
      p_vdataarr(i_var)%c_units         = p_iovars(i_var)%c_units
      p_vdataarr(i_var)%c_location      = 'node'
      p_vdataarr(i_var)%i_datatype      = 1
    END DO

    WRITE(c_tmp, '(A, A10)') TRIM(c_mesh), '_elmtindex'
    p_vdataarr(i_niovars+1)%c_varname       = TRIM(c_tmp)
    p_vdataarr(i_niovars+1)%c_long_name     = 'element index'
    p_vdataarr(i_niovars+1)%c_standard_name = 'element_index' ! not CF conform!
    p_vdataarr(i_niovars+1)%c_units         = 'None'
    p_vdataarr(i_niovars+1)%c_location      = 'face'
    p_vdataarr(i_niovars+1)%i_datatype      = 0

    WRITE(c_tmp, '(A, A10)') TRIM(c_mesh), '_nodeindex'
    p_vdataarr(i_niovars+2)%c_varname       = TRIM(c_tmp)
    p_vdataarr(i_niovars+2)%c_long_name     = 'node index'
    p_vdataarr(i_niovars+2)%c_standard_name = 'node_index' ! not CF conform!
    p_vdataarr(i_niovars+2)%c_units         = 'None'
    p_vdataarr(i_niovars+2)%c_location      = 'node'
    p_vdataarr(i_niovars+2)%i_datatype      = 0

    WRITE(c_tmp, '(A, A6)') TRIM(c_mesh), '_level'
    p_vdataarr(i_niovars+3)%c_varname       = TRIM(c_tmp)
    p_vdataarr(i_niovars+3)%c_long_name     = 'grid level'
    p_vdataarr(i_niovars+3)%c_standard_name = 'grid_level' ! not CF conform!
    p_vdataarr(i_niovars+3)%c_units         = 'None'
    p_vdataarr(i_niovars+3)%c_location      = 'face'
    p_vdataarr(i_niovars+3)%i_datatype      = 0

!--- allocate workspace
    i_eunknowns = GRID_femtypes%p_type(FEM_DG)%sig%i_unknowns

!--- calculate interpolation matrix for subtriangulation (first routine call only)
    IF (l_subtriangulation) THEN

      IF (i_nsubpts_edge == 1) CALL grid_error(c_error='[plot_netcdf]: SUBTRIANG_PTS has to be greater than 1')
      i_nsubpts = (i_nsubpts_edge+1) * i_nsubpts_edge / 2
      i_nsubelmts = (i_nsubpts_edge-1)**2

      ALLOCATE(i_subtriang(i_nsubelmts, GRID_elementnodes), &
               r_bary(GRID_elementnodes, i_nsubpts), &
               r_V(i_nsubpts, i_eunknowns), &
               r_psi(i_nsubpts, i_eunknowns), stat=i_alct)
      IF (i_alct /= 0) CALL grid_error(c_error='[plot_netcdf]: could not allocate data arrays')

      i_sub = 0
      DO i_cnt1 = 0, i_nsubpts_edge-1
        DO i_cnt2 = 0, i_nsubpts_edge-1-i_cnt1
          i_sub = i_sub + 1
          r_bary(:,i_sub) = 1.0_GRID_SR / REAL(i_nsubpts_edge-1, GRID_SR) * &
            [REAL(i_cnt1, GRID_SR), REAL(i_cnt2, GRID_SR), REAL(i_nsubpts_edge-1-i_cnt1-i_cnt2, GRID_SR)]
        END DO !i_cnt1
      END DO !i_cnt2

!--- determine subtriangulation indexing
      i_sub = 1
      DO i_cnt1 = 1, i_nsubpts_edge-1
        i_subtriang(i_sub,:) = [i_cnt1+i_nsubpts_edge, i_cnt1+1, i_cnt1]
        i_sub   = i_sub + 1
        i_shift = 0
        DO i_cnt2 = 1, i_cnt1-1
          i_subtriang(i_sub  ,:) = [i_cnt1+i_nsubpts_edge-i_cnt2+1+i_shift, &
                                    i_cnt1+i_shift, &
                                    i_cnt1+i_nsubpts_edge-i_cnt2+i_shift]
          i_subtriang(i_sub+1,:) = [i_cnt1+2*i_nsubpts_edge-2*i_cnt2+i_shift, &
                                    i_cnt1+i_nsubpts_edge-i_cnt2+1+i_shift, &
                                    i_cnt1+i_nsubpts_edge-i_cnt2+i_shift]
          i_sub   = i_sub + 2
          i_shift = i_shift + i_nsubpts_edge - i_cnt2
        END DO !i_cnt1
      END DO !i_cnt2

      CALL Vandermonde2D(GRID_femtypes%p_type(FEM_DG)%sig%i_degree, r_bary, r_V)
      r_psi = MATMUL(r_V, r_einvvander)

      DO i_cnt1 = 1, i_nsubpts
        DO i_cnt2 = 1, i_eunknowns
          IF (ABS(r_psi(i_cnt1,i_cnt2)) < GRID_EPS) r_psi(i_cnt1,i_cnt2) = 0.0_GRID_SR
        END DO !i_cnt2
      END DO !i_cnt1

      DEALLOCATE(r_bary, r_V)
      l_subtriangulation = .FALSE.
    END IF

    ALLOCATE(i_eltidx(i_numelmt*i_nsubelmts), &
             i_nodeidx(i_numelmt*GRID_elementnodes*i_nsubelmts), &
             i_eltlev(i_numelmt*i_nsubelmts), &
             i_crseltlev(i_numelmt), &
             r_tmpdat(i_niovars, i_numelmt*i_eunknowns), &
             r_vdata(i_niovars, i_numelmt*GRID_elementnodes*i_nsubelmts), &
             i_elmtnodesdg(GRID_elementnodes, i_numelmt*i_nsubelmts), &
             r_coonoddg(GRID_dimension, i_numelmt*GRID_elementnodes*i_nsubelmts), stat=i_alct)
    IF (i_alct /= 0) CALL grid_error(c_error='[plot_netcdf]: could not allocate data arrays')

!--- get arrays with some mesh data
    CALL grid_getinfo(p_ghand, l_finelevel=.TRUE., l_relative=.TRUE., &
                      i_elementlevel=i_crseltlev)

!--- construct temporary data array
    r_tmpdat(1:i_nprogvars,:) = r_Q
    r_tmpdat(i_nprogvars+1:i_nprogvars+i_nsrcterms,:) = r_S
    IF (l_solplot) r_tmpdat(i_nprogvars+i_nsrcterms+1:i_niovars,:) = r_sol
    DO i_var = 1, i_niovars
      r_tmpdat(i_var,:) = r_tmpdat(i_var,:) * p_iovars(i_var)%r_factor
    END DO

!--- compute remeshed data
    IF (i_eunknowns > 1_GRID_SI) THEN
      DO i_elmt = 1, i_numelmt
        DO i_sub = 1, i_nsubelmts
          i_currelmt                  = (i_elmt-1)*i_nsubelmts + i_sub
          i_eltidx(i_currelmt)        = i_elmt
          i_ndindices                 = (i_currelmt-1)*3 + [1, 2, 3]
          i_elmtnodesdg(:,i_currelmt) = i_ndindices
          i_nodeidx(i_ndindices)      = i_elmtnodes(:,i_elmt)
          r_coonoddg(:,i_ndindices)   = MATMUL(r_coodof(:,i_elmtdofs(:,i_elmt)), &
                                        TRANSPOSE(r_psi(i_subtriang(i_sub,:),:)))
          r_vdata(:,i_ndindices)      = MATMUL(r_tmpdat(:,i_elmtdofs(:,i_elmt)), &
                                        TRANSPOSE(r_psi(i_subtriang(i_sub,:),:)))
          i_eltlev(i_currelmt)        = i_crseltlev(i_elmt)
        END DO
      END DO
    ELSE
      DO i_elmt = 1, i_numelmt
        DO i_sub = 1, i_nsubelmts
          i_currelmt                  = (i_elmt-1)*i_nsubelmts + i_sub
          i_eltidx(i_currelmt)        = i_elmt
          i_ndindices                 = (i_currelmt-1)*3 + [1, 2, 3]
          i_elmtnodesdg(:,i_currelmt) = i_ndindices
          i_nodeidx(i_ndindices)      = i_elmtnodes(:,i_elmt)
          r_coonoddg(:,i_ndindices)   = r_coonod(:,i_elmtnodes(:,i_elmt))
          r_vdata(:,i_ndindices)      = MATMUL(r_tmpdat(:,i_elmtdofs(:,i_elmt)), &
                                        TRANSPOSE(r_psi(i_subtriang(i_sub,:),:)))
          i_eltlev(i_currelmt)        = i_crseltlev(i_elmt)
        END DO
      END DO
    END IF
      
!--- fill structure with data
    DO i_var = 1, i_niovars
      p_vdataarr(i_var)%p_rvardata => r_vdata(i_var,:)
    END DO
    p_vdataarr(i_niovars+1)%p_ivardata => i_eltidx
    p_vdataarr(i_niovars+2)%p_ivardata => i_nodeidx
    p_vdataarr(i_niovars+3)%p_ivardata => i_eltlev

!--- create the file for timestep data
    CALL ncugrid_createmesh(c_file, c_mesh, i_numelmt*GRID_elementnodes*i_nsubelmts, i_numelmt*i_nsubelmts, &
                            i_elmtnodesdg, r_coonoddg, c_title, TRIM(GRID_parameters%author_affil2))

!--- write data to file
    CALL ncugrid_putvariablearray(c_file, c_mesh, p_vdataarr, r_time)

!--- some examples for ncugrid_putvariable:
!     WRITE(c_tmp, '(A, A9)') TRIM(c_mesh), 'depthXX'
!     CALL ncugrid_putvariable(c_file, c_tmp, 'fluid depth', 'fluid_depth', 'm', &
!                              c_mesh, 'node', r_vdata(1,:), r_time)
!
!     WRITE(c_tmp, '(A, A8)') TRIM(c_mesh), 'levelXX'
!     CALL ncugrid_putvariable(c_file, c_tmp, 'grid level', 'grid_level', 'none', &
!                              c_mesh, 'face', i_eltlev, r_time)

!--- deallocate data arrays
    DO i_var=1,SIZE(p_vdataarr)
      NULLIFY(p_vdataarr(i_var)%p_ivardata)
      NULLIFY(p_vdataarr(i_var)%p_rvardata)
    END DO
    DEALLOCATE(p_vdataarr, r_tmpdat, r_vdata, i_eltlev, i_crseltlev, &
               i_elmtnodesdg, r_coonoddg, i_eltidx, i_nodeidx)

  END SUBROUTINE plot_netcdf

!*******************************************************************************
END MODULE IO_netcdfplot
