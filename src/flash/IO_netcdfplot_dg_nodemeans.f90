!*******************************************************************************
!
!> @file  IO_netcdfplot_dg_nodemeans.f90
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

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: plot_netcdf

  CONTAINS

!*******************************************************************************
! DESCRIPTION of [SUBROUTINE plot_netcdf]:
!> @brief creates output in NetCDF file format using the UGRID conventions.
!>        In this version the nodal mean values are computed and written
!>        as (continuous) nodal data.
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
!> @param[in]         r_Q
!> @param[in]         r_S
!> @param[in]         r_sol
!> @param[in]         r_einvvander
!> @param[in]         i_nsubpts_edge
!> @param[in]         l_solplot
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
      i_alct, i_node, i_patchelem, i_patchcnt, i_pth, i_eunknowns, i_niovars
    INTEGER, SAVE                                               :: i_timecount = 0
    CHARACTER (len=32)                                          :: c_file, c_mesh
    CHARACTER (len=64)                                          :: c_title, c_tmp
    INTEGER (KIND = GRID_SI), DIMENSION(:), ALLOCATABLE, TARGET :: i_eltlev
    TYPE (ncugrid_vardatatype), DIMENSION(:), ALLOCATABLE       :: p_vdataarr
    REAL (KIND = GRID_SR), DIMENSION(:,:), ALLOCATABLE, TARGET  :: r_vdata
    REAL (KIND = GRID_SR), DIMENSION(:,:), ALLOCATABLE          :: r_tmpdat
    REAL (KIND = GRID_SR), DIMENSION(:), ALLOCATABLE            :: r_vmean
    INTEGER(KIND = GRID_SI), DIMENSION(:), ALLOCATABLE          :: i_ndmaster
    INTEGER(KIND = GRID_SI), DIMENSION(:,:), ALLOCATABLE        :: i_ndpatch

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

    ALLOCATE(p_vdataarr(i_niovars+1), stat=i_alct)
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

    WRITE(c_tmp, '(A, A6)') TRIM(c_mesh), '_level'
    p_vdataarr(i_niovars+1)%c_varname       = TRIM(c_tmp)
    p_vdataarr(i_niovars+1)%c_long_name     = 'grid level'
    p_vdataarr(i_niovars+1)%c_standard_name = 'grid_level' ! not CF conform!
    p_vdataarr(i_niovars+1)%c_units         = 'None'
    p_vdataarr(i_niovars+1)%c_location      = 'face'
    p_vdataarr(i_niovars+1)%i_datatype      = 0

!--- allocate workspace
    i_eunknowns = GRID_femtypes%p_type(FEM_DG)%sig%i_unknowns

    ALLOCATE(r_tmpdat(i_niovars, i_numelmt*i_eunknowns), &
             r_vdata(i_niovars, i_numnode), &
             r_vmean(i_niovars), &
             i_eltlev(i_numelmt), &
             i_ndmaster(i_numnode), &
             i_ndpatch(GRID_patchelements, i_numnode), stat=i_alct)
    IF (i_alct /= 0) CALL grid_error(c_error='[plot_netcdf]: could not allocate data arrays')

!--- get arrays with some mesh data
    CALL grid_getinfo(p_ghand, l_finelevel=.TRUE., l_relative=.TRUE., &
                      i_elementlevel=i_eltlev, i_nodepatch=i_ndpatch, &
                      i_nodemaster=i_ndmaster)

!--- construct temporary data array
    r_tmpdat(1:i_nprogvars,:) = r_Q
    r_tmpdat(i_nprogvars+1:i_nprogvars+i_nsrcterms,:) = r_S
    IF (l_solplot) r_tmpdat(i_nprogvars+i_nsrcterms+1:i_niovars,:) = r_sol
    DO i_var = 1, i_niovars
      r_tmpdat(i_var,:) = r_tmpdat(i_var,:) * p_iovars(i_var)%r_factor
    END DO

!--- compute average over each node
    DO i_node = 1,i_numnode
      i_patchcnt  = 1_GRID_SI
      i_patchelem = i_ndpatch(i_patchcnt, i_node)
      r_vmean     = 0.0_GRID_SR

      patch_loop: DO WHILE(i_patchelem /= 0)
        DO i_pth=1, GRID_elementnodes
          IF(i_ndmaster(i_node) == i_ndmaster(i_elmtnodes(i_pth, i_patchelem))) THEN
            DO i_var = 1, i_niovars
              r_vmean(i_var) = r_vmean(i_var) + &
              r_tmpdat(i_var,i_elmtdofs(i_pth,i_patchelem))
            END DO
            EXIT
          END IF
        END DO

        ! loop increment
        i_patchcnt = i_patchcnt + 1_GRID_SI
        IF(i_patchcnt > GRID_patchelements) THEN
          CALL grid_error(1,'[plot_netcdf]: node patch overflow')
          EXIT patch_loop
        END IF
        i_patchelem = i_ndpatch(i_patchcnt, i_node)
      END DO patch_loop

      DO i_var = 1, i_niovars
        r_vdata(i_var,i_node) = r_vmean(i_var) / REAL(i_patchcnt-1, GRID_SR)
      END DO
    END DO

    DO i_node = 1,i_numnode
      r_vdata(:,i_node) = r_vdata(:,i_ndmaster(i_node))
    END DO

!--- fill structure with data
    DO i_var = 1, i_niovars
      p_vdataarr(i_var)%p_rvardata => r_vdata(i_var,:)
    END DO
    p_vdataarr(i_niovars+1)%p_ivardata => i_eltlev

!--- create the file for timestep data
    CALL ncugrid_createmesh(c_file, c_mesh, i_numnode, i_numelmt, &
                            i_elmtnodes, r_coonod, c_title, TRIM(GRID_parameters%author_affil2), &
                            i_numedge=i_numedge, i_edgenodes=i_edgenodes)

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
    DEALLOCATE(p_vdataarr, r_tmpdat, r_vdata, r_vmean, i_eltlev)

  END SUBROUTINE plot_netcdf

!*******************************************************************************
END MODULE IO_netcdfplot
