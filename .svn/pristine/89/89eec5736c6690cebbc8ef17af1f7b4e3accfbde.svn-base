!*******************************************************************************
!
!> @file  TEST_eqsource_readcoards.f90
!> @brief program to test IO of COARDS-files with the MISC_eqsource modules
!>        (this is part 2 of 2 and reads the netCDF files)
!
!*******************************************************************************
PROGRAM TEST_eqsource_readcoards

  USE GRID_api
  USE IO_ncugrid
  USE MISC_eqsource

  IMPLICIT NONE

  INTEGER (KIND = GRID_SI)                                  :: FEM_DG, VAR1D_DISP
  INTEGER (KIND = GRID_SI)                                  :: i_alct, i_count, &
    i_faceunknowns, i_numdofs, i_numelmt, i_numnode, i_reflevel
  INTEGER, DIMENSION(GRID_elementnodes)                     :: i_ndindices
  INTEGER (KIND = GRID_SI), DIMENSION(:,:), ALLOCATABLE     :: i_elmtnodes, i_elmtdofs, &
    i_elmtnodesdg
  REAL (KIND = GRID_SR), DIMENSION(:), ALLOCATABLE          :: r_uplift, r_vdata
  REAL (KIND = GRID_SR), DIMENSION(:,:), ALLOCATABLE        :: r_coonod, r_coodof, r_coonoddg

!--- set maximum refinement level
  i_reflevel = 14

!--- set FLASH description in global datastruct
  GRID_parameters%program_name = 'TEST_eqsource_readcoards'

!--- register the fem type
  FEM_DG = grid_registerfemtype('dg_linear.ftf')

!--- initialize grid generator
  CALL grid_initialize(i_logging=7, i_initvals=1_GRID_SI)

!--- register variable
  VAR1D_DISP = grid_registerfemvar(FEM_DG, i_length=1)

!--- initialize grid parameters
  CALL grid_setparameter(p_grid, i_coarselevel=i_reflevel, i_finelevel=i_reflevel)

!--- create initial triangulation
  CALL grid_createinitial(p_grid, c_filename='Triang-Sumatra04fine.dat')

!--- compute some constants
  i_faceunknowns = GRID_femtypes%p_type(FEM_DG)%sig%i_unknowns
  i_numelmt      = p_grid(i_timeplus)%i_enumfine
  i_numnode      = p_grid(i_timeplus)%i_nnumber
  i_numdofs      = i_numelmt*i_faceunknowns

!--- allocate workspace
  ALLOCATE(i_elmtnodes(3,i_numelmt), &
           i_elmtdofs(i_faceunknowns, i_numelmt), &
           i_elmtnodesdg(GRID_elementnodes, i_numelmt), &
           r_coonod(GRID_dimension,i_numnode), &
           r_coodof(GRID_dimension,i_numdofs), &
           r_coonoddg(GRID_dimension, i_numelmt*GRID_elementnodes), &
           r_uplift(i_numdofs), r_vdata(i_numelmt*GRID_elementnodes), stat=i_alct)
  IF(i_alct /= 0) &
    CALL grid_error(c_error='[TEST_eqsource_readcoards]: Could not allocate workspace')

!--- get grid information
  CALL grid_getinfo(p_grid(i_timeplus), i_femtype=FEM_DG, l_relative=.TRUE., &
                    l_finelevel=.TRUE., i_elementnodes=i_elmtnodes, &
                    i_elementdofs=i_elmtdofs, r_nodecoordinates=r_coonod, &
                    r_dofcoordinates=r_coodof)

!--- load fault data and compute generated uplift on the grid
  CALL initialize_source('TEST_dz_coards.nc')

  DO i_count = 1, i_numdofs
    r_uplift(i_count) = compute_uplift(r_coodof(:,i_count), 0.0_GRID_SR)
  END DO

!--- compute remeshed data for NetCDF output
  DO i_count = 1, i_numelmt
    i_ndindices               = (i_count-1)*3 + [1, 2, 3]
    i_elmtnodesdg(:,i_count)   = i_ndindices
    r_coonoddg(:,i_ndindices) = r_coonod(:,i_elmtnodes(:,i_count))
    r_vdata(i_ndindices)      = r_uplift(i_elmtdofs(:,i_count))
  END DO

!--- create the the NetCDF file
!#ifdef USE_NETCDF
  CALL ncugrid_createmesh('TEST_eqsource_readcoards.nc', 'Mesh2', i_numelmt*GRID_elementnodes, i_numelmt, &
                          i_elmtnodesdg, r_coonoddg, 'none', 'none')

  CALL ncugrid_putvariable('TEST_eqsource_readcoards.nc', 'Mesh2_disp', 'displacement', 'displacement', 'm', &
                           'Mesh2', 'node', r_vdata, 0.0_GRID_SR)

!#else
!  CALL grid_error(c_error='[TEST_eqsource_readcoards]: NetCDF library not installed - cannot write data!')
!#endif

!--- free memory, terminate grid, etc.
  DEALLOCATE(i_elmtnodes, i_elmtdofs, i_elmtnodesdg, r_coonod, r_coodof, &
             r_coonoddg, r_uplift, r_vdata)

  CALL grid_terminate

  WRITE(*,'(A)') 'Output has been written. Now terminating.'
  STOP

END PROGRAM TEST_eqsource_readcoards
