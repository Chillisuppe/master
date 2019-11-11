!*******************************************************************************
!
!> @file  IO_paraplot_dg.F90
!> @brief contains module IO_paraplot_dg
!
!*******************************************************************************
! MODULE DESCRIPTION:
!> @brief creates output in ASCII vtu file format
!
MODULE IO_paraplot_dg
  USE GRID_api
  USE FLASH_parameters
  USE DG_equation, ONLY : p_iovars, i_nprogvars, i_nsrcterms, VAR1D_PARACHECK

  PRIVATE

  INTEGER, SAVE :: i_timecounter= 0
  PUBLIC        :: plot_para_dg

  CONTAINS
!*******************************************************************************
! DESCRIPTION of [SUBROUTINE plot_para_dg]:
!> @brief create output in paraview compatible file format
!>
!> ON OUTPUT: a number of paraview compatible ascii files for each plotted time step
!>
!> @todo compute element averages such that they are also valid for higher order elements
!>
!> @param[in]         p_ghand         grid handling data structure
!> @param[in]         p_param         control structure for global parameters
!> @param[in]         i_time          time stamp for file naming (opt.)
!> @param[in]         r_time          current model time
!> @param[in]         i_faceunknowns  number of DOFs per element
!> @param[in]         i_degree        polynomial degree of element
!> @param[in]         i_numelmt       total number of elements
!> @param[in]         i_numnode       total number of nodes
!> @param[in]         i_elmtnodes     node indices of each element; DIMENSION(3, i_numelmt)
!> @param[in]         i_elmtdofs      DOF indices corresponding to each element; DIMENSION(i_faceunknowns, i_numelmt)
!> @param[in]         r_coonod        coordinates for each grid node
!> @param[in]         r_coodof        coordinates for each dof in the grid
!> @param[in]         r_Q             discrete solution vector
!> @param[in]         r_S             discrete fields used in source terms
!> @param[in]         r_sol           discrete representation of exact solution vector
!> @param[in]         l_solplot       switch for output of exact solution
!
  SUBROUTINE plot_para_dg(p_ghand, p_param, i_time, r_time, i_faceunknowns, &
                          i_degree, i_numelmt, i_numnode, i_elmtnodes, i_elmtdofs, &
                          r_coonod, r_coodof, r_Q, r_S, r_sol, l_solplot)

    IMPLICIT NONE

    TYPE (grid_handle),                       INTENT(IN)        :: p_ghand
    TYPE (control_struct),                    INTENT(IN)        :: p_param
    INTEGER (KIND = GRID_SI), OPTIONAL,       INTENT(IN)        :: i_time
    REAL (KIND = GRID_SR),                    INTENT(IN)        :: r_time
    INTEGER (KIND = GRID_SI),                 INTENT(IN)        :: i_faceunknowns, i_degree, &
                                                                   i_numelmt, i_numnode
    INTEGER (KIND = GRID_SI), DIMENSION(:,:), INTENT(IN)        :: i_elmtnodes, i_elmtdofs
    REAL (KIND = GRID_SR),    DIMENSION(:,:), INTENT(IN)        :: r_coonod, r_coodof
    REAL (KIND = GRID_SR),    DIMENSION(:,:), INTENT(IN)        :: r_Q, r_S
    REAL (KIND = GRID_SR),    DIMENSION(:,:), INTENT(IN)        :: r_sol
    LOGICAL,                                  INTENT(IN)        :: l_solplot

!--- local declarations
    INTEGER (KIND = GRID_SI)                                    :: i_iounit, i_elmt, &
                                                                   i_alct, i_numdofs, i_tcnt, &
                                                                   i_fst, &
                                                                   i_femtyp, i_dof, i_var, i_niovars
    INTEGER (KIND = GRID_SI), DIMENSION(i_faceunknowns)         :: i_edofsort
    CHARACTER (len=32)                                          :: c_file
    REAL (KIND = GRID_SR), DIMENSION(:,:), ALLOCATABLE          :: r_tmpdat

!--- file handling (open)

    IF(PRESENT(i_time)) THEN
      i_tcnt = i_time
    ELSE
      i_tcnt = i_timecounter
    END IF
    i_timecounter = i_timecounter+1

    IF (l_solplot) THEN
      i_niovars = 2 * i_nprogvars + i_nsrcterms
    ELSE
      i_niovars = i_nprogvars + i_nsrcterms
    END IF

    WRITE(c_file, '(A, I8.8, A4)') TRIM(GRID_parameters%program_name), i_tcnt, '.vtu'

    OPEN(NEWUNIT=i_iounit, file=c_file, form='formatted', iostat=i_fst)

    file_notopen: IF(i_fst /= 0) THEN
      write(0,*) 'ERROR: Filename: ', c_file
      IF(GRID_parameters%iolog > 0) &
        write(GRID_parameters%iolog,*) 'ERROR: Filename: ', c_file
      CALL grid_error(c_error='[plot_netcdf]: could not read input file')
    ELSE file_notopen
      IF(GRID_parameters%iolog > 0) THEN
        write(GRID_parameters%iolog,*) 'INFO: Filename: ', c_file, ' opened on unit: ', i_iounit
      END IF
    END IF file_notopen

!--- allocate workspace
    i_numdofs   = p_ghand%i_unknowns(FEM_DG) !number of dofs

    ALLOCATE(r_tmpdat(i_niovars, i_numdofs), stat=i_alct)
    IF (i_alct /= 0) CALL grid_error(c_error='[plot_netcdf]: could not allocate data arrays')

!--- construct temporary data array
    r_tmpdat(1:i_nprogvars,:) = r_Q
    r_tmpdat(i_nprogvars+1:i_nprogvars+i_nsrcterms,:) = r_S
    IF (l_solplot) r_tmpdat(i_nprogvars+i_nsrcterms+1:i_niovars,:) = r_sol
    DO i_var = 1, i_niovars
      r_tmpdat(i_var,:) = r_tmpdat(i_var,:) * p_iovars(i_var)%r_factor
    END DO

!--- check for consistency
    i_femtyp = grid_femvarquery(VAR1D_PARACHECK)
    IF(i_femtyp /= FEM_DG) CALL grid_error(c_error='[plot_para_dg]: Inconsistency with fem types')

!--- write HEADER
    WRITE(i_iounit,'(A)') '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'
    WRITE(i_iounit,'(A)') '<UnstructuredGrid>'
    WRITE(i_iounit,'(A, I8, A, I8, A)') '<Piece NumberOfPoints="',i_numdofs,'" NumberOfCells="',i_numelmt,'">'

!--- write POINTS
    WRITE(i_iounit,'(A)') '<Points>'
    IF (GRID_DIMENSION == 2) THEN
      WRITE(i_iounit,1030) GRID_DIMENSION+1
      DO i_dof=1,i_numdofs
       WRITE(i_iounit,*) r_coodof(:,i_dof), 0.0
      END DO
      WRITE(i_iounit,1031)
    ELSE IF (GRID_DIMENSION == 3) THEN
      WRITE(i_iounit,1030) GRID_DIMENSION
      DO i_dof=1,i_numdofs
       WRITE(i_iounit,*) r_coodof(:,i_dof)
      END DO
      WRITE(i_iounit,1031)
    END IF

    WRITE(i_iounit,'(A)') '</Points>'

!--- write CELLS
    WRITE(i_iounit,'(A)') '<Cells>'

    SELECT CASE(i_degree)
      CASE (1)
        i_edofsort = [1,2,3]
      CASE (2)
        i_edofsort = [1,2,3,4,6,5]
      CASE DEFAULT
        i_edofsort = 1
    END SELECT

    WRITE(i_iounit,1033) 'connectivity'
    DO i_elmt=1,i_numelmt
      WRITE(i_iounit,*) i_elmtdofs(i_edofsort,i_elmt)-1     !statt enodes
    END DO
    WRITE(i_iounit,1031)

    WRITE(i_iounit,1033) 'offsets'
    DO i_elmt=1,i_numelmt
      WRITE(i_iounit,*) i_faceunknowns*i_elmt
    END DO
    WRITE(i_iounit,1031)

    WRITE(i_iounit,1034) 'types' !(types=5 is triangle, types=22 is quadratic triangle)
    SELECT CASE(i_degree)
      CASE (1)
        DO i_elmt=1,i_numelmt
          WRITE(i_iounit,*) '5'
        END DO
      CASE (2)
        DO i_elmt=1,i_numelmt
          WRITE(i_iounit,*) '22'
        END DO
    END SELECT
    WRITE(i_iounit,1031)
    WRITE(i_iounit,'(A)') '</Cells>'

!--- write POINTDATA
    WRITE(i_iounit,'(A)') '<PointData>'
    WRITE(i_iounit,1033) 'nodeindex'
    DO i_dof=1,i_numdofs
      WRITE(i_iounit,*) i_dof
    END DO
    WRITE(i_iounit,1031)

    DO i_var = 1, i_niovars
      WRITE(i_iounit,1036) TRIM(p_iovars(i_var)%c_varname)
      DO i_dof=1,i_numdofs
        WRITE(i_iounit,*) r_tmpdat(i_var,i_dof)
      END DO
      WRITE(i_iounit,1031)
    END DO !i_var

    WRITE(i_iounit,'(A)') '</PointData>'

!--- write CELLDATA
    WRITE(i_iounit,'(A)') '<CellData>'
    WRITE(i_iounit,1033) 'elementindex'
    DO i_elmt=1,i_numelmt
      WRITE(i_iounit,*) i_elmt
    END DO
    WRITE(i_iounit,1031)

    WRITE(i_iounit,'(A)') '</CellData>'

!--- write FOOTER
    WRITE(i_iounit,'(A)') '</Piece>'
    WRITE(i_iounit,'(A)') '</UnstructuredGrid>'
    WRITE(i_iounit,'(A)') '</VTKFile>'

!--- close file
    CLOSE(i_iounit)
    IF(GRID_parameters%iolog > 0) &
      WRITE(GRID_parameters%iolog,*) 'INFO: Closed file on unit: ', i_iounit

    DEALLOCATE(r_tmpdat)

1030    FORMAT('<DataArray type="Float32" NumberOfComponents="', i2, '" Format="ascii">')
1032    FORMAT('<DataArray type="Float32" Name="', A, '" Format="ascii">')
1033    FORMAT('<DataArray type="Int32" Name="', A, '" Format="ascii">')
1036    FORMAT('<DataArray type="Float32" Name="', A, '" Format="ascii">')
1035    FORMAT('<DataArray type="Float32" NumberOfComponents="3" Name="', A, '" Format="ascii">')
1034    FORMAT('<DataArray type="UInt8" Name="', A, '" Format="ascii">')
1031    FORMAT('</DataArray>')

  END SUBROUTINE plot_para_dg

!*******************************************************************************
END MODULE IO_paraplot_dg
