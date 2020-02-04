!*****************************************************************
!
!> @file MISC_diag_mesh_efficiency.F90
!> @brief contains module MISC_diag
!
!*****************************************************************
!
! VERSION(S):
! 1. original version            j. behrens    8/2009
! 2. Adapted for StormFlash2d    n. beisiegel  7/2014
!
!*****************************************************************
! MODULE DESCRIPTION:
!> contains the diagnostics function
!
! Verändere nur SUBROUTINE diag_diagnostics
!
MODULE MISC_diag

  USE GRID_api
  USE FLASH_parameters
  USE DG_equation
  USE DG_initial
  USE MISC_utils
  USE MISC_timing

  INTEGER (KIND = GRID_SI)                              :: i_numdofs, i_iodiag, i_file
  REAL (KIND = GRID_SR)                                 :: r_mass0

  PRIVATE
  PUBLIC  :: diag_initialize, diag_quit, diag_syncdisk, diag_diagnostics

  CONTAINS
!*******************************************************************************
! DESCRIPTION of [SUBROUTINE diag_initialize]:
!> @brief initializes diagnostics to the computed solution
!>
!> @param[in]       p_ghand         grid handling data structure
!> @param[in,out]   p_param         control structure for global parameters
!> @param[in]       i_faceunknowns  number of DOFs per element
!> @param[in]       r_einvvand      inverse Vandermonde matrix
!
  SUBROUTINE diag_initialize(p_ghand, p_param, i_faceunknowns, r_einvvand)

    IMPLICIT NONE

    TYPE (grid_handle), INTENT(in)                                :: p_ghand
    TYPE (control_struct), INTENT(inout)                          :: p_param
    INTEGER (KIND = GRID_SI), INTENT(in)                          :: i_faceunknowns
    REAL (KIND = GRID_SR), DIMENSION(:,:), INTENT(in)             :: r_einvvand

!--- local declarations
    CHARACTER (len=32)                                            :: c_file
    INTEGER (KIND = GRID_SI)                                      :: i_fst, i_tmp

    i_numdofs = i_faceunknowns

!--- compute initial mass
    r_mass0 = grid_dg_globalnorm_l1(p_ghand, VAR1D_HZERO)

!--- open file for diagnostic output
    i_tmp = p_param%num%i_experiment
    WRITE(c_file,'(A, A6, I4.4)') TRIM(GRID_parameters%program_name), '_diag.', i_tmp
    OPEN(NEWUNIT=i_iodiag, file=c_file, action='write', form='formatted', iostat=i_fst)

    not_opened: IF(i_fst /= 0) THEN
      CALL grid_error(c_error='[diag_initialize]: could not open general diagnostics file')
    END IF not_opened

    IF(GRID_parameters%iolog > 0) &
      WRITE(GRID_parameters%iolog,*) 'INFO: Filename: ', c_file, ' opened on unit: ', i_iodiag

!--- print some version info and header
    WRITE(i_iodiag,1100) GRID_parameters%program_name, GRID_parameters%version, &
                         GRID_parameters%subversion, GRID_parameters%patchversion, &
                         p_param%num%i_crslevel, p_param%num%i_reflevel

    WRITE(i_iodiag,'(A8 ,1x)', advance='no') 'timestep            '
    WRITE(i_iodiag,'(A14,1x)', advance='no') 'modeltime           '
    WRITE(i_iodiag,'(A14,1x)', advance='no') 'cpu time            '
    WRITE(i_iodiag,'(A14,1x)', advance='no') 'cfl                 '
    WRITE(i_iodiag,'(A8 ,1x)', advance='no') 'elements            '
    WRITE(i_iodiag,'(A8 ,1x)', advance='no') 'edges               '
    WRITE(i_iodiag,'(A8 ,1x)', advance='no') 'nodes               '
    WRITE(i_iodiag,'(A14,1x)', advance='no') 'min.edge len.       '
    WRITE(i_iodiag,'(A14,1x)', advance='no') 'max.edge len.       '
    WRITE(i_iodiag,'(A20,1x)', advance='no') 'relative mass       '
    WRITE(i_iodiag,'(A14,1x)', advance='no') 'L2 err h            '
    WRITE(i_iodiag,'(A14,1x)', advance='no') 'L1 err h            '
    WRITE(i_iodiag,'(A14,1x)', advance='no') 'Linf err h          '
    WRITE(i_iodiag,'(A14,1x)', advance='no') 'L2 err hu           '
    WRITE(i_iodiag,'(A14,1x)', advance='no') 'L1 err hu           '
    WRITE(i_iodiag,'(A14,1x)', advance='no') 'Linf err hu         '
    WRITE(i_iodiag,*)

 1100 FORMAT('************************************************************', &
             '************************************************************', &
             '************************************************************', &
             '*********',/ &
             '***************** PROGRAM: ',a15,62x,'**********************', &
             '***************************************************************',/ &
             '***************** VERSION: ',i2.2,'.',i2.2,'.',i2.2,69x,'***', &
             '************************************************************', &
             '**********************',/ &
             '***************** Diagnostic output ',68x,'*****************', &
             '************************************************************', &
             '********',/ &
             '***************** Refinement level  ',i6,i6, 56x '**********', &
             '************************************************************', &
             '***************',/ &
             '************************************************************', &
             '************************************************************', &
             '*********', &
             '************************************************************')

  END SUBROUTINE diag_initialize

!*******************************************************************************
! DESCRIPTION of [SUBROUTINE diag_quit]:
!> @brief exits diagnostics to the computed solution
!>
!> @param[in]       p_ghand     grid handling data structure
!> @param[in,out]   p_param     control structure for global parameters
!
  SUBROUTINE diag_quit(p_ghand, p_param)

    IMPLICIT NONE

    TYPE (grid_handle), INTENT(in)                                :: p_ghand
    TYPE (control_struct), INTENT(inout)                          :: p_param

!--- write closing line to diagnostic output file
    WRITE(i_iodiag, 1200)

!--- close diagnostic output file
    CLOSE(i_iodiag)
    IF(GRID_parameters%iolog > 0) &
      WRITE(GRID_parameters%iolog,*) 'INFO: Closed file on unit: ', i_iodiag

 1200 FORMAT('************************************************', &
             '************************************************', &
             '************************************************', &
             '*********************************************')

  END SUBROUTINE diag_quit

!*******************************************************************************
! DESCRIPTION of [SUBROUTINE diag_syncdisk]:
!> @brief synchronizes diagnostics file to disk
!
  SUBROUTINE diag_syncdisk

    IMPLICIT NONE

!--- local declarations
    INTEGER (KIND = GRID_SI)                                      :: i_fst

!--- flush and sync
    FLUSH(i_iodiag)
    i_fst = fsync(FNUM(i_iodiag))

    IF (i_fst /= 0) CALL grid_error(c_error='[diag_syncdisk]: Error calling FSYNC')

  END SUBROUTINE diag_syncdisk

!*******************************************************************************
! DESCRIPTION of [SUBROUTINE diag_diagnostics]:
!> @brief applies diagnostics to the computed solution
!>
!> @param[in]       p_ghand     grid handling data structure
!> @param[in,out]   p_param     control structure for global parameters
!> @param[in]       p_tinfo     timestep info structure
!> @param[in]       i_elmtnodes     global element-node relation matrix
!> @param[in]       i_elmtdofs      global element-dof relation matrix
!> @param           r_coonod        coordinates for each grid node
!> @param           r_coodof        coordinates for each dof in the grid
!> @param           r_Q             discrete solution vector
!> @param           r_S             discrete fields used in source terms
!>
!> @note ON OUTPUT: plot output
!
  SUBROUTINE diag_diagnostics(p_ghand, p_param, p_tinfo, i_elmtnodes, i_elmtdofs, &
                              r_coonod, r_coodof, r_Q, r_S)

    IMPLICIT NONE

    TYPE (grid_handle), INTENT(inout)                             :: p_ghand
    TYPE (control_struct), INTENT(inout)                          :: p_param
    TYPE (rt_info), INTENT(in)                                    :: p_tinfo
    INTEGER (KIND = GRID_SI), INTENT(in), DIMENSION(:,:)          :: i_elmtnodes, i_elmtdofs
    REAL (KIND = GRID_SR), INTENT(in), DIMENSION(:,:)             :: r_coonod, r_coodof, r_Q, r_S

!--- local declarations
    INTEGER (KIND = GRID_SI)                                      :: i_alct, i_numelmt, &
      i_equadpts, i_elmt, i_quad
    REAL (KIND = GRID_SR)                                         :: r_time, r_relmass, &
      r_medln, r_mxedln, r_err_h_l2, r_err_mom_l2, r_err_h_linf, r_err_mom_linf, &
      r_err_h_l1, r_err_mom_l1
    REAL (KIND = GRID_SR), DIMENSION(:), ALLOCATABLE      :: r_diam
    REAL (KIND = GRID_SR), DIMENSION(:,:), ALLOCATABLE    :: r_epsi, r_sol
    REAL (KIND = GRID_SR), DIMENSION(:,:,:), ALLOCATABLE  :: r_solq

!--- initialize some values
    r_time     = p_tinfo%r_modeltime
    i_numelmt  = p_ghand%i_enumfine
    i_equadpts = GRID_femtypes%p_type(FEM_DG)%sig%i_equadpts

!--- allocate workspace
    ALLOCATE(r_epsi(i_numdofs, i_equadpts), &
             r_diam(i_numelmt))

    r_epsi = GRID_femtypes%p_type(FEM_DG)%sig%r_epsiquad

!--- get grid information
    CALL grid_getinfo(p_ghand, i_femtype=FEM_DG, r_elementwidth=r_diam)

!--- get minimum edge length
    r_medln  = MINVAL(r_diam)
    r_mxedln = MAXVAL(r_diam)
!    CALL grid_edgelength(p_ghand, r_min=r_medln, r_max=r_mxedln)

!--- compute mass (integral over fluid depth) and relative mass
    r_relmass = grid_dg_globalnorm_l1(p_ghand, VAR1D_HZERO) / r_mass0

!   Hier Dimension angepasst 3 -> i_nprogvars
    ALLOCATE(r_sol(i_nprogvars,p_ghand%i_enumfine*i_numdofs), &
             r_solq(i_nprogvars,p_ghand%i_enumfine,i_equadpts))

    ! compute analytical solution at quadrature points
    CALL dg_elmt_solution(r_coodof, r_time, i_numelmt, i_elmtdofs, &
                          r_coonod, i_elmtnodes, r_sol)
    DO i_elmt=1, i_numelmt
      r_solq(:,i_elmt,:) = MATMUL(r_sol(:,i_elmtdofs(:,i_elmt)), r_epsi) * GRID_GRAV
    END DO

    r_err_h_l2     = grid_dg_globalnorm_l2(p_ghand, VAR1D_HZERO, r_solq(1,:,:)) / &
                     grid_dg_globalnorm_l2(p_ghand, VAR1D_HZERO)
    r_err_h_l1     = grid_dg_globalnorm_l1(p_ghand, VAR1D_HZERO, r_solq(1,:,:)) / &
                     grid_dg_globalnorm_l1(p_ghand, VAR1D_HZERO)
    r_err_h_linf   = grid_dg_globalnorm_linf(p_ghand, VAR1D_HZERO, r_solq(1,:,:)) / &
                     grid_dg_globalnorm_linf(p_ghand, VAR1D_HZERO)

    r_err_mom_l2   = grid_dg_globalnorm_l2(p_ghand, VAR3D_MZERO, r_solq(2,:,:)) / &
                     grid_dg_globalnorm_l2(p_ghand, VAR3D_MZERO) + &
                     grid_dg_globalnorm_l2(p_ghand, VAR3D_MZERO+1, r_solq(3,:,:)) / &
                     grid_dg_globalnorm_l2(p_ghand, VAR3D_MZERO+1) + &
                     grid_dg_globalnorm_l2(p_ghand, VAR3D_MZERO+2, r_solq(4,:,:)) / &
                     grid_dg_globalnorm_l2(p_ghand, VAR3D_MZERO+2)
    r_err_mom_l1   = grid_dg_globalnorm_l1(p_ghand, VAR3D_MZERO, r_solq(2,:,:)) / &
                     grid_dg_globalnorm_l1(p_ghand, VAR3D_MZERO) + &
                     grid_dg_globalnorm_l1(p_ghand, VAR3D_MZERO+1, r_solq(3,:,:)) / &
                     grid_dg_globalnorm_l1(p_ghand, VAR3D_MZERO+1) + &
                     grid_dg_globalnorm_l1(p_ghand, VAR3D_MZERO+2, r_solq(4,:,:)) / &
                     grid_dg_globalnorm_l1(p_ghand, VAR3D_MZERO+2)
    r_err_mom_linf = grid_dg_globalnorm_linf(p_ghand, VAR3D_MZERO, r_solq(2,:,:)) / &
                     grid_dg_globalnorm_linf(p_ghand, VAR3D_MZERO) + &
                     grid_dg_globalnorm_linf(p_ghand, VAR3D_MZERO+1, r_solq(3,:,:)) / &
                     grid_dg_globalnorm_linf(p_ghand, VAR3D_MZERO+1) + &
                     grid_dg_globalnorm_linf(p_ghand, VAR3D_MZERO+2, r_solq(4,:,:)) / &
                     grid_dg_globalnorm_linf(p_ghand, VAR3D_MZERO+2)

    DEALLOCATE(r_sol, r_solq)

!--- print results to diagnostics file
    WRITE(i_iodiag,'(i8,1x)'    , advance='no') p_tinfo%i_step
    WRITE(i_iodiag,'(e14.7,1x)' , advance='no') r_time
    WRITE(i_iodiag,'(e14.7,1x)' , advance='no') p_tinfo%r_cputime
    WRITE(i_iodiag,'(e14.7,1x)' , advance='no') p_tinfo%r_cflnumber
    WRITE(i_iodiag,'(i8,1x)'    , advance='no') p_ghand%i_enumfine
    WRITE(i_iodiag,'(i8,1x)'    , advance='no') p_ghand%i_gnumfine
    WRITE(i_iodiag,'(i8,1x)'    , advance='no') p_ghand%i_nnumber
    WRITE(i_iodiag,'(e14.7,1x)' , advance='no') r_medln
    WRITE(i_iodiag,'(e14.7,1x)' , advance='no') r_mxedln
    WRITE(i_iodiag,'(e20.13,1x)', advance='no') r_relmass
    WRITE(i_iodiag,'(e14.7,1x)' , advance='no') r_err_h_l2
    WRITE(i_iodiag,'(e14.7,1x)' , advance='no') r_err_h_l1
    WRITE(i_iodiag,'(e14.7,1x)' , advance='no') r_err_h_linf
    WRITE(i_iodiag,'(e14.7,1x)' , advance='no') r_err_mom_l2
    WRITE(i_iodiag,'(e14.7,1x)' , advance='no') r_err_mom_l1
    WRITE(i_iodiag,'(e14.7,1x)' , advance='no') r_err_mom_linf
    WRITE(i_iodiag,*)

!--- deallocate workspace
    DEALLOCATE(r_epsi, r_diam)

  END SUBROUTINE diag_diagnostics

END MODULE MISC_diag
