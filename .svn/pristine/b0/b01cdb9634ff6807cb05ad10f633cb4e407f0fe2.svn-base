!*******************************************************************************
!
!> @file  DG_equation_swesrcbathy.F90
!> @brief contains module DG_equation
!
!*******************************************************************************
!
! VERSION(s):
! 1. first version                            s.beckers    07/2013
! 2. include dummies for DWR                  s.beckers    06/2014
!
!*******************************************************************************
! MODULE DESCRIPTION:
!> @brief Defines all quantities with respect to the equations used, i.e.,
!>        in this case the shallow water equations. Included source terms are
!>        bottom topography (explicit) and Manning friction (split implicit via
!>        pre_timestepping).
!
MODULE DG_equation

  USE GRID_api
  USE FLASH_parameters
  USE IO_equation, ONLY : p_equationparam

  PRIVATE
  PUBLIC  i_nprogvars, i_nsrcterms, i_valQ, i_valS, &
          VAR2D_MZERO, VAR1D_HZERO, VAR1D_HCNST, VAR1D_PARACHECK, &
          velocity, equation_initialize, divflux, flux, fluxGrav, source, &
          io_initialize, pre_timestepping, p_iovars

  INTEGER (KIND = GRID_SI), PARAMETER               :: i_nprogvars = 3 !< number of unknowns
  INTEGER (KIND = GRID_SI), PARAMETER               :: i_nsrcterms = 1 !< number of grid source term variables

  INTEGER (KIND = GRID_SI), DIMENSION(i_nprogvars)  :: i_valQ !< array of index handles for registered unknowns
  INTEGER (KIND = GRID_SI), DIMENSION(i_nsrcterms)  :: i_valS !< array of index handles for registered source term variables

!--- constants for the DG variables
  INTEGER (KIND = GRID_SI)                          :: VAR2D_MZERO     !< momentum variable
  INTEGER (KIND = GRID_SI)                          :: VAR1D_HZERO     !< height variable
  INTEGER (KIND = GRID_SI)                          :: VAR1D_HCNST     !< bathy variable
  INTEGER (KIND = GRID_SI)                          :: VAR1D_PARACHECK !< variable for consistency check in io_paraplot

!--- array for generic IO
  TYPE (io_vars), DIMENSION(2 * i_nprogvars + i_nsrcterms)  :: p_iovars

  INTERFACE velocity
    MODULE PROCEDURE velocity_scalar, velocity_array
  END INTERFACE

  CONTAINS

!*******************************************************************************
  FUNCTION velocity_scalar(r_h, r_hu) RESULT(r_u)

    IMPLICIT NONE

    REAL (KIND = GRID_SR), INTENT(in)               :: r_h, r_hu
    REAL (KIND = GRID_SR)                           :: r_u

    r_u = MERGE(r_hu/r_h, 0.0_GRID_SR, r_h > p_equationparam%r_wettol)

  END FUNCTION velocity_scalar

!*******************************************************************************
  FUNCTION velocity_array(r_h, r_hu) RESULT(r_u)

    IMPLICIT NONE

    REAL (KIND = GRID_SR), DIMENSION(:), INTENT(in) :: r_h, r_hu
    REAL (KIND = GRID_SR), DIMENSION(SIZE(r_h))     :: r_u

    r_u = MERGE(r_hu/r_h, 0.0_GRID_SR, r_h > p_equationparam%r_wettol)

  END FUNCTION velocity_array

!*******************************************************************************
! DESCRIPTION of [SUBROUTINE equation_initialize]:
!> @brief initializes grid variables for unknowns and source term variables
!>
!> @param[in]       i_femtype   index handle for fem type
!
  SUBROUTINE equation_initialize(i_femtype)

    IMPLICIT NONE

    INTEGER (KIND = GRID_SI), INTENT(in)            :: i_femtype

    VAR2D_MZERO = grid_registerfemvar(i_femtype, i_length=2)
    VAR1D_HZERO = grid_registerfemvar(i_femtype)
    VAR1D_HCNST = grid_registerfemvar(i_femtype)

    i_valQ = (/ VAR1D_HZERO, VAR2D_MZERO, VAR2D_MZERO+1 /)
    i_valS = (/ VAR1D_HCNST /)

  END SUBROUTINE equation_initialize

!*******************************************************************************
! DESCRIPTION of [FUNCTION flux]:
!> @brief Routine for flux computation. Has to be defined by the user.
!>
!> @param[in]       r_Q         state vector (fluid depth, momentum in x-, momentum in y-direction)
!> @param[in]       r_S         source terms
!> @return                      return value
!
  FUNCTION flux(r_Q, r_S, r_iswet) RESULT(r_flux)

    IMPLICIT NONE

    REAL (KIND = GRID_SR), DIMENSION(i_nprogvars), INTENT(in) :: r_Q
    REAL (KIND = GRID_SR), DIMENSION(i_nsrcterms), INTENT(in), OPTIONAL :: r_S
    REAL (KIND = GRID_SR), INTENT(in), OPTIONAL               :: r_iswet
    REAL (KIND = GRID_SR), DIMENSION(i_nprogvars,2)           :: r_flux

    REAL (KIND = GRID_SR)                                     :: r_u, r_v

    r_u = velocity(r_Q(1), r_Q(2))
    r_v = velocity(r_Q(1), r_Q(3))

    r_flux(1,1) = r_Q(2)
    r_flux(1,2) = r_Q(3)

    r_flux(2,1) = r_Q(2)*r_u
    r_flux(2,2) = r_Q(2)*r_v

    r_flux(3,1) = r_Q(3)*r_u
    r_flux(3,2) = r_Q(3)*r_v

    IF (.NOT. PRESENT(r_iswet) .OR. r_iswet==1.0) THEN
      r_flux(2,1) = r_flux(2,1) + 0.5_GRID_SR*r_Q(1)**2
      r_flux(3,2) = r_flux(3,2) + 0.5_GRID_SR*r_Q(1)**2
    END IF

  END FUNCTION flux

!*******************************************************************************
  FUNCTION fluxGrav(r_h) RESULT(r_flux)

    IMPLICIT NONE

    REAL (KIND = GRID_SR), INTENT(in)               :: r_h
    REAL (KIND = GRID_SR), DIMENSION(i_nprogvars,2) :: r_flux

    r_flux = 0.0_GRID_SR

    r_flux(2,1) = 0.5_GRID_SR*r_h**2
    r_flux(3,2) = 0.5_GRID_SR*r_h**2

  END FUNCTION fluxGrav

!*******************************************************************************
! DESCRIPTION of [SUBROUTINE source]:
!> @brief
!
!> @param         r_rhs
!> @param         i_faceunknowns
!> @param         i_equadpts
!> @param         r_eqwei
!> @param         r_eMinvpsi
!> @param         r_epsi
!> @param         r_ddxq
!> @param         r_ddyq
!> @param         r_lapq
!> @param         r_Qelmt
!> @param         r_Selmt
!> @param         r_coo
!> @param         r_time
!> @param         r_iswet
!
  SUBROUTINE source(r_rhs, i_faceunknowns, i_equadpts, r_eqwei, r_eMinvpsi, r_epsi, &
                    r_ddxq, r_ddyq, r_lapq, r_Qelmt, r_Selmt, r_coo, r_time, r_iswet)

    IMPLICIT NONE

    REAL (KIND = GRID_SR), DIMENSION(:,:), INTENT(inout)    :: r_rhs
    INTEGER (KIND = GRID_SI), INTENT(in)                    :: i_faceunknowns, i_equadpts
    REAL (KIND = GRID_SR), DIMENSION(:), INTENT(in)         :: r_eqwei
    REAL (KIND = GRID_SR), DIMENSION(:,:), INTENT(in)       :: r_eMinvpsi, &
                                                               r_epsi, r_ddxq, r_ddyq, r_lapq, &
                                                               r_Qelmt, r_Selmt, r_coo
    REAL (KIND = GRID_SR), INTENT(in)                       :: r_time, r_iswet

    INTEGER (KIND = GRID_SI)                                :: i_quad, i_dof
    REAL (KIND = GRID_SR), DIMENSION(i_equadpts)            :: r_h_q, r_hu_q, r_hv_q, &
                                                               r_cor_q, r_laphu_q, r_laphv_q
    REAL (KIND = GRID_SR), DIMENSION(i_equadpts,2)          :: r_gradb_q, r_tau_q
    REAL (KIND = GRID_SR), DIMENSION(i_nprogvars)           :: r_source

!--- interpolate to quadrature points and compute derivatives
    r_h_q          = MATMUL(r_Qelmt(1,:), r_epsi)
    r_hu_q         = MATMUL(r_Qelmt(2,:), r_epsi)
    r_hv_q         = MATMUL(r_Qelmt(3,:), r_epsi)
    r_gradb_q(:,1) = MATMUL(r_Selmt(1,:), r_ddxq)
    r_gradb_q(:,2) = MATMUL(r_Selmt(1,:), r_ddyq)
    r_laphu_q      = MATMUL(r_Qelmt(2,:), r_lapq)
    r_laphv_q      = MATMUL(r_Qelmt(3,:), r_lapq)

!--- quadrature loop
    elmt_quad_loop: DO i_quad=1, i_equadpts

!--- compute the source terms
      r_source(1) = 0.0_GRID_SR

      r_source(2) = &
        + r_iswet * r_h_q(i_quad) * r_gradb_q(i_quad,1)

      r_source(3) = &
        + r_iswet * r_h_q(i_quad) * r_gradb_q(i_quad,2)

!--- multiply with r_eqwei*Minv*psi for each dof
      elmt_dof_loop: DO i_dof=1,i_faceunknowns
        r_rhs(:,i_dof) = r_rhs(:,i_dof) - r_eqwei(i_quad)*r_eMinvpsi(i_dof,i_quad)*r_source
      END DO elmt_dof_loop

    END DO elmt_quad_loop

  END SUBROUTINE source

!*******************************************************************************
! DESCRIPTION of [FUNCTION divflux]:
!> @brief computes the divergence and source terms of the equation in consistent form
!
!> @param[in]       r_Q         state vector (fluid depth, momentum in x-, momentum in y-direction)
!> @param[in]       r_dQdx
!> @param[in]       r_dQdy
!> @param[in]       r_wet
!> @return                      return value
!
  FUNCTION divflux(r_Q, r_dQdx, r_dQdy, r_iswet) RESULT(r_Div)

    IMPLICIT NONE

    REAL (KIND = GRID_SR), DIMENSION(i_nprogvars), INTENT(in) :: r_Q, r_dQdx, r_dQdy
    REAL (KIND = GRID_SR), INTENT(in), OPTIONAL               :: r_iswet
    REAL (KIND = GRID_SR), DIMENSION(i_nprogvars)             :: r_Div

!--- local declarations
    REAL (KIND = GRID_SR)                                     :: r_u, r_v, r_wet
    REAL (KIND = GRID_SR), DIMENSION(i_nprogvars,2)           :: r_div_tmp

    IF (PRESENT(r_iswet)) THEN
      r_wet = r_iswet
    ELSE
      r_wet = 1.0_GRID_SR
    END IF

    r_u = velocity(r_Q(1), r_Q(2))
    r_v = velocity(r_Q(1), r_Q(3))

    r_div_tmp(1,1) = r_dQdx(2)
    r_div_tmp(1,2) = r_dQdy(3)

    r_div_tmp(2,1) = - r_u*r_u * r_dQdx(1) + 2.0_GRID_SR * r_dQdx(2) * r_u &
                     + r_wet*r_Q(1)*r_dQdx(1)
    r_div_tmp(2,2) = - r_u*r_v * r_dQdy(1) + r_dQdy(2)*r_v + r_dQdy(3)*r_u

    r_div_tmp(3,1) = - r_u*r_v * r_dQdx(1) + r_dQdx(2)*r_v + r_dQdx(3)*r_u
    r_div_tmp(3,2) = - r_v*r_v * r_dQdy(1) + 2.0_GRID_SR * r_dQdy(3) * r_v &
                     + r_wet*r_Q(1)*r_dQdy(1)

!--- compute full RHS
    r_Div(:) = r_div_tmp(:,1) + r_div_tmp(:,2)

  END FUNCTION divflux

!*******************************************************************************
! DESCRIPTION of [SUBROUTINE pre_timestepping]:
!> @brief pre-timestepping routine for, e.g., split implicit discretization
!>        of friction
!
!> @param         r_Q
!> @param         r_dt
!> @param         i_numelmt
!> @param         i_faceunknowns
!> @param         i_elementdofs
!
  SUBROUTINE pre_timestepping(r_Q, r_dt, i_numelmt, i_faceunknowns, i_elementdofs)

    IMPLICIT NONE

    REAL (KIND = GRID_SR), DIMENSION(:,:)                   :: r_Q
    REAL (KIND = GRID_SR)                                   :: r_dt
    INTEGER (KIND = GRID_SI)                                :: i_numelmt, i_faceunknowns
    INTEGER (KIND = GRID_SI), DIMENSION(:,:)                :: i_elementdofs

    INTEGER (KIND = GRID_SI)                                :: i_elmt, i_cnt
    REAL (KIND = GRID_SR)                                   :: r_n, r_mabs, r_fac, &
                                                               r_Sf, r_D

    r_fac = -GRID_GRAV**(4.0_GRID_SR/3.0_GRID_SR) * p_equationparam%r_gamma**2
    elmt_loop: DO i_elmt = 1,i_numelmt
      DO i_cnt = 1, i_faceunknowns
        r_mabs = SQRT(r_Q(2,i_elementdofs(i_cnt,i_elmt))**2 + r_Q(3,i_elementdofs(i_cnt,i_elmt))**2)
        IF ((r_Q(1,i_elementdofs(i_cnt,i_elmt)) > p_equationparam%r_wettol) .AND. &
            (r_mabs > 1.e-8)) THEN
          r_Sf = r_fac * r_Q(2,i_elementdofs(i_cnt,i_elmt)) * &
                 r_mabs / r_Q(1,i_elementdofs(i_cnt,i_elmt))**(7.0_GRID_SR/3.0_GRID_SR)
          r_D  = 1.0_GRID_SR - r_dt * r_fac * &
                 (2.0_GRID_SR*r_Q(2,i_elementdofs(i_cnt,i_elmt))**2 + r_Q(3,i_elementdofs(i_cnt,i_elmt))**2) / &
                 (r_Q(1,i_elementdofs(i_cnt,i_elmt))**(7.0_GRID_SR/3.0_GRID_SR) * r_mabs)
          IF (ABS(r_dt*r_Sf/r_D) < ABS(r_Q(2,i_elementdofs(i_cnt,i_elmt)))) THEN
            r_Q(2,i_elementdofs(i_cnt,i_elmt)) = r_Q(2,i_elementdofs(i_cnt,i_elmt)) + r_dt*r_Sf/r_D
          ELSE
            r_Q(2,i_elementdofs(i_cnt,i_elmt)) = 0.0_GRID_SR
          END IF

          r_Sf = r_fac * r_Q(3,i_elementdofs(i_cnt,i_elmt)) * &
                 r_mabs / r_Q(1,i_elementdofs(i_cnt,i_elmt))**(7.0_GRID_SR/3.0_GRID_SR)
          r_D  = 1.0_GRID_SR - r_dt * r_fac * &
                 (r_Q(2,i_elementdofs(i_cnt,i_elmt))**2 + 2.0_GRID_SR*r_Q(3,i_elementdofs(i_cnt,i_elmt))**2) / &
                 (r_Q(1,i_elementdofs(i_cnt,i_elmt))**(7.0_GRID_SR/3.0_GRID_SR) * r_mabs)

          IF (ABS(r_dt*r_Sf/r_D) < ABS(r_Q(3,i_elementdofs(i_cnt,i_elmt)))) THEN
            r_Q(3,i_elementdofs(i_cnt,i_elmt)) = r_Q(3,i_elementdofs(i_cnt,i_elmt)) + r_dt*r_Sf/r_D
          ELSE
            r_Q(3,i_elementdofs(i_cnt,i_elmt)) = 0.0_GRID_SR
          END IF
        END IF
      END DO
    END DO elmt_loop

  END SUBROUTINE pre_timestepping

!*******************************************************************************
! Only routines for IO follow
!*******************************************************************************
! DESCRIPTION of [SUBROUTINE io_initialize]:
!> @brief initializes general IO variable names and units for IO file output
!>
!
  SUBROUTINE io_initialize()

    IMPLICIT NONE

    VAR1D_PARACHECK = VAR1D_HZERO

!--- depth
    p_iovars(1)%r_factor        = 1.0_GRID_SR / GRID_GRAV
    p_iovars(1)%c_varname       = 'depth'
    p_iovars(1)%c_long_name     = 'fluid depth'
    p_iovars(1)%c_standard_name = 'fluid_depth' ! not CF conform!
    p_iovars(1)%c_units         = 'm'

!--- x-momentum
    p_iovars(2)%r_factor        = 1.0_GRID_SR / GRID_GRAV
    p_iovars(2)%c_varname       = 'm_x'
    p_iovars(2)%c_long_name     = 'momentum in x-direction'
    p_iovars(2)%c_standard_name = 'momentum_in_x-direction' ! not CF conform!
    p_iovars(2)%c_units         = 'm**2/s'

!--- y-momentum
    p_iovars(3)%r_factor        = 1.0_GRID_SR / GRID_GRAV
    p_iovars(3)%c_varname       = 'm_y'
    p_iovars(3)%c_long_name     = 'momentum in y-direction'
    p_iovars(3)%c_standard_name = 'momentum_in_y-direction' ! not CF conform!
    p_iovars(3)%c_units         = 'm**2/s'

!--- bathymetry
    p_iovars(4)%r_factor        = 1.0_GRID_SR / GRID_GRAV
    p_iovars(4)%c_varname       = 'bathy'
    p_iovars(4)%c_long_name     = 'bathymetry'
    p_iovars(4)%c_standard_name = 'bathymetry' ! not CF conform!
    p_iovars(4)%c_units         = 'm'

!--- exact depth
    p_iovars(5)%r_factor        = 1.0_GRID_SR
    p_iovars(5)%c_varname       = 'exact_depth'
    p_iovars(5)%c_long_name     = 'exact fluid depth'
    p_iovars(5)%c_standard_name = 'exact_fluid depth' ! not CF conform!
    p_iovars(5)%c_units         = 'm'

!--- exact x-momentum
    p_iovars(6)%r_factor        = 1.0_GRID_SR
    p_iovars(6)%c_varname       = 'exact_m_x'
    p_iovars(6)%c_long_name     = 'exact momentum in x-direction'
    p_iovars(6)%c_standard_name = 'exact_momentum_in_x-direction' ! not CF conform!
    p_iovars(6)%c_units         = 'm**2/s'

!--- exact y-momentum
    p_iovars(7)%r_factor        = 1.0_GRID_SR
    p_iovars(7)%c_varname       = 'exact_m_y'
    p_iovars(7)%c_long_name     = 'exact momentum in y-direction'
    p_iovars(7)%c_standard_name = 'exact_momentum_in_y-direction' ! not CF conform!
    p_iovars(7)%c_units         = 'm**2/s'

  END SUBROUTINE io_initialize

!*******************************************************************************
END MODULE DG_equation
