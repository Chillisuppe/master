!*******************************************************************************
!
!> @file  DG_equation_TE_notWB.F90
!> @brief contains module DG_equation
!
!*******************************************************************************
! MODULE DESCRIPTION:
!> @brief Defines equation specific parts of the scheme, such as prognostic
!>        variables, flux function, source terms, equation specific parameters
!>        This equation set uses total energy (w/o gravitational energy) as energy representation
!>        classical non-balanced discretization of gravity source term
!>        In this case, the source term yields instabilities since the regimes
!>        the momentum does not balance out. The additional momentum introduced
!>        by the source term \f$ S(q) = - \rho \cdot GRAV \f$ deviates from the pressure
!>        gradient.
!
MODULE DG_equation

  USE GRID_api
  USE FLASH_parameters
  USE IO_equation, ONLY : p_equationparam

  PRIVATE
  PUBLIC  i_nprogvars, i_nsrcterms, i_valQ, i_valS, &
          VAR1D_RHO, VAR2D_MOM, VAR1D_NRG, VAR1D_PARACHECK, &
          velocity, equation_initialize, divflux, pressure, &
          flux, source, pre_timestepping, PrimitiveToProgvars, &
          io_initialize, p_iovars

  INTEGER (KIND = GRID_SI), PARAMETER               :: i_nprogvars = 4 !< number of unknowns
  INTEGER (KIND = GRID_SI), PARAMETER               :: i_nsrcterms = 0 !< number of grid source term variables

  INTEGER (KIND = GRID_SI), DIMENSION(i_nprogvars)  :: i_valQ   !< array of index handles for registered unknowns
  INTEGER (KIND = GRID_SI), DIMENSION(i_nsrcterms)  :: i_valS   !< array of index handles for registered source term variables

!--- constants for the DG variables
  INTEGER (KIND = GRID_SI)                          :: VAR1D_RHO       !< density variable
  INTEGER (KIND = GRID_SI)                          :: VAR2D_MOM       !< momentum variable
  INTEGER (KIND = GRID_SI)                          :: VAR1D_NRG       !< potential energy variable
  INTEGER (KIND = GRID_SI)                          :: VAR1D_PARACHECK !< variable for consistency check in io_paraplot

!--- array for generic IO
  TYPE (io_vars), DIMENSION(2 * i_nprogvars + i_nsrcterms)  :: p_iovars

  INTERFACE velocity
    MODULE PROCEDURE velocity_scalar, velocity_array
  END INTERFACE

  CONTAINS

!*******************************************************************************
  FUNCTION velocity_scalar(r_rho, r_rhou) RESULT(r_u)

    IMPLICIT NONE

    REAL (KIND = GRID_SR), INTENT(in)               :: r_rho, r_rhou
    REAL (KIND = GRID_SR)                           :: r_u

    r_u = r_rhou/r_rho

  END FUNCTION velocity_scalar

!*******************************************************************************
  FUNCTION velocity_array(r_rho, r_rhou) RESULT(r_u)

    IMPLICIT NONE

    REAL (KIND = GRID_SR), DIMENSION(:), INTENT(in) :: r_rho, r_rhou
    REAL (KIND = GRID_SR), DIMENSION(SIZE(r_rho))   :: r_u

    r_u = r_rhou/r_rho

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

!--- total variables
    VAR1D_RHO = grid_registerfemvar(i_femtype)
    VAR2D_MOM = grid_registerfemvar(i_femtype, i_length=2)
    VAR1D_NRG = grid_registerfemvar(i_femtype)

    i_valQ = (/ VAR1D_RHO, VAR2D_MOM, VAR2D_MOM+1, VAR1D_NRG /)
!     i_valS = (/ /)

  END SUBROUTINE equation_initialize

!*******************************************************************************
! DESCRIPTION of [FUNCTION flux]:
!> @brief Routine for flux computation. Has to be defined by the user.
!>
!> @param[in]       r_Q         state vector (rho, rho u, rho v, rho e)
!> @param[in]       r_S         source terms
!> @return                      return value
!
  FUNCTION flux(r_Q, r_S) RESULT(r_flux)

    IMPLICIT NONE

    REAL (KIND = GRID_SR), DIMENSION(i_nprogvars), INTENT(in) :: r_Q
    REAL (KIND = GRID_SR), DIMENSION(i_nsrcterms), INTENT(in), OPTIONAL :: r_S
    REAL (KIND = GRID_SR), DIMENSION(i_nprogvars,2)           :: r_flux

    REAL (KIND = GRID_SR)                                     :: r_u, r_v, r_p

    r_u = velocity(r_Q(1), r_Q(2))
    r_v = velocity(r_Q(1), r_Q(3))
    r_p = pressure(r_Q)

    r_flux(:,:)  = 0.0_GRID_SR

    r_flux(1,1)  = r_Q(2)
    r_flux(1,2)  = r_Q(3)

    r_flux(2,1)  = r_Q(2)*r_u + r_p
    r_flux(2,2)  = r_Q(2)*r_v

    r_flux(3,1)  = r_Q(3)*r_u
    r_flux(3,2)  = r_Q(3)*r_v + r_p

    r_flux(4,1) = (r_Q(4) + r_p) * r_u
    r_flux(4,2) = (r_Q(4) + r_p) * r_v

  END FUNCTION flux

!*******************************************************************************
! DESCRIPTION of [SUBROUTINE source]:
!> @brief
!
!> @param         r_rhs
!> @param         i_faceunknowns
!> @param         i_equadpts
!> @param         r_eqwei
!> @param         r_eMinvpsi
!> @param         r_Qelmt       state vector (rho, rho u, rho v, rho theta)
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
    REAL (KIND = GRID_SR), DIMENSION(i_equadpts)            :: r_rho_q, r_rhou_q, r_rhov_q, &
     r_u_q, r_v_q, r_nrg_q, r_dpdy_q
    REAL (KIND = GRID_SR), DIMENSION(i_nprogvars)           :: r_source

!--- interpolate to quadrature points and compute derivatives
    r_rho_q  = MATMUL(r_Qelmt(1,:), r_epsi)

!--- calculate velocity
    r_v_q = velocity(r_rho_q, r_rhov_q)

!--- quadrature loop
    elmt_quad_loop: DO i_quad=1, i_equadpts

!--- compute the source terms
      r_source(:) = 0.0_GRID_SR
      r_source(3) = r_rho_q(i_quad) * GRID_GRAV * p_equationparam%r_gravityswitch
      r_source(4) = r_rho_q(i_quad) * GRID_GRAV * r_v_q(i_quad) * p_equationparam%r_gravityswitch

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
!> @param[in]       r_Q         state vector (rho, rho u, rho v, rho e)
!> @param[in]       r_dQdx
!> @param[in]       r_dQdy
!> @return                      return value
!
  FUNCTION divflux(r_Q, r_dQdx, r_dQdy) RESULT(r_Div)

    IMPLICIT NONE

    REAL (KIND = GRID_SR), DIMENSION(i_nprogvars), INTENT(IN) :: r_Q, r_dQdx, r_dQdy
    REAL (KIND = GRID_SR), DIMENSION(i_nprogvars)             :: r_Div

!--- local declarations
    REAL (KIND = GRID_SR)                                     :: r_u, r_v, r_p
    REAL (KIND = GRID_SR), DIMENSION(2)                       :: r_gradp
    REAL (KIND = GRID_SR), DIMENSION(i_nprogvars,2)           :: r_div_tmp

    r_u = velocity(r_Q(1), r_Q(2))
    r_v = velocity(r_Q(1), r_Q(3))

    r_p = pressure(r_Q)
    r_gradp = gradpressure(r_Q, r_dQdx, r_dQdy)

    r_div_tmp(:,:)  = 0.0_GRID_SR

    r_div_tmp(1,1)  = r_dQdx(2)
    r_div_tmp(1,2)  = r_dQdy(3)

    r_div_tmp(2,1)  = - r_u*r_u * r_dQdx(1) + 2.0 * r_dQdx(2) * r_u + r_gradp(1)
    r_div_tmp(2,2)  = - r_u*r_v * r_dQdy(1) + r_dQdy(2)*r_v + r_dQdy(3)*r_u

    r_div_tmp(3,1)  = - r_u*r_v * r_dQdx(1) + r_dQdx(2)*r_v + r_dQdx(3)*r_u
    r_div_tmp(3,2)  = - r_v*r_v * r_dQdy(1) + 2.0 * r_dQdy(3) * r_v + r_gradp(2)

    r_div_tmp(4,1) = r_u*(r_dQdx(4)+r_gradp(1)) + (r_Q(4)+r_p)/r_Q(1)*(r_dQdx(2)-r_u*r_dQdx(1))
    r_div_tmp(4,2) = r_v*(r_dQdy(4)+r_gradp(2)) + (r_Q(4)+r_p)/r_Q(1)*(r_dQdy(3)-r_v*r_dQdy(1))

!--- compute full RHS
    r_Div(:) = r_div_tmp(:,1) + r_div_tmp(:,2)

  END FUNCTION divflux

!*******************************************************************************
! DESCRIPTION of [FUNCTION pressure]:
!> @brief computes the pressure
!
!> @param[in]       r_Q         state vector (rho, rho u, rho v, rho e)
!> @return                      return value
!
  FUNCTION pressure(r_Q) RESULT(r_pressure)

    IMPLICIT NONE

    REAL (KIND = GRID_SR), DIMENSION(i_nprogvars), INTENT(IN) :: r_Q
    REAL (KIND = GRID_SR)                                     :: r_pressure, r_u, r_v

    r_u = velocity(r_Q(1), r_Q(2))
    r_v = velocity(r_Q(1), r_Q(3))

    r_pressure = (p_equationparam%r_gamma - 1.0_GRID_SR) * &
                 (r_Q(4) - 0.5_GRID_SR * (r_u**2+r_v**2)*r_Q(1))

  END FUNCTION pressure

!*******************************************************************************
! DESCRIPTION of [FUNCTION gradpressure]:
!> @brief computes the pressure gradient
!
!> @param[in]       r_Q             state vector \f$(\rho, \rho u, \rho v, \rho e)\f$
!> @param[in]       r_dQdx
!> @param[in]       r_dQdy
!> @return                          return value
!
  FUNCTION gradpressure(r_Q, r_dQdx, r_dQdy) RESULT(r_gradpressure)

    IMPLICIT NONE

    REAL (KIND = GRID_SR), DIMENSION(i_nprogvars), INTENT(IN) :: r_Q, r_dQdx, r_dQdy
    REAL (KIND = GRID_SR), DIMENSION(2)                       :: r_gradpressure

!--- local declarations
    REAL (KIND = GRID_SR) :: r_u, r_v, r_velaux, r_coeff

    r_u = velocity(r_Q(1), r_Q(2))
    r_v = velocity(r_Q(1), r_Q(3))

    r_velaux = 0.5_GRID_SR * (r_u**2 + r_v**2)
    r_coeff  = p_equationparam%r_gamma - 1.0_GRID_SR

    r_gradpressure(1) = r_coeff * (r_dQdx(4) - r_u * r_dQdx(2) - r_v * r_dQdx(3) + r_velaux * r_dQdx(1))
    r_gradpressure(2) = r_coeff * (r_dQdy(4) - r_u * r_dQdy(2) - r_v * r_dQdy(3) + r_velaux * r_dQdy(1))

  END FUNCTION gradpressure

!*******************************************************************************
! DESCRIPTION of [SUBROUTINE PrimitiveToProgvars]:
!> @brief convertes from primite to prognostic variables
!>
!> @param[in]       r_Qprim     vector of primitive variables (rho, u, v, pressure)
!> @param[out]      r_Q         state vector of primary variables
!
  SUBROUTINE PrimitiveToProgvars(r_Qprim, r_Q)

    IMPLICIT NONE

    REAL (KIND = GRID_SR), DIMENSION(:), INTENT(IN)  :: r_Qprim
    REAL (KIND = GRID_SR), DIMENSION(:), INTENT(OUT) :: r_Q

!--- calculate prognostic variables
!--- the total energy formula is rho*e = rho*c_v*exner*pot_temp + rho/2*u*u (see Chandrashekar & Zenk 2016)
    r_Q(1) = r_Qprim(1)
    r_Q(2) = r_Qprim(1) * r_Qprim(2)
    r_Q(3) = r_Qprim(1) * r_Qprim(3)
    r_Q(4) = r_Qprim(4) / (p_equationparam%r_gamma-1.0_GRID_SR) + &
             0.5_GRID_SR * DOT_PRODUCT(r_Qprim(2:3), r_Qprim(2:3)) * r_Q(1)

  END SUBROUTINE PrimitiveToProgvars

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

    VAR1D_PARACHECK = VAR1D_RHO

!--- density
    p_iovars(1)%r_factor        = 1.0_GRID_SR
    p_iovars(1)%c_varname       = 'density'
    p_iovars(1)%c_long_name     = 'fluid density'
    p_iovars(1)%c_standard_name = 'fluid_density' ! not CF conform!
    p_iovars(1)%c_units         = 'kg/m^3'

!--- x-momentum
    p_iovars(2)%r_factor        = 1.0_GRID_SR
    p_iovars(2)%c_varname       = 'm_x'
    p_iovars(2)%c_long_name     = 'momentum in x-direction'
    p_iovars(2)%c_standard_name = 'momentum_in_x-direction' ! not CF conform!
    p_iovars(2)%c_units         = 'm**2/s'

!--- y-momentum
    p_iovars(3)%r_factor        = 1.0_GRID_SR
    p_iovars(3)%c_varname       = 'm_y'
    p_iovars(3)%c_long_name     = 'momentum in y-direction'
    p_iovars(3)%c_standard_name = 'momentum_in_y-direction' ! not CF conform!
    p_iovars(3)%c_units         = 'm**2/s'

!--- energy
    p_iovars(4)%r_factor        = 1.0_GRID_SR
    p_iovars(4)%c_varname       = 'energy'
    p_iovars(4)%c_long_name     = 'energy'
    p_iovars(4)%c_standard_name = 'energy' ! not CF conform!
    p_iovars(4)%c_units         = 'J'

!--- exact depth
    p_iovars(5)%r_factor        = 1.0_GRID_SR
    p_iovars(5)%c_varname       = 'exact_density'
    p_iovars(5)%c_long_name     = 'exact fluid density'
    p_iovars(5)%c_standard_name = 'exact_fluid_density' ! not CF conform!
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

!--- exact energy
    p_iovars(8)%r_factor        = 1.0_GRID_SR
    p_iovars(8)%c_varname       = 'exact_energy'
    p_iovars(8)%c_long_name     = 'exact energy'
    p_iovars(8)%c_standard_name = 'exact_energy' ! not CF conform!
    p_iovars(8)%c_units         = 'J'


  END SUBROUTINE io_initialize

!*******************************************************************************
END MODULE DG_equation
