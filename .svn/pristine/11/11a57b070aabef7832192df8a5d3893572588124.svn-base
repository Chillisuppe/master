!*******************************************************************************
! DESCRIPTION:
!> @brief Correction of the flux term according to a paper by Audusse et al.
!>
!> This will only take effect for discontinuous bathymetry and will prevent
!> the sea surface height from artificially becoming negative.
!
!> @param[in,out]
!> @param[in]     p_param
!>
!*******************************************************************************
MODULE DG_fluxcorrection

  USE GRID_api
  USE FLASH_parameters
  USE DG_utils
  USE DG_riemann_solver

  PRIVATE
  PUBLIC  :: flux_correction

 CONTAINS

!*******************************************************************************
  SUBROUTINE flux_correction(r_ssh_star_e, r_ssh_star_k, &
                             r_bathy_e, r_bathy_k, r_u_e, r_u_k, r_v_e, r_v_k, &
                             r_hu_e, r_hu_k, r_hv_e, r_hv_k, &
                             r_height_e, r_height_k, r_ssh_e, r_ssh_k, &
                             r_normal, r_Fstrong_e, r_Fstrong_k)

    IMPLICIT NONE

    REAL(KIND=GRID_SR), INTENT(in)                  :: r_ssh_star_e, r_ssh_star_k, r_bathy_e, r_bathy_k, &
                                                       r_u_e, r_u_k, r_v_e, r_v_k, &
                                                       r_hu_e, r_hu_k, r_hv_e, r_hv_k, &
                                                       r_height_e, r_height_k, r_ssh_e, r_ssh_k
    REAL(KIND=GRID_SR)                              :: r_height_star_e, r_height_star_k, &
                                                       r_hu_star_e, r_hu_star_k, r_hv_star_e, r_hv_star_k
    REAL(KIND=GRID_SR), DIMENSION(:), INTENT(in)    :: r_normal
    REAL(KIND=GRID_SR), DIMENSION(:), INTENT(inout) :: r_Fstrong_e, r_Fstrong_k
    REAL (KIND = GRID_SR), DIMENSION(3)             :: r_flux_e, r_Fstar
    REAL (KIND = GRID_SR), DIMENSION(3,2)           :: r_F_e, r_flux_corr_e, r_flux_corr_k

! Evaluate numerical flux
    r_F_e    = flux(r_ssh_star_e, r_bathy_e, r_hu_e, r_hv_e)
    r_flux_e = matmul(r_F_e, r_normal)

! Compute star quantities
    r_height_star_e = r_ssh_star_e+r_bathy_e
    r_height_star_k = r_ssh_star_k+r_bathy_k

    r_hu_star_e = r_height_star_e*r_u_e
    r_hu_star_k = r_height_star_k*r_u_k
    r_hv_star_e = r_height_star_e*r_v_e
    r_hv_star_k = r_height_star_k*r_v_k

! Compute star flux
    r_Fstar = riemannsolver(r_ssh_star_e, r_ssh_star_k, r_bathy_e, r_bathy_k, &
                            r_u_e, r_u_k, r_v_e, r_v_k, r_hu_star_e, r_hu_star_k, &
                            r_hv_star_e, r_hv_star_k, r_normal)

! Using the correction as in Xing and Shu, 2010  !! speed correct for correction
! Initalize correction terms
    r_flux_corr_e    = 0._GRID_SR
    r_flux_corr_k    = 0._GRID_SR

! Compute the correction terms
    r_flux_corr_e(2,1) = r_height_e**2 - r_height_star_e**2
    r_flux_corr_e(3,2) = r_flux_corr_e(2,1)
    r_flux_corr_k(2,1) = r_height_k**2 - r_height_star_k**2
    r_flux_corr_k(3,2) = r_flux_corr_k(2,1)

! Rotate to normal components
    r_flux_corr_e(2,1) = 0.5*r_normal(1)*r_flux_corr_e(2,1)
    r_flux_corr_e(3,2) = 0.5*r_normal(2)*r_flux_corr_e(3,2)
    r_flux_corr_k(2,1) = 0.5*r_normal(1)*r_flux_corr_k(2,1)
    r_flux_corr_k(3,2) = 0.5*r_normal(2)*r_flux_corr_k(3,2)

    r_Fstrong_e  = r_Fstar + matmul(r_flux_corr_e,r_normal)
    r_Fstrong_k  = r_Fstar + matmul(r_flux_corr_k,r_normal)

 END SUBROUTINE flux_correction

!*******************************************************************************
END MODULE DG_fluxcorrection
