!*******************************************************************************
!
!> @file  MISC_bathy_etopo.f90
!> @brief contains module MISC_bathy
!
!*******************************************************************************
!
! VERSION(S):
!  1. original version                j. behrens    01/2005
!  2. this is for TsunaFlash          j. behrens    07/2006
!  3. homogenized with GRD bathy      j. behrens    12/2011
!
!*******************************************************************************
! MODULE DESCRIPTION:
!> @brief reads bathymetry data from ETOPO data set
!
MODULE MISC_bathy

  USE GRID_api

!--- global constants
  LOGICAL                                             :: l_bathyinit = .FALSE.
  INTEGER (KIND = GRID_SI), PARAMETER                 :: i_etop = 4321_GRID_SI
  INTEGER (KIND = GRID_SI), PARAMETER                 :: j_etop = 2161_GRID_SI
  REAL (KIND = GRID_SR), DIMENSION(i_etop,j_etop)     :: r_etopo5

  PRIVATE
  PUBLIC  :: bathy_initialize, bathy, bathy_inquire

  CONTAINS

!*******************************************************************************
! DESCRIPTION of [SUBROUTINE bathy_initialize]:
!> @brief initializes bathymetry data.
!>
!> @param[in]     c_filename       uplift data input file
!>
!
  SUBROUTINE bathy_initialize(c_filename)

    IMPLICIT NONE

    CHARACTER (len=*), INTENT(in)                         :: c_filename

!--- local declarations
    INTEGER*2, DIMENSION(:,:), ALLOCATABLE                :: i_rawetopo5
    INTEGER (KIND = GRID_SI)                              :: i_alct, i_unit, i_cnt, j_cnt, i_fst
    CHARACTER (len = 80)                                  :: c_dummy

!--- allocate raw data array
    ALLOCATE(i_rawetopo5(i_etop,j_etop), stat=i_alct)
    IF(i_alct /= 0) &
      CALL grid_error(c_error='[bathy_initialize]: could not allocate raw data field')

!--- open topography file
    i_unit = 21
    open(i_unit, file=c_filename, form='formatted', iostat=i_fst)
    IF(i_fst /= 0) &
      CALL grid_error(c_error='[bathy_initialize]: could open file topography file')
    IF(GRID_parameters%iolog > 0) &
      write(GRID_parameters%iolog,*) 'INFO: Filename: ', c_filename, ' opened on unit: ', i_unit

!--- read raw data
    DO j_cnt=1,8
      read(i_unit,*) c_dummy
    END DO
    DO j_cnt=1,j_etop
      read(i_unit,'(4321I6)') (i_rawetopo5(i_cnt,j_cnt), i_cnt=1,i_etop)
    END DO

!--- close topography file
    close(i_unit)
    IF(GRID_parameters%iolog > 0) &
      write(GRID_parameters%iolog,*) 'INFO: Closed file on unit: ', i_unit

!--- convert data
    r_etopo5 = REAL(i_rawetopo5, GRID_SR)

!--- deallocate raw array
    DEALLOCATE(i_rawetopo5)

    l_bathyinit = .TRUE.

  END SUBROUTINE bathy_initialize

!*******************************************************************************
! DESCRIPTION of [FUNCTION bathy]:
!> @brief calculates a tracer/constituent initial function
!>
!> @param         r_coos      coordinate array
!> @return                    value at postition r_coo
!
  FUNCTION bathy(r_coos) RESULT (r_val)

    IMPLICIT NONE

    REAL (KIND = GRID_SR), DIMENSION(GRID_dimension), INTENT(in)  :: r_coos
    REAL (KIND = GRID_SR)                                         :: r_val

!--- local declarations
    INTEGER (KIND = GRID_SI)                                      :: i_pos, j_pos, i_tmp, j_tmp

!--- check if bathymetry data is loaded
    IF(.NOT. l_bathyinit) THEN
      CALL grid_error(c_error='[bathy]: bathymetry must be initialized before!')
    END IF

!--- we assume spherical coordinate and derive the nearest
!           integer position from that

!     i_tmp = INT(((r_coos(1)+GRID_PI)/(2._GRID_SR* GRID_PI)) * &
!                 360._GRID_SR * 12._GRID_SR, GRID_SI) + 1_GRID_SI
!     j_tmp = INT(((-r_coos(2)+0.5_GRID_SR*GRID_PI)/(GRID_PI)) * &
!                 180._GRID_SR * 12._GRID_SR, GRID_SI) + 1_GRID_SI
    i_tmp = INT((r_coos(1) * 12._GRID_SR), GRID_SI) + 5_GRID_SI
    j_tmp = INT(((-r_coos(2) + 90._GRID_SR) * 12._GRID_SR), GRID_SI) - 8_GRID_SI
    j_pos = min(j_tmp, j_etop)
    i_pos = min(i_tmp, i_etop)

!--- now the value
    r_val = abs(r_etopo5(i_pos,j_pos))
    IF(r_val < 10.0_GRID_SR) r_val = 10.0_GRID_SR

  END FUNCTION bathy

!*******************************************************************************
! DESCRIPTION of [SUBROUTINE bathy_inquire]:
!> @brief inquires domain size from given bathymetry data file
!
!> @param         r_east      eastern boundary
!> @param         r_west      western boundary
!> @param         r_north     northern boundary
!> @param         r_south     southern boundary
!> @param         i_xsize     number of nodes in x
!> @param         i_ysize     number of nodes in y
!> @param         r_xres      resolution in x
!> @param         r_yres      resolution in y
!>
!> @note  all output parameters are optional
!
  SUBROUTINE bathy_inquire(r_east, r_west, r_north, r_south, i_xsize, i_ysize, r_xres, r_yres)

    IMPLICIT NONE

    REAL (KIND = GRID_SR), OPTIONAL     :: r_east
    REAL (KIND = GRID_SR), OPTIONAL     :: r_west
    REAL (KIND = GRID_SR), OPTIONAL     :: r_north
    REAL (KIND = GRID_SR), OPTIONAL     :: r_south
    REAL (KIND = GRID_SR), OPTIONAL     :: r_xres
    REAL (KIND = GRID_SR), OPTIONAL     :: r_yres
    INTEGER (KIND = GRID_SI), OPTIONAL  :: i_xsize
    INTEGER (KIND = GRID_SI), OPTIONAL  :: i_ysize

!--- local declarations

!--- check if bathymetry data is loaded
    IF(.NOT. l_bathyinit) THEN
      CALL grid_error(c_error='[bathy_inquire]: bathymetry must be initialized before!')
    END IF

!--- fill optional output parameters

    IF(present(r_east))  r_east  = 0.0_GRID_SR
    IF(present(r_west))  r_west  = 360.0_GRID_SR
    IF(present(r_north)) r_north = 180.0_GRID_SR
    IF(present(r_south)) r_south = 0.0_GRID_SR
    IF(present(r_xres))  r_xres  = 360.0_GRID_SR / (i_etop-1)
    IF(present(r_yres))  r_yres  = 180.0_GRID_SR / (j_etop-1)
    IF(present(i_xsize)) i_xsize = i_etop
    IF(present(i_ysize)) i_ysize = j_etop

  END SUBROUTINE bathy_inquire

!*******************************************************************************
END MODULE MISC_bathy
