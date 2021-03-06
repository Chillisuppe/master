!*******************************************************************************
!
!> @file  MISC_bathy_ncfile.f90
!> @brief contains module MISC_bathy
!
!*******************************************************************************
!
! VERSION(S):
!  1. original version                j. behrens    01/2005
!  2. this is for TsunaFlash          j. behrens    07/2006
!  3. include NetCDF file read        j. behrens    12/2011
!
!*******************************************************************************
! MODULE DESCRIPTION:
!> @brief read bathymetry data from NetCDF file. This can be:
!>         - a GRD (NetCDF) file as being retrieved from
!>           http://www.marine-geo.org/portals/gmrt/
!>         - a NetCDF file in COARD convention (note: such a file must have a
!>           global attribute "Conventions" with value "COARDS")
!
MODULE MISC_bathy

  USE netcdf
  USE GRID_api

!--- global constants
  LOGICAL                                             :: l_bathyinit = .FALSE.
  INTEGER (KIND = GRID_SI), SAVE                      :: i_nlonpts, i_nlatpts
  REAL (KIND = GRID_SR), DIMENSION(:), ALLOCATABLE    :: r_lon, r_lat
  REAL (KIND = GRID_SR), DIMENSION(:,:), ALLOCATABLE  :: r_topo

  PRIVATE
  PUBLIC  :: bathy_initialize, bathy, bathy_inquire, get_shorepoints

  CONTAINS
!*******************************************************************************
! DESCRIPTION of [SUBROUTINE bathy_initialize]:
!> @brief initializes bathymetry data.
!>        Reads bathymetry/topography data from a NetCDF formatted GRD or COARDS file
!>
!> @param[in]     c_filename       uplift data input file
!> @param[in]     i_bathysmooth    no. of iteration for next neighbor smoothing
!>
!
  SUBROUTINE bathy_initialize(c_filename, i_bathysmooth)

    IMPLICIT NONE

!    INCLUDE "netcdf.inc"

    CHARACTER (len=*), INTENT(in)                         :: c_filename
    INTEGER (KIND = GRID_SI), INTENT(in), OPTIONAL        :: i_bathysmooth

!--- local declarations
    INTEGER (KIND = GRID_SI)                              :: i_alct, i_cnt, j_cnt, &
                                                             k_cnt, i_case, i_bsmooth
    INTEGER (KIND = GRID_SI)                              :: i_fileid, i_dimid, i_varid
    INTEGER (KIND = GRID_SI)                              :: i_ncstat         ! NetCDF library return status
    INTEGER (KIND = GRID_SI)                              :: i_lenconventions
    CHARACTER (LEN=:), ALLOCATABLE                        :: c_conventions

!--- grd specific
    INTEGER (KIND = GRID_SI)                              :: i_bathydim, i_bathysize1D    ! NetCDF dimensions
    INTEGER (KIND = GRID_SI)                              :: i_bathyX, i_bathyY, i_bathyZ ! NetCDF variable IDs
    INTEGER (KIND = GRID_SI)                              :: i_bathydim2D     ! NetCDF variable IDs
    INTEGER (KIND = GRID_SI), DIMENSION(GRID_dimension)   :: i_bathydim2Dvals ! NetCDF variable
    REAL (KIND = GRID_SR), DIMENSION(:), ALLOCATABLE      :: r_rawbathy
    REAL (KIND = GRID_SR)                                 :: r_dlon, r_dlat
    REAL (KIND = GRID_SR), DIMENSION(GRID_dimension)      :: r_bathyXrange, r_bathyYrange

    IF(present(i_bathysmooth)) THEN
      i_bsmooth = i_bathysmooth
    ELSE
      i_bsmooth = 0
    END IF

!--- open bathymetry file
    i_ncstat = nf90_open(c_filename, NF90_NOWRITE, i_fileid)
    IF(i_ncstat /= NF90_NOERR) &
      CALL grid_error(c_error='[bathy_initialize]: could not open bathy data file')
    IF(GRID_parameters%iolog > 0) &
      write(GRID_parameters%iolog,*) 'INFO: Opened bathymetry data file ', c_filename

!--- determine file format
    i_ncstat = nf90_inquire_attribute(i_fileid, nf90_global, "Conventions", len=i_lenconventions)
    IF(i_ncstat /= NF90_NOERR) THEN
      i_case = 0
    ELSE
      ALLOCATE(character(len=i_lenconventions) :: c_conventions, stat=i_alct)
      IF(i_alct /= 0) &
      CALL grid_error(c_error='[bathy_initialize]: could not allocate data field')

      i_ncstat = nf90_get_att(i_fileid, nf90_global, "Conventions", c_conventions)
      IF (c_conventions(1:6)=='COARDS') THEN
        i_case = 1
      ELSE
        i_case = 0
      END iF
    END IF

!--- just read input for bathy_inquire
    SELECT CASE (i_case)
!--- read grd data
    CASE(0)
!--- read x and y range of data
      i_ncstat = nf90_inq_varid(i_fileid, 'x_range', i_bathyX)
      IF(i_ncstat /= NF90_NOERR) &
        CALL grid_error(c_error='[bathy_initialize]: could not determine x range varid')
      i_ncstat = nf90_inq_varid(i_fileid, 'y_range', i_bathyY)
      IF(i_ncstat /= NF90_NOERR) &
        CALL grid_error(c_error='[bathy_initialize]: could not determine y range varid')
      i_ncstat = nf90_get_var(i_fileid, i_bathyX, r_bathyXrange)
      IF(i_ncstat /= NF90_NOERR) &
        CALL grid_error(c_error='[bathy_initialize]: could not read bathy x range data')
      i_ncstat = nf90_get_var(i_fileid, i_bathyY, r_bathyYrange)
      IF(i_ncstat /= NF90_NOERR) &
        CALL grid_error(c_error='[bathy_initialize]: could not read bathy y range data')

!--- read 2D dimensions
      i_ncstat = nf90_inq_varid(i_fileid, 'dimension', i_bathydim2D)
      IF(i_ncstat /= NF90_NOERR) &
        CALL grid_error(c_error='[bathy_initialize]: could not determine dimension varid')
      i_ncstat = nf90_get_var(i_fileid, i_bathydim2D, i_bathydim2Dvals)
      IF(i_ncstat /= NF90_NOERR) &
        CALL grid_error(c_error='[bathy_initialize]: could not read 2D dimension data')

      i_nlonpts = i_bathydim2Dvals(1)
      i_nlatpts = i_bathydim2Dvals(2)

      !--- allocate lat/lon arrays
      ALLOCATE(r_lon(i_nlonpts), r_lat(i_nlatpts), stat=i_alct)
      IF(i_alct /= 0) &
        CALL grid_error(c_error='[bathy_initialize]: could not allocate data field')

      r_dlon = (r_bathyXrange(2)-r_bathyXrange(1)) / (i_nlonpts-1)
      r_dlat = (r_bathyYrange(2)-r_bathyYrange(1)) / (i_nlatpts-1)

!--- calculate lat and lon
      DO i_cnt = 1, i_nlonpts
        r_lon(i_cnt) = r_bathyXrange(1) + (i_cnt-1) * r_dlon
      END DO
      DO j_cnt = 1, i_nlatpts
        r_lat(j_cnt) = r_bathyYrange(1) + (j_cnt-1) * r_dlat
      END DO

!--- read coards data
    CASE(1)
!--- determine array sizes.
      i_ncstat = nf90_inq_dimid(i_fileid, 'lon', i_dimid)
      IF(i_ncstat /= NF90_NOERR) &
        CALL grid_error(c_error='[bathy_initialize]: could not identify x dimension')
      i_ncstat = nf90_inquire_dimension(i_fileid, i_dimid, len=i_nlonpts)
      IF(i_ncstat /= NF90_NOERR) &
        CALL grid_error(c_error='[bathy_initialize]: could not read x dimension')
      i_ncstat = nf90_inq_dimid(i_fileid, 'lat', i_dimid)
      IF(i_ncstat /= NF90_NOERR) &
        CALL grid_error(c_error='[bathy_initialize]: could not identify y dimension')
      i_ncstat = nf90_inquire_dimension(i_fileid, i_dimid, len=i_nlatpts)
      IF(i_ncstat /= NF90_NOERR) &
        CALL grid_error(c_error='[bathy_initialize]: could not read y dimension')

!--- allocate data array
      ALLOCATE(r_lon(i_nlonpts), r_lat(i_nlatpts), stat=i_alct)
      IF(i_alct /= 0) &
        CALL grid_error(c_error='[bathy_initialize]: could not allocate data field')

!--- read data arrays
      i_ncstat = nf90_inq_varid(i_fileid, 'lon', i_varid)
      IF(i_ncstat /= NF90_NOERR) &
        CALL grid_error(c_error='[bathy_initialize]: could not determine x varid')
      i_ncstat = nf90_get_var(i_fileid, i_varid, r_lon)
      IF(i_ncstat /= NF90_NOERR) &
        CALL grid_error(c_error='[bathy_initialize]: could not read x vector')

      i_ncstat = nf90_inq_varid(i_fileid, 'lat', i_varid)
      IF(i_ncstat /= NF90_NOERR) &
        CALL grid_error(c_error='[bathy_initialize]: could not determine y varid')
      i_ncstat = nf90_get_var(i_fileid, i_varid, r_lat)
      IF(i_ncstat /= NF90_NOERR) &
        CALL grid_error(c_error='[bathy_initialize]: could not read y vector')

    END SELECT

!--- allocate data array
    ALLOCATE(r_topo(i_nlonpts, i_nlatpts), stat=i_alct)
    IF(i_alct /= 0) &
      CALL grid_error(c_error='[bathy_initialize]: could not allocate raw data field')

    SELECT CASE(i_case)
!--- read grd data
    CASE(0)
!--- determine array sizes.
      i_ncstat = nf90_inq_dimid(i_fileid, 'side', i_dimid)
      IF(i_ncstat /= NF90_NOERR) &
        CALL grid_error(c_error='[bathy_initialize]: could not identify dimension')
      i_ncstat = nf90_inquire_dimension(i_fileid, i_dimid, len=i_bathydim)
      IF(i_ncstat /= NF90_NOERR) &
        CALL grid_error(c_error='[bathy_initialize]: could not read dimension')
      IF(i_bathydim /= 2_GRID_SI) &
        CALL grid_error(c_error='[bathy_initialize]: dimension mismatch')

      i_ncstat = nf90_inq_dimid(i_fileid, 'xysize', i_dimid)
      IF(i_ncstat /= NF90_NOERR) &
        CALL grid_error(c_error='[bathy_initialize]: could not read array size')
      i_ncstat= nf90_inquire_dimension(i_fileid, i_dimid, len=i_bathysize1D)
      IF(i_ncstat /= NF90_NOERR) &
        CALL grid_error(c_error='[bathy_initialize]: could not read array size')

!--- allocate raw data array
      ALLOCATE(r_rawbathy(i_bathysize1D), stat=i_alct)
      IF(i_alct /= 0) &
        CALL grid_error(c_error='[bathy_initialize]: could not allocate raw data field')

!--- read raw data array
      i_ncstat = nf90_inq_varid(i_fileid, 'z', i_bathyZ)
      IF(i_ncstat /= NF90_NOERR) &
        CALL grid_error(c_error='[bathy_initialize]: could not determine bathy z varid')
      i_ncstat = nf90_get_var(i_fileid, i_bathyZ, r_rawbathy)
      IF(i_ncstat /= NF90_NOERR) &
        CALL grid_error(c_error='[bathy_initialize]: could not read bathy z data')

!--- reshape data into 2D data array
!           invert sign, since bathymetry should appear positive, and elevation negative!
      k_cnt = 1
      DO j_cnt= i_bathydim2Dvals(2), 1, -1
        DO i_cnt= 1, i_bathydim2Dvals(1) !, 1, -1
          r_topo(i_cnt,j_cnt) = r_rawbathy(k_cnt)
          k_cnt               = k_cnt + 1
        END DO
      END DO

!--- deallocate raw array
      DEALLOCATE(r_rawbathy)

!--- read coards data
    CASE(1)
!--- read topography data
      i_ncstat = nf90_inq_varid(i_fileid, 'altitude', i_varid)
      IF(i_ncstat /= NF90_NOERR) &
        CALL grid_error(c_error='[bathy_initialize]: could not determine uplift varid')
      i_ncstat = nf90_get_var(i_fileid, i_varid, r_topo)
      IF(i_ncstat /= NF90_NOERR) &
        CALL grid_error(c_error='[bathy_initialize]: could not read uplift array')

    END SELECT

!--- close topography file
    i_ncstat = nf90_close(i_fileid)
    IF(i_ncstat /= NF90_NOERR) &
      CALL grid_error(c_error='[bathy_initialize]: could not close bathy data file')
    IF(GRID_parameters%iolog > 0) &
      write(GRID_parameters%iolog,*) 'INFO: Closed bathymetry data file'

!--- smooth bathymetry for better numerical stability...
    IF(i_bsmooth > 0_GRID_SI) THEN
      CALL smooth_bathymetry(i_bsmooth)
    END IF

    l_bathyinit = .TRUE.

  END SUBROUTINE bathy_initialize

!*******************************************************************************
! DESCRIPTION of [FUNCTION bathy]:
!> @brief calculates a bathymetry/topography height value
!>
!> @param         r_coos      coordinate array
!> @return                    value at postition r_coo
!
  FUNCTION bathy(r_coos) RESULT (r_val)

    IMPLICIT NONE

    REAL (KIND = GRID_SR), DIMENSION(GRID_dimension)  :: r_coos
    REAL (KIND = GRID_SR)                             :: r_val

!--- local declarations
    INTEGER (KIND = GRID_SI)                          :: i_xind, i_yind, i_cnt
    REAL (KIND = GRID_SR)                             :: r_ax, r_ay, r_xx1, r_xx2

!--- check if bathymetry data is loaded
    IF(.NOT. l_bathyinit) THEN
      CALL grid_error(c_error='[bathy]: bathymetry must be initialized before!')
    END IF

!--- interpolate bathy to given x, y pair
    inside: IF((r_coos(1) >= r_lon(1) .AND. r_coos(1) <= &
                r_lon(i_nlonpts)) .AND. &
               (r_coos(2) >= r_lat(1) .AND. r_coos(2) <= &
                r_lat(i_nlatpts))) THEN

!--- determine grid box surrounding r_coos and calculate weights for linear interpolation
!--- for longitude
      DO i_cnt = 1, i_nlonpts
        IF (r_lon(i_cnt+1)>=r_coos(1)) THEN
          i_xind = i_cnt
          EXIT
        END IF
      END DO
      IF (r_coos(1) > r_lon(i_nlonpts)) THEN
        r_ax   = 1.0_GRID_SR
        i_xind = i_nlonpts-1
      ELSE
        r_ax  = (r_coos(1) - r_lon(i_xind)) / (r_lon(i_xind+1) - r_lon(i_xind))
      END IF
!--- for latitude
      DO i_cnt = 1, i_nlatpts
        IF (r_lat(i_cnt+1)>=r_coos(2)) THEN
          i_yind = i_cnt
          EXIT
        END IF
      END DO
      IF (r_coos(2) > r_lat(i_nlatpts)) THEN
        r_ay   = 1.0_GRID_SR
        i_yind = i_nlonpts-1
      ELSE
        r_ay  = (r_coos(2) - r_lat(i_yind)) / (r_lat(i_yind+1) - r_lat(i_yind))
      END IF

      r_xx1 = (1._GRID_SR-r_ax) * r_topo(i_xind,i_yind  ) + r_ax * r_topo(i_xind+1,i_yind  )
      r_xx2 = (1._GRID_SR-r_ax) * r_topo(i_xind,i_yind+1) + r_ax * r_topo(i_xind+1,i_yind+1)

      r_val = (1._GRID_SR-r_ay) * r_xx1 + r_ay * r_xx2
    ELSE inside
      CALL grid_error(c_error='[bathy]: input coordinates are not inside of bathymetry range')
    END IF inside

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

!    INCLUDE "netcdf.inc"

    REAL (KIND = GRID_SR), OPTIONAL                     :: r_east
    REAL (KIND = GRID_SR), OPTIONAL                     :: r_west
    REAL (KIND = GRID_SR), OPTIONAL                     :: r_north
    REAL (KIND = GRID_SR), OPTIONAL                     :: r_south
    REAL (KIND = GRID_SR), OPTIONAL                     :: r_xres
    REAL (KIND = GRID_SR), OPTIONAL                     :: r_yres
    INTEGER (KIND = GRID_SI), OPTIONAL                  :: i_xsize
    INTEGER (KIND = GRID_SI), OPTIONAL                  :: i_ysize

!--- local declarations

!--- check if bathymetry data is loaded
    IF(.NOT. l_bathyinit) THEN
      CALL grid_error(c_error='[bathy_inquire]: bathymetry must be initialized before!')
    END IF

!--- fill optional output parameters
    IF(present(r_east))  r_east  = r_lon(i_nlonpts)
    IF(present(r_west))  r_west  = r_lon(1)
    IF(present(r_north)) r_north = r_lat(i_nlatpts)
    IF(present(r_south)) r_south = r_lat(1)
    IF(present(r_xres))  r_xres  = (r_lon(i_nlonpts) - r_lon(1)) / REAL(i_nlonpts, GRID_SR)
    IF(present(r_yres))  r_yres  = (r_lat(i_nlatpts) - r_lat(1)) / REAL(i_nlatpts, GRID_SR)
    IF(present(i_xsize)) i_xsize = i_nlonpts
    IF(present(i_ysize)) i_ysize = i_nlatpts

  END SUBROUTINE bathy_inquire

!*******************************************************************************
! DESCRIPTION of [FUNCTION get_shorepoints]:
!> @brief returns coordinates of nearshore points in bathymetry data
!>
!> @param[in]     r_shoreband       depth values to define band near shore
!> @param[in]     r_shorepoints     coordinates of nearshore points
!
  SUBROUTINE get_shorepoints(r_shoreband, r_shorepoints)

    IMPLICIT NONE

    REAL (KIND = GRID_SR), DIMENSION(GRID_dimension), INTENT(in)    :: r_shoreband
    REAL (KIND = GRID_SR), DIMENSION(:,:), ALLOCATABLE, INTENT(out) :: r_shorepoints

!--- local declarations
    INTEGER (KIND = GRID_SI)                                        :: i_nshorepts, i_shorept, &
                                                                       i_cnt, j_cnt

!--- check if bathymetry data is loaded
    IF(.NOT. l_bathyinit) THEN
      CALL grid_error(c_error='[get_shorepoints]: bathymetry must be initialized before!')
    END IF

    i_nshorepts = COUNT((r_shoreband(2) <= r_topo) .AND. (r_topo <= r_shoreband(1)))

    ALLOCATE(r_shorepoints(GRID_dimension, i_nshorepts))

    i_shorept = 1_GRID_SI
    DO i_cnt=1, i_nlonpts
      DO j_cnt=1, i_nlatpts
        IF ((r_shoreband(2) <= r_topo(i_cnt, j_cnt)) .AND. (r_topo(i_cnt, j_cnt) <= r_shoreband(1))) THEN
          r_shorepoints(1, i_shorept) = r_lon(i_cnt)
          r_shorepoints(2, i_shorept) = r_lat(j_cnt)
          i_shorept = i_shorept + 1_GRID_SI
        END IF
      END DO
    END DO

  END SUBROUTINE get_shorepoints

!*******************************************************************************
! DESCRIPTION of [SUBROUTINE smooth_bathymetry]:
!> @brief smoothes bathymetry field
!>
!> @param         i_iter      no. of iteration for next neighbor smoothing
!>
!> @note  other fields are used from module-global memory
!
  SUBROUTINE smooth_bathymetry(i_iter)

    IMPLICIT NONE

    INTEGER (KIND = GRID_SI)                            :: i_iter

!--- local declarations
    INTEGER (KIND = GRID_SI)                            :: i_cnt, j_cnt, i_it, i_alct
    REAL (KIND = GRID_SR)                               :: r_weig
    REAL (KIND = GRID_SR), DIMENSION(:,:), ALLOCATABLE  :: r_topotmp

    ALLOCATE(r_topotmp(i_nlonpts,i_nlatpts), stat=i_alct)
    IF(i_alct /= 0) &
      CALL grid_error(c_error='[smooth_bathymetry]: could not allocate 2D data field')

!--- use a five-point stencil to smooth the bathymetry field
!           repeat this iter times
    r_weig = 1._GRID_SR/8._GRID_SR
!    r_weig = 1._GRID_SR/4._GRID_SR
    iter_loop: DO i_it=1_GRID_SI,i_iter
      r_topotmp = r_topo
      DO j_cnt= 2_GRID_SI, i_nlatpts-1_GRID_SI
        DO i_cnt= 2_GRID_SI, i_nlonpts-1_GRID_SI
          r_topo(i_cnt,j_cnt) = r_weig*(4._GRID_SR*r_topotmp(i_cnt,j_cnt)+ &
            r_topotmp(i_cnt+1_GRID_SI,j_cnt)+r_topotmp(i_cnt-1_GRID_SI,j_cnt)+ &
            r_topotmp(i_cnt,j_cnt+1_GRID_SI)+r_topotmp(i_cnt,j_cnt-1_GRID_SI))
!          r_topo(i_cnt,j_cnt)= r_weig*( &
!            r_topotmp(i_cnt+1_GRID_SI,j_cnt)+r_topotmp(i_cnt-1_GRID_SI,j_cnt)+ &
!            r_topotmp(i_cnt,j_cnt+1_GRID_SI)+r_topotmp(i_cnt,j_cnt-1_GRID_SI))
        END DO
      END DO
    END DO iter_loop

!--- deallocate array
    DEALLOCATE(r_topotmp)

  END SUBROUTINE smooth_bathymetry

!*******************************************************************************
END MODULE MISC_bathy
