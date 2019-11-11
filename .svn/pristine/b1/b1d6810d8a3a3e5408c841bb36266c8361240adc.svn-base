!*******************************************************************************
!
!> @file MISC_eqsource_coards.f90
!> @brief contains module MISC_eqsource
!
!*******************************************************************************
!
! VERSION(S):
!  1. original version                  m. baensch    09/2016
!
!*******************************************************************************
! MODULE DESCRIPTION:
!> @brief reads pre-calculated uplift from a COARDS-file

MODULE MISC_eqsource

  USE GRID_api
  USE netcdf

!--- global constants
    INTEGER (KIND = GRID_SI)                              :: i_nlonpts, i_nlatpts, i_ntimepts
    REAL (KIND = GRID_SR), DIMENSION(:,:,:), ALLOCATABLE  :: r_Z
    REAL (KIND = GRID_SR), DIMENSION(:), ALLOCATABLE      :: r_lat
    REAL (KIND = GRID_SR), DIMENSION(:), ALLOCATABLE      :: r_lon
    REAL (KIND = GRID_SR), DIMENSION(:), ALLOCATABLE      :: r_time

  PRIVATE
  PUBLIC  :: initialize_source, compute_uplift

  CONTAINS

!*******************************************************************************
! DESCRIPTION of [SUBROUTINE initialize_source]:
!> @brief initializes the source and uplift
!>
!> @param[in]     c_filename       uplift data input file
!> @param[in]     r_dummy1         dummy
!> @param[in]     r_dummy2         dummy
!> @param[in]     i_dummy1         dummy
!> @param[in]     i_dummy2         dummy
!> @param[in]     i_dummy3         dummy
!
!> @note the only actual input is the filename but all dummies are needed to keep
!>       the initialize_source routine consistent with other modules
!
  SUBROUTINE initialize_source(c_filename, r_dummy1, r_dummy2, &
                               i_dummy1, i_dummy2, i_dummy3)

    IMPLICIT NONE

    CHARACTER (len=*), INTENT(in)                             :: c_filename
    INTEGER (KIND = GRID_SI), INTENT(IN), OPTIONAL            :: i_dummy3
    INTEGER (KIND = GRID_SI), INTENT(IN), OPTIONAL            :: i_dummy1, i_dummy2
    REAL (KIND = GRID_SR), DIMENSION(2), INTENT(IN), OPTIONAL :: r_dummy1, r_dummy2

!--- local declarations
    INTEGER (KIND = GRID_SI)                            :: i_alct, i_cnt, j_cnt, k_cnt
    INTEGER (KIND = GRID_SI)                            :: i_fileid
    INTEGER (KIND = GRID_SI)                            :: i_dimid, i_varid
    INTEGER (KIND = GRID_SI)                            :: i_ncstat         ! NetCDF library return status

    WRITE(*,'(A)') 'Initializing source. Please wait.'

    IF (PRESENT(r_dummy1) .OR. PRESENT(r_dummy2) .OR. &
        PRESENT(i_dummy1) .OR. PRESENT(i_dummy2) .OR. &
        PRESENT(i_dummy3)) THEN

      IF(GRID_parameters%iolog > 0) &
        WRITE(GRID_parameters%iolog,*) 'INFO: Optional input in MISC_eqsource is set but ignored.'
    END IF


!--- open uplift file
    i_ncstat = nf90_open(c_filename, NF90_NOWRITE, i_fileid)
    IF(i_ncstat /= NF90_NOERR) &
      CALL grid_error(c_error='[initialize_source]: could not open uplift data file')
    IF(GRID_parameters%iolog > 0) &
      WRITE(GRID_parameters%iolog,*) 'INFO: Opened uplift data file for reading: ', c_filename

!--- determine array sizes.
    i_ncstat = nf90_inq_dimid(i_fileid, 'x', i_dimid)
    IF(i_ncstat /= NF90_NOERR) &
      CALL grid_error(c_error='[initialize_source]: could not identify x dimension')
    i_ncstat = nf90_inquire_dimension(i_fileid, i_dimid, len=i_nlonpts)
    IF(i_ncstat /= NF90_NOERR) &
      CALL grid_error(c_error='[initialize_source]: could not read x dimension')
    i_ncstat = nf90_inq_dimid(i_fileid, 'y', i_dimid)
    IF(i_ncstat /= NF90_NOERR) &
      CALL grid_error(c_error='[initialize_source]: could not identify y dimension')
    i_ncstat = nf90_inquire_dimension(i_fileid, i_dimid, len=i_nlatpts)
    IF(i_ncstat /= NF90_NOERR) &
      CALL grid_error(c_error='[initialize_source]: could not read y dimension')
    i_ncstat = nf90_inq_dimid(i_fileid, 'time', i_dimid)
    IF(i_ncstat /= NF90_NOERR) &
      CALL grid_error(c_error='[initialize_source]: could not identify time dimension')
    i_ncstat = nf90_inquire_dimension(i_fileid, i_dimid, len=i_ntimepts)
    IF(i_ncstat /= NF90_NOERR) &
      CALL grid_error(c_error='[initialize_source]: could not read time dimension')

!--- allocate data array
    ALLOCATE(r_lon(i_nlonpts), r_lat(i_nlatpts), r_time(i_ntimepts), &
             r_Z(i_nlonpts, i_nlatpts, i_ntimepts), stat=i_alct)
    IF(i_alct /= 0) &
      CALL grid_error(c_error='[initialize_source]: could not allocate data field')

!--- read data arrays
    i_ncstat = nf90_inq_varid(i_fileid, 'x', i_varid)
    IF(i_ncstat /= NF90_NOERR) &
      CALL grid_error(c_error='[initialize_source]: could not determine x varid')
    i_ncstat = nf90_get_var(i_fileid, i_varid, r_lon)
    IF(i_ncstat /= NF90_NOERR) &
      CALL grid_error(c_error='[initialize_source]: could not read x vector')

    i_ncstat = nf90_inq_varid(i_fileid, 'y', i_varid)
    IF(i_ncstat /= NF90_NOERR) &
      CALL grid_error(c_error='[initialize_source]: could not determine y varid')
    i_ncstat = nf90_get_var(i_fileid, i_varid, r_lat)
    IF(i_ncstat /= NF90_NOERR) &
      CALL grid_error(c_error='[initialize_source]: could not read y vector')

    i_ncstat = nf90_inq_varid(i_fileid, 'time', i_varid)
    IF(i_ncstat /= NF90_NOERR) &
      CALL grid_error(c_error='[initialize_source]: could not determine time varid')
    i_ncstat = nf90_get_var(i_fileid, i_varid, r_time)
    IF(i_ncstat /= NF90_NOERR) &
      CALL grid_error(c_error='[initialize_source]: could not read time vector')

    i_ncstat = nf90_inq_varid(i_fileid, 'dz', i_varid)
    IF(i_ncstat /= NF90_NOERR) &
      CALL grid_error(c_error='[initialize_source]: could not determine uplift varid')
    i_ncstat = nf90_get_var(i_fileid, i_varid, r_Z)
    IF(i_ncstat /= NF90_NOERR) &
      CALL grid_error(c_error='[initialize_source]: could not read uplift array')

!--- close uplift file
    i_ncstat = nf90_close(i_fileid)
    IF(i_ncstat /= NF90_NOERR) &
      CALL grid_error(c_error='[initialize_source]: could not close uplift data file')
    IF(GRID_parameters%iolog > 0) &
      WRITE(GRID_parameters%iolog,*) 'INFO: Closed uplift data file: ', c_filename

  END SUBROUTINE initialize_source

!*******************************************************************************
! DESCRIPTION of [FUNCTION compute_uplift]:
!> @brief computes the interpolated uplift at the given coordinates
!>
!> @param[in]     r_coord           coordinate values for uplift calculation
!> @param[in]     r_timeuplift      time at which uplift is needed
!> @return                          calculated uplift
!
  FUNCTION compute_uplift(r_coord, r_timeuplift) RESULT (r_val)

    IMPLICIT NONE

    REAL (KIND = GRID_SR), DIMENSION(GRID_dimension)      :: r_coord
    REAL (KIND = GRID_SR)                                 :: r_timeuplift
    REAL (KIND = GRID_SR)                                 :: r_val

!--- local declarations
    INTEGER (KIND = GRID_SI)                              :: i_xind, i_yind, i_t, i_timept, i_xy, i_cnt
    REAL (KIND = GRID_SR)                                 :: r_ax, r_ay, r_xx1, r_xx2, r_at

!--- interpolate uplift to given x, y pair
    inside: IF((r_coord(1) >= r_lon(1) .AND. r_coord(1) <= &
                r_lon(i_nlonpts)) .AND. &
               (r_coord(2) >= r_lat(1) .AND. r_coord(2) <= &
                r_lat(i_nlatpts)) .AND. &
               (r_timeuplift >= r_time(1))) THEN

!--- determine grid box surrounding r_coord and calculate weights for linear interpolation
!--- for loongitude
      DO i_cnt = 1, i_nlonpts
        IF (r_lon(i_cnt+1)>=r_coord(1)) THEN
          i_xind = i_cnt
          EXIT
        END IF
      END DO
      IF (r_coord(1) > r_lon(i_nlonpts)) THEN
        r_ax   = 1.0_GRID_SR
        i_xind = i_nlonpts-1
      ELSE
        r_ax  = (r_coord(1) - r_lon(i_xind)) / (r_lon(i_xind+1) - &
                 r_lon(i_xind))
      END IF
!--- for latitude
      DO i_cnt = 1, i_nlatpts
        IF (r_lat(i_cnt+1)>=r_coord(2)) THEN
          i_yind = i_cnt
          EXIT
        END IF
      END DO
      IF (r_coord(2) > r_lat(i_nlatpts)) THEN
        r_ay   = 1.0_GRID_SR
        i_yind = i_nlonpts-1
      ELSE
        r_ay  = (r_coord(2) - r_lat(i_yind)) / (r_lat(i_yind+1) - &
                r_lat(i_yind))
      END IF

      IF (i_ntimepts > 1) THEN
!--- determine time snapshot to be used
        DO i_timept=1,i_ntimepts-1
          IF (r_time(i_timept+1)>=r_timeuplift) THEN
            i_t = i_timept
            EXIT
          END IF
        END DO
        IF (r_timeuplift > r_time(i_ntimepts)) THEN
          r_at = 1.0_GRID_SR
          i_t  = i_ntimepts-1
        ELSE
          r_at = (r_timeuplift - r_time(i_t)) / (r_time(i_t+1) - r_time(i_t))
        END IF

!--- interpolate linearly along x, y and t directions
        r_xx1 = (1._GRID_SR-r_at) * ((1._GRID_SR-r_ax) * r_Z(i_xind,i_yind  ,i_t  ) + r_ax * r_Z(i_xind+1,i_yind  ,i_t  )) + &
                r_at              * ((1._GRID_SR-r_ax) * r_Z(i_xind,i_yind  ,i_t+1) + r_ax * r_Z(i_xind+1,i_yind  ,i_t+1))

        r_xx2 = (1._GRID_SR-r_at) * ((1._GRID_SR-r_ax) * r_Z(i_xind,i_yind+1,i_t  ) + r_ax * r_Z(i_xind+1,i_yind+1,i_t  )) + &
                r_at              * ((1._GRID_SR-r_ax) * r_Z(i_xind,i_yind+1,i_t+1) + r_ax * r_Z(i_xind+1,i_yind+1,i_t+1))
      ELSE
        r_xx1 = (1._GRID_SR-r_ax) * r_Z(i_xind,i_yind  ,1) + r_ax * r_Z(i_xind+1,i_yind  ,1)
        r_xx2 = (1._GRID_SR-r_ax) * r_Z(i_xind,i_yind+1,1) + r_ax * r_Z(i_xind+1,i_yind+1,1)
      END IF

      r_val = (1._GRID_SR-r_ay) * r_xx1 + r_ay * r_xx2
    ELSE inside
      r_val = 0.0_GRID_SR
    END IF inside

  END FUNCTION compute_uplift

!*******************************************************************************
END MODULE MISC_eqsource
