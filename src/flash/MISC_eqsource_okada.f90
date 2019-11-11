!*******************************************************************************
!
!> @file MISC_eqsource_okada.f90
!> @brief contains module MISC_eqsource
!
!
!*******************************************************************************
!
! VERSION(S):
!  1. original version                       m. baensch    06/2016
!  2. changed structure for input file       m. baensch    07/2016
!  3. added option to write uplift file      m. baensch    09/2016
!
!*******************************************************************************
! MODULE DESCRIPTION:
!> @brief calculates resulting seafloor displacement computed by the Okada (1985)
!>        model from given finite fault data
!>
!>  @note the input file must have the following structure:
!>  Header/comment lines are denoted by a #
!>  First line has to be NUM_FAULTS followed by a line defining the number of
!>  subfaults.
!>  Preceding the data, there has to be a line with FAULT_PARAMETERS indicating
!>  that the data follows in the next line.
!>  The data then has to yield NUM_FAULTS rows with i_faultdim entries (usually
!>  10 or 12).
!>
!>  The columns have to be aranged in the following way (12 inputs):
!>   1) latitude in deg,
!>   2) longitude in deg,
!>   3) depth in m,
!>   4) slip in m,
!>   5) rake in deg,
!>   6) strike in deg,
!>   7) dip in deg,
!>   8) intial rupture time in s,
!>   9) starting rise time in s,
!>  10) ending rise time in s,
!>  11) length of rupture in m and
!>  12) width of rupture in m
!>
!>  As of yet, initial rupture time and rising times have not been implemented
!>  (07/2016)
!
MODULE MISC_eqsource

  USE GRID_api
  USE netcdf

!--- global constants
    REAL (KIND = GRID_SR), PARAMETER                      :: r_earth    = 6.3675e6_GRID_SR
    REAL (KIND = GRID_SR), PARAMETER                      :: r_poisson  = 0.25_GRID_SR ! Poisson ratio for Okada
    INTEGER (KIND = GRID_SI), PARAMETER                   :: i_faultdim = 12
    INTEGER (KIND = GRID_SI)                              :: i_nlonpts, i_nlatpts, i_ntimepts
    REAL (KIND = GRID_SR)                                 :: r_dlon, r_dlat
    REAL (KIND = GRID_SR), DIMENSION(:,:,:), ALLOCATABLE  :: r_Z
    REAL (KIND = GRID_SR), DIMENSION(:), ALLOCATABLE      :: r_lat
    REAL (KIND = GRID_SR), DIMENSION(:), ALLOCATABLE      :: r_lon
    REAL (KIND = GRID_SR), DIMENSION(:), ALLOCATABLE      :: r_time

! DESCRIPTION of [TYPE fault]:
!> @brief type definition for a fault description
!>
  TYPE fault
    INTEGER (KIND = GRID_SI)                        :: i_nsubfaults     !< number of subfaults
    REAL (KIND = GRID_SR), DIMENSION(:,:), POINTER  :: r_subfaults      !< subfault parameters
  END TYPE fault

  PRIVATE
  PUBLIC  :: initialize_source, compute_uplift, write_uplift_coards

  CONTAINS

!*******************************************************************************
! DESCRIPTION of [SUBROUTINE initialize_source]:
!> @brief initializes the source and uplift
!>
!> @param[in]     c_filename       uplift data input file
!> @param[in]     r_lonrange       lon. range in which the uplift is calculated
!> @param[in]     r_latrange       lat. range in which the uplift is calculated
!> @param[in]     i_nlonpoints     number of longitude points
!> @param[in]     i_nlatpoints     number of latitude points
!> @param[in]     i_boundwidth     size of boundary, where source is linearly phased out
!
!> @note for the time being, we do not consider time-dependent uplifts!
!
  SUBROUTINE initialize_source(c_filename, r_lonrange, r_latrange, &
                               i_nlonpoints, i_nlatpoints, i_boundwidth)

    IMPLICIT NONE

    CHARACTER (len=*), INTENT(in)                         :: c_filename
    INTEGER (KIND = GRID_SI), INTENT(IN)                  :: i_boundwidth
    INTEGER (KIND = GRID_SI), INTENT(IN)                  :: i_nlonpoints, i_nlatpoints
    REAL (KIND = GRID_SR), DIMENSION(2), INTENT(IN)       :: r_lonrange, r_latrange

!--- local declarations
    TYPE (fault)                                          :: p_fault
    INTEGER (KIND = GRID_SI)                              :: i_iofil
    INTEGER (KIND = GRID_SI), PARAMETER                   :: i_ioerr = 0
    INTEGER (KIND = GRID_SI)                              :: i_iost, i_ioend, i_alct
    INTEGER (KIND = GRID_SI)                              :: i_cnt, j_cnt
    INTEGER (KIND = GRID_SI), DIMENSION(1)                :: i_pos
    CHARACTER (len=80)                                    :: c_filrow
    REAL (KIND = GRID_SR), DIMENSION(:,:), POINTER        :: r_subfaults_tmp
    INTEGER (KIND = GRID_SI), DIMENSION(:), ALLOCATABLE   :: i_faulttimeidx, i_timesorted
    REAL (KIND = GRID_SR), DIMENSION(:), ALLOCATABLE      :: r_timeunsorted
    LOGICAL, DIMENSION(:), ALLOCATABLE                    :: l_positions

    WRITE(*,'(A)') 'Initializing source. Please wait.'

    i_nlonpts = i_nlonpoints
    i_nlatpts = i_nlatpoints

!--- load fault
!--- open file in subfault format for uplift parameters
    open(NEWUNIT=i_iofil, file=c_filename, status='OLD', action='READ', iostat=i_iost)
    file_notopen: IF(i_iost /= 0) THEN
      write(i_ioerr,*) 'ERROR: Filename: ', c_filename
      IF(GRID_parameters%iolog > 0) &
        write(GRID_parameters%iolog,*) 'ERROR: Filename: ', c_filename
      CALL grid_error(c_error='[initialize_source]: Could not open file for uplift parameters')
    ELSE file_notopen
      IF(GRID_parameters%iolog > 0) THEN
        write(GRID_parameters%iolog,*) 'INFO: Filename: ', c_filename, ' opened on unit: ', i_iofil
      END IF
    END IF file_notopen

    p_fault%i_nsubfaults = 0

!--- read line by line
    read_loop: DO
      read(i_iofil,'(A)',iostat=i_ioend) c_filrow

!--- if file ended
      file_end: IF(i_ioend /= 0) THEN
        close(i_iofil)
        IF(GRID_parameters%iolog > 0) &
          write(GRID_parameters%iolog,*) 'INFO: Closed file on unit: ', i_iofil
        EXIT
      ELSE file_end

!--- decide what to DO with line according to first character
        comment_line: IF(c_filrow(1:1) == '#' .OR. c_filrow(1:1) == '!') THEN
          CYCLE read_loop
        ELSE IF (TRIM(c_filrow) == 'NUM_FAULTS') THEN
          READ(i_iofil,*) p_fault%i_nsubfaults
          ALLOCATE(r_subfaults_tmp(i_faultdim, p_fault%i_nsubfaults), &
                   r_timeunsorted(p_fault%i_nsubfaults), &
                   i_faulttimeidx(p_fault%i_nsubfaults), STAT=i_alct)
            IF(i_alct /= 0) CALL grid_error(c_error='[initialize_source]: could not allocate faults array etc.')
            i_faulttimeidx = 0

        ELSE IF (TRIM(c_filrow) == 'FAULT_PARAMETERS' .AND. p_fault%i_nsubfaults>0) THEN
          i_ntimepts = 0
          DO i_cnt=1,p_fault%i_nsubfaults
            READ(i_iofil,*) r_subfaults_tmp(:,i_cnt)

!--- determine time slices
            DO j_cnt = 1, i_ntimepts
              IF (r_subfaults_tmp(8,i_cnt) == r_timeunsorted(j_cnt)) THEN
                i_faulttimeidx(i_cnt) = j_cnt
                EXIT
              END IF
            END DO ! j_cnt
            IF (i_faulttimeidx(i_cnt) == 0) THEN
              i_ntimepts                 = i_ntimepts + 1
              r_timeunsorted(i_ntimepts) = r_subfaults_tmp(8,i_cnt)
              i_faulttimeidx(i_cnt)      = i_ntimepts
            END IF
          END DO ! i_cnt
        ELSE comment_line
          CALL grid_error(c_error='[initialize_source]: Dont know what to do with this parameter file')
        END IF comment_line
      END IF file_end
    END DO read_loop

    p_fault%r_subfaults => r_subfaults_tmp
    NULLIFY(r_subfaults_tmp)

!--- compute order of time slices
    ALLOCATE(i_timesorted(i_ntimepts), r_time(i_ntimepts), &
             l_positions(i_ntimepts), STAT=i_alct)
    IF(i_alct /= 0) CALL grid_error(c_error='[initialize_source]: could not allocate time array')
    l_positions = .TRUE.

    DO i_cnt = 1, i_ntimepts
      i_pos = MINLOC(r_timeunsorted(1:i_ntimepts), l_positions)
      i_timesorted(i_pos(1)) = i_cnt
      r_time(i_cnt)       = r_timeunsorted(i_pos(1))
      l_positions(i_pos(1))  = .FALSE.
    END DO ! i_cnt

!    WRITE(*,'(A, I6)') '#times: ', i_ntimepts
!    WRITE(*,'(A)') 'times unsorted:'
!    WRITE(*,'(10F12.5)') r_timeunsorted(1:i_ntimepts)
!    WRITE(*,'(A)') 'times sorted  :'
!    WRITE(*,'(10F12.5)') r_time
!    WRITE(*,'(A)') 'i_timesorted:'
!    WRITE(*,'(20I5)') i_timesorted
!    WRITE(*,'(A)') 'fault timeidx:'
!    WRITE(*,'(20I5)') i_faulttimeidx

!--- allocate rupture time array
    ALLOCATE(r_Z(i_nlonpts,i_nlatpts,i_ntimepts), r_lat(i_nlatpts),&
             r_lon(i_nlonpts), STAT=i_alct)
    IF(i_alct /= 0) &
    CALL grid_error(c_error='[initialize_source]: could not allocate uplift array')

!--- compute resulting uplift on uniform grid now
!--- initialize lat/lon grid
    r_dlon = (r_lonrange(2)-r_lonrange(1)) /(i_nlonpts-1) ! longitude resolution in minutes
    r_dlat = (r_latrange(2)-r_latrange(1)) /(i_nlatpts-1) ! latitude resolution in minutes

    DO i_cnt = 1, i_nlatpts
      r_lat(i_cnt) = r_latrange(1) + REAL(i_cnt-1_GRID_SI, GRID_SR) * r_dlat
    END DO
    DO i_cnt = 1, i_nlonpts
      r_lon(i_cnt) = r_lonrange(1) + REAL(i_cnt-1_GRID_SI, GRID_SR) * r_dlon
    END DO

    r_Z = 0.0_GRID_SR

    IF(GRID_parameters%iolog > 0) &
      WRITE(GRID_parameters%iolog,*) 'INFO: Computing Okada model with ', p_fault%i_nsubfaults, 'subfaults.'

    subfault_loop: DO i_cnt= 1, p_fault%i_nsubfaults
      WRITE(*,'(A, I3, A, I3)') &
        'Computing contributon from ', i_cnt, 'th subfault, out of ', p_fault%i_nsubfaults
!       WRITE(*,'(A,I6)') 'time slab: ', i_timesorted(i_faulttimeidx(i_cnt))
      CALL add_uplift(p_fault%r_subfaults(:,i_cnt), i_timesorted(i_faulttimeidx(i_cnt)))
    END DO subfault_loop
    WRITE(*,*)

    DEALLOCATE(r_timeunsorted, i_timesorted, i_faulttimeidx, l_positions)

    IF(GRID_parameters%iolog > 0) &
      WRITE(GRID_parameters%iolog,*) 'INFO: Okada model computed.'

!--- smooth uplift at boundaries
    DO i_cnt = 1, i_boundwidth
      DO j_cnt = i_cnt, i_nlatpts+1-i_cnt
        r_Z(i_cnt            ,j_cnt, :) = REAL(i_cnt,GRID_SR)/REAL(i_boundwidth+1,GRID_SR) * r_Z(i_cnt            ,j_cnt,:)
        r_Z(i_nlonpts+1-i_cnt,j_cnt, :) = REAL(i_cnt,GRID_SR)/REAL(i_boundwidth+1,GRID_SR) * r_Z(i_nlonpts+1-i_cnt,j_cnt,:)
      END DO
      DO j_cnt = i_cnt+1, i_nlonpts-i_cnt
        r_Z(j_cnt,i_cnt            , :) = REAL(i_cnt,GRID_SR)/REAL(i_boundwidth+1,GRID_SR) * r_Z(j_cnt,i_cnt,:)
        r_Z(j_cnt,i_nlatpts+1-i_cnt, :) = REAL(i_cnt,GRID_SR)/REAL(i_boundwidth+1,GRID_SR) * r_Z(j_cnt,i_nlatpts+1-i_cnt,:)
      END DO
    END DO

! ! DEBUG:
!     OPEN(unit=13, file='displ.txt', status='NEW', action='WRITE')
!     DO i_cnt=1,i_nlat
!       WRITE(13,*) (r_Z(j_cnt,i_cnt,1), j_cnt=1,i_nlon)
!     END DO
!     CLOSE(13)

  END SUBROUTINE initialize_source

!*******************************************************************************
! DESCRIPTION of [FUNCTION compute_uplift]:
!> @brief computes the uplift and interpolates according to the Okada model (1985)
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
    INTEGER (KIND = GRID_SI)                              :: i_xind, i_yind, i_t, i_timept
    REAL (KIND = GRID_SR)                                 :: r_ys, r_xw, &
      r_ax, r_ay, r_xx1, r_xx2, r_at

!--- interpolate Okada uplift to given x,y pair
    inside: IF((r_coord(1) >= r_lon(1) .AND. r_coord(1) <= r_lon(i_nlonpts)) .AND. &
               (r_coord(2) >= r_lat(1) .AND. r_coord(2) <= r_lat(i_nlatpts)) .AND. &
               (r_timeuplift >= r_time(1))) THEN

!--- determine grid box surounding r_coord
      i_xind = INT((r_coord(1)-r_lon(1))/r_dlon, GRID_SI) + 1_GRID_SI
      i_yind = INT((r_coord(2)-r_lat(1))/r_dlat, GRID_SI) + 1_GRID_SI

      IF (r_coord(1) == r_lon(i_nlonpts)) THEN
        i_xind = i_xind - 1
      END IF
      IF (r_coord(2) == r_lat(i_nlatpts)) THEN
        i_yind = i_yind - 1
      END IF

      r_xw = r_lon(i_xind)
      r_ys = r_lat(i_yind)

!--- compute weights for linear interpolation
      r_ax = (r_coord(1) - r_xw) / r_dlon
      r_ay = (r_coord(2) - r_ys) / r_dlat

!--- determine time snapshot to be used
      DO i_timept=1,i_ntimepts
        IF (r_time(i_timept)<=r_timeuplift) THEN
          i_t = i_timept
        ELSE
          EXIT
        END IF
      END DO

!--- interpolate linearly along x, y and t directions
      r_xx1 = (1._GRID_SR-r_ax) * r_Z(i_xind,i_yind  ,i_t) + r_ax * r_Z(i_xind+1,i_yind  ,i_t)
      r_xx2 = (1._GRID_SR-r_ax) * r_Z(i_xind,i_yind+1,i_t) + r_ax * r_Z(i_xind+1,i_yind+1,i_t)

      r_val = (1._GRID_SR-r_ay) * r_xx1 + r_ay * r_xx2
    ELSE inside
      r_val = 0.0_GRID_SR
    END IF inside

  END FUNCTION compute_uplift

!*******************************************************************************
! DESCRIPTION of [SUBROUTINE write_uplift_coards]:
!> @brief writes the uplift into a COARDS file
!>
!> @param[in]     c_filename       uplift data input file
!> @param[in]     c_tunits         units for time coordinate
!> @param[in]     c_xunits         units for x coordinate
!> @param[in]     c_yunits         units for y coordinate
!> @param[in]     c_zunits         units for uplift
!
  SUBROUTINE write_uplift_coards(c_filename, c_tunits, c_xunits, c_yunits, c_zunits)

    IMPLICIT NONE

!--- includes
    INCLUDE "netcdf.inc"

    CHARACTER (LEN = *), INTENT(IN)                   :: c_filename
    CHARACTER (LEN = *), INTENT(IN), OPTIONAL         :: c_tunits
    CHARACTER (LEN = *), INTENT(IN), OPTIONAL         :: c_xunits
    CHARACTER (LEN = *), INTENT(IN), OPTIONAL         :: c_yunits
    CHARACTER (LEN = *), INTENT(IN), OPTIONAL         :: c_zunits

!--- local declarations
    INTEGER                                           :: i_ncid, i_ncstatus, &
      i_dimid_x, i_dimid_y, i_dimid_time, i_varid_x, i_varid_y, i_varid_time, i_varid_dz
    INTEGER, DIMENSION(8)                             :: i_dateandtime
    CHARACTER (LEN = 43)                              :: c_title = 'computed output for uplift from Okada model'
    CHARACTER (LEN = 64)                              :: c_tmp

!--- create the file for timestep data
    i_ncstatus = nf90_create(TRIM(c_filename), NF90_CLOBBER, i_ncid)
    IF (i_ncstatus /= NF90_NOERR) CALL netcdf_error(i_ncstatus)
    IF(GRID_parameters%iolog > 0) &
      WRITE(GRID_parameters%iolog,*) 'INFO: Created uplift data file for writing: ', c_filename

!--- put global attributes into the file
    i_ncstatus = nf90_put_att(i_ncid, NF90_GLOBAL, 'title', c_title)
    IF (i_ncstatus /= NF90_NOERR) CALL netcdf_error(i_ncstatus)
    i_ncstatus = nf90_put_att(i_ncid, NF90_GLOBAL, 'Conventions', 'COARDS')
    IF (i_ncstatus /= NF90_NOERR) CALL netcdf_error(i_ncstatus)

!--- define dimensions
    i_ncstatus = nf90_def_dim(i_ncid, 'x', i_nlonpts, i_dimid_x)
    IF (i_ncstatus /= NF90_NOERR) CALL netcdf_error(i_ncstatus)
    i_ncstatus = nf90_def_dim(i_ncid, 'y', i_nlatpts, i_dimid_y)
    IF (i_ncstatus /= NF90_NOERR) CALL netcdf_error(i_ncstatus)
    i_ncstatus = nf90_def_dim(i_ncid, 'time', i_ntimepts, i_dimid_time)
    IF (i_ncstatus /= NF90_NOERR) CALL netcdf_error(i_ncstatus)

!--- define data variables and write variable attribute data into netCDF file
    i_ncstatus = nf90_def_var(i_ncid, 'x', NF90_REAL, i_dimid_x, i_varid_x)
    IF (i_ncstatus /= NF90_NOERR) CALL netcdf_error(i_ncstatus)
    IF (PRESENT(c_xunits)) THEN
      i_ncstatus = nf90_put_att(i_ncid, i_varid_x, 'units', c_xunits)
    ELSE
      i_ncstatus = nf90_put_att(i_ncid, i_varid_x, 'units', 'degrees east')
    END IF
    IF (i_ncstatus /= NF90_NOERR) CALL netcdf_error(i_ncstatus)

    i_ncstatus = nf90_def_var(i_ncid, 'y', NF90_REAL,  i_dimid_y, i_varid_y)
    IF (i_ncstatus /= NF90_NOERR) CALL netcdf_error(i_ncstatus)
    IF (PRESENT(c_yunits)) THEN
      i_ncstatus = nf90_put_att(i_ncid, i_varid_y, 'units', c_yunits)
    ELSE
      i_ncstatus = nf90_put_att(i_ncid, i_varid_y, 'units', 'degrees north')
    END IF
    IF (i_ncstatus /= NF90_NOERR) CALL netcdf_error(i_ncstatus)

    i_ncstatus = nf90_def_var(i_ncid, 'time', NF90_REAL, i_dimid_time, i_varid_time)
    IF (i_ncstatus /= NF90_NOERR) CALL netcdf_error(i_ncstatus)
    IF (PRESENT(c_tunits)) THEN
      i_ncstatus = nf90_put_att(i_ncid, i_varid_time, 'units', c_tunits)
    ELSE
      i_ncstatus = nf90_put_att(i_ncid, i_varid_time, 'units', 'seconds since initial rupture')
    END IF
    IF (i_ncstatus /= NF90_NOERR) CALL netcdf_error(i_ncstatus)
    i_ncstatus = nf90_def_var(i_ncid, 'dz', NF90_REAL, (/ i_dimid_x, i_dimid_y, i_dimid_time /), i_varid_dz)
    IF (i_ncstatus /= NF90_NOERR) CALL netcdf_error(i_ncstatus)
    IF (PRESENT(c_zunits)) THEN
      i_ncstatus = nf90_put_att(i_ncid, i_varid_dz, 'units', c_zunits)
    ELSE
      i_ncstatus = nf90_put_att(i_ncid, i_varid_dz, 'units', 'm')
    END IF
    IF (i_ncstatus /= NF90_NOERR) CALL netcdf_error(i_ncstatus)

!--- write creation and modification date into netCDF file
    CALL DATE_AND_TIME(VALUES=i_dateandtime)
    WRITE(c_tmp, '(I4.4, A, I2.2, A, I2.2, A, I2.2, A, I2.2, A, I2.2, A, I3.2, A, I2.2)') &
      i_dateandtime(1), '-', i_dateandtime(2), '-', i_dateandtime(3), ' ', &
      i_dateandtime(5), ':', i_dateandtime(6), ':', i_dateandtime(7), ' ', &
      i_dateandtime(4)/60, ':', MOD(i_dateandtime(4),60)
    i_ncstatus = nf90_put_att(i_ncid, NF90_GLOBAL, 'creation date', c_tmp)
    IF (i_ncstatus /= NF90_NOERR) CALL netcdf_error(i_ncstatus)

!--- put the file into data mode
    i_ncstatus = nf90_enddef(i_ncid)
    IF (i_ncstatus /= NF90_NOERR) CALL netcdf_error(i_ncstatus)

!--- write variable data into netCDF file
    i_ncstatus = nf90_put_var(i_ncid, i_varid_x, r_lon)
    IF (i_ncstatus /= NF90_NOERR) CALL netcdf_error(i_ncstatus)
    i_ncstatus = nf90_put_var(i_ncid, i_varid_y, r_lat)
    IF (i_ncstatus /= NF90_NOERR) CALL netcdf_error(i_ncstatus)
    i_ncstatus = nf90_put_var(i_ncid, i_varid_time, r_time)
    IF (i_ncstatus /= NF90_NOERR) CALL netcdf_error(i_ncstatus)
    i_ncstatus = nf90_put_var(i_ncid, i_varid_dz, r_Z)
    IF (i_ncstatus /= NF90_NOERR) CALL netcdf_error(i_ncstatus)

!--- close NetCDF file
    i_ncstatus = nf90_close(i_ncid)
    IF (i_ncstatus /= NF90_NOERR) CALL netcdf_error(i_ncstatus)
    IF(GRID_parameters%iolog > 0) &
      WRITE(GRID_parameters%iolog,*) 'INFO: Closed uplift data file: ', c_filename

  END SUBROUTINE write_uplift_coards

!*******************************************************************************
! DESCRIPTION of [SUBROUTINE add_uplift]:
!> @brief Adds the uplift to the initialized source
!
!> @param[in]     r_subfaultparam   fault data from input file
!> @param[in]     i_timeslice       time slice
!
  SUBROUTINE add_uplift(r_subfaultparam, i_timeslice)

    IMPLICIT NONE

    REAL(KIND = GRID_SR), DIMENSION(:), INTENT(IN)      :: r_subfaultparam
    INTEGER (KIND = GRID_SI),           INTENT(IN)      :: i_timeslice

!--- local declarations
    REAL(KIND = GRID_SR)                                :: DEG2RAD, LAT2METER
    REAL(KIND = GRID_SR)                                :: r_lat0, r_lon0, r_slip, r_length, &
      r_width, r_strike, r_dip, r_depth, r_rake, TRU, TRI, MO
    REAL(KIND = GRID_SR)                                :: r_angstrike, r_angdip, r_angrake
    REAL(KIND = GRID_SR)                                :: r_latb, r_lonb, r_depb
    INTEGER (KIND = GRID_SI)                            :: i_cntlon, i_cntlat
    REAL(KIND = GRID_SR)                                :: r_dstrike, r_ddip
    REAL(KIND = GRID_SR)                                :: r_xx, r_yy, r_x1, r_x2, r_p, r_q
    REAL(KIND = GRID_SR)                                :: F1, F2, F3, F4, G1, G2, G3, G4
    REAL(KIND = GRID_SR)                                :: r_us, r_ud

    DEG2RAD   = GRID_PI / 180.0_GRID_SR
    LAT2METER = r_earth * DEG2RAD

    r_lat0   = r_subfaultparam( 1) ! origin of plate (latitude) at centroid
    r_lon0   = r_subfaultparam( 2) ! origin of plate (longitude) at centroid
    r_depth  = r_subfaultparam( 3) ! depth at centroid
    r_slip   = r_subfaultparam( 4) ! slip
    r_rake   = r_subfaultparam( 5) ! rake angle
    r_strike = r_subfaultparam( 6) ! strike angle
    r_dip    = r_subfaultparam( 7) ! dip angle
    TRU      = r_subfaultparam( 8) !
    TRI      = r_subfaultparam( 9) !
    MO       = r_subfaultparam(10) !
    r_length = r_subfaultparam(11) ! length
    r_width  = r_subfaultparam(12) ! width

!--- convert angles to radians:
    r_angstrike = DEG2RAD * r_strike
    r_angdip    = DEG2RAD * r_dip
    r_angrake   = DEG2RAD * r_rake

!--- Okada model assumes origin is at bottom center
    r_lonb = r_lon0  + 0.5_GRID_SR*r_width*COS(r_angdip)*COS(r_angstrike) / (LAT2METER*COS(DEG2RAD*r_lat0))
    r_latb = r_lat0  - 0.5_GRID_SR*r_width*COS(r_angdip)*SIN(r_angstrike) / (LAT2METER)
    r_depb = r_depth + 0.5_GRID_SR*r_width*SIN(r_angdip)

!--- displacement in direction of strike and dip
    r_dstrike = r_slip * COS(r_angrake)
    r_ddip    = r_slip * SIN(r_angrake)

    DO i_cntlon = 1, i_nlonpts
      DO i_cntlat = 1, i_nlatpts

!--- convert distance (r_lon,r_lat)-(r_lonb,r_latb) from degrees to meters
        r_xx = LAT2METER * (r_lon(i_cntlon) - r_lonb) * COS(DEG2RAD * r_lat(i_cntlat))
        r_yy = LAT2METER * (r_lat(i_cntlat) - r_latb)

!--- convert to distance along strike (x1) and dip (x2)
        r_x1 = r_xx * SIN(r_angstrike) + r_yy * COS(r_angstrike)
        r_x2 = r_xx * COS(r_angstrike) - r_yy * SIN(r_angstrike)

!--- in Okada's paper, x2 is distance up the fault plane, not down dip
        r_x2 = -r_x2

        r_p = r_x2 * COS(r_angdip) + r_depb * SIN(r_angdip)
        r_q = r_x2 * SIN(r_angdip) - r_depb * COS(r_angdip)

        F1 = strike_slip(r_x1+r_length/2.0_GRID_SR, r_p        , r_angdip, r_q)
        F2 = strike_slip(r_x1+r_length/2.0_GRID_SR, r_p-r_width, r_angdip, r_q)
        F3 = strike_slip(r_x1-r_length/2.0_GRID_SR, r_p        , r_angdip, r_q)
        F4 = strike_slip(r_x1-r_length/2.0_GRID_SR, r_p-r_width, r_angdip, r_q)
        G1 = dip_slip(   r_x1+r_length/2.0_GRID_SR, r_p        , r_angdip, r_q)
        G2 = dip_slip(   r_x1+r_length/2.0_GRID_SR, r_p-r_width, r_angdip, r_q)
        G3 = dip_slip(   r_x1-r_length/2.0_GRID_SR, r_p        , r_angdip, r_q)
        G4 = dip_slip(   r_x1-r_length/2.0_GRID_SR, r_p-r_width, r_angdip, r_q)

        r_us = (F1 - F2 - F3 + F4) * r_dstrike
        r_ud = (G1 - G2 - G3 + G4) * r_ddip

        r_Z(i_cntlon,i_cntlat,i_timeslice:i_ntimepts) = &
        r_Z(i_cntlon,i_cntlat,i_timeslice:i_ntimepts) + r_us + r_ud

      END DO
    END DO

  END SUBROUTINE add_uplift

!*******************************************************************************
!   Used for Okada's model
!   Methods from Yoshimitsu Okada (1985)
!*******************************************************************************
! DESCRIPTION of [FUNCTION strike_slip]:
!> @brief Calculates the strike slip used in Okada's model (1985)
!
!> @param[in]     y1              distance along strike
!> @param[in]     y2              distance up the fault plane
!> @param[in]     ang_dip         dip of fault
!> @param[in]     q               depth up the fault plane
!> @return                        displacement in z direction
!
  FUNCTION strike_slip(y1, y2, ang_dip, q) RESULT(f)

    IMPLICIT NONE

    REAL(KIND = GRID_SR), INTENT(IN)            :: y1, y2, ang_dip, q
    REAL(KIND = GRID_SR)                        :: f

    REAL(KIND = GRID_SR)                        :: sn, cs, d_bar, r, xx, a4

    sn    = SIN(ang_dip)
    cs    = COS(ang_dip)
    d_bar = y2*sn - q*cs
    r     = SQRT(y1**2 + y2**2 + q**2)
    xx    = SQRT(y1**2 + q**2)
    a4    = 2.0_GRID_SR*r_poisson/cs*(LOG(r+d_bar) - sn*LOG(r+y2))
    f     = -(d_bar*q/r/(r+y2) + q*sn/(r+y2) + a4*sn)/(2.0_GRID_SR*GRID_PI)

  END FUNCTION strike_slip

!*******************************************************************************
! DESCRIPTION of [FUNCTION dip_slip]:
!> @brief Calculates the dip slip used in Okada's model (1985)
!
!> @param[in]     y1              distance along strike
!> @param[in]     y2              distance up the fault plane
!> @param[in]     ang_dip         dip of fault
!> @param[in]     q               depth up the fault plane
!> @return                        displacement in z direction
!
  FUNCTION dip_slip(y1, y2, ang_dip, q) RESULT(f)

    IMPLICIT NONE

    REAL(KIND = GRID_SR), INTENT(IN)            :: y1, y2, ang_dip, q
    REAL(KIND = GRID_SR)                        :: f

    REAL(KIND = GRID_SR)                        :: sn, cs, d_bar, r, xx, a5

    sn = SIN(ang_dip)
    cs = COS(ang_dip)

    d_bar = y2*sn - q*cs;
    r     = SQRT(y1**2 + y2**2 + q**2)
    xx    = SQRT(y1**2 + q**2)
    a5    = 4.0_GRID_SR*r_poisson/cs*ATAN((y2*(xx+q*cs)+xx*(r+xx)*sn)/y1/(r+xx)/cs)
    f     = -(d_bar*q/r/(r+y1) + sn*ATAN(y1*y2/q/r) - a5*sn*cs)/(2.0_GRID_SR*GRID_PI)

  END FUNCTION dip_slip

!*******************************************************************************
! DESCRIPTION of [SUBROUTINE netcdf_error]:
!> @brief prints error message
!>
!> @param[in]       status
!
  SUBROUTINE netcdf_error(status)

    IMPLICIT NONE

!--- includes
    INCLUDE "netcdf.inc"

    INTEGER, INTENT (IN)            :: status

    IF (status /= NF90_NOERR) THEN
      print *, '[MISC_eqsource_okada:netcdf_error]: ', TRIM(nf90_strerror(status))
      STOP "Stopped"
    END IF

  END SUBROUTINE netcdf_error

!*******************************************************************************
END MODULE MISC_eqsource
