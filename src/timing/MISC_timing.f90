!*****************************************************************
!
! MODULE NAME:
!   MISC_timing
! FUNCTION:
!   miscellaneaus utility routines
! CONTAINS:
!   second, ticktock, stop_watch
!-----------------------------------------------------------------
!
! NAME:
!   stop_watch_init
! FUNCTION:
!   initialize time measurement
! SYNTAX:
!   CALL stop_watch_init(int, char, rtinfo)
! ON INPUT:
!   i_number:  number of timers in timer structure  INTEGER
!   c_strings: strings describing the code segments CHARACTER*16
!   p_tsinfo:  runtime info data structure          TYPE (rt_info)
! ON OUTPUT:
!   p_tsinfo: runtime info data structure (updated) TYPE (rt_info)
! CALLS:
!
! COMMENTS:
!
!-----------------------------------------------------------------
!
! NAME:
!   stop_watch
! FUNCTION:
!   measure time for a code segment and put info in runtime structure
! SYNTAX:
!   CALL stop_watch(char, int, rtinfo)
! ON INPUT:
!   c_action: start or stop the watch               CHARACTER*5
!   i_ident:  identification for the timer          INTEGER
!   p_tsinfo: runtime info data structure           TYPE (rt_info)
! ON OUTPUT:
!   p_tsinfo: runtime info data structure (updated) TYPE (rt_info)
! CALLS:
!
! COMMENTS:
!
!-----------------------------------------------------------------
! NAME:
!   second
! FUNCTION:
!   calculate actual second for timing
! SYNTAX:
!   real= second()
! ON INPUT:
!
! ON OUTPUT:
!   second: actual second measured from 1/1/1996  real
! CALLS:
!
! COMMENTS:
!   if no system clock is available, result is negative
!
!-----------------------------------------------------------------
! NAME:
!   ticktock
! FUNCTION:
!   fortran 90 portable timing routine
! SYNTAX:
!   real= ticktock(overhead= real, start= real)
! ON INPUT:
!   overhead: if present, overhead is subtracted  real
!   start:    first call to tick                  real
! ON OUTPUT:
!   overhead: if present, overhead is computed    real
!   ticktock: time of program segment             real
! CALLS:
!
! COMMENTS:
!   this routine has different functionalities:
! on 1st call:
!         r_sta = ticktock(overhead= r_ovh)
!     or
!         r_sta = ticktock()
!     gives back the starting tick mark for timing in r_sta and, if
!     requested, the overhead for calling the timing routine in r_ovh
! on second call:
!         r_tim= ticktock(start=r_sta)
!     or
!         r_tim= ticktock(start=r_sta, overhead= r_ovh)
!     gives back the time in [s] for the program segment in r_tim,
!     additionally, the overhead given by r_ovh is subtracted from total
!     time
!
! the convention is:
! noting is given:              first call
! only overhead is given:       first call
!       only start is given           following call
! overhead and start are given: following call
!
!-----------------------------------------------------------------
! PUBLIC:
!   stop_watch, stop_watch_init
! COMMENTS:
!
! USES:
!   MISC_globalparam
! LIBRARIES:
!
! REFERENCES:
!
! VERSION(S):
!  1. original version              j. behrens      07/96
!  2. stop_watch changed            j. behrens      01/97
!     stop_watch_init added
!  3. changed for self consistency  j. behrens      01/97
!  4. fixed ticktock                j. behrens      03/98
!  5. changed to be completely
!     self-consistent               j. behrens      12/2000
!
!*****************************************************************
MODULE MISC_timing
  INTEGER, PARAMETER :: i_ioerr= 0
  PRIVATE

!---------- structure for timing

  INTEGER, PARAMETER   :: DEF_timings= 10 ! number of different timers

  TYPE time_info
    CHARACTER (len=16) :: c_tim
    REAL               :: r_tim
    REAL               :: r_lap
  END TYPE time_info
  TYPE sw_info
    INTEGER                                  :: i_num
    TYPE (time_info), DIMENSION(DEF_timings) :: p_tim
  END TYPE sw_info

  PUBLIC :: stop_watch, stop_watch_init, time_info, sw_info, DEF_timings
  CONTAINS
!*****************************************************************
  SUBROUTINE stop_watch_init(i_number, c_strings, p_tsinfo)

!---------- local declarations

    IMPLICIT NONE
    INTEGER, INTENT(in)                                  :: i_number
    CHARACTER (len=16), INTENT(in), DIMENSION(i_number)  :: c_strings
    TYPE (sw_info), INTENT(inout)                        :: p_tsinfo
    INTEGER                                              :: i_cnt

!---------- initialize

    p_tsinfo%i_num = 0

!---------- consistency check

    check_con: IF(i_number > DEF_timings) THEN
      write(i_ioerr,*) 'TIMING: Initialization of timers omitted several'
      RETURN
    ELSE check_con

!---------- set number of timers in this structure

      p_tsinfo%i_num = min(DEF_timings, i_number)

!---------- initialize info structure

      main_loop: DO i_cnt=1, p_tsinfo%i_num
        p_tsinfo%p_tim(i_cnt)%c_tim= c_strings(i_cnt)
        p_tsinfo%p_tim(i_cnt)%r_tim= 0.0
        p_tsinfo%p_tim(i_cnt)%r_lap= 0.0
      END DO main_loop
    END IF check_con

    RETURN
  END SUBROUTINE stop_watch_init

!*****************************************************************
  SUBROUTINE stop_watch(c_action, i_ident, p_tsinfo)

!---------- local declarations

    IMPLICIT NONE
    CHARACTER (len= 5), INTENT(in) :: c_action
    INTEGER, INTENT(in)            :: i_ident
    TYPE (sw_info), INTENT(inout)  :: p_tsinfo
    REAL                           :: r_tmp

!---------- consistency check

    check_con: IF(i_ident > p_tsinfo%i_num) THEN
      write(i_ioerr,*) 'TIMING: Identification not correct, nothing done'
      RETURN
    END IF check_con

!---------- action start:

    start_stop: IF(c_action(1:4) == 'star') THEN
      p_tsinfo%p_tim(i_ident)%r_lap= ticktock()

!---------- action stop:

    ELSE IF(c_action(1:4) == 'stop') THEN start_stop
      r_tmp= ticktock(start= p_tsinfo%p_tim(i_ident)%r_lap)
      p_tsinfo%p_tim(i_ident)%r_tim= p_tsinfo%p_tim(i_ident)%r_tim+ r_tmp

!---------- not supported

    ELSE start_stop
      write(i_ioerr,*) 'TIMING: No supported action given: ', c_action
    END IF start_stop

    RETURN
  END SUBROUTINE stop_watch

!*****************************************************************
  FUNCTION second() RESULT (times)

!---------- local declarations

    IMPLICIT NONE

    REAL               :: times
!     CHARACTER (len=8)  :: a_date
    INTEGER, PARAMETER :: i_start= 19960101
    INTEGER, PARAMETER :: i_dsecs= 86400
    INTEGER, DIMENSION(10) :: i_time
    REAL, PARAMETER       :: r_mil=1/1000.
    REAL, PARAMETER       :: r_sec=1.
    REAL, PARAMETER       :: r_min=60.*r_sec
    REAL, PARAMETER       :: r_hor=60.*r_min
    REAL, PARAMETER       :: r_day=24.*r_hor
!     REAL, PARAMETER       :: r_yea=365.*r_day
!     INTEGER            :: i_day, i_mult
    INTEGER            :: i_cnt, i_rte
    INTEGER            :: i_total

!---------- CALL date_and_time for day info

! CALL date_and_time(date= a_date)
! date_ok: IF(a_date /= '        ') THEN
!   read(a_date,*) i_day
!   i_mult= abs(i_day- i_start)
! ELSE date_ok
!   i_mult= 1
! END IF date_ok

!---------- CALL date_and_time for easy time access

    CALL date_and_time(values=i_time)
    date_ok: IF(i_time(1) /= -HUGE(0)) THEN
      times= r_day* i_time(3)+ &
             r_hor* i_time(5)+ &
             r_min* i_time(6)+ &
             r_sec* i_time(7)+ &
             r_mil* i_time(8)
    ELSE date_ok
      times= -1.0
    END IF date_ok

!---------- CALL system_clock for maximum precision clock

! CALL system_clock(count= i_cnt, count_rate= i_rte)
! clock_ok: IF(i_rte /= 0) THEN
!   i_total= i_cnt+ i_mult* i_dsecs
!   i_total= i_cnt
!   times= float(i_total)/ float(i_rte)
! ELSE clock_ok
!   times= -1.0
! END IF clock_ok

    RETURN

  END FUNCTION second
!*****************************************************************
  FUNCTION ticktock(start, overhead) RESULT (timing)

!---------- local declarations

    IMPLICIT NONE

    REAL, OPTIONAL, INTENT(inout) :: overhead
    REAL, OPTIONAL, INTENT(in)    :: start
    REAL                          :: timing
    INTEGER, PARAMETER            :: i_loop= 10
    INTEGER                       :: i
    REAL                          :: r_tm, r_t1, r_tsum

!---------- initialize

    timing=0.0

!---------- check presence of input variables and act accordingly

!---------- second call to ticktock
    start_present: IF(present(start)) THEN
      overhead_present: IF(present(overhead)) THEN
        r_tm= second()- start- overhead
      ELSE overhead_present
        r_tm= second()- start
      END IF overhead_present

!---------- first call to ticktock
    ELSE start_present
      overhead_requested: IF(present(overhead)) THEN
        r_tsum= 0.0
        DO i=2, i_loop
          r_t1= second()
          r_tm= second()
          r_tsum= r_tsum+ (r_tm- r_t1)
        END DO
        overhead= r_tsum/float(i_loop)
        r_tm= second()
      ELSE overhead_requested
        r_tm= second()
      END IF overhead_requested
    END IF start_present

!---------- result

    timing= r_tm

    RETURN
  END FUNCTION ticktock
END MODULE MISC_timing
