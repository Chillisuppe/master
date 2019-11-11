!*****************************************************************
!
! MODULE NAME:
!       MISC_system
! FUNCTION:
!       contains system dependent functions (e.g. command line io)
! CONTAINS:
!-----------------------------------------------------------------
!
! NAME:
!       sys_getcmdargs
! FUNCTION:
!       same as getarg (read command line option)
! SYNTAX:
!       call sys_getcmdargs(int, char)
! ON INPUT:
!       i_pos:  position of command line option         integer
! ON OUTPUT:
!       c_opt:  option read from command line           character string
! CALLS:
!
! COMMENTS:
!       i_pos=0 returns the command itself
!-----------------------------------------------------------------
!
! NAME:
!       sys_numcmdargs
! FUNCTION:
!       same as iargc (tell number of command line options)
! SYNTAX:
!       int= sys_numcmdargs
! ON INPUT:
!
! ON OUTPUT:
!       i_argnum: number of command line items          integer
! CALLS:
!
! COMMENTS:
!       i_argnum excludes the command itself, thus i_argnum=0
!       indicates that there are no command line options given.
!-----------------------------------------------------------------
!
! PUBLIC:
!       sys_getcmdargs, sys_numcmdargs
! COMMENTS:
!       this is for "standard" compilers, supporting iargc and getarg
! USES:
!
! LIBRARIES:
!
! REFERENCES:
!
! VERSION(S):
!       1. original version     j. behrens      1/97
!
!*****************************************************************
MODULE MISC_system
  PRIVATE
  PUBLIC :: sys_getcmdargs, sys_numcmdargs

  CONTAINS
!*****************************************************************
  SUBROUTINE sys_getcmdargs(i_pos, c_opt, i_len)

    IMPLICIT NONE

!---------- local declarations

    INTEGER, INTENT(in)             :: i_pos
    CHARACTER (len=80), INTENT(out) :: c_opt
    INTEGER, OPTIONAL               :: i_len
    EXTERNAL GETARG

!         CALL GETARG(i_pos, c_opt)
    CALL GET_COMMAND_ARGUMENT(i_pos, c_opt)

    RETURN
  END SUBROUTINE sys_getcmdargs

!*****************************************************************
  FUNCTION sys_numcmdargs() RESULT(i_argnum)

!---------- local declarations

    IMPLICIT NONE

    INTEGER :: i_argnum
    EXTERNAL   IARGC
    INTEGER :: iargc

!         i_argnum= IARGC()
    i_argnum=COMMAND_ARGUMENT_COUNT()

    RETURN
  END FUNCTION sys_numcmdargs
END MODULE MISC_system
