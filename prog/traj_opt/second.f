C+----------------------------------------------------------------------
C
      SUBROUTINE SECOND (CPUSEC)
C
C PURPOSE:
C     Measures CPU time in same form as the CRAY utility.
C
C ARGUMENTS:
C    ARG      DIM  TYPE I/O/S DESCRIPTION
C    CPUSEC    -    R     O   CPU time used so far, in seconds.
C
C METHOD:
C     IRIS FORTRAN requires use of the C language's intrinsic clock
C     function, which returns CPU time used since the FIRST call to it.
C     However, because FORTRAN cannot access the intrinsic, we call a
C     C routine (iclock) that in turns calls clock.
C     Conversion to seconds is done here.
C
C KNOWN BUG:
C     The count of microseconds wraps around after about 36 minutes
C     for the original iclock form => the mclock form should be better.
C
C USAGE:
C        CALL SECOND (TIME1)
C        ::::::::::::::::::
C        CALL SECOND (TIME2)
C        TOTALT = TIME2 - TIME1
C
C ENVIRONMENT:  SGI IRIS, FORTRAN 77
C
C HISTORY:
C     08/31/82   Dan McKernan    VAX/VMS version.
C     05/30/90   David Saunders  IRIS version (iclock).
C     06/05/90   Dexter Hermstad IRIS version (continued).
C     03/21/97   D.Saunders      IRIS mclock form obtained from J.Reuther.
C     12/04/02       "           Intel's mclock seems to count millisecs.
C-----------------------------------------------------------------------

C     Arguments:

      REAL CPUSEC

C     C routine to call a system utility:

C*****INTEGER iclock
      INTEGER mclock

C     Execution:

C*****CPUSEC = REAL (iclock()) * 1.E-6
C*****CPUSEC = REAL (mclock()) * 0.01  ! SGI
      CPUSEC = REAL (mclock()) * 0.001 ! Intel (?)

      END
