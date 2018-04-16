C+----------------------------------------------------------------------
C
      SUBROUTINE LOGAX (LUNOUT, N, X, AXLENG, ORIGIN, ENDPT, CYCLE,
     >   ERROR)
C
C
C     Description and usage:
C
C          This is a skeleton version of a (tentatively) planned log axis
C        scaling routine to be modeled after LINAX.  The present version just
C        does some error checking, with axis scaling by DISSPLA.  We assume
C        that the axis is IN-creasing.
C
C           Note that at present (Oct. 85), we do all the error handling at
C        this level, so the ERROR flag is just decoration for now.  (The
C        output values are set, but need not be acted upon.)
C
C
C     Parameters:
C
C        Name    Dimension  Type  I/O/S  Description
C        LUNOUT              I    I      Unit number for error reporting.
C
C        N                   I    I      Total number of data points.
C
C        X           N       R    I      Packed array of data for one axis.
C
C        AXLENG              R    I      Physical axis length, in inches
C                                        (assumed positive).
C
C        ORIGIN              R    I/O    Requested origin point (left or
C                                        lower end of the axis).  Must be
C                                        greater than zero.  Use flag value
C                                        999.0 to request auto-scaling.
C
C        ENDPT               R    I/O    Requested end point (right or
C                                        upper end of the axis).  Must be
C                                        greater than zero.  Use flag value
C                                        999.0 to request auto-scaling.
C
C        CYCLE               R      O    Log cycle length determined by
C                                        DISSPLA (inches/cycle).
C
C        ERROR               I      O    Trouble flag. >>> Not really used <<<
C
C                                           0:  Normal termination.
C
C                                          +1:  Data range is zero.
C                                          +2:  Input value of ORIGIN negative.
C                                          +3:  Input value of ENDPT negative.
C
C                                          -1:  Number of data points (N)
C                                               is non-positive.  (Bad
C                                               programmer input - abort.)
C                                          -2:  Axis length is not positive.
C                                               (Programming error - abort.)
C
C
C     Environment:  Digital VAX-11/780 VMS Version 4.1 (FORTRAN 77).
C                   ISSCO Graphics DISSPLA Version 9.2
C
C
C     Notes:
C
C        (1)  IMPLICIT NONE is non-standard.
C
C
C     Author:  Robert Kennelly, Informatics General Corporation.
C
C
C     Development history:
C
C         2 Aug. 1985    RAK    Initial design and coding.
C        26 Sep. 1985    RAK    Import error checking and DISSPLA scaling
C                               from old QUICK.  (We don't try to resolve
C                               conflicts - just use DISSPLA's results.)
C         9 Oct. 1985    RAK    Error-check N, ORIGIN, ENDPT, and force the
C                               X array to be > TINY (don't use ABS!).
C                               Report problems on LUNOUT directly - we
C                               handle everything here for now (i.e. the
C                               ERROR output is not really required yet).
C         8 Jan. 1985    RAK    Increased exponents of HUGE and TINY to 38,
C                               essentially the VAX limit for REALs.
C        13 Apr. 1987    RAK    Added internal loop-back to check for >1
C                               input error.  This is NOT the final answer.
C                               Hard STOP if AXLENG is not positive.
C                               Added some protection to the protection
C                               of DISSPLA's log axis scaling routine.
C        16 July 1987    RAK    Re-ordered the test for too-small CYCLE
C                               length so that we test the result from
C                               DISSPLA routine ALGPLT, not just the input.
C
C-----------------------------------------------------------------------


C     Declarations.
C     -------------

      IMPLICIT NONE

C     Constants.

      REAL
     >   FLAG, HUGE, ZERO, TEN, TINY
      PARAMETER
     >   (FLAG = 999.0,                        ! Enables automatic scaling
     >    HUGE = 1.0E+38,                      ! Overflow protection
     >    ZERO = 0.0E+0,
     >    TEN  = 1.0E+1,
     >    TINY = 1.0E-37)                      ! Underflow protection

C     Variables.

      INTEGER
     >   ERROR, I, LUNOUT, N
      REAL
     >   AXLENG, CYCLE, ENDPT, LOGMAX, LOGMIN, ORIGIN, POWER, X (N)


C     Execution.
C     ----------

    1 CONTINUE
      ERROR = 0

      IF (N .LE. 0)         ERROR = -1
      IF (AXLENG .LE. ZERO) ERROR = -2
      IF (ORIGIN .LE. ZERO) ERROR = +2
      IF (ENDPT  .LE. ZERO) ERROR = +3

C     In this (preliminary) version, just report/fix what we can within this
C     routine.  It might be better to separate analysis and error handling
C     eventually, as in LINAX.

C     Negative ERROR's are bugs - the routine needs to assume that N and
C     AXLENG have already been checked.

      IF (ERROR .EQ. -1) THEN
         STOP 'LOGAX:  Fatal array dimension error!'
      ELSE IF (ERROR .EQ. -2) THEN
         STOP 'LOGAX:  Fatal axis length error!'

C     Recoverable errors.

      ELSE IF (ERROR .EQ. +2) THEN
         WRITE (LUNOUT, 1000) 'origin', ORIGIN
         ORIGIN = FLAG
      ELSE IF (ERROR .EQ. +3) THEN
         WRITE (LUNOUT, 1000) 'end point', ENDPT
         ENDPT = FLAG
      END IF

C     Loop back for re-check in case of multiple errors.  (If error checking
C     is ever moved out, this would be accomplished by re-CALLing LOGAX.)

      IF (ERROR .NE. 0) GO TO 1

C     Ensure that the data to be plotted is greater than zero.

      DO 10, I = 1, N
         X (I) = MAX (X (I), TINY)
   10 CONTINUE

C     Determine raw data limits.

      IF (ORIGIN .EQ. FLAG) THEN
         ORIGIN = HUGE
         DO 20, I = 1, N
            ORIGIN = MIN (X (I), ORIGIN)
   20    CONTINUE
      END IF

      IF (ENDPT .EQ. FLAG) THEN
         ENDPT = TINY
         DO 30, I = 1, N
            ENDPT = MAX (X (I), ENDPT)
   30    CONTINUE
      END IF

C     Are the axis limits usable?  Must guarantee positivity, and create
C     a new ORIGIN if data range is zero.

      ORIGIN = MAX (ORIGIN, TINY)
      ENDPT  = MAX (ENDPT, TINY)

      IF (ORIGIN .EQ. ENDPT) THEN
         WRITE (LUNOUT, 1000) 'range', ZERO
         ERROR = +1
         ENDPT = TEN ** (LOG10 (ORIGIN) + 1)
      END IF

C     Let DISSPLA choose the endpoints and stepsize.  ORIGIN and ENDPT are
C     treated as extrema of the data to be accommodated.  For now, ignore any
C     conflicts between ALGPLT-computed and user-input ORIGIN.

      CALL ALGPLT (ORIGIN, ENDPT, AXLENG, ORIGIN, CYCLE)


C     Protect DISSPLA from a too-small cycle length.
C     ----------------------------------------------

C     This nonsense wouldn't be necessary if DISSPLA's log-axis drawing
C     routine were more robust!

      IF (CYCLE .LE. 0.4000) THEN

C        There is not enough room to include all of the data, so cut off the
C        plot at the smallest integral power of ten which will still fit.
C        Protect LOG10 by bounding POWER to avoid underflow.

         LOGMIN = LOG10 (ORIGIN)
         LOGMAX = LOG10 (ENDPT)

         IF (ORIGIN .LT. ENDPT) THEN

C           Normal (increasing) axis.

            POWER  = MAX (TEN ** (LOGMAX - (AXLENG / 0.4001)),
     >         LOG10 (TINY))
            LOGMIN = LOG10 (POWER)
            ORIGIN = TEN ** (INT (LOGMIN) + 1)
            ENDPT  = MAX (ENDPT, TEN * ORIGIN)
         ELSE

C           Inverted.  "MAX" here is really the "lower" limit.

            POWER  = MAX (TEN ** (LOGMIN - (AXLENG / 0.4001)),
     >         LOG10 (TINY))
            LOGMAX = LOG10 (POWER)
            ENDPT  = TEN ** (INT (LOGMAX) + 1)
            ORIGIN = MAX (ORIGIN, TEN * ENDPT)
         END IF

C        Try again for an acceptable CYCLE > .40 inches.

         CALL ALGPLT (ORIGIN, ENDPT, AXLENG, ORIGIN, CYCLE)
      END IF


C     Termination.
C     ------------

      RETURN


C     Formats.
C     --------

 1000 FORMAT (' LOGAX:  Warning - bad log axis ', A, ' = ', E10.3, '.'/
     >   9X, 'Automatic scaling will be used.'//)

      END
