C+----------------------------------------------------------------------
C
      SUBROUTINE LINAX (N, X, MXNSTP, ORIGIN, ENDPT, STEP, ERROR)
C
C     One-liner:  Linear plot axis scale selection, with error checking.
C     ----------
C
C     Description and usage:
C     ----------------------
C
C        Supervises scale selection and error checking for a linear plot
C     axis.  As commonly required by layout constraints, a maximum on the
C     number of scale divisions is imposed.  By default, coordinates
C     increase from ORIGIN to ENDPT, but this is not required.  Any of
C     the inputs ORIGIN, ENDPT and STEP may be specified by the user -
C     the rest will be automatically chosen as follows:  STEP is taken
C     from an array of "nice" values, using as many divisions as possible,
C     but preferring values which divide ORIGIN or ENDPT evenly if they
C     have already been supplied.  The computed values for ORIGIN and
C     ENDPT delimit the smallest range which encompasses the data and is
C     consistent with STEP.  This range will include at least some of the
C     data if possible, but there are cases where no sensible choice is
C     available.  (For example, when the data are all negative, and
C     ORIGIN = 0.0, STEP = 1.0, and ENDPT is free, there is nothing to
C     be done.  Such errors are flagged, as are programming blunders like
C     non-positive array bounds.  See Parameters and Notes for details.
C
C        LINAX, and co-routines EXPRESS, GETSTP and ROUND, were written
C     for use with QPLOT, developed in the Aerodynamics Division of Ames
C     Research Center.  They are based in (small) part on SCALE1 by
C     C. R. Lewart of Bell Labs (see Bibliography).
C
C     Parameters:
C     -----------
C
C     Name     Dimension  Type  I/O/S  Description
C     N                    I    I      Number of data points.
C
C     X           N        R    I      Plot data for one axis.
C
C     MXNSTP               I    I      Maximum number of axis divisions
C                                      permitted.
C
C     ORIGIN               R    I/O    Plot origin (left or lower corner
C                                      of plot on the page).  The axis
C                                      numbering is assumed to begin with
C                                      ORIGIN, which will be calculated
C                                      if a special "flag" value is input.
C                                      See Notes.
C
C     ENDPT                R    I/O    The "other end" of the plot axis,
C                                      usually right or upper corner. It
C                                      will be calculated if a special
C                                     "flag" value is input. See Notes.
C
C     STEP                 R    I/O    Axis labeling increment.  The first
C                                      division will be at ORIGIN + STEP.
C                                      It will be calculated if a special
C                                      "flag" value is input.  If the flag
C                                      is negative (-999), the axis will
C                                      be decreasing.  See Notes.
C
C     ERROR                I      O    Error return flag.
C                                        0:  Normal termination
C
C                                       +1:  Axis length is zero
C                                       +2:  Step length zero or too big
C                                       +3:  Step has wrong sign
C                                            (Bad user input - could just
C                                            reset to defaults and retry)
C
C                                       -1:  Number of data points N
C                                            is non-positive
C                                       -2:  Maximum number of intervals
C                                            MXNSTP is non-positive
C                                            (Bad programmer input -
C                                            abort)
C
C     External references:
C     --------------------
C
C     GETSTP   Axis stepsize selection from an array of candidates.
C     ROUND    Rounds a quotient "up" or "down" in magnitude.
C
C     Environment:  Digital VAX-11/780 VMS FORTRAN (FORTRAN 77)
C     ------------  Apple Macintosh, Absoft MacFORTRAN/020 v2.3
C
C     Notes:
C     ------
C
C     (1)  IMPLICIT NONE is non-standard.
C
C     (2)  A "flag" value input for ORIGIN, ENDPT, or STEP specifies that
C          these be calculated from the data (the default mode).  The
C          present version uses FLAG = 999.0, specified in a PARAMETER
C          statement below.  STEP = -FLAG has special significance: the
C          axis values will decrease from ORIGIN to ENDPT.
C
C     (3)  Error handling:  we trap anything impossible, and many of the
C          unlikely combinations, returning values of ERROR greater than
C          zero for "user" errors (where no sensible choice is available),
C          and less than zero for "programmer" errors, where no action
C          at all is possible - see Parameters section for details.  In
C          the case of user errors, it may suffice to call LINAX again
C          with FLAG in place of the offending quantities.  Two classes
C          of possible "errors" on the part of the user are permitted:
C
C               (a) if an input value of STEP is too small, the number
C                   of scale divisions may exceed MXNSTP
C
C               (b) if all three of the inputs are supplied by the user,
C                   the data may lie outside the range of the plot
C
C          These oddities were intentionally NOT guarded against, as
C          there may be situations where they are intended, and the
C          results are in any event easily recognized as unusual.
C
C          LINCHK, a companion routine available from the author, may
C          be used to interpret the error flag and take recovery action
C          if possible.  It was written for an application where it was
C          necessary to "correct" the user's input and try again.
C
C     (4)  In calculating nice values for both ORIGIN and ENDPT (the
C          last of the four main cases handled), the rounding is tricky:
C          we round ORIGIN "down" (toward zero) when
C
C               (a) the axis is increasing (i.e., STEP > 0)   -and-
C               (b) the initial estimate ORIGIN / STEP > 0
C
C          or when
C
C               (a') STEP < 0            -and-
C               (b') ORIGIN / STEP < 0
C
C          and "up" otherwise so that the data range is covered.  ENDPT
C          is handled analogously, mutatis mutandis. (It took me awhile,
C          too!)
C
C     (5)  DELTA, the rounding grace margin, should be larger than
C          machine epsilon.  It should be small enough that the effect
C          on the plot appearance is not significant.  (The size of the
C          error depends on the plot scale and size, but DELTA should
C          probably be less than the relative error due to hardware
C          precision, <precision, inches> / <axis length, inches>.)  A
C          typical value is 10**-6 < DELTA < 10**-4.  Keep in mind that
C          although this approach preserves the appearance of the data,
C          small "overhanging" quantities may seem to have been lost, e.g.,
C          the value -1.0E-6 won't show up on an axis which also includes
C          data of order +1.0E+2.
C
C     Bibliography:
C     -------------
C
C     (1)  LEWART, C. R.  Algorithms SCALE1, SCALE2, and SCALE3 for
C             Determination of Scales on Computer Generated Plots.
C             CACM, Vol. 16, No. 10, Oct. 1973 (Algorithm 463).
C
C     Author:  Robert Kennelly, Informatics General Corporation
C     -------
C
C     Development history:
C     --------------------
C
C     18 Apr. 1985    RAK    Initial design and coding.
C     15 July 1985    RAK    Use a variable for FIXED (STEP).  Modify test
C                            for when ENDPT = XMIN (or ORIGIN = XMAX) to
C                            take care of STEP < 0 case.  Eliminated
C                            use of FLAG to signal GETSTP that TARGET is
C                            not used (0.0 suffices now).
C     16 July 1985    RAK    Use variables for FIXED (ORIGIN) and (ENDPT).
C     18 July 1985    RAK    When ORIGIN = 0.0, and ENDPT fixed as well,
C                            use ENDPT as step selection target.
C     19 Aug. 1985    RAK    Added 2.5 to NICE array.  Use BIGNEG and
C                            BIGPOS in place of +/- HUGE constant.
C                            Removed "dead" variables MULT and ORDER.
C     21 Aug. 1985    RAK    Added ERROR = +3 output to allow recovery
C                            from a step with wrong sign.
C     23 Aug. 1985    RAK    Must iterate to be sure the final number
C                            of steps does not exceed MXNSTP when
C                            both endpoints are free.
C     15 Oct. 1985    RAK    STEP = -FLAG means decreasing axis.
C     15 Apr. 1987    RAK    Added special handling for MXNSTP = 1 case.
C     24 Apr. 1987    RAK    Changed some multiplies to SIGN comparision
C                            to reduce possibility of overflow.
C      5 Oct. 1987    RAK    Modify search for min and max to eliminate
C                            need for BIGNEG and BIGPOS. Reduce DELTA
C                            to 1.0E-6.
C
C-----------------------------------------------------------------------

C     Declarations.
C     -------------

      IMPLICIT NONE

C     Constants.

      INTEGER
     >   MXNICE
      PARAMETER
     >   (MXNICE = 5)

      REAL
     >   DELTA, FLAG, ZERO, ONE, HALF
      PARAMETER
     >   (DELTA  = 1.0E-6,
     >    FLAG   = 999.0E+0,
     >    ZERO   = 0.0E+0,
     >    ONE    = 1.0E+0,
     >    HALF   = 0.5E+0)

C     Arguments.

      INTEGER
     >   ERROR, MXNSTP, N
      REAL
     >   ENDPT, ORIGIN, STEP, X (N)

C     Local variables.

      LOGICAL
     >   FIXEND, FIXORG, FIXSTP, FLIP, VALID
      INTEGER
     >   I
      REAL
     >   FLOOR, NICE (MXNICE), RAW, REALMX, ROUND, TARGET, UPDOWN,
     >   XMAX, XMIN

C     Storage.

      DATA
     >   NICE /1.0E+0, 2.0E+0, 2.5E+0, 5.0E+0, 10.0E+0/

C     Execution.
C     ----------

C     Check the input parameters for nonsense - note that at most one
C     error will be flagged.

      ERROR = 0

      IF (N .LE. 0)       ERROR = -1
      IF (MXNSTP .LE. 0)  ERROR = -2
      IF (STEP .EQ. ZERO) ERROR = +2

      IF (ERROR .NE. 0)  GO TO 990

C     Initialization.

      REALMX = REAL (MXNSTP)

      FLIP = (STEP .EQ. -FLAG)
      IF (FLIP) STEP = ABS (STEP)

      FIXORG = (ORIGIN .NE. FLAG)
      FIXEND = (ENDPT .NE. FLAG)
      FIXSTP = (STEP .NE. FLAG)

C     Select a case (one of four).
C     ----------------------------

      IF (FIXORG .AND. FIXEND) THEN

C        Both axis limits are fixed.
C        ---------------------------

         RAW = ENDPT - ORIGIN
         IF (RAW .EQ. ZERO) THEN
            ERROR = +1
            GO TO 990
         END IF

         IF (FIXSTP) THEN

C           Error-check the supplied step.  Note that we do NOT check whether
C           this value will yield too many intervals (> MXNSTP).

            RAW = RAW / STEP
            IF (RAW .LT. ZERO) THEN

C              Input STEP must point in the right direction, i.e.,
C              sign (RAW) = sign (STEP), or RAW / STEP > 0.

               ERROR = +3
               GO TO 990
            ELSE IF (RAW .LT. ONE) THEN

C              Input STEP must not exceed the axis length.

               ERROR = +2
               GO TO 990
            END IF
         ELSE

C           Choose STEP.  (We prefer one which divides ORIGIN or ENDPT.)

            TARGET = ORIGIN
            IF (TARGET .EQ. ZERO) TARGET = ENDPT
            STEP = ZERO
            CALL GETSTP (MXNICE, NICE, REALMX, TARGET, RAW, STEP)
         END IF
      ELSE

C        We will be needing one or both of these extrema later.

         XMAX = X (1)
         XMIN = X (1)
         DO 10, I = 1,  N
            XMAX = MAX (X (I), XMAX)
            XMIN = MIN (X (I), XMIN)
   10    CONTINUE

         IF (FIXORG) THEN

C           Choose ENDPT.
C           -------------

C           Estimate ENDPT using maximum, but if we would miss the data
C           entirely or if STEP indicates a reversed axis, choose minimum.

            ENDPT = XMAX
            IF (ORIGIN .GE. XMAX .OR. STEP .LT. ZERO) ENDPT = XMIN

            IF (ENDPT .EQ. ORIGIN) THEN

C              The data are apparently stacked up above the ORIGIN, so
C              just produce a plot which will at least exhibit the data.

               IF (FIXSTP) THEN
                  ENDPT = ORIGIN + STEP
               ELSE
                  ENDPT = ORIGIN + ONE
               END IF
            END IF

            RAW = ENDPT - ORIGIN
            IF (FIXSTP) THEN

C              Is STEP pointing in the right direction?

               IF (SIGN (ONE, RAW) .NE. SIGN (ONE, STEP)) THEN
                  ENDPT = FLAG
                  ERROR = +3
                  GO TO 990
               END IF
            ELSE

C              Choose STEP.

               STEP = ZERO
               CALL GETSTP (MXNICE, NICE, REALMX, ORIGIN, RAW, STEP)
            END IF

C           Pick ENDPT so that there are an even number of intervals
C           of length STEP, beginning with ORIGIN.

            ENDPT = ORIGIN + ROUND (RAW, STEP, DELTA, +ONE) * STEP
         ELSE IF (FIXEND) THEN

C           Choose ORIGIN.
C           --------------

C           Estimate ORIGIN using minimum, but if we would miss the data
C           entirely or if STEP indicates a reversed axis, choose maximum.

            ORIGIN = XMIN
            IF (ENDPT .LE. XMIN .OR. STEP .LT. ZERO) ORIGIN = XMAX

            IF (ENDPT .EQ. ORIGIN) THEN

C              The data are apparently stacked up above the ENDPT, so
C              just produce a plot which will at least exhibit the data.

               IF (FIXSTP) THEN
                  ORIGIN = ENDPT - STEP
               ELSE
                  ORIGIN = ENDPT - ONE
               END IF
            END IF

            RAW = ENDPT - ORIGIN
            IF (FIXSTP) THEN

C              Is STEP pointing in the right direction?

               IF (SIGN (ONE, RAW) .NE. SIGN (ONE, STEP)) THEN
                  ORIGIN = FLAG
                  ERROR  = +3
                  GO TO 990
               END IF
            ELSE

C              Choose STEP.

               STEP = ZERO
               CALL GETSTP (MXNICE, NICE, REALMX, ENDPT, RAW, STEP)
            END IF

C           Pick ORIGIN so that there are an even number of intervals
C           of length STEP, ending with ENDPT.

            ORIGIN = ENDPT - ROUND (RAW, STEP, DELTA, +ONE) * STEP
         ELSE

C           Choose both ORIGIN and ENDPT.
C           -----------------------------

            IF (XMAX .EQ. XMIN) THEN

C              Just produce a plot which will at least exhibit the data.

               IF (FIXSTP) THEN
                  XMAX = XMAX + STEP
               ELSE
                  XMAX = XMIN + ONE
               END IF
            END IF

C           May have to loop until the number of steps is less than MXNSTP
C           since moving both ends can result in violation of the bound.
C           (Relevant only when XMIN, XMAX, and XSTEP are all free.)

            FLOOR = ZERO
   20       CONTINUE
               ORIGIN = XMIN
               ENDPT  = XMAX
               VALID  = .TRUE.

               IF (FIXSTP) THEN

C                 Turn the axis around if necessary.

                  IF (STEP .LT. ZERO) THEN
                     ORIGIN = XMAX
                     ENDPT  = XMIN
                  END IF
               ELSE

C                 Choose STEP.  Note that TARGET = ZERO is used to disable
C                 optional choice of a step which divides ORIGIN or ENDPT.

                  TARGET = ZERO
                  RAW    = ENDPT - ORIGIN
                  STEP   = FLOOR
                  CALL GETSTP (MXNICE, NICE, REALMX, TARGET, RAW, STEP)
               END IF

C              Compute tentative axis endpoints, and set the flag which will
C              force iteration if they are not acceptable.

               IF ((MXNSTP .EQ. 1) .AND.
     >             (.NOT.FIXSTP)   .AND.
     >             (SIGN (ONE, ORIGIN) .NE. SIGN (ONE, ENDPT))) THEN

C                 A special case (alas!) since rounding ORIGIN and ENDPT
C                 both away from zero always yields at least two intervals.
C                 Try arranging the endpoints symmetrically about zero.

                  ORIGIN = SIGN (HALF * STEP, ORIGIN)
                  ENDPT  = -ORIGIN

                  VALID  = (MAX (ABS (XMIN), ABS (XMAX)) .LE.
     >               ABS (ORIGIN))
               ELSE

C                 Pick ORIGIN and ENDPT so that each is an even multiple of
C                 STEP.  (See Notes for discussion of rounding direction.)

                  UPDOWN = -(SIGN (ONE, ORIGIN) * SIGN (ONE, STEP))
                  ORIGIN = ROUND (ORIGIN, STEP, DELTA, UPDOWN) * STEP

                  UPDOWN = SIGN (ONE, ENDPT) * SIGN (ONE, STEP)
                  ENDPT  = ROUND (ENDPT, STEP, DELTA, UPDOWN) * STEP

                  VALID = (FIXSTP) .OR.
     >               (NINT ((ABS (ENDPT - ORIGIN)) / STEP) .LE. MXNSTP)
               END IF

               IF (.NOT.VALID) THEN

C                 The chosen STEP was not large enough (too many steps or
C                 did not cover the data), so increase the minimum stepsize 
C                 and try again.

                  FLOOR = STEP
                  GO TO 20

               END IF

         END IF
      END IF

C     If input STEP was -FLAG, then we want a decreasing axis.

      IF (FLIP .AND. STEP .GT. ZERO) THEN
         XMIN   = MIN (ORIGIN, ENDPT)
         ORIGIN = MAX (ORIGIN, ENDPT)
         ENDPT  = XMIN
         STEP   = -STEP
      END IF

C     Termination.
C     ------------

  990 CONTINUE
      RETURN
      END
C+----------------------------------------------------------------------
C
      SUBROUTINE EXPRESS (X, MULT, EXPON)
C
C     One-liner:  Find "scientific notation" multiplier and exponent.
C     ----------
C
C     Description and usage:
C     ----------------------
C
C        EXPRESS (carefully) breaks up a floating point quantity into
C     a multiplier and exponent, where the multiplier is guaranteed to
C     lie in the range 1.0 <= MULT < 10.0 despite the hazards of finite-
C     precision arithmetic.  It was written to protect GETSTP, an axis
C     scaling routine whose iterative alogorithm requires the multiplier
C     to lie within a known range.
C
C     Parameters:
C     -----------
C
C     Name      Dimension  Type  I/O/S  Description
C     X                     R    I      Floating point number to be
C                                       re-expressed in multiplier
C                                       and exponent notation.
C
C     MULT                  R      O    Multiplier in the range
C                                       [1.0, 10.0); may share
C                                       storage with X.
C
C     EXPON                 I      O    Exponent such that
C                                       X = MULT * (10 ** EXPON).
C
C
C     Environment:  Digital VAX-11/780, VMS FORTRAN (FORTRAN 77)
C     ------------  Apple Macintosh, Absoft MacFORTRAN/020 v2.3
C
C     Notes:
C     ------
C
C        (1)  IMPLICIT NONE and eight character variable names are
C             non-standard.
C
C     Author:  Robert Kennelly, Sterling Federal Systems
C     -------
C
C     Development history:
C     --------------------
C
C      5 Oct. 1987    RAK    Original design and coding.
C
C-----------------------------------------------------------------------

C     Declarations.
C     -------------

      IMPLICIT NONE

C     Constants.

      REAL
     >   ZERO, ONE, TEN
      PARAMETER
     >  (ZERO = 0.0E+0,
     >   ONE  = 1.0E+0,
     >   TEN  = 1.0E+1)

C     Arguments.

      INTEGER
     >   EXPON
      REAL
     >   MULT, X

C     Execution.
C     ----------

      IF (X .EQ. ZERO) THEN

C        Go home early!

         MULT  = ZERO
         EXPON = 0
         GO TO 990
      END IF

C     EXPON is the exponent of the largest power of TEN less than or
C     equal to X. In "scientific" notation:
C
C        X = MULT * TEN ** EXPON, where 1.0 <= MULT < 10.0
C
C     Note that we have to fudge for small X since INT rounds up toward
C     zero when LOG is negative (ABS (X) < 1.0), and we want to round down.

      EXPON = INT (LOG10 (ABS (X)))
      IF (ABS (X) .LT. ONE) EXPON = EXPON - 1

C     Map X into the interval [1.0, 10.0), with some insurance.

      MULT = X * TEN ** (-EXPON)

      IF (MULT .GE. TEN) THEN
         MULT  = MULT / TEN
         EXPON = EXPON + 1
      ELSE IF (MULT .LT. ONE) THEN
         MULT  = MULT * TEN
         EXPON = EXPON - 1
      END IF

C     Termination.
C     ------------

  990 CONTINUE
      RETURN
      END
C+----------------------------------------------------------------------
C
      SUBROUTINE GETSTP (MXNICE, NICE, REALMX, TARGET, RAW, STEP)
C
C     One-liner:  Axis stepsize selection from an array of candidates.
C     ----------
C
C     Description and usage:
C     ----------------------
C
C        Computes a "nice" stepsize for linear plot axes, given an
C     estimate of the axis length, a bound on the number of intervals,
C     and an array of candidates in the range [1.0, 10.0].  The smallest
C     acceptable step is used, i.e. the largest feasible number of axis
C     divisions.  If requested, an attempt is made to choose a step which
C     divides a specifed "target" evenly.  Thus if one end of the aixs is
C     known, the step may be tailored to it (if possible) so that zero
C     will fall on a scale division.
C
C        The value of STEP must be preset to a lower bound on the final
C     stepsize.  This bound will normally be zero, but may be used to force
C     selection of a larger step if a previous call resulted in too many
C     axis intervals after the endpoints were rounded off.
C
C        GETSTP was written to be called by LINAX, but can also be used
C     by itself for simple plot setup.  Error checking is expected to
C     have been performed by the calling routine (see LINAX, for example).
C
C     Parameters:
C     -----------
C
C     Name     Dimension  Type  I/O/S  Description
C     MXNICE               I    I      Number of candidate interval
C                                      lengths, must be >= 1.
C
C     NICE      MXNICE     R    I      Array of interval lengths, in
C                                      the range [1.0, 10.0].  Note that
C                                      10.0 must be included.
C
C     REALMX               R    I      Maximum number of steps permitted.
C                                      Expected to be a whole number,
C                                      greater than zero.
C
C     TARGET               R    I      STEP should divide TARGET evenly
C                                      if possible.  Note that this step
C                                      is skipped (perhaps because not
C                                      required) if 0.0 is input.
C
C     RAW                  R    I      Lower bound on data range.  Must
C                                      be non-zero.
C
C     STEP                 R    I/O    On entry, a lower bound on the
C                                      desired stepsize (final STEP will
C                                      be strictly greater in magnitude
C                                      than the initial value).
C
C                                      The output is a "nice" interval
C                                      scaled so that the entire axis can
C                                      be covered by no more than REALMX
C                                      divisions, and larger in magnitude
C                                      than input STEP.
C
C     External references:
C     --------------------
C
C     EXPRESS  Find "scientific notation" multiplier and exponent.
C
C     Environment:  Digital VAX-11/780 VMS FORTRAN (FORTRAN 77)
C     ------------  Apple Macintosh, Absoft MacFORTRAN/020 v2.3
C
C     Notes:
C     ------
C
C     (1)  IMPLICIT NONE is non-standard.
C
C     (2)  The array of candidates must include 10.0, with additional
C          intermediate points up to the user - a typical array might
C          include 1.0, 2.0, and 5.0 as well.
C
C     (3)  EPS is a machine-dependent constant:  it should be set a
C          little larger than the precision level of floating point
C          calculations.
C
C     (4)  If it is not required that STEP divide TARGET evenly, use zero
C          for TARGET.  In this case, the first acceptable element of
C          NICE will be used.
C
C     (5)  No error checking is performed on REALMX and RAW, which must
C          be non-zero.  This utility is intended for step length
C          calculation after higher levels take care of error detection
C          and recovery.
C
C     Author:  Robert Kennelly, Informatics General Corporation
C     -------
C
C     Development history:
C     --------------------
C
C     23 May  1985    RAK    Initial design and coding.
C      8 July 1985    RAK    Changed test in 10 loop: accept NICE(I)
C                            if equal to STEP.
C     15 July 1985    RAK    Have to pre-scale TARGET to determine if
C                            STEP can divide TARGET evenly.  Test for
C                            TARGET = 0.0 to skip this step (eliminate
C                            use of FLAG).
C     16 July 1985    RAK    Use SCALE to save arithmetic result.
C     23 Aug. 1985    RAK    Added provision for input STEP used as
C                            lower bound (helps LINAX iterate).
C     10 July 1987    RAK    Must enter 10-loop with I = 0 (not 1),
C                            else STEP = NICE (1) is not recognized.
C      6 Oct. 1987    RAK    To avoid trouble due to limited precision,
C                            routine EXPRESS now guarantees that STEP
C                            lies in [1.0, 10.0).  Meaning of SCALE
C                            reversed, plus other minor cleanup.
C
C-----------------------------------------------------------------------

C     Declarations.
C     -------------

      IMPLICIT NONE

C     Constants.

      REAL
     >   EPS, ZERO, ONE, TEN
      PARAMETER
     >   (EPS   = 1.0E-6,
     >    ZERO  = 0.0E+0,
     >    ONE   = 1.0E+0,
     >    TEN   = 1.0E+1)

C     Arguments.

      INTEGER
     >   MXNICE
      REAL
     >   NICE (MXNICE), RAW, REALMX, STEP, TARGET

C     Local variables.

      INTEGER
     >   I, J
      REAL
     >   SCALE, TEST

C     Execution.
C     ----------

C     Lower bound for STEP is based on division of the raw interval length
C     into the maximum number of pieces, or slightly larger than the input
C     value (if non-zero).

      STEP = MAX (ABS (RAW / REALMX), ABS (STEP * (ONE + EPS)))

C     Re-express STEP as a multiplier in the range [1.0, 10.0), with
C     scratch variable I used for the corresponding exponent.

      CALL EXPRESS (STEP, STEP, I)
      SCALE = TEN ** I

C     Find the smallest element of NICE >= the raw step.

      I = 0
   10 CONTINUE
         I = I + 1
         IF (NICE (I) .LT. STEP .AND. I .LT. MXNICE) GO TO 10


C     Look for a larger element of NICE which divides TARGET evenly.

      IF (TARGET .NE. ZERO) THEN

C        The target is pre-scaled for comparison with NICE (cheaper than
C        scaling each element of NICE).

         TEST = TARGET / SCALE
         J    = I - 1
   20    CONTINUE
            IF (J .LT. MXNICE) THEN
               J = J + 1

C              Is the scaled target (nearly) a multiple of NICE?

               IF (ABS (MOD (TEST, NICE (J))) .GT. EPS) THEN
                  GO TO 20

               ELSE

C                 Found one!  Save the result and drop through.

                  I = J
               END IF
            END IF

      END IF

C     Map the best STEP back to the real world.

      STEP = SIGN (NICE (I), RAW) * SCALE

C     Termination.
C     ------------

      RETURN
      END
C+----------------------------------------------------------------------
C
      FUNCTION ROUND (NUM, DENOM, DELTA, UPDOWN)
C
C     One-liner:  Rounds a quotient "up" or "down" in magnitude.
C     ----------
C
C     Description and usage:
C     ----------------------
C
C        Permits modular treatment of rounding a quotient to an integer.
C     This single routine handles symmetric rounding either "up" or
C     "down" in magnitude, with provision for a small grace margin to
C     account for roundoff error or for plotting devices with limited
C     precision.  Input and output is in the form of REAL variables to
C     facilitate floating point computations.
C
C        ROUND was written for use with LINAX, a linear plot axis setup
C     routine, but may find use in other applications. It should improve
C     readability over the usual truncate-and-adjust method in cases
C     where it is necessary to handle different signs.
C
C     Parameters:
C     -----------
C
C     Name     Dimension  Type  I/O/S  Description
C     ROUND                R      O    Quotient NUM / DENOM, rounded as
C                                      specified by flag UPDOWN. (FUNCTION
C                                      value.)
C
C     NUM                  R    I      Numerator of quotient to be rounded.
C
C     DENOM                R    I      Denominator of quotient. Assumed
C                                      to be non-zero.
C
C     DELTA                R    I      Absolute rounding tolerance - if
C                                      the quotient NUM / DENOM is within
C                                      DELTA of an integer, then that
C                                      value is returned even if it rounds
C                                      in the "wrong" direction.
C
C     UPDOWN               R    I      Specifies which direction to round
C                                      the magnitude of the quotient:
C                                       > 0:  "up" (away from zero)
C                                       < 0:  "down" (toward zero)
C                                       = 0:  plain truncation
C
C     Environment:  Digital VAX-11/780 VMS FORTRAN (FORTRAN 77)
C     ------------  Apple Macintosh, Absoft MacFORTRAN/020 v2.3
C
C     Notes:
C     ------
C
C     (1)  IMPLICIT NONE is non-standard.
C
C     (2)  If UPDOWN is zero, the truncated result of (NUM / DENOM) is
C          returned.   This may be useful in some cases, but a direct
C          call to AINT is generally more appropriate.
C
C     (3)  We don't check whether DENOM is zero (a tentative design
C          choice).  In the graphics applications for which ROUND was
C          originally designed, such error checking is better done at
C          a higher level.
C
C     (4)  DELTA may be zero if no grace margin is desired, i.e. strict
C          rounding up or down in magnitude.  It should normally be at
C          least as big as machine epsilon.
C
C     Author:  Robert Kennelly, Informatics General Corporation
C     -------
C
C     Development history:
C     --------------------
C
C      3 June 1985    RAK    Initial design and coding.
C     28 June 1985    RAK    SIGN of rounding correction taken from
C                            RAW, not ROUND.
C      7 Oct. 1987    RAK    Corrected header - ROUND is merely the
C                            truncated quotient when UPDOWN is zero, not
C                            the nearest integer. Other minor prettifying.
C
C-----------------------------------------------------------------------

C     Declarations.
C     -------------

      IMPLICIT NONE

C     Constants.

      REAL
     >   ZERO, ONE
      PARAMETER
     >   (ZERO  = 0.0E+0,
     >    ONE   = 1.0E+0)

C     Arguments.

      REAL
     >   DELTA, DENOM, NUM, RAW, ROUND, UPDOWN

C     Local variables.

      REAL
     >   TEST

C     Execution.
C     ----------

      RAW   = NUM / DENOM
      ROUND = AINT (RAW)

C     Round up or down (unless flag UPDOWN is zero).  If we're very close
C     (2.0 + epsilon, say, when rounding up), chose NEAREST integer (2.0).

      IF (UPDOWN .NE. ZERO) THEN
         TEST = DELTA
         IF (UPDOWN .LT. ZERO) TEST = ONE - DELTA

C        The change in magnitude resulting from the previous "downward"
C        round (truncation) is used to check whether we need to either
C
C           (a) add an additional amount (yielding an "upward" round), or
C           (b) correct a truncated value (if a value rounded "down" was
C               close to an integer).

         IF (ABS (ROUND - RAW) .GT. TEST)
     >      ROUND = ROUND + SIGN (ONE, RAW)
      END IF

C     Termination.
C     ------------

      RETURN
      END
