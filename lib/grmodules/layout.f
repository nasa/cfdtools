C+----------------------------------------------------------------------
C
      SUBROUTINE LAYOUT (LUNOUT, N, X, Y, XMIN, XMAX, XSTEP, YMIN, YMAX,
     >   YSTEP, HEIGHT, WIDTH, DEFHGT, DEFWID, PLOT)
C
C     One-liner:  Detailed plot scaling/sizing, with conflict resolution
C     ----------
C
C     Description and usage:
C     ----------------------
C
C        Detailed data plot layout is handled in this one module in order
C     to hide the messy details from the rest of an application.  LAYOUT
C     was originally written for use by QPLOT, which is based on the DISSPLA
C     proprietary graphics package.  This routine is more powerful (for
C     linear axes, at least) than the setup utilities provided by DISSPLA,
C     so it has been re-written for general use.  Its functions include:
C
C        (1) defaulting,
C        (2) input decoding,
C        (3) conflict resolution, and
C        (4) error checking/correction
C
C     of the axis scaling, physical size, and plot type.  Dictionary
C     look-ups permit abbreviations and synonyms for the string type
C     inputs.  A fairly powerful linear axis scaling routine will handle
C     nearly any combination of inputs, with (optional) error reporting
C     on the few cases where no sensible interpretation of the input data
C     is possible.  Some error checking is applied to logarithmic axes,
C     while polar plots are handled in a rudimentary fashion simply by
C     passing the data straight through.  (These limitations on log and
C     polar plots are a feature which should be improved if LAYOUT is to
C     be used with graphics packages other than DISSPLA.)
C
C        Arguments for which the application expects LAYOUT to determine
C     values should be FLAGged as 999. on input.
C
C        If your data is not in single X and Y arrays, or if you want to
C     suggest values for minimum and maximum, pass in arrays X and Y of
C     length 2 containing the desired extremes.
C
C        See the Notes below for some restrictions on LAYOUT's generality.
C     The QPLOT header contains detailed information on that program's
C     control inputs.  Dictionary routine LOOKUP is part of a small set
C     of FORTRAN 77 input utilities available from the author.
C
C     Arguments:
C     ----------
C
C     Name    Dimension  Type  I/O/S  Description
C     LUNOUT              I    I      Logical unit number for diagnostic
C                                     output.  Use < 0 to suppress output.
C
C     N                   I    I      Number of (X,Y) pairs to be plotted.
C
C     X         N         R    I      Packed array of abcissas.
C
C     Y         N         R    I      Packed array of ordinates.
C
C     XMIN                R    I/O    Left-hand endpoint of the horizontal
C                                     axis, in plot units. Use XMIN = FLAG
C                                     to indicate that a value is to be
C                                     calculated.
C
C     XMAX                R    I/O    Right-hand endpoint of the horizontal
C                                     plot axis.  Use FLAG to indicate
C                                     that a value is to be calculated.
C
C     XSTEP               R    I/O    Length of axis divisions on the
C                                     horizontal axis.  Use FLAG to indicate
C                                     that a value is to be calculated.
C                                     Negative FLAG means auto-scale with
C                                     decreasing axis.
C
C     YMIN                R    I/O    Left-hand endpoint of the vertical
C                                     plot axis.  Use FLAG to indicate
C                                     that a value is to be calculated.
C
C     YMAX                R    I/O    Right-hand endpoint of the vertical
C                                     plot axis.  Use FLAG to indicate
C                                     that a value is to be calculated.
C
C     YSTEP               R    I/O    Length of axis divisions on the
C                                     vertical axis.  Use FLAG to indicate
C                                     that a value is to be calculated.
C                                     Negative FLAG means auto-scale with
C                                     decreasing axis.
C
C     HEIGHT              R    I/O    User-requested length of the vertical
C                                     axis, in physical units.  This value
C                                     will be overridden if too small (see
C                                     parameter DENSITY = max. no. of scale
C                                     divisions per unit), or in case of
C                                     conflicts with equally-scaled plot
C                                     axes. If HEIGHT = FLAG, then the
C                                     default value is used (see DEFHGT).
C
C     WIDTH               R    I/O    User-requested length of the horizontal
C                                     axis, in physical units.  See HEIGHT.
C
C     DEFHGT              R    I      Default length of the vertical axis,
C                                     in physical units.
C
C     DEFWID              R    I      Default length of the horizontal axis,
C                                     in physical units.
C
C     PLOT                C    I/O    Plot type identifier.  After checking
C                                     dictionary, converted to one of the
C                                     standard forms:  'LINEAR', 'LOGLOG',
C                                     'LOGX', 'LOGY', 'POLAR', and 'SCALE',
C                                     whose meanings should be evident. The
C                                     default is 'LINEAR'.  Note that 'SCALE'
C                                     is changed to 'LINEAR' after the axes
C                                     have been set up successfully.
C
C     External references:
C     --------------------
C
C     FIXWINDO  Determines X range left by Y windowing, and vice versa
C     LINAX     Chooses (and error checks) axis limits and step sizes
C     LINCHK    Analyzes errors reported by LINAX, "corrects" some of them
C     LOGAX     Skeleton routine only, for now (uses DISSPLA)
C     LOOKUP    Dictionary lookup, with abbreviations and synonyms
C     UPCASE    Used to ensure PLOT is in upper case
C
C     Environment:  Digital VAX-11/780 VMS FORTRAN (FORTRAN 77)
C     ------------
C
C     Notes:
C     ------
C
C     (1)  IMPLICIT NONE, 8-character symbolic names, and use of "!" for
C          comments is not standard FORTRAN 77.
C
C     (2)  The SCALE option cannot, in general, respect all user choices
C          of X- and Y-axis length.  We have chosen to base the plot
C          layout on the physical length input for that axis whose length
C          is input, or the Y-axis in case of conflict.  There is no guar-
C          antee that the resulting layout will fit on the eventual output
C          medium!
C
C     (3)  As above, a user-supplied value for XSTEP will be ignored if it
C          conflicts with YSTEP in those cases (the majority) where the
C          Y-axis information takes precedence.
C
C     (4)  Only rudimentary range checking is applied to HEIGHT and WIDTH:
C          they must be at least one minimum-grid-space.  See parameter
C          DENSITY, below.  It may be desirable to add upper bounds as well.
C
C     (5)  For log axes, we re-use the STEP variables for passing DISSPLA's
C          CYCLE (inches/cycle) quantities.  This may be inappropriate (or
C          merely meaningless) for use with other graphics packages.
C
C     (6)  The FLAG value used here for defaulting (+/- 999) is arbitrary,
C          but must be consistent with the other routines.  The idea is NOT
C          to force the user to use FLAG to obtain defaults, but as a means
C          of communication between routines:  the calling program should
C          just preset the variable to FLAG, then read in the user's value
C          and pass it to LAYOUT.
C
C     (7)  For SCALE plots, we use the same STEP on both axes if possible,
C          but if the derived axis' STEP is too big we have to re-try.  The
C          result will be equally scaled, but may not look it.
C
C     (8)  The physical units used are expected to be inches, as reflected
C          in choice of DENSITY (equivalent to 3/4" inch per scale division)
C          and GRACE (calculated scaled axis lengths may exceed the default
C          values by 1").  These two constants will need to be changed if
C          another system of measurement is to be used.
C
C     History:
C     --------
C
C     11 Apr. 1985    RAK     Initial design and coding.
C      6 Aug. 1985    RAK     SCALE option adapted from EQUSTP/GRNICE
C                             approach used in earlier QPLOT.
C     12 Aug. 1985    RAK     Added LINCHK to decode error messages from
C                             LINAX (modular, saves space).  Moved
C                             assignment of default axis lengths down
C                             from QPLOT level.
C     21 Aug. 1985    RAK     Tie minimum HEIGHT and WIDTH to DENSITY.
C                             Print informational message for SCALE plots
C                             telling which axis controls the layout.
C                             Round down (INT) to ensure that the maximum
C                             grid density does not exceed DENSITY.  Added
C                             ORIENTation parameter.  Default HEIGHT and
C                             WIDTH depend on ORIENT.
C      9 Sep. 1985    RAK     Decreased DENSITY.  Shuffled diagnostic
C                             WRITEs.  Permit XSTEP to govern Y-axis in
C                             SCALE plots if YSTEP is not input, but
C                             retain preference for positive steps.
C     27 Sep. 1985    RAK     Dropped 'CONVERGENCE' plot option (just
C                             point back to 'LOGY').  Re-set 'SCALE' plot
C                             type to 'LINEAR' when through with setup.
C                             'POLAR' plots are passed through with no
C                             checking, and only radial axis is set up.
C      9 Oct. 1985    RAK     Nasty special case:  if MIN, MAX, and STEP
C                             are all specified for SCALE plot, but
C                             STEP is too big, then default STEP and try
C                             again (since we now know that MXNSTP can
C                             be computed).  Protect all MXNSTPs with
C                             MAX (*, 1), esp. for POLAR.
C     16 Oct. 1985    RAK     Try to retain the newly-significant sign
C                             when a STEP input equals -FLAG (indicating
C                             a decreasing axis).  Orientation dictionary
C                             was disoriented (LANDSCAPE out of place).
C     12 Mar. 1987    RAK     Renamed LAYOUT and (decreed to be general-
C                             purpose) moved to graphics library. Changed
C                             orientation terminology to standard
C                             PORTRAIT (default) and LANDSCAPE.
C     24 Mar. 1987  R.Langhi  Removed orientation option, and replaced
C                             it with default axis lengths input by the
C                             calling routine.
C     23 Apr. 1987    RAK     Revised error message spacing (new name).
C                             POLAR step lengths are now rounded values.
C     06 Oct. 1991 D.Saunders Introduced FIXWINDO to solve long-standing
C                             difficulty caused by partially-specified
C                             plot windows, where the global data extremes
C                             are not always good choices - an indicated
C                             window in X, say, may exclude quite a range
C                             of Y, and vice versa.
C
C     Author:  Robert Kennelly, NASA Ames Research Center, Mt. View, CA.
C     -------
C-----------------------------------------------------------------------

C     Declarations.
C     -------------

      IMPLICIT NONE

C     Constants.

      LOGICAL
     >   ORDER
      PARAMETER
     >  (ORDER = .TRUE.)

      INTEGER
     >   LENDIC, NPLDIC
      PARAMETER
     >  (LENDIC = 14,                 ! Length of plot dictionary entries
     >   NPLDIC = 10)                 ! Number of dictionary entries

      REAL
     >   DENSITY, FLAG, GRACE, MINLEN
      PARAMETER
     >  (DENSITY = 4.0 / 3.0,         ! Max. number of steps/inch
     >   FLAG    = 999.0,             ! Special value for indicating defaults
     >   GRACE   = 1.0,               ! Default SCALE plot "expansion" limit
     >   MINLEN  = 1.0 / DENSITY)     ! Min. length of an axis step

C     Arguments.

      INTEGER
     >   LUNOUT, N
      REAL
     >   DEFHGT, DEFWID, HEIGHT, WIDTH, X (N), XMAX, XMIN, XSTEP, Y (N),
     >   YMAX, YMIN, YSTEP
      CHARACTER
     >   PLOT * (*)

C     Local variables.

      LOGICAL
     >   CHECK
      INTEGER
     >   I, ITEM, MXNSTP, OOPS
      REAL
     >   EDGES (4), URXMAX, URXMIN, URXSTP, URYMAX, URYMIN, URYSTP
      CHARACTER
     >   PLDIC (NPLDIC) * (LENDIC)

C     Storage.

      DATA
     >  (PLDIC (I), I = 1, NPLDIC)
     >     /'DEFAULT=LINEAR',
     >      'LINEAR',
     >      'LOGLOG',
     >      'LOGX',
     >      'LOGY',
     >      'POLAR',
     >      'SCALE',
     >      'XLOG=LOGX',
     >      'XY=LINEAR',
     >      'YLOG=LOGY'/

C     Execution.
C     ----------

      HEIGHT = MAX (ABS (HEIGHT), MINLEN)
      WIDTH  = MAX (ABS (WIDTH),  MINLEN)

C     Fix up plot type flag.

      CALL UPCASE (PLOT)
      CALL LOOKUP (NPLDIC, PLDIC, ORDER, PLOT, ITEM)
      IF (ITEM .LE. 0) THEN
         IF (LUNOUT .GT. 0) WRITE (LUNOUT, 1000) ITEM
         PLOT = 'LINEAR'
      END IF

C     Save inputs so that the plot scale may be re-computed later.  (Prefix
C     "UR" means "original" in German, if that helps...)

      URXMIN = XMIN
      URXMAX = XMAX
      URXSTP = XSTEP

      URYMIN = YMIN
      URYMAX = YMAX
      URYSTP = YSTEP

C     Handle the problem of partially-specified windows: using global data
C     extrema for the unspecified parts of the plot window is not good
C     enough because the given clipping in X (say) may exclude some of the
C     Y range, and vice versa.

      CALL FIXWINDO (N, X, Y, URXMIN, URXMAX, URYMIN, URYMAX, FLAG,
     >               EDGES)

C     EDGES (1:2) can now replace X (*) in LINAX calls, and
C     EDGES (3:4) can now replace Y (*) similarly.


C     The 'SCALE' plot option can fail:  with default settings, a derived
C     X-axis can be too long, i.e., WIDTH much greater than DEFWID.  This
C     loop allows a second try, with the Y-axis derived from the X-axis.

   10 CONTINUE
         CHECK = .FALSE.

         IF (PLOT .EQ. 'SCALE') THEN

C           Find step size and axis lengths for equal axis scaling.
C           -------------------------------------------------------

C           Axis lengths and stepsizes are controlled by X-axis input only
C           if the user explicitly requests them, and does NOT explicitly
C           request a corresponding Y-axis value.

            IF (WIDTH .NE. FLAG .AND. HEIGHT .EQ. FLAG) THEN

C              Use X-axis length to determine required length of Y-axis.

               MXNSTP = MAX (INT (WIDTH * DENSITY), 1)
               IF (ABS (YSTEP) .NE. FLAG) XSTEP = ABS (YSTEP)
               IF (LUNOUT .GT. 0) WRITE (LUNOUT, 1010) 'X'
   20          CONTINUE
                  CALL LINAX (2, EDGES (1), MXNSTP, XMIN, XMAX, XSTEP,
     >                        OOPS)
                  IF (OOPS .NE. 0) THEN
                     CALL LINCHK (LUNOUT, OOPS, 'X', XMIN, XMAX, XSTEP)
                     GO TO 20

                  END IF

C              Since the two axes usually use the same step, specifying the
C              maximum number of steps for the derived axis may not be
C              required.  But if LINAX returns ERROR = +2, then the step
C              derived above is too big and we will have to try again, using
C              computed MXNSTP (YMAX and YMIN must have been fixed).

               MXNSTP = 999
               YSTEP  = SIGN (XSTEP, YSTEP)
   30          CONTINUE
                  CALL LINAX (2, EDGES (3), MXNSTP, YMIN, YMAX, YSTEP,
     >                        OOPS)
                  HEIGHT = WIDTH * ABS ((YMAX - YMIN) / (XMAX - XMIN))
                  IF (OOPS .NE. 0) THEN
                     CALL LINCHK (LUNOUT, OOPS, 'Y', YMIN, YMAX, YSTEP)
                     IF (OOPS .EQ. +2)
     >                  MXNSTP = MAX (INT (HEIGHT * DENSITY), 1)
                     GO TO 30

                  END IF

            ELSE

C              Use Y-axis length to determine required length of X-axis.

               IF (HEIGHT .EQ. FLAG) THEN
                  HEIGHT = DEFHGT
                  IF (WIDTH .EQ. FLAG) CHECK = .TRUE.
               END IF

               MXNSTP = MAX (INT (HEIGHT * DENSITY), 1)
               IF (ABS (YSTEP) .EQ. FLAG) YSTEP = SIGN (XSTEP, YSTEP)
               IF (LUNOUT .GT. 0) WRITE (LUNOUT, 1010) 'Y'
   40          CONTINUE
                  CALL LINAX (2, EDGES (3), MXNSTP, YMIN, YMAX, YSTEP,
     >                        OOPS)
                  IF (OOPS .NE. 0) THEN
                     CALL LINCHK (LUNOUT, OOPS, 'Y', YMIN, YMAX, YSTEP)
                     GO TO 40

                  END IF

C              See comment above on re-try's when a user-specified SCALE
C              plot must resort to different STEPs on X- and Y-axes.

               MXNSTP = 999
               XSTEP  = SIGN (YSTEP, XSTEP)
   50          CONTINUE
                  CALL LINAX (2, EDGES (1), MXNSTP, XMIN, XMAX, XSTEP,
     >                        OOPS)
                  WIDTH = HEIGHT * ABS ((XMAX - XMIN) / (YMAX - YMIN))
                  IF (OOPS .NE. 0) THEN
                     CALL LINCHK (LUNOUT, OOPS, 'X', XMIN, XMAX, XSTEP)
                     IF (OOPS .EQ. +2)
     >                  MXNSTP = MAX (INT (WIDTH * DENSITY), 1)
                     GO TO 50

                  END IF

            END IF         
         ELSE

C           Treat the X- and Y-axes individually (unequal scaling).
C           -------------------------------------------------------

C           Set axis sizes to default values if requested.

            IF (WIDTH  .EQ. FLAG) WIDTH  = DEFWID
            IF (HEIGHT .EQ. FLAG) HEIGHT = DEFHGT

            IF (PLOT .EQ. 'LINEAR' .OR. PLOT .EQ. 'LOGY') THEN

C              Linear X-axis.

               MXNSTP = MAX (INT (WIDTH * DENSITY), 1)
   60          CONTINUE
                  CALL LINAX (2, EDGES (1), MXNSTP, XMIN, XMAX, XSTEP,
     >                        OOPS)
                  IF (OOPS .NE. 0) THEN
                     CALL LINCHK (LUNOUT, OOPS, 'X', XMIN, XMAX, XSTEP)
                     GO TO 60

                  END IF

            ELSE IF (PLOT .EQ. 'LOGX' .OR. PLOT .EQ. 'LOGLOG') THEN

C              Logarithmic X-axis.

               CALL LOGAX (LUNOUT, 2, EDGES (1), WIDTH, XMIN, XMAX,
     >                     XSTEP, OOPS)
            END IF

            IF (PLOT .EQ. 'LINEAR' .OR. PLOT .EQ. 'LOGX') THEN

C              Linear Y-axis.

               MXNSTP = MAX (INT (HEIGHT * DENSITY), 1)
   70          CONTINUE
                  CALL LINAX (2, EDGES (3), MXNSTP, YMIN, YMAX, YSTEP,
     >                        OOPS)
                  IF (OOPS .NE. 0) THEN
                     CALL LINCHK (LUNOUT, OOPS, 'Y', YMIN, YMAX, YSTEP)
                     GO TO 70

                  END IF

            ELSE IF (PLOT .EQ. 'LOGY' .OR. PLOT .EQ. 'LOGLOG') THEN

C              Logarithmic Y-axis.

               CALL LOGAX (LUNOUT, 2, EDGES (3), HEIGHT, YMIN, YMAX,
     >            YSTEP, OOPS)

            ELSE IF (PLOT .EQ. 'POLAR') THEN

C              Polar diagram of X (THETA, in radians) vs. Y (RADIUS).

C              Required step parameter is units/inch, based on an axis
C              of length = (half of the minimum dimension).  THIS IS AN
C              ODDBALL CASE - DISSPLA always uses one step/inch for POLAR.
C              NOTE: since YMAX will be recomputed later by DISSPLA, the
C              final plot may not fill the available space exactly.

               MXNSTP = MAX (INT (0.50 * MIN (HEIGHT, WIDTH)), 1)

C              Default minimum radius is zero - and note that we can't
C              specify a decreasing R-axis using DISSPLA's POLAR routine.

               IF (YMIN .EQ. FLAG) YMIN = 0.0
               YMAX  = FLAG
               YSTEP = ABS (YSTEP)
   80          CONTINUE
                  CALL LINAX (2, EDGES (3), MXNSTP, YMIN, YMAX, YSTEP,
     >                           OOPS)
                  IF (OOPS .NE. 0) THEN
                     CALL LINCHK (LUNOUT, OOPS, 'R', YMIN, YMAX, YSTEP)
                     GO TO 80

                  END IF
            END IF
         END IF

         IF (CHECK .AND. WIDTH .GT. DEFWID + GRACE) THEN

C           A default 'SCALE' plot layout has failed - reset inputs,
C           specifying WIDTH explicitly so that the X-axis will dominate
C           on retry.  (Only one retry will be attempted.)

            HEIGHT = FLAG
            WIDTH  = DEFWID

            XMIN   = URXMIN
            XMAX   = URXMAX
            XSTEP  = URXSTP

            YMIN   = URYMIN
            YMAX   = URYMAX
            YSTEP  = URYSTP

            IF (LUNOUT .GT. 0) WRITE (LUNOUT, 1020)
            GO TO 10

         END IF

C     Once the axes have been set up, a SCALE plot is the same as LINEAR.

      IF (PLOT .EQ. 'SCALE') PLOT = 'LINEAR'

C     Termination.
C     ------------

      RETURN

C     Formats.
C     --------

 1000 FORMAT (' LAYOUT:  Warning - the plot type was improperly ',
     >   'specified.  (ITEM = ', I3, ')'/
     >   10X, 'The default will be used (LINEAR).')
 1010 FORMAT (' LAYOUT:  Axis lengths and scaling will be derived from ',
     >   'the ', A, '-axis information.')
 1020 FORMAT (' LAYOUT:  The derived X-axis length was much longer ',
     >   'than the default width.'/
     >   10X, 'The layout will be repeated with derived Y-axis length.')

      END
C+----------------------------------------------------------------------
C
      SUBROUTINE LINCHK (LUNOUT, ERROR, LABEL, AXMIN, AXMAX, AXSTEP)
C
C     One-liner:  Error checking/recovery for axis scaling routine LINAX
C     ----------
C
C     Description and usage:
C     ----------------------
C
C        Modularizes the error checking left out of LINAX.  Normal usage,
C     if error checking/correction is desired, is to call LINCHK after
C     LINAX if ERROR is non-zero, and then call LINAX again using the
C     patched axis input values (see LAYOUT for an example).  A hard STOP
C     is provided as a debugging aid if LINAX was called improperly and
C     no sensible recovery is possible.
C
C        LINCHK was designed as part of QPLOT, a general purpose plotting
C     package developed at NASA - Ames Research Center.
C
C     Arguments:
C     ----------
C
C     Name    Dimension  Type  I/O/S  Description
C     LUNOUT              I    I      Logical unit number for error
C                                     messages (line printer format).
C                                     Use < 0 to suppress output.
C
C     ERROR               I    I      Error code from LINAX to be
C                                     analyzed.
C
C     LABEL     *         C    I      Name of the axis being set up, e.g.
C                                     'X', for error messages.
C
C     AXMIN               R    I/O    Axis ORIGIN from LINAX.
C
C     AXMAX               R    I/O    Axis ENDPT from LINAX.
C
C     AXSTEP              R    I/O    Stepsize from LINAX.
C
C
C     Environment:  Digital VAX-11/780 VMS FORTRAN (FORTRAN 77)
C     ------------
C
C     Notes:
C     ------
C
C        (1)  IMPLICIT NONE is non-standard.
C
C
C     History:
C     --------
C
C        12 Aug. 1985    RAK    Initial design and coding.
C        14 Aug. 1985    RAK    For ERROR=+2 (bad AXSTEP), just retry
C                               with the sign reversed.
C        21 Aug. 1985    RAK    Had to distinguish between bad sign (+3)
C                               and bad length of AXSTEP (+2).  Reordered
C                               calling sequence.
C        24 Dec. 1985    RAK    Length of LABEL is now as-passed.  Try
C                               SP sign control editing in error message.
C        27 Mar. 1987    RAK    Negative LUNOUT suppresses messages. Use
C                               A, A in FORMAT to avoid concatenation.
C
C     Author:  Robert Kennelly, NASA Ames Research Center, Mt. View, CA.
C     -------
C-----------------------------------------------------------------------

C     Declarations.
C     -------------

      IMPLICIT NONE

C     Constants.

      REAL
     >   FLAG
      PARAMETER
     >  (FLAG = 999.0)

C     Arguments.

      INTEGER
     >   ERROR, LUNOUT
      REAL
     >   AXMIN, AXMAX, AXSTEP
      CHARACTER
     >   LABEL * (*)

C     Execution.
C     ----------

      IF (ERROR .EQ. +1) THEN

C        Axis length was zero!

         IF (LUNOUT .GT. 0) WRITE (LUNOUT, 1000) LABEL, '-axis', ERROR
         AXMIN = FLAG
         AXMAX = FLAG

      ELSE IF (ERROR .EQ. +2) THEN

C        Input AXSTEP was zero, or too big.

         IF (LUNOUT .GT. 0) WRITE (LUNOUT, 1000) LABEL, 'STEP', ERROR
         AXSTEP = FLAG

      ELSE IF (ERROR .EQ. +3) THEN

C        The step had the wrong sign - we'll "correct" it.

         IF (LUNOUT .GT. 0) WRITE (LUNOUT, 1000) LABEL, 'STEP', ERROR
         AXSTEP = -AXSTEP

      ELSE IF (ERROR .LT. 0) THEN

C        Programming error:  the number of data points or maximum number
C        of steps passed to LINAX was not positive.

         IF (LUNOUT .GT. 0) WRITE (LUNOUT, 1010) LABEL
         STOP 'LINCHK:  Fatal error!'

      END IF

C     Termination.
C     ------------

      RETURN

C     Formats.
C     --------

 1000 FORMAT (' LINCHK:  Warning - bad ', A, A, ' information.  ',
     >   '(ERROR = ', SP, I2, ')'/
     >   10X, 'Automatic scaling will be used.')
 1010 FORMAT (' LINCHK:  Error - bad call to LINAX for ', A, '-axis.')

      END
