C+------------------------------------------------------------------------------
C
      SUBROUTINE POLYLINE (N, X, Y, LINE, THICK, SYMBOL, COLOR, IER)
C
C One-liner:  
C ----------
C     Connect-the-dots with given pattern, thickness, symbol, color (level 3)
C
C Description and usage:
C ----------------------
C
C        POLYLINE is a high-level "curve" drawing routine with a generic
C     programmer interface.  It isolates at least some of the details of
C     the specific graphics library in use.  (Plot set-up functions still
C     have to be translated elsewhere when switching to another library.)
C
C        This version is for the DISSPLA package (originally developed
C     for the PLOTCL application and since adapted for QPLOT and FMAP).
C
C        Note that any effect POLYLINE might have on the plot legend (via
C     its call to CURVE) is up to the application.  The original DISSPLA
C     approach would initiate this by calling LEGLIN, LINEST, and LINES prior
C     to any call to POLYLINE.  The related LEGEND utility (an alternative
C     to DISSPLA's original LEGEND) eliminates any "hidden" connections between
C     POLYLINE/CURVE and the legend, so that awkward work-arounds for (say)
C     certain curve fitting methods are no longer needed.
C
C Arguments:
C ----------
C
C     Name    Dimension  Type  I/O/S  Description
C     N                   I    I      Number of points to be plotted.
C
C     X          N        R    I      Array of abscissas.
C
C     Y          N        R    I      Array of ordinates.
C
C     LINE                I    I      Line type code:
C                                        1 = connected symbols
C                                        2 = symbols only
C                                        3 = solid line
C                                        4 = dots
C                                        5 = dashes
C                                        6 = chain dots
C                                        7 = chain dashes
C                                        8 = long dashes
C                                        9 = thick solid line (.02 inch)
C
C     Notes:  1)  LINE < 0 indicates a closed curve.  E.g.  Line = -4 is a
C                 closed dotted line.
C             2)  Thickness factor is ignored if the line type code = 2
C                 (symbols only), or if the line type code = 9 (thick line 
C                 = .02 inch).  (LINE = 9 was retained for upward compatibility;
C                 LINE = 3 combined with THICK = 2 gives the same result.)
C
C     THICK               I    I      Line thickness factor (where pen width 
C                                     = THICK x .01 inch).  Note:  If THICK = 0,
C                                     DISSPLA default thickness is assumed.
C
C     SYMBOL              I    I      Symbol type.  Legal values are
C                                     [0,18], with automatic wraparound
C                                     if > 18.  See DISSPLA manual for
C                                     key.  Use negative for no symbol.
C
C     COLOR               R    I      Indicates the color to be
C                                     used for the curve:
C                                        0. = white
C                                        1. = black
C                                        2. = magenta
C                                        3. = red
C                                        4. = yellow
C                                        5. = green
C                                        6. = cyan
C                                        7. = blue
C
C                                      with colors between the six primaries 
C                                      specified using the real numbers in
C                                      [2., 7.].  (E.g. 3.5 = orange, 4.5 =
C                                      yellow-green, etc.)
C
C     IER                 I      O    Error flag:
C                                        0 = no problems
C                                        1 = N is not positive
C                                        2 = no line or symbol asked for
C
C Environment:  Digital VAX-11/780 VMS Version 4.4 (FORTRAN 77)
C ------------  Computer Associates DISSPLA Version 10.0
C
C Notes:
C ------
C     (1)  IMPLICIT NONE and 8-character symbolic names are non-standard.
C
C     (2)  See the DISSPLA manual for the order of the plotting symbols
C          used.  The first few are: Square, Octagon, Triangle ("up"),
C          Plus, X, Diamond, and Triangle ("down").
C
C Author:  Robert Kennelly, Sterling Software/NASA-Ames
C -------
C
C History:
C --------
C     08/27/87   R.A.Kennelly   Design and coding, based on PLCURV, fragments
C                               of QUICK, and incorporating pieces of MYSPEC.
C     06/23/89   M.D.Wong,      Header clarified in light of SMDLIB version;
C                D.A.Saunders   chaindot, chaindash redefined as for longdash.
C                               Symbols may now be plotted in combination with
C                               non-solid line patterns.
C     08/29/89   M.D.Wong       Updated to run on DISSPLA version 11.0. (%REF
C                               taken out of call to SETCLR.)
C     12/06/89   M.D.Wong       Ideas from FMCURV included here for contour
C                               plotting.  Line thicknesses may now be varied
C                               via the addition of THICK to the argument
C                               list.  Curve closing option added.  COLOR made 
C                               type real so that line colors may be specified 
C                               between the primaries.
C     02/05/90  D.A.Saunders    Use of RLVEC to close a curve was bad - does
C                               not acknowledge line-type or thickness.  Chose
C                               to call CURVE a second time rather than require
C                               a duplicate of the first point.  This could
C                               impact use of DISSPLA's legend scheme - use
C                               LEGEND from this library instead.
C-------------------------------------------------------------------------------


C     Declarations.
C     -------------

      IMPLICIT NONE

C     Arguments.

      INTEGER
     >   IER, LINE, N, SYMBOL, THICK
      REAL
     >   COLOR, X (N), Y (N)

C     Local Constants.

      REAL
     >   FACTOR, XTHK, TLENG
      PARAMETER
     >  (FACTOR = 1.0005,  ! Used to separate coincident points for DISSPLA. 
     >   XTHK = .01,       ! Line thickness multiplier in inches.
     >   TLENG = 0.2)      ! Overall line pattern length in inches.

C     Local variables.  

      INTEGER
     >   ICOLOR, LTYP, TYPE 
      REAL
     >   AINT, HUE, RATCHDSH (4), RATCHDOT (4), RATLODSH (2), SAT, X2,Y2
      CHARACTER
     >   COLORAY (0:1) * 5

C     Storage.

      DATA AINT, SAT /1., 1./              ! Fix intensity and saturation for
                                           ! colors between magenta and blue.
      DATA COLORAY   /'WHITE', 'BLACK'/

      DATA RATCHDSH  /1.7, 0.3, 0.7, 0.3/  ! Set line patterns.
      DATA RATCHDOT  /2.3, 0.3, 0.1, 0.3/
      DATA RATLODSH  /0.7, 0.3/

      SAVE
     >   AINT, COLORAY, RATCHDSH, RATCHDOT, RATLODSH, SAT


C     Execution.
C     ----------

      IER = 0

C     Simple-minded error checking.  These are not necessarily fatal, so
C     just return with error flag set.

      IF (N .LE. 0) THEN

C        No plot data?

         IER = 1
      ELSE IF (LINE .EQ. 2 .AND. SYMBOL .LT. 0) THEN

C        Invisible line!?

         IER = 2
      END IF

      IF (IER .NE. 0) GO TO 999

C     Preset to defaults.  Note that 'WHITE' actually shows up as black
C     in hardcopy.  Resetting 'DOT' takes care of all the patterns.

      CALL RESET ('SETCLR')
      CALL RESET ('THKCRV')
      CALL RESET ('DOT')


C     Set line color.
C     ---------------

      IF (COLOR .GE. 0.0 .AND. COLOR .LE. 7.0) THEN
         IF (COLOR .GE. 2.0) THEN   ! Assign HSI colors from magenta to blue.
            HUE = .5 * (COLOR - 1.0)
            CALL HWHSI (HUE, SAT, AINT)
         ELSE                       ! Handle black or white separately.
            ICOLOR = COLOR
            CALL SETCLR (COLORAY (ICOLOR))
         END IF
      END IF

C     Set line pattern.
C     -----------------

      IF (THICK .NE. 0) CALL THKCRV (THICK * XTHK)

      LTYP = ABS (LINE) 
      IF (LTYP .EQ. 4) THEN
         CALL DOT
      ELSE IF (LTYP .EQ. 5) THEN
         CALL DASH
      ELSE IF (LTYP .EQ. 6) THEN

C        Chaindot.  Don't use DISSPLA's CHNDOT (too long a pattern).
C        Define a "custom" line type - see DISSPLA manual, Part B, Section 9.
C        Alternating marks and spaces are prescribed by the ratios of the
C        segments to one another and by their total length.

         CALL MRSCOD (3.0 * TLENG, 4, RATCHDOT)

      ELSE IF (LTYP .EQ. 7) THEN

C        Chaindash.  Don't use DISSPLA's CHNDSH (too long a pattern).

         CALL MRSCOD (3.0 * TLENG, 4, RATCHDSH)

      ELSE IF (LTYP .EQ. 8) THEN   ! Longdash.  No DISSPLA utility available.

         CALL MRSCOD (TLENG, 2, RATLODSH)

      ELSE IF (LTYP .EQ. 9) THEN   ! Use default width for thick line.

         CALL THKCRV (2. * XTHK)

      END IF

C     Set point marker type.
C     ----------------------

C     Note that MARKER will look for a custom symbol if called with a
C     negative argument, so we have to test first.

      IF (SYMBOL .GE. 0) THEN 

         CALL MARKER (SYMBOL)

C        Set DISSPLA's TYPE flag for symbols and/or lines.

         IF (LTYP .NE. 2) THEN  ! Line plus symbols.
            TYPE = +1
         ELSE                   ! Symbols alone.
            TYPE = -1
         END IF
      ELSE                      ! Line only
         TYPE = 0
      END IF


C     Plot the curve, and close it if appropriate.
C     --------------------------------------------

C     DISSPLA is known to fail if the curve to be plotted is more
C     than a single thickness and of zero length.  Guard against
C     this by moving these points apart by some satisfactory distance.
C     Note:  Handle the N=2 case only - no attempt here to handle three or
C     more coincident points, since this involves calculating arc length.

      X2 = X (2)
      IF (N .EQ. 2) THEN
         IF (X (1) .EQ. X2 .AND. Y (1) .EQ. Y (2)) THEN
            X (2) = X2 * FACTOR
            IF (X2 .EQ. 0.) X(2) = 1.E-6
         END IF
      END IF

      CALL CURVE (X, Y, N, TYPE)

      IF (LINE .LT. 0 .AND. LTYP .NE. 2) THEN        ! Don't use RLVEC - it
CCCCC    CALL RLVEC ( X (N), Y (N), X (1), Y (1), 0) ! ignores thickness/type
         X (2) = X (N)
         Y2 = Y (2)
         Y (2) = Y (N)
         CALL CURVE (X, Y, 2, TYPE)     ! Preferable to requiring (N+1)th pt.,
         Y (2) = Y2                     ! but could affect DISSPLA legends -
      END IF                            ! use LEGEND from this library instead.
         
      X (2) = X2

C     Reset everything.
C     -----------------

      CALL RESET ('SETCLR')
      CALL RESET ('THKCRV')
      CALL RESET ('DOT')


C     Termination.
C     ------------

  999 RETURN
      END
