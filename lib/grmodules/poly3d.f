C+------------------------------------------------------------------------------
C
      SUBROUTINE POLY3D (N, X, Y, Z, LINE, THICK, SYMBOL, COLOR, IER)
C
C     Purpose: 3-space version of POLYLINE, q.v.
C
C     Arguments:
C
C     Name    Dimension  Type  I/O/S  Description
C     N                   I    I      Number of points to be plotted.
C
C     X,         N        R    I      Cartesian coordinates of
C     Y,         N        R    I      points to be plotted.
C     Z          N        R    I
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
C Environment:  VAX/VMS; FORTRAN 77; CA-DISSPLA Version 11.0
C
C Notes:
C     (1)  IMPLICIT NONE and 8-character symbolic names are non-standard.
C
C     (2)  See the DISSPLA manual for the order of the plotting symbols
C          used.  The first few are: Square, Octagon, Triangle ("up"),
C          Plus, X, Diamond, and Triangle ("down").
C
C Author:  Robert Kennelly, NASA Ames
C
C History:
C     08/27/87   R.A.Kennelly   Original POLYLINE (2-space).
C     02/18/91   D.A.Saunders   POLY3D adapted from POLYLINE.
C
C-------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments.

      INTEGER
     >   IER, LINE, N, SYMBOL, THICK
      REAL
     >   COLOR, X (N), Y (N), Z (N)

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
     >   AINT, HUE, RATCHDSH (4), RATCHDOT (4), RATLODSH (2), SAT,
     >   X2, Y2, Z2
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
C     in some hardcopy (but not for POstScript).  Resetting 'DOT' takes
C     care of all the patterns.

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

      CALL CURV3D (X, Y, Z, N, TYPE)

      IF (LINE .LT. 0 .AND. LTYP .NE. 2) THEN ! Don't use RLVEC3 - it ignores
CCCCC    CALL RLVEC3 ( X (N), Y (N), Z (N), X (1), Y (1), Z (1), 0) ! thick/type
         X (2) = X (N)
         Y2 = Y (2)
         Y (2) = Y (N)
         Z2 = Z (2)
         Z (2) = Z (N)
         CALL CURV3D (X, Y, Z, 2, TYPE) ! Preferable to requiring (N+1)th pt.,
         Y (2) = Y2                     ! but could affect DISSPLA legends -
         Z (2) = Z2
      END IF                            ! use LEGEND from this library instead.
         
      X (2) = X2

C     Reset everything.

      CALL RESET ('SETCLR')
      CALL RESET ('THKCRV')
      CALL RESET ('DOT')

  999 RETURN
      END
