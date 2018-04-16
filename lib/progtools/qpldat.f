C+----------------------------------------------------------------------
C
      SUBROUTINE QPLDAT
C
C PURPOSE: QPLDAT writes QPLOTable data to the indicated file, including
C          titles, axis labels and most of the optional namelist inputs
C          accepted by QPLOT.
C
C METHOD:  Separate entry points are used to distinguish between new
C          frame information (plot titles, etc. - entry QPLFRM), new
C          curve information (namelist options  - entry QPLCRV), and
C          strictly numerical values (possibly points to be inserted
C          in an existing curve - entry QPLNXY).
C
C          Floating point arguments are checked against a flag (also
C          an argument) to suppress any unwanted variables in the
C          namelist.  Blanks should be passed for any of the character
C          variables not wanted in the namelist.
C
C ARGUMENTS FOR ENTRY QPLFRM:
C
C    ARG       DIM     TYPE I/O/S DESCRIPTION
C
C  LUN          -        I    I   Logical unit number for QPLOT file.
C  TITLE        -      C*(*)  I   Title for a plot frame; QPLOT has a
C                                 limit of 80 characters.
C  SBTITL       -      C*(*)  I   Subtitle for plot frame.
C  XLABEL,      -      C*(*)  I   Axis labels (also up to 80 chars.).
C  YLABEL
C
C ARGUMENTS FOR ENTRY QPLCRV:
C
C    ARG       DIM     TYPE I/O/S DESCRIPTION
C
C  LUN          -        I    I   Logical unit number for QPLOT file.
C  NPTS         -        I    I   Length of data arrays. NPTS = 0 permits
C                                 just the namelist options to be written.
C  X,Y        NPTS       R    I   Data representing one curve.
C  ENDSTR       -        C    I   Either 'END CURVE', 'END FRAME', or ' '.
C  OMIT         -        R    I   Flag for omitting floating point
C                                 variables from the namelist. Input
C                                 arguments equal to OMIT will not be
C                                 written to the output file.
C
C  Note: The remaining arguments are variables for QPLOT's optional namelist.
C
C  PLOT         -      C*(*)  I   Plot type. Options include 'DEFAULT' =
C                                 'LINEAR', 'SCALE', 'LOGX', 'LOGY',
C                                 'LOGLOG', and 'POLAR'.
C  YMIN,YMAX,   -        R    I   Bounds used for windowing or suppress-
C  XMIN,XMAX                      ing points out of desired range.
C  WIDTH        -        R    I   Desired length of X axis in inches.
C  HEIGHT       -        R    I   Desired length of Y axis in inches.
C  LEGEND       -      C*(*)  I   Text for legend entry for this curve.
C  LINE         -      C*(*)  I   Line type for this curve. Use 'DEFAULT'
C                                 or 'CONNECT' for symbols and lines. Other
C                                 choices are 'SYMBOLS', 'SOLID', 'DASH',
C                                 'DOT', 'THICK', 'CHAINDASH', 'CHAINDOT',
C                                 OR 'LONGDASH'.
C
C ARGUMENTS FOR ENTRY QPLNXY:
C
C    ARG       DIM     TYPE I/O/S DESCRIPTION
C
C  LUN          -        I    I   Logical unit number for QPLOT file
C  NPTS         -        I    I   Length of data arrays. Must be > 0.
C  X,Y        NPTS       R    I   Data points to be included.
C  ENDSTR       -        C    I   Either 'END CURVE', 'END FRAME', or ' '.
C
C ERROR HANDLING:  None.
C
C STANDARDS VIOLATIONS:  Multiple entry points make sense here.
C
C ENVIRONMENT:  Fortran 90
C
C HISTORY:
C
C  10/30/83  DAS  Initial implementation (point suppression done here).
C  01/31/84  LJC  Added writing of optional namelist, leaving QPLOT to
C                 suppress points out of desired range.
C  04/12/84  DAS  Made use of QPLOT's legend capability.
C  10/23/84  LJC  Added more namelist options and QPLNXY entry point.
C  06/20/86  DAS  ENDSTR*(*) instead of ENDSTR*9.
C  12/10/86  DAS  X(*), Y(*) instead of X(NPTS), ... for the NPTS=0
C                 case of QPLNXY.
C  06/18/91  DAS  Kludge for missing ORIENT option: if WIDTH is too
C                 big for portrait but small enough for landscape
C                 (8.5 x 11 in both cases), insert ORIENT='LANDSCAPE'.
C  10/18/99  DAS  SGI list-directed I/O started putting out things like
C                 2*0.E+0, which cannot be parsed by QPLOT.
C
C AUTHOR: David Saunders, Sterling Software, Palo Alto, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT  NONE

C ... Arguments:

      INTEGER
     >   LUN, NPTS

      REAL
     >   HEIGHT, OMIT, WIDTH, X (*), XMAX, XMIN, Y (*), YAXIS, YMAX,
     >   YMIN

      CHARACTER
     >   ENDSTR * (*), LEGEND * (*), LINE * (*), PLOT * (*),
     >   SBTITL * (*), TITLE * (*), XLABEL * (*), YLABEL * (*)

C ... Local constants:

      CHARACTER * 1, PARAMETER ::
     >   BLANK = ' '

C ... Local variables:

      INTEGER
     >   I

      REAL
     >   EPS

C-----------------------------------------------------------------------

      ENTRY QPLFRM (LUN, TITLE, SBTITL, XLABEL, YLABEL)

C-----------------------------------------------------------------------

C ... Write the title and axis labels for a plot frame.

      WRITE (LUN, 100) TITLE
      WRITE (LUN, 100) SBTITL
      WRITE (LUN, 100) XLABEL
      WRITE (LUN, 100) YLABEL
      GO TO 99


C-----------------------------------------------------------------------

      ENTRY QPLCRV (LUN, NPTS, X, Y, ENDSTR, OMIT, PLOT, XMIN, XMAX,
     >              YMIN, YMAX, WIDTH, HEIGHT, LEGEND, LINE)

C-----------------------------------------------------------------------

C ... Write options in namelist format:

      WRITE (LUN, 100) ' $OPTIONS'

      IF (PLOT  /= BLANK) WRITE (LUN, 120) 'PLOT', PLOT
      IF (XMIN  /=  OMIT) WRITE (LUN, 110) 'XMIN', XMIN
      IF (XMAX  /=  OMIT) WRITE (LUN, 110) 'XMAX', XMAX
      IF (YMIN  /=  OMIT) WRITE (LUN, 110) 'YMIN', YMIN
      IF (YMAX  /=  OMIT) WRITE (LUN, 110) 'YMAX', YMAX
      IF (WIDTH /=  OMIT) THEN
         WRITE (LUN, 110) 'WIDTH',  WIDTH
         IF (WIDTH > 8.0 .AND. WIDTH < 11.0) THEN
            WRITE (LUN, 120) 'ORIENT', 'LANDSCAPE'
         END IF
      END IF
      IF (HEIGHT/=  OMIT) WRITE (LUN, 110) 'HEIGHT', HEIGHT
      IF (LEGEND/= BLANK) WRITE (LUN, 120) 'LEGEND', LEGEND
      IF (LINE  /= BLANK) WRITE (LUN, 120) 'LINE',   LINE

      WRITE (LUN, 100) ' $END'
 
C-----------------------------------------------------------------------

      ENTRY QPLNXY (LUN, NPTS, X, Y, ENDSTR)

C-----------------------------------------------------------------------

C ... Save data points:

      IF (NPTS > 0) THEN

         EPS = EPSILON (EPS) ! Avoid 2*0.E+0, etc., from SGI list-directed write

         IF (EPS > 1.E-10) THEN ! Assume 32-bit compile

            WRITE (LUN, '(1P, 2E15.7)') (X (I), Y (I), I = 1, NPTS)

         ELSE

            WRITE (LUN, '(1P, 2E23.15)') (X (I), Y (I), I = 1, NPTS)

         END IF

      END IF

      IF (ENDSTR /= BLANK) WRITE (LUN, 100) ENDSTR

   99 RETURN

C ... Formats:

  100 FORMAT (A)
  110 FORMAT (1X, A, ' = ', 1P, E13.6, ',')
  120 FORMAT (1X, A, ' = ''', A, ''',')

      END SUBROUTINE QPLDAT
