C+------------------------------------------------------------------------------
C
      SUBROUTINE LEGEND (TEXT, LENTXT, LINE1, LINE2, DY, XCORNR, YCORNR,
     >                   LEGINF, HITE, SEG, GAP, BOX)
C
C ONE-LINER:  Portable plot legend utility (Level 3)
C
C PURPOSE:
C
C        LEGEND composes a plot legend and displays it at a specified position.
C     If the number of legend entries is large and the texts are suitably short,
C     LEGEND may be called more than once to position different sections of the
C     legend alongside each other.  (Such positioning is up to the calling
C     program.)  An optional box may be drawn around the legend.
C
C        LEGEND was developed as an alternative to existing DISSPLA and SMDLIB
C     schemes for the following reasons:
C
C     (1) rather than underlining the text with the line pattern, the following
C         form is provided:           ---o---  Description
C     (2) CHARACTER-type text is expected instead of an awkward packed-integer
C         array data structure;
C     (3) any implicit connection between curve-drawing and updating the legend
C         is avoided, since this has required work-arounds in the past.
C
C METHOD:
C
C        The symbol/line drawing is handled separately from the writing the
C     text as a left-justified block.  The latter task was modularized when
C     it was found to be a potentially reusable function (subroutine LBTEXT)
C     while the former is done with the standard curve drawing routine,
C     POLYLINE, after appropriate conversion of units from inches to data
C     units via XINVRS and YINVRS.
C
C        Character height and line spacing are assumed constant for each call
C     to LEGEND.  The text color is the calling program's current color.
C
C FURTHER USAGE NOTES:
C
C        The (XCORNR, YCORNR) coordinates are in inches relative to the
C     current plot origin.  Carefully-centered (or otherwise positioned)
C     legends almost certainly require something like the following in the
C     calling program:
C
C        LWIDTH = XDIMTB (TEXT, LENTXT, LINE1, LINE2) + SEG + GAP
C        XCORNR = 0.5 * (WIDTH - LWIDTH)
C
C        WARNING:  If the legend is to be drawn outside the plotting area,
C     a suitable grace margin should be defined before the call to LEGEND.
C
C ARGUMENTS:
C
C     ARG        DIM     TYPE  I/O/S  DESCRIPTION
C     TEXT    (*) * (*)   C        I  Character array containing legend text.
C                                     (Terminating '$'s may be required,
C                                     depending on the usage of LENTXT (*).)
C                                     Elements LINE1:LINE2 will be displayed in
C                                     the current color (as opposed to the color
C                                     associated with each line/symbol).
C     LENTXT     (*)      I        I  Array containing number of characters
C                                     in text line(s), if known.  (No trailing
C                                     '$'s are needed in this case.)  Otherwise,
C                                     pass 100 as the first (and only) value,
C                                     and all lines of text will be self counted
C                                     (requiring the trailing '$'s).
C     LINE1       -       I      I    First line of legend text to display.
C     LINE2       -       I      I    Last line of legend text to display.
C     DY          -       R      I    Space between lines in inches.
C     XCORNR      -       R      I    Horizontal distance of lower left hand
C                                     corner of legend from origin in inches.
C     YCORNR      -       R      I    Vertical distance of lower left hand
C                                     corner of legend from origin in inches.
C     LEGINF    (3, *)    I      I    Integer array containing codes for
C                                     (1) line type, (2) symbol, and (3) color.
C                                     See POLYLINE for the code definitions.
C     HITE        -       R      I    Text character height in inches.
C     SEG         -       R      I    Length of legend line segment in inches.
C                                     0.54" is a good value: it suppresses the
C                                     last space of several interrupted line
C                                     patterns - see POLYLINE.
C     GAP         -       R      I    Length of gap between legend line segment
C                                     and text in inches. 0.25" is typical.
C     BOX         -       L      I    .TRUE. means draw a box around legend.
C
C ERROR HANDLING:  There is no check for calling LEGEND at the wrong level.
C
C PROCEDURES:
C     LBTEXT    Writes left justified block of text
C     XDIMTB    Finds maximum width of a text block. Used here only if a
C               BOX has been requested.
C     XINVRS    Converts location of a point from inches to data units
C     YINVRS       "              "      "          "         "
C     POLYLINE  Draws a curve using lines and/or symbols
C
C ENVIRONMENT:  VAX/VMS; FORTRAN 77; DISSPLA or SMDLIB
C
C HISTORY:
C     04/17/89    M.D.Wong     Initial implementation for QPLOT translation
C                              from DISSPLA to SMDLIB, using internal COMMONs
C                              and DSDRAW.
C     06/22/89    M.D.Wong     Modularized line/symbol drawing into POLYLINE,
C                              following introduction of conversion utilities
C                              XINVRS and YINVRS.  Eliminated internal COMMONs
C                              at the expense of not checking calling level and
C                              not restoring original line type in order to
C                              combine DISSPLA and SMDLIB versions as one.
C     01/30/90    M.D.Wong     Variable COLOR now passed as REAL to POLYLINE.
C     02/26/90    M.D.Wong     Added LENTXT to argument list to avoid reliance
C                              on DISSPLA dependent self counting option.
C     04/12/90    M.D.Wong     Modified LBTEXT to enable 100 to be passed as
C                              the first and only element of the LENTXT array
C                              to enable self-counting option.  No changes
C                              were necessary to LEGEND.
C     03/27/91   D.A.Saunders  Clarified usage in light of QPLOT revisions.
C
C AUTHOR:  Michael Wong, Sterling Software, Palo Alto, CA.
C
C-------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments.
C     ----------

      CHARACTER  TEXT (*) * (*)
      INTEGER    LENTXT (*), LEGINF (3, *), LINE1, LINE2
      REAL       DY, HITE, XCORNR, YCORNR
      LOGICAL    BOX

C     Local Constants.
C     ----------------

      REAL       HALF
      PARAMETER (HALF = 0.5)

C     Local Variables.
C     ----------------

      INTEGER    COLOR, I, IER, LINE, SYMBOL
      REAL       GAP, MAR, SEG, X (5), XBLEN, XPOS, Y (5), YBLEN, YPOS

C     Procedures.
C     -----------

      REAL       XDIMTB, XINVRS, YINVRS
      EXTERNAL   XDIMTB, XINVRS, YINVRS, LBTEXT, POLYLINE

C     Execution.
C     ----------


C     Display the legend text as a left-justified block in current color.

      CALL LBTEXT (TEXT, LENTXT, LINE1, LINE2, DY, XCORNR + SEG + GAP,
     >             YCORNR, HITE)


C     Draw a symbol and/or line segment for each legend entry.
C     POLYLINE applies if positions in inches are converted to data units.
C     Any grace margin needed for a legend outside the plot area is left
C     to the application program, since there is no way of restoring the
C     input setting here if it is temporarily adjusted.
C
C     Note that the polar plot case (not to mention general rotated coordinate
C     systems) forces all conversions from inches to stay inside the loop
C     (where some would be constant otherwise).


C     Start with the upper left corner.

      YPOS = YCORNR + (LINE2 - LINE1) * (HITE + DY) + HALF * HITE

      DO 200, I = LINE1, LINE2

         X (1)  = XINVRS (XCORNR, YPOS)
         Y (1)  = YINVRS (XCORNR, YPOS)
         X (2)  = XINVRS (XCORNR + SEG, YPOS)
         Y (2)  = YINVRS (XCORNR + SEG, YPOS)
         LINE   = LEGINF (1, I)
         SYMBOL = LEGINF (2, I)
         COLOR  = LEGINF (3, I)

         IF (LINE .NE. 2) THEN

C           Draw the line segment.

            CALL POLYLINE (2, X, Y, LINE, 0, -1, REAL (COLOR), IER)

         END IF

         IF (SYMBOL .GT. -1) THEN

C           Draw the symbol in the middle of the segment.

            X (1) = XINVRS (XCORNR + HALF * SEG, YPOS)
            Y (1) = YINVRS (XCORNR + HALF * SEG, YPOS)

            CALL POLYLINE (1, X (1), Y (1), 2, 0, SYMBOL, REAL (COLOR),
     >                     IER)

         END IF

         YPOS = YPOS - (HITE + DY)

  200 CONTINUE


      IF (BOX) THEN

C        Pick a reasonable box margin.

         MAR  = .6 * GAP
         XPOS = XCORNR - MAR
         YPOS = YCORNR - MAR
         MAR  = MAR + MAR

C        Calculate the box dimensions (inches).

         XBLEN = XDIMTB (TEXT, LENTXT, LINE1, LINE2) + SEG + GAP + MAR
         YBLEN = REAL (LINE2 - LINE1) * (HITE + DY) + HITE + MAR

C        Draw the box starting from the lower left corner.

         X (1) = XINVRS (XPOS, YPOS)
         Y (1) = YINVRS (XPOS, YPOS)
         X (2) = XINVRS (XPOS + XBLEN, YPOS)
         Y (2) = YINVRS (XPOS + XBLEN, YPOS)
         X (3) = XINVRS (XPOS + XBLEN, YPOS + YBLEN)
         Y (3) = YINVRS (XPOS + XBLEN, YPOS + YBLEN)
         X (4) = XINVRS (XPOS, YPOS + YBLEN)
         Y (4) = YINVRS (XPOS, YPOS + YBLEN)
         X (5) = X (1)
         Y (5) = Y (1)

         CALL POLYLINE (5, X, Y, 3, 0, -1, 0, IER)

      END IF

      RETURN
      END
