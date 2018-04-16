C+----------------------------------------------------------------------
C
      SUBROUTINE ANNOTATE (NUMNOTE, NOTES, XNOTES, YNOTES, XARROWS,
     >                     YARROWS, NOTEBOX, FLAG, LUNERR)
C
C ONE-LINER:  Annotate a plot with or without arrows (DISSPLA level 2, 3)
C
C PURPOSE:
C
C        ANNOTATE adds NUMNOTE text strings to a plot, with the lower
C     left point of string I defined by (XNOTES (I), YNOTES (I)) in
C     data units.  If (XARROWS (I), YARROWS (I)) is also defined
C     (= (FLAG, FLAG) if undefined), an arrow is drawn from the CENTER
C     of the text string rectangle.  Each arrow is clipped by a border
C     of width .05" around the annotation.  The area within this border
C     may be outlined and/or blanked (to prevent grid lines, etc. from
C     intruding).
C
C METHOD:
C
C        DISSPLA's MESSAG draws text at a point specified in inches,
C     and we need to convert from data units if an arrow is involved,
C     so choose MESSAG over RLMESS, which could draw the text too.
C     Clipping of the optional arrows is done explicitly in case no
C     blanking is requested.  Of many possible arrow types, the one
C     defined by the integer 1001 in DISSPLA's terminology for its
C     VECTOR utility is the only one provided here (i.e. single arrow
C     head, shaded, with modestly sharp point).  DISSPLA's BLREC blanks
C     each rectangle once the text is drawn, if blanking is requested.
C
C        The current character height (retrieved via WHATIS for the
C     rectangle calculations) applies to both the height of the text
C     and the length of the arrow head.
C
C ERROR HANDLING:
C
C        If either coordinate for an annotation is undefined (= FLAG),
C     an error message is written to LUNERR, and the annotation is skipped.
C
C ENVIRONMENT:
C     VAX/VMS, FORTRAN 77, with:
C     > IMPLICIT NONE
C     > Trailing ! comments
C     > Names up to 8 characters
C
C HISTORY:
C     09/27/91  D.A.Saunders  Initial design and code, for QPLOT.
C     10/19/91    "     "     Achieving precisely horizontal or vertical
C                             arrows was awkward originally.  Now,
C                             YARROWS (I) = YNOTES (I) is taken to mean
C                             "horizontal" (whether the text is blank or
C                             not).  Likewise, XARROWS (I) = XNOTES (I)
C                             means "vertical" now, from center of text.
C
C AUTHOR: David Saunders, Sterling Software/NASA Ames, Mountain View, CA.
C
C ----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments.

      INTEGER   NUMNOTE         ! (I) Number of annotations.  NUMNOTE >= 0.
      CHARACTER NOTES (*) * (*) ! (I) Annotation text strings.  If blank,
                                !        an arrow may still be drawn between
                                !        the specified points.
      REAL      XNOTES (*)      ! (I) Coordinates of bottom left corner of
      REAL      YNOTES (*)      !        annotation text, in data units.
      REAL      XARROWS (*)     ! (I) Coordinates of tips of arrow heads,
      REAL      YARROWS (*)     !        in data units.  An arrow is suppressed
                                !        if either coordinate equals FLAG.
                                !        To ensure a HORIZONTAL arrow, set
                                !        YARROWS (I) = YNOTES (I) (whether
                                !        NOTES (I) is blank or not).  Likewise,
                                !        set XARROWS (I) = XNOTES (I) to achieve
                                !        a precisely VERTICAL arrow from the
                                !        center of text.
      INTEGER   NOTEBOX         ! (I) Controls boxes around annotations:
                                !        0 means no frame and no blanking;
                                !        1 means no frame/rectangle is blanked;
                                !        2 means box is framed and blanked;
                                !       >2 means frame is thickened outwards:
                                !          (NOTEBOX - 1) scales the pen width.
      REAL      FLAG            ! (I) Flag for undefined coordinates. E.g. 999.
      INTEGER   LUNERR          ! (I) Logical unit number for error messages.

C     Procedures.

      REAL      XMESS
      EXTERNAL  BLREC           ! (DISSPLA) Blanks a non-tilted rectangle.
      EXTERNAL  MESSAG          ! (DISSPLA) Draws text at (X, Y) inches.
      EXTERNAL  SCAN2           ! Determines lengths of annotation strings.
      EXTERNAL  TOINCH          ! (DISPPLA) Converts from data units to inches.
      EXTERNAL  VECTOR          ! (DISSPLA) Draws an arrow.
      EXTERNAL  WHATIS          ! (DISSPLA) Used to determine HEIGHT setting.
      EXTERNAL  XMESS           ! (DISSPLA) Length of text string in inches.

C-----------------------------------------------------------------------

C     Local constants.

      INTEGER   IVEC
      REAL      BORDER, EPS, HALF
      CHARACTER BLANK * 1
      PARAMETER
     >  (BORDER = 0.05,         ! Width of borders around text (inches)
     >   EPS    = .01,          ! Need to avoid TAN (90 deg.)
     >   HALF   = 0.5E+0,
     >   IVEC   = 1001,         ! See IVEC arrow definition in Ch. 5, Sec. 2.1
     >   BLANK  = ' ')

C     Local variables.

      INTEGER
     >   FIRST, I, LAST, LENNOTE, MARK
      REAL
     >   FRM, HIGH, PHI, PI, TANTH, THETA, XA, XC, XE, XL, XR,
     >   YA, YB, YC, YE, YT
      LOGICAL
     >   HORIZ, VERTI, VERTICAL

C     Execution.

      LENNOTE = LEN (NOTES (1))
      PI = 4.0 * ATAN2 (1.0, 1.0)
      FRM = REAL (1 - NOTEBOX)   ! 0 means no frame, but text area is blanked.
                                 ! -ve causes frame thickening toward OUTside.

      DO 500, I = 1, NUMNOTE

         XL = XNOTES (I)
         YB = YNOTES (I)

         IF (XL .EQ. FLAG .OR. YB .EQ. FLAG) THEN

            WRITE (LUNERR, 1010) I

         ELSE

            XA = XARROWS (I)
            YA = YARROWS (I)
            HORIZ = YA .EQ. YB
            VERTI = XA .EQ. XL

            CALL TOINCH (XL, YB, XL, YB) ! Convert lower left coords. of string
                                         ! to inches, as needed for box & arrow

C           Draw the annotation unless it is blank.

            FIRST = 1
            LAST = LENNOTE
            CALL SCAN2 (NOTES (I), BLANK, FIRST, LAST, MARK)

            IF (LAST .GT. 0) THEN
               CALL MESSAG (NOTES (I) (1 : LAST), LAST, XL, YB)
            END IF

C           Determine text box border for blanking purposes, even if
C           no arrow is drawn.  XL, YB are modified for non-blank text, so:

            XE = XL                 ! End of arrow, if arrow is
            YE = YB                 ! requested but string is blank

            XR = XL + XMESS (NOTES (I), LAST) + BORDER
            XL = XL - BORDER

            CALL WHATIS ('HEIGHT', BLANK, MARK, HIGH, MARK)

            YT = YB + HIGH + BORDER
            YB = YB - BORDER

C           Draw any arrow (explicitly clipped) before any blanking of the text.

            IF (XA .NE. FLAG .AND. YA .NE. FLAG) THEN

               CALL TOINCH (XA, YA, XA, YA)

               IF (LAST .GT. 0) THEN

C                 The arrow (suitably clipped) will emanate from the center
C                 (XC, YC) of the rectangle (with border) containing the string.

                  XC = (XL + XR) * HALF
                  YC = (YB + YT) * HALF

                  IF (HORIZ) YA = YC               ! Fudge - otherwise the user
                  IF (VERTI) XA = XC               ! has to mix data units and
                                                   ! plot inches to estimate the
                                                   ! center of the text box.

                  PHI = ATAN2 (YT - YC, XR - XC)   ! Angle of line from center
                                                   ! to top right of box in
                                                   ! the range (0, PI/2 - del) 

                  THETA = ATAN2 (YA - YC, XA - XC) ! Angle of arrow in the
                                                   ! range [-PI, PI].

                  VERTICAL = ABS (ABS (THETA) - PI * HALF) .LT. EPS
                  IF (VERTICAL) THEN
                     XE = XC
                  ELSE
                     TANTH = (YA - YC) / (XA - XC)
                  END IF

C                 Determine which side of the rectangle clips the arrow,
C                 and adjust the end of the arrow (XE, YE) accordingly.

                  IF (ABS (THETA) .LE. PHI) THEN             ! Right side clips.
                     XE = XR
                     YE = YC + (XR - XC) * TANTH
                  ELSE IF (THETA .GT. PHI .AND.
     >                     THETA .LE. PI - PHI) THEN         ! Top clips.
                     IF (.NOT. VERTICAL) XE = XC + (YT - YC) / TANTH
                     YE = YT
                  ELSE IF (ABS (PI - THETA) .LE. PHI) THEN   ! Left side clips.
                     XE = XL
                     YE = YC + (XL - XC) * TANTH
                  ELSE                                       ! Bottom clips.
                     IF (.NOT. VERTICAL) XE = XC + (YB - YC) / TANTH
                     YE = YB
                  END IF
               END IF

C              Draw the arrow unless it is absurdly short.

               IF ((XA - XE) ** 2 + (YA - YE) ** 2 .GT. HIGH ** 2)
     >            CALL VECTOR (XE, YE, XA, YA, IVEC)

            END IF

C           Blank and/or outline the text area?

            IF (FRM .LE. 0.0) CALL BLREC (XL, YB, XR - XL, YT - YB, FRM)

         END IF

  500 CONTINUE

      RETURN

C     Formats.

 1010 FORMAT ('0*** ANNOTATE: Undefined coordinate(s) for annotation #',
     >        I3, '.  Proceeding...')

      END
