C=======================================================================
C
C     This collection of PostScript plot utilities appear to have been
C     translated from some earlier graphics package. They have also been
C     modified in ways that seem inconsistent with the original intent.
C     The present design and implementation still leave a lot to be
C     desired. However, extensive use of these utilities has been made
C     in the SYN87 wing/body design code for producing basic plots
C     directly, and they could prove helpful to other applications.
C
C     Written <when?> by <who?> for Prof. Antony Jameson, Princeton U.
C
C     Pre-7/94  J.J.Reuther   Various mods.  Some hard-coding.
C     07/29/94  D.A.Saunders  Cleaned it up. Introduced short
C                             forms for lineto, moveto, etc.
C     08/18/94  JJR/DAS       Added solid/dot/dash definitions.
C     Sept. 95  DAS           Introduced SYMTYPE (& LINETYPE)
C                             in an effort to reduce verbosity.
C                             MLABEL renamed as CSTRING; FONT
C                             made an argument in a few places.
C     01/30/96  JJR/DAS       Pages numbered to speed previewing.
C     11/20/96     "          Added simple-minded RGBCOLOR.
C     11/29/96  DAS           Added CTEXT.
C
C     Present Support:  David Saunders, Sterling Software/NASA Ames
C
C=======================================================================
C
      SUBROUTINE INITPL (XINCH, YINCH, SCALE, IPAGE)
C
C     INITPL outputs a PostScript prologue.  The application needs to
C     include COMMON /PSC/ JWRIT and initialize JWRIT to a plot file
C     logical unit number. Otherwise, the utilities are argument-driven.
C     XINCH and YINCH specify the (0,0) origin on the plot page (lower
C     left hand corner.  SCALE will probably be 1.  IPAGE needs to be
C     zeroed at the start of plotting to a given file; it is incremented
C     here.  Use INITPL and ENDPLT to start and end each plot page.
C
C=======================================================================

      COMMON /MIC/ NEN
      COMMON /PSC/ JWRIT

      NEN = 0               ! See PLOT
      WRITE (JWRIT, '(A)')
     >   '%!PS-Adobe-2.0'   ! Ghostview needs PS-Adobe- to do paging
      IF (IPAGE .EQ. 0) THEN
         WRITE (JWRIT, '(A)') '%%Pages: 999'  ! Ghostview uses EOF
      END IF
      WRITE (JWRIT, '(A)')
     >   '/inch {72 mul} def',
     >   '/l  /lineto  load def',
     >   '/m  /moveto  load def',
     >   '/r  /rotate  load def',
     >   '/rm /rmoveto load def',
     >   '/sh /show    load def',
     >   '/st /stroke  load def',
     >   '/solid     {[]            0 setdash} def',
     >   '/dot       {[0.001 0.01]  0 setdash} def',
     >   '/dash      {[0.05]        0 setdash} def',
     >   '/longdash  {[0.10]        0 setdash} def',
     >   '0 setlinewidth'

      IPAGE = IPAGE + 1
      WRITE (JWRIT, 1) XINCH, YINCH, SCALE, SCALE, IPAGE

      RETURN

    1 FORMAT (2(F5.2,' inch'),' translate ',
     >        2(F5.2,' inch'),' scale', /, '%%Page:', I3)
      END
C=======================================================================
C
      SUBROUTINE XAXIS (XMIN, XMAX, XINT, XORIG, YORIG, XLONG, ITIC,
     >                  ILABEL, XTITLE, NTITLE, FONT, ITICSIZE, IEND)
C
C     XAXIS draws a horizontal axis and returns a scale factor for
C     plotting variables on this axis.  (? No such factor is returned now.)
C     Absolute size of tics, labels and titles is set by SIZE
C     (normalized NDC coordinates).     (? The SIZE argument was not being
C     used until it was changed to refer to the tic size (of all things).)
C     Relative size of tics, labels and titles is set by internal
C     parameters which may be adjusted.
C
C     XMIN,  XMAX,  XINT  : min, max & intervals of x along axis
C     XORIG, YORIG, XLONG : origin & length of x axis in window coordinates
C     ITIC = -1  ----->     tics on negative y side of axis
C     ITIC =  0  ----->     no tics
C     ITIC =  1  ----->     tics on positive y side of axis
C     ITIC =  2  ----->     tics on both sides of axis
C     ILABEL=-1  ----->     x-coord labels on negative y side of axis
C     ILABEL= 0  ----->     no x-coord labels
C     ILABEL= 1  ----->     x-coord labels on positive y side of axis
C     ILABEL=+/-2  --->     same as +/-1 except x-coord labels are integers
C     XTITLE              : Title for x axis
C     NTITLE              : Number of characters in Title (Omit title if =0)
C     FONT                : String such as Times_Roman or Symbol
C     ITICSIZE            : Size of tics.  E.g. ITICSIZE = 8 gives 8/72"
C     IEND = -2  ----->     omit all x-coord labels except at ends of axis
C     IEND = -1  ----->     omit x-coord label at axis origin
C     IEND =  0  ----->     write all x-coord labels
C     IEND =  1  ----->     omit x-coord label at end of axis
C     IEND =  2  ----->     omit x-coord label at both ends of axis
C
C=======================================================================

      CHARACTER XTITLE * (*), FONT * (*)
      CHARACTER FLABEL * 7

      TICSIZ = REAL (ITICSIZE) / 72.
      CHASIZ = 1.25 * TICSIZ
      DIST   = 3.50 * TICSIZ   ! distance between labels and axis
      IF (ILABEL .LT. 0) DIST = -DIST + TICSIZ
      RELSIZ = 1.1             ! size of axis title versus tic labels

C     Draw the axis:

      CALL MOVEAB (XORIG, YORIG)

      XEND = XORIG + XLONG

      CALL LINEAB (XEND, YORIG)

      IF (ITIC .EQ. 0) GO TO 30

C     Draw the tics and (possibly) label them:

      IF (ITIC .GT. 0) THEN
         TIC = TICSIZ
      ELSE
         TIC = -TICSIZ
      END IF

      XPOS  = XORIG
      XSCAL = XLONG / (XMAX - XMIN)
      DX    = XINT * XSCAL
      N     = INT (XLONG / DX + 1.01)

      CALL SYMTYPE (FONT, CHASIZ)

      DO I = 1, N

         CALL MOVEAB (XPOS, YORIG)
         CALL LINEAB (XPOS, YORIG + TIC)
         IF (ITIC .EQ. 2) CALL LINEAB (XPOS, YORIG - TIC)

C        Move to center location of tic label and write x-value:

         IF (ILABEL .EQ. 0) GO TO 19
         IF (I .EQ. 1  .AND.  IEND .EQ.-1) GO TO 19
         IF (I .EQ. 1  .AND.  IEND .EQ. 2) GO TO 19
         IF (I .EQ. N  .AND.  IEND .GT. 0) GO TO 19
         IF (I .GT. 1  .AND.  I .LT. N  .AND. IEND .EQ. -2) GO TO 19

         CALL MOVEAB (XPOS, YORIG + DIST)

         XLABL = XMIN + REAL (I-1) * XINT
         IF (ABS (ILABEL) .EQ. 2) THEN
             IXLABL = INT (XLABL)
             WRITE (FLABEL, '(I5)') IXLABL
             NCHAR  = 5
         ELSE
             WRITE (FLABEL, '(F7.2)') XLABL
             NCHAR  = 7
         END IF

         CALL CSTRING (FLABEL, ' ', CHASIZ, NCHAR, 0., 1)

   19    XPOS = XPOS + DX
      END DO

   30 CONTINUE

C     Axis title:

      IF (NTITLE .GT. 0) THEN
         CALL MOVEAB (0.5 * (XORIG + XEND), YORIG + (1.+ RELSIZ) * DIST)
         CALL CSTRING (XTITLE, FONT, CHASIZ * RELSIZ, NTITLE, 0., 1)
      END IF

      RETURN
      END
C=======================================================================
C
      SUBROUTINE YAXIS (YMIN, YMAX, YINT, XORIG, YORIG, YLONG, ITIC,
     >                  ILABEL, YTITLE, NTITLE, FONT, ITICSIZE, IEND)
C
C     YAXIS draws a vertical axis and returns a scale factor for
C     plotting variables on this axis.  (? No such factor is returned now.)
C     Absolute size of tics, labels and titles is set by SIZE
C     (normalized NDC coordinates).     (? The SIZE argument was not being
C     used until it was changed to refer to the tic size (of all things).)
C     Relative size of tics, labels and titles is set by internal
C     parameters which may be adjusted.
C
C     YMIN,  YMAX,  YINT  : min, max & intervals of y along axis
C     XORIG, YORIG, YLONG : origin & length of y axis in window coordinates
C     ITIC = -1  ----->     tics on negative x side of axis
C     ITIC =  0  ----->     no tics
C     ITIC =  1  ----->     tics on positive x side of axis
C     ITIC =  2  ----->     tics on both sides of axis
C     ILABEL=-1  ----->     y-coord labels on negative x side of axis
C     ILABEL= 0  ----->     no y-coord labels
C     ILABEL= 1  ----->     y-coord labels on positive x side of axis
C     ILABEL=+/-2  --->     same as +/-1 except y-coord labels are integers
C     YTITLE              : Title for y axis
C     NTITLE              : Number of characters in Title (Omit title if =0)
C     FONT                : String such as Times_Roman or Symbol
C     ITICSIZE            : Size of tics.  E.g. ITICSIZE = 8 gives 8/72"
C     IEND = -2  ----->     omit all y-coord labels except at ends of axis
C     IEND = -1  ----->     omit y-coord label at axis origin
C     IEND =  0  ----->     write all y-coord labels
C     IEND =  1  ----->     omit y-coord label at end of axis
C     IEND =  2  ----->     omit y-coord label at both ends of axis
C
C=======================================================================

      CHARACTER YTITLE * (*), FONT * (*)
      CHARACTER FLABEL * 7

      TICSIZ = REAL (ITICSIZE) / 72.
      CHASIZ = 1.25 * TICSIZ
      DIST   = 3.50 * TICSIZ   ! distance between labels and axis
      IF (ILABEL .LT. 0) DIST = -DIST + TICSIZ
      RELSIZ = 1.1             ! size of axis title versus tic labels

C     Draw the axis:

      CALL MOVEAB (XORIG, YORIG)

      YEND = YORIG + YLONG

      CALL LINEAB (XORIG, YEND)

      IF (ITIC .EQ. 0) GO TO 30

C     Draw the tics and (possibly) label them:

      IF (ITIC .GT. 0) THEN
         TIC = TICSIZ
      ELSE
         TIC = -TICSIZ
      END IF

      YPOS  = YORIG
      YSCAL = YLONG / (YMAX - YMIN)
      DY    = YINT * YSCAL
      N     = INT (YLONG / DY + 1.01)

      CALL SYMTYPE (FONT, CHASIZ)

      DO I = 1, N

         CALL MOVEAB (XORIG, YPOS)
         CALL LINEAB (XORIG + TIC, YPOS)
         IF (ITIC .EQ. 2) CALL LINEAB (XORIG - TIC, YPOS)

C        Move to center location of tic label and write x-value:

         IF (ILABEL .EQ. 0) GO TO 19
         IF (I .EQ. 1  .AND.  IEND .EQ.-1) GO TO 19
         IF (I .EQ. 1  .AND.  IEND .EQ. 2) GO TO 19
         IF (I .EQ. N  .AND.  IEND .GT. 0) GO TO 19
         IF (I .GT. 1  .AND.  I .LT. N  .AND. IEND .EQ. -2) GO TO 19

         CALL MOVEAB (XORIG + DIST, YPOS)

         YLABL = YMIN + REAL (I-1) * YINT
         IF (ABS (ILABEL) .EQ. 2) THEN
             IYLABL = INT (YLABL)
             WRITE (FLABEL, '(I5)') IYLABL
             NCHAR = 5
         ELSE
             WRITE (FLABEL, '(F7.2)') YLABL
             NCHAR = 7
         END IF

         CALL CSTRING (FLABEL, ' ', CHASIZ, NCHAR, 90., 1)

   19    YPOS = YPOS + DY
      END DO

   30 CONTINUE

C     Axis title:

      IF (NTITLE .GT. 0) THEN
         CALL MOVEAB (XORIG + (1.+ RELSIZ) * DIST, (YORIG + YEND) * 0.5)
         CALL CSTRING (YTITLE, FONT, CHASIZ * RELSIZ, NTITLE, 90., 1)
      END IF

      RETURN
      END
C=======================================================================
C
      SUBROUTINE CSTRING (STRING, FONT, CHASIZ, NCHAR, ANGLE, ICENTR)
C
C     CSTRING writes a text string centered about or beginning at the
C     current position.
C
C     FONT  :  String such as Times_Roman or Symbol, ignored if blank
C     CHASIZ:  Character font size
C     NCHAR :  Number of characters
C     ANGLE :  Angle of rotation from horizontal axis in degrees
C              about the current position
C     ICENTR:  Center string about current position if ICENTR = 1,
C              otherwise write string beginning at current position.
C
C=======================================================================

      COMMON /PSC/ JWRIT

      CHARACTER  STRING * (*), FONT * (*)

      IF (FONT .NE. ' ') WRITE (JWRIT, 1) FONT, CHASIZ
      IF (ANGLE .NE. 0.) WRITE (JWRIT, 2) ANGLE
      IF (ICENTR .EQ. 1) WRITE (JWRIT, 3) STRING (1:NCHAR)
      WRITE (JWRIT, 4) STRING (1:NCHAR)
      IF (ANGLE .NE. 0.) WRITE (JWRIT, 2) -ANGLE

      RETURN

    1 FORMAT ('/', A, ' findfont ',F9.4,' scalefont setfont')
    2 FORMAT (F9.4,' r')      ! rotate
    3 FORMAT ('(',A,') stringwidth pop -.5 mul 0. rm')  ! rmoveto
    4 FORMAT ('(',A,') sh')   ! show

      END
C=======================================================================
C
      SUBROUTINE CTEXT (XO, YO, STRING, NCHAR)
C
C     CTEXT centers a text string about the indicated position, using
C     the current font and character size (no rotations).
C
C     XO, YO:  Text is centered at XO, and sits on YO
C     STRING:  Text to be printed
C     NCHAR :  Number of characters
C
C     11/29/96  D.A.Saunders  Adaptation of CSTRING & SYMBOL.
C=======================================================================

      COMMON /PSC/ JWRIT

      CHARACTER  STRING * (*)

      WRITE (JWRIT, 1) XO, YO
      WRITE (JWRIT, 2) STRING (1:NCHAR)
      WRITE (JWRIT, 3) STRING (1:NCHAR)

      RETURN

    1 FORMAT (2F9.4, ' m')    ! moveto
    2 FORMAT ('(',A,') stringwidth pop -.5 mul 0. rm')  ! rmoveto
    3 FORMAT ('(',A,') sh')   ! show

      END
C=======================================================================
C
      SUBROUTINE SYMTYPE (FONT, CHASIZ)
C
C     SYMTYPE establishes a character string font and size, so that
C     SYMBOL no longer does it for every string.
C
C     FONT  :  String such as Times_Roman or Symbol, ignored if blank
C     CHASIZ:  Font size
C=======================================================================

      CHARACTER FONT * (*)
      REAL      CHASIZ

      COMMON /PSC/ JWRIT

      WRITE (JWRIT, 1) FONT, CHASIZ

      RETURN

    1 FORMAT ('/', A, ' findfont ',F9.4,' scalefont setfont')

      END
C=======================================================================
C
      SUBROUTINE SYMBOL (XO, YO, STRING, NCHAR, ANGLE)
C
C     SYMBOL writes a text string at (XO,YO), centered if NCHAR = 1,
C     with optional rotation.  Use it in conjunction with SYMTYPE.
C     ANGLE :  Angle of rotation in degrees, counterclockwise
C
C=======================================================================

      CHARACTER STRING * (*)

      COMMON /PSC/ JWRIT

      WRITE (JWRIT, 1) XO, YO
      IF (ANGLE .NE. 0.) WRITE (JWRIT, 2) ANGLE
      IF (NCHAR .EQ. 1)  WRITE (JWRIT, 3) STRING (1:1)
      WRITE (JWRIT, 4) STRING (1:NCHAR)
      IF (ANGLE .NE. 0.) WRITE (JWRIT, 2) -ANGLE

      RETURN

    1 FORMAT (2F9.4, ' m')    ! moveto
    2 FORMAT (F9.4, ' r')     ! rotate
    3 FORMAT ('(', A, ') stringwidth pop -.5 mul dup rm')  ! center
    4 FORMAT ('(', A, ') sh') ! show

      END
C=======================================================================
C
      SUBROUTINE LINETYPE (TYPE)
C
C     Line-type utility.
C
C=======================================================================

      CHARACTER TYPE * (*)

      COMMON /PSC/ JWRIT

      WRITE (JWRIT, '(A,A)') 'st ', TYPE  ! st = stroke

      RETURN
      END
C=======================================================================
C
      SUBROUTINE RGBCOLOR (COLOR)
C
C     RGB color utility (primary colors + black only).
C
C=======================================================================

      CHARACTER COLOR * (*)

      COMMON /PSC/ JWRIT

      IF (COLOR .EQ. 'black') THEN
         CALL LINETYPE (' ')  ! Force a stroke before reverting to black
         WRITE (JWRIT, '(A)') '0. setgray'  ! 1. = white
      ELSE
         RED   = 0.
         GREEN = 0.
         BLUE  = 0.
         IF      (COLOR(1:1) .EQ. 'r') THEN
            RED = 1.
         ELSE IF (COLOR(1:1) .EQ. 'b') THEN
            BLUE = 1.
         ELSE IF (COLOR(1:1) .EQ. 'g') THEN
            GREEN = 1.
         END IF
      END IF

      WRITE (JWRIT, '(3F3.0, A)') RED, BLUE, GREEN, ' setrgbcolor'

      RETURN
      END
C=======================================================================
C
      SUBROUTINE PLOT (X, Y, IPEN)
C
C     Move/draw utility. IPEN = 2 means line to; IPEN = 3 means move to.
C
C=======================================================================

      COMMON /MIC/ NEN
      COMMON /PSC/ JWRIT

      NEN = NEN + 1
      IF (IPEN .EQ. 3) THEN ! Move
         WRITE (JWRIT, 1) X, Y
      ELSE ! IPEN .EQ. 2 (Draw) is the only other possibility
         WRITE (JWRIT, 2) X, Y
         IF (NEN .GE. 500) THEN
            NEN = 0
            WRITE (JWRIT, 3) 
            WRITE (JWRIT, 1) X, Y
         END IF
      END IF

      RETURN

    1 FORMAT (2F9.4, ' m')  ! moveto
    2 FORMAT (2F9.4, ' l')  ! lineto
    3 FORMAT ('st')         ! stroke

      END
C=======================================================================
C
      SUBROUTINE MOVEAB (X, Y)
C
C=======================================================================

      COMMON /PSC/ JWRIT

      WRITE (JWRIT, 1) X, Y

      RETURN

    1 FORMAT (2F9.4, ' m')  ! moveto

      END
C=======================================================================
C
      SUBROUTINE LINEAB (X, Y)
C
C=======================================================================

      COMMON /PSC/ JWRIT

      WRITE (JWRIT, 2) X, Y

      RETURN

    2 FORMAT (2F9.4, ' l')  ! lineto

      END
C=======================================================================
C
      SUBROUTINE ENDPLT
C
C=======================================================================

      COMMON /PSC/ JWRIT

      WRITE (JWRIT, '(A)') 'stroke', 'showpage'

      RETURN
      END
