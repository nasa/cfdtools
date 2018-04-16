C+----------------------------------------------------------------------
C
      SUBROUTINE CONT2D (IDIM, JDIM, F, X, Y, IS, IE, JS, JE, NDIM,
     >                   ACONT, NCNDIM, XCONT, YCONT, NADDIM, NAD,
     >                   IADIM, IA)
C
C     ACRONYM: CONTour line(s) for 2-Dimensional data
C              ----                - -
C
C     DESCRIPTION:
C
C        CONT2D was adapted from similar routines in Versions 2 and 3 of
C     Program PLOT3D by Pieter G. Buning, NASA/Ames Research Center.  It
C     calculates the contour line(s) for the given value of the function
C     F in the specified portion of a 2-D grid,  which may be singly- or
C     doubly-dimensioned.  It handles just one contour level per call to
C     minimize scratch storage requirements. If the workspace dimensions
C     are too small, partial contour information will be returned with a
C     warning message on logical unit 6.
C
C
C     ARGUMENTS:
C
C     ARG    DIM TYPE I/O/S DESCRIPTION
C     IDIM,    -    I   I   Declared dimensions of F(*,*) in the
C      JDIM                 calling routine.
C     F   IDIM,JDIM R   I   Function to be contoured.
C     X,       *    R   I   Coordinates matching F(I,J) are:
C      Y                    (  X(I), Y(J) ) if NDIM = 1;
C                           (X(I,J),Y(I,J)) if NDIM = 2.
C     IS,IE,   -    I   I   Define a window within F(*,*) for the
C      JS,JE                portion to be contoured.
C     NDIM     -    I   I   Number of dimensions of X & Y (1 or 2).
C     ACONT    -    R   I   Contour level to be traced.
C     NCNDIM   -    I   I   Maximum number of points provided for
C                           in XCONT(*) & YCONT(*).
C     XCONT, NCNDIM R   O   Packed X & Y coordinates for line segment(s)
C      YCONT                representing function value ACONT.
C     NADDIM   -    I   I   Dimension of contour line pointer array.
C                           NADDIM >= 1 + maximum number of line segments
C                           expected for one contour level.
C     NAD    NADDIM I   O   NAD(1) gives the number of line segments found.
C                           NAD(I) points to the start of the Ith segment
C                           in XCONT(*) & YCONT(*).  See usage outline.
C     IADIM    -    I   I   Dimension of scratch array.
C     IA     IADIM  I   S   Scratch array.
C
C
C     USAGE OUTLINE:
C
C        Apply CONT2D in a loop over multiple contour levels as follows:
C
C        DO LEV = 1, NLEVS
C
C           ACONT = RLEV(LEV)
C
C           CALL CONT2D (..., ACONT, NCNDIM, XCONT, YCONT, NADDIM, NAD, ...)
C
C           NLINES = NAD(1)  ! May be 0
C           NAD(1) = 1       ! So that the following works for LINE = 1
C
C           DO LINE = 1, NLINES
C              I1 = NAD(LINE)
C              NPTS = NAD(LINE+1) - I1  ! CONT2D defines NAD(NLINES+1)
C
C              CALL PLCURV (NPTS, XCONT(I1), YCONT(I1), ...)  ! Or whatever
C
C           END DO  ! Next line segment at this contour level
C
C        END DO  ! Next contour level
C
C
C     PROCEDURES:  None
C
C     ENVIRONMENT:  FORTRAN 77
C
C     HISTORY:
C
C     Jul. 86  R.Lefkowitz Adapted for use as RA Division Graphics Library
C                          module with addition of option for singly- or
C                          doubly-dimensioned grid coordinate arrays.
C     Mar. 92  D.Hermstad  Set NAD(1) (# of lines) to 0 before testing for
C                          degenerate case (IS same as IE or JS same as JE).
C     Nov. 95  D.Saunders  Some polishing for a flow solver application.
C
C     AUTHOR:  Pieter G. Buning, NASA Ames Research Center
C
C-----------------------------------------------------------------------

C     Arguments:

      INTEGER
     >   IDIM, JDIM, IS, IE, JS, JE, NDIM, NCNDIM, NADDIM, NAD(NADDIM),
     >   IADIM, IA(IADIM)
      REAL
     >   F(IDIM,JDIM), X(*), Y(*), ACONT, XCONT(NCNDIM), YCONT(NCNDIM)

C     Statement function:

      IDEX(I,J) = I + (J-1)*IDIM

C     Execution:

      NAD(1) = 0
      IF (IS .EQ. IE .OR. JS .EQ. JE) GO TO 999

      NAD(1) = 1
      NLINEP = 2

C     Zero the marker array.

      DO 10 IIA = 1, IADIM
         IA(IIA) = 0
   10 CONTINUE

C     LNSTRT = 1 means we're about to start a new line.

      LNSTRT = 1
      IGNEXT = 0

C     Flag points in IA where the the function increases through the contour
C     level, not including the boundaries.  This is so we have a list of at
C     least one point on each contour line that doesn't intersect a boundary.
C     IAE = the number of points where F(I,J) increases through ACONT between
C     two mesh points.
C     Set IMB = 1 when first finding F(I,J) < ACONT. Then set it back to zero
C     and set IAE = IAE+1 when F(I,J) increases through ACONT.

      IAE = 0
      DO 40 J = JS+1, JE-1
         IMB = 0
         DO 30 I = IS, IE
            IF (F(I,J) .LE. ACONT) THEN
               IMB = 1
            ELSE IF (IMB .EQ. 1) THEN
               IAE = IAE+1
               IF (IAE .GT. IADIM) THEN
C                 Warning - IA array full.
                  WRITE (6, 1000) 'IA', IADIM
                  GO TO 60
               END IF
               IA(IAE) = 1000*I+J
               IMB = 0
            END IF
   30    CONTINUE
   40 CONTINUE

C     Search along the boundaries for contour line starts.
C     IMA tells which boundary of the region we're on.

   60 CONTINUE
      IMA  = 1
      IMB  = 0
      IBEG = IS-1
      JBEG = JS

   70    IBEG = IBEG+1
         IF (IBEG .EQ. IE)  IMA = 2
         GO TO 90
   75    JBEG = JBEG+1
         IF (JBEG .EQ. JE)  IMA = 3
         GO TO 90
   80    IBEG = IBEG-1
         IF (IBEG .EQ. IS)  IMA = 4
         GO TO 90
   85    JBEG = JBEG-1
         IF (JBEG .EQ. JS)  IMA = 5
   90    IF (F(IBEG,JBEG) .GT. ACONT) GO TO 100
         IMB = 1

   95    GO TO (70, 75, 80, 85, 750), IMA

  100    IF (IMB .NE. 1) GO TO 95

C        Got a start point.

         IMB = 0
         I   = IBEG
         J   = JBEG
         FIJ = F(IBEG,JBEG)

C        Round the corner if necessary.

         GO TO (240, 110, 115, 120, 360), IMA

  110    IF (J .NE. JS) GO TO 280
         GO TO 240
  115    IF (I .NE. IE) GO TO 320
         GO TO 280
  120    IF (J .NE. JE) GO TO 360
         GO TO 320

C        This is how we start in the middle of the region, using IA,
C        where IIA points to the first non-zero element of IA.

  190    I = IA(IIA)/1000
         J = IA(IIA)-1000*I
         FIJ = F(I,J)
         IA(IIA) = 0
         GO TO 240
C                                                                         4
C        Look in all directions to see which way the contour line went: 1-|-3
C                                                                         2
  220    J = J+1
  240    I = I-1
         IF (I .LT. IS) GO TO 700
         IDIR = 1
         IF (F(I,J) .LE. ACONT) GO TO 380
         FIJ = F(I,J)
         GO TO 280
  260    I = I-1
  280    J = J-1
         IF (J .LT. JS) GO TO 700
         IDIR = 2
         IF (F(I,J) .LE. ACONT) GO TO 500
         FIJ = F(I,J)
         GO TO 320
  300    J = J-1
  320    I = I+1
         IF (I .GT. IE) GO TO 700
         IDIR = 3
         IF (F(I,J) .LE. ACONT) GO TO 500
         FIJ = F(I,J)
         GO TO 360
  340    I = I+1
  360    J = J+1
         IDIR = 4
         IF (J .GT. JE) GO TO 700
         IF (F(I,J) .LE. ACONT) GO TO 500
         FIJ = F(I,J)
         GO TO 240

C        Wipe this point out of IA if it's in the list.

  380    IF (IAE .EQ. 0) GO TO 500
         IJ = 1000*I+J+1000
         DO 400 IIA = 1,IAE
            IF (IA(IIA) .EQ. IJ) THEN
               IA(IIA) = 0
               GO TO 500
            END IF
  400    CONTINUE

C        Linear interpolation of (X,Y) coordinates:

  500    XYF = (ACONT - F(I,J)) / (FIJ - F(I,J))

C        This tests for a contour point coinciding with a grid point.
C        In this case the contour routine comes up with the same physical
C        coordinate twice.  If we don't trap it, it can (in some cases
C        significantly) increase the number of points in a contour line.
C        Also, if this happens on the first point in a line, the second point
C        could be misinterpreted as the end of a (circling) contour line.

         IF (XYF .EQ. 0.) IGNEXT = IGNEXT+1

         IF (NDIM .EQ. 1) THEN
            GO TO (550, 560, 570, 580), IDIR

  550       WXX = X(I) + XYF*(X(I+1) - X(I))
            WYY = Y(J)
            GO TO 655
  560       WXX = X(I)
            WYY = Y(J) + XYF*(Y(J+1) - Y(J))
            GO TO 655
  570       WXX = X(I) + XYF*(X(I-1) - X(I))
            WYY = Y(J)
            GO TO 655
  580       WXX = X(I)
            WYY = Y(J) + XYF*(Y(J-1) - Y(J))
C***        GO TO 655
          ELSE
            IJ = IDEX(I,J)
            GO TO (610, 620, 630, 640), IDIR

  610       IJINC = IDEX(I+1,J)
            GO TO 650
  620       IJINC = IDEX(I,J+1)
            GO TO 650
  630       IJINC = IDEX(I-1,J)
            GO TO 650
  640       IJINC = IDEX(I,J-1)
  650       CONTINUE
            WXX = X(IJ) + XYF*(X(IJINC) - X(IJ))
            WYY = Y(IJ) + XYF*(Y(IJINC) - Y(IJ))
         END IF

C        Figure out what to do with this point.

  655    IF (LNSTRT .NE. 1) GO TO 660

C        This is the first point in a contour line.

         NP = 1
         NAD(NLINEP) = NAD(NLINEP-1)
         XCONT(NAD(NLINEP)) = WXX
         YCONT(NAD(NLINEP)) = WYY
         WX = WXX
         WY = WYY
         LNSTRT = 0
         GO TO 670

C        Add a point to this line.  Check for duplicate point first.

  660    CONTINUE
         IF (IGNEXT .EQ. 2) THEN
            IF (WXX .EQ. XCONT(NP) .AND. WYY .EQ. YCONT(NP)) THEN
               IGNEXT = 0
               GO TO 670
            ELSE
               IGNEXT = 1
            END IF
         END IF

         NAD(NLINEP) = NAD(NLINEP) + 1
         IF (NAD(NLINEP) .GT. NCNDIM) THEN
C           Warning - XCONT array full.
            WRITE (6, 1000) 'XCONT', NCNDIM
            GO TO 900
         END IF
         NP = NP + 1
         XCONT(NAD(NLINEP)) = WXX
         YCONT(NAD(NLINEP)) = WYY

C        Check to see if we're back to the initial point.

         IF (WXX .NE. WX) GO TO 670
         IF (WYY .EQ. WY) GO TO 700

C        Nope.  Search for the next point on this line.

  670    CONTINUE
         GO TO (340, 220, 260, 300), IDIR

C        This is the end of a contour line, which is reached by going to the
C        boundary of the mesh or by completing a closed curve interior contour.
C        After this we'll start a new line.

  700    LNSTRT = 1
         IGNEXT = 0
         NAD(NLINEP) = NAD(NLINEP) + 1
         IF (NAD(NLINEP) .GT. NCNDIM) THEN
C           Warning - XCONT array full.
            WRITE (6, 1000) 'XCONT', NCNDIM
            GO TO 900
         END IF

         IF (NP .GT. 1) THEN
            NLINEP = NLINEP + 1
            IF (NLINEP .GT. NADDIM) THEN
C              Warning - NAD array full.
               WRITE (6, 1000) 'NAD', NADDIM
               GO TO 900
            END IF
         END IF

C        If we're not done looking along the boundaries, check for more.

         IF (IMA .NE. 5) GO TO 95

C        Otherwise, get the next start out of IA.

  750    CONTINUE
         IF (IAE .NE. 0) THEN
            DO 800 IIA = 1, IAE
               IF (IA(IIA) .NE. 0) GO TO 190
  800       CONTINUE
         END IF

C        Loop back for the next contour level.

  900 CONTINUE
      NAD(1) = NLINEP-2

  999 RETURN

C     Formats:

 1000 FORMAT (/, ' *** Warning: Scratch array ', A, ' full in contour',
     >   ' routine CONT2D.', /,
     >   ' Picture may be incomplete.  Array dimension: ', I5)

      END
