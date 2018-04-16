C+**********************************************************************
C
      SUBROUTINE CONT3D( IDIM, JDIM, F, X, Y, Z, IS, IE, JS, JE, ACONT,
     +                   NCNDIM, XCONT, YCONT, ZCONT, NADDIM, NAD,
     +                   IADIM, IA )
C
C ACRONYM: CONTours 3-D
C          ----     - -
C
C PURPOSE:
C   CONT3D  was  adapted from  similar routines in  Versions 2 and 3 of
C   Program PLOT3D  by  Pieter G. Buning, NASA/Ames Research Center. It
C   calculates  contour lines  for the function  F  for  a single plane
C   orientation per call on a grid defined by arrays X, Y and Z. In the
C   case of  an IJ ( or XY ) orientation  this should be handled easily
C   by the calling program.  For IK and JK ( or XZ/YZ ) orientations it
C   will  probably  be necessary to transpose  the original 3D grid and
C   function arrays into the required 2D arrays.   CONT3D  handles just
C   one contour level per call to minimize scratch storage requirements.
C   In the event that scratch  or output  array  dimensions are too low,
C   partial contour information will be returned with a warning message.
C
C
C ARGUMENTS:
C     ARG   TYPE I/O/S  DIM     DESCRIPTION
C    IDIM     I    I     -      Effective row dimension of F(*,*)
C                               in the routine setting up this array.
C    JDIM     I    I     -      Effective column dimension of F(*,*)
C                               in the routine setting up this array.
C    F        R    I IDIM,JDIM  Function to be contoured.
C    X,Y,Z    R    I IDIM,JDIM  The coordinates matching F(I,J).
C    IS,IE,   I    I     -      Define a window within the F array for
C     JS,JE                     the portion to be contoured.
C    ACONT    R    I     -      Contour level to be plotted.
C    NCNDIM   I    I     -      Maximum number of points to be returned
C                               in XCONT, YCONT, ZCONT.
C    XCONT,   R    O   NCNDIM   On return hold the X,Y and Z coordinates
C    YCONT,ZCONT                for plotting contour lines.
C    NADDIM   I    I     -      Dimension of contour line pointer array.
C                               Should be large enough to accomodate all
C                               contour lines expected from a single call.
C    NAD      I    O   NADDIM   NAD(1) gives the number of contour lines.
C                               NAD(I) points to the start of the Ith
C                               contour line in arrays XCONT,YCONT,ZCONT.
C    IADIM    I    I     -      Dimension of scratch array.
C    IA       I    S   IADIM    Scratch array.
C
C EXTERNAL REFERENCES: None
C
C KNOWN MACHINE DEPENDENCIES: None
C
C STANDARDS VIOLATIONS: None.
C
C ENVIRONMENT:  VAX/VMS FORTRAN 77
C
C DEVELOPMENT HISTORY:
C     DATE   INITIALS   DESCRIPTION
C    Jul. 86   RCL      Adapted for use as RA Division Graphics Library
C                       module.
C    Mar. 92   DLH      Set NAD(1) (# of lines) to 0 before testing for
C			degenerate case (IS same as IE or JS same as JE).
C
C AUTHOR: Pieter G. Buning, NASA Ames Research Center
C
C-**********************************************************************
C
C
      INTEGER    IA(IADIM), NAD(NADDIM)
      REAL       XCONT(NCNDIM), YCONT(NCNDIM), ZCONT(NCNDIM),
     +           X(IDIM,JDIM), Y(IDIM,JDIM), Z(IDIM,JDIM), F(IDIM,JDIM)
C
C
      NAD(1) = 0
      IF (IS.EQ.IE .OR. JS.EQ.JE) RETURN
C
      NAD(1) = 1
      NLINEP = 2
C
C   Zero the marker array.
C
      DO 10 IIA= 1,IADIM
         IA(IIA)= 0
 10   CONTINUE
C
C   LNSTRT=1 means we're about to start a new line.
C
      LNSTRT= 1
      IGNEXT= 0
C
C
C   Flag points in IA where the the function increases through the contour
C   level, not including the boundaries.  This is so we have a list of at
C   least one point on each contour line that doesn't intersect a boundary.
C   IAE = the number of points where F(I,J) increases through ACONT between
C   two mesh points.
C   Set IMB = 1 when first finding F(I,J) < ACONT. Then set it back to zero
C   and set IAE = IAE+1 when F(I,J) increases through ACONT.
C
         IAE= 0
         DO 40 J= JS+1,JE-1
            IMB  = 0
            IAEND= IAE
            DO 40 I= IS,IE
               IF (F(I,J).LE.ACONT) THEN
                  IMB    = 1
               ELSE IF (IMB.EQ.1) THEN
                  IAE    = IAE+1
                  IF (IAE.GT.IADIM) THEN
C                    Warning - IA array full.
                     WRITE(6,2110) IADIM
                     GOTO 60
                  END IF
                  IA(IAE)= 1000*I+J
                  IMB    = 0
C
               ENDIF
 40         CONTINUE
C
C   Search along the boundaries for contour line starts.  IMA tells which
C   boundary of the region we're on.
C
 60      CONTINUE
         IMA = 1
         IMB = 0
         IBEG= IS-1
         JBEG= JS
C
 70      IBEG= IBEG+1
         IF (IBEG.EQ.IE)  IMA= 2
         GOTO 90
 75      JBEG= JBEG+1
         IF (JBEG.EQ.JE)  IMA= 3
         GOTO 90
 80      IBEG= IBEG-1
         IF (IBEG.EQ.IS)  IMA= 4
         GOTO 90
 85      JBEG= JBEG-1
         IF (JBEG.EQ.JS)  IMA= 5
 90      IF (F(IBEG,JBEG).GT.ACONT) GOTO 100
         IMB= 1
 95      GOTO (70,75,80,85,750),IMA
 100     IF (IMB.NE.1) GOTO 95
C
C   Got a start point.
C
         IMB= 0
         I  = IBEG
         J  = JBEG
         FIJ= F(IBEG,JBEG)
C
C   Round the corner if necessary.
C
         GOTO (240,110,115,120,360),IMA
 110     IF (J.NE.JS) GOTO 280
         GOTO 240
 115     IF (I.NE.IE) GOTO 320
         GOTO 280
 120     IF (J.NE.JE) GOTO 360
         GOTO 320
C
C   This is how we start in the middle of the region, using IA, where
C   IIA points to the first non-zero element of IA.
C
 190     I      = IA(IIA)/1000
         J      = IA(IIA)-1000*I
         FIJ    = F(I,J)
         IA(IIA)= 0
         GOTO 240
C                                                                       4
C   Look different directions to see which way the contour line went: 1-|-3
C                                                                       2
 220     J   = J+1
 240     I   = I-1
         IF (I.LT.IS) GOTO 700
         IDIR= 1
         IF (F(I,J).LE.ACONT) GOTO 380
         FIJ = F(I,J)
         GOTO 280
 260     I   = I-1
 280     J   = J-1
         IF (J.LT.JS) GOTO 700
         IDIR= 2
         IF (F(I,J).LE.ACONT) GOTO 500
         FIJ = F(I,J)
         GOTO 320
 300     J   = J-1
 320     I   = I+1
         IF (I.GT.IE) GOTO 700
         IDIR= 3
         IF (F(I,J).LE.ACONT) GOTO 500
         FIJ = F(I,J)
         GOTO 360
 340     I   = I+1
 360     J   = J+1
         IDIR= 4
         IF (J.GT.JE) GOTO 700
         IF (F(I,J).LE.ACONT) GOTO 500
         FIJ = F(I,J)
         GOTO 240
C
C   Wipe this point out of IA if it's in the list.
C
 380     IF (IAE.EQ.0) GOTO 500
         IJ= 1000*I+J+1000
         DO 400 IIA= 1,IAE
            IF (IA(IIA).EQ.IJ) THEN
               IA(IIA)= 0
               GOTO 500
            ENDIF
 400     CONTINUE
C
C   Do interpolation for X,Y,Z coordinates.
C
 500     XYF= (ACONT-F(I,J))/(FIJ-F(I,J))
C
C   This tests for a contour point coinciding with a grid point.  In this case
C   the contour routine comes up with the same physical coordinate twice.  If
C   If we don't trap it, it can (in some cases significantly) increase the
C   number of points in a contour line.  Also, if this happens on the first
C   point in a line, the second point could be misinterpreted as the end of a
C   (circling) contour line.
C
         IF (XYF.EQ.0.) IGNEXT= IGNEXT+1

         GOTO (610,620,630,640),IDIR

 610     WXX= X(I,J)+XYF*(X(I+1,J)-X(I,J))
         WYY= Y(I,J)+XYF*(Y(I+1,J)-Y(I,J))
         WZZ= Z(I,J)+XYF*(Z(I+1,J)-Z(I,J))
         GOTO 650

 620     WXX= X(I,J)+XYF*(X(I,J+1)-X(I,J))
         WYY= Y(I,J)+XYF*(Y(I,J+1)-Y(I,J))
         WZZ= Z(I,J)+XYF*(Z(I,J+1)-Z(I,J))
         GOTO 650

 630     WXX= X(I,J)+XYF*(X(I-1,J)-X(I,J))
         WYY= Y(I,J)+XYF*(Y(I-1,J)-Y(I,J))
         WZZ= Z(I,J)+XYF*(Z(I-1,J)-Z(I,J))
         GOTO 650

 640     WXX= X(I,J)+XYF*(X(I,J-1)-X(I,J))
         WYY= Y(I,J)+XYF*(Y(I,J-1)-Y(I,J))
         WZZ= Z(I,J)+XYF*(Z(I,J-1)-Z(I,J))

 650     CONTINUE

C   Figure out what to do with this point.
C
 655     IF (LNSTRT.NE.1) GOTO 660
C
C   This is the first point in a contour line.
C
         NP                = 1
         NAD(NLINEP)       = NAD(NLINEP-1)
         XCONT(NAD(NLINEP))= WXX
         YCONT(NAD(NLINEP))= WYY
         ZCONT(NAD(NLINEP))= WZZ
         WX    = WXX
         WY    = WYY
         WZ    = WZZ
         LNSTRT= 0
         GOTO 670
C
C   Add a point to this line.  Check for duplicate point first.
C
 660     CONTINUE
         IF (IGNEXT.EQ.2) THEN
            IF (WXX.EQ.XCONT(NP) .AND. WYY.EQ.YCONT(NP)
     +         .AND. WZZ.EQ.ZCONT(NP) ) THEN
               IGNEXT= 0
               GOTO 670
            ELSE
               IGNEXT= 1
            ENDIF
         ENDIF
C
         NAD(NLINEP) = NAD(NLINEP) + 1
         IF ( NAD(NLINEP).GT.NCNDIM ) THEN
C           Warning - CONT arrays full.
            WRITE(6,2210) NCNDIM
            GOTO 1100
         END IF
         NP = NP + 1
         XCONT(NAD(NLINEP))= WXX
         YCONT(NAD(NLINEP))= WYY
         ZCONT(NAD(NLINEP))= WZZ
C
C
C   Check to see if we're back to the initial point.
C
         IF (WXX.NE.WX) GOTO 670
         IF (WYY.NE.WY) GOTO 670
         IF (WZZ.EQ.WZ) GOTO 700
C
C   Nope.  Search for the next point on this line.
C
 670     CONTINUE
         GOTO (340,220,260,300),IDIR
C
C   This is the end of a contour line, which is reached by going to the
C   boundary of the mesh or by completing a closed curve interior contour.
C   After this we'll start a new line.
C
 700     LNSTRT= 1
         IGNEXT= 0
         NAD(NLINEP) = NAD(NLINEP) + 1
         IF ( NAD(NLINEP).GT.NCNDIM ) THEN
C           Warning - CONT arrays full.
            WRITE(6,2210) NCNDIM
            GOTO 1100
         END IF

         IF ( NP .GT.1 ) THEN
            NLINEP = NLINEP + 1
            IF ( NLINEP.GT.NADDIM) THEN
C              Warning - NAD array full.
               WRITE(6,2310) NADDIM
               GOTO 1100
            END IF
         END IF
C
C   If we're not done looking along the boundaries, go look there some more.
C
         IF (IMA.NE.5) GOTO 95
C
C   Otherwise, get the next start out of IA.
C
 750     CONTINUE
         IF (IAE.NE.0) THEN
            DO 800 IIA= 1,IAE
               IF (IA(IIA).NE.0) GOTO 190
 800        CONTINUE
         ENDIF
C
C   Loop back for the next contour level.
C
 1100 CONTINUE
      NAD(1)= NLINEP-2
      RETURN
C
C
 2110 FORMAT('  Warning - Scratch array IA full in contour routine ',
     C 'CONT3D.'/
     C       '  Picture may be incomplete.  Array was dimensioned ',
     C I5,'.')
 2210 FORMAT('  Warning - Contour line arrays XCONT, YCONT, ZCONT ',
     C 'full in contour routine CONT3D.'/
     C       '  Picture may be incomplete.Arrays were dimensioned ',
     C I5,'.')
 2310 FORMAT('  Warning - Contour line pointer array NAD full in ',
     C 'contour routine CONT3D.'/
     C       '  Picture may be incomplete.  Array was dimensioned ',
     C I5,'.')
      END
