C+----------------------------------------------------------------------
C
      SUBROUTINE PLTGEO (DRAW, MXGPTS, XGEOM, YGEOM, NGFIG, NGPTS,
     >                   GCLIP, GCOLOR, GLINES, GTHICK, HEIGHT, WIDTH,
     >                   MXGWRK, XWORK, YWORK, LUNERR)

C  PURPOSE:
C     PLTGEO draws multiple curves on a plot with optional blanking of
C     the region interior or exterior to each.  The blanking is handled
C     by DISSPLA software, with the requirement that PLTGEO be called
C     after all the plot axis drawing and scaling are complete, but
C     before any of the curves to be tested for clipping are drawn.
C
C     Since the boundaries of the (optionally-clipped) regions may
C     well need to be drawn LAST so as not to be obscured by other
C     curves such as contours or grid lines, while setting up for
C     clipping must be done FIRST, we have a conflict.  Solution:
C     CALL PLTGEO first with DRAW = .FALSE., meaning just set up for
C     clipping (if any); then call PLTGEO last with DRAW = .TRUE.,
C     meaning do the drawing, but not the clipping set-up.  In fact,
C     suppress any clipping before doing the drawing, just to be sure
C     of getting thick lines properly (educated guess).
C
C  ARGUMENTS:
C   VAR      DIM     TYPE    I/O/S   DESCRIPTION
C  DRAW       -       L      I       DRAW=F means set up blanking (only, if any)
C                                    DRAW=T means draw outlines (only).
C                                    Normally, PLTGEO must be called twice.
C  MXGPTS     -       I      I       Maximum total number of geom. coordinates.
C  XGEOM    MXGPTS    R      I       Array containing abscissas for geometry(s).
C  YGEOM    MXGPTS    R      I       Array containing ordinates for geometry(s).
C  NGFIG      -       I      I       Number of geometries to be processed.
C  NGPTS    NGFIG     I      I       Number of data points per geometry.
C  GCLIP    NGFIG     C      I       Clipping type for geometry:
C                                    'NOCLIP', 'CLIPIN', or 'CLIPOUT'.
C  GCOLOR   NGFIG     R      I       Color index for each geometry.
C                                      0.   ...  Current default
C                                      1.   ...  Black
C                                      2.   ...  Magenta
C                                      3.   ...  Red
C                                      4.   ...  Yellow
C                                      5.   ...  Green
C                                      6.   ...  Cyan
C                                      7.   ...  Blue
C                                      8.   ...  White
C                                    with colors between the six primaries
C                                    specified using real numbers in [2.,7.]
C                                    (E.g. orange = 3.5, yellow-green = 4.5,
C                                    etc.)
C  GLINES   NGFIG     I      I       Line type code - see POLYLINE.
C  GTHICK   NGFIG     I      I       Multiplier of default geometry line
C                                    thickness.
C  HEIGHT     -       R      I       Length of vertical axis in inches.
C  WIDTH      -       R      I       Length of horizontal axis in inches.
C  MXGWRK     -       I      I       Dimension of X and Y work arrays.
C  XWORK,  MXGWRK     R      I/S     Work arrays used to define blanked areas,
C  YWORK                             where MXGWRK >= 2*MXGPTS + 8.
C  LUNERR     -       I      I       Logical unit for error messages.
C
C  PROCEDURES:
C   POLYLINE Draws a curve of specified line type, thickness, and color.
C   PLPOLY   Draws outline of a closed curve and blanks interior.
C   XPOSN    DISSPLA routine that converts data units to inches.
C   YPOSN          "                 "                  "
C
C  ENVIRONMENT:
C   VAX/VMS, FORTRAN 77, CA-DISSPLA V11-9003
C   IMPLICIT NONE is non-standard.
C
C  HISTORY:
C    9 Sep. 86   RCL   Original design and implementation.
C    6 Dec. 89   MDW   Replaced call to FMCURV with call to generic
C                      curve plotting routine POLYLINE.
C   17 Jan. 90   MDW   Some minor details.  END DO's eliminated.
C   13 Jan. 91   DAS   COLOR was INTEGER by mistake.  Also FRAME.
C                      Made drawing and setup for blanking mutually
C                      exclusive - always need separate calls now.
C                      Other refinements.
C   28 Dec. 91   DAS   MXGWRK size requirement was described wrongly.
C                      Why does the FIRST element of a geometry being
C                      clipped "out" have be "left"most?  (I.e. why
C                      the clumsy in-place shift if it's not?)
C
C  AUTHOR:    Rosalie Lefkowitz, Sterling Software, Palo Alto, CA
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C  *  Arguments:

      INTEGER
     >   NGFIG, GLINES (NGFIG), GTHICK (NGFIG), LUNERR, MXGPTS, MXGWRK, 
     >   NGPTS (NGFIG)
      REAL
     >   GCOLOR (NGFIG), HEIGHT, WIDTH, XGEOM (MXGPTS), XWORK (MXGWRK),
     >   YGEOM (MXGPTS), YWORK (MXGWRK)
      CHARACTER
     >   GCLIP (NGFIG) * (*)
      LOGICAL
     >   DRAW

C  *  Local constants:

      INTEGER
     >   MXGFIG
      PARAMETER
     >  (MXGFIG = 10)  ! Some reasonable number needed for local blanking IDs

C  *  Local variables:

      INTEGER
     >   I, ID (MXGFIG), IER, IFIRST, IG, ILAST, NPTS, OFFSET
      REAL
     >   COLOR, FRAME
      SAVE
     >   ID

C  *  Procedures:

      EXTERNAL
     >   XPOSN, YPOSN
      REAL
     >   XPOSN, YPOSN

C  *  Execution:


      IF (NGFIG .GT. MXGFIG) GO TO 920
      IF (MXGWRK .LT. 2*MXGPTS + 8) GO TO 940

      FRAME = 0.   ! No frame
      OFFSET = 0

      DO 500, IG = 1, NGFIG

C  *     Copy the curve into work arrays:

         IF (IG .GT. 1) OFFSET = OFFSET + NGPTS (IG-1)
         NPTS = NGPTS (IG)
         DO 100, I = 1, NPTS
            XWORK (I) = XGEOM (I+OFFSET)
            YWORK (I) = YGEOM (I+OFFSET)
 100     CONTINUE
      
         IF (DRAW) THEN

C  *        Draw the outline (only).
C           Suppress blanking in case it helps thick outlines:

            CALL BLOFF (ID (IG))     ! ID set up on a prior call to PLTGEO.

            COLOR = GCOLOR (IG)
            IF (GCOLOR (IG) .EQ. 8) COLOR = 0.

            CALL POLYLINE (NPTS, XWORK, YWORK, GLINES (IG),
     >                     GTHICK (IG), -1, COLOR, IER)

            GO TO 500

         END IF

C  *     Else we're done if no clipping has been requested:

         IF (GCLIP (IG) .EQ. 'NOCLIP') GO TO 500

C  *     Convert coordinates to inches for blanking software:

         DO 200, I = 1, NPTS
            XWORK (I) = XPOSN (XGEOM (I+OFFSET), YGEOM (I+OFFSET))
            YWORK (I) = YPOSN (XGEOM (I+OFFSET), YGEOM (I+OFFSET))
  200    CONTINUE

         IF (GCLIP (IG) .EQ. 'CLIPIN') THEN

             CALL BLPOLY (XWORK, YWORK, NPTS, FRAME)

          ELSE IF (GCLIP (IG) .EQ. 'CLIPOUT') THEN

C  *         Extend the "inch" arrays to the plot boundaries with a
C            connection at the leftmost point. The calling program
C            must provide the required space for this extension.
C            <I don't understand this part.  DAS, 12/28/91>

             IFIRST = 1
             DO 300, I = 2, NPTS
                IF (XWORK (I) .LT. XWORK (IFIRST)) IFIRST = I
  300        CONTINUE
             IF (IFIRST .EQ. 1) THEN
                ILAST = NPTS
             ELSE
                ILAST = IFIRST - 1
                DO 400, I = 1, ILAST
                   XWORK (NPTS + I) = XWORK (I)
                   YWORK (NPTS + I) = YWORK (I)
  400           CONTINUE
                ILAST = ILAST + NPTS
             END IF

             XWORK (ILAST+1) = 0.0
             XWORK (ILAST+2) = 0.0
             XWORK (ILAST+3) = WIDTH
             XWORK (ILAST+4) = WIDTH
             XWORK (ILAST+5) = 0.0
             XWORK (ILAST+6) = 0.0
             XWORK (ILAST+7) = XWORK (ILAST)
             YWORK (ILAST+1) = YWORK (ILAST)
             YWORK (ILAST+2) = HEIGHT
             YWORK (ILAST+3) = HEIGHT
             YWORK (ILAST+4) = 0.0
             YWORK (ILAST+5) = 0.0
             YWORK (ILAST+6) = YWORK (ILAST)
             YWORK (ILAST+7) = YWORK (ILAST)

             NPTS = NPTS + 7     

C  *         Set up blanked region:

             CALL BLPOLY (XWORK (IFIRST), YWORK (IFIRST), NPTS, FRAME)

C  *         Get an identifier and then turn off blanking temporarily:

             CALL BLKEY (ID (IG))
             CALL BLOFF (ID (IG))

         END IF

 500  CONTINUE

C  *  Turn all blanking on again:

      DO 600, IG = 1, NGFIG
         CALL BLON (ID (IG))
 600  CONTINUE

      GO TO 999


C  *  Error handling:

 920  WRITE (LUNERR, '(1X, A)')
     >   'PLTGEO: Abnormal termination - too many geometries requested.'
      GO TO 990

 940  WRITE (LUNERR, '(1X, A)')
     >   'PLTGEO: Abnormal termination - too little workspace provided.'
 990  STOP ' '

 999  RETURN
      END         
