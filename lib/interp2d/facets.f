C+---------------------------------------------------------------------
C
      SUBROUTINE FACETS ( XD, YD, ZD, ND, SCALED, XI, YI, ZI, NI,
     +                    IWK, IER )
C
C  PURPOSE: FACETS performs bivariate interpolation within irregularly
C           spaced data by a form of linear interpolation.  It is most
C           appropriate for modest-size datasets where the function is
C           not particularly smooth, meaning interpolation with a more
C           elaborate model may behave erratically.   (This  algorithm
C           cannot give function values outside the range of the given
C           function values,  in the same sense that linear interpola-
C           tion is "safe" in one dimension.)
C
C  METHOD:  For each given point requiring interpolation,  the nearest
C           three data points which are not essentially collinear, and
C           which surround the given point, are determined by exhaust-
C           ive search.  These define a facet or plane, which is eval-
C           uated at the given point.
C
C           More precisely:
C
C           > Scale one of the coordinates of the data points (X, say)
C             so that the data range is square, and "distances" become
C             meaningful.   (The scaling may be suppressed if the user
C             prefers to do his own scaling.)
C
C           > For each given point, find the 3 nearest points. This is
C             done by maintaining a list of points (0 or more) to skip
C             in the next search.  (Is there a better way?)
C
C           > Reject the third point if it does not include the target
C             point in the triangle it forms with the first two.  Just
C             add it to the list of points to ignore,  locate the next
C             closest point, check for enclosure, and so on.   (Reject
C             it also if the points are essentially collinear,   since
C             the formula for evaluating the plane through the  points
C             breaks down in this case.)
C
C           > Eventually, either three enclosing points are found else
C             the target point must have been outside the data range -
C             handled indirectly by a test for skipping ALL points  in
C             a search, rather than a more expensive test for  out-of-
C             range at the outset because this won't happen in  normal
C             usage.
C
C           The method may fail if the data points  are not distinct -
C           no attempt is made to detect this case.
C
C  PARAMETERS:
C       ARG    TYPE  I/O/S   DIM     DESCRIPTION
C        XD     R      I      ND     Input X locations of data points.
C        YD     R      I      ND     Input Y locations of data points.
C        ZD     R      I      ND     Input function values at data points.
C        ND     I      I      -      Number of data points.
C      SCALED   L      I      -      SCALED = T means do not scale and
C                                    unscale here; SCALED = F means the
C                                    X values are scaled to the Y range
C                                    so that distances are meaningful.
C                                    The Xs are unscaled before returning.
C        XI     R      I      NI     X locations of points at which to 
C                                    approximate.
C        YI     R      I      NI     Y locations of points at which to 
C                                    approximate.
C        ZI     R      O      NI     Approximated values at the points
C                                    (XI(I),YI(I)).  Values of 9999. mean
C                                    failure - see IER below.
C        NI     I      I      -      Number of interpolation points.
C        IWK    I      S      ND     Used for list of points to ignore.
C        IER    I      O      -      Error return code. 0 means no error.
C                                    IER = 99 means a target point was
C                                    outside the data range (or all data
C                                    points are essentially collinear) -
C                                    ZI(*) is set to 9999. in these cases.
C
C  PARAMETER CONSTANTS:
C    TOL      Tolerance for determining collinearity of 3 points -
C             basically a limit on the area of the triangle they make.
C
C  EXTERNAL REFERENCES:
C    NEARPT   Identifies nearest point to given point, with option to
C             ignore a list of points.
C
C  ENVIRONMENT: VAX/VMS FORTRAN
C
C  AUTHOR:  David Saunders, Sterling Software/Informatics, Palo Alto, CA
C
C  HISTORY:
C     02/13/86    DAS    Original design and coding (using much the same
C                        argument list as HARDY, TPSPLN, and 2-D Akima).
C     02/14/86    DAS    Nearest three points on their own do not work
C                        well - forcing enclosure makes the difference.
C     03/10/86    DAS    Set ZI(*) to 9999. for points where algorithm
C                        breaks down, but keep going.
C     03/21/86    DAS    Made scaling of the Xs optional - SMOOTH2D
C                        needs to do the scaling for other algorithms.
C
C-----------------------------------------------------------------------
C
      IMPLICIT  NONE

C ... Arguments:

      INTEGER   IER, ND, NI, IWK(ND)
      REAL      XD(ND), YD(ND), ZD(ND), XI(NI), YI(NI), ZI(NI)
      LOGICAL   SCALED

C ... Local constants:

      REAL      TOL
      PARAMETER ( TOL = 1.E-6 )

C ... Local variables:

      INTEGER   I, I1, I2, I3, NSKIP
      REAL      DET, DX1, DX2, DX3, DY1, DY2, DY3, DZ2, DZ3,
     +          SCALE, TOLER, T1, T2, T3, XMIN, XMAX, YMIN, YMAX,
     +          X0, Y0

C ... Execution:

      IER = 0

      IF ( .NOT. SCALED ) THEN

C ...    Scale the data so that distances make sense:

         XMIN = 1.E+36
         YMIN = XMIN
         XMAX = -XMIN
         YMAX = XMAX

         DO 200 I = 1, ND
            XMIN = MIN ( XD(I), XMIN )
            XMAX = MAX ( XD(I), XMAX )
            YMIN = MIN ( YD(I), YMIN )
            YMAX = MAX ( YD(I), YMAX )
  200    CONTINUE

         SCALE = ( YMAX - YMIN ) / ( XMAX - XMIN )
         TOLER = TOL * MAX ( ABS ( YMIN ), ABS ( YMAX) )

         DO 300 I = 1, ND
            XD(I) = XD(I) * SCALE
  300    CONTINUE

         DO 350 I = 1, NI
            XI(I) = XI(I) * SCALE
  350    CONTINUE

      END IF

C ... Handle each point to be interpolated one at a time:

      DO 600 I = 1, NI

         X0 = XI(I)
         Y0 = YI(I)

C ...    Locate nearest data point:

         NSKIP = 0

         CALL NEARPT ( X0, Y0, ND, XD, YD, NSKIP, IWK, I1 )

C ...    Locate second nearest point:

         NSKIP = 1
         IWK(1) = I1

         CALL NEARPT ( X0, Y0, ND, XD, YD, NSKIP, IWK, I2 )

C ...    Locate third nearest point:

         NSKIP = 2
         IWK(2) = I2

  400    CALL NEARPT ( X0, Y0, ND, XD, YD, NSKIP, IWK, I3 )

         IF ( I3 .EQ. 0 ) THEN
C  ...      NEARPT found all points were to be skipped - failure.

            ZI(I) = 9999.
            IER = 99
            GO TO 600
         END IF

C ...    I3 will have to be rejected if it is collinear with I1, I2:

         DX1 = X0  - XD(I1)
         DX2 = XD(I2) - XD(I1)
         DX3 = XD(I3) - XD(I1)
         DY1 = Y0  - YD(I1)
         DY2 = YD(I2) - YD(I1)
         DY3 = YD(I3) - YD(I1)
         DZ2 = ZD(I2) - ZD(I1)
         DZ3 = ZD(I3) - ZD(I1)

         DET = XD(I1)*YD(I2) - YD(I1)*XD(I2) + DX2*YD(I3) - DY2*XD(I3)

C ...    Reject point I3 also if it does not enclose target point.
C        The following is derived from the formulas in FUNCTION INMESH,
C        which deals with quadrilaterals rather than triangles:

         T1 = DX1 * DY2 - DY1 * DX2
         T2 = ( X0 - XD(I2) ) * ( YD(I3) - YD(I2) ) -
     +        ( Y0 - YD(I2) ) * ( XD(I3) - XD(I2) )
         T3 = ( Y0 - YD(I3) ) * DX3 - ( X0 - XD(I3) ) * DY3
         T3 = SIGN ( 1.E+0, T3 )
         T2 = SIGN ( 1.E+0, T2 )
         T1 = SIGN ( 1.E+0, T1 )

         IF ( ABS ( DET ) .GT. TOLER  .AND.
     +        T1 .EQ. T2  .AND.  T2. EQ. T3 ) THEN

C ...       Evaluate the plane facet defined by points I1, I2, I3:

            ZI(I) = ZD(I1) +
     +        ( DY1*( DX2*DZ3 - DX3*DZ2 ) - DX1*( DY2*DZ3 - DY3*DZ2 ) )/
     +        ( DX2*DY3 - DX3*DY2 )

         ELSE
C ...       Add point I3 to the list of points to skip, and try again:

            NSKIP = NSKIP + 1
            IWK(NSKIP) = I3
            GO TO 400
         END IF

  600 CONTINUE

      IF ( .NOT. SCALED ) THEN

C ...    Unscale the scaled coordinates:

         SCALE = 1.E+0 / SCALE
         DO 700 I = 1, ND
            XD(I) = XD(I) * SCALE
  700    CONTINUE

         DO 750 I = 1, NI
            XI(I) = XI(I) * SCALE
  750    CONTINUE

      END IF

  999 RETURN
      END
