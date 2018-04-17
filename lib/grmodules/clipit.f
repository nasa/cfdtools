C+----------------------------------------------------------------------
      LOGICAL FUNCTION CLIPIT (XG, YG, NGPTS, ICLIP, NGEOM, X, Y)
C  PURPOSE -
C    Clip concave polygon[s].  Return value is .TRUE. if point (X,Y) lies
C    inside/outside the polygon[s] formed by XG, YG., else .FALSE.  Inside/
C    outside testing is controlled by array ICLIP.
C
C  PARAMETERS -
C    ARG  TYPE I/O/S  DIM         DESCRIPTION
C   XG     R     I    (1)      Coordinates of polygon[s].  Start point
C   YG     R     I    (1)      should not be specified as endpoint.
C   NGPTS  I     I  (NGEOM)    Numbers of points in polygons
C   ICLIP  I     I  (NGEOM)    Inside/outside control code for each polygon
C				 2 = clip outside polygon, quit on point outside
C				 1 = clip outside polygon, but continue
C				 0 = ignore this polygon
C				-1 = clip inside  polygon, but continue
C				-2 = clip inside  polygon, quit on point inside
C   NGEOM  I     I     -       Number of polygons in XG, YG
C   X      R     I     -       Coordinates of point to test
C   Y      R     I     -
C				EXAMPLE - to show all points inside either of
C				 two squares:
C				 define first  square in (XG(I),YG(I),I=1,4)
C				 define second square in (XG(I),YG(I),I=5,8)
C				 NGPTS(1)= 4   NGPTS(2)= 4
C				 ICLIP(1)= 1   ICLIP(2)= 1
C				 NGEOM= 2
C				 test various X,Y, displaying if CLIPIT= .TRUE.
C
C  COMMONS USED -
C   /CLIPC/	Values computed on first call to CLIPIT, saved for later calls
C    VAR  TYPE I/O/S  DIM   OFFSET  DESCRIPTION
C   IFLAG  I    I/O    -            Flag for whether to load this common block
C				     0 = compute and load
C				     1 = do not compute and load
C				     set to 1 on output
C   GXMIN  R    I/O    -            Minimum X from XG array
C   GXMAX  R    I/O    -            Maximum X from XG array
C   GYMIN  R    I/O    -            Minimum Y from YG array
C   GYMAX  R    I/O    -            Maximum Y from YG array
C   PREFX  R    I/O    -            Coordinates of a point outside of the
C   PREFY  R    I/O    -             rectangle formed by GXMIN, GXMAX, etc.
C
C  EXTERNAL REFERENCES -
C    SEGINT (in this source file)
C
C  AUTHOR - Dexter L. Hermstad		Informatics General Corporation
C		11/16/82		Palo Alto, California
C
C  REVISIONS -
C    DATE      PERSON  STATEMENT OF CHANGES
C    06/29/83	DLH	Checks for end point of a geometry equal to start
C			 point.  Also optimized code.
C    04/27/83	DLH	Corrected bypass of geometry that is for display only
C    11/16/82	DLH	Original coding
C    11/16/82	DLH	Design
C
C  NOTES -
C
C    This algorithm determines whether a point satisfies the clipping
C    condition based on the number of times a line segment from the point
C    to a point outside the polygon crosses that polygon.  If the segment
C    crosses an odd number of times, then the point is inside.  The test
C    is repeated for each polygon.
C-----------------------------------------------------------------------
C
      DIMENSION XG(1), YG(1)
      INTEGER NGPTS(NGEOM), ICLIP(NGEOM)
C     ... CLIP1 is set to true or false for each polygon in the geometry.
C     ... CLIPIT is set .FALSE. for any failure of an AND clipping
C         (ICLIP(J)= 2 or -2).
C     ... CLIPIT is set .TRUE.  for any success of an OR  clipping
C         (ICLIP(J)= 1 or -1).
C     ... TEST is set   .TRUE.  after at least one polygon is tested
C
      LOGICAL CLIP1, TEST
C     INTEGER*2 IRCARY(100)
C
      COMMON /CLIPC/ IFLAG, GXMIN, GXMAX, GYMIN, GYMAX, PREFX, PREFY
C     DATA LUNDE/3/
C
      IF (IFLAG.NE.0) GO TO 200
         IFLAG= 1		! No need to compute after this time
C        ... Get total number of points in all polygons (XG, YG)
         NGALL= 0
         DO 100 I=1,NGEOM
            NGALL= NGALL+ NGPTS(I)
  100    CONTINUE
C        ... Compute minima and maxima of data ranges
         GXMIN= XG(1)
         GYMIN= YG(1)
         GXMAX= GXMIN
         GYMAX= GYMIN
         DO 150 I=2,NGALL
            IF (XG(I).LT.GXMIN) GXMIN= XG(I)
            IF (XG(I).GT.GXMAX) GXMAX= XG(I)
            IF (YG(I).LT.GYMIN) GYMIN= YG(I)
            IF (YG(I).GT.GYMAX) GYMAX= YG(I)
  150    CONTINUE
C     ... Place the reference point well outside a rectangle containing data
      PREFX= GXMAX+ 20.
      PREFY= GYMAX+ 20.
C
  200 CONTINUE
      CLIPIT= .TRUE.		! Set true in case all polygons have ICLIP=0
      TEST= .FALSE.		! No polygons tested so far
      IPT= 1			! Point to beginning of first polygon
C
C     ... For each polygon in the array
      DO 700 J=1,NGEOM
         CLIP1= .TRUE.		! Assume point is valid
         NPTS= IABS(NGPTS(J))
C        ... Do not do any clipping if this polygon is for display only
         IF (ICLIP(J).EQ.0) GO TO 680
C        ... If start point is same as end point do not consider it
         MPTS= NPTS
         IF (XG(IPT).EQ.XG(IPT+NPTS-1) .AND. YG(IPT).EQ.YG(IPT+NPTS-1))
     1       MPTS= MPTS- 1
C        ... Now assume point not valid - set it before test of first polygon
         IF (.NOT.TEST) CLIPIT= .FALSE.
         TEST= .TRUE.
         IMARK= 0
         X3= PREFX
         Y3= PREFY
C        ... Loop back here to fudge the reference point if
C            a strange intersection is found.
  490    CONTINUE
         NCROSS= 0
C        IF (IMARK.NE.0) WRITE (LUNDE,6490) X3, Y3,
C    1      (IRCARY(L), L=1,NPTS)
         IF (IMARK.EQ.1) X3= X3+ 0.928
         IF (IMARK.EQ.2) Y3= Y3+ 0.456
         IF (IMARK.EQ.3) CLIP1= .FALSE.
         IF (IMARK.EQ.3) GO TO 580	! Give up after three tries
C        DO 495 I=1,NPTS
C           IRCARY(I)= 9
C 495    CONTINUE
C
C        ... For every segment on the polygon
C            Determine if and how a segment from the reference point to
C            the test point intersects the polygon segment.
         DO 500 I=IPT,IPT+MPTS-1
            I1= I+ 1
            IF (I.EQ.IPT+MPTS-1) I1= IPT! For last segment use first vertex
            CALL SEGINT (XG(I), YG(I), XG(I1), YG(I1), X3, Y3, X, Y,IRC)
C           IRCARY(I)= IRC
            IF (IRC.EQ.2) NCROSS= NCROSS+ 1! IRC of 2 means another intersection
            IF (IRC.EQ.0 .OR. IRC.EQ.2) GO TO 500
            IF (IRC.EQ.1) GO TO 550	! IRC of 1 means point on segment
C              ... Alter reference point for strange intersections
               IMARK= IMARK+ 1
               GO TO 490
  500    CONTINUE
C
C        ... If segment crossed polygon an even number of times, point is out
C        *** NOTE: Use of IAND is much faster than MOD.  If your FORTRAN
C        ***       does not support IAND, then use the MOD statement or
C        ***       whatever function is equivalent to IAND.
CC         IF (MOD(NCROSS,2).EQ.0) CLIP1= .FALSE.
         IF (IAND(NCROSS,1).EQ.0) CLIP1= .FALSE.
C        ... But reverse sense of answer if clipping inside
         IF (ICLIP(J).LT.0) CLIP1= .NOT.CLIP1
         GO TO 580
C        ... Point is on polygon
  550    CONTINUE
         CLIP1=  .TRUE.
  580    CONTINUE
C        WRITE (LUNDE,6580) X, Y, CLIP1, (IRCARY(L),L=1,MPTS)
         IF (CLIP1) CLIPIT= .TRUE.
         IF (.NOT.CLIP1.AND. IABS(ICLIP(J)).EQ.2) CLIPIT= .FALSE.
         IF (.NOT.CLIP1.AND. IABS(ICLIP(J)).EQ.2) GO TO 800
C
  680    CONTINUE
         IPT= IPT+ NPTS		! Point to next polygon
  700 CONTINUE
C
  800 CONTINUE
C     ... CLIPIT is .TRUE. if no geometry was tested (all ICLIPs=0).
C         ELSE based on results of OR clipping and AND clipping
      RETURN
C6490 FORMAT (2F8.3,' ?',<NGEOM>I2)
C6580 FORMAT (2F8.3,L2,<NGEOM>I2)
      END
C+----------------------------------------------------------------------
      SUBROUTINE SEGINT (X1, Y1, X2, Y2, X3, Y3, X4, Y4, IRC)
C  PURPOSE -
C    Determine the type of intersection of segment (X1,Y1) (X2,Y2) and
C    segment (X3,Y3) (X4,Y4).  It is assumed that (X3,Y3) does not lie
C    on segment (X1,Y1) (X2,Y2).
C
C  PARAMETERS -
C    ARG  TYPE I/O/S  DIM         DESCRIPTION
C   X1     R     I     -       X- and Y- coordinates of first point
C   Y1     R     I     -
C   X2     R     I     -       X- and Y- coordinates of second point
C   Y2     R     I     -
C   X3     R     I     -       X- and Y- coordinates of third point
C   Y3     R     I     -
C   X4     R     I     -       X- and Y- coordinates of fourth point
C   Y4     R     I     -
C   IRC    I     O     -       Return code
C				0 = segment P1P2 does not intersect segment P3P4
C				1 = point P4 on segment P1P2
C				2 = intersection at point (not endpoint) on P1P2
C				3 = intersection at P1
C				4 = intersection at P2
C				5 = segment P1P2 lies on segment P3P4
C
C  COMMONS USED - (None)
C
C  EXTERNAL REFERENCES - (None)
C
C  AUTHOR - Dexter L. Hermstad		Informatics General Corporation
C		11/16/82		Palo Alto, California
C
C  REVISIONS -
C    DATE      PERSON  STATEMENT OF CHANGES
C    05/10/89	DLH	Epsilon changed from 1.E-3 to 1.E-6*(dist(p1,p2))
C    11/16/82	DLH	Original coding
C    11/16/82	DLH	Design
C
C  NOTES -
C-----------------------------------------------------------------------
C
      LOGICAL CL
      DATA EPSNRM/1.E-6/
C     ... Statement function to determine closeness
      CL(V1, V2)= ABS(V2-V1).LT.EPS
C
      DY= Y2- Y1
      DX= X2- X1
C     ... Rotate both segments so that segment P1P2 points in +X direction
      IF (.NOT.(DY.NE.0.)) GO TO 300
C        ... Find angle - use 90 or 270 degrees if P1P2 has infinite slope
         IF (DX.NE.0.) ANGLE= ATAN2(DY, DX)
         IF (DX.EQ.0.) ANGLE= 3.1415926/2.
         IF (DX.EQ.0. .AND. DY.LT.0.) ANGLE= -ANGLE
         COSANG= COS(ANGLE)
         SINANG= SIN(ANGLE)
         TX1= X1* COSANG+ Y1* SINANG
         TY1= -X1* SINANG+ Y1* COSANG
         TX2= X2* COSANG+ Y2* SINANG
         TY2= -X2* SINANG+ Y2* COSANG
         TX3= X3* COSANG+ Y3* SINANG
         TY3= -X3* SINANG+ Y3* COSANG
         TX4= X4* COSANG+ Y4* SINANG
         TY4= -X4* SINANG+ Y4* COSANG
         GO TO 400
C     ... Simple rotation of either 0 or 180 degrees if P1P2 has 0 slope
  300 CONTINUE
         SIGN= 1.
         IF (DX.LT.0.) SIGN= -SIGN
         TX1= SIGN* X1
         TY1= SIGN* Y1
         TX2= SIGN* X2
         TY2= SIGN* Y2
         TX3= SIGN* X3
         TY3= SIGN* Y3
         TX4= SIGN* X4
         TY4= SIGN* Y4
  400 CONTINUE
C     ... END rotate both segments so that segment P1P2 points in +X direction
C
C     ... Translate both segments so that segment P1P2 lies on X- axis
      IF (TY1.EQ.0.) GO TO 500
         TY4= TY4- TY1
         TY3= TY3- TY1
         TY2= TY2- TY1
         TY1= 0.
  500 CONTINUE
C
C     ... Translate both segments so that segment P1P2 starts at origin
      IF (TX1.EQ.0.) GO TO 600
         TX4= TX4- TX1
         TX3= TX3- TX1
         TX2= TX2- TX1
         TX1= 0.
  600 CONTINUE
C
C     ... Tranformations have placed point P1 at (0,0), point P2 at
C         (+X,0), and segment P3P4 wherever
C
      IRC= 0
C     ... Test for obvious misses
C     ... P3P4 above P1P2?
      IF (TY3.GT.0. .AND. TY4.GT.0.) GO TO 800
C     ... P3P4 below P1P2?
      IF (TY3.LT.0. .AND. TY4.LT.0.) GO TO 800
C     ... P3P4 to right of P1P2?
      IF (TX3.GT.TX2 .AND. TX4.GT.TX2) GO TO 800
C     ... P3P4 to left of P1P2?
      IF (TX3.LT.0. .AND. TX4.LT.0.) GO TO 800
C     ... END test for obvious misses
C
C     ... Test for point on segment P1P2
C     ... Y must be near 0.
      EPS= EPSNRM* TX2
      IF (.NOT.CL(TY4, 0.)) GO TO 720
C        ... X must be in range
         IF (TX4.LT.0. .OR. TX4.GT.TX2) GO TO 710
            IRC= 1
            GO TO 800
C     ... Point not on segment P1P2
C
C     ... Test for all of segment P1P2 contained in segment P3P4
  710 CONTINUE
C     ... Endpoints must be near X-axis (know P4 is near X-axis)
      IF (.NOT.CL(TY3, 0.)) GO TO 720
C        ... Signs of X-coordinates of endpoints must be different
         IF ((TX3.LT.0.   .AND. TX4.LT.0.) .OR.
     1       (TX3.GT.TX2 .AND. TX4.GT.TX2)) GO TO 800
            IRC= 5
            GO TO 800
C     ... Segment P1P2 not contained in segment P3P4
C
C     ... Find X-intercept of line segment (it must have an intercept
C         because cases where both Y values are alike were caught in
C         either the obvious miss section or the section after 710)
  720 CONTINUE
      XINT= (TX3*TY4 - TX4*TY3)/ (TY4-TY3)
C     ... Special return codes for intersection at an endpoint of P1P2
      IF (CL(XINT, 0.)) IRC= 3
      IF (CL(XINT, TX2)) IRC= 4
C     ... If no special code set, then indicate an interior intersection
      IF (IRC.EQ.0 .AND. XINT.GT.0. .AND. XINT.LT.TX2) IRC= 2
  800 CONTINUE
      RETURN
      END
