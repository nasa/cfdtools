      SUBROUTINE DTRHOC (IOPT, NDIM, P0, P1, P2, RHO, S, C, IER)
C
C     CONSTRUCT A SPLINE ARRAY (C) DESCRIBING A RHO CONIC, GIVEN A
C     SHOULDER POINT OR A RHO VALUE.
C
C     GIVEN THREE UNIQUE POINTS P0, P1 AND P2, LET M BE THE MIDPOINT
C     BETWEEN P0 AND P2.  THE SHOULDER POINT S IS THE POINT ON THE
C     CONIC BETWEEN M AND P0.  RHO IS THE DESIRED RATIO OF THE DIS-
C     TANCE BETWEEN M AND S AND THE DISTANCE BETWEEN M AND P1.
C
C
C     USAGE:
C
C         DOUBLE PRECISION P0(NDIM), P1(NDIM), P2(NDIM), C(NC), S(NDIM)
C         DOUBLE PRECISION RHO
C         CALL DTRHOC (IOPT, NDIM, P0, P1, P2, RHO, S, C, IER)
C
C             WHERE NC IS 11 + 3 * (NDIM+1)
C
C     INPUT:
C
C         IOPT        .EQ. 0 : RHO GIVEN, COMPUTE SHOULDER POINT (S) AND
C                              GENERATE CONIC
C                     .NE. 0 : SHOULDER POINT (S) GIVEN, COMPUTE RHO AND
C                              GENERATE CONIC
C
C         NDIM        DIMENSION OF POINTS (2 OR 3)
C
C         P0, P1, P2  COORDINATES OF CONTROL POINTS
C
C         RHO         IF IOPT .EQ. 0, RHO IS THE RATIO
C                             (0 .LT. RHO .LT. 1)
C
C         S           IF IOPT .NE. 0, S IS THE SHOULDER POINT.
C                     S MUST LIE ON THE SEGMENT BETWEEN P1 AND THE
C                     MIDPOINT OF THE P0 P2 LINE SEGMENT.
C
C     OUTPUT:
C
C         S           IF IOPT .EQ. 0, S IS THE COMPUTED SHOULDER POINT.
C
C         RHO         IF IOPT .NE. 0, RHO IS THE COMPUTED RATIO.
C
C         C           THE SPLINE DEFINITION ARRAY.
C
C         IER         SUCCESS/ERROR CODE.  FOR IER < 0, RESULTS HAVE 
C                     NOT BEEN COMPUTED AND HSRHOC HAS SET C(1) = -1
C
C                     IER =  0    NO ERRORS DETECTED
C                     IER = -1    NDIM .LT. 2 OR NDIM .GT. 3
C                     IER = -2    RHO .LE. 0 OR RHO .GE. 1
C                     IER = -3    THE THREE CONTROL POINTS ARE NOT
C                                 UNIQUE.
C                     IER = -4    S DOES NOT LIE BETWEEN P1 AND THE
C                                 MIDPOINT BETWEEN P0 AND P2.
C
C     WRITTEN BY
C                 DEBORAH PARSONS
C                 JULY 27, 1989
C
      EXTERNAL DTMCON
      DOUBLE PRECISION P0(*), P1(*), P2(*), C(*), S(*)
      DOUBLE PRECISION RHO
      INTEGER ISUM, I
      DOUBLE PRECISION EPS, DTMCON, DIV, D1, D2, D3, M(3)
C
      CHARACTER*8 SUBNAM
      DATA SUBNAM /'DTRHOC  '/
C
      IER = 0
      EPS = DTMCON(6)
C
C     DO SOME ERROR CHECKING
C
      IF ((NDIM .LT. 2) .OR. (NDIM .GT. 3)) THEN
          IER = -1
          GOTO 9000
      ENDIF
C
      IF ((IOPT .EQ. 0) .AND.
     *    ((RHO .LE. 0.0) .OR. (RHO .GE. 1.0))) THEN
          IER = -2
          GOTO 9000
      ENDIF
C
C     CHECK TO SEE IF CONTROL POINTS ARE UNIQUE
C
      ISUM = 0
      DO 10 I = 1, NDIM
          IF ((ABS(P0(I)-P1(I)) .LE. EPS) .AND.
     *       (ABS(P1(I)-P2(I)) .LE. EPS))        ISUM = ISUM+1
   10 CONTINUE
C
C     RETURN ERROR IF ALL COORDINATES ARE THE SAME
C
      IF (ISUM .EQ. 3) THEN
          IER = -3
          GOTO 9000
      ENDIF
C
C     CHECK TO SEE WHETHER S LIES ON THE BISECTOR
C
C     FIND THE MID-POINT BETWEEN P0 AND P2
C
      IF (IOPT .NE. 0) THEN
          DO 20 I = 1, NDIM
              M(I) = (P0(I) + P2(I)) / 2.0
   20     CONTINUE
C
C         COMPUTE THE 3 DISTANCES (THESE MAKE COMPUTATION OF RHO TRIVIAL)
C
          IF (NDIM .EQ. 2) THEN
              D1 = SQRT((P1(1)-M(1))**2+(P1(2)-M(2))**2)
              D2 = SQRT((S(1)-M(1))**2+(S(2)-M(2))**2)
              D3 = SQRT((S(1)-P1(1))**2+(S(2)-P1(2))**2)
          ELSE
              D1 = SQRT((P1(1)-M(1))**2+(P1(2)-M(2))**2+(P1(3)-M(3))**2)
              D2 = SQRT((S(1)-M(1))**2+(S(2)-M(2))**2+(S(3)-M(3))**2)
              D3 = SQRT((S(1)-P1(1))**2+(S(2)-P1(2))**2+(S(3)-P1(3))**2)
          ENDIF
          IF ((ABS(D1-(D2+D3)) .GT. EPS) .OR. (ABS(D1) .LE. EPS)) THEN
              IER = -4
              GOTO 9000
          ENDIF
      ENDIF
C
C     FIND RHO, OR FIND SHOULDER POINT?
C
      IF (IOPT .EQ. 0) THEN
C
C         FIND SHOULDER POINT
C
          DIV = (1.0-RHO) / 2.0
          DO 110 I = 1, NDIM
              S(I) = DIV * P0(I) + RHO * P1(I) + DIV * P2(I)
  110     CONTINUE
      ELSE
C
C         FIND RHO
C
          RHO = D2 / D1
      ENDIF
C
C     COMPUTE THE C SPLINE ARRAY
C
      C(1)  = 1
      C(3)  = 3
      C(4)  = 3
      C(5)  = 3
      C(6)  = 0
      C(7)  = 0
      C(8)  = 0
      C(9)  = 1
      C(10) = 1
      C(11) = 1
      IF (ABS(RHO-0.5) .GT. EPS) THEN
C
C         THE C SPLINE IS RATIONAL
C
          C(2)  = -(NDIM+1)
          DO 100 I = 1, NDIM
              C(12+(I-1)*3) = (1.0-RHO) * P0(I)
              C(13+(I-1)*3) = RHO * P1(I)
              C(14+(I-1)*3) = (1.0-RHO) * P2(I)
  100     CONTINUE
          C(11+3*NDIM+1) = 1.0-RHO
          C(11+3*NDIM+2) = RHO
          C(11+3*NDIM+3) = 1.0-RHO
      ELSE
C
C         THE C SPLINE ARRAY DOES NOT NEED TO BE RATIONAL
C
          C(2)  = NDIM
          DO 200 I = 1, NDIM
              C(12+(I-1)*3) = P0(I)
              C(13+(I-1)*3) = P1(I)
              C(14+(I-1)*3) = P2(I)
  200     CONTINUE
      ENDIF
C
      RETURN
 9000 CONTINUE
      CALL DTERR (1, SUBNAM, IER, 0)
      C(1) = -1.0
      RETURN
      END
