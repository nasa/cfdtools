      SUBROUTINE DTTSCL (NDIM, S, M, IER)
C
C     CREATE A TRANSFORMATION MATRIX TO SCALE AN OBJECT.
C
C
C     USAGE
C
C         DOUBLE PRECISION S(NDIM), M(NDIM+1, NDIM+1)
C         CALL DTTSCL (NDIM, S, M, IER)
C
C
C     INPUT
C
C         NDIM    THE DIMENSION OF THE OBJECTS TO BE SCALED.
C                 NDIM = 2 FOR PLANAR, OR NDIM = 3 FOR 3-D SPACE.
C
C         S       THE SCALE AMOUNTS:  SCALE BY A FACTOR OF S(1)
C                 ALONG THE X-AXIS; SCALE BY A FACTOR OF S(2) ALONG
C                 THE Y-AXIS; AND (IF NDIM = 3) SCALE BY A FACTOR OF
C                 S(3) ALONG THE Z-AXIS.
C
C
C     OUTPUT
C
C         M       THE TRANSFORMATION MATRIX IN HOMOGENEOUS COORDINATES.
C
C                         SX  0   0   0
C                 M =     0   SY  0   0
C                         0   0   SZ  0
C                         0   0   0   1
C
C                 TO SCALE THE POINT (X,Y,Z), MULTIPLY
C
C                     SX  0   0   0          X
C                     0   SY  0   0    .     Y
C                     0   0   SZ  0          Z
C                     0   0   0   1          1
C
C                          SX * X
C                 TO GET   SY * Y
C                          SZ * Z
C                             1
C
C         IER     SUCCESS/ERROR CODE.
C                 FOR IER < 0, DTTSCL HAS SET M(1) TO DTMCON (1).
C
C                 IER =  0    NO ERRORS DETECTED
C
C                 IER = -1    NDIM < 2 OR NDIM > 3
C
C
C     WRITTEN BY
C         DEBORAH PARSONS
C         MAY 2, 1989
C
      EXTERNAL DTMCON
      INTEGER I, J, NDIM, IER
      DOUBLE PRECISION S(*), M(NDIM+1, *), DTMCON
C
      CHARACTER*8 SUBNAM
      DATA SUBNAM /'DTTSCL  '/
C
C     ERROR CHECKING
C
      IER = 0
      IF ((NDIM .LT. 2) .OR. (NDIM .GT. 3)) THEN
          IER = -1
          GOTO 9900
      ENDIF
C
C     CREATE THE SCALING MATRIX
C
C             SX  0   0   0
C     M =     0   SY  0   0
C             0   0   SZ  0
C             0   0   0   1
C
      DO 10 I = 1, NDIM+1
          DO 10 J = 1, NDIM+1
              IF (I .EQ. J) THEN
                  IF (I .LE. NDIM) THEN
                      M(I,J) = S(I)
                  ELSE
                      M(I,J) = 1.0
                  ENDIF
              ELSE
                  M(I,J) = 0.0
              ENDIF
   10     CONTINUE
C
 9900 IF (IER .LT. 0) THEN
          M(1,1) = DTMCON(1)
          CALL DTERR (1, SUBNAM, IER, 0)
      ENDIF
      RETURN
      END
