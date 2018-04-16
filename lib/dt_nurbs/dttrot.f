      SUBROUTINE DTTROT (NDIM, R, M, IER)
C
C     CREATE A TRANSFORMATION MATRIX TO ROTATE AN OBJECT.
C
C
C     USAGE
C
C         DOUBLE PRECISION R(NDIM), M(NDIM+1, NDIM+1)
C         CALL DTTROT (NDIM, R, M, IER)
C
C
C     INPUT
C
C         NDIM    THE DIMENSION OF THE OBJECTS TO BE ROTATED.
C                 NDIM = 2 FOR PLANAR, OR NDIM = 3 FOR 3-D SPACE.
C
C         R       THE ROTATION ANGLES (DEGREES).  ALL ROTATIONS ARE
C                 COUNTER-CLOCKWISE.
C
C                 IF NDIM = 3, R(1) IS THE ANGLE OF ROTATION ABOUT THE
C                 X-AXIS; R(2) IS THE ANGLE OF ROTATION ABOUT THE
C                 Y-AXIS; AND R(3) IS THE ANGLE OF ROTATION ABOUT THE
C                 Z-AXIS.
C
C                 IF NDIM = 2, ONLY R(1) IS USED: R(1) IS THE ANGLE OF
C                 ROTATION ABOUT THE ORIGIN.
C
C
C     OUTPUT
C
C         M       THE TRANSFORMATION MATRIX IN HOMOGENEOUS COORDINATES.
C
C                         1           0           0           0
C                 RX =    0           COS R(1)   -SIN R(1)    0
C                         0           SIN R(1)    COS R(1)    0
C                         0           0           0           1
C
C                         COS R(2)    0           SIN R(2)    0
C                 RY =    0           1           0           0
C                        -SIN R(2)    0           COS R(2)    0
C                         0           0           0           1
C
C                         COS R(3)   -SIN R(3)    0           0
C                 RZ =    SIN R(3)    COS R(3)    0           0
C                         0           0           1           0
C                         0           0           0           1
C
C                 M = RZ * RY * RX
C
C                (SORT OF REVERSED, SO THAT M * P COMES OUT TO BE
C                 [RZ*[RY*[RX*P]]] AND THE ROTATION IS ABOUT THE
C                 X-AXIS, THEN THE Y-AXIS, THEN THE Z-AXIS, AS
C                 ADVERTISED.)
C
C
C                 TO ROTATE THE POINT P = (X,Y,Z), MULTIPLY
C
C                                X
C                     M    .     Y
C                                Z
C                                1
C
C                          X '
C                 TO GET   Y '
C                          Z '
C                           1
C
C         IER     SUCCESS/ERROR CODE.
C                 FOR IER < 0, DTTROT HAS SET M(1) TO DTMCON (1).
C
C                 IER =  0    NO ERRORS DETECTED
C
C                 IER = -1    NDIM < 2 OR NDIM > 3
C
C
C     WRITTEN BY
C         DEBORAH PARSONS
C         MAY 4, 1989
C
      EXTERNAL DTMCON, DTMMPS
      INTEGER I, J, NDIM, IER
      DOUBLE PRECISION R(*), M(NDIM+1,*), DTMCON
      DOUBLE PRECISION T1(4,4), T2(4,4), D2R
C
      CHARACTER*8 SUBNAM
      DATA SUBNAM /'DTTROT  '/
C
C     ERROR CHECKING
C
      IER = 0
      IF ((NDIM .LT. 2) .OR. (NDIM .GT. 3)) THEN
          IER = -1
          GOTO 9900
      ENDIF
C
C     INITIALIZE DEGREES TO RADIANS CONVERSION
C
      D2R = DTMCON(15)
C
C     INITIALIZE THE ROTATION MATRIX TO THE IDENTITY MATRIX
C
      DO 10 I = 1, NDIM+1
          DO 10 J = 1, NDIM + 1
              IF (I .EQ. J) THEN
                  M(I,J) = 1.0
              ELSE
                  M(I,J) = 0.0
              ENDIF
              T1(I,J) = M(I,J)
   10 CONTINUE
C
C     CREATE THE ROTATION MATRIX
C
      IF (NDIM .EQ. 2) THEN
C
C     ---  TWO-DIMENSIONAL (PLANAR)  ---
C
          M(1,1) =  COS(D2R*R(1))
          M(1,2) = -SIN(D2R*R(1))
          M(2,1) =  SIN(D2R*R(1))
          M(2,2) =  COS(D2R*R(1))
      ELSE
C
C     ---  THREE-DIMENSIONAL  ---
C
C
C     MULTIPLY M BY RZ (THE ROTATION ABOUT THE Z-AXIS)
C
C             COS R(3)   -SIN R(3)    0           0
C     RZ =    SIN R(3)    COS R(3)    0           0
C             0           0           1           0
C             0           0           0           1
C
          IF (R(3) .NE. 0.0) THEN
              M(1,1) =  COS(D2R*R(3))
              M(1,2) = -SIN(D2R*R(3))
              M(2,1) =  SIN(D2R*R(3))
              M(2,2) =  COS(D2R*R(3))
              DO 20 I = 1, 4
                  DO 20 J = 1, 4
                      T1(I,J) = M(I,J)
   20         CONTINUE
          ENDIF
C
C     MULTIPLY M BY T2 = RY (THE ROTATION ABOUT THE Y-AXIS)
C
C             COS R(2)    0           SIN R(2)    0
C     RY =    0           1           0           0
C            -SIN R(2)    0           COS R(2)    0
C             0           0           0           1
C
          IF (R(2) .NE. 0.0) THEN
              DO 30 I = 1, 4
                  DO 30 J = 1, NDIM + 1
                      IF (I .EQ. J) THEN
                          T2(I,J) = 1.0
                      ELSE
                          T2(I,J) = 0.0
                      ENDIF
   30         CONTINUE
              T2(1,1) =  COS(D2R*R(2))
              T2(1,3) =  SIN(D2R*R(2))
              T2(3,1) = -SIN(D2R*R(2))
              T2(3,3) =  COS(D2R*R(2))
              CALL DTMMPS (1, 4, 4, 4, T1, 4, T2, 4, M, 4, IER)
              DO 40 I = 1, 4
                  DO 40 J = 1, 4
                      T1(I,J) = M(I,J)
   40         CONTINUE
          ENDIF
C
C     MULTIPLY M BY T2 = RX (THE ROTATION ABOUT THE X-AXIS)
C
C             1           0           0           0
C     RX =    0           COS R(1)   -SIN R(1)    0
C             0           SIN R(1)    COS R(1)    0
C             0           0           0           1
C
              DO 50 I = 1, 4
                  DO 50 J = 1, NDIM + 1
                      IF (I .EQ. J) THEN
                          T2(I,J) = 1.0
                      ELSE
                          T2(I,J) = 0.0
                      ENDIF
   50         CONTINUE
          IF (R(1) .NE. 0.0) THEN
              T2(2,2) =  COS(D2R*R(1))
              T2(2,3) = -SIN(D2R*R(1))
              T2(3,2) =  SIN(D2R*R(1))
              T2(3,3) =  COS(D2R*R(1))
              CALL DTMMPS (1, 4, 4, 4, T1, 4, T2, 4, M, 4, IER)
          ENDIF
      ENDIF
C
 9900 IF (IER .LT. 0) THEN
          M(1,1) = DTMCON(1)
          CALL DTERR (1, SUBNAM, IER, 0)
      ENDIF
      RETURN
      END
