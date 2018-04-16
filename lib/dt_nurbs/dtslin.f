      SUBROUTINE DTSLIN (NDIM, P, HLEN, V, C, IER)
C
C     GENERATE A LINE IN SPLINE (C-ARRAY) FORM, GIVEN A CENTER POINT,
C     HALF-LENGTH, AND DIRECTION VECTOR (SLOPE).
C
C
C     USAGE
C
C         DOUBLE PRECISION P(NDIM), HLEN, V(NDIM), C(NC)
C         CALL DTSLIN (NDIM, P, HLEN, V, C, IER)
C
C         WHERE NC = 13 IF NDIM = 2
C            OR NC = 15 IF NDIM = 3
C
C
C     INPUT
C
C         NDIM    THE DIMENSION OF THE LINE.  NDIM = 2 (PLANAR) OR
C                 NDIM = 3 (3-SPACE).
C
C         P       THE CENTER POINT OF THE LINE.
C
C         HLEN    THE HALF-LENGTH OF THE LINE.  THE RESULTING LINE
C                 WILL HAVE LENGTH = 2 * HLEN.
C
C         V       A VECTOR DESCRIBING THE SLOPE (DIRECTION) OF THE
C                 LINE.  V MUST HAVE AT LEAST ONE NON-ZERO ELEMENT.
C
C
C     OUTPUT
C
C         C       THE SPLINE ARRAY.
C
C         IER     SUCCESS/ERROR CODE.
C                 FOR IER < 0, DTSLIN HAS SET C(1) = -1.
C
C                 IER =  0    NO ERRORS DETECTED.
C
C                 IER = -1    NDIM < 2 OR NDIM > 3.
C
C                 IER = -2    HLEN .LE. 0.
C
C                 IER = -3    V(I) = 0 FOR ALL I, I=1..NDIM.
C
C
C     WRITTEN BY
C         DEBORAH PARSONS
C         MAY 10, 1989
C
      INTEGER I, NDIM, IER
      DOUBLE PRECISION P(*), HLEN, V(*), C(*), X(2), Y(2), Z(2)
C
      CHARACTER*8 SUBNAM
      DATA SUBNAM /'DTSLIN  '/
C
C     ERROR CHECKING
C
      IER = 0
      IF ((NDIM .LT. 2) .OR. (NDIM .GT. 3)) THEN
          IER = -1
          GOTO 9900
      ENDIF
C
      IF (HLEN .LE. 0.0) THEN
          IER = -2
          GOTO 9900
      ENDIF
C
      IF (V(1) .NE. 0.0) THEN
C
C     SOLVE FOR X, THEN Y AND Z IN TERMS OF X, IF POSSIBLE
C
C         X = P(1) +/- HLEN/(SQRT((V(2)/V(1))**2 + (V(3)/V(1))**2 + 1))
C         Y = (V(2)/V(1)) * (X - P(1)) + P(2)
C         Z = (V(3)/V(1)) * (X - P(1)) + P(3)
C
          X(1) = P(1) - HLEN/(SQRT((V(2)/V(1))**2 + (V(3)/V(1))**2 + 1))
          X(2) = P(1) + HLEN/(SQRT((V(2)/V(1))**2 + (V(3)/V(1))**2 + 1))
          DO 10 I = 1, 2
              Y(I) = (V(2)/V(1)) * (X(I) - P(1)) + P(2)
              IF (NDIM .EQ. 3) THEN
                  Z(I) = (V(3)/V(1)) * (X(I) - P(1)) + P(3)
              ENDIF
   10     CONTINUE
      ELSE IF (V(2) .NE. 0.0) THEN
C
C     SOLVE FOR Y, THEN X AND Z IN TERMS OF Y, IF POSSIBLE
C
C         Y = P(2) +/- HLEN/(SQRT((V(1)/V(2))**2 + (V(3)/V(2))**2 + 1))
C         X = (V(1)/V(2)) * (Y - P(2)) + P(1)
C         Z = (V(3)/V(2)) * (Y - P(2)) + P(3)
C
          Y(1) = P(2) - HLEN/(SQRT((V(1)/V(2))**2 + (V(3)/V(2))**2 + 1))
          Y(2) = P(2) + HLEN/(SQRT((V(1)/V(2))**2 + (V(3)/V(2))**2 + 1))
          DO 20 I = 1, 2
              X(I) = (V(1)/V(2)) * (Y(I) - P(2)) + P(1)
              IF (NDIM .EQ. 3) THEN
                  Z(I) = (V(3)/V(2)) * (Y(I) - P(2)) + P(3)
              ENDIF
   20     CONTINUE
      ELSE IF ((NDIM .EQ. 3) .AND. (V(3) .NE. 0.0)) THEN
C
C     SOLVE FOR Z, THEN X AND Y IN TERMS OF Z, IF POSSIBLE
C
C         Z = P(3) +/- HLEN/(SQRT((V(1)/V(3))**2 + (V(2)/V(3))**2 + 1))
C         X = (V(1)/V(3)) * (Z - P(3)) + P(1)
C         Y = (V(2)/V(3)) * (Z - P(3)) + P(2)
C
          Z(1) = P(3) - HLEN/(SQRT((V(1)/V(3))**2 + (V(2)/V(3))**2 + 1))
          Z(2) = P(3) + HLEN/(SQRT((V(1)/V(3))**2 + (V(2)/V(3))**2 + 1))
          DO 30 I = 1, 2
              X(I) = (V(1)/V(3)) * (Z(I) - P(3)) + P(1)
              Y(I) = (V(2)/V(3)) * (Z(I) - P(3)) + P(2)
   30     CONTINUE
      ELSE
C
C     ERROR -- ALL ELEMENTS OF V ARE ZERO
C
          IER = -3
          GOTO 9900
      ENDIF
C
C     LOAD THE C ARRAY
C
      C(1) = 1.0
      C(2) = NDIM
      C(3) = 2.0
      C(4) = 2.0
      C(5) = 2.0
      C(6) = 0.0
      C(7) = 0.0
      C(8) = 1.0
      C(9) = 1.0
      C(10) = X(1)
      C(11) = X(2)
      C(12) = Y(1)
      C(13) = Y(2)
      IF (NDIM .EQ. 3) THEN
          C(14) = Z(1)
          C(15) = Z(2)
      ENDIF
C
 9900 IF (IER .LT. 0) THEN
          C(1) = -1.0
          CALL DTERR (1, SUBNAM, IER, 0)
      ENDIF
      RETURN
      END
