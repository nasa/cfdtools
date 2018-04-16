      SUBROUTINE DTOSL1(NRNG,KORD,NT,T,NU,U,
     1                  CPTIN,INC1,NDIM1,ICC,HOLD,IHOLD,
     2                  CPTOUT,INC2,NDIM2,IER)
C+
C  PURPOSE:  COMPUTES B-SPLINE COEFFICIENTS OF A SPLINE CURVE
C            WITH RESPECT TO AN ENRICHED KNOT SET.
C
C  METHOD:  OSLO ALGORITHM
C
C  INPUTS:  NRNG    DIMENSION OF CURVE.
C           KORD    ORDER OF CURVE.
C           NT      NUMBER OF KNOTS IN INPUT CURVE,
C                   INCLUDING MULTIPLICITIES.
C           T       KNOT SET OF INPUT CURVE.
C           NU      NUMBER OF KNOTS IN NEW KNOT SET,
C                   INCLUDING MULTIPLICITIES.
C           U       NEW KNOT SET.
C                   EACH KNOT IN T(*) MUST BE INCLUDED IN U(*),
C                   WITH EQUAL OR GREATER MULTIPLICITY.
C           CPTIN   B-SPLINE COEFFICIENTS OF INPUT CURVE.
C           INC1    INCREMENT PARAMETER FOR CPTIN ARRAY,
C                   INC1 .GE. 1 (NOT CHECKED).
C           NDIM1   FIRST DIMENSION OF CPTIN ARRAY,
C                   NDIM1 .GE. NT - KORD (NOT CHECKED).
C           ICC     INITAL CALL/CONTINUE PARAMETER.
C                   MUST BE PASSED AS A VARIABLE.
C           INC2    INCREMENT PARAMETER FOR CPTOUT ARRAY,
C                   INC2 .GE. 1 (NOT CHECKED).
C           NDIM2   FIRST DIMENSION OF CPTOUT ARRAY,
C                   NDIM2 .GE. NU - KORD (NOT CHECKED).
C  WORKING STORAGE:
C
C           HOLD    DOUBLE PRECISION  ARRAY OF LENGTH .GE.
C                   KORD * (NU - KORD).
C           IHOLD   INTEGER ARRAY OF LENTH .GE. (NU - KORD).
C
C  OUTPUT:  ICC     INITIAL CALL/CONTINUE PARAMETER.
C                   IF IER .EQ. 0, ICC IS SET TO ABS OF INPUT VALUE.
C                   IF IER .NE. 0, ICC IS SET TO ZERO.
C
C           CPTOUT  B-SPLINE COEFFICIENTS OF THE OUTPUT CURVE.
C
C           IER     SUCCESS/ERROR CODE.
C                   IER = 0    SUCCESS.
C                   IER = -11  ICC = 0 ON INPUT.
C                   IER = -100 ERROR IN DTSPV1,
C                              USUALLY A PROBLEM WITH THE OLD OR NEW
C                              KNOT VECTORS.
C
C  DATE:  5-DEC-84
C         8-JAN-85  SUBROUTINE NAMES REVISED
C-
      DOUBLE PRECISION T(*),U(*)
      DOUBLE PRECISION CPTIN(INC1,NDIM1,*),CPTOUT(INC2,NDIM2,*)
      DOUBLE PRECISION HOLD(KORD,*)
      INTEGER IHOLD(*)
C
      EXTERNAL DTSPV2
C
      DOUBLE PRECISION SUM,ZERO,ONE
      DATA ZERO,ONE / 0.0D0, 1.0D0 /
C
      IF(ICC .EQ. 0) THEN
        IER = -11
        GO TO 9900
      END IF
C
      NCPOUT = NU - KORD
C
C  FOR INITIAL CALL, CALULATE IHOLD(*) AND HOLD(*,*).
C
C  IHOLD(J) = THE INDEX MU SUCH THAT T(MU) .LE. U(J) .LT. T(MU+1).
C
C  HOLD(*,J) = ALPHA(MU-KORD+1..MU,J) WHERE ALPHA(*,*) IS SUCH THAT
C                 CPTOUT(J,*) = SUM ( ALPHA(I,J)*CPTIN(I,*) ).
C
      IF(ICC .LT. 0) THEN
        MU = KORD
C
        DO 500 J=1,NCPOUT
          CALL DTSPV2(U(J),KORD,T,1,NT,MU,MU,IER)
          IF(IER .NE. 0) THEN
            IER = -100
            GO TO 9900
          END IF
          IHOLD(J) = MU
C
          HOLD(KORD,J) = ONE
          DO 300 K=2,KORD
            HOLD(KORD-K+1,J) = (T(MU+1)-U(J+K-1))/(T(MU+1)-T(MU-K+2))
     1                       * HOLD(KORD-K+2,J)
            DO 100 I=MU-K+2,MU-1
              HOLD(KORD-MU+I,J) = (U(J+K-1)-T(I))/(T(I+K-1)-T(I))
     1                          * HOLD(KORD-MU+I,J)
     2                          + (T(I+K)-U(J+K-1))/(T(I+K)-T(I+1))
     3                         * HOLD(KORD-MU+I+1,J)
  100       CONTINUE
            HOLD(KORD,J) = (U(J+K-1)-T(MU))/(T(MU+K-1)-T(MU))
     1                   * HOLD(KORD,J)
  300     CONTINUE
  500   CONTINUE
      END IF
C
C  IN CASE (ICC .LT 0) OR (ICC .GT. 0),
C  COMPUTE NEW COEFFICIENTS BY MULTIPLYING WITH SPARSE MATRIX
C  STORED IN HOLD(*,*) AND IHOLD(*).
C
      DO 900 K=1,NRNG
        DO 800 J=1,NCPOUT
          SUM = ZERO
          MU = IHOLD(J)
          DO 700 I=1,KORD
            SUM = SUM + HOLD(I,J) * CPTIN(1,MU-KORD+I,K)
  700     CONTINUE
          CPTOUT(1,J,K) = SUM
  800   CONTINUE
  900 CONTINUE
C
C  NORMAL RETURN.
C
      IER = 0
      ICC = ABS(ICC)
      RETURN
C
C  ERROR RETURN.
C
 9900 CONTINUE
      ICC = 0
      RETURN
      END
