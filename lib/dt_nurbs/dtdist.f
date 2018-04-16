C***** DTDIST *****
C   FIND (LOCALLY) MINIMUM DISTANCE BETWEEN TWO CURVES OR SURFACES
C     AND A PAIR OF POINTS AT WHICH THIS DISTANCE IS ACHIEVED.  
C     WHETHER THE LOCAL MINIMUM IS ALSO THE GLOBAL MINIMUM DEPENDS ON
C     THE GIVEN INITIAL POINTS AND THE SHAPES OF THE CURVES OR SURFACES.
C
C     STARTING WITH A GIVEN POINT ON EACH SPLINE, ALTERNATELY MOVE ONE,
C     THEN THE OTHER, POINT UNTIL NO FURTHER DECREASE IN DISTANCE CAN BE 
C     ACHIEVED.  USES GAUSS'S METHOD WITH MARQUART MODIFICATIONS TO FIND 
C     NEXT (CLOSER) POINT.  IF AN INITIAL POINT'S PARAMETER IS OUTSIDE THE
C     VALID PARAMETER RANGE, THE NEAREST VALID PARAMETER VALUE IS TAKEN.
C
C   INPUT:
C     C1      SPLINE DATA ARRAY FOR FIRST CURVE OR SURFACE.
C     C2      SPLINE DATA ARRAY FOR SECOND CURVE OR SURFACE.
C     P1      PARAMETER VALUE(S) FOR INITIAL POINT ON C1.
C     P2      PARAMETER VALUE(S) FOR INITIAL POINT ON C2.
C     EPS     CONVERGENCE EPSILONS:
C             EPS(1)  CONTROLS ABSOLUTE DISTANCE.
C             EPS(2)  CONTROLS RELATIVE CHANGE IN DISTANCE.
C             EPS(3)  CONTROLS ABSOLUTE CHANGE IN ANSWER.
C             EPS(4)  CONTROLS RELATIVE CHANGE IN ANSWER.
C     GAMMA   INITIAL MARQUART FACTOR.  IF THE ROUTINE FAILS TO CONVERGE,
C             DUE TO A SADDLE POINT, YOU MAY TRY A MUCH LARGER GAMMA.
C     FNU     GAMMA KICK FACTOR.
C     MXITER  MAXIMUM NUMBER OF ITERATIONS ALLOWED.
C     MXHALF  MAXIMUM NUMBER OF HALVINGS ALLOWED.
C     MXMARQ  MAXIMUM NUMBER OF MARQUART MODIFICATIONS ALLOWED.
C   WORKSPACE:
C     WORK    WORK ARRAY OF LENGTH NWORK.
C     NWORK   LENGTH OF WORK ARRAY.  iT IS SUFFICIENT FOR
C             NWORK >= 5*KMAX**2 + 2*KMAX + 10,
C             WHERE KMAX IS THE LARGEST ORDER FOUND AMONG ALL PARAMETERS
C             OF BOTH SPLINES.
C   OUTPUT:
C     Q1      PARAMETER VALUE(S) FOR FINAL POINT ON C1.
C     Q2      PARAMETER VALUE(S) FOR FINAL POINT ON C2.
C     D       THE DISTANCE BETWEEN THE TWO CURVES OR SURFACES AT Q1 
C             AND Q2.
C     NITER   NUMBER OF ITERATIONS USED.
C     INFO    CONVERGENCE INFORMATION FLAGS FOR LAST MOVEMENTS OF THE 
C             FINAL POINTS ON THE TWO SPLINES.
C             INFO = 0  SUCCESSFUL MOVE TO CLOSER POINT
C             INFO = 1  INITIAL DISTANCE LESS THAN EPS(1).
C             INFO = 2  DISTANCE LESS THAN EPS(1).
C             INFO = 3  RELATIVE CHANGE IN DISTANCE LESS THAN EPS(2).
C             INFO = 4  ABSOLUTE CHANGE IN ANSWER LESS THAN EPS(3).
C             INFO = 5  RELATIVE CHANGE IN ANSWER LESS THAN EPS(4).
C             INFO = 6  SPLINE PERPENDICULAR TO RAY TO BETWEEN POINTS.
C             INFO = 7  HALVING AND MARQUART LOOPS GIVE NO IMPROVEMENT
C                       (A CHANGE IN THE EPS VALUES MAY BE IN ORDER).
C             INFO = 8  SOLUTION CONSTRAINED BY BOUNDARY.
C     IER     ERROR FLAG.  IF IER < 0, Q1(1) IS SET TO DTMCON(1) (CLOBBER 
C             CONSTANT).
C     
C   CALLS:
C     DCOPY
C     DTDIS1
C     DTERR
C     DTGET
C     DTMCON
C     DTNPDR
C     DTSPDR
C*****
      SUBROUTINE DTDIST (C1, C2, P1, P2, EPS, GAMMA, FNU, MXITER, 
     +    MXHALF, MXMARQ, WORK, NWORK, Q1, Q2, D, NITER, INFO, IER)
      INTEGER MXITER, MXHALF, MXMARQ, NWORK, NITER, INFO(2), IER
      DOUBLE PRECISION C1(*), C2(*), P1(*), P2(*), EPS(4), GAMMA, FNU
      DOUBLE PRECISION WORK(NWORK), Q1(*), Q2(*), D
C
C   LOCAL VARIABLES
      DOUBLE PRECISION A(3,2), ATA(3), ATR(2), DIAG(2), PBAR(2), 
     +    PNEW(2), R(3), RBAR(3), RNEW(3), TEST(2), XBASE(3,2), 
     +    X(3), XBAR(3), XNEW(3)
      DOUBLE PRECISION PLB(2,2), PUB(2,2), P(2,2), V(4,0:1)
      DOUBLE PRECISION GAMMAD, TEMP1, TEMP2, DET, HMAX, STEP, H, HBAR, 
     +    DBAR, DNEW, DOLD
      INTEGER JER, NDIM, I, N, NMARQ, NHALF, IPAR(2), MRAW(2), 
     +    KORD(2,2), NCOEF(2,2), IDER(3), NEED, KMAX, MOV, NON
      CHARACTER*8 SUBNAM
C
      DOUBLE PRECISION DTMCON
C
      DATA SUBNAM /'DTDIST'/
      DATA IDER /1, 0, 1/
C
C   INITIALIZATION
      IER = 0
      JER = 0
      INFO(1) = 0
      INFO(2) = 0
      CALL DCOPY (6, 0.0D0, 0, XBASE, 1)
      R(3) = 0.0D0
      Q1(1) = DTMCON(1)
C
C   GET BOTH SETS OF SPLINE PARAMETERS
      CALL DTGET (C1, .TRUE., 2, IPAR(1), MRAW(1), NDIM, KORD(1,1),
     +    NCOEF(1,1), PLB(1,1), PUB(1,1), IER)
      IF (IER .NE. 0) GO TO 9900
      IF (IPAR(1) .GT. 2) GO TO 9001
      IF (NDIM .LT. 2 .OR. NDIM .GT. 3) GO TO 9002
      CALL DTGET (C2, .TRUE., 2, IPAR(2), MRAW(2), J, KORD(1,2),
     +    NCOEF(1,2), PLB(1,2), PUB(1,2), IER)
      IF (IER .NE. 0) GO TO 9900
      IF (IPAR(2) .GT. 2) GO TO 9001
      IF (J .LT. 2 .OR. J .GT. 3) GO TO 9002
      IF (NDIM .NE. J) GO TO 9007
C
C   CHECK FOR SUFFICIENT WORKSPACE
      NEED = 0
      DO 10 MOV=1,2
        IF (IPAR(MOV) .EQ. 1) THEN
          NEED = MAX (NEED, KORD(1,MOV)*(KORD(1,MOV) + 4))
        ELSE
          KMAX = MAX (KORD(1,MOV), KORD(2,MOV))
          J = ABS (MRAW(MOV))
          NEED = MAX (NEED, KMAX*KMAX*(J+1) + 2*(KMAX+J))
        END IF
  10  CONTINUE
      IF (NWORK .LT. NEED) GO TO 9006
C
C   INITIALIZE STARTING POSITION DATA
      P(1,1) = P1(1)
      IF (IPAR(1) .EQ. 2) P(2,1) = P1(2)
      P(1,2) = P2(1)
      IF (IPAR(2) .EQ. 2) P(2,2) = P2(2)
      DO 30 MOV=1,2
        DO 20 I=1,IPAR(MOV)
          P(I,MOV) = MAX ( MIN (P(I,MOV), PUB(I,MOV)), PLB(I,MOV))
  20    CONTINUE
  30  CONTINUE
      CALL DTDIS1 (C1, C2, 2, P(1,2), XBASE, NDIM, WORK, NWORK, 
     +    XBASE(1,2), R, D, IER)
      IF (IER .NE. 0) GO TO 9099
      CALL DTDIS1 (C1, C2, 1, P(1,1), XBASE, NDIM, WORK, NWORK, 
     +    XBASE(1,1), R, D, IER)
      IF (IER .NE. 0) GO TO 9099
C
C   CHECK FOR TERMINATING FIRST GUESS
      IF (D .LT. EPS(1)) THEN
        INFO(1) = 1
        INFO(2) = 1
        GO TO 9000
      END IF
C
C   START MAIN LOOP (ONE MOVEMENT EACH POINT)
      NITER = 1
 100  CONTINUE
        DO 1010 MOV = 1,2
          NON = 3 - MOV
          INFO(MOV) = 0
          GAMMAD = GAMMA
          CALL DCOPY (6, 0.0D0, 0, A, 1)
C         GET FIRST PARTIAL DERIVATIVES AT MOVING POINT
          IF (IPAR(MOV) .EQ. 1) THEN
            IF (MOV .EQ. 1) THEN
              CALL DTSPDR (P(1,MOV), 1, C1, WORK, NWORK, V, 4, IER)
            ELSE
              CALL DTSPDR (P(1,MOV), 1, C2, WORK, NWORK, V, 4, IER)
            END IF
            IF (IER .LT. 0) GO TO 9099
            DO 110 I=1,NDIM
              A(I,1) = V(I,1)
 110        CONTINUE
          ELSE
            IF (MOV .EQ. 1) THEN
              CALL DTNPDR (P(1,MOV), 1, IDER(1), C1, WORK, NWORK, 
     +            A(1,1), IER)
              IF (IER .LT. 0) GO TO 9099
              CALL DTNPDR (P(1,MOV), 1, IDER(2), C1, WORK, NWORK, 
     +            A(1,2), IER)
              IF (IER .LT. 0) GO TO 9099
            ELSE
              CALL DTNPDR (P(1,MOV), 1, IDER(1), C2, WORK, NWORK, 
     +            A(1,1), IER)
              IF (IER .LT. 0) GO TO 9099
              CALL DTNPDR (P(1,MOV), 1, IDER(2), C2, WORK, NWORK, 
     +            A(1,2), IER)
              IF (IER .LT. 0) GO TO 9099
            END IF
          END IF
C  **
C  **     FORM NORMAL EQUATIONS
C  **
          ATA(1) = A(1,1)**2+A(2,1)**2+A(3,1)**2
          ATA(2) = A(1,1)*A(1,2)+A(2,1)*A(2,2)+A(3,1)*A(3,2)
          ATA(3) = A(1,2)**2+A(2,2)**2+A(3,2)**2
C
          ATR(1) = A(1,1)*R(1)+A(2,1)*R(2)+A(3,1)*R(3)
          ATR(2) = A(1,2)*R(1)+A(2,2)*R(2)+A(3,2)*R(3)
C  **
C  **     SAVE DIAGONAL ELEMENTS OF ATA IN DIAG
C  **
          DIAG(1) = ATA(1)
          DIAG(2) = ATA(3)
          NMARQ = 0
C
C       START LOOP FOR MARQUART MODIFICATIONS
  130     CONTINUE
C  **
C  **       DO MARQUART MODIFICATION
C  **
            ATA(1) = DIAG(1)+GAMMAD*DIAG(1)
            ATA(3) = DIAG(2)+GAMMAD*DIAG(2)
C  **
C  **       SOLVE MODIFIED NORMAL EQUATIONS FOR TEST DIRECTION
C  **       STORE TEST DIRECTION IN TEST(1) AND TEST(2)
C  **
C
C         UNSTRUCTURED GOTO TARGET (RECOMPUTE AFTER BDRY CONSTRAINT
C         IMPOSED)
  140       CONTINUE
            IF (IPAR(MOV) .EQ. 2) THEN
              DET = ATA(1)*ATA(3)-ATA(2)**2
              IF ( DET .LE. 0.0D0 ) GO TO 9013
              TEST(1) = (ATA(3)*ATR(1)-ATA(2)*ATR(2))/DET
              TEST(2) = (ATA(1)*ATR(2)-ATA(2)*ATR(1))/DET
            ELSE
              TEST(1)=ATR(1)/ATA(1)
              TEST(2)=0.0D0
            ENDIF
C  **
C  **       TEST FOR A NULL PERTURBATION
C  **
            IF ( ABS(TEST(1))+ABS(TEST(2)) .EQ. 0.0D0 ) THEN
              INFO(MOV) = 6
              IF (INFO(NON) .NE. 0) GO TO 9000
              GO TO 1000
            ENDIF
C  **
C  **       FIND MAXIMUM STEP
C  **
            HMAX = 2.0D0
            J = 0
            DO 180 I = 1,IPAR(MOV)
              IF (TEST(I) .GT. 0.0D0) THEN
                IF (P(I,MOV) .NE. PLB(I,MOV)) THEN
                  HMAX = MIN (HMAX, (P(I,MOV) - PLB(I,MOV))/TEST(I))
                ELSE
                  J = I
                  TEST(I) = 0.0D0
                END IF
              ELSE IF (TEST(I) .LT. 0.0D0) THEN
                IF (P(I,MOV) .NE. PUB(I,MOV)) THEN
                  HMAX = MIN (HMAX, (P(I,MOV) - PUB(I,MOV))/TEST(I))
                ELSE
                  J = I
                  TEST(I) = 0.0D0
                END IF
              END IF
  180       CONTINUE
C  **
C  **       CHECK FOR ABSOLUTE MOVEMENT
C  **
            STEP = ABS(TEST(1))+ABS(TEST(2))
            IF ( STEP .EQ. 0.0D0 ) THEN
              INFO(MOV) = 8
              IF (INFO(NON) .NE. 0) GO TO 9000
              GOTO 1000
            ENDIF
C  **
C  **       FORCE ONE-DIMENSIONAL OPTIMIZATION IF CONSTRAINED BY BOUNDARY
C  **
            IF ( J .NE. 0 ) THEN
              ATA(2) = 0.0D0
              ATR(J) = 0.0D0
              GO TO 140
C             RECOMPUTE AFTER BOUNDARY CONSTRAINT IMPOSED (UNSTRUC. GOTO)
            END IF
C  **
C  **       COMPUTE TEST POINT AND RESIDUALS
C  **
            NHALF = 0
            H = MIN (HMAX,1.0D0)
            HBAR = H
C
C         START HALVING LOOP
  200       CONTINUE
              DO 210 I = 1, IPAR(MOV)
                PNEW(I) = MAX (MIN (P(I,MOV)-H*TEST(I), 
     +              PUB(I,MOV)), PLB(I,MOV))
  210         CONTINUE
              CALL DTDIS1 (C1, C2, MOV, PNEW, XBASE, NDIM, 
     +            WORK, NWORK, XNEW, RNEW, DNEW, IER)
              IF (IER .LT. 0) GO TO 9099
C  **
C  **         APPLY PARABOLIC INTERPOLATION TO (HOPEFULLY) IMPROVE
C  **         THE ESTIMATED PARAMETRIC LOCATION OF THE CLOSEST POINT
C  **
              TEMP1 = H*(ATR(1)*TEST(1)+ATR(2)*TEST(2))
              IF ( TEMP1 .GT. 0.0D0 ) THEN
                TEMP2 = TEMP1+TEMP1+DNEW-D
                IF ( TEMP2 .GT. 0.0D0 ) THEN
                  HBAR = MIN(H*TEMP1/TEMP2,HMAX)
                  IF ( HBAR .NE. H ) THEN
                    DO 220 I = 1, IPAR(MOV)
                      PBAR(I) = MAX (MIN (P(I,MOV)-HBAR*TEST(I), 
     +                    PUB(I,MOV)), PLB(I,MOV))
  220               CONTINUE
                    CALL DTDIS1 (C1, C2, MOV, PBAR, XBASE, NDIM, 
     +                  WORK, NWORK, XBAR, RBAR, DBAR, IER)
                    IF (IER .LT. 0) GO TO 9099
                    IF ( DBAR .LT. DNEW ) THEN
                      DNEW = DBAR
                      PNEW(1) = PBAR(1)
                      PNEW(2) = PBAR(2)
                      RNEW(1) = RBAR(1)
                      RNEW(2) = RBAR(2)
                      RNEW(3) = RBAR(3)
                      XNEW(1) = XBAR(1)
                      XNEW(2) = XBAR(2)
                      XNEW(3) = XBAR(3)
                    END IF
                  END IF
                END IF
              END IF
C
C           EXIT HALVING AND MARQUART LOOPS IF NEW POINT IS CLOSER
              IF (DNEW .LT. D) GO TO 230
C  **
C  **         SOLUTION DID NOT IMPROVE - TRY HALVING H
C  **
              H = 0.5D0 * MIN (HBAR, H)
              NHALF = NHALF+1
              IF (NHALF .LE. MXHALF ) GO TO 200
C
C           END OF HALVING LOOP
C  **
C  **       ROTATE TOWARDS STEEPEST DESCENT TO TRY TO IMPROVE (MARQUART)
C  **
            GAMMAD = FNU*GAMMAD
            NMARQ = NMARQ+1
            IF (NMARQ .LT. MXMARQ) GO TO 130
C
C         END MARQUART MODIFICATIONS LOOP
          INFO(MOV) = 7
          IF (INFO(NON) .NE. 0) GO TO 9000
          GO TO 1000
C  **
C  **     SOLUTION IS GETTING BETTER.  ACCEPT NEW POINT.
C  **
  230     CONTINUE
          INFO(MOV) = 0
          P(1,MOV) = PNEW(1)
          P(2,MOV) = PNEW(2)
          R(1) = RNEW(1)
          R(2) = RNEW(2)
          R(3) = RNEW(3)
          XBASE(1,MOV) = XNEW(1)
          XBASE(2,MOV) = XNEW(2)
          XBASE(3,MOV) = XNEW(3)
          HMAX = DIM(HMAX,H)
          DOLD = D
          D = DNEW
C  **
C  **     TEST FOR CONVERGENCE CASES 2, 3, 4, OR 5
C  **
          IF (D .LE. EPS(1)) THEN
            INFO(MOV) = 2
            INFO(NON) = 2
            GO TO 9000
          ENDIF
C         IF (DOLD-D .LE. EPS(2)*DOLD) THEN
C           INFO(MOV) = 3
C           IF (INFO(NON) .NE. 0) GO TO 9000
C         ENDIF
          IF (H*STEP .LE. EPS(3)) THEN
            INFO(MOV) = 4
            IF (INFO(NON) .NE. 0) GO TO 9000
          ENDIF
          IF (H*STEP .LE. EPS(4)*(P(1,MOV)+P(2,MOV))) THEN
            INFO(MOV) = 5
            IF (INFO(NON) .NE. 0) GO TO 9000
          ENDIF
 1000     CONTINUE
C         SWITCH TO MOVING OTHER POINT (REVERSE DIFFERENCE VECTOR)
          R(1) = - R(1)
          R(2) = - R(2)
          R(3) = - R(3)
 1010   CONTINUE
        NITER = NITER + 1
        IF (NITER .LE. MXITER) GO TO 100
C     END MAIN LOOP
      IER = 1
C   NORMAL EXIT (MORE OR LESS)
 9000 CONTINUE
      D = SQRT (D)
      Q1(1) = P(1,1)
      IF (IPAR(1) .EQ. 2) Q1(2) = P(2,1)
      Q2(1) = P(1,2)
      IF (IPAR(2) .EQ. 2) Q2(2) = P(2,2)
      IF (IER .GT. 0) CALL DTERR (0, SUBNAM, IER, 0)
      RETURN
C
C   ERROR EXITS
C   
 9001 IER = -1
      GO TO 9900
 9002 IER = -2
      GO TO 9900
 9006 IER = -6
      CALL DTERR (2, SUBNAM, IER, NEED)
      RETURN
 9007 IER = -7
      GO TO 9900
 9013 IER = -13
      GO TO 9900
 9099 IER = -99
      CALL DTERR (4, SUBNAM, IER, 0)
      RETURN
 9900 CONTINUE
      CALL DTERR (1, SUBNAM, IER, 0)
      RETURN
      END
C
C
