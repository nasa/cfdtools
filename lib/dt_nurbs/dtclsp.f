      SUBROUTINE DTCLSP (C, PT, PIN, EPS, GAMMA, FNU, MXITER,
     *                   MXHALF, MXMARQ, WORK, NWORK, P, D, NITER, 
     *                   INFO, IER)
C+
C  **     COMPUTE THE POINT ON A SURFACE CLOSEST TO A GIVEN POINT.
C  **
C  **     USES GAUSS'S METHOD WITH MARQUART MODIFICATIONS
C  **
C  **     DAVE FERGUSON   9/79
C  **     DAVE DODSON     1/80
C  **     DAVE DODSON     3/80
C  **     DAVE DODSON     8/80
C  **     A.K. JONES      8/87    CONVERTED TO REAL*8
C  **     DEBORAH PARSONS 8/89    CONVERTED TO STANDARD FORM
C  **
C  **     INPUT:  C      SPLINE ARRAY DEFINING THE SPLINE OR SURFACE
C  **             PT     GIVEN POINT IN 3-DIMENSIONAL SPACE
C  **             PIN    INITIAL GUESS IN PARAMETER SPACE; LENGTH = C(1)
C  **             EPS(1)  CONTROLS ABSOLUTE DISTANCE
C  **                 2   CONTROLS RELATIVE CHANGE IN DISTANCE
C  **                 3   CONTROLS ABSOLUTE CHANGE IN ANSWER
C  **                 4   CONTROLS RELATIVE CHANGE IN ANSWER
C  **             GAMMA   INITIAL MARQUART FACTOR
C  **             FNU     GAMMA KICK FACTOR
C  **             MXITER  MAXIMUM NUMBER OF ITERATIONS ALLOWED
C  **             MXHALF  MAXIMUM NUMBER OF HALFINGS ALLOWED
C  **             MXMARQ  MAXIMUM NUMBER OF MARQUART MODIFICATIONS ALLOWED
C  **
C  **     OUTPUT: P      THE POINT ON THE SURFACE CLOSEST TO PT
C  **                    IN PARAMETER SPACE; LENGTH = C(1)
C  **             D      THE DISTANCE BETWEEN PT AND SURFACE
C  **             NITER  NUMBER OF ITERATIONS
C  **             INFO   CONVERGENCE INFORMATION FLAG.
C  **                 1  INITIAL DISTANCE LESS THAN EPS(1)
C  **                 2  DISTANCE LESS THAN EPS(1)
C  **                 3  RELATIVE CHANGE IN DISTANCE LESS THAN EPS(2)
C  **                 4  ABSOLUTE CHANGE IN ANSWER LESS THAN EPS(3)
C  **                 5  RELATIVE CHANGE IN ANSWER LESS THAN EPS(4)
C  **                 6  SURFACE PERPENDICULAR TO RAY TO PT
C  **                 7  FURTHER ITERATIONS GIVE NO IMPROVEMENT
C  **                    (A CHANGE IN THE EPS VALUES MAY BE IN ORDER)
C  **                 8  SOLUTION CONSTRAINED BY BOUNDARY
C  **             IER    ERROR FLAG.  IF IER .LT. 0, P(1) IS SET TO
C  **                    DTMCON(1)--CLOBBER CONSTANT.
C  **                 1  NO CONVERGENCE - MAXIMUM NUMBER OF ITERATIONS
C  **                 0  NO COMPUTATIONAL PROBLEMS
C  **                -1  C(1) < 1 OR C(1) > 2
C  **                -2  C(2) .NE. 3
C  **                -3  C(3) OR C(4) < 1 OR > DTMCON(11)
C  **                -4  INVALID NO. OF B-SPLINE COEFFICIENTS
C  **                -5  INVALID KNOT SEQUENCE
C  **                -6  INSUFFICIENT WORKING STORAGE
C  **               -13  ATA MATRIX NOT POSITIVE DEFINITE
C  **               -50  NEGATIVE ERROR FLAG RETURNED BY DTSPVL
C  **               -60  NEGATIVE ERROR FLAG RETURNED BY DTNPDR OR 
C  **                    DTSPDR
C  **
C-
      EXTERNAL DTCLS1, DTNPDR, DTSPDR, DTMCON, DTSCHK
      EXTERNAL DCOPY, DTERR, DTERPT
C
      INTEGER MXITER, MXHALF, MXMARQ, NWORK, IER, INFO
      INTEGER JER, NITER, NDIM, I, J, NMARQ, NHALF
      INTEGER NDIMV, IP1, IP2, IPAR, IDER(2), NEED, KMAX
      PARAMETER (NDIMV=4)
      DOUBLE PRECISION C(*), PT(*),PIN(2),P(*), WORK(*)
      DOUBLE PRECISION A(3,2),ATA(3),ATR(2),DIAG(2), EPS(4)
      DOUBLE PRECISION GAMMA, FNU, D
      DOUBLE PRECISION PBAR(2),PNEW(2),R(3),RBAR(3),RNEW(3),TEST(2)
      DOUBLE PRECISION GAMMAD, TEMP1, TEMP2, DET, HMAX, STEP, H
      DOUBLE PRECISION HBAR, DBAR, DNEW, DOLD
      DOUBLE PRECISION V(NDIMV,0:1), DTMCON
      DOUBLE PRECISION LKNOT(2), RKNOT(2)
      CHARACTER*8 SUBNAM
      SUBNAM = 'DTCLSP  '
C
C     ERROR CHECKING
C
      IER  = 0
      INFO = 0
      JER  = 0
C
      IPAR = INT(C(1))
      IF ((IPAR .LT. 1) .OR. (IPAR .GT. 2)) THEN
          IER = -1
          GOTO 9900
      ENDIF
C
      NDIM = INT(C(2))
      IF (NDIM .LT. 0) THEN
          NDIM = -1-NDIM
      ENDIF
      IF ((NDIM .LT. 2) .OR. (NDIM .GT. 3)) THEN
          IER = -2
          GOTO 9900
      ENDIF
C
      CALL DTERPT(0)
      CALL DTSCHK (C, IER)
      CALL DTERPT(1)
      IF (IER .NE. 0) THEN
          GOTO 9900
      ENDIF
C
      IF (IPAR .EQ. 1) THEN
          NEED = C(3)*C(3) + 4*C(3)
      ELSE
          KMAX = MAX(C(3),C(4))
          NEED = KMAX**2*ABS(C(2)) + 2*KMAX + KMAX**2 + 2*ABS(C(2))
      ENDIF
      IF (NWORK .LT. NEED) THEN
          IER = -6
          GOTO 9900
      ENDIF
C
C     EXTRACT THE KNOTS FROM THE SPLINE ARRAY
C
      IP1 = 3*C(1)+2
      IP2 = IP1+C(3)+C(INT(C(1))+3)
      DO 5 I = 1, IPAR
          LKNOT(I) = C(IP1+1)
          RKNOT(I) = C(IP2)
          IP1 = IP2
          IP2 = IP1+C(4)+C(INT(C(1))+4)
    5 CONTINUE
C  **
C  **     INITIALIZE
C  **
      GAMMAD = GAMMA
      NITER = 1
      DO 10 I = 1, IPAR
          P(I) = MAX(MIN(PIN(I),RKNOT(I)),LKNOT(I))
   10 CONTINUE
C  **
C  **     COMPUTE INITIAL RESIDUAL
C  **
      CALL DTCLS1 (C, P, PT, NDIM, WORK, NWORK, R, JER)
      IF ( JER .LT. 0 ) THEN
          IER = -50
          GOTO 9900
      ENDIF
      D = R(1)**2+R(2)**2+R(3)**2
C  **
C  **     TEST FOR INITIAL CONVERGENCE
C  **
      IF ( D .LT. EPS(1) ) THEN
          INFO = 1
          GOTO 9900
      ENDIF
C  **
C  **     MAIN LOOP
C  **
  100 CONTINUE
C  **
C  **     FIND THE PARTIAL DERIVATIVES
C  **
      IF (IPAR .EQ. 1) THEN
          CALL DTERPT(0)
          CALL DTSPDR (P, 1, C, WORK, NWORK, V, NDIMV, JER)
          CALL DTERPT(1)
          IF (JER .NE. 0) THEN
              IER = -60
              GOTO 9900
          ENDIF
          DO 110 I = 1, NDIM
              A(I,1) = V(I,1)
              A(I,2) = 0.0D0
  110     CONTINUE
      ELSE
          IDER(1) = 1
          IDER(2) = 0
          CALL DCOPY (6, 0.0D0, 0, A(1,1), 1)
          CALL DTERPT(0)
          CALL DTNPDR (P, 1, IDER, C, WORK, NWORK, A(1,1), JER)
          CALL DTERPT(1)
          IF (JER .NE. 0) THEN
              IER = -60
              GOTO 9900
          ENDIF
          IDER(1) = 0
          IDER(2) = 1
          CALL DTERPT(0)
          CALL DTNPDR (P, 1, IDER, C, WORK, NWORK, A(1,2), JER)
          CALL DTERPT(1)
          IF (JER .NE. 0) THEN
              IER = -60
              GOTO 9900
          ENDIF
      ENDIF
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
C  **
C  **     DO MARQUART MODIFICATION
C  **
  130 CONTINUE
      ATA(1) = DIAG(1)+GAMMAD*DIAG(1)
      ATA(3) = DIAG(2)+GAMMAD*DIAG(2)
C  **
C  **     SOLVE MODIFIED NORMAL EQUATIONS FOR TEST DIRECTION
C  **     STORE TEST DIRECTION IN TEST(1) AND TEST(2)
C  **
  140 CONTINUE
      IF (IPAR .EQ. 2) THEN
          DET = ATA(1)*ATA(3)-ATA(2)**2
          IF ( DET .LE. 0.0D0 ) THEN
              IER = -13
              GOTO 9900
          ENDIF
          TEST(1) = (ATA(3)*ATR(1)-ATA(2)*ATR(2))/DET
          TEST(2) = (ATA(1)*ATR(2)-ATA(2)*ATR(1))/DET
      ELSE
          TEST(1)=ATR(1)/ATA(1)
          TEST(2)=0.0D0
      ENDIF
C  **
C  **     TEST FOR A NULL PERTURBATION
C  **
      IF ( ABS(TEST(1))+ABS(TEST(2)) .EQ. 0.0D0 ) THEN
          INFO = 6
          GOTO 9900
      ENDIF
C  **
C  **     FIND MAXIMUM STEP
C  **
      HMAX = 2.0D0
      J = 0
      DO 180 I = 1,IPAR
        IF ( TEST(I) .EQ. 0.0D0 ) GO TO 180
          IF ( TEST(I) .LT. 0.0D0 ) GO TO 160
            IF ( P(I) .NE. LKNOT(I) ) GO TO 150
              J = I
              P(I) = LKNOT(I)
              TEST(I) = 0.0D0
              GO TO 180
  150       HMAX = MIN((P(I)-LKNOT(I))/TEST(I),HMAX)
            GO TO 180
  160     IF ( P(I) .NE. RKNOT(I) ) GO TO 170
            J = I
            P(I) = RKNOT(I)
            TEST(I) = 0.0D0
            GO TO 180
  170     HMAX = MIN((P(I)-RKNOT(I))/TEST(I),HMAX)
  180   CONTINUE
C  **
C  **     CHECK FOR ABSOLUTE MOVEMENT
C  **
      STEP = ABS(TEST(1))+ABS(TEST(2))
      IF ( STEP .EQ. 0.0D0 ) THEN
          INFO = 8
          GOTO 9900
      ENDIF
C  **
C  **     FORCE ONE-DIMENSIONAL OPTIMIZATION IF CONSTRAINED BY BOUNDARY
C  **
      IF ( J .EQ. 0 ) GO TO 190
      ATA(2) = 0.0D0
      ATR(J) = 0.0D0
      GO TO 140
C  **
C  **     COMPUTE TEST POINT AND RESIDUALS
C  **
  190 NHALF = 0
      H = MIN(HMAX,1.0D0)
      HBAR = H
  200 DO 205 I = 1, IPAR
          PNEW(I) = MAX(MIN(P(I)-H*TEST(I),RKNOT(I)),LKNOT(I))
  205 CONTINUE
      CALL DTCLS1 (C, PNEW, PT, NDIM, WORK, NWORK, RNEW, JER)
      IF ( JER .LT. 0 ) THEN
          IER = -50
          GOTO 9900
      ENDIF
      DNEW = RNEW(1)**2+RNEW(2)**2+RNEW(3)**2
C  **
C  **     APPLY PARABOLIC INTERPOLATION TO (HOPEFULLY) IMPROVE
C  **     THE ESTIMATED PARAMETRIC LOCATION OF THE CLOSEST POINT
C  **
      TEMP1 = H*(ATR(1)*TEST(1)+ATR(2)*TEST(2))
      IF ( TEMP1 .LE. 0.0D0 ) GO TO 210
      TEMP2 = TEMP1+TEMP1+DNEW-D
      IF ( TEMP2 .LE. 0.0D0 ) GO TO 210
      HBAR = MIN(H*TEMP1/TEMP2,HMAX)
      IF ( HBAR .EQ. H ) GO TO 210
      DO 206 I = 1, IPAR
          PBAR(I) = MAX(MIN(P(I)-HBAR*TEST(I),RKNOT(I)),LKNOT(I))
  206 CONTINUE
      CALL DTCLS1 (C, PBAR, PT, NDIM, WORK, NWORK, RBAR, JER)
      IF ( JER .LT. 0 ) THEN
          IER = -50
          GOTO 9900
      ENDIF
      DBAR = RBAR(1)**2+RBAR(2)**2+RBAR(3)**2
      IF ( DBAR .GE. DNEW ) GO TO 210
      DNEW = DBAR
      PNEW(1) = PBAR(1)
      PNEW(2) = PBAR(2)
      RNEW(1) = RBAR(1)
      RNEW(2) = RBAR(2)
      RNEW(3) = RBAR(3)
C  **
C  **     TEST FOR IMPROVEMENT
C  **
  210 IF ( DNEW .LT. D ) GO TO 230
C  **
C  **     SOLUTION DID NOT IMPROVE - TRY HALVING H
C  **
      IF ( NHALF .EQ. 0 ) GAMMAD = FNU*GAMMAD
      NHALF = NHALF+1
      IF ( NHALF .GT. MXHALF ) GO TO 220
      H = 0.5*MIN(HBAR,H)
      GO TO 200
C  **
C  **     ROTATE TOWARDS STEEPEST DESCENT TO TRY TO IMPROVE
C  **
  220 NMARQ = NMARQ+1
      IF ( NMARQ .LT. MXMARQ ) GO TO 130
      INFO = 7
      GOTO 9900
C  **
C  **     SOLUTION IS GETTING BETTER
C  **
  230 P(1) = PNEW(1)
      IF (IPAR .EQ. 2) P(2) = PNEW(2)
      R(1) = RNEW(1)
      R(2) = RNEW(2)
      R(3) = RNEW(3)
      HMAX = DIM(HMAX,H)
      DOLD = D
      D = DNEW
C  **
C  **     TEST FOR CONVERGENCE CASES 2, 3, 4, OR 5
C  **
      IF ( D .LE. EPS(1) ) THEN
          INFO = 2
          GOTO 9900
      ENDIF
C      IF ( DOLD-D .LE. EPS(2)*DOLD ) THEN
C          INFO = 3
C          GOTO 9900
C      ENDIF
      IF ( H*STEP .LE. EPS(3) ) THEN
          INFO = 4
          GOTO 9900
      ENDIF
      IF ( H*STEP .LE. EPS(4)*(ABS(P(1))+ABS(P(IPAR))) ) THEN
          INFO = 5
          GOTO 9900
      ENDIF
C  **
C  **     CONVERGENCE HAS NOT OCCURRED
C  **     TRY TO ROTATE BACK TOWARD GAUSS'S METHOD
C  **
      IF ( NHALF .EQ. 0  .AND.  HBAR .GT. 1.5*H ) GAMMAD = GAMMAD/FNU
C  **
C  **     ITERATE MXITER TIMES
C  **
      NITER = NITER+1
      IF ( NITER .LE. MXITER ) GO TO 100
      IER = 1
      GOTO 9900
C  **
C  **     
C  **     PRINT ERROR MESSAGES AND QUIT
C  **
C  **
 9900 CONTINUE
      D = SQRT(D)
      IF (IER .LT. 0) THEN
          P(1) = DTMCON(1)
          IF (IER .EQ. -6) THEN
              CALL DTERR (2, SUBNAM, IER, NEED)
          ELSE
              CALL DTERR (1, SUBNAM, IER, 0)
          ENDIF
      ELSE IF (IER .GT. 0) THEN
          CALL DTERR (0, SUBNAM, IER, 0)
      ENDIF
C  **
      RETURN
      END

      SUBROUTINE DTCLS1 (C, P, PT, NDIM, WORK, NWORK, RDIFF, JER)
C
C     FIND THE DIFFERENCE BETWEEN THE CURRENT POINT (PARAMETRIC 
C     FORM) AND THE DESIRED POINT (3-DIMENSIONAL COORDINATES).
C
      EXTERNAL DTNPVL
C
      DOUBLE PRECISION C(*), P(*), PT(*), RDIFF(*), WORK(*)
      INTEGER I, JER, NWORK, NDIM
C
      CALL DCOPY (3, 0.0D0, 0, RDIFF(1), 1)
      CALL DTERPT(0)
      CALL DTNPVL (P, 1, C, WORK, NWORK, RDIFF, JER)
      CALL DTERPT(1)
C
      DO 10 I = 1, NDIM
          RDIFF(I) = RDIFF(I) - PT(I)
   10 CONTINUE
C
      END

