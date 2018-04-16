      SUBROUTINE DTSCR3 (NDIM, P1, P2, P3, IOPT, WORK, NWORK, C, IER)
C
C     GENERATE THE B-SPLINE REPRESENTATION OF A CIRCLE, GIVEN THREE
C     POINTS ON THE CIRCLE.
C
C
C     USAGE
C
C     DOUBLE PRECISION P1(NDIM), P2(NDIM), P3(NDIM), WORK(NWORK), C(MC)
C     CALL DTSCR3 (P1, P2, P3, NDIM, IOPT, WORK, NWORK, C, IER)
C
C     WHERE MC >= 53
C
C
C     INPUT
C
C     NDIM        THE DIMENSION OF THE POINTS.  NDIM = 2 OR 3.
C
C     P1,P2,P3    THREE DISTINCT POINTS WHICH LIE ON THE DESIRED CIRCLE.
C
C     IOPT        FULL/PARTIAL CIRCLE OPTION.
C
C                 IOPT = 0 :  PARTIAL CIRCLE, STARTING AT P1, THROUGH
C                             P2, AND ENDING AT P3.
C
C                 IOPT <> 0:  FULL CIRCLE, STARTING AT P1, THROUGH P2,
C                             THEN P3, AND ENDING AT P1.
C
C
C     WORKING STORAGE
C
C     WORK
C
C     NWORK
C
C
C     OUTPUT
C
C     C           THE B-SPLINE REPRESENTATION OF THE CIRCLE.
C
C     IER         SUCCESS/ERROR CODE.  IF IER < 0, THEN C(1) HAS
C                 BEEN SET TO -1.
C
C                 = 0     SUCCESS.  RESULTS RETURNED IN C.
C
C                 = -1    NDIM < 2 OR NDIM > 3.
C
C                 = -2    P1, P2, P3 NOT DISTINCT.
C
C                 = -3    INSUFFICIENT WORKING STORAGE.  THE
C                         RECOMMENDED SIZE OF NWORK IS DISPLAYED IN
C                         THE WRITTEN ERROR MESSAGE.
C
C
C     WRITTEN BY
C         DEBORAH PARSONS
C         AUGUST 16, 1989
C
C
      EXTERNAL  DTMCON, DTSCRC, DTERR, DTERPT, DDOT
C
      INTEGER NDIM, IER, NWORK, NEED, I
      DOUBLE PRECISION DDOT
      DOUBLE PRECISION P1(*), P2(*), P3(*), C(*), WORK(*)
      DOUBLE PRECISION A(3), B(3), CTR(3), N(3), ADOTB
      DOUBLE PRECISION LNASQ, LNBSQ, LNNSQ, COEF1, COEF2, DENOM
C
      CHARACTER*8 SUBNAM
      DATA SUBNAM /'DTSCR3  '/
C
C     ERROR CHECKING
C
      IER = 0
C
      IF ((NDIM .LT. 2) .OR. (NDIM .GT. 3)) THEN
          IER = -1
          GOTO 9900
      ENDIF
C
      NEED = 273
      IF (NWORK .LT. NEED) THEN
          IER = -3
          GOTO 9900
      ENDIF

C
C     FIND THE NORMAL TO THE CIRCLE:
C
C     A IS THE VECTOR FROM P1 TO P2
C     B IS THE VECTOR FROM P1 TO P3
C     N IS THE NORMAL (N = A X B)
C
      DO 10 I = 1, NDIM
          A(I) = P2(I) - P1(I)
          B(I) = P3(I) - P1(I)
   10 CONTINUE
C
      IF (NDIM .EQ. 2) THEN
          A(3) = 0.0D0
          B(3) = 0.0D0
      ENDIF
C
      N(1) = A(2)*B(3) - A(3)*B(2)
      N(2) = A(3)*B(1) - A(1)*B(3)
      N(3) = A(1)*B(2) - A(2)*B(1)
C
C     FIND THE CENTER OF THE CIRCLE
C
      LNASQ = A(1)**2 + A(2)**2 + A(3)**2
      LNBSQ = B(1)**2 + B(2)**2 + B(3)**2
      LNNSQ = N(1)**2 + N(2)**2 + N(3)**2
      IF (LNNSQ .EQ. 0.0) THEN
          IER = -2
          GOTO 9900
      ENDIF
      ADOTB = DDOT (3, A, 1, B, 1)
      COEF1 = LNBSQ*(LNASQ-ADOTB)
      COEF2 = LNASQ*(LNBSQ-ADOTB)
      DENOM = LNNSQ*2.0D0
      DO 20 I = 1, 3
          CTR(I) = P1(I) + (COEF1*A(I) + COEF2*B(I)) / DENOM
   20 CONTINUE
C
C     CALL DTSCRC TO GENERATE THE CIRCLE
C
      CALL DTERPT (0)
      IF (IOPT .EQ. 0) THEN
          CALL DTSCRC (NDIM, CTR, N, P1, P3, WORK, NWORK, C, IER)
      ELSE
          CALL DTSCRC (NDIM, CTR, N, P1, P1, WORK, NWORK, C, IER)
      ENDIF
      if (ier .ne. 0) write (6,*) 'dtscrc returned ier = ',ier
      CALL DTERPT (1)
C
C     ERROR HANDLING
C
 9900 IF (IER .NE. 0) THEN
          C(1) = -1.0D0
          IF (IER .EQ. -3) THEN
              CALL DTERR (2, SUBNAM, IER, NEED)
          ELSE
              CALL DTERR (1, SUBNAM, IER, 0)
          ENDIF
      ENDIF
C
      RETURN
      END
