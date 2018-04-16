      SUBROUTINE DTSPAD(C,IDOM,IWORK,NWORK,COUT,IER)
C
C      COMPUTE PARTIAL DERIVATIVE OF SPLINE C
C
      DOUBLE PRECISION C(*),COUT(*)
      INTEGER IDOM,IWORK(*),NWORK,IER
C
C      INPUT VARIABLES:
C            C      INPUT SPLINE ARRAY
C            IDOM      INDEPENDENT VARIABLE WITH RESPECT TO WHICH
C                  DERIVATIVES ARE TAKEN
C            IWORK      WORK ARRAY
C            NWORK      DIMENSION OF WORK ARRAY
C            COUT      OUTPUT SPLINE ARRAY
C            IER      SUCCESS/ERROR CODE
C  **                 =    0 SUCCESS
C  **                 =   -1 KORD TOO SMALL
C  **                 =   -3 NWORK TOO SMALL
C  **                 =   -8 INVALID XKNT ARRAY
C  **                 =  -14 IDOM OUT OF RANGE
C  **                 =  -15 TOO FEW KNOTS
C  **                 =  -51 INVALID C(1)
C  **                 =  -52 INVALID C(2)
C  **
C  **  SUBPROGRAMS REFERENCED
C  **              DCOPY    DT1SVY   DTERPT   DTERR
C
C
      INTEGER MODE,NEED
      INTEGER KCOEFS,KORD,KSTEP,NCOEF1,NCOEF2,NCOEFS,NDEP,NDOM,
     1        NKNOTS,NMIDOM,ICOUNT,IER1,MC
      INTEGER I,IKORD,INDX,J,K,K1,K2,N,NCOUT,NSKNTS,NSTART
C
      DOUBLE PRECISION DENOMN
      EXTERNAL DCOPY,DT1SVY,DTERPT,DTERR
C
      CHARACTER*8 SUBNAM
      PARAMETER (SUBNAM = 'DTSPAD')
C
C      SET UP AND ERROR CHECKING
C
      IER = 0
      MODE = 1
      NEED = 0
C
C      NDOM = NUMBER OF INDEPENDENT VARIABLES
C      NDEP = NUMBER OF DEPENDENT VARIABLES
C
      NDOM = MAX(C(1),1.0D0)
      NDEP = C(2)
      IF (NWORK .LT. 2*NDOM) THEN
          IER = -3
          MODE = 2
          NEED = 2*NDOM
          GO TO 900
      END IF
      CALL DTERPT(0)
      CALL DT1SVY(C,MC,IWORK,IWORK,IWORK,NCOEF1,IER1)
      CALL DTERPT(1)
      IF (IER1 .LT. 0) THEN
          IER = IER1
          GO TO 900
      END IF
      IF ((IDOM .GT. NDOM) .OR. (IDOM .LE. 0)) THEN
          IER = -14
          GO TO 900
      END IF
C
C      SET UP COUT ARRAY
C
C      KSTEP = THE STEP BETWEEN COEFFICIENTS FOR IDOM
C      NMIDOM = THE NUMBER COEFFICIENTS TO COMPUTE FOR EACH IDOM
C             VALUE
C      IWORK(I+NDOM) = NUMBER OF COEFFICIENTS FOR ITH VARIABLE
C
      N = 3*NDOM + 2
      KSTEP = 1
      DO 20 I = 1,IDOM-1
          IKORD = C(I+2)
          IWORK(I+NDOM) = C(I+NDOM+2)
          N = N + IKORD + IWORK(I+NDOM)
          KSTEP = KSTEP * IWORK(I+NDOM)
  20      CONTINUE
      CALL DCOPY(N,C,1,COUT,1)
C
C      NSKNTS POINTS TO THE START OF THE KNOTS FOR IDOM
C      NSTART RECORDS POSITIONS IN THE C VECTOR
C      NCOUT POINTS TO LOCATIONS IN COUT
C      KCOEFS = OLD NUMBER OF COEFFICIENTS FOR IDOM
C      KORD = NEW ORDER FOR IDOM
C
      KCOEFS = C(IDOM+NDOM+2)
      KORD = C(IDOM+2) - 1.0D0
      COUT(2+IDOM) = KORD
      NSTART = N + 2
      NSKNTS = N + 1
      NCOUT = NSKNTS
      NMIDOM = KSTEP
C
C      NKNOTS = THE NEW NUMBER OF KNOTS FOR IDOM
C
      NKNOTS = 0
      DO 40 I = 1,KCOEFS - 1
          IF (C(NSTART) .LT. C(NSTART+KORD)) THEN
              COUT(NCOUT) = C(NSTART)
              NCOUT = NCOUT + 1
              NKNOTS = NKNOTS + 1
          END IF
          NSTART = NSTART + 1
  40  CONTINUE
      DO 50 I = 1,KORD
          COUT(NCOUT) = C(NSTART)
          NCOUT = NCOUT + 1
          NKNOTS = NKNOTS + 1
          NSTART = NSTART + 1
  50  CONTINUE
      NSTART = NSTART + 1
C
C      NCOEFS = THE NEW NUMBER OF COEFFICIENTS
C
      NCOEFS = NKNOTS - KORD
      COUT(2+NDOM+IDOM) = NCOEFS
      NCOEF2 = KSTEP * NCOEFS
      IF (IER .LT. 0) GO TO 900
      N = 0
      DO 60 I = IDOM+1,NDOM
          IKORD = C(I+2)
          IWORK(I+NDOM) = C(I+NDOM+2)
          N = N + IKORD + IWORK(I+NDOM)
          NMIDOM = NMIDOM * IWORK(I+NDOM)
          NCOEF2 = NCOEF2 * IWORK(I+NDOM)
  60  CONTINUE
      CALL DCOPY(N,C(NSTART),1,COUT(NCOUT),1)
      NCOUT = NCOUT + N
      NSTART = NSTART + N
C
C      COMPUTE NEW COEFICIENTS
C
      INDX = 1
      IF (INDX .EQ. IDOM) INDX = 2
      IWORK(IDOM) = -1
      ICOUNT = 0
      DO 300 I = 1,NCOEFS
  310     ICOUNT = ICOUNT + 1
          DENOMN = C(NSKNTS+ICOUNT+KORD) - C(NSKNTS+ICOUNT)
          IF (DENOMN .LE. 0.0D0) GO TO 310
          IWORK(IDOM) = IWORK(IDOM) + 1
          DENOMN = COUT(2+IDOM)/DENOMN
          DO 100 J = 1,NDOM
              IF (J .NE. IDOM) IWORK(J) = 0
 100      CONTINUE
          IWORK(INDX) = -1
          DO 301 J = 1,NMIDOM
C
C       COMPUTE INDICES
C
              IF (NDOM .GT. 1) THEN
                   K = INDX
 120               IWORK(K) = IWORK(K) + 1
                   IF (IWORK(K) .GE. IWORK(K+NDOM)) THEN
                       IWORK(K) = 0
                       K = K + 1
                       IF (K .EQ. IDOM) K = K + 1
                       GO TO 120
                   END IF
               END IF
               K1 = IWORK(NDOM)
               IF (NDOM .EQ. IDOM) K1 = 0
               DO 140 K = NDOM-1,IDOM+1,-1
                   K1 = K1*IWORK(K+NDOM) + IWORK(K)
 140           CONTINUE
               K2 = K1*NCOEFS + IWORK(IDOM)
               K1 = K1*KCOEFS + ICOUNT - 1
               DO 160 K = IDOM-1,1,-1
                   K1 = K1*IWORK(K+NDOM) + IWORK(K)
                   K2 = K2*IWORK(K+NDOM) + IWORK(K)
 160           CONTINUE
               K1 = K1 + NSTART
               K2 = K2 + NCOUT
               DO 301 L2 = 1,NDEP
                   COUT(K2) = (C(K1+KSTEP) - C(K1))*DENOMN
                   K2 = K2 + NCOEF2
                   K1 = K1 + NCOEF1
 301       CONTINUE
 300   CONTINUE
 900   IF (IER .LT. 0) THEN
          CALL DTERR(MODE,SUBNAM,IER,NEED)
          COUT(1) = -1.0D0
      END IF
      RETURN
      END
