      SUBROUTINE DTSOFF (CIN, DIST, TOLABS, WORK, NWORK, COUT, IER)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC
CC      THE SUBROUTINE DTSOFF IS USED TO DETERMINE OFFSETS TO PLANAR
CC      CURVES.  THE CALLING SEQUENCE HAS THE FOLLOWING FORM:
CC
CC      CALL DTSOFF (CIN, DIST, TOLABS, WORK, NWORK, COUT, IER)
CC
CC      WHERE THE PARAMETERS MEAN
CC
CC      CIN    - THE INPUT SPLINE ARRAY
CC      DIST   - THE EUCLIDEAN DISTANCE AT WHICH TO PLACE THE OFFSET
CC      TOLABS - THE ABSOLUTE TOLERANCE WHICH THE OUTPUT SPLINE MUST
CC               MODEL THE OFFSET CURVE
CC      WORK   - A WORK ARRAY OF LENGTH NWORK
CC      NWORK  - THE LENGTH OF THE WORK ARRAY
CC      COUT   - THE OUTPUT SPLINE ARRAY CONTAINING OFFSET APPROX.
CC      IER    - THE ERROR CONTROL FLAG
CC               IER = 0     NO ERRORS DETECTED
CC               IER = -1    CIN(1) <> 1
CC               IER = -2    CIN(2) <> 2 OR CIN(2) <> -3
CC               IER = -3    INVALID SPLINE ORDER;
CC                           CIN(1) < 1 OR CIN(1) > DTMCON(11)
CC               IER = -4    INVALID NUMBER OF b-SPLINE COEFFICIENTS
CC               IER = -5    INVALID KNOT SEQUENCE
CC               IER = -6    TOLABS <= 0.0
CC               IER = -7    INSUFFICIENT WORKING STORAGE
CC               IER = -100  UNEXPECTED ERROR RETURN FROM LOWER LEVEL
CC
CC      THOMAS GRANDINE
CC      JUNE, 1989
CC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      DOUBLE PRECISION CIN(*), DIST, WORK(*), COUT(*), V(2,2), SCALE
      DOUBLE PRECISION ERROR, MERROR, EST(2), TRUE(2), TOLABS, DTMCON
      INTEGER NWORK, IER, NKNOTS, NDEG, NPTS, NLEFT, ICC, IFAIL, IS
      INTEGER IX, IY, IZ, IW, NORD, NADD, NEED, NC
      CHARACTER*8 SUBNAM
      DATA SUBNAM /'DTSOFF  '/

C----- CHECK THE INPUT FOR VALIDITY

      CALL DTERPT (0)
      CALL DTSCHK (CIN, IER)
      IF (IER .NE. 0) GO TO 9900
      IF (TOLABS .LE. 0.0D0) THEN
        IER = -6
        GO TO 9900
        ENDIF

C----- DETERMINE IF ENOUGH WORKING STORAGE EXISTS FOR FIRST PASS

      NORD = CIN(3)
      NDEG = NORD - 1
      NKNOTS = NORD + CIN(4)
      NC = 2 * NKNOTS + NORD
      NEED = NC * (NORD + 6) + MAX (NC, 3 * NORD - 2) + 7
      IF (NEED .GT. NWORK) THEN
        IER = -7
        GO TO 9900
        ENDIF

C----- DETERMINE THE NEW KNOT SET

      CALL DCOPY (NORD, CIN(6), 1, COUT(6), 1)
      IX = 6 + NDEG
      IY = IX
      IS = NKNOTS - 2 * NORD
      DO 10 IZ = 1,IS
        IX = IX + 1
        IY = IY + 1
        COUT(IY) = CIN(IX)
        IF (CIN(IX+1) .NE. CIN(IX)) THEN
          IY = IY + 1
          COUT(IY) = CIN(IX)
          ENDIF
10      CONTINUE
      IX = IX + 1
      IY = IY + 1
      CALL DCOPY (NORD, CIN(IX), 1, COUT(IY), 1)
      IY = IY + NDEG

C----- DETERMINE THE POINTS TO EVALUATE AT

      NKNOTS = IY - 7
      DO 120 IX = 0,NKNOTS-NDEG
        WORK(2*IX+1) = 0.0D0
        DO 110 IY = 1,NDEG
          WORK(2*IX+1) = WORK(2*IX+1) + COUT(6+IX+IY)
110       CONTINUE
        WORK(2*IX+1) = WORK(2*IX+1) / NDEG
120     CONTINUE
      DO 130 IX = 1,NKNOTS-NDEG
        WORK(2*IX) = 0.5D0 * (WORK(2*IX-1) + WORK(2*IX+1))
130     CONTINUE
      NPTS = 2 * (NKNOTS - NDEG) + 1
      NLEFT = NWORK - 3 * NPTS

C----- NOW EVALUATE POINTS ON THE OFFSET

      DO 210 IX = 1,NPTS
        CALL DTSPDR (WORK(IX), 1, CIN, WORK(3*NPTS+1), NLEFT, V, 2, IER)
        IF (IER .LT. 0) THEN
          IER = -100
          GO TO 9900
          ENDIF
        SCALE = SQRT (V(1,2) ** 2 + V(2,2) ** 2)
        IF (SCALE .NE. 0.0D0) THEN
          WORK(NPTS+IX) = V(1,1) - DIST * V(2,2) / SCALE
          WORK(2*NPTS+IX) = V(2,1) + DIST * V(1,2) / SCALE
          ENDIF
210     CONTINUE

C----- NOW FIT THE OFFSET CURVE

      NKNOTS = NKNOTS - 2 * NDEG
      ICC = -1
      CALL DTPLFT (NPTS, WORK, WORK(NPTS+1), NDEG, ICC, 0, WORK,
     *             NKNOTS, COUT(7+NDEG), WORK(3*NPTS+1), NLEFT, COUT,
     *             IFAIL, IER)
      IF (IER .LT. 0) THEN
        IER = -100
        GO TO 9900
        ENDIF
      ICC = -1
      IS = 5 + COUT(3) + 2 * COUT(4)
      CALL DTPLFT (NPTS, WORK, WORK(2*NPTS+1), NDEG, ICC, 0, WORK,
     *             NKNOTS, COUT(7+NDEG), WORK(3*NPTS+1), NLEFT,
     *             COUT(IS+1), IFAIL, IER)
      IF (IER .LT. 0) THEN
        IER = -100
        GO TO 9900
        ENDIF
      COUT(2) = 2
      IY = COUT(4)
      IX = 5 + COUT(3) + COUT(4)
      CALL DCOPY (IY, COUT(IS+IX+1), 1, COUT(IS+1), 1)

C----- ATTEMPT TO DETERMINE HOW BAD THE SOLUTION IS

      MERROR = 0.0D0
      DO 310 IX = 2,NPTS
        SCALE = 0.5D0 * (WORK(IX-1) + WORK(IX))
        CALL DTSPVL (SCALE, COUT, WORK(3*NPTS+1), NLEFT, EST, IER)
        IF (IER .LT. 0) THEN
          IER = -100
          GO TO 9900
          ENDIF
        CALL DTSPDR (SCALE, 1, CIN, WORK(3*NPTS+1), NLEFT, V, 2, IER)
        IF (IER .LT. 0) THEN
          IER = -100
          GO TO 9900
          ENDIF
        SCALE = SQRT (V(1,2) ** 2 + V(2,2) ** 2)
        IF (SCALE .NE. 0.0D0) THEN
          TRUE(1) = V(1,1) - DIST * V(2,2) / SCALE
          TRUE(2) = V(2,1) + DIST * V(1,2) / SCALE
          ENDIF
        ERROR = SQRT ((EST(1) - TRUE(1)) ** 2 + (EST(2) - TRUE(2)) ** 2)
        MERROR = MAX (MERROR, ERROR)
310     CONTINUE

C----- NOW DETERMINE HOW MANY POINTS TO SPLIT THIS UP INTO

      NADD = (10.0D0 * MERROR / TOLABS) ** (1.0D0 / NORD)
      NKNOTS = NORD + COUT(4)
      CALL DCOPY (NKNOTS, COUT(6), 1, WORK, 1)
      IX = NDEG
      IY = 5 + NDEG
      IS = NKNOTS - 2 * NORD + 1
      DO 420 IZ = 1,IS
        IX = IX + 1
        IY = IY + 1
        COUT(IY) = WORK(IX)
        IF (WORK(IX+1) .NE. WORK(IX)) THEN
          DO 410 IW = 1,NADD
            IY = IY + 1
            COUT(IY) = (IW * WORK(IX+1) + (NADD + 1 - IW) * WORK(IX)) /
     *                 (NADD + 1)
410         CONTINUE
          ENDIF
420     CONTINUE
      IX = IX + 1
      IY = IY + 1
      CALL DCOPY (NORD, WORK(IX), 1, COUT(IY), 1)
      IY = IY + NDEG

C----- DETERMINE IF ENOUGH WORKING STORAGE EXISTS THIS TIME

      NKNOTS = IY - 7
      NC = 2 * NKNOTS + NORD
      NEED = NC * (NORD + 6) + MAX (NC, 3 * NORD - 2) + 7
      IF (NEED .GT. NWORK) THEN
        IER = -7
        GO TO 9900
        ENDIF

C----- DETERMINE THE ADDITIONAL POINTS TO EVALUATE AT

      DO 520 IX = 0,NKNOTS-NDEG
        WORK(2*IX+1) = 0.0D0
        DO 510 IY = 1,NDEG
          WORK(2*IX+1) = WORK(2*IX+1) + COUT(6+IX+IY)
510       CONTINUE
        WORK(2*IX+1) = WORK(2*IX+1) / NDEG
520     CONTINUE
      DO 530 IX = 1,NKNOTS-NDEG
        WORK(2*IX) = 0.5D0 * (WORK(2*IX-1) + WORK(2*IX+1))
530     CONTINUE
      NPTS = 2 * (NKNOTS - NDEG) + 1
      NLEFT = NWORK - 3 * NPTS

C----- NOW EVALUATE POINTS ON THE OFFSET

      DO 610 IX = 1,NPTS
        CALL DTSPDR (WORK(IX), 1, CIN, WORK(3*NPTS+1), NLEFT, V, 2, IER)
        IF (IER .LT. 0) THEN
          IER = -100
          GO TO 9900
          ENDIF
        SCALE = SQRT (V(1,2) ** 2 + V(2,2) ** 2)
        IF (SCALE .NE. 0.0D0) THEN
          WORK(NPTS+IX) = V(1,1) - DIST * V(2,2) / SCALE
          WORK(2*NPTS+IX) = V(2,1) + DIST * V(1,2) / SCALE
          ENDIF
610     CONTINUE

C----- NOW FIT THE OFFSET CURVE

      NKNOTS = NKNOTS - 2 * NDEG
      ICC = -1
      CALL DTPLFT (NPTS, WORK, WORK(NPTS+1), NDEG, ICC, 0, WORK,
     *             NKNOTS, COUT(7+NDEG), WORK(3*NPTS+1), NLEFT, COUT,
     *             IFAIL, IER)
      IF (IER .LT. 0) THEN
        IER = -100
        GO TO 9900
        ENDIF
      ICC = -1
      IS = 5 + COUT(3) + 2 * COUT(4)
      CALL DTPLFT (NPTS, WORK, WORK(2*NPTS+1), NDEG, ICC, 0, WORK,
     *             NKNOTS, COUT(7+NDEG), WORK(3*NPTS+1), NLEFT,
     *             COUT(IS+1), IFAIL, IER)
      IF (IER .LT. 0) THEN
        IER = -100
        GO TO 9900
        ENDIF
      COUT(2) = 2
      IY = COUT(4)
      IX = 5 + COUT(3) + COUT(4)
      CALL DCOPY (IY, COUT(IS+IX+1), 1, COUT(IS+1), 1)

C----- NOW PERFORM THE ERROR PROCESSING STEP

9900  CONTINUE
      CALL DTERPT (1)
      IF (IER .LT. 0) THEN
        COUT(1) = DTMCON (1)
        IF (IER .EQ. -7) THEN
          CALL DTERR (2, SUBNAM, IER, NEED)
        ELSE
          CALL DTERR (1, SUBNAM, IER, 0)
          ENDIF
        ENDIF
      END
