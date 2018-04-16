      SUBROUTINE DTSPJN (C1, C2, WORK, NWORK, COUT, IER)
C
C     SUBROUTINE TO JOIN TWO SPLINES TOGETHER
C
C
C     USAGE
C
C         DOUBLE PRECISION C1(MC1), C2(MC1), WORK(NWORK), COUT(MC3)
C         CALL DTSPJN (C1, C2, WORK, NWORK, COUT, IER)
C
C         MC3 <= 2*(MAX(C1(3),C2(3)) + 1)*MAX(MC1,MC2)
C
C
C     INPUT
C
C         C1, C2  THE TWO SPLINES TO BE JOINED TOGETHER.
C
C
C     WORKING STORAGE
C
C         WORK    ARRAY OF LENGTH NWORK
C
C         NWORK   LENGTH OF ARRAY WORK
C                 NWORK >= MAX(C1(3),C2(3)**2 + 4*MAX(C1(3),C2(3)) + MC3
C
C                          WHERE MC1, MC2 AND MC3 ARE THE LENGTHS OF
C                          C1, C2 AND COUT, RESPECTIVELY.
C
C
C     OUTPUT
C
C         COUT    THE RESULTING JOINED SPLINE ARRAY.
C
C         IER     SUCCESS/ERROR CODE
C                 FOR IER < 0, COUT(1) IS SET TO THE CLOBBER CONSTANT
C                 DTMCON(1).
C
C                 IER = 0     NO ERRORS DETECTED.
C
C                 IER = -1    C1(1) .NE. 1
C
C                 IER = -2    INVALID NUMBER OF DIMENSIONS (DEPENDENT
C                             VARIABLES) IN C1:  C1(2) = 0 OR
C                             C1(2) > DTMCON(11) OR C1(2) = -1 (IN-
C                             DICATING A RATIONAL SPLINE WITH 0 DEPEN-
C                             DENT VARIABLES) OR C1(2) < -DTMCON(11).
C
C                 IER = -3    INVALID SPLINE ORDER IN C1: C1(I) < 1
C                             OR C1(I) > DTMCON(11), FOR SOME
C                             I = 3, ..., C1(1).
C
C                 IER = -4    INVALID NUMBER OF B-SPLINE COEFFICIENTS
C                             FOR EACH OF C1(1) INDEPENDENT VARIABLES
C                             IN C1: THE NUMBER OF COEFFICIENTS IS
C                             LESS THAN THE ORDER OF THE SPLINES, OR
C                             GREATER THAN DTMCON(11).
C
C                 IER = -5    INVALID KNOT SEQUENCE IN C1: FIRST AND
C                             LAST KNOTS DO NOT HAVE MULTIPLICITY EQUAL
C                             TO THE ORDER OF THE SPLINE, OR INTERIOR
C                             KNOTS ARE NOT IN ASCENDING ORDER.
C
C                 IER = -6    C2(1) .NE. 1
C
C                 IER = -7    INVALID NUMBER OF DIMENSIONS (DEPENDENT
C                             VARIABLES) IN C2:  C2(2) = 0 OR
C                             C2(2) > DTMCON(11) OR C2(2) = -1 (IN-
C                             DICATING A RATIONAL SPLINE WITH 0 DEPEN-
C                             DENT VARIABLES) OR C2(2) < -DTMCON(11).
C
C                 IER = -8    INVALID SPLINE ORDER IN C2: C2(I) < 1
C                             OR C2(I) > DTMCON(11), FOR SOME
C                             I = 3, ..., C2(1).
C
C                 IER = -9    INVALID NUMBER OF B-SPLINE COEFFICIENTS
C                             FOR EACH OF C2(1) INDEPENDENT VARIABLES
C                             IN C2: THE NUMBER OF COEFFICIENTS IS
C                             LESS THAN THE ORDER OF THE SPLINES, OR
C                             GREATER THAN DTMCON(11).
C
C                 IER = -10   INVALID KNOT SEQUENCE IN C2: FIRST AND
C                             LAST KNOTS DO NOT HAVE MULTIPLICITY EQUAL
C                             TO THE ORDER OF THE SPLINE, OR INTERIOR
C                             KNOTS ARE NOT IN ASCENDING ORDER.
C
C                 IER = -11   INSUFFIENT WORKING STORAGE (NWORK).
C                             THE RECOMMENDED SIZE OF NWORK IS GIVEN
C                             BY THE PRINTED MESSAGE.
C
C                 IER = -12   ABS(C1(2)) OR ABS(C2(2)) > 10
C
C                 IER = -13   THE NUMBER OF DEPENDENT VARIABLES DO NOT
C                             MATCH.
C
C                 IER = -100  NUMERICAL ERROR IN UPDATING THE DEGREE
C                             DUE TO A PROBLEM WITH THE KNOT SET.
C
C                 IER = -150  ZERO WEIGHT COEFFICIENT IN RATIONAL
C                             SPLINE.
C
C                 IER = -200  UNEXPECTED ERROR IN A LOWER LEVEL ROUTINE.
C
C
C
C     WRITTEN BY TOM GRANDINE
C
C     CONVERTED TO WORK WITH RATIONAL SPLINES BY
C         DEBORAH PARSONS
C         SEPTEMBER 13, 1989
C
C     ALTERED TO COMPUTE THE CORRECT AMOUNT OF STORAGE NEEDED
C         ERIC BRECHNER
C         APRIL, 1990
C
      EXTERNAL DTMCON, DTERPT
C
      INTEGER DEG1, DEG2, NCOUT, NWORK, IER, IX, NDEP, KNOT1, KNOT2  
      INTEGER ADDKNT, ICOEF1, ICOEF2, ICOEF, KNOTS, NEED, ICERR
      INTEGER M, ILEN1, ILEN2, NWORKT, IPC1, IPC2, N, K, I, M1, M2
      INTEGER MC1, MC2
      DOUBLE PRECISION C1(*), C2(*), WORK(*), COUT(*), TRIGHT, TLEFT
      DOUBLE PRECISION DTMCON, VAL1(10,2), VAL2(10,2), ERROR, FIX
      CHARACTER*8 SUBNAM
      LOGICAL CON0, CON1, RATNL1, RATNL2
C
      DATA SUBNAM /'DTSPJN  '/
C
      IER = 0
      CALL DTERPT(0)
C
C     CHECK VALIDITY OF INPUT SPLINES
C
      IF (INT(C1(1)) .NE. 1) THEN
          IER = -1
          GOTO 9900
      ENDIF
      IF (INT(C2(1)) .NE. 1) THEN
          IER = -6
          GOTO 9900
      ENDIF
C
      CALL DTSCHK (C1, IER)
      IF (IER .NE. 0) THEN
          GOTO 9900
      ENDIF
C
      CALL DTSCHK (C2, IER)
      IF (IER .NE. 0) THEN
          IER = IER - 5
          GOTO 9900
      ENDIF
C
C     CHECK FOR NUMBER OF DEPENDENT VARIABLES WITHIN LIMITS
C
      IF (INT(ABS(C1(2))) .GT. 10) THEN
          IER = -12
          GOTO 9900
      ENDIF
C
      IF (INT(ABS(C2(2))) .GT. 10) THEN
          IER = -12
          GOTO 9900
      ENDIF
C
C     CHECK FOR RATIONAL SPLINES
C
      RATNL1 = (C1(2) .LT. 0)
      RATNL2 = (C2(2) .LT. 0)
C
C     CHECK RATIONAL SPLINES FOR VALID WEIGHTS
C
      IF (RATNL1) THEN
          M = INT(ABS(C1(2)))
          K = INT(C1(3))
          N = INT(C1(4))
          ERROR = 1.0
          DO 10 I = 1, N
              ERROR = ERROR * C1(5+K+M*N+I)
   10     CONTINUE
          IF (ERROR .EQ. 0.0) THEN
              IER = -150
              GOTO 9900
          ENDIF
      ENDIF
C
      IF (RATNL2) THEN
          M = INT(ABS(C2(2)))
          K = INT(C2(3))
          N = INT(C2(4))
          ERROR = 1.0
          DO 20 I = 1, N
              ERROR = ERROR * C2(5+K+M*N+I)
   20     CONTINUE
          IF (ERROR .EQ. 0.0) THEN
              IER = -150
              GOTO 9900
          ENDIF
      ENDIF
C
C     FIND THE LENGTH OF THE FIRST SPLINE
C
      M = ABS(INT(C1(2)))
C
C     IF THE SECOND SPLINE IS RATIONAL AND THIS ONE IS NOT,
C     THIS ONE WILL HAVE TO BE LENGTHENED TO SIMULATE A RATIONAL
C
      IF (RATNL2 .AND. .NOT. RATNL1) M = M+1
      DEG1 = INT(C1(3))
      MC1 = 5 + INT(C1(3)) + (M+1)*INT(C1(4))
C
C     FIND THE LENGTH OF THE SECOND SPLINE
C
      M = ABS(INT(C2(2)))
C
C     IF THE FIRST SPLINE IS RATIONAL AND THIS ONE IS NOT,
C     THIS ONE WILL HAVE TO BE LENGTHENED TO SIMULATE A RATIONAL
C
      IF (RATNL1 .AND. .NOT. RATNL2) M = M+1
      DEG2 = INT(C2(3))
      MC2 = 5 + INT(C2(3)) + (M+1)*INT(C2(4))
C
C     ALLOW FOR PROCESSING EXPANSION
C
      IF (DEG1 .LT. DEG2) THEN
          ILEN1 = (DEG2 - DEG1 + 1)*MC1
          ILEN2 = MC2
      ELSE
          ILEN1 = MC1
          ILEN2 = (DEG1 - DEG2 + 1)*MC2
      ENDIF
C
C     CHECK NWORK
C
      NWORKT = MAX(C1(3),C2(3))**2 + 4*MAX(C1(3),C2(3))
      IPC1 = NWORKT + 1
      IPC2 = IPC1 + ILEN1
      NEED = IPC2 + ILEN2 - 1
      IF (NEED .GT. NWORK) THEN
          IER = -11
          GOTO 9900
      ENDIF
C
C     CHECK THAT THE NUMBER OF DEPENDENT VARIABLES MATCHES
C
      M1 = INT(C1(2))
      IF (M1 .LT. 0) M1 = -1-M1
      M2 = INT(C2(2))
      IF (M2 .LT. 0) M2 = -1-M2
      IF (M1 .NE. M2) THEN
          IER = -13
          GOTO 9900
      ENDIF
C
C     COPY THE SPLINES INTO THE WORK ARRAY
C
C
C     IF ONE OF THE SPLINES WAS RATIONAL, THE OTHER WILL HAVE TO BE
C     CONVERTED TO RATIONAL, AS WELL
C
      IF (RATNL2 .AND. .NOT. RATNL1) THEN
          N = INT(C1(4))
          CALL DCOPY (MC1-N, C1(1), 1, WORK(IPC1), 1)
C
C         LENGTHEN THE 1ST SPLINE, SETTING THE WEIGHT COEFFICIENTS TO 1
C
          WORK(IPC1+1) = WORK(IPC1+1) + 1
          CALL DCOPY (N, 1.0D0, 0, WORK(IPC1+MC1-N), 1)
      ELSE
          CALL DCOPY (MC1, C1(1), 1, WORK(IPC1), 1)
          IF (RATNL1) WORK(IPC1+1) = ABS(WORK(IPC1+1))
      ENDIF
C
      IF (RATNL1 .AND. .NOT. RATNL2) THEN
          N = INT(C2(4))
          CALL DCOPY (MC2-N, C2(1), 1, WORK(IPC2), 1)
C
C         LENGTHEN THE 2ND SPLINE, SETTING THE WEIGHT COEFFICIENTS TO 1
C
          WORK(IPC2+1) = WORK(IPC2+1) + 1
          CALL DCOPY (N, 1.0D0, 0, WORK(IPC2+MC2-N), 1)
      ELSE
          CALL DCOPY (MC2, C2(1), 1, WORK(IPC2), 1)
          IF (RATNL2) WORK(IPC2+1) = ABS(WORK(IPC2+1))
      ENDIF
C
C----- MAKE BOTH SPLINES HAVE THE SAME DEGREE
C
  100 CONTINUE
C
      DEG1 = WORK(IPC1+2)
      DEG2 = WORK(IPC2+2)
      IF (DEG1 .LT. DEG2) THEN
          NCOUT = ILEN1
          CALL DTUPDG (WORK(IPC1), 1, NCOUT, WORK(1), NWORKT, COUT, IER)
          IF (IER .NE. 0) THEN
              IER = -200
              GOTO 9900
          ENDIF
          CALL DCOPY (NCOUT, COUT, 1, WORK(IPC1), 1)
          GO TO 100
      ENDIF
C
      IF (DEG1 .GT. DEG2) THEN
          NCOUT = ILEN2
          CALL DTUPDG (WORK(IPC2), 1, NCOUT, WORK(1), NWORKT, COUT, IER)
          IF (IER .NE. 0) THEN
              IER = -200
              GOTO 9900
          ENDIF
          CALL DCOPY (NCOUT, COUT, 1, WORK(IPC2), 1)
          GO TO 100
      ENDIF
C
C     IF ONE OF THE SPLINES WAS RATIONAL, MAKE SURE THAT THE
C     WEIGHT COEFFICIENTS MATCH UP
C
      IF (RATNL1 .OR. RATNL2) THEN
          M = INT(WORK(IPC1+1))
          K = INT(WORK(IPC1+2))
          N = INT(WORK(IPC1+3))
          ICOEF1 = WORK(IPC1+4+K+(M+1)*N)
C
          M = INT(WORK(IPC2+1))
          K = INT(WORK(IPC2+2))
          N = INT(WORK(IPC2+3))
          ICOEF2 = WORK(IPC2+5+K+M*N)
C
          IF (ICOEF1 .NE. ICOEF2) THEN
              IF (ICOEF2 .EQ. 0.0) THEN
                  IER = -150
                  GOTO 9900
              ENDIF
              FIX = DBLE(ICOEF1)/DBLE(ICOEF2)
              M = INT(WORK(IPC2+1))
              K = INT(WORK(IPC2+2))
              N = INT(WORK(IPC2+3))
              DO 500 I = (5+K+N), (4+K+(M+1)*N)
                  WORK(IPC2+I) = WORK(IPC2+I)*FIX
  500         CONTINUE
          ENDIF
      ENDIF
C
C----- NOW ESTABLISH THE OUTPUT SPLINE PARAMETERS
C
      NDEP    = WORK(IPC1+1)
      COUT(1) = 1.0D0
      IF (RATNL1 .OR. RATNL2) THEN
          COUT(2) = -NDEP
      ELSE
          COUT(2) = NDEP
      ENDIF
      COUT(3) = DEG1
      COUT(4) = WORK(IPC1+3) + WORK(IPC2+3)
      KNOT1   = WORK(IPC1+2) + WORK(IPC1+3)
      KNOT2   = WORK(IPC2+2) + WORK(IPC2+3)
      CALL DCOPY (KNOT1, WORK(IPC1+5), 1, COUT(6), 1)
C
C----- DETERMINE THE SMOOTHNESS OF THE JOIN
C
      IX = 5 + WORK(IPC1+2) + WORK(IPC1+3)
      TRIGHT = WORK(IPC1+IX-1)
      TLEFT = WORK(IPC2+5)
      CALL DTSPDR (TRIGHT, 1, WORK(IPC1), WORK, NWORK, VAL1, 10, IER)
      IF (IER .NE. 0) THEN
          IER = -200
          GOTO 9900
      ENDIF
      CALL DTSPDR (TLEFT,  1, WORK(IPC2), WORK, NWORK, VAL2, 10, IER)
      IF (IER .NE. 0) THEN
          IER = -200
          GOTO 9900
      ENDIF
      ERROR = 0.0D0
      DO 610 IX = 1, NDEP
          ERROR = ERROR + (VAL1(IX,1) - VAL2(IX,1)) ** 2
  610 CONTINUE
      ERROR = SQRT (ERROR)
      CON0 = (ERROR .LT. SQRT (DTMCON (6)))
      ERROR = 0.0D0
      DO 620 IX = 1, NDEP
          ERROR = ERROR + VAL1(IX,2) ** 2
  620 CONTINUE
      ERROR = SQRT (ERROR)
      IF (ERROR .NE. 0.0) THEN
          DO 630 IX = 1, NDEP
              VAL1(IX,1) = VAL1(IX,2) / ERROR
  630     CONTINUE
      ENDIF
      ERROR = 0.0D0
      DO 640 IX = 1, NDEP
          ERROR = ERROR + VAL2(IX,2) ** 2
  640 CONTINUE
      ERROR = SQRT (ERROR)
      IF (ERROR .NE. 0.0) THEN
          DO 650 IX = 1, NDEP
              VAL2(IX,1) = VAL2(IX,2) / ERROR
  650     CONTINUE
      ENDIF
      ERROR = 0.0D0
      DO 660 IX = 1, NDEP
          ERROR = ERROR + (VAL1(IX,1) - VAL2(IX,1)) ** 2
  660 CONTINUE
      ERROR = SQRT (ERROR)
      CON1 = ((ERROR .LT. SQRT (DTMCON (6))) .AND. CON0)
C
C----- NOW DETERMINE THE KNOT SCALINGS FOR THE SECOND SPLINE
C
      ERROR = 1.0D0
      IF (CON1) THEN
          TLEFT = WORK(IPC2+5)
          TRIGHT = WORK(IPC2+4+KNOT2)
          ADDKNT = WORK(IPC2+2) + 1
          ERROR = 0.0D0
          ICERR = 0
          DO 710 IX = 1, NDEP
              IF (VAL1(IX,1) * VAL2(IX,2) .NE. 0.0D0) THEN
                  ERROR = ERROR + VAL1(IX,2) * VAL2(IX,1) /
     *                            (VAL1(IX,1) * VAL2(IX,2))
              ELSE
                  ICERR = ICERR + 1
              ENDIF
  710     CONTINUE
          ERROR = ERROR / ( NDEP - ICERR)
          DO 720 IX = ADDKNT,KNOT2
              COUT(4+KNOT1+IX-ADDKNT) = WORK(IPC1+4+KNOT1) +
     *                      (WORK(IPC2+4+IX) - TLEFT) / ERROR
  720     CONTINUE
          COUT(4) = WORK(IPC1+3) + WORK(IPC2+3) - 2
          ICOEF = COUT(4)
          KNOTS = COUT(3) + COUT(4)
          ICOEF1 = WORK(IPC1+3)
          ICOEF2 = WORK(IPC2+3)
          DO 730 IX = 1, NDEP
              CALL DCOPY (ICOEF1, WORK(IPC1+5+KNOT1+(IX-1)*ICOEF1), 1,
     *                    COUT(6+KNOTS+(IX-1)*ICOEF), 1)
              CALL DCOPY (ICOEF2-1, WORK(IPC2+6+KNOT2+(IX-1)*ICOEF2), 1,
     *                    COUT(5+KNOTS+(IX-1)*ICOEF+ICOEF1), 1)
  730     CONTINUE
          GOTO 9900
      ENDIF
C
C----- NOW COPY STUFF OVER IF THE SPLINE IS JUST CONTINUOUS
C
      IF (CON0) THEN
          TLEFT = WORK(IPC2+5)
          TRIGHT = WORK(IPC2+4+KNOT2)
          ADDKNT = WORK(IPC2+2) + 1
          DO 810 IX = ADDKNT,KNOT2
              COUT(5+KNOT1+IX-ADDKNT) = WORK(IPC1+4+KNOT1) +
     *                                 (WORK(IPC2+4+IX) - TLEFT)
  810     CONTINUE
          COUT(4) = WORK(IPC1+3) + WORK(IPC2+3) - 1
          ICOEF = COUT(4)
          KNOTS = COUT(3) + COUT(4)
          ICOEF1 = WORK(IPC1+3)
          ICOEF2 = WORK(IPC2+3)
          DO 820 IX = 1, NDEP
              CALL DCOPY (ICOEF1, WORK(IPC1+5+KNOT1+(IX-1)*ICOEF1), 1,
     *                COUT(6+KNOTS+(IX-1)*ICOEF), 1)
              CALL DCOPY (ICOEF2-1, WORK(IPC2+6+KNOT2+(IX-1)*ICOEF2), 1,
     *                COUT(6+KNOTS+(IX-1)*ICOEF+ICOEF1), 1)
  820     CONTINUE
          GOTO 9900
      ENDIF
C
C----- NOW JUST COPY OVER THE STUFF
C
      TLEFT = WORK(IPC2+5)
      TRIGHT = WORK(IPC2+4+KNOT2)
      ADDKNT = WORK(IPC2+2) + 1
      DO 910 IX = ADDKNT,KNOT2
          COUT(6+KNOT1+IX-ADDKNT) = WORK(IPC1+4+KNOT1) +
     *                             (WORK(IPC2+4+IX) - TLEFT)
  910 CONTINUE
      COUT(4) = WORK(IPC1+3) + WORK(IPC2+3)
      ICOEF = COUT(4)
      KNOTS = COUT(3) + COUT(4)
      ICOEF1 = WORK(IPC1+3)
      ICOEF2 = WORK(IPC2+3)
      DO 920 IX = 1, NDEP
          CALL DCOPY (ICOEF1, WORK(IPC1+5+KNOT1+(IX-1)*ICOEF1), 1,
     *                COUT(6+KNOTS+(IX-1)*ICOEF), 1)
          CALL DCOPY (ICOEF2, WORK(IPC2+5+KNOT2+(IX-1)*ICOEF2), 1,
     *                COUT(6+KNOTS+(IX-1)*ICOEF+ICOEF1), 1)
  920 CONTINUE
C
C     ERROR HANDLING
C
 9900 CONTINUE
      CALL DTERPT(1)
      IF (IER .LT. 0) THEN
          COUT(1) = DTMCON(1)
          IF (IER .EQ. -11) THEN
              CALL DTERR (2, SUBNAM, IER, NEED)
          ELSE
              CALL DTERR (1, SUBNAM, IER, 0)
          ENDIF
      ENDIF
C
      RETURN
      END
