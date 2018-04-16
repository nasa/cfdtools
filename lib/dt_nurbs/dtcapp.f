      SUBROUTINE DTCAPP (CCOEF, PT1, PT2, NORMAL, TOL, WORK, NWORK,
     *                   NEED, CVEC, IER)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC
CC      DTCAPP APPROXIMATES CONIC SECTIONS BY CUBIC SPLINES TO A
CC      USER-SPECIFIED TOLERANCE.  A TRIGONOMETRIC PARAMETRIZATION
CC      IS USED FOR THE CUBIC SPLINE (WHICH WILL, IN GENERAL, BE
CC      DIFFERENT FROM THE RATIONAL PARAMETRIZATION).  THE CALLING
CC      SEQUENCE FOR THE SUBROUTINE IS
CC
CC      CALL DTCAPP (CCOEF, PT1, PT2, NORMAL, TOL, WORK, NWORK, NEED,
CC     *             CVEC, IER)
CC
CC      WHERE THE PARAMETERS HAVE THE FOLLOWING SIGNIFICANCE:
CC
CC      CCOEF  - THE CONIC COEFFICIENTS, I.E. THE VALUES A, B, C,
CC               D, E, AND F SUCH THAT
CC                 2           2
CC               AX  + BXY + CY  + DX + EY + F = 0.
CC      PT1    - THE INITIAL POINT ON THE CONIC CURVE
CC      PT2    - THE FINAL POINT ON THE CONIC CURVE
CC      NORMAL - THE NORMAL VALUE FOR RIGHT-HAND, LEFT-HAND SENSE
CC               NORMAL = +1 TRAVERSE THE CONIC COUNTERCLOCKWISE
CC               NORMAL = -1 TRAVERSE THE CONIC CLOCKWISE
CC      TOL    - THE TOLERANCE TO WHICH THE APPROXIMATION SHOULD BE
CC               ACCURATE.
CC      WORK   - A WORK ARRAY OF LENGTH NWORK.
CC      NWORK  - THE LENGTH OF THE WORK ARRAY.
CC      NEED   - IF IER = -6, THEN NEED CONTAINS WORK ARRAY SPACE NEEDED.
CC      CVEC   - THE SPLINE DEFINITION ARRAY FOR THE RESULTING
CC               VECTOR.
CC      IER    - THE SUCCESS / ERROR FLAG.
CC               IER = 0; SUCCESS, RESULTS COMPUTED.
CC               IER = -1; POINTS DO NOT LIE ON CONIC.
CC               IER = -2; TOL <= 0.
CC               IER = -3; NORMAL <> 1 OR NORMAL <> -1.
CC               IER = -4; COEFFICIENTS DESCRIBE THE EMPTY SET.
CC               IER = -5; POINTS ARE ON DIFFERENT BRANCHES.
CC               IER = -6; NOT ENOUGH WORKING STORAGE.
CC               IER = -100; UNEXPECTED ERROR FROM LOWER LEVEL.
CC
CC      THOMAS GRANDINE
CC      SEPTEMBER, 1989
CC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      DOUBLE PRECISION CCOEF(*), PT1(*), PT2(*), TOL, WORK(*)
      DOUBLE PRECISION CVEC(*), ERR1, ERR2, TYPE, THETA, A, C, D, E, F
      DOUBLE PRECISION SINE, COSINE, RPT1(2), RPT2(2)
      INTEGER NWORK, IER, PREAMB, NORMAL, NC, IX, NEED
      CHARACTER*8 SUBNAM
      DATA SUBNAM /'DTCAPP  '/

C----- CHECK INPUT ERRORS

      IF (TOL .LE. 0.0D0) THEN
        IER = -2
        CALL DTERR (1, SUBNAM, IER, 0)
        CVEC(1) = -1.0D0
        RETURN
        ENDIF
      IF ((NORMAL .NE. -1) .AND. (NORMAL .NE. 1)) THEN
        IER = -3
        CALL DTERR (1, SUBNAM, IER, 0)
        CVEC(1) = -1.0D0
        RETURN
        ENDIF

C----- VERIFY THAT THE TWO POINTS LIE ON THE CONIC

      ERR1 = CCOEF(1) * PT1(1) ** 2 + CCOEF(2) * PT1(1) * PT1(2) +
     *       CCOEF(3) * PT1(2) ** 2 + CCOEF(4) * PT1(1) +
     *       CCOEF(5) * PT1(2) + CCOEF(6)
      ERR2 = CCOEF(1) * PT2(1) ** 2 + CCOEF(2) * PT2(1) * PT2(2) +
     *       CCOEF(3) * PT2(2) ** 2 + CCOEF(4) * PT2(1) +
     *       CCOEF(5) * PT2(2) + CCOEF(6)
      IF ((ABS (ERR1) .GT. TOL) .OR. (ABS (ERR2) .GT. TOL)) THEN
        IER = -1
        CALL DTERR (1, SUBNAM, IER, 0)
        CVEC(1) = -1.0D0
        RETURN
        ENDIF

C----- DETERMINE THE TYPE OF CONIC TO APPROXIMATE

      IF (CCOEF(2) .EQ. 0.0D0) THEN
        THETA = 0.0D0
      ELSE
        THETA = ATAN2 (-CCOEF(2), CCOEF(1) - CCOEF(3)) / 2.0D0
        ENDIF
      SINE = SIN (THETA)
      COSINE = COS (THETA)
      A = CCOEF(1) * COSINE ** 2 - CCOEF(2) * SINE * COSINE +
     *    CCOEF(3) * SINE ** 2
      C = CCOEF(1) * SINE ** 2 + CCOEF(2) * SINE * COSINE +
     *    CCOEF(3) * COSINE ** 2
      D = CCOEF(4) * COSINE - CCOEF(5) * SINE
      E = CCOEF(4) * SINE + CCOEF(5) * COSINE     
      F = CCOEF(6)
      RPT1(1) = COSINE * PT1(1) - SINE * PT1(2)
      RPT1(2) = SINE * PT1(1) + COSINE * PT1(2)
      RPT2(1) = COSINE * PT2(1) - SINE * PT2(2)
      RPT2(2) = SINE * PT2(1) + COSINE * PT2(2)
      TYPE = A * C

C----- DO THE CORRECT THING FOR EACH OF THE THREE CASES

      IER = 0
      CALL DTERPT (0)
      IF (TYPE .GT. TOL) THEN
        CALL DTCAP1 (A, C, D, E, F, RPT1, RPT2, NORMAL, TOL,
     *               WORK, NWORK, NEED, CVEC, IER)
        ENDIF
      IF (ABS (TYPE) .LT. TOL) THEN
        CALL DTCAP2 (A, C, D, E, F, RPT1, RPT2, TOL, CVEC, IER)
        NEED = 0
        ENDIF
      IF (TYPE .LT. -TOL) THEN
        CALL DTCAP3 (A, C, D, E, F, RPT1, RPT2, TOL,
     *               WORK, NWORK, NEED, CVEC, IER)
        ENDIF
      CALL DTERPT (1)
      IF (IER .LT. 0) THEN
        CVEC(1) = -1.0D0
        IF (IER .EQ. -100) THEN
          CALL DTERR (5, SUBNAM, IER, 0)
          RETURN
          ENDIF
        IF (IER .EQ. -6) THEN
          CALL DTERR (2, SUBNAM, IER, NEED)
          RETURN
          ENDIF
        CALL DTERR (1, SUBNAM, IER, 0)
        RETURN
        ENDIF

C----- MAP THE ROTATED SPLINE BACK INTO POSITION

      PREAMB = 5 + CVEC(3) + CVEC(4)
      NC = CVEC(4)
      DO 110 IX = 1,NC
        RPT1(1) = CVEC(PREAMB+IX)
        RPT1(2) = CVEC(PREAMB+NC+IX)
        CVEC(PREAMB+IX) = COSINE * RPT1(1) + SINE * RPT1(2)
        CVEC(PREAMB+NC+IX) = COSINE * RPT1(2) - SINE * RPT1(1)
110     CONTINUE
      END

C----- THE FOLLOWING IS THE SUBROUTINE FOR APPROXIMATING ELLIPSES

      SUBROUTINE DTCAP1 (A, C, D, E, F, PT1, PT2, NORMAL, TOL,
     *                   WORK, NWORK, NEED, CVEC, IER)

      DOUBLE PRECISION A, C, D, E, F, PT1(*), PT2(*), TOL, WORK(*)
      DOUBLE PRECISION CVEC(*), VAL(2), ALFA, BETA, RADIUS, CENTER(2)
      DOUBLE PRECISION NUMER, DENOM, SPARM, EPARM, DTMCON, RPERT, PARM
      INTEGER NWORK, IER, NPCS, PREAMB, NCOEFS, NKNOTS, NCOL, NORMAL
      INTEGER NPTS, IL, IX, IAMAT, IRHS1, IRHS2, IRSD, ILEFT, NLEFT
      INTEGER NEED

C----- DETERMINE THE BASIC CIRCLE PARAMETERS

      ALFA = D ** 2 / (4.0D0 * A ** 2) +
     *       E ** 2 / (4.0D0 * A * C) - F / A
      IF (ALFA .LT. 0.0D0) THEN
        IER = -4
        RETURN
      ELSE
        ALFA = SQRT (ALFA)
        ENDIF
      BETA = D ** 2 / (4.0D0 * A * C) +
     *       E ** 2 / (4.0D0 * C ** 2) - F / C
      IF (BETA .LT. 0.0D0) THEN
        IER = -4
        RETURN
      ELSE
        BETA = SQRT (BETA)
        ENDIF
      RADIUS = MAX (ALFA, BETA)
      CENTER(1) = 0.5D0 * D / A
      CENTER(2) = 0.5D0 * E / C

C----- DETERMINE START AND END PARAMETERS

      NUMER = ALFA * NORMAL * (PT1(2) + CENTER(2))
      DENOM = BETA * (PT1(1) + CENTER(1))
      SPARM = ATAN2 (NUMER, DENOM)
      NUMER = ALFA * NORMAL * (PT2(2) + CENTER(2))
      DENOM = BETA * (PT2(1) + CENTER(1))
      EPARM = ATAN2 (NUMER, DENOM)
      IF (EPARM .LT. SPARM) EPARM = EPARM + 2.0D0 * DTMCON (12)
      NPCS = (EPARM - SPARM) * (0.002D0 * RADIUS / TOL) ** 0.25D0 + 1

C----- NOW INITIALIZE THE SPLINE ARRAY

      CVEC(1) = 1.0D0
      CVEC(2) = 2.0D0
      CVEC(3) = 4.0D0
      NCOEFS = 3 + NPCS
      CVEC(4) = NCOEFS
      NKNOTS = NCOEFS + 4
      CVEC(5) = 0.0D0
      CVEC(6) = SPARM
      CVEC(7) = SPARM
      CVEC(8) = SPARM
      CVEC(9) = SPARM
      CVEC(9+NPCS) = EPARM
      CVEC(10+NPCS) = EPARM
      CVEC(11+NPCS) = EPARM
      CVEC(12+NPCS) = EPARM
      DO 10 IX = 1,NPCS-1
        CVEC(9+IX) = SPARM + IX * (EPARM - SPARM) / NPCS
10      CONTINUE
      PREAMB = 12 + NPCS
      CALL DCOPY (2 * NCOEFS, 0.0D0, 0, CVEC(PREAMB+1), 1)
      CVEC(PREAMB+1) = PT1(1)
      CVEC(PREAMB+2) = PT1(1) + (CVEC(7) - CVEC(10)) * ALFA *
     *                 SIN (SPARM) / 3.0D0
      CVEC(PREAMB+NCOEFS-1) = PT2(1) + (CVEC(PREAMB-1) - CVEC(PREAMB-4))
     *                        * ALFA * SIN (EPARM) / 3.0D0
      CVEC(PREAMB+NCOEFS) = PT2(1)
      CVEC(PREAMB+NCOEFS+1) = PT1(2)
      CVEC(PREAMB+NCOEFS+2) = PT1(2) + (CVEC(10) - CVEC(7)) * BETA *
     *                        NORMAL * COS (SPARM) / 3.0D0
      CVEC(PREAMB+2*NCOEFS-1) = PT2(2) - (CVEC(PREAMB-1) -
     *                          CVEC(PREAMB-4)) * BETA * NORMAL *
     *                          COS (EPARM) / 3.0D0
      CVEC(PREAMB+2*NCOEFS) = PT2(2)
      IF (NCOEFS .EQ. 4) RETURN

C----- CHOP UP THE WORKSPACE INTO THE APPROPRIATE NUMBER OF PIECES

      NPTS = 2 * NCOEFS
      NCOL = NCOEFS - 4
      IAMAT = 1
      IRHS1 = IAMAT + NPTS * NCOL
      IRHS2 = IRHS1 + NPTS
      IRSD = IRHS2 + NPTS
      ILEFT = IRSD + NPTS
      NLEFT = NWORK - ILEFT + 1
      NEED = MAX (16, NCOL+2) - NLEFT + NWORK
      IF (NEED .LT. NWORK) THEN
        IER = -6
        RETURN
        ENDIF

C----- NOW DETERMINE THE POINTS ON THE CURVE AND THE FITTING MATRIX

      DO 110 IX = 1,NPTS
        PARM = SPARM + IX * (EPARM - SPARM) / (NPTS + 1)
        CALL DTBSPL (CVEC(6), NKNOTS, PARM, IL, 4, 0, NPTS, WORK(ILEFT),
     *               NLEFT, WORK(IRSD), IER)
        IF (IER .LT. 0) THEN
          IER = -100
          RETURN
          ENDIF
        CALL DCOPY (NCOL, WORK(IRSD+2), 1, WORK(IX), NPTS)
        WORK(IRHS1+IX-1) = ALFA * COS (PARM) - CENTER(1) 
        WORK(IRHS2+IX-1) = BETA * NORMAL * SIN (PARM) - CENTER(2)
        CALL DTSPVL (PARM, CVEC, WORK(ILEFT), NLEFT, VAL, IER)
        IF (IER .LT. 0) THEN
          IER = -100
          RETURN
          ENDIF
        WORK(IRHS1+IX-1) = WORK(IRHS1+IX-1) - VAL(1)
        WORK(IRHS2+IX-1) = WORK(IRHS2+IX-1) - VAL(2)
110     CONTINUE

C----- COMPUTE THE FACTORIZATION

      CALL DTLSLE (WORK(IAMAT), NPTS, NPTS, NCOL, WORK(IRHS1),
     *             WORK(IRSD), WORK(ILEFT), NLEFT, RPERT, IER)
      IF (IER .LT. 0) THEN
        IER = -100
        RETURN
        ENDIF
      CALL DCOPY (NCOL, WORK(IRHS1), 1, CVEC(PREAMB+3), 1)
      CALL DTLSSL (WORK(IAMAT), NPTS, NPTS, NCOL, WORK(IRHS2),
     *             WORK(IRSD), WORK(ILEFT), NLEFT, RPERT, IER)
      IF (IER .LT. 0) THEN
        IER = -100
        RETURN
        ENDIF
      CALL DCOPY (NCOL, WORK(IRHS2), 1, CVEC(PREAMB+NCOEFS+3), 1)
      END

C----- THIS SUBROUTINE GENERATES THE PARABOLIC STUFF

      SUBROUTINE DTCAP2 (A, C, D, E, F, PT1, PT2, TOL, CVEC, IER)

      DOUBLE PRECISION A, C, D, E, F, PT1(*), PT2(*), TOL, CVEC(*)
      INTEGER IER

C----- INITIALIZE THE SPLINE ARRAY

      CVEC(1) = 1.0D0
      CVEC(2) = 2.0D0
      CVEC(3) = 3.0D0
      CVEC(4) = 3.0D0
      CVEC(5) = 0.0D0
      CVEC(6) = 0.0D0
      CVEC(7) = 0.0D0
      CVEC(8) = 0.0D0
      CVEC(9) = 1.0D0
      CVEC(10) = 1.0D0
      CVEC(11) = 1.0D0
      CVEC(12) = PT1(1)
      CVEC(14) = PT2(1)
      CVEC(15) = PT1(2)
      CVEC(17) = PT2(2)

C----- FILL IN THE REST OF THIS ARRAY

      IF (ABS (C) .LT. TOL) THEN
        CVEC(13) = 0.5D0 * (PT1(1) + PT2(1))
        IF (ABS (E) .LT. TOL) THEN
          IF (ABS (PT1(1) - PT2(1)) .GT. TOL) THEN
            IER = -5
            RETURN
          ELSE
            CVEC(16) = 0.5D0 * (PT1(2) + PT2(2))
            ENDIF
        ELSE
          CVEC(16) = - (A * PT1(1) * PT2(1) + 0.5D0 * D *
     *                  (PT1(1) + PT2(1)) + F) / E
          ENDIF
        RETURN
        ENDIF
      CVEC(16) = 0.5D0 * (PT1(2) + PT2(2))
      IF (ABS (D) .LT. TOL) THEN
        IF (ABS (PT1(2) - PT2(2)) .GT. TOL) THEN
          IER = -5
          RETURN
        ELSE
          CVEC(13) = 0.5D0 * (PT1(1) + PT2(1))
          ENDIF
      ELSE
        CVEC(13) = - (C * PT1(2) * PT2(2) + 0.5D0 * E *
     *                (PT1(2) + PT2(2)) + F) / D
        ENDIF
      END

C----- THE FOLLOWING IS THE SUBROUTINE FOR APPROXIMATING HYPERBOLAS

      SUBROUTINE DTCAP3 (A, C, D, E, F, PT1, PT2, TOL,
     *                   WORK, NWORK, NEED, CVEC, IER)

      DOUBLE PRECISION A, C, D, E, F, PT1(*), PT2(*), TOL, WORK(*)
      DOUBLE PRECISION CVEC(*), VAL(2), ALFA, BETA, RADIUS, CENTER(2)
      DOUBLE PRECISION DPROD, SPARM, EPARM, DTMCON, RPERT, PARM
      INTEGER NWORK, IER, NPCS, PREAMB, NCOEFS, NKNOTS, NCOL, NEED
      INTEGER NPTS, IL, IX, IAMAT, IRHS1, IRHS2, IRSD, ILEFT, NLEFT

C----- DETERMINE THE BASIC CIRCLE PARAMETERS

      ALFA = D ** 2 / (4.0D0 * A ** 2) +
     *       E ** 2 / (4.0D0 * A * C) - F / A
      BETA = D ** 2 / (4.0D0 * A * C) +
     *       E ** 2 / (4.0D0 * C ** 2) - F / C
      IF (ALFA .GT. 0.0D0) THEN
        ALFA = SQRT (ALFA)
        BETA = SQRT (-A / C) * ALFA
      ELSE
        BETA = SQRT (BETA)
        ALFA = SQRT (-C / A) * BETA
        ENDIF
      RADIUS = MAX (ALFA, BETA)
      CENTER(1) = 0.5D0 * D / A
      CENTER(2) = 0.5D0 * E / C

C----- IF PAIR OF CROSSING LINES, THEN FILL IN THE RESULT

      CVEC(1) = 1.0D0
      CVEC(2) = 2.0D0
      CVEC(5) = 0.0D0
      IF (ALFA .LT. 10.0D0 * SQRT (DTMCON (5))) THEN
        CVEC(3) = 2.0D0
        CVEC(6) = 0.0D0
        CVEC(7) = 0.0D0
        CVEC(9) = 1.0D0
        DPROD = ABS ((PT1(1) + CENTER(1)) * (PT2(1) + CENTER(1)) +
     *               (PT1(2) + CENTER(2)) * (PT2(2) + CENTER(2)))
        PARM = SQRT ((PT1(1) + CENTER(1)) ** 2 +
     *               (PT1(2) + CENTER(2)) ** 2) *
     *         SQRT ((PT2(1) + CENTER(1)) ** 2 +
     *               (PT2(2) + CENTER(2)) ** 2)
        IF (ABS (DPROD - PARM) .LT. TOL) THEN
          CVEC(4) = 2.0D0
          CVEC(8) = 1.0D0
          CVEC(10) = PT1(1)
          CVEC(11) = PT2(1)
          CVEC(12) = PT1(2)
          CVEC(13) = PT2(2)
        ELSE
          CVEC(4) = 3.0D0
          CVEC(8) = 0.5D0
          CVEC(10) = 1.0D0
          CVEC(11) = PT1(1)
          CVEC(12) = -CENTER(1)
          CVEC(13) = PT2(1)
          CVEC(14) = PT1(2)
          CVEC(15) = -CENTER(2)
          CVEC(16) = PT2(2)
          ENDIF
        RETURN
        ENDIF

C----- DETERMINE START AND END PARAMETERS

      IF (A .GT. 0.0D0) THEN
        PARM = PT1(2) + CENTER(2)
        SPARM = LOG ((PARM + SQRT (PARM ** 2 + BETA ** 2)) / BETA)
        PARM = PT2(2) + CENTER(2)
        EPARM = LOG ((PARM + SQRT (PARM ** 2 + BETA ** 2)) / BETA)
        IF (EPARM .LT. SPARM) THEN
          BETA = -BETA
          SPARM = -SPARM
          EPARM = -EPARM
          ENDIF
      ELSE
        PARM = PT1(1) + CENTER(1)
        SPARM = LOG ((PARM + SQRT (PARM ** 2 + ALFA ** 2)) / ALFA)
        PARM = PT2(1) + CENTER(1)
        EPARM = LOG ((PARM + SQRT (PARM ** 2 + ALFA ** 2)) / ALFA)
        IF (EPARM .LT. SPARM) THEN
          ALFA = -ALFA
          SPARM = -SPARM
          EPARM = -EPARM
          ENDIF
        ENDIF
      NPCS = (EPARM - SPARM) * (0.002D0 * RADIUS / TOL) ** 0.25D0 + 1

C----- WORK OUT THE SIGNS OF THE VARIOUS PIECES

      IF (A .GT. 0.0D0) THEN
        VAL(1) = PT1(1) - ALFA * COSH (SPARM) + CENTER(1)
        VAL(2) = PT1(2) - BETA * SINH (SPARM) + CENTER(2)
        IF (SQRT (VAL(1) ** 2 + VAL(2) ** 2) .GT. TOL) ALFA = -ALFA
        VAL(1) = PT2(1) - ALFA * COSH (EPARM) + CENTER(1)
        VAL(2) = PT2(2) - BETA * SINH (EPARM) + CENTER(2)
      ELSE
        VAL(1) = PT1(1) - ALFA * SINH (SPARM) + CENTER(1)
        VAL(2) = PT1(2) - BETA * COSH (SPARM) + CENTER(2)
        IF (SQRT (VAL(1) ** 2 + VAL(2) ** 2) .GT. TOL) BETA = -BETA
        VAL(1) = PT2(1) - ALFA * SINH (EPARM) + CENTER(1)
        VAL(2) = PT2(2) - BETA * COSH (EPARM) + CENTER(2)
        ENDIF
      IF (SQRT (VAL(1) ** 2 + VAL(2) ** 2) .GT. TOL) THEN
        IER = -5
        RETURN
        ENDIF

C----- NOW INITIALIZE THE SPLINE ARRAY

      CVEC(3) = 4.0D0
      NCOEFS = 3 + NPCS
      CVEC(4) = NCOEFS
      NKNOTS = NCOEFS + 4
      CVEC(6) = SPARM
      CVEC(7) = SPARM
      CVEC(8) = SPARM
      CVEC(9) = SPARM
      CVEC(9+NPCS) = EPARM
      CVEC(10+NPCS) = EPARM
      CVEC(11+NPCS) = EPARM
      CVEC(12+NPCS) = EPARM
      DO 10 IX = 1,NPCS-1
        CVEC(9+IX) = SPARM + IX * (EPARM - SPARM) / NPCS
10      CONTINUE
      PREAMB = 12 + NPCS
      CALL DCOPY (2 * NCOEFS, 0.0D0, 0, CVEC(PREAMB+1), 1)
      CVEC(PREAMB+1) = PT1(1)
      IF (A .GT. 0.0D0) THEN
        CVEC(PREAMB+2) = PT1(1) + (CVEC(10) - CVEC(7)) * ALFA *
     *                   SINH (SPARM) / 3.0D0
        CVEC(PREAMB+NCOEFS-1) = PT2(1) - (CVEC(PREAMB-1) -
     *                   CVEC(PREAMB-4)) * ALFA * SINH (EPARM) / 3.0D0
      ELSE
        CVEC(PREAMB+2) = PT1(1) + (CVEC(10) - CVEC(7)) * ALFA *
     *                   COSH (SPARM) / 3.0D0
        CVEC(PREAMB+NCOEFS-1) = PT2(1) - (CVEC(PREAMB-1) -
     *                   CVEC(PREAMB-4)) * ALFA * COSH (EPARM) / 3.0D0
        ENDIF
      CVEC(PREAMB+NCOEFS) = PT2(1)
      CVEC(PREAMB+NCOEFS+1) = PT1(2)
      IF (A. GT. 0.0D0) THEN
        CVEC(PREAMB+NCOEFS+2) = PT1(2) + (CVEC(10) - CVEC(7)) * BETA *
     *                          COSH (SPARM) / 3.0D0
        CVEC(PREAMB+2*NCOEFS-1) = PT2(2) - (CVEC(PREAMB-1) -
     *                            CVEC(PREAMB-4)) * BETA *
     *                            COSH (EPARM) / 3.0D0
      ELSE
        CVEC(PREAMB+NCOEFS+2) = PT1(2) + (CVEC(10) - CVEC(7)) * BETA *
     *                          SINH (SPARM) / 3.0D0
        CVEC(PREAMB+2*NCOEFS-1) = PT2(2) - (CVEC(PREAMB-1) -
     *                            CVEC(PREAMB-4)) * BETA *
     *                            SINH (EPARM) / 3.0D0
        ENDIF
      CVEC(PREAMB+2*NCOEFS) = PT2(2)
      IF (NCOEFS .EQ. 4) RETURN

C----- CHOP UP THE WORKSPACE INTO THE APPROPRIATE NUMBER OF PIECES

      NPTS = 2 * NCOEFS
      NCOL = NCOEFS - 4
      IAMAT = 1
      IRHS1 = IAMAT + NPTS * NCOL
      IRHS2 = IRHS1 + NPTS
      IRSD = IRHS2 + NPTS
      ILEFT = IRSD + NPTS
      NLEFT = NWORK - ILEFT + 1
      NEED = MAX (16, NCOL+2) - NLEFT + NWORK
      IF (NEED .LT. NWORK) THEN
        IER = -6
        RETURN
        ENDIF

C----- NOW DETERMINE THE POINTS ON THE CURVE AND THE FITTING MATRIX

      DO 110 IX = 1,NPTS
        PARM = SPARM + IX * (EPARM - SPARM) / (NPTS + 1)
        CALL DTBSPL (CVEC(6), NKNOTS, PARM, IL, 4, 0, NPTS, WORK(ILEFT),
     *               NLEFT, WORK(IRSD), IER)
        IF (IER .LT. 0) THEN
          IER = -100
          RETURN
          ENDIF
        CALL DCOPY (NCOL, WORK(IRSD+2), 1, WORK(IX), NPTS)
        IF (A. GT. 0.0D0) THEN
          WORK(IRHS1+IX-1) = ALFA * COSH (PARM) - CENTER(1) 
          WORK(IRHS2+IX-1) = BETA * SINH (PARM) - CENTER(2)
        ELSE
          WORK(IRHS1+IX-1) = ALFA * SINH (PARM) - CENTER(1) 
          WORK(IRHS2+IX-1) = BETA * COSH (PARM) - CENTER(2)
          ENDIF
        CALL DTSPVL (PARM, CVEC, WORK(ILEFT), NLEFT, VAL, IER)
        IF (IER .LT. 0) THEN
          IER = -100
          RETURN
          ENDIF
        WORK(IRHS1+IX-1) = WORK(IRHS1+IX-1) - VAL(1)
        WORK(IRHS2+IX-1) = WORK(IRHS2+IX-1) - VAL(2)
110     CONTINUE

C----- COMPUTE THE FACTORIZATION

      CALL DTLSLE (WORK(IAMAT), NPTS, NPTS, NCOL, WORK(IRHS1),
     *             WORK(IRSD), WORK(ILEFT), NLEFT, RPERT, IER)
      IF (IER .LT. 0) THEN
        IER = -100
        RETURN
        ENDIF
      CALL DCOPY (NCOL, WORK(IRHS1), 1, CVEC(PREAMB+3), 1)
      CALL DTLSSL (WORK(IAMAT), NPTS, NPTS, NCOL, WORK(IRHS2),
     *             WORK(IRSD), WORK(ILEFT), NLEFT, RPERT, IER)
      IF (IER .LT. 0) THEN
        IER = -100
        RETURN
        ENDIF
      CALL DCOPY (NCOL, WORK(IRHS2), 1, CVEC(PREAMB+NCOEFS+3), 1)
      END
