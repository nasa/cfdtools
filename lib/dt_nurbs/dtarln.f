        SUBROUTINE DTARLN(C, TLOW, TUP, WORK, NWORK, ARCLEN, IER)
        EXTERNAL DTVAL, DTMCON
        DOUBLE PRECISION DTMCON
        DOUBLE PRECISION C(*), TLOW, TUP, WORK(*), ARCLEN
        DOUBLE PRECISION A, B, DTVAL, DTVAL1, TOL, ERR, ONEMEP
        INTEGER NWORK, IER
        INTEGER N, M, K, NCOEF, ISNG, LIMIT, LEVEL
        INTEGER LC, I
        DATA ISNG, LIMIT, TOL/0, 12, 1.E-4/
C  **
C  **           SUBROUTINE FOR COMPUTING THE ARC-LENGTH OF A 
C  **           (RATIONAL) SPLINE CURVE
C  **
C  **           INPUT:  C       A (RATIONAL) SPLINE ARRAY
C  **                   TLOW    PARAMETER VALUE OF LOWER END POINT
C  **                   TUP     PARAMETER VALUE OF UPPER END POINT
C  **
C  **
C  **           WORKING
C  **           STORAGE WORK    ARRAY OF LENGTH NWORK
C  **                   NWORK   LENGTH OF ARRAY WORK, MUST BE 
C  **                           AT LEAST C(3)**2 + 3 * C(3)
C  **
C  **
C  **           OUTPUT: ARCLEN  ARC LENGTH - A DOUBLE PRECISION NUMBER
C  **                   IER     SUCCESS/FAILURE
C  **                           0       SUCCESS
C  **                          -51      FAILURE: SPLINE IS 
C  **                                           NOT A CURVE
C  **
C  **   ONEMEP IS ONE MINUS EPSILON, THE LARGEST NUMBER STRICTLY LESS THAN ONE
        ONEMEP = 1.0D0 - DTMCON(6)
C  **
C  **           CHECK TO SEE IF C IS A CURVE (N = 1)
C  **
        N = C(1)
        IF (N .EQ. 1) GO TO 10
        IER = -51
        GO TO 900
10      CONTINUE
C  **
C  **           SET ORDER (K), NO. OF COMPONENTS (DABS(M))
C  **           NO. OF SPLINE COEFFICIENTS (NCOEF) AND LENGTH
C  **           OF THE SPLINE ARRAY (LC).
C  **
        K = C(3)
        M = DABS(C(2))
        NCOEF = C(4)
        LC = 5 + NCOEF + K + M * NCOEF
C  **
C  **           BREAK UP THE SPLINE INTO CONTINUOUS SUBPIECES AND FIND
C  **           THE ARC LENGTH OF THE SUBPIECES.
C  **
        ARCLEN = 0.0D0
        A = MAX(TLOW, C(5+K))
        DO 30 I = 6 + K, 5 + NCOEF
          IF (C(I+K-2) .GT. C(I))GO TO 30
          B = ONEMEP*C(I)
          IF (B .LE. A)GO TO 30
          IF (B .GE. TUP)B = TUP
          CALL DTQUAD(DTVAL, C, A, B, TOL, ISNG, LIMIT, DTVAL1
     *          , LEVEL, WORK, NWORK, ERR, IER)
          ARCLEN = ARCLEN + DTVAL1
          IF (B .GE. TUP)GO TO 35
          A = C(I)
30      CONTINUE
        B = MIN(C(6 + NCOEF), TUP)
        IF (B .LE. A) GO TO 35
        CALL DTQUAD(DTVAL, C, A, B, TOL, ISNG, LIMIT, DTVAL1
     *          , LEVEL, WORK, NWORK, ERR, IER)
        ARCLEN = ARCLEN + DTVAL1
35      CONTINUE
900     RETURN
        END
