C<*>

      SUBROUTINE DTNPVL ( X, INCX, C, WORK, NWORK, V, IER )
C
C     MODIFIED BY DEBORAH PARSONS, JUNE 14, 1989
C         TO INCLUDE RATIONAL B-SPLINES
C
C     MODIFIED BY J. MANKE, 9/3/91, FOR RATIONAL B-SPLINES CERTIFICATION
C
C    ===================================================================
C
C    --------------
C    ... PARAMETERS
C    --------------
C
      INTEGER      INCX,     NWORK,    IER
C
      DOUBLE PRECISION       X(*),     C(*),     WORK(*),  V(*)
C
C    ---------------------
C    ...INTERNAL VARIABLES
C    ---------------------
C
      INTEGER      I1,       I2,       I3,       I4,       I5,
     1             I,        ILC,      INCC,     INPT,
     2             IOPT,     ISPAN,    KMAX,     KORD,     MODE,
     3             NBS,      NCOEF,    NDIM,     NDOM,     NEED,
     4             NRNG,     NZERO,    NRAT,     I6
C
      CHARACTER*8  SUBNAM
C
      DOUBLE PRECISION  EPS, DIV
C
      LOGICAL      RATNL
C
C     -------------
C     ... EXTERNALS
C     -------------
C
      DOUBLE PRECISION       DTMCON
C
      EXTERNAL     DTERR,    DTMCON,   DTNPV1,   DTNPV2,   DTSPV1
C
C     -------------------
C     ... SUBROUTINE NAME
C     -------------------
C
      DATA SUBNAM / 'DTNPVL  ' /
C
C     ==================================================================
C
C     ---------------------
C     ...INPUT ERROR CHECKS
C     ---------------------
C
      IER  = 0
      NEED = 0
      MODE = 1
      EPS = DTMCON(6)
C
C     --------------------
C     ... CHECK VALID NRNG
C     --------------------
C
      NRNG = INT ( C(2) )
C
      IF ( NRNG .LT. -1 ) THEN
          RATNL = .TRUE.
          NRAT  = - NRNG
          NRNG  = NRAT - 1
      ELSE IF (NRNG .GT. 0) THEN
          RATNL = .FALSE.
      ELSE
          IER  = -52
          GOTO 9000
      END IF
C
C     --------------------
C     ... CHECK VALID NDOM
C     --------------------
C
      NDOM = INT( C(1) )
C
      IF( NDOM .LT. 1 ) THEN
C
          IER  = -51
          GOTO 9000
C
      END IF
C
C     ---------------------------------------------
C     ... CHECK VALID ORDER AND NUMBER OF B-SPLINES
C     ---------------------------------------------
C
      IOPT  =  2
      INPT  = IOPT + NDOM
      KMAX  = -1
      NCOEF = 1
      NZERO = 1
      ILC   = 3 + 3 * NDOM
C
      DO 10 I = 1, NDOM
C
          KORD  = INT( C(IOPT+I) )
          NBS   = INT( C(INPT+I) )
C
          IF ( KORD .LE. 0 ) THEN
C
              IER  = -1
              GOTO 9000
C
          END IF
C
          IF ( NBS .LT. KORD ) THEN
C
              IER  = -6
              GOTO 9000
C
          END IF
C
          KMAX  = MAX0( KORD, KMAX )
          NCOEF = NCOEF * NBS
          NZERO = NZERO * KORD
          ILC   = ILC + NBS + KORD
C
 10   CONTINUE
C
C     ==================================================================
C     PROCESS B-SPLINE CASE
C     ==================================================================
C
      IF( .NOT.RATNL) THEN
C
C     ---------------------
C     ... CHECK VALID NWORK
C     ---------------------
C
          IF ( NDOM .EQ. 1 ) THEN
C
              NEED = 5 * KORD - 2
C
          ELSE
C
              NZERO = NZERO / INT ( C(3) )
              NEED  = NZERO * ( NRNG + 1 ) + 3 * KMAX + NDOM
              NDIM  = NZERO
C
          END IF
C
          IF( NWORK .LT. NEED ) THEN
C
              IER  = -3
              MODE = 2
              GOTO 9000
C
          END IF
C
C     ------------------------------------
C     ... CALL DTSPV1 IF UNIVARIATE SPLINE
C     ------------------------------------
C
          IF ( NDOM .EQ. 1 ) THEN
C
              ISPAN = INT( C(5) )
              ISPAN = MAX0 ( KORD, ISPAN )
              ISPAN = MIN0 ( NBS,  ISPAN )
              INCC  = 1
C
              CALL DTSPV1 ( X     , KORD  , C(6)  , INCC  , NBS   ,
     1                      C(ILC), INCC  , NBS   , NRNG  , ISPAN , 
     2                      WORK  , V     , IER   )
C
              IF ( IER .NE. 0 ) THEN
                  GOTO 9000
              END IF
C
              C(5) = DBLE( ISPAN )
C
              GOTO 100
C
          END IF
C
C     ------------------------------------------------------------------
C     ... CALL LOWER LEVEL ROUTINE FOR MULTIVARIATE SPLINE EVALUATION
C     ------------------------------------------------------------------
C
          I1 = 1
          I2 = I1 + NZERO
          I3 = I2 + NDOM
          I4 = I3 + NZERO * NRNG
          I5 = I4 + KMAX
C
          IF ( NDOM .EQ. 2 ) THEN
C
              CALL DTNPV1 ( X,        INCX,     NRNG,     C,
     1                      NCOEF,    C(ILC),   NDIM,
     2                      WORK(I1), WORK(I3), WORK(I4),
     3                      WORK(I5), V,        IER                )
C
          ELSE
C
              CALL DTNPV2 ( X,        INCX,     NDOM,     NRNG,    C,
     1                      NCOEF,    C(ILC),   NDIM,     NZERO,
     2                      WORK(I1), WORK(I2), WORK(I3), WORK(I4),
     3                      WORK(I5), V,        IER                 )
C
          END IF
C
          IF ( IER .NE. 0 ) THEN
              GOTO 9000
          END IF
C
 100      CONTINUE
C
C     ==================================================================
C     PROCESS RATIONAL B-SPLINE CASE
C     ==================================================================
C
      ELSE
C
C     ---------------------
C     ... CHECK VALID NWORK
C     ---------------------
C
          IF ( NDOM .EQ. 1 ) THEN
C
              NEED = 5 * KORD - 2 + NRAT
C
          ELSE
C
              NZERO = NZERO / INT ( C(3) )
              NEED  = NZERO * ( NRAT + 1 ) + 3 * KMAX + NDOM + NRAT
              NDIM  = NZERO
C
          END IF
C
          IF( NWORK .LT. NEED ) THEN
C
              IER  = -3
              MODE = 2
              GOTO 9000
C
          END IF
C
C     ------------------------------------
C     ... CALL DTSPV1 IF UNIVARIATE SPLINE
C     ------------------------------------
C
          IF ( NDOM .EQ. 1 ) THEN
C
              ISPAN = INT( C(5) )
              ISPAN = MAX0 ( KORD, ISPAN )
              ISPAN = MIN0 ( NBS,  ISPAN )
              INCC  = 1
C
              I1 = 1
              I6 = I1 + 5 * KORD - 2
C
              CALL DTSPV1 ( X       , KORD    , C(6) , INCC, NBS   ,
     1                      C(ILC)  , INCC    , NBS  , NRAT, ISPAN , 
     2                      WORK(I1), WORK(I6), IER  )
C
              IF ( IER .NE. 0 ) THEN
                  GOTO 9000
              END IF
C
              C(5) = DBLE( ISPAN )
C
              GOTO 200
C
          END IF
C
C     ------------------------------------------------------------------
C     ... CALL LOWER LEVEL ROUTINE FOR MULTIVARIATE SPLINE EVALUATION
C     ------------------------------------------------------------------
C
          I1 = 1
          I2 = I1 + NZERO
          I3 = I2 + NDOM
          I4 = I3 + NZERO * NRAT
          I5 = I4 + KMAX
          I6 = I5 + 2*KMAX
C
          IF ( NDOM .EQ. 2 ) THEN
C
              CALL DTNPV1 ( X,        INCX,     NRAT,     C,
     1                      NCOEF,    C(ILC),   NDIM,
     2                      WORK(I1), WORK(I3), WORK(I4),
     3                      WORK(I5), WORK(I6), IER                    )
C
          ELSE
C
              CALL DTNPV2 ( X,        INCX,     NDOM,     NRAT,     C,
     1                      NCOEF,    C(ILC),   NDIM,     NZERO,
     2                      WORK(I1), WORK(I2), WORK(I3), WORK(I4),
     3                      WORK(I5), WORK(I6), IER                    )
C
          END IF
C
          IF ( IER .NE. 0 ) THEN
              GOTO 9000
          END IF
C
 200      CONTINUE
C
          DIV = WORK(I6+NRAT-1)
          IF (ABS(DIV) .LE. EPS) THEN
              IER = -10
              GOTO 9000
          END IF
          DO 310 I = 1, NRNG
              V(I) = WORK(I6+I-1) / DIV
 310      CONTINUE
C
C     ==================================================================
C
      END IF
C
 9000 CONTINUE
      IF (IER .LT. 0) THEN
          V(1) = DTMCON(1)
          CALL DTERR (MODE, SUBNAM, IER, NEED)
      END IF
C
C     ==================================================================
C
C     ----------
C     ... RETURN
C     ----------
C
      RETURN
      END
