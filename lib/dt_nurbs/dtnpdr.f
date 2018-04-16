C<*>
      SUBROUTINE DTNPDR ( X, INCX, IDER, C, WORK, NWORK, V, IER )
C
C     MODIFIED BY DEBORAH PARSONS, JULY 10, 1989
C         TO INCLUDE RATIONAL SPLINES
C
C     MODIFIED BY J. MANKE, 9/3/91, FOR RATIONAL B-SPLINE CERTIFICATION
C
C    ===================================================================
C
C    --------------
C    ... PARAMETERS
C    --------------
C
      INTEGER      INCX,     IDER(*),  NWORK,    IER
C
      DOUBLE PRECISION       X(*),     C(*),     WORK(*),  V(*)
C
C    ---------------------
C    ...INTERNAL VARIABLES
C    ---------------------
C
      INTEGER      I1,       I2,       I3,       I4,       I5,
     1             I,        J,        ILC,      INCC,     INPT,
     2             IOPT,     ISPAN,    KMAX,     KORD,     MODE,
     3             NBS,      NCOEF,    NDIM,     NDOM,     NEED,
     4             NRNG,     NZERO,    K,        NRAT
C
      INTEGER MAXDR, LNAT, LNBSV, LNWK, LNV, IV, IAT, IBSVAL, IWK, 
     *        M, N, IJ, IQ, IPQ, IPVX, IPVW
C
      CHARACTER*8  SUBNAM
C
      DOUBLE PRECISION  EPS, SUM, DIV, TX(2), CJM, CQN
C
      LOGICAL      RATNL
C
C     -------------
C     ... EXTERNALS
C     -------------
C
      DOUBLE PRECISION       DTMCON, DTFCTL
C
      EXTERNAL     DTERR, DTMCON, DTNPD1, DTNPD2, DTSPD1, DTFCTL
C
C     -------------------
C     ... SUBROUTINE NAME
C     -------------------
C
      DATA SUBNAM / 'DTNPDR  ' /
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
      EPS = DTMCON(5)
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
          C(2)  = NRAT
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
C     ---------------------------------------------------------
C     ... CHECK VALID DERIVATIVE, ORDER AND NUMBER OF B-SPLINES
C     ---------------------------------------------------------
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
          IF ( IDER(I) .LT. 0 ) THEN
C
              IER = -53
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
              NEED = KORD * ( KORD + 3 + IDER(1) )
     1             + NRNG * ( IDER(1) + 1 )
C
          ELSE
C
              NZERO = NZERO / INT ( C(3) )
              NDIM  = NZERO
              NEED  = NZERO * ( NRNG + 1 ) + 2 * KMAX * KMAX + NDOM
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
C     ... CALL DTSPD1 IF UNIVARIATE SPLINE
C     ------------------------------------
C
          IF ( NDOM .EQ. 1 ) THEN
C
              ISPAN = INT( C(5) )
              ISPAN = MAX0 ( KORD, ISPAN )
              ISPAN = MIN0 ( NBS,  ISPAN )
              INCC  = 1
              I1    = 1
              I2    = I1 + NRNG * ( IDER(1) + 1 )
C
              CALL DTSPD1 ( X     , IDER(1), KORD  , C(6)  , INCC  ,
     1                      NBS   , C(ILC) , INCC  , NBS   , NRNG  ,
     2                      ISPAN , WORK(I2),WORK(I1)      , NRNG  ,
     3                      IER   )
C
              IF (IER .NE. 0) THEN
                  GOTO 9000
              ENDIF
C
              DO 110 I = 1, NRNG
                 V(I) = WORK(I+IDER(1)*NRNG)
 110          CONTINUE
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
          IF ( NDOM .EQ. 2 ) THEN
C
              MAXDR = MAX (IDER(1)+1, IDER(2)+1)
              KMAX  = MAX(C(3), C(4))
              LNAT  = KMAX * KMAX * C(2)
              LNBSV = KMAX * MAXDR
              LNWK  = KMAX * KMAX
              LNV   = (IDER(1)+1) * (IDER(2)+1) * C(2)
C
C             RE-CHECK THE WORKING STORAGE
C
              NEED  = LNAT+LNBSV+LNWK+LNV
              IF (NEED .GT. NWORK) THEN
                  IER = -3
                  MODE = 2
                  GOTO 9000
              ENDIF
              IV     = 1
              IAT    = IV+LNV
              IBSVAL = IAT+LNAT
              IWK    = IBSVAL + LNBSV
              TX(1)  = X(1)
              TX(2)  = X(1+INCX)
              M      = IDER(1)
              N      = IDER(2)
              CALL DTNPD1 (TX, M, N, C, WORK(IAT), KMAX, WORK(IBSVAL), 
     *                     WORK(IWK), WORK(IV), IER)
              IF (IER .NE. 0) THEN
                  GOTO 9000
              ENDIF
C
              DO 120 I = 1, NRNG
                  V(I) = WORK((I-1)*(M+1)*(N+1) + N*(M+1) + (M+1))
 120          CONTINUE
C
          ELSE
C
              I1 = 1
              I2 = I1 + NZERO
              I3 = I2 + NDOM
              I4 = I3 + NZERO * NRNG
              I5 = I4 + KMAX * KMAX
C
              CALL DTNPD2 ( X,        INCX,     IDER,     NDOM,
     1                      NRNG,     C,        NCOEF,    C(ILC),
     2                      NDIM,     NZERO,    WORK(I1), WORK(I2),
     3                      WORK(I3), KMAX,     WORK(I4), WORK(I5),
     4                      V,        IER                    )
C
          END IF
C
          IF ( IER .NE. 0 ) THEN
              GOTO 9000
          END IF
C
 200      CONTINUE
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
              NEED = KORD * ( KORD + 3 + IDER(1) )
     1             + NRAT * ( IDER(1) + 1 )
C
          ELSE
C
              NZERO = NZERO / INT ( C(3) )
              NDIM  = NZERO
              NEED  = NZERO * ( NRAT + 1 ) + 2 * KMAX * KMAX + NDOM
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
C     ... CALL DTSPD1 IF UNIVARIATE SPLINE
C     ------------------------------------
C
          IF ( NDOM .EQ. 1 ) THEN
C
              ISPAN = INT( C(5) )
              ISPAN = MAX0 ( KORD, ISPAN )
              ISPAN = MIN0 ( NBS,  ISPAN )
              INCC  = 1
              I1    = 1
              I2    = I1 + NRAT * ( IDER(1) + 1 )
C
              CALL DTSPD1 ( X     , IDER(1), KORD  , C(6)  , INCC  ,
     1                      NBS   , C(ILC) , INCC  , NBS   , NRAT  ,
     2                      ISPAN , WORK(I2),WORK(I1)      , NRAT  ,
     3                      IER   )
C
              IF (IER .NE. 0) THEN
                  GOTO 9000
              ENDIF
C
C             -- PRE-CHECK FOR DIVIDE-BY-ZERO
              DIV = WORK(NRAT)
              IF (ABS(DIV) .LE. EPS) THEN
                  IER = -10
                  GOTO 9000
              ENDIF
C             -- LOOP THROUGH VALUES
              DO 250 I = 1, NRNG
                  WORK(I) = WORK(I)/DIV
C                 -- INITIALIZE WORK
                  DO 210 J = 1, IDER(1)
                      WORK(I2-1+J) = 1.0
 210             C ONTINUE
C                 -- LOOP THROUGH DERIVATIVES
                  DO 240 J = 1, IDER(1)
                      SUM = 0.0
                      DO 220 K = 1, J
                          SUM = SUM + WORK(I2-1+K) *
     *                    WORK((J-K)*NRAT+I) * WORK(K*NRAT+NRAT)
 220                  CONTINUE
                      WORK(J*NRAT+I) = (WORK(J*NRAT+I) - SUM)/DIV
                      DO 230 K = 1, J
                          WORK(I2-1+K) = DBLE(J+1)/DBLE(J+1-K) *
     *                    WORK(I2-1+K)
 230                 CONTINUE
 240              CONTINUE
 250          CONTINUE
C
              DO 260 I = 1, NRNG
                 V(I) = WORK(I+IDER(1)*NRAT)
 260          CONTINUE
C
              C(5) = DBLE( ISPAN )
C
              GOTO 400
C
          END IF
C
C     ------------------------------------------------------------------
C     ... CALL LOWER LEVEL ROUTINE FOR MULTIVARIATE SPLINE EVALUATION
C     ------------------------------------------------------------------
C
          IF ( NDOM .EQ. 2 ) THEN
C
              MAXDR = MAX (IDER(1)+1, IDER(2)+1)
              KMAX  = MAX(C(3), C(4))
              LNAT  = KMAX * KMAX * C(2)
              LNBSV = KMAX * MAXDR
              LNWK  = KMAX * KMAX
              LNV   = (IDER(1)+1) * (IDER(2)+1) * C(2)
C
C             RE-CHECK THE WORKING STORAGE
C
              NEED  = LNAT+LNBSV+LNWK+LNV
              IF (NEED .GT. NWORK) THEN
                  IER = -3
                  MODE = 2
                  GOTO 9000
              ENDIF
              IV     = 1
              IAT    = IV+LNV
              IBSVAL = IAT+LNAT
              IWK    = IBSVAL + LNBSV
              TX(1)  = X(1)
              TX(2)  = X(1+INCX)
              M      = IDER(1)
              N      = IDER(2)
              CALL DTNPD1 (TX, M, N, C, WORK(IAT), KMAX, WORK(IBSVAL), 
     *                     WORK(IWK), WORK(IV), IER)
              IF (IER .NE. 0) THEN
                  GOTO 9000
              ENDIF
C
C             COMPUTE THE DERIVATIVES OF THE RATIONAL SPLINE FROM
C             THE PARTIAL DERIVATIVES CONTAINED IN THE DTNPD1 V ARRAY
C
C             ALGORITHM:  Q(I,J) IS THE PARTIAL WITH RESPECT TO U^I 
C                         AND V^J.
C
C                         CIJ IS J! / (I!(J-I)!)
C
C                         SO,
C
C           X(M,N) - (SUM(J=0,M) (SUM(Q=0,N) CJM*CQN*Q(M-J,N-Q)*W(J,Q)))
C  Q(M,N) = ------------------------------------------------------------
C                                           W
C             (J+Q .NE. 0)
C
              DO 350 K = 1, NRNG
C
                  IPVW = (NRAT-1)*(M+1)*(N+1)+1
                  DIV = WORK(IPVW)
                  IF (ABS(DIV) .LE. EPS) THEN
                      IER = -10
                      GOTO 9000
                  ENDIF
C
C                 COMPUTE WORK(LNV+1) = Q(0,0) = V(1,1,K) / V(1,1,NRAT)
C
                  IPQ  = LNV+1
                  IPVX = (K-1)*(M+1)*(N+1)+1
                  WORK(IPQ) = WORK(IPVX) / DIV
                  DO 340 I = 1, N+M
                      DO 330 J = 0, MIN(I,N)
                          IF ((I-J) .LE. M) THEN
                              SUM = 0.0
                              DO 320 IJ = 0, I-J
                                  DO 310 IQ = 0,J
                                      IF ((IJ+IQ) .NE. 0) THEN
                                          CJM = DTFCTL(I-J)/
     *                                       (DTFCTL(IJ)*DTFCTL(I-J-IJ))
                                          CQN = DTFCTL(J)/
     *                                       (DTFCTL(IQ)*DTFCTL(J-IQ))
                                          IPQ  = LNV+(J-IQ)*(M+1)+
     *                                           I-J-IJ+1
                                          IPVW = (NRAT-1)*(M+1)*(N+1)+
     *                                        IQ*(M+1)+IJ+1
                                          SUM = SUM + CJM * CQN *
     *                                       WORK(IPQ) * WORK(IPVW)
                                      ENDIF
 310                              CONTINUE
 320                          CONTINUE
                              IPQ  = LNV+(J)*(M+1)+I-J+1
                              IPVX = (K-1)*(M+1)*(N+1)+J*(M+1)+(I-J)+1
                              WORK(IPQ) = (WORK(IPVX)-SUM)/ DIV
                          ENDIF
 330                  CONTINUE
 340              CONTINUE
                  IPQ  = LNV+N*(M+1)+M+1
                  V(K) = WORK(IPQ)
 350          CONTINUE
C
          ELSE
C
              IER = -54
              GOTO 9000
C
          END IF
C
          IF ( IER .NE. 0 ) THEN
          GOTO 9000
          END IF
C
 400      CONTINUE
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
C     PUT C(2) BACK IF CHANGED
C
      IF (RATNL) THEN
          C(2) = - NRAT
      ENDIF
C
C     ==================================================================
C
C     ----------
C     ... RETURN
C     ----------
C
      RETURN
      END
