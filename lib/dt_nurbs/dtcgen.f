      SUBROUTINE DTCGEN(NDEG,X,N,Y,MLT,YLHS,NLHS,YRHS,NRHS,ICC,
     &                 HOLD,NHOLD,MAXC,C,NC,IER)
C
C=======================================================================
C
C  PURPOSE:
C          DTCGEN constructs the spline C array representing the
C          set of real data in the Y array.  This is a low level 
C          routine which assumes that argument checking has occured
C          in the calling routine.
C
C  METHOD:
C          Given a set of independent variable values X(I),I=1,...,N
C          and the corresponding dependent variables YLHS,Y,YRHS,
C          DTCGEN forms the band matrix of b-spline values at the 
C          data points.  The linear system of equations choosing the
C          b-spline coefficients to match the data is then factored
C          and solved.
C
C  USAGE:
C          DOUBLE PRECISION X(N),Y(MLT,N),YLHS(NLHS),YRHS(NRHS),
C                           HOLD(NHOLD),C(MAXC)
C          CALL DTCGEN(NDEG,X,N,Y,MLT,YLHS,NLHS,YRHS,NRHS,ICC,
C                     HOLD,NHOLD,MAXC,C,NC,IER)
C
C  INPUT:
C          NDEG    Degree of spline, NDEG .GE. 0.
C
C          X       Array of N values of the independent variable in
C                  ascending order.
C
C          N       Number of data points, N .GE. MAX( 2, (NDEG+1)/2 )
C
C          Y       2-dimensional array containing the values of the
C                  dependent variables at data points. 
C                  Y(I,J) is the (I-1)st derivative, at the J-th
C                  data point.
C
C          MLT     Number of known function and derivative values at
C                  all data points.
C
C          YLHS    Additional derivative values at left endpoint.
C
C          NLHS    Number of elements in array YLHS.
C
C          YRHS    Additional derivative values at right endpoint.
C
C          NRHS    Number of elements in array YRHS.
C
C          ICC     Initial/continue call code. ICC is used to
C                  communicate whether this is an initial or
C                  continue call.  See usage remarks for further
C                  explanation.
C
C                  ICC .LT. 0  Initial call to DTCGEN.
C
C                  ICC .GT. 0  Continue call to DTCGEN.
C                         = 1  Overwrite old coefficients.
C                        >= 2  Append new coefficients.
C
C           MAXC   Maximum storage allocated to C array.
C
C  STORAGE:
C           HOLD   Hold array of length NHOLD.  This array must be
C                  unchanged by the user between subsequent calls
C                  to DTCGEN for ICC .GT. 0.
C
C                                       LENGTH
C
C                1  LU decomp           NCOEF*LDA
C               I2  HOLD vector         2*NCOEF
C               I3  HDSPB1 work         2*KORD*KORD
C
C
C          NHOLD   The length of array HOLD,
C                  NHOLD .GE. 3*K*NCOEF + 2*K**2
C                  where NCOEF = (N-2)*MLT + K  and K = NDEG + 1.
C
C  OUTPUT:
C          ICC     Initial/continue call code:
C
C                  ICC .GT. 0  Normal return, ICC has been set to the
C                              absolute value of ICC on input.  DTCGEN
C                              may be called again for a different set
C                              of values for the dependent variable.
C
C                  ICC .EQ. 0  An error has been detected, IER .NE. 0,
C                              Check value of IER.
C
C          C       Array of elements containing the information needed
C                  to evaluate the spline.  See DTSPVL abstract for a
C                  description of the contents of C.
C
C          NC      Length of spline array.
C                  NC = 5 + (NDEP + 1)*NCOEF + K.
C                  where  NDEP is the number of dependent variables.
C
C          IER     Success/error code.  DTCGEN calls the standard error
C                  handler DTERR to print out error and warning messages
C                  for IER .NE. 0.  For IER .LT. 0, DTCGEN sets C(1) =
C                  -1.0.
C
C                  IER =  0   Success.
C
C                  IER =  1   Results computed but they are sensitive
C                             to inaccuracies in the entries of X and
C                             Y.
C
C                  IER =  2   Results computed but they are strongly
C                             sensitive to inaccuracies in the entries
C                             of X and Y.
C
C                  IER = -1   NDEG .LT. 0.
C
C                  IER = -2   N .LT. MAX( 2, (NDEG+1)/2 ).
C
C                  IER = -3   NHOLD too small, The number of elements
C                             needed is given by the printed error
C                             message.
C
C                  IER = -11  Invalid value for ICC. This may be caused
C                             by failing to reset ICC following an
C                             unsuccessful call to DTCGEN.
C
C                  IER = -12  A continuation call has been requested
C                             with an invalid spline vector.
C
C                  IER = -13  A continuation call has been requested
C                             for input which is inconsistent with the
C                             c-vector.  Check that the number of data
C                             points, multiplicities, and order have
C                             been kept the same.
C
C                  IER = -14  The user has not supplied enough storage
C                             for the spline vector, C.
C
C                  IER = -20  The coefficient matrix for the interpolation
C                             problem is singular.
C
C  USAGE REMARKS:
C          The initial/continue call code allows the user to 
C          efficiently compute several spline interpolants based on
C          the same independent variable values.  For the initial call,
C          ICC .LT. 0, the user supplies all input arguments.  For
C          subsequent calls, ICC .GT. 0, the user can redefine the
C          Y argument corresponding to the initial X argument.
C
C
C  OTHER DT ROUTINES CALLED
C
C       DTCMAT
C       DGBCO
C       DGBSL
C
C     ******************************************************************
C
      DOUBLE PRECISION ZERO, ONE, C0P5, C2P0, C3P0
      PARAMETER  (ZERO=0.D0,ONE=1.D0,C0P5=0.5D0,C2P0=2.D0,C3P0=3.D0)
C
C  Arguments:
C
      INTEGER  ICC, IER, N, NDEG, MLT, NLHS, NRHS, NHOLD, MAXC
      DOUBLE PRECISION  C(MAXC),X(*),Y(MLT,*),HOLD(NHOLD),YLHS(*),
     &                  YRHS(*)
C
C  Internal:
C
      INTEGER  I,I2,I3,IL,IR,J,KEVEN,KHALF,KORD,LDA,ML,NCOEF
      DOUBLE PRECISION  RCOND, DIGCAL, DIGMAX, DTMCON
C
C     ******************************************************************
C
C=======================================================================
C     Set default values.
C=======================================================================
C
      IER = 0
      NM2 = N - 2
      KORD = NDEG + 1
      KHALF = (KORD+1)/2
      KEVEN = 2*KHALF
      NCOEF = KEVEN
      IF (N.GT.2)  NCOEF = NCOEF + MLT*NM2
C
C=======================================================================
C    Check ICC = The continuation flag.
C=======================================================================
C
      IF (ICC.GT.0)
     &  THEN
C
C             Check for error from previous matrix factorization.
C
          IF (C(1) .LT. ZERO) THEN
            IER = -12
            RETURN
          ENDIF
C
C             Check that input is consistent with c-vector.
C
          IF (NCOEF .NE. NINT(C(4))) THEN
            IER = -13
            RETURN
          ENDIF
          GO TO 320
        ELSEIF (ICC.EQ.0) THEN
          IER = -11
          RETURN
      ENDIF
C
C=======================================================================
C    Set the C-vector.
C=======================================================================
C
      C(1) = ONE
      C(2) = ONE
      C(3) = DBLE( FLOAT( KORD ) )
      C(4) = DBLE( FLOAT( NCOEF ) )
      C(5) = C(3)
      IL = 5
C
C=======================================================================
C    Set knots.
C=======================================================================
C
C  Left endpoint:
C
      DO 220 I=1,KORD
        IL = IL + 1
        C(IL) = X(1)
  220 CONTINUE
C
C  Interior:
C
      IF (KORD.EQ.KEVEN)
     &  THEN
C
C             Knots are at data points.
C
          DO 240 I=1,NM2
            DO 235 J=1,MLT
              IL = IL + 1
              C(IL) = X(I+1)
  235       CONTINUE
  240     CONTINUE
C
        ELSE
C
C             Knots are midway between data points.
C
          DO 250 I=1,N-1
            DO 245 J=1,MLT
              IL = IL + 1
              C(IL) = C0P5*(X(I+1)+X(I))
  245       CONTINUE
  250     CONTINUE
      ENDIF
C
C  Right endpoint.
C
      DO 260 I=1,KORD
        IL = IL + 1
        C(IL) = X(N)
  260 CONTINUE
C
C=======================================================================
C    Call DTCMAT to define the interpolation matrix.
C=======================================================================
C
      ML = KORD - 1
      LDA = 3*ML + 1
      I2 = 1 + LDA*NCOEF
      I3 = I2 + KORD*KORD
      CALL DTCMAT(X,N,C(6),NCOEF,KORD,MLT,NLHS,NRHS,LDA,HOLD,
     +            HOLD(I3),HOLD(I2))
C
C=======================================================================
C    Factor the interpolation matrix using LINPACK band factor.
C=======================================================================
C
      I3 = I2 + NCOEF
      RCOND = ZERO
      CALL DGBCO(HOLD,LDA,NCOEF,ML,ML,HOLD(I2),RCOND,HOLD(I3))
C
C  Test for exact singularity.
C
      IF (RCOND.EQ.ZERO) THEN
        IER  = -20
        RETURN
      ENDIF
C
C  Check condition number.
C
      DIGMAX = -DLOG10( DTMCON(5) )
      DIGCAL = -C3P0*DLOG10( RCOND )
      MODE   = 0
      IF ( DIGCAL .GE. DIGMAX ) THEN
        IF ( DIGCAL .LT. C2P0*DIGMAX )
     &    THEN
C
C             Matrix is poorly conditioned.
C
            IER = 1
C
          ELSE
C
C             Matrix is badly conditioned.
C
            IER = 2
        ENDIF
      ENDIF
C
C=======================================================================
C    Check MAXC.
C=======================================================================
C
  320 CONTINUE
      NDPOLD = NINT(C(2))
      IF (ICC.LT.2) 
     &  THEN
          NDEP = NDPOLD
        ELSE
          NDEP = NDPOLD + 1
      ENDIF
      IL = 6 + KORD + NDEP*NCOEF
      NC = IL + NCOEF - 1
C
      IF (NC.GT.MAXC) THEN
        IER = -14
        RETURN
      ENDIF
      C(2) = DBLE( FLOAT(NDEP) )
C
C=======================================================================
C    Load right hand side vector into C.
C=======================================================================
C
C  Left boundary conditions.
C
      IR = IL - 1
      NLPT = MLT + NLHS
      KLHS = KHALF - NLPT
      DO 410 I=1,KLHS
        C(IR+I) = ZERO
  410 CONTINUE
      DO 415 I=KLHS+1,KLHS+MLT
        J = I - KLHS
        C(IR+I) = Y(J,1)
  415 CONTINUE
      DO 420 I=KLHS+MLT+1,KHALF
        J = I - KLHS - MLT
        C(IR+I) = YLHS(J)
  420 CONTINUE
      IR = IR + KHALF
C
C  Interior interpolation points.
C
      DO 460 I=2,N
        DO 450 J=1,MLT
          C(IR+J) = Y(J,I)
  450   CONTINUE
        IR = IR + MLT
  460 CONTINUE
C
C  Right boundary conditions.
C
      DO 480 I=1,NRHS
        C(IR+I) = YRHS(I)
  480 CONTINUE
      DO 485 I=NRHS+1,KHALF-MLT
        C(IR+I) = ZERO
  485 CONTINUE
C
C=======================================================================
C    Solve band matrix for B-spline coefficients.
C=======================================================================
C
      ML = KORD - 1
      LDA = 3*ML + 1
      I2 = 1 + LDA*NCOEF
      CALL DGBSL(HOLD,LDA,NCOEF,ML,ML,HOLD(I2),C(IL),0)
      ICC = IABS( ICC )
C
      RETURN
      END
