      SUBROUTINE DTHSI(N,X,Y,MLT,ICC,HOLD,NHOLD,MAXC,C,NC,IER)
C
C=======================================================================
C
C  PURPOSE:
C          DTHSI computes the Hermite spline interpolant to a
C          set of real data.  It is assumed that the data at each
C          point consists of function values and the first MLT-1
C          derivatives.
C
C  METHOD:
C          Given a set of independent variable values X(I),I=1,...,N
C          and the corresponding dependent variables
C          (Y(1,I),...,Y(MLT,I)),I=1,...,N   DTHSI sets the input
C          parameters for DTCGEN that define a Hermite spline
C          interpolant to the data.
C
C  USAGE:
C          DOUBLE PRECISION X(N),Y(MLT,N),C(MAXC),HOLD(NHOLD)
C          CALL DTHSI(N,X,Y,MLT,ICC,HOLD,NHOLD,MAXC,C,NC,IER)
C
C  INPUT:
C          N       Number of data points, N .GE. MAX( 2, (NDEG+1)/2 ).
C
C          X       Array of N values of the independent variable in
C                  ascending order.
C
C          Y       2-dimensional array containing the values of the
C                  dependent variables.  Y(I,J) is the (I-1)st 
C                  derivative, at the J-th data point.
C
C          MLT     Number of known function and derivative values at
C                  all data points.
C
C          ICC     Initial/continue call code. ICC is used to 
C                  communicate whether this is an initial or continue
C                  call.  See Usage Remarks for further explanation.
C
C                  ICC .LT. 0  Initial call to DTHSI.
C
C                  ICC .GT. 0  Continue call to DTHSI.
C                      .EQ. 1  Overwrite previous coefficients with
C                              new values.
C                      .GT. 1  Append new coefficients to C vector.
C
C           MAXC   Maximum storage allocated to C array.
C
C  STORAGE:
C          HOLD    Hold array of length NHOLD.  This array must be
C                  unchanged by the user between subsequent calls
C                  to DTHSI for ICC .GT. 0.
C
C          NHOLD   The length of array HOLD,
C                  NHOLD .GE. NCOEF*(3*K-2) + 2*K**2 + 2*NCOEF 
C                  where NCOEF = (N-2)*MLT + K and K = NDEG + 1.
C
C  OUTPUT:
C          ICC     Initial/continue call code:
C
C                  ICC .GT. 0  Normal return, ICC has been set to the
C                              absolute value of ICC on input.  DTHSI
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
C                  where NDEP is the number of dependent variables.
C
C          IER     Success/error code.  DTHSI calls the standard error
C                  handler DTERR to print out error and warning messages
C                  for IER .NE. 0.  For IER .LT. 0, DTHSI sets C(1) =
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
C                  IER = -2   N .LT. MAX( 2, (NDEG + 1)/2 ).
C
C                  IER = -3   NHOLD too small, the number of elements
C                             needed is given by the printed error
C                             message.
C
C                  IER = -5   X(I) not in ascending order.
C
C                  IER = -11  Invalid value for ICC. This may be caused
C                             by failing to reset ICC following an
C                             unsuccessful call to DTHSI.
C
C                  IER = -12  A continuation call has been requested
C                             with an invalid spline vector.
C
C                  IER = -14  The user has not supplied enough storage
C                             for the spline vector, C.
C
C                  IER = -20  The coefficient matrix of the interpolation
C                             problem is singular.
C
C                  IER = -30  MLT.GT.KHALF
C                             where KHALF = 1 + NDEG/2.
C
C
C  USAGE REMARKS:
C          The initial/continue call code allows the user to
C          efficiently compute several spline interpolants based on the
C          same independent variable values.  For the initial call,
C          ICC .LT. 0, the user supplies all input arguments.  For
C          subsequent calls, ICC .GT. 0, the user can redefine the
C          Y argument corresponding to the initial X argument.
C
C
C  OTHER DT ROUTINES CALLED:
C
C       DTCGEN
C       DTERR
C
C
C     ******************************************************************
C
      CHARACTER*8 SUBNAM
      DOUBLE PRECISION ZERO, ONE
      PARAMETER  (ZERO=0.D0,ONE=1.D0)
      PARAMETER  (SUBNAM = 'DTHSI   ')
C
C  Arguments:
C
      INTEGER  N, ICC, NHOLD, MAXC, NC, IER
      DOUBLE PRECISION  C(MAXC), X(N), Y(MLT,N), HOLD(NHOLD)
C
C  Internal:
C
      INTEGER  I,K,MODE,NCOEF,NEED,NMIN,NDEG,NRHS,NLHS
      DOUBLE PRECISION  YDUM
C
C=======================================================================
C     Set default values.
C=======================================================================
C
      IER = 0
      NEED = 0
      NDEG = 2*MLT - 1
C
C=======================================================================
C    Check ICC = The continuation flag.
C=======================================================================
C
      IF (ICC.GT.0)
     &  THEN
          GO TO 70
        ELSEIF (ICC.EQ.0) THEN
          IER = -11
          CALL DTERR(MODE,SUBNAM,IER,NEED)
          RETURN
      ENDIF
      C(1) = -ONE
C
C=======================================================================
C    Check N.
C=======================================================================
C
      NMIN = ( NDEG + 1 ) / 2
      IF (N .LT. MAX0(NMIN,2)) THEN
        MODE = 1
        IER = -2
        CALL DTERR(MODE,SUBNAM,IER,NEED)
        RETURN
      ENDIF
C
C=======================================================================
C    Check NHOLD.
C=======================================================================
C
      K = NDEG + 1
      NCOEF = N + K - 2
      NEED =  NCOEF*(3*K - 2) + MAX0(2*K*K,2*NCOEF) + 4*K + 9
      IF (NHOLD .LT. NEED) THEN
        MODE = 2
        IER = -3
        CALL DTERR(MODE,SUBNAM,IER,NEED)
        RETURN
      ENDIF
C
C=======================================================================
C    Check that X(I) are strictly increasing.
C=======================================================================
C
      DO 60 I = 1, N-1
        IF ( X(I) .GE. X(I+1) ) THEN
          ICC  = 0
          MODE = 1
          IER = -5
          CALL DTERR(MODE,SUBNAM,IER,NEED)
          RETURN
        ENDIF
60    CONTINUE
C
C=======================================================================
C     Set remaining inputs and then call DTCGEN.
C=======================================================================
C
   70 CONTINUE
C
C             Set lhs conditions.
C
      NLHS = 0
C
C             Set rhs conditions.
C
      NRHS = 0
C
C             Construct spline array.
C
      CALL DTCGEN(NDEG,X,N,Y(1,1),MLT,YDUM,NLHS,YDUM,NRHS,ICC,
     &            HOLD,NHOLD,MAXC,C,NC,IER)
C
C=======================================================================
C    Check for process errors.
C=======================================================================
C
      IF (IER.NE.0) THEN
        IF (IER.GT.0)
     &    THEN
            MODE = 0
          ELSEIF (IER.EQ.-4) THEN
            MODE = 3
          ELSE
            MODE = 1
        ENDIF
        CALL DTERR(MODE,SUBNAM,IER,NEED)
      ENDIF
C
      RETURN
      END
