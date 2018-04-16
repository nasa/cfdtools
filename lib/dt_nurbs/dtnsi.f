      SUBROUTINE DTNSI(N,X,Y,NDEG,ICC,HOLD,NHOLD,MAXC,C,NC,IER)
C
C=======================================================================
C
C  PURPOSE:
C          DTNSI computes the natural spline interpolant to a set
C          of real data.
C
C  METHOD:
C          Given a set of independent variable values X(I), I=1,...,N
C          and the corresponding dependent variables  Y(I), I=1,...,N 
C          DTNSI sets the input parameters for DTCGEN that define a 
C          natural spline interpolant to the data.
C
C  USAGE:
C          DOUBLE PRECISION X(N),Y(N),C(MAXC),HOLD(NHOLD)
C	   INTEGER N, MAXC, NHOLD, NC, IER
C
C          CALL DTNSI(N,X,Y,NDEG,ICC,HOLD,NHOLD,MAXC,C,NC,IER)
C
C  INPUT:
C          N       Number of data points, N .GE. MAX( 2, (NDEG+1)/2 )
C
C          X       Array of N values of the independent variable in
C                  ascending order.
C
C          Y       Array containing the values of the dependent
C                  variable.
C
C          NDEG    Degree of spline, NDEG .GE. 0.
C
C          ICC     Initial/continue call code. ICC is used to
C                  communicate whether this is an initial or
C                  continue call.  See Usage Remarks for further
C                  explanation.
C
C                  ICC .LT. 0  initial call to DTNSI.
C
C                  ICC .GT. 0  continue call TO DTNSI.
C
C          MAXC    Maximum storage allocated to C array.
C
C  STORAGE:
C          HOLD    Hold array of length NHOLD.  This array must be
C                  unchanged by the user between subsequent calls
C                  to DTNSI for ICC .GT. 0.
C
C          NHOLD   The length of array, HOLD,
C                  NHOLD .GE. NCOEF*(3*K-2) + 2*K**2 + 2*NCOEF 
C                             + 4*K + 9
C                  where NCOEF = (N-2) + K and K = NDEG + 1.
C
C  OUTPUT:
C          ICC     Initial/continue call code:
C
C                  ICC .GT. 0  Normal return, ICC has been set to the
C                              absolute value of ICC on input.  DTNSI
C                              may be called again for a different set
C                              of values for the dependent variable.
C
C                  ICC .EQ. 0  An error has been detected, IER .NE. 0.
C                              Check value of IER.
C
C          NC      Length of spline array.
C                  NC = 5 + (NDEP + 1)*NCOEF + K.
C                  where NDEP is the number of dependent variables.
C
C          C       Array of elements containing the information needed
C                  to evaluate the spline.  See DTSPVL abstract for a
C                  description of the contents of C.
C
C          IER     Success/error code.  DTNSI calls the standard error
C                  handler DTERR to print out error and warning messages
C                  for IER .NE. 0. For IER .LT. 0, DTNSI sets C(1) =
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
C                  IER = -2   N .LT. MAX( 2, (NDEG + 1)/2 ).
C
C                  IER = -3   NHOLD too small, the number of elements
C                             needed is given by the printed
C                             error message.
C
C                  IER = -5   X(I) not in ascending order.
C
C                  IER = -11  Invalid value for ICC. This may be caused
C                             by failing to reset ICC following an
C                             unsuccessful call to DTNSI.
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
      PARAMETER  (SUBNAM = 'DTNSI   ')
C
C  Arguments:
C
      INTEGER  ICC, IER, N, NDEG, NHOLD, MAXC
      DOUBLE PRECISION  C(MAXC), X(N), Y(N), HOLD(NHOLD)
C
C  Internal:
C
      INTEGER  I,K,MODE,NCOEF,NEED,NMIN,NRHS,NLHS,MLT
      DOUBLE PRECISION  YDUM
C
C=======================================================================
C     Set default values.
C=======================================================================
C
      IER = 0
      NEED = 0
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
C    Check NDEG.
C=======================================================================
C
      IF( NDEG .LT. 1 ) THEN
        MODE = 1
        IER = -1
        CALL DTERR(MODE,SUBNAM,IER,NEED)
        RETURN
      ENDIF
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
C             Set interior conditions.
C
      MLT = 1
C
C             Set rhs conditions.
C
      NRHS = 0
C
C             Construct spline array.
C
      CALL DTCGEN(NDEG,X,N,Y,MLT,YDUM,NLHS,YDUM,NRHS,ICC,
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
