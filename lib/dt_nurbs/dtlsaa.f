      SUBROUTINE DTLSAA(NPTS,T,Y,NDEG,ICC,IWT,WHT,XKNOTS,NKNOTS,HOLD,
     +                  NHOLD,MAXC,C,NC,IFAIL,IER)
C
C=======================================================================
C
C  PURPOSE:
C          DTLSAA computes a weighted least squares spline fit to
C          user specified data.
C
C  METHOD:
C          The user supplies the abscissae, ordinates and weights
C          ( T(I), Y(I), W(I) ) I = 1,...,N, a set of internal knots 
C          X(1),...,X(M) and a degree NDEG.  This routine computes
C          the coefficients of a spline S(T) of degree NDEG (order
C          K = NDEG + 1) with internal knots X(1),...,X(M) which
C          minimizes the error defined by
C
C                N
C          E = SIGMA W(I)*( Y(I) - S(T(I)) )**2
C               I=1
C
C          Note that the internal knots need not be distinct.  Thus,
C          one can fit splines with user specified multiple knots.
C
C  USAGE:
C          DOUBLE PRECISION T(N),Y(N),WHT(N),XKNOTS(NKNOTS),C(MAXC),
C                           HOLD(NHOLD)
C          CALL DTLSAA(N,T,Y,NDEG,ICC,IWT,WHT,XKNOTS,NKNOTS,HOLD,
C                      NHOLD,MAXC,C,NC,IFAIL,IER)
C
C  INPUT:
C          NPTS    Number of datapoints, NPTS .GE. NKNOTS + K.
C
C          T       Array of NPTS values of the independent variable T
C                  in ascending order.
C
C          Y       Array of NPTS values of the dependent variable in
C                   correspondence with the T array.
C
C          NDEG    Degree of the spline, NDEG .GE. 0.
C
C          ICC     Initial/continue call code.  ICC is used to
C                  communicate whether this is an initial or
C                  continue call.  See usage remarks for further
C                  explanation.
C
C                  ICC .LT. 0    Initial call to DTLSAA.
C
C                  ICC .GT. 0    Continue call to DTLSAA.
C                      .EQ. 1    Overwrite previous coefficients with 
C                                new values.
C                      .GT. 1    Append new coefficients to C vector.
C
C          IWT     Weight option.
C
C                  IWT .EQ. 0  DTLSAA uses WHT(I) = 1.0 for all I and
C                              the WHT array is not used.
C
C                  IWT .NE. 0  DTLSAA uses weights provided in the
C                              WHT array.
C
C          WHT     Array of NPTS nonnegative weights in correspondence
C                  to the T array, WHT(I) .GE. 0. If IWT .EQ. 0, WHT
C                  is not used and it may be a dummy argument.  If
C                  IWT .NE. 0, at least NKNOTS + K of the WHT(I) must
C                  be positive, WHT(I) .GT. 0.
C
C          NKNOTS  Number of values in the XKNOTS array, NKNOTS .GE. 0.
C
C          XKNOTS  Array of NKNOTS values containing the internal
C                  knots in ascending order with multiplicity .LE. K.
C
C          MAXC    Storage allocated to spline c-vector.
C
C  STORAGE:
C          HOLD    Hold array of length NHOLD.  This array must not be
C                  changed by the user between subsequent calls to
C                  DTLSAA for ICC .GT. 0.  It contains information
C                  which must be preserved.
C                                        LENGTH
C                  1   Band matrix         K*NCOEF
C                  I2  Work area            3*K
C
C          NHOLD   Length of array HOLD,
C                  NHOLD .GE. NCOEF*K + 3*K
C                  where NCOEF = NKNOTS + K.
C
C  OUTPUT:
C          ICC     Initial/continue call code:
C
C                  ICC .GT. 0   Normal return, ICC has been set to the
C                               absolute value of ICC on input.  DTLSAA
C                               may be called again for a different
C                               set of dependent variable data.
C
C                  ICC .EQ. 0   An error has been detected, IER .NE. 0,
C                               Check value of IER.
C
C          C       Array of elements containing the information needed
C                  to evaluate the spline.  See DTSPVL abstract for a
C                  description of the contents of C. 
C
C          NC      Length of spline array.
C                  NC = 5 + (NDEP + 1)*NCOEF + K.
C                  where NDEP is number of dependent variables.
C
C          IFAIL   Knot interval where interlacing fails.  Used only if
C                  IER = -35.
C
C          IER     Success/error code. DTLSAA calls the standard error
C                  handler DTERR to print error and warning messages
C                  for IER .NE. 0.  For IER .LT. 0 DTLSAA sets C(1) =
C                  -1.0.
C
C                  IER = 0     Success.
C
C                  IER =  1    Results computed but they are sensitive
C                              to inaccuracies in the entries of T, Y,
C                              and XKNOTS.
C
C                  IER =  2    Results computed but they are strongly
C                              sensitive to inaccuracies in the entries
C                              of  T, Y, and XKNOTS.
C
C                  IER = -1    NDEG .LT. 0.
C
C                  IER = -2    NPTS .LT. NKNOTS + K or the number of
C                              positive weights less than NKNOTS + K.
C
C                  IER = -3    NHOLD too small. The number of elements
C                              needed is given by the printed error
C                              message.
C
C                  IER = -5    T(I) not in ascending order.
C
C                  IER = -6    NKNOTS .LT. 0.
C
C                  IER = -8    XKNOTS(I) not in ascending order or the
C                              multiplicity of a knot exceeds the order
C
C                  IER = -11   Invalid value of ICC.  This was probably
C                              caused by calling DTLSAA after an error
C                              occured without resetting ICC.
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
C                  IER = -20   The matrix representing the least square
C                              problem is singular.
C
C                  IER = -35   XKNOTS failed interlacing conditions.
C
C                  IER = -36   XKNOTS(I) .LT. T(1) or
C                              XKNOTS(I) .GT. T(NPTS)
C
C                  IER = -37   T(K) .LT. T(1) .OR. T(K) .GT. T(NPTS)
C                              for K = 2,...,NPTS-1.  If this error
C                              occurs on a subsequent call then
C                              the C or T vectors were probably
C                              altered between calls.
C
C                  IER = -100  Unexpected error return from DTILCK, See
C                              explanation of MODE = 5 error message
C                              from DTERR.
C
C  USAGE REMARK:
C          The initial/continue call code allows the user to efficiently
C          compute several spline approximations based on the same
C          independent variable values and weight conditions.  For the
C          initial call, ICC .LT. 0, the user supplies all input 
C          arguments.  For subsequent calls, ICC .GT. 0, the user can
C          redefine the Y argument corresponding to the initial T,
C          IWT and WHT arguments.  The T, IWT and WHT arguments must 
C          not be changed between subsequent calls to DTLSAA.
C
C
C  OTHER DT ROUTINES CALLED
C
C     DPBCO
C     DPBSL
C     DTLSA1
C     DTLSA2
C     DTILCK
C     DTERR
C
C=======================================================================
C
C  Parameters:
C
      CHARACTER*8  SUBNAM
      DOUBLE PRECISION ZERO, ONE, C2P0, C3P0
      PARAMETER  (ZERO=0.D0,ONE=1.D0,C2P0=2.D0,C3P0=3.D0,
     &            SUBNAM='DTLSA  ')
C
C  Arguments:
C
      INTEGER  ICC, IWT, NDEG, NKNOTS, NPTS, NHOLD, IFAIL, IER
      DOUBLE PRECISION  C(MAXC), T(NPTS), WHT(NPTS), XKNOTS(*),
     &                  Y(NPTS), HOLD(NHOLD)
C
C  Internal:
C
      DOUBLE PRECISION  DIGCAL, DIGMAX, DTMCON, RCOND
      INTEGER  INFO, I, I2, IL, K, KM1
      INTEGER  MODE, NACT, NCOEF, NEED, NTOTAL
      INTEGER  NDPOLD, NDEP
C
C=======================================================================
C     DEFINE INTEGER PARAMETERS
C=======================================================================
C
      IER = 0
      NEED = 0
      MODE = 1
      K = NDEG + 1
      KM1 = K - 1
      NCOEF = NKNOTS + K
      I2 = 1 +  K*NCOEF
C
C=======================================================================
C    Check ICC.
C=======================================================================
C
      IF (ICC.GT.0)
     &  THEN
C
C             Check for error from previous matrix factorization.
C
          IF (C(1) .LT. ZERO) THEN
            IER = -12
            CALL DTERR(MODE,SUBNAM,IER,NEED)
            RETURN
          ENDIF
C
C             Check that input is consistent with c-vector.
C
          IF (NCOEF .NE. NINT(C(4))) THEN
            IER = -13
            RETURN
          ENDIF
          GO TO 190
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
      IF (NDEG.LT.0) THEN
        IER = -1
        CALL DTERR(MODE,SUBNAM,IER,NEED)
        RETURN
      ENDIF
C
C=======================================================================
C    Check NPTS.
C=======================================================================
C
      IF (NPTS.LT.MAX0(NCOEF,2)) THEN
        IER = -2
        CALL DTERR(MODE,SUBNAM,IER,NEED)
        RETURN
      ENDIF
C
C=======================================================================
C    Check number of positive weights.
C=======================================================================
C
      IF (IWT.NE.0) THEN
        NACT = 0
        DO 40 I = 1, NPTS
          IF (WHT(I).GT.ZERO)  NACT = NACT + 1
40      CONTINUE
        IF (NACT.LT.MAX0(NCOEF,2)) THEN
          IER = -2
          CALL DTERR(MODE,SUBNAM,IER,NEED)
          RETURN
        ENDIF
      ENDIF
C
C=======================================================================
C    Check NHOLD.
C=======================================================================
C
      NEED = NCOEF*K + 3*K
      IF (NHOLD.LT.NEED) THEN
        IER = -3
        MODE = 2
        CALL DTERR(MODE,SUBNAM,IER,NEED)
        RETURN
      ENDIF
C
C=======================================================================
C    Check NKNOTS.
C=======================================================================
C
      IF (NKNOTS.LT.0) THEN
        IER = -6
        CALL DTERR(MODE,SUBNAM,IER,NEED)
        RETURN
      ENDIF
C
C=======================================================================
C     Check that T vector is increasing.
C=======================================================================
C
      DO 80 I = 2, NPTS
        IF (T(I) .LT. T(I-1) ) THEN
          IER = -5
          CALL DTERR(MODE,SUBNAM,IER,NEED)
          RETURN
        ENDIF
80    CONTINUE
C
C=======================================================================
C    Set the C-vector.
C=======================================================================
C
      C(2) = ONE
      C(3) = DBLE( FLOAT(K) )
      C(4) = DBLE( FLOAT(NCOEF) )
      C(5) = C(3)
      IL = 5
C
C=======================================================================
C    Set knots.
C=======================================================================
C
C  Left endpoint:
C
      DO 90 I = 1, K
        IL = IL + 1
        C(IL) = T(1)
90    CONTINUE
C
C  Interior:
C
      DO 110 I = 1, NKNOTS
          IF (T(1).LT.XKNOTS(I) .AND. XKNOTS(I).LT.T(NPTS))
     &      THEN
              IL = IL + 1
              C(IL) = XKNOTS(I)
            ELSE
              IER = -36
              CALL DTERR(MODE,SUBNAM,IER,NEED)
              RETURN
          ENDIF
110   CONTINUE
C
C  Right endpoint.
C
      DO 130 I = 1, K
        IL = IL + 1
        C(IL) = T(NPTS)
130   CONTINUE
C
C=======================================================================
C     Call DTILCK to check interlacing.
C=======================================================================
C
      NTOTAL = NCOEF + K
      CALL DTILCK(NPTS,T,K,IWT,WHT,NTOTAL,C(6),IFAIL,IER)
      IF (IER.NE.0) THEN
        IER = -100
        MODE = 5
        CALL DTERR(MODE,SUBNAM,IER,NEED)
        RETURN
      ENDIF
      IF (IFAIL.NE.0) THEN
        IF( IFAIL .LT. -100 ) IER = -8
        IF( IFAIL .GT. 0 ) IER = -35
        CALL DTERR(MODE,SUBNAM,IER,NEED)
        RETURN
      ENDIF
C
C=======================================================================
C     Call DTLSA1 to set up normal matrix.
C=======================================================================
C
      CALL DTLSA1(NPTS,T,Y,IWT,WHT,NCOEF,K,C(6),HOLD(I2),HOLD)
C
C=======================================================================
C     Factor the normal matrix via LINPACK.
C=======================================================================
C
      RCOND = ZERO
      CALL DPBCO(HOLD,K,NCOEF,KM1,RCOND,HOLD(I2),INFO)
C
C  Test for exact singularity.
C
      IF (RCOND.EQ.ZERO) THEN
        IER  = -20
        MODE = 3
        CALL DTERR(MODE,SUBNAM,IER,NEED)
        RETURN
      ENDIF
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
            CALL DTERR(MODE,SUBNAM,IER,NEED)
C
          ELSE
C
C             Matrix is badly conditioned.
C
            IER = 2
            CALL DTERR(MODE,SUBNAM,IER,NEED)
        ENDIF
      ENDIF
C
C=======================================================================
C    Check MAXC.
C=======================================================================
C
      C(1) = ONE
190   CONTINUE
      NDPOLD = NINT(C(2))
      IF (ICC.LT.2) 
     &  THEN
          NDEP = NDPOLD
        ELSE
          NDEP = NDPOLD + 1
      ENDIF
      IL = 6 + K + NDEP*NCOEF
      NC = IL + NCOEF - 1
C
      IF (NC.GT.MAXC) THEN
        IER = -14
        RETURN
      ENDIF
      C(2) = DBLE( FLOAT(NDEP) )
C
C=======================================================================
C    Define right hand side and solve.
C=======================================================================
C
      CALL DTLSA2(NPTS,T,Y,IWT,WHT,NCOEF,K,C(6),HOLD(I2),C(IL),IER)
C
      IF (IER.LT.0) THEN
        CALL DTERR(MODE,SUBNAM,IER,NEED)
        RETURN
      ENDIF
C
      CALL DPBSL(HOLD,K,NCOEF,KM1,C(IL))
      ICC = IABS( ICC )
C
      RETURN
      END
