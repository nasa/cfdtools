      SUBROUTINE DTLSA(NPTS,T,Y,NROWY,NCURV,NDEG,IWT,WHT,
     +                  XKNOTS,NKNOTS,HOLD,NHOLD,MAXC,C,NC,IFAIL,IER)
C
C=======================================================================
C
C  PURPOSE:
C          DTLSA computes a weighted least squares spline fit to
C          user specified curve data.  The number of dependent variables
C          is arbitrary.
C
C  METHOD:
C          The user supplies the abscissae, ordinates and weights
C          ( T(I), Y(I,J), W(I) ) I = 1,...,N, J = 1, NCURV, and a 
C          set of internal knots X(1),...,X(M) and a degree NDEG.  
C          This routine computes the coefficients of a spline S(T) 
C          of degree NDEG (order K = NDEG + 1) with internal knots 
C          X(1),...,X(M) which minimizes the error defined by
C
C                N
C          E = SIGMA W(I)*( Y(I) - S(T(I)) )**2
C               I=1
C
C          Note that the internal knots need not be distinct.  Thus,
C          one can fit splines with user specified multiple knots.
C
C  USAGE:
C          DOUBLE PRECISION T(N),Y(NROWY,NCURV),WHT(N)
C          DOUBLE PRECISION XKNOTS(NKNOTS),C(MAXC),HOLD(NHOLD)
C          CALL DTLSA(N,T,Y,NROWY,NCURV,NDEG,IWT,WHT,
C          XKNOTS,NKNOTS,HOLD,NHOLD,MAXC,C,NC,IFAIL,IER)
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
C          NROWY   Dimension constant for Y array.
C
C          NCURV   Number of dependent variables.
C
C          NDEG    Degree of the spline, NDEG .GE. 0.
C
C          IWT     Weight option.
C
C                  IWT .EQ. 0  DTLSA uses WHT(I) = 1.0 for all I and
C                              the WHT array is not used.
C
C                  IWT .NE. 0  DTLSA uses weights provided in the
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
C                  DTLSA for ICC .GT. 0.  It contains information
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
C          IER     Success/error code. DTLSA calls the standard error
C                  handler DTERR to print error and warning messages
C                  for IER .NE. 0.  For IER .LT. 0 DTLSA sets C(1) =
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
C                  IER = -11   Invalid value of ICC.  This error is 
C                              unexpected.  See DTLSAA documentation.
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
C
C  OTHER DT ROUTINES CALLED
C
C     DPBCO
C     DPBSL
C     DTLSAA
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
      EXTERNAL DTMCON
      DOUBLE PRECISION ZERO, ONE, C2P0, C3P0
      PARAMETER  (ZERO=0.D0,ONE=1.D0,C2P0=2.D0,C3P0=3.D0,
     &            SUBNAM='DTLSA  ')
C
C  Arguments:
C
      INTEGER  IWT, NDEG, NKNOTS, NPTS, NHOLD, IFAIL, IER
      DOUBLE PRECISION  C(MAXC), T(NPTS), WHT(NPTS), XKNOTS(*),
     &                  Y(NROWY,NCURV), HOLD(NHOLD)
C
C  Internal:
C
      DOUBLE PRECISION  DTMCON
      INTEGER  ICC, I, I2, K, KM1
      INTEGER  MODE, NCOEF, NEED
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
      ICC = -1
      DO 10 I = 1, NCURV
      CALL DTLSAA(NPTS,T,Y(1,I),NDEG,ICC,IWT,WHT,
     +                  XKNOTS,NKNOTS,HOLD,NHOLD,MAXC,C,NC,IFAIL,IER)
      ICC = 2
10    CONTINUE
      RETURN
      END
