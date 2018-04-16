C=======================================================================
C
      SUBROUTINE DTLSA (NPTS,T,Y,NROWY,NCURV,NDEG,IWT,WHT,IEND,
     &                  XKNOTS,NKNOTS,HOLD,NHOLD,MAXC,C,NC,IFAIL,IER)
C
C  PURPOSE:
C          DTLSA computes a weighted least squares spline fit to user-
C          specified curve data.  The number of dependent variables is
C          arbitrary.
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
C          where Y(I) - S(T(I)) means the Euclidean distance between
C          data point I (coordinates Y(I,*)) and the curve evaluated
C          at the parametric variable associated with that data point
C          (typically a cumulative chord length for lack of better
C          information).
C
C          Note that the internal knots need not be distinct.  Thus,
C          one can fit splines with user-specified multiple knots.
C
C          This version has an option to force the curve through the
C          first and last data points.  NDEG + 1 multiplicity for the
C          first and last knots is imposed here in either case - only
C          the INTERNAL knots are supplied by the user - but this
C          multiplicity alone does not constrain the end points of the
C          fitted curve.  If these are to match the end data points,
C          only the interior control points ("coefficients") should be
C          calculated, and this is done by working with only the interior
C          columns of the (rectangular) matrix and a suitably adjusted
C          right-hand side.  The effect of this on the Normal Equations
C          is handled by zeroing the first/last column elements of A so
C          the interior part of A'A is as it should be (even in packed
C          band form), and subtracting from b the first and last columns
C          of A multiplied by the constrained values of the first and
C          last coefficients, prior to forming A'b.
C          
C
C  USAGE:
C          DOUBLE PRECISION T(N),Y(NROWY,NCURV),WHT(N)
C          DOUBLE PRECISION XKNOTS(NKNOTS),C(MAXC),HOLD(NHOLD)
C          CALL DTLSA (N,T,Y,NROWY,NCURV,NDEG,IWT,WHT,IEND,
C          XKNOTS,NKNOTS,HOLD,NHOLD,MAXC,C,NC,IFAIL,IER)
C
C  INPUT:
C          NPTS    Number of data points, NPTS .GE. NKNOTS + K.
C
C          T       Array of NPTS values of the independent variable T
C                  in ascending order.
C
C          Y       Array of NPTS values of the dependent variable in
C                  correspondence with the T array.
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
C          IEND    End-point option.
C
C                  IEND .EQ. 0  means the end-points of the curve are
C                               not constrained - they are calculated
C                               along with the rest of the coefficients.
C
C                  IEND .NE. 0  means the end-points of the curve are
C                               constrained as the end data points.
C
C          NKNOTS  Number of values in the XKNOTS array, NKNOTS .GE. 0.
C
C          XKNOTS  Array of NKNOTS values containing the internal
C                  knots in ascending order with multiplicity .LE. K.
C                  (These are typically derived from the data points
C                  at the application level.  At least one data point
C                  should be contained by each knot span.)
C
C          MAXC    Storage allocated to spline c-vector.
C
C  STORAGE:
C          HOLD    Hold array of length NHOLD.
C
C          NHOLD   Length of array HOLD.  NHOLD .GE. NCOEF*(K + 1) + 2*K
C                  where NCOEF = NKNOTS + K.
C
C  OUTPUT:
C
C          C       Array of elements containing the information needed
C                  to evaluate the spline.  See DTSPVL abstract for a
C                  description of the contents of C. 
C
C          NC      Length of spline array.  NC = 5 + (NDEP + 1)*NCOEF + K
C                  where NDEP is number of dependent variables.
C
C          IFAIL   Knot interval where interlacing fails.  Used only if
C                  IER = -35.
C
C          IER     Success/error code as returned from DTLSAA, q.v.
C                  (Lower level routines use the standard error
C                  handler DTERR to print error and warning messages
C                  for IER .NE. 0.)  For IER .LT. 0 DTLSA sets C(1) =
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
C                  IER = -20   The matrix representing the least squares
C                              problem is singular.
C
C                  IER = -35   XKNOTS failed interlacing conditions.
C
C                  IER = -36   XKNOTS(I) .LT. T(1) or
C                              XKNOTS(I) .GT. T(NPTS)
C
C                  IER = -37   T(K) .LT. T(1) .OR. T(K) .GT. T(NPTS)
C                              for K = 2,...,NPTS-1.
C
C                  IER = -100  Unexpected error return from DTILCK. See
C                              explanation of MODE = 5 error message
C                              from DTERR.
C
C
C  OTHER DT ROUTINES CALLED
C
C     DTLSAA       Called once for each dimension.  (The dimensions are
C                  separable but use the same left-hand side matrix,
C                  which is factorized only once per call to DTLSA.)
C
C  HISTORY
C
C     198x         Boeing              Original code (author(s) not shown.
C                                      No end-point constraint option.
C
C     03/23/92     David Saunders      Introduced IEND to allow optional
C                  Sterling Software/  constraining of the curve end points.
C                  NASA Ames           Clarified the description somewhat.
C                                      NHOLD length NCOEF*K + 3*K seems
C                                      to be in error because DPBCO needs
C                                      a work vector of length NCOEF.
C
C=======================================================================

C     Arguments:

      INTEGER  NPTS, NROWY, NCURV, NDEG, IWT, IEND, NKNOTS, NHOLD, MAXC,
     &         NC, IFAIL, IER
      DOUBLE PRECISION  C(MAXC), T(NPTS), WHT(NPTS), XKNOTS(*),
     &                  Y(NROWY,NCURV), HOLD(NHOLD)

C     Local variables:

      INTEGER  I, ICC

C     Execution:

      ICC = -1
      DO 10, I = 1, NCURV
         CALL DTLSAA (NPTS,T,Y(1,I),NDEG,ICC,IWT,WHT,IEND,
     &                XKNOTS,NKNOTS,HOLD,NHOLD,MAXC,C,NC,IFAIL,IER)
         ICC = 2
   10 CONTINUE

      RETURN
      END
C=====================================================================
C
      SUBROUTINE DTLSA1 (NPTS,T,Y,IWT,WHT,IEND,NCOEF,KORD,XKNOTS,WORK,
     &                   ATA)

C     Purpose:
C
C     DTLSA1 sets up the LHS of the Normal Equations for the system
C     solved by DTLSAA.  LINPACK's symmetric positive definite banded
C     matrix storage convention is used.  This version provides for
C     constraining the end points of the fitted spline.
C
C     Method:
C
C     The underlying rectangular matrix A is irregularly block diagonal,
C     but A'A is banded with "NDEG" diagonals above the main diagonal.
C     Equations involving knot span [T(L), T(L+1)) contribute to the
C     elements I, J = L - K - 1 through L - 1 in A'A, where K = NDEG + 1
C     (for the unconstrained case).  Constraining the end-points means
C     only the interior control points ("coefficients") need to be
C     solved for.  This is achieved here by zeroing the first and last
C     columns of A, which corresponds to zeroing the first and last
C     rows and columns of A'A (apart from the main diagonal elements).
C     The system solved remains of order NCOEF in all cases.
C
C     03/24/92  D.A.Saunders  IEND argument & explanation introduced.
C
C=====================================================================

C     Arguments:

      INTEGER  IWT, IEND, KORD, NCOEF, NPTS
      DOUBLE PRECISION    ATA(KORD,*), WORK(*), T(NPTS), WHT(NPTS),
     &                    XKNOTS(*), Y(NPTS)

C     Local constants:

      PARAMETER  (ZERO=0.D0, ONE=1.D0)

C     Local variables:

      INTEGER  I, IK, IPOS, IR, I1, I2, J, JC, K
      DOUBLE PRECISION    W

C=====================================================================
C     Initialize.
C=====================================================================

      DO 20 J = 1, NCOEF
        DO 10 I = 1, KORD
          ATA(I,J) = ZERO
   10   CONTINUE
   20 CONTINUE

C=======================================================================
C     Loop through points defining contribution to ATA.
C=======================================================================

      I1 = 1 + KORD
      I2 = I1 + KORD
      IPOS = KORD
      DO 70 K = 1, NPTS

C        Determine position parameter.

   30    IF (T(K) .LE. XKNOTS(IPOS+1)) GO TO 40
         IPOS = IPOS + 1
         GO TO 30

C        Determine if point is to be used.

   40    W = ONE
         IF (IWT .NE. 0) THEN
            W = WHT(K)
            IF (W .LE. ZERO) GO TO 70
            W = W * W
         END IF

C        Evaluate the B-spline basis functions.

         IK = -1
         CALL DTBSP2 (XKNOTS,T(K),IPOS,IK,KORD,WORK(I1),WORK(I2),WORK)

C        Add contribution of this point.

         JC = IPOS - KORD
         DO 60 I = 1, KORD
            JC = JC + 1
            DO 50 J = 1, I
               IR = J - I + KORD
               ATA(IR,JC) = ATA(IR,JC) + WORK(J)*WORK(I)*W
   50       CONTINUE
   60    CONTINUE
   70 CONTINUE

      IF (IEND .NE. 0) THEN

C        The first and last coefficients are being constrained.
C        Zero the first and last columns of A, which leads to zeroing
C        these rows and columns of A'A.  The first and last main diagonal
C        elements are left alone - they just need to be nonzero.

         JC = 1
         DO 80, I = KORD - 1, 1, -1
            JC = JC + 1
            ATA(I,JC) = ZERO    ! First row (& column)
            ATA(I,NCOEF) = ZERO ! Last column (& row)
   80    CONTINUE

      END IF

C     <Should check that enough points are active.>

      RETURN
      END
C=====================================================================
C
      SUBROUTINE DTLSA2 (NPTS,T,Y,IWT,WHT,IEND,NCOEF,KORD,XKNOTS,WORK,
     &                   ATB,IER)

C     Purpose:
C
C     DTLSA2 sets up a RHS vector of the Normal Equations for the
C     system solved by DTLSAA.  See DTLSA1 for an explanation of the
C     constrained end-point case and the added IEND argument.
C
C=====================================================================

C     Arguments:

      INTEGER IWT, IEND, IER, KORD, NCOEF, NPTS
      DOUBLE PRECISION    ATB(*), WORK(*), T(*), WHT(*)
      DOUBLE PRECISION    XKNOTS(*), Y(*)

C     Local constants:

      PARAMETER  (ZERO=0.D0, ONE=1.D0)

C     Local variables:

      INTEGER  I, IK, IPOS, I1, I2, J, JC, K, NEND
      DOUBLE PRECISION    ADJUST, W, WW

C=====================================================================
C     Initialize.
C=====================================================================

      NEND  = NCOEF + 1
      DO 10 J = 1, NCOEF
         ATB(J) = ZERO
   10 CONTINUE

C=======================================================================
C     Loop through points defining contribution to A'b.
C=======================================================================

      I1 = 1 + KORD
      I2 = I1 + KORD
      IPOS = KORD
      DO 50 K = 1, NPTS

C        Verify T(K) is still in range.

         IF (XKNOTS(KORD).GT.T(K) .OR. T(K).GT.XKNOTS(NEND)) THEN
            IER = -37
            RETURN
         ENDIF

C        Determine position parameter.

   20    IF (T(K) .LE. XKNOTS(IPOS+1)) GO TO 30
         IPOS = IPOS + 1
         GO TO 20

C        Determine if point is to be used.

   30    W = ONE
         IF (IWT .NE. 0) THEN
            W = WHT(K)
            IF (W .LE. ZERO) GO TO 50
            W = W * W
         END IF

C        Evaluate the B-spline basis functions.

         IK = -1
         CALL DTBSP2 (XKNOTS,T(K),IPOS,IK,KORD,WORK(I1),WORK(I2),WORK)

         ADJUST = ZERO
         IF (IEND .NE. 0) THEN  ! First/last coefficients are constrained

C           Adjust the RHS, Y (not A'Y).  In the degenerate case, both
C           adjustments can affect the same element.

            IF (IPOS .EQ. KORD)  ADJUST = Y(1)*WORK(1)
            IF (IPOS .GE. NCOEF) ADJUST = Y(NPTS)*WORK(KORD) + ADJUST
         END IF

C        Add contribution of this point to A'b:

         JC = IPOS - KORD
         DO 40 I = 1, KORD
            JC = JC + 1
            ATB(JC) = ATB(JC) + WORK(I)*(Y(K) - ADJUST)*W
   40    CONTINUE
   50 CONTINUE

C     <Should check that enough points are active.>

      RETURN
      END
C=======================================================================
C
      SUBROUTINE DTLSAA (NPTS,T,Y,NDEG,ICC,IWT,WHT,IEND,XKNOTS,NKNOTS,
     &                   HOLD,NHOLD,MAXC,C,NC,IFAIL,IER)
C
C
C  PURPOSE:
C          DTLSAA computes a weighted least squares spline fit to
C          the coordinates in one dimension of a set of points
C          defining a curve.  It is called once per coordinate
C          dimension by applications such as DTLSA, since the
C          calculations for each dimension are separable.
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
C          This version has an option to force the curve through the
C          first and last data points.  NDEG + 1 multiplicity for the
C          first and last knots is imposed here in either case - only
C          the INTERNAL knots are supplied by the user - but this
C          multiplicity alone does not constrain the end points of the
C          fitted curve.  If these are to match the end data points,
C          only the interior control points ("coefficients") should be
C          calculated.  This is done by zeroing the first and last
C          columns of the (rectangular) matrix and suitably adjusting
C          the right-hand side.  See DTLSA for more details.
C
C  USAGE:
C          DOUBLE PRECISION T(N),Y(N),WHT(N),XKNOTS(NKNOTS),C(MAXC),
C                           HOLD(NHOLD)
C          CALL DTLSAA (N,T,Y,NDEG,ICC,IWT,WHT,IEND,XKNOTS,NKNOTS,HOLD,
C                       NHOLD,MAXC,C,NC,IFAIL,IER)
C
C  INPUT:
C          NPTS    Number of data points, NPTS .GE. NKNOTS + K.
C
C          T       Array of NPTS values of the independent variable T
C                  in ascending order.
C
C          Y       Array of NPTS values of the dependent variable in
C                  correspondence with the T array.
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
C          IEND    End-point option.
C
C                  IEND .EQ. 0  means the end-points of the curve are
C                               not constrained - they are calculated
C                               along with the rest of the coefficients.
C
C                  IEND .NE. 0  means the end-points of the curve are
C                               constrained as the end data points.
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
C                  which must be preserved:
C                                          LENGTH
C                  1   Band matrix         K*NCOEF
C                  I2  Work area           2*K + NCOEF
C
C          NHOLD   Length of array HOLD.  NHOLD .GE. NCOEF*(K + 1) + 2*K
C                  where NCOEF = NKNOTS + K.
C
C  OUTPUT:
C          ICC     Initial/continue call code:
C
C                  ICC .GT. 0   Normal return.  ICC has been set to the
C                               absolute value of ICC on input.  DTLSAA
C                               may be called again for a different
C                               set of dependent variable data.
C
C                  ICC .EQ. 0   An error has been detected.  IER .NE. 0.
C                               Check value of IER.
C
C          C       Array of elements containing the information needed
C                  to evaluate the spline.  See DTSPVL abstract for a
C                  description of the contents of C. 
C
C          NC      Length of spline array.  NC = 5 + (NDEP + 1)*NCOEF + K.
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
C  OTHER ROUTINES CALLED
C
C     DPBCO
C     DPBSL
C     DTLSA1
C     DTLSA2
C     DTILCK
C     DTERR
C
C  HISTORY
C
C     198x       Boeing              Original code (author(s) not shown.
C                                    No end-point constraint option.
C
C     03/23/92   David Saunders      Introduced IEND to allow optional
C                Sterling Software/  constraining of the curve end points.
C                NASA Ames           Clarified the description somewhat.
C                                    NHOLD prescription seemed in error
C                                    (see DTLSA).
C=======================================================================

C     Arguments:

      INTEGER  NPTS, NDEG, ICC, IWT, IEND, NKNOTS, NHOLD, MAXC, NC,
     &         IFAIL, IER
      DOUBLE PRECISION  C(MAXC), T(NPTS), WHT(NPTS), XKNOTS(*),
     &                  Y(NPTS), HOLD(NHOLD)

C     Local constants:

      CHARACTER*8  SUBNAM
      PARAMETER  (ZERO=0.D0, ONE=1.D0, TWO=2.D0, THREE=3.D0,
     &            SUBNAM='DTLSA  ')

C     Local variables:

      DOUBLE PRECISION  DIGCAL, DIGMAX, DTMCON, RCOND
      INTEGER  INFO, I, I2, IL, K
      INTEGER  MODE, NACT, NCOEF, NCOLS, NEED, NTOTAL

C=======================================================================
C     Define integer parameters.
C=======================================================================

      IER = 0
      NEED = 0
      MODE = 1
      K = NDEG + 1
      NCOEF = NKNOTS + K
      I2 = 1 +  K*NCOEF

C=======================================================================
C     Check ICC.
C=======================================================================

      IF (ICC.GT.0)
     &  THEN

C         Check for error from previous matrix factorization.

          IF (C(1) .LT. ZERO) THEN
            IER = -12
            CALL DTERR (MODE,SUBNAM,IER,NEED)
            RETURN
          ENDIF

C         Check that input is consistent with C-vector.

          IF (NCOEF .NE. NINT(C(4))) THEN
            IER = -13
            RETURN
          ENDIF
          GO TO 190

        ELSEIF (ICC.EQ.0) THEN
          IER = -11
          CALL DTERR (MODE,SUBNAM,IER,NEED)
          RETURN
      ENDIF
      C(1) = -ONE

C=======================================================================
C     Check other arguments.
C=======================================================================

      IF (NDEG.LT.0) THEN
        IER = -1
        CALL DTERR (MODE,SUBNAM,IER,NEED)
        RETURN
      ENDIF

      IF (NPTS.LT.MAX(NCOEF,2)) THEN
        IER = -2
        CALL DTERR (MODE,SUBNAM,IER,NEED)
        RETURN
      ENDIF

      IF (IWT.NE.0) THEN
        NACT = 0
        DO 40 I = 1, NPTS
          IF (WHT(I).GT.ZERO)  NACT = NACT + 1
40      CONTINUE
        IF (NACT.LT.MAX(NCOEF,2)) THEN
          IER = -2
          CALL DTERR (MODE,SUBNAM,IER,NEED)
          RETURN
        ENDIF
      ENDIF

      NEED = NCOEF*(K + 1) + 2*K
      IF (NHOLD.LT.NEED) THEN
        IER = -3
        MODE = 2
        CALL DTERR (MODE,SUBNAM,IER,NEED)
        RETURN
      ENDIF

      IF (NKNOTS.LT.0) THEN
        IER = -6
        CALL DTERR (MODE,SUBNAM,IER,NEED)
        RETURN
      ENDIF

C     Check that T vector is increasing.

      DO 80 I = 2, NPTS
        IF (T(I) .LT. T(I-1) ) THEN
          IER = -5
          CALL DTERR (MODE,SUBNAM,IER,NEED)
          RETURN
        ENDIF
   80 CONTINUE

C=======================================================================
C     Set the C-vector.
C=======================================================================

      C(2) = ONE
      C(3) = DBLE (K)
      C(4) = DBLE (NCOEF)
      C(5) = C(3)
      IL = 5

C=======================================================================
C     Set knots.
C=======================================================================

C     Left endpoint:

      DO 90 I = 1, K
        IL = IL + 1
        C(IL) = T(1)
   90 CONTINUE

C     Interior:

      DO 110 I = 1, NKNOTS
        IF (T(1).LT.XKNOTS(I) .AND. XKNOTS(I).LT.T(NPTS))
     &    THEN
            IL = IL + 1
            C(IL) = XKNOTS(I)
          ELSE
            IER = -36
            CALL DTERR (MODE,SUBNAM,IER,NEED)
            RETURN
        ENDIF
  110 CONTINUE

C     Right endpoint.

      DO 130 I = 1, K
        IL = IL + 1
        C(IL) = T(NPTS)
  130 CONTINUE

C=======================================================================
C     Check interlacing.
C=======================================================================

      NTOTAL = NCOEF + K
      CALL DTILCK (NPTS,T,K,IWT,WHT,NTOTAL,C(6),IFAIL,IER)
      IF (IER.NE.0) THEN
        IER = -100
        MODE = 5
        CALL DTERR (MODE,SUBNAM,IER,NEED)
        RETURN
      ENDIF
      IF (IFAIL.NE.0) THEN
        IF( IFAIL .LT. -100 ) IER = -8
        IF( IFAIL .GT. 0 ) IER = -35
        CALL DTERR (MODE,SUBNAM,IER,NEED)
        RETURN
      ENDIF

C=======================================================================
C     Set up the Normal Equations matrix.
C=======================================================================

C     If the end-points are constrained, the first and last columns
C     of (rectangular) A need to be suppressed (IEND .NE. 0).

      CALL DTLSA1 (NPTS,T,Y,IWT,WHT,IEND,NCOEF,K,C(6),HOLD(I2),HOLD)

C=======================================================================
C     Factor the normal matrix via LINPACK.
C=======================================================================

      RCOND = ZERO

      CALL DPBCO (HOLD, K, NCOEF, NDEG, RCOND, HOLD(I2), INFO)

C     Test for exact singularity.

      IF (RCOND.EQ.ZERO) THEN
        IER  = -20
        MODE = 3
        CALL DTERR (MODE,SUBNAM,IER,NEED)
        RETURN
      ENDIF

      DIGMAX = -DLOG10( DTMCON(5) )
      DIGCAL = -THREE*DLOG10( RCOND )
      MODE   = 0
      IF ( DIGCAL .GE. DIGMAX ) THEN
        IF ( DIGCAL .LT. TWO*DIGMAX ) THEN ! Matrix is poorly conditioned.
          IER = 1
          CALL DTERR (MODE,SUBNAM,IER,NEED)
        ELSE                               ! Matrix is badly conditioned.
          IER = 2
          CALL DTERR (MODE,SUBNAM,IER,NEED)
        ENDIF
      ENDIF

C=======================================================================
C     Check MAXC.
C=======================================================================

      C(1) = ONE

  190 CONTINUE
      NDPOLD = NINT (C(2))
      IF (ICC.LT.2) 
     &  THEN
          NDEP = NDPOLD
        ELSE
          NDEP = NDPOLD + 1
      ENDIF
      IL = 6 + K + NDEP*NCOEF
      NC = IL + NCOEF - 1

      IF (NC.GT.MAXC) THEN
        IER = -14
        RETURN
      ENDIF
      C(2) = DBLE (NDEP)

C=======================================================================
C     Define right hand side of the Normal Equations.
C=======================================================================

      CALL DTLSA2 (NPTS,T,Y,IWT,WHT,IEND,NCOEF,K,C(6),HOLD(I2),C(IL),
     >             IER)

      IF (IER.LT.0) THEN
        CALL DTERR (MODE,SUBNAM,IER,NEED)
        RETURN
      ENDIF

C=======================================================================
C     Solve the Normal Equations for this coordinate.
C=======================================================================

      CALL DPBSL (HOLD, K, NCOEF, NDEG, C(IL))
      ICC = ABS (ICC)

      IF (IEND .NE. 0) THEN

C       Constrain the end control points:

        C(IL) = Y(1)
        C(NC) = Y(NPTS)
      END IF

      RETURN
      END
