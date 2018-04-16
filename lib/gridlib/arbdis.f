C+------------------------------------------------------------------------------
C
      SUBROUTINE ARBDIS (NX, XA, XB, NC, XC, YC, START, LUNOUT, WORK,
     >                   X, IER)
C
C  ONELINER: ARBitrary DIStribution of points (1-D)
C
C  PURPOSE:  ARBDIS generates a 1-D distribution of points X(1:NX) in
C            the interval [XA, XB] with smoothly varying spacing pro-
C            portional in some sense to a given control  distribution
C            ("shape" function) defined at XA and XB (preferably) and
C            some number (possibly small) of interior points.   These
C            control points must have positive ordinates, but scaling
C            is unimportant. A uniform distribution would result from
C            any constant ordinates YC(*) (say 1.0 or some other pos-
C            itive  constant).   The smaller the (relative) value  of
C            YC(I) the smaller the spacing near X = XC(I) in the out-
C            put distribution X(*).   Note  that  "X" may well be the
C            usual approximation to arc length.
C
C  METHOD:   Let  S(X)  represent the control distribution - suitably
C            splined (shape function; monotonic spline is essential).
C            We want
C                        dX(I) ~ S( X(I) )        or (more precisely)
C                        dX(I) ~ S( XBAR(I) )     or
C                X(I+1) - X(I) = P * S( XBAR(I) )
C
C            for I = 1, ..., NX-1,  where XBAR(I) = (X(I) + X(I+1))/2
C            and P is an unknown constant of proportionality.
C
C            The NX = 4 case illustrates how this represents a system
C            of NX - 1 nonlinear equations:
C
C                F(1)  =  X(2) -  XA  -  P * S( XBAR(1) )  =  0
C                F(2)  =  X(3) - X(2) -  P * S( XBAR(2) )  =  0
C                F(3)  =   XB  - X(3) -  P * S( XBAR(3) )  =  0
C
C            We have three equations in three unknowns X(2), X(3), P.
C            Adding the equations gives a simple expression for P for
C            a given estimate of X(*):
C
C                P = (XB - XA) / Sum (I=1:NX-1) S( XBAR(I) )
C
C            The standard Newton iteration in-line can take advantage
C            of the sparse nature of the first derivatives of F(I)  -
C            they are +/-1 - 0.5 * P * S'( XBAR(I) ) or -S( XBAR(I) )
C            or zero. (A look at the piecewise cubic evaluated at the
C            mid-point between X(I) and X(I+1) shows that the partial
C            derivatives of S w.r.t. X(I) and w.r.t. X(I+1)  are both
C            just 0.5 * S'( XBAR(I) ) - a happy result.)
C                 
C            Provision is made for supplying a starting guess for the
C            desired X(I)s, but in practice it appears that a uniform
C            distribution is often adequate for the iteration to con-
C            verge. The steps taken at each iteration are safeguarded
C            to (almost) ensure that the Xs remain monotonic and that
C            the norm of the residual vector  F is decreasing.   How-
C            ever,  the step length can be halved only so many  times
C            without underflowing, especially in single precision, so
C            extreme cases may still fail to converge.
C
C  ARGUMENTS:
C    ARG    TYPE  I/O/S  DIM   DESCRIPTION
C     NX      I     I     -    Number of pts. in desired distribution; NX > 2
C     XA,     R     I     -    Coordinates of end pts. of desired distribution;
C      XB                      X(NX) = XB > XA = X(1)
C     NC      I     I     -    Number of control pts. supplied; NC > 1
C     XC      R     I    NC    Abscissas of control points in [XA, XB] and
C                              montonically increasing.  To avoid extrapolation,
C                              XC(1) = XA and XC(NC) = XB is recommended, but
C                              this is not essential.
C     YC      R     I    NC    Control values roughly proportional to the
C                              desired spacing of X(*) near XC(*).  Must be
C                              greater than zero, but the scaling is not
C                              important.  If all YC(*) = constant > 0.,
C                              X(*) will be uniform; smaller YC(I) gives
C                              smaller spacing near X = XC(I).
C    START   C*1    I     -    'A' for 'AUTO' means ARBDIS generates a uniform
C                              distribution as the starting-guess for X(*);
C                              'I' for 'INPUT' means X(*) is input with an
C                              estimate - not normally needed for the method
C                              to converge except on extreme cases.  Any input
C                              for START other than 'I' defaults to 'AUTO'.
C    LUNOUT   I     I     -    Logical unit for displaying convergence history.
C                              Set LUNOUT < 0 to suppress this output.
C    WORK     R     S    5*NX  Work-space used internally.  (This version
C                              avoids storing spline coefficients.)
C     X       R    I/O    NX   Desired distribution of abscissas.  See START
C                              if a starting guess is available.
C    IER      I     O     -    Error return code:
C                              IER = 0 means no error;
C                                  = 1 means NC < 2;
C                                  = 2 means NX < 3;
C                                  = 3 means not all YC(*)s are positive;
C                                  = 4 means START = 'I(nput)' but X(1) or
C                                      X(NX) does not match XA or XB properly;
C                                  = 5 means the Newton iteration diverged -
C                                      try an(other) input estimate for X(*).
C
C  PROCEDURES:
C     BINSLV     Solves the specialized linear system arising from the
C                particular nonlinear system produced by the ARBDIS algorithm
C     COPY       Used to keep copies of X and dX in case steps are cut back
C     LCSFIT     Local cubic spline fit & evaluation (S and S')
C
C  USAGE NOTES:
C
C     Note that there is no guarantee that final distribution X(*) will
C     include the XC(*) coordinates exactly.  But routine SMOOTHX is
C     available to perturb a given distribution so as to include any
C     specified points precisely.
C
C  ENVIRONMENT:  FORTRAN 77 with a few Fortran 90 upgrades.
C
C  HISTORY:
C  05/20/87    DAS    Initial implementation: simple-minded - gave the
C                     right general shape, but could not guarantee the
C                     precise location of dense regions, etc. - it used
C                     evaluations at the uniform Xs only - effectively
C                     just the starting guess for the eventual iteration.
C  07/10/87    DAS    Reformulated as the solution to a system of non-
C                     linear equations - unsymmetric form, with dX(I)
C                     proportional to S( X(I) ); demonstrated with general
C                     purpose nonlinear system solver from IMSL.
C  07/17/87    DAS    Took advantage of the very sparse Jacobian involved,
C                     solving the nonlinear system in-line by the usual
C                     Newton iteration.  Main worries: is the starting guess
C                     always good enough for convergence? and is the non-
C                     mid-point approximation serious?
C  07/27/87    DAS    Spline and derivatives are now evaluated at interval
C                     mid-points.
C  10/10/88    DAS    Safeguarded the iteration step lengths to (almost)
C                     ensure monotonic Xs and descending ||F||s.  (Still
C                     have to limit the number of times a step can be halved.)
C  08/25/89    DAS    Replaced MSFIT/CSDVAL with LCSFIT to avoid storing
C                     the spline coefficients.
C  09/05/89  CLH/DAS  IER = 0 was missing.
C  03/14/02    DAS    Allowed the search direction to reverse in the hope of
C                     achieving convergence for extreme cases; however, the
C                     need for X(*) to increase monotonically at every iteration
C                     must interfere with normal Newton iteration safeguarding.
C                     Also: checked for non-positive YC(*) values, changing
C                     IER usage a little.
C  06/29/06    DAS    Raised ITMAX from 600 to 2000 for difficult cases.
C  09/03/11     "     Todd White encountered what turned out to be an infinite
C                     loop introduced with the reverse-direction logic - fixed.
C
C  AUTHOR:     David Saunders, Sterling Software/NASA Ames, Moffett Field, CA
C
C-------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER
     >   IER, LUNOUT, NC, NX
      REAL
     >   WORK (*), X (NX), XA, XB, XC (NC), YC (NC)
      CHARACTER
     >   START * 1

C     Local constants:

      INTEGER, PARAMETER ::
     >   ITMAX = 2000     ! Safeguarding can mean many small steps

      REAL, PARAMETER ::
     >   BIG   = 1.E+30,
     >   HALF  = 0.5E+0,
     >   ONE   = 1.E+0,
     >   TOL   = 1.E-6,
     >   ZERO  = 0.E+0

C     Local variables:

      INTEGER
     >   I, ID, IF, IL, IP, ITER
      REAL
     >   ALPHA, FNORM, FNORM0, P, RANGE, SUMSI, TERM
      LOGICAL
     >   OVERSHOT, REVERSED

C     Procedures:

      EXTERNAL
     >   BINSLV, COPY, LCSFIT

C     Execution:

      IER = 0

      IF (NC < 2) THEN
         IER = 1
         GO TO 999
      END IF

      IF (NX <= 2) THEN
         IER = 2
         GO TO 999
      END IF

      DO I = 1, NC
         IF (YC (I) <= ZERO) THEN
            IER = 3
            GO TO 999
         END IF
      END DO

C     Partition the work-space:

                         ! 1 : NX - 1 for spline evaluated at mid-intervals
      IF = NX            ! Offset for function vector f (RHS) and soln. dX
      IP = IF + NX       !   "     "    "prime" (spline derivs. at mid-pts.)
      IL = IP + NX       !   "     "    lower diagonal of Jacobian matrix
      ID = IL + NX       !   "     "    diagonal of Jacobian
      RANGE = XB - XA

      IF (START /= 'I') THEN

C        Generate the uniform distribution as the initial estimate.

         TERM = RANGE / REAL (NX - 1)
         DO I = 1, NX - 1
            X (I) = XA + REAL (I - 1) * TERM
         END DO

      ELSE

C        An estimate is input by the caller in X (*).

         IF (X (1) /= XA .OR. X (NX) /= XB) THEN
            IER = 4
            GO TO 999
         END IF

      END IF

C     Perform a safeguarded Newton iteration for the nonlinear system:

      ITER = 0
      REVERSED = .FALSE.

  300 CONTINUE

C        Evaluate the (splined) shape function and its first derivative
C        at the midpoints of the intervals defined by the current X (*).
C        Overwrite the mid-point coordinates with the spline values.

         X (NX) = XB             ! Temporary usage; X (NX) = P below.

         DO I = 1, NX - 1
            WORK (I) = HALF * (X (I) + X (I + 1))
         END DO

         CALL LCSFIT (NC, XC, YC, .TRUE., 'M', NX - 1, WORK (1),
     >                WORK (1), WORK (IP + 1))

         IF (ITER == 0) THEN

C           We need a starting guess for constant of proportionality "P",
C           treated as the (NX-1)th unknown.  Use XB - XA = P * sum S(I):

            SUMSI = ZERO
            DO I = 1, NX - 1
               SUMSI = SUMSI + WORK (I)
            END DO

            P = RANGE / SUMSI

C           Also, initialize the convergence test and step-length safeguard.

            FNORM0 = BIG
            ALPHA = ZERO          ! For printout purposes only
         END IF

C        Find the norm of the residual vector, F, which should converge to 0.
C        Set up the right-hand side for the system J dX = -F in the process:

         FNORM = ZERO
         DO I = 1, NX - 1
            WORK (IF + I) = X (I) - X (I + 1) + P * WORK (I)   ! -F (I)
            FNORM = MAX (FNORM, ABS (WORK (IF + I)))
         END DO

         X (NX) = P

C        Halve the step length until || F || is reduced.
C        Note saving of old X and dX in WORK (*) below (after first iter.).

         IF (FNORM > FNORM0) THEN

            IF (ALPHA > TOL * RANGE) THEN  ! Presume TOL*RANGE > 2*MACHEPS
               ALPHA = HALF * ALPHA

               DO I = 2, NX
                  X (I) = WORK (IL + I - 1) + ALPHA * WORK (ID + I - 1)
               END DO
               P = ABS (X (NX))
               GO TO 300

            ELSE IF (.NOT. REVERSED) THEN ! Direction must have been uphill
               REVERSED = .TRUE.
               ALPHA = ONE

               DO I = 2, NX
                  WORK (ID + I - 1) = -WORK (ID + I - 1) ! Reverse the direction
                  X (I) = WORK (IL + I - 1) - WORK (ID + I - 1) ! * ONE
               END DO
               P = ABS (X (NX))
               GO TO 300
            END IF

            GO TO 800
         END IF

         IF (LUNOUT > 0) THEN
            WRITE (LUNOUT, '(A, 1X, I4, 3X, A, 1X, 1P, E11.3,
     >         3X, A, E11.3, 3X, A, E11.3)')
     >         ' ARBDIS iter.:', ITER, '||f(x)||:', FNORM,
     >         'step:', ALPHA, 'P:', P
         END IF

         IF (FNORM < TOL * RANGE) GO TO 900

         FNORM0 = FNORM

C        Set up the left-hand side of J dX = -F for the next iteration:

         DO I = 1, NX - 1
            WORK (I) = -WORK (I)                               ! U (I)
            TERM = (HALF * P) * WORK (IP + I)
            WORK (ID + I) =  ONE - TERM                        ! D (I)
            WORK (IL + I) = -ONE - TERM                        ! L (I)
         END DO

C        Solve the system.  Solution dX overwrites the RHS -F:

         CALL BINSLV (NX - 1, WORK (IL + 1), WORK (ID + 1),
     >                WORK (1), WORK (IF + 1))

C        Update the solution X.  Save previous iteration and new dX first.

         CALL COPY (NX - 1, X (2),         WORK (IL + 1))
         CALL COPY (NX - 1, WORK (IF + 1), WORK (ID + 1))

         ALPHA = ONE
         IF (ITER == 0) ALPHA = HALF         ! Conservative initial step

  600    CONTINUE                            ! Start of monotonicity loop

            DO I = 2, NX
               X (I) = WORK (IL + I - 1) + ALPHA * WORK (ID + I - 1)
            END DO
            P = X (NX)

C           Make sure abscissas remain ordered and P is positive:

            OVERSHOT = P <= ZERO
            X (NX) = XB

            DO I = 2, NX - 1
               IF (X (I - 1) >= X (I) .OR. X (I) >= X (I + 1))
     >            OVERSHOT = .TRUE.
            END DO

            IF (OVERSHOT) THEN
               IF (ALPHA > TOL * RANGE) THEN
                  ALPHA = HALF * ALPHA
                  GO TO 600
               END IF
               GO TO 800
            END IF

         REVERSED = .FALSE.
         ITER = ITER + 1
         IF (ITER < ITMAX)
     >GO TO 300

  800 IER = 5
      IF (LUNOUT > 0) WRITE (LUNOUT, '(/, A, I2)')
     >   ' ARBDIS: No convergence.  IER:', IER

  900 X (NX) = XB

  999 RETURN

      END SUBROUTINE ARBDIS
C+----------------------------------------------------------------------
C
      SUBROUTINE BINSLV (N, L, D, U, X)
C
C ACRONYM: BIdiagonal plus Nth column system SoLVer
C          --              -                 - --
C
C PURPOSE: BINSLV solves one system  A x = b,  where A has a structure
C          made up of a diagonal, a subdiagonal, and column N.  Such a
C          structure arises in the solution of the nonlinear system of
C          equations for the abscissas in a  1-D distribution  defined
C          by an arbitrary set of control points.   See routine ARBDIS
C          for its use.
C
C METHOD:  The LU factorization is
C
C               | D       U |     | D         | | 1       V |
C               | L D     U |     | L D       | |   1     V |
C               |   L D   U |  =  |   L D     | |     1   V |
C               |     L D U |     |     L D   | |       1 V |
C               |       L U |     |       L V | |         1 |
C
C          where the L, D, U and V elements are not necessarily equal.
C
C          The right-hand-side is overwritten with the solution.
C
C ARGUMENTS:
C   ARG  DIM TYPE I/O/S DESCRIPTION
C    N    -    I    I   Order of the system (>=3)
C    L    N    R   I/S  Off-diagonals of A in L(2:N); L(1) is unused
C    D    N    R   I/S  Diagonals of A in D(1:N-1); D(N) is workspace
C    U    N    R   I/S  Nth column elements in U(1:N)
C    X    N    R   I/O  Input with right-hand side;
C                       output with the solution
C
C ERROR HANDLING:  None.  Pivoting should be unnecessary (not proven).
C
C ENVIRONMENT:  FORTRAN 77 with a few Fortran 90 upgrades.
C
C HISTORY:
C   07/15/87   DAS    Initial implementation (unit diagonal case).
C   07/27/87   DAS    Generalized for case where ARBDIS uses mid-points
C                     of the intervals rather than left-hand-end points.
C
C AUTHOR: David Saunders, Sterling Software/NASA Ames, Moffett Field, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER
     >   N
      REAL
     >   L (N), D (N), U (N), X (N)

C     Local constants:

      REAL, PARAMETER ::
     >   ONE = 1.E+0

C     Local variables:

      INTEGER
     >   I
      REAL
     >   DINV

C     Execution:

C     Overwrite U (*) with V (*) from the sparse upper triangular factor.
C     Solve the L y = b part while we're at it:

      DINV = ONE / D (1)
      U (1) = U (1) * DINV
      X (1) = X (1) * DINV
      D (N) = ONE

      DO I = 2, N
         DINV = ONE / D (I)
         U (I) = (U (I) - L (I) * U (I - 1)) * DINV
         X (I) = (X (I) - L (I) * X (I - 1)) * DINV
      END DO

      X (N) = X (N) / U (N)

C     Solve U x = y:

      DO I = 1, N - 1
         X (I) = X (I) - U (I) * X (N)
      END DO

      END SUBROUTINE BINSLV
