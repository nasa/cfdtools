C+------------------------------------------------------------------------------
C
      SUBROUTINE LSFIT2 (N, X, F, DEGREE, NEARPTS, G, IER)
C
C  ACRONYM:  Local Least Squares method for smoothing 1-D data
C
C  DESCRIPTION:
C
C        For given data points X(I), F(I), I = 1 : N, LSFIT2 calculates
C     smoothed function values G(I) via local weighted least squares fits
C     of polynomials.  The abscissas should be monotonically increasing
C     or decreasing for the distance-based weighting to work properly.
C
C  ARGUMENTS:
C     NAME   DIM TYPE I/O/S DESCRIPTION
C     N           I     I   Number of data points;  N >= NEARPTS >= 2
C     X       N   R     I   Abscissas, monotonic increasing or decreasing
C     F       N   R     I   Ordinates (presumably noisy)
C     DEGREE      I     I   Degree of polynomials for local weighted least
C                           squares fits;  DEGREE = 1 (more smoothing)
C                           through 4 (less smoothing);  2 or 3 are normally
C                           the best choices; actual degree used may be
C                           less than the input DEGREE if N is small.
C     NEARPTS     I     I   Number of points in the neighborhood of each
C                           data point to use for each weighted local fit;
C                           normally DEGREE + 2 <= M <= DEGREE + 4 where
C                           M = NEARPTS, unless N is small, when M may be
C                           adjusted so as not to exceed N.  NEARPTS <= 10.
C     G       N   R     O   Smoothed function values corresponding to X(I)
C     IER     -   I     O   0: No errors
C                           1: N < 2                 (too few data points)
C                           2: DEGREE < 1            (should we allow 0 here?)
C                           3: NEARPTS < DEGREE + 2  (no smoothing is feasible)
C
C  PROCEDURES:
C
C     PNWFIT  Weighted least squares polynomial fit
C     PNEVAL  Evaluation of polynomial and its derivatives
C
C  REFERENCE:
C
C     Foley, Thomas A.  "Scattered data interpolation and approximation
C     with error bounds, Computer Aided Geometric Design 3" (1986) p.163
C
C  ALGORITHM OUTLINE:
C
C     >  For I = 1: N
C        >  Identify m points nearest to X(I)  (m = NEARPTS).
C        >  For each point k of these m points, calculate a weight
C               w(k) = (1 - s(k)/s)  +  c
C           where s(k) = | X - X(I) | for kth neighbor X,  s = max s(k),
C           and  c ~ 1 / m serves to ensure that all m points contribute
C           non-trivially.
C        >  Solve the weighted least squares polynomial problem.
C        >  Evaluate the polynomial at X(I), giving G(I).
C
C     >  Degenerate cases are handled as far as possible by using
C            reduced values of DEGREE and NEARPTS if necessary.
C
C  ENVIRONMENT:  Fortran 90
C
C  HISTORY:
C
C     11/18/87-  D. Saunders  Initial implementation of LSFIT1, a 1-D analog
C     12/04/90                of the 2-D algorithm referenced above.
C                             LSFIT1 has the option to spline the G(I) values.
C     01/29/04    "    "      Fortran 90 upgrade of LSFIT1's method 'G' (only).
C
C  AUTHOR:  David Saunders, ELORET/NASA Ames Research Center, CA.
C
C-------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) ::
     >   N, DEGREE, NEARPTS

      INTEGER, INTENT (OUT) ::
     >   IER

      REAL, DIMENSION (N), INTENT (IN) ::
     >   X, F

      REAL, DIMENSION (N), INTENT (OUT) ::
     >   G

C     Local constants:

      REAL, PARAMETER ::
     >   ONE   = 1.E+0,
     >   THREE = 3.E+0

C     Local variables:

      INTEGER
     >   DEG, I, J, K, K1, K2, LENWORK, M

      REAL
     >   RMSDEV, RSMAX, WMIN, XI,
     >   PNCOEFS(0:DEGREE),         ! Polynomial coefs. for a local fit
     >   WEIGHT(NEARPTS),
     >   WORK((DEGREE + 3)*NEARPTS) ! Work-space for solving each
                                    ! local overdetermined system
C     Procedures:

      EXTERNAL
     >   PNEVAL, PNWFIT

C     Execution:

      IER = 0
      IF (N < 2)                IER = 1
      IF (DEGREE < 1)           IER = 2  ! 0 works but is rather extreme.
      IF (NEARPTS < DEGREE + 2) IER = 3  ! No smoothing is feasible.
      IF (IER /= 0) GO TO 999

      LENWORK = (DEGREE + 3)*NEARPTS

C     Allow for small N:

      DEG  = MIN (DEGREE, N - 1)
      M    = MIN (N, NEARPTS)
      WMIN = ONE / MAX (REAL (M), THREE)  ! Should this constant be tuned a bit?

C     Perform local least squares smoothing in neighborhood of each data pt.:

      K1 = 1                    ! Define the nearest M points to X (I).
      K2 = M                    ! These are they for X (1).
      XI = X (1)

      DO I = 1, N

         RSMAX = ONE / MAX (ABS (X(K1) - XI), ABS (X(K2) - XI))
         J = 0

         DO K = K1, K2
            J = J + 1
            WEIGHT(J) = (ONE - RSMAX * ABS (X(K) - XI)) + WMIN
         END DO

C        Solve the weighted least squares polynomial problem:

         CALL PNWFIT (M, X(K1), F(K1), WEIGHT, DEG, LENWORK, WORK,
     >                PNCOEFS, RMSDEV, IER)

C        Evaluate the polynomial and its derivative(s) at XI:

         CALL PNEVAL (DEG, PNCOEFS, 1, XI, G(I))

C        Adjust for the next polynomial fit:

         XI = X(MIN (I+1, N))          ! Won't be used if I is already N.

         DO WHILE (K2 < N)
            IF (ABS (X(K2 + 1) - XI) < ABS (X(K1) - XI)) THEN
               K1 = K1 + 1
               K2 = K2 + 1
            ELSE
               EXIT
            END IF
         END DO

      END DO ! Next local fit

  999 RETURN

      END SUBROUTINE LSFIT2
