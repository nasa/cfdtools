C+----------------------------------------------------------------------
C
      SUBROUTINE GEODIS (XA, XB, N, D1, CLUS, X, LUNERR, IER)
C
C  ACRONYM:  GEOmetric-type DIStribution
C            ---            ---
C  PURPOSE:
C        GEODIS generates abscissas on the interval [XA, XB] with the
C     (possibly modified) geometric relationship indicated below, where
C     dX (I) = X (I+1) - X (I)):
C
C     dX (I) = dX (I-1) * (1 + M * ((I-1) / (N-2)) ** CLUS)  (I = 2:N-1)
C
C     with     dX (1) = D1 = X (2) - XA
C                N-1
C     and   XA + Sum dX(I) = XB
C                I=1
C
C     The multiplier M is computed here to give the desired sum of the dXs.
C     The input exponent CLUS may be zero (giving true geometric spacing),
C     although CONDIS with M = 2 and XM = XA + D1 would be more efficient
C     in this case.
C
C        The spacing may be increasing or decreasing, depending on how big
C     D1 is compared with (XB - XA) / (N - 1).  This may be the main virtue
C     of GEODIS over EXPDIS2.  See also GEODIS2 for a form of two-sided
C     geometric clustering.
C
C  METHOD:
C        A bisection search is considered adequate, since use of a general
C     zero-finding utility would be rather cumbersome here, while a Newton
C     iteration is out of the question for lack of an analytic 1st derivative.
C     X (N) - XB is forced arbitrarily close to zero.  Convergence criteria
C     are tight because the last interval may be the smallest, yet it in
C     effect absorbs all of the final error when X (N) is forced to be XB
C     exactly.
C
C  ARGUMENTS:
C    NAME    DIM   TYPE I/O/S DESCRIPTION
C   XA        -     R     I   Desired X (1); passing X (1) here is safe
C   XB        -     R     I   Desired X (N); passing X (N) here is safe
C                             Note that XA = 0 and XB = 1 can provide one
C                             normalized distribution for multiple reuse.
C   N         -     I     I   Number of points; N > 2
C   D1        -     R     I   Desired size of first interval, X (2) - X (1);
C                             may be smaller or larger than uniform dX
C   CLUS      -     R     I   CLUS = 0. gives true geometric spacing;
C                                  > 0. amplifies the bunching;
C                                  < 0. reduces the bunching;
C                                  +1.0 is extreme; << -1.0 can fail;
C                                  rule of thumb: [-1.0 : +0.5]
C   X         N     R     O   Output distribution, with X (1) = XA and
C                             X (N) = XB, and X (2) - X (1) = D1.
C   LUNERR    -     I     I   LUNERR > 0 shows the iteration history;
C                             LUNERR < 0 suppresses it (but any error
C                                        message goes to -LUNERR)
C   IER       -     I     O   IER = 0 means no problems;
C                                   1 means the iteration failed - must
C                                     have been bad inputs
C
C  ERROR HANDLING:
C     See LUNERR and IER.  Failure handling is left to the calling program,
C     which can presumably reprompt for parameters and try again.
C
C  ENVIRONMENT:  VAX/VMS, FORTRAN 77 with IMPLICIT NONE extension
C
C  DEVELOPMENT HISTORY:
C     DATE   INITIALS   DESCRIPTION 
C   11/01/86    JEM     Initial design and code as CLUSGEO.
C   05/08/87  RGL/DAS   Description clarified (but still some puzzles).
C   03/11/88    DAS     Functionality simplified and clarified; a few
C                       refinements to the convergence test, etc.;
C                       safeguarded the initial bracketing iteration;
C                       renamed GEODIS (with GEODIS2 in mind).
C   03/13/89    DAS     Raised TOL from 1.E-06 to 3.E-06.
C
C  AUTHOR:  John E. Melton, NASA/Ames Research Ctr., Moffett Field, CA
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments.

      INTEGER
     >   IER, LUNERR, N
      REAL
     >   D1, CLUS, X (N), XA, XB

C     Local constants.

      INTEGER
     >   MXITER
      REAL
     >   ONE, TOL, ZERO
      PARAMETER
     >   (TOL    = 3.E-06,    ! Absolute test on dM and relative test on d F(M).
     >    MXITER = 30,        ! ~ 20 bisections typical.  Note: 2**-20 ~ 1.E-6.
     >    ZERO   = 0.E+0,     ! (MXITER also safeguards initial bracketing.)
     >    ONE    = 1.E+0)

C     Local variables.

      INTEGER
     >   I, ITER
      REAL
     >   DEL, FRAC, RML, RMM, RMOLD, RMR, RNM2, SCALE, XR, YL, YM, YR

C     Execution.

      IER = 0
      XR = XB                     ! Local copy in case X (N) is passed as XB.
      IF (XR .LE. XA .OR.
     >    D1 .LE. ZERO .OR.
     >    N .LE. 2) GO TO 800

      X (1) = XA
      X (2) = XA + D1
      RNM2 = ONE / FLOAT (N - 2)
      SCALE = MAX (ABS (XA), ABS (XR))   ! An end-point could be zero.

C     Set up the bisection search to find the correct multiplier:

      RML = ZERO                  ! M = 0  <->  uniform distribution
      YL = XA + D1 * (N - 1) - XR 

      IF ( YL .LE. ZERO ) THEN    ! M = 0 is a lower bound
         RMR = RNM2               ! Try M ~ 1/N for an upper bound
      ELSE                        ! M = 0 is an upper bound
         RMR = -RNM2
      END IF

      ITER = 0
   60 CONTINUE                    ! Indefinite loop to bracket solution
         DEL = D1

         DO 80 I = 2, N - 1
            FRAC = FLOAT (I - 1) * RNM2
            DEL = DEL * (ONE + RMR * FRAC ** CLUS)

C           DEL can overflow if M is too big (no guaranteed safeguard).
C           DEL can also go negative if M * FRAC ** CLUS < -1, which is
C           quite possible for CLUS < 0 and small I.  Therefore:

            IF (DEL .LE. ZERO) THEN  ! Bad news.
               IF (LUNERR .GT. 0) WRITE (LUNERR,
     >            '(//'' GEODIS:  Negative dX. M = '', E12.4)') RMR
               RMR = RMR * 0.1E+0    ! May help if this was first iteration
               GO TO 90
            END IF

            X (I+1) = X (I) + DEL
   80    CONTINUE

         YR = X(N) - XR
         IF (YR * YL .LE. ZERO) GO TO 100  ! Else double M and try again

            RML = RMR
            YL  = YR
            RMR = RMR + RMR

   90       ITER = ITER + 1
            IF (ITER .LT. MXITER ) GO TO 60
            GO TO 800

  100 IF (LUNERR .GT. 0)
     >   WRITE (LUNERR, '(//'' GEODIS:  Bisection Search''/'' ITER'',
     >   ''        RML           YL          RMM           YM'',
     >   ''          RMR           YR'')')

      RMOLD = RMR

      DO 120, ITER = 1, MXITER

         RMM = (RMR + RML) / 2
         DEL = D1

         DO 110, I = 2, N - 1
            FRAC = FLOAT (I - 1) * RNM2
            DEL = DEL * (ONE + RMM * FRAC ** CLUS)
            X (I+1) = X (I) + DEL
  110    CONTINUE

         YM = X (N) - XR

         IF (LUNERR .GT. 0) WRITE (LUNERR, '(1X, I2, 1P, 6E13.5)')
     >      ITER, RML, YL, RMM, YM, RMR, YR

C        M should converge to a value in ~(-1., 1.),
C        and the function of M should converge to zero:

         IF (ABS (YM) .LT. TOL * SCALE .AND.        ! X may not be ~ 1.
     >       ABS (RMM - RMOLD) .LT. TOL) GO TO 900  ! We know |M| < 1.

         RMOLD = RMM

         IF (YM * YL .GT. ZERO) THEN
            RML = RMM
            YL = YM
         ELSE
            RMR = RMM
            YR = YM
         END IF

  120 CONTINUE

C     Dropped through - bad news:

      WRITE (ABS (LUNERR), '(/'' GEODIS: Bisection failed.'')')

  800 WRITE (ABS (LUNERR),
     >   '('' Inputs to GEODIS must be bad:'', /,
     >     '' N = '', I5, /,
     >     '' XA, XB, D1 = '', 1P, 3E10.3, /,
     >     '' CLUS = '', E10.3)') N, XA, XR, D1, CLUS

      IER = 1

  900 X (N) = XR

      RETURN
      END
