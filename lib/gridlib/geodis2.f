C+----------------------------------------------------------------------
C
      SUBROUTINE GEODIS2 (XA, XB, N, D1, CLUS1, CLUS2, X, LUNERR, IER)
C
C  ACRONYM:  GEOmetric-type DIStribution, 2-sided
C            ---            ---           -
C  PURPOSE:
C        GEODIS2 generates abscissas on the interval [XA, XB] using two
C     (possibly modified) geometric distributions of the type generated
C     by GEODIS, q.v.  Its main raison d'etre is its potential for UN-
C     symmetric bunching at the end points, with reasonable blending in
C     the middle of what are really two distinct distributions.  It can
C     also produce unsymmetric bunching around the mid-point, but one
C     of DSTRIB's may be preferable for this case.
C
C  METHOD:
C        The two given clustering exponents are used to proportion the
C     total number of points N in some reasonable way between the two
C     halves of the interval, as N1 and N2 say.  Then N1, D1, and CLUS1
C     defines a distribution on the first half, and N2, the last interval
C     from the first half, and CLUS2 define the rest of the points.
C
C        Ideally, if CLUS1=CLUS2, the result should be symmetric.  But
C     the only obvious way of achieving symmetry is, for odd N, to have
C     a repeated interval either side of the mid-point.  (If N is even,
C     another kind of outer iteration would be needed to get the middle
C     interval to straddle the mid-point properly - not attempted here.)
C
C        In general, having the last interval of the first distribution
C     be the first interval of the second is more appealing than having
C     two equal intervals in the middle.  Thus the symmetric case above
C     is treated as an exception.
C
C        The algorithm for determining N1 and N2 is somewhat arbitrary,
C     being basically
C                      N1 = N/2 + A (C1 - C2) 
C                      N2 = N/2 - A (C1 - C2)
C     for suitable choice of the constant A.  A = N/10 seems reasonable,
C     meaning each 0.1 difference in the clustering exponents adds 1%
C     of the points to one half and subtracts 1% from the other half.
C
C  DISCLAIMERS:
C        Until a single formulation is found, some composite approach is
C     the best we can do.  The lack of symmetry mentioned above for even
C     N and CLUS1=CLUS2 is troubling, but likely applications will have
C     different clustering parameters.
C
C        Another asymmetry is apparent in the explicit control of the
C     extreme interval at one end but not at the other.  Again, this may
C     be of little consequence: consider the intended initial application
C     to airfoil sections where the two surfaces are treated separately,
C     higher density is needed at the leading edge rather than at the
C     trailing edge, and some control over continuity between surfaces
C     is needed at the leading edge but not at the trailing edge...
C
C  ARGUMENTS:
C    NAME    DIM   TYPE I/O/S DESCRIPTION
C   XA        -     R     I   Desired X (1); passing X (1) here is safe.
C   XB        -     R     I   Desired X (N); passing X (N) here is safe.
C                             Note that XA = 0 and XB = 1 can provide one
C                             normalized distribution for multiple reuse.
C   N         -     I     I   Total number of points.
C   D1        -     R     I   Desired size of first interval, X (2) - X (1);
C                             may be smaller or larger than uniform dX
C   CLUS1     -     R     I   Clustering parameters for the first half
C   CLUS2     -     R     I   and second half of [XA, XB].  The difference
C                             CLUS1-CLUS2 influences the proportion of
C                             points in each half.  As for GEODIS,
C                             CLUS = 0. gives true geometric spacing;
C                                  > 0. amplifies the bunching;
C                                  < 0. reduces the bunching;
C                                  +1.0 is extreme; << -1.0 can fail;
C                                  rule of thumb: [-1.0 : +0.5]
C   X         N     R     O   Output distribution, with X (1) = XA and
C                             X (N) = XB, X (2) - X (1) = D1, and
C                             X (M) = (XA + XB)/2 for some M.
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
C   03/17/88    DAS     Initial adaptation of GEODIS for 2-sided case.
C                       Supersedes original composite use of CLUSGEO by RGL.
C
C  AUTHOR:  David Saunders/Ron Langhi, Sterling Software, Palo Alto, CA
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments.

      INTEGER
     >   IER, LUNERR, N
      REAL
     >   D1, CLUS1, CLUS2, X (N), XA, XB

C     Local constants.

      REAL
     >   FACTOR, HALF
      PARAMETER
     >   (FACTOR = 10.E+0,  ! Reasonable fraction of N applied to C1-C2.
     >    HALF   = 0.5E+0)

C     Local variables.

      INTEGER
     >   M, N1

C     Execution.

      N1 = (N + 1) / 2 + NINT ((CLUS1 - CLUS2) * FLOAT (N) / FACTOR)

      CALL GEODIS (XA, (XA + XB) * HALF, N1, D1, CLUS1, X, LUNERR, IER)
      IF (IER .NE. 0) GO TO 99

      M = N1 - 1                       ! The two distributions normally overlap

      IF (CLUS1 .EQ. CLUS2) THEN
         IF (MOD (N, 2) .NE. 0) THEN   ! Special case of desirable symmetry:
            M = N1                     ! last interval is duplicated in middle
         END IF
      END IF

      CALL GEODIS (X (M), XB, N-M+1, X (N1) - X (N1-1), CLUS2, X (M),
     >   LUNERR, IER)

   99 RETURN
      END

