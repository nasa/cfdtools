C*******************************************************************************
C
      SUBROUTINE BLGRID (N, D1, D2, NBLAYER, RBLAYER, X, LUNERR, IER)
C
C        BLGRID generates a specialized 1-D grid intended for radial lines
C     starting in a boundary layer and ending in the far field.  The first
C     and last X values are assumed to be in X(1) and X(N).  Points 1
C     through NBLAYER have spacing that varies geometrically according to
C     the ratio RBLAYER (or uniformly if RBLAYER is 1.0).  Spacing is
C     Vinokur-type beyond that.  X is probably arc-length (zero to STOTAL),
C     but could be decreasing if desired.  D1 & D2 are both positive.
C
C     Use NBLAYER <= 2 to suppress the special treatment.
C     See HTDIS4 for further details, including the LUNERR and IER arguments.
C
C     07/08/97  DAS  Effort to emulate the single-piece C grid result for
C                    the two-halves 2-D case forward of the leading edge.
C     08/06/97   "   Generalized to allow application off the trailing edge.
C     08/22/98   "   In-lined the VINOKUR call because it was easy.
C
C*******************************************************************************

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN)    :: N, NBLAYER, LUNERR
      REAL,    INTENT (IN)    :: D1, D2, RBLAYER
      REAL,    INTENT (INOUT) :: X(N)
      INTEGER, INTENT (OUT)   :: IER

C     Local constants:

      REAL, PARAMETER :: ONE = 1., ZERO = 0.

C     Local variables:

      INTEGER  I, I1, NB
      REAL     DX, RANGE, XI1, XI2

C     Execution:

      DX = SIGN (D1, X(N) - X(1))
      NB = MAX (NBLAYER, 2) ! Force a single pass if special treatment
                            ! is suppressed
      DO I = 2, NB
         X(I) = X(I-1) + DX
         DX   = DX * RBLAYER
      END DO

      DX    = ABS (DX / RBLAYER)
      I1    = NB - 1
      XI1   = X(I1)
      XI2   = X(N)
      RANGE = XI2 - XI1

C     Overlap the two distributions by one cell:

C *** CALL VINOKUR (I1, N, DX, D2, X, LUNERR, IER) ! In-lined below

C     Normalized Vinokur distribution in ascending order:

      CALL HTDIS4 (.TRUE., ZERO, ONE, ABS (DX/RANGE), ABS (D2/RANGE),
     >             N - I1 + 1, X(I1), -LUNERR, IER)

C     Denormalize:

      DO I = I1, N - 1
         X(I) = XI1 + X(I) * RANGE
      END DO

      X(N) = XI2 ! Exactly

      END SUBROUTINE BLGRID
