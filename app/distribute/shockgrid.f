C*******************************************************************************
C
      SUBROUTINE SHOCKGRID (N, D1, DE, NMARGIN, SEDGE, S, LUNERR, IER)
C
C        SHOCKGRID generates a specialized 1-D grid intended for radial lines
C     starting in a boundary layer and ending just outside a shock envelope.
C     The first and last S values are assumed to be in S(1) and S(N).  NMARGIN
C     points are distributed beyond SEDGE (the given shock edge) with 1-sided
C     stretching.  More precisely, counting 2 points at SEDGE - DE and SEDGE
C     that overlap the 2-sided distribution inboard, the outer distribution has
C     NMARGIN + 2 points on the interval [SEDGE - DE, S(N)].
C
C     DE is probably derived from D1 at the wall.  For example,  DE = 250 D1.
C     Spacing is 2-sided Vinokur-type from S(1) to S(N-NMARGIN) using D1 and
C     DE with an overlap of any margin points by one cell.
C
C     NMARGIN = 0 can be used to indicate that no shock edge is present.
C     DE will then be the outermost increment, with no 1-sided stretching beyond
C     SEDGE = S(N).  However, it is not obvious how to blend such a degenerate
C     case (likely to be in a wake grid block) with adjacent radial lines that
C     do encounter the shock, so this option is probably not helpful.
C
C     See HTDIS4 and EXPDIS4 for further details, including the LUNERR and IER
C     arguments.
C
C     08/22/98  DAS  Last version of BLGRID from which SHOCKGRID was adapted.
C     07/01/06   "   Initial implementation for heat-shield applications.
C
C*******************************************************************************

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN)    :: N, NMARGIN, LUNERR
      REAL,    INTENT (IN)    :: D1, DE, SEDGE     ! See above
      REAL,    INTENT (INOUT) :: S(N)
      INTEGER, INTENT (OUT)   :: IER

C     Local constants:

      REAL, PARAMETER :: ONE = 1., ZERO = 0.

C     Local variables:

      INTEGER  I1, N1SIDED, N2SIDED
      REAL     BETA, SI1, SI2

C     Execution:

      SI1        = S(1)
      SI2        = SEDGE
      N2SIDED    = N - NMARGIN  ! # points inside the shock
      S(N2SIDED) = SEDGE

C     2-sided Vinokur distribution:

C *** CALL VINOKUR (1, N2SIDED, D1, DE, S, LUNERR, IER) ! In-lined below

      CALL HTDIS4 (.TRUE., SI1, SI2, D1, DE, N2SIDED, S, -LUNERR, IER)

C     Overlap two distributions by one cell?

      IF (NMARGIN > 0) THEN

         SI1 = SEDGE - DE
         SI2 = S(N)
         I1  = N - NMARGIN - 1
         N1SIDED = NMARGIN + 2

         CALL EXPDIS4 (3, SI1, SI2, DE, N1SIDED, S(I1), BETA, -LUNERR)

      END IF

      END SUBROUTINE SHOCKGRID
