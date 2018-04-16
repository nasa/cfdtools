C+------------------------------------------------------------------------------
C
      SUBROUTINE INTERVAL2 (NX, X, XFIND, LEFT)
C
C     One-liner: Ascending-data variant of INTERVAL (no ARROW argument)
C
C     Description and usage, adapted from INTERVAL:
C
C        INTERVAL2 locates the interval in a 1-D array X(:) containing
C        a target value, XFIND, using an "interpolation" search method.
C        X(:) is assumed to be ascending, so handling of descending data
C        as in INTERVAL is dispensed with here for a tiny gain in efficiency
C        where it might matter.  The normal case returns the "left-hand"
C        endpoint of the interval bracketing the target, but for the out-of-
C        range cases, the interval to be used is the first or last.
C
C        Diagrammatically, LEFT is returned as shown below for the normal
C        case (no extrapolation):
C
C           X (1)  ...   X (LEFT)   X (LEFT+1)   ...      X (NX)
C                                 ^
C                               XFIND
C
C        And for extrapolation:
C
C                      X (LEFT = 1)  ...   X (NX)
C              ^
C            XFIND
C
C        or,
C                  X (1)  ...   X (LEFT = NX-1)    X (NX)
C                                                            ^
C                                                          XFIND
C
C        Even if XFIND = X(NX), LEFT is returned as NX - 1, not NX.
C
C     Arguments:
C
C        Name  Dimension  Type  I/O/S  Description
C        NX                I    I      Number of points in array X; must
C                                      be >= 2 (no check performed).
C
C        X        NX       R    I      Array of points defining the set
C                                      of intervals to be examined. Only
C                                      the first NX-1 points are required.
C
C        XFIND             R    I      The point for which a bracketing
C                                      interval is sought.
C
C        LEFT              I    I/O    Input: guessed index of left-hand
C                                      endpoint of the interval containing
C                                      the specified point.
C
C                                      Output: index of the largest array
C                                      value <= specified point.
C                                      Special cases for data out of range:
C                                      return left endpoint of closest interval.
C                                      Thus, LEFT = 1 for XFIND < X (2), and
C                                      LEFT = NX-1 for XFIND >= X (NX-1).
C     Bibliography:
C
C        (1) Sedgewick, R.  Algorithms.  Reading: Addison-Wesley, 1983.
C            (Chap. 14)
C
C     History:
C
C        20 Oct. 1987    RAK    Original INTERVAL interpolation search method
C                               adapted from CSEVAL (RAK) and TABLE1 (DAS).
C        23 Dec. 2014    DAS    INTERVAL2 variant for ascending data only,
C                               prompted by possible heavy use in a flow solver.
C
C     Author:  Robert Kennelly and David Saunders, Sterling Software/NASA ARC.
C              (Saunders is now with ERC, Inc. at NASA Ames Research Center.)
C
C-------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN)    :: NX      ! Number of data points being searched
      REAL,    INTENT (IN)    :: X (NX)  ! Ascending/distinct abscissas
      REAL,    INTENT (IN)    :: XFIND   ! Target abscissa
      INTEGER, INTENT (INOUT) :: LEFT    ! Input with some starting guess;
                                         ! output with the left index of the
                                         ! interval containing the target
C     Local variables:

      INTEGER :: LENGTH, NXLESS1, RIGHT, TRIAL

C     Execution:

C     Simplify things by disposing of two important special cases so that
C     X (LEFT) and X (RIGHT) can really bracket XFIND. As a by-product,
C     this also takes care of the NX = 2, 3 cases.

      NXLESS1 = NX - 1

      IF (XFIND >= X (NXLESS1)) THEN
         LEFT = NXLESS1
         GO TO 990
      ELSE IF (XFIND < X (2)) THEN
         LEFT = 1
         GO TO 990
      END IF

C      ---------------------------------
C     |                                 |
C     |   X (2) <= XFIND < X (NX - 1)   |
C     |                                 |
C     |   NX > 3                        |
C     |                                 |
C      ---------------------------------

C     Adjust the pointers. We hope that the calling routine has provided
C     a reasonable guess (since it's probably working on an ordered array
C     of points to evaluate), but check anyway.

      LEFT = MIN (MAX (2, LEFT), NX - 2)

      IF (XFIND >= X (LEFT)) THEN
         IF (XFIND < X (LEFT + 1)) THEN

C           XFIND is in the original guessed-at interval.

            GO TO 990
         ELSE

C           We'll look farther to the right. Evidently LEFT was < NX - 2.

            RIGHT = NXLESS1
            LEFT  = LEFT + 1
         END IF
      ELSE

C        Look to the left of the guess. Evidently LEFT was > 2.

         RIGHT = LEFT
         LEFT  = 2
      END IF

C      ----------------------------------
C     |                                  |
C     |   2 <= LEFT < RIGHT <= NX - 1    |
C     |                                  |
C      ----------------------------------

C     The interval length must decrease each time through - terminate
C     when the correct interval is found or when the interval length
C     cannot be decreased.

   10 CONTINUE
         LENGTH = RIGHT - LEFT
         IF (LENGTH > 1) THEN

C           The trial value is a "linear" estimate of the left-hand endpoint
C           of the interval bracketing the target XFIND, with protection
C           against round-off (which can affect convergence).

            TRIAL = MIN (RIGHT - 1, LEFT + MAX (0, INT (REAL (LENGTH) *
     &         (XFIND - X (LEFT)) / (X (RIGHT) - X (LEFT)))))

C            ------------------------------------------
C           |                                          |
C           |   2 <= LEFT <= TRIAL < RIGHT <= NX - 1   |
C           |                                          |
C            ------------------------------------------

C           Adjust pointers. Increase LEFT or decrease RIGHT until done.

            IF (XFIND >= X (TRIAL + 1)) THEN
               LEFT  = TRIAL + 1
            ELSE IF (XFIND < X (TRIAL)) THEN
               RIGHT = TRIAL
            ELSE

C              We're done: XFIND is in the interval [X (TRIAL), X (TRIAL+1)).

               LEFT  = TRIAL
               GO TO 990
            END IF
            GO TO 10

         END IF

C     Termination.

  990 RETURN

      END SUBROUTINE INTERVAL2
