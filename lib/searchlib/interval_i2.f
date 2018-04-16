C+------------------------------------------------------------------------------
C
      SUBROUTINE INTERVAL_I2 (NLIST, LIST, LFIND, LEFT)
C
C     1-liner: Ascending-integer-data variant of INTERVAL_I (no ARROW argument)
C
C     Description and usage, adapted from INTERVAL/INTERVAL2/INTERVAL_I:
C
C        INTERVAL_I2 locates the interval in a 1-D array LIST(:) containing
C        a target value, LFIND, using an "interpolation" search method.
C        LIST(:) is assumed to be ascending, so handling of descending data
C        as in INTERVAL[_I] is dispensed with here for a tiny gain in efficiency
C        where it might matter.  The normal case returns the "left-hand"
C        endpoint of the interval bracketing the target, but for the out-of-
C        range cases, the interval to be used is the first or last.
C
C        Diagrammatically, LEFT is returned as shown below for the normal
C        case (no extrapolation):
C
C           LIST (1)  ...   LIST (LEFT)   LIST (LEFT+1)   ...      LIST (NLIST)
C                                       ^
C                                     LFIND
C
C        And for extrapolation:
C
C                      LIST (LEFT = 1)  ...   LIST (NLIST)
C              ^
C            LFIND
C
C        or,
C            LIST (1)  ...   LIST (LEFT = NLIST-1)    LIST (NLIST)
C                                                                        ^
C                                                                      LFIND
C
C        Even if LFIND = LIST(NLIST), LEFT is returned as NLIST - 1, not NLIST.
C
C     Arguments:
C
C        Name  Dimension  Type  I/O/S  Description
C        NLIST                I    I   Number of points in array LIST; must
C                                      be >= 2 (no check performed).
C
C        LIST        NLIST    R    I   Array of points defining the set
C                                      of intervals to be examined. Only
C                                      the first NLIST-1 points are required.
C
C        LFIND             R    I      The target integer for which a bracketing
C                                      interval is sought.
C
C        LEFT              I    I/O    Input: guessed index of left-hand
C                                      endpoint of the interval containing
C                                      the target integer.
C
C                                      Output: index of the largest array
C                                      value <= specified target.
C                                      Special cases for data out of range:
C                                      return left endpoint of closest interval.
C                                      Thus, LEFT = 1 if LFIND < LIST (2), and
C                                      LEFT = NLIST-1 if LFIND >= LIST(NLIST-1).
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
C        24 Dec. 2014    DAS    INTERVAL_I2 variant of INTERVAL2/INTERVAL_I,
C                               for integer data/ascending only.
C
C     Author:  Robert Kennelly and David Saunders, Sterling Software/NASA ARC.
C              (Saunders is now with ERC, Inc. at NASA Ames Research Center.)
C
C-------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN)    :: NLIST   ! Number of data points being searched
      INTEGER, INTENT (IN)    :: LIST (NLIST)  ! Ascending/distinct integers
      INTEGER, INTENT (IN)    :: LFIND   ! Target integer
      INTEGER, INTENT (INOUT) :: LEFT    ! Input with some starting guess;
                                         ! output with the left index of the
                                         ! interval containing the target
C     Local variables:

      INTEGER :: LENGTH, NLISTM1, RIGHT, TRIAL

C     Execution:

C     Simplify things by disposing of two important special cases so that
C     LIST (LEFT) and LIST (RIGHT) can really bracket LFIND. As a by-product,
C     this also takes care of the NLIST = 2, 3 cases.

      NLISTM1 = NLIST - 1

      IF (LFIND >= LIST (NLISTM1)) THEN
         LEFT = NLISTM1
         GO TO 990
      ELSE IF (LFIND < LIST (2)) THEN
         LEFT = 1
         GO TO 990
      END IF

C      ------------------------------------------
C     |                                          |
C     |   LIST (2) <= LFIND < LIST (NLIST - 1)   |
C     |                                          |
C     |   NLIST > 3                              |
C     |                                          |
C      ------------------------------------------

C     Adjust the pointers. We hope that the calling routine has provided
C     a reasonable guess (since it's probably working on an ordered array
C     of points to evaluate), but check anyway.

      LEFT = MIN (MAX (2, LEFT), NLIST - 2)

      IF (LFIND >= LIST (LEFT)) THEN
         IF (LFIND < LIST (LEFT + 1)) THEN

C           LFIND is in the original guessed-at interval.

            GO TO 990
         ELSE

C           We'll look farther to the right. Evidently LEFT was < NLIST - 2.

            RIGHT = NLISTM1
            LEFT  = LEFT + 1
         END IF
      ELSE

C        Look to the left of the guess. Evidently LEFT was > 2.

         RIGHT = LEFT
         LEFT  = 2
      END IF

C      -------------------------------------
C     |                                     |
C     |   2 <= LEFT < RIGHT <= NLIST - 1    |
C     |                                     |
C      -------------------------------------

C     The interval length must decrease each time through - terminate
C     when the correct interval is found or when the interval length
C     cannot be decreased.

   10 CONTINUE
         LENGTH = RIGHT - LEFT
         IF (LENGTH > 1) THEN

C           The trial value is a "linear" estimate of the left-hand endpoint
C           of the interval bracketing the target LFIND, with protection
C           against round-off (which can affect convergence).

            TRIAL = MIN (RIGHT - 1,
     &         LEFT + MAX (0, INT (REAL (LENGTH) *
     &                     REAL (LFIND - LIST (LEFT)) /
     &                     REAL (LIST (RIGHT) - LIST (LEFT)))))

C            ---------------------------------------------
C           |                                             |
C           |   2 <= LEFT <= TRIAL < RIGHT <= NLIST - 1   |
C           |                                             |
C            ---------------------------------------------

C           Adjust pointers. Increase LEFT or decrease RIGHT until done.

            IF (LFIND >= LIST (TRIAL + 1)) THEN
               LEFT  = TRIAL + 1
            ELSE IF (LFIND < LIST (TRIAL)) THEN
               RIGHT = TRIAL
            ELSE

C              We're done: LFIND is in [LIST (TRIAL), LIST (TRIAL+1)).

               LEFT  = TRIAL
               GO TO 990
            END IF
            GO TO 10

         END IF

C     Termination.

  990 RETURN

      END SUBROUTINE INTERVAL_I2
