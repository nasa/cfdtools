C+------------------------------------------------------------------------------
C
      SUBROUTINE INTERVAL_I (NLIST, LIST, LFIND, ARROW, LEFT)
C
C     One-liner: Interpolation search for interval containing an integer.
C     ----------
C
C     Description and usage:
C     ----------------------
C
C        This is the integer-list analog of subroutine INTERVAL, q.v.
C     It performs an efficient search for the interval in an ordered
C     LIST of integers that contains the target integer, LFIND.  The
C     LIST may be ascending or descending according to ARROW = +/-1.
C
C        The normal case returns the "left-hand" endpoint of the interval
C     bracketing the point, but for the out-of-range cases below or above
C     the range of the list, the interval to be used is the first or last.
C     Diagrammatically, LEFT is returned as shown below for the normal case
C     (no extrapolation):
C
C        LIST (1)  ...   LIST (LEFT)   LIST (LEFT+1)   ...      LIST (NLIST)
C                                    ^
C                                  LFIND
C
C     And for extrapolation:
C
C               LIST (LEFT = 1)  ...   LIST (NLIST)
C          ^
C        LFIND
C
C     or,
C               LIST (1)  ...   LIST (LEFT = NLIST-1)    LIST (NLIST)
C                                                                      ^
C                                                                    LFIND
C
C     If the target matches one of the LIST integers, the index of that
C     list element is returned as LEFT.  Thus, the condition for a bracket
C     of an interior point is:
C
C        LIST (LEFT) <= LFIND < LIST (LEFT+1)  if  ARROW = +1,  or
C        LIST (LEFT) >= LFIND > LIST (LEFT+1)  if  ARROW = -1.
C
C        This is a low-level routine with minimal error checking. The
C     calling program is assumed to have verified the following:
C
C     (1)  NLIST >= 2
C     (2)  LIST is strictly monotonic
C     (3)  ARROW = +1 or -1
C
C        The interpolation search method used was adapted from ideas in
C     Sedgewick's book referenced below.
C
C     Arguments:
C     ----------
C
C     Name  Dimension  Type  I/O/S  Description
C     NLIST             I    I      Number of points in array LIST; must
C                                   be >= 2 (no check performed).
C
C     LIST    NLIST     I    I      Array of integers defining the set
C                                   of intervals to be examined.
C
C     LFIND             I    I      The targetinteger for which a bracketing
C                                   interval is sought.
C
C     ARROW             I    I      Monotonicity indicator for input
C                                   array LIST:
C                                     -1  strictly decreasing
C                                      0  NOT ALLOWED!
C                                     +1  strictly increasing
C                                   Supplied by the calling routine for
C                                   reasons of speed (not checked).
C
C     LEFT              I    I/O    Input: guessed index of left-hand
C                                   endpoint of the interval containing
C                                   the target integer.
C
C                                   Output: index of the largest array
C                                   value <= specified point (if ARROW = +1).
C                                   Special case for data out of range:
C                                   return left endpoint of closest interval.
C                                   Thus, LEFT = 1 for LFIND < LIST (2), and
C                                   LEFT = NLIST-1 for LFIND >= LIST (NLIST-1).
C                                   (If ARROW = -1, reverse the inequalities.)
C
C     Notes:
C     ------
C
C     (1)  In speed-critical applications, it might be a good idea to build
C          this algorithm in-line since it is typically called many times
C          from within a loop. Another potential speed-up is removal of the
C          ARROW multiplies, which restricts the method to increasing data.
C          So far, the simplicity of separating out the messy search details
C          and the generality of bi-directional searching have outweighed
C          the modest speed penalty incurred.
C
C     Bibliography:
C     -------------
C
C     (1) Sedgewick, R.  Algorithms.  Reading: Addison-Wesley, 1983.
C            (Chap. 14)
C
C     Author:  Robert Kennelly and David Saunders, Sterling Federal Systems
C     -------
C
C     History:
C     --------
C
C     20 Oct. 1987    RAK     Interpolation search adapted (with mods.
C                             for bidirectional search and some minor
C                             repair) from CSEVAL (RAK) and TABLE1 (DAS).
C     08 Aug. 1988 D.Saunders Clarified descriptions of bracketing, where
C                             the inequalities depend upon ARROW.
C     02 June 2004     "      Integer version of INTERVAL.
C
C-------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments.

      INTEGER
     &   NLIST, LIST (NLIST), LFIND, ARROW, LEFT

C     Local variables.

      INTEGER
     &   LENGTH, LBYARROW, NXLESS1, RIGHT, TRIAL

C     Execution.
C     ----------

      LBYARROW = LFIND * ARROW

C     Simplify things by disposing of two important special cases so that
C     LIST (LEFT) and LIST (RIGHT) can really bracket LFIND. As a by-product,
C     this also takes care of the NLIST = 2, 3 cases.

      NXLESS1 = NLIST - 1

      IF (LBYARROW .GE. LIST (NXLESS1) * ARROW) THEN
         LEFT = NXLESS1
         GO TO 990
      ELSE IF (LBYARROW .LT. LIST (2) * ARROW) THEN
         LEFT = 1
         GO TO 990
      END IF

C      ------------------------------------------
C     |                                          |
C     |   LIST (2) <= LFIND < LIST (NLIST - 1)   |
C     |            - or -                        |
C     |   LIST (2) > LFIND >= LIST (NLIST - 1)   |
C     |                                          |
C     |   NLIST > 3                              |
C     |                                          |
C      ------------------------------------------

C     Adjust the pointers. We hope that the calling routine has provided
C     a reasonable guess (since it's probably working on an ordered array
C     of points to evaluate), but check anyway.

      LEFT = MIN (MAX (2, LEFT), NLIST - 2)

      IF (LBYARROW .GE. LIST (LEFT) * ARROW) THEN
         IF (LBYARROW .LT. LIST (LEFT + 1) * ARROW) THEN

C           LFIND is in the original guessed-at interval.

            GO TO 990
         ELSE

C           We'll look farther to the right. Evidently LEFT was < NLIST - 2.

            RIGHT = NXLESS1
            LEFT  = LEFT + 1
         END IF
      ELSE

C        Look to the left of the guess. Evidently LEFT was > 2.

         RIGHT = LEFT
         LEFT  = 2
      END IF

C      ------------------------------------
C     |                                    |
C     |   2 <= LEFT < RIGHT <= NLIST - 1   |
C     |                                    |
C      ------------------------------------

C     The interval length must decrease each time through - terminate
C     when the correct interval is found or when the interval length
C     cannot be decreased.

   10 CONTINUE
         LENGTH = RIGHT - LEFT
         IF (LENGTH .GT. 1) THEN

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

            IF (LBYARROW .GE. LIST (TRIAL + 1) * ARROW) THEN
               LEFT  = TRIAL + 1
            ELSE IF (LBYARROW .LT. LIST (TRIAL) * ARROW) THEN
               RIGHT = TRIAL
            ELSE

C              We're done: LFIND is in the interval [LIST (TRIAL), LIST (TRIAL+1)).

               LEFT  = TRIAL
               GO TO 990
            END IF
            GO TO 10

         END IF

C     Termination.
C     ------------

  990 CONTINUE
      RETURN
      END
