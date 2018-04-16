C-----------------------------------------------------------------------
C
      SUBROUTINE MERGER (NINSERT, RLISTIN, NSORTED, RLISTIO, TOL)
C
C ACRONYM:  MERGE two Real lists (one of which is sorted)
C           -----     -
C PURPOSE:
C        MERGER inserts the (not necessarily ordered) NINSERT real values
C     of RLISTIN (*) into the ordered list RLISTIO (*) such that the
C     output RLISTIO (*) remains ordered, with no duplicates.  Its length,
C     NSORTED, is updated accordingly.
C
C ARGUMENTS:
C     ARG       DIM    TYPE I/O/S DESCRIPTION
C    NINSERT             I  I     Number of values to merge;  NINSERT >= 1
C    RLISTIN  NINSERT    R  I     Values to merge (need not be ordered)
C    NSORTED             I  I/O   Input & output with length of RLISTIO (*)
C                                 NSORTED (I) >= 2;
C                                 NSORTED (O) <= NSORTED (I) + NINSERT
C    RLISTIO NSORTED (I) R  I/O   Input & output list, in either ascending
C            + NINSERT            or descending order, with no duplicates
C    TOL                 R  I     Tolerance used for duplicate comparisons;
C                                 TOL >= 0.
C METHOD:
C        Search utility INTERVAL finds each point of insertion efficiently;
C     it also handles ordering in either direction.  The utility UTCOPY
C     was originally used to make room for insertions by shifting values to
C     the "right" but Fortran 90 now provides for correct in-place shifts.
C     Avoiding duplicates is problematic; hence the TOLerance argument in
C     case exact equality is inappropriate.
C
C USAGE:
C        Simple overwriting of the nearest point may seem adequate for the
C     initial application of ensuring that data point abscissas are included
C     exactly when a curve fit is plotted.  But this can fail if data points
C     are very close together.  True insertion (as done here) should find
C     non-graphical applications too.
C
C ERROR HANDLING:  None.  Make sure the sorted array RLISTIO (*) has room.
C
C EXTERNAL REFERENCES:
C    INTERVAL   "Interpolation search" utility
C
C HISTORY:
C   02/09/89   D.Saunders  Initial implementation.
C   01/25/11     "    "    UTCOPY is best avoided; also eliminated the
C              ERC/ARC     kludge of using an extra element in RLISTIO.
C
C AUTHOR:  David Saunders, Sterling Software/NASA Ames Research Center, CA
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER
     >   NINSERT, NSORTED
      REAL
     >   RLISTIN (NINSERT), RLISTIO (NSORTED + NINSERT), TOL

C     Local variables:

      INTEGER
     >   I, LEFT, NSHIFT
      REAL
     >   ARROW, RIN

C     Procedures:

      EXTERNAL
     >   INTERVAL

C     Execution.

      ARROW = SIGN (1.E+0, RLISTIO (2) - RLISTIO (1))

      DO I = 1, NINSERT

         RIN = RLISTIN (I)

C        Locate the interval for insertion:

         CALL INTERVAL (NSORTED, RLISTIO, RIN, ARROW, LEFT)

C        LEFT cannot be NSORTED the way INTERVAL works, so:

         IF (RIN * ARROW >= RLISTIO (NSORTED) * ARROW) LEFT = LEFT + 1

         IF (ABS (RIN - RLISTIO (LEFT)) > TOL) THEN
            NSHIFT = NSORTED - LEFT
            IF (NSHIFT > 0) THEN

C              Make room by shifting elements to the "right":

               RLISTIO (LEFT + 2 : LEFT + 1 + NSHIFT) =
     >            RLISTIO (LEFT + 1 : LEFT + NSHIFT)

CCC            CALL UTCOPY (-NSHIFT, RLISTIO (LEFT + 1),
CCC  >            RLISTIO (LEFT + 2))
            END IF

            RLISTIO (LEFT + 1) = RIN
            NSORTED = NSORTED + 1
         END IF

      END DO

      END SUBROUTINE MERGER
