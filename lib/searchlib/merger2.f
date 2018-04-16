C-----------------------------------------------------------------------
C
      SUBROUTINE MERGER2 (NINSERT, XINSERT, YINSERT, NSORTED, XINOUT,
     >                    YINOUT, TOL)
C
C ACRONYM:  MERGE two Real (X, Y) datasets (one of which is sorted in X)
C           -----     -
C PURPOSE:
C        MERGER2 inserts the (not necessarily ordered) NINSERT real values
C     of XINSERT (*), YINSERT (*) into the ordered-in-X dataset XINOUT (*),
C     YINOUT (*) such that the output XINOUT (*) remains ordered, with no
C     duplicates.  Its length, NSORTED, is updated accordingly.
C
C ARGUMENTS:
C     ARG       DIM     TYPE I/O/S DESCRIPTION
C    NINSERT              I  I     Number of values to merge;  NINSERT >= 1
C    XINSERT, NINSERT     R  I     (X, Y) pairs to merge (need not be ordered)
C    YINSERT
C    NSORTED              I  I/O   Input & output with length of XINOUT (*)
C                                 NSORTED (I) >= 2;
C                                 NSORTED (O) <= NSORTED (I) + NINSERT
C    XINOUT   NSORTED (I) R  I/O   Input & output list, in either ascending
C           + NINSERT            or descending order, with no duplicates
C    YINOUT   NSORTED (I) R  I/O   Corresponding Y values
C           + NINSERT
C    TOL                  R  I     Tolerance used for duplicate comparisons;
C                                  on X; TOL >= 0.
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
C ERROR HANDLING:  None.  Make sure the sorted array XINOUT (*) has room.
C
C EXTERNAL REFERENCES:
C    INTERVAL   "Interpolation search" utility
C
C HISTORY:
C   02/09/89   D.Saunders  Initial implementation of MERGER (just Xs).
C   10/24/93     "    "    MERGER2 adapted to handle (X, Y) datasets.
C   01/25/11     "    "    UTCOPY is best avoided; also eliminated the
C              ERC/ARC     kludge of using an extra element in XINOUT.
C
C AUTHOR:  David Saunders, Sterling Software/NASA Ames Research Center, CA
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER
     >   NINSERT, NSORTED
      REAL
     >   XINSERT (NINSERT), YINSERT (NINSERT),
     >   XINOUT (NSORTED + NINSERT), YINOUT (NSORTED + NINSERT), TOL

C     Local variables:

      INTEGER
     >   I, LEFT, NSHIFT
      REAL
     >   ARROW, XIN

C     Procedures:

      EXTERNAL
     >   INTERVAL

C     Execution.

      ARROW = SIGN (1.E+0, XINOUT (2) - XINOUT (1))

      DO I = 1, NINSERT

         XIN = XINSERT (I)

C        Locate the interval for insertion:

         CALL INTERVAL (NSORTED, XINOUT, XIN, ARROW, LEFT)

C        LEFT cannot be NSORTED the way INTERVAL works, so:

         IF (XIN * ARROW >= XINOUT (NSORTED) * ARROW) LEFT = LEFT + 1

         IF (ABS (XIN - XINOUT (LEFT)) > TOL) THEN
            NSHIFT = NSORTED - LEFT
            IF (NSHIFT > 0) THEN

C              Make room by shifting elements to the "right":

               XINOUT (LEFT + 2 : LEFT + 1 + NSHIFT) =
     >            XINOUT (LEFT + 1 : LEFT + NSHIFT)
               YINOUT (LEFT + 2 : LEFT + 1 + NSHIFT) =
     >            YINOUT (LEFT + 1 : LEFT + NSHIFT)

CCC            CALL UTCOPY (-NSHIFT, XINOUT (LEFT + 1),
CCC  >            XINOUT (LEFT + 2))
            END IF

            XINOUT (LEFT + 1) = XIN
            YINOUT (LEFT + 1) = YINSERT (I)
            NSORTED = NSORTED + 1
         END IF

      END DO

      END SUBROUTINE MERGER2
