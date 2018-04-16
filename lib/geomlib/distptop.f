C+----------------------------------------------------------------------
C
      SUBROUTINE DISTPTOP (P0, P1, DIST)
C
C     One-liner:  Computes the distance between two points in space.
C     ----------
C
C     Description and usage:
C     ----------------------
C
C     DISTPTOP computes the Euclidean distance between two points, paying
C     attention to eliminating unnecessary overflow by prescaling the
C     terms to be squared beforehand.
C
C     Arguments:
C     ----------
C
C     Name     Dimension  Type  I/O/S  Description
C     P0       3           R    I      Field point.
C
C     P1       3           R    I      Reference point. (Of course, it
C                                      doesn't matter here which point
C                                      is which; the names were chosen
C                                      to match the other DIST routines.)
C
C     DIST                 R      O    (Scalar) distance between P0 and P1.
C
C     Environment:  Digital VAX-11/780 VMS FORTRAN 77 v4.7
C     ------------
C
C     Notes:
C     ------
C
C     (1)  IMPLICIT NONE and eight character symbols are not (yet) standard.
C
C     Author:  Robert Kennelly, Sterling Federal Systems/NASA-Ames
C     -------
C
C     Development history:
C     --------------------
C
C     18 Jan. 1989  RAK/DAS  Initial design and coding, based in part
C                            on 2-dim. function CHORD.
C
C-----------------------------------------------------------------------

C     Declarations.
C     -------------

      IMPLICIT NONE

C     Constants.

      REAL
     &   ZERO, ONE
      PARAMETER
     &  (ZERO = 0.0E+0,
     &   ONE  = 1.0E+0)

C     Arguments.

      REAL
     &   DIST, P0 (3), P1 (3)

C     Local variables.

      REAL
     &   DENOM, DX, DY, DZ, INVDENOM

C     Execution.
C     ----------

      DX    = ABS (P1 (1) - P0 (1))
      DY    = ABS (P1 (2) - P0 (2))
      DZ    = ABS (P1 (3) - P0 (3))
      DENOM = MAX (DX, DY, DZ)

      IF (DENOM .GT. ZERO) THEN

C        Divide the terms to be squared by the largest delta to reduce
C        the chance of overflow. (An "extra" divide below saves a lot
C        of program logic.)

         INVDENOM = ONE / DENOM
         DIST     = DENOM * SQRT ((DX * INVDENOM) ** 2 +
     &                            (DY * INVDENOM) ** 2 +
     &                            (DZ * INVDENOM) ** 2)
      ELSE

C        The points are identical.

         DIST = ZERO
      END IF

C     Termination.
C     ------------

      RETURN
      END

