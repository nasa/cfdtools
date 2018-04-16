C+----------------------------------------------------------------------
C
      SUBROUTINE SHIFTX ( NX, X, SHIFT )
C
C PURPOSE: SHIFTX offsets the elements of the given array by the given
C          amount, in place:
C
C                    X := X + SHIFT
C
C PARAMETERS:
C   ARG     DIM    TYPE  I/O/S DESCRIPTION 
C   NX       -       I     I   Number of elements to offset
C   X        *       R    I/O  Array to be shifted
C   SHIFT    -       R     I   Offset added to each element
C
C AUTHOR:  David Saunders, Sterling Software, Palo Alto, CA., June 1986
C
C-----------------------------------------------------------------------

      INTEGER   NX
      REAL      X(*), SHIFT

      INTEGER   I
      REAL      OFFSET

      OFFSET = SHIFT
      DO 20 I = 1, NX
         X(I) = X(I) + OFFSET
   20 CONTINUE

      RETURN
      END
