C+----------------------------------------------------------------------
C
      SUBROUTINE GETXFORM (A, B, P, Q, SCALE, SHIFT)
C
C PURPOSE:
C     GETXFORM calculates the coefficients of the linear transformation
C
C                        X' <-- SCALE * X + SHIFT
C
C     such that X in interval [A, B] is transformed to X' in [P, Q].
C
C     USESCALE can be used to apply (or reverse) the transformation.
C
C     GETXFORM is a variation of GETSCALE for just the 1-D case.
C     It was introduced when the design of XFORMX was found to be
C     improper for the case of applying the same transformation to
C     more than one dataset.
C
C ARGUMENTS:
C   ARG   DIM  TYPE  I/O/S  DESCRIPTION 
C   A      -    R    I      Bounds of interval [A, B] to be transformed
C   B      -    R    I
C   P      -    R    I      Bounds of desired interval [P, Q]
C   Q      -    R    I
C   SCALE  -    R      O    Coefficients of the transformation
C   SHIFT  -    R      O
C
C HISTORY: 01/17/91  DAS  Ideas from XFORMX and GETSCALE combined for
C                         greater flexibility than XFORMX provided.
C
C AUTHOR:  David Saunders, NASA Ames/Sterling Software, Palo Alto, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT  NONE

C     Arguments:

      REAL      A, B, P, Q, SCALE, SHIFT

      SCALE = (Q - P) / (B - A)
      SHIFT = (B * P - A * Q) / (B - A)

      RETURN
      END
