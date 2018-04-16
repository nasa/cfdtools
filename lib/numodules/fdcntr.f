C+----------------------------------------------------------------------
C
      SUBROUTINE FDCNTR (I, X, F, FP, FPP)
C
C  PURPOSE: FDCNTR returns central 3-point finite-difference estimates
C           of the 1st and 2nd derivatives at the point  (X(I), F(I)).
C           Use FD1SID for the end-point cases.
C
C  ARGS:    Obvious from PURPOSE.
C
C  METHOD:  FPP is computed first, in case only FP is desired, so that
C           the same item may be passed for both arguments.
C
C  HISTORY: 08/17/91  DAS  FDCNTR adapted from FD12K's in-line code,
C                          as for FD1SID, for the case of a single I
C                          at a time (which FD12K can't do).
C
C  AUTHOR:  David Saunders, Sterling Software/NASA Ames, Palo Alto, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT   NONE

C     Arguments:

      INTEGER    I
      REAL       X (*), F (*), FP, FPP

C     Local constants:

      REAL       ONE
      PARAMETER (ONE = 1.E+0)

C     Local variables:

      REAL     DEL1, DEL2, DIV, DX1, DX2, W

C     Execution:

      DX1  =  X (I) - X (I-1)
      DEL1 = (F (I) - F (I-1)) / DX1
      DX2  =  X (I+1) - X (I)
      DEL2 = (F (I+1) - F (I)) / DX2
      DIV  = ONE / (DX1 + DX2)
      W    = DX2 * DIV
      FPP  = (DEL2 - DEL1) * (DIV + DIV)
      FP   = W * DEL1 + (ONE - W) * DEL2

      RETURN
      END
