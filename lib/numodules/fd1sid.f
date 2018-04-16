C+----------------------------------------------------------------------
C
      SUBROUTINE FD1SID (I, INC, X, F, FP, FPP)
C
C  PURPOSE: FD1SID returns one-sided 3-point finite-difference estimates
C           of the 1st and 2nd derivatives at the point  ( X(I), F(I) ).
C           If INC = 1, points I, I+1, I+2 are used,  while if INC = -1,
C           points  I-2, I-1, I are used. The abscissas need not be uni-
C           formly spaced.  
C
C  ARGS:    Obvious from PURPOSE.
C
C  METHOD:  FPP is computed first,  in case only FP is desired,  so that
C           the same item may be passed for both arguments. The formula-
C           tion is similar to that of central differencing - see FD12K.
C
C  HISTORY: 12/27/85  DAS  Initial implementation  (prompted by the need
C                          to approximate an airfoil leading edge with a
C                          cubic having specified slope at one end).
C           08/21/89  DAS  Formulation revised as for centrals (FD12K).
C
C  AUTHOR:  David Saunders, Sterling Software/NASA Ames, Palo Alto, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT   NONE

C     Arguments:

      INTEGER    INC, I
      REAL       X (*), F (*), FP, FPP

C     Local constants:

      REAL       ONE
      PARAMETER (ONE = 1.E+0)

C     Local variables:

      INTEGER  I1, I2
      REAL     DEL1, DEL2, DIV, DX1, DX2, W

C     Execution:

C     The minus signs take care of themselves for the backward case.

      I1   = I  + INC
      I2   = I1 + INC
      DX1  =  X (I1) - X (I)
      DEL1 = (F (I1) - F (I)) / DX1
      DX2  =  X (I2) - X (I1)
      DEL2 = (F (I2) - F (I1)) / DX2
      DIV  = ONE / (DX1 + DX2)
      W    = -DX1 * DIV
      FPP  = (DEL2 - DEL1) * (DIV + DIV)
      FP   = W * DEL2 + (ONE - W) * DEL1

      RETURN
      END
