C+---------------------------------------------------------------------
C
      SUBROUTINE RVERSE ( NX, X1, X2 )
C
C  PURPOSE:  RVERSE reverses the order of the elements of X1(*), re-
C            turning the result in X2(*). In-place reordering is OK.
C
C  ARGUMENTS:
C  ARG  DIM   TYPE I/O/S DESCRIPTION
C  NX   -       I    I   No. of elements in given array (odd or even).
C  X1   NX      R    I   Array whose ordering is to be reversed.
C  X2   NX      R    O   Reordered version of X1(*).  The calling pro-
C                        gram may pass the same array for X1 & X2.
C
C  ENVIRONMENT:  FORTRAN IV
C
C  AUTHOR: David Saunders, Informatics, Palo Alto, CA.  (07/30/82)
C
C----------------------------------------------------------------------
C
      DIMENSION  X1(NX), X2(NX)
C
      NXBY2 = (NX+1)/2
C
      DO 20 I = 1, NXBY2
         XI     = X1(I)
         II     = NX + 1 - I
         X2(I)  = X1(II)
         X2(II) = XI
 20   CONTINUE
C
      RETURN
      END
