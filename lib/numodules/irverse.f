C+---------------------------------------------------------------------
C
      SUBROUTINE IRVERSE (NI, I1, I2)
C
C  PURPOSE:  IRVERSE reverses the order of the elements of I1(*), re-
C            turning the result in I2(*). In-place reordering is OK.
C
C  ARGUMENTS:
C  ARG  DIM   TYPE I/O/S DESCRIPTION
C  NI   -       I    I   No. of elements in given array (odd or even).
C  I1   NI      R    I   Array whose ordering is to be reversed.
C  I2   NI      R    O   Reordered version of I1(*).  The calling pro-
C                        gram may pass the same array for I1 & I2.
C
C  AUTHOR:   David Saunders, Informatics, Palo Alto, CA.  (07/30/82)
C
C  03/16/2017  D.A.Saunders  IRVERSE adapted from the ancient RVERSE.
C
C----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN)    :: NI
      INTEGER, INTENT (INOUT) :: I1(NI), I2(NI)

C     Local variables:

      INTEGER :: I, II, NIBY2

C     Execution:

      NIBY2 = (NI+1)/2

      DO I = 1, NIBY2
         II     = I1(I)
         II     = NI + 1 - I
         I2(I)  = I1(II)
         I2(II) = II
      END DO

      END SUBROUTINE IRVERSE
