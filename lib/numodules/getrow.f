C+---------------------------------------------------------------------
C
      SUBROUTINE GETROW (IDIM, I, J1, J2, A, ROW)
C
C  PURPOSE:
C     GETROW extracts the (indicated part of the) Ith row from a 2-D array.
C
C  ARGUMENTS:
C     ARG   TYPE I/O/S   DIM     DESCRIPTION
C     IDIM   I     I      -      Declared row dimension of A
C     I      I     I      -      Target row
C     J1,    I     I      -      First and last elements of target row
C     J2                         to be extracted
C     A      R     I  (IDIM, *)  Array from which elements of the Ith row
C                                are extracted
C     ROW    R     O      -      Elements J1 : J2 of row I of A are
C                                output in ROW (1 : J2 - J1 + 1)
C
C  ENVIRONMENT:  VAX/VMS; FORTRAN 77
C
C  HISTORY:
C   05/24/91  D.Saunders   Original implementation.
C
C  AUTHOR:   David Saunders, Sterling Software/NASA Ames, Palo Alto, CA.
C
C----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER
     >   IDIM, I, J1, J2
      REAL
     >   A (IDIM, *), ROW (*)

C     Local variables:

      INTEGER
     >   J

C     Execution:

      DO 10, J = J1, J2
         ROW (J - J1 + 1) = A (I, J)
   10 CONTINUE

      RETURN
      END
