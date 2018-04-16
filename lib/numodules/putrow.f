C+---------------------------------------------------------------------
C
      SUBROUTINE PUTROW (IDIM, I, J1, J2, A, ROW)
C
C  PURPOSE:
C     PUTROW updates the (indicated part of the) Ith row in a 2-D array.
C
C  ARGUMENTS:
C     ARG   TYPE I/O/S   DIM     DESCRIPTION
C     IDIM   I     I      -      Declared row dimension of A
C     I      I     I      -      Target row
C     J1,    I     I      -      First and last elements of target row
C     J2                         to be transferred
C     A      R     O  (IDIM, *)  Array for which elements of the Ith row
C                                are updated
C     ROW    R     O      -      Elements J1 : J2 of row I of A are
C                                updated from ROW (1 : J2 - J1 + 1)
C
C  ENVIRONMENT:  VAX/VMS; FORTRAN 77
C
C  HISTORY:
C   05/24/91  D.Saunders   Original implementation, to match GETROW.
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
         A (I, J) = ROW (J - J1 + 1)
   10 CONTINUE

      RETURN
      END
