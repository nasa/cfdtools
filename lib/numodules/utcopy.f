C***********************************************************************
C
      SUBROUTINE UTCOPY (NWORDS, RINPUT, ROUTPT)
C
C ACRONYM: UTility for COPYing 1-D arrays or making room for insertions.
C          --          ----
C PURPOSE: Copy one real array to another, or shift the elements of a
C          single array up or down, meaning the input and output arrays
C          can overlap.
C     
C NOTES:   NWORDS > 0  covers the normal situation,  where there is no
C          worry about overwriting needed values.  NWORDS < 0 was used
C          to indicate the more awkward situation  (requiring starting
C          from the high-subscript end) in order not to change the no.
C          of arguments and hence not affect existing applications. 
C
C ARGUMENTS:
C    ARG     DIM   TYPE I/O/S DESCRIPTION 
C   NWORDS    -      I    I   Indicates no. of real words to copy:
C                             NWORDS>0 means either the arrays are not
C                                      overlapping, as in
C                                      CALL UTCOPY (N, X, Y),
C                                      or elements are to be moved in
C                                      the smaller subscript direction:
C                                      CALL UTCOPY (N, X(3), X(1));
C                             NWORDS<0 means shift |NWORDS| words in
C                                      the larger subscript direction:
C                                      CALL UTCOPY (-N, X(1), X(3)).
C   RINPUT |NWORDS|  R   I/O  Array being copied.
C   ROUTPT |NWORDS|  R   I/O  Copy of RINPUT(*), possibly overlapping.
C
C HISTORY:
C     c. 1980  D.A.Saunders  Original utility for data acquisition.
C
C     01/24/11  D. Saunders  Handle real arrays only; the original
C                            was coded for integers, and advised
C                            doubling NWORDS if the arrays were
C                            8-byte reals - a really bad practice.
C
C AUTHOR:  David Saunders, Informatics General/NASA Ames Research Center
C
C***********************************************************************

      IMPLICIT NONE

C     Arguments:

      INTEGER :: NWORDS
      REAL    :: RINPUT(:), ROUTPT(:)

C     Local variables:

      INTEGER :: I, I1, I2, INC

C     Execution:

      IF (NWORDS > 0) THEN
         I1  =  1
         I2  =  NWORDS
         INC =  1
      ELSE
         I1  = -NWORDS
         I2  =  1
         INC = -1
      END IF

      DO I = I1, I2, INC
         ROUTPT(I) = RINPUT(I)
      END DO

      END SUBROUTINE UTCOPY
