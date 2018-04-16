C+----------------------------------------------------------------------
C
      SUBROUTINE GETBLOCK (IDIMA, JDIMA, KDIMA, A, I1,I2, J1,J2, K1,K2,
     >                     IDIMB, JDIMB, KDIMB, B)
C
C  PURPOSE:
C
C        GETBLOCK extracts the indicated part of array A into array B,
C     where both are considered as 3-dimensional.  The initial application
C     was to extracting I, J, or K planes from PLOT3D files.  (Treating a
C     plane as a degenerate case of a block simplifies the code greatly.)
C
C  ARGUMENTS:
C
C     ARG     TYPE I/O/S   DIM     DESCRIPTION
C
C     IDIMA,   I   I        -      Dimensions of A in calling program.
C     JDIMA,
C     KDIMA
C     A        R   I    As shown.  3-D array from which a block or plane
C                                  is to be extracted.
C     I1, I2,  I   I        -      First and last elements of A to be
C     J1, J2,                      extracted in each dimension.
C     K1, K2
C     IDIMB,   I   I        -      Dimensions of B in calling program.
C     JDIMB,
C     KDIMB
C     B        R     O  As shown.  Elements B (1:NI, 1:NJ, 1:NK) are
C                                  filled, where NI = I2 - I1 + 1, etc.
C
C  ERROR HANDLING:   None.
C
C  ENVIRONMENT:  VAX/VMS, FORTRAN 77 with minor extensions:
C                IMPLICIT NONE and 8-character module name
C
C  HISTORY:
C   12/29/91  D.A.Saunders  Original implementation, for FMAP.
C
C  AUTHOR:    David Saunders, Sterling Software/NASA Ames, Mt. View, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER
     >   IDIMA, JDIMA, KDIMA, I1, I2, J1, J2, K1, K2,
     >   IDIMB, JDIMB, KDIMB
      REAL
     >   A (IDIMA, JDIMA, KDIMA), B (IDIMB, JDIMB, KDIMB)

C     Local variables:

      INTEGER
     >   I, J, K, IDIF, JDIF, KDIF

C     Execution:

      IDIF = I1 - 1
      JDIF = J1 - 1
      KDIF = K1 - 1

      IF (IDIF + JDIF + KDIF .EQ. 0) THEN  ! The common case
         DO 30, K = 1, K2
            DO 20, J = 1, J2
               DO 10, I = 1, I2
                  B (I, J, K) = A (I, J, K)
   10          CONTINUE
   20       CONTINUE
   30    CONTINUE
      ELSE
         DO 60, K = K1, K2
            DO 50, J = J1, J2
               DO 40, I = I1, I2
                  B (I - IDIF, J - JDIF, K - KDIF) = A (I, J, K)
   40          CONTINUE
   50       CONTINUE
   60    CONTINUE
      END IF

      RETURN
      END
