C+----------------------------------------------------------------------
C
      SUBROUTINE CROSS (A, B, C)
C
C PURPOSE: CROSS calculates the cross product A x B = C for vectors in
C          three-space.  The argument description is obvious.
C
C METHOD:  C is a vector perpendicular to the plane of A and B and such
C          that A, B and C form a right-handed system.  By definition,
C
C                                   |  i   j   k |
C                       A x B   =   | A1  A2  A3 |
C                                   | B1  B2  B3 |
C
C          where i, j, k are orthogonal unit vectors.  The magnitude of
C          A x B equals the area of the parallelogram with sides A and B.
C
C AUTHOR:  D.Saunders/R.Kennelly, NASA Ames, Mt. View, CA.  09/20/92
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      REAL
     >   A (3), B (3), C (3)
      
C     Execution:

      C (1) =   A (2) * B (3) - A (3) * B (2)
      C (2) = -(A (1) * B (3) - A (3) * B (1))
      C (3) =   A (1) * B (2) - A (2) * B (1)

      RETURN
      END
