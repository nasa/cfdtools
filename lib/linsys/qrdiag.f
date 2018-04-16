C+----------------------------------------------------------------------
C
      SUBROUTINE QRDIAG ( A, B, C, R, S, N )
C
C ACRONYM: QR factorization & solution of a triDIAGonal system
C          --                                  ----
C
C PURPOSE: QRDIAG solves one tridiagonal system which need not be diagonally
C          dominant, using the QR factorization rather than the standard LU
C          decomposition (which would require pivoting for stability).
C
C METHOD:  The lower diagonal elements are eliminated by a sequence of symmetric
C          plane rotations involving the diagonal elements.  The initial B(1)
C          and transformed B(2):B(N) are assumed to be non-zero.  The right-hand
C          side is transformed along the way.  The resulting upper tridiagonal
C          system is then solved straightforwardly.
C
C          This version is intended for the single-right-hand-side case.
C
C ARGUMENTS:
C   ARG  DIM TYPE I/O/S DESCRIPTION
C    A    N    R   I/S  Input with the lower diagonal in A(2):A(N);
C                       A(1):A(N-1) are reused for 2nd upper diagonal fill in.
C    B    N    R   I/S  Input with the main diagonal; overwritten with the
C                       transformed diagonal.
C    C    N    R   I/S  Input with the upper diagonal in C(1):C(N-1), and
C                       overwritten with the transformed upper diagonal.
C                       C(N) is used to avoid special handling of the (N-1)th
C                       elimination.
C    R    N    R    I   Input with the right-hand side vector; left intact
C                       unless R and solution S are the same locations.
C    S    N    R    O   Output with the solution; may be the same locations
C                       as right-hand-side R.
C    N    -    I    I   Order of the system.  N >= 3.
C
C ERROR HANDLING:  None.  Divide by zero can mean the matrix is singular
C                  but is more likely to mean unwitting errors in the call
C                  to QRDIAG.
C
C ENVIRONMENT:  VAX/VMS FORTRAN 77
C
C DEVELOPMENT HISTORY:
C   01/04/89   DAS   Initial implementation.
C
C AUTHOR: David Saunders, Sterling Software, Palo Alto, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT  NONE

C     Arguments:

      INTEGER   N
      REAL      A(N), B(N), C(N), R(N), S(N)

C     Local variables:

      INTEGER   I
      REAL      A2, B1, B2, C1, C2, CS, SN, T0, S1, S2

C     Execution:

      S(1) = R(1)        ! Leaves RHS intact (unless S and R are the same).
      C(N) = R(1)        ! Avoids undefined C(N) when I = N-1.

      DO 20, I = 1, N - 1

C        Find the rotation that eliminates a(i+1) via (transformed) b(i):

         B1 = B(I)
         A2 = A(I+1)
         T0 = SQRT ( B1 ** 2 + A2 ** 2 )
         CS = B1 / T0
         SN = A2 / T0

C        Modify affected elements of the LHS ...

         B(I) = T0

         C1 = C(I)
         B2 = B(I+1)
         C(I) = CS * C1 + SN * B2
         B(I+1) = SN * C1 - CS * B2

         C2 = C(I+1)
         A(I) = SN * C2
         C(I+1) = -CS * C2

C        ... and of the RHS:

         S1 = S(I)
         S2 = R(I+1)
         S(I) = CS * S1 + SN * S2
         S(I+1) = SN * S1 - CS * S2

   20 CONTINUE

C     Solve the transformed, upper tridiagonal system:

      S(N) = S(N) / B(N)
      S(N-1) = ( S(N-1) - C(N-1) * S(N) ) / B(N-1)

      DO 30, I = N-2, 1, -1
         S(I) = ( S(I) - C(I) * S(I+1) - A(I) * S(I+2) ) / B(I)
   30 CONTINUE

      RETURN
      END
