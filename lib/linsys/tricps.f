C+----------------------------------------------------------------------
C
      SUBROUTINE TRICPS ( N, D, L, R, S )
C
C ACRONYM: TRIdiagonal system (Cyclic, Positive definite, Symmetric)
C          ---                 -       -                  -
C
C PURPOSE: TRICPS solves one system  A s = r,  where A is a symmetric,
C          diagonally-dominant, irreducible matrix with positive diag-
C          onal elements, non-zero sub-/super-diagonals,  and non-zero
C          lower-left/upper-right elements, and zeros elsewhere. It is
C          suited to problems involving cyclic/periodic boundary  con-
C          ditions.
C
C METHOD:  The LDL' form of the Cholesky factorization of  A  is used.
C          L  in this case is unit lower bidiagonal with its last  row
C          also filled in:
C
C                 x x     x       1         x         1 x     x
C                 x x x           x 1         x         1 x   x
C                   x x x     =     x 1         x         1 x x
C                     x x x           x 1         x         1 x
C                 x     x x       x x x x 1         x         1
C
C          It is assumed that only one right-hand side is involved for
C          the given A,  so no attempt is made to separate the factor-
C          ization from the solution using the triangular factors.
C
C ARGUMENTS:
C   ARG  DIM TYPE I/O/S DESCRIPTION
C    N    -    I    I   Order of the system (>=3)
C    D    N    R   I/S  Main diagonal elements of A
C    L    N    R   I/S  Off-diagonals of A in L(1:N-1); A(N,1) in L(N)
C    R    N    R    I   Right-hand side elements (see next item)
C    S    N    R    O   Solution vector; may be the same argument as R
C                       in which case R is destroyed.  D and L are de-
C                       stroyed in either case.
C
C ERROR HANDLING:  None. A divide by zero can mean the matrix is sing-
C                  ular but is more likely to mean unwitting errors in
C                  the call to TRICPS.
C
C ENVIRONMENT:  VAX/VMS FORTRAN 77
C
C DEVELOPMENT HISTORY:
C     DATE  INITIALS  DESCRIPTION
C   05/24/84   DAS    Initial design/coding, for a spline application.
C   05/22/86   RAK    Three divides by DI are worth replacing with one
C                     reciprocal and three multiplies in the main loop.
C
C AUTHOR: David Saunders, Informatics General, Palo Alto, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C ... Arguments:

      INTEGER  N
      REAL     D(N), L(N), R(N), S(N)

C ... Local variables:

      INTEGER  I
      REAL     DI, DINV, LIM1

C ... Executable statements:

C ... The off-diagonals are modified in place.   The modified diagonals
C     do not have to be stored for the back-solve if the matrix factors
C     are treated as (LD)L' and not L(DL').  This opens the way for the
C     filled-in row of factor L to overwrite the original diagonals. It
C     is then not possible to include the I = N-1 case in the following
C     loop ( D(N-1) being the exceptional element ). Repeating the code
C     was considered preferable to an IF test inside the loop.

      DINV = 1.E+0 / D(1)
      D(1) = L(N) * DINV
      S(1) = R(1) * DINV
      D(N) = D(N) - L(N) * D(1)
      S(N) = R(N) - L(N) * S(1)

      DO 20 I = 2, N-2
         LIM1   = L(I-1)
         L(I-1) = LIM1 * DINV
         DI     = D(I) - LIM1 * L(I-1)
         DINV   = 1.E+0 / DI
         D(I)   = - LIM1 * D(I-1) * DINV
         S(I)   = ( R(I) - LIM1 * S(I-1) ) * DINV
         D(N)   = D(N) - ( DI * D(I) ) * D(I)
         S(N)   = S(N) - ( DI * D(I) ) * S(I)
   20 CONTINUE

      I      = N-1
      LIM1   = L(I-1)
      L(I-1) = LIM1 * DINV
      DI     = D(I) - LIM1 * L(I-1)
      D(I)   = ( L(I) - LIM1 * D(I-1) ) / DI
      S(I)   = ( R(I) - LIM1 * S(I-1) ) / DI
      D(N)   = D(N) - ( DI * D(I) ) * D(I)
      S(N)   = ( S(N) - ( DI * D(I) ) * S(I) ) / D(N)
      L(I)   = 0.E+0

      DO 30 I = N-1, 1, -1
         S(I) = S(I) - L(I) * S(I+1) - D(I)*S(N)
   30 CONTINUE

      RETURN
      END
