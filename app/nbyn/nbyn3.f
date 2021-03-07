C+----------------------------------------------------------------------
C
      PROGRAM NBYN3
C
C     PURPOSE:
C
C        NBYN3 is a program intended for solving small square systems of
C     linear equations interactively (such as 3x3 systems). For a given
C     matrix, more than one right-hand-side is provided for (one RHS at
C     a time).
C
C     METHOD:
C
C        NBYN3 drives versions of DECOMP and SOLVE (LU decomposition and
C     solution with partial pivoting) without the condition number estimate.
C     See also NBYN and NBYN2.
C
C     HISTORY:
C
C     07/16/79  DAS  Initial implementation.
C     04/07/87  DAS  Tidied up for publication.
C     05/19/93  DAS  Converted from FORTRAN 66; I/O made more portable.
C     09/24/07  DAS  NBYN3 adapted to call Fortran 90 DECOMP & SOLVE.
C
C     AUTHOR:  David Saunders, Sterling Software/NASA Ames, Mt. View, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Constants:

      INTEGER
     >   LUNCRT, LUNKBD, NMAX
      PARAMETER
     >  (LUNCRT = 6, LUNKBD = 5, NMAX = 6)

C     Variables:

      INTEGER
     >   I, IP (NMAX), J, N
      DOUBLE PRECISION
     >   A (NMAX, NMAX), B (NMAX)
      CHARACTER
     >   ANS * 1

C     Execution:

   20 WRITE (LUNCRT, '(A, I2, A, $)')
     >   ' Enter order of square system (<=', NMAX,').  N = '

      READ  (LUNKBD, *, ERR=20, END=99) N
      IF (N .GT. NMAX .OR. N .LE. 1) GO TO 20

      WRITE (LUNCRT, '(A)') ' Enter N*N matrix by rows (free form):'

      DO 40 I = 1, N
         WRITE (LUNCRT, '(/, A, I2, A)') ' Row', I, ':'
         READ  (LUNKBD, *) (A (I, J), J = 1, N)
   40 CONTINUE

      WRITE (LUNCRT, '(A)') ' Entered matrix follows:'

      DO 50 I = 1, N
         WRITE (LUNCRT, '(1X, 1P, 5D24.14)') (A (I, J), J = 1, N)
   50 CONTINUE

C     Factorize the matrix:

      CALL DECOMP (N, NMAX, A, IP)

      WRITE (LUNCRT, '(A)')
      IF (IP(N) == 0) THEN
         WRITE (LUNCRT, '(A)') ' The system is singular.  Try again.'
         GO TO 20
      END IF

   60 CONTINUE
C     (Permit multiple right-hand-sides.)

      WRITE (LUNCRT, '(A)')
      WRITE (LUNCRT, '(A)') ' Enter RHS vector (as a row):'
      READ  (LUNKBD, *) (B (I), I = 1, N)

C     Solve the system:

      CALL SOLVE (N, NMAX, A, B, IP)

      WRITE (LUNCRT, '(A)') ' Solution follows:'
      WRITE (LUNCRT, '(1X, 1P, 5D24.14)') (B (I), I = 1, N)

   70 WRITE (LUNCRT, '(//, A, $)') ' Another RHS? (Y or N) '
      READ  (LUNKBD, '(A)', ERR=70, END=99) ANS

      IF (ANS .EQ. 'Y' .OR. ANS .EQ. 'y') GO TO 60

   99 CONTINUE

      END PROGRAM NBYN3
