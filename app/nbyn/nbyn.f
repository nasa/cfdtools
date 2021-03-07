C+----------------------------------------------------------------------
C
      PROGRAM NBYN
C
C
C     PURPOSE:
C
C        NBYN is a program intended for solving small square systems of
C     linear equations interactively (such as 3x3 systems). For a given
C     matrix, more than one right-hand-side is provided for (one RHS at
C     a time).  Condition number estimates are also provided.
C
C
C     METHOD:
C
C        NBYN drives versions of DECOMP and SOLVE (LU decomposition and
C     solution with partial pivoting) from Forsythe, et al. Double pre-
C     cision is used just to be on the safe side. Thus this program can
C     be useful in the presence of exact data,  although typical inputs
C     are expected to be significant to a few digits only, and the full
C     precision used in displaying the solution should  be  interpreted
C     appropriately.
C
C     ENVIRONMENT:  DEC VAX/VMS and FORTRAN 77.
C
C     HISTORY:
C
C     07/16/79  DAS  Initial implementation.
C     04/07/87  DAS  Tidied up for publication.
C     05/19/93  DAS  Converted from FORTRAN 66; I/O made more portable.
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
     >   A (NMAX, NMAX), B (NMAX), WORK (NMAX), COND
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

      CALL DECOMP (NMAX, N, A, COND, IP, WORK)

      WRITE (LUNCRT, '(A)')
      IF (COND .EQ. 1.D+32) THEN
         WRITE (LUNCRT, '(A)') ' The system is singular.  Try again.'
         GO TO 20

      ELSE
         WRITE (LUNCRT, '(A, 1P, D12.3)')
     >      ' Condition number estimate:', COND
         IF (COND + 1.D+0 .EQ. COND) THEN
            WRITE (LUNCRT, '(A)')
     >         ' WARNING: Matrix is singular to working precision.'
         END IF
      END IF

   60 CONTINUE
C     (Permit multiple right-hand-sides.)

      WRITE (LUNCRT, '(A)')
      WRITE (LUNCRT, '(A)') ' Enter RHS vector (as a row):'
      READ  (LUNKBD, *) (B (I), I = 1, N)

C     Solve the system:

      CALL SOLVE (NMAX, N, A, B, IP)

      WRITE (LUNCRT, '(A)') ' Solution follows:'
      WRITE (LUNCRT, '(1X, 1P, 5D24.14)') (B (I), I = 1, N)

   70 WRITE (LUNCRT, '(//, A, $)') ' Another RHS? (Y or N) '
      READ  (LUNKBD, '(A)', ERR=70, END=99) ANS

      IF (ANS .EQ. 'Y' .OR. ANS .EQ. 'y') GO TO 60

   99 STOP ' '
      END
