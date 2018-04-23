C+----------------------------------------------------------------------
C
      PROGRAM NBYN2
C
C     07/09/93  DAS  Version of NBYN for testing LUSOLVE.
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
     >   I, IER, J, N
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

      WRITE (LUNCRT, '(A)')
      WRITE (LUNCRT, '(A)') ' Enter RHS vector (as a row):'
      READ  (LUNKBD, *) (B (I), I = 1, N)

C     Solve the system:

      CALL LUSOLVE (N, NMAX, A, B, IER)

      WRITE (LUNCRT, '(A)')
      IF (IER .NE. 0) THEN
         WRITE (LUNCRT, '(A)') ' The system is singular.  Try again.'
         GO TO 20
      END IF

      WRITE (LUNCRT, '(A)') ' Solution follows:'
      WRITE (LUNCRT, '(1X, 1P, 5D24.14)') (B (I), I = 1, N)

   70 WRITE (LUNCRT, '(//, A, $)') ' Another system? (Y/N/^Z) '
      READ  (LUNKBD, '(A)', ERR=70, END=99) ANS

      IF (ANS .EQ. 'Y' .OR. ANS .EQ. 'y') THEN
         WRITE (LUNCRT, '(A)')
         GO TO 20
      END IF

   99 STOP ' '
      END
