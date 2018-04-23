C+----------------------------------------------------------------------
C
      PROGRAM SHIFTSCALE
C
C     PURPOSE:
C
C        SHIFTSCALE determines interactively the shift/scale factors which
C     transform interval [A, B] to [P, Q] (and back). It is simply a driver
C     for the subroutine GETXFORM, which uses the following formulae:
C
C        SHIFT = (B P - A Q) / (B - A)   and   SCALE = (Q - P) / (B - A)
C
C     SHIFTSCALE was written as an aid to using QPLOT's secondary Y-axis.
C
C     HISTORY: 10/05/91  DAS  Initial implementation.
C              01/24/03   "   Minor F90 upgrades.
C
C     AUTHOR:  David Saunders, ELORET/NASA Ames, Moffett Field, CA
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Constants:

      INTEGER, PARAMETER ::
     >   LUNCRT = 6, LUNKBD = 5

C     Variables:

      INTEGER
     >   N
      REAL
     >   A, B, P, Q, REALS (2), SCALE, SHIFT

C     Execution:

      WRITE (LUNCRT, '(/, A)')
     >  ' Determination of transformations between [A, B] and [P, Q]:'

      DO ! Until CR or EOF

         N = 2
         CALL RDREALS (LUNCRT, '$Enter [A, B] as a pair: ',
     >                 LUNKBD, N, REALS)
         IF (N <= 0) EXIT   ! 0=CR, -1=EOF

         A = REALS (1)
         B = REALS (2)

         N = 2
         CALL RDREALS (LUNCRT, '$Enter [P, Q] as a pair: ',
     >                 LUNKBD, N, REALS)
         IF (N <= 0) EXIT   ! 0=CR, -1=EOF

         P = REALS (1)
         Q = REALS (2)

         CALL GETXFORM (A, B, P, Q, SCALE, SHIFT)

         WRITE (LUNCRT, '('' For [A, B] --> [P, Q]:  SHIFT ='', G15.7,
     >      ''  SCALE ='', G15.7)') SHIFT, SCALE

         CALL GETXFORM (P, Q, A, B, SCALE, SHIFT)

         WRITE (LUNCRT, '('' For [P, Q] --> [A, B]:  SHIFT ='', G15.7,
     >      ''  SCALE ='', G15.7)') SHIFT, SCALE

      END DO

      END PROGRAM SHIFTSCALE
