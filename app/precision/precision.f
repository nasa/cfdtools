C+----------------------------------------------------------------------
C
      PROGRAM PRECISION
C
C  PURPOSE:
C        This program drives the function-precision estimation routine
C     OBJEPS for a given function with calling sequence of the form
C
C        CALL xxx (N, X, F)                    ! 'Q' for QNMDIF
C        CALL xxx (MODE, N, X, F, G, NSTATE)   ! 'N' for NPSOL/NZSOL
C
C     Actually, the more recent forms, QNMDIF2 and NPOPT, now both expect
C
C        CALL xxx (N, X, F, NCALL, FAIL)
C
C     now, and OBJEPS should be adjusted accordingly (no call for that
C     yet at NASA ARC).
C
C     See the optlib collection for OBJEPS.
C
C     This version examines the function F8PT3 from Example 8.3 of the
C     book PRACTICAL OPTIMIZATION.
C
C     HISTORY:
C     05/09/91  D.A.Saunders  Simple initial test program.
C     05/10/91    "     "     OBJEPS had FUNTYPE argument added.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Constants:

      INTEGER
     >   KMAX, LUNOUT, N
      PARAMETER
     >  (KMAX = 7, LUNOUT = 6, N = 1)

C     Variables:

      INTEGER
     >   I, IER
      REAL
     >   FEPS (N), H (N), X (N)
      EXTERNAL
     >   F8PT3

C     Execution:

      OPEN (UNIT=LUNOUT, STATUS='NEW', FILE='precision.out')
      WRITE (LUNOUT, '(1X, A)')
     >   'Precision estimation for a function of N variables',
     >   '--------------------------------------------------'

      H (1) = 1.E-2
      X (1) = 10.E+0 + 3.0 * H (1)   ! Book's table is not centered

      CALL OBJEPS (N, X, H, F8PT3, 'Q', KMAX, LUNOUT, FEPS, IER)

      WRITE (LUNOUT, '(A)') '1'
      H (1) = 1.E-3
      X (1) = 1.115E+0 + 3.0 * H (1)

      CALL OBJEPS (N, X, H, F8PT3, 'Q', KMAX, LUNOUT, FEPS, IER)

      STOP ' '

      END
