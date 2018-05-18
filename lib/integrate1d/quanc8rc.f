C+------------------------------------------------------------------------------
C
      SUBROUTINE QUANC8RC (XEVAL, FEVAL, A, B, ABSERR, RELERR, RESULT,
     >                     ERREST, NOFUN, FLAG, ISTAT)
C
C  ACRONYM: QUadrature/Adaptive/Newton-Cotes/8-panel/Reverse-Communication
C           --         -        -      -     -       -       -
C
C  PURPOSE:
C
C        QUANC8RC is a reverse-communication version of QUANC8 for estimating
C     the integral of a function of one variable on the interval [A, B] to a
C     specified tolerance using an adaptive scheme based upon the 8-panel
C     Newton-Cotes rule.
C
C        This version avoids the data-communication problems encountered with
C     functions which require more than a single argument X to be sufficiently
C     defined.  (QUANC8 forces the additional information to be passed via
C     a COMMON block, which is usually inconvenient.)
C
C
C  USAGE:
C
C        :   :
C        X = A       ! Initialize the integration
C        ISTAT = 0
C
C    100 CONTINUE
C
C           CALL FUN (X, F, ..., IER)      ! Or whatever it takes to evaluate
C                                          ! the function as F at X
C           IF (IER .NE. 0) THEN
C              <Handle an error in the function evaluation>
C           END IF
C
C           CALL QUANC8RC (X, F, A, B, ABSERR, RELERR, RESULT,
C       >                  ERREST, NOFUN, FLAG, ISTAT)
C
C           IF (ISTAT .GT. 0)
C       >GO TO 100
C
C        IF (ISTAT .NE. 0) THEN   ! Trouble
C           <Display a warning with FLAG = XXX.YYY
C            and handle the failure>
C        END IF
C        :    :
C
C
C  INPUT:
C
C     XEVAL   is not really an input, but it needs to be initialized to A
C             as shown above for all evaluations to be done in the same loop.
C     FEVAL   is input with the value of the function corresponding to XEVAL.
C     A       is the lower limit of integration.
C     B       is the upper limit of integration.  (B may be less than A.)
C     ABSERR  is an absolute error tolerance which should be non-negative.
C     RELERR  is a relative error tolerance which should be non-negative.
C     ISTAT   is input as 0 on the FIRST call and SHOULD NOT BE CHANGED by
C             the calling program before the integration is complete, since
C             a positive value output by the previous call is meaningful as
C             input on the present call if it is positive.
C
C  OUTPUT:
C
C     ISTAT   is output as 0 upon successful completion;
C             ISTAT < 0 indicates a failure - see FLAG;
C             ISTAT > 0 is convenient as an index into this routine's working
C             arrays of abscissas and ordinates: the index output is used
C             as input on the next call to indicate where to store FEVAL.
C     RESULT  is an approximation to the integral which should satisfy the
C             least stringent of the two error tolerances.
C     ERREST  is an estimate of the magnitude of the actual error.
C     NOFUN   is the number of function values used in calculation of RESULT.
C     FLAG    is a reliability indicator.  If FLAG is 0., then RESULT
C             probably satisfies the error tolerance.  If FLAG is XXX.YYY,
C             then XXX = the number of intervals which have not converged
C             and 0.YYY = the fraction of the interval left to do when the
C             (internal) limit on NOFUN was approached.
C
C
C  REFERENCE:
C
C     Forsythe, G., M. Malcolm, and C. Moler.  Computer Methods for
C     Mathematical Computations (Englewood Cliffs: Prentice-Hall, 1977).
C
C
C  HISTORY:
C
C     29 Apr. 1984   R.A.Kennelly,   Partial conversion of QUANC8 to
C                    NASA Ames       FORTRAN 77 using generic functions.
C     30 Apr. 1986   R.A.K.          Added IMPLICIT NONE.
C     21 May  1992   D.A.Saunders,   Adapted QUANC8 as QUANC8RC to avoid
C                    Sterling/       COMMON blocks.  Subscript 0 for X (*)
C                    NASA Ames       and F (*) simplifies the code.
C
C-------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER
     >   ISTAT, NOFUN
      REAL
     >   A, B, ABSERR, ERREST, FEVAL, FLAG, RELERR, RESULT, XEVAL

C     Local constants (not all of them):

      REAL
     >   HALF, ONE, ZERO
      PARAMETER
     >  (HALF = 0.5E+0, ONE = 1.E+0, ZERO = 0.E+0)

C     Local variables:

      INTEGER
     >   LEVMIN, LEVMAX, LEVOUT, NOMAX, NOFIN, LEV, NIM, I, J
      REAL
     >   W0, W1, W2, W3, W4, AREA, STONE, STEP, COR11, TEMP, QPREV,
     >   QNOW, QDIFF, QLEFT, ESTERR, TOLERR, QRIGHT (31),
     >   F (0 : 16), X (0 : 16), FSAVE (8, 30), XSAVE (8, 30)

C     Local storage:

      SAVE

C     Execution:


C     ***   STAGE 1 ***   GENERAL INITIALIZATION

      IF (ISTAT .EQ. 0) THEN

         LEVMIN = 1
         LEVMAX = 30
         LEVOUT = 6
         NOFUN  = 0
         NOMAX  = 5000
         NOFIN  = NOMAX - 8 * (LEVMAX - LEVOUT + 2 ** (LEVOUT + 1)) !NOFUN limit

         W0 =   3956.0E+0 / 14175.0E+0
         W1 =  23552.0E+0 / 14175.0E+0
         W2 =  -3712.0E+0 / 14175.0E+0
         W3 =  41984.0E+0 / 14175.0E+0
         W4 = -18160.0E+0 / 14175.0E+0

         FLAG   = ZERO
         RESULT = ZERO
         COR11  = ZERO
         ERREST = ZERO
         AREA   = ZERO

         IF (A .EQ. B) GO TO 99


C        ***   STAGE 2 ***   INITIALIZATION FOR FIRST INTERVAL

         LEV = 0
         NIM = 1
         STONE = (B - A) / 16.0E+0

         DO 20, J = 0, 16, 2
            X (J) = A + STONE * J
   20    CONTINUE
         X (16) = B     ! To avoid round-off

         QPREV  = ZERO
      END IF

      IF (MOD (ISTAT, 2) .EQ. 0) THEN   ! Fill in an even-numbered F (*)

         NOFUN = NOFUN + 1
         F (ISTAT) = FEVAL

         IF (ISTAT .LT. 16) THEN
            ISTAT = ISTAT + 2
            XEVAL = X (ISTAT)
            GO TO 99                    ! Go back for another function value
         ELSE
            ISTAT = -1                  ! Initialize the odd function evals.
         END IF

      END IF


C     ***   STAGE 3 ***   CENTRAL CALCULATION
C     Requires QPREV, X0, X2, X4, ..., X16, F0, F2, F4, ..., F16.
C     Calculates X1, X3, ..., X15, F1, F3, ..., F15, QLEFT, QRIGHT, QNOW,
C     QDIFF, and AREA.


   30 IF (ISTAT .EQ. -1) THEN

         DO 35, J = 1, 15, 2
            X (J) = (X (J - 1) + X (J + 1)) * HALF
   35    CONTINUE

         XEVAL = X (1)
         ISTAT = 1
         GO TO 99                    ! Go get F (1)

      END IF

      NOFUN = NOFUN + 1
      F (ISTAT) = FEVAL

      IF (ISTAT .LT. 15) THEN
         ISTAT = ISTAT + 2
         XEVAL = X (ISTAT)
         GO TO 99                    ! Go get another odd function value
      END IF

      STEP = (X (16) - X (0)) / 16.0E+0
      QLEFT = (W0 * (F (0) + F (8)) + W1 * (F (1) + F (7)) +
     >         W2 * (F (2) + F (6)) + W3 * (F (3) + F (5)) +
     >         W4 * F (4)) * STEP
      QRIGHT (LEV + 1) = (W0 * (F (8) + F (16)) + W1 * (F (9) + F (15))
     >                 +  W2 * (F (10)+ F (14)) + W3 * (F (11)+ F (13))
     >                 +  W4 * F (12)) * STEP
      QNOW = QLEFT + QRIGHT (LEV + 1)
      QDIFF = QNOW - QPREV
      AREA = AREA + QDIFF

C     ***   STAGE 4 *** INTERVAL CONVERGENCE TEST

      ESTERR = ABS (QDIFF) / 1023.0E+0
      TOLERR = MAX (ABSERR, RELERR * ABS (AREA)) * (STEP / STONE)
      IF (LEV .LT. LEVMIN) GO TO 50
      IF (LEV .GE. LEVMAX) GO TO 65
      IF (NOFUN .GT. NOFIN) GO TO 60
      IF (ESTERR .LE. TOLERR) GO TO 70

C     ***   STAGE 5   ***   NO CONVERGENCE
C     Locate next interval.

   50 NIM = 2 * NIM
      LEV = LEV + 1

C     Store right hand elements for future use:

      DO 52, I = 1, 8
         FSAVE (I, LEV) = F (I + 8)
         XSAVE (I, LEV) = X (I + 8)
   52 CONTINUE

C     Assemble left hand elements for immediate use:

      QPREV = QLEFT
      DO 55, I = 1, 8
         J = -I
         F (2 * J + 18) = F (J + 9)
         X (2 * J + 18) = X (J + 9)
   55 CONTINUE

      ISTAT = -1
      GO TO 30


C     ***   STAGE 6   ***   TROUBLE HANDLING
C     The number of function values is about to exceed limit.

   60 NOFIN = NOFIN + NOFIN
      LEVMAX = LEVOUT
      FLAG = FLAG + (B - X (0)) / (B - A)
      GO TO 70

C     CURRENT LEVEL IS LEVMAX.

   65 FLAG = FLAG + ONE

C     ***   STAGE 7   ***   INTERVAL CONVERGED
C     Add contributions into the running sums:

   70 RESULT = RESULT + QNOW
      ERREST = ERREST + ESTERR
      COR11  = COR11  + QDIFF / 1023.0E+0

C     Locate the next interval:

   75 IF (NIM .NE. 2 * (NIM / 2)) THEN
         NIM = NIM / 2
         LEV = LEV - 1
         GO TO 75
      END IF

      NIM = NIM + 1
      IF (LEV .LE. 0) GO TO 80

C     Assemble elements required for the next interval:

      QPREV = QRIGHT (LEV)
      X (0) = X (16)
      F (0) = F (16)
      DO 78, I = 1, 8
         F (2 * I) = FSAVE (I, LEV)
         X (2 * I) = XSAVE (I, LEV)
   78 CONTINUE

      ISTAT = -1
      GO TO 30

C     ***   STAGE 8   ***   FINALIZE AND RETURN

   80 RESULT = RESULT + COR11

C     Make sure ERREST is not less than roundoff level:

      IF (ERREST .NE. ZERO) THEN
   90    TEMP = ABS (RESULT) + ERREST
         IF (TEMP .EQ. ABS (RESULT)) THEN
            ERREST = ERREST + ERREST
            GO TO 90
         END IF
      END IF

C     Set final ISTAT according to original error flag:

      IF (FLAG .EQ. ZERO) THEN
         ISTAT = 0
      ELSE
         ISTAT = -2
      END IF

   99 RETURN
      END
