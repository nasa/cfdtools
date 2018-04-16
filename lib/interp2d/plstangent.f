C+------------------------------------------------------------------------------
C
      SUBROUTINE PLSTANGENT (N, X, Y, T, IL, IR, CALCT, SLOPE, TOL,
     >                       TSLOPE, IT, LUNOUT, IER)
C
C     ONE-LINER: Parametric local spline:  Find T for given slope dY/dX
C
C     DESCRIPTION:
C
C        PLSTANGENT estimates the value of the parametric variable T where
C     a 2-space curve is touched by a tangent line of specified slope.
C
C        Reverse-communication zero-finder ZERORC is used to solve the
C     nonlinear equation  F(T) = Y'(T) - SLOPE * X'(T) = 0.  The first
C     derivatives are provided by local spline utility LCSFIT.  This
C     approach was considered less cumbersome than possible use of second
C     derivatives X"(T) and Y"(T) via conventional spline utilities
C     CSFIT/CSEVAL, which would have allowed a safeguarded Newton iteration.
C
C     USAGE NOTES:
C        Large or infinite SLOPEs present difficulties for the intended
C     application to tangents near an airfoil leading edge.  The application
C     may be advised to switch the roles of X and Y if that provides more
C     manageable magnitudes to work with.
C
C     HISTORY:
C        10/15/02  D.A.Saunders  Initial adaptation of PLXCUT (T for given X).
C
C     AUTHOR:
C        David Saunders, ELORET/NASA Ames Research Center, Moffett Field, CA.
C
C ------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) ::
     >    N                 ! Number of points on curve; N >= 2

      REAL, INTENT (IN), DIMENSION (N) ::
     >    X, Y              ! Curve coordinates

      REAL, INTENT (INOUT) ::
     >   T (N)              ! Cumulative chord lengths, either input
                            ! or computed here - see CALCT

      INTEGER, INTENT (IN) ::
     >   IL, IR             ! Curve indices bracketing the solution

      LOGICAL, INTENT (IN) ::
     >   CALCT              ! .FALSE. means the calling program has
                            !         supplied T (*)

      REAL, INTENT (IN) ::
     >   SLOPE,             ! Slope of tangent; see USAGE NOTES
     >   TOL                ! Tolerance on || f || used relative to the
                            ! curve length; try 1.E-6 or 1.E-14 for 32-
                            ! or 64-bit arithmetic.

      REAL, INTENT (OUT) ::
     >   TSLOPE             ! Estimated T where the tangent touches the curve.

      INTEGER, INTENT (OUT) ::
     >   IT                 ! Actual index nearest to the point of
                            ! tangency in INTERVAL's sense (that
                            ! is, between 1 and N - 1 inclusive)

      INTEGER, INTENT (IN) ::
     >   LUNOUT             ! Logical unit for displaying iterations,
                            ! which are suppressed if LUNOUT < 0.
                            ! |LUNOUT| is used for error messages.

      INTEGER, INTENT (OUT) ::
     >   IER                !   0 means no problem was encountered;
                            ! < 0 means trouble - see ISTAT in ZERORC
C     Procedures:

      EXTERNAL
     >   CHORDS2D,          ! Chord-length utility
     >   INTERVAL,          ! Search utility, also used by LCSFIT
     >   LCSFIT,            ! Local spline utility; "loose" fits are used
     >   ZERORC             ! Reverse-communication no-derivatives zero-finder

C-------------------------------------------------------------------------------

C     Local constants.

      INTEGER, PARAMETER ::
     >   MAXFUN = 40        ! Limit on number of function evaluations

      REAL, PARAMETER ::
     >   BIG  = 1.E+25,
     >   HALF = 0.5E+0,
     >   ONE  = 1.0E+0,
     >   ZERO = 0.0E+0

      LOGICAL, PARAMETER ::
     >   NEW = .TRUE.       ! Because the data change with each LCSFIT call

      CHARACTER, PARAMETER ::
     >   METHOD * 1 = 'B',  ! "Bessel" = "loose" spline fits
     >   NAME * 10 = 'PLSTANGENT'

C     Local variables:

      INTEGER
     >   ISTAT, LUNERR, NUMFUN

      REAL
     >   FI, TI, TL, TR, XI, XP, YI, YP, TOLER, TTOTAL, HOLD (13)

C     Execution:

      LUNERR = ABS (LUNOUT)

C     Set up the cumulative chord lengths unless they are supplied
C     by the higher level for efficiency reasons:

      IF (CALCT) THEN ! Don't normalize

         CALL CHORDS2D (N, X, Y, .FALSE., TTOTAL, T)

      ELSE
         TTOTAL = T (N)
      END IF

      TL     = T (IL)
      TR     = T (IR)
      TOLER  = TOL * TTOTAL
      NUMFUN = MAXFUN
      ISTAT  = 2          ! Initialize the iteration

      DO WHILE (ISTAT > 0)

         CALL ZERORC (TL, TR, TI, FI, TOLER, NUMFUN, NAME,
     >                LUNOUT, HOLD, ISTAT)

         IF (ISTAT > 0) THEN  ! Evaluate the function

            CALL LCSFIT (N, T, X, NEW, METHOD, 1, TI, XI, XP)
            CALL LCSFIT (N, T, Y, NEW, METHOD, 1, TI, YI, YP)

            FI = YP - SLOPE * XP

         END IF
            
      END DO


C     Wrap up:

      IER = ISTAT

      IF (ISTAT >= -1) THEN ! 0 = converged); -1 = iteration limit - may be OK

         TSLOPE = TI

         CALL INTERVAL (N, T, TI, ONE, IT) ! Recover the "left" index IT

      ELSE ! ISTAT < -1) THEN      ! Fatal error

         WRITE (LUNERR, 1010) 'Bad return from ZERORC.  ISTAT:', ISTAT

      END IF

      RETURN


C     Formats:

 1010 FORMAT (/, ' PLSTANGENT: ', A, I3)

      END SUBROUTINE PLSTANGENT
