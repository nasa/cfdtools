C+------------------------------------------------------------------------------
C
      SUBROUTINE BSEVAL (MODE, C, NEW, NTARGET, TARGET, TLEFT, TRIGHT,
     >                   TRESULT, XYRESULT, LUNERR, IER)
C
C ONE-LINER: B-Spline curve EVALuation(s) (Ys for given Xs or Xs for given Ys)
C            - -            ----
C PURPOSE:
C
C        BSEVAL evaluates a given 2-space DT_NURBS-format B-spline curve
C     in one of two ways:  either Y for given X if MODE = 'YFORX', or X
C     for given Y if MODE = 'XFORY'.  More than one evaluation per call
C     is provided for (but see the next paragraph).  If NTARGET = 1, the
C     corresponding value of parametric variable T is returned as well
C     (in case some derivatives are required, for instance).  Otherwise,
C     TRESULT (*) is a probably-undesired argument, and output is suppressed.
C
C        The application is expected to know the range of the parametric
C     variable T within which the solution(s) lie.  This range must be
C     such that the target coordinate is monotonic in that range.
C     It may be advantageous to request just one evaluation per call if
C     the range of T can thereby be narrowed for successive calls.
C
C        A typical application would be the evaluation of both surfaces of
C     an airfoil at the same X abscissa(s).
C
C METHOD:
C
C        This is a zero-finding problem.  E.g., for a given Xtarg, we need
C     to find the T for which X(T) = Xtarg to within round-off error.
C     Then we return the corresponding T and Y(T).
C
C ENVIRONMENT:
C     VAX/VMS; FORTRAN 77, with...
C     > IMPLICIT NONE
C     > Trailing ! comments
C     > Names up to 8 characters
C
C HISTORY:
C     05/15/92  D.A.Saunders  Initial design and implementation, with
C                             estimation of the maximum thickness of an
C                             airfoil in mind.
C     04/02/93     "    "     Added ABSCISSA, etc., to error messages.
C     04/21/93     "    "     Suppressed return of TRESULT (*) if there
C                             is more than one such (prompted by BSAFY4X).
C     02/01/94       "        Dispensed with DOUBLE PRECISION.  Compile with
C                             /REAL_LENGTH=64 or -r8 switches as necessary.
C
C AUTHOR: David Saunders, Sterling Software/NASA Ames, Mt. View, CA.
C
C ------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments
C     ---------

      CHARACTER MODE * 5         ! (I)   Upper case switch:
                                 !       'XFORY' means evaluate X for given Y;
                                 !       'YFORX' means the converse.

      REAL     C (*)             ! (I)   DT_NURBS B-spline curve vector:
                                 !       C(1) = 1 (# parametric variables);
                                 !       C(2) = 2 (for X and Y);
                                 !       C(3) = polynomial order k (degree + 1);
                                 !       C(4) = # control points, N;
                                 !       C(5) = pointer for efficient evalns.;
                                 !       C(6 : 5 + N + k) = knot vector,
                                 !       followed by X coords. of control pts.,
                                 !       followed by Y coords. of control pts.

      LOGICAL NEW                ! (I)   .FALSE. means C (*) ISN'T a new spline.
                                 !       This may (some day) permit suppressing
                                 !       redundant error-checking of C(*).

      INTEGER NTARGET            ! (I)   No. of target abscissas requiring
                                 !       evaluations.

      REAL    TARGET (NTARGET)   ! (I)   Target value(s) of X or of Y, depending
                                 !       on MODE.

      REAL    TLEFT, TRIGHT      ! (I)   Left and right values of T bracketing
                                 !       the target abscissa(s).

      REAL    TRESULT (NTARGET)  ! (O)   Calculated value of T matching the
                                 !       given target abscissa and XYRESULT
                                 !       (suppressed if NTARGET > 1).
                           
      REAL    XYRESULT (NTARGET) ! (O)   Calculated X(s) or Y(s) matching the
                                 !       target Y(s) or X(s).

      INTEGER LUNERR             ! (I)   Logical unit for error messages.

      INTEGER IER                ! (O)   Success/error code:
                                 !        0 means no problem was encountered;
                                 !       12 means <to be completed>

C     Procedures
C     ----------

      EXTERNAL  DTSPVL           ! Evaluates a B-spline curve at the given T.
      EXTERNAL  ZERORC           ! Reverse-communication 1-D zero-finder (avoids
                                 ! common blocks)

C-------------------------------------------------------------------------------

C     Local constants
C     ---------------

      INTEGER
     >   MAXFUN, MAXK, NWORK
      REAL
     >   TOL
      PARAMETER
     >  (MAXFUN  = 50,             ! Limit on zero-finder's # evaluations
     >   MAXK    = 7,              ! Local storage handles polynomial deg. <= 6
     >   NWORK   = 5 * MAXK - 2,   ! This is why the degree is limited
     >   TOL     = 0.E+0)

C     Local variables
C     ---------------

      INTEGER
     >   I, I1, I2, ISTAT, NUMFUN
      REAL
     >   ABSCISSA, F, HOLD (13), T, WORK (NWORK), XY (2)
      LOGICAL
     >   SINGLET

C     Execution
C     ---------


C     Most of the error checking is left to the lower level utilities.

      IF (NEW) THEN
         IF (C (3) .GT. MAXK) THEN
            IER = 3
            WRITE (LUNERR, 1000) 'Local storage exceeded.'
            ISTAT = 0
            GO TO 900
         END IF

C ***    DTSPVL doesn't permit use of NEW at this stage.

      END IF


C     Distinguish the cases:

      IF (MODE (1 : 1) .EQ. 'X') THEN  ! Evaluate X for Y
         I1 = 2      ! Initially, find T for Y (2nd fn. of T)
      ELSE
         I1 = 1      ! ... else initially find T for X
      END IF
      I2 = 3 - I1   ! Points to the other coordinate once T is found


      SINGLET = NTARGET .EQ. 1

C     For each target abscissa ...

      DO 300, I = 1, NTARGET

         ABSCISSA = TARGET (I)

C        Find the value of T which matches the abscissa:

         ISTAT = 2          ! 2 initializes the zero-finding iteration
         NUMFUN = MAXFUN


  200    CONTINUE

            CALL ZERORC (TLEFT, TRIGHT, T, F, TOL, NUMFUN, 'BSEVAL',
     >                   -LUNERR, HOLD, ISTAT)

            IF (ISTAT .EQ. -1) THEN  ! Probable application error

               IER = 1
               WRITE (LUNERR, 1000) 'Iteration limit reached.',
     >            ABSCISSA, TLEFT, TRIGHT
               GO TO 900

            ELSE IF (ISTAT .LT. 0) THEN  ! Some other fatal error

               IER = 2
               WRITE (LUNERR, 1000) 'Fatal return from ZERORC.',
     >            ABSCISSA, TLEFT, TRIGHT
               GO TO 900

            ELSE IF (ISTAT .GT. 0) THEN  ! Evaluate the function

               CALL DTSPVL (T, C, WORK, NWORK, XY, IER)

            IF (IER .NE. 0) THEN
               WRITE (LUNERR, 1000) 'Fatal return from DTSPVL.',
     >            ABSCISSA, TLEFT, TRIGHT
               GO TO 900
            END IF

            F = XY (I1) - ABSCISSA
            GO TO 200

         ELSE  ! ISTAT = 0 (success).  Ensure X and Y are both at final T

            CALL DTSPVL (T, C, WORK, NWORK, XY, IER)

            IF (SINGLET) TRESULT (1) = T
            XYRESULT (I) = XY (I2)

         END IF

  300 CONTINUE   ! Next target, or done

      GO TO 999


C     Error handling:

  900 WRITE (LUNERR, 1010) ISTAT, IER


C     Termination:

  999 RETURN


C     Formats:

 1000 FORMAT (/, ' *** BSEVAL trouble: ', A,
     >        /, ' Abscissa, T(left), T(right):', /, 1P, 3E24.16)
 1010 FORMAT (' ISTAT:', I3, '  IER: ', I4)

      END
