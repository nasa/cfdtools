C+------------------------------------------------------------------------------
C
      SUBROUTINE CONAREA (C, AREA, LUNERR, IER)

C ONE-LINER: Area under a conic arc (rational quadratic Bezier curve)
C
C PURPOSE:
C
C        CONAREA is a specialized form of BSAREA, intended to be much more
C     efficient for the case of a rational quadratic Bezier curve in 2-space.
C     The intended application is to quarter-sections of bodies such as
C     aircraft fuselages.
C
C        The arc is expected in DT_NURBS format, C(*), and the full length
C     of the arc is integrated.  The special-special case of a parabola 
C     (non-rational) is handled here by making sure the weights are equal.
C
C METHOD:
C
C        The integral with respect to x of a rational quadratic y=y(x),
C     where x=x(u) and y=y(u) are given, is still intractable analytically,
C     so adaptive quadrature is used as for the general case.  This can
C     involve hundreds of function and derivative evaluations for a single
C     calculation.  Area ruling of a fuselage may require hundreds of these
C     area calculations - hence this special version.  As for BSAREA, we
C                                 x=b                  t(b)
C     use the fact that  Integral y(x) dx  =  Integral y(t) x'(t) dt.
C                                 x=a                  t(a)
C     General purpose adaptive quadrature routine QUANC8RC does the rest.
C     The relative error tolerance chosen here, SQRT (machine epsilon), is
C     a portable compromise: too small leads to too many function evaluations.
C
C ERROR HANDLING:  C(*) is assumed to be checked at a higher level.
C                  Failure of the adaptive quadrature is trapped - see IER.
C
C ENVIRONMENT:
C     VAX/VMS; FORTRAN 77, with...
C     > IMPLICIT NONE
C     > Trailing ! comments
C     > Names up to 8 characters
C
C HISTORY:
C     06/16/92  D.A.Saunders  Initial adaptation of BSAREA.
C     02/01/94       "        Dispensed with DOUBLE PRECISION.  Compile with
C                             /REAL_LENGTH=64 or -r8 switches as necessary.
C
C AUTHOR: David Saunders, Sterling Software/NASA Ames, Mt. View, CA.
C
C ------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments
C     ---------

      REAL    C (*)     ! (I)   DT_NURBS-format Bezier curve:
                        !       C(1) = 1 (# parametric variables);
                        !       C(2) =-3 (or 2) (2-space/(non)rational)
                        !       C(3) = 2 (for quadratic);
                        !       C(4) = 3 (# control points);
                        !       C(5) = 3 (pointer into knot vector);
                        !       C(6 :11) = knots {0,0,0,1,1,1} or equiv.
                        !       C(12:14) = (weighted) Xs of control pts.,
                        !       C(15:17) = (weighted) Ys of control pts.,
                        !       C(18:20) = weights {1,w,1} or equivalent
                        !                  (not reqd. if C(2)=2).

      REAL    AREA      ! (O)   Desired integral.

      INTEGER LUNERR    ! (I)   Logical unit for error messages.

      INTEGER IER       ! (O)   Success/error code:
                        !       0 means no problem was encountered;
                        !      <0 means trouble within QUANC8RC - see
                        !         its description of ISTAT (=IER here)
C     Procedures
C     ----------

      EXTERNAL DTMCON   ! Provides machine constants in portable form.
      REAL     DTMCON
      EXTERNAL QUANC8RC ! Adapative quadrature utility.

C-------------------------------------------------------------------------------

C     Local constants
C     ---------------

      REAL
     >   ONE, ZERO
      PARAMETER
     >  (ONE  = 1.E+0,
     >   ZERO = 0.E+0)       ! Used for the quadrature tolerance ABSERR

C     Local variables
C     ---------------

      INTEGER
     >   ISTAT, NUMFUN
      REAL
     >   ERREST, F, FLAG, RELERR,
     >   D, DP, ONEMU, U, UA, UB, W1, W2, W3, XN, XP, YN
      LOGICAL
     >   FIRST

C     Storage
C     -------

      DATA
     >   FIRST /.TRUE./
      SAVE
     >   FIRST

C     Execution
C     ---------


      IF (FIRST) THEN
         FIRST = .FALSE.
         RELERR = SQRT (DTMCON (6)) ! SQRT (MACHINE EPS).  Not obvious what else
      END IF

C     Handle the nonrational case:

      IF (C (2) .LT. ZERO) THEN
         W1 = C (18)                ! Probably 1., but may not be
         W2 = C (19)
         W3 = C (20)
      ELSE
         W1 = ONE
         W2 = ONE
         W3 = ONE
      END IF

C     Initialize the integration:

      UA = C (6)                    ! First knot, probably 0.
      U  = UA
      UB = C (11)                   ! Last knot, probably 1.
      ISTAT = 0

  100 CONTINUE

C        We need y(u) and x'(u) only (and x(u), but not y'(u)).
C        Let    x = xn / d, y = yn / d
C        Then  x' = xn' / d - xn.d' / d.d
C        Everything is interpreted as a linear interpolation:

         ONEMU = ONE - U
         XN = ONEMU * (ONEMU * C (12)   + U * C (13)) +
     >            U * (ONEMU * C (13)   + U * C (14))
         YN = ONEMU * (ONEMU * C (15)   + U * C (16)) +
     >            U * (ONEMU * C (16)   + U * C (17))
         D  = ONEMU * (ONEMU * W1       + U * W2)     +
     >            U * (ONEMU * W2       + U * W3)
         XP = ONEMU * (C (13) - C (12)) + U * (C (14) - C (13)) ! 2. incl. below
         DP = ONEMU * (W2     - W1)     + U * (W3     - W2)     ! "    "    "

         F = (YN / D ** 2) * (XP - XN * DP / D)

         CALL QUANC8RC (U, F, UA, UB, ZERO, RELERR, AREA,
     >                  ERREST, NUMFUN, FLAG, ISTAT)
         IF (ISTAT .GT. 0)
     >GO TO 100


      IF (ISTAT .EQ. 0) THEN  ! Success
         AREA = AREA + AREA   ! Since a factor of 2 was omitted on XP and DP
      ELSE
         WRITE (LUNERR, 1000) ISTAT, NUMFUN, FLAG
         IER = ISTAT
      END IF


      RETURN


C     Formats
C     -------

 1000 FORMAT (/, ' CONAREA: Error in QUANC8RC.  ISTAT: ', I3,
     >        '   NUMFUN: ', I4, '   FLAG: ', F7.3)

      END
