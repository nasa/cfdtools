C+------------------------------------------------------------------------------
C
      SUBROUTINE BSAREA (C, TA, TB, AREA, NWORK, WORK, LUNERR, IER)

C ONE-LINER: B-Spline curve: AREA under it
C            - -             ----
C PURPOSE:
C
C        BSAREA estimates the area under a NURBS curve on the interval
C     [XA, XB] - that is, the definite integral with respect to X.  The
C     parameter values corresponding to XA and XB are assumed to be known.
C
C        Note that DT_NURBS utility DTPLAR determines the area ENCLOSED by
C     a curve, which is expected to be a closed curve.  DTPLAR is therefore
C     unsuited to some shapes such as airfoils with thick trailing edges.
C
C        Note also that CONAREA is much more efficient for the special
C     case of a conic represented by a (rational) quadratic Bezier curve.
C
C METHOD:
C                                   x=b                  t(b)
C        Use the fact that Integral y(x) dx  =  Integral y(t) x'(t) dt.
C                                   x=a                  t(a)
C     General purpose adaptive quadrature routine QUANC8RC does the rest.
C     The relative error tolerance chosen here, SQRT (machine epsilon),
C     is a portable compromise: too small leads to too many function
C     evaluations.
C
C ENVIRONMENT:
C     VAX/VMS; FORTRAN 77, with...
C     > IMPLICIT NONE
C     > Trailing ! comments
C     > Names up to 8 characters
C
C HISTORY:
C     05/28/92  D.A.Saunders  Initial implementation, for airfoil applications.
C     06/16/92    "    "      DERIVS (3, 2), not (2, 2), for the rational case.
C     02/01/94       "        Dispensed with DOUBLE PRECISION.  Compile with
C                             /REAL_LENGTH=64 or -r8 switches as necessary.
C
C AUTHOR: David Saunders, Sterling Software/NASA Ames, Mt. View, CA.
C
C ------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments
C     ---------

      REAL    C (*)     ! (I)   DT_NURBS B-spline curve vector:
                        !       C(1) = 1 (# parametric variables);
                        !       C(2) = 2 (for X and Y);
                        !       C(3) = polynomial order k (4 for cubics);
                        !       C(4) = # control points, N;
                        !       C(5) = pointer for efficient evaluations;
                        !       C(6 : 5 + N + k) = knot vector,
                        !       followed by X coords. of control pts.,
                        !       followed by Y coords. of control pts.

      REAL    TA, TB    ! (I)   Parametric limits of integration.

      REAL    AREA      ! (O)   Desired integral.

      INTEGER NWORK     ! (I)   NWORK >= K * (K + 4) for polynomials of order K.

      REAL    WORK (NWORK)      ! (S)   Work-space for DTSPDR.

      INTEGER LUNERR    ! (I)   Logical unit for error messages.

      INTEGER IER       ! (O)   Success/error code:
                        !        0 means no problem was encountered;
                        !        ? <to be completed>

C     Procedures
C     ----------

      EXTERNAL DTMCON   ! Provides machine constants in portable form.
      REAL     DTMCON
      EXTERNAL DTSPDR   ! B-spline curve function/derivative evaluator.
      EXTERNAL QUANC8RC ! Adapative quadrature utility.

C-------------------------------------------------------------------------------

C     Local constants
C     ---------------

      REAL
     >   ABSERR
      PARAMETER
     >  (ABSERR = 0.E+0) ! One of the quadrature tolerances

C     Local variables
C     ---------------

      INTEGER
     >   ISTAT, NUMFUN
      REAL
     >   ERREST, F, FLAG, DERIVS (3, 2), RELERR, T
      LOGICAL
     >   FIRST

C     Storage
C     -------

      DATA
     >   FIRST /.TRUE./
      SAVE
     >   FIRST, RELERR

C     Execution
C     ---------


      IF (FIRST) THEN
         FIRST = .FALSE.
         RELERR = SQRT (DTMCON (6)) ! SQRT (MACHINE EPS).  Not obvious what else
      END IF

      ISTAT = 0                     ! Initialize the integration
      T = TA

  100 CONTINUE

C        Use DTSPDR to evaluate the spline and its 1st derivs. at T:

         CALL DTSPDR (T, 1, C, WORK, NWORK, DERIVS, 3, IER)

         IF (IER .NE. 0) THEN
            WRITE (LUNERR, 1010) IER
            GO TO 999
         END IF

         F = DERIVS (2, 1) * DERIVS (1, 2)

         CALL QUANC8RC (T, F, TA, TB, ABSERR, RELERR, AREA,
     >                  ERREST, NUMFUN, FLAG, ISTAT)
         IF (ISTAT .GT. 0)
     >GO TO 100


      IF (ISTAT .NE. 0) THEN
         WRITE (LUNERR, 1020) ISTAT, NUMFUN, FLAG
         IER = ISTAT
      END IF


  999 RETURN


C     Formats.

 1010 FORMAT (' BSAREA:  Error in DTSPDR.  IER: ', I4)
 1020 FORMAT (' BSAREA:  Error in QUANC8RC.  ISTAT: ', I4,
     >        '  NUMFUN: ', I4, '   FLAG: ', F7.3)

      END
