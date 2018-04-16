C+----------------------------------------------------------------------
C
      SUBROUTINE HTDIS (XA, XB, D1, N, X, LUNERR, IER)
C
C  ACRONYM: Hyperbolic Tangent DIStribution (1-sided)
C           -          -       ---
C
C  PURPOSE:
C        HTDIS generates abscissas on the interval [XA, XB] using a
C     hyperbolic tangent method, where dX(I) = X(I+1) - X(I), with
C                                            N-1
C        dX(1) = D1 = X(2) - XA   and   XA + Sum dX(I) = XB
C                                            I=1
C  METHOD:
C        The mathematics were evidently developed by Marcel Vinokur
C     (Sterling Software/NASA Ames).  The result can be shown to be
C     equivalent to that of MODE=3 for the earlier routine EXPDIS2.
C     The latter employs a zero-finder to solve the relevant nonlinear
C     equation, while HTDIS uses bisection to be more self-contained,
C     and the formulation involves a different (but presumably equivalent)
C     nonlinear equation.
C     
C        X(N) - XB is forced arbitrarily close to zero.  The convergence
C     criteria are tight because the last interval may be the smallest,
C     yet it in effect absorbs all of the final error when X(N) is forced
C     to be XB exactly.
C
C  ARGUMENTS:
C    NAME    DIM   TYPE I/O/S DESCRIPTION
C   XA        -     R     I   Desired X(1); passing X(1) here is safe
C   XB        -     R     I   Desired X(N); passing X(N) here is safe
C                             Note that XA = 0 and XB = 1 can provide one
C                             normalized distribution for multiple reuse.
C   D1        -     R     I   Desired size of first interval, X(2) - X(1);
C                             must be less than or equal to uniform dX
C   N         -     I     I   Number of points; N > 2
C   X         N     R     O   Output distribution, with X(1) = XA and
C                             X(N) = XB, and X(2) - X(1) = D1.
C   LUNERR    -     I     I   LUNERR > 0 shows the iteration history;
C                             LUNERR < 0 suppresses it (but any error
C                                        message goes to -LUNERR)
C   IER       -     I     O   IER = 0 means no problems;
C                                   1 means the requested D1 was greater
C                                     than the uniform value;
C                                   2 means the iteration failed - must
C                                     have been other bad inputs
C
C  ERROR HANDLING:
C     See LUNERR and IER.  Failure handling is left to the calling
C     program, which can presumably reprompt for parameters and try again.
C
C  ENVIRONMENT:  VAX/VMS, FORTRAN 77 with IMPLICIT NONE extension
C
C  HISTORY:
C   ????????  M.Vinokur  Algorithm.
C   03/01/89  J.Melton   Initial implementation.
C   03/10/89  R.Langhi   Functionality simplified and clarified; a few
C                        refinements to the convergence test, etc.;
C                        safeguarded the initial bracketing iteration;
C                        documentation patterned after GEODIS.
C   11/27/90  D.Saunders Moved HTDIS from a grid-generation application
C                        to GRIDLIB for general use. (This was not done
C                        at the time of HTDIS2 because it was thought
C                        that EXPDIS2 included the same functionality.
C                        However, the formulation is different and the
C                        bisection method is more self-contained than
C                        EXPDIS2's use of a general-purpose zero-finder.)
C                        Inverted the nonlinear equation to avoid large
C                        numbers, and safeguarded the DEL = 0 case which
C                        corresponds to a uniform distribution.
C                        Also: the internal arithmetic is now done in
C                        DOUBLE PRECISION, as for EXPDIS2.  It turns
C                        out that the numerics are still strangely
C                        imprecise:  D1 is only approximately achieved.
C   05/04/91    "    "   Moved HTDIS from GRIDLIB to OBSOLETE because
C                        EXPDIS2 is more precise.  HTDIS could be made
C                        more precise the way HTDIS2 has been (with an
C                        outer iteration and a Newton inner iteration)
C                        but EXPDIS2 is simpler in this case.
C                        
C  AUTHOR:  John E. Melton, NASA/Ames Research Ctr., Moffett Field, CA
C
C-----------------------------------------------------------------------

C     Declarations.
C     -------------

      IMPLICIT NONE

C ... Arguments:

      INTEGER
     >   N, LUNERR, IER
      REAL
     >   XA, XB, D1, X (N)

C ... Local constants:

      INTEGER
     >   MXITER
      DOUBLE PRECISION
     >   BIG, EPS, HALF, ONE, SMALL, TOL, ZERO
      PARAMETER
     >  (MXITER=50, BIG=20.D+0, EPS=1.D-7, HALF=0.5D+0, ONE=1.D+0,
     >   SMALL=1.D-10, TOL=1.D-10, ZERO=0.D+0)

C ... Local variables:

      INTEGER
     >   I, ITER, LUNMSG
      DOUBLE PRECISION
     >   ARG, B, DEL, DL, DR, FRAC, FUNDEL, RANGE, RFNM1, RTNHD2

C     Execution.
C     ----------

      IER = 0

C ... Normalize the input first and last spaces:

      RANGE = DBLE (XB - XA)
      RFNM1 = ONE / DBLE (N - 1)
      B = DBLE (N - 1) * DBLE (D1) / RANGE

C ... Safeguard undefined cases.  TOL is too small in the following test.

      IF (ABS (B - ONE) .LT. EPS) THEN  ! B=1 ==> DEL=0 (uniform)
         
         DO 5, I = 1, N-1
            X (I) = XA + REAL (I-1) * D1
    5    CONTINUE
         GO TO 40

      ELSE IF (B .GT. ONE) THEN  ! No solution to  DEL / SINH (DEL) = B
         IER = 1
         GO TO 99

      END IF

C ... Use bisection to solve Eq. 8 on p. 229, vol. II of EAGLE manual:

      IF (LUNERR .GT. 0)
     >   WRITE (LUNERR,'('' ITER    DEL          F(DEL)       TOL'')')

      ITER = 0
      DL = SMALL
      DR = BIG

   10 CONTINUE
         DEL = HALF * (DL + DR)
         FUNDEL = DEL / SINH (DEL) - B

         IF (LUNERR .GT. 0)
     >      WRITE (LUNERR, '(I5, 1P, 3D13.5)') ITER, DEL, FUNDEL, TOL

         IF (ABS (FUNDEL) .LT. TOL) GO TO 20

         IF (FUNDEL .LT. ZERO) THEN
            DR = DEL
         ELSE
            DL = DEL
         ENDIF

         ITER = ITER + 1
         IF (ITER .LT. MXITER)
     >GO TO 10

      LUNMSG = ABS (LUNERR)
      WRITE (LUNMSG,'(/'' *** HTDIS:  Bisection failed ***''/
     >                 '' Number of iterations:'', I4/
     >                 '' DEL, F(DEL), TOL:'', 1P, 3D14.6/
     >                 '' No distribution computed.''/)')
     >                 MXITER, DEL, FUNDEL, TOL
      IER = 2
      GO TO 99

   20 CONTINUE

C ... Compute distribution, now that DEL has been found.
C ... Use Eq. 13 on p. 230, vol II of the EAGLE manual:

      RTNHD2 = ONE / TANH (HALF * DEL)

C*      X (1) = 0.E+0
C*      DO 30, I = 1, N-1
C*         ARG  = HALF * DEL * (REAL (I-1) * RFNM1 - ONE)
C*         X(I) = ONE + TANH (ARG) * RTNHD2
C*   30 CONTINUE      
C*      X(N) = ONE

      X (1) = XA
      DO 30, I = 2, N-1
         ARG  = HALF * DEL * (DBLE (I-1) * RFNM1 - ONE)
         FRAC = ONE + TANH (ARG) * RTNHD2
         X (I) = XA + REAL (FRAC * RANGE)
   30 CONTINUE      

   40 X (N) = XB

   99 RETURN
      END
