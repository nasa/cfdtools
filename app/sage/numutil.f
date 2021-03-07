C+---------------------------------------------------------------------
C
      SUBROUTINE CSDVAL (NX, X, Y, NU, U, B, C, D, S, SP, SPP)
C
C     One-liner: Cubic spline evaluation on an array, with derivatives.
C     ----------
C
C     Purpose:
C     --------
C
C        CSDVAL evaluates the 1-dimensional cubic spline function
C
C     S(U) = Y(I) + B(I)*(U - X(I)) + C(I)*(U - X(I))**2 + D(I)*(U - X(I))**3
C
C     along with its first and second derivatives
C
C     SP(U)  =      B(I)          + 2*C(I)*(U - X(I))  + 3*D(I)*(U - X(I))**2
C     SPP(U) =    2*C(I)          + 6*D(I)*(U - X(I))
C
C     using Horner's Rule for polynomials.  If U < X(1), then I = 1 is used;
C     if U > X(NX) then I = NX-1 is used. The data must be monotonic, but
C     may be increasing or decreasing.  (The coefficients associated with
C     the NXth point are not required.)  Most of the effort is in deciding
C     which set of spline coefficients to use - an "interpolation search"
C     is employed which should be highly efficient if the knots are fairly
C     uniformly spaced (commonly the case).
C
C        CSDVAL was developed as a companion routine to CSFIT, but it should
C     work with any method of generating cubic spline coefficients.  CSEVAL
C     should be used if derivatives are not required.  See also NOTES.
C
C     Arguments:
C     ----------
C
C     Name     Dimension  Type  I/O/S  Description
C     NX                   I    I      Number of data points defining the
C                                      spline ("knots"); >= 2.
C
C     X        NX          R    I      Abscissas of knots. Must be monotone
C                                      increasing or decreasing (consistent
C                                      with the data used to fit the spline).
C
C     Y        NX          R    I      Ordinates of knots.
C
C     NU                   I    I      Number of points at which to evaluate
C                                      the spline.
C
C     U        NU          R    I      Abscissas at which to evaluate spline.
C
C     B,C,D    NX          R    I      Spline coefficients, e.g. as computed
C                                      by CSFIT.
C
C     S       NU           R      O    Spline values at points U(I).
C
C     SP      NU           R      O    1st derivative values at points U(I).
C
C     SPP     NU           R      O    2nd derivative values at points U(I).
C
C     External modules:
C     -----------------
C
C     INTERVAL  Interpolation search for interval containing a point.
C
C     Environment:  Digital VAX-11/780 VMS FORTRAN 77 V4.7
C     ------------
C
C     Notes:
C     ------
C
C     (1)  IMPLICIT NONE and eight character symbols are not (yet) standard.
C
C     (2)  If derivatives are not desired, use CSEVAL instead. Or, take
C          advantage of the fact that the returned values are computed in
C          the order  SPP, SP, S so that if only first derivative is needed,
C          for example, the same array can be passed for SPP as for SP.
C
C     Bibliography:
C     -------------
C
C     (1)  Forsythe, G., M. Malcolm, and C. Moler.  Computer Methods for
C             Mathematical Computations.  Englewood Cliffs: Prentice-Hall,
C             1977.  (Chapter 4)
C
C     Author:  Robert Kennelly, Informatics General Corp.
C     -------
C
C     History:
C     --------
C
C      8/27/84    RGL    Adapted from CSEVAL to provide derivatives.
C     10/20/87    RAK    Abstract the interval search (with revisions) as
C                        a separate module.
C     06/09/91    DAS    Ensured S (U) = Y (NX) if U = X (NX).  (This
C                        is not guaranteed by CSDVAL's use of the cubic
C                        in the (NX-1)th interval for such a U, except
C                        if the arithmetic is exact, which it isn't.)
C
C----------------------------------------------------------------------

C     Declarations.
C     -------------

      IMPLICIT NONE

C     Constants.

      REAL
     &   ONE, TWO, THREE, SIX
      PARAMETER
     &  (ONE   = 1.0E+0,
     &   TWO   = 2.0E+0,
     &   THREE = 3.0E+0,
     &   SIX   = 6.0E+0)

C     Arguments.

      INTEGER
     &   NU, NX
      REAL
     &   B (NX), C (NX), D (NX), S (NU), SP (NU), SPP (NU), U (NU),
     &   X (NX), Y (NX)

C     Local variables.

      INTEGER
     &   IU, LEFT
      REAL
     &   ARROW, DX, XRIGHT

C     Execution.
C     ----------

      ARROW = SIGN (ONE, X (2) - X (1))
      XRIGHT = X (NX)
      LEFT = 1

      DO 10, IU = 1, NU

C        Search for the best "left-hand" endpoint to use for interpolation.

         CALL INTERVAL (NX, X, U (IU), ARROW, LEFT)

C        Evaluate the spline.  Note that if U is off-scale on the left,
C        the coefficients of the first interval are used for extrapolation,
C        while if U is greater than X(NX), the (NX-1)th set is used.
C        This means when U = X (NX), roundoff could mean S (U) is not
C        exactly Y (NX).  Hence the test for this special case.
C        Note also the order, in case y" is not really needed - see NOTES.

         DX = U (IU) - X (LEFT)

         SPP (IU) = TWO * C (LEFT) +
     &       (DX * (SIX * D (LEFT)))

         SP  (IU) = B (LEFT) +
     &       (DX * (TWO * C (LEFT) +
     &       (DX * (THREE * D (LEFT)))))

         S   (IU) = Y (LEFT) +
     &       (DX * (B (LEFT) +
     &       (DX * (C (LEFT) +
     &       (DX * (D (LEFT)))))))

         IF (U (IU) .EQ. XRIGHT) S (IU) = Y (NX)
   10 CONTINUE

C     Termination.
C     ------------

      END SUBROUTINE CSDVAL

C+---------------------------------------------------------------------
C
      SUBROUTINE CSEVAL (NDATA, X, Y, NU, U, B, C, D, S)
C
C     One-liner: Cubic spline evaluation on an array.
C     ----------
C
C     Description and usage:
C     ----------------------
C
C        CSEVAL evaluates the 1-dimensional cubic spline function
C
C     S(U) = Y(I) + B(I)*(U - X(I)) + C(I)*(U - X(I))**2 + D(I)*(U - X(I))**3
C
C     using Horner's Rule for polynomials. Normally, if U < X(1), then I = 1
C     is used, and if U > X(NX) then I = NX - 1 is used.  However, this version
C     treats periodic data properly - see the usage of NDATA.
C
C        The data must be monotonic, but may be increasing or decreasing.
C     (The coefficients associated with the NXth point are not required.)
C     Most of the effort is in deciding which set of spline coefficients to
C     use.  An "interpolation" search is employed which should be highly
C     efficient if the knots are fairly uniformly spaced (commonly the case).
C
C        CSEVAL was developed as a companion routine to CSFIT, but it should
C     work with any method of generating cubic spline coefficients.  CSDVAL
C     should be used if derivatives of the spline are required.
C
C     Arguments:
C     ----------
C
C     Name     Dimension  Type  I/O/S  Description
C     NDATA                I    I      NX = |NDATA| is the number of data points
C                                      defining the spline ("knots");
C                                      use NDATA = -NX to signify the periodic
C                                      data case; NX >= 2.
C
C     X        NX          R    I      Abscissas of knots. Must be monotone
C                                      increasing or decreasing (consistent
C                                      with the data used to fit the spline).
C
C     Y        NX          R    I      Ordinates of knots.
C
C     NU                   I    I      Number of points at which to evaluate
C                                      the spline.
C
C     U        NU          R    I      Abscissas at which to evaluate spline.
C
C     B,C,D    NX          R    I      Spline coefficients, e.g. as computed
C                                      by CSFIT.
C
C     S       NU           R      O    Spline values at points U(I).
C
C     External modules:
C     -----------------
C
C        INTERVAL   Interpolation search for interval containing a point.
C
C     Environment:  FORTRAN 77 + IMPLICIT NONE
C     ------------
C
C     Bibliography:
C     -------------
C
C     (1)  Forsythe, G., M. Malcolm, and C. Moler.  Computer Methods for
C             Mathematical Computations.  Englewood Cliffs: Prentice-Hall,
C             1977.  (Chapter 4)
C
C     Author:  Robert Kennelly, Informatics General Corp.
C     -------
C
C     History:
C     --------
C
C      7/27/84    RAK    Original design and coding, based in part on
C                        SEVAL from Forsythe, Malcolm, and Moler, and
C                        using a search method adapted from Sedgewick.
C                        An earlier rewrite of SEVAL by Trosin included
C                        conversion from a FUNCTION to a SUBROUTINE to
C                        allow more than one evaluation per call.
C     10/19/87    RAK    Abstracted the interval search (with revisions)
C                        as a separate module.
C     06/09/91    DAS    Ensured S (U) = Y (NX) if U = X (NX).  (This
C                        is not guaranteed by CSEVAL's use of the cubic
C                        in the (NX-1)th interval for such a U, except
C                        if the arithmetic is exact, which it isn't.)
C     08/16/97     "     Handled the periodic case via NX < 0 kludge.
C
C----------------------------------------------------------------------

C     Declarations.
C     -------------

      IMPLICIT NONE

C     Constants.

      REAL
     &   ONE
      PARAMETER
     &  (ONE = 1.0E+0)

C     Arguments.

      INTEGER
     &   NDATA, NU
      REAL
     &   B (*), C (*), D (*), S (NU), U (NU), X (*), Y (*)

C     Local variables.

      INTEGER
     &   IU, LEFT, NX
      REAL
     &   ARROW, DX, PERIOD, XEVAL, XLEFT, XRIGHT

C     Execution.
C     ----------

      ARROW  = SIGN (ONE, X (2) - X (1))
      NX     = ABS (NDATA)
      XRIGHT = X (NX)
      LEFT   = 1

      IF (NDATA .GT. 0) THEN  ! Non-periodic case

         DO IU = 1, NU

            XEVAL = U (IU)

C           Search for the best "left-hand" endpoint to use for interpolation.

            CALL INTERVAL (NX, X, XEVAL, ARROW, LEFT)

C           Evaluate the spline.  Note that if U is off-scale on the left,
C           the coefficients of the first interval are used for extrapolation,
C           while if U is greater than X(NX), the (NX-1)th set is used.
C           This means when U = X (NX), roundoff could mean S (U) is not
C           exactly Y (NX).  Hence the test for this special case.

            DX = XEVAL - X (LEFT)

            S (IU) = Y (LEFT) +
     &        (DX * (B (LEFT) +
     &        (DX * (C (LEFT) +
     &        (DX * (D (LEFT)))))))

            IF (XEVAL .EQ. XRIGHT) S (IU) = Y (NX)

         END DO

      ELSE  ! Periodic case

         IF (ARROW .EQ. ONE) THEN
            XLEFT = X (1)
         ELSE
            XLEFT = XRIGHT
            XRIGHT = X (1)
         END IF

         PERIOD = XRIGHT - XLEFT  ! Definitely positive

         DO IU = 1, NU

            XEVAL = U (IU)

            IF (XEVAL .LT. XLEFT) THEN

               XEVAL = XRIGHT - MOD (XLEFT - XEVAL, PERIOD)

            ELSE IF (XEVAL .GT. XRIGHT) THEN

               XEVAL = XLEFT + MOD (XEVAL - XRIGHT, PERIOD)

            END IF

            CALL INTERVAL (NX, X, XEVAL, ARROW, LEFT)

            DX = XEVAL - X (LEFT)

            S (IU) = Y (LEFT) +
     &        (DX * (B (LEFT) +
     &        (DX * (C (LEFT) +
     &        (DX * (D (LEFT)))))))

            IF (XEVAL .EQ. X (NX)) S (IU) = Y (NX)

         END DO

      END IF

C     Termination.
C     ------------

      END SUBROUTINE CSEVAL

C+----------------------------------------------------------------------
C
      SUBROUTINE CSFIT ( N, X, Y, IENDL, DERIVL, IENDR, DERIVR,
     +                   B, C, D, IER )
C
C  ACRONYM:  Cubic Spline FIT
C            -     -      ---
C  PURPOSE:
C    CSFIT computes the coefficients B(I), C(I), and D(I), I=1,2,...,N
C    for the conventional cubic spline defined by
C
C    S(X) = Y(I) + B(I)*(X-X(I)) + C(I)*(X-X(I))**2 + D(I)*(X-X(I))**3
C
C    for X(I) <= X <= X(I+1).  The spline interpolates the data points
C    (X(I), Y(I)) (also called the knots) and has continuous first and
C    second derivatives.
C
C    This implementation contains a number of enhancements to the sub-
C    routine SPLINE of Forsythe, Malcolm, and Moler (ref. below):
C
C       1. Several boundary conditions are provided, including period-
C          icity and the familiar "natural" spline (S"(X)=0 at the end
C          points). Different conditions may be specified at each end.
C          
C       2. Degenerate cases (N=2 and 3) are handled as far as possible
C          depending on the end-conditions specified - all of them are
C          applicable if N is at least 4. (A lot of algebra was needed
C          to derive the code for these cases.)
C
C       3. An Nth set of coefficients is computed,  corresponding to X
C          beyond X(N).  (This was needed by one version of CSEVAL and
C          is retained because some applications expect it.)
C
C       4. The search in CSEVAL now permits decreasing abscissas;  the
C          algebra in CSFIT applies without change in this case.
C
C  METHOD:
C    (0)  Represent the spline in the ith interval by
C
C         s(x) = w*y(i+1) + (1 - w)*y(i) + h(i)**2 *
C                [(w**3 - w)*sigma(i+1) + ((1 - w)**3 - w)*sigma(i)],
C
C         where h(i) = x(i+1) - x(i)  and  w = (x - x(i)) / h(i).
C
C         (as on p.71 of Forsythe, Malcolm, and Moler).
C    (1)  Check for and process any degenerate case of the cubic spline.
C         (See description for N under ARGUMENTS.)
C    (2)  Construct the bulk of the symmetric tridiagonal system
C         involved (common to all cases).
C    (3)  Apply the indicated end conditions for this system.
C    (4)  Solve the system for the sigmas.
C    (5)  Derive the corresponding polynomial coefficients.
C
C  ARGUMENTS:
C    NAME   DIM TYPE I/O/S DESCRIPTION
C    N       -   I     I   Number of data points, or knots.  Minimum N
C                          depends on the spline to be constructed--
C                          > 3: All possible end conditions are feasible;
C                          = 3: If IENDL=2, then IENDR must also be 2;
C                          = 2: The same order derivative must be used
C                          at both ends in all cases; if third deriva-
C                          tive is used, the same value also must be
C                          applied at each end.
C    X       N   R     I   Abscissas of the knots, strictly increasing
C                          or strictly decreasing.
C    Y       N   R     I   Ordinates of the knots.  Y(N)=Y(1) for the
C                          cyclic case (IENDL=4).
C    IENDL   -   I     I   = 0: The 3rd derivative of the spline at the
C                               left endpoint is to match the 3rd deriv-
C                               ative of the cubic passing through the
C                               first 4 data points; if N=3 and IENDL=
C                               IENDR=0, both parts of the spline become
C                               the same quadratic defined by the pts.;
C                          = 1: The 1st derivative of the spline at the
C                               left endpoint is to be the given DERIVL;
C                          = 2: The 2nd derivative of the spline at the
C                               left endpoint is to be the given DERIVL;
C                          = 3: The 3rd derivative of the spline at the
C                               left endpoint is to be the given DERIVL;
C                          = 4: The spline and its first 3 derivatives
C                               at X(N) are to match the values at X(1).
C                               This is the cyclic, or periodic, case.
C    DERIVL  -   R     I   Value of derivative used if IENDL=1, 2, or 3;
C                          ignored if IENDL=0 or 4.
C    IENDR,  -             As for IENDL, DERIVL, but pertaining to the
C    DERIVR                right endpoint; ignored if IENDL=4.
C    B,C,D   N   R     O   Spline coefficients (see PURPOSE and NOTES).
C    IER     -   I     O   =0: No errors were detected.
C                          =1: Too few data points; N < 2.
C                          =2: Degenerate cases of N=2 for all end con-
C                              ditions and of N=3 with second derivative
C                              end condition require that the same order
C                              of derivative is applied at both ends.
C                          =3: Degenerate case where N=2 with third der-
C                              ivative applied must have the same value
C                              for the derivative at each endpoint.
C                          =4: Cyclic mode -- Y(N) not equal to Y(1).
C                          =5: Non-cyclic mode -- IENDL or IENDR is
C                              beyond the range [0,4].
C
C  EXTERNAL REFERENCES:
C    NAME     DESCRIPTION
C    TRICPS   Solves a cyclic, positive-definite, symmetric tridiagonal
C             system of equations.  Also used here for true tridiagonal
C             cases.
C
C  SYSTEM DEPENDENCIES:
C    (1)  IMPLICIT NONE is an extension of FORTRAN 77.
C
C  NOTES:
C    (1)  Y(I) = S(X(I))
C         B(I) = S'(X(I))
C         C(I) = S''(X(I))/2
C         D(I) = S'''(X(I))/6
C
C    (2)  To evaluate the spline, use subroutine CSEVAL.
C
C  ENVIRONMENT:  VAX, CRAY -- FORTRAN 77
C
C  BIBLIOGRAPHY:
C    Forsythe, Malcolm, and Moler, Computer Methods for Mathematical
C       Computations. Englewood Cliffs, NJ:  Prentice-Hall, 1977.
C
C  DEVELOPMENT HISTORY:
C    DATE   INITIALS   DESCRIPTION 
C   c. 1977  F,M,M     Original implementation (IENDL=IENDR=0 case only)
C   06/25/84 RGL/RCL/  Now provides for multiple end conditions, includ-
C            RAK/DAS   ing cyclic functions.
C   10/23/87 DAS       Revised PURPOSE description (prompted by CSEVAL's
C                      new ability to deal with descending abscissas).
C
C  AUTHORS:  Forsythe, Malcolm, and Moler, c. 1977
C            Informatics, Palo Alto, CA., June 1984
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C ... Arguments:

      INTEGER  N, IENDL, IENDR, IER
      REAL     X(N), Y(N), B(N), C(N), D(N)

C ... Local constants:

      REAL      ZERO, ONE, TWO, THREE, FOUR, FIVE, SIX
      PARAMETER ( ZERO=0.E+0, ONE=1.E+0, TWO=2.E+0, THREE=3.E+0,
     +            FOUR=4.E+0, FIVE=5.E+0, SIX=6.E+0 )

C ... Local variables:

      INTEGER  I, ISYS, NSYS
      REAL     DERIVL, DERIVR, DYDX1, DYDX2, DYDX3, H1, H2, H3

C ... Procedures:

      EXTERNAL  TRICPS

C ... Execution:

C ... Check for too few data points:

      IER = 0
      IF ( N.LT.2 ) THEN
         IER = 1
         GO TO 999
      END IF

C ... N will commonly be more than enough, so:

      IF ( N.GT.4 ) GO TO 10

C ... Handle degenerate cases:

      IF ( N.EQ.2 ) THEN
         IF ( IENDL.NE.IENDR ) THEN
            IER = 2
            GO TO 999
         END IF

         H1    = X(2) - X(1)
         DYDX1 = (Y(2) - Y(1)) / H1

         IF ( IENDL.EQ.0 ) THEN
            B(1) = DYDX1
            C(1) = ZERO
            D(1) = ZERO
         ELSE IF ( IENDL.EQ.1 ) THEN
            B(1) = DERIVL
            C(1) = (DYDX1 * THREE - DERIVL * TWO - DERIVR) / H1
            D(1) = (DERIVL + DERIVR - DYDX1 * TWO) / H1**2
         ELSE IF ( IENDL.EQ.2 ) THEN
            B(1) = DYDX1 - (DERIVL * FOUR - DERIVR) * (H1 / SIX)
            C(1) = DERIVL / TWO
            D(1) = (DERIVR - DERIVL) / (H1 * SIX)
         ELSE IF ( IENDL.EQ.3 ) THEN
            IF ( DERIVL.EQ.DERIVR ) THEN
               B(1) = (Y(2) - Y(1)) / H1 - DERIVL / SIX * H1**2
               C(1) = ZERO
               D(1) = DERIVL / SIX
            ELSE
               IER = 3
               GO TO 999
            END IF
         ELSE IF ( IENDL.EQ.4 ) THEN
            IF ( Y(N).NE.Y(1) ) THEN
               IER = 4
               GO TO 999
            END IF
            B(1) = ZERO
            C(1) = ZERO
            D(1) = ZERO
         END IF
         GO TO 50

      ELSE IF ( N.EQ.3 ) THEN
         H1  = X(2) - X(1)
         H2  = X(3) - X(2)
         DYDX1 = (Y(2) - Y(1)) / H1
         DYDX2 = (Y(3) - Y(2)) / H2

         IF ( IENDL.EQ.2 .OR. IENDR.EQ.2 ) THEN
            IF ( IENDL.NE.IENDR ) THEN
               IER = 2
               GO TO 999
            END IF

            C(1) = DERIVL / SIX
            C(2) = (DYDX2 - DYDX1 - (DERIVL * H1 + DERIVR * H2) /
     +             SIX) / ((H1 + H2) * TWO)
            C(3) = DERIVR / SIX
            GO TO 30
         ELSE IF ( IENDL.EQ.4 ) THEN
            IF ( Y(N).NE.Y(1) ) THEN
               IER = 4
               GO TO 999
            END IF
            C(1) = (DYDX1 - DYDX2) / (H1 + H2)
            C(2) =-C(1)
            C(3) = C(1)
            GO TO 30
         END IF

      ELSE IF ( N.EQ.4 ) THEN
         IF ( IENDL.EQ.2 .AND. IENDR.EQ.2 ) THEN
            H1  = X(2) - X(1)
            H2  = X(3) - X(2)
            H3  = X(4) - X(3)
            DYDX1 = (Y(2) - Y(1)) / H1
            DYDX2 = (Y(3) - Y(2)) / H2
            DYDX3 = (Y(4) - Y(3)) / H3

            C(1)  = DERIVL / SIX
            C(2)  = (DYDX3 - DYDX2 - DERIVR * H3 / SIX -
     +              (H3 / H2 + ONE) * (DYDX2 - DYDX1 -
     +              DERIVL * H1 / SIX) * TWO) /
     +              (H2 - (H2 + H3) *
     +              (H1 / H2 + ONE) * FOUR)
            C(3)  = (DYDX2 - DYDX1 - DERIVL * H1 / SIX -
     +              (H1 / H2 + ONE) * (DYDX3 - DYDX2 -
     +              DERIVR * H3 / SIX) * TWO) /
     +              (H2 - (H2 + H3) *
     +              (H1 / H2 + ONE) * FOUR)
            C(4)  = DERIVR / SIX
            GO TO 30
         END IF
      END IF

   10 CONTINUE

C ... Set up the bulk of the tridiagonal system (common to all cases).
C     B = diagonal, D = offdiagonal, C = right hand side.

      D(1) = X(2) - X(1)

      DO 20 I = 2, N-1
         D(I) = X(I+1) - X(I)
         B(I) = ( X(I+1) - X(I-1) ) * TWO
         C(I) = ( Y(I+1) - Y(I) ) / D(I) -
     +          ( Y(I) - Y(I-1) ) / ( X(I) - X(I-1) )
   20 CONTINUE

C ... Now for the boundary conditions.  Cyclic case ignores IENDR.

      ISYS = 1
      IF ( IENDL.NE.4 ) THEN

C ...    The two ends of the spline are independent.

         NSYS = N
         D(N) = ZERO

C ...    Set the left end condition:

         IF ( IENDL.EQ.0 ) THEN

C ...       Use divided differences for the 3rd derivative of the
C           cubic through the first 4 points:

            B(1) = -D(1)
            C(1) = ZERO
            IF ( N.GT.3 )
     +         C(1) = (C(3)/(X(4)-X(2)) - C(2)/(X(3)-X(1)))*
     +                D(1)**2/(X(4)-X(1))
         ELSE IF ( IENDL.EQ.1 ) THEN
            B(1) = D(1) * TWO
            C(1) = - ( DERIVL - ( Y(2) - Y(1) )/D(1) )
         ELSE IF ( IENDL.EQ.2 ) THEN
            NSYS = NSYS - 1
            ISYS = 2
            C(1) = DERIVL/SIX
            C(2) = C(2) - D(1) * C(1)
         ELSE IF ( IENDL.EQ.3 ) THEN
            B(1) = -D(1)
            C(1) = ( DERIVL/SIX ) * D(1)**2
         ELSE
            IER  = 5
            GO TO 999
         END IF

C ...    Set the right end condition similarly:

         IF ( IENDR.EQ.0 ) THEN
            B(N) = -D(N-1)
            C(N) = ZERO
            IF ( N.GT.3 )
     +         C(N) = (C(N-2)/(X(N-1)-X(N-3)) - C(N-1)/(X(N)-X(N-2)))*
     +                D(N-1)**2/(X(N)-X(N-3))
         ELSE IF ( IENDR.EQ.1 ) THEN
            B(N) = D(N-1) * TWO
            C(N) = DERIVR - (Y(N)-Y(N-1))/D(N-1)
         ELSE IF ( IENDR.EQ.2 ) THEN
            NSYS   = NSYS - 1
            C(N)   = DERIVR/SIX
            C(N-1) = C(N-1) - D(N-1) * C(N)
            D(N-1) = ZERO
         ELSE IF ( IENDR.EQ.3 ) THEN
            B(N) = -D(N-1)
            C(N) = -(DERIVR/SIX) * D(N-1)**2
         ELSE
            IER  = 5
            GO TO 999
         END IF
      ELSE

C ...    Cyclic boundary conditions, IENDL = 4.  Sigma(N) = sigma(1) is
C        used to keep the system symmetric, giving one equation less:

         IF ( Y(N).NE.Y(1) ) THEN
            IER = 4
            GO TO 999
         END IF

         NSYS = N-1
         B(1) = (D(1) + D(NSYS)) * TWO
         C(1) = (Y(2) - Y(1)) / D(1) - (Y(N) - Y(NSYS)) / D(NSYS)

      END IF

C ... Solve the tridiagonal system for the sigmas (returned in C(*)):

      CALL TRICPS ( NSYS, B(ISYS), D(ISYS), C(ISYS), C(ISYS) )

      IF ( IENDL.EQ.4 ) C(N) = C(1)

C ... Derive the spline coefficients (noting that the DXs are lost):

   30 DO 40 I = 1, N-1
         H1   = X(I+1) - X(I)
         B(I) = (Y(I+1) - Y(I)) / H1 - (C(I+1) + C(I)*TWO) * H1
         D(I) = (C(I+1) - C(I)) / H1
         C(I) = C(I) * THREE
   40 CONTINUE

C ... Set coefficients for rightmost knot (and beyond).

   50 IF ( IENDL.NE.4 ) THEN

C ...    1st deriv. at X(N) uses right-sided deriv. for X in X(N-1):X(N)

         H1   = X(N) - X(N-1)
         B(N) = B(N-1) + C(N-1) * H1 * TWO + D(N-1) * H1**2 * THREE
         C(N) = C(N-1) + D(N-1) * H1 * THREE
         D(N) = D(N-1)
      ELSE

C ...    Use periodicity for IENDL=4 case:

         B(N) = B(1)
         C(N) = C(1)
         D(N) = D(1)
      END IF

  999 RETURN

      END SUBROUTINE CSFIT

C+----------------------------------------------------------------------
C
      SUBROUTINE HTDIS4 (DSINPUT, XA, XB, D1, D2, N, X, LUNOUT, IER)
C
C  ACRONYM:  Hyperbolic Tangent-type DIStribution, 2-sided, 4th version
C            -          -            ---                    -
C  PURPOSE
C  -------
C
C     HTDIS4 is the version of HTDIS2 intended for 64-bit systems.
C     It differs only in its internal procedure, HTDIS3, which need
C     not use DOUBLE PRECISION the way it does for HTDIS2.
C
C     The remaining documentation is that of HTDIS2.
C
C     HTDIS2 generates abscissas on the interval [XA, XB] using the
C     method of Marcel Vinokur to achieve asymmetric bunching.  The
C     input D1 and D2 may be either the desired initial and final
C     INCREMENTS ("deltas") if DSINPUT is .TRUE., or the desired
C     end-point SLOPES of the (normalized) stretching function if
C     DSINPUT is .FALSE.
C
C  METHOD
C  ------
C
C     The mathematics appear in the paper
C
C        "On One-Dimensional Stretching Functions Functions for Finite-
C        Difference Calculations" by Marcel Vinokur,
C        Journal of Computational Physics, Vol. 50, No. 2, May 1983.
C
C     The paper develops criteria for stretching functions that optimize
C     in some sense the truncation errors inherent in finite-difference
C     approximations to derivatives.  The simplest "universal function"
C     (of which a scaled portion is employed) satisfying the criteria is
C     w = tan z, where  z = x + iy is complex.
C
C     The analysis uses the quantity B = SQRT (S0 * S1), where S0 and S1
C     are dimensionless end-point SLOPES of the stretching function.
C     I.e.,  S0 = dXI / dT  at  T = 0,  and  S1 = dXI / dT  at  T = 1,
C     where XI and T are normalized variables and XI = XI (T) is the
C     stretching function.
C
C     For the cases of usual interest, B > 1 and z is purely imaginary,
C     leading to relationships involving sinh and tanh - hence the name
C     HTDIS2.  However, for B < 1, z is real and the analogous solution
C     involves sine and tangent.  Furthermore, for B ~ 1, both formula-
C     tions break down, so a third relation is used involving B - 1.
C
C     In this implementation, a Newton iteration is used to solve
C
C        SINH (DEL Y) / (DEL Y) = B   if  B > 1, or
C
C        SIN (DEL X) / (DEL X)  = B   if  B < 1.
C
C     Then for each XI = (I - 1) / (N - 1) an ANTISYMMETRIC function is
C     calculated:
C
C        U = 0.5 + 0.5 TANH [DEL Y * (XI - 0.5)] / TANH (0.5 * DEL Y) or
C
C        U = 0.5 + 0.5 TAN  [DEL X * (XI - 0.5)] / TAN  (0.5 * DEL X)
C
C     from which the unsymmetric general case is derived as
C
C        T = U / [A + (1 - A) * U]    where  A = SQRT (S0 / S1).
C
C     For B ~ 1. (actually | B - 1. | < 1.E-5), the approximation is
C
C        U = XI * [1. + 2. * (B - 1.) * (XI - 0.5) * (1. - XI)]
C
C     with no nonlinear equation to solve.  B = 1 if the distribution
C     is uniform, and this is handled as a special case if D1 = D2.
C
C  CLARIFICATION OF WHAT D1 AND D2 REALLY MEAN
C  -------------------------------------------
C
C     Confusion has arisen from the tendency to specify first and last
C     INTERVALS (D1 and D2), while the method calls for specifying SLOPES
C     at the stretching function's end-points.  In fact, one is related
C     to the reciprocal of the other.  We have, for the end-point slope,
C
C        Slope ~ [XI(2) - XI(1)] / [T(2) - T(1)]
C              = [(2 - 1) / (N - 1)] / dT(1)
C              = [1 / (N - 1)] * [1 / dT(1)]
C
C     Vinokur's use of the end slopes was intended to enable blending
C     of two distributions end-on-end smoothly.  However, there may be
C     situations where specifying and achieving precise increments is
C     preferred (for instance if it is desired to prepare a table
C     showing the behavior of the distribution in terms of largest and
C     smallest increments for varying N or for fixed N and varying D1, D2).
C
C     Therefore, this version has an option to perform an outer iteration
C     (via the secant method) in order to achieve specified increments
C     more precisely.  DSINPUT = .TRUE. invokes the outer iteration.
C
C     [This outer iteration really has to solve a nonlinear equation
C     in TWO variables - the D1, D2 arguments to lower-level routine
C     HTDIS3.  But independent parallel iterations with the secant
C     (non-derivative) method - one for adjusting D1 and one for
C     adjusting D2 - appear to be reliable in this peculiar circumstance.
C     Extreme cases can cause these iterations to converge before
C     the normal tolerances are met, but the results should be close
C     anyway.  IER = 3 in these cases.]
C
C  ARGUMENTS
C  ---------
C
C   NAME   DIM   TYPE I/O/S DESCRIPTION
C
C   DSINPUT -     L     I   .TRUE. means D1 and D2 are the desired first
C                                  and last increments;
C                           .FALSE. means they are the desired end-point
C                                  slopes of the normalized stretching
C                                  function.
C   XA      -     R     I   Desired X(1); passing X(1) here is safe.
C   XB      -     R     I   Desired X(N); passing X(N) here is safe.
C                           Note that XA = 0. and XB = 1. can provide one
C                           normalized distribution for multiple reuse,
C                           in some cases.  XA > XB is required.
C   D1,     -     R     I   See the DSINPUT description.  D1 and D2 should
C   D2                      be positive in either mode.  If they represent
C                           slopes, they refer to the end-slopes of the
C                           curve formed by plotting (I - 1) / (N - 1)
C                           (vertical axis) vs. (X (I) - XA) / (XB - XA)
C                           (horizontal axis).  This curve is independent
C                           of N for given slopes.  A larger slope means
C                           a higher density of points at the corresponding
C                           end, in normalized space.
C   N       -     I     I   Number of points; N > 3.
C   X       N     R     O   Reqd. distribn., with X(1) = XA, X(N) = XB.
C   LUNOUT  -     I     I   LUNOUT > 0 shows the iteration history;
C                           LUNOUT < 0 suppresses it (but any error
C                                      message goes to |LUNOUT|).
C   IER     -     I     O   IER = 0 means no problems;
C                                 1 means bad combination of D1 and D2,
C                                   such as D1 + D2 > XB - XA;
C                                 2 means the inner Newton iteration
C                                   failed - must have been bad inputs.
C                                   X (*) is not usable in this case;
C                                 3 means the outer iteration by the
C                                   secant method did not meet the
C                                   normal tolerances, presumably
C                                   because the case was extreme,
C                                   but things stopped changing so
C                                   it did the best it could.
C                                   X (*) should still be usable;
C                                 4 means the outer iteration did not
C                                   converge after 20 passes.  X (*)
C                                   is not usable in this case.
C  ERROR HANDLING:
C
C     See LUNOUT and IER.  Failure handling is left to the calling program,
C     which can proceed or reprompt for parameters and try again.
C
C  INTERNAL PROCEDURE:  HTDIS3 (the original HTDIS2 with additional arguments
C                       for performing an efficient outer iteration).
C
C  ENVIRONMENT:  FORTRAN 90
C
C  HISTORY:
C  c. 1980   M.Vinokur   Original analysis.
C  04/01/89  J.E.Melton  Initial implementation of HTDIS2 (bisection
C                        used for the inner iteration).
C  05/01/89  R.G.Langhi  Introduced error handling; patterned after
C                        one-sided routine HTDIS.
C  08/10/89  D.Saunders  Internal arithmetic is now in DOUBLE PRECISION
C                        to help application to Navier-Stokes grids.
C                        Results prove to be very similar, with the
C                        bisection just less likely to fail.  It seems
C                        great precision in DEL does not help:  X(2) and
C                        X(N-1) are still strangely imprecise except for
C                        the uniform case (where DEL in SINGLE and DOUBLE
C                        are quite different but give similar X(I)s).
C  11/30/90  D.Saunders  Safeguarded the B <= 1 cases, which have no soln.
C  04/23/91   "    "     Introduced an outer iteration for precise D1, D2:
C                        pushed the main algorithm down a level with added
C                        arguments for more efficient solution of the
C                        nonlinear equation.  6 steps and a fixed starting
C                        guess for DEL seem to suffice for likely cases.
C  04/28/91   "    "     Introduced a Newton inner iteration.  (Bisection
C                        hardly benefits from good starting guesses.)
C  06/04/91   "    "     Laborious explanation of the misconception that
C                        had arisen from the apparently imprecise results.
C                        Retained the more-precise option by setting the
C                        number of outer iterations to 1 or 6 according
C                        to the signs of the input D1, D2.
C  08/20/91   "    "     Incorporated the B < 1 and B ~ 1 cases.
C  08/28/91   "    "     The small-correction outer iteration was found to
C                        diverge on extreme cases with small N.  Replaced
C                        it with the secant method (up to 20 iterations
C                        for each of D1 and D2 in parallel).  IER = 3 or
C                        4 are new possibilities.
C  09/10/91   "    "     Introduced the DSINPUT argument to clarify the
C                        "deltas" or "slopes" options.
C  03/03/95   "    "     Encountered a near-singular case (D1 * D2 ~ Du**2)
C                        during iteration, which improperly terminated.
C                        EPS = 1.E-5 instead of 1.E-3 in HTDIS3 helps, and
C                        the outer iteration may now continue for CASE = 3.
C                        Improved the second estimates derived from the
C                        first call to HTDIS3.
C  03/03/95   "    "     HTDIS4 variation for 64-bit systems.
C  08/07/99   "    "     Fortran 90 upgrade: HTDIS3 is internal now, using
C                        just 2 arguments.  The other 10 become global.
C  01/23/05   "    "     One case out of nearly 28,000 hit the 30-iteration
C                        limit in HTDIS3.  It needed 36 iterations, so the
C                        limit is now 50.  Also, lowering EPS from 1.E-5 to
C                        1.E-6 gave identical results without encountering
C                        the case 3 situation that led to 36 iterations.
C
C  AUTHORS: Marcel Vinokur, NASA/Ames Research Center, Moffett Field, CA
C           and John Melton, Ron Langhi, David Saunders
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER
     >   IER, LUNOUT, N
      REAL
     >   D1, D2, XA, XB, X (N)
      LOGICAL
     >   DSINPUT

C     Local constants:

      INTEGER, PARAMETER ::
     >   MXITER = 20      ! Max. number of secant method iterations

      REAL, PARAMETER ::
     >   DLOW   = 1.E-2,  ! Smallest reduction in trial D1, D2 inputs to HTDIS3.
     >   ONE    = 1.E+0,
     >   THREE  = 3.E+0,
     >   TOL    = 1.E-4,  ! Relative tolerance on D1, D2
     >   ZERO   = 0.E+0

C     Local variables:

      INTEGER
     >   I

      REAL
     >   DEL0, DELSOLN, DMIN1, DMIN2, DX1 (3), DX2 (3), F1 (2), F2 (2),
     >   R, TOL1, TOL2, X1TARG, X2TARG

      LOGICAL
     >   CLOSE, TRIAL

C     Execution:
C     ----------

      TRIAL   = DSINPUT        ! Forces iteration to match Ds precisely
      DX1 (1) = D1
      DX2 (1) = D2

      IF (TRIAL) THEN
         TOL1   = DX1 (1) * TOL
         DMIN1  = DX1 (1) * DLOW
         X1TARG = XA + DX1 (1)      ! Target X (2)

         TOL2   = DX2 (1) * TOL
         DMIN2  = DX2 (1) * DLOW
         X2TARG = XB - DX2 (1)      ! Target X (N-1)

         IF (X1TARG > X2TARG) THEN  ! Bad D1 and D2
            IER = 1
            GO TO 99
         END IF
      END IF

      DEL0 = THREE  ! B > 1 gives DEL ~ 0.1 to 10 (increasingly non-uniform);
                    ! B < 1 gives DEL in the range [0+, PI-].
      CLOSE = .FALSE.     ! Used to detect "converged but above tolerance"

C     The first call to the basic Vinokur method is all we need if we are
C     matching slopes, else we need it to initialize the outer iteration:

C*****CALL HTDIS3 (DSINPUT, XA, XB, DX1 (1), DX2 (1), N, X, TRIAL, DEL0,
C    >             DELSOLN, LUNOUT, IER)

      CALL HTDIS3 (DX1 (1), DX2 (1)) ! Internal procedure now

      IF (IER /= 0) GO TO 99

      IF (.NOT. DSINPUT) GO TO 99    ! Slopes are matched without iteration.

      IF (DELSOLN == ZERO) GO TO 99  ! DEL = 0. signals uniform distribution
                                     ! which is specially handled.


C     Otherwise apply two secant iterations in parallel - one for D1 and
C     one for D2 - to solve what is really a nonlinear equation in TWO
C     variables.  Finish the first function evaluation:

      F1 (1) = X (2) - X1TARG      
      F2 (1) = X2TARG - X (N - 1)

C     Now for the second evaluation.  Attempt a better estimate.

      DX1 (2) = DX1 (1) ** 2 / (X (2) - XA)
      DX2 (2) = DX2 (1) ** 2 / (XB - X (N - 1))
      DEL0 = DELSOLN

      CALL HTDIS3 (DX1 (2), DX2 (2))

      IF (IER /= 0) GO TO 99

      F1 (2) = X (2) - X1TARG      
      F2 (2) = X2TARG - X (N - 1)

C     Here is the rest of the two secant iterations in parallel:

      DO I = 1, MXITER

         DEL0 = DELSOLN

         R = F1 (1) - F1 (2)
         IF (R /= ZERO) THEN
            R = F1 (1) / R
            DX1 (3) = R * DX1 (2) + (ONE - R) * DX1 (1)
            DX1 (3) = MAX (DX1 (3), DMIN1)
         ELSE
            DX1 (3) = DX1 (2)
         END IF

         R = F2 (1) - F2 (2)
         IF (R /= ZERO) THEN
            R = F2 (1) / R
            DX2 (3) = R * DX2 (2) + (ONE - R) * DX2 (1)
            DX2 (3) = MAX (DX2 (3), DMIN2)
         ELSE
            DX2 (3) = DX2 (2)
         END IF

         CALL HTDIS3 (DX1 (3), DX2 (3))

         IF (IER /= 0) GO TO 99

         IF (.NOT. TRIAL) THEN 
            IF (CLOSE) IER = 3
            GO TO 99
         END IF

C        Normal convergence test:
C        Reset TRIAL to get the full distribution on a final pass.

         F1 (1) = F1 (2)
         F1 (2) = X (2) - X1TARG      
         F2 (1) = F2 (2)
         F2 (2) = X2TARG - X (N - 1)

         IF (ABS (F1 (2)) < TOL1 .AND. ABS (F2 (2)) < TOL2) THEN
            TRIAL = .FALSE.
         ELSE IF (DX1 (3) == DX1 (2) .AND. DX2 (3) == DX2 (2)) THEN
            CLOSE = .TRUE.
            TRIAL = .FALSE.
         END IF

         DX1 (1) = DX1 (2)
         DX1 (2) = DX1 (3)
         DX2 (1) = DX2 (2)
         DX2 (2) = DX2 (3)

      END DO

      IER = 4  ! No convergence - X (*) is not usable.

   99 RETURN   ! From HTDIS4


C     Internal procedure for HTDIS4:
C     ------------------------------

      CONTAINS

!        -----------------------------------------------------------------------
!
         SUBROUTINE HTDIS3 (D1, D2)  ! The essence of the Vinokur method
!
!        -----------------------------------------------------------------------
!
!        ACRONYM:  Hyperbolic Tangent-type DIStribution, 2-sided (3rd version)
!                  -          -            ---                    -
!
!        HTDIS3 is the essence of the original HTDIS2.  It allows the specified
!        initial and final increments to be obtained with minimal work during
!        the intermediate iterations.
!
!        ADDITIONAL VARIABLES (See HTDIS4 for the others):
!
!        VAR      TYPE I/O/S DESCRIPTION
!
!        TRIAL     L     I   TRIAL = .TRUE. suppresses calculation of the
!                            interior X (*), since only X (2) and X (N-1)
!                            are being made use of at the higher level.
!        DEL0      R     I   Initial guess for the DELSOLN which solves
!                            the relevant nonlinear equation.  If B > 1,
!                            DEL is somewhere around 1. to 10. (larger
!                            for more nonuniform distributions); 1.E-10
!                            and 25. are used as bounds if necessary.
!                            If B < 1, the bounds are 1.E-10 and PI - eps.
!        DELSOLN   R     O   This solution is returned for reuse in
!                            a possible next call.
!
!        -----------------------------------------------------------------------

!        Arguments:

         REAL
     >      D1, D2

!        Local constants:

         INTEGER, PARAMETER ::
     >      MXITER = 50

         REAL, PARAMETER ::
     >      EPS = 1.E-6, HALF = 0.5E+0, ONE = 1.E+0,
     >      PIMINUS = 3.14159E+0, TOL = 1.E-12, TWO = 2.E+0,
     >      ZERO = 0.E+0

!        Local variables:

         INTEGER
     >      CASE, I, INC, ITER, LUNMSG

         REAL
     >      A, B, DEL, DELHI, DELLO, DELINV, DUNIFM, EXPDEL, FACTOR,
     >      FPRIME, FRACT, FUNDEL, RANGE, RFNM1, SLOPE0, SLOPE1, U, XI


!        Execution:
!        ----------

         IER = 0

!        Normalize the input control parameters:

         RANGE  = XB - XA
         RFNM1  = ONE / REAL (N - 1)
         DUNIFM = RANGE * RFNM1
         SLOPE0 = D1
         SLOPE1 = D2

         IF (DSINPUT) THEN
            SLOPE0 = DUNIFM / SLOPE0
            SLOPE1 = DUNIFM / SLOPE1
         END IF

         A = SQRT (SLOPE0 / SLOPE1)
         B = SQRT (SLOPE0 * SLOPE1)

!        B distinguishes the cases.

         IF (B == ONE) THEN

            IF (SLOPE0 == SLOPE1) THEN  ! D1 = D2 = uniform DX
               DEL = ZERO               ! Signal no need to iterate
               DO I = 1, N - 1
                  X (I) = XA + REAL (I-1) * DUNIFM
               END DO
               GO TO 40
            ELSE
               CASE = 3                 ! Singular case, but not uniform
               DEL = -ONE               ! Something other than zero
            END IF

         ELSE IF (B > ONE + EPS) THEN   ! Need to solve SINH (DEL) / DEL = B

            CASE = 1
            DELHI = 25.E+0
            DELLO = 1.E-10

         ELSE IF (B < ONE - EPS) THEN   ! Need to solve SIN (DEL) / DEL = B

            CASE = 2
            DELHI = PIMINUS
            DELLO = EPS

         ELSE                           ! Simpler approximation

            CASE = 3
            DEL = -ONE

         END IF


         IF (CASE /= 3) THEN

!           Newton iteration for the two main cases:

            IF (LUNOUT > 0) WRITE (LUNOUT, 1001) CASE

            DEL = DEL0

            DO ITER = 1, MXITER

               IF (CASE == 1) THEN  ! B > 1 + EPS

!                 Solve   SINH (DEL) / DEL = B  in the form
!
!                 F (DEL) = SINH (DEL) / (DEL * B) - 1 = 0  (better scaled?)
!                 The exponential form is used since EXP (DEL) and its
!                 reciprocal should be more efficient than SINH and COSH.

                  EXPDEL = EXP (DEL)
                  DELINV = ONE / DEL
                  FRACT  = HALF * DELINV / B
                  FUNDEL = FRACT * (EXPDEL - ONE / EXPDEL) - ONE
                  FPRIME = FRACT * ((ONE - DELINV) * EXPDEL +
     >                              (ONE + DELINV) / EXPDEL)

               ELSE  ! CASE = 2; B < 1 - EPS
         
!                 Solve   SIN (DEL) / DEL = B  in the form
!
!                 F (DEL) = SIN (DEL) / (DEL * B) - 1 = 0  (better scaled?)

                  FUNDEL = SIN (DEL) / (DEL * B) - ONE
                  FPRIME = (COS (DEL) / B - FUNDEL - ONE) / DEL

               END IF

               IF (LUNOUT > 0)
     >            WRITE (LUNOUT, 1002) ITER, DEL, FUNDEL, FPRIME, TOL

               IF (ABS (FUNDEL) < TOL) GO TO 20

               DEL = MIN (MAX (DEL - FUNDEL / FPRIME, DELLO), DELHI)

            END DO

            LUNMSG = ABS (LUNOUT)
            WRITE (LUNMSG, 1003) MXITER, DEL, FUNDEL, TOL
            IER = 2
            GO TO 99

   20       CONTINUE

         END IF

!        Calculate the antisymmetric function, from which the desired
!        stretching function is derived by another simple transformation.
!        The three cases are becoming cumbersome...

         IF (CASE == 1) THEN

            FACTOR = ONE / TANH (HALF * DEL)

         ELSE IF (CASE == 2) THEN

            FACTOR = ONE / TAN (HALF * DEL)

         ELSE

            FACTOR = TWO * (B - ONE)

         END IF

         X (1) = XA

         IF (TRIAL) THEN       ! Only the 1st & last increments are looked at
            INC = MAX (1, N - 3)
         ELSE
            INC = 1
         END IF

         DO I = 2, N - 1, INC

            XI = REAL (I-1) * RFNM1 - HALF  ! Subtracting 0.5 helps a bit

            IF (CASE == 1) THEN

               U = HALF * (ONE + TANH (DEL * XI) * FACTOR)

            ELSE IF (CASE == 2) THEN

               U = HALF * (ONE + TAN (DEL * XI) * FACTOR)

            ELSE

               U = (XI + HALF) * (ONE + FACTOR * XI * (HALF - XI))

            END IF

            FRACT = U / (A + (ONE - A) * U)
            X (I) = XA + FRACT * RANGE

         END DO

   40    X (N) = XB
         DELSOLN = DEL

         IF (LUNOUT > 0)
     >      WRITE (LUNOUT, 1004) B, X (2) - X (1), X (N) - X (N - 1)

   99    RETURN

 1001    FORMAT (/, ' ITER (Case = ', I1, ')    DEL', 7X, 'F(DEL)', 11X,
     >           'F''', 10X, 'TOL')
 1002    FORMAT (1X, I2, 1P, E20.12, 3E13.5)
 1003    FORMAT (/' *** HTDIS4:  Newton iteration failed ***',
     >           /' Number of iterations:', I4,
     >           /' DEL, F(DEL), TOL:', 1P, 3E14.6,
     >           /' No distribution computed.'/)
 1004    FORMAT (/, ' B:', 1P, E14.6, '  dX(1):', E14.6, '  dX(N-1):',
     >           E14.6)

         END SUBROUTINE HTDIS3

      END SUBROUTINE HTDIS4

C+----------------------------------------------------------------------
C
      SUBROUTINE INTERVAL (NX, X, XFIND, ARROW, LEFT)
C
C     One-liner: Interpolation search for interval containing a point.
C     ----------
C
C     Description and usage:
C     ----------------------
C
C        Written primarily for interval-based interpolations such as
C     piecewise linear or cubic spline, INTERVAL performs a search to
C     locate the best interval for evaluating the interpolant at a
C     given point. The normal case returns the "left-hand" endpoint of
C     the interval bracketing the point, but for the out-of-range cases
C     below or above the range of the knots, the interval to be used is
C     the first or last. The array of knots must be monotonic, either
C     increasing or decreasing. Diagrammatically, LEFT is returned as
C     shown below for the normal case (no extrapolation):
C
C          X (1)  ...   X (LEFT)   X (LEFT+1)   ...      X (NX)
C                                ^
C                              XFIND
C
C     And for extrapolation:
C
C                     X (LEFT = 1)  ...   X (NX)
C             ^
C           XFIND
C
C     or,
C                X (1)  ...   X (LEFT = NX-1)    X (NX)
C                                                           ^
C                                                         XFIND
C
C     If the point to be bracketed (XFIND) matches one of the knots, the
C     index of that knot is returned as LEFT, i.e., the condition for a
C     bracket of an interior point is:
C
C        X (LEFT) <= XFIND < X (LEFT+1)  if  ARROW = +1.0,  or
C        X (LEFT) >= XFIND > X (LEFT+1)  if  ARROW = -1.0.
C
C        This is a low-level routine with minimal error checking. The
C     calling program is assumed to have verified the following:
C
C     (1)  NX >= 2
C     (2)  X strictly monotonic
C     (3)  ARROW = +1.0 or -1.0
C
C     Subroutine PROTECT is available from the author for easily checking
C     conditions (2) and (3). LEFT is verified on input, but efficiency in
C     loops will benefit from passing the best estimate available, usually
C     just the result of the last call.
C
C        INTERVAL was originally written for use with CSEVAL and TABLE1.
C     The interpolation search was adapted from ideas in Sedgewick's book
C     referenced below.
C
C     Arguments:
C     ----------
C
C     Name  Dimension  Type  I/O/S  Description
C     NX                I    I      Number of points in array X; must
C                                   be >= 2 (no check performed).
C
C     X        NX       R    I      Array of points defining the set
C                                   of intervals to be examined. Only
C                                   the first NX-1 points are required.
C
C     XFIND             R    I      The point for which a bracketing
C                                   interval is sought.
C
C     ARROW             R    I      Monotonicity indicator for input
C                                   array X:
C                                     -1.0  strictly decreasing
C                                      0.0  NOT ALLOWED!
C                                     +1.0  strictly increasing
C                                   Supplied by the calling routine for
C                                   reasons of speed (not checked).
C
C     LEFT              I    I/O    Input: guessed index of left-hand
C                                   endpoint of the interval containing
C                                   the specified point.
C
C                                   Output: index of the largest array
C                                   value <= specified point (if ARROW=+1.0).
C                                   Special case for data out of range:
C                                   return left endpoint of closest interval.
C                                   Thus, LEFT = 1 for XFIND < X (2), and
C                                   LEFT = NX-1 for XFIND >= X (NX-1).
C                                   (If ARROW=-1.0, reverse the inequalities.)
C
C     Environment:  Digital VAX-11/780, VMS FORTRAN
C     ------------  Apple Macintosh, Absoft MacFORTRAN/020 v2.3
C
C     Notes:
C     ------
C
C     (1)  IMPLICIT NONE and eight character symbols are not (yet) standard.
C
C     (2)  In speed-critical applications, it might be a good idea to build
C          this algorithm in-line since it is typically called many times
C          from within a loop. Another potential speed-up is removal of the
C          ARROW multiplies, which restricts the method to increasing data.
C          So far, the simplicity of separating out the messy search details
C          and the generality of bi-directional searching have outweighed
C          the modest speed penalty incurred.
C
C     Bibliography:
C     -------------
C
C     (1) Sedgewick, R.  Algorithms.  Reading: Addison-Wesley, 1983.
C            (Chap. 14)
C
C     Author:  Robert Kennelly and David Saunders, Sterling Federal Systems
C     -------
C
C     Development history:
C     --------------------
C
C     20 Oct. 1987    RAK    Interpolation search adapted (with mods.
C                            for bidirectional search and some minor
C                            repair) from CSEVAL (RAK) and TABLE1 (DAS).
C     08 Aug. 1988    DAS    Clarified descriptions of bracketing, where
C                            the inequalities depend upon ARROW.
C
C-----------------------------------------------------------------------

C     Declarations.
C     -------------

      IMPLICIT NONE

C     Constants.

      REAL
     &   ONE
      PARAMETER
     &  (ONE = 1.0E+0)

C     Arguments.

      INTEGER
     &   LEFT, NX
      REAL
     &   ARROW, X (NX), XFIND

C     Local variables.

      INTEGER
     &   LENGTH, NXLESS1, RIGHT, TRIAL
      REAL
     &   XBYARROW

C     Execution.
C     ----------

      XBYARROW = XFIND * ARROW

C     Simplify things by disposing of two important special cases so that
C     X (LEFT) and X (RIGHT) can really bracket XFIND. As a by-product,
C     this also takes care of the NX = 2, 3 cases.

      NXLESS1 = NX - 1

      IF (XBYARROW .GE. X (NXLESS1) * ARROW) THEN
         LEFT = NXLESS1
         GO TO 990
      ELSE IF (XBYARROW .LT. X (2) * ARROW) THEN
         LEFT = 1
         GO TO 990
      END IF

C      ---------------------------------
C     |                                 |
C     |   X (2) <= XFIND < X (NX - 1)   |
C     |            - or -               |
C     |   X (2) > XFIND >= X (NX - 1)   |
C     |                                 |
C     |   NX > 3                        |
C     |                                 |
C      ---------------------------------

C     Adjust the pointers. We hope that the calling routine has provided
C     a reasonable guess (since it's probably working on an ordered array
C     of points to evaluate), but check anyway.

      LEFT = MIN (MAX (2, LEFT), NX - 2)

      IF (XBYARROW .GE. X (LEFT) * ARROW) THEN
         IF (XBYARROW .LT. X (LEFT + 1) * ARROW) THEN

C           XFIND is in the original guessed-at interval.

            GO TO 990
         ELSE

C           We'll look farther to the right. Evidently LEFT was < NX - 2.

            RIGHT = NXLESS1
            LEFT  = LEFT + 1
         END IF
      ELSE

C        Look to the left of the guess. Evidently LEFT was > 2.

         RIGHT = LEFT
         LEFT  = 2
      END IF

C      ----------------------------------
C     |                                  |
C     |   2 <= LEFT < RIGHT <= NX - 1    |
C     |                                  |
C      ----------------------------------

C     The interval length must decrease each time through - terminate
C     when the correct interval is found or when the interval length
C     cannot be decreased.

   10 CONTINUE
         LENGTH = RIGHT - LEFT
         IF (LENGTH .GT. 1) THEN

C           The trial value is a "linear" estimate of the left-hand endpoint
C           of the interval bracketing the target XFIND, with protection
C           against round-off (which can affect convergence).

            TRIAL = MIN (RIGHT - 1, LEFT + MAX (0, INT (REAL (LENGTH) *
     &         (XFIND - X (LEFT)) / (X (RIGHT) - X (LEFT)))))

C            ------------------------------------------
C           |                                          |
C           |   2 <= LEFT <= TRIAL < RIGHT <= NX - 1   |
C           |                                          |
C            ------------------------------------------

C           Adjust pointers. Increase LEFT or decrease RIGHT until done.

            IF (XBYARROW .GE. X (TRIAL + 1) * ARROW) THEN
               LEFT  = TRIAL + 1
            ELSE IF (XBYARROW .LT. X (TRIAL) * ARROW) THEN
               RIGHT = TRIAL
            ELSE

C              We're done: XFIND is in the interval [X (TRIAL), X (TRIAL+1)).

               LEFT  = TRIAL
               GO TO 990
            END IF
            GO TO 10

         END IF

C     Termination.
C     ------------

  990 CONTINUE

      END SUBROUTINE INTERVAL

C+----------------------------------------------------------------------
C
      SUBROUTINE TRICPS ( N, D, L, R, S )
C
C ACRONYM: TRIdiagonal system (Cyclic, Positive definite, Symmetric)
C          ---                 -       -                  -
C
C PURPOSE: TRICPS solves one system  A s = r,  where A is a symmetric,
C          diagonally-dominant, irreducible matrix with positive diag-
C          onal elements, non-zero sub-/super-diagonals,  and non-zero
C          lower-left/upper-right elements, and zeros elsewhere. It is
C          suited to problems involving cyclic/periodic boundary  con-
C          ditions.
C
C METHOD:  The LDL' form of the Cholesky factorization of  A  is used.
C          L  in this case is unit lower bidiagonal with its last  row
C          also filled in:
C
C                 x x     x       1         x         1 x     x
C                 x x x           x 1         x         1 x   x
C                   x x x     =     x 1         x         1 x x
C                     x x x           x 1         x         1 x
C                 x     x x       x x x x 1         x         1
C
C          It is assumed that only one right-hand side is involved for
C          the given A,  so no attempt is made to separate the factor-
C          ization from the solution using the triangular factors.
C
C ARGUMENTS:
C   ARG  DIM TYPE I/O/S DESCRIPTION
C    N    -    I    I   Order of the system (>=3)
C    D    N    R   I/S  Main diagonal elements of A
C    L    N    R   I/S  Off-diagonals of A in L(1:N-1); A(N,1) in L(N)
C    R    N    R    I   Right-hand side elements (see next item)
C    S    N    R    O   Solution vector; may be the same argument as R
C                       in which case R is destroyed.  D and L are de-
C                       stroyed in either case.
C
C ERROR HANDLING:  None. A divide by zero can mean the matrix is sing-
C                  ular but is more likely to mean unwitting errors in
C                  the call to TRICPS.
C
C ENVIRONMENT:  VAX/VMS FORTRAN 77
C
C HISTORY:
C
C   05/24/84   DAS   Initial design/coding, for a spline application.
C   05/22/86   RAK   Three divides by DI are worth replacing with one
C                    reciprocal and three multiplies in the main loop.
C
C AUTHOR: David Saunders, Informatics General, Palo Alto, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER  N
      REAL     D(N), L(N), R(N), S(N)

C     Local variables:

      INTEGER  I
      REAL     DI, DINV, LIM1

C     Execution:

C     The off-diagonals are modified in place.   The modified diagonals
C     do not have to be stored for the back-solve if the matrix factors
C     are treated as (LD)L' and not L(DL').  This opens the way for the
C     filled-in row of factor L to overwrite the original diagonals. It
C     is then not possible to include the I = N-1 case in the following
C     loop ( D(N-1) being the exceptional element ). Repeating the code
C     was considered preferable to an IF test inside the loop.

      DINV = 1.E+0 / D(1)
      D(1) = L(N) * DINV
      S(1) = R(1) * DINV
      D(N) = D(N) - L(N) * D(1)
      S(N) = R(N) - L(N) * S(1)

      DO 20 I = 2, N - 2
         LIM1   =  L(I-1)
         L(I-1) =  LIM1 * DINV
         DI     =  D(I) - LIM1 * L(I-1)
         DINV   =  1.E+0 / DI
         D(I)   = -LIM1 * D(I-1) * DINV
         S(I)   = (R(I) - LIM1 * S(I-1)) * DINV
         D(N)   =  D(N) - (DI * D(I)) * D(I)
         S(N)   =  S(N) - (DI * D(I)) * S(I)
   20 CONTINUE

      I      =  N - 1
      LIM1   =  L(I-1)
      L(I-1) =  LIM1 * DINV
      DI     =  D(I) - LIM1 * L(I-1)
      D(I)   = (L(I) - LIM1 * D(I-1)) / DI
      S(I)   = (R(I) - LIM1 * S(I-1)) / DI
      D(N)   =  D(N) - (DI * D(I)) * D(I)
      S(N)   = (S(N) - (DI * D(I)) * S(I)) / D(N)
      L(I)   =  0.E+0

      DO 30 I = N - 1, 1, -1
         S(I) = S(I) - L(I) * S(I+1) - D(I)*S(N)
   30 CONTINUE

      END SUBROUTINE TRICPS

C+----------------------------------------------------------------------
C
      SUBROUTINE VINOKUR (I1, I2, D1, D2, X, LUNOUT, IER)
C
C     VINOKUR imposes the indicated 2-sided Vinokur distribution on the
C     interior points of the grid line between I1 and I2. The end points
C     are assumed to be in place, and may be in either order: this
C     routine simplifies use of the underlying HTDIS4 for the case where
C     X(I2) is less than X(I1). D1 and D2 should be positive, full scale.
C     See HTDIS4 for further details, including LUNOUT & IER arguments.
C
C     07/07/97  DAS  2-sided variant of IC1 from SYN87.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER  I1, I2, LUNOUT, IER
      REAL     D1, D2, X(1:I2)

C     Local constants:

      REAL, PARAMETER :: ONE = 1., ZERO = 0.

C     Local variables:

      INTEGER  I
      REAL     RANGE, XI1, XI2

C     Execution:

      XI1 = X(I1)
      XI2 = X(I2)
      RANGE = XI2 - XI1

C     Normalized Vinokur distribution in ascending order:

      CALL HTDIS4 (.TRUE., ZERO, ONE, ABS (D1/RANGE), ABS (D2/RANGE),
     >             I2 - I1 + 1, X(I1), -LUNOUT, IER)

C     Denormalize:

      DO I = I1, I2
         X(I) = XI1 + X(I) * RANGE
      END DO

      X(I2) = XI2  ! Exactly

      END SUBROUTINE VINOKUR
