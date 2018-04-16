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
      END
