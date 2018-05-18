C+----------------------------------------------------------------------
C
      SUBROUTINE CSQUAD (N, X, F, OFFSET, C, FINT)
C
C
C     Description and usage:
C
C           CSQUAD (Cubic Spline QUADrature) is a specialized routine for
C        integrating a function using its cubic spline representation.  The
C        (definite) integral up to each knot is returned, with an optional
C        constant added to the value of the integral array at the first point.
C        Output array FINT may share storage with spline coefficient array C.
C
C           This version uses a compact, but tricky, form of the integral
C        taken from Forsythe, et al., which is faster than the naive
C        integration formula used in the original.
C
C
C     Arguments:
C
C        Name    Dimension  Type  I/O/S  Description
C         N                  I    I      Number of knots.
C         X        N         R    I      Abscissas (spline knots Xi).
C         F        N         R    I      Function values.
C         OFFSET             R    I      Value assigned to the first element
C                                        of FINT (will be added to all).
C         C        N         R    I      Array of coefficients of (X - Xi)**2
C                                        from spline fit (= F"/2 at X = Xi).
C         FINT     N         R      O    Integral of F up to each knot.
C                                        FINT (1) = OFFSET; FINT (N) =
C                                        OFFSET + integral from X (1) to X (N).
C
C
C     Standards violations:  IMPLICIT NONE is non-standard.
C
C
C     Bibliography:
C
C        (1)  Forsythe, G., M. Malcolm, and C. Moler.  Computer
C                Methods for Mathematical Computations (Englewood Cliffs:
C                Prentice-Hall, 1977).  Chapter 5.
C
C
C     Development environments:  Digital VAX-11/780  VMS/V4.1   FORTRAN
C                                Cray X-MP/48        COS/V1.14  CFT/V1.11
C
C
C     Author:  Robert Kennelly, Sterling Software, Palo Alto, CA.
C
C
C     History:
C
C         6 July 1984    RAK    Initial design/coding as INTEG in FLO6QNM.
C        10 July 1984    RAK    Modified to permit shared storage,
C                               added OFFSET parameter.
C        11 Feb. 1986    RAK    Use IMPLICIT NONE, and declare all
C                               constants locally.
C        20 Mar. 1986    RAK    Switch to cute FM&M form of integral.
C         9 Sep. 1986    RAK    Renamed CSQUAD for general use.
C
C-----------------------------------------------------------------------


C     Declarations.
C     -------------

      IMPLICIT NONE

C     Constants.

      REAL
     >   HALF, TWELTH
      PARAMETER
     >   (HALF   = 1.0E+0 / 2.0E+0,
     >    TWELTH = 1.0E+0 / 12.0E+0)

C     Local variables.

      INTEGER
     >   I, N
      REAL
     >   C (N), F (N), FINT (N), H, OFFSET, X (N)


C     Execution.
C     ----------

C     Compute integral on each subinterval.  Work backwards to avoid
C     overwriting the C array if it was also passed in as FINT.

      DO 10, I = N - 1, 1, -1
         H = X (I + 1) - X (I)
         FINT (I + 1) = H      * ((F (I) + F (I + 1)) * HALF -
     >                  H ** 2 *  (C (I) + C (I + 1)) * TWELTH)
   10 CONTINUE

C     Sum the contributions from each subinterval.

      FINT (1) = OFFSET
      DO 20, I = 2, N
         FINT (I) = FINT (I) + FINT (I - 1)
   20 CONTINUE


C     Termination.
C     ------------

      RETURN
      END
