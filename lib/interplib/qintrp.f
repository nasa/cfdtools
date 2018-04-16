C+----------------------------------------------------------------------
C
      SUBROUTINE QINTRP (NX, X, XL, YL, SL, XR, YR, SR, Y)
C
C  ONE-LINER:  Quadratic INTeRPolation given two points and one slope
C
C  PURPOSE:
C
C        QINTRP performs a form of quadratic interpolation, given two
C     points and one slope.  For given X(s), it returns corresponding
C     Y(s) along with the slope at the second end point (for possible
C     reuse in an adjacent interval).
C
C        The quadratic is of this form:  Y  =  A X ** 2  +  B X  +  C
C
C        QINTRP was prompted by an application involving nonlinear
C     spanwise lofting of wing sections.
C
C  ARGUMENTS:
C   ARG   DIM   TYPE I/O/S DESCRIPTION
C   NX     -      I    I   No. of interpolated elements required (>=1).
C   X      NX     R    I   Abscissa(s) at which to interpolate.
C   XL     -      R    I   Abscissa of "left" end point.
C   YL     -      R    I   Ordinate  "    "    "    "  .
C   SL     -      R    I   Slope     "    "    "    "  .
C   XR     -      R    I   Abscissa of "right" end point.
C   YR     -      R    I   Ordinate  "    "    "    "  .
C   SR     -      R    O   Slope calculated for right end point.
C   Y      NX     R    O   Interpolated value(s).
C    
C  ENVIRONMENT: FORTRAN 77
C
C  HISTORY:  10/16/92  DAS  Initial implementation.
C
C  AUTHOR:   David Saunders, Sterling Software/NASA Ames, Mt. View, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER
     >   NX
      REAL
     >   X (NX), XL, YL, SL, XR, YR, SR, Y (NX)

C     Local variables:

      INTEGER
     >   I
      REAL
     >   A, B, C, DX, DY, S

C     Execution:

C     The following expressions are obtained by solving a 3 x 3 system:

      DX = XR - XL
      DY = YR - YL
      S  = DY / DX
      SR = S + S - SL

      A  = (S - SL) / DX
      B  = (SL * (XR + XL) - S * (XL + XL)) / DX
      C  = YL - XL * (B + XL * A)

      DO I = 1, NX
         Y (I) = C + X (I) * (B + X (I) * A)
      END DO

      RETURN
      END
