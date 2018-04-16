C+----------------------------------------------------------------------
C
      SUBROUTINE CUBINTRP (NX, X, XL, YL, SL, XR, YR, SR, Y)
C
C  ONE-LINER:  CUBic INTeRPolation given two points and two slopes
C
C  PURPOSE:
C
C        CUBINTRP performs a form of cubic interpolation, given two
C     points and two slopes.  For given X(s), it returns corresponding
C     Y(s).  The relevant formula is:
C
C        Y (X)  =  A U ** 3  +  B U ** 2  +  SL U  +  YL
C
C     where  U = X - XL.  If  H = XR - XL  and  SM = (YR - YL) / H  then
C
C            A = (SL - 2 SM + SR) / H   and
C            B = (SM - SL) / H  -  A
C
C     The coefficients are NOT returned.
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
C   SR     -      R    I   Slope     "    "    "    "
C   Y      NX     R    O   Interpolated value(s).
C    
C  ENVIRONMENT: FORTRAN 77
C
C  HISTORY:  05/24/93  DAS  Adapted from QINTRP.
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

C     Local constants:

      REAL
     >   ONE
      PARAMETER
     >  (ONE = 1.E+0)

C     Local variables:

      INTEGER
     >   I
      REAL
     >   A, B, C, H, SM, U

C     Execution:

      H = ONE / (XR - XL)
      SM = (YR - YL) * H
      A = (SL - SM - SM + SR) * H
      B = (SM - SL) * H - A
      A = A * H

      DO I = 1, NX
         U = X (I) - XL
         Y (I) = YL + U * (SL + U * (B + U * A))
      END DO

      RETURN
      END
