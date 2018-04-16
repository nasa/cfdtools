C+------------------------------------------------------------------------------
C
      SUBROUTINE INVERSE_PARABOLA (NX, X, XL, YL, XR, YR, Y)
C
C  ONE-LINER:  Inverse parabolic interpolation given two points
C
C  PURPOSE:
C
C        INVERSE_PARABOLA performs inverse parabolic interpolation of
C     y vs. x given two points on the upper arm of a parabola of the form
C
C                          y ** 2  =  a x  +  b
C
C     for one or more points.  The target points are intended to be between
C     the two given points, although X > XR is OK, and V <= X < XL is also
C     OK, where V is at the vertex.  V is expressible as a function of the
C     control points, but no test is made for avoiding the negative square
C     root that would be implied if X < V.
C
C         This utility was prompted by an aerodynamic boundary layer
C     application, and adapted from the earlier QINTRP.
C
C  ARGUMENTS:
C     ARG   DIM   TYPE I/O/S DESCRIPTION
C     NX     -      I    I   No. of interpolated elements required (>=1)
C     X      NX     R    I   Abscissa(s) at which to interpolate; X >= vertex
C     XL     -      R    I   Abscissa of "left" control point
C     YL     -      R    I   Ordinate  "    "    "    "   "
C     XR     -      R    I   Abscissa of "right" control point
C     YR     -      R    I   Ordinate  "    "    "    "   "
C     Y      NX     R    O   Interpolated value(s)
C    
C  HISTORY:  10/31/02  DAS  Initial adaptation of QINTRP.
C
C  AUTHOR:   David Saunders, ELORET/NASA Ames Research Center, Moffett Field, CA
C
C-------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN)  ::  NX
      REAL,    INTENT (IN)  ::  X (NX), XL, YL, XR, YR
      REAL,    INTENT (OUT) ::  Y (NX)

C     Local constant:

      REAL, PARAMETER :: ZERO = 0.

C     Local variables:

      INTEGER :: I
      REAL    :: A, B, DX

C     Execution:

C     The following expressions are obtained by solving a 2 x 2 system:

      DX =  XR - XL
      B  =  YL ** 2            ! This is correct if XL = 0.
      A  = (YR ** 2 - B) / DX
      IF (XL /= ZERO) B  = (B * XR - YR ** 2 * XL) / DX

      DO I = 1, NX
         Y (I) = SQRT (A * X (I) + B)
      END DO

      END SUBROUTINE INVERSE_PARABOLA
