C+------------------------------------------------------------------------------
C
      SUBROUTINE NULINE3D (I1, I2, X0, Y0, Z0, X, Y, Z)
C
C  ONE-LINER: Adjust interior points of a 3-space line given new end points
C
C  DESCRIPTION:
C
C        NULINE3D perturbs the interior points of a 3-space curve given
C     perturbed end points.  More precisely, any pair of points I1, I2
C     may be the controlling points; each point between these two is moved
C     according to an arc-length-based combination of the distances between
C     the original and new points corresponding to I1 and I2.  The new
C     coordinates of these two points should be input in the desired output
C     curve.  In-place adjustment is NOT permitted.
C
C  ENVIRONMENT:  Fortran 90
C
C  HISTORY:
C
C     03/04/94  DAS  Implemented for James Reuther's grid perturbation scheme.
C     01/09/95   "   Renamed NULINE3D from WARP3D since WARP3D now applies
C                    to all points of a 3-space grid block; I1, I2, not 1, N.
C     03/16/98   "   Eliminated ARC argument - make it an automatic array.
C
C  AUTHOR:  David Saunders, Sterling Software/NASA Ames, Mt. View, CA.
C
C ------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER I1, I2                  !  I  Points given as perturbed
      REAL    X0 (*), Y0 (*), Z0 (*)  !  I  Original curve coordinates
      REAL    X (*), Y (*), Z (*)     ! I/O Desired curve, with points
                                      !     I1, I2 input

C-------------------------------------------------------------------------------

C     Local constants:

      REAL, PARAMETER :: ONE = 1., ZERO = 0.

C     Local variables:

      INTEGER I
      REAL    ARC(I1:I2), DX1, DX2, DY1, DY2, DZ1, DZ2, SINV, W1, W2

C     Execution:

      ARC (I1) = ZERO
      DO I = I1 + 1, I2
         ARC (I) = ARC (I - 1) + SQRT (
     >      (X0 (I) - X0 (I - 1)) ** 2 + (Y0 (I) - Y0 (I - 1)) ** 2 +
     >      (Z0 (I) - Z0 (I - 1)) ** 2)
      END DO

      DX1 = X (I1) - X0 (I1)
      DX2 = X (I2) - X0 (I2)
      DY1 = Y (I1) - Y0 (I1)
      DY2 = Y (I2) - Y0 (I2)
      DZ1 = Z (I1) - Z0 (I1)
      DZ2 = Z (I2) - Z0 (I2)
      SINV = ONE / ARC (I2)

      DO I = I1 + 1, I2 - 1
         W2 = ARC (I) * SINV
         W1 = ONE - W2
         X (I) = X0 (I) + W1 * DX1 + W2 * DX2
         Y (I) = Y0 (I) + W1 * DY1 + W2 * DY2
         Z (I) = Z0 (I) + W1 * DZ1 + W2 * DZ2
      END DO

      END SUBROUTINE NULINE3D
