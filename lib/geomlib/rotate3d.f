C+----------------------------------------------------------------------
C
      SUBROUTINE ROTATE3D (N, X, Y, Z, ANGLE, PX, PY, PZ, QX, QY, QZ)
C
C     ROTATE3D rotates the given point(s) (X, Y, Z) about the straight
C     line joining points P (Px, Py, Pz) and Q (Qx, Qy, Qz) by the
C     indicated angle (degrees, right-hand rule).  Input coordinates
C     are overwritten.
C
C     HISTORY: 06/14/00  DAS  Adaptation of ROTATE2D and ROTATE program.
C
C     AUTHOR: David Saunders, ELORET/NASA Ames, Moffett Field, CA.
C
C-----------------------------------------------------------------------

      USE TRIGD

C     Arguments:

      INTEGER, INTENT (IN) ::
     >   N                    ! Number of points to rotate
      REAL, INTENT (INOUT) ::
     >   X (N), Y (N), Z (N)  ! Point coordinates, overwritten
      REAL, INTENT (IN) ::
     >   ANGLE,               ! Angle in degrees; RH rule (thumb P -> Q)
     >   PX, PY, PZ,          ! Points defining the rotation axis
     >   QX, QY, QZ

C     Local constants:

      REAL, PARAMETER ::
     >   ONE = 1.0, ZERO = 0.0

C     Local variables:

      INTEGER
     >   I, J, K
      REAL
     >   CTH, STH, VTH, D, VEC (3), VECROT (3), ROT (3, 3)

C     Execution:

      CTH = COSD (ANGLE)
      STH = SIND (ANGLE)
      VTH = ONE - CTH

C     Derive a unit vector translated to the origin from the vector PQ:

      VEC (1) = QX - PX
      VEC (2) = QY - PY
      VEC (3) = QZ - PZ
      D = SQRT (VEC (1) ** 2 + VEC (2) ** 2 + VEC (3) ** 2)
      VEC = VEC / D

C     Rotation matrix ("Introduction to Robotics" by John Craig, 1986):

      ROT (1, 1) = VEC (1) * VEC (1) * VTH + CTH
      ROT (1, 2) = VEC (1) * VEC (2) * VTH - VEC (3) * STH
      ROT (1, 3) = VEC (1) * VEC (3) * VTH + VEC (2) * STH

      ROT (2, 1) = VEC (1) * VEC (2) * VTH + VEC (3) * STH
      ROT (2, 2) = VEC (2) * VEC (2) * VTH + CTH
      ROT (2, 3) = VEC (3) * VEC (2) * VTH - VEC (1) * STH

      ROT (3, 1) = VEC (1) * VEC (3) * VTH - VEC (2) * STH
      ROT (3, 2) = VEC (2) * VEC (3) * VTH + VEC (1) * STH
      ROT (3, 3) = VEC (3) * VEC (3) * VTH + CTH

C     Translate each point by P's coordinates, rotate about the unit axis
C     at the origin, then translate back:

      DO K = 1, N

         VEC (1) = X (K) - PX
         VEC (2) = Y (K) - PY
         VEC (3) = Z (K) - PZ

         DO I = 1, 3
            VECROT (I) = ZERO
            DO J = 1, 3
               VECROT (I) = VECROT (I) + ROT (I, J) * VEC (J)
            END DO
         END DO

         X (K) = VECROT (1) + PX
         Y (K) = VECROT (2) + PY
         Z (K) = VECROT (3) + PZ

      END DO

      END SUBROUTINE ROTATE3D
