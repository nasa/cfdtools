C+------------------------------------------------------------------------------
C
      SUBROUTINE WARP3D (IMIN, IMAX, JMIN, JMAX, KMIN, KMAX,
     >                   I1, I2, J1, J2, K1, K2, X0, Y0, Z0, S0,
     >                   DFACEI, DFACEJ, DFACEK, X, Y, Z)
C
C  ONE-LINER: Perturb interior of a 3-space grid given new faces
C
C  DESCRIPTION:
C
C        WARP3D perturbs the interior points of a 3-space grid block
C     given the original grid and perturbed faces. All six faces are
C     assumed to be perturbed (though some may of course be fixed).
C     If all edges are unperturbed, considerably less work is needed,
C     but here some corners and/or edges are assumed to have been
C     perturbed, requiring a three-stage algorithm. Stage 1 accounts
C     for corner motion, stage 2 for edge motion, and stage 3 matches
C     the desired faces.
C
C        The perturbed faces should be input as faces of the desired
C     output grid.  The original relative arc-length increments in
C     each index direction should also be input.  See PARAMXYZ for
C     setting them up in preparation for multiple perturbations.
C
C  ENVIRONMENT:  FORTRAN 77 with minor extensions
C
C  HISTORY:
C
C     03/04/94  DAS/JJR  WARP3D written to perturb a single 3-space line,
C                        for James's wing/body design code.
C     12/07/94   "   "   Rewritten to apply the idea at a higher level,
C                        for efficient perturbation of a whole block.
C     02/01/96    DAS    Use on a subgrid with full-grid S0 must transform S0.
C     02/10/96     "     It's now a 3-stage algorithm (corners, edges, faces)
C                        where originally the first 2 stages were blurred.
C     02/12/96     "     Storing all the intermediate face perturbations allows
C                        a single pass through the volume points.
C     03/29/96  JJR/DAS  Stage 2 is done more consistently now by distinguishing
C                        the two types of perturbation affecting each face.
C                        This required adding DFACE*(*,*,*,*,4).
C
C  AUTHOR:  David Saunders/James Reuther, NASA Ames, Mt. View, CA.
C
C ------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments.

      INTEGER IMIN, IMAX, JMIN, JMAX, KMIN, KMAX     ! I Grid array dimensions.

      INTEGER I1, I2, J1, J2, K1, K2                 ! I Define active volume.

      REAL    X0 (IMIN:IMAX, JMIN:JMAX, KMIN:KMAX),  ! I Original grid coords.
     >        Y0 (IMIN:IMAX, JMIN:JMAX, KMIN:KMAX),
     >        Z0 (IMIN:IMAX, JMIN:JMAX, KMIN:KMAX)

      REAL    S0 (IMIN:IMAX, JMIN:JMAX, KMIN:KMAX, 3)! I Relative arc-lengths
                                                     !   for the I,J,K lines
                                                     !   see PARAMXYZ.  If from
                                                     !   a grid larger than the
                                                     !   active subgrid, S0 is
                                                     !   transformed here.

      REAL    DFACEI (3, JMIN:JMAX, KMIN:KMAX, 2, 4),! S For face perturbations:
     >        DFACEJ (3, IMIN:IMAX, KMIN:KMAX, 2, 4),!   DFACEI(1:3,*,*,1,M) =
     >        DFACEK (3, IMIN:IMAX, JMIN:JMAX, 2, 4) !   dX,dY,dZ on the I = I1
                                                     !   face for stage M, etc.;
                                                     !   M = 1, 2a, 2b, 3.

      REAL    X (IMIN:IMAX, JMIN:JMAX, KMIN:KMAX),   !I/O Grid coords.: input
     >        Y (IMIN:IMAX, JMIN:JMAX, KMIN:KMAX),   !   with faces perturbed;
     >        Z (IMIN:IMAX, JMIN:JMAX, KMIN:KMAX)    !   output fully perturbed.

C-------------------------------------------------------------------------------

C     Local constants.

      REAL    EPS, ONE

      PARAMETER (EPS = 1.E-8, ONE = 1.E+0) ! EPS safeguards a divide by zero -
                                           ! presumably only if result is zero.
C     Local variables.

      INTEGER I, J, K, L
      REAL    DEL1, DEL2, DELI, DELJ, DELK,
     >        DELIJ, DELIK, DELJI, DELJK, DELKI, DELKJ,
     >        WTI1, WTI2, WTJ1, WTJ2, WTK1, WTK2


C     Execution.
C     ----------

C     Set up the face perturbations for all three stages:

C     Stage 1:  Corner motion (only). This stage is reusable by WARPQ3D.
C     --------

C     I = I1 and I2 intermediate faces.

      CALL DELQ3D (IMIN, IMAX, JMIN, JMAX, KMIN, KMAX,
     >             I1, I1, J1, J2, K1, K2, X0, Y0, Z0, S0,
     >             DFACEI (1,JMIN,KMIN,1,1), DFACEJ, DFACEK, X, Y, Z)

      CALL DELQ3D (IMIN, IMAX, JMIN, JMAX, KMIN, KMAX,
     >             I2, I2, J1, J2, K1, K2, X0, Y0, Z0, S0,
     >             DFACEI (1,JMIN,KMIN,2,1), DFACEJ, DFACEK, X, Y, Z)

C     J = J1 and J2 intermediate faces.

      CALL DELQ3D (IMIN, IMAX, JMIN, JMAX, KMIN, KMAX,
     >             I1, I2, J1, J1, K1, K2, X0, Y0, Z0, S0,
     >             DFACEI, DFACEJ (1,IMIN,KMIN,1,1), DFACEK, X, Y, Z)

      CALL DELQ3D (IMIN, IMAX, JMIN, JMAX, KMIN, KMAX,
     >             I1, I2, J2, J2, K1, K2, X0, Y0, Z0, S0,
     >             DFACEI, DFACEJ (1,IMIN,KMIN,2,1), DFACEK, X, Y, Z)

C     K = K1 and K2 intermediate faces.

      CALL DELQ3D (IMIN, IMAX, JMIN, JMAX, KMIN, KMAX,
     >             I1, I2, J1, J2, K1, K1, X0, Y0, Z0, S0,
     >             DFACEI, DFACEJ, DFACEK (1,IMIN,JMIN,1,1), X, Y, Z)

      CALL DELQ3D (IMIN, IMAX, JMIN, JMAX, KMIN, KMAX,
     >             I1, I2, J1, J2, K2, K2, X0, Y0, Z0, S0,
     >             DFACEI, DFACEJ, DFACEK (1,IMIN,JMIN,2,1), X, Y, Z)


C     Stage 2:  Handle edge motion from above interim edges to final edges.
C     --------

C     James's insight here: consider the 4 faces affected by the motion of
C     4 edges in a given (index) direction; furthermore, keep the two
C     directions of the resulting face perturbations separated for proper
C     weighted combination during the final interpolation into the interior.
C     The 3 index directions are then independent of each other (added).

C     I1 and I2 faces:

      L = 1
      DO I = I1, I2, I2 - I1

C        K = K1 and K2 edge perturbations (J direction edges):

         DO J = J1 + 1, J2 - 1
            DFACEI(1,J,K1,L,2) = X(I,J,K1)-X0(I,J,K1)-DFACEI(1,J,K1,L,1)
            DFACEI(2,J,K1,L,2) = Y(I,J,K1)-Y0(I,J,K1)-DFACEI(2,J,K1,L,1)
            DFACEI(3,J,K1,L,2) = Z(I,J,K1)-Z0(I,J,K1)-DFACEI(3,J,K1,L,1)
            DFACEI(1,J,K2,L,2) = X(I,J,K2)-X0(I,J,K2)-DFACEI(1,J,K2,L,1)
            DFACEI(2,J,K2,L,2) = Y(I,J,K2)-Y0(I,J,K2)-DFACEI(2,J,K2,L,1)
            DFACEI(3,J,K2,L,2) = Z(I,J,K2)-Z0(I,J,K2)-DFACEI(3,J,K2,L,1)
         END DO

C        J = J1 and J2 edge perturbations (K direction edges):

         DO K = K1 + 1, K2 - 1
            DFACEI(1,J1,K,L,3) = X(I,J1,K)-X0(I,J1,K)-DFACEI(1,J1,K,L,1)
            DFACEI(2,J1,K,L,3) = Y(I,J1,K)-Y0(I,J1,K)-DFACEI(2,J1,K,L,1)
            DFACEI(3,J1,K,L,3) = Z(I,J1,K)-Z0(I,J1,K)-DFACEI(3,J1,K,L,1)
            DFACEI(1,J2,K,L,3) = X(I,J2,K)-X0(I,J2,K)-DFACEI(1,J2,K,L,1)
            DFACEI(2,J2,K,L,3) = Y(I,J2,K)-Y0(I,J2,K)-DFACEI(2,J2,K,L,1)
            DFACEI(3,J2,K,L,3) = Z(I,J2,K)-Z0(I,J2,K)-DFACEI(3,J2,K,L,1)
         END DO

C        Interpolate stage 2 interior points for this I face, keeping the
C        two index directions separated.

         DO K = K1 + 1, K2 - 1
            DO J = J1 + 1, J2 - 1
               WTJ2 = (S0 (I, J,  K, 2) - S0 (I, J1, K, 2)) /
     >                (S0 (I, J2, K, 2) - S0 (I, J1, K, 2))
               WTJ1 = ONE - WTJ2
               WTK2 = (S0 (I, J, K,  3) - S0 (I, J, K1, 3)) /
     >                (S0 (I, J, K2, 3) - S0 (I, J, K1, 3))
               WTK1 = ONE - WTK2

               DFACEI(1,J,K,L,3) = WTJ1*DFACEI(1,J1,K,L,3) +
     >                             WTJ2*DFACEI(1,J2,K,L,3)
               DFACEI(2,J,K,L,3) = WTJ1*DFACEI(2,J1,K,L,3) +
     >                             WTJ2*DFACEI(2,J2,K,L,3)
               DFACEI(3,J,K,L,3) = WTJ1*DFACEI(3,J1,K,L,3) +
     >                             WTJ2*DFACEI(3,J2,K,L,3)

               DFACEI(1,J,K,L,2) = WTK1*DFACEI(1,J,K1,L,2) +
     >                             WTK2*DFACEI(1,J,K2,L,2)
               DFACEI(2,J,K,L,2) = WTK1*DFACEI(2,J,K1,L,2) +
     >                             WTK2*DFACEI(2,J,K2,L,2)
               DFACEI(3,J,K,L,2) = WTK1*DFACEI(3,J,K1,L,2) +
     >                             WTK2*DFACEI(3,J,K2,L,2)
            END DO
         END DO
         L = 2
      END DO

C     J1 and J2 faces, stage 2:

      L = 1
      DO J = J1, J2, J2 - J1

C        K = K1 and K2 edge perturbations (I direction):

         DO I = I1 + 1, I2-1
            DFACEJ(1,I,K1,L,2) = X(I,J,K1)-X0(I,J,K1)-DFACEJ(1,I,K1,L,1)
            DFACEJ(2,I,K1,L,2) = Y(I,J,K1)-Y0(I,J,K1)-DFACEJ(2,I,K1,L,1)
            DFACEJ(3,I,K1,L,2) = Z(I,J,K1)-Z0(I,J,K1)-DFACEJ(3,I,K1,L,1)
            DFACEJ(1,I,K2,L,2) = X(I,J,K2)-X0(I,J,K2)-DFACEJ(1,I,K2,L,1)
            DFACEJ(2,I,K2,L,2) = Y(I,J,K2)-Y0(I,J,K2)-DFACEJ(2,I,K2,L,1)
            DFACEJ(3,I,K2,L,2) = Z(I,J,K2)-Z0(I,J,K2)-DFACEJ(3,I,K2,L,1)
         END DO

C        I = I1 and I2 edge perturbations (K direction):

         DO K = K1 + 1, K2 - 1
            DFACEJ(1,I1,K,L,3) = X(I1,J,K)-X0(I1,J,K)-DFACEJ(1,I1,K,L,1)
            DFACEJ(2,I1,K,L,3) = Y(I1,J,K)-Y0(I1,J,K)-DFACEJ(2,I1,K,L,1)
            DFACEJ(3,I1,K,L,3) = Z(I1,J,K)-Z0(I1,J,K)-DFACEJ(3,I1,K,L,1)
            DFACEJ(1,I2,K,L,3) = X(I2,J,K)-X0(I2,J,K)-DFACEJ(1,I2,K,L,1)
            DFACEJ(2,I2,K,L,3) = Y(I2,J,K)-Y0(I2,J,K)-DFACEJ(2,I2,K,L,1)
            DFACEJ(3,I2,K,L,3) = Z(I2,J,K)-Z0(I2,J,K)-DFACEJ(3,I2,K,L,1)
         END DO

C        Interpolate stage 2 interior points for this J face.

         DO K = K1 + 1, K2 - 1
            DO I = I1 + 1, I2 - 1
               WTI2 = (S0 (I,  J, K, 1) - S0 (I1, J, K, 1)) /
     >                (S0 (I2, J, K, 1) - S0 (I1, J, K, 1))
               WTI1 = ONE - WTI2
               WTK2 = (S0 (I, J, K,  3) - S0 (I, J, K1, 3)) /
     >                (S0 (I, J, K2, 3) - S0 (I, J, K1, 3))
               WTK1 = ONE - WTK2

               DFACEJ(1,I,K,L,3) = WTI1*DFACEJ(1,I1,K,L,3) +
     >                             WTI2*DFACEJ(1,I2,K,L,3)
               DFACEJ(2,I,K,L,3) = WTI1*DFACEJ(2,I1,K,L,3) +
     >                             WTI2*DFACEJ(2,I2,K,L,3)
               DFACEJ(3,I,K,L,3) = WTI1*DFACEJ(3,I1,K,L,3) +
     >                             WTI2*DFACEJ(3,I2,K,L,3)

               DFACEJ(1,I,K,L,2) = WTK1*DFACEJ(1,I,K1,L,2) +
     >                             WTK2*DFACEJ(1,I,K2,L,2)
               DFACEJ(2,I,K,L,2) = WTK1*DFACEJ(2,I,K1,L,2) +
     >                             WTK2*DFACEJ(2,I,K2,L,2)
               DFACEJ(3,I,K,L,2) = WTK1*DFACEJ(3,I,K1,L,2) +
     >                             WTK2*DFACEJ(3,I,K2,L,2)
            END DO
         END DO
         L = 2
      END DO

C     K1 and K2 faces, stage 2:

      L = 1
      DO K = K1, K2, K2 - K1

C        J = J1 and J2 edge perturbations (I direction):

         DO I = I1 + 1, I2 - 1
            DFACEK(1,I,J1,L,2) = X(I,J1,K)-X0(I,J1,K)-DFACEK(1,I,J1,L,1)
            DFACEK(2,I,J1,L,2) = Y(I,J1,K)-Y0(I,J1,K)-DFACEK(2,I,J1,L,1)
            DFACEK(3,I,J1,L,2) = Z(I,J1,K)-Z0(I,J1,K)-DFACEK(3,I,J1,L,1)
            DFACEK(1,I,J2,L,2) = X(I,J2,K)-X0(I,J2,K)-DFACEK(1,I,J2,L,1)
            DFACEK(2,I,J2,L,2) = Y(I,J2,K)-Y0(I,J2,K)-DFACEK(2,I,J2,L,1)
            DFACEK(3,I,J2,L,2) = Z(I,J2,K)-Z0(I,J2,K)-DFACEK(3,I,J2,L,1)
         END DO

C        I = I1 and I2 edge perturbations (J direction):

         DO J = J1 + 1, J2 - 1
            DFACEK(1,I1,J,L,3) = X(I1,J,K)-X0(I1,J,K)-DFACEK(1,I1,J,L,1)
            DFACEK(2,I1,J,L,3) = Y(I1,J,K)-Y0(I1,J,K)-DFACEK(2,I1,J,L,1)
            DFACEK(3,I1,J,L,3) = Z(I1,J,K)-Z0(I1,J,K)-DFACEK(3,I1,J,L,1)
            DFACEK(1,I2,J,L,3) = X(I2,J,K)-X0(I2,J,K)-DFACEK(1,I2,J,L,1)
            DFACEK(2,I2,J,L,3) = Y(I2,J,K)-Y0(I2,J,K)-DFACEK(2,I2,J,L,1)
            DFACEK(3,I2,J,L,3) = Z(I2,J,K)-Z0(I2,J,K)-DFACEK(3,I2,J,L,1)
         END DO

C        Interpolate stage 2 interior points for this K face.

         DO J = J1 + 1, J2 - 1
            DO I = I1 + 1, I2 - 1
               WTI2 = (S0 (I,  J, K, 1) - S0 (I1, J, K, 1)) /
     >                (S0 (I2, J, K, 1) - S0 (I1, J, K, 1))
               WTI1 = ONE - WTI2
               WTJ2 = (S0 (I, J,  K, 2) - S0 (I, J1, K, 2)) /
     >                (S0 (I, J2, K, 2) - S0 (I, J1, K, 2))
               WTJ1 = ONE - WTJ2

               DFACEK(1,I,J,L,3) = WTI1*DFACEK(1,I1,J,L,3) +
     >                             WTI2*DFACEK(1,I2,J,L,3)
               DFACEK(2,I,J,L,3) = WTI1*DFACEK(2,I1,J,L,3) +
     >                             WTI2*DFACEK(2,I2,J,L,3)
               DFACEK(3,I,J,L,3) = WTI1*DFACEK(3,I1,J,L,3) +
     >                             WTI2*DFACEK(3,I2,J,L,3)

               DFACEK(1,I,J,L,2) = WTJ1*DFACEK(1,I,J1,L,2) +
     >                             WTJ2*DFACEK(1,I,J2,L,2)
               DFACEK(2,I,J,L,2) = WTJ1*DFACEK(2,I,J1,L,2) +
     >                             WTJ2*DFACEK(2,I,J2,L,2)
               DFACEK(3,I,J,L,2) = WTJ1*DFACEK(3,I,J1,L,2) +
     >                             WTJ2*DFACEK(3,I,J2,L,2)
            END DO
         END DO
         L = 2
      END DO


C     Stage 3:  Handle face motion from above interim faces to final faces.
C     --------

      DO K = K1 + 1, K2 - 1
         DO J = J1 + 1, J2 - 1
            DFACEI(1,J,K,1,4) = X(I1,J,K) - X0(I1,J,K) -
     >         DFACEI(1,J,K,1,1) - DFACEI(1,J,K,1,2) - DFACEI(1,J,K,1,3)
            DFACEI(2,J,K,1,4) = Y(I1,J,K) - Y0(I1,J,K) -
     >         DFACEI(2,J,K,1,1) - DFACEI(2,J,K,1,2) - DFACEI(2,J,K,1,3)
            DFACEI(3,J,K,1,4) = Z(I1,J,K) - Z0(I1,J,K) -
     >         DFACEI(3,J,K,1,1) - DFACEI(3,J,K,1,2) - DFACEI(3,J,K,1,3)
            DFACEI(1,J,K,2,4) = X(I2,J,K) - X0(I2,J,K) -
     >         DFACEI(1,J,K,2,1) - DFACEI(1,J,K,2,2) - DFACEI(1,J,K,2,3)
            DFACEI(2,J,K,2,4) = Y(I2,J,K) - Y0(I2,J,K) -
     >         DFACEI(2,J,K,2,1) - DFACEI(2,J,K,2,2) - DFACEI(2,J,K,2,3)
            DFACEI(3,J,K,2,4) = Z(I2,J,K) - Z0(I2,J,K) -
     >         DFACEI(3,J,K,2,1) - DFACEI(3,J,K,2,2) - DFACEI(3,J,K,2,3)
         END DO

         DO I = I1 + 1, I2 - 1
            DFACEJ(1,I,K,1,4) = X(I,J1,K) - X0(I,J1,K) -
     >         DFACEJ(1,I,K,1,1) - DFACEJ(1,I,K,1,2) - DFACEJ(1,I,K,1,3)
            DFACEJ(2,I,K,1,4) = Y(I,J1,K) - Y0(I,J1,K) -
     >         DFACEJ(2,I,K,1,1) - DFACEJ(2,I,K,1,2) - DFACEJ(2,I,K,1,3)
            DFACEJ(3,I,K,1,4) = Z(I,J1,K) - Z0(I,J1,K) -
     >         DFACEJ(3,I,K,1,1) - DFACEJ(3,I,K,1,2) - DFACEJ(3,I,K,1,3)
            DFACEJ(1,I,K,2,4) = X(I,J2,K) - X0(I,J2,K) -
     >         DFACEJ(1,I,K,2,1) - DFACEJ(1,I,K,2,2) - DFACEJ(1,I,K,2,3)
            DFACEJ(2,I,K,2,4) = Y(I,J2,K) - Y0(I,J2,K) -
     >         DFACEJ(2,I,K,2,1) - DFACEJ(2,I,K,2,2) - DFACEJ(2,I,K,2,3)
            DFACEJ(3,I,K,2,4) = Z(I,J2,K) - Z0(I,J2,K) -
     >         DFACEJ(3,I,K,2,1) - DFACEJ(3,I,K,2,2) - DFACEJ(3,I,K,2,3)
         END DO
      END DO

      DO J = J1 + 1, J2 - 1
         DO I = I1 + 1, I2 - 1
            DFACEK(1,I,J,1,4) = X(I,J,K1) - X0(I,J,K1) -
     >         DFACEK(1,I,J,1,1) - DFACEK(1,I,J,1,2) - DFACEK(1,I,J,1,3)
            DFACEK(2,I,J,1,4) = Y(I,J,K1) - Y0(I,J,K1) -
     >         DFACEK(2,I,J,1,1) - DFACEK(2,I,J,1,2) - DFACEK(2,I,J,1,3)
            DFACEK(3,I,J,1,4) = Z(I,J,K1) - Z0(I,J,K1) -
     >         DFACEK(3,I,J,1,1) - DFACEK(3,I,J,1,2) - DFACEK(3,I,J,1,3)
            DFACEK(1,I,J,2,4) = X(I,J,K2) - X0(I,J,K2) -
     >         DFACEK(1,I,J,2,1) - DFACEK(1,I,J,2,2) - DFACEK(1,I,J,2,3)
            DFACEK(2,I,J,2,4) = Y(I,J,K2) - Y0(I,J,K2) -
     >         DFACEK(2,I,J,2,1) - DFACEK(2,I,J,2,2) - DFACEK(2,I,J,2,3)
            DFACEK(3,I,J,2,4) = Z(I,J,K2) - Z0(I,J,K2) -
     >         DFACEK(3,I,J,2,1) - DFACEK(3,I,J,2,2) - DFACEK(3,I,J,2,3)
         END DO
      END DO


C     Perturb the interior volume points.
C     All stages are performed at once via interpolation from the face
C     perturbations stored for each stage.  Note that the three stages
C     accumulate the contributions from the three subscript directions
C     with varying degrees of independence.

      DO K = K1 + 1, K2 - 1
         DO J = J1 + 1, J2 - 1
            DO I = I1 + 1, I2 - 1
               WTI2 = (S0 (I,  J, K, 1) - S0 (I1, J, K, 1)) /
     >                (S0 (I2, J, K, 1) - S0 (I1, J, K, 1))
               WTI1 = ONE - WTI2
               WTJ2 = (S0 (I, J,  K, 2) - S0 (I, J1, K, 2)) /
     >                (S0 (I, J2, K, 2) - S0 (I, J1, K, 2))
               WTJ1 = ONE - WTJ2
               WTK2 = (S0 (I, J, K,  3) - S0 (I, J, K1, 3)) /
     >                (S0 (I, J, K2, 3) - S0 (I, J, K1, 3))
               WTK1 = ONE - WTK2

C              Stage 1:

               DELI = WTI1*DFACEI (1,J,K,1,1) + WTI2*DFACEI (1,J,K,2,1)
               DELJ = WTJ1*DFACEJ (1,I,K,1,1) + WTJ2*DFACEJ (1,I,K,2,1)
               DELK = WTK1*DFACEK (1,I,J,1,1) + WTK2*DFACEK (1,I,J,2,1)

               DEL1 = (ABS(DELI)*DELI + ABS(DELJ)*DELJ + ABS(DELK)*DELK)
     >                / MAX (ABS (DELI) + ABS (DELJ) + ABS (DELK), EPS)

C              Stage 2:

               DELIJ= WTJ1*DFACEJ (1,I,K,1,2) + WTJ2*DFACEJ (1,I,K,2,2)
               DELIK= WTK1*DFACEK (1,I,J,1,2) + WTK2*DFACEK (1,I,J,2,2)
               DELI = (ABS (DELIJ) * DELIJ + ABS (DELIK) * DELIK) /
     >                 MAX (ABS (DELIJ) + ABS (DELIK), EPS)

               DELJI= WTI1*DFACEI (1,J,K,1,2) + WTI2*DFACEI (1,J,K,2,2)
               DELJK= WTK1*DFACEK (1,I,J,1,3) + WTK2*DFACEK (1,I,J,2,3)
               DELJ = (ABS (DELJI) * DELJI + ABS (DELJK) * DELJK) /
     >                 MAX (ABS (DELJI) + ABS (DELJK), EPS)

               DELKI= WTI1*DFACEI (1,J,K,1,3) + WTI2*DFACEI (1,J,K,2,3)
               DELKJ= WTJ1*DFACEJ (1,I,K,1,3) + WTJ2*DFACEJ (1,I,K,2,3)
               DELK = (ABS (DELKI) * DELKI + ABS (DELKJ) * DELKJ) /
     >                 MAX (ABS (DELKI) + ABS (DELKJ), EPS)

               DEL2 = DELI + DELJ + DELK

C              Stage 3:

               DELI = WTI1*DFACEI (1,J,K,1,4) + WTI2*DFACEI (1,J,K,2,4)
               DELJ = WTJ1*DFACEJ (1,I,K,1,4) + WTJ2*DFACEJ (1,I,K,2,4)
               DELK = WTK1*DFACEK (1,I,J,1,4) + WTK2*DFACEK (1,I,J,2,4)

               X (I, J, K) = X0 (I, J, K) + DEL1 + DEL2 + DELI+DELJ+DELK

C              Repeat all stages for Y perturbations:

               DELI = WTI1*DFACEI (2,J,K,1,1) + WTI2*DFACEI (2,J,K,2,1)
               DELJ = WTJ1*DFACEJ (2,I,K,1,1) + WTJ2*DFACEJ (2,I,K,2,1)
               DELK = WTK1*DFACEK (2,I,J,1,1) + WTK2*DFACEK (2,I,J,2,1)

               DEL1 = (ABS(DELI)*DELI + ABS(DELJ)*DELJ + ABS(DELK)*DELK)
     >                / MAX (ABS (DELI) + ABS (DELJ) + ABS (DELK), EPS)

               DELIJ= WTJ1*DFACEJ (2,I,K,1,2) + WTJ2*DFACEJ (2,I,K,2,2)
               DELIK= WTK1*DFACEK (2,I,J,1,2) + WTK2*DFACEK (2,I,J,2,2)
               DELI = (ABS (DELIJ) * DELIJ + ABS (DELIK) * DELIK) /
     >                 MAX (ABS (DELIJ) + ABS (DELIK), EPS)

               DELJI= WTI1*DFACEI (2,J,K,1,2) + WTI2*DFACEI (2,J,K,2,2)
               DELJK= WTK1*DFACEK (2,I,J,1,3) + WTK2*DFACEK (2,I,J,2,3)
               DELJ = (ABS (DELJI) * DELJI + ABS (DELJK) * DELJK) /
     >                 MAX (ABS (DELJI) + ABS (DELJK), EPS)

               DELKI= WTI1*DFACEI (2,J,K,1,3) + WTI2*DFACEI (2,J,K,2,3)
               DELKJ= WTJ1*DFACEJ (2,I,K,1,3) + WTJ2*DFACEJ (2,I,K,2,3)
               DELK = (ABS (DELKI) * DELKI + ABS (DELKJ) * DELKJ) /
     >                 MAX (ABS (DELKI) + ABS (DELKJ), EPS)

               DEL2 = DELI + DELJ + DELK

               DELI = WTI1*DFACEI (2,J,K,1,4) + WTI2*DFACEI (2,J,K,2,4)
               DELJ = WTJ1*DFACEJ (2,I,K,1,4) + WTJ2*DFACEJ (2,I,K,2,4)
               DELK = WTK1*DFACEK (2,I,J,1,4) + WTK2*DFACEK (2,I,J,2,4)

               Y (I, J, K) = Y0 (I, J, K) + DEL1 + DEL2 + DELI+DELJ+DELK

C              ... and for Z perturbations:

               DELI = WTI1*DFACEI (3,J,K,1,1) + WTI2*DFACEI (3,J,K,2,1)
               DELJ = WTJ1*DFACEJ (3,I,K,1,1) + WTJ2*DFACEJ (3,I,K,2,1)
               DELK = WTK1*DFACEK (3,I,J,1,1) + WTK2*DFACEK (3,I,J,2,1)

               DEL1 = (ABS(DELI)*DELI + ABS(DELJ)*DELJ + ABS(DELK)*DELK)
     >                / MAX (ABS (DELI) + ABS (DELJ) + ABS (DELK), EPS)

               DELIJ= WTJ1*DFACEJ (3,I,K,1,2) + WTJ2*DFACEJ (3,I,K,2,2)
               DELIK= WTK1*DFACEK (3,I,J,1,2) + WTK2*DFACEK (3,I,J,2,2)
               DELI = (ABS (DELIJ) * DELIJ + ABS (DELIK) * DELIK) /
     >                 MAX (ABS (DELIJ) + ABS (DELIK), EPS)

               DELJI= WTI1*DFACEI (3,J,K,1,2) + WTI2*DFACEI (3,J,K,2,2)
               DELJK= WTK1*DFACEK (3,I,J,1,3) + WTK2*DFACEK (3,I,J,2,3)
               DELJ = (ABS (DELJI) * DELJI + ABS (DELJK) * DELJK) /
     >                 MAX (ABS (DELJI) + ABS (DELJK), EPS)

               DELKI= WTI1*DFACEI (3,J,K,1,3) + WTI2*DFACEI (3,J,K,2,3)
               DELKJ= WTJ1*DFACEJ (3,I,K,1,3) + WTJ2*DFACEJ (3,I,K,2,3)
               DELK = (ABS (DELKI) * DELKI + ABS (DELKJ) * DELKJ) /
     >                 MAX (ABS (DELKI) + ABS (DELKJ), EPS)

               DEL2 = DELI + DELJ + DELK

               DELI = WTI1*DFACEI (3,J,K,1,4) + WTI2*DFACEI (3,J,K,2,4)
               DELJ = WTJ1*DFACEJ (3,I,K,1,4) + WTJ2*DFACEJ (3,I,K,2,4)
               DELK = WTK1*DFACEK (3,I,J,1,4) + WTK2*DFACEK (3,I,J,2,4)

               Z (I, J, K) = Z0 (I, J, K) + DEL1 + DEL2 + DELI+DELJ+DELK
            END DO
         END DO
      END DO

      RETURN
      END
