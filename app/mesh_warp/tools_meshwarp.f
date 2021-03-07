C+-------------------------------------------------------------------------

      SUBROUTINE WARPBLK_GM (IFACEPTB, IEDGEPTB, IL, JL, KL, XYZ0, XYZ)

C     ******************************************************************
C     *   WARPBLK completes the perturbation of one block of a multi-  *
C     *   block grid structure given a base grid and all fixed and all *
C     *   explicitly-perturbed faces in-place for the desired block.   *
C     *                                                                *
C     *   Ancillary routines:   PARAM3DM, WARPQ3DM, DELQ3DM            *
C     *                                                                *
C     *   11/29/95  D.Saunders  Specialized adaptation of WARP3D and   *
C     *                         ancillary routines for FLO107-MB.      *
C     *   01/26/96  DAS/JJR     Edges affected implicitly by corner    *
C     *                         motion had been overlooked.            *
C     *   04/04/96     "        Algorithm is now 3-stage/1-pass, not   *
C     *                         2-stage/2-pass as it was originally.   *
C     *   05/24/96     "        Juan Alonso debugged it for us.        *
C     *   06/??/96    JJR       Implicit edge motion is specified at   *
C     *                         the higher level now via IEDGEPTB(*).  *
C     *   06/19/96    DAS       Handled degenerate edges in PARAM3DM.  *
C     *   10/99    M.Rimlinger  GMORPH/MESHWARP version does not need  *
C     *                         halos; allocate s0, dfacei,j,k instead *
C     *                         of passing them as arguments; killed   *
C     *                         ig, ig, and ncall.                     *
C     *                                                                *
C     *   David Saunders/James Reuther, NASA Ames Research Center, CA. *
C     ******************************************************************

      IMPLICIT NONE

C     Arguments:

      INTEGER   IFACEPTB(6)                 ! I  Controls for faces 1 - 6:
                                            !    0 means no perturbation;
                                            !    1 means implicit perturbation;
                                            !    2 means explicit perturbation
      INTEGER   IEDGEPTB(12)                ! I  Controls for edges 1 - 12:
                                            !    0 means no perturbation;
                                            !    1 means implicit perturbation;
                                            !    2 means explicit perturbation
      INTEGER   NCALL                       ! I  First call if NCALL = -3
      INTEGER   IG, IGO                     ! I  Design point i.d.
      INTEGER   IL, JL, KL                  ! I  Block dimensions
      REAL      XYZ0(3,IL,JL,KL)            ! I  Base block
      REAL      S0(3,IL,JL,KL)              !I/O Base normalized arc-lengths;
      REAL      DFACEI(3,JL,KL,2,4)         ! S  For face perturbations; e.g.,
      REAL      DFACEJ(3,IL,KL,2,4)         !    DFACEI(1:3,JL,KL,1,M) =
      REAL      DFACEK(3,IL,JL,2,4)         !    dX, dY, dZ on the I=1 face for
                                            !    stage M, etc.; M = 1,2a 2b,3.
                                            !    each can be (3,MDM,MDM,2,4) if
                                            !    MDM=MAX(IL,JL,KL) (all blocks)
      REAL      XYZ(3,IL,JL,KL)             !I/O Perturbed block

C     Local constants:

      REAL, PARAMETER ::
     >          EPS = 1.E-14, ONE = 1.E+0,  ! EPS safeguards a divide by zero -
     >          ZERO = 0.E+0                ! presumably only if result is zero.

C     Local variables:

      INTEGER   I, I1, I2, J, J1, J2, K, K1, K2, L, LC, LE, M
      REAL      DEL1, DEL2, DELI, DELJ, DELK,
     >          DELIJ, DELIK, DELJI, DELJK, DELKI, DELKJ,
     >          WTI1, WTI2, WTJ1, WTJ2, WTK1, WTK2
      REAL      DX(8), DY(8), DZ(8)         ! Corner perturbations

C     Execution:

C     Parameterize the block volume on the first call:

C     IF (NCALL <= -3 .AND. IG /= IGO)
C    >   CALL PARAM3DM (IL, JL, KL, XYZ0, S0)

      CALL PARAM3DM_GM (IL, JL, KL, XYZ0, S0)

C     Calculate corner point motion:

      LC = 1
      DO K = 1, KL, KL - 1
         DO J = 1, JL, JL - 1
            DO I = 1, IL, IL - 1
               DX(LC) = XYZ(1,I,J,K) - XYZ0(1,I,J,K)
               DY(LC) = XYZ(2,I,J,K) - XYZ0(2,I,J,K)
               DZ(LC) = XYZ(3,I,J,K) - XYZ0(3,I,J,K)
               LC = LC + 1
            END DO
         END DO
      END DO

C     Perturb implicit block edges.

C     Edges in the I direction:

      LE = 1
      LC = 1

      DO K = 1, KL, KL - 1
         DO J = 1, JL, JL - 1
            IF (IEDGEPTB(LE) == 1) THEN
               DO I = 2, IL - 1
                  WTI2 = S0(1,I,J,K)
                  WTI1 = ONE - WTI2
                  XYZ(1,I,J,K)=WTI1*DX(LC) +WTI2*DX(LC+1) +XYZ0(1,I,J,K)
                  XYZ(2,I,J,K)=WTI1*DY(LC) +WTI2*DY(LC+1) +XYZ0(2,I,J,K)
                  XYZ(3,I,J,K)=WTI1*DZ(LC) +WTI2*DZ(LC+1) +XYZ0(3,I,J,K)
               END DO
            END IF
            LE = LE + 1
            LC = LC + 2
         END DO
      END DO

C     Edges in the J direction:

      LE = 5
      LC = 1

      DO K = 1, KL, KL - 1
         DO I = 1, IL, IL - 1
            IF (IEDGEPTB(LE) == 1) THEN
               DO J = 2, JL - 1
                  WTJ2 = S0(2,I,J,K)
                  WTJ1 = ONE - WTJ2
                  XYZ(1,I,J,K)=WTJ1*DX(LC) +WTJ2*DX(LC+2) +XYZ0(1,I,J,K)
                  XYZ(2,I,J,K)=WTJ1*DY(LC) +WTJ2*DY(LC+2) +XYZ0(2,I,J,K)
                  XYZ(3,I,J,K)=WTJ1*DZ(LC) +WTJ2*DZ(LC+2) +XYZ0(3,I,J,K)
               END DO
            END IF
            LE = LE + 1
            LC = LC + 1
         END DO
         LC = LC + 2
      END DO

C     Edges in the K direction:

      LE = 9
      LC = 1

      DO J = 1, JL, JL - 1
         DO I = 1, IL, IL - 1
            IF (IEDGEPTB(LE) == 1 ) THEN
               DO K = 2, KL - 1
                  WTK2 = S0(3,I,J,K)
                  WTK1 = ONE - WTK2
                  XYZ(1,I,J,K)=WTK1*DX(LC) +WTK2*DX(LC+4) +XYZ0(1,I,J,K)
                  XYZ(2,I,J,K)=WTK1*DY(LC) +WTK2*DY(LC+4) +XYZ0(2,I,J,K)
                  XYZ(3,I,J,K)=WTK1*DZ(LC) +WTK2*DZ(LC+4) +XYZ0(3,I,J,K)
               END DO
            END IF
            LE = LE + 1
            LC = LC + 1
         END DO
      END DO

C     Perturb block faces implicitly affected by explicit edge changes:

      DO M = 1, 6

         IF (IFACEPTB(M) == 1) THEN

            I1 = 1   ! WARPQ3DM expects one pair of indices to be
            I2 = IL  ! equal to define the face
            J1 = 1
            J2 = JL
            K1 = 1
            K2 = KL

            IF (M == 1) THEN
               I2 = 1
            ELSE IF (M == 2) THEN
               I1 = IL
            ELSE IF (M == 3) THEN
               J2 = 1
            ELSE IF (M == 4) THEN
               J1 = JL
            ELSE IF (M == 5) THEN
               K2 = 1
            ELSE !  (M == 6)
               K1 = KL
            END IF

            CALL WARPQ3DM_GM (IL, JL, KL, I1, I2, J1, J2, K1, K2,
     >                        XYZ0, S0, DFACEI, DFACEJ, DFACEK, XYZ)
         END IF
      END DO


C     Perturb the volume grid as though all faces have changed.
C     The face perturbations are stored for all three stages so that
C     the volume points can be perturbed in a single pass.


C     Stage 1:  Corner motion (only). This stage is reusable by WARPQ3DM.
C     --------

C     I = 1 and IL intermediate faces.

      CALL DELQ3DM_GM (IL, JL, KL,  1,  1, 1, JL, 1, KL, XYZ0, S0,
     >                 DFACEI(1,1,1,1,1), DFACEJ, DFACEK, XYZ)

      CALL DELQ3DM_GM (IL, JL, KL, IL, IL, 1, JL, 1, KL, XYZ0, S0,
     >                 DFACEI(1,1,1,2,1), DFACEJ, DFACEK, XYZ)

C     J = 1 and JL intermediate faces.

      CALL DELQ3DM_GM (IL, JL, KL, 1, IL,  1,  1, 1, KL, XYZ0, S0,
     >                 DFACEI, DFACEJ(1,1,1,1,1), DFACEK, XYZ)

      CALL DELQ3DM_GM (IL, JL, KL, 1, IL, JL, JL, 1, KL, XYZ0, S0,
     >                 DFACEI, DFACEJ(1,1,1,2,1), DFACEK, XYZ)

C     K = 1 and KL intermediate faces.

      CALL DELQ3DM_GM (IL, JL, KL, 1, IL, 1, JL,  1,  1, XYZ0, S0,
     >                 DFACEI, DFACEJ, DFACEK(1,1,1,1,1), XYZ)

      CALL DELQ3DM_GM (IL, JL, KL, 1, IL, 1, JL, KL, KL, XYZ0, S0,
     >                 DFACEI, DFACEJ, DFACEK(1,1,1,2,1), XYZ)


C     Stage 2:  Handle edge motion from above interim edges to final edges.
C     --------

C     James's insight here: consider the 4 faces affected by the motion of
C     4 edges in a given (index) direction; furthermore, keep the two
C     directions of the resulting face perturbations separated for proper
C     weighted combination during the final interpolation into the interior.
C     The 3 index directions are then independent of each other (added).

C     I = 1 and IL faces:

      L = 1
      DO I = 1, IL, IL - 1

C        K = 1 and KL edge perturbations (J direction edges):

         DO J = 2, JL - 1
            DFACEI(1,J, 1,L,2) = XYZ(1,I,J, 1) - XYZ0(1,I,J, 1) -
     >      DFACEI(1,J, 1,L,1)
            DFACEI(2,J, 1,L,2) = XYZ(2,I,J, 1) - XYZ0(2,I,J, 1) -
     >      DFACEI(2,J, 1,L,1)
            DFACEI(3,J, 1,L,2) = XYZ(3,I,J, 1) - XYZ0(3,I,J, 1) -
     >      DFACEI(3,J, 1,L,1)
            DFACEI(1,J,KL,L,2) = XYZ(1,I,J,KL) - XYZ0(1,I,J,KL) -
     >      DFACEI(1,J,KL,L,1)
            DFACEI(2,J,KL,L,2) = XYZ(2,I,J,KL) - XYZ0(2,I,J,KL) -
     >      DFACEI(2,J,KL,L,1)
            DFACEI(3,J,KL,L,2) = XYZ(3,I,J,KL) - XYZ0(3,I,J,KL) -
     >      DFACEI(3,J,KL,L,1)
         END DO

C        J = 1 and JL edge perturbations (K direction edges):

         DO K = 2, KL - 1
            DFACEI(1, 1,K,L,3) = XYZ(1,I, 1,K) - XYZ0(1,I, 1,K) -
     >      DFACEI(1, 1,K,L,1)
            DFACEI(2, 1,K,L,3) = XYZ(2,I, 1,K) - XYZ0(2,I, 1,K) -
     >      DFACEI(2, 1,K,L,1)
            DFACEI(3, 1,K,L,3) = XYZ(3,I, 1,K) - XYZ0(3,I, 1,K) -
     >      DFACEI(3, 1,K,L,1)
            DFACEI(1,JL,K,L,3) = XYZ(1,I,JL,K) - XYZ0(1,I,JL,K) -
     >      DFACEI(1,JL,K,L,1)
            DFACEI(2,JL,K,L,3) = XYZ(2,I,JL,K) - XYZ0(2,I,JL,K) -
     >      DFACEI(2,JL,K,L,1)
            DFACEI(3,JL,K,L,3) = XYZ(3,I,JL,K) - XYZ0(3,I,JL,K) -
     >      DFACEI(3,JL,K,L,1)
         END DO

C        Interpolate stage 2 interior points for this I face, keeping the
C        two index directions separated.

         DO K = 2, KL - 1
            DO J = 2, JL - 1
               WTJ2 = S0(2,I,J,K)
               WTJ1 = ONE - WTJ2
               WTK2 = S0(3,I,J,K)
               WTK1 = ONE - WTK2

               DFACEI(1,J,K,L,3) = WTJ1*DFACEI(1, 1,K,L,3) +
     >                             WTJ2*DFACEI(1,JL,K,L,3)
               DFACEI(2,J,K,L,3) = WTJ1*DFACEI(2, 1,K,L,3) +
     >                             WTJ2*DFACEI(2,JL,K,L,3)
               DFACEI(3,J,K,L,3) = WTJ1*DFACEI(3, 1,K,L,3) +
     >                             WTJ2*DFACEI(3,JL,K,L,3)

               DFACEI(1,J,K,L,2) = WTK1*DFACEI(1,J, 1,L,2) +
     >                             WTK2*DFACEI(1,J,KL,L,2)
               DFACEI(2,J,K,L,2) = WTK1*DFACEI(2,J, 1,L,2) +
     >                             WTK2*DFACEI(2,J,KL,L,2)
               DFACEI(3,J,K,L,2) = WTK1*DFACEI(3,J, 1,L,2) +
     >                             WTK2*DFACEI(3,J,KL,L,2)
            END DO
         END DO
         L = 2
      END DO

C     J = 1 and JL faces, stage 2:

      L = 1
      DO J = 1, JL, JL - 1

C        K = 1 and KL edge perturbations (I direction):

         DO I = 2, IL - 1
            DFACEJ(1,I, 1,L,2) = XYZ(1,I,J, 1) - XYZ0(1,I,J, 1) -
     >      DFACEJ(1,I, 1,L,1)
            DFACEJ(2,I, 1,L,2) = XYZ(2,I,J, 1) - XYZ0(2,I,J, 1) -
     >      DFACEJ(2,I, 1,L,1)
            DFACEJ(3,I, 1,L,2) = XYZ(3,I,J, 1) - XYZ0(3,I,J, 1) -
     >      DFACEJ(3,I, 1,L,1)
            DFACEJ(1,I,KL,L,2) = XYZ(1,I,J,KL) - XYZ0(1,I,J,KL) -
     >      DFACEJ(1,I,KL,L,1)
            DFACEJ(2,I,KL,L,2) = XYZ(2,I,J,KL) - XYZ0(2,I,J,KL) -
     >      DFACEJ(2,I,KL,L,1)
            DFACEJ(3,I,KL,L,2) = XYZ(3,I,J,KL) - XYZ0(3,I,J,KL) -
     >      DFACEJ(3,I,KL,L,1)
         END DO

C        I = 1 and IL edge perturbations (K direction):

         DO K = 2, KL - 1
            DFACEJ(1, 1,K,L,3) = XYZ(1, 1,J,K) - XYZ0(1, 1,J,K) -
     >      DFACEJ(1, 1,K,L,1)
            DFACEJ(2, 1,K,L,3) = XYZ(2, 1,J,K) - XYZ0(2, 1,J,K) -
     >      DFACEJ(2, 1,K,L,1)
            DFACEJ(3, 1,K,L,3) = XYZ(3, 1,J,K) - XYZ0(3, 1,J,K) -
     >      DFACEJ(3, 1,K,L,1)
            DFACEJ(1,IL,K,L,3) = XYZ(1,IL,J,K) - XYZ0(1,IL,J,K) -
     >      DFACEJ(1,IL,K,L,1)
            DFACEJ(2,IL,K,L,3) = XYZ(2,IL,J,K) - XYZ0(2,IL,J,K) -
     >      DFACEJ(2,IL,K,L,1)
            DFACEJ(3,IL,K,L,3) = XYZ(3,IL,J,K) - XYZ0(3,IL,J,K) -
     >      DFACEJ(3,IL,K,L,1)
         END DO

C        Interpolate stage 2 interior points for this J face.

         DO K = 2, KL - 1
            DO I = 2, IL - 1
               WTI2 = S0(1,I,J,K)
               WTI1 = ONE - WTI2
               WTK2 = S0(3,I,J,K)
               WTK1 = ONE - WTK2

               DFACEJ(1,I,K,L,3) = WTI1*DFACEJ(1, 1,K,L,3) +
     >                             WTI2*DFACEJ(1,IL,K,L,3)
               DFACEJ(2,I,K,L,3) = WTI1*DFACEJ(2, 1,K,L,3) +
     >                             WTI2*DFACEJ(2,IL,K,L,3)
               DFACEJ(3,I,K,L,3) = WTI1*DFACEJ(3, 1,K,L,3) +
     >                             WTI2*DFACEJ(3,IL,K,L,3)

               DFACEJ(1,I,K,L,2) = WTK1*DFACEJ(1,I, 1,L,2) +
     >                             WTK2*DFACEJ(1,I,KL,L,2)
               DFACEJ(2,I,K,L,2) = WTK1*DFACEJ(2,I, 1,L,2) +
     >                             WTK2*DFACEJ(2,I,KL,L,2)
               DFACEJ(3,I,K,L,2) = WTK1*DFACEJ(3,I, 1,L,2) +
     >                             WTK2*DFACEJ(3,I,KL,L,2)
            END DO
         END DO
         L = 2
      END DO

C     K = 1 and KL faces, stage 2:

      L = 1
      DO K = 1, KL, KL - 1

C        J = 1 and JL edge perturbations (I direction):

         DO I = 2, IL - 1
            DFACEK(1,I, 1,L,2) = XYZ(1,I, 1,K) - XYZ0(1,I, 1,K) -
     >      DFACEK(1,I, 1,L,1)
            DFACEK(2,I, 1,L,2) = XYZ(2,I, 1,K) - XYZ0(2,I, 1,K) -
     >      DFACEK(2,I, 1,L,1)
            DFACEK(3,I, 1,L,2) = XYZ(3,I, 1,K) - XYZ0(3,I, 1,K) -
     >      DFACEK(3,I, 1,L,1)
            DFACEK(1,I,JL,L,2) = XYZ(1,I,JL,K) - XYZ0(1,I,JL,K) -
     >      DFACEK(1,I,JL,L,1)
            DFACEK(2,I,JL,L,2) = XYZ(2,I,JL,K) - XYZ0(2,I,JL,K) -
     >      DFACEK(2,I,JL,L,1)
            DFACEK(3,I,JL,L,2) = XYZ(3,I,JL,K) - XYZ0(3,I,JL,K) -
     >      DFACEK(3,I,JL,L,1)
         END DO

C        I = 1 and IL edge perturbations (J direction):

         DO J = 2, JL - 1
            DFACEK(1, 1,J,L,3) = XYZ(1, 1,J,K) - XYZ0(1, 1,J,K) -
     >      DFACEK(1, 1,J,L,1)
            DFACEK(2, 1,J,L,3) = XYZ(2, 1,J,K) - XYZ0(2, 1,J,K) -
     >      DFACEK(2, 1,J,L,1)
            DFACEK(3, 1,J,L,3) = XYZ(3, 1,J,K) - XYZ0(3, 1,J,K) -
     >      DFACEK(3, 1,J,L,1)
            DFACEK(1,IL,J,L,3) = XYZ(1,IL,J,K) - XYZ0(1,IL,J,K) -
     >      DFACEK(1,IL,J,L,1)
            DFACEK(2,IL,J,L,3) = XYZ(2,IL,J,K) - XYZ0(2,IL,J,K) -
     >      DFACEK(2,IL,J,L,1)
            DFACEK(3,IL,J,L,3) = XYZ(3,IL,J,K) - XYZ0(3,IL,J,K) -
     >      DFACEK(3,IL,J,L,1)
         END DO

C        Interpolate stage 2 interior points for this K face.

         DO J = 2, JL - 1
            DO I = 2, IL - 1
               WTI2 = S0(1,I,J,K)
               WTI1 = ONE - WTI2
               WTJ2 = S0(2,I,J,K)
               WTJ1 = ONE - WTJ2

               DFACEK(1,I,J,L,3) = WTI1*DFACEK(1, 1,J,L,3) +
     >                             WTI2*DFACEK(1,IL,J,L,3)
               DFACEK(2,I,J,L,3) = WTI1*DFACEK(2, 1,J,L,3) +
     >                             WTI2*DFACEK(2,IL,J,L,3)
               DFACEK(3,I,J,L,3) = WTI1*DFACEK(3, 1,J,L,3) +
     >                             WTI2*DFACEK(3,IL,J,L,3)

               DFACEK(1,I,J,L,2) = WTJ1*DFACEK(1,I, 1,L,2) +
     >                             WTJ2*DFACEK(1,I,JL,L,2)
               DFACEK(2,I,J,L,2) = WTJ1*DFACEK(2,I, 1,L,2) +
     >                             WTJ2*DFACEK(2,I,JL,L,2)
               DFACEK(3,I,J,L,2) = WTJ1*DFACEK(3,I, 1,L,2) +
     >                             WTJ2*DFACEK(3,I,JL,L,2)
            END DO
         END DO
         L = 2
      END DO


C     Stage 3:  Handle face motion from above interim faces to final faces.
C     --------

      DO K = 2, KL - 1
         DO J = 2, JL - 1
            DFACEI(1,J,K,1,4) = XYZ(1, 1,J,K) - XYZ0(1, 1,J,K) -
     >         DFACEI(1,J,K,1,1) - DFACEI(1,J,K,1,2) - DFACEI(1,J,K,1,3)
            DFACEI(2,J,K,1,4) = XYZ(2, 1,J,K) - XYZ0(2, 1,J,K) -
     >         DFACEI(2,J,K,1,1) - DFACEI(2,J,K,1,2) - DFACEI(2,J,K,1,3)
            DFACEI(3,J,K,1,4) = XYZ(3, 1,J,K) - XYZ0(3, 1,J,K) -
     >         DFACEI(3,J,K,1,1) - DFACEI(3,J,K,1,2) - DFACEI(3,J,K,1,3)
            DFACEI(1,J,K,2,4) = XYZ(1,IL,J,K) - XYZ0(1,IL,J,K) -
     >         DFACEI(1,J,K,2,1) - DFACEI(1,J,K,2,2) - DFACEI(1,J,K,2,3)
            DFACEI(2,J,K,2,4) = XYZ(2,IL,J,K) - XYZ0(2,IL,J,K) -
     >         DFACEI(2,J,K,2,1) - DFACEI(2,J,K,2,2) - DFACEI(2,J,K,2,3)
            DFACEI(3,J,K,2,4) = XYZ(3,IL,J,K) - XYZ0(3,IL,J,K) -
     >         DFACEI(3,J,K,2,1) - DFACEI(3,J,K,2,2) - DFACEI(3,J,K,2,3)
         END DO

         DO I = 2, IL - 1
            DFACEJ(1,I,K,1,4) = XYZ(1,I, 1,K) - XYZ0(1,I, 1,K) -
     >         DFACEJ(1,I,K,1,1) - DFACEJ(1,I,K,1,2) - DFACEJ(1,I,K,1,3)
            DFACEJ(2,I,K,1,4) = XYZ(2,I, 1,K) - XYZ0(2,I, 1,K) -
     >         DFACEJ(2,I,K,1,1) - DFACEJ(2,I,K,1,2) - DFACEJ(2,I,K,1,3)
            DFACEJ(3,I,K,1,4) = XYZ(3,I, 1,K) - XYZ0(3,I, 1,K) -
     >         DFACEJ(3,I,K,1,1) - DFACEJ(3,I,K,1,2) - DFACEJ(3,I,K,1,3)
            DFACEJ(1,I,K,2,4) = XYZ(1,I,JL,K) - XYZ0(1,I,JL,K) -
     >         DFACEJ(1,I,K,2,1) - DFACEJ(1,I,K,2,2) - DFACEJ(1,I,K,2,3)
            DFACEJ(2,I,K,2,4) = XYZ(2,I,JL,K) - XYZ0(2,I,JL,K) -
     >         DFACEJ(2,I,K,2,1) - DFACEJ(2,I,K,2,2) - DFACEJ(2,I,K,2,3)
            DFACEJ(3,I,K,2,4) = XYZ(3,I,JL,K) - XYZ0(3,I,JL,K) -
     >         DFACEJ(3,I,K,2,1) - DFACEJ(3,I,K,2,2) - DFACEJ(3,I,K,2,3)
         END DO
      END DO

      DO J = 2, JL - 1
         DO I = 2, IL - 1
            DFACEK(1,I,J,1,4) = XYZ(1,I,J, 1) - XYZ0(1,I,J, 1) -
     >         DFACEK(1,I,J,1,1) - DFACEK(1,I,J,1,2) - DFACEK(1,I,J,1,3)
            DFACEK(2,I,J,1,4) = XYZ(2,I,J, 1) - XYZ0(2,I,J, 1) -
     >         DFACEK(2,I,J,1,1) - DFACEK(2,I,J,1,2) - DFACEK(2,I,J,1,3)
            DFACEK(3,I,J,1,4) = XYZ(3,I,J, 1) - XYZ0(3,I,J, 1) -
     >         DFACEK(3,I,J,1,1) - DFACEK(3,I,J,1,2) - DFACEK(3,I,J,1,3)
            DFACEK(1,I,J,2,4) = XYZ(1,I,J,KL) - XYZ0(1,I,J,KL) -
     >         DFACEK(1,I,J,2,1) - DFACEK(1,I,J,2,2) - DFACEK(1,I,J,2,3)
            DFACEK(2,I,J,2,4) = XYZ(2,I,J,KL) - XYZ0(2,I,J,KL) -
     >         DFACEK(2,I,J,2,1) - DFACEK(2,I,J,2,2) - DFACEK(2,I,J,2,3)
            DFACEK(3,I,J,2,4) = XYZ(3,I,J,KL) - XYZ0(3,I,J,KL) -
     >         DFACEK(3,I,J,2,1) - DFACEK(3,I,J,2,2) - DFACEK(3,I,J,2,3)
         END DO
      END DO


C     Perturb the interior volume points.
C     All stages are performed at once via interpolation from the face
C     perturbations stored for each stage.  Note that the three stages
C     accumulate the contributions from the three subscript directions
C     with varying degrees of independence.

      DO K = 2, KL - 1
         DO J = 2, JL - 1
            DO I = 2, IL - 1
               WTI2 = S0(1,I,J,K)
               WTI1 = ONE - WTI2
               WTJ2 = S0(2,I,J,K)
               WTJ1 = ONE - WTJ2
               WTK2 = S0(3,I,J,K)
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

               XYZ(1,I,J,K) = XYZ0(1,I,J,K) + DEL1 +DEL2 +DELI+DELJ+DELK

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

               XYZ(2,I,J,K) = XYZ0(2,I,J,K) + DEL1 +DEL2 +DELI+DELJ+DELK

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

               XYZ(3,I,J,K) = XYZ0(3,I,J,K) + DEL1 +DEL2 +DELI+DELJ+DELK
            END DO
         END DO
      END DO

      END SUBROUTINE WARPBLK_GM

C+--------------------------------------------------------------------------

      SUBROUTINE PARAM3DM_GM (IL, JL, KL, XYZ, S)

C     ******************************************************************
C     *   PARAM3DM parameterizes the volume of one block of a multi-   *
C     *   block grid structure by setting up the normalized arc-length *
C     *   increments in all three index directions.                    *
C     *                                                                *
C     *   11/29/95  D.Saunders  Adaptation of PARAMXYZ for specialized *
C     *                         WARP-BLK used by FLO107-MB.            *
C     *   06/19/96      "       Allow for degenerate edges.            *
C     *                                                                *
C     *   10/99    M.Rimlinger  Removed halos for GMORPH/MESHWARP app. *
C     *                                                                *
C     *   David Saunders/James Reuther, NASA Ames Research Center, CA. *
C     ******************************************************************

      IMPLICIT NONE

C     Arguments:

      INTEGER   IL, JL, KL      ! I Grid array dimensions.
      REAL      XYZ(3,IL,JL,KL) ! I Grid coordinates
      REAL      S(3,IL,JL,KL)   ! O Normalized arc-lengths:
                                !   S(1,1,J,K) = 0.,
                                !   S(2,I,1,K) = 0.,
                                !   S(3,I,J,1) = 0.,
                                !   S(1,IL,J,K) = 1.,etc.
C     Local constants:

      REAL, PARAMETER :: ONE = 1.E+0, ZERO = 0.E+0

C     Local variables:

      INTEGER   I, J, K

C     Local functions:

      REAL      DELI, DELJ, DELK

      DELI(I,J,K) = SQRT ((XYZ(1,I,J,K) - XYZ(1,I-1,J,K)) ** 2 +
     >                    (XYZ(2,I,J,K) - XYZ(2,I-1,J,K)) ** 2 +
     >                    (XYZ(3,I,J,K) - XYZ(3,I-1,J,K)) ** 2)

      DELJ(I,J,K) = SQRT ((XYZ(1,I,J,K) - XYZ(1,I,J-1,K)) ** 2 +
     >                    (XYZ(2,I,J,K) - XYZ(2,I,J-1,K)) ** 2 +
     >                    (XYZ(3,I,J,K) - XYZ(3,I,J-1,K)) ** 2)

      DELK(I,J,K) = SQRT ((XYZ(1,I,J,K) - XYZ(1,I,J,K-1)) ** 2 +
     >                    (XYZ(2,I,J,K) - XYZ(2,I,J,K-1)) ** 2 +
     >                    (XYZ(3,I,J,K) - XYZ(3,I,J,K-1)) ** 2)

C     Execution:
C     ----------

C     Zero the three low-end faces (or edges if one plane is specified).

      DO K = 1, KL
         DO J = 1, JL
            S(1,1,J,K) = ZERO
         END DO

         DO I = 1, IL
            S(2,I,1,K) = ZERO
         END DO
      END DO

      DO J = 1, JL
         DO I = 1, IL
            S(3,I,J,1) = ZERO
         END DO
      END DO

C     Set up the low-end edge lines because they are missed by the
C     following loops over most of the low-end faces:

      DO I = 2, IL
         S(1,I,1,1) = S(1,I-1,1,1) + DELI(I,1,1)
      END DO

      DO J = 2, JL
         S(2,1,J,1) = S(2,1,J-1,1) + DELJ(1,J,1)
      END DO

      DO K = 2, KL
         S(3,1,1,K) = S(3,1,1,K-1) + DELK(1,1,K)
      END DO

C     Set up the rest of the low-end face lines because they are
C     missed by the the main loop over most of the volume.

      DO K = 2, KL
         DO J = 2, JL
            S(2,1,J,K) = S(2,1,J-1,K) + DELJ(1,J,K)
            S(3,1,J,K) = S(3,1,J,K-1) + DELK(1,J,K)
         END DO
         DO I = 2, IL
            S(1,I,1,K) = S(1,I-1,1,K) + DELI(I,1,K)
            S(3,I,1,K) = S(3,I,1,K-1) + DELK(I,1,K)
         END DO
      END DO

      DO J = 2, JL
         DO I = 2, IL
            S(1,I,J,1) = S(1,I-1,J,1) + DELI(I,J,1)
            S(2,I,J,1) = S(2,I,J-1,1) + DELJ(I,J,1)
         END DO
      END DO

C     Traverse the block just once for all lines except those within
C     the low-end faces.

      DO K = 2, KL
         DO J = 2, JL
            DO I = 2, IL
               S(1,I,J,K) = S(1,I-1,J,K) + DELI(I,J,K)
               S(2,I,J,K) = S(2,I,J-1,K) + DELJ(I,J,K)
               S(3,I,J,K) = S(3,I,J,K-1) + DELK(I,J,K)
            END DO
         END DO
      END DO

C     Normalizing requires another pass through the volume.
C     Handle lines of zero length first by inserting uniform
C     distributions.  Then the standard normalization can be
C     applied safely everywhere.

      DO K = 1, KL

C        Zero-length lines in the I direction?

         DO J = 1, JL
            IF (S(1,IL,J,K) == ZERO) THEN
               DO I = 2, IL
                  S(1,I,J,K) = I - 1
               END DO
            END IF
         END DO

C        Zero-length lines in the J direction?

         DO I = 1, IL
            IF (S(2,I,JL,K) == ZERO) THEN
               DO J = 2, JL
                  S(2,I,J,K) = J - 1
               END DO
            END IF
         END DO
      END DO

C     Zero-length lines in the K direction?

      DO J = 1, JL
         DO I = 1, IL
            IF (S(3,I,J,KL) == ZERO) THEN
               DO K = 2, KL
                  S(3,I,J,K) = K - 1
               END DO
            END IF
         END DO
      END DO

C     Normalize:

      DO K = 1, KL
         DO J = 1, JL
            DO I = 1, IL
               S(1,I,J,K) = S(1,I,J,K) / S(1,IL,J,K)
               S(2,I,J,K) = S(2,I,J,K) / S(2,I,JL,K)
               S(3,I,J,K) = S(3,I,J,K) / S(3,I,J,KL)
            END DO
         END DO
      END DO

C     Finally, precise 1s for the three high-end faces:

      DO K = 1, KL
         DO J = 1, JL
            S(1,IL,J,K) = ONE
         END DO

         DO I = 1, IL
            S(2,I,JL,K) = ONE
         END DO
      END DO

      DO J = 1, JL
         DO I = 1, IL
            S(3,I,J,KL) = ONE
         END DO
      END DO

      END SUBROUTINE PARAM3DM_GM

C+--------------------------------------------------------------------------

      SUBROUTINE WARPQ3DM_GM (IL, JL, KL, I1, I2, J1, J2, K1, K2,
     >                        XYZ0, S0, DFACEI, DFACEJ, DFACEK, XYZ)

C     ******************************************************************
C     *   WARPQ3DM perturbs the interior of one face of one block of a *
C     *   multiblock grid structure given perturbed edges of that face,*
C     *   which is indicated by one pair of equal index arguments.     *
C     *   (E.g.: I2 = I1 means a face in the J/K subspace.)            *
C     *                                                                *
C     *   The two-stage algorithm uses an intermediate perturbation to *
C     *   account for any corner motion, involving edges derived from  *
C     *   the original edges, then a second perturbation to account    *
C     *   for any differences between the intermediate edges and the   *
C     *   specified new edges.                                         *
C     *                                                                *
C     *   The perturbed edges should be input as edges of the desired  *
C     *   output face.  The original relative arc-length increments in *
C     *   each index direction should also be input.  See PARAM3DM for *
C     *   setting them up in preparation for multiple perturbations.   *
C     *                                                                *
C     *   11/29/95  D.Saunders  Adaptation of WARPQ3D for specialized  *
C     *                         WARP-BLK used by FLO107-MB.            *
C     *   04/04/96      "       DELQ3DM does stage 1 only now (all     *
C     *                         that WARP-BLK needs).                  *
C     *   10/99    M.Rimlinger  Removed halos for GMORPH/MESHWARP.     *
C     *                                                                *
C     *   David Saunders/James Reuther, NASA Ames Research Center, CA. *
C     ******************************************************************

      IMPLICIT NONE

C     Arguments:

      INTEGER   IL, JL, KL             ! I  Grid array dimensions.
      INTEGER   I1, I2, J1, J2, K1, K2 ! I  Define active face,
                                       !    one pair being equal.
      REAL      XYZ0(3,IL,JL,KL)       ! I  Base face coordinates in
                                       !    appropriate places
      REAL      S0(3,IL,JL,KL)         ! I  Base normalized arc-lengths
      REAL      DFACEI(3,JL,KL)        ! S  For face perturbations; e.g.,
      REAL      DFACEJ(3,IL,KL)        !    DFACEI(1:3,1:JL,1:KL) =
      REAL      DFACEK(3,IL,JL)        !    dX, dY, dZ for an I face, etc.
      REAL      XYZ(3,IL,JL,KL)        !I/O Grid coordinates: new edges of
                                       !    a face in; full new face out
C     Local constants:

      REAL, PARAMETER :: ONE = 1.E+0

C     Local variables:

      INTEGER   I, J, K
      REAL      DELI, DELJ, DELK, WTI1, WTI2, WTJ1, WTJ2, WTK1, WTK2

C     Execution:
C     ----------

C     Stage 1:
C     Handle any corner motion by generating an intermediate face with
C     the final corners but otherwise derived from the original edges.
C     Actually, just set up the appropriate face perturbations.

      CALL DELQ3DM_GM (IL, JL, KL, I1, I2, J1, J2, K1, K2, XYZ0, S0,
     >                 DFACEI, DFACEJ, DFACEK, XYZ)


C     Stage 2:
C     Set up the perturbations from the intermediate edges to the final
C     edges, then interpolate them into the interior points.

      IF (I1 == I2) THEN          ! I plane case:

         I = I1

C        J = 1 and JL edge perturbations:

         DO K = 2, KL - 1
            DFACEI(1, 1,K) = XYZ(1,I, 1,K)-XYZ0(1,I, 1,K)-DFACEI(1, 1,K)
            DFACEI(2, 1,K) = XYZ(2,I, 1,K)-XYZ0(2,I, 1,K)-DFACEI(2, 1,K)
            DFACEI(3, 1,K) = XYZ(3,I, 1,K)-XYZ0(3,I, 1,K)-DFACEI(3, 1,K)
            DFACEI(1,JL,K) = XYZ(1,I,JL,K)-XYZ0(1,I,JL,K)-DFACEI(1,JL,K)
            DFACEI(2,JL,K) = XYZ(2,I,JL,K)-XYZ0(2,I,JL,K)-DFACEI(2,JL,K)
            DFACEI(3,JL,K) = XYZ(3,I,JL,K)-XYZ0(3,I,JL,K)-DFACEI(3,JL,K)
         END DO

C        K = 1 and KL edge perturbations:

         DO J = 2, JL - 1
            DFACEI(1,J, 1) = XYZ(1,I,J, 1)-XYZ0(1,I,J, 1)-DFACEI(1,J, 1)
            DFACEI(2,J, 1) = XYZ(2,I,J, 1)-XYZ0(2,I,J, 1)-DFACEI(2,J, 1)
            DFACEI(3,J, 1) = XYZ(3,I,J, 1)-XYZ0(3,I,J, 1)-DFACEI(3,J, 1)
            DFACEI(1,J,KL) = XYZ(1,I,J,KL)-XYZ0(1,I,J,KL)-DFACEI(1,J,KL)
            DFACEI(2,J,KL) = XYZ(2,I,J,KL)-XYZ0(2,I,J,KL)-DFACEI(2,J,KL)
            DFACEI(3,J,KL) = XYZ(3,I,J,KL)-XYZ0(3,I,J,KL)-DFACEI(3,J,KL)
         END DO

C        Interior points: accumulate the (independent) contributions.

         DO K = 2, KL - 1
            DO J = 2, JL - 1
               WTJ2 = S0(2,I,J,K)
               WTJ1 = ONE - WTJ2
               WTK2 = S0(3,I,J,K)
               WTK1 = ONE - WTK2

               DELJ = WTJ1 * DFACEI(1,1,K) + WTJ2 * DFACEI(1,JL,K)
               DELK = WTK1 * DFACEI(1,J,1) + WTK2 * DFACEI(1,J,KL)

               XYZ(1,I,J,K) = XYZ0(1,I,J,K) + DFACEI(1,J,K) + DELJ +DELK

               DELJ = WTJ1 * DFACEI(2,1,K) + WTJ2 * DFACEI(2,JL,K)
               DELK = WTK1 * DFACEI(2,J,1) + WTK2 * DFACEI(2,J,KL)

               XYZ(2,I,J,K) = XYZ0(2,I,J,K) + DFACEI(2,J,K) + DELJ +DELK

               DELJ = WTJ1 * DFACEI(3,1,K) + WTJ2 * DFACEI(3,JL,K)
               DELK = WTK1 * DFACEI(3,J,1) + WTK2 * DFACEI(3,J,KL)

               XYZ(3,I,J,K) = XYZ0(3,I,J,K) + DFACEI(3,J,K) + DELJ +DELK
            END DO
         END DO

      ELSE IF (J1 == J2) THEN     ! J plane case:

         J = J1

C        I = 1 and IL edge perturbations:

         DO K = 2, KL - 1
            DFACEJ(1, 1,K) = XYZ(1, 1,J,K)-XYZ0(1, 1,J,K)-DFACEJ(1, 1,K)
            DFACEJ(2, 1,K) = XYZ(2, 1,J,K)-XYZ0(2, 1,J,K)-DFACEJ(2, 1,K)
            DFACEJ(3, 1,K) = XYZ(3, 1,J,K)-XYZ0(3, 1,J,K)-DFACEJ(3, 1,K)
            DFACEJ(1,IL,K) = XYZ(1,IL,J,K)-XYZ0(1,IL,J,K)-DFACEJ(1,IL,K)
            DFACEJ(2,IL,K) = XYZ(2,IL,J,K)-XYZ0(2,IL,J,K)-DFACEJ(2,IL,K)
            DFACEJ(3,IL,K) = XYZ(3,IL,J,K)-XYZ0(3,IL,J,K)-DFACEJ(3,IL,K)
         END DO

C        K = 1 and KL edge perturbations:

         DO I = 2, IL - 1
            DFACEJ(1,I, 1) = XYZ(1,I,J, 1)-XYZ0(1,I,J, 1)-DFACEJ(1,I, 1)
            DFACEJ(2,I, 1) = XYZ(2,I,J, 1)-XYZ0(2,I,J, 1)-DFACEJ(2,I, 1)
            DFACEJ(3,I, 1) = XYZ(3,I,J, 1)-XYZ0(3,I,J, 1)-DFACEJ(3,I, 1)
            DFACEJ(1,I,KL) = XYZ(1,I,J,KL)-XYZ0(1,I,J,KL)-DFACEJ(1,I,KL)
            DFACEJ(2,I,KL) = XYZ(2,I,J,KL)-XYZ0(2,I,J,KL)-DFACEJ(2,I,KL)
            DFACEJ(3,I,KL) = XYZ(3,I,J,KL)-XYZ0(3,I,J,KL)-DFACEJ(3,I,KL)
         END DO

C        Interior points: accumulate the (independent) contributions.

         DO K = 2, KL - 1
            DO I = 2, IL - 1
               WTI2 = S0(1,I,J,K)
               WTI1 = ONE - WTI2
               WTK2 = S0(3,I,J,K)
               WTK1 = ONE - WTK2

               DELI = WTI1 * DFACEJ(1,1,K) + WTI2 * DFACEJ(1,IL,K)
               DELK = WTK1 * DFACEJ(1,I,1) + WTK2 * DFACEJ(1,I,KL)

               XYZ(1,I,J,K) = XYZ0(1,I,J,K) + DFACEJ(1,I,K) + DELI +DELK

               DELI = WTI1 * DFACEJ(2,1,K) + WTI2 * DFACEJ(2,IL,K)
               DELK = WTK1 * DFACEJ(2,I,1) + WTK2 * DFACEJ(2,I,KL)

               XYZ(2,I,J,K) = XYZ0(2,I,J,K) + DFACEJ(2,I,K) + DELI +DELK

               DELI = WTI1 * DFACEJ(3,1,K) + WTI2 * DFACEJ(3,IL,K)
               DELK = WTK1 * DFACEJ(3,I,1) + WTK2 * DFACEJ(3,I,KL)

               XYZ(3,I,J,K) = XYZ0(3,I,J,K) + DFACEJ(3,I,K) + DELI +DELK
            END DO
         END DO

      ELSE IF (K1 == K2) THEN     ! K plane case:

         K = K1

C        I = 1 and IL edge perturbations:

         DO J = 2, JL - 1
            DFACEK(1, 1,J) = XYZ(1, 1,J,K)-XYZ0(1, 1,J,K)-DFACEK(1, 1,J)
            DFACEK(2, 1,J) = XYZ(2, 1,J,K)-XYZ0(2, 1,J,K)-DFACEK(2, 1,J)
            DFACEK(3, 1,J) = XYZ(3, 1,J,K)-XYZ0(3, 1,J,K)-DFACEK(3, 1,J)
            DFACEK(1,IL,J) = XYZ(1,IL,J,K)-XYZ0(1,IL,J,K)-DFACEK(1,IL,J)
            DFACEK(2,IL,J) = XYZ(2,IL,J,K)-XYZ0(2,IL,J,K)-DFACEK(2,IL,J)
            DFACEK(3,IL,J) = XYZ(3,IL,J,K)-XYZ0(3,IL,J,K)-DFACEK(3,IL,J)
         END DO

C        J = 1 and JL edge perturbations:

         DO I = 2, IL - 1
            DFACEK(1,I, 1) = XYZ(1,I, 1,K)-XYZ0(1,I, 1,K)-DFACEK(1,I, 1)
            DFACEK(2,I, 1) = XYZ(2,I, 1,K)-XYZ0(2,I, 1,K)-DFACEK(2,I, 1)
            DFACEK(3,I, 1) = XYZ(3,I, 1,K)-XYZ0(3,I, 1,K)-DFACEK(3,I, 1)
            DFACEK(1,I,JL) = XYZ(1,I,JL,K)-XYZ0(1,I,JL,K)-DFACEK(1,I,JL)
            DFACEK(2,I,JL) = XYZ(2,I,JL,K)-XYZ0(2,I,JL,K)-DFACEK(2,I,JL)
            DFACEK(3,I,JL) = XYZ(3,I,JL,K)-XYZ0(3,I,JL,K)-DFACEK(3,I,JL)
         END DO

C        Interior points: accumulate the (independent) contributions.

         DO J = 2, JL - 1
            DO I = 2, IL - 1
               WTI2 = S0(1,I,J,K)
               WTI1 = ONE - WTI2
               WTJ2 = S0(2,I,J,K)
               WTJ1 = ONE - WTJ2

               DELI = WTI1 * DFACEK(1,1,J) + WTI2 * DFACEK(1,IL,J)
               DELJ = WTJ1 * DFACEK(1,I,1) + WTJ2 * DFACEK(1,I,JL)

               XYZ(1,I,J,K) = XYZ0(1,I,J,K) + DFACEK(1,I,J) + DELI +DELJ

               DELI = WTI1 * DFACEK(2,1,J) + WTI2 * DFACEK(2,IL,J)
               DELJ = WTJ1 * DFACEK(2,I,1) + WTJ2 * DFACEK(2,I,JL)

               XYZ(2,I,J,K) = XYZ0(2,I,J,K) + DFACEK(2,I,J) + DELI +DELJ

               DELI = WTI1 * DFACEK(3,1,J) + WTI2 * DFACEK(3,IL,J)
               DELJ = WTJ1 * DFACEK(3,I,1) + WTJ2 * DFACEK(3,I,JL)

               XYZ(3,I,J,K) = XYZ0(3,I,J,K) + DFACEK(3,I,J) + DELI +DELJ
            END DO
         END DO

      END IF

      END SUBROUTINE WARPQ3DM_GM

C+--------------------------------------------------------------------------

      SUBROUTINE DELQ3DM_GM (IL, JL, KL, I1, I2, J1, J2, K1, K2,
     >                       XYZ0, S0,
     >                       DFACEI, DFACEJ, DFACEK, XYZ)

C     ******************************************************************
C     *   DELQ3DM performs stage 1 of the WARPQ3DM 3-space surface     *
C     *   grid perturbation in a form which is reusable by WARP-BLK    *
C     *   It returns face perturbations rather than perturbed face     *
C     *   coordinates.  The three cases of a block face are handled    *
C     *   here by three similar code sections.  Special handling of    *
C     *   fixed corners is avoided to keep the bulk down.              *
C     *                                                                *
C     *   11/29/95  D.Saunders  Adaptation of DELQ3D for specialized   *
C     *                         WARP-BLK and WARPQ3DM used by          *
C     *                         FLO107-MB.                             *
C     *   04/04/96      "       DELQ3DM does only stage 1 now.         *
C     *   10/99    M.Rimlinger  Removed halo cells for GMORPH/MESHWARP.*
C     *                                                                *
C     *   David Saunders/James Reuther, NASA Ames Research Center, CA. *
C     ******************************************************************

      IMPLICIT NONE

C     Arguments:

      INTEGER   IL, JL, KL             ! I Grid array dimensions.
      INTEGER   I1, I2, J1, J2, K1, K2 ! I Define active face,
                                       !   one pair being equal.
      REAL      XYZ0(3,IL,JL,KL)       ! I Base face coordinates in
                                       !   appropriate places
      REAL      S0(3,IL,JL,KL)         ! I Base normalized arc-lengths
      REAL      DFACEI(3,JL,KL)        ! O Reqd. face perturbations:
      REAL      DFACEJ(3,IL,KL)        !   DFACEI(1:3,1:JL,1:KL) =
      REAL      DFACEK(3,IL,JL)        !   dX, dY, dZ for an I face, etc.
      REAL      XYZ(3,IL,JL,KL)        ! I Grid coordinates: new edges of a
                                       !   face input; unchanged on output
C     Local constants:

      REAL, PARAMETER ::
     >          EPS = 1.E-14,          ! EPS safeguards a divide by zero -
     >          ONE = 1.E+0            ! presumably only if result is zero.

C     Local variables:

      INTEGER   I, J, K
      REAL      DELI, DELJ, DELK, WTI1, WTI2, WTJ1, WTJ2, WTK1, WTK2

C     Execution:
C     ----------

      IF (I1 == I2) THEN

C        I plane case:
C        -------------

         I = I1

C        Set up the corner perturbations:

         DO K = 1, KL, KL - 1
            DO J = 1, JL, JL - 1
               DFACEI(1,J,K) = XYZ(1,I,J,K) - XYZ0(1,I,J,K)
               DFACEI(2,J,K) = XYZ(2,I,J,K) - XYZ0(2,I,J,K)
               DFACEI(3,J,K) = XYZ(3,I,J,K) - XYZ0(3,I,J,K)
            END DO
         END DO

C        Set up intermediate edge perturbations corresponding to the
C        final corners but otherwise derived from the original edges.

         DO J = 1, JL, JL - 1
            DO K = 2, KL - 1
               WTK2 = S0(3,I,J,K)
               WTK1 = ONE - WTK2
               DFACEI(1,J,K) = WTK1 * DFACEI(1,J, 1) +
     >                         WTK2 * DFACEI(1,J,KL)
               DFACEI(2,J,K) = WTK1 * DFACEI(2,J, 1) +
     >                         WTK2 * DFACEI(2,J,KL)
               DFACEI(3,J,K) = WTK1 * DFACEI(3,J, 1) +
     >                         WTK2 * DFACEI(3,J,KL)
            END DO
         END DO

         DO K = 1, KL, KL - 1
            DO J = 2, JL - 1
               WTJ2 = S0(2,I,J,K)
               WTJ1 = ONE - WTJ2
               DFACEI(1,J,K) = WTJ1 * DFACEI(1, 1,K) +
     >                         WTJ2 * DFACEI(1,JL,K)
               DFACEI(2,J,K) = WTJ1 * DFACEI(2, 1,K) +
     >                         WTJ2 * DFACEI(2,JL,K)
               DFACEI(3,J,K) = WTJ1 * DFACEI(3, 1,K) +
     >                         WTJ2 * DFACEI(3,JL,K)
            END DO
         END DO

C        Interpolate the intermediate perturbations of interior points.
C        The contributions from each pair of edges are not independent.

         DO K = 2, KL - 1
            DO J = 2, JL - 1
               WTJ2 = S0(2,I,J,K)
               WTJ1 = ONE - WTJ2
               WTK2 = S0(3,I,J,K)
               WTK1 = ONE - WTK2

               DELJ = WTJ1 * DFACEI(1,1,K) + WTJ2 * DFACEI(1,JL,K)
               DELK = WTK1 * DFACEI(1,J,1) + WTK2 * DFACEI(1,J,KL)

               DFACEI(1,J,K) = (ABS (DELJ) * DELJ + ABS (DELK) * DELK) /
     >                          MAX (ABS (DELJ) + ABS (DELK), EPS)

               DELJ = WTJ1 * DFACEI(2,1,K) + WTJ2 * DFACEI(2,JL,K)
               DELK = WTK1 * DFACEI(2,J,1) + WTK2 * DFACEI(2,J,KL)

               DFACEI(2,J,K) = (ABS (DELJ) * DELJ + ABS (DELK) * DELK) /
     >                          MAX (ABS (DELJ) + ABS (DELK), EPS)

               DELJ = WTJ1 * DFACEI(3,1,K) + WTJ2 * DFACEI(3,JL,K)
               DELK = WTK1 * DFACEI(3,J,1) + WTK2 * DFACEI(3,J,KL)

               DFACEI(3,J,K) = (ABS (DELJ) * DELJ + ABS (DELK) * DELK) /
     >                          MAX (ABS (DELJ) + ABS (DELK), EPS)
            END DO
         END DO

      ELSE IF (J1 == J2) THEN

C        J plane case:
C        -------------

         J = J1

C        Corner perturbations:

         DO K = 1, KL, KL - 1
            DO I = 1, IL, IL - 1
               DFACEJ(1,I,K) = XYZ(1,I,J,K) - XYZ0(1,I,J,K)
               DFACEJ(2,I,K) = XYZ(2,I,J,K) - XYZ0(2,I,J,K)
               DFACEJ(3,I,K) = XYZ(3,I,J,K) - XYZ0(3,I,J,K)
            END DO
         END DO

C        Intermediate edge perturbations:

         DO I = 1, IL, IL - 1
            DO K = 2, KL - 1
               WTK2 = S0(3,I,J,K)
               WTK1 = ONE - WTK2
               DFACEJ(1,I,K) = WTK1 * DFACEJ(1,I, 1) +
     >                         WTK2 * DFACEJ(1,I,KL)
               DFACEJ(2,I,K) = WTK1 * DFACEJ(2,I, 1) +
     >                         WTK2 * DFACEJ(2,I,KL)
               DFACEJ(3,I,K) = WTK1 * DFACEJ(3,I, 1) +
     >                         WTK2 * DFACEJ(3,I,KL)
            END DO
         END DO

         DO K = 1, KL, KL - 1
            DO I = 2, IL - 1
               WTI2 = S0(1,I,J,K)
               WTI1 = ONE - WTI2
               DFACEJ(1,I,K) = WTI1 * DFACEJ(1, 1,K) +
     >                         WTI2 * DFACEJ(1,IL,K)
               DFACEJ(2,I,K) = WTI1 * DFACEJ(2, 1,K) +
     >                         WTI2 * DFACEJ(2,IL,K)
               DFACEJ(3,I,K) = WTI1 * DFACEJ(3, 1,K) +
     >                         WTI2 * DFACEJ(3,IL,K)
            END DO
         END DO

C        Intermediate perturbations of interior points:

         DO K = 2, KL - 1
            DO I = 2, IL - 1
               WTI2 = S0(1,I,J,K)
               WTI1 = ONE - WTI2
               WTK2 = S0(3,I,J,K)
               WTK1 = ONE - WTK2

               DELI = WTI1 * DFACEJ(1,1,K) + WTI2 * DFACEJ(1,IL,K)
               DELK = WTK1 * DFACEJ(1,I,1) + WTK2 * DFACEJ(1,I,KL)

               DFACEJ(1,I,K) = (ABS (DELI) * DELI + ABS (DELK) * DELK) /
     >                          MAX (ABS (DELI) + ABS (DELK), EPS)

               DELI = WTI1 * DFACEJ(2,1,K) + WTI2 * DFACEJ(2,IL,K)
               DELK = WTK1 * DFACEJ(2,I,1) + WTK2 * DFACEJ(2,I,KL)

               DFACEJ(2,I,K) = (ABS (DELI) * DELI + ABS (DELK) * DELK) /
     >                          MAX (ABS (DELI) + ABS (DELK), EPS)

               DELI = WTI1 * DFACEJ(3,1,K) + WTI2 * DFACEJ(3,IL,K)
               DELK = WTK1 * DFACEJ(3,I,1) + WTK2 * DFACEJ(3,I,KL)

               DFACEJ(3,I,K) = (ABS (DELI) * DELI + ABS (DELK) * DELK) /
     >                          MAX (ABS (DELI) + ABS (DELK), EPS)
            END DO
         END DO

      ELSE IF (K1 == K2) THEN

C        K plane case:
C        -------------

         K = K1

C        Corner perturbations:

         DO J = 1, JL, JL - 1
            DO I = 1, IL, IL - 1
               DFACEK(1,I,J) = XYZ(1,I,J,K) - XYZ0(1,I,J,K)
               DFACEK(2,I,J) = XYZ(2,I,J,K) - XYZ0(2,I,J,K)
               DFACEK(3,I,J) = XYZ(3,I,J,K) - XYZ0(3,I,J,K)
            END DO
         END DO

C        Intermediate edge perturbations:

         DO I = 1, IL, IL - 1
            DO J = 2, JL - 1
               WTJ2 = S0 (2,I,J,K)
               WTJ1 = ONE - WTJ2
               DFACEK(1,I,J) = WTJ1 * DFACEK(1,I, 1) +
     >                         WTJ2 * DFACEK(1,I,JL)
               DFACEK(2,I,J) = WTJ1 * DFACEK(2,I, 1) +
     >                         WTJ2 * DFACEK(2,I,JL)
               DFACEK(3,I,J) = WTJ1 * DFACEK(3,I, 1) +
     >                         WTJ2 * DFACEK(3,I,JL)
            END DO
         END DO

         DO J = 1, JL, JL - 1
            DO I = 2, IL - 1
               WTI2 = S0(1,I,J,K)
               WTI1 = ONE - WTI2
               DFACEK(1,I,J) = WTI1 * DFACEK(1, 1,J) +
     >                         WTI2 * DFACEK(1,IL,J)
               DFACEK(2,I,J) = WTI1 * DFACEK(2, 1,J) +
     >                         WTI2 * DFACEK(2,IL,J)
               DFACEK(3,I,J) = WTI1 * DFACEK(3, 1,J) +
     >                         WTI2 * DFACEK(3,IL,J)
            END DO
         END DO

C        Intermediate perturbations of interior points:

         DO J = 2, JL - 1
            DO I = 2, IL - 1
               WTI2 = S0(1,I,J,K)
               WTI1 = ONE - WTI2
               WTJ2 = S0(2,I,J,K)
               WTJ1 = ONE - WTJ2

               DELI = WTI1 * DFACEK(1,1,J) + WTI2 * DFACEK(1,IL,J)
               DELJ = WTJ1 * DFACEK(1,I,1) + WTJ2 * DFACEK(1,I,JL)

               DFACEK(1,I,J) = (ABS (DELI) * DELI + ABS (DELJ) * DELJ) /
     >                          MAX (ABS (DELI) + ABS (DELJ), EPS)

               DELI = WTI1 * DFACEK(2,1,J) + WTI2 * DFACEK(2,IL,J)
               DELJ = WTJ1 * DFACEK(2,I,1) + WTJ2 * DFACEK(2,I,JL)

               DFACEK(2,I,J) = (ABS (DELI) * DELI + ABS (DELJ) * DELJ) /
     >                          MAX (ABS (DELI) + ABS (DELJ), EPS)

               DELI = WTI1 * DFACEK(3,1,J) + WTI2 * DFACEK(3,IL,J)
               DELJ = WTJ1 * DFACEK(3,I,1) + WTJ2 * DFACEK(3,I,JL)

               DFACEK(3,I,J) = (ABS (DELI) * DELI + ABS (DELJ) * DELJ) /
     >                          MAX (ABS (DELI) + ABS (DELJ), EPS)
            END DO
         END DO

      END IF

      END SUBROUTINE DELQ3DM_GM

C+--------------------------------------------------------------------------

      SUBROUTINE WARPFACE (DXC, DXEJ, DXEK, JL, KL, DFACE, U0, V0)

C     ******************************************************************
C     *   WARPFACE controls the correction of design faces by using    *
C     *   updated corner and edge information.                         *
C     *                                                                *
C     *   ??? What's the history, etc. ???                             *
C     ******************************************************************

      IMPLICIT NONE

C     Arguments:

      INTEGER   JL, KL                     ! I  Face dimensions
      REAL      DXC(3,4)                   ! I  Corner perturbations
      REAL      DXEJ(3,2,JL), DXEK(3,2,KL),
     >          DXCJ(3,2,JL), DXCK(3,2,KL)
      REAL      DFACE(3,JL,KL)             !I/O Face perturbations
      REAL      U0(JL,KL), V0(JL,KL)       ! I  Base normalized arc-lengths;
                                           !    halo allows use of XYZ pointers

C     Local constants:

      REAL, PARAMETER ::
     >          EPS = 1.E-14, ONE = 1.E+0, ! EPS safeguards a divide by zero -
     >          ZERO = 0.E+0               ! presumably only if result is zero.

C     Local variables:

      INTEGER   J, K, LC, LE
      REAL      DELJ, DELK, WTJ1, WTJ2, WTK1, WTK2
      REAL      TXEJ(3,2,JL), TXEK(3,2,KL)

C     Execution:

C     Perturb implicit block edges based on correct corner and edge values.

C     NOTE:  Corner values are assumed to be correct.  Edge values are NOT
C     assumed to be correct, and thus are implicitely corrected before they
C     are added into the face array.

      LC = 1
      LE = 1

      DO K = 1, KL, KL - 1
         DO J = 1, JL

            WTJ2 = U0(J,K)
            WTJ1 = ONE - WTJ2

            DXCJ(1,LE,J) = WTJ1*DXC(1,LC) +
     >                     WTJ2*DXC(1,LC+1) - DFACE(1,J,K)
            DXCJ(2,LE,J) = WTJ1*DXC(2,LC) +
     >                     WTJ2*DXC(2,LC+1) - DFACE(2,J,K)
            DXCJ(3,LE,J) = WTJ1*DXC(3,LC) +
     >                     WTJ2*DXC(3,LC+1) - DFACE(3,J,K)

            TXEJ(1,LE,J) = DXEJ(1,LE,J) -
     >                     (WTJ1*DXEJ(1,LE,1) + WTJ2*DXEJ(1,LE,JL))
            TXEJ(2,LE,J) = DXEJ(2,LE,J) -
     >                     (WTJ1*DXEJ(2,LE,1) + WTJ2*DXEJ(2,LE,JL))
            TXEJ(3,LE,J) = DXEJ(3,LE,J) -
     >                     (WTJ1*DXEJ(3,LE,1) + WTJ2*DXEJ(3,LE,JL))

         END DO
         LC = LC + 2
         LE = LE + 1
      END DO

      LC = 1
      LE = 1

      DO J = 1, JL, JL - 1
         DO K = 1, KL

            WTK2 = V0(J,K)
            WTK1 = ONE - WTK2

            DXCK(1,LE,K) = WTK1*DXC(1,LC) +
     >                     WTK2*DXC(1,LC+2) - DFACE(1,J,K)
            DXCK(2,LE,K) = WTK1*DXC(2,LC) +
     >                     WTK2*DXC(2,LC+2) - DFACE(2,J,K)
            DXCK(3,LE,K) = WTK1*DXC(3,LC) +
     >                     WTK2*DXC(3,LC+2) - DFACE(3,J,K)

            TXEK(1,LE,K) = DXEK(1,LE,K) -
     >                     (WTK1*DXEK(1,LE,1) + WTK2*DXEK(1,LE,KL))
            TXEK(2,LE,K) = DXEK(2,LE,K) -
     >                     (WTK1*DXEK(2,LE,1) + WTK2*DXEK(2,LE,KL))
            TXEK(3,LE,K) = DXEK(3,LE,K) -
     >                     (WTK1*DXEK(3,LE,1) + WTK2*DXEK(3,LE,KL))

         END DO
         LC = LC + 1
         LE = LE + 1
      END DO

C     Fix up the borders:

      DO J = 1, JL
         DFACE(1,J,1)  = DXCJ(1,1,J) + TXEJ(1,1,J) + DFACE(1,J,1)
         DFACE(2,J,1)  = DXCJ(2,1,J) + TXEJ(2,1,J) + DFACE(2,J,1)
         DFACE(3,J,1)  = DXCJ(3,1,J) + TXEJ(3,1,J) + DFACE(3,J,1)
         DFACE(1,J,KL) = DXCJ(1,2,J) + TXEJ(1,2,J) + DFACE(1,J,KL)
         DFACE(2,J,KL) = DXCJ(2,2,J) + TXEJ(2,2,J) + DFACE(2,J,KL)
         DFACE(3,J,KL) = DXCJ(3,2,J) + TXEJ(3,2,J) + DFACE(3,J,KL)
      END DO

      DO K = 2, KL - 1
         DFACE(1,1,K)  = DXCK(1,1,K) + TXEK(1,1,K) + DFACE(1,1,K)
         DFACE(2,1,K)  = DXCK(2,1,K) + TXEK(2,1,K) + DFACE(2,1,K)
         DFACE(3,1,K)  = DXCK(3,1,K) + TXEK(3,1,K) + DFACE(3,1,K)
         DFACE(1,JL,K) = DXCK(1,2,K) + TXEK(1,2,K) + DFACE(1,JL,K)
         DFACE(2,JL,K) = DXCK(2,2,K) + TXEK(2,2,K) + DFACE(2,JL,K)
         DFACE(3,JL,K) = DXCK(3,2,K) + TXEK(3,2,K) + DFACE(3,JL,K)
      END DO

C     Perturb block faces implicitly affected by explicit edge changes:

      DO K = 2, KL - 1
      DO J = 2, JL - 1
         WTJ2 = U0(J,K)
         WTJ1 = ONE - WTJ2
         WTK2 = V0(J,K)
         WTK1 = ONE - WTK2

         DELJ = WTJ1 * DXCK(1,1,K) + WTJ2 * DXCK(1,2,K)
         DELK = WTK1 * DXCJ(1,1,J) + WTK2 * DXCJ(1,2,J)

         DFACE(1,J,K) = (ABS (DELJ) * DELJ + ABS (DELK) * DELK) /
     >                   MAX (ABS (DELJ) + ABS (DELK), EPS)
     >                 + DFACE(1,J,K)

         DELJ = WTJ1 * DXCK(2,1,K) + WTJ2 * DXCK(2,2,K)
         DELK = WTK1 * DXCJ(2,1,J) + WTK2 * DXCJ(2,2,J)

         DFACE(2,J,K) = (ABS (DELJ) * DELJ + ABS (DELK) * DELK) /
     >                   MAX (ABS (DELJ) + ABS (DELK), EPS)
     >                 + DFACE(2,J,K)

         DELJ = WTJ1 * DXCK(3,1,K) + WTJ2 * DXCK(3,2,K)
         DELK = WTK1 * DXCJ(3,1,J) + WTK2 * DXCJ(3,2,J)

         DFACE(3,J,K) = (ABS (DELJ) * DELJ + ABS (DELK) * DELK) /
     >                   MAX (ABS (DELJ) + ABS (DELK), EPS)
     >                 + DFACE(3,J,K)

         DELJ = WTJ1 * TXEK(1,1,K) + WTJ2 * TXEK(1,2,K)
         DELK = WTK1 * TXEJ(1,1,J) + WTK2 * TXEJ(1,2,J)

         DFACE(1,J,K) = DFACE(1,J,K) + DELJ + DELK

         DELJ = WTJ1 * TXEK(2,1,K) + WTJ2 * TXEK(2,2,K)
         DELK = WTK1 * TXEJ(2,1,J) + WTK2 * TXEJ(2,2,J)

         DFACE(2,J,K) = DFACE(2,J,K) + DELJ + DELK

         DELJ = WTJ1 * TXEK(3,1,K) + WTJ2 * TXEK(3,2,K)
         DELK = WTK1 * TXEJ(3,1,J) + WTK2 * TXEJ(3,2,J)

         DFACE(3,J,K) = DFACE(3,J,K) + DELJ + DELK

      END DO
      END DO

      END SUBROUTINE WARPFACE

C+------------------------------------------------------------------------------
C
      SUBROUTINE BILINT (IDIM, JDIM, X, Y, XTARGET, YTARGET,
     >                   ITARGET, JTARGET, EPS, P, Q, STATUS)
C
C ACRONYM: BILinear INTerpolation (quasi-rectangular data in 2 dimensions)
C          ---      ---
C
C DESCRIPTION:
C
C         Interpolation within a quasi-rectangular (sub)mesh involves:
C     (1) locating the cell (if any) containing the current target point, and
C     (2) deriving a function value from the values of the function at that
C     cell's vertices (and there may be more than one function).
C
C         This routine serves both requirements: if the cell indicated by (i,j)
C     contains the target point (x,y), an interpolated function value is given
C     by a generalization of the 4-point bilinear formula for rectangular data:
C
C         F(x,y) = (1-p)(1-q)Fi,j + p(1-q)Fi+1,j + (1-p)qFi,j+1 + pqFi+1,j+1
C
C     for some p and q in [0, 1].  If not, a direction to move in is indicated.
C
C         Two nonlinear equations in p and q are obtained by applying the above
C     to X and Y:
C
C     f1 = (1-p)(1-q)Xi,j + p(1-q)Xi+1,j + (1-p)qXi,j+1 + pqXi+1,j+1 - Xtarg = 0
C     f2 = (1-p)(1-q)Yi,j + p(1-q)Yi+1,j + (1-p)qYi,j+1 + pqYi+1,j+1 - Ytarg = 0
C
C         Solution by a Newton iteration generalizes nicely to the trilinear
C     case.  Starting guesses of 0.5 should lead to rapid convergence if the
C     cell is the correct one.  Iterates outside [-4, +4] terminate the
C     calculation immediately to save time and guard against overflow: as with
C     convergence to a solution outside [0, 1], a direction for moving i and j
C     is implied.
C
C         Singularity in the Jacobian matrix is probably an indication of
C     singularity in the mesh.  Proceeding to the next cell in the higher
C     level's search strategy is probably appropriate.  Failure with all
C     cells means either the target point is outside the mesh or the mesh
C     is corrupted.  Handling this is application-dependent.
C
C         For a comprehensive treatment of the 3D case, see the INT3D program
C     in Pieter Buning's PLOT3D Tools collection at NASA Ames.  It supports
C     PLOT3D's IBLANK feature.  This routine is an attempt to provide a more
C     general-purpose lower level module without too much loss in efficiency,
C     and with the 3D analogue in mind.
C
C         Since the functions associated with the mesh may be stored either as
C     F(i,j,1:n) PLOT3D-style or as F(1:n,i,j) and similar uncertainty applies
C     to storing the interpolated function values, the evaluations of F must
C     be left to the higher level.
C
C ARGUMENTS:
C
C     ARG    DIM   TYPE I/O/S DESCRIPTION
C     IDIM,         I     I   Dimensions of X and Y in the calling program
C     JDIM                    (all cells not necessarily active)
C     X,  IDIM,JDIM R     I   Coordinates of the quasi-rectangular (sub)mesh;
C     Y                       those for the target cell are assumed meaningful
C     XTARGET,      R     I   Coordinates of the target point
C     YTARGET
C     ITARGET,      I     I   These indicate the cell to be processed on this
C     JTARGET                 call.  See STATUS.
C     EPS           R     I   Tolerance used to determine whether the computed
C                             P & Q are effectively within [0, 1], and also to
C                             help terminate the Newton iteration.  Suggestion:
C                             1.E-6 or 1.E-12 for 32- or 64-bit arithmetic.
C     P,            R     O   Interpolation coefficients as indicated above.
C     Q                       See STATUS.
C     STATUS        I     O   0:  ITARGET & JTARGET define the enclosing cell
C                                 and P & Q are inside [0-EPS, 1+EPS] and are
C                                 usable as bilinear interplation coefficients.
C                             1:  P & Q are outside [0, 1] and their values
C                                 suggest how to adjust ITARGET & JTARGET on
C                                 the next call.  E.g.:
C                                 IF (P .GT. 1.) ITARGET = MIN (ITARGET + 1, I2)
C                                 Avoiding an infinite loop may be tricky.
C                            -1:  No P & Q could be calculated because of
C                                 matrix singularity.  See discussion above.
C
C PROCEDURES:
C
C     LUSOLVE    Square system solution by LU decomposition with pivoting
C                (now in-line for efficiency)
C
C ERROR HANDLING:
C
C     None, for efficiency.  ITARGET & JTARGET are assumed to be valid.
C
C HISTORY:
C
C   07/16/93   DAS    Initial implementation, with 3D version in mind,
C                     following discussion with James Reuther and Scott
C                     Thomas, and perusal of Pieter Buning's INT3D.
C   11/24/93    "     NI, NJ arguments replaced by JDIM argument.
C                     All are actually redundant, but submesh usage is clearer.
C   11/01/96    "     In-lined LUSOLVE.
C   12/08/97    "     FNORM wasn't being tested upon hitting MXITER, yet
C                     p and q could both be in [0, 1].  Non-convex cells can
C                     still meet both convergence tests with non-unique
C                     solutions - no work-around is provided here. MXITER
C                     has been raised from 4 (really 5) to 8 in view of
C                     the missing step-halving and possibly tiny EPS.
C   01/08/98    "     FNORM was unintentionally being updated as ||correction||.
C   01/10/98    "     Application to unnormalized data requires a reference
C                     length, but we can't change the arguments at this stage,
C                     so a combination of current cell size and target (X,Y)
C                     units is used to determine a residual scale factor.
C
C AUTHOR: David Saunders, Sterling Software/NASA Ames, Moffett Field, CA.
C
C-------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER
     >   IDIM, JDIM, ITARGET, JTARGET, STATUS
      REAL
     >   X (IDIM, JDIM), Y (IDIM, JDIM), XTARGET, YTARGET, EPS, P, Q

C     Local constants:

      INTEGER
     >   MAXITER
      REAL
     >   HALF, TOOBIG, TOOLOW, ZERO
      PARAMETER
     >  (MAXITER = 8, HALF = 0.5E+0, TOOBIG = 4.E+0, TOOLOW = -3.E+0,
     >   ZERO = 0.E+0)

C     Local variables:

      INTEGER
     >   I, J, K, M
      REAL
     >   AJAC (2, 2), DP (2), FNORM, PK, QK, RTOL, T,
     >   XP, YP, XQ, YQ, XPQ, YPQ, XT, YT

C     Procedures:

C**** EXTERNAL
C****>   LUSOLVE


C     Execution:

      STATUS = 1      ! Assume most calls are searching for the right cell
      I = ITARGET
      J = JTARGET
      XT = X (I, J) - XTARGET
      YT = Y (I, J) - YTARGET
      XQ = X (I, J + 1) - X (I, J)
      YQ = Y (I, J + 1) - Y (I, J)
      XP = X (I + 1, J) - X (I, J)
      YP = Y (I + 1, J) - Y (I, J)
      XPQ = X (I + 1, J + 1) - X (I, J + 1) - XP
      YPQ = Y (I + 1, J + 1) - Y (I, J + 1) - YP
      RTOL = EPS * MAX (ABS(XP), ABS(YP), ABS (XTARGET), ABS (YTARGET))

C     For any usable solution to P & Q, the Newton iteration should not
C     need safeguarding - hence no step-halving inner iteration is included
C     to ensure that ||f|| decreases the way one would in general.

      K = 0
      PK = HALF
      QK = HALF

  300 CONTINUE

         AJAC (1, 1) = XP + QK * XPQ      ! J (i, j) = partial df(i)/dp(j)
         AJAC (2, 1) = YP + QK * YPQ
         AJAC (1, 2) = XQ + PK * XPQ
         AJAC (2, 2) = YQ + PK * YPQ

         DP (1) = PK * AJAC (1, 1) + QK * XQ + XT     ! f(1) and f (2)
         DP (2) = PK * AJAC (2, 1) + QK * YQ + YT

         FNORM = MAX (ABS (DP (1)), ABS (DP (2)))


C*****   CALL LUSOLVE (2, 2, AJAC, DP, IER)       ! Solve J dp = f

C        Cray inlining doesn't seem to help, so do the 2x2 case explicitly
C        to avoid an external reference and to avoid degenerate loops:

C        Perform the LU factorization, solving L y = b as we go:

C        Determine pivot element:

         M = 1
         T = AJAC (1, 1)
         IF (ABS (AJAC (2, 1)) .GT. ABS (T)) THEN
            AJAC (1, 1) = AJAC (2, 1)
            AJAC (2, 1) = T
            M = 2
         END IF

         IF (AJAC (1, 1) .EQ. ZERO) THEN  ! Singular matrix
            STATUS = -1
            GO TO 999
         END IF

C        Store the multiplier:

         AJAC (2, 1) = -AJAC (2, 1) / AJAC (1, 1)

C        Apply the multiplier to current submatrix and RHS:

         T = AJAC (M, 2)
         AJAC (M, 2) = AJAC (1, 2)
         AJAC (1, 2) = T
         AJAC (2, 2) = AJAC (2, 2) + AJAC (2, 1) * T

         T = DP (M)
         DP (M) = DP (1)
         DP (1) = T
         DP (2) = DP (2) + AJAC (2, 1) * T

         IF (AJAC (2, 2) .EQ. ZERO) THEN
            STATUS = -1       ! Singular matrix
            GO TO 999
         END IF

C        Back substitution (solution of U x = y):

         DP (2) = DP (2) / AJAC (2, 2)
         DP (1) = (DP (1) - AJAC (1, 2) * DP (2)) / AJAC (1, 1)


         PK = PK - DP (1)     ! New solution
         QK = QK - DP (2)

C        Assume that more often than not, we're searching for the right I, J
C        rather than being in the eventual enclosing cell. Exit early if we can.

         IF (MAX (PK, QK) .GT. TOOBIG) GO TO 900
         IF (MIN (PK, QK) .LT. TOOLOW) GO TO 900

         IF (FNORM .GT. RTOL) THEN
            K = K + 1
            IF (K .LT. MAXITER) GO TO 300
         ELSE
C           Drop through - either converged for the correct cell, converged for
C           the wrong cell (a nearby one?), or unconverged.  As for a diverged
C           case, the current solution in the latter two cases should still
C           point to a better I and J for the next call.

CCC      K = K + 1

            IF (MAX (ABS (PK - HALF), ABS (QK - HALF)) .LT. HALF + EPS)
     >         STATUS = 0


         END IF

  900 P = PK
      Q = QK


  999 RETURN
      END

C+------------------------------------------------------------------------------
C
      SUBROUTINE PARAM2D (IDIM, JDIM, I1, I2, J1, J2, X, Y, Z, U, V)
C
C ONE-LINER: PARAMeterization in 2-space of an XYZ surface
C            -----               -
C PURPOSE:
C
C        PARAM2D parameterizes the given (sub)mesh on a surface using
C     the usual chord-length approximation to arc-length between grid
C     points.  The arc length of each row and column is generally
C     normalized to 1, but this version can suppress normalization.
C     U(I1,J1) should be input as 0. or -999. respectively.
C     This version also handles degenerate lines by inserting uniform
C     u or v for the normalized case.
C
C        It is hoped that tying the dimensions of U and V to those of
C     X and Y (to keep the argument list modest) is not unwise, though
C     it could be if small submeshes are commonly involved.
C
C HISTORY:
C
C     11/19/93  DAS  Initial implementation.
C     12/03/93   "   Normalized each row and column.
C     03/20/95   "   Suppress normalization kludge: U(I1,J1) = -999.
C     02/06/98   "   Handled degenerate lines by inserting uniform u or v.
C
C AUTHOR: David Saunders, Sterling Software/NASA Ames, Mt. View, CA.
C
C ------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER IDIM, JDIM       ! (I) Grid dimensions in the calling program.

      INTEGER I1, I2, J1, J2   ! (I) Define the submesh to be parameterized.

      REAL    X (IDIM, JDIM),  ! (I) The surface grid coordinates.
     >        Y (IDIM, JDIM),
     >        Z (IDIM, JDIM)

      REAL    U (IDIM, JDIM),  ! (O) The chord-length-based parameterization:
     >        V (IDIM, JDIM)   !     U (I, J) = SUM (K=2:I) dS (K) / S (J) where
                               !     dS (K) ** 2 = (X (K, J) - X (K-1, J) ** 2 +
                               !                   (Y (K, J) - Y (K-1, J) ** 2 +
                               !                   (Z (K, J) - Z (K-1, J) ** 2
                               !     and S (J) is the total length of row J;
                               !     similarly for V (I, J).
                               ! (I) KLUDGE: If U(I1,J1) = -999. on input,
                               !     the normalization is suppressed.
C-------------------------------------------------------------------------------

C     Local constants:

      REAL, PARAMETER :: EPS = 1.E-6, FLAG = -999., ONE = 1., ZERO = 0.

C     Local variables:

      INTEGER   I, J
      REAL      RLENGTH
      LOGICAL   NORM

C     Execution:

      NORM = U (I1, J1) /= FLAG

      DO J = J1, J2

         U (I1, J) = ZERO

         DO I = I1 + 1, I2
            U (I, J) = U (I - 1, J)  +  SQRT (
     >         (X (I, J) - X (I - 1, J)) ** 2 +
     >         (Y (I, J) - Y (I - 1, J)) ** 2 +
     >         (Z (I, J) - Z (I - 1, J)) ** 2)
         END DO

         IF (NORM) THEN

            IF (U (I2, J) > EPS) THEN
               RLENGTH = ONE / U (I2, J)
               DO I = I1 + 1, I2 - 1
                  U (I, J) = U (I, J) * RLENGTH
               END DO
            ELSE
               RLENGTH = ONE / REAL (I2 - I1)
               DO I = I1 + 1, I2 - 1
                  U (I, J) = REAL (I - I1) * RLENGTH
               END DO
            END IF
            U (I2, J) = ONE

         END IF

      END DO


      DO I = I1, I2

         V (I, J1) = ZERO

         DO J = J1 + 1, J2
            V (I, J) = V (I, J - 1)  +  SQRT (
     >         (X (I, J) - X (I, J - 1)) ** 2 +
     >         (Y (I, J) - Y (I, J - 1)) ** 2 +
     >         (Z (I, J) - Z (I, J - 1)) ** 2)
         END DO

         IF (NORM) THEN

            IF (V (I, J2) .GT. EPS) THEN
               RLENGTH = ONE / V (I, J2)
               DO J = J1 + 1, J2 - 1
                  V (I, J) = V (I, J) * RLENGTH
               END DO
            ELSE
               RLENGTH = ONE / REAL (J2 - J1)
               DO J = J1 + 1, J2 - 1
                  V (I, J) = REAL (J - J1) * RLENGTH
               END DO
            END IF
            V (I, J2) = ONE

         END IF

      END DO

      END SUBROUTINE PARAM2D

C+------------------------------------------------------------------------------
C
      SUBROUTINE RIPPLE2D (MDIM, NDIM, M1, M2, N1, N2, X, Y, XP, YP,
     >                     M, N, EPS, P, Q, STATUS)
C
C     ONE-LINER:  "RIPPLE" search of a 2D quasi-rectangular grid
C                  ------              --
C     PURPOSE:
C
C        RIPPLE2D performs a search within a quasi-rectangular 2D grid for
C     the cell containing a target point by using a "ripple" pattern.  The
C     search should be efficient for typical applications where the specified
C     initial cell is a good starting guess (probably the result of a previous
C     call).  It also returns the corresponding bilinear interpolation
C     coefficients (though they may not apply to, say, bicubic applications).
C
C     METHOD:
C
C        Points are tested around the starting guess in rectangular "layers"
C     of increasing size.  If the cell defined by (M, N) does not contain
C     the target point, then the next search domain would be the surrounding
C     grid layer defined by the corner cells (M + 1, N + 1), (M + 1, N - 1),
C     (M - 1, N - 1), and (M - 1, N + 1).  If the cell is not found then the
C     next domain would be the layer defined by (M + 2, N + 2), (M + 2, N - 2),
C     (M - 2, N - 2), and (M - 2, N + 2), and so on.  The indices for each of
C     the corners are incremented or decremented until they reach the indicated
C     grid limits (M1 <= I <= M2 - 1 and N1 <= J <= N2 - 1) after which they
C     are held constant.  The search stops when a match is found or after all
C     corners of the layer have reached the corners of the (sub)grid, by which
C     time all specified cells will have been tested.
C
C        Within the rectangular layer, the corner points and the points interior
C     to them, here called edges, are treated separately for the sake of
C     symmetry.  An edge is no longer tested after the boundary of a grid has
C     been encountered by it.  Both corners and edges are cycled through in
C     clockwise directions as follows (where L is the current layer or level):
C
C
C          (M - L, N + L)                           (M + L, N + L)
C
C                          4 +. . . . . . . . .+ 1
C                            .        4        .  \
C                            .                 .   \
C                            .                 .    Corner
C                            .                 .
C                            . 3      +      1<------Edge
C                            .         (M,N)   .
C                            .                 .
C                            .                 .
C                            .        2        .
C                          3 +. . . . . . . . .+ 2
C
C          (M - L, N - L)                           (M + L, N - L)
C
C
C        The subroutine BILINT is used to perform the test.  BILINT returns
C     bilinear interpolation coefficients which are not used here but are
C     returned to the higher level.
C
C     ARGUMENTS:
C
C     NAME    DIM    TYPE  I/O/S     DESCRIPTION
C     ----    ---    ----  -----     -----------
C     MDIM,    -       I   I         Dimensions of the grid arrays in the
C     NDIM                           calling program
C     M1,      -       I   I         Indices specifying the (sub)grid to be
C     M2,                            searched as points (M1:M2, N1:N2)
C     N1,
C     N2
C     X,   (MDIM,NDIM) R   I         Grid coordinates
C     Y
C     XP,      -       R   I         Coordinates of the target point
C     YP
C     M,       -       I   I O       Input indices of starting guess and
C     N                              output indices of "lower left corner"
C                                    of the enclosing cell (or nearest to it)
C     EPS      -       R   I         Search tolerance used by BILINT, q.v.
C     P,       -       R   O         Bilinear interpolation coefficients;
C     Q                              see BILINT for their usage
C     STATUS   -       I   O         0 = match found;
C                                    1 = failed to find match; M, N, P, Q
C                                        represent the smallest mismatch
C                                    (Note that for the sake of efficiency, the
C                                    matrix singularity condition checked for
C                                    by BILINT is not tested for here.)
C
C     PROCEDURES:
C
C        BILINT     Performs test for enclosing cell and returns P, Q
C
C     ENVIRONMENT:
C
C        FORTRAN 77 with minor extensions
C
C     HISTORY:
C
C        11/10/93  M.Wong      Original design and coding.
C        11/24/93  D.Saunders  Arguments changed to allow search of sub-
C                              grids, etc.
C        11/28/93      "       Ensured "lower left corner" indices (not
C                              "upper" or "right" corner corresponding to
C                              P or Q = 1.), as for 1-D utility INTERVAL.
C        08/15/94      "       Egregious omission rectified: P, Q are now
C                              returned for bilinear applications (not re-
C                              quired for the original bicubic applications).
C                              Also: success on the initial BILINT call was
C                              skipping the possible adjustment of 11/28/93.
C        01/23/96      "       Returned the nearest M, N, P, Q to treat possible
C                              difficulties near the edges of (sub)regions.
C                              This can help establish an appropriate EPS too.
C        01/10/98      "       The 11/28/93 adjustment is dubious - should at
C                              least update P or Q via BILINT.
C     AUTHOR:
C
C        Michael D. Wong, Sterling Software, NASA/Ames, Moffett Field, CA
C
C-------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments.

      INTEGER  MDIM, NDIM, M1, M2, N1, N2, M, N, STATUS
      REAL     X (MDIM, NDIM), Y (MDIM, NDIM), XP, YP, EPS, P, Q

C     Local constants.

      REAL       ONE, ZERO
      PARAMETER (ONE = 1.E+0, ZERO = 0.E+0)

C     Local variables.

      INTEGER    I, J, L, ICORNER, M1P1, M2M1, M2M2, N1P1, N2M1, N2M2,
     >           NEDGE, ISTART, IEND, JEND, JSTART, MBEST, NBEST
      REAL       EBEST, ERR, PBEST, PERR, QBEST, QERR
      LOGICAL    ITEST, JTEST, TEST

C     Execution.


C     Begin with the initial guess.

      CALL BILINT (MDIM, NDIM, X, Y, XP, YP, M, N, EPS, P, Q, STATUS)

      IF (STATUS .EQ. 0) GO TO 90    ! M, N, P, Q may need adjusting

      IF (P .LT. ZERO) THEN
         PERR = -P
      ELSE IF (P .GT. ONE) THEN
         PERR = P - ONE
      ELSE
         PERR = ZERO
      END IF

      IF (Q .LT. ZERO) THEN
         QERR = -Q
      ELSE IF (Q .GT. ONE) THEN
         QERR = Q - ONE
      ELSE
         QERR = ZERO
      END IF

      EBEST = PERR + QERR
      PBEST = P
      QBEST = Q
      MBEST = M
      NBEST = N

C     Enter the rectangular layer search pattern.

      M1P1 = M1 + 1
      M2M1 = M2 - 1
      M2M2 = M2M1 - 1
      N1P1 = N1 + 1
      N2M1 = N2 - 1
      N2M2 = N2M1 - 1

      L = 0
   10 CONTINUE

         L = L + 1  ! Increment "ripple" level.

C        Define and test the corner points.

         DO ICORNER = 1, 4

            ITEST = .TRUE.
            JTEST = .TRUE.

            IF (ICORNER .EQ. 1) THEN

               I = M + L
               J = N + L

               IF (I .GT. M2M1) THEN
                  I     = M2M1
                  ITEST = .FALSE.
               END IF

               IF (J .GT. N2M1) THEN
                  J     = N2M1
                  JTEST = .FALSE.
               END IF

            ELSE IF (ICORNER .EQ. 2) THEN

               I = M + L
               J = N - L

               IF (I .GT. M2M1) THEN
                  I     = M2M1
                  ITEST = .FALSE.
               END IF

               IF (J .LT. N1) THEN
                  J     = N1
                  JTEST = .FALSE.
               END IF

            ELSE IF (ICORNER .EQ. 3) THEN

               I = M - L
               J = N - L

               IF (I .LT. M1) THEN
                  I     = M1
                  ITEST = .FALSE.
               END IF

               IF (J .LT. N1) THEN
                  J     = N1
                  JTEST = .FALSE.
               END IF

            ELSE  ! ICORNER = 4

               I = M - L
               J = N + L

               IF (I .LT. M1) THEN
                  I     = M1
                  ITEST = .FALSE.
               END IF

               IF (J .GT. N2M1) THEN
                  J     = N2M1
                  JTEST = .FALSE.
               END IF

            END IF

C           Avoid repeated tests.

            IF (ITEST .OR. JTEST) THEN

               CALL BILINT (MDIM, NDIM, X, Y, XP, YP, I, J, EPS,
     >                      P, Q, STATUS)

               IF (STATUS .EQ. 0) THEN
                  M = I
                  N = J
                  GO TO 90
               END IF

               IF (P .LT. ZERO) THEN
                  PERR = -P
               ELSE IF (P .GT. ONE) THEN
                  PERR = P - ONE
               ELSE
                  PERR = ZERO
               END IF

               IF (Q .LT. ZERO) THEN
                  QERR = -Q
               ELSE IF (Q .GT. ONE) THEN
                  QERR = Q - ONE
               ELSE
                  QERR = ZERO
               END IF

               ERR = PERR + QERR
               IF (ERR .LT. EBEST) THEN
                  EBEST = ERR
                  PBEST = P
                  QBEST = Q
                  MBEST = I
                  NBEST = J
               END IF
            END IF
         END DO

C        Cycle through interior points of the connector lines one time only.

         DO NEDGE = 1, 4

            IF (NEDGE .EQ. 1) THEN

               ISTART = M + L
               IEND   = ISTART
               JSTART = N - L + 1
               JEND   = N + L - 1
               TEST   = IEND .LE. M2M1

            ELSE IF (NEDGE .EQ. 2) THEN

               ISTART = M - L + 1
               IEND   = M + L - 1
               JSTART = N - L
               JEND   = JSTART
               TEST   = JSTART .GE. N1

            ELSE IF (NEDGE .EQ. 3) THEN

               ISTART = M - L
               IEND   = ISTART
               JSTART = N - L + 1
               JEND   = N + L - 1
               TEST   = ISTART .GE. M1

            ELSE ! Edge 4.

               ISTART = M - L + 1
               IEND   = M + L - 1
               JSTART = N + L
               JEND   = JSTART
               TEST   = JEND .LE. N2M1

            END IF

            IF (TEST) THEN

C              Limit the search to the interior points of the edge.

               IF (ISTART .NE. IEND) THEN

                  IF (ISTART .LT. M1P1) ISTART = M1P1
                  IF (IEND   .GT. M2M2) IEND   = M2M2

               ELSE IF (JSTART .NE. JEND) THEN

                  IF (JSTART .LT. N1P1) JSTART = N1P1
                  IF (JEND   .GT. N2M2) JEND   = N2M2

               END IF

               DO I = ISTART, IEND
                  DO J = JSTART, JEND

                     CALL BILINT (MDIM, NDIM, X, Y, XP, YP, I, J,
     >                            EPS, P, Q, STATUS)

                     IF (STATUS .EQ. 0) THEN
                        M = I
                        N = J
                        GO TO 90
                     END IF

                     IF (P .LT. ZERO) THEN
                        PERR = -P
                     ELSE IF (P .GT. ONE) THEN
                        PERR = P - ONE
                     ELSE
                        PERR = ZERO
                     END IF

                     IF (Q .LT. ZERO) THEN
                        QERR = -Q
                     ELSE IF (Q .GT. ONE) THEN
                        QERR = Q - ONE
                     ELSE
                        QERR = ZERO
                     END IF

                     ERR = PERR + QERR
                     IF (ERR .LT. EBEST) THEN
                        EBEST = ERR
                        PBEST = P
                        QBEST = Q
                        MBEST = I
                        NBEST = J
                     END IF
                  END DO
               END DO
            END IF
         END DO

C        Check to see if the entire grid has been tested.

         IF (M - L .GT. M1 .OR. M + L .LT. M2M1 .OR.
     >       N - L .GT. N1 .OR. N + L .LT. N2M1)
     >GO TO 10  ! Next search level.


C     Dropped through - entire (sub)grid was tested without success.
C     Return the nearest result - may be adequate near a boundary:

      M = MBEST
      N = NBEST
      P = PBEST
      Q = QBEST


   90 CONTINUE

C     Ensure "lower left" as for INTERVAL's 1-D search:

      IF (P .GE. ONE) THEN
         IF (M + 1 .LT. M2) THEN
            M = M + 1
            CALL BILINT (MDIM, NDIM, X, Y, XP, YP, M, N, EPS, P, Q,
     >                   STATUS)
         END IF
      END IF

      IF (Q .GE. ONE) THEN
         IF (N + 1 .LT. N2) THEN
            N = N + 1
            CALL BILINT (MDIM, NDIM, X, Y, XP, YP, M, N, EPS, P, Q,
     >                   STATUS)
         END IF
      END IF

      RETURN
      END

C+------------------------------------------------------------------------------
C
      SUBROUTINE LOCASE (STRING)
C
C PURPOSE:  LOCASE changes all upper case letters in the given
C           character string to lower case.
C
C METHOD:   Each character in STRING is treated in turn.  The intrinsic
C           function INDEX effectively allows a table lookup, with
C           the local strings LOW and UPP acting as two tables.
C           This method avoids the use of CHAR and ICHAR, which appear
C           to behave differently on ASCII and EBCDIC machines.
C
C ARGUMENTS:
C    ARG    TYPE I/O/S DESCRIPTION
C  STRING    C   I/O   Character string possibly containing some
C                      uppercase letters on input;
C                      strictly lowercase letters on output with no
C                      change to any non-alphabetic characters.
C
C HISTORY:
C
C  09/10/1985  M. Saunders   Rewrite of version by Hooper/Kennelly, 1983.
C  July 1999   M. Rimlinger  Adaptation of existing UPCASE.
C
C AUTHOR (UPCASE): Michael Saunders, Systems Optimization Lab., Stanford U.
C
C-------------------------------------------------------------------------------

C   Arguments:

      CHARACTER, INTENT (INOUT) :: STRING * (*)

C   Local constants:

      CHARACTER, PARAMETER :: LOW*26 = 'abcdefghijklmnopqrstuvwxyz',
     .                        UPP*26 = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
C   Local variables:

      INTEGER    I, J
      CHARACTER  C*1

C   Execution:

      DO J = 1, LEN (STRING)
         C = STRING(J:J)
         IF (C >= 'A' .AND. C <= 'Z') THEN
            I = INDEX (UPP, C)
            IF (I > 0) STRING(J:J) = LOW(I:I)
         END IF
      END DO

      END SUBROUTINE LOCASE

C+----------------------------------------------------------------------
C
      SUBROUTINE SECOND (CPUSEC)
C
C PURPOSE:
C     Measures CPU time in same form as the CRAY utility.
C
C ARGUMENTS:
C    ARG      DIM  TYPE I/O/S DESCRIPTION
C    CPUSEC    -    R     O   CPU time used so far, in seconds.
C
C METHOD:
C     IRIS FORTRAN requires use of the C language's intrinsic clock
C     function, which returns CPU time used since the FIRST call to it.
C     However, because FORTRAN cannot access the intrinsic, we call a
C     C routine (iclock) that in turns calls clock.
C     Conversion to seconds is done here.
C
C KNOWN BUG:
C     The count of microseconds wraps around after about 36 minutes
C     for the original iclock form => the mclock form should be better.
C
C USAGE:
C        CALL SECOND (TIME1)
C        ::::::::::::::::::
C        CALL SECOND (TIME2)
C        TOTALT = TIME2 - TIME1
C
C ENVIRONMENT:  SGI IRIS, FORTRAN 77
C
C HISTORY:
C     08/31/82   Dan McKernan    VAX/VMS version.
C     05/30/90   David Saunders  IRIS version (iclock).
C     06/05/90   Dexter Hermstad IRIS version (continued).
C     03/21/97   D.Saunders      IRIS mclock form obtained from J.Reuther.
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      REAL CPUSEC

C     C routine to call a system utility:

C*****INTEGER iclock
      INTEGER mclock

C     Execution:

C*****CPUSEC = FLOAT (iclock()) * 1.E-6
      CPUSEC = FLOAT (mclock()) * 0.01

      END
