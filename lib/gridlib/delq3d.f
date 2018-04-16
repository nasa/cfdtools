C+------------------------------------------------------------------------------
C
      SUBROUTINE DELQ3D (IMIN, IMAX, JMIN, JMAX, KMIN, KMAX,
     >                   I1, I2, J1, J2, K1, K2, X0, Y0, Z0, S0,
     >                   DFACEI, DFACEJ, DFACEK, X, Y, Z)
C
C  ONE-LINER: 3-space surface perturbations (not XYZs) given new corners
C
C  DESCRIPTION:
C
C        DELQ3D performs stage 1 of the WARPQ3D 3-space surface grid
C     perturbation in a form which is reusable by WARP3D.  It returns
C     face perturbations rather than perturbed face coordinates, to
C     avoid work-space problems in WARP3D's first stage.  See WARP2D
C     and WARPQ3D for further details of these two-stage algorithms.
C
C        The three cases of a block face are handled here by three similar
C     code sections.  Eliminating WARP2D's special handling of the case of
C     fixed corners helps keep the bulk down.
C
C  HISTORY:
C
C     12/20/94  DAS/JJR  Adaptation of WARPQ3D for WARP3D.
C     12/23/94    DAS    DELQ3D is now used by WARPQ3D.
C     02/08/96     "     DELQ3D does only stage 1 now. Use on a subgrid with
C                        full-grid S0 requires transforming S0 here.
C
C  AUTHOR:  David Saunders/James Reuther, NASA Ames, Mt. View, CA.
C
C ------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments.

      INTEGER IMIN, IMAX, JMIN, JMAX, KMIN, KMAX      ! I Grid array dimensions.

      INTEGER I1, I2, J1, J2, K1, K2                  ! I Define active face,
                                                      !   one pair being equal.

      REAL    X0 (IMIN:IMAX, JMIN:JMAX, KMIN:KMAX),   ! I Original face coords.
     >        Y0 (IMIN:IMAX, JMIN:JMAX, KMIN:KMAX),   !   in appropriate places.
     >        Z0 (IMIN:IMAX, JMIN:JMAX, KMIN:KMAX)

      REAL    S0 (IMIN:IMAX, JMIN:JMAX, KMIN:KMAX, 3) ! I Relative arc-lengths
                                                      !   for the I,J,K lines -
                                                      !   see PARAMXYZ.  If from
                                                      !   a grid larger than the
                                                      !   active subgrid, S0 is
                                                      !   transformed here.

      REAL    DFACEI (3, JMIN:JMAX, KMIN:KMAX),    ! O Reqd. face perturbations:
     >        DFACEJ (3, IMIN:IMAX, KMIN:KMAX),    !   DFACEI (1:3,J1:J2,K1:K2)
     >        DFACEK (3, IMIN:IMAX, JMIN:JMAX)     !   = dX,dY,dZ on an I face
                                                   !   of a 3-space grid block.

      REAL    X (IMIN:IMAX, JMIN:JMAX, KMIN:KMAX), ! I Grid coordinates:
     >        Y (IMIN:IMAX, JMIN:JMAX, KMIN:KMAX), !   new edges of a face in;
     >        Z (IMIN:IMAX, JMIN:JMAX, KMIN:KMAX)  !   unchanged on output.

C-------------------------------------------------------------------------------

C     Local constants.

      REAL    EPS, ONE

      PARAMETER (EPS = 1.E-8, ONE = 1.E+0) ! EPS safeguards a divide by zero -
                                           ! presumably only if result is zero.
C     Local variables.

      INTEGER I, J, K
      REAL    DELI, DELJ, DELK, S0I1, S0J1, S0K1, SRANGEI, SRANGEJ,
     >        SRANGEK, WTI1, WTI2, WTJ1, WTJ2, WTK1, WTK2

C     Execution.
C     ----------


      IF (I1 .EQ. I2) THEN

C        I plane case:
C        -------------

         I = I1

C        Set up the corner perturbations:

         DO K = K1, K2, K2 - K1
            DO J = J1, J2, J2 - J1
               DFACEI (1, J, K) = X (I, J, K) - X0 (I, J, K)
               DFACEI (2, J, K) = Y (I, J, K) - Y0 (I, J, K)
               DFACEI (3, J, K) = Z (I, J, K) - Z0 (I, J, K)
            END DO
         END DO

C        Set up intermediate edge perturbations corresponding to the
C        final corners but otherwise derived from the original edges.

         DO J = J1, J2, J2 - J1
            S0K1 = S0 (I, J, K1, 3)
            SRANGEK = ONE / (S0 (I, J, K2, 3) - S0K1)
            DO K = K1 + 1, K2 - 1
               WTK2 = (S0 (I, J, K, 3) - S0K1) * SRANGEK
               WTK1 = ONE - WTK2
               DFACEI (1, J, K) = WTK1 * DFACEI (1, J, K1) +
     >                            WTK2 * DFACEI (1, J, K2)
               DFACEI (2, J, K) = WTK1 * DFACEI (2, J, K1) +
     >                            WTK2 * DFACEI (2, J, K2)
               DFACEI (3, J, K) = WTK1 * DFACEI (3, J, K1) +
     >                            WTK2 * DFACEI (3, J, K2)
            END DO
         END DO

         DO K = K1, K2, K2 - K1
            S0J1 = S0 (I, J1, K, 2)
            SRANGEJ = ONE / (S0 (I, J2, K, 2) - S0J1)
            DO J = J1 + 1, J2 - 1
               WTJ2 = (S0 (I, J, K, 2) - S0J1) * SRANGEJ
               WTJ1 = ONE - WTJ2
               DFACEI (1, J, K) = WTJ1 * DFACEI (1, J1, K) +
     >                            WTJ2 * DFACEI (1, J2, K)
               DFACEI (2, J, K) = WTJ1 * DFACEI (2, J1, K) +
     >                            WTJ2 * DFACEI (2, J2, K)
               DFACEI (3, J, K) = WTJ1 * DFACEI (3, J1, K) +
     >                            WTJ2 * DFACEI (3, J2, K)
            END DO
         END DO

C        Interpolate the intermediate perturbations of interior points.
C        The contributions from each pair of edges are not independent.

         DO K = K1 + 1, K2 - 1
            DO J = J1 + 1, J2 - 1
               WTJ2 = (S0 (I, J,  K, 2) - S0 (I, J1, K, 2)) /
     >                (S0 (I, J2, K, 2) - S0 (I, J1, K, 2))
               WTJ1 = ONE - WTJ2
               WTK2 = (S0 (I, J, K,  3) - S0 (I, J, K1, 3)) /
     >                (S0 (I, J, K2, 3) - S0 (I, J, K1, 3))
               WTK1 = ONE - WTK2

               DELJ = WTJ1 * DFACEI (1, J1, K) + WTJ2 * DFACEI (1, J2,K)
               DELK = WTK1 * DFACEI (1, J, K1) + WTK2 * DFACEI (1, J,K2)

               DFACEI (1, J, K) = (ABS (DELJ) *DELJ + ABS (DELK) *DELK)/
     >                             MAX (ABS (DELJ) + ABS (DELK), EPS)

               DELJ = WTJ1 * DFACEI (2, J1, K) + WTJ2 * DFACEI (2, J2,K)
               DELK = WTK1 * DFACEI (2, J, K1) + WTK2 * DFACEI (2, J,K2)

               DFACEI (2, J, K) = (ABS (DELJ) *DELJ + ABS (DELK) *DELK)/
     >                             MAX (ABS (DELJ) + ABS (DELK), EPS)

               DELJ = WTJ1 * DFACEI (3, J1, K) + WTJ2 * DFACEI (3, J2,K)
               DELK = WTK1 * DFACEI (3, J, K1) + WTK2 * DFACEI (3, J,K2)

               DFACEI (3, J, K) = (ABS (DELJ) *DELJ + ABS (DELK) *DELK)/
     >                             MAX (ABS (DELJ) + ABS (DELK), EPS)
            END DO
         END DO

      ELSE IF (J1 .EQ. J2) THEN

C        J plane case:
C        -------------

         J = J1

C        Corner perturbations:

         DO K = K1, K2, K2 - K1
            DO I = I1, I2, I2 - I1
               DFACEJ (1, I, K) = X (I, J, K) - X0 (I, J, K)
               DFACEJ (2, I, K) = Y (I, J, K) - Y0 (I, J, K)
               DFACEJ (3, I, K) = Z (I, J, K) - Z0 (I, J, K)
            END DO
         END DO

C        Intermediate edge perturbations:

         DO I = I1, I2, I2 - I1
            S0K1 = S0 (I, J, K1, 3)
            SRANGEK = ONE / (S0 (I, J, K2, 3) - S0K1)
            DO K = K1 + 1, K2 - 1
               WTK2 = (S0 (I, J, K, 3) - S0K1) * SRANGEK
               WTK1 = ONE - WTK2
               DFACEJ (1, I, K) = WTK1 * DFACEJ (1, I, K1) +
     >                            WTK2 * DFACEJ (1, I, K2)
               DFACEJ (2, I, K) = WTK1 * DFACEJ (2, I, K1) +
     >                            WTK2 * DFACEJ (2, I, K2)
               DFACEJ (3, I, K) = WTK1 * DFACEJ (3, I, K1) +
     >                            WTK2 * DFACEJ (3, I, K2)
            END DO
         END DO

         DO K = K1, K2, K2 - K1
            S0I1 = S0 (I1, J, K, 1)
            SRANGEI = ONE / (S0 (I2, J, K, 1) - S0I1)
            DO I = I1 + 1, I2 - 1
               WTI2 = (S0 (I, J, K, 1) - S0I1) * SRANGEI
               WTI1 = ONE - WTI2
               DFACEJ (1, I, K) = WTI1 * DFACEJ (1, I1, K) +
     >                            WTI2 * DFACEJ (1, I2, K)
               DFACEJ (2, I, K) = WTI1 * DFACEJ (2, I1, K) +
     >                            WTI2 * DFACEJ (2, I2, K)
               DFACEJ (3, I, K) = WTI1 * DFACEJ (3, I1, K) +
     >                            WTI2 * DFACEJ (3, I2, K)
            END DO
         END DO

C        Intermediate perturbations of interior points:

         DO K = K1 + 1, K2 - 1
            DO I = I1 + 1, I2 - 1
               WTI2 = (S0 (I,  J, K, 1) - S0 (I1, J, K, 1)) /
     >                (S0 (I2, J, K, 1) - S0 (I1, J, K, 1))
               WTI1 = ONE - WTI2
               WTK2 = (S0 (I, J, K,  3) - S0 (I, J, K1, 3)) /
     >                (S0 (I, J, K2, 3) - S0 (I, J, K1, 3))
               WTK1 = ONE - WTK2

               DELI = WTI1 * DFACEJ (1, I1, K) + WTI2 * DFACEJ (1, I2,K)
               DELK = WTK1 * DFACEJ (1, I, K1) + WTK2 * DFACEJ (1, I,K2)

               DFACEJ (1, I, K) = (ABS (DELI) *DELI + ABS (DELK) *DELK)/
     >                             MAX (ABS (DELI) + ABS (DELK), EPS)

               DELI = WTI1 * DFACEJ (2, I1, K) + WTI2 * DFACEJ (2, I2,K)
               DELK = WTK1 * DFACEJ (2, I, K1) + WTK2 * DFACEJ (2, I,K2)

               DFACEJ (2, I, K) = (ABS (DELI) *DELI + ABS (DELK) *DELK)/
     >                             MAX (ABS (DELI) + ABS (DELK), EPS)

               DELI = WTI1 * DFACEJ (3, I1, K) + WTI2 * DFACEJ (3, I2,K)
               DELK = WTK1 * DFACEJ (3, I, K1) + WTK2 * DFACEJ (3, I,K2)

               DFACEJ (3, I, K) = (ABS (DELI) *DELI + ABS (DELK) *DELK)/
     >                             MAX (ABS (DELI) + ABS (DELK), EPS)
            END DO
         END DO

      ELSE IF (K1 .EQ. K2) THEN

C        K plane case:
C        -------------

         K = K1

C        Corner perturbations:

         DO J = J1, J2, J2 - J1
            DO I = I1, I2, I2 - I1
               DFACEK (1, I, J) = X (I, J, K) - X0 (I, J, K)
               DFACEK (2, I, J) = Y (I, J, K) - Y0 (I, J, K)
               DFACEK (3, I, J) = Z (I, J, K) - Z0 (I, J, K)
            END DO
         END DO

C        Intermediate edge perturbations:

         DO I = I1, I2, I2 - I1
            S0J1 = S0 (I, J1, K, 2)
            SRANGEJ = ONE / (S0 (I, J2, K, 2) - S0J1)
            DO J = J1 + 1, J2 - 1
               WTJ2 = (S0 (I, J, K, 2) - S0J1) * SRANGEJ
               WTJ1 = ONE - WTJ2
               DFACEK (1, I, J) = WTJ1 * DFACEK (1, I, J1) +
     >                            WTJ2 * DFACEK (1, I, J2)
               DFACEK (2, I, J) = WTJ1 * DFACEK (2, I, J1) +
     >                            WTJ2 * DFACEK (2, I, J2)
               DFACEK (3, I, J) = WTJ1 * DFACEK (3, I, J1) +
     >                            WTJ2 * DFACEK (3, I, J2)
            END DO
         END DO

         DO J = J1, J2, J2 - J1
            S0I1 = S0 (I1, J, K, 1)
            SRANGEI = ONE / (S0 (I2, J, K, 1) - S0I1)
            DO I = I1 + 1, I2 - 1
               WTI2 = (S0 (I, J, K, 1) - S0I1) * SRANGEI
               WTI1 = ONE - WTI2
               DFACEK (1, I, J) = WTI1 * DFACEK (1, I1, J) +
     >                            WTI2 * DFACEK (1, I2, J)
               DFACEK (2, I, J) = WTI1 * DFACEK (2, I1, J) +
     >                            WTI2 * DFACEK (2, I2, J)
               DFACEK (3, I, J) = WTI1 * DFACEK (3, I1, J) +
     >                            WTI2 * DFACEK (3, I2, J)
            END DO
         END DO

C        Intermediate perturbations of interior points:

         DO J = J1 + 1, J2 - 1
            DO I = I1 + 1, I2 - 1
               WTI2 = (S0 (I,  J, K, 1) - S0 (I1, J, K, 1)) /
     >                (S0 (I2, J, K, 1) - S0 (I1, J, K, 1))
               WTI1 = ONE - WTI2
               WTJ2 = (S0 (I, J,  K, 2) - S0 (I, J1, K, 2)) /
     >                (S0 (I, J2, K, 2) - S0 (I, J1, K, 2))
               WTJ1 = ONE - WTJ2

               DELI = WTI1 * DFACEK (1, I1, J) + WTI2 * DFACEK (1, I2,J)
               DELJ = WTJ1 * DFACEK (1, I, J1) + WTJ2 * DFACEK (1, I,J2)

               DFACEK (1, I, J) = (ABS (DELI) *DELI + ABS (DELJ) *DELJ)/
     >                             MAX (ABS (DELI) + ABS (DELJ), EPS)

               DELI = WTI1 * DFACEK (2, I1, J) + WTI2 * DFACEK (2, I2,J)
               DELJ = WTJ1 * DFACEK (2, I, J1) + WTJ2 * DFACEK (2, I,J2)

               DFACEK (2, I, J) = (ABS (DELI) *DELI + ABS (DELJ) *DELJ)/
     >                             MAX (ABS (DELI) + ABS (DELJ), EPS)

               DELI = WTI1 * DFACEK (3, I1, J) + WTI2 * DFACEK (3, I2,J)
               DELJ = WTJ1 * DFACEK (3, I, J1) + WTJ2 * DFACEK (3, I,J2)

               DFACEK (3, I, J) = (ABS (DELI) *DELI + ABS (DELJ) *DELJ)/
     >                             MAX (ABS (DELI) + ABS (DELJ), EPS)
            END DO
         END DO

      END IF

      RETURN
      END
