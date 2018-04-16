!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      SUBROUTINE WARPQ3D2 (IMIN, IMAX, JMIN, JMAX, KMIN, KMAX,                 &
                           I1, I2, J1, J2, K1, K2, X0, Y0, Z0, S0, X, Y, Z)
!
!  DESCRIPTION:
!
!        WARPQ3D2 is a Fortran 90 version of WARPQ3D, from which the following
!     description is retained.
!
!        WARPQ3D is the quasi-3D form of WARP2D.  It perturbs the
!     interior of one face of a 3-space grid block given perturbed
!     edges of that face, which is indicated by one pair of equal index
!     arguments.  (E.g.: I2 = I1 means a face in the J/K subspace.)
!
!        The two-stage algorithm uses an intermediate perturbation to
!     account for any corner motion, involving edges derived from the
!     original edges, then a second perturbation to account for any
!     differences between the intermediate edges and the specified
!     new edges.
!
!        The perturbed edges should be input as edges of the desired
!     output face.  The original relative arc-length increments in
!     each index direction should also be input.  See PARAMXYZ for
!     setting them up in preparation for multiple perturbations.
!
!  HISTORY:
!
!     12/16/94  DAS/JDR  Quasi-3D analogue of WARP2D's two-stage algorithm,
!                        with stage 1 of WARP3D's algorithm in mind.
!     12/23/94    DAS    Lots of code replaced by DELQ3D, which is the form
!                        of WARPQ3D needed by WARP3D.
!     02/08/96     "     DELQ3D does stage 1 only now (all that WARP3D needs).
!                        S0(*) values are renormalized before use, in case
!                        the routine is being applied to a subgrid.
!     12/12/09     "     Fortran 90 version with DFACE* work-space local now,
!                        and renamed from WARPQ3D to WARPQ3D2.
!
!  AUTHOR:  David Saunders/James Reuther, NASA Ames Research Center, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IMPLICIT NONE

!     Arguments:

      INTEGER, INTENT (IN) :: &
         IMIN, IMAX, JMIN, JMAX, KMIN, KMAX      ! Grid array dimensions.

      INTEGER, INTENT (IN) :: &
         I1, I2, J1, J2, K1, K2                  ! Define the active face,
                                                 ! with one pair equal.

      REAL, INTENT (IN), DIMENSION (IMIN:IMAX, JMIN:JMAX, KMIN:KMAX) :: &
         X0, Y0, Z0                              ! Original face coordinates

      REAL, INTENT (IN), DIMENSION (IMIN:IMAX, JMIN:JMAX, KMIN:KMAX, 3) :: & 
         S0                                      ! Relative arc-lengths for the
                                                 ! active face (and probably
                                                 ! the underlying 3D array)

      REAL, INTENT (INOUT), DIMENSION (IMIN:IMAX, JMIN:JMAX, KMIN:KMAX) :: &
         X, Y, Z                                 ! Grid coordinates: new edges
                                                 ! of a face in; full face out
!     Local constants:

      REAL, PARAMETER :: &
         ONE = 1.E+0

!     Local variables:

      INTEGER :: &
         I, J, K

      REAL :: &
         DELI, DELJ, DELK, WTI1, WTI2, WTJ1, WTJ2, WTK1, WTK2

      REAL, ALLOCATABLE, DIMENSION (:,:,:) :: &  ! For face perturbations:
         DFACEI, DFACEJ, DFACEK                  ! DFACEI (1:3, J1:J2, K:K2) =
                                                 ! dX,dY,dZ on an I face, etc.
!     Execution:
!     ----------

      IF (I1 == I2) THEN  ! I plane case:

         ALLOCATE (DFACEI (3, JMIN:JMAX, KMIN:KMAX), &  ! Avoid unallocated
                   DFACEJ (1, 1, 1), DFACEK (1, 1, 1))  ! arguments

!        Stage 1:
!        Handle any corner motion by generating an intermediate face with
!        the final corners but otherwise derived from the original edges.
!        Actually, just set up the appropriate face perturbations.

         CALL DELQ3D (IMIN, IMAX, JMIN, JMAX, KMIN, KMAX,     &
                      I1, I2, J1, J2, K1, K2, X0, Y0, Z0, S0, &
                      DFACEI, DFACEJ, DFACEK, X, Y, Z)

!        Stage 2:
!        Set up the perturbations from the intermediate edges to the final
!        edges, then interpolate them into the interior points.

         I = I1

!        J = J1 and J2 edge perturbations:

         DO K = K1 + 1, K2 - 1
            DFACEI (1,J1,K) = X (I,J1,K) - X0 (I,J1,K) - DFACEI (1,J1,K)
            DFACEI (2,J1,K) = Y (I,J1,K) - Y0 (I,J1,K) - DFACEI (2,J1,K)
            DFACEI (3,J1,K) = Z (I,J1,K) - Z0 (I,J1,K) - DFACEI (3,J1,K)
            DFACEI (1,J2,K) = X (I,J2,K) - X0 (I,J2,K) - DFACEI (1,J2,K)
            DFACEI (2,J2,K) = Y (I,J2,K) - Y0 (I,J2,K) - DFACEI (2,J2,K)
            DFACEI (3,J2,K) = Z (I,J2,K) - Z0 (I,J2,K) - DFACEI (3,J2,K)
         END DO

!        K = K1 and K2 edge perturbations:

         DO J = J1 + 1, J2 - 1
            DFACEI (1,J,K1) = X (I,J,K1) - X0 (I,J,K1) - DFACEI (1,J,K1)
            DFACEI (2,J,K1) = Y (I,J,K1) - Y0 (I,J,K1) - DFACEI (2,J,K1)
            DFACEI (3,J,K1) = Z (I,J,K1) - Z0 (I,J,K1) - DFACEI (3,J,K1)
            DFACEI (1,J,K2) = X (I,J,K2) - X0 (I,J,K2) - DFACEI (1,J,K2)
            DFACEI (2,J,K2) = Y (I,J,K2) - Y0 (I,J,K2) - DFACEI (2,J,K2)
            DFACEI (3,J,K2) = Z (I,J,K2) - Z0 (I,J,K2) - DFACEI (3,J,K2)
         END DO

!        Interior points: accumulate the (independent) contributions.

         DO K = K1 + 1, K2 - 1
            DO J = J1 + 1, J2 - 1
               WTJ2 = (S0 (I, J,  K, 2) - S0 (I, J1, K, 2)) / &
                      (S0 (I, J2, K, 2) - S0 (I, J1, K, 2))
               WTJ1 = ONE - WTJ2
               WTK2 = (S0 (I, J, K,  3) - S0 (I, J, K1, 3)) / &
                      (S0 (I, J, K2, 3) - S0 (I, J, K1, 3))
               WTK1 = ONE - WTK2

               DELJ = WTJ1 * DFACEI (1, J1, K) + WTJ2 * DFACEI (1, J2,K)
               DELK = WTK1 * DFACEI (1, J, K1) + WTK2 * DFACEI (1, J,K2)

               X (I, J, K) = X0 (I, J, K) + DFACEI (1,J,K) + DELJ + DELK

               DELJ = WTJ1 * DFACEI (2, J1, K) + WTJ2 * DFACEI (2, J2,K)
               DELK = WTK1 * DFACEI (2, J, K1) + WTK2 * DFACEI (2, J,K2)

               Y (I, J, K) = Y0 (I, J, K) + DFACEI (2,J,K) + DELJ + DELK

               DELJ = WTJ1 * DFACEI (3, J1, K) + WTJ2 * DFACEI (3, J2,K)
               DELK = WTK1 * DFACEI (3, J, K1) + WTK2 * DFACEI (3, J,K2)

               Z (I, J, K) = Z0 (I, J, K) + DFACEI (3,J,K) + DELJ + DELK
            END DO
         END DO

      ELSE IF (J1 == J2) THEN  ! J plane case:

         ALLOCATE (DFACEJ (3, IMIN:IMAX, KMIN:KMAX), &
                   DFACEI (1, 1, 1), DFACEK (1, 1, 1))

!        Stage 1:
!        Handle any corner motion by generating an intermediate face with
!        the final corners but otherwise derived from the original edges.
!        Actually, just set up the appropriate face perturbations.

         CALL DELQ3D (IMIN, IMAX, JMIN, JMAX, KMIN, KMAX,     &
                      I1, I2, J1, J2, K1, K2, X0, Y0, Z0, S0, &
                      DFACEI, DFACEJ, DFACEK, X, Y, Z)

!        Stage 2:
!        Set up the perturbations from the intermediate edges to the final
!        edges, then interpolate them into the interior points.

         J = J1

!        I = I1 and I2 edge perturbations:

         DO K = K1 + 1, K2 - 1
            DFACEJ (1,I1,K) = X (I1,J,K) - X0 (I1,J,K) - DFACEJ (1,I1,K)
            DFACEJ (2,I1,K) = Y (I1,J,K) - Y0 (I1,J,K) - DFACEJ (2,I1,K)
            DFACEJ (3,I1,K) = Z (I1,J,K) - Z0 (I1,J,K) - DFACEJ (3,I1,K)
            DFACEJ (1,I2,K) = X (I2,J,K) - X0 (I2,J,K) - DFACEJ (1,I2,K)
            DFACEJ (2,I2,K) = Y (I2,J,K) - Y0 (I2,J,K) - DFACEJ (2,I2,K)
            DFACEJ (3,I2,K) = Z (I2,J,K) - Z0 (I2,J,K) - DFACEJ (3,I2,K)
         END DO

!        K = K1 and K2 edge perturbations:

         DO I = I1 + 1, I2 - 1
            DFACEJ (1,I,K1) = X (I,J,K1) - X0 (I,J,K1) - DFACEJ (1,I,K1)
            DFACEJ (2,I,K1) = Y (I,J,K1) - Y0 (I,J,K1) - DFACEJ (2,I,K1)
            DFACEJ (3,I,K1) = Z (I,J,K1) - Z0 (I,J,K1) - DFACEJ (3,I,K1)
            DFACEJ (1,I,K2) = X (I,J,K2) - X0 (I,J,K2) - DFACEJ (1,I,K2)
            DFACEJ (2,I,K2) = Y (I,J,K2) - Y0 (I,J,K2) - DFACEJ (2,I,K2)
            DFACEJ (3,I,K2) = Z (I,J,K2) - Z0 (I,J,K2) - DFACEJ (3,I,K2)
         END DO

!        Interior points: accumulate the (independent) contributions.

         DO K = K1 + 1, K2 - 1
            DO I = I1 + 1, I2 - 1
               WTI2 = (S0 (I,  J, K, 1) - S0 (I1, J, K, 1)) / &
                      (S0 (I2, J, K, 1) - S0 (I1, J, K, 1))
               WTI1 = ONE - WTI2
               WTK2 = (S0 (I, J, K,  3) - S0 (I, J, K1, 3)) / &
                      (S0 (I, J, K2, 3) - S0 (I, J, K1, 3))
               WTK1 = ONE - WTK2

               DELI = WTI1 * DFACEJ (1, I1, K) + WTI2 * DFACEJ (1, I2,K)
               DELK = WTK1 * DFACEJ (1, I, K1) + WTK2 * DFACEJ (1, I,K2)

               X (I, J, K) = X0 (I, J, K) + DFACEJ (1,I,K) + DELI + DELK

               DELI = WTI1 * DFACEJ (2, I1, K) + WTI2 * DFACEJ (2, I2,K)
               DELK = WTK1 * DFACEJ (2, I, K1) + WTK2 * DFACEJ (2, I,K2)

               Y (I, J, K) = Y0 (I, J, K) + DFACEJ (2,I,K) + DELI + DELK

               DELI = WTI1 * DFACEJ (3, I1, K) + WTI2 * DFACEJ (3, I2,K)
               DELK = WTK1 * DFACEJ (3, I, K1) + WTK2 * DFACEJ (3, I,K2)

               Z (I, J, K) = Z0 (I, J, K) + DFACEJ (3,I,K) + DELI + DELK
            END DO
         END DO

      ELSE ! IF (K1 == K2) THEN  ! K plane case:

         ALLOCATE (DFACEK (3, IMIN:IMAX, JMIN:JMAX), &
                   DFACEI (1, 1, 1), DFACEJ (1, 1, 1))

!        Stage 1:
!        Handle any corner motion by generating an intermediate face with
!        the final corners but otherwise derived from the original edges.
!        Actually, just set up the appropriate face perturbations.

         CALL DELQ3D (IMIN, IMAX, JMIN, JMAX, KMIN, KMAX,     &
                      I1, I2, J1, J2, K1, K2, X0, Y0, Z0, S0, &
                      DFACEI, DFACEJ, DFACEK, X, Y, Z)

!        Stage 2:
!        Set up the perturbations from the intermediate edges to the final
!        edges, then interpolate them into the interior points.

         K = K1

!        I = I1 and I2 edge perturbations:

         DO J = J1 + 1, J2 - 1
            DFACEK (1,I1,J) = X (I1,J,K) - X0 (I1,J,K) - DFACEK (1,I1,J)
            DFACEK (2,I1,J) = Y (I1,J,K) - Y0 (I1,J,K) - DFACEK (2,I1,J)
            DFACEK (3,I1,J) = Z (I1,J,K) - Z0 (I1,J,K) - DFACEK (3,I1,J)
            DFACEK (1,I2,J) = X (I2,J,K) - X0 (I2,J,K) - DFACEK (1,I2,J)
            DFACEK (2,I2,J) = Y (I2,J,K) - Y0 (I2,J,K) - DFACEK (2,I2,J)
            DFACEK (3,I2,J) = Z (I2,J,K) - Z0 (I2,J,K) - DFACEK (3,I2,J)
         END DO

!        J = J1 and J2 edge perturbations:

         DO I = I1 + 1, I2 - 1
            DFACEK (1,I,J1) = X (I,J1,K) - X0 (I,J1,K) - DFACEK (1,I,J1)
            DFACEK (2,I,J1) = Y (I,J1,K) - Y0 (I,J1,K) - DFACEK (2,I,J1)
            DFACEK (3,I,J1) = Z (I,J1,K) - Z0 (I,J1,K) - DFACEK (3,I,J1)
            DFACEK (1,I,J2) = X (I,J2,K) - X0 (I,J2,K) - DFACEK (1,I,J2)
            DFACEK (2,I,J2) = Y (I,J2,K) - Y0 (I,J2,K) - DFACEK (2,I,J2)
            DFACEK (3,I,J2) = Z (I,J2,K) - Z0 (I,J2,K) - DFACEK (3,I,J2)
         END DO

!        Interior points: accumulate the (independent) contributions.

         DO J = J1 + 1, J2 - 1
            DO I = I1 + 1, I2 - 1
               WTI2 = (S0 (I,  J, K, 1) - S0 (I1, J, K, 1)) / &
                      (S0 (I2, J, K, 1) - S0 (I1, J, K, 1))
               WTI1 = ONE - WTI2
               WTJ2 = (S0 (I, J,  K, 2) - S0 (I, J1, K, 2)) / &
                      (S0 (I, J2, K, 2) - S0 (I, J1, K, 2))
               WTJ1 = ONE - WTJ2

               DELI = WTI1 * DFACEK (1, I1, J) + WTI2 * DFACEK (1, I2,J)
               DELJ = WTJ1 * DFACEK (1, I, J1) + WTJ2 * DFACEK (1, I,J2)

               X (I, J, K) = X0 (I, J, K) + DFACEK (1,I,J) + DELI + DELJ

               DELI = WTI1 * DFACEK (2, I1, J) + WTI2 * DFACEK (2, I2,J)
               DELJ = WTJ1 * DFACEK (2, I, J1) + WTJ2 * DFACEK (2, I,J2)

               Y (I, J, K) = Y0 (I, J, K) + DFACEK (2,I,J) + DELI + DELJ

               DELI = WTI1 * DFACEK (3, I1, J) + WTI2 * DFACEK (3, I2,J)
               DELJ = WTJ1 * DFACEK (3, I, J1) + WTJ2 * DFACEK (3, I,J2)

               Z (I, J, K) = Z0 (I, J, K) + DFACEK (3,I,J) + DELI + DELJ
            END DO
         END DO

      END IF

      DEALLOCATE (DFACEI, DFACEJ, DFACEK)

      END SUBROUTINE WARPQ3D2
