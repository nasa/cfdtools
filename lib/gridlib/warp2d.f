C+------------------------------------------------------------------------------
C
      SUBROUTINE WARP2D (IMIN, IMAX, JMIN, JMAX, I1, I2, J1, J2,
     >                   X0, Y0, S0, DEDGEI, DEDGEJ, X, Y)
C
C  ONE-LINER: Perturb interior of a plane grid given new edges
C
C  DESCRIPTION:
C
C        WARP2D perturbs the interior points of a 2-space grid block
C     given the original grid and perturbed edges. All four edges are
C     assumed to be perturbed (though some may of course be fixed).
C     If all corners are found to be unperturbed, considerably less
C     work is performed.  In general, an intermediate perturbation is
C     required to account for the corner motion, using edges derived
C     from the original edges, then a second perturbation accounts
C     for differences between the intermediate edges and the specified
C     new edges.
C
C        The perturbed edges should be input as edges of the desired
C     output grid.  The original relative arc-length increments in
C     each index direction should also be input.  See PARAMXY for
C     setting them up in preparation for multiple perturbations.
C
C  ENVIRONMENT:  FORTRAN 77 with minor extensions
C
C  HISTORY:
C
C     03/04/94  DAS/JJR  WARP2D written to perturb a single 2-space line,
C                        for James's airfoil design code.
C     12/13/94   "   "   WARP2D adapted from generalized WARP3D to test
C                        variations of the algorithm.
C     12/14/94   "   "   Two-stage algorithm devised to handle corner motion.
C     02/01/96    DAS    Use on a subgrid with full-grid S0 must transform S0.
C
C  AUTHOR:  David Saunders/James Reuther, NASA Ames, Mt. View, CA.
C
C ------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments.

      INTEGER IMIN, IMAX, JMIN, JMAX       !  I  Grid array dimensions.

      INTEGER I1, I2, J1, J2               !  I  Active area is (I1:I2, J1:J2).

      REAL    X0 (IMIN:IMAX, JMIN:JMAX),   !  I  Original grid coordinates.
     >        Y0 (IMIN:IMAX, JMIN:JMAX)

      REAL    S0 (IMIN:IMAX, JMIN:JMAX, 2) !  I  Relative arc-length increments
                                           !     for all lines in the I and J
                                           !     directions - see PARAMXY.
                                           !     If from a grid larger than the
                                           !     active subgrid, S0 is trans-
                                           !     formed here.

      REAL    DEDGEI (2, JMIN:JMAX, 2),    !  S  Storage for edge perturbations:
     >        DEDGEJ (2, IMIN:IMAX, 2)     !     DEDGEI (1:2, J1:J2, 1) = dX, dY
                                           !     along the I = I1 edge, etc.

      REAL    X (IMIN:IMAX, JMIN:JMAX),    ! I/O Grid coordinates: input with
     >        Y (IMIN:IMAX, JMIN:JMAX)     !     the edges perturbed; output
                                           !     fully perturbed.
C-------------------------------------------------------------------------------

C     Local constants.

      REAL    EPS, ONE

      PARAMETER (EPS = 1.E-8, ONE = 1.E+0) ! EPS safeguards a divide by zero -
                                           ! presumably only if result is zero.
C     Local variables.

      INTEGER I, J
      REAL    DELI, DELJ, WTI1, WTI2, WTJ1, WTJ2

      LOGICAL SAME

C     Execution.
C     ----------

C     Have any corners moved?  Set up the perturbations in available
C     storage because they are needed anyway if the corners have moved.

      SAME = .TRUE.
      I = I1

      DO J = 1, 2
         DEDGEI (1, J1, J) = X (I, J1) - X0 (I, J1)
         DEDGEI (2, J1, J) = Y (I, J1) - Y0 (I, J1)
         DEDGEI (1, J2, J) = X (I, J2) - X0 (I, J2)
         DEDGEI (2, J2, J) = Y (I, J2) - Y0 (I, J2)

         SAME = DEDGEI (1, J1, J) .EQ. 0. .AND. 
     >          DEDGEI (2, J1, J) .EQ. 0. .AND.
     >          DEDGEI (1, J2, J) .EQ. 0. .AND. 
     >          DEDGEI (2, J2, J) .EQ. 0. .AND. SAME
         I = I2
      END DO


      IF (.NOT. SAME) THEN  ! Corners have moved

C        Set up intermediate edges with the final corners but
C        otherwise derived from the original edges.  Actually,
C        all we need are the intermediate edge perturbations.

C        I = I1 and I2 intermediate edges:

         DO J = J1 + 1, J2 - 1
            WTJ2 = (S0 (I1, J,  2) - S0 (I1, J1, 2)) /
     >             (S0 (I1, J2, 2) - S0 (I1, J1, 2))
            WTJ1 = ONE - WTJ2
            DEDGEI (1, J, 1) = WTJ1 * DEDGEI (1, J1, 1) +
     >                         WTJ2 * DEDGEI (1, J2, 1)
            DEDGEI (2, J, 1) = WTJ1 * DEDGEI (2, J1, 1) +
     >                         WTJ2 * DEDGEI (2, J2, 1)

            WTJ2 = (S0 (I2, J,  2) - S0 (I2, J1, 2)) /
     >             (S0 (I2, J2, 2) - S0 (I2, J1, 2))
            WTJ1 = ONE - WTJ2
            DEDGEI (1, J, 2) = WTJ1 * DEDGEI (1, J1, 2) +
     >                         WTJ2 * DEDGEI (1, J2, 2)
            DEDGEI (2, J, 2) = WTJ1 * DEDGEI (2, J1, 2) +
     >                         WTJ2 * DEDGEI (2, J2, 2)
         END DO

C        J = J1 and J2 intermediate edges:

         DO I = I1 + 1, I2 - 1
            WTI2 = (S0 (I,  J1, 1) - S0 (I1, J1, 1)) /
     >             (S0 (I2, J1, 1) - S0 (I1, J1, 1))
            WTI1 = ONE - WTI2
            DEDGEJ (1, I, 1) = WTI1 * DEDGEI (1, J1, 1) +
     >                         WTI2 * DEDGEI (1, J1, 2)
            DEDGEJ (2, I, 1) = WTI1 * DEDGEI (2, J1, 1) +
     >                         WTI2 * DEDGEI (2, J1, 2)

            WTI2 = (S0 (I,  J2, 1) - S0 (I1, J2, 1)) /
     >             (S0 (I2, J2, 1) - S0 (I1, J2, 1))
            WTI1 = ONE - WTI2
            DEDGEJ (1, I, 2) = WTI1 * DEDGEI (1, J2, 1) +
     >                         WTI2 * DEDGEI (1, J2, 2)
            DEDGEJ (2, I, 2) = WTI1 * DEDGEI (2, J2, 1) +
     >                         WTI2 * DEDGEI (2, J2, 2)
         END DO


C        Perturb the interior based on the original interior and
C        compatible intermediate edges connecting the moved corners.
C        The contributions from each index direction are NOT
C        independent, so they are combined as a weighted average.

         DO J = J1 + 1, J2 - 1
            DO I = I1 + 1, I2 - 1
               WTI2 = (S0 (I,  J, 1) - S0 (I1, J, 1)) /
     >                (S0 (I2, J, 1) - S0 (I1, J, 1))
               WTI1 = ONE - WTI2
               WTJ2 = (S0 (I, J,  2) - S0 (I, J1, 2)) /
     >                (S0 (I, J2, 2) - S0 (I, J1, 2))
               WTJ1 = ONE - WTJ2

               DELI = WTI1 * DEDGEI (1, J, 1) + WTI2 * DEDGEI (1, J, 2)
               DELJ = WTJ1 * DEDGEJ (1, I, 1) + WTJ2 * DEDGEJ (1, I, 2)

               X (I, J) = X0 (I, J) +
     >            (ABS (DELI) * DELI + ABS (DELJ) * DELJ) /
     >             MAX (ABS (DELI) + ABS (DELJ), EPS)

               DELI = WTI1 * DEDGEI (2, J, 1) + WTI2 * DEDGEI (2, J, 2)
               DELJ = WTJ1 * DEDGEJ (2, I, 1) + WTJ2 * DEDGEJ (2, I, 2)

               Y (I, J) = Y0 (I, J) +
     >            (ABS (DELI) * DELI + ABS (DELJ) * DELJ) /
     >             MAX (ABS (DELI) + ABS (DELJ), EPS)
            END DO
         END DO

C        Set up the second-phase edge perturbations from the
C        intermediate edges to the final edges:

C        I = I1 and I2 edges:

         DO J = J1 + 1, J2 - 1
            DEDGEI (1, J, 1) = X (I1, J) - X0 (I1, J) - DEDGEI (1, J, 1)
            DEDGEI (2, J, 1) = Y (I1, J) - Y0 (I1, J) - DEDGEI (2, J, 1)
            DEDGEI (1, J, 2) = X (I2, J) - X0 (I2, J) - DEDGEI (1, J, 2)
            DEDGEI (2, J, 2) = Y (I2, J) - Y0 (I2, J) - DEDGEI (2, J, 2)
         END DO

C        J = J1 and J2 edges:

         DO I = I1 + 1, I2 - 1
            DEDGEJ (1, I, 1) = X (I, J1) - X0 (I, J1) - DEDGEJ (1, I, 1)
            DEDGEJ (2, I, 1) = Y (I, J1) - Y0 (I, J1) - DEDGEJ (2, I, 1)
            DEDGEJ (1, I, 2) = X (I, J2) - X0 (I, J2) - DEDGEJ (1, I, 2)
            DEDGEJ (2, I, 2) = Y (I, J2) - Y0 (I, J2) - DEDGEJ (2, I, 2)
         END DO

      ELSE

C        Just set up the edge perturbations for a single phase.

C        I = I1 and I2 edges:

         DO J = J1 + 1, J2 - 1
            DEDGEI (1, J, 1) = X (I1, J) - X0 (I1, J)
            DEDGEI (2, J, 1) = Y (I1, J) - Y0 (I1, J)
            DEDGEI (1, J, 2) = X (I2, J) - X0 (I2, J)
            DEDGEI (2, J, 2) = Y (I2, J) - Y0 (I2, J)
         END DO

C        J = J1 and J2 edges:

         DO I = I1 + 1, I2 - 1
            DEDGEJ (1, I, 1) = X (I, J1) - X0 (I, J1)
            DEDGEJ (2, I, 1) = Y (I, J1) - Y0 (I, J1)
            DEDGEJ (1, I, 2) = X (I, J2) - X0 (I, J2)
            DEDGEJ (2, I, 2) = Y (I, J2) - Y0 (I, J2)
         END DO

C        Transfer the original interior to avoid duplicating code below:

         DO J = J1 + 1, J2 - 1
            DO I = I1 + 1, I2 - 1
               X (I, J) = X0 (I, J)
               Y (I, J) = Y0 (I, J)
            END DO
         END DO

      END IF


C     Update the interior with contributions from all edges.  The two
C     directions are independent, so just accumulate the contributions.

      DO J = J1 + 1, J2 - 1
         DO I = I1 + 1, I2 - 1
            WTI2 = (S0 (I,  J, 1) - S0 (I1, J, 1)) /
     >             (S0 (I2, J, 1) - S0 (I1, J, 1))
            WTI1 = ONE - WTI2
            WTJ2 = (S0 (I, J,  2) - S0 (I, J1, 2)) /
     >             (S0 (I, J2, 2) - S0 (I, J1, 2))
            WTJ1 = ONE - WTJ2

            DELI = WTI1 * DEDGEI (1, J, 1) + WTI2 * DEDGEI (1, J, 2)
            DELJ = WTJ1 * DEDGEJ (1, I, 1) + WTJ2 * DEDGEJ (1, I, 2)

            X (I, J) = X (I, J) + DELI + DELJ

            DELI = WTI1 * DEDGEI (2, J, 1) + WTI2 * DEDGEI (2, J, 2)
            DELJ = WTJ1 * DEDGEJ (2, I, 1) + WTJ2 * DEDGEJ (2, I, 2)

            Y (I, J) = Y (I, J) + DELI + DELJ
         END DO
      END DO

      RETURN
      END
