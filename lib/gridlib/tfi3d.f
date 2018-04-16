C+------------------------------------------------------------------------------
C
      SUBROUTINE TFI3D (IDIM, JDIM, KDIM, I1, I2, J1, J2, K1, K2,
     >                  X, Y, Z, S, BF, PRINT, ITMAX, DMAX)
C
C     Purpose:
C
C        TFI3D generates a volume grid given the six boundary surface
C     meshes, using transfinite interpolation with Soni-type optimal
C     blending functions.
C
C     Method:
C
C        The method follows Soni's paper "Two- and Three-Dimensional Grid
C     Generation for Internal Flow Applications of CFD", AIAA 85-1526.
C     The final 3-stage algorithm is best documented in the Projectors
C     section of "Numerical Grid Generation, Foundations and Applications"
C     by Thompson, Warsi, and Mastin.
C
C     History:
C
C     09/05/96     SJE      Initial implementation.
C     09/10-13/96  SJE/DAS  Various efficiency refinements, e.g.,
C                           single pass through the volume in TFINT3F.
C                           Tried interpolating normalized face arc-
C                           lengths into the volume first then using
C                           them as blending functions for the volume
C                           X,Y,Zs, but lines may cross.  Would
C                           interpolating UNnormalized face arc-lengths
C                           into the interior first be better?
C     05/22/98     DAS      A bad subscript in initialization of the
C                           volume blending functions showed up for a
C                           Navier-Stokes grid; Fortran 90 updates.
C     05/27/98      "       ITMAX = 0 case works now.  True face arc
C                           lengths (normalized) are used now.  One or
C                           two more iterations are needed as a result
C                           but it's cleaner than use of edges only.
C     05/28/98     SJE/DAS  Another index error was in the iteration,
C                           though the effect has been minor.  Revised
C                           nomenclature reveals no further glitches.
C                           Behavior is still poor for a subblock above
C                           a wing wake when the fuselage is present,
C                           but the GRIDGEN result is equally poor too.
C                           Too large an increment off the body crown/
C                           keel is a contributing factor, but results
C                           are still not well understood.
C
C     Author:  Stephen J. Edwards, NASA Ames, Mtn. View, CA.
C ------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) ::
     >   IDIM, JDIM, KDIM,       ! Max. block dims. in the calling program
     >   I1, I2, J1, J2, K1, K2, ! Active block boundary indices
     >   ITMAX                   ! Max. optimal blending iterations (e.g., 10)

      REAL, DIMENSION (IDIM,JDIM,KDIM), INTENT (INOUT) ::
     >   X, Y, Z                 ! Input with initial (sub)surface grids;
                                 ! (sub)volume grid filled upon return 

      REAL, INTENT (OUT) ::
     >   S(IDIM,JDIM,KDIM,3),    ! Normalized arcs in each direction (scratch)
     >   BF(IDIM,JDIM,KDIM,3)    ! 3D blending functions

      LOGICAL, INTENT (IN) ::
     >   PRINT                   ! .TRUE. shows the optimal blending function
                                 ! iterations on unit 6

      REAL, INTENT (IN) ::
     >   DMAX                    ! Convergence tol. for the Jacobi iteration
                                 ! on optimal volume blending fns., which are
                                 ! in [0,1]; derive it from the smallest
                                 ! spacing as a fraction of the data range

C-------------------------------------------------------------------------------

C     Local constants:

      REAL, PARAMETER :: HALF = 0.5, ONE = 1., ZERO = 0.

C     Local variables:

      INTEGER
     >   I, J, K, IT
      REAL
     >   DJ1, DK1, DV1,
     >   DK2, DI2, DV2,
     >   DI3, DJ3, DV3, P, Q, R, PM1, QM1, RM1

C     Execution:

C     Calculate normalized arc length distributions on each face:

      CALL PARAMXYZ (1, IDIM, 1, JDIM, 1, KDIM, I1, I2, J1, J2,
     >               K1, K1, X, Y, Z, S)

      CALL PARAMXYZ (1, IDIM, 1, JDIM, 1, KDIM, I1, I2, J1, J2,
     >               K2, K2, X, Y, Z, S)

      CALL PARAMXYZ (1, IDIM, 1, JDIM, 1, KDIM, I1, I2, J1, J1,
     >               K1, K2, X, Y, Z, S)

      CALL PARAMXYZ (1, IDIM, 1, JDIM, 1, KDIM, I1, I2, J2, J2,
     >               K1, K2, X, Y, Z, S)

      CALL PARAMXYZ (1, IDIM, 1, JDIM, 1, KDIM, I1, I1, J1, J2,
     >               K1, K2, X, Y, Z, S)

      CALL PARAMXYZ (1, IDIM, 1, JDIM, 1, KDIM, I2, I2, J1, J2,
     >               K1, K2, X, Y, Z, S)

C     PARAMXYZ puts in 1s for the degenerate directions, so:

      S(I1,J1:J2,K1:K2,1) = ZERO
      S(I1:I2,J1,K1:K2,2) = ZERO
      S(I1:I2,J1:J2,K1,3) = ZERO

C     Initialize the 3D volume blending functions by applying the 2D algorithm
C     to I, J, K planes.  A sweep in one direction gives estimates for 2 of the
C     3 blending functions at a point. For efficiency, the 3 sweeps and the
C     averaging are performed together.

      DO K = K1 + 1, K2 - 1

         DO J = J1 + 1, J2 - 1

            DJ3 = S(I2,J,K,2) - S(I1,J,K,2) ! Between I edges/J dirn./K plane
            DK2 = S(I2,J,K,3) - S(I1,J,K,3) !         I       K       J

            DO I = I1 + 1, I2 - 1

               DK1 = S(I,J2,K,3) - S(I,J1,K,3)
               DI3 = S(I,J2,K,1) - S(I,J1,K,1)
               DI2 = S(I,J,K2,1) - S(I,J,K1,1)
               DJ1 = S(I,J,K2,2) - S(I,J,K1,2)

               DV1 = ONE / (ONE - DJ1*DK1)
               DV2 = ONE / (ONE - DK2*DI2)
               DV3 = ONE / (ONE - DI3*DJ3)

               S(I,J,K,1) = ((S(I,J,K1,1) + S(I1,J,K,3)*DI2)*DV2  +
     >                       (S(I,J1,K,1) + S(I1,J,K,2)*DI3)*DV3) * HALF
               S(I,J,K,2) = ((S(I,J,K1,2) + S(I,J1,K,3)*DJ1)*DV1  +
     >                       (S(I1,J,K,2) + S(I,J1,K,1)*DJ3)*DV3) * HALF
               S(I,J,K,3) = ((S(I,J1,K,3) + S(I,J,K1,2)*DK1)*DV1  +
     >                       (S(I1,J,K,3) + S(I,J,K1,1)*DK2)*DV2) * HALF
            END DO
         END DO
      END DO


C     Jacobi iteration to satisfy Soni's equations A.5 - A.7.
C     BF & S(,,,1:3) end up with the optimal volume blending functions:

      IF (PRINT) THEN
         I = I2 - I1 + 1
         J = J2 - J1 + 1
         K = K2 - K1 + 1
         WRITE (6,1000) I, J, K, ITMAX, DMAX
      END IF

      BF(I1:I2,J1:J2,K1:K2,1:3) = S(I1:I2,J1:J2,K1:K2,1:3)

      DO IT = 1, ITMAX ! Or until convergence

         DI2 = ZERO
         DJ1 = ZERO
         DK1 = ZERO

         DO K = K1, K2
            DO J = J1, J2
               DO I = I1, I2

                  P = BF(I,J,K,1)
                  Q = BF(I,J,K,2)
                  R = BF(I,J,K,3)
                  PM1 = ONE - P
                  QM1 = ONE - Q
                  RM1 = ONE - R

                  S(I,J,K,1) = RM1*(QM1*BF(I,J1,K1,1) + Q*BF(I,J2,K1,1))
     >                         + R*(QM1*BF(I,J1,K2,1) + Q*BF(I,J2,K2,1))

                  S(I,J,K,2) = PM1*(RM1*BF(I1,J,K1,2) + R*BF(I1,J,K2,2))
     >                         + P*(RM1*BF(I2,J,K1,2) + R*BF(I2,J,K2,2))

                  S(I,J,K,3) = QM1*(PM1*BF(I1,J1,K,3) + P*BF(I2,J1,K,3))
     >                         + Q*(PM1*BF(I1,J2,K,3) + P*BF(I2,J2,K,3))

                  DI2 = MAX (ABS (BF(I,J,K,1) - S(I,J,K,1)), DI2)
                  DJ1 = MAX (ABS (BF(I,J,K,2) - S(I,J,K,2)), DJ1)
                  DK1 = MAX (ABS (BF(I,J,K,3) - S(I,J,K,3)), DK1)

                  BF(I,J,K,1:3) = S(I,J,K,1:3) ! More cache efficient than
                                               ! after the loops, and
               END DO                          ! equivalent at convergence
            END DO
         END DO

         IF (PRINT) WRITE (6, '(I5, 3F15.8)') IT, DI2, DJ1, DK1

         IF (MAX (DI2, DJ1, DK1) < DMAX) THEN
            EXIT
         ELSE IF (IT == ITMAX) THEN
            IF (.NOT. PRINT) WRITE (6, '(/, A, 3E12.5)')
     >         ' TFI3D iteration limit. Final maximum corrections:',
     >         DI2, DJ1, DK1
         END IF

      END DO ! Next iteration


C     The 3D interpolation steps appear on p. 326 of Thompson, Warsi, Mastin.
C     The S arguments are more than enough work-space.

      CALL TFINT3F (1, IDIM, 1, JDIM, 1, KDIM, I1, I2, J1, J2, K1, K2,
     >              BF, S(1,1,1,1), S(1,1,1,2), S(1,1,1,3), X, Y, Z)

      RETURN

 1000 FORMAT (/, ' TFI3D mesh size:', I5, '  x', I4, '  x', I4,
     >   '   Optimal blending iterations:', I4,
     >   '   Max. correction tolerance:', F10.6, //,
     >   ' ITER', 10X, 'DIMAX', 10X, 'DJMAX', 10X, 'DKMAX')

      END SUBROUTINE TFI3D
