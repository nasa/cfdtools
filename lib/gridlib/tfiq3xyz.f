C+------------------------------------------------------------------------------
C
      SUBROUTINE TFIQ3XYZ (IDIM, JDIM, KDIM, I1, I2, J1, J2, K1, K2,
     >                     X, Y, Z, ARCS)
C
C     Purpose:
C
C        TFIQ3XYZ fills the interior of one (sub)-plane of a 3-space grid
C     block given the four edges of the plane, using transfinite interpolation
C     with Soni-type blending functions.
C
C        This is more general than the earlier TFIQ3D, which is limited to
C     "K" planes only.  It also allows for input of the normalized edge arc
C     lengths in case they are already available.
C
C     History:
C
C     06/15/97   D.Saunders   Initial adaptation of parts of TFI3D & TFIQ3D.
C
C ------------------------------------------------------------------------------

C     Arguments:

      IMPLICIT NONE

      INTEGER
     >   IDIM, JDIM, KDIM,      ! I   Max. block dims. in the calling program
     >   I1, I2, J1, J2, K1, K2 ! I   Active sub-plane indices, two being equal

      REAL
     >   X(IDIM,JDIM,KDIM),     ! I/O Grid coordinates, input at specified
     >   Y(IDIM,JDIM,KDIM),     !     edges, output with corresponding interior
     >   Z(IDIM,JDIM,KDIM),     !     filled
     >   ARCS(IDIM,JDIM,KDIM,3) ! I/O Normalized arcs in each direction, needed
                                !     on the relevant edges;
                                !     ARCS(I1,J1,K1,n) < 0. on input forces
                                !     computing them here, where n = 2, 3, 1
                                !     for an I, J, K plane, respectively
                                !     (NOT 1, 2, 3 because that would imply
                                !     testing in the degenerate direction)
C-------------------------------------------------------------------------------

C     Local constants:

      REAL
     >   ONE, ZERO
      PARAMETER
     >  (ONE = 1., ZERO = 0.)

C     Local variables:

      INTEGER
     >   I1P1, I2M1, J1P1, J2M1, K1P1, K2M1, I, J, K
      REAL
     >   DI, DJ, DK, DV, F1, F2, R1MF1, R1MF2

C     Execution:

      I1P1 = I1 + 1
      I2M1 = I2 - 1
      J1P1 = J1 + 1
      J2M1 = J2 - 1
      K1P1 = K1 + 1
      K2M1 = K2 - 1

      IF (I1 .EQ. I2) THEN  ! I plane case

         I = I1
         IF (ARCS(I,J1,K1,2) .LT. ZERO) THEN  ! Note 2, not 1

C           Calculate normalized arc length distributions on each edge:

            CALL PARAMXYZ (1, IDIM, 1, JDIM, 1, KDIM, I, I, J1, J2,
     >                     K1, K1, X, Y, Z, ARCS)     !  K1 edge

            CALL PARAMXYZ (1, IDIM, 1, JDIM, 1, KDIM, I, I, J1, J2,
     >                     K2, K2, X, Y, Z, ARCS)     !  K2 edge

            CALL PARAMXYZ (1, IDIM, 1, JDIM, 1, KDIM, I, I, J1, J1,
     >                     K1, K2, X, Y, Z, ARCS)     !  J1 edge

            CALL PARAMXYZ (1, IDIM, 1, JDIM, 1, KDIM, I, I, J2, J2,
     >                     K1, K2, X, Y, Z, ARCS)     !  J2 edge
         END IF

         DO K = K1P1, K2M1
            DK = ARCS(I,J2,K,3) - ARCS(I,J1,K,3)
            DO J = J1P1, J2M1
               DJ = ARCS(I,J,K2,2) - ARCS(I,J,K1,2)
               DV = ONE / (ONE - DJ*DK)

               F1 = (ARCS(I,J,K1,2) + ARCS(I,J1,K,3)*DJ) * DV
               F2 = (ARCS(I,J1,K,3) + ARCS(I,J,K1,2)*DK) * DV

               R1MF1 = ONE - F1
               R1MF2 = ONE - F2

               X(I,J,K) = R1MF1 * X(I,J1,K) + F1 * X(I,J2,K)
     >                  + R1MF2 * X(I,J,K1) + F2 * X(I,J,K2)
     >                  - R1MF2 * (R1MF1 * X(I,J1,K1) + F1 * X(I,J2,K1))
     >                  - F2    * (R1MF1 * X(I,J1,K2) + F1 * X(I,J2,K2))
 
               Y(I,J,K) = R1MF1 * Y(I,J1,K) + F1 * Y(I,J2,K)
     >                  + R1MF2 * Y(I,J,K1) + F2 * Y(I,J,K2)
     >                  - R1MF2 * (R1MF1 * Y(I,J1,K1) + F1 * Y(I,J2,K1))
     >                  - F2    * (R1MF1 * Y(I,J1,K2) + F1 * Y(I,J2,K2))
 
               Z(I,J,K) = R1MF1 * Z(I,J1,K) + F1 * Z(I,J2,K)
     >                  + R1MF2 * Z(I,J,K1) + F2 * Z(I,J,K2)
     >                  - R1MF2 * (R1MF1 * Z(I,J1,K1) + F1 * Z(I,J2,K1))
     >                  - F2    * (R1MF1 * Z(I,J1,K2) + F1 * Z(I,J2,K2))
            END DO
         END DO

      ELSE IF (J1 .EQ. J2) THEN  ! J plane case

         J = J1
         IF (ARCS(I1,J,K1,3) .LT. ZERO) THEN

C           Calculate normalized arc length distributions on each edge:

            CALL PARAMXYZ (1, IDIM, 1, JDIM, 1, KDIM, I1, I2, J, J,
     >                     K1, K1, X, Y, Z, ARCS)     !  K1 edge

            CALL PARAMXYZ (1, IDIM, 1, JDIM, 1, KDIM, I1, I2, J, J,
     >                     K2, K2, X, Y, Z, ARCS)     !  K2 edge

            CALL PARAMXYZ (1, IDIM, 1, JDIM, 1, KDIM, I1, I1, J, J,
     >                     K1, K2, X, Y, Z, ARCS)     !  I1 edge

            CALL PARAMXYZ (1, IDIM, 1, JDIM, 1, KDIM, I2, I2, J, J,
     >                     K1, K2, X, Y, Z, ARCS)     !  I2 edge
         END IF

         DO K = K1P1, K2M1
            DK = ARCS(I2,J,K,3) - ARCS(I1,J,K,3)
            DO I = I1P1, I2M1
               DI = ARCS(I,J,K2,1) - ARCS(I,J,K1,1)
               DV = ONE / (ONE - DI*DK)

               F1 = (ARCS(I,J,K1,1) + ARCS(I1,J,K,3)*DI) * DV
               F2 = (ARCS(I1,J,K,3) + ARCS(I,J,K1,1)*DK) * DV

               R1MF1 = ONE - F1
               R1MF2 = ONE - F2

               X(I,J,K) = R1MF1 * X(I1,J,K) + F1 * X(I2,J,K)
     >                  + R1MF2 * X(I,J,K1) + F2 * X(I,J,K2)
     >                  - R1MF2 * (R1MF1 * X(I1,J,K1) + F1 * X(I2,J,K1))
     >                  - F2    * (R1MF1 * X(I1,J,K2) + F1 * X(I2,J,K2))
 
               Y(I,J,K) = R1MF1 * Y(I1,J,K) + F1 * Y(I2,J,K)
     >                  + R1MF2 * Y(I,J,K1) + F2 * Y(I,J,K2)
     >                  - R1MF2 * (R1MF1 * Y(I1,J,K1) + F1 * Y(I2,J,K1))
     >                  - F2    * (R1MF1 * Y(I1,J,K2) + F1 * Y(I2,J,K2))
 
               Z(I,J,K) = R1MF1 * Z(I1,J,K) + F1 * Z(I2,J,K)
     >                  + R1MF2 * Z(I,J,K1) + F2 * Z(I,J,K2)
     >                  - R1MF2 * (R1MF1 * Z(I1,J,K1) + F1 * Z(I2,J,K1))
     >                  - F2    * (R1MF1 * Z(I1,J,K2) + F1 * Z(I2,J,K2))
            END DO
         END DO

      ELSE IF (K1 .EQ. K2) THEN  ! K plane case

         K = K1
         IF (ARCS(I1,J1,K,1) .LT. ZERO) THEN

C           Calculate normalized arc length distributions on each edge:

            CALL PARAMXYZ (1, IDIM, 1, JDIM, 1, KDIM, I1, I2, J1, J1,
     >                     K, K, X, Y, Z, ARCS)       !  J1 edge

            CALL PARAMXYZ (1, IDIM, 1, JDIM, 1, KDIM, I1, I2, J2, J2,
     >                     K, K, X, Y, Z, ARCS)       !  J2 edge

            CALL PARAMXYZ (1, IDIM, 1, JDIM, 1, KDIM, I1, I1, J1, J2,
     >                     K, K, X, Y, Z, ARCS)       !  I1 edge

            CALL PARAMXYZ (1, IDIM, 1, JDIM, 1, KDIM, I2, I2, J1, J2,
     >                     K, K, X, Y, Z, ARCS)       !  I2 edge
         END IF

         DO J = J1P1, J2M1
            DJ = ARCS(I2,J,K,2) - ARCS(I1,J,K,2)
            DO I = I1P1, I2M1
               DI = ARCS(I,J2,K,1) - ARCS(I,J1,K,1)
               DV = ONE / (ONE - DI*DJ)

               F1 = (ARCS(I,J1,K,1) + ARCS(I1,J,K,2)*DI) * DV
               F2 = (ARCS(I1,J,K,2) + ARCS(I,J1,K,1)*DJ) * DV

               R1MF1 = ONE - F1
               R1MF2 = ONE - F2

               X(I,J,K) = R1MF1 * X(I1,J,K) + F1 * X(I2,J,K)
     >                  + R1MF2 * X(I,J1,K) + F2 * X(I,J2,K)
     >                  - R1MF2 * (R1MF1 * X(I1,J1,K) + F1 * X(I2,J1,K))
     >                  - F2    * (R1MF1 * X(I1,J2,K) + F1 * X(I2,J2,K))
 
               Y(I,J,K) = R1MF1 * Y(I1,J,K) + F1 * Y(I2,J,K)
     >                  + R1MF2 * Y(I,J1,K) + F2 * Y(I,J2,K)
     >                  - R1MF2 * (R1MF1 * Y(I1,J1,K) + F1 * Y(I2,J1,K))
     >                  - F2    * (R1MF1 * Y(I1,J2,K) + F1 * Y(I2,J2,K))
 
               Z(I,J,K) = R1MF1 * Z(I1,J,K) + F1 * Z(I2,J,K)
     >                  + R1MF2 * Z(I,J1,K) + F2 * Z(I,J2,K)
     >                  - R1MF2 * (R1MF1 * Z(I1,J1,K) + F1 * Z(I2,J1,K))
     >                  - F2    * (R1MF1 * Z(I1,J2,K) + F1 * Z(I2,J2,K))
             END DO
         END DO

      END IF

      RETURN
      END
