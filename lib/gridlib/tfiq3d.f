C+----------------------------------------------------------------------
C
      SUBROUTINE TFIQ3D (IDIM, I1, I2, J1, J2, X, Y, Z, ARCS)
C
C  ACRONYM:  TransFinite Interpolation (Quasi-3D grid generation)
C            -    -      -              -     --
C  PURPOSE:
C
C        TFIQ3D generates interior points of a regular surface mesh
C     in 3-space given the points along all four edges, using a form
C     of transfinite interpolation.  This is the analogue of TFI2D, q.v.
C
C  HISTORY:
C
C   04/20/88  Ron Langhi      Initial extension of TFI2D.
C   09/08/99  David Saunders  Analogue of Fortran 90 version of TFI2D.
C
C  AUTHOR:  John E. Melton, NASA/Ames Research Center, Moffett Field, CA
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) ::
     >   IDIM,              ! "Row" dimension in calling program.
     >   I1, I2,            ! I1:I2 is the active "row" size, where
     >                      ! 1 <= I1 and I1 + 2 <= I2 <= IDIM.
     >   J1, J2             ! J1:J2 is the active "column" size, where
                            ! 1 <= J1 and J1 + 2 <= J2.

      REAL, INTENT (INOUT), DIMENSION (IDIM,J2) ::
     >   X, Y, Z            ! Grid coordinates for (sub)surface in 3-space:
     >                      ! Input:  desired points along boundaries
                            !         defined by I1, I2, J1, J2;
                            ! Output: interior grid points, with the
                            !         boundary points (and beyond) unchanged.
                            ! Any one of the active edges may be degenerate.

      REAL, INTENT (OUT), TARGET ::
     >   ARCS (2*(I2 + J2)) ! Work-space for normalized boundary arc lengths.

C     Local constants:

      REAL, PARAMETER ::
     >   ONE = 1., TOL = 0.00001, ZERO = 0.

C     Local variables:

      INTEGER
     >   I, J, NI, NJ

      REAL
     >   DI, DJ, DV, F1, F2, R1MF1, R1MF2, RARC

      REAL, POINTER, DIMENSION (:) ::
     >   SI1, SI2, SJ1, SJ2

      LOGICAL
     >   FIRST_COLLAPSED

C     Execution:

C     Partition the work-space:

      J = 2 * J2
      I = J + I2
      SI1 => ARCS(J1:J2)    ! SI1(1:J2-J1+1) <-> ARCS(J1:J2), etc.
      SI2 => ARCS(J2+1:J)
      SJ1 => ARCS(J+1:I)
      SJ2 => ARCS(I+1:I+I2)

C     Compute normalized arc lengths along active boundaries:

      SI1(1) = ZERO
      SI2(1) = ZERO
      NJ = 1

      DO J = J1 + 1, J2
         NJ = NJ + 1
         SI1(NJ) = SI1(NJ-1) + SQRT
     >      ((X(I1,J) - X(I1,J-1)) ** 2 +
     >       (Y(I1,J) - Y(I1,J-1)) ** 2 +
     >       (Z(I1,J) - Z(I1,J-1)) ** 2)
         SI2(NJ) = SI2(NJ-1) + SQRT
     >      ((X(I2,J) - X(I2,J-1)) ** 2 +
     >       (Y(I2,J) - Y(I2,J-1)) ** 2 +
     >       (Z(I2,J) - Z(I2,J-1)) ** 2)
      END DO

      SJ1(1) = ZERO
      SJ2(1) = ZERO
      NI = 1

      DO I = I1 + 1, I2
         NI = NI + 1
         SJ1(NI) = SJ1(NI-1) + SQRT
     >      ((X(I,J1) - X(I-1,J1)) ** 2 +
     >       (Y(I,J1) - Y(I-1,J1)) ** 2 +
     >       (Z(I,J1) - Z(I-1,J1)) ** 2)
         SJ2(NI) = SJ2(NI-1) + SQRT
     >      ((X(I,J2) - X(I-1,J2)) ** 2 +
     >       (Y(I,J2) - Y(I-1,J2)) ** 2 +
     >       (Z(I,J2) - Z(I-1,J2)) ** 2)
      END DO

C     Normalize if possible, else copy the opposite edge after normalizing.
C     Assume only 1 of the 4 edges can degenerate:

      FIRST_COLLAPSED = SI1(NJ) < TOL

      IF (.NOT. FIRST_COLLAPSED) THEN
         RARC    = ONE / SI1(NJ)
         SI1     = SI1 * RARC
         SI1(NJ) = ONE ! Make sure of it
      END IF

      IF (SI2(NJ) > TOL) THEN
         RARC    = ONE / SI2(NJ)
         SI2     = SI2 * RARC
         SI2(NJ) = ONE
      ELSE
         SI2 = SI1
      END IF

      IF (FIRST_COLLAPSED) SI1 = SI2

      FIRST_COLLAPSED = SJ1(NI) < TOL

      IF (.NOT. FIRST_COLLAPSED) THEN
         RARC    = ONE / SJ1(NI)
         SJ1     = SJ1 * RARC
         SJ1(NI) = ONE
      END IF

      IF (SJ2(NI) > TOL) THEN
         RARC    = ONE / SJ2(NI)
         SJ2     = SJ2 * RARC
         SJ2(NI) = ONE
      ELSE
         SJ2 = SJ1
      END IF

      IF (FIRST_COLLAPSED) SJ1 = SJ2

C     Interpolate from the edges to the interior:

      NJ = 1

      DO J = J1 + 1, J2 - 1

         NJ = NJ + 1
         DJ = SI2(NJ) - SI1(NJ)
         NI = 1

         DO I = I1 + 1, I2 - 1

            NI = NI + 1
            DI = SJ2(NI) - SJ1(NI)
            DV = ONE / (ONE - DI * DJ)
            F1 = (SJ1(NI) + SI1(NJ) * DI) * DV
            F2 = (SI1(NJ) + SJ1(NI) * DJ) * DV

            R1MF1 = ONE - F1
            R1MF2 = ONE - F2

            X(I,J) =   R1MF1 * X(I1,J) + F1 * X(I2,J)
     >               + R1MF2 * X(I,J1) + F2 * X(I,J2)
     >               - R1MF2 * (R1MF1 * X(I1,J1) + F1 * X(I2,J1))
     >               - F2    * (R1MF1 * X(I1,J2) + F1 * X(I2,J2))

            Y(I,J) =   R1MF1 * Y(I1,J) + F1 * Y(I2,J)
     >               + R1MF2 * Y(I,J1) + F2 * Y(I,J2)
     >               - R1MF2 * (R1MF1 * Y(I1,J1) + F1 * Y(I2,J1))
     >               - F2    * (R1MF1 * Y(I1,J2) + F1 * Y(I2,J2))

            Z(I,J) =   R1MF1 * Z(I1,J) + F1 * Z(I2,J)
     >               + R1MF2 * Z(I,J1) + F2 * Z(I,J2)
     >               - R1MF2 * (R1MF1 * Z(I1,J1) + F1 * Z(I2,J1))
     >               - F2    * (R1MF1 * Z(I1,J2) + F1 * Z(I2,J2))
         END DO

      END DO

      END SUBROUTINE TFIQ3D
