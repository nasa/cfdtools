C+----------------------------------------------------------------------
C
      SUBROUTINE TFI2D (IDIM, I1, I2, J1, J2, X, Y, ARCS)
C
C  ACRONYM:  TransFinite Interpolation (2D grid generation)
C            -    -      -              --
C  PURPOSE:
C
C        TFI2D generates interior points of a regular mesh in 2-space
C     given the points along all four edges, using a form of transfinite
C     interpolation.  This version uses Soni blending functions.
C     It provides for filling any subregion of the plane.  Any one of
C     the four edges may be collapsed.  (Originally, only the I1 or J1
C     edges were allowed to be degenerate.)  Normalized arcs from the
C     opposite edge are used in such cases.
C
C  METHOD:
C
C        References for transfinite interpolation are given below.
C     "Transfinite" refers to matching of a function on a boundary at
C     a nondenumerable number of points.  The algebra used here may be
C     summarized as follows:
C
C        A typical point (I,J) defines a 4-point stencil of boundary
C     points:  (I1,J), (I2,J), and (I,J1), (I,J2).  Weights p, 1-p are
C     applied to the first pair, and q, 1-q to the second pair, where
C     p = p(I) and q = q(J) are appropriate functions of arc-lengths
C     along opposite boundaries.  A 4-point tensor term is also formed
C     as for standard 2-D linear interpolation, using the "corners"
C     with p, q, 1-p, and 1-q.  The final combination is of the form
C     F(p) + F(q) - F(p,q)  where F is X or Y.
C
C  REFERENCES: "Numerical Grid Generation: Foundations and Applications"
C              Thompson, Warsi, and Martin (1985)
C              "Two- and Three-dimensional Grid Generation for Internal
C              Flow Applications of CFD", AIAA 85-1526
C
C  HISTORY:
C
C   03/01/88    JEM      Initial implementation for (1:N,1:M) and two
C                        cases of degeneracy (bottom or left).
C   04/01/88  RGL/DAS    Provided for in-place filling of (I1:I2,J1:J2);
C                        single work-space argument to simplify calling
C                        program.
C   05/11/95    DAS      Factorized main expressions a little more.
C   05/19/95    DAS      Work-space offsets were wrong if I1 or J1 > 1.
C   09/03/96  S.Edwards  Switched from linear blending functions to the
C                        Soni blending.
C   06/21/97    DAS      DI, DJ definitions now match Soni (not -DI, -DJ).
C   09/07/99     "       Mark Rimlinger suggested taking advantage of
C                        Fortran 90 to improve the nomenclature, & handle
C                        all 4 possible degeneracies instead of just 2,
C                        while retaining the same calling sequence to
C                        avoid affecting numerous existing applications.
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
     >   X, Y               ! Grid coordinates for 2D (sub)region:
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
     >      ((X(I1,J) - X(I1,J-1)) ** 2 + (Y(I1,J) - Y(I1,J-1)) ** 2)
         SI2(NJ) = SI2(NJ-1) + SQRT
     >      ((X(I2,J) - X(I2,J-1)) ** 2 + (Y(I2,J) - Y(I2,J-1)) ** 2)
      END DO

      SJ1(1) = ZERO
      SJ2(1) = ZERO
      NI = 1

      DO I = I1 + 1, I2
         NI = NI + 1
         SJ1(NI) = SJ1(NI-1) + SQRT
     >      ((X(I,J1) - X(I-1,J1)) ** 2 + (Y(I,J1) - Y(I-1,J1)) ** 2)
         SJ2(NI) = SJ2(NI-1) + SQRT
     >      ((X(I,J2) - X(I-1,J2)) ** 2 + (Y(I,J2) - Y(I-1,J2)) ** 2)
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
         END DO

      END DO

      END SUBROUTINE TFI2D
