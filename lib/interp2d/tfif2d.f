C+----------------------------------------------------------------------
C
      SUBROUTINE TFIF2D (NF, IDIM, I1, I2, J1, J2, X, Y, F)
C
C  ACRONYM:  TransFinite Interpolation of Function(s) in 2 Dimensions
C            -    -      -                -              - -
C  PURPOSE:
C
C        TFIF2D generates interior points of a regular mesh in 2-space
C     given the points along all four edges, using a form of transfinite
C     interpolation, and interpolates one or more functions other than
C     (x,y) as well (unlike TFI2D, from which it has been adapted).
C
C     [Referring to TFI2D:]  This version uses Soni blending functions.
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
C   03/01/88  J.Melton   Initial implementation for (1:N,1:M) and two
C                        cases of degeneracy (bottom or left).
C   04/01/88  R.Langhi   Provided for in-place filling of (I1:I2,J1:J2);
C             D.Saunders single work-space argument to simplify calling
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
C   02/27/07     "       TFI2DF adapted from TFI2D to deal with one or
C                        more functions defined on the edges of a 2-D
C                        structured grid block.  Note that if Z were
C                        passed as a function, results would not be the
C                        same as those from the TFIQ3D surface grid
C                        utility because the arc lengths used here
C                        involve only X and Y.
C   03/07/07     "       Use of just X and Y for the arc lengths ignores
C                        what TFIQ3D would do if F were Z.  Scaling F is
C                        desirable but arbitrary.  (Too big would swamp
C                        X and Y in the quasi-arc-lengths.)  Choose to
C                        transform each F to the average data range of
C                        X and Y.  Dispense with the ARCS argument,
C                        since we now need other local work space.
C
C  AUTHOR:  John E. Melton, NASA/Ames Research Center, Moffett Field, CA
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) ::
     >   NF,                ! # functions to be interpolated in addition
     >                      ! to X and Y; NF >= 1
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

      REAL, INTENT (INOUT), DIMENSION (NF,IDIM,J2) ::
     >   F                  ! Additional functions defined on the edges as
                            ! for X and Y

C     Local constants:

      REAL, PARAMETER ::
     >   HALF = 0.5, ONE = 1., TOL = 0.00001, ZERO = 0.

C     Local variables:

      INTEGER
     >   I, J, N, NI, NJ

      REAL
     >   DI, DJ, DV, F1, F2, R1MF1, R1MF2, RARC,
     >   FMIN, FMAX, FRANGE, FSCALE, FSHIFT, XMAX, XMIN, YMAX, YMIN

      REAL, POINTER, DIMENSION (:) ::
     >   SI1, SI2, SJ1, SJ2

      REAL, ALLOCATABLE, TARGET, DIMENSION (:) ::
     >   ARCS

      REAL, ALLOCATABLE, DIMENSION (:,:) ::
     >   FNORM

      LOGICAL
     >   FIRST_COLLAPSED

C     Execution:

      ALLOCATE (ARCS (2*(I2 + J2)), FNORM(IDIM,J2))

C     Partition the work-space:

      J = 2 * J2
      I = J + I2
      SI1 => ARCS(J1:J2)    ! SI1(1:J2-J1+1) <-> ARCS(J1:J2), etc.
      SI2 => ARCS(J2+1:J)
      SJ1 => ARCS(J+1:I)
      SJ2 => ARCS(I+1:I+I2)

C     Determine the data ranges for X and Y:

      XMIN = X(I1,J1);  XMAX = XMIN
      YMIN = Y(I1,J1);  YMAX = YMIN

      DO J = J1, J2, J2 - J1
         DO I = I1, I2, I2 - I1
            XMIN = MIN (XMIN, X(I,J))
            XMAX = MAX (XMAX, X(I,J))
            YMIN = MIN (YMIN, Y(I,J))
            YMAX = MAX (YMAX, Y(I,J))
         END DO
      END DO

      YMAX = ((XMAX - XMIN) + (YMAX - YMIN)) * HALF

C     The "arc" lengths differ for each function:

      DO N = 1, NF

         FNORM(:,J1) = F(N,:,J1)  ! Transfer given edge values
         FNORM(:,J2) = F(N,:,J2)
         FNORM(I1,:) = F(N,I1,:)
         FNORM(I2,:) = F(N,I2,:)

         FMIN = FNORM(I1,J1);  FMAX = FMIN

         DO J = J1, J2
            DO I = I1, I2
               IF (I == I1 .OR. I == I2 .OR. J == J1 .OR. J == J2) THEN
                  FMIN = MIN (FMIN, FNORM(I,J))
                  FMAX = MAX (FMAX, FNORM(I,J))
               END IF
            END DO
         END DO

         FRANGE = FMAX - FMIN

         IF (FRANGE < TOL) THEN
             FSCALE = ONE
         ELSE
             FSCALE = YMAX / FRANGE
         END IF

C        Transform this F edge values to [0, average of X and Y data ranges]:

         DO J = J1, J2
            DO I = I1, I2
               IF (I == I1 .OR. I == I2 .OR. J == J1 .OR. J == J2) THEN
                  FNORM(I,J) = (FNORM(I,J) - FMIN) * FSCALE
               END IF
            END DO
         END DO

C        Compute normalized arc lengths along active boundaries:

         SI1(1) = ZERO
         SI2(1) = ZERO
         NJ = 1

         DO J = J1 + 1, J2
            NJ = NJ + 1
            SI1(NJ) = SI1(NJ-1) + SQRT
     >         ((X(I1,J) - X(I1,J-1)) ** 2 + (Y(I1,J) - Y(I1,J-1)) ** 2
     >          + (FNORM(I1,J) - FNORM(I1,J-1)) ** 2)
            SI2(NJ) = SI2(NJ-1) + SQRT
     >         ((X(I2,J) - X(I2,J-1)) ** 2 + (Y(I2,J) - Y(I2,J-1)) ** 2
     >          + (FNORM(I2,J) - FNORM(I2,J-1)) ** 2)
         END DO

         SJ1(1) = ZERO
         SJ2(1) = ZERO
         NI = 1

         DO I = I1 + 1, I2
            NI = NI + 1
            SJ1(NI) = SJ1(NI-1) + SQRT
     >         ((X(I,J1) - X(I-1,J1)) ** 2 + (Y(I,J1) - Y(I-1,J1)) ** 2
     >          + (FNORM(I,J1) - FNORM(I-1,J1)) ** 2)
            SJ2(NI) = SJ2(NI-1) + SQRT
     >         ((X(I,J2) - X(I-1,J2)) ** 2 + (Y(I,J2) - Y(I-1,J2)) ** 2
     >          + (FNORM(I,J2) - FNORM(I-1,J2)) ** 2)
         END DO

C        Normalize if possible, else copy the opposite edge after normalizing.
C        Assume only 1 of the 4 edges can degenerate:

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

C        Interpolate from the edges to the interior:

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
     >                  + R1MF2 * X(I,J1) + F2 * X(I,J2)
     >                  - R1MF2 * (R1MF1 * X(I1,J1) + F1 * X(I2,J1))
     >                  - F2    * (R1MF1 * X(I1,J2) + F1 * X(I2,J2))

               Y(I,J) =   R1MF1 * Y(I1,J) + F1 * Y(I2,J)
     >                  + R1MF2 * Y(I,J1) + F2 * Y(I,J2)
     >                  - R1MF2 * (R1MF1 * Y(I1,J1) + F1 * Y(I2,J1))
     >                  - F2    * (R1MF1 * Y(I1,J2) + F1 * Y(I2,J2))

             FNORM(I,J) = R1MF1 * FNORM(I1,J) + F1 * FNORM(I2,J)
     >                  + R1MF2 * FNORM(I,J1) + F2 * FNORM(I,J2)
     >                  - R1MF2 * (R1MF1*FNORM(I1,J1) + F1*FNORM(I2,J1))
     >                  - F2    * (R1MF1*FNORM(I1,J2) + F1*FNORM(I2,J2))
            END DO

         END DO

C        Denormalize F:

         FSCALE = ONE / FSCALE
         F(N,I1+1:I2-1,J1+1:J2-1) = FNORM(I1+1:I2-1,J1+1:J2-1) * FSCALE
     >                            + FMIN
      END DO ! Next function

      DEALLOCATE (ARCS, FNORM)

      END SUBROUTINE TFIF2D
