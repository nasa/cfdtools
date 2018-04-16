C+------------------------------------------------------------------------------
C
      SUBROUTINE MORPH_LINE_3D (I1, I2, S1, S2, X0, Y0, Z0, X, Y, Z)
C
C  ONE-LINER: Variant of NULINE3D to control end points AND end-point slopes
C
C  DESCRIPTION:
C
C        MORPH_LINE_3D perturbs the interior points of a 3-space curve given
C     perturbed end points and desired end-point slopes.  More precisely, any
C     pair of points I1, I2 may be the controlling points; each point between
C     these two is moved according to an arc-length-based combination of the
C     distances between the original and new points corresponding to I1 and I2.
C
C        This first stage (that of NULINE3D) is further adjusted to obtain the
C     specified end-point slopes in a similar arc-length-based way which turns
C     out to be overdetermined:  N - 1 intervals in which 2-point derivative
C     estimates are specified but only N - 2 interior points to play with.
C
C        The new coordinates of the two end points should be input in the
C     desired output curve.
C
C  ENVIRONMENT:  Fortran 90
C
C  HISTORY:
C
C     03/16/98  DAS  Original NULINE3D, used here as the first stage.
C     02/18/04   "   MORPH_LINE_3D applies a linear least squares technique
C                    to achieve some degree of slope control at the cost of
C                    losing the original relative spacing.  Therefore ...
C     02/24/04   "   Interpolate back to the original relative spacings, and
C                    iterate to optimize some measure of goodness.
C     02/26/04   "   Replaced dense HDECOM/HSOLVE pair with specialized sparse
C                    pair BIDECOM/BISOLVE, and weighted the end-point slopes.
C                    Excessive weighting can lose slopes in the middle.
C                    Suppress the redistribution of the stage 1 result as
C                    hardly worth the extra arithmetic.
C     06/14/05   "   TOL should be sqrt (eps) for a minimum, not eps.
C
C  AUTHOR:  David Saunders, ELORET/NASA Ames Research Center, Moffet Field, CA.
C
C ------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) ::
     >   I1, I2                ! Indices of given perturbed "end" points
      REAL, INTENT (IN) ::
     >   S1(3), S2(3)          ! Unit tangents desired at points I1, I2
      REAL, INTENT (IN) ::
     >   X0(*), Y0(*), Z0(*)   ! Original curve coordinates
      REAL, INTENT (INOUT) ::
     >   X(*), Y(*), Z(*)      ! Desired curve, with points I1, I2 input

C-------------------------------------------------------------------------------

C     Local constants:

      INTEGER, PARAMETER ::
     >   LUNOUT = -6           ! Suppress FMINRC iterations; diagnostics to +6
      REAL, PARAMETER ::
     >   ONE = 1., ZERO = 0.
      CHARACTER, PARAMETER ::
     >   METHOD * 1 = 'B',     ! Plain Hermite cubics of X, Y, Z vs. arc
     >   CALLER * 8 = 'MORPH_3D'
      LOGICAL, PARAMETER ::
     >   CLOSED = .FALSE.

C     Local variables:

      INTEGER
     >   I, IEVAL, ISTAT, J, LUNERR, M, N, NDATA, NUMFUN
      REAL, DIMENSION (I1:I2) ::
     >   ARC, ARC0, B1, B2, B3, C, D, S, U, W, XI, YI, ZI
      REAL
     >   AV, BV, DX1, DX2, DY1, DY2, DZ1, DZ2, OBJ, T, TOL, TLENGTH1,
     >   V, W1, W2, DERIVS(3)
      LOGICAL
     >   FIRST, NEW

C     Procedures:

      EXTERNAL
     >   BIDECOM,  ! Specialized linear least squares pair
     >   BISOLVE,
     >   FMINRC,   ! Reverse-communication 1-D minimizer
     >   PLSCRV3D  ! Parametric spline interpolation utility

C     Execution:

CCC   write (6, '(/, (a, 3f12.7))') ' s1: ', s1, ' s2: ', s2

      M = I2 - I1;  N = M - 1;  NDATA = M + 1
      DERIVS(1) = -999.  ! Suppresses derivative outputs from interpolations

C     First, apply the simple end-point-location algorithm of NULINE3D:

      ARC0(I1) = ZERO
      DO I = I1 + 1, I2
         ARC0(I) = ARC0(I - 1) + SQRT (
     >     (X0(I) - X0(I - 1)) ** 2 +(Y0(I) - Y0(I - 1)) ** 2 +
     >     (Z0(I) - Z0(I - 1)) ** 2)
      END DO

      T = ONE / ARC0(I2)
      ARC0(I1:I2) = ARC0(I1:I2) * T ! Normalize original arc lengths

      XI(I1) = X(I1)    ! End points of coordinates to be iterated on
      YI(I1) = Y(I1)
      ZI(I1) = Z(I1)
      XI(I2) = X(I2)
      YI(I2) = Y(I2)
      ZI(I2) = Z(I2)

      DX1 = X(I1) - X0(I1)
      DX2 = X(I2) - X0(I2)
      DY1 = Y(I1) - Y0(I1)
      DY2 = Y(I2) - Y0(I2)
      DZ1 = Z(I1) - Z0(I1)
      DZ2 = Z(I2) - Z0(I2)

      DO I = I1 + 1, I2 - 1
         W2 = ARC0(I)
         W1 = ONE - W2
         XI(I) = X0(I) + W1 * DX1 + W2 * DX2
         YI(I) = Y0(I) + W1 * DY1 + W2 * DY2
         ZI(I) = Z0(I) + W1 * DZ1 + W2 * DZ2
      END DO

C     This stage 1 result now needs to be adjusted to obtain the specified
C     end-point slopes.
C     The interim relative spacing is NOT necessarily the same as for the
C     original curve - only if it's a straight line, we believe.
C     Therefore, update the current relative spacing and total length:

      ARC(I1) = ZERO
      DO I = I1 + 1, I2
         ARC(I) = ARC(I - 1) + SQRT (
     >      (XI(I) - XI(I - 1)) ** 2 + (YI(I) - YI(I - 1)) ** 2 +
     >      (ZI(I) - ZI(I - 1)) ** 2)
      END DO

      TLENGTH1 = ARC(I2)
      T = ONE / TLENGTH1

      ARC(I1:I2) = ARC(I1:I2) * T ! Normalized interim arc lengths

CCC   write (6, '(a, f12.7)') ' Stage 1 total length:', tlength1
CCC   write (10, '(a)') ' Stage 1 arc lengths:'
CCC   write (10, '(i3, 2f12.8)') (i, arc0(i), arc(i), i = i1, i2)

C     Impose the original relative spacing precisely on the stage 1 result:
C     NO - this doesn't seem to improve things significantly.

CCC   NEW = .TRUE.
CCC   IEVAL = 1

CCC   DO I = I1 + 1, I2 - 1
CCC      CALL PLSCRV3D (NDATA, XI(I1), YI(I1), ZI(I1), ARC(I1),
CCC  >                  METHOD, NEW, CLOSED, ARC0(I), IEVAL,
CCC  >                  X(I), Y(I), Z(I), DERIVS)
CCC      NEW = .FALSE.
CCC   END DO

CCC   DO I = I1 + 1, I2 - 1 ! Can't do this till PLSCRV3D is done
CCC      XI(I) = X(I)
CCC      YI(I) = Y(I)
CCC      ZI(I) = Z(I)
CCC   END DO
 
C     Copying ARC0(*) here is very close to recalculating arcs:

CCC   DO I = I1 + 1, I2
CCC      ARC(I) = ARC0(I)
CCC   END DO

C     Begin stage 2, where we attempt to impose slopes in ALL intervals,
C     including the first and last interval.  This tends to conflict with
C     imposing the original relative spacing, but we recover that at the
C     end at the expense of giving up a little on the end-point slopes.
C     Two-point derivatives keep the overdetermined system very simple.

      V = ONE ! Reasonable estimate for the unknown multiplier of the
              ! stage 1 total length that gives the final total length
      T = V * TLENGTH1

C     A minimum can be found to within sqrt (epsilon), but avoid the sqrt:

      IF (EPSILON (TOL) < 1.E-10) THEN
         TOL = 1.E-8
      ELSE
         TOL = 1.E-4
      END IF

      FIRST = .TRUE.
      ISTAT = 2        ! Initialize the minimization
      AV = V * 0.5
      BV = V + V
      NUMFUN = 30      ! Limit; FMINRC typically takes about 6 iterations
      LUNERR = ABS (LUNOUT)

   10 CONTINUE

         CALL FMINRC (AV, BV, V, OBJ, TOL, NUMFUN, CALLER, LUNOUT,
     >                ISTAT)

         IF (ISTAT < -1) THEN

            WRITE (LUNERR, '(/, 2A)') CALLER, ': FMINRC fatal error'
            STOP

         ELSE IF (ISTAT < 0) THEN ! Iteration limit; may be usable

            WRITE (LUNERR, '(/, 2A)') CALLER, ': Iteration limit.'

         ELSE IF (ISTAT > 0) THEN ! Evaluate the objective function

            CALL OBJECTIVE
            GO TO 10

         ELSE ! ISTAT = 0 (success).

         END IF

C     Ensure that everything matches V(best), not V(last):

      CALL OBJECTIVE

C     Calculate the optimized arc lengths:

      DO I = I1 + 1, I2
         ARC(I) = ARC(I - 1) + SQRT (
     >      (XI(I) - XI(I - 1)) ** 2 + (YI(I) - YI(I - 1)) ** 2 +
     >      (ZI(I) - ZI(I - 1)) ** 2)
      END DO

CCC   write (6, '(a, f12.7)') ' Current total length:', arc(i2)

      T = ONE / ARC(I2)

      ARC(I1:I2) = T * ARC(I1:I2) ! Renormalize

CCC   write (10, '(a)') ' Penultimate arc lengths:'
CCC   write (10, '(i3, 2f12.8)') (i, arc0(i), arc(i), i = i1, i2)

C     Recover the original relative spacing precisely:

      NEW = .TRUE.
      IEVAL = 1

      DO I = I1 + 1, I2 - 1
         CALL PLSCRV3D (NDATA, XI(I1), YI(I1), ZI(I1), ARC(I1),
     >                  METHOD, NEW, CLOSED, ARC0(I), IEVAL,
     >                  X(I), Y(I), Z(I), DERIVS)
         NEW = .FALSE.
      END DO

CCC   write (6, '(a, 2f12.7)') ' Target slopes: ',
CCC  >   s1(2)/s1(1), s2(2)/s2(1)
CCC   write (6, '(a, 2f12.7)') ' Final  slopes: ',
CCC  >   (y(i1+1) - y(i1)) / (x(i1+1) - x(i1)),
CCC  >   (y(i2) - y(i2-1)) / (x(i2) - x(i2-1))

C     Done.

C     Internal procedure for 2nd stage of 2-stage MORPH_LINE_3D algorithm:

      CONTAINS

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         SUBROUTINE OBJECTIVE

!        Evaluate a function of variable V being minimized to achieve specified
!        [end-point] slopes as well as possible.  V is a multiplier on the stage
!        1 total arc length that produces the final unknown arc length.
!        Actually, the optimal V turns out to differ, and it's not clear why.
!        Including a term aimed at preserving the original relative spacing as
!        well appears to be counter-productive.

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!        Local constants:

         REAL, PARAMETER ::
     >      HALF   = 0.5,
     >      WRANGE = 1.5 ! Highest : lowest weight = 1 : 1 / WRANGE
                         ! with parabolic variation from end to end; 10 too big!
!        Local variables:

         REAL
     >      A, B, DS, RESIDUAL1, RESIDUAL2, RESIDUAL3

!        Execution:

         IF (FIRST) THEN ! Factorize the matrix once only

            FIRST = .FALSE.

            B = ONE / WRANGE   ! Coeffs. of parabola on [0, 1] touching B at 0.5
            A = 4. * (ONE - B)

            DO I = I1, I2 - 1  ! Avoid starting W at I1+1; only the RH sides do
               W2 = (ARC(I) + ARC(I+1)) * HALF  ! Preserve any data symmetry
               W(I) = A * (W2 - HALF) ** 2 + B
            END DO

            CALL BIDECOM (NDATA, W, D, U, C, S) ! Upper bidiagonal + transforms

         END IF

C        The three coordinates present three right-hand-sides:

         T = V * TLENGTH1 ! Using V rather than T as the variable eases
                          ! bracketing the solution for any curve
         DO I = I1 + 1, I2
            W2 = (ARC(I-1) + ARC(I)) * HALF ! Normalized
            W1 = ONE - W2
            DS = (ARC(I) - ARC(I-1)) * T    ! Unnormalized
            B1(I) = DS * (W1 * S1(1) + W2 * S2(1))
            B2(I) = DS * (W1 * S1(2) + W2 * S2(2))
            B3(I) = DS * (W1 * S1(3) + W2 * S2(3))
         END DO

         I = I2
         B1(I) = B1(I) - XI(I)
         B2(I) = B2(I) - YI(I)
         B3(I) = B3(I) - ZI(I)

         I = I1 + 1
         B1(I) = B1(I) + XI(I1)
         B2(I) = B2(I) + YI(I1)
         B3(I) = B3(I) + ZI(I1)
 
         CALL BISOLVE (NDATA, W, D, U, C, S, B1(I), RESIDUAL1)
         CALL BISOLVE (NDATA, W, D, U, C, S, B2(I), RESIDUAL2)
         CALL BISOLVE (NDATA, W, D, U, C, S, B3(I), RESIDUAL3)

         OBJ = RESIDUAL1 + RESIDUAL2 + RESIDUAL3

         DO I = I1 + 1, I2 - 1
            XI(I) = B1(I)
            YI(I) = B2(I)
            ZI(I) = B3(I)
         END DO

CCC      write (6, '(a, 2f12.7)') ' Target   slopes: ',
CCC  >      s1(2)/s1(1), s2(2)/s2(1)
CCC      write (6, '(a, 2f12.7)') ' Iterated slopes: ',
CCC  >      (yi(i1+1) - yi(i1)) / (xi(i1+1) - xi(i1)),
CCC  >      (yi(i2) - yi(i2-1)) / (xi(i2) - xi(i2-1))

CCC      write (6, '(a, f16.13, 1p, 4e19.11)')
CCC  >      ' V/R123/OBJ: ',
CCC  >      v, residual1, residual2, residual3, obj

CCC      write (11, '(a, f16.13, 1p, 4e19.11)')
CCC  >      ' V/R123/OBJ: ',
CCC  >      v, residual1, residual2, residual3, obj

         END SUBROUTINE OBJECTIVE

      END SUBROUTINE MORPH_LINE_3D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine bidecom (n, w, d, u, c, s)
!
!        BIDECOM and BISOLVE treat a specialized weighted linear least squares
!     problem arising from perturbing a grid line in 3-space while controlling
!     the slopes at its new end points (and everywhere in between).  For n grid
!     points, n - 1 equations express the desired 2-point gradients but there
!     are only n - 2 interior grid points to adjust once the end point locations
!     have been imposed.
!
!        If the slopes in the n - 1 interior intervals are weighted by weights
!     w(1), w(2), ..., w(n-1), the relevant (n-1) x (n-2) matrix is:
!
!                       |  w1                              |
!                       | -w2  w2                          |
!                       |     -w3  w3                      |
!                       !         -w4  w4                  |
!                       :              :    :              |
!                       :                   :    :         |
!                       |                        :  w(n-2) |
!                       |                          -w(n-1) |
!
!        Applying a sequence of Givens matrices (symmetric plane rotations)
!     converts this to an upper bidiagonal system defined by vectors d(*) and
!     u(*).  BISOLVE applies the sequence (saved here in vectors c(*) and s(*))
!     to a given right-hand-side then solves the upper triangular system (or
!     rather the upper (n-2) x (n-2) portion of it) to produce the desired
!     linear least squares solution.  The minimized residual (2-norm squared)
!     is the square of the last element of the transformed right-hand side.
!
!     02/25/04  DAS  Initial implementation for MORPH_LINE_3D.
!
!     Author:  David Saunders, ELORET/NASA Ames Research Center, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

!     Arguments:

      integer, intent (in) :: n   ! Underlying number of grid points defining
                                  ! n - 1 intervals and an (n-1) x (n-2) matrix.
      real, intent (in) :: w(n-1) ! Weights scaling the matrix as shown.

      real, intent (out), dimension (n-1) :: d, u, c, s  ! See outline above.

!     Local variables:

      integer i

      real    gi, ci, si, wi

!     Execution:

      d(1) = w(1)   ! Initialize the forward pass

      do i = 2, n - 1

!        Construct the orthogonal matrix that zeros out -w(i), update the
!        diagonal above it and alongside it, and assign the corresponding
!        triangle factor element above the diagonal:

         wi = w(i)
         gi = sqrt (d(i-1)**2 + wi**2);  ! "Gamma"
         ci = d(i-1) / gi;               c(i) = ci
         si =   -wi  / gi;               s(i) = si
         d(i-1) =  gi
         u(i-1) =  si * wi
         d(i)   = -ci * wi

      end do

      end subroutine bidecom

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine bisolve (n, w, d, u, c, s, b, rsq)
!
!        BISOLVE completes the solution of a specialized weighted linear least
!     squares problem for a given right-hand-side vector b.  See BIDECOM for
!     further details.
!
!     02/26/04  DAS  Initial implementation for MORPH_LINE_3D.
!
!     Author:  David Saunders, ELORET/NASA Ames Research Center, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

!     Arguments:

      integer, intent (in) :: n   ! Underlying number of grid points defining
                                  ! n - 1 intervals and an (n-1) x (n-2) matrix
                                  ! which has been triangularized by BIDECOM
      real, intent (in) :: w(n-1) ! Weights for scaling the RHS vector b

      real, intent (in), dimension (n-1) :: d, u, c, s  ! Factorization info.
                                                        ! from BIDECOM
      real, intent (inout) :: b(n-1) ! Input with unweighted RHS;
                                     ! output with least squares solution
      real, intent (out) :: rsq      ! Minimized sum of squares 

!     Local variables:

      integer i

      real    bim1, ci, si

!     Execution:

!     Apply the weights to the RHS:

      do i = 1, n - 1
         b(i) = w(i) * b(i)
      end do

!     Transform the RHS proper:

      do i = 2, n - 1
         bim1   = b(i-1)
         b(i-1) = c(i) * bim1 + s(i) * b(i)
         b(i)   = s(i) * bim1 - c(i) * b(i)
      end do

      rsq = b(n - 1)**2

!     Solve the bidiagonal upper triangular system:

      i = n - 2
      b(i) = b(i) / d(i)

      do i = n - 3, 1, -1
         b(i) = (b(i) - u(i) * b(i+1)) / d(i)
      end do

      end subroutine bisolve
