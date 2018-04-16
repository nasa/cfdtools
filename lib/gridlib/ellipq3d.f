C+------------------------------------------------------------------------------
C
      SUBROUTINE ELLIPQ3D (IDIM, JDIM, I1, I2, J1, J2, X, Y, Z, ARCS,
     >                     PRINT, BGMODE, FGMODE, SPMODE,
     >                     ITMAX, ITFLOAT, ITFREEZE, JLAYER,
     >                     CONVG, CONVMIN, DMAX, OMEGA,
     >                     POWERI, POWERJ, EXPI1, EXPI2, EXPJ1, EXPJ2,
     >                     FGLIMIT, URFG, URFLOAT)
C
C     ELLIPQ3D smooths the interior of the specified 3-D surface (sub)grid
C     with control of spacing and edge orthogonality based on the formulation
C     of GRIDGEN, where P and Q are also referred to as PHI and PSI:
C
C     DEL^2 Xi = P(Xi,Eta) = p1(Eta) exp(-aj1 Xi) + p2(Eta) exp(-aj2(Ximax-Xi))
C                          + p3(Xi) exp(-ai1 Eta) + p4(Xi) exp(-ai2(Etamax-Eta))
C     DEL^2 Eta= Q(Xi,Eta) = q1(Eta) exp(-aj1 Xi) + q2 ...
C
C     The right-hand-side expressions are referred to as "foreground" forcing
C     functions.  Only one of the four components can be significant at a given
C     point, or perhaps two near a corner.  Thomas-Middlecoff-type "background"
C     forcing functions for interior spacing control are (optionally) blended
C     in as the foreground terms are decayed away from the edges.
C
C     Poisson's equation is solved via SLOR.  The equations above are solved
C     in parametric space U and V, where U and V are index- (not arc-) based.
C     The metric coefficients are determined in X/Y/Z space, and lagged by one
C     iteration.  See page 5-22 of GRIDGEN Volume I: The Final Report.
C
C     The initial spacing off a given edge may either be derived from the end-
C     point increments at adjoining edges (only), or it may preserve the input
C     first increments along the entire edge.
C
C     References:  Steinbrenner et. all, GRIDGEN final report ( July 1990)
C                  Sorenson (May 1980) NASA Tech. Mem. 81198
C                  Thompson, Warsi, Mastin (1985), p. 231
C
C     NOTE:  The initial surface grid must represent the underlying surface
C            adequately, since the smoothed grid is obtained by parametric
C            bilinear interpolation (evaluation) within the original grid.
C
C     Typical inputs for half of an Euler-type C-grid:
C
C     BGMODE = 'YYYY', FGMODE = 'YYYY', SPMODE = 'YYYY',
C     ITMAX = 200, ITFLOAT = 0, ITFREEZE = 100, JLAYER = 4,
C     CONVG = 5., CONVMIN = 0., DMAX = 0.001, OMEGA = 1.2,
C     POWERI = 1., POWERJ = .5, EXPI1 = .4, EXPI2 = .4, EXPJ1 = .4, EXPJ2 = .5,
C     FGLIMIT = 1., URFG = .1, URFLOAT = .1
C
C     ???????? Lockheed   Original WBGRID implementation (TTM2D).
C     09/13/94 D.Saunders Replaced TRIB with TRID2R.
C     10/03/94   "   "    Introduced IDIM, KDIM, I1, K1;
C                         Dubious PHI1(ILE) = 0. suppressed.
C     02/07/95   "   "    Avoid tests for coincident boundary pts;
C                         replace IFs with MAX (MIN ( to vectorize.
C     05/12/95   "   "    K, Z replaced by J, Y to match TTM3D changes.
C     05/21/95   "   "    Arg. ILE /= 999 forces PHI = 0. at I = ILE.
C                         [ILE has since been eliminated, but it may be
C                         helpful on single-block C grids.]
C                         2D source limiters are now separate from 3D's;
C                         application at the boundaries suffices.
C     6/95  DAS/J.Reuther Control of edge orthogonality introduced from
C                         Thompson, et al. Right-handedness is assumed.
C     06/23/95    DAS     Arguments BGMODE, FGMODE simplify application
C                         to various situations - preferred to making
C                         all the control inputs arguments.  E.g., a
C                         degenerate C grid must have FRCPJ1 = 1., but
C                         this is poor for regular airfoil C grids.
C                         BGMODE = 'YYYY' activates edges I1, I2, J1, J2;
C                         FGMODE = 'NNNN' disables orthogonality control;
C                         and so on.
C     7/95-8/95 S.Edwards Unsatisfactory orthogonality control replaced
C                         by a Sorenson-type implementation, along with
C                         an option to float the outer J2 boundary.
C     9/95     D.Saunders TTM2D renamed as ELLIP2D; some clean-up; fully
C                         argument-driven apart from much local work-space.
C     7/96      SJE/DAS   Some glitches fixed; arc-length-based interpolation
C                         of corner increments into edge interiors; nonlinear
C                         arc-length-based interpolation of edge PHIs & PSIs
C                         into the interior; in-line tridiagonal solutions.
C     8/96       SJE      Adapted ELLIP2D to handle 3-D surfaces using index-
C                         based parametric UV space.
C     6/97      SJE/DAS   Off-diagonal terms subtracted in main loop;
C                         background edge PHI and PSI values smoothed,
C                         requiring corner values to be generated;
C                         removed floating of J2 boundary.
C     8/97       "   "    POWERI < 0. extension needed for symmetric C grids.
C     4/27/98     DAS     Fortran 90 allows use of automatic arrays.
C ------------------------------------------------------------------------------

C     Arguments:

      INTEGER    IDIM, JDIM       ! I   Calling program dimensions for X, Y, Z
      INTEGER    I1, I2, J1, J2   ! I   Active element range for X, Y, Z;
                                  !     I2, J2 <= IDIM, JDIM
      REAL       X(IDIM,JDIM),    ! I/O Grid coordinates, input with initial
     >           Y(IDIM,JDIM),    !     estimates
     >           Z(IDIM,JDIM)     !
                                  !
      REAL     ARCS(IDIM,JDIM,1,3)! S   Work-space for relative arc lengths;
                                  !     dimensions must match X, Y, Z in
                                  !     the call to PARAMXYZ
      LOGICAL    PRINT            ! I   .TRUE. shows iterations on unit 6
                                  !
      CHARACTER  BGMODE * 4       ! I   4 upper case characters referring to
                                  !     edges I1, I2, J1, J2 in that order;
                                  !     'Y' activates "background" control
                                  !     of interior spacing from that edge.
      CHARACTER  FGMODE * 4       ! I   'Y' activates "foreground" control
                                  !     of orthogonality at that edge.
      CHARACTER  SPMODE * 4       ! I   'Y' means interior initial spacing
                                  !     increments are controlled from the
                                  !     corner increments - possibly non-
                                  !     linearly - see POWERI, POWERJ;
                                  !     'N' means starting guess values are
                                  !     retained for the entire boundary
      INTEGER    ITMAX            ! I   Limit on no. of smoothing iterations
      INTEGER    ITFLOAT          ! I   Retained from ELLIP2D but inactive
      INTEGER    ITFREEZE         ! I   "  "
      INTEGER    JLAYER           ! I   "  "
      REAL       CONVG            ! I   Orders of magnitude reduction sought
                                  !     in the maximum change in dU or dV
      REAL       CONVMIN          ! I   CONVMIN > 0. suppresses foreground
                                  !     terms until some degree of convergence
                                  !     is achieved
      REAL       DMAX             ! I   Allows early termination if dU and dV
                                  !     are all below this value
      REAL       OMEGA            ! I   SLOR over-relaxation parameter
      REAL       POWERI, POWERJ   ! I   1. means linear combination of opposite
                                  !        edge PHI/PSI for background control;
                                  !   > 1. emphasizes the lower edge;
                                  !   < 1.     "       "  upper edge;
                                  !     use POWERI/J < 0. to switch meaning of
                                  !     "lower" and "upper" - e.g., POWERI < 0.
                                  !     for the upper half of a typical C grid
      REAL       EXPI1, EXPI2,    ! I   Decay factors for foreground terms;
     >           EXPJ1, EXPJ2     !     larger means more rapid decay
      REAL       FGLIMIT          ! I   Limiters on the foreground terms;
                                  !     FGLIMIT / (1+edge curvature) is used
      REAL       URFG             ! I   Under-relaxation factors applied to
                                  !     the changes in the foreground terms
      REAL       URFLOAT          ! I   Under-relaxes the boundary floating

C-------------------------------------------------------------------------------

C     Local constants:

      REAL, PARAMETER ::
     >   EPS = 1.E-16, ! Safeguards 0/0
     >   FIVE = 5., FOUR = 4., HALF = 0.5, ONE = 1., TWO = 2., ZERO = 0.

      CHARACTER * 1, PARAMETER ::
     >   YES = 'Y', NO = 'N'

C     Automatic arrays:

      REAL, DIMENSION (I1:I2) ::
     >   A, B, C, DECAYI1, DECAYI2, FGMAXJ1, FGMAXJ2, P1, P2, Q1, Q2,
     >   PHIJ1, PHIJ2, RHS1, RHS2, RHS3, SEDGEJ1, SEDGEJ2, SJ1, SJ2,
     >   XXJ1, YXJ1, ZXJ1, XXXJ1, YXXJ1, ZXXJ1,
     >   XEJ1, YEJ1, ZEJ1, XEXJ1, YEXJ1, ZEXJ1,
     >   XXJ2, YXJ2, ZXJ2, XXXJ2, YXXJ2, ZXXJ2,
     >   XEJ2, YEJ2, ZEJ2, XEXJ2, YEXJ2, ZEXJ2

      REAL, DIMENSION (J1:J2) ::
     >   DECAYJ1, DECAYJ2, FGMAXI1, FGMAXI2, P3, P4, Q3, Q4,
     >   PSII1, PSII2, SEDGEI1, SEDGEI2, SI1, SI2,
     >   XEI1, YEI1, ZEI1, XEEI1, YEEI1, ZEEI1,
     >   XXI1, YXI1, ZXI1, XEXI1, YEXI1, ZEXI1,
     >   XEI2, YEI2, ZEI2, XEEI2, YEEI2, ZEEI2,
     >   XXI2, YXI2, ZXI2, XEXI2, YEXI2, ZEXI2

      REAL, DIMENSION (I1:I2,J1:J2) ::
     >   AJAC, ALP, BET, GAM, U, V, XORIG, YORIG, ZORIG

C     Execution:

      IF (ITMAX == 0) GO TO 999

C     Activate the requested "background" spacing control functions.
C     ("FRC" originally distinguished 2D from 3D "SRC" source terms.)

      FRCSI1 = ONE
      IF (BGMODE(1:1) /= YES) FRCSI1 = ZERO ! S = "PSI"
      FRCSI2 = ONE
      IF (BGMODE(2:2) /= YES) FRCSI2 = ZERO
      FRCPJ1 = ONE
      IF (BGMODE(3:3) /= YES) FRCPJ1 = ZERO ! P = "PHI"
      FRCPJ2 = ONE
      IF (BGMODE(4:4) /= YES) FRCPJ2 = ZERO

      I1P1 = I1 + 1
      I2M1 = I2 - 1
      J1P1 = J1 + 1
      J2M1 = J2 - 1

      RI2MI1 = ONE / REAL (I2 - I1)
      RJ2MJ1 = ONE / REAL (J2 - J1)

C     If the edges need independent limiting and under-relaxing, make
C     all of the following arguments (as they were originally):

      FGLIMI1 = FGLIMIT
      FGLIMI2 = FGLIMIT
      FGLIMJ1 = FGLIMIT
      FGLIMJ2 = FGLIMIT

      URI1 = URFG
      URI2 = URFG
      URJ1 = URFG
      URJ2 = URFG

C     Starting-guess normalized arc lengths for this sub-grid:

      CALL PARAMXYZ (1, IDIM, 1, JDIM, 1, 1,
     >               I1, I2,  J1, J2,  1, 1, X, Y, Z, ARCS)

C     Edge arc lengths are needed for explicit smoothing of background fns.

      SEDGEJ1 = ARCS(I1:I2,J1,1,1)
      SEDGEJ2 = ARCS(I1:I2,J2,1,1)

      SEDGEI1 = ARCS(I1,J1:J2,1,2)
      SEDGEI2 = ARCS(I2,J1:J2,1,2)

C     Heuristic foreground decay terms:

      I1POWER = NINT (-LOG10 (MIN (ARCS(I1P1,J1,1,1),
     >                             ARCS(I1P1,J2,1,1))))
      I1POWER = MAX ((I1POWER + 1)/ 2, 3)
      I2POWER = NINT (-LOG10 (ONE - MAX (ARCS(I2M1,J1,1,1),
     >                                   ARCS(I2M1,J2,1,1))))
      I2POWER = MAX ((I2POWER + 1)/ 2, 3)

      J1POWER = NINT (-LOG10 (MIN (ARCS(I1,J1P1,1,2),
     >                             ARCS(I2,J1P1,1,2))))
      J1POWER = MAX ((J1POWER + 1)/ 2, 3)
      J2POWER = NINT (-LOG10 (ONE - MAX (ARCS(I1,J2M1,1,2),
     >                                   ARCS(I2,J2M1,1,2))))
      J2POWER = MAX ((J2POWER + 1)/ 2, 3)

C     Store the exponential decay coefficients to avoid repeated EXPs:

      IF (FGMODE(1:1) == NO) THEN ! No foreground

         FRCFI1  = ZERO
         DECAYI1 = ZERO

      ELSE ! (FGMODE(1:1) == YES) ! Index-based decay

         FRCFI1 = ONE
         DO I = I1, I2
            DECAYI1(I) = (RI2MI1 * REAL (I2 - I)) ** I1POWER
     >             * EXP (-EXPI1 * REAL (I - I1))
         END DO

      END IF

      IF (FGMODE(2:2) == NO) THEN

         FRCFI2  = ZERO
         DECAYI2 = ZERO

      ELSE

         FRCFI2 = ONE
         DO I = I1, I2
            DECAYI2(I) = (RI2MI1 * REAL (I - I1)) ** I2POWER
     >             * EXP (-EXPI2 * REAL (I2 - I))
         END DO

      END IF

      IF (FGMODE(3:3) == NO) THEN

         FRCFJ1  = ZERO
         DECAYJ1 = ZERO

      ELSE

         FRCFJ1 = ONE
         DO J = J1, J2
            DECAYJ1(J) = (RJ2MJ1 * REAL (J2 - J)) ** J1POWER
     >             * EXP (-EXPJ1 * REAL (J - J1))
         END DO

      END IF

      IF (FGMODE(4:4) == NO) THEN

         FRCFJ2  = ZERO
         DECAYJ2 = ZERO

      ELSE

         FRCFJ2 = ONE
         DO J = J1, J2
            DECAYJ2(J) = (RJ2MJ1 * REAL (J - J1)) ** J2POWER
     >             * EXP (-EXPJ2 * REAL (J2 - J))
         END DO

      END IF

C     Moderating the INTERIOR relative arc-lengths (POWER* < 1.)
C     helps background control if the curvature at one boundary
C     greatly exceeds the curvature at the opposite boundary.

      IF (ABS (POWERI) /= ONE) THEN
         IF (POWERI == HALF) THEN
            DO J = J1P1, J2M1
               DO I = I1P1, I2M1
                  ARCS(I,J,1,1) = SQRT (ARCS(I,J,1,1))
               END DO
            END DO
         ELSE IF (POWERI == -HALF) THEN
            DO J = J1P1, J2M1
               DO I = I1P1, I2M1
                  ARCS(I,J,1,1) = ONE - SQRT (ONE - ARCS(I,J,1,1))
               END DO
            END DO
         ELSE IF (POWERI > ZERO) THEN
            DO J = J1P1, J2M1
               DO I = I1P1, I2M1
                  ARCS(I,J,1,1) = ARCS(I,J,1,1) ** POWERI
               END DO
            END DO
         ELSE ! POWERI < ZERO
            DO J = J1P1, J2M1
               DO I = I1P1, I2M1
                  ARCS(I,J,1,1) = ONE - (ONE-ARCS(I,J,1,1)) ** (-POWERI)
               END DO
            END DO
         END IF
      END IF

      IF (ABS (POWERJ) /= ONE) THEN
         IF (POWERJ == HALF) THEN
            DO J = J1P1, J2M1
               DO I = I1P1, I2M1
                  ARCS(I,J,1,2) = SQRT (ARCS(I,J,1,2))
               END DO
            END DO
         ELSE IF (POWERJ == -HALF) THEN
            DO J = J1P1, J2M1
               DO I = I1P1, I2M1
                  ARCS(I,J,1,2) = ONE - SQRT (ONE - ARCS(I,J,1,2))
               END DO
            END DO
         ELSE IF (POWERJ > ZERO) THEN
            DO J = J1P1, J2M1
               DO I = I1P1, I2M1
                  ARCS(I,J,1,2) = ARCS(I,J,1,2) ** POWERJ
               END DO
            END DO
         ELSE ! POWERJ < ZERO
            DO J = J1P1, J2M1
               DO I = I1P1, I2M1
                  ARCS(I,J,1,2) = ONE - (ONE-ARCS(I,J,1,2)) ** (-POWERJ)
               END DO
            END DO
         END IF
      END IF

C     D1MN = distance increment from corner (M,N) along edge M:

      D1I1J1 = SQRT ((X(I1,J1P1) - X(I1,J1))**2 +
     >               (Y(I1,J1P1) - Y(I1,J1))**2 +
     >               (Z(I1,J1P1) - Z(I1,J1))**2)
      D1I2J1 = SQRT ((X(I2,J1P1) - X(I2,J1))**2 +
     >               (Y(I2,J1P1) - Y(I2,J1))**2 +
     >               (Z(I2,J1P1) - Z(I2,J1))**2)
      D1I1J2 = SQRT ((X(I1,J2M1) - X(I1,J2))**2 +
     >               (Y(I1,J2M1) - Y(I1,J2))**2 +
     >               (Z(I1,J2M1) - Z(I1,J2))**2)
      D1I2J2 = SQRT ((X(I2,J2M1) - X(I2,J2))**2 +
     >               (Y(I2,J2M1) - Y(I2,J2))**2 +
     >               (Z(I2,J2M1) - Z(I2,J2))**2)

      D1J1I1 = SQRT ((X(I1P1,J1) - X(I1,J1))**2 +
     >               (Y(I1P1,J1) - Y(I1,J1))**2 +
     >               (Z(I1P1,J1) - Z(I1,J1))**2)
      D1J2I1 = SQRT ((X(I1P1,J2) - X(I1,J2))**2 +
     >               (Y(I1P1,J2) - Y(I1,J2))**2 +
     >               (Z(I1P1,J2) - Z(I1,J2))**2)
      D1J1I2 = SQRT ((X(I2M1,J1) - X(I2,J1))**2 +
     >               (Y(I2M1,J1) - Y(I2,J1))**2 +
     >               (Z(I2M1,J1) - Z(I2,J1))**2)
      D1J2I2 = SQRT ((X(I2M1,J2) - X(I2,J2))**2 +
     >               (Y(I2M1,J2) - Y(I2,J2))**2 +
     >               (Z(I2M1,J2) - Z(I2,J2))**2)

C     Is off-boundary spacing controlled by end-pt. edge increments only (Y)
C     or preserved from the starting guess for the entire boundary?

      IF (SPMODE(1:1) == YES) THEN
         DO J = J1, J2
            RJ = ARCS(I1,J,1,2)
            SI1(J) = (ONE - RJ) * D1J1I1 + RJ * D1J2I1
         END DO
      ELSE
         DO J = J1, J2
            SI1(J) = SQRT ((X(I1P1,J) - X(I1,J))**2 +
     >                     (Y(I1P1,J) - Y(I1,J))**2 +
     >                     (Z(I1P1,J) - Z(I1,J))**2)
         END DO
      END IF

      IF (SPMODE(2:2) == YES) THEN
         DO J = J1, J2
            RJ = ARCS(I2,J,1,2)
            SI2(J) = (ONE - RJ) * D1J1I2 + RJ * D1J2I2
         END DO
      ELSE
         DO J = J1, J2
            SI2(J) = SQRT ((X(I2,J) - X(I2M1,J))**2 +
     >                     (Y(I2,J) - Y(I2M1,J))**2 +
     >                     (Z(I2,J) - Z(I2M1,J))**2)
         END DO
      END IF

      IF (SPMODE(3:3) == YES) THEN
         DO I = I1, I2
            RI = ARCS(I,J1,1,1)
            SJ1(I) = (ONE - RI) * D1I1J1 + RI * D1I2J1
         END DO
      ELSE
         DO I = I1, I2
            SJ1(I) = SQRT ((X(I,J1P1) - X(I,J1))**2 +
     >                     (Y(I,J1P1) - Y(I,J1))**2 +
     >                     (Z(I,J1P1) - Z(I,J1))**2)
         END DO
      END IF

      IF (SPMODE(4:4) == YES) THEN
         DO I = I1, I2
            RI = ARCS(I,J2,1,1)
            SJ2(I) = (ONE - RI) * D1I1J2 + RI * D1I2J2
         END DO
      ELSE
         DO I = I1, I2
            SJ2(I) = SQRT ((X(I,J2) - X(I,J2M1))**2 +
     >                     (Y(I,J2) - Y(I,J2M1))**2 +
     >                     (Z(I,J2) - Z(I,J2M1))**2)
         END DO
      END IF

C     Initialize U and V to the grid indices:

      DO J = J1, J2
         DO I = I1, I2
            U(I,J) = REAL (I)
            V(I,J) = REAL (J)
            XORIG(I,J) = X(I,J)
            YORIG(I,J) = Y(I,J)
            ZORIG(I,J) = Z(I,J)
         END DO
      END DO

C     Set up the J1/J2 edge values of the background contributions to PHI,
C     and the foreground terms involving edge-wise derivatives.
C     The corner PHI and PSI values are needed only for explicit smoothing
C     of the full interior edge distributions.

      I = I1
      J = J1

      XXJ1(I)  = -1.5*X(I,J) + TWO*X(I1P1,J) - HALF*X(I1P1+1,J)
      YXJ1(I)  = -1.5*Y(I,J) + TWO*Y(I1P1,J) - HALF*Y(I1P1+1,J)
      ZXJ1(I)  = -1.5*Z(I,J) + TWO*Z(I1P1,J) - HALF*Z(I1P1+1,J)

      GAM(I,J) = XXJ1(I)**2 + YXJ1(I)**2 + ZXJ1(I)**2

C     The 3-point formula seems to lead to discontinuous PHI at the corner:

      XXXJ1(I) = TWO*X(I,J)-FIVE*X(I1P1,J)+FOUR*X(I1P1+1,J) -X(I1P1+2,J)
      YXXJ1(I) = TWO*Y(I,J)-FIVE*Y(I1P1,J)+FOUR*Y(I1P1+1,J) -Y(I1P1+2,J)
      ZXXJ1(I) = TWO*Z(I,J)-FIVE*Z(I1P1,J)+FOUR*Z(I1P1+1,J) -Z(I1P1+2,J)

      PHI = -(XXJ1(I) * XXXJ1(I) + YXJ1(I) * YXXJ1(I) +
     >        ZXJ1(I) * ZXXJ1(I)) / (GAM(I,J) + EPS)

      PHIJ1(I) = MIN (MAX (PHI, -FRCPJ1), FRCPJ1)

      XEJ1(I)  = -1.5*X(I,J) + TWO*X(I,J1P1) - HALF*X(I,J1P1+1)
      YEJ1(I)  = -1.5*Y(I,J) + TWO*Y(I,J1P1) - HALF*Y(I,J1P1+1)
      ZEJ1(I)  = -1.5*Z(I,J) + TWO*Z(I,J1P1) - HALF*Z(I,J1P1+1)

      J = J2

      XXJ2(I)  = -1.5*X(I,J) + TWO*X(I1P1,J) - HALF*X(I1P1+1,J)
      YXJ2(I)  = -1.5*Y(I,J) + TWO*Y(I1P1,J) - HALF*Y(I1P1+1,J)
      ZXJ2(I)  = -1.5*Z(I,J) + TWO*Z(I1P1,J) - HALF*Z(I1P1+1,J)

      GAM(I,J) = XXJ2(I)**2 + YXJ2(I)**2 + ZXJ2(I)**2

      XXXJ2(I) = TWO*X(I,J)-FIVE*X(I1P1,J)+FOUR*X(I1P1+1,J) -X(I1P1+2,J)
      YXXJ2(I) = TWO*Y(I,J)-FIVE*Y(I1P1,J)+FOUR*Y(I1P1+1,J) -Y(I1P1+2,J)
      ZXXJ2(I) = TWO*Z(I,J)-FIVE*Z(I1P1,J)+FOUR*Z(I1P1+1,J) -Z(I1P1+2,J)

      PHI = -(XXJ2(I) * XXXJ2(I) + YXJ2(I) * YXXJ2(I) +
     >        ZXJ2(I) * ZXXJ2(I)) / (GAM(I,J) + EPS)

      PHIJ2(I) = MIN (MAX (PHI, -FRCPJ2), FRCPJ2)

      XEJ2(I)  = 1.5*X(I,J) - TWO*X(I,J2M1) + HALF*X(I,J2M1-1)
      YEJ2(I)  = 1.5*Y(I,J) - TWO*Y(I,J2M1) + HALF*Y(I,J2M1-1)
      ZEJ2(I)  = 1.5*Z(I,J) - TWO*Z(I,J2M1) + HALF*Z(I,J2M1-1)

      DO I = I1P1, I2M1
         XXJ1(I) = (X(I+1,J1) - X(I-1,J1)) * HALF
         YXJ1(I) = (Y(I+1,J1) - Y(I-1,J1)) * HALF
         ZXJ1(I) = (Z(I+1,J1) - Z(I-1,J1)) * HALF

         XN = -1.5*X(I,J1) + TWO*X(I,J1P1) - HALF*X(I,J1P1+1)
         YN = -1.5*Y(I,J1) + TWO*Y(I,J1P1) - HALF*Y(I,J1P1+1)
         ZN = -1.5*Z(I,J1) + TWO*Z(I,J1P1) - HALF*Z(I,J1P1+1)

         GAM13 = YXJ1(I) * ZN - ZXJ1(I) * YN
         GAM23 = XN * ZXJ1(I) - ZN * XXJ1(I)
         GAM33 = XXJ1(I) * YN - YXJ1(I) * XN

C        Surface normals:

         RALP33 = ONE / SQRT (GAM13**2 + GAM23**2 + GAM33**2 + EPS)
         ORTHO1 = GAM13*RALP33
         ORTHO2 = GAM23*RALP33
         ORTHO3 = GAM33*RALP33

         GAM12 = ZXJ1(I)*ORTHO2 - YXJ1(I)*ORTHO3
         GAM22 = XXJ1(I)*ORTHO3 - ZXJ1(I)*ORTHO1
         GAM32 = YXJ1(I)*ORTHO1 - XXJ1(I)*ORTHO2

         DQ = ONE / SQRT (GAM12**2 + GAM22**2 + GAM32**2 + EPS)
         XEJ1(I) = (DQ * SJ1(I)) * GAM12
         YEJ1(I) = (DQ * SJ1(I)) * GAM22
         ZEJ1(I) = (DQ * SJ1(I)) * GAM32

         ALP(I,J1) = XEJ1(I)**2 + YEJ1(I)**2 + ZEJ1(I)**2
         BET(I,J1) = XEJ1(I)*XXJ1(I) + YEJ1(I)*YXJ1(I) +
     >               ZEJ1(I)*ZXJ1(I)
         GAM(I,J1) = XXJ1(I)**2 + YXJ1(I)**2 + ZXJ1(I)**2

         AJAC(I,J1) = SQRT (ALP(I,J1)*GAM(I,J1) -
     >                      BET(I,J1)*BET(I,J1)) + EPS

         XXXJ1(I) = X(I+1,J1) - TWO*X(I,J1) + X(I-1,J1)
         YXXJ1(I) = Y(I+1,J1) - TWO*Y(I,J1) + Y(I-1,J1)
         ZXXJ1(I) = Z(I+1,J1) - TWO*Z(I,J1) + Z(I-1,J1)

         CURV = SQRT (((XXJ1(I)*YXXJ1(I) - XXXJ1(I)*YXJ1(I))**2 +
     >                 (YXJ1(I)*ZXXJ1(I) - YXXJ1(I)*ZXJ1(I))**2 +
     >                 (ZXJ1(I)*XXXJ1(I) - ZXXJ1(I)*XXJ1(I))**2) /
     >                (GAM(I,J1) + EPS))
         FGMAXJ1(I) = MIN (ONE, FGLIMJ1 / (ONE + CURV))

         PHI = -(XXJ1(I) * XXXJ1(I) + YXJ1(I) * YXXJ1(I) +
     >           ZXJ1(I) * ZXXJ1(I)) / (GAM(I,J1) + EPS)

         PHIJ1(I) = MIN (MAX (PHI, -FRCPJ1), FRCPJ1) ! Background PHI at J1 edge

C        And for the J2 edge ...

         XXJ2(I) = (X(I+1,J2) - X(I-1,J2)) * HALF
         YXJ2(I) = (Y(I+1,J2) - Y(I-1,J2)) * HALF
         ZXJ2(I) = (Z(I+1,J2) - Z(I-1,J2)) * HALF

         XN = 1.5*X(I,J2) - TWO*X(I,J2M1) + HALF*X(I,J2M1-1)
         YN = 1.5*Y(I,J2) - TWO*Y(I,J2M1) + HALF*Y(I,J2M1-1)
         ZN = 1.5*Z(I,J2) - TWO*Z(I,J2M1) + HALF*Z(I,J2M1-1)

         GAM13 = YXJ2(I) * ZN - ZXJ2(I) * YN
         GAM23 = XN * ZXJ2(I) - ZN * XXJ2(I)
         GAM33 = XXJ2(I) * YN - YXJ2(I) * XN

         RALP33 = ONE / SQRT (GAM13**2 + GAM23**2 + GAM33*2 + EPS)
         ORTHO1 = GAM13*RALP33
         ORTHO2 = GAM23*RALP33
         ORTHO3 = GAM33*RALP33

         GAM12 = ZXJ2(I)*ORTHO2 - YXJ2(I)*ORTHO3
         GAM22 = XXJ2(I)*ORTHO3 - ZXJ2(I)*ORTHO1
         GAM32 = YXJ2(I)*ORTHO1 - XXJ2(I)*ORTHO2

         DQ = ONE / SQRT (GAM12**2 + GAM22**2 + GAM32**2 + EPS)

         XEJ2(I) = (DQ * SJ2(I)) * GAM12
         YEJ2(I) = (DQ * SJ2(I)) * GAM22
         ZEJ2(I) = (DQ * SJ2(I)) * GAM32

         ALP(I,J2) = XEJ2(I)**2 + YEJ2(I)**2 + ZEJ2(I)**2
         BET(I,J2) = XEJ2(I)*XXJ2(I) + YEJ2(I)*YXJ2(I) +
     >               ZEJ2(I)*ZXJ2(I)
         GAM(I,J2) = XXJ2(I)**2 + YXJ2(I)**2 + ZXJ2(I)**2

         AJAC(I,J2) = SQRT (ALP(I,J2)*GAM(I,J2) -
     >                      BET(I,J2)*BET(I,J2)) + EPS

         XXXJ2(I) = X(I+1,J2) - TWO*X(I,J2) + X(I-1,J2)
         YXXJ2(I) = Y(I+1,J2) - TWO*Y(I,J2) + Y(I-1,J2)
         ZXXJ2(I) = Z(I+1,J2) - TWO*Z(I,J2) + Z(I-1,J2)

         CURV = SQRT (((XXJ2(I)*YXXJ2(I) - XXXJ2(I)*YXJ2(I))**2 +
     >                 (YXJ2(I)*ZXXJ2(I) - YXXJ2(I)*ZXJ2(I))**2 +
     >                 (ZXJ2(I)*XXXJ2(I) - ZXXJ2(I)*XXJ2(I))**2) /
     >                (GAM(I,J2) + EPS))
         FGMAXJ2(I) = MIN (ONE, FGLIMJ2 / (ONE + CURV))

         PHI = -(XXJ2(I) * XXXJ2(I) + YXJ2(I) * YXXJ2(I) +
     >           ZXJ2(I) * ZXXJ2(I)) / (GAM(I,J2) + EPS)

         PHIJ2(I) = MIN (MAX (PHI, -FRCPJ2), FRCPJ2) ! Background PHI at J2 edge
      END DO

      I = I2
      J = J1

      XXJ1(I)  = 1.5*X(I,J) - TWO*X(I2M1,J) + HALF*X(I2M1-1,J)
      YXJ1(I)  = 1.5*Y(I,J) - TWO*Y(I2M1,J) + HALF*Y(I2M1-1,J)
      ZXJ1(I)  = 1.5*Z(I,J) - TWO*Z(I2M1,J) + HALF*Z(I2M1-1,J)

      GAM(I,J) = XXJ1(I)**2 + YXJ1(I)**2 + ZXJ1(I)**2

      XXXJ1(I) = TWO*X(I,J)-FIVE*X(I2M1,J)+FOUR*X(I2M1-1,J) -X(I2M1-2,J)
      YXXJ1(I) = TWO*Y(I,J)-FIVE*Y(I2M1,J)+FOUR*Y(I2M1-1,J) -Y(I2M1-2,J)
      ZXXJ1(I) = TWO*Z(I,J)-FIVE*Z(I2M1,J)+FOUR*Z(I2M1-1,J) -Z(I2M1-2,J)

      PHI = -(XXJ1(I) * XXXJ1(I) + YXJ1(I) * YXXJ1(I) +
     >        ZXJ1(I) * ZXXJ1(I)) / (GAM(I,J) + EPS)

      PHIJ1(I) = MIN (MAX (PHI, -FRCPJ1), FRCPJ1)

      XEJ1(I)  = -1.5*X(I,J) + TWO*X(I,J1P1) - HALF*X(I,J1P1+1)
      YEJ1(I)  = -1.5*Y(I,J) + TWO*Y(I,J1P1) - HALF*Y(I,J1P1+1)
      ZEJ1(I)  = -1.5*Z(I,J) + TWO*Z(I,J1P1) - HALF*Z(I,J1P1+1)

      J = J2

      XXJ2(I)  =  1.5*X(I,J) - TWO*X(I2M1,J) + HALF*X(I2M1-1,J)
      YXJ2(I)  =  1.5*Y(I,J) - TWO*Y(I2M1,J) + HALF*Y(I2M1-1,J)
      ZXJ2(I)  =  1.5*Z(I,J) - TWO*Z(I2M1,J) + HALF*Z(I2M1-1,J)

      GAM(I,J) = XXJ2(I)**2 + YXJ2(I)**2 + ZXJ2(I)**2

      XXXJ2(I) = TWO*X(I,J)-FIVE*X(I2M1,J)+FOUR*X(I2M1-1,J) -X(I2M1-2,J)
      YXXJ2(I) = TWO*Y(I,J)-FIVE*Y(I2M1,J)+FOUR*Y(I2M1-1,J) -Y(I2M1-2,J)
      ZXXJ2(I) = TWO*Z(I,J)-FIVE*Z(I2M1,J)+FOUR*Z(I2M1-1,J) -Z(I2M1-2,J)

      PHI = -(XXJ2(I) * XXXJ2(I) + YXJ2(I) * YXXJ2(I) +
     >        ZXJ2(I) * ZXXJ2(I)) / (GAM(I,J) + EPS)

      PHIJ2(I) = MIN (MAX (PHI, -FRCPJ2), FRCPJ2)

      XEJ2(I)  = 1.5*X(I,J) - TWO*X(I,J2M1) + HALF*X(I,J2M1-1)
      YEJ2(I)  = 1.5*Y(I,J) - TWO*Y(I,J2M1) + HALF*Y(I,J2M1-1)
      ZEJ2(I)  = 1.5*Z(I,J) - TWO*Z(I,J2M1) + HALF*Z(I,J2M1-1)

C     Cross-derivatives:

      DO I = I1P1, I2M1
         XEXJ1(I) = (XEJ1(I+1) - XEJ1(I-1)) * HALF
         YEXJ1(I) = (YEJ1(I+1) - YEJ1(I-1)) * HALF
         ZEXJ1(I) = (ZEJ1(I+1) - ZEJ1(I-1)) * HALF
         XEXJ2(I) = (XEJ2(I+1) - XEJ2(I-1)) * HALF
         YEXJ2(I) = (YEJ2(I+1) - YEJ2(I-1)) * HALF
         ZEXJ2(I) = (ZEJ2(I+1) - ZEJ2(I-1)) * HALF
      END DO

C     Likewise for the I1 and I2 edges ...

      I = I1
      J = J1

      XEI1(J) = -1.5*X(I,J) + TWO*X(I,J1P1) - HALF*X(I,J1P1+1)
      YEI1(J) = -1.5*Y(I,J) + TWO*Y(I,J1P1) - HALF*Y(I,J1P1+1)
      ZEI1(J) = -1.5*Z(I,J) + TWO*Z(I,J1P1) - HALF*Z(I,J1P1+1)

      ALP(I,J) = XEI1(J)**2 + YEI1(J)**2 + ZEI1(J)**2

      XEEI1(J) = TWO*X(I,J)-FIVE*X(I,J1P1)+FOUR*X(I,J1P1+1) -X(I,J1P1+2)
      YEEI1(J) = TWO*Y(I,J)-FIVE*Y(I,J1P1)+FOUR*Y(I,J1P1+1) -Y(I,J1P1+2)
      ZEEI1(J) = TWO*Z(I,J)-FIVE*Z(I,J1P1)+FOUR*Z(I,J1P1+1) -Z(I,J1P1+2)

      PSI = -(XEI1(J)*XEEI1(J) + YEI1(J)*YEEI1(J) +
     >        ZEI1(J)*ZEEI1(J)) / (ALP(I,J) + EPS)
      PSII1(J) = MIN (MAX (PSI, -FRCSI1), FRCSI1)

      XXI1(J) = -1.5*X(I,J) + TWO*X(I1P1,J) - HALF*X(I1P1+1,J)
      YXI1(J) = -1.5*Y(I,J) + TWO*Y(I1P1,J) - HALF*Y(I1P1+1,J)
      ZXI1(J) = -1.5*Z(I,J) + TWO*Z(I1P1,J) - HALF*Z(I1P1+1,J)

      I = I2

      XEI2(J) = -1.5*X(I,J) + TWO*X(I,J1P1) - HALF*X(I,J1P1+1)
      YEI2(J) = -1.5*Y(I,J) + TWO*Y(I,J1P1) - HALF*Y(I,J1P1+1)
      ZEI2(J) = -1.5*Z(I,J) + TWO*Z(I,J1P1) - HALF*Z(I,J1P1+1)

      ALP(I,J) = XEI2(J)**2 + YEI2(J)**2 + ZEI2(J)**2

      XEEI2(J) = TWO*X(I,J)-FIVE*X(I,J1P1)+FOUR*X(I,J1P1+1) -X(I,J1P1+2)
      YEEI2(J) = TWO*Y(I,J)-FIVE*Y(I,J1P1)+FOUR*Y(I,J1P1+1) -Y(I,J1P1+2)
      ZEEI2(J) = TWO*Z(I,J)-FIVE*Z(I,J1P1)+FOUR*Z(I,J1P1+1) -Z(I,J1P1+2)

      PSI = -(XEI2(J)*XEEI2(J) + YEI2(J)*YEEI2(J) +
     >        ZEI2(J)*ZEEI2(J)) / (ALP(I,J) + EPS)
      PSII2(J) = MIN (MAX (PSI, -FRCSI2), FRCSI2)

      XXI2(J)  = 1.5*X(I,J) - TWO*X(I2M1,J) + HALF*X(I2M1-1,J)
      YXI2(J)  = 1.5*Y(I,J) - TWO*Y(I2M1,J) + HALF*Y(I2M1-1,J)
      ZXI2(J)  = 1.5*Z(I,J) - TWO*Z(I2M1,J) + HALF*Z(I2M1-1,J)

      DO J = J1P1, J2M1
         XEI1(J) = (X(I1,J+1) - X(I1,J-1)) * HALF
         YEI1(J) = (Y(I1,J+1) - Y(I1,J-1)) * HALF
         ZEI1(J) = (Z(I1,J+1) - Z(I1,J-1)) * HALF
         XX = -1.5*X(I1,J) + TWO*X(I1P1,J) - HALF*X(I1P1+1,J)
         YX = -1.5*Y(I1,J) + TWO*Y(I1P1,J) - HALF*Y(I1P1+1,J)
         ZX = -1.5*Z(I1,J) + TWO*Z(I1P1,J) - HALF*Z(I1P1+1,J)

         GAM13 = YX * ZEI1(J) - ZX * YEI1(J)
         GAM23 = XEI1(J) * ZX - ZEI1(J) * XX
         GAM33 = XX * YEI1(J) - YX * XEI1(J)

         RALP33 = ONE / SQRT (GAM13**2 + GAM23**2 + GAM33**2 + EPS)
         ORTHO1 = GAM13*RALP33
         ORTHO2 = GAM23*RALP33
         ORTHO3 = GAM33*RALP33

         GAM11 = YEI1(J)*ORTHO3 - ZEI1(J)*ORTHO2
         GAM21 = ZEI1(J)*ORTHO1 - XEI1(J)*ORTHO3
         GAM31 = XEI1(J)*ORTHO2 - YEI1(J)*ORTHO1

         DQ = ONE / SQRT (GAM11**2 + GAM21**2 + GAM31**2 + EPS)

         XXI1(J) = (DQ * SI1(J)) * GAM11
         YXI1(J) = (DQ * SI1(J)) * GAM21
         ZXI1(J) = (DQ * SI1(J)) * GAM31

         ALP(I1,J) = XEI1(J)**2 + YEI1(J)**2 + ZEI1(J)**2
         BET(I1,J) = XEI1(J)*XXI1(J) + YEI1(J)*YXI1(J) +
     >               ZEI1(J)*ZXI1(J)
         GAM(I1,J) = XXI1(J)**2 + YXI1(J)**2 + ZXI1(J)**2

         AJAC(I1,J) = SQRT (ALP(I1,J)*GAM(I1,J) -
     >                      BET(I1,J)*BET(I1,J)) + EPS

         XEEI1(J) = X(I1,J+1) - TWO*X(I1,J) + X(I1,J-1)
         YEEI1(J) = Y(I1,J+1) - TWO*Y(I1,J) + Y(I1,J-1)
         ZEEI1(J) = Z(I1,J+1) - TWO*Z(I1,J) + Z(I1,J-1)

         CURV = SQRT (((XEI1(J)*YEEI1(J) - XEEI1(J)*YEI1(J))**2 +
     >                 (YEI1(J)*ZEEI1(J) - YEEI1(J)*ZEI1(J))**2 +
     >                 (ZEI1(J)*XEEI1(J) - ZEEI1(J)*XEI1(J))**2) /
     >                (ALP(I1,J) + EPS))

         FGMAXI1(J) = MIN (ONE, FGLIMI1 / (ONE + CURV))

         PSI = -(XEI1(J)*XEEI1(J) + YEI1(J)*YEEI1(J) +
     >           ZEI1(J)*ZEEI1(J)) / (ALP(I1,J) + EPS)
         PSII1(J) = MIN (MAX (PSI, -FRCSI1), FRCSI1) ! Background PSI at I1 edge

C        And for the I2 edge ...

         XEI2(J) = (X(I2,J+1) - X(I2,J-1)) * HALF
         YEI2(J) = (Y(I2,J+1) - Y(I2,J-1)) * HALF
         ZEI2(J) = (Z(I2,J+1) - Z(I2,J-1)) * HALF
         XX =  1.5*X(I2,J) - TWO*X(I2M1,J) + HALF*X(I2M1-1,J)
         YX =  1.5*Y(I2,J) - TWO*Y(I2M1,J) + HALF*Y(I2M1-1,J)
         ZX =  1.5*Z(I2,J) - TWO*Z(I2M1,J) + HALF*Z(I2M1-1,J)

         GAM13 = YX * ZEI2(J) - ZX * YEI2(J)
         GAM23 = XEI2(J) * ZX - ZEI2(J) * XX
         GAM33 = XX * YEI2(J) - YX * XEI2(J)

         RALP33 = ONE / SQRT (GAM13**2 + GAM23**2 + GAM33**2 + EPS)
         ORTHO1 = GAM13*RALP33
         ORTHO2 = GAM23*RALP33
         ORTHO3 = GAM33*RALP33

         GAM11 = YEI2(J)*ORTHO3 - ZEI2(J)*ORTHO2
         GAM21 = ZEI2(J)*ORTHO1 - XEI2(J)*ORTHO3
         GAM31 = XEI2(J)*ORTHO2 - YEI2(J)*ORTHO1

         DQ = ONE / SQRT (GAM11**2 + GAM21**2 + GAM31**2 + EPS)

         XXI2(J) = (DQ * SI2(J)) * GAM11
         YXI2(J) = (DQ * SI2(J)) * GAM21
         ZXI2(J) = (DQ * SI2(J)) * GAM31

         ALP(I2,J) = XEI2(J)**2 + YEI2(J)**2 + ZEI2(J)**2
         BET(I2,J) = XEI2(J)*XXI2(J) + YEI2(J)*YXI2(J) +
     >               ZEI2(J)*ZXI2(J)
         GAM(I2,J) = XXI2(J)**2 + YXI2(J)**2 + ZXI2(J)**2

         AJAC(I2,J) = SQRT (ALP(I2,J)*GAM(I2,J) -
     >                      BET(I2,J)*BET(I2,J)) + EPS

         XEEI2(J) = X(I2,J+1) - TWO*X(I2,J) + X(I2,J-1)
         YEEI2(J) = Y(I2,J+1) - TWO*Y(I2,J) + Y(I2,J-1)
         ZEEI2(J) = Z(I2,J+1) - TWO*Z(I2,J) + Z(I2,J-1)

         CURV = SQRT (((XEI2(J)*YEEI2(J) - XEEI2(J)*YEI2(J))**2 +
     >                 (YEI2(J)*ZEEI2(J) - YEEI2(J)*ZEI2(J))**2 +
     >                 (ZEI2(J)*XEEI2(J) - ZEEI2(J)*XEI2(J))**2) /
     >                (ALP(I2,J) + EPS))
         FGMAXI2(J) = MIN (ONE, FGLIMI2 / (ONE + CURV))

         PSI = -(XEI2(J)*XEEI2(J) + YEI2(J)*YEEI2(J) +
     >           ZEI2(J)*ZEEI2(J)) / (ALP(I2,J) + EPS)
         PSII2(J) = MIN (MAX (PSI, -FRCSI2), FRCSI2) ! Background PSI at I2 edge
      END DO

      I = I1
      J = J2

      XEI1(J) = 1.5*X(I,J) - TWO*X(I,J2M1) + HALF*X(I,J2M1-1)
      YEI1(J) = 1.5*Y(I,J) - TWO*Y(I,J2M1) + HALF*Y(I,J2M1-1)
      ZEI1(J) = 1.5*Z(I,J) - TWO*Z(I,J2M1) + HALF*Z(I,J2M1-1)

      ALP(I,J) = XEI1(J)**2 + YEI1(J)**2 + ZEI1(J)**2

      XEEI1(J) = TWO*X(I,J)-FIVE*X(I,J2M1)+FOUR*X(I,J2M1-1) -X(I,J2M1-2)
      YEEI1(J) = TWO*Y(I,J)-FIVE*Y(I,J2M1)+FOUR*Y(I,J2M1-1) -Y(I,J2M1-2)
      ZEEI1(J) = TWO*Z(I,J)-FIVE*Z(I,J2M1)+FOUR*Z(I,J2M1-1) -Z(I,J2M1-2)

      PSI = -(XEI1(J)*XEEI1(J) + YEI1(J)*YEEI1(J) +
     >        ZEI1(J)*ZEEI1(J)) / (ALP(I,J) + EPS)
      PSII1(J) = MIN (MAX (PSI, -FRCSI1), FRCSI1)

      XXI1(J) = -1.5*X(I,J) + TWO*X(I1P1,J) - HALF*X(I1P1+1,J)
      YXI1(J) = -1.5*Y(I,J) + TWO*Y(I1P1,J) - HALF*Y(I1P1+1,J)
      ZXI1(J) = -1.5*Z(I,J) + TWO*Z(I1P1,J) - HALF*Z(I1P1+1,J)

      I = I2

      XEI2(J)  = 1.5*X(I,J) - TWO*X(I,J2M1) + HALF*X(I,J2M1-1)
      YEI2(J)  = 1.5*Y(I,J) - TWO*Y(I,J2M1) + HALF*Y(I,J2M1-1)
      ZEI2(J)  = 1.5*Z(I,J) - TWO*Z(I,J2M1) + HALF*Z(I,J2M1-1)

      ALP(I,J) = XEI2(J)**2 + YEI2(J)**2 + ZEI2(J)**2

      XEEI2(J) = TWO*X(I,J)-FIVE*X(I,J2M1)+FOUR*X(I,J2M1-1) -X(I,J2M1-2)
      YEEI2(J) = TWO*Y(I,J)-FIVE*Y(I,J2M1)+FOUR*Y(I,J2M1-1) -Y(I,J2M1-2)
      ZEEI2(J) = TWO*Z(I,J)-FIVE*Z(I,J2M1)+FOUR*Z(I,J2M1-1) -Z(I,J2M1-2)

      PSI = -(XEI2(J)*XEEI2(J) + YEI2(J)*YEEI2(J) +
     >        ZEI2(J)*ZEEI2(J)) / (ALP(I,J) + EPS)
      PSII2(J) = MIN (MAX (PSI, -FRCSI2), FRCSI2)

      XXI2(J)  = 1.5*X(I,J) - TWO*X(I2M1,J) + HALF*X(I2M1-1,J)
      YXI2(J)  = 1.5*Y(I,J) - TWO*Y(I2M1,J) + HALF*Y(I2M1-1,J)
      ZXI2(J)  = 1.5*Z(I,J) - TWO*Z(I2M1,J) + HALF*Z(I2M1-1,J)

C     Cross-derivatives:

      DO J = J1P1, J2M1
         XEXI1(J) = (XXI1(J+1) - XXI1(J-1)) * HALF
         YEXI1(J) = (YXI1(J+1) - YXI1(J-1)) * HALF
         ZEXI1(J) = (ZXI1(J+1) - ZXI1(J-1)) * HALF
         XEXI2(J) = (XXI2(J+1) - XXI2(J-1)) * HALF
         YEXI2(J) = (YXI2(J+1) - YXI2(J-1)) * HALF
         ZEXI2(J) = (ZXI2(J+1) - ZXI2(J-1)) * HALF
      END DO

C     Smooth spikes or jumps in the background control functions to
C     reduce irregularities in the smoothed grid:

      CALL SMOOTH1D (I1P1, I2M1, SEDGEJ1(I1P1), PHIJ1(I1P1))
      CALL SMOOTH1D (I1P1, I2M1, SEDGEJ2(I1P1), PHIJ2(I1P1))
      CALL SMOOTH1D (J1P1, J2M1, SEDGEI1(J1P1), PSII1(J1P1))
      CALL SMOOTH1D (J1P1, J2M1, SEDGEI2(J1P1), PSII2(J1P1))


C     *********************
C     Begin SLOR iteration:
C     *********************

      P1 = ZERO ! Foreground parts of PHI, PSI at J1 edge
      Q1 = ZERO
      P2 = ZERO !  "   "   "   "    "    "    "   J2  "
      Q2 = ZERO
      P3 = ZERO
      Q3 = ZERO
      P4 = ZERO
      Q4 = ZERO

      ROMEGA = ONE / OMEGA
C**** NUMI   = I2M1 - I1P1 + 1
      CONV   = ZERO

      IF (PRINT) WRITE (6,1001) I1, I2, J1, J2, ITMAX, CONVG

      DO IT = 1, ITMAX ! Or until convergence

         DUMAX = ZERO
         DVMAX = ZERO
         DUSUM = ZERO
         DVSUM = ZERO

C        Set up Alpha, Gamma, and Beta for the RHS:

         DO I = I1P1, I2M1
            J  = J1
            XX = (X(I+1,J) - X(I-1,J)) * HALF
            YX = (Y(I+1,J) - Y(I-1,J)) * HALF
            ZX = (Z(I+1,J) - Z(I-1,J)) * HALF
            XN = -1.5*X(I,J) + TWO*X(I,J+1) - HALF*X(I,J+2)
            YN = -1.5*Y(I,J) + TWO*Y(I,J+1) - HALF*Y(I,J+2)
            ZN = -1.5*Z(I,J) + TWO*Z(I,J+1) - HALF*Z(I,J+2)

            ALP(I,J) = XN**2 + YN**2 + ZN**2
            BET(I,J) = XN*XX + YN*YX + ZN*ZX
            GAM(I,J) = XX**2 + YX**2 + ZX**2

            AJAC(I,J) = SQRT (ALP(I,J)*GAM(I,J) -
     >                        BET(I,J)*BET(I,J)) + EPS
            J  = J2
            XX = (X(I+1,J) - X(I-1,J)) * HALF
            YX = (Y(I+1,J) - Y(I-1,J)) * HALF
            ZX = (Z(I+1,J) - Z(I-1,J)) * HALF
            XN = 1.5*X(I,J) - TWO*X(I,J-1) + HALF*X(I,J-2)
            YN = 1.5*Y(I,J) - TWO*Y(I,J-1) + HALF*Y(I,J-2)
            ZN = 1.5*Z(I,J) - TWO*Z(I,J-1) + HALF*Z(I,J-2)

            ALP(I,J) = XN**2 + YN**2 + ZN**2
            BET(I,J) = XN*XX + YN*YX + ZN*ZX
            GAM(I,J) = XX**2 + YX**2 + ZX**2

            AJAC(I,J) = SQRT (ALP(I,J)*GAM(I,J) -
     >                        BET(I,J)*BET(I,J)) + EPS
         END DO

         DO J = J1P1, J2M1
            I  = I1
            XX = -1.5*X(I,J) + TWO*X(I+1,J) - HALF*X(I+2,J)
            YX = -1.5*Y(I,J) + TWO*Y(I+1,J) - HALF*Y(I+2,J)
            ZX = -1.5*Z(I,J) + TWO*Z(I+1,J) - HALF*Z(I+2,J)
            XN = (X(I,J+1) - X(I,J-1)) * HALF
            YN = (Y(I,J+1) - Y(I,J-1)) * HALF
            ZN = (Z(I,J+1) - Z(I,J-1)) * HALF

            ALP(I,J) = XN**2 + YN**2 + ZN**2
            BET(I,J) = XN*XX + YN*YX + ZN*ZX
            GAM(I,J) = XX**2 + YX**2 + ZX**2

            AJAC(I,J) = SQRT (ALP(I,J)*GAM(I,J) -
     >                        BET(I,J)*BET(I,J)) + EPS

            DO I = I1P1, I2M1
               XX = (X(I+1,J) - X(I-1,J)) * HALF
               YX = (Y(I+1,J) - Y(I-1,J)) * HALF
               ZX = (Z(I+1,J) - Z(I-1,J)) * HALF
               XN = (X(I,J+1) - X(I,J-1)) * HALF
               YN = (Y(I,J+1) - Y(I,J-1)) * HALF
               ZN = (Z(I,J+1) - Z(I,J-1)) * HALF

               ALP(I,J) = XN**2 + YN**2 + ZN**2
               BET(I,J) = XN*XX + YN*YX + ZN*ZX
               GAM(I,J) = XX**2 + YX**2 + ZX**2

               AJAC(I,J) = SQRT (ALP(I,J)*GAM(I,J) -
     >                           BET(I,J)*BET(I,J)) + EPS
            END DO

            I  = I2
            XX = 1.5*X(I,J) - TWO*X(I-1,J) + HALF*X(I-2,J)
            YX = 1.5*Y(I,J) - TWO*Y(I-1,J) + HALF*Y(I-2,J)
            ZX = 1.5*Z(I,J) - TWO*Z(I-1,J) + HALF*Z(I-2,J)
            XN = (X(I,J+1) - X(I,J-1)) * HALF
            YN = (Y(I,J+1) - Y(I,J-1)) * HALF
            ZN = (Z(I,J+1) - Z(I,J-1)) * HALF

            ALP(I,J) = XN**2 + YN**2 + ZN**2
            BET(I,J) = XN*XX + YN*YX + ZN*ZX
            GAM(I,J) = XX**2 + YX**2 + ZX**2

            AJAC(I,J) = SQRT (ALP(I,J)*GAM(I,J) -
     >                        BET(I,J)*BET(I,J)) + EPS
         END DO

C        CONVMIN allows for some degree of convergence before introducing the
C        foreground terms.

         IF (CONV >= CONVMIN) THEN

C           Add in foreground terms at the edges involving normal 2nd derivs.:

            IF (FRCFJ1 /= ZERO) THEN
               J = J1
               DO I = I1P1, I2M1
                  S1 = -3.5*X(I,J) + FOUR*X(I,J1P1) - HALF*X(I,J+2) -
     >                  3.0*XEJ1(I)
                  S2 = -3.5*Y(I,J) + FOUR*Y(I,J1P1) - HALF*Y(I,J+2) -
     >                  3.0*YEJ1(I)
                  S3 = -3.5*Z(I,J) + FOUR*Z(I,J1P1) - HALF*Z(I,J+2) -
     >                  3.0*ZEJ1(I)

                  R1 = ALP(I,J)*XXXJ1(I) - TWO*BET(I,J)*XEXJ1(I) +
     >                 GAM(I,J)*S1
                  R2 = ALP(I,J)*YXXJ1(I) - TWO*BET(I,J)*YEXJ1(I) +
     >                 GAM(I,J)*S2
                  R3 = ALP(I,J)*ZXXJ1(I) - TWO*BET(I,J)*ZEXJ1(I) +
     >                 GAM(I,J)*S3

C                 Limit the foreground magnitudes, and under-relax the changes:

                  T1   = R1*XXJ1(I) + R2*YXJ1(I) + R3*ZXJ1(I)
                  T2   = R1*XEJ1(I) + R2*YEJ1(I) + R3*ZEJ1(I)
                  RJAC = ONE / AJAC(I,J)**2

                  DELP = RJAC * (T2 * BET(I,J) / ALP(I,J) - T1) - P1(I)
                  DELQ = RJAC * (T1 * BET(I,J) / GAM(I,J) - T2) - Q1(I)

                  DELP = SIGN (MIN (URJ1 * ABS (DELP),
     >                   FGMAXJ1(I) * MAX (ABS (P1(I)), ONE)), DELP)
                  DELQ = SIGN (MIN (URJ1 * ABS (DELQ),
     >                   FGMAXJ1(I) * MAX (ABS (Q1(I)), ONE)), DELQ)

                  P1(I) = P1(I) + DELP
                  Q1(I) = Q1(I) + DELQ
               END DO
            END IF

            IF (FRCFJ2 /= ZERO) THEN
               J = J2
               DO I = I1P1, I2M1
                  S1 = -3.5*X(I,J) + FOUR*X(I,J2M1) - HALF*X(I,J-2) +
     >                  3.0*XEJ2(I)
                  S2 = -3.5*Y(I,J) + FOUR*Y(I,J2M1) - HALF*Y(I,J-2) +
     >                  3.0*YEJ2(I)
                  S3 = -3.5*Z(I,J) + FOUR*Z(I,J2M1) - HALF*Z(I,J-2) +
     >                  3.0*ZEJ2(I)

                  R1 = ALP(I,J)*XXXJ2(I) - TWO*BET(I,J)*XEXJ2(I) +
     >                 GAM(I,J)*S1
                  R2 = ALP(I,J)*YXXJ2(I) - TWO*BET(I,J)*YEXJ2(I) +
     >                 GAM(I,J)*S2
                  R3 = ALP(I,J)*ZXXJ2(I) - TWO*BET(I,J)*ZEXJ2(I) +
     >                 GAM(I,J)*S3

C                 Limit the foreground magnitudes, and under-relax the changes:

                  T1 = R1*XXJ2(I) + R2*YXJ2(I) + R3*ZXJ2(I)
                  T2 = R1*XEJ2(I) + R2*YEJ2(I) + R3*ZEJ2(I)
                  RJAC = ONE / AJAC(I,J)**2

                  DELP = RJAC * (T2 * BET(I,J) / ALP(I,J) - T1) - P2(I)
                  DELQ = RJAC * (T1 * BET(I,J) / GAM(I,J) - T2) - Q2(I)

                  DELP = SIGN (MIN (URJ2 * ABS (DELP),
     >                   FGMAXJ2(I) * MAX (ABS (P2(I)), ONE)), DELP)
                  DELQ = SIGN (MIN (URJ2 * ABS (DELQ),
     >                   FGMAXJ2(I) * MAX (ABS (Q2(I)), ONE)), DELQ)

                  P2(I) = P2(I) + DELP
                  Q2(I) = Q2(I) + DELQ
               END DO
            END IF

            IF (FRCFI1 /= ZERO) THEN
               I = I1
               DO J = J1P1, J2M1
                  S1 = -3.5*X(I,J) + FOUR*X(I1P1,J) - HALF*X(I+2,J) -
     >                  3.0*XXI1(J)
                  S2 = -3.5*Y(I,J) + FOUR*Y(I1P1,J) - HALF*Y(I+2,J) -
     >                  3.0*YXI1(J)
                  S3 = -3.5*Z(I,J) + FOUR*Z(I1P1,J) - HALF*Z(I+2,J) -
     >                  3.0*ZXI1(J)

                  R1 = ALP(I,J)*S1 - TWO*BET(I,J)*XEXI1(J) +
     >                 GAM(I,J)*XEEI1(J)
                  R2 = ALP(I,J)*S2 - TWO*BET(I,J)*YEXI1(J) +
     >                 GAM(I,J)*YEEI1(J)
                  R3 = ALP(I,J)*S3 - TWO*BET(I,J)*ZEXI1(J) +
     >                 GAM(I,J)*ZEEI1(J)

C                 Limit the foreground magnitudes, and under-relax the changes:

                  T1 = R1*XXI1(J) + R2*YXI1(J) + R3*ZXI1(J)
                  T2 = R1*XEI1(J) + R2*YEI1(J) + R3*ZEI1(J)
                  RJAC = ONE / AJAC(I,J)**2

                  DELP = RJAC * (T2 * BET(I,J) / ALP(I,J) - T1) - P3(J)
                  DELQ = RJAC * (T1 * BET(I,J) / GAM(I,J) - T2) - Q3(J)

                  DELP = SIGN (MIN (URI1 * ABS (DELP),
     >                   FGMAXI1(J) * MAX (ABS (P3(J)), ONE)), DELP)
                  DELQ = SIGN (MIN (URI1 * ABS (DELQ),
     >                   FGMAXI1(J) * MAX (ABS (Q3(J)), ONE)), DELQ)

                  P3(J) = P3(J) + DELP
                  Q3(J) = Q3(J) + DELQ
               END DO
            END IF

            IF (FRCFI2 /= ZERO) THEN
               I = I2
               DO J = J1P1, J2M1
                  S1 = -3.5*X(I,J) + FOUR*X(I2M1,J) - HALF*X(I-2,J) +
     >                  3.0*XXI2(J)
                  S2 = -3.5*Y(I,J) + FOUR*Y(I2M1,J) - HALF*Y(I-2,J) +
     >                  3.0*YXI2(J)
                  S3 = -3.5*Z(I,J) + FOUR*Z(I2M1,J) - HALF*Z(I-2,J) +
     >                  3.0*ZXI2(J)

                  R1 = ALP(I,J)*S1 - TWO*BET(I,J)*XEXI2(J) +
     >                 GAM(I,J)*XEEI2(J)
                  R2 = ALP(I,J)*S2 - TWO*BET(I,J)*YEXI2(J) +
     >                 GAM(I,J)*YEEI2(J)
                  R3 = ALP(I,J)*S3 - TWO*BET(I,J)*ZEXI2(J) +
     >                 GAM(I,J)*ZEEI2(J)

C                 Limit the foreground magnitudes, and under-relax the changes:

                  T1 = R1*XXI2(J) + R2*YXI2(J) + R3*ZXI2(J)
                  T2 = R1*XEI2(J) + R2*YEI2(J) + R3*ZEI2(J)
                  RJAC = ONE / AJAC(I,J)**2

                  DELP = RJAC * (T2 * BET(I,J) / ALP(I,J) - T1) - P4(J)
                  DELQ = RJAC * (T1 * BET(I,J) / GAM(I,J) - T2) - Q4(J)

                  DELP = SIGN (MIN (URI2 * ABS (DELP),
     >                   FGMAXI2(J) * MAX (ABS (P4(J)), ONE)), DELP)
                  DELQ = SIGN (MIN (URI2 * ABS (DELQ),
     >                   FGMAXI2(J) * MAX (ABS (Q4(J)), ONE)), DELQ)

                  P4(J) = P4(J) + DELP
                  Q4(J) = Q4(J) + DELQ
               END DO
            END IF
         END IF

C        Interior points by lines:

         DO 80 J = J1P1, J2M1
            FGJ1 = DECAYJ1(J)
            FGJ2 = DECAYJ2(J)
            DECAYJ = MAX (FGJ1, FGJ2)

            DO 60 I = I1P1, I2M1
               FGI1 = DECAYI1(I)
               FGI2 = DECAYI2(I)
               DECAYI = MAX (FGI1, FGI2)
               DECAY  = ONE - MAX (DECAYI, DECAYJ)

               RJ  = ARCS(I,J,1,2)
               PHI = (ONE - RJ)*PHIJ1(I) + RJ*PHIJ2(I) ! Background PHI

C              Foreground PHI: Only 1 edge contributes much, or 2 at a corner.

               DENOM = ONE / (DECAYJ + DECAYI + EPS)
               PHF = ((P1(I)*FGJ1 + P2(I)*FGJ2)*DECAYJ +
     >                (P3(J)*FGI1 + P4(J)*FGI2)*DECAYI) * DENOM

C              Blend foreground and background PHI:

               PHI = (DECAY * PHI + PHF) * ALP(I,J)
               SGNPHI = SIGN (ONE, PHI)

               RI  = ARCS(I,J,1,1)
               PSI = (ONE - RI)*PSII1(J) + RI*PSII2(J)          ! Bg. PSI
               PSF = ((Q1(I)*FGJ1 + Q2(I)*FGJ2)*DECAYJ +
     >                (Q3(J)*FGI1 + Q4(J)*FGI2)*DECAYI) * DENOM ! Fg. PSI

C              Blend foreground and background PSI:

               PSI = (DECAY * PSI + PSF) * GAM(I,J)
               SGNPSI = SIGN (ONE, PSI)

C              RHS terms used to keep grid on surface:

               WIPJ = HALF*(ALP(I+1,J)/AJAC(I+1,J) +
     >                      ALP(I,J)  /AJAC(I,J))
               WIJ = -HALF*(ALP(I+1,J)/AJAC(I+1,J) +
     >                  TWO*ALP(I,J)  /AJAC(I,J)   +
     >                      ALP(I-1,J)/AJAC(I-1,J) +
     >                      GAM(I,J+1)/AJAC(I,J+1) +
     >                  TWO*GAM(I,J)  /AJAC(I,J)   +
     >                      GAM(I,J-1)/AJAC(I,J-1))
               WIMJ = HALF*(ALP(I-1,J)/AJAC(I-1,J) +
     >                      ALP(I,J)  /AJAC(I,J))
               WIJP = HALF*(GAM(I,J+1)/AJAC(I,J+1) +
     >                      GAM(I,J)  /AJAC(I,J))
               WIJM = HALF*(GAM(I,J-1)/AJAC(I,J-1) +
     >                      GAM(I,J)  /AJAC(I,J))
               WIPJP =-.25*(BET(I+1,J)/AJAC(I+1,J) +
     >                      BET(I,J+1)/AJAC(I,J+1))
               WIPJM = .25*(BET(I+1,J)/AJAC(I+1,J) +
     >                      BET(I,J-1)/AJAC(I,J-1))
               WIMJP = .25*(BET(I-1,J)/AJAC(I-1,J) +
     >                      BET(I,J+1)/AJAC(I,J+1))
               WIMJM =-.25*(BET(I-1,J)/AJAC(I-1,J) +
     >                      BET(I,J-1)/AJAC(I,J-1))

C              Derivatives in the U and V directions:

               UXX =  U(I+1,J) - TWO*U(I,J) + U(I-1,J)
               VXX =  V(I+1,J) - TWO*V(I,J) + V(I-1,J)
               UXN = (U(I+1,J+1) - U(I-1,J+1) -
     >                U(I+1,J-1) + U(I-1,J-1)) * HALF ! 0.5 rather than 0.25
               VXN = (V(I+1,J+1) - V(I-1,J+1) -       ! allows factors of TWO
     >                V(I+1,J-1) + V(I-1,J-1)) * HALF ! to be removed below
               UNN =  U(I,J+1) - TWO*U(I,J) + U(I,J-1)
               VNN =  V(I,J+1) - TWO*V(I,J) + V(I,J-1)

               RHS1(I) = -ALP(I,J)*UXX + BET(I,J)*UXN -
     >                    GAM(I,J)*UNN - HALF *
     >            (PHI*((SGNPHI-ONE)*U(I-1,J) - TWO*SGNPHI*U(I,J) +
     >                  (SGNPHI+ONE)*U(I+1,J)) +
     >             PSI*((SGNPSI-ONE)*U(I,J-1) - TWO*SGNPSI*U(I,J) +
     >                  (SGNPSI+ONE)*U(I,J+1))) + AJAC(I,J)*
     >           (WIPJ*U(I+1,J)+WIPJP*U(I+1,J+1)+WIPJM*U(I+1,J-1) +
     >             WIJ*U(I,J)  +WIJP *U(I,J+1)  +WIJM *U(I,J-1)   +
     >            WIMJ*U(I-1,J)+WIMJP*U(I-1,J+1)+WIMJM*U(I-1,J-1))

               RHS2(I) = -ALP(I,J)*VXX + BET(I,J)*VXN -
     >                    GAM(I,J)*VNN - HALF *
     >            (PHI*((SGNPHI-ONE)*V(I-1,J) - TWO*SGNPHI*V(I,J) +
     >                  (SGNPHI+ONE)*V(I+1,J)) +
     >             PSI*((SGNPSI-ONE)*V(I,J-1) - TWO*SGNPSI*V(I,J) +
     >                  (SGNPSI+ONE)*V(I,J+1))) + AJAC(I,J)*
     >           (WIPJ*V(I+1,J)+WIPJP*V(I+1,J+1)+WIPJM*V(I+1,J-1) +
     >             WIJ*V(I,J)  +WIJP *V(I,J+1)  +WIJM *V(I,J-1)   +
     >            WIMJ*V(I-1,J)+WIMJP*V(I-1,J+1)+WIMJM*V(I-1,J-1))

C              The AJAC(I,J)*W terms are included in the subdiagonals
C              to maintain diagonal dominance. This is equivalent to
C              lagging rather than linearizing the Wij*Uij term.

               A(I) = (ALP(I,J) + (HALF*PHI)*(SGNPHI-ONE) -
     >                 AJAC(I,J)*WIMJ)*ROMEGA
               C(I) = (ALP(I,J) + (HALF*PHI)*(SGNPHI+ONE) -
     >                 AJAC(I,J)*WIPJ)*ROMEGA
               B(I) = (-TWO*(ALP(I,J) + GAM(I,J)) - (PHI*SGNPHI +
     >                 PSI*SGNPSI))*ROMEGA

   60       CONTINUE ! Next I


C           Solve two tridiagonal systems with a common LHS.
C           In-lining the solver improves efficiency:

C****       CALL TRID2R (NUMI, A(I1P1), B(I1P1), C(I1P1),
C****>                   RHS1(I1P1), RHS2(I1P1))

            W = ONE / B(I1P1)
            RHS1(I1P1) = W * RHS1(I1P1)
            RHS2(I1P1) = W * RHS2(I1P1)

            DO I = I1P1 + 1, I2M1
               B(I-1) = C(I-1) * W
               W = ONE / (B(I) - A(I) * B(I-1))
               RHS1(I) = W * (RHS1(I) - A(I) * RHS1(I-1))
               RHS2(I) = W * (RHS2(I) - A(I) * RHS2(I-1))
            END DO

            DO I = I2M1 - 1, I1P1, -1
               RHS1(I) = RHS1(I) - B(I) * RHS1(I+1)
               RHS2(I) = RHS2(I) - B(I) * RHS2(I+1)
            END DO

C           Determine the largest corrections (doesn't vectorize):

            DO I = I1P1, I2M1
               IF (DUMAX < ABS (RHS1(I))) THEN
                   DUMAX = ABS (RHS1(I))
                   IU = I
                   JU = J
               END IF
               IF (DVMAX < ABS (RHS2(I))) THEN
                   DVMAX = ABS (RHS2(I))
                   IV = I
                   JV = J
               END IF
            END DO

C           Update the grid interior (vectorizes):

            DO I = I1P1, I2M1
               U(I,J) = U(I,J) + RHS1(I)
               DUSUM  = DUSUM  + ABS (RHS1(I))
               V(I,J) = V(I,J) + RHS2(I)
               DVSUM  = DVSUM  + ABS (RHS2(I))
            END DO

   80    CONTINUE ! Next J

C        Bilinear interpolation of the original surface mesh at the new (U,V)s:

         DO J = J1P1, J2M1
            DO I = I1P1, I2M1
               IINT = MAX (I1, MIN (I2M1, INT (U(I,J))))
               PINT = U(I,J) - REAL (IINT)
               PM1  = ONE - PINT
               JINT = MAX (J1, MIN (J2M1, INT (V(I,J))))
               QINT = V(I,J) - REAL (JINT)
               QM1  = ONE - QINT
               X(I,J) = QM1  * (PM1  * XORIG(IINT,JINT)    +
     >                          PINT * XORIG(IINT+1,JINT)) +
     >                  QINT * (PM1  * XORIG(IINT,JINT+1)  +
     >                          PINT * XORIG(IINT+1,JINT+1))
               Y(I,J) = QM1  * (PM1  * YORIG(IINT,JINT)    +
     >                          PINT * YORIG(IINT+1,JINT)) +
     >                  QINT * (PM1  * YORIG(IINT,JINT+1)  +
     >                          PINT * YORIG(IINT+1,JINT+1))
               Z(I,J) = QM1  * (PM1  * ZORIG(IINT,JINT)    +
     >                          PINT * ZORIG(IINT+1,JINT)) +
     >                  QINT * (PM1  * ZORIG(IINT,JINT+1)  +
     >                          PINT * ZORIG(IINT+1,JINT+1))
            END DO
         END DO

         IF (PRINT) WRITE (6, 1002)
     >      IT, DUSUM, DUMAX, IU, JU, DVSUM, DVMAX, IV, JV

         IF (IT <= 1) THEN
            DUM1 = DUMAX
            DVM1 = DVMAX
         END IF

         CONV = LOG10 (MAX (DUM1 / DUMAX, DVM1 / DVMAX))

         IF (CONV > CONVG .OR. MAX (DUMAX, DVMAX) < DMAX) EXIT

      END DO ! Next SLOR iteration

  999 RETURN

 1001 FORMAT (/, ' I1, I2:', 2I4, '  J1, J2:', 2I4,
     >   '   Iteration limit:', I5,
     >   '   Orders of magnitude:', F5.2, //, ' ITER',
     >   12X, 'DUSUM', 12X, 'DUMAX    I    J',
     >   12X, 'DVSUM', 12X, 'DVMAX    I    J', /)
 1002 FORMAT (I5, 2(1P, 2E17.7, 2I5))

      END SUBROUTINE ELLIPQ3D
