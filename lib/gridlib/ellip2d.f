C+------------------------------------------------------------------------------
C
      SUBROUTINE ELLIP2D (IDIM, JDIM, I1, I2, J1, J2, X, Y, ARCS, PRINT,
     >                    BGMODE, FGMODE, SPMODE, ITMAX, ITFLOAT,
     >                    ITFREEZE, JLAYER, CONVG, CONVMIN, DMAX, OMEGA,
     >                    POWERI, POWERJ, EXPI1, EXPI2, EXPJ1, EXPJ2,
     >                    FGLIMIT, URFG, URFLOAT)
C
C     ELLIP2D smooths the interior of the specified 2-D (sub)grid with
C     control of spacing and edge orthogonality based on the formulation
C     of Sorenson, where P and Q are also referred to as PHI and PSI:
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
C     Poisson's equation is solved via SLOR.
C
C     The initial spacing off a given edge may either be derived from the end-
C     point increments at adjoining edges (only), or it may preserve the input
C     first increments along the entire edge.  Provision is also made for
C     floating the outer boundary distribution (J = J2 only at present).
C
C     This is a major revision of the routine TTM2D found in the WBGRID wing/
C     body grid generator. Initial efforts to introduce control of orthogonality
C     as described by Thompson, et al. have been replaced with the supposedly
C     equivalent scheme described in the GRAPE user guide (Sorenson).
C
C     References:  Sorenson (May 1980) NASA Tech. Mem. 81198
C                  Thompson, Warsi, Mastin (1985), p. 231
C
C     Typical inputs for half of an Euler-type C-grid on an airfoil:
C
C     BGMODE = 'YYYY', FGMODE = 'YYYY', SPMODE = 'YYYY',
C     ITMAX = 200, ITFLOAT = 0, ITFREEZE = 100, JLAYER = 4,
C     CONVG = 6., CONVMIN = 0., DMAX = 0.0001, OMEGA = 1.2,
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
C     8/96       "   "    FGLIMIT and URFG arguments serve for all edges now;
C                         local FGLIMI1, etc. could easily be made arguments
C                         again if necessary; ARCS(*,*,2) has to be an argument
C                         with dimensions of X & Y to suit the PARAMXY call.
C     6-7/97     "   "    Control of orthogonality may now be arc-length-based.
C     8/97       "   "    POWERI < 0. extension needed for symmetric C grids.
C     04/27/98    DAS     Fortran 90 allows use of automatic arrays and
C                         internal procedures.
C     05/08/98     "      Applied SMOOTH1D to the background edge distributions.
C     09/08/99     "      Mark Rimlinger found NEW was not initialized.
C     10/04/99     "      IMPLICIT NONE belatedly used to help trace a failure
C                         on a Linux system.
C ------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) ::
     >   IDIM, JDIM,        ! Calling program dimensions for X, Y, & ARCS
     >   I1, I2, J1, J2     ! Active index range for X, Y, & ARCS;
                            ! I2, J2 <= IDIM, JDIM

      REAL, INTENT (INOUT), DIMENSION (IDIM, JDIM) ::
     >   X, Y               ! Grid coordinates, input with initial estimates
                                  !
      REAL, INTENT (OUT) ::
     >   ARCS(IDIM,JDIM,2)  ! Work-space for relative arc lengths;
                            ! dimensions must match X, Y in PARAMXY

      LOGICAL, INTENT (IN) ::
     >   PRINT              ! .TRUE. shows iterations on unit 6

      CHARACTER * 4, INTENT (IN) ::
     >   BGMODE,            ! 4 upper case characters referring to edges
                            ! I1, I2, J1, J2 in that order;
                            ! 'Y' activates "background" control of
                            !     interior spacing from that edge;
                            ! 'N' suppresses background control
     >   FGMODE,            ! 'N' suppresses "foreground" control of
                            !     orthogonality at that edge;
                            ! 'Y' activates the index-based scheme for
                            !     exponential decay;
                            ! 'A' activates the arc-lgth-based scheme,
                            !     meaning foreground term coefs. decay
                            !     exponentially from 1.0 at the edge to 0.1 at
                            !     the fractional distance input as EXPI1, etc.
     >   SPMODE             ! 'Y' means interior initial spacing increments
                            !     are controlled from the corner increments;
                            ! 'N' means starting guess values are retained
                            !     for the entire boundary

      INTEGER, INTENT (IN) ::
     >   ITMAX,             ! Limit on number of smoothing iterations
     >   ITFLOAT,           ! 0 turns off the option to float the outer edge
                            !   (J2 only at present)
     >   ITFREEZE,          ! An edge is updated every ITFLOAT iterations until
                            ! it is frozen for full convergence after ITFREEZE
                            ! iterations
     >   JLAYER             ! If ITFLOAT > 0, the relative I distribution along
                            ! the line J = J2 - JLAYER is imposed at J = J2

      REAL, INTENT (IN) ::
     >   CONVG,             ! Orders of magnitude reduction sought in the
                            ! maximum change in dX or dY
     >   CONVMIN,           ! CONVMIN > 0. suppresses foreground terms until
                            ! some degree of convergence is achieved
     >   DMAX,              ! Allows early termination if dX and dY are all
                            ! below this value
     >   OMEGA,             ! SLOR over-relaxation parameter
     >   POWERI, POWERJ,    !   1. means linear combination of opposite edge
                            !      PHI/PSI for background control;
                            ! > 1. emphasizes the lower edge;
                            ! < 1.     "       "  upper edge;
                            !      use POWERI/J < 0. to switch meaning of
                            !      "lower" and "upper" - e.g., POWERI < 0.
                            !      for the upper half of a typical C grid
     >   EXPI1, EXPI2,      ! Decay factors for foreground terms;
     >   EXPJ1, EXPJ2,      ! if FGMODE(n:n) = 'Y', larger means more rapid
                            ! decay - try 0.40 - 0.50;
                            ! if FGMODE(n:n) = 'A', input the fraction of arc
                            ! length from the edge at which the foreground
                            ! coef. should decay to 0.1
     >   FGLIMIT,           ! Limiter on the foreground terms;
                            ! FGLIMIT / (1 + edge curvature) is used
     >   URFG,              ! Under-relaxation factor applied to the changes
                            ! in the foreground terms
     >   URFLOAT            ! Under-relaxes the boundary floating

C-------------------------------------------------------------------------------

C     Local constants:

      REAL, PARAMETER ::
     >   EPS = 1.E-16,     ! Safeguards 0/0
     >   RLN = 2.302585,   ! exp(-RLN) = 0.1
     >   FOUR = 4., HALF = 0.5, ONE = 1., TWO = 2., ZERO = 0.

      LOGICAL, PARAMETER ::
     >   NEW = .TRUE.

      CHARACTER * 1, PARAMETER ::
     >   MONO = 'M', NO = 'N', YES = 'Y'

C     Local variables:

      INTEGER
     >   I,       I1P1,    I2M1,    IT,      IX,      IY,
     >   J,       J1P1,    J2M1,    JX,      JY,      JIN,
     >   I1POWER, I2POWER, J1POWER, J2POWER, NCOUNT,  NEDGE,   NUMI

      REAL
     >   AJ2,     ALP,     BET,     CONV,    CURV,
     >   D1I1J1,  D1I2J1,  D1I1J2,  D1I2J2,
     >   D1J1I1,  D1J2I1,  D1J1I2,  D1J2I2,
     >   DECAY,   DECAYI,  DECAYJ,  DELP,    DELQ,    DENOM,   DQ,
     >   DXMAX,   DXM1,    DXSUM,   DYMAX,   DYM1,    DYSUM,
     >   FGI1,    FGI2,    FGJ1,    FGJ2,
     >   FGLIMI1, FGLIMI2, FGLIMJ1, FGLIMJ2,
     >   FRCFI1,  FRCFI2,  FRCFJ1,  FRCFJ2, 
     >   FRCSI1,  FRCSI2,  FRCPJ1,  FRCPJ2,
     >   GAM,     PHF,     PHI,     PSF,     PSI,
     >   R1,      R2,      RI,      RI2MI1,  RJ,      RJ2MJ1,  RJAC,
     >   ROMEGA,  RSEDGE,  RSINNER, S1,      S2,      SGNPHI,  SGNPSI,
     >   URI1,    URI2,    URJ1,    URJ2,    W,
     >   XN, XNN, XX, XXN, XXX, YN, YNN, YX, YXN, YXX

C     Automatic arrays:

      REAL, DIMENSION (I1:I2) ::
     >   A, B, C, FGMAXJ1, FGMAXJ2, P1, P2, Q1, Q2, PHIJ1, PHIJ2,
     >   RHS1, RHS2, SEDGEJ1, SEDGEJ2, SINNER, SJ1, SJ2, XNEW, YNEW,
     >   XXJ1, YXJ1, XXXJ1, YXXJ1, XEJ1, YEJ1, XEXJ1, YEXJ1,
     >   XXJ2, YXJ2, XXXJ2, YXXJ2, XEJ2, YEJ2, XEXJ2, YEXJ2

      REAL, DIMENSION (J1:J2) ::
     >   FGMAXI1, FGMAXI2, P3, P4, Q3, Q4, PSII1, PSII2,
     >   SEDGEI1, SEDGEI2, SI1, SI2,
     >   XEI1, YEI1, XEEI1, YEEI1, XXI1, YXI1, XEXI1, YEXI1,
     >   XEI2, YEI2, XEEI2, YEEI2, XXI2, YXI2, XEXI2, YEXI2

      REAL, DIMENSION (I1:I2, J1:J2) ::
     >   DECAYI1, DECAYI2, DECAYJ1, DECAYJ2

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

      CALL PARAMXY (1, IDIM, 1, JDIM, I1, I2, J1, J2, X, Y, ARCS)

C     Edge arc lengths are needed for explicit smoothing of background fns.

      SEDGEJ1(I1:I2) = ARCS(I1:I2,J1,1)
      SEDGEJ2(I1:I2) = ARCS(I1:I2,J2,1)

      SEDGEI1(J1:J2) = ARCS(I1,J1:J2,2)
      SEDGEI2(J1:J2) = ARCS(I2,J1:J2,2)

C     Store the foreground decay coefficients:

      I1POWER = NINT (-LOG10 (MIN (ARCS(I1P1,J1,1), ARCS(I1P1,J2,1))))
      I1POWER = MAX ((I1POWER + 1)/ 2, 3)
      I2POWER = NINT (-LOG10 (ONE - MAX (ARCS(I2M1,J1,1),
     >                                   ARCS(I2M1,J2,1))))
      I2POWER = MAX ((I2POWER + 1)/ 2, 3)

      J1POWER = NINT (-LOG10 (MIN (ARCS(I1,J1P1,2), ARCS(I2,J1P1,2))))
      J1POWER = MAX ((J1POWER + 1)/ 2, 3)
      J2POWER = NINT (-LOG10 (ONE - MAX (ARCS(I1,J2M1,2),
     >                                   ARCS(I2,J2M1,2))))
      J2POWER = MAX ((J2POWER + 1)/ 2, 3)

      IF (FGMODE(1:1) == NO) THEN ! No foreground

         FRCFI1  = ZERO
         DECAYI1 = ZERO

      ELSE IF (FGMODE(1:1) == YES) THEN ! Index-based decay

         FRCFI1 = ONE

         DO I = I1, I2
            DECAYI1(I,J1) = (RI2MI1 * REAL (I2 - I)) ** I1POWER
     >                * EXP (-EXPI1 * REAL (I - I1))
         END DO

         DO J = J1P1, J2
            DECAYI1(I1:I2,J) = DECAYI1(I1:I2,J1)
         END DO

      ELSE ! 'A' for arc-length-based decay

         FRCFI1 = ONE
         DECAY  = RLN / EXPI1

         DO J = J1, J2
            DO I = I1, I2
               XX = DECAY * ARCS(I,J,1)
               DECAYI1(I,J) = (RI2MI1 * REAL (I2 - I)) ** I1POWER
     >                      * EXP (-XX)
            END DO
         END DO

      END IF

      IF (FGMODE(2:2) == NO) THEN

         FRCFI2  = ZERO
         DECAYI2 = ZERO

      ELSE IF (FGMODE(2:2) == YES) THEN

         FRCFI2 = ONE

         DO I = I1, I2
            DECAYI2(I,J1) = (RI2MI1 * REAL (I - I1)) ** I2POWER
     >                * EXP (-EXPI2 * REAL (I2 - I))
         END DO

         DO J = J1P1, J2
            DECAYI2(I1:I2,J) = DECAYI2(I1:I2,J1)
         END DO

      ELSE

         FRCFI2 = ONE
         DECAY  = RLN / EXPI2

         DO J = J1, J2
            DO I = I1, I2
               XX = DECAY * (ONE - ARCS(I,J,1))
               DECAYI2(I,J) = (RI2MI1 * REAL (I - I1)) ** I2POWER
     >                      * EXP (-XX)
            END DO
         END DO

      END IF

      IF (FGMODE(3:3) == NO) THEN

         FRCFJ1  = ZERO
         DECAYJ1 = ZERO

      ELSE IF (FGMODE(3:3) == YES) THEN

         FRCFJ1 = ONE

         DO J = J1, J2
            DECAYJ1(I1,J) = (RJ2MJ1 * REAL (J2 - J)) ** J1POWER
     >                * EXP (-EXPJ1 * REAL (J - J1))
         END DO

         DO I = I1P1, I2
            DECAYJ1(I,J1:J2) = DECAYJ1(I1,J1:J2)
         END DO

      ELSE

         FRCFJ1 = ONE
         DECAY  = RLN / EXPJ1

         DO J = J1, J2
            DO I = I1, I2
               XX = DECAY * ARCS(I,J,2)
               DECAYJ1(I,J) = (RJ2MJ1 * REAL (J2 - J)) ** J1POWER
     >                      * EXP(-XX)
            END DO
         END DO

      END IF

      IF (FGMODE(4:4) == NO) THEN

         FRCFJ2  = ZERO
         DECAYJ2 = ZERO

      ELSE IF (FGMODE(4:4) == YES) THEN

         FRCFJ2 = ONE

         DO J = J1, J2
            DECAYJ2(I1,J) = (RJ2MJ1 * REAL (J - J1)) ** J2POWER
     >                * EXP (-EXPJ2 * REAL (J2 - J))
         END DO

         DO I = I1P1, I2
            DECAYJ2(I,J1:J2) = DECAYJ2(I1,J1:J2)
         END DO

      ELSE

         FRCFJ2 = ONE
         DECAY  = RLN / EXPJ2

         DO J = J1, J2
            DO I = I1, I2
               XX = DECAY * (ONE - ARCS(I,J,2))
               DECAYJ2(I,J) = (RJ2MJ1 * REAL (J - J1)) ** J2POWER
     >                      * EXP (-XX)
            END DO
         END DO

      END IF


C     Moderating the INTERIOR relative arc-lengths (POWER* < 1.)
C     helps background control if the curvature at one boundary
C     greatly exceeds the curvature at the opposite boundary.

      CALL MODERATE_ARCS ()


C     D1MN = distance increment from corner (M,N) along edge M:

      D1I1J1 = SQRT ((X(I1,J1P1)-X(I1,J1))**2 +(Y(I1,J1P1)-Y(I1,J1))**2)
      D1I2J1 = SQRT ((X(I2,J1P1)-X(I2,J1))**2 +(Y(I2,J1P1)-Y(I2,J1))**2)
      D1I1J2 = SQRT ((X(I1,J2M1)-X(I1,J2))**2 +(Y(I1,J2M1)-Y(I1,J2))**2)
      D1I2J2 = SQRT ((X(I2,J2M1)-X(I2,J2))**2 +(Y(I2,J2M1)-Y(I2,J2))**2)

      D1J1I1 = SQRT ((X(I1P1,J1)-X(I1,J1))**2 +(Y(I1P1,J1)-Y(I1,J1))**2)
      D1J2I1 = SQRT ((X(I1P1,J2)-X(I1,J2))**2 +(Y(I1P1,J2)-Y(I1,J2))**2)
      D1J1I2 = SQRT ((X(I2M1,J1)-X(I2,J1))**2 +(Y(I2M1,J1)-Y(I2,J1))**2)
      D1J2I2 = SQRT ((X(I2M1,J2)-X(I2,J2))**2 +(Y(I2M1,J2)-Y(I2,J2))**2)


C     Set up control of initial spacing off each edge at each edge point:

      CALL EDGE_SPACING (.TRUE.)


C     Set up the J1/J2 edge values of the background contributions to PHI,
C     & the foreground terms involving edge-wise derivatives.

      XEJ1(I1) = -1.5*X(I1,J1) + TWO*X(I1,J1P1) - HALF*X(I1,J1P1+1)
      YEJ1(I1) = -1.5*Y(I1,J1) + TWO*Y(I1,J1P1) - HALF*Y(I1,J1P1+1)
      XEJ2(I1) =  1.5*X(I1,J2) - TWO*X(I1,J2M1) + HALF*X(I1,J2M1-1)
      YEJ2(I1) =  1.5*Y(I1,J2) - TWO*Y(I1,J2M1) + HALF*Y(I1,J2M1-1)

      P1(I1P1:I2M1) = ZERO ! Parts of PHI, PSI at J1 edge
      Q1(I1P1:I2M1) = ZERO !
      P2(I1P1:I2M1) = ZERO !   "    "    "    "   J2  "
      Q2(I1P1:I2M1) = ZERO

      DO I = I1P1, I2M1
         XXJ1(I) = (X(I+1,J1) - X(I-1,J1)) * HALF
         YXJ1(I) = (Y(I+1,J1) - Y(I-1,J1)) * HALF

         DQ = ONE / SQRT (XXJ1(I)**2 + YXJ1(I)**2 + EPS)
         XEJ1(I) = -(DQ * SJ1(I)) * YXJ1(I)
         YEJ1(I) =  (DQ * SJ1(I)) * XXJ1(I)

         XXXJ1(I) = X(I+1,J1) - TWO*X(I,J1) + X(I-1,J1)
         YXXJ1(I) = Y(I+1,J1) - TWO*Y(I,J1) + Y(I-1,J1)

C        Handle edge discontinuities or high curvature via foreground limiting:

         CURV = (XXJ1(I) * YXXJ1(I) - YXJ1(I) * XXXJ1(I)) * DQ**3
         FGMAXJ1(I) = MIN (ONE, FGLIMJ1 / (ONE + ABS (CURV)))

         PHI = -(XXJ1(I) * XXXJ1(I) + YXJ1(I) * YXXJ1(I)) * DQ**2
         PHIJ1(I) = MIN (MAX (PHI, -FRCPJ1), FRCPJ1) ! Background PHI at J1 edge

C        And for the J2 edge ...

         XXJ2(I) = (X(I+1,J2) - X(I-1,J2)) * HALF
         YXJ2(I) = (Y(I+1,J2) - Y(I-1,J2)) * HALF

         DQ = ONE / SQRT (XXJ2(I)**2 + YXJ2(I)**2 + EPS)
         XEJ2(I) = -(DQ * SJ2(I)) * YXJ2(I)
         YEJ2(I) =  (DQ * SJ2(I)) * XXJ2(I)

         XXXJ2(I) = X(I+1,J2) - TWO*X(I,J2) + X(I-1,J2)
         YXXJ2(I) = Y(I+1,J2) - TWO*Y(I,J2) + Y(I-1,J2)

         CURV = (XXJ2(I) * YXXJ2(I) - YXJ2(I) * XXXJ2(I)) * DQ**3
         FGMAXJ2(I) = MIN (ONE, FGLIMJ2 / (ONE + ABS (CURV)))

         PHI = -(XXJ2(I) * XXXJ2(I) + YXJ2(I) * YXXJ2(I)) * DQ**2
         PHIJ2(I) = MIN (MAX (PHI, -FRCPJ2), FRCPJ2) ! Background PHI at J2 edge
      END DO

      XEJ1(I2) = -1.5*X(I2,J1) + TWO*X(I2,J1P1) - HALF*X(I2,J1P1+1)
      YEJ1(I2) = -1.5*Y(I2,J1) + TWO*Y(I2,J1P1) - HALF*Y(I2,J1P1+1)
      XEJ2(I2) =  1.5*X(I2,J2) - TWO*X(I2,J2M1) + HALF*X(I2,J2M1-1)
      YEJ2(I2) =  1.5*Y(I2,J2) - TWO*Y(I2,J2M1) + HALF*Y(I2,J2M1-1)

C     Cross-derivatives:

      DO I = I1P1, I2M1
         XEXJ1(I) = (XEJ1(I+1) - XEJ1(I-1)) * HALF
         YEXJ1(I) = (YEJ1(I+1) - YEJ1(I-1)) * HALF
         XEXJ2(I) = (XEJ2(I+1) - XEJ2(I-1)) * HALF
         YEXJ2(I) = (YEJ2(I+1) - YEJ2(I-1)) * HALF
      END DO

C     Likewise for the I1 and I2 edges ...

      XXI1(J1) = -1.5*X(I1,J1) + TWO*X(I1P1,J1) - HALF*X(I1P1+1,J1)
      YXI1(J1) = -1.5*Y(I1,J1) + TWO*Y(I1P1,J1) - HALF*Y(I1P1+1,J1)
      XXI2(J1) =  1.5*X(I2,J1) - TWO*X(I2M1,J1) + HALF*X(I2M1-1,J1)
      YXI2(J1) =  1.5*Y(I2,J1) - TWO*Y(I2M1,J1) + HALF*Y(I2M1-1,J1)

      P3(J1P1:J2M1) = ZERO
      Q3(J1P1:J2M1) = ZERO
      P4(J1P1:J2M1) = ZERO
      Q4(J1P1:J2M1) = ZERO

      DO J = J1P1, J2M1
         XEI1(J) = (X(I1,J+1) - X(I1,J-1)) * HALF
         YEI1(J) = (Y(I1,J+1) - Y(I1,J-1)) * HALF

         DQ = ONE / SQRT (XEI1(J)**2 + YEI1(J)**2 + EPS)
         XXI1(J) =  (DQ * SI1(J)) * YEI1(J)
         YXI1(J) = -(DQ * SI1(J)) * XEI1(J)

         XEEI1(J) = X(I1,J+1) - TWO*X(I1,J) + X(I1,J-1)
         YEEI1(J) = Y(I1,J+1) - TWO*Y(I1,J) + Y(I1,J-1)

         CURV = (XEI1(J)*YEEI1(J) - YEI1(J)*XEEI1(J)) * DQ**3
         FGMAXI1(J) = MIN (ONE, FGLIMI1 / (ONE + ABS (CURV)))

         PSI = -(XEI1(J)*XEEI1(J) + YEI1(J)*YEEI1(J)) * DQ**2
         PSII1(J) = MIN (MAX (PSI*FRCSI1, -ONE), ONE) ! Backgd. PSI at I1 edge

C        And for the I2 edge ...

         XEI2(J) = (X(I2,J+1) - X(I2,J-1)) * HALF
         YEI2(J) = (Y(I2,J+1) - Y(I2,J-1)) * HALF

         DQ = ONE / SQRT (XEI2(J)**2 + YEI2(J)**2 + EPS)
         XXI2(J) =  (DQ * SI2(J)) * YEI2(J)
         YXI2(J) = -(DQ * SI2(J)) * XEI2(J)

         XEEI2(J) = X(I2,J+1) - TWO*X(I2,J) + X(I2,J-1)
         YEEI2(J) = Y(I2,J+1) - TWO*Y(I2,J) + Y(I2,J-1)

         CURV = (XEI2(J)*YEEI2(J) - YEI2(J)*XEEI2(J)) * DQ**3
         FGMAXI2(J) = MIN (ONE, FGLIMI2 / (ONE + ABS (CURV)))
         PSI = -(XEI2(J)*XEEI2(J) + YEI2(J)*YEEI2(J)) * DQ**2
         PSII2(J) = MIN (MAX (PSI*FRCSI2, -ONE), ONE) ! Backgd. PSI at I2 edge
      END DO

      XXI1(J2) = -1.5*X(I1,J2) + TWO*X(I1P1,J2) - HALF*X(I1P1+1,J2)
      YXI1(J2) = -1.5*Y(I1,J2) + TWO*Y(I1P1,J2) - HALF*Y(I1P1+1,J2)
      XXI2(J2) =  1.5*X(I2,J2) - TWO*X(I2M1,J2) + HALF*X(I2M1-1,J2)
      YXI2(J2) =  1.5*Y(I2,J2) - TWO*Y(I2M1,J2) + HALF*Y(I2M1-1,J2)

C     Cross-derivatives:

      DO J = J1P1, J2M1
         XEXI1(J) = (XXI1(J+1) - XXI1(J-1)) * HALF
         YEXI1(J) = (YXI1(J+1) - YXI1(J-1)) * HALF
         XEXI2(J) = (XXI2(J+1) - XXI2(J-1)) * HALF
         YEXI2(J) = (YXI2(J+1) - YXI2(J-1)) * HALF
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

      ROMEGA = ONE / OMEGA
      NUMI   = I2M1 - I1P1 + 1
      NEDGE  = NUMI + 2
      JIN    = J2 - JLAYER
      NCOUNT = ITFLOAT ! NCOUNT counts backwards
      CONV   = ZERO

      IF (PRINT) WRITE (6,1001) I1, I2, J1, J2, ITMAX, CONVG

      DO IT = 1, ITMAX ! Or until convergence

         DXMAX = ZERO
         DYMAX = ZERO
         DXSUM = ZERO
         DYSUM = ZERO

         IF (IT == ITFREEZE) NCOUNT = 0 ! Stops floating of boundary

C        Float the J2 boundary?  (Other boundaries not implemented yet.)

         IF (NCOUNT == 1) THEN

            SEDGEJ2(I1) = ZERO
            SINNER(I1)  = ZERO
            DO I = I1P1, I2
               SEDGEJ2(I) = SEDGEJ2(I-1) + SQRT
     >            ((X(I,J2) - X(I-1,J2))**2 + (Y(I,J2) - Y(I-1,J2))**2)
               SINNER(I) = SINNER(I-1) + SQRT
     >            ((X(I,JIN)-X(I-1,JIN))**2 + (Y(I,JIN)-Y(I-1,JIN))**2)
            END DO

C           Impose the distribution at J = J2 - JLAYER on the J2 edge:

            RSEDGE  = URFLOAT / SEDGEJ2(I2) ! Under-relax the changes
            RSINNER = URFLOAT / SINNER(I2)

            DO I = I1P1, I2M1
               SINNER(I) = SEDGEJ2(I) + SEDGEJ2(I2) *
     >           (SINNER(I) * RSINNER - SEDGEJ2(I) * RSEDGE)
            END DO

            CALL LCSFIT (NEDGE, SEDGEJ2(I1), X(I1,J2), NEW, MONO,
     >                   NUMI, SINNER(I1P1), XNEW(I1P1), XNEW(I1P1))

            CALL LCSFIT (NEDGE, SEDGEJ2(I1), Y(I1,J2), NEW, MONO,
     >                   NUMI, SINNER(I1P1), YNEW(I1P1), YNEW(I1P1))

            X(I1P1:I2M1,J2) = XNEW(I1P1:I2M1)
            Y(I1P1:I2M1,J2) = YNEW(I1P1:I2M1)

C           Updating all normalized arc lengths shouldn't hurt and may help:

            CALL PARAMXY (1, IDIM, 1, JDIM, I1, I2, J1, J2, X, Y, ARCS)

            SEDGEJ2(I1:I2) = ARCS(I1:I2,J2,1)

C           Adjust the arc lengths as done initially:

            CALL MODERATE_ARCS ()

            D1J2I1 = SQRT ((X(I1P1,J2) - X(I1,J2))**2 +
     >                     (Y(I1P1,J2) - Y(I1,J2))**2)
            D1J2I2 = SQRT ((X(I2M1,J2) - X(I2,J2))**2 +
     >                     (Y(I2M1,J2) - Y(I2,J2))**2)

C           Update any corner-based off-the-edge spacing controls:

            CALL EDGE_SPACING (.FALSE.)

C           Update the J2 edge derivatives:

            DO I = I1P1, I2M1
               XXJ2(I) = (X(I+1,J2) - X(I-1,J2)) * HALF
               YXJ2(I) = (Y(I+1,J2) - Y(I-1,J2)) * HALF

               DQ = ONE / SQRT (XXJ2(I)**2 + YXJ2(I)**2 + EPS)

               XEJ2(I) = -(DQ * SJ2(I)) * YXJ2(I)
               YEJ2(I) =  (DQ * SJ2(I)) * XXJ2(I)

               XXXJ2(I) = X(I+1,J2) - TWO*X(I,J2) + X(I-1,J2)
               YXXJ2(I) = Y(I+1,J2) - TWO*Y(I,J2) + Y(I-1,J2)

               CURV = (XXJ2(I) * YXXJ2(I) - YXJ2(I) * XXXJ2(I)) * DQ**3
               FGMAXJ2(I) = MIN (ONE, FGLIMJ2 / (ONE + ABS (CURV)))

               PHI = -(XXJ2(I) * XXXJ2(I) + YXJ2(I) * YXXJ2(I)) * DQ**2
               PHIJ2(I) = MIN (MAX (PHI, -FRCPJ2), FRCPJ2)
            END DO

            CALL SMOOTH1D (I1P1, I2M1, SEDGEJ2(I1P1), PHIJ2(I1P1))

            DO I = I1P1, I2M1
               XEXJ2(I) = (XEJ2(I+1) - XEJ2(I-1)) * HALF
               YEXJ2(I) = (YEJ2(I+1) - YEJ2(I-1)) * HALF
            END DO

            NCOUNT = ITFLOAT
         END IF

         NCOUNT = NCOUNT - 1

C        CONVMIN allows for some degree of convergence before introducing the
C        foreground terms.  *** Is it worth retaining???

         IF (CONV >= CONVMIN) THEN

C           Add in foreground terms at the edges involving normal 2nd derivs.:

            IF (FRCFJ1 /= ZERO) THEN
               DO I = I1P1, I2M1
                  S1 = -3.5*X(I,J1) + FOUR*X(I,J1P1) - HALF*X(I,J1P1+1)
     >                 -3.0*XEJ1(I)
                  S2 = -3.5*Y(I,J1) + FOUR*Y(I,J1P1) - HALF*Y(I,J1P1+1)
     >                 -3.0*YEJ1(I)

                  ALP  =        XEJ1(I)*XEJ1(I) + YEJ1(I)*YEJ1(I)
                  BET  =        XXJ1(I)*XEJ1(I) + YXJ1(I)*YEJ1(I)
                  GAM  =        XXJ1(I)*XXJ1(I) + YXJ1(I)*YXJ1(I)
                  RJAC = ONE / (XXJ1(I)*YEJ1(I) - XEJ1(I)*YXJ1(I) + EPS)

                  R1 =-(ALP*XXXJ1(I) -TWO*BET*XEXJ1(I) +GAM*S1)*RJAC**2
                  R2 =-(ALP*YXXJ1(I) -TWO*BET*YEXJ1(I) +GAM*S2)*RJAC**2

C                 Limit the foreground magnitudes, and under-relax the changes:

                  DELP = ( YEJ1(I)*R1 - XEJ1(I)*R2) * RJAC - P1(I)
                  DELQ = (-YXJ1(I)*R1 + XXJ1(I)*R2) * RJAC - Q1(I)

                  DELP = SIGN (MIN (URJ1 * ABS (DELP),
     >                   FGMAXJ1(I) * MAX (ABS (P1(I)), ONE)), DELP)
                  DELQ = SIGN (MIN (URJ1 * ABS (DELQ),
     >                   FGMAXJ1(I) * MAX (ABS (Q1(I)), ONE)), DELQ)

                  P1(I) = P1(I) + DELP
                  Q1(I) = Q1(I) + DELQ
               END DO
            END IF

            IF (FRCFJ2 /= ZERO) THEN
               DO I = I1P1, I2M1
                  S1 = -3.5*X(I,J2) + FOUR*X(I,J2M1) - HALF*X(I,J2M1-1)
     >                 +3.0*XEJ2(I)
                  S2 = -3.5*Y(I,J2) + FOUR*Y(I,J2M1) - HALF*Y(I,J2M1-1)
     >                 +3.0*YEJ2(I)

                  ALP  =        XEJ2(I)*XEJ2(I) + YEJ2(I)*YEJ2(I)
                  BET  =        XXJ2(I)*XEJ2(I) + YXJ2(I)*YEJ2(I)
                  GAM  =        XXJ2(I)*XXJ2(I) + YXJ2(I)*YXJ2(I)
                  RJAC = ONE / (XXJ2(I)*YEJ2(I) - XEJ2(I)*YXJ2(I) + EPS)

                  R1 =-(ALP*XXXJ2(I) -TWO*BET*XEXJ2(I) +GAM*S1)*RJAC**2
                  R2 =-(ALP*YXXJ2(I) -TWO*BET*YEXJ2(I) +GAM*S2)*RJAC**2

                  DELP = ( YEJ2(I)*R1 - XEJ2(I)*R2) * RJAC - P2(I)
                  DELQ = (-YXJ2(I)*R1 + XXJ2(I)*R2) * RJAC - Q2(I)

                  DELP = SIGN (MIN (URJ2 * ABS (DELP),
     >                   FGMAXJ2(I) * MAX (ABS (P2(I)), ONE)), DELP)
                  DELQ = SIGN (MIN (URJ2 * ABS (DELQ),
     >                   FGMAXJ2(I) * MAX (ABS (Q2(I)), ONE)), DELQ)

                  P2(I) = P2(I) + DELP
                  Q2(I) = Q2(I) + DELQ
               END DO
            END IF

            IF (FRCFI1 /= ZERO) THEN
               DO J = J1P1, J2M1
                  S1 = -3.5*X(I1,J) + FOUR*X(I1P1,J) - HALF*X(I1P1+1,J)
     >                 -3.0*XXI1(J)
                  S2 = -3.5*Y(I1,J) + FOUR*Y(I1P1,J) - HALF*Y(I1P1+1,J)
     >                 -3.0*YXI1(J)

                  ALP  =        XEI1(J)*XEI1(J) + YEI1(J)*YEI1(J)
                  BET  =        XXI1(J)*XEI1(J) + YXI1(J)*YEI1(J)
                  GAM  =        XXI1(J)*XXI1(J) + YXI1(J)*YXI1(J)
                  RJAC = ONE / (XXI1(J)*YEI1(J) - XEI1(J)*YXI1(J) + EPS)

                  R1 =-(ALP*S1 -TWO*BET*XEXI1(J) +GAM*XEEI1(J))*RJAC**2
                  R2 =-(ALP*S2 -TWO*BET*YEXI1(J) +GAM*YEEI1(J))*RJAC**2

                  DELP = ( YEI1(J)*R1 - XEI1(J)*R2) * RJAC - P3(J)
                  DELQ = (-YXI1(J)*R1 + XXI1(J)*R2) * RJAC - Q3(J)

                  DELP = SIGN (MIN (URI1 * ABS (DELP),
     >                   FGMAXI1(J) * MAX (ABS (P3(J)), ONE)), DELP)
                  DELQ = SIGN (MIN (URI1 * ABS (DELQ),
     >                   FGMAXI1(J) * MAX (ABS (Q3(J)), ONE)), DELQ)

                  P3(J) = P3(J) + DELP
                  Q3(J) = Q3(J) + DELQ
               END DO
            END IF

            IF (FRCFI2 /= ZERO) THEN
               DO J = J1P1, J2M1
                  S1 = -3.5*X(I2,J) + FOUR*X(I2M1,J) - HALF*X(I2M1-1,J)
     >                 +3.0*XXI2(J)
                  S2 = -3.5*Y(I2,J) + FOUR*Y(I2M1,J) - HALF*Y(I2M1-1,J)
     >                 +3.0*YXI2(J)

                  ALP  =        XEI2(J)*XEI2(J) + YEI2(J)*YEI2(J)
                  BET  =        XXI2(J)*XEI2(J) + YXI2(J)*YEI2(J)
                  GAM  =        XXI2(J)*XXI2(J) + YXI2(J)*YXI2(J)
                  RJAC = ONE / (XXI2(J)*YEI2(J) - XEI2(J)*YXI2(J) + EPS)

                  R1 =-(ALP*S1 -TWO*BET*XEXI2(J) +GAM*XEEI2(J))*RJAC**2
                  R2 =-(ALP*S2 -TWO*BET*YEXI2(J) +GAM*YEEI2(J))*RJAC**2

                  DELP = ( YEI2(J)*R1 - XEI2(J)*R2) * RJAC - P4(J)
                  DELQ = (-YXI2(J)*R1 + XXI2(J)*R2) * RJAC - Q4(J)

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

         DO J = J1P1, J2M1

            DO I = I1P1, I2M1
               FGJ1 = DECAYJ1(I,J)
               FGJ2 = DECAYJ2(I,J)
               DECAYJ = MAX (FGJ1, FGJ2)

               FGI1 = DECAYI1(I,J)
               FGI2 = DECAYI2(I,J)
               DECAYI = MAX (FGI1, FGI2)
               DECAY  = ONE - MAX (DECAYI, DECAYJ)

               XX  = (X(I+1,J) - X(I-1,J)) * HALF
               YX  = (Y(I+1,J) - Y(I-1,J)) * HALF
               XN  = (X(I,J+1) - X(I,J-1)) * HALF
               YN  = (Y(I,J+1) - Y(I,J-1)) * HALF
               XXX =  X(I+1,J) - TWO*X(I,J) + X(I-1,J)
               YXX =  Y(I+1,J) - TWO*Y(I,J) + Y(I-1,J)
               XXN = (X(I+1,J+1) - X(I-1,J+1) -
     >                X(I+1,J-1) + X(I-1,J-1)) * HALF ! 0.5 rather than 0.25
               YXN = (Y(I+1,J+1) - Y(I-1,J+1) -       ! allows factors of 2.
     >                Y(I+1,J-1) + Y(I-1,J-1)) * HALF ! to be removed below
               XNN =  X(I,J+1) - TWO*X(I,J) + X(I,J-1)
               YNN =  Y(I,J+1) - TWO*Y(I,J) + Y(I,J-1)
               ALP =  XN*XN + YN*YN
               BET =  XX*XN + YX*YN
               GAM =  XX*XX + YX*YX
               AJ2 = (XX*YN - YX*XN)**2

               RJ  = ARCS(I,J,2)
               PHI = ((ONE - RJ)*PHIJ1(I) + RJ*PHIJ2(I)) * ALP ! Background PHI

C              Foreground PHI: Only 1 edge contributes much, or 2 at a corner.

               DENOM = ONE / (DECAYI + DECAYJ + EPS)
               PHF = ((P1(I)*FGJ1 + P2(I)*FGJ2)*DECAYJ +
     >                (P3(J)*FGI1 + P4(J)*FGI2)*DECAYI) * DENOM

C              Blend foreground and background PHI:

               PHI = DECAY * PHI + PHF * AJ2
               SGNPHI = SIGN (ONE, PHI)

               RI  = ARCS(I,J,1)
               PSI = ((ONE - RI)*PSII1(J) + RI*PSII2(J)) * GAM  ! Background PSI
               PSF = ((Q1(I)*FGJ1 + Q2(I)*FGJ2)*DECAYJ +
     >                (Q3(J)*FGI1 + Q4(J)*FGI2)*DECAYI) * DENOM ! Fg. PSI

C              Blend foreground and background PSI:

               PSI = DECAY * PSI + PSF * AJ2
               SGNPSI = SIGN (ONE, PSI)

               RHS1(I) = -ALP*XXX + BET*XXN - GAM*XNN - HALF*
     >            (PHI*((SGNPHI-ONE)*X(I-1,J) - TWO*SGNPHI*X(I,J) +
     >                  (SGNPHI+ONE)*X(I+1,J)) +
     >             PSI*((SGNPSI-ONE)*X(I,J-1) - TWO*SGNPSI*X(I,J) +
     >                  (SGNPSI+ONE)*X(I,J+1)))

               RHS2(I) = -ALP*YXX + BET*YXN - GAM*YNN - HALF*
     >            (PHI*((SGNPHI-ONE)*Y(I-1,J) - TWO*SGNPHI*Y(I,J) +
     >                  (SGNPHI+ONE)*Y(I+1,J)) +
     >             PSI*((SGNPSI-ONE)*Y(I,J-1) - TWO*SGNPSI*Y(I,J) +
     >                  (SGNPSI+ONE)*Y(I,J+1)))

               A(I) = (ALP + (HALF*PHI)*(SGNPHI-ONE))*ROMEGA
               C(I) = (ALP + (HALF*PHI)*(SGNPHI+ONE))*ROMEGA
               B(I) = (-TWO*(ALP + GAM) - (PHI*SGNPHI + PSI*SGNPSI))*
     >                ROMEGA

            END DO ! Next I


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
               IF (DXMAX < ABS (RHS1(I))) THEN
                   DXMAX = ABS (RHS1(I))
                   IX = I
                   JX = J
               END IF
               IF (DYMAX < ABS (RHS2(I))) THEN
                   DYMAX = ABS (RHS2(I))
                   IY = I
                   JY = J
               END IF
            END DO

C           Update the grid interior (vectorizes):

            DO I = I1P1, I2M1
               X(I,J) = X(I,J) + RHS1(I)
               DXSUM  = DXSUM  + ABS (RHS1(I))
               Y(I,J) = Y(I,J) + RHS2(I)
               DYSUM  = DYSUM  + ABS (RHS2(I))
            END DO

         END DO ! Next J

         IF (PRINT) WRITE (6, 1002)
     >      IT, DXSUM, DXMAX, IX, JX, DYSUM, DYMAX, IY, JY

         IF (IT <= 1) THEN
            DXM1 = DXMAX
            DYM1 = DYMAX
         END IF

         CONV = LOG10 (MAX (DXM1 / DXMAX, DYM1 / DYMAX))

         IF (CONV > CONVG .OR. MAX (DXMAX, DYMAX) < DMAX) EXIT

      END DO ! Next SLOR iteration

  999 RETURN


C     Formats:

 1001 FORMAT (/, ' I1, I2:', 2I4, '  J1, J2:', 2I4,
     >   '   Iteration limit:', I5,
     >   '   Orders of magnitude:', F5.2, //, ' ITER',
     >   12X, 'DXSUM', 12X, 'DXMAX    I    J',
     >   12X, 'DYSUM', 12X, 'DYMAX    I    J', /)
 1002 FORMAT (I5, 2(1P, 2E17.7, 2I5))


C     ELLIP2D internal procedures:
C     ----------------------------

      CONTAINS

      SUBROUTINE MODERATE_ARCS ()

!     Moderating the INTERIOR relative arc-lengths (POWER* < 1.)
!     helps background control if the curvature at one boundary
!     greatly exceeds the curvature at the opposite boundary.

      IF (ABS (POWERI) /= ONE) THEN
         IF (POWERI == HALF) THEN
            DO J = J1P1, J2M1
               DO I = I1P1, I2M1
                  ARCS(I,J,1) = SQRT (ARCS(I,J,1))
               END DO
            END DO
         ELSE IF (POWERI == -HALF) THEN
            DO J = J1P1, J2M1
               DO I = I1P1, I2M1
                  ARCS(I,J,1) = ONE - SQRT (ONE - ARCS(I,J,1))
               END DO
            END DO
         ELSE IF (POWERI > ZERO) THEN
            DO J = J1P1, J2M1
               DO I = I1P1, I2M1
                  ARCS(I,J,1) = ARCS(I,J,1) ** POWERI
               END DO
            END DO
         ELSE ! POWERI < ZERO
            DO J = J1P1, J2M1
               DO I = I1P1, I2M1
                  ARCS(I,J,1) = ONE - (ONE - ARCS(I,J,1)) ** (-POWERI)
               END DO
            END DO
         END IF
      END IF

      IF (ABS (POWERJ) /= ONE) THEN
         IF (POWERJ == HALF) THEN
            DO J = J1P1, J2M1
               DO I = I1P1, I2M1
                  ARCS(I,J,2) = SQRT (ARCS(I,J,2))
               END DO
            END DO
         ELSE IF (POWERJ == -HALF) THEN
            DO J = J1P1, J2M1
               DO I = I1P1, I2M1
                  ARCS(I,J,2) = ONE - SQRT (ONE - ARCS(I,J,2))
               END DO
            END DO
         ELSE IF (POWERJ > ZERO) THEN
            DO J = J1P1, J2M1
               DO I = I1P1, I2M1
                  ARCS(I,J,2) = ARCS(I,J,2) ** POWERJ
               END DO
            END DO
         ELSE ! POWERJ < ZERO
            DO J = J1P1, J2M1
               DO I = I1P1, I2M1
                  ARCS(I,J,2) = ONE - (ONE - ARCS(I,J,2)) ** (-POWERJ)
               END DO
            END DO
         END IF
      END IF

      END SUBROUTINE MODERATE_ARCS


      SUBROUTINE EDGE_SPACING (FIRST)

      LOGICAL FIRST ! If SPMODE = N, don't update starting guess controls 

!     Is off-boundary spacing controlled by end-pt. edge increments only (Y)
!     or preserved from the starting guess for the entire boundary?

      IF (SPMODE(1:1) == YES) THEN
         DO J = J1, J2
            RJ = ARCS(I1,J,2)
            SI1(J) = (ONE - RJ) * D1J1I1 + RJ * D1J2I1
         END DO
      ELSE IF (FIRST) THEN
         DO J = J1, J2
            SI1(J) = SQRT ((X(I1P1,J) - X(I1,J))**2 +
     >                     (Y(I1P1,J) - Y(I1,J))**2)
         END DO
      END IF

      IF (SPMODE(2:2) == YES) THEN
         DO J = J1, J2
            RJ = ARCS(I2,J,2)
            SI2(J) = (ONE - RJ) * D1J1I2 + RJ * D1J2I2
         END DO
      ELSE IF (FIRST) THEN
         DO J = J1, J2
            SI2(J) = SQRT ((X(I2,J) - X(I2M1,J))**2 +
     >                     (Y(I2,J) - Y(I2M1,J))**2)
         END DO
      END IF

      IF (SPMODE(3:3) == YES) THEN
         DO I = I1, I2
            RI = ARCS(I,J1,1)
            SJ1(I) = (ONE - RI) * D1I1J1 + RI * D1I2J1
         END DO
      ELSE IF (FIRST) THEN
         DO I = I1, I2
            SJ1(I) = SQRT ((X(I,J1P1) - X(I,J1))**2 +
     >                     (Y(I,J1P1) - Y(I,J1))**2)
         END DO
      END IF

      IF (SPMODE(4:4) == YES) THEN
         DO I = I1, I2
            RI = ARCS(I,J2,1)
            SJ2(I) = (ONE - RI) * D1I1J2 + RI * D1I2J2
         END DO
      ELSE IF (FIRST) THEN
         DO I = I1, I2
            SJ2(I) = SQRT ((X(I,J2) - X(I,J2M1))**2 +
     >                     (Y(I,J2) - Y(I,J2M1))**2)
         END DO
      END IF

      END SUBROUTINE EDGE_SPACING

      END SUBROUTINE ELLIP2D
