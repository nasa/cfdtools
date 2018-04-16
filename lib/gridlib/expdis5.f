C-------------------------------------------------------------------------------
C
      SUBROUTINE EXPDIS5 (MODE, XA, XB, DX, N, X, LUNOUT)
C
C  DESCRIPTION:
C
C        EXPDIS5 is a variant of EXPDIS2/4 that safeguards all possibilities.
C     Instead of assuming the input DX is the smallest output interval, it
C     assumes that DX is desired at the end XA or XB according to MODE.  This
C     means that the spacing may increase or decrease from the indicated end,
C     depending on how DX compares with the uniform spacing, DU.  If DX is too
C     close to DU from below, a geometric distribution is generated, and this
C     may be combined in a weighted average with the intended distribution for
C     a small range of DX such as [0.9 DU, 0.95 DU].  See B1 below for the
C     actual blend window.  [0.94 DU, 0.97 DU] was settled on at the time of
C     implementing the blending, along with quadratic weighting.  Variations
C     of the window and the type of weighting all produce irregular changes in
C     the spacing at the opposite end for regular changes at the indicated end,
C     but this is believed to be inevitable given the different character of
C     the two types of distribution involved.  Note that if either iteration
C     fails, a uniform distribution is resorted to, but it is virtually certain
C     that neither can fail for the safeguarded cases permitted.
C
C        The original EXPDIS* routine was so-named before it was realized
C     that its stretching is equivalent to the one-sided Vinokur stretching.
C
C        Common usage performs one-sided stretching first with this routine,
C    followed by two-sided Vinokur stretching (see HTDIS4 or VINOKUR) with
C    some fraction of the resulting DX2 (e.g., DX2 * 0.2) as the spacing to
C    use at the other end of interval [XA, XB].  [Afterthought:  Is picking
C    that fraction well any easier than defining the second end spacing well
C    directly?]
C
C        Note that the EXPDIS4 meanings of MODE = 1 and 2 have been changed to
C     the more natural order, and the symmetric case has been dropped, since
C     two-sided stretching is the right choice for the same spacing at each end.
C
C        MODE = 1 means impose DX [or possibly DU] at the XA end;
C        MODE = 2   "     "     "     "     "     "     " XB end.
C
C  FURTHER NOTES:
C
C        EXPDIS2/4/5 apply the analysis of EXPDIS to solve the inverse problem:
C     "What BETA gives a specified DX?" as opposed to "What DX does a given
C     BETA produce?" by using a zero-finder to solve the relevant nonlinear
C     equation in BETA.  However:
C
C        As indicated above, the Vinokur method's zero-finding iteration fails
C     if DX is too close to DU.  In fact, it fails if DX > DU - epsilon for some
C     small epsilon.  The method of EXPDIS2/4 requires the spacing to INCREASE
C     sufficiently from the target DX at the indicated end of [XA, XB].
C
C        Sufficient conditions for the existence of a solution for the one-sided
C     Vinokur stretching have not been determined, except that BETA = 1. cor-
C     sponds to zero spacing at the low end, while uniform spacing corresponds
C     to infinite BETA.  The search interval for BETA is [1.00000001, 10.], 
C     which is believed to cover all likely applications.
C
C        As already stated, where the Vinokur stretching is guaranteed or likely
C     to fail, a geometric distribution is employed instead, with blending of
C     the two if DX is in a small window below DU.
C
C  ARGUMENTS:
C     ARG    TYPE  I/O/S  DIM   DESCRIPTION
C     MODE     I     I     -    MODE = 1 means apply DX at the XA end;
C                               MODE = 2 means apply DX at the XB end; see above
C     XA       R     I     -    Left- and right-hand ends of the total interval
C     XB       R     I     -    to be distributed.
C     DX       R     I     -    Desired length of the spacing at the indicated
C                               end of [XA, XB] (but see above).
C     N        I     I     -    Required number of points; N > 1.
C     X        R     O     N    1-D coordinates varying smoothly from XA to XB,
C                               with steady stretching in either direction (or
C                               strictly uniform if DX is too close to uniform)
C     LUNOUT   I     I     -    LUNOUT < 0 suppresses iteration history printing
C
C  ERROR HANDLING:
C
C        Bad inputs for XA, XB, N, and DX (such as XA > XB) are not checked for.
C     A zero finder failure leads to returning the uniform distribution, with no
C     error message, unlike EXPDIS4, which also proceeds with BETA = 1.00000001.
C
C  PROCEDURES:
C
C     ZERORC   General-purpose reverse-communication 1-D zero-finder
C
C  HISTORY:
C
C  08/21/85    DAS    Original implementation of EXPDIS (after a formula
C                     found in the SCRAM2D flow solver of A.J.Kumar).
C  02/25/88     "     EXPDIS2 derived from the analysis of EXPDIS.
C  08/09/89     "     Calculation of BETA is in DOUBLE PRECISION now to
C                     help requirement for initial dX of 1.E-7 on [0,1].
C                     If ZEROIN fails, keep going with very small BETA.
C  05/04/91     "     Description updated slightly in light of HTDIS2, etc.
C  12/01/94     "     Substituted ZERORC for ZEROIN to avoid COMMON and to
C                     avoid having both forms in one application.
C   "   "       "     EXPDIS4 adapted from EXPDIS2 for 64-bit systems.
C  08/12/09     "     EXPDIS5 adapted from EXPDIS4 to safeguard all likely
C                     combinations of XA, XB, N, and DX.
C  08/13/09     "     Where the spacing needs to decrease, geometric variation
C                     is preferable to an inverted Vinokur-type distribution.
C                     Matt Bartkowitz recommended this for cases where Navier-
C                     Stokes-type stretching at the end away from DX is hardly
C                     what is wanted.  Geometric keeps the variation as smooth
C                     as possible.  In extreme cases, relaxing the specified DX
C                     may be better yet, but that still has to be done at a
C                     higher level.
C  10/14/09     "     Todd White suggested blending the two types of point
C                     distribution somehow so that as radial grid lines in a
C                     grid block vary across the block, results near where
C                     DX ~ DU vary more smoothly.  This version does that.
C
C  AUTHOR: David Saunders, Sterling Software/NASA Ames Research Center, CA
C          (now with ELORET Corporation at ARC).
C
C-------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN)  :: LUNOUT, MODE, N
      REAL,    INTENT (IN)  :: DX, XA, XB
      REAL,    INTENT (OUT) :: X (N)

C     Local constants:

      INTEGER, PARAMETER ::
     >   BLEND_E = 2,             ! Exponent for blending function
     >   MAXFUN  = 35             ! Max. # zero-finder function evaluations

      REAL, PARAMETER ::
     >   HALF    = 5.E-1,
     >   ONE     = 1.E+0,
     >   TWO     = 2.E+0,
     >   ZERO    = 0.E+0,
     >   BETAMAX = 10.E+0,        ! Search interval for the value of BETA
     >   BETAMIN = 1.00000001E+0, ! which produces the desired DX
     >   B1      = 9.4E-1,        ! Blend region is [DU*B1, DU*B2]
     >   B2      = (B1 + ONE)*HALF,
     >   TOL     = 0.E+0          ! Ask ZERORC for full precision

      CHARACTER, PARAMETER :: SUBNAME *7 = 'EXPDIS5'

C     Local variables:

      INTEGER ::
     >   I, ISTAT, NUMFUN
      REAL ::
     >   B, BETA, DELTA, DSB1, DSB2, DSBM, DU, DXNORM, FUN, FI, GAMMA,
     >   HALF_BLEND, HOLD (13), POWER, RANGE, REALN, RNM1, WG, WV
      LOGICAL ::
     >   NEED_GEOMETRIC, NEED_VINOKUR, UNIFORM

      REAL, ALLOCATABLE :: DELTAX (:), XG (:)

C     Procedures:

      EXTERNAL :: GEODIS, ZERORC

C     Execution:

      RANGE = XB - XA
      REALN = REAL (N)
      DU    = ONE / (REALN - ONE)
      POWER = (REALN - TWO) * DU
      DU    = RANGE * DU            ! Uniform interval, full scale
      DSB1  = DU * B1               ! Pure Vinokur if DS < DSB1
      DSB2  = DU * B2               ! Pure geometric if DS > DSB2

      NEED_GEOMETRIC = DX > DSB1    ! Blend if both of these are true
      NEED_VINOKUR   = DX < DSB2
      UNIFORM        = .FALSE.

      IF (NEED_VINOKUR) THEN        ! Spacing grows sufficiently from DX end

         DELTA   = DX
         ISTAT   = 2                ! Initialize zero-finder to solve for BETA
         NUMFUN  = MAXFUN

   10    CONTINUE

            CALL ZERORC (BETAMIN, BETAMAX, BETA, FUN, TOL, NUMFUN,
     >                   SUBNAME, LUNOUT, HOLD, ISTAT)

            IF (ISTAT > 0) THEN                 ! Evaluate the function

               GAMMA = ((BETA + ONE) / (BETA - ONE)) ** POWER
               FUN = (ONE - BETA * (GAMMA - ONE) / (GAMMA + ONE)) *
     >               RANGE - DELTA
               GO TO 10

            ELSE IF (ISTAT < 0) THEN            ! -1 => MAXFUN reached; usable

               IF (ISTAT < -1) UNIFORM = .TRUE. ! -2, -3 => bad inputs, but
                                                ! -4 should mean DX too near DU
CC             LUNERR = ABS (LUNOUT)            ! No - we don't want print for
CC             WRITE (LUNERR, 90) ISTAT, XA, XB, DX, N, MODE  ! ~trouble cases

C*****      ELSE ISTAT == 0                     ! A zero has been found for Beta

            END IF

      END IF

      IF (NEED_GEOMETRIC) THEN  ! Spacing decreases or is too close to uniform

         ALLOCATE (XG (N))

         CALL GEODIS (XA, XB, N, DX, ZERO, XG, LUNOUT, ISTAT)  ! For DX at XA
                                   ! ZERO power gives pure geometric spacing

         UNIFORM = ISTAT /= 0  ! Should never happen, but use uniform if it does

      END IF

      IF (UNIFORM) THEN

         DO I = 1, N
            X (I) = XA + DU * REAL (I - 1)
         END DO

      ELSE  ! Err on the most likely side of pure Vinokur

         IF (NEED_VINOKUR) THEN  ! Spacing increases from specified end

C           Note that for F(I) = A ** P(I), where A is independent of I,
C           it is cheaper to use the form  F(I) = EXP (LOG(A) * P(I)):

            B = LOG ((BETA + ONE) / (BETA - ONE))
            RNM1 = B / REAL (N - 1)

            IF (MODE == 1) THEN   ! DX specified at XA

               DO I = 1, N
                  FI = EXP (RNM1 *  REAL (N - I))
                  X (I) = RANGE * (ONE - BETA * (FI - ONE) / (FI + ONE))
     >                    + XA
               END DO

            ELSE       ! MODE = 2;  DX specified at XB

               DO I = 1, N
                  FI = EXP (RNM1 * REAL (I - 1))
                  X (I) = (RANGE * BETA) * (FI - ONE) / (FI + ONE) + XA
               END DO

            END IF

         END IF

         IF (NEED_GEOMETRIC) THEN ! Check for DX at XB, not XA

            IF (MODE /= 1) THEN   ! Reverse GEODIS distribution for DX at XB

               ALLOCATE (DELTAX (N))

               DO I = 1, N - 1
                  DELTAX (I) = XB - XG (I)  ! Distance from XB
               END DO

               DO I = 2, N
                  XG (I) = XA + DELTAX (N + 1 - I)
               END DO

               DEALLOCATE (DELTAX)

            END IF

            IF (.NOT. NEED_VINOKUR) THEN

               X (:) = XG (:)

            ELSE  ! We need to blend the two distributions

               DSBM       = HALF * (DSB1 + DSB2)  ! Mid-point of blend interval
               HALF_BLEND = HALF * (DSB2 - DSB1)  ! Half width of blend interval

               IF (DX <= DSBM) THEN
                  DXNORM = (DX - DSB1) / HALF_BLEND
                  WG = HALF * DXNORM ** BLEND_E
               ELSE  ! Inverted/reflected form of DXNORM ** BLEND_E
                  DXNORM = (DX - DSBM) / HALF_BLEND
                  WG = ONE - HALF * (ONE - DXNORM) ** BLEND_E
               END IF

               WV = ONE - WG

CCC            write (*, '(a, 2f12.8)') 'Blending with WV, WG = ', WV,WG

               X (:) = WV * X (:) + WG * XG (:)

            END IF

            DEALLOCATE (XG)

         END IF

      END IF

      X (1) = XA  ! Make sure of it
      X (N) = XB

C     Formats:

CC 90 FORMAT (///, ' EXPDIS5:  Bad return from ZERORC. ISTAT = ', I2,
CC   >        /, ' XA, XB, DX =', 1P, 3E15.6,
CC   >        /, ' N = ', I4, '   MODE = ', I1,
CC   >        /, ' Proceeding with geometric distribution.')

      END SUBROUTINE EXPDIS5
