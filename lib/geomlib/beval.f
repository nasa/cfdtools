C+----------------------------------------------------------------------
C
      SUBROUTINE BEVAL (BNAME, NP, P, ADD, NX, X, FX)
C
C PURPOSE:  BEVAL evaluates the indicated "bump" (shape function) at the
C           given abscissas, which are assumed to be normalized.  As the
C           bump names suggest,  these  functions  are commonly used for
C           perturbing airfoils.   They originated as the  "Hicks-Henne"
C           shape functions.  Some of them are also handy for generating
C           distributions of nonuniform weighting factors.
C
C           An option is provided to add the evaluations to the existing
C           values in the output array rather than just return values of
C           the current function.  This can save work-space in a calling
C           routine that is accumulating the effect of several bumps.
C
C METHOD:   The bump is selected by name rather than code number.   This
C           poses a problem of variable-length names,  resolved here  by
C           dealing with only the first four characters  -  an arbitrary
C           choice.   UPPER case is assumed  -  it seemed unnecessary to
C           handle lower case as well. Similarly, there is no attempt to
C           handle  unnormalized abscissas.
C
C           It was considered too inefficient, at this level, to attempt
C           handling of a given bump function's parameters by name. Thus
C           the ordering of the elements of P(*) IS important.
C
C           The option to ADD rather than just evaluate  is  implemented
C           by ALWAYS adding, but zeroing out first if necessary.   This
C           relieves the calling program from doing the zeroing.
C
C ARGUMENTS:
C    ARG    DIM  TYPE  I/O/S DESCRIPTION
C   BNAME    -   C*(*)   I   Name of desired  bump function.   Only
C                            the first 4 characters are  looked at,
C                            and they must be UPPER case.   Look at
C                            the code for valid names.
C    NP      -     I     I   Number of parameters (other than X) in
C                            the selected bump expression.
C    P      NP     R     I   The given values of the parameters, in
C                            a definite order.
C    ADD     -     L     I   ADD = .TRUE. means the evaluations for
C                            each X(I) are added into FX(I);
C                            ADD = .FALSE. means  FX(I)  is  zeroed
C                            out before this addition.
C    NX      -     I     I   The number of abscissas where the sel-
C                            ected function is to be evaluated.
C    X      NX     R     I   Abscissas in the range [0,1]. The case
C                            of simple scaling is an exception: the
C                            inputs and outputs involve ordinates.
C    FX     NX     R    I/O  The desired bump function values (pos-
C                            sibly added to input values; see ADD).
C
C NOTES:
C   *   The ordering of the shape function parameters is such that
C       the LAST one (P (NP)) is usually a multiplicative factor.
C
C HISTORY:
C   09/15/83   DAS   Initial design as FUNCTION BUMP.
C   01/14/84   DAS   Made simple scaling additive rather than multipli-
C                    cative: S*Y = Y + (S-1)*Y (simplifies usage).
C   02/22/84   DAS   Eliminated SQRT and SIN bump (redundant).
C   04/09/84   DAS   Adapted as SUBROUTINE BEVAL to get the loop inside
C                    rather than outside. Switched to selection by name
C                    rather than by type code; provided for  adding  as
C                    well as just evaluating; Wagner functions in-line.
C   07/17/86   DAS   Added simple "RAMP" function option.   Made Wagner
C                    functions the first choice.
C   08/11/86   DAS   Added "FLAP" and "SLAT" options (simple shearing).
C   07/26/89   RAK   ALOG changed to generic LOG.
C   11/06/93   RAK   "DROOP" function needed (1 - X) factor.
C   06/19/96   DAS   "EXP" function now has value 1. at specified X, as
C                    recommended by James Reuther;  added the symmetric
C                    forms of the modified sine function (SIN1 & SIN2),
C                    requiring use of leading 4 characters, not 3.
C   12/18/96   DAS   Added SINF, COSL, COSR, LCOS, & RCOS functions.
C   06/03/97   DAS   RADDEG name was misleading - changed it to DEGRAD.
C   10/20/99    "    Added SIN3 (fore/aft symmetry via SIN1 on [0,1])
C                    and   SIN4 (  "   "   "   "   "   "   " each half).
C
C AUTHORS: Leslie Collins, Robert Kennelly, David Saunders (Sterling);
C          Ray Hicks (NASA Ames Research Center)
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER
     >   NP, NX
      REAL
     >   FX (NX), P (NP), X (NX)
      LOGICAL
     >   ADD
      CHARACTER
     >   BNAME * (*)

C     Local constants:

      REAL, PARAMETER ::
     >   DEGRAD =  0.017453292519943E+0,
     >   PI     =  3.141592653589793E+0,
     >   PIBY2  =  1.570796326794897E+0,
     >   PT5LOG = -0.693147180559945E+0,
     >   RPI    =  0.318309886183790E+0,
     >   HALF   =  5.E-1,
     >   ONE    =  1.E+0,
     >   TWO    =  2.E+0,
     >   ZERO   =  0.E+0

C     Local variables:

      INTEGER
     >   I
      REAL
     >   AEXP, BEXP, CENTER, CENTER2, N, ONEMC, RN, POWER, POWERL,
     >   POWERR, SINXI, TANGNT, THETA, XI
      CHARACTER
     >   KEY * 4

C     Statement functions:

      REAL
     >   EBUMP, SBUMP, PWR, WIDTH, XNRM

      EBUMP (WIDTH, PWR, XNRM) =
     >   XNRM ** PWR * (ONE - XNRM) * EXP (-WIDTH * XNRM)

      SBUMP (WIDTH, PWR, XNRM) =
     >   (MAX (SIN (PI * XNRM ** PWR), ZERO)) ** WIDTH

C     Execution:

C     Check for just evaluating, not adding:

      IF (.NOT. ADD) THEN
         DO I = 1, NX
            FX (I) = ZERO
         END DO
      END IF

C     Avoid comparison of different-length strings:

      KEY = BNAME (1:4)

      IF (KEY == 'WAGN') THEN   ! Wagner functions.

C        Reference:  Ramamoorthy, P. and Padmavathi, K.  "Airfoil Design
C        by Optimization" in J. Aircraft, Vol. 14 (Feb. 1977), 219-221.

         N = P (1)
         RN = ONE / N

         IF (N == ONE) THEN
            DO I = 1, NX
               THETA = TWO * ASIN (SQRT (X (I)))
               FX (I) = FX (I) + ((THETA + SIN (THETA)) * RPI -
     >            (SIN (HALF * THETA)) ** 2) * P (2)
            END DO
         ELSE
            DO I = 1, NX
               THETA = TWO * ASIN (SQRT (X (I)))
               FX (I) = FX (I) + ((SIN (N * THETA) * RN  +
     >            SIN ((N-ONE) * THETA)) * RPI) * P (2)
            END DO
         END IF

      ELSE IF (KEY == 'SINE') THEN  ! Modified "SINE" (unsymmetric):
                                      ! P1 = "center", P2 = "width"
         POWER = PT5LOG / LOG (P (1))
         DO I = 1, NX
            FX (I) = P (3) * SBUMP (P (2), POWER, X (I)) + FX (I)
         END DO

      ELSE IF (KEY == 'SINF') THEN  ! Flipped "SINE" (unsymmetric):
                                      ! P1 = "center", P2 = "width"
         POWER = PT5LOG / LOG (ONE - P (1))
         DO I = 1, NX
            FX (I) = P (3) * SBUMP (P (2), POWER, ONE - X (I)) + FX (I)
         END DO

      ELSE IF (KEY == 'SIN1') THEN  ! Modified "SINE" (symmetric):
                                      ! P1 = "center", P2 = "width"
         CENTER = P (1)
         IF (CENTER <= HALF) THEN
            POWER = PT5LOG / LOG (CENTER)
            DO I = 1, NX
               FX (I) = P (3) * SBUMP (P (2), POWER, X (I)) + FX (I)
            END DO
         ELSE
            POWER = PT5LOG / LOG (ONE - CENTER)
            DO I = 1, NX
               FX (I) = P (3) * SBUMP (P (2), POWER, ONE - X (I)) +
     >                  FX (I)
            END DO
         END IF

      ELSE IF (KEY == 'SIN2') THEN  ! Modified "SINE" (symmetric):
                                      ! P1 = "center", P2 = "width"
         CENTER = P (1)
         IF (CENTER <= HALF) THEN
            POWER = PT5LOG / LOG (ONE - CENTER)
            DO I = 1, NX
               FX (I) = P (3) * SBUMP (P (2), POWER, ONE - X (I)) +
     >                  FX (I)
            END DO
         ELSE
            POWER = PT5LOG / LOG (CENTER)
            DO I = 1, NX
               FX (I) = P (3) * SBUMP (P (2), POWER, X (I)) + FX (I)
            END DO
         END IF

      ELSE IF (KEY == 'SIN3') THEN ! Fore/aft symmetry via SIN1 on [0, 1]

C        Assume [XA, XB] = [0, 1] and CENTER in (0, 0.5]:

         CENTER = P (1)
         POWER  = PT5LOG / LOG (CENTER)

         DO I = 1, NX
            XI = X (I)
            IF (XI > HALF) XI = ONE - XI
            FX (I) = P (3) * SBUMP (P (2), POWER, XI) + FX (I)
         END DO

      ELSE IF (KEY == 'SIN4') THEN ! Fore/aft symmetry via SIN1 on each half

C        Assume [XA, XB] = [0, 0.5] and CENTER in (0, 1).  Left half:

         CENTER  = P (1)
         POWERL  = PT5LOG / LOG (CENTER)
         CENTER2 = ONE - CENTER
         POWERR  = PT5LOG / LOG (CENTER2)

         IF (CENTER <= HALF) THEN
            DO I = 1, NX
               XI = MAX (ZERO, MIN (X (I) * TWO, ONE))
               FX (I) = P (3) * SBUMP (P (2), POWERL, XI) + FX (I)
            END DO
         ELSE
            DO I = 1, NX
               XI = MAX (ZERO, MIN (X (I) * TWO, ONE))
               FX (I) = P (3) * SBUMP (P (2), POWERR, ONE - XI) + FX (I)
            END DO
         END IF

C        Now the [0.5, 0.1] half with center <-- 1 - center

         IF (CENTER2 <= HALF) THEN
            DO I = 1, NX
               XI = MAX (ZERO, MIN ((X (I) - HALF) * TWO, ONE))
               FX (I) = P (3) * SBUMP (P (2), POWERR, XI) + FX (I)
            END DO
         ELSE
            DO I = 1, NX
               XI = MAX (ZERO, MIN ((X (I) - HALF) * TWO, ONE))
               FX (I) = P (3) * SBUMP (P (2), POWERL, ONE - XI) + FX (I)
            END DO
         END IF

      ELSE IF (KEY == 'COSL') THEN  ! 1/4 cosine (or sine), peak at left

C        SIN is used instead of COS because in the alternative SFEVAL form,
C        where COSL was first installed, PI is adjusted so that SIN (PI) > 0.

         DO I = 1, NX
            SINXI = MAX (SIN (PIBY2 * (X (I) + ONE)), ZERO)
            FX (I) = P (2) * SINXI ** P (1) + FX (I)
         END DO

      ELSE IF (KEY == 'COSR') THEN  ! 1/4 cosine (or sine), peak at right

C        SIN is used instead of COS for consistency with COSL.

         DO I = 1, NX
            SINXI = MAX (SIN (PIBY2 * X (I)), ZERO)
            FX (I) = P (2) * SINXI ** P (1) + FX (I)
         END DO

      ELSE IF (KEY == 'LCOS') THEN  ! Inverted 1/4 (co)sine, peak at left

         DO I = 1, NX
            SINXI = MAX (SIN (PIBY2 * X (I)), ZERO)
            FX (I) = P (2) * (ONE - SINXI ** P (1)) + FX (I)
         END DO

      ELSE IF (KEY == 'RCOS') THEN  ! Inverted 1/4 (co)sine, peak at right

         DO I = 1, NX
            SINXI = MAX (SIN (PIBY2 * (X (I) + ONE)), ZERO)
            FX (I) = P (2) * (ONE - SINXI ** P (1)) + FX (I)
         END DO

      ELSE IF (KEY == 'EXPO') THEN  ! "EXPONENTIAL" (peak height = 1.):
                                      ! P1 = "center", P2 = "width"
         ONEMC = ONE - P (1)
         AEXP  = P (1) * (ONE + ONEMC * P (2)) / ONEMC
         BEXP  = P (3) / EBUMP (P (2), AEXP, P (1))

         DO I = 1, NX
            FX (I) = BEXP * EBUMP (P (2), AEXP, X (I))  +  FX (I)
         END DO

      ELSE IF (KEY == 'DROO') THEN  ! "DROOP":  P1 = "width"

         DO I = 1, NX
            FX (I) = ((ONE - X (I)) * EXP (-P (1) * X (I))) * P (2)  +
     >               FX (I)
         END DO

      ELSE IF (KEY == 'LEAD') THEN  ! "LEADING"-edge:  P1 = "power"

         DO I = 1, NX
            FX (I) = ((ONE - X (I)) ** P (1)) * P (2)  +  FX (I)
         END DO

      ELSE IF (KEY == 'TRAI') THEN  ! "TRAILING"-edge:  P1 = "power"

         DO I = 1, NX
            FX (I) = (X (I) ** P (1)) * P (2)  +  FX (I)
         END DO

      ELSE IF (KEY == 'FLAP') THEN

C        "FLAP"-like function (shearing transformation only):
C        P (1) is the hinge point (fraction of chord from leading edge);
C        P (2) is the flap angle in DEGREEs (positive is "down").

         TANGNT = TAN (DEGRAD * P (2))
         DO I = 1, NX
            IF (X (I) > P (1))
     >         FX (I) = FX (I) - (X (I) - P (1)) * TANGNT
         END DO

      ELSE IF (KEY == 'SLAT') THEN

C        "SLAT"-like function (shearing transformation only):
C        P (1) is the hinge point (fraction of chord from leading edge);
C        P (2) is the flap angle in DEGREEs (positive is "down").

         TANGNT = TAN (DEGRAD * P (2))
         DO I = 1, NX
            IF (X (I) < P (1))
     >         FX (I) = FX (I) - (P (1) - X (I)) * TANGNT
         END DO

      ELSE IF (KEY == 'RAMP') THEN  ! "RAMP":  Y = P(1) * X

C        Commonly used in conjunction with Wagner functions.

         DO I = 1, NX
            FX (I) = P (1) * X (I)  +  FX (I)
         END DO

      ELSE IF (KEY == 'SCAL') THEN

C        Simple scaling - X(I) are probably ordinates.  Note that the
C        effect is arranged to be  additive,  not multiplicative,  so
C        that scaling can be treated just as for the  other functions
C        when it is included in a group for perturbing purposes.

         DO I = 1, NX
            FX (I) = X (I) * (P (1) - ONE)  +  FX (I)
         END DO

      ELSE
         WRITE (*, '(/, A)')  ' BEVAL: Illegal bump name.'
         STOP
      END IF

      END SUBROUTINE BEVAL
