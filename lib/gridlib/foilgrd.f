C+----------------------------------------------------------------------
C
      SUBROUTINE FOILGRD (N, XMIN, XMAX, WTLIN, WTQUAD, WTSIN, WTCOS, X)
C
C  ONE-LINER:  1-D grid distribution tailored for airfoils
C
C  DESCRIPTION:
C
C        FOILGRD generates N abscissas in the range [XMIN, XMAX] using a
C     combination of linear, quadratic, sine, and cosine shape functions.
C     The four weights should add to 1. and each should be in [0., 1.].
C     Values suggested for a typical airfoil:  0.04, 0.0, 0.3, 0.66
C
C  ARGUMENTS:
C    ARG    TYPE  I/O/S  DIM   DESCRIPTION
C    N        I     I     -    |N| = desired number of abscissas;
C                              N < 0 can be used to reverse the bunching
C                              of the N > 0 case, to suit the lower
C                              surface of a wrap-around airfoil
C  XMIN,XMAX  R     I     -    First and last values to include in X (*)
C   WTLIN     R     I     -    Weight applied to the linear component
C   WTQUAD    R     I     -    Weight applied to the quadratic term
C   WTSIN     R     I     -    Weight applied to the sinusoidal term
C   WTCOS     R     I     -    Weight applied to the cosine term
C    X        R     O     N    Desired distribution of abscissas;
C                              X (1) = XMIN and X (|N|) = XMAX
C
C  08/23/95  JJR  Adaptation of FOILGRID (sine + quadratic).
C  08/29/95  DAS  Tidied up the nomenclature.
C  07/21/97   "   Implemented N < 0 option to simplify usage.
C
C  AUTHOR:  James Reuther, Ames Research Center, Mt. View, CA
C
C----------------------------------------------------------------------

C     Arguments:

      INTEGER    N
      REAL       XMIN, XMAX, WTLIN, WTQUAD, WTSIN, WTCOS, X (*)

C     Local constants:

      REAL       HALF, ONE, PI

      PARAMETER (HALF = 0.5,
     >           ONE  = 1.0,
     >           PI   = 180.)

C     Local variables:

      INTEGER    I, J, NPTS
      REAL       RANGE, RN, SNORM, WL, WQ, WS, WC, XNORM
      REAL       TOTAL

C     Intrinsics:

      REAL       COSD

C     Execution.

      TOTAL = ONE / (WTLIN + WTQUAD + WTSIN + WTCOS)
      WL    = WTLIN * TOTAL
      WQ    = WTQUAD* TOTAL
      WS    = WTSIN * TOTAL
      WC    = WTCOS * TOTAL
      RANGE = XMAX - XMIN
      NPTS  = ABS (N)
      RN    = ONE / REAL (NPTS - 1)
      X (1) = XMIN

      IF (N .GT. 0) THEN  ! Ascending order

         DO I = 2, NPTS - 1
            XNORM = REAL (I - 1) * RN
            SNORM = WL * XNORM +
     >              WQ * XNORM ** 2 +
     >              WS * (ONE - COSD (XNORM * PI * HALF)) +
     >              WC * (ONE - COSD (XNORM * PI)) * HALF
            X (I) = XMIN + SNORM * RANGE
         END DO

      ELSE  ! N < 0 - reverse the order

         J = NPTS - 1
         DO I = 2, NPTS - 1
            J = J - 1
            XNORM = ONE - REAL (J) * RN
            SNORM = WL * XNORM +
     >              WQ * XNORM ** 2 +
     >              WS * (ONE - COSD (XNORM * PI * HALF)) +
     >              WC * (ONE - COSD (XNORM * PI)) * HALF
            X (I) = XMIN + SNORM * RANGE
         END DO

      END IF

      X (NPTS) = XMAX

      RETURN
      END
