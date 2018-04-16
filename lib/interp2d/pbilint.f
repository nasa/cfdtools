C+------------------------------------------------------------------------------
C
      SUBROUTINE PBILINT (MODE, IDIM, JDIM, I1, I2, J1, J2, XS, YS, ZS,
     >                    US, VS, UTARG, VTARG, IS, JS, EPS, P, Q,
     >                    XINT, YINT, ZINT, XYZU, XYZV, IER)
C
C     Parametric bilinear interpolation within a regular 3-space surface mesh.
C     Arguments are much as for PLBICUBE, q.v.  New arguments P, Q are outputs.
C
C     Jan. 95   DAS  Initial implementation (no derivatives).
C     08/06/99   "   Analytic derivatives can come from the BILINT formulation.
C                    Also, return the relevant (p,q) that is sometimes needed.
C
C     Author:  David Saunders  Sterling Software/NASA Ames, Mt. View, CA.
C
C-------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) ::
     >   MODE, IDIM, JDIM, I1, I2, J1, J2

      REAL, INTENT (IN), DIMENSION (IDIM, JDIM) ::
     >   XS, YS, ZS, US, VS 

      REAL, INTENT (IN) ::
     >   UTARG, VTARG, EPS

      INTEGER, INTENT (INOUT) ::
     >   IS, JS

      REAL, INTENT (OUT) ::
     >   P, Q, XINT, YINT, ZINT, XYZU (3), XYZV (3)

      INTEGER, INTENT (OUT) ::
     >   IER

C     Local constants:

      REAL, PARAMETER ::
     >   ONE = 1.

C     Local variables:

      REAL
     >   A (2, 2), AJ (2, 2), DP (2), FP, FQ, FPQ, PM1, QM1

C     Procedures:

      EXTERNAL
     >   RIPPLE2D

C     Execution:
      
C     Locate the enclosing cell:

      CALL RIPPLE2D (IDIM, JDIM, I1, I2, J1, J2, US, VS,
     >               UTARG, VTARG, IS, JS, EPS, P, Q, IER)

C     Ignore IER - the nearest cell is always returned.

      PM1  = ONE - P
      QM1  = ONE - Q

      XINT = QM1 * (PM1 * XS (IS, JS)   + P * XS (IS+1, JS)) +
     >         Q * (PM1 * XS (IS, JS+1) + P * XS (IS+1, JS+1))
      YINT = QM1 * (PM1 * YS (IS, JS)   + P * YS (IS+1, JS)) +
     >         Q * (PM1 * YS (IS, JS+1) + P * YS (IS+1, JS+1))
      ZINT = QM1 * (PM1 * ZS (IS, JS)   + P * ZS (IS+1, JS)) +
     >         Q * (PM1 * ZS (IS, JS+1) + P * ZS (IS+1, JS+1))

C     Interpolated derivatives as well?

      IF (MODE /= 0) THEN

C        The BILINT formulation provides analytic partial dx/dp, etc.
C        The chain rule converts p/q derivatives to u/v derivatives.

C        Set up the 2 x 2 LHS matrix common to derivatives of x, y, and z.
C        The algebra saves a few operations over the more obvious form.

         FP  = US (IS + 1, JS) - US (IS, JS)
         FQ  = US (IS, JS + 1) - US (IS, JS)
         FPQ = US (IS + 1, JS + 1) - US (IS, JS + 1) - FP

         AJ (1, 1) = FP + Q * FPQ ! du/dp
         AJ (2, 1) = FQ + P * FPQ ! du/dq

         FP  = VS (IS + 1, JS) - VS (IS, JS)
         FQ  = VS (IS, JS + 1) - VS (IS, JS)
         FPQ = VS (IS + 1, JS + 1) - VS (IS, JS + 1) - FP

         AJ (1, 2) = FP + Q * FPQ ! dv/dp
         AJ (2, 2) = FQ + P * FPQ ! dv/dq

C        For X ...

         FP  = XS (IS + 1, JS) - XS (IS, JS)
         FQ  = XS (IS, JS + 1) - XS (IS, JS)
         FPQ = XS (IS + 1, JS + 1) - XS (IS, JS + 1) - FP

         DP (1) = FP + Q * FPQ ! dx/dp
         DP (2) = FQ + P * FPQ ! dx/dq

         A = AJ ! Could factorize once & solve 3 times, but LUSOLVE
                ! is in use elsewhere, while DECOMP & SOLVE are not.

         CALL LUSOLVE (2, 2, A, DP, IER) ! Solve J dfuv = dfpq

         XYZU (1) = DP (1) ! dx/du
         XYZV (1) = DP (2) ! dx/dv

C        For Y ...

         FP  = YS (IS + 1, JS) - YS (IS, JS)
         FQ  = YS (IS, JS + 1) - YS (IS, JS)
         FPQ = YS (IS + 1, JS + 1) - YS (IS, JS + 1) - FP

         DP (1) = FP + Q * FPQ ! dy/dp
         DP (2) = FQ + P * FPQ ! dy/dq

         A = AJ

         CALL LUSOLVE (2, 2, A, DP, IER)

         XYZU (2) = DP (1) ! dy/du
         XYZV (2) = DP (2) ! dy/dv

C        For Z ...

         FP  = ZS (IS + 1, JS) - ZS (IS, JS)
         FQ  = ZS (IS, JS + 1) - ZS (IS, JS)
         FPQ = ZS (IS + 1, JS + 1) - ZS (IS, JS + 1) - FP

         DP (1) = FP + Q * FPQ ! dz/dp
         DP (2) = FQ + P * FPQ ! dz/dq

         A = AJ

         CALL LUSOLVE (2, 2, A, DP, IER)

         XYZU (3) = DP (1) ! dz/du
         XYZV (3) = DP (2) ! dz/dv

      END IF

      END SUBROUTINE PBILINT
