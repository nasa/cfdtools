C+------------------------------------------------------------------------------
C
      SUBROUTINE PLINCUB (IDIM, JDIM, I1, I2, J1, J2, XS, YS, ZS,
     >                    US, VS, UTARG, VTARG, LINUV, CMETHOD, CLOSED,
     >                    IS, JS, EPS, XINT, YINT, ZINT, IER)
C
C     PLINCUB performs parametric interpolation within a regular 3-space surface
C     mesh for the case where linear interpolation is required in one of the
C     (u,v) directions and cubic in the other, as is typical for wing surface
C     grids.  The (u,v) cells should be ~parallelograms for best results.
C
C     LINUV   = 1 or 2 means the 1-D cubics are in the V or U directions resp.;
C     CMETHOD = 'M' or 'B' for monotonic ("tight") or "loose" ("Bessel")
C               controls the 1-D cubic results from LCSFIT.
C
C     The other arguments are as for PLBICUBE and PBILINT, q.v.
C
C     12/17/97   DAS   Adaptation of PBILINT for wing surfaces (Fortran 90).
C     01/08/98   DAS   Removed a LINXYZ argument - strictly 1-D methods now.
C
C     Author:  David Saunders  Sterling Software/NASA Ames, Moffett Field, CA.
C
C-------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER   IDIM, JDIM, I1, I2, J1, J2, IS, JS, LINUV, IER
      REAL      XS (IDIM, JDIM), YS (IDIM, JDIM), ZS (IDIM, JDIM),
     >          US (IDIM, JDIM), VS (IDIM, JDIM),
     >          UTARG, VTARG, EPS, XINT, YINT, ZINT
      CHARACTER CMETHOD * 1
      LOGICAL   CLOSED

C     Local constants:

      REAL,    PARAMETER :: DERIVS = -999., ONE = 1.
      LOGICAL, PARAMETER :: NEW = .TRUE.

C     Local variables:

      INTEGER   I, IA, IB, IEVAL, J, JA, JB, N
      REAL      X(4), Y(4), Z(4), T(4), P, PM1, Q, QM1, TEVAL

C     Procedures:

      EXTERNAL  PLSCRV3D, RIPPLE2D

C     Execution:
      
C     Locate the enclosing cell:

      CALL RIPPLE2D (IDIM, JDIM, I1, I2, J1, J2, US, VS,
     >               UTARG, VTARG, IS, JS, EPS, P, Q, IER)

C     Ignore IER - the nearest cell is always returned.

      IF (LINUV .EQ. 2) THEN  ! Linear in the V (J) direction

         QM1 = ONE - Q
         IA  = MAX (I1, IS - 1)  ! PLSCRV3D uses a 4-pt. method
         IB  = MIN (I2, IS + 2)
         N   = 0

         DO I = IA, IB
            N = N + 1
            X(N) = QM1 * XS(I,JS) + Q * XS(I,JS+1)
            Y(N) = QM1 * YS(I,JS) + Q * YS(I,JS+1)
            Z(N) = QM1 * ZS(I,JS) + Q * ZS(I,JS+1)
            T(N) = QM1 * US(I,JS) + Q * US(I,JS+1)
         END DO

         IEVAL = IS - IA + 1
         TEVAL = UTARG

      ELSE  ! LINUV = 1 meaning 1-D cubics in the 2nd (V/J) direction

         PM1 = ONE - P
         JA  = MAX (J1, JS - 1)
         JB  = MIN (J2, JS + 2)
         N   = 0

         DO J = JA, JB
            N = N + 1
            X(N) = PM1 * XS(IS,J) + P * XS(IS+1,J)
            Y(N) = PM1 * YS(IS,J) + P * YS(IS+1,J)
            Z(N) = PM1 * ZS(IS,J) + P * ZS(IS+1,J)
            T(N) = PM1 * VS(IS,J) + P * VS(IS+1,J)
         END DO

         IEVAL = JS - JA + 1
         TEVAL = VTARG

      END IF

      CALL PLSCRV3D (N, X, Y, Z, T, CMETHOD, NEW, CLOSED, TEVAL, IEVAL,
     >               XINT, YINT, ZINT, DERIVS)

      RETURN
      END
