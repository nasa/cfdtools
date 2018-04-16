C+----------------------------------------------------------------------
C
      SUBROUTINE FD12K (N, X, F, FP, FPP, FK)
C
C  PURPOSE: FD12K returns estimates of 1st and 2nd derivatives and of
C           curvature, by finite differencing, for each of the points
C           (X(I), F(I)), I = 1 : N.  The abscissas are assumed to be
C           nonuniform, and they must be monotonic.
C
C           This routine combines calls to FDCNTR, FD1SID, FDCURV for
C           the common case of needing results for N >= 2 points.
C
C           If (say) curvature is wanted at a single point only, call
C           either FDCNTR or FD1SID and FDCURV directly.
C
C  INPUTS:  X(*) & F(*) are N-vectors defining some curve in 2-space.
C           For N > 2, the 3-pt formulas are used for all I (with the
C                      one-sided forms used at the end-points).
C           For N = 2, the 2-pt formulas are used.
C
C  OUTPUTS: FP, FPP, FK are N-vectors representing 1st and 2nd deriv-
C           atives and curvature respectively.  These are assigned in
C           reverse order (FK, FPP, FP) so that a call such as
C
C                     CALL FD12K (N, X, Y, YP, YP, YP)
C
C           can be used if just 1st derivatives are desired, to avoid
C           declaring storage for FPP and FK. (Similarly for the case
C           when 1st and 2nd derivatives are desired but curvature is
C           not. The unnecessary arithmetic in these cases is consid-
C           ered preferable to another argument and extra logic.)
C
C  METHOD:  Central differencing is used at all interior points, with
C           one-sided 3-point formulas used at each end-point.
C
C           The curvature formula is safeguarded against overflow  in
C           the presence of large 1st derivatives.  The basic formula
C           used here is:
C
C               FK (I) = FPP (I) / (1 + FP(I) ** 2) ** 3/2
C
C           Note that if X is not necessarily monotonic, curvature is
C           defined as
C
C               CURVATURE = (X" ** 2  +  Y" ** 2) ** 1/2
C
C           where " means 2nd derivative with respect to  arc-length.
C           See modules CURV2D and CURV3D for these parametric cases.
C
C  NOTES:   1. Finite differencing errors can be large if the delta-X
C              values are too small,  especially if the precision  in
C              the function values is less than full.
C           2. Nevertheless, finite differences have been observed to
C              behave better than the analytic derivatives of splines
C              in airfoil geometry applications.
C
C  EXTERNALS:
C           FDCNTR modularizes the central 3-point formulas for first
C                  and second derivatives by finite differences.
C           FDCURV modularizes the curvature formula (safe-guarded).
C           FD1SID modularizes the 1-sided forward and backward 3-pt.
C                  formulas for first and second derivatives.
C
C  HISTORY:
C           09/15/83   DAS   Initial implementation (interior pts only).
C           12/27/85   DAS   End points are handled by FD1SID now.
C           09/18/87   DAS   The N=2 case is handled now.
C           08/21/89   DAS   Formulation revised to use separate dF/dX
C                            terms instead of a common denominator.
C           08/17/91   DAS   Introduced FDCNTR and FDCURV when it was
C                            found that FD12K did not lend itself to
C                            application to one point at a time.
C
C  AUTHOR:  David Saunders, Sterling Software/NASA Ames, Palo Alto, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT   NONE

C     Arguments:

      INTEGER    N
      REAL       X (N), F (N), FP (N), FPP (N), FK (N)

C     Local constants:

      REAL       ZERO
      PARAMETER (ZERO = 0.E+0)

C     Local variables:

      INTEGER    I
      REAL       FPI, FPPI

C     Procedures:

      EXTERNAL   FDCNTR, FDCURV, FD1SID

C     Execution:

C     Assign values in the order  curvature, f", f'  so that the
C     user can pass the same array for, say, curvature and f" if
C     curvature is not wanted:

      IF (N .EQ. 2) THEN

         FK (1)  = ZERO
         FPP (1) = ZERO
         FP (1)  = (F (2) - F (1)) / (X (2) - X (1))
         FK (2)  = ZERO
         FPP (2) = ZERO
         FP (2)  = FP (1)

      ELSE

C        Forward 3-pt. differencing for the first point:

         CALL FD1SID (1, 1, X, F, FPI, FPPI)
         CALL FDCURV (FPI, FPPI, FK (1))
         FPP (1) = FPPI
         FP (1)  = FPI
      
C        Central 3-pt. differencing for the bulk of the points:

         DO 20 I = 2, N-1
            CALL FDCNTR (I, X, F, FPI, FPPI)
            CALL FDCURV (FPI, FPPI, FK (I))
            FPP (I) = FPPI
            FP (I)  = FPI
   20    CONTINUE

C        Backward differencing for the last point:

         CALL FD1SID (N, -1, X, F, FPI, FPPI)
         CALL FDCURV (FPI, FPPI, FK (N))
         FPP (N) = FPPI
         FP (N)  = FPI
      END IF

      RETURN
      END
