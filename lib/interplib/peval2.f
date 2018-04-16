C+----------------------------------------------------------------------
C
      SUBROUTINE PEVAL2 (MODE, NT, T, NGEOM, XGEOM, YGEOM, TGEOM,
     +                   PSCOFS, X, Y)
C
C  ACRONYM:    Parametric spline EVALuation (2-D) for given T(s)
C              -                 ----        -
C  PURPOSE:    PEVAL2 evaluates the 2-D parametric spline previously
C              fitted by PFIT2.  For the given value(s) of parameter
C              T, it returns values of X and Y.  T(*) may be entered
C              as a normalized distribution, depending on MODE.
C
C  ARGUMENTS:
C   ARG    TYPE  I/O/S   DIM   DESCRIPTION
C   MODE    I      I      -    0 means T(*) is not normalized;
C                              1 means a normalized or relative
C                              distribution is supplied, so evaluations
C                              are at T = T(I) * TGEOM(NGEOM).
C                              For storage efficiency, T(*) is scaled
C                              in-place then renormalized.
C   NT      I      I      -    Number of evaluations required.  NT>=1.
C   T       R      I      NT   Values of the parametric variable,
C                              possibly normalized, for which X, Y values
C                              are desired.  See MODE.
C   NGEOM   R      I      -    Number of data points fitted by PFIT2.
C   XGEOM,  R      I    NGEOM  The data point coordinates.
C    YGEOM
C   TGEOM   R      I    NGEOM  Corresponding parametric variable values
C                              from PFIT2.
C   PSCOFS  R      I   NGEOM*6 Spline coefs. calculated by PFIT2.
C   X, Y    R      O      NT   Evaluated coordinates.
C
C  PROCEDURES:
C     CSEVAL   Evaluates a 1-dimensional cubic spline
C     SSCAL    BLAS utility
C
C  ENVIRONMENT:  FORTRAN 77
C
C  HISTORY:
C     03/29/94  D.A.Saunders  Adapted from PEVAL3.
C
C  AUTHOR:  David Saunders, Sterling Software/NASA Ames, Mt. View, CA.
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER    MODE, NT, NGEOM
      REAL       XGEOM (NGEOM), YGEOM (NGEOM), TGEOM (NGEOM),
     >           PSCOFS (NGEOM, 3, 2), T (NT), X (NT), Y (NT)

C     Local constants:

      REAL       ONE
      PARAMETER (ONE = 1.E+0)

C     Procedures:

      EXTERNAL   CSEVAL, SSCAL

C     Execution:

C     Is T (*) normalized?

      IF (MODE .GT. 0) CALL SSCAL (NT, TGEOM (NGEOM), T, 1)

C     Evaluate spline XGEOM vs. TGEOM at all values of T; repeat for Y.

      CALL CSEVAL (NGEOM, TGEOM, XGEOM, NT, T, PSCOFS (1, 1, 1),
     >             PSCOFS (1, 2, 1), PSCOFS (1, 3, 1), X)

      CALL CSEVAL (NGEOM, TGEOM, YGEOM, NT, T, PSCOFS (1, 1, 2),
     >             PSCOFS (1, 2, 2), PSCOFS (1, 3, 2), Y)

      IF (MODE .GT. 0) THEN
         CALL SSCAL (NT, ONE / TGEOM (NGEOM), T, 1)
         T (NT) = ONE
      END IF

      RETURN
      END
