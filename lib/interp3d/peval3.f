C+---------------------------------------------------------------------
C
      SUBROUTINE PEVAL3 ( NT, T, NGEOM, XGEOM, YGEOM, ZGEOM, TGEOM,
     +                    PSCOFS, X, Y, Z )
C
C  ACRONYM:    Parametric spline EVALuation (3-D) for given T(s)
C              -                 ----        -
C  PURPOSE:    PEVAL3 evaluates the 3-D parametric spline previously
C              fitted by PFIT3.  For the given value(s) of parameter
C              T, it returns values of X, Y, and Z.
C
C  ARGUMENTS:
C   ARG    TYPE  I/O/S   DIM   DESCRIPTION
C   NT      I      I      -    Number of evaluations required.  NT>=1.
C   T       R      I      NT   Arbitrary values of parametric variable,
C                              for which X, Y, Z values are desired.
C   NGEOM   R      I      -    Number of data points fitted by PFIT3.
C   XGEOM,  R      I    NGEOM  Discrete points defining curve fitted
C    YGEOM,                    by PFIT3.
C    ZGEOM
C   TGEOM   R      I    NGEOM  Values of parametric variable from PFIT3.
C   PSCOFS  R      I   NGEOM*9 Spline coefs. calculated by PFIT3.
C   X,      R      O      NT   Evaluated coordinates.
C    Y,
C    Z
C
C  ERROR HANDLING: None.
C
C  EXTERNAL REFERENCES:
C   CSEVAL - Evaluates a 1-dimensional cubic spline
C
C  ENVIRONMENT:  VAX/VMS FORTRAN
C
C  AUTHOR:   David Saunders, Sterling Software, Palo Alto, CA.
C
C  DEVELOPMENT HISTORY:
C    DATE   INITIALS  DESCRIPTION
C  09/18/86   DAS     Adapted from PSTVAL. No attempt at derivatives.
C
C----------------------------------------------------------------------

      IMPLICIT  NONE

C  *  Arguments:

      INTEGER   NT, NGEOM
      REAL      XGEOM(NGEOM), YGEOM(NGEOM), ZGEOM(NGEOM), TGEOM(NGEOM),
     +          PSCOFS(NGEOM,3,3), T(NT), X(NT), Y(NT), Z(NT)

C  *  Procedures:

      EXTERNAL  CSEVAL

C  *  Execution:

C  *  Evaluate spline XGEOM vs. TGEOM at all values of T; repeat for Y, Z:

      CALL CSEVAL ( NGEOM, TGEOM, XGEOM, NT, T, PSCOFS(1,1,1),
     +              PSCOFS(1,2,1), PSCOFS(1,3,1), X )

      CALL CSEVAL ( NGEOM, TGEOM, YGEOM, NT, T, PSCOFS(1,1,2),
     +              PSCOFS(1,2,2), PSCOFS(1,3,2), Y )

      CALL CSEVAL ( NGEOM, TGEOM, ZGEOM, NT, T, PSCOFS(1,1,3),
     +              PSCOFS(1,2,3), PSCOFS(1,3,3), Z )

      RETURN
      END
