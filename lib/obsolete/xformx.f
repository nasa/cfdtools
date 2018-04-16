C+----------------------------------------------------------------------
C
      SUBROUTINE XFORMX ( NX, XOLD, ANEW, BNEW, XNEW )
C
C PURPOSE:
C     XFORMX is intended to apply the general linear transformation
C
C                        x' := m x + c
C
C     where m, c are such that interval [a, b] is transformed to [p, q].
C
C     In this implementation (prompted by a need to transform arbitrary
C     units to the interval [0, 2*pi]), the elements of vector XOLD (*)
C     are transformed to XNEW (*), and [a, b] is assumed to be defined
C     by [XOLD (1), XOLD (N)], thereby eliminating two arguments.
C
C ARGUMENTS:
C   ARG     DIM    TYPE  I/O/S DESCRIPTION 
C   NX       -       I     I   Number of elements to be transformed (>1)
C   XOLD     NX      R     I   Array to be transformed to [ANEW, BNEW].
C                              XOLD (1) and XOLD (N) define [a, b].
C   ANEW     -       R     I   Bounds of desired transformed interval
C   BNEW     -       R     I
C   XNEW     NX      R     O   Transformed array - may be same locations
C                              as XOLD (*)
C
C AUTHOR:  David Saunders, Sterling Software, Palo Alto, CA., Dec. 1986
C
C-----------------------------------------------------------------------

      INTEGER   NX
      REAL      ANEW, BNEW, XNEW (NX), XOLD (NX)

      INTEGER   I
      REAL      C1, C2

C     Ensure that the end-points are transformed exactly (no round-off).

      C2 = XOLD (NX) - XOLD (1)
      C1 = ( BNEW - ANEW ) / C2
      C2 = ( XOLD (NX) * ANEW - XOLD (1) * BNEW ) / C2

      XNEW (1) = ANEW
      DO 20, I = 2, NX - 1
         XNEW (I) = C1 * XOLD (I) + C2
   20 CONTINUE
      XNEW (NX) = BNEW

      RETURN
      END
