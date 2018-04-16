C+----------------------------------------------------------------------
C
      SUBROUTINE SUMX (N, X, ISTEP, SUM)
C
C PURPOSE: SUMX calculates the sum of N elements of a single precision
C          array X(*).  The elements must be some constant step apart.
C          Element 1 is taken as the first to be included.
C
C METHOD:  A local double precision accumulator is used prior to return-
C          ing the result in single precision form.  No attempt is made
C          to provide for proper summing of elements in subsets via more
C          than one call to SUMX.  (IMSL's extended precision routines
C          (Chapter V) do provide such capability but require one call
C          per element of the sum.)
C
C    Usage:
C
C          Elements 1, 4, 7, ... :  CALL SUMX (N, X, 3, SUM1)
C
C          Elements 2, 4, 6, ... :  CALL SUMX (N, X (2), 2, SUM2)
C
C PARAMETERS:
C   ARG     DIM    TYPE  I/O/S DESCRIPTION 
C   N        -       I     I   Number of elements to sum
C   X        *       R     I   Elements 1, 1+ISTEP, 1+2*ISTEP, ...
C   ISTEP    -       I     I   of X (*) are summed
C   SUM      -       R     O   Sum of specified elements
C
C ENVIRONMENT:  VAX/VMS, FORTRAN 77.
C
C HISTORY:
C   08/19/87  DAS  Adapted from MEANX, partly to avoid MEANX's recursive
C                  approach to preserving precision, which involves N
C                  divides.  But MEANX could be translated to double
C                  precision where SUMX is really strictly for single
C                  precision data.
C
C AUTHOR:  David Saunders, Sterling Software, Palo Alto, CA.
C
C-----------------------------------------------------------------------

C     Arguments:

      INTEGER
     >   ISTEP, N
      REAL
     >   X (*), SUM
      
C     Local variables:

      INTEGER
     >   I
      DOUBLE PRECISION
     >   ACCUM

C     Execution:

      ACCUM = 0.D+0
      DO 20 I = 1, 1 + (N - 1) * ISTEP, ISTEP
         ACCUM = ACCUM + DBLE (X (I))
   20 CONTINUE

      SUM = REAL (ACCUM)

      RETURN
      END
