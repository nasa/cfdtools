C+----------------------------------------------------------------------
C
      SUBROUTINE MEANX ( N, X, ISTEP, XMEAN )
C
C PURPOSE: MEANX calculates the average of N elements of a real array.
C          The elements must be some constant step >=1 apart.  Element 1
C          is taken as the first to be included.
C
C METHOD:  To avoid accumulating large sums, use the recurrence relation
C
C               XMEAN(I) = XMEAN(I-1) + ( X(I) - XMEAN(I-1) ) / I
C
C    Usage:
C
C          Elements 1, 4, 7, ... :  CALL MEANX ( N, X, 3, XMEAN )
C
C          Elements 2, 4, 6, ... :  CALL MEANX ( N, X(2), 2, XMEAN )
C
C PARAMETERS:
C   ARG     DIM    TYPE  I/O/S DESCRIPTION 
C   N        -       I     I   Number of elements to average
C   X        *       R     I   Elements 1, 1+ISTEP, 1+2*ISTEP, ...
C   ISTEP    -       I     I   of X(*) are averaged
C   XMEAN    -       R     O   Average of specified elements
C
C ENVIRONMENT:  VAX/VMS, FORTRAN 77.
C
C AUTHOR:  David Saunders, Sterling Software, Palo Alto, CA.
C          June 1986, adapted from version for integer data.
C
C-----------------------------------------------------------------------

C ... Arguments:

      INTEGER   ISTEP, N
      REAL      X(*), XMEAN
      
C ... Local variables:

      INTEGER   I
      REAL      AVG

C ... Intrinsics:

      INTRINSIC REAL

C ... Execution:

      AVG = 0.E+0
      DO 20 I = 1, 1 + ( N-1 ) * ISTEP, ISTEP
         AVG = AVG  +  ( X(I) - AVG ) / REAL ( I )
   20 CONTINUE

      XMEAN = AVG
      RETURN
      END
