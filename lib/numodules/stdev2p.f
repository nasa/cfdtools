C+----------------------------------------------------------------------
C
      SUBROUTINE STDEV2P ( N, X, STRIDE, MEAN, STDEV )
C
C ACRONYM: STandard DEViation for real data; 2-Pass/2-Precision method
C          --       ---                      - -    - -
C
C PURPOSE: STDEV2P returns the mean and standard deviation of N elements
C          of a real array X(*).  The elements must be some constant 
C          stride apart.  Element 1 is taken as the first to be included.
C
C          STDEV2P and STDEV1P (1-pass/1-precision) are recommended as 
C          the best of numerous possibilities that have been examined.
C          See STDEV1P for a complete discussion of various methods.  See 
C          also STDEVI2 for the case of INTEGER*2 data.
C
C USAGE:
C
C          Elements 1, 4, 7, ... CALL STDEV2P ( N, X, 3, XMEAN, XSD )
C
C          Elements I, I+M, I+2M, ... CALL STDEV2P ( N, X(I), M, XBAR, SD )
C          (not compatible with FORTRAN 8X)
C
C ARGUMENTS:
C    ARG     DIM   TYPE I/O/S DESCRIPTION 
C    N        -     I     I   Number of elements in sums involved; N > 0.
C    X        *     R     I   Elements 1, 1+STRIDE, 1+2*STRIDE, ... of
C                             X(*) are processed.
C    STRIDE   -     I     I   See description of argument X.
C    MEAN     -     R     O   Mean of of specified elements.
C    STDEV    -     R     O   Standard deviation, where SD**2 is taken to
C                             be the sum of [ ( X(I) - XMEAN )**2 ] / N and
C                             not   "     "     "      "      "     / N-1
C    
C METHOD:  STDEV2P implements the best known compromise assuming that internal
C          use of double precision is viable.  Its two passes - one for the
C          sum of the Xs, and one for the sum of squares of the deviations from
C          the mean - are efficient and accurate for any likely application.
C
C          See STDEV1P for a complete discussion of various methods and 
C          recommendations for their use.
C     
C ERROR HANDLING:  None.
C
C ENVIRONMENT:  Digital VAX/VMS FORTRAN
C               IMPLICIT NONE and the 7-character name are non-standard.
C
C DEVELOPMENT HISTORY:
C       DATE   INITIALS   DESCRIPTION 
C     05/25/88   CLH      Initial implementation based on testing of
C                         various methods (see STDEV1P).
C
C AUTHOR: Charles L. Hooper, Sterling Software, Palo Alto, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT  NONE

C --- Arguments:

      INTEGER   N, STRIDE
      REAL      X(*), MEAN, STDEV

C --- Local constants:

      DOUBLE PRECISION  ZERO
      PARAMETER ( ZERO = 0.0E+00 )

C --- Local variables:

      INTEGER   I
      DOUBLE PRECISION  AV2, MEAND

C --- Initialize accumulators:

      MEAND = ZERO 
      AV2   = ZERO

C --- Find the sum of the elements:

      DO 100 I = 1, 1+(N-1)*STRIDE, STRIDE
         MEAND = MEAND + DBLE ( X(I) )
  100 CONTINUE
     
C --- Compute the mean:

      MEAND = MEAND / DBLE ( N )

C --- Sum the squares of the deviations from the mean:

      DO 200 I = 1, 1+(N-1)*STRIDE, STRIDE
         AV2 = AV2 + ( DBLE ( X(I) ) - MEAND )**2 
  200 CONTINUE

      MEAN  = SNGL ( MEAND )
      STDEV = SNGL ( SQRT ( MAX ( AV2 / DBLE ( N ), ZERO )))

      RETURN
      END
