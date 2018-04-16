C+----------------------------------------------------------------------
C
      SUBROUTINE STDEV1P (N, X, STRIDE, MEAN, SD)
C
C ACRONYM: STandard DEViation for real data; 1-Pass/1-Precision method
C          --       ---                      - -    - -
C PURPOSE:
C        STDEV1P returns the mean and standard deviation of N elements
C     of a real array X(*).  The elements must be some constant stride
C     apart.  Element 1 is taken as the first to be included.
C
C        STDEV1P and STDEV2P are recommended as the best of numerous
C     possibilities that have been examined, as discussed below.  See
C     also STDEVI2 for the case of INTEGER*2 data.
C
C     Usage:
C
C        Elements 1, 4, 7, ...   CALL STDEV1P (N, X, 3, XMEAN, XSD)
C
C        Elements I, I+M, I+2M, ...  "   "  (N, X (I), M, XBAR, SD)
C        (as for STDEVI2, but not compatible with FORTRAN 8X, however...)
C
C ARGUMENTS:
C   ARG     DIM  TYPE  I/O/S  DESCRIPTION 
C   N        -      I    I    Number of elements in sums involved; N > 0
C   X        *      R    I    Elements 1, 1+STRIDE, 1+2*STRIDE, ... of
C   STRIDE   -      I    I       X (*) are processed
C   MEAN     -      R    O    Mean of specified elements
C   SD       -      R    O    Standard deviation, where SD**2 is taken to
C                                be  ( <Sum> (X (I) - MEAN)**2 ) / N
C                                not     "    "    "    "    "   / (N-1)
C
C METHOD:
C        STDEV1P implements the best known method for computing standard
C     deviation (and mean) given the restrictions of using a single pass
C     and not using any double precision.  It proves to be about 30% slower
C     than the more obvious two-pass/double-precision-accumulation method
C     of STDEV2P, but has the virtue of permitting a double-precision analog
C     if such is needed (where quadruple precision would be needed by other
C     methods).
C     
C        Since the best compromise between speed and accuracy has been unclear
C     in the past, the various possibilites are now summarized:
C
C     One-pass methods:
C
C     (1) Sum the elements (SUM1) and the squares of the elements (SUM2)
C         in the obvious way, and use
C            MEAN  = SUM1 / N
C            SD**2 = SUM2 / N - MEAN**2
C         This involves the least arithmetic but can lead to large sums
C         and resulting truncation errors.  It is not recommended even if
C         double precision accumulators are used (except on I*2 data).
C
C     (2) Use "running averages" for the means of both X and X**2:
C            MEAN1 (I) = MEAN1 (I-1) + ( X (I) - MEAN1 (I-1) ) / I    and
C            MEAN2 (I) = MEAN2 (I-1) + ( X (I)**2 - MEAN2 (I-1) ) / I  (I=1:N)
C         Then   SD**2 = MEAN2 - MEAN1**2.
C         This approach avoids the SUM of squares, but still involves squaring
C         the data and requires safeguarding of the final square root.  It may
C         suit applications involving MANY means and std. devs., where double
C         precision use of (1) would double the storage over (2) in single.
C
C     (3) Use the above recurrence for the mean, and the following recurrence
C         derived from the form  N * SD**2 = <sum> ( X - MEAN )**2
C         for the standard deviation:
C            SSQ (I) = SSQ (I-1) + I*(I-1)*(( X (I) - MEAN (I-1)) / I )**2
C         which is certainly non-negative.  Then SD = SQRT (SSQ (N) / N).
C         This avoids squares of the data (though the mean may be zero), at
C         some cost in arithmetic.  This is the method of STDEV1P, using
C         one-precision arithmetic throughout.
C         
C     Two-pass methods:
C
C     (1) Compute MEAN by summing the Xs using double precision accumulation.
C         Compute the standard deviation of (X - MEAN) in a second pass, again
C         using double precision accumulation for the sum of squares.  This
C         is the straightforward method of STDEV2P.  It proves to be both fast
C         and accurate.  Note that if MEAN=0 then subtracting the mean cannot
C         help.  But it is safe to say that normally, subtracting the mean
C         DOES help, especially if the standard deviation is small compared
C         with the mean.
C
C     (2) Other two-pass methods, such as simple summation for MEAN and the
C         simple recursion for the sum of (X - MEAN)**2, both using double
C         precision as far as possible, prove to be merely slower with no
C         significant gain in accuracy except possibly for extreme cases of
C         huge data ranges (standard deviations) and N in the millions.
C
C     Recommendations (following exhaustive testing by Charles Hooper):
C
C     STDEV2P is both accurate and efficient if internal use of double precision
C             is viable.  Two passes are affected by paging on virtual memory
C             machines sooner (for smaller N) than is a single pass (but this
C             method still wins by about 10%).
C     STDEV1P is accurate for any reasonable application, 30% less efficient
C             than STDEV2P, and retained because it lends itself to a double-
C             precision-data analog.
C
C ENVIRONMENT:  VAX/VMS, FORTRAN 77
C               IMPLICIT NONE and the 7-character name are nonstandard
C
C HISTORY:
C   04/13/88  DAS  Initial implementation of method (3) above, in preference
C                  to earlier implementation of one-pass method (2) (SCSTDV),
C                  and in the hope of clarifying the possibilities.  See also
C                  STDEVI2 for an instance of one-pass method (1).
C   05/23/88  DAS  Documentation updated in the light of comparisons made by
C                  Charles Hooper (recommending STDEV2P and STDEV1P in that
C                  order for normal applications).
C   11/11/88  DAS  Updated details of methods (1), (2) and (3) descriptions.
C                  (Method (2) seems a fair compromise for updating whole
C                  arrays of means and standard deviations.)
C
C AUTHOR:  David Saunders, Sterling Software, Palo Alto, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT   NONE

C     Arguments:

      INTEGER    N, STRIDE
      REAL       X (*), MEAN, SD
      
C     Local variables:

      INTEGER    I 
      REAL       XBAR, SSQ, TERM

C     Constants:

      REAL       ZERO
      PARAMETER (ZERO = 0.E+0)

C     Execution:

      XBAR = ZERO
      SSQ  = ZERO

      DO 20, I = 1, 1 + (N - 1) * STRIDE, STRIDE
         TERM = (X (I) - XBAR) / REAL (I)            ! Permits just 1 divide
         XBAR = XBAR + TERM

C*****   SSQ  = SSQ  + TERM ** 2 * REAL (I * (I-1))  ! Bad...integer overflow

         SSQ  = SSQ  + (TERM * REAL (I)) * (TERM * REAL (I-1))  ! May be better
                                                     ! than TERM**2*R(I)*R(I-1)
   20 CONTINUE

      MEAN = XBAR
      SD   = SQRT (SSQ / REAL (N))

      RETURN
      END
