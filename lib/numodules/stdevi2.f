C+----------------------------------------------------------------------
C
      SUBROUTINE STDEVI2 (N, COUNTS, STRIDE, MEAN, STDEV)
C
C ACRONYM: STandard DEViation for I*2 data
C          --       ---           - -
C PURPOSE:
C     STDEVI2 returns the mean and standard deviation of N elements of
C     an INTEGER*2 array COUNTS(*).  The elements must be some constant
C     stride apart.  Element 1 is taken as the first to be included.
C
C METHOD:
C     One pass through the data is the main reason for STDEVI2.  This
C     requires squaring each count as opposed to its deviation from the
C     mean, but use of double precision accumulators (~16 significant
C     decimal digits) means no truncation with I*2 data unless (in the
C     worst case),
C
C             N * 32767 ** 2  ~  10 ** 16    (round-off in 16th place)
C             N  ~  10 ** 7
C
C     meaning up to 10 megawords of data are handled safely.  This should
C     accommodate typical applications to analog-to-digital data.  (The
C     2-byte integer data from typical A/D converters prompted this
C     utility.)  The routine aborts with a diagnostic if N > 10 million.
C
C     For N = 1 million and STRIDE = 1, STDEVI2 uses ~15 CPU secs. on a
C     VAX-11/780 with 8 megabytes and a 400-page working set (and about
C     0.1 sec. for the more likely N = 10,000 case).
C
C     For extremely large N (> 10 million), a two-pass method working with
C     data offset by the mean may be preferable will be slower.
C
C     Usage:
C
C        Elements 1, 4, 7, ...   CALL STDEVI2 (N, IBUF, 3, R4MEAN, R4SD)
C
C        Elements I, I+M, I+2M, ...   "   "   (N, IBUF (I), M, XBAR, SD)
C
C ARGUMENTS:
C   ARG     DIM    TYPE  I/O/S DESCRIPTION 
C   N        -      I*4    I   Number of elements in each sum; N > 0
C   COUNTS   *      I*2    I   Elements 1, 1+STRIDE, 1+2*STRIDE, ... of
C   STRIDE   -      I*4    I     COUNTS (*) are processed
C   MEAN     -      R*4    O   Mean of specified elements
C   STDEV    -      R*4    O   Standard deviation s, where s**2 is given
C                                by  ( Sum (counts**2) ) / N  -  mean**2
C                                rather than the more proper use of N-1.
C
C ENVIRONMENT:  VAX/VMS, FORTRAN 77
C               INTEGER*2, INTEGER*4, REAL*4, REAL*8 are used for clarity.
C
C HISTORY:
C   01/14/88  DAS  Initial implementation for INTEGER*2 data and one pass.
C   02/10/88  DAS  Safeguarded square root because of round-off even if
C                  N is not huge.
C   04/18/88  DAS  Slip-up in the final steps - should be double precision
C                  throughout; raised cut-off to 10 million from 1 million
C                  after comparisons with REAL*16 version by Charley Hooper.
C
C AUTHOR:  David Saunders, Sterling Software, Palo Alto, CA.
C
C-----------------------------------------------------------------------

C     Arguments:

      INTEGER*2  COUNTS (*)
      INTEGER*4  N, STRIDE
      REAL*4     MEAN, STDEV
      
C     Local variables:

      INTEGER*4  I 
      REAL*8     SUM1, SUM2  ! Avoid other locals: registers should do it

C     Constants:

      REAL*8     ZERO
      PARAMETER (ZERO = 0.D+0)

C     Execution:

      IF (N .GT. 10000000)
     >   STOP 'STDEVI2: N > 10 million. Beware of St.Dev.!'

      SUM1 = ZERO
      SUM2 = ZERO
      DO 20, I = 1, 1 + (N - 1) * STRIDE, STRIDE
         SUM1 = SUM1 + DBLE (COUNTS (I))
         SUM2 = SUM2 + DBLE (COUNTS (I)) ** 2
   20 CONTINUE

C     Stay in double precision as long as possible:

      SUM1 = SUM1 / DBLE (N)
      MEAN = REAL (SUM1)
      STDEV = REAL (SQRT (MAX (ZERO, SUM2 / DBLE (N) - SUM1 ** 2)))

      RETURN
      END
