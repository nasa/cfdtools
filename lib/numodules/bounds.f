C+---------------------------------------------------------------------
C
      SUBROUTINE BOUNDS ( M, N, MDIM, F, FMIN, FMAX )
C
C  PURPOSE:
C     BOUNDS determines the maximum and minimum values contained in the
C     input array F, which is treated as 2-dimensional.  This version
C     provides for finding extremes across multiple arrays (multiple
C     calls, one per array).
C
C  ARGUMENTS:
C     ARG   TYPE I/O/S DIM    DESCRIPTION
C     M,N    I     I    -     This routine checks elements F(1:M,1:N).
C                             N=1 for a 1-D array.
C     MDIM   I     I    -     Effective row dimension of F in the
C                             routine that sets up this array.
C     F      R     I  MDIM,N  Array to be scanned for min and max.
C     FMIN,  R    I/O   -     Minimum and maximum values found (so far).
C     FMAX
C            WARNING:  Do not initialize FMIN or FMAX with extreme
C                      values, since they may never be reset because
C                      the code takes advantage of the fact that a new
C                      maximum cannot also be a new minimum.
C
C            Sample usage (1-D case):
C
C            FMIN = F1(1)
C            FMAX = FMIN
C            CALL BOUNDS ( NPTS, 1, MXPTS, F1, FMIN, FMAX )
C            CALL BOUNDS ( NPTS, 1, MXPTS, F2, FMIN, FMAX )
C
C            Then FMIN, FMAX are the extreme values across 1-D arrays F1 & F2.
C
C  ENVIRONMENT:  VAX, CRAY -- FORTRAN 77
C
C  HISTORY:
C   01/04/82    PJT    Original design and coding.
C   11/23/82    RGL    Extended to handle 2D arrays as well as vectors.
C   02/11/87    DAS    Provided for usage across multiple arrays.
C   10/27/88    DAS    Clarified usage for the 1-D array case.  Revised
C                      code to take advantage of fact that a new maximum
C                      cannot also be a new minimum.  (But use of MIN/MAX
C                      intrinsics may be better on vector machines.)
C   06/12/89    DAS    Aaargh!  Minimum is never set if inputs are
C                      monotonically increasing.  Therefore, it must be
C                      input as a legitimate value and not some extreme.
C                      Chose to warn user in header rather than change
C                      the code, but is the efficiency worth the risk?
C
C  AUTHOR:   Jeff Trosin, Informatics Inc.
C
C----------------------------------------------------------------------

      INTEGER   M, N, MDIM
      REAL      F(MDIM,N), FMIN, FMAX
      INTEGER   I, J

      DO 60, J = 1, N
         DO 50, I = 1, M
            IF ( F(I,J) .GT. FMAX ) THEN
               FMAX = F(I,J)
            ELSE IF ( F(I,J) .LT. FMIN ) THEN
               FMIN = F(I,J)
            END IF
   50    CONTINUE
   60 CONTINUE

      RETURN
      END
