C+----------------------------------------------------------------------
C
      SUBROUTINE ZEROS (NX, X, F, MXNZS, NZS, XZS, DCREAS, IER)
C
C PURPOSE: ZEROS  locates, in crude fashion,  all of the zeros (if any)
C          of the given one-variable function represented  by  discrete
C          data points in some range.
C
C METHOD:  Search the discretized data from left to right, noting where
C          the function changes sign.   Once a zero has been bracketed,
C          use just linear interpolation to estimate it.  Indicate whe-
C          ther the function was increasing or decreasing, for the case
C          where it is really the maxima/minima of some other  function
C          that is being sought.
C
C NOTES:   The finer the discretization, the more precise (and more ex-
C          pensive) this approach is.  But traditional zero finders ex-
C          pect the function to have opposite signs at each end of  the
C          given interval, with no attempt to identify more than 1 zero
C          per call.  It is not clear how to use such routines here.
C
C ARGUMENTS:
C    ARG    DIM  TYPE I/O/S DESCRIPTION
C    NX      -     I    I   Number of discrete data points
C    X      NX     R    I   Abscissas of data points (monotonic, incr.)
C    F      NX     R    I   Function values corresponding to X(*)
C    MXNZS   -     I    I   Max. room provided in NZS(*), DCREAS(*)
C    NZS     -     I    O   Number of zeros identified in [X(1),X(NX)]
C    XZS    NZS    R    O   Abscissas of zeros found
C    DCREAS NZS    L    O   DCREAS(J)=T means F(X) is decreasing around
C                           the Jth zero found.   If F(X) is really the
C                           1st derivative of, say,  G(X),  then XZS(J)
C                           estimates the location of a maximum of G(X)
C                           else it represents a minimum.    (Points of
C                           inflection in G(X) will not be noticed.)
C    IER     -     I    O   IER=0 means no problems encountered;
C                           IER=2 means MXNZS was not big enough.
C
C ENVIRONMENT: FORTRAN 77
C
C HISTORY:  06/27/83   DAS  Initial design and code.
C
C AUTHOR: David Saunders, Informatics General, Palo Alto, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER  NX, MXNZS, NZS, IER
      REAL     F (NX), X (NX), XZS (MXNZS)
      LOGICAL  DCREAS (MXNZS)
        
C     Local variables:

      INTEGER  I
      REAL     SIGNL, SIGNR

C     Execution:

      NZS = 0
      SIGNL = SIGN (1.E+0, F (1))

      DO 200, I = 2, NX
         SIGNR = SIGN (1.E+0, F (I))
         IF (SIGNR .NE. SIGNL) THEN
            NZS = NZS + 1
            IF (NZS .GT. MXNZS) GO TO 900

            XZS (NZS) = X (I-1) - F (I-1) * (X (I) - X (I-1)) /
     >                                      (F (I) - F (I-1))
            DCREAS (NZS) = SIGNL .GT. 0.E+0
            SIGNL = SIGNR
         END IF
  200 CONTINUE

      IER = 0
      GO TO 999

  900 IER = 2

  999 RETURN
      END
