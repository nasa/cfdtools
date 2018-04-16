C+----------------------------------------------------------------------
C
      SUBROUTINE LINT2D (MAXF, MAXX, MAXY, NF, NX, NY, F, X, Y,
     >                   NTARGS, XTARGS, YTARG, FINTERP)
C
C ACRONYM: Linear INTerpolation in 2 Dimensions (multi-function table)
C          -      ---              - -
C PURPOSE:
C        LINT2D returns linearly-interpolated (or extrapolated) values
C     from a given 2-D table, which should be rectangular (NX entries
C     in each row and NY entries in each column, with all coordinates
C     defined by two vectors).  More than one function may be associated
C     with each point of the table, and the table coordinates may be
C     monotonic increasing or decreasing (and not necessarily uniform).
C
C        This version is somewhat optimized for interpolating a new ROW
C     of values - new columns would need to be done one call per row.
C     (This compromise is considered better than handling just one point
C     per call, while passing in a rectangular set of target points is
C     considered too cumbersome.)
C
C METHOD:
C        The efficient search offered by subroutine INTERVAL is used in
C     each direction.  By allowing for interpolating a row at a time,
C     much of the repeated searching implied by interpolating from one
C     rectangular mesh to another is eliminated.  Providing for multiple
C     functions likewise can eliminate much repeated searching.
C
C ARGUMENTS:
C    ARG    DIM   TYPE I/O/S DESCRIPTION
C   MAXF     -      I    I   Max. # functions provided for in F(*,i,j)
C   MAXX     -      I    I   Max. # X-coords.   "   "   "   " F(n,*,j)
C   MAXY     -      I    I   Max. # Y-coords.   "   "   "   " F(n,i,*)
C                            (MAXY is redundant but retained for symmetry.)
C    NF      -      I    I   Actual number of functions at each pt. (NF>=1)
C    NX      -      I    I   Actual number of columns in table (NX>=2)
C    NY      -      I    I   Actual number of rows    in table (NY>=2)
C    F  MAXF,MAXX,MAXY   I   F(1:NF,1:NX,1:NY) are the meaningful function
C                            values in the table
C    X      NX      R    I   Table coordinates in the 1st dimension
C    Y      NY      R    I   Table coordinates in the 2nd dimension
C   NTARGS   -      I    I   Number of points in a "row" of target points;
C                            NTARGS >= 1
C   XTARGS  NTARGS  R    I   X-coordinates of target points
C   YTARG    -      R    I   Y-coordinate of (all) target points
C   FINTERP MAXF,NTARGS  O   One "row" of interpolated values is returned
C                            in FINTERP(1:NF,1:NTARGS)
C                     
C USAGE NOTES:
C
C   There are at least two ways to extrapolate with LINT2D:
C
C      (1) Linear extrapolation is obtained just as for linear interpolation
C          if a target point is outside the table and LINT2D is used normally.
C
C      (2) A restricted form of extrapolation can be obtained simply if the
C          target point is adjusted as follows (suggested by Charles Hooper):
C
C                XT = MAX (X (1), MIN (X (NX), XTARG))
C                YT = MAX (Y (1), MIN (Y (NY), YTARG))
C
C                CALL LINT2D (MAXF, MAXX, MAXY, NF, NX, NY, F, X, Y,
C               >             1, XT, YT, FINTERP)
C
C          where the single-target point case with ascending table coordinates
C          is shown for simplicity.  This has the effect of extending table
C          function values at the edges of the table into regions outside the
C          table - a desirable choice for some applications.
C          
C
C EXTERNAL REFERENCES:
C  INTERVAL   1-D search utility
C
C ERROR HANDLING:
C     None, for reasons of efficiency.  The calling program should ensure
C     that NX >=2, etc., and that the table coordinates are monotonic.
C     Since extrapolation is permitted, there is nothing else to consider.
C
C ENVIRONMENT:  VAX/VMS FORTRAN 77
C
C DEVELOPMENT HISTORY:
C     DATE  INITIALS  DESCRIPTION
C   10/29/87   DAS    Initial design and code, adapted from TABLE2.
C   11/10/86   DAS    Replaced flexible IGNORE(1:MAXF) argument with
C                     simpler NF argument - probably sufficient.
C
C AUTHOR: David Saunders, Sterling Software, Palo Alto, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER
     >   MAXF, MAXX, MAXY, NF, NTARGS, NX, NY
      REAL
     >   F (MAXF, MAXX, MAXY), FINTERP (MAXF, NTARGS), X (NX),
     >   XTARGS (NTARGS), Y (NY), YTARG

C     Local constants:

      REAL
     >   ONE
      PARAMETER
     >  (ONE = 1.E+0)

C     Local variables:

      INTEGER
     >   I, ITARG, J, N
      REAL
     >   P, PM1, Q

C     Procedures:

      EXTERNAL
     >   INTERVAL

C     Execution:

C     Just one target Y to deal with:

      J = 1
      CALL INTERVAL (NY, Y, YTARG, SIGN (ONE, Y (2) - Y (1)), J)

      Q = (YTARG - Y (J)) / (Y (J+1) - Y (J))

C     One or more target (X, Y)s with this Y coordinate:

      PM1 = SIGN (ONE, X (2) - X (1))
      I = 1

      DO 30, ITARG = 1, NTARGS

         CALL INTERVAL (NX, X, XTARGS (ITARG), PM1, I)

         P = (XTARGS (ITARG) - X (I)) / (X (I+1) - X (I))

C        One or more functions at this target (X, Y):

         DO 20, N = 1, NF
            FINTERP (N, ITARG) = (ONE - P) * (ONE - Q) * F (N, I, J)   +
     >                                   P * (ONE - Q) * F (N, I+1, J) +
     >                           (ONE - P) * Q         * F (N, I, J+1) +
     >                                   P * Q         * F (N, I+1, J+1)
   20    CONTINUE
   30 CONTINUE
      
      RETURN         ! End of LINT2D
      END
