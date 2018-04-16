C+----------------------------------------------------------------------
C
      SUBROUTINE LINT1D ( MXFUNC, NFUNC, NTABLE, XTABLE, FTABLE, NX, 
     >                    X, FINTERP )
C
C ACRONYM: Linear INTerpolation/extrapolation (1D) (multi-function table)
C          -      ---                          --
C
C PURPOSE: For each X(I), LINT1D returns linearly-interpolated function
C          values FINTERP(1:NFUNC,I) from the given 1D table whose abscissas 
C          must be monotone increasing or monotone decreasing.
C
C          This version permits extrapolation.
C
C          LINT1D is a multi-function adaptation of LINE1D, which is the
C          subroutine form of FUNCTION TABLE1.  Each has its place, although
C          LINE1D could well be dispensed with.
C
C METHOD:  Nothing more than linear interpolation (2-point formula)
C          is attempted, partly because multi-point formulas demand
C          special treatment at the ends of the table, and partly to
C          permit "safe" extrapolation, which is often due to round-
C          off or is otherwise intended.  Furthermore, piecewise line
C          segments (as the table may well represent) MUST be treat-
C          ed linearly.
C
C ARGUMENTS:
C    ARG    DIM   TYPE I/O/S DESCRIPTION
C  MXFUNC    -      I    I   Max. no. of function values at each point 
C                            in the table; MXFUNC >= NFUNC.
C  NFUNC     -      I    I   Requested no. of function values at each 
C                            point in the table; NFUNC >= 1 assumed.
C                            This allows for interpolating fewer than
C                            MXFUNC function values.
C  NTABLE    -      I    I   Length of the table; NTABLE >= 1 assumed.
C  XTABLE  NTABLE   R    I   Abscissas in table, assumed monotonic -
C                            increasing or decreasing.
C  FTABLE MXFUNC,NTABLE   
C                   R    I   Ordinates in table.
C    NX      -      I    I   Number of points in table at which 
C                            interpolations are requested.
C    X       NX     R    I   Abscissas at which interpolation is 
C                            requested.
C  FINTERP MXFUNC,NX  
C                   R    O   Interpolated/extrapolated function values
C                            returned in FINTERP(1:NFUNC,1:NX).
C
C EXTERNAL REFERENCES:
C    NAME     DESCRIPTION
C  INTERVAL   Interpolation search for interval containing a point.
C
C ERROR HANDLING:
C     None - checking for monotonicity in XTABLE(*) would compromise
C     efficiency too much.  FORTRAN handles NTABLE <= 0.
C
C ENVIRONMENT:  VAX/VMS FORTRAN 77
C
C DEVELOPMENT HISTORY:
C     DATE  INITIALS  DESCRIPTION
C   06-06-88   CLH    Adapted from LINE1D (David Saunders) to interpolate
C                     one or more functions at each point in the line.
C
C AUTHOR: Charles L. Hooper, Sterling Software, Palo Alto, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT    NONE

C --- Arguments:

      INTEGER     MXFUNC, NFUNC, NTABLE, NX
      REAL        XTABLE(NTABLE), X(NX), FTABLE(MXFUNC,NTABLE), 
     >            FINTERP(MXFUNC,NX)

C --- Local variables:

      INTEGER     I, LEFT, N
      REAL        PM1, RATIO

C --- Execution:

      IF ( NTABLE.EQ.1 ) THEN

C ---    Only one table entry:

         DO 100 N=1,NFUNC
            FINTERP(N,1) = FTABLE(N,1)
  100    CONTINUE

      ELSE

         PM1  = SIGN ( 1.E+0, XTABLE(2) - XTABLE(1) )
         LEFT = 1

         DO 300 I=1,NX

C ---       Search for the best "left-hand" endpoint to use for 
C           interpolation:

            CALL INTERVAL ( NTABLE, XTABLE, X(I), PM1, LEFT )

C ---       Linear interpolation or extrapolation.  Order does not matter:

            RATIO = ( X(I) - XTABLE(LEFT) ) /
     >         ( XTABLE(LEFT+1) - XTABLE(LEFT) )

            DO 200 N=1,NFUNC
               FINTERP(N,I) = FTABLE(N,LEFT) + 
     >            ( FTABLE(N,LEFT+1) - FTABLE(N,LEFT) ) * RATIO
  200       CONTINUE
  300    CONTINUE

      END IF

  999 RETURN
      END
