C+----------------------------------------------------------------------
C
      SUBROUTINE PLEVAL ( ND, NPTS, ARC, X, Y, Z,
     >                    NEVAL, ARCEVAL, XEVAL, YEVAL, ZEVAL, IER )
C
C Acronym:  Parametric Linear curve EVALuation using arc length 
C
C Purpose:
C   Given a curve that has been linearly parameterized in arc length,
C   PLEVAL evaluates the Cartesian coordinates for the input arc length
C   parameter values.  The evaluation may be for 1, 2, or 3 dimensions.
C   Routine PLFIT can be used to parameterize the curve.
C
C Notes:
C   For 1D and 2D cases, arguments Y and/or Z do not apply.  The calling
C   program should pass dummies, or pass X more than once.
C 
C Method:
C   For each arc length value in ARCEVAL() to be evaluated:
C     - find the elements in the curve parameter array, ARC(), that 
C       bracket the value.  Sequential searches in either direction
C       are used: start at the low end for the first point; start from
C       the previous interval thereafter.
C     - use the relative position in the bracketing interval to interpolate
C       X,Y,ZEVAL() from the X,Y,Z() values of the bracketing points
C
C Arguments:
C    NAME      DIM   TYPE I/O/S DESCRIPTION
C   ND          -     I     I   Number of dimensions to interpolate
C   NPTS        -     I     I   Number of points in input curve
C   ARC       NPTS    R     I   Arc lengths of points on input curve
C   X,        NPTS    R     I   + Coordinates
C   Y,          "     "     "   | of the
C   Z           "     "     "   + input curve.
C   NEVAL       -     I     I   Number of points to be evaluated
C   ARCEVAL   NEVAL   R     I   Arc lengths of evaluation points
C   XEVAL,    NEVAL   R     O   + Coordinates
C   YEVAL,      "     "     "   | of evaluation
C   ZEVAL       "     "     "   + points
C   IER         -     I     O   Error code:
C                               IER = 0: successful completion
C                               IER =-1: bad number of dimensions
C                               IER > 0: ARCEVAL(IER) is out-of-range
C                                        of ARC(1)...ARC(NPTS).
C                                        X,Y,ZEVAL(I) undefined for
C                                        I .GE. IER.
C
C System Dependencies:
C   IMPLICIT NONE is a VAX extension.
C   Variable names longer than 6 characters are non-standard.
C
C Environment:  VAX/VMS, FORTRAN 77
C
C Author: Ron Langhi, Sterling Software, Palo Alto, CA
C 
C Modification history:
C   Dec 86  RGL  Original design and coding
C   Jan 87  dbs  Minor mods to allow for more general usage
C 
C-----------------------------------------------------------------------

      IMPLICIT NONE

C ... Arguments.

      INTEGER ND, NPTS, NEVAL, IER
      REAL    ARC(NPTS), X(NPTS), Y(NPTS), Z(NPTS),
     >        ARCEVAL(NEVAL), XEVAL(NEVAL), YEVAL(NEVAL), ZEVAL(NEVAL)

C ... Local variables.

      INTEGER N ,LOWER ,UPPER
      REAL    RATIO

C ... Execution.

      IF ( ND .GT. 3 .OR. ND .LT. 1 ) THEN
         IER = -1
         GOTO 999
      ENDIF

C ... First point uses simple search; thereafter, assume input values of
C     ARCEVAL() are often ordered, so searching up or down from previous
C     interval should be quite efficient.

      IER = 0
      LOWER = 1
      UPPER = 2

C ... For each arc length to be evaluated:

      DO 100, N = 1, NEVAL

C ...    If the arc length to evaluate is not on the curve, abort:

         IF ( ARCEVAL(N) .LT. ARC(1)  .OR.
     &        ARCEVAL(N) .GT. ARC(NPTS)   ) THEN
            IER = N
            GOTO 999
         ENDIF

C ...    Find the segment of the curve this arc length is in by linear
C        search in the proper direction from the current segment.
C        First check the current segment to see if a search is needed.

         IF ( ARCEVAL(N) .GT. ARC(UPPER) ) THEN
C ...       Search the segments above the current one.
 10         UPPER = UPPER + 1
            IF ( ARC(UPPER) .LT. ARCEVAL(N) ) GOTO 10
            LOWER = UPPER - 1

         ELSE IF ( ARCEVAL(N) .LT. ARC(LOWER) ) THEN
C ...       Search the segments below the current on.
 20         LOWER = LOWER - 1
            IF (ARC(LOWER) .GT. ARCEVAL(N) ) GOTO 20
            UPPER = LOWER + 1

C        {ELSE}
C           {The EVAL point is in this segment, so do nothing.}

         ENDIF

C ...    Use the relative position as a ratio to compute the X,Y,Z values.

         RATIO = ( ARCEVAL(N) - ARC(LOWER) )
     &         / ( ARC(UPPER) - ARC(LOWER) )

         XEVAL(N) = X(LOWER) + RATIO * ( X(UPPER) - X(LOWER) )

         IF( ND .GT. 1 ) YEVAL(N) = Y(LOWER)
     &                            + RATIO * ( Y(UPPER) - Y(LOWER) )

         IF( ND .GT. 2 ) ZEVAL(N) = Z(LOWER)
     &                            + RATIO * ( Z(UPPER) - Z(LOWER) )

 100  CONTINUE

 999  RETURN
      END
