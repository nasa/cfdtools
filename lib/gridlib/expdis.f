C+---------------------------------------------------------------------
C
      SUBROUTINE EXPDIS ( MODE, BETA, N, RATIOS )
C
C  PURPOSE:  EXPDIS (EXPonential DIStribution) generates N values in
C            the range [0,1] which vary more and more nonlinearly as
C            parameter BETA tends towards 1 from above.  Such a dis-
C            tribution may be generated once then used repeatedly to
C            fill in 2- or 3-D meshes as follows:
C
C            X(I)  =  Xmin  +  (Xmax - Xmin) * RATIOS(I)     I = 1:N
C            
C            The bunching may be one-sided (MODE = 1 or MODE = 3) or
C            two-sided, in which case the distribution is symmetric.
C
C  METHOD:   EXPDIS modularizes formulas found in SCRAM2D/SCRAM3D by
C            Ajay Kumar. The algebra is obscure - these are not ord-
C            inary exponential distributions. It is not obvious that
C            the two-sided bunching is symmetric.  However, RATIO(I)
C            = 1 - RATIO(N+1-I) in observed results.  The code takes
C            advantage of this.
C
C  ARGUMENTS:
C    ARG    TYPE  I/O/S  DIM   DESCRIPTION
C    MODE     I     I     -    MODE=1 gives bunching towards 1.
C                              MODE=2 gives bunching at 0. and 1.
C                              MODE=3 gives bunching towards 0.
C    BETA     R     I     -    BETA > 1 controls bunching at the end
C                              points.  Typical values:  1.02 - 1.2.
C                              Closer to 1. gives more bunching.
C    N        I     I     -    Required number of points. N > 1.
C    RATIOS   R     O     N    Values varying nonuniformly from 0 to
C                              1, inclusive.
C
C  ERROR HANDLING:
C    Avoid a logical unit number for diagnostic output by using
C    STOP 'EXPDIS: .....' if BETA <= 1.0, or N <= 1.
C
C  ENVIRONMENT:  VAX/VMS FORTRAN 77
C
C  HISTORY:
C    08/21/85    DAS    Original design, from Kumar's code (2-sided case).
C    08/23/85    DAS    Introduced MODE to provide one-sided case too.
C    05/19/86    dbs    Added MODE=3 for bunching towards 0.
C    11/11/87    dbs    Fix for MODE=3 to calculate values in reverse order.
C    05/27/89    DAS    X(I) = XA + (XB - XA) * R(I) form is better (above).
C    06/15/89     "     Small revision involving RNM1 was wrong for MODE=2.
C
C  AUTHOR: David Saunders, Sterling Software, CA; Ajay Kumar, NASA/Langley
C
C----------------------------------------------------------------------

C ... Arguments:

      INTEGER    MODE, N
      REAL       BETA, RATIOS(N)

C ... Local constants:

      REAL       HALF, ONE
      PARAMETER (HALF = 0.5E+0, ONE = 1.E+0)

C ... Local variables:

      INTEGER    I
      REAL       B, FI, RNM1

C ... System functions:

      INTRINSIC  EXP, FLOAT, LOG

C ... Execution:

      IF ( BETA .LE. ONE  .OR.  N .LE. 1 ) STOP 'EXPDIS: Bad BETA or N'

C ... Note that for F(I) = A ** P(I), where A is independent of I,
C     it is cheaper to use the form  F(I) = EXP ( LOG(A) * P(I) ):

      B = LOG ( ( BETA + ONE ) / ( BETA - ONE ) )
      RNM1 = B / FLOAT ( N - 1 )

      IF ( MODE .EQ. 2 ) THEN

C ...    Symmetric bunching towards both 0 and 1:

         DO 200 I = 1, (N+1)/2
            FI = EXP ( RNM1 * FLOAT ( 2*I-N-1 ) )
            RATIOS(I) = ( ONE + BETA * ( FI-ONE ) / ( FI+ONE ) ) * HALF
            RATIOS(N+1-I) = ONE - RATIOS(I)
  200    CONTINUE

      ELSE IF ( MODE .EQ. 1 )THEN

C ...    Bunch towards 1:

         DO 100 I = 1, N
            FI = EXP ( RNM1 * FLOAT (I-1) )
            RATIOS(I) = BETA * ( FI-ONE ) / ( FI+ONE )
  100    CONTINUE


      ELSE IF ( MODE .EQ. 3 )THEN

C ...    Bunch towards 0 by subtracting the MODE = 1 ratios from 1 and
C        calculating them in reverse order:

         DO 300 I = 1 ,N
            FI = EXP ( RNM1 * FLOAT (N-I) )
            RATIOS(I) = ONE - BETA * ( FI-ONE ) / ( FI+ONE )
  300    CONTINUE

      END IF

      RETURN
      END
