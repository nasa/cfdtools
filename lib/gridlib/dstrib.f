C+---------------------------------------------------------------------
C
      SUBROUTINE DSTRIB (MODE, NP, P, NPTS, XMIN, XMAX, X)
C
C  PURPOSE:  DSTRIB generates a distribution of points in the range
C            [XMIN,XMAX].  A choice of distributions is offered via
C            argument MODE.  Most of the choices require additional
C            parameters,  unlike those offered by subroutine XGRID,
C            from which DSTRIB was derived.
C
C            (Two reasons to leave XGRID intact: (1) It is not easy
C            to track down all existing usages to change the calls;
C            (2) XGRID avoids fractional powers (i.e., logarithms);
C            DSTRIB does fractional exponentiation even if the  ex-
C            ponent is 1.)
C            
C  METHOD:   The uniform distribution is retained for completeness.
C            MODES 1 and 2 are generalizations of XGRID's sinusoid-
C            al distributions, providing for powers of cosine other
C            than 1.   MODE 3 provides for bunching at an arbitrary
C            internal point only,  or at an internal point and each
C            end as well.
C 
C  NOTES:    The sinusoidal distributions (MODE = -1,1,2,3)  have a
C            common feature of giving greater bunching in the  same
C            regions as "WIDTH" parameter P(1)=1.0 gives if P(1) is
C            LESS than 1.0 (and greater than zero - i.e., fraction-
C            al power as opposed to squaring, cubing, etc.).   When
C            P(1) > 1. is used,  results are less easily described,
C            and likely to be other than what was intended.  Brief-
C            ly, the higher this exponent, the less the bunching at
C            the regions given by P(1)=1.,  and the more the bunch-
C            ing towards either XMAX (MODE=1) or the midpoint (MODE
C            =2) or the "CENTER" (P(2), MODE=3).
C
C  ARGUMENTS:
C    ARG      DIM  TYPE  I/O/S  DESCRIPTION
C    MODE      -     I     I    MODE controls type of distribution,
C                               as described below.
C    NP        -     I     I    Number of input parameters P(*). In
C                               the case of MODE=3, NP should be 3,
C                               not 2, to allow for returning info.
C    P        NP     R     I    Parameters depending on MODE below.
C    NPTS      -     I     I    No. of pts in desired distribution.
C  XMIN,XMAX   -     R     I    First and last values to include.
C    X       NPTS    R     O    Desired distribution of points.
C
C  USAGE:
C       MODE = -1 means same as for MODE = 1 except the bunching is
C                 is towards XMAX instead of XMIN.
C       MODE = 0  means uniform distribution. Use NP=1; P not used.
C       MODE = 1  means sinusoidal bunching towards XMIN if P(1) is
C                 in the range (0.,1.].   P(1) > 1. gives bunchings
C                 at both XMIN and XMAX that differ from each other
C                 in general, and are probably not what you want.
C       MODE = 2  means sinusoidal bunchings of the  same  type  at
C                 both XMIN and XMAX if P(1) is in (0.,1.].   These
C                 bunchings are reduced, and additional bunching is
C                 given at the mid-point, if P(1) > 1. is used. The
C                 distribution is symmetric about the center point.
C       MODE = 3  means sinusoidal bunching around the internal pt.
C                 defined by P(2), where XMIN < P(2) < XMAX.   P(1)
C                 controls bunching as for MODE=2,  with additional
C                 bunching at the ends of the interval if P(1) > 1.
C
C                 MODE 3 NOTES:
C
C              1. P(3) is RETURNED with  the corresponding interior
C                 index M (as a REAL) such that X(M) = P(2).
C
C              2. If X = P(2) is to be included exactly (as assumed
C                 here),  ensuring smoothness in  dX is nontrivial.
C                 DSTRIB does NOT attempt to ensure this.   It does
C                 no more than determine the  proportion  of points
C                 to be placed on either side of the interior point
C                 then generate two similar distributions.  Routine
C                 SMOOTHX can be used to smooth the result if reqd.
C
C  ENVIRONMENT:  VAX/VMS FORTRAN 77
C
C  HISTORY:
C  05/24/85    DAS    Added NP, P(*) arguments to handle generation
C                     of internally-bunched distributions (MODE=3);
C                     generalized XGRID's options (MODE=1,2) in the
C                     process.  Left XGRID intact.
C  06/25/85    DAS    Revised MODE=3: straightforward sine function
C                     used now, where the modified sine function of
C                     subroutine BEVAL had been used originally. It
C                     turned out that the latter gives a  different
C                     type of distribution at the end  points  that
C                     is hard to predict and is therefore unwise.
C  01/22/87    DAS    Patch in MODE=3 to handle symmetric case like
C                     P(2)=0.5 on [0.,1.] with odd N properly. (For
C                     even N, there is no symmetric solution here.)
C                     Also: had to return subscript of interior pt.
C                     Also: introduced MODE=-1 for completeness.
C  08/14/89    DAS    Return MODE=3's interior point index in P(3),
C                     not NP.
C
C  AUTHOR:  David Saunders, Sterling Software/NASA Ames, Palo Alto, CA.
C
C----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments.

      INTEGER
     >   MODE, NP, NPTS
      REAL
     >   P (NP), X (NPTS), XMAX, XMIN

C     Local constants.

      REAL
     >   EPS, HALF, ONE
      PARAMETER
     >  (EPS = 1.E-6, HALF = 5.E-1, ONE = 1.E+0)

C     Local variables.

      INTEGER
     >   I, MIDDLE, NPTS1, NPTS2
      REAL
     >   DTHETA, DX, PIBY2, RANGE, YNORM

C     Execution.


      PIBY2  = ASIN (ONE)
      DTHETA = PIBY2 / REAL (NPTS - 1)
      RANGE  = XMAX - XMIN


      IF (MODE .EQ. 0) THEN

C  *     Uniform distribution:

         DX = RANGE / REAL (NPTS - 1)

         DO 50 I = 2, NPTS - 1
            X (I) = REAL (I - 1) * DX + XMIN
   50    CONTINUE


      ELSE IF (MODE .EQ. 1) THEN

C  *     Sinusoidal bunching at the low end if 0.0 < P (1) <= 1.0;
C        unsymmetric bunching at both ends if P (1) > 1.0:

         DO 100 I = 2, NPTS - 1
            X (I) = XMIN +
     +         RANGE * (ONE - (COS (REAL (I-1) * DTHETA)) ** P (1))
  100    CONTINUE


      ELSE IF (MODE .EQ. -1) THEN

C  *     Sinusoidal bunching at the high end if 0.0 < P (1) <= 1.0;
C        unsymmetric bunching at both ends if P (1) > 1.0:

         DO 150 I = 2, NPTS - 1
            X (I) = XMIN + RANGE * (SIN (REAL (I-1) * DTHETA)) ** P (1)
  150    CONTINUE


      ELSE IF (MODE .EQ. 2) THEN

C  *     Symmetric sinusoidal bunching at both ends if 0. < P (1) <= 1.;
C        additional bunching in the middle if P (1) > 1.0:

         DTHETA = DTHETA + DTHETA
         RANGE  = RANGE * HALF
         MIDDLE = NPTS / 2

         DO 200 I = 2, MIDDLE
            YNORM = ONE - (COS (REAL (I-1) * DTHETA)) ** P (1)
            X (I)        = XMIN + RANGE * YNORM
            X (NPTS+1-I) = XMAX - RANGE * YNORM
  200    CONTINUE

C  *     The following avoids difficulties with cosine at what is
C        supposed to be exactly PI/2 but may be slightly greater,
C        giving a negative cosine, when NPTS is odd:

         IF (MIDDLE + MIDDLE .LT. NPTS) X (MIDDLE+1) = XMIN + RANGE


      ELSE IF (MODE .EQ. 3) THEN

C  *     Modified sinusoidal bunching around internal point P(2) if
C        0.0 < P(1) <= 1.0; bunching also towards the end-points if
C        P(1) > 1.0.    The following gives proportional numbers of
C        points in the two intervals involved.    For instance,  if
C        NPTS = 100 and CENTER = XMIN + 0.7 * RANGE, then 70 of the
C        points should be in [XMIN,CENTER] and 30 in [CENTER,XMAX],
C        with bunching towards CENTER (and possibly towards the end
C        points XMIN and XMAX if P(1) > 1.).
C
C        See NOTE above about lack of continuity in dX across P(2).
C
C        Patch for symmetric case with odd N:  we need to round up,
C        but e.g. for N = 101 on [0.,1.] with roundoff we might get
C        NINT (50.49999...) = 50, not 51 - hence the EPS.
C
C        Also:  we need to return the subscript of the interior pt.
C        NP was the only available argument - caller beware.

         NPTS1 = NINT (REAL (NPTS) * (P (2) - XMIN) / RANGE + EPS)
         DTHETA = PIBY2 / (NPTS1 - 1)
         RANGE = P (2) - XMIN

         DO 300 I = 2, NPTS1 - 1
            X (I) = XMIN + RANGE * (SIN (REAL (I-1) * DTHETA)) ** P (1)
  300    CONTINUE

         P (3) = REAL (NPTS1)
         X (NPTS1) = P (2)
         NPTS2 = NPTS - NPTS1 + 1
         DTHETA = PIBY2 / REAL (NPTS2 - 1)
         RANGE = XMAX - P (2)

         DO 400 I = 2, NPTS2 - 1
            X (I+NPTS1-1) = P (2) + RANGE *
     +         (ONE - (SIN (PIBY2 + REAL (I-1) * DTHETA)) ** P (1))
  400    CONTINUE
         
      END IF


C  *  Ensure precise values at end-points:

      X (1) = XMIN
      X (NPTS) = XMAX

      RETURN
      END
