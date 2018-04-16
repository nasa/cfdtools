C+----------------------------------------------------------------------
C
      SUBROUTINE ASRCH3D (N, X, Y, Z, T, LEFT, TLEFT, RIGHT, TRIGHT,
     >   NEW)
C
C     One-liner: Interpolation search along a parametric curve in 3-space
C     ----------
C
C     Description and usage:
C     ----------------------
C
C        ASRCH3D is the X-Y-Z analog of ARCSRCH, from which the following
C     description is retained:
C
C        ARCSRCH modularizes the arclength interval search originally
C     written for PLSFIT (Parametric Local Spline FIT). An arclength
C     interval is sought which contains some specified value, T, or
C     which is at least the nearest interval if T is off-scale.  The
C     condition for a bracket is: TLEFT <= T <= T (LEFT+1). A logical
C     flag, NEW, must be set .TRUE. on the first call for a new set of
C     data, and the total arclength must be supplied as the initial value
C     of input variable TRIGHT. Subsequent calls with the same data make
C     use of local SAVEd variables to avoid unnecessary recalculation.
C
C        There is minimal error checking in this low-level routine. The
C     calling program is assumed to have verified that N >= 2. Efficiency
C     will benefit from passing the best estimate available, usually just
C     the result of the last call.
C
C        This is not a fully independent utility - a number of quantities
C     are shared between PLSFIT and ARCSRCH. Speed is thereby enhanced at
C     the expense of generality and ease-of-use. With some care, ARCSRCH
C     could be transplanted into another setting (see PLSFIT for usage).
C
C        The interpolation search was adapted from ideas in Sedgewick's
C     book, referenced below.
C
C     Arguments:
C     ----------
C
C     Name     Dimension  Type  I/O/S  Description
C     N                    I    I      Number of points in arrays X, Y, Z;
C                                      must be >= 2; no check performed.
C
C     X, Y, Z     N        R    I      Array of distinct points defining
C                                      the set of intervals to be examined;
C                                      not checked.
C
C     T                    R    I      The distance along the chord for
C                                      which a bracketing interval is sought.
C                                      Normally in the range zero to total
C                                      chord length, but may lie outside.
C
C     NOTE: If NEW = .TRUE., only the value of TRIGHT needs to be supplied.
C
C     LEFT                 I    I/O    Input: estimate of index of left
C                                      endpoint of the interval containing
C                                      the specified point. Must be in the
C                                      range [1, N-1]; not error checked.
C                                      LEFT < RIGHT is assumed.
C
C                                      Output: index of the largest array
C                                      value <=, or sometimes <, specified
C                                      point - see Notes. Special case for
C                                      data out of range: returns left
C                                      endpoint of closest interval.
C
C     TLEFT                R    I/O    Arclength up to index LEFT.
C
C     RIGHT                I    I/O    Input: estimate of index of right
C                                      endpoint of the interval containing
C                                      the specified point. Must be in the
C                                      range [2, N]; not error checked.
C                                      RIGHT > LEFT is assumed.
C
C                                      Output: index of the largest array
C                                      value >, or sometimes >=, specified
C                                      point - see Notes. Special case for
C                                      data out of range: return right
C                                      endpoint of closest interval.
C
C     TRIGHT               R    I/O    Arclength up to index RIGHT. NOTE:
C                                      when NEW = .TRUE., TRIGHT must be
C                                      the total arclength of the curve.
C                                      This trick permits PLSFIT3D/ASRCH3D
C                                      to avoid unnecessary recalculation.
C
C     NEW                  L    I/O    Must be set .TRUE. when ASRCH3D is
C                                      first called for any dataset.
C
C                                      Set to .FALSE. on output.
C
C     Environment:  Digital VAX-11/780, VMS FORTRAN
C     ------------
C
C     Notes:
C     ------
C
C     (1)  IMPLICIT NONE and 8-character symbols are not (yet) standard.
C
C     (2)  The arclength calculations are approximations based on straight
C          line segments between the data points.
C
C     (3)  The algorithm is designed to return bracket endpoints such that
C          TLEFT <= T < TRIGHT = T (LEFT+1), with strict inequality on the
C          right, but roundoff errors in the chord calculations sometimes
C          interfere. We protect the routine from failures such as collapse
C          of the interval to zero, i.e., LEFT = RIGHT, by checking ahead
C          when the endpoints are adjusted during the main loop, and simply
C          accept the occasional odd result. Note also that the arclengths
C          of the interval's endpoints should be regarded as being slightly
C          fuzzy, perhaps a few parts in 10**7 for typical single precision.
C          Neither condition is likely to be a problem in typical graphics
C          applications, and the calling routine can easily check the result
C          if it is critical.
C
C     Bibliography:
C     -------------
C
C     (1) Sedgewick, R.  Algorithms.  Reading: Addison-Wesley, 1983.
C            (Chap. 14)
C
C     Author:  Robert Kennelly, Sterling Software/NASA Ames, Palo Alto, CA.
C     -------
C
C     History:
C     --------
C
C     25 Mar. 1988    RAK    ARCSRCH adapted from INTERVAL and PLSFIT.
C     16 Aug. 1989  DAS/RAK  ASRCH3D adapted from ARCSRCH.
C
C-----------------------------------------------------------------------

C     Declarations.
C     -------------

      IMPLICIT NONE

C     Constants.

      REAL
     >   ZERO
      PARAMETER
     >  (ZERO = 0.0E+0)

C     Arguments.

      LOGICAL
     >   NEW
      INTEGER
     >   LEFT, N, RIGHT
      REAL
     >   TLEFT, TRIGHT, X (N), Y (N), Z (N), T

C     Local variables.

      INTEGER
     >   LENGTH, NLESS1, TRIAL
      REAL
     >   LSENTRY, RSENTRY, TEMP1, TEMP2, URTRIGHT

C     Storage.

      SAVE
     >   NLESS1, LSENTRY, RSENTRY, URTRIGHT

C     Procedures.

      REAL
     >   CHORD3D
      EXTERNAL
     >   CHORD3D

C     Execution.
C     ----------

      IF (NEW) THEN

C        Initialization for new data set. Note special use of TRIGHT here,
C        to avoid having to repeat the relatively expensive summation over
C        the entire curve.

         NEW      = .FALSE.
         URTRIGHT = TRIGHT
         NLESS1   = N - 1

         LSENTRY  = CHORD3D (X, Y, Z, 1, 2)
         RSENTRY  = TRIGHT - CHORD3D (X, Y, Z, NLESS1, N)

C        The following will be appropriate only if the "simplification"
C        below doesn't apply, but we initialize everything here to get it
C        over with.

         LEFT   = 2
         TLEFT  = LSENTRY
         RIGHT  = NLESS1
         TRIGHT = RSENTRY
      END IF

C     Simplify things by disposing of two important special cases so that
C     TLEFT and TRIGHT can really bracket T. As a byproduct, this also
C     takes care of the N = 2, 3 cases (one or two intervals).

      IF (T .LT. LSENTRY) THEN       ! First interval applies
         LEFT   = 1
         TLEFT  = ZERO
         RIGHT  = 2
         TRIGHT = LSENTRY
         GO TO 990
      ELSE IF (T .GE. RSENTRY) THEN  ! Last interval applies
         LEFT   = NLESS1
         TLEFT  = RSENTRY
         RIGHT  = N
         TRIGHT = URTRIGHT
         GO TO 990
      END IF

C      ----------------------------------------------------------------
C     |                                                                |
C     |   N > 3                                                        |
C     |                                                                |
C     |   1 <= LEFT < RIGHT <= N                                       |
C     |                                                                |
C     |   LSENTRY <= T < RSENTRY                                       |
C     |                                                                |
C      ----------------------------------------------------------------

C     Refine bracket estimate, checking in particular whether the current
C     values are already correct.

      IF (T .GE. TLEFT) THEN
         IF (T .LT. TRIGHT) THEN

C           T is somewhere in the original interval - are we done?

            IF (RIGHT - LEFT .LE. 1) GO TO 990
         ELSE

C           We'll look farther to the right. Evidently LEFT was < N - 2.

            LEFT   = RIGHT
            TLEFT  = TRIGHT
            RIGHT  = NLESS1
            TRIGHT = RSENTRY
         END IF
      ELSE

C        Look to the left of the guess. Evidently LEFT was > 2.

         RIGHT  = LEFT
         TRIGHT = TLEFT
         LEFT   = 2
         TLEFT  = LSENTRY
      END IF

C      ----------------------------------------------------------------
C     |                                                                |
C     |   2 <= LEFT < RIGHT <= N - 1                                   |
C     |                                                                |
C     |   LSENTRY <= TLEFT <= T < TRIGHT <= RSENTRY                    |
C     |                                                                |
C      ----------------------------------------------------------------

C     The interval length must decrease each search iteration. Terminate
C     when the interval length drops to 1.

   10 CONTINUE
         LENGTH = RIGHT - LEFT
         IF (LENGTH .GT. 1) THEN

C           The trial value is a "linear" estimate of the left-hand
C           endpoint of the interval bracketing the target T, with
C           protection against round-off (which can affect convergence).

            TRIAL = MIN (RIGHT - 1, LEFT + MAX (0, INT (REAL (LENGTH) *
     >         (T - TLEFT) / (TRIGHT - TLEFT))))

C            ----------------------------------------------------------
C           |                                                          |
C           |   2 <= LEFT <= TRIAL < RIGHT <= N - 1                    |
C           |                                                          |
C            ----------------------------------------------------------

C           Compute the cumulative arclength up to TRIAL as cheaply as
C           possible using previous results. (The search runs equally well
C           for an increasing or decreasing sequence of T's.)

            IF ((TRIAL - LEFT) .LE. (RIGHT - TRIAL)) THEN
               TEMP1 = TLEFT + CHORD3D (X, Y, Z, LEFT, TRIAL)
            ELSE
               TEMP1 = TRIGHT - CHORD3D (X, Y, Z, TRIAL, RIGHT)
            END IF

C           Similar trick for arclength up to TRIAL + 1 since we may well
C           already know it, e.g., just before termination.

            IF (RIGHT .EQ. TRIAL + 1) THEN
               TEMP2 = TRIGHT
            ELSE
               TEMP2 = TEMP1 + CHORD3D (X, Y, Z, TRIAL, TRIAL + 1)
            END IF

C           Adjust the endpoints to reduce LENGTH as much as possible, but
C           not less than 1.

            IF (T .GE. TEMP2) THEN

C              Increase LEFT carefully to avoid overshoot (which can
C              occur due to roundoff error in chord calculations).

               IF (TRIAL + 1 .LT. RIGHT) THEN
                  LEFT   = TRIAL + 1
                  TLEFT  = TEMP2
               ELSE
                  LEFT   = RIGHT - 1
                  TLEFT  = TRIGHT - CHORD3D (X, Y, Z, LEFT, RIGHT)
               END IF
            ELSE IF (T .LT. TEMP1) THEN

C              Decrease RIGHT.

               IF (TRIAL .GT. LEFT) THEN
                  RIGHT  = TRIAL
                  TRIGHT = TEMP1
               ELSE
                  RIGHT  = LEFT + 1
                  TRIGHT = TLEFT + CHORD3D (X, Y, Z, LEFT, RIGHT)
               END IF
            ELSE

C              Adjust both LEFT and RIGHT. We're done since T is in the
C              interval [T (TRIAL), T (TRIAL+1)), but defer termination
C              until LENGTH is tested at the top of the loop.

               LEFT   = TRIAL
               TLEFT  = TEMP1
               RIGHT  = TRIAL + 1
               TRIGHT = TEMP2
            END IF
            GO TO 10

         END IF

C     Termination.
C     ------------

  990 CONTINUE
      RETURN
      END
