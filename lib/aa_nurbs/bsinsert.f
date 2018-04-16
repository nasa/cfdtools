C+------------------------------------------------------------------------------
C
      SUBROUTINE BSINSERT (MAXC, CIN, TINSERT, MULT, LASTCPT, COUT, IER)
C
C     One-liner:
C
C        Inserts a knot with specified multiplicity for a B-spline curve
C
C     Description:
C
C           BSINSERT modifies a B-spline curve representation by inserting
C        the given knot with specified multiplicity and adjusting the control
C        points.  (Each singular insertion adds a control point and affects
C        neighboring control points; the curve shape is unchanged.)  If the
C        multiplicity equals the degree, then the last control point added
C        lies on the curve at the specified knot.  This last index is returned
C        to assist applications requiring fixed interior points on the curve.
C
C           Pre-existing knots of the same value are checked for so that the
C        final multiplicity cannot exceed the degree of the curve.
C
C           The output curve (in DT_NURBS form) may overwrite the input curve.
C        
C     Method:
C
C           To facilitate in-place insertion, only one at a time is
C        attempted via a lower level routine not intended for other uses.
C        (Minor local storage is still needed there.)  For each insertion,
C        p new control points k-p+1 : k are computed in place of p-1 old
C        ones (one more each time).  Here, the new knot is in the interval
C        [k, k+1) and p is the degree.
C
C           The basic algorithm is on p. 113 of the lecture notes by Tiller:
C
C           Q(i) = alpha(i) P(i) + (1 - alpha(i)) P(i-1)
C
C        where
C                    = 1                                   0 <= i <= k-p
C           alpha(i) = (TINSERT - T(i)) / (T(i+p) - T(i))  k-p < i < k+1
C                    = 0                                   i >= k+1
C
C     Environment:
C
C        FORTRAN 77 with minor extensions.
C
C     History:
C
C        25 Jan. 1994   DAS   Initial implementation, with fixing the leading
C                             edge of an airfoil in mind.
C
C     Author:  David Saunders, Sterling Software/NASA Ames, Mt. View, CA
C
C ------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER MAXC      ! Length of COUT (*)
      REAL    CIN (*)   ! Input B-spline curve in DT_NURBS form
      REAL    TINSERT   ! Knot value to be inserted
      INTEGER MULT      ! Number of repeated insertions.  If the given knot
     >                  ! already has multiplicity M, the final multiplicity
     >                  ! will be the smaller of M + MULT and degree P.
      INTEGER LASTCPT   ! Index of the last control point inserted (which is
     >                  ! on the curve if the final multiplicity is P)
      REAL    COUT (*)  ! Output curve, which may be passed as the same array
                        ! as the input curve
      INTEGER IER       ! IER =   0 means no error detected;
                        !     =  -1 means output C-array is too short
                        !     = -50 means TINSERT was out of range.
 
C-------------------------------------------------------------------------------

C     Local variables:

      INTEGER I, J, K, IX1, KORDER, M, NCPTS, NDIM, NKNOTS, NINSERT
      LOGICAL RATIONAL

C     Execution:

      NDIM   = INT (CIN (2))        ! 2 if X/Y/nonrational, -3 if rational, etc.
      KORDER = INT (CIN (3))        ! Degree + 1
      NCPTS  = INT (CIN (4))        ! # control pts.
      NKNOTS = NCPTS + KORDER       ! Initial # knots

      RATIONAL = NDIM .LT. 0
      IF (RATIONAL) THEN
         NDIM = -NDIM
         STOP 'BSINSERT: Rational not implemented yet.'
      END IF

C     Check that the output array is long enough.

      IF (MAXC .LT. 5 + KORDER + (NDIM + 1) * (NCPTS + MULT)) THEN
         IER = -1
         GO TO 999
      END IF

C     Locate the knot interval, K, in [KORDER, NKNOTS - KORDER]:

      CALL DTSPV2 (TINSERT, KORDER, CIN (6), 1, NKNOTS, 1, K, IER)
      IF (IER .NE. 0) GO TO 999

C     Check for any existing multiplicity, M.
C     DTSPV2 doesn't seem to return K = NKNOTS - NDEGREE if TINSERT
C     is equal to the last (repeated) knot.  Therefore, use a REAL test.

C*****IF (K .GT. NKNOTS - KORDER) THEN  ! Right-end special case
      IF (TINSERT .EQ. CIN (5 + NKNOTS)) THEN
         M = KORDER
      ELSE          ! This should be OK for left end as well as interior
         J = K
  200    IF (TINSERT .EQ. CIN (5 + J)) THEN
            J = J - 1
            IF (J .GT. 1) GO TO 200
         END IF
         M = K - J
      END IF

C     The (interior) multiplicity can be the degree at most:

      NINSERT = MAX (MIN (MULT, KORDER - 1 - M), 0)
      LASTCPT = NINSERT  ! 0 signals no change to the B-spline representation

C     Copy the whole thing in case NINSERT = 0 and so that we're
C     definitely updating in-place.

      DO I = 1, 5 + NKNOTS + NCPTS * NDIM
         COUT (I) = CIN (I)
      END DO

C     Insert just one knot at a time to minimize local storage, at the
C     expense of efficiency.  We need data structures!

      DO I = 1, NINSERT
         IX1 = 6 + NKNOTS     ! Index of first control pt. X in old curve

         CALL BSINSRT1 (NDIM, KORDER, NCPTS, NKNOTS, TINSERT, COUT (6),
     >                  K, COUT (IX1), COUT (IX1 + 1))

         NKNOTS = NKNOTS + 1
         NCPTS = NCPTS + 1
         COUT (4) = NCPTS  ! Inside the loop so LASTCPT stays 0 if NINSERT = 0
         LASTCPT = K       ! On the curve if multiplicity = degree
      END DO

C     Termination:

  999 RETURN
      END
C+------------------------------------------------------------------------------
C
      SUBROUTINE BSINSRT1 (NDIM, KORDER, NCPTS, NKNOTS, TINSERT, TKNOTS,
     >                     K, CIN, COUT)
C
C     Purpose:
C
C        Ancillary routine for BSINSERT, not intended to be used elsewhere.
C     Provides structure for the control point coordinates in the DT_NURBS
C     format, thus simplifying in-place insertion of one knot.
C
C-------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER NDIM      ! Number of dimensions (+ve here)
      INTEGER KORDER    ! Degree + 1
      INTEGER NCPTS     ! Input no. of control points
      INTEGER NKNOTS    !   "   "   "  knots
      REAL    TINSERT   ! Knot value being inserted as the (K+1)st
      REAL    TKNOTS (*)! Knots (before and after)
      INTEGER K         ! TINSERT lies in knot interval [K, K+1)
      REAL    CIN (NCPTS, NDIM)       ! Input control points
      REAL    COUT (NCPTS + 1, NDIM)  ! Output control points
                                      ! (same C array at higher level)

C     Local constants:

      REAL       ONE
      PARAMETER (ONE = 1.E+0)

C     Local variables:

      INTEGER I, IDIM, NDEG
      REAL    ALPHA

C     Execution:

      NDEG = KORDER - 1

C     Shift control point coordinates up, in-place.
C     Reverse order avoids overwriting.  Make room for one new control pt.

      DO IDIM = NDIM, 1, -1
         DO I = NCPTS, K - NDEG, -1
            COUT (I + 1, IDIM) = CIN (I, IDIM)
         END DO
         DO I = K - NDEG, 1, -1
            COUT (I, IDIM) = CIN (I, IDIM)
         END DO
      END DO

C     Apply the insertion algorithm to each coordinate:

      DO I = K - NDEG + 1, K
         ALPHA = (TINSERT - TKNOTS (I)) /
     >           (TKNOTS (I + NDEG) - TKNOTS (I))
         DO IDIM = 1, NDIM
            COUT (I, IDIM) = ALPHA * COUT (I + 1, IDIM) +
     >                       (ONE - ALPHA) * COUT (I, IDIM)
         END DO
      END DO

C     Update the knots.  Leave their number to be updated at the
C     higher level.  (Ditto for the number of control points.)

      DO I = NKNOTS, K + 1, -1
         TKNOTS (I + 1) = TKNOTS (I)
      END DO

      TKNOTS (K + 1) = TINSERT

      RETURN
      END
