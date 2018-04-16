C+----------------------------------------------------------------------
C
      SUBROUTINE PLSINT2D (MXXD, MXYD, NXD, NYD, I1D, I2D, XD, YD, ZD,
     >                     NLE, XLE, YLE, NTE, XTE, YTE, ROUND,
     >                     LEDGMTHD, TEDGMTHD, ROWMTHD, COLMTHD,
     >                     MXXE, MXYE, I1E, I2E, J1E, J2E, XE, YE, ZE)
C
C ACRONYM: Parametric Local Spline INTerpolation in 2 Dimensions
C          -          -     -      ---              - -
C PURPOSE:
C
C     PLSINT2D performs a specialized kind of parametric interpolation
C     for functions of two variables using 1-D local cubic splines.
C     It is the parametric form of LCSINT2D: the latter requires that
C     Z is strictly a single-valued function of X and Y, suited to
C     wing surface applications where the leading is sharp, whereas
C     PLSINT2D permits proper handling of rounded leading or trailing
C     edges by treating upper/lower surfaces together parametrically.
C
C     This version allows for "rows" to be independently round or
C     sharp at an interior I1D or I2D.
C
C     PLSINT2D should be called once per subarea with monotonic Xs
C     (and Ys), but the interpolation can still involve curve fits
C     which wrap around non-sharp edges of the subarea.
C
C     The initial application to wing surface grid generation determines
C     its restrictions:
C
C     > The input data consists of at least two "rows" (e.g. straight
C       sections of a wing) with constant Y within a given row, and the
C       same number of points along all rows (preferably the same relative
C       distribution).  The full surface should be contained in each
C       row (e.g. lower surface trailing to leading edge followed by the
C       upper surface leading to trailing edge).
C
C     > The parametric capability is in the X or I direction only.
C       (Providing it for the Y or J direction as well would not be
C       very meaningful given the way the surface is defined by rows.)
C
C     > Specifying the rows as "closed" curves (ROWMTHD = 'C') is not yet
C       supported: it corresponds to rounded trailing edges, which are
C       not envisaged for the wing surface grid generation application.
C       But there are a couple of ways such capability could be included.
C
C     > The points requiring evaluation are in pseudo-rectangular form
C       but are otherwise not assumed to have any structure.  (They may
C       have been generated as interior points via transfinite interpol-
C       ation from boundary points.)
C
C     > The "edges" of the region of interpolation may be defined by the
C       indicated end points of the sub-rows, but provision is made for
C       entering independent edge definitions as would be useful if an
C       edge had curved segments in the XY plane.  (We don't want to
C       require an input row at every edge point if the edge curvature
C       varies a lot: better to define the edge as densely as needed for
C       piecewise linear or cubic interpolation via LCSFIT, regardless
C       of the number of input rows.  This is a compromise, short of
C       providing for an external subroutine to define the edges to any
C       precision desired.)
C
C METHOD:
C
C     Use of local cubic splines avoids storing any spline coefficients.
C     First derivative continuity along interpolation lines is achieved;
C     second derivative continuity is NOT guaranteed, but for sufficiently
C     well-defined input data, this should not be a problem.
C
C     At most 4 points in a line define the spline in that neighborhood
C     (3 points at a boundary).  The units in all directions are assumed
C     to be reasonably scaled.  The steps:
C
C     > For each target point (XE, YE), determine J1 and J2 defining the
C       rows to use for splining across the rows in this neighborhood.
C       In general YE will lie in the middle interval defined by 4 rows
C       J1 : J2.
C
C     > Determine the RELATIVE location of XE along the (imaginary) row
C       through YE.  Local spline interpolation may be used to find the
C       end-points of this imaginary row, but strict linear interpolation
C       is an option as it may be more correct for some applications,
C       such as wings with cranked leading edges.
C
C     > For rows J = J1 : J2, interpolate at the same relative X.
C
C     > Use these 4 (or 3 or 2) interpolated points to interpolate at YE.
C
C ARGUMENTS:
C    ARG    DIM   TYPE I/O/S DESCRIPTION
C   MXXD,    -      I    I   Max. # points provided for in the data arrays.
C   MXYD                     (MXYD is superfluous but retained for symmetry.)
C   NXD      -      I    I   Actual number of points along each row.  NXD >= 2.
C   NYD      -      I    I   Actual number of rows.  NYD >= 2.
C   I1D,     -      I    I   Indices defining strictly monotonic Xs in the
C   I2D                      data.  1 <= I1D < I2D <= NXD.  Interpolations
C                            near an interior I1D or I2D can involve data
C                            beyond those indices, since they are done
C                            parametrically as needed for rounded edges.
C                            I1D, I2D must define true edges if NLE, NTE = 0.
C   XD  MXXD,MXYD   R    I   XD (1:NXD, 1:NYD) are the meaningful X coordinates
C                            supplied by rows; XD (I1D:I2D, 1:NYD) must
C                            have strictly increasing or decreasing Xs.
C   YD     MXYD     R    I   YD (1:NYD) are the Y coords. of the rows.
C   ZD  MXXD,MXYD   R    I   Z coordinates of the rows, as for XD (*, *).
C   NLE      -      I    I   Number of independent "leading" edge points
C                            supplied.  NLE=0 means the end-points of
C                            the rows are treated as defining the edge
C                            (at I1 or I2, whichever has smaller X).
C   XLE,    NLE     R    I   Leading edge definition (if NLE > 0).  (There is
C   YLE                      no convenient way to provide ZLE (*) as well.)
C   NTE      -      I    I   As above, for the "trailing" edge.  NTE=0 means
C   XTE,                     that the I=NXD end-points are used instead of
C   YTE,                     XTE (*) and YTE (*).
C   ROUND   NYD     L    I   ROUND (J) = .TRUE. causes interpolations along
C                            "row" J to be done parametrically.
C   LEDGMTHD -     C*1   I   Interpolation method to be used along the
C                            "leading" edge when determining the relative X
C                            position of each target point.  Options are
C                            those of the LCSFIT utility, q.v.:
C                            'L' means Linear;
C                            'M' means Monotonic local cubic spline;
C                            'B' means non-monotonic "Bessel"-type local
C                                cubic spline (looser fit).
C   LEDGMTHD -     C*1   I   Interpolation method to be used along the
C                            "trailing" edge, as above (LCSFIT control).
C   ROWMTHD  -     C*1   I   Method to be used for the interior interpolations
C                            along rows.  Options are those of PLSINTRP, with
C                            'B' (for Bessel-type loose fits) most suited
C                            to rounded edges.  'C' for closed-curve rows
C                            is not yet supported.
C   COLMTHD  -     C*1   I   Method to be used for the interior interpolations
C                            across rows, as above (LCSFIT control).
C   MXXE,    -      I    I   Max. # points provided for in the output arrays.
C   MXYE                     (MXYE is superfluous but retained for symmetry.)
C   I1E,     -      I    I   First/last indices in XE (*, *) and YE (*, *) of
C   I2E,                     points at which to interpolate.  (Two pointers
C   J1E,                     rather than 1 and NXE, say, allow singular points
C   J2E                      on the boundaries to be avoided.)  All of the
C                            target points should lie within the data region
C                            defined by I1D, I2D and 1, NYD.
C   XE  MXXE,MXYE   R    I   XE (I1E:I2E, J1E:J2E) are the meaningful X
C                            coordinates of the points of interpolation.
C   YE  MXYE,MXYE   R    I   Y coordinates as for XE.
C   ZE  MXXE,MXYE   R    O   Desired interpolated Z coordinates as for XE.
C                     
C USAGE NOTES:
C
C     >  Since LCSFIT and PLSINTRP perform no error trapping, nor does this
C        routine.  (Both may extrapolate without warning, but in the present
C        case, this is likely to be unintended.)
C
C     >  One "function" (Z) per call is assumed.  Since the local spline
C        method of computing a handful of spline coefficients on the fly
C        for each evaluation is reasonably efficient, applications with
C        more than one quantity per (X,Y) would not be TOO inefficient.
C
C     >  For strictly sharp edges, use LCSINT2D instead.
C
C EXTERNAL REFERENCES:
C
C     INTERVAL   1-D search utility
C     LCSFIT     1-D local cubic spline utility (Y vs. X)
C     PLSINTRP   Parametric form of LCSFIT (X vs. T, Y vs. T)
C
C ENVIRONMENT:  VAX/VMS FORTRAN 77
C
C HISTORY:
C   06/21/91   DAS   Initial adaptation of LCSINT2D for application to
C                    a double-delta wing with rounded leading edge.
C   06/30/91   DAS   Allowed for individual rows to be rounded or not,
C                    by adding the ROUND (*) argument.
C
C AUTHOR: David Saunders, Sterling Software, Palo Alto, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER
     >   MXXD, MXXE, MXYD, MXYE, NLE, NTE, NXD, NYD, I1D, I2D,
     >   I1E, I2E, J1E, J2E
      REAL
     >   XD (MXXD, MXYD), XE (MXXE, MXYE), XLE (*), XTE (*),
     >   YD (MXYD),       YE (MXXE, MXYE), YLE (*), YTE (*),
     >   ZD (MXXD, MXYD), ZE (MXXE, MXYE)
      LOGICAL
     >   ROUND (NYD)
      CHARACTER
     >   LEDGMTHD * 1, TEDGMTHD * 1, ROWMTHD * 1, COLMTHD * 1
     
C     Local constants:

      REAL
     >   ONE
      LOGICAL
     >   NEW
      PARAMETER
     >  (ONE = 1.E+0, NEW = .TRUE.)

C     Local variables:

      INTEGER
     >   I, ILE, ITE, J, JT, J1, J2, K, NFIT
      REAL
     >   F (4), RATIO, SIGNY, X, XT (4), X1, X2, Y

C     Procedures:

      EXTERNAL
     >   INTERVAL, LCSFIT, PLSINTRP

C     Execution:

      SIGNY = SIGN (ONE, YD (2) - YD (1))
      IF (XD (I1D, 1) .LT. XD (I2D, 1)) THEN
         ILE = I1D
         ITE = I2D
      ELSE
         ILE = I2D
         ITE = I1D
      END IF

C     For each target point ...

      DO 700, J = J1E, J2E
         DO 600, I = I1E, I2E
            X = XE (I, J)
            Y = YE (I, J)

C           Find the rows bracketing the point:

            J1 = 1
            CALL INTERVAL (NYD, YD, Y, SIGNY, J1)

C           LCSFIT uses at most 4 pts, and only 3 at a boundary.
C           NYD = 2 is also possible.

            J2 = MIN (J1 + 2, NYD)
            J1 = MAX (J1 - 1, 1)
            NFIT = J2 - J1 + 1

C           Determine the RELATIVE location of XE along the (imaginary) row
C           through Y.

            IF (NLE .GT. 0) THEN ! Interpolate the independent leading edge data

               CALL LCSFIT (NLE, YLE, XLE, NEW, LEDGMTHD, 1, Y, X1, X1)

            ELSE           ! Pick off end-pts. of the rows to define the edge.
               JT = J1
               DO 300, K = 1, NFIT
                  F (K) = XD (ILE, JT)
                  JT = JT + 1
  300          CONTINUE

               CALL LCSFIT (NFIT, YD (J1), F, NEW, LEDGMTHD, 1, Y, X1,
     >                      X1)
            END IF

C           ... and now the opposite edge:

            IF (NTE .GT. 0) THEN

               CALL LCSFIT (NTE, YTE, XTE, NEW, TEDGMTHD, 1, Y, X2, X2)

            ELSE
               JT = J1
               DO 400, K = 1, NFIT
                  F (K) = XD (ITE, JT)
                  JT = JT + 1
  400          CONTINUE

               CALL LCSFIT (NFIT, YD (J1), F, NEW, TEDGMTHD, 1, Y, X2,
     >                      X2)
            END IF

            RATIO = (X - X1) / (X2 - X1)
            
C           For rows J = J1 : J2, interpolate at the same relative X.
C           Do it parametrically if the row is rounded:

            K = 1
            DO 500, JT = J1, J2
               XT (K) = (XD (ITE, JT) - XD (ILE, JT)) * RATIO +
     >                   XD (ILE, JT)

               IF (ROUND (JT)) THEN

                  CALL PLSINTRP (NXD, XD (1, JT), ZD (1, JT), I1D, I2D,
     >                           ROWMTHD, 1, XT (K), F (K), F (K))
               ELSE

                  CALL LCSFIT (I2D-I1D+1, XD (I1D, JT), ZD (I1D, JT),
     >                         NEW, ROWMTHD, 1, XT (K), F (K), F (K))
               END IF

               K = K + 1
  500       CONTINUE

C           Use these 4 (or 3 or 2) interpolated points to interpolate at Y:

            CALL LCSFIT (NFIT, YD (J1), F, NEW, COLMTHD, 1, Y,
     >                   ZE (I, J), ZE (I, J))
  600    CONTINUE
  700 CONTINUE
      
      RETURN         ! End of LCSINT2D
      END
