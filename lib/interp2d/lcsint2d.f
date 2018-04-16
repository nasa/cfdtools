C+----------------------------------------------------------------------
C
      SUBROUTINE LCSINT2D (MXXD, MXYD, NXD, NYD, XD, YD, ZD,
     >                     NLE, XLE, YLE, NTE, XTE, YTE,
     >                     LEDGMTHD, TEDGMTHD, ROWMTHD, COLMTHD,
     >                     MXXE, MXYE, I1E, I2E, J1E, J2E, XE, YE, ZE)
C
C ACRONYM: Local Cubic Spline INTerpolation in 2 Dimensions
C          -     -     -      ---              - -
C PURPOSE:
C
C     LCSINT2D performs a specialized kind of interpolation in 2-space
C     using 1-D local cubic spline methods.  Its initial application to
C     wing surface grid generation determines its restrictions:
C
C     > The input data consists of at least two "rows" (e.g. straight
C       sections of a wing) with constant Y within a given row, and the
C       same number of points along all rows (preferably the same relative
C       distribution).
C
C     > The points requiring evaluation are in pseudo-rectangular form
C       but are otherwise not assumed to have any structure.  (They
C       may have been generated as interior points via transfinite
C       interpolation from boundary points.)
C
C     > The "edges" of the region of interpolation may be defined by the
C       end points of the input rows, but provision is made for entering
C       independent edge definitions as would be useful if an edge had
C       curved segments.  (We don't want to require an input row at every
C       edge point if the edge curvature varies a lot: better to define
C       the edge as densely as needed for piecewise linear or cubic
C       interpolation via LCSFIT, regardless of the number of input rows.
C       This is a compromise, short of providing for an external subroutine
C       to define the edges to any precision desired.)
C
C METHOD:
C
C     Use of local cubic splines avoids storing spline coefficients for
C     each row.  First derivative continuity along interpolation lines is
C     guaranteed; second derivative continuity is NOT guaranteed, but for
C     sufficiently well-defined input data, this should not be a problem.
C     At most 4 points in a line define the spline in that neighborhood
C     (3 points at a boundary).  The units in all directions are assumed
C     to be reasonably scaled.  The steps:
C
C     > For each target point (XE, YE), determine J1 and J2 defining the
C       rows to use for splining across the rows in this neighborhood.
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
C   XD  MXXD,MXYD   R    I   XD (1:NXD, 1:NYD) are the meaningful X coordinates
C                            supplied by rows.
C   YD     MXYD     R    I   YD (1:NYD) are the Y coords. of the rows.
C   ZD  MXXD,MXYD   R    I   Z coordinates of the rows, as for XD (*, *).
C   NLE      -      I    I   Number of independent "leading" edge points
C                            supplied.  NLE=0 means the I=1 end-points of
C                            the rows are treated as defining the edge.
C   XLE,    NLE     R    I   Leading edge definition (if NLE > 0).  (There is
C   YLE                      no convenient way to provide ZLE (*) as well.)
C   NTE      -      I    I   As above, for the "trailing" edge.  NTE=0 means
C   XTE,                     that the I=NXD end-points are used instead of
C   YTE,                     XTE (*) and YTE (*).
C   LEDGMTHD -     C*1   I   Interpolation method to be used along the
C                            "leading" edge when determining the relative X
C                            position of each target point.  Options are
C                            those of the LCSFIT utility, q.v.:
C                            'L' means Linear;
C                            'M' means Monotonic local cubic spline;
C                            'B' means non-monotonic "Bessel"-type local
C                                cubic spline (looser fit).
C   LEDGMTHD -     C*1   I   Interpolation method to be used along the
C                            "trailing" edge, as above.
C   ROWMTHD  -     C*1   I   Method to be used for the interior interpolations
C                            along rows, as above.
C   COLMTHD  -     C*1   I   Method to be used for the interior interpolations
C                            across rows, as above.
C   MXXE,    -      I    I   Max. # points provided for in the output arrays.
C   MXYE                     (MXYE is superfluous but retained for symmetry.)
C   I1E,     -      I    I   First/last indices in XE (*, *) and YE (*, *) of
C   I2E,                     points at which to interpolate.  (Two pointers
C   J1E,                     rather than 1 and NXE, say, allow singular points
C   J2E                      on the boundaries to be avoided.)
C   XE  MXXE,MXYE   R    I   XE (I1E:I2E, J1E:J2E) are the meaningful X
C                            coordinates of the points of interpolation.
C   YE  MXYE,MXYE   R    I   Y coordinates as for XE.
C   ZE  MXXE,MXYE   R    O   Desired interpolated Z coordinates as for XE.
C                     
C USAGE NOTES:
C
C     >  Since LCSFIT performs no error trapping, nor does this routine.
C        (Both may extrapolate without warning, but in the present case,
C        this is likely to be unintended.)
C
C     >  One "function" (Z) per call is assumed.  Since the local spline
C        method of computing a handful of spline coefficients on the fly
C        for each evaluation is reasonably efficient, applications with
C        more than one quantity per (X,Y) would not be TOO inefficient.
C
C     >  Application to wings with rounded leading edges is not intended.
C        (Proper handling of the leading edge region would involve treating
C        both surfaces as one, and using parametric interpolation.)
C
C EXTERNAL REFERENCES:
C     INTERVAL   1-D search utility
C     LCSFIT     1-D local cubic spline utility
C
C ENVIRONMENT:  VAX/VMS FORTRAN 77
C
C HISTORY:
C   04/08/91   DAS   Initial implementation for a double-delta wing
C                    surface grid application.
C   04/15/91    "    Introduced the separate leading/trailing edge
C                    definition.  (John Melton had proposed an external
C                    routine; Jeff Trosin suggested the discrete point
C                    compromise.)
C   04/17/91    "    Changed arguments NXE, NYE, to I1E, I2E, etc., to
C                    overcome a singular point problem in the evaluation
C                    arrays.  Introduced ROWMTHD, COLMTHD in place of
C                    a common SURFMTHD.
C   04/18/91    "    Introduced LEDGMTHD, TEDGMTHD in place of EDGEMTHD.
C
C AUTHOR: David Saunders, Sterling Software, Palo Alto, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER
     >   MXXD, MXXE, MXYD, MXYE, NLE, NTE, NXD, NYD, I1E, I2E, J1E, J2E
      REAL
     >   XD (MXXD, MXYD), XE (MXXE, MXYE), XLE (*), XTE (*),
     >   YD (MXYD),       YE (MXXE, MXYE), YLE (*), YTE (*),
     >   ZD (MXXD, MXYD), ZE (MXXE, MXYE)
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
     >   I, J, JT, J1, J2, K, NFIT
      REAL
     >   F (4), RATIO, SIGNY, X, XT (4), X1, X2, Y

C     Procedures:

      EXTERNAL
     >   INTERVAL, LCSFIT

C     Execution:

      SIGNY = SIGN (ONE, YD (2) - YD (1))

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
                  F (K) = XD (1, JT)
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
                  F (K) = XD (NXD, JT)
                  JT = JT + 1
  400          CONTINUE

               CALL LCSFIT (NFIT, YD (J1), F, NEW, TEDGMTHD, 1, Y, X2,
     >                      X2)
            END IF

            RATIO = (X - X1) / (X2 - X1)
            
C           For rows J = J1 : J2, interpolate at the same relative X:

            K = 1
            DO 500, JT = J1, J2
               XT (K) = XD (1, JT) + RATIO * (XD (NXD, JT) - XD (1, JT))

               CALL LCSFIT (NXD, XD (1, JT), ZD (1, JT), NEW,
     >                      ROWMTHD, 1, XT (K), F (K), F (K))
               K = K + 1
  500       CONTINUE

C           Use these 4 (or 3 or 2) interpolated points to interpolate at Y:

            CALL LCSFIT (NFIT, YD (J1), F, NEW, COLMTHD, 1, Y,
     >                   ZE (I, J), ZE (I, J))
  600    CONTINUE
  700 CONTINUE
      
      RETURN         ! End of LCSINT2D
      END
