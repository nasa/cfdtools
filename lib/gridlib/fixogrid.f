C+----------------------------------------------------------------------
C
      SUBROUTINE FIXOGRID (JDIM, KDIM, JMAX, KMAX, KFINAL, Y, Z, W,
     >                     YC, ZC, A, B, DR1, DR2, METHOD)
C
C  ONE-LINER: Redistribute/stretch the "spokes" of one (half) O-grid.
C
C  PURPOSE:
C
C        FIXOGRID redistributes the radial lines of (half of) an O-grid
C     in two dimensions.  (The third dimension coordinates are presumably
C     determined by the body surface values, and are ignored here.)
C
C        The outer/upper boundary is made elliptical in the process, but
C     the inner/lower boundary is preserved.  The number of points along
C     each line may be changed if desired, and precise first and last
C     radial increments are established (see METHOD).
C
C        The input coordinates are overwritten upon return.  
C
C  METHOD/USAGE:
C
C        The input grid is assumed to be in reasonable shape, to the
C     extent that any spline extrapolation of radial lines to the
C     outer boundary will produce acceptable results.  (The application
C     may force LINEAR extrapolation if necessary, by overwriting the
C     last point of each radial line with one sufficiently far away that
C     this routine will always interpolate rather than extrapolate.
C     But extrapolation should normally be safe.)
C
C        With flow-solver applications in mind, the specified initial
C     increment DR1 is applied at the surface for all radial lines.
C     The last increment, on the other hand, varies to match the ellipse
C     defined by RA, RB, with DR2 assumed to apply to the major axis.
C
C        No assumption about the orientation of the input grid is
C     needed.  The ellipse is taken to be aligned with the coordinate
C     axes in standard form, and for each radial line, the point of
C     intersection with the ellipse is initially found by solving the
C     equation
C
C        F(T) = ((Y(T) - YC) / A) ** 2 + ((Z(T) - ZC) / B) ** 2 - 1 = 0
C
C     for T using local spline techniques.  The interval [0, T] is then
C     redistributed according to DR1, the J-dependent function of DR2,
C     and KFINAL using the method of Vinokur in subroutine HTDIS2.
C
C  ARGUMENTS:
C
C     ARG     DIM   TYPE I/O/S DESCRIPTION
C
C     JDIM,    -     I   I     Dimensions of input/output grid arrays.
C     KDIM
C     JMAX,    -     I   I     Y/Z (1:JMAX, 1:KMAX) contain the input
C                              coordinates.
C     KFINAL   -     I   I     Y/Z (1:JMAX, 1:KFINAL) contain the output
C                              coordinates.  The K = 1 line is preserved.
C     Y, (JDIM,KDIM) R   I/O   Grid coordinates as indicated above.
C     Z
C     W   (5*KDIM)   R      S  Work-space for Y,Z,T vectors along a J-line,
C                              and for redistributed Ts and Ys/Zs.
C     YC,      -     R   I     Center of ellipse defining the outer
C     ZC                       boundary to be imposed.
C     A,       -     R   I     Ellipse's Y/Z (major/minor) axis lengths.
C     B
C     DR1      -     R   I     First radial increment applied to each J line.
C     DR2      -     R   I     Last radial increment applied to major axis
C                              of ellipse, and adjusted according to the new
C                              outer boundary at each J if B and A differ.
C     METHOD        C*1  I     Type of piecewise fit to be used by the local
C                              method of LCSFIT for X = X(T) and Y = Y(T):
C                              'B' means "Bessel" (= loose fit, normally OK);
C                              'M' means Monotonic between data points;
C                              'L' means Linear.
C
C  PROCEDURES:
C     CHORD      Chord-length utility
C     ECOMPASS   Utility for identifying where a spoke meets the ellipse
C     GETROW     Extracts a row from a 2-D array
C     HTDIS2     1-D grid utility (Vinokur's method)
C     LCSFIT     Local cubic spline method
C     PUTROW     Inserts a row in a 2-D array
C
C  ENVIRONMENT:
C     VAX/VMS; FORTRAN 77 with ...
C     > IMPLICIT NONE
C     > Trailing ! comments
C     > Names up to 8 characters
C
C  HISTORY:
C     07/26/91  D.A.Saunders  Initial design and code, for HYPEROH.
C     08/01/91    "    "      Added METHOD argument: 'M' or 'L' may
C                             help difficult cases.
C     09/15/91    "    "      HTDIS2 has an additional argument now.
C
C  AUTHOR: David Saunders, Sterling Software/NASA Ames, Palo Alto, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments.

      INTEGER
     >   JDIM, KDIM, JMAX, KMAX, KFINAL
      REAL
     >   Y (JDIM, KDIM), Z (JDIM, KDIM), W (5 * KDIM), YC, ZC,
     >   A, B, DR1, DR2
      CHARACTER
     >   METHOD * 1

C     Local constants.

      REAL
     >   ONE, ZERO
      LOGICAL
     >   NEW

      PARAMETER
     >  (NEW    = .TRUE.,
     >   ONE    = 1.E+0,
     >   ZERO   = 0.E+0)

C     Local variables.

      INTEGER
     >   IER, IR, IT, IY, IZ, J, K, KINDEX
      REAL
     >   DR2J, RATIO, TC, TINT, YINT, ZINT

C     Procedures.

      EXTERNAL
     >   CHORD, ECOMPASS, GETROW, HTDIS2, PUTROW
      REAL
     >   CHORD


C     Execution.

      TC = B / A       ! "Thickness/chord" ratio of ellipse
      IY = KMAX + 1    ! Work-space pointers for input Y, Z along a J-line.
      IZ = KMAX + IY   ! (The initial Ts start in W (1).)
      IT = KMAX + IZ   ! Redistributed Ts
      IR = KFINAL + IT ! Results (reused for Y and Z)

      W (1) = ZERO

C     For each J-line...

      DO 500, J = 1, JMAX

C        Extract the J-line as a vector:

         CALL GETROW (JDIM, J, 1, KMAX, Y, W (IY))
         CALL GETROW (JDIM, J, 1, KMAX, Z, W (IZ))

C        Set up the original cumulative chord lengths:

         DO 200, K = 2, KMAX
            W (K) = CHORD (W (IY), W (IZ), K - 1, K) + W (K - 1)
  200    CONTINUE

C        Locate the intersection of this J-line with the ellipse.

         CALL ECOMPASS (YC, ZC, A, B, KMAX, W (IY), W (IZ), W (1),
     >                  YINT, ZINT, TINT, KINDEX, IER)

         IF (IER .NE. 0 .OR.           ! No intersection.
     >       KINDEX .EQ. 0) GO TO 910  ! Off the end to the "left" (T < 0).

C        Adjust the outer radial increment according to the position
C        of the boundary point along the major axis:

         RATIO = ABS (YINT - YC) / A     ! Would something else be better?
         DR2J = (RATIO + (ONE - RATIO) * TC) * DR2

C        Redistribute the point distances along the J-line.

         CALL HTDIS2 (.TRUE., ZERO, TINT, DR1, DR2J, KFINAL, W (IT), -6,
     >                IER)
         IF (IER .NE. 0) GO TO 920

C        Interpolate Y and Z at the redistributed T values, and overwrite
C        the input J-line.

C        (The original PLSFIT is appropriate here, but LCSFIT is used
C        already by ECOMPASS, so...)

         CALL LCSFIT (KMAX, W (1), W (IY), NEW, METHOD, KFINAL, W (IT),
     >                W (IR), W (IR))

         CALL PUTROW (JDIM, J, 1, KFINAL, Y, W (IR))

         CALL LCSFIT (KMAX, W (1), W (IZ), NEW, METHOD, KFINAL, W (IT),
     >                W (IR), W (IR))

         CALL PUTROW (JDIM, J, 1, KFINAL, Z, W (IR))

  500 CONTINUE

      RETURN

        
C     Error handling.

  910 STOP 'FIXOGRID: Bad return from ECOMPASS'

  920 STOP 'FIXOGRID: Bad return from HTDIS2'

      END
