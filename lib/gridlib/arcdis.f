C+------------------------------------------------------------------------------
C
      SUBROUTINE ARCDIS (NEWDIS, QUERY, LUNCRT, LUNKBD,
     >                   DMETHOD, PMETHOD, MODE, IPARAM, RPARAM,
     >                   NWK, WK, CLOSED, NGEOM, I1GEOM, I2GEOM,
     >                   XGEOM, YGEOM, ZGEOM, TNORM,
     >                   NEVAL, XEVAL, YEVAL, ZEVAL, TEVAL, IER)
C
C  ONE-LINER:  Redistribution of points along (part of) an arc in 3-space.
C
C  PURPOSE:
C
C        ARCDIS is a high-level routine for applying the available 1-D grid
C     point distribution utilities and a choice of parametric fits to the
C     problem of redistributing points along arcs in 3-space.  It operates
C     on one (sub)arc per call but provides for applying similar distributions
C     efficiently to many (sub)arcs in a row (as in the case of surface grid
C     generation for the many stations of a wing).  It packages these
C     capabilities together with the appropriate data scaling that has been
C     found to be advisable.
C
C  METHOD:
C
C        ARCDIS wraps carefully scaled 3-space arclength calculations around
C     the 1-D distribution options offered by DISTRIB.  The user is referred
C     to DISTRIB for some of the argument descriptions because the two have
C     so much in common (and for maintenance reasons).
C
C        Provision has been made for repeated application of the current
C     distribution method without repetitious prompting.  Thus all prompts
C     are either enabled or disabled via argument QUERY, and relevant
C     parameters are assumed ready for (re)use if QUERY is false.  These
C     parameters may be determined at run time by a call to ARCDIS with
C     QUERY true, but they may also be hard-coded at compile time.
C
C        Sample application (run time control during first of many calls):
C
C        ::::::::::::::
C        MODE = 'R'            ! "Relative" means normalized Ts can be reused.
C        NEWDIS = .TRUE.       ! This argument is redundant if MODE = 'A'.
C        QUERY = .TRUE.        ! Prompting or not is independent of controlling
C                              ! generation of TNORM(*) on first call only.
C        DO xx K = 1, NSTATION
C           <set up the geometry data for this station>
C
C           CALL ARCDIS (NEWDIS, QUERY, ..., MODE,         )
C           NEWDIS = .FALSE.
C           QUERY = .FALSE.
C    xx  CONTINUE
C        ::::::::::::::
C
C        Some methods require additional work-space, but a given application
C     may have no intention of using certain methods.  Therefore the calling
C     program should supply only enough for the methods that will be selected.
C
C        Error handling is left to that of DISTRIB:  if a method fails, the
C     user is prompted as to whether to retry (go back to the menu), proceed
C     anyway, (since some iterative methods may still have given a usable
C     result), or quit.
C
C        PROGRAMMER NOTE:  Scaling and unscaling are performed in-place,
C     with care to preserve original end points, but some round-off may show
C     up in the interior points of the original X-Y-Z coordinates.  The WK(*)
C     argument can be used for temporary copies should this become a problem.
C
C  ARGUMENTS:
C     ARG      DIM   TYPE   I/O/S   DESCRIPTION
C     NEWDIS    -      L    I       Used to suppress recalculation of relative
C                                   distribution; ignored if MODE = 'A'.
C     QUERY     -      L    I       Used to suppress prompts.  Also:
C                                   QUERY=T: relevant choices for input/output
C                                            arguments are output.
C                                   QUERY=F: above choices are assumed input.
C     LUNCRT    -      I    I       Logical unit for prompts/diagnostics.
C     LUNKBD    -      I    I          "     "    "  keyboard entries.
C     DMETHOD   -      I    I/O     Distribution menu item number - see summary
C                                   in DISTRIB.
C     PMETHOD   -     C*1   I       Parametric fit method, presently limited
C                                   to those offered by PLSFIT3D:
C                                   'B' for BESSEL-type (loose fit) cubics;
C                                   'M' for MONOTONIC (tight fit) cubics;
C                                   'L' for piecewise linear.
C     MODE      -     C*1   I       'R' for RELATIVE means specifications for
C                                       certain distributions (such as the
C                                       initial increment) apply to the
C                                       NORMALIZED arclength, i.e. [0, 1].
C                                   'A' for ABSOLUTE means the specifications
C                                       apply to the UNNORMALIZED arclength.
C                                       Some methods (e.g. simple sinusoidal
C                                       distributions) are independent of MODE.
C     IPARAM    *      I    I/O     INTEGER parameter(s) used by some methods.
C                                   Usually no more than 1 - see DISTRIB.
C     RPARAM    *      R    I/O     REAL parameters used by certain methods.
C                                   Usually 1 to 3 values - see DISTRIB for
C                                   significance of -ve dX values (percentages).
C     NWK       -      I    I       Length of WK (*) needed.  NWK = 0 may be OK.
C                                   But NWK should be at least 1 for methods
C                                   using dX > 0. values in RPARAM (*).
C     WK        *      R        S   Work-space needed by some methods - see
C                                   DISTRIB.  Also see PROGRAMMER NOTE above.
C     CLOSED    -      L    I       CLOSED=.TRUE. means the input curve wraps
C                                   around smoothly on itself so periodic
C                                   boundary conditions are used by PLSFIT3D.
C                                   The first and last data points must agree.
C     NGEOM     -      I    I       Number of input data points.
C     I1GEOM,   -      I    I       Indices defining part of *GEOM (*) to be re-
C     I2GEOM                        distributed. 1 <= I1GEOM <= I2GEOM <= NGEOM.
C                                   (Curve fits may include data outside the
C                                   I1GEOM/I2GEOM range, for proper smoothness.
C                                   If an end-point is to treated as "sharp",
C                                   I1GEOM &/or I2GEOM should be 1 &/or NGEOM.)
C     XGEOM,  NGEOM    R    I       Input data coordinates.  See CLOSED.
C     YGEOM,
C     ZGEOM
C     TNORM   NGEOM    R    I/O     Normalized arclengths in the range [0, 1].
C                                   Applies only if MODE='R':
C                                   output (and used) if NEWDIS=.TRUE.;
C                                   input from a previous call if NEWDIS=.FALSE.
C     NEVAL     -      I    I/O     Number of interpolated points.
C     XEVAL,  NEVAL    R      O     Interpolated coordinates.
C     YEVAL,
C     ZEVAL
C     TEVAL   NEVAL    R        S   Workspace for the interpolated arclengths.
C     IER       -      I      O     IER = 0 means no error.
C                                   IER > 0 means some control parameter to
C                                           ARCDIS, PLSFIT3D, or DISTRIB was
C                                           bad.  Diagnostics are given.
C                                   IER =-1 means quit in DISTRIB (probably
C                                           abnormal in this application).
C
C  PROCEDURES:
C     CHORD3D   Chord-length utility.
C     DISTRIB   Drives the available 1-D distribution utilities.
C     GETSCALE, Scaling utilities. Their option used here transforms the coord.
C     USESCALE  with greatest range to [0,1], and the others to preserve shape.
C     PLSFIT3D  3-space parametric local spline, also used to find arc lengths.
C
C  FILES USED:
C     See LUNCRT and LUNKBD arguments.
C
C  ENVIRONMENT:
C     VAX/VMS; FORTRAN 77 with extensions:
C     >  IMPLICIT NONE
C     >  A few names longer than 6 characters
C     >  ! comments
C
C  HISTORY:
C     03/16/89  R.G.Langhi    Original REDIS utility had the main ideas of
C                             ARCDIS, but was not reusable.
C     08/09/89  D.A.Saunders  Wrapped 3-space arc redistribution around the
C                             DISTRIB utility (possibly more trouble than it
C                             was worth).
C     04/23/91   "    "       It WAS worthwhile!  (First application: a
C                             delta wing surface grid generator.)  However:
C                             grubby problem: any RPARAM (*) values repre-
C                             senting absolute distances (dX(1), etc.) have
C                             to be transformed as for the geometry data.
C                             (Relative distances are properly handled
C                             within DISTRIB.)  Therefore: pass the appropriate
C                             scale factor to DISTRIB via WK (1).
C     06/25/91   "    "       Introduced I1/I2GEOM for proper parametric
C                             interpolation of sub-curves (possibly rounded
C                             at I1GEOM or I2GEOM and therefore involving
C                             data points beyond the indicated range).
C                             Had to use CHORD3D directly rather than via
C                             PLSFIT3D if a sub-curve is being interpolated.
C
C  AUTHOR:   David Saunders, Sterling Software/NASA Ames, Palo Alto, CA.
C
C-------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER
     >   DMETHOD, IER, IPARAM (*), I1GEOM, I2GEOM, LUNCRT, LUNKBD,
     >   NEVAL, NGEOM, NWK
      REAL
     >   RPARAM (*), WK (*),
     >   TEVAL (NEVAL), TNORM (NGEOM), XEVAL (NEVAL), XGEOM (NGEOM),
     >   YEVAL (NEVAL), YGEOM (NGEOM), ZEVAL (NEVAL), ZGEOM (NGEOM)
      LOGICAL
     >   CLOSED, NEWDIS, QUERY
      CHARACTER
     >   MODE * 1, PMETHOD * 1

C     Local constants:

      INTEGER
     >   NDIM
      REAL
     >   ONE, ZERO
      CHARACTER
     >   DENORM * 1, DUNUSED * 1, GEOM * 1, NORM * 1, START * 1,
     >   SUBNAME * 6

      PARAMETER
     >  (DENORM  = 'D',          ! "Denormalize"
     >   DUNUSED = 'U',          ! Distribution method argument for PLSFIT3D,
                                 ! unused because its NEVAL arg. is always 1.
     >   GEOM    = 'G',          ! Geometric (shape-preserving) scaling
     >   NDIM    = 3,
     >   NORM    = 'N',          ! "Normalize"
     >   ONE     = 1.E+0,
     >   START   = 'A',          ! Auto-start for some iterative methods
     >   SUBNAME = 'ARCDIS',
     >   ZERO    = 0.E+0)

C     Local variables:

      INTEGER
     >   I, NEVALLOC
      REAL
     >   SCALE (3), SHIFT (3), TI1, TI2, XI1, YI1, ZI1, XI2, YI2, ZI2,
     >   X1, Y1, Z1, XN, YN, ZN
      LOGICAL
     >   NEWLOC
      CHARACTER
     >   DETAILS * 1      ! Dummy

C     Procedures:

      REAL
     >   CHORD3D
      EXTERNAL
     >   CHORD3D, DISTRIB, GETSCALE, USESCALE, PLSFIT3D

C     Execution:

C     Check for errors not handled by DISTRIB or PLSFIT3D:

      IER = 0
      IF (PMETHOD .NE. 'B' .AND.
     >    PMETHOD .NE. 'M' .AND.
     >    PMETHOD .NE. 'L')  GO TO 910
      IF (MODE .NE. 'A' .AND.
     >    MODE .NE. 'R')     GO TO 920
      IF (NEVAL .LE. 0)      GO TO 930

C     Scale the largest-spanning coordinate to [0, 1] and the others
C     in such a way that the shape is preserved:

      X1  = XGEOM (1)       ! Preserving end-points exactly is all we
      Y1  = YGEOM (1)       ! can reasonably do to avoid round-off during
      Z1  = ZGEOM (1)       ! in-place scaling/unscaling, short of making
      XI1 = XGEOM (I1GEOM)  ! temporary copies in WK (*).
      YI1 = YGEOM (I1GEOM)
      ZI1 = ZGEOM (I1GEOM)
      XI2 = XGEOM (I2GEOM)
      YI2 = YGEOM (I2GEOM)
      ZI2 = ZGEOM (I2GEOM)
      XN  = XGEOM (NGEOM)  
      YN  = YGEOM (NGEOM)
      ZN  = ZGEOM (NGEOM)

C     Have to work with the whole geometry because PLSFIT3D may go beyond
C     the range of I1/I2GEOM.

      CALL GETSCALE (GEOM, NDIM, NGEOM, XGEOM, YGEOM, ZGEOM,
     >               SCALE, SHIFT, IER)

      CALL USESCALE (NORM, NDIM, NGEOM, XGEOM, YGEOM, ZGEOM,
     >               SCALE, SHIFT, IER)

C     Determine the (scaled) sub-curve arc-length range.  If it really
C     is the whole curve, we should make use of PLSFIT3D's unavoidable
C     calculation of the total arc-length on its first call:

      NEVALLOC = 1
      NEWLOC = .TRUE.

      IF (I1GEOM .NE. 1 .OR. I2GEOM .NE. NGEOM) THEN
         TI1 = CHORD3D (XGEOM, YGEOM, ZGEOM, 1, I1GEOM)        ! Possibly zero.
         TI2 = CHORD3D (XGEOM, YGEOM, ZGEOM, I1GEOM, I2GEOM) + TI1
      ELSE
         TI1 = -ONE                          ! Tells PLSFIT3D to reset these
         TI2 = -ONE                          ! to zero and total arc-length

         CALL PLSFIT3D (NGEOM, XGEOM, YGEOM, ZGEOM, TI1, TI2,
     >                  NEVALLOC, XEVAL (I), YEVAL (I), ZEVAL (I),
     >                  NEWLOC, CLOSED, PMETHOD, DUNUSED, IER)
         IF (IER .NE. 0) GO TO 940

         NEWLOC = .FALSE.                    ! To avoid recalculating TI2.
      END IF

C     Generate the specified 1-D distribution of points along the arc.

C     Awkwardness:
C     Any RPARAM (*) values > 0. representing absolute distances need to be
C     scaled too.  Percentage distances (RPARAM (*) < 0.) are OK.
C     Best to confine this grubbiness to DISTRIB in case other methods are
C     added there.  Therefore, use its WK (*) argument to pass the appropriate
C     scale factor, in case scaling is needed.

      IF (MODE .EQ. 'A') THEN

C        "Absolute" case: no option to reuse a normalized distribution.
C        Generate the distribution TEVAL (*) on [TI1, TI2]:

         WK (1) = SCALE (1)   ! All SCALE (*) elements are the same here.

         CALL DISTRIB (DMETHOD, QUERY, NEVAL, TI1, TI2, LUNCRT,
     >                 LUNKBD, IPARAM, RPARAM, START, NWK, WK, TEVAL,
     >                 DETAILS, IER)
         IF (IER .NE. 0) GO TO 950

      ELSE

C        "Relative" distribution - may be able to reuse TNORM (*) on [0, 1]:

         IF (NEWDIS) THEN

C           Any non-relative distances (RPARAM (*) > 0) would have to
C           apply to [0, 1], so no scaling is required:

            WK (1) = ZERO
            CALL DISTRIB (DMETHOD, QUERY, NEVAL, ZERO, ONE, LUNCRT,
     >                    LUNKBD, IPARAM, RPARAM, START, NWK, WK, TNORM,
     >                    DETAILS, IER)
            IF (IER .NE. 0) GO TO 950

         END IF

         DO 300, I = 1, NEVAL
            TEVAL (I) = TI1 + TNORM (I) * (TI2 - TI1)
  300    CONTINUE

      END IF

C     Evaluate the interpolated curve coordinates one point at a time to
C     avoid PLSFIT3D's graphical prejudice towards internally-generated Ts:

      DO 400, I = 1, NEVAL
        CALL PLSFIT3D (NGEOM, XGEOM, YGEOM, ZGEOM, TEVAL (I), TEVAL (I),
     >                 NEVALLOC, XEVAL (I), YEVAL (I), ZEVAL (I),
     >                 NEWLOC, CLOSED, PMETHOD, DUNUSED, IER)
        IF (IER .NE. 0) GO TO 940
        IF (NEWLOC) NEWLOC = .FALSE.
  400 CONTINUE

C     Transform the interpolated values to the original units:

      CALL USESCALE (DENORM, NDIM, NEVAL, XEVAL, YEVAL, ZEVAL, SCALE,
     >               SHIFT, IER)

C     Restore the input data to the original units also:

      CALL USESCALE (DENORM, NDIM, NGEOM, XGEOM, YGEOM, ZGEOM, SCALE,
     >               SHIFT, IER)

      XGEOM (1) = X1
      YGEOM (1) = Y1
      ZGEOM (1) = Z1
      XGEOM (I1GEOM) = XI1
      YGEOM (I1GEOM) = YI1
      ZGEOM (I1GEOM) = ZI1
      XGEOM (I2GEOM) = XI2
      YGEOM (I2GEOM) = YI2
      ZGEOM (I2GEOM) = ZI2
      XGEOM (NGEOM)  = XN
      YGEOM (NGEOM)  = YN
      ZGEOM (NGEOM)  = ZN

      XEVAL (1) = XI1
      YEVAL (1) = YI1
      ZEVAL (1) = ZI1
      XEVAL (NEVAL) = XI2
      YEVAL (NEVAL) = YI2
      ZEVAL (NEVAL) = ZI2

      GO TO 999


C     Error handling:

  910 WRITE (LUNCRT, 1003) SUBNAME, ': Invalid PMETHOD: ', PMETHOD
      GO TO 990

  920 WRITE (LUNCRT, 1003) SUBNAME, ': Invalid MODE: ', MODE
      GO TO 990

  930 WRITE (LUNCRT, 1004) SUBNAME, ': Invalid NEVAL: ', NEVAL
      GO TO 990

  940 WRITE (LUNCRT, 1004)
     >   SUBNAME, ': Bad return from PLSFIT3D. IER: ', IER
      GO TO 990

  950 WRITE (LUNCRT, 1004)
     >   SUBNAME, ': Abnormal return from DISTRIB. IER: ', IER
      IF (IER .EQ. -1) GO TO 999   ! Quit; may want to retry with other params.

  990 IER = 1

  999 RETURN

C     Reusable formats:

 1003 FORMAT (A, A, A)
 1004 FORMAT (A, A, I6)

      END     ! End of ARCDIS
