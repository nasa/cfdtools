C-------------------------------------------------------------------------------
C
      PROGRAM CURVATURE3D
C
C     Purpose:
C
C        CURVATURE3D serves two purposes:  The first purpose is to estimate
C     the curvature, k(s), at each point along a given curve in 3-space
C     using parametric methods, where s is the arc length along the curve.
C     Four methods are provided for estimating the necessary derivatives:
C     the finite-differencing method, the conventional cubic spline,
C     the quasi-Hermite cubic spline, and the monotonic cubic spline.  Both
C     the quasi-Hermite and monotonic splines prove poor as their lack of
C     continuous 2nd derivatives would suggest, but they are retained for
C     illustrative purposes.  The results are saved in QPLOTable form:
C
C                  y vs. x,   z vs. x,   z vs. y,
C                  s vs. x,   s vs. y,   s vs. z,   or
C                  k vs. s,
C      
C     any of which may be suppressed.
C
C        The second purpose of CURVATURE3D is to serve as a driving program
C     for the subroutine CURVDIS3D, which is presented as an alternative to
C     estimating and plotting local curvatures.  This option redistributes
C     the points of the input curve with local curve spacing more or less
C     inversely proportional to the local curvature for any specified number
C     of output points.  Three scaling options are provided since curvature
C     magnitudes are sensitive to the units and CURVDIS3D expects data of
C     order 1.
C
C        This program began as a pedagogical exercise practicing Sterling
C     Software's programming style.  It has evolved into a convenient utility
C     for displaying curvature or redistributing points for (say) airframe
C     sections.
C
C     Method:
C
C        The program reads points from a file with a "SMOOTH2D" format:
C     [title and number of points in records 1 and 2, now optional, followed by]
C     columns of (x,y,z) data; any further columns are ignored.
C     
C        As the first step in displaying curvature, "PLFIT" computes the
C     cumulative chord lengths, which serve as the parametric variable, s.
C     Derivatives with respect to s may be estimated in several ways.  For
C     the conventional cubic spline option, "CSFIT" calculates the spline
C     coefficients - first for x vs. s, second for y vs. s, and then z vs. s.
C     After each "CSFIT" call, "CSDVAL" computes the 1st and 2nd derivatives
C     (dx/ds, d2x/ds2, dy/ds, d2y/ds2, dz/ds, and d2z/ds2) at each point.
C     The magnitude of the local curvature is then determined from:
C
C                    |k(s)| = | X' *  X" | / | X'| ** 3
C
C     where:       X' = (x', y', z') and X" = (x", y", z"),
C
C                     x' = dx / ds,     x" = d2x / ds2
C                     y' = dy / ds,     y" = d2y / ds2
C                     z' = dz / ds,     z' = d2z / ds2
C 
C
C     For the quasi-Hermite spline fit, the same steps are taken with "QHSFIT"
C     used instead of "CSFIT."  And for the monotonic cubic spline, "MSFIT" is
C     used in place of "CSFIT."  For the finite differencing method, values of
C     x vs. s, y vs. s and z vs. s are given to "FD12K", which generates the
C     derivatives directly.  Then k(s) is evaluated as before.
C
C        To redistribute points according to the local curvature, the user
C     may first normalize the input coordinates such that: 1) x, y and z are
C     in the range [0,1],  2) the coordinates are around the range [0,1], but
C     the geometric shape is retained, or 3) the coordinates are normalized
C     in precisely the same way as data from another geometric curve.  Or, the
C     user may elect not to rescale the data at all.   After the optional
C     rescaling, CURVDIS3D distributes the specified number of points along
C     the curve, with a fractional exponent provided to vary the clustering
C     effects.  Results are scaled back to the original units if necessary and
C     saved in plottable form.
C
C     External References:
C
C        ALPHA         Distinguishes numeric and alphanumeric input
C        CSFIT         Computes conventional cubic spline coefficients
C        CSDVAL        Evaluates 1-D cubic spline function and derivatives
C        CURVDIS3D     Redistributes points such that the local spacing
C                         is related to the local curvature
C        CURV3D        Modular form of the 3-space curvature formula
C        FD12K         Estimates derivatives through finite differencing
C        GETSCALE      Determines data normalization parameters
C        MSFIT         Computes monotonic cubic spline coefficients
C        OPENER        Prompts for and opens files
C        PLFIT         Calculates the chord lengths
C        QHSFIT        Computes quasi-Hermite spline coefficients
C        READER        Prompting utility
C        TOGGLE        Turns on and off a number of options
C        USESCALE      Normalizes or denormalizes 1, 2, or 3 dimensional data
C
C     History:
C
C        08/16/88   BAN/DAS   Initial design & implementation of CURVATURE
C                             for XY data
C        09/14/88   BAN       Installed subroutine CURVDIS
C        11/02/88   DAS       A few refinements in readiness for an XYZ analog
C        12/02/88   MDW       Modularized curve normalization procedures
C        12/08/88   MDW       Adapted CURVATURE3D from 2D program CURVATURE
C        08/21/91   PJT       Increased MXPTS from 350 to 4000.  First
C                             tabulation of redistributed points now contains
C                             all three coordinates.
C        08/27/91   DAS       Installed CURV3D in place of in-line code.
C        10/30/91   DAS       CURV3D now has the loop over N inside it.
C                             The PLFIT call had NDIM=2 where it should be 3.
C                             Replaced IQHSCU with QHSFIT.
C        09/02/08   DAS       Carriage control is a thing of the past now.
C                             Make the header records optional.
C        03/06/13   DAS       8 (not 7) significant digits match CAPSULE_GRID.
C        11/29/13   DAS       CURVDIS3D has been brought up to date with
C                             CURVDIS (smoothing options; vertex option).
C
C     Authors:  Michael Wong, NASA Ames/Sterling Software, Palo Alto, CA.
C               Brian Nishida
C               David Saunders   
C               
C-------------------------------------------------------------------------------

      IMPLICIT NONE


C     Local Constants:
C     ----------------

      INTEGER, PARAMETER ::
     >   LUNCRT = 6,
     >   LUNKBD = 5,
     >   LUNRD  = 1,                  ! Initial data & optional reference data
     >   LUNPLT = 2,                  ! QPLOTable results
     >   NDIM   = 3,                  ! Number of dimensions
     >   MXPTS  = 8000                ! Max. no. of points

      CHARACTER, PARAMETER :: BLANK*1 = ' '

      REAL, PARAMETER :: HALF = 0.5E+0, ONE = 1.E+0


C     Local Variables:
C     ----------------

      INTEGER    IENDL, IENDR, IER, IOS, ISMOOTH,
     >           METHOD, N, NNEW,
     >           NPTS, NPTREF, NSELEC
      REAL       B (MXPTS), C (MXPTS), D (MXPTS), ! Spline coefficents for CSFIT
     >           DX (MXPTS), DY (MXPTS),     ! 1st derivatives
     >           DZ (MXPTS),
     >           DXX (MXPTS), DYY (MXPTS),   ! 2nd derivatives
     >           DZZ (MXPTS),
     >           POWER,                      ! Controls clustering effects
     >           K (MXPTS),                  ! Local curvatures
     >           LX (MXPTS), LY (MXPTS),     ! Spline function evaluations (not
     >           LZ (MXPTS),                 !    used--just need derivatives)
     >           S (MXPTS),                  ! Data arc lengths
     >           SNEW (MXPTS),               ! Redistributed arc lengths
     >           STOTAL,
     >           XMAX, XMIN, YMAX, YMIN,
     >           ZMAX, ZMIN,
     >           X (MXPTS), Y (MXPTS),       ! x, y and z coordinates
     >           Z (MXPTS),
     >           XNEW (MXPTS), YNEW (MXPTS), ! Redistributed x, y and z values
     >           ZNEW (MXPTS),
     >           XREF (MXPTS), YREF (MXPTS), ! Reference curve data for 
     >           ZREF (MXPTS),               ! computing reference scaling 
                                             ! parameters
     >           SCALE (NDIM),               ! Scale factores and shifts
     >           SHIFT (NDIM)
      LOGICAL    ALPHA, CR, CYCLIC, EOF,
     >           ONOFF (7),                  ! Toggles for yes/no menu options
     >           ORIGINAL, DOSCALE, YES
      CHARACTER  PLOT (7) * 16,              ! Plotting options menu
     >           FILE1 * 60, FILE2 * 60,
     >           TITLE * 60, LEGEND * 60


C     Procedures:
C     -----------

      EXTERNAL
     >   ALPHA


C     Storage:
C     --------

      DATA FILE2
     >   /'curvature.plt'/
      DATA ONOFF
     >   /6 * .FALSE., .TRUE./
      DATA PLOT
     >   /'1: Plot Y vs. X?',
     >    '2: Plot Z vs. X?',
     >    '3: Plot Z vs. Y?',
     >    '4: Plot S vs. X?',
     >    '5: Plot S vs. Y?',
     >    '6: Plot S vs. Z?',
     >    '7: Plot K vs. S?'/


C     Execution:
C     ----------

      WRITE (LUNCRT, 1001)
      WRITE (LUNCRT, 1004)
     >   ' CURVATURE3D displays curvature properties for', 
     >   ' one (x,y,z) dataset, or',
     >   ' generates a new curvature-related point distribution.'


C     Open the input and output files

      CALL OPENER (LUNCRT, 'Enter input curve filename: ',
     >   LUNKBD, FILE1, LUNRD, 'OLD:FORMATTED')

      WRITE (LUNCRT, 1001) BLANK

      CALL OPENER (LUNCRT, 'Enter the filename of the ' //
     >   'QPLOTable results (<cr> = curvature.plt): ',
     >   LUNKBD, FILE2, LUNPLT, 'UNKNOWN:FORMATTED')

      WRITE (LUNCRT, 1001) BLANK

C     Read the input curve data.  (Separate READs ignore additional columns.)

      READ (LUNRD, 1001) TITLE

      IF (.NOT. ALPHA (TITLE)) THEN  ! Assume just numeric input
         TITLE = FILE1
         REWIND (LUNRD)
         NPTS = 0
         IOS = 0
         DO WHILE (IOS == 0) ! Read till EOF
            NPTS = NPTS + 1
            IF (NPTS .GT. MXPTS) GO TO 900
            READ (LUNRD, *, IOSTAT=IOS) X (NPTS), Y (NPTS), Z (NPTS)
         END DO
         NPTS = NPTS - 1
      ELSE
         READ (LUNRD, *) NPTS
         IF (NPTS .GT. MXPTS) GO TO 900
         DO N = 1, NPTS
            READ (LUNRD, *) X (N), Y (N), Z (N)
         END DO
      END IF

      CLOSE (LUNRD)
      YES = .FALSE.
      WRITE (LUNCRT, 1001) ' Would you rather redistribute points than'

      CALL READY (LUNCRT,
     >   'display curvature properties? <cr> = no: ',
     >   LUNKBD, YES, CR, EOF)

      IF (YES) GO TO 500
      WRITE (LUNCRT, 1001) BLANK
      CYCLIC = .FALSE.

      CALL READY (LUNCRT, 
     >   'Is the input curve cyclic (y/n)? <cr> = no: ',
     >   LUNKBD, CYCLIC, CR, EOF)

      WRITE (LUNCRT, 1001) BLANK

C     Compute the cumulative chord lengths:

      CALL PLFIT (3, NPTS, X, Y, Z, S, IER)
      IF (IER .NE. 0) GO TO 910


C     ***** Plotting Options *****

      CALL TOGGLE (7, PLOT, ONOFF, LUNCRT, LUNKBD)

      IF (ONOFF (1)) THEN
         WRITE (LUNPLT, 1001)
     >      TRIM (TITLE),
     >      'Input Curve: Y vs. X',
     >      'X Coordinate',
     >      'Y Coordinate'
         WRITE (LUNPLT, 1005) (X (N), Y (N), N = 1, NPTS)
         WRITE (LUNPLT, 1001) 'END FRAME'

      END IF


      IF (ONOFF (2)) THEN
         WRITE (LUNPLT, 1001)
     >      TRIM (TITLE),
     >      'Input Curve: Z vs. X',
     >      'X Coordinate', 
     >      'Z Coordinate'
         WRITE (LUNPLT, 1005) (X (N), Z (N), N = 1, NPTS)
         WRITE (LUNPLT, 1001) 'END FRAME'

      END IF


      IF (ONOFF (3)) THEN
         WRITE (LUNPLT, 1001)
     >      TRIM (TITLE),
     >      'Input Curve: Z vs. Y',
     >      'Y Coordinate',
     >      'Z Coordinate'
         WRITE (LUNPLT, 1005) (Y (N), Z (N), N = 1, NPTS)
         WRITE (LUNPLT, 1001) 'END FRAME'

      END IF


      IF (ONOFF (4)) THEN
         WRITE (LUNPLT, 1001)
     >      TRIM (TITLE),
     >      'X vs. S',
     >      'S, Cumulative Chord Length',
     >      'X Coordinate'
         WRITE (LUNPLT, 1005) (S (N), X (N), N = 1, NPTS)
         WRITE (LUNPLT, 1001) 'END FRAME'

      END IF


      IF (ONOFF (5)) THEN
         WRITE (LUNPLT, 1001)
     >      TRIM (TITLE),
     >      'Y vs. S',
     >      'S, Cumulative Chord Length',
     >      'Y Coordinate'
         WRITE (LUNPLT, 1005) (S (N), Y (N), N = 1, NPTS)
         WRITE (LUNPLT, 1001) 'END FRAME'

      END IF

      IF (ONOFF (6)) THEN
         WRITE (LUNPLT, 1001)
     >      TRIM (TITLE),
     >      'Z vs. S',
     >      'S, Cumulative Chord Length',
     >      'Z Coordinate'
         WRITE (LUNPLT, 1005) (S (N), Z (N), N = 1, NPTS)
         WRITE (LUNPLT, 1001) 'END FRAME'

      END IF
  
C     Curvature is normally requested, but . . .
      IF (.NOT. ONOFF (7)) GO TO 800
      


C     ***** Loop Over the Derivative Estimation Choices ****

      WRITE (LUNPLT, 1001)
     >   TRIM (TITLE),
     >   '[M6]k[M0] vs. S',
     >   'S, cumulative chord length',
     >   '[M6]k[M0], local curvature'
      WRITE (LUNCRT, 1001) ' ',
     >   ' Select a method of derivative estimation: '

  300 CONTINUE

      WRITE (LUNCRT, 1001)
     >   '        1) Finite differencing (best)',
     >   '        2) Conventional cubic spline',
     >   '        3) Quasi-Hermite spline',
     >   '        4) Monotonic local cubic spline'

      CR = .FALSE.
      CALL READI (LUNCRT, 'or <cr> to exit: ', LUNKBD, METHOD, CR, EOF)
      IF (CR .OR. EOF) GO TO 800

      WRITE (LUNCRT, 1001) BLANK


      IF (METHOD .EQ. 1) THEN

         LEGEND = 'Finite Differencing'

C        Calculate d2x/ds2
         CALL FD12K (NPTS, S, X, DX, DXX, DXX)

C        Calculate d2y/ds2
         CALL FD12K (NPTS, S, Y, DY, DYY, DYY)

C        Calculate d2z/dz2
         CALL FD12K (NPTS, S, Z, DZ, DZZ, DZZ)


      ELSE IF (METHOD .EQ. 2) THEN

         LEGEND = 'Standard Cubic Spline'

         IF (CYCLIC) THEN
            IENDL = 4
            IENDR = 4

         ELSE
            IENDL = 0
            IENDR = 0

         END IF

C        Calculate the spline coefficients for x(s)
         CALL CSFIT (NPTS, S, X, IENDL, 0, IENDR, 0, B, C, D, IER)
         IF (IER .NE .0) GO TO 920

C        Compute d2x/ds2
         CALL CSDVAL (NPTS, S, X, NPTS, S, B, C, D, LX, DX, DXX)

C        Calculate the spline coefficients for y(s)
         CALL CSFIT (NPTS, S, Y, IENDL, 0, IENDR, 0, B, C, D, IER)
         IF (IER .NE. 0) GO TO 930

C        Compute d2y/ds2
         CALL CSDVAL (NPTS, S, Y, NPTS, S, B, C, D, LY, DY, DYY)

C        Calculate the spline coefficients for z(s)
         CALL CSFIT (NPTS, S, Z, IENDL, 0, IENDR, 0, B, C, D, IER)
         IF (IER .NE. 0) GO TO 931

C        Compute d2z/dz2
         CALL CSDVAL (NPTS, S, Z, NPTS, S, B, C, D, LZ, DZ, DZZ)


      ELSE IF (METHOD .EQ. 3) THEN

         LEGEND = 'Quasi-Hermite Spline'

         WRITE (LUNCRT, 1003)

C        Calculate the spline coefficients for x(s)
         CALL QHSFIT (NPTS, S, X, 'P', B, C, D, IER)
         IF (IER .NE. 0) GO TO 932

C        Calculate d2x/ds2
         CALL CSDVAL (NPTS, S, X, NPTS, S, B, C, D, LX, DX, DXX)         

C        Calculate the spline coefficients for y(s)
         CALL QHSFIT (NPTS, S, Y, 'P', B, C, D, IER)
         IF (IER .NE. 0) GO TO 934

C        Calculate d2y/ds2
         CALL CSDVAL (NPTS, S, Y, NPTS, S, B, C, D, LY, DY, DYY)

C        Calculate the spline coefficients for z(s)
         CALL QHSFIT (NPTS, S, Z, 'P', B, C, D, IER)
         IF (IER .NE. 0) GO TO 936

C        Calculate d2z/ds2
         CALL CSDVAL (NPTS, S, Z, NPTS, S, B, C, D, LZ, DZ, DZZ)


      ELSE IF (METHOD .EQ. 4) THEN

         LEGEND = 'Monotonic Cubic Spline'

         WRITE (LUNCRT, 1003)

C        Find the spline coefficients for x(s)
         CALL MSFIT (NPTS, S, X, CYCLIC, B, C, D, IER)
         IF (IER .NE. 0) GO TO 940

C        Calculate d2x/ds2
         CALL CSDVAL (NPTS, S, X, NPTS, S, B, C, D, LX, DX, DXX)
 
C        Find the spline coefficients for y(s)
         CALL MSFIT (NPTS, S, Y, CYCLIC, B, C, D, IER)
         IF (IER .NE. 0) GO TO 950

C        Calculate d2y/ds2
         CALL CSDVAL (NPTS, S, Y, NPTS, S, B, C, D, LY, DY, DYY)

C        Find the spline coefficients for z(s)
         CALL MSFIT (NPTS, S, Z, CYCLIC, B, C, D, IER)
         IF (IER .NE. 0) GO TO 955

C        Calculate d2z/ds2
         CALL CSDVAL (NPTS, S, Z, NPTS, S, B, C, D, LZ, DZ, DZZ)


      ELSE
         WRITE (LUNCRT, 1001)
     >      ' That is not a valid choice. Please try again.', BLANK
         GO TO 300

      END IF


C     Compute the local curvatures:

      CALL CURV3D (NPTS, DX, DXX, DY, DYY, DZ, DZZ, K)


      WRITE (LUNPLT, 1001) ' $OPTIONS'
      WRITE (LUNPLT, 1002) 'LEGEND', 'Derivatives via ' // TRIM (LEGEND)
      WRITE (LUNPLT, 1001) ' $END'
      WRITE (LUNPLT, 1005) (S (N), K (N), N = 1, NPTS)
      WRITE (LUNPLT, 1001) 'END CURVE'
      WRITE (LUNCRT, 1001)
     >   ' You may repeat curvature calculations with another method:'
      GO TO 300



C     ***** Redistribute Points According to Local Curvatures *****

  500 CONTINUE
  
      DOSCALE = .TRUE.
      WRITE (LUNCRT, 1001) BLANK

      CALL READY (LUNCRT,
     >   'Would you like results to be plotted to scale? <cr> = yes: ',
     >   LUNKBD, DOSCALE, CR, EOF)

      ORIGINAL = .FALSE.
      WRITE (LUNCRT, 1001) BLANK

      CALL READY (LUNCRT,
     >   'Do you want to display the original data? <cr> = no: ',
     >   LUNKBD, ORIGINAL, CR, EOF)

      WRITE (LUNCRT, 1007)
     >   ' Number of points found in the input data: ', NPTS
      NNEW = NPTS
      WRITE (LUNCRT, 1001) 
     >   ' Enter the number of points desired for redistribution.'

      CALL READI (LUNCRT, '<cr> means the same number: ',
     >            LUNKBD, NNEW, CR, EOF)

      WRITE (LUNCRT, 1001) BLANK,
     >  ' CURVDIS3D can smooth the curvature-based shape fn. and/or',
     >  ' the redistributed arc lengths.  ISMOOTH options:',
     >  '    0: no smoothing',
     >  '    1: shape function only',
     >  '    2: new arc lengths only',
     >  '    3: both smoothings'

      ISMOOTH = 3
      CALL READI (LUNCRT,
     >    'Enter a smoothing choice [<cr> = 3 => both]: ',
     >    LUNKBD, ISMOOTH, CR, EOF)

      IF (ISMOOTH /= 0) THEN
         YES = .FALSE.
         CALL READY (LUNCRT,
     >      'Resolve any sharp vertices? [<cr> = no]: ',
     >      LUNKBD, YES, CR, EOF)
         IF (YES) ISMOOTH = -ISMOOTH
      END IF

C     Save the original dataset before any scaling is done on it:

      IF (ORIGINAL) THEN
         WRITE (LUNPLT, 1001)
     >      TRIM (TITLE),
     >      'Input coordinates',
     >      'X',
     >      'Y'

         IF (DOSCALE) 
     >      WRITE (LUNPLT, 1001)
     >      ' $OPTIONS',
     >      ' PLOT = ''SCALE'',',
     >      ' $END'

         WRITE (LUNPLT, 1009) (X (N), Y (N), N = 1, NPTS)
         WRITE (LUNPLT, 1001) 'END FRAME'

         WRITE (LUNPLT, 1001)
     >      TRIM (TITLE),
     >      'Input Coordinates',
     >      'X',
     >      'Z'
         WRITE (LUNPLT, 1009) (X (N), Z (N), N = 1, NPTS)  
         WRITE (LUNPLT, 1001) 'END FRAME'

         WRITE (LUNPLT, 1001)
     >      TRIM (TITLE),
     >      'Input Coordinates',
     >      'Y',
     >      'Z'
         WRITE (LUNPLT, 1009) (Y (N), Z (N), N = 1, NPTS)
         WRITE (LUNPLT, 1001) 'END FRAME'
    
      END IF


C     The curvature-based calculations by CURVDIS3D behave best if the
C     dataset is normalized to the region of the unit circle (roughly).
C     Since shape should be preserved for geometric data, there is more
C     than one possibility:

  520 NSELEC = 3

      WRITE (LUNCRT, 1001) ' ',
     >   ' Normalizing the data is recommended for best results.' //
     >   ' Options:',
     >   '   1) No scaling',
     >   '   2) X, Y, Z normalized to [0., 1.] independently',
     >   '   3) Geometric scaling (preserves shape)',
     >   '   4) Scaling reusing the parameters of another dataset'

      CALL READI (LUNCRT, 'Pick one [<cr> = 3 = preserve shape]: ',
     >           LUNKBD, NSELEC, CR, EOF)

      IF (NSELEC .LT. 1 .OR. NSELEC .GT. 4) THEN
         GO TO 520

      ELSE IF (NSELEC .EQ. 2) THEN   
         CALL GETSCALE ('N', NDIM, NPTS, X, Y, Z, SCALE, SHIFT, IER)
                     ! "Normal" or "independent" normalization

      ELSE IF (NSELEC .EQ. 3) THEN    
         CALL GETSCALE ('G', NDIM, NPTS, X, Y, Z, SCALE, SHIFT, IER)
                     ! Geometric normalization

      ELSE IF (NSELEC .EQ. 4) THEN

C        Read a reference dataset to be used for the scaling parameters:

         CALL OPENER (LUNCRT,
     >      ' Enter the filename of the data whose scaling' //
     >      ' parameters you wish to use: ', LUNKBD, FILE1,
     >      LUNRD, 'OLD, FORMATTED')

         READ (LUNRD, *)  
         READ (LUNRD, *) NPTREF
         DO 540 N=1, NPTREF
            READ (LUNRD, *) XREF (N), YREF (N), ZREF (N)
  540    CONTINUE

         CALL GETSCALE ('G', NDIM, NPTREF, XREF, YREF, ZREF, SCALE,
     >      SHIFT, IER)

      END IF   

      IF (NSELEC .NE. 1) THEN         ! Do the normalization
         CALL USESCALE ('N', NDIM, NPTS, X, Y, Z, SCALE, SHIFT, IER)
      END IF


C     Generate the curvature-based distribution:

      WRITE (LUNCRT, 1004) 
     >   ' Enter an exponent in the range [0., 1.] such as 0.5;',
     >   ' 0. => uniform spacing; 1. maximizes curvature effects.'


C     Start of possible-retry loop:

  610 POWER = 0.5E+0

      CALL READR (LUNCRT, 'Exponent choice [<cr> = 0.5; ^D = quit]: ',
     >           LUNKBD, POWER, CR, EOF)

      IF (EOF) GO TO 999
      WRITE (LUNCRT, 1001) BLANK

 
C     Redistribute points according to the local curvature:
C     To enable a look at the redistributed arc lengths, it is
C     convenient to do the calculation twice, with only the second
C     call completing the interpolation of X, Y, Z.

      CALL CURVDIS3D (NPTS, X, Y, Z, NNEW, POWER, ISMOOTH,
     >                LUNCRT, .TRUE., SNEW, SNEW, SNEW, IER)

      write (luncrt, '(a, i3)') 'CURVDIS3D (1) IER:', ier
      IF (IER .NE. 0) THEN
         WRITE (LUNCRT, 1006) 'CURVDIS3D', IER
         GO TO 610
      END IF

      STOTAL = SNEW (NNEW)   ! Normalize the redistributed arcs, for
      SNEW   = SNEW / STOTAL ! possible comparison with DISTRIBUTE results

      CALL CURVDIS3D (NPTS, X, Y, Z, NNEW, POWER, ISMOOTH,
     >                LUNCRT, .FALSE., XNEW, YNEW, ZNEW, IER)

      write (luncrt, '(a, i3)') 'CURVDIS3D (2) IER:', ier

C     Transform (denormalize) the new coordinates back to the original units:

      IF (NSELEC .NE. 1) THEN
         CALL USESCALE ('D', NDIM, NNEW, XNEW, YNEW, ZNEW, SCALE, SHIFT,
     >                  IER)
      END IF


C     Save results:

      WRITE (LUNPLT, 1001)
     >   TRIM (TITLE)
      WRITE (LUNPLT, 1008)
     >   'Redistributed coordinates.  Exponent:', POWER
      WRITE (LUNPLT, 1001)
     >   'X',
     >   'Y'
      IF (DOSCALE)
     >   WRITE (LUNPLT, 1001)
     >   ' $OPTIONS',
     >   ' PLOT = ''SCALE'',',
     >   ' $END'

C     The first tabulation contains all three coordinates, for possible
C     use in gridding applications (e.g. on a CAD system)...

      WRITE (LUNPLT, 1010) (XNEW (N), YNEW (N), ZNEW (N), N = 1, NNEW)
      WRITE (LUNPLT, 1001) 'END FRAME'

      WRITE (LUNPLT, 1001)
     >   TRIM (TITLE) 
      WRITE (LUNPLT, 1008)
     >   'Redistributed coordinates.  Exponent:', POWER
      WRITE (LUNPLT, 1001)
     >   'X',
     >   'Z'
      WRITE (LUNPLT, 1009) (XNEW (N), ZNEW (N), N = 1, NNEW)
      WRITE (LUNPLT, 1001) 'END FRAME'

      WRITE (LUNPLT, 1001)
     >   TRIM (TITLE)
      WRITE (LUNPLT, 1008)
     >   'Redistributed coordinates.  Exponent:', POWER
      WRITE (LUNPLT, 1001)
     >   'Y',
     >   'Z'
      WRITE (LUNPLT, 1009) (YNEW (N), ZNEW (N), N = 1, NNEW)
      WRITE (LUNPLT, 1001) 'END FRAME'

      WRITE (LUNPLT, 1001)
     >   TRIM (TITLE),
     >   'Normalized spacing distribution',
     >   'Snorm',
     >   'dSnorm'
      WRITE (LUNPLT, 1005) ((SNEW (N + 1) + SNEW (N)) * HALF,
     >                       SNEW (N + 1) - SNEW (N), N = 1, NNEW - 1)
      WRITE (LUNPLT, 1001) 'END FRAME'

C     ***** End Redistribution of Points *****

  800 WRITE (LUNCRT, 1004) '   QPLOTable results:  ', FILE2
      GO TO 999


C     Error Handling
C     --------------

  900 WRITE (LUNCRT, 1004) ' ',
     >   ' NPTS exceeds the maximum # points allowed for.  ',
     >   ' Please change the parameter "MXPTS."'
      GO TO 999

  910 WRITE (LUNCRT, 1006)  'PLFIT', IER
      GO TO 999

  920 WRITE (LUNCRT, 1006)  'CSFIT for X vs. S', IER
      GO TO 999

  930 WRITE (LUNCRT, 1006)  'CSFIT for Y vs. S', IER
      GO TO 999

  931 WRITE (LUNCRT, 1006)  'CSFIT for Z vs. S', IER
      GO TO 999

  932 WRITE (LUNCRT, 1006) 'QHSFIT for X vs. S', IER
      GO TO 999

  934 WRITE (LUNCRT, 1006) 'QHSFIT for Y vs. S', IER
      GO TO 999

  936 WRITE (LUNCRT, 1006) 'QHSFIT for Z vs. S', IER
      GO TO 999

  940 WRITE (LUNCRT, 1006)  'MSFIT for X vs. S', IER
      GO TO 999

  950 WRITE (LUNCRT, 1006)  'MSFIT for Y vs. S', IER
      GO TO 999

  955 WRITE (LUNCRT, 1006)  'MSFIT for Z vs. S', IER

  999 CONTINUE

C     Formats
C     -------

 1001 FORMAT (A)
 1002 FORMAT (1X, A, ' = ''', A, ''',' )
 1003 FORMAT (/, ' Warning: 2nd derivatives (and hence curvatures) ',
     >   'may be discontinuous.')
 1004 FORMAT (A, A)
 1005 FORMAT (2ES15.7)
 1006 FORMAT (//, ' Error return from ', A, ': IER = ', I3)
 1007 FORMAT (A, I3)
 1008 FORMAT (A, F8.5)
 1009 FORMAT (2ES15.7)
 1010 FORMAT (3ES15.7)
 1011 FORMAT (4ES15.7)

      END PROGRAM CURVATURE3D
