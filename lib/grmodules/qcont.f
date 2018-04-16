C+**********************************************************************
C
      SUBROUTINE QCONT ( NDIM, IDIM, FUNCT, FMIN, FMAX, XGRID,
     >                   YGRID, XMIN, XMAX, YMIN, YMAX, IMIN,
     >                   IMAX, JMIN, JMAX, NLEVS, RLEVS, LINTYP,
     >                   LINTHK, MXCOLR, LINCOL, LPRLEV, ILEGND,
     >                   HTLGND, ANGLBL, MAXPTS, MXNLEV, RWORK, 
     >                   MXLSEG, IWORK, LUNERR )
C
C ACRONYM: Quick CONTour routine
C          -     ----
C PURPOSE:
C    QCONT is a general-purpose routine for identifying and plotting
C    all requested contour levels within a given region. The function
C    values within the region must be represented as a (pseudo-)rect-
C    angular array.  The corresponding coordinate arrays may be 1- or
C    2-dimensional, according to whether the region is truly or only
C    pseudo-rectangular.
C
C    Provision is also made for direct labeling of contour level
C    indices.
C
C    QCONT does not handle plotting of the boundary of the contoured
C    region, or of axes, etc.   See modules PLTGEO, POLYLINE, etc., for
C    assembling complete composite plots.
C
C METHOD:
C    Assume the following prior set-up for this (sub)plot:
C      *  physical origin is set;
C      *  axes are scaled and drawn.
C
C    Set character height and angle for labels on contour lines.
C    If requested, get contour levels by calling CONSCL.
C    If requested, set default line type and line color values.
C    FOR each requested contour level DO
C       Obtain contour data by calling CONT2D;
C       IF there is at least 1 contour line THEN
C          FOR each contour line at this level DO
C             Plot line;
C             Plot level number level at mid-pt. of line.
C
C NOTES:
C    QCONT requires MXNLEV words of storage for the arrays RLEVS, 
C    LINTYP, LINTHK and LINCOL because the number of contour levels
C    will not always be known in advance.  QCONT will issue warning
C    messages and return partial plot data if work space is inadequate. 
C    The legend subroutine is called outside QCONT to avoid its use
C    in a 2D plot projected in 3D space.
C
C ARGUMENTS:
C     ARG   TYPE I/O/S  DIM     DESCRIPTION
C    NDIM     I    I     -      No. of dimensions of XGRID and YGRID.
C    IDIM     I    I     -      Effective row dimension of FUNCT(*,*)
C                               in the routine setting up this array.
C    FUNCT    R    I IDIM,JMAX  Function values to be contoured are in
C                               FUNCT(IMIN:IMAX,JMIN:JMAX).
C    FMIN,    R    I     -      FMIN and FMAX are arbitrary maximum
C     FMAX                      and minimum function data levels. They
C                               are used only when NLEVS <=0 requesting
C                               automatic selection of contour levels.
C    XGRID,   R    I  IDIM,*    The coords. matching FUNCT(I,J) are:
C    YGRID                      (  XGRID(I), YGRID(J) ) if NDIM=1;
C                               (XGRID(I,J),YGRID(I,J)) if NDIM=2.
C                               (* = 1 or JMAX.)
C  XMIN,XMAX, R    I     -      Data range of plot axes to enable 
C   YMIN,YMAX                   clipping of contour labels. DISSPLA
C                               handles clipping of curves.
C    IMIN,IMAX,          -      Define a window within the FUNCT array
C    JMIN,JMAX                  for the portion to be contoured.
C    NLEVS    I   I/O    -      NLEVS>0 means plot NLEVS user-specified
C                               contour levels.
C                               NLEVS=0 means plot contour levels from an
C                               array with increment RLEVS(1) and maximum
C                               and minimum from requested data window.
C                               NLEVS<0 means means plot approximately
C                               ABS(NLEVS) nicely scaled contour levels.
C                               On output NLEVS is the number of levels
C                               plotted.
C    RLEVS    R   I/S  MXNLEV   Array of contour levels. Calling program
C                               must provide at least MXNLEV words even if
C                               the levels are not user-specified. RLEVS
C                               is not modified on return.
C                               RLEVS(1) = desired contour level increment
C                               if NLEVS=0.
C    LINTYP   I    I   MXNLEV   The line type corresponding to the current 
C                               contour level.  (Note: Use of POLYLINE
C                               necessitates that the 1st contour level be 
C                               of line type 3.)
C                               Contour level      Line Type      Pattern
C                               -------------      ---------      -------
C                                     1              3             Solid
C                                     2              4             Dotted
C                                     3              5             Dashed
C                                     4              6             Chain-dotted
C                                     5              7             Chain-dashed
C                                     6              8             Long-dashed
C                               LINTYP(1)=99 means plot all contour
C                               lines solid with default thickness and 
C                               overrides any values of LINTHK.
C    LINTHK   I    I   MXNLEV   Defines the thickness of the line drawn
C                               in multiples of .01 inch. 1:9 is valid.
C    LINCOL   R    I   MXNLEV   LINCOL(I) specifies the color of the
C                               Ith contour line  plotted,  but if
C                               LINCOL(1)=99.,  QCONT will set a default
C                               color range from blue to magenta with
C                               increasing level.
C                               Any other number outside the range 1.- 8.
C                               will be ignored by POLYLINE, which will
C                               use the output device default color.
C                                  1.   ...  Black
C                                  2.   ...  Magenta
C                                  3.   ...  Red
C                                  4.   ...  Yellow
C                                  5.   ...  Green
C                                  6.   ...  Cyan
C                                  7.   ...  Blue
C                                  8.   ...  White
C                               with colors between the six primaries
C                               specified using real numbers in [2.,7.]
C                               (E.g. orange = 3.5, yellow-green = 4.5, etc.)
C    LPRLEV   I    O   MXNLEV   LPRLEV returns the number of line segments
C                               found per contour level as required by
C                               the legend drawing routine.
C    ILEGND   I    I     -      ILEGND=0 means lines are labeled with
C                               level numbers. (A legend may be provided
C                               outside this routine.)
C    ANGLBL   R    I     -      Angle (degrees) for contour labels.
C    MAXPTS   I    I     -      Used for dimensions of the floating-
C                               point work-space arrays needed by QCONT
C                               and CONT2D.  (Try MAXPTS = 1000.) 
C    MXNLEV    I    I     -     MXNLEV is used by QCONT as the dimension of 
C                               RLEVS, LINCOL, LINTHK, and LINTYP 
C                               and as the maximum number of contour levels
C                               plotted.
C    RWORK    R    S  MAXPTS*4  RWORK(1:2*MAXPTS) holds the abscissas and
C                               ordinates of contour lines returned by
C                               CONT2D;
C                               RWORK(2*MAXPTS+1:3*MAXPTS) is scratch for
C                               CONT2D;
C                               RWORK(3*MAXPTS+1:4*MAXPTS) is scratch for
C                               QCONT.
C    MXLSEG    I    I    -      Dimension of integer work-space arrays
C                               needed by CONT2D involving expected max-
C                               imum number of separate lines for  any
C                               given contour level. 
C    IWORK    I    S  MXLSEG    IWORK(1) gives the number of contour
C                               lines found in one call to CONT2D;
C                               IWORK(I) for I=2,MXLSEG points to the
C                               start of the Ith contour line returned
C                               in RWORK (reused for each call to 
C                               CONT2D).
C   LUNERR     I    I    -      Logical unit for error messages.
C
C
C EXTERNAL REFERENCES:
C    CONT2D    Prepares data for contour plotting.
C    CONSCL    Finds "nice" scaling of a number of values over a given
C              data range.
C    POLYLINE  Plots curve of indicated type/thickness/color.
C    CLRNGE    Sets color indices for contour lines.
C    Certain DISSPLA routines for line labels, legend and color.
C
C KNOWN MACHINE DEPENDENCIES:  DISSPLA graphics software. Color output
C    depends on plotting hardware. Requests for color are ignored when
C    color is  not available.   Interpretation  of  color selection is
C    device dependent.
C
C STANDARDS VIOLATIONS: None.
C
C ENVIRONMENT:  VAX/VMS, FORTRAN 77
C
C HISTORY:
C  31 Aug.1987  RCL       Added argument MXLSEG to simplify changes in
C                         dimension of the pointer array for individual
C                         line segments for each contour level. Control
C                         is now entirley in the calling routine.
C  24 Mar.1987  RCL       Color range now calculated with real numbers to
C                         allow specification of intermediate colors.
C  3 Sept.1986  RCL       This is a new version of the original PLLEVS,
C                         with revised handling of line type and thick-
C                         ness and an option for user input of contour
C                         level max and min.
C  6 Dec. 1989  MDW       Replaced call to FMCURV with call to curve utility
C                         POLYLINE.  
C 18 Jan. 1990  MDW       Added LUNERR to argument list.  Some minor changes to
C                         array declarations.
C
C AUTHORS: Jeff Trosin, Carol Green, Ron Langhi, David Saunders,
C          Rosalie Lefkowitz, Informatics General Corp., Palo Alto, CA.
C
C-**********************************************************************

      IMPLICIT NONE

C  *  Arguments.

      INTEGER    IDIM, ILEGND, IMAX, IMIN, JMAX, JMIN, LUNERR, MAXPTS, 
     >           MXCOLOR, MXLSEG, MXNLEV, NDIM, NLEVS
      INTEGER
     >           LINTYP (*) , LINTHK (*), LPRLEV (MXNLEV), 
     >           MXCOLR, IWORK (MXLSEG)

      REAL       ANGLBL, FMAX, FMIN, HTLGND, XMAX, XMIN, YMAX, YMIN

      REAL       FUNCT(IDIM,JMAX), XGRID(IDIM,*), YGRID(IDIM,*),
     >           RLEVS (*), LINCOL (*), RWORK(MAXPTS*4)

C  *  Local variables.

      INTEGER    I, I1, IADIM, ICOLRY, IER, IIA, IX, IXCON, IY, IYCON,
     >           LEV, LINE, LTYP, MXCONT, NCONT, NLINES, NPTS, THK

      REAL       ACONT, CINC, COLOR, XLBL, YLBL

C  *  Execution.

      CALL HEIGHT ( HTLGND )
      CALL ANGLE  ( ANGLBL )
      IXCON  = 1
      IYCON  = MAXPTS + 1
      IIA    = MAXPTS + IYCON
      ICOLRY = MAXPTS + IIA
      IADIM  = MAXPTS
      MXCONT = MIN(MXNLEV,MAXPTS)
      NCONT  = ABS(NLEVS)
      IF ( NCONT.GT.MXCONT ) THEN
         NCONT = MXCONT
         WRITE (LUNERR,2000) MXCONT
      END IF

      IF ( NLEVS.LE.0 ) THEN

C  *     Find contour levels for user specified increment, RLEVS(1),
C        when NLEVS=0 or for user specified number of levels. The
C        number of levels is returned in NCONT.

         CINC = RLEVS(1)
        
C  *     CONSCL finds the actual levels.

         CALL CONSCL (NLEVS, FMIN, FMAX, NCONT, CINC)

         NCONT = MIN(NCONT,MXCONT)       
         NLEVS = NCONT
         DO 40 I=1,NCONT
            RLEVS(I) = FMIN + (I-1)*CINC         
 40      CONTINUE

      END IF

C  *  Now that the number of contour levels is fixed, we can set colors
C     over the "reverse" range 7.-2., from blue to magenta. Assume a
C     maximum of MXCOLR colors.

      IF ( LINCOL(1).EQ.99. )
     >   CALL CLRNGE ( MXCOLR, RWORK(ICOLRY), NCONT, LINCOL )

C  *  Process one contour level at a time.

      DO 200 LEV = 1, NCONT

C  *     Find all contour lines at this level.

         ACONT = RLEVS(LEV)

         CALL CONT2D ( IDIM, JMAX, FUNCT, XGRID, YGRID, IMIN, IMAX,
     +                 JMIN, JMAX, NDIM,  ACONT, MAXPTS, RWORK(IXCON),
     +                 RWORK(IYCON), MXLSEG, IWORK, IADIM,
     +                 RWORK(IIA) )

C  *     Process contour line(s) on this level.

         NLINES = IWORK(1)
         LPRLEV(LEV) = NLINES

         IF ( NLINES.GT.0 ) THEN
            IWORK(1) = 1
            IF ( LINTYP(1).NE.99 ) THEN
                  LTYP = LINTYP(LEV)
                  THK  = LINTHK(LEV)
               ELSE
                  LTYP = 3
                  THK  = 1
            END IF

            COLOR = LINCOL(LEV)
            IF (LINCOL (LEV) .EQ. 8) COLOR = 0

C  *        IWORK(I) points to the start of the Ith line. (Note
C           that IWORK(NLINES+1) is defined in CONT2D even though
C           the line itself does not exist.)  Draw contour line(s)
C           with labels if requested (no symbols).

            DO 150 LINE = 1, NLINES
               I1   = IWORK(LINE)
               NPTS = IWORK( LINE + 1 ) - I1
               IX   =   IXCON - 1 + I1
               IY   =   IYCON - 1 + I1
               CALL POLYLINE (NPTS, RWORK(IX), RWORK(IY),
     >                        LTYP, THK, -1, COLOR, IER)

               IF ( ILEGND.EQ.0 ) THEN
                  XLBL = RWORK( IX + NPTS/2 )
                  YLBL = RWORK( IY + NPTS/2 )
                  IF ( XLBL.GE.XMIN .AND. XLBL.LE.XMAX .AND.
     >               YLBL.GE.YMIN .AND. YLBL.LE.YMAX )
     >               CALL RLINT ( LEV, XLBL, YLBL )
               END IF
  150       CONTINUE

         END IF
  200 CONTINUE


      RETURN
 2000 FORMAT('  QCONT: WARNING!  Number of contour levels limited to ',
     >    I4,' as set by MXNLEV')
      END
