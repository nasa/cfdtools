C+**********************************************************************
C
      SUBROUTINE PLLEVS ( NDIM, JDIM, FUNCT, XGRID, YGRID, JMIN, JMAX,
     +                    KMIN, KMAX, NLEVS, RLEVS, LINTYP, LINCOL,
     +                    ILEGND, XLEGND, YLEGND, TXLGND, HTLGND,
     +                    NDGLBL, ANGLBL, ISIZ1, ISIZ2, RWORK, IWORK )
C
C ACRONYM: PLot contour LEVelS
C          --           ---  -
C PURPOSE:
C    PLLEVS is a general-purpose routine for identifying and plotting
C    all requested contour levels within a given region. The function
C    values within the region must be represented as a (pseudo) rect-
C    angular array.  The corresponding coordinate arrays may be 1- or
C    2-dimensional,  according to whether the region is truly or only
C    pseudo-rectangular.
C
C    Provision is also made for an optional legend,  direct  labeling 
C    of contours  or  no labeling.  (The legend type  is inextricably
C    tied in with actual plotting of  the contours because of the way
C    the midpoints of the curves are labeled.)  Use of color in draw-
C    ing contour lines may obviate the need for labels.
C
C    PLLEVS does not handle plotting of the boundary of the contoured
C    region, or of axes, etc.   See modules PLCURV, PLGEOM, etc., for
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
C             Plot level number or actual level at mid-pt. of line.
C    IF lines labeled with level numbers THEN
C       Plot legend at indicated location.
C
C NOTES:
C    There is  a   certain ambiguity about the meaning of ISIZ1 and
C    ISIZ2 as used here.   They are being retained in order to make
C    the 1986 update of PLLEVS transparent  to  current users.  The
C    program will issue  warning messages  and  return partial plot
C    data if work space is inadequate. This should not be a problem
C    since the  RWORK  requirements have dropped  from  6*ISIZ1  to
C    4*ISIZ1.   However, PLLEVS now requires ISIZ2 words of storage
C    for the arrays RLEVS, LINTYP and  LINCOL because the number of
C    contour levels will not always be known in advance.
C         
C    A good deal of confusion has arisen over whether the 2D arrays
C    have to be stored one way or the other.  The answer is that it
C    doesn't matter!   As long as XGRID(J,K) and YGRID(J,K) are the
C    coordinates to be shown horizontal  and  vertical respectively
C    on a (normal) plot,  and  they  do  correspond  to  the  value
C    FUNCT(J,K), and not FUNCT(K,J), say,  then either way of stor-
C    ing these arrays is correct, and the user should choose which-
C    ever way is more convenient.
C
C ARGUMENTS:
C     ARG   TYPE I/O/S  DIM     DESCRIPTION
C    NDIM     I    I     -      No. of dimensions of XGRID and YGRID.
C    JDIM     I    I     -      Effective row dimension of FUNCT(*,*)
C                               in the routine setting up this array.
C    FUNCT    R    I JDIM,KMAX  Function values to be contoured are in
C                               FUNCT(JMIN:JMAX,KMIN:KMAX).
C    XGRID,   R    I  JDIM,*    The coords. matching FUNCT(J,K) are:
C    YGRID                      (  XGRID(J), YGRID(K) ) if NDIM=1;
C                               (XGRID(J,K),YGRID(J,K)) if NDIM=2.
C                               (* = 1 or KMAX.)
C    JMIN,JMAX,          -      Define a window within the FUNCT array
C    KMIN,KMAX                  for the portion to be contoured.
C    NLEVS    I    I     -      NLEVS>0 means plot NLEVS user-specified
C                               contour levels.
C                               NLEVS=0 means plot contour levels from an
C                               array with increment RLEVS(1) and maximum
C                               and minimum from requested data window.
C                               NLEVS<0 means means plot approximately
C                               ABS(NLEVS) nicely scaled contour levels.
C    RLEVS    R   I/S  ISIZ2    Array of contour levels. Calling program
C                               must provide at least ISIZ2 words even if
C                               the levels are not user-specified. RLEVS
C                               is not modified on return.
C                               RLEVS(1) = desired contour level increment
C                               if NLEVS=0.
C    LINTYP   I    I   ISIZ2    LINTYP(I) defines the type and thick-
C                               ness of the Ith contour line plotted.
C                               The first digit of each  defines  the
C                               thickness of the line drawn in multi-
C                               ples of .01 inch.  1:9 is valid.  The
C                               second digit determines the line type:
C                                   1   ...  Solid
C                                   2   ...  Dashed
C                                   3   ...  Dotted
C                                   4   ...  Chain-dashed
C                                   5   ...  Chain-dotted
C                               LINTYP(1)=99 means plot all contour
C                               lines solid with default thickness.
C    LINCOL   I    I   ISIZ2    LINCOL(I) specifies the color of the
C                               Ith contour line  plotted,  but if
C                               LINCOL(1)=99,  PLLEVS will set a default
C                               color range from blue to magenta with
C                               increasing level.
C                               Any other number outside the range 1-8
C                               will be ignored by PLCURV, which will
C                               use the output device default color.
C                                   1   ...  Black
C                                   2   ...  Magenta
C                                   3   ...  Red
C                                   4   ...  Yellow
C                                   5   ...  Green
C                                   6   ...  Cyan
C                                   7   ...  Blue
C                                   8   ...  White
C    ILEGND   I    I     -      ILEGND=0 means contour lines are lab-
C                               eled directly (no legend).
C                               ILEGND>0 means lines are labeled with
C                               level numbers and a legend is provided.
C                               ILEGND<0 means no legend or labels will
C                               be drawn.
C    XLEGND   R    I     -      Horizontal position of legend.
C    YLEGND   R    I     -      Vertical position of legend.
C    TXLGND  C*10  I     -      Legend text ( 10 packed characters ),
C                               heading column of actual contour vals.
C    HTLGND   R    I     -      Height of legend, and also of the lab-
C                               els on the contour lines.
C    NDGLBL   I    I     -      No. of digits after the decimal point
C                               to show on contour labels.
C    ANGLBL   R    I     -      Angle ( degrees ) for contour labels.
C    ISIZ1    I    I     -      Used for dimensions of the floating-
C                               point work-space arrays needed by PLLEVS
C                               and CONT2D. ( Try ISIZ1 = 1000. )
C    ISIZ2    I    I     -      Dimension of integer work-space arrays
C                               needed by CONT2D involving expected max-
C                               imum number of separate lines for  any
C                               given contour level. (Try ISIZ2 = 30.)
C                               It is also used by PLLEVS as the dimen-
C                               sion of RLEVS, LINCOL, and LINTYP  and
C                               maximum number of contour levels plotted.
C    RWORK    R    S  ISIZ1*4   RWORK(1:2*ISIZ1) holds the abscissas and
C                               ordinates of contour lines returned by
C                               CONT2D;
C                               RWORK(2*ISIZ1+1:3*ISIZ1) is scratch for
C                               CONT2D;
C                               RWORK(3*ISIZ1+1:4*ISIZ1) is scratch for
C                               PLLEVS.
C    IWORK    I    S  ISIZ2*2   IWORK(1) gives the number of contour
C                               lines found in one call to CONT2D;
C                               IWORK(I) for I=2,ISIZ2 points to the
C                               start of the Ith contour line returned
C                               in RWORK ( reused for each call to 
C                               CONT2D ); the remainder of IWORK is
C                               used by PLLEVS to keep track of the
C                               number of contour lines per contour
C                               level as required for PLCLEG.
C
C
C EXTERNAL REFERENCES:
C    NAME      DESCRIPTION
C    CONT2D    Prepares data for contour plotting.
C    CONSCL    Finds "nice" scaling of a number of values over a given
C              data range.
C    PLCURV    Plots curve of indicated type/thickness/color.
C    PLCLEG    Plots a legend for the contour lines.
C    Also:     DISSPLA routines for line labels, legend and color.
C
C KNOWN SYSTEM DEPENDENCIES:
C    > DISSPLA graphics software.
C    > Color output depends on plotting hardware.
C    > Requests for color are ignored when color is not available.
C    > Interpretation of color selection is device dependent.
C
C STANDARDS VIOLATIONS: None.
C
C ENVIRONMENT:  VAX/VMS, CRAY-1S -- FORTRAN 77
C
C DEVELOPMENT HISTORY:
C     DATE   INITIALS   DESCRIPTION 
C    Jul. 86    RCL     Changed the contour line finder from PLCONTOUR
C                       to CONT2D (algorithm from PLOT3D).
C                       Added optional automatic selection of "nicely"
C                       scaled contour levels; selection of contour
C                       levels from input contour increment; and default
C                       line colors and type.
C    Feb. 85    RCL     Added color capability and option for no labels.
C    Nov. 82    DAS     Mnemonics and work-space allocation revised;
C                       now one of the "building blocks" for composite
C                       plots.  See also PLCURV, PLFREE, etc.
C    Sep. 82    RGL     Adapted from part of original KONPLT and trans-
C                       lated to FORTRAN 77.
C    1980-81  PJT/CJG   Original coding for PLTSRF application.
C
C AUTHORS: Jeff Trosin, Carol Green, Ron Langhi, David Saunders,
C          Rosalie Lefkowitz, Informatics General Corp., Palo Alto, CA.
C
C-**********************************************************************

      CHARACTER  TXLGND*10

      DIMENSION  FUNCT(JDIM,KMAX), XGRID(JDIM,1), YGRID(JDIM,1),
     +           RLEVS(1), LINTYP(1), LINCOL(1),
     +           RWORK(ISIZ1*4), IWORK(ISIZ2*2)


      CALL HEIGHT ( HTLGND )
      CALL ANGLE  ( ANGLBL )

      IXCON  = 1
      IYCON  = ISIZ1 + 1
      IIA    = ISIZ1 + IYCON
      IRLEVS = ISIZ1 + IIA
      IADIM  = ISIZ1
      INPRLV = ISIZ2 + 1
      MXCONT = MIN(ISIZ2,ISIZ1)
      NCONT  = ABS(NLEVS)
      IF ( NCONT.GT.MXCONT ) THEN
         NCONT = MXCONT
         WRITE (6,2000) MXCONT
      END IF

C  *  Save original RLEVS.

      DO 20 I=1,MXCONT
         RWORK(IRLEVS+I-1) = RLEVS(I)
 20   CONTINUE

      IF ( NLEVS.LE.0 ) THEN

C  *     Find contour levels for user specified increment, RLEVS(1),
C        when NLEVS=0 or for user specified number of levels. The
C        number of levels is returned in NCONT.

         CINC = RLEVS(1)
        
C        Find the range (minimum and maximum) of the function
C        in the region JMIN to JMAX, KMIN to KMAX.

         FMIN = FUNCT(JMIN,KMIN)
         FMAX = FMIN
         DO 35 K=KMIN,KMAX
            DO 30 J=JMIN,JMAX
               FMIN = AMIN1(FMIN,FUNCT(J,K))
               FMAX = AMAX1(FMAX,FUNCT(J,K))
 30         CONTINUE
 35      CONTINUE

C  *     CONSCL finds the actual levels.

         CALL CONSCL(NLEVS,FMIN,FMAX,NCONT,CINC)

         NCONT = MIN(NCONT,MXCONT)       
         DO 40 I=1,NCONT
            RLEVS(I) = FMIN + (I-1)*CINC         
 40      CONTINUE

      END IF

C  *  Now that the number of contour levels is fixed, we can set default
C     colors over the "reverse" range 7-2, from blue to magenta.

      IF ( LINCOL(1).EQ.99 ) THEN
         CINC = 6./REAL(NCONT)
         DO 60 LEV=1,NCONT
            LINCOL(NCONT-LEV+1) = 2 + REAL(LEV-1)*CINC
 60      CONTINUE
      END IF

C  *  Process one contour level at a time.

      DO 200 LEV = 1, NCONT

C  *     Find all contour lines at this level.

         ACONT       = RLEVS(LEV)

         CALL CONT2D ( JDIM, KMAX, FUNCT, XGRID, YGRID, JMIN, JMAX,
     +                 KMIN, KMAX, NDIM,  ACONT, ISIZ1, RWORK(IXCON),
     +                 RWORK(IYCON), ISIZ2, IWORK, IADIM,
     +                 RWORK(IIA) )

C  *     Process contour line(s) on this level.

         NLINES = IWORK(1)
         IWORK(INPRLV-1+LEV) = NLINES

         IF ( NLINES.GT.0 ) THEN
            IWORK(1) = 1
            IF ( LINTYP(1).NE.99 ) THEN
               LTYP = LINTYP(LEV)
            ELSE
               LTYP = 11
            END IF
            ITHK = LTYP/10
            LTYP = LTYP - 10*ITHK
            ICOLOR = LINCOL(LEV)

C  *        IWORK(I) points to the start of the Ith line. (Note
C           that IWORK(NLINES+1) is defined in CONT2D even though
C           the line itself does not exist.)  Draw contour line(s)
C           with labels if requested (no symbols).

            DO 150 LINE = 1, NLINES
               I1   = IWORK(LINE)
               NPTS = IWORK( LINE + 1 ) - I1
               IX =   IXCON - 1 + I1
               IY =   IYCON - 1 + I1

               CALL PLCURV ( NPTS, RWORK(IX), RWORK(IY),
     +                       LTYP, ITHK*.01, 0, 0, ICOLOR )

               XLBL = RWORK( IX + NPTS/2 )
               YLBL = RWORK( IY + NPTS/2 )

               IF ( ILEGND.GT.0 ) THEN
                  CALL RLINT  ( LEV, XLBL, YLBL )
               ELSE IF ( ILEGND.EQ.0 ) THEN
                  CALL RLREAL ( ACONT, NDGLBL, XLBL, YLBL )
               END IF
  150       CONTINUE

         END IF
  200 CONTINUE

C  *  Print legend if desired.

      IF ( ILEGND.GT.0 )
     +   CALL PLCLEG ( NCONT, RLEVS, IWORK(INPRLV), NDGLBL, 
     +                 XLEGND, YLEGND, TXLGND, HTLGND )

C  *  Return original RLEVS.

      DO 240 I=1,MXCONT
          RLEVS(I) = RWORK(IRLEVS+I-1)
 240  CONTINUE


      RETURN
 2000 FORMAT('  PLLEVS: WARNING!  Number of contour levels limited to ',
     >   I4,' as set by ISIZ2')
      END
