C+------------------------------------------------------------------------------
C
      PROGRAM QPLOT
C
C
C     Description and usage:
C
C           This program is intended for "generic" DISSPLA plots with
C        minimal user involvement.  The number of curves per plot is
C        (almost) unlimited, and any number of plot frames may be produced
C        in one run if the appropriate separator is used between the data
C        sets.  The plot parameters are fully defaulted - a minimal input
C        file consists of one or two columns of numbers and will produce a
C        reasonable-looking plot of a single curve.  Multiple curves or
C        frames, titles and legends, scaling, grids, colored lines, data
C        transformations, etc. are all possible but strictly optional.
C        Blank lines anywhere in the input file are ignored.
C
C           The input file consists of the following:
C
C        (1)  A plot title, subtitle, labels for the X and Y axes, and
C             virtually any number of caption lines, each on a separate
C             input line.  Three text lines imply no subtitle.
C             Text lines after the 4th are treated as captions.
C             Blank text strings may be specified in quotes or by "" or ''.
C
C        (2)  An optional namelist denoted OPTIONS, which may be used
C             to override the default plot size, axis scaling, color,
C             and line type; add legend entries, create a background grid,
C             or specify the columns from which the data is to be read.
C             Comments may also be included, but will not appear in the
C             plot.  The namelist must precede the data for each curve.
C
C        (3)  The data to be plotted, in multi-column format.  If more
C             than one curve is to be plotted, each set of data points
C             must be followed by an OPTIONS namelist or a line beginning
C             with "END CURVE."  The columns of data may be separated by
C             blanks, tabs, commas, colons, or equal signs (see SCANNR
C             header for details).
C
C        (4)  The final curve of a frame must be terminated by a line
C             beginning with "END FRAME" if more data is to follow.
C
C           This version permits (nested) indirect data files at the CURVE
C        level - see "Notes" below.
C
C           This version is set up for a maximum of 5,000 curves.  The number
C        of points per curve is arbitrary, but the total number of points
C        per frame must not exceed 500,000.  There is no limit on the number
C        of frames.  The legend may have up to 30 entries, but 20 or fewer
C        are better suited to the default height and width of the plot area.
C
C           This version provides a choice of plot devices (screen and/or
C        metafile).  The file qplot.config in a standard location defines
C        the devices supported at a given site.  Plots may be previewed
C        frame-by-frame with optional metafile output.  A "plot-the-rest"
C        option (all remaining frames to the metafile) is also provided.
C
C           The plot device(s) may be specified on the command line, as may
C        the input data file name.  This version also permits specification
C        of an output metafile name from the command line.  Examples:
C
C           qplot /tek /ps xy.dat     ! A qplot.ps PostScript file is output
C           qplot -tek -ps xy.dat xy  ! An   xy.ps   "    "    "    "    "
C           qplot -ps qwerty.plt same ! A qwerty.ps  "    "    "    "    "
C
C        (Use of "same" in the third example saves retyping awkward names.)
C
C           At NASA-Ames, a number of companion utilities are available:
C
C        (1)  Multiple QPLOT input files may be merged using QMERGE.
C        (2)  Programs SMOOTH, for data fitting, and PROFILE, for airfoil
C             geometry manipulation, produce output files in QPLOT format.
C
C
C     Input data:
C
C        Any input line (including any title line) may contain a trailing 
C        comment starting with '!'.  Blank or comment-only lines are ignored.
C
C        Name     Type   Description
C        TITLES    C     Plot title, axis labels and captions.  Max. 132 
C                        characters each. The number of characters which will
C                        fit on a line is less than this, but the higher limit
C                        permits inclusion of special instructions if needed.
C                        Interpretation of the input depends on the number of
C                        lines of text found, as follows:
C
C                           0 lines  =>  no titles, labels or captions
C                           1 line   =>  single title
C                           2 lines  =>  title and subtitle
C                           3 lines  =>  title and X, Y axis labels
C                           4 lines  =>  title, subtitle, and X, Y axis labels
C                           5 lines  =>  title, subtitle, X, Y axis labels and
C                               .        one caption line
C                               .
C                               .
C                           n lines  =>  title, subtitle, X, Y axis labels and
C                                        n - 4 caption lines (n = 30 currently)
C
C                        If IDENT = .TRUE., two of the caption lines are used
C                        for a blank line and a line with file name and date.
C
C                        A significant blank title line may be indicated with
C                        a quoted blank or '' or "".  Similarly, numeric looking
C                        text should be distinguished from the data proper using
C                        quotes.  Embedded square brackets '[' and ']' have 
C                        special significance for alternative character sets
C                        or other text manipulations, although these escape
C                        characters may be changed.  See Notes.
C
C        X, Y      R     Plot data arranged in vertical columns.  Default
C                        correspondences are X = <col. 1>, Y = <col. 2>.
C
C
C     NAMELIST input data:
C
C        /OPTIONS/       (NOTE:  This NAMELIST may be freely omitted!)
C
C        Name      Type  Description
C
C        ANGLE      R    Angle in degrees to rotate current curve through
C                        about (CENTERX, CENTERY).  Positive is counter-
C                        clockwise.  See also RADIANS.  A rotation is
C                        applied BEFORE any shifting/scaling.
C                        Synonym:  DEGREES   Default: 0 degrees
C
C        CENTERX,   R    Center of rotation for rotating current curve.
C        CENTERY         Synonyms:  XCENTER, YCENTER.   Default:  (0, 0)
C
C        CLASS      C    A string containing a special label for the plot 
C                        which appears at the top and bottom of the page
C                        in large print.  It was inspired by the need to 
C                        label classified documents.
C
C        COLOR      C    Choice of color for individual curves, from among
C                        'WHITE' (default), 'MAGENTA', 'RED', 'YELLOW',
C                        'GREEN', 'CYAN', 'BLUE', and 'BLACK' (usually not
C                        appropriate).  The various colors may (but generally
C                        should not be) mixed with the different symbol and
C                        line types.
C
C        COMMENT    C    A comment string may be used for file documentation.
C                        The string does not appear on the plot.
C
C        DEBUG      L    A logical flag which controls printing of diagnostic
C                        messages from DISSPLA (as well as job information
C                        along the left side of the plot).  Normally .FALSE.
C
C        ESCAPE    C*2   A pair of characters used to turn on and off DISSPLA's
C                        "instruction" alphabet for math. symbols, subscripts,
C                        etc.  The default is '[]' (reset for each frame).
C
C        EXPLAN     C    An optional subtitle for CLASS.
C                        
C        FIT        C    Plotted points may be connected in a variety of
C                        ways, specified by input FIT for each curve.  The
C                        default is 'LINEAR'; the other choices are
C                        'TIGHT' and 'LOOSE' parametric piecewise cubics
C                        using locally-developed routine PLSFIT.  For the
C                        special case of a smoothly-closed curve, e.g. for
C                        geometrical figures, use FIT = 'CLOSED'.
C
C        FORMAT     I    Plot format:
C                        2   Box with tick marks on all four sides.
C                        1   Ticks extending from the lower left corner
C                            of the box along two axes (default).
C                        0   No box and no axes.
C                       -1   No box; axes extending from lower left.  
C                       -2   No box; axes crossing at or near the origins.
C                        (Enabling the secondary Y-axis forces FORMAT=2.)
C
C        GRID       I    If GRID > 0, number of grid lines per scale division;
C                        if < 0, number of tick marks per scale division. 
C                        Default is GRID = -1.  (Note:  If GRID is specified
C                        as 0, it is assigned the default value.)
C
C        HEIGHT,    R    Nominal plot dimensions, in inches.  Default is
C          WIDTH         5" x 5" (portrait) or 3.75" x 7.5" (landscape).
C                        This is small enough that a viewgraph of a typical
C                        plot with titles, legends, and a few captions
C                        will fit a projector comfortably.
C
C        IDENT      L    .TRUE. causes QPLOT to insert and extra caption
C                        line (separated by a blank line from data captions)
C                        showing  Data file: xxx.xxx   Plot date: xx-xxx-xx.
C                        .FALSE. is the default.
C
C        LEGEND     C    Optional legend string for each curve.  Max. 132
C                        characters.  At most 30 legends will appear in the 
C                        plot.  A "blank" legend such as ' ' will appear 
C                        with a blank text field.
C
C        LINE       C    Optional flags to indicate the line type to be
C                        used.  The default is a solid line with symbols
C                        for each curve.  Other options are:  'SYMBOLS',
C                        'SOLID', 'DOTS', 'DASHES', LONGDASHES', 'CHAINDOTS',
C                        'CHAINDASHES', and 'THICK'.  (Default 'CONNECTED'
C                        may, but need not, be specified.)
C
C        NOTE,      C    Annotation text, placed according to XNOTE, YNOTE
C          NOTEn         (for the common case of a single annotation) or
C                        XNOTEn, YNOTEn where n = 1, 2, 3, ..., 9 (up to
C                        9 more annotations per X/YSHIFT/SCALE).  E.g.:
C                           NOTE = 'Freestream velocity' or
C                           NOTE1 = 'Point of inflection'
C                        The maximum total number of annotations per frame is
C                        30.  All strings may be up to 132 characters long.
C
C        NOTEBOX    I    Controls the boxes around any annotation strings
C                        entered via NOTE[n].  If boxes are drawn, they
C                        allow a border of .05" around the text, which is
C                        .12" high.  It may be desirable to prevent other
C                        plot elements from overwriting the text area by
C                        "blanking" it.
C                           0 means no frames and no blanking;
C                           1 means no frames but boxes are blanked [default];
C                           2 means the boxes are framed and blanked;
C                          >2 thickens the frames: a factor of NOTEBOX - 1
C                             is applied to the pen width.
C
C        ORIENT     C    Plot orientation.  For 8.5" x 11" paper,
C                          'PORTRAIT' = X-axis along 8.5" edge (vertical);
C                          'LANDSCAPE' = X-axis along 11" edge (horizontal).
C                        (Default is 'PORTRAIT')
C
C        PLOT       C    Plot type:
C                          'LINEAR' = linear scales on both axes (default);
C                          'SCALE'  = equal scaling on linear axes;
C                          'LOGX'   = log X-, linear Y-axis;
C                          'LOGY'   = log Y-, linear X-axis;
C                          'LOGLOG' = log scales on both axes;
C                          'POLAR'  = THETA vs. RADIUS polar diagram.
C
C        RADIANS    R    Angle of rotation for a curve, in radians.  See
C                        ANGLE above.
C
C        SCALEX,         Linear transformation coefficients for the abscissas:
C          SHIFTX        X (plotted) = SCALEX * X (input) + SHIFTX.
C                        Scaling/shifting is applied AFTER any rotation.
C                        Synonyms:  XSCALE, XSHIFT   Defaults: 1. and 0.
C
C        SCALEY,         Linear transformation coefficients for the ordinates,
C          SHIFTY        as for SCALEX, SCALEY above.
C
C        SPLINE     C    Spline type (parametric or nonparametric):
C                          'PARAMETRIC' = X vs. arc & Y vs. arc (default);
C                          'STANDARD'   = Y vs. X (X increasing or decreasing);
C                          'REVERSED'   = X vs. Y (Y increasing or decreasing);
C                        used in conjunction with FIT='LOOSE' or 'TIGHT' and
C                        PLOT='LINEAR' or 'SCALE' (not log. or polar plots;
C                        parametric fits are performed in normalized space;
C                        nonparametric fits are performed in data space.
C
C        STYLE      C    Overrides default text style for an entire frame,
C                        primarily for "final" publication or presentation
C                        charts using shaded characters (SLOW-BULKY-OBSCURE!)
C                        Only a selection of DISSPLA's options are listed
C                        here.  (The rest are available under their normal
C                        names - see the DISSPLA manual, if you must...)
C                          'DEFAULT'    = double-stroke serif;
C                          'LIGHT'      = single-stroke sans serif;
C                          'STICK'      = rudimentary, fast, single-stroke;
C                          'TIMESROMAN' = like newspaper text (shaded);
C                          'FUTURA'     = very plain sans serif (shaded);
C                          'HELVETICA'  = sans serif (shaded).
C                   
C        SYMBOL     I    Curve marking symbol number, an integer between
C                        0 and 18 (higher values wrap, -ve means use default
C                        as determined by line type, thus for special line
C                        types -1 means NO symbol, and for connected points
C                        and symbols alone -1 means cycle through the list).
C                        See QPLOT Cheat Sheat or the DISSPLA manual for
C                        the key.
C
C        XARROW,    R    (X, Y) coordinates (in data units) of the tip of an
C          YARROW,       arrow drawn from the CENTER of the text block defined
C          XARROWn,      by the corresponding X/Y/NOTE[n].  (The arrow is
C          YARROWn       clipped by a .05" border around the text.)  Scaling
C                        and shifting applies as for X/YNOTE[n], q.v.
C                        N.B. Set YARROW[n] = YNOTE[n] if a precisely
C                        HORIZONTAL arrow is desired (whether NOTE[n] is blank
C                        or not).  Likewise, set XARROW[n] = XNOTE[n] if a
C                        precisely VERTICAL arrow is desired.  It will be
C                        realigned with the center of the text box.
C
C        XCOL,      I    Column numbers from which X and Y data are to be
C          YCOL          read.  Defaults are XCOL = 1, YCOL = 2.  If XCOL or
C                        YCOL are zero, the ordinals i = 1,2,3, ... will be
C                        used.
C
C        XMAX,      R    Optional X axis scaling limits.  The input values
C          XMIN          will be used without rounding further. By default,
C                        "nice" axis limits will be calculated.
C
C        XSTEP      R    X-axis step size.  Default is auto-scaling, with
C                        increasing values.  Use -999.0 for decreasing.
C
C        XNOTE,     R    (X, Y) coordinates (in data units) of the bottom
C          YNOTE,        left corner of the text entered as NOTE or NOTEn, q.v.
C          XNOTEn,       Any X/YNOTE[n] entered with an XSCALE or YSCALE or
C          YNOTEn        XSHIFT or YSHIFT will be adjusted accordingly.
C                        See also X/YARROW[n].
C
C        XNUMBERS,  C    X- and Y-axis numbering type.  Default is 'REAL';
C          YNUMBERS      use 'INTEGER' or 'WHOLE' to suppress decimal points.
C
C        YMAX,      R    Optional Y axis scaling limits.  Default is automatic
C          YMIN          scaling, as for the X axis.  Choose YMIN > YMAX for
C                        an inverted axis.
C
C        YSTEP      R    Y-axis step size.  Default is auto-scaling, with
C                        increasing values.  Use -999.0 for decreasing.
C
C        Y2LABEL    C    Label for secondary Y-axis if any (right side of
C                        plot box).
C
C        Y2SCALE,   R    Scale and shift factors applied to the tick marks
C          Y2SHIFT       of the secondary Y-axis.  Y2SCALE = Y2SHIFT = 0.
C                        suppresses numbering of these ticks, and this is
C                        the default.
C
C     Input format:
C
C        Sample Plot Title                        <-- The first eight lines of
C        Optional Subtitle                            text are optional, but
C        X (horizontal) Axis Label                    must be in order.
C        Y (vertical) Axis Label
C        Caption 1
C        Caption 2
C        Caption 3
C        .....
C
C         $OPTIONS                                <-- The NAMELIST is optional.
C         XSTEP = 0.5,                                It must begin in col. 2.
C         YMIN = 0.0,
C         YMAX = 10.0,
C         LINE = 'SOLID',                         <-- Note that CHARACTER type
C         $END                                        inputs via NAMELIST must
C                                                     be enclosed in quotes.
C        6,3
C        4,5                                      <-- Columns of data may be
C        8,9                                          separated by spaces,
C        5,8                                          tabs, commas, colons, or
C        2,4                                          equal signs.
C        1,3
C        END CURVE                                <-- "End-of-curve" mark (not
C                                                     required if next curve
C         $OPTIONS                                    has NAMELIST).
C         XCOL = 3,
C         LEGEND = 'Curve #2',
C         LINE = 'LONGD',                         <-- String flags may be
C         $END                                        abbreviated.
C        .546   .53472    .234
C        .123   .76598    .345                    <-- Multiple column input.
C
C           (etc.)
C
C        @file-name                               <-- Optional indirect file
C                                                     containing complete
C                                                     curves.  A namelist
C                                                     immediately preceding
C                                                     @file-name (or EOF)
C                                                     applies to the next
C                                                     curve read.
C           (etc.)
C
C        .23646  5.465457  1.46474
C        END CURVE                                <-- If "END FRAME" is used
C                                                     here instead, the entire
C           (etc.)                                    data set input format
C                                                     may be repeated.
C
C
C     External files:
C
C        Unit  I/O/S  Description
C         2      O    Diagnostics file.
C         3      O    Plot metafile.
C         5    I      Keyboard inputs.
C         6      O    Screen prompts and plot previews.
C         7+   I      Input data file(s).
C
C        
C     External references:
C
C        ALPHA    Distinguishes alphabetic from purely numeric data.
C        DONEPL   Closes plot file.  (DISSPLA)
C        GETPRM   Gets command line parameter string.
C        EXIT     Stops execution.  See Notes below.  (VMS)
C        LAYOUT   Sets default plot parameters, error-checks user inputs.
C        LOOKUP   Dictionary lookup, with abbreviations and synonyms.
C        NLCHEK   Identifies legal or suspicious NAMELISTs.
C        PARSENAM Used to derive an output metafile name from the input file.
C        PLOTDEV  Determines available plot devices; selects; initializes.
C        QUICK    Plot driver which calls DISSPLA routines.
C        RESCALE  Linear transformation of data packed in an array.
C        READCOLS Reads a "curve" to be plotted, in multi-column format.
C        ROTATE   Rotates one packed curve about an arbitrary point.
C        SCAN2    Identifies significant characters.
C        UPCASE   Converts a string to uppercase.
C        WRAPPER  Checks/corrects match of first and last points of curve.
C
C
C     Environment:  FORTRAN 77
C                   DEC/OpenVMS; SGI/IRIX; CA-DISSPLA Version 11.5-9312
C
C
C     Notes:
C
C        (1)  IMPLICIT NONE, 8-character symbolic names, and use of "!"
C             for comments are non-standard (until Fortran 90).
C
C        (2)  NAMELIST input is non-standard.  The VMS NAMELIST symbol is
C             $ or &.  All NAMELIST input must begin in column 2.  Spaces,
C             tabs, or commas are valid separators.  The data may appear
C             in any order.  It is not necessary to assign values to all
C             of the NAMELIST input variables.
C
C        (3)  CALL EXIT (N) is VMS-specific.  A parameter value of 1 indicates
C             successful completion, while 3 is defined here to mean abnormal
C             termination.  These values are available to the system after
C             the end of the run.  <Disabled in multi-plot-device version.>
C
C        (4)  The OPEN statement for the input data file contains the
C             VMS-specific parameter READONLY, which permits input across
C             different accounts.  <IRIX too - see the OPENER utility.>
C
C        (5)  A large variety of type sizes and styles are available.  See
C             Chapter 6, Section 11 of the DISSPLA manual for details on the
C             "instruction string" commands.  The escape to the special
C             command set is defaulted to the character '[' and the return is
C             signaled by ']'.  These characters can therefore only appear
C             in a plotted text string if ESCAPE is used to switch the escape
C             characters.
C
C        (6)  Drastic re-sizing (using HEIGHT and WIDTH) may require changes
C             in character sizes or spacing.  This may be done by embedding
C             the appropriate directives in the input strings, or more
C             generally, by modifying a copy of plot routine QUICK and re-
C             linking.  See command procedure QPLOTLNK.COM for an example.
C
C        (7)  Unfortunately, any use of embedded directives in a string
C             disables the PostScript output option's hardware character
C             capability for that string: DISSPLA reverts to stroke text.
C             To avoid a mixture of text strings in PostScript output,
C             specify "pss" rather than "ps" and all text will be stroked.
C
C        (8)  Any of the CHARACTER type flags (LINE, ORIENT, and PLOT) may
C             be abbreviated to their shortest unambiguous equivalent.  There
C             are synonyms available for many of them.
C
C        (9)  For POLAR diagrams, angle THETA must be in column 1, or given
C             by THCOL in $OPTIONS (or just use XCOL).  The angular unit is
C             radians.  RADIUS is in column 2, or as specified by RCOL (same
C             as YCOL).  A non-zero radial origin may be input with RMIN
C             (= YMIN).  The other scaling inputs are ignored, as are the
C             X- and Y-axis labels.
C
C       (10)  This version supports indirection: any data line after the titles
C             containing @file-name means QPLOT opens and reads the indicated
C             file, which should NOT START with more titles and should not end
C             with a partial curve, but is otherwise unrestricted.
C
C             There are two likely uses of indirect files:
C
C             o  for overlaying a given self-contained curve repeatedly on other
C                curves which are varying with the plot.
C
C             o  for plotting more than one pair of columns from the same file
C                without having to replicate the file.
C
C             Both usages (particularly the second) demand the following:  a
C             namelist followed IMMEDIATELY by @file-name (or by end-of-file or
C             by another namelist) applies to the next curve read.  (Normally,
C             such read exceptions cause the current curve and its controls to
C             be processed, followed by reinitializing of the curve controls
C             prior to dealing with the exception.  But "no-data-yet" means
C             "don't (re)reinitialize.")
C
C
C     Author:  Robert Kennelly
C              Mail Stop 227-6
C              NASA-Ames Research Center
C              Moffett Field, CA  94035
C
C
C     Development history:
C
C        14 May  1982 R.Kennelly Original design and coding.
C        29 Sep. 1982    RAK     Drives multiple-plot version of QUICK.
C        29 Nov. 1982 D.McKernan Allows multiple curves per plotframe.
C        31 Dec. 1982    RAK     Restructured for proper termination.
C        17 Feb. 1983    RAK     Modified to permit axis scaling and line
C                                type inputs.
C        28 Feb. 1983    RAK     Generalized TYPE input.
C        23 Mar. 1983    RAK     Fixed input bug - was skipping too many
C                                lines when OPTIONS was absent.
C         2 June 1983    RAK     Added explicit OPEN statements (except for
C                                unit 9, which is handled by DIP).
C        13 June 1983    RAK     Added provision for multiple frames per run.
C        24 June 1983    RAK     Added a BACKSPACE for when NAMELIST is
C                                absent (bug), corrected defaulting, spruced
C                                up READ protection.
C         3 Nov. 1983    RAK     Updated IOCHEK parameter list (CONVER).
C        13 Feb. 1984    RAK     Extensive revisions.  Added legends, grid,
C                                multiple column input, more line type options,
C                                etc.  More error recovery is done in RDCOLS.
C                                New packed data structure for X, Y arrays.
C        23 Feb. 1984    RAK     Added subtitles, with all titles now optional.
C                                COMMENTs now permitted in the namelist.
C        13 Mar. 1984    RAK     Added namelist inputs HEIGHT and WIDTH.
C                                Eliminated ECHOing of input data.  Legend
C                                strings are packed in a CHARACTER variable.
C        23 Apr. 1984    RAK     Corrected counting of titles when blank found.
C                                Any blank line (anywhere) is now ignored.
C                                A number of QMOD variables are declared and
C                                included in the namelist for compatibility.
C        11 May  1984  R.Langhi  Added capability for a logarithmic vertical
C                                axis in the form of namelist option PLOT,
C                                derived from COVPLT.
C        18 June 1984  RAK/DAS   Minor tidying up.  PLOT is converted to
C                                upper case and checked for legality before
C                                calling QUICK.  Added protection against
C                                invoking the plot routine with no data.
C                                Reordered namelist data packing to allow
C                                curves with no data to be ignored.  RDCOLS
C                                arguments reordered.
C        23 Jan. 1985  L.Collins Added 'SCALE' to test for plot types.
C        11 Apr. 1985    RAK     Added SETUP between QPLOT and QUICK,
C                                including keyword handling via LOOKUP. Use
C                                TOTAL rather than MAXPT for QUICK. Drop
C                                SCALE array in favor of XMIN, XMAX, etc.
C        12 Aug. 1985    RAK     Moved default assignment of HEIGHT and WIDTH
C                                down to SETUP.
C        21 Aug. 1985    RAK     ORIENTation input added.
C        10 Sep. 1985    RAK     Added dictionary lookup for line types.
C        26 Sep. 1985    RAK     New line types: 'LONGDASHES' and 'FAT'.
C                                Line types now use plural form (e.g "dots"),
C                                but LOOKUP still accepts singular.  QMOD
C                                namelist updated.
C         4 Oct. 1985    RAK     Symbol type may be specified, by number
C                                using DISSPLA sequence.  LINDAT array is
C                                now INTEGER, using unsorted dictionary.
C                                Eliminated IOCHEK.  Moved normal termination
C                                to the end.  Only the first MAXLEG curves
C                                may now have a legend entry.  Use EQUIVALENCE
C                                to provide synonyms for several inputs.
C        15 Oct. 1985    RAK     Added description of -FLAG for STEP inputs.
C        25 Aug. 1986 R.Lefkowitz Added captions option.
C        16 Oct. 1986    RAK     Changed default GRID to 0 (box around plot).
C                                Added COLOR option.
C        22 Oct. 1986    RAK     The LAST New Option: STYLE, for setting
C                                new default lettering style.  Improved
C                                recovery from dictionary failures: set
C                                string to 'DEFAULT' and re-CALL, so that
C                                defaulting is handled within dictionary.
C        03 Mar. 1987    RAK     REALLY The Last New Option: a fit option
C                                using PLSFIT.  Added X/YINC synonyms for
C                                X/YSTEP.  Dropped QMOD compatibility.
C                                Removed obsolete LEGS variable. SETUP
C                                renamed LAYOUT, and revised so that the
C                                orientation decoding must be handled here.
C         8 Apr. 1987    RAK     Finished the above mods.  Added WRAPPER
C                                to verify/repair the wrap-around data
C                                format required by the CLOSED option.
C                                Repaired too-many-legends warning bug.
C                                Added DEBUG flag to control DISSPLA msgs.
C        17 Apr. 1987    RAK     LDEBUG passed to QUICK as logical unit
C                                number for messages, and as control.
C                                Permit a special line type to have a
C                                symbol as well.  Negative SYMDAT now means
C                                no symbol to be plotted.
C        21 Apr. 1987    RAK     Added SCALEX,Y and SHIFTX,Y inputs.
C         7 May  1987    RAK     Added some equivalent input names to
C                                $OPTIONS: XCOLUMN, YCOLUMN, YBOTTOM.
C        18 May  1987    RAK     Fiddled with message formats.
C        31 July 1987    RAK     Updated header, added some synonyms,
C                                tightened test for 'DEFAULT' fit.
C         9 Oct. 1987    RAK     Added some SHIFT and SCALE synonyms to
C                                $OPTIONS, handling of flakey namelists
C                                while reading titles (lame omission in
C                                previous versions). NLCHEK now uses
C                                GETLINE, which soaks up blanks and
C                                comments. Grouped CALL EXIT (3)s at end.
C                                Test on FRAME when EOF during titles
C                                repaired (was .GT. 0). NLGRIPE writes
C                                namelist tutorial to LWRITE.
C         5 Apr. 1988    RAK     Pass MAXPT to QUICK to permit use of
C                                X, Y arrays as workspace.
C         8 Mar. 1989  M.D.Wong  Added namelist options for printing special 
C                                headings at the top and bottom of the page.
C                                Minor modifications made to spacing of
C                                captions and X axis label.
C        12 May  1989    MDW     Eliminated bug induced by addition of special 
C                                headings.  (MAXTIT should have been left at
C                                8, not raised to 10.)  Also, length of text
C                                strings increased from 81 to 133 to handle
C                                more in-line editing commands.
C        Apr/May 1989    MDW     Translated from DISSPLA to SMDLIB graphics pkg.
C           June 1989    MDW     Ideas from the SMDLIB translation retrofitted
C                                in DISSPLA version: DISSPLA's underlined legend
C                                format replaced by same-line format (using a
C                                new LEGEND - not DISSPLA's utility of that
C                                name); all curve-drawing via POLYLINE.
C        12 July 1989    MDW     Namelist parameters WIDE and HIGH now specify
C                                the lengths of the x and y axes in inches in
C                                QMS plots. Default sizes have changed slightly.
C        28 Aug. 1989    MDW     Updated for run on DISSPLA version 11.0.  %REF
C                                is no longer used to pass character variables
C                                into DISSPLA routines.  Other minor changes
C                                were made to adapt to new version.
C        28 Dec. 1989    MDW     Untitled plots re-enabled for DISSPLA 11.0 
C                                version.  Enabled drawing of tick marks around 
C                                box using keyword GRID.  Introduced keyword 
C                                FORMAT to eliminate kludge used for specifying 
C                                plot format. READCOLS (formerly RDCOLS) revised
C                                to handle some exceptions at a higher level.  
C                                Two bug fixes in the process:  infinite loop 
C                                caused by flakey NAMELIST handling; and proper 
C                                handling of single column dataset even when
C                                YCOL is not 1. 
C        20 Feb. 1990    MDW     Introduced "#" in place of "$" to force blank
C                                titles or legends.  Eliminated use of "$" as a 
C                                string termination character. 
C         3 Mar. 1990    MDW     NAMELIST dictionary added for use by READCOLS.
C        15 Mar. 1990    MDW     Enabled reading of quoted numeric values and 
C                                blanks intended as titles through use of 
C                                STRIPPER.  (Replaces use of '#'.)
C        29 Oct. 1990    MDW     Increased maximum number of points from 10,000
C                                to 20,000.
C        10 Feb. 1991 D.Saunders Provided for choice of output devices via
C                                PLOTDEV (and optional command-line params.).
C                                QPLOT.COM (dealing with input and output files)
C                                is redundant now.  The frame-by-frame preview
C                                option is implemented as in program FMAP.
C                                PostScript from DISSPLA is usable at last,
C                                although only the FIRST frame of multi-frame
C                                plots prints on LaserWriters (%%EOFs at fault).
C                                Default color must be BLACK, not white, for
C                                PostScript. "qplot.out" now opened as 'unknown'
C                                for Unix compatibility.
C        22 Feb. 1991   "    "   Scale factors for true inches are no longer
C                                hard-coded: PLOTDEV handles them now.  Two
C                                lists of devices and the indirect DEVINDEX
C                                argument make life difficult - had to go
C                                to 2-D arrays for the lists (3-D for the
C                                hardware parameters).  Sorry!
C        25 Apr. 1991   "    "   Had to UPCASE any device type on the command
C                                line for Unix purposes.
C        23 May  1991   "    "   Arranged to force all strings to be done with
C                                stroke characters in the PostScript option
C                                if specified by /pss or -pss on the command
C                                line, as may be preferred if SOME strings have
C                                embedded directives.
C        31 May  1991   "    "   Restored system-dependent calls to EXIT at
C                                the end for use with a command file.
C        12 June 1991   "    "   Cray compiler doesn't like 'NML =' in read.
C        22 Sep. 1991   "    "   Rethought positioning of the plot on the page
C                                (allowing for height of legend, caption, etc.);
C                                provided for second Y axis: same tick marks
C                                but values scaled/shifted by user inputs;
C                                provided for essentially any number of caption
C                                lines; introduced XNOTE*, YNOTE*, NOTE*.
C        24 Sep. 1991  DAS/RAK   Introduced IDENT option.  Excess caption
C                                lines are now properly skipped (not fatal).
C        26 Sep. 1991    DAS     QUICK no longer has any legend limit, so the
C                                max. here has been raised from 20 to 30.
C                                QUICK's legend-positioning has been improved
C                                to suppress unused line-segment space, etc.
C        27 Sep. 1991  DAS/DBS   Applied X/YSHIFT/SCALE to X/YNOTE* entered
C                                in the same namelist; introduced X/YARROW[n].
C        04 Oct. 1991    DAS     Changed STRIPPER to return zero-length where
C                                previously it forced a minimum length of 1;
C                                made NOTEBOX = 1 the default, not 0.
C        09 Oct. 1991  RAK/DAS   Set default landscape height back to 4.5", not
C                                3.75"; DON'T suppress blank caption line if
C                                file/date identifier is the only caption line.
C        26 Oct. 1991    DAS     Rob found the identifier caption could show up
C                                as a title or label.  We ensure NUMTIT >= 4
C                                first now.  Couldn't resist: added ESCAPE and
C                                X/YNUMBERS options while I was at it.
C        17 Nov. 1991    DAS     Arranged for data-file indirection at the
C                                CURVE level (possibly nested).  First, reduced
C                                the functionality of READCOLS with a possible
C                                true keyword control scheme in mind: it now
C                                returns to the calling program upon finding
C                                any nonnumeric first token, without trying
C                                to identify it.  (At present, legal ones
C                                are END CURVE, END FRAME, @file, or namelist.)
C                                All back-spacing has been eliminated except
C                                for the unavoidable one prior to reading a
C                                namelist.  NLCHEK's functionality is also
C                                reduced: it no longer reads data - it just
C                                checks the buffer, and is called in two places,
C                                not three.  Eliminating namelist altogether
C                                should now be straightforward when the time
C                                comes.
C        19 Nov. 1991    DAS     Clarified likely usages of indirect files,
C                                requiring suppression of the curve control
C                                reinitialization in the special case of
C                                indirection (or EOF or another namelist)
C                                encountered IMMEDIATELY following a namelist.
C        01 Jan. 1992    DAS     Enabled alternative output metafile names.
C        08 Jan. 1992    DAS     Modularized some code as utility PARSENAM.
C        04 May  1993    DAS     Added option to rotate any curve about a point.
C        14 Sep. 1997    DAS     HEIGHT & WIDTH need restoring to nominal values
C                                after a preview before replotting to a metafile
C                                with different hardware calibration factors.
C                                WORSE: the factors should be applied AFTER
C                                calling LAYOUT, not before, because LAYOUT can
C                                update WIDTH for scale plots.  Also: no need
C                                to call LAYOUT again when replotting.
C                                If XCOL or YCOL = 0 for all curves on a frame,
C                                default to integer numbering.
C                                DATE_AND_TIME avoids the year 2000 problem.
C        17 Sep. 1997    DAS     Max. # curves: 5,000; max # points: 500,000.
C        29 July 2002    DAS     Moved open of qplot.out on LWRITE=2 to the top
C                                for compatibility with DISSPLA's SETDEV in
C                                PLOTDEV.  However, DONEPL still writes to
C                                fort.2 instead of qplot.out if DEBUG=.TRUE.
C        31 July, 2002   DAS     Introduced SPLINE keyword to allow use of
C                                nonparametric splines (Y vs. X or X vs. Y).
C
C-------------------------------------------------------------------------------


C     Declarations.
C     -------------

      IMPLICIT NONE

C     Local constants.

      INTEGER
     >   LENDIC, LKEYBD, LMETA, LREAD1, LSCREEN, LWRITE, MAXCUR, MAXDEV,
     >   MAXLEG, MAXLEN, MAXLIST, MAXNAME, MAXNOTE, MAXPT, MAXTIT,
     >   NCLDIC, NFDIC, NLNDIC, NMENU, NNMLDIC, NNUMDIC, NORDIC,
     >   NSPLDIC, NSTDIC

      PARAMETER
     >  (LENDIC = 20,             ! Length (max.) of entries, all dictionaries
     >   LKEYBD = 5,              ! Logical unit number for keyboard input
     >   LMETA  = 3,              ! Logical unit number for plot metafile output
     >   LREAD1 = 7,              ! Logical unit number for base data file
     >   LSCREEN= 6,              ! Logical unit for prompts and plot previews
     >   LWRITE = 2,              ! Logical unit number for diagnostic messages
     >   MAXCUR = 5000,           ! Max. number of curves
     >   MAXDEV = 2,              ! Max. number of active output devices
     >   MAXLEG = 30,             ! Max. number of legend strings
     >   MAXLEN = 132,            ! Max. length of titles, subtitles, etc.
     >   MAXLIST= 5,              ! Max. number of screen OR metafile options
     >   MAXNAME= 20,             ! Max. length of output metafile name
     >   MAXNOTE= 30,             ! Max. number of annotated points per frame
     >   MAXPT  = 200000,         ! Max. number of data points
     >   MAXTIT = 30,             ! Room for 2 titles, 2 labels, and the rest
                                  ! captions, 2 of which may be the optional
                                  ! blank line plus file/date identifier
     >   NCLDIC = 12,             ! Number of entries in COLOR dictionary
     >   NFDIC  = 8,              ! Ditto for FIT
     >   NLNDIC = 16,             ! Ditto for LINE types
     >   NMENU  = 4,              ! Ditto for plotting option menu
     >   NNMLDIC= 2,              ! Ditto for NAMELIST keywords
     >   NNUMDIC= 3,              ! Ditto for axis number types
     >   NORDIC = 7,              ! Ditto for plot ORIENTation
     >   NSPLDIC= 4,              ! Ditto for SPLINE types
     >   NSTDIC = 20)             ! Ditto for text STYLEs

      REAL
     >   FLAG, HLAND, HPORT, WLAND, WPORT, ONE, ZERO
      PARAMETER
     >  (FLAG   = 999.E+0,        ! Preset default value for input REALs
     >   HLAND  = 4.5E+0,         ! Default plot frame height & width,
     >   WLAND  = 7.5E+0,         ! landscape & portrait, in inches
     >   HPORT  = 5.0E+0,
     >   WPORT  = 5.0E+0,
     >   ONE    = 1.E+0,
     >   ZERO   = 0.E+0)

      CHARACTER
     >   BLANK, CFLAG * 7, EXCLAM, NAMLIST * 7, QUOTES * 2
      PARAMETER
     >  (BLANK  = ' ',
     >   CFLAG  = 'DEFAULT',
     >   EXCLAM = '!',
     >   NAMLIST= 'OPTIONS',
     >   QUOTES = '''"')

C     Variables.

      INTEGER
     >   COLDAT (MAXCUR), CURVE, DEFSYM, DEVINDEX (MAXDEV), DEVINDX,
     >   DEVLIST (MAXLIST, 2), ENDSTATE, FIRST, FIRST2, FIRST3,
     >   FITDAT (MAXCUR), FORMAT, FRAME, GRID, I, IDATE_TIME(8),
     >   IER, IOS, ITEM, LAST, LAST2, LAST3, LDEBUG, LEGINF (3, MAXLEG),
     >   LEG, LENBUF, LENLEG (MAXLEG), LENM, LENPRM, LENTIT (MAXTIT+2),
     >   LINDAT (MAXCUR), LREAD, MARK, NDEV (MAXDEV), NOTEBOX,
     >   NPT (MAXCUR), NTOKENS, NUMNOTE, NUMTIT, OLDCUR, POINTER, RCOL,
     >   SPLDAT (MAXCUR), SYMBOL, SYMDAT (MAXCUR), THCOL, TOTAL,
     >   XCOL, XCOLUMN, YCOL, YCOLUMN
      REAL
     >   ANGLE, CENTERX, CENTERY, DEFHGT, DEFWID, DEG2RAD, DEGREES,
     >   HEIGHT, HSAVE, HWPARAMS (MAXLIST, 2, 2), RADIANS,
     >   SCALEX, SCALEY, SHIFTX, SHIFTY, RMIN, RSCALE, RSHIFT, SCALETH,
     >   SCALER, SHIFTR, SHIFTTH, RSTEP, THSCALE, THSHIFT, WIDTH, WSAVE,
     >   X (MAXPT), XCENTER, XINC, XLEFT, XMAX, XMIN, XRIGHT, XSCALE,
     >   XSCLFR, XSHIFT, XSTEP,
     >   XARROW,  XARROW1, XARROW2, XARROW3, XARROW4, XARROW5, XARROW6,
     >   XARROW7, XARROW8, XARROW9, XARROWQ (0 : 9), XARROWS (MAXNOTE),
     >   XNOTE,  XNOTE1, XNOTE2, XNOTE3, XNOTE4, XNOTE5, XNOTE6, XNOTE7,
     >   XNOTE8, XNOTE9, XNOTEQ (0 : 9), XNOTES (MAXNOTE),
     >   Y (MAXPT), Y2SCALE, Y2SHIFT, YBOT, YBOTTOM, YCENTER, YINC,
     >   YMAX, YMIN, YSCALE, YSCLFR, YSHIFT, YSTEP, YTOP,
     >   YARROW,  YARROW1, YARROW2, YARROW3, YARROW4, YARROW5, YARROW6,
     >   YARROW7, YARROW8, YARROW9, YARROWQ (0 : 9), YARROWS (MAXNOTE),
     >   YNOTE,  YNOTE1, YNOTE2, YNOTE3, YNOTE4, YNOTE5, YNOTE6, YNOTE7,
     >   YNOTE8, YNOTE9, YNOTEQ (0 : 9), YNOTES (MAXNOTE)
      LOGICAL
     >   ALL, ALLXCOL0, ALLYCOL0, ALPHA, BADDEV, CR, DATANAME, DEBUG,
     >   ENDFRAME, EOF, ERROR, FINIS, FLAKEY, GRAPHICS, IDENT, INTEGERX,
     >   INTEGERY, META, NEWCNTRL, NEWFILE, NOFIT, NOSPLINE, OUTNAME,
     >   PRESENT, PREVIEW, QUIT, REPLOT, SCREEN, STROKE, YES
      CHARACTER
     >   BUFFERS (2) * (MAXLEN),          CLASS * (MAXLEN),
     >   CDATE_TIME (3) * 10,
     >   COLDIC (NCLDIC) * (LENDIC),      COLOR * (LENDIC),
     >   COMMENT,                         CTEMP * 5,
     >   DATAFILE * (MAXLEN),             DEVTYPE (MAXDEV) * 3,
     >   ESCAPE * 2,                      EXPLAN * (MAXLEN),
     >   FILENAME * 9,
     >   FIT * 8,                         FITDIC (NFDIC) * (LENDIC),
     >   LEGEND * (MAXLEN),               LEGTXT (MAXLEG) * (MAXLEN),
     >   LINDIC (NLNDIC) * (LENDIC),      LINE * (LENDIC),
     >   METAFILE * (MAXNAME),            NMLDIC (NNMLDIC) * (LENDIC),
     >   NOTE  * (MAXLEN),
     >   NOTE1 * (MAXLEN),                NOTE2 * (MAXLEN),
     >   NOTE3 * (MAXLEN),                NOTE4 * (MAXLEN),
     >   NOTE5 * (MAXLEN),                NOTE6 * (MAXLEN),
     >   NOTE7 * (MAXLEN),                NOTE8 * (MAXLEN),
     >   NOTE9 * (MAXLEN),                NOTEQ (0 : 9) * (MAXLEN),
     >   NOTES (MAXNOTE) * (MAXLEN),
     >   NUMDIC (NNUMDIC) * 7,
     >   OPT * 1,                         ORIENT * 9,
     >   ORDIC (NORDIC) * (LENDIC),       PARAMS * 80,
     >   PLOT * 12,                       PLOTMENU (NMENU) * 42,
     >   PROGNAME * 5,                    PROMPT * 43,
     >   SEPS * 2,
     >   SPLDIC (NSPLDIC) * (LENDIC),     SPLINE * (LENDIC),
     >   STYDIC (NSTDIC) * (LENDIC),      STYLE * (LENDIC),
     >   TITLES (MAXTIT+2) * (MAXLEN),
     >   TOKEN (2) * 8,      ! Just enough for 'ENDCURVE' at this stage.
     >   TYPE * 3,                        TYPELIST (MAXLIST, 2) * 3,
     >   XNUMBERS * 7,                    YNUMBERS * 7,
     >   Y2LABEL * (MAXLEN)

C     Storage.  The growing list of EQUIVALENCEs, together with "redundant"
C     entries in the NAMELIST, provides a rudimentary synonym capability for
C     the variables in $OPTIONS.  CLASS and EXPLAN are included here because
C     they are read with $OPTIONS but passed as part of the TITLES array.
C     The NOTE* equivalences are forced by systems which do not support
C     CHARACTER arrays in namelists.

      EQUIVALENCE
     >   (ANGLE,   DEGREES),
     >   (CENTERX, XCENTER),
     >   (CENTERY, YCENTER),
     >   (CLASS,   TITLES (MAXTIT + 1)),
     >   (EXPLAN,  TITLES (MAXTIT + 2)),
     >   (NOTE,    NOTEQ (0)),
     >   (NOTE1,   NOTEQ (1)), (NOTE2, NOTEQ (2)), (NOTE3, NOTEQ (3)),
     >   (NOTE4,   NOTEQ (4)), (NOTE5, NOTEQ (5)), (NOTE6, NOTEQ (6)),
     >   (NOTE7,   NOTEQ (7)), (NOTE8, NOTEQ (8)), (NOTE9, NOTEQ (9)),
     >   (RCOL,    YCOL),
     >   (RMIN,    YMIN),
     >   (RSCALE,  YSCALE),
     >   (RSHIFT,  YSHIFT),
     >   (RSTEP,   YSTEP),
     >   (SCALER,  RSCALE),
     >   (SCALETH, THSCALE),
     >   (SCALETH, XSCALE),
     >   (SCALEX,  XSCALE),
     >   (SCALEY,  YSCALE),
     >   (SHIFTR,  RSHIFT),
     >   (SHIFTTH, THSHIFT),
     >   (SHIFTTH, XSHIFT),
     >   (SHIFTX,  XSHIFT),
     >   (SHIFTY,  YSHIFT),
     >   (THCOL,   XCOL),
     >   (XARROW,  XARROWQ (0)),
     >   (XARROW1, XARROWQ (1)), (XARROW2, XARROWQ (2)),
     >   (XARROW3, XARROWQ (3)), (XARROW4, XARROWQ (4)),
     >   (XARROW5, XARROWQ (5)), (XARROW6, XARROWQ (6)),
     >   (XARROW7, XARROWQ (7)), (XARROW8, XARROWQ (8)),
     >   (XARROW9, XARROWQ (9)),
     >   (XCOLUMN, XCOL),
     >   (XLEFT,   XMIN),
     >   (XRIGHT,  XMAX),
     >   (XSTEP,   XINC),
     >   (XNOTE,   XNOTEQ (0)),
     >   (XNOTE1,  XNOTEQ (1)), (XNOTE2, XNOTEQ (2)),
     >   (XNOTE3,  XNOTEQ (3)), (XNOTE4, XNOTEQ (4)),
     >   (XNOTE5,  XNOTEQ (5)), (XNOTE6, XNOTEQ (6)),
     >   (XNOTE7,  XNOTEQ (7)), (XNOTE8, XNOTEQ (8)),
     >   (XNOTE9,  XNOTEQ (9)),
     >   (YARROW,  YARROWQ (0)),
     >   (YARROW1, YARROWQ (1)), (YARROW2, YARROWQ (2)),
     >   (YARROW3, YARROWQ (3)), (YARROW4, YARROWQ (4)),
     >   (YARROW5, YARROWQ (5)), (YARROW6, YARROWQ (6)),
     >   (YARROW7, YARROWQ (7)), (YARROW8, YARROWQ (8)),
     >   (YARROW9, YARROWQ (9)),
     >   (YCOLUMN, YCOL),
     >   (YBOT,    YMIN),
     >   (YBOTTOM, YMIN),
     >   (YSTEP,   YINC),
     >   (YTOP,    YMAX),
     >   (YNOTE,   YNOTEQ (0)),
     >   (YNOTE1,  YNOTEQ (1)), (YNOTE2, YNOTEQ (2)),
     >   (YNOTE3,  YNOTEQ (3)), (YNOTE4, YNOTEQ (4)),
     >   (YNOTE5,  YNOTEQ (5)), (YNOTE6, YNOTEQ (6)),
     >   (YNOTE7,  YNOTEQ (7)), (YNOTE8, YNOTEQ (8)),
     >   (YNOTE9,  YNOTEQ (9))

      NAMELIST /OPTIONS/
     >   ANGLE, CENTERX, CENTERY, CLASS, COLOR, COMMENT, DEBUG, DEGREES,
     >   ESCAPE, EXPLAN, FIT, FORMAT, GRID, HEIGHT, IDENT, LEGEND, LINE,
     >   NOTE, NOTE1, NOTE2, NOTE3, NOTE4, NOTE5, NOTE6, NOTE7, NOTE8,
     >   NOTE9, NOTEBOX, ORIENT, PLOT, RADIANS, RCOL, RMIN, RSCALE,
     >   RSHIFT, RSTEP, SCALER, SCALETH, SCALEX, SCALEY, SHIFTR,
     >   SHIFTTH, SHIFTX, SHIFTY, SPLINE, STYLE, SYMBOL, THCOL, THSCALE,
     >   THSHIFT, WIDTH,
     >   XARROW,  XARROW1, XARROW2, XARROW3, XARROW4, XARROW5, XARROW6,
     >   XARROW7, XARROW8, XARROW9,
     >   XCENTER, XCOL, XCOLUMN, XINC, XLEFT, XMAX, XMIN,
     >   XNOTE,  XNOTE1, XNOTE2, XNOTE3, XNOTE4, XNOTE5, XNOTE6, XNOTE7,
     >   XNOTE8, XNOTE9, XNUMBERS, XRIGHT, XSCALE, XSHIFT, XSTEP,
     >   YARROW,  YARROW1, YARROW2, YARROW3, YARROW4, YARROW5, YARROW6,
     >   YARROW7, YARROW8, YARROW9,
     >   Y2LABEL, Y2SCALE, Y2SHIFT, YBOT, YBOTTOM, YCENTER, YCOL,
     >   YCOLUMN, YINC, YMAX, YMIN,
     >   YNOTE,  YNOTE1, YNOTE2, YNOTE3, YNOTE4, YNOTE5, YNOTE6, YNOTE7,
     >   YNOTE8, YNOTE9, YNUMBERS, YSCALE, YSHIFT, YSTEP, YTOP

      DATA
     >   PROGNAME /'qplot'/,
     >   PROMPT   /'          *** OPTIONS:  Frame number nn ***'/

      DATA (PLOTMENU (I), I = 1, NMENU)
     >  /' P - Preview the plot (metafile optional).',
     >   ' M - Metafile output exclusively.         ',
     >   ' S - Skip to the next frame.              ',
     >   ' Q - Quit from QPLOT.                     '/

      DATA (COLDIC (I), I = 1, NCLDIC)
     >  /'BLACK',
     >   'MAGENTA',
     >   'RED',
     >   'YELLOW',
     >   'GREEN',
     >   'CYAN',
     >   'BLUE',
     >   'WHITE',
     >   'DEFAULT=BLACK',            ! Needed for PostScript; OK for Tek, DIP
     >   'NEUTRAL=WHITE',
     >   'INVISIBLE=BLACK',          ! ? Not for PS
     >   'NONE=BLACK'/               !   "    "

      DATA (FITDIC (I), I = 1, NFDIC)
     >  /'LINEAR',
     >   'MONOTONE',
     >   'BESSEL',
     >   'CLOSED',
     >   'DEFAULT=LINEAR',
     >   'LOOSE=BESSEL',
     >   'STRAIGHT=LINEAR',
     >   'TIGHT=MONOTONE'/

      DATA (LINDIC (I), I = 1, NLNDIC)
     >  /'CONNECTED',
     >   'SYMBOLS',
     >   'SOLID',
     >   'DOTS',
     >   'DASHES',
     >   'CHAINDOTS',
     >   'CHAINDASHES',
     >   'LONGDASHES',
     >   'THICK',
     >   'BOTH=CONNECTED',
     >   'DEFAULT=CONNECTED',
     >   'FAT=THICK',
     >   'NONE=SYMBOLS',
     >   'NULL=SYMBOLS',
     >   'STRAIGHT=CONNECTED',
     >   'SHORTDASH=DASHES'/

      DATA (NMLDIC (I), I = 1, NNMLDIC)   ! This will become a control
     >  /'$OPTIONS',                      ! keyword dictionary some day.
     >   '&OPTIONS'/

      DATA (NUMDIC (I), I = 1, NNUMDIC)
     >  /'REAL   ',
     >   'INTEGER',
     >   'WHOLE  '/

      DATA (ORDIC (I), I = 1, NORDIC)
     >  /'COMIC=PORTRAIT',
     >   'DEFAULT=PORTRAIT',
     >   'HORIZONTAL=LANDSCAPE',
     >   'LANDSCAPE',
     >   'MOVIE=LANDSCAPE',
     >   'PORTRAIT',
     >   'VERTICAL=PORTRAIT'/
        
      DATA (SPLDIC (I), I = 1, NSPLDIC)
     >  /'STANDARD',
     >   'REVERSED',
     >   'PARAMETRIC',
     >   'DEFAULT=PARAMETRIC'/

      DATA (STYDIC (I), I = 1, NSTDIC)
     >  /'DISSPLA',
     >   'CARTOG',
     >   'SIMPLX',
     >   'SCMPLX',
     >   'COMPLX',
     >   'DUPLX',
     >   'TRIPLX',
     >   'GOTHIC',
     >   'FUTURA',
     >   'SERIF',
     >   'FASHON',
     >   'LOGO1',
     >   'SWISSL',
     >   'SWISSM',
     >   'SWISSB',
     >   'DEFAULT=SCMPLX',
     >   'STICK=DISSPLA',
     >   'LIGHT=SIMPLX',
     >   'TIMESROMAN=SERIF',
     >   'HELVETICA=SWISSL'/


C     Execution.
C     ----------

      LREAD    = LREAD1        ! Base LUN for input data (add 1 for each @)
      ENDSTATE = 0             ! For normal completion
      GRAPHICS = .FALSE.       ! No graphics done yet
      SEPS = BLANK // CHAR (9) ! Blank & tab are the only insignificant chars.

C     Open the diagnostics file.  Unix forces overwriting of any existing one.

      FILENAME = 'qplot.out'
      CALL OPENER (LSCREEN, BLANK, LKEYBD, FILENAME, LWRITE, 'UNKNOWN')
      WRITE (LWRITE, 1000)

C     Obtain the list of preview (screen) devices supported at this site for
C     this application (via the file PROGNAME.config in a standard location).

      NDEV (1) = MAXLIST
      CALL PLOTDEV ('P', PROGNAME, LSCREEN, LKEYBD, LMETA, BLANK,
     >   NDEV (1), DEVLIST (1, 1), TYPELIST (1, 1), HWPARAMS (1, 1, 1),
     >   DEVINDEX)

C     Obtain the list of metafiles supported.

      NDEV (2) = MAXLIST
      CALL PLOTDEV ('M', PROGNAME, LSCREEN, LKEYBD, LMETA, BLANK,
     >   NDEV (2), DEVLIST (1, 2), TYPELIST (1, 2), HWPARAMS (1, 1, 2),
     >   DEVINDEX)


C     Check the command line arguments for plot device(s) and in/out file names.
C     --------------------------------------------------------------------------

C     (This clutter should probably be moved out of the main program some day.)

      CALL GETPRM (PARAMS, LENPRM)

      DATANAME = .FALSE.
      OUTNAME  = .FALSE.
      SCREEN   = .FALSE.
      META     = .FALSE.
      BADDEV   = .FALSE.
      MARK     = -1
   10 CONTINUE
         FIRST = MARK + 2
         IF (FIRST .LE. LENPRM) THEN    ! Isolate a parameter

            CALL SCAN2 (PARAMS, SEPS, FIRST, LENPRM, MARK)

            IF (PARAMS (FIRST : FIRST) .NE. '/' .AND.
     >          PARAMS (FIRST : FIRST) .NE. '-') THEN

               IF (.NOT. DATANAME) THEN  ! Must be the input datafile name
                  DATANAME = .TRUE.
                  DATAFILE = PARAMS (FIRST : MARK)

               ELSE                      ! Must be the metafile name or 'same'.
                  OUTNAME  = .TRUE.      ! Can't really trap an unintended
                  METAFILE = PARAMS (FIRST : MARK)  ! extension because its
                  LENM = MARK - FIRST + 1           ! OK under Unix and there
               END IF                    ! are too many possible metafiles.

            ELSE                         ! Check for a screen device
               TYPE = PARAMS (FIRST + 1 : MARK)
               CALL UPCASE (TYPE)
               CALL LOOKUP (NDEV (1), TYPELIST (1, 1), .FALSE., TYPE, I)
               IF (I .GT. 0) THEN
                  DEVINDEX (1) = I
                  DEVTYPE (1) = TYPE
                  SCREEN = .TRUE.
               ELSE                      ! Check for a metafile type
		  STROKE = TYPE (1 : 3) .EQ. 'PSS' ! Kludge for forcing stroke
		  IF (STROKE) TYPE (3 : 3) = BLANK ! characters for all text

                  CALL LOOKUP (NDEV (2), TYPELIST (1, 2), .FALSE., TYPE, 
     >               I)
                  IF (I .GT. 0) THEN
                     DEVINDEX (2) = I
                     DEVTYPE (2) = TYPE
                     META = .TRUE.
                  ELSE
                     WRITE (LSCREEN, 1020)
     >               ' >> Invalid or ambiguous plot device.  Reenter.'
                     BADDEV = .TRUE.
                  END IF
               END IF
            END IF

            GO TO 10                            ! Look for another parameter

         END IF        

          
C     Open the data file.  (May still need to prompt for its name.)
C     -------------------------------------------------------------

      IF (DATANAME) THEN
         CALL OPENER (LSCREEN, BLANK, LKEYBD, DATAFILE, LREAD, 'OLD')
      ELSE
         WRITE (LSCREEN, 1020)
         DATAFILE = 'qplot.dat'
         CALL OPENER (LSCREEN, 'Data file? (<CR> = qplot.dat) ',
     >      LKEYBD, DATAFILE, LREAD, 'OLD')
      END IF


C     Parse the data file name for possible use in a metafile name or an
C     identifying caption.  Any subdirectory or version number is suppressed.
C     End-of-extension is not needed here, so make it end-of-version (if any).

      CALL PARSENAM (DATAFILE, FIRST2, LAST2, FIRST3, LAST3) ! Delimits name/ext

      LAST3 = MAX (LAST3, INDEX (DATAFILE, BLANK) - 1)   ! Used below for ident.


C     Default the metafile name if necessary.  A prompt to override it would
C     be tiresome for those who fail to use command line arguments, so don't
C     bother.

      IF (.NOT. OUTNAME .OR. LAST2 .EQ. 0) THEN
         METAFILE = PROGNAME    ! I.e., 'qplot'
         LENM = 5
      ELSE          ! Check for 'same' entered as the metafile name
         CTEMP = METAFILE
         CALL UPCASE (CTEMP)
         IF (CTEMP .EQ. 'SAME ') THEN
            METAFILE = DATAFILE (FIRST2 : LAST2)
            LENM = LAST2 - FIRST2 + 1
         END IF
      END IF


C     Prompt only if neither preview device nor metafile type was specified
C     on the command line, or if a bad device was encountered.

      IF ((.NOT. SCREEN .AND. .NOT. META) .OR. BADDEV) THEN
         IF (.NOT. SCREEN) THEN
            CALL PLOTDEV ('S', PROGNAME, LSCREEN, LKEYBD, LMETA,
     >         'Preview device selection: ("None" for metafile only.)',
     >         NDEV (1), DEVLIST (1, 1), TYPELIST (1, 1),
     >         HWPARAMS (1, 1, 1), DEVINDEX (1))
            SCREEN = DEVINDEX (1) .NE. 0
         END IF

         IF (.NOT. META) THEN
            CALL PLOTDEV ('S', METAFILE (1:LENM), LSCREEN, LKEYBD,
     >         LMETA, 'Metafile selection:', NDEV (2), DEVLIST (1, 2),
     >         TYPELIST (1, 2), HWPARAMS (1, 1, 2), DEVINDEX (2))
            META = DEVINDEX (2) .NE. 0
         END IF

         IF (.NOT. SCREEN .AND. .NOT. META) GO TO 985    ! Early exit.
      END IF


C     Read and plot data frame-by-frame.
C     ----------------------------------

      ALL     = META .AND. .NOT. SCREEN  ! T = do rest of frames (no prompts)
      PREVIEW = SCREEN
      DEBUG   = .FALSE.          ! Leave it set once it is set
      IDENT   = .FALSE.          ! Ditto
      FINIS   = .FALSE.          ! Not EOF yet
      DEG2RAD = ASIN (ONE) / 90.
      FRAME   = 0

  100 CONTINUE
         FRAME = FRAME + 1
         IF (FINIS) GO TO 990   ! Normal end where EOF is already known of

         ENDFRAME = .FALSE.

C        (Re)set default parameters for each frame.

         NUMTIT   = 0
         NUMNOTE  = 0
         LEG      = 0
         DEFSYM   = 0
         NOTEBOX  = 1
         FORMAT   = 1
         GRID     = -1
         HEIGHT   = FLAG
         WIDTH    = FLAG
         XMIN     = FLAG
         XMAX     = FLAG
         XSTEP    = FLAG
         YMIN     = FLAG
         YMAX     = FLAG
         YSTEP    = FLAG
         STYLE    = CFLAG
         ORIENT   = CFLAG
         PLOT     = CFLAG
         ESCAPE   = '[]'
         CLASS    = BLANK
         EXPLAN   = BLANK
         Y2LABEL  = BLANK
         Y2SCALE  = ZERO
         Y2SHIFT  = ZERO
         XNUMBERS = CFLAG
         YNUMBERS = CFLAG

         DO 110, I = 1, MAXTIT + 2  ! CLASS and EXPLAN are not affected by this.
            LENTIT (I) = 0
  110    CONTINUE

C        Look for title, subtitle, axis labels and captions.
C        ---------------------------------------------------

  200    CONTINUE

            CALL GETLINE (LREAD, EXCLAM, BUFFERS (1), LENBUF, IOS)

            IF (IOS .LT. 0) THEN         ! EOF
               IF (FRAME .GT. 1 .AND. NUMTIT .EQ. 0) THEN

C                 Terminate normally.

                  GO TO 990
               ELSE                      ! EOF, but no plottable data

                  WRITE (LWRITE, 1010) FRAME
                  GO TO 980
               END IF

            ELSE IF (IOS .NE. 0) THEN    ! System-dependent read error

               WRITE (LWRITE, 1015) IOS
               GO TO 980
            END IF

            IF (LENBUF .EQ. 0) GO TO 200   ! Empty lines are ignored


C           Check for the optional NAMELIST.

            CALL NLCHEK (BUFFERS (1), LENBUF, NAMLIST, PRESENT, FLAKEY)

            IF (FLAKEY) THEN

C              A namelist appears to be present, but the format is wrong
C              (name is incorrect, or doesn't begin in column 2).

               WRITE (LWRITE, 1050) CURVE + 1, FRAME
               CALL NLGRIPE (LWRITE, NAMLIST)
               GO TO 980

            ELSE IF (PRESENT) THEN   ! Do the unavoidable backspace at the read

!!!!!!!!!!!!   BACKSPACE (LREAD)

            ELSE                     ! Must be text or numerical data.

               IF (ALPHA (BUFFERS (1) (1 : LENBUF))) THEN

C                 It's text.  Add it to the list of titles and look for more,
C                 unless it's an indirect file, which can be handled after
C                 the call to READCOLS.

                  FIRST = 1                  
                  CALL SCAN2 (BUFFERS (1), SEPS, FIRST, LENBUF, MARK)

                  IF (BUFFERS (1) (FIRST : FIRST) .NE. '@') THEN

                     IF (NUMTIT .LT. MAXTIT) THEN

C                       Strip quotes from numeric data or blanks intended to be
C                       used as titles or labels.

                        CALL STRIPPER (BUFFERS (1), LENBUF, QUOTES)

                        NUMTIT = NUMTIT + 1
                        LENTIT (NUMTIT) = LENBUF
                        TITLES (NUMTIT) (1 : LENBUF) =
     >                      BUFFERS (1) (1 : LENBUF)
                     ELSE
                        WRITE (LWRITE, 1090) MAXTIT
                     END IF

                     GO TO 200      ! Go back for the next input line.
                  ELSE
                     ! It's an indirect file name.  Pass it on ...
                  END IF

               ELSE
C                 It's numeric.  Pass it on ...
               END IF

            END IF               ! End of check for namelist.

C        End of loop over frame-level text lines.


C        Read curves to be plotted until the plot frame is complete.
C        -----------------------------------------------------------

         CURVE = 0
         ALLXCOL0 = .TRUE.  ! To enable defaulting to integer axes
         ALLYCOL0 = .TRUE.

  400    CONTINUE

C           Re(set) default plot parameters for each curve.

            COLOR   = CFLAG
            FIT     = CFLAG
            SPLINE  = CFLAG
            LINE    = CFLAG
            LEGEND  = BLANK
            SCALEX  = ONE
            SCALEY  = ONE
            SHIFTX  = ZERO
            SHIFTY  = ZERO
            CENTERX = ZERO
            CENTERY = ZERO
            DEGREES = ZERO
            RADIANS = ZERO
            SYMBOL  = -1
            XCOL    = 1
            YCOL    = 2

            DO 410, I = 0, 9
               NOTEQ (I)   = BLANK
               XNOTEQ (I)  = FLAG
               YNOTEQ (I)  = FLAG
               XARROWQ (I) = FLAG
               YARROWQ (I) = FLAG
  410       CONTINUE


C           Any namelist has already been detected at the end of the titles
C           or at the end of a curve.
           
            IF (PRESENT) THEN

               PRESENT = .FALSE.
               NEWCNTRL = .TRUE.
               LENBUF = 0       ! BUFFERS (1) is not valid after the next read.

               BACKSPACE (LREAD)

C              Read the namelist.
C              ------------------

C*****         READ (LREAD, OPTIONS, IOSTAT = IOS)  ! Not allowed on Cray
               IOS = -1
               READ (LREAD, OPTIONS, ERR = 420)
C*****         READ (LREAD, NML=OPTIONS, ERR = 420) ! Reqd. on IRIS; bad on Cray
               IOS = 0
  420          IF (IOS .NE. 0) THEN
                  WRITE (LWRITE, 1050) CURVE + 1, FRAME
                  CALL NLGRIPE (LWRITE, NAMLIST)
                  GO TO 980
               END IF

            ELSE
               NEWCNTRL = .FALSE.
            END IF


  430       CONTINUE

C           Read a curve.
C           -------------

C           Note that the packing details are (mostly) hidden.  The curve
C           counter will only be incremented if some plottable data was found.

            OLDCUR = CURVE

            CALL READCOLS (LREAD, LWRITE, XCOL, YCOL, MAXPT, NPT,
     >         MAXCUR, CURVE, X, Y, BUFFERS, LENBUF, FINIS, IER)

            ERROR = IER .NE. 0
            NEWFILE = .FALSE.

            IF (IER .EQ. +2) THEN

C              Illegal input data found.
             
               WRITE (LWRITE, 1080)
               WRITE (LWRITE, 1030) BUFFERS (1) (1:60)
               WRITE (LWRITE, 1040) BUFFERS (2) (1:60)
               WRITE (LWRITE, 1020)

            ELSE IF (IER .EQ. +1) THEN

C              READCOLS encountered a curve-ending keyword.
C              Check to see if it's part of a namelist.

  440          CALL NLCHEK (BUFFERS (1), LENBUF, NAMLIST, PRESENT,
     >                      FLAKEY)

               IF (PRESENT) THEN       ! Process it as part of next curve.

                  ERROR = .FALSE.

               ELSE IF (FLAKEY) THEN   ! Gripe. (ERROR = .TRUE. ends frame.)

                  CALL NLGRIPE (LWRITE, NAMLIST)

               ELSE

C                 We need to tokenize the keyword(s).

                  NTOKENS = 2
                  CALL TOKENS (BUFFERS (1) (1 : LENBUF), NTOKENS, TOKEN)

                  IF (TOKEN (1) (1:3) .EQ. 'END') THEN
                     IF (TOKEN (2) (1:5) .EQ. 'CURVE' .OR.
     >                   TOKEN (1) (4:8) .EQ. 'CURVE') THEN

                        ERROR = .FALSE.
                        LENBUF = 0

                     ELSE IF (TOKEN (2) (1:5) .EQ. 'FRAME' .OR.
     >                        TOKEN (1) (4:8) .EQ. 'FRAME') THEN
                        ENDFRAME = .TRUE.
                        ERROR = .FALSE.
                        LENBUF = 0
                     END IF

                  ELSE IF (TOKEN (1) (1:1) .EQ. '@') THEN

C                     Go to a(nother) level of indirection.

                      I = INDEX (BUFFERS (1) (1 : LENBUF), '@') + 1

                      LREAD = LREAD + 1
                      CALL OPENER (LSCREEN, BLANK, LKEYBD,
     >                   BUFFERS (1) (I : LENBUF), LREAD, 'OLD')
                      LENBUF = 0
                      ERROR = .FALSE.
                      NEWFILE = .TRUE. ! Needed for special case (NEWCNTRL=T).
                  ELSE
C                    Not a namelist, END CURVE, END FRAME, or file name.
C                    Handle true keyword controls here some day.

                     WRITE (LWRITE, 1020) ' Bad keyword encountered.'
                     WRITE (LWRITE, 1030) BUFFERS (1) (1:60)
                     LENBUF = 0
                  END IF

               END IF

C              Special case where we DON'T reinitialize curve controls:

               IF (NEWCNTRL .AND. NEWFILE) THEN   ! No data yet
                  NEWCNTRL = .FALSE.
                  GO TO 430
               END IF

            END IF

C           Back up a level of indirection?

            IF (FINIS) THEN                ! EOF for current data file
               IF (LREAD .GT. LREAD1) THEN
                  CLOSE (UNIT = LREAD)
                  LREAD = LREAD - 1
                  FINIS = .FALSE.
               END IF
               LENBUF = 0
            END IF


C           If a curve was accepted, examine or apply related controls.
C           -----------------------------------------------------------

            IF (CURVE .GT. OLDCUR) THEN

               IF (XCOL .NE. 0) ALLXCOL0 = .FALSE.
               IF (YCOL .NE. 0) ALLYCOL0 = .FALSE.

C              Apply any specified rotation.

               IF (RADIANS .EQ. ZERO) RADIANS = DEGREES * DEG2RAD

               IF (RADIANS .NE. ZERO)
     >            CALL ROTATE (MAXPT, NPT, CURVE, X, Y, RADIANS,
     >                         CENTERX, CENTERY)

C              Apply specified linear transformations.

               IF (SCALEX .NE. ONE .OR. SHIFTX .NE. ZERO)
     >            CALL RESCALE (MAXPT, NPT, CURVE, X, SCALEX, SHIFTX)

               IF (SCALEY .NE. ONE .OR. SHIFTY .NE. ZERO)
     >            CALL RESCALE (MAXPT, NPT, CURVE, Y, SCALEY, SHIFTY)

C              Convert color string to uppercase and decode.  Note that
C              we only attempt to call LOOKUP twice here (and in what
C              follows) to eliminate even the possibility of looping.

               CALL UPCASE (COLOR)
               CALL LOOKUP (NCLDIC, COLDIC, .FALSE., COLOR, ITEM)
               IF (ITEM .LE. 0) THEN
                  WRITE (LWRITE, 1070) CURVE, FRAME
                  COLOR = CFLAG
                  CALL LOOKUP (NCLDIC, COLDIC, .FALSE., COLOR, ITEM)
               END IF
               IF (ITEM .EQ. 8) ITEM = 0
               COLDAT (CURVE) = ITEM

C              Convert fit type string to uppercase and decode.  Use local
C              flag NOFIT to remember if the original request was just for
C              the DEFAULT as opposed to LINEAR (to be used below to resolve
C              possible LINE conflicts).

               CALL UPCASE (FIT)
               NOFIT = FIT .EQ. CFLAG
               CALL LOOKUP (NFDIC, FITDIC, .FALSE., FIT, ITEM)
               IF (ITEM .LE. 0) THEN
                  WRITE (LWRITE, 1100) CURVE, FRAME
                  NOFIT = .TRUE.
                  FIT = CFLAG
                  CALL LOOKUP (NFDIC, FITDIC, .FALSE., FIT, ITEM)
               END IF
               FITDAT (CURVE) = ITEM

               IF (FIT .EQ. 'CLOSED') THEN
                  IF (NPT (CURVE) .GT. 2) THEN

C                    Protect the CLOSED option:  PLSFIT (called by QUICK)
C                    requires that the first and last points match.

                     CALL WRAPPER (MAXPT, NPT, CURVE, X, Y, IER)
                     IF (IER .NE. 0) THEN

C                       The data was not in wrap-around form and could not be
C                       patched, so gripe and revert to LINEAR interpolation.

                        WRITE (LWRITE, 1110) CURVE, FRAME
                        FITDAT (CURVE) = 1
                     END IF
                  ELSE ! IF (NPT (CURVE) .LE. 2) THEN

C                    Can't really "close" a curve with only two points! Just
C                    specify LINEAR and continue.

                     FITDAT (CURVE) = 1
                  END IF
               END IF

C              Convert spline type string to uppercase and decode.
C              Use local flag NOSPLINE to distinguish between an explicit
C              entry and the default (no entry).

               CALL UPCASE (SPLINE)
               NOSPLINE = SPLINE .EQ. CFLAG

               CALL LOOKUP (NSPLDIC, SPLDIC, .FALSE., SPLINE, ITEM)
               IF (ITEM .LE. 0) THEN
                  WRITE (LWRITE, 1105) CURVE, FRAME
                  SPLINE = CFLAG
                  CALL LOOKUP (NSPLDIC, SPLDIC, .FALSE., SPLINE, ITEM)
               ELSE IF (LINDAT (CURVE) .EQ. 1) THEN ! 1 = LINEAR
                  IF (.NOT. NOSPLINE) WRITE (LWRITE, 1107) CURVE, FRAME
               ELSE IF (LINDAT (CURVE) .EQ. 4) THEN ! 4 = CLOSED
                  IF (ITEM .LT. 3) THEN
                     ITEM = 3 ! Only parametric is appropriate for closed crvs.
                     WRITE (LWRITE, 1109) CURVE, FRAME
                  END IF
               END IF
               SPLDAT (CURVE) = ITEM

C              Convert line type string to uppercase and decode.  Note that
C              if a fit was explicitly requested above, then symbols-only
C              doesn't make sense here so we issue a message and override.

               CALL UPCASE (LINE)
               CALL LOOKUP (NLNDIC, LINDIC, .FALSE., LINE, ITEM)
               IF (ITEM .LE. 0) THEN
                  WRITE (LWRITE, 1120) CURVE, FRAME
                  LINE = CFLAG
                  CALL LOOKUP (NLNDIC, LINDIC, .FALSE., LINE, ITEM)
               ELSE IF (LINE .EQ. 'SYMBOLS' .AND. .NOT. NOFIT) THEN
                  WRITE (LWRITE, 1125) CURVE, FRAME
                  LINE = CFLAG
                  CALL LOOKUP (NLNDIC, LINDIC, .FALSE., LINE, ITEM)
               END IF
               LINDAT (CURVE) = ITEM

C              Save symbol type.  Negative values mean no symbol specifically
C              requested - the interpretation of this depends on the line
C              chosen.  For simple lines-connected-by-symbols and for symbols
C              alone, we cycle through the available types (DISSPLA wraps on
C              SYMBOL > 18), and increment the default only if it was used.
C              For the special line types, the default is NO symbol, in which
C              case the negative value is just passed along.

               IF (LINE .EQ. 'CONNECTED' .OR. LINE .EQ. 'SYMBOLS') THEN
                  IF (SYMBOL .LT. 0) THEN
                     SYMBOL = DEFSYM
                     DEFSYM = DEFSYM + 1
                  END IF
               END IF
               SYMDAT (CURVE) = SYMBOL

               IF (LEG .LT. MAXLEG .AND. LEGEND .NE. BLANK) THEN

C                 Pack the string into LEGTXT and keep track of length.  Assign 
C                 current linetype, color and symbol to parallel arrays. 

                  LEG = LEG + 1
                  FIRST = 1
                  LAST  = MAXLEN
                  CALL SCAN2 (LEGEND, SEPS, FIRST, LAST, MARK)
                  IF (LAST .NE. 0 .AND. 
     >               LEGEND (FIRST : LAST) .NE. BLANK) THEN
                     LEGTXT (LEG) (1 : LAST) = LEGEND
		     IF (STROKE) THEN
			LAST = MIN (LAST + 2, MAXLEN)
			LEGTXT (LEG) (LAST - 1 : LAST) = ESCAPE
                     END IF
                     LENLEG (LEG) = LAST
                  ELSE
                     LEGTXT (LEG) = BLANK
                     LENLEG (LEG) = 1
                  END IF
                    
                  LEGINF (1, LEG) = LINDAT (CURVE)
                  LEGINF (2, LEG) = SYMDAT (CURVE)
                  LEGINF (3, LEG) = COLDAT (CURVE)
               ELSE IF (LEGEND .NE. BLANK) THEN

C                 A significant legend string will have to be ignored.

                  WRITE (LWRITE, 1060) MAXLEG, CURVE
               END IF
            END IF

 
C           Entry of other curve-level controls may also require action.
C           ------------------------------------------------------------

            IF (FORMAT .LT. -2 .OR. FORMAT .GT. 2) THEN

C              FORMAT is out of range.  Default value is assumed. 

               FORMAT = 1
               WRITE (LWRITE, 1190)
            END IF

            IF (GRID .EQ. 0) THEN

C              Assign default value if GRID is specified as 0.

               GRID = -1
               WRITE (LWRITE, 1200)
            END IF

            IF (FORMAT .EQ. -1 .OR. FORMAT .EQ. -2) THEN
  
C              For axis plots, tick marks are to be drawn in spite of request
C              for grid lines.
    
               GRID = -ABS (GRID)
               WRITE (LWRITE, 1210)
            END IF

C           If annotations were entered, look for scaling their coordinates.

            DO 490, I = 0, 9

C              Allow blank text in case just an arrow is desired.

               IF (XNOTEQ (I) .NE. FLAG .AND. YNOTEQ (I) .NE. FLAG) THEN
                  NUMNOTE = NUMNOTE + 1
                  IF (NUMNOTE .LT. MAXNOTE) THEN
                     NOTES  (NUMNOTE) = NOTEQ  (I)
                     XNOTES (NUMNOTE) = XNOTEQ (I) * XSCALE + XSHIFT
                     YNOTES (NUMNOTE) = YNOTEQ (I) * YSCALE + YSHIFT
                     XARROWS(NUMNOTE) = XARROWQ(I)
                     IF (XARROWQ (I) .NE. FLAG) XARROWS (NUMNOTE) =
     >                   XARROWQ (I) * XSCALE + XSHIFT
                     YARROWS(NUMNOTE) = YARROWQ(I)
                     IF (YARROWQ (I) .NE. FLAG) YARROWS (NUMNOTE) =
     >                   YARROWQ (I) * YSCALE + YSHIFT
                  ELSE
                     WRITE (LWRITE, 1095) MAXNOTE
                  END IF
               END IF
  490       CONTINUE


C           Another curve?

            IF (.NOT. (FINIS .OR. ENDFRAME .OR. ERROR)) GO TO 400


C        Handle the frame-level inputs.
C        ------------------------------

         IF (CURVE .EQ. 0) THEN         ! Too cumbersome not to jump out early.
            WRITE (LWRITE, 1160) FRAME
            IF (ERROR) THEN             ! Quit because of error reading data.
               WRITE (LWRITE, 1170) CURVE, FRAME
               GO TO 980
            ELSE                        ! Keep going - may be more frames.
               GO TO 100
            END IF
         END IF

C        Convert lettering style to uppercase and decode.

         CALL UPCASE (STYLE)
         CALL LOOKUP (NSTDIC, STYDIC, .FALSE., STYLE, ITEM)
         IF (ITEM .LE. 0) THEN
            WRITE (LWRITE, 1140) FRAME
            STYLE = CFLAG
            CALL LOOKUP (NSTDIC, STYDIC, .FALSE., STYLE, ITEM)
         END IF

C        Convert orientation flag to uppercase, decode, and set default plot
C        height and width.

         CALL UPCASE (ORIENT)
         CALL LOOKUP (NORDIC, ORDIC, .TRUE., ORIENT, ITEM)
         IF (ITEM .LE. 0) THEN
            WRITE (LWRITE, 1130) FRAME
            ORIENT = CFLAG
            CALL LOOKUP (NORDIC, ORDIC, .TRUE., ORIENT, ITEM)
         END IF

C        Was integer axis-numbering requested?  Can it be the default?

         IF (XNUMBERS .EQ. CFLAG) THEN  ! It was defaulted
            IF (ALLXCOL0) THEN
               ITEM = 2
            ELSE
               ITEM = 1
            END IF
         ELSE  ! It wasn't defaulted, so ignore the possibility of all XCOLs = 0
            CALL UPCASE (XNUMBERS)
            CALL LOOKUP (NNUMDIC, NUMDIC, .FALSE., XNUMBERS, ITEM)
            IF (ITEM .LE. 0) THEN
               WRITE (LWRITE, 1145) FRAME
               ITEM = 1
            END IF
         END IF
         INTEGERX = ITEM .NE. 1

         IF (YNUMBERS .EQ. CFLAG) THEN
            IF (ALLYCOL0) THEN
               ITEM = 2
            ELSE
               ITEM = 1
            END IF
         ELSE
            CALL UPCASE (YNUMBERS)
            CALL LOOKUP (NNUMDIC, NUMDIC, .FALSE., YNUMBERS, ITEM)
            IF (ITEM .LE. 0) THEN
               WRITE (LWRITE, 1145) FRAME
               ITEM = 1
            END IF
         END IF
         INTEGERY = ITEM .NE. 1

C        Set unit number for DISSPLA messages.  Negative leaves only fatal
C        messages enabled.

         IF (DEBUG) THEN
            LDEBUG = +LWRITE
         ELSE
            LDEBUG = -LWRITE
         END IF

         TOTAL = 0
         DO 500, I = 1, CURVE
            TOTAL = TOTAL + NPT (I)
  500    CONTINUE


         IF (NUMTIT .EQ. 3) THEN

C           Special case:  no subtitle was present - null it out and shift
C           the axis labels "down" by one position.

            NUMTIT = 4
            TITLES (4) = TITLES (3)
            TITLES (3) = TITLES (2)
            TITLES (2) = BLANK
            LENTIT (4) = LENTIT (3)
            LENTIT (3) = LENTIT (2)
            LENTIT (2) = 0
         END IF

C        Insert file/date identifier as 2 extra caption lines?

         IF (IDENT) THEN
            IF (NUMTIT .LT. MAXTIT - 1) THEN

C              Make sure the identifier doesn't show up as a title or label.

               DO 510, I = NUMTIT + 1, 4
                  NUMTIT = I
                  TITLES (I) = BLANK
                  LENTIT (I) = 0
  510          CONTINUE

               NUMTIT = NUMTIT + 1
               TITLES (NUMTIT) = BLANK  ! Good idea even if no other captions.
               LENTIT (NUMTIT) = 0
               NUMTIT = NUMTIT + 1

               TITLES (NUMTIT) = 'File name: ' //
     >                           DATAFILE (FIRST2 : LAST3)
               LAST = LAST3 - FIRST2 + 18
               I = LAST + 9
               TITLES (NUMTIT) (LAST : I) = 'Plot date:'
               I = I + 1

C*****         CALL DATE (TITLES (NUMTIT) (I + 1 :))  ! Only 2-digit year

               CALL DATE_AND_TIME (CDATE_TIME (1), CDATE_TIME (2),
     >                             CDATE_TIME (3), IDATE_TIME)

C              DATE_AND_TIME returns CCYYMMDD, hhmmss.sss, +/-hhmm, & 8 integers

               LAST = I + 18 - MAXLEN  ! Avoid going off the end
               IF (LAST .GT. 0) I = I - LAST

               TITLES (NUMTIT) (I + 1 :) = 'CCYY-MM-DD @ hh:mm'
               TITLES (NUMTIT) (I + 1 : I + 4)  = CDATE_TIME (1) (1 : 4)
               TITLES (NUMTIT) (I + 6 : I + 7)  = CDATE_TIME (1) (5 : 6)
               TITLES (NUMTIT) (I + 9 : I + 10) = CDATE_TIME (1) (7 : 8)
               TITLES (NUMTIT) (I + 14: I + 15) = CDATE_TIME (2) (1 : 2)
               TITLES (NUMTIT) (I + 17: I + 18) = CDATE_TIME (2) (3 : 4)

               LAST = MAXLEN
               CALL SCAN2 (TITLES (NUMTIT), BLANK, FIRST, LAST, MARK)
               LENTIT (NUMTIT) = LAST
            ELSE
               WRITE (LWRITE, 1090) MAXTIT
            END IF
         END IF


C        If classification titles were entered, determine their length,
C        and move them up if necessary.
 
         DO 520, I = MAXTIT + 1, MAXTIT + 2
            IF (TITLES (I) .NE. BLANK) THEN
               FIRST = 1
               LAST  = MAXLEN
               CALL SCAN2 (TITLES (I), SEPS, FIRST, LAST, MARK)
               LENTIT (I) = LAST

               IF (NUMTIT .LT. MAXTIT) THEN
                  TITLES (NUMTIT + I - MAXTIT) = TITLES (I)
                  LENTIT (NUMTIT + I - MAXTIT) = LENTIT (I)
               END IF
            END IF
  520    CONTINUE


C        Appending the instruction alphabet escape characters is a kludge
C        for forcing DISSPLA to use stroke characters in each text string.

	 IF (STROKE) THEN
	    DO 530, I = 1, NUMTIT
               IF (LENTIT (I) .GT. 0) THEN
                  LENTIT (I) = MIN (LENTIT (I) + 2, MAXLEN)
	          TITLES (I) (LENTIT (I) - 1 : LENTIT (I)) = ESCAPE
               END IF
  530       CONTINUE
	 END IF


C        Determine output device(s) for this frame.
C        ------------------------------------------

C        If not sending all frames to metafile only (yet), then prompt.

         IF (.NOT. ALL) THEN

            OPT = 'P'
            IF (FRAME .NE. 1) THEN
               WRITE (PROMPT (38:39), '(I2)') FRAME
               WRITE (LSCREEN, 1020) PROMPT, PLOTMENU
               CALL READC (LSCREEN,
     >            'Choose option (default = P(review)): ',
     >            LKEYBD, OPT, CR, QUIT)
            END IF

            IF (OPT .EQ. 'M') THEN     ! All output to metafile.
               PREVIEW = .FALSE.
               IF (META) THEN
                  ALL = .TRUE.
               ELSE
                  WRITE (LSCREEN, 1020) ' Warning: No metafile is open.'
                  CALL PLOTDEV ('B', METAFILE (1:LENM), LSCREEN, LKEYBD,
     >               LMETA, 'Metafile selection: ', NDEV (2),
     >               DEVLIST (1, 2), TYPELIST (1, 2),
     >               HWPARAMS (1, 1, 2), DEVINDEX (2))
                  META = DEVINDEX (2) .NE. 0
                  ALL = META
               END IF

            ELSE IF (OPT .EQ. 'S') THEN     ! Skip to next frame.
               GO TO 100

            ELSE IF (OPT .EQ. 'Q') THEN     ! End session.
               GO TO 990
            END IF
         END IF


C        Begin plotting.
C        ---------------

         GRAPHICS = .TRUE.        ! Helps proper clean-up at the end.
         IF (PREVIEW) THEN        ! Not all to just metafile (yet).
            DEVINDX = 1
         ELSE
            DEVINDX = 2
         END IF

C        Start of loop over up to two output devices (screen and/or metafile):

         IF (ORIENT .EQ. 'PORTRAIT') THEN
            DEFHGT = HPORT
            DEFWID = WPORT
        ELSE           ! Landscape
            DEFHGT = HLAND
            DEFWID = WLAND
         END IF

         REPLOT = .FALSE.
  600    CONTINUE

C           Initialize output device.

            POINTER = DEVINDEX (DEVINDX)
            TYPE = TYPELIST (POINTER, DEVINDX)
            XSCLFR = HWPARAMS (POINTER, 1, DEVINDX)
            YSCLFR = HWPARAMS (POINTER, 2, DEVINDX)

            CALL PLOTDEV ('I', METAFILE (1:LENM), LSCREEN, LKEYBD,
     >         LMETA, BLANK, NDEV (DEVINDX), DEVLIST (1, DEVINDX),
     >         TYPELIST, HWPARAMS, POINTER)

            IF (.NOT. REPLOT) THEN

C              Determine unspecified plot frame parameters.

               CALL LAYOUT (LWRITE, TOTAL, X, Y, XMIN, XMAX, XSTEP,
     >            YMIN, YMAX, YSTEP, HEIGHT, WIDTH, DEFHGT, DEFWID,
     >            PLOT)

               HSAVE = HEIGHT
               WSAVE = WIDTH

            ELSE  ! No need to repeat LAYOUT
               HEIGHT = HSAVE
               WIDTH  = WSAVE
            END IF

C           Apply calibration factors for drawing axes in exact inches 
C           on the current output device (from qplot.config file).

            IF (ORIENT .EQ. 'PORTRAIT') THEN
               HEIGHT = HEIGHT * YSCLFR 
               WIDTH  = WIDTH  * XSCLFR
            ELSE
               HEIGHT = HEIGHT * XSCLFR
               WIDTH  = WIDTH  * YSCLFR
            END IF

C           ---------------------------------------
C           Generate a DISSPLA plot for this frame.
C           ---------------------------------------

            CALL QUICK
     >        (FRAME, NUMTIT, TITLES, LENTIT, MAXPT, CURVE, NPT, X, Y,
     >         XMIN, XMAX, XSTEP, YMIN, YMAX, YSTEP, Y2LABEL, Y2SCALE,
     >         Y2SHIFT, HEIGHT, WIDTH, ORIENT, PLOT, STYLE, FORMAT,
     >         GRID, FITDAT, SPLDAT, LINDAT, SYMDAT, COLDAT, LEG,LEGTXT,
     >         LENLEG, LEGINF, NUMNOTE, NOTEBOX, NOTES, XNOTES, YNOTES,
     >         XARROWS, YARROWS, INTEGERX, INTEGERY, ESCAPE, LDEBUG)

            WRITE (LWRITE, 1150) CURVE, FRAME

            IF (ERROR) THEN             ! Quit because of error reading data.
               WRITE (LWRITE, 1170) CURVE, FRAME
               GO TO 980
            END IF

            IF (DEVINDX .EQ. 1) THEN

C              Clear plot from screen and home the cursor in VT100 mode.

               IF (TYPE .EQ. 'TEK') THEN
                  CALL SCLEAN (3)
                  CALL SCLEAN (1)
!!!            ELSE IF (TYPE .EQ. 'VT') THEN   ! What?
!!!            ELSE IF (TYPE .EQ. 'X' ) THEN   ! What?
               END IF

               IF (META) THEN
                  YES = .TRUE.
                  CALL READY (LSCREEN, 'Do you wish to add this plot' //
     >               ' to the metafile? (<CR>=Y) ',
     >               LKEYBD, YES, CR, QUIT)
                  IF (QUIT) GO TO 980
                  IF (YES) THEN
                     DEVINDX = 2      ! Make the device the metafile.
                     REPLOT = .TRUE.  ! Avoid redoing the layout.
                     GO TO 600
                  END IF
               END IF
            END IF

C        Loop back for another plot frame.
C        ---------------------------------

         GO TO 100


C     Abnormal termination.
C     ---------------------

  980 CONTINUE
      ENDSTATE = 3
      GO TO 990

C     Very early exit (nothing in qplot.out).
C     ---------------------------------------

  985 CONTINUE
      ENDSTATE = 1
      GO TO 990

C     Normal termination.   ENDSTATE = 0 from initialization.
C     -------------------------------------------------------

  990 CONTINUE
      IF (GRAPHICS) THEN
         CALL DONEPL     ! Sign off the plotting device.

         IF (SCREEN) THEN   ! Clean up and home the cursor in VT100 mode.
            TYPE = DEVTYPE (1)
            IF (TYPE .EQ. 'TEK') THEN
               CALL SCLEAN (3)
               CALL SCLEAN (1)
!!!         ELSE IF (TYPE .EQ. 'VT') THEN   ! What?
!!!         ELSE IF (TYPE .EQ. 'X' ) THEN   ! What?
            END IF
         END IF
      END IF

      IF (ENDSTATE .EQ. 0) THEN         ! Normal end
         WRITE (LWRITE, 1180) FRAME - 1
         CALL EXIT (1)
      ELSE IF (ENDSTATE .EQ. 1) THEN    ! Early exit - nothing to do?
      ELSE IF (ENDSTATE .EQ. 3) THEN    ! Point user to diagnostics file.
         WRITE (LSCREEN, 1220)
         CALL EXIT (3)
      END IF

C*****STOP ' '  ! Avoid machine dependencies - final END suffices.

C     Formats.
C     --------

 1000 FORMAT
     >  (' -------------- QPLOT (July 2002 version) -------------'/
     >   '                  run-time diagnostics'//
     >   ' Use $OPTIONS DEBUG = .TRUE. $END for DISSPLA messages.'//)
 1010 FORMAT (/' QPLOT:  Aborting with no numeric data found in frame ',
     >   I3)
 1015 FORMAT (/' QPLOT:  System read error - aborting.  IOS:', I5)
 1020 FORMAT (A)
 1030 FORMAT (10X, '>>', A60)
 1040 FORMAT (10X, A60)
 1050 FORMAT (/' QPLOT:   Abnormal termination due to error while ',
     >   'reading $OPTIONS for curve '/
     >   10X, 'number ', I3, ' in frame number ', I2, '.  ',
     >   'No curve(s) were plotted in this frame.'//)
 1060 FORMAT (/' QPLOT:   Warning - only the first ', I2, ' curves may',
     >   ' have a legend.'/
     >   10X, 'The legend for curve number ', I3, ' will be ignored.')
 1070 FORMAT (/' QPLOT:   Warning - color was improperly ',
     >   'specified for curve '/
     >   10X, 'number ', I3, ' in frame number ', I3, '.  ',
     >   'The default will be used.')
 1080 FORMAT (/' QPLOT:   Trouble - illegal input data beginning '/
     >   10X, 'with:')
 1090 FORMAT (/' QPLOT:   Warning - caption line limit exceeded.  The ',
     >   'limit is ', I2/
     >   10X, 'including file/date id. if specified.  Excess captions ',
     >   'are ignored.')
 1095 FORMAT (/' QPLOT:   Warning - annotation/arrow limit exceeded.  ',
     >   'The limit is ', I2/
     >   10X, 'Excess annotations are ignored.')
 1100 FORMAT (/' QPLOT:   Warning - fit was improperly ',
     >   'specified for curve '/
     >   10X, 'number ', I3, ' in frame number ', I3, '.  ',
     >   'The default will be used.')
 1105 FORMAT (/' QPLOT:   Warning - spline type was improperly ',
     >   'specified for curve '/
     >   10X, 'number ', I3, ' in frame number ', I3, '.  ',
     >   'The default will be used.')
 1107 FORMAT (/' QPLOT:   Warning - piecewise linear fit ',
     >   'indicated for curve '/
     >   10X, 'number ', I3, ' in frame number ', I3, '.  ',
     >   'Specified spline fit will be suppressed.')
 1109 FORMAT (/' QPLOT:   Warning - closed fit type ',
     >   'specified for curve '/
     >   10X, 'number ', I3, ' in frame number ', I3, '.  ',
     >   'Switching to parametric spline.')
 1110 FORMAT (/' QPLOT:   Warning - insufficient storage to close ',
     >   'the data for curve '/
     >   10X, 'number ', I3, ' in frame number ', I3, '.  ',
     >   'Linear interpolation will be used.')
 1120 FORMAT (/' QPLOT:   Warning - line type was improperly ',
     >   'specified for curve '/
     >   10X, 'number ', I3, ' in frame number ', I3, '.  ',
     >   'The default will be used.')
 1125 FORMAT (/' QPLOT:   Warning - LINE = ''SYMBOLS'' (only) ',
     >   'conflicts with curve fit specified '/
     >   10X, 'for curve number ', I3, ' in frame number ', I3, '.'/,
     >   10X, 'The curve fit will be shown with symbols.')
 1130 FORMAT (/' QPLOT:   Warning - plot orientation was improperly ',
     >   'specified for '/
     >   10X, 'frame number ', I3, '.  The default will be used.')
 1140 FORMAT (/' QPLOT:   Warning - lettering style was improperly ',
     >   'specified for '/
     >   10X, 'frame number ', I3, '.  The default will be used.')
 1145 FORMAT (/' QPLOT:   Warning - axis numbering was improperly ',
     >   'specified for '/
     >   10X, 'frame number ', I3, '.  The default will be used.')
 1150 FORMAT (/' QPLOT:   ', I3, ' curve(s) plotted in frame number ',
     >   I3, '.'//)
 1160 FORMAT (/' QPLOT:   Warning - frame number ', I3, ' is empty.'//)
 1170 FORMAT (/' QPLOT:   Abnormal termination - the last frame ',
     >   'plotted may '/
     >   10X, 'be incomplete; ', I3, ' curve(s) plotted in frame ',
     >   'number ', I3, '.'//)
 1180 FORMAT (/' QPLOT:   Normal termination. ', I3, ' frame(s) ',
     >   'were plotted.'//)
 1190 FORMAT (/' QPLOT:   Warning - FORMAT is out of range.  The ',
     >   'default, FORMAT = 1, '
     > / 10X, 'is assumed.'/)
 1200 FORMAT (/' QPLOT:   Warning - GRID = 0 is not valid.  Changed ',
     >   'to default, GRID = -1.')
 1210 FORMAT (/' QPLOT:   Warning - grid lines requested where tick ',
     >   'marks appropriate.  '
     > / 10X, 'Tick marks will be used.')
 1220 FORMAT (/' Diagnostics are in the file qplot.out.'//)

      END
C+----------------------------------------------------------------------
C
      SUBROUTINE ANNOTATE (NUMNOTE, NOTES, XNOTES, YNOTES, XARROWS,
     >                     YARROWS, NOTEBOX, FLAG, LUNERR)
C
C ONE-LINER:  Annotate a plot with or without arrows (DISSPLA level 2, 3)
C
C PURPOSE:
C
C        ANNOTATE adds NUMNOTE text strings to a plot, with the lower
C     left point of string I defined by (XNOTES (I), YNOTES (I)) in
C     data units.  If (XARROWS (I), YARROWS (I)) is also defined
C     (= (FLAG, FLAG) if undefined), an arrow is drawn from the CENTER
C     of the text string rectangle.  Each arrow is clipped by a border
C     of width .05" around the annotation.  The area within this border
C     may be outlined and/or blanked (to prevent grid lines, etc. from
C     intruding).
C
C METHOD:
C
C        DISSPLA's MESSAG draws text at a point specified in inches,
C     and we need to convert from data units if an arrow is involved,
C     so choose MESSAG over RLMESS, which could draw the text too.
C     Clipping of the optional arrows is done explicitly in case no
C     blanking is requested.  Of many possible arrow types, the one
C     defined by the integer 1001 in DISSPLA's terminology for its
C     VECTOR utility is the only one provided here (i.e. single arrow
C     head, shaded, with modestly sharp point).  DISSPLA's BLREC blanks
C     each rectangle once the text is drawn, if blanking is requested.
C
C        The current character height (retrieved via WHATIS for the
C     rectangle calculations) applies to both the height of the text
C     and the length of the arrow head.
C
C ERROR HANDLING:
C
C        If either coordinate for an annotation is undefined (= FLAG),
C     an error message is written to LUNERR, and the annotation is skipped.
C
C ENVIRONMENT:
C     VAX/VMS, FORTRAN 77, with:
C     > IMPLICIT NONE
C     > Trailing ! comments
C     > Names up to 8 characters
C
C HISTORY:
C     09/27/91  D.A.Saunders  Initial design and code, for QPLOT.
C     10/19/91    "     "     Achieving precisely horizontal or vertical
C                             arrows was awkward originally.  Now,
C                             YARROWS (I) = YNOTES (I) is taken to mean
C                             "horizontal" (whether the text is blank or
C                             not).  Likewise, XARROWS (I) = XNOTES (I)
C                             means "vertical" now, from center of text.
C
C AUTHOR: David Saunders, Sterling Software/NASA Ames, Mountain View, CA.
C
C ----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments.

      INTEGER   NUMNOTE         ! (I) Number of annotations.  NUMNOTE >= 0.
      CHARACTER NOTES (*) * (*) ! (I) Annotation text strings.  If blank,
                                !        an arrow may still be drawn between
                                !        the specified points.
      REAL      XNOTES (*)      ! (I) Coordinates of bottom left corner of
      REAL      YNOTES (*)      !        annotation text, in data units.
      REAL      XARROWS (*)     ! (I) Coordinates of tips of arrow heads,
      REAL      YARROWS (*)     !        in data units.  An arrow is suppressed
                                !        if either coordinate equals FLAG.
                                !        To ensure a HORIZONTAL arrow, set
                                !        YARROWS (I) = YNOTES (I) (whether
                                !        NOTES (I) is blank or not).  Likewise,
                                !        set XARROWS (I) = XNOTES (I) to achieve
                                !        a precisely VERTICAL arrow from the
                                !        center of text.
      INTEGER   NOTEBOX         ! (I) Controls boxes around annotations:
                                !        0 means no frame and no blanking;
                                !        1 means no frame/rectangle is blanked;
                                !        2 means box is framed and blanked;
                                !       >2 means frame is thickened outwards:
                                !          (NOTEBOX - 1) scales the pen width.
      REAL      FLAG            ! (I) Flag for undefined coordinates. E.g. 999.
      INTEGER   LUNERR          ! (I) Logical unit number for error messages.

C     Procedures.

      REAL      XMESS
      EXTERNAL  BLREC           ! (DISSPLA) Blanks a non-tilted rectangle.
      EXTERNAL  MESSAG          ! (DISSPLA) Draws text at (X, Y) inches.
      EXTERNAL  SCAN2           ! Determines lengths of annotation strings.
      EXTERNAL  TOINCH          ! (DISPPLA) Converts from data units to inches.
      EXTERNAL  VECTOR          ! (DISSPLA) Draws an arrow.
      EXTERNAL  WHATIS          ! (DISSPLA) Used to determine HEIGHT setting.
      EXTERNAL  XMESS           ! (DISSPLA) Length of text string in inches.

C-----------------------------------------------------------------------

C     Local constants.

      INTEGER   IVEC
      REAL      BORDER, EPS, HALF
      CHARACTER BLANK * 1
      PARAMETER
     >  (BORDER = 0.05,         ! Width of borders around text (inches)
     >   EPS    = .01,          ! Need to avoid TAN (90 deg.)
     >   HALF   = 0.5E+0,
     >   IVEC   = 1001,         ! See IVEC arrow definition in Ch. 5, Sec. 2.1
     >   BLANK  = ' ')

C     Local variables.

      INTEGER
     >   FIRST, I, LAST, LENNOTE, MARK
      REAL
     >   FRM, HIGH, PHI, PI, TANTH, THETA, XA, XC, XE, XL, XR,
     >   YA, YB, YC, YE, YT
      LOGICAL
     >   HORIZ, VERTI, VERTICAL

C     Execution.

      LENNOTE = LEN (NOTES (1))
      PI = 4.0 * ATAN2 (1.0, 1.0)
      FRM = REAL (1 - NOTEBOX)   ! 0 means no frame, but text area is blanked.
                                 ! -ve causes frame thickening toward OUTside.

      DO 500, I = 1, NUMNOTE

         XL = XNOTES (I)
         YB = YNOTES (I)

         IF (XL .EQ. FLAG .OR. YB .EQ. FLAG) THEN

            WRITE (LUNERR, 1010) I

         ELSE

            XA = XARROWS (I)
            YA = YARROWS (I)
            HORIZ = YA .EQ. YB
            VERTI = XA .EQ. XL

            CALL TOINCH (XL, YB, XL, YB) ! Convert lower left coords. of string
                                         ! to inches, as needed for box & arrow

C           Draw the annotation unless it is blank.

            FIRST = 1
            LAST = LENNOTE
            CALL SCAN2 (NOTES (I), BLANK, FIRST, LAST, MARK)

            IF (LAST .GT. 0) THEN
               CALL MESSAG (NOTES (I) (1 : LAST), LAST, XL, YB)
            END IF

C           Determine text box border for blanking purposes, even if
C           no arrow is drawn.  XL, YB are modified for non-blank text, so:

            XE = XL                 ! End of arrow, if arrow is
            YE = YB                 ! requested but string is blank

            XR = XL + XMESS (NOTES (I), LAST) + BORDER
            XL = XL - BORDER

            CALL WHATIS ('HEIGHT', BLANK, MARK, HIGH, MARK)

            YT = YB + HIGH + BORDER
            YB = YB - BORDER

C           Draw any arrow (explicitly clipped) before any blanking of the text.

            IF (XA .NE. FLAG .AND. YA .NE. FLAG) THEN

               CALL TOINCH (XA, YA, XA, YA)

               IF (LAST .GT. 0) THEN

C                 The arrow (suitably clipped) will emanate from the center
C                 (XC, YC) of the rectangle (with border) containing the string.

                  XC = (XL + XR) * HALF
                  YC = (YB + YT) * HALF

                  IF (HORIZ) YA = YC               ! Fudge - otherwise the user
                  IF (VERTI) XA = XC               ! has to mix data units and
                                                   ! plot inches to estimate the
                                                   ! center of the text box.

                  PHI = ATAN2 (YT - YC, XR - XC)   ! Angle of line from center
                                                   ! to top right of box in
                                                   ! the range (0, PI/2 - del) 

                  THETA = ATAN2 (YA - YC, XA - XC) ! Angle of arrow in the
                                                   ! range [-PI, PI].

                  VERTICAL = ABS (ABS (THETA) - PI * HALF) .LT. EPS
                  IF (VERTICAL) THEN
                     XE = XC
                  ELSE
                     TANTH = (YA - YC) / (XA - XC)
                  END IF

C                 Determine which side of the rectangle clips the arrow,
C                 and adjust the end of the arrow (XE, YE) accordingly.

                  IF (ABS (THETA) .LE. PHI) THEN             ! Right side clips.
                     XE = XR
                     YE = YC + (XR - XC) * TANTH
                  ELSE IF (THETA .GT. PHI .AND.
     >                     THETA .LE. PI - PHI) THEN         ! Top clips.
                     IF (.NOT. VERTICAL) XE = XC + (YT - YC) / TANTH
                     YE = YT
                  ELSE IF (ABS (PI - THETA) .LE. PHI) THEN   ! Left side clips.
                     XE = XL
                     YE = YC + (XL - XC) * TANTH
                  ELSE                                       ! Bottom clips.
                     IF (.NOT. VERTICAL) XE = XC + (YB - YC) / TANTH
                     YE = YB
                  END IF
               END IF

C              Draw the arrow unless it is absurdly short.

               IF ((XA - XE) ** 2 + (YA - YE) ** 2 .GT. HIGH ** 2)
     >            CALL VECTOR (XE, YE, XA, YA, IVEC)

            END IF

C           Blank and/or outline the text area?

            IF (FRM .LE. 0.0) CALL BLREC (XL, YB, XR - XL, YT - YB, FRM)

         END IF

  500 CONTINUE

      RETURN

C     Formats.

 1010 FORMAT ('0*** ANNOTATE: Undefined coordinate(s) for annotation #',
     >        I3, '.  Proceeding...')

      END
C+----------------------------------------------------------------------
C
      FUNCTION ALPHA (STRING)
C
C
C     Description and usage:
C
C           A simple(-minded) test for numeric data is implemented by
C        searching an input string for disqualifying characters.  Most,
C        but not all, disqualifiers are included in the definition of
C        statement function LETTER - included are some punctuation marks
C        and the alphabet, except for the letters D and E which may be
C        part of a single or double precision quantity.  Some insurance
C        is provided by requiring that any numeric string have at least
C        one digit.  Note that a few ambiguities remain:
C
C           (a)  A string might have the form of numerical data but be
C                intended as text - no general test can hope to detect
C                such cases.
C        
C           (b)  This routine does not check for correctness of the
C                data format.  A meaningless string such as 'E1E2E3'
C                will not be tagged as ALPHA, but is not numeric either.
C
C        Despite these weaknesses, this method should work in the vast
C        majority of ordinary cases.
C
C           This routine was written for use by QPLOT but may find other
C        applications where loosely-formatted input data is required.
C
C
C     Parameters:
C
C        Name    Dimension  Type  I/O/S  Description
C        ALPHA               L      O    Set .TRUE. if STRING is not
C                                        numerical data.
C        STRING              C    I      Input data to be tested.
C
C
C     Environment:  Digital VAX-11/780 VMS FORTRAN (FORTRAN 77).
C
C
C     Notes:
C
C        (1)  IMPLICIT NONE is non-standard.
C
C        (2)  The test and conversion of lowercase letters is meaningless
C             for systems which do not recognize lowercase - it may be
C             deleted.
C
C        (3)  The ASCII character set and collating sequence is assumed in
C             the choice of some symbols used and in the use of comparison
C             functions .LE. and .GE. (rather than the more general lexical
C             comparison library functions).
C
C        (4)  For QPLOT, any title line may be guaranteed to be alphabetic
C             by enclosing it in matching quotes (single or double).  The
C             quotes are stripped out by STRIPPER.
C
C        (5)  COMPLEX data with parentheses will look alphabetic.
C
C
C     Development history:
C
C        23 Feb. 1984      RAK     Initial design and coding.
C        27 Feb. 1984      RAK     Added more punctuation, including '$',
C                                  to the list of disqualifiers.
C        14 June 1984      RAK     Amended header.
C        16 Nov. 1991  D.Saunders  Lone dates like 11/16/91 are now deemed
C                                  alphabetic.  (Added '/' to LETTER;
C                                  also added ATOM .GE. '{'.)
C
C     Author:  Robert Kennelly, Informatics General Corporation.
C
C-----------------------------------------------------------------------


C     Declarations.
C     -------------

      IMPLICIT NONE

C     Arguments.

      LOGICAL
     >   ALPHA
      CHARACTER
     >   STRING * (*)

C     Local variables.

      LOGICAL
     >   DIGIT, LETTER, MAYBE
      INTEGER
     >   I, LENGTH
      CHARACTER
     >   ATOM * 1

C     Statement functions.

      LETTER (ATOM) = ATOM .GE. ':' .AND. ATOM .LE. '_' .AND.
     >                ATOM .NE. 'D' .AND. ATOM .NE. 'E' .OR.
     >                ATOM .GE. '!' .AND. ATOM .LE. '*' .OR.
     >                ATOM .GE. '{'  .OR. ATOM .EQ. '/'

      DIGIT (ATOM) =  ATOM .GE. '0' .AND. ATOM .LE. '9'


C     Execution.
C     ----------

      LENGTH = LEN (STRING)
      MAYBE = .FALSE.

      I = 0
   10 CONTINUE
         I = I + 1
         ATOM = STRING (I : I)

C        Convert any lowercase letters to uppercase.

         IF (ATOM .GE. 'a' .AND. ATOM .LE. 'z')
     >      ATOM = CHAR (ICHAR (ATOM) - 32)

C        The presence of any "letter" means the string is NOT numeric.  For
C        insurance, keep track of whether any honest digits have been found.

         ALPHA = LETTER (ATOM)
         MAYBE = MAYBE .OR. DIGIT (ATOM)
         IF (.NOT.ALPHA .AND. I .LT. LENGTH) GO TO 10

      ALPHA = ALPHA .OR. .NOT.MAYBE


C     Termination.
C     ------------

      RETURN
      END
C+----------------------------------------------------------------------
C
      FUNCTION BESSEL (J, H, DEL)
C
C     One-liner: First derivative using central 3-point formula
C     ----------
C
C     Description and usage:
C     ----------------------
C
C        Computes a first derivative approximation using the central
C     3-point formula.  The data must be in the form of arrays containing
C     finite difference interval lengths and 2-point forward difference
C     derivatives.  BESSEL is intended to be used by PLSFIT for determin-
C     ing end conditions on an interval for (non-monotonic) interpolation
C     by piecewise cubics.  See the PLSFIT header for more details.
C
C     Arguments:
C     ----------
C
C     Name    Type/Dimension  I/O/S  Description
C     J       I               I      Indicates at which end of the
C                                    interval the derivative is to be
C                                    estimated. J = 0 means left-hand
C                                    side, J = 1 means right.
C
C     H       R (-1:1)        I      Array of interval lengths. The 0th
C                                    element is the length of the interval
C                                    on which the cubic is to be deter-
C                                    mined.
C
C     DEL     R (-1:1)        I      Array of derivative estimates. The
C                                    0th element is the forward difference
C                                    derivative over the interval on which
C                                    the cubic is to be determined.
C                                     
C     BESSEL  R                 O    The function value is the adjusted
C                                    derivative.
C
C     Notes:
C     ------
C
C     (1)  IMPLICIT NONE is non-standard.
C
C     Author:  Robert Kennelly, Sterling Federal Systems/NASA-Ames
C     -------
C
C     Development history:
C     --------------------
C
C     18 Feb. 1987    RAK    Initial design and coding.
C
C-----------------------------------------------------------------------

C     Declarations.
C     -------------

      IMPLICIT NONE

C     Constants.

      REAL
     &   ONE
      PARAMETER
     &  (ONE = 1.0E+0)

C     Arguments.

      INTEGER
     &   J
      REAL
     &   H (-1:1), DEL (-1:1), BESSEL

C     Local variables.

      REAL
     &   WEIGHT

C     Execution.
C     ----------

C     Estimate first derivative on left (J = 0) or right side (J = 1) of
C     an interval.

      WEIGHT = H (J) / (H (J) + H (J - 1))
      BESSEL = WEIGHT * DEL (J - 1) + (ONE - WEIGHT) * DEL (J)

C     Termination.
C     ------------

      RETURN
      END
C+----------------------------------------------------------------------
C
      FUNCTION BRODLIE (J, H, DEL)
C
C     One-liner: First derivative, adjusted for monotonicity
C     ----------
C
C     Description and usage:
C     ----------------------
C
C        BRODLIE is intended to be used by PLSFIT for determining end
C     conditions on an interval for monotonic interpolation by piecewise
C     cubics. The data must be in the form of arrays containing finite
C     difference interval lengths and 2-point forward difference deriva-
C     tives. See the PLSFIT header for more details.
C
C        The method is due to Brodlie, Butland, Carlson, and Fritsch,
C     as referenced in the PLSFIT header.
C
C     Arguments:
C     ----------
C
C     Name    Type/Dimension  I/O/S  Description
C     J       I               I      Indicates at which end of the
C                                    interval the derivative is to be
C                                    estimated. J = 0 means left-hand
C                                    side, J = 1 means right.
C
C     H       R (-1:1)        I      Array of interval lengths. The 0th
C                                    element is the length of the interval
C                                    on which the cubic is to be deter-
C                                    mined.
C
C     DEL     R (-1:1)        I      Array of derivative estimates. The
C                                    0th element is the forward difference
C                                    derivative over the interval on which
C                                    the cubic is to be determined.
C                                     
C     BRODLIE R                 O    The function value is the adjusted
C                                    derivative.
C
C     Notes:
C     ------
C
C     (1)  IMPLICIT NONE and 8-character symbolic names are non-standard.
C
C     Author:  Robert Kennelly, Sterling Federal Systems/NASA-Ames
C     -------
C
C     Development history:
C     --------------------
C
C     18 Feb. 1987    RAK    Initial design and coding.
C
C-----------------------------------------------------------------------

C     Declarations.
C     -------------

      IMPLICIT NONE

C     Constants.

      REAL
     &   ZERO, ONE, THIRD
      PARAMETER
     &  (ZERO   = 0.0E+0,
     &   ONE    = 1.0E+0,
     &   THIRD  = ONE / 3.0E+0)

C     Arguments.

      INTEGER
     &   J
      REAL
     &   BRODLIE, H (-1:1), DEL (-1:1)

C     Local variables.

      REAL
     &   ALPHA

C     Execution.
C     ----------

C     Compare the algebraic signs of the two DEL's.  Have to test that
C     at least one is positive to avoid a zero denominator (this fancy
C     test permits one term to be zero, but the answer below is zero
C     anyway in these cases).  The trick is to work around the SIGN
C     function, which returns positive even if its 2nd argument is zero.

      IF (SIGN (ONE, -DEL (J - 1)) .NE. SIGN (ONE, DEL (J))) THEN

C        Form "weighted harmonic mean" of the two finite-difference
C        derivative approximations.  Note that we try to avoid overflow
C        by not multiplying them together directly.

         ALPHA   = THIRD * (ONE + H (J) / (H (J - 1) + H (J)))
         BRODLIE = DEL (J - 1) * (DEL (J) / (ALPHA * DEL (J) +
     &      (ONE - ALPHA) * DEL (J - 1)))
      ELSE

C        The signs differ, so make this point a local extremum.

         BRODLIE = ZERO
      END IF

C     Termination.
C     ------------

      RETURN
      END
C+----------------------------------------------------------------------
C
      FUNCTION BUTLAND (J, H, DEL)
C
C     One-liner: First derivative, non-central 3-point formula, adjusted
C     ----------
C
C     Description and usage:
C     ----------------------
C
C        Computes a first derivative approximation for PLSFIT over an
C     interval at a data boundary, using a modified forward or backward
C     3-point formula.  The data must be in the form of arrays containing
C     finite difference interval lengths and 2-point forward difference
C     derivatives, and the differencing direction is controlled by a flag.
C     See PLSFIT for more details, or THREEPT for the pure 3-pt. formula.
C
C        The "shape preserving adjustments" are from PCHIP, a monotone
C     piecewise cubic interpolation package by F. N. Fritsch.
C
C     Arguments:
C     ----------
C
C     Name    Type/Dimension  I/O/S  Description
C     J       I               I      Indicates at which end of the
C                                    interval the derivative is to be
C                                    estimated. J = 0 means left-hand
C                                    side, J = 1 means right. 
C
C     H       R (-1:1)        I      Array of interval lengths. The 0th
C                                    element is the length of the interval
C                                    on which the cubic is to be deter-
C                                    mined.
C
C     DEL     R (-1:1)        I      Array of derivative estimates. The
C                                    0th element is the forward difference
C                                    derivative over the interval on which
C                                    the cubic is to be determined.
C                                     
C     BUTLAND R                 O    The function value is the adjusted
C                                    derivative.
C
C     Environment:  VAX/VMS; FORTRAN 77
C     ------------
C
C     IMPLICIT NONE and 8-character symbolic names are non-standard.
C
C     Author:  Robert Kennelly, Sterling Federal Systems/NASA-Ames
C     -------
C
C     History:
C     --------
C
C     18 Feb. 1987    RAK    Initial design and coding, as THREEPT.
C     20 June 1991    DAS    Monotonic form renamed BUTLAND; THREEPT
C                            is now the pure 3-point formula.
C
C-----------------------------------------------------------------------

C     Declarations.
C     -------------

      IMPLICIT NONE

C     Arguments.

      INTEGER
     &   J
      REAL
     &   H (-1:1), DEL (-1:1), BUTLAND

C     Local constants.

      REAL
     &   ZERO, ONE, THREE
      PARAMETER
     &  (ZERO  = 0.0E+0,
     &   ONE   = 1.0E+0,
     &   THREE = 3.0E+0)

C     Local variables.

      INTEGER
     &   STEP
      REAL
     &   DMAX, WEIGHT

C     Execution.
C     ----------

C     Estimate first derivative on a left-hand boundary using a 3-point
C     forward difference (STEP = +1), or with a backward difference for
C     the right-hand boundary (STEP = -1).

      STEP = 1 - J - J   ! J here is consistent with related modules.

C     In {H, DEL} form, the derivative looks like a weighted average.

      WEIGHT  = -H (0) / (H (0) + H (STEP))
      BUTLAND = WEIGHT * DEL (STEP) + (ONE - WEIGHT) * DEL (0)

C     Shape-preserving adjustments.  Note that we try to avoid overflow
C     by not multiplying quantities directly.

      IF (SIGN (ONE, BUTLAND) .NE. SIGN (ONE, DEL (0)) .OR.
     &   DEL (0) .EQ. ZERO) THEN

C        Defer to the estimate closest to the boundary.

         BUTLAND = ZERO
      ELSE IF (SIGN (ONE, DEL (0)) .NE. SIGN (ONE, DEL (STEP))) THEN

C        If the monotonicity switches, may need to bound the estimate.

         DMAX = THREE * DEL (0)
         IF (ABS (BUTLAND) .GT. ABS (DMAX)) BUTLAND = DMAX
      END IF

C     Termination.
C     ------------

      RETURN
      END
C+----------------------------------------------------------------------
C
      FUNCTION CHORD (X, Y, I1, I2)
C
C     One-liner: Summed chord-lengths for X-Y curve over range of indices
C     ----------
C
C     Description and usage:
C     ----------------------
C
C        Computes the sum of the Euclidean distance between adjacent
C     points in a curve represented by two arrays. The calculation is
C     rearranged so as to avoid (almost) all chance of overflow, and
C     any unnecessary loss of precision when one component of the
C     distance from one point to the next is small relative to the other.
C     The calling routine must supply beginning and ending indices for
C     the summation (this is intended to facilitate operations with
C     packed data arrays). The result does not depend on the order of
C     I1 and I2.
C
C        CHORD was originally written for use with PLSFIT, which performs
C     parametric cubic interpolation with cumulative chord length as the
C     curve parameter. In use, it is a good idea to try to use all
C     available information to avoid (expensively) re-calculating the
C     lengths of the same intervals over and over; CHORD should be
C     thought of as providing length increments.
C
C     Arguments:
C     ----------
C
C     Name    Type/Dimension  I/O/S  Description
C     X         R (*)         I      Array of abscissas.
C
C     Y         R (*)         I      Array of ordinates.
C
C     I1,I2                   I      Indices for summation. The loop
C                                    runs from MIN(I1,I2) to MAX(I1,I2)
C                                    so the result is independent of
C                                    order.
C
C     CHORD   R                 O    Function value is the sum of the
C                                    chord lengths along the curve
C                                    defined by the X and Y arrays
C                                    between indices I1 and I2.
C
C     Notes:
C     ------
C
C     (1)  IMPLICIT NONE is non-standard.
C
C     Author:  Robert Kennelly, Sterling Federal Systems/NASA-Ames
C     -------
C
C     Development history:
C     --------------------
C
C     18 Feb. 1987    RAK    Initial design and coding.
C      8 Apr. 1988  RAK/DAS  Reformulated to reduce chance of overflow
C                            or unnecessary underflow. Result is not
C                            dependent on order of I1, I2.
C
C-----------------------------------------------------------------------

C     Declarations.
C     -------------

      IMPLICIT NONE

C     Constants.

      REAL
     &   ZERO, ONE
      PARAMETER
     &  (ZERO = 0.0E+0,
     &   ONE  = 1.0E+0)

C     Arguments.

      INTEGER
     &   I1, I2
      REAL
     &   CHORD, X (*), Y (*)

C     Local variables.

      INTEGER
     &   I
      REAL
     &   DENOM, DX, DY

C     Execution.
C     ----------

      CHORD = ZERO

      DO 10, I = MIN (I1, I2), MAX (I1, I2) - 1
         DX    = ABS (X (I + 1) - X (I))
         DY    = ABS (Y (I + 1) - Y (I))
         DENOM = MAX (DX, DY)
         IF (DENOM .GT. ZERO) CHORD = CHORD +
     &      DENOM * SQRT (ONE + (MIN (DX, DY) / DENOM) ** 2)
   10 CONTINUE

C     Termination.
C     ------------

      RETURN
      END

C+------------------------------------------------------------------------------
C
      SUBROUTINE CLASSIF (ORIENT, CLASS, EXPLAN, LENTXT)
C
C     One liner:  Writes classification headings on plots using DISSPLA.
C
C     Purpose:    CLASSIF is intended to write special headings for plots
C                 in either landscape or portrait mode.  The format consists
C                 of both a classification status in large print and an
C                 explanation, both to appear at the top and bottom of
C                 the page.  The headings are suitable for 8.5 x 11" paper.
C
C     Method:     The titles are written using DISSPLA routines.  They
C                 are centered, and positioned at hard-coded distances
C                 from the top and bottom of the page.  Text heights are
C                 also hard-coded.
C                 
C                 CAUTION:  Care must be taken since a new subplot is
C                 defined by the routine, and the page size, plot area and
C                 origin are changed.  CLASSIF should be used after all
C                 other plotting is done.
C
C     Arguments: 
C        ARG         DIM    TYPE   I/O/S    DESCRIPTION
C      ORIENT       LENDIC   C*1     I     'Landscape' or 'Portrait' mode
C      CLASS          *       C      I      Classification title
C      EXPLAN         *       C      I      Explanation title
C      LENTXT         2       I      I      Lengths of titles. Pass 100 to
C                                           invoke self-counting option.
C
C     Procedures (all DISSPLA utilities):
C         AREA2D    Sets the plot area in inches
C         BSHIFT    Used with (0.,0.) to ensure origin is at edge of page
C         ENDGR     Terminates previous subplot
C         HEIGHT    Sets character height
C         MESSAG    Prints text string at specified position
C         PAGE      Sets the page size
C         PHYSOR    Sets the origin
C         XMESS     Returns the length of a character string
C
C     Environment:  VAX/VMS; FORTRAN 77; DISSPLA V11-9003
C
C     History:  03/02/89  M.D.Wong     Initial design and coding.
C               08/27/89   "   "       Updated for DISSPLA V11.0.
C               02/22/90   "   "       Added LENTXT to argument list.
C               08/20/90  D.A.Saunders Recovered DISSPLA version from
C                                      SMDLIB version accidentally installed.
C               03/26/91   "   "       PostScript/Landscape combination was
C                                      giving trouble.  Page size cannot
C                                      exceed 8.5 x 11 it seems.  Use of
C                                      BSHIFT instead (in QPLOT and here)
C                                      got around the problem.
C
C     Author:  Michael Wong, NASA Ames/Sterling Software, Palo Alto, CA.
C
C-------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments.

      INTEGER
     >   LENTXT (2)
      CHARACTER
     >   ORIENT * 1, CLASS * (*), EXPLAN * (*)

C     Local constants.

      REAL
     >   HALF, HTIT, HEXP, SAFETY, SPACE
      PARAMETER
     >  (HALF = 0.5, HTIT = .15, HEXP =.12, SAFETY = .20, SPACE = .10)
      CHARACTER
     >   BLANK * 1
      PARAMETER
     >  (BLANK = ' ')

C     Local variables.

      REAL
     >   BOTTOM, TOP, XCENTER, XPOSN, YPOSN

C     Procedures.

      EXTERNAL
     >   AREA2D, BSHIFT, ENDGR, HEIGHT, MESSAG, PAGE, PHYSOR, XMESS
      REAL
     >   XMESS

C     Execution.

      IF (CLASS .EQ. BLANK .AND. EXPLAN .EQ. BLANK) GO TO 99
   
C     Terminate current subplot.

      CALL ENDGR (0)

C     Set page values according to landscape or portrait mode.

      IF (ORIENT .EQ. 'L') THEN
         CALL BSHIFT (0., 0.)  ! Reset any earlier fudge used to lower the plot
         CALL PAGE (11.0, 8.5)
         XCENTER = 5.75        ! 1/4" fudge for nicer centering of typical plots
         BOTTOM = 0.125        ! 1/8" fudge needed for landscape/Postscript
         TOP = 8.5 - BOTTOM    ! (SAFETY below isn't enough)
      ELSE
         CALL PAGE (8.5, 11.0)      
         XCENTER = 4.50
         BOTTOM = 0.0
         TOP = 11.0
      END IF

C     Define physical origin and new subplot plot area.  Is 6.5 arbitrary?

      CALL PHYSOR (0., 0.)
      CALL AREA2D (6.5, 6.5)

C     Write classification titles in large letters.

      IF (CLASS .NE. BLANK) THEN

         CALL HEIGHT (HTIT)

C        Center first title at top of page, then at bottom of page.
  
         XPOSN = XCENTER - HALF * XMESS (CLASS, LENTXT (1))
         YPOSN = TOP - SAFETY - HTIT
         CALL MESSAG (CLASS, LENTXT (1), XPOSN, YPOSN)

         YPOSN = BOTTOM + SAFETY
         CALL MESSAG (CLASS, LENTXT (1), XPOSN, YPOSN)
      END IF

C     Write explanation titles in smaller letters. 

      IF (EXPLAN .NE. BLANK) THEN

         CALL HEIGHT (HEXP)

         XPOSN = XCENTER - HALF * XMESS (EXPLAN, LENTXT (2))
         YPOSN = TOP - SAFETY - HTIT - SPACE - HEXP
         CALL MESSAG (EXPLAN, LENTXT (2), XPOSN, YPOSN)

         YPOSN = BOTTOM + SAFETY + HTIT + SPACE
         CALL MESSAG (EXPLAN, LENTXT (2), XPOSN, YPOSN)
      END IF

   99 RETURN
      END
C+----------------------------------------------------------------------
C
      SUBROUTINE FIXWINDO (N, X, Y, XMIN, XMAX, YMIN, YMAX, FLAG, EDGES)
C
C ONE-LINER:  FIX up a partially-defined plotting WINDOw
C
C PURPOSE:
C
C        FIXWINDO deals with the problem of determining a well-chosen
C     plotting window which is only partially defined to begin with.
C     Use of the global data limits is not smart enough: the window may
C     be smaller than this because only the data satisfying the partial
C     specification should be eligible, not all of it.
C
C METHOD:
C        Very little is done if the entire window is already specified.
C     Otherwise, the global data limits are determined.  If neither Y
C     limit is specified, then the X window cannot be modified, and
C     this case is handled for efficiency reasons: the X limits are set
C     to the global data extrema for X unless they are specified.
C     The case of neither X limit's being specified is handled analogously.
C
C        If either or both Y limits are specified, the X window is reduced
C     to the extremes of the data points within the Y window, and vice
C     versa.  These cases require more work.
C
C        The strategy is for LAYOUT (or equivalent) to use EDGES (1:2)
C     and EDGES (3:4) in its calls to LINAX instead of the full X (*)
C     and Y (*), thus leaving LINAX intact.  Modularizing the determin-
C     ation of EDGES (*) here leaves LAYOUT (already non-trivial)
C     largely unchanged.
C
C ERROR HANDLING:
C
C        If a specified window edge (in X, say) excludes all the data,
C     the edge(s) being determined (Ymin and/or Ymax) are assigned the
C     the appropriate global data value rather than values suggesting
C     a very narrow window.  This situation is assumed to be handled at
C     a higher level.
C
C ENVIRONMENT:
C     VAX/VMS; FORTRAN 77, with...
C     > IMPLICIT NONE
C     > Trailing ! comments
C     > Names up to 8 characters
C
C HISTORY:
C     10/06/91  D.A.Saunders   Initial implementation to solve a LAYOUT
C                              problem with partially-windowed data.
C
C AUTHOR: David Saunders, Sterling Software/NASA Ames, Mountain View, CA.
C
C ----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments.

      INTEGER   N         ! (I) Number of data points. N >= 1.
      REAL      X (N)     ! (I) Given abscissas ...
      REAL      Y (N)     ! (I) ... and ordinates.
      REAL      XMIN      ! (I) Partial (?) definition of plotting window
      REAL      XMAX      ! (I)    as input to (say) LAYOUT, q.v.
      REAL      YMIN      ! (I)
      REAL      YMAX      ! (I)
      REAL      FLAG      ! (I) Flag for indicating that XMIN, XMAX, etc.
                          !        is not yet specified.
      REAL      EDGES (4) ! (O) Xmin, Xmax, Ymin, Ymax determined as
                          !        defining a well-chosen plot window.

C-----------------------------------------------------------------------

C     Local variables.

      INTEGER   I, I1
      REAL      BOUNDHI, BOUNDLO, WMAX, WMIN, XHI, XLO, YHI, YLO
      LOGICAL   HARDER

C     Execution.

      EDGES (1) = XMIN
      EDGES (2) = XMAX
      EDGES (3) = YMIN
      EDGES (4) = YMAX

      HARDER = .FALSE.
      DO I = 1, 4
         IF (EDGES (I) .EQ. FLAG) HARDER = .TRUE.
      END DO

      IF (HARDER) THEN

C        Find the global data range, which may or may not be narrowed.

         XLO = X (1)
         XHI = XLO
         YLO = Y (1)
         YHI = YLO
         DO 200, I = 2, N
            IF (X (I) .GT. XHI) THEN   ! It can't also be a new minimum...
               XHI = X (I)
            ELSE IF (X (I) .LT. XLO) THEN
               XLO = X (I)
            END IF
            IF (Y (I) .GT. YHI) THEN
               YHI = Y (I)
            ELSE IF (Y (I) .LT. YLO) THEN
               YLO = Y (I)
            END IF
  200    CONTINUE

C        Consider the case of specifying YMIN:
C        XMIN and XMAX should allow for only those points I with
C        Y (I) >= YMIN, not for the full X range.
C        But if neither Y limit were specified, then global data limits
C        for X will do where an X limit was not specified.

         BOUNDLO = EDGES (3)    ! YMIN
         IF (BOUNDLO .EQ. FLAG) BOUNDLO = YLO
         BOUNDHI = EDGES (4)    ! YMAX
         IF (BOUNDHI .EQ. FLAG) BOUNDHI = YHI

         IF (BOUNDLO .GT. YLO .OR. BOUNDHI .LT. YHI) THEN

            ! The Y range was restricted.  The X range may be affected too.

            WMIN = FLAG            ! No valid values can be guaranteed
            WMAX = FLAG
            DO 300, I = 1, N
               I1 = I
               IF (BOUNDLO .LE. Y (I) .AND. Y (I) .LE. BOUNDHI) THEN
                  WMIN = X (I)
                  WMAX = WMIN
                  GO TO 310
               END IF
  300       CONTINUE
  310       DO 320, I = I1 + 1, N
               IF (BOUNDLO .LE. Y (I) .AND. Y (I) .LE. BOUNDHI) THEN
                  WMIN = MIN (WMIN, X (I))
                  WMAX = MAX (WMAX, X (I))
               END IF
  320       CONTINUE

C           Only assign reduced limits if they weren't specified.

            IF (EDGES (1) .EQ. FLAG) EDGES (1) = WMIN  ! May still be undefined.
            IF (EDGES (2) .EQ. FLAG) EDGES (2) = WMAX

         END IF


C        Analogous treatment of Y where X range may or may not be limited.

         BOUNDLO = EDGES (1)
         IF (BOUNDLO .EQ. FLAG) BOUNDLO = XLO
         BOUNDHI = EDGES (2)
         IF (BOUNDHI .EQ. FLAG) BOUNDHI = XHI

         IF (BOUNDLO .GT. XLO .OR. BOUNDHI .LT. XHI) THEN

            ! The X range was restricted.  The Y range may be affected too.

            WMIN = FLAG
            WMAX = FLAG

            DO 500, I = 1, N
               I1 = I
               IF (BOUNDLO .LE. X (I) .AND. X (I) .LE. BOUNDHI) THEN
                  WMIN = Y (I)
                  WMAX = WMIN
                  GO TO 510
               END IF
  500       CONTINUE
  510       DO 520, I = I1 + 1, N
               IF (BOUNDLO .LE. X (I) .AND. X (I) .LE. BOUNDHI) THEN
                  WMIN = MIN (WMIN, Y (I))
                  WMAX = MAX (WMAX, Y (I))
               END IF
  520       CONTINUE

            IF (EDGES (3) .EQ. FLAG) EDGES (3) = WMIN
            IF (EDGES (4) .EQ. FLAG) EDGES (4) = WMAX

         END IF

C        Finally, ensure that all limits are defined.

         IF (EDGES (1) .EQ. FLAG) EDGES (1) = XLO
         IF (EDGES (2) .EQ. FLAG) EDGES (2) = XHI
         IF (EDGES (3) .EQ. FLAG) EDGES (3) = YLO
         IF (EDGES (4) .EQ. FLAG) EDGES (4) = YHI

      END IF

      RETURN
      END
C+----------------------------------------------------------------------
C
      SUBROUTINE GETLINE (LUN, COMMENT, LINE, LAST, IOS)
C
C  One-liner:  Low-level data reader; suppresses trailing comments/blanks
C
C  Description and usage:
C
C        GETLINE is intended to be a standard low level text input utility,
C     providing a uniform input format which permits free use of blank
C     lines, comment lines, trailing comments, and "commented-out" lines.
C     It reads a record from logical unit LUN and returns the "significant"
C     portion in LINE, with only trailing COMMENTs, blanks, or tabs after
C     LAST.  It returns LAST = 0 if the line is effectively empty.
C
C        Double-COMMENTs are replaced with single ones so that COMMENT may
C     still be used in a string if required.  SPECIAL CASE: if COMMENT is
C     the FIRST significant character, the line is considered empty.  This
C     covers the common case of "commenting-out" lines with more than one
C     COMMENT character, as in  !!!! 0.800000   1.
C
C  Arguments:
C
C     Name    Type/Dimension  I/O/S  Description
C
C     LUN     I               I      Logical unit for reading text.
C
C     COMMENT C*1             I      Character signalling the end of the
C                                    "significant" part of a text string
C                                    and the beginning of comments. If
C                                    blank, this feature is disabled and
C                                    much less scanning is needed.
C
C     LINE    C*(*)             O/S  Buffer for reading one record;
C                                    returned with "significant" text in
C                                    LINE (1 : LAST) unless LAST is zero.
C
C     LAST    I                 O    Index of last significant character
C                                    in LINE.  If LINE is null or all
C                                    blanks, tabs, and/or comment, LAST
C                                    is returned as zero.
C
C     IOS     I                 O    Error status of read;  = 0 means no
C                                    read error, negative is EOF, and the
C                                    meaning of positive IOS is system-
C                                    dependent.
C
C  Environment:  Digital VAX-11/780 VMS FORTRAN (FORTRAN 77)
C
C  Notes:
C
C     (1)  IMPLICIT NONE and 7-character names are non-standard.
C
C     (2)  The Q edit descriptor for reading a text string to find the
C          actual length of the record is a VAX extension.  See the
C          commented-out section of code below for a standard approach.
C
C     (3)  Note that GETLINE does NOT keep reading if it encounters an
C          empty line.  Part of the purpose of GETLINE is to avoid this
C          very drawback of FORTRAN's list-directed I/O.
C
C  Authors:  Ronald Langhi/Robert Kennelly,  Sterling Software, Palo Alto, CA
C
C  History:
C
C     10 Mar. 1987  RGL/RAK  Initial design and code as GETSTRING.
C      4 Apr. 1987    RAK    Cosmetics, documentation revised.
C      8 May  1987    DAS    Name changed to GETLINE.
C     22 Apr. 1988  DAS/RAK  Trim trailing tabs as well as the blanks.
C      7 July 1990  DAS/RAK  Handle special case of commenting-out lines
C                            with more than one leading COMMENT.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Constants.

      CHARACTER
     >   BLANK * 1
      PARAMETER
     >   (BLANK = ' ')

C     Arguments.

      INTEGER
     >   LUN, LAST, IOS
      CHARACTER
     >   COMMENT * 1, LINE * (*)

C     Local variables.

      INTEGER
     >   I, J, N              ! N can be LAST throughout, but may be inefficient
      LOGICAL
     >   SUPPRESD
      CHARACTER
     >   HTAB * 1

C     Execution.
C     ----------

C     Read the next line, using the non-standard Q edit descriptor
C     to determine the actual length of the record being read.

      HTAB = CHAR (9)   ! For Absoft FORTRAN on 68000 machines

      READ (LUN, '(Q, A)', IOSTAT = IOS) N, LINE
      IF (IOS .NE. 0) GO TO 99

      IF (N .EQ. 0) THEN
         GO TO 99
      ELSE
         N = MIN (N, LEN (LINE))
      END IF

C     Alternatively, read the text using standard FORTRAN 77.
C
C     READ (LUN, '(A)', IOSTAT = IOS) LINE
C     IF (IOS .NE. 0) GO TO 99
C
C     N = LEN (LINE)

C     Perform search for comments?
C     ----------------------------

      IF (COMMENT .NE. BLANK) THEN

C        Examine one character at a time in LINE (1 : N) for a
C        transition between significant characters and comments.

         I = 0
   10    CONTINUE
            I = I + 1
            IF (LINE (I : I) .EQ. COMMENT) THEN

C              Handle a special case here rather than impact all lines.
C              If the FIRST significant character is a COMMENT, consider
C              the line suppressed.  This covers the common case where
C              more than one COMMENT character is used to "comment-out"
C              an input line.
C
C              Search for preceding significant text:

               SUPPRESD = .TRUE.
               DO 20, J = 1, I-1
                  IF (LINE (J : J) .NE. BLANK .AND.
     >                LINE (J : J) .NE. HTAB) SUPPRESD = .FALSE.
   20          CONTINUE

               IF (SUPPRESD) THEN
                  N = 0
                  GO TO 40       ! For consistency; GO TO 99 is more direct.
               ELSE IF (I .EQ. N) THEN

C                 Single comment at the very end of the text - done.

                  N = I - 1
                  GO TO 40
               ELSE IF (LINE (I+1 : I+1) .EQ. COMMENT) THEN

C                 Double comment - strip the first one from the text,
C                 decrement N, and keep searching.  The DO-loop is
C                 clunky but required by the standard.

                  DO 30, J = I, N - 1
                     LINE(J : J) = LINE(J+1 : J+1)
   30             CONTINUE
                  N = N - 1
               ELSE

C                 Single comment embedded in the text - done.

                  N = I - 1
                  GO TO 40
               END IF
            END IF
            IF (I .LT. N) GO TO 10

   40    CONTINUE
      END IF

C     Remove trailing blanks and tabs.
C     --------------------------------

C     We're taking advantage of the fact that if the DO-loop runs to
C     completion, then N = 0 on exit.

      DO 50, N = N, 1, -1
         IF (LINE (N : N) .NE. BLANK .AND.
     >       LINE (N : N) .NE. HTAB) GO TO 99
   50 CONTINUE

C     Termination.
C     ------------

   99 LAST = N
      RETURN
      END
C+----------------------------------------------------------------------
C
      SUBROUTINE GETPRM (PARAM, LENGTH)
C
C  ACRONYM:  GET command line PaRaMeters
C            ---              - - -
C  PURPOSE:
C     Get the system command line and return into a string the whole
C     command from after the first blank string to end of line.
C     Example:
C
C        > PROG FILENAME.TYP
C
C        GETPRM returns 'FILENAME.TYP' in PARAM (1:LENGTH) where LENGTH=12.
C
C  ARGUMENTS:
C     ARG   TYPE I/O/S  DIM    DESCRIPTION
C     PARAM  C*    O     -     Parameter string.  See NOTES below.
C     LENGTH I     O     -     Length of string.  0 means no parameter.
C
C  COMMONS USED:  None.
C
C  EXTERNAL REFERENCES:
C     IARGC   Unix system routine to get # of command line args
C     GETARG  Unix system routine to get one argument
C
C  ENVIRONMENT:  VAX/Ultrix x.x, FORTRAN 77; IRIS/IRIX, f77
C
C  AUTHORS: Dexter L. Hermstad		Informatics General Corporation
C		07/02/81		Palo Alto, California
C           Daniel S. McCoy		Informatics, Inc.
C
C  REVISIONS:
C     12/10/91   DAS   Dexter pointed out that NP = GETARG (I, ARG) gives
C                      NP = LEN (ARG), not the number of characters in the
C                      Ith argument.
C     09/20/91   DAS   Matched the VAX version, which now allows PARAM
C                      to have variable length.  Changed CALL GETARG ...
C                      to NP = GETARG ...  Added error handling.
C     07/03/90   DLH   Initialize PARAM to ' ' before any processing, to
C                      prevent returning meaningless stray characters.
C     04/06/89   DLH   Version with PARAM as character variable
C     03/01/89   DLH   Ultrix version
C     07/02/81   DLH   Original coding - adapted from PDP-11 version
C     03/??/81   DSM   Design
C
C  NOTES:
C     > Multiple blanks and tabs are compressed to single blanks.
C       If the command line exceeds the input string, execution halts
C       with a message to unit 6.
C     > New applications would be wiser to use Unix's getarg and iargc
C       utilities directly.
C-----------------------------------------------------------------------

C     Arguments:

      CHARACTER  PARAM * (*)
      INTEGER    LENGTH

C     Local constants:

      INTEGER    MXBUF
      CHARACTER  BLANK * 1
      PARAMETER (BLANK = ' ',
     >           MXBUF = 80)   ! Max. length of any one argument

C     Local variables:

      INTEGER    I, LENMX, NP
      INTEGER*4  IARGC, NARG   ! Must be longword, even if compiled /NOI4
      CHARACTER  ARG * (MXBUF)

C ... Get the command line by appending each argument.

      LENMX = LEN (PARAM)
      PARAM = BLANK
      NARG = IARGC ()
      LENGTH = 0
      DO 100, I = 1, NARG
         IF (I .GT. 1) THEN
            LENGTH = LENGTH + 1
            IF (LENGTH .GT. LENMX) GO TO 900
            PARAM (LENGTH:LENGTH) = BLANK
         END IF
         CALL GETARG (I, ARG)
         NP = INDEX (ARG, BLANK) - 1
         IF (LENGTH + NP .GT. LENMX) GO TO 900
         PARAM (LENGTH + 1:LENGTH + NP) = ARG (1:NP)
         LENGTH = LENGTH + NP
  100 CONTINUE
      GO TO 999

  900 WRITE (6, '(A)')
     >   ' GETPRM: Command line argument(s) exceeded input buffer.',
     >   BLANK     ! Because IRIS FORTRAN's "stop" clobbers last line.
      STOP ' '

  999 RETURN
      END
C+----------------------------------------------------------------------
C
      SUBROUTINE INTERVAL (NX, X, XFIND, ARROW, LEFT)
C
C     One-liner: Interpolation search for interval containing a point.
C     ----------
C
C     Description and usage:
C     ----------------------
C
C        Written primarily for interval-based interpolations such as
C     piecewise linear or cubic spline, INTERVAL performs a search to
C     locate the best interval for evaluating the interpolant at a
C     given point. The normal case returns the "left-hand" endpoint of
C     the interval bracketing the point, but for the out-of-range cases
C     below or above the range of the knots, the interval to be used is
C     the first or last. The array of knots must be monotonic, either
C     increasing or decreasing. Diagrammatically, LEFT is returned as
C     shown below for the normal case (no extrapolation):
C
C          X (1)  ...   X (LEFT)   X (LEFT+1)   ...      X (NX)
C                                ^
C                              XFIND
C
C     And for extrapolation:
C
C                     X (LEFT = 1)  ...   X (NX)
C             ^
C           XFIND
C
C     or,
C                X (1)  ...   X (LEFT = NX-1)    X (NX)
C                                                           ^
C                                                         XFIND
C
C     If the point to be bracketed (XFIND) matches one of the knots, the
C     index of that knot is returned as LEFT, i.e., the condition for a
C     bracket of an interior point is:
C
C        X (LEFT) <= XFIND < X (LEFT+1)  if  ARROW = +1.0,  or
C        X (LEFT) >= XFIND > X (LEFT+1)  if  ARROW = -1.0.
C
C        This is a low-level routine with minimal error checking. The
C     calling program is assumed to have verified the following:
C
C     (1)  NX >= 2
C     (2)  X strictly monotonic
C     (3)  ARROW = +1.0 or -1.0
C
C     Subroutine PROTECT is available from the author for easily checking
C     conditions (2) and (3). LEFT is verified on input, but efficiency in
C     loops will benefit from passing the best estimate available, usually
C     just the result of the last call.
C
C        INTERVAL was originally written for use with CSEVAL and TABLE1.
C     The interpolation search was adapted from ideas in Sedgewick's book
C     referenced below.
C
C     Arguments:
C     ----------
C
C     Name  Dimension  Type  I/O/S  Description
C     NX                I    I      Number of points in array X; must
C                                   be >= 2 (no check performed).
C
C     X        NX       R    I      Array of points defining the set
C                                   of intervals to be examined. Only
C                                   the first NX-1 points are required.
C
C     XFIND             R    I      The point for which a bracketing
C                                   interval is sought.
C
C     ARROW             R    I      Monotonicity indicator for input
C                                   array X:
C                                     -1.0  strictly decreasing
C                                      0.0  NOT ALLOWED!
C                                     +1.0  strictly increasing
C                                   Supplied by the calling routine for
C                                   reasons of speed (not checked).
C
C     LEFT              I    I/O    Input: guessed index of left-hand
C                                   endpoint of the interval containing
C                                   the specified point.
C
C                                   Output: index of the largest array
C                                   value <= specified point (if ARROW=+1.0).
C                                   Special case for data out of range:
C                                   return left endpoint of closest interval.
C                                   Thus, LEFT = 1 for XFIND < X (2), and
C                                   LEFT = NX-1 for XFIND >= X (NX-1).
C                                   (If ARROW=-1.0, reverse the inequalities.)
C
C     Environment:  Digital VAX-11/780, VMS FORTRAN
C     ------------  Apple Macintosh, Absoft MacFORTRAN/020 v2.3
C
C     Notes:
C     ------
C
C     (1)  IMPLICIT NONE and eight character symbols are not (yet) standard.
C
C     (2)  In speed-critical applications, it might be a good idea to build
C          this algorithm in-line since it is typically called many times
C          from within a loop. Another potential speed-up is removal of the
C          ARROW multiplies, which restricts the method to increasing data.
C          So far, the simplicity of separating out the messy search details
C          and the generality of bi-directional searching have outweighed
C          the modest speed penalty incurred.
C
C     Bibliography:
C     -------------
C
C     (1) Sedgewick, R.  Algorithms.  Reading: Addison-Wesley, 1983.
C            (Chap. 14)
C
C     Author:  Robert Kennelly and David Saunders, Sterling Federal Systems
C     -------
C
C     Development history:
C     --------------------
C
C     20 Oct. 1987    RAK    Interpolation search adapted (with mods.
C                            for bidirectional search and some minor
C                            repair) from CSEVAL (RAK) and TABLE1 (DAS).
C     08 Aug. 1988    DAS    Clarified descriptions of bracketing, where
C                            the inequalities depend upon ARROW.
C
C-----------------------------------------------------------------------

C     Declarations.
C     -------------

      IMPLICIT NONE

C     Constants.

      REAL
     &   ONE
      PARAMETER
     &  (ONE = 1.0E+0)

C     Arguments.

      INTEGER
     &   LEFT, NX
      REAL
     &   ARROW, X (NX), XFIND

C     Local variables.

      INTEGER
     &   LENGTH, NXLESS1, RIGHT, TRIAL
      REAL
     &   XBYARROW

C     Execution.
C     ----------

      XBYARROW = XFIND * ARROW

C     Simplify things by disposing of two important special cases so that
C     X (LEFT) and X (RIGHT) can really bracket XFIND. As a by-product,
C     this also takes care of the NX = 2, 3 cases.

      NXLESS1 = NX - 1

      IF (XBYARROW .GE. X (NXLESS1) * ARROW) THEN
         LEFT = NXLESS1
         GO TO 990
      ELSE IF (XBYARROW .LT. X (2) * ARROW) THEN
         LEFT = 1
         GO TO 990
      END IF

C      ---------------------------------
C     |                                 |
C     |   X (2) <= XFIND < X (NX - 1)   |
C     |            - or -               |
C     |   X (2) > XFIND >= X (NX - 1)   |
C     |                                 |
C     |   NX > 3                        |
C     |                                 |
C      ---------------------------------

C     Adjust the pointers. We hope that the calling routine has provided
C     a reasonable guess (since it's probably working on an ordered array
C     of points to evaluate), but check anyway.

      LEFT = MIN (MAX (2, LEFT), NX - 2)

      IF (XBYARROW .GE. X (LEFT) * ARROW) THEN
         IF (XBYARROW .LT. X (LEFT + 1) * ARROW) THEN

C           XFIND is in the original guessed-at interval.

            GO TO 990
         ELSE

C           We'll look farther to the right. Evidently LEFT was < NX - 2.

            RIGHT = NXLESS1
            LEFT  = LEFT + 1
         END IF
      ELSE

C        Look to the left of the guess. Evidently LEFT was > 2.

         RIGHT = LEFT
         LEFT  = 2
      END IF

C      ----------------------------------
C     |                                  |
C     |   2 <= LEFT < RIGHT <= NX - 1    |
C     |                                  |
C      ----------------------------------

C     The interval length must decrease each time through - terminate
C     when the correct interval is found or when the interval length
C     cannot be decreased.

   10 CONTINUE
         LENGTH = RIGHT - LEFT
         IF (LENGTH .GT. 1) THEN

C           The trial value is a "linear" estimate of the left-hand endpoint
C           of the interval bracketing the target XFIND, with protection
C           against round-off (which can affect convergence).

            TRIAL = MIN (RIGHT - 1, LEFT + MAX (0, INT (REAL (LENGTH) *
     &         (XFIND - X (LEFT)) / (X (RIGHT) - X (LEFT)))))

C            ------------------------------------------
C           |                                          |
C           |   2 <= LEFT <= TRIAL < RIGHT <= NX - 1   |
C           |                                          |
C            ------------------------------------------

C           Adjust pointers. Increase LEFT or decrease RIGHT until done.

            IF (XBYARROW .GE. X (TRIAL + 1) * ARROW) THEN
               LEFT  = TRIAL + 1
            ELSE IF (XBYARROW .LT. X (TRIAL) * ARROW) THEN
               RIGHT = TRIAL
            ELSE

C              We're done: XFIND is in the interval [X (TRIAL), X (TRIAL+1)).

               LEFT  = TRIAL
               GO TO 990
            END IF
            GO TO 10

         END IF

C     Termination.
C     ------------

  990 CONTINUE
      RETURN
      END
C+----------------------------------------------------------------------
C
      SUBROUTINE IOCHEK (IOS, OK, EOF, LDSYN, CONVER)
C
C
C     Description:
C
C           IOCHEK handles the machine dependent details of interpreting
C        the IO status flag set during READ statements in FORTRAN 77.  The 
C        LOGICAL flags returned may be more convenient to deal with in
C        the calling program than the integer values of the status flag.
C
C
C     Arguments:
C
C        Name    Dimension   Type   I/O/S   Description
C        IOS                  I     I       IO status flag from a READ.
C        OK                   L       O     Set .TRUE. if no READ errors.
C        EOF                  L       O     Set .TRUE. if end-of-file was
C                                           encountered.
C        LDSYN                L       O     Set .TRUE. if an input list-
C                                           directed I/O syntax error
C                                           occurred.
C                                           [Seems impossible on the IRIS 4D.]
C        CONVER               L       O     Set .TRUE. if an input conver-
C                                           sion error occurred.
C
C
C     Environment:  IRIS 4D; FORTRAN 77
C
C
C     Notes:
C
C        (1)  Consult local system documentation for appropriate values for
C             the return status flag.  These values are MACHINE DEPENDENT.
C
C        (2)  IMPLICIT NONE is non-standard.
C
C
C     Author:  Robert Kennelly, Informatics General Corporation.
C
C
C     History:
C
C        18 Dec. 1982    RAK    Initial coding.
C        27 Dec. 1982    RAK    Cleaned up initialization.
C        13 June 1983    RAK    Minor rearrangement, used LDSYN instead
C                               of (incorrect) CONVER error flag name.
C        29 June 1983    CLH    Added CONVER.
C        02 Feb. 1990    DAS    IRIS 4D version.  Unknown code for ILDSYN.
C
C-----------------------------------------------------------------------


C     Specifications.
C     ---------------

      IMPLICIT NONE


C     Parameters.

      INTEGER
     >   IOS, IOK, IEOF, ILDSYN, ICNVER
      PARAMETER
     >  (IOK = 0, IEOF = -1, ILDSYN = 59, ICNVER = 115)
C                            59 is VAX value - always get 0 on IRIS 4D.

C     Variables.

      LOGICAL
     >   OK, EOF, LDSYN, CONVER


C     Control.
C     --------

C     Set output flags.

      OK = IOS .EQ. IOK
      EOF = IOS .LE. IEOF
      LDSYN = IOS .EQ. ILDSYN
      CONVER = IOS .EQ. ICNVER

      RETURN
      END
C+----------------------------------------------------------------------
C
      SUBROUTINE LAYOUT (LUNOUT, N, X, Y, XMIN, XMAX, XSTEP, YMIN, YMAX,
     >   YSTEP, HEIGHT, WIDTH, DEFHGT, DEFWID, PLOT)
C
C     One-liner:  Detailed plot scaling/sizing, with conflict resolution
C     ----------
C
C     Description and usage:
C     ----------------------
C
C        Detailed data plot layout is handled in this one module in order
C     to hide the messy details from the rest of an application.  LAYOUT
C     was originally written for use by QPLOT, which is based on the DISSPLA
C     proprietary graphics package.  This routine is more powerful (for
C     linear axes, at least) than the setup utilities provided by DISSPLA,
C     so it has been re-written for general use.  Its functions include:
C
C        (1) defaulting,
C        (2) input decoding,
C        (3) conflict resolution, and
C        (4) error checking/correction
C
C     of the axis scaling, physical size, and plot type.  Dictionary
C     look-ups permit abbreviations and synonyms for the string type
C     inputs.  A fairly powerful linear axis scaling routine will handle
C     nearly any combination of inputs, with (optional) error reporting
C     on the few cases where no sensible interpretation of the input data
C     is possible.  Some error checking is applied to logarithmic axes,
C     while polar plots are handled in a rudimentary fashion simply by
C     passing the data straight through.  (These limitations on log and
C     polar plots are a feature which should be improved if LAYOUT is to
C     be used with graphics packages other than DISSPLA.)
C
C        Arguments for which the application expects LAYOUT to determine
C     values should be FLAGged as 999. on input.
C
C        If your data is not in single X and Y arrays, or if you want to
C     suggest values for minimum and maximum, pass in arrays X and Y of
C     length 2 containing the desired extremes.
C
C        See the Notes below for some restrictions on LAYOUT's generality.
C     The QPLOT header contains detailed information on that program's
C     control inputs.  Dictionary routine LOOKUP is part of a small set
C     of FORTRAN 77 input utilities available from the author.
C
C     Arguments:
C     ----------
C
C     Name    Dimension  Type  I/O/S  Description
C     LUNOUT              I    I      Logical unit number for diagnostic
C                                     output.  Use < 0 to suppress output.
C
C     N                   I    I      Number of (X,Y) pairs to be plotted.
C
C     X         N         R    I      Packed array of abcissas.
C
C     Y         N         R    I      Packed array of ordinates.
C
C     XMIN                R    I/O    Left-hand endpoint of the horizontal
C                                     axis, in plot units. Use XMIN = FLAG
C                                     to indicate that a value is to be
C                                     calculated.
C
C     XMAX                R    I/O    Right-hand endpoint of the horizontal
C                                     plot axis.  Use FLAG to indicate
C                                     that a value is to be calculated.
C
C     XSTEP               R    I/O    Length of axis divisions on the
C                                     horizontal axis.  Use FLAG to indicate
C                                     that a value is to be calculated.
C                                     Negative FLAG means auto-scale with
C                                     decreasing axis.
C
C     YMIN                R    I/O    Left-hand endpoint of the vertical
C                                     plot axis.  Use FLAG to indicate
C                                     that a value is to be calculated.
C
C     YMAX                R    I/O    Right-hand endpoint of the vertical
C                                     plot axis.  Use FLAG to indicate
C                                     that a value is to be calculated.
C
C     YSTEP               R    I/O    Length of axis divisions on the
C                                     vertical axis.  Use FLAG to indicate
C                                     that a value is to be calculated.
C                                     Negative FLAG means auto-scale with
C                                     decreasing axis.
C
C     HEIGHT              R    I/O    User-requested length of the vertical
C                                     axis, in physical units.  This value
C                                     will be overridden if too small (see
C                                     parameter DENSITY = max. no. of scale
C                                     divisions per unit), or in case of
C                                     conflicts with equally-scaled plot
C                                     axes. If HEIGHT = FLAG, then the
C                                     default value is used (see DEFHGT).
C
C     WIDTH               R    I/O    User-requested length of the horizontal
C                                     axis, in physical units.  See HEIGHT.
C
C     DEFHGT              R    I      Default length of the vertical axis,
C                                     in physical units.
C
C     DEFWID              R    I      Default length of the horizontal axis,
C                                     in physical units.
C
C     PLOT                C    I/O    Plot type identifier.  After checking
C                                     dictionary, converted to one of the
C                                     standard forms:  'LINEAR', 'LOGLOG',
C                                     'LOGX', 'LOGY', 'POLAR', and 'SCALE',
C                                     whose meanings should be evident. The
C                                     default is 'LINEAR'.  Note that 'SCALE'
C                                     is changed to 'LINEAR' after the axes
C                                     have been set up successfully.
C
C     External references:
C     --------------------
C
C     FIXWINDO  Determines X range left by Y windowing, and vice versa
C     LINAX     Chooses (and error checks) axis limits and step sizes
C     LINCHK    Analyzes errors reported by LINAX, "corrects" some of them
C     LOGAX     Skeleton routine only, for now (uses DISSPLA)
C     LOOKUP    Dictionary lookup, with abbreviations and synonyms
C     UPCASE    Used to ensure PLOT is in upper case
C
C     Environment:  Digital VAX-11/780 VMS FORTRAN (FORTRAN 77)
C     ------------
C
C     Notes:
C     ------
C
C     (1)  IMPLICIT NONE, 8-character symbolic names, and use of "!" for
C          comments is not standard FORTRAN 77.
C
C     (2)  The SCALE option cannot, in general, respect all user choices
C          of X- and Y-axis length.  We have chosen to base the plot
C          layout on the physical length input for that axis whose length
C          is input, or the Y-axis in case of conflict.  There is no guar-
C          antee that the resulting layout will fit on the eventual output
C          medium!
C
C     (3)  As above, a user-supplied value for XSTEP will be ignored if it
C          conflicts with YSTEP in those cases (the majority) where the
C          Y-axis information takes precedence.
C
C     (4)  Only rudimentary range checking is applied to HEIGHT and WIDTH:
C          they must be at least one minimum-grid-space.  See parameter
C          DENSITY, below.  It may be desirable to add upper bounds as well.
C
C     (5)  For log axes, we re-use the STEP variables for passing DISSPLA's
C          CYCLE (inches/cycle) quantities.  This may be inappropriate (or
C          merely meaningless) for use with other graphics packages.
C
C     (6)  The FLAG value used here for defaulting (+/- 999) is arbitrary,
C          but must be consistent with the other routines.  The idea is NOT
C          to force the user to use FLAG to obtain defaults, but as a means
C          of communication between routines:  the calling program should
C          just preset the variable to FLAG, then read in the user's value
C          and pass it to LAYOUT.
C
C     (7)  For SCALE plots, we use the same STEP on both axes if possible,
C          but if the derived axis' STEP is too big we have to re-try.  The
C          result will be equally scaled, but may not look it.
C
C     (8)  The physical units used are expected to be inches, as reflected
C          in choice of DENSITY (equivalent to 3/4" inch per scale division)
C          and GRACE (calculated scaled axis lengths may exceed the default
C          values by 1").  These two constants will need to be changed if
C          another system of measurement is to be used.
C
C     History:
C     --------
C
C     11 Apr. 1985    RAK     Initial design and coding.
C      6 Aug. 1985    RAK     SCALE option adapted from EQUSTP/GRNICE
C                             approach used in earlier QPLOT.
C     12 Aug. 1985    RAK     Added LINCHK to decode error messages from
C                             LINAX (modular, saves space).  Moved
C                             assignment of default axis lengths down
C                             from QPLOT level.
C     21 Aug. 1985    RAK     Tie minimum HEIGHT and WIDTH to DENSITY.
C                             Print informational message for SCALE plots
C                             telling which axis controls the layout.
C                             Round down (INT) to ensure that the maximum
C                             grid density does not exceed DENSITY.  Added
C                             ORIENTation parameter.  Default HEIGHT and
C                             WIDTH depend on ORIENT.
C      9 Sep. 1985    RAK     Decreased DENSITY.  Shuffled diagnostic
C                             WRITEs.  Permit XSTEP to govern Y-axis in
C                             SCALE plots if YSTEP is not input, but
C                             retain preference for positive steps.
C     27 Sep. 1985    RAK     Dropped 'CONVERGENCE' plot option (just
C                             point back to 'LOGY').  Re-set 'SCALE' plot
C                             type to 'LINEAR' when through with setup.
C                             'POLAR' plots are passed through with no
C                             checking, and only radial axis is set up.
C      9 Oct. 1985    RAK     Nasty special case:  if MIN, MAX, and STEP
C                             are all specified for SCALE plot, but
C                             STEP is too big, then default STEP and try
C                             again (since we now know that MXNSTP can
C                             be computed).  Protect all MXNSTPs with
C                             MAX (*, 1), esp. for POLAR.
C     16 Oct. 1985    RAK     Try to retain the newly-significant sign
C                             when a STEP input equals -FLAG (indicating
C                             a decreasing axis).  Orientation dictionary
C                             was disoriented (LANDSCAPE out of place).
C     12 Mar. 1987    RAK     Renamed LAYOUT and (decreed to be general-
C                             purpose) moved to graphics library. Changed
C                             orientation terminology to standard
C                             PORTRAIT (default) and LANDSCAPE.
C     24 Mar. 1987  R.Langhi  Removed orientation option, and replaced
C                             it with default axis lengths input by the
C                             calling routine.
C     23 Apr. 1987    RAK     Revised error message spacing (new name).
C                             POLAR step lengths are now rounded values.
C     06 Oct. 1991 D.Saunders Introduced FIXWINDO to solve long-standing
C                             difficulty caused by partially-specified
C                             plot windows, where the global data extremes
C                             are not always good choices - an indicated
C                             window in X, say, may exclude quite a range
C                             of Y, and vice versa.
C
C     Author:  Robert Kennelly, NASA Ames Research Center, Mt. View, CA.
C     -------
C-----------------------------------------------------------------------

C     Declarations.
C     -------------

      IMPLICIT NONE

C     Constants.

      LOGICAL
     >   ORDER
      PARAMETER
     >  (ORDER = .TRUE.)

      INTEGER
     >   LENDIC, NPLDIC
      PARAMETER
     >  (LENDIC = 14,                 ! Length of plot dictionary entries
     >   NPLDIC = 10)                 ! Number of dictionary entries

      REAL
     >   DENSITY, FLAG, GRACE, MINLEN
      PARAMETER
     >  (DENSITY = 4.0 / 3.0,         ! Max. number of steps/inch
     >   FLAG    = 999.0,             ! Special value for indicating defaults
     >   GRACE   = 1.0,               ! Default SCALE plot "expansion" limit
     >   MINLEN  = 1.0 / DENSITY)     ! Min. length of an axis step

C     Arguments.

      INTEGER
     >   LUNOUT, N
      REAL
     >   DEFHGT, DEFWID, HEIGHT, WIDTH, X (N), XMAX, XMIN, XSTEP, Y (N),
     >   YMAX, YMIN, YSTEP
      CHARACTER
     >   PLOT * (*)

C     Local variables.

      LOGICAL
     >   CHECK
      INTEGER
     >   I, ITEM, MXNSTP, OOPS
      REAL
     >   EDGES (4), URXMAX, URXMIN, URXSTP, URYMAX, URYMIN, URYSTP
      CHARACTER
     >   PLDIC (NPLDIC) * (LENDIC)

C     Storage.

      DATA
     >  (PLDIC (I), I = 1, NPLDIC)
     >     /'DEFAULT=LINEAR',
     >      'LINEAR',
     >      'LOGLOG',
     >      'LOGX',
     >      'LOGY',
     >      'POLAR',
     >      'SCALE',
     >      'XLOG=LOGX',
     >      'XY=LINEAR',
     >      'YLOG=LOGY'/

C     Execution.
C     ----------

      HEIGHT = MAX (ABS (HEIGHT), MINLEN)
      WIDTH  = MAX (ABS (WIDTH),  MINLEN)

C     Fix up plot type flag.

      CALL UPCASE (PLOT)
      CALL LOOKUP (NPLDIC, PLDIC, ORDER, PLOT, ITEM)
      IF (ITEM .LE. 0) THEN
         IF (LUNOUT .GT. 0) WRITE (LUNOUT, 1000) ITEM
         PLOT = 'LINEAR'
      END IF

C     Save inputs so that the plot scale may be re-computed later.  (Prefix
C     "UR" means "original" in German, if that helps...)

      URXMIN = XMIN
      URXMAX = XMAX
      URXSTP = XSTEP

      URYMIN = YMIN
      URYMAX = YMAX
      URYSTP = YSTEP

C     Handle the problem of partially-specified windows: using global data
C     extrema for the unspecified parts of the plot window is not good
C     enough because the given clipping in X (say) may exclude some of the
C     Y range, and vice versa.

      CALL FIXWINDO (N, X, Y, URXMIN, URXMAX, URYMIN, URYMAX, FLAG,
     >               EDGES)

C     EDGES (1:2) can now replace X (*) in LINAX calls, and
C     EDGES (3:4) can now replace Y (*) similarly.


C     The 'SCALE' plot option can fail:  with default settings, a derived
C     X-axis can be too long, i.e., WIDTH much greater than DEFWID.  This
C     loop allows a second try, with the Y-axis derived from the X-axis.

   10 CONTINUE
         CHECK = .FALSE.

         IF (PLOT .EQ. 'SCALE') THEN

C           Find step size and axis lengths for equal axis scaling.
C           -------------------------------------------------------

C           Axis lengths and stepsizes are controlled by X-axis input only
C           if the user explicitly requests them, and does NOT explicitly
C           request a corresponding Y-axis value.

            IF (WIDTH .NE. FLAG .AND. HEIGHT .EQ. FLAG) THEN

C              Use X-axis length to determine required length of Y-axis.

               MXNSTP = MAX (INT (WIDTH * DENSITY), 1)
               IF (ABS (YSTEP) .NE. FLAG) XSTEP = ABS (YSTEP)
               IF (LUNOUT .GT. 0) WRITE (LUNOUT, 1010) 'X'
   20          CONTINUE
                  CALL LINAX (2, EDGES (1), MXNSTP, XMIN, XMAX, XSTEP,
     >                        OOPS)
                  IF (OOPS .NE. 0) THEN
                     CALL LINCHK (LUNOUT, OOPS, 'X', XMIN, XMAX, XSTEP)
                     GO TO 20

                  END IF

C              Since the two axes usually use the same step, specifying the
C              maximum number of steps for the derived axis may not be
C              required.  But if LINAX returns ERROR = +2, then the step
C              derived above is too big and we will have to try again, using
C              computed MXNSTP (YMAX and YMIN must have been fixed).

               MXNSTP = 999
               YSTEP  = SIGN (XSTEP, YSTEP)
   30          CONTINUE
                  CALL LINAX (2, EDGES (3), MXNSTP, YMIN, YMAX, YSTEP,
     >                        OOPS)
                  HEIGHT = WIDTH * ABS ((YMAX - YMIN) / (XMAX - XMIN))
                  IF (OOPS .NE. 0) THEN
                     CALL LINCHK (LUNOUT, OOPS, 'Y', YMIN, YMAX, YSTEP)
                     IF (OOPS .EQ. +2)
     >                  MXNSTP = MAX (INT (HEIGHT * DENSITY), 1)
                     GO TO 30

                  END IF

            ELSE

C              Use Y-axis length to determine required length of X-axis.

               IF (HEIGHT .EQ. FLAG) THEN
                  HEIGHT = DEFHGT
                  IF (WIDTH .EQ. FLAG) CHECK = .TRUE.
               END IF

               MXNSTP = MAX (INT (HEIGHT * DENSITY), 1)
               IF (ABS (YSTEP) .EQ. FLAG) YSTEP = SIGN (XSTEP, YSTEP)
               IF (LUNOUT .GT. 0) WRITE (LUNOUT, 1010) 'Y'
   40          CONTINUE
                  CALL LINAX (2, EDGES (3), MXNSTP, YMIN, YMAX, YSTEP,
     >                        OOPS)
                  IF (OOPS .NE. 0) THEN
                     CALL LINCHK (LUNOUT, OOPS, 'Y', YMIN, YMAX, YSTEP)
                     GO TO 40

                  END IF

C              See comment above on re-try's when a user-specified SCALE
C              plot must resort to different STEPs on X- and Y-axes.

               MXNSTP = 999
               XSTEP  = SIGN (YSTEP, XSTEP)
   50          CONTINUE
                  CALL LINAX (2, EDGES (1), MXNSTP, XMIN, XMAX, XSTEP,
     >                        OOPS)
                  WIDTH = HEIGHT * ABS ((XMAX - XMIN) / (YMAX - YMIN))
                  IF (OOPS .NE. 0) THEN
                     CALL LINCHK (LUNOUT, OOPS, 'X', XMIN, XMAX, XSTEP)
                     IF (OOPS .EQ. +2)
     >                  MXNSTP = MAX (INT (WIDTH * DENSITY), 1)
                     GO TO 50

                  END IF

            END IF         
         ELSE

C           Treat the X- and Y-axes individually (unequal scaling).
C           -------------------------------------------------------

C           Set axis sizes to default values if requested.

            IF (WIDTH  .EQ. FLAG) WIDTH  = DEFWID
            IF (HEIGHT .EQ. FLAG) HEIGHT = DEFHGT

            IF (PLOT .EQ. 'LINEAR' .OR. PLOT .EQ. 'LOGY') THEN

C              Linear X-axis.

               MXNSTP = MAX (INT (WIDTH * DENSITY), 1)
   60          CONTINUE
                  CALL LINAX (2, EDGES (1), MXNSTP, XMIN, XMAX, XSTEP,
     >                        OOPS)
                  IF (OOPS .NE. 0) THEN
                     CALL LINCHK (LUNOUT, OOPS, 'X', XMIN, XMAX, XSTEP)
                     GO TO 60

                  END IF

            ELSE IF (PLOT .EQ. 'LOGX' .OR. PLOT .EQ. 'LOGLOG') THEN

C              Logarithmic X-axis.

               CALL LOGAX (LUNOUT, 2, EDGES (1), WIDTH, XMIN, XMAX,
     >                     XSTEP, OOPS)
            END IF

            IF (PLOT .EQ. 'LINEAR' .OR. PLOT .EQ. 'LOGX') THEN

C              Linear Y-axis.

               MXNSTP = MAX (INT (HEIGHT * DENSITY), 1)
   70          CONTINUE
                  CALL LINAX (2, EDGES (3), MXNSTP, YMIN, YMAX, YSTEP,
     >                        OOPS)
                  IF (OOPS .NE. 0) THEN
                     CALL LINCHK (LUNOUT, OOPS, 'Y', YMIN, YMAX, YSTEP)
                     GO TO 70

                  END IF

            ELSE IF (PLOT .EQ. 'LOGY' .OR. PLOT .EQ. 'LOGLOG') THEN

C              Logarithmic Y-axis.

               CALL LOGAX (LUNOUT, 2, EDGES (3), HEIGHT, YMIN, YMAX,
     >            YSTEP, OOPS)

            ELSE IF (PLOT .EQ. 'POLAR') THEN

C              Polar diagram of X (THETA, in radians) vs. Y (RADIUS).

C              Required step parameter is units/inch, based on an axis
C              of length = (half of the minimum dimension).  THIS IS AN
C              ODDBALL CASE - DISSPLA always uses one step/inch for POLAR.
C              NOTE: since YMAX will be recomputed later by DISSPLA, the
C              final plot may not fill the available space exactly.

               MXNSTP = MAX (INT (0.50 * MIN (HEIGHT, WIDTH)), 1)

C              Default minimum radius is zero - and note that we can't
C              specify a decreasing R-axis using DISSPLA's POLAR routine.

               IF (YMIN .EQ. FLAG) YMIN = 0.0
               YMAX  = FLAG
               YSTEP = ABS (YSTEP)
   80          CONTINUE
                  CALL LINAX (2, EDGES (3), MXNSTP, YMIN, YMAX, YSTEP,
     >                           OOPS)
                  IF (OOPS .NE. 0) THEN
                     CALL LINCHK (LUNOUT, OOPS, 'R', YMIN, YMAX, YSTEP)
                     GO TO 80

                  END IF
            END IF
         END IF

         IF (CHECK .AND. WIDTH .GT. DEFWID + GRACE) THEN

C           A default 'SCALE' plot layout has failed - reset inputs,
C           specifying WIDTH explicitly so that the X-axis will dominate
C           on retry.  (Only one retry will be attempted.)

            HEIGHT = FLAG
            WIDTH  = DEFWID

            XMIN   = URXMIN
            XMAX   = URXMAX
            XSTEP  = URXSTP

            YMIN   = URYMIN
            YMAX   = URYMAX
            YSTEP  = URYSTP

            IF (LUNOUT .GT. 0) WRITE (LUNOUT, 1020)
            GO TO 10

         END IF

C     Once the axes have been set up, a SCALE plot is the same as LINEAR.

      IF (PLOT .EQ. 'SCALE') PLOT = 'LINEAR'

C     Termination.
C     ------------

      RETURN

C     Formats.
C     --------

 1000 FORMAT (' LAYOUT:  Warning - the plot type was improperly ',
     >   'specified.  (ITEM = ', I3, ')'/
     >   10X, 'The default will be used (LINEAR).')
 1010 FORMAT (' LAYOUT:  Axis lengths and scaling will be derived from',
     >   ' the ', A, '-axis information.')
 1020 FORMAT (' LAYOUT:  The derived X-axis length was much longer ',
     >   'than the default width.'/
     >   10X, 'The layout will be repeated with derived Y-axis length.')

      END
C+----------------------------------------------------------------------
C
      SUBROUTINE LINCHK (LUNOUT, ERROR, LABEL, AXMIN, AXMAX, AXSTEP)
C
C     One-liner:  Error checking/recovery for axis scaling routine LINAX
C     ----------
C
C     Description and usage:
C     ----------------------
C
C        Modularizes the error checking left out of LINAX.  Normal usage,
C     if error checking/correction is desired, is to call LINCHK after
C     LINAX if ERROR is non-zero, and then call LINAX again using the
C     patched axis input values (see LAYOUT for an example).  A hard STOP
C     is provided as a debugging aid if LINAX was called improperly and
C     no sensible recovery is possible.
C
C        LINCHK was designed as part of QPLOT, a general purpose plotting
C     package developed at NASA - Ames Research Center.
C
C     Arguments:
C     ----------
C
C     Name    Dimension  Type  I/O/S  Description
C     LUNOUT              I    I      Logical unit number for error
C                                     messages (line printer format).
C                                     Use < 0 to suppress output.
C
C     ERROR               I    I      Error code from LINAX to be
C                                     analyzed.
C
C     LABEL     *         C    I      Name of the axis being set up, e.g.
C                                     'X', for error messages.
C
C     AXMIN               R    I/O    Axis ORIGIN from LINAX.
C
C     AXMAX               R    I/O    Axis ENDPT from LINAX.
C
C     AXSTEP              R    I/O    Stepsize from LINAX.
C
C
C     Environment:  Digital VAX-11/780 VMS FORTRAN (FORTRAN 77)
C     ------------
C
C     Notes:
C     ------
C
C        (1)  IMPLICIT NONE is non-standard.
C
C
C     History:
C     --------
C
C        12 Aug. 1985    RAK    Initial design and coding.
C        14 Aug. 1985    RAK    For ERROR=+2 (bad AXSTEP), just retry
C                               with the sign reversed.
C        21 Aug. 1985    RAK    Had to distinguish between bad sign (+3)
C                               and bad length of AXSTEP (+2).  Reordered
C                               calling sequence.
C        24 Dec. 1985    RAK    Length of LABEL is now as-passed.  Try
C                               SP sign control editing in error message.
C        27 Mar. 1987    RAK    Negative LUNOUT suppresses messages. Use
C                               A, A in FORMAT to avoid concatenation.
C
C     Author:  Robert Kennelly, NASA Ames Research Center, Mt. View, CA.
C     -------
C-----------------------------------------------------------------------

C     Declarations.
C     -------------

      IMPLICIT NONE

C     Constants.

      REAL
     >   FLAG
      PARAMETER
     >  (FLAG = 999.0)

C     Arguments.

      INTEGER
     >   ERROR, LUNOUT
      REAL
     >   AXMIN, AXMAX, AXSTEP
      CHARACTER
     >   LABEL * (*)

C     Execution.
C     ----------

      IF (ERROR .EQ. +1) THEN

C        Axis length was zero!

         IF (LUNOUT .GT. 0) WRITE (LUNOUT, 1000) LABEL, '-axis', ERROR
         AXMIN = FLAG
         AXMAX = FLAG

      ELSE IF (ERROR .EQ. +2) THEN

C        Input AXSTEP was zero, or too big.

         IF (LUNOUT .GT. 0) WRITE (LUNOUT, 1000) LABEL, 'STEP', ERROR
         AXSTEP = FLAG

      ELSE IF (ERROR .EQ. +3) THEN

C        The step had the wrong sign - we'll "correct" it.

         IF (LUNOUT .GT. 0) WRITE (LUNOUT, 1000) LABEL, 'STEP', ERROR
         AXSTEP = -AXSTEP

      ELSE IF (ERROR .LT. 0) THEN

C        Programming error:  the number of data points or maximum number
C        of steps passed to LINAX was not positive.

         IF (LUNOUT .GT. 0) WRITE (LUNOUT, 1010) LABEL
         STOP 'LINCHK:  Fatal error!'

      END IF

C     Termination.
C     ------------

      RETURN

C     Formats.
C     --------

 1000 FORMAT (' LINCHK:  Warning - bad ', A, A, ' information.  ',
     >   '(ERROR = ', SP, I2, ')'/
     >   10X, 'Automatic scaling will be used.')
 1010 FORMAT (' LINCHK:  Error - bad call to LINAX for ', A, '-axis.')

      END
C+------------------------------------------------------------------------------
C
      SUBROUTINE LBTEXT (TEXT, LENTXT, LINE1, LINE2, DELTAY, XCORNR, 
     >                   YCORNR, HITE)
C
C ONE LINER:
C     Displays a left-justified text block using DISSPLA utilities. (Level 2-3)
C
C PURPOSE:
C     LBTEXT composes a left-justified text block and displays it at a specified
C     position.  Character height and line spacing are assumed constant within
C     the indicated lines of text.
C
C     LBTEXT was developed as part of the TXTLEG alternative to the original
C     SMDLIB legend utility (TXTBLK).  This is the DISSPLA version. 
C
C ARGUMENTS:
C      ARG       DIM     TYPE    I/O/S     DESCRIPTION
C     TEXT    (*) * (*)   C        I       Character array containing text.
C                                          (Terminating '$'s may be required,
C                                          depending on the usage of 
C                                          LENTXT (*).)  Elements LINE1:LINE2 
C                                          will be displayed in the current 
C                                          color (as opposed to the color 
C                                          associated with each line/symbol).
C     LENTXT     (*)      I        I       Array containing number of characters
C                                          in text line(s), if known.  (No 
C                                          trailing '$'s are needed in this 
C                                          case.)  Otherwise, pass 100 as the 
C                                          first (and only) value, and all lines
C                                          of text will be self counted and all
C                                          lines of text will be self counted
C                                          (requiring the trailing '$'s).
C     LINE1       -       I        I       First line of text to display.
C     LINE2       -       I        I       Last line of text to display.
C     DELTAY      -       R        I       Space between lines in inches.
C     XCORNR      -       R        I       Horizontal position of lower left
C                                          hand corner of text block in inches.
C     YCORNR      -       R        I       Vertical position of lower left hand
C                                          corner of text block in inches.
C     HITE        -       R        I       Text character height in inches.
C
C
C ERROR HANDLING:  No check if LBTEXT is called at the wrong level.
C
C EXTERNAL REFERENCES:  DISSPLA utilities.
C
C ENVIRONMENT:  VAX/VMS; FORTRAN 77
C
C HISTORY:
C     04/17/89    M.D. Wong      Initial design and coding for SMDLIB.
C     06/06/89    M.D. Wong      Adapted from SMDLIB to DISSPLA.  
C     08/29/89    M.D. Wong      Updated to run on DISSPLA version 11.0
C                                (%REF taken out of call to LBTEXT.)
C     02/26/90    M.D. Wong      LENTXT added to argument list to avoid 
C                                reliance on DISSPLA's self counting option.
C     04/12/90    M.D. Wong      Enabled 100 to be passed as first and only
C                                element of the LENTXT array to envoke the
C                                self-counting option.
C
C AUTHOR:  Michael Wong, Sterling Software, Palo Alto, CA.
C
C-------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments
C     ---------

      CHARACTER  TEXT (*) * (*)
      INTEGER    LENTXT (*), LINE1, LINE2
      REAL       DELTAY, HITE, XCORNR, YCORNR

C     Local Variables
C     ---------------

      INTEGER    I, LENGTH
      REAL       GAP, Y 

C     Execution
C     ---------

C     Find distance between lines.

      GAP = DELTAY + HITE 

C     Start with upper left hand corner of text block.

      Y = YCORNR + (LINE2 - LINE1) * GAP

      DO 200 I = LINE1, LINE2
         IF (LENTXT (1) .EQ. 100) THEN
            LENGTH = LENTXT (1)
         ELSE
            LENGTH = LENTXT (I)
         END IF
         CALL MESSAG (TEXT (I), LENGTH, XCORNR, Y) 
         Y = Y - GAP
  200 CONTINUE

      RETURN
      END
C+----------------------------------------------------------------------
C
      SUBROUTINE LCSFIT (NDATA, X, Y, NEW, METHOD, NEVAL, XEVAL, YEVAL,
     &   YPEVAL)
C
C     Two-liner:  Storage-efficient local cubic spline fit (2-space)
C     ----------  (monotonic and piecewise linear options too)
C
C     Description and usage:
C     ----------------------
C
C        LCSFIT is the non-parametric analog of PLSFIT (parametric).
C     It is intended for spline applications which do not require the
C     spline coefficients as output.  It is efficient for repeated
C     calls with the same data, so repeated use with NEVAL = 1 may be
C     preferable to storing vectors of results.
C
C        LCSFIT offers monotonic spline and piecewise linear options
C     also.  And it returns an interpolated first derivative along
C     with the function value.  (The second derivative is omitted
C     because Y" is not guaranteed to be continuous by local methods.)
C
C        See PLSFIT for more details on local methods.  As with most
C     numerical methods, scaling of the data to the unit interval (and
C     unscaling of the result) is recommended to avoid unnecessary
C     effects attributable to the data units.  Utilities GETSCALE and
C     USESCALE from the present authors are appropriate.  The data
C     abscissas should be distinct and either ascending or descending.
C     PROTECT is available to check this.  Extrapolation is permitted
C     (mainly in case of round-off; it is normally inadvisable).
C
C        The CSFIT/CSEVAL or CSDVAL pair are probably preferable if
C     efficiency is not an issue, since CSFIT gives Y" continuity.
C
C     Arguments:
C     ----------
C
C     Name    Type/Dimension  I/O/S  Description
C     NDATA   I               I      Length of X, Y input data arrays.
C
C     X,      R (NDATA)       I      Input data coordinates.  The Xs
C     Y                              must be distinct and monotonic,
C                                    either ascending or descending.
C                                    (No check here.) 
C
C     NEW     L               I      If control flag NEW is .TRUE., the
C                                    search for a bracket starts from
C                                    scratch, otherwise locally-saved
C                                    search and fit information will be
C                                    assumed to be correct. If calling
C                                    LCSFIT from within a loop, set
C                                    NEW = .FALSE. after the first call.
C
C     METHOD   C*1            I      (Uppercase) Type of fit to be used:
C                                    'M' means Monotonic piecewise cubics;
C                                    'B' means non-monotonic "Bessel"-type
C                                        piecewise cubics (looser fit);
C                                    'L' means piecewise Linear fit;
C                                    'C' means Cyclic (periodic) end
C                                        conditions: loose fit assumed.
C
C     NEVAL   I               I      Number of interpolations requested.
C                                    NEVAL >= 1.  One call per result
C                                    (NEVAL = 1) may save storage, and is
C                                    not too inefficient as long as NEW
C                                    is set to .FALSE. after the first.
C
C     XEVAL   R (NEVAL)       I      Abscissa(s) to interpolate to.  These
C                                    are normally in the data range, but
C                                    extrapolation - probably due to
C                                    round-off - is not prevented.
C
C     YEVAL   R (NEVAL)       O      Interpolated function value(s).
C
C     YPEVAL  R (NEVAL)       O      Interpolated 1st derivative value(s).
C                                    Pass the same storage as for YEVAL
C                                    if no derivatives are required.
C
C     Significant local variables:
C     ----------------------------
C
C     MEMORY         Indicates that coefficients are correct for the
C                    current point.
C
C     H, DEL         Delta X and forward difference derivative arrays.
C
C     B, C, D        Coefficients of cubic on the bracketing interval.
C
C     Procedures:
C     -----------
C
C     INTERVAL  1-D "interpolation" search.
C     BESSEL    First derivative (central 3-point formula).
C     BRODLIE   First derivative (central), adjusted for monotonicity.
C     BUTLAND   First derivative (non-central), adjusted for monotonicity.
C     THREEPT   First derivative (non-central 3-point formula).
C
C     Environment:  FORTRAN 90
C     ------------
C
C     Error handling:  None
C     ---------------
C
C     Notes:
C     ------
C
C     (1)  Since many of the calculations must be repeated at both ends
C          of an interval, the various finite difference quantities used
C          are stored as arrays. The following "map" of a typical interior
C          interval and its neighbors should help in understanding the
C          notation.  The local array indices are all numbered relative
C          to the left-hand end of the interval which brackets the point
C          to be evaluated.
C
C                                  LEFT       RIGHT
C
C          Point         -1          0         +1          +2
C
C          Data           x -------- x -------- x --------- x
C
C          Interval      -1          0         +1
C
C
C     Author: Robert Kennelly, Sterling Software/NASA Ames  (PLSFIT)
C     -------
C
C     History:
C     --------
C
C     27 Feb. 1987  R.A.Kennelly  Initial implementation of PLSFIT.
C     23 Aug. 1989  D.A.Saunders  LCSFIT adapted as non-parametric form,
C                                 for embedding in other utilities where
C                                 minimizing work-space is desirable.
C     20 June 1991    "    "      THREEPT (monotonic) renamed BUTLAND;
C                                 THREEPT (pure 3-pt. formula) now used
C                                 for nonmonotonic end-point handling;
C                                 METHOD='C' case belatedly added, as
C                                 needed by PLSINTRP for closed curves.
C     23 July 1991    "    "      The tests for being in the same interval
C                                 as before were not allowing for the
C                                 descending-Xs case.
C     06 May  1998    "    "      Minor Fortran 90 updates.
C     09 Oct  2008    "    "      Revert to FORTRAN 77 for g77 compiler.
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER    NDATA, NEVAL
      REAL       X (NDATA), Y (NDATA), XEVAL (NEVAL)
      REAL       YEVAL (NEVAL), YPEVAL (NEVAL)
      LOGICAL    NEW
      CHARACTER  METHOD * 1

C     Local constants:

      REAL       ZERO, ONE, TWO, THREE
      PARAMETER (ZERO = 0., ONE = 1., TWO = 2., THREE = 3.)

C     Local variables:

      LOGICAL
     &   CYCLIC, LINEAR, MEMORY, MONO
      INTEGER
     &   IEVAL, J, K, LEFT, RIGHT
      REAL
     &   ARROW, B (0:1), C, DELY (-1:1), D, DX, H (-1:1), XBYARROW, XE

C     Procedures:

      REAL      BESSEL, BRODLIE, BUTLAND, THREEPT
      EXTERNAL  BESSEL, BRODLIE, BUTLAND, THREEPT

C     Storage:

      SAVE
     &   ARROW, B, C, D, LEFT, RIGHT

C     Execution:
C     ----------

      MONO   = METHOD == 'M'
      CYCLIC = METHOD == 'C'
      LINEAR = METHOD == 'L'

      IF (CYCLIC) THEN
         IF (Y (NDATA) /= Y (1)) STOP 'LCSFIT: End points must match.'
      END IF

C     Initialize search or avoid it if possible:

      IF (NEW) THEN
         MEMORY = .FALSE.
         ARROW  = SIGN (ONE, X (2) - X (1))
         LEFT   = 1
      END IF

      IEVAL = 1
      XE = XEVAL (1)
      XBYARROW = XE * ARROW

      IF (.NOT. NEW) THEN
      
C        We can save a lot of time when LCSFIT is being called from within
C        a loop by setting MEMORY if possible. The out-of-range checking
C        relies on the fact that RIGHT = LEFT + 1 after the last return.
C        Cater to the more likely case of XE in the previous, interior
C        interval.

         MEMORY = XBYARROW >= X (LEFT)  * ARROW .AND.
     &            XBYARROW <  X (RIGHT) * ARROW

         IF (.NOT. MEMORY) THEN
            MEMORY =
     &         LEFT  == 1     .AND. XBYARROW <  X (RIGHT) * ARROW
     &         .OR.
     &         RIGHT == NDATA .AND. XBYARROW >= X (LEFT)  * ARROW
         END IF

      END IF

      IF (MEMORY) GO TO 70 ! Skip the bulk of the computation

C     Loop over evaluation points requiring a new search:
C     ---------------------------------------------------

   10 CONTINUE

C        Interpolation search for bracketing interval:
C        ---------------------------------------------

         CALL INTERVAL (NDATA, X, XE, ARROW, LEFT)

         RIGHT = LEFT + 1

C         -------------------------------------------
C        |                                           |
C        |   1 <= LEFT < RIGHT = LEFT + 1 <= NDATA   |
C        |                                           |
C         -------------------------------------------

C        Compute derivatives by finite-differences:
C        ------------------------------------------

         IF (NDATA > 2 .AND. .NOT. LINEAR) THEN

C           Interval and derivative approximations:
C           ---------------------------------------

C           The following duplicates more code than PLSFIT's approach,
C           but eliminates some indirection - no need to wrap-around here.
C           Handle the end conditions first to minimize testing LEFT, RIGHT.

            IF (LEFT == 1) THEN

               H (0) = X (2) - X (1)
               DELY (0) = (Y (2) - Y (1)) / H (0)
               H (1) = X (3) - X (2)
               DELY (1) = (Y (3) - Y (2)) / H (1)

               IF (CYCLIC) THEN ! Loose fit assumed
                  H (-1) = X (NDATA) - X (NDATA - 1)
                  DELY (-1) = (Y (NDATA) - Y (NDATA - 1)) / H (-1)
                  B (0) = BESSEL (0, H, DELY)
                  B (1) = BESSEL (1, H, DELY)
               ELSE
                  IF (MONO) THEN
                     B (0) = BUTLAND (0, H, DELY)
                     B (1) = BRODLIE (1, H, DELY)
                  ELSE
                     B (0) = THREEPT (0, H, DELY)
                     B (1) = BESSEL  (1, H, DELY)
                  END IF
               END IF

            ELSE IF (RIGHT == NDATA) THEN

               H(-1) = X (LEFT) - X (LEFT-1)
               DELY(-1) = (Y (LEFT) - Y (LEFT-1)) / H (-1)
               H (0) = X (RIGHT) - X (LEFT)
               DELY (0) = (Y (RIGHT) - Y (LEFT))  / H (0)

               IF (CYCLIC) THEN
                  H (1) = X (2) - X (1)
                  DELY (1) = (Y (2) - Y (1)) / H (1)
                  B (0) = BESSEL (0, H, DELY)
                  B (1) = BESSEL (1, H, DELY)
               ELSE

                  IF (MONO) THEN
                     B (0) = BRODLIE (0, H, DELY)
                     B (1) = BUTLAND (1, H, DELY)
                  ELSE
                     B (0) = BESSEL  (0, H, DELY)
                     B (1) = THREEPT (1, H, DELY)
                  END IF
               END IF

            ELSE

               K = LEFT
               DO J = -1, +1
                  H (J)    =  X (K) - X (K-1)
                  DELY (J) = (Y (K) - Y (K-1)) / H (J)
                  K = K + 1
               END DO

C              Select interpolation scheme:
C              ----------------------------

C              Compute (possibly adjusted) first derivatives at both
C              left- and right-hand endpoints of the interval.

               IF (MONO) THEN

C                 Monotone - use Brodlie modification of Butland's
C                 formula to adjust the derivatives at the knots.

                  B (0) = BRODLIE (0, H, DELY)
                  B (1) = BRODLIE (1, H, DELY)

               ELSE ! METHOD = 'B'

C                 Bessel - use central difference formula at the knots.

                  B (0) = BESSEL (0, H, DELY)
                  B (1) = BESSEL (1, H, DELY)

               END IF

            END IF

C           Compute the remaining cubic coefficients.

            C = (THREE * DELY (0) - TWO * B (0) - B (1)) / H (0)
            D = ( -TWO * DELY (0) + B (0) + B (1)) / H (0) ** 2

         ELSE ! NDATA = 2 .OR. METHOD = 'L'

C           Degenerate case (linear).
C           -------------------------

            B (0) = (Y (RIGHT) - Y (LEFT)) / (X (RIGHT) - X (LEFT))
            C     = ZERO
            D     = ZERO

         END IF

C        Evaluate the cubic (derivative first in case only YEVAL is reqd.):
C        ------------------------------------------------------------------

   70    CONTINUE ! Start of same-interval loop inside new-interval loop

            DX = XE - X (LEFT)
            YPEVAL (IEVAL) = B (0) + DX * (TWO * C + DX * THREE * D)
            YEVAL  (IEVAL) = Y (LEFT) + DX * (B (0) + DX * (C + DX * D))

C           The next evaluation point may be in the same interval:
C           ------------------------------------------------------

            IF (IEVAL < NEVAL) THEN ! Skips this if NEVAL = 1

               IEVAL = IEVAL + 1
               XE = XEVAL (IEVAL)
               XBYARROW = XE * ARROW
               IF (XBYARROW >= X (LEFT)  * ARROW  .AND.
     &             XBYARROW <  X (RIGHT) * ARROW) GO TO 70
            
               GO TO 10 ! Else much more work required.

            END IF

C     Termination:
C     ------------

      RETURN

      END SUBROUTINE LCSFIT
C+------------------------------------------------------------------------------
C
      SUBROUTINE LEGENDZ (TEXT, LENTXT, LINE1, LINE2, DY,
     >                    XCORNR, YCORNR, LEGINF, HITE, SEG, GAP, BOX)
C
C ONE-LINER:  Portable plot legend utility (Level 3)
C
C PURPOSE:
C
C        LEGENDZ composes a plot legend and displays it at a specified position.
C     If the number of legend entries is large and the texts are suitably short,
C     LEGENDZ may be called more than once to position different sections of the
C     legend alongside each other.  (Such positioning is up to the calling
C     program.)  An optional box may be drawn around the legend.
C
C        LEGENDZ was developed as an alternative to existing DISSPLA and SMDLIB
C     schemes for the following reasons:
C
C     (1) rather than underlining the text with the line pattern, the following
C         form is provided:           ---o---  Description
C     (2) CHARACTER-type text is expected instead of an awkward packed-integer
C         array data structure;
C     (3) any implicit connection between curve-drawing and updating the legend
C         is avoided, since this has required work-arounds in the past.
C
C METHOD:
C
C        The symbol/line drawing is handled separately from the writing the
C     text as a left-justified block.  The latter task was modularized when
C     it was found to be a potentially reusable function (subroutine LBTEXT)
C     while the former is done with the standard curve drawing routine,
C     POLYLINE, after appropriate conversion of units from inches to data
C     units via XINVRS and YINVRS.
C
C        Character height and line spacing are assumed constant for each call
C     to LEGENDZ.  The text color is the calling program's current color.
C
C FURTHER USAGE NOTES:
C
C        The (XCORNR, YCORNR) coordinates are in inches relative to the
C     current plot origin.  Carefully-centered (or otherwise positioned)
C     legends almost certainly require something like the following in the
C     calling program:
C
C        LWIDTH = XDIMTB (TEXT, LENTXT, LINE1, LINE2) + SEG + GAP
C        XCORNR = 0.5 * (WIDTH - LWIDTH)
C
C        WARNING:  If the legend is to be drawn outside the plotting area,
C     a suitable grace margin should be defined before the call to LEGENDZ.
C
C ARGUMENTS:
C
C     ARG        DIM     TYPE  I/O/S  DESCRIPTION
C     TEXT    (*) * (*)   C        I  Character array containing legend text.
C                                     (Terminating '$'s may be required,
C                                     depending on the usage of LENTXT (*).)
C                                     Elements LINE1:LINE2 will be displayed in
C                                     the current color (as opposed to the color
C                                     associated with each line/symbol).
C     LENTXT     (*)      I        I  Array containing number of characters
C                                     in text line(s), if known.  (No trailing
C                                     '$'s are needed in this case.)  Otherwise,
C                                     pass 100 as the first (and only) value,
C                                     and all lines of text will be self counted
C                                     (requiring the trailing '$'s).
C     LINE1       -       I      I    First line of legend text to display.
C     LINE2       -       I      I    Last line of legend text to display.
C     DY          -       R      I    Space between lines in inches.
C     XCORNR      -       R      I    Horizontal distance of lower left hand
C                                     corner of legend from origin in inches.
C     YCORNR      -       R      I    Vertical distance of lower left hand
C                                     corner of legend from origin in inches.
C     LEGINF    (3, *)    I      I    Integer array containing codes for
C                                     (1) line type, (2) symbol, and (3) color.
C                                     See POLYLINE for the code definitions.
C     HITE        -       R      I    Text character height in inches.
C     SEG         -       R      I    Length of legend line segment in inches.
C                                     0.54" is a good value: it suppresses the
C                                     last space of several interrupted line
C                                     patterns - see POLYLINE.
C     GAP         -       R      I    Length of gap between legend line segment
C                                     and text in inches. 0.25" is typical.
C     BOX         -       L      I    .TRUE. means draw a box around legend.
C
C ERROR HANDLING:  There is no check for calling LEGEND at the wrong level.
C
C PROCEDURES:
C     LBTEXT    Writes left justified block of text
C     XDIMTB    Finds maximum width of a text block. Used here only if a
C               BOX has been requested.
C     XINVRS    Converts location of a point from inches to data units
C     YINVRS       "              "      "          "         "
C     POLYLINE  Draws a curve using lines and/or symbols
C
C ENVIRONMENT:  VAX/VMS; FORTRAN 77; DISSPLA or SMDLIB
C
C HISTORY:
C     04/17/89    M.D.Wong     Initial implementation for QPLOT translation
C                              from DISSPLA to SMDLIB, using internal COMMONs
C                              and DSDRAW.
C     06/22/89    M.D.Wong     Modularized line/symbol drawing into POLYLINE,
C                              following introduction of conversion utilities
C                              XINVRS and YINVRS.  Eliminated internal COMMONs
C                              at the expense of not checking calling level and
C                              not restoring original line type in order to
C                              combine DISSPLA and SMDLIB versions as one.
C     01/30/90    M.D.Wong     Variable COLOR now passed as REAL to POLYLINE.
C     02/26/90    M.D.Wong     Added LENTXT to argument list to avoid reliance
C                              on DISSPLA dependent self counting option.
C     04/12/90    M.D.Wong     Modified LBTEXT to enable 100 to be passed as
C                              the first and only element of the LENTXT array
C                              to enable self-counting option.  No changes
C                              were necessary to LEGEND.
C     03/27/91   D.A.Saunders  Clarified usage in light of QPLOT revisions.
C     04/07/06     "    "      Changed LEGEND name to LEGENDZ to avoid a
C                              conflict with the DISSPLA version, as pointed
C                              out by Bahman Zohuri of GAE.
C
C AUTHOR:  Michael Wong, Sterling Software, Palo Alto, CA.
C
C-------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments.
C     ----------

      CHARACTER  TEXT (*) * (*)
      INTEGER    LENTXT (*), LEGINF (3, *), LINE1, LINE2
      REAL       DY, HITE, XCORNR, YCORNR
      LOGICAL    BOX

C     Local Constants.
C     ----------------

      REAL       HALF, ZERO
      PARAMETER (HALF = 0.5, ZERO = 0.0)

C     Local Variables.
C     ----------------

      INTEGER    COLOR, I, IER, LINE, SYMBOL
      REAL       GAP, MAR, SEG, X (5), XBLEN, XPOS, Y (5), YBLEN, YPOS

C     Procedures.
C     -----------

      REAL       XDIMTB, XINVRS, YINVRS
      EXTERNAL   XDIMTB, XINVRS, YINVRS, LBTEXT, POLYLINE

C     Execution.
C     ----------


C     Display the legend text as a left-justified block in current color.

      CALL LBTEXT (TEXT, LENTXT, LINE1, LINE2, DY, XCORNR + SEG + GAP,
     >             YCORNR, HITE)


C     Draw a symbol and/or line segment for each legend entry.
C     POLYLINE applies if positions in inches are converted to data units.
C     Any grace margin needed for a legend outside the plot area is left
C     to the application program, since there is no way of restoring the
C     input setting here if it is temporarily adjusted.
C
C     Note that the polar plot case (not to mention general rotated coordinate
C     systems) forces all conversions from inches to stay inside the loop
C     (where some would be constant otherwise).


C     Start with the upper left corner.

      YPOS = YCORNR + (LINE2 - LINE1) * (HITE + DY) + HALF * HITE

      DO 200, I = LINE1, LINE2

         X (1)  = XINVRS (XCORNR, YPOS)
         Y (1)  = YINVRS (XCORNR, YPOS)
         X (2)  = XINVRS (XCORNR + SEG, YPOS)
         Y (2)  = YINVRS (XCORNR + SEG, YPOS)
         LINE   = LEGINF (1, I)
         SYMBOL = LEGINF (2, I)
         COLOR  = LEGINF (3, I)

         IF (LINE .NE. 2) THEN

C           Draw the line segment.

            CALL POLYLINE (2, X, Y, LINE, 0, -1, REAL (COLOR), IER)

         END IF

         IF (SYMBOL .GT. -1) THEN

C           Draw the symbol in the middle of the segment.

            X (1) = XINVRS (XCORNR + HALF * SEG, YPOS)
            Y (1) = YINVRS (XCORNR + HALF * SEG, YPOS)

            CALL POLYLINE (1, X (1), Y (1), 2, 0, SYMBOL, REAL (COLOR),
     >                     IER)

         END IF

         YPOS = YPOS - (HITE + DY)

  200 CONTINUE


      IF (BOX) THEN

C        Pick a reasonable box margin.

         MAR  = .6 * GAP
         XPOS = XCORNR - MAR
         YPOS = YCORNR - MAR
         MAR  = MAR + MAR

C        Calculate the box dimensions (inches).

         XBLEN = XDIMTB (TEXT, LENTXT, LINE1, LINE2) + SEG + GAP + MAR
         YBLEN = REAL (LINE2 - LINE1) * (HITE + DY) + HITE + MAR

C        Draw the box starting from the lower left corner.

         X (1) = XINVRS (XPOS, YPOS)
         Y (1) = YINVRS (XPOS, YPOS)
         X (2) = XINVRS (XPOS + XBLEN, YPOS)
         Y (2) = YINVRS (XPOS + XBLEN, YPOS)
         X (3) = XINVRS (XPOS + XBLEN, YPOS + YBLEN)
         Y (3) = YINVRS (XPOS + XBLEN, YPOS + YBLEN)
         X (4) = XINVRS (XPOS, YPOS + YBLEN)
         Y (4) = YINVRS (XPOS, YPOS + YBLEN)
         X (5) = X (1)
         Y (5) = Y (1)

         CALL POLYLINE (5, X, Y, 3, 0, -1, ZERO, IER)

      END IF

      RETURN
      END
C+----------------------------------------------------------------------
C
      SUBROUTINE LINAX (N, X, MXNSTP, ORIGIN, ENDPT, STEP, ERROR)
C
C     One-liner:  Linear plot axis scale selection, with error checking.
C     ----------
C
C     Description and usage:
C     ----------------------
C
C        Supervises scale selection and error checking for a linear plot
C     axis.  As commonly required by layout constraints, a maximum on the
C     number of scale divisions is imposed.  By default, coordinates
C     increase from ORIGIN to ENDPT, but this is not required.  Any of
C     the inputs ORIGIN, ENDPT and STEP may be specified by the user -
C     the rest will be automatically chosen as follows:  STEP is taken
C     from an array of "nice" values, using as many divisions as possible,
C     but preferring values which divide ORIGIN or ENDPT evenly if they
C     have already been supplied.  The computed values for ORIGIN and
C     ENDPT delimit the smallest range which encompasses the data and is
C     consistent with STEP.  This range will include at least some of the
C     data if possible, but there are cases where no sensible choice is
C     available.  (For example, when the data are all negative, and
C     ORIGIN = 0.0, STEP = 1.0, and ENDPT is free, there is nothing to
C     be done.  Such errors are flagged, as are programming blunders like
C     non-positive array bounds.  See Parameters and Notes for details.
C
C        LINAX, and co-routines EXPRESS, GETSTP and ROUND, were written
C     for use with QPLOT, developed in the Aerodynamics Division of Ames
C     Research Center.  They are based in (small) part on SCALE1 by
C     C. R. Lewart of Bell Labs (see Bibliography).
C
C     Parameters:
C     -----------
C
C     Name     Dimension  Type  I/O/S  Description
C     N                    I    I      Number of data points.
C
C     X           N        R    I      Plot data for one axis.
C
C     MXNSTP               I    I      Maximum number of axis divisions
C                                      permitted.
C
C     ORIGIN               R    I/O    Plot origin (left or lower corner
C                                      of plot on the page).  The axis
C                                      numbering is assumed to begin with
C                                      ORIGIN, which will be calculated
C                                      if a special "flag" value is input.
C                                      See Notes.
C
C     ENDPT                R    I/O    The "other end" of the plot axis,
C                                      usually right or upper corner. It
C                                      will be calculated if a special
C                                     "flag" value is input. See Notes.
C
C     STEP                 R    I/O    Axis labeling increment.  The first
C                                      division will be at ORIGIN + STEP.
C                                      It will be calculated if a special
C                                      "flag" value is input.  If the flag
C                                      is negative (-999), the axis will
C                                      be decreasing.  See Notes.
C
C     ERROR                I      O    Error return flag.
C                                        0:  Normal termination
C
C                                       +1:  Axis length is zero
C                                       +2:  Step length zero or too big
C                                       +3:  Step has wrong sign
C                                            (Bad user input - could just
C                                            reset to defaults and retry)
C
C                                       -1:  Number of data points N
C                                            is non-positive
C                                       -2:  Maximum number of intervals
C                                            MXNSTP is non-positive
C                                            (Bad programmer input -
C                                            abort)
C
C     External references:
C     --------------------
C
C     GETSTP   Axis stepsize selection from an array of candidates.
C     ROUND    Rounds a quotient "up" or "down" in magnitude.
C
C     Environment:  Digital VAX-11/780 VMS FORTRAN (FORTRAN 77)
C     ------------  Apple Macintosh, Absoft MacFORTRAN/020 v2.3
C
C     Notes:
C     ------
C
C     (1)  IMPLICIT NONE is non-standard.
C
C     (2)  A "flag" value input for ORIGIN, ENDPT, or STEP specifies that
C          these be calculated from the data (the default mode).  The
C          present version uses FLAG = 999.0, specified in a PARAMETER
C          statement below.  STEP = -FLAG has special significance: the
C          axis values will decrease from ORIGIN to ENDPT.
C
C     (3)  Error handling:  we trap anything impossible, and many of the
C          unlikely combinations, returning values of ERROR greater than
C          zero for "user" errors (where no sensible choice is available),
C          and less than zero for "programmer" errors, where no action
C          at all is possible - see Parameters section for details.  In
C          the case of user errors, it may suffice to call LINAX again
C          with FLAG in place of the offending quantities.  Two classes
C          of possible "errors" on the part of the user are permitted:
C
C               (a) if an input value of STEP is too small, the number
C                   of scale divisions may exceed MXNSTP
C
C               (b) if all three of the inputs are supplied by the user,
C                   the data may lie outside the range of the plot
C
C          These oddities were intentionally NOT guarded against, as
C          there may be situations where they are intended, and the
C          results are in any event easily recognized as unusual.
C
C          LINCHK, a companion routine available from the author, may
C          be used to interpret the error flag and take recovery action
C          if possible.  It was written for an application where it was
C          necessary to "correct" the user's input and try again.
C
C     (4)  In calculating nice values for both ORIGIN and ENDPT (the
C          last of the four main cases handled), the rounding is tricky:
C          we round ORIGIN "down" (toward zero) when
C
C               (a) the axis is increasing (i.e., STEP > 0)   -and-
C               (b) the initial estimate ORIGIN / STEP > 0
C
C          or when
C
C               (a') STEP < 0            -and-
C               (b') ORIGIN / STEP < 0
C
C          and "up" otherwise so that the data range is covered.  ENDPT
C          is handled analogously, mutatis mutandis. (It took me awhile,
C          too!)
C
C     (5)  DELTA, the rounding grace margin, should be larger than
C          machine epsilon.  It should be small enough that the effect
C          on the plot appearance is not significant.  (The size of the
C          error depends on the plot scale and size, but DELTA should
C          probably be less than the relative error due to hardware
C          precision, <precision, inches> / <axis length, inches>.)  A
C          typical value is 10**-6 < DELTA < 10**-4.  Keep in mind that
C          although this approach preserves the appearance of the data,
C          small "overhanging" quantities may seem to have been lost, e.g.,
C          the value -1.0E-6 won't show up on an axis which also includes
C          data of order +1.0E+2.
C
C     Bibliography:
C     -------------
C
C     (1)  LEWART, C. R.  Algorithms SCALE1, SCALE2, and SCALE3 for
C             Determination of Scales on Computer Generated Plots.
C             CACM, Vol. 16, No. 10, Oct. 1973 (Algorithm 463).
C
C     Author:  Robert Kennelly, Informatics General Corporation
C     -------
C
C     Development history:
C     --------------------
C
C     18 Apr. 1985    RAK    Initial design and coding.
C     15 July 1985    RAK    Use a variable for FIXED (STEP).  Modify test
C                            for when ENDPT = XMIN (or ORIGIN = XMAX) to
C                            take care of STEP < 0 case.  Eliminated
C                            use of FLAG to signal GETSTP that TARGET is
C                            not used (0.0 suffices now).
C     16 July 1985    RAK    Use variables for FIXED (ORIGIN) and (ENDPT).
C     18 July 1985    RAK    When ORIGIN = 0.0, and ENDPT fixed as well,
C                            use ENDPT as step selection target.
C     19 Aug. 1985    RAK    Added 2.5 to NICE array.  Use BIGNEG and
C                            BIGPOS in place of +/- HUGE constant.
C                            Removed "dead" variables MULT and ORDER.
C     21 Aug. 1985    RAK    Added ERROR = +3 output to allow recovery
C                            from a step with wrong sign.
C     23 Aug. 1985    RAK    Must iterate to be sure the final number
C                            of steps does not exceed MXNSTP when
C                            both endpoints are free.
C     15 Oct. 1985    RAK    STEP = -FLAG means decreasing axis.
C     15 Apr. 1987    RAK    Added special handling for MXNSTP = 1 case.
C     24 Apr. 1987    RAK    Changed some multiplies to SIGN comparision
C                            to reduce possibility of overflow.
C      5 Oct. 1987    RAK    Modify search for min and max to eliminate
C                            need for BIGNEG and BIGPOS. Reduce DELTA
C                            to 1.0E-6.
C
C-----------------------------------------------------------------------

C     Declarations.
C     -------------

      IMPLICIT NONE

C     Constants.

      INTEGER
     >   MXNICE
      PARAMETER
     >   (MXNICE = 5)

      REAL
     >   DELTA, FLAG, ZERO, ONE, HALF
      PARAMETER
     >   (DELTA  = 1.0E-6,
     >    FLAG   = 999.0E+0,
     >    ZERO   = 0.0E+0,
     >    ONE    = 1.0E+0,
     >    HALF   = 0.5E+0)

C     Arguments.

      INTEGER
     >   ERROR, MXNSTP, N
      REAL
     >   ENDPT, ORIGIN, STEP, X (N)

C     Local variables.

      LOGICAL
     >   FIXEND, FIXORG, FIXSTP, FLIP, VALID
      INTEGER
     >   I
      REAL
     >   FLOOR, NICE (MXNICE), RAW, REALMX, ROUND, TARGET, UPDOWN,
     >   XMAX, XMIN

C     Storage.

      DATA
     >   NICE /1.0E+0, 2.0E+0, 2.5E+0, 5.0E+0, 10.0E+0/

C     Execution.
C     ----------

C     Check the input parameters for nonsense - note that at most one
C     error will be flagged.

      ERROR = 0

      IF (N .LE. 0)       ERROR = -1
      IF (MXNSTP .LE. 0)  ERROR = -2
      IF (STEP .EQ. ZERO) ERROR = +2

      IF (ERROR .NE. 0)  GO TO 990

C     Initialization.

      REALMX = REAL (MXNSTP)

      FLIP = (STEP .EQ. -FLAG)
      IF (FLIP) STEP = ABS (STEP)

      FIXORG = (ORIGIN .NE. FLAG)
      FIXEND = (ENDPT .NE. FLAG)
      FIXSTP = (STEP .NE. FLAG)

C     Select a case (one of four).
C     ----------------------------

      IF (FIXORG .AND. FIXEND) THEN

C        Both axis limits are fixed.
C        ---------------------------

         RAW = ENDPT - ORIGIN
         IF (RAW .EQ. ZERO) THEN
            ERROR = +1
            GO TO 990
         END IF

         IF (FIXSTP) THEN

C           Error-check the supplied step.  Note that we do NOT check whether
C           this value will yield too many intervals (> MXNSTP).

            RAW = RAW / STEP
            IF (RAW .LT. ZERO) THEN

C              Input STEP must point in the right direction, i.e.,
C              sign (RAW) = sign (STEP), or RAW / STEP > 0.

               ERROR = +3
               GO TO 990
            ELSE IF (RAW .LT. ONE) THEN

C              Input STEP must not exceed the axis length.

               ERROR = +2
               GO TO 990
            END IF
         ELSE

C           Choose STEP.  (We prefer one which divides ORIGIN or ENDPT.)

            TARGET = ORIGIN
            IF (TARGET .EQ. ZERO) TARGET = ENDPT
            STEP = ZERO
            CALL GETSTP (MXNICE, NICE, REALMX, TARGET, RAW, STEP)
         END IF
      ELSE

C        We will be needing one or both of these extrema later.

         XMAX = X (1)
         XMIN = X (1)
         DO 10, I = 1,  N
            XMAX = MAX (X (I), XMAX)
            XMIN = MIN (X (I), XMIN)
   10    CONTINUE

         IF (FIXORG) THEN

C           Choose ENDPT.
C           -------------

C           Estimate ENDPT using maximum, but if we would miss the data
C           entirely or if STEP indicates a reversed axis, choose minimum.

            ENDPT = XMAX
            IF (ORIGIN .GE. XMAX .OR. STEP .LT. ZERO) ENDPT = XMIN

            IF (ENDPT .EQ. ORIGIN) THEN

C              The data are apparently stacked up above the ORIGIN, so
C              just produce a plot which will at least exhibit the data.

               IF (FIXSTP) THEN
                  ENDPT = ORIGIN + STEP
               ELSE
                  ENDPT = ORIGIN + ONE
               END IF
            END IF

            RAW = ENDPT - ORIGIN
            IF (FIXSTP) THEN

C              Is STEP pointing in the right direction?

               IF (SIGN (ONE, RAW) .NE. SIGN (ONE, STEP)) THEN
                  ENDPT = FLAG
                  ERROR = +3
                  GO TO 990
               END IF
            ELSE

C              Choose STEP.

               STEP = ZERO
               CALL GETSTP (MXNICE, NICE, REALMX, ORIGIN, RAW, STEP)
            END IF

C           Pick ENDPT so that there are an even number of intervals
C           of length STEP, beginning with ORIGIN.

            ENDPT = ORIGIN + ROUND (RAW, STEP, DELTA, +ONE) * STEP
         ELSE IF (FIXEND) THEN

C           Choose ORIGIN.
C           --------------

C           Estimate ORIGIN using minimum, but if we would miss the data
C           entirely or if STEP indicates a reversed axis, choose maximum.

            ORIGIN = XMIN
            IF (ENDPT .LE. XMIN .OR. STEP .LT. ZERO) ORIGIN = XMAX

            IF (ENDPT .EQ. ORIGIN) THEN

C              The data are apparently stacked up above the ENDPT, so
C              just produce a plot which will at least exhibit the data.

               IF (FIXSTP) THEN
                  ORIGIN = ENDPT - STEP
               ELSE
                  ORIGIN = ENDPT - ONE
               END IF
            END IF

            RAW = ENDPT - ORIGIN
            IF (FIXSTP) THEN

C              Is STEP pointing in the right direction?

               IF (SIGN (ONE, RAW) .NE. SIGN (ONE, STEP)) THEN
                  ORIGIN = FLAG
                  ERROR  = +3
                  GO TO 990
               END IF
            ELSE

C              Choose STEP.

               STEP = ZERO
               CALL GETSTP (MXNICE, NICE, REALMX, ENDPT, RAW, STEP)
            END IF

C           Pick ORIGIN so that there are an even number of intervals
C           of length STEP, ending with ENDPT.

            ORIGIN = ENDPT - ROUND (RAW, STEP, DELTA, +ONE) * STEP
         ELSE

C           Choose both ORIGIN and ENDPT.
C           -----------------------------

            IF (XMAX .EQ. XMIN) THEN

C              Just produce a plot which will at least exhibit the data.

               IF (FIXSTP) THEN
                  XMAX = XMAX + STEP
               ELSE
                  XMAX = XMIN + ONE
               END IF
            END IF

C           May have to loop until the number of steps is less than MXNSTP
C           since moving both ends can result in violation of the bound.
C           (Relevant only when XMIN, XMAX, and XSTEP are all free.)

            FLOOR = ZERO
   20       CONTINUE
               ORIGIN = XMIN
               ENDPT  = XMAX
               VALID  = .TRUE.

               IF (FIXSTP) THEN

C                 Turn the axis around if necessary.

                  IF (STEP .LT. ZERO) THEN
                     ORIGIN = XMAX
                     ENDPT  = XMIN
                  END IF
               ELSE

C                 Choose STEP.  Note that TARGET = ZERO is used to disable
C                 optional choice of a step which divides ORIGIN or ENDPT.

                  TARGET = ZERO
                  RAW    = ENDPT - ORIGIN
                  STEP   = FLOOR
                  CALL GETSTP (MXNICE, NICE, REALMX, TARGET, RAW, STEP)
               END IF

C              Compute tentative axis endpoints, and set the flag which will
C              force iteration if they are not acceptable.

               IF ((MXNSTP .EQ. 1) .AND.
     >             (.NOT.FIXSTP)   .AND.
     >             (SIGN (ONE, ORIGIN) .NE. SIGN (ONE, ENDPT))) THEN

C                 A special case (alas!) since rounding ORIGIN and ENDPT
C                 both away from zero always yields at least two intervals.
C                 Try arranging the endpoints symmetrically about zero.

                  ORIGIN = SIGN (HALF * STEP, ORIGIN)
                  ENDPT  = -ORIGIN

                  VALID  = (MAX (ABS (XMIN), ABS (XMAX)) .LE.
     >               ABS (ORIGIN))
               ELSE

C                 Pick ORIGIN and ENDPT so that each is an even multiple of
C                 STEP.  (See Notes for discussion of rounding direction.)

                  UPDOWN = -(SIGN (ONE, ORIGIN) * SIGN (ONE, STEP))
                  ORIGIN = ROUND (ORIGIN, STEP, DELTA, UPDOWN) * STEP

                  UPDOWN = SIGN (ONE, ENDPT) * SIGN (ONE, STEP)
                  ENDPT  = ROUND (ENDPT, STEP, DELTA, UPDOWN) * STEP

                  VALID = (FIXSTP) .OR.
     >               (NINT ((ABS (ENDPT - ORIGIN)) / STEP) .LE. MXNSTP)
               END IF

               IF (.NOT.VALID) THEN

C                 The chosen STEP was not large enough (too many steps or
C                 did not cover the data), so increase the minimum stepsize 
C                 and try again.

                  FLOOR = STEP
                  GO TO 20

               END IF

         END IF
      END IF

C     If input STEP was -FLAG, then we want a decreasing axis.

      IF (FLIP .AND. STEP .GT. ZERO) THEN
         XMIN   = MIN (ORIGIN, ENDPT)
         ORIGIN = MAX (ORIGIN, ENDPT)
         ENDPT  = XMIN
         STEP   = -STEP
      END IF

C     Termination.
C     ------------

  990 CONTINUE
      RETURN
      END
C+----------------------------------------------------------------------
C
      SUBROUTINE EXPRESS (X, MULT, EXPON)
C
C     One-liner:  Find "scientific notation" multiplier and exponent.
C     ----------
C
C     Description and usage:
C     ----------------------
C
C        EXPRESS (carefully) breaks up a floating point quantity into
C     a multiplier and exponent, where the multiplier is guaranteed to
C     lie in the range 1.0 <= MULT < 10.0 despite the hazards of finite-
C     precision arithmetic.  It was written to protect GETSTP, an axis
C     scaling routine whose iterative alogorithm requires the multiplier
C     to lie within a known range.
C
C     Parameters:
C     -----------
C
C     Name      Dimension  Type  I/O/S  Description
C     X                     R    I      Floating point number to be
C                                       re-expressed in multiplier
C                                       and exponent notation.
C
C     MULT                  R      O    Multiplier in the range
C                                       [1.0, 10.0); may share
C                                       storage with X.
C
C     EXPON                 I      O    Exponent such that
C                                       X = MULT * (10 ** EXPON).
C
C
C     Environment:  Digital VAX-11/780, VMS FORTRAN (FORTRAN 77)
C     ------------  Apple Macintosh, Absoft MacFORTRAN/020 v2.3
C
C     Notes:
C     ------
C
C        (1)  IMPLICIT NONE and eight character variable names are
C             non-standard.
C
C     Author:  Robert Kennelly, Sterling Federal Systems
C     -------
C
C     Development history:
C     --------------------
C
C      5 Oct. 1987    RAK    Original design and coding.
C
C-----------------------------------------------------------------------

C     Declarations.
C     -------------

      IMPLICIT NONE

C     Constants.

      REAL
     >   ZERO, ONE, TEN
      PARAMETER
     >  (ZERO = 0.0E+0,
     >   ONE  = 1.0E+0,
     >   TEN  = 1.0E+1)

C     Arguments.

      INTEGER
     >   EXPON
      REAL
     >   MULT, X

C     Execution.
C     ----------

      IF (X .EQ. ZERO) THEN

C        Go home early!

         MULT  = ZERO
         EXPON = 0
         GO TO 990
      END IF

C     EXPON is the exponent of the largest power of TEN less than or
C     equal to X. In "scientific" notation:
C
C        X = MULT * TEN ** EXPON, where 1.0 <= MULT < 10.0
C
C     Note that we have to fudge for small X since INT rounds up toward
C     zero when LOG is negative (ABS (X) < 1.0), and we want to round down.

      EXPON = INT (LOG10 (ABS (X)))
      IF (ABS (X) .LT. ONE) EXPON = EXPON - 1

C     Map X into the interval [1.0, 10.0), with some insurance.

      MULT = X * TEN ** (-EXPON)

      IF (MULT .GE. TEN) THEN
         MULT  = MULT / TEN
         EXPON = EXPON + 1
      ELSE IF (MULT .LT. ONE) THEN
         MULT  = MULT * TEN
         EXPON = EXPON - 1
      END IF

C     Termination.
C     ------------

  990 CONTINUE
      RETURN
      END
C+----------------------------------------------------------------------
C
      SUBROUTINE GETSTP (MXNICE, NICE, REALMX, TARGET, RAW, STEP)
C
C     One-liner:  Axis stepsize selection from an array of candidates.
C     ----------
C
C     Description and usage:
C     ----------------------
C
C        Computes a "nice" stepsize for linear plot axes, given an
C     estimate of the axis length, a bound on the number of intervals,
C     and an array of candidates in the range [1.0, 10.0].  The smallest
C     acceptable step is used, i.e. the largest feasible number of axis
C     divisions.  If requested, an attempt is made to choose a step which
C     divides a specifed "target" evenly.  Thus if one end of the aixs is
C     known, the step may be tailored to it (if possible) so that zero
C     will fall on a scale division.
C
C        The value of STEP must be preset to a lower bound on the final
C     stepsize.  This bound will normally be zero, but may be used to force
C     selection of a larger step if a previous call resulted in too many
C     axis intervals after the endpoints were rounded off.
C
C        GETSTP was written to be called by LINAX, but can also be used
C     by itself for simple plot setup.  Error checking is expected to
C     have been performed by the calling routine (see LINAX, for example).
C
C     Parameters:
C     -----------
C
C     Name     Dimension  Type  I/O/S  Description
C     MXNICE               I    I      Number of candidate interval
C                                      lengths, must be >= 1.
C
C     NICE      MXNICE     R    I      Array of interval lengths, in
C                                      the range [1.0, 10.0].  Note that
C                                      10.0 must be included.
C
C     REALMX               R    I      Maximum number of steps permitted.
C                                      Expected to be a whole number,
C                                      greater than zero.
C
C     TARGET               R    I      STEP should divide TARGET evenly
C                                      if possible.  Note that this step
C                                      is skipped (perhaps because not
C                                      required) if 0.0 is input.
C
C     RAW                  R    I      Lower bound on data range.  Must
C                                      be non-zero.
C
C     STEP                 R    I/O    On entry, a lower bound on the
C                                      desired stepsize (final STEP will
C                                      be strictly greater in magnitude
C                                      than the initial value).
C
C                                      The output is a "nice" interval
C                                      scaled so that the entire axis can
C                                      be covered by no more than REALMX
C                                      divisions, and larger in magnitude
C                                      than input STEP.
C
C     External references:
C     --------------------
C
C     EXPRESS  Find "scientific notation" multiplier and exponent.
C
C     Environment:  Digital VAX-11/780 VMS FORTRAN (FORTRAN 77)
C     ------------  Apple Macintosh, Absoft MacFORTRAN/020 v2.3
C
C     Notes:
C     ------
C
C     (1)  IMPLICIT NONE is non-standard.
C
C     (2)  The array of candidates must include 10.0, with additional
C          intermediate points up to the user - a typical array might
C          include 1.0, 2.0, and 5.0 as well.
C
C     (3)  EPS is a machine-dependent constant:  it should be set a
C          little larger than the precision level of floating point
C          calculations.
C
C     (4)  If it is not required that STEP divide TARGET evenly, use zero
C          for TARGET.  In this case, the first acceptable element of
C          NICE will be used.
C
C     (5)  No error checking is performed on REALMX and RAW, which must
C          be non-zero.  This utility is intended for step length
C          calculation after higher levels take care of error detection
C          and recovery.
C
C     Author:  Robert Kennelly, Informatics General Corporation
C     -------
C
C     Development history:
C     --------------------
C
C     23 May  1985    RAK    Initial design and coding.
C      8 July 1985    RAK    Changed test in 10 loop: accept NICE(I)
C                            if equal to STEP.
C     15 July 1985    RAK    Have to pre-scale TARGET to determine if
C                            STEP can divide TARGET evenly.  Test for
C                            TARGET = 0.0 to skip this step (eliminate
C                            use of FLAG).
C     16 July 1985    RAK    Use SCALE to save arithmetic result.
C     23 Aug. 1985    RAK    Added provision for input STEP used as
C                            lower bound (helps LINAX iterate).
C     10 July 1987    RAK    Must enter 10-loop with I = 0 (not 1),
C                            else STEP = NICE (1) is not recognized.
C      6 Oct. 1987    RAK    To avoid trouble due to limited precision,
C                            routine EXPRESS now guarantees that STEP
C                            lies in [1.0, 10.0).  Meaning of SCALE
C                            reversed, plus other minor cleanup.
C
C-----------------------------------------------------------------------

C     Declarations.
C     -------------

      IMPLICIT NONE

C     Constants.

      REAL
     >   EPS, ZERO, ONE, TEN
      PARAMETER
     >   (EPS   = 1.0E-6,
     >    ZERO  = 0.0E+0,
     >    ONE   = 1.0E+0,
     >    TEN   = 1.0E+1)

C     Arguments.

      INTEGER
     >   MXNICE
      REAL
     >   NICE (MXNICE), RAW, REALMX, STEP, TARGET

C     Local variables.

      INTEGER
     >   I, J
      REAL
     >   SCALE, TEST

C     Execution.
C     ----------

C     Lower bound for STEP is based on division of the raw interval length
C     into the maximum number of pieces, or slightly larger than the input
C     value (if non-zero).

      STEP = MAX (ABS (RAW / REALMX), ABS (STEP * (ONE + EPS)))

C     Re-express STEP as a multiplier in the range [1.0, 10.0), with
C     scratch variable I used for the corresponding exponent.

      CALL EXPRESS (STEP, STEP, I)
      SCALE = TEN ** I

C     Find the smallest element of NICE >= the raw step.

      I = 0
   10 CONTINUE
         I = I + 1
         IF (NICE (I) .LT. STEP .AND. I .LT. MXNICE) GO TO 10


C     Look for a larger element of NICE which divides TARGET evenly.

      IF (TARGET .NE. ZERO) THEN

C        The target is pre-scaled for comparison with NICE (cheaper than
C        scaling each element of NICE).

         TEST = TARGET / SCALE
         J    = I - 1
   20    CONTINUE
            IF (J .LT. MXNICE) THEN
               J = J + 1

C              Is the scaled target (nearly) a multiple of NICE?

               IF (ABS (MOD (TEST, NICE (J))) .GT. EPS) THEN
                  GO TO 20

               ELSE

C                 Found one!  Save the result and drop through.

                  I = J
               END IF
            END IF

      END IF

C     Map the best STEP back to the real world.

      STEP = SIGN (NICE (I), RAW) * SCALE

C     Termination.
C     ------------

      RETURN
      END
C+----------------------------------------------------------------------
C
      FUNCTION ROUND (NUM, DENOM, DELTA, UPDOWN)
C
C     One-liner:  Rounds a quotient "up" or "down" in magnitude.
C     ----------
C
C     Description and usage:
C     ----------------------
C
C        Permits modular treatment of rounding a quotient to an integer.
C     This single routine handles symmetric rounding either "up" or
C     "down" in magnitude, with provision for a small grace margin to
C     account for roundoff error or for plotting devices with limited
C     precision.  Input and output is in the form of REAL variables to
C     facilitate floating point computations.
C
C        ROUND was written for use with LINAX, a linear plot axis setup
C     routine, but may find use in other applications. It should improve
C     readability over the usual truncate-and-adjust method in cases
C     where it is necessary to handle different signs.
C
C     Parameters:
C     -----------
C
C     Name     Dimension  Type  I/O/S  Description
C     ROUND                R      O    Quotient NUM / DENOM, rounded as
C                                      specified by flag UPDOWN. (FUNCTION
C                                      value.)
C
C     NUM                  R    I      Numerator of quotient to be rounded.
C
C     DENOM                R    I      Denominator of quotient. Assumed
C                                      to be non-zero.
C
C     DELTA                R    I      Absolute rounding tolerance - if
C                                      the quotient NUM / DENOM is within
C                                      DELTA of an integer, then that
C                                      value is returned even if it rounds
C                                      in the "wrong" direction.
C
C     UPDOWN               R    I      Specifies which direction to round
C                                      the magnitude of the quotient:
C                                       > 0:  "up" (away from zero)
C                                       < 0:  "down" (toward zero)
C                                       = 0:  plain truncation
C
C     Environment:  Digital VAX-11/780 VMS FORTRAN (FORTRAN 77)
C     ------------  Apple Macintosh, Absoft MacFORTRAN/020 v2.3
C
C     Notes:
C     ------
C
C     (1)  IMPLICIT NONE is non-standard.
C
C     (2)  If UPDOWN is zero, the truncated result of (NUM / DENOM) is
C          returned.   This may be useful in some cases, but a direct
C          call to AINT is generally more appropriate.
C
C     (3)  We don't check whether DENOM is zero (a tentative design
C          choice).  In the graphics applications for which ROUND was
C          originally designed, such error checking is better done at
C          a higher level.
C
C     (4)  DELTA may be zero if no grace margin is desired, i.e. strict
C          rounding up or down in magnitude.  It should normally be at
C          least as big as machine epsilon.
C
C     Author:  Robert Kennelly, Informatics General Corporation
C     -------
C
C     Development history:
C     --------------------
C
C      3 June 1985    RAK    Initial design and coding.
C     28 June 1985    RAK    SIGN of rounding correction taken from
C                            RAW, not ROUND.
C      7 Oct. 1987    RAK    Corrected header - ROUND is merely the
C                            truncated quotient when UPDOWN is zero, not
C                            the nearest integer. Other minor prettifying.
C
C-----------------------------------------------------------------------

C     Declarations.
C     -------------

      IMPLICIT NONE

C     Constants.

      REAL
     >   ZERO, ONE
      PARAMETER
     >   (ZERO  = 0.0E+0,
     >    ONE   = 1.0E+0)

C     Arguments.

      REAL
     >   DELTA, DENOM, NUM, RAW, ROUND, UPDOWN

C     Local variables.

      REAL
     >   TEST

C     Execution.
C     ----------

      RAW   = NUM / DENOM
      ROUND = AINT (RAW)

C     Round up or down (unless flag UPDOWN is zero).  If we're very close
C     (2.0 + epsilon, say, when rounding up), chose NEAREST integer (2.0).

      IF (UPDOWN .NE. ZERO) THEN
         TEST = DELTA
         IF (UPDOWN .LT. ZERO) TEST = ONE - DELTA

C        The change in magnitude resulting from the previous "downward"
C        round (truncation) is used to check whether we need to either
C
C           (a) add an additional amount (yielding an "upward" round), or
C           (b) correct a truncated value (if a value rounded "down" was
C               close to an integer).

         IF (ABS (ROUND - RAW) .GT. TEST)
     >      ROUND = ROUND + SIGN (ONE, RAW)
      END IF

C     Termination.
C     ------------

      RETURN
      END
C+----------------------------------------------------------------------
C
      SUBROUTINE LOGAX (LUNOUT, N, X, AXLENG, ORIGIN, ENDPT, CYCLE,
     >   ERROR)
C
C
C     Description and usage:
C
C          This is a skeleton version of a (tentatively) planned log axis
C        scaling routine to be modeled after LINAX.  The present version just
C        does some error checking, with axis scaling by DISSPLA.  We assume
C        that the axis is IN-creasing.
C
C           Note that at present (Oct. 85), we do all the error handling at
C        this level, so the ERROR flag is just decoration for now.  (The
C        output values are set, but need not be acted upon.)
C
C
C     Parameters:
C
C        Name    Dimension  Type  I/O/S  Description
C        LUNOUT              I    I      Unit number for error reporting.
C
C        N                   I    I      Total number of data points.
C
C        X           N       R    I      Packed array of data for one axis.
C
C        AXLENG              R    I      Physical axis length, in inches
C                                        (assumed positive).
C
C        ORIGIN              R    I/O    Requested origin point (left or
C                                        lower end of the axis).  Must be
C                                        greater than zero.  Use flag value
C                                        999.0 to request auto-scaling.
C
C        ENDPT               R    I/O    Requested end point (right or
C                                        upper end of the axis).  Must be
C                                        greater than zero.  Use flag value
C                                        999.0 to request auto-scaling.
C
C        CYCLE               R      O    Log cycle length determined by
C                                        DISSPLA (inches/cycle).
C
C        ERROR               I      O    Trouble flag. >>> Not really used <<<
C
C                                           0:  Normal termination.
C
C                                          +1:  Data range is zero.
C                                          +2:  Input value of ORIGIN negative.
C                                          +3:  Input value of ENDPT negative.
C
C                                          -1:  Number of data points (N)
C                                               is non-positive.  (Bad
C                                               programmer input - abort.)
C                                          -2:  Axis length is not positive.
C                                               (Programming error - abort.)
C
C
C     Environment:  Digital VAX-11/780 VMS Version 4.1 (FORTRAN 77).
C                   ISSCO Graphics DISSPLA Version 9.2
C
C
C     Notes:
C
C        (1)  IMPLICIT NONE is non-standard.
C
C
C     Author:  Robert Kennelly, Informatics General Corporation.
C
C
C     Development history:
C
C         2 Aug. 1985    RAK    Initial design and coding.
C        26 Sep. 1985    RAK    Import error checking and DISSPLA scaling
C                               from old QUICK.  (We don't try to resolve
C                               conflicts - just use DISSPLA's results.)
C         9 Oct. 1985    RAK    Error-check N, ORIGIN, ENDPT, and force the
C                               X array to be > TINY (don't use ABS!).
C                               Report problems on LUNOUT directly - we
C                               handle everything here for now (i.e. the
C                               ERROR output is not really required yet).
C         8 Jan. 1985    RAK    Increased exponents of HUGE and TINY to 38,
C                               essentially the VAX limit for REALs.
C        13 Apr. 1987    RAK    Added internal loop-back to check for >1
C                               input error.  This is NOT the final answer.
C                               Hard STOP if AXLENG is not positive.
C                               Added some protection to the protection
C                               of DISSPLA's log axis scaling routine.
C        16 July 1987    RAK    Re-ordered the test for too-small CYCLE
C                               length so that we test the result from
C                               DISSPLA routine ALGPLT, not just the input.
C
C-----------------------------------------------------------------------


C     Declarations.
C     -------------

      IMPLICIT NONE

C     Constants.

      REAL
     >   FLAG, HUGE, ZERO, TEN, TINY
      PARAMETER
     >   (FLAG = 999.0,                        ! Enables automatic scaling
     >    HUGE = 1.0E+38,                      ! Overflow protection
     >    ZERO = 0.0E+0,
     >    TEN  = 1.0E+1,
     >    TINY = 1.0E-38)                      ! Underflow protection

C     Variables.

      INTEGER
     >   ERROR, I, LUNOUT, N
      REAL
     >   AXLENG, CYCLE, ENDPT, LOGMAX, LOGMIN, ORIGIN, POWER, X (N)


C     Execution.
C     ----------

    1 CONTINUE
      ERROR = 0

      IF (N .LE. 0)         ERROR = -1
      IF (AXLENG .LE. ZERO) ERROR = -2
      IF (ORIGIN .LE. ZERO) ERROR = +2
      IF (ENDPT  .LE. ZERO) ERROR = +3

C     In this (preliminary) version, just report/fix what we can within this
C     routine.  It might be better to separate analysis and error handling
C     eventually, as in LINAX.

C     Negative ERROR's are bugs - the routine needs to assume that N and
C     AXLENG have already been checked.

      IF (ERROR .EQ. -1) THEN
         STOP 'LOGAX:  Fatal array dimension error!'
      ELSE IF (ERROR .EQ. -2) THEN
         STOP 'LOGAX:  Fatal axis length error!'

C     Recoverable errors.

      ELSE IF (ERROR .EQ. +2) THEN
         WRITE (LUNOUT, 1000) 'origin', ORIGIN
         ORIGIN = FLAG
      ELSE IF (ERROR .EQ. +3) THEN
         WRITE (LUNOUT, 1000) 'end point', ENDPT
         ENDPT = FLAG
      END IF

C     Loop back for re-check in case of multiple errors.  (If error checking
C     is ever moved out, this would be accomplished by re-CALLing LOGAX.)

      IF (ERROR .NE. 0) GO TO 1

C     Ensure that the data to be plotted is greater than zero.

      DO 10, I = 1, N
         X (I) = MAX (X (I), TINY)
   10 CONTINUE

C     Determine raw data limits.

      IF (ORIGIN .EQ. FLAG) THEN
         ORIGIN = HUGE
         DO 20, I = 1, N
            ORIGIN = MIN (X (I), ORIGIN)
   20    CONTINUE
      END IF

      IF (ENDPT .EQ. FLAG) THEN
         ENDPT = TINY
         DO 30, I = 1, N
            ENDPT = MAX (X (I), ENDPT)
   30    CONTINUE
      END IF

C     Are the axis limits usable?  Must guarantee positivity, and create
C     a new ORIGIN if data range is zero.

      ORIGIN = MAX (ORIGIN, TINY)
      ENDPT  = MAX (ENDPT, TINY)

      IF (ORIGIN .EQ. ENDPT) THEN
         WRITE (LUNOUT, 1000) 'range', ZERO
         ERROR = +1
         ENDPT = TEN ** (LOG10 (ORIGIN) + 1)
      END IF

C     Let DISSPLA choose the endpoints and stepsize.  ORIGIN and ENDPT are
C     treated as extrema of the data to be accommodated.  For now, ignore any
C     conflicts between ALGPLT-computed and user-input ORIGIN.

      CALL ALGPLT (ORIGIN, ENDPT, AXLENG, ORIGIN, CYCLE)


C     Protect DISSPLA from a too-small cycle length.
C     ----------------------------------------------

C     This nonsense wouldn't be necessary if DISSPLA's log-axis drawing
C     routine were more robust!

      IF (CYCLE .LE. 0.4000) THEN

C        There is not enough room to include all of the data, so cut off the
C        plot at the smallest integral power of ten which will still fit.
C        Protect LOG10 by bounding POWER to avoid underflow.

         LOGMIN = LOG10 (ORIGIN)
         LOGMAX = LOG10 (ENDPT)

         IF (ORIGIN .LT. ENDPT) THEN

C           Normal (increasing) axis.

            POWER  = MAX (TEN ** (LOGMAX - (AXLENG / 0.4001)),
     >         LOG10 (TINY))
            LOGMIN = LOG10 (POWER)
            ORIGIN = TEN ** (INT (LOGMIN) + 1)
            ENDPT  = MAX (ENDPT, TEN * ORIGIN)
         ELSE

C           Inverted.  "MAX" here is really the "lower" limit.

            POWER  = MAX (TEN ** (LOGMIN - (AXLENG / 0.4001)),
     >         LOG10 (TINY))
            LOGMAX = LOG10 (POWER)
            ENDPT  = TEN ** (INT (LOGMAX) + 1)
            ORIGIN = MAX (ORIGIN, TEN * ENDPT)
         END IF

C        Try again for an acceptable CYCLE > .40 inches.

         CALL ALGPLT (ORIGIN, ENDPT, AXLENG, ORIGIN, CYCLE)
      END IF


C     Termination.
C     ------------

      RETURN


C     Formats.
C     --------

 1000 FORMAT (' LOGAX:  Warning - bad log axis ', A, ' = ', E10.3, '.'/
     >   9X, 'Automatic scaling will be used.'//)

      END
C+----------------------------------------------------------------------
C
      SUBROUTINE LOOKUP (NDICT, DICTRY, ALPHA, KEY, ENTRY)
C
C
C     Description and usage:
C
C           Performs dictionary lookups.  A pointer is returned if a
C        match is found between the input key and the corresponding
C        initial characters of one of the elements of the dictionary.
C        If a "synonym" has been provided for an entry, the search is
C        continued until a match to a primary dictionary entry is found.
C        Cases of no match, or multiple matches, are also provided for.
C
C           Dictionary entries must be left-justified, and may be alphabetized
C        for faster searches.  Secondary entries, if any, are composed of
C        two words separated by one or more characters such as blank, tab,
C        comma, colon, or equal sign which are treated as non-significant
C        by SCANNR.  The first entry of each such pair serves as a synonym
C        for the second, more fundamental keyword.
C
C           The ordered search stops after the section of the dictionary
C        having the same first letters as the key has been checked, or
C        after a specified number of entries have been examined.  A special
C        dictionary entry, the vertical bar '|', will also terminate the
C        search.  This will speed things up if an appropriate dictionary
C        length parameter cannot be determined.  Both types of search are
C        sequential.  See "Notes" below for some suggestions if efficiency
C        is an issue.
C
C
C     Arguments:
C
C        Name    Dimension  Type  I/O/S  Description
C        NDICT               I    I      Number of dictionary entries to be
C                                        examined.
C        DICTRY  NDICT       C    I      Array of dictionary entries,
C                                        left-justified in their fields.
C                                        May be alphabetized for efficiency,
C                                        in which case ALPHA should be .TRUE.
C                                        Entries with synonyms are of the form
C                                        'ENTRY:SYNONYM', where 'SYNONYM'
C                                        is a more fundamental entry in the
C                                        same dictionary.  NOTE: Don't build
C                                        "circular" dictionaries!
C        ALPHA               L    I      Indicates whether the dictionary
C                                        is in alphabetical order, in which
C                                        case the search can be terminated
C                                        sooner.
C        KEY                 C    I/O    String to be compared against the
C                                        dictionary.  Abbreviations are legal
C                                        provided they correspond to a unique
C                                        entry in the dictionary.  KEY is
C                                        replaced on termination by its most
C                                        fundamental equivalent dictionary
C                                        entry (uppercase, left-justified) if
C                                        a match was found.
C        ENTRY               I      O    Dictionary pointer.  If > 0, it
C                                        indicates which entry matched KEY.
C                                        In case of trouble, a negative value
C                                        means that a UNIQUE match was not
C                                        found - the absolute value of ENTRY
C                                        points to the second dictionary entry
C                                        which matched KEY.  Zero means that
C                                        NO match could be found.  ENTRY
C                                        always refers to the last search
C                                        performed - in searching a chain of
C                                        synonyms a non-positive value will be
C                                        returned if there is any break, even
C                                        if the original input key was found.
C
C
C     External references:
C
C        Name    Description
C        SCANNR  Finds first and last significant characters.
C
C
C     Environment:  Digital VAX-11/780 VMS FORTRAN (FORTRAN 77).
C
C
C     Notes:
C
C        (1)  IMPLICIT NONE is non-standard.
C
C        (2)  We have assumed that the dictionary is not too big.  If
C             many searches are to be done or if the dictionary has more
C             than a dozen or so entries, it may be advantageous to build
C             an index array of pointers to the beginning of the section
C             of the dictionary containing each letter, then pass in the
C             portion of the dictionary beginning with DICTRY (INDEX).
C             (This won't generally work for dictionaries with synonyms.)
C             For very large problems, a completely different approach may
C             be advisable, e.g. a binary search for ordered dictionaries.
C
C        (3)  LOOKUP is case sensitive.  In most applications it will be
C             necessary to use an uppercase dictionary, and to convert the
C             input key to uppercase before calling LOOKUP.  Companion
C             routines TOKENS and PAIRS, available from the author, already
C             take care of this.
C
C        (4)  The key need not be left-justified.  Any leading (or
C             trailing) characters which are "non-significant" to SCANNR
C             will be ignored.  These include blanks, horizontal tabs,
C             commas, colons, and equal signs.  See SCANNR for details.
C
C        (5)  The ASCII collating sequence for character data is assumed.
C             (Note that this means the numerals precede the alphabet, unlike
C             common practice!)  On some machines, it may be necessary to
C             use the FORTRAN lexical library routines to force use of the
C             ASCII sequence.
C
C        (6)  Parameter NUMSIG sets a limit on the length of significant
C             dictionary entries.  Special applications may require that this
C             be increased.  (It is 16 in the present version.)
C
C        (7)  No protection against "circular" dictionaries is provided: don't
C             claim that A is B, and that B is A.  All synonym chains must
C             terminate!  Other potential errors not checked for include
C             duplicate or mis-ordered entries.
C
C        (8)  The handling of ambiguities introduces some ambiguity:
C
C                ALPHA = .TRUE.  A potential problem, when one entry
C                                looks like an abbreviation for another
C                                (eg. does 'A' match 'A' or 'AB'?) was
C                                resolved by dropping out of the search
C                                immediately when an "exact" match is found.
C
C                ALPHA = .FALSE. The programmer must ensure that the above
C                                situation does not arise: each dictionary
C                                entry must be recognizable, at least when
C                                specified to full length.  Otherwise, the
C                                result of a search will depend on the
C                                order of entries.
C
C
C     Author:  Robert Kennelly, Informatics General Corporation.
C
C
C     Development history:
C
C        24 Feb. 1984  RAK/DAS  Initial design and coding.
C        25 Feb. 1984    RAK    Combined the two searches by suitable
C                               choice of terminator FLAG.
C        28 Feb. 1984    RAK    Optional synonyms in dictionary, no
C                               longer update KEY.
C        29 Mar. 1984    RAK    Put back replacement of KEY by its
C                               corresponding entry.
C        21 June 1984    RAK    Corrected bug in error handling for cases
C                               where no match was found.
C        23 Apr. 1985    RAK    Introduced test for exact matches, which
C                               permits use of dictionary entries which
C                               would appear to be ambiguous (for ordered
C                               case).  Return -I to point to the entry
C                               which appeared ambiguous (had been -1).
C                               Repaired loop termination - had to use
C                               equal length strings or risk quitting too
C                               soon when one entry is an abbreviation
C                               for another.  Eliminated HIT, reduced
C                               NUMSIG to 16.
C        16 May  1988    DAS    Had to use LEN in definition of LAST to
C                               suit revised SCANNR; early termination if
C                               KEY length exceeds LEN (DICTRY (1)).
C
C-----------------------------------------------------------------------


C     Declarations.
C     -------------

      IMPLICIT NONE

C     Arguments.

      INTEGER
     >   ENTRY, NDICT
      LOGICAL
     >   ALPHA
      CHARACTER
     >   DICTRY (NDICT) * (*), KEY * (*)

C     Local constants.

      INTEGER
     >   NUMSIG
      CHARACTER
     >   BLANK, CURLY
      PARAMETER
     >   (BLANK = ' ', CURLY = '{', NUMSIG = 16)

C     Local variables.

      INTEGER
     >   FIRST, I, LAST, LENDIC, LENGTH, MARK
      CHARACTER
     >   FLAG * (NUMSIG), TARGET * (NUMSIG)

C     Procedures.

      EXTERNAL
     >   SCANNR


C     Execution.
C     ----------

      ENTRY = 0

C     Isolate the significant portion of the input key (if any).

      FIRST = 1
      LAST = LEN (KEY)
      CALL SCANNR (KEY, FIRST, LAST, MARK)
      IF (MARK .EQ. 0) GO TO 99

C     Can't hope to find a match if the key is longer than dictionary entries.

      LENGTH = MARK - FIRST + 1
      LENDIC = LEN (DICTRY (1))
      IF (LENGTH .GT. LENDIC) GO TO 99


C     The search starts with the input key, but may be repeated if that
C     target is just a synonym for a more fundamental dictionary entry.
C     NUMSIG = LEN (TARGET) is assumed to be plenty big enough.

      TARGET = KEY (FIRST:MARK)

   10 CONTINUE

C        Select search strategy by cunning choice of termination test
C        flag.  The left curly bracket follows all the alphabetic
C        characters in the ASCII collating sequence, but precedes the
C        vertical bar.

         IF (ALPHA) THEN
            FLAG = TARGET
         ELSE
            FLAG = CURLY
         END IF


C        Perform search.
C        ---------------

         I = 0
   20    CONTINUE
            I = I + 1
            IF (TARGET (1:LENGTH) .EQ. DICTRY (I) (1:LENGTH)) THEN
               IF (ENTRY .EQ. 0) THEN

C                 First "hit" - must still guard against ambiguities
C                 by searching until we've gone beyond the key (ordered
C                 dictionary), or until the end-of-dictionary mark is
C                 reached (exhaustive search).

                  ENTRY = I

C                 Special handling if match is exact - terminate search.
C                 We thus avoid confusion if one dictionary entry looks
C                 like an abbreviation of another.  This fix won't
C                 generally work for un-ordered dictionaries!

                  FIRST = 1
                  LAST = LENDIC
                  CALL SCANNR (DICTRY (ENTRY), FIRST, LAST, MARK)
                  IF (MARK .EQ. LENGTH) I = NDICT
               ELSE


C                 Oops - two hits!  Abnormal termination.
C                 ---------------------------------------

                  ENTRY = -I
                  GO TO 99
               END IF
            END IF

C           Check whether we've gone past the appropriate section of the 
C           dictionary.  The test on the index provides insurance and an
C           optional means for limiting the extent of the search.

            IF (DICTRY (I) (1:LENGTH) .LE. FLAG .AND. I .LT. NDICT)
     >   GO TO 20


C        Check for a synonym.
C        --------------------

         IF (ENTRY .GT. 0) THEN

C           Look for a second entry "behind" the first entry.  (FIRST
C           and MARK were determined above when the hit was detected.)

            FIRST = MARK + 2
            CALL SCANNR (DICTRY (ENTRY), FIRST, LAST, MARK)
            IF (MARK .GT. 0) THEN

C              Reset the target and dictionary pointer and repeat the
C              search for the synonym instead of the original key.

               TARGET = DICTRY (ENTRY) (FIRST:MARK)
               LENGTH = MARK - FIRST + 1
               ENTRY = 0
               GO TO 10

            ELSE

C              Expand the key to the full dictionary entry as a possible aid
C              to the calling program (which may prefer to avoid dealing with
C              entry numbers).

               KEY = DICTRY (ENTRY)

            END IF
         END IF


C     Normal termination.
C     -------------------

   99 RETURN
      END
C+----------------------------------------------------------------------
C
      SUBROUTINE MERGER (NINSERT, RLISTIN, NSORTED, RLISTIO, TOL)
C
C ACRONYM:  MERGE two Real lists (one of which is sorted)
C           -----     -
C PURPOSE:
C        MERGER inserts the (not necessarily ordered) NINSERT real values
C     of RLISTIN (*) into the ordered list RLISTIO (*) such that the
C     output RLISTIO (*) remains ordered, with no duplicates.  Its length,
C     NSORTED, is updated accordingly.
C
C ARGUMENTS:
C     ARG       DIM    TYPE I/O/S DESCRIPTION
C    NINSERT             I  I     Number of values to merge;  NINSERT >= 1
C    RLISTIN  NINSERT    R  I     Values to merge (need not be ordered)
C    NSORTED             I  I/O   Input & output with length of RLISTIO (*)
C                                 NSORTED (I) >= 2;
C                                 NSORTED (O) <= NSORTED (I) + NINSERT
C    RLISTIO NSORTED (I) R  I/O   Input & output list, in either ascending
C            + NINSERT            or descending order, with no duplicates;
C            + 1                  the + 1 allows adding an element temporarily
C                                 to handle INTERVAL's LEFT < NX property
C    TOL                 R  I     Tolerance used for duplicate comparisons;
C                                 TOL >= 0.
C
C METHOD:
C        Search utility INTERVAL finds each point of insertion efficiently;
C     it also handles ordering in either direction.  UTCOPY "makes room" by
C     shifting points to the "right."  Avoiding duplicates is problematic;
C     hence the TOLerance argument in case exact equality is inappropriate.
C
C USAGE:
C        Simple overwriting of the nearest point may seem adequate for the
C     initial application of ensuring that data point abscissas are included
C     exactly when a curve fit is plotted.  But this can fail if data points
C     are very close together.  True insertion (as done here) should find
C     non-graphical applications too.
C
C ERROR HANDLING:  None.  Make sure the sorted array RLISTIO (*) has room.
C
C EXTERNAL REFERENCES:
C    INTERVAL   "Interpolation search" utility
C    UTCOPY     Utility designed for moving array elements "up" or "down"
C
C ENVIRONMENT:  VAX/VMS FORTRAN 77
C
C HISTORY:
C   02/09/89   D.Saunders  Initial implementation.
C
C AUTHOR:  David Saunders, Sterling Software, Palo Alto, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments.

      INTEGER
     >   NINSERT, NSORTED
      REAL
     >   RLISTIN (NINSERT), RLISTIO (NSORTED + NINSERT + 1), TOL

C     Local variables.

      INTEGER
     >   I, LEFT, NSHIFT
      REAL
     >   ARROW

C     Procedures.

      EXTERNAL
     >   INTERVAL, UTCOPY

C     Execution.

      ARROW = SIGN (1.E+0, RLISTIO (2) - RLISTIO (1))

      DO 200, I = 1, NINSERT

C        Since INTERVAL never points to the last point, inserting a
C        temporary extra one saves some logic:

         RLISTIO (NSORTED + 1) = RLISTIO (NSORTED) +
     >      (RLISTIO (NSORTED) - RLISTIO (NSORTED - 1))

C        Locate the interval for insertion:

         CALL INTERVAL (NSORTED + 1, RLISTIO, RLISTIN (I), ARROW, LEFT)

         IF (ABS (RLISTIN (I) - RLISTIO (LEFT)) .GT. TOL) THEN
            NSHIFT = NSORTED - LEFT
            IF (NSHIFT .GT. 0) THEN

C              Make room by shifting elements to the "right":

               NSHIFT = -NSHIFT

               CALL UTCOPY (NSHIFT, RLISTIO (LEFT + 1),
     >            RLISTIO (LEFT + 2))
            END IF

            RLISTIO (LEFT + 1) = RLISTIN (I)
            NSORTED = NSORTED + 1
         END IF

  200 CONTINUE

      RETURN
      END
C+----------------------------------------------------------------------
C
      SUBROUTINE MYSTYL (STYLE)
C
C
C     Description and usage:
C
C           This is boiler-plate code pulled out of QUICK to reduce its
C        already imposing bulk.  Although modeled on MYSPEC, which provides
C        "memory" for line types and colors, it need not have been coded
C        as a separate module.  As in the rest of QUICK, default selection
C        is assumed to have been handled at a higher level.
C
C           Several of the styles are part of the "shaded characters" option
C        which is not part of basic DISSPLA, and hence may not be available
C        on all systems.
C
C           MYSTYL was written for use in QPLOT, an easy-to-use scientific
C        plotting program, but may be useful in other DISSPLA-based
C        applications.
C
C
C     Arguments:
C
C        Name    Dimension  Type  I/O/S  Description
C        STYLE               C    I      Name of one of DISSPLA's lettering
C                                        styles.  Assumed to have been
C                                        previously checked: must be left-
C                                        justified, uppercase, with no
C                                        abbreviations or synonyms. See
C                                        below for legal values. NOTE: to
C                                        reset to default, use 'DISSPLA'.
C
C
C     Environment:  Digital VAX-11/780 VMS FORTRAN.
C
C
C     Notes:
C
C        (1)  IMPLICIT NONE is non-standard.
C
C
C     Author:  Robert Kennelly, NASA-Ames/Sterling Software.
C
C
C     Development history:
C
C        22 Oct. 1986    RAK    Initial design and coding.
C         8 Sep. 1987    RAK    Permit 'DISSPLA' as input, signifying reset.
C
C-----------------------------------------------------------------------


C     Declarations.
C     -------------

      IMPLICIT NONE

C     Arguments.

      CHARACTER
     >   STYLE * (*)


C     Execution.
C     ----------

      IF (STYLE .EQ. 'DISSPLA') THEN

C        Reset to DISSPLA's default style (simple stick figures).

         CALL RESET ('CARTOG')
      ELSE IF (STYLE .EQ. 'CARTOG') THEN
         CALL CARTOG
      ELSE IF (STYLE .EQ. 'SIMPLX') THEN
         CALL SIMPLX
      ELSE IF (STYLE .EQ. 'SCMPLX') THEN
         CALL SCMPLX
      ELSE IF (STYLE .EQ. 'COMPLX') THEN
         CALL COMPLX
      ELSE IF (STYLE .EQ. 'DUPLX') THEN
         CALL DUPLX
      ELSE IF (STYLE .EQ. 'TRIPLX') THEN
         CALL TRIPLX
      ELSE IF (STYLE .EQ. 'GOTHIC') THEN
         CALL GOTHIC


C     The remaining styles are optional "shaded character" types.
C     -----------------------------------------------------------

      ELSE IF (STYLE .EQ. 'FUTURA') THEN
         CALL FUTURA
      ELSE IF (STYLE .EQ. 'SERIF') THEN
         CALL SERIF
      ELSE IF (STYLE .EQ. 'FASHON') THEN
         CALL FASHON
      ELSE IF (STYLE .EQ. 'LOGO1') THEN
         CALL LOGO1
      ELSE IF (STYLE .EQ. 'SWISSL') THEN
         CALL SWISSL
      ELSE IF (STYLE .EQ. 'SWISSM') THEN
         CALL SWISSM
      ELSE IF (STYLE .EQ. 'SWISSB') THEN
         CALL SWISSB
      ELSE

C        Coding error - should be FIXED, not just recovered-from.

         STOP 'MYSTYL:  ERROR - unrecognized STYLE!'
      END IF


C     Termination.
C     ------------

      RETURN
      END
C+----------------------------------------------------------------------
C
      SUBROUTINE NLCHEK (LINE, LAST, NAME, PRESENT, FLAKEY)
C
C
C     One-liner:  Identifies legal or suspicious namelists.
C
C
C     Description:
C
C           This utility checks for the start of a namelist in the given
C        string.  The following simple errors are checked-for to trap
C        some bad input:
C
C           (a)  "$" or "&" in column 2, followed by something other than
C                <name>,
C
C           (b)  line beginning with <name> (possibly preceded by "$"
C                or "&", but not in column 2).
C
C        (We don't try to identify namelists which are mis-spelled AND
C        mis-positioned.)
C
C           NLCHEK was originally written for QPLOT but may find use in
C        other contexts where a flexible input format is desired.
C
C           This version does NOT read data (unlike earlier versions).
C
C
C     Parameters:
C
C        Name    Dimension  Type  I/O/S  Description
C        LINE                C    I      Line to be checked. Its case is
C                                        insignificant.
C        LAST                I    I      Saves a LEN (LINE).
C        NAME      *         C    I      The namelist name sought (without
C                                        the namelist symbol, upper case).
C        PRESENT             L      O    Indicates that the namelist is
C                                        actually present.
C        FLAKEY              L      O    Something like a namelist is
C                                        present, but doesn't have correct
C                                        form.
C
C
C     External References:
C
C        UPCASE   Converts a string to uppercase.
C
C
C     Environment:  Digital VAX-11/780 VMS FORTRAN (FORTRAN 77)
C                   SGI IRIS 4D, IRIX 3.3, f77
C                   Cray Y-MP, UNICOS, cft77
C
C     Notes:
C
C        (1)  IMPLICIT NONE and eight character variable names are not (yet)
C             standard FORTRAN.  Use of CHAR function in a PARAMETER statement
C             is not standard and has been removed from this version.
C
C        (2)  NAMELIST is non-standard.  This version checks for the two
C             VAX-legal NAMELIST symbols $ and &, which must be in column 2.
C             The NAMELIST name must be followed by a blank or a tab. Other
C             implementations may differ.
C
C
C     History:
C
C        10 June 1983    RAK    Initial design and coding.
C         3 Jan. 1984    RAK    Skips blank lines, handles positioning of
C                               input file.
C         6 Feb. 1984    RAK    Added EOF output flag.
C         9 Oct. 1987    RAK    Use GETLINE to filter out blank and comment
C                               lines (flagged by exclamation mark). Revised
C                               error categories. TEST lengthened to
C                               permit longer names. Hard STOPs for blank
C                               or too-long namelist name (debugging).
C        28 Dec. 1988    RAK    Permit the NAMELIST name to be followed by
C                               a tab, since the VAX doesn't mind.
C        10 June 1991    DAS    CHAR (9) moved from PARAMETER statement
C                               for Cray purposes.
C        16 Nov. 1991    DAS    LINE is an input now - no reading.
C
C     Author:  Robert Kennelly, Informatics General Corporation.
C
C-----------------------------------------------------------------------


C     Declarations.
C     -------------

      IMPLICIT NONE

C     Arguments.

      INTEGER
     >   LAST
      CHARACTER
     >   LINE * (*), NAME * (*)
      LOGICAL
     >   FLAKEY, PRESENT

C     Local constants.

      CHARACTER
     >   AMPERS, BLANK, DOLLAR, EXCLAM
      PARAMETER
     >  (AMPERS = '&',
     >   BLANK  = ' ', 
     >   DOLLAR = '$',
     >   EXCLAM = '!')

C     Local variables.

      INTEGER
     >   FIRST, IOS, LENGTH, MARK
      CHARACTER
     >   HT * 1, TEST * 8


C     Execution.
C     ----------

      HT = CHAR (9)
      PRESENT = .FALSE.
      FLAKEY = .FALSE.

      IF (LINE (1:2) .EQ. BLANK // DOLLAR .OR.
     >    LINE (1:2) .EQ. BLANK // AMPERS) THEN

C        A namelist symbol is in the right place. Look for the correct
C        name followed by a blank or a tab.

         LENGTH = LEN (NAME)
         TEST = LINE (3:LENGTH + 2)
         CALL UPCASE (TEST)

         IF (TEST (1:LENGTH) .EQ. NAME .AND.
     >      (LINE (LENGTH + 3:LENGTH + 3) .EQ. BLANK .OR.
     >       LINE (LENGTH + 3:LENGTH + 3) .EQ. HT)) THEN

C           The namelist begins properly.

            PRESENT = .TRUE.
         ELSE

C           Looks like a garbled namelist - the symbol is in the correct
C           place but it's not followed by the requested name.

            FLAKEY = .TRUE.
         END IF
      ELSE

C        Check for the name, with or without the symbol, as first token
C        on the line. We strip off leading $ or &, if present, since the
C        subsequent test on the rest of the token is all that counts here.

         FIRST = 1
         CALL SCANNR (LINE, FIRST, LAST, MARK)

         IF (LINE (FIRST:FIRST) .EQ. DOLLAR .OR.
     >       LINE (FIRST:FIRST) .EQ. AMPERS) FIRST = FIRST + 1

         TEST = LINE (FIRST:MARK)
         CALL UPCASE (TEST)
         IF (TEST (1:LENGTH) .EQ. NAME) FLAKEY = .TRUE.
      END IF


C     Termination.
C     ------------

      RETURN
      END
C+----------------------------------------------------------------------
C
      SUBROUTINE NLGRIPE (LUN, NAME)
C
C
C     Description and usage:
C
C        This is a not-entirely-serious modular lecture on proper namelist
C     format for naive VAX users. The idea is to be able to call up a block
C     of expanatory text on a standard topic without cluttering the calling
C     routine (especially since use of namelists is on the way out). Thus it
C     should be easy to remove NLGRIPE when appropriate.
C
C        NLGRIPE was written for use with QPLOT, which for (somewhat dubious)
C     historical reasons supports a namelist for plotting options, causing
C     much grief to neophytes who've never seen a namelist before.
C
C
C     Arguments:
C
C        Name    Dimension  Type  I/O/S  Description
C        LUN                 I    I      Logical unit number to receive the
C                                        lecture on namelist format. Must be
C                                        > 0, or the (non-fatal) message will
C                                        be skipped (for compatibility with
C                                        other message schemes).
C                                        
C        NAME      *         C    I      The namelist name to be used as
C                                        an example.
C
C
C     External files:
C
C        Unit    I/O/S  Description
C        LUN       O    Receives the text of the sermon.
C
C
C     Environment:  Digital VAX-11/780 VMS FORTRAN
C
C
C     Notes:
C
C        (1)  IMPLICIT NONE and eight character variable names are not (yet)
C             standard FORTRAN.
C
C
C     Author:  Robert Kennelly, NASA-Ames/Sterling Federal Systems
C
C
C     Development history:
C
C         9 Oct. 1987    RAK    Initial design and coding.
C
C-----------------------------------------------------------------------


C     Declarations.
C     -------------

      IMPLICIT NONE

C     Arguments.

      INTEGER
     >   LUN
      CHARACTER
     >   NAME * (*)


C     Execution.
C     ----------

      IF (LUN .LT. 0) GO TO 990

      WRITE (LUN, 1000)
     >   ' ',
     >   '          Summary of namelist format rules',
     >   ' ',
     >   'A FORTRAN namelist is introduced by a special symbol ($)',
     >   'in the second column of the input file, followed by a name'
      WRITE (LUN, 1010)
     >   '(in this case, ', NAME, '), in upper or lowercase. Data items'
      WRITE (LUN, 1000)
     >   'follow in the form KEYWORD = VALUE, with separating commas.',
     >   'Except for the first line, the text may begin in any column',
     >   'beyond the first (the first is ignored), and is terminated',
     >   'by the special string $END. The various keywords may be in',
     >   'any order, and need not all be present. Misspelled keywords',
     >   'are a common problem. Note that CHARACTER type values must',
     >   'be enclosed in single quote marks.',
     >   ' '
      WRITE (LUN, 1000)
     >   'A simple example follows. Variable N is integer, STRING is',
     >   'CHARACTER type. Note the $ in column two!',
     >   ' '
      WRITE (LUN, 1010)
     >   ' $', NAME, ' N = 123,'
      WRITE (LUN, 1000)
     >   ' STRING = ''Some text, enclosed by single quotes'',',
     >   ' $END',
     >   ' ',
     >   '^ Column one is blank.',
     >   ' '


C     Termination.
C     ------------

  990 CONTINUE
      RETURN


C     Formats.
C     --------

 1000 FORMAT (9X, A)
 1010 FORMAT (9X, A, A, A)

      END
C+----------------------------------------------------------------------
C
      FUNCTION NUMBER( STRING )
C
C     Acronym:
C
C        if ( NUMBER( string ) ) then
C           string very likely represents a number - integer or real
C
C     Description and usage:
C
C        A simple(-minded) test for numeric data is implemented by
C        searching the input string for legitimate characters:
C                digits 0 to 9, D, E, -, + and .
C        Insurance is provided by requiring that a numeric string
C        have at least one digit, at most one D, E or .
C        and at most two -s or +s.  Note that a few ambiguities remain:
C
C           (a)  A string might have the form of numeric data but be
C                intended as text.  No general test can hope to detect
C                such cases.
C
C           (b)  There is no check for correctness of the data format.
C                For example a meaningless string such as 'E1.+2-'
C                will be accepted as numeric.
C
C        Despite these weaknesses, the method should work in the
C        majority of cases.
C
C     Arguments:
C
C        Name    Dimension  Type  I/O/S  Description
C        NUMBER      -       L      O    Set .TRUE. if STRING appears
C                                        to be numerical data, else
C                                        set .FALSE.
C        STRING      *       C    I      Input data to be tested,
C                                        assumed to be in upper case.
C
C     Notes:
C
C        (1)  It is assumed that STRING has been extracted by
C             a "token" utility - hence the upper case assumption.
C
C        (2)  The scan of STRING stops at the first blank.
C
C        (3)  COMPLEX data with parentheses will not look numeric.
C
C     Environment:  ANSI FORTRAN 77.
C
C     Michael Saunders, Systems Optimization Lab., Stanford University.
C     12 Nov  1985    Initial design and coding, starting from the
C                     routine ALPHA from Informatics General, Inc.
C     23 May  1986    OPNUMB name changed to NUMBER, analogous to ALPHA,
C                     for use at NASA Ames - D. Saunders, Informatics.
C
C-----------------------------------------------------------------------

C  Arguments:

      LOGICAL          NUMBER
      CHARACTER*(*)    STRING

C  Local variables:

      LOGICAL         NUM
      INTEGER         J, LENGTH, NDIGIT, NEXP, NMINUS, NPLUS, NPOINT
      CHARACTER*1     ATOM

C  Executable statements:

      NDIGIT = 0
      NEXP   = 0
      NMINUS = 0
      NPLUS  = 0
      NPOINT = 0
      NUM    = .TRUE.
      LENGTH = LEN (STRING)
      J      = 0

   10    J    = J + 1
         ATOM = STRING (J:J)
         IF      (ATOM .GE. '0'  .AND.  ATOM .LE. '9') THEN
            NDIGIT = NDIGIT + 1
         ELSE IF (ATOM .EQ. 'D'  .OR.   ATOM .EQ. 'E') THEN
            NEXP   = NEXP   + 1
         ELSE IF (ATOM .EQ. '-') THEN
            NMINUS = NMINUS + 1
         ELSE IF (ATOM .EQ. '+') THEN
            NPLUS  = NPLUS  + 1
         ELSE IF (ATOM .EQ. '.') THEN
            NPOINT = NPOINT + 1
         ELSE IF (ATOM .EQ. ' ') THEN
            J      = LENGTH
         ELSE
            NUM    = .FALSE.
         END IF

         IF (NUM  .AND.  J .LT. LENGTH) GO TO 10

      NUMBER = NUM
     $         .AND.  NDIGIT .GE. 1
     $         .AND.  NEXP   .LE. 1
     $         .AND.  NMINUS .LE. 2
     $         .AND.  NPLUS  .LE. 2
     $         .AND.  NPOINT .LE. 1

      RETURN
      END
C+----------------------------------------------------------------------
C
      SUBROUTINE OPENER (LUNCRT, PROMPT, LUNKBD, FNAME, LUNFIL, FSTAT)
C
C
C  PURPOSE:
C
C        OPENER modularizes the common situation of prompting for a file
C     name, opening the file, and reprompting if the file is not found.
C     Isolating any system dependencies here (such as VAX/VMS's READONLY
C     extension) enhances the portability of typical applications.
C
C        This version is restricted to sequential files, either formatted
C     or unformatted, old or new, with a few of the other occasionally
C     desirable attributes (e.g. CARRIAGECONTROL='LIST') provided for too.
C
C        This version also permits proceeding anyway if a file specified
C     as old is not found.  This option required FSTAT to be used as an
C     OUTPUT as well as an input; all other options use it as input only.
C
C        System-dependent feature:
C
C        If there is trouble opening the file, the program user has the
C     option to execute an operating-system command - probably to look in
C     some directory to see what the filename should be - then try again.
C
C
C  ARGUMENTS:
C
C   NAME     DIM    TYPE I/O/S DESCRIPTION
C  LUNCRT     -      I     I   Logical unit for screen prompts.
C  PROMPT     *      C     I   Prompt (possibly indicating default filename).
C                              May be blank to suppress the prompt (but FNAME
C                              must be non-blank in this case).
C  LUNKBD     -      I     I   Logical unit for keyboard entries.
C  FNAME      *      C    I/O  Name of file to be opened - see METHOD/NOTES.
C                              May be blank to permit termination of an
C                              indefinite loop over file names (carriage
C                              return response to prompt, no open attempted,
C                              and a check in the calling program to see if
C                              FNAME is still blank).  A non-blank FNAME on
C                              input will be treated as the default file name.
C  LUNFIL     -      I     I   Logical unit for file to be opened.
C  FSTAT      *      C   I[/O] INPUT:  String of file status attributes,
C                              separated by commas, colons, or blanks.
C                              Unique abbreviations suffice, in upper or
C                              lower case.  The more common possibilities
C                              provided for appear below.
C                              OUTPUT:  None, with one exception: if the
C                              keyword 'IfPresent' is included in the input 
C                              FSTAT string, then FSTAT = 'MISSING' (upper
C                              case) is returned as output if the file is
C                              not found.  ('Old' is implied by 'IfPresent'
C                              here.)  FSTAT must be a character VARIABLE
C                              in this one case; a CONSTANT is fine in all
C                              other cases.
C                              
C  FSTAT examples:
C
C     'NEW, LIST, 160'      New formatted file with up to 160 characters per
C                           record, and implied (not explicit) carriage control
C
C     'old:binary:write'    Existing unformatted file where READONLY access (the
C                           default) is not enough
C
C     'IfPresent'           Suppresses reprompting if the specified file is
C                           not found.
C
C
C  FSTAT token summary:
C
C     (Meaningful combinations should be concatenated.)
C
C  Attribute        Description          Corresponding OPEN keyword   Default
C
C     'OLD'        File already exists            STATUS              'UNKNOWN'
C     'NEW'        New file desired                "  "                 "   "
C     'SCRATCH'    New file; deleted upon closing  "  "                 "   "
C
C     'BINARY'     Unformatted file                FORM              'FORMATTED'
C     'UNFORMATTED'     "      "                   "  "                 "   "
C
C     'LIST'       Single-spaced records   CARRIAGECONTROL 'FORTRAN' (formatted)
C     'NONE'       No implied carriage ctrl     "  "  "    Default for unfrmttd.
C
C     'WRITE'      Write (and read) access        READONLY           The default
C                  is 'READONLY' if STATUS='OLD', else it is 'WRITE'=read/write.
C
C     'nnn'        Max. record length             RECL         <System default> 
C
C     'IfPresent'  See FSTAT description and examples.
C
C  Further FSTAT Notes:
C
C     Integers nnn for RECL refer to bytes (characters) for formatted files
C     (where 132 is normally the upper limit), or longwords (4-byte units)
C     for unformatted files (where the default is system-dependent).
C
C     READONLY is normally needed for reading files of another owner.
C
C     The defaults are also legitimate (if redundant) input values.
C
C
C  METHOD:
C
C     (1) Default the file attributes.  (The READONLY one will be adjusted
C         later if the file status is new or unknown.)
C
C     (2) Decode the subset of attributes indicated in FSTAT, one token at
C         a time - just itemize the finite number of cases in a way that
C         could be extended.  Any unknown attributes are considered fatal
C         programmer errors - stop with a diagnostic.
C
C     (3) Prompt for the name of the file to be opened (unless the prompt
C         is blank - handy for opening files with fixed names).
C
C     (4) If a carriage return is entered, and FNAME was passed to OPENER
C         as all blanks, return immediately without opening a file.  The
C         calling program can detect this case by checking if FNAME is
C         still blank - handy for performing an indefinite loop over
C         multiple files.
C
C     (5) If a carriage return is entered, and FNAME was NOT all blanks,
C         then FNAME is assumed to represent a default filename (presumably
C         indicated as part of the prompt).
C
C     (6) If "EOF" is entered (^Z under VMS; ^D under Unix), this is assumed
C         to mean "STOP."  Perhaps this should be indicated to the calling
C         program via a returned argument value.  However, the original
C         design had the stop occurring here, and there is no good, upwardly
C         compatible way of changing that.
C
C     (7) If the keyword 'IfPresent' has been input in FSTAT, use INQUIRE to
C         detect existence, and return with FSTAT = 'MISSING' if appropriate.
C
C     (8) Attempt to open the file:
C
C         IF the file cannot be opened THEN
C            Inform the user.
C            Prompt for either the correct file or an operating
C               system command (probably for a directory listing).
C            IF the first character found is a '$' THEN
C               Execute the associated command (system-dependent)
C               Go back to (3)
C            ELSE
C               Go back to (8)
C            END IF
C         END IF
C
C  EXTERNAL REFERENCES:
C
C     LOOKUP        Dictionary lookup: permits abbreviations
C     NUMBER        Identifies RECL parameter
C     READER        Prompting utility
C     SCANNR        Identifies tokens in FSTAT string
C     SYSTEM        IRIS utility analogous to VMS's LIB$SPAWN
C     UPCASE        Upper-casifies FSTAT prior to SCANNR/LOOKUP
C
C  ERROR HANDLING:
C
C        Normally, if a file cannot be opened, a message to that effect is
C     sent to the screen.  The user then has three options:  re-enter the
C     file, which is consequently opened; type a command to list directory
C     contents, whereupon the prompt/response cycle is reinitiated; or quit.
C
C        Alternatively, a missing file may be handled by the calling program.
C     See the FSTAT description for more.
C
C        File attributes found in FSTAT (e.g. 'NULL' for the BLANK keyword)
C     which are not (yet) handled by OPENER are fatal: OPENER will STOP.
C
C  NOTES:
C
C        The declared length of the name of the file to be opened should be
C     enough to accommodate a possibly lengthy directory specification.  A
C     generous size, such as CHARACTER * 60, is suggested.
C
C  SYSTEM DEPENDENCIES:
C
C     (1) IMPLICIT NONE is nonstandard.
C     (2) READONLY and CARRIAGECONTROL keywords are VAX/VMS extensions.
C     (3) CALL SYSTEM is IRIS-specific.
C
C  ENVIRONMENT:  IRIS/IRIX, FORTRAN 77
C
C  HISTORY:
C
C  02/27/86  R.G.Langhi  Initial design and code (formatted files only).
C
C  05/24/86  D.Saunders  Provided for defaulting filename, and for
C                        returning with FNAME=' ' and no file opened
C                        to handle the indefinite-loop-over-files case;
C                        documentation clarified (sequential/formatted).
C
C  08/28/86  RGL         Added unformatted sequential capability.
C
C  05/28/87  RGL         Bug:  the default filename was lost after an
C                        erroneous response to the prompt; need to keep
C                        a copy locally.  On cancellation, FNAME is
C                        returned now with the default instead of blank.
C                        (Note:  The default filename MAY be blank.)
C
C  05/10/88  DAS         Generalized FSTAT usage to indicate more than
C                        just old/new and formatted/unformatted, giving
C                        CARRIAGECONTROL, READONLY, and RECL control too;
C                        suppressed prompt if it is blank, so that OPENER
C                        can still be used with fixed file names.
C
C  02/15/90  DAS         Generalized FSTAT further to permit the calling
C                        program to proceed without the specified (old)
C                        file (as opposed to OPENER's reprompting, which
C                        remains the usual option).
C                        Also: changed from READC to READS to accommodate
C                        the case-sensitive Unix community.
C
C  02/22/90  DAS         IRIS 4D version: suppressed LIB$SPAWN feature
C                        (is there an IRIX equivalent?); ^D instead of ^Z
C                        indicated in the prompts.
C
C  03/14/90  DLH/DAS     Dexter found "call system" as a substitute for
C                        LIB$SPAWN.
C
C  08/28/90  DAS         Unformatted files needed CC='NONE'.
C
C  06/10/91  DAS         NUMBER was declared INTEGER - should be LOGICAL.
C
C  07/29/02  DAS         READONLY and CARRIAGECONTROL are not standard F90.
C
C  AUTHOR:  Ronald Langhi, NASA Ames/Sterling Software, Palo Alto, CA
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER
     >   LUNCRT, LUNKBD, LUNFIL
      CHARACTER
     >   FNAME * (*), FSTAT * (*), PROMPT * (*)

C     Local constants:

      INTEGER
     >   MXCHAR, MXWORD
      CHARACTER
     >   BLANK * 1

      PARAMETER
     >  (BLANK  = ' ',
     >   MXCHAR = 11,
     >   MXWORD = 13)

C     Local variables:

      INTEGER
     >   ENTRY, FIRST, LAST, MARK, RECLEN
      CHARACTER
     >   ATTRIB * 80, CC * 7, DICT (MXWORD) * (MXCHAR), FORIG * 80,
     >   FORM * 11, KEY * (MXCHAR), STATUS * 7
      LOGICAL
     >   CR, DFRECL, ENQUIRE, EOF, PRESENT, RDONLY

C     Procedures:

      LOGICAL
     >   NUMBER
      EXTERNAL
     >   LOOKUP, NUMBER, SCANNR, READS, SYSTEM, UPCASE

C     Storage:

      DATA
     >   DICT
     >      /'BINARY', 'FORMATTED', 'FORTRAN', 'LIST', 'NEW', 'NONE',
     >       'OLD', 'READONLY', 'SCRATCH', 'UNFORMATTED', 'UNKNOWN',
     >       'WRITE', 'IFPRESENT'/
C     The dictionary should be in upper case.  It need not be alphabetized.


C     Execution:

C     Save the default input filename in case of an erroneous response 
C     to the prompt:

      FORIG = FNAME

C     Default the file attributes so that corresponding lookup hits can
C     be ignored (and input string FSTAT can be short):

      FORM    = DICT (2)
      CC      = DICT (3)
      STATUS  = DICT (11)
      RDONLY  = .TRUE.
      DFRECL  = .TRUE.
      ENQUIRE = .FALSE.

C     Ensure that the attributes text is upper case:

      ATTRIB = FSTAT
      CALL UPCASE (ATTRIB)

C     Start of loop over tokens in attributes text:

      FIRST = 1
      LAST = LEN (FSTAT)
   50 CONTINUE

         CALL SCANNR (ATTRIB, FIRST, LAST, MARK)
         IF (MARK .EQ. 0) GO TO 70

         KEY = ATTRIB (FIRST : MARK)
         CALL LOOKUP (MXWORD, DICT, .FALSE., KEY, ENTRY)

         IF (ENTRY .GT. 0) THEN

C           There is only a modest number of possibilities.
C           Avoid replicating the dictionary text by working with
C           subscripts rather than text.  Adding keywords at the
C           end of the dictionary will not affect this code, which
C           is why LOOKUP's non-alphabetic option is used above.

            IF (ENTRY .EQ. 1 .OR. ENTRY .EQ. 10) THEN
               FORM = DICT (10)
               CC = DICT (6)
            END IF
            IF (ENTRY .EQ. 4 .OR.
     >          ENTRY .EQ. 6)    CC = DICT (ENTRY)
            IF (ENTRY .EQ. 12)   RDONLY = .FALSE.
            IF (ENTRY .EQ. 5 .OR.
     >          ENTRY .EQ. 7 .OR.
     >          ENTRY .EQ. 9)    STATUS = DICT (ENTRY)
            IF (ENTRY .EQ. 13)   ENQUIRE = .TRUE.

         ELSE IF (NUMBER (ATTRIB (FIRST : MARK))) THEN
            DFRECL = .FALSE.
            READ (ATTRIB (FIRST : MARK), '(BN, I11)') RECLEN
         ELSE
            GO TO 810
         END IF

         FIRST = MARK + 2
         IF (FIRST .LE. LAST)
     >GO TO 50

   70 CONTINUE

C     Check for inconsistencies:

      IF (STATUS .NE. DICT (7)) RDONLY = .FALSE.


  100 CONTINUE

C     Start of loop over retries at opening the file:

      IF (PROMPT .EQ. BLANK) GO TO 300

C     The prompt should indicate any default file name here.
C     Use READS instead of the original READC, for Unix reasons.

      CALL READS (LUNCRT, PROMPT, LUNKBD, FNAME, CR, EOF)

  200 CONTINUE

      IF (CR .AND. FNAME .EQ. BLANK)  GO TO 999
      IF (EOF) GO TO 800

  300 CONTINUE

C     Try to open the file, but test for its existence first in one case:

      IF (ENQUIRE) THEN
         INQUIRE (FILE=FNAME, EXIST=PRESENT)
         IF (.NOT. PRESENT) THEN
            FSTAT = 'MISSING'
            GO TO 999
         END IF
      END IF

C     Oddball READONLY keyword, and uncertain system default for RECL
C     force four possibilities here:

      IF (RDONLY .AND. DFRECL) THEN

         OPEN (UNIT=LUNFIL, FILE=FNAME, FORM=FORM, STATUS=STATUS,
     >      ERR=400)
CCCC >      CARRIAGECONTROL=CC, READONLY, ERR=400)

      ELSE IF (RDONLY .AND. .NOT.DFRECL) THEN

         OPEN (UNIT=LUNFIL, FILE=FNAME, FORM=FORM, STATUS=STATUS,
     >      RECL=RECLEN, ERR=400)
CCCC >      CARRIAGECONTROL=CC, READONLY, RECL=RECLEN, ERR=400)

      ELSE IF (.NOT.RDONLY .AND. DFRECL) THEN

         OPEN (UNIT=LUNFIL, FILE=FNAME, FORM=FORM, STATUS=STATUS,
     >      ERR=400)
CCCC >      CARRIAGECONTROL=CC, ERR=400)

      ELSE
C****    IF (.NOT.RDONLY .AND. .NOT.DFRECL) THEN

         OPEN (UNIT=LUNFIL, FILE=FNAME, FORM=FORM, STATUS=STATUS,
     >      RECL=RECLEN, ERR=400)
CCCC >      CARRIAGECONTROL=CC, RECL=RECLEN, ERR=400)

      END IF

      GO TO 999


C     OPEN error handling:

  400 CONTINUE

      WRITE (LUNCRT, 1000) 'Error opening following file:', FNAME
      FNAME = FORIG

      IF (FNAME .EQ. BLANK) THEN
         CALL READS (LUNCRT,
     >      'Try again, look in directory ("$ls ..."), cancel open ' //
     >      '(CR), or stop (^D): ', LUNKBD, FNAME, CR, EOF)
      ELSE
         CALL READS (LUNCRT,
     >      'Try again, look in directory ("$ls ..."), open default ' //
     >      '(CR), or stop (^D): ', LUNKBD, FNAME, CR, EOF)
      END IF

C     Either another attempt at the file name was entered ...

      IF (FNAME (1 : 1) .NE. '$') GO TO 200

C     ... or an operating-system command was requested (system-dependent) ...

      CALL SYSTEM (FNAME (2 :))

C     ... and start over:

      FNAME = FORIG
      GO TO 100


  800 WRITE (LUNCRT, 1000) 'Stopping as requested.'
      GO TO 990

  810 WRITE (LUNCRT, 1000) 'OPENER: Bad file attribute in FSTAT:',
     >   FSTAT (FIRST : MARK), 'Aborting.'
C**** GO TO 990

  990 STOP ' '

  999 RETURN

C     Formats:

 1000 FORMAT (1X, A)

      END
C+----------------------------------------------------------------------
C
      SUBROUTINE PARSENAM (STRING, NAMFIRST, NAMLAST, EXTFIRST, EXTLAST)
C
C     One-liner:
C
C        Isolate file name and extension if any (VMS or Unix).
C
C     Description and usage:
C
C           PARSNAME identifies the basic file name and extension (if any)
C        within a string containing a file specification.  Since it uses a
C        backward search, it may apply to VAX/VMS file names or to Unix
C        file names, and is thus portable.  STRING (NAMFIRST : NAMLAST) is
C        determined as the basic file name; STRING (EXTFIRST : EXTLAST) is
C        determined as the file extension if one is present.  For example:
C
C           dev:[dir.subdir]some_file.type;2    (VMS) and
C           /d1/dx/whatever/some.file.type      (Unix)
C                           ^       ^ ^  ^
C        would produce the same pointers as indicated.  Note that the LAST
C        "extension" is isolated under Unix, while ;n version numbers are
C        not valid for Unix.
C
C           In order for STRING (NAMFIRST : EXTLAST) always to delimit the
C        file name apart from directory and version number even if the
C        basic name or an extension is missing, it was necessary to use
C        NAMLAST = 0 and EXTFIRST = 0 as the signals for missing elements
C        - user beware.  (If BOTH are zero, at most a directory was present.)
C
C           Any device/directory or version number can thus be deduced if
C        desired, but the initial intent here is to modularize extraction
C        of the basic name for the common case of deriving a related file
C        name.  Delimination of a possible file extension, and of the full
C        name other than directory and version number, were included in
C        case they prove handy.
C
C     Arguments:
C
C        Name      Type  I/O/S  Description
C        STRING     C    I      STRING (1 : LEN (STRING)) will be scanned
C                               for a file name and extension.
C
C        NAMFIRST,  I      O    First & last characters of basic file name.
C          NAMLAST              NAMLAST = 0 signals a missing basic name as
C                               for /dir/.alias (say).  NAMFIRST = 6 here.
C
C        EXTFIRST,  I      O    First & last characters of any extension.
C          EXTLAST              EXTFIRST = 0 signals no extension, as for
C                               /dir/mydata or mydata. (say).  EXTLAST is
C                               11 or 12 in these cases respectively.
C
C     Procedures:
C
C        SCAN3      Backward-scan utility
C
C
C     Environment:  VAX/VMS or Unix, FORTRAN 77, with
C                   IMPLICIT NONE, 8-character names, and ! comments
C
C     History:
C
C        06 Jan. 1991   DAS   Initial implementation.
C
C     Author:  David Saunders, NASA Ames/Sterling Software, Mt. View, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments.

      CHARACTER
     >   STRING * (*)
      INTEGER
     >   NAMFIRST, NAMLAST, EXTFIRST, EXTLAST

C     Local variables.

      INTEGER
     >   FIRST, LAST, MARK

C     Procedures.

      EXTERNAL
     >   SCAN3

C     Execution.

C     First, dismiss any device/directory spec.
C     Allow for trailing blanks, but not undefined characters.

      LAST = LEN (STRING)

      CALL SCAN3 (STRING, ' ]/', 1, LAST, MARK)  ! Backward search

      NAMFIRST = 1
      NAMLAST = 0
      EXTFIRST = 0
      EXTLAST = LAST

      IF (MARK .EQ. 0) GO TO 99          ! Remove the most degenerate case

C     Now dismiss any version number.

      FIRST = MARK                       ! Start of significant text
      NAMFIRST = FIRST

      CALL SCAN3 (STRING, ';', FIRST, LAST, MARK)

      IF (MARK .EQ. 0) GO TO 99          ! ; is a grubby possibility
      IF (MARK .EQ. FIRST + 1) GO TO 99  ! ;n is too

      IF (MARK .GT. FIRST + 1) THEN      ! More than just a version
         LAST = MARK - 2
         EXTLAST = LAST
      END IF                             ! Else MARK = FIRST (no version #)

C     Directory and version have been cleared out of the way.
C     Now distinguish basic name from extension.  Either may be missing.

      IF (STRING (LAST : LAST) .EQ. '.') THEN  ! No extension proper
         IF (LAST .GT. FIRST) THEN
            NAMLAST = LAST - 1
         END IF                             ! And we're done either way

      ELSE                               ! Search backwards for '.'

         CALL SCAN3 (STRING, '.', FIRST, LAST, MARK)

         IF (MARK .GT. FIRST) THEN          ! Valid extension
            EXTFIRST = MARK                 ! Normal case
            IF (STRING (FIRST : FIRST) .NE. '.') THEN
               NAMLAST = MARK - 2              ! Normal case
            END IF                             ! Else name is missing
         ELSE                               ! Name only
            NAMLAST = LAST                  ! Since extension is missing
         END IF

      END IF
        
   99 RETURN
      END
C+------------------------------------------------------------------------------
C
      SUBROUTINE PLOTDEV (MODE, NAME, LUNCRT, LUNKBD, LUNDEV, EXPLAIN,
     >   NLIST, LIST, TYPELIST, PARMLIST, DEVINDEX)
C
C
C PURPOSE:
C
C        PLOTDEV prompts for and/or initializes the current graphics output
C     device, for applications using the CA-DISSPLA graphics package.  In
C     the case of metafile output, the basic file name is application-dependent
C     and hence indicated as an argument.  The appropriate file extension is
C     appended here.
C
C        This version is more elaborate than the original SMDLIB version in
C     order to permit switching devices repeatedly, as needed for flexible
C     previewing and hard-copying of multiple plot frames in a single run.
C
C        This version also provides a means of determining the graphics
C     output devices available to a given application at a given site.
C     Doing this via a configuration file in a choice of standard locations
C     is more easily extensible than hard-coding the lists of screen
C     devices and metafiles in the application.  Device attributes such
C     as X, Y scale factors needed for true inches are also handled here.
C
C
C USAGE:
C
C        Suppression of devices and selection of default devices (the first
C     in the preview/metafile lists) is determined by the configuration file
C     associated with each application.  The format of this file is indicated
C     under the DEVINDEX argument description below.
C
C        PLOTDEV's various modes are of two basic kinds: one kind sets up
C     the available list(s) of devices; the second kind picks one device
C     from a list and/or (re)initializes the device.
C
C        There are also two kinds of ways an application might do the
C     switching between previewing and hard-copy: the simple approach
C     provides for previewing one frame at a time until the user either
C     quits or requests hard-copy, in which case any input data files are
C     rewound and all frames are replotted to the metafile; the more
C     elaborate approach avoids rewinding the data by providing for
C     plotting of each frame to either one or two output devices, with
C     a "hard-copy-all-the-rest" option.  Thus an application may have
C     as few as two calls to PLOTDEV, or as many as five.
C
C        When initializing a device, the first use of PLOTDEV raises the
C     DISSPLA level from 0 to 1.  (An exception is for CGM metafile
C     initialization, which leaves the level at 0.)  In this case if
C     PLOTDEV is called at DISSPLA level 1, the level will remain at 1.
C
C
C ARGUMENTS:
C
C     ARG    TYPE/DIM   I/O/S  DESCRIPTION
C
C     MODE     C*1      I      Controls the mode of operation as follows:
C
C                |             'P': Read list of PREVIEW devices supported
C                |                  for this application at this site, from
C         Read one device           the file NAME.config in a standard place.
C       list.  NAME applies    'M': Read list of METAFILES supported from the
C       to the config. file.        configuration file.
C                |             'A': Reads list of ALL active devices from the
C                |                  configuration file (as needed by simpler
C                |                  applications).
C
C                |             'S': SELECT output device (preview or metafile)
C                |                  via a prompt; do not initialize the device.
C    Select and/or init. one   'I': INITIALIZE device according to the input
C    device. NAME applies to        device code index, without prompting.
C      the metafile if any.    'B': BOTH prompt for and initialize a device
C                |                  from the given LIST.  DEVINDEX is returned
C                |                  (see its description).
C                              
C     NAME     C*(*)    I      Name used to establish the input configuration
C                              file name and the output metafile name(s).
C                              The application program's name may serve in
C                              both cases, or the output file name may differ
C                              if NAME is used as indicated alongside the
C                              description of MODE.  Example: if NAME = 'qplot'
C                              and MODE = 'P', 'M', or 'A', the configuration
C                              file must be qplot.config in one of the standard
C                              locations shown below.
C                              For the other class of MODEs, NAME = 'qplot'
C                              would define the output PostScript file as
C                              qplot.ps.  But NAME here might instead be derived
C                              from the name of the application program's input
C                              data file, or from a user prompt.
C                               
C     LUNCRT     I      I      Unit number for screen output
C     LUNKBD     I      I      Unit number for keyboard input
C     LUNDEV     I      I      Unit number for device output (or config. file)
C
C     EXPLAIN  C*(*)    I      Explanatory text prefacing menu defined by LIST
C                              (applicable to MODE = 'S' or 'B').
C
C     NLIST      I      I/O    No. of device codes in LIST.  For MODE = 'P',
C                              'M' or 'A', NLIST should be input with the
C                              size of the array passed as LIST(*); this will
C                              be updated with the number of entries in the
C                              returned list.  For MODE = 'S' or 'I' or 'B',
C                              input the NLIST found from a prior call.
C
C     LIST   I(NLIST)   I/O    Subset of device codes defining desired menu.
C                              The first element should be the default code.
C                              Output for MODE = 'P', 'M', or 'A'; input for
C                              'S', 'I', or 'B'.  See DEVINDEX for options.
C
C   TYPELIST C(NLIST)*3   O    3-letter mnemonic device code corresponding to
C                              LIST(*) from config. file.  Output as for LIST.
C
C   PARMLIST R(NLIST,*)   O    Parameters associated with the output device.
C                              An arbitrary number per device is handled
C                              by scanning the config. file until a non-numeric
C                              token is encountered.  So far:
C                              PARMLIST (1) = X scale factor for true inches;
C                              PARMLIST (2) = Y   "     "     "     "     "
C
C   DEVINDEX     I      I/O    Device code index into the LIST arrays
C                              (input if MODE = 'I'; output if MODE = 'S' or
C                              'B'; ignored otherwise).  See definitions below.
C
C                              DEVINDEX = 0 output means no device was selected.
C
C     The configuration file format is as follows, with undesired options
C     commented out via '!'.  This should match the DATA statement below
C     (too much trouble to avoid the latter).
C
C     Standard locations searched    Samples
C
C     graph$:[NAME]                  graph$:[qplot]qplot.config
C     sys$system:                    sys$system:qplot.config
C     /usr/local/src/graphics/NAME   /usr/local/src/graphics/qplot/qplot.config
C
C     ! qplot.config file on node RALph
C     !
C     ! Type  Dev.code  Xscale  Yscale  ! Device/metafile
C     ! -------------------------------------------------
C     TEK       1       1.0     1.0     ! Tektronix 4014, default preview device
C     TEK       2       1.0     1.0     ! Tektronix 4105/4107
C     TEK       3       1.0     1.0     ! Tektronix 4109
C     ! TEK     4       1.0     1.0     ! Tektronix 4115
C     VT        5       1.0     1.0     ! VT240
C     ! LN      6       1.0     1.0     ! LN03 Plus
C     PS        7    0.996933  1.00719  ! PostScript, def. metaf.; Apple scales
C     CGM       8       1.0     1.0     ! CGM
C     ! DIS     9       1.0     1.0     ! DISSPOP
C     DIP      10       1.052   1.043   ! DIP    (QMS 1200 scale factors)
C     ! X      11       1.0     1.0     ! X Window System
C     ! <to be extended>
C
C EXTERNAL REFERENCES:
C     GETLINE    Gets one line of text, skipping over '!' lines.
C     NUMBER     Identifies (non)numeric tokens.
C     READER     Prompting utility.
C     SCANNR     Used here to trim trailing blanks.
C     UPCASE     Uppercase utility.
C     Also:      Numerous CA-DISSPLA device interfaces.
C
C ENVIRONMENT:   VAX/VMS, FORTRAN 77;  also IRIX
C                CA-DISSPLA, Version 11.0-9003
C
C HISTORY:
C
C     05/15/89   PLOTDEV (formerly SETDEV) was adapted from various ideas found
C     DAS/MDW    in SMDLIB and demo programs to solve the following problems:
C                1: IDEV had to be returned so that the screen could be cleared
C                   for SOME devices.
C                2: Several metafiles had to be handled, and these should be
C                   named after the application program.
C                3: A mechanism was needed for suppressing devices known to
C                   SMDLIB but irrelevant to the environment without losing
C                   information.  A hard-coded mapping was used as a reasonably
C                   maintainable compromise.
C
C     10/27/89   Adapted SMDLIB version to initialize devices supported by the
C     MDW        DISSPLA graphics package, including metafiles such as DIP, CGM,
C                PostScript and DISSPOP plus the option of choosing additional
C                Tektronix devices through further prompting.
C                Introduced argument MODE to separate prompting and initializ-
C                ation, and arguments LIST, NLIST to suppress devices from the
C                menu at the application level.
C
C     01/04/90   Deleted mythical GKS option.  Stripped trailing blanks from
C     MDW        device prompt.
C
C     01/23/90   Enabled hardware characters for PostScript and CGM output.
C     MDW/DAS    EXPLAIN argument found desirable to specify terminals and
C                metafiles in separate calls to PLOTDEV.
C
C     11/20/90   A correction.  Now passing 60 instead of 100 through the
C     MDW        argument list in CGM metafile calls.
C
C     01/23/91   Added configuration file scheme to keep changes to the
C     DAS        available devices from affecting applications.  Here was the
C                logical place, at some cost in explaining the different modes.
C                The argument list was reordered, and the two TYPE arguments
C                were added in the hope of keeping hard-coding of device codes
C                in application programs to a minimum.
C
C     02/19/91   Added PARMLIST argument to handle the X, Y scale factor
C     DAS        problem; changed IDEV to DEVINDEX (the relevant index into what
C                are now THREE lists); eliminated corresponding TYPE argument.
C
C     03/04/91   PostScript's pen-width was not being initialized when /PS was
C     DAS/DBS    specified on QPLOT's command line (to suppress prompting).
C                Would passing it IN via PARMLIST (3) (say) make sense?
C                For now, a data statement keeps the application simpler, but
C                forces it to invoke the menu to override the default.
C
C     03/22/91   The prompt for type of CGM file was in the wrong place, and
C     DAS        the <NAME>.cgm file was not being set up properly.
C
C     04/30/91   Had to introduce UPCASE when tried on the IRIS, which also
C     DAS        didn't like concatenation in a passed CHARACTER argument
C                (.pop and .dip file names).
C                Lots of trouble with the output PostScript file name on the
C                IRIS.  (Frames after the first went to fort.3, where LUNDEV
C                was 3.)  Had to avoid any kind of initialization if the
C                device did not change from the previous call (MODE = 'I').
C
C     05/31/91   DIP and DISSPOP metafile names were using self-counting 100s
C     DAS        where NCHARS should have been.
C
C     06/12/91   Eliminated DIP for Cray version (%REF (STRING) won't compile).
C                "Standard" location for *.config file had to be changed.
C                READONLY keyword suppressed in the OPEN for *.config.
C                IOMGR's weak way of dealing with file names is word-length
C                dependent.  LENWORD=4 or 8 and A4 or A8 are still hard-coded.
C
C     01/04/92   Metafile name(s) may not be related to the application program
C     DAS        name - made this clear in describing the two uses of NAME.
C
C AUTHORS:       David Saunders/Michael Wong
C                NASA Ames/Sterling Software, Palo Alto, CA.
C
C-------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:
C     ----------

      INTEGER
     >   DEVINDEX, LUNCRT, LUNKBD, LUNDEV, NLIST, LIST (*)

      REAL
     >   PARMLIST (NLIST, *)

      CHARACTER
     >   EXPLAIN * (*), NAME * (*), MODE * 1, TYPELIST (*) * 3

C     Local constants:
C     ----------------

      INTEGER
     >   LENSTR, LENWORD, MAXDEV, MAXTXT
      CHARACTER
     >   BLANK * 1
      PARAMETER
     >  (BLANK  = ' ',
     >   LENSTR = 60,   ! 60 >= MAXTXT + 27  (See prompt for menu item.)
                        ! NOTE: STRING is also used for the PostScript file
                        ! name set up, and thus should be at least 4*15 chars.
     >   LENWORD= 4,    ! Wordlength in bytes (8 for CRAY, 4 for DEC, SGI)
     >   MAXDEV = 10,
     >   MAXTXT = 30)   ! Length of each menu item

C     Local variables:
C     ----------------

      INTEGER
     >   FIRST, I, IDEV, IDEVLAST, IOS, ITEM, LAST, MARK, NCHARS,
     >   NWORDS, PSCALLS, IBUFF (16)
      REAL
     >   PWIDTH, XPAGE, YPAGE
      LOGICAL
     >   CR, EOF, PRE
      CHARACTER
     >   ALLDEV (MAXDEV) * (MAXTXT), CGMTYPE * 1, IFMT * 4, RFMT * 10,
     >   STRING * (LENSTR), TYP * 3

C     Procedures:
C     -----------

      LOGICAL
     >   NUMBER
      EXTERNAL
     >   NUMBER

C     Storage:
C     --------

      DATA
     +   ALLDEV /
     1   'Tektronix (monochrome) (4014) ',
     2   'Tektronix (color) (4105, 4107)',
     3   'Tek 4109                      ',
     4   'Tek 4115                      ',
     5   'VT240                         ',
     6   'LN03 Plus                     ',
     7   'PostScript metafile           ',
     8   'CGM metafile                  ',
     9   'DISSPOP metafile              ',
     +   'DIP metafile                  '/

      DATA
     >   CGMTYPE  /'B'/,
     >   IDEVLAST /-1/,     ! Initialize it to an invalid device code
     >   PSCALLS  /0/,
     >   PWIDTH   /.007/,
     >   IFMT     /'(In)'/,
     >   RFMT     /'(BN,Fnn.0)'/

      SAVE
     >   ALLDEV, CGMTYPE, IDEVLAST, IFMT, PSCALLS, PWIDTH, RFMT

C     Execution:
C     ----------

C     Either process the configuration file ...
C     -----------------------------------------

      IF (MODE .EQ. 'P' .OR. MODE .EQ. 'M' .OR. MODE .EQ. 'A') THEN

C        Look for the file NAME.config in one of several standard locations:

        DO 100, I = 1, 2
            IF (I .EQ. 1) STRING = '/sde1-fs/saunders/' // NAME // '/'
            IF (I .EQ. 2) STRING = ' '
            LAST = INDEX (STRING, BLANK)
            STRING (LAST :) = NAME // '.config'
            LAST = LAST + LEN (NAME) + 6
            OPEN (UNIT = LUNDEV, FILE = STRING (1 : LAST),
C****>            STATUS = 'OLD', IOSTAT = IOS, READONLY)
     >            STATUS = 'OLD', IOSTAT = IOS)
            IF (IOS .EQ. 0) GO TO 120          ! Success
  100    CONTINUE

         WRITE (LUNCRT, '(/, A, A, A)')
     >     ' PLOTDEV: Unable to open ', NAME, '.config file. Aborting.',
     >     BLANK
         GO TO 990

  120    CONTINUE

C        Read the config. file till EOF.  There are so few choices for
C        preview (screen) devices that they are itemized here rather
C        requiring another column (P or M) to distinguish them from metafiles.

         ITEM = 0
  140    CONTINUE

            CALL GETLINE (LUNDEV, '!', STRING, LAST, IOS)

            IF (IOS .GT. 0) THEN
               WRITE (LUNCRT, '(A)')
     >            ' PLOTDEV: Error reading config. file.  Aborting.',
     >            BLANK
               GO TO 990
            END IF

            IF (IOS  .LT. 0) GO TO 160         ! EOF
            IF (LAST .EQ. 0) GO TO 140         ! Empty or suppressed line

C           Isolate the first token:
            FIRST = 1
            CALL SCANNR (STRING, FIRST, LAST, MARK)
            TYP = STRING (FIRST : MARK)
            CALL UPCASE (TYP)
            PRE = TYP .EQ. 'TEK' .OR.
     >            TYP .EQ. 'VT'  .OR.
     >            TYP .EQ. 'X'

            IF (MODE .EQ. 'A' .OR.
     >         (MODE .EQ. 'P' .AND. PRE) .OR.
     >         (MODE .EQ. 'M' .AND. .NOT. PRE)) THEN    ! Add it to the list.
               ITEM = ITEM + 1
               TYPELIST (ITEM) = TYP

C              The second token is the integer device code:

               FIRST = MARK + 2
               CALL SCANNR (STRING, FIRST, LAST, MARK)
               WRITE (IFMT (3 : 3), '(I1)') MARK - FIRST + 1
               READ (STRING (FIRST : MARK), IFMT) LIST (ITEM)

C              Remaining numeric tokens are device parameters.
C              E.g.: X, Y scale factors for true inches.

               I = 0
  150          CONTINUE

                  FIRST = MARK + 2
                  CALL SCANNR (STRING, FIRST, LAST, MARK)
                  IF (NUMBER (STRING (FIRST : MARK))) THEN
                     I = I + 1
                     WRITE (RFMT (6 : 7), '(I2)') MARK - FIRST + 1
                     READ (STRING (FIRST : MARK), RFMT)
     >                  PARMLIST (ITEM, I)
                     GO TO 150    ! No need to check for exceeding LAST.
                  END IF

            END IF
            IF (ITEM .LT. NLIST)
     >   GO TO 140

  160    NLIST = ITEM
         CLOSE (UNIT = LUNDEV)
      END IF                  ! List is set up - return.


C     ... or select and/or initialize a device:
C     -----------------------------------------

  400 IF (MODE .EQ. 'S' .OR. MODE .EQ. 'B') THEN

C        Prompt for and assign code for device or metafile.

         WRITE (LUNCRT, 1001) BLANK, EXPLAIN, BLANK
         WRITE (LUNCRT, 1002)
     >      (ITEM, ALLDEV (LIST (ITEM)), ITEM = 1, NLIST),
     >      NLIST + 1, '(or EOF) None of the above'
         WRITE (LUNCRT, 1001)

C        Strip trailing blanks from description of default device (LIST(1)):

         FIRST = 1
         LAST  = LEN (ALLDEV (LIST (1)))
         CALL SCANNR (ALLDEV (LIST (1)), FIRST, LAST, MARK)

         STRING (1:24) = 'Item number? (Def = 1 = '
         STRING (25:24+LAST) = ALLDEV (LIST (1)) (1:LAST)
         STRING (25+LAST:27+LAST) = '): '
         ITEM = 1
         CALL READI (LUNCRT, STRING (1:27+LAST), LUNKBD, ITEM, CR, EOF)

         IF (EOF .OR. ITEM .EQ. NLIST + 1) THEN  ! Quit.
            DEVINDEX = 0
            GO TO 999
         END IF
         IF ( ITEM .LE. 0 .OR. ITEM .GT. NLIST) THEN ! Must have been a mistake.
            GO TO 400
         END IF
         DEVINDEX = ITEM

         IF (TYPELIST (ITEM) .EQ. 'PS') THEN
            PWIDTH = 7.
            CALL READR (LUNCRT,
     >         'Pen width in thousandths of an inch? (<CR> = 7): ',
     >         LUNKBD, PWIDTH, CR, EOF)
            PWIDTH = PWIDTH * .001

         ELSE IF (TYPELIST (ITEM) .EQ. 'CGM') THEN
            CALL READC (LUNCRT, 'CGM metafile type? ' //
     >         '(B(inary), C(har.), T(ext); <CR>=B): ',
     >         LUNKBD, CGMTYPE, CR, EOF)
         END IF

      END IF

      IF (MODE .EQ. 'I' .OR. MODE .EQ. 'B') THEN

         IDEV = LIST (DEVINDEX)
         IF (IDEV .EQ. IDEVLAST) GO TO 999  ! No need to do anything.
                                            ! Too much indenting to avoid GO TO.
         IDEVLAST = IDEV

C        Initialize the selected device or metafile.
C        First, TURN OFF hardware characters for most devices, and
C        set the device configuration to primary I/O (normal default):

         IF (IDEV .NE. 7 .AND. IDEV .NE. 8) THEN
            CALL RESET ('HWCHAR')
            CALL IOMGR (0, -102)
         END IF

         IF (IDEV .EQ. 1) THEN                         ! Tektronix 4014 mono

            CALL TK4014 (960, 0)

         ELSE IF (IDEV .EQ. 2) THEN                    ! Tektronix 4105/7 color

            CALL PTK41

         ELSE IF (IDEV .EQ. 3) THEN                    ! Tektronix 4109

            CALL TK41DO (1, 4109)

         ELSE IF (IDEV .EQ. 4) THEN                    ! Tektronix 4115

            CALL TK41DO (1, 4115)

         ELSE IF (IDEV .EQ. 5) THEN                    ! VT240

            CALL VT240

         ELSE IF (IDEV .EQ. 6) THEN                    ! LN03 Plus

            CALL LN01TK (LUNDEV)

         ELSE IF (IDEV .EQ. 7) THEN                    ! PostScript

            PSCALLS = PSCALLS + 1
            IBUFF (1) = 16            ! Buffer length
            CALL IOMGR (IBUFF, -1)    ! Initialize I/O system

C           Set the I/O configuration:

            CALL IOMGR (5, -102)      ! Direct data file output

C           Set up new user-supplied file name in IBUFF (1 : 15) (awkward):

            STRING = NAME
            NCHARS = MIN (LEN (NAME) + 3, 15 * LENWORD)
            STRING (NCHARS - 2 : NCHARS) = '.ps'
            NWORDS = (NCHARS + LENWORD - 1) / LENWORD

	    DO 500, I = 1, 16
               IBUFF (I) = 0
  500       CONTINUE

            DO 510, I = 1, NWORDS
               READ (STRING ((I - 1) * LENWORD + 1 : I * LENWORD),
     >            '(A4)') IBUFF (I)
  510       CONTINUE

            IBUFF (16) = NCHARS
            CALL IOMGR (IBUFF, -103)

C           Set the file mode:

            IF (PSCALLS .EQ. 1) THEN
               IBUFF (1) = 3    ! No overwrite (= new version where applicable)
            ELSE
               IBUFF (1) = 0    ! Append
            END IF

            CALL IOMGR (IBUFF, -104)

C           Set the logical unit number:

            CALL IOMGR (LUNDEV, -110)

C           Set the file carriage control and padding:

            IBUFF (1) = 0             ! Carriage control off
            IBUFF (2) = 0             ! Shouldn't be needed (cc type = space)
            IBUFF (3) = 0             ! Pad with spaces, not nulls
            CALL IOMGR (IBUFF, -111)

C           Plot area:                ! Switch between portrait and landscape
                                      ! at the application level?
            XPAGE = 7.99              ! Recommended for Apple LaserWriter.
            YPAGE = 10.78

            CALL PSCRPT (XPAGE, YPAGE, PWIDTH)

            IF (PSCALLS .EQ. 1) THEN  ! Reset to append, else we get a new
               IBUFF (1) = 0          ! PostScript file for each plot on VAX.
               CALL IOMGR (IBUFF, -104)
            END IF

C           Specify hardware characters for all available alphabets,
C           with reasonable tolerances:

            CALL HWCHAR (20., 20., 90., 1, 'HEBREW')

         ELSE IF (IDEV .EQ. 8) THEN                    ! CGM metafile

            STRING = NAME
            NCHARS = LEN (NAME) + 4
            STRING (NCHARS - 3 : NCHARS) = '.cgm'

            IF (CGMTYPE .EQ. 'B') CALL CGMBO (STRING, NCHARS, 0)
            IF (CGMTYPE .EQ. 'C') CALL CGMCO (STRING, NCHARS, 0)
            IF (CGMTYPE .EQ. 'T') CALL CGMTO (STRING, NCHARS, 0)

            CALL HWCHAR (20., 20., 90., 1, 'HEBREW')

         ELSE IF (IDEV .EQ. 9) THEN                    ! DISSPOP metafile

            STRING = NAME
            NCHARS = LEN (NAME) + 4
            STRING (NCHARS - 3 : NCHARS) = '.pop'
            CALL COMPRS
            CALL SETCPR (LUNDEV, 0, 0, 0)
            CALL POPNAM (STRING, NCHARS)

         ELSE IF (IDEV .EQ. 10) THEN                   ! DIP metafile

            STRING = NAME
            NCHARS = LEN (NAME) + 4
            STRING (NCHARS - 3 : NCHARS) = '.dip'
C*****      CALL DIP (LUNDEV, %REF (STRING), NCHARS)
            WRITE (LUNCRT, 1001) 'PLOTDEV: DIP is no longer supported.'
            GO TO 990

         END IF
      END IF

      GO TO 999

  990 STOP ' '

  999 RETURN

C     Formats:
C     --------

 1001 FORMAT (1X, A)
 1002 FORMAT (3X, I2, ': ', A)

      END
C+----------------------------------------------------------------------
C
      SUBROUTINE PLSFIT (NDATA, X, Y, TBEGIN, TEND, NEVAL, XEVAL, YEVAL,
     &   NEW, CLOSED, METHOD, DISTRIB, IER)
C
C     One-liner:  Storage-efficient parametric piecewise cubic fit (2-space)
C     ----------
C
C     Description and usage:
C     ----------------------
C
C        PLSFIT (Parametric Local Spline FIT) is intended for general
C     curve drawing applications where it is not necessary to provide
C     spline coefficients as output. It is also suited to interpolating
C     the same data repeatedly at maximum efficiency (as in interactive
C     graphics, perhaps). Application to grid generation is another
C     possibility (as in redistributing X and Y vs. T where preferred
C     values of T are determined externally).
C
C        "Local" piecewise methods typically produce continuity in the
C     function and its first derivative across the data points, but the
C     second derivatives are not guaranteed to be continuous.  With
C     reasonably smooth data, this limitation should not be a problem.
C
C        Since it avoids storing spline coefficients, PLSFIT requires
C     very little working storage, yet its speed is on a par with other
C     methods in a single pass. The routine was designed for evaluation
C     of either many points at once (preferred) or one point at a time.
C     Computational efficiency is achieved by retaining some internal
C     quantities between calls for re-use where appropriate.
C
C        The primary method employed, references (1)-(3), uses monotone
C     piecewise cubic interpolation for X and Y vs. T, so the resulting
C     curve follows the data closely. In fact, the fit is sometimes too
C     tight, especially where X or Y appears constant over a region
C     which is intended to be curved. (A circle with unfortunate choice
C     of data points can develop flat spots.)  A partial solution is
C     offered by the "Bessel" option, which produces a more rounded
C     interpretation of the data. Both techniques are based on matching
C     the first (but not second) derivatives of the cubics on adjacent
C     intervals. The structure of PLSFIT permits easy addition of other
C     local, 3-point methods for these derivatives.
C
C        Output arrays XEVAL and YEVAL are filled by sampling the fit at
C     intervals along the arc's length, which is estimated by summing
C     the chords between adjacent points. For curve drawing (NEVAL large),
C     this version of PLSFIT offers only equal increments of arclength,
C     but other options may be added in the future (a method based on
C     curvature, for example). Other applications may be require calling
C     the subroutine repeatedly with NEVAL = 1 to evaluate the fit at
C     externally-generated values of T.  The evaluation can proceed
C     efficiently "backwards" (decreasing T) if desired.
C
C
C     Curve-drawing example:
C     ----------------------
C
C        PLSFIT usage, with extra-careful prescaling to avoid possible
C     overflow or underflow in chord length calculations, might look like
C     this for a curve-drawing application:
C
C        CALL RANGER (..., X, XSCALE, XSHIFT)   ! Scan input arrays
C        CALL RANGER (..., Y, YSCALE, YSHIFT)
C
C        CALL RESCALE (..., X, XSCALE, XSHIFT)  ! Transform data to [0, 1]
C        CALL RESCALE (..., Y, YSCALE, YSHIFT)
C
C        CALL PROTECT (..., DISTINCT)           ! Check for zero chords
C        IF (.NOT.DISTINCT) GO TO 900
C
C        NEW     = .TRUE.                       ! Initialize PLSFIT parameters
C        CLOSED  = .FALSE.
C        DISTRIB = 'U'                          ! Uniform (upper case)
C        METHOD  = 'M'                          ! Monotone ( "     " )
C        NFIT    = 500                          ! Some reasonable number
C        TBEGIN  = -1.0                         ! Interpolate entire curve
C        TEND    = -1.0
C
C        CALL PLSFIT (NXY, X, Y, TBEGIN, TEND, NFIT, XFIT, YFIT, NEW,
C       &   CLOSED, METHOD, DISTRIB, IER)
C        IF (IER .NE. 0) GO TO 910
C
C        Restore result using the inverse transformation:
C
C        CALL RESCALE (..., XFIT, ONE / XSCALE, -XSHIFT / XSCALE)
C        CALL RESCALE (..., YFIT, ONE / YSCALE, -YSHIFT / YSCALE)
C        (May need to restore the data also.)
C
C     Utilities RANGER, RESCALE, and PROTECT are available from the author.
C     (More recent utilities GETSCALE and USESCALE may be preferable now.)
C     They should only be needed if the magnitude and spacing of the data
C     points are unsuitable for ordinary single-precision calculations (all
C     too frequently the case, alas!). CHORD, which is used repeatedly by
C     PLSFIT, may also be invoked by the calling routine to provide initial
C     and final arclengths associated with particular array indices:
C
C        TBEGIN = CHORD (X, Y, 1, IBEGIN)
C        TEND   = TBEGIN + CHORD (X, Y, IBEGIN, IEND)
C
C     If the entire curve is to be interpolated at internally-generated
C     Ts, it is more efficient to pass -1.0s (indicating default choice
C     of arclength) and let PLSFIT compute (and re-use) the total arclength.
C
C     
C     Grid generation example:
C     ------------------------
C
C        A grid generation application may need to compute total arclength,
C     redistribute the intermediate values of T, and interpolate to these
C     by calling PLSFIT once per value of T.  Again, using PLSFIT to do
C     the initial calculation of TEND is preferable to direct use of CHORD,
C     because PLSFIT needs to do other initialization on its first call
C     if unnecessary recalculation of TEND is to be avoided:
C
C        Determine total arclength:
C
C        NEW     = .TRUE.                  ! Initialize PLSFIT parameters;
C        NEVAL   = 1                       ! force one pass through the
C                                          ! search to initialize it too.
C        TBEGIN  = -1.0                    ! Causes TEND to be returned;
C        TEND    = -1.0                    ! TBEGIN returns as zero.
C        ..............
C
C        CALL PLSFIT (NXY, X, Y, TBEGIN, TEND, NEVAL, XEVAL, YEVAL, NEW,
C       &   CLOSED, METHOD, DISTRIB, IER)
C        IF (IER .NE. 0) GO TO 910
C
C        Redistribute the total arclength:
C
C        CALL DISTRIB (..., NEVAL, ZERO, TEND, ..., T, ...)
C
C        Evaluate the redistributed X and Y values one pair at a time
C        to avoid PLSFIT's internally-generated Ts:
C
C        NEW = .FALSE.
C        DO 300, I = 1, NEVAL
C           CALL PLSFIT (NXY, X, Y, T (I), T (I), 1, XEVAL (I), ...)
C    300 CONTINUE
C        ..............
C
C
C     Arguments:
C     ----------
C
C     Name    Type/Dimension  I/O/S  Description
C     NDATA   I               I      Length of X, Y input data arrays.
C
C     X,      R (NDATA)       I      Input data coordinates. Successive
C     Y                              data points must be distinct in the
C                                    sense that DX ** 2 + DY ** 2 > 0, but
C                                    no check is performed at this level.
C                                    If CLOSED, then first and last points
C                                    must agree (not checked here).
C
C     TBEGIN  R               I/O    Specifies the arclength at which to
C                                    start interpolation. If NEVAL = 1,
C                                    TBEGIN is the only point at which
C                                    the fit will be evaluated. If both
C                                    TBEGIN and TEND are negative, a
C                                    default value of zero will be used
C                                    for TBEGIN.  In this case, TBEGIN
C                                    is returned as zero and TEND is
C                                    returned as the total arclength.
C
C     TEND    R               I/O    The arclength at which to end the
C                                    interpolation, unless both TBEGIN
C                                    and TEND are negative, in which case
C                                    the total chord length will be used
C                                    for TEND. If NEVAL = 2, TEND is the
C                                    second (and last) at which the fit
C                                    will be evaluated.  See TBEGIN and
C                                    NEVAL.
C
C     NEVAL   I               I      Number of points at which the fit is
C                                    to be evaluated.  NEVAL >= 1.
C                                    NEVAL > 1 means PLSFIT generates its
C                                    own values of T according to DISTRIB;
C                                    NEVAL = 1 with TBEGIN, TEND < 0. enables
C                                    PLSFIT to calculate total arclength
C                                    while initializing itself for further
C                                    interpolation at externally-supplied
C                                    values of T.  See the grid-generation
C                                    example above.
C
C     XEVAL,  R (NEVAL)         O    Interpolated coordinates.
C     YEVAL
C
C     NEW     L               I      If control flag NEW is .TRUE., the
C                                    search for a bracket starts from
C                                    scratch, otherwise locally-saved
C                                    search and fit information will be
C                                    assumed to be correct. If calling
C                                    PLSFIT from within a loop, set
C                                    NEW = .FALSE. after the first call.
C
C     CLOSED   L              I      Logical flag indicating that periodic
C                                    boundary conditions are to be used.
C                                    (The curve should wrap around smoothly
C                                    on itself.). The calling routine must
C                                    ensure that the ends agree.
C
C     METHOD   C*1            I      The type of fit to be used:
C                                    'M' means monotonic piecewise cubics;
C                                    'B' means non-monotonic "Bessel"-type
C                                    piecewise cubics (looser fit).
C                                    METHOD must be uppercase.
C
C     DISTRIB  C*1            I      Controls the distribution of the
C                                    interpolated points along the curve.
C                                    Irrelevant if NEVAL = 1, which is
C                                    appropriate for externally-generated
C                                    values of arclength for interpolation.
C                                    The only choice at present is 'U' for
C                                    uniform arclength increments. DISTRIB
C                                    must be uppercase.
C
C     IER      I                O    Error flag.  The possibilities are:
C                                       0:  No problems
C                                      +1:  NDATA < 2 (all cases), or
C                                           NDATA < 3 (periodic case)
C                                      +2:  Ends don't match (periodic)
C                                      +3:  Bad METHOD requested
C                                      +4:  Bad DISTRIB requested
C
C     Significant local variables:
C     ----------------------------
C
C     MEMORY         Indicates that search information and cubic coefficients
C                    are correct for the current point.
C
C     NEWFLAG        A local copy of NEW which is set .FALSE. during the
C                    first entry to the interval search routine, ARCSRCH,
C                    which SAVEs several local variables from one call to
C                    the next.
C
C     LEFT, RIGHT    Current endpoints of the bracketing interval.
C     TLEFT, TRIGHT  Corresponding cumulative arclengths.
C
C     IND            Index array which keeps track of the endpoints of
C                    the intervals preceding and following the bracketing
C                    interval, with wrap-around at the extremes.
C
C     H, DEL         Arclength and forward difference derivative arrays.
C
C     BX, CX, DX     Coefficients of X-cubic on the bracketing interval.
C     BY, CY, DY     Ditto, for Y-cubic.
C
C     Procedures:
C     -----------
C
C     ARCSRCH   Interpolation search along chord of a parametric curve.
C     BESSEL    First derivative using central 3-point formula.
C     BRODLIE   First derivative (central), adjusted for monotonicity.
C     BUTLAND   First derivative (non-central), adjusted for monotonicity.
C     CHORD     Summed chord-lengths for X-Y curve over range of indices.
C     THREEPT   First derivative using non-central 3-point formula.
C
C     Environment:  Digital VAX-11/780, VMS FORTRAN 77
C     ------------  Apple Macintosh, Absoft MacFORTRAN/020
C
C     Notes:
C     ------
C
C     (1)  IMPLICIT NONE, 8-character symbolic names, and "!" as comment
C          character are not (yet) standard.
C
C     (2)  Since many of the calculations must be repeated at both ends
C          of an interval, and for both X and Y data, the various finite
C          difference quantities used are stored as arrays. The following
C          "map" of a typical interior interval and its neighbors should
C          help in understanding the notation.  The local array indices
C          are all numbered relative to the left-hand end of the interval
C          which brackets the point to be evaluated.
C
C                                  LEFT       RIGHT
C
C          Point         -1          0         +1          +2
C
C          Data           x -------- x -------- x --------- x
C
C          Interval      -1          0         +1
C
C     (3)  Within PLSFIT, we have allowed for extrapolation to simplify
C          things, but it's not intended for general use (especially for
C          the periodic case where we don't attempt to wrap around past
C          the last interval).
C
C     (4)  Error-checking philosophy: several cheap tests are performed
C          upon entry (see IER description, above), but the user must
C          guarantee that successive points are distinct to avoid division
C          by arclength intervals of length zero. Such a test is best
C          performed at a higher level to avoid the considerable overhead
C          of scanning both data arrays on each entry.
C
C     (5)  Note that if some function of the data is to be plotted,
C          e.g., logarithm, then a smoother curve will be obtained by
C          passing transformed data, fitting at equal increments, then 
C          reverse-transforming the result.
C
C     Bibliography:
C     -------------
C
C     (1)  Brodlie, K. W.  A Review of Methods for Curve and Function
C             Drawing, in Mathematical Methods in Computer Graphics and
C             Design, ed. K. W. Brodlie.  London: Academic Press, 1980.
C             Pp. 1-37.  (See esp. the discussion, pp. 33-37)
C
C     (2)  Fritsch, F. N., and J. Butland.  A Method for Constructing
C             Local Monotone Piecewise Cubic Interpolants.  SIAM J. Sci.
C             Stat. Comput., Vol. 5, No. 2 (June 1984).  Pp. 300- 304.
C
C     (3)  Fritsch, F. N., and R. E. Carlson.  Monotone Piecewise Cubic
C             Interpolation.  SIAM J. Num. Anal., Vol. 17, No. 2
C             (April 1980).  Pp. 238-246.
C
C     (4)  Sedgewick, R.  Algorithms.  Reading: Addison-Wesley, 1983.
C             Chap. 14.  (Interpolation search)
C
C     Author:  Robert Kennelly, ex Sterling Software, now RAC Branch
C     -------
C              Mail Stop 227-2
C              NASA-Ames Research Center
C              Moffett Field, CA  94035
C
C              Phone (415) 604-5860
C
C     History:
C     --------
C
C     27 Feb. 1987    RAK    Initial design and coding.
C     16 Apr. 1987    RAK    Repaired handling of the degenerate linear
C                            case - DELX and DELY were not being set.
C     22 Apr. 1987    RAK    Fixed typo in linear case (DY not set).
C     25 Mar. 1988    RAK    Modularized ARCSRCH, with repairs (now works
C                            bidirectionally). Both TBEGIN and TEND have
C                            to be negative to signal use of defaults
C                            (zero and total chord length). Pass total
C                            arclength to ARCSRCH on first entry to
C                            avoid redundant calculation.
C     12 Apr. 1988    RAK    Set MEMORY = .TRUE. if possible when PLSFIT
C                            is called repeatedly with the same data.
C     23 Aug. 1989    DAS    Application to grid generation revealed that
C                            NEVAL = 0 (for arclength calculation) didn't
C                            initialize the search properly - had to go
C                            to NEVAL = 1.  Descriptions above revised
C                            to cover both graphics and grid generation.
C     20 June 1991    DAS    THREEPT was renamed BUTLAND (monotonic);
C                            THREEPT (pure 3-pt. formula) is now used
C                            for the loose fit at the boundaries.
C                            
C-----------------------------------------------------------------------

C     Declarations.
C     -------------

      IMPLICIT NONE

C     Arguments.

      LOGICAL
     &   CLOSED, NEW
      INTEGER
     &   IER, NDATA, NEVAL
      REAL
     &   TBEGIN, TEND, X (NDATA), XEVAL (NEVAL), Y (NDATA),
     &   YEVAL (NEVAL)
      CHARACTER
     &   DISTRIB * 1, METHOD * 1

C     Local constants.

      REAL
     &   ZERO, ONE, TWO, THREE
      PARAMETER
     &  (ZERO  = 0.0E+0,
     &   ONE   = 1.0E+0,
     &   TWO   = 2.0E+0,
     &   THREE = 3.0E+0)

C     Local variables.

      LOGICAL
     &   MEMORY, MONO, NEWFLAG
      INTEGER
     &   IEVAL, IND (-1:2), J, LEFT, RIGHT
      REAL
     &   BX (0:1), BY (0:1), CX, CY, DT, DELX (-1:1), DELY (-1:1),
     &   DX, DY, H (-1:1), RH, TEVAL, TINC, TLEFT, TRIGHT, TTOTAL

C     Procedures.

      REAL
     &   BESSEL, BRODLIE, BUTLAND, CHORD, THREEPT
      EXTERNAL
     &   BESSEL, BRODLIE, BUTLAND, CHORD, THREEPT

C     Storage.

      SAVE
     &   BX, BY, CX, CY, DX, DY, LEFT, RIGHT, TLEFT, TRIGHT, TTOTAL

C     Error checking and initialization.
C     ----------------------------------

      IER = 0
      IF (NDATA .LT. 2) THEN
         IER = +1
         GO TO 990
      END IF

      IF (CLOSED) THEN

C        At least three points are required, and the first and last
C        points must be identical.

         IF (NDATA .LT. 3) IER = +1
         IF (X (1) .NE. X (NDATA) .OR. Y (1) .NE. Y (NDATA)) IER = +2
      END IF

C     We'll need the total arclength to initialize the search efficiently.

      IF (NEW) THEN
         TTOTAL = CHORD (X, Y, 1, NDATA)
         TRIGHT = TTOTAL
      END IF

C     The whole curve will be interpolated by default if the calling
C     routine passes negative quantities for both TBEGIN and TEND. 

      IF (TBEGIN .LT. ZERO .AND. TEND .LT. ZERO) THEN
         TBEGIN = ZERO
         TEND   = TTOTAL
      END IF

C     Check the requested interpolation method and point distribution.

      MONO = METHOD .EQ. 'M'
      IF (.NOT. MONO .AND. METHOD .NE. 'B') IER = +3

      IF (NEVAL .GT. 1) THEN           ! NEVAL = 1 for T supplied externally.
         IF (DISTRIB .EQ. 'U') THEN    ! Uniform values of T generated here.
            TINC = (TEND - TBEGIN) / REAL (NEVAL - 1)
         ELSE
            IER = +4
         END IF
      END IF

C     Bail out if any error was found. (Only the last found is reported.)
C     Original test on NEVAL = 0 here was invalid - have to go through ARCSRCH.

      IF (IER .NE. 0) GO TO 990


C     Initialize bracket quantities. Note that when NEW = .TRUE. the
C     MEMORY flag is never set, thus ARCSRCH (below) gets initialized
C     properly, as long as it is entered at least once (NEVAL >= 1).

      TEVAL   = TBEGIN
      NEWFLAG = NEW

      IF (NEWFLAG) THEN
         MEMORY = .FALSE.
      ELSE
C        We can save a lot of time when PLSFIT is being called from within
C        a loop by setting MEMORY if possible. The out-of-range checking
C        relies on the fact that RIGHT = LEFT + 1 after the last return.
C        Cater to the more likely case of TEVAL in the previous, interior
C        interval.

         MEMORY = (TEVAL .GE. TLEFT) .AND. (TEVAL .LT. TRIGHT)
         IF (.NOT. MEMORY) THEN
            MEMORY = (LEFT  .EQ. 1)     .AND. (TEVAL .LT. TRIGHT) .OR.
     &               (RIGHT .EQ. NDATA) .AND. (TEVAL .GE. TLEFT)
         END IF
      END IF

C     Loop over evaluation points.
C     ----------------------------

      IEVAL = 0
   10 CONTINUE
         IEVAL = IEVAL + 1

         IF (MEMORY) GO TO 70        ! Skip the bulk of the computation.


C        Interpolation search for bracketing interval.
C        ---------------------------------------------

         CALL ARCSRCH (NDATA, X, Y, TEVAL, LEFT, TLEFT, RIGHT,
     &      TRIGHT, NEWFLAG)

C         -------------------------------------------------------------
C        |                                                             |
C        |   1 <= LEFT < RIGHT = LEFT + 1 <= NDATA                     |
C        |                                                             |
C         -------------------------------------------------------------

C        Compute derivatives by finite-differences.
C        ------------------------------------------

         IF (NDATA .GT. 2) THEN

C           Three cases are handled together: (1) both endpoints are
C           interior to the data range, or the interval is at the
C           beginning/end of the range with either (2) periodic,
C           or (3) free end conditions.  For special cases (2) and (3),
C           the initial calculations will be overridden below.

            IND (-1) = LEFT - 1
            IND ( 0) = LEFT
            IND (+1) = RIGHT
            IND (+2) = RIGHT + 1

C           Patch index array for periodic case (2). (Later, the free end
C           case will overridden again, following derivative calculation.)

            IF (LEFT .EQ. 1) THEN

C              Left side wrap-around boundary condition.

               IND (-1) = NDATA - 1
            ELSE IF (RIGHT .EQ. NDATA) THEN

C              Right side.

               IND (2) = 2
            END IF

C           Interval and derivative approximations.
C           ---------------------------------------

C           Eliminate possible division by zero due to cancellation while
C           subtracting by computing the chord from LEFT to RIGHT explicitly.

            H (-1) = CHORD (X, Y, IND (-1), IND (-1) + 1)
            H ( 0) = CHORD (X, Y, LEFT, RIGHT)
            H (+1) = CHORD (X, Y, IND (2) - 1, IND (2))

            DO 40, J = -1, +1
               RH       = ONE / H (J)
               DELX (J) = (X (IND (J + 1)) - X (IND (J))) * RH
               DELY (J) = (Y (IND (J + 1)) - Y (IND (J))) * RH
   40       CONTINUE

C           Select interpolation scheme.
C           ----------------------------

C           Compute adjusted X and Y derivatives at both left- and
C           right-hand endpoints of the interval.

            IF (MONO) THEN

C              Monotone - use Brodlie modification of Butland's
C              formula to adjust the derivatives at the knots.

               DO 50, J = 0, +1
                  BX (J) = BRODLIE (J, H, DELX)
                  BY (J) = BRODLIE (J, H, DELY)
  50           CONTINUE

            ELSE                             ! IF (METHOD .EQ. 'B') THEN

C              Bessel - use central difference formula at the knots.

               DO 60, J = 0, +1
                  BX (J) = BESSEL (J, H, DELX)
                  BY (J) = BESSEL (J, H, DELY)
  60           CONTINUE
            END IF

C           Patch initial/final derivatives if not periodic, case (3).

            IF (.NOT. CLOSED) THEN
               IF (LEFT .EQ. 1) THEN
                  IF (.NOT. MONO) THEN
                     BX (0) = THREEPT (0, H, DELX)
                     BY (0) = THREEPT (0, H, DELY)
                  ELSE
                     BX (0) = BUTLAND (0, H, DELX)
                     BY (0) = BUTLAND (0, H, DELY)
                  END IF
               ELSE IF (RIGHT .EQ. NDATA) THEN
                  IF (.NOT. MONO) THEN
                     BX (1) = THREEPT (1, H, DELX)
                     BY (1) = THREEPT (1, H, DELY)
                  ELSE
                     BX (1) = BUTLAND (1, H, DELX)
                     BY (1) = BUTLAND (1, H, DELY)
                  END IF
               END IF
            END IF

C           Compute the remaining cubic coefficients for X and Y relative
C           to the left-hand endpoint.

            RH = ONE / H (0)
            CX = (THREE * DELX (0) - TWO * BX (0) - BX (1)) * RH
            CY = (THREE * DELY (0) - TWO * BY (0) - BY (1)) * RH
            DX = (-TWO * DELX (0) + BX (0) + BX (1)) * RH ** 2
            DY = (-TWO * DELY (0) + BY (0) + BY (1)) * RH ** 2

         ELSE                                   ! IF (NDATA .EQ. 2) THEN

C           Degenerate case (linear).
C           -------------------------

            H (0)  = TRIGHT - TLEFT
            RH     = ONE / H (0)
            BX (0) = (X (RIGHT) - X (LEFT)) * RH
            BY (0) = (Y (RIGHT) - Y (LEFT)) * RH
            CX     = ZERO
            CY     = ZERO
            DX     = ZERO
            DY     = ZERO
         END IF

C        Evaluate the cubics for XEVAL and YEVAL.
C        ----------------------------------------

   70    CONTINUE

         DT            = TEVAL - TLEFT
         XEVAL (IEVAL) = X (LEFT) + DT * (BX (0) + DT * (CX + DT * DX))
         YEVAL (IEVAL) = Y (LEFT) + DT * (BY (0) + DT * (CY + DT * DY))

C        Choose next evaluation point and loop back.
C        -------------------------------------------

         IF (IEVAL .LT. NEVAL) THEN           ! Skips this if NEVAL = 1.

C           Uniform spacing - this fast, simple update should suffice
C           for graphics despite roundoff accumulation, especially since
C           extrapolation is allowed. (NOTE: Add test for DISTRIB here if
C           other distribution options are added.)

            TEVAL  = TEVAL + TINC
            MEMORY = (TEVAL .LT. TRIGHT) .AND. (TEVAL .GE. TLEFT)
            GO TO 10

         END IF

C     Termination.
C     ------------

  990 CONTINUE
      RETURN
      END
C+----------------------------------------------------------------------
C
      SUBROUTINE ARCSRCH (N, X, Y, T, LEFT, TLEFT, RIGHT, TRIGHT, NEW)
C
C     One-liner: Interpolation search along chord of a parametric curve.
C     ----------
C
C     Description and usage:
C     ----------------------
C
C        A modular implementation of the arclength interval search
C     originally written for PLSFIT (Parametric Local Spline FIT). An
C     arclength interval is sought which contains some specified value,
C     T, or which is at least the nearest interval if T is off-scale.
C     The condition for a bracket is: TLEFT <= T <= T (LEFT+1). A logical
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
C     N                    I    I      Number of points in arrays X, Y;
C                                      must be >= 2; no check performed.
C
C     X, Y     N           R    I      Array of distinct points defining
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
C                                      This trick permits PLSFIT & ARCSRCH
C                                      to avoid unnecessary recalculation.
C
C     NEW                  L    I/O    Must be set .TRUE. when ARCSRCH is
C                                      first called for any dataset.
C
C                                      Set to .FALSE. on output.
C
C     Environment:  Digital VAX-11/780, VMS FORTRAN
C     ------------  Apple Macintosh, Absoft MacFORTRAN/020 v2.3
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
C     25 Mar. 1988    RAK    Adapted from INTERVAL and PLSFIT.
C     24 Aug. 1988    RAK    Added protection against roundoff-induced
C                            failure in main loop (LEFT = RIGHT). Revised
C                            loop for single termination (pure WHILE).
C                            Take advantage if TEMP2 is already known.
C                            Retain URTRIGHT instead of re-calculating.
C     10 July 1989    DAS    URTRIGHT was integer by mistake.
C
C-----------------------------------------------------------------------

C     Declarations.
C     -------------

      IMPLICIT NONE

C     Constants.

      REAL
     &   ZERO
      PARAMETER
     &  (ZERO = 0.0E+0)

C     Arguments.

      LOGICAL
     &   NEW
      INTEGER
     &   LEFT, N, RIGHT
      REAL
     &   TLEFT, TRIGHT, X (N), Y (N), T

C     Local variables.

      INTEGER
     &   LENGTH, NLESS1, TRIAL
      REAL
     &   LSENTRY, RSENTRY, TEMP1, TEMP2, URTRIGHT

C     Storage.

      SAVE
     &   NLESS1, LSENTRY, RSENTRY, URTRIGHT

C     Procedures.

      REAL
     &   CHORD
      EXTERNAL
     &   CHORD

C     Execution.
C     ----------

      IF (NEW) THEN

C        Initialization for new data set. Note special use of TRIGHT here,
C        to avoid having to repeat the relatively expensive summation over
C        the entire curve.

         NEW      = .FALSE.
         URTRIGHT = TRIGHT
         NLESS1   = N - 1

         LSENTRY  = CHORD (X, Y, 1, 2)
         RSENTRY  = TRIGHT - CHORD (X, Y, NLESS1, N)

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
     &         (T - TLEFT) / (TRIGHT - TLEFT))))

C            ----------------------------------------------------------
C           |                                                          |
C           |   2 <= LEFT <= TRIAL < RIGHT <= N - 1                    |
C           |                                                          |
C            ----------------------------------------------------------

C           Compute the cumulative arclength up to TRIAL as cheaply as
C           possible using previous results. (The search runs equally well
C           for an increasing or decreasing sequence of T's.)

            IF ((TRIAL - LEFT) .LE. (RIGHT - TRIAL)) THEN
               TEMP1 = TLEFT + CHORD (X, Y, LEFT, TRIAL)
            ELSE
               TEMP1 = TRIGHT - CHORD (X, Y, TRIAL, RIGHT)
            END IF

C           Similar trick for arclength up to TRIAL + 1 since we may well
C           already know it, e.g., just before termination.

            IF (RIGHT .EQ. TRIAL + 1) THEN
               TEMP2 = TRIGHT
            ELSE
               TEMP2 = TEMP1 + CHORD (X, Y, TRIAL, TRIAL + 1)
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
                  TLEFT  = TRIGHT - CHORD (X, Y, LEFT, RIGHT)
               END IF
            ELSE IF (T .LT. TEMP1) THEN

C              Decrease RIGHT.

               IF (TRIAL .GT. LEFT) THEN
                  RIGHT  = TRIAL
                  TRIGHT = TEMP1
               ELSE
                  RIGHT  = LEFT + 1
                  TRIGHT = TLEFT + CHORD (X, Y, LEFT, RIGHT)
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
C+------------------------------------------------------------------------------
C
      SUBROUTINE POLYLINE (N, X, Y, LINE, THICK, SYMBOL, COLOR, IER)
C
C One-liner:  
C ----------
C     Connect-the-dots with given pattern, thickness, symbol, color (level 3)
C
C Description and usage:
C ----------------------
C
C        POLYLINE is a high-level "curve" drawing routine with a generic
C     programmer interface.  It isolates at least some of the details of
C     the specific graphics library in use.  (Plot set-up functions still
C     have to be translated elsewhere when switching to another library.)
C
C        This version is for the DISSPLA package (originally developed
C     for the PLOTCL application and since adapted for QPLOT and FMAP).
C
C        Note that any effect POLYLINE might have on the plot legend (via
C     its call to CURVE) is up to the application.  The original DISSPLA
C     approach would initiate this by calling LEGLIN, LINEST, and LINES prior
C     to any call to POLYLINE.  The related LEGEND utility (an alternative
C     to DISSPLA's original LEGEND) eliminates any "hidden" connections between
C     POLYLINE/CURVE and the legend, so that awkward work-arounds for (say)
C     certain curve fitting methods are no longer needed.
C
C Arguments:
C ----------
C
C     Name    Dimension  Type  I/O/S  Description
C     N                   I    I      Number of points to be plotted.
C
C     X          N        R    I      Array of abscissas.
C
C     Y          N        R    I      Array of ordinates.
C
C     LINE                I    I      Line type code:
C                                        1 = connected symbols
C                                        2 = symbols only
C                                        3 = solid line
C                                        4 = dots
C                                        5 = dashes
C                                        6 = chain dots
C                                        7 = chain dashes
C                                        8 = long dashes
C                                        9 = thick solid line (.02 inch)
C
C     Notes:  1)  LINE < 0 indicates a closed curve.  E.g.  Line = -4 is a
C                 closed dotted line.
C             2)  Thickness factor is ignored if the line type code = 2
C                 (symbols only), or if the line type code = 9 (thick line 
C                 = .02 inch).  (LINE = 9 was retained for upward compatibility;
C                 LINE = 3 combined with THICK = 2 gives the same result.)
C
C     THICK               I    I      Line thickness factor (where pen width 
C                                     = THICK x .01 inch).  Note:  If THICK = 0,
C                                     DISSPLA default thickness is assumed.
C
C     SYMBOL              I    I      Symbol type.  Legal values are
C                                     [0,18], with automatic wraparound
C                                     if > 18.  See DISSPLA manual for
C                                     key.  Use negative for no symbol.
C
C     COLOR               R    I      Indicates the color to be
C                                     used for the curve:
C                                        0. = white
C                                        1. = black
C                                        2. = magenta
C                                        3. = red
C                                        4. = yellow
C                                        5. = green
C                                        6. = cyan
C                                        7. = blue
C
C                                      with colors between the six primaries 
C                                      specified using the real numbers in
C                                      [2., 7.].  (E.g. 3.5 = orange, 4.5 =
C                                      yellow-green, etc.)
C
C     IER                 I      O    Error flag:
C                                        0 = no problems
C                                        1 = N is not positive
C                                        2 = no line or symbol asked for
C
C Environment:  Digital VAX-11/780 VMS Version 4.4 (FORTRAN 77)
C ------------  Computer Associates DISSPLA Version 10.0
C
C Notes:
C ------
C     (1)  IMPLICIT NONE and 8-character symbolic names are non-standard.
C
C     (2)  See the DISSPLA manual for the order of the plotting symbols
C          used.  The first few are: Square, Octagon, Triangle ("up"),
C          Plus, X, Diamond, and Triangle ("down").
C
C Author:  Robert Kennelly, Sterling Software/NASA-Ames
C -------
C
C History:
C --------
C     08/27/87   R.A.Kennelly   Design and coding, based on PLCURV, fragments
C                               of QUICK, and incorporating pieces of MYSPEC.
C     06/23/89   M.D.Wong,      Header clarified in light of SMDLIB version;
C                D.A.Saunders   chaindot, chaindash redefined as for longdash.
C                               Symbols may now be plotted in combination with
C                               non-solid line patterns.
C     08/29/89   M.D.Wong       Updated to run on DISSPLA version 11.0. (%REF
C                               taken out of call to SETCLR.)
C     12/06/89   M.D.Wong       Ideas from FMCURV included here for contour
C                               plotting.  Line thicknesses may now be varied
C                               via the addition of THICK to the argument
C                               list.  Curve closing option added.  COLOR made 
C                               type real so that line colors may be specified 
C                               between the primaries.
C     02/05/90  D.A.Saunders    Use of RLVEC to close a curve was bad - does
C                               not acknowledge line-type or thickness.  Chose
C                               to call CURVE a second time rather than require
C                               a duplicate of the first point.  This could
C                               impact use of DISSPLA's legend scheme - use
C                               LEGEND from this library instead.
C-------------------------------------------------------------------------------


C     Declarations.
C     -------------

      IMPLICIT NONE

C     Arguments.

      INTEGER
     >   IER, LINE, N, SYMBOL, THICK
      REAL
     >   COLOR, X (N), Y (N)

C     Local Constants.

      REAL
     >   FACTOR, XTHK, TLENG
      PARAMETER
     >  (FACTOR = 1.0005,  ! Used to separate coincident points for DISSPLA. 
     >   XTHK = .01,       ! Line thickness multiplier in inches.
     >   TLENG = 0.2)      ! Overall line pattern length in inches.

C     Local variables.  

      INTEGER
     >   ICOLOR, LTYP, TYPE 
      REAL
     >   AINT, HUE, RATCHDSH (4), RATCHDOT (4), RATLODSH (2), SAT, X2,Y2
      CHARACTER
     >   COLORAY (0:1) * 5

C     Storage.

      DATA AINT, SAT /1., 1./              ! Fix intensity and saturation for
                                           ! colors between magenta and blue.
      DATA COLORAY   /'WHITE', 'BLACK'/

      DATA RATCHDSH  /1.7, 0.3, 0.7, 0.3/  ! Set line patterns.
      DATA RATCHDOT  /2.3, 0.3, 0.1, 0.3/
      DATA RATLODSH  /0.7, 0.3/

      SAVE
     >   AINT, COLORAY, RATCHDSH, RATCHDOT, RATLODSH, SAT


C     Execution.
C     ----------

      IER = 0

C     Simple-minded error checking.  These are not necessarily fatal, so
C     just return with error flag set.

      IF (N .LE. 0) THEN

C        No plot data?

         IER = 1
      ELSE IF (LINE .EQ. 2 .AND. SYMBOL .LT. 0) THEN

C        Invisible line!?

         IER = 2
      END IF

      IF (IER .NE. 0) GO TO 999

C     Preset to defaults.  Note that 'WHITE' actually shows up as black
C     in hardcopy.  Resetting 'DOT' takes care of all the patterns.

      CALL RESET ('SETCLR')
      CALL RESET ('THKCRV')
      CALL RESET ('DOT')


C     Set line color.
C     ---------------

      IF (COLOR .GE. 0.0 .AND. COLOR .LE. 7.0) THEN
         IF (COLOR .GE. 2.0) THEN   ! Assign HSI colors from magenta to blue.
            HUE = .5 * (COLOR - 1.0)
            CALL HWHSI (HUE, SAT, AINT)
         ELSE                       ! Handle black or white separately.
            ICOLOR = COLOR
            CALL SETCLR (COLORAY (ICOLOR))
         END IF
      END IF

C     Set line pattern.
C     -----------------

      IF (THICK .NE. 0) CALL THKCRV (THICK * XTHK)

      LTYP = ABS (LINE) 
      IF (LTYP .EQ. 4) THEN
         CALL DOT
      ELSE IF (LTYP .EQ. 5) THEN
         CALL DASH
      ELSE IF (LTYP .EQ. 6) THEN

C        Chaindot.  Don't use DISSPLA's CHNDOT (too long a pattern).
C        Define a "custom" line type - see DISSPLA manual, Part B, Section 9.
C        Alternating marks and spaces are prescribed by the ratios of the
C        segments to one another and by their total length.

         CALL MRSCOD (3.0 * TLENG, 4, RATCHDOT)

      ELSE IF (LTYP .EQ. 7) THEN

C        Chaindash.  Don't use DISSPLA's CHNDSH (too long a pattern).

         CALL MRSCOD (3.0 * TLENG, 4, RATCHDSH)

      ELSE IF (LTYP .EQ. 8) THEN   ! Longdash.  No DISSPLA utility available.

         CALL MRSCOD (TLENG, 2, RATLODSH)

      ELSE IF (LTYP .EQ. 9) THEN   ! Use default width for thick line.

         CALL THKCRV (2. * XTHK)

      END IF

C     Set point marker type.
C     ----------------------

C     Note that MARKER will look for a custom symbol if called with a
C     negative argument, so we have to test first.

      IF (SYMBOL .GE. 0) THEN 

         CALL MARKER (SYMBOL)

C        Set DISSPLA's TYPE flag for symbols and/or lines.

         IF (LTYP .NE. 2) THEN  ! Line plus symbols.
            TYPE = +1
         ELSE                   ! Symbols alone.
            TYPE = -1
         END IF
      ELSE                      ! Line only
         TYPE = 0
      END IF


C     Plot the curve, and close it if appropriate.
C     --------------------------------------------

C     DISSPLA is known to fail if the curve to be plotted is more
C     than a single thickness and of zero length.  Guard against
C     this by moving these points apart by some satisfactory distance.
C     Note:  Handle the N=2 case only - no attempt here to handle three or
C     more coincident points, since this involves calculating arc length.

      X2 = X (2)
      IF (N .EQ. 2) THEN
         IF (X (1) .EQ. X2 .AND. Y (1) .EQ. Y (2)) THEN
            X (2) = X2 * FACTOR
            IF (X2 .EQ. 0.) X(2) = 1.E-6
         END IF
      END IF

      CALL CURVE (X, Y, N, TYPE)

      IF (LINE .LT. 0 .AND. LTYP .NE. 2) THEN        ! Don't use RLVEC - it
CCCCC    CALL RLVEC ( X (N), Y (N), X (1), Y (1), 0) ! ignores thickness/type
         X (2) = X (N)
         Y2 = Y (2)
         Y (2) = Y (N)
         CALL CURVE (X, Y, 2, TYPE)     ! Preferable to requiring (N+1)th pt.,
         Y (2) = Y2                     ! but could affect DISSPLA legends -
      END IF                            ! use LEGEND from this library instead.
         
      X (2) = X2

C     Reset everything.
C     -----------------

      CALL RESET ('SETCLR')
      CALL RESET ('THKCRV')
      CALL RESET ('DOT')


C     Termination.
C     ------------

  999 RETURN
      END
C+----------------------------------------------------------------------
C
      SUBROUTINE PROTECT (NX, X, Y, ARROW, DISTINCT)
C
C     One-liner: Test data for monotonicity and distinctness
C     ----------
C
C     Description and usage:
C     ----------------------
C
C        PROTECT scans a curve represented by arrays X and Y to check
C     strict monotonicity (increasing or decreasing) of the X data and
C     to verify that no two successive points match. One or the other
C     of these tests is required by conventional or parametric spline
C     fitting methods. Since it is often inefficient to perform these
C     tests in the fitting routine, which may be called more than once
C     with the same data or with data known to be good, the present
C     modular approach was chosen.
C
C        The initial application is program SMOOTH, where a single set
C     of data points may be fit by several different techniques - the
C     idea was to check each dataset just once as it was read in, and
C     then test flags ARROW and DISTINCT as needed.
C
C     Arguments:
C     ----------
C
C     Name     Type/Dimension  I/O/S  Description
C     NX       I               I      Dimension of X and Y arrays. If
C                                     NX <= 1, then ARROW will be
C                                     0.0 and DISTINCT will be .FALSE.
C
C     X        R (NX)          I      Array of abscissas.
C
C     Y        R (NX)          I      Array of ordinates.
C
C     ARROW    R                 O    Monotonicity indicator:
C                                       -1.0  strictly decreasing
C                                        0.0  neither
C                                       +1.0  strictly increasing
C
C     DISTINCT L                 O    Indicates whether successive points
C                                     are distinct in both X and Y.
C
C     Notes:
C     ------
C
C     (1)  IMPLICIT NONE and 8-character symbolic names are non-standard.
C
C     (2)  There is no provision for roundoff error in the monotonicity
C          test, i.e., all interval lengths are compared to zero.
C
C     Author:  Robert Kennelly, Sterling Federal Systems/NASA-Ames
C     -------
C
C     Development history:
C     --------------------
C
C     22 Feb. 1987    RAK    Initial design and coding.
C     14 Apr. 1988    RAK    Corrected value returned by DISTINCT when
C                            NX <= 1. Revised comments.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Constants.

      REAL
     &   ZERO, ONE
      PARAMETER
     &  (ZERO   = 0.0E+0,
     &   ONE    = 1.0E+0)

C     Arguments.

      LOGICAL
     &   DISTINCT
      INTEGER
     &   NX
      REAL
     &   ARROW, X (NX), Y (NX)

C     Local variables.

      LOGICAL
     &   FIRST
      INTEGER
     &   J
      REAL
     &   DX

C     Execution.
C     ----------

C     Set the output flags and bail out if the input data is trivial
C     or improper.

      ARROW    = ZERO
      DISTINCT = .FALSE.
      IF (NX .LE. 1) GO TO 990

C     Reset the direction flag according to the first interval, and reset
C     distinctness flag prior to testing.

      DX = X (2) - X (1)
      IF (DX .NE. ZERO) THEN
         FIRST = (DX .GT. ZERO)
         ARROW = SIGN (ONE, DX)
      END IF

      DISTINCT = .TRUE.

C     The approach is to try to set ARROW = ZERO and DISTINCT = .FALSE.
C     as we go, and quit as soon as possible. No harm is done if ARROW is
C     set repeatedly - it's not worth testing for.

      DO 10, J = 1, NX - 1
         DX = X (J + 1) - X (J)
         IF (DX .NE. ZERO) THEN

C           Compare the sign of the increment in this interval to that of
C           the first interval.

            IF ((DX .GT. ZERO) .NEQV. FIRST) ARROW = ZERO
         ELSE

C           If a pair of X's match, the data is not strictly monotonic. We
C           must still check whether the Y's are distinct.

            ARROW = ZERO
            IF (Y (J + 1) - Y (J) .EQ. ZERO) THEN

C              The data is neither monotonic nor distinct - time to die.

               DISTINCT = .FALSE.
               GO TO 990

            END IF
         END  IF
   10 CONTINUE

C     Termination.
C     ------------

  990 CONTINUE
      RETURN
      END
C+------------------------------------------------------------------------------
C
      SUBROUTINE QUICK
     >  (FRAME, NUMTIT, TITLES, LENTIT, MAXPT, CURVES, POINTS, X, Y,
     >   XMIN, XMAX, XSTEP, YMIN, YMAX, YSTEP, Y2LABEL, Y2SCALE,
     >   Y2SHIFT, HIGH, WIDE, ORIENT, PLOT, STYLE, FORMDAT, GRIDAT,
     >   FITDAT, SPLDAT, LINDAT, SYMDAT, COLDAT, LEGS, LEGTXT, LENLEG,
     >   LEGINF, NUMNOTE, NOTEBOX, NOTES, XNOTES, YNOTES, XARROWS,
     >   YARROWS, INTEGERX, INTEGERY, ESCAPE, LDEBUG)
C
C
C     Description and usage:
C
C           This routine produces simple plots using DISSPLA.  For stand-
C        alone use, QUICK is called by PROGRAM QPLOT, but it may also be
C        used as a general purpose utility in other applications.  (Since
C        all defaulting/error checking/etc. is assumed to have already been
C        done, the present version is not as appropriate for this as in the
C        past.)  Many curves may be plotted on one set of axes - the data
C        are packed in arrays X, Y.  The number of points to be plotted is
C        unlimited.  A number of options regarding minimum and maximum values,
C        step sizes, plot dimensions and orientation, color, line types, and
C        background grid are available.  Several fitting options are available
C        for parametric or nonparametric spline interpolation.
C
C           The plot will be centered on a medium 11" high (portrait) or 8.5"
C        high (landscape) with unlimited width (e.g. Versatec or roll-fed pen
C        plotter).  If the width parameter is small enough, the plot will be
C        centered horizontally on a single sheet of paper.  Axis numbering
C        is controlled by the MIN, MAX, and STEP inputs.  QPLOT uses
C        LAYOUT/LINAX/etc. for determination of "nice" linear values, and
C        LOGAX for log plots.  Polar diagrams are set up by DISSPLA.
C
C           The title, subtitle, and legend are centered at the top of the
C        page, and the axis labels are automatically centered on their axes.
C        The legend, if present, will be in one- or two-column format
C        depending on the number of entries and the widths of the two halves.
C        Caption(s) if present, will be left-justified below the plot.  The
C        plotted data points may be marked with symbols, colored, and/or
C        connected by one of several line types.
C
C           Plotting must be initiated and terminated from the calling
C        routine as follows:
C
C               :
C              CALL PLOTDEV ( ... )                  <-- Plot device
C               :                                        initialization.
C               :
C              CALL LAYOUT ( ... )                   <-- Helps determine
C                                                        good axis params.
C
C              CALL QUICK (NTITLES, ... )            <-- May be repeated.
C               :                                        Each call produces
C               :                                        one plot frame.
C              CALL DONEPL
C               :
C
C     Arguments:
C
C        Name     Dimension   Type   I/O/S   Description
C
C        FRAME                 I     I       Number of plot frame (page) for
C                                            use in diagnostic messages
C
C        NUMTIT                I     I       Number of title/subtitle/x-label/
C                                            y-label/caption strings.  The
C                                            special "classification" strings
C                                            are assumed to be in positions
C                                            NUMTIT+1 and NUMTIT+2.
C
C        TITLES (NUMTIT+2)*(*) C     I       Plot title, subtitle, axis labels,
C                                            captions, etc.  The size and style
C                                            of the lettering may be changed
C                                            using embedded DISSPLA commands,
C                                            for which the escape characters are
C                                            given by ESCAPE (below).
C                                            LENTIT (I) = 0 suppresses the
C                                            Ith text string.
C
C        LENTIT (NUMTIT+2)     I     I       Array containing number of the
C                                            characters in the TITLES.  (Pass
C                                            100 to invoke DISSPLA self counting
C                                            option.)
C
C        MAXPT                 I     I       Length of the packed data arrays
C                                            as declared in calling program.
C                                            (Helps decide whether X,Y arrays
C                                            can be used as workspace.)
C
C        CURVES                I     I       Number of curves to be plotted.
C
C        POINTS   CURVES       I     I       Array containing number of points
C                                            to be plotted for each curve.
C
C        X        MAXPT        R     I/O     Packed array of abscissas.
C                                            (Angle THETA for POLAR plots.)
C                                            NOTE: Scaled in place for fitted
C                                            curves if there is insufficient
C                                            space at end of array, then
C                                            restored - expect small changes.
C
C        Y        MAXPT        R     I/O     Packed array of ordinates.
C                                            (RADIUS for POLAR plots.) See
C                                            NOTE above.
C
C        XMIN,                 R     I       Axis limits corresponding to
C          XMAX,                             min-X, max-X, min-Y, and max-Y
C          YMIN,                             (Only YMIN, the minimum radius,
C          YMAX                              is meaningful for POLAR plots,
C                                            and an axis will only be drawn
C                                            if a label is requested.)  XMIN
C                                            or YMIN = "Origin" (in data units)
C                                            in log axes.
C
C        XSTEP,                R     I       Interval lengths on X- and Y-axes.
C          YSTEP                             (Used for CYCLE lengths in inches
C                                            with LOG plots; only YSTEP=B radial
C                                            units/inch is meaningful for
C                                            POLAR plots.)
C
C        Y2LABEL             C*(*)   I       Label for second Y axis if any.
C
C        Y2SCALE,              R     I       Scale and shift factors to apply
C          Y2SHIFT                           to the left hand tick values to
C                                            give the right hand tick values.
C                                            If both are zero, or if the scale is
C                                            1. and the shift is 0., the second Y
C                                            axis is suppressed.
C
C        HIGH,                 R     I       Nominal plot dimensions in inches.
C          WIDE
C
C        ORIENT     *          C     I       Plot orientation flag determines
C                                            which edge of a standard page
C                                            the X-axis lies along:
C                                              'LANDSCAPE' = 11" edge;
C                                              'PORTRAIT' = 8-1/2" edge.
C
C        PLOT       *          C     I       Plot type descriptor:
C                                              'LINEAR' = linear plot;
C                                              'LOGLOG' = log vs. log;
C                                              'LOGX' = linear vs. log;
C                                              'LOGY' = log vs. linear;
C                                              'POLAR' = THETA vs. R.
C
C        STYLE      *          C     I       Default lettering style for this
C                                            frame.  All the DISSPLA options
C                                            are available, including shaded
C                                            characters. See MYSTYL, and the
C                                            DISSPLA manual, for details. Use
C                                            string 'DISSPLA' for the default
C                                            (simple stick characters).
C
C        FORMDAT               I     I       Plot format:
C                                            +2   Tick marks all around the
C                                                 plot.
C                                            +1   Box with tick marks extending
C                                                 from lower left corner along
C                                                 both axes (Default).
C                                             0   No box or axes.
C                                            -1   Axes only (crossing at lower
C                                                 left hand corner).
C                                            -2   Force axis crossing to be
C                                                 at or near zero.
C
C        GRIDAT                I     I       Number of grid lines per scale
C                                            division if > 0.  If < 0, number of
C                                            tick marks per scale division.
C                                            (Note: Guarding against GRID = 0
C                                            is left to the application.)
C
C        FITDAT   CURVES       I     I       Type of fit used to interpolate
C                                            the data points:
C                                              1 = linear;
C                                              2 = monotone (tight);
C                                              3 = bessel (loose);
C                                              4 = closed (same as loose with
C                                                  wrap-around end conditions).
C                                            The calling routine is expected
C                                            to have checked for wrap-around
C                                            in case 4.  The fit will be
C                                            skipped if the number of points
C                                            is less than 2 (cases 2 and 3) or
C                                            3 (case 4).
C
C        SPLDAT   CURVES       I     I       Type of spline used in conjunction
C                                            with FITDAT = 2 or 3:
C                                              1 = standard (Y vs. X);
C                                              2 = reversed (X vs. Y);
C                                              3 = parametric (X & Y vs. arc);
C                                            The spline will be suppressed if
C                                            the abscissas are not monotonically
C                                            increasing or decreasing.
C
C        LINDAT   CURVES       I     I       Line type codes for each curve.
C                                            See POLYLINE for the definitions.
C                                            Not used explicitly here, except in
C                                            one error check for symbols only.
C
C        SYMDAT   CURVES       I     I       Symbol type codes (applying only
C                                            to line types which use symbols).
C                                            See POLYLINE for the definitions.
C                                            Not used explicitly here, except
C                                            in one test for negativity, which
C                                            means no symbols.
C
C        COLDAT   CURVES       I     I       Color codes for each curve.
C                                            Not used explicitly here - see
C                                            POLYLINE for the definitions.
C
C        LEGS                  I     I       Number of legend entries.  The
C                                            application imposes some limit.
C                                            LEGS=0 is OK.
C
C        LEGTXT   (LEGS)*(*)   C     I       Legend text entries, if any.
C                                            No trailing '$'s needed unless
C                                            LENLEG (*) = 100.
C
C        LENLEG     (*)        I     I       Array containing number of
C                                            characters in legend line(s).
C                                            (Pass 100 to invoke DISSPLA self
C                                            counting option.)
C
C        LEGINF   (3, LEGS)    I     I       Legend info. apart from the text:
C                                            (1, *) = line type code;
C                                            (2, *) = symbol code;
C                                            (3, *) = color code.
C                                            Not used explicitly here - see
C                                            POLYLINE for the definitions.
C
C        NUMNOTE               I     I       (Max.) number of annotated points.
C
C        NOTEBOX               I     I       Controls the boxes around any
C                                            annotation strings in NOTES (*).
C                                            If boxes are drawn, they allow
C                                            a border of .05" around the text,
C                                            which is .12" high.
C                                            0 means no box frames and no
C                                              blanking;
C                                            1 means no frames but the boxes are
C                                              blanked;
C                                            2 means the boxes are framed and
C                                              blanked;
C                                           >2 thickens the frames: a factor of
C                                               NOTEBOX - 1 is applied to the
C                                               pen width.
C
C        NOTES  (NUMNOTE)*(*)  C     I       Text for annotated points.
C                                            May be blank for arrow-only.
C
C        XNOTES,  (NUMNOTE)    R     I       Coordinates (in data units) at
C          YNOTES                            which to place (non-blank)
C                                            annotation strings. Bottom left
C                                            corners of text go here (or
C                                            tails of arrows if text is blank).
C
C        XARROWS, (NUMNOTE)    R     I       Coordinates (in data units) of the
C          YARROWS                           tips of arrows drawn from the
C                                            CENTER of the text blocks defined
C                                            by the corresponding elements of
C                                            XNOTES, YNOTES, and NOTES.  Each
C                                            arrow is clipped by a .05" border
C                                            around the text.
C
C                                            N.B.: Set YARROWS (I) = YNOTES (I)
C                                            if a precisely HORIZONTAL arrow is
C                                            desired (whether NOTES (I) is blank
C                                            or not).  Likewise, set XARROWS (I)
C                                            = XNOTES (I) if an exactly VERTICAL
C                                            arrow is desired.  It will be
C                                            realigned with the center of the
C                                            text box.
C
C        INTEGERX,             L     I       Axis numbering styles.
C          INTEGERY                          .TRUE. means no decimal points;
C                                            .FALSE. means real numbers.
C
C        ESCAPE               C*2    I       Escape characters for switching
C                                            to DISSPLA's "instruction" alpha-
C                                            bet and back.  E.g.: '[]'
C
C        LDEBUG                I     I       Logical unit number for error
C                                            messages from DISSPLA or curve
C                                            fit routine PLSFIT. Also controls
C                                            whether the page-edge job
C                                            information string is printed.
C                                            The possibilities are:
C                                               >0: ordinary unit number;
C                                               =0: suppress all output;
C                                               <0: print only fatal messages
C                                                   on unit ABS (LDEBUG).
C
C
C     External references:
C
C        CHORD    Summed chord lengths for X-Y curve over range of indices.
C        CLASSIF  Writes special headings on plots.
C        LCSFIT   Local cubic spline utility (Y vs. X).
C        LEGEND   Portable legend utility (not DISSPLA's).
C        MERGER   Data-merging utility.
C        MYSTYL   Sets default lettering style.
C        PLSFIT   Storage-efficient parametric local cubic spline fit.
C        POLYLINE Connect-the-data-points utility.
C        PROTECT  Test data for monotonicity and distinctness.
C        RANGER   Prepare for transformation of data to interval [0, 1].
C        RESCALE  Linear transformation of data packed in an array.
C        SCAN2    Used to find length of annotation text lines.
C        WINDOW   Find a segment of curve (X, Y) enclosed by a rectangle.
C
C     Environment:  Digital VAX-11/780 VMS Version 4.1 (FORTRAN 77)
C                   Computer Associates CA-DISSPLA Version 11.0-9003
C
C     NOTES:
C
C        (1)  IMPLICIT NONE, 8-character symbolic names, and "!" as comment
C             character are not (yet) standard. CALL EXIT (3) is used to
C             inform the operating system of abnormal termination (VAXish).
C
C        (2)  The titles or legends may be blank.  This is best achieved
C             by entering zero for their lengths.
C
C        (3)  See the DISSPLA manual for the order of the plotting symbols
C             used.  The first few are: Square, Octagon, Triangle ("up"),
C             +, X, Diamond, and Triangle ("down").
C
C        (4)  A large variety of type sizes and styles are available. 
C
C             See Vol. 1, Ch. 6, Sec. 11 of the DISSPLA manual for details on
C             the 'INSTRUCTION' alphabetic commands. The escape to the special
C             command set is ESCAPE (1:1) (e.g. '[') and the return is
C             ESCAPE (2:2) (e.g. ']').  These characters therefore cannot
C             be plotted in a text string.  Taking them to be '[' and ']',
C             for instance, the string '[M6]a[M0]' will create a Greek alpha
C             by switching to the 'MATHEMATIC' character set, invoking the
C             appropriate character, then returning to the 'STANDARD' alphabet.
C             Alternative fonts are similarly specified, e.g. [F1] for
C             CARTOGraphic and [F3] for SCMPLX (the QUICK default). (Remember
C             QUICK's STYLE argument is a global control though.) Other
C             commands of interest are:  E for superscripts, L for subscripts,
C             and Hr for multiplying the character height by factor r. There
C             are LOTS more.  See the available tables for details.
C
C        (5)  Drastic re-sizing (using HIGH and WIDE) may require changes in
C             character sizes or spacing. This may be done by embedding the
C             appropriate directives in the input strings, or more generally,
C             by modifying the calls to HEIGHT in this routine.
C
C        (6)  User-defined line types and legend handling are simpler now,
C             at the expense of the application's having to update the legend 
C             entries explicitly (as opposed to implicit updating with each call
C             to CURVE).
C
C        (7)  The parametric interpolation schemes are local, requiring no
C             significant storage except for a buffer array to hold the
C             results, so there is no limit to the number of points which
C             can be interpolated (except for max. no. of points handled
C             by QUICK). We use a simple-minded approach to filling in the
C             interpolated points: a fixed number of points are distributed
C             uniformly along the part of the curve which appears to lie
C             within the plot window.
C
C
C     Author:  Robert Kennelly, Sterling Software/NASA Ames, Moffett Field, CA
C
C
C     Development history:
C
C        14 May  1982    RAK    Original design and coding.
C        29 Sep. 1982    RAK    Modified to plot multiple curves.
C        18 Dec. 1982    RAK    Repaired scaling loop indexing.
C        17 Feb. 1983    RAK    Modified to permit axis limit and line
C                               type inputs.
C        28 Feb. 1983    RAK    Generalized TYPE input.
C        31 Mar. 1983    RAK    SCALE is now left unchanged.
C         3 Nov. 1983    RAK    Shrink the plot (the better to make
C                               viewgraphs), eliminate border, add
C                               edge labels, change font.
C        13 Feb. 1984    RAK    Extensive revisions for legends, grids,
C                               various line types. Use packed data arrays.
C                               Layout dimensions changed to PARAMETERs.
C                               Input text is now CHARACTER type.
C        23 Feb. 1984    RAK    Added subtitle option using STORY.
C         2 Mar. 1984    RAK    Tied CROSS to GRID. Title position depends
C                               on number of legend entries.
C        13 Mar. 1984    RAK    Added 'INSTRUCTION' alphabet, input SIZE,
C                               modified legend format, reversed tick marks,
C                               moved labels and modified layout details,
C                               used default physical origin, suppressed
C                               marginal run data, packed legend data. Whew!
C        23 Apr. 1984    RAK    Reinstated marginal run data, with running
C                               frame numbers.
C        11 May  1984    RGL    Added logarithmic Y scale vs. linear X
C                               plot option; diagnostic now written onto
C                               the otherwise blank plot area when name-
C                               list option PLOT incorrectly designated.
C        18 June 1984    RAK    Moved error checking on PLOT to higher level.
C                               Added protection for axis limits.
C        15 Aug. 1984    RGL    For convergence plots, guarded against
C                               too many cycles for plotting on the
C                               logarithmic Y axis and calculated an
C                               interval size along the X axis in order
C                               to obtain true integer labelling.
C        21 Sep. 1984    LJC    Removed limits on WIDTH option for plots
C                               more than one page wide.
C        27 Sep. 1984    RGL    Revised the two-column legend layout to
C                               center both columns together as if one
C                               extra-long legend over the entire plot,
C                               in favor of the previous centering of
C                               each column over a separate half of plot.
C        23 Jan. 1985    LJC    Added 'SCALE' (equal scaling) to PLOT options.
C        11 Apr. 1985    RAK    Use passed-in TOTAL rather than MAXPT. Use
C                               XMIN, XMAX, etc. rather than SCALE(I), and
C                               drop SIZE(I) in favor of HIGH, WIDE.
C        21 Aug. 1985    RAK    Try out an ORIENTation input. Page size,
C                               and hence physical origin as set by AREA2D,
C                               depend on orientation now. Must now
C                               set PLOT = 'LINEAR' from calling routine -
C                               defaulting no longer done at this level.
C        11 Sep. 1985    RAK    Use 'CONNECT' for default line type, set by
C                               QPLOT using dictionary routine. Revert to
C                               old method of centering two-column legends.
C                               (Previous method, 27 Sep. 1984, was buggy.)
C        13 Sep. 1985    RAK    Confine legend to center of plot, even if
C                               WIDE is larger than a page (see COLWID).
C        27 Sep. 1985    RAK    Support for log plots (X- and/or Y-axis).
C                               PLOT type string must be spelled out. Drop
C                               the 'CONVERGENCE' plot option (use 'LOGY'
C                               for almost the same result). Added new
C                               line type 'LONGDASHES' ("hand-made"). Drop
C                               'SCALE' type (at this level) since it is
C                               now controlled by calling routine. New
C                               'POLAR' plot type, with radial origin YMIN.
C         3 Oct. 1985    RAK    Refined PACK dimension estimate. Only look
C                               at first MAXLEG curves (we couldn't really
C                               label curves no. 51, and up, before). Pass
C                               an array of (INTEGER) symbol type flags.
C                               Line type array is now INTEGER. Use special
C                               SPCMOD/MYSPEC to save line types (required
C                               for custom line types or color).
C         9 Oct. 1985    RAK    Protect YMIN in LOGLOG since it may really
C                               apply only to a linear scale to be spliced
C                               in later (e.g., LOGX plots).
C        21 Oct. 1985    RAK    Re-repaired the previous bug: we actually
C                               have to handle each case separately, because
C                               GRAF and LOGLOG are too fragile.
C        25 Aug. 1986    RCL    Added four title lines to be used to draw
C                               captions below the plot.
C        16 Oct. 1986    RAK    Distinguish between GRIDAT = -1, which
C                               now means cross axes in lower left but
C                               don't draw the rest of the box, and -2
C                               which means cross axes at origin (old -1).
C                               Added color. MYSPEC now (re)sets line type
C                               for each curve.
C        22 Oct. 1986    RAK    Added call to SHDCHR to enable filled,
C                               shaded fonts. Input parameter STYLE permits
C                               the calling routine to specify the default
C                               lettering style for the whole frame, set
C                               in routine MYSTYL.
C         6 Mar. 1987    RAK    Added fitting options. To accommodate
C                               changes in QPLOT & LAYOUT, translate
C                               from PORTRAIT/LANDSCAPE terminology to
C                               DISSPLA's peculiar COMIC/MOVIE terms.
C         8 Apr. 1987    RAK    Use WINDOW to fit each visible curve
C                               segment separately. For now, just use
C                               fixed number of points per segment, at
C                               at equal intervals of arclength (kludge).
C                               Check for longest legend to aid layout.
C                               Added DEBUG to optionally suppress plot
C                               summary and "out-of-range" warnings.
C        17 Apr. 1987    RAK    Re-arrange the curve drawing so that
C                               dummy curves are drawn before LEGEND,
C                               followed later by both linear and fitted
C                               curves. Use LDEBUG for logical unit no.
C                               Pass dummy values to GRAF when X- or
C                               Y-axis will be overridden anyway. Permit
C                               a user-requested symbol along with a
C                               special line type. The DUMMY points must
C                               be positive so that (dummy) LOG plots work.
C                               Tie POLAR axis drawing to existence of a
C                               requested axis label. Increase NFIT to 500.
C         8 May  1987    RAK    Draw both POLAR axes if either is labeled.
C                               Add some defensive RESETs for line type.
C        13 May  1987    RAK    Reduce magnitude of X/YDUMMY to protect
C                               fragile CURVE routine (NOT MY FAULT!).
C        12 Oct. 1987    RAK    Use RESET ('SETCLR') instead of switching
C                               to 'WHITE' - neater.
C        31 Mar. 1988    RAK    Rescale data fed to PLSFIT to reduce
C                               chance of overflow. Added PROTECT for
C                               final check of data to be plotted. Pass
C                               in MAXPT, compute TOTAL so that any extra
C                               workspace can be used while fitting. Use
C                               logical unit LDEBUG to issue warning if
C                               a curve fit had to be skipped. Choose
C                               X/YDUMMY based on data, but plot off the
C                               page by shifting origin (avoids overflow
C                               in CURVE).
C        17 Aug. 1988    RAK    Repaired handling of XMIN, XMAX, etc. Was
C                               performing LOG inside loop, which could
C                               result in LOG (-ve) on second fitted curve.
C        28 Dec. 1988    RAK    Fixed not one, but two typos in the test
C                               for distinctness before data is pushed
C                               upstairs for rescaling! Was comparing Y
C                               to X, and also requiring both X's and Y's
C                               to be different. Must have been a Monday.
C        20 Mar. 1989  M.D.Wong Added optional classification headings to
C                               plot.  Made minor adjustments to spacing of
C                               captions and X axis label to accomodate new
C                               features.
C        27 Mar. 1989    MDW    Adapted to use SMDLIB in place of DISSPLA.
C                               (Includes simplified legend handling.)
C           June 1989    MDW    Ideas from DISSPLA and SMDLIB versions merged
C                               into new DISSPLA version.  This includes use
C                               of new LEGEND (not DISSPLA's) and POLYLINE.
C           Aug. 1989    MDW    Updated to run on DISSPLA version 11.0.
C                               Took out all %REF's.  Substituted BLANK for
C                               DOLLAR in some subroutine calls to enable
C                               printing of tick labels.
C        17 Jan. 1990    MDW    Enabled drawing of tick marks around the plot
C                               by re-thinking GRIDDAT and adding FORMDAT to
C                               specify plot format.
C        20 Feb. 1990    MDW    TITLES = ' ' now forces a blank to be printed,
C                               and TITLES (2) = 'NONE' suppresses subtitle.
C                               Replaced use of STORY to print title and
C                               subtitle with MESSAG.  Number of characters in
C                               text strings passed through argument list to
C                               avoid reliance on DISSPLA's self-counting.
C         4 Mar. 1991  DAS/RAK  Raised legend text height from .08" to .09"
C                               to suit PostScript output better.
C        25 Mar. 1991    "      Made positioning of the halves of a split legend
C                               more dynamic: first, try to center each half
C                               over half the plot, but with not too big a
C                               space between; otherwise, treat as a unit with
C                               minimum space between and with grace margins of
C                               .75" & .25" left and right.  Any overflow will
C                               be to the right.  Vertical gap between entries
C                               is now 0.9*normal for 8 or more entries/column.
C                               CALL LINESP (1.8) has been removed (not used).
C                               Caption positioning is now less obscure.
C                               NOTE: All text positioning is in NOMINAL inches.
C        22 Sep. 1991    "      Rethought positioning of the plot based on the
C                               dimensions of all of its components (where no
C                               explicit call to PHYSOR was made before);
C                               allowed for any number of caption lines;
C                               introduced second Y axis and annotations.
C        26 Sep. 1991    "      Suppressed width of (missing) line segment in
C                               full columns of symbols-only legend entries;
C                               improved 1-/2-column legend decision-making,
C                               with no hard limit at this level on how many
C                               legend entries there may be.  There must be
C                               at least 4 entries before two columns are
C                               considered, though.
C        27 Sep. 1991  DAS/DBS  Introduced optional arrows for the annotations
C                               (which may be blank).
C        04 Oct. 1991    DAS    LENTIT (*) = 0 is legal anywhere.  (No need
C                               for 'NONE' in subtitle to suppress it.)
C                               FORMDAT = 0 suppresses axes from polar plots
C                               as for other plots.
C        26 Oct. 1991    DAS    Introduced ESCAPE and INTEGERX/Y arguments.
C        25 Nov. 1991  DAS/RAK  Rob points out that a binding margin of 0.4"
C                               moves the typical default-size plot too far
C                               to the right.  (The eye tends to see the box,
C                               not the Y-axis label.)  Reduced binding margins
C                               to zero in all cases as a compromise.
C        29 Apr. 1992    DAS    Y2SCALE default should be 1, not 0.  Rather than
C                               change the user guide, force 1 if 0 is entered.
C                               (It's used to scale YSTEP, which mustn't be 0.)
C        15 Sep. 1997    DAS    Y2SCALE = 0. should be changed to 1. BEFORE the
C                               test, which should be against 1., not 0.
C        31 July 2002    DAS    Added SPLINE argument for nonparametric options;
C                               local FRAME counter did not match that of the
C                               calling program if the plot is previewed before
C                               being sent to a metafile - make it an argument.
C-------------------------------------------------------------------------------


      IMPLICIT NONE

C     Arguments.
C     ----------

C     (CURVES must be declared ahead of COLDAT.)

      INTEGER
     >   FRAME, CURVES, COLDAT (CURVES), LENTIT (*), LENLEG (*),
     >   FITDAT (CURVES), SPLDAT (CURVES), FORMDAT, GRIDAT, LDEBUG,
     >   LEGS, LINDAT (CURVES), LEGINF (3, *), MAXPT, NOTEBOX, NUMNOTE,
     >   NUMTIT, POINTS (CURVES), SYMDAT (CURVES)
      REAL
     >   HIGH, WIDE, X (MAXPT), XMAX, XMIN, XARROWS (*), XNOTES (*),
     >   XSTEP, Y2SCALE, Y2SHIFT, Y (MAXPT), YARROWS (*), YMAX, YMIN,
     >   YNOTES (*), YSTEP
      LOGICAL
     >   INTEGERX, INTEGERY
      CHARACTER
     >   ESCAPE * 2, LEGTXT (*) * (*), NOTES (*) * (*), ORIENT * (*),
     >   PLOT * (*), STYLE * (*), TITLES (*) * (*), Y2LABEL * (*)

C     Local constants.
C     ----------------

      INTEGER
     >   MAXCUR, NFIT
      PARAMETER
     >  (MAXCUR = 200,           ! Max. number of curves
     >   NFIT   = 500)           ! Number of fitted points/curve segment

      REAL
     >   FLAG, GAP, GAPMAX, GAPMIN, GLEGND, GRACEL, HALF, HCAPTN,
     >   HLEGND, HSTITL, HTITL, OFFSET, ONE, SEG, SGAP, SYM, TENTH, TWO,
     >   ZERO

      PARAMETER
     >  (FLAG   = 999.,          ! "Undefined" value used by ANNOTATE to tell
                                 ! if an arrow has been requested.
     >   GAP    = 0.25,          ! Legend/box/title space; see X-lab./captn. too
     >   GAPMAX = 0.50,          ! Max. gap between two halves of 2-col. legend
     >   GAPMIN = 0.25,          ! Min.  "     "
     >   GLEGND = 0.09,          ! Gap between entries = HLEGND unless 8 or more
     >   GRACEL = 1.E+0,         ! X grace margin for legend, split 3/4 & 1/4
     >   HALF   = 0.5E+0,
     >   HCAPTN = 0.12,
     >   HLEGND = 0.09,
     >   HSTITL = 0.12,
     >   HTITL  = 0.15,
     >   OFFSET = 0.0,
     >   ONE    = 1.0E+0,
     >   SEG    = 0.54,          ! Length of line segment in legend entries
     >   SGAP   = 0.20,          ! Gap between line segment and legend text
     >   SYM    = 0.09,          ! Empirical width of symbols
     >   TENTH  = 0.1,
     >   TWO    = 2.0E+0,
     >   ZERO   = 0.0E+0)

C     NOTE: SEG = 0.54 rather than 0.60 suppresses the last space of
C     ----- the chaindot/chaindash/longdash patterns in the legend.
C           (3 * 0.2" for the whole pattern; 0.3 * 0.2" for the gaps.)

      CHARACTER
     >   BESSEL * 1, BLANK * 1, MONOTONE * 1
      PARAMETER
     >  (BESSEL = 'B', BLANK  = ' ', MONOTONE = 'M')

C     Local variables.
C     ----------------

      INTEGER
     >   BASE, BASETEMP, FIRST, I, IBEGIN, IEND, IER, ITEMP, J,
     >   JJ, LAST, LENNOTE, LENY2, MARK, NLEG1, NSORTED, TOTAL
      REAL
     >   ARROW, CORNER, DELTA, DELTAY, DERIVS (NFIT), FRACTION,
     >   GUTTERX, GUTTERY, HALFW,
     >   HEIGHTCA, HEIGHTLE, HEIGHTPA, HEIGHTPL, HEIGHTTI, HEIGHTXL,
     >   ORIGINX, ORIGINY, SPACE, SUBGAP, TBEGIN, TEND, WIDTH1, WIDTH2,
     >   WIDTHLE, WIDTHPA, WIDTHPL, WIDTHSE, WIDTHSE1, WIDTHSE2,
     >   WIDTHY1, WIDTHY2,
     >   XCAPTN, XFIT (NFIT), XLEG1, XLEG2, XMAXWIND, XMINWIND, XSCALE,
     >   XSHIFT,
     >   YCAPTN, YFIT (NFIT), YLEG1, YLEG2, YMAXWIND, YMINWIND, YSCALE,
     >   YSHIFT
      LOGICAL
     >   CLOSED, DISTINCT, NEW, PUSHDATA, SECONDY, TWOCOL
      CHARACTER
     >   DISTRIB * 8, METHOD * 1

C     Procedures.
C     -----------

      REAL
     >   CHORD, XMESS, XDIMTB
      EXTERNAL
     >   CHORD, CLASSIF, LCSFIT, LEGEND, MERGER, MYSTYL, PLSFIT,
     >   POLYLINE, PROTECT, RANGER, RESCALE, SCAN2, WINDOW, XDIMTB,
     >   XMESS

C     Execution.
C     ----------

C     Basic DISSPLA page setup.

      CALL RESET ('ALL')

      IF (LDEBUG .LE. 0) THEN
         CALL NOCHEK
         CALL SETDEV (0, 0)
      ELSE
         CALL SETDEV (ABS (LDEBUG), ABS (LDEBUG))
      END IF

      IF (ORIENT .EQ. 'LANDSCAPE') THEN
         CALL HWROT ('MOVIE')
         WIDTHPA  = 11.0     ! Modified below
         HEIGHTPA = 8.5
         GUTTERX  = ZERO
         GUTTERY  = 0.0      ! Binding margin (top) (disabled now)
      ELSE                   ! Portrait mode
         CALL HWROT ('COMIC')
         WIDTHPA  = 8.5      ! Modified below
         HEIGHTPA = 11.0
         GUTTERX  = 0.0      ! Binding margin (left) (disabled)
         GUTTERY  = 0.0      ! Could fudge a bit here (-0.1?) for extra at bot.
      END IF

      WIDTHPA = MAX (WIDTHPA, WIDTHPL + GUTTERX)  ! Adjusts only wide plots

      CALL HWSCAL ('NONE')
      CALL NOBRDR

C     Not clear how to ensure seeing DISSPLA job information on the edge of
C     the page.  (It's sometimes in the plot but off the page.)

      CALL PAGE (WIDTHPA, HEIGHTPA)

C     Select alphabets and style (before XDIMTB is called as part of
C     dynamically positioning the plot).

      CALL MX1ALF ('STANDARD', ESCAPE (2:2))
      CALL MX2ALF ('INSTRUCTION', ESCAPE (1:1))
      CALL SHDCHR (90.0, 1, .005, 1)

      IF (STYLE .NE. 'DISSPLA') CALL MYSTYL (STYLE)

C     Some more initialization.

C     Y2SCALE = 1. would've been a better default, but we need to guard
C     against 0. anyway since it's used to scale YSTEP, which can't be 0.

      IF (Y2SCALE .EQ. ZERO) Y2SCALE = ONE
      SECONDY = Y2SCALE .NE. ONE .OR. Y2SHIFT .NE. ZERO

      IF (SECONDY) THEN
         FIRST = 1         ! Find no. of characters in secondary Y label, if any
         LENY2 = LEN (Y2LABEL)
         CALL SCAN2 (Y2LABEL, BLANK, FIRST, LENY2, MARK)
      ELSE
         LENY2 = 0
      END IF


C     Determine the dimensions in inches of all plot components.
C     ----------------------------------------------------------

C     First, the title(s).

      HEIGHTTI = ZERO
      IF (LENTIT (1) .GT. 0) HEIGHTTI = HTITL
      IF (LENTIT (2) .GT. 0) HEIGHTTI = HEIGHTTI + 1.8 * HSTITL

C     Now the legend, if any.

      IF (LEGS .GT. 0) THEN

C        Do some of the figuring for the 1- or 2-column legend question.
C        Use one column if there are only a few entries, or if first and
C        last half of the entries are too long to fit side by side.

         NLEG1 = (LEGS + 1) / 2

C        Kludge: PostScript option fails in XDIMTB unless at level 2 or 3

         CALL AREA2D (WIDE, HIGH)      ! Moves to level 2
         CALL HEIGHT (HLEGND)

         WIDTH1 = XDIMTB (LEGTXT, LENLEG, 1, NLEG1)
         WIDTH2 = XDIMTB (LEGTXT, LENLEG, NLEG1 + 1, LEGS)

         CALL RESET ('HEIGHT')
         CALL ENDGR (-1)

C        Check for suppressing line segments if columns are symbols-only.

         WIDTHSE1 = SYM
         DO 2, I = 1, NLEG1
            IF (LEGINF (1, I) .NE. 2) WIDTHSE1 = SEG
    2    CONTINUE

         WIDTHSE2 = SYM
         DO 4, I = NLEG1 + 1, LEGS
            IF (LEGINF (1, I) .NE. 2) WIDTHSE2 = SEG
    4    CONTINUE

         WIDTHSE = MAX (WIDTHSE1, WIDTHSE2) ! In case we go to 1 column
         WIDTHLE = MAX (WIDTH1, WIDTH2) + WIDTHSE + SGAP
         WIDTH1  = WIDTH1 + WIDTHSE1 + SGAP
         WIDTH2  = WIDTH2 + WIDTHSE2 + SGAP

         TWOCOL = (WIDTH1 + GAPMIN + WIDTH2 .LE. WIDE + GRACEL) .AND.
     >            (LEGS .GT. 3)
         I = LEGS
         IF (TWOCOL) I = NLEG1
         DELTAY = GLEGND
         IF (I .GE. 8) DELTAY = 0.9 * GLEGND
         HEIGHTLE = REAL (I) * (HLEGND + DELTAY) - DELTAY + GAP + GAP

      ELSE  ! No legend: leave somewhat less than 2*GAP between plot and titles.

         HEIGHTLE = GAP + HALF * GAP

      END IF   

C     Now the X axis labeling.

      IF (LENTIT (3) .EQ. 0) THEN  ! A bit more than the numerical labels
         HEIGHTXL = 1.5 * GAP      ! in case there is a caption
      ELSE
         HEIGHTXL = 2.3 * GAP      ! Lower than DISSPLA's default
      END IF

C     Now the captions.

      ITEMP = NUMTIT - 4
      IF (ITEMP .EQ. 0) THEN
         HEIGHTCA = ZERO
      ELSE
         HEIGHTCA = GAP + REAL (ITEMP) * (1.8 * HCAPTN) - 0.8 * HCAPTN
      END IF

C     Now the width of the Y axis label(s).

      IF (LENTIT (4) .EQ. 0) THEN
         WIDTHY1 = 0.4     ! Empirical width of DISSPLA's numerical Y labels
      ELSE
         WIDTHY1 = 0.765   ! 0.625 + 0.14 (Label offset + height of label)
      END IF

      IF (SECONDY) THEN
         IF (LENY2 .EQ. 0) THEN
            WIDTHY2 = 0.4  ! As for main Y axis
         ELSE
            WIDTHY2 = 0.765
         END IF
      ELSE
         WIDTHY2 = ZERO
      END IF

C     Now we know the plot size:

      WIDTHPL = WIDTHY1 + WIDE + WIDTHY2
      HEIGHTPL = HEIGHTTI + HEIGHTLE + HIGH + HEIGHTXL + HEIGHTCA

C     Finally, the plot origin can be placed explicitly (lower left corner).

      ORIGINX = (( WIDTHPA - GUTTERX) - WIDTHPL ) * HALF + GUTTERX +
     >          WIDTHY1
      ORIGINY = ((HEIGHTPA - GUTTERY) - HEIGHTPL) * HALF +
     >          HEIGHTCA + HEIGHTXL

      CALL PHYSOR (ORIGINX, ORIGINY)
      CALL AREA2D (WIDE, HIGH)


C     Tick marks?  Grid?
C     ------------------

      IF (GRIDAT .LE. -1) THEN
         CALL XREVTK   ! Turn the tick marks inwards since there is no grid.
         CALL YREVTK
         CALL XTICKS (ABS (GRIDAT))   ! Display ticks. 
         CALL YTICKS (ABS (GRIDAT))
      ELSE
         CALL XTICKS (0)   ! We don't need tick marks if there is a grid.
         CALL YTICKS (0)
      END IF

C     Y-axis numbering is horizontal.

      CALL YAXANG (ZERO)

C     Integer numbering if possible?

      IF (INTEGERX) CALL XINTAX
      IF (INTEGERY) CALL YINTAX

C     Force axes to cross at zero, if appropriate.

      IF (FORMDAT .EQ. -2) CALL CROSS

C     Enable labeling and drawing of both axes, except as modified below.

      IF (FORMDAT .NE. 0 .AND. PLOT .NE. 'POLAR') THEN
         CALL XNAME (BLANK, 1)
         CALL YNAME (BLANK, 1)
      END IF


C     Draw the axes.
C     --------------

      IF (PLOT .EQ. 'LINEAR') THEN

C        Linear X- and Y-axes.

C        Draw any secondary axes first, because if they're rescaled
C        they affect the data drawn later.

         IF (FORMDAT .EQ. 2 .OR. SECONDY) THEN

            CALL XNONUM
            CALL YNONUM

C           Kludge: must bring DISSPLA to level 3 via dummy call to GRAF.

            CALL GRAF (ZERO, ONE, ONE, ZERO, ONE, ONE)

C           Draw an opposing unlabeled linear X-axis with ticks below it.

            CALL RESET ('XREVTK')
            CALL XGRAXS (XMIN, XSTEP, XMAX, WIDE, BLANK, 1, 
     >                   OFFSET, OFFSET + HIGH)     
            CALL XREVTK
            CALL RESET ('XNONUM')

            IF (.NOT. SECONDY) THEN

C              Draw an opposing unlabeled linear Y-axis with ticks to the left.

               CALL RESET ('YREVTK')
               CALL YGRAXS (YMIN, YSTEP, YMAX, HIGH, BLANK, 1, 
     >                      OFFSET + WIDE, OFFSET)     
               CALL YREVTK
               CALL RESET ('YNONUM')

            ELSE

C              Draw a secondary Y-axis with the same tick marks but different
C              tick values determined by user-input scale and shift.
C              (Example: units = ft. on primary axis; metres on secondary axis)

               CALL RESET ('YNONUM')
               CALL YGRAXS (YMIN * Y2SCALE + Y2SHIFT, YSTEP * Y2SCALE,
     >                      YMAX * Y2SCALE + Y2SHIFT, HIGH, Y2LABEL,
     >                      -LENY2, OFFSET + WIDE, OFFSET)     

            END IF

            CALL ENDGR (-1)             ! Can't call next GRAF from level 3
            CALL AREA2D (WIDE, HIGH)    ! Move back from level 1 to level 2
         END IF

C        Primary axes:

         CALL GRAF (XMIN, XSTEP, XMAX, YMIN, YSTEP, YMAX)

      ELSE IF (PLOT .EQ. 'LOGX') THEN

C        Start with dummy linear X-axis, then override it.

         CALL XNAME (BLANK, 0)
CCC      CALL GRAF (ZERO, 'SCALE', ONE, YMIN, YSTEP, YMAX)  ! Why 'SCALE'?
         CALL GRAF (ZERO,  TENTH,  ONE, YMIN, YSTEP, YMAX)

         CALL XNAME (BLANK, +1)
         CALL XLGAXS (XMIN, XSTEP, WIDE, BLANK, 1, OFFSET, OFFSET)

         IF (FORMDAT .EQ. 2) THEN

C           Add opposing log X-axis.

            CALL XNONUM
            CALL RESET ('XREVTK')
            CALL XLGAXS (XMIN, XSTEP, WIDE, BLANK, 1, OFFSET, 
     >                   OFFSET + HIGH)

C           Then add opposing linear Y-axis.

            CALL YNONUM
            CALL RESET ('YREVTK')
            CALL YGRAXS (YMIN, YSTEP, YMAX, HIGH, BLANK, 1, OFFSET + 
     >                   WIDE, OFFSET)     
         END IF

      ELSE IF (PLOT .EQ. 'LOGY') THEN

C        Start with dummy linear Y-axis, then override it.

         CALL YNAME (BLANK, 0)
CCC      CALL GRAF (XMIN, XSTEP, XMAX, ZERO, 'SCALE', ONE)  ! Why 'SCALE'?
         CALL GRAF (XMIN, XSTEP, XMAX, ZERO,  TENTH,  ONE)

         CALL YNAME (BLANK, +1)
         CALL YLGAXS (YMIN, YSTEP, HIGH, BLANK, 1, OFFSET, OFFSET)

         IF (FORMDAT .EQ. 2) THEN

C           Add opposing log Y-axis.

            CALL YNONUM
            CALL RESET ('YREVTK')
            CALL YLGAXS (YMIN, YSTEP, HIGH, BLANK, 1, OFFSET + WIDE, 
     >                   OFFSET)

C           Then add opposing linear X-axis.

            CALL XNONUM
            CALL RESET ('XREVTK')
            CALL XGRAXS (XMIN, XSTEP, XMAX, WIDE, BLANK, 1, OFFSET, 
     >                   OFFSET + HIGH)     
         END IF

      ELSE IF (PLOT .EQ. 'LOGLOG') THEN

C        Both axes are logarithmic.

         CALL LOGLOG (XMIN, XSTEP, YMIN, YSTEP)
    
         IF (FORMDAT .EQ. 2) THEN

C           Draw an opposing set of unlabeled logarithmic axes.

            CALL XNONUM
            CALL YNONUM
            CALL RESET ('XREVTK')
            CALL RESET ('YREVTK')
            CALL XLGAXS (XMIN, XSTEP, WIDE, BLANK, 1, OFFSET, 
     >                   OFFSET + HIGH) 
            CALL YLGAXS (YMIN, YSTEP, HIGH, BLANK, 1, OFFSET + 
     >                   WIDE, OFFSET)  
         END IF

      ELSE IF (PLOT .EQ. 'POLAR') THEN

C        Polar diagram. Since the axes (with ticks and numbering) can
C        get in the way of the data, only draw them if at least one axis
C        label was explicitly requested.

         IF (FORMDAT .NE. 0 .AND.
     >       LENTIT (3) .GT. 0 .OR. LENTIT (4) .GT. 0) THEN
             CALL XNAME (TITLES (3), LENTIT (3))
             CALL YNAME (TITLES (4), LENTIT (4))
         END IF

         CALL POLORG (YMIN)
         CALL POLAR (ONE, YSTEP, HALF * WIDE, HALF * HIGH)

      END IF


C     Label the primary axes explicitly.
C     ----------------------------------

      IF (PLOT .NE. 'POLAR') THEN

C        The X-axis label will be put in by brute force because the default
C        positioning chosen by XNAME is unaesthetically close to the axis
C        numbering.

         IF (LENTIT (3) .GT. 0) THEN
            CORNER = HALF * (WIDE - XMESS (TITLES (3), LENTIT (3)))
            CALL MESSAG (TITLES (3), LENTIT (3), CORNER, -HEIGHTXL)
         END IF

C        The label for the Y-axis is handled similarly so that the same
C        generous string length permitted by MESSAG can be used for both axes.
C        (YNAME has a limit of 60 characters (VAX) or 72 (Cray).)

         IF (LENTIT (4) .GT. 0) THEN
            CALL ANGLE (90.0)
            CORNER = HALF * (HIGH - XMESS (TITLES (4), LENTIT (4)))
            CALL MESSAG (TITLES (4), LENTIT (4), -0.625, CORNER)
            CALL RESET ('ANGLE')
         END IF
      END IF


C     Annotations should be drawn early in case their areas are blanked.
C     ------------------------------------------------------------------

      CALL HEIGHT (HSTITL)

      CALL ANNOTATE (NUMNOTE, NOTES, XNOTES, YNOTES, XARROWS, YARROWS,
     >               NOTEBOX, FLAG, ABS (LDEBUG))


C     Mark the grid (if any).
C     -----------------------

      IF (FORMDAT .GT. 0) THEN
         CALL DOT
         IF (GRIDAT .LE. -1) THEN
            CALL GRID (0, 0)
         ELSE
            CALL GRID (GRIDAT, GRIDAT)
         END IF
         CALL RESET ('DOT')
      END IF


C     Finish the legend set-up and draw it (one- or two-column format).
C     -----------------------------------------------------------------

      IF (LEGS .GT. 0) THEN

C        Expand the plotting area to allow for the drawing of legend
C        lines and symbols via POLYLINE/CURVE.

         CALL GRACE (3.0)
         CALL HEIGHT (HLEGND)
         YLEG1 = HIGH + GAP

C        Some of the legend placement figuring got moved upstairs as part
C        of the overall plot placement.  No need to redo it.

         IF (.NOT. TWOCOL) THEN   ! Single-column legend.

            XLEG1 = HALF * (WIDE - WIDTHLE)  ! Corner, in inches from origin
            IF (XLEG1 .LT. ZERO) THEN
               FRACTION = HALF
               IF (LENY2 .EQ. 0) FRACTION = 0.75
               XLEG1 = -FRACTION * MIN (-XLEG1, GRACEL)
            END IF

            CALL LEGENDZ (LEGTXT, LENLEG, 1, LEGS, DELTAY, XLEG1, YLEG1,
     >                    LEGINF, HLEGND, WIDTHSE, SGAP, .FALSE.)
         ELSE

C           Two columns with tops at the same level.

            YLEG2 = YLEG1 + REAL (MOD (LEGS, 2) * 2) * HLEGND
            HALFW = HALF * (WIDE - GAPMIN)

            IF (WIDTH1 .LE. HALFW .AND. WIDTH2 .LE. HALFW) THEN

C              Center the two halves above the two halves of the box ...

               XLEG1 = HALF * (HALFW - WIDTH1)
               XLEG2 = HALF * (HALFW - WIDTH2) + HALFW + GAPMIN

C              ... but don't allow too big a space between the columns.

               SPACE = XLEG2 - XLEG1 - WIDTH1
               SPACE = MAX (SPACE - GAPMAX, ZERO) * HALF
               XLEG1 = XLEG1 + SPACE
               XLEG2 = XLEG2 - SPACE

            ELSE                

C              One half is wider than half the box (or both halves are wider).
C              Try to center the two halves as one whole.

               XLEG1 = HALF * (WIDE - (WIDTH1 + GAPMIN + WIDTH2))

               IF (XLEG1 .LT. ZERO) THEN
                  FRACTION = HALF
                  IF (LENY2 .EQ. 0) FRACTION = 0.75
                  XLEG1 = -FRACTION * MIN (-XLEG1, GRACEL)
               END IF
               XLEG2 = XLEG1 + WIDTH1 + GAPMIN  ! Overflow will be to the right.
            END IF

            CALL LEGENDZ (LEGTXT, LENLEG, 1, NLEG1, DELTAY, XLEG1,
     >                    YLEG1, LEGINF, HLEGND, WIDTHSE1, SGAP,
     >                    .FALSE.)

            CALL LEGENDZ (LEGTXT, LENLEG, NLEG1 + 1, LEGS, DELTAY,
     >                    XLEG2, YLEG2, LEGINF, HLEGND, WIDTHSE2, SGAP,
     >                    .FALSE.)
         END IF
      END IF

C     We're done with the legend, so reset parameters.

      CALL GRACE (ZERO)
      CALL RESET ('DOT')
      CALL RESET ('THKCRV')
      CALL RESET ('SETCLR')


C     Plot the data.
C     --------------

      TOTAL = 0
      DO 6, I = 1, CURVES
         TOTAL = TOTAL + POINTS (I)
    6 CONTINUE

C     Compute transformed window edges which may be needed for fitted
C     curves in the loop below.

      IF (PLOT .EQ. 'LOGX' .OR. PLOT .EQ. 'LOGLOG') THEN
         XMINWIND = LOG (XMIN)    
         XMAXWIND = LOG (XMAX)
      ELSE
         XMINWIND = XMIN
         XMAXWIND = XMAX
      END IF

      IF (PLOT .EQ. 'LOGY' .OR. PLOT .EQ. 'LOGLOG') THEN
         YMINWIND = LOG (YMIN)    
         YMAXWIND = LOG (YMAX)
      ELSE
         YMINWIND = YMIN
         YMAXWIND = YMAX
      END IF

C     Loop over curves to interpolate and/or plot them.
C     Parametric splines are calculated in normalized space (~ plot box);
C     nonparametric splines are calculated in data space so as not to be
C     misleading when plotted.
C
C     The "TEMP" quantities are used during the parametric spline
C     portion of the loop to manipulate scaled/transformed data which
C     will usually have been moved to scratch space at the end of the
C     data arrays.

      BASE = 1

      DO 200, I = 1, CURVES

         BASETEMP = BASE
         ITEMP    = I

         IF (FITDAT (I) .EQ. 1 .OR. LINDAT (I) .EQ. 2) THEN

C           Straight line segments, or symbols only (no spline fit)
C           -------------------------------------------------------

            CALL POLYLINE (POINTS (I), X (BASE), Y (BASE), LINDAT (I),
     >                     0, SYMDAT (I), REAL (COLDAT (I)), IER)

         ELSE IF (SPLDAT (I) .EQ. 1) THEN ! Spline of Y vs. X

C           Nonparametric spline fit in data space if Xs are monotonic
C           ----------------------------------------------------------

            CALL PROTECT (POINTS (I), X (BASE), Y (BASE),
     >                    ARROW, DISTINCT)

            IF (ARROW .EQ. ZERO) THEN
               IF (LDEBUG .NE. 0) WRITE (ABS (LDEBUG), 1000)
     >            'abscissas', I, FRAME

               IF (SYMDAT (I) .LT. 0) THEN ! No symbols; force piecewise linear

                  CALL POLYLINE (POINTS (I), X (BASE), Y (BASE),
     >                           LINDAT (I), 0, SYMDAT (I),
     >                           REAL (COLDAT (I)), IER)
               END IF

               GO TO 20
            END IF

            IF (FITDAT (I) .EQ. 2) THEN
               METHOD = MONOTONE
            ELSE IF (FITDAT (I) .EQ. 3) THEN
               METHOD = BESSEL
            END IF

            JJ = BASE + POINTS (I) - 1
            DELTA = (X (JJ) - X (BASE)) / REAL (NFIT - 1)

            DO 10, J = 1, NFIT - 1
               XFIT (J) = X (BASE) + REAL (J-1) * DELTA
   10       CONTINUE
            XFIT (NFIT) = X (JJ)

            NEW = .TRUE.

C           Merge the data points into the evaluation abscissas if there's
C           room off the end of the data proper:

            IF (POINTS (I) + NFIT .GE. MAXPT - TOTAL) THEN ! Forget it

               CALL LCSFIT (POINTS (I), X (BASE), Y (BASE), NEW, METHOD,
     >                      NFIT, XFIT, YFIT, DERIVS)

               CALL POLYLINE (NFIT, XFIT, YFIT, LINDAT (I), 0, -1,
     >                        REAL (COLDAT (I)), IER)
            ELSE ! Merge

               BASETEMP = TOTAL + 1    ! Move the uniform fit Xs to where
               JJ = BASETEMP           ! there's room to insert the data Xs.
               DO 15, J = 1, NFIT
                  X (JJ) = XFIT (J)
                  JJ     = JJ + 1
   15          CONTINUE

               NSORTED = NFIT

               CALL MERGER (POINTS (I), X (BASE), NSORTED,
     >                      X (BASETEMP), ZERO)

               CALL LCSFIT (POINTS (I), X (BASE), Y (BASE), NEW, METHOD,
     >                      NSORTED, X (BASETEMP), Y (BASETEMP),
     >                      Y (BASETEMP))  ! No distinct derivs. array may
                                           ! cause trouble on some systems

               CALL POLYLINE (NSORTED, X (BASETEMP), Y (BASETEMP),
     >                        LINDAT (I), 0, -1, REAL (COLDAT (I)), IER)
            END IF

C           Draw in the original points as symbols on top of the fitted
C           curve, if required.

   20       IF (SYMDAT (I) .GE. 0) THEN
               CALL POLYLINE (POINTS (I), X (BASE), Y (BASE), 2, 0,
     >                        SYMDAT (I), REAL (COLDAT (I)), IER)
            END IF

         ELSE IF (SPLDAT (I) .EQ. 2) THEN ! Spline of X vs. Y

C           Nonparametric spline fit in data space if Ys are monotonic
C           ----------------------------------------------------------

            CALL PROTECT (POINTS (I), Y (BASE), X (BASE),
     >                    ARROW, DISTINCT)

            IF (ARROW .EQ. ZERO) THEN
               IF (LDEBUG .NE. 0) WRITE (ABS (LDEBUG), 1000)
     >            'ordinates', I, FRAME

               IF (SYMDAT (I) .LT. 0) THEN ! No symbols; force piecewise linear

                  CALL POLYLINE (POINTS (I), X (BASE), Y (BASE),
     >                           LINDAT (I), 0, SYMDAT (I),
     >                           REAL (COLDAT (I)), IER)
               END IF

               GO TO 40
            END IF

            IF (FITDAT (I) .EQ. 2) THEN
               METHOD = MONOTONE
            ELSE IF (FITDAT (I) .EQ. 3) THEN
               METHOD = BESSEL
            END IF

            JJ = BASE + POINTS (I) - 1
            DELTA = (Y (JJ) - Y (BASE)) / REAL (NFIT - 1)

            DO 30, J = 1, NFIT - 1
               YFIT (J) = Y (BASE) + REAL (J-1) * DELTA
   30       CONTINUE
            YFIT (NFIT) = Y (JJ)

            NEW = .TRUE.

C           Merge the data points into the evaluation abscissas if there's
C           room off the end of the data proper:

            IF (POINTS (I) + NFIT .GE. MAXPT - TOTAL) THEN ! Forget it

               CALL LCSFIT (POINTS (I), Y (BASE), X (BASE), NEW, METHOD,
     >                      NFIT, YFIT, XFIT, DERIVS)

               CALL POLYLINE (NFIT, XFIT, YFIT, LINDAT (I), 0, -1,
     >                        REAL (COLDAT (I)), IER)
            ELSE ! Merge

               BASETEMP = TOTAL + 1    ! Move the uniform fit Ys to where
               JJ = BASETEMP           ! there's room to insert the data Ys.
               DO 35, J = 1, NFIT
                  Y (JJ) = YFIT (J)
                  JJ     = JJ + 1
   35          CONTINUE

               NSORTED = NFIT

               CALL MERGER (POINTS (I), Y (BASE), NSORTED,
     >                      Y (BASETEMP), ZERO)

               CALL LCSFIT (POINTS (I), Y (BASE), X (BASE), NEW, METHOD,
     >                      NSORTED, Y (BASETEMP), X (BASETEMP),
     >                      X (BASETEMP))  ! No distinct derivs. array may
                                           ! cause trouble on some systems

               CALL POLYLINE (NSORTED, X (BASETEMP), Y (BASETEMP),
     >                        LINDAT (I), 0, -1, REAL (COLDAT (I)), IER)
            END IF

C           Draw in the original points as symbols on top of the fitted
C           curve, if required.

   40       IF (SYMDAT (I) .GE. 0) THEN
               CALL POLYLINE (POINTS (I), X (BASE), Y (BASE), 2, 0,
     >                        SYMDAT (I), REAL (COLDAT (I)), IER)
            END IF

         ELSE ! SPLDAT (I) = 3 (parametric spline)

C           Piecewise-cubic interpolated curve (parametric in normalized space)
C           -------------------------------------------------------------------

C           Push the data upstairs if there's room for it.

            PUSHDATA = (POINTS (I) .LE. MAXPT - TOTAL) .AND.
     >                 (CURVES .LT. MAXCUR)

            IF (PUSHDATA) THEN
               JJ = TOTAL
               DO 50, J = BASE, BASE + POINTS (I) - 2
                  IF (X (J + 1) .NE. X (J) .OR.
     >                Y (J + 1) .NE. Y (J)) THEN

C                    Move the distinct points up into scratch space (simply
C                    overwrite any data previously moved).

                     JJ     = JJ + 1
                     X (JJ) = X (J)
                     Y (JJ) = Y (J)
                  END IF
   50          CONTINUE
               JJ     = JJ + 1
               X (JJ) = X (J)
               Y (JJ) = Y (J)

C              Set up pointers to the "new" data.

               BASETEMP       = TOTAL + 1
               ITEMP          = CURVES + 1
               POINTS (ITEMP) = JJ - TOTAL
            END IF

C           If the scaling of the two axes is very different, it helps to
C           temporarily rescale the data before interpolating. (The problem
C           arises when a vertical curve segment, say, appears to have zero
C           arclength because the chords of the neighboring intervals are
C           dominated by the data along the horizontal axis. Think about it.)

            IF (PLOT .EQ. 'LOGX' .OR. PLOT .EQ. 'LOGLOG') THEN
               DO 60, J = BASETEMP, BASETEMP + POINTS (ITEMP) - 1
                  X (J) = LOG (X (J))
   60          CONTINUE
            END IF

            IF (PLOT .EQ. 'LOGY' .OR. PLOT .EQ. 'LOGLOG') THEN
               DO 70, J = BASETEMP, BASETEMP + POINTS (ITEMP) - 1
                  Y (J) = LOG (Y (J))
   70          CONTINUE
            END IF

C           Determine normalization factors for X and Y, and apply them:

            CALL RANGER  (MAXPT, POINTS, ITEMP, X, XSCALE, XSHIFT)
            CALL RESCALE (MAXPT, POINTS, ITEMP, X, XSCALE, XSHIFT)

            CALL RANGER  (MAXPT, POINTS, ITEMP, Y, YSCALE, YSHIFT)
            CALL RESCALE (MAXPT, POINTS, ITEMP, Y, YSCALE, YSHIFT)

C           Check whether the scaled data points are still sufficiently
C           distinct for parametric interpolation.

            CALL PROTECT (POINTS (ITEMP), X (BASETEMP), Y (BASETEMP),
     >                    ARROW, DISTINCT)

            IF (.NOT.DISTINCT) THEN
               IF (LDEBUG .NE. 0) WRITE (ABS (LDEBUG), 1005) I, FRAME

               IF (SYMDAT (I) .LT. 0) THEN ! No symbols; force piecewise linear

                  IF (PUSHDATA) THEN ! Original data has not been scaled

                     CALL POLYLINE (POINTS (I), X (BASE), Y (BASE),
     >                              LINDAT (I), 0, SYMDAT (I),
     >                              REAL (COLDAT (I)), IER)
                  END IF

               END IF

               GO TO 110
            END IF


C           Loop over the segments of each curve (if any) within the window.
C           ----------------------------------------------------------------

            IBEGIN = BASETEMP
            IEND   = BASETEMP + POINTS (ITEMP) - 1
   80       CONTINUE

C              Since we let DISSPLA handle POLAR plots, with (almost) no
C              control over data range, the segment-finding is only applied
C              to the other plot types. Data limits must be scaled to agree
C              with (temporarily) rescaled plot data.

               IF (PLOT .NE. 'POLAR') THEN
                  CALL WINDOW (BASETEMP + POINTS (ITEMP) - 1, X, Y,
     >                         XMINWIND * XSCALE + XSHIFT,
     >                         XMAXWIND * XSCALE + XSHIFT,
     >                         YMINWIND * YSCALE + YSHIFT,
     >                         YMAXWIND * YSCALE + YSHIFT,
     >                         IBEGIN, IEND, IER)

C                 In case there was not enough data, or if none of the points
C                 lie inside the plot window, we're done fitting this curve.

                  IF (IER .NE. 0) GO TO 110
               END IF

C              Convert index limits to arclength for PLSFIT.

               TBEGIN  = CHORD (X, Y, BASETEMP, IBEGIN)
               TEND    = TBEGIN + CHORD (X, Y, IBEGIN, IEND)

               NEW     = .TRUE.
               CLOSED  = .FALSE.
               DISTRIB = 'UNIFORM'

               IF (FITDAT (I) .EQ. 2) THEN
                  METHOD = MONOTONE
               ELSE IF (FITDAT (I) .EQ. 3) THEN
                  METHOD = BESSEL
               ELSE IF (FITDAT (I) .EQ. 4) THEN

C                 This is QPLOT's "closed" case, for geometrical shapes.

                  CLOSED = .TRUE.
                  METHOD = BESSEL
               END IF

C              Piecewise cubic local spline interpolation.

               CALL PLSFIT (POINTS (ITEMP), X (BASETEMP), Y (BASETEMP),
     >                      TBEGIN, TEND, NFIT, XFIT, YFIT, NEW, CLOSED,
     >                      METHOD, DISTRIB, IER)

               IF (IER .EQ. 0) THEN

C                 Apply the inverse of the original data scaling to the
C                 interpolated points if necessary, then plot them. (We're
C                 cheating and passing scalar NFIT as an array.)

                  CALL RESCALE (NFIT, NFIT, 1, XFIT, ONE / XSCALE,
     >                          -XSHIFT / XSCALE)
                  CALL RESCALE (NFIT, NFIT, 1, YFIT, ONE / YSCALE,
     >                          -YSHIFT / YSCALE)

                  IF (PLOT .EQ. 'LOGX' .OR. PLOT .EQ. 'LOGLOG') THEN
                     DO 90, J = 1, NFIT
                        XFIT (J) = EXP (XFIT (J))
   90                CONTINUE
                  END IF

                  IF (PLOT .EQ. 'LOGY' .OR. PLOT .EQ. 'LOGLOG') THEN
                     DO 100, J = 1, NFIT
                        YFIT (J) = EXP (YFIT (J))
  100                CONTINUE
                  END IF

                  CALL POLYLINE (NFIT, XFIT, YFIT, LINDAT (I), 0, -1,
     >                           REAL (COLDAT (I)), IER)

               ELSE IF (IER .EQ. +1) THEN

C                 Insufficient number of data points - just skip the plotting.

                  CONTINUE
               ELSE                              ! IF (IER .GT. +1) THEN

C                 Somebody goofed - hit the brakes. All of the inputs to
C                 PLSFIT should have been verified prior to the call.

                  IF (LDEBUG .NE. 0) WRITE (ABS (LDEBUG), 1010) I, FRAME
                  CALL EXIT (3)
               END IF

C              Update pointers and loop back to look for another segment.
C              (We merely restore IEND, since it was changed by WINDOW.)

               IBEGIN = IEND
               IEND   = BASETEMP + POINTS (ITEMP) - 1
               IF (IBEGIN .LT. IEND) GO TO 80  ! Next segment.

  110       CONTINUE  ! End of loop over line segments.

C           Restore original data scaling if it was done in place.

            IF (.NOT. PUSHDATA) THEN
               CALL RESCALE (TOTAL, POINTS, I, X, ONE / XSCALE,
     >                       -XSHIFT / XSCALE)
               CALL RESCALE (TOTAL, POINTS, I, Y, ONE / YSCALE,
     >                       -YSHIFT / YSCALE)

               IF (PLOT .EQ. 'LOGX' .OR. PLOT .EQ. 'LOGLOG') THEN
                  DO 120, J = BASE, BASE + POINTS (I) - 1
                     X (J) = EXP (X (J))
  120             CONTINUE
               END IF

               IF (PLOT .EQ. 'LOGY' .OR. PLOT .EQ. 'LOGLOG') THEN
                  DO 130, J = BASE, BASE + POINTS (I) - 1
                     Y (J) = EXP (Y (J))
  130             CONTINUE
               END IF
            END IF

C           Draw in the original points as symbols on top of the fitted
C           curve, if required.

            IF (SYMDAT (I) .GE. 0) THEN
               CALL POLYLINE (POINTS (I), X (BASE), Y (BASE), 2, 0,
     >                        SYMDAT (I), REAL (COLDAT (I)), IER)
            END IF

         END IF

C        Update the index to the beginning of the next curve, then loop
C        back to plot it.

         BASE = BASE + POINTS (I)

  200 CONTINUE                                   ! End loop over curves.


C     Display title(s), caption(s), and classification headings.
C     ----------------------------------------------------------

C     Write optional subtitle below the main title.

      IF (LENTIT (2) .GT. 0) THEN
         CALL HEIGHT (HSTITL)
         CORNER = HALF * (WIDE - XMESS (TITLES (2), LENTIT (2)))
         CALL MESSAG (TITLES (2), LENTIT(2), CORNER, HIGH + HEIGHTLE)
         SUBGAP = 1.8 * HSTITL
      ELSE
         SUBGAP = ZERO
      END IF

      CALL HEIGHT (HTITL)
      CORNER = HALF * (WIDE - XMESS (TITLES (1), LENTIT (1)))
      CALL MESSAG (TITLES (1), LENTIT(1), CORNER,
     >             HIGH + HEIGHTLE + SUBGAP)

C     Captions:  MESSAG is used instead of LSTORY so that text may indented.

      CALL HEIGHT (HCAPTN)
      XCAPTN = ZERO
      YCAPTN = -(HEIGHTXL + GAP + HCAPTN)   ! See X-label positioning
      DO 300, I = 5, NUMTIT
         CALL MESSAG (TITLES (I), LENTIT (I), XCAPTN, YCAPTN)
         YCAPTN = YCAPTN - 1.8 * HCAPTN
  300 CONTINUE

C     Classification headings, if any, are in fixed locations top and bottom.

      CALL CLASSIF (ORIENT, TITLES (NUMTIT+1), TITLES (NUMTIT+2), 
     >              LENTIT (NUMTIT+1))

C     Termination.
C     ------------

      IF (LDEBUG .GT. 0) THEN

C        Close the frame and add edge-label information string.

         CALL ENDPL (FRAME)
      ELSE
         CALL ENDPL (0)
      END IF

      RETURN


C     Formats.
C     --------

 1000 FORMAT (/' QUICK:  Warning - ', A, ' are not monotonically ',
     >   'increasing or decreasing'/
     >   9X, 'for curve number ' , I3, ' in frame number ', I3, '.'/
     >   9X, 'Spline fit downgraded to piecewise linear.')
 1005 FORMAT (/' QUICK:  Warning - data points are insufficiently ',
     >   'distinct for parametric spline fit'/
     >   9X, 'for curve number ' , I3, ' in frame number ', I3, '.'/
     >   9X, 'Interpolation will be skipped.')
 1010 FORMAT (/' QUICK:  Abnormal termination - programmer error ',
     >   'in call to PLSFIT for'/
     >   9X, 'curve number ' , I3, ' in frame number ', I3, '.')

      END
C+----------------------------------------------------------------------
C
      SUBROUTINE RANGER (MAXPT, NPT, CURVE, X, SCALE, SHIFT)
C
C
C     One-liner: Prepare for transformation of data array to interval [0, 1].
C
C
C     Description and usage:
C
C           Prepare for a linear transformation of data stored in a packed
C        array. The upper bound (XUB) and lower bound (XLB) of the data are
C        used to determine the parameters required to map the data to the
C        interval [0, 1], with protection against all array elements being
C        identical. Since the data array may be large, we check whether
C        anything really needs to be done before proceeding. The eventual
C        transformation is assumed to have the form:
C
C                    X (output) = SCALE * X (input) + SHIFT
C
C        Hard STOPs are provided if MAXPT or the NPT array are unreasonable.
C
C           RANGER was written for use with RESCALE in QPLOT, a general-
C        purpose data graphics program developed at NASA-Ames. It could
C        find other applications, but may not pay unless the data is 
C        actually packed.
C
C
C     Arguments:
C
C        Name    Dimension  Type  I/O/S  Description
C        MAXPT              I     I      Declared size of data array in
C                                        calling routine (required for
C                                        error checking).
C
C        NPT     CURVE      I     I      Point counts for each curve. Simple
C                                        tests are performed to verify that
C                                        the computed data limits within
C                                        the X array are sensible, i.e.
C                                        > 0 and <= MAXPT.
C
C        CURVE              I     I      Index of curve to be checked.
C
C        X       MAXPT      R     I      Array to be examined.
C
C        SCALE              R       O    Scale factor: 1 / (XUB - XLB).
C                                        If XUB = XLB, then SCALE is set
C                                        to 1 / XUB unless XUB = 0, in
C                                        which case SCALE = 1.
C
C        SHIFT              R       O    Offset: -XLB / (XUB - XLB).
C                                        If XUB = XLB, then SHIFT is set
C                                        to 0 unless XUB = 0, in which
C                                        case SHIFT = 1.
C
C
C     Environment:  Digital VAX-11/780 VMS FORTRAN 77
C
C
C     Notes:
C
C        (1)  IMPLICIT NONE and 8-character symbolic names are non-standard.
C
C        (2)  The protection for cases where all data points are identical
C             is such that SCALE and SHIFT are chosen to map the data to
C             [1, 1] when XLB (lower bound) = XUB (upper bound).
C
C
C     Author:  Robert Kennelly, NASA-Ames/Sterling Software
C
C
C     Development history:
C
C        31 Mar. 1988    RAK    Initial design and coding.
C
C-----------------------------------------------------------------------


C     Declarations.
C     -------------

      IMPLICIT NONE

C     Constants.

      REAL
     >   ZERO, ONE
      PARAMETER
     >  (ZERO = 0.0E+0,
     >   ONE  = 1.0E+0)

C     Arguments.

      INTEGER
     >   CURVE, MAXPT, NPT (CURVE)
      REAL
     >   SCALE, SHIFT, X (MAXPT)

C     Local variables.

      INTEGER
     >   FIRST, I, LAST
      REAL
     >   XLB, XUB


C     Execution.
C     ----------

      IF (MAXPT .LE. 0) STOP 'RANGER: Array dimension error!'

C     Count the points in each curve up to this one.

      LAST = 0
      DO 10, I = 1, CURVE
         LAST = LAST + NPT (I)
   10 CONTINUE

      IF (LAST .LE. 0)     STOP 'RANGER: No data!'
      IF (LAST .GT. MAXPT) STOP 'RANGER: (Apparently) too much data!'

C     Set index for the begining of this curve's data points.

      FIRST = LAST - NPT (CURVE) + 1

C     Scan the data for upper and lower bound.

      XUB = X (FIRST)
      XLB = X (FIRST)

      DO 20, I = FIRST, LAST
         XUB = MAX (X (I), XUB)
         XLB = MIN (X (I), XLB)
   20 CONTINUE

C     Compute transformation parameters.

      IF (XUB .NE. XLB) THEN
         SCALE = ONE / (XUB - XLB)
         SHIFT = -XLB * SCALE
      ELSE                                      ! IF (XUB .EQ. XLB) THEN

C        Fudge the transformation so we will map to [1, 1].

         IF (XUB .NE. ZERO) THEN
            SCALE = ONE / XUB
            SHIFT = ZERO
         ELSE
            SCALE = ONE
            SHIFT = ONE
         END IF
      END IF

C     Termination.
C     ------------

  990 CONTINUE
      RETURN
      END
C+----------------------------------------------------------------------
C
      SUBROUTINE READCOLS (LREAD, LWRITE, XCOL, YCOL, MAXPT, POINTS,
     >   MAXCUR, CURVE, X, Y, BUFFERS, LAST, FINIS, ERROR)
C
C     One-liner:   READ two COLumnS of numerical data
C
C     Description:
C
C           READCOLS reads one set of X, Y pairs until either a curve-ending
C        keyword or EOF is encountered.  To read several sets of data, this
C        routine must be called from within a loop which should terminate
C        normally when flag FINIS is .TRUE. (signalling end of file for the
C        current data file).  The other return conditions are many - see the
C        ERROR argument description and "Notes" below.
C
C           READCOLS was originally written for use by QPLOT, but it may
C        also be useful in other applications involving multi-column files.
C
C           This version returns whenever the first token on a line is
C        found to be non-numeric (presumably a keyword).  All handling of
C        keywords is left to the calling program.  (END CURVE and END
C        FRAME were originally handled in READCOLS, but this was eliminated
C        with a view to implementing a keyword control scheme as opposed to
C        the original namelist scheme of QPLOT.)  This version also avoids
C        the backspacing that namelist handling had demanded.
C
C
C     Arguments:
C
C        Name    Dimension  Type  I/O/S  Description
C        LREAD               I    I      Logical unit number for input.
C        LWRITE              I    I      Logical unit number for output.
C                                        Negative suppresses warnings but
C                                        fatal errors are written to unit
C                                        ABS (LWRITE) unless LWRITE = 0
C                                        in which case NO output is written.
C        XCOL                I    I      Input column number from which X
C                                        data is to be read.  If XCOL=0,
C                                        then X(I) = I for I=1,2,3,....
C        YCOL                I    I/O    Input column number from which Y
C                                        data is to be read.  If YCOL=0,
C                                        then Y(I) = I for I=1,2,3,....
C        MAXPT               I    I      Maximum total number of points to
C                                        be read.
C        POINTS  MAXCUR      I    I/O    Number of points in each curve.
C                                        This array indexes the X and Y
C                                        data arrays.
C        MAXCUR              I    I      Maximum number of curves.
C        CURVE               I    I/O    Current number of curves.  Enter 0 for
C                                        first call.
C        X       MAXPT       R    I      Packed array of abscissas.
C        Y       MAXPT       R    I      Packed array of ordinates.
C        BUFFERS   2       C*(*)  I/O    IF LAST > 0, BUFFERS (1) (1:LAST)
C                                        contains the last data record read
C                                        (both on input and on output).
C                                        On output BUFFERS (2) stores formatted
C                                        point and curve information, which may
C                                        be used for upper level error handling.
C        LAST                I    I/O    See BUFFERS.  LAST = 0 means the
C                                        buffer is not meaningful.
C        FINIS               L      O    End-of-input flag (legal end of file).
C        ERROR               I      O    Error flag which indicates an input 
C                                        error was encountered. It assumes a
C                                        value as follows:
C                                        ERROR = 0  Normal (including FINIS=T).
C                                               -1  Non-correctable error
C                                                   reported (some software
C                                                   limit has been exceeded).
C                                               +1  Keyword caused termination 
C                                                   of data read.
C                                               +2  Invalid numeric data read.
C
C     Error Handling:
C
C           A message is printed on unit LWRITE for any error condition
C        (if LWRITE > 0), and on unit ABS (LWRITE) for fatal errors.
C        This permits the calling program to suppress warning messages.
C        Processing halts for READ errors judged fatal.  See LWRITE and
C        ERROR descriptions above, and NOTES below.
C
C
C     External References:
C
C        ALPHA   Identifies non-numeric characters.
C        EXIT    Stop and return status flag to system (VAX/VMS FORTRAN).
C        IOCHEK  Returns LOGICAL flags describing READ status.
C        GETLINE Reads one line, suppressing trailing blanks/comments
C        TOKENS  Finds tokens in text string and returns them in an array.
C
C
C     External Files:
C
C        Unit    I/O/S  Description
C        LREAD   I      Input file 
C        LWRITE    O    Output file ("line printer").
C
C
C     Environment:  Digital VAX-11/780 VMS FORTRAN (FORTRAN 77).
C
C
C     Notes:
C
C        (1)  IMPLICIT NONE is non-standard.
C
C        (2)  Some symbols are up to eight characters long.
C
C        (3)  CALL EXIT (N) is VAX-specific.  It causes program termination,
C             closes all files, and returns control to the operating system.
C             Its argument is then available to the system.  N = 3 is defined
C             here to mean abnormal termination due to fatal read error.
C
C        (4)  Since this routine is intended to be embedded within a loop
C             over sets of data (several frames each consisting of several
C             curves, in plotting terminology), there are two sets of logical
C             control flags floating around.  It is important to keep them
C             straight:
C
C                FINIS  Indicates a normal EOF for the current data file.
C                ERROR  Tells the calling program to use what it has (if
C                       CURVE > 0) and then quit if ERROR .NE. 0 or +1.
C                DONE   Indicates when the internal loop over points is
C                       complete. (Local)
C                FATAL  General bail-out for unrecoverable read errors.
C                       (Local)
C
C        (5)  Some of the details concerning the data packing are handled
C             at this level in an attempt to simplify the calling routine
C             (cf. BASE, MAXCUR).
C
C
C     History:
C
C        31 Dec. 1982    RAK    Initial design and coding.
C        21 Feb. 1983    RAK    "End-of-curve" mark may begin in any column.
C        13 June 1983    RAK    Added "end-of-frame" input option.  FINIS
C                               output replaces EOF flag.
C        21 July 1983    RAK    Updated call to IOCHEK (added CONVER).
C        23 Jan. 1984    RAK    Extensively revised, but old data sets are
C                               upward compatible except that "END CURVE"
C                               and "END FRAME" marks are less free.
C                               Multiple-column input format and the packed
C                               data structure for the X and Y arrays are
C                               the significant changes.  A bug in the
C                               handling of the MAXCURth curve was repaired.
C        18 June 1984    RAK    Reduced BACKSPACE-ing by using TOKENS and BN
C                               format (much faster now on large files).
C                               List-directed READs eliminated.  Maximum
C                               number of columns and data field width
C                               changed.  Some GO TOs used to reduce level
C                               of nesting.  Dropped warning on blank lines.
C                               Reordered arguments.
C         5 Feb. 1985    RAK    Moved initialization of POINT to beginning
C                               of routine (exit processing was sometimes
C                               wrong when array bounds were violated).
C        18 May  1987  RAK/DAS  GETLINE introduced for more efficient reads
C                               (VAX-dependent, but easily modified if
C                               required).  Revise FORMATs to use /' ' for
C                               pushing a blank line (not 0).  Print
C                               warnings only if LWRITE > 0 and fatal
C                               error messages only if LWRITE <> 0.
C        31 July 1987    RAK    Increased MAXCOL to 30.  Pass EXCLAM to
C                               GETLINE as comment character.
C         3 Aug. 1987    RAK    Increased MAXTOK from 24 to 40 (enough
C                               for Cray DOUBLE PRECISION, plus a little).
C         2 Jan. 1990  MDW/RAK  Interpretation of unexpected alpha or numeric 
C                               input passed to to a higher level.  Ensured the
C                               reading of a single column of Y-data even when
C                               YCOL set incorrectly. 
C         5 Mar. '90  M.D.Wong  Keyword dictionary now passed as argument to
C                               terminate curve read.
C        16 Nov. '91 D.Saunders Major revision with indirection and a keyword
C                               input scheme in mind:
C                               1. Avoid all backspacing (except at a higher
C                                  level for a namelist) by passing a meaningful
C                                  BUFFER (1) (1:LAST) in (and out).
C                               2. Reduce end-of-curve handling at this level
C                                  to simply identifying a nonnumeric 1st token.
C                                  The dictionary argument is thus redundant.
C                                  (Using LOOKUP on the first token of every
C                                  line was inefficient anyway: ALPHA can detect
C                                  a nonnumeric string more cheaply.)
C                               3. Handling of arbitrary text in column 1 with
C                                  XCOL, YCOL > 1 appears impossible for a
C                                  keyword scheme - no way of telling if the
C                                  text is a misspelled keyword.  However,
C                                  the special case of *** in column 1 is
C                                  allowed here.
C
C     Author:  Robert Kennelly, NASA Ames/Sterling Software, Palo Alto, CA.
C
C-----------------------------------------------------------------------


C     Declarations.
C     -------------

      IMPLICIT NONE

C     Arguments.

      LOGICAL
     >   FINIS
      INTEGER
     >   LAST, LREAD, LWRITE, MAXPT, MAXCUR, POINTS (MAXCUR), XCOL, YCOL
      REAL
     >   X (MAXPT), Y (MAXPT)
      CHARACTER
     >   BUFFERS (2) * (*)

C     Local constants.

      INTEGER
     >   MAXCOL, MAXTOK
      CHARACTER
     >   EXCLAM * 1
      PARAMETER
     >  (EXCLAM = '!',
     >   MAXCOL = 30,   ! Highest accessible column number.
     >   MAXTOK = 40)   ! Largest data field width which may be read.

C     Local variables.

      INTEGER
     >   BASE, CURVE, ERROR, INSERT, IOS, NUMCOL, POINT, XYMAX
      LOGICAL
     >   ALPHA, CONVER, DONE, EOF, FATAL, LDSYN, NUMERIC, OK
      CHARACTER
     >   CONCAT * (2 * MAXTOK), LIST (0:MAXCOL) * (MAXTOK), XYFMT * 11

C     Procedures.

      EXTERNAL
     >   ALPHA, GETLINE, IOCHEK, TOKENS

C     Storage.

      DATA
     >   XYFMT /'(BN,2Fnn.0)'/
      SAVE
     >   BASE, XYFMT


C     Execution.
C     ----------

      ERROR = 0
      POINT = 0
      CURVE = CURVE + 1
      FINIS = .FALSE.
      FATAL = .FALSE.

C     Variable format in portable form:

      WRITE (XYFMT (7:8), '(I2)') MAXTOK

C     (Re)set local pointer to beginning of storage for this curve.

      IF (CURVE .EQ. 1) THEN
         BASE = 1
      ELSE
         BASE = BASE + POINTS (CURVE - 1)
      END IF

C     We can't check for too many curves (or pts.) till we're sure we have data.

C     Check specified column numbers.

      XYMAX = MAX (XCOL, YCOL)
      IF (XYMAX .GT. MAXCOL) THEN

C        Too many columns requested for present array dimensions.

         IF (LWRITE .GT. 0) THEN
            WRITE (LWRITE, 1030)
            WRITE (LWRITE, 1040) XCOL, YCOL
         END IF

         ERROR = -1
         GO TO 800
      END IF

      IF (XCOL .LE. 0 .AND. YCOL .LE. 0) THEN

C        Nonsense input values.

         IF (LWRITE .GT. 0) THEN
            WRITE (LWRITE, 1050) CURVE
            WRITE (LWRITE, 1040) XCOL, YCOL
         END IF

         ERROR = -1
         GO TO 800
      END IF


C     Look for a curve.
C     -----------------

      DONE = .FALSE.
  200 CONTINUE
         POINT = POINT + 1

         IF (POINT .GT. 1 .OR. LAST .EQ. 0) THEN   ! Read another line.

            CALL GETLINE (LREAD, EXCLAM, BUFFERS (1), LAST, IOS)
            CALL IOCHEK (IOS, OK, EOF, LDSYN, CONVER)

         ELSE             ! The calling program passes a valid buffer initially.
            OK = .TRUE.   ! This is simpler than avoiding the next test.
         END IF

         IF (.NOT. OK) THEN
            POINT = POINT - 1
            DONE = .TRUE.

            IF (EOF) THEN     ! We're done (EOF can terminate last curve).
               FINIS = .TRUE.
            ELSE              ! Fatal read error.
               FATAL = .TRUE.
            END IF

         ELSE

C           Analyze the input line.
C           -----------------------

            IF (LAST .EQ. 0) THEN   ! Skip a blank line.  
               POINT = POINT - 1
               GO TO 200
            END IF

            NUMCOL = XYMAX
            CALL TOKENS (BUFFERS (1) (1:LAST), NUMCOL, LIST (1))

            IF (ALPHA (LIST (1) (1:3))) THEN   ! First 3 chars. are not numeric.

C              The token is MOST LIKELY a control keyword.  Return to caller.

C              SPECIAL CASE:  Allow a bunch of ***s in column 1 which will
C                             be ignored if XCOL, YCOL are not 1.

               NUMERIC = LIST (1) (1:3) .EQ. '***'

               IF (.NOT. NUMERIC) THEN   ! Signal some kind of keyword
                  ERROR = +1
                  POINT = POINT - 1               
                  DONE  = .TRUE.
               END IF
C
C              DISTURBING REALIZATION:
C
C              The possibility of having arbitrary text in column 1 along with
C              other columns of numeric data specified by XCOL and YCOL cannot
C              be handled here.  Even at the higher level, there is no way of
C              distinguishing this case from the case of a misspelled keyword
C              followed by (possibly numeric) keyword value(s).
              

            ELSE    ! The first token cannot be a keyword.

               NUMERIC = .TRUE.

            END IF

            IF (NUMERIC) THEN

C              Are there enough columns for the numeric data?

               IF (NUMCOL .LT. XYMAX) THEN
              
                  IF (XCOL .EQ. 0 .AND. NUMCOL .EQ. 1) THEN

C                    Special case for a single column of Y-data where
C                    QPLOT's default YCOL = 2 is corrected here.  [Note
C                    that a single-column row among 2-column rows in
C                    this case would give spurious results.  DAS, 11/16/91]

                     YCOL = 1
                  ELSE
                     IF (LWRITE .GT. 0) THEN
                        WRITE (LWRITE, 1060) 
     >                     CURVE, NUMCOL, BUFFERS (1) (1:60)
                        WRITE (LWRITE, 1020) POINT, CURVE
                     END IF
                     POINT = POINT - 1
                     ERROR = -1
                     DONE = .TRUE.
                  END IF
               END IF
            END IF


            IF (.NOT. DONE) THEN

               IF (POINT .EQ. 1) THEN         ! Don't test curve # every pt. #.
                  IF (CURVE .GT. MAXCUR) THEN

C                     We've reached the curve limit but not the end of the data.

                      IF (LWRITE .GT. 0) WRITE (LWRITE, 1010) CURVE
                      POINT = POINT - 1
                      ERROR = -1
                      GO TO 800
                  END IF
               END IF

C              Read the data internally in floating point format.
C              --------------------------------------------------

               INSERT = BASE + POINT - 1
               IF (INSERT .LE. MAXPT) THEN

C                 One or the other of XCOL, YCOL may be zero.

                  IF (XCOL * YCOL .EQ. 0)
     >               WRITE (LIST (0), XYFMT) FLOAT (POINT)

C                 "Read" both data points at once to save error handling.

                  CONCAT = LIST (XCOL) // LIST (YCOL)
                  READ (CONCAT, XYFMT, IOSTAT = IOS)
     >               X (INSERT), Y (INSERT)
                  CALL IOCHEK (IOS, OK, EOF, LDSYN, CONVER)
                  IF (.NOT. OK) THEN
                     IF (CONVER) THEN
                        ERROR = +2
                        WRITE (BUFFERS (2), 1025) POINT, CURVE
                     ELSE
                        FATAL = .TRUE.
                     END IF
                     POINT = POINT - 1
                     DONE = .TRUE.
                  END IF
               ELSE

C                 We have run out of array space.

                  IF (LWRITE .GT. 0) THEN
                     WRITE (LWRITE, 1070)
                     WRITE (LWRITE, 1020) POINT, CURVE
                  END IF
                  POINT = POINT - 1
                  ERROR = -1
                  DONE = .TRUE.
               END IF
            END IF
         END IF

         IF (.NOT. DONE) GO TO 200     ! Look for more data


C     Termination.
C     ------------

  800 CONTINUE
      IF (FATAL) THEN

C        Abnormal termination due to fatal read error - bail out!

         IF (LWRITE .NE. 0) THEN
            WRITE (ABS (LWRITE), 1080) IOS
            WRITE (ABS (LWRITE), 1020) POINT + 1, CURVE
         END IF
         CALL EXIT (3)
      ELSE

C        Check whether last attempt actually read any data and adjust
C        counters accordingly.

         IF (POINT .GT. 0) THEN
            POINTS (CURVE) = POINT
         ELSE
            CURVE = CURVE - 1
            IF (CURVE .GT. 0) BASE = BASE - POINTS (CURVE)
         END IF
      END IF

      RETURN


C     Formats.
C     --------

 1010 FORMAT (/' READCOLS: Trouble - too many curves for array size.'/
     >   11X, '(curve number ', I5, ')')
 1020 FORMAT (11X, '(point number ', I5, ', curve number ', I5, ')')
 1025 FORMAT ('(point number ', I5, ', curve number ', I5, ')')
 1030 FORMAT (/' READCOLS: Trouble - too many data columns requested.')
 1040 FORMAT (11X, '(XCOL = ', I3, ', YCOL = ', I3, ')')
 1050 FORMAT (/' READCOLS: Trouble - illegal input column data for ',
     >   'curve number ', I5)
 1060 FORMAT (/' READCOLS: Trouble - not enough columns of input data ',
     >   'for curve'/
     >   11X, 'number ', I5, '.  Only ', I3, ' column(s) were found !'/
     >   11X, 'The questionable data began with:'/
     >   11X, '>>', A, '<<')
 1070 FORMAT (/' READCOLS: Trouble - too many data points.  Terminate ',
     >   'this frame.')
 1080 FORMAT (/' READCOLS: Trouble - fatal READ error, IO status  = ',
     >   I6, '.')

      END
C+----------------------------------------------------------------------
C
      SUBROUTINE READER
C
C
C     Description and usage:
C
C           READER provides a set of utilities offering some of the
C        flexibility of list-directed input to interactive FORTRAN programs,
C        while permitting use of <CR> (i.e. Carriage Return) as a legal,
C        identifiable response.  This is impossible with a list-directed
C        read, which ignores <CR> and continues to wait for input.  In
C        case of bad input (presumably due to user error), a warning is
C        displayed and the prompt repeated.  "Quit" (as opposed to "default")
C        is provided for by the "End-of-file" argument EOF.  (See Notes.)
C        If CR or EOF is .TRUE., the argument value returned is unchanged,
C        so these flags should be checked by the calling routine.  (However,
C        assigning a default value before calling READER can save checking
C        the CR argument.)
C
C           Several ENTRY points are provided for reading INTEGER, REAL,
C        DOUBLE PRECISION, STRING, CHARACTER, or YES/NO data.  The names
C        are READx, where x = I, R, D, S, C, or Y.  READC returns all the
C        non-blank characters found, converted to upper case if alphabetic
C        and packed from the left into the output string with blanks removed.
C        The READS option merely returns the string as entered from the
C        keyboard, without modification.  READC is normally appropriate for
C        entering a single item or token (such as an identifier or a file
C        name), while READS is appropriate for literal strings such as plot
C        titles.  HOWEVER: since Unix file names are case sensitive, READS
C        is actually the better choice for file name prompts on Unix systems.
C
C           All calling sequences are identical except for the type of the
C        value to be returned.  Note that the "top" of the subroutine (READER)
C        is NOT a legitimate entry point!
C
C           Sample usage:
C
C           :      :                              (^D under Unix)
C       210 NPTS = 100                              |
C           CALL READI (LUNCRT,                     |
C          >   'Enter number of points.  <CR>=100; ^Z=quit: ',
C          >   LUNKBD, NPTS, CR, EOF)
C           IF (EOF) GO TO 999
C           IF (NPTS .LT. 1 .OR. NPTS .GT. MXPTS) GO TO 210
C           :      :
C
C
C     Arguments:
C
C        Name    Dimension  Type  I/O/S  Description
C        SCREEN              I    I      Logical unit number to which
C                                        the prompt is written (often 6).
C        PROMPT    *         C    I      Prompt string to be written to
C                                        unit SCREEN.
C        KEYBRD              I    I      Logical unit number from which
C                                        data is to be read (often 5).
C
C
C   ---> Only ONE of the following six output choices is to be used:
C
C        INT4                I      O    INTEGER quantity to be returned.
C        REAL4               R      O    REAL quantity to be returned.
C        REAL8               D      O    DOUBLE PRECISION quantity to be
C                                        returned.
C        CHARS               C      O    A CHARACTER string, converted to
C                                        upper case if alphabetic, packed
C                                        and left-justified.
C        STRING              C      O    A CHARACTER string, as entered.
C        YES                 L      O    Logical flag set .TRUE. for a 'YES'
C                                        response (first character 'Y' or 'y')
C                                        and .FALSE. for 'N' or 'n'.  Other
C                                        inputs force a reprompt.
C        CR                  L      O    CR = .TRUE. if a null value (i.e.
C                                        Carriage Return only) was read,
C                                        else CR = .FALSE.
C        EOF                 L      O    EOF = .TRUE. if "End-of-File"
C                                        was detected (^Z under VMS; ^D
C                                        under Unix), else EOF = .FALSE.
C                                        
C
C     External devices:  See arguments SCREEN and KEYBRD.
C
C
C     External procedures:  UPCASE is used by the READC and READY options
C                           in place of the original in-line code - see Notes.
C
C
C     Environment:  Digital VAX-11/780 VMS FORTRAN 77
C                   Also:   VAX/ULTRIX and Silicon Graphics IRIS 4D
C
C
C     Method:
C
C           The user is prompted for input, with the cursor left at the
C        end of the prompt if there is room for the expected response.  The
C        reply is read into a buffer as a character string and, except for
C        the READC and READS entries, left-justified with blanks squeezed
C        out and re-read internally in the requested format.  For READC, the
C        same repacking is done, then the result is converted to upper case.
C        No repacking or conversion is done in the case of READS.
C
C
C     Notes:
C
C        (1)  IMPLICIT NONE is non-standard.
C
C        (2)  Conversions to upper case were originally done in-line by
C             adding the appropriate offset (= ICHAR ('A') - ICHAR ('a')).
C             However this has been shown to fail on some systems and has
C             been replaced by use of the UPCASE utility, which cannot fail.
C
C        (3)  The length of the string returned by READC and READS is set by
C             the calling program - excess characters read will be truncated
C             from the right.  No error flag is set for this condition.  The
C             maximum length for CHARS and STRING is MAXBUF.
C
C        (4)  The '$' carriage control character used for short prompts is
C             a VAX extension.  It starts output at the beginning of the next
C             line, and suppresses the carriage return at the end of the line.
C             This feature may have to be modified or eliminated on other
C             systems.  Indeed, the ('$', A) form was found not to work on
C             the IRIS 4D, and has been replaced by the alternative DEC form
C             (A, $) which is compatible with the IRIS 4D.
C
C
C     Author:  Robert Kennelly, NASA Ames/Sterling Software, Palo Alto, CA.
C
C
C     Development history:
C
C        27 Sep. 1983    RAK    Initial design and coding.
C         6 Oct. 1983    RAK    Extended to multiple data types.
C         2 Nov. 1983  RAK/LJC  Added READS for strings without modification.
C        16 Mar. 1984  RAK/LJC  Included "BN" specifier in formats to avoid
C                               unneccessary packing.
C        14 June 1984  DAS/RAK  Corrected header, TEST always defined, added
C                               CC, eliminated redundant trap on CASE, used
C                               shorter flags, e.g. 'R' instead of 'REAL'.
C         7 Sep. 1984    RAK    Leave STRING unchanged when CR entered,
C                               thus consistent with the other modes.
C        19 Oct. 1984    RAK    Initialize BUFFER prior to READ.
C        14 Dec. 1984    RAK    Make sure length of CHARS is not exceeded
C                               during repacking/case-conversion.
C        19 Dec. 1985  RAK/DBS  Addition of READY entry, and general 
C                               streamlining of code.
C        30 Dec. 1985  RAK/DAS  Edited header.
C        09 Sep. 1988    DAS    Sample usage added above; other cosmetics.
C        05 May  1989  DAS/RGL  Unix- and VMS-compatible now: ^Z/^D usage
C                               documented; STOP 'READER' message <=8 chars.;
C                               OFFSET not defined as a PARAMETER constant.
C                               (UNICOS displays '$' in column 1 at time of
C                               writing, but this may go away...)
C        01 Feb. 1990    DAS    '$' carriage control changed to (A,$) form
C                               (see Notes above; thanks to Scott Thomas).
C                               READS recommended over READC for file names
C                               on Unix systems.  Conversions to upper case
C                               now done via UPCASE - see Notes above.
C
C-----------------------------------------------------------------------


C     DECLARATIONS.
C     -------------

      IMPLICIT NONE

C     Local constants.  Note that the field lengths in the format strings
C     are hardcoded to MAXBUF.

      INTEGER
     >   MAXBUF, MAXLIN
      CHARACTER
     >   BLANK, FORMI * 9, FORMR * 11
      PARAMETER
     >  (BLANK  = ' ',
     >   FORMR  = '(BN, F80.0)',
     >   FORMI  = '(BN, I80)',
     >   MAXBUF = 80,
     >   MAXLIN = 80)

C     Arguments.

      INTEGER
     >   INT4, KEYBRD, SCREEN
      REAL
     >   REAL4
      DOUBLE PRECISION
     >   REAL8
      CHARACTER
     >   CHARS * (*), PROMPT * (*), STRING * (*)
      LOGICAL
     >   CR, EOF

C     Local variables.

      LOGICAL
     >   YES
      INTEGER
     >   FILL, SCAN, STATUS, TEST
      CHARACTER
     >   BUFFER * (MAXBUF), CASE, FORMP * 8

C     External references.

      EXTERNAL
     >   UPCASE

C     Storage.

      DATA
     >   FORMP /'(1X,A  )'/        ! Valid blanks embedded so that the (A,$)
                                   ! form may be commented out below if reqd.


C     EXECUTION.
C     ----------

      STOP 'READER'                ! Illegal entry point.

      ENTRY READI (SCREEN, PROMPT, KEYBRD, INT4, CR, EOF)
         CASE = 'I'
         GO TO 10

      ENTRY READR (SCREEN, PROMPT, KEYBRD, REAL4, CR, EOF)
         CASE = 'R'
         GO TO 10

      ENTRY READD (SCREEN, PROMPT, KEYBRD, REAL8, CR, EOF)
         CASE = 'D'
         GO TO 10

      ENTRY READS (SCREEN, PROMPT, KEYBRD, STRING, CR, EOF)
         CASE = 'S'
         TEST = LEN (STRING)
         GO TO 20

      ENTRY READC (SCREEN, PROMPT, KEYBRD, CHARS, CR, EOF)
         CASE = 'C'
         TEST = LEN (CHARS)
         GO TO 20

      ENTRY READY (SCREEN, PROMPT, KEYBRD, YES, CR, EOF)
         CASE = 'Y'
         TEST = 5
         GO TO 20

   10 CONTINUE

C     Set a typical maximum length for response to the non-character
C     entries.  This value and the length of the prompt string are used
C     below to determine where the cursor is to be positioned after
C     prompting.

      TEST = 16

   20 CONTINUE

C     Set format for displaying prompt.  The IF determines whether to leave
C     the cursor at the end of the prompt line (VAX, IRIS 4D - may be deleted).

      IF (TEST + LEN (PROMPT) .LE. MAXLIN) THEN
         FORMP (6 : 7) = ',$'
      ELSE
         FORMP (6 : 7) = BLANK
      END IF

C     Prompt user for input.
C     ----------------------

   30 CONTINUE
         BUFFER = BLANK

C        Issue prompt, then read up to MAXBUF characters from the keyboard.

         WRITE (SCREEN, FORMP) PROMPT
         READ (KEYBRD, '(A)', IOSTAT = STATUS) BUFFER

C        Re-try on errors, exit on End-of-File, or continue.

         IF (STATUS .GT. 0) THEN
            WRITE (SCREEN, 1000)
            GO TO 30

         ELSE IF (STATUS .LT. 0) THEN
            EOF = .TRUE.
            CR = .FALSE.
            GO TO 990

         ELSE
            EOF = .FALSE.
            CR = (BUFFER .EQ. BLANK)
         END IF


C        Process the data in the buffer, if any.
C        ---------------------------------------

         IF (.NOT. CR) THEN

            IF (CASE .EQ. 'I') THEN
               READ (BUFFER, FORMI, IOSTAT = STATUS) INT4

            ELSE IF (CASE .EQ. 'R') THEN
               READ (BUFFER, FORMR, IOSTAT = STATUS) REAL4

            ELSE IF (CASE .EQ. 'D') THEN
               READ (BUFFER, FORMR, IOSTAT = STATUS) REAL8

            ELSE IF (CASE .EQ. 'S') THEN

C              For literal strings, copy the buffered input directly to
C              the output variable.  If CR is true, we leave STRING alone
C              so that any default value set in the calling routine will
C              remain intact.  No error checking is required.

               STRING = BUFFER

            ELSE 

C              CASE is either 'C' or 'Y' - scan the input buffer, starting
C              from the left, to pack the data and count the non-blanks.

               FILL = 0
               DO 40 SCAN = 1, MAXBUF
                  IF (BUFFER (SCAN : SCAN) .NE. BLANK) THEN
                     FILL = FILL + 1
                     BUFFER (FILL : FILL) = BUFFER (SCAN : SCAN)
                  END IF
   40          CONTINUE

C              Convert to upper case.

               CALL UPCASE (BUFFER (1 : FILL))

               IF (CASE .EQ. 'C') THEN

                  CHARS = BUFFER (1 : FILL)

               ELSE

C                 CASE must be 'Y' - all we need is the first character.

                  IF (BUFFER (1 : 1) .EQ. 'Y') THEN
                     YES = .TRUE.
                  ELSE IF (BUFFER (1 : 1) .EQ. 'N') THEN
                     YES = .FALSE.
                  ELSE

C                    Buffer value is invalid - generate an ersatz input error
C                    to force a re-prompt.

                     STATUS = 999
                  END IF
               END IF
            END IF 

C           Re-try on errors, exit on End-of-File, else continue.
C           (Dropping through for 'C', 'S', 'Y' cases is OK and saves code.)

            EOF = (STATUS .LT. 0)
            IF (STATUS .GT. 0) THEN
               WRITE (SCREEN, 1000)
               GO TO 30

            END IF

      END IF


C     TERMINATION.
C     ------------

  990 RETURN


C     FORMATS.
C     --------

 1000 FORMAT (' Input error!  Please try again.')

      END
C+----------------------------------------------------------------------
C
      SUBROUTINE RESCALE (MAXPT, NPT, CURVE, X, SCALE, SHIFT)
C
C
C     One-liner:  Linear transformation of data packed in an array
C
C
C     Description and usage:
C
C           Another hide-the-ugly-details utility which applies a linear
C        transformation to data stored in a packed array.  Since the data
C        array may be large, we check whether anything really needs to be
C        done before proceding.  The order of operations is as follows:
C
C                    X (output) = SCALE * X (input) + SHIFT
C
C        Hard stops are provided if MAXPT or the NPT array are unreasonable.
C
C           RESCALE was written for use by QPLOT, a general-purpose data
C        graphics program developed at NASA-Ames.  It could find other
C        applications, but may not pay unless the data is actually packed.
C
C
C     Arguments:
C
C        Name    Dimension  Type  I/O/S  Description
C        MAXPT              I     I      Declared size of data array in
C                                        calling routine (required for
C                                        error checking).
C
C        NPT     CURVE      I     I      Point counts for each curve. Simple
C                                        tests are performed to verify that
C                                        the computed data limits within
C                                        the X array are sensible, i.e.
C                                        > 0 and <= MAXPT.
C
C        CURVE              I     I      Index of curve to be checked.
C
C        X       MAXPT      R     I/O    Array to be transformed.
C
C        SCALE              R     I      Scale factor.  If SCALE = 1.0 and
C                                        SHIFT = 0.0, just return.
C
C        SHIFT              R     I      Offset.
C
C
C     Environment:  Digital VAX-11/780 VMS FORTRAN 77
C
C
C     Notes:
C
C        (1)  IMPLICIT NONE and 8-character symbolic names are non-standard.
C
C
C     Author:  Robert Kennelly, NASA-Ames/Sterling Software
C
C
C     Development history:
C
C        21 Apr. 1987    RAK    Initial design and coding.
C
C-----------------------------------------------------------------------


C     Declarations.
C     -------------

      IMPLICIT NONE

C     Arguments.

      INTEGER
     >   CURVE, MAXPT, NPT (CURVE)
      REAL
     >   SCALE, SHIFT, X (MAXPT)

C     Local variables.

      INTEGER
     >   FIRST, I, LAST


C     Execution.
C     ----------

      IF (MAXPT .LE. 0) STOP 'RESCALE:  Array dimension error!'
      IF (SCALE .EQ. 1.0 .AND. SHIFT .EQ. 0.0) GO TO 990

C     Count the points in each curve up to this one.

      LAST = 0
      DO 10, I = 1, CURVE
         LAST = LAST + NPT (I)
   10 CONTINUE
      IF (LAST .LE. 0) STOP 'RESCALE:  No data!'
      IF (LAST .GT. MAXPT) STOP 'RESCALE:  (Apparently) too much data!'

C     Set index for the begining of this curve's data points.

      FIRST = LAST - NPT (CURVE) + 1

C     Transform the data (at least one of SCALE or SHIFT was found to be
C     active above).

      DO 20, I = FIRST, LAST
         X (I) = SCALE * X (I) + SHIFT
   20 CONTINUE


C     Termination.
C     ------------

  990 CONTINUE
      RETURN
      END
C+----------------------------------------------------------------------
C
      SUBROUTINE ROTATE (MAXPT, NPT, CURVE, X, Y, RADIANS, CENTERX,
     >                   CENTERY)
C
C
C     One-liner:  Rotation of data packed in an array about a point
C
C
C     Description and usage:
C
C           Another hide-the-ugly-details utility which rotates one curve
C        of data stored in a packed array about an arbitrary point, in place.
C
C           This analogue of RESCALE (q.v.) was forced for use in QPLOT
C        because of the scheme used for identifying the start and end of
C        a given curve.
C
C           Error checking is omitted because the application should have
C        done it already.
C
C
C     Arguments:
C
C        Name    Dimension  Type  I/O/S  Description
C        MAXPT              I     I      Declared size of data array in
C                                        calling routine (required for
C                                        dimensioning purposes only).
C
C        NPT     CURVE      I     I      Point counts for each curve.
C
C        CURVE              I     I      Index of curve to be rotated.
C
C        X       MAXPT      R     I/O    Arrays to be transformed.
C        Y       MAXPT
C
C        RADIANS            R     I      Angle of rotation.  Positive is
C                                        counterclockwise.  Radians are
C                                        expected because QPLOT allows
C                                        either degrees or radians, and
C                                        any conversion here could mean
C                                        two approximations where none was
C                                        needed.
C
C        CENTERX            R     I      Center of rotation.
C        CENTERY
C
C
C     Environment:  Digital VAX/VMS, FORTRAN 77 + ...
C                   IMPLICIT NONE
C                   Names up to 8 characters
C
C     History:
C
C        4 Apr. 1993    DAS    Adaptation of ROTATE2D necessitated by
C                              QPLOT's packed-curve scheme.
C
C     Author:  David Saunders, Sterling Software, NASA-Ames, Mt. View, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments.

      INTEGER
     >   CURVE, MAXPT, NPT (CURVE)
      REAL
     >   CENTERX, CENTERY, RADIANS, X (MAXPT), Y (MAXPT)

C     Local variables.

      INTEGER
     >   FIRST, I, LAST
      REAL
     >   CI, SI, DX, DY

C     Execution.

      IF (RADIANS .EQ. 0.0) GO TO 99

C     Count the points in each curve up to this one.

      LAST = 0
      DO 10, I = 1, CURVE
         LAST = LAST + NPT (I)
   10 CONTINUE

C     Set index for the begining of this curve's data points.

      FIRST = LAST - NPT (CURVE) + 1

C     Transform the coordinates:

      CI = COS (RADIANS)
      SI = SIN (RADIANS)

      DO 20, I = FIRST, LAST
         DX = X (I) - CENTERX
         DY = Y (I) - CENTERY
         X (I) = DX * CI - DY * SI + CENTERX
         Y (I) = DX * SI + DY * CI + CENTERY
   20 CONTINUE

   99 RETURN
      END
C+----------------------------------------------------------------------
C
      SUBROUTINE SCANNR (STRING, FIRST, LAST, MARK)
C
C
C     Description and usage:
C
C           Looks for significant fields ("tokens") in a string, where the
C        fields are of arbitrary length and separated by blanks, tabs, commas,
C        colons, or equal signs.  The position of the end of the first token
C        is also returned so that this routine may be conveniently used within
C        a loop to process an entire line of text.
C
C           The procedure examines a substring, STRING (FIRST : LAST), which
C        may of course be the entire string (in which case just call SCANNR
C        with FIRST = 1 and LAST = LEN (STRING) ).  The indices returned
C        are relative to STRING itself, not the substring.
C
C
C     Arguments:
C
C        Name    Dimension  Type  I/O/S  Description
C        STRING    *         C    I      Text string containing data to be
C                                        scanned.
C        FIRST               I    I/O    Index of first character of interest
C                                        in STRING.  FIRST >= 1 is assumed.
C                                        Output is index of the beginning of
C                                        the first token, or 0 if no token
C                                        was found.
C        LAST                I    I/O    Index of last character of interest
C                                        in STRING.  LAST <= LEN (STRING) is
C                                        assumed.  Output is index of the end
C                                        of the last non-separator, or 0 if no
C                                        token was found.
C        MARK                I      O    Points to the end of the first token
C                                        found in the specified portion of
C                                        STRING.  Use FIRST = MARK + 2 for
C                                        further searches in the same string.
C                                        MARK is set to 0 if no token was found.
C
C
C     Environment:  Digital VAX-11/780, VMS, FORTRAN 77
C                   Also: SGI IRIS 4D, IRIX 3.3, f77
C                         Cray Y-MP, UNICOS, cft77
C     Notes:
C
C        (1)  IMPLICIT NONE is non-standard.  Constant HT (Tab) is defined
C             in a standard way since the CHAR function is not permitted
C             in a PARAMETER declaration (OK on VAX, though).  For Absoft
C             FORTRAN 77 on 68000 machines, use HT = 9.
C
C        (2)  The pseudo-recursive structure of the original version has been
C             abandoned because the VAX compiler treated the SOLID statement
C             function as an in-line subroutine, with substantial penalty
C             in speed (factor of almost 3!).  The single-loop form used
C             later was almost as fast (especially for lines with only a few
C             tokens), and was more compact and easy to change if a different
C             set of delimiters was required.  However, ...
C
C        (3)  This version adjusts LAST using a BACKward search, which is
C             degenerate for all tokens after the first.  Repeated forward
C             scanning to locate LAST (in a calling program's loop over tokens)
C             amounts to an O(4N**2) operation count for N tokens in a string
C             versus the O(2N) that it should be.  (The 2 is in there because
C             the non-token fields are as numerous as the tokens and take a
C             similar time to scan.)  The price paid is some repeated code.
C
C        (4)  FIRST >= 1 and LAST <= LEN (STRING) are assumed upon entry,
C             because these are the obvious bounds needed by the calling
C             program on the first scan of a given string - no need to
C             evaluate LEN (STRING) here with every call, while use of
C             MIN and MAX tends to disguise programming errors.
C
C        (5)  Entering with FIRST > LAST is an eventual normal consequence
C             of setting FIRST = MARK + 2 for further scans of the same string.
C             The zero DO loop iteration count feature of the language is
C             taken advantage of here - the result is simply to drop
C             out the bottom with FIRST = LAST = MARK = 0.
C
C        (6)  Zeroing all of FIRST, LAST, and MARK (rather than just MARK)
C             when no token is found may, in retrospect, be undesirable in
C             that information is lost (esp. LAST), but there are too many
C             applications to change this now.
C
C        (7)  The variety of separators recognized limits the usefulness of
C             this routine somewhat.  The intent is to facilitate handling
C             such tokens as keywords or numerical values.  In other
C             applications, it might be necessary for ALL printing characters
C             to be significant.  Use SCAN2 (from the same author) if user-
C             supplied delimiters are appropriate.
C
C        (8)  Note that "null" characters are not treated here.  This should
C             not be a problem in that FORTRAN READs of short records, like
C             assignment of short strings, pads with blanks.
C
C
C     Author:  Robert Kennelly, Informatics General Corporation.
C
C
C     History:
C
C        29 Dec. 1984    RAK    Initial design and coding, (very) loosely
C                               based on SCAN_STRING by Ralph Carmichael.
C        25 Feb. 1984    RAK    Added ':' and '=' to list of separators.
C        16 Apr. 1985    RAK    Defined SOLID in terms of variable DUMMY
C                               (previous re-use of STRING was ambiguous).
C         3 Mar. 1986    RAK    Restructured, without SOLID.  Checks for
C                               "state transitions" while executing a
C                               single forward DO loop.  Protect against
C                               funny input (FIRST > LAST).
C        13 May  1988    DAS    Introduced backward search for LAST;
C                               eliminated MIN, MAX and LEN.  (See Notes.)
C        10 June 1991    DAS    Cray's cft77 doesn't allow CHAR (9) in a
C                               PARAMETER statement.
C
C-----------------------------------------------------------------------


C     Declarations.
C     -------------

      IMPLICIT NONE

C     Arguments.

      INTEGER
     >   FIRST, LAST, MARK
      CHARACTER
     >   STRING * (*)

C     Local constants.

      CHARACTER
     >   BLANK, EQUAL, COLON, COMMA
      PARAMETER
     >   (BLANK = ' ',
     >    EQUAL = '=',
     >    COLON = ':',
     >    COMMA = ',')

C     Local variables.

      INTEGER
     >   HEAD, I, J, TAIL
      LOGICAL
     >   FOUND
      CHARACTER
     >   HT * 1

C     Execution.
C     ----------

      HT = CHAR (9)
      HEAD = FIRST
      TAIL = LAST

      FIRST = 0
      LAST = 0
      MARK = 0
      FOUND = .FALSE.

C     Look at each character in STRING (HEAD : TAIL) from left to right
C     until a token is found.  (Drop through if LAST > FIRST.)

      DO 30, I = HEAD, TAIL

C        Is the character a separator?

         IF (STRING (I : I) .EQ. BLANK) GO TO 20
         IF (STRING (I : I) .EQ. HT   ) GO TO 20
         IF (STRING (I : I) .EQ. COMMA) GO TO 20
         IF (STRING (I : I) .EQ. EQUAL) GO TO 20
         IF (STRING (I : I) .EQ. COLON) GO TO 20

C           Not a separator.  Check for the beginning of the FIRST token.

            IF (.NOT. FOUND) THEN
               FIRST = I
               FOUND = .TRUE.
            END IF
            GO TO 30

   20    CONTINUE

C           We found a separator.

            IF (FOUND) THEN

C              We just passed the "trailing edge" of the first token.

               MARK = I - 1
               GO TO 40
            END IF
   30 CONTINUE

C     We reached the last character while still seeking the first token.
C     Either this is part of the first token, or MARK and LAST are still zero.

      IF (FOUND) THEN
         MARK = TAIL
         LAST = TAIL
      END IF
      GO TO 99


   40 CONTINUE

C     We found the first token but haven't reached the end yet to adjust LAST.

      DO 60, I = TAIL, MARK + 1, -1

C        Keep searching backwards as long as we have a separator.

         IF (STRING (I : I) .EQ. BLANK) GO TO 60
         IF (STRING (I : I) .EQ. HT   ) GO TO 60
         IF (STRING (I : I) .EQ. COMMA) GO TO 60
         IF (STRING (I : I) .EQ. EQUAL) GO TO 60
         IF (STRING (I : I) .EQ. COLON) GO TO 60

C           Not a separator, so we've found LAST.

            LAST = I
            GO TO 99
   60 CONTINUE

C     Dropped through, so there were no more tokens.

      LAST = MARK
     
C     Termination.
C     ------------

   99 RETURN
      END
C+----------------------------------------------------------------------
C
      SUBROUTINE SCAN2 (STRING, SEPS, FIRST, LAST, MARK)
C
C
C     Description and usage:
C
C           Looks for non-blank fields ("tokens") in a string, where the
C        fields are of arbitrary length and separated by any of a set of
C        user-specified separators (e.g., blanks or commas).  The position
C        of the end of the first token is also returned so that this
C        routine may be conveniently used within a loop to process an
C        entire line of text.
C
C           The procedure examines a substring, STRING (FIRST : LAST), which
C        may of course be the entire string (in which case just call SCAN2
C        with FIRST = 1 and LAST = LEN (STRING) ).  The indices returned
C        are relative to STRING itself, not the substring.
C
C
C     Arguments:
C
C        Name    Dimension  Type  I/O/S  Description
C        STRING    *         C    I      Text string containing data to be
C                                        scanned.
C        SEPS      *         C    I      String containing the separators.
C                                        Each character in SEPS counts as a
C                                        token delimiter.
C        FIRST               I    I/O    Index of first character of interest
C                                        in STRING.  FIRST >= 1 is assumed.
C                                        Output is index of the beginning of
C                                        the first non-separator, or 0 if no
C                                        token was found.
C        LAST                I    I/O    Index of last character of interest
C                                        in STRING.  LAST <= LEN (STRING) is
C                                        assumed.  Output is index of the end
C                                        of the last non-separator, or 0 if no
C                                        token was found.
C        MARK                I      O    Points to the end of the first token
C                                        found in the specified portion of
C                                        STRING.  Use FIRST = MARK + 2 for
C                                        further searches in the same string.
C                                        MARK is set to 0 if no token was found.
C
C
C     Environment:  Digital VAX-11/780 VMS FORTRAN (FORTRAN 77).
C
C
C     Notes:
C
C        (1)  IMPLICIT NONE is non-standard.
C
C        (2)  FIRST >= 1 and LAST <= LEN (STRING) are assumed upon entry,
C             because these are the obvious bounds needed by the calling
C             program on the first scan of a given string - no need to
C             evaluate LEN (STRING) here with every call, while use of
C             MIN and MAX tends to disguise programming errors.
C
C        (3)  Entering with FIRST > LAST is an eventual normal consequence
C             of setting FIRST = MARK + 2 for further scans of the same string.
C             The zero DO loop iteration count feature of the language is
C             taken advantage of here - the result is simply to drop
C             out the bottom with FIRST = LAST = MARK = 0.
C
C        (4)  This version adjusts LAST using a BACKward search, which is
C             degenerate for all tokens after the first.  (Doing it in
C             the one forward loop means unnecessary repeated tokenizing
C             to find the end.)
C
C        (5)  Zeroing all of FIRST, LAST, and MARK (rather than just MARK)
C             when no token is found may, in retrospect, be undesirable in
C             that information is lost (esp. LAST), but there are too many
C             applications to change this now.
C
C
C     Author:  Robert Kennelly, Informatics General Corporation.
C
C
C     Development history:
C
C         4 Mar 1986    RAK    Variation of SCANNR, which is hard-coded
C                              (for historical reasons and a modest speed
C                              advantage) for separators BLANK, TAB,
C                              COMMA, COLON, and EQUAL.
C
C         5 May  1988   DAS    Reverse search used to find LAST; MAX, MIN,
C                              and LEN also eliminated (see Notes).
C
C-----------------------------------------------------------------------


C     Declarations.
C     -------------

      IMPLICIT NONE

C     Arguments.

      INTEGER
     >   FIRST, LAST, MARK
      CHARACTER
     >   SEPS * (*), STRING * (*)

C     Local variables.

      INTEGER
     >   HEAD, I, J, NSEPS, TAIL
      LOGICAL
     >   FOUND


C     Execution.
C     ----------

      NSEPS = LEN (SEPS)
      HEAD = FIRST
      TAIL = LAST

      FIRST = 0
      LAST = 0
      MARK = 0
      FOUND = .FALSE.

C     Look at each character in STRING (HEAD : TAIL) from left to right
C     until a token is found.  (Drop through if LAST > FIRST.)

      DO 30, I = HEAD, TAIL

C        Is the character a separator?

         DO 10, J = 1, NSEPS
            IF (STRING (I : I) .EQ. SEPS (J : J)) GO TO 20
   10    CONTINUE

C           Not a separator.  Check for the beginning of the FIRST token.

            IF (.NOT. FOUND) THEN
               FIRST = I
               FOUND = .TRUE.
            END IF
            GO TO 30

   20    CONTINUE

C           We found a separator.

            IF (FOUND) THEN

C              We just passed the "trailing edge" of the first token.

               MARK = I - 1
               GO TO 40
            END IF
   30 CONTINUE

C     We reached the last character while still seeking the first token.
C     Either this is part of the first token, or MARK and LAST are still zero.

      IF (FOUND) THEN
         MARK = TAIL
         LAST = TAIL
      END IF
      GO TO 99


   40 CONTINUE

C     We found the first token but haven't reached the end yet to adjust LAST.

      DO 60, I = TAIL, MARK + 1, -1

C        Keep searching as long as we have a separator.

         DO 50, J = 1, NSEPS
            IF (STRING (I : I) .EQ. SEPS (J : J)) GO TO 60
   50    CONTINUE

C           Not a separator, so we've found LAST.

            LAST = I
            GO TO 99
   60 CONTINUE

C     Dropped through, so there were no more tokens.

      LAST = MARK
     
C     Termination.
C     ------------

   99 RETURN
      END
C+----------------------------------------------------------------------
C
      SUBROUTINE SCAN3 (STRING, SEPS, FIRST, LAST, MARK)
C
C     One-liner:
C
C        Backward-search variant of SCAN2 - find last token in a (sub)string
C
C     Description and usage:
C
C           Looks for the LAST non-blank field ("token") in a (sub)string by
C        searching from right to left.  The fields are of arbitrary length
C        and separated by any of a set of user-specified separators (e.g.,
C        blanks or commas).  This routine may be conveniently used within
C        a loop to process an entire line of text BACKwards.  (However, the
C        need for "TOKEN3" to go with SCAN3 appears limited - just use SCAN3
C        on the substring defined by FIRST and LAST = MARK-2 if you want further
C        tokens in reverse order.)
C
C           To clarify:  Given STRING and positions FIRST, LAST,
C
C        SCAN3 locates RIGHT-most token as STRING (MARK : LAST)  (LAST updated);
C        SCAN2    "    LEFT   "    "    "  STRING (FIRST : MARK) (FIRST  "   ").
C
C           In both cases, MARK = 0 if no token was found.
C
C     Arguments:
C
C        Name    Dimension  Type  I/O/S  Description
C        STRING    *         C    I      Text string containing data to be
C                                        scanned.
C        SEPS      *         C    I      String containing desired separators.
C                                        Each character in SEPS is treated
C                                        as a possible token delimiter.
C        FIRST               I    I      Index of first character of STRING
C                                        of interest.  Not updated here.
C        LAST                I    I/O    Input as last character of STRING
C                                        of interest.  Output is right-most
C                                        end of right-most token found in
C                                        specified (sub)string.  Unchanged if
C                                        no token was found (MARK = 0).
C        MARK                I      O    Points to left-most character of the
C                                        first token found searching BACKward
C                                        in the specified substring. MARK = 0
C                                        means that no token was found.
C
C
C     Environment:  VAX/VMS FORTRAN 77
C
C
C     Notes:
C
C        (1)  IMPLICIT NONE is non-standard.
C
C        (2)  Error checking:  None.  The calling program is expected
C             to set FIRST and LAST correctly on the first call for a
C             given STRING (probably at 1 and LEN (STRING)) and update
C             LAST as MARK-2 for further searches of the same string.
C
C
C     Author:  Robert Kennelly/David Saunders, Sterling Software, Palo Alto.
C
C
C     Development history:
C
C      4 Mar. 1986   RAK   SCAN2 derived from SCANNR.
C     21 Apr. 1987   DAS   SCAN3 derived from SCAN2 for the backward search
C                          case.  Simplified with one-or-two-trailing tokens
C                          per line in mind.
C
C-----------------------------------------------------------------------


      IMPLICIT NONE

C     Arguments.

      INTEGER
     >   FIRST, LAST, MARK
      CHARACTER
     >   SEPS * (*), STRING * (*)

C     Local variables.

      INTEGER
     >   I, J, TAIL
      LOGICAL
     >   FOUND


C     Execution.
C     ----------

      MARK = 0
      TAIL = LAST
      FOUND = .FALSE.

C     Look at each character in STRING (FIRST : TAIL) from right to left
C     until a token is found or until FIRST is encountered.  Drop through
C     with MARK at 0 if TAIL < FIRST.

      DO 30, I = TAIL, FIRST, -1

C        Is the character a separator?

         DO 10, J = 1, LEN (SEPS)
            IF (STRING (I : I) .EQ. SEPS (J : J)) GO TO 20
   10    CONTINUE

C           Not a separator.  Check for the beginning of the right-most token.

            IF (.NOT. FOUND) THEN
               LAST = I
               FOUND = .TRUE.
            END IF
            GO TO 30

   20    CONTINUE

C           We found a separator.

            IF (FOUND) THEN

C              We just passed the "leading edge" of the desired token.
C              Bail out with a successful search.

               MARK = I + 1
               GO TO 99
            END IF
   30 CONTINUE

C     We reached the FIRST character.  Either it is part of a token, or
C     MARK is still zero.

      IF (FOUND) MARK = FIRST


C     Termination.

   99 RETURN
      END
C+----------------------------------------------------------------------
C
      SUBROUTINE SCLEAN (I)
C
C  PURPOSE:  SCLEAN writes escape sequences which erase the screen of a
C            DEI VT640  "green screen"  graphics terminal. The graphics
C            display  can  be left  in  place or  erased along with the
C            VT100  display.  In both cases  the  terminal is placed in
C            VT100  mode.   A  separate call can  be used to return the
C            terminal to alpha mode.
C
C            <This is obscure.  Did Rosalie mean the I=2 case, which
C            returns the terminal to vector graphics mode?  DAS>
C
C            Logical unit 6 is hard-coded.
C
C  ARGUMENTS:
C   VAR  TYPE  DIM  I/O/S  DESCRIPTION
C    I    I     -     I    1 means  go from  alpha mode  to transparent
C                            mode,  erase all of the VT100 display, but
C                            leave the graphics display with the cursor
C                            set to home.
C                          2 means return to TEK 4014 vector mode.
C                          3 means wipe screen and return to transparent
C                            (i.e. VT100) mode.
C
C  ENVIRONMENT:  DEC VAX/VMS; FORTRAN; VT640 and Selanar termiinals.
C         Also:  SGI IRIS/IRIX; GraphOn GO-200 terminals.
C
C  HISTORY:
C    Oct. 1983    RCL     Original design and coding.
C    Dec. 1985    RCL     Changed SCLEAN(3) to be compatible with 
C                         Selanar terminals as well as VT640's.
C    May 1 '91 D.Saunders Writing integer ASCII codes didn't work on the
C                         IRIS - use CHAR(n) instead of n in the I/O list.
C
C  AUTHOR: Rosalie Lefkowitz, Informatics General Corp.
C
C-----------------------------------------------------------------------

C  *  Arguments:

      INTEGER I

C  *  Execution:

      IF (I .EQ. 1) THEN

C  *     Go from alpha mode to transparent (VT100) mode, erase all of the
C        VT100 display (leaving the graphics display), and home the cursor.
C        The sequence is:       ESC, ESC, FF, CAN, ESC, [, 2, J, ESC, [, H

         WRITE (6, 1001) CHAR (27), CHAR (27), CHAR (12), CHAR (24),
     >      CHAR (27), CHAR (91), CHAR (50), CHAR (74), CHAR (27),
     >      CHAR (91), CHAR (72)

      ELSE IF (I .EQ. 2) THEN

C  *     Return to TEK 4014 vector mode.      
C                        ESC,  ESC, GS

         WRITE (6, 1001) CHAR (27), CHAR (27), CHAR (29)

      ELSE IF (I .EQ. 3) THEN

C  *     Clear graphics screen and return to transparent (VT100) mode.
C                        GS, ESC, FF, CAN

         WRITE (6, 1001) CHAR (29), CHAR (27), CHAR (12), CHAR (24)
      END IF

 999  RETURN

1001  FORMAT ('+', 11A1)
      END
C+----------------------------------------------------------------------
C
      SUBROUTINE STRIPPER (STRING, LENSTR, DELIM)
C
C
C     Description and usage:
C
C        STRIPPER scans a string for enclosing delimiters and passes back
C     the contents (possibly left-shifted one character) along with its
C     length (possibly reduced by 2).  The length of the string on input is
C     assumed to be positive (see the argument description).  If a pair of
C     delimiters is stripped out, the last two characters of the original
C     string are blanked.  Any blanks or tabs in front of the first delimiter
C     are ignored (and retained).  A returned length of zero is possible.
C
C
C        Examples:
C
C        (1)          '12345'           is returned as
C                     12345             (with length reduced by 2).
C                          ^^
C
C        (2)          " "               (a quoted blank) becomes
C                                       (3 blanks with length returned as 1).
C                     ^^^
C
C        (3)          ""                becomes
C                                       (2 blanks but length is returned as 0).
C                     ^^
C
C           STRIPPER was prompted by the need to distinguish blank or
C        apparently numeric text from missing or truly numeric data by "quoting"
C        the text.  The quotes need to be checked for and eliminated if present,
C        with a left shift in place, and a length reduced by 2 (normally).
C
C
C     Arguments:
C
C        Name    Dimension  Type  I/O/S  Description
C        STRING      *       C    I/O    Text string to be searched.  If found,
C                                        one pair of delimiters is eliminated
C                                        on output, and remaining characters
C                                        are left-shifted one position.
C        LENSTR              I    I/O    Length of text string.  LENSTR > 0 is
C                                        assumed as input since the string's
C                                        last significant character should have
C                                        been determined by a prior call to
C                                        GETLINE, SCANx, or equivalent.  If
C                                        enclosing delimiters are found, LENSTR
C                                        is reduced by two (meaning LENSTR = 0
C                                        is a possible output).
C        DELIM       *       C    I      String containing delimiters to be
C                                        searched for, probably ' or " or both.
C                                        If a matching pair is found, the search
C                                        terminates - no further delimiters are
C                                        considered.  A valid pair consists of
C                                        the character STRING (LENSTR:LENSTR)
C                                        and the first SIGNIFICANT character in
C                                        STRING.  (Leading blanks or tabs are
C                                        ignored.)
C
C     External References:
C
C        SCAN2   Used here to find the first significant character.
C
C
C     Environment:  Digital VAX-11/785 VMS FORTRAN (FORTRAN 77).
C                   IMPLICIT NONE is non-standard.
C
C     History:
C
C        14 Mar 1990    R.A.Kennelly    Initial design.
C        15 Mar 1990    M.D.Wong        Initial implementation.
C        02 Oct 1991    D.A.Saunders    Allowed LENSTR=0 upon return,
C                                       where 1 was forced originally.
C
C     Author:  Michael Wong, NASA Ames/Sterling Software, Palo Alto, CA
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments.

      INTEGER
     >   LENSTR
      CHARACTER
     >   DELIM * (*), STRING * (*)

C     Local constants.

      CHARACTER
     >   BLANK
      PARAMETER
     >  (BLANK = ' ')

C     Local variables.

      CHARACTER
     >   SEPS * 2
      INTEGER
     >   HEAD, I, J, MARK, NDELIM, TAIL

C     Execution.


      HEAD = 1
      TAIL = LENSTR
      SEPS = BLANK // CHAR (9)      ! Tab and blank are considered
                                    ! insignificant as leading characters.
      NDELIM = LEN (DELIM)

C     Compare each element in DELIM with first and last significant characters
C     in STRING.  The search is terminated if they both match.

      DO 20, J = 1, NDELIM
         IF (STRING (TAIL : TAIL) .EQ. DELIM (J : J)) THEN

C           The last character is a possible delimiter, so continue search
C           with the first significant character, located by SCAN2.

            CALL SCAN2 (STRING, SEPS, HEAD, TAIL, MARK)

            IF (STRING (HEAD : HEAD) .EQ. DELIM (J : J) .AND.
     >         HEAD .NE. TAIL) THEN        ! Guard against single delimiter.

C              Adjust string length, allowing 0 (where 1 was once forced).

               LENSTR = TAIL - 2

C              Shift left one character, in place, and blank out last two
C              characters.  (The loop does nothing in the degenerate case.)

               DO 10, I = HEAD, LENSTR
                  STRING (I : I) = STRING (I + 1 : I + 1)
   10          CONTINUE
               STRING (TAIL - 1 : TAIL) = BLANK

               GO TO 99
            END IF
         END IF
   20 CONTINUE

   99 RETURN
      END
C+----------------------------------------------------------------------
C
      FUNCTION THREEPT (J, H, DEL)
C
C     One-liner: First derivative, non-central 3-point formula
C     ----------
C
C     Description and usage:
C     ----------------------
C
C        Computes a first derivative approximation for PLSFIT over an
C     interval at a data boundary, using a forward or backward 3-point
C     formula.  The data must be in the form of arrays containing finite
C     difference interval lengths and 2-point forward difference deriva-
C     tives, and the differencing direction is controlled by a flag. See
C     PLSFIT for more details.
C
C        See module BUTLAND for a version with "shape-preserving"
C     adjustments.
C
C     Arguments:
C     ----------
C
C     Name    Type/Dimension  I/O/S  Description
C     J       I               I      Indicates at which end of the
C                                    interval the derivative is to be
C                                    estimated. J = 0 means left-hand
C                                    side, J = 1 means right. 
C
C     H       R (-1:1)        I      Array of interval lengths. The 0th
C                                    element is the length of the interval
C                                    on which the cubic is to be deter-
C                                    mined.
C
C     DEL     R (-1:1)        I      Array of derivative estimates. The
C                                    0th element is the forward difference
C                                    derivative over the interval on which
C                                    the cubic is to be determined.
C                                     
C     THREEPT R                 O    The function value is the derivative.
C
C     Environment:  VAX/VMS; FORTRAN 77
C     ------------
C
C     IMPLICIT NONE and 8-character symbolic names are non-standard.
C
C     Author:  Robert Kennelly, Sterling Federal Systems/NASA-Ames
C     -------
C
C     History:
C     --------
C
C     18 Feb. 1987    RAK    Initial design and coding.
C     06 June 1991    DAS    Original THREEPT renamed BUTLAND; THREEPT
C                            now gives unmodified 1-sided 3-pt. results.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments.

      INTEGER
     &   J
      REAL
     &   H (-1:1), DEL (-1:1), THREEPT

C     Local constants.

      REAL
     &   ONE
      PARAMETER
     &  (ONE = 1.0E+0)

C     Local variables.

      INTEGER
     &   STEP
      REAL
     &   WEIGHT

C     Execution.
C     ----------

C     Estimate first derivative on a left-hand boundary using a 3-point
C     forward difference (STEP = +1), or with a backward difference for
C     the right-hand boundary (STEP = -1).

      STEP = 1 - J - J   ! J here is consistent with related modules.

C     In {H, DEL} form, the derivative looks like a weighted average.

      WEIGHT  = -H (0) / (H (0) + H (STEP))
      THREEPT = WEIGHT * DEL (STEP) + (ONE - WEIGHT) * DEL (0)

C     Termination.
C     ------------

      RETURN
      END
C+----------------------------------------------------------------------
C
      SUBROUTINE TOKENS (STRING, NUMBER, LIST)
C
C
C     Description and usage:
C
C           An aid to parsing input data.  The individual "tokens" in a
C        character string are isolated, converted to uppercase, and stored
C        in an array.  Here, a token is a group of significant, contiguous
C        characters.  The following are NON-significant, and hence may
C        serve as separators:  blanks, horizontal tabs, commas, colons,
C        and equal signs.  See SCANNR for details.  (Use TOKEN2 and SCAN2
C        if you need a variable list of delimiters.)
C
C           Processing continues until the requested number of tokens have
C        been found or the end of the input string is reached.
C
C
C     Parameters:
C
C        Name    Dimension  Type  I/O/S  Description
C        STRING              C    I      Input character string to be analyzed.
C        NUMBER              I    I/O    Number of tokens requested (input) and
C                                        found (output).
C        LIST    NUMBER      C      O    Array of tokens, changed to uppercase.
C
C
C     External references:
C
C        Name    Description
C        SCANNR  Finds positions of the first and last significant characters.
C        UPCASE  Converts a string to uppercase.
C
C
C     Environment:  Digital VAX-11/780 VMS FORTRAN (FORTRAN 77).
C
C
C     Notes:
C
C        (1)  IMPLICIT NONE is non-standard.
C
C
C     Author:  Robert Kennelly, Informatics General Corporation.
C
C
C     Development history:
C
C        16 Jan. 1984    RAK    Initial design and coding.
C        16 Mar. 1984    RAK    Revised header to reflect full list of
C                               separators, repaired faulty WHILE clause
C                               in "10" loop.
C        18 Sep. 1984    RAK    Change elements of LIST to uppercase one
C                               at a time, leaving STRING unchanged.
C        12 Mar. 1986    RAK    Cross-referenced TOKEN2 variation.
C
C-----------------------------------------------------------------------


C     Variable declarations.
C     ----------------------

      IMPLICIT NONE

C     Parameters.

      CHARACTER
     >   BLANK
      PARAMETER
     >   (BLANK = ' ')

C     Variables.

      INTEGER
     >   COUNT, FIRST, I, LAST, MARK, NUMBER
      CHARACTER
     >   STRING * (*), LIST (NUMBER) * (*)

C     Procedures.

      EXTERNAL
     >   UPCASE, SCANNR


C     Executable statements.
C     ----------------------

C     WHILE there are tokens to find, loop UNTIL enough have been found.

      FIRST = 1
      LAST = LEN (STRING)

      COUNT = 0
   10 CONTINUE

C        Get delimiting indices of next token, if any.

         CALL SCANNR (STRING, FIRST, LAST, MARK)
         IF (LAST .GT. 0) THEN
            COUNT = COUNT + 1

C           Pass token to output string array, then change case.

            LIST (COUNT) = STRING (FIRST : MARK)
            CALL UPCASE (LIST (COUNT))
            FIRST = MARK + 2
            IF (COUNT .LT. NUMBER) GO TO 10

         END IF


C     Fill the rest of LIST with blanks and set NUMBER for output.

      DO 20 I = COUNT + 1, NUMBER
         LIST (I) = BLANK
   20 CONTINUE

      NUMBER = COUNT


C     Termination.
C     ------------

      RETURN
      END
C+----------------------------------------------------------------------
C
      SUBROUTINE UPCASE( STRING )
C
C PURPOSE:  UPCASE changes all lower case letters in the given
C           character string to upper case.
C
C METHOD:   Each character in STRING is treated in turn.  The intrinsic
C           function INDEX effectively allows a table lookup, with
C           the local strings LOW and UPP acting as two tables.
C           This method avoids the use of CHAR and ICHAR, which appear
C           to behave differently on ASCII and EBCDIC machines.
C
C ARGUMENTS
C    ARG       DIM     TYPE I/O/S DESCRIPTION
C  STRING       *       C   I/O   Character string possibly containing
C                                 some lower-case letters on input;
C                                 strictly upper-case letters on output
C                                 with no change to any non-alphabetic
C                                 characters.
C
C EXTERNAL REFERENCES:
C  LEN    - Returns the declared length of a CHARACTER variable.
C  INDEX  - Returns the position of second string within first.
C
C ENVIRONMENT:  ANSI FORTRAN 77
C
C AUTHOR: Michael Saunders, Systems Optimization Lab., Stanford, 9/10/85
C         (rewrite of version by Hooper/Kennelly, Informatics, 1983)
C
C-----------------------------------------------------------------------

      CHARACTER      STRING * (*)
      INTEGER        I, J
      CHARACTER      C*1, LOW*26, UPP*26
      DATA           LOW /'abcdefghijklmnopqrstuvwxyz'/,
     $               UPP /'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/

      DO 10 J = 1, LEN(STRING)
         C = STRING(J:J)
         IF (C .GE. 'a'  .AND.  C .LE. 'z') THEN
            I = INDEX( LOW, C )
            IF (I .GT. 0) STRING(J:J) = UPP(I:I)
         END IF
   10 CONTINUE

      RETURN
      END
C+**********************************************************************
C
      SUBROUTINE UTCOPY ( NWORDS, RINPUT, OUTPUT )
C
C ACRONYM: UTility for COPYing arrays.
C          --          ----
C PURPOSE: To copy one array to another.  Also appropriate to shifting
C          the elements of an array up or down, meaning that the input
C          and output arrays can overlap.  Can be used for integers or
C          floating point data -- just be sure to get the number of 4-
C          byte words correct.  [Later: this is no longer a good idea;
C          confine usage to REAL words of whatever precision.]
C     
C NOTES:   NWORDS > 0  covers the normal situation,  where there is no
C          worry about overwriting needed values.  NWORDS < 0 was used
C          to indicate the more awkward situation  (requiring starting
C          from the high-subscript end) in order not to change the no.
C          of arguments and hence not affect existing applications. 
C
C PARAMETERS:
C    ARG     DIM   TYPE I/O/S DESCRIPTION 
C   NWORDS    -      I    I   Indicates no. of REAL words to copy:
C                             NWORDS>0 means either the arrays are not
C                                      overlapping, as in
C                                      CALL UTCOPY ( N, X, Y ),
C                                      or elements are to be moved in
C                                      the smaller subscript direction:
C                                      CALL UTCOPY ( N, X(3), X(1) );
C                             NWORDS<0 means shift ABS(NWORDS) words in
C                                      the larger subscript direction:
C                                      CALL UTCOPY ( -N, X(1), X(3) ).
C   INPUT  |NWORDS|  R    I   Array being copied.
C   OUTPUT |NWORDS|  R    O   Copy of RINPUT(*), possibly overlapping.
C    
C HISTORY:  c. 1980   DAS  Original implementation.
C           04/07/06   "   Make the arrays REAL, not INTEGER.
C
C AUTHOR:  David Saunders, Informatics General.
C
C-**********************************************************************

      INTEGER  NWORDS
      REAL     RINPUT(*), OUTPUT(*)

      IF ( NWORDS.LT.0 ) GO TO 20
            I1  =  1
            I2  =  NWORDS
            INC =  1
            GO TO 30
 20      CONTINUE
            I1  = -NWORDS
            I2  =  1
            INC = -1

 30      DO 40 I = I1, I2, INC
            OUTPUT(I) = RINPUT(I)
 40      CONTINUE

      RETURN
      END
C+----------------------------------------------------------------------
C
      SUBROUTINE WINDOW (NDATA, X, Y, X1, X2, Y1, Y2, NBEGIN, NFINAL,
     >   IER)
C
C     One-liner:  Find a segment of curve (X, Y) enclosed by a rectangle
C     ----------
C
C     Description and usage:
C     ----------------------
C
C        WINDOW is a simple method for finding the portions of an (X, Y)
C     curve which lie inside a rectangular window, for the purpose of
C     plotting parametric curves.  We wish to determine which points must
C     be connected in order to save time by not fitting the portions of
C     the curve which will not be visible.  The beginning index returned
C     is thus either the first one (if the initial point lies inside the
C     window) or indicates the point immediately preceding the first point
C     which does lie inside.  The final index is set analogously.  For a
C     segment to be identified as interior, at least one point (X, Y) must
C     lie inside the window.
C
C        Indices NBEGIN and NFINAL must be initialized by the calling
C     routine to delimit the range of interest.  WINDOW resets these on
C     return, and the user must update them before calling WINDOW again.
C     An error flag is set if the number of data points supplied is
C     less than 2, or if no portion of the curve is visible.
C
C     Arguments:
C     ----------
C
C     Name     Type/Dimension  I/O/S  Description
C     NDATA    I               I      Length of the X, Y data arrays.
C
C     X        R (NDATA)       I      Array of input points' first
C                                     components (abscissas).
C
C     Y        R (NDATA)       I      Second components (ordinates).
C
C     X1,X2    R               I      "Side" coordinates of the rectangular
C                                     window of interest.  The values need
C                                     not be ordered - each point's "X"
C                                     coordinate will be checked to see if
C                                     it is between X1 and X2.
C
C     Y1,Y2    R               I      "Top" and "bottom" coordinates of the
C                                     window.
C
C     NBEGIN   I               I/O    Set by caller to index of first
C                                     point to be examined (usually
C                                     just 1).  On output, NBEGIN is
C                                     the index of the first point before
C                                     the window, i.e.,  point NBEGIN + 1
C                                     is inside.  NBEGIN = 1 if the
C                                     curve originates inside.
C
C     NFINAL   I               I/O    Set by caller to index of last
C                                     point to be examined (usually
C                                     just NDATA).  On output, NFINAL is
C                                     the index of the first point beyond
C                                     the window, i.e.,  point NFINAL - 1
C                                     is inside.  NFINAL = NDATA if the
C                                     curve terminates inside.
C
C     IER      I                 O    Error flag:
C                                       0 = normal return (a non-empty
C                                           segment was found);
C                                       1 = NDATA < 2, or NFINAL - NBEGIN
C                                           < 1 (either a blunder, or
C                                           we've run out of data);
C                                       2 = no segment lies within the
C                                           window.
C
C     Significant local variables:
C     ----------------------------
C
C     CURROUT  Status of current point, .TRUE. means outside the window
C     PREVOUT  Status of previous point
C     EDGES    Counts edge-crossings
C     OUTSIDE  Statement function which returns status of a point
C
C     Environment:  Digital VAX-11/780 VMS FORTRAN 77
C     ------------  Apple Macintosh Absoft MacFORTRAN/020
C
C     Notes:
C     ------
C
C     (1)  IMPLICIT NONE, 8-character symbolic names, and use of "!" for
C          comments are non-standard (until FORTRAN-8X).
C
C     (2)  For the time being, we do NOT attempt to handle the case
C          where two adjacent points in the data arrays are both outside
C          the window, but where the fitted curve loops inside.  We also
C          ignore the possiblity that a portion of the fitted curve
C          loops outside and then back - such a segment is considered
C          to be entirely enclosed by the window.  The rationale is that
C          the cost of examining such cases outweighs their benefit,
C          since the first case is presumably very rare, and the second
C          merely means that the interior parts of the reentrant curve
C          segment are fitted with somewhat fewer points than they would
C          be if each portion were fitted separately.  Time will tell if
C          this is acceptable in practice - a next level of sophisti-
C          cation might be to implement a LINEAR clipping algorithm and
C          check for intersection of the successive chords with the
C          window.  (See, for example CLIPPER in Foley & van Dam.)
C
C     Author:  Robert Kennelly, NASA-Ames/Sterling Federal Systems
C     -------
C
C     Development history:
C     --------------------
C
C     12 Mar. 1987    RAK    Initial design and coding.
C     20 Apr. 1987    RAK    Generalized - X1,X2 and Y1,Y2 pairs need not
C                            be ordered.
C
C-----------------------------------------------------------------------

C     Declarations.
C     -------------

      IMPLICIT NONE

C     Arguments.

      INTEGER
     >   IER, NBEGIN, NDATA, NFINAL
      REAL
     >   X1, X2, X (NDATA), Y1, Y2, Y (NDATA)

C     Local variables.

      LOGICAL
     >   CURROUT, PREVOUT
      INTEGER
     >   EDGES, N
      REAL
     >   XLEFT, XRIGHT, YBOTTOM, YTOP

C     Procedures.

      LOGICAL
     >   OUTSIDE
      REAL
     >   XDUMMY, YDUMMY

      OUTSIDE (XDUMMY, YDUMMY) = (XDUMMY .LT. XLEFT .OR.
     >                            XDUMMY .GT. XRIGHT .OR.
     >                            YDUMMY .LT. YBOTTOM .OR.
     >                            YDUMMY .GT. YTOP)

C     Execution.
C     ----------

      IER    = 0
      NBEGIN = MAX (1, NBEGIN)
      NFINAL = MIN (NDATA, NFINAL)

      IF (NDATA .LT. 2 .OR. (NFINAL - NBEGIN) .LT. 0) THEN

C        There is not enought data to draw a line through, or there are
C        no points to be checked.

         IER = +1
         GO TO 990
      END IF

C      ----------------------------------------------------------------
C     |                                                                |
C     |  1 <= NBEGIN <= NFINAL <= NDATA                                |
C     |  2 <= NDATA                                                    |
C     |                                                                |
C      ----------------------------------------------------------------

      XLEFT   = MIN (X1, X2)
      XRIGHT  = MAX (X1, X2)
      YBOTTOM = MIN (Y1, Y2)
      YTOP    = MAX (Y1, Y2)

C     Loop over the data (if there is more than one point), looking for
C     edge crossings.

      CURROUT = OUTSIDE (X (NBEGIN), Y (NBEGIN))
      EDGES   = 0

      DO 10, N = NBEGIN + 1, NFINAL

C        Save the status of the previous point.

         PREVOUT = CURROUT
         CURROUT = OUTSIDE (X (N), Y (N))
         IF (CURROUT .NEQV. PREVOUT) THEN

C           An edge has been crossed, do the bookkeeping.

            EDGES = EDGES + 1
            IF (CURROUT) THEN

C              It's the end of the line as soon as a transition is made
C              from inside to outside.  Set "right-hand" index and quit.

               NFINAL = MIN (NDATA, N + 1)
               GO TO 990
            ELSE

C              We just entered the visible region from outside.  Set
C              "left-hand" index and continue looping to look for an
C              exit transition.

               NBEGIN = N - 1
            END IF
         END IF
   10 CONTINUE

C     Clean up some special cases.

      IF (EDGES .EQ. 0 .AND. CURROUT) THEN

C        Evidently the entire curve lies outside the window.

         IER = +2
         GO TO 990
      ELSE IF (EDGES .EQ. 1 .AND. .NOT.CURROUT) THEN

C        The final point is inside.

         NFINAL = NDATA
      END IF

C     Termination.
C     ------------

  990 CONTINUE
      RETURN
      END
C+----------------------------------------------------------------------
C
      SUBROUTINE WRAPPER (MAXPT, NPT, CURVE, X, Y, IER)
C
C
C     One-liner:  Checks/corrects match of first and last points of curve
C
C
C     Description and usage:
C
C           A hide-the-ugly-details utility which checks for proper wrap-
C        around data format for X, Y curves stored in a packed array.  For
C        fitting with parametric piecewise cubics (PLSFIT), the first and
C        last points must agree.  WRAPPER either verifies that this is true
C        or inserts a point (if possible).  An error flag is set if there
C        is no more room in the data arrays, in which case the calling
C        routine will just have to cope (for example, by ignoring the
C        request for periodic end conditions).
C
C           WRAPPER was written for use by QPLOT, a general-purpose data
C        graphics program developed at NASA-Ames.  It could find other
C        application, but may not pay unless the data is actually packed.
C
C
C     Arguments:
C
C        Name    Dimension  Type  I/O/S  Description
C        MAXPT              I     I      Declared size of data arrays in
C                                        calling routine (required for
C                                        error checking).
C        NPT     CURVE      I     I      Point counts for each curve.
C        CURVE              I     I      Index of curve to be checked.
C        X       MAXPT      R     I/O    Array of abscissas to be checked
C                                        for wrap-around condition.  If
C                                        storage permits (determined from
C                                        MAXPT), an extra point will be
C                                        added to the end of the data if
C                                        required.
C        Y       MAXPT      R     I/O    Array of ordinates (see X).
C        IER                I       O    Error flag.  0 means all's well,
C                                        and 1 signals that the first and
C                                        last points didn't match and that
C                                        no recovery was possible.
C
C
C     Environment:  Digital VAX-11/780 VMS FORTRAN 77
C
C
C     Notes:
C
C        (1)  IMPLICIT NONE and 8-character symbolic names are non-standard.
C
C
C     Author:  Robert Kennelly, NASA-Ames/Sterling Software
C
C
C     Development history:
C
C         7 Apr. 1987    RAK    Initial design and coding.
C
C-----------------------------------------------------------------------


C     Declarations.
C     -------------

      IMPLICIT NONE

C     Arguments.

      INTEGER
     >   CURVE, IER, MAXPT, NPT (CURVE)
      REAL
     >   X (MAXPT), Y (MAXPT)

C     Local variables.

      INTEGER
     >   FIRST, I, LAST


C     Execution.
C     ----------

      IER = 0
      IF (MAXPT .LE. 0) STOP 'WRAPPER:  Array dimension error!'

C     Count the points in each curve up to this one.

      LAST = 0
      DO 10, I = 1, CURVE
         LAST = LAST + NPT (I)
   10 CONTINUE
      IF (LAST .LE. 0) STOP 'WRAPPER:  No data!'

C     Set index for the begining of this curve's data points.

      FIRST = LAST - NPT (CURVE) + 1

C     Does the data wrap around?

      IF (X (FIRST) .NE. X (LAST) .OR. Y (FIRST) .NE. Y (LAST)) THEN
         IF (LAST .LT. MAXPT) THEN

C           Patch up the data arrays and adjust the points array.

            LAST = LAST + 1
            X (LAST) = X (FIRST)
            Y (LAST) = Y (FIRST)
            NPT (CURVE) = NPT (CURVE) + 1
         ELSE

C           Oops - there is no room for the required point to be added!

            IER = +1
         END IF
      END IF


C     Termination.
C     ------------

      RETURN
      END

C+------------------------------------------------------------------------------
C
      REAL FUNCTION XDIMTB (TEXT, LENTXT, LINE1, LINE2)
C
C ONE LINER:  Calculates length of a text block in inches. (Level 1 - 3)
C
C PURPOSE:
C    XDIMTB finds the length in inches of an imaginary rectangle enclosing a
C    block of text.  This routine may be used as an alternative to the DISSPLA
C    subroutine XSTORY, for handling data in character array format.
C
C    The character size and font are assumed to be set outside the routine
C    and hence constant within the indicated lines of text.
C
C    XDIMTB was introduced as when it was realized that the SMDLIB utility
C    TXTLEG needed such a function, both prior to calling TXTLEG (for 
C    positioning the legend) and within TXTLEG (for the optional box).  The
C    translation of TXTLEG to DISSPLA required the similar translation of 
C    XDIMTB.
C
C ARGUMENTS:
C     ARG      DIM     TYPE      I/O/S      DESCRIPTION
C     TEXT    (*) * (*)   C        I       Character array containing text.
C                                          (Terminating '$'s may be required,
C                                          depending on the usage of 
C                                          LENTXT (*).)  Elements LINE1:LINE2 
C                                          will be displayed in the current 
C                                          color (as opposed to the color 
C                                          associated with each line/symbol).
C     LENTXT     (*)      I        I       Array containing number of characters
C                                          in text line(s), if known.  (No 
C                                          trailing '$'s are needed in this 
C                                          case.)  Otherwise, pass 100 as the 
C                                          first (and only) value, and all lines
C                                          of text will be self counted and all
C                                          lines of text will be self counted
C                                          (requiring the trailing '$'s).
C    LINE1      -        I         I        First line of text to process.
C    LINE2      -        I         I        Last line of text to process.
C
C METHOD:
C    The size of the block is the length of the longest line in the block.
C    DISSPLA's XMESS utility, with its self counting option, is used for
C    finding the length of a string.
C
C ERROR HANDLING:  No error message if called at wrong level.
C
C EXTERNAL REFERENCES:
C    XMESS      Finds length of text string in inches
C
C ENVIRONMENT:  VMS/VAX; FORTRAN 77
C
C HISTORY:
C    04/19/89    M.D. Wong        Initial design and implementation for SMDLIB
C    06/07/89    M.D. Wong        Adapted for DISSPLA from SMDLIB.  Error 
C                                 handling removed (since internal common blocks
C                                 are avoided).
C    08/29/89    M.D. Wong        Updated to run on DISSPLA version 11.0.
C                                 (%REF taken out of call to function XMESS.)
C    02/26/90    M.D. Wong        Added LENTXT to argument list to avoid 
C                                 reliance on DISSPLA self counting option.
C    04/12/90    M.D. Wong        Enabled 100 to be passed as first and only
C                                 element of the LENTXT array to enable self 
C                                 counting option.
C
C AUTHOR:  Michael Wong, Sterling Software, Palo Alto, CA
C
C-------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments
C     ---------

      CHARACTER TEXT (*) * (*)
      INTEGER   LENTXT (*), LINE1, LINE2

C     Local variables
C     ---------------

      INTEGER   I, LENGTH

C     Procedures
C     ----------

      REAL      XMESS
      EXTERNAL  XMESS

C     Execution
C     ---------

      XDIMTB = 0.
      DO 100, I = LINE1, LINE2
         IF (LENTXT (1) .EQ. 100) THEN
            LENGTH = LENTXT (1)
         ELSE
            LENGTH = LENTXT (I)
         END IF
         XDIMTB = MAX (XDIMTB, XMESS (TEXT (I), LENGTH))
  100 CONTINUE

      RETURN
      END
