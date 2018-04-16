C+------------------------------------------------------------------------------
C
      SUBROUTINE DISTRIB (METHOD, QUERY, NX, XA, XB, LUNCRT, LUNKBD,
     >   IPARAM, RPARAM, START, NWK, WK, X, DETAILS, IER)
C
C  PURPOSE
C  -------
C
C        DISTRIB is a high-level routine for driving the available 1-D grid
C     point distribution utilities.  It was modularized from the bulk of
C     the in-line code of program DISTRIBUTE when application of a choice
C     of methods to distributions along arcs in 3-space was required.
C     (See ARCDIS for the latter application.)
C
C  METHOD
C  ------
C
C        Provision must be made for repeated application of the current
C     method with the current parameters (as in the ARCDIS application)
C     and for showing details of the method (as needed for the plots of
C     DISTRIBUTE results).  The approach taken is to have all prompts
C     either enabled or disabled (argument QUERY) and to go with current
C     values in METHOD, IPARAM (*) and RPARAM (*) if QUERY is false.
C
C        Some methods require additional work-space, but a given application
C     may have no intention of using certain methods.  Therefore the calling
C     program should supply only enough for the methods that will be selected.
C
C        Error handling is problematic.  If a method fails, the user is
C     prompted as to whether to retry (go back to the menu), proceed anyway,
C     (since some iterative methods may still have given a usable result),
C     or quit.
C
C  ARGUMENTS
C  ---------
C
C     ARG     DIM    TYPE   I/O/S   DESCRIPTION
C
C     METHOD   -       I    I/O     Menu item number - see menu below.
C                                   If QUERY=T: METHOD is input with a default;
C                                               output with chosen item #.
C                                   If QUERY=F: METHOD is input with item # from
C                                               an earlier call; unchanged on
C                                               output.
C     QUERY    -       L    I       Used to suppress prompts - see METHOD. Also:
C                                   If QUERY=T: relevant IPARAM(*) and RPARAM(*)
C                                               choices are output.
C                                   If QUERY=F: relevant IPARAM(*) and RPARAM(*)
C                                               values are assumed to be input.
C     NX       -       I    I       Requested number of points in distribution.
C     XA       -       R    I       [XA, XB] is the interval being distributed;
C     XB                            normally XA < XB, but XA > XB may be OK.
C     LUNCRT   -       I    I       Logical unit for prompts/diagnostics.
C     LUNKBD   -       I    I          "     "    "  keyboard entries.
C     IPARAM   *       I    I/O     INTEGER parameter(s) used by some methods.
C                                   Usually no more than 1 - see summary below.
C                                   (Normally unchanged upon return if QUERY=F.)
C     RPARAM   *       R    I/O     REAL parameters used by certain methods.
C                                   0 to 4 values - see summary below.
C                                   Any NEGATIVE value intended as a "delta X"
C                                   indicates a PERCENTAGE of (XB - XA) rather
C                                   than an absolute increment.
C                                   (Normally unchanged upon return if QUERY=F.)
C     START    -      C*1   I       Provides for entering a starting guess
C                                   in X (*) for certain iterative methods.
C                                   'A' for 'AUTO' means the method should
C                                       self-start (typically with uniform Xs);
C                                   'I' for 'INPUT' means a starting guess is
C                                       supplied to DISTRIB in X (1:NX).
C     NWK      -       I    I       Amount of work-space supplied for methods
C                                   that need it - see summary below.  NWK=0 is
C                                   OK if WK (*) is unused.  NWK=1 for several
C                                   methods involving dX > 0 in RPARAM (*).
C     WK       *       R    I / S   Work-space needed by some methods.  Length
C                                   should be at least NWK if it is needed,
C                                   else any valid location passed here will do.
C                                   For some methods, WK (1) should be input
C                                   with a scale factor to divide dX > 0 values
C                                   by in RPARAM (*).
C     X       NX       R    I/O     Output with desired distribution on [XA,XB].
C                                   Maybe input with a starting guess for
C                                   iterative methods.  See START.
C     DETAILS  -  C*(LENMENU) O     Output with details of chosen method;
C                                   ignored if QUERY=F.  Normally just the
C                                   appropriate menu line, but methods which
C                                   use further parameters will show these
C                                   following the appropriate subroutine name.
C                                   LENMENU is defined below.
C     IER      -       I      O     IER = 0 means no error;
C                                   IER =-1 means quit (normal or abnormal,
C                                           depending on the application);
C                                   IER =+1 means METHOD was invalid;
C                                   IER = 2 means NWK was too small;
C                                   IER = 3 means the specified method failed.
C
C  CONTROL PARAMETER & WORK-SPACE REQUIREMENTS
C  -------------------------------------------
C
C      METHOD   Description                               IPARAM/RPARAM/WK usage
C                                                         -      -
C  (1) DSTRIB:  Uniform distribution                         No parameters
C  (2) DSTRIB:  Simple sinusoidal bunching at end point(s)   I(1): mode 
C  (3) DSTRIB:  Modified sinusoidal bunching at end pt(s)    I(I): mode
C                                                            R(1): exponent
C  (4) DSTRIB:  Modified sinusoidal bunching at interior pt  I(1): mode
C                                                            R(1): exponent
C                                    Kludge to return the    R(2): interior pt.
C                                    corresponding index --> R(3): index (outp.)
C  (5) EXPDIS:  Exponential-type distribution (given Beta)   I(1): mode
C                                                            R(1): Beta (> 1.)
C  (6) EXPDIS4: Exponential-type (given smallest interval)   I(1): mode
C                                                            R(1): dX at XA | XB
C                                                            R(2): corresp. Beta
C  (7) CONDIS:  Exponential-type (constrained interior pt M) I(1): subscript M
C                                                            R(1): interior X(M)
C  (8) ARBDIS:  Arbitrary-shape distribution                 No parameters.
C                                              NWK = 5 * NX + 2 * # control pts.
C  (9) GEODIS:  Geometric-type (given initial dX; 1-sided)   R(1): dX (left)
C                                                            R(2): exponent
C (10) GEODIS2: Geometric-type (given initial dX; 2-sided)   R(1): dX (left)
C                                                            R(2): expt. (left)
C                                                            R(3): expt. (right)
C (11) HTDIS4:  Hyperbolic tangent method (end-pt. slopes)   R(1): slope (left)
C                                                            R(2): slope (right)
C (12) HTDIS4:  Hyperbolic tangent method (dX at each end)   R(1): dX (left)
C                                                            R(2): dX (right)
C
C (13) FOILGRID: Hybrid cosine/quadratic suited to airfoils; R(1): Sine term wt.
C                heuristic adjustment at leading edge              in [0.,1.]
C
C (14) FOILGRD:  Hybrid linear/quadratic/sine/cosine suited  R(1): Linear wt.
C                to airfoils; all weights in [0.,1.]         R(2): Quadratic wt.
C                                                            R(3): Sine wt.
C                                                            R(4): Cosine wt.
C
C (15) BLGRID:   Hybrid geometric/Vinokur suited for         I(1): NBLAYER [20]
C                boundary layers                             R(1): dX (left)
C                                                            R(2): dX (right)
C                                                            R(3): RBLAYER [1.1]
C (16) STRETCH:  One-sided, given smallest interval and      R(1): dX (left)
C                an estimate of stretching parameter Beta    R(2): corresp. Beta
C
C (17) SHOCKGRID:  2-sided + 1-sided combination intended    I(1): NMARGIN [3]
C                  to resolve a boundary layer and a shock   R(1): D1 (left)
C                  with NMARGIN points outside the shock     R(2): DE (sh. edge)
C                  location defined by SEDGE                 R(3): SEDGE
C
C (18) EXPDIS5:  One-sided, but possibly decreasing spacing  I(1): mode
C                                                            R(1): dX at XA | XB
C  PROCEDURES
C  ----------
C
C     See above for the distribution utilities.  Also:
C     READER   Prompting utility
C     SELECT   Main menu utility
C
C  FILES USED
C  ----------
C
C     See LUNCRT and LUNKBD arguments.
C
C  ENVIRONMENT
C  -----------
C
C     FORTRAN 77 with extensions:
C     >  IMPLICIT NONE
C     >  A few names longer than 6 characters
C     >  ! comments
C
C  HISTORY
C  -------
C
C     03/16/89  R.Langhi    Original reuse of pieces of DISTRIBUTE in REDIS,
C                           but this aggravated the maintenance problem.
C     08/09/89  D.Saunders  Modularized main part of DISTRIBUTE for use by both
C                           it and REDIS (renamed ARCDIS) to save maintenance.
C                           (Possibly more trouble than it was worth.)
C     11/19/90   "    "     HTDIS installed (only HTDIS2 originally, because
C                           EXPDIS2 offers the same functionality; but HTDIS
C                           has its place (more self-contained)).
C     04/16/91   "    "     Arranged for -dX entries to mean +ve percentages.
C     04/23/91   "    "     Arranged for scaling +dX values (dividing by
C                           WK (1)), as required by the ARCDIS application.
C     05/04/91   "    "     Had to shorten the legend info. from GEODIS2.
C     06/04/91   "    "     HTDIS2 has a "precise dXs" or "match slopes"
C                           choice, wedged in via -dX ... awkward.
C     09/10/91   "    "     Added second menu item for the "match slopes" case
C                           of HTDIS2, whose DSINPUT arg. replaces -dX kludge.
C                           Removed HTDIS since HTDIS2 (or EXPDIS2) is better.
C     09/16/91   "    "     When a method failed, but was probably close,
C                           and QUERY was true, and "proceed" was requested,
C                           DETAILS was not being set up, causing DISTRIBUTE
C                           driving program to fail when adding to DETAILS.
C     05/13/93   "    "     Added THREEINC (composite use of HTDIS2 to control
C                           the max./min interior dX as well as the end dXs.
C     02/27/95   "    "     Safeguard for EXPDIS2 should divide by NX-1, not NX.
C     03/30/95   "    "     Eliminated THREEINC (it never worked); installed
C                           FOILGRID.
C     08/29/96   "    "     Installed FOILGRD.
C     09/16/96   "    "     Installed BLGRID.
C     03/22/04   "    "     Installed STRETCH.  Use EXPDIS4 and HTDIS4 and
C                           compile and link with 64-bit precision (only).
C     07/03/06   "    "     Installed SHOCKGRID.
C     08/03/09   "    "     Installed EXPDIS5.
C
C  AUTHOR:   David Saunders, Sterling Software/NASA Ames, Mt. View, CA.
C
C-------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER
     >   IER, IPARAM (*), LUNCRT, LUNKBD, METHOD, NWK, NX
      REAL
     >   RPARAM (*), X (NX), XA, XB, WK (*)
      LOGICAL
     >   QUERY
      CHARACTER
     >   START * 1,
     >   DETAILS * (*)        ! Should be same length as MENU entries below.

C     Local constants:

      INTEGER
     >   LENMENU, MXMENU
      REAL
     >   HALF, ONE, PERCENT, TWOTON, ZERO
      CHARACTER
     >   BLANK * 1, SUBNAME * 7

      PARAMETER
     >  (BLANK   = ' ',
     >   HALF    = 0.5E+0,
     >   LENMENU = 60,   ! Applies to both MENU and DETAILS
     >   MXMENU  = 18,   ! Max. number of distribution options handled
     >   ONE     = 1.E+0,
     >   PERCENT = 1.E-2,
     >   SUBNAME = 'DISTRIB',
     >   TWOTON  = 200.E+0,
     >   ZERO    = 0.E+0)

C     Local variables:

      INTEGER
     >   I, LUNOUT, M, MODE, NC, NEEDED
      REAL
     >   ADJUST, BETA, DX, DX2, PCTRANGE, POWER, POWER2, RANGE, SCALE,
     >   XM, XTEMP, YTEMP
      LOGICAL
     >   AUTO, DEFAULT, DSINPUT, QUERYLOC, QUIT, YES
      CHARACTER
     >   ANSWER*1, DISNAME*8, STARTLOC*1, MENU (MXMENU) * (LENMENU)

C     Procedures:

      EXTERNAL
     >   ARBDIS, CONDIS, DSTRIB, EXPDIS, EXPDIS4, EXPDIS5, FOILGRD,
     >   FOILGRID, GEODIS, GEODIS2, HTDIS4, READC, READI, READR, READY,
     >   SELECT, SHOCKGRID, STRETCH

C     Storage:

      DATA
     >   MENU
     >     /' (1) DSTRIB:  Uniform distribution',
     >      ' (2) DSTRIB:  Simple sinusoidal bunching at end point(s)',
     >      ' (3) DSTRIB:  Modified sinusoidal bunching at end pt(s)',
     >      ' (4) DSTRIB:  Modified sinusoidal bunching at interior pt',
     >      ' (5) EXPDIS:  Exponential-type distribution (given Beta)',
     >      ' (6) EXPDIS4: Exponential-type (1st or last dX specified)',
     >      ' (7) CONDIS:  Exponential-type (constrained interior pt.)',
     >      ' (8) ARBDIS:  Arbitrary-shape distribution',
     >      ' (9) GEODIS:  Geometric-type (given initial dX; 1-sided)',
     >      '(10) GEODIS2: Geometric-type (given initial dX; 2-sided)',
     >      '(11) HTDIS4:  Vinokur''s TANH-based method (end slopes)',
     >      '(12) HTDIS4:  Vinokur''s method iterated for precise dXs',
     >      '(13) FOILGRID Sine/linear/quadratic, suited to airfoils',
     >      '(14) FOILGRD  Linear/quadratic/sine/cosine, for airfoils',
     >      '(15) BLGRID:  Geometric/Vinokur, suited to boundary layer',
     >      '(16) STRETCH: One-sided (initial dX; Beta guess updated)',
     >      '(17) SHOCKGRD Two-sided + 1-sided for b.layer + shock',
     >      '(18) EXPDIS5: EXPDIS4 variant to allow decreasing spacing'
     >     /

C                 ^ Code below assumes method names start with 6th character,
C                 while DISNAME*8 above limits the name of the distribution.


C     Execution:

      IF (METHOD .LT. 1 .OR. METHOD .GT. MXMENU) GO TO 900

      SCALE = ONE       ! ARCDIS (higher level) needs to scale the dXs
      IF (NWK .GT. 0) THEN    ! Kludge to avoid changing calling sequence
         IF (WK (1) .GT. ZERO) SCALE = ONE / WK (1)
      END IF

      RANGE = XB - XA
      PCTRANGE = -PERCENT * RANGE
      DISNAME = MENU (METHOD) (6 :)
      QUERYLOC = QUERY

  200 IF (QUERYLOC) THEN

         CALL SELECT ('Select a distribution by number or name.',
     >      MXMENU, MENU, .FALSE., LUNCRT, LUNKBD, METHOD, DISNAME,
     >      QUIT)
         IF (QUIT) GO TO 990

         DETAILS = MENU (METHOD) (6 :)   ! For methods with no further details

         IF ((METHOD .GE. 5 .AND. METHOD .LT. 13) .OR.
     >        METHOD .EQ. 18) THEN
            YES = .FALSE.
            CALL READY (LUNCRT,
     >         'Display iteration history? (Y/N; <CR>=No) ',
     >         LUNKBD, YES, DEFAULT, QUIT)
            IF (YES) THEN
               LUNOUT = LUNCRT
            ELSE
               LUNOUT = -LUNCRT
            END IF
         END IF
      ELSE
         LUNOUT = -LUNCRT                ! Suppress iteration history
      END IF

C     The following cases contain duplicate code in an attempt to
C     keep each case self-contained...

      GO TO (410, 420, 430, 440, 450, 460, 470, 480, 490, 500, 520,
     >       520, 530, 540, 550, 560, 570, 580) METHOD

  410    CONTINUE    ! DSTRIB: Uniform distribution

            CALL DSTRIB (0, 1, RPARAM, NX, XA, XB, X)
            GO TO 800

  420    CONTINUE    ! DSTRIB: Simple sinusoidal distribution - bunch at end(s)

            IF (QUERYLOC) THEN
               MODE = 1
               ANSWER = 'L'
               CALL READC (LUNCRT,
     >      'Bunch at low end? high end? or both? [L/H/B; <CR> = Low] ',
     >            LUNKBD, ANSWER, DEFAULT, QUIT)
               IF (ANSWER .EQ. 'H') MODE = -1
               IF (ANSWER .EQ. 'B') MODE = 2
               IPARAM (1) = MODE
            END IF

            CALL DSTRIB (IPARAM, 1, ONE, NX, XA, XB, X)
            GO TO 800

  430    CONTINUE    ! DSTRIB: Modified sinusoidal - bunch at end point(s)

            IF (QUERYLOC) THEN
               MODE = 1
               ANSWER = 'L'
               CALL READC (LUNCRT,
     >      'Bunch at low end? high end? or both? [L/H/B; <CR> = Low] ',
     >            LUNKBD, ANSWER, DEFAULT, QUIT)
               IF (ANSWER .EQ. 'H') MODE = -1
               IF (ANSWER .EQ. 'B') MODE = 2
               IPARAM (1) = MODE

               WRITE (LUNCRT, 1010)
               POWER = HALF
               CALL READR (LUNCRT,
     >            'Enter fractional exponent > 0. Default is 0.5: ',
     >            LUNKBD, POWER, DEFAULT, QUIT)
               IF (POWER .LE. ZERO) GO TO 430
               RPARAM (1) = POWER
               DETAILS = DISNAME
               WRITE (DETAILS (11 : 34), '(A, I2, A, F5.3)')
     >            'Mode = ', MODE, '  Power = ', POWER
            END IF

            CALL DSTRIB (IPARAM, 1, RPARAM, NX, XA, XB, X)
            GO TO 800

  440    CONTINUE    ! DSTRIB: Modified sinusoidal - bunch at interior pt.

            IF (QUERYLOC) THEN
               RPARAM (2) = (XA + XB) * HALF
               CALL READR (LUNCRT,
     >            'Interior pt. at which to bunch? [<CR>=mid-pt.] ',
     >            LUNKBD, RPARAM (2), DEFAULT, QUIT)
               IF (RPARAM (2) .LE. XA .OR. RPARAM (2) .GE. XB) GO TO 440

               WRITE (LUNCRT, 1010)
               IPARAM (1) = 3
               RPARAM (1) = ONE
               CALL READR (LUNCRT,
     >            'Enter fractional exponent > 0. Default is 1.0: ',
     >            LUNKBD, RPARAM (1), DEFAULT, QUIT)
               IF (RPARAM (1) .LE. ZERO) GO TO 440

               DETAILS = DISNAME
               WRITE (DETAILS (11 : 33), '(A, I1, A, F5.3)')
     >            'Mode = ', MODE, '  Power = ', RPARAM (1)
            END IF

            CALL DSTRIB (IPARAM, 2, RPARAM, NX, XA, XB, X)
                                  ! RPARAM (3) is returned with the index
                                  ! of the interior pt., as a REAL (kludge)
            GO TO 800

  450    CONTINUE    ! EXPDIS: Exponential-type distribution (given Beta)

            IF (QUERYLOC) THEN
               MODE = 3
               ANSWER = 'L'
               CALL READC (LUNCRT,
     >      'Bunch at low end? high end? or both? [L/H/B; <CR> = Low] ',
     >            LUNKBD, ANSWER, DEFAULT, QUIT)
               IF (ANSWER .EQ. 'H') MODE = 1
               IF (ANSWER .EQ. 'B') MODE = 2
               IPARAM (1) = MODE

  453          BETA = 1.01E+0
               CALL READR (LUNCRT,
     >            'Enter stretching parameter Beta > 1.0.  [1.01]: ',
     >            LUNKBD, BETA, DEFAULT, QUIT)
               IF (BETA .LE. ONE) GO TO 453
               RPARAM (1) = BETA
            END IF

            CALL EXPDIS (IPARAM, RPARAM, NX, X)   ! X comes back in [0., 1.]

            DO 455, I = 1, NX
               X (I) = XA + RANGE * X (I)
  455       CONTINUE

            IF (QUERY) THEN
               DX = X (2) - X (1)
               IF (IPARAM (1) .EQ. 1) DX = X (NX) - X (NX-1)
               DETAILS = DISNAME
               WRITE (DETAILS (11 : 44), '(A, F9.6, A, 1P, E10.3)')
     >            'Beta =', RPARAM (1), '  dX(1) =', X (2) - X (1)
            END IF
            GO TO 800

  460    CONTINUE    ! EXPDIS4: Exponential-type (given smallest interval)

            IF (QUERYLOC) THEN
               MODE = 3
               ANSWER = 'L'
            
               CALL READC (LUNCRT,
     >      'Bunch at low end? high end? or both? [L/H/B; <CR> = Low] ',
     >            LUNKBD, ANSWER, DEFAULT, QUIT)
               IF (ANSWER .EQ. 'H') MODE = 1
               IF (ANSWER .EQ. 'B') MODE = 2
               IPARAM (1) = MODE

  465          DX = RANGE / TWOTON     ! DX << range / (N - 1)

               WRITE (LUNCRT, 1002)
     >            ' Desired smallest interval? (dX > 0 is absolute;',
     >            ' dX < 0 is relative)'
               CALL READR (LUNCRT,
     >           'dX = -r means r% of range; default = 0.5% of range: ',
     >            LUNKBD, DX, DEFAULT,QUIT)
               IF (DX .EQ. ZERO .OR. DX .GE. RANGE / (NX - 1)) GO TO 465
               RPARAM (1) = DX
            END IF

            IF (RPARAM (1) .LT. ZERO) THEN
               ADJUST = PCTRANGE
            ELSE
               ADJUST = SCALE
            END IF

            DX = RPARAM (1) * ADJUST

            CALL EXPDIS4 (IPARAM, XA, XB, DX, NX, X, BETA, LUNCRT)

            RPARAM (2) = BETA           ! Output
            IF (QUERY) THEN
               DETAILS = DISNAME
               IF (MODE /= 1) THEN
                  WRITE (DETAILS (11:44), '(A, 1P, E10.3, A, 0P, F9.6)')
     >               'dX(1) =', X (2) - X (1), '  Beta =', BETA
               ELSE
                  WRITE (DETAILS (11:44), '(A, 1P, E10.3, A, 0P, F9.6)')
     >               'dX(N) =', X (NX) - X (NX - 1), '  Beta =', BETA
               END IF
            END IF
            GO TO 800

  470    CONTINUE    ! CONDIS: Exponential-type (constrained interior pt.)

            IF (QUERYLOC) THEN
  473          XM = (XA + XB) * HALF       ! No particularly good default
               CALL READR (LUNCRT,
     >            'Interior X to be included in the distribution? ',
     >            LUNKBD, XM, DEFAULT, QUIT)
               IF (XM .LE. XA .OR. XM .GE. XB) GO TO 473

  475          M = (NX + 1) / 2
               CALL READI (LUNCRT, 'Corresponding desired subscript? ',
     >            LUNKBD, M, DEFAULT, QUIT)
               IF (M .LE. 1 .OR. M .GE. NX) GO TO 475

               IPARAM (1) = M
               RPARAM (1) = XM
            END IF

            CALL CONDIS (NX, IPARAM, XA, RPARAM, XB, X)

            IF (QUERY) THEN
               DETAILS = DISNAME
               WRITE (DETAILS (11 : 35), '(A, 1P, E10.3, A, I4)')
     >            'X(M) =', X (M), '  M =', IPARAM (1)
            END IF
            GO TO 800

  480    CONTINUE    ! ARBDIS: Arbitrary-shape distribution

C           <Presently arranged for interactive use only.  Entering
C           the control distribution is problematic otherwise.>

            WRITE (LUNCRT, '(/, (A))')
     >   ' Enter nominal control points (XC,YC) on [XA, XB], monotonic',
     >         ' and preferably including end points XA and XB.',
     >         ' YC > 0 is proportional to the desired spacing at XC.',
     >         ' YC set all the same <-> uniform distribution;',
     >         ' smaller YC produces smaller spacing (bunching).',
     >         ' CTRL/Z means no more.',
     >         BLANK

  482       NC = 0
  484       CONTINUE
               WRITE (LUNCRT, 1003, ADVANCE='NO')
     >            ' Enter (XC,YC) for control pt. # ', NC + 1, ': '
               READ (LUNKBD, *, ERR=484, END=485) XTEMP, YTEMP

               IF (XTEMP .LT. XA .OR. XTEMP .GT. XB) THEN
                  WRITE (LUNCRT, '(/, A)') ' Bad abscissa - try again.'
                  GO TO 484
               END IF

C              Using a single work array is a little inconvenient.
C              Put the YCs at the end of WK(*) till we know how many:

               NC = NC + 1
               WK (NC) = XTEMP
               WK (NWK + 1 - NC) = YTEMP
               IF (NC .LT. NWK / 2)
     >      GO TO 484

  485       CONTINUE
            IF (NC .LT. 2) THEN
               WRITE (LUNCRT, '(A)')
     >            ' You need at least 2 control points.  Start again.'
               GO TO 482
            END IF

C           Move the YC values to follow the XCs:

            DO 486, I = 1, NC
               WK (NC + I) = WK (NWK + 1 - I)
  486       CONTINUE

C           Provide for hand-entering a starting solution (test purposes).
C           START = 'A' already provides for the application to supply it.

            STARTLOC = START
            IF (STARTLOC .NE. 'I') THEN
               IF (NX .LE. 20) THEN      ! Don't hand enter too many...
                  AUTO = .TRUE.
                  WRITE (LUNCRT, '(/, A)')
     >    ' You may hand-enter a starting guess for X(*) or default it.'
                  CALL READY (LUNCRT,
     >          'Do you want the iteration to self-start?  [<CR>=Yes] ',
     >               LUNKBD, AUTO, DEFAULT, QUIT)

                  IF (AUTO) THEN
                     STARTLOC = 'A'
                  ELSE
                     STARTLOC = 'I'
                     WRITE (LUNCRT, '(/, (A))')
     >                  ' X(1) and X(N) should match XA and XB.', BLANK
                     I = 1
  488                CONTINUE
                        WRITE (LUNCRT, 1003, ADVANCE='NO')
     >                     ' Enter X(', I, '): '
                        READ (LUNKBD, *, ERR=488) X (I)
                        I = I + 1
                        IF (I .LE. NX)
     >               GO TO 488
                  END IF
               END IF
            END IF
      
            NEEDED = 2 * NC + 5 * NX
            IF (NWK .LT. NEEDED) GO TO 910  ! Not enough work-space

C           Generate distribution defined by the control points and NX:

            CALL ARBDIS (NX, XA, XB, NC, WK (1), WK (NC + 1), STARTLOC,
     >          LUNOUT, WK (NC + NC + 1), X, IER)
            IF (IER .NE. 0) GO TO 920
            GO TO 800

  490    CONTINUE    ! GEODIS: Modified geometric distribution (1-sided)

            IF (QUERYLOC) THEN
               DX = RANGE / TWOTON
               WRITE (LUNCRT, 1002)
     >            ' First increment? (dX > 0 is absolute;',
     >            ' dX < 0 is relative)'
               CALL READR (LUNCRT,
     >           'dX = -r means r% of range; default = 0.5% of range: ',
     >            LUNKBD, DX, DEFAULT,QUIT)
               IF (DX .EQ. ZERO .OR. DX .GE. RANGE) GO TO 490
               RPARAM (1) = DX

               WRITE (LUNCRT, 1020)
               POWER = ZERO
               CALL READR (LUNCRT,
     >            'Clustering exponent?  (<CR> gives 0.) ',
     >            LUNKBD, POWER, DEFAULT, QUIT)

               RPARAM (1) = DX
               RPARAM (2) = POWER
            END IF

            IF (RPARAM (1) .LT. ZERO) THEN
               ADJUST = PCTRANGE
            ELSE
               ADJUST = SCALE
            END IF

            DX = RPARAM (1) * ADJUST

            CALL GEODIS (XA, XB, NX, DX, RPARAM (2), X, LUNOUT, IER)
  
            IF (QUERY) THEN
               DETAILS = DISNAME
               WRITE (DETAILS (11 : 42), '(A, 1P, E10.3, A, 0P, F6.2)')
     >            'dX(1) =', X (2) - X (1), '  Power =', RPARAM (2)
            END IF

            IF (IER .NE. 0) GO TO 920
            GO TO 800

  500    CONTINUE    ! GEODIS2: Modified geometric distribution (2-sided)

            IF (QUERYLOC) THEN
               DX = RANGE / TWOTON
               WRITE (LUNCRT, 1002)
     >            ' First increment? (dX > 0 is absolute;',
     >            ' dX < 0 is relative)'
               CALL READR (LUNCRT,
     >           'dX = -r means r% of range; default = 0.5% of range: ',
     >            LUNKBD, DX, DEFAULT,QUIT)
               IF (DX .EQ. ZERO .OR. DX .GE. RANGE * HALF) GO TO 500

               WRITE (LUNCRT, 1020)
               POWER = ZERO
               CALL READR (LUNCRT, 'First clustering exponent?  ',
     >            LUNKBD, POWER, DEFAULT, QUIT)

               POWER2 = ZERO
               CALL READR (LUNCRT, 'Second clustering exponent? ',
     >            LUNKBD, POWER2, DEFAULT, QUIT)

               RPARAM (1) = DX
               RPARAM (2) = POWER
               RPARAM (3) = POWER2
            END IF

            IF (RPARAM (1) .LT. ZERO) THEN
               ADJUST = PCTRANGE
            ELSE
               ADJUST = SCALE
            END IF

            DX = RPARAM (1) * ADJUST

            CALL GEODIS2 (XA, XB, NX, DX, RPARAM (2), RPARAM (3),
     >                    X, LUNOUT, IER)

            IF (QUERY) THEN
               DETAILS = DISNAME
               WRITE (DETAILS (11 : 46),
     >            '(A, 1P, E10.3, A, 0P, 2F6.2)')
     >            'dX(1) =', X (2) - X (1), '  Pwrs:', RPARAM (2),
     >            RPARAM (3)
            END IF

            IF (IER .NE. 0) GO TO 920
            GO TO 800

  520    CONTINUE    ! HTDIS4: Vinokur's TANH-based method (end slopes or
                     !         iterated for precise end dXs).

            DSINPUT = METHOD .EQ. 12    ! Else METHOD = 11

            IF (QUERYLOC) THEN
               IF (.NOT. DSINPUT) THEN  ! Specify end slopes
                  DX = ZERO
                  CALL READR (LUNCRT,
     >               'Stretching function slope at lower end? ',
     >               LUNKBD, DX, DEFAULT, QUIT)
                  IF (DX .EQ. ZERO) GO TO 520

                  DX2 = ZERO
                  CALL READR (LUNCRT,
     >               'Stretching function slope at upper end? ',
     >               LUNKBD, DX2, DEFAULT, QUIT)
                  IF (DX2 .EQ. ZERO) GO TO 520

               ELSE                     ! Specify end dXs
                  DX = RANGE / TWOTON
                  WRITE (LUNCRT, 1002)
     >               ' First increment? (dX > 0 is absolute;',
     >               ' dX < 0 is relative)'
                  CALL READR (LUNCRT,
     >           'dX = -r means r% of range; default = 0.5% of range: ',
     >               LUNKBD, DX, DEFAULT,QUIT)
                  IF (DX .EQ. ZERO .OR. DX .GE. RANGE * HALF) GO TO 520

                  DX2 = DX
                  CALL READR (LUNCRT,
     >               'Last increment?   Default = same as first: ',
     >               LUNKBD, DX2, DEFAULT, QUIT)
                  IF (DX2 .EQ. ZERO .OR. DX2 .GE. RANGE * HALF) GOTO 520

               END IF

               RPARAM (1) = DX
               RPARAM (2) = DX2
            END IF

            IF (RPARAM (1) .LT. ZERO) THEN
               ADJUST = PCTRANGE
            ELSE
               ADJUST = SCALE
            END IF

            DX = RPARAM (1) * ADJUST

            IF (RPARAM (2) .LT. ZERO) THEN
               ADJUST = PCTRANGE
            ELSE
               ADJUST = SCALE
            END IF

            DX2 = RPARAM (2) * ADJUST

            CALL HTDIS4 (DSINPUT, XA, XB, DX, DX2, NX, X, LUNOUT, IER)

            IF (QUERY) THEN
               DETAILS = DISNAME
               WRITE (DETAILS (11 : 48), '(A, 1P, E10.3, A, E10.3)')
     >            'dX(1) =', X (2) - X (1), '  dX(N-1) =',
     >            X (NX) - X (NX-1)
            END IF

            IF (IER .NE. 0) GO TO 920
            GO TO 800

  530    CONTINUE    ! FOILGRID: Hybrid cosine/quadratic, suited to airfoils

            IF (QUERYLOC) THEN
               DX = 0.7
               CALL READR (LUNCRT,
     >           'Weight on cosine term?  ([0., 1.]; <CR> = 0.70) ',
     >            LUNKBD, DX, DEFAULT,QUIT)
               IF (DX .LT. ZERO .OR. DX .GT. ONE) GO TO 530

               RPARAM (1) = DX
            END IF

            CALL FOILGRID (NX, XA, XB, RPARAM, X)

            IF (QUERY) THEN
               DETAILS = DISNAME
               WRITE (DETAILS (9 : 32), '(A, F4.2)')
     >            ': Cosine weight = ', RPARAM (1)
            END IF

            GO TO 800

  540    CONTINUE    ! FOILGRD: Linear/quadratic/sine/cosine, suited to airfoils

            IF (QUERYLOC) THEN
               WRITE (LUNCRT, '(A)', ADVANCE='NO')
     >         ' Linear/quad/sin/cos weights? [E.g., .04, .0, .3, .66] '
               READ (LUNKBD, *)
     >            RPARAM (1), RPARAM (2), RPARAM (3), RPARAM (4)
            END IF

            CALL FOILGRD (NX, XA, XB, RPARAM (1), RPARAM (2),
     >                    RPARAM (3), RPARAM (4), X)

            IF (QUERY) THEN
               DETAILS = DISNAME
               WRITE (DETAILS (9 : 45), '(A, 4F5.2)')
     >            ': L/Q/S/C wts. = ', RPARAM (1), RPARAM (2),
     >            RPARAM (3), RPARAM (4)
            END IF

            GO TO 800

  550    CONTINUE    ! BLGRID: Geometric/Vinokur suited to boundary layers

            IF (QUERYLOC) THEN
  555          M = 20
               CALL READI (LUNCRT, 'NBLAYER? ', LUNKBD, M, DEFAULT,
     >                     QUIT)
               IF (M .LT. 0 .OR. M .GE. NX) GO TO 555

               IPARAM (1) = M

               WRITE (LUNCRT, '(A)', ADVANCE='NO')
     >         ' D1, D2, and RBLAYER? [E.g., 1.E-6, 0.1, 1.1] '
               READ (LUNKBD, *)
     >            RPARAM (1), RPARAM (2), RPARAM (3)
            END IF

            X (1)  = XA ! Assumed to be in place by BLGRID
            X (NX) = XB

            CALL BLGRID (NX, RPARAM (1), RPARAM (2), IPARAM (1),
     >                   RPARAM (3), X, LUNCRT, IER)

            IF (IER .NE. 0) GO TO 920

            IF (QUERY) THEN
               DETAILS = DISNAME
               WRITE (DETAILS (9 :), '(A, I2, A, F5.2)')
     >            ': NBL =', IPARAM (1), ', RBL =', RPARAM (3)
            END IF

            GO TO 800

  560    CONTINUE    ! STRETCH: Exponential-type (given smallest interval)

            IF (QUERYLOC) THEN

  563          DX = RANGE / TWOTON     ! DX << range / (N - 1)

               WRITE (LUNCRT, 1002)
     >            ' Desired smallest interval? (dX > 0 is absolute;',
     >            ' dX < 0 is relative)'
               CALL READR (LUNCRT,
     >           'dX = -r means r% of range; default = 0.5% of range: ',
     >            LUNKBD, DX, DEFAULT,QUIT)
               IF (DX .EQ. ZERO .OR. DX .GE. RANGE / (NX - 1)) GO TO 563
               RPARAM (1) = DX

  566          BETA = 1.001E+0
               WRITE (LUNCRT, '(2A)')
     >            ' Suggestion for stretching parameter Beta: ',
     >            ' 1 + dX(1) / (XB - XA)'
               CALL READR (LUNCRT,
     >            'Enter starting guess for Beta > 1.0.  [1.001]: ',
     >            LUNKBD, BETA, DEFAULT, QUIT)
               IF (BETA .LE. ONE) GO TO 566
               RPARAM (2) = BETA
            END IF

            IF (RPARAM (1) .LT. ZERO) THEN
               ADJUST = PCTRANGE
            ELSE
               ADJUST = SCALE
            END IF

            DX = (RPARAM (1) * ADJUST) / RANGE

            CALL STRETCH (NX, DX, RPARAM (2), X)  ! Assumes [0, 1] interval

            BETA = RPARAM (2)           ! Output

            DO I = 1, NX
               X (I) = XA + RANGE * X (I)
            END DO

            IF (QUERY) THEN
               DETAILS = DISNAME
               WRITE (DETAILS (11 : 44), '(A, 1P, E10.3, A, 0P, F9.6)')
     >            'dX(1) =', X (2) - X (1), '  Beta =', BETA
            END IF
            GO TO 800

  570    CONTINUE    ! SHOCKGRID: 2-sided + 1-sided Vinokur for b.layer + shock

            IF (QUERYLOC) THEN
  575          M = 3
               CALL READI (LUNCRT, 'NMARGIN? ', LUNKBD, M, DEFAULT,
     >                     QUIT)
               IF (M .LT. 0 .OR. M .GE. NX) GO TO 575

               IPARAM (1) = M

               WRITE (LUNCRT, '(A)', ADVANCE='NO')
     >         ' D1, DE, and SEDGE? [E.g., 1.E-6, 2.5E-4, 0.98] '
               READ (LUNKBD, *)
     >            RPARAM (1), RPARAM (2), RPARAM (3)
            END IF

            X (1)  = XA ! Assumed to be in place by SHOCKGRID
            X (NX) = XB

            CALL SHOCKGRID (NX, RPARAM (1), RPARAM (2), IPARAM (1),
     >                      RPARAM (3), X, LUNCRT, IER)

            IF (IER .NE. 0) GO TO 920

            IF (QUERY) THEN
               DETAILS = DISNAME
               WRITE (DETAILS (9 :), '(A, I2, A, F5.2)')
     >            ': NMARGIN =', IPARAM (1), ', SEDGE =', RPARAM (3)
            END IF

            GO TO 800

  580    CONTINUE    ! EXPDIS5: Exponential-type (given smallest interval)

            IF (QUERYLOC) THEN
               MODE = 1
               ANSWER = 'L'

               CALL READC (LUNCRT,
     >            'Input dX at low end? high end? [L/H; <CR> = Low] ',
     >            LUNKBD, ANSWER, DEFAULT, QUIT)
               IF (ANSWER .EQ. 'L') MODE = 1
               IF (ANSWER .EQ. 'H') MODE = 2
               IPARAM (1) = MODE

  585          DX = RANGE / TWOTON     ! DX << range / (N - 1)

               WRITE (LUNCRT, 1002)
     >            ' Desired end interval? (dX > 0 is absolute;',
     >            ' dX < 0 is relative)'
               CALL READR (LUNCRT,
     >           'dX = -r means r% of range; default = 0.5% of range: ',
     >            LUNKBD, DX, DEFAULT,QUIT)
               IF (DX .EQ. ZERO .OR. DX .GE. RANGE) GO TO 585
               RPARAM (1) = DX
            END IF

            IF (RPARAM (1) .LT. ZERO) THEN
               ADJUST = PCTRANGE
            ELSE
               ADJUST = SCALE
            END IF

            DX = RPARAM (1) * ADJUST

            CALL EXPDIS5 (IPARAM, XA, XB, DX, NX, X, LUNOUT)

            IF (QUERY) THEN
               DETAILS = DISNAME
               IF (MODE == 1) THEN
                  WRITE (DETAILS (11 : 27), '(A, 1P, E10.3)')
     >               'dX(1) =', X (2) - X (1)
               ELSE
                  WRITE (DETAILS (11 : 27), '(A, 1P, E10.3)')
     >               'dX(N) =', X (NX) - X (NX - 1)
               END IF
            END IF

            GO TO 800


  800    CONTINUE

C     End of case statement over available methods.

      IER = 0
      GO TO 999


C     Error handling:

  900 WRITE (LUNCRT, 1005) SUBNAME, ': Invalid METHOD: ', METHOD
      IER = 1
      GO TO 980

  910 WRITE (LUNCRT, 1004)
     >   SUBNAME, ': Not enough work-space for ', DISNAME,
     >   'Needed: ', NEEDED, '  Provided: ', NWK
      IER = 2
      GO TO 980

  920 WRITE (LUNCRT, 1004)
     >   SUBNAME, ': Bad return from ', DISNAME,
     >   'Error code returned: ', IER
      IER = 3
      GO TO 980


  980 ANSWER = 'R'
      CALL READC (LUNCRT,
     >   'Retry?  Proceed anyway?  Or Quit?  [R/P/Q; <CR> = Retry] ',
     >   LUNKBD, ANSWER, DEFAULT, QUIT)
      QUERYLOC = .TRUE.
      IF (ANSWER .EQ. 'R') GO TO 200
      IF (ANSWER .EQ. 'P') GO TO 800     ! Else drop through (Quit).


  990 IER = -1          ! Quit from menu - may be normal (DISTRIBUTE) or
      GO TO 999         ! abnormal (other applications).

  999 RETURN

C     Formats:

 1002 FORMAT (A, A)
 1003 FORMAT (A, I2, A)
 1004 FORMAT (1X, A, A, A, /, 1X, (A, I4))
 1005 FORMAT (A, A, I2)
 1010 FORMAT (/, ' Exponents in the range (0., 1.) give greater',
     >        ' bunching at the end(s) than the', /,
     >        ' simple 1. gives, and no additional bunching.', /,
     >        ' Powers > 1. also give bunching at the other end',
     >        ' or mid-pt (probably undesired).')
 1020 FORMAT (/, ' True geometric variation corresponds to a zero',
     >        ' clustering exponent.', /,
     >        ' > 0 amplifies bunching; < 0 reduces it.',
     >        ' Rule of thumb:  -1.0 through +0.5.')

      END SUBROUTINE DISTRIB
