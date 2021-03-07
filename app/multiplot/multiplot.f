C*******************************************************************************
C
      PROGRAM MULTIPLOT
C
C     Description:
C     ------------
C
C     MULTIPLOT produces PostScript plots of scalar surface data from a
C     pair of "XYZ" and "Q" files containing one or more surface datasets
C     in PLOT3D /mgrid format.  Plots from a second pair may be overlaid.
C     Structured target Cp distributions may also be overlaid, and some
C     provision is made for displaying wind tunnel data (no interpolation).
C
C     The plots are produced at the indicated cuts (constant X and/or Z),
C     one plot page per cut.  (Use PLOT3D or equivalent to obtain plots along
C     grid lines.)
C
C     Where relevant, total and sectional force coefficients are calculated
C     and displayed on each plot page.
C
C     This version has an option to save (X,Y,Cp) and (Z,Y,Cp) data points
C     from the cuts for plotting by some other means.  The output format
C     matches that of the optional input wind tunnel data, but with Y also.
C
C     (X,Y,Z) Convention:  That of the Jameson flow codes (right-handed):
C     -------------------
C
C        X is streamwise  (+ve pointing downstream)
C        Y is vertical    (+ve pointing up)
C        Z is spanwise    (+ve pointing from root to tip of left wing)
C
C     Input data format:  PLOT3D /mgrid /whole /[un]formatted grid and Q files
C     ------------------
C
C        NGRID                                  NGRID
C        IDIM1 JDIM1 1 [IDIM2 JDIM2 1] ...      IDIM1 JDIM1 1 ...
C        All Xs for surface # 1                 FSMACH  ALPHA  RE  TIME
C        All Ys  "     "      "                 All Q1s for surface # 1
C        All Zs  "     "      "                 All Q2s  "     "      "
C       [All Xs for surface # 2]                 :   :   :   :   :   :
C         :   :   :   :   :   :                  :   :   :   :   :   :
C         :   :   :   :   :   :                [FSMACH  ALPHA  RE  TIME]
C         :   :   :   :   :   :                  :   :   :   :   :   :
C
C     Target Cp data format (ASCII):
C     ------------------------------
C
C        IDIM JDIM
C        x  y  z  Cp
C        :  :  :  :
C        :  :  :  :
C
C     Wind tunnel Cp data format:  ASCII; separate files for each cut direction
C     ---------------------------
C
C        Title                       and/or    Title
C        NROWS  (of pressure taps)             NROWS  (of pressure taps)
C        7      (say)                          7
C        NPTS  X                               NPTS  Z
C        34    88.5                            34    88.9
C        Z     CP                              X     CP
C        0.000 -0.1234                         123.0 -0.1234
C        1.234 -0.1234                         123.5 -0.1245
C        :     :                               :     :
C        :     :                               :     :
C        NPTS  X                               NPTS  Z
C        m     x                               :     :
C        Z     CP
C        z     p
C        :     :
C
C     Control inputs:
C     ---------------
C
C        multiplot.inp      See a sample.
C
C     Outputs:
C     --------
C
C        xcuts.ps           PostScript plots of cuts at specified Xs
C        zcuts.ps             "    "    "    "    "    "    "     Zs
C
C        Cp may be suppressed if the Q file name is 'none' or 'NONE'.
C        See also the SHOWX and SHOWF inputs associated with each block.
C
C        xcuts.dat          (Z,Y,Cp) points of cuts at specified Xs
C        zcuts.dat          (X,Y,Cp)    "    "    "    "    "    Zs
C
C        The format is analogous to the above wind tunnel data format.
C        Cp may be suppressed if the Q file name is 'none' or 'NONE'.
C        See also the SHOWX and SHOWF inputs associated with each block.
C
C     Environment:          Fortran 90
C     ------------
C
C     History:
C     --------
C
C     11/22/96  D.Saunders  Initial design around SYN87 plot utilities,
C                           for SYN87-MB purposes (FORTRAN 77).
C     12/03/97      "       Provided for target Cps.
C     03/28/98      "       Large surface patches can produce skewed (u,v)s.
C                           Packing the data allows handling more patches.
C                           Went to Fortran 90 (barely).
C     03/31/98      "       Abandoned the contouring approach to avoid skewed
C                           (u,v) problems.  CUTSURF from Scott Thomas adapts
C                           his unstructured surface grid techniques.
C     04/10/98      "       CUTSURF & CUTORDER were revised to assemble the
C                           line segments the way the original CONT2D did.
C     04/16/98      "       Provided for wind tunnel data.
C     01/22/99      "       Enabled suppressing the Q-file function for
C                           specified blocks but retaining the corresponding
C                           geometry in the plot; switched from "blank"ing to
C                           "display"ing specified blocks; provided for
C                           reversing the sign of coefs. computed for each
C                           block; added a legend line for the primary data
C                           (where the plot title sufficed originally);
C                           enabled shifting the default legend position.
C     04/08/99      "       Entering 'none' or 'NONE' for the *.q file name
C                           enables slicing just the (x,y,z)s (kludge).
C     05/13/99      "       Option to save sliced data for plotting elsewhere.
C
C     Author:  David Saunders, Sterling Software/Raytheon/NASA Ames, CA
C
C*******************************************************************************

      IMPLICIT NONE

C     Constants:

      INTEGER, PARAMETER ::
     >   IFUNCP = 114,   ! PLOT3D's function code for Cp
     >   LUNCRT = 6,
     >   LUNCTL = 1,     ! Control inputs
     >   LUNDAT = 2,     ! XYZ and Q file(s) and target Cps and wind tunnel data
     >   LUNKBD = 5,
     >   LUNPLT = 3,     ! Output PostScript plot file(s)
     >   LUNSAV = 4,     ! Cut data saved for plotting elsewhere if specified
     >   MAXBLK = 30,    ! Max. # blocks (surfaces) per set
     >   MAXCPS = 50000, ! Max. # target Cp points
     >   MAXCUT = 201,   ! Max. # plot cuts in either direction
     >   MAXROW = 30,    ! Max. # rows of WT Cp data in either direction
     >   MAXPTS = 200000,! Max. # pts. total in the current or overlaid datasets
     >   MAXSEG = 2000,  ! Max. # 2-pt. line segments handled per plot cut
     >   MAXTAP = 200    ! Max. # pressure taps per row of WT Cp data

      REAL, PARAMETER ::
     >   DEFAULT= 999.,  ! Inputs read as DEFAULT are reset to sensible values
     >   EPSWT  = 0.001, ! Fraction of wind tunnel data range for cut-match tol.
     >   ONE    = 1.,
     >   ZERO   = 0.

      CHARACTER, PARAMETER ::
     >   NONEL * 4 = 'none',
     >   NONEU * 4 = 'NONE'

C     Variables:

      INTEGER
     >   MF(MAXBLK),  MQ(MAXBLK),  MX(MAXBLK),
     >   MFO(MAXBLK), MQO(MAXBLK), MXO(MAXBLK),
     >   NI(MAXBLK),  NIO(MAXBLK), NIQ(MAXBLK),
     >   NJ(MAXBLK),  NJO(MAXBLK), NJQ(MAXBLK),
     >   NTAPS(MAXROW),
     >   I, I1, I2, IFUNPLOT3D, IPAGE, J, JROW, LUNSAVE, MAXI, MAXJ,
     >   N, NBLK, NBLKO, NBLKQ, NFIRST, NFIRSTO, NIJ,
     >   NITARG, NJTARG, NLAST, NLASTO, NPTS, NROWS, NXCUTS, NZCUTS
      REAL
     >   FUN(MAXPTS),  Q(MAXPTS*5),  XYZ(MAXPTS*3),
     >   FUNO(MAXPTS), QO(MAXPTS*5), XYZO(MAXPTS*3),
     >   SIGNBLK(MAXBLK), SIGNBLKO(MAXBLK),
     >   CPTARG(MAXCPS),   CPWT(MAXTAP,MAXROW),
     >   XYZTARG(MAXCPS*3), XWT(MAXTAP,MAXROW), ZWT(MAXROW),
     >   SEG(4,2*MAXSEG), XCUTS(MAXCUT), ZCUTS(MAXCUT),
     >   FSMACH,  ALPHA,  RE,  REY,   TIME,  CREF,  SREF,
     >   FSMACHO, ALPHAO, REO, REYO,  TIMEO, CREFO, SREFO,
     >   CL,  CD,  CY,  CLBLK, CDBLK, CYBLK,
     >   CLO, CDO, CYO, CPC,   COSA,  SINA, TOLWT, WTMIN, WTMAX,
     >   XLEFT, XRIGHT, XSTEP, XSHIFTLEG, FBOT, FTOP, FSTEP,
     >   ZLEFT, ZRIGHT, ZSTEP, YSHIFTLEG, YMIN, YSCALEVX, YSCALEVZ
      LOGICAL
     >   SHOWF(MAXBLK), SHOWX(MAXBLK), SHOWFO(MAXBLK), SHOWXO(MAXBLK),
     >   FINDYMIN, FIRST, FORMATTED, LAST, LOCALCFS, OVERLAY, LINES,
     >   PLTXCUTS, PLTZCUTS, SAVECUTS, SECCFS, SECCFSO, TARGETS,
     >   TARGSYM, TUNNELX, TUNNELZ, XYZONLY, XAXES, YAXES, ZAXES
      CHARACTER
     >   FILENAME*48, FLABEL*2, FORM*11, LEGENDC*48, LEGENDO*48,
     >   LEGENDW*48, TITLE*80, WTXFILE*48, WTZFILE*48

C     Procedures:

      REAL
     >   COSD, SIND

C     Execution:

C     ********************** Read control inputs ***********************

C     Control inputs are read as they are needed, as opposed to reading
C     them all ahead of time.  This keeps down the count of variables
C     such as FILENAME and LUNDAT.  [Potentially a bad design decision.]

      OPEN (UNIT=LUNCTL, FILE='multiplot.inp', STATUS='OLD')

      READ (LUNCTL, '(A)') TITLE
      READ (LUNCTL, *) FORMATTED
      IF (FORMATTED) THEN
         FORM = 'FORMATTED'
      ELSE
         FORM = 'UNFORMATTED'
      END IF

      READ (LUNCTL, *) FILENAME ! Current multiblock grid file

      OPEN (UNIT=LUNDAT, FILE=FILENAME, STATUS='OLD', FORM=FORM)

      IF (FORMATTED) THEN
         READ (LUNDAT, *) NBLK
      ELSE
         READ (LUNDAT) NBLK
      END IF

      IF (NBLK > MAXBLK) THEN
         WRITE (LUNCRT, '(1X, A, I6)')
     >      'Too many blocks in XYZ file:', NBLK, 'Limit:', MAXBLK
         GO TO 999
      END IF

      IF (FORMATTED) THEN
         READ (LUNDAT, *) (NI(N), NJ(N), I, N = 1, NBLK)
      ELSE
         READ (LUNDAT)    (NI(N), NJ(N), I, N = 1, NBLK)
      END IF

      NPTS = 0
      DO N = 1, NBLK
         NPTS  = NI(N) * NJ(N) + NPTS
      END DO

      WRITE (LUNCRT, '(/, A)') ' Surface Block Dimensions'
      WRITE (LUNCRT, '(I3, 2I5)') (N, NI(N), NJ(N), N = 1, NBLK)
      WRITE (LUNCRT, '(/, A, I8)') ' Total Number of Points:', NPTS

      IF (NPTS > MAXPTS) THEN
         WRITE (LUNCRT, '(A, I8)') ' Too many points. Limit:', MAXPTS
         GO TO 999
      END IF

      MF(1) = 1
      MX(1) = 1
      MQ(1) = 1

      DO N = 1, NBLK - 1
         NIJ = NI(N) * NJ(N)
         MF(N+1) = MF(N) + NIJ     ! Index of 1st point in each F block
         MX(N+1) = MX(N) + NIJ * 3 ! ... and each XYZ block
         MQ(N+1) = MQ(N) + NIJ * 5 ! ... and each Q block
      END DO

      DO N = 1, NBLK
         I1 = MX(N)
         I2 = I1 + 3 * NI(N) * NJ(N) - 1
         IF (FORMATTED) THEN
            READ (LUNDAT, *) XYZ(I1:I2)
         ELSE
            READ (LUNDAT)    XYZ(I1:I2)
         END IF
      END DO
            
      CLOSE (LUNDAT)

      READ (LUNCTL, *) FILENAME ! Current multiblock Q file

      XYZONLY = FILENAME(1:4) == NONEL .OR.
     >          FILENAME(1:4) == NONEU

      IF (.NOT. XYZONLY) THEN

         OPEN (UNIT=LUNDAT, FILE=FILENAME, STATUS='OLD', FORM=FORM)

         IF (FORMATTED) THEN
            READ (LUNDAT, *) NBLKQ
         ELSE
            READ (LUNDAT) NBLKQ
         END IF

         IF (NBLK /= NBLKQ) THEN
            WRITE (LUNCRT, '(1X, A, 2I6)')
     >         'Mismatched XYZ & Q file block counts:', NBLK, NBLKQ
            GO TO 999
         END IF

         IF (FORMATTED) THEN
            READ (LUNDAT, *) (NIQ(N), NJQ(N), I, N = 1, NBLK)
         ELSE
            READ (LUNDAT)    (NIQ(N), NJQ(N), I, N = 1, NBLK)
         END IF

         DO N = 1, NBLK

            IF (NI(N) /= NIQ(N) .OR. NJ(N) /= NJQ(N)) THEN
               WRITE (LUNCRT, '(1X, A, I6)')
     >            'Mismatched XYZ & Q block sizes.  Block:', N,
     >            'NIX:', NI(N),  'NJX:', NJ(N),
     >            'NIQ:', NIQ(N), 'NJQ:', NJQ(N)
               GO TO 999
            END IF

            I1 = MQ(N)
            I2 = I1 + 5 * NI(N) * NJ(N) - 1

            IF (FORMATTED) THEN
               READ (LUNDAT, *) FSMACH, ALPHA, RE, TIME
               READ (LUNDAT, *) Q(I1:I2)
            ELSE
               READ (LUNDAT) FSMACH, ALPHA, RE, TIME
               READ (LUNDAT) Q(I1:I2)
            END IF

         END DO
            
         CLOSE (LUNDAT)

      ELSE
         FSMACH = ZERO
         ALPHA  = ZERO
         RE     = ZERO
         TIME   = ZERO
      END IF

      READ (LUNCTL, *) LEGENDC       ! Legend for current solution
      READ (LUNCTL, *) CREF, SREF    ! Reference chord & area
      READ (LUNCTL, *) REY           ! Reynolds #; 999. = use data value
      READ (LUNCTL, *) SECCFS        ! Section coefs.? (Z cuts only)
      READ (LUNCTL, *) SAVECUTS      ! Save cuts for plotting elsewhere?

      IF (REY /= DEFAULT) RE = REY   ! Control input may be 0. to mean inviscid

      DO N = 1, NBLK                 ! Show Q and/or geometry for which blocks?
         READ (LUNCTL, *) I, SHOWX(N), SHOWF(N), SIGNBLK(N)  ! & reverse normal?
         IF (XYZONLY) SHOWF(N) = .FALSE.
      END DO

      READ (LUNCTL, *) OVERLAY       ! Overlay a second solution?
      READ (LUNCTL, *) FILENAME      ! Overlay grid file name, used or not

      IF (OVERLAY) THEN

         OPEN (UNIT=LUNDAT, FILE=FILENAME, STATUS='OLD', FORM=FORM)

         IF (FORMATTED) THEN
            READ (LUNDAT, *) NBLKO
         ELSE
            READ (LUNDAT) NBLKO
         END IF

         IF (NBLKO > MAXBLK) THEN
            WRITE (LUNCRT, '(1X, A, I6)')
     >         'Too many blocks in overlay X file:', NBLKO,
     >         'Limit:', MAXBLK
            GO TO 999
         END IF

         IF (FORMATTED) THEN
            READ (LUNDAT, *) (NIO(N), NJO(N), I, N = 1, NBLKO)
         ELSE
            READ (LUNDAT)    (NIO(N), NJO(N), I, N = 1, NBLKO)
         END IF

         NPTS = 0
         DO N = 1, NBLKO
            NPTS = NIO(N) * NJO(N) + NPTS
         END DO

         WRITE (LUNCRT, '(/, A)') ' Overlay Block Dimensions'
         WRITE (LUNCRT, '(I3, 2I5)') (N, NIO(N), NJO(N), N = 1, NBLKO)
         WRITE (LUNCRT, '(/, A, I8)') ' Total Number of Points:', NPTS

         IF (NPTS > MAXPTS) THEN
            WRITE (LUNCRT, '(A, I8)') ' Too many points. Limit:', MAXPTS
            GO TO 999
         END IF

         MFO(1) = 1
         MXO(1) = 1
         MQO(1) = 1

         DO N = 1, NBLKO - 1
            NIJ = NIO(N) * NJO(N)
            MFO(N+1) = MFO(N) + NIJ     ! 1st point in overlay F block
            MXO(N+1) = MXO(N) + NIJ * 3 ! ... and each overlay XYZ block
            MQO(N+1) = MQO(N) + NIJ * 5 ! ... and each overlay Q block
         END DO

         DO N = 1, NBLKO
            I1 = MXO(N)
            I2 = I1 + 3 * NIO(N) * NJO(N) - 1
            IF (FORMATTED) THEN
               READ (LUNDAT, *) XYZO(I1:I2)
            ELSE
               READ (LUNDAT)    XYZO(I1:I2)
            END IF
         END DO
            
         CLOSE (LUNDAT)

      ELSE
         NBLKO = 0
      END IF

      READ (LUNCTL, *) FILENAME ! Overlay Q file

      IF (XYZONLY) THEN
         FSMACHO = ZERO
         ALPHAO  = ZERO
         REO     = ZERO
         TIMEO   = ZERO
      ELSE

         IF (OVERLAY) THEN

            OPEN (UNIT=LUNDAT, FILE=FILENAME, STATUS='OLD', FORM=FORM)

            IF (FORMATTED) THEN
               READ (LUNDAT, *) NBLKQ
            ELSE
               READ (LUNDAT) NBLKQ
            END IF

            IF (NBLKO /= NBLKQ) THEN
               WRITE (LUNCRT, '(1X, A, 2I6)')
     >          'Mismatched XYZ & Q overlay block counts:', NBLKO, NBLKQ
               GO TO 999
            END IF

            IF (FORMATTED) THEN
               READ (LUNDAT, *) (NIQ(N), NJQ(N), I, N = 1, NBLKO)
            ELSE
               READ (LUNDAT)    (NIQ(N), NJQ(N), I, N = 1, NBLKO)
            END IF

            DO N = 1, NBLKO

               IF (NIO(N) /= NIQ(N) .OR. NJO(N) /= NJQ(N)) THEN
                  WRITE (LUNCRT, '(1X, A, I6)')
     >             'Mismatched overlay XYZ & Q block sizes.  Block:', N,
     >             'NIO:', NIO(N), 'NJO:', NJO(N),
     >             'NIQ:', NIQ(N), 'NJQ:', NJQ(N)
                  GO TO 999
               END IF

               I1 = MQO(N)
               I2 = I1 + 5 * NIO(N) * NJO(N) - 1

               IF (FORMATTED) THEN
                  READ (LUNDAT, *) FSMACHO, ALPHAO, REO, TIMEO
                  READ (LUNDAT, *) QO(I1:I2)
               ELSE
                  READ (LUNDAT) FSMACHO, ALPHAO, REO, TIMEO
                  READ (LUNDAT) QO(I1:I2)
               END IF

            END DO
            
            CLOSE (LUNDAT)

         END IF

      END IF

      READ (LUNCTL, *) LEGENDO       ! Legend for overlaid solution
      READ (LUNCTL, *) CREFO, SREFO  ! Reference chord & area for overlaid soln.
      READ (LUNCTL, *) REYO          ! Reynolds #; 999. = use data value
      READ (LUNCTL, *) SECCFSO       ! Overlaid section coefs.? (Z cuts only)

      IF (REYO /= DEFAULT) REO = REYO! Control may be 0. to mean inviscid

C     Overlaid surface suppression is handled last in case it is irrelevant.

      READ (LUNCTL, *) TARGETS, TARGSYM ! Target Cps?  Symbols?
      READ (LUNCTL, *) FILENAME      ! Whether it'll be used or not
      READ (LUNCTL, *) TUNNELX       ! Wind tunnel Cps for X cuts?
      READ (LUNCTL, *) WTXFILE       ! Whether it'll be used or not
      READ (LUNCTL, *) TUNNELZ       ! Wind tunnel Cps for Z cuts?
      READ (LUNCTL, *) WTZFILE       ! Whether it'll be used or not
      READ (LUNCTL, *) LEGENDW       ! Legend for wind tunnel data
      READ (LUNCTL, *) IFUNPLOT3D    ! Scalar plot quantity; 114 = Cp, etc.

      TARGETS = TARGETS .AND. IFUNPLOT3D == IFUNCP
      TUNNELX = TUNNELX .AND. IFUNPLOT3D == IFUNCP
      TUNNELZ = TUNNELZ .AND. IFUNPLOT3D == IFUNCP

      READ (LUNCTL, *) PLTXCUTS      ! Allows retention of suppressed cuts
      READ (LUNCTL, *) PLTZCUTS
      READ (LUNCTL, *) NXCUTS

      IF (NXCUTS > MAXCUT) THEN
         WRITE (LUNCRT, '(1X, A, I6)') 'Too many X cuts:', NXCUTS,
     >      'Limit:', MAXCUT
         GO TO 999 ! More trouble than it's worth skipping over the excess
      END IF

      DO N = 1, NXCUTS
         READ (LUNCTL, *) I, XCUTS(N)
      END DO

      READ (LUNCTL, *) NZCUTS

      IF (NZCUTS > MAXCUT) THEN
         WRITE (LUNCRT, '(1X, A, I6)') 'Too many Z cuts:', NZCUTS,
     >      'Limit:', MAXCUT
         GO TO 999
      END IF

      DO N = 1, NZCUTS
         READ (LUNCTL, *) I, ZCUTS(N)
      END DO

      READ (LUNCTL, *) XAXES, YAXES,  ZAXES
      READ (LUNCTL, *) XLEFT, XRIGHT, XSTEP
      READ (LUNCTL, *) FBOT,  FTOP,   FSTEP
      READ (LUNCTL, *) ZLEFT, ZRIGHT, ZSTEP
      READ (LUNCTL, *) YMIN
      READ (LUNCTL, *) YSCALEVX, YSCALEVZ
      READ (LUNCTL, *) XSHIFTLEG, YSHIFTLEG
      
      IF (OVERLAY) THEN ! Which part of the overlaid data to show?
         DO N = 1, NBLKO
            READ (LUNCTL, *) I, SHOWXO(N), SHOWFO(N), SIGNBLKO(N)
            IF (XYZONLY) SHOWFO(N) = .FALSE.
         END DO
      END IF

      CLOSE (LUNCTL)

C     ********************* End of control inputs **********************


      IF (TARGETS) THEN

         OPEN (UNIT=LUNDAT, FILE=FILENAME, STATUS='OLD')

         READ (LUNDAT, *) NITARG, NJTARG
         NIJ = NITARG * NJTARG

         IF (NIJ > MAXCPS) THEN
            WRITE (LUNCRT, '(1X, A, 2I7)')
     >         'Too many target Cps: ', NITARG, NJTARG,
     >         'Total and limit:     ', NIJ, MAXCPS
            GO TO 999
         END IF

         READ (LUNDAT, *) (XYZTARG(I), XYZTARG(I+NIJ),
     >                     XYZTARG(I+2*NIJ), CPTARG(I), I = 1, NIJ)
         CLOSE (LUNDAT)

      END IF


C     ************* Force integrations for current surfaces ************

      FINDYMIN = YMIN == DEFAULT
      IF (FINDYMIN) YMIN = HUGE (YMIN)

      COSA = COSD (ALPHA)
      SINA = SIND (ALPHA)
      CL = ZERO
      CD = ZERO
      CY = ZERO

      DO N = 1, NBLK

         MAXI = NI(N)
         MAXJ = NJ(N)
         I1   = MQ(N)
         I2   = MF(N)
         NIJ  = MAXI * MAXJ

         IF (XYZONLY) THEN

            FUN(I2:I2 + NIJ - 1) = ZERO ! PLCUT does range calculations

         ELSE

C           Evaluate Cps on this surface:

            CALL FPLOT3D (IFUNCP, FSMACH, ALPHA, RE, MAXI, MAXJ,
     >                    1, 1, MAXI, 1, MAXJ, 1, 1, Q(I1), FUN(I2))

C           Forces on this surface:

            I1 = MX(N)

            CALL FORCES (COSA, SINA, SREF, MAXI, MAXJ, 1, MAXI, 1, MAXJ,
     >                   XYZ(I1), FUN(I2), CLBLK, CDBLK, CYBLK)

            CL = CL + CLBLK * SIGNBLK(N)
            CD = CD + CDBLK * SIGNBLK(N)
            CY = CY + CYBLK * SIGNBLK(N)

         END IF

         IF (FINDYMIN) THEN
            I1  = MX(N)
            DO I = I1 + NIJ, I1 + 2 * NIJ - 1
               YMIN = MIN (XYZ(I), YMIN)
            END DO
         END IF

      END DO

      IF (FINDYMIN) WRITE (LUNCRT, '(/, A, F10.2)')
     >   ' Minimum Y found: ', YMIN


C     ********************* Prepare data for plots *********************

C     Evaluate the specified scalar function:

      DO N = 1, NBLK

         IF (SHOWF(N)) THEN

            IF (IFUNPLOT3D /= IFUNCP) THEN

               MAXI = NI(N)
               MAXJ = NJ(N)
               I1   = MQ(N)
               I2   = MF(N)

               CALL FPLOT3D (IFUNPLOT3D, FSMACH, ALPHA, RE, MAXI, MAXJ,
     >                       1, 1, MAXI, 1, MAXJ, 1, 1, Q(I1), FUN(I2))
            END IF

         END IF

      END DO


C     Likewise for an overlaid solution?

      IF (OVERLAY) THEN

         COSA = COSD (ALPHAO)
         SINA = SIND (ALPHAO)
         CLO = ZERO
         CDO = ZERO
         CYO = ZERO

         DO N = 1, NBLKO

            MAXI = NIO(N)
            MAXJ = NJO(N)
            I1   = MQO(N)
            I2   = MFO(N)
            NIJ  = MAXI * MAXJ

            IF (XYZONLY) THEN

               FUNO(I2:I2 + NIJ - 1) = ZERO ! PLCUT does range calculations

            ELSE

C              Evaluate Cps on this surface:

               CALL FPLOT3D (IFUNCP, FSMACHO, ALPHAO, REO, MAXI, MAXJ,
     >                       1, 1, MAXI, 1, MAXJ, 1, 1, QO(I1),FUNO(I2))

C              Forces on this surface:

               I1 = MXO(N)

               CALL FORCES (COSA, SINA, SREFO, MAXI, MAXJ, 1, MAXI, 1,
     >                      MAXJ, XYZO(I1), FUNO(I2), CLBLK,CDBLK,CYBLK)

               CLO = CLO + CLBLK * SIGNBLKO(N)
               CDO = CDO + CDBLK * SIGNBLKO(N)
               CYO = CYO + CYBLK * SIGNBLKO(N)

            END IF

         END DO

C        Prepare overlaid surfaces for plotting:

         DO N = 1, NBLKO

            IF (SHOWFO(N)) THEN

               IF (IFUNPLOT3D /= IFUNCP) THEN

                  MAXI = NIO(N)
                  MAXJ = NJO(N)
                  I1   = MQO(N)
                  I2   = MFO(N)

                  CALL FPLOT3D (IFUNPLOT3D, FSMACHO, ALPHAO, REO,
     >                          MAXI, MAXJ, 1, 1, MAXI, 1, MAXJ,
     >                          1, 1, QO(I1), FUNO(I2))
               END IF

            END IF

         END DO

      END IF


C     ********************* Plots at specified Xs **********************

      CPC = -999.
      FLABEL = 'Cp'

      IF (IFUNPLOT3D == IFUNCP) THEN
         IF (FSMACH > ZERO)
     >      CPC = ((SQRT ((5. + FSMACH**2) / 6.))**7 - 1.) /
     >            (0.7*FSMACH**2)
      ELSE
         FLABEL = '??'
      END IF

      DO N = 1, NBLK
         IF (SHOWX(N)) EXIT
      END DO

      NFIRST = N

      DO N = NBLK, 1, -1
         IF (SHOWX(N)) EXIT
      END DO

      NLAST = N

      IF (NBLKO > 0) THEN
         DO N = 1, NBLKO
            IF (SHOWXO(N)) EXIT
         END DO

         NFIRSTO = N

         DO N = NBLKO, 1, -1
            IF (SHOWXO(N)) EXIT
         END DO

         NLASTO = N
      ELSE
         NFIRSTO = 1
         NLASTO  = 0
      END IF


      IF (NXCUTS > 0 .AND. PLTXCUTS) THEN

         OPEN (UNIT=LUNPLT, FILE='xcuts.ps', STATUS='UNKNOWN')

         IF (SAVECUTS) THEN

            OPEN (UNIT=LUNSAV, FILE='xcuts.dat', STATUS='UNKNOWN')

            I = LEN_TRIM (TITLE)
            WRITE (LUNSAV, '(A)') TITLE(1:I), 'Number of X cuts'
            WRITE (LUNSAV, '(I3)') NXCUTS
            LUNSAVE = LUNSAV
         ELSE
            LUNSAVE = 0
         END IF

         IF (TUNNELX) THEN

            OPEN (UNIT=LUNDAT, FILE=WTXFILE, STATUS='OLD')

            READ (LUNDAT, *) ! Skip title
            READ (LUNDAT, *)
            READ (LUNDAT, *) NROWS

            IF (NROWS > MAXROW) THEN
               WRITE (LUNCRT, '(1X, A, 2I7)')
     >          'Too many wind tunnel tap rows (X cuts):', NROWS, MAXROW
               GO TO 999
            END IF

            WTMIN = HUGE (WTMIN) * 0.5
            WTMAX = -WTMIN

            DO J = 1, NROWS

               READ (LUNDAT, *)
               READ (LUNDAT, *) NTAPS(J), ZWT(J)

               IF (NTAPS(J) > MAXTAP) THEN
                  WRITE (LUNCRT, '(A, I3, /, A, 2I7)')
     >               ' Too many wind tunnel taps.  Row # (X cuts):', J,
     >               ' # taps and limit:', NTAPS(J), MAXTAP
                  GO TO 999
               END IF

               READ (LUNDAT, *)
               READ (LUNDAT, *) (XWT(I,J), CPWT(I,J), I = 1, NTAPS(J))

               DO I = 1, NTAPS(J)
                  WTMIN = MIN (WTMIN, XWT(I,J))
                  WTMAX = MAX (WTMAX, XWT(I,J))
               END DO

            END DO

            TOLWT = (WTMAX - WTMIN) * EPSWT

            CLOSE (LUNDAT)

         END IF

         LOCALCFS = .FALSE. ! Section coefs. aren't practical in general
         IPAGE = 0

         DO I = 1, NXCUTS

C           Current solution:

            LINES = .TRUE.

C           Check for wind tunnel Cps in order to end the plot if none:

            JROW = 0
            IF (TUNNELX) THEN

               DO J = 1, NROWS
                  IF (ABS (XCUTS(I) - ZWT(J)) < TOLWT) THEN
                     JROW = J
                     EXIT
                  END IF
               END DO

            END IF

            DO N = NFIRST, NLAST

               FIRST = N == NFIRST
               LAST  = N == NLAST .AND.
     >                 .NOT. (OVERLAY .OR. TARGETS .OR. JROW > 0)

               IF (SHOWX(N)) THEN

                  MAXI = NI(N)
                  MAXJ = NJ(N)
                  NIJ  = MAXI * MAXJ
                  I1   = MX(N)
                  I2   = MF(N)

                  CALL GRAPHF
     >              ('C', SHOWF(N), MAXI, MAXJ, MAXSEG,
     >               XYZ(I1+2*NIJ), XYZ(I1+NIJ), XYZ(I1),
     >               FUN(I2), SEG, 1, MAXI, 1, MAXJ,
     >               XCUTS(I), ZAXES, YAXES, ZLEFT, ZRIGHT, ZSTEP,
     >               FBOT, FTOP, FSTEP, YMIN, YSCALEVZ, 'Z',
     >               FLABEL, TITLE, LEGENDC, LEGENDO, LEGENDW,
     >               XSHIFTLEG, YSHIFTLEG, FSMACH, ALPHA, RE,
     >               CL, CD, CY, CPC, LOCALCFS, SIGNBLK(N), CREF,
     >               LINES, IPAGE, FIRST, LAST, LUNPLT, LUNSAVE)

               END IF

            END DO

C           Overlaid solution?

            DO N = NFIRSTO, NLASTO

               FIRST = N == NFIRSTO
               LAST  = N == NLASTO .AND. .NOT. (TARGETS .OR. JROW > 0)

               IF (SHOWXO(N)) THEN

                  MAXI = NIO(N)
                  MAXJ = NJO(N)
                  NIJ  = MAXI * MAXJ
                  I1   = MXO(N)
                  I2   = MFO(N)

                  CALL GRAPHF
     >              ('O', SHOWFO(N), MAXI, MAXJ, MAXSEG,
     >               XYZO(I1+2*NIJ), XYZO(I1+NIJ), XYZO(I1),
     >               FUNO(I2), SEG, 1, MAXI, 1, MAXJ,
     >               XCUTS(I), ZAXES, YAXES, ZLEFT, ZRIGHT, ZSTEP,
     >               FBOT, FTOP, FSTEP, YMIN, YSCALEVZ, 'Z',
     >               FLABEL, TITLE, LEGENDC, LEGENDO, LEGENDW,
     >               XSHIFTLEG, YSHIFTLEG, FSMACHO, ALPHAO, REO,
     >               CLO, CDO, CYO, CPC, LOCALCFS, SIGNBLKO(N), CREFO,
     >               LINES, IPAGE, FIRST, LAST, LUNPLT, LUNSAVE)

               END IF

            END DO

C           Target Cps?

            IF (TARGETS) THEN

               LINES = .NOT. TARGSYM
               FIRST = .TRUE.
               LAST  = JROW == 0
               NIJ   = NITARG * NJTARG

               CALL GRAPHF
     >           ('T', .TRUE., NITARG, NJTARG, MAXSEG,
     >            XYZTARG(1+2*NIJ), XYZTARG(1+NIJ), XYZTARG(1),
     >            CPTARG, SEG, 1, NITARG, 1, NJTARG,
     >            XCUTS(I), ZAXES, YAXES, ZLEFT, ZRIGHT, ZSTEP,
     >            FBOT, FTOP, FSTEP, YMIN, YSCALEVZ, 'Z',
     >            FLABEL, TITLE, LEGENDC, LEGENDO, LEGENDW,
     >            XSHIFTLEG, YSHIFTLEG, FSMACHO, ALPHAO, REO,
     >            CLO, CDO, CYO, CPC, LOCALCFS, ONE, CREFO, LINES,
     >            IPAGE, FIRST, LAST, LUNPLT, LUNSAVE)

            END IF

C           Wind tunnel Cps?

            IF (JROW > 0) THEN

               LINES = .FALSE.
               FIRST = .TRUE.
               LAST  = .TRUE. ! Terminates the plot

               CALL GRAPHF
     >           ('W', .TRUE., MAXTAP, MAXROW, MAXSEG,
     >            XWT, XWT, XWT, CPWT, SEG, 1,
     >            NTAPS(JROW), JROW, JROW, XCUTS(I),
     >            ZAXES, YAXES, ZLEFT, ZRIGHT, ZSTEP, FBOT, FTOP, FSTEP,
     >            YMIN, YSCALEVZ, 'Z', FLABEL, TITLE, LEGENDC, LEGENDO,
     >            LEGENDW, XSHIFTLEG, YSHIFTLEG, FSMACHO, ALPHAO, REO,
     >            CLO, CDO, CYO, CPC, LOCALCFS, ONE, CREFO, LINES,
     >            IPAGE, FIRST, LAST, LUNPLT, LUNSAVE)

            END IF

         END DO

         CLOSE (LUNPLT)
         WRITE (LUNCRT, '(A)') ' Cross-stream plots completed: xcuts.ps'

         IF (SAVECUTS) THEN
            CLOSE (LUNPLT)
            WRITE (LUNCRT, '(A)')
     >         ' Plottable results also saved: xcuts.dat'
         END IF

      END IF


C     ********************* Plots at specified Zs **********************

      IF (NZCUTS > 0 .AND. PLTZCUTS) THEN

         OPEN (UNIT=LUNPLT, FILE='zcuts.ps', STATUS='UNKNOWN')

         IF (SAVECUTS) THEN

            OPEN (UNIT=LUNSAV, FILE='zcuts.dat', STATUS='UNKNOWN')

            I = LEN_TRIM (TITLE)
            WRITE (LUNSAV, '(A)') TITLE(1:I), 'Number of Z cuts'
            WRITE (LUNSAV, '(I3)') NZCUTS
            LUNSAVE = LUNSAV
         ELSE
            LUNSAVE = 0
         END IF

         IF (TUNNELZ) THEN

            OPEN (UNIT=LUNDAT, FILE=WTZFILE, STATUS='OLD')

            READ (LUNDAT, *) ! Skip title
            READ (LUNDAT, *)
            READ (LUNDAT, *) NROWS

            IF (NROWS > MAXROW) THEN
               WRITE (LUNCRT, '(1X, A, 2I7)')
     >          'Too many wind tunnel tap rows (Z cuts):', NROWS, MAXROW
               GO TO 999
            END IF

            WTMIN = HUGE (WTMIN) * 0.5
            WTMAX = -WTMIN

            DO J = 1, NROWS

               READ (LUNDAT, *)
               READ (LUNDAT, *) NTAPS(J), ZWT(J)

               IF (NTAPS(J) > MAXTAP) THEN
                  WRITE (LUNCRT, '(A, I3, /, A, 2I7)')
     >               ' Too many wind tunnel taps.  Row # (Z cuts):', J,
     >               ' # taps and limit:', NTAPS(J), MAXTAP
                  GO TO 999
               END IF

               READ (LUNDAT, *)
               READ (LUNDAT, *) (XWT(I,J), CPWT(I,J), I = 1, NTAPS(J))

               DO I = 1, NTAPS(J)
                  WTMIN = MIN (WTMIN, XWT(I,J))
                  WTMAX = MAX (WTMAX, XWT(I,J))
               END DO

            END DO

            CLOSE (LUNDAT)

            TOLWT = (WTMAX - WTMIN) * EPSWT

         END IF

         IPAGE = 0

         DO J = 1, NZCUTS

C           Current solution:

            LINES = .TRUE.

C           Check for wind tunnel Cps in order to end the plot if none:

            JROW = 0
            IF (TUNNELZ) THEN

               DO I = 1, NROWS
                  IF (ABS (ZCUTS(J) - ZWT(I)) < TOLWT) THEN
                     JROW = I
                     EXIT
                  END IF
               END DO

            END IF

            DO N = NFIRST, NLAST

               FIRST = N == NFIRST
               LAST  = N == NLAST .AND.
     >                 .NOT. (OVERLAY .OR. TARGETS .OR. JROW > 0)

               IF (SHOWX(N)) THEN

                  MAXI = NI(N)
                  MAXJ = NJ(N)
                  NIJ  = MAXI * MAXJ
                  I1   = MX(N)
                  I2   = MF(N)

                  CALL GRAPHF
     >              ('C', SHOWF(N), MAXI, MAXJ, MAXSEG,
     >               XYZ(I1), XYZ(I1+NIJ), XYZ(I1+2*NIJ),
     >               FUN(I2), SEG, 1, MAXI, 1, MAXJ,
     >               ZCUTS(J), XAXES, YAXES, XLEFT, XRIGHT, XSTEP,
     >               FBOT, FTOP, FSTEP, YMIN, YSCALEVX, 'X/C',
     >               FLABEL, TITLE, LEGENDC, LEGENDO, LEGENDW,
     >               XSHIFTLEG, YSHIFTLEG, FSMACH, ALPHA, RE,
     >               CL, CD, CY, CPC, SECCFS, SIGNBLK(N), CREF,
     >               LINES, IPAGE, FIRST, LAST, LUNPLT, LUNSAVE)

               END IF

            END DO

C           Overlaid solution?

            DO N = NFIRSTO, NLASTO

               FIRST = N == NFIRSTO
               LAST  = N == NLASTO .AND. .NOT. (TARGETS .OR. JROW > 0)

               IF (SHOWXO(N)) THEN

                  MAXI = NIO(N)
                  MAXJ = NJO(N)
                  NIJ  = MAXI * MAXJ
                  I1   = MXO(N)
                  I2   = MFO(N)

                  CALL GRAPHF
     >              ('O', SHOWFO(N), MAXI, MAXJ, MAXSEG,
     >               XYZO(I1), XYZO(I1+NIJ), XYZO(I1+2*NIJ),
     >               FUNO(I2), SEG, 1, MAXI, 1, MAXJ,
     >               ZCUTS(J), XAXES, YAXES, XLEFT, XRIGHT, XSTEP,
     >               FBOT, FTOP, FSTEP, YMIN, YSCALEVX, 'X/C',
     >               FLABEL, TITLE, LEGENDC, LEGENDO, LEGENDW,
     >               XSHIFTLEG, YSHIFTLEG, FSMACHO, ALPHAO, REO,
     >               CLO, CDO, CYO, CPC, SECCFSO, SIGNBLKO(N), CREFO,
     >               LINES, IPAGE, FIRST, LAST, LUNPLT, LUNSAVE)

               END IF

            END DO

C           Target Cps?

            IF (TARGETS) THEN

               LINES = .NOT. TARGSYM
               FIRST = .TRUE.
               LAST  = JROW == 0
               NIJ   = NITARG * NJTARG

               CALL GRAPHF
     >           ('T', .TRUE., NITARG, NJTARG, MAXSEG,
     >            XYZTARG(1), XYZTARG(1+NIJ), XYZTARG(1+2*NIJ),
     >            CPTARG, SEG, 1, NITARG, 1, NJTARG,
     >            ZCUTS(J), XAXES, YAXES, XLEFT, XRIGHT, XSTEP,
     >            FBOT, FTOP, FSTEP, YMIN, YSCALEVX, 'X/C',
     >            FLABEL, TITLE, LEGENDC, LEGENDO, LEGENDW,
     >            XSHIFTLEG, YSHIFTLEG, FSMACHO, ALPHAO, REO,
     >            CLO, CDO, CYO, CPC, .FALSE., CREFO, ONE, LINES,
     >            IPAGE, FIRST, LAST, LUNPLT, LUNSAVE)

            END IF

C           Wind tunnel Cps?

            IF (JROW > 0) THEN

               LINES = .FALSE.
               FIRST = .TRUE.
               LAST  = .TRUE. ! Terminates the plot

               CALL GRAPHF
     >           ('W', .TRUE., MAXTAP, MAXROW, MAXSEG,
     >            XWT, XWT, XWT, CPWT, SEG, 1,
     >            NTAPS(JROW), JROW, JROW, ZCUTS(J),
     >            XAXES, YAXES, XLEFT, XRIGHT, XSTEP, FBOT, FTOP, FSTEP,
     >            YMIN, YSCALEVX, 'X', FLABEL, TITLE, LEGENDC, LEGENDO,
     >            LEGENDW, XSHIFTLEG, YSHIFTLEG, FSMACHO, ALPHAO, REO,
     >            CLO, CDO, CYO, CPC, .FALSE., ONE, CREFO, LINES,
     >            IPAGE, FIRST, LAST, LUNPLT, LUNSAVE)

            END IF

         END DO

         CLOSE (LUNPLT)
         WRITE (LUNCRT, '(A)') '    Chordwise plots completed: zcuts.ps'

         IF (SAVECUTS) THEN
            CLOSE (LUNPLT)
            WRITE (LUNCRT, '(A)')
     >         ' Plottable results also saved: zcuts.dat'
         END IF

      END IF


  999 CONTINUE

      END PROGRAM MULTIPLOT


C***********************************************************************
C
      SUBROUTINE FORCES (COSA, SINA, SREF, MAXI, MAXJ, I1, I2, J1, J2,
     >                   XYZ, CP, CL, CD, CY)
C
C     FORCES integrates lift, drag, and yaw coefficients for one surface
C     given Cps at the grid points (right-handed coordinate system).
C
C     11/23/96  D.A.Saunders  Adaptation of SYN87 code for MULTIPLOT.
C
C***********************************************************************

      IMPLICIT NONE

C     Arguments:

      INTEGER
     >   MAXI, MAXJ, I1, I2, J1, J2
      REAL
     >   COSA, SINA, SREF, XYZ(MAXI,MAXJ,3), CP(MAXI,MAXJ), CL, CD, CY

C     Local variables:

      INTEGER
     >   I, J, L, M
      REAL
     >   AX, AY, AZ, CT

C     Execution:

      CL = 0.
      CD = 0.
      CY = 0.

      DO J = J1 + 1, J2
         M = J - 1
         DO I = I1 + 1, I2
            L  = I - 1
            AX = (XYZ(I,J,2) - XYZ(L,M,2)) * (XYZ(L,J,3) - XYZ(I,M,3)) -
     >           (XYZ(I,J,3) - XYZ(L,M,3)) * (XYZ(L,J,2) - XYZ(I,M,2))
            AY = (XYZ(I,J,3) - XYZ(L,M,3)) * (XYZ(L,J,1) - XYZ(I,M,1)) -
     >           (XYZ(I,J,1) - XYZ(L,M,1)) * (XYZ(L,J,3) - XYZ(I,M,3))
            AZ = (XYZ(I,J,1) - XYZ(L,M,1)) * (XYZ(L,J,2) - XYZ(I,M,2)) -
     >           (XYZ(I,J,2) - XYZ(L,M,2)) * (XYZ(L,J,1) - XYZ(I,M,1))
            CT = -0.125 * (CP(L,M) + CP(I,M) + CP(L,J) + CP(I,J))
            CD = CT * AX + CD
            CL = CT * AY + CL
            CY = CT * AZ + CY
         END DO
      END DO

      CL = CL / SREF
      CD = CD / SREF
      CY = CY / SREF
      CT = CL * COSA - CD * SINA
      CD = CD * COSA + CL * SINA
      CL = CT

      END SUBROUTINE FORCES


C***********************************************************************
C
      SUBROUTINE FPLOT3D (IFUN, FSMACH, ALPHA, REYNO, MAXI, MAXJ, MAXK,
     >                    I1, I2, J1, J2, K1, K2, Q, FUN)
C
C     FPLOT3D evaluates the scalar function defined by IFUN, which
C     should match the values known to PLOT3D.
C
C     11/22/96  D.A.Saunders  Initial implementation for MULTIPLOT.
C
C***********************************************************************

      IMPLICIT NONE

C     Arguments:

      INTEGER
     >   IFUN, MAXI, MAXJ, MAXK, I1, I2, J1, J2, K1, K2
      REAL
     >   FSMACH, ALPHA, REYNO, Q(MAXI,MAXJ,MAXK,5), FUN(MAXI,MAXJ,MAXK)

C     Local constants:

      REAL, PARAMETER ::
     >   CINF    = 1.,         ! Freestream speed of sound
     >   GAMMA   = 1.4,
     >   GAMMAM1 = GAMMA - 1.,
     >   PINF    = 1./ GAMMA,  ! Freestream pressure
     >   RHOINF  = 1.          ! Freestream density

C     Local variables:

      INTEGER
     >   I, J, K
      REAL
     >   DENOM, P, RHO

C     Execution:

      IF (IFUN == 114) THEN ! Cp:  Reference: PLOT3D Manual p.130

         DENOM = 2./ (RHOINF * (FSMACH * CINF) ** 2)

         DO K = K1, K2
            DO J = J1, J2
               DO I = I1, I2
                  RHO = SIGN (MAX (1.E-9, ABS (Q(I,J,K,1))), Q(I,J,K,1))
                  P = (Q(I,J,K,2)**2 + Q(I,J,K,3)**2 + Q(I,J,K,4)**2) *
     >                (0.5 / RHO)
                  P = GAMMAM1 * DIM (Q(I,J,K,5), P) ! Pressure >= 0.
                  FUN(I,J,K) = (P - PINF) * DENOM
               END DO
            END DO
         END DO

      ELSE
         WRITE (6, '(/, A, I6)') ' Unknown function code:', IFUN
         STOP
      END IF

      END SUBROUTINE FPLOT3D


C***********************************************************************
C
      SUBROUTINE GRAPHF (CASE, SHOWF, MAXI, MAXJ, MAXSEG, X, Y, Z, F,
     >                   SEG, I1, I2, J1, J2, ZSLICE, XAXES, YAXES,
     >                   XLEFT, XRIGHT, XSTEP, FBOT, FTOP, FSTEP, YMIN,
     >                   YSCALE, XLABEL, FLABEL, TITLE, LEGENDC,
     >                   LEGENDO, LEGENDW, XSHIFTLEG, YSHIFTLEG,
     >                   FSMACH, ALPHA, RE, CL, CD, CY, CPC, SECCFS,
     >                   SIGNBLK, CREF, LINES, IPAGE, FIRST, LAST,
     >                   LUNPLT, LUNSAVE)
C
C     GRAPHF prepares & plots one cut of a scalar function on a surface
C     at a specified "Z" (which may be Z or X).  The function F will be
C     Cp normally, but some effort has been made to allow for other Fs.
C     GRAPHF may be called more than once for a given plot page, (a) if
C     there are two or more surfaces (see FIRST and LAST usage), or (b)
C     to overlay another result and/or target Cp data and/or wind tunnel
C     data (see CASE usage).
C
C     11/27/96  D.A.Saunders  Adaptation of SYN87 modules for MULTIPLOT,
C                             pushing the generality of PLCUT one level
C                             higher.
C     03/31/98      "         Eliminated (u,v) approach in favor of
C                             unstructured surface techniques (CUTSURF).
C     04/16/98      "         Handled wind tunnel data.
C     01/22/99      "         Provided for suppressing the function,
C                             moving the legend, etc.
C
C***********************************************************************

      IMPLICIT NONE

C     Arguments:

      INTEGER
     >   MAXI, MAXJ, MAXSEG, I1, I2, J1, J2, IPAGE, LUNPLT, LUNSAVE
      REAL
     >   X(MAXI,MAXJ), Y(MAXI,MAXJ), Z(MAXI,MAXJ), F(MAXI,MAXJ),
     >   SEG(4,2*MAXSEG)
      REAL
     >   ZSLICE, XLEFT, XRIGHT, XSTEP, FBOT, FTOP, FSTEP,
     >   YMIN, YSCALE, XSHIFTLEG, YSHIFTLEG, FSMACH, ALPHA,
     >   RE, CL, CD, CY, CPC, SIGNBLK, CREF
      LOGICAL
     >   SHOWF, XAXES, YAXES, SECCFS, LINES, FIRST, LAST
      CHARACTER
     >   CASE*1, LINTYP*8, XLABEL*(*), FLABEL*(*), TITLE*(*),
     >   LEGENDC*(*), LEGENDO*(*), LEGENDW*(*)

C     Global variables:

      INTEGER      JWRIT
      COMMON /PSC/ JWRIT

C     Local constants:

      INTEGER, PARAMETER ::
     >   NADDIM  = 20      ! 1 + max. # line segments expected per cut

      REAL, PARAMETER ::   ! These are inches from the origin of non-offset axes
     >   AXSHIFT = -0.40,  ! Axes are offset down and left
     >   DY      =  0.02,  ! Raises legend symbols a shade
     >   DYLEG   =  0.20,  ! Distance between legend lines
     >   DYTITLE =  0.20,  ! Distance between caption lines
     >   XLEG1   =  1.70,  ! Default left end of legend symbol/line segment
     >   XLEG2   =  1.95,  ! Default right ...
     >   XLEG3   =  2.20,  ! Default left end of legend text
     >   YLEG1   =  2.20,  ! Default vertical location of first legend + DYLEG
     >   YLEGY   =  1.20,  ! Default vertical location of Y stretch factor text
     >   XLENGTH =  5.0,   ! Axis lengths in inches
     >   YLENGTH =  8.0,
     >   YTITLE  = -1.4,   ! Main title location below end of Y axis
     >   HALF    =  0.5,
     >   ONE     =  1.0,
     >   ZERO    =  0.

      CHARACTER, PARAMETER ::
     >   CURRENT = 'C', OVERLAY = 'O', TARGET = 'T', TUNNEL = 'W',
     >   DASH*4  = 'dash', LONG*8 = 'longdash', SOLID*5 = 'solid',
     >   SYMFONT*6 = 'Symbol', TXTFONT*11 = 'Times-Roman'

C     Local variables:

      INTEGER
     >   NAD(NADDIM), I, L, N, NCRV, NPTS
      REAL
     >   CENTERS(10), XCRV(MAXSEG), YCRV(MAXSEG), ZCRV(MAXSEG),
     >   CPCRV(MAXSEG),
     >   ANGLE, CDLOC, CLLOC, CMLOC, CLOAD, DLOAD, CPA, CHORD, DCD, DCL,
     >   FORIGIN, FSCALE, SA, CA, XA, YA, XLE, YLE, XMIN, XOFF,
     >   XPT, YPT, XSCALE, YLEG, YSHIFT
      LOGICAL
     >   DOLEGENDO, DOLEGENDT, DOYSCALE, LOCALCFS
      CHARACTER
     >   B*28, CUTTYPE*1, SYMB*1, XLAB*3

C     Procedures:

      REAL
     >   COSD, SIND

C     Storage:

      DATA       ! Heuristic placement of caption values
     >   CENTERS /-0.2, 0.3, 0.9, 1.6, 2.2, 2.9, 3.6, 4.2, 4.8, 5.4/
      SAVE
     >   CENTERS, DOLEGENDO, DOLEGENDT, FORIGIN, FSCALE, LOCALCFS,
     >   XMIN, XOFF, XSCALE, YLEG, YSHIFT

C     Execution:

      JWRIT = LUNPLT

      IF (FIRST) THEN
         DOLEGENDO = CASE == OVERLAY .AND. LEGENDO /= ' '
         DOLEGENDT = CASE == TARGET
         LOCALCFS  = SECCFS .AND. SHOWF
      END IF

      IF (FIRST .AND. CASE == CURRENT) THEN ! First block of primary dataset

         YLEG     = YLEG1 + YSHIFTLEG
         DOYSCALE = .TRUE.
         XSCALE   = XLENGTH / (XRIGHT - XLEFT)
         XMIN     = XLEFT         ! Because PLCUT changes its XMIN argument
         XLAB     = XLABEL

         IF (XLAB == 'Z') THEN
            CUTTYPE = 'X'
            YSHIFT  = 0.25        ! Inches, applied to all geometry points
            L = 1                 ! Label length
         ELSE  ! XLABEL = 'X/C'
            CUTTYPE = 'Z'
            YSHIFT  = ZERO
            IF (XLEFT == ZERO .AND. XRIGHT == ONE) THEN
               XSCALE = -XLENGTH  ! PLCUT adjusts it; SAVE keeps it for overlay
               L = 3
            ELSE                  ! Avoid 'X/C' when wing surface is suppressed
               L = 1              ! XLAB = 'X'
               LOCALCFS = .FALSE. ! Suppress section coefs. if X is not X/C
            END IF
         END IF

         CALL INITPL (2., 2.4, ONE, IPAGE) ! (2,2.4)" <-> plotted (X,Y) = (0,0)

         IF (XAXES) THEN

            CALL XAXIS (XLEFT, XRIGHT, XSTEP, ZERO, AXSHIFT, XLENGTH,
     >                  -1, -1, XLAB, L, TXTFONT, 8, 0)
         END IF

         IF (YAXES) THEN
            L = LEN (FLABEL)
            CALL YAXIS (FBOT, FTOP, FSTEP, AXSHIFT, ZERO, YLENGTH,
     >                  -1, -1, FLABEL, L, TXTFONT, 8, 0)
         END IF

         FSCALE  = YLENGTH / (FTOP - FBOT)
         FORIGIN = -FBOT * FSCALE

         IF (CPC <= FBOT .AND. CPC >= FTOP) THEN
            XPT = AXSHIFT
            YPT = FORIGIN + FSCALE * CPC
            CALL SYMTYPE (TXTFONT, 0.12)
            CALL SYMBOL (XPT, YPT, '--- Cp* ---', 11, ZERO)
         END IF

C        Center the title:

         CALL SYMTYPE (TXTFONT, 0.20)

         DO I = LEN (TITLE), 1, -1
            IF (TITLE(I:I) /= ' ') EXIT
         END DO

         XPT = XLENGTH * HALF
         YPT = YTITLE
         CALL CTEXT (XPT, YPT, TITLE(1:I),  I)

C        Six or ten caption headers depending on type of cut:

         CALL SYMTYPE (TXTFONT, 0.15)
         YPT  = YPT - DYTITLE
         XOFF = 1.1 ! Offset to ~ center just six quantities
         IF (LOCALCFS) XOFF = -0.15

         CALL CTEXT (CENTERS(1) + XOFF, YPT, 'Mach',  4)
         CALL CTEXT (CENTERS(2) + XOFF, YPT, 'Alpha', 5)
         CALL CTEXT (CENTERS(3) + XOFF, YPT, 'Re',    2)
         CALL CTEXT (CENTERS(4) + XOFF, YPT, 'CL',    2)
         CALL CTEXT (CENTERS(5) + XOFF, YPT, 'CD',    2)
         CALL CTEXT (CENTERS(6) + XOFF, YPT, CUTTYPE, 1)

         IF (LOCALCFS) THEN
            CALL CTEXT (CENTERS(7)  + XOFF, YPT, 'Cl',   2)
            CALL CTEXT (CENTERS(8)  + XOFF, YPT, 'Cd',   2)
            CALL CTEXT (CENTERS(9)  + XOFF, YPT, 'Cm',   2)
            CALL CTEXT (CENTERS(10) + XOFF, YPT, 'Load', 4)
         END IF

         N = LEN_TRIM (LEGENDC)

         IF (N > 0) THEN
            YLEG = YLEG - DYLEG

            CALL LINETYPE (SOLID)
            CALL PLOT (XLEG1 + XSHIFTLEG, YLEG, 3) ! moveto
            CALL PLOT (XLEG2 + XSHIFTLEG, YLEG, 2) ! lineto
            CALL SYMTYPE (TXTFONT, 0.18)
            CALL SYMBOL (XLEG3 + XSHIFTLEG, YLEG, LEGENDC, N, ZERO)
         END IF

      END IF

C     First 6 caption values?

      IF (FIRST .AND. CASE /= TARGET .AND. CASE /= TUNNEL) THEN

         YPT = YTITLE - DYTITLE - DYTITLE
         IF (CASE == OVERLAY) YPT = YPT - DYTITLE

         CALL SYMTYPE (TXTFONT, 0.15)

         WRITE (B(1:5), '(F5.3)') FSMACH
         CALL CTEXT (CENTERS(1) + XOFF, YPT, B, 5)

         WRITE (B(1:6), '(F6.3)') ALPHA
         CALL CTEXT (CENTERS(2) + XOFF, YPT, B, 6)

         WRITE (B(1:9), '(1P, E9.2)') RE
         CALL CTEXT (CENTERS(3) + XOFF, YPT, B, 9)

         WRITE (B(1:8), '(F8.5)') CL
         CALL CTEXT (CENTERS(4) + XOFF, YPT, B, 8)

         WRITE (B(1:8), '(F8.5)') CD
         CALL CTEXT (CENTERS(5) + XOFF, YPT, B, 8)

         WRITE (B(1:8), '(F8.3)') ZSLICE
         CALL CTEXT (CENTERS(6) + XOFF, YPT, B, 8)

      END IF


C     Slice & plot the surface data:

      IF (CASE == CURRENT) THEN
         LINTYP = SOLID
         SYMB   = '+'
      ELSE IF (CASE == OVERLAY) THEN
         LINTYP = LONG
         SYMB   = 'o'
      ELSE IF (CASE == TARGET) THEN
         LINTYP = DASH
         SYMB   = '*'
      ELSE   ! CASE == 'W' for wind tunnel data
         SYMB   = '+'
      END IF

      IF (LINES) THEN
         ANGLE = -999. ! Angle < 0. is the switch within PLCUT
      ELSE ! symbols
         ANGLE = ZERO
      END IF

      CALL PLCUT (CASE, SHOWF, ZSLICE, MAXI, MAXJ, X, Y, Z, F,
     >            I1, I2, J1, J2, MAXSEG, SEG, NADDIM, NAD,
     >            XCRV, YCRV, ZCRV, CPCRV,
     >            XMIN, YMIN, XSCALE, YSHIFT, YSCALE, LINTYP,
     >            FORIGIN, FSCALE, SYMB, ANGLE, LUNSAVE)


      IF (NAD(1) > 0) THEN

         CALL LINETYPE (LINTYP)

         IF (DOYSCALE) THEN
            DOYSCALE = .FALSE.
            IF (YSCALE /= ONE) THEN
               CALL SYMTYPE (TXTFONT, 0.18)
               WRITE (B(1:28), '(A, F6.1)')
     >            '     Y stretch factor:', YSCALE
               XPT = XLEG2 + XSHIFTLEG
               YPT = YLEGY + YSHIFTLEG
               CALL SYMBOL (XPT, YPT, B, 28, ZERO)
            END IF
         END IF

         IF (DOLEGENDO) THEN

            DOLEGENDO = .FALSE.
            YLEG = YLEG - DYLEG

            CALL PLOT (XLEG1 + XSHIFTLEG, YLEG, 3) ! moveto
            CALL PLOT (XLEG2 + XSHIFTLEG, YLEG, 2) ! lineto
            CALL SYMTYPE (TXTFONT, 0.18)

            N = LEN_TRIM (LEGENDO)
            CALL SYMBOL (XLEG3 + XSHIFTLEG, YLEG, LEGENDO, N, ZERO)

            IF (.NOT. LINES) THEN
               CALL SYMTYPE (SYMFONT, 0.18)
               CALL SYMBOL (XLEG1 + XSHIFTLEG, YLEG + DY, 'o ', 2, ZERO)
            END IF          ! Left justify 'o'

         ELSE IF (DOLEGENDT) THEN

            DOLEGENDT = .FALSE.
            YLEG = YLEG - DYLEG

            CALL SYMTYPE (TXTFONT, 0.18)
            CALL SYMBOL (XLEG2 + XSHIFTLEG, YLEG, '     Target',11,ZERO)
            CALL PLOT (XLEG1 + XSHIFTLEG, YLEG, 3)
            CALL PLOT (XLEG2 + XSHIFTLEG, YLEG, 2)

            IF (.NOT. LINES) THEN
               CALL SYMTYPE (SYMFONT, 0.18)
               CALL SYMBOL (XLEG1 + XSHIFTLEG, YLEG + DY, '* ', 2, ZERO)
            END IF

         ELSE IF (CASE == TUNNEL) THEN

            YLEG = YLEG - DYLEG
            N = LEN_TRIM (LEGENDW)

            CALL SYMTYPE (TXTFONT, 0.18)
            CALL SYMBOL (XLEG3 + XSHIFTLEG, YLEG, LEGENDW, N, ZERO)
            CALL SYMTYPE (SYMFONT, 0.18)
            XPT = XLEG1 + XSHIFTLEG
            CALL SYMBOL (XPT, YLEG + DY, '+ ', 2, ZERO) ! Left justify '+'

         END IF

      END IF


C     Section force coefficients?

      IF (NAD(1) == 1 .AND. LOCALCFS .AND. SHOWF) THEN ! Wing sec. (presumably)

         NPTS = NAD(2) - 1
         XLE  = XCRV(1)
         L = 1

         DO I = 2, NPTS
            IF (XCRV(I) < XLE) THEN
               XLE = XCRV(I)
               L = I
            END IF
         END DO

         YLE = YCRV(L)
         CHORD = MAX (XCRV(1), XCRV(NPTS)) - XLE

         CLLOC = ZERO
         CDLOC = ZERO
         CMLOC = ZERO

         DO I = 2, NPTS
            CPA = (CPCRV(I) + CPCRV(I-1)) * HALF
            DCL =-(XCRV(I)  -  XCRV(I-1)) * CPA
            XA  = (XCRV(I)  +  XCRV(I-1)) * HALF
            DCD = (YCRV(I)  -  YCRV(I-1)) * CPA
            YA  = (YCRV(I)  +  YCRV(I-1)) * HALF
            CLLOC = CLLOC + DCL
            CDLOC = CDLOC + DCD
            CMLOC = CMLOC + DCD*(YA - YLE) - DCL*(XA - XLE)
         END DO

         CLLOC = CLLOC / CHORD
         CDLOC = CDLOC / CHORD
         CMLOC = CMLOC / CHORD ** 2
         SA    = SIND (ALPHA)
         CA    = COSD (ALPHA)
         XA    = CLLOC * CA - CDLOC * SA
         CDLOC = CLLOC * SA + CDLOC * CA
         CLLOC = XA
         CLOAD = (CHORD / CREF) * CLLOC
         DLOAD = (CHORD / CREF) * CDLOC

         IF (SIGNBLK /= ONE) THEN
            CLLOC = -CLLOC
            CDLOC = -CDLOC
            CMLOC = -CMLOC
            CLOAD = -CLOAD
            DLOAD = -DLOAD
         END IF

         CALL SYMTYPE (TXTFONT, 0.15)

         YPT = YTITLE - DYTITLE - DYTITLE
         IF (CASE == OVERLAY) YPT = YPT - DYTITLE

         WRITE (B(1:8), '(F8.5)') CLLOC
         CALL CTEXT (CENTERS(7) + XOFF, YPT, B, 8)

         WRITE (B(1:8), '(F8.5)') CDLOC
         CALL CTEXT (CENTERS(8) + XOFF, YPT, B, 8)

         WRITE (B(1:8), '(F8.5)') CMLOC
         CALL CTEXT (CENTERS(9) + XOFF, YPT, B, 8)

         WRITE (B(1:8), '(F8.5)') CLOAD
         CALL CTEXT (CENTERS(10) + XOFF, YPT, B, 8)

      END IF

      IF (LAST) CALL ENDPLT

      END SUBROUTINE GRAPHF

C***********************************************************************
C
      SUBROUTINE PLCUT (CASE, SHOWF, ZSLICE, IDIM, JDIM, X, Y, Z, CP,
     >                  I1, I2, J1, J2, MAXSEG, SEG, NADDIM, NAD,
     >                  XCRV, YCRV, ZCRV, CPCRV,
     >                  XMIN, YMIN, SCALE, YSHIFT, YSCALE,
     >                  LINTYP, CPORG, CPSCL, SYMB, ANGLE, LUNSAVE)

C     PLCUT prepares & plots a slice of surface data consisting of
C     regular (X,Y,Z) coordinates, and a single scalar function, CP.
C     ZSLICE is considered to be the value of the THIRD coordinate
C     specified (Z at this level). The plotting is interleaved with
C     the data preparation partly because the original contouring
C     approach could produce an indefinite number of line segments,
C     and partly to modularize the plotting (Y and CP vs. X at this
C     level) for repeated application to various datasets and sub-surfaces.
C
C     03/08/96 D.Saunders Generalized parts of original GRAPHX.
C     03/11/96     "      SCALE < 0. causes XY-scaling here;
C                         ANGLE < 0. means lines for CP, not symbols
C     07/26/96     "      YSCALE allows magnifying Y coordinates.
C     03/31/98     "      Switched from CONT2D to CUTSURF approach.
C     04/16/98     "      Wind tunnel data required CASE argument.
C     08/21/98     "      CUTSURF & CUTORDER were upgraded to Fortran 90.
C     01/22/99     "      SHOWF argument allows suppressing the function.
C     05/13/99     "      LUNSAVE /= 0 allows saving cut pts. to a file.
C
C***********************************************************************

      IMPLICIT NONE

C     Arguments:

      CHARACTER
     >   CASE*1
      LOGICAL
     >   SHOWF
      INTEGER
     >   IDIM, JDIM, I1, I2, J1, J2, MAXSEG, LUNSAVE
      REAL
     >   X(IDIM,JDIM), Y(IDIM,JDIM), Z(IDIM,JDIM), CP(IDIM,JDIM),
     >   SEG(4,2*MAXSEG)
      INTEGER
     >   NADDIM, NAD(NADDIM)
      REAL
     >   XCRV(MAXSEG), YCRV(MAXSEG), ZCRV(MAXSEG), CPCRV(MAXSEG),
     >   ZSLICE, XMIN, YMIN, SCALE, YSHIFT, YSCALE, CPORG, CPSCL, ANGLE
      CHARACTER
     >   LINTYP*(*), SYMB*1

C     Local constants:

      REAL, PARAMETER ::
     >   EPS = 2.E-7,  ! Fraction of data range for CUTSURF
     >   ONE = 1., ZERO = 0.
      CHARACTER, PARAMETER ::
     >   CURRENT*1 = 'C', SYMFONT*6 = 'Symbol', TUNNEL*1  = 'W'

C     Automatic arrays:

      INTEGER    NDGCUT(2,MAXSEG) ! Pointers used by CUTSURF & CUTORDER

C     Local variables:

      INTEGER    IER, L, L1, L2, LINE, NCRV, NNODE, NPTS, NSEG
      REAL       DATAMAX, DATAMIN, DRANGE, EPSCP, EPSXYZ, P(3), S(3),
     >           SCALEY, X0, XMAX, XPT, Y0, YPT

C     Execution:

      IF (CASE /= TUNNEL) THEN

C        Determine the data range to provide CUTSURF with tolerances:

         DRANGE  = ZERO
         DATAMAX = X(I1,J1) ! Necessary initialization
         DATAMIN = DATAMAX

         CALL BOUNDS (I2, J2, IDIM, X, DATAMIN, DATAMAX)
         DRANGE = MAX (DRANGE, DATAMAX - DATAMIN)

         CALL BOUNDS (I2, J2, IDIM, Y, DATAMIN, DATAMAX)
         DRANGE = MAX (DRANGE, DATAMAX - DATAMIN)

         CALL BOUNDS (I2, J2, IDIM, Z, DATAMIN, DATAMAX)
         DRANGE = MAX (DRANGE, DATAMAX - DATAMIN)

         EPSXYZ = EPS * DRANGE

         DATAMAX = CP(I1,J1)
         DATAMIN = DATAMAX

         CALL BOUNDS (I2, J2, IDIM, CP, DATAMIN, DATAMAX)

         EPSCP = EPS * (DATAMAX - DATAMIN)

C        Determine the 2-pt. line segments forming the cut at Z = ZSLICE:

         P(1) = ZERO   ! Point
         P(2) = ZERO
         P(3) = ZSLICE

         S(1) = ZERO   ! Vector from the point defining the cutting plane
         S(2) = ZERO
         S(3) = ONE

         NSEG  = 0     ! Reqd. on input
         NNODE = 0

         CALL CUTSURF (X, Y, Z, CP, I1, I2, J1, J2, IDIM, JDIM,
     >                 SEG, NNODE, 2*MAXSEG, NDGCUT, NSEG, MAXSEG,
     >                 P, S, EPSXYZ, EPSCP, IER)

         IF (IER /= 0) THEN
            WRITE (6, '(/,A,I2,A,F10.3,A,I5,/,A)')
     >         ' PLCUT: IER from CUTSURF =', IER,
     >         '   Cut =', ZSLICE, '   NSEG =', NSEG, ' Proceeding.'
         END IF

         NAD(1) = NSEG ! May be zero; caller may need to know

         IF (NSEG == 0) GO TO 999 ! Avoids indenting


C        Convert the 2-pt. line segments to one or more contiguous curves:

         CALL CUTORDER (NAD, NADDIM, XCRV, YCRV, ZCRV, CPCRV,
     >                  MAXSEG, SEG, NNODE, 2*MAXSEG, NDGCUT,
     >                  NSEG, MAXSEG, EPSXYZ, EPSCP, IER)

         IF (IER /= 0) THEN
            WRITE (6, '(/,A,I2,A,F10.3,A,I5,/,A)')
     >         ' PLCUT: IER from CUTORDER =', IER, '   Cut =',
     >         ZSLICE, '   NSEG =', NSEG, ' Aborting this cut.'
            NAD(1) = 0
            GO TO 999
         END IF

         NCRV = NAD(1)

      ELSE ! Wind tunnel data case - simply copy the data as one curve

         NCRV   = 1
         L2     = I2 - I1 + 1
         NAD(2) = L2 + 1

         XCRV (1:L2) = X (I1:I2, J1)
         CPCRV(1:L2) = CP(I1:I2, J1)
         
      END IF


C     Self-scaling option:

      IF (SCALE < ZERO) THEN

         XMAX = XCRV(1)
         XMIN = XMAX
         YMIN = YCRV(1)
         L2 = NAD(NCRV+1) - 1

         DO L = 2, L2
            XMAX = MAX (XCRV(L), XMAX)
            XMIN = MIN (XCRV(L), XMIN)
            YMIN = MIN (YCRV(L), YMIN)
         END DO

         SCALE = -SCALE / (XMAX - XMIN)

      END IF

      SCALEY = SCALE * YSCALE


C     Plot the distinct curves:

      CALL LINETYPE (LINTYP)

      NAD(1) = 1 ! So the following works for line 1

      DO LINE = 1, NCRV
         L1 = NAD(LINE)
         L2 = NAD(LINE+1) - 1

         IF (CASE /= TUNNEL) THEN

C           Plot the geometry:

            XPT = (XCRV(L1) - XMIN) * SCALE
            YPT = (YCRV(L1) - YMIN) * SCALEY + YSHIFT

            CALL PLOT (XPT, YPT, 3) ! moveto

            DO L = L1 + 1, L2
               XPT = (XCRV(L) - XMIN) * SCALE
               YPT = (YCRV(L) - YMIN) * SCALEY + YSHIFT
               CALL PLOT (XPT, YPT, 2) ! lineto
            END DO

         END IF

C        Plot the scalar function?

         IF (SHOWF) THEN

            IF (ANGLE >= ZERO) THEN ! Use symbols

               CALL SYMTYPE (SYMFONT, 0.18)

               DO L = L1, L2
                  XPT = (XCRV(L) - XMIN) * SCALE
                  YPT = CPORG + CPSCL * CPCRV(L)
                  CALL SYMBOL (XPT, YPT, SYMB, 1, ANGLE)
               END DO

            ELSE

               XPT = (XCRV(L1) - XMIN) * SCALE
               YPT = CPORG + CPSCL * CPCRV(L1)

               CALL PLOT (XPT, YPT, 3) ! moveto

               DO L = L1 + 1, L2
                  XPT = (XCRV(L) - XMIN) * SCALE
                  YPT = CPORG + CPSCL * CPCRV(L)
                  CALL PLOT (XPT, YPT, 2) ! lineto
               END DO

            END IF

         END IF

C        Save the cut points for plotting elsewhere?

         IF (LUNSAVE /= 0) THEN

            IF (CASE == CURRENT) THEN

               NPTS = L2 - L1 + 1
               WRITE (LUNSAVE, '(A, /, I4, F12.6)')
     >            'NPTS      ZSLICE', NPTS, ZSLICE

               IF (SHOWF) THEN
                  WRITE (LUNSAVE, '(A)')
     >               '       X           Y           Cp'
                  WRITE (LUNSAVE, '(3F12.6)')
     >               (XCRV(L), YCRV(L), CPCRV(L), L = L1, L2)
               ELSE
                  WRITE (LUNSAVE, '(A)')
     >               '       X           Y'
                  WRITE (LUNSAVE, '(2F12.6)')
     >               (XCRV(L), YCRV(L), L = L1, L2)
               END IF

            END IF

         END IF

      END DO ! Next line segment

      NAD(1) = NCRV

  999 RETURN

      END SUBROUTINE PLCUT
