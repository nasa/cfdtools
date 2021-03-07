!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      PROGRAM GEN1D
!
!  PURPOSE:
!
!     GEN1D is intended to generate data that can be used for evaluating 1-D
!     interpolation/smoothing routines.  It provides a menu of likely functions
!     of one variable (easily modified as necessary).  It also provides the
!     option to evaluate the specified function at a self-generated set of grid
!     points or at a set read from a file.  Results are written to 'gen1d.dat'
!     in the following (Tecplotable) format:
!
!                      # Title with function details
!                      # NPTS
!                      X(1)     Y(1)
!                      X(2)     Y(2)
!                       :        :
!                      X(NPTS)  Y(NPTS)
!
!     If a grid is read from a file, its format is the same except that it has
!     1 column, not 2.  (Columns after the first will be ignored.)
!
!     Provision is made for parametric functions in the form X = X(T), Y = Y(T).
!     Here, reading a T distribution is an unlikely requirement and is NOT
!     provided for.
!
!     An option for applying noise to the Y values (or to both X and Y for the
!     parametric functions) is also provided.
!
!  METHOD:
!
!        Select function of X from menu.
!        Prompt for type of grid to use:
!        IF <self-generated> THEN
!           Prompt for XMIN, XMAX, DX.
!           Generate points X(*).
!        ELSE
!           Prompt for grid file name and read X(*).
!        END IF
!        Prompt for additional constants and evaluate function at all pts.
!        Write them to the output file.
!
!  FILES USED:

!    LUNCRT     O     For prompts and error messages
!    LUNKBD     I     For user responses to prompts
!    LUNX       I     For input X coordinates (if any)
!    LUNXY      O     For generated (X,Y) values
!
!  PROCEDURES:
!
!    GAUSS   NUMODULES  Gaussian random deviate generator (for noise)
!    QUINTIC INTERPLIB  Generates the quintic between a (x,y,y',y") pair
!    READER  PROGTOOLS  Prompting utility
!
!  HISTORY:
!
!    DAS   10/17/86   Adapted from GEN2D
!    DAS   02/18/87   Provided for X=X(T), Y=Y(T) parametric functions
!    DAS   11/20/87   Installed option to apply noise to the data
!    DAS   06/05/91   Added sin(x)/x and sinh(x)/x options.
!    DAS   03/21/92   Added parametric form of the parabola.
!    DAS   01/09/93   Made the cubic a quartic.
!    DAS   08/22/97   Added Hermite cubic on [0, 1]
!    DAS   09/27/00   Replaced IMSL's GGNML with GAUSS/RAN3 found on the web.
!    DAS   10/26/01   Added QUINTIC.
!    DAS   10/29/03   Slight changes for 64-bit version.
!    DAS   07/19/04   Use of $ carriage control was obsolete.
!    DAS   10/27/10   Added f(x) = Aexp(Bx) + C option; Fortran 90 upgrades.
!    DAS   11/15/10   Generalized previous option to handle Aexp(B/x) + C too.
!    DAS   11/22/10   Added Ax**n*exp(Bx) and power law option:  Ax**B + C.
!    DAS   05/24/11   Outputs are 64-bit precision now, not single precision.
!    DAS   10/06/11   Added sin**n(x.pi/2) for x in [0, 1] & 2 catenary options.
!    DAS   10/13/11   Replaced explicit discretization of a cosh (x/a) - a with
!                     a call to catenary grid (left, right or both halves)
!                     unless a desired parameter a is input.
!    DAS   05/16/14   Added Gaussian and normal distribution options, along with
!                     the error function-related options.
!
!  AUTHOR:  David Saunders, Sterling Software/ELORET/ERC at NASA Ames.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IMPLICIT NONE

!     Local constants:

      INTEGER, PARAMETER :: &
         LUNX=2, LUNXY=4, LUNKBD=5, LUNCRT=6, MXMENU=20, MXPTS=50000

      REAL, PARAMETER :: &
         ONE = 1.E+0, HALF = 0.5E+0, ZERO = 0.E+0

!     Local variables:

      INTEGER :: &
         I, IER, IFUN, IOS, ISEED, MODE, N, NPTS

      REAL :: &
         COEFS(5), X(MXPTS), Y(MXPTS), &
         A, B, C, D, DNOISE, E, DEGRAD, DT, DX, D1, D2, H, MU, PERCENT, PI, &
         SDEV, SIGMA, T, TMAX, X1, X2, XMAX, XMIN, XR, Y1, Y2, YMAX, YP1, YP2, &
         YPP1, YPP2

      LOGICAL :: &
         ARRHENIUS, CR, EOF, PARAMETRIC, XINDEGREES, YES

      CHARACTER :: &
         ANSWER*1, DATAFILE*50, DETAILS*70, MENU(0:MXMENU)*64

!     Procedures:

      REAL :: &
         ERF, ERFC, ERFCX, GAUSS
      EXTERNAL :: &
         GAUSS, GAUSSIAN, NORMAL_DISTRIBUTION, QUINTIC, &
         READC, READI, READR, READY

!     Storage:

      DATA   MENU   / &
         '  0:  Help                                                      ',   &
         '  1:  f(x) = A x**2 + B x + C                                   ',   &
         '  2:  f(x) = A x**4 + B x**3 + C x**2 + D x + E                 ',   &
         '  3:  f(x) = A cos x + B sin x + C cos 2x + D sin 2x + E        ',   &
         '  4:  x = A t**2, y = 2 A t    (parametric parabola)            ',   &
         '  5:  x = A cos (t), y = B sin (t)    (circle/ellipse)          ',   &
         '  6:  f(x) = sin x / x                                          ',   &
         '  7:  f(x) = sinh x / x                                         ',   &
         '  8:  Hermite cubic on [0,1]                                    ',   &
         '  9:  Quintic between (x,y,yp,ypp) pair                         ',   &
         ' 10:  f(x) = A exp(Bx) + C     or    A exp(B/x) + C             ',   &
         ' 11:  f(x) = A x**B exp(Cx)    or    A x**B exp(C/x)            ',   &
         ' 12:  f(x) = A x**B + C  power law                              ',   &
         ' 13:  f(x) = A (sin**n (x pi/2) - B) for x in [0, 1]            ',   &
         ' 14:  f(x) = A cosh (x/A) - A or find A for given semispan/defl.',   &
         ' 15:  f(x) = erf(x)   (error function)                          ',   &
         ' 16:  f(x) = erfc(x)  (complementary error function)            ',   &
         ' 17:  f(x) = erfcx(x) (scaled complementary error function)     ',   &
         ' 18:  f(x) = a exp(-.5((x - b)/c)**2) + d  (Gaussian function)  ',   &
         ' 19:  f(x) = normal distribution for given mu and sigma         ',   &
         ' 20:  f(x) = ? (any ideas?)                                     '/

!     Execution:

      PI = 2.E+0 * ASIN (ONE)
      DEGRAD = PI / 180.E+0
      DETAILS = ' '
      WRITE (LUNCRT, 1000) &
         ' GEN1D generates a dataset from a choice of 1-D functions.', &
         ' Results are in Tecplottable column format.'

!     Open file for output results:

      OPEN (UNIT=LUNXY, FILE='gen1d.dat', STATUS='UNKNOWN', ERR=910)

      WRITE (LUNCRT, 1001) MENU

!     Select a function:

      IFUN = 0
      DO WHILE (IFUN <= 0 .OR. IFUN >= MXMENU)
         WRITE (LUNCRT, 1001)
         CALL READI (LUNCRT, &
            'Select a function, 0 for help, or <CR> when through: ', &
            LUNKBD, IFUN, CR, EOF)
         IF (CR .OR. EOF) GO TO 999
      END DO

!     The following may be expanded to ranges of IFUN some day:

      XINDEGREES = IFUN == 3
      PARAMETRIC = IFUN == 4 .OR. IFUN == 5

      IF (PARAMETRIC) THEN ! Assume just a number of points is required as input

         NPTS = 0
         DO WHILE (NPTS <= 1  .OR.  NPTS > MXPTS)
            NPTS = 101
            CALL READI (LUNCRT, &
               'How many points do you want?  Default is 101: ', &
               LUNKBD, NPTS, CR, EOF)
         END DO

      ELSE

!        Prompt for type of grid to evaluate function on:

         YES = .TRUE.
         CALL READY (LUNCRT, &
            'Generate a grid? (Y); or read it from a file? (N)? [Y] ', &
            LUNKBD, YES, CR, EOF)

         IF (YES) THEN

            IF (XINDEGREES) THEN
               WRITE (LUNCRT, 1002) 'X is assumed to be in degrees.'
            END IF

            IOS = -1;  NPTS = 0
            DO WHILE (IOS /= 0 .OR. NPTS > MXPTS .OR. NPTS < 1)
               WRITE (LUNCRT, 1002, ADVANCE='NO') 'Enter XMIN, XMAX, DX: '
               READ (LUNKBD, *, IOSTAT=IOS) XMIN, XMAX, DX
               NPTS = 0
               IF (IOS == 0) THEN
                  NPTS = ABS ((XMAX - XMIN + 0.5E+0 * DX) / DX) + 1
                  IF (NPTS > MXPTS) WRITE (LUNCRT, 1004) &
                     ' Too many points requested. Limit: ', MXPTS
               END IF
            END DO

            DO I = 1, NPTS - 1
               X(I) = XMIN + REAL (I-1) * DX
            END DO
            X(NPTS) = XMAX

         ELSE ! Read X coordinates from a file

            IOS = -1
            DO WHILE (IOS /= 0)
               DATAFILE = 'xgrid.dat'
               CALL READC (LUNCRT, &
                  'Enter the X grid file name (<CR>="xgrid.dat"): ', &
                  LUNKBD, DATAFILE, CR, EOF)
               IF (EOF) GO TO 999

               OPEN (UNIT=LUNX, FILE=DATAFILE, STATUS='OLD', IOSTAT=IOS)
            END DO

            READ (LUNX, 1001, ERR=940)      ! Skip the title
            READ (LUNX, *, ERR=950) NPTS
            IF (NPTS < 1  .OR.  NPTS > MXPTS) GO TO 950

            DO I = 1, NPTS                  ! Ignore any column > 1
               READ (LUNX, *, ERR=960) X(I)
            END DO

         END IF

      END IF  ! End parametric/non-parametric branch

      WRITE (LUNCRT, 1001) MENU(IFUN), ' '

  500 CONTINUE

      SELECT CASE (IFUN)

         CASE (1)  ! F(X) = A X**2 + B X + C

            WRITE (LUNCRT, 1002, ADVANCE='NO') 'Enter A, B, C: '
            READ  (LUNKBD, *, ERR=500) A, B, C
            WRITE (DETAILS, '(''A, B, C: '', 1P, 3E12.4)') A, B, C

            DO I = 1, NPTS
               Y(I) = (A*X(I) + B) * X(I) + C
            END DO

         CASE (2)  ! F(X) = A X**4 + B X**3 + C X**2 + D X + E

            WRITE (LUNCRT, 1002, ADVANCE='NO') 'Enter A, B, C, D, E: '
            READ  (LUNKBD, *, ERR=500) A, B, C, D, E
            WRITE (DETAILS, '(''Coefs: '', 1P, 5E12.4)') A,B,C,D,E

            DO I = 1, NPTS
               Y(I) = (((A*X(I) + B) * X(I) + C) * X(I) + D) * X(I) + E
            END DO

         CASE (3)  ! F(X) = A COS X + B SIN X + C COS 2X + D SIN 2X + E

            WRITE (LUNCRT, 1002, ADVANCE='NO') 'Enter A, B, C, D, E: '
            READ  (LUNKBD, *, ERR=500) A, B, C, D, E
            WRITE (DETAILS, '(''A,B,C,D,E: '', 5F10.4)') A,B,C,D,E

            DO I = 1, NPTS
               XR = X(I) * DEGRAD
               Y(I) = A * COS(XR) + B * SIN(XR) + C * COS(XR+XR) + &
                      D * SIN(XR+XR) + E
            END DO

         CASE (4)  ! X = A T**2, Y = 2 A T    (parametric parabola)

            WRITE (LUNCRT, 1002, ADVANCE='NO') 'Enter A: '
            READ  (LUNKBD, *, ERR=500) A
            WRITE (DETAILS, '(''A: '', 1P, E12.4)') A

            WRITE (LUNCRT, 1002, ADVANCE='NO') 'Enter |T (max)|: '
            READ  (LUNKBD, *, ERR=500) TMAX
            DT = (TMAX + TMAX) / REAL (NPTS - 1)
            DO I = 1, (NPTS + 1) / 2
               T = TMAX - REAL (I - 1) * DT
               X(I) = A * T * T
               Y(I) = (A + A) * T
               X(NPTS + 1 - I) =  X(I)
               Y(NPTS + 1 - I) = -Y(I)
            END DO
            IF (MOD (NPTS, 2) /= 0) THEN
               I = (NPTS + 1) / 2
               X(I) = ZERO
               Y(I) = ZERO
            END IF

         CASE (5)  ! X = A COS (T), Y = B SIN (T)    (circle/ellipse)

            WRITE (LUNCRT, 1002, ADVANCE='NO') 'Enter A, B: '
            READ  (LUNKBD, *, ERR=500) A, B
            WRITE (DETAILS, '(''A, B: '', 1P, 2E12.4)') A, B

            DT = (PI + PI) / (NPTS - 1)
            DO I = 1, NPTS
               T = REAL (I - 1) * DT
               X(I) = A * COS (T)
               Y(I) = B * SIN (T)
            END DO
            X(NPTS) = A
            Y(NPTS) = ZERO

         CASE (6)  ! F(X) = SIN X / X

            DO I = 1, NPTS
               IF (X(I) == ZERO) THEN
                  Y(I) = ONE
               ELSE
                  Y(I) = SIN (X(I)) / X(I)
               END IF
            END DO

         CASE (7)  ! F(X) = SINH X / X

            DO I = 1, NPTS
               IF (X(I) == ZERO) THEN
                  Y(I) = ONE
               ELSE
                  Y(I) = SINH (X(I)) / X(I)
               END IF
            END DO

         CASE (8)  ! Hermite cubic on [0, 1]

            WRITE (LUNCRT, 1002, ADVANCE='NO') &
               'Enter Y1, Y2, Deriv1, & Deriv2: '
            READ  (LUNKBD, *, ERR=500) Y1, Y2, D1, D2
            WRITE (DETAILS, '(''Y1, Y2, D1, D2: '', 1P, 4E12.4)') &
                                Y1, Y2, D1, D2

            CALL HERMITE (Y1, Y2, D1, D2, NPTS, X, Y)

         CASE (9)  ! Quintic between (x,y,y',y") pairs

            WRITE (LUNCRT, 1002, ADVANCE='NO') 'Enter X1, Y1, YP1, YPP1: '
            READ  (LUNKBD, *, ERR=500) X1, Y1, YP1, YPP1
            WRITE (LUNCRT, 1002, ADVANCE='NO') 'Enter X2, Y2, YP2, YPP2: '
            READ  (LUNKBD, *, ERR=500) X2, Y2, YP2, YPP2
            WRITE (DETAILS, '(''YP1, Y"1, YP2, Y"2: '', 1P, 4E12.4)') &
                                YP1, YPP1, YP2, YPP2

            CALL QUINTIC (X1, Y1, YP1, YPP1, X2, Y2, YP2, YPP2, COEFS)

            DO I = 1, NPTS
              H    = X(I) - X1
              Y(I) = Y1 + H * (COEFS(1) + H * (COEFS(2) + H * &
                              (COEFS(3) + H * (COEFS(4) + H * COEFS(5)))))
            END DO

         CASE (10)  ! F(X) = A e**(BX)  +  C  or  A e**(B/X) + C

            WRITE (LUNCRT, 1002, ADVANCE='NO') 'exp(Bx) or exp(B/x)? [1|2]: '
            READ  (LUNKBD, *, ERR=500) MODE
            WRITE (LUNCRT, 1002, ADVANCE='NO') 'Enter A, B, C: '
            READ  (LUNKBD, *, ERR=500) A, B, C
            WRITE (DETAILS, '(''A, B, C: '', 1P, 3E12.4)') A, B, C

            IF (MODE == 1) THEN
               DO I = 1, NPTS
                  Y(I) = A * EXP (B*X(I)) + C
               END DO
            ELSE
               DO I = 1, NPTS
                  Y(I) = A * EXP (B/X(I)) + C
               END DO
            END IF

         CASE (11)  ! F(X) = A X**B EXP (CX)  or  A X**B EXP (C/X)

            WRITE (LUNCRT, 1002, ADVANCE='NO') 'exp(Cx) or exp(C/x)? [1|2]: '
            READ  (LUNKBD, *, ERR=500) MODE
            WRITE (LUNCRT, 1002, ADVANCE='NO') 'Enter A, B, C: '
            READ  (LUNKBD, *, ERR=500) A, B, C
            WRITE (DETAILS, '(''A, B, C: '', 1P, 3E12.4)') A, B, C

            IF (MODE == 1) THEN
               DO I = 1, NPTS
                  Y(I) = A * X(I)**B * EXP (C*X(I))
               END DO
            ELSE
               DO I = 1, NPTS
                  Y(I) = A * X(I)**B * EXP (C/X(I))
               END DO
            END IF

         CASE (12)  ! F(X) = A X**B  +  C

            WRITE (LUNCRT, 1002, ADVANCE='NO') 'Enter A, B, C: '
            READ  (LUNKBD, *, ERR=500) A, B, C
            WRITE (DETAILS, '(''A, B, C: '', 1P, 3E12.4)') A, B, C

            DO I = 1, NPTS
               Y(I) = A * X(I)**B + C
            END DO

         CASE (13)  ! F(X) = SIN**N (X PI/2) for X in [0, 1]

            WRITE (LUNCRT, 1002, ADVANCE='NO') 'Enter N, A, B: '
            READ  (LUNKBD, *, ERR=500) N, A, B
            WRITE (DETAILS, '(''N, A, B: '', I2, 2F10.6)') N, A, B

            DO I = 1, NPTS
               Y(I) = A * ((SIN (HALF*PI*X(I)))**N - B)
            END DO

         CASE (14)  ! F(X) = A COSH (X/A) - A  (Catenary)

            WRITE (LUNCRT, 1002, ADVANCE='NO') &
               'Given parameter a or given semispan+deflection? [a|s]: '
            READ  (LUNKBD, *) ANSWER
            YES = ANSWER /= 'a' .AND. ANSWER /= 'A'

            IF (YES) THEN

               IF (X(NPTS) == ZERO) THEN
                  MODE = 1
               ELSE IF (X(1) == ZERO) THEN
                  MODE = 2
               ELSE
                  MODE = 3  ! X(1) = -X(NPTS) assumed
               END IF

               WRITE (LUNCRT, 1002, ADVANCE='NO') &
                  'Enter semispan, maximum droop, and high-end y: '
               READ  (LUNKBD, *, ERR=500) B, D, YMAX

               CALL CATENARY_GRID (B, D, YMAX, MODE, NPTS, X, Y, A, IER)

               IF (IER /= 0) GO TO 999

               WRITE (DETAILS, '(''B, D, A: '', 2F10.6, 1P, E16.8)') B, D, A

            ELSE

               WRITE (LUNCRT, 1002, ADVANCE='NO') 'Enter A: '
               READ  (LUNKBD, *, ERR=500) A
               WRITE (DETAILS, '(''A: '', F10.6)') A

               DO I = 1, NPTS
                  Y(I) = A * COSH (X(I) / A) - A
               END DO

            END IF

         CASE (15)  ! F(X) = ERF (X)  (error function)

            DO I = 1, NPTS
               Y(I) = ERF (X(I))
            END DO

         CASE (16)  ! F(X) = ERFC (X)  (complementary error function)

            DO I = 1, NPTS
               Y(I) = ERFC (X(I))
            END DO

         CASE (17)  ! F(X) = ERFCX (X)  (scaled complementary error function)

            DO I = 1, NPTS
               Y(I) = ERFCX (X(I))
            END DO

         CASE (18)  ! F(X) = Gaussian function

            WRITE (LUNCRT, 1002, ADVANCE='NO') 'Enter a, b, c, d: '
            READ  (LUNKBD, *, ERR=500) A, B, C, D

            CALL GAUSSIAN (A, B, C, D, NPTS, X, Y)

         CASE (19)  ! F(X) = normal distribution

            WRITE (LUNCRT, 1002, ADVANCE='NO') 'Enter mean & std. deviation: '
            READ  (LUNKBD, *, ERR=500) MU, SIGMA

            CALL NORMAL_DISTRIBUTION (MU, SIGMA, NPTS, X, Y)

         CASE DEFAULT

!!! Can't   WRITE (LUNCRT, 1002) 'Invalid menu item - try again.'

      END SELECT

  800 CONTINUE  ! End of case statement over functions.

      YES = .FALSE.
      CALL READY (LUNCRT, &
         'Apply Gaussian type noise? (Yes/No; (<CR>=No) ', &
         LUNKBD, YES, CR, EOF)

      IF (YES) THEN

         PERCENT = ONE
         CALL READR (LUNCRT, &
            '% of data range to use as st.dev. of noise? (<CR>=1%) ', &
            LUNKBD, PERCENT, CR, EOF)
         PERCENT = PERCENT * 0.01E+0

!        Deal with Y noise first:

         XMAX = Y (1)
         XMIN = Y (1)
         DO I = 1, NPTS
            XMAX = MAX (XMAX, Y(I))
            XMIN = MIN (XMIN, Y(I))
         END DO

         SDEV = (XMAX - XMIN) * PERCENT
         IF (XMAX == XMIN) SDEV = PERCENT      ! In case Y = const.

!        Generate pseudo-normal (0, 1) deviates:

         ISEED = -123457  ! < 0 initializes a sequence

         DO I = 1, NPTS
            DNOISE = GAUSS (ISEED)
            Y(I) = Y(I) + SDEV * DNOISE ! Mean noise is zero.
         END DO

         IF (PARAMETRIC) THEN

!           Apply noise to Xs as well:

            XMAX = X(1)
            XMIN = X(1)
            DO I = 2, NPTS
               XMAX = MAX (XMAX, X(I))
               XMIN = MIN (XMIN, X(I))
            END DO

            SDEV = (XMAX - XMIN) * PERCENT

            DO I = 1, NPTS
               DNOISE = GAUSS (ISEED)
               X(I) = X(I) + SDEV * DNOISE
            END DO

         END IF
      END IF

!     Write results to the output file:

      WRITE (LUNXY, 1003) '# ', MENU(IFUN) (7:)
      WRITE (LUNXY, 1005) NPTS, DETAILS
      WRITE (LUNXY, 1016) (X(I), Y(I), I = 1, NPTS)

      WRITE (LUNCRT, 1002) '*** (X,Y) pairs are in gen1d.dat ***'
      GO TO 999

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!                           Error handling:
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  910 WRITE (LUNCRT, 1000) ' Error opening output file.  Privilege?'
      GO TO 999
  940 WRITE (LUNCRT, 1000) ' Error reading the X grid title.'
      GO TO 999
  950 WRITE (LUNCRT, 1000) ' Error in the no. of X points read:'
      WRITE (LUNCRT, 1004) ' Minimum required =', 1,     &
                           ' Maximum allowed  =', MXPTS, &
                           ' Number input     =', NPTS
      GO TO 999
  960 WRITE (LUNCRT, 1000) ' Error reading X coordinates.'
      GO TO 999

!     Formats:

 1000 FORMAT (/, (A))
 1001 FORMAT (A)
 1002 FORMAT (1X, A)
 1003 FORMAT (2A)
 1004 FORMAT (A, I4)
 1005 FORMAT ('# ', I5, ' ! ', A)
 1016 FORMAT (1X, 1P, 2E24.15)

  999 CONTINUE  ! Avoid machine-dependent STOPs
      
      END PROGRAM GEN1D
