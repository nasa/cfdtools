C+----------------------------------------------------------------------
C
      SUBROUTINE VECFIT (NXY, X, Y, LUNCRT, LUNKBD, LUNVEC, LUNLOG,
     >                   WK, MAXC, NC, C, YFIT, RMSDEV, IER)
C
C  ACRONYM: VECtor FIT: linear combination of discrete-data vectors
C           ---    ---
C  PURPOSE: VECFIT fits a linear combination of vectors FJ, J = 1:NVEC,
C           to the given vector Y, where all vectors are defined at the
C           same abscissas represented by vector X. This capability was
C           prompted by the need to determine  empirical  relationships
C           between various distributions related to airfoil geometries
C           and surface flow quantities, for inverse design methods. It
C           may find other uses.
C
C           VECFIT also evaluates the fit, for reasons described below.
C
C  METHOD:  VECFIT is interactive, so that an indefinite number of dis-
C           crete datasets can be read from disk rather than passed  as
C           arguments.  The only main arguments are the data NXY, X, Y,
C           for which a linear least squares-type fit is required, plus
C           workspace.  The overdetermined problem to be solved is:
C
C           |  F1(1)  F2(1)  F3(1) ..... FN(1)  1  | | C1 |    | Y(1) |
C           |  F1(2)  F2(2)  F3(2) ..... FN(2)  1  | | C2 |    | Y(2) |
C           |  F1(3)  F2(3)  F3(3) ..... FN(3)  1  | | C3 |    | Y(3) |
C           |    :      :      :           :    :  | |  : |    |  :   |
C           |    :      :      :           :    :  | |  : | ~= |  :   |
C           |    :      :      :           :    :  | |  : |    |  :   |
C           |  F1(N)  F2(N)  F3(N) ..... FN(N)  1  | | CN |    | Y(N) |
C           |    :      :      :           :    :  |           |  :   |
C           |    :      :      :           :    :  |           |  :   |
C           |    :      :      :           :    :  |           |  :   |
C           |  F1(M)  F2(M)  F3(M)       FN(M)  1  |           | Y(M) |
C
C           where  M = NXY;
C                  N = NC.
C
C           Suppressing the last coefficient shown (and its column of 1
C           elements) is an option.  The degenerate case of calculating
C           the best shift of F1(*) to match Y(*) is also treated.
C
C           The discrete vectors may come from one file or more.  These
C           files must have the following format:
C
C           TITLE         <for all datasets in this file>
C           M             <# points in one dataset, plus optional text>
C           X(1)   F(1)   <read as a pair, list-directed>
C           X(2)   F(2)   <any further columns are ignored>
C           X(3)   F(3)
C            :      :
C            :      :
C           X(M)   F(M)   <this may be the last record in the file>
C           M             <same M; record is included for optional text>
C           X(1)   G(1)
C           X(2)   G(2)
C            :      :     <and so on until EOF>
C
C           Each file is processed until EOF,  another is prompted for,
C           and so on, until an empty file name is detected.  Each file
C           is also echoed to the log file in the order read.
C
C           Evaluating the best fit at each X is problematical, because
C           the matrix above is destroyed during factorization  by  the
C           general purpose linear least squares solver.   The approach
C           taken is to require about twice the work-space  that  would
C           otherwise be needed - to keep the vectors around for apply-
C           ing the best coefficients after the latter are calculated.
C
C  ARGUMENTS:
C    ARG    TYPE  I/O/S   DIM    DESCRIPTION
C    NXY     I      I      -     Number of data points being fitted.
C    X       R      I     NXY    Abscissas for data points.  Must match
C                                abscissas for each discrete vector.
C    Y       R      I     NXY    Ordinates of data points.
C    LUNCRT  I      I      -     Logical unit for prompts.
C    LUNKBD  I      I      -     Logical unit for responses.
C    LUNVEC  I      I      -     Logical unit for reading vector data.
C    LUNLOG  I      I      -     Logical unit for writing log info.
C    WK      R      S    NXY *   Used for storing vectors read, and for
C                     (2*MAXC+2) setting up overdetermined system. (The
C                                data read must be kept for  evaluating
C                                the computed fit at each X (I).)
C    MAXC    I      I      -     Maximum no. of coefficients allowed.
C    NC      I      O      -     Actual no. of coefficients calculated.
C    C       R      O     NC     Computed coefs. Look at log file to
C                                interpret C(NC).
C    YFIT    R      O     NXY    YFIT(I) = SUM (J=1:NC) C(J) * FJ(X(I))
C                                where the last FJ may be 1 (constant).
C    RMSDEV  R      O      -     SQRT (Optimal sum of squares / NXY)
C    IER     I      O      -     0 means no errors;
C                                1 means failure (change data/problem).
C                                Diagnostics are written to the screen.
C
C  EXTERNAL REFERENCES:
C    HDESOL     Linear least squares by Householder decomposition.
C    LSTFIL     Echos input file to log file.
C    OPENER     File opening utility.
C    READER     Prompting utility.
C
C  ENVIRONMENT:  VAX/VMS FORTRAN;
C                IMPLICIT NONE is nonstandard.
C
C  HISTORY:
C  05/24/86  D.A.Saunders  Initial implementation, adapted from WAGFIT.
C  07/18/89    "     "     Cosmetic changes prior to Mac translation.
C  11/25/02    "     "     Added degenerate case of shifting one vector.
C
C  AUTHOR:  David Saunders, Sterling Software/NASA Ames, Palo Alto, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C ... Arguments:

      INTEGER
     >   IER, LUNCRT, LUNKBD, LUNVEC, LUNLOG, MAXC, NC, NXY
      REAL
     >   C (MAXC), RMSDEV, X (NXY), Y (NXY), YFIT (NXY),
     >   WK (NXY * (2*MAXC + 2))

C ... Local constants:

      REAL
     >   ONE, TOL, ZERO
      CHARACTER
     >   BLANK*1, OLD*3
      PARAMETER
     >  (BLANK = ' ',
     >   OLD   = 'OLD',
     >   ONE   = 1.E+0,
     >   TOL   = 1.E-05,
     >   ZERO  = 0.E+0)

C ... Local variables:

      INTEGER
     >   I, IW, J, MBYN, NPTS
      REAL
     >   DEV, SUM, XI
      LOGICAL
     >   CR, CTERM, EOF, SHIFTONLY
      CHARACTER
     >   BUFFER*80, FILNAM*50

C ... Execution:

C ... Most of the error checking must be done as files are read in:

      IER = 1
      IF (NXY .LT. 1) GO TO 910

C ... Work-space usage:
C     The first NC "columns" of length NXY in WK (*) are used to
C     read the data and keep them around for evaluating the fit.
C     (NC is not known until an empty file name is detected.)
C     The next NC columns are a repeat of the first, representing
C     the left-hand-side of the system to be solved.  The right-
C     hand-side follows, and a scratch vector follows that.

      NC = 0

      WRITE (LUNLOG, 1010) 'Data being fitted follow.'
      WRITE (LUNLOG, 1005) 'Number of points:', NXY
      WRITE (LUNLOG, 1015) (X (I), Y (I), I = 1, NXY)
      WRITE (LUNLOG, 1010) 'Discrete vectors used in fit follow. '//
     >                     '(One or more per file read.)'

C ... While FILNAM is not returned as blank ...

  200 CONTINUE
         WRITE (LUNCRT, 1005)
     >      'Number of discrete vectors set up so far:', NC

         FILNAM = BLANK
         CALL OPENER (LUNCRT, 'Enter file name for next vector(s): ',
     >                LUNKBD, FILNAM, LUNVEC, OLD)
         IF (FILNAM .EQ. BLANK) GO TO 300

C ...    Show each file in the run log for future reference:

         WRITE (LUNLOG, 1005) 'Input vector file: ' // FILNAM
         CALL LSTFIL (LUNVEC, LUNLOG, BUFFER)

C ...    Skip the title record:

         REWIND LUNVEC
         READ (LUNVEC, 1015)

C ...    Repeat until end of file ...

  250    CONTINUE

            READ (LUNVEC, *, END=200, ERR=920) NPTS
            IF (NPTS .NE. NXY) GO TO 930

C ...       Read the function values into the next column of WK (*).
C           Ignore extra columns in the input file:
            
            NC = NC + 1
            IF (NC .GT. MAXC) GO TO 940

            IW = (NC - 1) * NPTS

            DO 260 I = 1, NPTS
               IW = IW + 1
               READ (LUNVEC, *, ERR=950) XI, WK (IW)

C ...          Checking for exact match in X is not practical:

               IF (XI * X (I) .EQ. ZERO) THEN
                  IF (XI .NE. X (I)) GO TO 960
               ELSE
                  IF (ABS ((XI - X (I)) / XI) .GT. TOL) GO TO 960
               END IF

  260       CONTINUE

            GO TO 250

C ...    End repeat until EOF

C ... End while file name not blank

  300 CONTINUE

C     No more vectors.
C     Check for the degenerate case of simply shifting one vector to
C     match the given data vector best:

      IF (NC .EQ. 1) THEN

         SHIFTONLY = .FALSE.
         CALL READY (LUNCRT,
     >               'Calculate the best shift (only)? (Y/N; CR=N) ',
     >               LUNKBD, SHIFTONLY, CR, EOF)

         IF (SHIFTONLY) THEN ! Just average the difference between elements

            SUM    = ZERO
            RMSDEV = ZERO

            DO I = 1, NPTS
               DEV    = Y (I) - WK (I)
               SUM    = DEV + SUM
               RMSDEV = DEV * DEV + RMSDEV
            END DO

            C(1)   = SUM / REAL (NPTS)
            RMSDEV = SQRT (RMSDEV / REAL (NPTS))

            DO I = 1, NPTS ! Calculate the shifted Ys
               YFIT (I) = WK (I) + C (1)
            END DO

            IER = 0
            GO TO 999 ! Avoid indenting the general method

         END IF

      END IF

C ... What about a constant term?

      CTERM = .TRUE.
      CALL READY (LUNCRT,
     >   'Do you want a constant term in the fit? (Y/N; CR=Y) ',
     >   LUNKBD, CTERM, CR, EOF)

      IF (CTERM) THEN
         NC = NC + 1
         IF (NC .GT. MAXC) GO TO 940

         I = IW
         DO 320 IW = I + 1, I + NPTS
            WK (IW) = ONE
  320    CONTINUE
      END IF

      IF (NC .GT. NPTS) GO TO 970

C ... Copy these columns as the left-hand-side in the second
C     NPTS*NC segment of WK (*):

      MBYN = NPTS * NC

      DO 340 I = 1, MBYN
         WK (I + MBYN) = WK (I)
  340 CONTINUE

C ... And the right-hand-side:

      IW = 2 * MBYN
      DO 360 I = 1, NPTS
         WK (I + IW) = Y (I)
  360 CONTINUE
      IW = IW + NPTS

C ... Solve the overdetermined system by Householder factorization.
C     The required coefs. are returned starting in WK (IW+1):

      CALL HDESOL (NPTS, NPTS, NC+1, WK (MBYN+1), WK (IW+1), RMSDEV)

      DO 400 I = 1, NC
         C (I) = WK (IW + I)
  400 CONTINUE

      RMSDEV = SQRT (RMSDEV / REAL (NPTS))

C ... Finally, evaluate the fit using the saved data in the first
C     NC columns of WK (*):

      DO 440 I = 1, NPTS
         SUM = ZERO
         IW = I
         DO 420 J = 1, NC
            SUM = WK (IW) * C (J) + SUM
            IW = IW + NPTS
  420    CONTINUE
         YFIT (I) = SUM
  440 CONTINUE

      IER = 0
      GO TO 999

C ... Error handling:

  910 WRITE (LUNCRT, 1010) 'Bad number of data points. NXY =', NXY
      GO TO 999

  920 WRITE (LUNCRT, 1010) 'Error reading NPTS - quitting.'
      GO TO 999

  930 WRITE (LUNCRT, 1010) 'Mismatch: NPTS found =', NPTS,
     >                     'Number NXY expected  =', NXY
      GO TO 999

  940 WRITE (LUNCRT, 1010) 'Too many coefficients requested:  ', NC,
     >                     'Limit provided by calling program:', MAXC
      GO TO 999

  950 WRITE (LUNCRT, 1010) 'Error reading (X,Y) in dataset no.', NC
      GO TO 999

  960 WRITE (LUNCRT, 1010) 'Mismatch in abscissas. Dataset no. =', NC,
     >                     '                       Element no. =', I
      GO TO 999

  970 WRITE (LUNCRT, 1010) 'Too many coefficients requested:', NC,
     >                     'Number of data points is only:  ', NPTS
      GO TO 999

  999 RETURN

C ... FORMATs:

 1005 FORMAT (/, 1X, A, I5)
 1010 FORMAT (/, 1X, 'VECFIT:  ', A, I4)
 1015 FORMAT (1X, 1P, 2E15.6)
      END
