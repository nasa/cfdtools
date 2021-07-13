C+----------------------------------------------------------------------
C
      PROGRAM RESHAPE
C
C  PURPOSE:
C
C     RESHAPE is a utility to transform (X,Y) data in one or more of
C     a variety of ways such as scaling, rotating, reversing order, etc.
C     Since X and Y may be in specified columns other than 1 and 2,
C     RESHAPE can also serve the purpose of extracting two columns of
C     data from many without other changes.
C
C  METHOD:
C
C     > Deal with just one dataset at a time.  The formats chosen are
C       those of program SMOOTH:
C
C          [TITLE]             <Optional title, up to 80 characters>
C          [N]                 <Optional no. of pts - EOF also serves>
C                              <Blank lines are ignored>
C          X (1)  Y (1)
C           :     :            <X and Y may be extracted from specified
C           :     :             columns if there are more than two>
C           :     :
C          ! X     Y           <! suppresses points or adds comments>
C          X (N)  Y (N)
C
C     > Transformations are done in-place (Y and/or X).
C
C     > "Undo last" and "start over" operations are done with spare copies.
C
C  PROCEDURES:
C
C     ALPHA     Distinguishes text and numeric data
C     GETLINE   Gets next significant line
C     OPENER    File opening utility
C     RDLIST    Gets an indefinite number of integers
C     RDXYZ     Gets 2 or 3 columns from many
C     READER    Prompting utility (entry pts. READI, READR, etc.)
C
C  ENVIRONMENT: Fortran 90
C
C  HISTORY:
C
C     08/29/86   DAS   Initial implementation (in haste) -
C                      simple Y-translation or scaling.
C     09/03/86   DAS   Added rotation and X-translation/scaling.
C     10/21/86   DAS   Using EOF to determine number of points turns
C                      out to be inconvenient for other utilities.
C                      Expect N to be with the data.  (Changed later.)
C     08/19/88   DAS   Added option to reverse the order 1:N.
C     10/21/88   DAS   Bug fix for reverse-order option; added
C                      "start over", "switch X and Y", and
C                      "reflect" options.
C     10/13/92   DAS   Handled multi-column files with ! comments or
C                      blank lines ignored.  N and title are now optional.
C     01/03/98   DAS   Provided for rotating Y about (Zc=0, Yc).
C     05/20/99   DAS   Minor Fortran 90 changes.
C     08/30/10   DAS   The advance='no' prompting was misbehaving; added a
C                      "Done" item to the menu.
C     05/24/11   DAS   64-bit precision outputs now, not single precision.
C     07/10/21   DAS   In order to test the revised CHANGEN2D, install
C                      it as one more option here, as first done for testing
C                      CHANGEN via RESHAPE3D.  Add saving of before and
C                      after cell growth rates for this option.
C
C  AUTHOR:  David Saunders, Sterling Software/NASA Ames, Mt. View, CA.
C           Later with ELORET, Inc. and AMA, Inc. at NASA ARC.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Constants:

      INTEGER, PARAMETER ::
     >   LUNCRT = 6,
     >   LUNKBD = 5,
     >   LUNIN  = 7,
     >   LUNOUT = 8,
     >   MXMENU = 13,
     >   MXPTS  = 200000,
     >   HALFM  = (MXMENU + 4) / 2 ! 4 here allows for -2, -1, 0, and MXMENU+1

C     Variables:

      INTEGER ::
     >   CHOICE, COLUMNS (2), I, IER, IOS, J, LAST, N, NCOL, NNEW
      REAL ::
     >   ANGLE, CI, GROWTH, P, Q, SCALE, SHIFT, SI, TEMP, TOTAL, XP, YP
      REAL, DIMENSION (MXPTS) ::
     >   X, Y, XLAST, YLAST, XORIG, YORIG
      REAL, ALLOCATABLE, DIMENSION (:) ::
     >   ARC, XNEW, YNEW
      LOGICAL ::
     >   CR, EOF, FINIS, YESTITLE
      CHARACTER ::
     >   DATASET*128, MENU (-2 : MXMENU + 1)*30, METHOD*1, TITLE*80

C     Procedures:

      LOGICAL ::
     >   ALPHA
      EXTERNAL ::
     >   ALPHA, CHANGEN2D, CHORDS2D, GETLINE, OPENER, RDLIST, RDXYZ,
     >   READI, READR, READS

C     Storage:

      DATA MENU /
     >   '  -2: Start over',
     >   '  -1: Undo last transformation',
     >   '   0: Review data',
     >   '   1: Translate Y',
     >   '   2: Scale Y',
     >   '   3: Translate X',
     >   '   4: Scale X',
     >   '   5: Rotate (X,Y) about (p,q)',
     >   '   6: Reflect about the X-axis',
     >   '   7: Reflect about the Y-axis',
     >   '   8: Reverse the order (1:N)',
     >   '   9: Switch X and Y',
     >   '  10: Rotate Y about (Zc=0,Yc)',
     >   '  11: Scale X & Y the same way',
     >   '  12: Change N; same rel. sp.',
     >   '  13: Done',
     >   '                              '/ ! Last blank eases display
                                           ! of menu as two columns.
C     Execution:

      WRITE (LUNCRT, 1000)
     >   'RESHAPE applies simple transformations to an (X,Y) dataset.',
     >   'Input and output files are in [program] SMOOTH format(s).',
     >   ' '

C     :::::::::::::
C     Get the data:
C     :::::::::::::

      CALL OPENER (LUNCRT, 'Dataset name: ',
     >             LUNKBD, DATASET, LUNIN, 'OLD')

C     Look for the optional title:

  120 CALL GETLINE (LUNIN, '!', TITLE, LAST, IOS)

      IF (IOS  /= 0) GO TO 900
      IF (LAST == 0) GO TO 120

C     If we have strictly numeric data, assume there is no title:

      YESTITLE = ALPHA (TITLE)
      IF (.NOT. YESTITLE) THEN
         REWIND LUNIN         ! RDXYZ should have an optional input buffer
      END IF

C     Too awkward to avoid a prompt here (read/tokenize/backspace):

      COLUMNS (1) = 1
      COLUMNS (2) = 2
      NCOL = 2
      CALL RDLIST (LUNCRT, 'X and Y column numbers?  [1, 2]: ',
     >             LUNKBD, NCOL, COLUMNS)
      IF (NCOL < 0) GO TO 999           ! EOF = Quit

      WRITE (LUNCRT, 1001)

C     Read the data proper:

      CALL RDXYZ (2, LUNIN, LUNCRT, COLUMNS, MXPTS, N, X, Y,
     >            Y, FINIS, IER)
      IF (IER /= 0) GO TO 900  ! RDXYZ explains any errors
      IF (N == 0)   GO TO 901  ! Doesn't make sense in this application

C     Keep some copies around:

      DO I = 1, N
         XORIG (I) = X (I)
         YORIG (I) = Y (I)
         XLAST (I) = X (I)
         YLAST (I) = Y (I)
      END DO


C     :::::::::::::::::::::::::::::::::::::::::::::::::::
C     Loop over possibly several transformations per run:
C     :::::::::::::::::::::::::::::::::::::::::::::::::::

  200 CONTINUE

         WRITE (LUNCRT, 1004)                          ! Basically, I=1:MXMENU/2
     >      (MENU (I), MENU(I + HALFM), I = -2, HALFM - 3)

         IF (MOD (MXMENU, 2) == 0) WRITE (LUNCRT, 1001)

  210    CALL READI (LUNCRT, 'Pick one. EOF (^Z or ^D) means no more. ',
     >               LUNKBD, CHOICE, CR, EOF)
         IF (EOF) GO TO 800
         IF (CR)  GO TO 210

         IF (CHOICE > 0 .OR. CHOICE == -2) THEN   ! Save current values:
            DO I = 1, N
               XLAST (I) = X (I)
               YLAST (I) = Y (I)
            END DO
         END IF

         IF (CHOICE == -2) THEN      ! "Start over" from scratch:
            DO I = 1, N
               X (I) = XORIG (I)
               Y (I) = YORIG (I)
            END DO

         ELSE IF (CHOICE == -1) THEN ! "Undo" previous operation:

            DO I = 1, N
               X (I) = XLAST (I)
               Y (I) = YLAST (I)
            END DO

         ELSE IF (CHOICE == 0) THEN  ! "Review": Display data pairs.

            WRITE (LUNCRT, '(/, (I6, 2ES16.7))')
     >         (I, X(I), Y(I), I = 1, N)

         ELSE IF (CHOICE == 1) THEN  ! "Translate Y":

            CALL READR (LUNCRT, '   Enter Y shift (+ or -): ',
     >                  LUNKBD, SHIFT, CR, EOF)
            IF (CR .OR. EOF) GO TO 210

            DO I = 1, N
               Y (I) = Y (I) + SHIFT
            END DO

         ELSE IF (CHOICE == 2) THEN  ! "Scale Y":

            CALL READR (LUNCRT, '   Enter Y scale (+ or -): ',
     >                  LUNKBD, SCALE, CR, EOF)
            IF (CR .OR. EOF) GO TO 210

            DO I = 1, N
               Y (I) = Y (I) * SCALE
            END DO

         ELSE IF (CHOICE == 3) THEN  ! "Translate X":

            CALL READR (LUNCRT, '   Enter X shift (+ or -): ',
     >                  LUNKBD, SHIFT, CR, EOF)
            IF (CR .OR. EOF) GO TO 210

            DO I = 1, N
               X (I) = X (I) + SHIFT
            END DO

         ELSE IF (CHOICE == 4) THEN  ! "Scale X":

            CALL READR (LUNCRT, '   Enter X scale (+ or -): ',
     >                  LUNKBD, SCALE, CR, EOF)
            IF (CR .OR. EOF) GO TO 210

            DO I = 1, N
               X (I) = X (I) * SCALE
            END DO

         ELSE IF (CHOICE == 5) THEN  ! "Rotate (x,y) about (p,q)":

            CALL READR (LUNCRT,
     >         '   Enter rotation in degrees (+ve is anticlockwise): ',
     >                  LUNKBD, ANGLE, CR, EOF)
            IF (CR .OR. EOF) GO TO 210

  350       WRITE (LUNCRT, 1001, ADVANCE='NO')
     >         '    Enter center of rotation (Xc,Yc) as two values: '
            READ (LUNKBD, *, ERR=350, END=210) P, Q

            ANGLE = ANGLE * ASIN (1.E+0) / 90.E+0
            CI = COS (ANGLE)
            SI = SIN (ANGLE)

            DO I = 1, N
               XP = X (I) - P
               YP = Y (I) - Q
               X (I) = XP * CI - YP * SI + P
               Y (I) = XP * SI + YP * CI + Q
            END DO

         ELSE IF (CHOICE == 6) THEN  ! "Reflect Y about the X-axis":

            DO I = 1, N
               Y (I) = -Y (I)
            END DO

         ELSE IF (CHOICE == 7) THEN  ! "Reflect X about the Y-axis":

            DO I = 1, N
               X (I) = -X (I)
            END DO

         ELSE IF (CHOICE == 8) THEN  ! "Reverse the order":

            DO I = 1, (N + 1) / 2
               J = N + 1 - I
               TEMP  = X (I)
               X (I) = X (J)
               X (J) = TEMP
               TEMP  = Y (I)
               Y (I) = Y (J)
               Y (J) = TEMP
            END DO

         ELSE IF (CHOICE == 9) THEN  ! "Switch X and Y":

            DO I = 1, N
               TEMP = X (I)
               X (I) = Y (I)
               Y (I) = TEMP
            END DO

         ELSE IF (CHOICE == 10) THEN  ! "Rotate Y about (Z=0,Yc)":

            CALL READR (LUNCRT,
     >         '   Enter rotation in degrees (+ve is anticlockwise): ',
     >                  LUNKBD, ANGLE, CR, EOF)
            IF (CR .OR. EOF) GO TO 210

            CALL READR (LUNCRT,
     >         '   Enter Yc for center of rotation (Zc=0,Yc): ',
     >                  LUNKBD, Q, CR, EOF)
            IF (CR .OR. EOF) GO TO 210

            ANGLE = ANGLE * ASIN (1.E+0) / 90.E+0
            CI = COS (ANGLE)

            DO I = 1, N
               Y (I) = (Y (I) - Q) * CI + Q
            END DO

         ELSE IF (CHOICE == 11) THEN  ! "Scale X & Y the same way":

            CALL READR (LUNCRT, '   Enter scale (+ or -): ',
     >                  LUNKBD, SCALE, CR, EOF)
            IF (CR .OR. EOF) GO TO 210

            DO I = 1, N
               X (I) = X (I) * SCALE
               Y (I) = Y (I) * SCALE
            END DO

         ELSE IF (CHOICE == 12) THEN  ! "Change N; same relative spacing":

            WRITE (LUNCRT, '(A, I7)') '   Current number of points: ', N
            CALL READI (LUNCRT,        '  Desired number of points: ',
     >                  LUNKBD, NNEW, CR, EOF)
            IF (CR .OR. EOF) GO TO 210

            WRITE (LUNCRT, '(A)') '   Interpolation methods:',
     >                         '   M (monotonic), B (loose), L (Linear)'
            CALL READC (LUNCRT, '  Interpolation choice: ',
     >                  LUNKBD, METHOD, CR, EOF)
            IF (CR .OR. EOF) GO TO 210

C           Save the current growth rates:

            OPEN (LUNOUT, FILE='growth-rates-before.dat',
     >            STATUS='UNKNOWN')

            ALLOCATE (ARC(N))

            CALL CHORDS2D (N, X, Y, .FALSE., TOTAL, ARC)

            DO I = 2, N-1
               GROWTH = (ARC(I+1) - ARC(I)) / (ARC(I) - ARC(I-1))
               WRITE (LUNOUT, '(2ES16.8)') ARC(I), GROWTH
            END DO
            CLOSE (LUNOUT)

            DEALLOCATE (ARC)
            ALLOCATE (XNEW(NNEW), YNEW(NNEW))

            CALL CHANGEN2D (1, N, X, Y, 1, NNEW, XNEW, YNEW, METHOD)

            N = NNEW
            X(1:N) = XNEW(:)
            Y(1:N) = YNEW(:)

            DEALLOCATE (XNEW, YNEW)

            OPEN (LUNOUT, FILE='growth-rates-after.dat',
     >            STATUS='UNKNOWN')

            ALLOCATE (ARC(N))

            CALL CHORDS2D (N, X, Y, .FALSE., TOTAL, ARC)

            DO I = 2, N-1
               GROWTH = (ARC(I+1) - ARC(I)) / (ARC(I) - ARC(I-1))
               WRITE (LUNOUT, '(2ES16.8)') ARC(I), GROWTH
            END DO
            CLOSE (LUNOUT)
            DEALLOCATE (ARC)

         ELSE IF (CHOICE == 13) THEN  ! Done (works better than ^D)

            GO TO 800
         ELSE
            GO TO 210
         END IF

      GO TO 200

C     :::::::::::::
C     Save results:
C     :::::::::::::

  800 WRITE (LUNCRT, 1001)

      CALL OPENER (LUNCRT, 'Output file name?  EOF = quit: ',
     >   LUNKBD, DATASET, LUNOUT, 'NEW')

      WRITE (LUNCRT, 1001)

      IF (YESTITLE) THEN
         CALL READS (LUNCRT, 'Output title line? <CR>=same: ',
     >               LUNKBD, TITLE, CR, EOF)

         I = LEN_TRIM (TITLE)
         WRITE (LUNOUT, 1001, ERR=902) TITLE(1:I)
         WRITE (LUNOUT, 1003, ERR=902) N
      END IF

      WRITE (LUNOUT, '(2ES24.15)', ERR=902) (X(I), Y(I), I = 1, N)

      GO TO 999


C     Error handling:

  900 WRITE (LUNCRT, 1000) 'RESHAPE: Error reading the data.'
      GO TO 999
  901 WRITE (LUNCRT, 1000) 'RESHAPE: No data points found.'
      GO TO 999
  902 WRITE (LUNCRT, 1000) 'RESHAPE: Error saving results.'

  999 CONTINUE  ! Avoid system-dependencies with STOP


C     Formats:

 1000 FORMAT (/, (1X, A))
 1001 FORMAT (A)
 1002 FORMAT (A, I3)
 1003 FORMAT (I5)
 1004 FORMAT (/, (A, 10X, A))

      END PROGRAM RESHAPE
