C+----------------------------------------------------------------------
C
      PROGRAM SORTER
C
C  PURPOSE:
C
C        SORTER is a simple driver for testing or applying the available
C     sorting algorithms.  It is NOT intended to rival powerful utilities
C     such as VMS's SORT.  It was prompted by the need to investigate an
C     apparent error in one of the heap-sort utilities.  Incorporation of
C     further algorithms as they come along will be straightforward.
C
C        This version treats numerical data only (not character data).
C
C  METHOD:
C
C     > One or more columns of data (integer or real) are read from a
C       formatted file.  The lines are then sorted by the specified column
C       and output in the specified order (ascending or descending).
C
C     > The data format is compatible with other utilities ("SMOOTH" format):
C
C          TITLE
C          X (1)  Y (1)  Z (1)  ...  <Any number of columns>
C          X (2)  Y (2)  Z (2)  ...
C           :     :
C           :     :                  <Read as text tokens; converted to reals>
C           :     :
C          X (N)  Y (N)  Z (N)  ...  <EOF signals end of data>
C
C  PROCEDURES:
C
C     COUNTR   Used to determine how many columns are in the data file
C     HSORTRI  Heap sort utility
C     OPENER   File opening utility
C     RDXYZ    Convenient for picking out up to 3 columns from many
C     READER   Prompting utility
C     XSORT    Bubble sort utility
C
C  ENVIRONMENT:
C
C     VAX/VMS, FORTRAN 77, with the following exceptions:
C
C     > IMPLICIT NONE
C     > Variable names up to 8 characters
C     > Trailing ! comments
C
C  HISTORY:
C
C     12/18/91  D.A.Saunders  Initial implementation to investigate a
C                             reported failure in HSORTRI.  Used RDXYZ
C                             (up to 3 columns) to start with - only
C                             generalize it if the need arises.
C     06/28/06    "     "     HSORTRI overwrites its sorting column (RLIST
C                             argument), while other utilities do not.
C                             This was not properly handled by POINTERS.
C
C  AUTHOR:  David Saunders, Sterling Software/NASA Ames, Mt. View, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Local constants:

      INTEGER
     >   LUNCRT, LUNKBD, LUNIN, LUNOUT, MXCOLS, MXMENU, MXPTS
      CHARACTER
     >   BLANK * 1
      PARAMETER
     >  (BLANK  = ' ',
     >   LUNCRT = 6,
     >   LUNKBD = 5,
     >   LUNIN  = 7,
     >   LUNOUT = 8,
     >   MXCOLS = 3,
     >   MXMENU = 2,
     >   MXPTS  = 5000)

C     Local variables:

      INTEGER
     >   CHOICE, COLUMNS (3), I, I1, I2, IER, ILIST (MXPTS), INC,
     >   J, K, N, NUMTOK, SORTCOL
      REAL
     >   ROW (MXCOLS), X (MXPTS, MXCOLS)
      LOGICAL
     >   ASCENDING, DEFAULT, EOF, FINIS, INPLACE, OK, POINTERS, QUIT
      CHARACTER
     >   BUFFER * 80, DATASET * 60, MENU (MXMENU) * 48, TITLE * 80

C     Procedures:

      EXTERNAL
     >   COUNTR, HSORTRI, OPENER, RDXYZ, READI, READS, READY, XSORT

C     Data:

      DATA
     >   MENU / '   1. HSORTRI: Heap sort of Reals via Indices   ',
     >          '   2. XSORT:   Bubble sort of Reals via Indices '/

C     Execution:
C     ----------

      WRITE (LUNCRT, '(/, (A))')
     >   ' SORTER drives available sorting algorithms.',
     >   ' Input and output files are in [program] SMOOTH format',
     >   ' (up to 3 columns of data in this version).',
     >   BLANK

C     Get the data:

      CALL OPENER (LUNCRT, 'Enter dataset name: ',
     >   LUNKBD, DATASET, LUNIN, 'OLD')

      READ (LUNIN, 1001) TITLE

C     Count the entries on the first data line, then backspace:

      CALL COUNTR (LUNIN, BUFFER, NUMTOK, OK, EOF)
      IF (EOF .OR. .NOT. OK) GO TO 901
      IF (NUMTOK .GT. MXCOLS) GO TO 902

      DO 150, I = 1, NUMTOK   ! RDXYZ allows up to 3 from 20, but assume
         COLUMNS (I) = I      ! the first 1, 2, or 3 for now
  150 CONTINUE

  160 SORTCOL = 1
      IF (NUMTOK .GT. 1) THEN
         CALL READI (LUNCRT, 'Sort field number?  (<CR>=1): ',
     >               LUNKBD, SORTCOL, DEFAULT, QUIT)
         IF (QUIT) GO TO 999
         IF (SORTCOL .LT. 1 .OR. SORTCOL .GT. NUMTOK) GO TO 160
      END IF

      ASCENDING = .TRUE.
      CALL READY (LUNCRT,
     >   'Sort in ascending order?  (<CR>=Yes; No=descending): ',
     >   LUNKBD, ASCENDING, DEFAULT, QUIT)
      IF (QUIT) GO TO 999

C     Read the data into memory for outputting in the specified order.

      CALL RDXYZ (NUMTOK, LUNIN, LUNCRT, COLUMNS, MXPTS, N,
     >            X (1, 1), X (1, 2), X (1, 3), FINIS, IER)
      IF (IER .NE. 0) GO TO 999     ! RDXYZ gives its own diagnostics
      IF (N .EQ. 0) GO TO 999

C     Initialize the pointers (done internally by some methods):

      DO 180, I = 1, N
         ILIST (I) = I
  180 CONTINUE


C     Select a sorting method (one method per run):
C     ---------------------------------------------

      WRITE (LUNCRT, '(/, (A))') MENU, BLANK

  210 CHOICE = 1
      CALL READI (LUNCRT, 'Pick a method. <CR=1> ',
     >            LUNKBD, CHOICE, DEFAULT, QUIT)
      IF (QUIT) GO TO 999
      IF (CHOICE .LT. 1 .OR. CHOICE .GT. MXMENU) GO TO 210


      IF (CHOICE .EQ. 1) THEN                   ! Heap sort

         CALL HSORTRI (X (1, SORTCOL), ILIST, N)
         POINTERS = .TRUE.                      ! HSORTRI reorders RLIST
         INPLACE  = .TRUE.                      ! in-place

      ELSE IF (CHOICE .EQ. 2) THEN              ! Bubble sort

         CALL XSORT (N, X (1, SORTCOL), ILIST)
         POINTERS = .TRUE.
         INPLACE  = .FALSE.

      END IF


C     Save results:

  800 CALL OPENER (LUNCRT,
     >   'Enter output file name.  <CR>=same as input; EOF=quit. ',
     >   LUNKBD, DATASET, LUNOUT, 'UNKNOWN')

      CALL READS (LUNCRT, 'Enter output title line. <CR>=same. ',
     >   LUNKBD, TITLE, DEFAULT, QUIT)

      WRITE (LUNOUT, 1001, ERR=903) TITLE

      IF (ASCENDING) THEN
         I1  = 1
         I2  = N
         INC = 1
      ELSE
         I1  = N
         I2  = 1
         INC = -1
      END IF

      DO I = I1, I2, INC
         IF (POINTERS) THEN
            J = ILIST (I)
         ELSE
            J = I
         END IF
         DO K = 1, NUMTOK
            ROW (K) = X (J, K)
         END DO
         IF (INPLACE) ROW (SORTCOL) = X (I, SORTCOL)
         WRITE (LUNOUT, *, ERR=903) (ROW (K), K = 1, NUMTOK)
      END DO

      GO TO 999


C     Error handling:

  901 STOP 'Error counting fields on first data line.'
  902 STOP 'Cannot handle this many fields per line yet.'
  903 STOP 'Error saving results.'

  999 STOP ' '

C     Formats:

 1001 FORMAT (A)
 1002 FORMAT (A, I3)
 1003 FORMAT (I5)

      END
