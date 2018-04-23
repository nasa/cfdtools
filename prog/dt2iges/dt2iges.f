C+------------------------------------------------------------------------------
C
      PROGRAM DT2IGES
C
C     Purpose:
C
C           DT2IGES converts one file containing one NURBS curve or surface
C        between the David Taylor Research Center (DTRC) C-array format and
C        the Initial Graphics Exchange Specification (IGES) format (either
C        direction). 
C
C     Method:
C
C           The direction is determined from the first line of the input file:
C        if columns 73:80 contain S0000001 it is assumed to be an IGES file.
C
C     Notes:
C
C           The decision to process just one B-spline in a single file per run
C        was influenced by the DT routines DTRDCA, DTWRCA (now merged as DTRDWR)
C        which originally included opening and closing of the DT file by name.
C        This was incompatible with the < in > out paradigm, so DTRDWR now
C        assumes the file has been opened by the application.  But it still does
C        not deal properly with more than one spline per file, at least for
C        reading, since it reads to EOF unless the length is known ahead of
C        time.
C
C           In the curve case, DTWRIG always writes a 3-space rational B-spline
C        representation (with all Z coordinates 0. and all weights 1. if the
C        curve is really 2-space and nonrational as normally output by
C        BSPROFILE).  Therefore, this program checks for reading such an IGES
C        file and converts it to the 2-space/nonrational DT representation if
C        possible, for BSPROFILE purposes.
C
C           The following is an indication of the size of a DT C-array:
C
C              Curve:   # control pts. * NDIM + # knots + 5
C                       where NDIM = 2 for non-rational (X, Y) curve;
C                                  = 3 for rational (X, Y) curve, etc.
C
C              Surface: # U control pts. * # V control pts. * NDIM +
C                       # U knots + # V knots + 8
C
C           Conservatively large array dimensions are used here - let the
C        error handling in the subroutines catch any over-runs.
C
C     Procedures:
C
C        DTRDIG   DT_NURBS utility for converting the CHARACTER image of
C                 of an IGES-format NURBS curve or surface to C-array form.
C
C        DTRDWR   Supplement to the DT_NURBS library for reading or writing
C                 a file containing one C-array form of a curve or surface.
C
C        DTWRIG   Analog of DTRDIG (C-array to IGES character array).
C
C     Environment:
C
C        SGI IRIS/IRIX, DEC AXP/OpenVMS, FORTRAN 77 + minor extensions
C
C     History:
C
C        04/28/92   D.Saunders,  Initial implementation following the
C                   Sterling     decision not to encumber design codes
C                   Software/    with DTRDIG and DTWRIG.  (DTRDWR is
C                   NASA Ames    much more compact and should suffice.)
C
C        05/08/92   D.Saunders   Added a title to the DTRDWR file; inserted
C                                it into and extracted it from the IGES file.
C
C        07/27/97       "        Gave up on SYS$ERROR - IRIX f90 won't compile
C                                it.  This meant abandoning the < in > out
C                                paradigm.  Kept it simple by using READER.
C
C        07/28/97       "        If a rational 3-space IGES curve is really a
C                                2-space nonrational curve, the DT form is
C                                written accordingly for BSPROFILE purposes.
C                                Likewise, a 2-space nonrational DT curve is
C                                converted to 3-space/rational IGES form.
C
C-------------------------------------------------------------------------------

      IMPLICIT NONE

C     Local constants:

      INTEGER
     >   LUNCRT, LUNKBD, LUNIN, LUNOUT, MAXC, MAXIG
      DOUBLE PRECISION
     >   ONE, ZERO
      CHARACTER * 8
     >   TITLETAG
      PARAMETER
     >  (LUNCRT   = 6,          ! Screen
     >   LUNKBD   = 5,          ! Keyboard
     >   LUNIN    = 1,          ! Input file
     >   LUNOUT   = 2,          ! Output file
     >   MAXC     = 20000,      ! Limit on the length of a C array
     >   MAXIG    = 2000,       ! Limit on the number of lines in the IGES file
     >   ONE      = 1.D+0,
     >   TITLETAG = 'S0000001', ! IGES code for a title (not optional here)
     >   ZERO     = 0.D+0)

C     Local variables:

      INTEGER
     >   I, IER, IX, IY, IZ, IW, NC, NCPTS, NDEG, NDEP, NINDEP, NLINES
      DOUBLE PRECISION
     >   C (MAXC), Z
      LOGICAL
     >   DEFAULT, QUIT, RATIONAL, ZCONSTANT
      CHARACTER
     >   FILENAME * 48, IGES (MAXIG) * 80, TITLE * 80

C     Procedures:

      EXTERNAL
     >   DTRDIG, DTRDWR, DTWRIG, READS

C     Execution:

  100 CALL READS (LUNCRT, 'Input file?  ', LUNKBD, FILENAME, DEFAULT,
     >            QUIT)
      IF (QUIT) GO TO 999
      IF (DEFAULT) GO TO 100

      OPEN (UNIT=LUNIN, FILE=FILENAME, STATUS='OLD', ERR=100)

  150 CALL READS (LUNCRT, 'Output file? ', LUNKBD, FILENAME, DEFAULT,
     >            QUIT)
      IF (QUIT) GO TO 999
      IF (DEFAULT) GO TO 150

      OPEN (UNIT=LUNOUT, FILE=FILENAME, STATUS='UNKNOWN', ERR=150)

C     Read the first line of the input file as character data.

      READ (LUNIN, '(A)') TITLE
      REWIND LUNIN

      IF (TITLE (73 : 80) .EQ. TITLETAG) THEN

C        Treat the input file as an IGES file.
C        Read all lines into memory (one B-spline expected).

         READ (LUNIN, '(A)', END=200) IGES

C        Convert these lines to a C-array:

  200    CALL DTRDIG (IGES, MAXC, C, NC, IER)

         IF (IER .NE. 0) THEN
            WRITE (LUNCRT, '(/, A, I6)')
     >         ' *** Bad return from DTRDIG.  IER: ', IER,
     >         '     Length of C-array needed if IER = -3: ', NC
            GO TO 999
         END IF

C        Do we a curve that is really 2-space/nonrational?

         NINDEP = NINT (C (1))  ! # independent variables
         NDEP   = NINT (C (2))  ! # dependent variables

         IF (NINDEP .EQ. 1 .AND. NDEP .EQ. -4) THEN  ! 3-space, rational curve

            NDEG  = NINT (C (3)) - 1
            NCPTS = NINT (C (4))
            IX = 6 + NDEG + NCPTS  ! Offset of X coordinates of control pts.
            IY = IX + NCPTS
            IZ = IY + NCPTS
            IW = IZ + NCPTS
            Z  = C (1 + IZ)
            ZCONSTANT = .TRUE.
            RATIONAL  = .FALSE.

            DO I = 1, NCPTS
               IF (C (I + IZ) .NE. Z) ZCONSTANT = .FALSE.
               IF (C (I + IW) .NE. ONE) RATIONAL = .TRUE.
            END DO

            IF (ZCONSTANT .AND. .NOT. RATIONAL) THEN
               C (2) = 2
               NC = NC - NCPTS - NCPTS

               WRITE (LUNCRT, '(A, E15.6)')
     >            ' Curve is 2-space & nonrational.  Constant Z: ', Z
            END IF

         END IF

C        Titles saved by DTWRIG are worthless:

         WRITE (LUNCRT, '(1X, A)') 'Title found:', TITLE (1 : 72)
         CALL READS (LUNCRT, 'Preferred title? (<CR> = leave it):',
     >               LUNKBD, TITLE (1 : 72), DEFAULT, QUIT)
         IF (QUIT) GO TO 999

C        Write the curve or surface in DT format:

         CALL DTRDWR ('W', LUNOUT, TITLE (1 : 72), NC, C, NC, IER)

         IF (IER .NE. 0) THEN
            WRITE (LUNCRT, '(/, A, I6)')
     >         ' *** Bad return from DTRDWR (writing).  IER: ', IER
            GO TO 999
         END IF

      ELSE

C        Read a titled C-array from the input file:

         CALL DTRDWR ('R', LUNIN, TITLE, MAXC, C, NC, IER)

         IF (IER .NE. 0) THEN
            WRITE (LUNCRT, '(/, A, I6)')
     >         ' *** Bad return from DTRDWR (reading).  IER: ', IER
            GO TO 999
         END IF

C        Do we have a curve that is 2-space/nonrational?

         NINDEP = NINT (C (1))  ! # independent variables
         NDEP   = NINT (C (2))  ! # dependent variables

         IF (NINDEP .EQ. 1 .AND. NDEP .EQ. 2) THEN

            Z = ZERO
            CALL READD (LUNCRT,
     >         'Z coordinate for this curve? (<CR> = 0.>) ',
     >         LUNKBD, Z, DEFAULT, QUIT)
            IF (QUIT) GO TO 999

            C (2) = -4  ! 3-space, rational
            NDEG  = NINT (C (3)) - 1
            NCPTS = NINT (C (4))

            NC = NC + NCPTS + NCPTS

            IF (NC .GT. MAXC) THEN
               WRITE (LUNCRT, '(A, 2I6)')
     >            ' No room to add Zs & weights. NC & MAXC: ', NC, MAXC
               GO TO 999
            END IF

            IX = 6 + NDEG + NCPTS  ! Offset of X coordinates of control pts.
            IY = IX + NCPTS
            IZ = IY + NCPTS
            IW = IZ + NCPTS

            DO I = 1, NCPTS
               C (I + IZ) = Z
               C (I + IW) = ONE
            END DO

         END IF

         IGES (1) = TITLE (1 : 72) // TITLETAG

C        Convert the C-array to a character array image of an IGES file:

         CALL DTWRIG (C, MAXIG, IGES, NLINES, IER)

         IF (IER .NE. 0) THEN
            WRITE (LUNCRT, '(/, A, I6)')
     >         ' *** Bad return from DTWRIG.  IER: ', IER,
     >         '     Length of IGES array needed if IER = -3: ', NLINES
            GO TO 999
         END IF

C        Write the IGES file:

         WRITE (LUNOUT, '(A)') (IGES (I), I = 1, NLINES)

      END IF

  999 CONTINUE  ! Avoid STOP problems

      END
