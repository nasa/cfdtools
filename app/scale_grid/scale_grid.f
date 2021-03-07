C*******************************************************************************
C
      PROGRAM SCALE_GRID
C
C     Description:
C
C        SCALE_GRID applies one or more transformations (shift, scale, rotate,
C     etc.) to all blocks of the indicated grid, which may be in any of the
C     standard formats.
C
C     Procedures:
C
C        CFD_IO_PACKAGE   I/O utilities for CFD files (Mark Rimlinger)
C        READER           Prompting utility
C        ROTATE2D         2-space rotation utility
C
C     History:
C
C        02/18/00  D. Saunders  Initial adaptation of RESHAPE3D, which
C                               operates on 3-column formatted datasets.
C        09/09/06  "      "     Added explicit inches <-> meters options.
C
C     Author:  David Saunders, Raytheon/NASA Ames, Moffett Field, CA.
C
C*******************************************************************************

C     Modules:

      USE CFD_IO_PACKAGE

C     Constants:

      IMPLICIT NONE

      INTEGER, PARAMETER ::
     >   LUNKBD = 5,
     >   LUNCRT = 6,
     >   MXMENU = 20,
     >   HALFM  = (MXMENU + 4) / 2 ! 4 here allows for -2, -1, 0, and MXMENU + 1

      REAL, PARAMETER ::
     >   INCHES_TO_METERS = 0.0254,
     >   METERS_TO_INCHES = 1.0 / INCHES_TO_METERS

C     Variables:

      INTEGER
     >   CHOICE, I, I1, I2, IBLOCK, IOS, N, NBLOCK

      LOGICAL
     >   CR, EOF

      CHARACTER
     >   ANSWER * 1, MENU (-2 : MXMENU + 1) * 35

      REAL, POINTER ::
     >   X(:,:)

      INTEGER, POINTER ::
     >   IJK(:,:), IP(:)

C     Composite data types:

      TYPE (IOMESHFILETYPE), POINTER ::
     >   INFILE, OUTFILE

C     Data:

      DATA MENU /
     >   '  -2: Start over',
     >   '  -1: Undo last transformation',
     >   '   0: Review data',
     >   '   1: Translate X',
     >   '   2: Scale X',
     >   '   3: Translate Y',
     >   '   4: Scale Y',
     >   '   5: Translate Z',
     >   '   6: Scale Z',
     >   '   7: Reflect Z about the XY-plane',
     >   '   8: Reflect X about the YZ-plane',
     >   '   9: Reflect Y about the ZX-plane',
     >   '  10: Reverse the order (1:N)',
     >   '  11: Switch X & Y',
     >   '  12: Switch Y & Z',
     >   '  13: Switch Z & X',
     >   '  14: Scale X & Y & Z the same way',
     >   '  15: Rotate (X,Y) about (Xc,Yc)',
     >   '  16: Rotate (Y,Z) about (Yc,Zc)',
     >   '  17: Rotate (Z,X) about (Zc,Xc)',
     >   '  18: Convert inches to meters',
     >   '  19: Convert meters to inches',
     >   '  99: Done',
     >   '                                   '/ ! Last ' ' eases display
                                                ! of menu as two columns.

C     Execution:

      WRITE (LUNCRT, '(A)') ' ',
     >   ' SCALE_GRID applies simple transformations to an XYZ grid.',
     >   ' Be sure to specify contiguous Xs, contiguous Ys, etc., for',
     >   ' the internal storage mode (2).'

C     Read the input grid.  File attributes will be prompted for:

      CALL IO_R_STR_X_P3D (X, NBLOCK, IJK, IP, -1, INFILE)

C     Display the block dimensions to reassure the user:

      WRITE (LUNCRT, '(A)')
      WRITE (LUNCRT, '(I6, 2X, 3I5)')
     >    (IBLOCK, IJK(1:3,IBLOCK), IBLOCK = 1, NBLOCK)

C     Loop over possibly several transformations per run:
C     ---------------------------------------------------

  200 CONTINUE

         WRITE (LUNCRT, '(/, (A, 10X, A))')            ! Basically, I=1:MXMENU/2
     >      (MENU (I), MENU(I + HALFM), I = -2, HALFM - 3)

         IF (MOD (MXMENU, 2) == 0) WRITE (LUNCRT, '(A)')

  210    CALL READI (LUNCRT, 'Make a choice. EOF (^D) means no more. ',
     >      LUNKBD, CHOICE, CR, EOF)

         IF (CHOICE == 99) GO TO 800
         IF (EOF) GO TO 800
         IF (CR)  GO TO 210
         IF (CHOICE < -2 .OR. CHOICE >= MXMENU) GO TO 210

         IF (CHOICE > 0 .OR. CHOICE == -2) THEN   ! Save current values:

C ***       XLAST = X ! Leave hooks in case the memory usage is worth it

         END IF

         IF (CHOICE == -2) THEN      ! "Start over" from scratch:

C ***       X = XORIG

            WRITE (LUNCRT, '(/, A)')
     >         ' Starting over is not implemented.  Sorry.'

         ELSE IF (CHOICE == -1) THEN ! "Undo" previous operation:

C ***       X = XLAST
        
            WRITE (LUNCRT, '(/, A)')
     >         ' The undo operation is not implemented.  Sorry.'

         ELSE IF (CHOICE == 0) THEN  ! "Review": Display the data.

            WRITE (LUNCRT, '(/, A, /)') ' First 10 points of block 1:'
            WRITE (LUNCRT, '(1X, 1P, 3E16.7)')
     >         (X(I,1), X(I,2), X(I,3), I = 1, 10)

         ELSE

C           Process grid blocks one at a time:

            DO IBLOCK = 1, NBLOCK

               I1 = IP(IBLOCK)        ! First word of this block
               I2 = IP(IBLOCK+1) - 1
               N  = I2 - I1 + 1       ! # pts in block

               CR  = .FALSE. ! TRANSFORM may want to return to the menu
               EOF = .FALSE.

C              Internal procedure similar to RESHAPE3D:

               CALL TRANSFORM (IBLOCK, N, X(I1:I2,1), X(I1:I2,2),
     >                         X(I1:I2,3))

               IF (CR .OR. EOF) EXIT

            END DO ! Next block

         END IF

CCC   GO TO 200
      GO TO 210  ! Repeating the menu is rather redundant

  800 CONTINUE

C     Save transformed grid?

      CALL READC (LUNCRT,
     >   'Save results?  Y or <CR>=Yes; ^D = No): ',
     >   LUNKBD, ANSWER, CR, EOF)

      IF (.NOT. EOF) THEN

         IF (.NOT. CR) THEN

            CR = ANSWER == 'Y' .OR. ANSWER == 'y'

            IF (.NOT. CR) THEN

               IF (ANSWER /= 'N' .AND. ANSWER /= 'n') GO TO 800

            END IF

         END IF

         IF (CR) THEN

             CALL IO_W_STR_X_P3D (X, NBLOCK, IJK, IP, -1, OUTFILE)

         END IF

      END IF
 
C *** STOP ! Avoid system dependencies.


C     SCALE_GRID internal procedure:

      CONTAINS

!        ------------------------------------------
         SUBROUTINE TRANSFORM (IBLOCK, N, X, Y, Z)
!        ------------------------------------------

!        The X, Y, Z arguments simplify reuse of code from RESHAPE3D.
!        Suppress local prompts if IBLOCK is not 1.
!        CHOICE, CR, EOF, and the LUNs are inherited from the calling program.

!        Arguments:

         INTEGER
     >      IBLOCK, N
         REAL, DIMENSION (N) ::
     >      X, Y, Z

!        Local variables:

         INTEGER
     >      I, J
         REAL
     >      TEMP
         REAL, SAVE ::
     >      ANGLE, P, Q, SCALE, SHIFT
         LOGICAL
     >      PROMPT

!        Execution:

         PROMPT = IBLOCK == 1

         IF (CHOICE == 1) THEN  ! "Translate X":

            IF (PROMPT) THEN
               CALL READR (LUNCRT, '   Enter X shift (+ or -): ',
     >            LUNKBD, SHIFT, CR, EOF)
               IF (CR .OR. EOF) GO TO 999
            END IF

            DO I = 1, N
               X (I) = X (I) + SHIFT
            END DO

         ELSE IF (CHOICE == 2) THEN  ! "Scale X":

            IF (PROMPT) THEN
               CALL READR (LUNCRT, '   Enter X scale (+ or -): ',
     >            LUNKBD, SCALE, CR, EOF)
               IF (CR .OR. EOF) GO TO 999
            END IF

            DO I = 1, N
               X (I) = X (I) * SCALE
            END DO

         ELSE IF (CHOICE == 3) THEN  ! "Translate Y":

            IF (PROMPT) THEN
               CALL READR (LUNCRT, '   Enter Y shift (+ or -): ',
     >            LUNKBD, SHIFT, CR, EOF)
               IF (CR .OR. EOF) GO TO 999
            END IF

            DO I = 1, N
               Y (I) = Y (I) + SHIFT
            END DO

         ELSE IF (CHOICE == 4) THEN  ! "Scale Y":

            IF (PROMPT) THEN
               CALL READR (LUNCRT, '   Enter Y scale (+ or -): ',
     >            LUNKBD, SCALE, CR, EOF)
               IF (CR .OR. EOF) GO TO 999
            END IF

            DO I = 1, N
               Y (I) = Y (I) * SCALE
            END DO

         ELSE IF (CHOICE == 5) THEN  ! "Translate Z"

            IF (PROMPT) THEN
               CALL READR (LUNCRT, '   Enter Z shift (+ or -): ',
     >            LUNKBD, SHIFT, CR, EOF)
               IF (CR .OR. EOF) GO TO 999
            END IF

            DO I = 1, N
               Z (I) = Z (I) + SHIFT
            END DO

         ELSE IF (CHOICE == 6) THEN  ! "Scale Z"

            IF (PROMPT) THEN
               CALL READR (LUNCRT, '   Enter Z scale (+ or -): ',
     >            LUNKBD, SCALE, CR, EOF)
               IF (CR .OR. EOF) GO TO 999
            END IF

            DO I = 1, N
               Z (I) = Z (I) * SCALE
            END DO

         ELSE IF (CHOICE == 7) THEN  ! "Reflect Z about the XY-plane":

            DO I = 1, N
               Z (I) = -Z (I)
            END DO

         ELSE IF (CHOICE == 8) THEN  ! "Reflect X about the YZ-plane":

            DO I = 1, N
               X (I) = -X (I)
            END DO

         ELSE IF (CHOICE == 9) THEN  ! "Reflect Y about the ZX-plane":

            DO I = 1, N
               Y (I) = -Y (I)
            END DO

         ELSE IF (CHOICE == 10) THEN  ! "Reverse the order":

            DO I = 1, (N + 1) / 2
               J = N + 1 - I
               TEMP  = X (I)
               X (I) = X (J)
               X (J) = TEMP
               TEMP  = Y (I)
               Y (I) = Y (J)
               Y (J) = TEMP
               TEMP  = Z (I)
               Z (I) = Z (J)
               Z (J) = TEMP
            END DO

         ELSE IF (CHOICE == 11) THEN  ! "Switch X and Y":

            DO I = 1, N
               TEMP = X (I)
               X (I) = Y (I)
               Y (I) = TEMP
            END DO

         ELSE IF (CHOICE == 12) THEN  ! "Switch Y and Z":

            DO I = 1, N
               TEMP = Y (I)
               Y (I) = Z (I)
               Z (I) = TEMP
            END DO

         ELSE IF (CHOICE == 13) THEN  ! "Switch Z and X":

            DO I = 1, N
               TEMP = X (I)
               X (I) = Z (I)
               Z (I) = TEMP
            END DO
        
         ELSE IF (CHOICE == 14) THEN  ! "Scale X & Y & Z"

            IF (PROMPT) THEN
               CALL READR (LUNCRT, '   Enter scale (+ or -): ',
     >            LUNKBD, SCALE, CR, EOF)
               IF (CR .OR. EOF) GO TO 999
            END IF

            DO I = 1, N
               X (I) = X (I) * SCALE
               Y (I) = Y (I) * SCALE
               Z (I) = Z (I) * SCALE
            END DO

         ELSE IF (CHOICE == 15) THEN  ! "Rotate (X,Y)":

            IF (PROMPT) THEN
               CALL READR (LUNCRT, '   Degrees anticlockwise: ',
     >            LUNKBD, ANGLE, CR, EOF)
               IF (CR .OR. EOF) GO TO 999

  150          WRITE (LUNCRT, '(A)', ADVANCE='NO')'  (Center (Xc, Yc): '
               READ  (LUNKBD, *, ERR=150) P, Q
               WRITE (LUNCRT, '(A)')
            END IF

            CALL ROTATE2D (N, X, Y, ANGLE, P, Q)

         ELSE IF (CHOICE == 16) THEN  ! "Rotate (Y,Z)":

            IF (PROMPT) THEN
               CALL READR (LUNCRT, '   Degrees anticlockwise: ',
     >            LUNKBD, ANGLE, CR, EOF)
               IF (CR .OR. EOF) GO TO 999

  160          WRITE (LUNCRT, '(A)', ADVANCE='NO')'  (Center (Yc, Zc): '
               READ  (LUNKBD, *, ERR=160) P, Q
               WRITE (LUNCRT, '(A)')
            END IF

            CALL ROTATE2D (N, Y, Z, ANGLE, P, Q)

         ELSE IF (CHOICE == 17) THEN  ! "Rotate (Z,X)":

            IF (PROMPT) THEN
               CALL READR (LUNCRT, '   Degrees anticlockwise: ',
     >            LUNKBD, ANGLE, CR, EOF)
               IF (CR .OR. EOF) GO TO 999

  170          WRITE (LUNCRT, '(A)', ADVANCE='NO')'  (Center (Zc, Xc): '
               READ  (LUNKBD, *, ERR=170) P, Q
               WRITE (LUNCRT, '(A)')
            END IF

            CALL ROTATE2D (N, Z, X, ANGLE, P, Q)

         ELSE IF (CHOICE == 18) THEN  ! "Convert inches to meters"

            DO I = 1, N
               X (I) = X (I) * INCHES_TO_METERS
               Y (I) = Y (I) * INCHES_TO_METERS
               Z (I) = Z (I) * INCHES_TO_METERS
            END DO

         ELSE IF (CHOICE == 19) THEN  ! "Convert meters to inches"

            DO I = 1, N
               X (I) = X (I) * METERS_TO_INCHES
               Y (I) = Y (I) * METERS_TO_INCHES
               Z (I) = Z (I) * METERS_TO_INCHES
            END DO

         END IF

  999    RETURN

         END SUBROUTINE TRANSFORM

      END PROGRAM SCALE_GRID
