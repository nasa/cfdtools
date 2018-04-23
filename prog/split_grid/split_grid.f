C-------------------------------------------------------------------------------
C
      PROGRAM SPLIT_GRID
C
C  Description:
C
C     SPLIT_GRID splits a 3-D multiblock grid into another multiblock form, with
C     options to reverse any of the index directions, change the i, j, k order,
C     and/or suppress input blocks.  An accompanying PLOT3D-type function file
C     may be processed in the same run.
C
C     A preliminary option is provided to scan a grid file to be split so as
C     to generate a template control file for the splitting.
C
C  Sample input control file (split.inp):
C
C     SPLIT_GRID controls  ! Input blocks in order; omitted blocks suppressed
C     mygrid.g             ! Input grid file name
C     T                    ! T = formatted; F = unformatted
C     myfunction.f         ! Input function file name, or 'none'
C     T                    ! T = formatted; F = unformatted
C     splitgrid.g          ! Output grid file name
C     T                    ! T = formatted; F = unformatted
C     splitfunction.f      ! Output function file name (if present)
C     T                    ! T = formatted; F = unformatted
C     block   isplits jsplits ksplits   output ijk order; splits apply to input
C     1       6       2       2         1 2 3
C     i:      1,49 49,81 81,113 145,113 177,145 225,177
C     j:      1,49 49,97
C     k:      1,17 17,33
C     block   isplits jsplits ksplits   ijk order
C     2       6       1       1         2 1 3
C     i:      1,49 49,81 81,113 145,113 177,145 225,177
C     j:      1,97
C     k:      33,1
C     block   isplits jsplits ksplits   ijk order
C     4       1       1       1         1 2 3
C     i:      1,225
C     j:      1,97
C     k:      1,33
C
C  Procedures:
C
C     XYZQ_IO package  I/O utilities for CFD files (David Saunders)
C
C  History:
C
C     ??/??/??   J.Reuther/  Original SPLIT_XYZ.
C                M.Rimlinger
C
C     10/25/99-  D.Saunders  Incorporated CFD_IO_PACKAGE; provided for
C     11/10/99               reversing the order in i, j, k directions;
C                            provided for permuting the i, j, k order;
C                            more compact control file + template option.
C     02/24/00     "   "     A short control file needed the final loop over
C                            included blocks, not over all blocks.
C     03/05/04     "   "     Set up of IJKO(1:3,*) wasn't right!
C     04/12/04     "   "     Replaced defunct SYN107-MB restart file option
C                            with PLOT3D function file option.
C     04/27/04     "   "     The function file option didn't work.
C                            Replaced CFD_IO_PACKAGE with XYZQ_IO package.
C
C  Author:  David Saunders, ELORET/NASA Ames Research Center, Moffett Field, CA
C
C-------------------------------------------------------------------------------

C     Modules:

      USE GRID_BLOCK_STRUCTURE
      USE XYZQ_IO_MODULE

      IMPLICIT NONE

C     Constants:

      INTEGER, PARAMETER ::
     >   LUNCTL = 1, LUNXI = 2, LUNXO = 3, LUNKBD = 5, LUNCRT = 6,
     >   LUNFI  = 7, LUNFO = 8

      LOGICAL, PARAMETER ::
     >   TRUE = .TRUE.

      CHARACTER, PARAMETER ::
     >   FORMAT * 11 = 'unformatted',
     >   NONE   * 4  = 'none'

C     Variables:

      INTEGER
     >   I3, IB, IBLOCK, IOS, ISPLIT,
     >   I, I1, I2, IINC, IS, ISPLITS,
     >   J, J1, J2, JINC, JS, JSPLITS,
     >   K, K1, K2, KINC, KS, KSPLITS,
     >   NB, NBSPLIT, NBLOCKI, NBLOCKO, NF, NPTS, NSPLITS

      LOGICAL
     >   FFILE, FORMATTED_FI, FORMATTED_FO, FORMATTED_XI, FORMATTED_XO,
     >   XFILE

      CHARACTER
     >   ANSWER * 1, FILENAME * 80, LABEL * 2

      CHARACTER ::
     >   CIJK (3) * 1 = (/'i', 'j', 'k'/)

      INTEGER, POINTER ::
     >   L, M, N

      INTEGER, POINTER, DIMENSION (:) ::
     >   ISPLIT1, ISPLIT2, JSPLIT1, JSPLIT2, KSPLIT1, KSPLIT2

      INTEGER, TARGET ::
     >   II, JJ, KK

      INTEGER, TARGET, DIMENSION (3) ::
     >   IJKORDER

      INTEGER, POINTER, DIMENSION (:,:) ::
     >   IJKI, IJKO

      TYPE (GRID_TYPE), POINTER, DIMENSION (:) ::
     >   GRIDI, GRIDO

C     Execution:
C     ----------

      WRITE (LUNCRT, '(/, A)', ADVANCE='NO')
     >   ' Generate a template control file? (Y/N): '
      READ (LUNKBD, '(A)') ANSWER

      IF (ANSWER == 'Y' .or. ANSWER == 'y') THEN

         CALL TEMPLATE () ! Internal procedure

         GO TO 999 ! Done

      END IF


      OPEN (LUNCTL, FILE='split.inp', STATUS='OLD', IOSTAT=IOS)

      IF (IOS /= 0) THEN
         WRITE (LUNCRT, '(/, A)')
     >      ' Unable to open control file split.inp.'
         GO TO 999
      END IF

      READ (LUNCTL, *) ! Skip header
      READ (LUNCTL, *) FILENAME
      READ (LUNCTL, *) FORMATTED_XI

      XFILE = FILENAME(1:4) /= none

      IF (.NOT. XFILE) THEN
         WRITE (LUNCRT, '(/, 2A)') ' Processing just a function file',
     >      ' is disallowed by the XYZQ_IO package.'
         GO TO 999
      END IF

C     Read the active headers in preparation for processing one block at a
C     time.  There is no need to store all input or output blocks.

      I1 = 1;  IF (FORMATTED_XI) I1 = 3

      OPEN (LUNXI, FILE=FILENAME, STATUS='OLD', FORM=FORMAT(I1:11),
     >      IOSTAT=IOS)

      IF (IOS /= 0) THEN
         WRITE (LUNCRT, '(/, 2A)') ' Unable to open the input grid: ',
     >      FILENAME(1:LEN_TRIM(FILENAME))
         GO TO 999
      END IF

      CALL XYZ_HEADER_IO (1, LUNXI, FORMATTED_XI, NBLOCKI, GRIDI, IOS)

      IF (IOS /= 0) THEN
         WRITE (LUNCRT, '(/, A)')
     >      ' Trouble allocating or reading the input grid header.'
         GO TO 999
      END IF

      READ (LUNCTL, *) FILENAME
      READ (LUNCTL, *) FORMATTED_FI

      FFILE = FILENAME(1:4) /= none

      IF (FFILE) THEN

         I1 = 1;  IF (FORMATTED_FI) I1 = 3

         OPEN (LUNFI, FILE=FILENAME, STATUS='OLD', FORM=FORMAT(I1:11),
     >         IOSTAT=IOS)

         IF (IOS /= 0) THEN
            WRITE (LUNCRT, '(/, 2A)')
     >         ' Unable to open the input function file: ',
     >         FILENAME(1:LEN_TRIM(FILENAME))
            GO TO 999
         END IF

         CALL Q_HEADER_IO (1, LUNFI, FORMATTED_FI, NBLOCKI, NF, GRIDI,
     >                     IOS)
         IF (IOS /= 0) THEN
            WRITE (LUNCRT, '(/, A)')
     >         ' Trouble reading input function file header.'
            GO TO 999
         END IF

         IF (GRIDI(1)%MI /= GRIDI(1)%NI) THEN
            WRITE (LUNCRT, '(/, A)')
     >         ' The function data must be at the grid vertices.'
            GO TO 999
         END IF

      END IF

      READ (LUNCTL, *) FILENAME
      READ (LUNCTL, *) FORMATTED_XO

      I1 = 1;  IF (FORMATTED_XO) I1 = 3

      OPEN (LUNXO, FILE=FILENAME, STATUS='NEW', FORM=FORMAT(I1:11),
     >      IOSTAT=IOS)

      IF (IOS /= 0) THEN
         WRITE (LUNCRT, '(/, 2A)') ' Unable to open output grid: ',
     >      FILENAME(1:LEN_TRIM(FILENAME))
         GO TO 999
      END IF

      READ (LUNCTL, *) FILENAME
      READ (LUNCTL, *) FORMATTED_FO

      IF (FFILE) THEN

         I1 = 1;  IF (FORMATTED_FO) I1 = 3

         OPEN (LUNFO, FILE=FILENAME, STATUS='NEW',
     >         FORM=FORMAT(I1:11), IOSTAT=IOS)

         IF (IOS /= 0) THEN
            WRITE (LUNCRT, '(/, 2A)')
     >         ' Unable to open output function file: ',
     >         FILENAME(1:LEN_TRIM(FILENAME))
            GO TO 999
         END IF

      END IF

C     It's easiest to retain some of the CFD_IO package's variables:

      ALLOCATE (IJKI(3,NBLOCKI))

      DO IB = 1, NBLOCKI
         IJKI(1,IB) = GRIDI(IB)%NI
         IJKI(2,IB) = GRIDI(IB)%NJ
         IJKI(3,IB) = GRIDI(IB)%NK
      END DO

C     A first scan of the control file determines the output number of blocks:

      NBSPLIT = 0 ! # input blocks to be split (<= NBLOCKI)
      NBLOCKO = 0 ! # output blocks
      NSPLITS = 1 ! Most splits of any block in any direction

      DO ! Until EOF

         CALL SKIP (1)     ! Internal procedure

         IF (IOS < 0) EXIT ! Normal EOF

         IF (IOS /= 0) THEN
            WRITE (LUNCRT, '(/, A, I6)')
     >         ' Error skipping split definition header; block #:', IB
            GO TO 999
         END IF

         READ (LUNCTL, *, IOSTAT=IOS) IB, ISPLITS, JSPLITS, KSPLITS

         IF (IOS /= 0) THEN
            WRITE (LUNCRT, '(/, A, I6)')
     >         ' Error reading split definitions; block #:', IB
            GO TO 999
         END IF

         NBSPLIT = NBSPLIT + 1

         IF (NBSPLIT > NBLOCKI) THEN
            WRITE (LUNCRT, '(/, A, /, (A, I6))')
     >         ' *** Control file has too many sets of splits.',
     >         ' # blocks in X/F file:  ', NBLOCKI,
     >         ' # sets of block splits:', NBSPLIT
            GO TO 999
         END IF

         NBLOCKO = NBLOCKO   +   ISPLITS * JSPLITS * KSPLITS
         NSPLITS = MAX (NSPLITS, ISPLITS,  JSPLITS,  KSPLITS)

         CALL SKIP (3)

      END DO

      REWIND (LUNCTL)
      CALL SKIP (9)

      ALLOCATE (IJKO(3,NBLOCKO),
     >          ISPLIT1(NSPLITS), ISPLIT2(NSPLITS), JSPLIT1(NSPLITS),
     >          JSPLIT2(NSPLITS), KSPLIT1(NSPLITS), KSPLIT2(NSPLITS))


C     A second scan of the control file allows set-up of the output header info:

      NB = 0 ! Output block number

      DO ISPLIT = 1, NBSPLIT

         CALL SKIP (1)

         READ (LUNCTL, *, IOSTAT=IOS) IB, ISPLITS, JSPLITS, KSPLITS,
     >                                IJKORDER
         IF (IOS /= 0) THEN
            WRITE (LUNCRT, '(/, A, I6)')
     >         ' Error reading ijk order; input block #:', IB
            GO TO 999
         END IF

         READ (LUNCTL, *, IOSTAT=IOS) ! Splits refer to input i,j,k order
     >      LABEL, (ISPLIT1(IS), ISPLIT2(IS), IS = 1, ISPLITS),
     >      LABEL, (JSPLIT1(JS), JSPLIT2(JS), JS = 1, JSPLITS),
     >      LABEL, (KSPLIT1(KS), KSPLIT2(KS), KS = 1, KSPLITS)

         IF (IOS /= 0) THEN
            WRITE (LUNCRT, '(/, A, I5)')
     >         ' Error reading split indices; input block #:', IB
            GO TO 999
         END IF

C        Protect the user from obviously wrong split indices:

         DO IS = 1, ISPLITS
            IF (MAX (ISPLIT1(IS), ISPLIT2(IS)) > IJKI(1,IB)) THEN
               I = 1
               GO TO 800
            END IF
         END DO

         DO JS = 1, JSPLITS
            IF (MAX (JSPLIT1(JS), JSPLIT2(JS)) > IJKI(2,IB)) THEN
               I = 2
               GO TO 800
            END IF
         END DO

         DO KS = 1, KSPLITS
            IF (MAX (KSPLIT1(KS), KSPLIT2(KS)) > IJKI(3,IB)) THEN
               I = 3
               GO TO 800
            END IF
         END DO

C        Handle changing the order of i, j, k during splitting:

CCCC     L => IJKORDER(1)   ! No!  Assigning IJKO(L/M/N,NB) doesn't work
CCCC     M => IJKORDER(2)
CCCC     N => IJKORDER(3)

         DO I = 1, 3
            J = IJKORDER(I)
            IF (J == 1) THEN
               I1 = I
            ELSE IF (J == 2) THEN
               I2 = I
            ELSE IF (J == 3) THEN
               I3 = I
            ELSE
               GO TO 900
            END IF
         END DO

         DO KS = 1, KSPLITS
            DO JS = 1, JSPLITS
               DO IS = 1, ISPLITS
                  NB = NB + 1
                  IJKO(I1,NB) = ABS (ISPLIT2(IS) - ISPLIT1(IS)) + 1
                  IJKO(I2,NB) = ABS (JSPLIT2(JS) - JSPLIT1(JS)) + 1
                  IJKO(I3,NB) = ABS (KSPLIT2(KS) - KSPLIT1(KS)) + 1
               END DO
            END DO
         END DO

      END DO ! Next block to be split

      REWIND (LUNCTL)
      CALL SKIP (9)


C     Write the header records for the active output file(s):

      ALLOCATE (GRIDO(NBLOCKO))

      DO IB = 1, NBLOCKO
         GRIDO(IB)%NI = IJKO(1,IB)
         GRIDO(IB)%NJ = IJKO(2,IB)
         GRIDO(IB)%NK = IJKO(3,IB)
      END DO

      CALL XYZ_HEADER_IO (2, LUNXO, FORMATTED_XO, NBLOCKO, GRIDO,
     >                    IOS)
      IF (IOS /= 0) THEN
         WRITE (LUNCRT, '(/, A)')
     >      ' Trouble writing the output grid header.'
         GO TO 999
      END IF

      IF (FFILE) THEN

         DO IB = 1, NBLOCKO
            GRIDO(IB)%MI = IJKO(1,IB)
            GRIDO(IB)%MJ = IJKO(2,IB)
            GRIDO(IB)%MK = IJKO(3,IB)
         END DO

         CALL Q_HEADER_IO (2, LUNFO, FORMATTED_FO, NBLOCKO, NF, GRIDO,
     >                     IOS)
         IF (IOS /= 0) THEN
            WRITE (LUNCRT, '(/, A)')
     >         ' Trouble writing the output function file header.'
            GO TO 999
         END IF

      END IF


C     Process one input block at a time for all active file types:

      IBLOCK = 1 ! Input  block number
      NB     = 0 ! Output block number

      DO ISPLIT = 1, NBSPLIT

         CALL SKIP (1)

         READ (LUNCTL, *) IB, ISPLITS, JSPLITS, KSPLITS, IJKORDER

         READ (LUNCTL, *)
     >      LABEL, (ISPLIT1(IS), ISPLIT2(IS), IS = 1, ISPLITS),
     >      LABEL, (JSPLIT1(JS), JSPLIT2(JS), JS = 1, JSPLITS),
     >      LABEL, (KSPLIT1(KS), KSPLIT2(KS), KS = 1, KSPLITS)

C        Find the indicated input block:

         DO WHILE (IBLOCK <= IB) ! Read a block (no way of skipping)

            CALL XYZ_ALLOCATE (GRIDI(IBLOCK), IOS)
            IF (IOS /= 0) GO TO 999

            NPTS = IJKI(1,IBLOCK) * IJKI(2,IBLOCK) * IJKI(3,IBLOCK)

            CALL XYZ_BLOCK_IO (1, LUNXI, FORMATTED_XI, NPTS,
     >                         GRIDI(IBLOCK)%X, GRIDI(IBLOCK)%Y,
     >                         GRIDI(IBLOCK)%Z, IOS)
            IF (IOS /= 0) GO TO 999

            IF (IBLOCK < IB) THEN
               DEALLOCATE (GRIDI(IBLOCK)%X, GRIDI(IBLOCK)%Y,
     >                     GRIDI(IBLOCK)%Z)
            END IF

            IF (FFILE) THEN

               CALL Q_ALLOCATE (GRIDI(IBLOCK), NF, IOS)
               IF (IOS /= 0) GO TO 999

               CALL Q_BLOCK_IO (1, LUNFI, FORMATTED_FI, NF,
     >                          GRIDI(IBLOCK)%MI, GRIDI(IBLOCK)%MJ,
     >                          GRIDI(IBLOCK)%MK, GRIDI(IBLOCK)%Q, IOS)
               IF (IOS /= 0) GO TO 999

               IF (IBLOCK < IB) DEALLOCATE (GRIDI(IBLOCK)%Q)

            END IF

            IBLOCK = IBLOCK + 1

         END DO

C        Arrange for permuting i,j,k in the output block:

         I = IJKORDER(1)
         IF (I == 1) THEN
            L => II
         ELSE IF (I == 2) THEN
            L => JJ
         ELSE ! I == 3
            L => KK
         END IF

         J = IJKORDER(2)
         IF (J == 1) THEN
            M => II
         ELSE IF (J == 2) THEN
            M => JJ
         ELSE ! J == 3
            M => KK
         END IF

         K = IJKORDER(3)
         IF (K == 1) THEN
            N => II
         ELSE IF (K == 2) THEN
            N => JJ
         ELSE ! K == 3
            N => KK
         END IF

         DO KS = 1, KSPLITS

            K1 = KSPLIT1(KS)
            K2 = KSPLIT2(KS)
            KINC = SIGN (1, K2 - K1)

            DO JS = 1, JSPLITS

               J1 = JSPLIT1(JS)
               J2 = JSPLIT2(JS)
               JINC = SIGN (1, J2 - J1)

               DO IS = 1, ISPLITS

                  I1 = ISPLIT1(IS)
                  I2 = ISPLIT2(IS)
                  IINC = SIGN (1, I2 - I1)

                  NB = NB + 1

                  CALL XYZ_ALLOCATE (GRIDO(NB), IOS)
                  IF (IOS /= 0) GO TO 999

                  KK = 0
                  DO K = K1, K2, KINC
                     KK = KK + 1
                     JJ = 0
                     DO J = J1, J2, JINC
                        JJ = JJ + 1
                        II = 0
                        DO I = I1, I2, IINC
                           II = II + 1
                           GRIDO(NB)%X(L,M,N) = GRIDI(IB)%X(I,J,K)
                           GRIDO(NB)%Y(L,M,N) = GRIDI(IB)%Y(I,J,K)
                           GRIDO(NB)%Z(L,M,N) = GRIDI(IB)%Z(I,J,K)
                        END DO
                     END DO
                  END DO

                  NPTS = IJKO(1,NB) * IJKO(2,NB) * IJKO(3,NB)

                  CALL XYZ_BLOCK_IO (2, LUNXO, FORMATTED_XO, NPTS,
     >                               GRIDO(NB)%X, GRIDO(NB)%Y,
     >                               GRIDO(NB)%Z, IOS)
                  IF (IOS /= 0) GO TO 999

                  DEALLOCATE (GRIDO(NB)%X, GRIDO(NB)%Y, GRIDO(NB)%Z)

                  IF (FFILE) THEN

                     CALL Q_ALLOCATE (GRIDO(NB), NF, IOS)

                     IF (IOS /= 0) GO TO 999

                     KK = 0
                     DO K = K1, K2, KINC
                        KK = KK + 1
                        JJ = 0
                        DO J = J1, J2, JINC
                           JJ = JJ + 1
                           II = 0
                           DO I = I1, I2, IINC
                             II = II + 1
                             GRIDO(NB)%Q(:,L,M,N) = GRIDI(IB)%Q(:,I,J,K)
                           END DO
                        END DO
                     END DO

                     CALL Q_BLOCK_IO (2, LUNFO, FORMATTED_FO, NF,
     >                                GRIDO(NB)%MI, GRIDO(NB)%MJ,
     >                                GRIDO(NB)%MK, GRIDO(NB)%Q, IOS)
                     IF (IOS /= 0) GO TO 999

                     DEALLOCATE (GRIDO(NB)%Q)

                  END IF

               END DO ! Next I split

            END DO ! Next J split

         END DO ! Next K split

         DEALLOCATE (GRIDI(IB)%X, GRIDI(IB)%Y, GRIDI(IB)%Z)

         IF (FFILE) DEALLOCATE (GRIDI(IB)%Q)

      END DO ! Next input block

      GO TO 999


C     Error handling:

  800 WRITE (LUNCRT, '(/, 1X, A1, A, I5)')
     >   CIJK(I), ' split indices are bad; input block #:', IB
      GO TO 999

  900 WRITE (LUNCRT, '(/, (A, 3I5))')
     >   ' Bad IJK order: ', IJKORDER, ' Input block #:', IB

  999 CONTINUE ! Avoid system dependencies with "STOP"


C     Internal procedures for SPLIT_GRID:

      CONTAINS

!        ---------------------------------------------------------
         SUBROUTINE SKIP (NLINES) ! Skip control file header lines
!        ---------------------------------------------------------

         INTEGER L, NLINES

         DO L = 1, NLINES
            READ (LUNCTL, *, IOSTAT=IOS)
         END DO

         END SUBROUTINE SKIP

!        ----------------------------------------------------------
         SUBROUTINE TEMPLATE () ! Generate a template for split.inp
!        ----------------------------------------------------------

         CHARACTER, PARAMETER ::
     >      ORDER * 14 = '         1 2 3'

         CHARACTER ::
     >      CFORMAT * 16 = '(I1, I8, 2I8, A)'

!        Execution:

         OPEN (LUNCTL, FILE='split.inp', STATUS='NEW', IOSTAT=IOS)

         IF (IOS /= 0) THEN
            WRITE (LUNCRT, '(/, A)')
     >         ' Unable to open new control file split.inp.'
            GO TO 999
         END IF

!        Prompt for a grid file and open it:

         CALL FILE_PROMPT (LUNXI, 0, 'grid to be split', 'OLD', TRUE,
     >                     FORMATTED_XI, IOS)

         IF (IOS /= 0) GO TO 999

!        Allocate an array of blocks, and read the header info:

         CALL XYZ_HEADER_IO (1, LUNXI, FORMATTED_XI, NBLOCKI, GRIDI,
     >                       IOS)
         IF (IOS /= 0) GO TO 999

         WRITE (LUNCTL, '(2A)')
     >      'SPLIT_GRID controls      ',
     >      '! Input blocks in order; omitted blocks suppressed'
         WRITE (LUNCTL, '(A)')
     >      'mygrid.g                 ! Input grid file name',
     >      'T                        ! T = formatted; F = unformatted',
     >      'none                     ! Input function file, or "none"',
     >      'T                        ! T = formatted; F = unformatted',
     >      'splitgrid.g              ! Output grid file name',
     >      'T                        ! T = formatted; F = unformatted',
     >      'splitfunction.f          ! Output func. file (if present)',
     >      'T                        ! T = formatted; F = unformatted'

         DO IB = 1, NBLOCKI

            WRITE (LUNCTL, '(A)')
     >         'block   isplits jsplits ksplits   ijk order'

            IF (IB > 9 .AND. IB <= 99) THEN
               CFORMAT(2:7) = 'I2, I7'
            ELSE IF (IB > 99) THEN
               CFORMAT(2:7) = 'I3, I6'
            END IF

            WRITE (LUNCTL, CFORMAT) IB, 1, 1, 1, ORDER

            WRITE (LUNCTL, '(A1, '':      1,'', I3)')
     >         CIJK(1), GRIDI(IB)%NI,
     >         CIJK(2), GRIDI(IB)%NJ,
     >         CIJK(3), GRIDI(IB)%NK

         END DO

         WRITE (LUNCRT, '(A)') ' Template generated: split.inp.'

  999    CONTINUE

         END SUBROUTINE TEMPLATE

      END PROGRAM SPLIT_GRID
