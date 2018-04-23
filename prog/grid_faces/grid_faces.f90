!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      PROGRAM GRID_FACES
!
!     Description:
!
!        GRID_FACES extracts the indicated block faces from the input grid 
!     and/or function file.  All files are 3-space PLOT3D multiblock type.
!
!     Control file:
!
!        'grid_faces.inp' should contain the input and output file specs.
!     followed by a list of (block, face) pairs, where faces 1 - 6 mean 
!     imin, imax, jmin, jmax, kmin, kmax respectively.  E.g.:
!
!     GRID_FACES control file
!     ---------- Input grid file ----------
!     volume.g   ! or none
!     T          ! Formatted? [T|F]
!     ---------- Input function file ------
!     volume.f   ! or none
!     T          ! Formatted? [T|F]
!     ---------- Output grid file ---------
!     surface.g  ! or none
!     T          ! Formatted? [T|F]
!     ---------- Output function file -----
!     surface.f  ! or none
!     T          ! Formatted? [T|F]
!     ---------- Block/face list ----------
!        171 3   ! Read one per line to EOF
!        172 3  
!        173 3
!         :  :
!        229 5
!        230 5
!
!     Procedures:
!
!        XYZQ_IO          PLOT3D-type I/O package
!        COPY_FACE        Wouldn't work as an internal procedure
!
!     History:
!
!        06/19/03  D. Saunders  Initial adaptation of SCALE_GRID.
!        06/01/04       "       Converting face # to dimension # was wrong!
!        12/14/04       "       Switched from CFD_IO_PACKAGE to XYZQ_IO package
!                               in order to treat function files as for grids.
!
!     Author:  David Saunders, Eloret/NASA Ames, Moffett Field, CA.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Modules:

      USE GRID_BLOCK_STRUCTURE
      USE XYZQ_IO_MODULE

      IMPLICIT NONE

!     Constants:

      INTEGER, PARAMETER :: &
         LUNCTL  = 1,       &
         LUNFIN  = 2,       &
         LUNGIN  = 3,       &
         LUNKBD  = 5,       &
         LUNCRT  = 6,       &
         LUNFOUT = 7,       &
         LUNGOUT = 8

      CHARACTER, PARAMETER ::         &
         FORMAT * 11 = 'UNFORMATTED', &
         NONE   *  4 = 'none'

!     Variables:

      INTEGER :: &
         I, I1, I2, IB, ID, IF, IOS, L1, N, NBLOCK_IN, NBLOCK_OUT, NBLOCK_RD, &
         NPTS, NVARS

      LOGICAL :: &
         FORMATTED_FIN, FORMATTED_FOUT, FORMATTED_GIN, FORMATTED_GOUT,        &
         FUNC, GRID

      CHARACTER :: &
         FILENAME * 80

      INTEGER, ALLOCATABLE, DIMENSION (:) :: &
         IBLOCK, IFACE

      TYPE (GRID_TYPE), POINTER, DIMENSION (:) :: &
         X_IN, X_OUT

!     Execution:

      OPEN (LUNCTL, FILE='grid_faces.inp', STATUS='OLD', IOSTAT=IOS)

      IF (IOS /= 0) THEN
         WRITE (LUNCRT, '(/, A)') &
            ' Unable to open grid_faces.inp control file.'
         GO TO 999
      END IF

      READ (LUNCTL, *)
      READ (LUNCTL, *) ! Input volume grid file
      READ (LUNCTL, *) FILENAME
      GRID = FILENAME /= NONE
      READ (LUNCTL, *) FORMATTED_GIN
      L1 = 1;  IF (FORMATTED_GIN) L1 = 3

      IF (GRID) THEN
         OPEN (LUNGIN, FILE=FILENAME, STATUS='OLD', FORM=FORMAT(L1:11),       &
               IOSTAT=IOS)
         IF (IOS /= 0) THEN
            WRITE (LUNCRT, '(/, 2A)') ' Unable to open input grid: ',         &
               TRIM (FILENAME)
            GO TO 999
         END IF
      END IF

      READ (LUNCTL, *) ! Input volume function file
      READ (LUNCTL, *) FILENAME
      FUNC = FILENAME /= NONE
      READ (LUNCTL, *) FORMATTED_FIN
      L1 = 1;  IF (FORMATTED_FIN) L1 = 3

      IF (FUNC) THEN
         OPEN (LUNFIN, FILE=FILENAME, STATUS='OLD', FORM=FORMAT(L1:11),       &
               IOSTAT=IOS)
         IF (IOS /= 0) THEN
            WRITE (LUNCRT, '(/, 2A)') ' Unable to open input function file: ',&
               TRIM (FILENAME)
            GO TO 999
         END IF
      END IF

      READ (LUNCTL, *) ! Output surface grid file
      READ (LUNCTL, *) FILENAME
      READ (LUNCTL, *) FORMATTED_GOUT
      L1 = 1;  IF (FORMATTED_GOUT) L1 = 3

      IF (GRID) THEN
         OPEN (LUNGOUT, FILE=FILENAME, &    ! Avoid an a2ps line wrap
               STATUS='UNKNOWN', FORM=FORMAT(L1:11), IOSTAT=IOS)
         IF (IOS /= 0) THEN
            WRITE (LUNCRT, '(/, 2A)') ' Unable to open output grid: ',        &
               TRIM (FILENAME)
            GO TO 999
         END IF
      END IF

      READ (LUNCTL, *) ! Output surface function file
      READ (LUNCTL, *) FILENAME
      READ (LUNCTL, *) FORMATTED_FOUT
      L1 = 1;  IF (FORMATTED_FOUT) L1 = 3

      IF (FUNC) THEN
         OPEN (LUNFOUT, FILE=FILENAME, &
               STATUS='UNKNOWN', FORM=FORMAT(L1:11), IOSTAT=IOS)
         IF (IOS /= 0) THEN
            WRITE (LUNCRT, '(/, 2A)') ' Cannot open output function file: ',  &
               TRIM (FILENAME)
            GO TO 999
         END IF
      END IF

      READ (LUNCTL, *) ! Header for list of block/face pairs

      NBLOCK_OUT = 0   ! We have to count them before we can store them

      DO ! Until EOF
         READ (LUNCTL, *, IOSTAT=IOS) I1, I2
         IF (IOS == 0) THEN
            NBLOCK_OUT = NBLOCK_OUT + 1
            IF (I2 < 1 .OR. I2 > 6) THEN
               WRITE (LUNCRT, '(A, I6, A, I6)') &
                  ' Bad face number = ', I2, ',  line #: ', NBLOCK_OUT
                  GO TO 999
            END IF
         ELSE IF (IOS < 0) THEN
            EXIT
         ELSE
            WRITE (LUNCRT, '(A, I5)') &
               ' Bad integer in control file.  Line #: ', NBLOCK_OUT + 1
            GO TO 999
         END IF
      END DO

      WRITE (LUNCRT, '(/, A, I8)') ' # faces being extracted:', NBLOCK_OUT

      ALLOCATE (IBLOCK(NBLOCK_OUT), IFACE(NBLOCK_OUT))

      REWIND (LUNCTL)

      DO I = 1, 14
         READ (LUNCTL, *) ! Skip to start of list
      END DO

      NBLOCK_RD = 0
      DO IB = 1, NBLOCK_OUT
         READ (LUNCTL, *) IBLOCK(IB), IFACE(IB)
         NBLOCK_RD = MAX (IBLOCK(IB), NBLOCK_RD) ! Avoid reading unneeded blocks
      END DO

      CLOSE (LUNCTL)

      WRITE (LUNCRT, '(A, I5)') ' # blocks that will be read:', NBLOCK_RD

!     Read a grid if it is present.  The fact that it might NOT be present
!     (function file input only) conflicts with an assumption in XYZQ_IO about
!     always reading a grid before a function file, so we have to work a bit
!     harder than if we assume we always have a grid.

      IF (GRID) THEN

!        Read the number of input blocks, allocate them, and read dimensions:

         CALL XYZ_HEADER_IO (1, LUNGIN, FORMATTED_GIN, NBLOCK_IN, X_IN, IOS)

         WRITE (LUNCRT, '(A, I5)') ' # blocks found:', NBLOCK_IN

         IF (IOS /= 0) THEN
            WRITE (LUNCRT, '(/, A)') ' Trouble reading grid header.'
            GO TO 999
         END IF

         DO IB = 1, NBLOCK_RD ! <= NBLOCK_IN

            CALL XYZ_ALLOCATE (X_IN(IB), IOS)

            IF (IOS /= 0) THEN
               WRITE (LUNCRT, '(/, A, I5)') &
                  ' Trouble allocating grid (x,y,z)s.  Block #:', IB
               GO TO 999
            END IF

            NPTS = X_IN(IB)%NI * X_IN(IB)%NJ * X_IN(IB)%NK

            CALL XYZ_BLOCK_IO (1, LUNGIN, FORMATTED_GIN, NPTS,                &
                               X_IN(IB)%X, X_IN(IB)%Y, X_IN(IB)%Z, IOS)
            IF (IOS /= 0) THEN
               WRITE (LUNCRT, '(/, A, I5)') &
                  ' Trouble reading grid (x,y,z)s.  Block #:', IB
               GO TO 999
            END IF

         END DO ! Next input grid block

         CLOSE (LUNGIN)

      END IF

!     Read a function file?

      IF (FUNC) THEN

         IF (GRID) THEN ! Blocks have already been allocated

            CALL Q_HEADER_IO (1, LUNFIN, FORMATTED_FIN, NBLOCK_IN, NVARS,     &
                              X_IN, IOS)
            IF (IOS /= 0) THEN
               WRITE (LUNCRT, '(/, A)') ' Trouble reading function file header.'
               GO TO 999
            END IF

         ELSE ! Read the number of blocks, allocate them, and read dimensions:

            IF (FORMATTED_FIN) THEN
               READ (LUNFIN, *, IOSTAT=IOS) NBLOCK_IN
            ELSE
               READ (LUNFIN,    IOSTAT=IOS) NBLOCK_IN
            END IF

            IF (IOS /= 0) THEN
               WRITE (LUNCRT, '(/, A)') ' Trouble reading # func. file blocks.'
               GO TO 999
            END IF

            REWIND (LUNFIN)

            ALLOCATE (X_IN(NBLOCK_IN))

            CALL Q_HEADER_IO (1, LUNFIN, FORMATTED_FIN, NBLOCK_IN, NVARS,     &
                              X_IN, IOS)
            IF (IOS /= 0) THEN
               WRITE (LUNCRT, '(/, A)') ' Trouble reading function file header.'
               GO TO 999
            END IF

         END IF

         DO IB = 1, NBLOCK_RD ! <= NBLOCK_IN

            CALL Q_ALLOCATE (X_IN(IB), NVARS, IOS)

            IF (IOS /= 0) THEN
               WRITE (LUNCRT, '(/, A, I5)') &
                  ' Trouble allocating function file variables.  Block #:', IB
               GO TO 999
            END IF

            CALL Q_BLOCK_IO (1, LUNFIN, FORMATTED_FIN, NVARS,                 &
                             X_IN(IB)%MI, X_IN(IB)%MJ, X_IN(IB)%MK,           &
                             X_IN(IB)%Q, IOS)
            IF (IOS /= 0) THEN
               WRITE (LUNCRT, '(/, A, I5)') &
                  ' Trouble reading function variables.  Block #:', IB
               GO TO 999
            END IF

         END DO ! Next input function block

         CLOSE (LUNFIN)

      END IF

!     Display the block dimensions to reassure the user:

      WRITE (LUNCRT, '(A)')
      IF (GRID) THEN
         WRITE (LUNCRT, '(4(I6, 2X, 3I5))') &
           (IB, X_IN(IB)%NI, X_IN(IB)%NJ, X_IN(IB)%NK, IB = 1, NBLOCK_IN)
      ELSE
         WRITE (LUNCRT, '(4(I6, 2X, 4I5))') &
           (IB, X_IN(IB)%MI, X_IN(IB)%MJ, X_IN(IB)%MK, NVARS, IB = 1, NBLOCK_IN)
      END IF

!     Transcribe the face dimensions, etc., to the output file arrays:

      ALLOCATE (X_OUT(NBLOCK_OUT))

      DO I = 1, NBLOCK_OUT

         IB = IBLOCK(I)
         IF = IFACE(I)
!!!!!!   ID = MOD (IF + 1, 2) ! Convert face # to dimension 1, 2, or 3
         ID = (IF + 1) / 2    ! Convert face # to dimension 1, 2, or 3
         I1 = MOD (ID, 3) + 1 ! 1st dimension of face IF
         I2 = MOD (I1, 3) + 1 ! 2nd ....................

!        I1 and I2 come out cyclic for ID = 1, 2, 3.
!        I.e., they are 2,3 and 3,1 and 1,2 respectively, meaning
!        nj,nk and nk,ni and ni,nj respectively.

         IF (ID /= 2) THEN ! I or K face
            SELECT CASE (I1)
               CASE (1)
                  IF (GRID) X_OUT(I)%NI = X_IN(IB)%NI
                  IF (FUNC) X_OUT(I)%MI = X_IN(IB)%MI
               CASE (2)
                  IF (GRID) X_OUT(I)%NI = X_IN(IB)%NJ
                  IF (FUNC) X_OUT(I)%MI = X_IN(IB)%MJ
               CASE (3)
                  IF (GRID) X_OUT(I)%NI = X_IN(IB)%NK
                  IF (FUNC) X_OUT(I)%MI = X_IN(IB)%MK
            END SELECT
            SELECT CASE (I2)
               CASE (1)
                  IF (GRID) X_OUT(I)%NJ = X_IN(IB)%NI
                  IF (FUNC) X_OUT(I)%MJ = X_IN(IB)%MI
               CASE (2)
                  IF (GRID) X_OUT(I)%NJ = X_IN(IB)%NJ
                  IF (FUNC) X_OUT(I)%MJ = X_IN(IB)%MJ
               CASE (3)
                  IF (GRID) X_OUT(I)%NJ = X_IN(IB)%NK
                  IF (FUNC) X_OUT(I)%MJ = X_IN(IB)%MK
            END SELECT
         ELSE ! J face - be consistent with UV_MAP's "extractHyperfaces"
            SELECT CASE (I2)
               CASE (1)
                  IF (GRID) X_OUT(I)%NI = X_IN(IB)%NI
                  IF (FUNC) X_OUT(I)%MI = X_IN(IB)%MI
               CASE (2)
                  IF (GRID) X_OUT(I)%NI = X_IN(IB)%NJ
                  IF (FUNC) X_OUT(I)%MI = X_IN(IB)%MJ
               CASE (3)
                  IF (GRID) X_OUT(I)%NI = X_IN(IB)%NK
                  IF (FUNC) X_OUT(I)%MI = X_IN(IB)%MK
            END SELECT
            SELECT CASE (I1)
               CASE (1)
                  IF (GRID) X_OUT(I)%NJ = X_IN(IB)%NI
                  IF (FUNC) X_OUT(I)%MJ = X_IN(IB)%MI
               CASE (2)
                  IF (GRID) X_OUT(I)%NJ = X_IN(IB)%NJ
                  IF (FUNC) X_OUT(I)%MJ = X_IN(IB)%MJ
               CASE (3)
                  IF (GRID) X_OUT(I)%NJ = X_IN(IB)%NK
                  IF (FUNC) X_OUT(I)%MJ = X_IN(IB)%MK
            END SELECT
         END IF

         IF (GRID) X_OUT(I)%NK = 1
         IF (FUNC) X_OUT(I)%MK = 1

      END DO

!     Write the output file header(s):

      IF (GRID) THEN

         CALL XYZ_HEADER_IO (2, LUNGOUT, FORMATTED_GOUT, NBLOCK_OUT, X_OUT,   &
                             IOS)
         IF (IOS /= 0) THEN
            WRITE (LUNCRT, '(/, A)') ' Trouble writing output grid header.'
            GO TO 999
         END IF

      END IF

      IF (FUNC) THEN

         CALL Q_HEADER_IO (2, LUNFOUT, FORMATTED_FOUT, NBLOCK_OUT, NVARS,     &
                           X_OUT, IOS)
         IF (IOS /= 0) THEN
            WRITE (LUNCRT, '(/, A)') ' Trouble writing output fun. file header.'
            GO TO 999
         END IF

      END IF

!     Extract the specified faces in the indicated order and write them:

      DO I = 1, NBLOCK_OUT

         IB = IBLOCK(I)

         IF (GRID) THEN

            CALL XYZ_ALLOCATE (X_OUT(I), IOS)

            CALL COPY_GRID_FACE (X_IN(IB), IFACE(I), X_OUT(I))

            NPTS = X_OUT(I)%NI * X_OUT(I)%NJ ! X_OUT(I)%NK = 1

            CALL XYZ_BLOCK_IO (2, LUNGOUT, FORMATTED_GOUT, NPTS,              &
                               X_OUT(I)%X, X_OUT(I)%Y, X_OUT(I)%Z, IOS)
            IF (IOS /= 0) THEN
               WRITE (LUNCRT, '(/, A, I5)') &
                  ' Trouble writing grid face.  Output block #:', I
               GO TO 999
            END IF

            DEALLOCATE (X_OUT(I)%X, X_OUT(I)%Y, X_OUT(I)%Z)

         END IF

         IF (FUNC) THEN

            CALL Q_ALLOCATE (X_OUT(I), NVARS, IOS)

            CALL COPY_FUNC_FACE (X_IN(IB), IFACE(I), X_OUT(I))

            CALL Q_BLOCK_IO (2, LUNFOUT, FORMATTED_FOUT, NVARS,               &
                             X_OUT(I)%MI, X_OUT(I)%MJ, X_OUT(I)%MK,           &
                             X_OUT(I)%Q, IOS)
            IF (IOS /= 0) THEN
               WRITE (LUNCRT, '(/, A, I5)') &
                  ' Trouble writing function face.  Output block #:', I
               GO TO 999
            END IF

            DEALLOCATE (X_OUT(I)%Q)

         END IF

      END DO

  999 CONTINUE
 
! *** STOP ! Avoid system dependencies.

      END PROGRAM GRID_FACES

!     --------------------------------------------------------------------------

      SUBROUTINE COPY_GRID_FACE (X_IN, IFACE, X_OUT)

!     Transcribe the indicated face of the input block to the output block.
!     Note that J faces are left-handed as in UV_MAP's "extractHyperfaces".

!     --------------------------------------------------------------------------

      USE GRID_BLOCK_STRUCTURE

      IMPLICIT NONE

!     Arguments:

      TYPE (GRID_TYPE), INTENT (IN)    :: X_IN  ! Input block
      INTEGER, INTENT (IN)             :: IFACE ! Block face to copy
      TYPE (GRID_TYPE), INTENT (INOUT) :: X_OUT ! Output block, dimensions set

!     Local variables:

      INTEGER :: NI_IN, NJ_IN, NK_IN

!     Execution:

      NI_IN = X_IN%NI;  NJ_IN = X_IN%NJ;  NK_IN = X_IN%NK

      SELECT CASE (IFACE)

         CASE (1)

            X_OUT%X(:,:,1) = X_IN%X (1, 1:NJ_IN, 1:NK_IN)
            X_OUT%Y(:,:,1) = X_IN%Y (1, 1:NJ_IN, 1:NK_IN)
            X_OUT%Z(:,:,1) = X_IN%Z (1, 1:NJ_IN, 1:NK_IN)

         CASE (2)

            X_OUT%X(:,:,1) = X_IN%X (NI_IN, 1:NJ_IN, 1:NK_IN)
            X_OUT%Y(:,:,1) = X_IN%Y (NI_IN, 1:NJ_IN, 1:NK_IN)
            X_OUT%Z(:,:,1) = X_IN%Z (NI_IN, 1:NJ_IN, 1:NK_IN)

         CASE (3)

            X_OUT%X(:,:,1) = X_IN%X (1:NI_IN, 1, 1:NK_IN)  ! Left-handed
            X_OUT%Y(:,:,1) = X_IN%Y (1:NI_IN, 1, 1:NK_IN)  ! as in UV_MAP
            X_OUT%Z(:,:,1) = X_IN%Z (1:NI_IN, 1, 1:NK_IN)

         CASE (4)

            X_OUT%X(:,:,1) = X_IN%X (1:NI_IN, NJ_IN, 1:NK_IN)  ! Likewise
            X_OUT%Y(:,:,1) = X_IN%Y (1:NI_IN, NJ_IN, 1:NK_IN)
            X_OUT%Z(:,:,1) = X_IN%Z (1:NI_IN, NJ_IN, 1:NK_IN)

         CASE (5)

            X_OUT%X(:,:,1) = X_IN%X (1:NI_IN, 1:NJ_IN, 1)
            X_OUT%Y(:,:,1) = X_IN%Y (1:NI_IN, 1:NJ_IN, 1)
            X_OUT%Z(:,:,1) = X_IN%Z (1:NI_IN, 1:NJ_IN, 1)

         CASE (6)

            X_OUT%X(:,:,1) = X_IN%X (1:NI_IN, 1:NJ_IN, NK_IN)
            X_OUT%Y(:,:,1) = X_IN%Y (1:NI_IN, 1:NJ_IN, NK_IN)
            X_OUT%Z(:,:,1) = X_IN%Z (1:NI_IN, 1:NJ_IN, NK_IN)

      END SELECT

      END SUBROUTINE COPY_GRID_FACE

!     --------------------------------------------------------------------------

      SUBROUTINE COPY_FUNC_FACE (X_IN, IFACE, X_OUT)

!     Transcribe the indicated face of the input block to the output block.
!     Note that J faces are left-handed as in UV_MAP's "extractHyperfaces".

!     --------------------------------------------------------------------------

      USE GRID_BLOCK_STRUCTURE

      IMPLICIT NONE

!     Arguments:

      TYPE (GRID_TYPE), INTENT (IN)    :: X_IN  ! Input block
      INTEGER, INTENT (IN)             :: IFACE ! Block face to copy
      TYPE (GRID_TYPE), INTENT (INOUT) :: X_OUT ! Output block, dimensions set

!     Local variables:

      INTEGER :: MI_IN, MJ_IN, MK_IN

!     Execution:

      MI_IN = X_IN%MI;  MJ_IN = X_IN%MJ;  MK_IN = X_IN%MK

      SELECT CASE (IFACE)

         CASE (1)

            X_OUT%Q(:,:,:,1) = X_IN%Q (:, 1, 1:MJ_IN, 1:MK_IN)

         CASE (2)

            X_OUT%Q(:,:,:,1) = X_IN%Q (:, MI_IN, 1:MJ_IN, 1:MK_IN)

         CASE (3)

            X_OUT%Q(:,:,:,1) = X_IN%Q (:, 1:MI_IN, 1, 1:MK_IN) ! Left-handed

         CASE (4)

            X_OUT%Q(:,:,:,1) = X_IN%Q (:, 1:MI_IN, MJ_IN, 1:MK_IN) ! Likewise

         CASE (5)

            X_OUT%Q(:,:,:,1) = X_IN%Q (:, 1:MI_IN, 1:MJ_IN, 1)

         CASE (6)

            X_OUT%Q(:,:,:,1) = X_IN%Q (:, 1:MI_IN, 1:MJ_IN, MK_IN)

      END SELECT

      END SUBROUTINE COPY_FUNC_FACE
