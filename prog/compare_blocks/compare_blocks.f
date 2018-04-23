C***********************************************************************
C
      PROGRAM COMPARE_BLOCKS
C
C     COMPARE_BLOCKS reads two PLOT3D volume grids and locates the
C     largest point difference in the specified subscript range(s).
C     The grids are assumed to be /3D, /MGRID and /FORMATTED or /IEEE_D.
C
C     07/21/03  D.A.Saunders   Adaptation of COMPARE_PATCHES.
C
C     David Saunders, ELORET/NASA Ames Research Center, Moffett Fld., CA
C
C***********************************************************************

      IMPLICIT NONE

C     Constants:

      INTEGER, PARAMETER ::
     >   LUN1    = 1,
     >   LUN2    = 2,
     >   LUNCRT  = 6,
     >   LUNKBD  = 5,
     >   LUNOUT  = 3

      REAL, PARAMETER ::
     >   ZERO = 0.

      CHARACTER, PARAMETER ::
     >   FORM * 11 = 'UNFORMATTED'

C     Allocatable variables:

      INTEGER, ALLOCATABLE, DIMENSION (:) ::
     >   NI, NJ, NK, NI2, NJ2, NK2

      INTEGER, ALLOCATABLE, DIMENSION (:) ::
     >   IBLANK, IP

      REAL, ALLOCATABLE, DIMENSION (:) ::
     >   X1, Y1, Z1, X2, Y2, Z2

C     Local variables:

      INTEGER
     >   I, J, K, L, I1, I2, J1, J2, K1, K2, L1, L2,
     >   IOSTAT, IMAX, JMAX, KMAX, LMAX,
     >   N, NBLANK, NBLOCKS, NCOMPARE, NPOINTS
      REAL
     >   D, DMAX, DRMS, DX, DY, DZ, DXMAX, DYMAX, DZMAX
      LOGICAL
     >   BLANKED, FORMATTED, PRINT
      CHARACTER
     >   ANSWER * 1, FILENAME * 64

C     Execution:

      WRITE (LUNCRT, '(A)', ADVANCE='NO')
     >   ' Formatted or unformatted? [f/u]: '
      READ  (LUNKBD, '(A)') ANSWER
      FORMATTED = ANSWER == 'f' .OR. ANSWER == 'F'

      WRITE (LUNCRT, '(A)', ADVANCE='NO')
     >   ' Blanked or unblanked?     [b/u]: '
      READ  (LUNKBD, '(A)') ANSWER
      BLANKED   = ANSWER == 'b' .OR. ANSWER == 'B'

      WRITE (LUNCRT, '(A)', ADVANCE='NO') ' 1st volume grid file: '
      READ  (LUNKBD, '(A)') FILENAME

      K1 = 1
      IF (FORMATTED) K1 = 3

      OPEN (UNIT=LUN1, FILE=FILENAME, STATUS='OLD', FORM=FORM(K1:11))

      WRITE (LUNCRT, '(A)', ADVANCE='NO') ' 2nd volume grid file: '
      READ  (LUNKBD, '(A)') FILENAME

      IF (FORMATTED) THEN
         READ (LUN1, *, IOSTAT=IOSTAT) NBLOCKS
      ELSE
         READ (LUN1,    IOSTAT=IOSTAT) NBLOCKS
      END IF

      IF (IOSTAT /= 0) THEN
         WRITE (LUNCRT, '(/, A, I6)')
     >      ' Error reading NBLOCKS for grid 1.  IOSTAT: ', IOSTAT
         GO TO 999
      END IF

      ALLOCATE (NI(NBLOCKS),  NJ(NBLOCKS),  NK(NBLOCKS),
     >          NI2(NBLOCKS), NJ2(NBLOCKS), NK2(NBLOCKS))

      IF (FORMATTED) THEN
         READ (LUN1, *, IOSTAT=IOSTAT)
     >      (NI(L), NJ(L), NK(L), L = 1, NBLOCKS)
      ELSE
         READ (LUN1,    IOSTAT=IOSTAT)
     >      (NI(L), NJ(L), NK(L), L = 1, NBLOCKS)
      END IF

      IF (IOSTAT /= 0) THEN
         WRITE (LUNCRT, '(/, A, I6)')
     >      ' Error reading block dimensions, grid 1.  IOSTAT: ', IOSTAT
         GO TO 999
      END IF

      WRITE (LUNCRT, '(/, A)') ' Grid 1 block dimensions: '
      WRITE (LUNCRT, '(1X, I4, 3I7)')
     >   (L, NI(L), NJ(L), NK(L), L = 1, NBLOCKS)

C     Pack the blocks in memory.  Determine the start of each block:

      ALLOCATE (IP(NBLOCKS + 1))

      IP(1) = 1
      DO L = 1, NBLOCKS
         NPOINTS = NI(L) * NJ(L) * NK(L)
         IP(L+1) = IP(L) + NPOINTS
      END DO

      NPOINTS = IP(NBLOCKS+1) - 1

      ALLOCATE (X1(NPOINTS), Y1(NPOINTS), Z1(NPOINTS),
     >          X2(NPOINTS), Y2(NPOINTS), Z2(NPOINTS), IBLANK(NPOINTS))

      NBLANK = 0
      IF (BLANKED) NBLANK = NPOINTS

      DO L = 1, NBLOCKS

         I1 = IP(L)
         NPOINTS = IP(L+1) - I1

         CALL RD_BLK (LUN1, FORMATTED, NPOINTS, NBLANK,
     >                X1(I1), Y1(I1), Z1(I1), IBLANK(I1))

         IF (IOSTAT /= 0) THEN
            WRITE (LUNCRT, '(/, A, I5, A, I5)')
     >         ' Trouble reading grid 1.  Block #:', L,
     >         '.  IOSTAT:', IOSTAT
            GO TO 999
         END IF

      END DO

      CLOSE (LUN1)

      OPEN (UNIT=LUN2, FILE=FILENAME, STATUS='OLD', FORM=FORM(K1:11))

      IF (FORMATTED) THEN
         READ (LUN2, *, IOSTAT=IOSTAT) L
      ELSE
         READ (LUN2,    IOSTAT=IOSTAT) L
      END IF

      IF (IOSTAT /= 0) THEN
         WRITE (LUNCRT, '(/, A, I6)')
     >      ' Error reading NBLOCKS for grid 2.  IOSTAT: ', IOSTAT
         GO TO 999
      END IF

      IF (L /= NBLOCKS) THEN
         WRITE (LUNCRT, '(A, /, I8)') ' # blocks in grid 2:', L
         GO TO 999
      END IF

      IF (FORMATTED) THEN
         READ (LUN2, *, IOSTAT=IOSTAT)
     >      (NI2(L), NJ2(L), NK2(L), L = 1, NBLOCKS)
      ELSE
         READ (LUN2,    IOSTAT=IOSTAT)
     >      (NI2(L), NJ2(L), NK2(L), L = 1, NBLOCKS)
      END IF

      IF (IOSTAT /= 0) THEN
         WRITE (LUNCRT, '(/, A, I6)')
     >      ' Error reading block sizes for grid 2.  IOSTAT: ', IOSTAT
         GO TO 999
      END IF

      DO L = 1, NBLOCKS

         IF (NI2(L) /= NI(L) .OR. NJ2(L) /= NJ(L) .OR.
     >       NK2(L) /= NK(L)) THEN
            WRITE (LUNCRT, '(A, I3, A, (/, 3I8))')
     >         ' Mismatched dimensions for block # ', L, ':',
     >         NI(L), NJ(L), NK(L), NI2(L), NJ2(L), NK2(L) 
            GO TO 999
         END IF

         I1 = IP(L)
         NPOINTS = IP(L+1) - I1

         CALL RD_BLK (LUN2, FORMATTED, NPOINTS, NBLANK,
     >                X2(I1), Y2(I1), Z2(I1), IBLANK(I1))

         IF (IOSTAT /= 0) THEN
            WRITE (LUNCRT, '(/, A, I5, A, I5)')
     >         ' Trouble reading grid 2.  Block #:', L,
     >         '.  IOSTAT:', IOSTAT
            GO TO 999
         END IF

      END DO

      CLOSE (LUN2)

      OPEN (UNIT=LUNOUT, FILE='compare.out', STATUS='UNKNOWN',
     >      FORM='FORMATTED')

C     Loop over blocks to compare:

  100 CONTINUE

         WRITE (LUNCRT, '(A)', ADVANCE='NO')
     >      ' First and last blocks to be compared: '

         READ  (LUNKBD, *, END=999) L1, L2
         IF (L1 < 1 .OR. L1 > NBLOCKS) GO TO 100
         IF (L2 < 1 .OR. L2 > NBLOCKS) GO TO 100

         IF (L1 == L2) THEN
            L = L1
            WRITE (LUNCRT, '(A, 3I5)')
     >         ' Dimensions: ', NI(L), NJ(L), NK(L)
  102       WRITE (LUNCRT, '(A)', ADVANCE='NO')
     >         ' (I,J,K) range to compare (3 pairs): '

            READ  (LUNKBD, *, END=999) I1, I2, J1, J2, K1, K2
            IF (I1 < 1 .OR. I1 > NI(L)) GO TO 102
            IF (I2 < 1 .OR. I2 > NI(L)) GO TO 102
            IF (J1 < 1 .OR. J1 > NJ(L)) GO TO 102
            IF (J2 < 1 .OR. J2 > NJ(L)) GO TO 102
            IF (K1 < 1 .OR. K1 > NK(L)) GO TO 102
            IF (K2 < 1 .OR. K2 > NK(L)) GO TO 102

            NCOMPARE = (I2 - I1 + 1) * (J2 - J1 + 1) * (K2 - K1 + 1)
            PRINT = NCOMPARE <= 1000

            IF (NCOMPARE == 1) THEN
               N = IP(L) + I1 - 1 + (J1 - 1) * NI(L) +
     >             (K1 - 1) * NI(L) * NJ(L)
               WRITE (LUNCRT, '(A, 3F12.6)')
     >            ' X/Y/Z1: ', X1(N), Y1(N), Z1(N),
     >            ' X/Y/Z2: ', X2(N), Y2(N), Z2(N)
            END IF

         ELSE
            PRINT = .FALSE.
         END IF

         IF (PRINT) WRITE (LUNOUT, '(//, A, /)')
     >      '    #    I    J    K          dX          dY          dZ'

         DMAX = -1.
         DRMS = ZERO
         NCOMPARE = 0

         DO L = L1, L2

            IF (L1 /= L2) THEN
               I1 = 1
               I2 = NI(L)
               J1 = 1
               J2 = NJ(L)
               K1 = 1
               K2 = NK(L)
            END IF

            DO K = K1, K2
               DO J = J1, J2
                  N = IP(L) + I1 - 1 + (J - 1) * NI(L) +
     >                        (K - 1)  * NI(L) * NJ(L)
                  DO I = I1, I2
                     DX = X2(N) - X1(N)
                     DY = Y2(N) - Y1(N)
                     DZ = Z2(N) - Z1(N)
                     N  = N + 1

                     IF (PRINT) WRITE (LUNOUT, '(4I5,3F12.7)')
     >                  L, I, J, K, DX, DY, DZ

                     D = DX**2 + DY**2 + DZ**2
                     DRMS = DRMS + D
                     NCOMPARE = NCOMPARE + 1

                     IF (D > DMAX) THEN
                        DMAX  = D
                        DXMAX = DX
                        DYMAX = DY
                        DZMAX = DZ
                        IMAX  = I
                        JMAX  = J
                        KMAX  = K
                        LMAX  = L
                     END IF
                  END DO
               END DO
            END DO
         END DO

         DMAX = SQRT (DMAX)
         DRMS = SQRT (DRMS / REAL (NCOMPARE))

         WRITE (LUNCRT, '(A,1P,E12.5,0P,A,3I4,A,I4,A,3F12.6)')
     >      ' Max. difference:', DMAX, '  (i,j,k): ', IMAX, JMAX, KMAX,
     >      '  block:', LMAX, '  dX/Y/Z:', DXMAX, DYMAX, DZMAX,
     >      ' RMS  difference:', DRMS
      GO TO 100

  999 CONTINUE

      CONTAINS ! Internal procedure for COMPARE_BLOCKS

!       ---------------------------------------------------------------
        subroutine RD_BLK (lun, formatted, npoints, nblank,
     .                     x, y, z, iblank)
!       ---------------------------------------------------------------

!       RD_BLK reads one block of [un]blanked data in PLOT3D order on disk,
!       formatted or not, efficiently by avoiding subscripting.
!       It is an adaptation of the IO_R_P3D utility from CFD_IO_PACKAGE.

!       Arguments:

        integer, intent (in)  :: lun                ! Logical unit to read from
        logical, intent (in)  :: formatted          ! T or F
        integer, intent (in)  :: npoints,           ! # points in this block
     .                           nblank             ! # iblanks (possibly 0)
        real,    intent (out), dimension (npoints) ::
     .                           x, y, z            ! Packed block (x,y,z)s ...
        integer, intent (out) :: iblank(nblank)     ! ... and its blanking

!       Execution:

        iostat = 1

        if (formatted) then
          if (nblank == 0) then
            read (lun, *, err=999, iostat=iostat) x, y, z
          else
            read (lun, *, err=999, iostat=iostat) x, y, z, iblank
          end if
        else
          if (nblank == 0) then
            read (lun,    err=999, iostat=iostat) x, y, z
          else
            read (lun,    err=999, iostat=iostat) x, y, z, iblank
          end if
        end if

        iostat = 0

  999   return

        end subroutine RD_BLK

      END PROGRAM COMPARE_BLOCKS
