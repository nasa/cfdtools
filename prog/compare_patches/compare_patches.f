C***********************************************************************
C
      PROGRAM COMPARE_PATCHES
C
C     COMPARE_PATCHES reads 2 PLOT3D surface grids and locates the
C     largest point difference in the specified subscript range(s).
C     The grids are assumed to be /3D /MGRID /FORMATTED or /UNF.
C
C     04/11/98  D.A.Saunders   Adaptation of single-block COMPARE_GRIDS,
C                              for the AEROSURF repeat checking.
C     02/06/99      "          Cases with more than 50 blocks can't be
C                              handled by PLOT3D, so use this program to
C                              display actual coordinates 1 pt. at time.
C     05/05/99  Justin Jagger  Dynamic workspace version.
C     11/26/99  David Saunders Unformatted option, and 2D/blanked options
C                              for UV_MAP outputs.
C     12/15/99    "     "      Display mean deviation as well as max.
C     10/26/11    "     "      The 2D formatted case had a read error.
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

C     Allocatable Variables:

      INTEGER, ALLOCATABLE, DIMENSION (:) ::
     >   NI, NJ, NI2, NJ2

      INTEGER, ALLOCATABLE, DIMENSION (:,:,:) ::
     >   IBLANK

      REAL, ALLOCATABLE, DIMENSION (:,:,:) ::
     >   X1, Y1, Z1, X2, Y2, Z2

C     Local variables:

      INTEGER
     >   I, J, K, I1, I2, J1, J2, K1, K2, IMAX, JMAX, KMAX,
     >   NCOMPARE, NPATCH, IDIM, JDIM
      REAL
     >   D, DMAX, DRMS, DX, DY, DZ, DXMAX, DYMAX, DZMAX
      LOGICAL
     >   BLANKED, FORMATTED, PRINT, TWOD
      CHARACTER
     >   ANSWER * 1, FILENAME * 48

C     Execution:

      WRITE (LUNCRT, '(A)', ADVANCE='NO')
     >   ' Formatted or unformatted? [f/u]: '
      READ  (LUNKBD, '(A)') ANSWER
      FORMATTED = ANSWER == 'f' .OR. ANSWER == 'F'

      WRITE (LUNCRT, '(A)', ADVANCE='NO')
     >   ' 2D or 3D surface patches? [2/3]: '
      READ  (LUNKBD, *) K
      TWOD = K == 2

      IF (TWOD .AND. .NOT. FORMATTED) THEN
         WRITE (LUNCRT, '(A)', ADVANCE='NO')
     >      ' Blanked or unblanked? [b/u]: '
         READ  (LUNKBD, '(A)') ANSWER
         BLANKED = ANSWER == 'b' .OR. ANSWER == 'B'
      END IF

      WRITE (LUNCRT, '(A)', ADVANCE='NO') ' 1st patch file: '
      READ  (LUNKBD, '(A)') FILENAME

      IF (FORMATTED) THEN
         OPEN (UNIT=LUN1, FILE=FILENAME, STATUS='OLD', FORM='FORMATTED')
      ELSE
         OPEN (UNIT=LUN1, FILE=FILENAME, STATUS='OLD',
     >         FORM='UNFORMATTED')
      END IF

      WRITE (LUNCRT, '(A)', ADVANCE='NO') ' 2nd patch file: '
      READ  (LUNKBD, '(A)') FILENAME

      IF (FORMATTED) THEN
         OPEN (UNIT=LUN2, FILE=FILENAME, STATUS='OLD', FORM='FORMATTED')
      ELSE
         OPEN (UNIT=LUN2, FILE=FILENAME, STATUS='OLD',
     >         FORM='UNFORMATTED')
      END IF

      IF (FORMATTED) THEN
         READ (LUN1, *) NPATCH
      ELSE
         READ (LUN1)    NPATCH
      END IF

      ALLOCATE (NI(NPATCH), NJ(NPATCH), NI2(NPATCH), NJ2(NPATCH))

      IF (FORMATTED) THEN
         IF (TWOD) THEN
            READ (LUN1, *) (NI(K), NJ(K),    K = 1, NPATCH)
         ELSE
            READ (LUN1, *) (NI(K), NJ(K), I, K = 1, NPATCH)
         END IF
      ELSE
         IF (TWOD) THEN
            READ (LUN1) (NI(K), NJ(K),    K = 1, NPATCH)
         ELSE
            READ (LUN1) (NI(K), NJ(K), I, K = 1, NPATCH)
         END IF
      END IF

      WRITE (LUNCRT, '(/, A)') ' Grid 1 patch dimensions: '
      WRITE (LUNCRT, '(1X, I4, 2I7)') (K, NI(K), NJ(K), K = 1, NPATCH)

C     We're avoiding packing the blocks, so make all patches
C     as big as the biggest one.

      IDIM = 0
      JDIM = 0
      DO K = 1, NPATCH
         IDIM = MAX (IDIM, NI(K))
         JDIM = MAX (JDIM, NJ(K))
      END DO

      ALLOCATE (X1(IDIM, JDIM, NPATCH), Y1(IDIM, JDIM, NPATCH),
     >          Z1(IDIM, JDIM, NPATCH), X2(IDIM, JDIM, NPATCH),
     >          Y2(IDIM, JDIM, NPATCH), Z2(IDIM, JDIM, NPATCH))

      IF (BLANKED) ALLOCATE (IBLANK(IDIM, JDIM, NPATCH))
      
      DO K = 1, NPATCH
         IF (FORMATTED) THEN
            IF (TWOD) THEN
               READ (LUN1, *)
     >            ((X1(I,J,K), I = 1, NI(K)), J = 1, NJ(K)),
     >            ((Y1(I,J,K), I = 1, NI(K)), J = 1, NJ(K))
            ELSE
               READ (LUN1, *)
     >            ((X1(I,J,K), I = 1, NI(K)), J = 1, NJ(K)),
     >            ((Y1(I,J,K), I = 1, NI(K)), J = 1, NJ(K)),
     >            ((Z1(I,J,K), I = 1, NI(K)), J = 1, NJ(K))
            END IF
         ELSE
            IF (TWOD) THEN
               IF (BLANKED) THEN
                  READ (LUN1)
     >               ((X1(I,J,K), I = 1, NI(K)), J = 1, NJ(K)),
     >               ((Y1(I,J,K), I = 1, NI(K)), J = 1, NJ(K)),
     >               ((IBLANK(I,J,K), I = 1, NI(K)), J = 1, NJ(K))
               ELSE
                  READ (LUN1)
     >               ((X1(I,J,K), I = 1, NI(K)), J = 1, NJ(K)),
     >               ((Y1(I,J,K), I = 1, NI(K)), J = 1, NJ(K))
               END IF
            ELSE
               READ (LUN1)
     >            ((X1(I,J,K), I = 1, NI(K)), J = 1, NJ(K)),
     >            ((Y1(I,J,K), I = 1, NI(K)), J = 1, NJ(K)),
     >            ((Z1(I,J,K), I = 1, NI(K)), J = 1, NJ(K))
            END IF
         END IF
      END DO

      CLOSE (LUN1)

      IF (FORMATTED) THEN
         READ (LUN2, *) K
      ELSE
         READ (LUN2   ) K
      END IF

      IF (K /= NPATCH) THEN
         WRITE (LUNCRT, '(A, /, I8)') ' # blocks in grid 2:', K
         GO TO 999
      END IF

      IF (FORMATTED) THEN
         IF (TWOD) THEN
            READ (LUN2, *) (NI2(K), NJ2(K),    K = 1, NPATCH)
         ELSE
            READ (LUN2, *) (NI2(K), NJ2(K), I, K = 1, NPATCH)
         END IF
      ELSE
         IF (TWOD) THEN
            READ (LUN2) (NI2(K), NJ2(K),    K = 1, NPATCH)
         ELSE
            READ (LUN2) (NI2(K), NJ2(K), I, K = 1, NPATCH)
         END IF
      END IF

      DO K = 1, NPATCH
         IF (NI2(K) /= NI(K) .OR. NJ2(K) /= NJ(K)) THEN
            WRITE (LUNCRT, '(A, I2, (/, 2I8))')
     >         ' Mismatched dimensions for patch # ', K,
     >         NI(K), NJ(K), NI2(K), NJ2(K)
            GO TO 999
         END IF

         IF (FORMATTED) THEN
            IF (TWOD) THEN
               READ (LUN2, *)
     >            ((X2(I,J,K), I = 1, NI(K)), J = 1, NJ(K)),
     >            ((Y2(I,J,K), I = 1, NI(K)), J = 1, NJ(K))
            ELSE
               READ (LUN2, *)
     >            ((X2(I,J,K), I = 1, NI(K)), J = 1, NJ(K)),
     >            ((Y2(I,J,K), I = 1, NI(K)), J = 1, NJ(K)),
     >            ((Z2(I,J,K), I = 1, NI(K)), J = 1, NJ(K))
            END IF
         ELSE
            IF (TWOD) THEN
               IF (BLANKED) THEN
                  READ (LUN2)
     >               ((X2(I,J,K), I = 1, NI(K)), J = 1, NJ(K)),
     >               ((Y2(I,J,K), I = 1, NI(K)), J = 1, NJ(K)),
     >               ((IBLANK(I,J,K), I = 1, NI(K)), J = 1, NJ(K))
               ELSE
                  READ (LUN2)
     >               ((X2(I,J,K), I = 1, NI(K)), J = 1, NJ(K)),
     >               ((Y2(I,J,K), I = 1, NI(K)), J = 1, NJ(K))
               END IF
            ELSE
               READ (LUN2)
     >            ((X2(I,J,K), I = 1, NI(K)), J = 1, NJ(K)),
     >            ((Y2(I,J,K), I = 1, NI(K)), J = 1, NJ(K)),
     >            ((Z2(I,J,K), I = 1, NI(K)), J = 1, NJ(K))
            END IF
         END IF

      END DO

      CLOSE (LUN2)

      IF (TWOD) THEN
         Z1 = ZERO
         Z2 = ZERO
      END IF

      OPEN (UNIT=LUNOUT, FILE='compare.out', STATUS='UNKNOWN',
     >      FORM='FORMATTED')

C     Loop over patches to compare

  100 CONTINUE

         WRITE (LUNCRT, '(A)', ADVANCE='NO')
     >      ' First and last patches to be compared: '

         READ  (LUNKBD, *, END=999) K1, K2
         IF (K1 < 1 .OR. K1 > NPATCH) GO TO 100
         IF (K2 < 1 .OR. K2 > NPATCH) GO TO 100

         IF (K1 == K2) THEN
            K = K1
            WRITE (LUNCRT, '(A, 2I5)') ' Dimensions: ', NI(K), NJ(K)
  102       WRITE (LUNCRT, '(A)', ADVANCE='NO')
     >         ' (I,J) range to compare (2 pairs): '

            READ  (LUNKBD, *, END=999) I1, I2, J1, J2
            IF (I1 < 1 .OR. I1 > NI(K)) GO TO 102
            IF (I2 < 1 .OR. I2 > NI(K)) GO TO 102
            IF (J1 < 1 .OR. J1 > NJ(K)) GO TO 102
            IF (J2 < 1 .OR. J2 > NJ(K)) GO TO 102

            NCOMPARE = (I2 - I1 + 1) * (J2 - J1 + 1)
            PRINT = NCOMPARE <= 1000

            IF (NCOMPARE == 1) THEN
               WRITE (LUNCRT, '(A, 3F12.6)')
     >            ' X/Y/Z1: ', X1(I1,J1,K), Y1(I1,J1,K), Z1(I1,J1,K),
     >            ' X/Y/Z2: ', X2(I1,J1,K), Y2(I1,J1,K), Z2(I1,J1,K)
            END IF

         ELSE
            PRINT = .FALSE.
         END IF

         IF (PRINT) WRITE (LUNOUT, '(//, A, /)')
     >      '    #    I    J          dX          dY          dZ'

         DMAX = -1.
         DRMS = ZERO
         NCOMPARE = 0

         DO K = K1, K2

            IF (K1 /= K2) THEN
               I1 = 1
               I2 = NI(K)
               J1 = 1
               J2 = NJ(K)
            END IF

            DO J = J1, J2
               DO I = I1, I2
                  DX = X2(I,J,K) - X1(I,J,K)
                  DY = Y2(I,J,K) - Y1(I,J,K)
                  DZ = Z2(I,J,K) - Z1(I,J,K)
                  IF (PRINT)
     >               WRITE (LUNOUT, '(3I5,3F12.7)') K, I, J, DX, DY, DZ
                  D    = DX**2 + DY**2 + DZ**2
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
                  END IF
               END DO
            END DO
         END DO

         DMAX = SQRT (DMAX)
         DRMS = SQRT (DRMS / REAL (NCOMPARE))

         WRITE (LUNCRT, '(A,1P,E12.5,0P,A,2I4,A,I3,A,3F12.6)')
     >      ' Max. difference:', DMAX, '  (i,j): ', IMAX, JMAX,
     >      '  patch:', KMAX, '  dX/Y/Z:', DXMAX, DYMAX, DZMAX,
     >      ' RMS  difference:', DRMS
      GO TO 100

  999 CONTINUE

      END PROGRAM COMPARE_PATCHES
