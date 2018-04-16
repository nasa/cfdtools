      SUBROUTINE D2UNLD (CMEM, IMEM, DMEM, IDE, IER)

C   Unlock the double precision space belonging to the data entity
C     pointed to by IDE.
C
C   DYNAMIC MEMORY INPUT-OUTPUT:
C     CMEM    Dynamically managed character data
C     IMEM    Dynamically managed integer data
C     DMEM    Dynamically managed double precision data
C   INPUT:
C     IDE     Pointer to data entity to be unlocked
C   OUTPUT:
C     IER     Returned error value.  Zero implies no errors.  Otherwise,
C                  1 =  Entity was already unlocked
C                 -1 =  Dynamic memory is corrupt or uninitialized.
C                 -2 =  Null pointer (IDE = 0)
C                 -3 =  Garbage value in pointer (points outside valid
C                       range in IMEM)
C                 -4 =  Either a deleted entity or a garbage value in
C                       pointer
C                 -5 =  Pointer refers to a deleted entity
C
C     CALLS:
C        D0PTR
C
C     HISTORY:
C        10Mar92  D. Parsons  created
C        13Apr92  D. Parsons  Add check of dynamic memory.
C*****
C
C   Long Name Alias:
C     ENTRY D2_UNLOCK_DMEM (CMEM, IMEM, DMEM, IDE, IER)

      EXTERNAL D2TSLD
      LOGICAL  D2TSLD

      CHARACTER CMEM*(*)
      INTEGER IMEM(*)
      DOUBLE PRECISION DMEM(*)
      INTEGER IDE, IER

      INTEGER IHDR, KHDR, NHDR, ITYP, LOWC, LOWI, LOWD, I

      CHARACTER*6 SUBNAM
      DATA SUBNAM /'D2UNLD'/

C     ------

C     Call the general pointer-check utility routine

      CALL D0PTR (CMEM, IMEM, DMEM, IDE, 0, 0, IHDR, ITYP, IER)
      IF (IER .NE. 0) GOTO 9900

C     Is entity already unlocked?

      IF (.NOT. D2TSLD (CMEM, IMEM, DMEM, IDE)) THEN
         IER = +1
         GOTO 9990
      ENDIF

C     Unlock double precision data space

      CALL D0UNLK (CMEM, IMEM, DMEM, IHDR, 3, LOWD)
      LOWC = IMEM(4)
      LOWI = IMEM(5)


C     Did the unlocking expose other, previously deleted entities to
C     at least partial garbage collection?

      NHDR = IMEM(IHDR+7)/256
      IF ((LOWD .LT. IMEM(6)) .AND. (NHDR .NE. 0)) THEN
        DO 110 I=1,1+(IMEM(2)-IMEM(16))/8
          KHDR = NHDR
          NHDR = IMEM(KHDR+7)/256
          IF (IMEM(KHDR) .LT. 0) CALL D0ERAS (CMEM, IMEM, DMEM, KHDR)
          IF (NHDR .EQ. 0) GO TO 120
          IF (IMEM(NHDR+1) .LT. LOWC .AND. IMEM(NHDR+2) .LT. LOWI
     +        .AND. IMEM(NHDR+3) .LT. LOWD) GO TO 120
  110   CONTINUE

C       The DO loop bound is chosen so that it CANNOT terminate here.
        STOP 'Dynamic memory fatally corrupted - D2UNLD'

  120   CONTINUE
      END IF

      IER = 0

      RETURN

C   Error reporting section

 9900 CALL DTERR (1, SUBNAM, IER, 0)
      RETURN

 9990 CALL DTERR (0, SUBNAM, IER, 0)
      RETURN

      END
