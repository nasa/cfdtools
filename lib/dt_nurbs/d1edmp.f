      SUBROUTINE D1EDMP (CMEM, IMEM, DMEM, MXC, MXI, MXD, LIM, LOWHDR,
     +    LUNIT, IHDR)

C   Dump one entity, assuming as little as possible beyond the
C   correctness of the input parameters.
C
C   DYNAMIC MEMORY INPUT-OUTPUT:
C     CMEM    Dynamically managed character data
C     IMEM    Dynamically managed integer data
C     DMEM    Dynamically managed double precision data
C   INPUT:
C     MXC     Correct length of CMEM
C     MXI     Correct length of IMEM
C     MXD     Correct length of DMEM
C     LIM     Maximum number of lines to print of one data space
C     LOWHDR  Lowest header index
C     LUNIT   Logical unit used for dump file (assumed already open)
C     IHDR    Index of header block for data entity to dump
C   OUTPUT:
C     Additional lines to dump file on logical unit LUNIT
C
C   CALLS:    D0CDMP, D0IDMP
C   HISTORY:
C     10Feb92 P Kraushar    Created.
C     26Mar92 D. Parsons    Make the minimum number of lines=3
C*****
C
      CHARACTER CMEM*(*)
      INTEGER IMEM(*)
      DOUBLE PRECISION DMEM(*)
      INTEGER MXC, MXI, MXD, LIM, LIMX, LOWHDR, LUNIT, IHDR

      INTEGER I1, I2, KREL, KABS, KEND, KLEN, LCNT, LDO, I, J
      LOGICAL ESKIP, MSKIP, SPECAL
      CHARACTER CLASS*9, OUTINT(5)*12, OUTSTR*56, OUTCTL*56
      INTEGER MEMIN(3)

      DATA MEMIN / 2, 20, 2 /
      DATA OUTSTR, OUTCTL /2*
     +    '|          |          |          |          |          |'/
C****
      LIMX = LIM
      IF (LIMX .LE. 0) LIMX = 20
      LIMX = MAX (LIMX,3)
      IF (IHDR .LT. 1 .OR. IHDR .GT. MXI-8) THEN
        WRITE (LUNIT, '(A)') ' Header block location not all in IMEM'
        RETURN
      END IF
      IF (IHDR .LT. LOWHDR .OR. MOD(MXI-IHDR,8) .NE. 0)
     +    WRITE (LUNIT, '(A)') ' Header block location not valid'
C   Classify header block
      IF (IMEM(IHDR) .GT. 0) THEN
        CLASS = 'ACTIVE'
      ELSE IF (IMEM(IHDR+4) .NE. 0 .OR. IMEM(IHDR+5) .NE. 0
     +    .OR. IMEM(IHDR+6) .NE. 0) THEN
        CLASS = 'DELETED'
      ELSE
        CLASS = 'FREE HDR'
      END IF
C   Dump header block to dump file
      I1 = IMEM(IHDR)/256
      I2 = MOD(ABS(IMEM(IHDR)),256)
      WRITE (LUNIT,  5) IHDR
      WRITE (LUNIT, 10) IMEM(IHDR), I1, I2, CLASS
      WRITE (LUNIT, 21) IMEM(IHDR+1)
      WRITE (LUNIT, 22) IMEM(IHDR+2)
      WRITE (LUNIT, 23) IMEM(IHDR+3)
      WRITE (LUNIT, 24) IMEM(IHDR+4)
      WRITE (LUNIT, 25) IMEM(IHDR+5)
      WRITE (LUNIT, 26) IMEM(IHDR+6)
      I1 = IMEM(IHDR+7)/256
      I2 = MOD(ABS(IMEM(IHDR+7)),256)
      WRITE (LUNIT, 30) IMEM(IHDR+7), I1, I2

C   Dump character data space to file
      IF (IMEM(IHDR+4) .NE. 0) THEN
        KABS = IMEM(IHDR+1)
        KLEN = ABS(IMEM(IHDR+4))
        IF (KABS .GE. MXC .OR. KABS+KLEN .LE. MEMIN(1)) THEN
          WRITE (LUNIT, '(A)') ' CMEM data space totally out of range'
        ELSE
          IF (KABS .LT. MEMIN(1)) THEN
            WRITE (LUNIT, '(A)') ' CMEM data space begins out of range'
            KREL = MEMIN(1) - KABS
            KABS = KABS + KREL
            KLEN = KLEN - KREL
            KREL = 1 + KREL
          ELSE
            WRITE (LUNIT, '(A)') ' CMEM'
            KREL = 1
          END IF
          ESKIP = (KABS + KLEN .GT. MXC)
          IF (ESKIP) KLEN = MXC - KABS
          KEND = KABS + KLEN - 1
          LCNT = 1 + (KLEN - 1)/50
          MSKIP = (LCNT .GT. LIMX)
          LDO = MIN (LCNT - 2, LIMX - 2)
   41     FORMAT(13X,A56)
          DO 100 I=1,LDO
            CALL D0CDMP (CMEM(KABS:KABS+49), OUTSTR, OUTCTL, SPECAL)
            WRITE (LUNIT, 40) KABS, KREL, OUTSTR
            IF (SPECAL) WRITE (LUNIT, 41) OUTCTL
            KABS = KABS + 50
            KREL = KREL + 50
  100     CONTINUE
          IF (MSKIP) THEN
            WRITE (LUNIT, '(A)') ' ... middle section skipped'
            KABS = KABS + 50*(LCNT-LIMX)
            KREL = KREL + 50*(LCNT-LIMX)
          END IF
          IF (LCNT .GT. 1) THEN
            CALL D0CDMP (CMEM(KABS:KABS+49), OUTSTR, OUTCTL, SPECAL)
            WRITE (LUNIT, 40) KABS, KREL, OUTSTR
            IF (SPECAL) WRITE (LUNIT, 41) OUTCTL
            KABS = KABS + 50
            KREL = KREL + 50
          END IF
          CALL D0CDMP (CMEM(KABS:KEND), OUTSTR, OUTCTL, SPECAL)
          WRITE (LUNIT, 40) KABS, KREL, OUTSTR
          IF (SPECAL) WRITE (LUNIT, 41) OUTCTL
          IF (ESKIP) WRITE (LUNIT, '(A)')
     +         ' ... data space ends out of range'
        END IF
      END IF
C   Dump integer data space to file
      IF (IMEM(IHDR+5) .NE. 0) THEN
        KABS = IMEM(IHDR+2)
        KLEN = ABS(IMEM(IHDR+5))
        IF (KABS .GE. MXI .OR. KABS+KLEN .LE. MEMIN(2)) THEN
          WRITE (LUNIT, '(A)') ' IMEM data space totally out of range'
        ELSE
          IF (KABS .LT. MEMIN(2)) THEN
            WRITE (LUNIT, '(A)') ' IMEM data space begins out of range'
            KREL = MEMIN(2) - KABS
            KABS = KABS + KREL
            KLEN = KLEN - KREL
            KREL = 1 + KREL
          ELSE
            WRITE (LUNIT, '(A)') ' IMEM'
            KREL = 1
          END IF
          ESKIP = (KABS + KLEN .GT. MXI)
          IF (ESKIP) KLEN = MXI - KABS
          KEND = KABS + KLEN - 1
          LCNT = 1 + (KLEN - 1)/5
          MSKIP = (LCNT .GT. LIMX)
          LDO = MIN (LCNT - 2, LIMX - 2)
          DO 210 I=1,LDO
            DO 200 J=1,5
              CALL D0IDMP (IMEM(KABS+J-1), LOWHDR, MXI, OUTINT(J))
  200       CONTINUE
            WRITE (LUNIT, 50) KABS, KREL, OUTINT
            KABS = KABS + 5
            KREL = KREL + 5
  210     CONTINUE
          IF (MSKIP) THEN
            WRITE (LUNIT, '(A)') ' ... middle section skipped'
            KABS = KABS + 5*(LCNT-LIMX)
            KREL = KREL + 5*(LCNT-LIMX)
          END IF
          IF (LCNT .GT. 1) THEN
            DO 220 J=1,5
              CALL D0IDMP (IMEM(KABS+J-1), LOWHDR, MXI, OUTINT(J))
  220       CONTINUE
            WRITE (LUNIT, 50) KABS, KREL, OUTINT
            KABS = KABS + 5
            KREL = KREL + 5
          END IF
          LDO = KEND - KABS + 1
          DO 230 J=1,LDO
            CALL D0IDMP (IMEM(KABS+J-1), LOWHDR, MXI, OUTINT(J))
  230     CONTINUE
          WRITE (LUNIT, 50) KABS, KREL, (OUTINT(J), J=1,LDO)
          IF (ESKIP) WRITE (LUNIT, '(A)')
     +         ' ... data space ends out of range'
        END IF
      END IF
C   Dump double precision data space to file
      IF (IMEM(IHDR+6) .NE. 0) THEN
        KABS = IMEM(IHDR+3)
        KLEN = ABS(IMEM(IHDR+6))
        IF (KABS .GE. MXD .OR. KABS+KLEN .LE. MEMIN(3)) THEN
          WRITE (LUNIT, '(A)') ' DMEM data space totally out of range'
        ELSE
          IF (KABS .LT. MEMIN(3)) THEN
            WRITE (LUNIT, '(A)') ' DMEM data space begins out of range'
            KREL = MEMIN(3) - KABS
            KABS = KABS + KREL
            KLEN = KLEN - KREL
            KREL = 1 + KREL
          ELSE
            WRITE (LUNIT, '(A)') ' DMEM'
            KREL = 1
          END IF
          ESKIP = (KABS + KLEN .GT. MXD)
          IF (ESKIP) KLEN = MXD - KABS
          KEND = KABS + KLEN - 1
          LCNT = 1 + (KLEN - 1)/3
          MSKIP = (LCNT .GT. LIMX)
          LDO = MIN (LCNT - 2, LIMX - 2)
          DO 300 I=1,LDO
            WRITE (LUNIT, 60) KABS, KREL, (DMEM(J), J=KABS,KABS+2)
            KABS = KABS + 3
            KREL = KREL + 3
  300     CONTINUE
          IF (MSKIP) THEN
            WRITE (LUNIT, '(A)') ' ... middle section skipped'
            KABS = KABS + 3*(LCNT-LIMX)
            KREL = KREL + 3*(LCNT-LIMX)
          END IF
          IF (LCNT .GT. 1) THEN
            WRITE (LUNIT, 60) KABS, KREL, (DMEM(J), J=KABS,KABS+2)
            KABS = KABS + 3
            KREL = KREL + 3
          END IF
          LDO = KEND - KABS + 1
          WRITE (LUNIT, 60) KABS, KREL, (DMEM(KABS+J-1), J=1,LDO)
          IF (ESKIP) WRITE (LUNIT, '(A)')
     +         ' ... data space ends out of range'
        END IF
      END IF

      RETURN

    5 FORMAT(' IHDR = ',I13)
   10 FORMAT(' IDE:                 IMEM(IHDR+0) =',
     +      I13,' = ',I8,'#',I3.3,2X,A9)
   21 FORMAT(' Pointer to CMEM:     IMEM(IHDR+1) =',I13)
   22 FORMAT(' Pointer to IMEM:     IMEM(IHDR+2) =',I13)
   23 FORMAT(' Pointer to DMEM:     IMEM(IHDR+3) =',I13)
   24 FORMAT(' Length of Character: IMEM(IHDR+4) =',I13)
   25 FORMAT(' Length of Integer:   IMEM(IHDR+5) =',I13)
   26 FORMAT(' Length of Double:    IMEM(IHDR+6) =',I13)
   30 FORMAT(' Link to next header: IMEM(IHDR+7) =',
     +      I13,' = ',I8,'#',I3.3)
   40 FORMAT(I7,I5,1X,A56)
   50 FORMAT(I7,I5,1X,5(1X,A12))
   60 FORMAT(I7,I5,1X,5(1X,G20.12))

      END
