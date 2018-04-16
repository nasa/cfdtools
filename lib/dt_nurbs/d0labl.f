      SUBROUTINE D0LABL (CMEM, IMEM, DMEM, IHDR, LABEL, LENLBL)

C     Generate a label from IGES Entity directory entries 18 and 19.
C     This subroutine used during translations from IGES to MEM objects.
C
C   DYNAMIC MEMORY INPUT-OUTPUT
C     CMEM    Dynamically managed character data
C     IMEM    Dynamically managed integer data
C     DMEM    Dynamically managed double precision data
C   INPUT:
C     IHDR    Header index of IGES Entity entity
C   OUTPUT:
C     LABEL   Destination string for label (maximum label is 16 characters)
C     LENLBL  Length of label generated
C
C   CALLS:    none
C
C   HISTORY:
C     24Aug93  P Kraushar   Created
C*****
      CHARACTER CMEM*(*), LABEL*(*)
      INTEGER IMEM(*), IHDR, LENLBL
      DOUBLE PRECISION DMEM(*)

      CHARACTER LABL*8
      INTEGER IXI, IXCB, IXCE, MAXLBL, I
C****
      MAXLBL = LEN (LABEL)
C     Locate IGES directory entry 18, "Entity Label"
      IXCB = IMEM(IHDR+1)+16
      IXCE = IXCB + 7
C     Squeeze out trailing and leading spaces
C     START LOOP
  100    CONTINUE
         IF (CMEM(IXCE:IXCE) .NE. ' ') GOTO 110
         IXCE = IXCE - 1
         IF (IXCE .GE. IXCB) GOTO 100
C     END LOOP
  110 CONTINUE
C     START LOOP
  120    CONTINUE
         IF (IXCB .GT. IXCE) GOTO 130
         IF (CMEM(IXCB:IXCB) .NE. ' ') GOTO 130
         IXCB = IXCB + 1
         GOTO 120
C     END LOOP
  130 CONTINUE
C     Locate IGES directory entry 19, "Entity Subscript"
      IXI = IMEM(IHDR+2) + 18
      IF (IMEM(IXI) .EQ. 0) THEN
C        No subscript
         LENLBL = MIN (IXCE-IXCB+1, MAXLBL)
         IF (IXCB .LE. IXCE) THEN
            LABEL = CMEM(IXCB:IXCE)
         ELSE
            LABEL = ' '
         ENDIF
      ELSE
C        Non-zero subscript
         WRITE (LABL, '(I8)') IMEM(IXI)
         DO 140 I=8,1,-1
            IF (LABL(I:I) .EQ. ' ') THEN
               IXI = I
               LABL(I:I) = '_'
               GOTO 150
            ENDIF
  140    CONTINUE
C        Eight-digit subscript!
         IXI = 1
  150    CONTINUE
         LENLBL = MIN (IXCE-IXCB+10-IXI, MAXLBL)
         IF (IXCB .LE. IXCE) THEN
            LABEL = CMEM(IXCB:IXCE) // LABL(IXI:8)
         ELSE
            LABEL = LABL(IXI:8)
         ENDIF
      ENDIF
      RETURN
      END
