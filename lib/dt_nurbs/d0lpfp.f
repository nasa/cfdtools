      SUBROUTINE D0LPFP (CMEM, IMEM, DMEM, IDE, LABEL, ITS, IXLP,
     +      MAXEG, IHDR, IER)

C     PURPOSE:
C        Fetch parameters from a loop entity
C
C
C     DYNAMIC MEMORY INPUT-OUTPUT:
C        CMEM     Dynamically managed character data
C        IMEM     Dynamically managed integer data
C        DMEM     Dynamically managed double precision data
C
C     INPUT:
C        IDE      Entity ID
C
C     OUTPUT:
C        LABEL    The first character of the label
C        ITS      Pointer to trimmed surface entity
C        IXLP     Index into ITS
C        MAXEG    Maximum edge pointers
C        IHDR     Index of header block corresponding to IPTR.
C        IER      Error flag.
C                 -1  = IDE does not point to a loop entity
C                 -999= Should not occur
C
C     CALLS:
C        D0PTR
C
C     HISTORY:
C        19May92  D. Parsons   Created.
C
C     ------

      CHARACTER         CMEM*(*)
      INTEGER           IMEM(*)
      DOUBLE PRECISION  DMEM(*)

C     Get the Reserved Entity Type ID Numbers

C =====================================================================
C
C     INCLUDE FILE DITYID.INC
C
C     Reserved Entity Type ID Numbers
C
C     HISTORY
C        12May92     D. Parsons     created
C
C                          Loop of a Trimmed Surface
      INTEGER     ENTLP
      PARAMETER  (ENTLP = 242)

C =====================================================================

      INTEGER     IDE, IXLP, MAXEG, IHDR, IER
      CHARACTER   LABEL
      INTEGER     JTYP, ITYP, IERX, IA(3), LENC

C     ------

      IER = 0

      LABEL = ' '
      ITS   = 0
      IXLP  = 0
      MAXEG = 0
      IHDR  = 0

C     Call the general pointer-check utility routine

      JTYP = ENTLP

      CALL D0PTR (CMEM, IMEM, DMEM, IDE, JTYP, 0, IHDR, ITYP, IERX)

      IF (IERX .NE. 0) THEN
         IER = -1
         GOTO 9999
      ENDIF

C     Extract the information

      CALL D2FELC (CMEM, IMEM, DMEM, IDE, LENC, IERX)
      IF (IERX .NE. 0) THEN
         IER = -999
         GOTO 9999
      ENDIF

C     Fetch the first character of the label

      IF (LENC .GT. 0) THEN
         CALL D2FEEC (CMEM, IMEM, DMEM, IDE, 1, LABEL, IERX)
         IF (IERX .NE. 0) THEN
            IER = -999
            GOTO 9999
         ENDIF
      ELSE
         LABEL = ' '
      ENDIF

      CALL D2FEAI (CMEM, IMEM, DMEM, IDE, 1, 3, IA, IERX)
      IF (IERX .NE. 0) THEN
         IER = -999
         GOTO 9999
      ENDIF

      ITS   = IA(1)
      IXLP  = IA(2)
      MAXEG = IA(3)

 9999 CONTINUE

      RETURN
      END
