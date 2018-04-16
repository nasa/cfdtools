      SUBROUTINE D0JSFP (CMEM, IMEM, DMEM, IDE, LABEL, ICLS, 
     +      MAXTS, IHDR, IER)

C     PURPOSE:
C        Fetch parameters from a joined surface entity
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
C        ICLS     Closed flag
C                    -1 Closed (shell)
C                    +1 Open
C                    0  Unknown
C        MAXTS    Maximum trimmed surface pointers
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
C                          Joined Surface
      INTEGER     ENTJS
      PARAMETER  (ENTJS = 240)

C =====================================================================

      INTEGER     IDE, ICLS, MAXTS, IHDR, IER
      CHARACTER   LABEL
      INTEGER     JTYP, ITYP, IERX, IA(2), LENC

C     ------

      IER = 0

      LABEL = ' '
      ICLS  = 0
      MAXTS = 0
      IHDR  = 0

C     Call the general pointer-check utility routine

      JTYP = ENTJS

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

      CALL D2FEAI (CMEM, IMEM, DMEM, IDE, 1, 2, IA, IERX)
      IF (IERX .NE. 0) THEN
         IER = -999
         GOTO 9999
      ENDIF

      ICLS  = IA(1)
      MAXTS = IA(2)

 9999 CONTINUE

      RETURN
      END
