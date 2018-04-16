      SUBROUTINE D0TSFP (CMEM, IMEM, DMEM, IDE, LABEL, IBFS, IJS, IXTS,
     +      IORIEN, MAXLP, IHDR, IER)

C     PURPOSE:
C        Fetch parameters from a trimmed surface entity
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
C        IBFS     Pointer to survace B-spline entity
C        IJS      Pointer to Joined Surface entity
C        IXTS     Index of Trimmed surfaces in IJS
C        IORIEN   Orientation
C        MAXLP    Maximum number of Loop entities
C        IHDR     Index of header block corresponding to IPTR.
C        IER      Error flag.
C                 -1  = IDE does not point to a trimmed surface entity
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
C                          Trimmed Surface of a Joined Surface
      INTEGER     ENTTS
      PARAMETER  (ENTTS = 241)

C =====================================================================

      INTEGER     IDE, IBFS, IJS, IXTS, IORIEN, MAXLP, IHDR, IER
      CHARACTER   LABEL
      INTEGER     JTYP, ITYP, IERX, IA(5), LENC

C     ------

      IER = 0

      LABEL    = ' '
      IBFS     = 0
      IJS      = 0
      IXTS     = 0
      IORIEN   = 0
      MAXLP    = 0
      IHDR     = 0

C     Call the general pointer-check utility routine

      JTYP = ENTTS

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

      CALL D2FEAI (CMEM, IMEM, DMEM, IDE, 1, 5, IA, IERX)
      IF (IERX .NE. 0) THEN
         IER = -999
         GOTO 9999
      ENDIF

      IBFS   = IA(1)
      IJS    = IA(2)
      IXTS   = IA(3)
      IORIEN = IA(4)
      MAXLP  = IA(5)

 9999 CONTINUE

      RETURN
      END
