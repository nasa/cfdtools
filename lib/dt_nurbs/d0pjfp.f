      SUBROUTINE D0PJFP (CMEM, IMEM, DMEM, IDE, LABEL, IEG, ITS,
     +      T, U, V, X, Y, Z, IHDR, IER)

C     PURPOSE:
C        Fetch parameters from a point on a joined surface entity
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
C        IEG      Pointer to Edge entity
C        ITS      Pointer to Trimmed Surface Entity
C        T        Parameter value for Edge
C        U,V      Parameter values for Trimmed Surface
C        X,Y,Z    3-Space coordinates
C        IHDR     Index of header block corresponding to IPTR.
C        IER      Error flag.
C                 -1  = IDE does not point to a point on a joined
C                    surface entity
C                 -999= Should not occur
C
C     CALLS:
C        D0PTR
C        D2FELC
C        D2FEEC
C        D2FEAI
C        D2FEAD
C
C     HISTORY:
C        10Jun92  D. Parsons   Created.
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
C                          Point on a Joined Surface
      INTEGER     ENTPJ
      PARAMETER  (ENTPJ = 239)

C =====================================================================

      INTEGER           IDE, IEG, ITS, IHDR, IER
      DOUBLE PRECISION  T, U, V, X, Y, Z
      CHARACTER         LABEL
      INTEGER           JTYP, ITYP, IERX, LENC
      INTEGER           IA(2)
      DOUBLE PRECISION  DA(6)

C     ------

      IER = 0

      LABEL = ' '
      IEG   = 0
      ITS   = 0
      IHDR  = 0
      T     = 0.0D0
      U     = 0.0D0
      V     = 0.0D0
      X     = 0.0D0
      Y     = 0.0D0
      Z     = 0.0D0

C     Call the general pointer-check utility routine

      JTYP = ENTPJ

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

C     Fetch the integer data

      CALL D2FEAI (CMEM, IMEM, DMEM, IDE, 1, 2, IA, IERX)
      IF (IERX .NE. 0) THEN
         IER = -999
         GOTO 9999
      ENDIF

      IEG = IA(1)
      ITS = IA(2)

C     Fetch the double precision data

      CALL D2FEAD (CMEM, IMEM, DMEM, IDE, 1, 6, DA, IERX)
      IF (IERX .NE. 0) THEN
         IER = -999
         GOTO 9999
      ENDIF

      T = DA(1)
      U = DA(2)
      V = DA(3)
      X = DA(4)
      Y = DA(5)
      Z = DA(6)

 9999 CONTINUE

      RETURN
      END
