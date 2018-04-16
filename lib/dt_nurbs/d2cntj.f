      SUBROUTINE D2CNTJ (CMEM, IMEM, DMEM, IXTS, ITS, IJS, IER)

C     PURPOSE:
C        ConNect Trimmed Surface to Joined Surface.  Record
C        Trimmed Surface ITS as the IXTSth Trimmed Surface in
C        Joined Surface IJS.
C
C        If ITS is zero, the IJS(IXTS) is set to zero.
C
C        If ITS has a non-zero Joined Surface pointer, or IJS already
C        has a Trimmed Surface defined in its IXTSth index, an attempt
C        will be made to remove those connections.
C
C     DYNAMIC MEMORY INPUT-OUTPUT:
C        CMEM     Dynamically managed character data
C        IMEM     Dynamically managed integer data
C        DMEM     Dynamically managed double precision data
C
C     INPUT:
C        IXTS     Index value of the Trimmed Surface.
C        ITS      Pointer to the Trimmed Surface.
C        IJS      Pointer to the Joined Surface.
C
C     OUTPUT:
C        IER      Error flag.
C                 -1  = Dynamic memory is corrupt or uninitialized
C                 -2  = IJS does not point to an Trimmed Surface entity
C                 -3  = IXTS not between 1 and MAXTS
C                 -4  = ITS does not point to a Trimmed Surface entity
C                 -999= Unexpected error in a lower level routine
C
C     CALLS:
C        D1MBAD
C        D0TSFP
C        D0JSFP
C        D2FEEI
C        D2STEI
C        DTERR
C
C     HISTORY:
C        09Jul92  D. Parsons  Created.
C        21Jul92  D. Parsons  Added zeroing of old connections.
C
C     ------

C     Long name alias:
C        ENTRY D2_CONNECT_TRIMMED_JOINED (CMEM, IMEM, DMEM, IXTS, ITS,
C    +         IJS, IER)

      EXTERNAL D1MBAD
      LOGICAL  D1MBAD

      CHARACTER         CMEM*(*)
      INTEGER           IMEM(*)
      DOUBLE PRECISION  DMEM(*)

      INTEGER           IXTS, ITS, IJS, IER
      CHARACTER         LABELX
      INTEGER           IDMY, MAXTS, ITSX, IXTSX, IJSX, IERX

      CHARACTER*6 SUBNAM
      DATA SUBNAM /'D2CNTJ'/

C     ------

      IER = 0

      IF (CMEM(1:1) .NE. 'M' .AND. CMEM(1:1) .NE. 'L') THEN
         IF (D1MBAD (CMEM, IMEM, DMEM, 0, 0, 0, IERX)) THEN
            IER = -1
            GOTO 9000
         ENDIF
      ENDIF

C     Verify the validity of the Joined Surface and get MAXTS.  Try
C     to disconnect any old connections

      IF (IJS .NE. 0) THEN
         CALL D0JSFP (CMEM, IMEM, DMEM, IJS, LABELX, IDMY,
     +         MAXTS, IDMY, IERX)
         IF (IERX .NE. 0) THEN
            IER = -2
            GOTO 9000
         ENDIF
         IF ((IXTS .LT. 1) .OR. (IXTS .GT. MAXTS)) THEN
            IER = -3
            GOTO 9000
         ENDIF

C        Look for old connections

         CALL D2FEEI (CMEM, IMEM, DMEM, IJS, 2+IXTS, ITSX, IERX)
         IF (IERX .NE. 0) THEN
            IER = -999
            GOTO 9000
         ENDIF
         IF (ITSX .NE. 0) THEN
            CALL D0TSFP (CMEM, IMEM, DMEM, ITSX, LABELX, IDMY, IJSX, 
     +            IDMY, IDMY, IDMY, IDMY, IERX)
            IF (IJSX .EQ. IJS) THEN
C              Zero the Joined Surface pointer
               CALL D2STEI (CMEM, IMEM, DMEM, 0, 2, ITSX, IERX)
               IF (IERX .NE. 0) THEN
                  IER = -999
                  GOTO 9000
               ENDIF
C              Zero the index
               CALL D2STEI (CMEM, IMEM, DMEM, 0, 3, ITSX, IERX)
               IF (IERX .NE. 0) THEN
                  IER = -999
                  GOTO 9000
               ENDIF
            ENDIF
         ENDIF
      ENDIF

      IF (ITS .NE. 0) THEN

C        Verify the validity of the Trimmed Surface

         CALL D0TSFP (CMEM, IMEM, DMEM, ITS, LABELX, IDMY, IJSX, IXTSX,
     +      IDMY, IDMY, IDMY, IERX)
         IF (IERX .NE. 0) THEN
            IER = -4
            GOTO 9000
         ENDIF

C        If the Trimmed Surface is already pointing to a Joined Surface,
C        try to disconnect the link

         IF (IJSX .NE. 0) THEN
            CALL D0JSFP (CMEM, IMEM, DMEM, IJSX, LABELX, IDMY, IDMY,
     +            IDMY, IERX)
            IF (IERX .EQ. 0) THEN
               CALL D2FEEI (CMEM, IMEM, DMEM, IJSX, 2+IXTSX, ITSX, IERX)
               IF (ITSX .EQ. ITS)
     +               CALL D2STEI (CMEM, IMEM, DMEM, 0, 2+IXTSX,
     +               IJSX, IERX)
            ENDIF
         ENDIF
      ENDIF

C     Replace the IXTSth Trimmed Surface in IJS with ITS

      IF (IJS .NE. 0) THEN
         CALL D2STEI (CMEM, IMEM, DMEM, ITS, 2+IXTS, IJS, IERX)
         IF (IERX .NE. 0) THEN
            IER = -999
            GOTO 9000
         ENDIF
      ENDIF

C     Set the Trimmed Surface ITS to point back to this IJS

      IF (ITS .NE. 0) THEN
C        Store the Joined Surface pointer into ITS
         CALL D2STEI (CMEM, IMEM, DMEM, IJS, 2, ITS, IERX)
         IF (IERX .NE. 0) THEN
            IER = -999
            GOTO 9000
         ENDIF
C        Store the index IXTS pointer into ITS
         IF (IJS .NE. 0) THEN
            CALL D2STEI (CMEM, IMEM, DMEM, IXTS, 3, ITS, IERX)
         ELSE
C           Zero the index pointer when zeroing the IJS
            CALL D2STEI (CMEM, IMEM, DMEM, 0, 3, ITS, IERX)
         ENDIF
         IF (IERX .NE. 0) THEN
            IER = -999
            GOTO 9000
         ENDIF
      ENDIF

      GOTO 9999

 9000 CONTINUE

C     Error Handling

      IF (IER .EQ. -999) THEN
         CALL DTERR (5, SUBNAM, IER, 0)
      ELSE
         CALL DTERR (1, SUBNAM, IER, 0)
      ENDIF

 9999 CONTINUE

      RETURN
      END
